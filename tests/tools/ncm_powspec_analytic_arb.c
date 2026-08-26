/***************************************************************************
 *            ncm_powspec_analytic_arb.c
 *
 *  Tue August 26 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Certified reference values for #NcmPowspecAnalytic.
 *
 * Re-implements the closed form of ncm_powspec_analytic.c in Arb ball
 * arithmetic, independently of the library. That duplication is the point: an
 * independent re-derivation is what makes this a reference rather than a
 * tautology. Every value comes with a rigorous radius, so a test can assert
 * against a number whose accuracy is proved rather than assumed.
 *
 * Working precision is not fixed. Each entry is recomputed at doubling
 * precision until its relative radius falls below the target, which is what
 * makes the tool robust where a fixed-precision one would silently return
 * garbage -- see the growth factor below the Limber band in the k-integral
 * plan.
 *
 * Emits TSV on stdout: shape growth k z value radius prec.
 *
 * Build (standalone, when meson has not been told about FLINT):
 *   gcc -O2 -o powspec_arb ncm_powspec_analytic_arb.c \
 *       $(pkg-config --cflags --libs flint) -lm
 */

/* arb.h sits at the top level in some packagings and under flint/ in others;
 * the top-level meson.build probes for both. __has_include keeps this file
 * compilable standalone too, where that probe's defines are absent. */
#if defined (__has_include )
#if __has_include (<flint/arb.h>)
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#else
#include <arb.h>
#include <arb_hypgeom.h>
#endif
#else
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BBKS_C1 2.34
#define BBKS_C2 3.89
#define BBKS_C3 16.1
#define BBKS_C4 5.46
#define BBKS_C5 6.71

typedef enum
{
  SHAPE_POWER_LAW = 0,
  SHAPE_BBKS,
  SHAPE_RATIONAL,
} Shape;

typedef enum
{
  GROWTH_NONE = 0,
  GROWTH_LCDM,
  GROWTH_RATIONAL,
} Growth;

typedef struct
{
  Shape shape;
  Growth growth;
  double amplitude;
  double n_s;
  double k_eq;
  double Omega_m;
  double a_t;
  double bao_amplitude;
  double bao_rd;
  double bao_sigma;
} Par;

/* T(k). */
static void
transfer (arb_t res, const Par *p, const arb_t k, slong prec)
{
  arb_t q, t, u, M, L;

  arb_init (q);
  arb_init (t);
  arb_init (u);
  arb_init (M);
  arb_init (L);

  arb_set_d (t, p->k_eq);
  arb_div (q, k, t, prec); /* q = k / k_eq */

  switch (p->shape)
  {
    case SHAPE_POWER_LAW:
      arb_one (res);
      break;

    case SHAPE_BBKS:
      /* M = 1 + C2 q + (C3 q)^2 + (C4 q)^3 + (C5 q)^4, by Horner in q. */
      arb_set_d (M, BBKS_C5 * BBKS_C5 * BBKS_C5 * BBKS_C5);
      arb_mul (M, M, q, prec);
      arb_set_d (t, BBKS_C4 * BBKS_C4 * BBKS_C4);
      arb_add (M, M, t, prec);
      arb_mul (M, M, q, prec);
      arb_set_d (t, BBKS_C3 * BBKS_C3);
      arb_add (M, M, t, prec);
      arb_mul (M, M, q, prec);
      arb_set_d (t, BBKS_C2);
      arb_add (M, M, t, prec);
      arb_mul (M, M, q, prec);
      arb_add_si (M, M, 1, prec);

      /* L = log1p (C1 q) / (C1 q), with the removable 0/0 at q = 0. */
      arb_set_d (t, BBKS_C1);
      arb_mul (u, q, t, prec);

      if (arb_is_zero (u))
      {
        arb_one (L);
      }
      else
      {
        arb_log1p (L, u, prec);
        arb_div (L, L, u, prec);
      }

      arb_set_d (t, -0.25);
      arb_pow (M, M, t, prec);
      arb_mul (res, L, M, prec);
      break;

    case SHAPE_RATIONAL:
      arb_sqr (t, q, prec);
      arb_add_si (t, t, 1, prec);
      arb_inv (res, t, prec);
      break;

    default:
      fprintf (stderr, "transfer: bad shape\n");
      exit (1);
  }

  arb_clear (q);
  arb_clear (t);
  arb_clear (u);
  arb_clear (M);
  arb_clear (L);
}

/* B(k) = 1 + a sinc (k r_d) exp (-(k sigma_d)^2), one when a = 0. */
static void
bao (arb_t res, const Par *p, const arb_t k, slong prec)
{
  arb_t x, s, d, t;

  if (p->bao_amplitude == 0.0)
  {
    arb_one (res);

    return;
  }

  arb_init (x);
  arb_init (s);
  arb_init (d);
  arb_init (t);

  arb_set_d (t, p->bao_rd);
  arb_mul (x, k, t, prec);

  if (arb_is_zero (x))
  {
    arb_one (s);
  }
  else
  {
    arb_sin (s, x, prec);
    arb_div (s, s, x, prec);
  }

  arb_set_d (t, p->bao_sigma);
  arb_mul (d, k, t, prec);
  arb_sqr (d, d, prec);
  arb_neg (d, d);
  arb_exp (d, d, prec);

  arb_mul (s, s, d, prec);
  arb_set_d (t, p->bao_amplitude);
  arb_mul (s, s, t, prec);
  arb_add_si (res, s, 1, prec);

  arb_clear (x);
  arb_clear (s);
  arb_clear (d);
  arb_clear (t);
}

/*
 * Unnormalized growth. For LCDM this is
 *
 *   D(a) = a 2F1 (1/3, 1; 11/6; -(1 - Om)/Om a^3),
 *
 * evaluated directly: unlike GSL, Arb has no |z| < 1 restriction, so the
 * Pfaff transformation the library uses is not needed here. Comparing the two
 * therefore also checks that transformation.
 */
static void
growth_raw (arb_t res, const Par *p, const arb_t a, slong prec)
{
  arb_t x, t, A, B, C;

  arb_init (x);
  arb_init (t);
  arb_init (A);
  arb_init (B);
  arb_init (C);

  switch (p->growth)
  {
    case GROWTH_NONE:
      arb_one (res);
      break;

    case GROWTH_LCDM:
      arb_pow_ui (x, a, 3, prec);
      arb_set_d (t, -(1.0 - p->Omega_m) / p->Omega_m);
      arb_mul (x, x, t, prec);

      arb_set_ui (A, 1);
      arb_div_ui (A, A, 3, prec); /* 1/3  */
      arb_one (B);                /* 1    */
      arb_set_ui (C, 11);
      arb_div_ui (C, C, 6, prec); /* 11/6 */

      arb_hypgeom_2f1 (res, A, B, C, x, 0, prec);
      arb_mul (res, res, a, prec);
      break;

    case GROWTH_RATIONAL:
      arb_set_d (t, p->a_t);
      arb_div (x, a, t, prec);
      arb_pow_ui (x, x, 3, prec);
      arb_add_si (x, x, 1, prec);
      arb_set_d (t, -1.0 / 3.0);
      arb_pow (x, x, t, prec);
      arb_mul (res, a, x, prec);
      break;

    default:
      fprintf (stderr, "growth_raw: bad growth\n");
      exit (1);
  }

  arb_clear (x);
  arb_clear (t);
  arb_clear (A);
  arb_clear (B);
  arb_clear (C);
}

/* P(k, z) = A k^n_s T(k)^2 B(k) D(z)^2, with D(0) = 1. */
static void
powspec (arb_t res, const Par *p, double kd, double zd, slong prec)
{
  arb_t k, a, one, T, Bk, D, Dn, t;

  arb_init (k);
  arb_init (a);
  arb_init (one);
  arb_init (T);
  arb_init (Bk);
  arb_init (D);
  arb_init (Dn);
  arb_init (t);

  arb_set_d (k, kd);
  arb_set_d (t, 1.0 + zd);
  arb_one (a);
  arb_div (a, a, t, prec);
  arb_one (one);

  transfer (T, p, k, prec);
  bao (Bk, p, k, prec);
  growth_raw (D, p, a, prec);
  growth_raw (Dn, p, one, prec);
  arb_div (D, D, Dn, prec);

  arb_set_d (t, p->n_s);
  arb_pow (res, k, t, prec);
  arb_sqr (T, T, prec);
  arb_mul (res, res, T, prec);
  arb_mul (res, res, Bk, prec);
  arb_sqr (D, D, prec);
  arb_mul (res, res, D, prec);
  arb_set_d (t, p->amplitude);
  arb_mul (res, res, t, prec);

  arb_clear (k);
  arb_clear (a);
  arb_clear (one);
  arb_clear (T);
  arb_clear (Bk);
  arb_clear (D);
  arb_clear (Dn);
  arb_clear (t);
}

/*
 * Recompute at doubling precision until the relative radius is below
 * @target_rel. Returns the precision that succeeded.
 */
static slong
powspec_certified (arb_t res, const Par *p, double k, double z, double target_rel)
{
  slong prec;

  for (prec = 128; prec <= 8192; prec *= 2)
  {
    double r, m;

    powspec (res, p, k, z, prec);

    if (arb_is_zero (res))
      return prec;

    r = mag_get_d (arb_radref (res));
    m = arf_get_d (arb_midref (res), ARF_RND_NEAR);

    if ((m != 0.0) && (r / fabs (m) < target_rel))
      return prec;
  }

  fprintf (stderr, "powspec_certified: did not reach %g at k=%g z=%g\n", target_rel, k, z);
  exit (1);
}

static Shape
parse_shape (const char *s)
{
  if (strcmp (s, "power_law") == 0)
    return SHAPE_POWER_LAW;

  if (strcmp (s, "bbks") == 0)
    return SHAPE_BBKS;

  if (strcmp (s, "rational") == 0)
    return SHAPE_RATIONAL;

  fprintf (stderr, "unknown shape '%s'\n", s);
  exit (1);
}

static Growth
parse_growth (const char *s)
{
  if (strcmp (s, "none") == 0)
    return GROWTH_NONE;

  if (strcmp (s, "lcdm") == 0)
    return GROWTH_LCDM;

  if (strcmp (s, "rational") == 0)
    return GROWTH_RATIONAL;

  fprintf (stderr, "unknown growth '%s'\n", s);
  exit (1);
}

int
main (int argc, char **argv)
{
  /* Defaults match ncm_powspec_analytic.c's property defaults. */
  Par p = {
    .shape     = SHAPE_BBKS, .growth = GROWTH_LCDM,
    .amplitude = 2.08e7, .n_s = 0.9875, .k_eq = 0.10594, .Omega_m = 0.3,
    .a_t       = 1.0, .bao_amplitude = 0.0, .bao_rd = 147.0, .bao_sigma = 10.0
  };
  double target_rel = 1.0e-25;
  double lk_min = -5.0, lk_max = 1.0;
  int n_k     = 61;
  double zs[] = { 0.0, 0.1, 0.5, 1.0, 3.0, 10.0, 20.0 };
  int n_z     = (int) (sizeof (zs) / sizeof (zs[0]));
  int i, j;

  for (i = 1; i < argc; i++)
  {
    if (strncmp (argv[i], "--shape=", 8) == 0)
    {
      p.shape = parse_shape (argv[i] + 8);
    }
    else if (strncmp (argv[i], "--growth=", 9) == 0)
    {
      p.growth = parse_growth (argv[i] + 9);
    }
    else if (strncmp (argv[i], "--amplitude=", 12) == 0)
    {
      p.amplitude = atof (argv[i] + 12);
    }
    else if (strncmp (argv[i], "--n-s=", 6) == 0)
    {
      p.n_s = atof (argv[i] + 6);
    }
    else if (strncmp (argv[i], "--k-eq=", 7) == 0)
    {
      p.k_eq = atof (argv[i] + 7);
    }
    else if (strncmp (argv[i], "--Omega-m=", 10) == 0)
    {
      p.Omega_m = atof (argv[i] + 10);
    }
    else if (strncmp (argv[i], "--a-t=", 6) == 0)
    {
      p.a_t = atof (argv[i] + 6);
    }
    else if (strncmp (argv[i], "--bao-amplitude=", 16) == 0)
    {
      p.bao_amplitude = atof (argv[i] + 16);
    }
    else if (strncmp (argv[i], "--n-k=", 6) == 0)
    {
      n_k = atoi (argv[i] + 6);
    }
    else if (strncmp (argv[i], "--target-rel=", 13) == 0)
    {
      target_rel = atof (argv[i] + 13);
    }
    else
    {
      fprintf (stderr, "unknown argument '%s'\n", argv[i]);
      exit (1);
    }
  }

  printf ("# shape\tgrowth\tk\tz\tvalue\tradius\tprec\n");

  for (i = 0; i < n_k; i++)
  {
    const double k = pow (10.0, lk_min + (lk_max - lk_min) * i / (double) (n_k - 1));

    for (j = 0; j < n_z; j++)
    {
      arb_t res;
      slong prec;
      char *s;

      arb_init (res);
      prec = powspec_certified (res, &p, k, zs[j], target_rel);
      s    = arb_get_str (res, 30, ARB_STR_NO_RADIUS);

      printf ("%d\t%d\t%.17g\t%.17g\t%s\t%.6e\t%ld\n",
              (int) p.shape, (int) p.growth, k, zs[j], s,
              mag_get_d (arb_radref (res)), (long) prec);

      flint_free (s);
      arb_clear (res);
    }
  }

  flint_cleanup ();

  return 0;
}


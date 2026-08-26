/***************************************************************************
 *            nc_xcor_kernel_analytic_arb.c
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
 * Certified reference values for the radial integral of the analytic xcor
 * windows,
 *
 *   I_ell(k) = int W(chi) j_ell(k chi) dchi ,   W normalized to int W dchi = 1
 *
 * over the window's own truncated support, with chi in Mpc.
 *
 * Re-implements each window of numcosmo/nc/xcor/tests/nc_xcor_kernel_analytic_*.c
 * in Arb ball arithmetic, independently of the library. The normalization is
 * *not* copied from the library's closed form -- it is computed here as a
 * certified integral of the unnormalized shape, so comparing against the
 * library checks its normalization too.
 *
 * Working precision is not fixed: every entry is recomputed at doubling
 * precision until its relative radius clears the target. That is what makes
 * the tool usable in the deep sub-band regime, where a fixed-precision
 * reference returns a plausible wrong number instead of reporting failure.
 *
 * Emits TSV on stdout: shape ell k value radius prec.
 *
 * Build (standalone):
 *   gcc -O2 -o xcor_window_arb nc_xcor_kernel_analytic_arb.c \
 *       $(pkg-config --cflags --libs flint) -lm
 */

#if defined (__has_include )
#if __has_include (<flint/acb.h>)
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_calc.h>
#else
#include <acb.h>
#include <acb_hypgeom.h>
#include <acb_calc.h>
#endif
#else
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_calc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_BUMPS 8
#define MAX_PANELS 8

typedef enum
{
  SHAPE_GAUSS = 0,
  SHAPE_TOPHAT,
  SHAPE_TOPHAT_SMOOTH,
  SHAPE_STUDENT_T,
  SHAPE_POWER_EXP,
  SHAPE_LENSING,
  SHAPE_MULTI,
  SHAPE_LEN,
} Shape;

static const char *shape_names[SHAPE_LEN] = {
  "gauss", "tophat", "tophat_smooth", "student_t", "power_exp", "lensing", "multi"
};

typedef struct
{
  Shape shape;
  long ell;
  double k;

  /* gauss / student_t / multi */
  double chi_mean, chi_sigma, n_sigma;
  /* tophat / tophat_smooth / power_exp */
  double chi_lower, chi_upper;
  /* student_t */
  double chi_scale, nu, n_scale;
  /* power_exp */
  double alpha, beta;
  /* lensing */
  double chi_source_lower, chi_source_upper;
  /* multi */
  int n_bumps;
  double mu[MAX_BUMPS], sg[MAX_BUMPS], wt[MAX_BUMPS];

  /* support, filled by shape_support () */
  double chi_min, chi_max;
  /* internal breakpoints where the shape is not analytic */
  int n_break;
  double breaks[MAX_PANELS];

  /* set while integrating: 0 = bare window, 1 = window times j_ell */
  int with_bessel;
  /* 0 = emit I_ell(k) for k on stdin, 1 = emit W(chi) for chi on stdin */
  int window_mode;
} Par;

/*
 * The unnormalized window, on a complex ball.
 *
 * Returns non-analytic (an indeterminate ball) when @order is 1 and the input
 * ball may cross a branch cut, which is how acb_calc_integrate is told to
 * subdivide rather than trust the value. Only `lensing` and `power_exp` have a
 * cut at all -- both take a log or a real power of chi, analytic for Re chi > 0
 * -- so the guard is a positivity test on the real part.
 */
static int
window_u (acb_t out, const acb_t chi, Par *p, slong order, slong prec)
{
  acb_t t, u, v, w;
  int i;

  acb_init (t);
  acb_init (u);
  acb_init (v);
  acb_init (w);

  /* Analyticity guard for the shapes that take log/pow of chi. */
  if ((order > 0) &&
      ((p->shape == SHAPE_LENSING) || (p->shape == SHAPE_POWER_EXP)))
  {
    arb_t re;

    arb_init (re);
    acb_get_real (re, chi);

    if (!arb_is_positive (re))
    {
      acb_indeterminate (out);
      arb_clear (re);
      acb_clear (t);
      acb_clear (u);
      acb_clear (v);
      acb_clear (w);

      return 0;
    }

    arb_clear (re);
  }

  switch (p->shape)
  {
    case SHAPE_GAUSS:
      acb_set_d (t, p->chi_mean);
      acb_sub (u, chi, t, prec);
      acb_set_d (t, p->chi_sigma);
      acb_div (u, u, t, prec);
      acb_sqr (u, u, prec);
      acb_div_si (u, u, -2, prec);
      acb_exp (out, u, prec);
      break;

    case SHAPE_TOPHAT:
      acb_one (out);
      break;

    case SHAPE_TOPHAT_SMOOTH:
      /* erf ((chi_upper - chi)/(sqrt2 sigma)) - erf ((chi_lower - chi)/(sqrt2 sigma)) */
      acb_set_d (w, M_SQRT2 * p->chi_sigma);

      acb_set_d (t, p->chi_upper);
      acb_sub (u, t, chi, prec);
      acb_div (u, u, w, prec);
      acb_hypgeom_erf (u, u, prec);

      acb_set_d (t, p->chi_lower);
      acb_sub (v, t, chi, prec);
      acb_div (v, v, w, prec);
      acb_hypgeom_erf (v, v, prec);

      acb_sub (out, u, v, prec);
      break;

    case SHAPE_STUDENT_T:
      /* (1 + t^2/nu)^(-(nu+1)/2),  t = (chi - mean)/scale */
      acb_set_d (t, p->chi_mean);
      acb_sub (u, chi, t, prec);
      acb_set_d (t, p->chi_scale);
      acb_div (u, u, t, prec);
      acb_sqr (u, u, prec);
      acb_set_d (t, p->nu);
      acb_div (u, u, t, prec);
      acb_add_si (u, u, 1, prec);
      acb_set_d (v, -0.5 * (p->nu + 1.0));
      acb_pow (out, u, v, prec);
      break;

    case SHAPE_POWER_EXP:
      /* x^alpha exp (-x^beta),  x = chi / chi_scale */
      acb_set_d (t, p->chi_scale);
      acb_div (u, chi, t, prec);
      acb_set_d (v, p->alpha);
      acb_pow (w, u, v, prec);
      acb_set_d (v, p->beta);
      acb_pow (u, u, v, prec);
      acb_neg (u, u);
      acb_exp (u, u, prec);
      acb_mul (out, w, u, prec);
      break;

    case SHAPE_LENSING:
    {
      /* chi ((b - m) - chi log (b/m)) / (b - a), with m = max (chi, a).
       * The two branches meet at chi = a, where the shape has a kink; the
       * caller panels on it, so each call sees only one branch. */
      const double a = p->chi_source_lower;
      const double b = p->chi_source_upper;
      arb_t re, av;
      int upper;

      arb_init (re);
      arb_init (av);
      acb_get_real (re, chi);
      arb_set_d (av, a);
      upper = arb_gt (re, av);
      arb_clear (re);
      arb_clear (av);

      if (upper)
        acb_set (v, chi);  /* m = chi */
      else
        acb_set_d (v, a);  /* m = a   */

      acb_set_d (t, b);
      acb_div (u, t, v, prec); /* b / m */
      acb_log (u, u, prec);
      acb_mul (u, u, chi, prec); /* chi log (b/m) */

      acb_sub (w, t, v, prec); /* b - m */
      acb_sub (w, w, u, prec);
      acb_mul (w, w, chi, prec);

      acb_set_d (t, b - a);
      acb_div (out, w, t, prec);
      break;
    }

    case SHAPE_MULTI:
      acb_zero (out);

      for (i = 0; i < p->n_bumps; i++)
      {
        acb_set_d (t, p->mu[i]);
        acb_sub (u, chi, t, prec);
        acb_set_d (t, p->sg[i]);
        acb_div (u, u, t, prec);
        acb_sqr (u, u, prec);
        acb_div_si (u, u, -2, prec);
        acb_exp (u, u, prec);
        acb_set_d (t, p->wt[i]);
        acb_mul (u, u, t, prec);
        acb_add (out, out, u, prec);
      }

      break;

    default:
      fprintf (stderr, "window_u: bad shape\n");
      exit (1);
  }

  acb_clear (t);
  acb_clear (u);
  acb_clear (v);
  acb_clear (w);

  return 0;
}

/* The integrand handed to acb_calc_integrate: W_u, optionally times j_ell. */
static int
integrand (acb_ptr out, const acb_t chi, void *param, slong order, slong prec)
{
  Par *p = (Par *) param;
  acb_t W, z, nu, J, t;

  if (order > 1)
    return 0;

  acb_init (W);
  window_u (W, chi, p, order, prec);

  if (!p->with_bessel)
  {
    acb_set (out, W);
    acb_clear (W);

    return 0;
  }

  acb_init (z);
  acb_init (nu);
  acb_init (J);
  acb_init (t);

  acb_set_d (t, p->k);
  acb_mul (z, chi, t, prec);

  /*
   * j_ell (z) = sqrt (pi) / 2^(ell+1) * z^ell * 0F1~ (; ell + 3/2; -z^2/4),
   * with 0F1~ the regularized confluent limit.
   *
   * Not the textbook sqrt (pi / 2z) J_{ell+1/2} (z): that form has a removable
   * singularity at z = 0, and a window whose support reaches the observer
   * (`multi` clamps to chi_min = 0) puts z = 0 inside the domain, where the
   * ball arithmetic sees 1/0 and subdivides without end. This form is entire.
   */
  acb_set_si (nu, 2 * p->ell + 3);
  acb_div_si (nu, nu, 2, prec); /* ell + 3/2 */
  acb_sqr (J, z, prec);
  acb_div_si (J, J, -4, prec); /* -z^2/4 */
  acb_hypgeom_0f1 (J, nu, J, 1, prec);
  acb_pow_ui (t, z, (ulong) p->ell, prec);
  acb_mul (J, J, t, prec);
  acb_const_pi (t, prec);
  acb_sqrt (t, t, prec);
  acb_mul (J, J, t, prec);
  acb_mul_2exp_si (J, J, -(p->ell + 1));

  acb_mul (out, W, J, prec);

  acb_clear (W);
  acb_clear (z);
  acb_clear (nu);
  acb_clear (J);
  acb_clear (t);

  return 0;
}

/* Sum of the certified integral over the shape's panels, at one precision. */
static void
integrate_panels (acb_t res, Par *p, slong prec)
{
  double edges[MAX_PANELS + 2];
  acb_calc_integrate_opt_t opt;
  acb_t A, B, part;
  mag_t tol;
  int n = 0, i;

  edges[n++] = p->chi_min;

  for (i = 0; i < p->n_break; i++)
    if ((p->breaks[i] > p->chi_min) && (p->breaks[i] < p->chi_max))
      edges[n++] = p->breaks[i];

  edges[n++] = p->chi_max;

  acb_calc_integrate_opt_init (opt);
  acb_init (A);
  acb_init (B);
  acb_init (part);
  mag_init (tol);
  mag_set_ui_2exp_si (tol, 1, -prec);

  acb_zero (res);

  for (i = 0; i + 1 < n; i++)
  {
    acb_set_d (A, edges[i]);
    acb_set_d (B, edges[i + 1]);
    acb_calc_integrate (part, integrand, p, A, B, prec, tol, opt, prec);
    acb_add (res, res, part, prec);
  }

  acb_clear (A);
  acb_clear (B);
  acb_clear (part);
  mag_clear (tol);
}

/* Recompute at doubling precision until the relative radius clears @target. */
static slong
certified (acb_t res, Par *p, double target)
{
  slong prec;

  for (prec = 128; prec <= 8192; prec *= 2)
  {
    double r, m;

    integrate_panels (res, p, prec);

    if (acb_is_zero (res))
      return prec;

    r = mag_get_d (arb_radref (acb_realref (res)));
    m = arf_get_d (arb_midref (acb_realref (res)), ARF_RND_NEAR);

    if ((m != 0.0) && (r / fabs (m) < target))
      return prec;
  }

  fprintf (stderr, "certified: did not reach %g (shape %s ell %ld k %g)\n",
           target, shape_names[p->shape], p->ell, p->k);
  exit (1);
}

/* Support and internal breakpoints, matching the library's constructors. */
static void
shape_support (Par *p)
{
  int i;

  p->n_break = 0;

  switch (p->shape)
  {
    case SHAPE_GAUSS:
      p->chi_min = fmax (0.0, p->chi_mean - p->n_sigma * p->chi_sigma);
      p->chi_max = p->chi_mean + p->n_sigma * p->chi_sigma;
      break;

    case SHAPE_TOPHAT:
      p->chi_min = p->chi_lower;
      p->chi_max = p->chi_upper;
      break;

    case SHAPE_TOPHAT_SMOOTH:
      p->chi_min = fmax (0.0, p->chi_lower - p->n_sigma * p->chi_sigma);
      p->chi_max = p->chi_upper + p->n_sigma * p->chi_sigma;
      break;

    case SHAPE_STUDENT_T:
      p->chi_min = fmax (0.0, p->chi_mean - p->n_scale * p->chi_scale);
      p->chi_max = p->chi_mean + p->n_scale * p->chi_scale;
      break;

    case SHAPE_POWER_EXP:
      p->chi_min = p->chi_lower;
      p->chi_max = p->chi_upper;
      break;

    case SHAPE_LENSING:
      p->chi_min              = p->chi_lower;
      p->chi_max              = p->chi_source_upper;
      p->breaks[p->n_break++] = p->chi_source_lower;
      break;

    case SHAPE_MULTI:
      p->chi_min = p->mu[0] - p->n_sigma * p->sg[0];
      p->chi_max = p->mu[0] + p->n_sigma * p->sg[0];

      for (i = 1; i < p->n_bumps; i++)
      {
        p->chi_min = fmin (p->chi_min, p->mu[i] - p->n_sigma * p->sg[i]);
        p->chi_max = fmax (p->chi_max, p->mu[i] + p->n_sigma * p->sg[i]);
      }

      p->chi_min = fmax (0.0, p->chi_min);
      break;

    default:
      fprintf (stderr, "shape_support: bad shape\n");
      exit (1);
  }
}

static Shape
parse_shape (const char *s)
{
  int i;

  for (i = 0; i < SHAPE_LEN; i++)
    if (strcmp (s, shape_names[i]) == 0)
      return (Shape) i;

  fprintf (stderr, "unknown shape '%s'\n", s);
  exit (1);
}

static void
parse_list (const char *s, double *out, int *n)
{
  char *copy = strdup (s);
  char *tok  = strtok (copy, ",");

  *n = 0;

  while ((tok != NULL) && (*n < MAX_BUMPS))
  {
    out[(*n)++] = atof (tok);
    tok         = strtok (NULL, ",");
  }

  free (copy);
}

int
main (int argc, char **argv)
{
  Par p;
  double target = 1.0e-25;
  long ell      = 2;
  double k;
  int i;

  memset (&p, 0, sizeof (Par));
  p.shape   = SHAPE_GAUSS;
  p.n_sigma = 4.0;
  p.n_scale = 6.0;

  for (i = 1; i < argc; i++)
  {
    const char *a = argv[i];

    if (strncmp (a, "--shape=", 8) == 0)
    {
      p.shape = parse_shape (a + 8);
    }
    else if (strncmp (a, "--ell=", 6) == 0)
    {
      ell = atol (a + 6);
    }
    else if (strncmp (a, "--chi-mean=", 11) == 0)
    {
      p.chi_mean = atof (a + 11);
    }
    else if (strncmp (a, "--chi-sigma=", 12) == 0)
    {
      p.chi_sigma = atof (a + 12);
    }
    else if (strncmp (a, "--n-sigma=", 10) == 0)
    {
      p.n_sigma = atof (a + 10);
    }
    else if (strncmp (a, "--chi-lower=", 12) == 0)
    {
      p.chi_lower = atof (a + 12);
    }
    else if (strncmp (a, "--chi-upper=", 12) == 0)
    {
      p.chi_upper = atof (a + 12);
    }
    else if (strncmp (a, "--chi-scale=", 12) == 0)
    {
      p.chi_scale = atof (a + 12);
    }
    else if (strncmp (a, "--nu=", 5) == 0)
    {
      p.nu = atof (a + 5);
    }
    else if (strncmp (a, "--n-scale=", 10) == 0)
    {
      p.n_scale = atof (a + 10);
    }
    else if (strncmp (a, "--alpha=", 8) == 0)
    {
      p.alpha = atof (a + 8);
    }
    else if (strncmp (a, "--beta=", 7) == 0)
    {
      p.beta = atof (a + 7);
    }
    else if (strncmp (a, "--chi-source-lower=", 19) == 0)
    {
      p.chi_source_lower = atof (a + 19);
    }
    else if (strncmp (a, "--chi-source-upper=", 19) == 0)
    {
      p.chi_source_upper = atof (a + 19);
    }
    else if (strncmp (a, "--mu=", 5) == 0)
    {
      parse_list (a + 5, p.mu, &p.n_bumps);
    }
    else if (strncmp (a, "--sigma=", 8) == 0)
    {
      int m;

      parse_list (a + 8, p.sg, &m);
    }
    else if (strncmp (a, "--weight=", 9) == 0)
    {
      int m;

      parse_list (a + 9, p.wt, &m);
    }
    else if (strncmp (a, "--target-rel=", 13) == 0)
    {
      target = atof (a + 13);
    }
    else if (strcmp (a, "--window") == 0)
    {
      p.window_mode = 1;
    }
    else
    {
      fprintf (stderr, "unknown argument '%s'\n", a);
      exit (1);
    }
  }

  p.ell = ell;
  shape_support (&p);

  /* The normalization is measured, not copied from the library. */
  {
    acb_t nrm;
    char *s;

    acb_init (nrm);
    p.with_bessel = 0;
    certified (nrm, &p, target);
    s = arb_get_str (acb_realref (nrm), 30, ARB_STR_NO_RADIUS);
    printf ("# shape=%s ell=%ld support=[%.17g,%.17g] norm=%s\n",
            shape_names[p.shape], p.ell, p.chi_min, p.chi_max, s);
    flint_free (s);

    if (p.window_mode)
    {
      /* W(chi) = W_u(chi) / norm, for chi values on stdin. The window is a
       * closed form, so only the normalization carries a radius. */
      double chi;

      printf ("# shape\tchi\tW\n");

      while (scanf ("%lf", &chi) == 1)
      {
        acb_t x, Wu, Wn;

        acb_init (x);
        acb_init (Wu);
        acb_init (Wn);
        acb_set_d (x, chi);

        if ((chi < p.chi_min) || (chi > p.chi_max))
        {
          acb_zero (Wn);
        }
        else
        {
          window_u (Wu, x, &p, 0, 256);
          acb_div (Wn, Wu, nrm, 256);
        }

        s = arb_get_str (acb_realref (Wn), 30, ARB_STR_NO_RADIUS);
        printf ("%s\t%.17g\t%s\n", shape_names[p.shape], chi, s);
        flint_free (s);
        acb_clear (x);
        acb_clear (Wu);
        acb_clear (Wn);
      }

      acb_clear (nrm);
      flint_cleanup ();

      return 0;
    }

    printf ("# shape\tell\tk\tvalue\tradius\tprec\n");
    p.with_bessel = 1;

    while (scanf ("%lf", &k) == 1)
    {
      acb_t res, val;
      slong prec;

      p.k = k;
      acb_init (res);
      acb_init (val);
      prec = certified (res, &p, target);
      acb_div (val, res, nrm, prec);
      s = arb_get_str (acb_realref (val), 30, ARB_STR_NO_RADIUS);
      printf ("%s\t%ld\t%.17g\t%s\t%.4e\t%ld\n",
              shape_names[p.shape], p.ell, k, s,
              mag_get_d (arb_radref (acb_realref (val))), (long) prec);
      flint_free (s);
      acb_clear (res);
      acb_clear (val);
    }

    acb_clear (nrm);
  }

  flint_cleanup ();

  return 0;
}


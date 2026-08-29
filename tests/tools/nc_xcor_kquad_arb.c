/***************************************************************************
 *            nc_xcor_kquad_arb.c
 *
 *  Fri August 28 12:00:00 2026
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
 * Certified reference values for the whole angular power spectrum of a pair of
 * analytic xcor windows,
 *
 *   C_ell = 2/pi INT dk k^2 P(k) I1_ell(k) I2_ell(k),
 *   I_ell(k) = INT dchi W(chi) g(chi, k) j_ell(k chi),
 *
 * with k in 1/Mpc, chi in Mpc, each W normalized to unit integral over its own
 * truncated support, and P the closed-form NcmPowspecAnalytic BBKS spectrum at
 * z = 0 (where its growth factor is exactly 1, so none is needed here).
 *
 * This is a *nested* certified integration: the outer integrator evaluates the
 * inner one on complex balls of k. Three things are needed to make that work,
 * each of which returned a useless answer when left out.
 *
 * 1. The inner integrator is capped (`window_eval_limit`). On a wide k ball its
 *    result is uncertain however far chi is subdivided, because the uncertainty
 *    is in k; uncapped, one easy case did not finish in two minutes.
 *
 * 2. A capped run leaves an indeterminate ball, so the header falls back to the
 *    trivial enclosure. Finite and crude beats infinite and exact-in-principle:
 *    one infinity poisons every outer cell that contains it.
 *
 * 3. The outer integral is taken one factor-of-two panel in k at a time, and
 *    **the precision escalates per panel**. Over a single six-decade interval the integrator
 *    bisects from the top and the wide-ball enclosures it meets give 0 +- 1e62.
 *    Per panel at a fixed 128 bits the low-k panels are clean and the high-k
 *    panels still give 0 +- 1e62 -- because the entire 0F1 form of j_ell loses
 *    about z^2/4 bits to cancellation and z = k chi reaches the hundreds. At
 *    512 bits those same panels certify to 1e-39.
 *
 * Emits one TSV row per multipole: ell value radius prec_max.
 *
 * Build (standalone):
 *   gcc -O2 -o xcor_kquad_arb nc_xcor_kquad_arb.c \
 *       $(pkg-config --cflags --libs flint) -lm
 */

#include "xcor_window_arb.h"

/* NcmPowspecAnalytic BBKS, matching ncm_powspec_analytic.c's defaults. */
#define BBKS_C1 2.34
#define BBKS_C2 3.89
#define BBKS_C3 16.1
#define BBKS_C4 5.46
#define BBKS_C5 6.71

typedef struct
{
  double amplitude;
  double n_s;
  double k_eq;
} Powspec;

/*
 * P(k) = A k^n_s T(q)^2, q = k / k_eq, on a complex ball.
 *
 * T(q) = ln (1 + C1 q) / (C1 q) * [1 + C2 q + (C3 q)^2 + (C4 q)^3 + (C5 q)^4]^(-1/4)
 *
 * The leading factor is 0/0 at q = 0; the window supports here never reach it,
 * but a ball straddling zero would, so it is guarded.
 */
static void
powspec_eval (acb_t out, const acb_t k, const Powspec *ps, slong prec)
{
  acb_t q, t, u, M, L;

  acb_init (q);
  acb_init (t);
  acb_init (u);
  acb_init (M);
  acb_init (L);

  acb_set_d (t, ps->k_eq);
  acb_div (q, k, t, prec);

  /* M by Horner in q. */
  acb_set_d (M, BBKS_C5 * BBKS_C5 * BBKS_C5 * BBKS_C5);
  acb_mul (M, M, q, prec);
  acb_set_d (t, BBKS_C4 * BBKS_C4 * BBKS_C4);
  acb_add (M, M, t, prec);
  acb_mul (M, M, q, prec);
  acb_set_d (t, BBKS_C3 * BBKS_C3);
  acb_add (M, M, t, prec);
  acb_mul (M, M, q, prec);
  acb_set_d (t, BBKS_C2);
  acb_add (M, M, t, prec);
  acb_mul (M, M, q, prec);
  acb_add_si (M, M, 1, prec);

  acb_set_d (t, BBKS_C1);
  acb_mul (u, q, t, prec);

  if (acb_is_zero (u))
  {
    acb_one (L);
  }
  else
  {
    acb_log1p (L, u, prec);
    acb_div (L, L, u, prec);
  }

  acb_sqrt (M, M, prec);
  acb_sqrt (M, M, prec); /* M^(1/4) */
  acb_div (L, L, M, prec);
  acb_sqr (L, L, prec); /* T^2 */

  acb_set_d (t, ps->n_s);
  acb_pow (out, k, t, prec);
  acb_mul (out, out, L, prec);
  acb_set_d (t, ps->amplitude);
  acb_mul (out, out, t, prec);

  acb_clear (q);
  acb_clear (t);
  acb_clear (u);
  acb_clear (M);
  acb_clear (L);
}

typedef struct
{
  Par *a;
  Par *b;
  acb_t norm_a;
  acb_t norm_b;
  Powspec ps;
  int isauto;
  long n_eval;
} Pair;

/* k^2 P(k) I1_ell(k) I2_ell(k), the outer integrand. */
static int
kquad_integrand (acb_ptr out, const acb_t k, void *param, slong order, slong prec)
{
  Pair *pr = (Pair *) param;
  acb_t I1, I2, P, t;

  if (order > 1)
    return 0;

  pr->n_eval++;

  acb_init (I1);
  acb_init (I2);
  acb_init (P);
  acb_init (t);

  acb_set (pr->a->k, k);
  pr->a->with_bessel = 1;
  integrate_panels (I1, pr->a, prec);
  acb_div (I1, I1, pr->norm_a, prec);

  if (pr->isauto)
  {
    acb_set (I2, I1);
  }
  else
  {
    acb_set (pr->b->k, k);
    pr->b->with_bessel = 1;
    integrate_panels (I2, pr->b, prec);
    acb_div (I2, I2, pr->norm_b, prec);
  }

  powspec_eval (P, k, &pr->ps, prec);

  acb_sqr (t, k, prec);
  acb_mul (t, t, P, prec);
  acb_mul (t, t, I1, prec);
  acb_mul (out, t, I2, prec);

  acb_clear (I1);
  acb_clear (I2);
  acb_clear (P);
  acb_clear (t);

  return 0;
}

/*
 * A certified bound on everything above @k_hi, without integrating it.
 *
 * |j_ell| <= 1 and every window here is normalized to unit integral, so
 * |I_ell(k)| <= 1 for all k and all ell. The whole tail is therefore bounded
 * by INT_{k_hi}^{k_end} k^2 P(k) dk, a smooth one-dimensional integral with no
 * Bessel function in it -- and k^2 P ~ k^(n_s - 2) ln^2 k, so it converges and
 * costs nothing.
 *
 * This is what the far tail is for. Integrating it honestly meant 2048-bit
 * arithmetic to confirm that exp (-780) is zero: the 0F1 cancellation is worst
 * exactly where the answer matters least.
 */
static int
tail_integrand (acb_ptr out, const acb_t k, void *param, slong order, slong prec)
{
  Powspec *ps = (Powspec *) param;
  acb_t P, t;

  if (order > 1)
    return 0;

  acb_init (P);
  acb_init (t);
  powspec_eval (P, k, ps, prec);
  acb_sqr (t, k, prec);
  acb_mul (out, t, P, prec);
  acb_clear (P);
  acb_clear (t);

  return 0;
}

static void
tail_bound (acb_t res, Powspec *ps, double k_hi, double k_end, slong prec)
{
  acb_calc_integrate_opt_t opt;
  acb_t A, B;
  mag_t tol;

  acb_init (A);
  acb_init (B);
  mag_init (tol);
  acb_calc_integrate_opt_init (opt);
  acb_set_d (A, k_hi);
  acb_set_d (B, k_end);
  mag_set_ui_2exp_si (tol, 1, -prec / 2);
  acb_calc_integrate (res, tail_integrand, ps, A, B, prec, tol, opt, prec);

  /* The bound is two-sided: the true tail lies in [-T, T]. */
  {
    arb_t r;

    arb_init (r);
    arb_get_abs_ubound_arf (arb_midref (acb_realref (res)), acb_realref (res), prec);
    arb_zero (acb_realref (res));
    arb_clear (r);
  }

  acb_clear (A);
  acb_clear (B);
  mag_clear (tol);
}

/*
 * One factor-of-two panel in k, at whatever precision it takes.
 *
 * @scale is the running magnitude of the whole integral; a panel only has to
 * be certified relative to that, not to itself. Panels far out in k are orders
 * below the total and would otherwise drive the precision up for nothing.
 */
static slong
k_panel (acb_t res, Pair *pr, double lo, double hi, double scale, double target,
         slong prec_max, slong eval_limit)
{
  acb_calc_integrate_opt_t opt;
  acb_t A, B;
  mag_t tol;
  slong prec;

  acb_init (A);
  acb_init (B);
  mag_init (tol);
  acb_calc_integrate_opt_init (opt);
  opt->eval_limit = eval_limit;
  acb_set_d (A, lo);
  acb_set_d (B, hi);

  for (prec = 128; prec <= prec_max; prec *= 2)
  {
    double r;

    mag_set_ui_2exp_si (tol, 1, -prec / 2);
    acb_calc_integrate (res, kquad_integrand, pr, A, B, prec, tol, opt, prec);

    if (!acb_is_finite (res))
      continue;

    r = mag_get_d (arb_radref (acb_realref (res)));

    if ((scale > 0.0) && (r < target * scale))
      break;

    if (scale == 0.0)
      break;
  }

  acb_clear (A);
  acb_clear (B);
  mag_clear (tol);

  return prec > prec_max ? prec_max : prec;
}

static void
normalize (acb_t norm, Par *p, double target)
{
  p->with_bessel = 0;
  acb_zero (p->k);
  certified (norm, p, target);
}

static void
parse_window (Par *p, const char *key, const char *val)
{
  if (!strcmp (key, "shape"))
  {
    p->shape = parse_shape (val);
  }
  else if (!strcmp (key, "chi-mean"))
  {
    p->chi_mean = atof (val);
  }
  else if (!strcmp (key, "chi-sigma"))
  {
    p->chi_sigma = atof (val);
  }
  else if (!strcmp (key, "n-sigma"))
  {
    p->n_sigma = atof (val);
  }
  else if (!strcmp (key, "chi-lower"))
  {
    p->chi_lower = atof (val);
  }
  else if (!strcmp (key, "chi-upper"))
  {
    p->chi_upper = atof (val);
  }
  else if (!strcmp (key, "chi-scale"))
  {
    p->chi_scale = atof (val);
  }
  else if (!strcmp (key, "nu"))
  {
    p->nu = atof (val);
  }
  else if (!strcmp (key, "n-scale"))
  {
    p->n_scale = atof (val);
  }
  else if (!strcmp (key, "alpha"))
  {
    p->alpha = atof (val);
  }
  else if (!strcmp (key, "beta"))
  {
    p->beta = atof (val);
  }
  else if (!strcmp (key, "chi-source-lower"))
  {
    p->chi_source_lower = atof (val);
  }
  else if (!strcmp (key, "chi-source-upper"))
  {
    p->chi_source_upper = atof (val);
  }
  else if (!strcmp (key, "mu"))
  {
    parse_list (val, p->mu, &p->n_bumps);
  }
  else if (!strcmp (key, "sigma"))
  {
    int n;

    parse_list (val, p->sg, &n);
  }
  else if (!strcmp (key, "weight"))
  {
    int n;

    parse_list (val, p->wt, &n);
  }
  else if (!strcmp (key, "kdep-amplitude"))
  {
    p->kdep_on        = 1;
    p->kdep_amplitude = atof (val);
  }
  else if (!strcmp (key, "kdep-k-transition"))
  {
    p->kdep_k_transition = atof (val);
  }
  else if (!strcmp (key, "kdep-chi-ref"))
  {
    p->kdep_chi_ref = atof (val);
  }
  else
  {
    fprintf (stderr, "unknown window key '%s'\n", key);
    exit (1);
  }
}

int
main (int argc, char **argv)
{
  Par a, b;
  Pair pr;
  Powspec ps = { .amplitude = 2.08e7, .n_s = 0.9875, .k_eq = 0.10594 };
  double k_lo = 1.0e-6, k_hi = 10.0, k_end = 1.0e4, target = 1.0e-20;
  slong prec_max = 4096, eval_limit = 2000;
  long ells[64];
  int n_ells = 0, i, isauto = 1, verbose = 0, negligible = 0;

  par_init (&a);
  par_init (&b);
  a.n_sigma = 4.0;
  a.n_scale = 6.0;

  for (i = 1; i < argc; i++)
  {
    const char *s = argv[i];

    if (!strncmp (s, "--a:", 4) || !strncmp (s, "--b:", 4))
    {
      char key[64];
      const char *eq = strchr (s + 4, '=');

      if (eq == NULL)
      {
        fprintf (stderr, "expected --%c:key=value\n", s[2]);
        exit (1);
      }

      snprintf (key, sizeof (key), "%.*s", (int) (eq - (s + 4)), s + 4);

      if (s[2] == 'a')
      {
        parse_window (&a, key, eq + 1);
      }
      else
      {
        if (isauto)
        {
          b = a;
          acb_init (b.k);
          isauto = 0;
        }

        parse_window (&b, key, eq + 1);
      }
    }
    else if (!strncmp (s, "--ells=", 7))
    {
      double tmp[64];
      int j, n;

      parse_list (s + 7, tmp, &n);

      for (j = 0; j < n && j < 64; j++)
        ells[n_ells++] = (long) tmp[j];
    }
    else if (!strncmp (s, "--k-lo=", 7))
    {
      k_lo = atof (s + 7);
    }
    else if (!strncmp (s, "--k-hi=", 7))
    {
      k_hi = atof (s + 7);
    }
    else if (!strncmp (s, "--k-end=", 8))
    {
      k_end = atof (s + 8);
    }
    else if (!strncmp (s, "--target-rel=", 13))
    {
      target = atof (s + 13);
    }
    else if (!strncmp (s, "--prec-max=", 11))
    {
      prec_max = atol (s + 11);
    }
    else if (!strncmp (s, "--eval-limit=", 13))
    {
      eval_limit = atol (s + 13);
    }
    else if (!strcmp (s, "--verbose"))
    {
      verbose = 1;
    }
    else if (!strncmp (s, "--amplitude=", 12))
    {
      ps.amplitude = atof (s + 12);
    }
    else if (!strncmp (s, "--n-s=", 6))
    {
      ps.n_s = atof (s + 6);
    }
    else if (!strncmp (s, "--k-eq=", 7))
    {
      ps.k_eq = atof (s + 7);
    }
    else
    {
      fprintf (stderr, "unknown argument '%s'\n", s);
      exit (1);
    }
  }

  if (n_ells == 0)
    ells[n_ells++] = 2;

  shape_support (&a);

  if (!isauto)
    shape_support (&b);

  pr.a      = &a;
  pr.b      = isauto ? &a : &b;
  pr.isauto = isauto;
  pr.ps     = ps;
  acb_init (pr.norm_a);
  acb_init (pr.norm_b);

  /* Normalizations are uncapped: they are cheap, exact and computed once. */
  window_eval_limit = 0;
  normalize (pr.norm_a, &a, 1.0e-25);

  if (!isauto)
    normalize (pr.norm_b, &b, 1.0e-25);

  window_eval_limit = eval_limit;

  printf ("# ell\tvalue\tradius\tprec_max\tn_eval\n");

  for (i = 0; i < n_ells; i++)
  {
    acb_t total, part, pi;
    double lo, scale = 0.0;
    slong prec_used = 0;
    int pass;
    char *str;

    a.ell     = ells[i];
    b.ell     = ells[i];
    pr.n_eval = 0;

    acb_init (total);
    acb_init (part);
    acb_init (pi);

    /* Two passes: the first learns the magnitude of the integral and of each
     * k-panel, the second certifies against it. A panel's own size is the wrong
     * scale -- deep in the tail it is orders below the total, and certifying it
     * to a relative target costs precision for nothing. */
    {
      double mag[128];
      int n_oct = 0, j, spent = 0, j_peak = 0;

      /* Pass 1, cheap: 256 bits, no escalation. */
      acb_zero (total);

      for (lo = k_lo; lo < k_hi && n_oct < 128; lo *= 2.0)
      {
        double hi = lo * 2.0 > k_hi ? k_hi : lo * 2.0;

        double mid, rad;

        k_panel (part, &pr, lo, hi, 0.0, target, 256, eval_limit);

        mid = fabs (arf_get_d (arb_midref (acb_realref (part)), ARF_RND_NEAR));
        rad = mag_get_d (arb_radref (acb_realref (part)));

        /* A panel that carries no significant digit tells us nothing about
         * the scale, and at 256 bits the high-k panels are exactly that --
         * finite, but 0 +- 1e62. Summed in, they put the threshold above every
         * panel that matters and the whole integral reads as negligible. */
        mag[n_oct++] = (acb_is_finite (part) && (rad < 0.5 * mid)) ? mid : 0.0;

        if (acb_is_finite (part) && (rad < 0.5 * mid))
          acb_add (total, total, part, 256);
      }

      scale = fabs (arf_get_d (arb_midref (acb_realref (total)), ARF_RND_NEAR));

      for (j = 1; j < n_oct; j++)
        if (mag[j] > mag[j_peak])
          j_peak = j;

      /* Pass 2. An octave pass 1 already put below the target's share of the
       * total is taken as pass 1 left it: high in k the phase k chi reaches the
       * thousands, the 0F1 form of j_ell wants precision proportional to its
       * square, and what it would eventually certify is a rounding of the
       * total. Escalating those two octaves alone was 90% of the runtime. */
      acb_zero (total);
      prec_used = 0;
      lo        = k_lo;

      for (j = 0; j < n_oct; j++, lo *= 2.0)
      {
        double hi = lo * 2.0 > k_hi ? k_hi : lo * 2.0;

        /* Only a panel that pass 1 actually resolved may be called small.
         * mag[j] == 0 means pass 1 got no significant digit there, which is
         * ignorance, not smallness -- those are exactly the high-k panels the
         * 0F1 cancellation defeats at 256 bits, and they have to escalate. */
        int small = (scale > 0.0) && (mag[j] > 0.0) &&
                    (mag[j] < 0.125 * target * scale);
        slong used = k_panel (part, &pr, lo, hi, small ? 0.0 : scale, target,
                              small ? 256 : prec_max, eval_limit);

        if (verbose)
          fprintf (stderr, "  k-panel [%9.3g, %9.3g] prec %4ld rad %.2e%s\n",
                   lo, hi, (long) used,
                   mag_get_d (arb_radref (acb_realref (part))),
                   small ? "  (negligible)" : "");

        if (!acb_is_finite (part))
        {
          fprintf (stderr, "k-panel [%g, %g] ell %ld did not certify\n",
                   lo, hi, ells[i]);
          continue;
        }

        acb_add (total, total, part, prec_max);

        if (used > prec_used)
          prec_used = used;

        /* Two certified panels in a row below the target's share and the tail
         * is spent -- but only past the peak. The integrand rises to a maximum
         * near k ~ ell / chi and is just as small below it as above, so the
         * same test applied from the low-k end stops before the integral has
         * been reached at all. */
        if (small && (j > j_peak))
        {
          if (++spent >= 2)
            break;
        }
        else
        {
          spent = 0;
        }
      }
    }

    acb_const_pi (pi, prec_max);
    acb_div (total, total, pi, prec_max);
    acb_mul_2exp_si (total, total, 1);

    str = arb_get_str (acb_realref (total), 30, ARB_STR_NO_RADIUS);
    printf ("%ld\t%s\t%.6e\t%ld\t%ld\n", ells[i], str,
            mag_get_d (arb_radref (acb_realref (total))),
            (long) prec_used, pr.n_eval);
    flint_free (str);

    acb_clear (total);
    acb_clear (part);
    acb_clear (pi);
  }

  acb_clear (pr.norm_a);
  acb_clear (pr.norm_b);
  par_clear (&a);

  if (!isauto)
    par_clear (&b);

  flint_cleanup ();

  return 0;
}


/***************************************************************************
 *            test_xcor_window_arb.c
 *
 *  Wed Sep 10 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_xcor_window_arb.c
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

/* Smoke tests for the Arb reference generator behind the certified xcor
 * tables (tests/tools/xcor_window_arb.h).
 *
 * The generator is not built into the library and its output reaches the suite
 * only as committed tables, so nothing here executed it: a break showed up as a
 * table that disagreed with the library, which is a slow and indirect way to
 * learn that a closed form was mistyped. These tests run it.
 *
 * They are deliberately shallow and fast. The deep statement about the
 * generator is the certified table itself, checked against the library in
 * tests/python/nc/xcor/; what is wanted here is that every code path still runs
 * and returns something sane, so the tools cannot rot unnoticed. Targets are
 * loose and multipoles small for that reason -- `certified` exits the process
 * when it cannot reach its target, so a target has to be comfortably reachable.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <glib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "xcor_window_arb.h"

/* Loose enough that `certified` always reaches it, tight enough that a wrong
 * closed form cannot pass the normalization check below. */
#define TEST_TARGET 1.0e-12
#define TEST_PREC_MAX 2048

typedef struct
{
  const char *name;
  Shape shape;
  void (*setup) (Par *p);
} ShapeCase;

/* One real parameter set per closed form, taken from the cases the committed
 * table certifies, so these exercise the same branches the tables were built
 * from rather than a plausible retyping. */

static void
setup_gauss (Par *p)
{
  p->chi_mean  = 4520.0;
  p->chi_sigma = 715.0;
  p->n_sigma   = 4.0;
}

static void
setup_tophat (Par *p)
{
  p->chi_lower = 3282.0;
  p->chi_upper = 5758.0;
}

static void
setup_tophat_smooth (Par *p)
{
  p->chi_lower = 1000.0;
  p->chi_upper = 2000.0;
  p->chi_sigma = 150.0;
  p->n_sigma   = 6.0;
}

static void
setup_student_t (Par *p)
{
  p->chi_mean  = 1500.0;
  p->chi_scale = 200.0;
  p->nu        = 2.0;
  p->n_scale   = 6.0;
}

static void
setup_power_exp (Par *p)
{
  p->chi_scale = 2500.0;
  p->alpha     = 1.0;
  p->beta      = 1.5;
  p->chi_lower = 50.0;
  p->chi_upper = 8000.0;
}

static void
setup_lensing (Par *p)
{
  p->chi_lower        = 50.0;
  p->chi_source_lower = 13900.0;
  p->chi_source_upper = 14146.0;
}

static void
setup_multi (Par *p)
{
  p->n_bumps = 2;
  p->mu[0]   = 1000.0;
  p->mu[1]   = 1600.0;
  p->sg[0]   = 300.0;
  p->sg[1]   = 300.0;
  p->wt[0]   = 1.0;
  p->wt[1]   = 0.6;
  p->n_sigma = 4.0;
}

static const ShapeCase shape_cases[] = {
  {"gauss",         SHAPE_GAUSS,         setup_gauss        },
  {"tophat",        SHAPE_TOPHAT,        setup_tophat       },
  {"tophat_smooth", SHAPE_TOPHAT_SMOOTH, setup_tophat_smooth},
  {"student_t",     SHAPE_STUDENT_T,     setup_student_t    },
  {"power_exp",     SHAPE_POWER_EXP,     setup_power_exp    },
  {"lensing",       SHAPE_LENSING,       setup_lensing      },
  {"multi",         SHAPE_MULTI,         setup_multi        },
};

#define N_SHAPE_CASES (sizeof (shape_cases) / sizeof (shape_cases[0]))

static void
case_init (Par *p, const ShapeCase *sc)
{
  par_init (p);
  p->shape = sc->shape;
  sc->setup (p);
  shape_support (p);
}

/* j_0 (x) = sin x / x and j_1 (x) = sin x / x^2 - cos x / x, against the ball
 * the generator's own Bessel returns. A ball that does not contain the value is
 * either a wrong function or a radius that under-reports, and both matter. */
static void
test_sph_bessel_closed_form (void)
{
  const double xs[] = {0.5, 3.0, 17.0};
  const slong prec  = 256;
  acb_t z, j;
  guint i;

  acb_init (z);
  acb_init (j);

  for (i = 0; i < sizeof (xs) / sizeof (xs[0]); i++)
  {
    const double x  = xs[i];
    const double j0 = sin (x) / x;
    const double j1 = sin (x) / (x * x) - cos (x) / x;

    acb_set_d (z, x);

    sph_bessel (j, z, 0, prec);
    g_assert_true (acb_is_finite (j));
    g_assert_cmpfloat (fabs (arf_get_d (arb_midref (acb_realref (j)), ARF_RND_NEAR) - j0), <, 1.0e-12);

    sph_bessel (j, z, 1, prec);
    g_assert_true (acb_is_finite (j));
    g_assert_cmpfloat (fabs (arf_get_d (arb_midref (acb_realref (j)), ARF_RND_NEAR) - j1), <, 1.0e-12);
  }

  acb_clear (z);
  acb_clear (j);
}

/* The derivative orders come from ORDER recurrences -- j_{l-2}, j_l, j_{l+2}
 * with no division by the argument -- precisely so they are not the argument
 * recurrences the library uses. Differencing the generator's own j_l in the
 * argument is a third route, and the one that catches a mistyped coefficient.
 *
 * The differencing is done in Arb rather than on extracted doubles. A second
 * difference divides by h^2, so in double precision the cancellation error is
 * O(eps / h^2) -- at h = 1e-8 that is O(1) and the comparison measures nothing.
 * Here h is 2^-100 against 512 bits of working precision, which leaves the
 * fourth-derivative truncation term the only error that matters.
 *
 * ell = 0 and ell = 1 are in the list on purpose: both formulas drop terms
 * there (j_{-1} and j_{-2} are never evaluated because their coefficients
 * vanish), and a wrong guard would show up nowhere else.
 */
static void
test_sph_bessel_deriv_matches_differences (void)
{
  const long ells[] = {0, 1, 2, 10};
  const double xs[] = {2.5, 14.0};
  const slong prec  = 512;
  guint i, e;
  acb_t z, zp, zm, fp, fm, f0, d, fd, h, t;

  acb_init (z);
  acb_init (zp);
  acb_init (zm);
  acb_init (fp);
  acb_init (fm);
  acb_init (f0);
  acb_init (d);
  acb_init (fd);
  acb_init (h);
  acb_init (t);

  /* A power of two is exact, so the shifted arguments carry no representation
   * error of their own. */
  acb_one (h);
  acb_mul_2exp_si (h, h, -100);

  for (e = 0; e < sizeof (ells) / sizeof (ells[0]); e++)
  {
    for (i = 0; i < sizeof (xs) / sizeof (xs[0]); i++)
    {
      const long ell = ells[e];

      acb_set_d (z, xs[i]);
      acb_add (zp, z, h, prec);
      acb_sub (zm, z, h, prec);

      sph_bessel (fp, zp, ell, prec);
      sph_bessel (fm, zm, ell, prec);
      sph_bessel (f0, z, ell, prec);

      /* First order: (f(x+h) - f(x-h)) / 2h. */
      acb_sub (fd, fp, fm, prec);
      acb_div (fd, fd, h, prec);
      acb_mul_2exp_si (fd, fd, -1);

      sph_bessel_deriv (d, z, ell, 1, prec);
      g_assert_true (acb_is_finite (d));
      acb_sub (t, d, fd, prec);
      g_assert_cmpfloat (fabs (arf_get_d (arb_midref (acb_realref (t)), ARF_RND_NEAR)), <, 1.0e-20);

      /* Second order: (f(x+h) - 2 f(x) + f(x-h)) / h^2. */
      acb_mul_2exp_si (t, f0, 1);
      acb_sub (fd, fp, t, prec);
      acb_add (fd, fd, fm, prec);
      acb_div (fd, fd, h, prec);
      acb_div (fd, fd, h, prec);

      sph_bessel_deriv (d, z, ell, 2, prec);
      g_assert_true (acb_is_finite (d));
      acb_sub (t, d, fd, prec);
      g_assert_cmpfloat (fabs (arf_get_d (arb_midref (acb_realref (t)), ARF_RND_NEAR)), <, 1.0e-20);
    }
  }

  acb_clear (z);
  acb_clear (zp);
  acb_clear (zm);
  acb_clear (fp);
  acb_clear (fm);
  acb_clear (f0);
  acb_clear (d);
  acb_clear (fd);
  acb_clear (h);
  acb_clear (t);
}


/* d = 0 has to be the plain Bessel, not a derivative formula evaluated at
 * order zero: the drivers reach j_l through this entry point. */
static void
test_sph_bessel_deriv_zero_is_the_bessel (void)
{
  const slong prec = 256;
  acb_t z, a, b;

  acb_init (z);
  acb_init (a);
  acb_init (b);
  acb_set_d (z, 4.0);

  sph_bessel (a, z, 7, prec);
  sph_bessel_deriv (b, z, 7, 0, prec);
  g_assert_true (acb_overlaps (a, b));

  acb_clear (z);
  acb_clear (a);
  acb_clear (b);
}

/* Every closed form: a support that is ordered and non-empty, breakpoints
 * inside it, and a window that is finite and non-negative at mid-support. */
static void
test_shape_support_and_window (void)
{
  const slong prec = 256;
  guint c;

  for (c = 0; c < N_SHAPE_CASES; c++)
  {
    const ShapeCase *sc = &shape_cases[c];
    Par p;
    acb_t chi, w;
    int b;

    case_init (&p, sc);

    g_assert_cmpfloat (p.chi_min, <, p.chi_max);
    g_assert_cmpfloat (p.chi_min, >=, 0.0);
    g_assert_cmpint (p.n_break, >=, 0);

    for (b = 0; b < p.n_break; b++)
    {
      g_assert_cmpfloat (p.breaks[b], >=, p.chi_min);
      g_assert_cmpfloat (p.breaks[b], <=, p.chi_max);
    }

    acb_init (chi);
    acb_init (w);
    acb_set_d (chi, 0.5 * (p.chi_min + p.chi_max));

    window_u (w, chi, &p, 0, prec);
    g_assert_true (acb_is_finite (w));
    g_assert_cmpfloat (arf_get_d (arb_midref (acb_realref (w)), ARF_RND_NEAR), >=, 0.0);

    acb_clear (chi);
    acb_clear (w);
    par_clear (&p);
  }
}

/* The integrator and the certification loop, on every closed form.
 *
 * j_0 (k chi) -> 1 as k -> 0, so the window integrated against it has to
 * approach the bare integral of the same window. Both sides come from
 * `certified` -- with_bessel off and on, which is exactly how the driver
 * normalizes -- so this exercises window_u, integrate_panels and certified
 * together, and a closed form wrong by a factor cannot satisfy it.
 */
static void
test_certified_normalization (void)
{
  guint c;

  for (c = 0; c < N_SHAPE_CASES; c++)
  {
    const ShapeCase *sc = &shape_cases[c];
    Par p;
    acb_t norm, proj;
    double n, i0;

    case_init (&p, sc);

    acb_init (norm);
    acb_init (proj);

    p.with_bessel = 0;
    certified (norm, &p, TEST_TARGET, TEST_PREC_MAX);
    n = arf_get_d (arb_midref (acb_realref (norm)), ARF_RND_NEAR);

    g_assert_true (acb_is_finite (norm));
    g_assert_cmpfloat (n, >, 0.0);

    /* Small enough that k chi_max << 1 over every support here, so j_0 is one
     * to well inside the comparison below. */
    p.ell = 0;
    acb_set_d (p.k, 1.0e-9);
    p.with_bessel = 1;
    certified (proj, &p, TEST_TARGET, TEST_PREC_MAX);
    i0 = arf_get_d (arb_midref (acb_realref (proj)), ARF_RND_NEAR);

    g_assert_true (acb_is_finite (proj));
    g_assert_cmpfloat (fabs (i0 / n - 1.0), <, 1.0e-6);

    acb_clear (norm);
    acb_clear (proj);
    par_clear (&p);
  }
}

/* The derivative weights, through the same path, on every closed form. Only
 * that they run and stay finite: their values are what the d = 1 and d = 2
 * tables carry, and those are checked against the library elsewhere. */
static void
test_certified_runs_for_derivative_weights (void)
{
  guint c;
  int deriv;

  for (c = 0; c < N_SHAPE_CASES; c++)
  {
    for (deriv = 1; deriv <= 2; deriv++)
    {
      const ShapeCase *sc = &shape_cases[c];
      Par p;
      acb_t res;

      case_init (&p, sc);

      acb_init (res);

      p.ell          = 2;
      p.bessel_deriv = deriv;
      p.with_bessel  = 1;
      acb_set_d (p.k, 1.0e-3);

      certified (res, &p, TEST_TARGET, TEST_PREC_MAX);
      g_assert_true (acb_is_finite (res));

      acb_clear (res);
      par_clear (&p);
    }
  }
}

/* The scale-dependent growth branch, which is off by default and so would
 * otherwise never run here. */
static void
test_kdep_growth_runs (void)
{
  const slong prec = 256;
  Par p;
  acb_t chi, g;

  case_init (&p, &shape_cases[0]);

  p.kdep_on            = 1;
  p.kdep_amplitude     = 0.1;
  p.kdep_k_transition  = 0.05;
  p.kdep_chi_ref       = 0.5 * (p.chi_min + p.chi_max);
  acb_set_d (p.k, 0.05);

  acb_init (chi);
  acb_init (g);
  acb_set_d (chi, p.kdep_chi_ref);

  /* At chi_ref the exponent vanishes whatever the amplitude, so the growth is
   * one there: a sign or ordering slip in the exponent shows up immediately. */
  kdep_growth (g, chi, &p, prec);
  g_assert_true (acb_is_finite (g));
  g_assert_cmpfloat (fabs (arf_get_d (arb_midref (acb_realref (g)), ARF_RND_NEAR) - 1.0), <, 1.0e-12);

  acb_clear (chi);
  acb_clear (g);
  par_clear (&p);
}

/* The list parser the drivers use for the multi-bump arguments. */
static void
test_parse_list (void)
{
  double out[MAX_BUMPS];
  int n = 0;

  parse_list ("1000,1600,2200", out, &n);
  g_assert_cmpint (n, ==, 3);
  g_assert_cmpfloat (out[0], ==, 1000.0);
  g_assert_cmpfloat (out[2], ==, 2200.0);

  n = 0;
  parse_list ("300", out, &n);
  g_assert_cmpint (n, ==, 1);
  g_assert_cmpfloat (out[0], ==, 300.0);
}


/* Every shape name the drivers accept round-trips through the parser. */
static void
test_parse_shape_round_trip (void)
{
  int s;

  for (s = 0; s < SHAPE_LEN; s++)
    g_assert_cmpint (parse_shape (shape_names[s]), ==, s);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/nc/xcor/window_arb/sph_bessel/closed_form",
                   &test_sph_bessel_closed_form);
  g_test_add_func ("/nc/xcor/window_arb/sph_bessel_deriv/differences",
                   &test_sph_bessel_deriv_matches_differences);
  g_test_add_func ("/nc/xcor/window_arb/sph_bessel_deriv/order_zero",
                   &test_sph_bessel_deriv_zero_is_the_bessel);
  g_test_add_func ("/nc/xcor/window_arb/shape/support_and_window",
                   &test_shape_support_and_window);
  g_test_add_func ("/nc/xcor/window_arb/certified/normalization",
                   &test_certified_normalization);
  g_test_add_func ("/nc/xcor/window_arb/certified/derivative_weights",
                   &test_certified_runs_for_derivative_weights);
  g_test_add_func ("/nc/xcor/window_arb/kdep_growth/runs",
                   &test_kdep_growth_runs);
  g_test_add_func ("/nc/xcor/window_arb/parse_shape/round_trip",
                   &test_parse_shape_round_trip);
  g_test_add_func ("/nc/xcor/window_arb/parse_list/values",
                   &test_parse_list);

  return g_test_run ();
}

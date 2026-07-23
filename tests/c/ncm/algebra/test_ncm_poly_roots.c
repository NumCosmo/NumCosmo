/***************************************************************************
 *            test_ncm_poly_roots.c
 *
 *  Mon Jul 6 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

/* Independent tests for NcmPolyRoots, deliberately with no reference to the
 * weak-lensing problem it was extracted from -- this module is pure
 * polynomial-root-finding, cross-checked here against gsl_poly_complex_solve
 * (trustworthy for correctness; it was only ever avoided in the hot loop
 * this was extracted from for performance reasons, not correctness). */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_poly.h>

static void
_sort (gdouble *x, gint n)
{
  gint i, j;

  for (i = 1; i < n; i++)
  {
    const gdouble key = x[i];

    for (j = i - 1; (j >= 0) && (x[j] > key); j--)
      x[j + 1] = x[j];

    x[j + 1] = key;
  }
}

/* Ground truth via GSL's general (eigenvalue-based) solver, ascending
 * coefficients, degree = n_coeffs-1, leading coefficient assumed nonzero. */
static gint
_gsl_real_roots (const gdouble *a, gint n_coeffs, gdouble *roots)
{
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc (n_coeffs);
  gdouble z[8];
  gint status = gsl_poly_complex_solve (a, n_coeffs, w, z);
  gint n      = 0;
  gint k;

  gsl_poly_complex_workspace_free (w);
  g_assert_cmpint (status, ==, GSL_SUCCESS);

  for (k = 0; k < n_coeffs - 1; k++)
  {
    if (fabs (z[2 * k + 1]) < 1.0e-7 * (fabs (z[2 * k]) + 1.0))
      roots[n++] = z[2 * k];
  }

  _sort (roots, n);

  return n;
}

static void
test_known_quadratic (void)
{
  /* (x-1)(x-2) = x^2 - 3x + 2 */
  const gdouble a[3] = {2.0, -3.0, 1.0};
  gdouble roots[2];
  gint n = ncm_poly_roots_real_quadratic (a, roots);

  g_assert_cmpint (n, ==, 2);
  _sort (roots, n);
  g_assert_cmpfloat (fabs (roots[0] - 1.0), <, 1.0e-12);
  g_assert_cmpfloat (fabs (roots[1] - 2.0), <, 1.0e-12);
}

static void
test_known_quadratic_no_real_roots (void)
{
  /* x^2 + 1 = 0 */
  const gdouble a[3] = {1.0, 0.0, 1.0};
  gdouble roots[2];
  gint n = ncm_poly_roots_real_quadratic (a, roots);

  g_assert_cmpint (n, ==, 0);
}

static void
test_known_cubic (void)
{
  /* (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 */
  const gdouble a[4] = {-6.0, 11.0, -6.0, 1.0};
  gdouble roots[3];
  gint n = ncm_poly_roots_real_cubic (a, roots);

  g_assert_cmpint (n, ==, 3);
  _sort (roots, n);
  g_assert_cmpfloat (fabs (roots[0] - 1.0), <, 1.0e-11);
  g_assert_cmpfloat (fabs (roots[1] - 2.0), <, 1.0e-11);
  g_assert_cmpfloat (fabs (roots[2] - 3.0), <, 1.0e-11);
}

static void
test_known_quartic_all_real (void)
{
  /* (x+2)(x-1)(x-2)(x-5) = x^4-6x^3+x^2+24x-20 */
  const gdouble a[5] = {-20.0, 24.0, 1.0, -6.0, 1.0};
  gdouble roots[4];
  gint n = ncm_poly_roots_real_quartic (a, roots);

  g_assert_cmpint (n, ==, 4);
  _sort (roots, n);
  g_assert_cmpfloat (fabs (roots[0] - (-2.0)), <, 1.0e-9);
  g_assert_cmpfloat (fabs (roots[1] - 1.0), <, 1.0e-9);
  g_assert_cmpfloat (fabs (roots[2] - 2.0), <, 1.0e-9);
  g_assert_cmpfloat (fabs (roots[3] - 5.0), <, 1.0e-9);
}

static void
test_known_quartic_two_real_two_complex (void)
{
  /* (x-1)(x-2)(x^2+1) = x^4-3x^3+3x^2-3x+2 */
  const gdouble a[5] = {2.0, -3.0, 3.0, -3.0, 1.0};
  gdouble roots[4];
  gint n = ncm_poly_roots_real_quartic (a, roots);

  g_assert_cmpint (n, ==, 2);
  _sort (roots, n);
  g_assert_cmpfloat (fabs (roots[0] - 1.0), <, 1.0e-9);
  g_assert_cmpfloat (fabs (roots[1] - 2.0), <, 1.0e-9);
}

static void
test_quartic_or_lower_degeneracy (void)
{
  gdouble roots[4];
  gint n;

  /* a4 = 0: genuinely a cubic, (x-1)(x-2)(x-3) */
  {
    const gdouble a[5] = {-6.0, 11.0, -6.0, 1.0, 0.0};

    n = ncm_poly_roots_real_quartic_or_lower (a, roots);
    g_assert_cmpint (n, ==, 3);
    _sort (roots, n);
    g_assert_cmpfloat (fabs (roots[0] - 1.0), <, 1.0e-11);
    g_assert_cmpfloat (fabs (roots[1] - 2.0), <, 1.0e-11);
    g_assert_cmpfloat (fabs (roots[2] - 3.0), <, 1.0e-11);
  }

  /* a4 = a3 = 0: genuinely a quadratic, (x-1)(x-2) */
  {
    const gdouble a[5] = {2.0, -3.0, 1.0, 0.0, 0.0};

    n = ncm_poly_roots_real_quartic_or_lower (a, roots);
    g_assert_cmpint (n, ==, 2);
    _sort (roots, n);
    g_assert_cmpfloat (fabs (roots[0] - 1.0), <, 1.0e-12);
    g_assert_cmpfloat (fabs (roots[1] - 2.0), <, 1.0e-12);
  }

  /* a4 = a3 = a2 = 0: genuinely linear, x - 5 */
  {
    const gdouble a[5] = {-5.0, 1.0, 0.0, 0.0, 0.0};

    n = ncm_poly_roots_real_quartic_or_lower (a, roots);
    g_assert_cmpint (n, ==, 1);
    g_assert_cmpfloat (fabs (roots[0] - 5.0), <, 1.0e-12);
  }
}

/* Random-trial cross-checks against GSL's general solver: root count and
 * values must agree, across many random cubics/quartics with genuinely
 * random (not hand-picked) coefficients. */
static void
test_random_cubic_matches_gsl (void)
{
  GRand *rng = g_rand_new_with_seed (20260706);
  guint trial;

  for (trial = 0; trial < 2000; trial++)
  {
    gdouble a[4];
    gdouble mine[3], ref[3];
    gint n_mine, n_ref, k;
    gint i;

    for (i = 0; i < 4; i++)
      a[i] = g_rand_double_range (rng, -3.0, 3.0);

    if (fabs (a[3]) < 1.0e-2)
      continue;  /* keep it a genuine cubic */

    n_mine = ncm_poly_roots_real_cubic (a, mine);
    n_ref  = _gsl_real_roots (a, 4, ref);

    g_assert_cmpint (n_mine, ==, n_ref);
    _sort (mine, n_mine);

    for (k = 0; k < n_mine; k++)
      g_assert_cmpfloat (fabs (mine[k] - ref[k]), <, 1.0e-8 * (1.0 + fabs (ref[k])));
  }

  g_rand_free (rng);
}

static void
test_random_quartic_matches_gsl (void)
{
  GRand *rng = g_rand_new_with_seed (20260706);
  guint trial;

  for (trial = 0; trial < 2000; trial++)
  {
    gdouble a[5];
    gdouble mine[4], ref[4];
    gint n_mine, n_ref, k;
    gint i;

    for (i = 0; i < 5; i++)
      a[i] = g_rand_double_range (rng, -3.0, 3.0);

    if (fabs (a[4]) < 1.0e-2)
      continue;  /* keep it a genuine quartic */

    n_mine = ncm_poly_roots_real_quartic (a, mine);
    n_ref  = _gsl_real_roots (a, 5, ref);

    g_assert_cmpint (n_mine, ==, n_ref);
    _sort (mine, n_mine);

    for (k = 0; k < n_mine; k++)
      g_assert_cmpfloat (fabs (mine[k] - ref[k]), <, 1.0e-7 * (1.0 + fabs (ref[k])));
  }

  g_rand_free (rng);
}

/* Random-trial cross-checks against GSL's general solver now using
 *  random coefficients between -3000 and 3000 */
static void
test_random_cubic_matches_gsl_large (void)
{
  GRand *rng = g_rand_new_with_seed (20260706);
  guint trial;

  for (trial = 0; trial < 2000; trial++)
  {
    gdouble a[4];
    gdouble mine[3], ref[3];
    gint n_mine, n_ref, k;
    gint i;

    for (i = 0; i < 4; i++)
      a[i] = g_rand_double_range (rng, -3000.0, 3000.0);

    if (fabs (a[3]) < 1.0e-2)
      continue;  /* keep it a genuine cubic */

    n_mine = ncm_poly_roots_real_cubic (a, mine);
    n_ref  = _gsl_real_roots (a, 4, ref);

    g_assert_cmpint (n_mine, ==, n_ref);
    _sort (mine, n_mine);

    for (k = 0; k < n_mine; k++)
      g_assert_cmpfloat (fabs (mine[k] - ref[k]), <, 1.0e-8 * (1.0 + fabs (ref[k])));
  }

  g_rand_free (rng);
}

static void
test_random_quartic_matches_gsl_large (void)
{
  GRand *rng = g_rand_new_with_seed (20260706);
  guint trial;

  for (trial = 0; trial < 2000; trial++)
  {
    gdouble a[5];
    gdouble mine[4], ref[4];
    gint n_mine, n_ref, k;
    gint i;

    for (i = 0; i < 5; i++)
      a[i] = g_rand_double_range (rng, -3000.0, 3000.0);

    if (fabs (a[4]) < 1.0e-2)
      continue;  /* keep it a genuine quartic */

    n_mine = ncm_poly_roots_real_quartic (a, mine);
    n_ref  = _gsl_real_roots (a, 5, ref);

    g_assert_cmpint (n_mine, ==, n_ref);
    _sort (mine, n_mine);

    for (k = 0; k < n_mine; k++)
      g_assert_cmpfloat (fabs (mine[k] - ref[k]), <, 1.0e-7 * (1.0 + fabs (ref[k])));
  }

  g_rand_free (rng);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/poly_roots/known_quadratic", &test_known_quadratic);
  g_test_add_func ("/ncm/poly_roots/known_quadratic_no_real_roots", &test_known_quadratic_no_real_roots);
  g_test_add_func ("/ncm/poly_roots/known_cubic", &test_known_cubic);
  g_test_add_func ("/ncm/poly_roots/known_quartic_all_real", &test_known_quartic_all_real);
  g_test_add_func ("/ncm/poly_roots/known_quartic_two_real_two_complex", &test_known_quartic_two_real_two_complex);
  g_test_add_func ("/ncm/poly_roots/quartic_or_lower_degeneracy", &test_quartic_or_lower_degeneracy);
  g_test_add_func ("/ncm/poly_roots/random_cubic_matches_gsl", &test_random_cubic_matches_gsl);
  g_test_add_func ("/ncm/poly_roots/random_quartic_matches_gsl", &test_random_quartic_matches_gsl);
  g_test_add_func ("/ncm/poly_roots/random_cubic_matches_gsl/large", &test_random_cubic_matches_gsl_large);
  g_test_add_func ("/ncm/poly_roots/random_quartic_matches_gsl/large", &test_random_quartic_matches_gsl_large);

  return g_test_run ();
}


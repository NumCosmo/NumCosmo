/***************************************************************************
 *            test_ncm_spline_vec.c
 *
 *  Sat Mar 15 19:53:22 2026
 *  Copyright  2026
 *  Sandro Dias Pinto Vitenti
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_sf.h>

typedef struct _TestNcmSplineVec
{
  NcmSpline *s_base;
  NcmSplineVec *sv;
  NcmVector *xv;
  NcmMatrix *ym;
  guint nknots;
  guint nvec;
  gdouble xi;
  gdouble xf;
} TestNcmSplineVec;

static void
test_ncm_spline_vec_new (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);

  g_assert_true (NCM_IS_SPLINE_VEC (sv));
  g_assert_cmpuint (ncm_spline_vec_get_len (sv), ==, test->nvec);
  g_assert_true (ncm_spline_vec_is_init (sv));

  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_new_gpa (TestNcmSplineVec *test, gconstpointer pdata)
{
  GPtrArray *yv = g_ptr_array_new ();
  guint i;

  for (i = 0; i < test->nvec; i++)
  {
    NcmVector *yv_i = ncm_matrix_get_row (test->ym, i);

    g_ptr_array_add (yv, yv_i);
  }

  NcmSplineVec *sv = ncm_spline_vec_new_gpa (test->s_base, test->xv, yv, TRUE);

  g_assert_true (NCM_IS_SPLINE_VEC (sv));
  g_assert_cmpuint (ncm_spline_vec_get_len (sv), ==, test->nvec);
  g_assert_true (ncm_spline_vec_is_init (sv));

  for (i = 0; i < yv->len; i++)
    ncm_vector_free (g_ptr_array_index (yv, i));

  g_ptr_array_unref (yv);
  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_get_nknots (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  guint i;

  /* Check that get_nknots returns the correct value */
  g_assert_cmpuint (ncm_spline_vec_get_nknots (sv), ==, test->nknots);

  /* Verify it matches the underlying spline length */
  for (i = 0; i < test->nvec; i++)
  {
    NcmSpline *s_i = ncm_spline_vec_peek_spline (sv, i);

    g_assert_cmpuint (ncm_spline_vec_get_nknots (sv), ==, ncm_spline_get_len (s_i));
  }

  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_eval (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  NcmVector *res   = ncm_vector_new (test->nvec);
  const gdouble x  = (test->xi + test->xf) * 0.5;
  guint i;

  ncm_spline_vec_eval (sv, x, res);

  /* Check that each component matches the individual spline evaluation */
  for (i = 0; i < test->nvec; i++)
  {
    NcmSpline *s_i         = ncm_spline_vec_peek_spline (sv, i);
    const gdouble expected = ncm_spline_eval (s_i, x);
    const gdouble computed = ncm_vector_get (res, i);

    g_assert_cmpfloat (fabs (computed - expected), <, 1e-15);
  }

  ncm_vector_free (res);
  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_deriv (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  NcmVector *res   = ncm_vector_new (test->nvec);
  const gdouble x  = (test->xi + test->xf) * 0.5;
  guint i;

  ncm_spline_vec_deriv (sv, x, res);

  /* Check that each component matches the individual spline derivative */
  for (i = 0; i < test->nvec; i++)
  {
    NcmSpline *s_i         = ncm_spline_vec_peek_spline (sv, i);
    const gdouble expected = ncm_spline_eval_deriv (s_i, x);
    const gdouble computed = ncm_vector_get (res, i);

    g_assert_cmpfloat (fabs (computed - expected), <, 1e-15);
  }

  ncm_vector_free (res);
  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_integ (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  NcmVector *res   = ncm_vector_new (test->nvec);
  const gdouble xi = test->xi + (test->xf - test->xi) * 0.25;
  const gdouble xf = test->xi + (test->xf - test->xi) * 0.75;
  guint i;

  ncm_spline_vec_integ (sv, xi, xf, res);

  /* Check that each component matches the individual spline integration */
  for (i = 0; i < test->nvec; i++)
  {
    NcmSpline *s_i         = ncm_spline_vec_peek_spline (sv, i);
    const gdouble expected = ncm_spline_eval_integ (s_i, xi, xf);
    const gdouble computed = ncm_vector_get (res, i);

    g_assert_cmpfloat (fabs (computed - expected), <, 1e-15);
  }

  ncm_vector_free (res);
  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_set (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, FALSE);

  g_assert_false (ncm_spline_vec_is_init (sv));

  ncm_spline_vec_prepare (sv);

  g_assert_true (ncm_spline_vec_is_init (sv));

  /* Reset with new data */
  ncm_spline_vec_set (sv, test->xv, test->ym, TRUE);

  g_assert_true (ncm_spline_vec_is_init (sv));
  g_assert_cmpuint (ncm_spline_vec_get_len (sv), ==, test->nvec);

  ncm_spline_vec_free (sv);
}

static void
test_ncm_spline_vec_free (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv     = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  NcmSplineVec *sv_ref = ncm_spline_vec_ref (sv);

  ncm_spline_vec_clear (&sv);
  g_assert_null (sv);
  g_assert_nonnull (sv_ref);

  ncm_spline_vec_free (sv_ref);
}

/* Against the closed forms, through both the NcmVector and the GArray variants. x^2 and
 * x^3 are reproduced by the not-a-knot cubic to rounding; the cos(x) bounds are 10 times
 * the measured errors on these 100 knots. */
static void
test_ncm_spline_vec_values (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv   = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  NcmVector *res     = ncm_vector_new (test->nvec);
  GArray *res_a      = NULL;
  const gdouble xs[] = {0.37, 5.0, 9.93};
  const gdouble ab[] = {2.0, 7.0, 0.37, 9.93};
  guint j, k;

  for (j = 0; j < G_N_ELEMENTS (xs); j++)
  {
    const gdouble x = xs[j];

    ncm_spline_vec_eval (sv, x, res);
    ncm_spline_vec_eval_array (sv, x, &res_a);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 0), ==, x * x, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 1), ==, x * x * x, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 2), ==, cos (x), 0.0, 1.0e-5);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (g_array_index (res_a, gdouble, k), ==, ncm_vector_get (res, k));

    ncm_spline_vec_deriv (sv, x, res);
    ncm_spline_vec_deriv_array (sv, x, &res_a);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 0), ==, 2.0 * x, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 1), ==, 3.0 * x * x, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 2), ==, -sin (x), 0.0, 1.0e-3);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (g_array_index (res_a, gdouble, k), ==, ncm_vector_get (res, k));
  }

  for (j = 0; j < G_N_ELEMENTS (ab); j += 2)
  {
    const gdouble a = ab[j];
    const gdouble b = ab[j + 1];

    ncm_spline_vec_integ (sv, a, b, res);
    ncm_spline_vec_integ_array (sv, a, b, &res_a);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 0), ==, (b * b * b - a * a * a) / 3.0, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 1), ==, (gsl_pow_4 (b) - gsl_pow_4 (a)) / 4.0, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (res, 2), ==, sin (b) - sin (a), 0.0, 1.0e-6);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (g_array_index (res_a, gdouble, k), ==, ncm_vector_get (res, k));

    /* Reversed limits change the sign */
    ncm_spline_vec_integ_array (sv, b, a, &res_a);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (g_array_index (res_a, gdouble, k), ==, -ncm_vector_get (res, k));
  }

  g_array_unref (res_a);
  ncm_vector_free (res);
  ncm_spline_vec_free (sv);
}

/* The matrix and the GPtrArray constructors build the same components. */
static void
test_ncm_spline_vec_gpa_equiv (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSplineVec *sv = ncm_spline_vec_new (test->s_base, test->xv, test->ym, TRUE);
  GPtrArray *yv    = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  NcmSplineVec *sv_gpa;
  guint i, k;

  for (k = 0; k < test->nvec; k++)
    g_ptr_array_add (yv, ncm_matrix_get_row (test->ym, k));

  sv_gpa = ncm_spline_vec_new_gpa (test->s_base, test->xv, yv, TRUE);

  for (i = 0; i < 20; i++)
  {
    const gdouble x = test->xi + (test->xf - test->xi) * i / 19.0;

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (ncm_spline_eval (ncm_spline_vec_peek_spline (sv_gpa, k), x), ==,
                         ncm_spline_eval (ncm_spline_vec_peek_spline (sv, k), x));
  }

  g_ptr_array_unref (yv);
  ncm_spline_vec_free (sv_gpa);
  ncm_spline_vec_free (sv);
}

/* A component type that does not use the interval index evaluates, differentiates and
 * integrates without it. */
static void
test_ncm_spline_vec_bspline (TestNcmSplineVec *test, gconstpointer pdata)
{
  NcmSpline *s_bs    = NCM_SPLINE (ncm_spline_bspline_new (8));
  NcmSplineVec *sv   = ncm_spline_vec_new (s_bs, test->xv, test->ym, TRUE);
  NcmVector *res     = ncm_vector_new (test->nvec);
  const gdouble xs[] = {0.37, 5.0, 9.93};
  guint j, k;

  for (j = 0; j < G_N_ELEMENTS (xs); j++)
  {
    ncm_spline_vec_eval (sv, xs[j], res);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (ncm_vector_get (res, k), ==, ncm_spline_eval (ncm_spline_vec_peek_spline (sv, k), xs[j]));

    ncm_spline_vec_deriv (sv, xs[j], res);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (ncm_vector_get (res, k), ==, ncm_spline_eval_deriv (ncm_spline_vec_peek_spline (sv, k), xs[j]));

    ncm_spline_vec_integ (sv, 1.0, xs[j], res);

    for (k = 0; k < test->nvec; k++)
      g_assert_cmpfloat (ncm_vector_get (res, k), ==, ncm_spline_eval_integ (ncm_spline_vec_peek_spline (sv, k), 1.0, xs[j]));
  }

  ncm_vector_free (res);
  ncm_spline_vec_free (sv);
  ncm_spline_free (s_bs);
}

static void
test_ncm_spline_vec_empty (void)
{
  g_test_trap_subprocess ("/ncm/spline_vec/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*no components*");
}

static void
test_ncm_spline_vec_empty_subprocess (void)
{
  NcmSpline *s  = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmVector *xv = ncm_vector_new (10);
  GPtrArray *yv = g_ptr_array_new ();

  ncm_vector_set_all (xv, 0.0);
  ncm_spline_vec_new_gpa (s, xv, yv, TRUE);
}

static void
test_ncm_spline_vec_setup (TestNcmSplineVec *test, gconstpointer pdata)
{
  guint i;

  test->nknots = 100;
  test->nvec   = 3;
  test->xi     = 0.0;
  test->xf     = 10.0;

  test->s_base = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  test->xv     = ncm_vector_new (test->nknots);
  test->ym     = ncm_matrix_new (test->nvec, test->nknots);

  for (i = 0; i < test->nknots; i++)
  {
    const gdouble x = test->xi + (test->xf - test->xi) * i / (test->nknots - 1.0);

    ncm_vector_set (test->xv, i, x);

    /* Initialize y values for each component with different functions */
    ncm_matrix_set (test->ym, 0, i, x * x);     /* y0 = x^2 */
    ncm_matrix_set (test->ym, 1, i, x * x * x); /* y1 = x^3 */
    ncm_matrix_set (test->ym, 2, i, cos (x));   /* y2 = cos(x) */
  }
}

static void
test_ncm_spline_vec_teardown (TestNcmSplineVec *test, gconstpointer pdata)
{
  ncm_spline_free (test->s_base);
  ncm_vector_free (test->xv);
  ncm_matrix_free (test->ym);
}

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/spline_vec/new", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_new,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/new_gpa", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_new_gpa,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/get_nknots", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_get_nknots,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/eval", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_eval,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/deriv", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_deriv,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/integ", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_integ,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/set", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_set,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/free", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_free,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/values", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_values,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/gpa_equiv", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_gpa_equiv,
              &test_ncm_spline_vec_teardown);

  g_test_add ("/ncm/spline_vec/bspline", TestNcmSplineVec, NULL,
              &test_ncm_spline_vec_setup,
              &test_ncm_spline_vec_bspline,
              &test_ncm_spline_vec_teardown);

  g_test_add_func ("/ncm/spline_vec/empty", &test_ncm_spline_vec_empty);
  g_test_add_func ("/ncm/spline_vec/empty/subprocess", &test_ncm_spline_vec_empty_subprocess);

  g_test_run ();

  return 0;
}


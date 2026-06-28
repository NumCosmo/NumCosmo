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

  g_test_run ();

  return 0;
}


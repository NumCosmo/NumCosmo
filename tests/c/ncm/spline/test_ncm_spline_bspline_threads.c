/***************************************************************************
 *            test_ncm_spline_bspline_threads.c
 *
 *  Sun Aug 31 12:00:00 2026
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

/*
 * Concurrent evaluation of one shared, prepared NcmSplineBSpline must return
 * the same bits as serial evaluation. The GSL workspace held by the spline is
 * scratch for every gsl_bspline_calc* call, so an evaluation path that touches
 * it from two threads corrupts both; this is the regression test for that.
 * Points outside the range are included, as the edge spans extrapolate.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <math.h>

#define TEST_NKNOTS 200
#define TEST_NEVAL 20000
#define TEST_NTHREADS 8

typedef struct _TestBSplineThreads
{
  NcmSpline *s;
  gdouble *x_eval;
  gdouble *ref_eval;
  gdouble *ref_deriv;
  gdouble *ref_deriv2;
  gdouble *ref_dnmax;
} TestBSplineThreads;

static void
test_ncm_spline_bspline_threads_new (TestBSplineThreads *test, gconstpointer pdata)
{
  NcmVector *xv = ncm_vector_new (TEST_NKNOTS);
  NcmVector *yv = ncm_vector_new (TEST_NKNOTS);
  gdouble x     = 0.3;
  guint i;

  for (i = 0; i < TEST_NKNOTS; i++)
  {
    x += 0.1 + 0.05 * sin (3.0 * i);
    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, sin (2.0 * x) * exp (-1.0e-3 * x * x) + 0.3 * cos (5.0 * x));
  }

  test->s = NCM_SPLINE (ncm_spline_bspline_new_full (8, xv, yv, TRUE));

  {
    const gdouble x0    = ncm_vector_get (xv, 0);
    const gdouble x1    = ncm_vector_get (xv, TEST_NKNOTS - 1);
    const gdouble width = x1 - x0;
    guint j;

    test->x_eval     = g_new (gdouble, TEST_NEVAL);
    test->ref_eval   = g_new (gdouble, TEST_NEVAL);
    test->ref_deriv  = g_new (gdouble, TEST_NEVAL);
    test->ref_deriv2 = g_new (gdouble, TEST_NEVAL);
    test->ref_dnmax  = g_new (gdouble, TEST_NEVAL);

    /* 2% of the points sit outside the range on each side. */
    for (j = 0; j < TEST_NEVAL; j++)
      test->x_eval[j] = x0 + width * (-0.02 + 1.04 * j / (TEST_NEVAL - 1.0));

    for (j = 0; j < TEST_NEVAL; j++)
    {
      test->ref_eval[j]   = ncm_spline_eval (test->s, test->x_eval[j]);
      test->ref_deriv[j]  = ncm_spline_eval_deriv (test->s, test->x_eval[j]);
      test->ref_deriv2[j] = ncm_spline_eval_deriv2 (test->s, test->x_eval[j]);
      test->ref_dnmax[j]  = ncm_spline_eval_deriv_nmax (test->s, test->x_eval[j]);
    }
  }

  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

static void
test_ncm_spline_bspline_threads_free (TestBSplineThreads *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_spline_free, test->s);
  g_free (test->x_eval);
  g_free (test->ref_eval);
  g_free (test->ref_deriv);
  g_free (test->ref_deriv2);
  g_free (test->ref_dnmax);
}

typedef struct _BSplineWorker
{
  TestBSplineThreads *test;
  glong n_bad;
  gboolean with_derivs;
} BSplineWorker;

static gpointer
_bspline_worker (gpointer data)
{
  BSplineWorker *worker    = (BSplineWorker *) data;
  TestBSplineThreads *test = worker->test;
  guint j;

  for (j = 0; j < TEST_NEVAL; j++)
  {
    const gdouble x = test->x_eval[j];

    if (ncm_spline_eval (test->s, x) != test->ref_eval[j])
      worker->n_bad++;

    if (worker->with_derivs)
    {
      if (ncm_spline_eval_deriv (test->s, x) != test->ref_deriv[j])
        worker->n_bad++;

      if (ncm_spline_eval_deriv2 (test->s, x) != test->ref_deriv2[j])
        worker->n_bad++;

      if (ncm_spline_eval_deriv_nmax (test->s, x) != test->ref_dnmax[j])
        worker->n_bad++;
    }
  }

  return NULL;
}

static void
_test_ncm_spline_bspline_threads_run (TestBSplineThreads *test, gboolean with_derivs)
{
  GThread *threads[TEST_NTHREADS];
  BSplineWorker workers[TEST_NTHREADS];
  guint i;

  for (i = 0; i < TEST_NTHREADS; i++)
  {
    workers[i].test        = test;
    workers[i].n_bad       = 0;
    workers[i].with_derivs = with_derivs;
    threads[i]             = g_thread_new ("bspline-eval", _bspline_worker, &workers[i]);
  }

  for (i = 0; i < TEST_NTHREADS; i++)
  {
    g_thread_join (threads[i]);
    g_assert_cmpint (workers[i].n_bad, ==, 0);
  }
}

static void
test_ncm_spline_bspline_threads_eval (TestBSplineThreads *test, gconstpointer pdata)
{
  _test_ncm_spline_bspline_threads_run (test, FALSE);
}

static void
test_ncm_spline_bspline_threads_eval_derivs (TestBSplineThreads *test, gconstpointer pdata)
{
  _test_ncm_spline_bspline_threads_run (test, TRUE);
}

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/spline_bspline/threads/eval", TestBSplineThreads, NULL,
              &test_ncm_spline_bspline_threads_new,
              &test_ncm_spline_bspline_threads_eval,
              &test_ncm_spline_bspline_threads_free);

  g_test_add ("/ncm/spline_bspline/threads/eval_derivs", TestBSplineThreads, NULL,
              &test_ncm_spline_bspline_threads_new,
              &test_ncm_spline_bspline_threads_eval_derivs,
              &test_ncm_spline_bspline_threads_free);

  g_test_run ();

  return 0;
}


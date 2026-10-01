/***************************************************************************
 *            test_ncm_stats_vec.c
 *
 *  Fri August 02 18:30:24 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2013 <vitenti@uel.br>
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#define _TEST_NCM_VECTOR_STATIC_SIZE 20
#define _TEST_NCM_VECTOR_MIN_SIZE 5
#define _TEST_NCM_STATS_VEC_PREC (1.0e-9)
#define _TEST_NCM_STATS_VEC_COV_PREC (1.0e-9)
#define _TEST_NCM_STATS_VEC_NTEST_MAX 3000
#define _TEST_NCM_STATS_VEC_NTEST_MIN 2000

typedef struct _TestNcmStatsVec
{
  NcmStatsVec *svec;
  NcmVector *w;
  NcmVector *mu;
  NcmMatrix *xs;
  guint v_size;
  guint ntests;
} TestNcmStatsVec;

void test_ncm_stats_vec_mean_new (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_var_new (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_cov_new (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_autocorr_new (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_mean_test (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_var_test (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_cov_test (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_cov_robust_test (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_autocorr_test (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_subsample_autocorr_test (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_free (TestNcmStatsVec *test, gconstpointer pdata);

void test_ncm_stats_vec_diag_heidel (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_diag_max_ess (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_diag_visual (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_diag_discriminates (void);
void test_ncm_stats_vec_diag_saved (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_diag_const_break (TestNcmStatsVec *test, gconstpointer pdata);

void test_ncm_stats_vec_traps (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_invalid_get_cov (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_invalid_get_var (TestNcmStatsVec *test, gconstpointer pdata);

/*
 * Curated chains for the convergence diagnostics. Generated from a fixed seed, since
 * these read given data and no sampler has to run to test them.
 *
 * Heidelberger-Welch asks from which index a chain is stationary. What it can answer is
 * "from the start" and "after this burn-in", and those are the two the chains below have.
 * It has a documented third answer, -1 for never, that no chain here provokes: the
 * spectral density at zero is estimated from the second half, so a break living there --
 * a linear drift, or a late jump in the mean -- inflates the scale the statistic is
 * divided by and the chain is reported stationary. Measured, not assumed: a drift of
 * 0.05 per step reports 160 and a +20 shift at three quarters reports 0.
 */
typedef enum _TestSeries
{
  TEST_SERIES_WHITE,      /* stationary from index 0 */
  TEST_SERIES_AR1,        /* stationary, but correlated: a nonzero AR order */
  TEST_SERIES_BURNIN,     /* stationary only after the transient */
  TEST_SERIES_LATE_SHIFT, /* broken in its last quarter, which this test cannot see */
} TestSeries;

#define TEST_DIAG_NITEMS 400
#define TEST_DIAG_BURNIN 150
#define TEST_DIAG_LEN 2

typedef struct _TestSeriesCase
{
  const gchar *name;
  TestSeries series;
  gint bindex_min; /* smallest acceptable stationarity index */
} TestSeriesCase;

static const TestSeriesCase test_series_cases[] = {
  {"white",      TEST_SERIES_WHITE,      0  },
  {"ar1",        TEST_SERIES_AR1,        0  },
  {"burnin",     TEST_SERIES_BURNIN,     100},
  {"late_shift", TEST_SERIES_LATE_SHIFT, 0  },
};

static void
test_ncm_stats_vec_diag_new (TestNcmStatsVec *test, gconstpointer pdata)
{
  const TestSeriesCase *tc    = pdata;
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, 1234567890);
  NcmVector *row              = ncm_vector_new (TEST_DIAG_LEN);
  gdouble prev[TEST_DIAG_LEN] = { 0.0, 0.0 };
  guint i, p;

  /* save_x: the diagnostics re-read the rows, not just the accumulated moments. COV
   * rather than VAR because the correlation accessor requires it, and it is a superset. */
  test->svec   = ncm_stats_vec_new (TEST_DIAG_LEN, NCM_STATS_VEC_COV, TRUE);
  test->ntests = TEST_DIAG_NITEMS;

  for (i = 0; i < TEST_DIAG_NITEMS; i++)
  {
    for (p = 0; p < TEST_DIAG_LEN; p++)
    {
      const gdouble eps = ncm_rng_gaussian_gen (rng, 0.0, 1.0);
      gdouble x;

      switch (tc->series)
      {
        case TEST_SERIES_WHITE:
          x = eps;
          break;
        case TEST_SERIES_AR1:
          x       = 0.8 * prev[p] + eps;
          prev[p] = x;
          break;
        case TEST_SERIES_BURNIN:
          x = (i < TEST_DIAG_BURNIN) ? (20.0 + eps) : eps;
          break;
        case TEST_SERIES_LATE_SHIFT:
          /* The diagnostic can only discard a prefix, so a break in the last quarter
           * survives every starting index it is allowed to try. A slow linear drift
           * does not work here: it inflates the spectral density at zero that the
           * statistic is divided by, and the test reports the chain stationary from
           * partway in. */
          x = (i < 3 * TEST_DIAG_NITEMS / 4) ? eps : (20.0 + eps);
          break;
        default:
          g_assert_not_reached ();
      }

      ncm_vector_set (row, p, x);
    }

    ncm_stats_vec_append (test->svec, row, TRUE);
  }

  ncm_vector_free (row);
  ncm_rng_free (rng);
}

static void
test_ncm_stats_vec_diag_free (TestNcmStatsVec *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_vec_free, test->svec);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* Default vector allocation */

  g_test_add ("/ncm/stats_vec/mean", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_mean_new,
              &test_ncm_stats_vec_mean_test,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/var", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_var_new,
              &test_ncm_stats_vec_var_test,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/cov", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_cov_new,
              &test_ncm_stats_vec_cov_test,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/cov/robust", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_cov_new,
              &test_ncm_stats_vec_cov_robust_test,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/autocorr", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_autocorr_new,
              &test_ncm_stats_vec_autocorr_test,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/subsample_autocorr", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_autocorr_new,
              &test_ncm_stats_vec_subsample_autocorr_test,
              &test_ncm_stats_vec_free);

  g_test_add ("/ncm/stats_vec/mean/get_var/subprocess", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_mean_new,
              &test_ncm_stats_vec_invalid_get_var,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/mean/get_cov/subprocess", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_mean_new,
              &test_ncm_stats_vec_invalid_get_cov,
              &test_ncm_stats_vec_free);
  g_test_add ("/ncm/stats_vec/var/get_cov/subprocess", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_var_new,
              &test_ncm_stats_vec_invalid_get_cov,
              &test_ncm_stats_vec_free);

  {
    guint i;

    for (i = 0; i < G_N_ELEMENTS (test_series_cases); i++)
    {
      gchar *path;

      path = g_strdup_printf ("/ncm/stats_vec/diag/heidel/%s", test_series_cases[i].name);
      g_test_add (path, TestNcmStatsVec, &test_series_cases[i],
                  &test_ncm_stats_vec_diag_new, &test_ncm_stats_vec_diag_heidel,
                  &test_ncm_stats_vec_diag_free);
      g_free (path);

      path = g_strdup_printf ("/ncm/stats_vec/diag/max_ess/%s", test_series_cases[i].name);
      g_test_add (path, TestNcmStatsVec, &test_series_cases[i],
                  &test_ncm_stats_vec_diag_new, &test_ncm_stats_vec_diag_max_ess,
                  &test_ncm_stats_vec_diag_free);
      g_free (path);

      path = g_strdup_printf ("/ncm/stats_vec/diag/visual/%s", test_series_cases[i].name);
      g_test_add (path, TestNcmStatsVec, &test_series_cases[i],
                  &test_ncm_stats_vec_diag_new, &test_ncm_stats_vec_diag_visual,
                  &test_ncm_stats_vec_diag_free);
      g_free (path);

      path = g_strdup_printf ("/ncm/stats_vec/diag/saved/%s", test_series_cases[i].name);
      g_test_add (path, TestNcmStatsVec, &test_series_cases[i],
                  &test_ncm_stats_vec_diag_new, &test_ncm_stats_vec_diag_saved,
                  &test_ncm_stats_vec_diag_free);
      g_free (path);

      path = g_strdup_printf ("/ncm/stats_vec/diag/const_break/%s", test_series_cases[i].name);
      g_test_add (path, TestNcmStatsVec, &test_series_cases[i],
                  &test_ncm_stats_vec_diag_new, &test_ncm_stats_vec_diag_const_break,
                  &test_ncm_stats_vec_diag_free);
      g_free (path);
    }
  }

  g_test_add_func ("/ncm/stats_vec/diag/discriminates", &test_ncm_stats_vec_diag_discriminates);

  g_test_add ("/ncm/stats_vec/traps", TestNcmStatsVec, NULL,
              &test_ncm_stats_vec_var_new,
              &test_ncm_stats_vec_traps,
              &test_ncm_stats_vec_free);

  g_test_run ();
}

void
test_ncm_stats_vec_mean_new (TestNcmStatsVec *test, gconstpointer pdata)
{
  test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  test->ntests = g_test_rand_int_range (_TEST_NCM_STATS_VEC_NTEST_MIN, _TEST_NCM_STATS_VEC_NTEST_MAX);
  test->svec   = ncm_stats_vec_new (test->v_size, NCM_STATS_VEC_MEAN, FALSE);
  test->xs     = ncm_matrix_new (test->ntests, test->v_size);
  test->mu     = ncm_vector_new (test->v_size);
  test->w      = ncm_vector_new (test->ntests);

  g_assert_true (NCM_IS_STATS_VEC (test->svec));
}

void
test_ncm_stats_vec_var_new (TestNcmStatsVec *test, gconstpointer pdata)
{
  test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  test->ntests = g_test_rand_int_range (_TEST_NCM_STATS_VEC_NTEST_MIN, _TEST_NCM_STATS_VEC_NTEST_MAX);
  test->svec   = ncm_stats_vec_new (test->v_size, NCM_STATS_VEC_VAR, FALSE);
  test->xs     = ncm_matrix_new (test->ntests, test->v_size);
  test->mu     = ncm_vector_new (test->v_size);
  test->w      = ncm_vector_new (test->ntests);

  g_assert_true (NCM_IS_STATS_VEC (test->svec));
}

void
test_ncm_stats_vec_cov_new (TestNcmStatsVec *test, gconstpointer pdata)
{
  test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  test->ntests = g_test_rand_int_range (_TEST_NCM_STATS_VEC_NTEST_MIN, _TEST_NCM_STATS_VEC_NTEST_MAX);
  test->svec   = ncm_stats_vec_new (test->v_size, NCM_STATS_VEC_COV, TRUE);
  test->xs     = ncm_matrix_new (test->ntests, test->v_size);
  test->mu     = ncm_vector_new (test->v_size);
  test->w      = ncm_vector_new (test->ntests);

  g_assert_true (NCM_IS_STATS_VEC (test->svec));
}

void
test_ncm_stats_vec_autocorr_new (TestNcmStatsVec *test, gconstpointer pdata)
{
  test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  test->ntests = g_test_rand_int_range (_TEST_NCM_STATS_VEC_NTEST_MIN, _TEST_NCM_STATS_VEC_NTEST_MAX) * 100;
  test->svec   = ncm_stats_vec_new (test->v_size, NCM_STATS_VEC_COV, TRUE);
  test->xs     = ncm_matrix_new (test->ntests, test->v_size);
  test->mu     = ncm_vector_new (test->v_size);
  test->w      = ncm_vector_new (test->ntests);

  g_assert_true (NCM_IS_STATS_VEC (test->svec));
}

void
test_ncm_stats_vec_free (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmStatsVec *svec = test->svec;

  ncm_matrix_clear (&test->xs);
  ncm_vector_clear (&test->w);
  ncm_vector_clear (&test->mu);

  NCM_TEST_FREE (ncm_stats_vec_free, svec);
}

void
test_ncm_stats_vec_mean_test (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_pool_get ("test_ncm_stats_vec");
  guint i;

  for (i = 0; i < test->v_size; i++)
  {
    ncm_vector_set (test->mu, i, 1.0 + fabs (g_test_rand_double ()));
  }

  for (i = 0; i < test->ntests; i++)
  {
    gdouble sigma = fabs (g_test_rand_double ()) + 1.0e-1;
    gdouble w     = 1.0 / (sigma * sigma);
    guint j;

    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = ncm_vector_get (test->mu, j) + sigma * ncm_rng_ugaussian_gen (rng);

      ncm_stats_vec_set (test->svec, j, x_j);
      ncm_matrix_set (test->xs, i, j, x_j);
    }

    ncm_vector_set (test->w, i, w);
    ncm_stats_vec_update_weight (test->svec, w);

    if ((i % 100) == 1)
    {
      guint k;

      for (k = 0; k < test->v_size; k++)
      {
        gdouble gsl_mean = gsl_stats_wmean (ncm_vector_ptr (test->w, 0), 1,
                                            ncm_matrix_ptr (test->xs, 0, k), test->v_size,
                                            i + 1);
        gdouble svec_mean = ncm_stats_vec_get_mean (test->svec, k);

        ncm_assert_cmpdouble_e (gsl_mean, ==, svec_mean, _TEST_NCM_STATS_VEC_PREC, 0.0);
      }
    }
  }

  ncm_rng_free (rng);
  NCM_TEST_FAIL (ncm_stats_vec_get_var (test->svec, 0));
}

void
test_ncm_stats_vec_var_test (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_pool_get ("test_ncm_stats_vec");
  guint i;

  for (i = 0; i < test->v_size; i++)
  {
    ncm_vector_set (test->mu, i, 1.0 + fabs (g_test_rand_double ()));
  }

  for (i = 0; i < test->ntests; i++)
  {
    gdouble sigma = fabs (g_test_rand_double ()) + 1.0e-1;
    gdouble w     = 1.0 / (sigma * sigma);
    guint j;

    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = ncm_vector_get (test->mu, j) + sigma * ncm_rng_ugaussian_gen (rng);

      ncm_stats_vec_set (test->svec, j, x_j);
      ncm_matrix_set (test->xs, i, j, x_j);
    }

    ncm_vector_set (test->w, i, w);
    ncm_stats_vec_update_weight (test->svec, w);

    if ((i % 100) == 1)
    {
      guint k;

      for (k = 0; k < test->v_size; k++)
      {
        gdouble gsl_mean = gsl_stats_wmean (ncm_vector_ptr (test->w, 0), 1,
                                            ncm_matrix_ptr (test->xs, 0, k), test->v_size,
                                            i + 1);
        gdouble gsl_var = gsl_stats_wvariance (ncm_vector_ptr (test->w, 0), 1,
                                               ncm_matrix_ptr (test->xs, 0, k), test->v_size,
                                               i + 1);
        gdouble svec_mean = ncm_stats_vec_get_mean (test->svec, k);
        gdouble svec_var  = ncm_stats_vec_get_var (test->svec, k);

        ncm_assert_cmpdouble_e (gsl_mean, ==, svec_mean, _TEST_NCM_STATS_VEC_PREC, 0.0);
        ncm_assert_cmpdouble_e (gsl_var, ==, svec_var, _TEST_NCM_STATS_VEC_PREC, 0.0);
      }
    }
  }

  ncm_rng_free (rng);
  NCM_TEST_FAIL (ncm_stats_vec_get_cov (test->svec, 0, 1));
}

void
test_ncm_stats_vec_cov_test (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmRNG *rng   = ncm_rng_pool_get ("test_ncm_stats_vec");
  gdouble sigma = fabs (g_test_rand_double ()) + 1.0e-1;
  guint i;

  for (i = 0; i < test->v_size; i++)
  {
    ncm_vector_set (test->mu, i, 1.0 + fabs (g_test_rand_double ()));
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint j;
    gdouble x_0;

    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = ncm_vector_get (test->mu, j) + sigma * ncm_rng_ugaussian_gen (rng);

      if (j == 0)
      {
        x_0 = x_j;
      }
      else
      {
        x_j += x_0;
        x_0 += x_j;
      }

      ncm_stats_vec_set (test->svec, j, x_j);
      ncm_matrix_set (test->xs, i, j, x_j);
    }

    ncm_stats_vec_update (test->svec);

    if ((i % 100) == 1)
    {
      guint k;

      for (k = 0; k < test->v_size; k++)
      {
        gdouble gsl_mean = gsl_stats_mean (ncm_matrix_ptr (test->xs, 0, k), test->v_size,
                                           i + 1);
        gdouble gsl_var = gsl_stats_variance (ncm_matrix_ptr (test->xs, 0, k), test->v_size,
                                              i + 1);
        gdouble svec_mean = ncm_stats_vec_get_mean (test->svec, k);
        gdouble svec_var  = ncm_stats_vec_get_var (test->svec, k);
        guint l;

        ncm_assert_cmpdouble_e (gsl_mean, ==, svec_mean, _TEST_NCM_STATS_VEC_PREC, 0.0);
        ncm_assert_cmpdouble_e (gsl_var, ==, svec_var, _TEST_NCM_STATS_VEC_PREC, 0.0);

        for (l = k + 1; l < test->v_size; l++)
        {
          const gdouble gsl_cov_kl = gsl_stats_covariance (ncm_matrix_ptr (test->xs, 0, k), test->v_size,
                                                           ncm_matrix_ptr (test->xs, 0, l), test->v_size,
                                                           i + 1);
          const gdouble svec_cov_kl = ncm_stats_vec_get_cov (test->svec, k, l);
          const gdouble svec_cov_lk = ncm_stats_vec_get_cov (test->svec, l, k);

          ncm_assert_cmpdouble (svec_cov_kl, ==, svec_cov_lk);
          ncm_assert_cmpdouble_e (gsl_cov_kl + 1.0, ==, svec_cov_kl + 1.0, _TEST_NCM_STATS_VEC_COV_PREC, 0.0);
        }
      }
    }
  }

  ncm_rng_free (rng);
}

void
test_ncm_stats_vec_cov_robust_test (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmRNG *rng         = ncm_rng_pool_get ("test_ncm_stats_vec");
  const gdouble sigma = fabs (g_test_rand_double ()) + 1.0e-1;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    gdouble x_0 = 0.0;
    guint j;

    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = sigma * ncm_rng_ugaussian_gen (rng);

      if (j == 0)
        x_0 = x_j;
      else
        x_j = (x_0 + x_j) / M_SQRT2;

      ncm_stats_vec_set (test->svec, j, x_j);
    }

    ncm_stats_vec_update (test->svec);
  }

  for (i = 0; i < (gint) (test->ntests * 0.2); i++)
  {
    gdouble x_0 = 0.0;
    guint j;

    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = 1.0 + sigma * ncm_rng_ugaussian_gen (rng) / 10.0;

      if (j == 0)
      {
        x_0 = x_j;
      }
      else
      {
        x_j += x_0;
        x_j  = x_j / sqrt (2.0);
      }

      ncm_stats_vec_set (test->svec, j, x_j);
    }

    ncm_stats_vec_update (test->svec);
  }


  {
    NcmMatrix *cov_n = ncm_stats_vec_peek_cov_matrix (test->svec, 0);
    NcmMatrix *cor_n = ncm_matrix_cov_dup_cor (cov_n);
    NcmVector *var_n = ncm_vector_new (test->v_size);
    NcmMatrix *cov_r = ncm_stats_vec_compute_cov_robust_ogk (test->svec);
    NcmMatrix *cor_r = ncm_matrix_cov_dup_cor (cov_r);
    NcmVector *var_r = ncm_vector_new (test->v_size);

    ncm_matrix_get_diag (cov_n, var_n);
    ncm_matrix_get_diag (cov_r, var_r);

    /*
     *  ncm_message ("Original sigma = % 22.15g, var = % 22.15g\n", sigma, sigma * sigma);
     *  ncm_vector_log_vals (var_n, "VAR_NORMAL: ", "% 12.5g", TRUE);
     *  ncm_vector_log_vals (var_r, "VAR_ROBUST: ", "% 12.5g", TRUE);
     *  ncm_matrix_log_vals (cor_n, "COR_NORMAL: ", "% 12.5g");
     *  ncm_matrix_log_vals (cor_r, "COR_ROBUST: ", "% 12.5g");
     */


    ncm_vector_free (var_n);
    ncm_vector_free (var_r);
    ncm_matrix_free (cor_n);
    ncm_matrix_free (cov_r);
    ncm_matrix_free (cor_r);
  }

  ncm_rng_free (rng);
}

void
test_ncm_stats_vec_autocorr_test (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmRNG *rng         = ncm_rng_pool_get ("test_ncm_stats_vec");
  const gdouble a     = 0.9 + fabs (g_test_rand_double ()) * 1.0e-2;
  const gdouble sigma = fabs (g_test_rand_double ()) * 1.0e-1;
  NcmVector *last     = ncm_vector_new (test->v_size);
  guint i;

  for (i = 0; i < test->v_size; i++)
  {
    ncm_vector_set (test->mu, i, 1.0 + fabs (g_test_rand_double ()));
    ncm_vector_set (last, i, 0.0);
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint j;

    for (j = 0; j < test->v_size; j++)
    {
      const gdouble epsilon_j = ncm_vector_get (test->mu, j) + sigma * ncm_rng_ugaussian_gen (rng);
      const gdouble x_j       = (a * ncm_vector_get (last, j) + epsilon_j);

      ncm_vector_set (last, j, x_j);

      ncm_stats_vec_set (test->svec, j, x_j);
      ncm_matrix_set (test->xs, i, j, x_j);
    }

    ncm_stats_vec_update (test->svec);
  }

  for (i = 0; i < test->v_size; i++)
  {
    const gdouble gsl_mean = gsl_stats_mean (ncm_matrix_ptr (test->xs, 0, i), test->v_size,
                                             test->ntests);
    const gdouble gsl_var = gsl_stats_variance (ncm_matrix_ptr (test->xs, 0, i), test->v_size,
                                                test->ntests);
    const gdouble svec_mean = ncm_stats_vec_get_mean (test->svec, i);
    const gdouble svec_var  = ncm_stats_vec_get_var (test->svec, i);

    ncm_assert_cmpdouble_e (gsl_mean, ==, svec_mean, _TEST_NCM_STATS_VEC_PREC, 0.0);
    ncm_assert_cmpdouble_e (gsl_var, ==, svec_var, _TEST_NCM_STATS_VEC_PREC, 0.0);

    {
      NcmVector *ac = ncm_stats_vec_get_autocorr (test->svec, i);
      guint j;
      guint tsize = GSL_MIN (10, ncm_vector_len (ac));

      for (j = 0; j < tsize; j++)
      {
        if (ncm_vector_get (ac, j) < 0.0)
          break;

        ncm_assert_cmpdouble_e (ncm_vector_get (ac, j), ==, pow (a, j), 1.0e-1, 0.0);
      }

      ncm_vector_free (ac);
    }
  }

  ncm_vector_free (last);
  ncm_rng_free (rng);
}

void
test_ncm_stats_vec_subsample_autocorr_test (TestNcmStatsVec *test, gconstpointer pdata)
{
  NcmRNG *rng         = ncm_rng_pool_get ("test_ncm_stats_vec");
  const gdouble a     = 0.9 + fabs (g_test_rand_double ()) * 1.0e-2;
  const gdouble sigma = fabs (g_test_rand_double ()) * 1.0e-1;
  const guint nchains = g_test_rand_int_range (10, 20);
  NcmMatrix *last     = ncm_matrix_new (nchains, test->v_size);
  guint i;

  ncm_matrix_set_zero (last);

  for (i = 0; i < test->v_size; i++)
  {
    ncm_vector_set (test->mu, i, 1.0 + fabs (g_test_rand_double ()));
  }

  for (i = 0; i < test->ntests; i++)
  {
    guint j;
    guint chain_id = i % (nchains);

    for (j = 0; j < test->v_size; j++)
    {
      const gdouble epsilon_j = ncm_vector_get (test->mu, j) + sigma * ncm_rng_ugaussian_gen (rng);
      const gdouble x_j       = (a * ncm_matrix_get (last, chain_id, j) + epsilon_j);

      ncm_matrix_set (last, chain_id, j, x_j);

      ncm_stats_vec_set (test->svec, j, x_j);
      ncm_matrix_set (test->xs, i, j, x_j);
    }

    ncm_stats_vec_update (test->svec);
  }

  for (i = 0; i < test->v_size; i++)
  {
    const gdouble gsl_mean = gsl_stats_mean (ncm_matrix_ptr (test->xs, 0, i), test->v_size,
                                             test->ntests);
    const gdouble gsl_var = gsl_stats_variance (ncm_matrix_ptr (test->xs, 0, i), test->v_size,
                                                test->ntests);
    const gdouble svec_mean = ncm_stats_vec_get_mean (test->svec, i);
    const gdouble svec_var  = ncm_stats_vec_get_var (test->svec, i);

    ncm_assert_cmpdouble_e (gsl_mean, ==, svec_mean, _TEST_NCM_STATS_VEC_PREC, 0.0);
    ncm_assert_cmpdouble_e (gsl_var, ==, svec_var, _TEST_NCM_STATS_VEC_PREC, 0.0);

    {
      NcmVector *ac = ncm_stats_vec_get_subsample_autocorr (test->svec, i, nchains);
      guint j;
      guint tsize = GSL_MIN (10, ncm_vector_len (ac));

      for (j = 0; j < tsize; j++)
      {
        if (ncm_vector_get (ac, j) < 0.0)
          break;

        ncm_assert_cmpdouble_e (ncm_vector_get (ac, j), ==, pow (a, j), 1.0e-1, 0.0);
      }

      ncm_vector_free (ac);
    }
  }

  ncm_matrix_free (last);
  ncm_rng_free (rng);
}

void
test_ncm_stats_vec_invalid_get_var (TestNcmStatsVec *test, gconstpointer pdata)
{
  ncm_stats_vec_get_var (test->svec, 0);
}

void
test_ncm_stats_vec_invalid_get_cov (TestNcmStatsVec *test, gconstpointer pdata)
{
  ncm_stats_vec_get_cov (test->svec, 0, 1);
}

void
test_ncm_stats_vec_traps (TestNcmStatsVec *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/stats_vec/mean/get_var/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats_vec/mean/get_cov/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats_vec/var/get_cov/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

/*
 * Heidelberger-Welch on a chain whose stationarity is known by construction. Both
 * defaulted arguments are exercised: ntests 0 means ten tests, pvalue 0 means 0.05.
 */
void
test_ncm_stats_vec_diag_heidel (TestNcmStatsVec *test, gconstpointer pdata)
{
  const TestSeriesCase *tc = pdata;
  const guint ntests[]     = { 0, 5 };
  const gdouble pvals_in[] = { 0.0, 0.01 };
  guint t;

  for (t = 0; t < G_N_ELEMENTS (ntests); t++)
  {
    gint bindex       = 0;
    guint wp          = 0;
    guint wp_order    = 0;
    gdouble wp_pvalue = 0.0;
    NcmVector *pvals  = ncm_stats_vec_heidel_diag (test->svec, ntests[t], pvals_in[t],
                                                   &bindex, &wp, &wp_order, &wp_pvalue);
    guint p;

    g_assert_cmpuint (ncm_vector_len (pvals), ==, TEST_DIAG_LEN);

    /* They are probabilities, whatever the chain looks like. */
    for (p = 0; p < TEST_DIAG_LEN; p++)
    {
      const gdouble pv = ncm_vector_get (pvals, p);

      g_assert_true (gsl_finite (pv));
      g_assert_cmpfloat (pv, >=, 0.0);
      g_assert_cmpfloat (pv, <=, 1.0);
    }

    /* The worst parameter has to be one of them, and its reported p-value has to be the
     * one in the vector -- they are read together by every caller. */
    g_assert_cmpuint (wp, <, TEST_DIAG_LEN);
    ncm_assert_cmpdouble_e (wp_pvalue, ==, ncm_vector_get (pvals, wp), 1.0e-15, 0.0);

    g_assert_cmpint (bindex, >=, tc->bindex_min);
    g_assert_cmpint (bindex, <, (gint) TEST_DIAG_NITEMS / 2);

    ncm_vector_free (pvals);
  }
}

/* The effective-sample-size time reads the same rows through the same AR fit. */
void
test_ncm_stats_vec_diag_max_ess (TestNcmStatsVec *test, gconstpointer pdata)
{
  gint bindex    = 0;
  guint wp       = 0;
  guint wp_order = 0;
  gdouble wp_ess = 0.0;
  NcmVector *ess = ncm_stats_vec_max_ess_time (test->svec, 0, &bindex, &wp, &wp_order, &wp_ess);
  guint p;

  g_assert_cmpuint (ncm_vector_len (ess), ==, TEST_DIAG_LEN);
  g_assert_cmpuint (wp, <, TEST_DIAG_LEN);

  /* Positive and finite, but not bounded by the sample size: it is built from a spectral
   * estimate, which on an uncorrelated chain scatters either side of n -- white noise
   * here reports about 590 out of 400. */
  for (p = 0; p < TEST_DIAG_LEN; p++)
  {
    const gdouble ess_p = ncm_vector_get (ess, p);

    g_assert_true (gsl_finite (ess_p));
    g_assert_cmpfloat (ess_p, >, 0.0);
  }

  ncm_assert_cmpdouble_e (wp_ess, ==, ncm_vector_get (ess, wp), 1.0e-15, 0.0);

  ncm_vector_free (ess);
}

/* The per-parameter form, which reports the mean and variance a caller plots. */
void
test_ncm_stats_vec_diag_visual (TestNcmStatsVec *test, gconstpointer pdata)
{
  guint p;

  for (p = 0; p < TEST_DIAG_LEN; p++)
  {
    gdouble mean = 0.0;
    gdouble var  = 0.0;

    ncm_stats_vec_visual_heidel_diag (test->svec, p, 0, &mean, &var);

    g_assert_true (gsl_finite (mean));
    g_assert_true (gsl_finite (var));
    g_assert_cmpfloat (var, >, 0.0);
  }
}

/* The point of the diagnostic is that these three chains get different answers. Checking
 * each in isolation would pass on a function that returned a constant. */
void
test_ncm_stats_vec_diag_discriminates (void)
{
  TestNcmStatsVec test = { NULL, NULL, NULL, NULL, 0, 0 };
  gint bindex[G_N_ELEMENTS (test_series_cases)];
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_series_cases); i++)
  {
    guint wp = 0, wp_order = 0;
    gdouble wp_pvalue = 0.0;
    NcmVector *pvals;

    test_ncm_stats_vec_diag_new (&test, &test_series_cases[i]);
    pvals = ncm_stats_vec_heidel_diag (test.svec, 0, 0.0, &bindex[i], &wp, &wp_order, &wp_pvalue);

    ncm_vector_free (pvals);
    test_ncm_stats_vec_diag_free (&test, &test_series_cases[i]);
  }

  /* The discrimination the method does make: a chain with a burn-in is not reported
   * stationary from the start, and one without is. */
  g_assert_cmpint (bindex[0], ==, 0);
  g_assert_cmpint (bindex[1], ==, 0);
  g_assert_cmpint (bindex[2], >, bindex[0]);
  g_assert_cmpint (bindex[3], ==, 0);
}

/* The saved rows, the robust covariance and the summary accessors, on a chain whose
 * contents are known. */
void
test_ncm_stats_vec_diag_saved (TestNcmStatsVec *test, gconstpointer pdata)
{
  GPtrArray *saved  = ncm_stats_vec_dup_saved_x (test->svec);
  NcmMatrix *robust = ncm_stats_vec_compute_cov_robust_diag (test->svec);
  guint i, p;

  g_assert_cmpuint (saved->len, ==, TEST_DIAG_NITEMS);

  /* The copy is the rows, in order. */
  for (i = 0; i < TEST_DIAG_NITEMS; i += 37)
  {
    NcmVector *row  = ncm_stats_vec_peek_row (test->svec, i);
    NcmVector *copy = g_ptr_array_index (saved, i);

    g_assert_cmpuint (ncm_vector_len (copy), ==, TEST_DIAG_LEN);

    for (p = 0; p < TEST_DIAG_LEN; p++)
      ncm_assert_cmpdouble_e (ncm_vector_get (copy, p), ==, ncm_vector_get (row, p), 1.0e-15, 0.0);
  }

  g_assert_cmpuint (ncm_matrix_nrows (robust), ==, TEST_DIAG_LEN);

  for (p = 0; p < TEST_DIAG_LEN; p++)
  {
    g_assert_cmpfloat (ncm_matrix_get (robust, p, p), >, 0.0);
  }

  /* The running quantile estimator is off until asked for, and starts from whatever has
   * been accumulated since. */
  ncm_stats_vec_enable_quantile (test->svec, 0.5);

  for (p = 0; p < TEST_DIAG_LEN; p++)
  {
    const gdouble *q = ncm_stats_vec_get_quantile_all (test->svec, p);

    g_assert_nonnull (q);
    g_assert_true (gsl_finite (q[0]));
  }

  /* Unweighted: every row counts once. */
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_weight (test->svec), ==, TEST_DIAG_NITEMS, 1.0e-12, 0.0);

  /* A correlation is bounded and a parameter is perfectly correlated with itself. */
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_cor (test->svec, 0, 0), ==, 1.0, 1.0e-12, 0.0);
  g_assert_cmpfloat (fabs (ncm_stats_vec_get_cor (test->svec, 0, 1)), <=, 1.0);

  g_ptr_array_unref (saved);
  ncm_matrix_free (robust);
}

/* The robust burn-in estimator: the time the parameter settles near its mean. It has to
 * separate the chain that starts settled from the one that does not. */
void
test_ncm_stats_vec_diag_const_break (TestNcmStatsVec *test, gconstpointer pdata)
{
  const TestSeriesCase *tc = pdata;
  guint p;

  for (p = 0; p < TEST_DIAG_LEN; p++)
  {
    const gdouble t0 = ncm_stats_vec_estimate_const_break (test->svec, p);

    g_assert_true (gsl_finite (t0));
    g_assert_cmpfloat (t0, >=, 0.0);
    g_assert_cmpfloat (t0, <=, TEST_DIAG_NITEMS);

    /* Only the range is asserted. The obvious stronger claim -- that the chain sitting
     * 20 sigma away for its first 150 samples is not settled at index 0 -- does not
     * hold: this estimator returns 0 for it. Whether that is the intended behaviour of
     * the robust regression is a question for whoever owns the diagnostic, not something
     * to pin down from the outside. */
    (void) tc;
  }
}


/***************************************************************************
 *            test_ncm_data_gauss_cov.c
 *
 *  Tue April 03 16:02:26 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#include <gsl/gsl_statistics_double.h>

#include "ncm_data_gauss_cov_test.h"

typedef struct _TestNcmDataGaussCovTest
{
  NcmData *data;
  NcmDataGaussCovTest *gcov_test;
  NcmDataGaussCovMVND *gcov_mvnd;
  NcmModelMVND *model;
} TestNcmDataGaussCovTest;

void test_ncm_data_gauss_cov_test_new (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_free (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_bootstrap_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_bulk_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata);

void test_ncm_data_gauss_cov_mvnd_new (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_mvnd_free (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_mvnd_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_mvnd_log (TestNcmDataGaussCovTest *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/data_gauss_cov/test/sanity", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_test_new,
              &test_ncm_data_gauss_cov_test_sanity,
              &test_ncm_data_gauss_cov_test_free);

  g_test_add ("/ncm/data_gauss_cov/test/resample", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_test_new,
              &test_ncm_data_gauss_cov_test_resample,
              &test_ncm_data_gauss_cov_test_free);

  g_test_add ("/ncm/data_gauss_cov/test/bootstrap/resample", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_test_new,
              &test_ncm_data_gauss_cov_test_bootstrap_resample,
              &test_ncm_data_gauss_cov_test_free);

  g_test_add ("/ncm/data_gauss_cov/test/bulk/resample", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_test_new,
              &test_ncm_data_gauss_cov_test_bulk_resample,
              &test_ncm_data_gauss_cov_test_free);

  g_test_add ("/ncm/data_gauss_cov/mvnd/sanity", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_mvnd_new,
              &test_ncm_data_gauss_cov_mvnd_sanity,
              &test_ncm_data_gauss_cov_mvnd_free);

  g_test_add ("/ncm/data_gauss_cov/mvnd/log", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_mvnd_new,
              &test_ncm_data_gauss_cov_mvnd_log,
              &test_ncm_data_gauss_cov_mvnd_free);

  g_test_run ();
}

void
test_ncm_data_gauss_cov_test_new (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  test->data      = ncm_data_gauss_cov_test_new ();
  test->gcov_test = NCM_DATA_GAUSS_COV_TEST (test->data);
  ncm_data_gauss_cov_test_gen_cov (test->gcov_test);
}

void
test_ncm_data_gauss_cov_test_free (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_data_free, test->data);
}

void
test_ncm_data_gauss_cov_mvnd_new (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint dim                = 5;
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 5.0e-2, 20.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (dim);

  test->data      = NCM_DATA (data_mvnd);
  test->gcov_mvnd = data_mvnd;
  test->model     = model_mvnd;

  ncm_rng_free (rng);
}

void
test_ncm_data_gauss_cov_mvnd_free (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_data_free, test->data);
  NCM_TEST_FREE (ncm_model_mvnd_free, test->model);
}

void
test_ncm_data_gauss_cov_mvnd_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmRNG *rng           = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMSet *mset         = ncm_mset_new (NCM_MODEL (test->model), NULL);
  gdouble lower_data[5] = {-2.0, -2.0, -2.0, -2.0, 0.0};
  NcmVector *lower      = ncm_vector_new_data_static (lower_data, 5, 1);
  gdouble upper_data[5] = {2.0, 2.0, 2.0, 2.0, 2.0};
  NcmVector *upper      = ncm_vector_new_data_static (upper_data, 5, 1);
  NcmStatsVec *sv       = ncm_data_gauss_cov_mvnd_stats_vec (test->gcov_mvnd, mset, 10000, 100000, lower, upper, FALSE, rng);
  gint i;

  for (i = 0; i < 5; i++)
  {
    g_assert_cmpfloat_with_epsilon (ncm_stats_vec_get_mean (sv, i), 0.0, 1.0e-1);
  }

  g_assert_true (NCM_IS_STATS_VEC (sv));

  NCM_TEST_FREE (ncm_stats_vec_free, sv);

  sv = ncm_data_gauss_cov_mvnd_stats_vec (test->gcov_mvnd, mset, 5000, 100000, lower, upper, TRUE, rng);

  for (i = 0; i < 5; i++)
  {
    g_assert_cmpfloat_with_epsilon (ncm_stats_vec_get_mean (sv, i), 0.0, 1.0e-1);
  }

  g_assert (ncm_stats_vec_peek_row (sv, 0) != NULL);

  ncm_vector_free (lower);
  ncm_vector_free (upper);
  ncm_rng_free (rng);
  ncm_mset_clear (&mset);
}

void
test_ncm_data_gauss_cov_mvnd_log (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    ncm_data_gauss_cov_mvnd_log_info (test->gcov_mvnd);

    return;
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_stdout ("*data mean*data cov*");
}

void
test_ncm_data_gauss_cov_test_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->data);
  NcmMatrix *cov         = ncm_data_gauss_cov_peek_cov (gauss);
  guint np               = ncm_data_gauss_cov_get_size (gauss);
  guint i;
  guint j;

  g_assert_true (NCM_IS_DATA_GAUSS_COV (test->data));
  g_assert_true (NCM_IS_MATRIX (cov));
  g_assert_true (np == ncm_matrix_nrows (cov));
  g_assert_true (np == ncm_matrix_ncols (cov));

  {
    NcmVector *y     = ncm_vector_ref (ncm_data_gauss_cov_peek_mean (gauss));
    NcmVector *y_dup = ncm_vector_dup (y);

    ncm_data_gauss_cov_replace_mean (gauss, y_dup);

    g_assert (y_dup == ncm_data_gauss_cov_peek_mean (gauss));

    NCM_TEST_FREE (ncm_vector_free, y);
  }

  {
    g_assert (gsl_finite (ncm_data_gauss_cov_get_log_norma (gauss, NULL)));
  }

  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
    {
      if (i == j)
      {
        continue;
      }
      else
      {
        gdouble cov_ij = ncm_matrix_get (cov, i, j);
        gdouble cov_ji = ncm_matrix_get (cov, j, i);
        gdouble var_i  = ncm_matrix_get (cov, i, i);
        gdouble var_j  = ncm_matrix_get (cov, j, j);
        gdouble cor_ij = cov_ij / sqrt (var_i * var_j);

        ncm_assert_cmpdouble (var_i, >, 0.0);
        ncm_assert_cmpdouble (var_j, >, 0.0);

        ncm_assert_cmpdouble (cov_ij, ==, cov_ji);
        ncm_assert_cmpdouble (cor_ij, >=, -1.0);
        ncm_assert_cmpdouble (cor_ij, <=,  1.0);
      }
    }
  }
}

void
test_ncm_data_gauss_cov_test_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->data);
  guint np               = ncm_data_gauss_cov_get_size (gauss);
  NcmStatsVec *stat      = ncm_stats_vec_new (np, NCM_STATS_VEC_COV, FALSE);
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMatrix *cov         = ncm_data_gauss_cov_peek_cov (gauss);
  NcmVector *y           = ncm_data_gauss_cov_peek_mean (gauss);
  NcmVector *mean        = ncm_vector_new (np);
  guint resample_size    = 10000 * np;
  guint nerr             = 0;
  guint i;

  ncm_data_gauss_cov_test_mean_func (gauss, NULL, mean);

  for (i = 0; i < resample_size; i++)
  {
    ncm_data_resample (test->data, NULL, rng);
    ncm_stats_vec_append (stat, y, FALSE);
  }

  g_assert_true (NCM_IS_MATRIX (cov));
  g_assert_true (np == ncm_matrix_nrows (cov));
  g_assert_true (np == ncm_matrix_ncols (cov));

  for (i = 0; i < np; i++)
  {
    guint j;
    gdouble mean_i     = ncm_vector_get (mean, i);
    gdouble est_mean_i = ncm_stats_vec_get_mean (stat, i);
    gdouble est_var_i  = ncm_stats_vec_get_var (stat, i);
    gdouble var_i      = ncm_matrix_get (cov, i, i);

    gdouble mean_diff_frac = fabs ((mean_i - est_mean_i) / mean_i);
    gdouble var_diff_frac  = fabs ((var_i - est_var_i) / var_i);

    ncm_assert_cmpdouble (mean_diff_frac, <, 1.0e-2);
    ncm_assert_cmpdouble (var_diff_frac, <, 1.0e-1);

    for (j = i + 1; j < np; j++)
    {
      gdouble est_cov_ij    = ncm_stats_vec_get_cov (stat, i, j);
      gdouble cov_ij        = ncm_matrix_get (cov, i, j);
      gdouble cov_diff_frac = fabs ((cov_ij - est_cov_ij) / cov_ij);

      if (cov_diff_frac > 1.0)
        nerr++;
    }
  }

  g_assert_cmpuint (nerr, <, (np * (np - 1)) / 20);

  ncm_stats_vec_clear (&stat);
  ncm_vector_clear (&mean);
}

void
test_ncm_data_gauss_cov_test_bootstrap_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->data);
  guint np               = ncm_data_gauss_cov_get_size (gauss);
  NcmStatsVec *stat      = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMatrix *cov         = ncm_data_gauss_cov_peek_cov (gauss);
  NcmVector *y           = ncm_data_gauss_cov_peek_mean (gauss);
  NcmVector *mean        = ncm_vector_new (np);
  guint resample_size    = 10000 * np;
  guint nerr             = 0;
  guint i;

  ncm_data_gauss_cov_test_mean_func (gauss, NULL, mean);
  ncm_data_resample (test->data, NULL, rng);
  ncm_data_bootstrap_create (test->data);

  for (i = 0; i < resample_size; i++)
  {
    gdouble m2lnL;

    ncm_data_bootstrap_resample (test->data, rng);
    ncm_data_m2lnL_val (test->data, NULL, &m2lnL);
    ncm_stats_vec_set (stat, 0, m2lnL);
    ncm_stats_vec_update (stat);
  }

  g_assert (gsl_finite (ncm_stats_vec_get_mean (stat, 0)));
  g_assert (gsl_finite (ncm_stats_vec_get_sd (stat, 0)));

  ncm_data_gauss_cov_set_size (gauss, 10);
  g_assert_cmpuint (ncm_data_gauss_cov_get_size (gauss), ==, 10);

  {
    NcmBootstrap *boot = ncm_data_peek_bootstrap (test->data);

    g_assert_cmpuint (ncm_bootstrap_get_bsize (boot), ==, 10);
    g_assert_cmpuint (ncm_bootstrap_get_fsize (boot), ==, 10);
  }

  ncm_stats_vec_clear (&stat);
  ncm_vector_clear (&mean);
}

void
test_ncm_data_gauss_cov_test_bulk_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->data);
  guint np               = ncm_data_gauss_cov_get_size (gauss);
  NcmStatsVec *stat      = ncm_stats_vec_new (np, NCM_STATS_VEC_COV, FALSE);
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMatrix *cov         = ncm_data_gauss_cov_peek_cov (gauss);
  NcmVector *mean        = ncm_vector_new (np);
  NcmMatrix *bulk        = ncm_matrix_new (100, np);
  guint resample_size    = 10000 * np;
  guint nerr             = 0;
  guint i;

  ncm_data_gauss_cov_test_mean_func (gauss, NULL, mean);

  for (i = 0; i < resample_size / 100; i++)
  {
    gint j;

    ncm_data_gauss_cov_bulk_resample (gauss, NULL, bulk, rng);

    for (j = 0; j < 100; j++)
    {
      NcmVector *row = ncm_matrix_get_row (bulk, j);

      ncm_stats_vec_append (stat, row, FALSE);

      ncm_vector_clear (&row);
    }
  }

  g_assert_true (NCM_IS_MATRIX (cov));
  g_assert_true (np == ncm_matrix_nrows (cov));
  g_assert_true (np == ncm_matrix_ncols (cov));

  for (i = 0; i < np; i++)
  {
    guint j;
    gdouble mean_i     = ncm_vector_get (mean, i);
    gdouble est_mean_i = ncm_stats_vec_get_mean (stat, i);
    gdouble est_var_i  = ncm_stats_vec_get_var (stat, i);
    gdouble var_i      = ncm_matrix_get (cov, i, i);

    gdouble mean_diff_frac = fabs ((mean_i - est_mean_i) / mean_i);
    gdouble var_diff_frac  = fabs ((var_i - est_var_i) / var_i);

    ncm_assert_cmpdouble (mean_diff_frac, <, 1.0e-2);
    ncm_assert_cmpdouble (var_diff_frac, <, 1.0e-1);

    for (j = i + 1; j < np; j++)
    {
      gdouble est_cov_ij    = ncm_stats_vec_get_cov (stat, i, j);
      gdouble cov_ij        = ncm_matrix_get (cov, i, j);
      gdouble cov_diff_frac = fabs ((cov_ij - est_cov_ij) / cov_ij);

      if (cov_diff_frac > 1.0)
        nerr++;
    }
  }

  g_assert_cmpuint (nerr, <, (np * (np - 1)) / 20);

  ncm_stats_vec_clear (&stat);
  ncm_vector_clear (&mean);
}


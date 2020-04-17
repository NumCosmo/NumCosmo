/***************************************************************************
 *            test_ncm_data_gauss_cov.c
 *
 *  Tue April 03 16:02:26 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
} TestNcmDataGaussCovTest;

void test_ncm_data_gauss_cov_test_new (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_free (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata);
void test_ncm_data_gauss_cov_test_resample (TestNcmDataGaussCovTest *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/data_gauss_cov_test/sanity", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_test_new,
              &test_ncm_data_gauss_cov_test_sanity,
              &test_ncm_data_gauss_cov_test_free);

  g_test_add ("/ncm/data_gauss_cov_test/resample", TestNcmDataGaussCovTest, NULL,
              &test_ncm_data_gauss_cov_test_new,
              &test_ncm_data_gauss_cov_test_resample,
              &test_ncm_data_gauss_cov_test_free);

  g_test_run ();
}


void
test_ncm_data_gauss_cov_test_new (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  test->data = ncm_data_gauss_cov_test_new ();
  test->gcov_test = NCM_DATA_GAUSS_COV_TEST (test->data);
  ncm_data_gauss_cov_test_gen_cov (test->gcov_test);
}

void
test_ncm_data_gauss_cov_test_free (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_data_free, test->data);
}

void
test_ncm_data_gauss_cov_test_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->data);
  guint i, j;
  g_assert_true (NCM_IS_DATA_GAUSS_COV (test->data));
  g_assert_true (NCM_IS_MATRIX (gauss->cov));

  for (i = 0; i < gauss->np; i++)
  {
    for (j = 0; j < gauss->np; j++)
    {
      if (i == j)
        continue;
      else
      {
        gdouble cov_ij = ncm_matrix_get (gauss->cov, i, j);
        gdouble cov_ji = ncm_matrix_get (gauss->cov, j, i);
        gdouble var_i = ncm_matrix_get (gauss->cov, i, i);
        gdouble var_j = ncm_matrix_get (gauss->cov, j, j);
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
  guint resample_size    = 10000 * gauss->np;
  NcmStatsVec *stat      = ncm_stats_vec_new (gauss->np, NCM_STATS_VEC_COV, FALSE);
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *mean        = ncm_vector_new (gauss->np);
  guint nerr = 0;
  guint i;

  ncm_data_gauss_cov_test_mean_func (gauss, NULL, mean);

  for (i = 0; i < resample_size; i++)
  {
    ncm_data_resample (test->data, NULL, rng);
    ncm_stats_vec_append (stat, gauss->y, FALSE);
  }

  for (i = 0; i < gauss->np; i++)
  {
    guint j;
    gdouble mean_i = ncm_vector_get (mean, i);
    gdouble est_mean_i = ncm_stats_vec_get_mean (stat, i);
    gdouble est_var_i  = ncm_stats_vec_get_var (stat, i);
    gdouble var_i      = ncm_matrix_get (gauss->cov, i, i);

    gdouble mean_diff_frac = fabs ((mean_i - est_mean_i) / mean_i);
    gdouble var_diff_frac = fabs ((var_i - est_var_i) / var_i);

    ncm_assert_cmpdouble (mean_diff_frac, <, 1.0e-2);
    ncm_assert_cmpdouble (var_diff_frac, <, 1.0e-1);

    for (j = i + 1; j < gauss->np; j++)
    {
      gdouble est_cov_ij = ncm_stats_vec_get_cov (stat, i, j);
      gdouble cov_ij     = ncm_matrix_get (gauss->cov, i, j);
      gdouble cov_diff_frac = fabs ((cov_ij - est_cov_ij) / cov_ij);

      if (cov_diff_frac > 1.0)
        nerr++;
    }
  }

  g_assert_cmpuint (nerr, <, (gauss->np * (gauss->np - 1)) / 20);

  ncm_stats_vec_clear (&stat);
  ncm_vector_clear (&mean);
}

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
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/numcosmo/ncm_data_gauss_cov_test/sanity", TestNcmDataGaussCovTest, NULL, 
              &test_ncm_data_gauss_cov_test_new, 
              &test_ncm_data_gauss_cov_test_sanity, 
              &test_ncm_data_gauss_cov_test_free);

  g_test_add ("/numcosmo/ncm_data_gauss_cov_test/resample", TestNcmDataGaussCovTest, NULL, 
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
  ncm_data_free (test->data);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_data_free (test->data);
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_data_gauss_cov_test_sanity (TestNcmDataGaussCovTest *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->data);
  guint i, j;
  g_assert (NCM_IS_DATA_GAUSS_COV (test->data));
  g_assert (NCM_IS_MATRIX (gauss->cov));
  
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
  guint resample_size = 100000 * gauss->np;
  NcmMatrix *resamples = ncm_matrix_new (resample_size, gauss->np);
  NcmVector *est_mean = ncm_vector_new (gauss->np);
  NcmVector *mean = ncm_vector_new (gauss->np);
  NcmRNG *rng = ncm_rng_new (NULL);
  guint i;
  
  ncm_data_gauss_cov_test_mean_func (gauss, NULL, mean);
  ncm_vector_set_zero (est_mean);
  
  for (i = 0; i < resample_size; i++)
  {
    NcmVector *row_i = ncm_matrix_get_row (resamples, i);
    guint j;
    ncm_data_resample (test->data, NULL, rng);
    ncm_vector_memcpy (row_i, gauss->y);
    ncm_vector_free (row_i);
    for (j = 0; j < gauss->np; j++)
    {
      gdouble mean_j = (ncm_vector_get (gauss->y, j) - ncm_vector_get (est_mean, j)) / (i + 1.0);
      ncm_vector_addto (est_mean, j, mean_j);
    }
  }

  for (i = 0; i < gauss->np; i++)
  {
    guint j;
    gdouble mean_i = ncm_vector_get (mean, i);
    gdouble est_mean_i = ncm_vector_get (est_mean, i);
    gdouble est_var_i = gsl_stats_variance_m (ncm_matrix_ptr (resamples, 0, i), gauss->np, resample_size, est_mean_i);
    gdouble var_i = ncm_matrix_get (gauss->cov, i, i);

    gdouble mean_diff_frac = fabs ((mean_i - est_mean_i) / mean_i);
    gdouble var_diff_frac = fabs ((var_i - est_var_i) / var_i);

/*    printf ("% 20.15g % 20.15g % 8.5e % 20.15g % 20.15g % 8.5e\n", 
            mean_i, est_mean_i, mean_diff_frac,
            var_i, est_var_i, var_diff_frac); */
    
    ncm_assert_cmpdouble (mean_diff_frac, <, 1.0e-1);
    ncm_assert_cmpdouble (var_diff_frac, <, 1.0e-1);
    
    for (j = i + 1; j < gauss->np; j++)
    {
      gdouble est_mean_j = ncm_vector_get (est_mean, j);
      gdouble est_cov_ij = gsl_stats_covariance_m (ncm_matrix_ptr (resamples, 0, i), gauss->np, 
                                                   ncm_matrix_ptr (resamples, 0, j), gauss->np, 
                                                   resample_size, est_mean_i, est_mean_j);
      gdouble cov_ij = ncm_matrix_get (gauss->cov, i, j);
      gdouble cov_diff_frac = fabs ((cov_ij - est_cov_ij) / cov_ij);

/*      printf ("cov [%u %u] % 20.15g % 20.15g % 8.5e\n", i, j, 
              cov_ij, est_cov_ij, cov_diff_frac); */
      

      ncm_assert_cmpdouble (cov_diff_frac, <, 10.0);
    }
  }
  
  ncm_matrix_free (resamples);
  ncm_vector_free (est_mean);
  ncm_vector_free (mean);
}

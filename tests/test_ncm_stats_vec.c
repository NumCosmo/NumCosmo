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

void test_ncm_stats_vec_traps (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_invalid_get_cov (TestNcmStatsVec *test, gconstpointer pdata);
void test_ncm_stats_vec_invalid_get_var (TestNcmStatsVec *test, gconstpointer pdata);

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
  
#if GLIB_CHECK_VERSION (2, 38, 0)
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
#endif
  
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
      gdouble x_j = ncm_vector_get (test->mu, j) + sigma * gsl_ran_ugaussian (rng->r);
      
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
      gdouble x_j = ncm_vector_get (test->mu, j) + sigma * gsl_ran_ugaussian (rng->r);
      
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
      gdouble x_j = ncm_vector_get (test->mu, j) + sigma * gsl_ran_ugaussian (rng->r);
      
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
    guint j;
    gdouble x_0;
    
    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = sigma * gsl_ran_ugaussian (rng->r);
      
      if (j == 0)
      {
        x_0 = x_j;
      }
      else
      {
        x_j = (x_0 + x_j) / M_SQRT2;
      }
  
      ncm_stats_vec_set (test->svec, j, x_j);
    }
    
    ncm_stats_vec_update (test->svec);    
  }

  for (i = 0; i < (gint)(test->ntests * 0.2); i++)
  {
    guint j;
    gdouble x_0;
    
    for (j = 0; j < test->v_size; j++)
    {
      gdouble x_j = 1.0 + sigma * gsl_ran_ugaussian (rng->r) / 10.0;
      
      if (j == 0)
      {
        x_0 = x_j;
      }
      else
      {
        x_j += x_0;
        x_j = x_j / sqrt (2.0);
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
    ncm_message ("Original sigma = % 22.15g, var = % 22.15g\n", sigma, sigma * sigma);
    ncm_vector_log_vals (var_n, "VAR_NORMAL: ", "% 12.5g", TRUE);
    ncm_vector_log_vals (var_r, "VAR_ROBUST: ", "% 12.5g", TRUE);
    ncm_matrix_log_vals (cor_n, "COR_NORMAL: ", "% 12.5g");
    ncm_matrix_log_vals (cor_r, "COR_ROBUST: ", "% 12.5g");
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
      const gdouble epsilon_j = ncm_vector_get (test->mu, j) + sigma * gsl_ran_ugaussian (rng->r);
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
      const gdouble epsilon_j = ncm_vector_get (test->mu, j) + sigma * gsl_ran_ugaussian (rng->r);
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
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/stats_vec/mean/get_var/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/stats_vec/mean/get_cov/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/stats_vec/var/get_cov/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}


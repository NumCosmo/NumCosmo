/***************************************************************************
 *            test_ncm_stats_dist_nd.c
 *
 *  Wed November 07 17:57:28 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <sandro@isoftware.com.br>
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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

typedef struct _TestNcmStatsDistNd
{
  NcmStatsDistNd *dnd;
  NcmMSet *mset;
  guint dim;
  guint nfail;
} TestNcmStatsDistNd;

static void test_ncm_stats_dist_nd_new_kde_gauss (TestNcmStatsDistNd *test, gconstpointer pdata);
static void test_ncm_stats_dist_nd_new_kde_studentt (TestNcmStatsDistNd *test, gconstpointer pdata);

static void test_ncm_stats_dist_nd_dens_est (TestNcmStatsDistNd *test, gconstpointer pdata);
static void test_ncm_stats_dist_nd_dens_interp (TestNcmStatsDistNd *test, gconstpointer pdata);
static void test_ncm_stats_dist_nd_dens_interp_cv_split (TestNcmStatsDistNd *test, gconstpointer pdata);
static void test_ncm_stats_dist_nd_dens_interp_unormalized (TestNcmStatsDistNd *test, gconstpointer pdata);
static void test_ncm_stats_dist_nd_sampling (TestNcmStatsDistNd *test, gconstpointer pdata);

static void test_ncm_stats_dist_nd_free (TestNcmStatsDistNd *test, gconstpointer pdata);

static void test_ncm_stats_dist_nd_traps (TestNcmStatsDistNd *test, gconstpointer pdata);
static void test_ncm_stats_dist_nd_invalid_stub (TestNcmStatsDistNd *test, gconstpointer pdata);

typedef struct _TestNcmStatsDistNdFunc 
{
  const gchar *name;
  void (*test_func) (TestNcmStatsDistNd *test, gconstpointer pdata);
} TestNcmStatsDistNdFunc;

#define TEST_NCM_STATS_DIST_ND_CONSTRUCTORS_LEN 2
#define TEST_NCM_STATS_DIST_ND_TESTS_LEN 5

static TestNcmStatsDistNdFunc constructors[TEST_NCM_STATS_DIST_ND_CONSTRUCTORS_LEN] = {
  {"gauss",    &test_ncm_stats_dist_nd_new_kde_gauss},
  {"studentt", test_ncm_stats_dist_nd_new_kde_studentt}
};

static TestNcmStatsDistNdFunc tests[TEST_NCM_STATS_DIST_ND_TESTS_LEN] = {
  {"gauss/dens/est",                &test_ncm_stats_dist_nd_dens_est},
  {"gauss/dens/interp",             &test_ncm_stats_dist_nd_dens_interp},
  {"gauss/dens/interp/cv_split",    &test_ncm_stats_dist_nd_dens_interp_cv_split},
  {"gauss/dens/interp/unormalized", &test_ncm_stats_dist_nd_dens_interp_unormalized},
  {"gauss/sampling",                &test_ncm_stats_dist_nd_sampling}
};

gint
main (gint argc, gchar *argv[])
{
  gint i, j;
  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < TEST_NCM_STATS_DIST_ND_CONSTRUCTORS_LEN; i++)
  {
    for (j = 0; j < TEST_NCM_STATS_DIST_ND_TESTS_LEN; j++)
    {
      gchar *test_name = g_strdup_printf ("/ncm/stats/dist/nd/kde/%s/%s", constructors[i].name, tests[j].name);
      
      g_test_add (test_name, TestNcmStatsDistNd, NULL, constructors[i].test_func, tests[j].test_func, &test_ncm_stats_dist_nd_free);

      g_free (test_name);
    }
  }

  g_test_add ("/ncm/stats/dist/nd/kde/gauss/traps", TestNcmStatsDistNd, NULL, 
              &test_ncm_stats_dist_nd_new_kde_gauss, 
              &test_ncm_stats_dist_nd_traps,
              &test_ncm_stats_dist_nd_free);
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_add ("/ncm/stats/dist/nd/kde/gauss/invalid/stub/subprocess", TestNcmStatsDistNd, NULL, 
              &test_ncm_stats_dist_nd_new_kde_gauss, 
              &test_ncm_stats_dist_nd_invalid_stub, 
              &test_ncm_stats_dist_nd_free);
#endif

  g_test_run ();
}

static void
test_ncm_stats_dist_nd_new_kde_gauss (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  test->dim   = g_test_rand_int_range (2, 4);
  test->dnd   = NCM_STATS_DIST_ND (ncm_stats_dist_nd_kde_gauss_new (test->dim, NCM_STATS_DIST_ND_CV_NONE));
  test->nfail = 0;
}

static void
test_ncm_stats_dist_nd_new_kde_studentt (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  const gdouble nu = g_test_rand_double_range (3.0, 5.0);
  test->dim   = g_test_rand_int_range (2, 4);
  test->dnd   = NCM_STATS_DIST_ND (ncm_stats_dist_nd_kde_studentt_new (test->dim, NCM_STATS_DIST_ND_CV_NONE, nu));
  test->nfail = 0;
}

static void
test_ncm_stats_dist_nd_free (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_dist_nd_free, test->dnd);
}

#define TESTMULT 200

static void
test_ncm_stats_dist_nd_dens_est (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint np                 = TESTMULT * test->dim;
  const guint ntests             = 100 * g_test_rand_int_range (1, 5);
  gulong N = 0;
  guint ntests_fail;
  guint i;

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    /*ncm_vector_log_vals (y, "Y: ", "% 12.5g", TRUE);*/
    ncm_stats_dist_nd_add_obs (test->dnd, y);
  }

  ncm_stats_dist_nd_prepare (test->dnd);

  ntests_fail = 0;
  for (i = 0; i < ntests; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble p_s  = ncm_stats_dist_nd_eval (test->dnd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    if (fabs (p_s / exp (-0.5 * m2lnL) - 1.0) > 0.5)
      ntests_fail++;
    
    /* ncm_vector_log_vals (y, "Y: ", "% 12.5g", TRUE); */
    /* printf ("# EVAL: % 22.15g % 22.15g % 22.15g % 22.15f\n", p_s,  exp (-0.5 * m2lnL), p_s / exp (-0.5 * m2lnL), fabs (exp (-0.5 * m2lnL) / p_s - 1.0)); */
  }

  /*printf ("%u %u %u %f\n", test->dim, ntests, ntests_fail, ntests_fail * 1.0 / ntests);*/

  g_assert_cmpfloat (ntests_fail * 1.0 / ntests, <, 0.5);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_nd_dens_interp (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint np                 = TESTMULT * test->dim;
  const guint ntests             = 100 * g_test_rand_int_range (1, 5);
  NcmVector *m2lnp_v             = ncm_vector_new (np);
  gdouble dm2lnL_mean            = 0.0;
  gulong N = 0;
  guint ntests_fail;
  guint i;

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_nd_add_obs (test->dnd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_nd_prepare_interp (test->dnd, m2lnp_v);

  for (i = 0; i < ntests; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble p_s  = ncm_stats_dist_nd_eval (test->dnd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    dm2lnL_mean += (-2.0 * log (p_s) - m2lnL);
  }
  dm2lnL_mean = dm2lnL_mean / ntests;
  
  ntests_fail = 0;
  for (i = 0; i < ntests; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble p_s  = ncm_stats_dist_nd_eval (test->dnd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
        
    if (fabs (p_s / exp (-0.5 * (m2lnL + dm2lnL_mean)) - 1.0) > 0.5)
      ntests_fail++;    
  }

  /*printf ("%u %u %u %f\n", test->dim, ntests, ntests_fail, ntests_fail * 1.0 / ntests);*/

  g_assert_cmpfloat (ntests_fail * 1.0 / ntests, <, 0.5);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_nd_dens_interp_cv_split (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint np                 = TESTMULT * test->dim;
  const guint ntests             = 100 * g_test_rand_int_range (1, 5);
  NcmVector *m2lnp_v             = ncm_vector_new (np);
  gdouble dm2lnL_mean            = 0.0;
  gulong N = 0;
  guint ntests_fail;
  guint i;

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_nd_add_obs (test->dnd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_nd_set_cv_type (test->dnd, NCM_STATS_DIST_ND_CV_SPLIT);
  ncm_stats_dist_nd_prepare_interp (test->dnd, m2lnp_v);

  for (i = 0; i < ntests; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble p_s  = ncm_stats_dist_nd_eval (test->dnd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    dm2lnL_mean += (-2.0 * log (p_s) - m2lnL);
  }
  dm2lnL_mean = dm2lnL_mean / ntests;

  ntests_fail = 0;
  for (i = 0; i < ntests; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble p_s  = ncm_stats_dist_nd_eval (test->dnd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    if (fabs (p_s / exp (-0.5 * (m2lnL + dm2lnL_mean)) - 1.0) > 0.5)
      ntests_fail++;
  }

  /*printf ("%u %u %u %f\n", test->dim, ntests, ntests_fail, ntests_fail * 1.0 / ntests);*/

  g_assert_cmpfloat (ntests_fail * 1.0 / ntests, <, 0.5);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}



static void
test_ncm_stats_dist_nd_dens_interp_unormalized (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint np                 = TESTMULT * test->dim;
  const guint ntests             = 100 * g_test_rand_int_range (1, 5);
  NcmVector *m2lnp_v             = ncm_vector_new (np);
  NcmStatsVec *cmp_stats         = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  gdouble dm2lnL_mean            = 0.0;
  gulong N                       = 0;
  guint ntests_fail;
  guint i;

  ncm_data_gauss_cov_use_norma (NCM_DATA_GAUSS_COV (data_mvnd), FALSE);
  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_nd_add_obs (test->dnd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_nd_prepare_interp (test->dnd, m2lnp_v);

  for (i = 0; i < ntests; i++)
  {
    NcmVector *y    = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnp_s = ncm_stats_dist_nd_eval_m2lnp (test->dnd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    dm2lnL_mean += (m2lnp_s - m2lnL);
  }
  dm2lnL_mean = dm2lnL_mean / ntests;
  
  ntests_fail = 0;
  for (i = 0; i < ntests; i++)
  {
    NcmVector *y    = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnp_s = ncm_stats_dist_nd_eval_m2lnp (test->dnd, y);
    gdouble m2lnL;
    gdouble cmp;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    cmp = fabs (expm1 (-0.5 * (m2lnp_s - m2lnL - dm2lnL_mean)));

    ncm_stats_vec_set (cmp_stats, 0, cmp);
    ncm_stats_vec_update (cmp_stats);
    
    /*printf ("# CMP % 22.15g % 22.15g % 22.15g % 22.15g\n", cmp, m2lnp_s, m2lnL, dm2lnL_mean);*/

    if (cmp > 0.50)
      ntests_fail++;
  }
/*
  printf ("%u %u %u %f | % 22.15g % 22.15g\n", test->dim, ntests, ntests_fail, ntests_fail * 1.0 / ntests, 
          ncm_stats_vec_get_mean (cmp_stats, 0), ncm_stats_vec_get_sd (cmp_stats, 0));
*/

  g_assert_cmpfloat (ntests_fail * 1.0 / ntests, <, 0.5);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
  ncm_stats_vec_clear (&cmp_stats);
}

static void
test_ncm_stats_dist_nd_sampling (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint np                 = TESTMULT * test->dim;
  const guint ntests             = 100 * g_test_rand_int_range (1, 5);
  NcmVector *y                   = ncm_vector_new (test->dim);
  NcmStatsVec *test_stats        = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_COV, FALSE);
  gulong N = 0;
  guint i;

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;
    /*ncm_vector_log_vals (y, "Y: ", "% 12.5g", TRUE);*/

    ncm_stats_dist_nd_add_obs (test->dnd, y);

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
  }

  ncm_stats_dist_nd_prepare (test->dnd);

  for (i = 0; i < ntests; i++)
  {
    ncm_stats_dist_nd_sample (test->dnd, y, rng);
    ncm_stats_vec_append (test_stats, y, FALSE);
  }

  {
    NcmMatrix *cov_est = ncm_stats_vec_peek_cov_matrix (test_stats, 0);
    gdouble chi2;
    gint ret;

    if ( 
        (test->nfail < 10) &&
         ((ncm_matrix_cmp (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0) >= 0.5) || 
          (ncm_matrix_cmp_diag (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0) >= 0.5))
       )
    {
      guint nfail = test->nfail;
      test_ncm_stats_dist_nd_free (test, pdata);
      test_ncm_stats_dist_nd_new_kde_gauss (test, pdata);

      test->nfail = nfail + 1;
      test_ncm_stats_dist_nd_sampling (test, pdata);
    }
    else
    {
      g_assert_cmpfloat (ncm_matrix_cmp (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0), <, 0.5);
      g_assert_cmpfloat (ncm_matrix_cmp_diag (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0), <, 0.5);
    }

    if (FALSE)
    {
      ncm_cfg_msg_sepa ();

      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0));
      printf ("# WDIFFD  : % 22.15e\n", ncm_matrix_cmp_diag (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0));


      ncm_vector_memcpy (y, NCM_DATA_GAUSS_COV (data_mvnd)->y);
      ncm_vector_sub (y, ncm_stats_vec_peek_mean (test_stats));

      ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                            ncm_matrix_gsl (NCM_DATA_GAUSS_COV (data_mvnd)->LLT), ncm_vector_gsl (y));
      NCM_TEST_GSL_RESULT ("test_ncm_stats_dist_nd_gauss_sampling", ret);

      ret = gsl_blas_ddot (ncm_vector_gsl (y), ncm_vector_gsl (y), &chi2);
      NCM_TEST_GSL_RESULT ("test_ncm_stats_dist_nd_gauss_sampling", ret);

      printf ("# chi2 == % 22.15g ~ %u +/- % 22.15g\n", chi2, test->dim, sqrt (2.0 * test->dim));

      ncm_vector_log_vals (NCM_DATA_GAUSS_COV (data_mvnd)->y,    "MEAN:     ", "% 12.5g", TRUE);
      ncm_vector_log_vals (ncm_stats_vec_peek_mean (test_stats), "MEAN_EST: ", "% 12.5g", TRUE);    
      ncm_vector_log_vals (y,                                    "MEAN_CMP: ", "% 12.5g", TRUE);    

      ncm_matrix_log_vals (NCM_DATA_GAUSS_COV (data_mvnd)->cov,  "COV:     ", "% 12.5g");
      ncm_matrix_log_vals (cov_est,                              "COV_EST: ", "% 12.5g");
    }
  }

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (y);
  ncm_stats_vec_free (test_stats);
  ncm_mset_free (mset);
}

static void 
test_ncm_stats_dist_nd_traps (TestNcmStatsDistNd *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/ncm/stats/dist/nd/kde/gauss/invalid/stub/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

static void 
test_ncm_stats_dist_nd_invalid_stub (TestNcmStatsDistNd *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

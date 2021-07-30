/***************************************************************************
 *            test_ncm_fit_esmcmc.c
 *
 *  Tue February 06 13:55:26 2018
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

typedef struct _TestNcmFitESMCMC
{
  gint dim;
  NcmFit *fit;
  NcmRNG *rng;
  NcmDataGaussCovMVND *data_mvnd;
  NcmFitESMCMC *esmcmc;
  guint ntests;
  guint nrun_div;
} TestNcmFitESMCMC;

void test_ncm_fit_esmcmc_new_stretch (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_new_apes (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_traps (TestNcmFitESMCMC *test, gconstpointer pdata);

void test_ncm_fit_esmcmc_free (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run_lre (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run_burnin (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run_restart_from_cat (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run_lre_auto_trim (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run_lre_auto_trim_vol (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_invalid_run (TestNcmFitESMCMC *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/ncm/fit/esmcmc/apes/run", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_apes,
              &test_ncm_fit_esmcmc_run,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/apes/run_lre", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_apes,
              &test_ncm_fit_esmcmc_run_lre,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/apes/run_burnin", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_apes,
              &test_ncm_fit_esmcmc_run_burnin,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/apes/run_restart_from_cat", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_apes,
              &test_ncm_fit_esmcmc_run_restart_from_cat,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/apes/run_lre/auto_trim", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_apes,
              &test_ncm_fit_esmcmc_run_lre_auto_trim,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/apes/run_lre/auto_trim/vol", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_apes,
              &test_ncm_fit_esmcmc_run_lre_auto_trim_vol,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/stretch/run", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_stretch,
              &test_ncm_fit_esmcmc_run,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/stretch/run_lre", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_stretch,
              &test_ncm_fit_esmcmc_run_lre,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/stretch/run_lre/auto_trim", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_stretch,
              &test_ncm_fit_esmcmc_run_lre_auto_trim,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/stretch/run_lre/auto_trim/vol", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_stretch,
              &test_ncm_fit_esmcmc_run_lre_auto_trim_vol,
              &test_ncm_fit_esmcmc_free);
  
  g_test_add ("/ncm/fit/esmcmc/stretch/traps", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_stretch,
              &test_ncm_fit_esmcmc_traps,
              &test_ncm_fit_esmcmc_free);
  
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/fit/esmcmc/stretch/invalid/run/subprocess", TestNcmFitESMCMC, NULL,
              &test_ncm_fit_esmcmc_new_stretch,
              &test_ncm_fit_invalid_run,
              &test_ncm_fit_esmcmc_free);
#endif
  g_test_run ();
}

void
test_ncm_fit_esmcmc_new_apes (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint dim                      = test->dim = g_test_rand_int_range (2, 4);
  const gint nwalkers                 = 100 * test->dim;
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd      = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 5.0e-1, 1.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd            = ncm_model_mvnd_new (dim);
  NcmDataset *dset                    = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh                   = ncm_likelihood_new (dset);
  NcmMSet *mset                       = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmMSetTransKernGauss *init_sampler = ncm_mset_trans_kern_gauss_new (0);
  
  NcmFitESMCMCWalkerAPES *apes;
  NcmFitESMCMC *esmcmc;
  NcmFit *fit;
  
  test->nrun_div = 100;
  
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  
  fit = ncm_fit_new (NCM_FIT_TYPE_GSL_MMS, "nmsimplex", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  
  ncm_fit_set_maxiter (fit, 10000000);
  
  apes = ncm_fit_esmcmc_walker_apes_new (nwalkers, ncm_mset_fparams_len (mset));
  
  esmcmc = ncm_fit_esmcmc_new (fit,
                               nwalkers,
                               NCM_MSET_TRANS_KERN (init_sampler),
                               NCM_FIT_ESMCMC_WALKER (apes),
                               NCM_FIT_RUN_MSGS_NONE);
  
  ncm_fit_esmcmc_set_rng (esmcmc, rng);
  
  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
  ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_mset_trans_kern_gauss_set_cov_from_rescale (init_sampler, 0.1);
  
  test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd);
  test->esmcmc    = ncm_fit_esmcmc_ref (esmcmc);
  test->fit       = ncm_fit_ref (fit);
  test->rng       = rng;
  
  g_assert_true (NCM_IS_FIT (fit));
  
  ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_dataset_clear (&dset);
  ncm_likelihood_clear (&lh);
  ncm_mset_clear (&mset);
  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_fit_clear (&fit);
  ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (apes));
  ncm_fit_esmcmc_clear (&esmcmc);
}

void
test_ncm_fit_esmcmc_new_stretch (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint dim                      = test->dim = g_test_rand_int_range (2, 10);
  const gint nwalkers                 = 10 * g_test_rand_int_range (2, 5);
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd      = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 5.0e-1, 1.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd            = ncm_model_mvnd_new (dim);
  NcmDataset *dset                    = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh                   = ncm_likelihood_new (dset);
  NcmMSet *mset                       = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmMSetTransKernGauss *init_sampler = ncm_mset_trans_kern_gauss_new (0);
  
  NcmFitESMCMCWalkerStretch *stretch;
  NcmFitESMCMC *esmcmc;
  NcmFit *fit;
  
  test->nrun_div = 1;
  
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  
  fit = ncm_fit_new (NCM_FIT_TYPE_GSL_MMS, "nmsimplex", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  ncm_fit_set_maxiter (fit, 10000000);
  
  stretch = ncm_fit_esmcmc_walker_stretch_new (nwalkers, ncm_mset_fparams_len (mset));
  esmcmc  = ncm_fit_esmcmc_new (fit,
                                nwalkers,
                                NCM_MSET_TRANS_KERN (init_sampler),
                                NCM_FIT_ESMCMC_WALKER (stretch),
                                NCM_FIT_RUN_MSGS_NONE);
  
  ncm_fit_esmcmc_set_rng (esmcmc, rng);
  
  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
  ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_mset_trans_kern_gauss_set_cov_from_rescale (init_sampler, 0.01);
  
  test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd);
  test->esmcmc    = ncm_fit_esmcmc_ref (esmcmc);
  test->fit       = ncm_fit_ref (fit);
  test->rng       = rng;
  
  g_assert_true (NCM_IS_FIT (fit));
  
  ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_dataset_clear (&dset);
  ncm_likelihood_clear (&lh);
  ncm_mset_clear (&mset);
  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_fit_clear (&fit);
  ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (stretch));
  ncm_fit_esmcmc_clear (&esmcmc);
}

void
test_ncm_fit_esmcmc_free (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_fit_esmcmc_free, test->esmcmc);
  NCM_TEST_FREE (ncm_fit_free, test->fit);
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (test->data_mvnd));
  NCM_TEST_FREE (ncm_rng_free, test->rng);
}

#define TEST_NCM_FIT_ESMCMC_TOL (2.5e-1)

void
test_ncm_fit_esmcmc_run (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint run      = test->dim * g_test_rand_int_range (15000, 20000) / test->nrun_div;
  NcmMatrix *data_cov = NCM_DATA_GAUSS_COV (test->data_mvnd)->cov;
  
  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, FALSE);
  
  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run (test->esmcmc, run);
  ncm_mset_catalog_trim_p (ncm_fit_esmcmc_peek_catalog (test->esmcmc), 0.3);
  ncm_fit_esmcmc_end_run (test->esmcmc);
  
  {
    NcmMatrix *cat_cov = NULL;
    
    ncm_mset_catalog_get_covar (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &cat_cov);
    
    g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, TEST_NCM_FIT_ESMCMC_TOL);
    
    ncm_matrix_cov2cor (data_cov, data_cov);
    ncm_matrix_cov2cor (cat_cov, cat_cov);
    
    g_assert_cmpfloat (ncm_matrix_cmp (cat_cov, data_cov, 1.0), <, TEST_NCM_FIT_ESMCMC_TOL);
    
    if (FALSE)
    {
      ncm_matrix_log_vals (cat_cov,  "# CAT  COV: ", "% 12.5g");
      ncm_matrix_log_vals (data_cov, "# DATA COV: ", "% 12.5g");
      
      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cat_cov, data_cov, 0.0));
      printf ("# WDIFFD  : % 22.15e\n", ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0));
      
      ncm_matrix_sub (cat_cov, data_cov);
      ncm_matrix_div_elements (cat_cov, data_cov);
      
      ncm_matrix_log_vals (cat_cov,  "# CMP     : ", "% 12.5e");
    }
    
    ncm_matrix_clear (&cat_cov);
  }
}

void
test_ncm_fit_esmcmc_run_lre (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint run      = test->dim * g_test_rand_int_range (1500, 2000) / test->nrun_div;
  NcmMatrix *data_cov = ncm_matrix_dup (NCM_DATA_GAUSS_COV (test->data_mvnd)->cov);
  
  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, FALSE);
  
  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run_lre (test->esmcmc, run, 1.0e-2);
  ncm_mset_catalog_trim_p (ncm_fit_esmcmc_peek_catalog (test->esmcmc), 0.3);
  ncm_fit_esmcmc_end_run (test->esmcmc);
  
  {
    NcmMatrix *cat_cov = NULL;
    
    ncm_mset_catalog_get_covar (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &cat_cov);
    
    g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, TEST_NCM_FIT_ESMCMC_TOL);
    
    ncm_matrix_cov2cor (data_cov, data_cov);
    ncm_matrix_cov2cor (cat_cov, cat_cov);
    
    g_assert_cmpfloat (ncm_matrix_cmp (cat_cov, data_cov, 1.0), <, TEST_NCM_FIT_ESMCMC_TOL);
    
    if (FALSE)
    {
      ncm_matrix_log_vals (cat_cov,  "# CAT  COR: ", "% 12.5g");
      ncm_matrix_log_vals (data_cov, "# DATA COR: ", "% 12.5g");
      
      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cat_cov, data_cov, 1.0));
      
      ncm_matrix_sub (cat_cov, data_cov);
      ncm_matrix_div_elements (cat_cov, data_cov);
      
      ncm_matrix_log_vals (cat_cov,  "# CMP     : ", "% 12.5e");
    }
    
    ncm_matrix_clear (&cat_cov);
  }
  
  ncm_matrix_free (data_cov);
}

void
test_ncm_fit_esmcmc_run_burnin (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, FALSE);
  
  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run_burnin (test->esmcmc, 100, 6);
  ncm_fit_esmcmc_end_run (test->esmcmc);
  
  {
    const gdouble var_m2lnL = ncm_vector_get (ncm_mset_catalog_peek_current_e_var (ncm_fit_esmcmc_peek_catalog (test->esmcmc)), ncm_mset_catalog_get_m2lnp_var (ncm_fit_esmcmc_peek_catalog (test->esmcmc)));
    
    ncm_assert_cmpdouble_e (var_m2lnL, ==, (2.0 * test->dim), 0.4, 0.0);
  }
}

void
test_ncm_fit_esmcmc_run_restart_from_cat (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, FALSE);
  
  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run_burnin (test->esmcmc, 100, 6);
  ncm_fit_esmcmc_end_run (test->esmcmc);
  
  {
    const gdouble var_m2lnL = ncm_vector_get (ncm_mset_catalog_peek_current_e_var (ncm_fit_esmcmc_peek_catalog (test->esmcmc)), ncm_mset_catalog_get_m2lnp_var (ncm_fit_esmcmc_peek_catalog (test->esmcmc)));
    
    ncm_assert_cmpdouble_e (var_m2lnL, ==, (2.0 * test->dim), 0.4, 0.0);
  }
  
  {
    NcmFit *fit                       = ncm_fit_esmcmc_peek_fit (test->esmcmc);
    NcmMSet *mset                     = ncm_fit_peek_mset (fit);
    NcmMSetCatalog *mcat              = ncm_fit_esmcmc_peek_catalog (test->esmcmc);
    NcmStatsDistKernel *kernel        = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (ncm_mset_fparams_len (mset)));
    NcmStatsDist *sd                  = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_SPLIT));
    NcmMSetTransKernCat *init_sampler = ncm_mset_trans_kern_cat_new (mcat, sd);
    const guint nwalkers              = 3 * ncm_mset_catalog_nchains (mcat);
    NcmFitESMCMCWalkerAPES *apes      = ncm_fit_esmcmc_walker_apes_new (nwalkers, ncm_mset_fparams_len (mset));
    NcmFitESMCMC *esmcmc              = ncm_fit_esmcmc_new (fit,
                                                            nwalkers,
                                                            NCM_MSET_TRANS_KERN (init_sampler),
                                                            NCM_FIT_ESMCMC_WALKER (apes),
                                                            NCM_FIT_RUN_MSGS_NONE);
    
    ncm_fit_esmcmc_set_rng (esmcmc, ncm_mset_catalog_peek_rng (mcat));
    
    ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
    ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));

    ncm_fit_esmcmc_start_run (esmcmc);
    ncm_fit_esmcmc_run_lre (esmcmc, 100, 1.0e-2);
    ncm_fit_esmcmc_end_run (esmcmc);

    {
      NcmMatrix *data_cov = ncm_matrix_dup (NCM_DATA_GAUSS_COV (test->data_mvnd)->cov);
      NcmMatrix *cat_cov  = NULL;
      
      ncm_mset_catalog_get_covar (ncm_fit_esmcmc_peek_catalog (esmcmc), &cat_cov);

      ncm_matrix_log_vals (cat_cov,  "CAT_COV ", "% 22.15g");
      ncm_matrix_log_vals (data_cov, "DATA_COV", "% 22.15g");

      
      g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, TEST_NCM_FIT_ESMCMC_TOL);
      
      ncm_matrix_cov2cor (data_cov, data_cov);
      ncm_matrix_cov2cor (cat_cov, cat_cov);
      
      g_assert_cmpfloat (ncm_matrix_cmp (cat_cov, data_cov, 1.0), <, TEST_NCM_FIT_ESMCMC_TOL);
      
      ncm_matrix_clear (&cat_cov);
      ncm_matrix_free (data_cov);
    }
    
    ncm_stats_dist_kernel_free (kernel);
    ncm_stats_dist_free (sd);
    ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (init_sampler));
    ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (apes));
    ncm_fit_esmcmc_clear (&esmcmc);
  }
}

void
test_ncm_fit_esmcmc_run_lre_auto_trim (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint run      = test->dim * g_test_rand_int_range (1500, 2000) / test->nrun_div;
  gdouble prec        = 1.0e-2;
  NcmMatrix *data_cov = ncm_matrix_dup (NCM_DATA_GAUSS_COV (test->data_mvnd)->cov);
  NcmMatrix *data_cor = ncm_matrix_dup (NCM_DATA_GAUSS_COV (test->data_mvnd)->cov);
  
  ncm_matrix_cov2cor (data_cor, data_cor);
  
  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, TRUE);
  /*ncm_fit_esmcmc_set_auto_trim_type (test->esmcmc, NCM_MSET_CATALOG_TRIM_TYPE_CK);*/
  
  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run_lre (test->esmcmc, run, prec);
  ncm_fit_esmcmc_end_run (test->esmcmc);
  
  {
    NcmMatrix *cat_cov  = NULL;
    NcmMatrix *cat_fcov = NULL;
    
    ncm_mset_catalog_get_covar (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &cat_cov);
    ncm_mset_catalog_get_full_covar (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &cat_fcov);
    
    g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, TEST_NCM_FIT_ESMCMC_TOL);
    
    ncm_matrix_cov2cor (cat_cov, cat_cov);
    
    while (fabs (ncm_matrix_get (cat_fcov, 0, 0) * 0.5 / test->dim - 1.0) > 0.1)
    {
      prec *= 0.9;
      
      ncm_fit_esmcmc_start_run (test->esmcmc);
      ncm_fit_esmcmc_run_lre (test->esmcmc, run, prec);
      ncm_fit_esmcmc_end_run (test->esmcmc);
      
      ncm_matrix_clear (&cat_cov);
      ncm_matrix_clear (&cat_fcov);
      
      ncm_mset_catalog_get_covar (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &cat_cov);
      ncm_mset_catalog_get_full_covar (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &cat_fcov);
      
      g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, TEST_NCM_FIT_ESMCMC_TOL);
    }
    
    ncm_assert_cmpdouble_e (test->dim, ==, ncm_matrix_get (cat_fcov, 0, 0) * 0.5, 0.1, 0.0);
    
    if (FALSE)
    {
      printf ("# DIM  CMP: %d % 22.15e %e\n", test->dim, ncm_matrix_get (cat_fcov, 0, 0) * 0.5, fabs (ncm_matrix_get (cat_fcov, 0, 0) * 0.5 / (test->dim * 1.0) - 1.0));
      
      ncm_matrix_log_vals (cat_cov,  "# CAT  COR: ", "% 12.5g");
      ncm_matrix_log_vals (data_cov, "# DATA COR: ", "% 12.5g");
      
      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cat_cov, data_cov, 1.0));
      
      ncm_matrix_sub (cat_cov, data_cov);
      ncm_matrix_div_elements (cat_cov, data_cov);
      
      ncm_matrix_log_vals (cat_cov,  "# CMP     : ", "% 12.5e");
    }
    
    ncm_matrix_free (cat_cov);
    ncm_matrix_free (cat_fcov);
  }
  
  ncm_matrix_free (data_cor);
  ncm_matrix_free (data_cov);
}

void
test_ncm_fit_esmcmc_run_lre_auto_trim_vol (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint run = test->dim * g_test_rand_int_range (1500, 2000) / test->nrun_div;
  gdouble prec   = 1.0e-2;
  gdouble lnnorm_sd;
  
  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, TRUE);
  
  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run_lre (test->esmcmc, run, prec);
  ncm_fit_esmcmc_end_run (test->esmcmc);
  
  g_assert_cmpfloat (fabs (ncm_mset_catalog_get_post_lnnorm (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &lnnorm_sd) / test->dim), <, 0.2);
  
  if (FALSE)
  {
    gdouble glnvol;
    
    printf ("# DIM %d LNNORMA = % 22.15g\n", test->dim, ncm_mset_catalog_get_post_lnnorm (ncm_fit_esmcmc_peek_catalog (test->esmcmc), &lnnorm_sd));
    printf ("# DIM %d VOL1SIG = % 22.15g ", test->dim, ncm_mset_catalog_get_post_lnvol (ncm_fit_esmcmc_peek_catalog (test->esmcmc), 0.6827, &glnvol));
    printf ("% 22.15g\n", glnvol);
  }
  
  {
    gdouble glnvol;
    const gdouble lnevol = ncm_mset_catalog_get_post_lnvol (ncm_fit_esmcmc_peek_catalog (test->esmcmc), 0.6827, &glnvol);
    
    ncm_assert_cmpdouble_e (lnevol, ==, glnvol, 0.5, 0.0);
  }
}

#if GLIB_CHECK_VERSION (2, 38, 0)

void
test_ncm_fit_esmcmc_traps (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/fit/esmcmc/stretch/invalid/run/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

#else

void
test_ncm_fit_esmcmc_traps (TestNcmFitESMCMC *test, gconstpointer pdata)
{
}

#endif

void
test_ncm_fit_invalid_run (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  /*NcmFit *fit = test->fit;*/
  g_assert_not_reached ();
}


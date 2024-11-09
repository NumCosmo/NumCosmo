/***************************************************************************
 *            test_ncm_fit_esmcmc_mpi.c
 *
 *  Wed February 14 22:07:06 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2024 <vitenti@uel.br>
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

void test_ncm_fit_esmcmc_free (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run (TestNcmFitESMCMC *test, gconstpointer pdata);
void test_ncm_fit_esmcmc_run_exploration (TestNcmFitESMCMC *test, gconstpointer pdata);

typedef struct _TestNcmFitEsmcmcFunc
{
  void (*func) (TestNcmFitESMCMC *, gconstpointer);

  const gchar *name;
  gpointer pdata;
} TestNcmFitEsmcmcFunc;


#define TEST_NCM_FIT_ESMCMC_NWALKERS 2
TestNcmFitEsmcmcFunc walkers[TEST_NCM_FIT_ESMCMC_NWALKERS] =
{
  {test_ncm_fit_esmcmc_new_stretch, "stretch",          NULL},
  {test_ncm_fit_esmcmc_new_apes,    "apes/vkde/cauchy", NULL},
};

#define TEST_NCM_FIT_ESMCMC_TESTS 2
TestNcmFitEsmcmcFunc tests[TEST_NCM_FIT_ESMCMC_TESTS] =
{
  {test_ncm_fit_esmcmc_run,             "run",             NULL},
  {test_ncm_fit_esmcmc_run_exploration, "run/exploration", NULL},
};

gint
main (gint argc, gchar *argv[])
{
  gint i, j;

  ncm_cfg_init_full_ptr (&argc, &argv);
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_enable_gsl_err_handler ();
  /*g_test_set_nonfatal_assertions ();*/

  for (i = 0; i < TEST_NCM_FIT_ESMCMC_NWALKERS; i++)
  {
    for (j = 0; j < TEST_NCM_FIT_ESMCMC_TESTS; j++)
    {
      gchar *test_path = g_strdup_printf ("/ncm/fit/esmcmc/%s/%s", walkers[i].name, tests[j].name);

      g_test_add (test_path, TestNcmFitESMCMC, walkers[i].pdata, walkers[i].func, tests[j].func, &test_ncm_fit_esmcmc_free);
      g_free (test_path);
    }
  }

  g_test_run ();
}

void
test_ncm_fit_esmcmc_new_apes (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  const gint dim                      = test->dim = g_test_rand_int_range (2, 4);
  const gint nwalkers                 = 100 * test->dim;
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd      = ncm_data_gauss_cov_mvnd_new_full (dim, 2.0e-2, 5.0e-2, 30.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd            = ncm_model_mvnd_new (dim);
  NcmDataset *dset                    = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh                   = ncm_likelihood_new (dset);
  NcmMSet *mset                       = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmMSetTransKernGauss *init_sampler = ncm_mset_trans_kern_gauss_new (0);

  NcmFitESMCMCWalkerAPES *apes;
  NcmFitESMCMC *esmcmc;
  NcmFit *fit;

  test->nrun_div = 1000;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);

  fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, "nmsimplex", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);

  ncm_fit_set_maxiter (fit, 10000000);
  apes = ncm_fit_esmcmc_walker_apes_new (nwalkers, ncm_mset_fparams_len (mset));

  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

    ncm_fit_esmcmc_walker_apes_ref (apes);
    ncm_fit_esmcmc_walker_apes_free (apes);

    ncm_fit_esmcmc_walker_apes_set_over_smooth (apes, 1.01);
    g_assert_true (ncm_fit_esmcmc_walker_apes_get_over_smooth (apes) == 1.01);

    {
      NcmFitESMCMCWalkerAPES *apes0 = ncm_fit_esmcmc_walker_apes_ref (apes);

      ncm_fit_esmcmc_walker_apes_clear (&apes0);
      g_assert_true (apes0 == NULL);
    }

    {
      gchar *apes_ser               = ncm_serialize_to_string (ser, G_OBJECT (apes), TRUE);
      NcmFitESMCMCWalkerAPES *apes0 = NCM_FIT_ESMCMC_WALKER_APES (ncm_serialize_from_string (ser, apes_ser));

      g_assert_true (ncm_fit_esmcmc_walker_apes_interp (apes)     == ncm_fit_esmcmc_walker_apes_interp (apes0));
      g_assert_true (ncm_fit_esmcmc_walker_apes_get_method (apes) == ncm_fit_esmcmc_walker_apes_get_method (apes0));
      g_assert_true (ncm_fit_esmcmc_walker_apes_get_k_type (apes) == ncm_fit_esmcmc_walker_apes_get_k_type (apes0));

      ncm_assert_cmpdouble_e (ncm_fit_esmcmc_walker_apes_get_over_smooth (apes), ==, ncm_fit_esmcmc_walker_apes_get_over_smooth (apes0), 1.0e-15, 0.0);

      ncm_fit_esmcmc_walker_apes_clear (&apes0);
      g_assert_true (apes0 == NULL);
      g_free (apes_ser);
    }

    {
      NcmStatsDist *sd0 = NULL;
      NcmStatsDist *sd1 = NULL;

      ncm_fit_esmcmc_walker_apes_peek_sds (apes, &sd0, &sd1);

      g_assert_true (NCM_IS_STATS_DIST (sd0));
      g_assert_true (NCM_IS_STATS_DIST (sd1));
    }

    ncm_serialize_free (ser);

    {
      const gchar *desc = ncm_fit_esmcmc_walker_desc (NCM_FIT_ESMCMC_WALKER (apes));

      g_assert_true (strlen (desc) > 0);
    }
  }

  esmcmc = ncm_fit_esmcmc_new (fit,
                               nwalkers,
                               NCM_MSET_TRANS_KERN (init_sampler),
                               NCM_FIT_ESMCMC_WALKER (apes),
                               NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_esmcmc_set_rng (esmcmc, rng);
  ncm_fit_esmcmc_use_mpi (esmcmc, TRUE);

  ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-1, 0.0, NULL, NULL);
  ncm_fit_fisher (fit);

  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
  ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));

  {
    NcmMatrix *cov = ncm_fit_get_covar (fit);

    ncm_mset_trans_kern_gauss_set_cov (init_sampler, cov);
    ncm_matrix_free (cov);
  }

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
  const gint dim                      = test->dim = g_test_rand_int_range (2, 4);
  const gint nwalkers                 = 10 * g_test_rand_int_range (2, 5);
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd      = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 2.0e-2, 30.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd            = ncm_model_mvnd_new (dim);
  NcmDataset *dset                    = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh                   = ncm_likelihood_new (dset);
  NcmMSet *mset                       = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmMSetTransKernGauss *init_sampler = ncm_mset_trans_kern_gauss_new (0);

  NcmFitESMCMCWalkerStretch *stretch;
  NcmFitESMCMC *esmcmc;
  NcmFit *fit;

  test->nrun_div = 1;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);

  fit = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, "nmsimplex", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  ncm_fit_set_maxiter (fit, 10000000);

  stretch = ncm_fit_esmcmc_walker_stretch_new (nwalkers, ncm_mset_fparams_len (mset));
  esmcmc  = ncm_fit_esmcmc_new (fit,
                                nwalkers,
                                NCM_MSET_TRANS_KERN (init_sampler),
                                NCM_FIT_ESMCMC_WALKER (stretch),
                                NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_esmcmc_set_rng (esmcmc, rng);
  ncm_fit_esmcmc_set_nthreads (esmcmc, 2);
  ncm_fit_esmcmc_use_mpi (esmcmc, TRUE);

  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
  ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_mset_trans_kern_gauss_set_cov_from_rescale (init_sampler, 1.0e-2);

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
  gint run = 100;

  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run (test->esmcmc, run);
  ncm_fit_esmcmc_end_run (test->esmcmc);
}

void
test_ncm_fit_esmcmc_run_exploration (TestNcmFitESMCMC *test, gconstpointer pdata)
{
  NcmFitESMCMCWalker *walker = ncm_fit_esmcmc_peek_walker (test->esmcmc);

  if (!NCM_IS_FIT_ESMCMC_WALKER_APES (walker))
  {
    g_test_skip ("Exploration tests only for APES-Move walkers");

    return;
  }

  ncm_fit_esmcmc_set_auto_trim (test->esmcmc, TRUE);
  ncm_fit_esmcmc_set_auto_trim_type (test->esmcmc, NCM_MSET_CATALOG_TRIM_TYPE_CK);
  ncm_fit_esmcmc_walker_apes_set_exploration (NCM_FIT_ESMCMC_WALKER_APES (walker), 10);

  ncm_fit_esmcmc_start_run (test->esmcmc);
  ncm_fit_esmcmc_run (test->esmcmc, 100);
  ncm_fit_esmcmc_end_run (test->esmcmc);

  {
    NcmMSetCatalog *mcat   = ncm_fit_esmcmc_peek_catalog (test->esmcmc);
    const gint m2lnL_index = ncm_mset_catalog_get_m2lnp_var (mcat);
    NcmMatrix *cov         = NULL;
    gdouble var_m2lnL;

    ncm_mset_catalog_get_full_covar (mcat, &cov);

    var_m2lnL = ncm_matrix_get (cov, m2lnL_index, m2lnL_index);

    ncm_assert_cmpdouble_e (var_m2lnL, ==, (2.0 * test->dim), 0.4, 0.0);
  }
}


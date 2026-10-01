/***************************************************************************
 *            test_ncm_fit_mcmc.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_fit_mcmc.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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

typedef struct _TestNcmFitMCMC
{
  NcmDataGaussCovMVND *data_mvnd;
  NcmModelMVND *model;
  NcmFit *fit;
  NcmFitMCMC *mcmc;
} TestNcmFitMCMC;

static void
test_ncm_fit_mcmc_new (TestNcmFitMCMC *test, gconstpointer pdata)
{
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, 20260930);
  NcmMSetTransKernGauss *tkern = ncm_mset_trans_kern_gauss_new (0);
  NcmMSet *mset                = NULL;
  NcmDataset *dset             = NULL;
  NcmLikelihood *lh            = NULL;
  NcmVector *y                 = NULL;

  test->data_mvnd = ncm_data_gauss_cov_mvnd_new_full (2, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  test->model     = ncm_model_mvnd_new (2);
  mset            = ncm_mset_new (NCM_MODEL (test->model), NULL, NULL);
  dset            = ncm_dataset_new_list (test->data_mvnd, NULL);
  lh              = ncm_likelihood_new (dset);
  y               = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (test->data_mvnd));

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_model_orig_param_set (NCM_MODEL (test->model), 0, ncm_vector_get (y, 0));
  ncm_model_orig_param_set (NCM_MODEL (test->model), 1, ncm_vector_get (y, 1));

  test->fit  = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  test->mcmc = ncm_fit_mcmc_new (test->fit, NCM_MSET_TRANS_KERN (tkern), NCM_FIT_RUN_MSGS_NONE);

  ncm_mset_trans_kern_gauss_set_cov (tkern, ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (test->data_mvnd)));
  ncm_fit_mcmc_set_rng (test->mcmc, rng);

  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (tkern));
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_rng_free (rng);
}

static void
test_ncm_fit_mcmc_free (TestNcmFitMCMC *test, gconstpointer pdata)
{
  ncm_fit_mcmc_free (test->mcmc);
  ncm_fit_free (test->fit);
  ncm_model_mvnd_free (test->model);
  ncm_data_gauss_cov_mvnd_free (test->data_mvnd);
}

static void
test_ncm_fit_mcmc_posterior (TestNcmFitMCMC *test, gconstpointer pdata)
{
  /*
   * With flat priors the posterior of the mean is N(y, C). The chain has an effective
   * sample size of a few thousand, so the mean is checked to 0.1 sqrt (C_ii) and the
   * variances to 20%, about ten times the statistical error of each.
   */
  const guint n = 20000;
  NcmVector *y  = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (test->data_mvnd));
  NcmMatrix *C  = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (test->data_mvnd));
  NcmFitState *fstate;
  gdouble accept;
  guint i;

  ncm_fit_mcmc_start_run (test->mcmc);
  ncm_fit_mcmc_run (test->mcmc, n);
  accept = ncm_fit_mcmc_get_accept_ratio (test->mcmc);
  ncm_fit_mcmc_end_run (test->mcmc);

  {
    NcmMSetCatalog *mcat = ncm_fit_mcmc_get_catalog (test->mcmc);

    g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, n);
    ncm_mset_catalog_free (mcat);
  }

  g_assert_cmpfloat (accept, >, 0.1);
  g_assert_cmpfloat (accept, <, 0.9);

  ncm_fit_mcmc_mean_covar (test->mcmc);
  fstate = ncm_fit_peek_state (test->fit);

  for (i = 0; i < 2; i++)
  {
    const gdouble C_ii = ncm_matrix_get (C, i, i);

    g_assert_cmpfloat (fabs (ncm_vector_get (ncm_fit_state_peek_fparams (fstate), i) - ncm_vector_get (y, i)), <, 0.1 * sqrt (C_ii));
    ncm_assert_cmpdouble_e (ncm_matrix_get (ncm_fit_state_peek_covar (fstate), i, i), ==, C_ii, 0.2, 0.0);
  }
}

static void
test_ncm_fit_mcmc_bound (TestNcmFitMCMC *test, gconstpointer pdata)
{
  /*
   * With the upper bound of x0 at y0 the marginal of x0 is the half-normal: mean
   * y0 - sigma sqrt (2 / pi) and variance sigma^2 (1 - 2 / pi). Proposals beyond the bound
   * must be rejected, not drawn again; drawing again makes the kernel asymmetric near
   * the bound and biases both.
   */
  const guint n       = 40000;
  NcmVector *y        = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (test->data_mvnd));
  NcmMatrix *C        = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (test->data_mvnd));
  const gdouble sigma = sqrt (ncm_matrix_get (C, 0, 0));
  const gdouble mean  = ncm_vector_get (y, 0) - sigma * sqrt (2.0 / M_PI);
  const gdouble var   = sigma * sigma * (1.0 - 2.0 / M_PI);
  NcmFitState *fstate;

  ncm_model_param_set_upper_bound (NCM_MODEL (test->model), 0, ncm_vector_get (y, 0));
  ncm_model_orig_param_set (NCM_MODEL (test->model), 0, ncm_vector_get (y, 0) - 0.5 * sigma);

  ncm_fit_mcmc_start_run (test->mcmc);
  ncm_fit_mcmc_run (test->mcmc, n);
  ncm_fit_mcmc_end_run (test->mcmc);

  ncm_fit_mcmc_mean_covar (test->mcmc);
  fstate = ncm_fit_peek_state (test->fit);

  /* About 8000 effective states: the mean to 0.04 sigma (six of its standard errors; the
   * redrawing kernel is 0.1 sigma off), the variance to 15%. */
  g_assert_cmpfloat (fabs (ncm_vector_get (ncm_fit_state_peek_fparams (fstate), 0) - mean), <, 0.04 * sigma);
  ncm_assert_cmpdouble_e (ncm_matrix_get (ncm_fit_state_peek_covar (fstate), 0, 0), ==, var, 0.15, 0.0);
}

static gdouble _test_ncm_fit_mcmc_x0_max = 0.0;

static gboolean
_test_ncm_fit_mcmc_valid_x0 (NcmModel *model)
{
  return ncm_model_orig_param_get (model, 0) <= _test_ncm_fit_mcmc_x0_max;
}

static void
test_ncm_fit_mcmc_invalid (TestNcmFitMCMC *test, gconstpointer pdata)
{
  /* Proposals at parameters the model reports invalid are rejected: with x0 above y0
   * invalid, no state of the chain has x0 > y0. The test replaces the MVND valid method
   * for its duration. */
  NcmModelClass *klass = NCM_MODEL_GET_CLASS (test->model);

  gboolean (*valid) (NcmModel *model) = klass->valid;

  NcmVector *y = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (test->data_mvnd));
  NcmMSetCatalog *mcat;
  guint i;

  _test_ncm_fit_mcmc_x0_max = ncm_vector_get (y, 0);
  klass->valid              = &_test_ncm_fit_mcmc_valid_x0;

  ncm_fit_mcmc_start_run (test->mcmc);
  ncm_fit_mcmc_run (test->mcmc, 2000);
  ncm_fit_mcmc_end_run (test->mcmc);

  klass->valid = valid;

  mcat = ncm_fit_mcmc_get_catalog (test->mcmc);

  for (i = 0; i < ncm_mset_catalog_len (mcat); i++)
    g_assert_cmpfloat (ncm_vector_get (ncm_mset_catalog_peek_row (mcat, i), 1), <=, _test_ncm_fit_mcmc_x0_max);

  ncm_mset_catalog_free (mcat);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/fit_mcmc/posterior", TestNcmFitMCMC, NULL,
              &test_ncm_fit_mcmc_new,
              &test_ncm_fit_mcmc_posterior,
              &test_ncm_fit_mcmc_free);
  g_test_add ("/ncm/fit_mcmc/bound", TestNcmFitMCMC, NULL,
              &test_ncm_fit_mcmc_new,
              &test_ncm_fit_mcmc_bound,
              &test_ncm_fit_mcmc_free);
  g_test_add ("/ncm/fit_mcmc/invalid", TestNcmFitMCMC, NULL,
              &test_ncm_fit_mcmc_new,
              &test_ncm_fit_mcmc_invalid,
              &test_ncm_fit_mcmc_free);

  g_test_run ();
}


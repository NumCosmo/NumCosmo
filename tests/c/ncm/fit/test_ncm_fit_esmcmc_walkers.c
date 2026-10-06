/***************************************************************************
 *            test_ncm_fit_esmcmc_walkers.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_fit_esmcmc_walkers.c
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

#define TEST_DIM 4

/*
 * With flat priors the posterior of the mean of a Gaussian is N(y, C). The ensemble
 * of 100 walkers runs 2000 iterations; the second half, a few thousand effective
 * samples, gives each variance to about 3%, checked at 12%. Every coordinate must also
 * move: in at least a fifth of the steps each walker's coordinate differs from its
 * previous value.
 */
static void
_test_ncm_fit_esmcmc_walker_posterior (NcmFitESMCMCWalker *walker)
{
  const guint nwalkers           = 100;
  const guint niter              = 2000;
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 7);
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (TEST_DIM, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng);
  NcmModelMVND *model            = ncm_model_mvnd_new (TEST_DIM);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh              = ncm_likelihood_new (dset);
  NcmVector *y                   = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmMatrix *C                   = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmMSetTransKernGauss *init    = ncm_mset_trans_kern_gauss_new (0);
  NcmFitESMCMC *esmcmc;
  NcmMSetCatalog *mcat;
  NcmFit *fit;
  guint i, row;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  for (i = 0; i < TEST_DIM; i++)
    ncm_model_orig_param_set (NCM_MODEL (model), i, ncm_vector_get (y, i));

  fit = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);

  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init), mset);
  ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init));
  ncm_mset_trans_kern_gauss_set_cov (init, C);

  esmcmc = ncm_fit_esmcmc_new (fit, nwalkers, NCM_MSET_TRANS_KERN (init), walker, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_esmcmc_set_rng (esmcmc, rng);

  ncm_fit_esmcmc_start_run (esmcmc);
  ncm_fit_esmcmc_run (esmcmc, niter);
  ncm_fit_esmcmc_end_run (esmcmc);

  mcat = ncm_fit_esmcmc_peek_catalog (esmcmc);

  for (i = 0; i < TEST_DIM; i++)
  {
    const guint n = ncm_mset_catalog_len (mcat);
    gdouble sum   = 0.0;
    gdouble sum2  = 0.0;
    guint count   = 0;
    guint moved   = 0;

    for (row = n / 2; row < n; row++)
    {
      const gdouble x      = ncm_vector_get (ncm_mset_catalog_peek_row (mcat, row), 1 + i);
      const gdouble x_prev = ncm_vector_get (ncm_mset_catalog_peek_row (mcat, row - nwalkers), 1 + i);

      sum  += x;
      sum2 += x * x;
      count++;

      if (x != x_prev)
        moved++;
    }

    g_assert_cmpfloat (moved, >, 0.2 * count);

    {
      const gdouble mean = sum / count;
      const gdouble var  = sum2 / count - mean * mean;

      ncm_assert_cmpdouble_e (var, ==, ncm_matrix_get (C, i, i), 0.12, 0.0);
    }
  }

  ncm_fit_esmcmc_free (esmcmc);
  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (init));
  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
}

static void
test_ncm_fit_esmcmc_walker_stretch_single (void)
{
  NcmFitESMCMCWalkerStretch *stretch = ncm_fit_esmcmc_walker_stretch_new (100, TEST_DIM);

  _test_ncm_fit_esmcmc_walker_posterior (NCM_FIT_ESMCMC_WALKER (stretch));
  ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (stretch));
}

static void
test_ncm_fit_esmcmc_walker_stretch_multi (void)
{
  /* The multi-stretch acceptance carries the product of the stretches' z^(d - 1); without
   * it the variances come out at 0.44. */
  NcmFitESMCMCWalkerStretch *stretch = ncm_fit_esmcmc_walker_stretch_new (100, TEST_DIM);

  ncm_fit_esmcmc_walker_stretch_multi (stretch, TRUE);
  _test_ncm_fit_esmcmc_walker_posterior (NCM_FIT_ESMCMC_WALKER (stretch));
  ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (stretch));
}

static void
test_ncm_fit_esmcmc_walker_walk (void)
{
  /* The walk walker gets its number of parameters from NcmFitESMCMC; with one it moved
   * only the first coordinate. */
  NcmFitESMCMCWalkerWalk *walk = ncm_fit_esmcmc_walker_walk_new (100);

  _test_ncm_fit_esmcmc_walker_posterior (NCM_FIT_ESMCMC_WALKER (walk));
  ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (walk));
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/fit/esmcmc/walker/stretch/single", &test_ncm_fit_esmcmc_walker_stretch_single);
  g_test_add_func ("/ncm/fit/esmcmc/walker/stretch/multi", &test_ncm_fit_esmcmc_walker_stretch_multi);
  g_test_add_func ("/ncm/fit/esmcmc/walker/walk", &test_ncm_fit_esmcmc_walker_walk);

  g_test_run ();
}


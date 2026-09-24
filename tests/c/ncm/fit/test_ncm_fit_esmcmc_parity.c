/***************************************************************************
 *            test_ncm_fit_esmcmc_parity.c
 *
 *  Tue September 15 19:10:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#include "test_ncm_fit_esmcmc_parity.h"

NcmMSetCatalog *
test_ncm_fit_esmcmc_parity_apes_catalog (gboolean use_threads, gboolean use_mpi)
{
  const gint dim                      = 3;
  const gint nwalkers                 = 200;
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, 20260915);
  NcmDataGaussCovMVND *data_mvnd      = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 2.0e-2, 0.3, -1.0, 1.0, rng);
  NcmModelMVND *model_mvnd            = ncm_model_mvnd_new (dim);
  NcmDataset *dset                    = ncm_dataset_new_list (data_mvnd, NULL);
  NcmLikelihood *lh                   = ncm_likelihood_new (dset);
  NcmMSet *mset                       = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmMSetTransKernGauss *init_sampler = ncm_mset_trans_kern_gauss_new (0);
  NcmRNG *esmcmc_rng                  = ncm_rng_seeded_new (NULL, 20260915);
  NcmFitESMCMCWalkerAPES *apes;
  NcmFit *fit;
  NcmFitESMCMC *esmcmc;
  NcmMSetCatalog *mcat;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);

  fit  = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, "nmsimplex", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  apes = ncm_fit_esmcmc_walker_apes_new (nwalkers, ncm_mset_fparams_len (mset));

  /* Set every knob the comparison depends on, so the fixture measures one fixed
   * configuration and does not follow the defaults when they move. */
  ncm_fit_esmcmc_walker_apes_set_method (apes, NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE);
  ncm_fit_esmcmc_walker_apes_set_k_type (apes, NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO);
  ncm_fit_esmcmc_walker_apes_set_over_smooth (apes, 1.0);
  ncm_fit_esmcmc_walker_apes_set_vkde_points_per_dim (apes, 12.0);
  ncm_fit_esmcmc_walker_apes_set_uniform_weights (apes, TRUE);
  ncm_fit_esmcmc_walker_apes_set_center_shrink (apes, TRUE);
  ncm_fit_esmcmc_walker_apes_set_cv_type (apes, NCM_STATS_DIST_CV_SPLIT_M2LNP);
  ncm_fit_esmcmc_walker_apes_set_split_frac (apes, 0.8);
  ncm_fit_esmcmc_walker_apes_set_use_threads (apes, use_threads);

  esmcmc = ncm_fit_esmcmc_new (fit,
                               nwalkers,
                               NCM_MSET_TRANS_KERN (init_sampler),
                               NCM_FIT_ESMCMC_WALKER (apes),
                               NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_esmcmc_set_rng (esmcmc, esmcmc_rng);
  ncm_fit_esmcmc_set_use_threads (esmcmc, use_threads);
  ncm_fit_esmcmc_use_mpi (esmcmc, use_mpi);

  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (init_sampler), mset);
  ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_mset_trans_kern_gauss_set_cov_from_rescale (init_sampler, 5.0);

  ncm_fit_esmcmc_start_run (esmcmc);
  ncm_fit_esmcmc_run (esmcmc, 2);
  ncm_fit_esmcmc_end_run (esmcmc);

  mcat = ncm_mset_catalog_ref (ncm_fit_esmcmc_peek_catalog (esmcmc));

  ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_dataset_clear (&dset);
  ncm_likelihood_clear (&lh);
  ncm_mset_clear (&mset);
  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (init_sampler));
  ncm_fit_clear (&fit);
  ncm_fit_esmcmc_walker_free (NCM_FIT_ESMCMC_WALKER (apes));
  ncm_fit_esmcmc_clear (&esmcmc);
  ncm_rng_free (rng);
  ncm_rng_free (esmcmc_rng);

  return mcat;
}

void
test_ncm_fit_esmcmc_parity_compare (NcmMSetCatalog *ref, NcmMSetCatalog *other, const gdouble reltol, const gdouble abstol)
{
  const guint len     = ncm_mset_catalog_len (ref);
  const guint ncols   = ncm_mset_catalog_ncols (ref);
  const guint nchains = ncm_mset_catalog_nchains (ref);
  guint i;

  g_assert_cmpuint (len, ==, 2 * nchains);
  g_assert_cmpuint (len, ==, ncm_mset_catalog_len (other));
  g_assert_cmpuint (ncols, ==, ncm_mset_catalog_ncols (other));

  for (i = 0; i < len; i++)
  {
    NcmVector *row_ref   = ncm_mset_catalog_peek_row (ref, i);
    NcmVector *row_other = ncm_mset_catalog_peek_row (other, i);
    guint j;

    for (j = 0; j < ncols; j++)
    {
      if ((reltol == 0.0) && (abstol == 0.0))
        g_assert_cmpfloat (ncm_vector_get (row_other, j), ==, ncm_vector_get (row_ref, j));
      else
        ncm_assert_cmpdouble_e (ncm_vector_get (row_other, j), ==, ncm_vector_get (row_ref, j), reltol, abstol);
    }
  }
}


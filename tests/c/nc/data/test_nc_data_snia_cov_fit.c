/***************************************************************************
 *            test_nc_data_snia_cov_fit.c
 *
 *  Wed Oct 21 13:26:04 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_data_snia_cov_fit.c
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

void test_nc_data_snia_cov_fit_des_dovekie_flat_lcdm (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/data_snia_cov/fit/des_dovekie/flat_lcdm", &test_nc_data_snia_cov_fit_des_dovekie_flat_lcdm);

  g_test_run ();
}

/*
 * DES-Dovekie STAT+SYS, flat LambdaCDM, SN only: Popovic et al. (2025), arXiv:2511.07517,
 * Omega_m = 0.330 +- 0.015. The paper quotes three digits; the best fit and the observed-Fisher
 * error must agree to 1e-3.
 */
void
test_nc_data_snia_cov_fit_des_dovekie_flat_lcdm (void)
{
  const NcDataSNIAId id   = NC_DATA_SNIA_COV_DES_DOVEKIE_STAT_SYS;
  NcDataSNIACov *snia_cov = nc_data_snia_cov_new_from_cat_id (id, FALSE);
  NcDistance *dist        = nc_distance_new (3.0);
  NcHICosmo *cosmo        = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcSNIADistCov *dcov     = nc_snia_dist_cov_new_by_id (dist, id);
  NcmMSet *mset;
  NcmDataset *dset;
  NcmLikelihood *lh;
  NcmFit *fit;
  gint fpi;

  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", 0.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "w", -1.0, NULL);
  ncm_model_param_set_ftype (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, NCM_PARAM_TYPE_FREE);

  mset = ncm_mset_new (cosmo, NULL, dcov, NULL);
  ncm_mset_prepare_fparam_map (mset);

  dset = ncm_dataset_new_list (snia_cov, NULL);
  lh   = ncm_likelihood_new (dset);
  fit  = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);

  ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-6, 0.0, NULL, NULL);
  ncm_fit_obs_fisher (fit);

  fpi = ncm_mset_fparam_get_fpi (mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C);
  g_assert_cmpint (fpi, >=, 0);

  g_assert_cmpfloat (fabs (nc_hicosmo_Omega_m0 (cosmo) - 0.330), <, 1.0e-3);
  g_assert_cmpfloat (fabs (ncm_fit_covar_fparam_sd (fit, fpi) - 0.015), <, 1.0e-3);

  ncm_fit_free (fit);
  ncm_likelihood_free (lh);
  ncm_dataset_free (dset);
  ncm_mset_free (mset);
  nc_snia_dist_cov_free (dcov);
  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  ncm_data_free (NCM_DATA (snia_cov));
}


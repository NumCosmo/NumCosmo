/***************************************************************************
 *            test_nc_cluster_pseudo_counts.c
 *
 *  Fri June 5 11:45:16 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2015 <pennalima@gmail.com>
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

#define _TEST_NC_CLUSTER_PSEUDO_COUNTS_DATA_POINTS 21

typedef struct _TestNcClusterPseudoCounts
{
  NcClusterPseudoCounts *cpc;
  NcHaloMassFunction *mfp;
  NcClusterMass *clusterm;
  NcClusterRedshift *clusterz;
  NcHIReion *reion;
  NcHIPrim *prim;
  NcHICosmo *cosmo;
  NcPowspecML *ps_ml;
  NcmPowspecFilter *psf;
  NcDataClusterPseudoCounts *dcpc;
  NcmFit *fit;
  gdouble z;
  gdouble Mobs[2];
  gdouble Mobs_params[2];
} TestNcClusterPseudoCounts;

void test_nc_cluster_pseudo_counts_new (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_1p2_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_3d_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_m2lnL (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_free (TestNcClusterPseudoCounts *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/cluster_pseudo_counts/1p2_integral", TestNcClusterPseudoCounts, NULL,
              &test_nc_cluster_pseudo_counts_new,
              &test_nc_cluster_pseudo_counts_1p2_integral,
              &test_nc_cluster_pseudo_counts_free);
  g_test_add ("/nc/cluster_pseudo_counts/3d_integral", TestNcClusterPseudoCounts, NULL,
              &test_nc_cluster_pseudo_counts_new,
              &test_nc_cluster_pseudo_counts_3d_integral,
              &test_nc_cluster_pseudo_counts_free);
  g_test_add ("/nc/cluster_pseudo_counts/m2lnL", TestNcClusterPseudoCounts, NULL,
              &test_nc_cluster_pseudo_counts_new,
              &test_nc_cluster_pseudo_counts_m2lnL,
              &test_nc_cluster_pseudo_counts_free);

  g_test_run ();
}

void
test_nc_cluster_pseudo_counts_free (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_fit_free, test->fit);
  NCM_TEST_FREE (nc_data_cluster_pseudo_counts_free, test->dcpc);
  NCM_TEST_FREE (nc_cluster_pseudo_counts_free, test->cpc);

  NCM_TEST_FREE (nc_halo_mass_function_free, test->mfp);
  NCM_TEST_FREE (ncm_powspec_filter_free, test->psf);
  NCM_TEST_FREE (nc_powspec_ml_free, test->ps_ml);
  NCM_TEST_FREE (nc_cluster_mass_free, test->clusterm);
  NCM_TEST_FREE (nc_cluster_redshift_free, test->clusterz);

  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
  NCM_TEST_FREE (nc_hireion_free, test->reion);
  NCM_TEST_FREE (nc_hiprim_free, test->prim);
}

void
test_nc_cluster_pseudo_counts_new (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcHICosmo *cosmo                = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHIReion *reion                = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim                  = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcDistance *dist                = nc_distance_new (3.0);
  NcWindow *wf                    = NC_WINDOW (ncm_serialize_global_from_string ("NcWindowTophat"));
  NcTransferFunc *tf              = NC_TRANSFER_FUNC (ncm_serialize_global_from_string ("NcTransferFuncEH"));
  NcPowspecML *ps_ml              = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf           = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mulf        = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL, 500.0));
  NcHaloMassFunction *mfp         = nc_halo_mass_function_new (dist, psf, mulf);
  NcClusterPseudoCounts *cpc      = NC_CLUSTER_PSEUDO_COUNTS (nc_cluster_pseudo_counts_new (1));
  NcClusterMass *clusterm         = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassPlCL"));
  NcClusterRedshift *clusterz     = NC_CLUSTER_REDSHIFT (ncm_serialize_global_from_string ("NcClusterRedshiftNodist{'z-min':<0.1>, 'z-max':<1.0>}"));
  NcClusterAbundance *cad         = nc_cluster_abundance_new (mfp, NULL);
  NcDataClusterPseudoCounts *dcpc = nc_data_cluster_pseudo_counts_new (cad);
  NcmMSet *mset                   = ncm_mset_new (cosmo, clusterz, clusterm, cpc, NULL);
  NcmDataset *dset                = ncm_dataset_new ();
  NcmMatrix *m                    = ncm_matrix_new (1, 5);
  gdouble z                       = g_test_rand_double_range (0.188, 0.890);
  NcmLikelihood *lh;
  NcmFit *fit;
  gdouble m1, m2;

  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));

  m1                   = g_test_rand_double_range (1.235, 2.496); /* ln(M/M0), M0 = 10^14 h^-1 M_sun */
  m2                   = m1 + 0.4;
  test->Mobs[0]        = exp (m1) * 1.0e14;
  test->Mobs[1]        = exp (m2) * 1.0e14;
  test->Mobs_params[0] = 0.09 * test->Mobs[0];
  test->Mobs_params[1] = 0.12 * test->Mobs[1];

  g_assert_true (cpc != NULL);
  g_assert_true (clusterm != NULL);
  g_assert_true (clusterz != NULL);

  nc_halo_mass_function_set_area (mfp, 1.0);

  test->cosmo    = cosmo;
  test->reion    = reion;
  test->prim     = prim;
  test->ps_ml    = ps_ml;
  test->psf      = psf;
  test->cpc      = cpc;
  test->mfp      = mfp;
  test->clusterm = clusterm;
  test->clusterz = clusterz;
  test->z        = z;
  g_assert_true (NC_IS_CLUSTER_PSEUDO_COUNTS (cpc));

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       67.8);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.258);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  0.692);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.0482);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);

  ncm_model_orig_param_set (NCM_MODEL (test->prim), NC_HIPRIM_POWER_LAW_N_SA,       0.9608);
  ncm_model_orig_param_set (NCM_MODEL (test->prim), NC_HIPRIM_POWER_LAW_LN10E10ASA, 3.1);

  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Asz", 0.7, NULL);      /*1.0); */
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Bsz", 0.35, NULL);     /*0.2); */
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "sigma_sz", 0.1, NULL); /*vary until 0.005; */
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Al", 1.0, NULL);       /* vary until 0.001; */
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Bl", 0.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "sigma_l", 0.1, NULL); /* vary until 0.15); */
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "cor", 0.4, NULL);     /*0.4); */

  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 0, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 1, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 3, NCM_PARAM_TYPE_FREE);

  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "lnMCut",     33.5, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "sigma_Mcut",  0.05, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "zmin",        0.188, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "Deltaz",      0.72, NULL);

  ncm_matrix_set (m, 0, 0, test->z);
  ncm_matrix_set (m, 0, 1, test->Mobs[0]);
  ncm_matrix_set (m, 0, 2, test->Mobs[1]);
  ncm_matrix_set (m, 0, 3, test->Mobs_params[0]);
  ncm_matrix_set (m, 0, 4, test->Mobs_params[1]);

  nc_data_cluster_pseudo_counts_set_nclusters (dcpc, 1);
  nc_data_cluster_pseudo_counts_set_obs (dcpc, m);
  test->dcpc = dcpc;
  ncm_dataset_append_data (dset, NCM_DATA (test->dcpc));

  lh        = ncm_likelihood_new (dset);
  fit       = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  test->fit = fit;

  nc_cluster_abundance_prepare_if_needed (cad, cosmo, clusterz, clusterm);

  ncm_matrix_free (m);
  nc_distance_free (dist);
  nc_window_free (wf);
  nc_transfer_func_free (tf);
  nc_multiplicity_func_free (mulf);
  nc_cluster_abundance_free (cad);
  ncm_mset_free (mset);
  ncm_dataset_free (dset);
  ncm_likelihood_free (lh);
}

void
test_nc_cluster_pseudo_counts_1p2_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcHICosmo *cosmo           = test->cosmo;
  NcClusterMass *clusterm    = test->clusterm;
  NcClusterPseudoCounts *cpc = test->cpc;
  NcHaloMassFunction *mfp    = test->mfp;
  gdouble I1p2, I3d;

  I1p2 = nc_cluster_pseudo_counts_posterior_numerator (cpc, mfp, clusterm, cosmo, test->z, test->Mobs, test->Mobs_params);
  I3d  = nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, mfp, clusterm, cosmo, test->z, test->Mobs[0], test->Mobs[1], test->Mobs_params[0], test->Mobs_params[1]);

  ncm_assert_cmpdouble_e (I1p2, ==, I3d, 1.0e-2, 0.0);
}

void
test_nc_cluster_pseudo_counts_3d_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcHICosmo *cosmo           = test->cosmo;
  NcClusterMass *clusterm    = test->clusterm;
  NcClusterPseudoCounts *cpc = test->cpc;
  NcHaloMassFunction *mfp    = test->mfp;

  gdouble I1p2, I3d;

  I1p2 = nc_cluster_pseudo_counts_posterior_numerator (cpc, mfp, clusterm, cosmo, test->z, test->Mobs, test->Mobs_params);
  I3d  = nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, mfp, clusterm, cosmo, test->z, test->Mobs[0], test->Mobs[1], test->Mobs_params[0], test->Mobs_params[1]);

  ncm_assert_cmpdouble_e (I1p2, ==, I3d, 1.0e-2, 0.0);
}

void
test_nc_cluster_pseudo_counts_m2lnL (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMSet *mset = ncm_fit_peek_mset (test->fit);
  gdouble m2lnL;

  nc_data_cluster_pseudo_counts_init_from_sampling (test->dcpc, mset, rng, 100);

  ncm_fit_set_params_reltol (test->fit, 1.0e-5);
  ncm_fit_m2lnL_val (test->fit, &m2lnL);

  g_assert_true (gsl_finite (m2lnL));

/*
 *  ncm_fit_run (test->fit, NCM_FIT_RUN_MSGS_FULL);
 *  printf ("m2lnL = % 20.15g\n", m2lnL);
 */
  ncm_rng_clear (&rng);
}


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
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  NcMatterVar *vp;
  NcDataClusterPseudoCounts *dcpc;
  NcmFit *fit;
  gdouble z;
  gdouble Mobs[2];
  gdouble Mobs_params[2];
} TestNcClusterPseudoCounts;

void test_nc_cluster_pseudo_counts_new (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_1p2_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_3d_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_3d_integral_new_variables (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_m2lnL (TestNcClusterPseudoCounts *test, gconstpointer pdata);
void test_nc_cluster_pseudo_counts_free (TestNcClusterPseudoCounts *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  
  g_test_add ("/numcosmo/nc_cluster_pseudo_counts/1p2_integral", TestNcClusterPseudoCounts, NULL, 
              &test_nc_cluster_pseudo_counts_new, 
              &test_nc_cluster_pseudo_counts_1p2_integral, 
              &test_nc_cluster_pseudo_counts_free);
  /*g_test_add ("/numcosmo/nc_cluster_pseudo_counts/3d_integral", TestNcClusterPseudoCounts, NULL, 
              &test_nc_cluster_pseudo_counts_new, 
              &test_nc_cluster_pseudo_counts_3d_integral, 
              &test_nc_cluster_pseudo_counts_free);*/
  g_test_add ("/numcosmo/nc_cluster_pseudo_counts/3d_integral_new_variables", TestNcClusterPseudoCounts, NULL, 
              &test_nc_cluster_pseudo_counts_new, 
              &test_nc_cluster_pseudo_counts_3d_integral_new_variables, 
              &test_nc_cluster_pseudo_counts_free);
  g_test_add ("/numcosmo/nc_cluster_pseudo_counts/m2lnL", TestNcClusterPseudoCounts, NULL, 
              &test_nc_cluster_pseudo_counts_new, 
              &test_nc_cluster_pseudo_counts_m2lnL, 
              &test_nc_cluster_pseudo_counts_free);
 
  g_test_run ();
}

void _set_destroyed (gpointer b) { gboolean *destroyed = b; *destroyed = TRUE; }

void
test_nc_cluster_pseudo_counts_free (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{  
  NcClusterPseudoCounts *cpc = NC_CLUSTER_PSEUDO_COUNTS (test->cpc);
  gboolean destroyed = FALSE;
  
  nc_matter_var_free (test->vp);
  nc_hicosmo_free (test->cosmo);
  nc_cluster_mass_free (test->clusterm);
  nc_data_cluster_pseudo_counts_free (test->dcpc);
  ncm_fit_free (test->fit);
  
  g_object_set_data_full (G_OBJECT (cpc), "test-destroy", &destroyed, _set_destroyed);
  nc_cluster_pseudo_counts_free (cpc);
  g_assert (destroyed);
}

void
test_nc_cluster_pseudo_counts_new (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcHICosmo *cosmo                = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcDistance *dist                = nc_distance_new (3.0);
  NcWindow *wf                    = nc_window_new_from_name ("NcWindowTophat"); 
  NcTransferFunc *tf              = nc_transfer_func_new_from_name ("NcTransferFuncEH");
  NcMatterVar *vp                 = nc_matter_var_new (NC_MATTER_VAR_FFT, wf, tf);
  NcGrowthFunc *gf                = nc_growth_func_new ();
  NcMultiplicityFunc *mulf        = nc_multiplicity_func_new_from_name ("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}");
  NcMassFunction *mfp             = nc_mass_function_new (dist, vp, gf, mulf);
  NcClusterPseudoCounts *cpc      = NC_CLUSTER_PSEUDO_COUNTS (nc_cluster_pseudo_counts_new (mfp));
  NcClusterMass *clusterm         = NC_CLUSTER_MASS (nc_cluster_mass_new_from_name ("NcClusterMassPlCL"));
  NcDataClusterPseudoCounts *dcpc = nc_data_cluster_pseudo_counts_new ();
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (cosmo), NCM_MODEL (clusterm), NCM_MODEL (cpc));
  NcmDataset *dset                = ncm_dataset_new ();
  NcmLikelihood *lh;              
  NcmFit *fit;                    
  gdouble m1, m2;
  gdouble z                       = g_test_rand_double_range (0.188, 0.890);
  NcmMatrix *m                     = ncm_matrix_new (1, 5);
  
  m1 = g_test_rand_double_range (1.235, 2.496); /* ln(M/M0), M0 = 10^14 h^-1 M_sun */
  m2 = m1 + 0.4;
  test->Mobs[0]        = exp (m1) * 1.0e14;
  test->Mobs[1]        = exp (m2) * 1.0e14;
  test->Mobs_params[0] = 0.09 * test->Mobs[0];
  test->Mobs_params[1] = 0.12 * test->Mobs[1];
  
  g_assert (cpc != NULL);
  g_assert (clusterm != NULL);
  test->cosmo     = cosmo;
  test->vp        = vp;
  test->cpc       = NC_CLUSTER_PSEUDO_COUNTS (cpc);
  test->clusterm  = clusterm;
  test->z         = z;
  g_assert (NC_IS_CLUSTER_PSEUDO_COUNTS (cpc));

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0, 67.8);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C, 0.258);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0, 0.692);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B, 0.0482);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_SPECINDEX, 0.9608);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_SIGMA8, 0.826);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W, -1.0);

  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Asz", 0.7);  //1.0);
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Bsz", 0.35); //0.2);
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "sigma_sz", 0.005); //vary until 0.005;
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Al", 0.002); // vary until 0.001;
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "Bl", 0.0);
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "sigma_l", 0.15); // vary until 0.15);
  ncm_model_param_set_by_name (NCM_MODEL (clusterm), "cor", 0.4); //0.4);

  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 0, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 1, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 2, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 3, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 4, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 5, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (clusterm), 6, NCM_PARAM_TYPE_FREE);
  
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "lnMcut", 34.5);
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "sigma_Mcut", 0.05);
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "zmin", 0.188);
  ncm_model_param_set_by_name (NCM_MODEL (test->cpc), "Deltaz", 0.72);

  ncm_model_param_set_ftype (NCM_MODEL (test->cpc), 0, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (test->cpc), 1, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (test->cpc), 2, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (test->cpc), 3, NCM_PARAM_TYPE_FREE);

  ncm_matrix_set (m, 0, 0, test->z);
  ncm_matrix_set (m, 0, 1, test->Mobs[0]);
  ncm_matrix_set (m, 0, 2, test->Mobs[1]);
  ncm_matrix_set (m, 0, 3, test->Mobs_params[0]);
  ncm_matrix_set (m, 0, 4, test->Mobs_params[1]);
  nc_data_cluster_pseudo_counts_set_obs (dcpc, m);
  test->dcpc = dcpc;
  
  ncm_dataset_append_data (dset, NCM_DATA (test->dcpc));
  //printf ("rows = %d cols = %d np = %d dset = %d\n", ncm_matrix_nrows (m), ncm_matrix_ncols (m), test->dcpc->np, ncm_dataset_get_ndata (dset));
  
  lh        = ncm_likelihood_new (dset);
  fit       = ncm_fit_new (NCM_FIT_TYPE_NLOPT, "ln-neldermead", lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  test->fit = fit;
  
  nc_distance_free (dist);
  nc_window_free (wf);
  nc_transfer_func_free (tf);
  nc_growth_func_free (gf);
  nc_multiplicity_func_free (mulf);
  nc_mass_function_free (mfp);
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
 
  printf ("z = %.5g Msz = %.5g Ml = %.5g\n", test->z, test->Mobs[0], test->Mobs[1]);
  //ncm_model_params_log_all (NCM_MODEL (cosmo));
  //ncm_model_params_log_all (NCM_MODEL (clusterm));
  //ncm_model_params_log_all (NCM_MODEL (cpc));

  nc_matter_var_prepare (test->vp, test->cosmo);

  printf ("Integral 1p2\n");
  nc_cluster_pseudo_counts_posterior_numerator (cpc, clusterm, cosmo, test->z, test->Mobs, test->Mobs_params);  
}

void 
test_nc_cluster_pseudo_counts_3d_integral (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcHICosmo *cosmo           = test->cosmo;
  NcClusterMass *clusterm    = test->clusterm;
  NcClusterPseudoCounts *cpc = test->cpc;
 
  printf ("z = %.5g Msz = %.5g Ml = %.5g SDsz = %.5g SDl = %.5g\n", test->z, test->Mobs[0], test->Mobs[1], test->Mobs_params[0], test->Mobs_params[1]);
  //ncm_model_params_log_all (NCM_MODEL (cosmo));
  //ncm_model_params_log_all (NCM_MODEL (clusterm));
  //ncm_model_params_log_all (NCM_MODEL (cpc));

  nc_matter_var_prepare (test->vp, test->cosmo);

  printf ("Integral 3d\n");
  nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, clusterm, cosmo, test->z, test->Mobs, test->Mobs_params);
}

void 
test_nc_cluster_pseudo_counts_3d_integral_new_variables (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  NcHICosmo *cosmo           = test->cosmo;
  NcClusterMass *clusterm    = test->clusterm;
  NcClusterPseudoCounts *cpc = test->cpc;
 
  printf ("z = %.5g Msz = %.5g Ml = %.5g SDsz = %.5g SDl = %.5g\n", test->z, test->Mobs[0], test->Mobs[1], test->Mobs_params[0], test->Mobs_params[1]);
  //ncm_model_params_log_all (NCM_MODEL (cosmo));
  //ncm_model_params_log_all (NCM_MODEL (clusterm));
  //ncm_model_params_log_all (NCM_MODEL (cpc));

  nc_matter_var_prepare (test->vp, test->cosmo);

  printf ("Integral 3dnew variables\n");
  nc_cluster_pseudo_counts_posterior_numerator_plcl_new_variables (cpc, clusterm, cosmo, test->z, test->Mobs[0], test->Mobs[1], test->Mobs_params[0], test->Mobs_params[1]);
}

void 
test_nc_cluster_pseudo_counts_m2lnL (TestNcClusterPseudoCounts *test, gconstpointer pdata)
{
  gdouble m2lnL;
 
  printf ("z = %.5g Msz = %.5g Ml = %.5g SDsz = %.5g SDl = %.5g\n", test->z, test->Mobs[0], test->Mobs[1], test->Mobs_params[0], test->Mobs_params[1]);
  //ncm_model_params_log_all (NCM_MODEL (cosmo));
  //ncm_model_params_log_all (NCM_MODEL (clusterm));
  //ncm_model_params_log_all (NCM_MODEL (cpc));

  nc_matter_var_prepare (test->vp, test->cosmo);

  printf ("Test m2lnL\n");
  ncm_fit_set_params_reltol (test->fit, 1.0e-5);
  ncm_fit_m2lnL_val (test->fit, &m2lnL);
//  printf ("m2lnL = %.5g\n", m2lnL);
  //ncm_fit_log_info (test->fit);
  //ncm_fit_run (test->fit, NCM_FIT_RUN_MSGS_FULL); 
}


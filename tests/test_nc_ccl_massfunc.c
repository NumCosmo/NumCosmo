/***************************************************************************
 *            test_nc_ccl_massfunc.c
 *
 *  Fri March 08 10:09:40 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti and Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * test_nc_ccl_massfunc.c
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti and Mariana Penna Lima
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <ccl.h>

typedef struct _TestNcCCLMassFunc
{
  ccl_cosmology * ccl_cosmo;
  NcHICosmo *cosmo;
  NcGrowthFunc *gf;
  NcPowspecML *Pk;
  NcmPowspecFilter *psf;
	NcHaloMassFunction *hmf; 
} TestNcCCLMassFunc;

void test_nc_ccl_massfunc_new_model1_bbks (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model2_bbks (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model3_bbks (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model1_eh (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model2_eh (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model3_eh (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model1_cbe (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model2_cbe (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_new_model3_cbe (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_free (TestNcCCLMassFunc *test, gconstpointer pdata);

void test_nc_ccl_massfunc_cmp_m2r (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_cmp_sigma_R (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_cmp_sigma_M (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_cmp_hmf (TestNcCCLMassFunc *test, gconstpointer pdata);

void test_nc_ccl_massfunc_bbks_traps (TestNcCCLMassFunc *test, gconstpointer pdata);
void test_nc_ccl_massfunc_bbks_invalid_st (TestNcCCLMassFunc *test, gconstpointer pdata);

#define BBKS_TOLERANCE 1.0E-5
#define CTEST_DATA(name) typedef struct _##name##_data name##_data; struct _##name##_data
#define CTEST_SETUP(name) void name##_setup (name##_data* data)

typedef struct _ccl_params_data 
{
  double Omega_c;
  double Omega_b;
  double h;
  double A_s;
  double n_s;
  double sigma8;
  double Neff;
  double mnu;
  ccl_mnu_convention mnu_type;
  double Omega_v[5];
  double Omega_k[5];
  double w_0[5];
  double w_a[5];
} ccl_params_data;

void
ccl_params_data_init (ccl_params_data *data)
{
  data->Omega_c  = 0.25;
  data->Omega_b  = 0.05;
  data->h        = 0.7;
  data->A_s      = 2.1e-9;
  data->n_s      = 0.96;
  data->sigma8   = 0.8;
  data->Neff     = 3.046;
  data->mnu      = 0.0;
  data->mnu_type = ccl_mnu_sum;

  double Omega_v[5] = { 0.7,  0.7,  0.7,  0.65, 0.75};
  double w_0[5]     = {-1.0, -0.9, -0.9, -0.9, -0.9};
  double w_a[5]     = { 0.0,  0.0,  0.1,  0.1,  0.1};

  for (int i=0;i<5;i++) 
  {
    data->Omega_v[i] = Omega_v[i];
    data->w_0[i]     = w_0[i];
    data->w_a[i]     = w_a[i];
    data->Omega_k[i] = 1.0 - data->Omega_c - data->Omega_b - data->Omega_v[i];
  }
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

	g_test_add ("/nc/ccl/MassFunc/model1/bbks/cmp/m2r", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_bbks,
              &test_nc_ccl_massfunc_cmp_m2r,
              &test_nc_ccl_massfunc_free);
/*
  g_test_add ("/nc/ccl/MassFunc/model1/bbks/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_bbks,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/bbks/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_bbks,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/bbks/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_bbks,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model1/eh/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_eh,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/eh/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_eh,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/eh/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_eh,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model1/cbe/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_cbe,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/cbe/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_cbe,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/cbe/cmp/sigmaR", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_cbe,
              &test_nc_ccl_massfunc_cmp_sigma_R,
              &test_nc_ccl_massfunc_free);
*/
	g_test_add ("/nc/ccl/MassFunc/model1/bbks/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_bbks,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/bbks/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_bbks,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/bbks/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_bbks,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);
	g_test_add ("/nc/ccl/MassFunc/model1/eh/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_eh,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/eh/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_eh,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/eh/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_eh,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);
	g_test_add ("/nc/ccl/MassFunc/model1/cbe/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_cbe,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/cbe/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_cbe,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/cbe/cmp/sigmaM", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_cbe,
              &test_nc_ccl_massfunc_cmp_sigma_M,
              &test_nc_ccl_massfunc_free);
/*	
	g_test_add ("/nc/ccl/MassFunc/model1/eh/cmp/hmf", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_eh,
              &test_nc_ccl_massfunc_cmp_hmf,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model2/eh/cmp/hmf", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model2_eh,
              &test_nc_ccl_massfunc_cmp_hmf,
              &test_nc_ccl_massfunc_free);

  g_test_add ("/nc/ccl/MassFunc/model3/eh/cmp/hmf", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model3_eh,
              &test_nc_ccl_massfunc_cmp_hmf,
              &test_nc_ccl_massfunc_free);
	*/

  g_test_add ("/nc/ccl/MassFunc/model1/bbks/traps", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_bbks,
              &test_nc_ccl_massfunc_bbks_traps,
              &test_nc_ccl_massfunc_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/nc/ccl/MassFunc/model1/bbks/invalid/st/subprocess", TestNcCCLMassFunc, NULL,
              &test_nc_ccl_massfunc_new_model1_bbks,
              &test_nc_ccl_massfunc_bbks_invalid_st,
              &test_nc_ccl_massfunc_free);
#endif
  g_test_run ();
}

ccl_cosmology *
test_nc_create_ccl_cosmo (gint i_model, transfer_function_t pk_t)
{
  ccl_params_data *data    = g_new0 (ccl_params_data, 1);
  ccl_configuration config = default_config;
  gint status              = 0;
  ccl_parameters params;
  ccl_cosmology *cosmo;

  ccl_params_data_init (data);
  config.transfer_function_method = pk_t;
	config.mass_function_method     = 1;

  params = ccl_parameters_create (data->Omega_c, data->Omega_b, data->Omega_k[i_model - 1], data->Neff, &data->mnu, data->mnu_type, data->w_0[i_model - 1], data->w_a[i_model - 1], data->h, data->sigma8, data->n_s, -1, -1, -1, -1, NULL, NULL, &status);
  /*params.Omega_g = 0;*/
  /*printf ("# Setting % 22.15g % 22.15g\n", params.Omega_l, data->Omega_v[i_model - 1]);*/
  /*params.Omega_l = data->Omega_v[i_model - 1];*/
  /*params.sigma8  = data->sigma8;*/

  cosmo = ccl_cosmology_create (params, config);

  cosmo->spline_params.A_SPLINE_NLOG_PK = 200;
  cosmo->spline_params.A_SPLINE_NA_PK   = 200;
  /*cosmo->spline_params.N_K              = 600;*/

  g_free (data);

  return cosmo;
}

void 
test_nc_ccl_massfunc_create_BBKS (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  NcTransferFunc *tf       = nc_transfer_func_bbks_new ();
	NcDistance *dist         = nc_distance_new (6.0);
	NcMultiplicityFunc *mulf = nc_multiplicity_func_new_from_name ("NcMultiplicityFuncTinkerMean{'Delta':<200.0>}");

  test->Pk = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  test->gf = nc_powspec_ml_transfer_peek_gf (NC_POWSPEC_ML_TRANSFER (test->Pk));
	
  nc_transfer_func_bbks_set_type (NC_TRANSFER_FUNC_BBKS (tf), NC_TRANSFER_FUNC_BBKS_TYPE_BARYONS);

  ncm_powspec_require_zi (NCM_POWSPEC (test->Pk), 0.0);
  ncm_powspec_require_zf (NCM_POWSPEC (test->Pk), 6.0);
  ncm_powspec_require_kmin (NCM_POWSPEC (test->Pk), 1.0e-5);
  ncm_powspec_require_kmax (NCM_POWSPEC (test->Pk), 1.0e+3);

  ncm_powspec_prepare (NCM_POWSPEC (test->Pk), NCM_MODEL (test->cosmo));

  {
    NcmModel *mprim = ncm_model_peek_submodel_by_mid (NCM_MODEL (test->cosmo), nc_hiprim_id ());
    ncm_model_param_set_by_name (mprim, "ln10e10ASA", 
                                 ncm_model_param_get_by_name (mprim, "ln10e10ASA") 
                                 + 2.0 * log (test->ccl_cosmo->params.sigma8 / nc_powspec_ml_sigma_R (test->Pk, NCM_MODEL (test->cosmo), 1.0e-7, 0.0, 8.0 / nc_hicosmo_h (test->cosmo))));
  }

  ncm_powspec_prepare (NCM_POWSPEC (test->Pk), NCM_MODEL (test->cosmo));

  test->psf = ncm_powspec_filter_new (NCM_POWSPEC (test->Pk), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  ncm_powspec_filter_set_best_lnr0 (test->psf); 

	test->hmf = nc_halo_mass_function_new (dist, test->psf, mulf);

  ncm_powspec_filter_prepare (test->psf, NCM_MODEL (test->cosmo));
	nc_halo_mass_function_prepare (test->hmf, test->cosmo);

  nc_transfer_func_free (tf);
	nc_distance_free (dist);
	nc_multiplicity_func_free (mulf);
}

void 
test_nc_ccl_massfunc_create_EH (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  NcTransferFunc *tf = nc_transfer_func_eh_new ();
  NcDistance *dist         = nc_distance_new (6.0);
	NcMultiplicityFunc *mulf = nc_multiplicity_func_new_from_name ("NcMultiplicityFuncTinkerMean{'Delta':<200.0>}");
  //NcMultiplicityFunc *mulf = nc_multiplicity_func_new_from_name ("NcMultiplicityFuncTinkerMeanNormalized{'Delta':<200.0>}");
	
  nc_transfer_func_eh_set_CCL_comp (NC_TRANSFER_FUNC_EH (tf), TRUE);
  
  test->Pk = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  test->gf = nc_powspec_ml_transfer_peek_gf (NC_POWSPEC_ML_TRANSFER (test->Pk));

  ncm_powspec_require_zi (NCM_POWSPEC (test->Pk), 0.0);
  ncm_powspec_require_zf (NCM_POWSPEC (test->Pk), 6.0);
  ncm_powspec_require_kmin (NCM_POWSPEC (test->Pk), 1.0e-5);
  ncm_powspec_require_kmax (NCM_POWSPEC (test->Pk), 1.0e+3);

  ncm_powspec_prepare (NCM_POWSPEC (test->Pk), NCM_MODEL (test->cosmo));
  
  {
    NcmModel *mprim = ncm_model_peek_submodel_by_mid (NCM_MODEL (test->cosmo), nc_hiprim_id ());
    ncm_model_param_set_by_name (mprim, "ln10e10ASA", ncm_model_param_get_by_name (mprim, "ln10e10ASA") 
                                 + 2.0 * log (test->ccl_cosmo->params.sigma8 / nc_powspec_ml_sigma_R (test->Pk, NCM_MODEL (test->cosmo), 1.0e-7, 0.0, 8.0 / nc_hicosmo_h (test->cosmo))));
  }
  
  ncm_powspec_prepare (NCM_POWSPEC (test->Pk), NCM_MODEL (test->cosmo));

  test->psf = ncm_powspec_filter_new (NCM_POWSPEC (test->Pk), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  ncm_powspec_filter_set_best_lnr0 (test->psf);

	test->hmf = nc_halo_mass_function_new (dist, test->psf, mulf);

  ncm_powspec_filter_prepare (test->psf, NCM_MODEL (test->cosmo));
	nc_halo_mass_function_prepare (test->hmf, test->cosmo);

  nc_transfer_func_free (tf);
	nc_distance_free (dist);
	nc_multiplicity_func_free (mulf);
}

void 
test_nc_ccl_massfunc_create_cbe (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  NcCBE *cbe;
	NcDistance *dist         = nc_distance_new (6.0);
	NcMultiplicityFunc *mulf = nc_multiplicity_func_new_from_name ("NcMultiplicityFuncTinkerMean{'Delta':<200.0>}");

  test->Pk = NC_POWSPEC_ML (nc_powspec_ml_cbe_new ());
  cbe      = nc_powspec_ml_cbe_peek_cbe (NC_POWSPEC_ML_CBE (test->Pk));

  nc_cbe_use_ppf (cbe, TRUE);
  /*g_object_set (cbe, "verbosity", 1, NULL);*/

  ncm_powspec_require_zi (NCM_POWSPEC (test->Pk), 0.0);
  ncm_powspec_require_zf (NCM_POWSPEC (test->Pk), 6.0);
  ncm_powspec_require_kmin (NCM_POWSPEC (test->Pk), 1.0e-6);
  ncm_powspec_require_kmax (NCM_POWSPEC (test->Pk), 1.0e+3);

  ncm_powspec_prepare (NCM_POWSPEC (test->Pk), NCM_MODEL (test->cosmo));

  test->psf = ncm_powspec_filter_new (NCM_POWSPEC (test->Pk), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  ncm_powspec_filter_set_best_lnr0 (test->psf);

	test->hmf = nc_halo_mass_function_new (dist, test->psf, mulf);

  ncm_powspec_filter_prepare (test->psf, NCM_MODEL (test->cosmo));
	nc_halo_mass_function_prepare (test->hmf, test->cosmo);

	nc_distance_free (dist);
	nc_multiplicity_func_free (mulf);
}

void
test_nc_ccl_massfunc_new_model1_bbks (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, ccl_bbks);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_BBKS (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model2_bbks (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, ccl_bbks);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_BBKS (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model3_bbks (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, ccl_bbks);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_BBKS (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model1_eh (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, ccl_eisenstein_hu);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_EH (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model2_eh (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, ccl_eisenstein_hu);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_EH (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model3_eh (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, ccl_eisenstein_hu);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_EH (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model1_cbe (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (1, ccl_boltzmann_class);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_cbe (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model2_cbe (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (2, ccl_boltzmann_class);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_cbe (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_new_model3_cbe (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  test->ccl_cosmo = test_nc_create_ccl_cosmo (3, ccl_boltzmann_class);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  /*ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Tgamma0", 0.0);*/

  test_nc_ccl_massfunc_create_cbe (test, pdata);
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_massfunc_free (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
	NCM_TEST_FREE (nc_halo_mass_function_free, test->hmf);
  NCM_TEST_FREE (ncm_powspec_filter_free, test->psf);
  NCM_TEST_FREE (nc_powspec_ml_free, test->Pk);

  ccl_parameters_free (&test->ccl_cosmo->params);
  ccl_cosmology_free (test->ccl_cosmo);  
}

void
test_nc_ccl_massfunc_cmp_m2r (TestNcCCLMassFunc *test, gconstpointer pdata)
{
	const gint ntests   = 100;
	const gdouble tol   = 5.0e-10;
	  
  gint status = 0;
  gint i;

  for (i = 0; i < ntests; i++)
  {
	  const gdouble log10M  = 9.0 + (8.0 * i) / (ntests - 1.0) ;
    const gdouble M       = pow (10.0, log10M);
		const gdouble lnR     = nc_halo_mass_function_lnM_to_lnR (test->hmf, test->cosmo, log (M));
		const gdouble ncR     = exp (lnR);
		const gdouble cclR    = ccl_massfunc_m2r(test->ccl_cosmo, M, &status); 
			
		ncm_assert_cmpdouble_e (ncR, ==, cclR, tol, 0.0);
		//printf ("% 22.15e | % 22.15g % 22.15g %17.10e\n", M, ncR, cclR, fabs (cclR / ncR - 1.0));
  }
}

void
test_nc_ccl_massfunc_cmp_sigma_R (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  const gint ntests   = 200;
  const gdouble z_max = 6.0;
  const gdouble tol   = 5.0e-5;
  const gdouble Rmin  = 0.01;
  const gdouble Rmax  = 100.0;

  gint status = 0;
  gint i;

  for (i = 0; i < ntests; i++)
  {
    const gdouble z     = z_max / (1.0 * ntests) * i;
    const gdouble a     = 1.0 / (1.0 + z);
    gint j;
    
    for (j = 0; j < ntests; j++)
    {
      const gdouble R          = exp (log (Rmin) + log (Rmax / Rmin) * j / (ntests - 1.0));
      /* const gdouble sigmaR_int = nc_powspec_ml_sigma_R (test->Pk, NCM_MODEL (test->cosmo), 1.0e-5, z, R); */
      const gdouble sigmaR     = ncm_powspec_filter_eval_sigma (test->psf, z, R);
      const gdouble cclsigmaR  = ccl_sigmaR (test->ccl_cosmo, R, a, &status);

      ncm_assert_cmpdouble_e (sigmaR, ==, cclsigmaR, tol, 0.0);
      /* printf ("% 22.15g % 22.15g | % 22.15g % 22.15g % 22.15g % 17.10e %17.10e\n", 
              R, z, sigmaR_fft, sigmaR_int, cclsigmaR, fabs (sigmaR / sigmaR_int - 1.0), fabs (cclsigmaR / sigmaR_int - 1.0)); */
    }
  }
}

void
test_nc_ccl_massfunc_cmp_sigma_M (TestNcCCLMassFunc *test, gconstpointer pdata)
{
	const gdouble z_max = 6.0;
  const gdouble tol   = 5.0e-5;
	const gint skip     = 3;
  gint status = 0;
  gint i;

	ccl_cosmology_compute_sigma (test->ccl_cosmo, &status);
	{
	  const gint ntests   = test->ccl_cosmo->data.logsigma->size;

	  for (i = 0; i < ntests; i += skip)
  	{
		  const gdouble z     = z_max / (1.0 * ntests) * i;
		  const gdouble a     = 1.0 / (1.0 + z);
		  gint j;

		  for (j = 0; j < ntests; j += skip)
		  {
			  const gdouble log10M     = test->ccl_cosmo->data.logsigma->x[j];
			  const gdouble M          = pow (10.0, log10M);
			  const gdouble R          = ccl_massfunc_m2r(test->ccl_cosmo, M, &status);
		  	const gdouble sigmaR_fft = ncm_powspec_filter_eval_sigma (test->psf, z, R);
			  /* const gdouble sigmaR_int = nc_powspec_ml_sigma_R (test->Pk, NCM_MODEL (test->cosmo), 1.0e-5, z, R); */
		  	const gdouble cclsigmaM  = ccl_sigmaM (test->ccl_cosmo, M, a, &status);

	  		ncm_assert_cmpdouble_e (sigmaR_fft, ==, cclsigmaM, tol, 0.0);
  			/*printf ("% 22.15g % 22.15e % 22.15g | % 22.15g % 22.15g % 22.15g %17.10e %17.10e %17.10e\n", R, M, z, 
			            sigmaR_fft, sigmaR_int, cclsigmaM, 
		  	          fabs (sigmaR_fft / sigmaR_int - 1.0), 
	  							fabs (cclsigmaM / sigmaR_int - 1.0), 
  								fabs (cclsigmaM / cclsigmaR - 1.0)); */ 
		  }
	  }
  }
}

void
test_nc_ccl_massfunc_cmp_hmf (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  const gint ntests   = 100;
  const gdouble z_max = 2.0;
  const gdouble tol   = 5.0e-4;
  const gdouble Delta = 200.0;

  gint status = 0;
  gint i;

  for (i = 0; i < ntests; i++)
  {
    const gdouble z     = z_max / (1.0 * ntests) * i;
    const gdouble a     = 1.0 / (1.0 + z);
    gint j;
    
    for (j = 0; j < ntests; j++)
    {
      const gdouble M          = exp (log (10e12) + log (10e16 / 10e12) * j / (ntests - 1.0));
			const gdouble lnM        = log (M);
			const gdouble nchmf      = log (10.0) * nc_halo_mass_function_dn_dlnM (test->hmf, test->cosmo, lnM, z);
      const gdouble cclhmf     = ccl_massfunc (test->ccl_cosmo, M, a, Delta, &status);

			/* Absolute error 1.0e-30. The test does not pass with tol for few points where the hmf is irrelevant. */
      ncm_assert_cmpdouble_e (nchmf, ==, cclhmf, tol, 1.0e-25);fflush (stderr);
      //printf ("MF: % 22.15e % 22.15g | % 22.15g % 22.15g %17.10e\n", M, z, nchmf, cclhmf, fabs (cclhmf / nchmf - 1.0));fflush (stdout);
    }
  }
}

void
test_nc_ccl_massfunc_bbks_traps (TestNcCCLMassFunc *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/nc/ccl/MassFunc/model1/bbks/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_nc_ccl_massfunc_bbks_invalid_st (TestNcCCLMassFunc *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

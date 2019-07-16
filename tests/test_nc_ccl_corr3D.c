/***************************************************************************
 *            test_nc_ccl_Corr3D.c
 *
 *  Thu March 14 09:47:56 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * test_nc_ccl_Corr3D.c
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
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

typedef struct _TestNcCCLCorr3D
{
  ccl_cosmology *ccl_cosmo;
  NcHICosmo *cosmo;
  NcmPowspec *Pk;
  NcmPowspecCorr3d *psc;
} TestNcCCLCorr3D;

void test_nc_ccl_corr_3d_new (TestNcCCLCorr3D *test, gconstpointer pdata);
void test_nc_ccl_corr_3d_free (TestNcCCLCorr3D *test, gconstpointer pdata);

void test_nc_ccl_corr_3d_create_BBKS (TestNcCCLCorr3D *test, gconstpointer pdata);
void test_nc_ccl_corr_3d_create_EH (TestNcCCLCorr3D *test, gconstpointer pdata);
void test_nc_ccl_corr_3d_create_CBE (TestNcCCLCorr3D *test, gconstpointer pdata);

void test_nc_ccl_corr_3d_cmp_corr_3d (TestNcCCLCorr3D *test, gconstpointer pdata);

void test_nc_ccl_corr_3d_bbks_traps (TestNcCCLCorr3D *test, gconstpointer pdata);
void test_nc_ccl_corr_3d_bbks_invalid_st (TestNcCCLCorr3D *test, gconstpointer pdata);

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
  gdouble Omega_v[5] = { 0.7,  0.7,  0.7,  0.65, 0.75};
  gdouble w_0[5]     = {-1.0, -0.9, -0.9, -0.9, -0.9};
  gdouble w_a[5]     = { 0.0,  0.0,  0.1,  0.1,  0.1};
  gint i;

  data->Omega_c  = 0.25;
  data->Omega_b  = 0.05;
  data->h        = 0.7;
  data->A_s      = 2.1e-9;
  data->n_s      = 0.96;
  data->sigma8   = 0.8;
  data->Neff     = 3.046;
  data->mnu      = 0.0;
  data->mnu_type = ccl_mnu_sum;
  
  for (i = 0; i < 5; i++) 
  {
    data->Omega_v[i] = Omega_v[i];
    data->w_0[i]     = w_0[i];
    data->w_a[i]     = w_a[i];
    data->Omega_k[i] = 1.0 - data->Omega_c - data->Omega_b - data->Omega_v[i];
  }
}

typedef struct _TestNcCCLCorr3DPk
{
  const gchar *path;
  void (*create_pk) (TestNcCCLCorr3D *test, gconstpointer pdata);
  transfer_function_t pk_t;
} TestNcCCLCorr3DPk;

typedef struct _TestNcCCLCorr3DModel
{
  const gchar *path;
  gint model_i;
} TestNcCCLCorr3DModel;

typedef struct _TestNcCCLCorr3DCMP
{
  const gchar *path;
  void (*cmp) (TestNcCCLCorr3D *test, gconstpointer pdata);
} TestNcCCLCorr3DCMP;

typedef struct _TestNcCCLCorr3DData
{
  TestNcCCLCorr3DPk pk;
  TestNcCCLCorr3DModel model;
  TestNcCCLCorr3DCMP cmp;
} TestNcCCLCorr3DData;

#define TEST_NC_CCL_NLPK_NPKS 1
TestNcCCLCorr3DPk pks[1] = {
  {"cbe+Halofit", &test_nc_ccl_corr_3d_create_CBE, ccl_boltzmann_class},
};

#define TEST_NC_CCL_NLPK_NMODELS 1
TestNcCCLCorr3DModel models[1] = {
  {"model1", 1},
/*  {"model2", 2},*/
/*  {"model3", 3},*/
/*  {"model4", 4},*/
/*  {"model5", 5},*/
};

#define TEST_NC_CCL_NLPK_NCMPS 1
TestNcCCLCorr3DCMP cmps[1] = {
  {"Corr3D", &test_nc_ccl_corr_3d_cmp_corr_3d},
};

#define NTESTS ((TEST_NC_CCL_NLPK_NPKS)*(TEST_NC_CCL_NLPK_NMODELS)*(TEST_NC_CCL_NLPK_NCMPS))

gint
main (gint argc, gchar *argv[])
{
  TestNcCCLCorr3DData data[NTESTS];
  gint i, m;
 
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  m = 0;
  for (i = 0; i < TEST_NC_CCL_NLPK_NPKS; i++)
  {
    gint j;
    for (j = 0; j < TEST_NC_CCL_NLPK_NMODELS; j++)
    {
      gint k;
      for (k = 0; k < TEST_NC_CCL_NLPK_NCMPS; k++)
      {
        gchar *path   = g_strdup_printf ("/nc/ccl/Corr3D/%s/%s/%s", pks[i].path, models[j].path, cmps[k].path);
        data[m].pk    = pks[i];
        data[m].model = models[j];
        data[m].cmp   = cmps[k];

        g_test_add (path, TestNcCCLCorr3D, &data[m], &test_nc_ccl_corr_3d_new, cmps[k].cmp, &test_nc_ccl_corr_3d_free);

        g_free (path);
        m++;
      }
    }
  }

  g_test_add ("/nc/ccl/Pk/bbks/model1/traps", TestNcCCLCorr3D, &data[0],
              &test_nc_ccl_corr_3d_new,
              &test_nc_ccl_corr_3d_bbks_traps,
              &test_nc_ccl_corr_3d_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/nc/ccl/Pk/model1/bbks/invalid/st/subprocess", TestNcCCLCorr3D, &data[0],
              &test_nc_ccl_corr_3d_new,
              &test_nc_ccl_corr_3d_bbks_invalid_st,
              &test_nc_ccl_corr_3d_free);
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
  config.transfer_function_method     = pk_t;
  config.matter_power_spectrum_method = ccl_halofit;

  params = ccl_parameters_create (data->Omega_c, data->Omega_b, data->Omega_k[i_model - 1], data->Neff, &data->mnu, data->mnu_type, data->w_0[i_model - 1], data->w_a[i_model - 1], data->h, data->sigma8, data->n_s, -1, -1, -1, 0.0, 0.0, -1, NULL, NULL, &status);
  cosmo = ccl_cosmology_create (params, config);

  /*cosmo->spline_params.A_SPLINE_NLOG_PK = 200;*/
  /*cosmo->spline_params.A_SPLINE_NA_PK   = 200;*/
  /*cosmo->spline_params.N_K              = 100;*/
	/*cosmo->spline_params.K_MAX_SPLINE     = 1.0e3;*/

  g_free (data);

  return cosmo;
}

void 
test_nc_ccl_corr_3d_create_BBKS (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  NcTransferFunc *tf = nc_transfer_func_bbks_new ();
  NcPowspecML *ps_ml = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf)); 

  nc_transfer_func_bbks_set_type (NC_TRANSFER_FUNC_BBKS (tf), NC_TRANSFER_FUNC_BBKS_TYPE_BARYONS);

  ncm_powspec_require_zi (NCM_POWSPEC (ps_ml), 0.0);
  ncm_powspec_require_zf (NCM_POWSPEC (ps_ml), 6.0);
  ncm_powspec_require_kmin (NCM_POWSPEC (ps_ml), 1.0e-6);
  ncm_powspec_require_kmax (NCM_POWSPEC (ps_ml), 1.0e+3);
  ncm_powspec_prepare (NCM_POWSPEC (ps_ml), NCM_MODEL (test->cosmo));

  {
    NcmModel *mprim = ncm_model_peek_submodel_by_mid (NCM_MODEL (test->cosmo), nc_hiprim_id ());
    ncm_model_param_set_by_name (mprim, "ln10e10ASA", 
                                 ncm_model_param_get_by_name (mprim, "ln10e10ASA") 
                                 + 2.0 * log (test->ccl_cosmo->params.sigma8 / ncm_powspec_sigma_tophat_R (NCM_POWSPEC (ps_ml), NCM_MODEL (test->cosmo), 1.0e-7, 0.0, 8.0 / nc_hicosmo_h (test->cosmo))));
  }

  test->Pk = NCM_POWSPEC (nc_powspec_mnl_halofit_new (ps_ml, 6.0, 1.0e-5));

  ncm_powspec_require_zi (test->Pk, 0.0);
  ncm_powspec_require_zf (test->Pk, 6.0);
  ncm_powspec_require_kmin (test->Pk, 1.0e-6);
  ncm_powspec_require_kmax (test->Pk, 1.0e+3);

  ncm_powspec_prepare (test->Pk, NCM_MODEL (test->cosmo));

  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
}

void 
test_nc_ccl_corr_3d_create_EH (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  NcTransferFunc *tf = nc_transfer_func_eh_new ();
  NcPowspecML *ps_ml = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf)); 

  nc_transfer_func_eh_set_CCL_comp (NC_TRANSFER_FUNC_EH (tf), TRUE);
  
  ncm_powspec_require_zi (NCM_POWSPEC (ps_ml), 0.0);
  ncm_powspec_require_zf (NCM_POWSPEC (ps_ml), 6.0);
  ncm_powspec_require_kmin (NCM_POWSPEC (ps_ml), 1.0e-6);
  ncm_powspec_require_kmax (NCM_POWSPEC (ps_ml), 1.0e+3);
  ncm_powspec_prepare (NCM_POWSPEC (ps_ml), NCM_MODEL (test->cosmo));
  
  {
    NcmModel *mprim = ncm_model_peek_submodel_by_mid (NCM_MODEL (test->cosmo), nc_hiprim_id ());
    ncm_model_param_set_by_name (mprim, "ln10e10ASA", ncm_model_param_get_by_name (mprim, "ln10e10ASA") 
                                 + 2.0 * log (test->ccl_cosmo->params.sigma8 / ncm_powspec_sigma_tophat_R (NCM_POWSPEC (ps_ml), NCM_MODEL (test->cosmo), 1.0e-7, 0.0, 8.0 / nc_hicosmo_h (test->cosmo))));
  }
  
  test->Pk = NCM_POWSPEC (nc_powspec_mnl_halofit_new (ps_ml, 6.0, 1.0e-5));
  ncm_powspec_require_zi (test->Pk, 0.0);
  ncm_powspec_require_zf (test->Pk, 6.0);
  ncm_powspec_require_kmin (test->Pk, 1.0e-6);
  ncm_powspec_require_kmax (test->Pk, 1.0e+3);

  ncm_powspec_prepare (test->Pk, NCM_MODEL (test->cosmo));
  
  nc_transfer_func_free (tf);
}

void 
test_nc_ccl_corr_3d_create_CBE (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  NcPowspecML *ps_ml = NC_POWSPEC_ML (nc_powspec_ml_cbe_new ());
  NcCBE *cbe         = nc_powspec_ml_cbe_peek_cbe (NC_POWSPEC_ML_CBE (ps_ml));

	nc_powspec_ml_cbe_set_intern_k_min (NC_POWSPEC_ML_CBE (ps_ml), test->ccl_cosmo->spline_params.K_MIN);
	nc_powspec_ml_cbe_set_intern_k_max (NC_POWSPEC_ML_CBE (ps_ml), test->ccl_cosmo->spline_params.K_MAX_SPLINE);
  nc_cbe_use_ppf (cbe, TRUE);
  /*g_object_set (cbe, "verbosity", 1, NULL);*/

  test->Pk = NCM_POWSPEC (nc_powspec_mnl_halofit_new (ps_ml, 6.0, 1.0e-5));
  ncm_powspec_require_zi (test->Pk, 0.0);
  ncm_powspec_require_zf (test->Pk, 6.0);
  ncm_powspec_require_kmin (test->Pk, 1.0e-6);
  ncm_powspec_require_kmax (test->Pk, 1.0e+3);

  ncm_powspec_prepare (test->Pk, NCM_MODEL (test->cosmo));

  nc_powspec_ml_free (ps_ml);
}

void
test_nc_ccl_corr_3d_new (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  const TestNcCCLCorr3DData *data = pdata;

  test->ccl_cosmo = test_nc_create_ccl_cosmo (data->model.model_i, data->pk.pk_t);
  test->cosmo     = NC_HICOSMO (nc_hicosmo_de_cpl_new_from_ccl (&test->ccl_cosmo->params));

  data->pk.create_pk (test, pdata);

  test->psc = ncm_powspec_corr3d_new (test->Pk);
  
  ncm_powspec_corr3d_prepare (test->psc, NCM_MODEL (test->cosmo));
  
  g_assert (NC_IS_HICOSMO_DE_CPL (test->cosmo));
}

void
test_nc_ccl_corr_3d_free (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_powspec_corr3d_free, test->psc);
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
  NCM_TEST_FREE (ncm_powspec_free, test->Pk);

  ccl_parameters_free (&test->ccl_cosmo->params);
  ccl_cosmology_free (test->ccl_cosmo);  
}

#define NR 200

void
test_nc_ccl_corr_3d_cmp_corr_3d (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  const gint ntests      = 20;
  const gdouble z_max    = 2.0;
  const gdouble reltol   = 1.0e-3;
  const gdouble abstol   = 1.0e-4;
  gint status = 0;
  gdouble r[NR];
  gdouble xi[NR];
  gint i;
  
  for (i = 0; i < NR; i++)
    r[i] = exp (log (5.0e-1) + log (2.0e3) * i / (NR - 1.0));
  
  for (i = 0; i < ntests; i++)
  {
    const gdouble z = z_max / (1.0 * ntests) * i;
    const gdouble a = 1.0 / (1.0 + z);
    gint j;

    ccl_correlation_3d (test->ccl_cosmo, a, NR, r, xi, 0, NULL, &status);

    for (j = 0; j < NR; j++)
    {
      const gdouble r_j       = r[j];
      const gdouble Corr3D    = ncm_powspec_corr3d_eval_xi (test->psc, z, r_j);
      /*const gdouble Corr3DInt = ncm_powspec_corr3d (test->Pk, NCM_MODEL (test->cosmo), 1.0e-6, z, r_j);*/
      const gdouble cclCorr3D = xi[j];

      ncm_assert_cmpdouble_e (Corr3D, ==, cclCorr3D, reltol, abstol);
      /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %17.10e %17.10e %17.10e\n", r_j, z, Corr3DInt, Corr3D, cclCorr3D, Corr3D / Corr3DInt, cclCorr3D / Corr3DInt, fabs (Corr3D / Corr3DInt - 1.0), fabs (cclCorr3D / Corr3DInt - 1.0), fabs (cclCorr3D / Corr3D - 1.0));*/
    }
  }
}

void
test_nc_ccl_corr_3d_bbks_traps (TestNcCCLCorr3D *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/nc/ccl/Pk/model1/bbks/invalid/st/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_nc_ccl_corr_3d_bbks_invalid_st (TestNcCCLCorr3D *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

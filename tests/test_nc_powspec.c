/***************************************************************************
 *            test_nc_powspec.c
 *
 *  Mon Dec 12 08:30:12 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2022 <vitenti@uel.br>
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

typedef struct _TestNcPowspec TestNcPowspec;

struct _TestNcPowspec
{
  NcmPowspec *ps;
  NcmModel *model;
};

typedef void (*TestNcPowspecFunc) (TestNcPowspec *test, gconstpointer pdata);

void test_nc_powspec_ml_transfer_new_EH (TestNcPowspec *test, gconstpointer pdata);
void test_nc_powspec_ml_transfer_new_BBKS (TestNcPowspec *test, gconstpointer pdata);
void test_nc_powspec_ml_cbe_new (TestNcPowspec *test, gconstpointer pdata);
void test_nc_powspec_mnl_halofit_new (TestNcPowspec *test, gconstpointer pdata);

void test_nc_powspec_eval (TestNcPowspec *test, gconstpointer pdata);
void test_nc_powspec_filter_tophat (TestNcPowspec *test, gconstpointer pdata);
void test_nc_powspec_corr3d (TestNcPowspec *test, gconstpointer pdata);

void test_nc_powspec_free (TestNcPowspec *test, gconstpointer pdata);

typedef struct _TestNcPowspecFunc
{
  void (*func) (TestNcPowspec *, gconstpointer);
  const gchar *name;
  gpointer pdata;
} TestNcPowspecFuncData;

#define TEST_NC_POWSPECS_LEN 6
TestNcPowspecFuncData powspecs[TEST_NC_POWSPECS_LEN] =
{
  {test_nc_powspec_ml_transfer_new_EH,   "ml/transfer/EH",               NULL},
  {test_nc_powspec_ml_transfer_new_BBKS, "ml/transfer/BBKS",             NULL},
  {test_nc_powspec_ml_cbe_new,           "ml/cbe",                       NULL},
  {test_nc_powspec_mnl_halofit_new,      "mnl/halofit/ml/transfer/EH",   test_nc_powspec_ml_transfer_new_EH},
  {test_nc_powspec_mnl_halofit_new,      "mnl/halofit/ml/transfer/BBKS", test_nc_powspec_ml_transfer_new_BBKS},
  {test_nc_powspec_mnl_halofit_new,      "mnl/halofit/ml/cbe",           test_nc_powspec_ml_cbe_new},
};

#define TEST_NC_POWSPEC_TESTS 3
TestNcPowspecFuncData tests[TEST_NC_POWSPEC_TESTS] =
{
    {test_nc_powspec_eval,          "eval",          NULL},
    {test_nc_powspec_filter_tophat, "filter/tophat", NULL},
    {test_nc_powspec_corr3d,        "corr3d",        NULL},
};

gint
main (gint argc, gchar *argv[])
{
  gint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  for (i = 0; i < TEST_NC_POWSPECS_LEN; i++)
  {
    for (j = 0; j < TEST_NC_POWSPEC_TESTS; j++)
    {
      gchar *test_path = g_strdup_printf ("/nc/powspec/%s/%s", powspecs[i].name, tests[j].name);
      g_test_add (test_path, TestNcPowspec, powspecs[i].pdata, powspecs[i].func, tests[j].func,
          &test_nc_powspec_free);
      g_free (test_path);
    }
  }
  
  g_test_run ();
}

void
test_nc_powspec_ml_transfer_new_EH (TestNcPowspec *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcHIReion *reion            = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim              = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcTransferFunc *tf          = nc_transfer_func_new_from_name ("NcTransferFuncEH");
  NcPowspecML *ps_ml          = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
  
  test->model = NCM_MODEL (nc_hicosmo_ref (cosmo));
  test->ps    = NCM_POWSPEC (nc_powspec_ml_ref (ps_ml));
  
  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
  
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
}

void
test_nc_powspec_ml_transfer_new_BBKS (TestNcPowspec *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcHIReion *reion            = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim              = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcTransferFunc *tf          = nc_transfer_func_new_from_name ("NcTransferFuncBBKS");
  NcPowspecML *ps_ml          = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
  
  test->model = NCM_MODEL (nc_hicosmo_ref (cosmo));
  test->ps    = NCM_POWSPEC (nc_powspec_ml_ref (ps_ml));
  
  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
  
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
}

void
test_nc_powspec_ml_cbe_new (TestNcPowspec *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcHIReion *reion            = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim              = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcPowspecML *ps_ml          = NC_POWSPEC_ML (nc_powspec_ml_cbe_new ());

  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));

  test->model = NCM_MODEL (nc_hicosmo_ref (cosmo));
  test->ps    = NCM_POWSPEC (nc_powspec_ml_ref (ps_ml));

  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));

  nc_powspec_ml_free (ps_ml);
}

void
test_nc_powspec_mnl_halofit_new (TestNcPowspec *test, gconstpointer pdata)
{
  ((TestNcPowspecFunc)pdata)(test, NULL);
  {
    NcPowspecMNLHaloFit *ps_mln = nc_powspec_mnl_halofit_new (NC_POWSPEC_ML (test->ps), 3.0, 1.0e-4);

    ncm_powspec_free (test->ps);
    test->ps = NCM_POWSPEC (ps_mln);
  }
}

void
test_nc_powspec_eval (TestNcPowspec *test, gconstpointer pdata)
{
  const gdouble zi   = ncm_powspec_get_zi (test->ps);
  const gdouble zf   = ncm_powspec_get_zf (test->ps);
  const gdouble kmin = ncm_powspec_get_kmin (test->ps);
  const gdouble kmax = ncm_powspec_get_kmax (test->ps);
  NcmVector *kv      = NULL;
  NcmVector *Pkv     = NULL;
  guint Nz           = 0;
  guint Nk           = 0;
  gint i, j;

  g_assert_cmpfloat (zi, >=, 0.0);
  g_assert_cmpfloat (zi, <, zf);
  g_assert_cmpfloat (kmin, >, 0.0);
  g_assert_cmpfloat (kmin, <, kmax);

  ncm_powspec_prepare (test->ps, test->model);
  ncm_powspec_get_nknots (test->ps, &Nz, &Nk);

  g_assert_cmpuint (Nz, >, 0);
  g_assert_cmpuint (Nk, >, 0);

  kv  = ncm_vector_new (Nk * 10);
  Pkv = ncm_vector_new (Nk * 10);

  for (i = 0; i < Nk * 10; i++)
  {
    const gdouble lnk = log (kmin) + log (kmax / kmin) / (Nk * 10.0 - 1.0) * i;
    ncm_vector_set (kv, i, exp (lnk));
  }

  for (i = 0; i < Nz * 10; i++)
  {
    const gdouble z = zi + (zf - zi) / (Nz * 10.0 - 1.0) * i;
    ncm_powspec_eval_vec (test->ps, test->model, z, kv, Pkv);
    for (j = 0; j < Nk * 10; j++)
    {
      const gdouble k = ncm_vector_get (kv, j);

      ncm_assert_cmpdouble_e (ncm_vector_get (Pkv, j), ==, ncm_powspec_eval (test->ps, test->model, z, k), 1.0e-10, 0.0);
    }
  }

  {
    NcmSpline2d *Pks = ncm_powspec_get_spline_2d (test->ps, test->model);

    for (i = 0; i < Nz * 10; i++)
    {
      const gdouble z = zi + (zf - zi) / (Nz * 10.0 - 1.0) * i;
      ncm_powspec_eval_vec (test->ps, test->model, z, kv, Pkv);
      for (j = 0; j < Nk * 10; j++)
      {
        const gdouble k = ncm_vector_get (kv, j);

        ncm_assert_cmpdouble_e (ncm_vector_get (Pkv, j), ==, ncm_spline2d_eval (Pks, z, k),
            test->ps->reltol_spline * 10.0, 0.0);
      }
    }

    ncm_spline2d_free (Pks);
  }

  ncm_vector_free (kv);
  ncm_vector_free (Pkv);
}

void
test_nc_powspec_filter_tophat (TestNcPowspec *test, gconstpointer pdata)
{
  NcmPowspecFilter *psf = ncm_powspec_filter_new (NCM_POWSPEC (test->ps), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  const gdouble zi   = ncm_powspec_get_zi (test->ps);
  const gdouble zf   = ncm_powspec_get_zf (test->ps);
  const gdouble kmin = ncm_powspec_get_kmin (test->ps);
  const gdouble kmax = ncm_powspec_get_kmax (test->ps);

  ncm_powspec_filter_prepare (psf, test->model);

  g_assert_cmpfloat (zi, >=, 0.0);
  g_assert_cmpfloat (zi, <, zf);
  g_assert_cmpfloat (kmin, >, 0.0);
  g_assert_cmpfloat (kmin, <, kmax);

  {
    const gdouble r_min = ncm_powspec_filter_get_r_min (psf);
    const gdouble r_max = ncm_powspec_filter_get_r_max (psf);
    gint i, j;

    for (i = 0; i < 10; i++)
    {
      const gdouble z = zi + (zf - zi) / (100.0 - 1.0) * i;
      for (j = 0; j < 10; j++)
      {
        const gdouble lnR    = log (r_min) + log (r_max / r_min) / (100.0 - 1.0) * j;
        const gdouble R      = exp (lnR);
        const gdouble var0   = ncm_powspec_var_tophat_R (test->ps, test->model, psf->reltol, z, R);
        const gdouble sigma0 = ncm_powspec_sigma_tophat_R (test->ps, test->model, psf->reltol, z, R);
        const gdouble var1   = ncm_powspec_filter_eval_var_lnr (psf, z, lnR);
        const gdouble sigma1 = ncm_powspec_filter_eval_sigma_lnr (psf, z, lnR);

        if (var0 > 1.0e-4)
        {
          ncm_assert_cmpdouble_e (var0, ==, var1, psf->reltol * 10.0, 0.0);
          ncm_assert_cmpdouble_e (sigma0, ==, sigma1, psf->reltol * 10.0, 0.0);
        }
      }
    }
  }
  ncm_powspec_filter_free (psf);
}


void
test_nc_powspec_corr3d (TestNcPowspec *test, gconstpointer pdata)
{
  NcmPowspecCorr3d *psc = ncm_powspec_corr3d_new (NCM_POWSPEC (test->ps));
  const gdouble zi      = ncm_powspec_get_zi (test->ps);
  const gdouble zf      = ncm_powspec_get_zf (test->ps);
  const gdouble kmin    = ncm_powspec_get_kmin (test->ps);
  const gdouble kmax    = ncm_powspec_get_kmax (test->ps);

  ncm_powspec_corr3d_prepare (psc, test->model);

  g_assert_cmpfloat (zi, >=, 0.0);
  g_assert_cmpfloat (zi, <, zf);
  g_assert_cmpfloat (kmin, >, 0.0);
  g_assert_cmpfloat (kmin, <, kmax);

  {
    const gdouble r_min = ncm_powspec_corr3d_get_r_min (psc);
    const gdouble r_max = ncm_powspec_corr3d_get_r_max (psc);
    gint i, j;

    for (i = 0; i < 10; i++)
    {
      const gdouble z = zi + (zf - zi) / (100.0 - 1.0) * i;
      for (j = 0; j < 10; j++)
      {
        const gdouble lnR = log (r_min) + log (r_max / r_min) / (100.0 - 1.0) * j;
        const gdouble R   = exp (lnR);
        const gdouble xi0 = ncm_powspec_corr3d (test->ps, test->model, psc->reltol, z, R);
        const gdouble xi1 = ncm_powspec_corr3d_eval_xi_lnr (psc, z, lnR);

        ncm_assert_cmpdouble_e (xi0, ==, xi1, psc->reltol * 10.0, 0.0);
      }
    }
  }
  ncm_powspec_corr3d_free (psc);
}

void
test_nc_powspec_free (TestNcPowspec *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_powspec_free, test->ps);
  NCM_TEST_FREE (ncm_model_free, test->model);
}


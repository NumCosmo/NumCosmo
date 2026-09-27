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

void test_nc_powspec_ml_cbe_extrapolation (void);

typedef struct _TestNcPowspecFunc
{
  void (*func) (TestNcPowspec *, gconstpointer);

  const gchar *name;
  gpointer pdata;
} TestNcPowspecFuncData;

/* Built twice from this source (see tests/c/meson.build): with -DPOWSPEC_SPLIT_CBE only the
 * CLASS-backed (cbe) spectra, which are FAST here (the spectrum is splined once, ~4s); with
 * -DPOWSPEC_SPLIT_TRANSFER only the analytic transfer-function spectra (EH/BBKS + halofit +
 * corr3d), which are the SLOW half (>2 min, re-evaluated pointwise). Splitting lets the fast
 * cbe half stay on the fast lane (unit) and the heavy transfer half move to the coverage-only
 * acceptance tier. With neither macro all six run (local default). */
TestNcPowspecFuncData powspecs[] =
{
#ifndef POWSPEC_SPLIT_CBE
  {test_nc_powspec_ml_transfer_new_EH,   "ml/transfer/EH",               NULL},
  {test_nc_powspec_ml_transfer_new_BBKS, "ml/transfer/BBKS",             NULL},
  {test_nc_powspec_mnl_halofit_new,      "mnl/halofit/ml/transfer/EH",   test_nc_powspec_ml_transfer_new_EH},
  {test_nc_powspec_mnl_halofit_new,      "mnl/halofit/ml/transfer/BBKS", test_nc_powspec_ml_transfer_new_BBKS},
#endif
#ifndef POWSPEC_SPLIT_TRANSFER
  {test_nc_powspec_ml_cbe_new,           "ml/cbe",                       NULL},
  {test_nc_powspec_mnl_halofit_new,      "mnl/halofit/ml/cbe",           test_nc_powspec_ml_cbe_new},
#endif
};

#define TEST_NC_POWSPECS_LEN G_N_ELEMENTS (powspecs)

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
  guint i, j;

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

#ifndef POWSPEC_SPLIT_TRANSFER
  g_test_add_func ("/nc/powspec/ml/cbe/extrapolation", &test_nc_powspec_ml_cbe_extrapolation);
#endif

  g_test_run ();
}

void
test_nc_powspec_ml_transfer_new_EH (TestNcPowspec *test, gconstpointer pdata)
{
  NcHIReion *reion   = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim     = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcTransferFunc *tf = NC_TRANSFER_FUNC (ncm_serialize_global_from_string ("NcTransferFuncEH"));
  NcPowspecML *ps_ml = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));

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
  NcHIReion *reion   = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim     = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcTransferFunc *tf = NC_TRANSFER_FUNC (ncm_serialize_global_from_string ("NcTransferFuncBBKS"));
  NcPowspecML *ps_ml = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));

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
  NcHIReion *reion   = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim     = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcPowspecML *ps_ml = NC_POWSPEC_ML (nc_powspec_ml_cbe_new ());

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
  ((TestNcPowspecFunc) pdata)(test, NULL);
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
  guint i, j;

  g_assert_cmpfloat (zi, >=, 0.0);
  g_assert_cmpfloat (zi, <, zf);
  g_assert_cmpfloat (kmin, >, 0.0);
  g_assert_cmpfloat (kmin, <, kmax);

  ncm_powspec_prepare (test->ps, test->model);
  ncm_powspec_get_nknots (test->ps, &Nz, &Nk);

  g_assert_cmpuint (Nz, >, 0);
  g_assert_cmpuint (Nk, >, 0);

  Nk = MIN (Nk, 1000);
  Nz = MIN (Nz, 1000);

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
                                ncm_powspec_get_reltol_spline (test->ps) * 10.0, 0.0);
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
  const gdouble zi      = ncm_powspec_get_zi (test->ps);
  const gdouble zf      = ncm_powspec_get_zf (test->ps);
  const gdouble kmin    = ncm_powspec_get_kmin (test->ps);
  const gdouble kmax    = ncm_powspec_get_kmax (test->ps);
  const gdouble reltol  = ncm_powspec_filter_get_reltol (psf);

  ncm_powspec_filter_prepare (psf, test->model);

  g_assert_cmpfloat (zi, >=, 0.0);
  g_assert_cmpfloat (zi, <, zf);
  g_assert_cmpfloat (kmin, >, 0.0);
  g_assert_cmpfloat (kmin, <, kmax);

  /* ncm_powspec_var_tophat_R integrates the table only, while the filter continues it
   * beyond both ends; at R = 1 / k_max the continuation carries 16% of sigma^2
   * (Eisenstein-Hu) and at R = 1 / k_min 1% (BBKS), while two decades inside the two agree
   * to 1e-9. Compare there. */
  {
    const gdouble r_min = 100.0 * ncm_powspec_filter_get_r_min (psf);
    const gdouble r_max = ncm_powspec_filter_get_r_max (psf) / 100.0;
    gint i, j;

    for (i = 0; i < 10; i++)
    {
      const gdouble z = zi + (zf - zi) / (100.0 - 1.0) * i;

      for (j = 0; j < 10; j++)
      {
        const gdouble lnR    = log (r_min) + log (r_max / r_min) / (10.0 - 1.0) * j;
        const gdouble R      = exp (lnR);
        const gdouble var0   = ncm_powspec_var_tophat_R (test->ps, test->model, reltol, z, R);
        const gdouble sigma0 = ncm_powspec_sigma_tophat_R (test->ps, test->model, reltol, z, R);
        const gdouble var1   = ncm_powspec_filter_eval_var_lnr (psf, z, lnR);
        const gdouble sigma1 = ncm_powspec_filter_eval_sigma_lnr (psf, z, lnR);

        if (var0 > 1.0e-4)
        {
          ncm_assert_cmpdouble_e (var0, ==, var1, reltol * 10.0, 0.0);
          ncm_assert_cmpdouble_e (sigma0, ==, sigma1, reltol * 10.0, 0.0);
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
  const gdouble reltol  = ncm_powspec_corr3d_get_reltol (psc);

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
        const gdouble xi0 = ncm_powspec_corr3d (test->ps, test->model, reltol, z, R);
        const gdouble xi1 = ncm_powspec_corr3d_eval_xi_lnr (psc, z, lnR);

        ncm_assert_cmpdouble_e (xi0, ==, xi1, reltol * 10.0, 0.0);
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

/*
 * Outside the range CLASS computed, P / P_EH continues as a power law from the nearest edge,
 * with the mean slope of ln (P / P_EH) over the decade of computed modes next to it: P is
 * continuous at the edges, ln (P / P_EH) is linear in ln k beyond them, and the derivative
 * in z is that of the extrapolated P.
 */
void
test_nc_powspec_ml_cbe_extrapolation (void)
{
  NcHIReion *reion   = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim     = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcmModel *model    = NCM_MODEL (cosmo);
  NcPowspecMLCBE *ps = nc_powspec_ml_cbe_new ();
  NcmPowspec *cbe    = NCM_POWSPEC (ps);
  NcTransferFunc *tf = nc_transfer_func_eh_new ();
  NcmPowspec *eh     = NCM_POWSPEC (nc_powspec_ml_transfer_new (tf));
  NcmSpline2d *lnPk;
  NcmVector *lnk_v;
  gdouble lnk_edge[2];
  guint e, i;

  ncm_powspec_set_kmax (cbe, 1.0e3);
  ncm_powspec_set_kmax (eh, 1.0e3);
  ncm_powspec_prepare (cbe, model);
  ncm_powspec_prepare (eh, model);

  lnPk        = nc_cbe_get_matter_ps (nc_powspec_ml_cbe_peek_cbe (ps));
  lnk_v       = ncm_spline2d_peek_xv (lnPk);
  lnk_edge[0] = ncm_vector_get (lnk_v, 0);
  lnk_edge[1] = ncm_vector_get (lnk_v, ncm_vector_len (lnk_v) - 1);

  for (e = 0; e < 2; e++)
  {
    const gdouble lnk_e = lnk_edge[e];
    const gdouble lnk_1 = lnk_e + ((e == 0) ? M_LN10 : -M_LN10);
    const gdouble out   = (e == 0) ? -1.0 : 1.0;

    for (i = 0; i < 3; i++)
    {
      const gdouble z     = 0.1 + 0.75 * i;
      const gdouble lnr_e = ncm_spline2d_eval (lnPk, lnk_e, z) - log (ncm_powspec_eval (eh, model, z, exp (lnk_e)));
      const gdouble lnr_1 = ncm_spline2d_eval (lnPk, lnk_1, z) - log (ncm_powspec_eval (eh, model, z, exp (lnk_1)));
      const gdouble beta  = (lnr_e - lnr_1) / (lnk_e - lnk_1);
      guint j;

      ncm_assert_cmpdouble_e (ncm_powspec_eval (cbe, model, z, exp (lnk_e + out * 1.0e-10)), ==,
                              ncm_powspec_eval (cbe, model, z, exp (lnk_e - out * 1.0e-10)), 1.0e-8, 0.0);

      for (j = 1; j <= 4; j++)
      {
        const gdouble lnk  = lnk_e + out * 0.5 * j * M_LN10;
        const gdouble k    = exp (lnk);
        const gdouble P    = ncm_powspec_eval (cbe, model, z, k);
        const gdouble lnr  = log (P) - log (ncm_powspec_eval (eh, model, z, k));
        const gdouble dz   = 1.0e-4;
        const gdouble fd_z = (ncm_powspec_eval (cbe, model, z + dz, k) - ncm_powspec_eval (cbe, model, z - dz, k)) / (2.0 * dz);

        ncm_assert_cmpdouble_e (lnr, ==, lnr_e + beta * (lnk - lnk_e), 1.0e-10, 1.0e-10);
        ncm_assert_cmpdouble_e (ncm_powspec_deriv_z (cbe, model, z, k), ==, fd_z, 1.0e-7, 0.0);
      }
    }
  }

  ncm_spline2d_free (lnPk);
  ncm_powspec_free (eh);
  nc_transfer_func_free (tf);
  nc_powspec_ml_free (NC_POWSPEC_ML (ps));
  ncm_model_free (NCM_MODEL (cosmo));
  ncm_model_free (NCM_MODEL (reion));
  ncm_model_free (NCM_MODEL (prim));
}


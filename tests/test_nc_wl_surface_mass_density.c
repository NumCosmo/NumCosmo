/***************************************************************************
 *            test_nc_wl_surface_mass_density.c
 *
 *  Tue Feb 13 23:45:28 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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

typedef struct _TestNcWLSurfaceMassDensity
{
  NcWLSurfaceMassDensity *smd;
  NcHaloDensityProfile *dp;
  NcHaloMassSummary *hms;
  NcHICosmo *cosmo;
  gdouble R1, R2, R3;
  gdouble zs, zl, zc;
} TestNcWLSurfaceMassDensity;

void test_nc_wl_surface_mass_density_new (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_new_Okp (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_new_Okn (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_sigma (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_sigma_mean (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_sigma_critical (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_convergence (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_shear (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_reduced_shear (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_reduced_shear_array (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_reduced_shear_infinity (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);
void test_nc_wl_surface_mass_density_free (TestNcWLSurfaceMassDensity *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/wl_surface_mass_density/nfw/sigma", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_sigma,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/sigma_mean", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_sigma_mean,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/sigma_critical", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_sigma_critical,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/convergence", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_convergence,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/shear", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_shear,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/reduced_shear", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_reduced_shear,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/reduced_shear/array", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_reduced_shear_array,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/Okp/nfw/reduced_shear/array", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new_Okp,
              &test_nc_wl_surface_mass_density_reduced_shear_array,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/Okn/nfw/reduced_shear/array", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new_Okn,
              &test_nc_wl_surface_mass_density_reduced_shear_array,
              &test_nc_wl_surface_mass_density_free);
  g_test_add ("/nc/wl_surface_mass_density/nfw/reduced_shear_infinity", TestNcWLSurfaceMassDensity, NULL,
              &test_nc_wl_surface_mass_density_new,
              &test_nc_wl_surface_mass_density_reduced_shear_infinity,
              &test_nc_wl_surface_mass_density_free);

  g_test_run ();
}

void
test_nc_wl_surface_mass_density_free (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_halo_density_profile_free, test->dp);
  NCM_TEST_FREE (nc_halo_mass_summary_free, test->hms);
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
  NCM_TEST_FREE (nc_wl_surface_mass_density_free, test->smd);
}

void
test_nc_wl_surface_mass_density_new (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist            = nc_distance_new (3.0);
  NcHaloMassSummary *hms      = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);

  g_assert_true (smd != NULL);

  test->cosmo = cosmo;
  test->dp    = dp;
  test->hms   = hms;
  test->smd   = smd;
  test->R1    = 0.3; /* Mpc */
  test->R3    = 10.0;
  g_assert_true (NC_IS_HALO_CM_PARAM (hms));
  g_assert_true (NC_IS_HALO_DENSITY_PROFILE_NFW (dp));

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 0.0, NULL);

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 15.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "cDelta",  4.0, NULL);

  test->zc = 1.0;
  test->zl = 1.0;
  test->zs = 1.5;

  nc_distance_free (dist);

  nc_wl_surface_mass_density_prepare (smd, cosmo);
}

void
test_nc_wl_surface_mass_density_new_Okp (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist            = nc_distance_new (3.0);
  NcHaloMassSummary *hms      = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);

  g_assert_true (smd != NULL);

  test->cosmo = cosmo;
  test->hms   = hms;
  test->dp    = dp;
  test->smd   = smd;
  test->R1    = 0.3; /* Mpc */
  test->R3    = 10.0;
  g_assert_true (NC_IS_HALO_CM_PARAM (hms));
  g_assert_true (NC_IS_HALO_DENSITY_PROFILE_NFW (dp));

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 0.1, NULL);

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 15.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "cDelta",  4.0, NULL);

  test->zc = 1.0;
  test->zl = 1.0;
  test->zs = 1.5;

  nc_distance_free (dist);

  nc_wl_surface_mass_density_prepare (smd, cosmo);
}

void
test_nc_wl_surface_mass_density_new_Okn (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist            = nc_distance_new (3.0);
  NcHaloMassSummary *hms      = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);

  g_assert_true (smd != NULL);

  test->cosmo = cosmo;
  test->hms   = hms;
  test->dp    = dp;
  test->smd   = smd;
  test->R1    = 0.3; /* Mpc */
  test->R3    = 10.0;
  g_assert_true (NC_IS_HALO_CM_PARAM (hms));
  g_assert_true (NC_IS_HALO_DENSITY_PROFILE_NFW (dp));

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", -0.1, NULL);

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 15.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "cDelta",  4.0, NULL);

  test->zc = 1.0;
  test->zl = 1.0;
  test->zs = 1.5;

  nc_distance_free (dist);

  nc_wl_surface_mass_density_prepare (smd, cosmo);
}

void
test_nc_wl_surface_mass_density_sigma (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;

  test->R2 = nc_halo_density_profile_r_s (dp, cosmo, test->zc);

  gdouble sig1 = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, test->R1, test->zc);
  gdouble sig2 = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, test->R2, test->zc);
  gdouble sig3 = nc_wl_surface_mass_density_sigma (smd, dp, cosmo, test->R3, test->zc);

  ncm_assert_cmpdouble_e (sig1, ==, 6.349089e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (sig2, ==, 5.241061e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (sig3, ==, 1.861740e+12, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_sigma_mean (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;

  test->R2 = nc_halo_density_profile_r_s (dp, cosmo, test->zc);

  gdouble sig1 = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, test->R1, test->zc);
  gdouble sig2 = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, test->R2, test->zc);
  gdouble sig3 = nc_wl_surface_mass_density_sigma_mean (smd, dp, cosmo, test->R3, test->zc);

  ncm_assert_cmpdouble_e (sig1, ==, 1.116721e+15, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (sig2, ==, 9.649405e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (sig3, ==, 1.063058e+13, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_sigma_critical (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcWLSurfaceMassDensity *smd = test->smd;

  gdouble sig_crit = nc_wl_surface_mass_density_sigma_critical (smd, cosmo, test->zs, test->zl, test->zc);

  ncm_assert_cmpdouble_e (sig_crit, ==, 4.145043e+15, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_convergence (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;

  test->R2 = nc_halo_density_profile_r_s (dp, cosmo, test->zc);

  gdouble k1 = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, test->R1, test->zs, test->zl, test->zc);
  gdouble k2 = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, test->R2, test->zs, test->zl, test->zc);
  gdouble k3 = nc_wl_surface_mass_density_convergence (smd, dp, cosmo, test->R3, test->zs, test->zl, test->zc);

  ncm_assert_cmpdouble_e (k1, ==, 0.153173057,   1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k2, ==, 0.126441657,   1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k3, ==, 0.00044914867, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_shear (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;

  test->R2 = nc_halo_density_profile_r_s (dp, cosmo, test->zc);

  gdouble k1 = nc_wl_surface_mass_density_shear (smd, dp, cosmo, test->R1, test->zs, test->zl, test->zc);
  gdouble k2 = nc_wl_surface_mass_density_shear (smd, dp, cosmo, test->R2, test->zs, test->zl, test->zc);
  gdouble k3 = nc_wl_surface_mass_density_shear (smd, dp, cosmo, test->R3, test->zs, test->zl, test->zc);

  ncm_assert_cmpdouble_e (k1, ==, 0.1162381196,  1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k2, ==, 0.1063522168,  1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k3, ==, 0.00211550005, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_reduced_shear (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;

  test->R2 = nc_halo_density_profile_r_s (dp, cosmo, test->zc);

  gdouble k1 = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, test->R1, test->zs, test->zl, test->zc);
  gdouble k2 = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, test->R2, test->zs, test->zl, test->zc);
  gdouble k3 = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, test->R3, test->zs, test->zl, test->zc);

  ncm_assert_cmpdouble_e (k1, ==, 0.137263133389,  1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k2, ==, 0.121745980289,  1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k3, ==, 0.0021164506538, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_reduced_shear_infinity (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;

  test->R2 = nc_halo_density_profile_r_s (dp, cosmo, test->zc);

  gdouble k1 = nc_wl_surface_mass_density_reduced_shear_infinity (smd, dp, cosmo, test->R1, test->zs, test->zl, test->zc);
  gdouble k2 = nc_wl_surface_mass_density_reduced_shear_infinity (smd, dp, cosmo, test->R2, test->zs, test->zl, test->zc);
  gdouble k3 = nc_wl_surface_mass_density_reduced_shear_infinity (smd, dp, cosmo, test->R3, test->zs, test->zl, test->zc);

  ncm_assert_cmpdouble_e (k1, ==, 0.13726313,   1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k2, ==, 0.12174598,   1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (k3, ==, 0.0021164506, 1.0e-5, 0.0);
}

void
test_nc_wl_surface_mass_density_reduced_shear_array (TestNcWLSurfaceMassDensity *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcHaloDensityProfile *dp    = test->dp;
  NcWLSurfaceMassDensity *smd = test->smd;
  const gint nR               = g_test_rand_int_range (40, 50);
  const gint nzs              = g_test_rand_int_range (40, 50);
  GArray *R                   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), nR);
  GArray *zs                  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), nzs);
  gint i;

  g_array_set_size (R, nR);
  g_array_set_size (zs, nzs);

  for (i = 0; i < nR; i++)
    g_array_index (R, gdouble, i) = exp (g_test_rand_double_range (log (1.0e-2), log (1.0e2)));

  for (i = 0; i < nzs; i++)
    g_array_index (zs, gdouble, i) = g_test_rand_double_range (0.05, 1.5);

  {
    GArray *res = nc_wl_surface_mass_density_reduced_shear_array (smd, dp, cosmo, R, 1.0, 1.0, zs, test->zl, test->zc);

    for (i = 0; i < nR; i++)
    {
      const gdouble R_i = g_array_index (R, gdouble, i);
      gint j;

      for (j = 0; j < nzs; j++)
      {
        const gdouble zs_j = g_array_index (zs, gdouble, j);
        const gdouble rs_a = g_array_index (res, gdouble, nzs * i + j);
        const gdouble rs_d = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, R_i, zs_j, test->zl, test->zc);

        ncm_assert_cmpdouble_e (rs_a, ==, rs_d,  1.0e-8, 0.0);
      }
    }

    g_array_unref (res);
  }
}


/***************************************************************************
 *            test_nc_density_profile_nfw.c
 *
 *  Mon Feb 12 23:04:16 2017
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

typedef struct _TestNcDensityProfileNFW
{
  NcDensityProfileNFW *dp;
  NcHICosmo *cosmo;
  gdouble z;
  gdouble R1, R2, R3;
} TestNcDensityProfileNFW;

void test_nc_density_profile_nfw_new (TestNcDensityProfileNFW *test, gconstpointer pdata);
void test_nc_density_profile_nfw_eval_density (TestNcDensityProfileNFW *test, gconstpointer pdata);
void test_nc_density_profile_nfw_integral_density_los (TestNcDensityProfileNFW *test, gconstpointer pdata);
void test_nc_density_profile_nfw_integral_density_2d (TestNcDensityProfileNFW *test, gconstpointer pdata);
void test_nc_density_profile_nfw_scale_radius (TestNcDensityProfileNFW *test, gconstpointer pdata);
void test_nc_density_profile_nfw_free (TestNcDensityProfileNFW *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
	g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/density_profile_nfw/eval_density", TestNcDensityProfileNFW, NULL,
              &test_nc_density_profile_nfw_new,
              &test_nc_density_profile_nfw_eval_density,
              &test_nc_density_profile_nfw_free);
  g_test_add ("/nc/density_profile_nfw/integral_density_los", TestNcDensityProfileNFW, NULL,
              &test_nc_density_profile_nfw_new,
              &test_nc_density_profile_nfw_integral_density_los,
              &test_nc_density_profile_nfw_free);
  g_test_add ("/nc/density_profile_nfw/integral_density_2d", TestNcDensityProfileNFW, NULL,
              &test_nc_density_profile_nfw_new,
              &test_nc_density_profile_nfw_integral_density_2d,
              &test_nc_density_profile_nfw_free);
  g_test_add ("/nc/density_profile_nfw/scale_radius", TestNcDensityProfileNFW, NULL,
              &test_nc_density_profile_nfw_new,
              &test_nc_density_profile_nfw_scale_radius,
              &test_nc_density_profile_nfw_free);

  g_test_run ();
}

void
test_nc_density_profile_nfw_free (TestNcDensityProfileNFW *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_density_profile_free, NC_DENSITY_PROFILE (test->dp));
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
}

void
test_nc_density_profile_nfw_new (TestNcDensityProfileNFW *test, gconstpointer pdata)
{
  NcHICosmo *cosmo        = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcDistance *dist        = nc_distance_new (3.0);
  NcDensityProfileNFW *dp = NC_DENSITY_PROFILE_NFW (nc_density_profile_nfw_new ());
  
  g_assert (dp != NULL);

  test->cosmo = cosmo;
  test->dp    = dp;
  test->z     = 1.0;
  test->R1    = 0.3; /* Mpc */
  test->R3    = 10.0;
	
  g_assert (NC_IS_DENSITY_PROFILE_NFW (dp));

	ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
	
	nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo));
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 0.0);
  
  ncm_model_param_set_by_name (NCM_MODEL (test->dp), "MDelta",  1.0e15);
  ncm_model_param_set_by_name (NCM_MODEL (test->dp), "cDelta",  4.0);

  nc_distance_free (dist);
}

void
test_nc_density_profile_nfw_eval_density (TestNcDensityProfileNFW *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDensityProfileNFW *dp_nfw = test->dp;
	NcDensityProfile *dp        = NC_DENSITY_PROFILE (dp_nfw);

  test->R2 = nc_density_profile_scale_radius (dp, cosmo, test->z);

  gdouble rho1 = nc_density_profile_eval_density (dp, cosmo, test->R1, test->z);
  gdouble rho2 = nc_density_profile_eval_density (dp, cosmo, test->R2, test->z);
  gdouble rho3 = nc_density_profile_eval_density (dp, cosmo, test->R3, test->z);

  ncm_assert_cmpdouble_e (rho1, ==, 7.670479e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (rho2, ==, 5.557793e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (rho3, ==, 9.171099e+10, 1.0e-5, 0.0);
}

void
test_nc_density_profile_nfw_integral_density_los (TestNcDensityProfileNFW *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDensityProfileNFW *dp_nfw = test->dp;
	NcDensityProfile *dp        = NC_DENSITY_PROFILE (dp_nfw);
  
  test->R2 = nc_density_profile_scale_radius (dp, cosmo, test->z);

  gdouble I1 = nc_density_profile_integral_density_los (dp, cosmo, test->R1, test->z);
  gdouble I2 = nc_density_profile_integral_density_los (dp, cosmo, test->R2, test->z);
  gdouble I3 = nc_density_profile_integral_density_los (dp, cosmo, test->R3, test->z);

  ncm_assert_cmpdouble_e (I1, ==, 3.174544e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (I2, ==, 2.620530e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (I3, ==, 9.308702e+11, 1.0e-5, 0.0);
}

void
test_nc_density_profile_nfw_integral_density_2d (TestNcDensityProfileNFW *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDensityProfileNFW *dp_nfw = test->dp;
	NcDensityProfile *dp        = NC_DENSITY_PROFILE (dp_nfw);

	test->R2 = nc_density_profile_scale_radius (dp, cosmo, test->z);

  gdouble I1 = nc_density_profile_integral_density_2d (dp, cosmo, test->R1, test->z);
  gdouble I2 = nc_density_profile_integral_density_2d (dp, cosmo, test->R2, test->z);
  gdouble I3 = nc_density_profile_integral_density_2d (dp, cosmo, test->R3, test->z);

  ncm_assert_cmpdouble_e (I1, ==, 4.018472e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (I2, ==, 4.824703e+14, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (I3, ==, 4.250409e+15, 1.0e-5, 0.0);
}

void
test_nc_density_profile_nfw_scale_radius (TestNcDensityProfileNFW *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDensityProfileNFW *dp_nfw = test->dp;
	NcDensityProfile *dp        = NC_DENSITY_PROFILE (dp_nfw);
  
  gdouble rs = nc_density_profile_scale_radius (dp, cosmo, test->z);

  ncm_assert_cmpdouble_e (rs, ==, 0.353629, 1.0e-5, 0.0);
}

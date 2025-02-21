/***************************************************************************
 *            test_nc_distance.c
 *
 *  Tue Feb 20 16:05:43 2017
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

typedef struct _TestNcDistance
{
  NcDistance *dist;
  NcHICosmo *cosmo;
  gdouble z1, z2, z3;
  guint ntests;
} TestNcDistance;

void test_nc_distance_new (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_new_spherical (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_new_hyperbolic (TestNcDistance *test, gconstpointer pdata);

void test_nc_distance_new_no_lambda (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_new_no_lambda_spherical (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_new_no_lambda_hyperbolic (TestNcDistance *test, gconstpointer pdata);

void test_nc_distance_comoving (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_transverse (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_angular_diameter (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_DH_r (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_DA_r (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_Dt_r (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_comoving_z_to_infinity (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_transverse_z_to_infinity (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_comoving_z1_z2 (TestNcDistance *test, gconstpointer pdata);

void test_nc_distance_comoving_vector (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_transverse_vector (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_luminosity_vector (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_angular_diameter_vector (TestNcDistance *test, gconstpointer pdata);

void test_nc_distance_comoving_no_lambda (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_transverse_no_lambda (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_luminosity_no_lambda (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_angular_diameter_no_lambda (TestNcDistance *test, gconstpointer pdata);

void test_nc_distance_free (TestNcDistance *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/distance/comoving", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_comoving,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/transverse", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_transverse,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/angular_diameter", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_angular_diameter,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/DH_r", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_DH_r,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/DA_r", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_DA_r,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/Dt_r", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_Dt_r,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/comoving_z_to_infinity", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_comoving_z_to_infinity,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/transverse_to_infinity", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_transverse_z_to_infinity,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/comoving_z1_z2", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_comoving_z1_z2,
              &test_nc_distance_free);

  /*
   * Testing vectorized distances
   */
  g_test_add ("/nc/distance/flat/comoving_vector", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_comoving_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/spherical/comoving_vector", TestNcDistance, NULL,
              &test_nc_distance_new_spherical,
              &test_nc_distance_comoving_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/hyperbolic/comoving_vector", TestNcDistance, NULL,
              &test_nc_distance_new_hyperbolic,
              &test_nc_distance_comoving_vector,
              &test_nc_distance_free);

  g_test_add ("/nc/distance/flat/transverse_vector", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_transverse_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/spherical/transverse_vector", TestNcDistance, NULL,
              &test_nc_distance_new_spherical,
              &test_nc_distance_transverse_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/hyperbolic/transverse_vector", TestNcDistance, NULL,
              &test_nc_distance_new_hyperbolic,
              &test_nc_distance_transverse_vector,
              &test_nc_distance_free);

  g_test_add ("/nc/distance/flat/luminosity_vector", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_luminosity_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/spherical/luminosity_vector", TestNcDistance, NULL,
              &test_nc_distance_new_spherical,
              &test_nc_distance_luminosity_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/hyperbolic/luminosity_vector", TestNcDistance, NULL,
              &test_nc_distance_new_hyperbolic,
              &test_nc_distance_luminosity_vector,
              &test_nc_distance_free);

  g_test_add ("/nc/distance/flat/angular_diameter_vector", TestNcDistance, NULL,
              &test_nc_distance_new,
              &test_nc_distance_angular_diameter_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/spherical/angular_diameter_vector", TestNcDistance, NULL,
              &test_nc_distance_new_spherical,
              &test_nc_distance_angular_diameter_vector,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/hyperbolic/angular_diameter_vector", TestNcDistance, NULL,
              &test_nc_distance_new_hyperbolic,
              &test_nc_distance_angular_diameter_vector,
              &test_nc_distance_free);

  /*
   * Testing the no-lambda analytical solution
   */
  g_test_add ("/nc/distance/no_lambda/flat/comoving", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda,
              &test_nc_distance_comoving_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/flat/transverse", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda,
              &test_nc_distance_transverse_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/flat/luminosity", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda,
              &test_nc_distance_luminosity_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/flat/angular_diameter", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda,
              &test_nc_distance_angular_diameter_no_lambda,
              &test_nc_distance_free);

  g_test_add ("/nc/distance/no_lambda/spherical/comoving", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_spherical,
              &test_nc_distance_comoving_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/spherical/transverse", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_spherical,
              &test_nc_distance_transverse_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/spherical/luminosity", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_spherical,
              &test_nc_distance_luminosity_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/spherical/angular_diameter", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_spherical,
              &test_nc_distance_angular_diameter_no_lambda,
              &test_nc_distance_free);

  g_test_add ("/nc/distance/no_lambda/hyperbolic/comoving", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_hyperbolic,
              &test_nc_distance_comoving_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/hyperbolic/transverse", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_hyperbolic,
              &test_nc_distance_transverse_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/hyperbolic/luminosity", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_hyperbolic,
              &test_nc_distance_luminosity_no_lambda,
              &test_nc_distance_free);
  g_test_add ("/nc/distance/no_lambda/hyperbolic/angular_diameter", TestNcDistance, NULL,
              &test_nc_distance_new_no_lambda_hyperbolic,
              &test_nc_distance_angular_diameter_no_lambda,
              &test_nc_distance_free);

  g_test_run ();
}

void
test_nc_distance_free (TestNcDistance *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
  NCM_TEST_FREE (nc_distance_free, test->dist);
}

void
test_nc_distance_new (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (6.0);

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  test->cosmo  = cosmo;
  test->dist   = dist;
  test->z1     = 0.5;
  test->z2     = 2.5;
  test->z3     = 5.0;
  test->ntests = 10000;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 0.0, NULL);

  nc_distance_prepare (dist, cosmo);
}

void
test_nc_distance_new_spherical (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (6.0);

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  test->cosmo  = cosmo;
  test->dist   = dist;
  test->z1     = 0.5;
  test->z2     = 2.5;
  test->z3     = 5.0;
  test->ntests = 10000;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 1.0e-1, NULL);

  nc_distance_prepare (dist, cosmo);
}

void
test_nc_distance_new_hyperbolic (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (6.0);

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  test->cosmo  = cosmo;
  test->dist   = dist;
  test->z1     = 0.5;
  test->z2     = 2.5;
  test->z3     = 5.0;
  test->ntests = 10000;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", -1.0e-1, NULL);

  nc_distance_prepare (dist, cosmo);
}

void
test_nc_distance_new_no_lambda (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (6.0);

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  test->cosmo  = cosmo;
  test->dist   = dist;
  test->z1     = 0.5;
  test->z2     = 2.5;
  test->z3     = 5.0;
  test->ntests = 10000;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.05);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   1.0 - nc_hicosmo_Omega_t0 (test->cosmo));

  nc_distance_prepare (dist, cosmo);
}

void
test_nc_distance_new_no_lambda_spherical (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (6.0);

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  test->cosmo  = cosmo;
  test->dist   = dist;
  test->z1     = 0.5;
  test->z2     = 2.5;
  test->z3     = 5.0;
  test->ntests = 10000;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.05);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   1.0 - nc_hicosmo_Omega_t0 (test->cosmo) + 1.0e-1);

  nc_distance_prepare (dist, cosmo);
}

void
test_nc_distance_new_no_lambda_hyperbolic (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (6.0);

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  test->cosmo  = cosmo;
  test->dist   = dist;
  test->z1     = 0.5;
  test->z2     = 2.5;
  test->z3     = 5.0;
  test->ntests = 10000;

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.05);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   1.0 - nc_hicosmo_Omega_t0 (test->cosmo) - 1.0e-1);

  nc_distance_prepare (dist, cosmo);
}

void
test_nc_distance_comoving (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;

  gdouble d1 = nc_distance_comoving (dist, cosmo, test->z1);
  gdouble d2 = nc_distance_comoving (dist, cosmo, test->z2);
  gdouble d3 = nc_distance_comoving (dist, cosmo, test->z3);

  ncm_assert_cmpdouble_e (d1, ==, 0.440963874582, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d2, ==, 1.36043818854, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d3, ==, 1.81495826687, 1.0e-5, 0.0);
}

void
test_nc_distance_transverse (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;

  gdouble d1 = nc_distance_transverse (dist, cosmo, test->z1);
  gdouble d2 = nc_distance_transverse (dist, cosmo, test->z2);
  gdouble d3 = nc_distance_transverse (dist, cosmo, test->z3);

  ncm_assert_cmpdouble_e (d1, ==, 0.440963874582, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d2, ==, 1.36043818854, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d3, ==, 1.81495826687, 1.0e-5, 0.0);
}

void
test_nc_distance_angular_diameter (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;

  gdouble d1 = nc_distance_angular_diameter (dist, cosmo, test->z1);
  gdouble d2 = nc_distance_angular_diameter (dist, cosmo, test->z2);
  gdouble d3 = nc_distance_angular_diameter (dist, cosmo, test->z3);

  ncm_assert_cmpdouble_e (d1, ==, 0.293975916388, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d2, ==, 0.388696625296, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d3, ==, 0.302493044478, 1.0e-5, 0.0);
}

void
test_nc_distance_DH_r (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    const gdouble z    = g_test_rand_double_range (1.0e-1, 3.0);
    const gdouble DH_r = nc_distance_DH_r (dist, cosmo, z);
    const gdouble r_zd = nc_distance_r_zd (dist, cosmo);
    const gdouble E    = nc_hicosmo_E (cosmo, z);

    ncm_assert_cmpdouble_e ((1.0 / E) / r_zd, ==, DH_r, 1.0e-15, 0.0);
  }
}

void
test_nc_distance_DA_r (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    const gdouble z    = g_test_rand_double_range (1.0e-1, 3.0);
    const gdouble DA_r = nc_distance_DA_r (dist, cosmo, z);
    const gdouble r_zd = nc_distance_r_zd (dist, cosmo);
    const gdouble DA   = nc_distance_angular_diameter (dist, cosmo, z);

    ncm_assert_cmpdouble_e (DA / r_zd, ==, DA_r, 1.0e-15, 0.0);
  }
}

void
test_nc_distance_Dt_r (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    const gdouble z    = g_test_rand_double_range (1.0e-1, 3.0);
    const gdouble Dt_r = nc_distance_Dt_r (dist, cosmo, z);
    const gdouble r_zd = nc_distance_r_zd (dist, cosmo);
    const gdouble Dt   = nc_distance_transverse (dist, cosmo, z);

    ncm_assert_cmpdouble_e (Dt / r_zd, ==, Dt_r, 1.0e-15, 0.0);
  }
}

void
test_nc_distance_comoving_z_to_infinity (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;

  gdouble d1 = nc_distance_comoving_z_to_infinity (dist, cosmo, test->z1);
  gdouble d2 = nc_distance_comoving_z_to_infinity (dist, cosmo, test->z2);
  gdouble d3 = nc_distance_comoving_z_to_infinity (dist, cosmo, test->z3);

  ncm_assert_cmpdouble_e (d1, ==, 2.80328046107, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d2, ==, 1.88380614789, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d3, ==, 1.42928606871, 1.0e-5, 0.0);
}

void
test_nc_distance_transverse_z_to_infinity (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;

  gdouble d1 = nc_distance_transverse_z_to_infinity (dist, cosmo, test->z1);
  gdouble d2 = nc_distance_transverse_z_to_infinity (dist, cosmo, test->z2);
  gdouble d3 = nc_distance_transverse_z_to_infinity (dist, cosmo, test->z3);

  ncm_assert_cmpdouble_e (d1, ==, 2.80328046107, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d2, ==, 1.88380614789, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (d3, ==, 1.42928606871, 1.0e-5, 0.0);
}

void
test_nc_distance_comoving_z1_z2 (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  gdouble d1       = nc_distance_comoving (dist, cosmo, test->z1);
  gdouble d2       = nc_distance_comoving (dist, cosmo, test->z2);
  gdouble d3       = nc_distance_comoving (dist, cosmo, test->z3);
  gdouble d12      = nc_distance_comoving_z1_z2 (dist, cosmo, test->z1, test->z2);
  gdouble d13      = nc_distance_comoving_z1_z2 (dist, cosmo, test->z1, test->z3);
  gdouble d23      = nc_distance_comoving_z1_z2 (dist, cosmo, test->z2, test->z3);

  ncm_assert_cmpdouble_e (d12, ==, d2 - d1, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (d13, ==, d3 - d1, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (d23, ==, d3 - d2, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (d13, ==, d12 + d23, 1.0e-15, 0.0);

  {
    const gdouble delta = 1.0e-5;
    gdouble d1pdelta    = nc_distance_comoving (dist, cosmo, test->z1 + delta);
    gdouble d11pdelta   = nc_distance_comoving_z1_z2 (dist, cosmo, test->z1, test->z1 + delta);

    g_assert_cmpfloat (d1pdelta, >, d1);
    ncm_assert_cmpdouble_e (d1pdelta, ==, d1 + 1.0 / nc_hicosmo_E (cosmo, test->z1) * delta, 1.0e-10, 0.0);
    ncm_assert_cmpdouble_e (d11pdelta, ==, 1.0 / nc_hicosmo_E (cosmo, test->z1) * delta, 1.0e-5, 0.0);
  }
}

void
test_nc_distance_comoving_vector (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  GArray *z        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->ntests);
  GArray *d;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    gdouble zi = g_test_rand_double_range (1.0e-1, 8.0);

    g_array_append_val (z, zi);
  }

  for (dist->use_cache = 0; dist->use_cache < 2; dist->use_cache++)
  {
    d = nc_distance_comoving_vector (dist, cosmo, z);

    for (i = 0; i < z->len; i++)
    {
      gdouble zi = g_array_index (z, gdouble, i);
      gdouble di = g_array_index (d, gdouble, i);
      gdouble dc = nc_distance_comoving (dist, cosmo, zi);

      ncm_assert_cmpdouble_e (di, ==, dc, 1.0e-15, 0.0);
    }

    g_array_unref (d);
  }

  g_array_unref (z);
}

void
test_nc_distance_transverse_vector (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  GArray *z        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->ntests);
  GArray *d;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    gdouble zi = g_test_rand_double_range (1.0e-1, 8.0);

    g_array_append_val (z, zi);
  }

  for (dist->use_cache = 0; dist->use_cache < 2; dist->use_cache++)
  {
    d = nc_distance_transverse_vector (dist, cosmo, z);

    for (i = 0; i < z->len; i++)
    {
      gdouble zi = g_array_index (z, gdouble, i);
      gdouble di = g_array_index (d, gdouble, i);
      gdouble dc = nc_distance_transverse (dist, cosmo, zi);

      ncm_assert_cmpdouble_e (di, ==, dc, 1.0e-15, 0.0);
    }

    g_array_unref (d);
  }

  g_array_unref (z);
}

void
test_nc_distance_luminosity_vector (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  GArray *z        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->ntests);
  GArray *d;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    gdouble zi = g_test_rand_double_range (1.0e-1, 8.0);

    g_array_append_val (z, zi);
  }

  for (dist->use_cache = 0; dist->use_cache < 2; dist->use_cache++)
  {
    d = nc_distance_luminosity_vector (dist, cosmo, z);

    for (i = 0; i < z->len; i++)
    {
      gdouble zi = g_array_index (z, gdouble, i);
      gdouble di = g_array_index (d, gdouble, i);
      gdouble dc = nc_distance_luminosity (dist, cosmo, zi);

      ncm_assert_cmpdouble_e (di, ==, dc, 1.0e-15, 0.0);
    }

    g_array_unref (d);
  }

  g_array_unref (z);
}

void
test_nc_distance_angular_diameter_vector (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  GArray *z        = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->ntests);
  GArray *d;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    gdouble zi = g_test_rand_double_range (1.0e-1, 8.0);

    g_array_append_val (z, zi);
  }

  for (dist->use_cache = 0; dist->use_cache < 2; dist->use_cache++)
  {
    d = nc_distance_angular_diameter_vector (dist, cosmo, z);

    for (i = 0; i < z->len; i++)
    {
      gdouble zi = g_array_index (z, gdouble, i);
      gdouble di = g_array_index (d, gdouble, i);
      gdouble dc = nc_distance_angular_diameter (dist, cosmo, zi);

      ncm_assert_cmpdouble_e (di, ==, dc, 1.0e-15, 0.0);
    }

    g_array_unref (d);
  }

  g_array_unref (z);
}

static gdouble
_no_lambda_comoving (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble E2        = nc_hicosmo_E2 (cosmo, z);
  const gdouble E2Omega_r = nc_hicosmo_E2Omega_r (cosmo, z);

  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble Omega_r0 = nc_hicosmo_Omega_r0 (cosmo);

  const gdouble x = (sqrt (E2Omega_r) - sqrt (E2)) / (1.0 + z);
  const gdouble y = sqrt (Omega_r0) - 1.0;

  if (fabs (Omega_k0) < 1.0e-14)
  {
    return 2.0 / x - 2.0 / y;
  }
  else if (Omega_k0 < 0.0)
  {
    const gdouble atan_x = atan (x / sqrt (fabs (Omega_k0)));
    const gdouble atan_y = atan (y / sqrt (fabs (Omega_k0)));

    return -2.0 * (atan_x - atan_y) / sqrt (fabs (Omega_k0));
  }
  else
  {
    const gdouble xa  = x / sqrt (Omega_k0);
    const gdouble ya  = y / sqrt (Omega_k0);
    const gdouble arg = 2.0 * (xa - ya) / ((1.0 - xa) * (1.0 + ya));

    return log1p (arg) / sqrt (Omega_k0);
  }
}

void
test_nc_distance_comoving_no_lambda (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  guint i;

  for (i = 0; i < test->ntests; i++)
  {
    const gdouble z    = g_test_rand_double_range (1.0e-1, 8.0);
    const gdouble d    = nc_distance_comoving (dist, cosmo, z);
    const gdouble d_nl = _no_lambda_comoving (cosmo, z);

    ncm_assert_cmpdouble_e (d, ==, d_nl, 1.0e-10, 0.0);
  }
}

void
test_nc_distance_transverse_no_lambda (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDistance *dist            = test->dist;
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (k)
  {
    case 0:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z    = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d    = nc_distance_transverse (dist, cosmo, z);
        const gdouble d_nl = _no_lambda_comoving (cosmo, z);


        ncm_assert_cmpdouble_e (d, ==, d_nl, 1.0e-10, 0.0);
      }

      break;
    case -1:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_transverse (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dT_nl = sinh (sqrt_Omega_k0 * d_nl) / sqrt_Omega_k0;

        ncm_assert_cmpdouble_e (d, ==, dT_nl, 1.0e-10, 0.0);
      }

      break;
    case 1:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_transverse (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dT_nl = sin (sqrt_Omega_k0 * d_nl) / sqrt_Omega_k0;

        ncm_assert_cmpdouble_e (d, ==, dT_nl, 1.0e-10, 0.0);
      }

      break;
    default:
      g_assert_not_reached ();

      break;
  }
}

void
test_nc_distance_luminosity_no_lambda (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDistance *dist            = test->dist;
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (k)
  {
    case 0:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_luminosity (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dL_nl = (1.0 + z) * d_nl;

        ncm_assert_cmpdouble_e (d, ==, dL_nl, 1.0e-10, 0.0);
      }

      break;
    case -1:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_luminosity (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dL_nl = (1.0 + z) * sinh (sqrt_Omega_k0 * d_nl) / sqrt_Omega_k0;

        ncm_assert_cmpdouble_e (d, ==, dL_nl, 1.0e-10, 0.0);
      }

      break;
    case 1:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_luminosity (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dL_nl = (1.0 + z) * sin (sqrt_Omega_k0 * d_nl) / sqrt_Omega_k0;

        ncm_assert_cmpdouble_e (d, ==, dL_nl, 1.0e-10, 0.0);
      }

      break;
    default:
      g_assert_not_reached ();

      break;
  }
}

void
test_nc_distance_angular_diameter_no_lambda (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo            = test->cosmo;
  NcDistance *dist            = test->dist;
  const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
  const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
  guint i;

  switch (k)
  {
    case 0:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_angular_diameter (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dA_nl = d_nl / (1.0 + z);

        ncm_assert_cmpdouble_e (d, ==, dA_nl, 1.0e-10, 0.0);
      }

      break;
    case -1:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_angular_diameter (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dA_nl = sinh (sqrt_Omega_k0 * d_nl) / (1.0 + z) / sqrt_Omega_k0;

        ncm_assert_cmpdouble_e (d, ==, dA_nl, 1.0e-10, 0.0);
      }

      break;
    case 1:

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble z     = g_test_rand_double_range (1.0e-1, 8.0);
        const gdouble d     = nc_distance_angular_diameter (dist, cosmo, z);
        const gdouble d_nl  = _no_lambda_comoving (cosmo, z);
        const gdouble dA_nl = sin (sqrt_Omega_k0 * d_nl) / (1.0 + z) / sqrt_Omega_k0;

        ncm_assert_cmpdouble_e (d, ==, dA_nl, 1.0e-10, 0.0);
      }

      break;
    default:
      g_assert_not_reached ();

      break;
  }
}


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
} TestNcDistance;

void test_nc_distance_new (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_comoving (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_transverse (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_angular_diameter (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_DH_r (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_DA_r (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_Dt_r (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_comoving_z_to_infinity (TestNcDistance *test, gconstpointer pdata);
void test_nc_distance_transverse_z_to_infinity (TestNcDistance *test, gconstpointer pdata);
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
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");
  NcDistance *dist = nc_distance_new (6.0);
  
  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));
  
  test->cosmo = cosmo;
  test->dist  = dist;
  test->z1    = 0.5;
  test->z2    = 2.5;
  test->z3    = 5.0;
  
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo));
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 0.0);
  
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
test_nc_distance_DH_r (TestNcDistance *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NcDistance *dist = test->dist;
  gint i;

  for (i = 0; i < 1000; i++)
  {
    const gdouble z = g_test_rand_double_range (1.0e-1, 3.0);
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
  gint i;

  for (i = 0; i < 1000; i++)
  {
    const gdouble z = g_test_rand_double_range (1.0e-1, 3.0);
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
  gint i;

  for (i = 0; i < 1000; i++)
  {
    const gdouble z = g_test_rand_double_range (1.0e-1, 3.0);
    const gdouble Dt_r = nc_distance_Dt_r (dist, cosmo, z);
    const gdouble r_zd = nc_distance_r_zd (dist, cosmo);
    const gdouble Dt   = nc_distance_transverse (dist, cosmo, z);

    ncm_assert_cmpdouble_e (Dt / r_zd, ==, Dt_r, 1.0e-15, 0.0);
  }
}

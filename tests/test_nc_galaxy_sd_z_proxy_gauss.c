/***************************************************************************
 *            test_nc_galaxy_sd_z_proxy_gauss.c
 *
 *  mon May 06 23:27:20 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiolimadeoliveira@pm.me>
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

typedef struct _TestNcGalaxySDZProxyGauss
{
  NcGalaxySDZProxyGauss *sdzpg;
} TestNcGalaxySDZProxyGauss;


static void test_nc_galaxy_sd_z_proxy_gauss_new (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_free (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_basic (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_serialize (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_sigma (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_z_lim (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_true_z_lim (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);


static void test_nc_galaxy_sd_z_proxy_gauss_gen (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_gauss_integ (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/basic", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_basic,
              &test_nc_galaxy_sd_z_proxy_gauss_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/serialize", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_serialize,
              &test_nc_galaxy_sd_z_proxy_gauss_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/sigma", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_sigma,
              &test_nc_galaxy_sd_z_proxy_gauss_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/z_lim", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_z_lim,
              &test_nc_galaxy_sd_z_proxy_gauss_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/true_z_lim", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_true_z_lim,
              &test_nc_galaxy_sd_z_proxy_gauss_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/gen", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_gen,
              &test_nc_galaxy_sd_z_proxy_gauss_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/gauss/integ", TestNcGalaxySDZProxyGauss, NULL,
              &test_nc_galaxy_sd_z_proxy_gauss_new,
              &test_nc_galaxy_sd_z_proxy_gauss_integ,
              &test_nc_galaxy_sd_z_proxy_gauss_free);


  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_z_proxy_gauss_new (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = nc_galaxy_sd_z_proxy_gauss_new (1.0e-6, 10, 0.05);

  test->sdzpg = sdzpg;

  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (sdzpg));
}

static void
test_nc_galaxy_sd_z_proxy_gauss_free (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_z_proxy_gauss_free, test->sdzpg);
}

static void
test_nc_galaxy_sd_z_proxy_gauss_basic (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;
  NcGalaxySDZProxyGauss *sdzpg2;

  sdzpg2 = nc_galaxy_sd_z_proxy_gauss_ref (sdzpg);
  nc_galaxy_sd_z_proxy_gauss_clear (&sdzpg2);
  g_assert_true (sdzpg2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (sdzpg));
}

static void
test_nc_galaxy_sd_z_proxy_gauss_serialize (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;
  NcGalaxySDZProxyGauss *sdzpg_dup;
  NcmSerialize *ser;
  gchar *sdzpg_ser;
  gdouble z_min;
  gdouble z_max;
  gdouble z_min_dup;
  gdouble z_max_dup;
  gdouble true_z_min;
  gdouble true_z_max;
  gdouble true_z_min_dup;
  gdouble true_z_max_dup;

  ser       = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  sdzpg_ser = ncm_serialize_to_string (ser, G_OBJECT (sdzpg), TRUE);
  sdzpg_dup = NC_GALAXY_SD_Z_PROXY_GAUSS (ncm_serialize_from_string (ser, sdzpg_ser));

  nc_galaxy_sd_z_proxy_gauss_get_z_lim (sdzpg, &z_min, &z_max);
  nc_galaxy_sd_z_proxy_gauss_get_z_lim (sdzpg_dup, &z_min_dup, &z_max_dup);
  nc_galaxy_sd_z_proxy_get_true_z_lim (NC_GALAXY_SD_Z_PROXY (sdzpg), 0.4, &true_z_min, &true_z_max);
  nc_galaxy_sd_z_proxy_get_true_z_lim (NC_GALAXY_SD_Z_PROXY (sdzpg_dup), 0.4, &true_z_min_dup, &true_z_max_dup);

  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_gauss_get_sigma (sdzpg), ==, nc_galaxy_sd_z_proxy_gauss_get_sigma (sdzpg_dup));
  g_assert_cmpfloat (z_min, ==, z_min_dup);
  g_assert_cmpfloat (z_max, ==, z_max_dup);
  g_assert_cmpfloat (true_z_min, ==, true_z_min_dup);
  g_assert_cmpfloat (true_z_max, ==, true_z_max_dup);
  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_integ (NC_GALAXY_SD_Z_PROXY (sdzpg), 0.5, 0.5), ==, nc_galaxy_sd_z_proxy_integ (NC_GALAXY_SD_Z_PROXY (sdzpg_dup), 0.5, 0.5));

  ncm_serialize_free (ser);
  g_free (sdzpg_ser);
  nc_galaxy_sd_z_proxy_gauss_free (sdzpg_dup);
}

static void
test_nc_galaxy_sd_z_proxy_gauss_sigma (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;

  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_gauss_get_sigma (sdzpg), ==, 0.05);
}

static void
test_nc_galaxy_sd_z_proxy_gauss_z_lim (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;
  gdouble z_min;
  gdouble z_max;

  nc_galaxy_sd_z_proxy_gauss_get_z_lim (sdzpg, &z_min, &z_max);

  g_assert_cmpfloat (z_min, ==, 1.0e-6);
  g_assert_cmpfloat (z_max, ==, 10);
}

static void
test_nc_galaxy_sd_z_proxy_gauss_true_z_lim (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;
  gdouble true_z_min;
  gdouble true_z_max;

  nc_galaxy_sd_z_proxy_gauss_set_true_z_min (sdzpg, 1.0e-3);

  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_gauss_get_true_z_min (sdzpg), ==, 1.0e-3);

  nc_galaxy_sd_z_proxy_get_true_z_lim (NC_GALAXY_SD_Z_PROXY (sdzpg), 1.1, &true_z_min, &true_z_max);

  g_assert_cmpfloat (true_z_min, >, 1.0e-3);

  nc_galaxy_sd_z_proxy_gauss_set_true_z_min (sdzpg, 0.4);
  nc_galaxy_sd_z_proxy_get_true_z_lim (NC_GALAXY_SD_Z_PROXY (sdzpg), 1.0, &true_z_min, &true_z_max);

  g_assert_cmpfloat (true_z_min, ==, 0.4);
}

static void
test_nc_galaxy_sd_z_proxy_gauss_gen (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  gdouble zp;

  nc_galaxy_sd_z_proxy_gauss_set_z_lim (sdzpg, 0.0, 1.0);
  g_assert_true (!nc_galaxy_sd_z_proxy_gen (NC_GALAXY_SD_Z_PROXY (sdzpg), rng, 100, &zp));
  g_assert_true (zp >= 1.0);

  nc_galaxy_sd_z_proxy_gauss_set_z_lim (sdzpg, 0.5, 1.5);
  g_assert_true (nc_galaxy_sd_z_proxy_gen (NC_GALAXY_SD_Z_PROXY (sdzpg), rng, 1, &zp));
  g_assert_true (zp > 0.5 && zp < 1.5);

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_z_proxy_gauss_integ (TestNcGalaxySDZProxyGauss *test, gconstpointer pdata)
{
  NcGalaxySDZProxyGauss *sdzpg = test->sdzpg;
  gdouble z                    = 0.5;
  gdouble zp                   = 0.7;

  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_integ (NC_GALAXY_SD_Z_PROXY (sdzpg), z, zp), >, 0.0);
}


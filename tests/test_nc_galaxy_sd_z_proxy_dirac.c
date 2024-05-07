/***************************************************************************
 *            test_nc_data_cluster_wl.c
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

typedef struct _TestNcGalaxySDZProxyDirac
{
  NcGalaxySDZProxyDirac *sdzpd;
} TestNcGalaxySDZProxyDirac;


static void test_nc_galaxy_sd_z_proxy_dirac_new (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_dirac_free (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_dirac_gen (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata);

static void test_nc_galaxy_sd_z_proxy_dirac_integ (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_z_proxy/dirac/gen", TestNcGalaxySDZProxyDirac, NULL,
              &test_nc_galaxy_sd_z_proxy_dirac_new,
              &test_nc_galaxy_sd_z_proxy_dirac_gen,
              &test_nc_galaxy_sd_z_proxy_dirac_free);

  g_test_add ("/nc/galaxy_sd_z_proxy/dirac/integ", TestNcGalaxySDZProxyDirac, NULL,
              &test_nc_galaxy_sd_z_proxy_dirac_new,
              &test_nc_galaxy_sd_z_proxy_dirac_integ,
              &test_nc_galaxy_sd_z_proxy_dirac_free);


  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_z_proxy_dirac_new (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata)
{
  NcGalaxySDZProxyDirac *sdzpd = nc_galaxy_sd_z_proxy_dirac_new ();

  test->sdzpd = sdzpd;

  g_assert_true (NC_IS_GALAXY_SD_Z_PROXY_DIRAC (sdzpd));
}

static void
test_nc_galaxy_sd_z_proxy_dirac_free (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_z_proxy_dirac_free, test->sdzpd);
}

static void
test_nc_galaxy_sd_z_proxy_dirac_gen (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata)
{
  NcGalaxySDZProxyDirac *sdzpd = test->sdzpd;
  gdouble z                    = 0.5;
  gdouble zp                   = 0.0;

  g_assert_true (nc_galaxy_sd_z_proxy_gen (NC_GALAXY_SD_Z_PROXY (sdzpd), NULL, z, &zp));
  g_assert_cmpfloat (zp, ==, z);
}

static void
test_nc_galaxy_sd_z_proxy_dirac_integ (TestNcGalaxySDZProxyDirac *test, gconstpointer pdata)
{
  NcGalaxySDZProxyDirac *sdzpd = test->sdzpd;
  gdouble z                    = 0.5;
  gdouble zp                   = 0.0;

  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_integ (NC_GALAXY_SD_Z_PROXY (sdzpd), z, z), ==, 1.0);
  g_assert_cmpfloat (nc_galaxy_sd_z_proxy_integ (NC_GALAXY_SD_Z_PROXY (sdzpd), z, zp), ==, 0.0);
}


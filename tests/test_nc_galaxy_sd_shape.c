/***************************************************************************
 *            test_nc_galaxy_sd_shape.c
 *
 *  mon May 07 00:02:42 2024
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

typedef struct _TestNcGalaxySDShapeGauss
{
  NcGalaxySDShapeGauss *sdsg;
} TestNcGalaxySDShapeGauss;


static void test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_basic (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_get_sigma (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_gen (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_gen_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_integ_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_integ_optzs (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_integ_optzs_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_shape/gauss/basic", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_basic,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/serialize", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_serialize,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/get_sigma", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_get_sigma,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/gen", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_gen,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/gen_strong", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_gen_strong,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ_strong", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ_strong,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ_optzs", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ_optzs,
              &test_nc_galaxy_sd_shape_gauss_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ_optzs_strong", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ_optzs_strong,
              &test_nc_galaxy_sd_shape_gauss_free);


  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg = nc_galaxy_sd_shape_gauss_new ();

  nc_galaxy_sd_shape_gauss_set_sigma (sdsg, 0.05);

  test->sdsg = sdsg;

  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (sdsg));
}

static void
test_nc_galaxy_sd_shape_gauss_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_shape_gauss_free, test->sdsg);
}

static void
test_nc_galaxy_sd_shape_gauss_basic (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg = test->sdsg;
  NcGalaxySDShapeGauss *sdsg2;

  sdsg2 = nc_galaxy_sd_shape_gauss_ref (sdsg);
  nc_galaxy_sd_shape_gauss_clear (&sdsg2);
  g_assert_true (sdsg2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (sdsg));
}

static void
test_nc_galaxy_sd_shape_gauss_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 1.0;
  gdouble z                   = 0.6;
  gdouble et                  = 0.01;
  gdouble ex                  = 0.01;
  NcGalaxySDShapeGauss *sdsg_dup;
  NcmSerialize *ser;
  gchar *sdsg_ser;

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  sdsg_ser = ncm_serialize_to_string (ser, G_OBJECT (sdsg), TRUE);
  sdsg_dup = NC_GALAXY_SD_SHAPE_GAUSS (ncm_serialize_from_string (ser, sdsg_ser));

  g_assert_cmpfloat (nc_galaxy_sd_shape_gauss_get_sigma (sdsg), ==, nc_galaxy_sd_shape_gauss_get_sigma (sdsg_dup));
  g_assert_cmpfloat (nc_galaxy_sd_shape_integ (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, r, z, et, ex), ==, nc_galaxy_sd_shape_integ (NC_GALAXY_SD_SHAPE (sdsg_dup), cosmo, dp, smd, z_cluster, r, z, et, ex));

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  ncm_serialize_free (ser);
  g_free (sdsg_ser);
  nc_galaxy_sd_shape_gauss_free (sdsg_dup);
}

static void
test_nc_galaxy_sd_shape_gauss_get_sigma (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg = test->sdsg;

  g_assert_cmpfloat (nc_galaxy_sd_shape_gauss_get_sigma (sdsg), ==, 0.05);
}

static void
test_nc_galaxy_sd_shape_gauss_gen (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 1.0;
  gdouble z                   = 0.6;
  gdouble et                  = 0.0;
  gdouble ex                  = 0.0;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_galaxy_sd_shape_gen (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, rng, r, z, &et, &ex);

  g_assert_cmpfloat (et, <, nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cluster, z_cluster) + 5.0 * 0.05);
  g_assert_cmpfloat (et, >, nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cluster, z_cluster) - 5.0 * 0.05);
  g_assert_cmpfloat (ex, <, 5.0 * 0.05);
  g_assert_cmpfloat (ex, >, -5.0 * 0.05);

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_gauss_gen_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_lcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL, 200.0));
  NcDistance *dist            = nc_distance_new (10.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 0.1;
  gdouble z                   = 0.6;
  gdouble et                  = 0.0;
  gdouble ex                  = 0.0;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  ncm_model_param_set_by_name (NCM_MODEL (dp), "log10MDelta", log10 (1.0e16));
  nc_galaxy_sd_shape_gen (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, rng, r, z, &et, &ex);

  g_assert_cmpfloat (sqrt (et * et + ex * ex), <, 1.0);

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 1.0;
  gdouble z                   = 0.6;
  gdouble et                  = 0.0;
  gdouble ex                  = 0.0;
  gdouble res;

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  res = nc_galaxy_sd_shape_integ (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, r, z, et, ex);

  g_assert_cmpfloat (res, >, 0.0);

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}

static void
test_nc_galaxy_sd_shape_gauss_integ_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_lcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL, 200.0));
  NcDistance *dist            = nc_distance_new (10.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 0.1;
  gdouble z                   = 0.6;
  gdouble et                  = 0.0;
  gdouble ex                  = 0.0;
  gdouble res;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  ncm_model_param_set_by_name (NCM_MODEL (dp), "log10MDelta", log10 (1.0e16));

  res = nc_galaxy_sd_shape_integ (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, r, z, et, ex);

  g_assert_cmpfloat (res, >, 1.0e-110);

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}

static void
test_nc_galaxy_sd_shape_gauss_integ_optzs (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 1.0;
  gdouble z                   = 0.6;
  gdouble et                  = 0.0;
  gdouble ex                  = 0.0;
  gdouble res;

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  nc_galaxy_sd_shape_integ_optzs_prep (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, r);

  res = nc_galaxy_sd_shape_integ_optzs (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, z, et, ex);

  g_assert_cmpfloat (res, >, 0.0);

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}

static void
test_nc_galaxy_sd_shape_gauss_integ_optzs_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *sdsg  = test->sdsg;
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_lcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL, 200.0));
  NcDistance *dist            = nc_distance_new (10.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble r                   = 0.1;
  gdouble z                   = 0.6;
  gdouble et                  = 0.0;
  gdouble ex                  = 0.0;
  gdouble res;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  ncm_model_param_set_by_name (NCM_MODEL (dp), "log10MDelta", log10 (1.0e16));

  nc_galaxy_sd_shape_integ_optzs_prep (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, r);

  res = nc_galaxy_sd_shape_integ_optzs (NC_GALAXY_SD_SHAPE (sdsg), cosmo, dp, smd, z_cluster, z, et, ex);

  g_assert_cmpfloat (res, >, 1.0e-110);

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}


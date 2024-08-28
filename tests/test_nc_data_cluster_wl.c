/***************************************************************************
 *            test_nc_data_cluster_wl.c
 *
 *  thu December 07 13:54:21 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2023 <vitenti@uel.br>
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

typedef struct _TestNcDataClusterWL
{
  NcDataClusterWL *dcwl;
} TestNcDataClusterWL;


static void test_nc_data_cluster_wl_new_spec (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_new_gauss (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_r_cut (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_set_obs_spec (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_set_obs_gauss (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_gen_obs_spec (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_gen_obs_gauss (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_serialize_spec (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_serialize_gauss (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_gata_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_m2lnP_integ (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_fit (TestNcDataClusterWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  /* g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/serialize", TestNcDataClusterWL, NULL, */
  /*             &test_nc_data_cluster_wl_new_spec, */
  /*             &test_nc_data_cluster_wl_serialize_spec, */
  /*             &test_nc_data_cluster_wl_free); */

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/r_cut", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_r_cut,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/set_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_set_obs_spec,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_gen_obs_spec,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/m2lnP", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_gata_cluster_wl_m2lnP,
              &test_nc_data_cluster_wl_free);

  /* g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/serialize", TestNcDataClusterWL, NULL, */
  /*             &test_nc_data_cluster_wl_new_gauss, */
  /*             &test_nc_data_cluster_wl_serialize_gauss, */
  /*             &test_nc_data_cluster_wl_free); */

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/r_cut", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_r_cut,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/set_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_set_obs_gauss,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_gen_obs_gauss,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/m2lnP_integ", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_m2lnP_integ,
              &test_nc_data_cluster_wl_free);

  g_test_run ();

  return 0;
}

static void
test_nc_data_cluster_wl_new_spec (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (0.0, 5.0));
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true_dist));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new (s_dist, z_dist, p_dist);
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloPosition *hp                  = nc_halo_position_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  nc_galaxy_sd_shape_set_models (s_dist, cosmo, hp);

  test->dcwl = dcwl;

  nc_galaxy_sd_shape_free (s_dist);
  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_galaxy_sd_obs_redshift_free (z_dist);
  nc_galaxy_sd_position_free (p_dist);
  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_position_free (hp);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_new_gauss (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (0.0, 5.0));
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_gauss_new (z_true_dist));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new (s_dist, z_dist, p_dist);
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloPosition *hp                  = nc_halo_position_new (dist);

  nc_halo_position_prepare (hp, cosmo);
  nc_galaxy_sd_shape_set_models (s_dist, cosmo, hp);

  test->dcwl = dcwl;

  nc_galaxy_sd_shape_free (s_dist);
  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_galaxy_sd_obs_redshift_free (z_dist);
  nc_galaxy_sd_position_free (p_dist);
  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_position_free (hp);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_data_cluster_wl_free, test->dcwl);
}

static void
test_nc_data_cluster_wl_r_cut (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  guint ngals                 = 200;
  NcGalaxyWLObs *obs;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);
  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, hp, ngals, rng, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN);

  obs = nc_data_cluster_wl_peek_obs (test->dcwl);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra  = nc_galaxy_wl_obs_get (obs, "ra", i);
    gdouble dec = nc_galaxy_wl_obs_get (obs, "dec", i);
    gdouble theta, phi, r;

    nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);
    r = nc_halo_position_projected_radius (hp, cosmo, theta);

    g_assert_cmpfloat (r, >=, 0.3);
    g_assert_cmpfloat (r, <=, 3.0);
  }

  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
  nc_galaxy_wl_obs_free (obs);
}

static void
test_nc_data_cluster_wl_set_obs_spec (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint ngals        = 200;
  GStrv header       = g_strsplit ("ra dec z e1 e2 e_rms e_sigma", " ", -1);
  NcGalaxyWLObs *obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, ngals, header);
  NcGalaxyWLObs *obs2;
  guint i;

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra      = ncm_rng_uniform_gen (rng, 0.0, 2.0 * M_PI);
    gdouble dec     = ncm_rng_uniform_gen (rng, -M_PI / 2.0, M_PI / 2.0);
    gdouble z       = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    gdouble e1      = ncm_rng_uniform_gen (rng, -1.0, 1.0);
    gdouble e2      = ncm_rng_uniform_gen (rng, -1.0, 1.0);
    gdouble e_rms   = ncm_rng_uniform_gen (rng, 0.0, 0.3);
    gdouble e_sigma = ncm_rng_uniform_gen (rng, 0.0, 0.1);

    nc_galaxy_wl_obs_set (obs, "ra", i, ra);
    nc_galaxy_wl_obs_set (obs, "dec", i, dec);
    nc_galaxy_wl_obs_set (obs, "z", i, z);
    nc_galaxy_wl_obs_set (obs, "e1", i, e1);
    nc_galaxy_wl_obs_set (obs, "e2", i, e2);
    nc_galaxy_wl_obs_set (obs, "e_rms", i, e_rms);
    nc_galaxy_wl_obs_set (obs, "e_sigma", i, e_sigma);
  }

  nc_data_cluster_wl_set_obs (test->dcwl, obs);

  obs2 = nc_data_cluster_wl_peek_obs (test->dcwl);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra      = nc_galaxy_wl_obs_get (obs2, "ra", i);
    gdouble dec     = nc_galaxy_wl_obs_get (obs2, "dec", i);
    gdouble z       = nc_galaxy_wl_obs_get (obs2, "z", i);
    gdouble e1      = nc_galaxy_wl_obs_get (obs2, "e1", i);
    gdouble e2      = nc_galaxy_wl_obs_get (obs2, "e2", i);
    gdouble e_rms   = nc_galaxy_wl_obs_get (obs2, "e_rms", i);
    gdouble e_sigma = nc_galaxy_wl_obs_get (obs2, "e_sigma", i);

    g_assert_cmpfloat (ra, ==, nc_galaxy_wl_obs_get (obs, "ra", i));
    g_assert_cmpfloat (dec, ==, nc_galaxy_wl_obs_get (obs, "dec", i));
    g_assert_cmpfloat (z, ==, nc_galaxy_wl_obs_get (obs, "z", i));
    g_assert_cmpfloat (e1, ==, nc_galaxy_wl_obs_get (obs, "e1", i));
    g_assert_cmpfloat (e2, ==, nc_galaxy_wl_obs_get (obs, "e2", i));
    g_assert_cmpfloat (e_rms, ==, nc_galaxy_wl_obs_get (obs, "e_rms", i));
    g_assert_cmpfloat (e_sigma, ==, nc_galaxy_wl_obs_get (obs, "e_sigma", i));
  }

  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_wl_obs_free (obs2);
}

static void
test_nc_data_cluster_wl_set_obs_gauss (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint ngals        = 200;
  GStrv header       = g_strsplit ("ra dec zp zp_sigma e1 e2 e_rms e_sigma", " ", -1);
  NcGalaxyWLObs *obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, ngals, header);
  NcGalaxyWLObs *obs2;
  guint i;


  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra       = ncm_rng_uniform_gen (rng, 0.0, 2.0 * M_PI);
    gdouble dec      = ncm_rng_uniform_gen (rng, -M_PI / 2.0, M_PI / 2.0);
    gdouble zp       = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    gdouble zp_sigma = ncm_rng_uniform_gen (rng, 0.0, 0.05);
    gdouble e1       = ncm_rng_uniform_gen (rng, -1.0, 1.0);
    gdouble e2       = ncm_rng_uniform_gen (rng, -1.0, 1.0);
    gdouble e_rms    = ncm_rng_uniform_gen (rng, 0.0, 0.3);
    gdouble e_sigma  = ncm_rng_uniform_gen (rng, 0.0, 0.1);

    nc_galaxy_wl_obs_set (obs, "ra", i, ra);
    nc_galaxy_wl_obs_set (obs, "dec", i, dec);
    nc_galaxy_wl_obs_set (obs, "zp", i, zp);
    nc_galaxy_wl_obs_set (obs, "zp_sigma", i, zp_sigma);
    nc_galaxy_wl_obs_set (obs, "e1", i, e1);
    nc_galaxy_wl_obs_set (obs, "e2", i, e2);
    nc_galaxy_wl_obs_set (obs, "e_rms", i, e_rms);
    nc_galaxy_wl_obs_set (obs, "e_sigma", i, e_sigma);
  }

  nc_data_cluster_wl_set_obs (test->dcwl, obs);

  obs2 = nc_data_cluster_wl_peek_obs (test->dcwl);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra       = nc_galaxy_wl_obs_get (obs2, "ra", i);
    gdouble dec      = nc_galaxy_wl_obs_get (obs2, "dec", i);
    gdouble zp       = nc_galaxy_wl_obs_get (obs2, "zp", i);
    gdouble zp_sigma = nc_galaxy_wl_obs_get (obs2, "zp_sigma", i);
    gdouble e1       = nc_galaxy_wl_obs_get (obs2, "e1", i);
    gdouble e2       = nc_galaxy_wl_obs_get (obs2, "e2", i);
    gdouble e_rms    = nc_galaxy_wl_obs_get (obs2, "e_rms", i);
    gdouble e_sigma  = nc_galaxy_wl_obs_get (obs2, "e_sigma", i);

    g_assert_cmpfloat (ra, ==, nc_galaxy_wl_obs_get (obs, "ra", i));
    g_assert_cmpfloat (dec, ==, nc_galaxy_wl_obs_get (obs, "dec", i));
    g_assert_cmpfloat (zp, ==, nc_galaxy_wl_obs_get (obs, "zp", i));
    g_assert_cmpfloat (zp_sigma, ==, nc_galaxy_wl_obs_get (obs, "zp_sigma", i));
    g_assert_cmpfloat (e1, ==, nc_galaxy_wl_obs_get (obs, "e1", i));
    g_assert_cmpfloat (e2, ==, nc_galaxy_wl_obs_get (obs, "e2", i));
    g_assert_cmpfloat (e_rms, ==, nc_galaxy_wl_obs_get (obs, "e_rms", i));
    g_assert_cmpfloat (e_sigma, ==, nc_galaxy_wl_obs_get (obs, "e_sigma", i));
  }

  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_wl_obs_free (obs2);
}

static void
test_nc_data_cluster_wl_gen_obs_spec (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  guint ngals                 = 200;
  NcGalaxyWLObs *obs;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);
  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, hp, ngals, rng, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN);

  obs = nc_data_cluster_wl_peek_obs (test->dcwl);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra      = nc_galaxy_wl_obs_get (obs, "ra", i);
    gdouble dec     = nc_galaxy_wl_obs_get (obs, "dec", i);
    gdouble z       = nc_galaxy_wl_obs_get (obs, "z", i);
    gdouble e1      = nc_galaxy_wl_obs_get (obs, "e1", i);
    gdouble e2      = nc_galaxy_wl_obs_get (obs, "e2", i);
    gdouble e_rms   = nc_galaxy_wl_obs_get (obs, "e_rms", i);
    gdouble e_sigma = nc_galaxy_wl_obs_get (obs, "e_sigma", i);
    gdouble theta, phi, r;

    nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);
    r = nc_halo_position_projected_radius (hp, cosmo, theta);

    g_assert_cmpfloat (r, >=, 0.3);
    g_assert_cmpfloat (r, <=, 3.0);
    g_assert_cmpfloat (z, >=, 0.0);
    g_assert_cmpfloat (z, <=, 5.0);
    g_assert_cmpfloat (e1, >=, -5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e1, <=, 5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e2, >=, -5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e2, <=, 5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e_rms, ==, 0.3);
    g_assert_cmpfloat (e_sigma, ==, 0.1);
  }

  nc_galaxy_wl_obs_free (obs);
}

static void
test_nc_data_cluster_wl_gen_obs_gauss (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  guint ngals                 = 200;
  NcGalaxyWLObs *obs;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);
  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, hp, ngals, rng, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN);

  obs = nc_data_cluster_wl_peek_obs (test->dcwl);

  for (i = 0; i < ngals; i++)
  {
    gdouble ra       = nc_galaxy_wl_obs_get (obs, "ra", i);
    gdouble dec      = nc_galaxy_wl_obs_get (obs, "dec", i);
    gdouble zp       = nc_galaxy_wl_obs_get (obs, "zp", i);
    gdouble zp_sigma = nc_galaxy_wl_obs_get (obs, "zp_sigma", i);
    gdouble e1       = nc_galaxy_wl_obs_get (obs, "e1", i);
    gdouble e2       = nc_galaxy_wl_obs_get (obs, "e2", i);
    gdouble e_rms    = nc_galaxy_wl_obs_get (obs, "e_rms", i);
    gdouble e_sigma  = nc_galaxy_wl_obs_get (obs, "e_sigma", i);
    gdouble theta, phi, r;

    nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);
    r = nc_halo_position_projected_radius (hp, cosmo, theta);

    g_assert_cmpfloat (r, >=, 0.3);
    g_assert_cmpfloat (r, <=, 3.0);
    g_assert_cmpfloat (zp, >=, 0.0);
    g_assert_cmpfloat (zp, <=, 5.0);
    g_assert_cmpfloat (zp_sigma, ==, 0.05);
    g_assert_cmpfloat (e1, >=, -5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e1, <=, 5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e2, >=, -5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e2, <=, 5 * sqrt (e_rms * e_rms + e_sigma * e_sigma));
    g_assert_cmpfloat (e_rms, ==, 0.3);
    g_assert_cmpfloat (e_sigma, ==, 0.1);
  }

  nc_galaxy_wl_obs_free (obs);
}

static void
test_nc_data_cluster_wl_serialize_spec (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo             = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp     = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist             = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd  = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp           = nc_halo_position_new (dist);
  const guint ngals            = 200;
  NcmVector *m2lnP_int_gal     = ncm_vector_new (ngals);
  NcmVector *m2lnP_int_gal_dup = ncm_vector_new (ngals);
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  GVariant *dcwl_ser;
  NcDataClusterWL *dcwl_dup;
  guint i;

  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (hp), NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (hp), NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, 0.0);

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);
  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, hp, ngals, rng, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN);

  ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  dcwl_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->dcwl));
  dcwl_dup = NC_DATA_CLUSTER_WL (ncm_serialize_from_variant (ser, dcwl_ser));

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl_dup));

  nc_data_cluster_wl_eval_m2lnP (test->dcwl, cosmo, dp, smd, hp, m2lnP_int_gal);
  nc_data_cluster_wl_eval_m2lnP (dcwl_dup, cosmo, dp, smd, hp, m2lnP_int_gal_dup);

  for (i = 0; i < ngals; i++)
    g_assert_cmpfloat (ncm_vector_get (m2lnP_int_gal, i), ==, ncm_vector_get (m2lnP_int_gal_dup, i));

  ncm_vector_free (m2lnP_int_gal);
  ncm_vector_free (m2lnP_int_gal_dup);
  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
  ncm_serialize_free (ser);
  g_free (dcwl_ser);
  nc_data_cluster_wl_free (dcwl_dup);
}

static void
test_nc_data_cluster_wl_serialize_gauss (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo             = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp     = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist             = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd  = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp           = nc_halo_position_new (dist);
  const guint ngals            = 200;
  NcmVector *m2lnP_int_gal     = ncm_vector_new (ngals);
  NcmVector *m2lnP_int_gal_dup = ncm_vector_new (ngals);
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  GVariant *dcwl_ser;
  NcDataClusterWL *dcwl_dup;
  guint i;

  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (hp), NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (hp), NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, 0.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, 0.0);

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);
  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, hp, ngals, rng, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN);

  ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  dcwl_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->dcwl));
  dcwl_dup = NC_DATA_CLUSTER_WL (ncm_serialize_from_variant (ser, dcwl_ser));

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl_dup));

  nc_data_cluster_wl_eval_m2lnP_integ (test->dcwl, cosmo, dp, smd, hp, m2lnP_int_gal);
  nc_data_cluster_wl_eval_m2lnP_integ (dcwl_dup, cosmo, dp, smd, hp, m2lnP_int_gal_dup);

  for (i = 0; i < ngals; i++)
    g_assert_cmpfloat (ncm_vector_get (m2lnP_int_gal, i), ==, ncm_vector_get (m2lnP_int_gal_dup, i));

  ncm_vector_free (m2lnP_int_gal);
  ncm_vector_free (m2lnP_int_gal_dup);
  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
  ncm_serialize_free (ser);
  g_free (dcwl_ser);
  nc_data_cluster_wl_free (dcwl_dup);
}

static void
test_nc_gata_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata)
{
}

static void
test_nc_data_cluster_wl_m2lnP_integ (TestNcDataClusterWL *test, gconstpointer pdata)
{
}

static void
test_nc_data_cluster_wl_fit (TestNcDataClusterWL *test, gconstpointer pdata)
{
}

/* static void */
/* test_nc_data_cluster_wl_fit (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ()); */
/*   NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ()); */
/*   NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0)); */
/*   NcDistance *dist            = nc_distance_new (100.0); */
/*   NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist); */
/*   const guint ngals           = 200; */
/*   NcmMatrix *gal_obs          = NULL; */

/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0); */

/*   nc_wl_surface_mass_density_prepare (smd, cosmo); */

/*   nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, ngals, rng); */

/*   { */
/*     nc_data_cluster_wl_set_use_kde (test->dcwl, FALSE); */
/*     ncm_data_set_init (NCM_DATA (test->dcwl), TRUE); */

/*     { */
/*       NcmDataset *dataset       = ncm_dataset_new_list (test->dcwl, NULL); */
/*       NcmLikelihood *likelihood = ncm_likelihood_new (dataset); */
/*       NcmMSet *mset             = ncm_mset_new (cosmo, dp, smd, NULL); */
/*       NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", likelihood, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL); */

/*       ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-3, 0.0, NULL, NULL); */
/*       ncm_fit_obs_fisher (fit); */

/*       { */
/*         const gdouble c_bf = ncm_model_param_get (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA); */
/*         const gdouble m_bf = ncm_model_param_get (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA); */
/*         const gdouble c_sd = ncm_fit_covar_sd (fit, nc_halo_density_profile_id (), NC_HALO_DENSITY_PROFILE_C_DELTA); */
/*         const gdouble m_sd = ncm_fit_covar_sd (fit, nc_halo_density_profile_id (), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA); */

/*         g_assert_cmpfloat (c_bf, <, 4.0 + 5.0 * c_sd); */
/*         g_assert_cmpfloat (c_bf, >, 4.0 - 5.0 * c_sd); */

/*         g_assert_cmpfloat (m_bf, <, 14.0 + 5.0 * m_sd); */
/*         g_assert_cmpfloat (m_bf, >, 14.0 - 5.0 * m_sd); */
/*       } */

/*       ncm_dataset_free (dataset); */
/*       ncm_likelihood_free (likelihood); */
/*       ncm_mset_free (mset); */
/*       ncm_fit_free (fit); */
/*     } */
/*   } */

/*   ncm_rng_free (rng); */
/*   nc_hicosmo_free (cosmo); */
/*   nc_halo_density_profile_free (dp); */
/*   nc_distance_free (dist); */
/*   nc_wl_surface_mass_density_free (smd); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_kde_cmp (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ()); */
/*   NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ()); */
/*   NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0)); */
/*   NcDistance *dist            = nc_distance_new (100.0); */
/*   NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist); */
/*   const guint ngals           = 1000; */
/*   NcmStatsVec *m2lnP_stats    = ncm_stats_vec_new (3, NCM_STATS_VEC_COV, FALSE); */
/*   NcmVector *m2lnP_int_gal    = ncm_vector_new (ngals); */
/*   NcmVector *m2lnP_kde_gal    = ncm_vector_new (ngals); */
/*   NcmMatrix *gal_obs          = NULL; */
/*   NcmMatrix *obs; */
/*   guint i; */

/*   ncm_stats_vec_enable_quantile (m2lnP_stats, 0.5); */

/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0); */

/*   nc_wl_surface_mass_density_prepare (smd, cosmo); */

/*   nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, ngals, rng); */
/*   nc_data_cluster_wl_set_ndata (test->dcwl, 1000); */

/*   obs = nc_data_cluster_wl_peek_obs (test->dcwl); */

/*   nc_data_cluster_wl_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_int_gal); */
/*   nc_data_cluster_wl_kde_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_kde_gal); */

/*   for (i = 0; i < ngals; i++) */
/*   { */
/*     const gdouble m2lnP_int = ncm_vector_get (m2lnP_int_gal, i); */
/*     const gdouble m2lnP_kde = ncm_vector_get (m2lnP_kde_gal, i); */

/*     ncm_stats_vec_set (m2lnP_stats, 0, m2lnP_int); */
/*     ncm_stats_vec_set (m2lnP_stats, 1, m2lnP_kde); */
/*     ncm_stats_vec_set (m2lnP_stats, 2, m2lnP_int - m2lnP_kde); */

/*     ncm_stats_vec_update (m2lnP_stats); */
/*   } */

/*   { */
/*     const gdouble mean_diff = ncm_stats_vec_get_mean (m2lnP_stats, 2); */

/*     ncm_stats_vec_reset (m2lnP_stats, TRUE); */

/*     for (i = 0; i < ngals; i++) */
/*     { */
/*       const gdouble m2lnP_int = ncm_vector_get (m2lnP_int_gal, i) - mean_diff; */
/*       const gdouble m2lnP_kde = ncm_vector_get (m2lnP_kde_gal, i); */

/*       ncm_stats_vec_set (m2lnP_stats, 0, m2lnP_int); */
/*       ncm_stats_vec_set (m2lnP_stats, 1, m2lnP_kde); */
/*       ncm_stats_vec_set (m2lnP_stats, 2, fabs (m2lnP_int / m2lnP_kde - 1.0)); */

/*       ncm_stats_vec_update (m2lnP_stats); */
/*     } */
/*   } */

/*   g_assert_cmpfloat (ncm_stats_vec_get_quantile (m2lnP_stats, 2), <, 0.1); */

/*   ncm_vector_free (m2lnP_int_gal); */
/*   ncm_vector_free (m2lnP_kde_gal); */
/*   ncm_rng_free (rng); */
/*   nc_hicosmo_free (cosmo); */
/*   nc_halo_density_profile_free (dp); */
/*   nc_distance_free (dist); */
/*   nc_wl_surface_mass_density_free (smd); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_set_obs (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcmRNG *rng       = ncm_rng_seeded_new (NULL, g_test_rand_int ()); */
/*   const guint ngals = 1000; */
/*   NcmMatrix *obs    = ncm_matrix_new (ngals, 4); */
/*   NcmMatrix *dcwl_obs; */
/*   guint i; */

/*   for (i = 0; i < ngals; i++) */
/*   { */
/*     ncm_matrix_set (obs, i, 0, ncm_rng_uniform_gen_gen (rng, 0.0, 2.0)); */
/*     ncm_matrix_set (obs, i, 1, ncm_rng_uniform_gen_gen (rng, 0.0, 2.0)); */
/*     ncm_matrix_set (obs, i, 2, ncm_rng_uniform_gen_gen (rng, 0.0, 2.0)); */
/*     ncm_matrix_set (obs, i, 3, ncm_rng_uniform_gen_gen (rng, 0.0, 2.0)); */
/*   } */

/*   nc_data_cluster_wl_set_obs (test->dcwl, obs); */

/*   dcwl_obs = nc_data_cluster_wl_peek_obs (test->dcwl); */

/*   for (i = 0; i < ngals; i++) */
/*   { */
/*     const gdouble w = ncm_matrix_get (obs, i, 0); */
/*     const gdouble x = ncm_matrix_get (obs, i, 1); */
/*     const gdouble y = ncm_matrix_get (obs, i, 2); */
/*     const gdouble z = ncm_matrix_get (obs, i, 3); */

/*     g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 0), ==, w); */
/*     g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 1), ==, x); */
/*     g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 2), ==, y); */
/*     g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 3), ==, z); */
/*   } */

/*   ncm_rng_free (rng); */
/*   ncm_matrix_free (obs); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcmRNG *rng                  = ncm_rng_seeded_new (NULL, g_test_rand_int ()); */
/*   NcHICosmo *cosmo             = NC_HICOSMO (nc_hicosmo_de_xcdm_new ()); */
/*   NcHaloDensityProfile *dp     = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0)); */
/*   NcDistance *dist             = nc_distance_new (100.0); */
/*   NcWLSurfaceMassDensity *smd  = nc_wl_surface_mass_density_new (dist); */
/*   const guint ngals            = 200; */
/*   NcmMatrix *gal_obs           = NULL; */
/*   NcmVector *m2lnP_int_gal     = ncm_vector_new (ngals); */
/*   NcmVector *m2lnP_kde_gal     = ncm_vector_new (ngals); */
/*   NcmVector *m2lnP_int_gal_dup = ncm_vector_new (ngals); */
/*   NcmVector *m2lnP_kde_gal_dup = ncm_vector_new (ngals); */
/*   NcmSerialize *ser; */
/*   gchar *sd_ser; */
/*   NcDataClusterWL *sd_dup; */
/*   guint i; */

/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0); */

/*   nc_wl_surface_mass_density_prepare (smd, cosmo); */

/*   nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, ngals, rng); */
/*   nc_data_cluster_wl_set_prec (test->dcwl, 1.0e-3); */
/*   nc_data_cluster_wl_set_use_kde (test->dcwl, FALSE); */

/*   ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE); */
/*   sd_ser = ncm_serialize_to_string (ser, G_OBJECT (test->dcwl), TRUE); */
/*   sd_dup = NC_DATA_CLUSTER_WL (ncm_serialize_from_string (ser, sd_ser)); */

/*   nc_data_cluster_wl_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_int_gal); */
/*   nc_data_cluster_wl_eval_m2lnP (sd_dup, cosmo, dp, smd, m2lnP_int_gal_dup); */

/*   for (i = 0; i < ngals; i++) */
/*   { */
/*     const gdouble m2lnP_int     = ncm_vector_get (m2lnP_int_gal, i); */
/*     const gdouble m2lnP_int_dup = ncm_vector_get (m2lnP_int_gal_dup, i); */

/*     ncm_assert_cmpdouble_e (m2lnP_int, ==, m2lnP_int_dup, 1.0e-3, 0.0); */
/*   } */

/*   nc_data_cluster_wl_kde_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_kde_gal); */
/*   nc_data_cluster_wl_kde_eval_m2lnP (sd_dup, cosmo, dp, smd, m2lnP_kde_gal_dup); */

/*   for (i = 0; i < ngals; i++) */
/*   { */
/*     const gdouble m2lnP_kde     = ncm_vector_get (m2lnP_kde_gal, i); */
/*     const gdouble m2lnP_kde_dup = ncm_vector_get (m2lnP_kde_gal_dup, i); */

/*     ncm_assert_cmpdouble_e (m2lnP_kde, ==, m2lnP_kde_dup, 1.0e-3, 0.0); */
/*   } */

/*   ncm_vector_free (m2lnP_int_gal); */
/*   ncm_vector_free (m2lnP_kde_gal); */
/*   ncm_serialize_free (ser); */
/*   g_free (sd_ser); */
/*   nc_data_cluster_wl_free (sd_dup); */
/*   ncm_rng_free (rng); */
/*   nc_hicosmo_free (cosmo); */
/*   nc_halo_density_profile_free (dp); */
/*   nc_distance_free (dist); */
/*   nc_wl_surface_mass_density_free (smd); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_r_lim (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcDataClusterWL *dcwl = test->dcwl; */
/*   gdouble r_min; */
/*   gdouble r_max; */

/*   g_object_set (dcwl, "r-min", 0.0, "r-max", 2.0, NULL); */
/*   g_object_get (dcwl, "r-min", &r_min, "r-max", &r_max, NULL); */

/*   g_assert_cmpfloat (r_min, ==, 0.0); */
/*   g_assert_cmpfloat (r_max, ==, 2.0); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_peek_kde (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcDataClusterWL *dcwl       = test->dcwl; */
/*   NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ()); */
/*   NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0)); */
/*   NcDistance *dist            = nc_distance_new (100.0); */
/*   NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist); */
/*   NcmStatsDist *kde; */

/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0); */

/*   nc_wl_surface_mass_density_prepare (smd, cosmo); */
/*   nc_data_cluster_wl_set_use_kde (dcwl, TRUE); */
/*   nc_data_cluster_wl_set_ndata (dcwl, 1000); */
/*   nc_data_cluster_wl_prepare_kde (dcwl, cosmo, dp, smd); */

/*   kde = nc_data_cluster_wl_peek_kde (dcwl); */

/*   g_assert_true (NCM_IS_STATS_DIST (kde)); */

/*   nc_hicosmo_free (cosmo); */
/*   nc_halo_density_profile_free (dp); */
/*   nc_distance_free (dist); */
/*   nc_wl_surface_mass_density_free (smd); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_basic (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcDataClusterWL *dcwl = test->dcwl; */
/*   NcDataClusterWL *dcwl2; */

/*   dcwl2 = nc_data_cluster_wl_ref (dcwl); */
/*   nc_data_cluster_wl_clear (&dcwl2); */
/*   g_assert_true (dcwl2 == NULL); */

/*   g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl)); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_prepare_kde (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcDataClusterWL *dcwl       = test->dcwl; */
/*   NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ()); */
/*   NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0)); */
/*   NcDistance *dist            = nc_distance_new (100.0); */
/*   NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist); */
/*   NcmMatrix *obs              = ncm_matrix_new (1, 4); */
/*   gdouble cut_fraction; */

/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0); */

/*   ncm_matrix_set (obs, 0, 0, 0.0); */
/*   ncm_matrix_set (obs, 0, 1, 0.0); */
/*   ncm_matrix_set (obs, 0, 2, 0.0); */
/*   ncm_matrix_set (obs, 0, 3, 0.0); */

/*   nc_data_cluster_wl_set_obs (dcwl, obs); */

/*   nc_wl_surface_mass_density_prepare (smd, cosmo); */
/*   nc_data_cluster_wl_set_use_kde (dcwl, TRUE); */
/*   nc_data_cluster_wl_set_ndata (dcwl, 100); */
/*   nc_data_cluster_wl_set_cut (dcwl, 10.0, 12.0); */
/*   nc_data_cluster_wl_prepare_kde (dcwl, cosmo, dp, smd); */

/*   g_assert_true (gsl_isnan (nc_data_cluster_wl_kde_eval_m2lnP (dcwl, cosmo, dp, smd, NULL))); */

/*   nc_data_cluster_wl_set_cut (dcwl, 0.0, 10.0); */
/*   nc_data_cluster_wl_prepare_kde (dcwl, cosmo, dp, smd); */

/*   g_assert_true (!gsl_isnan (nc_data_cluster_wl_kde_eval_m2lnP (dcwl, cosmo, dp, smd, NULL))); */
/* } */

/* static void */
/* test_nc_data_cluster_wl_val_kde (TestNcDataClusterWL *test, gconstpointer pdata) */
/* { */
/*   NcDataClusterWL *dcwl_int   = test->dcwl; */
/*   NcDataClusterWL *dcwl_kde   = nc_data_cluster_wl_ref (dcwl_int); */
/*   NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ()); */
/*   NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0)); */
/*   NcDistance *dist            = nc_distance_new (100.0); */
/*   NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist); */
/*   NcmMSet *mset               = ncm_mset_new (cosmo, dp, smd, NULL); */
/*   NcmMatrix *obs              = ncm_matrix_new (1, 4); */
/*   gdouble m2lnL_int           = 0.0; */
/*   gdouble m2lnL_kde           = 0.0; */
/*   NcmStatsDist *kde; */
/*   guint i; */

/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0); */
/*   ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0); */

/*   nc_wl_surface_mass_density_prepare (smd, cosmo); */
/*   nc_data_cluster_wl_set_use_kde (dcwl_kde, TRUE); */
/*   nc_data_cluster_wl_set_ndata (dcwl_kde, 1000); */
/*   nc_data_cluster_wl_prepare_kde (dcwl_kde, cosmo, dp, smd); */

/*   ncm_matrix_set (obs, 0, 0, 0.0); */
/*   ncm_matrix_set (obs, 0, 1, 0.0); */
/*   ncm_matrix_set (obs, 0, 2, 0.0); */
/*   ncm_matrix_set (obs, 0, 3, 0.0); */

/*   nc_data_cluster_wl_set_obs (dcwl_kde, obs); */

/*   ncm_data_set_init (NCM_DATA (dcwl_int), TRUE); */
/*   ncm_data_set_init (NCM_DATA (dcwl_kde), TRUE); */

/*   ncm_data_m2lnL_val (NCM_DATA (dcwl_int), mset, &m2lnL_int); */
/*   ncm_data_m2lnL_val (NCM_DATA (dcwl_kde), mset, &m2lnL_kde); */

/*   g_assert_cmpfloat (m2lnL_int, ==, m2lnL_kde); */

/*   nc_data_cluster_wl_free (dcwl_kde); */
/*   nc_hicosmo_free (cosmo); */
/*   nc_halo_density_profile_free (dp); */
/*   nc_distance_free (dist); */
/*   nc_wl_surface_mass_density_free (smd); */
/*   ncm_mset_free (mset); */
/* } */


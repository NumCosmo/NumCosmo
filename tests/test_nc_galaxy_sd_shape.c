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
  NcHICosmo *cosmo;
  NcHaloDensityProfile *density_profile;
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcGalaxySDObsRedshift *galaxy_redshift;
  NcGalaxySDPosition *galaxy_position;
  NcGalaxySDShape *galaxy_shape;
  NcmMSet *mset;
} TestNcGalaxySDShapeGauss;


static void test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_gen (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_required_columns (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_shape/gauss/serialize", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_serialize,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/model_id", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_model_id,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/gen", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_gen,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/required_columns", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_required_columns,
              &test_nc_galaxy_sd_shape_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (0.0, 5.0));
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true_dist));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->cosmo = cosmo;

  test->density_profile      = dp;
  test->halo_position        = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_distance_free (dist);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (s_dist));
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (s_dist));
}

static void
test_nc_galaxy_sd_shape_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  nc_hicosmo_clear (&test->cosmo);
  nc_halo_density_profile_clear (&test->density_profile);
  nc_halo_position_clear (&test->halo_position);
  nc_wl_surface_mass_density_clear (&test->surface_mass_density);

  nc_galaxy_sd_obs_redshift_clear (&test->galaxy_redshift);
  nc_galaxy_sd_position_clear (&test->galaxy_position);

  ncm_mset_free (test->mset);

  NCM_TEST_FREE (nc_galaxy_sd_shape_free, test->galaxy_shape);
}

static void
test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxySDShape *sdsg_dup = NC_GALAXY_SD_SHAPE (ncm_serialize_dup_obj (ser, G_OBJECT (test->galaxy_shape)));

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (sdsg_dup));

  ncm_serialize_free (ser);
  nc_galaxy_sd_shape_free (sdsg_dup);
}

static void
test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmMSet *model_set  = ncm_mset_empty_new ();
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmModel *model_dup = ncm_model_dup (NCM_MODEL (test->galaxy_shape), ser);

  ncm_mset_set (model_set, model_dup, NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (ncm_mset_peek (model_set, nc_galaxy_sd_shape_id ())));

  ncm_model_free (model_dup);
  ncm_mset_free (model_set);
  ncm_serialize_free (ser);
}

static void
test_nc_galaxy_sd_shape_gauss_gen (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble sigma_obs           = 0.03;
  const guint ntest                 = 1000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    gdouble epsilon_1_out, epsilon_2_out, sigma_obs_out;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.2, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -15.0, 15.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -15.0, 15.0);
    nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, sigma_obs, rng);
    nc_galaxy_sd_shape_gauss_data_get (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &sigma_obs_out);

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);
  }
}

static void
test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble sigma_obs           = 0.03;
  const guint ntest                 = 1000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    gdouble epsilon_1_out, epsilon_2_out, sigma_obs_out;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.2, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -15.0, 15.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -15.0, 15.0);
    nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, sigma_obs, rng);

    nc_galaxy_sd_shape_gauss_data_get (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &sigma_obs_out);

    g_assert_cmpfloat (s_data->epsilon_int_1, >=, -1.2);
    g_assert_cmpfloat (s_data->epsilon_int_1, <=, 1.2);
    g_assert_cmpfloat (s_data->epsilon_int_2, >=, -1.2);
    g_assert_cmpfloat (s_data->epsilon_int_2, <=, 1.2);

    /* g_assert_cmpfloat (epsilon_1_out, >=, -1.2); */
    /* g_assert_cmpfloat (epsilon_1_out, <=, 1.2); */
    /* g_assert_cmpfloat (epsilon_2_out, >=, -1.2); */
    /* g_assert_cmpfloat (epsilon_2_out, <=, 1.2); */
    /* g_assert_cmpfloat (sigma_obs_out, ==, sigma_obs); */
  }
}

static void
test_nc_galaxy_sd_shape_gauss_required_columns (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  GList *columns                    = nc_galaxy_sd_shape_data_required_columns (s_data);

  g_assert_cmpuint (g_list_length (columns), ==, 8);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1);
  g_assert_cmpstr (g_list_nth_data (columns, 1), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2);
  g_assert_cmpstr (g_list_nth_data (columns, 2), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1);
  g_assert_cmpstr (g_list_nth_data (columns, 3), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2);
  g_assert_cmpstr (g_list_nth_data (columns, 4), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS);

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_shape_data_free (s_data);
}


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
  NcGalaxySDShape *gsds;
} TestNcGalaxySDShapeGauss;


static void test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_set_models (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_gauss_gen (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_gen_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_integ_optzs (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_prepare (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_get_header (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_gauss_get_vec_size (TestNcGalaxySDShapeGauss *test, gconstpointer pdata);

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

  g_test_add ("/nc/galaxy_sd_shape/gauss/set_models", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_set_models,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/gen", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_gen,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/gen_strong", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_gen_strong,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/integ_optzs", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_integ_optzs,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/prepare", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_prepare,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/get_header", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_get_header,
              &test_nc_galaxy_sd_shape_free);

  g_test_add ("/nc/galaxy_sd_shape/gauss/get_vec_size", TestNcGalaxySDShapeGauss, NULL,
              &test_nc_galaxy_sd_shape_gauss_new,
              &test_nc_galaxy_sd_shape_gauss_get_vec_size,
              &test_nc_galaxy_sd_shape_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_shape_gauss_new (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds = nc_galaxy_sd_shape_gauss_new ();

  test->gsds = NC_GALAXY_SD_SHAPE (gsds);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsds));
}

static void
test_nc_galaxy_sd_shape_free (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_shape_free, test->gsds);
}

static void
test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShape *gsds     = test->gsds;
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *sdsg_ser           = ncm_serialize_to_string (ser, G_OBJECT (gsds), TRUE);
  NcGalaxySDShape *sdsg_dup = NC_GALAXY_SD_SHAPE (ncm_serialize_from_string (ser, sdsg_ser));

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (sdsg_dup));

  ncm_serialize_free (ser);
  g_free (sdsg_ser);
  nc_galaxy_sd_shape_free (sdsg_dup);
}

static void
test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmMSet *model_set = ncm_mset_empty_new ();
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  ncm_mset_set (model_set, ncm_model_dup (NCM_MODEL (test->gsds), ser));

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (ncm_mset_peek (model_set, nc_galaxy_sd_shape_id ())));

  ncm_mset_free (model_set);
  ncm_serialize_free (ser);
}

static void
test_nc_galaxy_sd_shape_set_models (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist   = nc_distance_new (100.0);
  NcHaloPosition *hp = nc_halo_position_new (dist);

  g_assert_true (nc_galaxy_sd_shape_set_models (test->gsds, cosmo, hp));

  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_position_free (hp);
}

static void
test_nc_galaxy_sd_shape_gauss_gen (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds  = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  gdouble z_cluster           = ncm_model_param_get_by_name (NCM_MODEL (hp), "z");
  gdouble ra                  = g_test_rand_double_range (0.0, 1.0e-3);
  gdouble dec                 = g_test_rand_double_range (0.0, 1.0e-3);
  gdouble z                   = g_test_rand_double_range (0.4, 0.6);
  gdouble e1_var              = gsl_pow_2 (ncm_model_param_get_by_name (NCM_MODEL (gsds), "sigma"));
  gdouble e2_var              = gsl_pow_2 (ncm_model_param_get_by_name (NCM_MODEL (gsds), "sigma"));
  guint nruns                 = 10;
  guint ndata                 = 10000;
  gdouble e1_celestial, e2_celestial, e1_euclidean, e2_euclidean;
  gdouble e1_avg, e2_avg, et, theta, phi, r;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);

  r  = nc_halo_position_projected_radius (hp, cosmo, theta);
  et = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cluster, z_cluster);

  e1_avg = -et *cos (2.0 *phi);

  e2_avg = -et *sin (2.0 *phi);

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    guint j;

    for (j = 0; j < ndata; j++)
    {
      gint seed    = g_test_rand_int ();
      NcmRNG *rng1 = ncm_rng_seeded_new (NULL, seed);
      NcmRNG *rng2 = ncm_rng_seeded_new (NULL, seed);

      nc_galaxy_sd_shape_gen (NC_GALAXY_SD_SHAPE (gsds), cosmo, dp, smd, hp, rng1, NC_GALAXY_WL_OBS_COORD_CELESTIAL, ra, dec, z, &e1_celestial, &e2_celestial);
      nc_galaxy_sd_shape_gen (NC_GALAXY_SD_SHAPE (gsds), cosmo, dp, smd, hp, rng2, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, ra, dec, z, &e1_euclidean, &e2_euclidean);

      g_assert_cmpfloat (e1_celestial, ==, e1_euclidean);
      g_assert_cmpfloat (e2_celestial, ==, -e2_euclidean);

      ncm_stats_vec_set (pos_sample, 0, e1_euclidean);
      ncm_stats_vec_set (pos_sample, 1, e2_euclidean);

      ncm_stats_vec_update (pos_sample);
      ncm_rng_free (rng1);
      ncm_rng_free (rng2);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, e1_avg + 5 * sqrt (e1_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, e1_avg - 5 * sqrt (e1_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), <, e2_avg + 5 * sqrt (e2_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), >, e2_avg - 5 * sqrt (e2_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / e1_var - 1.0), <, 0.1);
    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 1) / e2_var - 1.0), <, 0.1);

    ncm_stats_vec_free (pos_sample);
  }

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_gauss_gen_strong (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds  = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_lcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_CRITICAL, 200.0));
  NcDistance *dist            = nc_distance_new (10.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  gdouble z_cluster           = 0.4;
  gdouble z                   = 0.6;
  gdouble ra                  = 0.0035;
  gdouble dec                 = 0.0035;
  guint nruns                 = 10;
  guint ndata                 = 10000;
  gdouble theta, phi, r;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  ncm_model_param_set_by_name (NCM_MODEL (dp), "log10MDelta", log10 (1.0e16));
  ncm_model_param_set_by_name (NCM_MODEL (hp), "z", z_cluster);

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
    guint j;

    for (j = 0; j < ndata; j++)
    {
      gdouble e1, e2;

      nc_galaxy_sd_shape_gen (NC_GALAXY_SD_SHAPE (gsds), cosmo, dp, smd, hp, rng, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, ra, dec, z, &e1, &e2);

      ncm_stats_vec_set (pos_sample, 0, e1 * e1 + e2 * e2);
      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, 1.0);

    ncm_stats_vec_free (pos_sample);
  }

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds  = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  guint nruns                 = 10000;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  for (i = 0; i < nruns; i++)
  {
    NcmVector *data  = ncm_vector_new (6);
    gdouble r        = g_test_rand_double_range (0.0, 1.0);
    gdouble phi      = g_test_rand_double_range (0.0, 2.0 * M_PI);
    gdouble e1       = g_test_rand_double_range (-1.0, 1.0);
    gdouble e1_sigma = g_test_rand_double_range (0.1, 0.5);
    gdouble e2       = g_test_rand_double_range (-1.0, 1.0);
    gdouble e2_sigma = g_test_rand_double_range (0.1, 0.5);
    gdouble z        = g_test_rand_double_range (0.0, 5.0);

    ncm_vector_set (data, 0, r);
    ncm_vector_set (data, 1, phi);
    ncm_vector_set (data, 2, e1);
    ncm_vector_set (data, 3, e1_sigma);
    ncm_vector_set (data, 4, e2);
    ncm_vector_set (data, 5, e2_sigma);

    g_assert_cmpfloat (nc_galaxy_sd_shape_integ (NC_GALAXY_SD_SHAPE (gsds), cosmo, dp, smd, hp, z, data), >=, 0.0);

    ncm_vector_free (data);
  }

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
}

static void
test_nc_galaxy_sd_shape_gauss_integ_optzs (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds  = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  guint nruns                 = 10000;
  guint i;

  nc_wl_surface_mass_density_prepare (smd, cosmo);
  nc_halo_position_prepare (hp, cosmo);

  for (i = 0; i < nruns; i++)
  {
    NcmVector *data  = ncm_vector_new (6);
    gdouble r        = g_test_rand_double_range (0.0, 1.0);
    gdouble phi      = g_test_rand_double_range (0.0, 2.0 * M_PI);
    gdouble e1       = g_test_rand_double_range (-1.0, 1.0);
    gdouble e1_sigma = g_test_rand_double_range (0.1, 0.5);
    gdouble e2       = g_test_rand_double_range (-1.0, 1.0);
    gdouble e2_sigma = g_test_rand_double_range (0.1, 0.5);
    gdouble z        = g_test_rand_double_range (0.0, 5.0);

    ncm_vector_set (data, 0, r);
    ncm_vector_set (data, 1, phi);
    ncm_vector_set (data, 2, e1);
    ncm_vector_set (data, 3, e1_sigma);
    ncm_vector_set (data, 4, e2);
    ncm_vector_set (data, 5, e2_sigma);

    nc_galaxy_sd_shape_integ_optzs_prep (NC_GALAXY_SD_SHAPE (gsds), cosmo, dp, smd, hp, data);

    g_assert_cmpfloat (nc_galaxy_sd_shape_integ_optzs (NC_GALAXY_SD_SHAPE (gsds), cosmo, dp, smd, hp, z, data), >=, 0.0);

    ncm_vector_free (data);
  }

  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
}

static void
test_nc_galaxy_sd_shape_gauss_prepare (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  NcHICosmo *cosmo           = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloPosition *hp         = nc_halo_position_new (nc_distance_new (10.0));
  NcGalaxyWLObsCoord coord   = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  NcmObjArray *data          = ncm_obj_array_new ();
  NcmObjArray *data_prep     = ncm_obj_array_new ();
  NcmObjArray *data_prep1    = ncm_obj_array_new ();
  NcmObjArray *data_prep2    = ncm_obj_array_new ();
  guint ndata                = 1000;
  guint i;

  for (i = 0; i < ndata; i++)
  {
    NcmVector *vec       = ncm_vector_new (6);
    NcmVector *vec_prep  = ncm_vector_new (6);
    NcmVector *vec_prep1 = ncm_vector_new (6);
    NcmVector *vec_prep2 = ncm_vector_new (6);
    gdouble ra           = g_test_rand_double_range (0.0, 1.0e-3);
    gdouble dec          = g_test_rand_double_range (0.0, 1.0e-3);
    gdouble e1           = g_test_rand_double_range (-1.0, 1.0);
    gdouble e1_sigma     = g_test_rand_double_range (0.1, 0.5);
    gdouble e2           = g_test_rand_double_range (-1.0, 1.0);
    gdouble e2_sigma     = g_test_rand_double_range (0.1, 0.5);

    ncm_vector_set (vec, 0, ra);
    ncm_vector_set (vec, 1, dec);
    ncm_vector_set (vec, 2, e1);
    ncm_vector_set (vec, 3, e1_sigma);
    ncm_vector_set (vec, 4, e2);
    ncm_vector_set (vec, 5, e2_sigma);

    ncm_vector_set (vec_prep, 0, 10000.0);
    ncm_vector_set (vec_prep, 1, 10000.0);
    ncm_vector_set (vec_prep, 2, 10000.0);
    ncm_vector_set (vec_prep, 3, 10000.0);
    ncm_vector_set (vec_prep, 4, 10000.0);
    ncm_vector_set (vec_prep, 5, 10000.0);

    ncm_vector_set (vec_prep1, 0, 10000.0);
    ncm_vector_set (vec_prep1, 1, 10000.0);
    ncm_vector_set (vec_prep1, 2, 10000.0);
    ncm_vector_set (vec_prep1, 3, 10000.0);
    ncm_vector_set (vec_prep1, 4, 10000.0);
    ncm_vector_set (vec_prep1, 5, 10000.0);

    ncm_vector_set (vec_prep2, 0, 10000.0);
    ncm_vector_set (vec_prep2, 1, 10000.0);
    ncm_vector_set (vec_prep2, 2, 10000.0);
    ncm_vector_set (vec_prep2, 3, 10000.0);
    ncm_vector_set (vec_prep2, 4, 10000.0);
    ncm_vector_set (vec_prep2, 5, 10000.0);

    ncm_obj_array_add (data, G_OBJECT (vec));
    ncm_obj_array_add (data_prep, G_OBJECT (vec_prep));
    ncm_obj_array_add (data_prep1, G_OBJECT (vec_prep1));
    ncm_obj_array_add (data_prep2, G_OBJECT (vec_prep2));
  }


  nc_halo_position_prepare (hp, cosmo);
  nc_galaxy_sd_shape_set_models (NC_GALAXY_SD_SHAPE (gsds), cosmo, hp);

  g_assert_true (nc_galaxy_sd_shape_prepare (NC_GALAXY_SD_SHAPE (gsds), cosmo, hp, coord, TRUE, data, data_prep));

  for (i = 0; i < ndata; i++)
  {
    NcmVector *vec = NCM_VECTOR (ncm_obj_array_get (data_prep, i));

    g_assert_cmpfloat (ncm_vector_get (vec, 0), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 1), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 2), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 3), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 4), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 5), !=, 10000.0);
  }

  g_assert_false (nc_galaxy_sd_shape_prepare (NC_GALAXY_SD_SHAPE (gsds), cosmo, hp, coord, FALSE, data, data_prep1));

  for (i = 0; i < ndata; i++)
  {
    NcmVector *vec = NCM_VECTOR (ncm_obj_array_get (data_prep1, i));

    g_assert_cmpfloat (ncm_vector_get (vec, 0), ==, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 1), ==, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 2), ==, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 3), ==, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 4), ==, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 5), ==, 10000.0);
  }

  ncm_model_param_set_by_name (NCM_MODEL (hp), "ra", 0.001);

  g_assert_true (nc_galaxy_sd_shape_prepare (NC_GALAXY_SD_SHAPE (gsds), cosmo, hp, coord, FALSE, data, data_prep2));

  for (i = 0; i < ndata; i++)
  {
    NcmVector *vec = NCM_VECTOR (ncm_obj_array_get (data_prep2, i));

    g_assert_cmpfloat (ncm_vector_get (vec, 0), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 1), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 2), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 3), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 4), !=, 10000.0);
    g_assert_cmpfloat (ncm_vector_get (vec, 5), !=, 10000.0);
  }
}

static void
test_nc_galaxy_sd_shape_gauss_get_header (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  GStrv header               = nc_galaxy_sd_shape_get_header (NC_GALAXY_SD_SHAPE (gsds));
  GStrv header_ctrl          = g_strsplit ("ra dec e1 e1_sigma e2 e2_sigma", " ", -1);

  g_assert_cmpstrv (header, header_ctrl);
}

static void
test_nc_galaxy_sd_shape_gauss_get_vec_size (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDShapeGauss *gsds = NC_GALAXY_SD_SHAPE_GAUSS (test->gsds);
  guint vec_size             = nc_galaxy_sd_shape_get_vec_size (NC_GALAXY_SD_SHAPE (gsds));

  g_assert_cmpuint (vec_size, ==, 6);
}


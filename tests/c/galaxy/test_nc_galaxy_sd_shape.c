/***************************************************************************
 *            test_nc_galaxy_sd_shape.c
 *
 *  mon May 07 00:02:42 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiolimadeoliveira@pm.me>
 * Copyright (C) Sandro Dias Pinto Vitenti 2024 <vitenti@uel.br>
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
#include <gsl/gsl_math.h>

typedef struct _TestNcGalaxySDShape
{
  NcHICosmo *cosmo;
  NcHaloMassSummary *hms;
  NcHaloDensityProfile *density_profile;
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcGalaxySDObsRedshift *galaxy_redshift;
  NcGalaxySDPosition *galaxy_position;
  NcGalaxySDShape *galaxy_shape;
  NcmMSet *mset;
  NcGalaxyWLObsEllipConv ell_conv;
  NcGalaxyWLObsCoord ell_coord;
  gdouble g_fiduc;
  glong seed;
  gdouble g_estimator;
} TestNcGalaxySDShape;


static void test_nc_galaxy_sd_shape_hsm_gauss_global_new (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_new (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_free (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_convert_coord (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord_noise (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_convert_coord_noise (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_gen (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_gen (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_integ (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_integ (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_stats (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_stats (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_estimate (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_estimate (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_required_columns (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_required_columns (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_data_setget (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_data_setget (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_strong_lensing (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_strong_lensing (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_at_nodes (TestNcGalaxySDShape *test, gconstpointer pdata);
static void test_nc_galaxy_sd_shape_hsm_gauss_at_nodes (TestNcGalaxySDShape *test, gconstpointer pdata);

static void test_nc_galaxy_sd_shape_hsm_gauss_global_sigma_conversions (void);

typedef struct _TestNcGalaxySDShapeTests
{
  gchar *test_name;

  void (*test_func) (TestNcGalaxySDShape *test, gconstpointer pdata);
} TestNcGalaxySDShapeTests;

typedef struct _TestEllDefinition
{
  NcGalaxyWLObsEllipConv ell_conv;
  NcGalaxyWLObsCoord ell_coord;
  gchar *ell_conv_name;
  gchar *ell_coord_name;
} TestEllDefinition;

gint
main (gint argc, gchar *argv[])
{
  TestNcGalaxySDShapeTests tests_gauss[11] = {
    {"serialize", &test_nc_galaxy_sd_shape_serialize},
    {"model_id", &test_nc_galaxy_sd_shape_model_id},
    {"coord", &test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord},
    {"gen", &test_nc_galaxy_sd_shape_hsm_gauss_global_gen},
    {"integ", &test_nc_galaxy_sd_shape_hsm_gauss_global_integ},
    {"stats", &test_nc_galaxy_sd_shape_hsm_gauss_global_stats},
    {"estimate", &test_nc_galaxy_sd_shape_hsm_gauss_global_estimate},
    {"required_columns", &test_nc_galaxy_sd_shape_hsm_gauss_global_required_columns},
    {"data_setget", &test_nc_galaxy_sd_shape_hsm_gauss_global_data_setget},
    {"strong_lensing", &test_nc_galaxy_sd_shape_hsm_gauss_global_strong_lensing},
    {"at_nodes", &test_nc_galaxy_sd_shape_hsm_gauss_global_at_nodes},
  };
  TestNcGalaxySDShapeTests tests_gauss_hsc[11] = {
    {"serialize", &test_nc_galaxy_sd_shape_serialize},
    {"model_id", &test_nc_galaxy_sd_shape_model_id},
    {"coord", &test_nc_galaxy_sd_shape_hsm_gauss_convert_coord},
    {"gen", &test_nc_galaxy_sd_shape_hsm_gauss_gen},
    {"integ", &test_nc_galaxy_sd_shape_hsm_gauss_integ},
    {"stats", &test_nc_galaxy_sd_shape_hsm_gauss_stats},
    {"estimate", &test_nc_galaxy_sd_shape_hsm_gauss_estimate},
    {"required_columns", &test_nc_galaxy_sd_shape_hsm_gauss_required_columns},
    {"data_setget", &test_nc_galaxy_sd_shape_hsm_gauss_data_setget},
    {"strong_lensing", &test_nc_galaxy_sd_shape_hsm_gauss_strong_lensing},
    {"at_nodes", &test_nc_galaxy_sd_shape_hsm_gauss_at_nodes},
  };

  gint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  TestEllDefinition ell_def_convert[2] = {
    {
      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE,
      NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
      "trace",
      "euclidean"
    },
    {
      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
      NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
      "trace_det",
      "euclidean"
    }
  };
  TestEllDefinition ell_defs[4] = {
    {
      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE,
      NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
      "trace",
      "euclidean"
    },
    {
      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
      NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
      "trace_det",
      "euclidean"
    },
    {
      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE,
      NC_GALAXY_WL_OBS_COORD_CELESTIAL,
      "trace",
      "celestial"
    },
    {
      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
      NC_GALAXY_WL_OBS_COORD_CELESTIAL,
      "trace_det",
      "celestial"
    }
  };

  for (i = 0; i < 2; i++)
  {
    TestEllDefinition *ell_def = &ell_def_convert[i];
    gchar *ell_conv_name       = ell_def->ell_conv_name;

    gchar *test_path_convert_gauss = g_strdup_printf (
      "/nc/galaxy_sd_shape/shapeHSM/gauss/global/%s/convert_coordinate_system", ell_conv_name
    );
    gchar *test_path_convert_gauss_hsc = g_strdup_printf (
      "/nc/galaxy_sd_shape/shapeHSM/gauss/%s/convert_coordinate_system", ell_conv_name
    );
    gchar *test_path_convert_gauss_noise = g_strdup_printf (
      "/nc/galaxy_sd_shape/shapeHSM/gauss/global/%s/convert_coordinate_system_noise", ell_conv_name
    );
    gchar *test_path_convert_gauss_hsc_noise = g_strdup_printf (
      "/nc/galaxy_sd_shape/shapeHSM/gauss/%s/convert_coordinate_system_noise", ell_conv_name
    );

    g_test_add (test_path_convert_gauss, TestNcGalaxySDShape, (gpointer) ell_def,
                &test_nc_galaxy_sd_shape_hsm_gauss_global_new,
                &test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord,
                &test_nc_galaxy_sd_shape_free);

    g_test_add (test_path_convert_gauss_hsc, TestNcGalaxySDShape, (gpointer) ell_def,
                &test_nc_galaxy_sd_shape_hsm_gauss_new,
                &test_nc_galaxy_sd_shape_hsm_gauss_convert_coord,
                &test_nc_galaxy_sd_shape_free);

    g_test_add (test_path_convert_gauss_noise, TestNcGalaxySDShape, (gpointer) ell_def,
                &test_nc_galaxy_sd_shape_hsm_gauss_global_new,
                &test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord_noise,
                &test_nc_galaxy_sd_shape_free);

    g_test_add (test_path_convert_gauss_hsc_noise, TestNcGalaxySDShape, (gpointer) ell_def,
                &test_nc_galaxy_sd_shape_hsm_gauss_new,
                &test_nc_galaxy_sd_shape_hsm_gauss_convert_coord_noise,
                &test_nc_galaxy_sd_shape_free);

    g_free (test_path_convert_gauss);
    g_free (test_path_convert_gauss_hsc);
    g_free (test_path_convert_gauss_noise);
    g_free (test_path_convert_gauss_hsc_noise);
  }

  for (i = 0; i < 4; i++)
  {
    TestEllDefinition *ell_def = &ell_defs[i];
    gchar *ell_conv_name       = ell_def->ell_conv_name;
    gchar *ell_coord_name      = ell_def->ell_coord_name;

    for (j = 0; j < 11; j++)
    {
      gchar *test_name = tests_gauss[j].test_name;

      void (*test_func) (TestNcGalaxySDShape *test, gconstpointer pdata) = tests_gauss[j].test_func;

      gchar *test_path;

      test_path = g_strdup_printf ("/nc/galaxy_sd_shape/shapeHSM/gauss/global/%s/%s/%s", ell_coord_name, ell_conv_name, test_name);
      g_test_add (test_path, TestNcGalaxySDShape, (gpointer) ell_def,
                  &test_nc_galaxy_sd_shape_hsm_gauss_global_new,
                  test_func,
                  &test_nc_galaxy_sd_shape_free);
      g_free (test_path);
    }

    for (j = 0; j < 11; j++)
    {
      gchar *test_name = tests_gauss_hsc[j].test_name;

      void (*test_func) (TestNcGalaxySDShape *test, gconstpointer pdata) = tests_gauss_hsc[j].test_func;

      gchar *test_path;

      test_path = g_strdup_printf ("/nc/galaxy_sd_shape/shapeHSM/gauss/%s/%s/%s", ell_coord_name, ell_conv_name, test_name);
      g_test_add (test_path, TestNcGalaxySDShape, (gpointer) ell_def,
                  &test_nc_galaxy_sd_shape_hsm_gauss_new,
                  test_func,
                  &test_nc_galaxy_sd_shape_free);
      g_free (test_path);
    }
  }

  g_test_add_func ("/nc/galaxy_sd_shape/gauss/sigma_conversion", &test_nc_galaxy_sd_shape_hsm_gauss_global_sigma_conversions);

  g_test_run ();

  return 0;
}

#define TEST_LOG10_MASS 14.0

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_new (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  TestEllDefinition *ell_def          = (TestEllDefinition *) pdata;
  NcGalaxyWLObsEllipConv ell_conv     = ell_def->ell_conv;
  NcGalaxyWLObsCoord ell_coord        = ell_def->ell_coord;
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_hsm_gauss_global_new (ell_conv));
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true_dist, 0.0, 2.0));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));

  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist            = nc_distance_new (100.0);
  NcHaloMassSummary *hms      = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->cosmo = cosmo;

  test->hms                  = hms;
  test->density_profile      = dp;
  test->halo_position        = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->ell_conv  = ell_conv;
  test->ell_coord = ell_coord;

  ncm_model_param_set_by_name (NCM_MODEL (hms), "log10MDelta", TEST_LOG10_MASS, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (s_dist), "sigma", 0.3, NULL);

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_distance_free (dist);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (s_dist));
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (s_dist));
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_new (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  TestEllDefinition *ell_def          = (TestEllDefinition *) pdata;
  NcGalaxyWLObsEllipConv ell_conv     = ell_def->ell_conv;
  NcGalaxyWLObsCoord ell_coord        = ell_def->ell_coord;
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_hsm_gauss_new (ell_conv));
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true_dist, 0.0, 2.0));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));

  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist            = nc_distance_new (100.0);
  NcHaloMassSummary *hms      = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->cosmo = cosmo;

  test->hms                  = hms;
  test->density_profile      = dp;
  test->halo_position        = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->ell_conv  = ell_conv;
  test->ell_coord = ell_coord;

  ncm_model_param_set_by_name (NCM_MODEL (hms), "log10MDelta", TEST_LOG10_MASS, NULL);

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_distance_free (dist);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (s_dist));
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS (s_dist));
}

static void
test_nc_galaxy_sd_shape_free (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  nc_hicosmo_clear (&test->cosmo);
  nc_halo_density_profile_clear (&test->density_profile);
  nc_halo_mass_summary_clear (&test->hms);
  nc_halo_position_clear (&test->halo_position);
  nc_wl_surface_mass_density_clear (&test->surface_mass_density);

  nc_galaxy_sd_obs_redshift_clear (&test->galaxy_redshift);
  nc_galaxy_sd_position_clear (&test->galaxy_position);

  ncm_mset_free (test->mset);

  NCM_TEST_FREE (nc_galaxy_sd_shape_free, test->galaxy_shape);
}

static void
test_nc_galaxy_sd_shape_serialize (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxySDShape *sdsg_dup = NC_GALAXY_SD_SHAPE (ncm_serialize_dup_obj (ser, G_OBJECT (test->galaxy_shape)));

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (sdsg_dup));
  g_assert_true (nc_galaxy_sd_shape_get_ellip_conv (test->galaxy_shape) == test->ell_conv);
  g_assert_true (nc_galaxy_sd_shape_get_ellip_conv (sdsg_dup) == test->ell_conv);

  ncm_serialize_free (ser);
  nc_galaxy_sd_shape_free (sdsg_dup);
}

static void
test_nc_galaxy_sd_shape_model_id (TestNcGalaxySDShape *test, gconstpointer pdata)
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
test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  gint seed                         = g_test_rand_int ();
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmRNG *rng_shape                 = ncm_rng_seeded_new (NULL, seed);
  NcmRNG *rng_shape2                = ncm_rng_seeded_new (NULL, seed);
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble std_noise           = 0.0;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble c1 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m  = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble e_s_1_celestial, e_s_2_celestial;
    gdouble e_s_1_euclidean, e_s_2_euclidean;
    gdouble e_o_1_celestial, e_o_2_celestial, std_noise_celestial, c1_euclidean, c2_euclidean, m_euclidean;
    gdouble e_o_1_euclidean, e_o_2_euclidean, std_noise_euclidean, c1_celestial, c2_celestial, m_celestial;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);

    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                             test->mset, s_data, std_noise, c1, c2, m,
                                             NC_GALAXY_WL_OBS_COORD_CELESTIAL,
                                             rng_shape);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, &e_o_1_celestial,
                                                  &e_o_2_celestial,
                                                  &std_noise_celestial,
                                                  &c1_celestial, &c2_celestial,
                                                  &m_celestial);

    e_s_1_celestial = s_data->epsilon_int_1;
    e_s_2_celestial = s_data->epsilon_int_2;

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = e_s_1_celestial + I * e_s_2_celestial;
      complex double e_o = e_s;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + (c1 + I * c2);

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));  /* LCOV_EXCL_LINE */
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c1 + I * c2;
      }

      ncm_assert_cmpdouble_e (creal (e_o), ==, e_o_1_celestial, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_o), ==, e_o_2_celestial, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (creal (e_s), ==, e_s_1_celestial, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_s), ==, e_s_2_celestial, 1e-10, 1e-10);
    }

    /* The sign of c2 is flipped when converting the coordinate system */
    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                             test->mset, s_data, std_noise, c1, -c2,
                                             m, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
                                             rng_shape2);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, &e_o_1_euclidean,
                                                  &e_o_2_euclidean,
                                                  &std_noise_euclidean,
                                                  &c1_euclidean, &c2_euclidean,
                                                  &m_euclidean);

    e_s_1_euclidean = s_data->epsilon_int_1;
    e_s_2_euclidean = s_data->epsilon_int_2;

    ncm_assert_cmpdouble_e (e_s_1_celestial, ==, e_s_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_s_2_celestial, ==, -e_s_2_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_1_celestial, ==, e_o_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_2_celestial, ==, -e_o_2_euclidean, 1e-10, 1e-10);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = e_s_1_euclidean + I * e_s_2_euclidean;
      complex double e_o = e_s;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      r   = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);
      phi = M_PI - phi;

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + (c1 - I * c2);

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));  /* LCOV_EXCL_LINE */
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c1 - I * c2;
      }

      ncm_assert_cmpdouble_e (creal (e_s), ==, e_s_1_euclidean, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_s), ==, e_s_2_euclidean, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (creal (e_o), ==, e_o_1_euclidean, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_o), ==, e_o_2_euclidean, 1e-10, 1e-10);
    }
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_rng_free (rng);
  ncm_rng_free (rng_shape);
  ncm_rng_free (rng_shape2);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_convert_coord (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  gint seed                         = g_test_rand_int ();
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmRNG *rng_shape                 = ncm_rng_seeded_new (NULL, seed);
  NcmRNG *rng_shape2                = ncm_rng_seeded_new (NULL, seed);
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble std_noise           = 0.0;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
    const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble e_s_1_celestial, e_s_2_celestial;
    gdouble e_s_1_euclidean, e_s_2_euclidean;
    gdouble e_o_1_celestial, e_o_2_celestial, std_shape_celestial;
    gdouble std_noise_celestial, c1_celestial, c2_celestial, m_celestial;
    gdouble e_o_1_euclidean, e_o_2_euclidean, std_shape_euclidean;
    gdouble std_noise_euclidean, c1_euclidean, c2_euclidean, m_euclidean;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                      test->mset, s_data, std_shape, std_noise, c1,
                                      c2, m, NC_GALAXY_WL_OBS_COORD_CELESTIAL,
                                      rng_shape);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, &e_o_1_celestial, &e_o_2_celestial,
                                           &std_shape_celestial, &std_noise_celestial,
                                           &c1_celestial, &c2_celestial, &m_celestial);

    e_s_1_celestial = s_data->epsilon_int_1;
    e_s_2_celestial = s_data->epsilon_int_2;

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = e_s_1_celestial + I * e_s_2_celestial;
      complex double e_o = e_s;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + (c1 + I * c2);

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));  /* LCOV_EXCL_LINE */
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c1 + I * c2;
      }

      ncm_assert_cmpdouble_e (creal (e_s), ==, e_s_1_celestial, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_s), ==, e_s_2_celestial, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (creal (e_o), ==, e_o_1_celestial, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_o), ==, e_o_2_celestial, 1e-10, 1e-10);
    }

    /* The sign of c2 is flipped when converting the coordinate system */
    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                      test->mset, s_data, std_shape, std_noise, c1,
                                      -c2, m, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
                                      rng_shape2);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, &e_o_1_euclidean, &e_o_2_euclidean,
                                           &std_shape_euclidean, &std_noise_euclidean,
                                           &c1_euclidean, &c2_euclidean, &m_euclidean);

    e_s_1_euclidean = s_data->epsilon_int_1;
    e_s_2_euclidean = s_data->epsilon_int_2;

    ncm_assert_cmpdouble_e (e_s_1_celestial, ==, e_s_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_s_2_celestial, ==, -e_s_2_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_1_celestial, ==, e_o_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_2_celestial, ==, -e_o_2_euclidean, 1e-10, 1e-10);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = e_s_1_euclidean + I * e_s_2_euclidean;
      complex double e_o = e_s;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      r   = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);
      phi = M_PI - phi;

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + (c1 - I * c2);

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));  /* LCOV_EXCL_LINE */
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c1 - I * c2;
      }

      ncm_assert_cmpdouble_e (creal (e_s), ==, e_s_1_euclidean, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_s), ==, e_s_2_euclidean, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (creal (e_o), ==, e_o_1_euclidean, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_o), ==, e_o_2_euclidean, 1e-10, 1e-10);
    }
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_rng_free (rng);
  ncm_rng_free (rng_shape);
  ncm_rng_free (rng_shape2);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_convert_coord_noise (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  gint seed                         = g_test_rand_int ();
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmRNG *rng_shape                 = ncm_rng_seeded_new (NULL, seed);
  NcmRNG *rng_shape2                = ncm_rng_seeded_new (NULL, seed);
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble std_noise           = 0.03;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble c1 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m  = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble e_s_1_celestial, e_s_2_celestial;
    gdouble e_s_1_euclidean, e_s_2_euclidean;
    gdouble e_o_1_celestial, e_o_2_celestial, std_noise_celestial, c1_celestial, c2_celestial, m_celestial;
    gdouble e_o_1_euclidean, e_o_2_euclidean, std_noise_euclidean, c1_euclidean, c2_euclidean, m_euclidean;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);

    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                             test->mset, s_data, std_noise, c1, c2, m,
                                             NC_GALAXY_WL_OBS_COORD_CELESTIAL,
                                             rng_shape);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, &e_o_1_celestial,
                                                  &e_o_2_celestial,
                                                  &std_noise_celestial,
                                                  &c1_celestial, &c2_celestial,
                                                  &m_celestial);

    e_s_1_celestial = s_data->epsilon_int_1;
    e_s_2_celestial = s_data->epsilon_int_2;

    /* The sign of c2 is flipped when converting the coordinate system */
    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                             test->mset, s_data, std_noise, c1, -c2,
                                             m, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
                                             rng_shape2);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, &e_o_1_euclidean,
                                                  &e_o_2_euclidean,
                                                  &std_noise_euclidean,
                                                  &c1_euclidean, &c2_euclidean,
                                                  &m_euclidean);

    e_s_1_euclidean = s_data->epsilon_int_1;
    e_s_2_euclidean = s_data->epsilon_int_2;

    ncm_assert_cmpdouble_e (e_s_1_celestial, ==, e_s_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_s_2_celestial, ==, -e_s_2_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_1_celestial, ==, e_o_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_2_celestial, ==, -e_o_2_euclidean, 1e-10, 1e-10);
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_rng_free (rng);
  ncm_rng_free (rng_shape);
  ncm_rng_free (rng_shape2);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_convert_coord_noise (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  gint seed                         = g_test_rand_int ();
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmRNG *rng_shape                 = ncm_rng_seeded_new (NULL, seed);
  NcmRNG *rng_shape2                = ncm_rng_seeded_new (NULL, seed);
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble std_noise           = 0.03;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
    const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble e_s_1_celestial, e_s_2_celestial;
    gdouble e_s_1_euclidean, e_s_2_euclidean;
    gdouble e_o_1_celestial, e_o_2_celestial, std_shape_celestial;
    gdouble std_noise_celestial, c1_celestial, c2_celestial, m_celestial;
    gdouble e_o_1_euclidean, e_o_2_euclidean, std_shape_euclidean;
    gdouble std_noise_euclidean, c1_euclidean, c2_euclidean, m_euclidean;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                      test->mset, s_data, std_shape, std_noise, c1,
                                      c2, m, NC_GALAXY_WL_OBS_COORD_CELESTIAL,
                                      rng_shape);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, &e_o_1_celestial, &e_o_2_celestial,
                                           &std_shape_celestial, &std_noise_celestial,
                                           &c1_celestial, &c2_celestial, &m_celestial);

    e_s_1_celestial = s_data->epsilon_int_1;
    e_s_2_celestial = s_data->epsilon_int_2;

    /* The sign of c2 is flipped when converting the coordinate system */
    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                      test->mset, s_data, std_shape, std_noise, c1,
                                      -c2, m, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
                                      rng_shape2);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, &e_o_1_euclidean, &e_o_2_euclidean,
                                           &std_shape_euclidean, &std_noise_euclidean,
                                           &c1_euclidean, &c2_euclidean, &m_euclidean);

    e_s_1_euclidean = s_data->epsilon_int_1;
    e_s_2_euclidean = s_data->epsilon_int_2;

    ncm_assert_cmpdouble_e (e_s_1_celestial, ==, e_s_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_s_2_celestial, ==, -e_s_2_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_1_celestial, ==, e_o_1_euclidean, 1e-10, 1e-10);
    ncm_assert_cmpdouble_e (e_o_2_celestial, ==, -e_o_2_euclidean, 1e-10, 1e-10);
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_rng_free (rng);
  ncm_rng_free (rng_shape);
  ncm_rng_free (rng_shape2);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_gen (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcmStatsVec *stats                = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
  const gdouble std_noise           = 0.03;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble c1 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m  = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble epsilon_1_out, epsilon_2_out, std_noise_out, c1_out, c2_out, m_out;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                             test->mset, s_data, std_noise, c1, c2, m,
                                             test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, &epsilon_1_out,
                                                  &epsilon_2_out, &std_noise_out,
                                                  &c1_out, &c2_out, &m_out);

    ncm_assert_cmpdouble_e (std_noise_out, ==, std_noise, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c1_out, ==, c1, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c2_out, ==, c2, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (m_out, ==, m, 1.0e-12, 0.0);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = (s_data->epsilon_int_1 + I * s_data->epsilon_int_2);
      complex double e_o = e_s;
      complex double c   = c1 + I * c2;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);

        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + c;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));  /* LCOV_EXCL_LINE */
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c;
      }

      ncm_stats_vec_set (stats, 0, epsilon_1_out - creal (e_o));
      ncm_stats_vec_set (stats, 1, epsilon_2_out - cimag (e_o));

      ncm_stats_vec_update (stats);
    }

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);
  }

  ncm_assert_cmpdouble_e (ncm_stats_vec_get_sd (stats, 0), ==, std_noise, 1.0e-1, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_sd (stats, 1), ==, std_noise, 1.0e-1, 0.0);

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_stats_vec_free (stats);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_gen (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcmStatsVec *stats                = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
  const gdouble std_noise           = 0.03;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
    const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble epsilon_1_out, epsilon_2_out, std_shape_out, std_noise_out, c1_out, c2_out, m_out;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                      test->mset, s_data, std_shape, std_noise, c1,
                                      c2, m, test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, &epsilon_1_out, &epsilon_2_out,
                                           &std_shape_out, &std_noise_out, &c1_out,
                                           &c2_out, &m_out);

    ncm_assert_cmpdouble_e (std_shape_out, ==, std_shape, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (std_noise_out, ==, std_noise, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c1_out, ==, c1, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c2_out, ==, c2, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (m_out, ==, m, 1.0e-12, 0.0);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = (s_data->epsilon_int_1 + I * s_data->epsilon_int_2);
      complex double e_o = e_s;
      complex double c   = c1 + I * c2;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);

        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + c;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));  /* LCOV_EXCL_LINE */
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c;
      }

      ncm_stats_vec_set (stats, 0, epsilon_1_out - creal (e_o));
      ncm_stats_vec_set (stats, 1, epsilon_2_out - cimag (e_o));

      ncm_stats_vec_update (stats);
    }

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);
  }

  ncm_assert_cmpdouble_e (ncm_stats_vec_get_sd (stats, 0), ==, std_noise, 1.0e-1, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_sd (stats, 1), ==, std_noise, 1.0e-1, 0.0);

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_stats_vec_free (stats);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_integ (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data      = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data         = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data            = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcmStatsVec *stats                     = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
  NcGalaxySDShapeIntegrand *integrand    = nc_galaxy_sd_shape_integ (test->galaxy_shape, FALSE);
  NcGalaxySDShapeIntegrand *ln_integrand = nc_galaxy_sd_shape_integ (test->galaxy_shape, TRUE);
  GPtrArray *data_array                  = g_ptr_array_new ();
  const gdouble std_noise                = 0.03;
  const gdouble sigma                    = ncm_model_orig_param_get (NCM_MODEL (test->galaxy_shape), NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SIGMA);
  const guint ntest                      = 10000;
  guint i;

  g_ptr_array_add (data_array, s_data);
  nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);
  nc_galaxy_sd_shape_integrand_prepare (ln_integrand, test->mset);

  for (i = 0; i < ntest; i++)
  {
    const gdouble c1 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m  = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble epsilon_1_out, epsilon_2_out, std_noise_out, c1_out, c2_out, m_out;
    gdouble int0, lnint0;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0e-2, 2.0e-2);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0e-2, 2.0e-2);

    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                             test->mset, s_data, std_noise, c1, c2, m,
                                             test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, &epsilon_1_out,
                                                  &epsilon_2_out, &std_noise_out,
                                                  &c1_out, &c2_out, &m_out);

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);

    ncm_assert_cmpdouble_e (std_noise_out, ==, std_noise, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c1_out, ==, c1, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c2_out, ==, c2, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (m_out, ==, m, 1.0e-12, 0.0);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);

    int0   = nc_galaxy_sd_shape_integrand_eval (integrand, z_data->z, s_data);
    lnint0 = nc_galaxy_sd_shape_integrand_eval (ln_integrand, z_data->z, s_data);

    ncm_assert_cmpdouble_e (log (int0), ==, lnint0, 1e-10, 1.0e-10);

    {
      const gdouble z_cl      = nc_halo_position_get_redshift (test->halo_position);
      complex double e_o      = epsilon_1_out + I * epsilon_2_out;
      complex double data_e_s = s_data->epsilon_int_1 + I * s_data->epsilon_int_2;
      complex double e_s      = 0.0;
      complex double g        = 0.0;
      gdouble gt              = 0.0;
      gdouble theta, phi, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
      }
      else
      {
        const gdouble step = exp ((z_data->z - z_cl) / 0.001);

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_cl * (1.0 + GSL_DBL_EPSILON), z_cl, z_cl);

        gt *= step;
      }

      g = gt * cexp (2.0 * I * phi);
      g = (1.0 + m) * g + (c1 + I * c2);

      switch (test->ell_conv)
      {
        case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

          if (gt > 1.0)
            e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));  /* LCOV_EXCL_LINE */
          else
            e_s = (e_o - g) / (1.0 - conj (g) * e_o);

          break;

        case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
          e_s = (e_o + g * (g * conj (e_o) - 2.0)) / (1.0 + g * conj (g) - 2.0 * creal (g * conj (e_o)));
          break;

        default: /* LCOV_EXCL_LINE */
          g_assert_not_reached ();
      }

      {
        const gdouble var_int   = gsl_pow_2 (sigma);
        const gdouble total_var = var_int + gsl_pow_2 (std_noise_out);
        const gdouble chi2_1    = gsl_pow_2 (creal (e_s)) / total_var;
        const gdouble chi2_2    = gsl_pow_2 (cimag (e_s)) / total_var;
        gdouble jac_num, jac_den, int1, lnint1;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
            jac_num = gsl_pow_2 (1.0 - g * conj (g));
            jac_den = gsl_pow_2 (1.0 - 2.0 * creal (conj (g) * e_o) + g * conj (g) * conj (e_o) * e_o);
            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            jac_num = gsl_pow_3 (1.0 - g * conj (g));
            jac_den = gsl_pow_3 (1.0 - 2.0 * creal (conj (g) * e_o) + g * conj (g));
            break;

          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }

        int1   = exp (-0.5 * (chi2_1 + chi2_2)) / (2.0 * M_PI * total_var) * (jac_num / jac_den);
        lnint1 = -0.5 * (chi2_1 + chi2_2) - log (2.0 * M_PI * total_var) + log (jac_num / jac_den);

        ncm_assert_cmpdouble_e (int0, ==, int1, 1e-10, 1.0e-10);
        ncm_assert_cmpdouble_e (lnint0, ==, lnint1, 1e-10, 1.0e-10);

        e_s      = e_s * cexp (-2.0 * I * phi);
        e_o      = e_o * cexp (-2.0 * I * phi);
        data_e_s = data_e_s * cexp (-2.0 * I * phi);

        /* For the trace convention e_o is a biased estimator for g. */
        if (test->ell_conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE)
          e_o = 0.5 * e_o / (1.0 - var_int);

        ncm_stats_vec_set (stats, 0, creal (e_s));
        ncm_stats_vec_set (stats, 1, cimag (e_s));
        ncm_stats_vec_set (stats, 2, creal (data_e_s));
        ncm_stats_vec_set (stats, 3, cimag (data_e_s));
        ncm_stats_vec_set (stats, 4, creal (e_o));
        ncm_stats_vec_set (stats, 5, cimag (e_o));
        ncm_stats_vec_set (stats, 6, gt);

        ncm_stats_vec_update (stats);
      }
    }
  }

  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 0), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 1), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 2), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 3), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 4), ==, ncm_stats_vec_get_mean (stats, 6), 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 5), ==, 0.0, 0.0, 2.0e-2);

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_sd_shape_integrand_free (integrand);
  ncm_stats_vec_free (stats);
  ncm_rng_free (rng);
  g_ptr_array_unref (data_array);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_integ (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data      = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data         = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data            = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcmStatsVec *stats                     = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
  NcGalaxySDShapeIntegrand *integrand    = nc_galaxy_sd_shape_integ (test->galaxy_shape, FALSE);
  NcGalaxySDShapeIntegrand *ln_integrand = nc_galaxy_sd_shape_integ (test->galaxy_shape, TRUE);
  GPtrArray *data_array                  = g_ptr_array_new ();
  const gdouble std_noise                = 0.03;
  const guint ntest                      = 10000;
  guint i;

  g_ptr_array_add (data_array, s_data);
  nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);
  nc_galaxy_sd_shape_integrand_prepare (ln_integrand, test->mset);

  for (i = 0; i < ntest; i++)
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
    const gdouble sigma     = nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (std_shape);
    const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble epsilon_1_out, epsilon_2_out, std_noise_out, std_shape_out, c1_out, c2_out, m_out;
    gdouble int0, lnint0;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.1, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0e-2, 2.0e-2);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0e-2, 2.0e-2);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                      test->mset, s_data, std_shape, std_noise, c1, c2,
                                      m, test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, &epsilon_1_out, &epsilon_2_out,
                                           &std_shape_out, &std_noise_out, &c1_out,
                                           &c2_out, &m_out);

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);

    ncm_assert_cmpdouble_e (std_shape_out, ==, std_shape, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (std_noise_out, ==, std_noise, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c1_out, ==, c1, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (c2_out, ==, c2, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (m_out, ==, m, 1.0e-12, 0.0);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);

    int0   = nc_galaxy_sd_shape_integrand_eval (integrand, z_data->z, s_data);
    lnint0 = nc_galaxy_sd_shape_integrand_eval (ln_integrand, z_data->z, s_data);

    ncm_assert_cmpdouble_e (log (int0), ==, lnint0, 1e-10, 1.0e-10);

    {
      const gdouble z_cl      = nc_halo_position_get_redshift (test->halo_position);
      complex double e_o      = epsilon_1_out + I * epsilon_2_out;
      complex double data_e_s = s_data->epsilon_int_1 + I * s_data->epsilon_int_2;
      complex double e_s      = 0.0;
      complex double g        = 0.0;
      gdouble gt              = 0.0;
      gdouble theta, phi, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
      }
      else
      {
        const gdouble step = exp ((z_data->z - z_cl) / 0.001);

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_cl * (1.0 + GSL_DBL_EPSILON), z_cl, z_cl);

        gt *= step;
      }

      g = gt * cexp (2.0 * I * phi);
      g = (1.0 + m) * g + (c1 + I * c2);

      switch (test->ell_conv)
      {
        case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

          if (gt > 1.0)
            e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
          else
            e_s = (e_o - g) / (1.0 - conj (g) * e_o);  /* LCOV_EXCL_LINE */

          break;

        case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
          e_s = (e_o + g * (g * conj (e_o) - 2.0)) / (1.0 + g * conj (g) - 2.0 * creal (g * conj (e_o)));
          break;

        default: /* LCOV_EXCL_LINE */
          g_assert_not_reached ();
      }

      {
        const gdouble var_int   = gsl_pow_2 (sigma);
        const gdouble total_var = var_int + gsl_pow_2 (std_noise_out);
        const gdouble chi2_1    = gsl_pow_2 (creal (e_s)) / total_var;
        const gdouble chi2_2    = gsl_pow_2 (cimag (e_s)) / total_var;
        gdouble jac_num, jac_den, int1, lnint1;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
            jac_num = gsl_pow_2 (1.0 - g * conj (g));
            jac_den = gsl_pow_2 (1.0 - 2.0 * creal (conj (g) * e_o) + g * conj (g) * conj (e_o) * e_o);
            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            jac_num = gsl_pow_3 (1.0 - g * conj (g));
            jac_den = gsl_pow_3 (1.0 - 2.0 * creal (conj (g) * e_o) + g * conj (g));
            break;

          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }

        int1   = exp (-0.5 * (chi2_1 + chi2_2)) / (2.0 * M_PI * total_var) * (jac_num / jac_den);
        lnint1 = -0.5 * (chi2_1 + chi2_2) - log (2.0 * M_PI * total_var) + log (jac_num / jac_den);

        ncm_assert_cmpdouble_e (int0, ==, int1, 1e-10, 1.0e-10);
        ncm_assert_cmpdouble_e (lnint0, ==, lnint1, 1e-10, 1.0e-10);

        e_s      = e_s * cexp (-2.0 * I * phi);
        e_o      = e_o * cexp (-2.0 * I * phi);
        data_e_s = data_e_s * cexp (-2.0 * I * phi);

        /* For the trace convention e_o is a biased estimator for g. */
        if (test->ell_conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE)
          e_o = 0.5 * e_o / (1.0 - var_int);

        ncm_stats_vec_set (stats, 0, creal (e_s));
        ncm_stats_vec_set (stats, 1, cimag (e_s));
        ncm_stats_vec_set (stats, 2, creal (data_e_s));
        ncm_stats_vec_set (stats, 3, cimag (data_e_s));
        ncm_stats_vec_set (stats, 4, creal (e_o));
        ncm_stats_vec_set (stats, 5, cimag (e_o));
        ncm_stats_vec_set (stats, 6, gt);

        ncm_stats_vec_update (stats);
      }
    }
  }

  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 0), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 1), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 2), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 3), ==, 0.0, 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 4), ==, ncm_stats_vec_get_mean (stats, 6), 0.0, 2.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 5), ==, 0.0, 0.0, 2.0e-2);

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_sd_shape_integrand_free (integrand);
  ncm_stats_vec_free (stats);
  ncm_rng_free (rng);
  g_ptr_array_unref (data_array);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_stats (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble std_noise           = 0.0;
  const guint ngals                 = 2000;
  GList *required_columns           = nc_galaxy_sd_shape_data_required_columns (s_data);
  GList *required_columns_iter      = required_columns;
  GStrvBuilder *builder             = g_strv_builder_new ();
  NcmStatsVec *stats                = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
  GStrv required_columns_strv;
  NcGalaxyWLObs *obs;
  gdouble radius;

  guint i;

  while (required_columns_iter)
  {
    g_strv_builder_add (builder, required_columns_iter->data);
    required_columns_iter = g_list_next (required_columns_iter);
  }

  required_columns_strv = g_strv_builder_end (builder);
  g_strv_builder_unref (builder);

  obs = nc_galaxy_wl_obs_new (
    nc_galaxy_sd_shape_get_ellip_conv (test->galaxy_shape),
    test->ell_coord,
    ngals,
    required_columns_strv);

  z_data->z   = 0.80;
  p_data->ra  = 0.01;
  p_data->dec = -0.03;

  for (i = 0; i < ngals; i++)
  {
    const gdouble c1 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m  = ncm_rng_uniform_gen (rng, -0.2, 0.2);

    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), test->mset, s_data, std_noise, c1, c2, m, test->ell_coord, rng);
    nc_galaxy_sd_shape_data_write_row (s_data, obs, i);
  }

  radius = nc_galaxy_sd_shape_data_get_radius (s_data);

  {
    NcDataClusterWL *dcwl = nc_data_cluster_wl_new ();
    NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    const gdouble z_cl    = nc_halo_position_get_redshift (test->halo_position);
    const gdouble true_gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                                      test->density_profile,
                                                                      test->cosmo,
                                                                      radius, z_data->z, z_cl, z_cl);

    nc_data_cluster_wl_set_obs (dcwl, obs);
    nc_data_cluster_wl_set_resample_flag (dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_SHAPE);
    ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);

    {
      NcmObjArray *data_array = nc_data_cluster_wl_peek_data_array (dcwl);
      NcmData *data           = NCM_DATA (dcwl);
      NcmDataset *dataset     = ncm_dataset_new_array (&data, 1);
      NcmLikelihood *like     = ncm_likelihood_new (dataset);
      NcmFit *fit             = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", like, test->mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
      NcmStatsVec *stats      = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
      guint nfits             = 20;

      for (i = 0; i < nfits; i++)
      {
        gdouble gt, hat_gt, hat_sigma_gt, hat_gx, hat_sigma_gx, hat_rho;

        ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", TEST_LOG10_MASS, NULL);
        ncm_data_resample (data, test->mset, rng);
        ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

        nc_galaxy_sd_shape_direct_estimate (test->galaxy_shape, test->mset, data_array, &hat_gt, &hat_gx, &hat_sigma_gt, &hat_sigma_gx, &hat_rho);

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       radius, z_data->z, z_cl, z_cl);

        ncm_stats_vec_set (stats, 0, ncm_mset_fparam_get (test->mset, 0));
        ncm_stats_vec_set (stats, 1, gt);
        ncm_stats_vec_set (stats, 2, hat_gt);
        ncm_stats_vec_set (stats, 3, hat_gx);
        ncm_stats_vec_set (stats, 4, hat_sigma_gt);
        ncm_stats_vec_set (stats, 5, hat_sigma_gx);
        ncm_stats_vec_set (stats, 6, hat_rho);

        ncm_stats_vec_update (stats);
      }

      {
        const gdouble estimated_log10M = ncm_stats_vec_get_mean (stats, 0);
        const gdouble sigma_log10M     = ncm_stats_vec_get_sd (stats, 0);

        g_assert (gsl_finite (estimated_log10M));
        g_assert (gsl_finite (sigma_log10M));

        ncm_assert_cmpdouble (estimated_log10M, >, TEST_LOG10_MASS - 5.0 * sigma_log10M / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_log10M, <, TEST_LOG10_MASS + 5.0 * sigma_log10M / sqrt (nfits));
      }

      {
        const gdouble estimated_gt = ncm_stats_vec_get_mean (stats, 1);
        const gdouble sigma_gt     = ncm_stats_vec_get_sd (stats, 1);

        g_assert (gsl_finite (estimated_gt));
        g_assert (gsl_finite (sigma_gt));

        ncm_assert_cmpdouble (estimated_gt, >, true_gt - 5.0 * sigma_gt / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_gt, <, true_gt + 5.0 * sigma_gt / sqrt (nfits));
      }

      {
        const gdouble estimated_hat_gt = ncm_stats_vec_get_mean (stats, 2);
        const gdouble sigma_hat_gt     = ncm_stats_vec_get_sd (stats, 2);

        g_assert (gsl_finite (estimated_hat_gt));
        g_assert (gsl_finite (sigma_hat_gt));

        ncm_assert_cmpdouble (estimated_hat_gt, >, true_gt - 5.0 * sigma_hat_gt / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_hat_gt, <, true_gt + 5.0 * sigma_hat_gt / sqrt (nfits));
      }

      {
        const gdouble estimated_hat_gx = ncm_stats_vec_get_mean (stats, 3);
        const gdouble sigma_hat_gx     = ncm_stats_vec_get_sd (stats, 3);

        g_assert (gsl_finite (estimated_hat_gx));
        g_assert (gsl_finite (sigma_hat_gx));

        ncm_assert_cmpdouble (estimated_hat_gx, >=, 0.0 - 5.0 * sigma_hat_gx / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_hat_gx, <=, 0.0 + 5.0 * sigma_hat_gx / sqrt (nfits));
      }

      ncm_data_free (data);
      ncm_likelihood_free (like);
      ncm_dataset_free (dataset);
      ncm_fit_free (fit);
    }

    ncm_serialize_free (ser);
  }

  g_list_free_full (required_columns, g_free);
  g_strfreev (required_columns_strv);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_rng_free (rng);
  nc_galaxy_wl_obs_free (obs);
  ncm_stats_vec_free (stats);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_stats (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble std_noise           = 0.0;
  const guint ngals                 = 2000;
  GList *required_columns           = nc_galaxy_sd_shape_data_required_columns (s_data);
  GList *required_columns_iter      = required_columns;
  GStrvBuilder *builder             = g_strv_builder_new ();
  NcmStatsVec *stats                = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
  GStrv required_columns_strv;
  NcGalaxyWLObs *obs;
  gdouble radius;

  guint i;

  while (required_columns_iter)
  {
    g_strv_builder_add (builder, required_columns_iter->data);
    required_columns_iter = g_list_next (required_columns_iter);
  }

  required_columns_strv = g_strv_builder_end (builder);
  g_strv_builder_unref (builder);

  obs = nc_galaxy_wl_obs_new (
    nc_galaxy_sd_shape_get_ellip_conv (test->galaxy_shape),
    test->ell_coord,
    ngals,
    required_columns_strv);

  z_data->z   = 0.80;
  p_data->ra  = 0.01;
  p_data->dec = -0.03;

  for (i = 0; i < ngals; i++)
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
    const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), test->mset, s_data,
                                      std_shape, std_noise, c1, c2, m,
                                      test->ell_coord, rng);
    nc_galaxy_sd_shape_data_write_row (s_data, obs, i);
  }

  radius = nc_galaxy_sd_shape_data_get_radius (s_data);

  {
    NcDataClusterWL *dcwl = nc_data_cluster_wl_new ();
    NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    const gdouble z_cl    = nc_halo_position_get_redshift (test->halo_position);
    const gdouble true_gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                                      test->density_profile,
                                                                      test->cosmo,
                                                                      radius, z_data->z, z_cl, z_cl);

    nc_data_cluster_wl_set_obs (dcwl, obs);
    nc_data_cluster_wl_set_resample_flag (dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_SHAPE);
    ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);

    {
      NcmObjArray *data_array = nc_data_cluster_wl_peek_data_array (dcwl);
      NcmData *data           = NCM_DATA (dcwl);
      NcmDataset *dataset     = ncm_dataset_new_array (&data, 1);
      NcmLikelihood *like     = ncm_likelihood_new (dataset);
      NcmFit *fit             = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", like, test->mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
      NcmStatsVec *stats      = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
      guint nfits             = 50;

      for (i = 0; i < nfits; i++)
      {
        gdouble gt, hat_gt, hat_sigma_gt, hat_gx, hat_sigma_gx, hat_rho;

        ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", TEST_LOG10_MASS, NULL);
        ncm_data_resample (data, test->mset, rng);
        ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

        nc_galaxy_sd_shape_direct_estimate (test->galaxy_shape, test->mset, data_array, &hat_gt, &hat_gx, &hat_sigma_gt, &hat_sigma_gx, &hat_rho);

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       radius, z_data->z, z_cl, z_cl);

        ncm_stats_vec_set (stats, 0, ncm_mset_fparam_get (test->mset, 0));
        ncm_stats_vec_set (stats, 1, gt);
        ncm_stats_vec_set (stats, 2, hat_gt);
        ncm_stats_vec_set (stats, 3, hat_gx);
        ncm_stats_vec_set (stats, 4, hat_sigma_gt);
        ncm_stats_vec_set (stats, 5, hat_sigma_gx);
        ncm_stats_vec_set (stats, 6, hat_rho);

        ncm_stats_vec_update (stats);
      }

      {
        const gdouble estimated_log10M = ncm_stats_vec_get_mean (stats, 0);
        const gdouble sigma_log10M     = ncm_stats_vec_get_sd (stats, 0);

        g_assert (gsl_finite (estimated_log10M));
        g_assert (gsl_finite (sigma_log10M));

        ncm_assert_cmpdouble (estimated_log10M, >, TEST_LOG10_MASS - 5.0 * sigma_log10M / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_log10M, <, TEST_LOG10_MASS + 5.0 * sigma_log10M / sqrt (nfits));
      }

      {
        const gdouble estimated_gt = ncm_stats_vec_get_mean (stats, 1);
        const gdouble sigma_gt     = ncm_stats_vec_get_sd (stats, 1);

        g_assert (gsl_finite (estimated_gt));
        g_assert (gsl_finite (sigma_gt));

        ncm_assert_cmpdouble (estimated_gt, >, true_gt - 5.0 * sigma_gt / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_gt, <, true_gt + 5.0 * sigma_gt / sqrt (nfits));
      }

      {
        const gdouble estimated_hat_gt = ncm_stats_vec_get_mean (stats, 2);
        const gdouble sigma_hat_gt     = ncm_stats_vec_get_sd (stats, 2);

        g_assert (gsl_finite (estimated_hat_gt));
        g_assert (gsl_finite (sigma_hat_gt));

        ncm_assert_cmpdouble (estimated_hat_gt, >, true_gt - 5.0 * sigma_hat_gt / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_hat_gt, <, true_gt + 5.0 * sigma_hat_gt / sqrt (nfits));
      }

      {
        const gdouble estimated_hat_gx = ncm_stats_vec_get_mean (stats, 3);
        const gdouble sigma_hat_gx     = ncm_stats_vec_get_sd (stats, 3);

        g_assert (gsl_finite (estimated_hat_gx));
        g_assert (gsl_finite (sigma_hat_gx));

        ncm_assert_cmpdouble (estimated_hat_gx, >=, 0.0 - 5.0 * sigma_hat_gx / sqrt (nfits));
        ncm_assert_cmpdouble (estimated_hat_gx, <=, 0.0 + 5.0 * sigma_hat_gx / sqrt (nfits));
      }

      ncm_data_free (data);
      ncm_likelihood_free (like);
      ncm_dataset_free (dataset);
      ncm_fit_free (fit);
      ncm_stats_vec_free (stats);
    }

    ncm_serialize_free (ser);
  }

  g_list_free_full (required_columns, g_free);
  g_strfreev (required_columns_strv);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  ncm_rng_free (rng);
  nc_galaxy_wl_obs_free (obs);
  ncm_stats_vec_free (stats);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_required_columns (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  GList *columns                    = nc_galaxy_sd_shape_data_required_columns (s_data);

  g_assert_cmpuint (g_list_length (columns), ==, 11);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1);
  g_assert_cmpstr (g_list_nth_data (columns, 1), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2);
  g_assert_cmpstr (g_list_nth_data (columns, 2), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_EPSILON_OBS_1);
  g_assert_cmpstr (g_list_nth_data (columns, 3), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_EPSILON_OBS_2);
  g_assert_cmpstr (g_list_nth_data (columns, 4), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_STD_NOISE);
  g_assert_cmpstr (g_list_nth_data (columns, 5), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C1);
  g_assert_cmpstr (g_list_nth_data (columns, 6), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C2);
  g_assert_cmpstr (g_list_nth_data (columns, 7), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_M);

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_required_columns (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  GList *columns                    = nc_galaxy_sd_shape_data_required_columns (s_data);

  g_assert_cmpuint (g_list_length (columns), ==, 12);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1);
  g_assert_cmpstr (g_list_nth_data (columns, 1), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2);
  g_assert_cmpstr (g_list_nth_data (columns, 2), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_1);
  g_assert_cmpstr (g_list_nth_data (columns, 3), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_2);
  g_assert_cmpstr (g_list_nth_data (columns, 4), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_SHAPE);
  g_assert_cmpstr (g_list_nth_data (columns, 5), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_NOISE);
  g_assert_cmpstr (g_list_nth_data (columns, 6), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C1);
  g_assert_cmpstr (g_list_nth_data (columns, 7), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C2);
  g_assert_cmpstr (g_list_nth_data (columns, 8), ==, NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_M);

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_data_setget (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble epsilon_1           = 0.1;
  const gdouble epsilon_2           = 0.2;
  const gdouble std_noise           = 0.03;
  const gdouble c1                  = 0.05;
  const gdouble c2                  = 0.06;
  const gdouble m                   = 0.07;

  nc_galaxy_sd_shape_hsm_gauss_global_data_set (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), s_data, epsilon_1, epsilon_2, std_noise, c1, c2, m);

  {
    gdouble epsilon_1_out, epsilon_2_out, std_noise_out, c1_out, c2_out, m_out;

    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &std_noise_out, &c1_out, &c2_out, &m_out);

    g_assert_cmpfloat (epsilon_1_out, ==, epsilon_1);
    g_assert_cmpfloat (epsilon_2_out, ==, epsilon_2);
    g_assert_cmpfloat (std_noise_out, ==, std_noise);
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_data_setget (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  const gdouble epsilon_1           = 0.1;
  const gdouble epsilon_2           = 0.2;
  const gdouble std_noise           = 0.03;
  const gdouble std_shape           = 0.04;
  const gdouble c1                  = 0.05;
  const gdouble c2                  = 0.06;
  const gdouble m                   = 0.07;

  nc_galaxy_sd_shape_hsm_gauss_data_set (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), s_data,
                                         epsilon_1, epsilon_2, std_shape, std_noise, c1, c2, m);

  {
    gdouble epsilon_1_out, epsilon_2_out, std_noise_out, std_shape_out, c1_out, c2_out, m_out;

    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), s_data,
                                           &epsilon_1_out, &epsilon_2_out, &std_shape_out, &std_noise_out, &c1_out, &c2_out, &m_out);

    g_assert_cmpfloat (epsilon_1_out, ==, epsilon_1);
    g_assert_cmpfloat (epsilon_2_out, ==, epsilon_2);
    g_assert_cmpfloat (std_noise_out, ==, std_noise);
    g_assert_cmpfloat (std_shape_out, ==, std_shape);
    g_assert_cmpfloat (c1_out, ==, c1);
    g_assert_cmpfloat (c2_out, ==, c2);
    g_assert_cmpfloat (m_out, ==, m);
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_strong_lensing (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data   = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data      = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data         = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcGalaxySDShapeIntegrand *integrand = nc_galaxy_sd_shape_integ (test->galaxy_shape, TRUE);
  GPtrArray *data_array               = g_ptr_array_new ();
  const gdouble std_noise             = 0.0;
  const gdouble sigma                 = ncm_model_orig_param_get (NCM_MODEL (test->galaxy_shape), NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SIGMA);
  const guint ntest                   = 10000;
  guint i;

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 15, NULL);

  g_ptr_array_add (data_array, s_data);
  nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);

  for (i = 0; i < ntest; i++)
  {
    const gdouble c1 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2 = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m  = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble epsilon_1_out, epsilon_2_out, std_noise_out, c1_out, c2_out, m_out, int0;

    z_data->z   = 4.2 + ncm_rng_uniform_gen (rng, -0.1, 0.1);
    p_data->ra  = 0.0001 + ncm_rng_uniform_gen (rng, -0.00001, 0.00001);
    p_data->dec = 0.0001 + ncm_rng_uniform_gen (rng, -0.00001, 0.00001);

    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), test->mset, s_data, std_noise, c1, c2, m, test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &std_noise_out, &c1_out, &c2_out, &m_out);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = (s_data->epsilon_int_1 + I * s_data->epsilon_int_2);
      complex double e_o = e_s;
      complex double c   = c1 + I * c2;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);

        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + c;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);  /* LCOV_EXCL_LINE */

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c; /* LCOV_EXCL_LINE */
      }

      ncm_assert_cmpdouble_e (creal (e_o), ==, epsilon_1_out, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_o), ==, epsilon_2_out, 1e-10, 1e-10);
    }

    nc_galaxy_sd_shape_hsm_gauss_global_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), test->mset, s_data, 0.03, c1, c2, m, test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_global_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &std_noise_out, &c1_out, &c2_out, &m_out);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);
    int0 = nc_galaxy_sd_shape_integrand_eval (integrand, z_data->z, s_data);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_o = epsilon_1_out + I * epsilon_2_out;
      complex double e_s = 0.0;
      complex double g   = 0.0;
      gdouble gt         = 0.0;
      gdouble theta, phi, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + (c1 + I * c2);

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (gt > 1.0)
              e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
            else
              e_s = (e_o - g) / (1.0 - conj (g) * e_o);  /* LCOV_EXCL_LINE */

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_s = (e_o + g * (g * conj (e_o) - 2.0)) / (1.0 + g * conj (g) - 2.0 * creal (g * conj (e_o)));
            break;

          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
        }
      }
      else
      {
        e_s = e_o - (c1 + I * c2); /* LCOV_EXCL_LINE */
      }

      {
        const gdouble var_int       = gsl_pow_2 (sigma);
        const gdouble total_var     = var_int + gsl_pow_2 (std_noise_out);
        const gdouble chi2_1        = gsl_pow_2 (creal (e_s)) / total_var;
        const gdouble chi2_2        = gsl_pow_2 (cimag (e_s)) / total_var;
        const complex double g_conj = conj (g);
        const gdouble abs_g2        = g * g_conj;
        gdouble jac_num, jac_den, m2ln_int1;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
            jac_num = gsl_pow_2 (1.0 - abs_g2);

            if (abs_g2 <= 1.0)
              jac_den = gsl_pow_2 (1.0 - 2.0 * creal (conj (g) * e_o) + abs_g2 * conj (e_o) * e_o);  /* LCOV_EXCL_LINE */
            else
              jac_den = gsl_pow_2 (conj (e_o - g) * (e_o - g));

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            jac_num = fabs (gsl_pow_3 (1.0 - abs_g2));
            jac_den = fabs (gsl_pow_3 (1.0 - 2.0 * creal (conj (g) * e_o) + abs_g2));
            break;

          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }

        m2ln_int1 = chi2_1 + chi2_2 + 2.0 * log (2.0 * M_PI * total_var) + 2.0 * log (jac_den / jac_num);

        ncm_assert_cmpdouble_e (-2.0 * int0, ==, m2ln_int1, 1.0e-10, 1.0e-10);
      }
    }
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_sd_shape_integrand_free (integrand);
  g_ptr_array_unref (data_array);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_strong_lensing (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data   = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data      = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data         = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcGalaxySDShapeIntegrand *integrand = nc_galaxy_sd_shape_integ (test->galaxy_shape, TRUE);
  GPtrArray *data_array               = g_ptr_array_new ();
  const gdouble std_noise             = 0.0;
  const guint ntest                   = 10000;
  guint i;

  ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", 15, NULL);

  g_ptr_array_add (data_array, s_data);
  nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);

  for (i = 0; i < ntest; i++)
  {
    const gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
    const gdouble sigma     = nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (std_shape);
    const gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
    const gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    gdouble epsilon_1_out, epsilon_2_out, std_shape_out, std_noise_out, c1_out, c2_out, m_out, int0;

    z_data->z   = 4.2 + ncm_rng_uniform_gen (rng, -0.1, 0.1);
    p_data->ra  = 0.0001 + ncm_rng_uniform_gen (rng, -0.00001, 0.00001);
    p_data->dec = 0.0001 + ncm_rng_uniform_gen (rng, -0.00001, 0.00001);

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), test->mset, s_data,
                                      std_shape, std_noise, c1, c2, m,
                                      test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), s_data,
                                           &epsilon_1_out, &epsilon_2_out, &std_shape_out, &std_noise_out, &c1_out, &c2_out, &m_out);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = (s_data->epsilon_int_1 + I * s_data->epsilon_int_2);
      complex double e_o = e_s;
      complex double c   = c1 + I * c2;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        complex double g;

        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);

        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + c;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (fabs (gt) > 1.0)
              e_o = (1.0 + g * conj (e_s)) / (conj (e_s) + conj (g));
            else
              e_o = (e_s + g) / (1.0 + conj (g) * e_s);  /* LCOV_EXCL_LINE */

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_o = (e_s + g * (g * conj (e_s) + 2.0)) / (1.0 + g * conj (g) + 2.0 * creal (g * conj (e_s)));
            break;
          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
            break;
        }
      }
      else
      {
        e_o = e_s + c; /* LCOV_EXCL_LINE */
      }

      ncm_assert_cmpdouble_e (creal (e_o), ==, epsilon_1_out, 1e-10, 1e-10);
      ncm_assert_cmpdouble_e (cimag (e_o), ==, epsilon_2_out, 1e-10, 1e-10);
    }

    nc_galaxy_sd_shape_hsm_gauss_gen (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), test->mset, s_data,
                                      std_shape, 0.03, c1, c2, m,
                                      test->ell_coord, rng);
    nc_galaxy_sd_shape_hsm_gauss_data_get (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape), s_data,
                                           &epsilon_1_out, &epsilon_2_out, &std_shape_out, &std_noise_out, &c1_out, &c2_out, &m_out);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);
    int0 = nc_galaxy_sd_shape_integrand_eval (integrand, z_data->z, s_data);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_o = epsilon_1_out + I * epsilon_2_out;
      complex double e_s = 0.0;
      complex double g   = 0.0;
      gdouble gt         = 0.0;
      gdouble theta, phi, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);

      if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
        phi = M_PI - phi;

      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);
        g = (1.0 + m) * g + (c1 + I * c2);

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:

            if (gt > 1.0)
              e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
            else
              e_s = (e_o - g) / (1.0 - conj (g) * e_o);  /* LCOV_EXCL_LINE */

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            e_s = (e_o + g * (g * conj (e_o) - 2.0)) / (1.0 + g * conj (g) - 2.0 * creal (g * conj (e_o)));
            break;

          default: /* LCOV_EXCL_LINE */
            g_assert_not_reached ();
        }
      }
      else
      {
        e_s = e_o - (c1 + I * c2); /* LCOV_EXCL_LINE */
      }

      {
        const gdouble var_int       = gsl_pow_2 (sigma);
        const gdouble total_var     = var_int + gsl_pow_2 (std_noise_out);
        const gdouble chi2_1        = gsl_pow_2 (creal (e_s)) / total_var;
        const gdouble chi2_2        = gsl_pow_2 (cimag (e_s)) / total_var;
        const complex double g_conj = conj (g);
        const gdouble abs_g2        = g * g_conj;
        gdouble jac_num, jac_den, m2ln_int1;

        switch (test->ell_conv)
        {
          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
            jac_num = gsl_pow_2 (1.0 - abs_g2);

            if (abs_g2 <= 1.0)
              jac_den = gsl_pow_2 (1.0 - 2.0 * creal (conj (g) * e_o) + abs_g2 * conj (e_o) * e_o);  /* LCOV_EXCL_LINE */
            else
              jac_den = gsl_pow_2 (conj (e_o - g) * (e_o - g));

            break;

          case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
            jac_num = fabs (gsl_pow_3 (1.0 - abs_g2));
            jac_den = fabs (gsl_pow_3 (1.0 - 2.0 * creal (conj (g) * e_o) + abs_g2));
            break;

          default:
            g_assert_not_reached (); /* LCOV_EXCL_LINE */
            break;
        }

        m2ln_int1 = chi2_1 + chi2_2 + 2.0 * log (2.0 * M_PI * total_var) + 2.0 * log (jac_den / jac_num);

        ncm_assert_cmpdouble_e (-2.0 * int0, ==, m2ln_int1, 1.0e-10, 1.0e-10);
      }
    }
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_sd_shape_integrand_free (integrand);
  g_ptr_array_unref (data_array);
  ncm_rng_free (rng);
}

/* Verifies that eval_at_nodes and integrand_eval produce identical values at each node.
 * Also checks that re-preparing after a parameter change yields different (updated) values. */
static void
_test_nc_galaxy_sd_shape_at_nodes_impl (
  TestNcGalaxySDShape *test,
  gconstpointer       pdata,
  gboolean            is_global
)
{
  NcmRNG *rng                   = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  GPtrArray *z_nodes_per_galaxy = g_ptr_array_new_with_free_func (
    (GDestroyNotify) ncm_vector_free
  );
  NcGalaxySDShapeIntegrand *integrand = nc_galaxy_sd_shape_integ (
    test->galaxy_shape,
    FALSE
  );
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (
    test->galaxy_redshift
  );
  NcGalaxySDPositionData *p_data = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data    = nc_galaxy_sd_shape_data_new (
    test->galaxy_shape,
    p_data
  );
  GPtrArray *data_array    = g_ptr_array_new ();
  const gdouble z_node_min = 0.1;
  const gdouble z_node_max = 1.5;
  const guint n_total      = 100;
  NcmVector *z_nodes       = ncm_vector_new (n_total);
  NcmVector *out           = ncm_vector_new (n_total);
  NcmVector *out_ref       = ncm_vector_new (n_total);
  NcmVector *out_alt       = ncm_vector_new (n_total);
  guint i;

  for (i = 0; i < n_total; i++)
  {
    const gdouble z_i = z_node_min + (z_node_max - z_node_min) * i / (n_total - 1);

    ncm_vector_set (z_nodes, i, z_i);
  }

  g_ptr_array_add (z_nodes_per_galaxy, z_nodes);
  nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);

  z_data->z   = g_test_rand_double_range (0.15, 1.4);
  p_data->ra  = g_test_rand_double_range (-1.0e-2, 1.0e-2);
  p_data->dec = g_test_rand_double_range (-1.0e-2, 1.0e-2);

  if (is_global)
    nc_galaxy_sd_shape_hsm_gauss_global_gen (
      NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
      test->mset,
      s_data,
      0.03,
      0.0,
      0.0,
      0.0,
      test->ell_coord,
      rng);
  else
    nc_galaxy_sd_shape_hsm_gauss_gen (
      NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
      test->mset,
      s_data,
      0.3,
      0.03,
      0.0,
      0.0,
      0.0,
      test->ell_coord,
      rng
    );

  g_ptr_array_add (data_array, s_data);
  nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);

  {
    for (i = 0; i < n_total; i++)
    {
      const gdouble z_i   = ncm_vector_get (z_nodes, i);
      const gdouble ref_i = nc_galaxy_sd_shape_integrand_eval (integrand, z_i, s_data);

      ncm_vector_fast_set (out_ref, i, ref_i);
    }

    nc_galaxy_sd_shape_prepare_at_nodes (
      test->galaxy_shape,
      test->mset,
      data_array,
      z_nodes_per_galaxy
    );
    nc_galaxy_sd_shape_eval_at_nodes (
      test->galaxy_shape,
      test->mset,
      s_data,
      z_nodes,
      out
    );

    for (i = 0; i < n_total; i++)
    {
      const gdouble ref_i = ncm_vector_fast_get (out_ref, i);
      const gdouble out_i = ncm_vector_fast_get (out, i);

      ncm_assert_cmpdouble_e (out_i, ==, ref_i, 1.0e-12, 0.0);
    }
  }

  /* Stale-cache guard: change log10MDelta and verify values change */
  {
    const gdouble mass_orig = ncm_model_param_get_by_name (NCM_MODEL (test->hms), "log10MDelta", NULL);
    const gdouble mass_alt  = mass_orig - 0.5;

    ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", mass_alt, NULL);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);
    nc_galaxy_sd_shape_prepare_at_nodes (test->galaxy_shape, test->mset, data_array, z_nodes_per_galaxy);
    nc_galaxy_sd_shape_eval_at_nodes (test->galaxy_shape, test->mset, s_data, z_nodes, out_alt);

    nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);

    for (i = 0; i < n_total; i++)
    {
      const gdouble z_i       = ncm_vector_get (z_nodes, i);
      const gdouble ref_i     = nc_galaxy_sd_shape_integrand_eval (integrand, z_i, s_data);
      const gdouble out_i     = ncm_vector_get (out, i);
      const gdouble out_alt_i = ncm_vector_get (out_alt, i);

      /* Values must have changed */
      if (z_i > nc_halo_position_get_redshift (test->halo_position))
        g_assert_cmpfloat (out_alt_i, !=, out_i);

      /* And must still agree with integrand_eval */
      ncm_assert_cmpdouble_e (out_alt_i, ==, ref_i, 1.0e-12, 0.0);
    }

    /* Restore */
    ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", mass_orig, NULL);
  }

  /* Stale-cache guard: change H0 and verify eval_at_nodes stays consistent with integrand_eval. */
  {
    const gdouble h0_orig = ncm_model_param_get_by_name (NCM_MODEL (test->cosmo), "H0", NULL);
    const gdouble h0_alt  = h0_orig + 5.0;

    ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "H0", h0_alt, NULL);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);
    nc_galaxy_sd_shape_prepare_at_nodes (test->galaxy_shape, test->mset, data_array, z_nodes_per_galaxy);
    nc_galaxy_sd_shape_eval_at_nodes (test->galaxy_shape, test->mset, s_data, z_nodes, out_alt);

    nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);

    for (i = 0; i < n_total; i++)
    {
      const gdouble z_i       = ncm_vector_get (z_nodes, i);
      const gdouble ref_i     = nc_galaxy_sd_shape_integrand_eval (integrand, z_i, s_data);
      const gdouble out_i     = ncm_vector_get (out, i);
      const gdouble out_alt_i = ncm_vector_get (out_alt, i);

      /* Values must have changed */
      if (z_i > nc_halo_position_get_redshift (test->halo_position))
        g_assert_cmpfloat (out_alt_i, !=, out_i);

      ncm_assert_cmpdouble_e (out_alt_i, ==, ref_i, 1.0e-6, 0.0);
    }

    /* Restore */
    ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "H0", h0_orig, NULL);
  }

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_sd_shape_integrand_free (integrand);
  ncm_vector_free (out);
  ncm_vector_free (out_ref);
  ncm_vector_free (out_alt);
  g_ptr_array_unref (data_array);
  g_ptr_array_unref (z_nodes_per_galaxy);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_at_nodes (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  _test_nc_galaxy_sd_shape_at_nodes_impl (test, pdata, TRUE);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_at_nodes (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  _test_nc_galaxy_sd_shape_at_nodes_impl (test, pdata, FALSE);
}

/*
 * Deterministic test of the direct shear estimator. The observed ellipticities
 * are constructed so that the per-galaxy rotated shape hat_g = e_o * exp (-2 i phi)
 * takes prescribed values whose sample mean is (INJ_GT, INJ_GX) and whose real and
 * imaginary parts share the *same* zero-mean offset sequence (so sd_re == sd_im
 * exactly). With no bias (m = c = 0) the estimator must return, for the responsivity
 * R = 1 - std_shape^2:
 *   - TRACE_DET (epsilon):    gt = INJ_GT,           gx = INJ_GX           (no R factor);
 *   - TRACE     (distortion): gt = 0.5 * INJ_GT / R, gx = 0.5 * INJ_GX / R (R on both).
 * In both conventions sigma_t / sigma_x == 1, since the responsivity (when present)
 * multiplies both components equally. This pins down the R placement that the
 * tangential/cross asymmetry bug got wrong.
 */
#define TEST_ESTIMATE_INJ_GT 0.2
#define TEST_ESTIMATE_INJ_GX 0.1

static void
test_nc_galaxy_sd_shape_estimate_check (TestNcGalaxySDShape *test, GPtrArray *data_array, const gdouble R)
{
  gdouble hat_gt, hat_gx, hat_sigma_t, hat_sigma_x, hat_rho;
  gdouble exp_gt, exp_gx;

  nc_galaxy_sd_shape_direct_estimate (test->galaxy_shape, test->mset, data_array,
                                      &hat_gt, &hat_gx, &hat_sigma_t, &hat_sigma_x, &hat_rho);

  switch (test->ell_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      exp_gt = TEST_ESTIMATE_INJ_GT;
      exp_gx = TEST_ESTIMATE_INJ_GX;
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      exp_gt = 0.5 * TEST_ESTIMATE_INJ_GT / R;
      exp_gx = 0.5 * TEST_ESTIMATE_INJ_GX / R;
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  g_assert_true (gsl_finite (hat_gt) && gsl_finite (hat_gx));
  g_assert_true (gsl_finite (hat_sigma_t) && gsl_finite (hat_sigma_x));
  g_assert_cmpfloat (hat_sigma_x, >, 0.0);

  ncm_assert_cmpdouble_e (hat_gt, ==, exp_gt, 1.0e-9, 0.0);
  ncm_assert_cmpdouble_e (hat_gx, ==, exp_gx, 1.0e-9, 0.0);
  /* Responsivity is isotropic: it cancels in the ratio of the two components. */
  ncm_assert_cmpdouble_e (hat_sigma_t / hat_sigma_x, ==, 1.0, 1.0e-9, 0.0);
}

/* hat_g for galaxy i: prescribed mean plus a deterministic symmetric (zero-mean) spread,
 * identical on the real and imaginary axes so that sd_re == sd_im. */
static complex double
test_nc_galaxy_sd_shape_estimate_hat_g (const guint i, const guint ngals, const gdouble std_shape)
{
  const gdouble off = std_shape * (2.0 * i / (ngals - 1.0) - 1.0);

  return (TEST_ESTIMATE_INJ_GT + off) + I * (TEST_ESTIMATE_INJ_GX + off);
}

/* e_o that yields the prescribed hat_g at the galaxy position, inverting the rotation
 * (and Euclidean flip) the estimator applies. */
static complex double
test_nc_galaxy_sd_shape_estimate_e_o (TestNcGalaxySDShape *test, const gdouble ra, const gdouble dec, const complex double hat_g)
{
  gdouble theta, phi;

  nc_halo_position_polar_angles (test->halo_position, ra, dec, &theta, &phi);

  if (test->ell_coord == NC_GALAXY_WL_OBS_COORD_EUCLIDEAN)
    phi = M_PI - phi;

  return hat_g * cexp (2.0 * I * phi);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_estimate (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  GPtrArray *data_array   = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_sd_shape_data_unref);
  const guint ngals       = 1000;
  const gdouble sigma     = ncm_model_orig_param_get (NCM_MODEL (test->galaxy_shape), NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SIGMA);
  const gdouble std_shape = nc_galaxy_sd_shape_hsm_gauss_global_std_shape_from_sigma (sigma);
  const gdouble R         = 1.0 - std_shape * std_shape;
  const gdouble ra        = 0.05;
  const gdouble dec       = -0.04;
  guint i;

  for (i = 0; i < ngals; i++)
  {
    NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
    NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
    NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
    const complex double hat_g        = test_nc_galaxy_sd_shape_estimate_hat_g (i, ngals, std_shape);
    complex double e_o;

    p_data->ra    = ra;
    p_data->dec   = dec;
    s_data->coord = test->ell_coord;

    e_o = test_nc_galaxy_sd_shape_estimate_e_o (test, ra, dec, hat_g);

    nc_galaxy_sd_shape_hsm_gauss_global_data_set (NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (test->galaxy_shape),
                                                  s_data, creal (e_o), cimag (e_o), 0.1, 0.0, 0.0, 0.0);

    g_ptr_array_add (data_array, s_data);
    nc_galaxy_sd_obs_redshift_data_unref (z_data);
    nc_galaxy_sd_position_data_unref (p_data);
  }

  test_nc_galaxy_sd_shape_estimate_check (test, data_array, R);

  g_ptr_array_unref (data_array);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_estimate (TestNcGalaxySDShape *test, gconstpointer pdata)
{
  GPtrArray *data_array   = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_sd_shape_data_unref);
  const guint ngals       = 1000;
  const gdouble std_shape = 0.3;
  const gdouble R         = 1.0 - std_shape * std_shape;
  const gdouble ra        = 0.05;
  const gdouble dec       = -0.04;
  guint i;

  for (i = 0; i < ngals; i++)
  {
    NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
    NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
    NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
    const complex double hat_g        = test_nc_galaxy_sd_shape_estimate_hat_g (i, ngals, std_shape);
    complex double e_o;

    p_data->ra    = ra;
    p_data->dec   = dec;
    s_data->coord = test->ell_coord;

    e_o = test_nc_galaxy_sd_shape_estimate_e_o (test, ra, dec, hat_g);

    nc_galaxy_sd_shape_hsm_gauss_data_set (NC_GALAXY_SD_SHAPE_HSM_GAUSS (test->galaxy_shape),
                                           s_data, creal (e_o), cimag (e_o), std_shape, 0.1, 0.0, 0.0, 0.0);

    g_ptr_array_add (data_array, s_data);
    nc_galaxy_sd_obs_redshift_data_unref (z_data);
    nc_galaxy_sd_position_data_unref (p_data);
  }

  test_nc_galaxy_sd_shape_estimate_check (test, data_array, R);

  g_ptr_array_unref (data_array);
}

static void
test_nc_galaxy_sd_shape_hsm_gauss_global_sigma_conversions (void)
{
  gint i;

  for (i = 0; i < 10000; i++)
  {
    const gdouble sigma      = g_test_rand_double_range (0.01, 2.0);
    const gdouble std_shape  = nc_galaxy_sd_shape_hsm_gauss_global_std_shape_from_sigma (sigma);
    const gdouble sigma2     = nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (std_shape);
    const gdouble std_shape2 = nc_galaxy_sd_shape_hsm_gauss_global_std_shape_from_sigma (sigma2);

    g_assert (gsl_finite (sigma));
    g_assert (gsl_finite (std_shape));
    ncm_assert_cmpdouble_e (sigma, ==, sigma2, 1e-12, 0.0);
    ncm_assert_cmpdouble_e (std_shape, ==, std_shape2, 1e-12, 0.0);
  }
}


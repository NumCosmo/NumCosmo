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
  NcHICosmo *cosmo;
  NcHaloMassSummary *hms;
  NcHaloDensityProfile *density_profile;
  NcHaloPosition *hp;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcGalaxySDObsRedshift *galaxy_redshift;
  NcGalaxySDPosition *galaxy_position;
  NcGalaxySDShape *galaxy_shape;
  NcmMSet *mset;
  NcGalaxyWLObsCoord ell_coord;
  NcGalaxyWLObsEllipConv ell_conv;
  gdouble min_radius;
  gdouble max_radius;
} TestNcDataClusterWL;

typedef struct _TestNcDataClusterWLTests
{
  char *test_name;

  void (*test_func) (TestNcDataClusterWL *test, gconstpointer pdata);
} TestNcDataClusterWLTests;

typedef struct _TestNcDataClusterWLTestsObj
{
  gchar *shape_name;
  gchar *redshift_name;
  gchar *ell_conv_name;
  gchar *ell_coord_name;
} TestNcDataClusterWLTestsObj;

static void test_nc_data_cluster_wl_new (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_gen (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_gen_obs (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_resample (TestNcDataClusterWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  TestNcDataClusterWLTests tests[5] = {
    {"gen_obs", &test_nc_data_cluster_wl_gen_obs},
    {"m2lnP", &test_nc_data_cluster_wl_m2lnP},
    {"serialize", &test_nc_data_cluster_wl_serialize},
    {"resample", &test_nc_data_cluster_wl_resample},
  };
  TestNcDataClusterWLTestsObj tests_obj[24] = {
    {"gauss", "spec", "trace", "celestial"},
    {"gauss", "spec", "trace", "euclidean"},
    {"gauss", "spec", "trace_det", "celestial"},
    {"gauss", "spec", "trace_det", "euclidean"},
    {"gauss", "gauss", "trace", "celestial"},
    {"gauss", "gauss", "trace", "euclidean"},
    {"gauss", "gauss", "trace_det", "celestial"},
    {"gauss", "gauss", "trace_det", "euclidean"},
    {"gauss", "pz", "trace", "celestial"},
    {"gauss", "pz", "trace", "euclidean"},
    {"gauss", "pz", "trace_det", "celestial"},
    {"gauss", "pz", "trace_det", "euclidean"},
    {"gauss_hsc", "spec", "trace", "celestial"},
    {"gauss_hsc", "spec", "trace", "euclidean"},
    {"gauss_hsc", "spec", "trace_det", "celestial"},
    {"gauss_hsc", "spec", "trace_det", "euclidean"},
    {"gauss_hsc", "gauss", "trace", "celestial"},
    {"gauss_hsc", "gauss", "trace", "euclidean"},
    {"gauss_hsc", "gauss", "trace_det", "celestial"},
    {"gauss_hsc", "gauss", "trace_det", "euclidean"},
    {"gauss_hsc", "pz", "trace", "celestial"},
    {"gauss_hsc", "pz", "trace", "euclidean"},
    {"gauss_hsc", "pz", "trace_det", "celestial"},
    {"gauss_hsc", "pz", "trace_det", "euclidean"}
  };
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < 24; i++)
  {
    for (j = 0; j < 4; j++)
    {
      gchar *test_name = g_strdup_printf ("/nc/data_cluster_wl/%s/%s/%s/%s/%s",
                                          tests_obj[i].shape_name,
                                          tests_obj[i].redshift_name,
                                          tests_obj[i].ell_conv_name,
                                          tests_obj[i].ell_coord_name,
                                          tests[j].test_name);

      g_test_add (test_name, TestNcDataClusterWL, (gconstpointer) & tests_obj[i],
                  &test_nc_data_cluster_wl_new,
                  tests[j].test_func,
                  test_nc_data_cluster_wl_free);
      g_free (test_name);
    }
  }

  g_test_run ();

  return 0;
}

static void
test_nc_data_cluster_wl_new (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloMassSummary *hms              = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);
  TestNcDataClusterWLTestsObj *obj    = (TestNcDataClusterWLTestsObj *) pdata;
  gchar *shape_name                   = obj->shape_name;
  gchar *redshift_name                = obj->redshift_name;
  gchar *ell_conv_name                = obj->ell_conv_name;
  gchar *ell_coord_name               = obj->ell_coord_name;
  gdouble min_radius                  = ncm_rng_uniform_gen (rng, 0.1, 0.5);
  gdouble max_radius                  = ncm_rng_uniform_gen (rng, 2.0, 5.0);
  gdouble ra                          = ncm_rng_uniform_gen (rng, -180, 180);
  gdouble dec                         = ncm_rng_uniform_gen (rng, -90, 90);
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDPosition *p_dist;
  NcGalaxySDObsRedshift *z_dist;
  NcGalaxySDShape *s_dist;
  NcGalaxyWLObsEllipConv ell_conv;
  NcGalaxyWLObsCoord ell_coord;

  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, ncm_rng_uniform_gen (rng, 13.0, 15.0));
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_RA, ra);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_DEC, dec);

  ncm_model_param_set_lower_bound (NCM_MODEL (hp), NC_HALO_POSITION_RA, -180.0);
  ncm_model_param_set_upper_bound (NCM_MODEL (hp), NC_HALO_POSITION_RA, 180.0);
  ncm_model_param_set_lower_bound (NCM_MODEL (hp), NC_HALO_POSITION_DEC, -90.0);
  ncm_model_param_set_upper_bound (NCM_MODEL (hp), NC_HALO_POSITION_DEC, 90.0);
  ncm_model_param_set_lower_bound (NCM_MODEL (hp), NC_HALO_POSITION_RA, ra - 0.025 / cos (dec * M_PI / 180.0));
  ncm_model_param_set_upper_bound (NCM_MODEL (hp), NC_HALO_POSITION_RA, ra + 0.025 / cos (dec * M_PI / 180.0));
  ncm_model_param_set_lower_bound (NCM_MODEL (hp), NC_HALO_POSITION_DEC, dec - 0.025);
  ncm_model_param_set_upper_bound (NCM_MODEL (hp), NC_HALO_POSITION_DEC, dec + 0.025);

  nc_halo_position_prepare (hp, cosmo);

  nc_data_cluster_wl_set_cut (dcwl, min_radius, max_radius);

  p_dist = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (ra - 0.2, ra + 0.2, dec - 0.2, dec + 0.2));

  if (g_strcmp0 (ell_conv_name, "trace") == 0)
    ell_conv = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;
  else if (g_strcmp0 (ell_conv_name, "trace_det") == 0)
    ell_conv = NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET;
  else
    g_error ("test_nc_data_cluster_wl_new_spec_gauss: unknown ellip_conv name.");

  if (g_strcmp0 (ell_coord_name, "celestial") == 0)
    ell_coord = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  else if (g_strcmp0 (ell_coord_name, "euclidean") == 0)
    ell_coord = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;
  else
    g_error ("test_nc_data_cluster_wl_new_spec_gauss: unknown ellip_conv name.");

  if (g_strcmp0 (shape_name, "gauss") == 0)
    s_dist = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new (ell_conv));
  else if (g_strcmp0 (shape_name, "gauss_hsc") == 0)
    s_dist = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_hsc_new (ell_conv));
  else
    g_error ("test_nc_data_cluster_wl_new_spec_gauss: unknown shape name.");

  if (g_strcmp0 (redshift_name, "spec") == 0)
    z_dist = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true_dist, 0.0, 5.0));
  else if (g_strcmp0 (redshift_name, "gauss") == 0)
    z_dist = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_gauss_new (z_true_dist, 0.1, 4.8));
  else if (g_strcmp0 (redshift_name, "pz") == 0)
    z_dist = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_pz_new ());
  else
    g_error ("test_nc_data_cluster_wl_new_spec_gauss: unknown redshift name.");

  test->dcwl  = dcwl;
  test->cosmo = cosmo;

  test->hms                  = hms;
  test->density_profile      = dp;
  test->hp                   = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  test->ell_conv  = ell_conv;
  test->ell_coord = ell_coord;

  test->min_radius = min_radius;
  test->max_radius = max_radius;

  test_nc_data_cluster_wl_gen (test, pdata);

  nc_galaxy_sd_true_redshift_free (z_true_dist);
}

static void
test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata)
{
  nc_hicosmo_clear (&test->cosmo);
  nc_halo_density_profile_clear (&test->density_profile);
  nc_halo_mass_summary_clear (&test->hms);
  nc_halo_position_clear (&test->hp);
  nc_wl_surface_mass_density_clear (&test->surface_mass_density);

  nc_galaxy_sd_obs_redshift_clear (&test->galaxy_redshift);
  nc_galaxy_sd_position_clear (&test->galaxy_position);
  nc_galaxy_sd_shape_clear (&test->galaxy_shape);

  ncm_mset_free (test->mset);

  NCM_TEST_FREE (nc_data_cluster_wl_free, test->dcwl);
}

static void
test_nc_data_cluster_wl_gen (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  GList *columns                    = nc_galaxy_sd_shape_data_required_columns (s_data);
  GList *l                          = columns;
  GStrvBuilder *builder             = g_strv_builder_new ();
  guint nrows                       = 1000;
  guint npoints                     = 20;
  gdouble z_min                     = 0.01;
  gdouble z_max                     = 5.0;
  NcGalaxyWLObs *obs;
  GStrv columns_strv;
  guint i;

  while (l)
  {
    g_strv_builder_add (builder, l->data);

    l = g_list_next (l);
  }

  columns_strv = g_strv_builder_end (builder);
  obs          = nc_galaxy_wl_obs_new (test->ell_conv, test->ell_coord, nrows, columns_strv);

  for (i = 0; i < nrows; i++)
  {
    if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
    {
      nc_galaxy_sd_obs_redshift_gauss_gen (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift), test->mset, z_data, 0.03, rng);
    }
    else if (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (test->galaxy_redshift))
    {
      nc_galaxy_sd_obs_redshift_spec_gen (NC_GALAXY_SD_OBS_REDSHIFT_SPEC (test->galaxy_redshift), test->mset, z_data, rng);
    }
    else if (NC_IS_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift))
    {
      NcmVector *xv = ncm_vector_new (npoints);
      NcmVector *yv = ncm_vector_new (npoints);
      gdouble z_avg = g_test_rand_double_range (z_min + 0.5, z_max - 0.5);
      gdouble z_sd  = 0.03 * (1.0 + z_avg);
      gdouble x_min = z_avg - 5.0 * z_sd;
      gdouble x_max = z_avg + 5.0 * z_sd;
      NcmSpline *pz;
      guint j;

      if (x_min < z_min)
        x_min = z_min;  /* LCOV_EXCL_LINE */

      if (x_max > z_max)
        x_max = z_max;  /* LCOV_EXCL_LINE */

      for (j = 0; j < npoints; j++)
      {
        gdouble x = x_min + (x_max - x_min) * j / ((gdouble) npoints - 1.0);
        gdouble y = exp (-0.5 * gsl_pow_2 ((x - z_avg) / z_sd)) / (sqrt (2.0 * M_PI) * z_sd);

        ncm_vector_fast_set (xv, j, x);
        ncm_vector_fast_set (yv, j, y);
      }

      pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

      nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift), z_data, pz);
      nc_galaxy_sd_obs_redshift_pz_prepare (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift), z_data);
      nc_galaxy_sd_obs_redshift_pz_gen (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift), test->mset, z_data, rng);

      ncm_vector_free (xv);
      ncm_vector_free (yv);
      ncm_spline_free (pz);
    }
    else
    {
      g_error ("test_nc_data_cluster_wl_gen: unknown redshift name.");
    }

    gdouble radius = 0.0;

    do {
      nc_galaxy_sd_position_gen (test->galaxy_position, p_data, rng);

      radius = nc_halo_position_projected_radius_from_ra_dec (test->hp, test->cosmo, p_data->ra, p_data->dec);
    } while (radius < test->min_radius || radius > test->max_radius);

    if (NC_IS_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape))
    {
      gdouble std_noise = ncm_rng_uniform_gen (rng, 0.0, 0.3);

      nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, std_noise, test->ell_coord, rng);
    }
    else if (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape))
    {
      gdouble c1        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
      gdouble c2        = ncm_rng_uniform_gen (rng, -0.01, 0.01);
      gdouble m         = ncm_rng_uniform_gen (rng, -0.2, 0.2);
      gdouble std_shape = ncm_rng_uniform_gen (rng, 0.15, 0.3);
      gdouble std_noise = ncm_rng_uniform_gen (rng, 0.0, 0.3);

      nc_galaxy_sd_shape_gauss_hsc_gen (NC_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape), test->mset, s_data, std_shape, std_noise, c1, c2, m, test->ell_coord, rng);
    }
    else
    {
      g_error ("test_nc_data_cluster_wl_gen: unknown shape name.");
    }

    nc_galaxy_sd_shape_data_write_row (s_data, obs, i);
  }

  nc_data_cluster_wl_set_obs (test->dcwl, obs);

  g_strfreev (columns_strv);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_wl_obs_free (obs);
  ncm_rng_free (rng);
  g_strv_builder_unref (builder);
  g_list_free (columns);
  g_list_free (l);
}

static void
test_nc_data_cluster_wl_gen_obs (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxyWLObs *obs;
  gdouble ngals;
  guint i;

  obs   = nc_data_cluster_wl_peek_obs (test->dcwl);
  ngals = nc_galaxy_wl_obs_len (obs);

  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_OBS_REDSHIFT_COL_Z));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_POSITION_COL_RA));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_POSITION_COL_DEC));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE));

  if (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape))
  {
    g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1));
    g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2));
    g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M));
    g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE));
  }

  if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
  {
    g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA));
    g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0));
  }

  for (i = 0; i < ngals; i++)
  {
    const gdouble z             = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
    const gdouble ra            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_RA, i);
    const gdouble dec           = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_DEC, i);
    const gdouble ra_cl         = ncm_model_param_get (NCM_MODEL (test->hp), NC_HALO_POSITION_RA);
    const gdouble dec_cl        = ncm_model_param_get (NCM_MODEL (test->hp), NC_HALO_POSITION_DEC);
    const gdouble epsilon_int_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
    const gdouble epsilon_int_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
    const gdouble std_noise     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);
    const gdouble var_obs       = std_noise * std_noise;
    gdouble sigma;

    if (NC_IS_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape))
      sigma = ncm_model_orig_param_get (NCM_MODEL (test->galaxy_shape), NC_GALAXY_SD_SHAPE_GAUSS_SIGMA);
    else
      sigma = nc_galaxy_sd_shape_gauss_sigma_from_std_shape (nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i));

    const gdouble var_int = sigma * sigma;
    const gdouble e_rms   = sqrt (var_int + var_obs);
    const gdouble radius  = nc_halo_position_projected_radius_from_ra_dec (test->hp, test->cosmo, ra, dec);

    g_assert_cmpfloat (ra, >=, ra_cl - 0.2);
    g_assert_cmpfloat (ra, <=, ra_cl + 0.2);
    g_assert_cmpfloat (dec, >=, dec_cl - 0.2);
    g_assert_cmpfloat (dec, <=, dec_cl + 0.2);
    g_assert_cmpfloat (z, >=, 0.0);
    g_assert_cmpfloat (z, <=, 10.0);
    g_assert_cmpfloat (epsilon_int_1, >=, -5.0 * e_rms);
    g_assert_cmpfloat (epsilon_int_1, <=, 5.0 * e_rms);
    g_assert_cmpfloat (epsilon_int_2, >=, -5.0 * e_rms);
    g_assert_cmpfloat (epsilon_int_2, <=, 5.0 * e_rms);
    g_assert_cmpfloat (radius, >=, test->min_radius);
    g_assert_cmpfloat (radius, <=, test->max_radius);
  }

  ncm_rng_free (rng);
}

static void
test_nc_data_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata)
{
  gdouble m2lnL_a, m2lnL_b;

  ncm_data_prepare (NCM_DATA (test->dcwl), test->mset);

  {
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_a);
    g_assert (gsl_finite (m2lnL_a));

    ncm_mset_param_set (test->mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, 0.3);

    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_b);
    g_assert (gsl_finite (m2lnL_a));

    ncm_assert_cmpdouble (m2lnL_a, !=, m2lnL_b);
  }

  {
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_a);
    g_assert (gsl_finite (m2lnL_a));

    ncm_mset_param_set (test->mset, nc_halo_position_id (), NC_HALO_POSITION_RA, 0.2);

    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_b);
    g_assert (gsl_finite (m2lnL_a));

    ncm_assert_cmpdouble (m2lnL_a, !=, m2lnL_b);
  }

  {
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_a);
    g_assert (gsl_finite (m2lnL_a));

    ncm_mset_param_set (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, log10 (2.123e14));

    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_b);
    g_assert (gsl_finite (m2lnL_a));

    ncm_assert_cmpdouble (m2lnL_a, !=, m2lnL_b);
  }

  {
    g_object_set (test->dcwl, "enable-parallel", FALSE, NULL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_a);
    g_assert (gsl_finite (m2lnL_a));

    g_object_set (test->dcwl, "enable-parallel", TRUE, NULL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_b);
    g_assert (gsl_finite (m2lnL_b));

    ncm_assert_cmpdouble_e (m2lnL_a, ==, m2lnL_b, 1.0e-11, 0.0);
  }
}

static void
test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

  {
    NcmMSet *mset_dup = ncm_mset_dup (test->mset, ser);
    NcmData *dcwl_dup = ncm_data_dup (NCM_DATA (test->dcwl), ser);
    gdouble m2lnL, m2lnL_dup;

    g_object_set (dcwl_dup, "enable-parallel", FALSE, NULL);
    g_object_set (test->dcwl, "enable-parallel", FALSE, NULL);

    g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl_dup));

    ncm_data_m2lnL_val (dcwl_dup, mset_dup, &m2lnL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_dup);

    ncm_assert_cmpdouble (m2lnL, ==, m2lnL_dup);

    ncm_mset_param_set (test->mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, 0.3);
    ncm_mset_param_set (mset_dup, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, 0.3);

    ncm_data_m2lnL_val (dcwl_dup, mset_dup, &m2lnL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_dup);

    ncm_assert_cmpdouble (m2lnL, ==, m2lnL_dup);

    ncm_mset_param_set (test->mset, nc_halo_position_id (), NC_HALO_POSITION_DEC, 0.1);
    ncm_mset_param_set (mset_dup, nc_halo_position_id (), NC_HALO_POSITION_DEC, 0.1);

    ncm_data_m2lnL_val (dcwl_dup, mset_dup, &m2lnL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_dup);

    ncm_assert_cmpdouble (m2lnL, ==, m2lnL_dup);

    ncm_mset_param_set (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, log10 (2.123e14));
    ncm_mset_param_set (mset_dup, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, log10 (2.123e14));

    ncm_data_m2lnL_val (dcwl_dup, mset_dup, &m2lnL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), test->mset, &m2lnL_dup);

    ncm_assert_cmpdouble (m2lnL, ==, m2lnL_dup);

    ncm_mset_param_set (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, log10 (1.123e14));
    ncm_mset_param_set (mset_dup, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, log10 (1.123e14));

    ncm_data_m2lnL_val (dcwl_dup, test->mset, &m2lnL);
    ncm_data_m2lnL_val (NCM_DATA (test->dcwl), mset_dup, &m2lnL_dup);

    ncm_assert_cmpdouble (m2lnL, ==, m2lnL_dup);

    g_object_set (dcwl_dup, "enable-parallel", TRUE, NULL);
    g_object_set (test->dcwl, "enable-parallel", TRUE, NULL);

    ncm_mset_free (mset_dup);
    ncm_data_free (dcwl_dup);
  }

  ncm_serialize_free (ser);
}

static void
test_nc_data_cluster_wl_resample (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmData *data                     = NCM_DATA (test->dcwl);
  NcmRNG *rng                       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxyWLObs *obs                = nc_data_cluster_wl_peek_obs (test->dcwl);
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  GStrv columns                     = nc_galaxy_wl_obs_peek_columns (obs);
  GList *l                          = nc_galaxy_sd_shape_data_required_columns (s_data);
  gdouble ngals                     = nc_galaxy_wl_obs_len (obs);
  NcGalaxyWLObs *obs_copy           = nc_galaxy_wl_obs_new (test->ell_conv, test->ell_coord, ngals, columns);
  NcGalaxyWLObs *obs2;
  guint i;

  while (l)
  {
    const gchar *col_name = l->data;

    for (i = 0; i < ngals; i++)
    {
      const gdouble val = nc_galaxy_wl_obs_get (obs, col_name, i);

      nc_galaxy_wl_obs_set (obs_copy, col_name, i, val);
    }

    l = g_list_next (l);
  }

  /* PARAM_FLAG = ALL */
  nc_data_cluster_wl_set_resample_flag (test->dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL);

  {
    ncm_data_resample (data, test->mset, rng);

    obs2 = nc_data_cluster_wl_peek_obs (test->dcwl);

    for (i = 0; i < ngals; i++)
    {
      const gdouble ra         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec        = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);
      const gdouble ra2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z2         = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);

      g_assert_cmpfloat (ra, !=, ra2);
      g_assert_cmpfloat (dec, !=, dec2);
      g_assert_cmpfloat (z, !=, z2);
      g_assert_cmpfloat (e1_int, !=, e1_int2);
      g_assert_cmpfloat (e2_int, !=, e2_int2);
      g_assert_cmpfloat (e1_obs, !=, e1_obs2);
      g_assert_cmpfloat (e2_obs, !=, e2_obs2);
      g_assert_cmpfloat (std_noise, ==, std_noise2);

      if (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape))
      {
        const gdouble std_shape  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);
        const gdouble std_shape2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m_2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);

        g_assert_cmpfloat (std_shape, ==, std_shape2);
        g_assert_cmpfloat (c1, ==, c1_2);
        g_assert_cmpfloat (c2, ==, c2_2);
        g_assert_cmpfloat (m, ==, m_2);
      }

      if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
      {
        const gdouble zp      = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma   = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma0  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);
        const gdouble zp2     = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma2  = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma02 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);

        g_assert_cmpfloat (zp, !=, zp2);
        g_assert_cmpfloat (sigma, !=, sigma2);
        g_assert_cmpfloat (sigma0, ==, sigma02);
      }
    }
  }

  l = nc_galaxy_sd_shape_data_required_columns (s_data);

  while (l)
  {
    const gchar *col_name = l->data;

    for (i = 0; i < ngals; i++)
    {
      const gdouble val = nc_galaxy_wl_obs_get (obs, col_name, i);

      nc_galaxy_wl_obs_set (obs_copy, col_name, i, val);
    }

    l = g_list_next (l);
  }

  /* PARAM_FLAG = REDSHIFT */
  nc_data_cluster_wl_set_resample_flag (test->dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_REDSHIFT);

  {
    ncm_data_resample (data, test->mset, rng);

    obs2 = nc_data_cluster_wl_peek_obs (test->dcwl);

    for (i = 0; i < ngals; i++)
    {
      const gdouble ra         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec        = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);
      const gdouble ra2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z2         = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);

      g_assert_cmpfloat (ra, ==, ra2);
      g_assert_cmpfloat (dec, ==, dec2);
      g_assert_cmpfloat (z, !=, z2);
      g_assert_cmpfloat (e1_int, !=, e1_int2);
      g_assert_cmpfloat (e2_int, !=, e2_int2);
      g_assert_cmpfloat (e1_obs, !=, e1_obs2);
      g_assert_cmpfloat (e2_obs, !=, e2_obs2);
      g_assert_cmpfloat (std_noise, ==, std_noise2);

      if (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape))
      {
        const gdouble std_shape  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);
        const gdouble std_shape2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m_2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);

        g_assert_cmpfloat (std_shape, ==, std_shape2);
        g_assert_cmpfloat (c1, ==, c1_2);
        g_assert_cmpfloat (c2, ==, c2_2);
        g_assert_cmpfloat (m, ==, m_2);
      }

      if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
      {
        const gdouble zp      = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma   = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma0  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);
        const gdouble zp2     = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma2  = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma02 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);

        g_assert_cmpfloat (zp, !=, zp2);
        g_assert_cmpfloat (sigma, !=, sigma2);
        g_assert_cmpfloat (sigma0, ==, sigma02);
      }
    }
  }

  l = nc_galaxy_sd_shape_data_required_columns (s_data);

  while (l)
  {
    const gchar *col_name = l->data;

    for (i = 0; i < ngals; i++)
    {
      const gdouble val = nc_galaxy_wl_obs_get (obs, col_name, i);

      nc_galaxy_wl_obs_set (obs_copy, col_name, i, val);
    }

    l = g_list_next (l);
  }

  /* PARAM_FLAG = POSITION */
  nc_data_cluster_wl_set_resample_flag (test->dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_POSITION);

  {
    ncm_data_resample (data, test->mset, rng);

    obs2 = nc_data_cluster_wl_peek_obs (test->dcwl);

    for (i = 0; i < ngals; i++)
    {
      const gdouble ra         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec        = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);
      const gdouble ra2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z2         = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);

      g_assert_cmpfloat (ra, !=, ra2);
      g_assert_cmpfloat (dec, !=, dec2);
      g_assert_cmpfloat (z, ==, z2);
      g_assert_cmpfloat (e1_int, !=, e1_int2);
      g_assert_cmpfloat (e2_int, !=, e2_int2);
      g_assert_cmpfloat (e1_obs, !=, e1_obs2);
      g_assert_cmpfloat (e2_obs, !=, e2_obs2);
      g_assert_cmpfloat (std_noise, ==, std_noise2);

      if (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape))
      {
        const gdouble std_shape  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);
        const gdouble std_shape2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m_2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);

        g_assert_cmpfloat (std_shape, ==, std_shape2);
        g_assert_cmpfloat (c1, ==, c1_2);
        g_assert_cmpfloat (c2, ==, c2_2);
        g_assert_cmpfloat (m, ==, m_2);
      }

      if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
      {
        const gdouble zp      = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma   = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma0  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);
        const gdouble zp2     = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma2  = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma02 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);

        g_assert_cmpfloat (zp, ==, zp2);
        g_assert_cmpfloat (sigma, ==, sigma2);
        g_assert_cmpfloat (sigma0, ==, sigma02);
      }
    }
  }

  l = nc_galaxy_sd_shape_data_required_columns (s_data);

  while (l)
  {
    const gchar *col_name = l->data;

    for (i = 0; i < ngals; i++)
    {
      const gdouble val = nc_galaxy_wl_obs_get (obs, col_name, i);

      nc_galaxy_wl_obs_set (obs_copy, col_name, i, val);
    }

    l = g_list_next (l);
  }

  /* PARAM_FLAG = SHAPE */
  nc_data_cluster_wl_set_resample_flag (test->dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_SHAPE);

  {
    ncm_data_resample (data, test->mset, rng);

    obs2 = nc_data_cluster_wl_peek_obs (test->dcwl);

    for (i = 0; i < ngals; i++)
    {
      const gdouble ra         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec        = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs     = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);
      const gdouble ra2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_RA, i);
      const gdouble dec2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_POSITION_COL_DEC, i);
      const gdouble z2         = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
      const gdouble e1_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
      const gdouble e2_int2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
      const gdouble e1_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1, i);
      const gdouble e2_obs2    = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2, i);
      const gdouble std_noise2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_COL_STD_NOISE, i);

      g_assert_cmpfloat (ra, ==, ra2);
      g_assert_cmpfloat (dec, ==, dec2);
      g_assert_cmpfloat (z, ==, z2);
      g_assert_cmpfloat (e1_int, !=, e1_int2);
      g_assert_cmpfloat (e2_int, !=, e2_int2);
      g_assert_cmpfloat (e1_obs, !=, e1_obs2);
      g_assert_cmpfloat (e2_obs, !=, e2_obs2);
      g_assert_cmpfloat (std_noise, ==, std_noise2);

      if (NC_IS_GALAXY_SD_SHAPE_GAUSS_HSC (test->galaxy_shape))
      {
        const gdouble std_shape  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2         = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m          = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);
        const gdouble std_shape2 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_STD_SHAPE, i);
        const gdouble c1_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C1, i);
        const gdouble c2_2       = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_C2, i);
        const gdouble m_2        = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_M, i);

        g_assert_cmpfloat (std_shape, ==, std_shape2);
        g_assert_cmpfloat (c1, ==, c1_2);
        g_assert_cmpfloat (c2, ==, c2_2);
        g_assert_cmpfloat (m, ==, m_2);
      }

      if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
      {
        const gdouble zp      = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma   = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma0  = nc_galaxy_wl_obs_get (obs_copy, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);
        const gdouble zp2     = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
        const gdouble sigma2  = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
        const gdouble sigma02 = nc_galaxy_wl_obs_get (obs2, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);

        g_assert_cmpfloat (zp, ==, zp2);
        g_assert_cmpfloat (sigma, ==, sigma2);
        g_assert_cmpfloat (sigma0, ==, sigma02);
      }
    }
  }

  ncm_rng_clear (&rng);
  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  g_list_free (l);
  nc_galaxy_wl_obs_free (obs_copy);
}


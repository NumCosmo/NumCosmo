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
} TestNcDataClusterWLTestsObj;

static void test_nc_data_cluster_wl_new (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_gen (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_monte_carlo (TestNcDataClusterWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  TestNcDataClusterWLTests tests[5] = {
    {"monte_carlo", &test_nc_data_cluster_wl_monte_carlo}
  };
  TestNcDataClusterWLTestsObj tests_obj[12] = {
    {"gauss", "spec", "trace"},
    {"gauss", "spec", "trace_det"},
    {"gauss", "gauss", "trace"},
    {"gauss", "gauss", "trace_det"},
    {"gauss", "pz", "trace"},
    {"gauss", "pz", "trace_det"},
    {"gauss_hsc", "spec", "trace"},
    {"gauss_hsc", "spec", "trace_det"},
    {"gauss_hsc", "gauss", "trace"},
    {"gauss_hsc", "gauss", "trace_det"},
    {"gauss_hsc", "pz", "trace"},
    {"gauss_hsc", "pz", "trace_det"},
  };
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < 12; i++)
  {
    for (j = 0; j < 1; j++)
    {
      gchar *test_name = g_strdup_printf ("/nc/data_cluster_wl_monte_carlo/%s/%s/%s",
                                          tests_obj[i].shape_name,
                                          tests_obj[i].redshift_name,
                                          tests_obj[i].ell_conv_name);

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
  gdouble min_radius                  = ncm_rng_uniform_gen (rng, 0.1, 0.5);
  gdouble max_radius                  = ncm_rng_uniform_gen (rng, 2.0, 5.0);
  gdouble ra                          = ncm_rng_uniform_gen (rng, -180, 180);
  gdouble dec                         = ncm_rng_uniform_gen (rng, -90, 90);
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDPosition *p_dist;
  NcGalaxySDObsRedshift *z_dist;
  NcGalaxySDShape *s_dist;
  NcGalaxyWLObsEllipConv ell_conv;

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
  test->ell_coord = NC_GALAXY_WL_OBS_COORD_CELESTIAL;

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
  guint nrows                       = 200;
  guint npoints                     = 20;
  gdouble z_min                     = 0.01;
  gdouble z_max                     = 5.0;
  NcGalaxyWLObs *obs;
  GStrv columns_strv;
  guint i;

  // if (NC_IS_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift))
  //   nrows = 200;
  // else if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
  //   nrows = 200;

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
test_nc_data_cluster_wl_monte_carlo (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmData *data       = NCM_DATA (test->dcwl);
  NcmDataset *dataset = ncm_dataset_new_array (&data, 1);
  NcmLikelihood *like = ncm_likelihood_new (dataset);
  NcmFit *fit         = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", like, test->mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmStatsVec *stats  = ncm_stats_vec_new (3, NCM_STATS_VEC_COV, FALSE);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint nfits         = 10;
  guint nruns         = 1;
  guint i, j;

  nc_data_cluster_wl_set_resample_flag (test->dcwl, NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL);
  g_object_set (test->dcwl, "enable-parallel", TRUE, NULL);

  /* Test mass fits */
  ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_mset_param_set_ftype (test->mset, nc_halo_position_id (), NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (test->mset, nc_halo_position_id (), NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FIXED);

  for (i = 0; i < nruns; i++)
  {
    gdouble log10M     = ncm_rng_uniform_gen (rng, 14.0, 16.0);
    gdouble ra         = ncm_rng_uniform_gen (rng, -180.0, 180.0);
    gdouble dec        = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    gdouble min_radius = ncm_rng_uniform_gen (rng, 0.1, 0.7);
    gdouble max_radius = ncm_rng_uniform_gen (rng, 1.0, 5.0);

    /* Set fiducial variables */
    ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", log10M, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", ra, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", dec, NULL);

    nc_galaxy_sd_position_set_ra_lim (test->galaxy_position, ra - 0.2 / cos (dec * M_PI / 180.0), ra + 0.2 / cos (dec * M_PI / 180.0));
    nc_galaxy_sd_position_set_dec_lim (test->galaxy_position, dec - 0.2, dec + 0.2);

    nc_data_cluster_wl_set_cut (test->dcwl, min_radius, max_radius);

    j = 0;

    while (j < nfits)
    {
      /* reset parameters before resample */
      ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", log10M, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", ra, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", dec, NULL);

      ncm_data_resample (data, test->mset, rng);
      ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

      gdouble param_fit = ncm_mset_fparam_get (test->mset, 0);

      if ((param_fit > 10.0) && (param_fit < 17.0))
      {
        ncm_stats_vec_set (stats, 0, ncm_mset_fparam_get (test->mset, 0));

        ncm_stats_vec_update (stats);
        j++;
      }
    }

    const gdouble mean_log10M = ncm_stats_vec_get_mean (stats, 0);
    const gdouble sd_log10M   = ncm_stats_vec_get_sd (stats, 0);

    ncm_assert_cmpdouble (mean_log10M, >, log10M - 8.0 * sd_log10M / sqrt (nfits));
    ncm_assert_cmpdouble (mean_log10M, <, log10M + 8.0 * sd_log10M / sqrt (nfits));
  }

  /* Test ra fits */
  ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (test->mset, nc_halo_position_id (), NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  ncm_mset_param_set_ftype (test->mset, nc_halo_position_id (), NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FIXED);

  for (i = 0; i < nruns; i++)
  {
    gdouble log10M     = ncm_rng_uniform_gen (rng, 14.0, 16.0);
    gdouble ra         = ncm_rng_uniform_gen (rng, -180.0, 180.0);
    gdouble dec        = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    gdouble min_radius = ncm_rng_uniform_gen (rng, 0.1, 0.7);
    gdouble max_radius = ncm_rng_uniform_gen (rng, 1.0, 5.0);

    /* Set fiducial variables */
    ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", log10M, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", ra, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", dec, NULL);

    ncm_model_param_set_lower_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_RA, -180.0);
    ncm_model_param_set_upper_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_RA, 180.0);

    ncm_model_param_set_lower_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_RA, ra - 0.008 / cos (dec * M_PI / 180.0));
    ncm_model_param_set_upper_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_RA, ra + 0.008 / cos (dec * M_PI / 180.0));

    nc_galaxy_sd_position_set_ra_lim (test->galaxy_position, ra - 0.2 / cos (dec * M_PI / 180.0), ra + 0.2 / cos (dec * M_PI / 180.0));
    nc_galaxy_sd_position_set_dec_lim (test->galaxy_position, dec - 0.2, dec + 0.2);

    nc_data_cluster_wl_set_cut (test->dcwl, min_radius, max_radius);

    for (j = 0; j < nfits; j++)
    {
      /* reset parameters before resample */
      ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", log10M, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", ra, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", dec, NULL);

      ncm_data_resample (data, test->mset, rng);
      ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

      ncm_stats_vec_set (stats, 0, ncm_mset_fparam_get (test->mset, 0));

      ncm_stats_vec_update (stats);
    }

    const gdouble mean_RA = ncm_stats_vec_get_mean (stats, 0);
    const gdouble sd_RA   = ncm_stats_vec_get_sd (stats, 0);

    ncm_assert_cmpdouble (mean_RA, >, ra - 6.0 * sd_RA / sqrt (nfits));
    ncm_assert_cmpdouble (mean_RA, <, ra + 6.0 * sd_RA / sqrt (nfits));
  }

  /* Test dec fits */
  ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), NC_HALO_CM_PARAM_LOG10M_DELTA, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (test->mset, nc_halo_position_id (), NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (test->mset, nc_halo_position_id (), NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);

  for (i = 0; i < nruns; i++)
  {
    gdouble log10M     = ncm_rng_uniform_gen (rng, 14.0, 16.0);
    gdouble ra         = ncm_rng_uniform_gen (rng, -180.0, 180.0);
    gdouble dec        = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    gdouble min_radius = ncm_rng_uniform_gen (rng, 0.1, 0.7);
    gdouble max_radius = ncm_rng_uniform_gen (rng, 1.0, 5.0);

    /* Set fiducial variables */
    ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", log10M, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", ra, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", dec, NULL);

    ncm_model_param_set_lower_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_DEC, -90.0);
    ncm_model_param_set_upper_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_DEC, 90.0);

    ncm_model_param_set_lower_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_DEC, dec - 0.008);
    ncm_model_param_set_upper_bound (NCM_MODEL (test->hp), NC_HALO_POSITION_DEC, dec + 0.008);

    nc_galaxy_sd_position_set_ra_lim (test->galaxy_position, ra - 0.2 / cos (dec * M_PI / 180.0), ra + 0.2 / cos (dec * M_PI / 180.0));
    nc_galaxy_sd_position_set_dec_lim (test->galaxy_position, dec - 0.2, dec + 0.2);

    nc_data_cluster_wl_set_cut (test->dcwl, min_radius, max_radius);

    for (j = 0; j < nfits; j++)
    {
      /* reset parameters before resample */
      ncm_model_param_set_by_name (NCM_MODEL (test->hms), "log10MDelta", log10M, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "ra", ra, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->hp), "dec", dec, NULL);

      ncm_data_resample (data, test->mset, rng);
      ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

      ncm_stats_vec_set (stats, 0, ncm_mset_fparam_get (test->mset, 0));

      ncm_stats_vec_update (stats);
    }

    const gdouble mean_DEC = ncm_stats_vec_get_mean (stats, 0);
    const gdouble sd_DEC   = ncm_stats_vec_get_sd (stats, 0);

    ncm_assert_cmpdouble (mean_DEC, >, dec - 6.0 * sd_DEC / sqrt (nfits));
    ncm_assert_cmpdouble (mean_DEC, <, dec + 6.0 * sd_DEC / sqrt (nfits));
  }

  ncm_dataset_clear (&dataset);
  ncm_likelihood_clear (&like);
  ncm_fit_clear (&fit);
  ncm_stats_vec_clear (&stats);
  ncm_rng_clear (&rng);
}


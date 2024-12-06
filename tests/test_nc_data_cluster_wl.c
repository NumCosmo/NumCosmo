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
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcGalaxySDObsRedshift *galaxy_redshift;
  NcGalaxySDPosition *galaxy_position;
  NcGalaxySDShape *galaxy_shape;
  NcmMSet *mset;
} TestNcDataClusterWL;


static void test_nc_data_cluster_wl_new_spec (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_new_gauss (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_new_pz (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_gen (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_ra_dec_limits (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_gen_obs (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/data_cluster_wl/spec/gen", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_gen,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/spec/ra_dec/limits", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_ra_dec_limits,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/spec/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_gen_obs,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/spec/m2lnP", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_m2lnP,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/spec/serialize", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_serialize,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/gauss/ra_dec/limits", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_ra_dec_limits,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/gauss/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_gen_obs,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/gauss/m2lnP", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_m2lnP,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/gauss/serialize", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_serialize,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/pz/ra_dec/limits", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_pz,
              &test_nc_data_cluster_wl_ra_dec_limits,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/pz/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_pz,
              &test_nc_data_cluster_wl_gen_obs,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/pz/m2lnP", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_pz,
              &test_nc_data_cluster_wl_m2lnP,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/pz/serialize", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_pz,
              &test_nc_data_cluster_wl_serialize,
              &test_nc_data_cluster_wl_free);

  g_test_run ();

  return 0;
}

static void
test_nc_data_cluster_wl_new_spec (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (z_true_dist));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloMassSummary *hms              = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->dcwl  = dcwl;
  test->cosmo = cosmo;

  test->hms                  = hms;
  test->density_profile      = dp;
  test->halo_position        = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_distance_free (dist);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_new_gauss (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_gauss_new (z_true_dist));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloMassSummary *hms              = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->dcwl  = dcwl;
  test->cosmo = cosmo;

  test->hms                  = hms;
  test->density_profile      = dp;
  test->halo_position        = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  nc_galaxy_sd_true_redshift_free (z_true_dist);
  nc_distance_free (dist);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_new_pz (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_pz_new ());
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloMassSummary *hms              = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->dcwl  = dcwl;
  test->cosmo = cosmo;

  test->hms                  = hms;
  test->density_profile      = dp;
  test->halo_position        = hp;
  test->surface_mass_density = smd;

  test->galaxy_redshift = z_dist;
  test->galaxy_position = p_dist;
  test->galaxy_shape    = s_dist;

  test->mset = ncm_mset_new (cosmo, NULL, dp, hp, smd, z_dist, p_dist, s_dist, NULL);

  nc_distance_free (dist);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata)
{
  nc_hicosmo_clear (&test->cosmo);
  nc_halo_density_profile_clear (&test->density_profile);
  nc_halo_mass_summary_clear (&test->hms);
  nc_halo_position_clear (&test->halo_position);
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
  NcGalaxyWLObs *obs;
  GStrv columns_strv;
  guint i;

  while (l)
  {
    g_strv_builder_add (builder, l->data);
    l = g_list_next (l);
  }

  columns_strv = g_strv_builder_end (builder);
  g_list_free_full (columns, g_free);
  obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, nrows, columns_strv);

  if (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift))
  {
    for (i = 0; i < nrows; i++)
    {
      nc_galaxy_sd_obs_redshift_gauss_gen (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->galaxy_redshift), test->mset, z_data, 0.03, rng);
      nc_galaxy_sd_position_flat_gen (NC_GALAXY_SD_POSITION_FLAT (test->galaxy_position), test->mset, p_data, rng);
      nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, 0.1, 0.1, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, rng);
      nc_galaxy_sd_shape_data_write_row (s_data, obs, i);
    }
  }
  else if (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (test->galaxy_redshift))
  {
    for (i = 0; i < nrows; i++)
    {
      nc_galaxy_sd_obs_redshift_spec_gen (NC_GALAXY_SD_OBS_REDSHIFT_SPEC (test->galaxy_redshift), test->mset, z_data, rng);
      nc_galaxy_sd_position_flat_gen (NC_GALAXY_SD_POSITION_FLAT (test->galaxy_position), test->mset, p_data, rng);
      nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, 0.1, 0.1, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, rng);
      nc_galaxy_sd_shape_data_write_row (s_data, obs, i);
    }
  }
  else if (NC_IS_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift))
  {
    guint nrows   = 100;
    guint npoints = 1000;
    gdouble z_min = 0.01;
    gdouble z_max = 10.0;
    gdouble z_avg;
    gdouble z_sd;

    obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, nrows, columns_strv);

    for (i = 0; i < nrows; i++)
    {
      NcmVector *xv = ncm_vector_new (npoints);
      NcmVector *yv = ncm_vector_new (npoints);
      NcmSpline *pz;
      guint j;

      z_avg = g_test_rand_double_range (z_min, z_max);
      z_sd  = 0.03 * (1.0 + z_avg);

      for (j = 0; j < npoints; j++)
      {
        gdouble x_min = z_avg - 5.0 * z_sd;

        if (x_min < z_min)
          x_min = z_min;

        gdouble x = x_min + 10.0 * z_sd * j / ((gdouble) npoints - 1.0);
        gdouble y = exp (-0.5 * gsl_pow_2 ((x - z_avg) / z_sd)) / (sqrt (2.0 * M_PI) * z_sd);

        ncm_vector_fast_set (xv, j, x);
        ncm_vector_fast_set (yv, j, y);
      }

      pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

      nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift), z_data, pz);
      nc_galaxy_sd_obs_redshift_pz_prepare (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift), z_data);

      nc_galaxy_sd_obs_redshift_pz_gen (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->galaxy_redshift), test->mset, z_data, rng);
      nc_galaxy_sd_position_flat_gen (NC_GALAXY_SD_POSITION_FLAT (test->galaxy_position), test->mset, p_data, rng);
      nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, 0.1, 0.1, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, rng);
      nc_galaxy_sd_shape_data_write_row (s_data, obs, i);

      ncm_vector_free (xv);
      ncm_vector_free (yv);
      ncm_spline_free (pz);
    }
  }
  else
  {
    g_error ("Unknown galaxy redshift distribution");
  }

  g_strfreev (columns_strv);

  nc_data_cluster_wl_set_obs (test->dcwl, obs);

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_wl_obs_free (obs);
  ncm_rng_free (rng);
  g_strv_builder_unref (builder);
}

static void
test_nc_data_cluster_wl_ra_dec_limits (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxyWLObs *obs;
  gdouble ra_min, ra_max, dec_min, dec_max;
  guint ngals;
  guint i;

  nc_galaxy_sd_position_get_ra_lim (test->galaxy_position, &ra_min, &ra_max);
  nc_galaxy_sd_position_get_dec_lim (test->galaxy_position, &dec_min, &dec_max);

  test_nc_data_cluster_wl_gen (test, pdata);

  obs   = nc_data_cluster_wl_peek_obs (test->dcwl);
  ngals = nc_galaxy_wl_obs_len (obs);

  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_POSITION_COL_RA));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_POSITION_COL_DEC));

  for (i = 0; i < ngals; i++)
  {
    gdouble ra  = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_RA, i);
    gdouble dec = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_DEC, i);
    gdouble theta, phi;

    nc_halo_position_polar_angles (test->halo_position, ra, dec, &theta, &phi);

    g_assert_cmpfloat (ra, >=, ra_min);
    g_assert_cmpfloat (ra, <=, ra_max);

    g_assert_cmpfloat (dec, >=, dec_min);
    g_assert_cmpfloat (dec, <=, dec_max);
  }
}

static void
test_nc_data_cluster_wl_gen_obs (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxyWLObs *obs;
  gdouble ngals;
  guint i;

  nc_data_cluster_wl_set_cut (test->dcwl, 0.3, 3.0);
  test_nc_data_cluster_wl_gen (test, pdata);

  obs   = nc_data_cluster_wl_peek_obs (test->dcwl);
  ngals = nc_galaxy_wl_obs_len (obs);

  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_OBS_REDSHIFT_COL_Z));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_POSITION_COL_RA));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_POSITION_COL_DEC));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS_1));
  g_assert_true (g_strv_contains ((const gchar * const *) nc_galaxy_wl_obs_peek_columns (obs), NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS_2));

  for (i = 0; i < ngals; i++)
  {
    const gdouble z             = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
    const gdouble ra            = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_RA, i);
    const gdouble dec           = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_DEC, i);
    const gdouble epsilon_int_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
    const gdouble epsilon_int_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);
    const gdouble sigma_obs_1   = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS_1, i);
    const gdouble sigma_obs_2   = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS_2, i);
    const gdouble var_obs_1     = sigma_obs_1 * sigma_obs_1;
    const gdouble var_obs_2     = sigma_obs_2 * sigma_obs_2;
    const gdouble sigma_int     = ncm_model_orig_param_get (NCM_MODEL (test->galaxy_shape), NC_GALAXY_SD_SHAPE_GAUSS_SIGMA_INT);
    const gdouble var_int       = sigma_int * sigma_int;
    const gdouble e_rms_1       = sqrt (var_int + var_obs_1);
    const gdouble e_rms_2       = sqrt (var_int + var_obs_2);
    gdouble theta, phi;

    nc_halo_position_polar_angles (test->halo_position, ra, dec, &theta, &phi);

    g_assert_cmpfloat (ra, >=, -0.2);
    g_assert_cmpfloat (ra, <=, 0.2);
    g_assert_cmpfloat (dec, >=, -0.2);
    g_assert_cmpfloat (dec, <=, 0.2);
    g_assert_cmpfloat (z, >=, 0.0);
    g_assert_cmpfloat (z, <=, 20.0);
    g_assert_cmpfloat (epsilon_int_1, >=, -5.0 * e_rms_1);
    g_assert_cmpfloat (epsilon_int_1, <=, 5.0 * e_rms_1);
    g_assert_cmpfloat (epsilon_int_2, >=, -5.0 * e_rms_2);
    g_assert_cmpfloat (epsilon_int_2, <=, 5.0 * e_rms_2);
    g_assert_cmpfloat (sigma_int, ==, 0.3);
    g_assert_cmpfloat (sigma_obs_1, ==, 0.1);
    g_assert_cmpfloat (sigma_obs_2, ==, 0.1);
  }

  ncm_rng_free (rng);
}

static void
test_nc_data_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata)
{
  gdouble m2lnL_a, m2lnL_b;

  test_nc_data_cluster_wl_gen (test, pdata);

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

    ncm_mset_param_set (test->mset, nc_halo_position_id (), NC_HALO_POSITION_DEC, 0.1);

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
}

static void
test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

  test_nc_data_cluster_wl_gen (test, pdata);

  {
    NcmMSet *mset_dup = ncm_mset_dup (test->mset, ser);
    NcmData *dcwl_dup = ncm_data_dup (NCM_DATA (test->dcwl), ser);
    gdouble m2lnL, m2lnL_dup;

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

    ncm_mset_free (mset_dup);
    ncm_data_free (dcwl_dup);
  }

  ncm_serialize_free (ser);
}


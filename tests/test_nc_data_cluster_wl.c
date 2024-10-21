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
static void test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_gen (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_ra_dec_limits (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_gen_obs (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_gata_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_m2lnP_integ (TestNcDataClusterWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/gen", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_gen,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/ra_dec/limits", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_ra_dec_limits,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_data_cluster_wl_gen_obs,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/spec/m2lnP", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_spec,
              &test_nc_gata_cluster_wl_m2lnP,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/r_cut", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_ra_dec_limits,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/numcosmo/data/nc_data_cluster_wl/gauss/gen_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_gauss,
              &test_nc_data_cluster_wl_gen_obs,
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
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->dcwl  = dcwl;
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

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_new_gauss (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDShape *s_dist             = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_gauss_new ());
  NcGalaxySDTrueRedshift *z_true_dist = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (0.0, 5.0));
  NcGalaxySDObsRedshift *z_dist       = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_gauss_new (z_true_dist));
  NcGalaxySDPosition *p_dist          = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (-0.2, 0.2, -0.2, 0.2));
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new ();
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                    = nc_distance_new (100.0);
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcHaloPosition *hp                  = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd         = nc_wl_surface_mass_density_new (dist);

  nc_halo_position_prepare (hp, cosmo);

  test->dcwl  = dcwl;
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

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata)
{
  nc_hicosmo_clear (&test->cosmo);
  nc_halo_density_profile_clear (&test->density_profile);
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
  else
  {
    g_error ("Unknown galaxy redshift distribution");
  }

  g_strfreev (columns_strv);

  nc_data_cluster_wl_set_obs (test->dcwl, obs);

  nc_galaxy_sd_obs_redshift_data_unref (z_data);
  nc_galaxy_sd_position_data_unref (p_data);
  nc_galaxy_sd_shape_data_unref (s_data);
  nc_galaxy_wl_obs_clear (&obs);
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
    g_assert_cmpfloat (z, <=, 5.0);
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
test_nc_gata_cluster_wl_m2lnP (TestNcDataClusterWL *test, gconstpointer pdata)
{
}

static void
test_nc_data_cluster_wl_m2lnP_integ (TestNcDataClusterWL *test, gconstpointer pdata)
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


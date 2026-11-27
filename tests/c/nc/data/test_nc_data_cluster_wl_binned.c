/***************************************************************************
 *            test_nc_data_cluster_wl_binned.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * test_nc_data_cluster_wl_binned.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 *
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
#include <complex.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_Z_CL (0.3)
#define TEST_N_GAL (400)

static const NcmCatalogColType test_double_cols[64] = {NCM_CATALOG_COL_TYPE_DOUBLE};

typedef struct _TestNcDataClusterWLBinned
{
  NcmMSet *mset;
  NcGalaxyWLObs *obs;
  GArray *edges;
} TestNcDataClusterWLBinned;

static void
test_nc_data_cluster_wl_binned_new (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcHICosmo *cosmo                 = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                 = nc_distance_new (100.0);
  NcHaloMassSummary *hms           = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp         = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp               = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd      = nc_wl_surface_mass_density_new (dist);
  NcGalaxyWLBinnedScatterGauss *sg = nc_galaxy_wl_binned_scatter_gauss_new ();
  NcmRNG *rng                      = ncm_rng_seeded_new (NULL, 42);
  const gchar *cols[]              = {"ra", "dec", "z", "epsilon_obs_1", "epsilon_obs_2", "std_noise", "e_rms", NULL};
  const gdouble edges[]            = {0.3, 0.6, 1.2, 2.4};
  guint i;

  ncm_model_param_set_by_name (NCM_MODEL (hms), "log10MDelta", 14.5, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (hp), "z", TEST_Z_CL, NULL);
  nc_halo_position_prepare (hp, cosmo);
  nc_wl_surface_mass_density_prepare (smd, cosmo);

  test->mset = ncm_mset_new (NCM_MODEL (cosmo), NULL, NCM_MODEL (dp), NCM_MODEL (hp), NCM_MODEL (smd), NCM_MODEL (sg), NULL);
  test->obs  = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_WL_ELLIPTICITY_FRAME_CELESTIAL, TEST_N_GAL, (GStrv) cols, test_double_cols, g_strv_length ((GStrv) cols));

  /* Noiseless galaxies: epsilon_obs is the reduced shear */
  for (i = 0; i < TEST_N_GAL; i++)
  {
    const gdouble ra  = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    const gdouble dec = ncm_rng_uniform_gen (rng, -0.2, 0.2);
    const gdouble z   = ncm_rng_uniform_gen (rng, 0.1, 1.5);
    gdouble theta, phi, radius, g_t = 0.0;
    complex double g;

    nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);
    radius = nc_halo_position_projected_radius (hp, cosmo, theta);

    if (z > TEST_Z_CL)
      g_t = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, radius, z, TEST_Z_CL, TEST_Z_CL);

    g = g_t * cexp (2.0 * I * phi);

    ncm_catalog_set (NCM_CATALOG (test->obs), "ra", i, ra, NULL);
    ncm_catalog_set (NCM_CATALOG (test->obs), "dec", i, dec, NULL);
    ncm_catalog_set (NCM_CATALOG (test->obs), "z", i, z, NULL);
    ncm_catalog_set (NCM_CATALOG (test->obs), "epsilon_obs_1", i, creal (g), NULL);
    ncm_catalog_set (NCM_CATALOG (test->obs), "epsilon_obs_2", i, cimag (g), NULL);
    ncm_catalog_set (NCM_CATALOG (test->obs), "std_noise", i, 0.1, NULL);
    ncm_catalog_set (NCM_CATALOG (test->obs), "e_rms", i, 0.35, NULL);
  }

  test->edges = g_array_new (FALSE, FALSE, sizeof (gdouble));
  g_array_append_vals (test->edges, edges, 4);

  ncm_rng_free (rng);
  nc_galaxy_wl_binned_scatter_gauss_free (sg);
  nc_wl_surface_mass_density_free (smd);
  nc_halo_position_free (hp);
  nc_halo_density_profile_free (dp);
  nc_halo_mass_summary_free (hms);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

static void
test_nc_data_cluster_wl_binned_free (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  g_array_unref (test->edges);
  nc_galaxy_wl_obs_free (test->obs);
  ncm_mset_free (test->mset);
}

static NcGalaxyWLBinnedEstimator *
_test_estimator (gconstpointer pdata)
{
  if (GPOINTER_TO_INT (pdata) == 0)
    return NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());
  else
    return NC_GALAXY_WL_BINNED_ESTIMATOR (g_object_new (NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR_GT,
                                                        "error-type", NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_SHAPE_NOISE,
                                                        NULL));
}

static void
test_nc_data_cluster_wl_binned_exact (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = _test_estimator (pdata);
  NcDataClusterWLBinned *dcwlb   = nc_data_cluster_wl_binned_new_from_obs (est, test->mset, test->obs, test->edges);
  NcGalaxyWLObsBinned *binned    = nc_data_cluster_wl_binned_peek_binned (dcwlb);
  const guint n_bins             = nc_galaxy_wl_obs_binned_len (binned);
  NcmVector *model               = ncm_vector_new (n_bins);
  gdouble m2lnL;
  guint b;

  g_assert_cmpuint (n_bins, ==, 3);
  g_assert_cmpuint (ncm_data_get_length (NCM_DATA (dcwlb)), ==, 3);
  g_assert_true (nc_galaxy_wl_obs_binned_has_members (binned));
  g_assert_true (nc_galaxy_wl_obs_binned_peek_source (binned) == test->obs);
  g_assert_cmpstr (nc_galaxy_wl_obs_binned_peek_estimator (binned), ==, nc_galaxy_wl_binned_estimator_get_name (est));
  ncm_assert_cmpdouble_e (nc_galaxy_wl_obs_binned_get_z_cluster (binned), ==, TEST_Z_CL, 1.0e-15, 0.0);

  nc_data_cluster_wl_binned_set_exact (dcwlb, TRUE);
  g_assert_true (nc_data_cluster_wl_binned_get_exact (dcwlb));
  nc_data_cluster_wl_binned_eval_model (dcwlb, test->mset, model);

  for (b = 0; b < n_bins; b++)
    ncm_assert_cmpdouble_e (ncm_vector_get (model, b), ==, nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, b), 1.0e-10, 0.0);

  ncm_data_m2lnL_val (NCM_DATA (dcwlb), test->mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, <, 1.0e-12);

  nc_data_cluster_wl_binned_set_exact (dcwlb, FALSE);
  ncm_data_m2lnL_val (NCM_DATA (dcwlb), test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));
  g_assert_cmpfloat (m2lnL, >, 0.0);

  {
    NcmSerialize *ser                = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcDataClusterWLBinned *dcwlb_dup = NC_DATA_CLUSTER_WL_BINNED (ncm_serialize_dup_obj (ser, G_OBJECT (dcwlb)));
    gdouble m2lnL_dup;

    ncm_data_m2lnL_val (NCM_DATA (dcwlb_dup), test->mset, &m2lnL_dup);
    ncm_assert_cmpdouble_e (m2lnL_dup, ==, m2lnL, 1.0e-14, 0.0);

    nc_data_cluster_wl_binned_free (dcwlb_dup);
    ncm_serialize_free (ser);
  }

  ncm_vector_free (model);
  nc_data_cluster_wl_binned_free (dcwlb);
  nc_galaxy_wl_binned_estimator_free (est);
}

static void
test_nc_data_cluster_wl_binned_resample (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = _test_estimator (pdata);
  NcDataClusterWLBinned *dcwlb   = nc_data_cluster_wl_binned_new_from_obs (est, test->mset, test->obs, test->edges);
  NcGalaxyWLObsBinned *binned    = nc_data_cluster_wl_binned_peek_binned (dcwlb);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 7);
  NcmVector *before              = nc_galaxy_wl_obs_binned_get_column (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE);
  NcmVector *after;
  gdouble m2lnL;

  ncm_data_resample (NCM_DATA (dcwlb), test->mset, rng);
  after = nc_galaxy_wl_obs_binned_get_column (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE);

  g_assert_cmpfloat (fabs (ncm_vector_get (after, 0) - ncm_vector_get (before, 0)), >, 0.0);

  ncm_data_m2lnL_val (NCM_DATA (dcwlb), test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));

  ncm_vector_free (before);
  ncm_vector_free (after);
  ncm_rng_free (rng);
  nc_data_cluster_wl_binned_free (dcwlb);
  nc_galaxy_wl_binned_estimator_free (est);
}

static void
test_nc_data_cluster_wl_binned_properties (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds = nc_galaxy_wl_binned_estimator_delta_sigma_new ();
  NcGalaxyWLBinnedEstimatorGt *est_gt         = nc_galaxy_wl_binned_estimator_gt_new ();
  NcGalaxyWLBinnedEstimator *est              = NC_GALAXY_WL_BINNED_ESTIMATOR (est_ds);
  NcGalaxyWLBinnedScatterGauss *sg            = nc_galaxy_wl_binned_scatter_gauss_new ();
  gboolean subtract, include_lnnorm, beta_weighted;
  gdouble e_rms, sigma_pop;
  gdouble pz_reltol;
  guint exact_intervals;
  gint error_type;

  g_object_set (est_ds, "subtract-additive-bias", FALSE, "e-rms", 0.3, "pz-reltol", 1.0e-4, "exact-intervals", 10, NULL);
  g_object_get (est_ds, "subtract-additive-bias", &subtract, "e-rms", &e_rms, "pz-reltol", &pz_reltol, "exact-intervals", &exact_intervals, NULL);
  g_assert_false (subtract);
  g_assert_cmpfloat (e_rms, ==, 0.3);
  g_assert_cmpfloat (pz_reltol, ==, 1.0e-4);
  g_assert_cmpuint (exact_intervals, ==, 10);
  g_assert_false (nc_galaxy_wl_binned_estimator_delta_sigma_get_subtract_additive_bias (est_ds));
  g_assert_cmpfloat (nc_galaxy_wl_binned_estimator_get_e_rms (est), ==, 0.3);
  g_assert_cmpfloat (nc_galaxy_wl_binned_estimator_get_pz_reltol (est), ==, 1.0e-4);
  g_assert_cmpuint (nc_galaxy_wl_binned_estimator_get_exact_intervals (est), ==, 10);

  g_object_set (est_gt, "sigma-pop", 0.2, "error-type", NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_SHAPE_NOISE, "beta-weighted", TRUE, NULL);
  g_object_get (est_gt, "sigma-pop", &sigma_pop, "error-type", &error_type, "beta-weighted", &beta_weighted, NULL);
  g_assert_true (beta_weighted);
  g_assert_true (nc_galaxy_wl_binned_estimator_gt_get_beta_weighted (est_gt));
  g_assert_cmpfloat (sigma_pop, ==, 0.2);
  g_assert_cmpint (error_type, ==, NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_SHAPE_NOISE);

  g_object_set (sg, "include-lnnorm", TRUE, NULL);
  g_object_get (sg, "include-lnnorm", &include_lnnorm, NULL);
  g_assert_true (include_lnnorm);

  {
    GList *cols = nc_galaxy_wl_binned_estimator_required_columns (est);

    g_assert_cmpuint (g_list_length (cols), ==, 5);
    g_list_free_full (cols, g_free);
  }

  {
    NcGalaxyWLObsBinned *binned  = nc_galaxy_wl_binned_estimator_build (est, test->mset, test->obs, test->edges);
    NcDataClusterWLBinned *dcwlb = g_object_new (NC_TYPE_DATA_CLUSTER_WL_BINNED,
                                                 "estimator", est,
                                                 "binned", binned,
                                                 "exact", TRUE,
                                                 NULL);
    NcGalaxyWLObsBinned *binned_get;
    NcGalaxyWLBinnedEstimator *est_get;
    gboolean exact;

    g_object_get (dcwlb, "binned", &binned_get, "estimator", &est_get, "exact", &exact, NULL);
    g_assert_true (binned_get == binned);
    g_assert_true (est_get == est);
    g_assert_true (exact);
    g_assert_true (nc_data_cluster_wl_binned_peek_estimator (dcwlb) == est);

    {
      NcGalaxyWLObs *source;
      NcmCatalog *members;
      NcmMatrix *cov;
      gdouble z_cluster;
      gchar *name;

      g_object_get (binned, "source", &source, "members", &members, "cov", &cov, "z-cluster", &z_cluster, "estimator", &name, NULL);
      g_assert_true (source == test->obs);
      g_assert_true (members == nc_galaxy_wl_obs_binned_peek_members (binned));
      g_assert_true (cov == nc_galaxy_wl_obs_binned_peek_cov (binned));
      g_assert_cmpfloat (z_cluster, ==, TEST_Z_CL);
      g_assert_cmpstr (name, ==, "delta_sigma");

      nc_galaxy_wl_obs_free (source);
      ncm_catalog_free (members);
      ncm_matrix_free (cov);
      g_free (name);
    }

    nc_galaxy_wl_obs_binned_set_members (binned, NULL, NULL);
    g_assert_false (nc_galaxy_wl_obs_binned_has_members (binned));

    nc_galaxy_wl_obs_binned_free (binned_get);
    nc_galaxy_wl_binned_estimator_free (est_get);
    nc_data_cluster_wl_binned_free (dcwlb);
    nc_galaxy_wl_obs_binned_free (binned);
  }

  nc_galaxy_wl_binned_scatter_gauss_free (sg);
  nc_galaxy_wl_binned_estimator_gt_free (est_gt);
  nc_galaxy_wl_binned_estimator_delta_sigma_free (est_ds);
}

static NcGalaxyWLObsBinned *
_test_external_profile (void)
{
  const gdouble r[]    = {0.5, 1.0, 2.0};
  const gdouble v[]    = {1.0e14, 5.0e13, 2.0e13};
  GArray *radius       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *value        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmMatrix *cov       = ncm_matrix_new (3, 3);
  const gchar *extra[] = {NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT, NULL};
  NcGalaxyWLObsBinned *binned;
  guint b;

  ncm_matrix_set_identity (cov);
  ncm_matrix_scale (cov, 1.0e26);

  g_array_append_vals (radius, r, 3);
  g_array_append_vals (value, v, 3);

  binned = nc_galaxy_wl_obs_binned_new_from_arrays (radius, value, cov, (GStrv) extra);

  for (b = 0; b < 3; b++)
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT, b, 2.0e-16);

  g_array_unref (radius);
  g_array_unref (value);
  ncm_matrix_free (cov);

  return binned;
}

static void
test_nc_data_cluster_wl_binned_external (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());
  NcGalaxyWLObsBinned *binned    = _test_external_profile ();
  NcDataClusterWLBinned *dcwlb   = nc_data_cluster_wl_binned_new (est, binned);
  gdouble m2lnL;

  g_assert_false (nc_galaxy_wl_obs_binned_has_members (binned));
  g_assert_null (nc_galaxy_wl_obs_binned_peek_source (binned));
  g_assert_true (gsl_isnan (nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_R_LO, 0)));

  ncm_data_m2lnL_val (NCM_DATA (dcwlb), test->mset, &m2lnL);
  g_assert_true (gsl_finite (m2lnL));

  nc_data_cluster_wl_binned_free (dcwlb);
  nc_galaxy_wl_obs_binned_free (binned);
  nc_galaxy_wl_binned_estimator_free (est);
}

static void
test_nc_data_cluster_wl_binned_pitfalls (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/nc/data_cluster_wl_binned/exact_external/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exact model needs a profile built from a per-galaxy catalog*");

  g_test_trap_subprocess ("/nc/data_cluster_wl_binned/edges/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*bin edges must increase*");

  g_test_trap_subprocess ("/nc/data_cluster_wl_binned/no_scatter/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*must hold a NcGalaxyWLBinnedScatter*");

  g_test_trap_subprocess ("/nc/data_cluster_wl_binned/no_e_rms/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*has no column `e_rms'*");

  g_test_trap_subprocess ("/nc/data_cluster_wl_binned/exact_no_members/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the profile has no members*");
}

static void
test_nc_data_cluster_wl_binned_exact_external_subprocess (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());
  NcGalaxyWLObsBinned *binned    = _test_external_profile ();
  NcDataClusterWLBinned *dcwlb   = nc_data_cluster_wl_binned_new (est, binned);

  nc_data_cluster_wl_binned_set_exact (dcwlb, TRUE);
}

static void
test_nc_data_cluster_wl_binned_edges_subprocess (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());

  g_array_index (test->edges, gdouble, 1) = 0.1;
  nc_galaxy_wl_binned_estimator_build (est, test->mset, test->obs, test->edges);
}

static void
test_nc_data_cluster_wl_binned_no_scatter_subprocess (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());
  NcDataClusterWLBinned *dcwlb   = nc_data_cluster_wl_binned_new_from_obs (est, test->mset, test->obs, test->edges);
  gdouble m2lnL;

  ncm_mset_remove (test->mset, nc_galaxy_wl_binned_scatter_id ());
  ncm_data_m2lnL_val (NCM_DATA (dcwlb), test->mset, &m2lnL);
}

static void
test_nc_data_cluster_wl_binned_no_e_rms_subprocess (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  const gchar *cols[]            = {"ra", "dec", "z", "epsilon_obs_1", "epsilon_obs_2", "std_noise", NULL};
  NcGalaxyWLObs *obs             = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_WL_ELLIPTICITY_FRAME_CELESTIAL, 2, (GStrv) cols, test_double_cols, g_strv_length ((GStrv) cols));
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());

  nc_galaxy_wl_binned_estimator_build (est, test->mset, obs, test->edges);
}

static void
test_nc_data_cluster_wl_binned_exact_no_members_subprocess (TestNcDataClusterWLBinned *test, gconstpointer pdata)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_delta_sigma_new ());
  NcGalaxyWLObsBinned *binned    = _test_external_profile ();
  NcmVector *model               = ncm_vector_new (3);

  nc_galaxy_wl_binned_estimator_eval_model_exact (est, test->mset, binned, model);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/data_cluster_wl_binned/delta_sigma/exact", TestNcDataClusterWLBinned, GINT_TO_POINTER (0),
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_exact, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/gt/exact", TestNcDataClusterWLBinned, GINT_TO_POINTER (1),
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_exact, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/delta_sigma/resample", TestNcDataClusterWLBinned, GINT_TO_POINTER (0),
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_resample, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/gt/resample", TestNcDataClusterWLBinned, GINT_TO_POINTER (1),
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_resample, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/properties", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_properties, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/external", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_external, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/pitfalls", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_pitfalls, &test_nc_data_cluster_wl_binned_free);

  g_test_add ("/nc/data_cluster_wl_binned/exact_external/subprocess", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_exact_external_subprocess, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/edges/subprocess", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_edges_subprocess, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/no_scatter/subprocess", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_no_scatter_subprocess, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/no_e_rms/subprocess", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_no_e_rms_subprocess, &test_nc_data_cluster_wl_binned_free);
  g_test_add ("/nc/data_cluster_wl_binned/exact_no_members/subprocess", TestNcDataClusterWLBinned, NULL,
              &test_nc_data_cluster_wl_binned_new, &test_nc_data_cluster_wl_binned_exact_no_members_subprocess, &test_nc_data_cluster_wl_binned_free);

  g_test_run ();
}


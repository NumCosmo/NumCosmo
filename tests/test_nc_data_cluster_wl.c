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


static void test_nc_data_cluster_wl_new_flat (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_new_lsst_srd (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata);

static void test_nc_data_cluster_wl_fit (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_kde_cmp (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_set_obs (TestNcDataClusterWL *test, gconstpointer pdata);
static void test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/data_cluster_wl/flat/fit", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_flat,
              &test_nc_data_cluster_wl_fit,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/flat/kde_cmp", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_flat,
              &test_nc_data_cluster_wl_kde_cmp,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/flat/set_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_flat,
              &test_nc_data_cluster_wl_set_obs,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/flat/serialize", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_flat,
              &test_nc_data_cluster_wl_serialize,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/lsst_srd/fit", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_lsst_srd,
              &test_nc_data_cluster_wl_fit,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/lsst_srd/kde_cmp", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_lsst_srd,
              &test_nc_data_cluster_wl_kde_cmp,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/lsst_srd/set_obs", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_lsst_srd,
              &test_nc_data_cluster_wl_set_obs,
              &test_nc_data_cluster_wl_free);

  g_test_add ("/nc/data_cluster_wl/lsst_srd/serialize", TestNcDataClusterWL, NULL,
              &test_nc_data_cluster_wl_new_lsst_srd,
              &test_nc_data_cluster_wl_serialize,
              &test_nc_data_cluster_wl_free);

  g_test_run ();

  return 0;
}

static void
test_nc_data_cluster_wl_new_flat (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDPositionFlat *gsdpf = nc_galaxy_sd_position_flat_new (0.0, 2.0, 0.4, 1.2);
  NcGalaxySDZProxyGauss *gsdzpg = nc_galaxy_sd_z_proxy_gauss_new (0.4, 1.0, 0.03);
  NcGalaxySDShapeGauss *gss     = nc_galaxy_sd_shape_gauss_new ();
  NcDataClusterWL *dcwl         = nc_data_cluster_wl_new (NC_GALAXY_SD_SHAPE (gss),
                                                          NC_GALAXY_SD_Z_PROXY (gsdzpg),
                                                          NC_GALAXY_SD_POSITION (gsdpf),
                                                          0.4);

  nc_galaxy_sd_shape_gauss_set_sigma (gss, 0.001);
  nc_data_cluster_wl_set_cut (dcwl, 0.0, 2.0);

  test->dcwl = dcwl;

  nc_galaxy_sd_position_flat_free (gsdpf);
  nc_galaxy_sd_z_proxy_gauss_free (gsdzpg);
  nc_galaxy_sd_shape_gauss_free (gss);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_new_lsst_srd (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst = nc_galaxy_sd_position_lsst_srd_new (0.0, 2.0, 0.4, 1.2);
  NcGalaxySDZProxyGauss *gsdzpg       = nc_galaxy_sd_z_proxy_gauss_new (0.4, 1.0, 0.03);
  NcGalaxySDShapeGauss *gss           = nc_galaxy_sd_shape_gauss_new ();
  NcDataClusterWL *dcwl               = nc_data_cluster_wl_new (NC_GALAXY_SD_SHAPE (gss),
                                                                NC_GALAXY_SD_Z_PROXY (gsdzpg),
                                                                NC_GALAXY_SD_POSITION (gsdplsst),
                                                                0.4);

  nc_galaxy_sd_shape_gauss_set_sigma (gss, 0.001);
  nc_data_cluster_wl_set_cut (dcwl, 0.0, 2.0);

  test->dcwl = dcwl;

  nc_galaxy_sd_position_lsst_srd_free (gsdplsst);
  nc_galaxy_sd_z_proxy_gauss_free (gsdzpg);
  nc_galaxy_sd_shape_gauss_free (gss);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));
}

static void
test_nc_data_cluster_wl_free (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_data_cluster_wl_free, test->dcwl);
}

static void
test_nc_data_cluster_wl_fit (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  const guint ngals           = 200;
  NcmMatrix *gal_obs          = NULL;

  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0);

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, ngals, rng);

  {
    nc_data_cluster_wl_set_use_kde (test->dcwl, FALSE);
    ncm_data_set_init (NCM_DATA (test->dcwl), TRUE);

    {
      NcmDataset *dataset       = ncm_dataset_new_list (test->dcwl, NULL);
      NcmLikelihood *likelihood = ncm_likelihood_new (dataset);
      NcmMSet *mset             = ncm_mset_new (cosmo, dp, smd, NULL);
      NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_NLOPT, "ln-neldermead", likelihood, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL);

      ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-3, 0.0, NULL, NULL);
      ncm_fit_obs_fisher (fit);

      {
        const gdouble c_bf = ncm_model_param_get (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA);
        const gdouble m_bf = ncm_model_param_get (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA);
        const gdouble c_sd = ncm_fit_covar_sd (fit, nc_halo_density_profile_id (), NC_HALO_DENSITY_PROFILE_C_DELTA);
        const gdouble m_sd = ncm_fit_covar_sd (fit, nc_halo_density_profile_id (), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA);

        g_assert_cmpfloat (c_bf, <, 4.0 + 5.0 * c_sd);
        g_assert_cmpfloat (c_bf, >, 4.0 - 5.0 * c_sd);

        g_assert_cmpfloat (m_bf, <, 14.0 + 5.0 * m_sd);
        g_assert_cmpfloat (m_bf, >, 14.0 - 5.0 * m_sd);
      }

      ncm_dataset_free (dataset);
      ncm_likelihood_free (likelihood);
      ncm_mset_free (mset);
      ncm_fit_free (fit);
    }
  }

  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}

static void
test_nc_data_cluster_wl_kde_cmp (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist            = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);
  const guint ngals           = 1000;
  NcmStatsVec *m2lnP_stats    = ncm_stats_vec_new (3, NCM_STATS_VEC_COV, FALSE);
  NcmVector *m2lnP_int_gal    = ncm_vector_new (ngals);
  NcmVector *m2lnP_kde_gal    = ncm_vector_new (ngals);
  NcmMatrix *gal_obs          = NULL;
  NcmMatrix *obs;
  guint i;

  ncm_stats_vec_enable_quantile (m2lnP_stats, 0.5);

  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0);

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, ngals, rng);
  nc_data_cluster_wl_set_ndata (test->dcwl, 1000);

  obs = nc_data_cluster_wl_peek_obs (test->dcwl);

  nc_data_cluster_wl_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_int_gal);
  nc_data_cluster_wl_kde_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_kde_gal);

  for (i = 0; i < ngals; i++)
  {
    const gdouble m2lnP_int = ncm_vector_get (m2lnP_int_gal, i);
    const gdouble m2lnP_kde = ncm_vector_get (m2lnP_kde_gal, i);

    ncm_stats_vec_set (m2lnP_stats, 0, m2lnP_int);
    ncm_stats_vec_set (m2lnP_stats, 1, m2lnP_kde);
    ncm_stats_vec_set (m2lnP_stats, 2, m2lnP_int - m2lnP_kde);

    ncm_stats_vec_update (m2lnP_stats);
  }

  {
    const gdouble mean_diff = ncm_stats_vec_get_mean (m2lnP_stats, 2);

    ncm_stats_vec_reset (m2lnP_stats, TRUE);

    for (i = 0; i < ngals; i++)
    {
      const gdouble m2lnP_int = ncm_vector_get (m2lnP_int_gal, i) - mean_diff;
      const gdouble m2lnP_kde = ncm_vector_get (m2lnP_kde_gal, i);

      ncm_stats_vec_set (m2lnP_stats, 0, m2lnP_int);
      ncm_stats_vec_set (m2lnP_stats, 1, m2lnP_kde);
      ncm_stats_vec_set (m2lnP_stats, 2, fabs (m2lnP_int / m2lnP_kde - 1.0));

      ncm_stats_vec_update (m2lnP_stats);
    }
  }

  /*
   *  printf ("mean: %g\n", ncm_stats_vec_get_mean (m2lnP_stats, 2));
   *  printf ("sd: %g\n", ncm_stats_vec_get_sd (m2lnP_stats, 2));
   *  printf ("q50: %g\n", ncm_stats_vec_get_quantile (m2lnP_stats, 2));
   */
  g_assert_cmpfloat (ncm_stats_vec_get_quantile (m2lnP_stats, 2), <, 0.1);

  ncm_vector_free (m2lnP_int_gal);
  ncm_vector_free (m2lnP_kde_gal);
  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}

static void
test_nc_data_cluster_wl_set_obs (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ngals = 1000;
  NcmMatrix *obs    = ncm_matrix_new (ngals, 4);
  NcmMatrix *dcwl_obs;
  guint i;

  for (i = 0; i < ngals; i++)
  {
    ncm_matrix_set (obs, i, 0, ncm_rng_uniform_gen (rng, 0.0, 2.0));
    ncm_matrix_set (obs, i, 1, ncm_rng_uniform_gen (rng, 0.0, 2.0));
    ncm_matrix_set (obs, i, 2, ncm_rng_uniform_gen (rng, 0.0, 2.0));
    ncm_matrix_set (obs, i, 3, ncm_rng_uniform_gen (rng, 0.0, 2.0));
  }

  nc_data_cluster_wl_set_obs (test->dcwl, obs);

  dcwl_obs = nc_data_cluster_wl_peek_obs (test->dcwl);

  for (i = 0; i < ngals; i++)
  {
    const gdouble w = ncm_matrix_get (obs, i, 0);
    const gdouble x = ncm_matrix_get (obs, i, 1);
    const gdouble y = ncm_matrix_get (obs, i, 2);
    const gdouble z = ncm_matrix_get (obs, i, 3);

    g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 0), ==, w);
    g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 1), ==, x);
    g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 2), ==, y);
    g_assert_cmpfloat (ncm_matrix_get (dcwl_obs, i, 3), ==, z);
  }

  ncm_rng_free (rng);
  ncm_matrix_free (obs);
}

static void
test_nc_data_cluster_wl_serialize (TestNcDataClusterWL *test, gconstpointer pdata)
{
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcHICosmo *cosmo             = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHaloDensityProfile *dp     = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcDistance *dist             = nc_distance_new (100.0);
  NcWLSurfaceMassDensity *smd  = nc_wl_surface_mass_density_new (dist);
  const guint ngals            = 200;
  NcmMatrix *gal_obs           = NULL;
  NcmVector *m2lnP_int_gal     = ncm_vector_new (ngals);
  NcmVector *m2lnP_kde_gal     = ncm_vector_new (ngals);
  NcmVector *m2lnP_int_gal_dup = ncm_vector_new (ngals);
  NcmVector *m2lnP_kde_gal_dup = ncm_vector_new (ngals);
  NcmSerialize *ser;
  gchar *sd_ser;
  NcDataClusterWL *sd_dup;
  guint i;

  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_C_DELTA, 4.0);
  ncm_model_param_set (NCM_MODEL (dp), NC_HALO_DENSITY_PROFILE_LOG10M_DELTA, 14.0);

  nc_wl_surface_mass_density_prepare (smd, cosmo);

  nc_data_cluster_wl_gen_obs (test->dcwl, cosmo, dp, smd, ngals, rng);
  nc_data_cluster_wl_set_prec (test->dcwl, 1.0e-3);
  nc_data_cluster_wl_set_use_kde (test->dcwl, FALSE);

  ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  sd_ser = ncm_serialize_to_string (ser, G_OBJECT (test->dcwl), TRUE);
  sd_dup = NC_DATA_CLUSTER_WL (ncm_serialize_from_string (ser, sd_ser));

  nc_data_cluster_wl_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_int_gal);
  nc_data_cluster_wl_eval_m2lnP (sd_dup, cosmo, dp, smd, m2lnP_int_gal_dup);

  for (i = 0; i < ngals; i++)
  {
    const gdouble m2lnP_int     = ncm_vector_get (m2lnP_int_gal, i);
    const gdouble m2lnP_int_dup = ncm_vector_get (m2lnP_int_gal_dup, i);

    ncm_assert_cmpdouble_e (m2lnP_int, ==, m2lnP_int_dup, 1.0e-3, 0.0);
  }

  nc_data_cluster_wl_kde_eval_m2lnP (test->dcwl, cosmo, dp, smd, m2lnP_kde_gal);
  nc_data_cluster_wl_kde_eval_m2lnP (sd_dup, cosmo, dp, smd, m2lnP_kde_gal_dup);

  for (i = 0; i < ngals; i++)
  {
    const gdouble m2lnP_kde     = ncm_vector_get (m2lnP_kde_gal, i);
    const gdouble m2lnP_kde_dup = ncm_vector_get (m2lnP_kde_gal_dup, i);

    ncm_assert_cmpdouble_e (m2lnP_kde, ==, m2lnP_kde_dup, 1.0e-3, 0.0);
  }

  ncm_vector_free (m2lnP_int_gal);
  ncm_vector_free (m2lnP_kde_gal);
  ncm_serialize_free (ser);
  g_free (sd_ser);
  nc_data_cluster_wl_free (sd_dup);
  ncm_rng_free (rng);
  nc_hicosmo_free (cosmo);
  nc_halo_density_profile_free (dp);
  nc_distance_free (dist);
  nc_wl_surface_mass_density_free (smd);
}


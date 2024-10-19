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

  NcHICosmo *cosmo            = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist            = nc_distance_new (100.0);
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (NC_HALO_DENSITY_PROFILE_MASS_DEF_MEAN, 200.0));
  NcHaloPosition *hp          = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd = nc_wl_surface_mass_density_new (dist);

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
  NcmStatsVec *stats                = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
  const gdouble sigma_obs           = 0.03;
  const guint ntest                 = 10000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    gdouble epsilon_1_out, epsilon_2_out, sigma_obs_out_1, sigma_obs_out_2;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.2, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, sigma_obs, sigma_obs, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, rng);
    nc_galaxy_sd_shape_gauss_data_get (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &sigma_obs_out_1, &sigma_obs_out_2);

    {
      const gdouble z_cl = nc_halo_position_get_redshift (test->halo_position);
      complex double e_s = (s_data->epsilon_int_1 + I * s_data->epsilon_int_2);
      complex double e_o = e_s;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);
      r   = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);
      e_s = e_s * cexp (-2.0 * I * phi);

      if (z_data->z > z_cl)
      {
        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);

        if (fabs (gt) > 1.0)
          e_o = (1.0 + gt * conj (e_s)) / (conj (e_s) + gt);
        else
          e_o = (e_s + gt) / (1.0 + gt * e_s);
      }

      e_s = e_s * cexp (2.0 * I * phi);
      e_o = e_o * cexp (2.0 * I * phi);

      ncm_stats_vec_set (stats, 0, epsilon_1_out - creal (e_o));
      ncm_stats_vec_set (stats, 1, epsilon_2_out - cimag (e_o));

      ncm_stats_vec_update (stats);
    }

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);
  }

  ncm_assert_cmpdouble_e (ncm_stats_vec_get_sd (stats, 0), ==, sigma_obs, 1.0e-1, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_sd (stats, 1), ==, sigma_obs, 1.0e-1, 0.0);

  ncm_stats_vec_free (stats);
}

static void
test_nc_galaxy_sd_shape_gauss_integ (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *z_data   = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data      = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data         = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  NcmStatsVec *stats                  = ncm_stats_vec_new (7, NCM_STATS_VEC_COV, FALSE);
  NcGalaxySDShapeIntegrand *integrand = nc_galaxy_sd_shape_integ (test->galaxy_shape);
  GPtrArray *data_array               = g_ptr_array_new ();
  const gdouble sigma_obs             = 0.03;
  const gdouble sigma_int             = ncm_model_orig_param_get (NCM_MODEL (test->galaxy_shape), NC_GALAXY_SD_SHAPE_GAUSS_SIGMA_INT);
  const guint ntest                   = 10000;
  guint i;

  g_ptr_array_add (data_array, s_data);
  nc_galaxy_sd_shape_integrand_prepare (integrand, test->mset);

  for (i = 0; i < ntest; i++)
  {
    gdouble epsilon_1_out, epsilon_2_out, sigma_obs_out_1, sigma_obs_out_2;
    gdouble int0;

    z_data->z   = ncm_rng_uniform_gen (rng, 0.2, 1.2);
    p_data->ra  = ncm_rng_uniform_gen (rng, -2.0e-2, 2.0e-2);
    p_data->dec = ncm_rng_uniform_gen (rng, -2.0e-2, 2.0e-2);
    nc_galaxy_sd_shape_gauss_gen (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), test->mset, s_data, sigma_obs, sigma_obs, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, rng);
    nc_galaxy_sd_shape_gauss_data_get (NC_GALAXY_SD_SHAPE_GAUSS (test->galaxy_shape), s_data, &epsilon_1_out, &epsilon_2_out, &sigma_obs_out_1, &sigma_obs_out_2);

    g_assert_cmpfloat (hypot (s_data->epsilon_int_1, s_data->epsilon_int_2), <, 1.0);
    g_assert_cmpfloat (hypot (epsilon_1_out, epsilon_2_out), <, 1.2);

    nc_galaxy_sd_shape_prepare_data_array (test->galaxy_shape, test->mset, data_array);
    int0 = nc_galaxy_sd_shape_integrand_eval (integrand, z_data->z, s_data);

    {
      const gdouble z_cl      = nc_halo_position_get_redshift (test->halo_position);
      complex double e_o      = epsilon_1_out + I * epsilon_2_out;
      complex double data_e_s = s_data->epsilon_int_1 + I * s_data->epsilon_int_2;
      complex double e_s      = 0.0;
      complex double g        = 0.0;
      gdouble theta, phi, gt, r;

      nc_halo_position_polar_angles (test->halo_position, p_data->ra, p_data->dec, &theta, &phi);
      r = nc_halo_position_projected_radius (test->halo_position, test->cosmo, theta);

      if (z_data->z > z_cl)
      {
        gt = nc_wl_surface_mass_density_reduced_shear (test->surface_mass_density,
                                                       test->density_profile,
                                                       test->cosmo,
                                                       r, z_data->z, z_cl, z_cl);
        g = gt * cexp (2.0 * I * phi);

        if (gt > 1.0)
          e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
        else
          e_s = (e_o - g) / (1.0 - conj (g) * e_o);
      }

      {
        const gdouble var_int     = gsl_pow_2 (sigma_int);
        const gdouble total_var_1 = var_int + gsl_pow_2 (sigma_obs_out_1);
        const gdouble total_var_2 = var_int + gsl_pow_2 (sigma_obs_out_2);
        const gdouble chi2_1      = gsl_pow_2 (creal (e_s)) / total_var_1;
        const gdouble chi2_2      = gsl_pow_2 (cimag (e_s)) / total_var_2;
        const gdouble m2ln_int1   = chi2_1 + chi2_2 + log (2.0 * M_PI * total_var_1 * total_var_2);

        ncm_assert_cmpdouble_e (-2.0 * log (int0), ==, m2ln_int1, 0.0, 1.0e-13);
      }

      e_s      = e_s * cexp (-2.0 * I * phi);
      e_o      = e_o * cexp (-2.0 * I * phi);
      data_e_s = data_e_s * cexp (-2.0 * I * phi);

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

  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 0), ==, 0.0, 0.0, 1.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 1), ==, 0.0, 0.0, 1.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 2), ==, 0.0, 0.0, 1.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 3), ==, 0.0, 0.0, 1.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 4), ==, ncm_stats_vec_get_mean (stats, 6), 0.0, 1.0e-2);
  ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (stats, 5), ==, 0.0, 0.0, 1.0e-2);

  nc_galaxy_sd_shape_integrand_free (integrand);
  ncm_stats_vec_free (stats);
}

static void
test_nc_galaxy_sd_shape_gauss_required_columns (TestNcGalaxySDShapeGauss *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (test->galaxy_redshift);
  NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (test->galaxy_position, z_data);
  NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (test->galaxy_shape, p_data);
  GList *columns                    = nc_galaxy_sd_shape_data_required_columns (s_data);

  g_assert_cmpuint (g_list_length (columns), ==, 9);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1);
  g_assert_cmpstr (g_list_nth_data (columns, 1), ==, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2);
  g_assert_cmpstr (g_list_nth_data (columns, 2), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_1);
  g_assert_cmpstr (g_list_nth_data (columns, 3), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_EPSILON_OBS_2);
  g_assert_cmpstr (g_list_nth_data (columns, 4), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS_1);
  g_assert_cmpstr (g_list_nth_data (columns, 5), ==, NC_GALAXY_SD_SHAPE_GAUSS_COL_SIGMA_OBS_2);

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_shape_data_free (s_data);
}


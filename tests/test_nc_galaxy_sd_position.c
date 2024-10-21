/***************************************************************************
 *            test_nc_galaxy_sd_position_flat.c
 *
 *  Mon February 27 11:01:27 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2023 <caiolimadeoliveira@pm.me>
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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>


typedef struct _TestNcGalaxySDPosition
{
  NcGalaxySDObsRedshift *gsdor;
  NcGalaxySDPosition *gsdp;
  NcmMSet *mset;
} TestNcGalaxySDPosition;

static void test_nc_galaxy_sd_position_flat_new (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_free (TestNcGalaxySDPosition *test, gconstpointer pdata);

static void test_nc_galaxy_sd_position_lim (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_serialize (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_model_id (TestNcGalaxySDPosition *test, gconstpointer pdata);

static void test_nc_galaxy_sd_position_flat_gen (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_flat_integ (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_flat_required_columns (TestNcGalaxySDPosition *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_position_flat/lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/serialize", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_serialize,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/model_id", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_model_id,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/gen", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_gen,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/integ", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_integ,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/required_columns", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_required_columns,
              &test_nc_galaxy_sd_position_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_position_flat_new (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  gdouble ra_min, ra_max, dec_min, dec_max;

  ra_min = g_test_rand_double_range (0.0, 360.0);

  do {
    ra_max = g_test_rand_double_range (0.0, 360.0);
  } while (ra_max <= ra_min);

  dec_min = g_test_rand_double_range (-90.0, 90.0);

  do {
    dec_max = g_test_rand_double_range (-90.0, 90.0);
  } while (dec_max <= dec_min);

  {
    NcGalaxySDTrueRedshift *gsdtr    = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (0.0, 1100.0));
    NcGalaxySDObsRedshift *gsdor     = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_spec_new (gsdtr));
    NcGalaxySDPositionFlat *gsdpflat = nc_galaxy_sd_position_flat_new (ra_min, ra_max, dec_min, dec_max);

    test->gsdp  = NC_GALAXY_SD_POSITION (gsdpflat);
    test->gsdor = gsdor;
    test->mset  = ncm_mset_empty_new ();

    g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpflat));

    nc_galaxy_sd_true_redshift_free (gsdtr);
  }
}

static void
test_nc_galaxy_sd_position_free (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  nc_galaxy_sd_obs_redshift_free (test->gsdor);
  ncm_mset_free (test->mset);
  NCM_TEST_FREE (nc_galaxy_sd_position_free, test->gsdp);
}

static void
test_nc_galaxy_sd_position_lim (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  guint i;

  for (i = 0; i < 10000; i++)
  {
    gdouble ra_min, ra_max, dec_min, dec_max;
    gdouble ra_min_peek, ra_max_peek, dec_min_peek, dec_max_peek;

    ra_min = g_test_rand_double_range (0.0, 360.0);

    do {
      ra_max = g_test_rand_double_range (0.0, 360.0);
    } while (ra_max <= ra_min);

    dec_min = g_test_rand_double_range (-90.0, 90.0);

    do {
      dec_max = g_test_rand_double_range (-90.0, 90.0);
    } while (dec_max <= dec_min);


    nc_galaxy_sd_position_set_ra_lim (test->gsdp, ra_min, ra_max);
    nc_galaxy_sd_position_set_dec_lim (test->gsdp, dec_min, dec_max);

    nc_galaxy_sd_position_get_ra_lim (test->gsdp, &ra_min_peek, &ra_max_peek);
    nc_galaxy_sd_position_get_dec_lim (test->gsdp, &dec_min_peek, &dec_max_peek);

    g_assert_cmpint (ra_min, ==, ra_min_peek);
    g_assert_cmpint (ra_max, ==, ra_max_peek);
    g_assert_cmpint (dec_min, ==, dec_min_peek);
    g_assert_cmpint (dec_max, ==, dec_max_peek);
  }
}

static void
test_nc_galaxy_sd_position_serialize (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *gsdp_ser              = ncm_serialize_to_string (ser, G_OBJECT (test->gsdp), TRUE);
  NcGalaxySDPosition *gsdp_dup = NC_GALAXY_SD_POSITION (ncm_serialize_from_string (ser, gsdp_ser));
  gdouble ra_min, ra_max, dec_min, dec_max;
  gdouble ra_min_dup, ra_max_dup, dec_min_dup, dec_max_dup;

  g_assert_true (NC_IS_GALAXY_SD_POSITION (gsdp_dup));

  nc_galaxy_sd_position_get_ra_lim (test->gsdp, &ra_min, &ra_max);
  nc_galaxy_sd_position_get_dec_lim (test->gsdp, &dec_min, &dec_max);
  nc_galaxy_sd_position_get_ra_lim (gsdp_dup, &ra_min_dup, &ra_max_dup);
  nc_galaxy_sd_position_get_dec_lim (gsdp_dup, &dec_min_dup, &dec_max_dup);

  g_assert_cmpfloat (ra_min, ==, ra_min_dup);
  g_assert_cmpfloat (ra_max, ==, ra_max_dup);
  g_assert_cmpfloat (dec_min, ==, dec_min_dup);
  g_assert_cmpfloat (dec_max, ==, dec_max_dup);

  ncm_serialize_free (ser);
  g_free (gsdp_ser);
  nc_galaxy_sd_position_free (gsdp_dup);
}

static void
test_nc_galaxy_sd_position_model_id (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcmMSet *model_set  = ncm_mset_empty_new ();
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmModel *model_dup = ncm_model_dup (NCM_MODEL (test->gsdp), ser);

  ncm_mset_set (model_set, model_dup, NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION (ncm_mset_peek (model_set, nc_galaxy_sd_position_id ())));

  ncm_model_free (model_dup);
  ncm_mset_free (model_set);
  ncm_serialize_free (ser);
}

static void
test_nc_galaxy_sd_position_flat_required_columns (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *sdz_data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDPositionData *gsdp_data   = nc_galaxy_sd_position_data_new (test->gsdp, sdz_data);
  GList *columns                      = nc_galaxy_sd_position_data_required_columns (gsdp_data);

  g_assert_cmpint (g_list_length (columns), ==, 3);

  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, NC_GALAXY_SD_POSITION_COL_RA);
  g_assert_cmpstr (g_list_nth_data (columns, 1), ==, NC_GALAXY_SD_POSITION_COL_DEC);
  g_assert_cmpstr (g_list_nth_data (columns, 2), ==, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z);

  nc_galaxy_sd_obs_redshift_data_unref (sdz_data);
  nc_galaxy_sd_position_data_unref (gsdp_data);
  g_list_free_full (columns, g_free);
}

static void
test_nc_galaxy_sd_position_flat_gen (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *sdz_data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDPositionData *gsdp_data   = nc_galaxy_sd_position_data_new (test->gsdp, sdz_data);
  const guint nruns                   = 10;
  const guint ndata                   = 10000;
  gdouble ra_min, ra_max, dec_min, dec_max;
  gdouble ra_avg, ra_var, dec_avg, dec_var;
  guint i;

  ra_min = g_test_rand_double_range (0.0, 360.0);

  do {
    ra_max = g_test_rand_double_range (0.0, 360.0);
  } while (ra_max <= ra_min);

  dec_min = g_test_rand_double_range (-90.0, 90.0);

  do {
    dec_max = g_test_rand_double_range (-90.0, 90.0);
  } while (dec_max <= dec_min);

  ra_avg = 0.5 * (ra_max + ra_min);
  ra_var = (ra_max - ra_min) * (ra_max - ra_min) / 12.0;

  {
    const gdouble delta_min  = ncm_c_degree_to_radian (dec_min);
    const gdouble delta_max  = ncm_c_degree_to_radian (dec_max);
    const gdouble delta_p    = 0.5 * (delta_max + delta_min);
    const gdouble delta_m    = 0.5 * (delta_max - delta_min);
    const gdouble delta_mean = delta_p - tan (delta_p) + delta_m * tan (delta_p) / tan (delta_m);
    const gdouble delta_var  = gsl_pow_2 (delta_m / sin (delta_m)) - 1.0 - gsl_pow_2 ((1.0 - delta_m / tan (delta_m)) / cos (delta_p));

    dec_avg = ncm_c_radian_to_degree (delta_mean);
    dec_var = ncm_c_radian_to_degree (ncm_c_radian_to_degree (delta_var));
  }

  nc_galaxy_sd_position_set_ra_lim (test->gsdp, ra_min, ra_max);
  nc_galaxy_sd_position_set_dec_lim (test->gsdp, dec_min, dec_max);

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    guint j;

    for (j = 0; j < ndata; j++)
    {
      gdouble gen_ra, gen_dec;

      nc_galaxy_sd_obs_redshift_spec_gen (NC_GALAXY_SD_OBS_REDSHIFT_SPEC (test->gsdor), test->mset, sdz_data, rng);
      nc_galaxy_sd_position_flat_gen (NC_GALAXY_SD_POSITION_FLAT (test->gsdp), test->mset, gsdp_data, rng);

      gen_ra  = gsdp_data->ra;
      gen_dec = gsdp_data->dec;

      g_assert_cmpfloat (gen_ra, >, ra_min);
      g_assert_cmpfloat (gen_ra, <, ra_max);
      g_assert_cmpfloat (gen_dec, >, dec_min);
      g_assert_cmpfloat (gen_dec, <, dec_max);

      ncm_stats_vec_set (pos_sample, 0, gen_ra);
      ncm_stats_vec_set (pos_sample, 1, gen_dec);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, ra_avg + 5.0 * sqrt (ra_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, ra_avg - 5.0 * sqrt (ra_var / ndata));

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), <, dec_avg + 5.0 * sqrt (dec_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), >, dec_avg - 5.0 * sqrt (dec_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / ra_var - 1.0), <, 1.0e-1);
    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 1) / dec_var - 1.0), <, 1.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  nc_galaxy_sd_obs_redshift_data_unref (sdz_data);
  nc_galaxy_sd_position_data_unref (gsdp_data);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_position_flat_integ (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcmRNG *rng                            = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *sdz_data    = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDPositionData *gsdp_data      = nc_galaxy_sd_position_data_new (test->gsdp, sdz_data);
  NcGalaxySDPositionIntegrand *integrand = nc_galaxy_sd_position_integ (test->gsdp);
  NcmMSet *mset                          = ncm_mset_empty_new ();
  const guint nruns                      = 10000;
  gdouble ra_min, ra_max, dec_min, dec_max;
  guint i;

  ra_min = g_test_rand_double_range (0.0, 360.0);

  do {
    ra_max = g_test_rand_double_range (0.0, 360.0);
  } while (ra_max <= ra_min);

  dec_min = g_test_rand_double_range (-90.0, 90.0);

  do {
    dec_max = g_test_rand_double_range (-90.0, 90.0);
  } while (dec_max <= dec_min);

  nc_galaxy_sd_position_set_ra_lim (test->gsdp, ra_min, ra_max);
  nc_galaxy_sd_position_set_dec_lim (test->gsdp, dec_min, dec_max);

  nc_galaxy_sd_position_integrand_prepare (integrand, mset);

  for (i = 0; i < nruns; i++)
  {
    gdouble ra  = g_test_rand_double_range (0.0, 360.0);
    gdouble dec = g_test_rand_double_range (-90.0, 90.0);

    gsdp_data->ra  = ra;
    gsdp_data->dec = dec;

    if ((ra < ra_min) || (ra > ra_max) || (dec < dec_min) || (dec > dec_max))
      g_assert_cmpfloat (nc_galaxy_sd_position_integrand_eval (integrand, gsdp_data), ==, 0.0);
    else
      g_assert_cmpfloat (nc_galaxy_sd_position_integrand_eval (integrand, gsdp_data), >, 0.0);
  }

  nc_galaxy_sd_obs_redshift_data_unref (sdz_data);
  nc_galaxy_sd_position_data_unref (gsdp_data);
  nc_galaxy_sd_position_integrand_free (integrand);
  ncm_rng_free (rng);
  ncm_mset_free (mset);
}


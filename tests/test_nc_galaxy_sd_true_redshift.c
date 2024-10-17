/***************************************************************************
 *            test_nc_galaxy_sd_true_redshift.c
 *
 *  Thu Aug 15 17:12:30 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <code.caio@limadeoliveira.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <code.caio@limadeoliveira.me>
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


typedef struct _TestNcGalaxySDTrueRedshift
{
  NcGalaxySDTrueRedshift *gsdtr;
} TestNcGalaxySDTrueRedshift;

static void test_nc_galaxy_sd_true_redshift_lsst_srd_new (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_true_redshift_free (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_true_redshift_lim (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_true_redshift_serialize (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_true_redshift_model_id (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_true_redshift_lsst_srd_gen (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_true_redshift_lsst_srd_integ (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd/lim", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_new,
              &test_nc_galaxy_sd_true_redshift_lim,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd/serialize", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_new,
              &test_nc_galaxy_sd_true_redshift_serialize,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd/model_id", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_new,
              &test_nc_galaxy_sd_true_redshift_model_id,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd/gen", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_new,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_gen,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd/integ", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_new,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_integ,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd_y10/lim", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new,
              &test_nc_galaxy_sd_true_redshift_lim,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd_y10/serialize", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new,
              &test_nc_galaxy_sd_true_redshift_serialize,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd_y10/model_id", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new,
              &test_nc_galaxy_sd_true_redshift_model_id,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd_y10/gen", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_gen,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_add ("/nc/galaxy_sd_true_redshift_lsst_srd_y10/integ", TestNcGalaxySDTrueRedshift, NULL,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new,
              &test_nc_galaxy_sd_true_redshift_lsst_srd_integ,
              &test_nc_galaxy_sd_true_redshift_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_true_redshift_lsst_srd_new (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  gdouble z_min, z_max;

  z_min = g_test_rand_double_range (0.0, 0.5);

  do {
    z_max = g_test_rand_double_range (0.3, 5.0);
  } while (z_max <= z_min);

  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = nc_galaxy_sd_true_redshift_lsst_srd_new (z_min, z_max);

  test->gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (gsdtrlsst);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtrlsst));
}

static void
test_nc_galaxy_sd_true_redshift_lsst_srd_y10_new (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  gdouble z_min, z_max;

  z_min = g_test_rand_double_range (0.0, 0.5);

  do {
    z_max = g_test_rand_double_range (0.3, 5.0);
  } while (z_max <= z_min);

  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlssty10 = nc_galaxy_sd_true_redshift_lsst_srd_new_y10 (z_min, z_max);

  test->gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (gsdtrlssty10);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtrlssty10));
}

static void
test_nc_galaxy_sd_true_redshift_free (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_free, test->gsdtr);
}

static void
test_nc_galaxy_sd_true_redshift_lim (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  gdouble z_min, z_max;
  gdouble z_min_peek, z_max_peek;

  z_min = g_test_rand_double_range (0.0, 0.5);

  do {
    z_max = g_test_rand_double_range (0.3, 5.0);
  } while (z_max <= z_min);

  nc_galaxy_sd_true_redshift_set_lim (test->gsdtr, z_min, z_max);


  nc_galaxy_sd_true_redshift_get_lim (test->gsdtr, &z_min_peek, &z_max_peek);

  g_assert_cmpfloat (z_min, ==, z_min_peek);
  g_assert_cmpfloat (z_max, ==, z_max_peek);
}

static void
test_nc_galaxy_sd_true_redshift_serialize (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  NcGalaxySDTrueRedshift *gsdtr     = test->gsdtr;
  NcmSerialize *ser                 = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *gsdtr_ser                  = ncm_serialize_to_string (ser, G_OBJECT (gsdtr), TRUE);
  NcGalaxySDTrueRedshift *gsdtr_dup = NC_GALAXY_SD_TRUE_REDSHIFT (ncm_serialize_from_string (ser, gsdtr_ser));
  gdouble z_min, z_max;
  gdouble z_min_dup, z_max_dup;

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT (gsdtr_dup));

  nc_galaxy_sd_true_redshift_get_lim (gsdtr, &z_min, &z_max);
  nc_galaxy_sd_true_redshift_get_lim (gsdtr_dup, &z_min_dup, &z_max_dup);

  g_assert_cmpfloat (z_min, ==, z_min_dup);
  g_assert_cmpfloat (z_max, ==, z_max_dup);

  ncm_serialize_clear (&ser);
  g_free (gsdtr_ser);
  nc_galaxy_sd_true_redshift_clear (&gsdtr_dup);
}

static void
test_nc_galaxy_sd_true_redshift_model_id (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  NcmMSet *model_set                = ncm_mset_empty_new ();
  NcmSerialize *ser                 = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmModel *model_dup               = ncm_model_dup (NCM_MODEL (test->gsdtr), ser);
  NcGalaxySDObsRedshiftSpec *gsdors = nc_galaxy_sd_obs_redshift_spec_new (NC_GALAXY_SD_TRUE_REDSHIFT (model_dup));

  ncm_mset_set (model_set, NCM_MODEL (gsdors), NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT (ncm_mset_peek (model_set, nc_galaxy_sd_true_redshift_id ())));

  ncm_model_free (model_dup);
  ncm_mset_clear (&model_set);
  ncm_serialize_clear (&ser);
  nc_galaxy_sd_obs_redshift_spec_free (gsdors);
}

static gdouble
_test_get_z_avg (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst)
{
  gdouble z_min, z_max;

  nc_galaxy_sd_true_redshift_get_lim (NC_GALAXY_SD_TRUE_REDSHIFT (gsdtrlsst), &z_min, &z_max);
  {
    const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "alpha", NULL);
    const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "beta", NULL);
    const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "z0", NULL);
    const gdouble y_lb  = pow (z_min / z0, alpha);
    const gdouble y_ub  = pow (z_max / z0, alpha);

    return z0 *
           (gsl_sf_gamma_inc_Q ((2.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((2.0 + beta) / alpha, y_ub)) /
           (gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_ub)) *
           gsl_sf_poch ((1.0 + beta) / alpha, 1.0 / alpha);
  }
}

static gdouble
_test_get_z_var (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst)
{
  gdouble z_min, z_max;

  nc_galaxy_sd_true_redshift_get_lim (NC_GALAXY_SD_TRUE_REDSHIFT (gsdtrlsst), &z_min, &z_max);
  {
    const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "alpha", NULL);
    const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "beta", NULL);
    const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "z0", NULL);
    const gdouble y_lb  = pow (z_min / z0, alpha);
    const gdouble y_ub  = pow (z_max / z0, alpha);

    return z0 * z0 *
           (gsl_sf_gamma_inc_Q ((3.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((3.0 + beta) / alpha, y_ub)) /
           (gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_ub)) *
           gsl_sf_poch ((1.0 + beta) / alpha, 2.0 / alpha) -
           gsl_pow_2 (_test_get_z_avg (gsdtrlsst));
  }
}

static void
test_nc_galaxy_sd_true_redshift_lsst_srd_gen (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (test->gsdtr);
  NcmRNG *rng                              = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  gdouble z_min, z_max;

  z_min = g_test_rand_double_range (0.0, 0.5);

  do {
    z_max = g_test_rand_double_range (0.3, 5.0);
  } while (z_max <= z_min);

  nc_galaxy_sd_true_redshift_set_lim (test->gsdtr, z_min, z_max);

  {
    const gdouble z_avg = _test_get_z_avg (gsdtrlsst);
    const gdouble z_var = _test_get_z_var (gsdtrlsst);
    const gint nruns    = 10;
    const gint ndata    = 10000;
    gint i;

    for (i = 0; i < nruns; i++)
    {
      NcmStatsVec *pos_sample = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
      gint j;

      for (j = 0; j < ndata; j++)
      {
        const gdouble gen_z = nc_galaxy_sd_true_redshift_gen (test->gsdtr, rng);

        ncm_stats_vec_set (pos_sample, 0, gen_z);

        ncm_stats_vec_update (pos_sample);
      }

      g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 5.0 * sqrt (z_var / ndata));
      g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 5.0 * sqrt (z_var / ndata));

      g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);

      ncm_stats_vec_free (pos_sample);
    }
  }

  ncm_rng_clear (&rng);
}

static void
test_nc_galaxy_sd_true_redshift_lsst_srd_integ (TestNcGalaxySDTrueRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint nruns = 10000;
  gdouble z_min, z_max;
  guint i;

  z_min = g_test_rand_double_range (0.0, 0.5);

  do {
    z_max = g_test_rand_double_range (0.3, 5.0);
  } while (z_max <= z_min);

  nc_galaxy_sd_true_redshift_set_lim (test->gsdtr, z_min, z_max);

  for (i = 0; i < nruns; i++)
  {
    gdouble z = g_test_rand_double_range (0.0, 5.0);

    if ((z < z_min) || (z > z_max))
      g_assert_cmpfloat (nc_galaxy_sd_true_redshift_integ (test->gsdtr, z), ==, 0.0);
    else
      g_assert_cmpfloat (nc_galaxy_sd_true_redshift_integ (test->gsdtr, z), >, 0.0);
  }

  ncm_rng_clear (&rng);
}


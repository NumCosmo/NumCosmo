/***************************************************************************
 *            test_nc_galaxy_sd_obs_redshift.c
 *
 *  Thu Aug 15 17:12:30 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <code.caio@limadeoliveira.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <code.caio@limadeoliveira.me>
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
#include "galaxy/nc_galaxy_sd_obs_redshift_gauss.c"

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>


typedef struct _TestNcGalaxySDObsRedshift
{
  NcGalaxySDObsRedshift *gsdor;
  NcmMSet *mset;
  gdouble z_min;
  gdouble z_max;
} TestNcGalaxySDObsRedshift;

static void test_nc_galaxy_sd_obs_redshift_spec_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_pz_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_free (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_serialize (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_model_id (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_spec_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_spec_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_spec_required_columns (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_spec_lim (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_gauss_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_required_columns (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_lim (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_pz_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_pz_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_pz_required_columns (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_pz_lim (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_gauss_data_setget (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_pz_data_setget (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/serialize", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_serialize,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/model_id", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_model_id,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/gen", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_spec_gen,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/integ", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_spec_integ,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/required_columns", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_spec_required_columns,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/lim", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_spec_lim,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/serialize", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_serialize,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/model_id", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_model_id,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/gen", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_gauss_gen,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/integ", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_gauss_integ,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/required_columns", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_gauss_required_columns,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/data/lim", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_gauss_lim,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/data/setget", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_gauss_data_setget,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/serialize", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_serialize,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/model_id", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_model_id,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/gen", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_pz_gen,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/integ", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_pz_integ,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/required_columns", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_pz_required_columns,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/data/lim", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_pz_lim,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_add ("/nc/galaxy_sd_obs_redshift/pz/data/setget", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_pz_new,
              &test_nc_galaxy_sd_obs_redshift_pz_data_setget,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_obs_redshift_spec_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  const gdouble z_min                  = 0.01;
  const gdouble z_max                  = 10.0;
  NcGalaxySDTrueRedshift *gsdtr        = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshiftSpec *gsdorspec = nc_galaxy_sd_obs_redshift_spec_new (gsdtr, z_min, z_max);

  test->gsdor = NC_GALAXY_SD_OBS_REDSHIFT (gsdorspec);
  test->z_min = z_min;
  test->z_max = z_max;
  test->mset  = ncm_mset_empty_new ();

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdorspec));

  nc_galaxy_sd_true_redshift_free (gsdtr);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  const gdouble z_min                    = 0.01;
  const gdouble z_max                    = 10.0;
  NcGalaxySDTrueRedshift *gsdtr          = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshiftGauss *gsdorgauss = nc_galaxy_sd_obs_redshift_gauss_new (gsdtr, z_min, z_max);

  test->gsdor = NC_GALAXY_SD_OBS_REDSHIFT (gsdorgauss);
  test->z_min = z_min;
  test->z_max = z_max;
  test->mset  = ncm_mset_empty_new ();

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

  nc_galaxy_sd_true_redshift_free (gsdtr);
}

static void
test_nc_galaxy_sd_obs_redshift_pz_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  const gdouble z_min              = 0.01;
  const gdouble z_max              = 10.0;
  NcGalaxySDObsRedshiftPz *gsdorpz = nc_galaxy_sd_obs_redshift_pz_new ();

  test->gsdor = NC_GALAXY_SD_OBS_REDSHIFT (gsdorpz);
  test->z_min = z_min;
  test->z_max = z_max;
  test->mset  = ncm_mset_empty_new ();

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_PZ (gsdorpz));
}

static void
test_nc_galaxy_sd_obs_redshift_free (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  ncm_mset_free (test->mset);
  NCM_TEST_FREE (nc_galaxy_sd_obs_redshift_free, test->gsdor);
}

static void
test_nc_galaxy_sd_obs_redshift_serialize (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshift *gsdor     = test->gsdor;
  NcmSerialize *ser                = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *gsdor_ser                 = ncm_serialize_to_string (ser, G_OBJECT (gsdor), TRUE);
  NcGalaxySDObsRedshift *gsdor_dup = NC_GALAXY_SD_OBS_REDSHIFT (ncm_serialize_from_string (ser, gsdor_ser));

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor_dup));

  g_free (gsdor_ser);
  ncm_serialize_free (ser);
  nc_galaxy_sd_obs_redshift_free (gsdor_dup);
}

static void
test_nc_galaxy_sd_obs_redshift_model_id (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmMSet *model_set  = ncm_mset_empty_new ();
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmModel *model_dup = ncm_model_dup (NCM_MODEL (test->gsdor), ser);

  ncm_mset_set (model_set, model_dup, NULL);

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT (ncm_mset_peek (model_set, nc_galaxy_sd_obs_redshift_id ())));

  ncm_model_free (model_dup);
  ncm_mset_free (model_set);
  ncm_serialize_free (ser);
}

static gdouble
_test_get_z_avg (TestNcGalaxySDObsRedshift *test)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = nc_galaxy_sd_true_redshift_lsst_srd_new ();

  const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "alpha", NULL);
  const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "beta", NULL);
  const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "z0", NULL);
  const gdouble y_lb  = pow (test->z_min / z0, alpha);
  const gdouble y_ub  = pow (test->z_max / z0, alpha);

  nc_galaxy_sd_true_redshift_lsst_srd_free (gsdtrlsst);

  return z0 *
         (gsl_sf_gamma_inc_Q ((2.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((2.0 + beta) / alpha, y_ub)) /
         (gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_ub)) *
         gsl_sf_poch ((1.0 + beta) / alpha, 1.0 / alpha);
}

static gdouble
_test_get_z_var (TestNcGalaxySDObsRedshift *test)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = nc_galaxy_sd_true_redshift_lsst_srd_new ();

  const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "alpha", NULL);
  const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "beta", NULL);
  const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "z0", NULL);
  const gdouble y_lb  = pow (test->z_min / z0, alpha);
  const gdouble y_ub  = pow (test->z_max / z0, alpha);

  nc_galaxy_sd_true_redshift_lsst_srd_free (gsdtrlsst);

  return z0 * z0 *
         (gsl_sf_gamma_inc_Q ((3.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((3.0 + beta) / alpha, y_ub)) /
         (gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_ub)) *
         gsl_sf_poch ((1.0 + beta) / alpha, 2.0 / alpha) -
         gsl_pow_2 (_test_get_z_avg (test));
}

static void
test_nc_galaxy_sd_obs_redshift_spec_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  const gdouble z_avg             = _test_get_z_avg (test);
  const gdouble z_var             = _test_get_z_var (test);
  const guint nruns               = 10;
  const guint ndata               = 10000;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
    guint j;

    nc_galaxy_sd_obs_redshift_prepare (test->gsdor, data);

    for (j = 0; j < ndata; j++)
    {
      nc_galaxy_sd_obs_redshift_spec_gen (NC_GALAXY_SD_OBS_REDSHIFT_SPEC (test->gsdor), test->mset, data, rng);

      g_assert_cmpfloat (data->z, >=, test->z_min);
      g_assert_cmpfloat (data->z, <=, test->z_max);

      ncm_stats_vec_set (pos_sample, 0, data->z);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 5.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 5.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
    guint j;

    nc_galaxy_sd_obs_redshift_prepare (test->gsdor, data);

    for (j = 0; j < ndata; j++)
    {
      gboolean result = nc_galaxy_sd_obs_redshift_spec_gen1 (NC_GALAXY_SD_OBS_REDSHIFT_SPEC (test->gsdor), test->mset, data, rng);

      if (result)
        g_assert_true ((data->z >= test->z_min) && (data->z <= test->z_max));
      else
        g_assert_true ((data->z <= test->z_min) || (data->z >= test->z_max));

      ncm_stats_vec_set (pos_sample, 0, data->z);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 5.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 5.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  const gdouble sigma             = 0.05;
  const gdouble z_avg             = _test_get_z_avg (test);
  const gdouble z_var             = sqrt (gsl_pow_4 (sigma * (1 + z_avg)) + gsl_pow_2 (_test_get_z_var (test)));
  const guint nruns               = 10;
  const guint ndata               = 10000;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    guint j;

    nc_galaxy_sd_obs_redshift_prepare (test->gsdor, data);

    for (j = 0; j < ndata; j++)
    {
      nc_galaxy_sd_obs_redshift_gauss_gen (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), test->mset, data, sigma, rng);

      g_assert_cmpfloat (ldata->zp, >=, test->z_min);
      g_assert_cmpfloat (ldata->zp, <=, test->z_max);

      ncm_stats_vec_set (pos_sample, 0, data->z);
      ncm_stats_vec_set (pos_sample, 1, ldata->zp);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 6.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 6.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), <, z_avg + 6.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), >, z_avg - 6.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);
    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 1) / z_var - 1.0), <, 2.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    guint j;

    nc_galaxy_sd_obs_redshift_prepare (test->gsdor, data);

    for (j = 0; j < ndata; j++)
    {
      gboolean result = nc_galaxy_sd_obs_redshift_gauss_gen1 (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), test->mset, data, sigma, rng);

      if (result)
        g_assert_true ((ldata->zp >= test->z_min) && (ldata->zp <= test->z_max));
      else
        g_assert_true ((ldata->zp <= test->z_min) || (ldata->zp >= test->z_max));

      ncm_stats_vec_set (pos_sample, 0, data->z);
      ncm_stats_vec_set (pos_sample, 1, ldata->zp);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 6.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 6.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), <, z_avg + 6.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), >, z_avg - 6.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);
    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 1) / z_var - 1.0), <, 2.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_obs_redshift_pz_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint nruns   = 10;
  const guint npoints = 10000;
  const guint ndata   = 10000;
  gdouble z_avg;
  gdouble z_sd;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
    NcmStatsVec *pos_sample         = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
    NcmVector *xv                   = ncm_vector_new (npoints);
    NcmVector *yv                   = ncm_vector_new (npoints);
    NcmSpline *pz;
    guint j;
    guint k;

    z_avg = g_test_rand_double_range (test->z_min, test->z_max);
    z_sd  = 0.03 * (1.0 + z_avg);

    for (j = 0; j < npoints; j++)
    {
      gdouble x = z_avg - 5.0 * z_sd + 10.0 * z_sd * j / ((gdouble) npoints - 1.0);
      gdouble y = exp (-0.5 * gsl_pow_2 ((x - z_avg) / z_sd)) / (sqrt (2.0 * M_PI) * z_sd);

      ncm_vector_set (xv, j, x);
      ncm_vector_set (yv, j, y);
    }

    pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

    nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), data, pz);
    nc_galaxy_sd_obs_redshift_pz_prepare (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), data);

    for (k = 0; k < ndata; k++)
    {
      nc_galaxy_sd_obs_redshift_pz_gen (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), test->mset, data, rng);
      ncm_stats_vec_set (pos_sample, 0, data->z);
      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 6.0 * z_sd / sqrt (npoints));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 6.0 * z_sd / sqrt (npoints));

    nc_galaxy_sd_obs_redshift_data_unref (data);
    ncm_stats_vec_free (pos_sample);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_spline_free (pz);
  }

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng                               = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *data           = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDObsRedshiftIntegrand *integrand = nc_galaxy_sd_obs_redshift_integ (test->gsdor);
  NcGalaxySDTrueRedshift *gsdtr             = NC_GALAXY_SD_TRUE_REDSHIFT (ncm_model_peek_submodel_by_mid (NCM_MODEL (test->gsdor), nc_galaxy_sd_true_redshift_id ()));
  NcmMSet *mset                             = ncm_mset_empty_new ();
  const gdouble zp                          = 0.4;
  const gdouble sigma0                      = 0.05;
  const guint nruns                         = 10000;
  guint i;

  nc_galaxy_sd_obs_redshift_gauss_set_use_true_z (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), TRUE);

  nc_galaxy_sd_obs_redshift_gauss_data_set (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), data, zp, sigma0, sigma0);
  nc_galaxy_sd_obs_redshift_integrand_prepare (integrand, mset);

  for (i = 0; i < nruns; i++)
  {
    const gdouble z       = g_test_rand_double_range (0.0, 5.0);
    const gdouble sigmaz  = sigma0 * (1.0 + z);
    const gdouble norm    = sqrt (2.0 * M_PI) * sigmaz * 0.5 * (1.0 + erf (z / (M_SQRT2 * sigmaz)));
    const gdouble int_z   = nc_galaxy_sd_true_redshift_integ (gsdtr, z);
    const gdouble int_zp  = exp (-0.5 * gsl_pow_2 ((zp - z) / sigmaz)) / norm;
    const gdouble control = int_z * int_zp;
    const gdouble res     = nc_galaxy_sd_obs_redshift_integrand_eval (integrand, z, data);

    ncm_assert_cmpdouble_e (control, ==, res, 1.0e-10, 0.0);
  }

  nc_galaxy_sd_obs_redshift_gauss_set_use_true_z (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), FALSE);

  for (i = 0; i < nruns; i++)
  {
    const gdouble z       = g_test_rand_double_range (0.0, 5.0);
    const gdouble norm    = sqrt (2.0 * M_PI) * sigma0;
    const gdouble chi2    = exp (-0.5 * gsl_pow_2 ((zp - z) / sigma0));
    const gdouble control = chi2 / norm;
    const gdouble res     = nc_galaxy_sd_obs_redshift_integrand_eval (integrand, z, data);

    ncm_assert_cmpdouble_e (control, ==, res, 1.0e-10, 0.0);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
  nc_galaxy_sd_obs_redshift_integrand_free (integrand);
  ncm_mset_free (mset);
  ncm_rng_clear (&rng);
}

static void
test_nc_galaxy_sd_obs_redshift_spec_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng                               = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *data           = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDObsRedshiftIntegrand *integrand = nc_galaxy_sd_obs_redshift_integ (test->gsdor);
  NcGalaxySDTrueRedshift *gsdtr             = NC_GALAXY_SD_TRUE_REDSHIFT (ncm_model_peek_submodel_by_mid (NCM_MODEL (test->gsdor), nc_galaxy_sd_true_redshift_id ()));
  NcmMSet *mset                             = ncm_mset_empty_new ();
  const guint nruns                         = 10000;
  guint i;

  nc_galaxy_sd_obs_redshift_integrand_prepare (integrand, mset);

  for (i = 0; i < nruns; i++)
  {
    const gdouble z          = g_test_rand_double_range (0.0, 5.0);
    const gdouble int_true_z = nc_galaxy_sd_true_redshift_integ (gsdtr, z);
    const gdouble int_obs_z  = nc_galaxy_sd_obs_redshift_integrand_eval (integrand, z, data);

    ncm_assert_cmpdouble_e (int_obs_z, ==, int_true_z, 1.0e-10, 0.0);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
  nc_galaxy_sd_obs_redshift_integrand_free (integrand);
  ncm_mset_free (mset);
  ncm_rng_clear (&rng);
}

static void
test_nc_galaxy_sd_obs_redshift_pz_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng                               = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxySDObsRedshiftData *data           = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  NcGalaxySDObsRedshiftIntegrand *integrand = nc_galaxy_sd_obs_redshift_integ (test->gsdor);
  NcmMSet *mset                             = ncm_mset_empty_new ();
  const guint nruns                         = 10000;
  const guint npoints                       = 10000;
  gdouble z_avg                             = g_test_rand_double_range (test->z_min, test->z_max);
  gdouble z_sd                              = 0.03 * (1.0 + z_avg);
  NcmVector *xv                             = ncm_vector_new (npoints);
  NcmVector *yv                             = ncm_vector_new (npoints);
  NcmSpline *pz;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z = z_avg - 5.0 * z_sd + 10.0 * z_sd * i / ((gdouble) npoints - 1.0);
    const gdouble f = exp (-0.5 * gsl_pow_2 ((z - z_avg) / z_sd)) / (sqrt (2.0 * M_PI) * z_sd);

    ncm_vector_set (xv, i, z);
    ncm_vector_set (yv, i, f);
  }

  pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), data, pz);
  nc_galaxy_sd_obs_redshift_integrand_prepare (integrand, mset);

  for (i = 0; i < nruns; i++)
  {
    const gdouble z     = g_test_rand_double_range (ncm_vector_get (xv, 0), ncm_vector_get (xv, npoints - 1));
    const gdouble f     = exp (-0.5 * gsl_pow_2 ((z - z_avg) / z_sd)) / (sqrt (2.0 * M_PI) * z_sd);
    const gdouble integ = nc_galaxy_sd_obs_redshift_integrand_eval (integrand, z, data);

    ncm_assert_cmpdouble_e (f, ==, integ, 1.0e-6, 0.0);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
  nc_galaxy_sd_obs_redshift_integrand_free (integrand);
  ncm_spline_clear (&pz);
  ncm_mset_free (mset);
  ncm_rng_clear (&rng);
  ncm_vector_clear (&xv);
  ncm_vector_clear (&yv);
}

static void
test_nc_galaxy_sd_obs_redshift_spec_required_columns (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  GList *columns                  = nc_galaxy_sd_obs_redshift_data_required_columns (data);


  g_assert_cmpuint (g_list_length (columns), ==, 1);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, "z");

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_required_columns (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  GList *columns                  = nc_galaxy_sd_obs_redshift_data_required_columns (data);

  g_assert_cmpuint (g_list_length (columns), ==, 4);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, "z");
  g_assert_cmpstr (g_list_nth_data (columns, 1), ==, "zp");
  g_assert_cmpstr (g_list_nth_data (columns, 2), ==, "sigma_0");
  g_assert_cmpstr (g_list_nth_data (columns, 3), ==, "sigma_z");

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_spec_lim (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  gdouble z_min;
  gdouble z_max;

  nc_galaxy_sd_obs_redshift_get_lim (test->gsdor, data, &z_min, &z_max);

  g_assert_cmpfloat (z_min, ==, test->z_min);
  g_assert_cmpfloat (z_max, ==, test->z_max);

  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_lim (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  gdouble z_min;
  gdouble z_max;

  nc_galaxy_sd_obs_redshift_get_lim (test->gsdor, data, &z_min, &z_max);

  g_assert_cmpfloat (z_min, ==, test->z_min);
  g_assert_cmpfloat (z_max, ==, test->z_max);

  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_pz_lim (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  gdouble npoints                 = 1000;
  NcmVector *xv                   = ncm_vector_new (npoints);
  NcmVector *yv                   = ncm_vector_new (npoints);
  NcmSpline *pz;
  gdouble z_min;
  gdouble z_max;
  gdouble z_avg;
  gdouble z_sd;
  guint i;

  z_avg = g_test_rand_double_range (test->z_min, test->z_max);
  z_sd  = 0.03 * (1.0 + z_avg);

  for (i = 0; i < npoints; i++)
  {
    gdouble z = z_avg - 5 * z_sd + 10 * z_sd * i / (npoints - 1);
    gdouble f = exp (-0.5 * gsl_pow_2 (z - _test_get_z_avg (test)) / _test_get_z_var (test)) / (sqrt (2.0 * M_PI) * sqrt (_test_get_z_var (test)));

    ncm_vector_set (xv, i, z);
    ncm_vector_set (yv, i, f);
  }

  pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), data, pz);

  nc_galaxy_sd_obs_redshift_get_lim (test->gsdor, data, &z_min, &z_max);

  g_assert_cmpfloat_with_epsilon (z_min, z_avg - 5 * z_sd, 1e-10);
  g_assert_cmpfloat_with_epsilon (z_max, z_avg + 5 * z_sd, 1e-10);

  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_pz_required_columns (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  GList *columns                  = nc_galaxy_sd_obs_redshift_data_required_columns (data);

  g_assert_cmpuint (g_list_length (columns), ==, 1);
  g_assert_cmpstr (g_list_nth_data (columns, 0), ==, "z");

  g_list_free_full (columns, g_free);
  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_data_setget (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  const guint ntests              = 1000;
  guint i;

  for (i = 0; i < ntests; i++)
  {
    const gdouble zp      = g_test_rand_double_range (0.0, 5.0);
    const gdouble sigma0  = g_test_rand_double_range (0.01, 1.1);
    const gdouble sigma_z = g_test_rand_double_range (0.01, 1.1);
    gdouble zp_out;
    gdouble sigma0_out;
    gdouble sigma_z_out;

    nc_galaxy_sd_obs_redshift_gauss_data_set (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), data, zp, sigma0, sigma_z);
    nc_galaxy_sd_obs_redshift_gauss_data_get (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (test->gsdor), data, &zp_out, &sigma0_out, &sigma_z_out);

    g_assert_cmpfloat (zp_out, ==, zp);
    g_assert_cmpfloat (sigma0_out, ==, sigma0);
    g_assert_cmpfloat (sigma_z_out, ==, sigma_z);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
}

static void
test_nc_galaxy_sd_obs_redshift_pz_data_setget (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcGalaxySDObsRedshiftData *data = nc_galaxy_sd_obs_redshift_data_new (test->gsdor);
  const guint ntests              = 10;
  const guint npoints             = 1000;
  const guint ndata               = 1000;
  gdouble z_avg                   = g_test_rand_double_range (test->z_min, test->z_max);
  gdouble z_sd                    = 0.03 * (1.0 + z_avg);
  guint i;

  for (i = 0; i < ntests; i++)
  {
    NcmVector *xv = ncm_vector_new (npoints);
    NcmVector *yv = ncm_vector_new (npoints);
    NcmSpline *pz;
    NcmSpline *pz_out;
    guint j;

    for (j = 0; j < npoints; j++)
    {
      gdouble z = z_avg - 5.0 * z_sd + 10.0 * z_sd * j / ((gdouble) npoints - 1.0);
      gdouble f = exp (-0.5 * gsl_pow_2 (z - z_avg) / z_sd) / (sqrt (2.0 * M_PI) * z_sd);

      ncm_vector_set (xv, j, z);
      ncm_vector_set (yv, j, f);
    }

    pz = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

    nc_galaxy_sd_obs_redshift_pz_data_set (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), data, pz);
    nc_galaxy_sd_obs_redshift_pz_data_get (NC_GALAXY_SD_OBS_REDSHIFT_PZ (test->gsdor), data, &pz_out);

    for (j = 0; j < ndata; j++)
    {
      gdouble z = g_test_rand_double_range (ncm_vector_get (xv, 0), ncm_vector_get (xv, npoints - 1));
      gdouble y;
      gdouble y_out;

      y     = ncm_spline_eval (pz, z);
      y_out = ncm_spline_eval (pz_out, z);

      g_assert_cmpfloat_with_epsilon (y, y_out, 1e-10);
    }

    ncm_vector_clear (&xv);
    ncm_vector_clear (&yv);
    ncm_spline_clear (&pz);
  }

  nc_galaxy_sd_obs_redshift_data_unref (data);
}


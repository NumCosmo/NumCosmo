/***************************************************************************
 *            test_nc_galaxy_sd_obs_redshift.c
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


typedef struct _TestNcGalaxySDObsRedshift
{
  NcGalaxySDObsRedshift *gsdor;
  gdouble z_min;
  gdouble z_max;
} TestNcGalaxySDObsRedshift;

static void test_nc_galaxy_sd_obs_redshift_spec_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_free (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_serialize (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_model_id (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

static void test_nc_galaxy_sd_obs_redshift_spec_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_spec_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_spec_get_header (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);
static void test_nc_galaxy_sd_obs_redshift_gauss_get_header (TestNcGalaxySDObsRedshift *test, gconstpointer pdata);

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

  g_test_add ("/nc/galaxy_sd_obs_redshift/spec/get_header", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_spec_new,
              &test_nc_galaxy_sd_obs_redshift_spec_get_header,
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

  g_test_add ("/nc/galaxy_sd_obs_redshift/gauss/get_header", TestNcGalaxySDObsRedshift, NULL,
              &test_nc_galaxy_sd_obs_redshift_gauss_new,
              &test_nc_galaxy_sd_obs_redshift_gauss_get_header,
              &test_nc_galaxy_sd_obs_redshift_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_obs_redshift_spec_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  gdouble z_min = 0.0;
  gdouble z_max = 1100.0;

  NcGalaxySDTrueRedshift *gsdtr        = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (z_min, z_max));
  NcGalaxySDObsRedshiftSpec *gsdorspec = nc_galaxy_sd_obs_redshift_spec_new (gsdtr);

  test->gsdor = NC_GALAXY_SD_OBS_REDSHIFT (gsdorspec);
  test->z_min = z_min;
  test->z_max = z_max;

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdorspec));

  nc_galaxy_sd_true_redshift_free (gsdtr);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_new (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  gdouble z_min = 0.0;
  gdouble z_max = 1100.0;

  NcGalaxySDTrueRedshift *gsdtr          = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new (z_min, z_max));
  NcGalaxySDObsRedshiftGauss *gsdorgauss = nc_galaxy_sd_obs_redshift_gauss_new (gsdtr);

  test->gsdor = NC_GALAXY_SD_OBS_REDSHIFT (gsdorgauss);
  test->z_min = z_min;
  test->z_max = z_max;

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

  nc_galaxy_sd_true_redshift_free (gsdtr);
}

static void
test_nc_galaxy_sd_obs_redshift_free (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
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

  ncm_serialize_free (ser);
  g_free (gsdor_ser);
  nc_galaxy_sd_obs_redshift_free (gsdor_dup);
}

static void
test_nc_galaxy_sd_obs_redshift_model_id (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmMSet *model_set = ncm_mset_empty_new ();
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  ncm_mset_set (model_set, ncm_model_dup (NCM_MODEL (test->gsdor), ser));

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT (ncm_mset_peek (model_set, nc_galaxy_sd_obs_redshift_id ())));

  ncm_mset_free (model_set);
  ncm_serialize_free (ser);
}

static gdouble
_test_get_z_avg (TestNcGalaxySDObsRedshift *test)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = nc_galaxy_sd_true_redshift_lsst_srd_new (test->z_min, test->z_max);

  const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "alpha");
  const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "beta");
  const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "z0");
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
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = nc_galaxy_sd_true_redshift_lsst_srd_new (test->z_min, test->z_max);

  const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "alpha");
  const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "beta");
  const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdtrlsst), "z0");
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
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *data     = ncm_vector_new (1);
  const gdouble z_avg = _test_get_z_avg (test);
  const gdouble z_var = _test_get_z_var (test);
  const guint nruns   = 10;
  const guint ndata   = 10000;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
    guint j;

    for (j = 0; j < ndata; j++)
    {
      const gdouble z = nc_galaxy_sd_obs_redshift_gen (test->gsdor, rng, data);

      ncm_stats_vec_set (pos_sample, 0, z);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 5.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 5.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  ncm_vector_free (data);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_gen (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *data     = ncm_vector_new (1);
  const gdouble z_avg = _test_get_z_avg (test);
  const gdouble z_var = sqrt (gsl_pow_4 (ncm_model_param_get_by_name (NCM_MODEL (test->gsdor), "sigma") * (1 + z_avg)) + gsl_pow_2 (_test_get_z_var (test)));
  const guint nruns   = 10;
  const guint ndata   = 10000;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (1, NCM_STATS_VEC_COV, FALSE);
    guint j;

    for (j = 0; j < ndata; j++)
    {
      const gdouble z = nc_galaxy_sd_obs_redshift_gen (test->gsdor, rng, data);

      ncm_stats_vec_set (pos_sample, 0, z);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 6.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 6.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  ncm_vector_free (data);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_obs_redshift_spec_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *data   = ncm_vector_new (1);
  const guint nruns = 10000;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    gdouble z = g_test_rand_double_range (0.0, 5.0);

    if ((z < test->z_min) || (z > test->z_max))
      g_assert_cmpfloat (nc_galaxy_sd_obs_redshift_integ (test->gsdor, z, data), ==, 0.0);
    else
      g_assert_cmpfloat (nc_galaxy_sd_obs_redshift_integ (test->gsdor, z, data), >=, 0.0);
  }

  ncm_vector_clear (&data);
  ncm_rng_clear (&rng);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_integ (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *data   = ncm_vector_new (2);
  const guint nruns = 10000;
  guint i;

  for (i = 0; i < nruns; i++)
  {
    gdouble z     = g_test_rand_double_range (0.0, 5.0);
    gdouble sigma = g_test_rand_double_range (0.01, 0.05);
    gdouble zp    = ncm_rng_gaussian_gen (rng, z, sigma);

    ncm_vector_set (data, 0, zp);
    ncm_vector_set (data, 1, sigma);

    if ((z < test->z_min) || (z > test->z_max))
      g_assert_cmpfloat (nc_galaxy_sd_obs_redshift_integ (test->gsdor, z, data), ==, 0.0);
    else
      g_assert_cmpfloat (nc_galaxy_sd_obs_redshift_integ (test->gsdor, z, data), >=, 0.0);
  }

  ncm_vector_clear (&data);
  ncm_rng_clear (&rng);
}

static void
test_nc_galaxy_sd_obs_redshift_spec_get_header (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  GStrv header      = nc_galaxy_sd_obs_redshift_get_header (test->gsdor);
  GStrv header_ctrl = g_strsplit ("z", " ", -1);

  g_assert_cmpstrv (header, header_ctrl);

  g_strfreev (header);
  g_strfreev (header_ctrl);
}

static void
test_nc_galaxy_sd_obs_redshift_gauss_get_header (TestNcGalaxySDObsRedshift *test, gconstpointer pdata)
{
  GStrv header      = nc_galaxy_sd_obs_redshift_get_header (test->gsdor);
  GStrv header_ctrl = g_strsplit ("zp zp_sigma", " ", -1);

  g_assert_cmpstrv (header, header_ctrl);

  g_strfreev (header);
  g_strfreev (header_ctrl);
}


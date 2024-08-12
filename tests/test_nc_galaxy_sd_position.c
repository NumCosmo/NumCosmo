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
  NcGalaxySDPosition *gsdp;
} TestNcGalaxySDPosition;

static void test_nc_galaxy_sd_position_flat_new (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_lsst_srd_new (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_lsst_srd_new_y10 (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_free (TestNcGalaxySDPosition *test, gconstpointer pdata);

static void test_nc_galaxy_sd_position_set_get_lim (TestNcGalaxySDPosition *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_gen_lim (TestNcGalaxySDPosition *test, gconstpointer pdata);

static void test_nc_galaxy_sd_position_flat_gen_dist (TestNcGalaxySDPosition *test, gconstpointer pdata);

static void test_nc_galaxy_sd_position_lsst_srd_gen_dist (TestNcGalaxySDPosition *test, gconstpointer pdata);

static void test_nc_galaxy_sd_position_model_id (TestNcGalaxySDPosition *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_position_flat/set_lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_set_get_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/gen_lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_gen_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/gen_dist", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_gen_dist,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_flat/model_id", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_model_id,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd/set_lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new,
              &test_nc_galaxy_sd_position_set_get_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd/gen_lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new,
              &test_nc_galaxy_sd_position_gen_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd/gen_dist", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new,
              &test_nc_galaxy_sd_position_lsst_srd_gen_dist,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd/model_id", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new,
              &test_nc_galaxy_sd_position_model_id,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd_y10/set_lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new_y10,
              &test_nc_galaxy_sd_position_set_get_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd_y10/gen_lim", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new_y10,
              &test_nc_galaxy_sd_position_gen_lim,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd_y10/gen_dist", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new_y10,
              &test_nc_galaxy_sd_position_lsst_srd_gen_dist,
              &test_nc_galaxy_sd_position_free);

  g_test_add ("/nc/galaxy_sd_position_lsst_srd_y10/model_id", TestNcGalaxySDPosition, NULL,
              &test_nc_galaxy_sd_position_lsst_srd_new_y10,
              &test_nc_galaxy_sd_position_model_id,
              &test_nc_galaxy_sd_position_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_sd_position_flat_new (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  const gdouble z_ll               = g_test_rand_double_range (0.1, 0.5);
  const gdouble z_ul               = g_test_rand_double_range (5, 1000);
  const gdouble r_ll               = g_test_rand_double_range (0.001, 0.1);
  const gdouble r_ul               = g_test_rand_double_range (1, 10);
  NcGalaxySDPositionFlat *gsdpflat = nc_galaxy_sd_position_flat_new (z_ll, z_ul, r_ll, r_ul);

  test->gsdp = NC_GALAXY_SD_POSITION (gsdpflat);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpflat));
}

static void
test_nc_galaxy_sd_position_lsst_srd_new (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  const gdouble z_ll                  = g_test_rand_double_range (0.1, 0.5);
  const gdouble z_ul                  = g_test_rand_double_range (5, 1000);
  const gdouble r_ll                  = g_test_rand_double_range (0.001, 0.1);
  const gdouble r_ul                  = g_test_rand_double_range (1, 10);
  NcGalaxySDPositionLSSTSRD *gsdplsst = nc_galaxy_sd_position_lsst_srd_new (z_ll, z_ul, r_ll, r_ul);

  test->gsdp = NC_GALAXY_SD_POSITION (gsdplsst);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_LSST_SRD (gsdplsst));
}

static void
test_nc_galaxy_sd_position_lsst_srd_new_y10 (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  const gdouble z_ll                     = g_test_rand_double_range (0.1, 0.5);
  const gdouble z_ul                     = g_test_rand_double_range (5, 1000);
  const gdouble r_ll                     = g_test_rand_double_range (0.001, 0.1);
  const gdouble r_ul                     = g_test_rand_double_range (1, 10);
  NcGalaxySDPositionLSSTSRD *gsdplssty10 = nc_galaxy_sd_position_lsst_srd_new_y10 (z_ll, z_ul, r_ll, r_ul);

  test->gsdp = NC_GALAXY_SD_POSITION (gsdplssty10);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_LSST_SRD (gsdplssty10));
}

static void
test_nc_galaxy_sd_position_free (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_position_free, test->gsdp);
}

static void
test_nc_galaxy_sd_position_set_get_lim (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  gint i;

  for (i = 0; i < 10000; i++)
  {
    const gdouble z_ll = g_test_rand_double_range (0.1, 0.5);
    const gdouble z_ul = g_test_rand_double_range (5, 1000);
    const gdouble r_ll = g_test_rand_double_range (0.001, 0.1);
    const gdouble r_ul = g_test_rand_double_range (1, 10);
    gdouble z_ll_peek, z_ul_peek;
    gdouble r_ll_peek, r_ul_peek;
    gint j;

    nc_galaxy_sd_position_set_z_lim (test->gsdp, z_ll, z_ul);
    nc_galaxy_sd_position_set_r_lim (test->gsdp, r_ll, r_ul);

    nc_galaxy_sd_position_get_z_lim (test->gsdp, &z_ll_peek, &z_ul_peek);
    nc_galaxy_sd_position_get_r_lim (test->gsdp, &r_ll_peek, &r_ul_peek);

    g_assert_cmpint (z_ll, ==, z_ll_peek);
    g_assert_cmpint (z_ul, ==, z_ul_peek);
  }
}

static void
test_nc_galaxy_sd_position_gen_lim (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const gint nlims = 1000;
  const gint nruns = 1000;
  gint i;

  for (i = 0; i < nlims; i++)
  {
    const gdouble z_ll = g_test_rand_double_range (0.1, 0.5);
    const gdouble z_ul = g_test_rand_double_range (5, 1000);
    const gdouble r_ll = g_test_rand_double_range (0.001, 0.1);
    const gdouble r_ul = g_test_rand_double_range (1, 10);
    gdouble z_ll_peek, z_ul_peek;
    gdouble r_ll_peek, r_ul_peek;
    gint j;

    nc_galaxy_sd_position_set_z_lim (test->gsdp, z_ll, z_ul);
    nc_galaxy_sd_position_set_r_lim (test->gsdp, r_ll, r_ul);

    nc_galaxy_sd_position_get_z_lim (test->gsdp, &z_ll_peek, &z_ul_peek);
    nc_galaxy_sd_position_get_r_lim (test->gsdp, &r_ll_peek, &r_ul_peek);

    for (j = 0; j < nruns; j++)
    {
      const gdouble gen_z = nc_galaxy_sd_position_gen_z (test->gsdp, rng);
      const gdouble gen_r = nc_galaxy_sd_position_gen_r (test->gsdp, rng);

      g_assert_cmpfloat (gen_z, >, z_ll_peek);
      g_assert_cmpfloat (gen_z, <, z_ul_peek);
      g_assert_cmpfloat (gen_r, >, r_ll_peek);
      g_assert_cmpfloat (gen_r, <, r_ul_peek);
    }
  }
}

static void
test_nc_galaxy_sd_position_flat_gen_dist (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcGalaxySDPositionFlat *gsdpflat = NC_GALAXY_SD_POSITION_FLAT (test->gsdp);
  NcmRNG *rng                      = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *z_lim                 = ncm_vector_new (2);
  NcmVector *r_lim                 = ncm_vector_new (2);
  const gdouble z_ll               = g_test_rand_double_range (0.1, 0.5);
  const gdouble z_ul               = g_test_rand_double_range (5.0, 1000.0);
  const gdouble r_ll               = g_test_rand_double_range (0.001, 0.1);
  const gdouble r_ul               = g_test_rand_double_range (1.0, 10.0);
  const gdouble z_avg              = (z_ul + z_ll) / 2.0;
  const gdouble r_avg              = 2.0 / 3.0 * (r_ll + r_ul - r_ll * r_ul / (r_ll + r_ul));
  const gdouble z_var              = 1.0 / 12.0 * (gsl_pow_2 (z_ul - z_ll));
  const gdouble r_var              = 1.0 / 18.0 * (gsl_pow_2 (r_ul - r_ll) * (gsl_pow_2 (r_ul + r_ll) + 2.0 * r_ul * r_ll) / gsl_pow_2 (r_ll + r_ul));
  const gint nruns                 = 10;
  const gint ndata                 = 10000;
  gint i;

  nc_galaxy_sd_position_set_z_lim (test->gsdp, z_ll, z_ul);
  nc_galaxy_sd_position_set_r_lim (test->gsdp, r_ll, r_ul);

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    gint j;

    for (j = 0; j < ndata; j++)
    {
      const gdouble gen_z = nc_galaxy_sd_position_gen_z (test->gsdp, rng);
      const gdouble gen_r = nc_galaxy_sd_position_gen_r (test->gsdp, rng);

      ncm_stats_vec_set (pos_sample, 0, gen_z);
      ncm_stats_vec_set (pos_sample, 1, gen_r);

      ncm_stats_vec_update (pos_sample);
    }

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 5.0 * sqrt (z_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 5.0 * sqrt (z_var / ndata));

    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), <, r_avg + 5.0 * sqrt (r_var / ndata));
    g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), >, r_avg - 5.0 * sqrt (r_var / ndata));

    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 1.0e-1);
    g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 1) / r_var - 1.0), <, 1.0e-1);

    ncm_stats_vec_free (pos_sample);
  }

  ncm_vector_free (z_lim);
  ncm_vector_free (r_lim);
  ncm_rng_free (rng);
}

static gdouble
_test_get_z_avr (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  gdouble z_lb, z_ub;

  nc_galaxy_sd_position_get_z_lim (NC_GALAXY_SD_POSITION (gsdplsst), &z_lb, &z_ub);
  {
    const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdplsst), "alpha");
    const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdplsst), "beta");
    const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdplsst), "z0");
    const gdouble y_lb  = pow (z_lb / z0, alpha);
    const gdouble y_ub  = pow (z_ub / z0, alpha);

    return z0 *
           (gsl_sf_gamma_inc_Q ((2.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((2.0 + beta) / alpha, y_ub)) /
           (gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_ub)) *
           gsl_sf_poch ((1.0 + beta) / alpha, 1.0 / alpha);
  }
}

static gdouble
_test_get_z_var (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  gdouble z_lb, z_ub;

  nc_galaxy_sd_position_get_z_lim (NC_GALAXY_SD_POSITION (gsdplsst), &z_lb, &z_ub);
  {
    const gdouble alpha = ncm_model_param_get_by_name (NCM_MODEL (gsdplsst), "alpha");
    const gdouble beta  = ncm_model_param_get_by_name (NCM_MODEL (gsdplsst), "beta");
    const gdouble z0    = ncm_model_param_get_by_name (NCM_MODEL (gsdplsst), "z0");
    const gdouble y_lb  = pow (z_lb / z0, alpha);
    const gdouble y_ub  = pow (z_ub / z0, alpha);

    return z0 * z0 *
           (gsl_sf_gamma_inc_Q ((3.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((3.0 + beta) / alpha, y_ub)) /
           (gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_lb) - gsl_sf_gamma_inc_Q ((1.0 + beta) / alpha, y_ub)) *
           gsl_sf_poch ((1.0 + beta) / alpha, 2.0 / alpha) -
           gsl_pow_2 (_test_get_z_avr (gsdplsst));
  }
}

static void
test_nc_galaxy_sd_position_lsst_srd_gen_dist (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst = NC_GALAXY_SD_POSITION_LSST_SRD (test->gsdp);
  NcmRNG *rng                         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const gdouble z_ll                  = g_test_rand_double_range (0.1, 0.5);
  const gdouble z_ul                  = g_test_rand_double_range (5.0, 1000.0);
  const gdouble r_ll                  = g_test_rand_double_range (0.001, 0.1);
  const gdouble r_ul                  = g_test_rand_double_range (1.0, 10.0);

  nc_galaxy_sd_position_set_z_lim (test->gsdp, z_ll, z_ul);
  nc_galaxy_sd_position_set_r_lim (test->gsdp, r_ll, r_ul);

  {
    const gdouble z_avg = _test_get_z_avr (gsdplsst);
    const gdouble r_avg = 2.0 / 3.0 * (r_ll + r_ul - r_ll * r_ul / (r_ll + r_ul));
    const gdouble z_var = _test_get_z_var (gsdplsst);
    const gdouble r_var = 1.0 / 18.0 * (gsl_pow_2 (r_ul - r_ll) * (gsl_pow_2 (r_ul + r_ll) + 2.0 * r_ul * r_ll) / gsl_pow_2 (r_ll + r_ul));
    const gint nruns    = 10;
    const gint ndata    = 10000;
    gint i;

    for (i = 0; i < nruns; i++)
    {
      NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
      gint j;

      for (j = 0; j < ndata; j++)
      {
        const gdouble gen_z = nc_galaxy_sd_position_gen_z (test->gsdp, rng);
        const gdouble gen_r = nc_galaxy_sd_position_gen_r (test->gsdp, rng);

        ncm_stats_vec_set (pos_sample, 0, gen_z);
        ncm_stats_vec_set (pos_sample, 1, gen_r);

        ncm_stats_vec_update (pos_sample);
      }

      g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), <, z_avg + 5.0 * sqrt (z_var / ndata));
      g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 0), >, z_avg - 5.0 * sqrt (z_var / ndata));
      g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), <, r_avg + 5.0 * sqrt (r_var / ndata));
      g_assert_cmpfloat (ncm_stats_vec_get_mean (pos_sample, 1), >, r_avg - 5.0 * sqrt (r_var / ndata));

      g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 0) / z_var - 1.0), <, 2.0e-1);
      g_assert_cmpfloat (fabs (ncm_stats_vec_get_var (pos_sample, 1) / r_var - 1.0), <, 1.0e-1);

      ncm_stats_vec_free (pos_sample);
    }
  }

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_sd_position_model_id (TestNcGalaxySDPosition *test, gconstpointer pdata)
{
  NcGalaxySDPosition *gsdp = test->gsdp;
  NcmMSet *model_set       = ncm_mset_empty_new ();
  NcmSerialize *ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  ncm_mset_set (model_set, ncm_model_dup (NCM_MODEL (gsdp), ser));

  g_assert_true (NC_IS_GALAXY_SD_POSITION (ncm_mset_peek (model_set, nc_galaxy_sd_position_id ())));
}


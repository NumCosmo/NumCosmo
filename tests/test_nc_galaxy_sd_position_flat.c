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


typedef struct _TestNcGalaxySDPositionFlat
{
  NcGalaxySDPositionFlat *gsdpflat;
} TestNcGalaxySDPositionFlat;

static void test_nc_galaxy_sd_position_flat_new (TestNcGalaxySDPositionFlat *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_flat_set_lim (TestNcGalaxySDPositionFlat *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_flat_gen_lim (TestNcGalaxySDPositionFlat *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_flat_gen_dist (TestNcGalaxySDPositionFlat *test, gconstpointer pdata);
static void test_nc_galaxy_sd_position_flat_free (TestNcGalaxySDPositionFlat *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

  g_test_add ("/nc/galaxy_sd_position_flat/set_lim", TestNcGalaxySDPositionFlat, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_set_lim,
              &test_nc_galaxy_sd_position_flat_free);

  g_test_add ("/nc/galaxy_sd_position_flat/gen_lim", TestNcGalaxySDPositionFlat, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_gen_lim,
              &test_nc_galaxy_sd_position_flat_free);

  g_test_add ("/nc/galaxy_sd_position_flat/gen_dist", TestNcGalaxySDPositionFlat, NULL,
              &test_nc_galaxy_sd_position_flat_new,
              &test_nc_galaxy_sd_position_flat_gen_dist,
              &test_nc_galaxy_sd_position_flat_free);

  g_test_run ();
}

static void
test_nc_galaxy_sd_position_flat_new (TestNcGalaxySDPositionFlat *test, gconstpointer pdata)
{
  NcGalaxySDPositionFlat *gsdpflat = nc_galaxy_sd_position_flat_new ();

  test->gsdpflat = nc_galaxy_sd_position_flat_new ();

  g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpflat));
}

static void
test_nc_galaxy_sd_position_flat_free (TestNcGalaxySDPositionFlat *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_sd_position_flat_free, test->gsdpflat);
}

static void
test_nc_galaxy_sd_position_flat_set_lim (TestNcGalaxySDPositionFlat *test, gconstpointer pdata)
{
  gint i;

  for (i = 0; i < 10000; i++)
  {
    gint j;
    NcmVector *z_lim   = ncm_vector_new (2);
    NcmVector *r_lim   = ncm_vector_new (2);
    const gdouble z_ll = g_test_rand_double_range (0.1, 0.5);
    const gdouble z_ul = g_test_rand_double_range (5, 1000);
    const gdouble r_ll = g_test_rand_double_range (0.001, 0.1);
    const gdouble r_ul = g_test_rand_double_range (1, 10);

    ncm_vector_set (z_lim, 0, z_ll);
    ncm_vector_set (z_lim, 1, z_ul);
    ncm_vector_set (r_lim, 0, r_ll);
    ncm_vector_set (r_lim, 1, r_ul);

    nc_galaxy_sd_position_flat_set_z_lim (test->gsdpflat, z_lim);
    nc_galaxy_sd_position_flat_set_r_lim (test->gsdpflat, r_lim);

    NcmVector *z_lim_peek = nc_galaxy_sd_position_flat_peek_z_lim (test->gsdpflat);
    NcmVector *r_lim_peek = nc_galaxy_sd_position_flat_peek_r_lim (test->gsdpflat);

    g_assert_cmpint (ncm_vector_len (z_lim), ==, ncm_vector_len (z_lim_peek));
    g_assert_cmpint (ncm_vector_len (r_lim), ==, ncm_vector_len (r_lim_peek));

    for (j = 0; j < 2; j++)
    {
      g_assert_cmpfloat (ncm_vector_get (z_lim_peek, j), ==, ncm_vector_get (z_lim, j));
      g_assert_cmpfloat (ncm_vector_get (r_lim_peek, j), ==, ncm_vector_get (r_lim, j));
    }
  }
}

static void
test_nc_galaxy_sd_position_flat_gen_lim (TestNcGalaxySDPositionFlat *test, gconstpointer pdata)
{
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const gint nlims = 1000;
  const gint nruns = 1000;
  gint i;

  for (i = 0; i < nlims; i++)
  {
    gint j;
    NcmVector *z_lim   = ncm_vector_new (2);
    NcmVector *r_lim   = ncm_vector_new (2);
    const gdouble z_ll = g_test_rand_double_range (0.1, 0.5);
    const gdouble z_ul = g_test_rand_double_range (5, 1000);
    const gdouble r_ll = g_test_rand_double_range (0.001, 0.1);
    const gdouble r_ul = g_test_rand_double_range (1, 10);

    ncm_vector_set (z_lim, 0, z_ll);
    ncm_vector_set (z_lim, 1, z_ul);
    ncm_vector_set (r_lim, 0, r_ll);
    ncm_vector_set (r_lim, 1, r_ul);

    nc_galaxy_sd_position_flat_set_z_lim (test->gsdpflat, z_lim);
    nc_galaxy_sd_position_flat_set_r_lim (test->gsdpflat, r_lim);

    NcmVector *z_lim_peek = nc_galaxy_sd_position_flat_peek_z_lim (test->gsdpflat);
    NcmVector *r_lim_peek = nc_galaxy_sd_position_flat_peek_r_lim (test->gsdpflat);

    for (j = 0; j < nruns; j++)
    {
      const gdouble gen_z = nc_galaxy_sd_position_gen_z (NC_GALAXY_SD_POSITION (test->gsdpflat), rng);
      const gdouble gen_r = nc_galaxy_sd_position_gen_r (NC_GALAXY_SD_POSITION (test->gsdpflat), rng);

      g_assert_cmpfloat (gen_z, >, ncm_vector_get (z_lim_peek, 0));
      g_assert_cmpfloat (gen_z, <, ncm_vector_get (z_lim_peek, 1));
      g_assert_cmpfloat (gen_r, >, ncm_vector_get (r_lim_peek, 0));
      g_assert_cmpfloat (gen_r, <, ncm_vector_get (r_lim_peek, 1));
    }
  }
}

static void
test_nc_galaxy_sd_position_flat_gen_dist (TestNcGalaxySDPositionFlat *test, gconstpointer pdata)
{
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmVector *z_lim    = ncm_vector_new (2);
  NcmVector *r_lim    = ncm_vector_new (2);
  const gdouble z_ll  = g_test_rand_double_range (0.1, 0.5);
  const gdouble z_ul  = g_test_rand_double_range (5.0, 1000.0);
  const gdouble r_ll  = g_test_rand_double_range (0.001, 0.1);
  const gdouble r_ul  = g_test_rand_double_range (1.0, 10.0);
  const gdouble z_avg = (z_ul + z_ll) / 2.0;
  const gdouble r_avg = 2.0 / 3.0 * (r_ll + r_ul - r_ll * r_ul / (r_ll + r_ul));
  const gdouble z_var = 1.0 / 12.0 * (gsl_pow_2 (z_ul - z_ll));
  const gdouble r_var = 1.0 / 18.0 * (gsl_pow_2 (r_ul - r_ll) * (gsl_pow_2 (r_ul + r_ll) + 2.0 * r_ul * r_ll) / gsl_pow_2 (r_ll + r_ul));
  const gint nruns    = 10;
  const gint ndata    = 10000;
  gint i;

  ncm_vector_set (z_lim, 0, z_ll);
  ncm_vector_set (z_lim, 1, z_ul);
  ncm_vector_set (r_lim, 0, r_ll);
  ncm_vector_set (r_lim, 1, r_ul);

  nc_galaxy_sd_position_flat_set_z_lim (test->gsdpflat, z_lim);
  nc_galaxy_sd_position_flat_set_r_lim (test->gsdpflat, r_lim);

  for (i = 0; i < nruns; i++)
  {
    NcmStatsVec *pos_sample = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    gint j;

    for (j = 0; j < ndata; j++)
    {
      const gdouble gen_z = nc_galaxy_sd_position_gen_z (NC_GALAXY_SD_POSITION (test->gsdpflat), rng);
      const gdouble gen_r = nc_galaxy_sd_position_gen_r (NC_GALAXY_SD_POSITION (test->gsdpflat), rng);

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


/***************************************************************************
 *            test_nc_galaxy_wl_obs.c
 *
 *  Tue Jul 30 23:54:39 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiooliveiraCode@proton.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiooliveiraCode@proton.me>
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

typedef struct _TestNcGalaxyWLObs
{
  NcGalaxyWLObs *obs;
} TestNcGalaxyWLObs;

static void test_nc_galaxy_wl_obs_new (TestNcGalaxyWLObs *fixture, gconstpointer user_data);
static void test_nc_galaxy_wl_obs_new_pz (TestNcGalaxyWLObs *fixture, gconstpointer user_data);
static void test_nc_galaxy_wl_obs_free (TestNcGalaxyWLObs *fixture, gconstpointer user_data);

static void test_nc_galaxy_wl_obs_data (TestNcGalaxyWLObs *fixture, gconstpointer user_data);
static void test_nc_galaxy_wl_obs_data_pz (TestNcGalaxyWLObs *fixture, gconstpointer user_data);
static void test_nc_galaxy_wl_obs_header (TestNcGalaxyWLObs *fixture, gconstpointer user_data);
static void test_nc_galaxy_wl_obs_header_pz (TestNcGalaxyWLObs *fixture, gconstpointer user_data);
static void test_nc_galaxy_wl_obs_coord (TestNcGalaxyWLObs *fixture, gconstpointer user_data);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/data", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_data,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/data_pz", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new_pz,
              &test_nc_galaxy_wl_obs_data_pz,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/header", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_header,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/header_pz", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new_pz,
              &test_nc_galaxy_wl_obs_header_pz,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/coord", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_coord,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/coord_pz", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new_pz,
              &test_nc_galaxy_wl_obs_coord,
              &test_nc_galaxy_wl_obs_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_wl_obs_new (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcGalaxyWLObsCoord coord = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  guint nrows              = 100;
  GStrv col_names          = g_strsplit ("ra dec z", " ", -1);

  test->obs = nc_galaxy_wl_obs_new (coord, nrows, col_names);
  g_assert_nonnull (test->obs);
  g_assert_true (NC_IS_GALAXY_WL_OBS (test->obs));
}

static void
test_nc_galaxy_wl_obs_new_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcGalaxyWLObsCoord coord = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  guint nrows              = 100;
  GStrv col_names          = g_strsplit ("ra dec pz", " ", -1);

  test->obs = nc_galaxy_wl_obs_new (coord, nrows, col_names);
  g_assert_nonnull (test->obs);
  g_assert_true (NC_IS_GALAXY_WL_OBS (test->obs));
}

static void
test_nc_galaxy_wl_obs_free (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_wl_obs_free, test->obs);
}

static void
test_nc_galaxy_wl_obs_data (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint i;
  gdouble val;

  for (i = 0; i < nc_galaxy_wl_obs_len (test->obs); i++)
  {
    val = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i), ==, val);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i), ==, val);

    val = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    nc_galaxy_wl_obs_set (test->obs, "z", i, val);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "z", i), ==, val);
  }
}

static void
test_nc_galaxy_wl_obs_data_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint i;

  for (i = 0; i < nc_galaxy_wl_obs_len (test->obs); i++)
  {
    NcmSplineCubicNotaknot *nak = ncm_spline_cubic_notaknot_new ();
    NcmVector *xv               = ncm_vector_new (100);
    NcmVector *yv               = ncm_vector_new (100);
    NcmSpline *pz;
    gdouble val;
    guint j;

    for (j = 0; j < 100; j++)
    {
      ncm_vector_set (xv, j, (gdouble) j);
      ncm_vector_set (yv, j, (gdouble) j * j);
    }

    pz = ncm_spline_new (NCM_SPLINE (nak), xv, yv, TRUE);

    val = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i), ==, val);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i), ==, val);

    val = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    nc_galaxy_wl_obs_set (test->obs, "z", i, val);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "z", i), ==, val);

    nc_galaxy_wl_obs_set_pz (test->obs, i, pz);
    g_assert_true (NCM_IS_SPLINE (nc_galaxy_wl_obs_peek_pz (test->obs, i)));

    for (j = 0; j < 100; j++)
    {
      g_assert_cmpfloat_with_epsilon (ncm_spline_eval (nc_galaxy_wl_obs_peek_pz (test->obs, i), j), j * j, 1e-10);
    }
  }
}

static void
test_nc_galaxy_wl_obs_header (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmVarDict *header = nc_galaxy_wl_obs_peek_header (test->obs);
  gboolean a;
  gint i;

  g_assert_nonnull (header);

  g_assert_true (ncm_var_dict_get_int (header, "ra", &i));
  g_assert_cmpint (i, ==, 0);

  g_assert_true (ncm_var_dict_get_int (header, "dec", &i));
  g_assert_cmpint (i, ==, 1);

  g_assert_true (ncm_var_dict_get_int (header, "z", &i));
  g_assert_cmpint (i, ==, 2);

  g_assert_false (ncm_var_dict_get_boolean (header, "pz", &a));
}

static void
test_nc_galaxy_wl_obs_header_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmVarDict *header = nc_galaxy_wl_obs_peek_header (test->obs);
  gboolean a;
  gint i;

  g_assert_nonnull (header);

  g_assert_true (ncm_var_dict_get_int (header, "ra", &i));
  g_assert_cmpint (i, ==, 0);

  g_assert_true (ncm_var_dict_get_int (header, "dec", &i));
  g_assert_cmpint (i, ==, 1);

  g_assert_true (ncm_var_dict_get_boolean (header, "pz", &a));
  g_assert_true (a);
}

static void
test_nc_galaxy_wl_obs_coord (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcGalaxyWLObsCoord coord = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;

  g_assert_cmpint (nc_galaxy_wl_obs_get_coord (test->obs), ==, NC_GALAXY_WL_OBS_COORD_CELESTIAL);
  nc_galaxy_wl_obs_set_coord (test->obs, coord);
  g_assert_cmpint (nc_galaxy_wl_obs_get_coord (test->obs), ==, coord);
}


/***************************************************************************
 *            test_nc_galaxy_wl_obs.c
 *
 *  Tue Jul 30 23:54:39 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiolimadeoliveira@pm.me>
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcGalaxyWLObs
{
  NcGalaxyWLObs *obs;
} TestNcGalaxyWLObs;

static void test_nc_galaxy_wl_obs_new (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_new_pz (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_free (TestNcGalaxyWLObs *test, gconstpointer pdata);

static void test_nc_galaxy_wl_obs_serialize (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_serialize_pz (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_setget (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_data (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_data_pz (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_header (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_header_pz (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_coord (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_bad_set (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_bad_set_pz (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_bad_get_pz (TestNcGalaxyWLObs *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/serialize", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_serialize,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/serialize_pz", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new_pz,
              &test_nc_galaxy_wl_obs_serialize_pz,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/setget", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_setget,
              &test_nc_galaxy_wl_obs_free);

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

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/bad_set", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_bad_set,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/bad_set_pz", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_bad_set_pz,
              &test_nc_galaxy_wl_obs_free);

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/bad_get_pz", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_bad_get_pz,
              &test_nc_galaxy_wl_obs_free);

  g_test_run ();

  return 0;
}

static void
test_nc_galaxy_wl_obs_new (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcGalaxyWLObsCoord coord = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  GStrv col_names          = g_strsplit ("ra dec z", " ", -1);
  guint nrows              = 100;

  test->obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, coord, nrows, col_names);
  g_assert_nonnull (test->obs);
  g_assert_true (NC_IS_GALAXY_WL_OBS (test->obs));

  g_strfreev (col_names);
}

static void
test_nc_galaxy_wl_obs_new_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcGalaxyWLObsCoord coord = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  GStrv col_names          = g_strsplit ("ra dec pz", " ", -1);
  guint nrows              = 100;

  test->obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, coord, nrows, col_names);
  g_assert_nonnull (test->obs);
  g_assert_true (NC_IS_GALAXY_WL_OBS (test->obs));

  g_strfreev (col_names);
}

static void
test_nc_galaxy_wl_obs_free (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_wl_obs_free, test->obs);
}

static void
test_nc_galaxy_wl_obs_serialize (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxyWLObs *obs_dup;
  NcmSerialize *ser;
  GVariant *obs_ser;
  gdouble val;
  guint i;
  guint j;
  guint k;

  for (i = 0; i < nc_galaxy_wl_obs_len (test->obs); i++)
  {
    val = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val);

    val = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    nc_galaxy_wl_obs_set (test->obs, "z", i, val);
  }

  ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->obs));
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));

  g_assert_nonnull (obs_dup);
  g_assert_true (NC_IS_GALAXY_WL_OBS (obs_dup));

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "ra", &j));
  g_assert_true (nc_galaxy_wl_obs_get_index (obs_dup, "ra", &k));
  g_assert_cmpint (j, ==, k);

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "dec", &j));
  g_assert_true (nc_galaxy_wl_obs_get_index (obs_dup, "dec", &k));
  g_assert_cmpint (j, ==, k);

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "z", &j));
  g_assert_true (nc_galaxy_wl_obs_get_index (obs_dup, "z", &k));
  g_assert_cmpint (j, ==, k);

  for (i = 0; i < nc_galaxy_wl_obs_len (test->obs); i++)
  {
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i), ==, nc_galaxy_wl_obs_get (obs_dup, "ra", i));
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i), ==, nc_galaxy_wl_obs_get (obs_dup, "dec", i));
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "z", i), ==, nc_galaxy_wl_obs_get (obs_dup, "z", i));
  }

  g_assert_cmpint (nc_galaxy_wl_obs_get_coord (test->obs), ==, nc_galaxy_wl_obs_get_coord (obs_dup));

  g_variant_unref (obs_ser);
  nc_galaxy_wl_obs_free (obs_dup);
  ncm_serialize_free (ser);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_serialize_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcGalaxyWLObs *obs_dup;
  NcmSerialize *ser;
  GVariant *obs_ser;
  guint i;
  guint j;
  guint k;

  for (i = 0; i < nc_galaxy_wl_obs_len (test->obs); i++)
  {
    NcmSplineCubicNotaknot *nak = ncm_spline_cubic_notaknot_new ();
    NcmVector *xv               = ncm_vector_new (100);
    NcmVector *yv               = ncm_vector_new (100);
    NcmSpline *pz;
    gdouble val;
    guint l;

    for (l = 0; l < 100; l++)
    {
      ncm_vector_set (xv, l, (gdouble) l);
      ncm_vector_set (yv, l, (gdouble) l * l);
    }

    pz = ncm_spline_new (NCM_SPLINE (nak), xv, yv, TRUE);

    val = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val);

    nc_galaxy_wl_obs_set_pz (test->obs, i, pz);

    ncm_spline_free (NCM_SPLINE (nak));
    ncm_spline_free (pz);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }

  ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->obs));
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));

  g_assert_nonnull (obs_dup);
  g_assert_true (NC_IS_GALAXY_WL_OBS (obs_dup));

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "ra", &j));
  g_assert_true (nc_galaxy_wl_obs_get_index (obs_dup, "ra", &k));
  g_assert_cmpint (j, ==, k);

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "dec", &j));
  g_assert_true (nc_galaxy_wl_obs_get_index (obs_dup, "dec", &k));
  g_assert_cmpint (j, ==, k);

  g_assert_cmpint (nc_galaxy_wl_obs_get_coord (test->obs), ==, nc_galaxy_wl_obs_get_coord (obs_dup));

  for (i = 0; i < nc_galaxy_wl_obs_len (test->obs); i++)
  {
    guint l;

    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i), ==, nc_galaxy_wl_obs_get (obs_dup, "ra", i));
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i), ==, nc_galaxy_wl_obs_get (obs_dup, "dec", i));

    for (l = 0; l < 100; l++)
    {
      g_assert_cmpfloat (ncm_spline_eval (nc_galaxy_wl_obs_peek_pz (test->obs, i), l), ==, ncm_spline_eval (nc_galaxy_wl_obs_peek_pz (obs_dup, i), l));
    }
  }

  g_variant_unref (obs_ser);
  ncm_serialize_free (ser);
  ncm_rng_free (rng);
  nc_galaxy_wl_obs_free (obs_dup);
}

static void
test_nc_galaxy_wl_obs_setget (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint ngals        = 200;
  GStrv header       = g_strsplit ("ra dec z e1 e2 e1_int e2_int e_rms e_sigma", " ", -1);
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObs *obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_GALAXY_WL_OBS_COORD_EUCLIDEAN, ngals, header);
  NcGalaxyWLObs *obs2;
  guint i;

  for (i = 0; i < ngals; i++)
  {
    gdouble ra      = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    gdouble dec     = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    gdouble z       = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    gdouble e1      = ncm_rng_uniform_gen (rng, -1.0, 1.0);
    gdouble e2      = ncm_rng_uniform_gen (rng, -1.0, 1.0);
    gdouble e_rms   = ncm_rng_uniform_gen (rng, 0.0, 0.3);
    gdouble e_sigma = ncm_rng_uniform_gen (rng, 0.0, 0.1);

    nc_galaxy_wl_obs_set (obs, "ra", i, ra);
    nc_galaxy_wl_obs_set (obs, "dec", i, dec);
    nc_galaxy_wl_obs_set (obs, "z", i, z);
    nc_galaxy_wl_obs_set (obs, "e1", i, e1);
    nc_galaxy_wl_obs_set (obs, "e2", i, e2);
    nc_galaxy_wl_obs_set (obs, "e1_int", i, e1);
    nc_galaxy_wl_obs_set (obs, "e2_int", i, e2);
    nc_galaxy_wl_obs_set (obs, "e_rms", i, e_rms);
    nc_galaxy_wl_obs_set (obs, "e_sigma", i, e_sigma);
  }

  obs2 = NC_GALAXY_WL_OBS (ncm_serialize_dup_obj (ser, G_OBJECT (obs)));

  for (i = 0; i < ngals; i++)
  {
    gdouble ra      = nc_galaxy_wl_obs_get (obs2, "ra", i);
    gdouble dec     = nc_galaxy_wl_obs_get (obs2, "dec", i);
    gdouble z       = nc_galaxy_wl_obs_get (obs2, "z", i);
    gdouble e1      = nc_galaxy_wl_obs_get (obs2, "e1", i);
    gdouble e2      = nc_galaxy_wl_obs_get (obs2, "e2", i);
    gdouble e_rms   = nc_galaxy_wl_obs_get (obs2, "e_rms", i);
    gdouble e_sigma = nc_galaxy_wl_obs_get (obs2, "e_sigma", i);

    g_assert_cmpfloat (ra, ==, nc_galaxy_wl_obs_get (obs, "ra", i));
    g_assert_cmpfloat (dec, ==, nc_galaxy_wl_obs_get (obs, "dec", i));
    g_assert_cmpfloat (z, ==, nc_galaxy_wl_obs_get (obs, "z", i));
    g_assert_cmpfloat (e1, ==, nc_galaxy_wl_obs_get (obs, "e1", i));
    g_assert_cmpfloat (e2, ==, nc_galaxy_wl_obs_get (obs, "e2", i));
    g_assert_cmpfloat (e_rms, ==, nc_galaxy_wl_obs_get (obs, "e_rms", i));
    g_assert_cmpfloat (e_sigma, ==, nc_galaxy_wl_obs_get (obs, "e_sigma", i));
  }

  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_wl_obs_free (obs2);
  ncm_serialize_free (ser);
  g_strfreev (header);
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_data (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  gdouble val;
  guint i;

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

  ncm_rng_free (rng);
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

    nc_galaxy_wl_obs_set_pz (test->obs, i, pz);
    g_assert_true (NCM_IS_SPLINE (nc_galaxy_wl_obs_peek_pz (test->obs, i)));

    for (j = 0; j < 100; j++)
    {
      g_assert_cmpfloat_with_epsilon (ncm_spline_eval (nc_galaxy_wl_obs_peek_pz (test->obs, i), j), (gdouble) j * j, 1e-10);
    }

    ncm_spline_free (NCM_SPLINE (nak));
    ncm_spline_free (pz);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_header (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  guint i;

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "ra", &i));
  g_assert_cmpint (i, ==, 0);

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "dec", &i));
  g_assert_cmpint (i, ==, 1);

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "z", &i));
  g_assert_cmpint (i, ==, 2);

  g_assert_false (nc_galaxy_wl_obs_get_index (test->obs, "pz", &i));
}

static void
test_nc_galaxy_wl_obs_header_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  guint i;

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "ra", &i));
  g_assert_cmpint (i, ==, 0);

  g_assert_true (nc_galaxy_wl_obs_get_index (test->obs, "dec", &i));
  g_assert_cmpint (i, ==, 1);
}

static void
test_nc_galaxy_wl_obs_coord (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcGalaxyWLObsCoord coord = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;

  g_assert_cmpint (nc_galaxy_wl_obs_get_coord (test->obs), ==, NC_GALAXY_WL_OBS_COORD_CELESTIAL);
  nc_galaxy_wl_obs_set_coord (test->obs, coord);
  g_assert_cmpint (nc_galaxy_wl_obs_get_coord (test->obs), ==, coord);
}

static void
test_nc_galaxy_wl_obs_bad_set_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());

  if (g_test_subprocess ())
    nc_galaxy_wl_obs_set (test->obs, "pz", ncm_rng_uniform_gen (rng, 0, 99), 0);

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_bad_set (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());

  if (g_test_subprocess ())
    nc_galaxy_wl_obs_set (test->obs, "a", ncm_rng_uniform_gen (rng, 0, 99), 0);

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_bad_get_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());

  if (g_test_subprocess ())
    nc_galaxy_wl_obs_get (test->obs, "pz", ncm_rng_uniform_gen (rng, 0, 99));

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  ncm_rng_free (rng);
}


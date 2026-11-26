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
#include <glib/gstdio.h>
#include <string.h>

/* The environment variable that gives the abort test of the data hash its catalog file */
#define TEST_NC_GALAXY_WL_OBS_FILE_ENV "NC_GALAXY_WL_OBS_TEST_FILE"

typedef struct _TestNcGalaxyWLObs
{
  NcGalaxyWLObs *obs;
} TestNcGalaxyWLObs;

static void _test_nc_galaxy_wl_obs_assert_column_not_found_error (GError **error, const gchar *col);
static gchar *_test_nc_galaxy_wl_obs_prepare_file (NcGalaxyWLObs *obs, const gchar *dir);

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
static void test_nc_galaxy_wl_obs_shape_interp (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_file (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_new_from_file (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_pz_spline (void);
static void test_nc_galaxy_wl_obs_pz_spline_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_type (void);
static void test_nc_galaxy_wl_obs_shape_interp_type_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_hash (void);
static void test_nc_galaxy_wl_obs_shape_interp_hash_subprocess (void);
static void test_nc_galaxy_wl_obs_cheb_grid (TestNcGalaxyWLObs *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_cheb_grid_errors (void);
static void test_nc_galaxy_wl_obs_cheb_grid_errors_subprocess (void);

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

  g_test_add ("/numcosmo/nc_galaxy_wl_obs/shape_interp", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_shape_interp,
              &test_nc_galaxy_wl_obs_free);
  g_test_add ("/numcosmo/nc_galaxy_wl_obs/shape_interp/file", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_shape_interp_file,
              &test_nc_galaxy_wl_obs_free);
  g_test_add ("/numcosmo/nc_galaxy_wl_obs/new_from_file", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new,
              &test_nc_galaxy_wl_obs_new_from_file,
              &test_nc_galaxy_wl_obs_free);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/pz_spline", &test_nc_galaxy_wl_obs_pz_spline);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/pz_spline/subprocess", &test_nc_galaxy_wl_obs_pz_spline_subprocess);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/shape_interp/type", &test_nc_galaxy_wl_obs_shape_interp_type);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/shape_interp/type/subprocess", &test_nc_galaxy_wl_obs_shape_interp_type_subprocess);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/shape_interp/hash", &test_nc_galaxy_wl_obs_shape_interp_hash);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/shape_interp/hash/subprocess", &test_nc_galaxy_wl_obs_shape_interp_hash_subprocess);
  g_test_add ("/numcosmo/nc_galaxy_wl_obs/cheb_grid", TestNcGalaxyWLObs, NULL,
              &test_nc_galaxy_wl_obs_new_pz,
              &test_nc_galaxy_wl_obs_cheb_grid,
              &test_nc_galaxy_wl_obs_free);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/cheb_grid/errors", &test_nc_galaxy_wl_obs_cheb_grid_errors);
  g_test_add_func ("/numcosmo/nc_galaxy_wl_obs/cheb_grid/errors/subprocess", &test_nc_galaxy_wl_obs_cheb_grid_errors_subprocess);

  g_test_run ();

  return 0;
}

static void
_test_nc_galaxy_wl_obs_assert_column_not_found_error (GError **error, const gchar *col)
{
  gchar *expected_message = g_strdup_printf ("Column '%s' not found.", col);

  g_assert_nonnull (*error);
  g_assert_error (*error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND);
  g_assert_cmpstr ((*error)->message, ==, expected_message);

  g_clear_error (error);
  g_assert_null (*error);
  g_free (expected_message);
}

/*
 * Sets two interpolants and a record on obs, writes the interpolants to dir/tables.gvar and
 * sets on obs, in their place, an object backed by that file. Returns the absolute name of
 * the file.
 */
static gchar *
_test_nc_galaxy_wl_obs_prepare_file (NcGalaxyWLObs *obs, const gchar *dir)
{
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectralTensor *spp_0   = ncm_spectral_tensor_new (2);
  NcmSpectralTensor *spp_3   = ncm_spectral_tensor_new (2);
  NcmObjDictStr *source      = ncm_obj_dict_str_new ();
  NcGalaxyShapePopGauss *pop = nc_galaxy_shape_pop_gauss_new ();
  gchar *filename            = g_build_filename (dir, "tables.gvar", NULL);
  NcGalaxyWLObsShapeInterp *si, *file_si;

  ncm_obj_dict_str_add (source, "shape-pop", G_OBJECT (pop));
  nc_galaxy_wl_obs_set (obs, "z", 0, 0.5, NULL);
  nc_galaxy_wl_obs_set_shape_interp (obs, 0, spp_0);
  nc_galaxy_wl_obs_set_shape_interp (obs, 3, spp_3);
  nc_galaxy_wl_obs_set_shape_interp_source (obs, source);

  g_object_get (obs, "shape-interp", &si, NULL);
  ncm_serialize_to_binfile (ser, G_OBJECT (si), filename);

  file_si = nc_galaxy_wl_obs_shape_interp_new_from_file (filename, nc_galaxy_wl_obs_shape_interp_peek_data_hash (si),
                                                         nc_galaxy_wl_obs_shape_interp_peek_source (si));
  g_object_set (obs, "shape-interp", file_si, NULL);

  nc_galaxy_wl_obs_shape_interp_free (file_si);
  nc_galaxy_wl_obs_shape_interp_free (si);
  nc_galaxy_shape_pop_gauss_free (pop);
  ncm_obj_dict_str_unref (source);
  ncm_spectral_tensor_free (spp_0);
  ncm_spectral_tensor_free (spp_3);
  ncm_serialize_free (ser);

  return filename;
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
  g_assert_true (NCM_IS_CATALOG (test->obs));

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
  g_assert_true (NCM_IS_CATALOG (test->obs));

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
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val, NULL);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val, NULL);

    val = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    nc_galaxy_wl_obs_set (test->obs, "z", i, val, NULL);
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
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i, NULL), ==, nc_galaxy_wl_obs_get (obs_dup, "ra", i, NULL));
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i, NULL), ==, nc_galaxy_wl_obs_get (obs_dup, "dec", i, NULL));
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "z", i, NULL), ==, nc_galaxy_wl_obs_get (obs_dup, "z", i, NULL));
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
    NcmVector *xv = ncm_vector_new (100);
    NcmVector *yv = ncm_vector_new (100);
    NcGalaxyPz *pz;
    gdouble val;
    guint l;

    for (l = 0; l < 100; l++)
    {
      ncm_vector_set (xv, l, (gdouble) l);
      ncm_vector_set (yv, l, (gdouble) l * l);
    }

    pz = NC_GALAXY_PZ (nc_galaxy_pz_bspline_sqrt_new_from_table (xv, yv));

    val = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val, NULL);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val, NULL);

    nc_galaxy_wl_obs_set_pz (test->obs, i, pz);

    nc_galaxy_pz_free (pz);
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

    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i, NULL), ==, nc_galaxy_wl_obs_get (obs_dup, "ra", i, NULL));
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i, NULL), ==, nc_galaxy_wl_obs_get (obs_dup, "dec", i, NULL));

    for (l = 0; l < 100; l++)
    {
      g_assert_cmpfloat (nc_galaxy_pz_eval_p (nc_galaxy_wl_obs_peek_pz (test->obs, i), l), ==, nc_galaxy_pz_eval_p (nc_galaxy_wl_obs_peek_pz (obs_dup, i), l));
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

    nc_galaxy_wl_obs_set (obs, "ra", i, ra, NULL);
    nc_galaxy_wl_obs_set (obs, "dec", i, dec, NULL);
    nc_galaxy_wl_obs_set (obs, "z", i, z, NULL);
    nc_galaxy_wl_obs_set (obs, "e1", i, e1, NULL);
    nc_galaxy_wl_obs_set (obs, "e2", i, e2, NULL);
    nc_galaxy_wl_obs_set (obs, "e1_int", i, e1, NULL);
    nc_galaxy_wl_obs_set (obs, "e2_int", i, e2, NULL);
    nc_galaxy_wl_obs_set (obs, "e_rms", i, e_rms, NULL);
    nc_galaxy_wl_obs_set (obs, "e_sigma", i, e_sigma, NULL);
  }

  obs2 = NC_GALAXY_WL_OBS (ncm_serialize_dup_obj (ser, G_OBJECT (obs)));

  for (i = 0; i < ngals; i++)
  {
    gdouble ra      = nc_galaxy_wl_obs_get (obs2, "ra", i, NULL);
    gdouble dec     = nc_galaxy_wl_obs_get (obs2, "dec", i, NULL);
    gdouble z       = nc_galaxy_wl_obs_get (obs2, "z", i, NULL);
    gdouble e1      = nc_galaxy_wl_obs_get (obs2, "e1", i, NULL);
    gdouble e2      = nc_galaxy_wl_obs_get (obs2, "e2", i, NULL);
    gdouble e_rms   = nc_galaxy_wl_obs_get (obs2, "e_rms", i, NULL);
    gdouble e_sigma = nc_galaxy_wl_obs_get (obs2, "e_sigma", i, NULL);

    g_assert_cmpfloat (ra, ==, nc_galaxy_wl_obs_get (obs, "ra", i, NULL));
    g_assert_cmpfloat (dec, ==, nc_galaxy_wl_obs_get (obs, "dec", i, NULL));
    g_assert_cmpfloat (z, ==, nc_galaxy_wl_obs_get (obs, "z", i, NULL));
    g_assert_cmpfloat (e1, ==, nc_galaxy_wl_obs_get (obs, "e1", i, NULL));
    g_assert_cmpfloat (e2, ==, nc_galaxy_wl_obs_get (obs, "e2", i, NULL));
    g_assert_cmpfloat (e_rms, ==, nc_galaxy_wl_obs_get (obs, "e_rms", i, NULL));
    g_assert_cmpfloat (e_sigma, ==, nc_galaxy_wl_obs_get (obs, "e_sigma", i, NULL));
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
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val, NULL);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i, NULL), ==, val);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val, NULL);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i, NULL), ==, val);

    val = ncm_rng_uniform_gen (rng, 0.0, 5.0);
    nc_galaxy_wl_obs_set (test->obs, "z", i, val, NULL);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "z", i, NULL), ==, val);
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
    NcmVector *xv = ncm_vector_new (100);
    NcmVector *yv = ncm_vector_new (100);
    NcGalaxyPz *pz;
    gdouble val;
    guint j;

    for (j = 0; j < 100; j++)
    {
      ncm_vector_set (xv, j, (gdouble) j);
      ncm_vector_set (yv, j, (gdouble) j * j);
    }

    pz = NC_GALAXY_PZ (nc_galaxy_pz_bspline_sqrt_new_from_table (xv, yv));

    val = ncm_rng_uniform_gen (rng, 0.0, 360.0);
    nc_galaxy_wl_obs_set (test->obs, "ra", i, val, NULL);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "ra", i, NULL), ==, val);

    val = ncm_rng_uniform_gen (rng, -90.0, 90.0);
    nc_galaxy_wl_obs_set (test->obs, "dec", i, val, NULL);
    g_assert_cmpfloat (nc_galaxy_wl_obs_get (test->obs, "dec", i, NULL), ==, val);

    nc_galaxy_wl_obs_set_pz (test->obs, i, pz);
    g_assert_true (NC_IS_GALAXY_PZ (nc_galaxy_wl_obs_peek_pz (test->obs, i)));

    for (j = 0; j < 100; j++)
    {
      g_assert_cmpfloat_with_epsilon (nc_galaxy_pz_eval_p (nc_galaxy_wl_obs_peek_pz (test->obs, i), j), (gdouble) j * j, 1e-10);
    }

    nc_galaxy_pz_free (pz);
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
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  GError *error = NULL;

  nc_galaxy_wl_obs_set (test->obs, "pz", ncm_rng_uniform_gen (rng, 0, 99), 0, &error);
  _test_nc_galaxy_wl_obs_assert_column_not_found_error (&error, "pz");

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_bad_set (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  GError *error = NULL;

  nc_galaxy_wl_obs_set (test->obs, "a", ncm_rng_uniform_gen (rng, 0, 99), 0, &error);
  _test_nc_galaxy_wl_obs_assert_column_not_found_error (&error, "a");

  ncm_rng_free (rng);
}

static void
test_nc_galaxy_wl_obs_bad_get_pz (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  GError *error = NULL;

  g_assert_true (isnan (nc_galaxy_wl_obs_get (test->obs, "pz", ncm_rng_uniform_gen (rng, 0, 99), &error)));
  _test_nc_galaxy_wl_obs_assert_column_not_found_error (&error, "pz");

  ncm_rng_free (rng);
}

/*
 * The interpolants are visible only with a record and while the rows are those of the
 * record; a serialized catalog carries them, and a stale one serializes none.
 */
static void
test_nc_galaxy_wl_obs_shape_interp (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectralTensor *spp_0   = ncm_spectral_tensor_new (2);
  NcmSpectralTensor *spp_3   = ncm_spectral_tensor_new (2);
  NcmObjDictStr *source      = ncm_obj_dict_str_new ();
  NcGalaxyShapePopGauss *pop = nc_galaxy_shape_pop_gauss_new ();
  NcGalaxyWLObsShapeInterp *prop;
  NcGalaxyWLObs *obs_dup;
  GVariant *obs_ser;

  ncm_obj_dict_str_add (source, "shape-pop", G_OBJECT (pop));
  nc_galaxy_wl_obs_set (test->obs, "z", 0, 0.5, NULL);

  nc_galaxy_wl_obs_set_shape_interp (test->obs, 0, spp_0);
  nc_galaxy_wl_obs_set_shape_interp (test->obs, 3, spp_3);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 0));
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp_source (test->obs));

  nc_galaxy_wl_obs_set_shape_interp_source (test->obs, source);
  g_assert_true (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 0) == spp_0);
  g_assert_true (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 3) == spp_3);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 1));
  g_assert_true (nc_galaxy_wl_obs_peek_shape_interp_source (test->obs) == source);

  /* The round trip keeps the interpolants and the record */
  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (test->obs));
  ncm_serialize_reset (ser, FALSE);
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));
  g_variant_unref (obs_ser);
  ncm_serialize_reset (ser, FALSE);

  g_assert_cmpfloat (nc_galaxy_wl_obs_get (obs_dup, "z", 0, NULL), ==, 0.5);
  g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 0)));
  g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 3)));
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 1));
  g_assert_true (NC_IS_GALAXY_SHAPE_POP_GAUSS (ncm_obj_dict_str_peek (nc_galaxy_wl_obs_peek_shape_interp_source (obs_dup), "shape-pop")));

  /* A change of the rows makes them stale, in both catalogs independently */
  nc_galaxy_wl_obs_set (obs_dup, "z", 1, 1.0, NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 0));
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp_source (obs_dup));
  g_assert_true (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 0) == spp_0);

  g_object_get (obs_dup, "shape-interp", &prop, NULL);
  g_assert_null (prop);

  /* A stale catalog serializes no interpolants */
  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (obs_dup));
  ncm_serialize_reset (ser, FALSE);
  nc_galaxy_wl_obs_free (obs_dup);
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));
  g_variant_unref (obs_ser);

  g_assert_cmpfloat (nc_galaxy_wl_obs_get (obs_dup, "z", 1, NULL), ==, 1.0);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 0));
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp_source (obs_dup));

  /* A new record makes the current rows the described ones */
  nc_galaxy_wl_obs_set (test->obs, "z", 2, 2.0, NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 0));
  nc_galaxy_wl_obs_set_shape_interp_source (test->obs, source);
  g_assert_true (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 0) == spp_0);

  nc_galaxy_wl_obs_free (obs_dup);
  nc_galaxy_shape_pop_gauss_free (pop);
  ncm_obj_dict_str_unref (source);
  ncm_spectral_tensor_free (spp_0);
  ncm_spectral_tensor_free (spp_3);
  ncm_serialize_free (ser);
}

/*
 * A catalog whose interpolants are backed by a file: they are read at the first use. Written
 * into a file of the same directory, the catalog carries the bare name of the interpolants
 * file and the record, not the interpolants, and reads back with the same interpolants; a
 * copy made without a file carries the absolute name; a change of the rows makes the
 * interpolants stale; setting an interpolant reads the file and keeps them in memory.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_file (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectralTensor *st = ncm_spectral_tensor_new (2);
  GError *error         = NULL;
  gchar *dir            = g_dir_make_tmp ("nc_galaxy_wl_obs_XXXXXX", &error);
  gchar *tables, *catalog, *data;
  NcGalaxyWLObsShapeInterp *si;
  NcGalaxyWLObs *obs_file, *obs_dup;
  GVariant *var, *props, *si_var, *si_props, *name, *val;
  gsize len;

  g_assert_no_error (error);
  tables  = _test_nc_galaxy_wl_obs_prepare_file (test->obs, dir);
  catalog = g_build_filename (dir, "catalog.gvar", NULL);

  g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 0)));
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (test->obs, 1));
  g_assert_true (NC_IS_GALAXY_SHAPE_POP_GAUSS (ncm_obj_dict_str_peek (nc_galaxy_wl_obs_peek_shape_interp_source (test->obs), "shape-pop")));

  /* Written next to the interpolants file: the bare name and the record, no interpolants */
  ncm_serialize_to_binfile (ser, G_OBJECT (test->obs), catalog);
  ncm_serialize_reset (ser, TRUE);
  g_assert_true (g_file_get_contents (catalog, &data, &len, &error));
  g_assert_no_error (error);
  var      = g_variant_ref_sink (g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), data, len, TRUE, g_free, data));
  props    = g_variant_get_child_value (var, 1);
  si_var   = g_variant_lookup_value (props, "shape-interp", G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE));
  si_props = g_variant_get_child_value (si_var, 1);
  name     = g_variant_lookup_value (si_props, "filename", G_VARIANT_TYPE_STRING);
  g_assert_cmpstr (g_variant_get_string (name, NULL), ==, "tables.gvar");
  val = g_variant_lookup_value (si_props, "interp", NULL);
  g_assert_null (val);
  val = g_variant_lookup_value (si_props, "source", NULL);
  g_assert_nonnull (val);
  g_variant_unref (val);

  /* Read back, the name resolves against the directory of the catalog file */
  obs_file = NC_GALAXY_WL_OBS (ncm_serialize_from_binfile (ser, catalog));
  ncm_serialize_reset (ser, TRUE);
  g_object_get (obs_file, "shape-interp", &si, NULL);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, tables);
  nc_galaxy_wl_obs_shape_interp_free (si);
  g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (obs_file, 3)));

  /* A copy made without a file carries the absolute name */
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_dup_obj (ser, G_OBJECT (obs_file)));
  ncm_serialize_reset (ser, TRUE);
  g_object_get (obs_dup, "shape-interp", &si, NULL);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, tables);
  nc_galaxy_wl_obs_shape_interp_free (si);
  g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 0)));

  /* A change of the rows makes the interpolants stale */
  nc_galaxy_wl_obs_set (obs_dup, "z", 1, 1.0, NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs_dup, 0));
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp_source (obs_dup));
  g_object_get (obs_dup, "shape-interp", &si, NULL);
  g_assert_null (si);

  /* Setting an interpolant reads the file and keeps the interpolants in memory */
  nc_galaxy_wl_obs_set_shape_interp (obs_file, 5, st);
  g_object_get (obs_file, "shape-interp", &si, NULL);
  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_filename (si));
  g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0, NULL)));
  g_assert_true (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 5, NULL) == st);
  nc_galaxy_wl_obs_shape_interp_free (si);

  g_assert_cmpint (g_remove (tables), ==, 0);
  g_assert_cmpint (g_remove (catalog), ==, 0);
  g_assert_cmpint (g_rmdir (dir), ==, 0);

  g_variant_unref (name);
  g_variant_unref (si_props);
  g_variant_unref (si_var);
  g_variant_unref (props);
  g_variant_unref (var);
  nc_galaxy_wl_obs_free (obs_dup);
  nc_galaxy_wl_obs_free (obs_file);
  ncm_spectral_tensor_free (st);
  ncm_serialize_free (ser);
  g_free (catalog);
  g_free (tables);
  g_free (dir);
}

/*
 * nc_galaxy_wl_obs_new_from_file() reads a catalog, and on a file of another type of object
 * sets NC_GALAXY_WL_OBS_ERROR_WRONG_TYPE.
 */
static void
test_nc_galaxy_wl_obs_new_from_file (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmVector *v      = ncm_vector_new (1);
  GError *error     = NULL;
  gchar *dir        = g_dir_make_tmp ("nc_galaxy_wl_obs_XXXXXX", &error);
  gchar *filename;
  NcGalaxyWLObs *obs;

  g_assert_no_error (error);
  filename = g_build_filename (dir, "catalog.gvar", NULL);

  nc_galaxy_wl_obs_set (test->obs, "z", 0, 0.5, NULL);
  ncm_serialize_to_binfile (ser, G_OBJECT (test->obs), filename);
  ncm_serialize_reset (ser, FALSE);

  obs = nc_galaxy_wl_obs_new_from_file (filename, NULL);
  g_assert_cmpfloat (nc_galaxy_wl_obs_get (obs, "z", 0, NULL), ==, 0.5);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs, 0));
  nc_galaxy_wl_obs_free (obs);

  ncm_vector_set (v, 0, 1.0);
  ncm_serialize_to_binfile (ser, G_OBJECT (v), filename);
  g_assert_null (nc_galaxy_wl_obs_new_from_file (filename, &error));
  g_assert_error (error, NC_GALAXY_WL_OBS_ERROR, NC_GALAXY_WL_OBS_ERROR_WRONG_TYPE);
  g_assert_nonnull (strstr (error->message, "holds a NcmVector, not a NcGalaxyWLObs"));
  g_clear_error (&error);

  g_assert_cmpint (g_remove (filename), ==, 0);
  g_assert_cmpint (g_rmdir (dir), ==, 0);

  ncm_vector_free (v);
  ncm_serialize_free (ser);
  g_free (filename);
  g_free (dir);
}

/*
 * A catalog written before NcGalaxyPz holds a NcmSpline in its pz property: building
 * it aborts with a message that names the update tool. This test is dropped in 1.0, with
 * the check it covers.
 */
static void
test_nc_galaxy_wl_obs_pz_spline (void)
{
  g_test_trap_subprocess ("/numcosmo/nc_galaxy_wl_obs/pz_spline/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the catalog holds the p(z) of row 3 as a NcmSplineCubicNotaknot*");
  g_test_trap_assert_stderr ("*tools/nc_galaxy_pz_update_catalog.py*");
}

static void
test_nc_galaxy_wl_obs_pz_spline_subprocess (void)
{
  GStrv col_names   = g_strsplit ("ra dec", " ", -1);
  NcmObjDictInt *pz = ncm_obj_dict_int_new ();
  NcmVector *xv     = ncm_vector_new (8);
  NcmVector *yv     = ncm_vector_new (8);
  NcmSpline *spline;
  guint j;

  for (j = 0; j < 8; j++)
  {
    ncm_vector_set (xv, j, 0.1 * j);
    ncm_vector_set (yv, j, 1.0 + j);
  }

  spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));
  ncm_obj_dict_int_add (pz, 3, G_OBJECT (spline));

  g_object_new (NC_TYPE_GALAXY_WL_OBS,
                "ellip-conv", NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                "columns", col_names,
                "coord", NC_GALAXY_WL_OBS_COORD_CELESTIAL,
                "len", 10,
                "pz", pz,
                NULL);
}

/*
 * An interpolant property holding an object that is not a NcmSpectralTensor: reading the
 * interpolant aborts with a message that names the row and the type.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_type (void)
{
  g_test_trap_subprocess ("/numcosmo/nc_galaxy_wl_obs/shape_interp/type/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the interpolant of row 2 is a NcmVector, not a NcmSpectralTensor*");
}

static void
test_nc_galaxy_wl_obs_shape_interp_type_subprocess (void)
{
  GStrv col_names       = g_strsplit ("ra dec", " ", -1);
  NcmObjDictInt *interp = ncm_obj_dict_int_new ();
  NcmObjDictStr *source = ncm_obj_dict_str_new ();
  NcmVector *v          = ncm_vector_new (1);
  NcGalaxyWLObsShapeInterp *si;
  NcGalaxyWLObs *obs;

  ncm_obj_dict_int_add (interp, 2, G_OBJECT (v));
  si = g_object_new (NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP, "interp", interp, NULL);

  obs = g_object_new (NC_TYPE_GALAXY_WL_OBS,
                      "ellip-conv", NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                      "columns", col_names,
                      "coord", NC_GALAXY_WL_OBS_COORD_CELESTIAL,
                      "len", 10,
                      NULL);

  g_object_set (obs, "shape-interp", si, NULL);
  nc_galaxy_wl_obs_set_shape_interp_source (obs, source);
  nc_galaxy_wl_obs_peek_shape_interp (obs, 2);
}

/*
 * A catalog read from a file whose rows change before the first use of its interpolants: the
 * comparison of the data hash aborts with a message that gives both hashes.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_hash (void)
{
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  GStrv col_names    = g_strsplit ("ra dec z", " ", -1);
  NcGalaxyWLObs *obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_GALAXY_WL_OBS_COORD_CELESTIAL, 10, col_names);
  GError *error      = NULL;
  gchar *dir         = g_dir_make_tmp ("nc_galaxy_wl_obs_XXXXXX", &error);
  gchar *tables, *catalog;

  g_assert_no_error (error);
  tables  = _test_nc_galaxy_wl_obs_prepare_file (obs, dir);
  catalog = g_build_filename (dir, "catalog.gvar", NULL);
  ncm_serialize_to_binfile (ser, G_OBJECT (obs), catalog);

  g_setenv (TEST_NC_GALAXY_WL_OBS_FILE_ENV, catalog, TRUE);
  g_test_trap_subprocess ("/numcosmo/nc_galaxy_wl_obs/shape_interp/hash/subprocess", 0, 0);
  g_unsetenv (TEST_NC_GALAXY_WL_OBS_FILE_ENV);

  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the interpolants describe the rows of data hash *, not these rows of data hash *");

  g_assert_cmpint (g_remove (tables), ==, 0);
  g_assert_cmpint (g_remove (catalog), ==, 0);
  g_assert_cmpint (g_rmdir (dir), ==, 0);

  nc_galaxy_wl_obs_free (obs);
  ncm_serialize_free (ser);
  g_strfreev (col_names);
  g_free (catalog);
  g_free (tables);
  g_free (dir);
}

static void
test_nc_galaxy_wl_obs_shape_interp_hash_subprocess (void)
{
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObs *obs = NC_GALAXY_WL_OBS (ncm_serialize_from_binfile (ser, g_getenv (TEST_NC_GALAXY_WL_OBS_FILE_ENV)));

  nc_galaxy_wl_obs_set (obs, "z", 1, 1.0, NULL);
  nc_galaxy_wl_obs_peek_shape_interp (obs, 0);
}

/* A grid of @obs as its rows are now: two weights and two redshift values per galaxy. */
static NcGalaxyWLObsChebGrid *
_test_nc_galaxy_wl_obs_cheb_grid_new (NcGalaxyWLObs *obs)
{
  gchar *hash                 = ncm_catalog_get_data_hash (NCM_CATALOG (obs));
  NcGalaxyWLObsChebGrid *grid = nc_galaxy_wl_obs_cheb_grid_new (nc_galaxy_wl_obs_len (obs), 2, 2, 0.3, "NcGalaxyRedshiftFactorSpline", hash);

  g_free (hash);

  return grid;
}

/*
 * The grid is kept while the rows and the p(z) are unchanged, travels with a serialized
 * catalog, and is dropped by a change of a row or of a p(z).
 */
static void
test_nc_galaxy_wl_obs_cheb_grid (TestNcGalaxyWLObs *test, gconstpointer pdata)
{
  NcmSerialize *ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObsChebGrid *grid = _test_nc_galaxy_wl_obs_cheb_grid_new (test->obs);
  NcGalaxyWLObs *dup;
  GVariant *var;

  g_assert_null (nc_galaxy_wl_obs_peek_cheb_grid (test->obs));

  nc_galaxy_wl_obs_set_cheb_grid (test->obs, grid);
  g_assert_true (nc_galaxy_wl_obs_peek_cheb_grid (test->obs) == grid);

  var = ncm_serialize_to_variant (ser, G_OBJECT (test->obs));
  dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, var));
  g_assert_nonnull (nc_galaxy_wl_obs_peek_cheb_grid (dup));
  g_assert_cmpuint (nc_galaxy_wl_obs_cheb_grid_len (nc_galaxy_wl_obs_peek_cheb_grid (dup)), ==, nc_galaxy_wl_obs_len (test->obs));

  nc_galaxy_wl_obs_set (dup, "ra", 0, 1.0 + nc_galaxy_wl_obs_get (dup, "ra", 0, NULL), NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_cheb_grid (dup));

  {
    NcmVector *zv = ncm_vector_new (8);
    NcmVector *pv = ncm_vector_new (8);
    NcGalaxyPzBSplineSqrt *pz;
    guint k;

    for (k = 0; k < 8; k++)
    {
      ncm_vector_set (zv, k, 0.1 * k);
      ncm_vector_set (pv, k, 1.0 + k);
    }

    pz = nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv);
    nc_galaxy_wl_obs_set_pz (test->obs, 0, NC_GALAXY_PZ (pz));
    g_assert_null (nc_galaxy_wl_obs_peek_cheb_grid (test->obs));

    nc_galaxy_pz_free (NC_GALAXY_PZ (pz));
    ncm_vector_free (zv);
    ncm_vector_free (pv);
  }

  nc_galaxy_wl_obs_set_cheb_grid (test->obs, grid);
  g_assert_nonnull (nc_galaxy_wl_obs_peek_cheb_grid (test->obs));
  nc_galaxy_wl_obs_set_cheb_grid (test->obs, NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_cheb_grid (test->obs));

  g_variant_unref (var);
  ncm_serialize_free (ser);
  nc_galaxy_wl_obs_cheb_grid_free (grid);
  NCM_TEST_FREE (nc_galaxy_wl_obs_free, dup);
}

static void
test_nc_galaxy_wl_obs_cheb_grid_errors (void)
{
  const gchar *cases[][2] = {
    {"len",  "*nc_galaxy_wl_obs_set_cheb_grid: the grid has 4 galaxies, the catalog 5*"},
    {"hash", "*the Chebyshev grid describes the rows of data hash other, not these rows of data hash *"},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    g_setenv (TEST_NC_GALAXY_WL_OBS_FILE_ENV, cases[i][0], TRUE);
    g_test_trap_subprocess ("/numcosmo/nc_galaxy_wl_obs/cheb_grid/errors/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (cases[i][1]);
  }

  g_unsetenv (TEST_NC_GALAXY_WL_OBS_FILE_ENV);
}

static void
test_nc_galaxy_wl_obs_cheb_grid_errors_subprocess (void)
{
  const gchar *which = g_getenv (TEST_NC_GALAXY_WL_OBS_FILE_ENV);
  GStrv col_names    = g_strsplit ("ra dec z", " ", -1);
  NcGalaxyWLObs *obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_GALAXY_WL_OBS_COORD_CELESTIAL, 5, col_names);

  if (g_strcmp0 (which, "len") == 0)
  {
    nc_galaxy_wl_obs_set_cheb_grid (obs, nc_galaxy_wl_obs_cheb_grid_new (4, 2, 2, 0.3, "NcGalaxyRedshiftFactorSpline", "other"));
  }
  else
  {
    g_object_set (obs, "cheb-grid", nc_galaxy_wl_obs_cheb_grid_new (5, 2, 2, 0.3, "NcGalaxyRedshiftFactorSpline", "other"), NULL);
    nc_galaxy_wl_obs_peek_cheb_grid (obs);
  }
}


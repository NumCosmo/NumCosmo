/***************************************************************************
 *            test_nc_galaxy_wl_obs_shape_interp.c
 *
 *  Wed November 25 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_wl_obs_shape_interp.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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

#include <glib.h>
#include <glib-object.h>
#include <glib/gstdio.h>

/* The environment variables that give the abort tests their files */
#define TEST_SHAPE_INTERP_FILE_ENV "NC_GALAXY_WL_OBS_SHAPE_INTERP_TEST_FILE"
#define TEST_SHAPE_INTERP_OUT_ENV "NC_GALAXY_WL_OBS_SHAPE_INTERP_TEST_OUT"

/* The data hash of the tests: any string, compared as written */
#define TEST_SHAPE_INTERP_HASH "0123abcd"

typedef struct _TestNcGalaxyWLObsShapeInterp
{
  NcGalaxyWLObsShapeInterp *si;
  NcmObjDictStr *source;
  gchar *dir;
  gchar *filename;
} TestNcGalaxyWLObsShapeInterp;

static void _test_nc_galaxy_wl_obs_shape_interp_trap (const gchar *filename, const gchar *out, const gchar *subprocess, const gchar *pattern);

static void test_nc_galaxy_wl_obs_shape_interp_new (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_free (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);

static void test_nc_galaxy_wl_obs_shape_interp_memory (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_file (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_lazy (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_write (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_traps (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_write_traps (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata);
static void test_nc_galaxy_wl_obs_shape_interp_other_dir_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_hash_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_type_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_refers_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_both_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_no_hash_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_writing_peek_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_writing_serialize_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_close_no_hash_subprocess (void);
static void test_nc_galaxy_wl_obs_shape_interp_close_not_writing_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_wl_obs_shape_interp/memory", TestNcGalaxyWLObsShapeInterp, NULL,
              &test_nc_galaxy_wl_obs_shape_interp_new,
              &test_nc_galaxy_wl_obs_shape_interp_memory,
              &test_nc_galaxy_wl_obs_shape_interp_free);
  g_test_add ("/nc/galaxy_wl_obs_shape_interp/file", TestNcGalaxyWLObsShapeInterp, NULL,
              &test_nc_galaxy_wl_obs_shape_interp_new,
              &test_nc_galaxy_wl_obs_shape_interp_file,
              &test_nc_galaxy_wl_obs_shape_interp_free);
  g_test_add ("/nc/galaxy_wl_obs_shape_interp/lazy", TestNcGalaxyWLObsShapeInterp, NULL,
              &test_nc_galaxy_wl_obs_shape_interp_new,
              &test_nc_galaxy_wl_obs_shape_interp_lazy,
              &test_nc_galaxy_wl_obs_shape_interp_free);
  g_test_add ("/nc/galaxy_wl_obs_shape_interp/write", TestNcGalaxyWLObsShapeInterp, NULL,
              &test_nc_galaxy_wl_obs_shape_interp_new,
              &test_nc_galaxy_wl_obs_shape_interp_write,
              &test_nc_galaxy_wl_obs_shape_interp_free);
  g_test_add ("/nc/galaxy_wl_obs_shape_interp/traps", TestNcGalaxyWLObsShapeInterp, NULL,
              &test_nc_galaxy_wl_obs_shape_interp_new,
              &test_nc_galaxy_wl_obs_shape_interp_traps,
              &test_nc_galaxy_wl_obs_shape_interp_free);
  g_test_add ("/nc/galaxy_wl_obs_shape_interp/write_traps", TestNcGalaxyWLObsShapeInterp, NULL,
              &test_nc_galaxy_wl_obs_shape_interp_new,
              &test_nc_galaxy_wl_obs_shape_interp_write_traps,
              &test_nc_galaxy_wl_obs_shape_interp_free);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/other_dir/subprocess", &test_nc_galaxy_wl_obs_shape_interp_other_dir_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/hash/subprocess", &test_nc_galaxy_wl_obs_shape_interp_hash_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/type/subprocess", &test_nc_galaxy_wl_obs_shape_interp_type_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/refers/subprocess", &test_nc_galaxy_wl_obs_shape_interp_refers_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/both/subprocess", &test_nc_galaxy_wl_obs_shape_interp_both_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/no_hash/subprocess", &test_nc_galaxy_wl_obs_shape_interp_no_hash_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/writing_peek/subprocess", &test_nc_galaxy_wl_obs_shape_interp_writing_peek_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/writing_serialize/subprocess", &test_nc_galaxy_wl_obs_shape_interp_writing_serialize_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/close_no_hash/subprocess", &test_nc_galaxy_wl_obs_shape_interp_close_no_hash_subprocess);
  g_test_add_func ("/nc/galaxy_wl_obs_shape_interp/close_not_writing/subprocess", &test_nc_galaxy_wl_obs_shape_interp_close_not_writing_subprocess);

  g_test_run ();

  return 0;
}

/*
 * Runs the test subprocess with filename and out in the environment, and asserts that it
 * aborts with a message that matches pattern.
 */
static void
_test_nc_galaxy_wl_obs_shape_interp_trap (const gchar *filename, const gchar *out, const gchar *subprocess, const gchar *pattern)
{
  g_setenv (TEST_SHAPE_INTERP_FILE_ENV, filename, TRUE);
  g_setenv (TEST_SHAPE_INTERP_OUT_ENV, (out != NULL) ? out : "", TRUE);
  g_test_trap_subprocess (subprocess, 0, 0);
  g_unsetenv (TEST_SHAPE_INTERP_FILE_ENV);
  g_unsetenv (TEST_SHAPE_INTERP_OUT_ENV);

  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr (pattern);
}

/*
 * An object in memory with interpolants of rows 0 and 3, a record and the data hash
 * TEST_SHAPE_INTERP_HASH, written to test->filename in a new directory.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_new (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyShapePopGauss *pop = nc_galaxy_shape_pop_gauss_new ();
  NcmSpectralTensor *spp_0   = ncm_spectral_tensor_new (2);
  NcmSpectralTensor *spp_3   = ncm_spectral_tensor_new (3);
  GError *error              = NULL;

  test->dir = g_dir_make_tmp ("nc_wl_obs_shape_interp_XXXXXX", &error);
  g_assert_no_error (error);
  test->filename = g_build_filename (test->dir, "tables.gvar", NULL);

  test->source = ncm_obj_dict_str_new ();
  ncm_obj_dict_str_add (test->source, "shape-pop", G_OBJECT (pop));

  test->si = nc_galaxy_wl_obs_shape_interp_new ();
  nc_galaxy_wl_obs_shape_interp_set_interp (test->si, 0, spp_0);
  nc_galaxy_wl_obs_shape_interp_set_interp (test->si, 3, spp_3);
  nc_galaxy_wl_obs_shape_interp_set_source (test->si, test->source);
  nc_galaxy_wl_obs_shape_interp_set_data_hash (test->si, TEST_SHAPE_INTERP_HASH);

  ncm_serialize_to_binfile (ser, G_OBJECT (test->si), test->filename);

  nc_galaxy_shape_pop_gauss_free (pop);
  ncm_spectral_tensor_free (spp_0);
  ncm_spectral_tensor_free (spp_3);
  ncm_serialize_free (ser);
}

static void
test_nc_galaxy_wl_obs_shape_interp_free (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  g_assert_cmpint (g_remove (test->filename), ==, 0);
  g_assert_cmpint (g_rmdir (test->dir), ==, 0);

  NCM_TEST_FREE (nc_galaxy_wl_obs_shape_interp_free, test->si);
  ncm_obj_dict_str_unref (test->source);
  g_free (test->filename);
  g_free (test->dir);
}

/*
 * In memory: the interpolants, the record and the data hash; no file; a copy made through
 * a variant holds the same; nc_galaxy_wl_obs_shape_interp_ref() and
 * nc_galaxy_wl_obs_shape_interp_clear() manage the reference count.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_memory (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObsShapeInterp *dup, *ref;

  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_filename (test->si));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (test->si, 3)), ==, 3);
  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_interp (test->si, 1));
  g_assert_true (nc_galaxy_wl_obs_shape_interp_peek_source (test->si) == test->source);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_data_hash (test->si), ==, TEST_SHAPE_INTERP_HASH);

  dup = NC_GALAXY_WL_OBS_SHAPE_INTERP (ncm_serialize_dup_obj (ser, G_OBJECT (test->si)));
  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_filename (dup));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (dup, 3)), ==, 3);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_data_hash (dup), ==, TEST_SHAPE_INTERP_HASH);
  g_assert_true (NC_IS_GALAXY_SHAPE_POP_GAUSS (ncm_obj_dict_str_peek (nc_galaxy_wl_obs_shape_interp_peek_source (dup), "shape-pop")));

  ref = nc_galaxy_wl_obs_shape_interp_ref (dup);
  nc_galaxy_wl_obs_shape_interp_free (dup);
  nc_galaxy_wl_obs_shape_interp_clear (&ref);
  g_assert_null (ref);

  ncm_serialize_free (ser);
}

/*
 * Backed by the file of the fixture: the name is absolute; the interpolants are read at the
 * first use; written into a file of the same directory the object carries the bare name and
 * reads back from there; written without a file it carries the absolute name; a relative
 * name given outside a file function resolves against the base directory; a change of the
 * data hash reads the file and keeps the interpolants in memory.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_file (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser             = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObsShapeInterp *si  = nc_galaxy_wl_obs_shape_interp_new_from_file (test->filename, TEST_SHAPE_INTERP_HASH, test->source);
  gchar *ref_name               = g_build_filename (test->dir, "ref.gvar", NULL);
  gchar *base_name              = g_build_filename (ncm_serialize_context_peek_dir (ncm_serialize_context_peek_current ()), "tables.gvar", NULL);
  NcGalaxyWLObsShapeInterp *dup = NULL;
  NcGalaxyWLObsShapeInterp *rel;
  GVariant *var, *props, *name;
  gchar *written;

  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, test->filename);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0)), ==, 2);
  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 1));

  /* Same directory: the bare name, and the copy reads from there */
  ncm_serialize_to_binfile (ser, G_OBJECT (si), ref_name);
  ncm_serialize_reset (ser, TRUE);
  dup = NC_GALAXY_WL_OBS_SHAPE_INTERP (ncm_serialize_from_binfile (ser, ref_name));
  ncm_serialize_reset (ser, TRUE);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (dup), ==, test->filename);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (dup, 3)), ==, 3);
  nc_galaxy_wl_obs_shape_interp_clear (&dup);

  {
    gchar *data;
    gsize len;

    g_assert_true (g_file_get_contents (ref_name, &data, &len, NULL));
    var   = g_variant_ref_sink (g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), data, len, TRUE, g_free, data));
    props = g_variant_get_child_value (var, 1);
    name  = g_variant_lookup_value (props, "filename", G_VARIANT_TYPE_STRING);
    g_assert_cmpstr (g_variant_get_string (name, NULL), ==, "tables.gvar");
    g_assert_null (g_variant_lookup_value (props, "interp", NULL));
    g_variant_unref (name);
    g_variant_unref (props);
    g_variant_unref (var);
  }

  /* Without a file: the absolute name */
  g_object_get (si, "filename", &written, NULL);
  g_assert_cmpstr (written, ==, test->filename);
  g_free (written);

  /* A relative name outside a file function: the base directory */
  rel = nc_galaxy_wl_obs_shape_interp_new_from_file ("tables.gvar", TEST_SHAPE_INTERP_HASH, test->source);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (rel), ==, base_name);
  nc_galaxy_wl_obs_shape_interp_free (rel);

  /* A change reads the file and keeps the interpolants in memory */
  nc_galaxy_wl_obs_shape_interp_set_data_hash (si, "ffff");
  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_filename (si));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 3)), ==, 3);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_data_hash (si), ==, "ffff");

  g_assert_cmpint (g_remove (ref_name), ==, 0);

  nc_galaxy_wl_obs_shape_interp_free (si);
  ncm_serialize_free (ser);
  g_free (base_name);
  g_free (ref_name);
}

/*
 * An object backed by a file that does not exist: its record, data hash and name are read
 * and it is serialized without reading the file.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_lazy (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *missing               = g_build_filename (test->dir, "missing.gvar", NULL);
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new_from_file (missing, TEST_SHAPE_INTERP_HASH, test->source);
  NcGalaxyWLObsShapeInterp *dup;

  g_assert_true (nc_galaxy_wl_obs_shape_interp_peek_source (si) == test->source);
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_data_hash (si), ==, TEST_SHAPE_INTERP_HASH);

  dup = NC_GALAXY_WL_OBS_SHAPE_INTERP (ncm_serialize_dup_obj (ser, G_OBJECT (si)));
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (dup), ==, missing);

  nc_galaxy_wl_obs_shape_interp_free (dup);
  nc_galaxy_wl_obs_shape_interp_free (si);
  ncm_serialize_free (ser);
  g_free (missing);
}

/*
 * Written through create_empty: the records file exists until close; the last interpolant
 * of a row set twice is written; after close the records file is gone and the object is
 * backed by the file, which another object reads with the same interpolants, record and
 * data hash. An object freed before close removes its records file and writes nothing.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_write (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  gchar *written_name          = g_build_filename (test->dir, "written.gvar", NULL);
  gchar *records_name          = g_strconcat (written_name, ".records", NULL);
  gchar *unclosed_name         = g_build_filename (test->dir, "unclosed.gvar", NULL);
  gchar *unclosed_records      = g_strconcat (unclosed_name, ".records", NULL);
  NcmSpectralTensor *spp_1     = ncm_spectral_tensor_new (1);
  NcmSpectralTensor *spp_2     = ncm_spectral_tensor_new (2);
  NcmSpectralTensor *spp_3     = ncm_spectral_tensor_new (3);
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_create_empty (written_name);
  NcGalaxyWLObsShapeInterp *reader;

  g_assert_true (g_file_test (records_name, G_FILE_TEST_EXISTS));
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, written_name);

  nc_galaxy_wl_obs_shape_interp_set_interp (si, 3, spp_1);
  nc_galaxy_wl_obs_shape_interp_set_interp (si, 0, spp_2);
  nc_galaxy_wl_obs_shape_interp_set_interp (si, 3, spp_3);
  nc_galaxy_wl_obs_shape_interp_set_source (si, test->source);
  nc_galaxy_wl_obs_shape_interp_set_data_hash (si, TEST_SHAPE_INTERP_HASH);
  nc_galaxy_wl_obs_shape_interp_close (si);

  g_assert_false (g_file_test (records_name, G_FILE_TEST_EXISTS));
  g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_filename (si), ==, written_name);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 3)), ==, 3);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0)), ==, 2);
  g_assert_null (nc_galaxy_wl_obs_shape_interp_peek_interp (si, 1));

  reader = nc_galaxy_wl_obs_shape_interp_new_from_file (written_name, TEST_SHAPE_INTERP_HASH, NULL);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (reader, 3)), ==, 3);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_shape_interp_peek_interp (reader, 0)), ==, 2);
  nc_galaxy_wl_obs_shape_interp_free (reader);

  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    GObject *obj      = ncm_serialize_from_binfile (ser, written_name);

    g_assert_true (NC_IS_GALAXY_WL_OBS_SHAPE_INTERP (obj));
    g_assert_cmpstr (nc_galaxy_wl_obs_shape_interp_peek_data_hash (NC_GALAXY_WL_OBS_SHAPE_INTERP (obj)), ==, TEST_SHAPE_INTERP_HASH);
    g_assert_true (NC_IS_GALAXY_SHAPE_POP_GAUSS (ncm_obj_dict_str_peek (nc_galaxy_wl_obs_shape_interp_peek_source (NC_GALAXY_WL_OBS_SHAPE_INTERP (obj)), "shape-pop")));
    g_object_unref (obj);
    ncm_serialize_free (ser);
  }

  /* Freed before close: nothing written */
  {
    NcGalaxyWLObsShapeInterp *unclosed = nc_galaxy_wl_obs_shape_interp_create_empty (unclosed_name);

    nc_galaxy_wl_obs_shape_interp_set_interp (unclosed, 0, spp_1);
    g_assert_true (g_file_test (unclosed_records, G_FILE_TEST_EXISTS));
    nc_galaxy_wl_obs_shape_interp_free (unclosed);
    g_assert_false (g_file_test (unclosed_records, G_FILE_TEST_EXISTS));
    g_assert_false (g_file_test (unclosed_name, G_FILE_TEST_EXISTS));
  }

  g_assert_cmpint (g_remove (written_name), ==, 0);

  nc_galaxy_wl_obs_shape_interp_free (si);
  ncm_spectral_tensor_free (spp_1);
  ncm_spectral_tensor_free (spp_2);
  ncm_spectral_tensor_free (spp_3);
  g_free (unclosed_records);
  g_free (unclosed_name);
  g_free (records_name);
  g_free (written_name);
}

/*
 * Each of these aborts with a message that names the cause: writing an object backed by a
 * file into a file of another directory; reading a file of another data hash, of another
 * type of object, or that refers to another file; constructing an object with both a file
 * and interpolants, or with a file and no data hash.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_traps (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  GError *error     = NULL;
  gchar *other_dir  = g_dir_make_tmp ("nc_wl_obs_shape_interp_XXXXXX", &error);
  gchar *out        = g_build_filename (other_dir, "ref.gvar", NULL);
  gchar *vec_name   = g_build_filename (test->dir, "vector.gvar", NULL);
  gchar *ref_name   = g_build_filename (test->dir, "ref.gvar", NULL);
  NcmVector *v      = ncm_vector_new (1);

  g_assert_no_error (error);

  _test_nc_galaxy_wl_obs_shape_interp_trap (test->filename, out, "/nc/galaxy_wl_obs_shape_interp/other_dir/subprocess",
                                            "*the interpolants are in *tables.gvar, which is not in *");
  _test_nc_galaxy_wl_obs_shape_interp_trap (test->filename, NULL, "/nc/galaxy_wl_obs_shape_interp/hash/subprocess",
                                            "*holds the interpolants of the data hash " TEST_SHAPE_INTERP_HASH ", not ffff*");

  ncm_vector_set (v, 0, 1.0);
  ncm_serialize_to_binfile (ser, G_OBJECT (v), vec_name);
  _test_nc_galaxy_wl_obs_shape_interp_trap (vec_name, NULL, "/nc/galaxy_wl_obs_shape_interp/type/subprocess",
                                            "*vector.gvar holds a NcmVector, not a NcGalaxyWLObsShapeInterp*");

  {
    NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new_from_file (test->filename, TEST_SHAPE_INTERP_HASH, test->source);

    ncm_serialize_to_binfile (ser, G_OBJECT (si), ref_name);
    nc_galaxy_wl_obs_shape_interp_free (si);
  }

  _test_nc_galaxy_wl_obs_shape_interp_trap (ref_name, NULL, "/nc/galaxy_wl_obs_shape_interp/refers/subprocess",
                                            "*ref.gvar refers to the file *tables.gvar instead of holding the interpolants*");
  _test_nc_galaxy_wl_obs_shape_interp_trap (test->filename, NULL, "/nc/galaxy_wl_obs_shape_interp/both/subprocess",
                                            "*an object backed by the file *tables.gvar also holds interpolants*");
  _test_nc_galaxy_wl_obs_shape_interp_trap (test->filename, NULL, "/nc/galaxy_wl_obs_shape_interp/no_hash/subprocess",
                                            "*an object backed by the file *tables.gvar has no data hash*");

  g_assert_cmpint (g_remove (vec_name), ==, 0);
  g_assert_cmpint (g_remove (ref_name), ==, 0);
  g_assert_cmpint (g_rmdir (other_dir), ==, 0);

  ncm_vector_free (v);
  ncm_serialize_free (ser);
  g_free (ref_name);
  g_free (vec_name);
  g_free (out);
  g_free (other_dir);
}

/*
 * Each of these aborts with a message that names the cause: reading an interpolant or
 * serializing an object being written; closing one without a data hash; closing an object
 * that is not being written. The parent removes the records file each subprocess leaves.
 */
static void
test_nc_galaxy_wl_obs_shape_interp_write_traps (TestNcGalaxyWLObsShapeInterp *test, gconstpointer pdata)
{
  gchar *writing_name = g_build_filename (test->dir, "writing.gvar", NULL);
  gchar *records_name = g_strconcat (writing_name, ".records", NULL);

  _test_nc_galaxy_wl_obs_shape_interp_trap (writing_name, NULL, "/nc/galaxy_wl_obs_shape_interp/writing_peek/subprocess",
                                            "*the file *writing.gvar is being written*");
  g_assert_cmpint (g_remove (records_name), ==, 0);
  _test_nc_galaxy_wl_obs_shape_interp_trap (writing_name, NULL, "/nc/galaxy_wl_obs_shape_interp/writing_serialize/subprocess",
                                            "*the file *writing.gvar is being written*");
  g_assert_cmpint (g_remove (records_name), ==, 0);
  _test_nc_galaxy_wl_obs_shape_interp_trap (writing_name, NULL, "/nc/galaxy_wl_obs_shape_interp/close_no_hash/subprocess",
                                            "*the object writing *writing.gvar has no data hash*");
  g_assert_cmpint (g_remove (records_name), ==, 0);
  _test_nc_galaxy_wl_obs_shape_interp_trap (writing_name, NULL, "/nc/galaxy_wl_obs_shape_interp/close_not_writing/subprocess",
                                            "*the object is not being written*");

  g_free (records_name);
  g_free (writing_name);
}

static void
test_nc_galaxy_wl_obs_shape_interp_other_dir_subprocess (void)
{
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new_from_file (g_getenv (TEST_SHAPE_INTERP_FILE_ENV), TEST_SHAPE_INTERP_HASH, NULL);

  ncm_serialize_to_binfile (ser, G_OBJECT (si), g_getenv (TEST_SHAPE_INTERP_OUT_ENV));
}

static void
test_nc_galaxy_wl_obs_shape_interp_hash_subprocess (void)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new_from_file (g_getenv (TEST_SHAPE_INTERP_FILE_ENV), "ffff", NULL);

  nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0);
}

static void
test_nc_galaxy_wl_obs_shape_interp_type_subprocess (void)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new_from_file (g_getenv (TEST_SHAPE_INTERP_FILE_ENV), TEST_SHAPE_INTERP_HASH, NULL);

  nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0);
}

static void
test_nc_galaxy_wl_obs_shape_interp_refers_subprocess (void)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new_from_file (g_getenv (TEST_SHAPE_INTERP_FILE_ENV), TEST_SHAPE_INTERP_HASH, NULL);

  nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0);
}

static void
test_nc_galaxy_wl_obs_shape_interp_both_subprocess (void)
{
  NcmObjDictInt *interp = ncm_obj_dict_int_new ();

  g_object_new (NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP,
                "filename", g_getenv (TEST_SHAPE_INTERP_FILE_ENV),
                "data-hash", TEST_SHAPE_INTERP_HASH,
                "interp", interp,
                NULL);
}

static void
test_nc_galaxy_wl_obs_shape_interp_no_hash_subprocess (void)
{
  g_object_new (NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP, "filename", g_getenv (TEST_SHAPE_INTERP_FILE_ENV), NULL);
}

static void
test_nc_galaxy_wl_obs_shape_interp_writing_peek_subprocess (void)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_create_empty (g_getenv (TEST_SHAPE_INTERP_FILE_ENV));

  nc_galaxy_wl_obs_shape_interp_peek_interp (si, 0);
}

static void
test_nc_galaxy_wl_obs_shape_interp_writing_serialize_subprocess (void)
{
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_create_empty (g_getenv (TEST_SHAPE_INTERP_FILE_ENV));

  nc_galaxy_wl_obs_shape_interp_set_data_hash (si, TEST_SHAPE_INTERP_HASH);
  ncm_serialize_to_variant (ser, G_OBJECT (si));
}

static void
test_nc_galaxy_wl_obs_shape_interp_close_no_hash_subprocess (void)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_create_empty (g_getenv (TEST_SHAPE_INTERP_FILE_ENV));

  nc_galaxy_wl_obs_shape_interp_close (si);
}

static void
test_nc_galaxy_wl_obs_shape_interp_close_not_writing_subprocess (void)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new ();

  nc_galaxy_wl_obs_shape_interp_close (si);
}


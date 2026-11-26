/***************************************************************************
 *            test_ncm_catalog.c
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_catalog.c
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <string.h>

typedef struct _TestNcmCatalog
{
  NcmCatalog *catalog;
  NcmCatalog *typed;
} TestNcmCatalog;

/* The largest object_id of the HSC S16A shape catalog, above 2^53. */
#define TEST_NCM_CATALOG_BIG_ID (G_GINT64_CONSTANT (38558648720430529))

static void _test_ncm_catalog_assert_column_not_found_error (GError **error, const gchar *col);
static void _test_ncm_catalog_assert_wrong_type_error (GError **error, const gchar *message);

void test_ncm_catalog_new (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_free (TestNcmCatalog *test, gconstpointer pdata);

void test_ncm_catalog_ref (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_index (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_set_get (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_types (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_columns (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_wrong_type (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_meta (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_serialize (TestNcmCatalog *test, gconstpointer pdata);

void test_ncm_catalog_invalid_get (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_invalid_set (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_invalid_col_type (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_invalid_row (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_data_prop (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_data_hash (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_data_hash_bytes (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_key (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_key_saved (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_errors (void);
void test_ncm_catalog_errors_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/catalog/ref", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_ref,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/index", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_index,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/set_get", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_set_get,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/types", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_types,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/columns", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_columns,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/wrong_type", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_wrong_type,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/meta", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_meta,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/serialize", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_serialize,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/invalid/get", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_invalid_get,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/invalid/set", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_invalid_set,
              &test_ncm_catalog_free);

  g_test_add ("/ncm/catalog/invalid/row", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_invalid_row,
              &test_ncm_catalog_free);
  g_test_add ("/ncm/catalog/data_prop", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_data_prop,
              &test_ncm_catalog_free);
  g_test_add ("/ncm/catalog/key", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_key,
              &test_ncm_catalog_free);
  g_test_add ("/ncm/catalog/key/saved", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_key_saved,
              &test_ncm_catalog_free);
  g_test_add ("/ncm/catalog/data_hash", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_data_hash,
              &test_ncm_catalog_free);
  g_test_add ("/ncm/catalog/data_hash/bytes", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_data_hash_bytes,
              &test_ncm_catalog_free);
  g_test_add_func ("/ncm/catalog/errors", &test_ncm_catalog_errors);
  g_test_add_func ("/ncm/catalog/errors/subprocess", &test_ncm_catalog_errors_subprocess);
  g_test_add ("/ncm/catalog/invalid/col_type", TestNcmCatalog, NULL,
              &test_ncm_catalog_new,
              &test_ncm_catalog_invalid_col_type,
              &test_ncm_catalog_free);

  g_test_run ();

  return 0;
}

static void
_test_ncm_catalog_assert_column_not_found_error (GError **error, const gchar *col)
{
  gchar *expected_message = g_strdup_printf ("Column '%s' not found.", col);

  g_assert_nonnull (*error);
  g_assert_error (*error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND);
  g_assert_cmpstr ((*error)->message, ==, expected_message);

  g_clear_error (error);
  g_assert_null (*error);
  g_free (expected_message);
}

static void
_test_ncm_catalog_assert_wrong_type_error (GError **error, const gchar *message)
{
  g_assert_error (*error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_WRONG_TYPE);
  g_assert_cmpstr ((*error)->message, ==, message);
  g_clear_error (error);
}

static void
_test_ncm_catalog_assert_row_error (GError **error, const gchar *col)
{
  gchar *expected_message = g_strdup_printf ("Row 4 of column '%s' requested, but the catalog has 4 rows.", col);

  g_assert_error (*error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE);
  g_assert_cmpstr ((*error)->message, ==, expected_message);
  g_clear_error (error);
  g_free (expected_message);
}

void
test_ncm_catalog_new (TestNcmCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[] = {"ra", "dec", "z", NULL};

  const gchar *typed_names[]            = {"ra", "id", "dec", "flag", NULL};
  const NcmCatalogColType typed_types[] = {
    NCM_CATALOG_COL_TYPE_DOUBLE,
    NCM_CATALOG_COL_TYPE_INT,
    NCM_CATALOG_COL_TYPE_DOUBLE,
    NCM_CATALOG_COL_TYPE_BOOL,
  };

  test->catalog = ncm_catalog_new (4, (GStrv) col_names);
  test->typed   = ncm_catalog_new_full (4, (GStrv) typed_names, typed_types, 4);

  g_assert_true (NCM_IS_CATALOG (test->catalog));
  g_assert_cmpuint (ncm_catalog_len (test->catalog), ==, 4);
  g_assert_cmpuint (ncm_catalog_ncols (test->catalog), ==, 3);
  g_assert_cmpuint (ncm_catalog_len (test->typed), ==, 4);
  g_assert_cmpuint (ncm_catalog_ncols (test->typed), ==, 4);
}

void
test_ncm_catalog_free (TestNcmCatalog *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_catalog_free, test->catalog);
  NCM_TEST_FREE (ncm_catalog_free, test->typed);
}

void
test_ncm_catalog_ref (TestNcmCatalog *test, gconstpointer pdata)
{
  NcmCatalog *catalog_ref = ncm_catalog_ref (test->catalog);

  g_assert_true (catalog_ref == test->catalog);

  ncm_catalog_clear (&catalog_ref);
  g_assert_null (catalog_ref);

  /* The original reference is still alive. */
  g_assert_true (NCM_IS_CATALOG (test->catalog));
}

void
test_ncm_catalog_index (TestNcmCatalog *test, gconstpointer pdata)
{
  guint i = G_MAXUINT;

  g_assert_true (ncm_catalog_has_column (test->catalog, "dec"));
  g_assert_false (ncm_catalog_has_column (test->catalog, "missing"));

  g_assert_true (ncm_catalog_get_index (test->catalog, "z", &i));
  g_assert_cmpuint (i, ==, 2);

  g_assert_false (ncm_catalog_get_index (test->catalog, "missing", &i));

  /* The catalog index counts every column, the data index only the DOUBLE ones. */
  g_assert_true (ncm_catalog_get_index (test->typed, "dec", &i));
  g_assert_cmpuint (i, ==, 2);
  g_assert_true (ncm_catalog_get_data_index (test->typed, "dec", &i));
  g_assert_cmpuint (i, ==, 1);
  g_assert_false (ncm_catalog_get_data_index (test->typed, "id", &i));
  g_assert_false (ncm_catalog_get_data_index (test->typed, "flag", &i));
  g_assert_false (ncm_catalog_get_data_index (test->typed, "missing", &i));
}

void
test_ncm_catalog_set_get (TestNcmCatalog *test, gconstpointer pdata)
{
  NcmMatrix *data;
  guint i;
  guint j;

  /* A new catalog is zero-initialized. */
  for (i = 0; i < ncm_catalog_len (test->catalog); i++)
  {
    g_assert_cmpfloat (ncm_catalog_get (test->catalog, "ra", i, NULL), ==, 0.0);
    g_assert_cmpfloat (ncm_catalog_get (test->catalog, "dec", i, NULL), ==, 0.0);
    g_assert_cmpfloat (ncm_catalog_get (test->catalog, "z", i, NULL), ==, 0.0);
  }

  ncm_catalog_set (test->catalog, "ra", 0, 1.5, NULL);
  ncm_catalog_set (test->catalog, "z", 3, -3.25, NULL);

  g_assert_cmpfloat (ncm_catalog_get (test->catalog, "ra", 0, NULL), ==, 1.5);
  g_assert_cmpfloat (ncm_catalog_get (test->catalog, "z", 3, NULL), ==, -3.25);
  g_assert_cmpfloat (ncm_catalog_get (test->catalog, "z", 0, NULL), ==, 0.0);

  /* The backing matrix reflects the same values. */
  data = ncm_catalog_peek_data (test->catalog);
  g_assert_cmpuint (ncm_matrix_nrows (data), ==, 4);
  g_assert_cmpuint (ncm_matrix_ncols (data), ==, 3);

  ncm_catalog_get_index (test->catalog, "ra", &j);
  g_assert_cmpfloat (ncm_matrix_get (data, 0, j), ==, 1.5);
  ncm_catalog_get_index (test->catalog, "z", &j);
  g_assert_cmpfloat (ncm_matrix_get (data, 3, j), ==, -3.25);
}

void
test_ncm_catalog_types (TestNcmCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[]            = {"id", "mass", "detected", NULL};
  const NcmCatalogColType col_types[] = {
    NCM_CATALOG_COL_TYPE_INT,
    NCM_CATALOG_COL_TYPE_DOUBLE,
    NCM_CATALOG_COL_TYPE_BOOL,
  };
  NcmCatalog *typed = ncm_catalog_new_full (2, (GStrv) col_names, col_types, 3);

  /* Plain new defaults every column to DOUBLE. */
  g_assert_cmpuint (ncm_catalog_get_col_type (test->catalog, "ra", NULL), ==, NCM_CATALOG_COL_TYPE_DOUBLE);

  /* new_full records the requested per-column types. */
  g_assert_cmpuint (ncm_catalog_get_col_type (typed, "id", NULL), ==, NCM_CATALOG_COL_TYPE_INT);
  g_assert_cmpuint (ncm_catalog_get_col_type (typed, "mass", NULL), ==, NCM_CATALOG_COL_TYPE_DOUBLE);
  g_assert_cmpuint (ncm_catalog_get_col_type (typed, "detected", NULL), ==, NCM_CATALOG_COL_TYPE_BOOL);

  /* Only the DOUBLE column is in the matrix. */
  g_assert_cmpuint (ncm_matrix_ncols (ncm_catalog_peek_data (typed)), ==, 1);
  g_assert_cmpuint (ncm_matrix_nrows (ncm_catalog_peek_data (typed)), ==, 2);

  /* INT values are exact above 2^53, where a double is not. */
  ncm_catalog_set_int (typed, "id", 0, TEST_NCM_CATALOG_BIG_ID, NULL);
  ncm_catalog_set_int (typed, "id", 1, TEST_NCM_CATALOG_BIG_ID + 1, NULL);
  ncm_catalog_set (typed, "mass", 1, 14.5, NULL);
  ncm_catalog_set_bool (typed, "detected", 0, TRUE, NULL);
  ncm_catalog_set_bool (typed, "detected", 1, FALSE, NULL);

  g_assert_cmpint (ncm_catalog_get_int (typed, "id", 0, NULL), ==, TEST_NCM_CATALOG_BIG_ID);
  g_assert_cmpint (ncm_catalog_get_int (typed, "id", 1, NULL), ==, TEST_NCM_CATALOG_BIG_ID + 1);
  g_assert_cmpfloat (ncm_catalog_get (typed, "mass", 1, NULL), ==, 14.5);
  g_assert_cmpfloat (ncm_matrix_get (ncm_catalog_peek_data (typed), 1, 0), ==, 14.5);
  g_assert_true (ncm_catalog_get_bool (typed, "detected", 0, NULL));
  g_assert_false (ncm_catalog_get_bool (typed, "detected", 1, NULL));

  NCM_TEST_FREE (ncm_catalog_free, typed);
}

/*
 * The column functions set and get every row of an INT or BOOL column; a length other
 * than the number of rows is an error, and a given output array is resized.
 */
void
test_ncm_catalog_columns (TestNcmCatalog *test, gconstpointer pdata)
{
  const gint64 ids[]     = {TEST_NCM_CATALOG_BIG_ID, -1, 0, G_MAXINT64};
  const gboolean flags[] = {TRUE, FALSE, FALSE, TRUE};
  GArray *in_ids         = g_array_new (FALSE, FALSE, sizeof (gint64));
  GArray *in_flags       = g_array_new (FALSE, FALSE, sizeof (gboolean));
  GArray *out_ids        = NULL;
  GArray *out_flags      = g_array_new (FALSE, TRUE, sizeof (gboolean));
  GError *error          = NULL;
  guint i;

  g_array_append_vals (in_ids, ids, 4);
  g_array_append_vals (in_flags, flags, 4);
  g_array_set_size (out_flags, 9);

  ncm_catalog_set_int_column (test->typed, "id", in_ids, &error);
  g_assert_no_error (error);
  ncm_catalog_set_bool_column (test->typed, "flag", in_flags, &error);
  g_assert_no_error (error);

  ncm_catalog_get_int_column (test->typed, "id", &out_ids, &error);
  g_assert_no_error (error);
  ncm_catalog_get_bool_column (test->typed, "flag", &out_flags, &error);
  g_assert_no_error (error);

  g_assert_cmpuint (out_ids->len, ==, 4);
  g_assert_cmpuint (out_flags->len, ==, 4);

  for (i = 0; i < 4; i++)
  {
    g_assert_cmpint (g_array_index (out_ids, gint64, i), ==, ids[i]);
    g_assert_cmpint (ncm_catalog_get_int (test->typed, "id", i, NULL), ==, ids[i]);
    g_assert_cmpint (g_array_index (out_flags, gboolean, i), ==, flags[i]);
    g_assert_cmpint (ncm_catalog_get_bool (test->typed, "flag", i, NULL), ==, flags[i]);
  }

  g_array_set_size (in_ids, 3);
  ncm_catalog_set_int_column (test->typed, "id", in_ids, &error);
  g_assert_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE);
  g_assert_cmpstr (error->message, ==, "ncm_catalog_set_int_column: 3 values for the 4 rows of column 'id'.");
  g_clear_error (&error);

  g_array_set_size (in_flags, 5);
  ncm_catalog_set_bool_column (test->typed, "flag", in_flags, &error);
  g_assert_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE);
  g_assert_cmpstr (error->message, ==, "ncm_catalog_set_bool_column: 5 values for the 4 rows of column 'flag'.");
  g_clear_error (&error);

  /* The failed writes left the columns unchanged. */
  g_assert_cmpint (ncm_catalog_get_int (test->typed, "id", 0, NULL), ==, TEST_NCM_CATALOG_BIG_ID);
  g_assert_true (ncm_catalog_get_bool (test->typed, "flag", 3, NULL));

  ncm_catalog_get_int_column (test->typed, "missing", &out_ids, &error);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");

  g_array_unref (in_ids);
  g_array_unref (in_flags);
  g_array_unref (out_ids);
  g_array_unref (out_flags);
}

/*
 * Every accessor acts only on columns of its type; a column of another type sets
 * NCM_CATALOG_ERROR_WRONG_TYPE and leaves the catalog unchanged.
 */
void
test_ncm_catalog_wrong_type (TestNcmCatalog *test, gconstpointer pdata)
{
  GArray *ints  = g_array_new (FALSE, TRUE, sizeof (gint64));
  GArray *bools = g_array_new (FALSE, TRUE, sizeof (gboolean));
  GError *error = NULL;
  guint64 key   = ncm_catalog_get_key (test->typed);
  gchar *hash   = ncm_catalog_get_data_hash (test->typed);
  gchar *hash_b;

  g_array_set_size (ints, 4);
  g_array_set_size (bools, 4);

  g_assert_true (isnan (ncm_catalog_get (test->typed, "id", 0, &error)));
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'id' is INT, not DOUBLE.");

  ncm_catalog_set (test->typed, "flag", 0, 1.0, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'flag' is BOOL, not DOUBLE.");

  g_assert_cmpint (ncm_catalog_get_int (test->typed, "ra", 0, &error), ==, 0);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'ra' is DOUBLE, not INT.");

  ncm_catalog_set_int (test->typed, "flag", 0, 1, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'flag' is BOOL, not INT.");

  g_assert_false (ncm_catalog_get_bool (test->typed, "id", 0, &error));
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'id' is INT, not BOOL.");

  ncm_catalog_set_bool (test->typed, "dec", 0, TRUE, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'dec' is DOUBLE, not BOOL.");

  ncm_catalog_set_int_column (test->typed, "ra", ints, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'ra' is DOUBLE, not INT.");

  ncm_catalog_get_int_column (test->typed, "flag", &ints, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'flag' is BOOL, not INT.");

  ncm_catalog_set_bool_column (test->typed, "id", bools, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'id' is INT, not BOOL.");

  ncm_catalog_get_bool_column (test->typed, "ra", &bools, &error);
  _test_ncm_catalog_assert_wrong_type_error (&error, "Column 'ra' is DOUBLE, not BOOL.");

  hash_b = ncm_catalog_get_data_hash (test->typed);
  g_assert_cmpuint (ncm_catalog_get_key (test->typed), ==, key);
  g_assert_cmpstr (hash_b, ==, hash);

  g_free (hash);
  g_free (hash_b);
  g_array_unref (ints);
  g_array_unref (bools);
}

void
test_ncm_catalog_meta (TestNcmCatalog *test, gconstpointer pdata)
{
  NcmVarDict *meta;
  NcmVarDict *got;
  gdouble val;

  /* A new catalog carries no metadata. */
  g_assert_null (ncm_catalog_peek_meta (test->catalog));

  meta = ncm_var_dict_new ();
  ncm_var_dict_set_double (meta, "cluster_z", 0.234);
  ncm_catalog_set_meta (test->catalog, meta);

  got = ncm_catalog_peek_meta (test->catalog);
  g_assert_true (got == meta);
  g_assert_true (ncm_var_dict_get_double (got, "cluster_z", &val));
  g_assert_cmpfloat (val, ==, 0.234);
  g_assert_false (ncm_var_dict_has_key (got, "missing"));

  /* Setting NULL clears it. */
  ncm_catalog_set_meta (test->catalog, NULL);
  g_assert_null (ncm_catalog_peek_meta (test->catalog));

  ncm_var_dict_unref (meta);
}

void
test_ncm_catalog_serialize (TestNcmCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[]            = {"id", "flag", "x", NULL};
  const NcmCatalogColType col_types[] = {
    NCM_CATALOG_COL_TYPE_INT,
    NCM_CATALOG_COL_TYPE_BOOL,
    NCM_CATALOG_COL_TYPE_DOUBLE,
  };
  NcmCatalog *typed = ncm_catalog_new_full (1, (GStrv) col_names, col_types, 3);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmVarDict *meta  = ncm_var_dict_new ();
  GVariant *var;
  NcmCatalog *dup;
  NcmVarDict *dup_meta;
  gdouble val;

  ncm_catalog_set_int (typed, "id", 0, TEST_NCM_CATALOG_BIG_ID, NULL);
  ncm_catalog_set_bool (typed, "flag", 0, TRUE, NULL);
  ncm_catalog_set (typed, "x", 0, -2.5, NULL);

  ncm_var_dict_set_double (meta, "cluster_ra", 30.4273);
  ncm_catalog_set_meta (typed, meta);

  var = ncm_serialize_to_variant (ser, G_OBJECT (typed));
  dup = NCM_CATALOG (ncm_serialize_from_variant (ser, var));

  g_assert_cmpuint (ncm_catalog_len (dup), ==, 1);
  g_assert_cmpuint (ncm_catalog_ncols (dup), ==, 3);
  g_assert_cmpuint (ncm_catalog_get_col_type (dup, "id", NULL), ==, NCM_CATALOG_COL_TYPE_INT);
  g_assert_cmpuint (ncm_catalog_get_col_type (dup, "flag", NULL), ==, NCM_CATALOG_COL_TYPE_BOOL);
  g_assert_cmpint (ncm_catalog_get_int (dup, "id", 0, NULL), ==, TEST_NCM_CATALOG_BIG_ID);
  g_assert_true (ncm_catalog_get_bool (dup, "flag", 0, NULL));
  g_assert_cmpfloat (ncm_catalog_get (dup, "x", 0, NULL), ==, -2.5);
  g_assert_cmpuint (ncm_matrix_ncols (ncm_catalog_peek_data (dup)), ==, 1);

  dup_meta = ncm_catalog_peek_meta (dup);
  g_assert_nonnull (dup_meta);
  g_assert_true (ncm_var_dict_get_double (dup_meta, "cluster_ra", &val));
  g_assert_cmpfloat (val, ==, 30.4273);

  g_variant_unref (var);
  ncm_serialize_free (ser);
  ncm_var_dict_unref (meta);
  NCM_TEST_FREE (ncm_catalog_free, typed);
  NCM_TEST_FREE (ncm_catalog_free, dup);
}

void
test_ncm_catalog_invalid_get (TestNcmCatalog *test, gconstpointer pdata)
{
  GError *error = NULL;

  g_assert_true (isnan (ncm_catalog_get (test->catalog, "missing", 0, &error)));
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");

  ncm_catalog_get_int (test->catalog, "missing", 0, &error);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");

  ncm_catalog_get_bool (test->catalog, "missing", 0, &error);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");
}

void
test_ncm_catalog_invalid_set (TestNcmCatalog *test, gconstpointer pdata)
{
  GError *error = NULL;

  ncm_catalog_set (test->catalog, "missing", 0, 1.0, &error);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");

  ncm_catalog_set_int (test->catalog, "missing", 0, 42, &error);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");

  ncm_catalog_set_bool (test->catalog, "missing", 0, TRUE, &error);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");
}

void
test_ncm_catalog_invalid_col_type (TestNcmCatalog *test, gconstpointer pdata)
{
  GError *error = NULL;

  g_assert_cmpint (ncm_catalog_get_col_type (test->catalog, "missing", &error), ==, NCM_CATALOG_COL_TYPE_INVALID);
  _test_ncm_catalog_assert_column_not_found_error (&error, "missing");
}

void
test_ncm_catalog_invalid_row (TestNcmCatalog *test, gconstpointer pdata)
{
  GError *error = NULL;

  /* The last row is fine, the next one is reported. */
  ncm_catalog_set (test->catalog, "z", 3, 1.5, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (ncm_catalog_get (test->catalog, "z", 3, &error), ==, 1.5);
  g_assert_no_error (error);

  g_assert_true (isnan (ncm_catalog_get (test->catalog, "z", 4, &error)));
  _test_ncm_catalog_assert_row_error (&error, "z");

  ncm_catalog_set (test->catalog, "z", 4, 1.0, &error);
  _test_ncm_catalog_assert_row_error (&error, "z");

  ncm_catalog_set_int (test->typed, "id", 3, 7, &error);
  g_assert_no_error (error);
  g_assert_cmpint (ncm_catalog_get_int (test->typed, "id", 3, &error), ==, 7);
  g_assert_no_error (error);

  g_assert_cmpint (ncm_catalog_get_int (test->typed, "id", 4, &error), ==, 0);
  _test_ncm_catalog_assert_row_error (&error, "id");

  ncm_catalog_set_int (test->typed, "id", 4, 1, &error);
  _test_ncm_catalog_assert_row_error (&error, "id");

  g_assert_false (ncm_catalog_get_bool (test->typed, "flag", 4, &error));
  _test_ncm_catalog_assert_row_error (&error, "flag");

  ncm_catalog_set_bool (test->typed, "flag", 4, TRUE, &error);
  _test_ncm_catalog_assert_row_error (&error, "flag");
}

void
test_ncm_catalog_data_prop (TestNcmCatalog *test, gconstpointer pdata)
{
  NcmMatrix *data = ncm_matrix_new (6, 3);

  /* A matrix with one column per name replaces the data and the length. */
  ncm_matrix_set_zero (data);
  ncm_matrix_set (data, 5, 2, 7.0);
  g_object_set (test->catalog, "data", data, NULL);

  g_assert_cmpuint (ncm_catalog_len (test->catalog), ==, 6);
  g_assert_cmpfloat (ncm_catalog_get (test->catalog, "z", 5, NULL), ==, 7.0);

  ncm_matrix_free (data);
}

/* A saved key is current until the next write; a pending one takes the current key at the first check. */
void
test_ncm_catalog_key_saved (TestNcmCatalog *test, gconstpointer pdata)
{
  guint64 saved   = ncm_catalog_get_key (test->catalog);
  guint64 pending = NCM_CATALOG_KEY_PENDING;

  g_assert_cmpuint (saved, !=, NCM_CATALOG_KEY_PENDING);
  g_assert_true (ncm_catalog_key_is_current (test->catalog, &saved));

  ncm_catalog_set (test->catalog, "ra", 0, 1.5, NULL);
  g_assert_false (ncm_catalog_key_is_current (test->catalog, &saved));
  g_assert_cmpuint (saved, !=, ncm_catalog_get_key (test->catalog));

  g_assert_true (ncm_catalog_key_is_current (test->catalog, &pending));
  g_assert_cmpuint (pending, ==, ncm_catalog_get_key (test->catalog));

  ncm_catalog_set (test->catalog, "ra", 1, 2.5, NULL);
  g_assert_false (ncm_catalog_key_is_current (test->catalog, &pending));
}

/* Every write through the setters, the data property or mark_outdated changes the key; reads and failed writes do not. */
void
test_ncm_catalog_key (TestNcmCatalog *test, gconstpointer pdata)
{
  NcmMatrix *data = ncm_matrix_new (2, 3);
  GError *error   = NULL;
  guint64 key     = ncm_catalog_get_key (test->catalog);

  ncm_catalog_get (test->catalog, "ra", 0, NULL);
  g_assert_cmpuint (ncm_catalog_get_key (test->catalog), ==, key);

  ncm_catalog_set (test->catalog, "ra", 0, 1.5, NULL);
  g_assert_cmpuint (ncm_catalog_get_key (test->catalog), >, key);
  key = ncm_catalog_get_key (test->catalog);

  {
    guint64 tkey = ncm_catalog_get_key (test->typed);
    GArray *ids  = g_array_new (FALSE, TRUE, sizeof (gint64));
    GArray *bv   = g_array_new (FALSE, TRUE, sizeof (gboolean));

    g_array_set_size (ids, 4);
    g_array_set_size (bv, 4);

    ncm_catalog_set_int (test->typed, "id", 1, 3, NULL);
    g_assert_cmpuint (ncm_catalog_get_key (test->typed), >, tkey);
    tkey = ncm_catalog_get_key (test->typed);

    ncm_catalog_set_bool (test->typed, "flag", 2, TRUE, NULL);
    g_assert_cmpuint (ncm_catalog_get_key (test->typed), >, tkey);
    tkey = ncm_catalog_get_key (test->typed);

    ncm_catalog_set_int_column (test->typed, "id", ids, NULL);
    g_assert_cmpuint (ncm_catalog_get_key (test->typed), >, tkey);
    tkey = ncm_catalog_get_key (test->typed);

    ncm_catalog_set_bool_column (test->typed, "flag", bv, NULL);
    g_assert_cmpuint (ncm_catalog_get_key (test->typed), >, tkey);
    tkey = ncm_catalog_get_key (test->typed);

    ncm_catalog_set_int (test->typed, "ra", 0, 1, &error);
    g_assert_nonnull (error);
    g_clear_error (&error);
    g_assert_cmpuint (ncm_catalog_get_key (test->typed), ==, tkey);

    g_array_unref (ids);
    g_array_unref (bv);
  }

  ncm_catalog_set (test->catalog, "nonexistent", 0, 1.0, &error);
  g_assert_nonnull (error);
  g_clear_error (&error);
  g_assert_cmpuint (ncm_catalog_get_key (test->catalog), ==, key);

  ncm_catalog_mark_outdated (test->catalog);
  g_assert_cmpuint (ncm_catalog_get_key (test->catalog), >, key);
  key = ncm_catalog_get_key (test->catalog);

  ncm_matrix_set_zero (data);
  g_object_set (test->catalog, "data", data, NULL);
  g_assert_cmpuint (ncm_catalog_get_key (test->catalog), >, key);

  ncm_matrix_free (data);
}

void
test_ncm_catalog_errors (void)
{
  const gchar *cases[][2] = {
    {"data_cols",  "*the catalog has 3 double columns, but the matrix has 2*"},
    {"data_null",  "*the data matrix of a catalog cannot be unset*"},
    {"old_layout", "*the data matrix has all 3 columns of the catalog*must be regenerated*"},
    {"rows",       "*the catalog has 4 rows and INT or BOOL columns, but the matrix has 5 rows*"},
    {"int_len",    "*the catalog has 4 rows and 1 INT columns, but the integer data has 3 values*"},
    {"bool_len",   "*the catalog has 4 rows and 1 BOOL columns, but the boolean data has 3 values*"},
    {"duplicate",  "*column `ra' appears more than once*"},
    {"types_len",  "*2 column types for 3 columns*"},
    {"types_val",  "*column `dec' has the invalid type 7*"},
    {"no_columns", "*a NcmCatalog needs its column names*"},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    g_setenv ("TEST_NCM_CATALOG_ERROR", cases[i][0], TRUE);
    g_test_trap_subprocess ("/ncm/catalog/errors/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (cases[i][1]);
  }

  g_unsetenv ("TEST_NCM_CATALOG_ERROR");
}

void
test_ncm_catalog_errors_subprocess (void)
{
  const gchar *which       = g_getenv ("TEST_NCM_CATALOG_ERROR");
  const gchar *col_names[] = {"ra", "dec", "z", NULL};

  if ((g_strcmp0 (which, "old_layout") == 0) || (g_strcmp0 (which, "rows") == 0) ||
      (g_strcmp0 (which, "int_len") == 0) || (g_strcmp0 (which, "bool_len") == 0))
  {
    const NcmCatalogColType int_types[]  = {NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_INT, NCM_CATALOG_COL_TYPE_DOUBLE};
    const NcmCatalogColType bool_types[] = {NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_BOOL, NCM_CATALOG_COL_TYPE_DOUBLE};
    const gint64 ints[]                  = {1, 2, 3};
    const guchar bools[]                 = {1, 0, 1};

    if (g_strcmp0 (which, "bool_len") == 0)
    {
      NcmCatalog *catalog = ncm_catalog_new_full (4, (GStrv) col_names, bool_types, 3);

      g_object_set (catalog, "bool-data", g_variant_new_fixed_array (G_VARIANT_TYPE_BOOLEAN, bools, 3, sizeof (guchar)), NULL);
    }
    else
    {
      NcmCatalog *catalog = ncm_catalog_new_full (4, (GStrv) col_names, int_types, 3);

      if (g_strcmp0 (which, "old_layout") == 0)
        g_object_set (catalog, "data", ncm_matrix_new (4, 3), NULL);
      else if (g_strcmp0 (which, "rows") == 0)
        g_object_set (catalog, "data", ncm_matrix_new (5, 2), NULL);
      else
        g_object_set (catalog, "int-data", g_variant_new_fixed_array (G_VARIANT_TYPE_INT64, ints, 3, sizeof (gint64)), NULL);
    }
  }
  else if (g_str_has_prefix (which, "data"))
  {
    NcmCatalog *catalog = ncm_catalog_new (4, (GStrv) col_names);

    if (g_strcmp0 (which, "data_cols") == 0)
      g_object_set (catalog, "data", ncm_matrix_new (4, 2), NULL);
    else
      g_object_set (catalog, "data", NULL, NULL);
  }
  else if (g_strcmp0 (which, "duplicate") == 0)
  {
    const gchar *dup_names[] = {"ra", "dec", "ra", NULL};

    ncm_catalog_new (4, (GStrv) dup_names);
  }
  else if (g_strcmp0 (which, "types_len") == 0)
  {
    const NcmCatalogColType types[] = {NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_INT};

    ncm_catalog_new_full (4, (GStrv) col_names, types, 2);
  }
  else if (g_strcmp0 (which, "types_val") == 0)
  {
    const NcmCatalogColType types[] = {NCM_CATALOG_COL_TYPE_DOUBLE, 7, NCM_CATALOG_COL_TYPE_INT};

    ncm_catalog_new_full (4, (GStrv) col_names, types, 3);
  }
  else
  {
    g_object_new (NCM_TYPE_CATALOG, NULL);
  }
}

/*
 * The data hash is 64 hexadecimal digits, equal for a serialized copy and for a change of
 * the metadata, and different after a change of one value or with other column names.
 */
void
test_ncm_catalog_data_hash (TestNcmCatalog *test, gconstpointer pdata)
{
  const gchar *other_names[] = {"ra", "dec", "zz", NULL};
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmCatalog *other          = ncm_catalog_new (4, (GStrv) other_names);
  NcmVarDict *meta           = ncm_var_dict_new ();
  gchar *hash, *hash_b;
  GVariant *var;
  NcmCatalog *dup;
  guint k;

  ncm_catalog_set (test->catalog, "z", 2, 0.5, NULL);
  hash = ncm_catalog_get_data_hash (test->catalog);
  g_assert_cmpuint (strlen (hash), ==, 64);

  for (k = 0; k < 64; k++)
    g_assert_true (g_ascii_isxdigit (hash[k]));

  var    = ncm_serialize_to_variant (ser, G_OBJECT (test->catalog));
  dup    = NCM_CATALOG (ncm_serialize_from_variant (ser, var));
  hash_b = ncm_catalog_get_data_hash (dup);
  g_assert_cmpstr (hash_b, ==, hash);
  g_free (hash_b);

  ncm_var_dict_set_double (meta, "cluster_ra", 30.4273);
  ncm_catalog_set_meta (dup, meta);
  hash_b = ncm_catalog_get_data_hash (dup);
  g_assert_cmpstr (hash_b, ==, hash);
  g_free (hash_b);

  ncm_catalog_set (dup, "z", 2, 0.5 + 1.0e-15, NULL);
  hash_b = ncm_catalog_get_data_hash (dup);
  g_assert_cmpstr (hash_b, !=, hash);
  g_free (hash_b);

  ncm_catalog_set (other, "zz", 2, 0.5, NULL);
  hash_b = ncm_catalog_get_data_hash (other);
  g_assert_cmpstr (hash_b, !=, hash);
  g_free (hash_b);

  g_free (hash);

  /* A change of an INT or of a BOOL value changes the hash. */
  hash = ncm_catalog_get_data_hash (test->typed);
  ncm_catalog_set_int (test->typed, "id", 3, TEST_NCM_CATALOG_BIG_ID, NULL);
  hash_b = ncm_catalog_get_data_hash (test->typed);
  g_assert_cmpstr (hash_b, !=, hash);
  g_free (hash);
  hash = hash_b;
  ncm_catalog_set_bool (test->typed, "flag", 0, TRUE, NULL);
  hash_b = ncm_catalog_get_data_hash (test->typed);
  g_assert_cmpstr (hash_b, !=, hash);
  g_free (hash_b);

  g_free (hash);
  g_variant_unref (var);
  ncm_var_dict_unref (meta);
  ncm_serialize_free (ser);
  NCM_TEST_FREE (ncm_catalog_free, dup);
  NCM_TEST_FREE (ncm_catalog_free, other);
}

/*
 * The data hash is the SHA-256 of the documented byte sequence: the names with their zero
 * byte, the types as 32-bit integers, the DOUBLE values row by row, then the INT and the
 * BOOL values row by row. For a catalog of DOUBLE columns only, it is the hash of the
 * names, the types and the rows of the matrix, as before INT and BOOL had their own
 * storage, so the interpolant files of existing catalogs keep matching them.
 */
void
test_ncm_catalog_data_hash_bytes (TestNcmCatalog *test, gconstpointer pdata)
{
  NcmCatalog *cats[2] = {test->catalog, test->typed};
  guint c;

  ncm_catalog_set (test->catalog, "dec", 1, 0.25, NULL);
  ncm_catalog_set (test->catalog, "z", 3, -1.5e3, NULL);
  ncm_catalog_set (test->typed, "dec", 2, 7.0, NULL);
  ncm_catalog_set_int (test->typed, "id", 1, TEST_NCM_CATALOG_BIG_ID, NULL);
  ncm_catalog_set_bool (test->typed, "flag", 3, TRUE, NULL);

  for (c = 0; c < 2; c++)
  {
    NcmCatalog *cat   = cats[c];
    GStrv cols        = ncm_catalog_peek_columns (cat);
    const guint ncols = ncm_catalog_ncols (cat);
    const guint len   = ncm_catalog_len (cat);
    GChecksum *sum    = g_checksum_new (G_CHECKSUM_SHA256);
    gchar *hash       = ncm_catalog_get_data_hash (cat);
    NcmCatalogColType t;
    guint i, j;

    for (j = 0; j < ncols; j++)
      g_checksum_update (sum, (const guchar *) cols[j], strlen (cols[j]) + 1);

    for (j = 0; j < ncols; j++)
    {
      const guint32 t32 = ncm_catalog_get_col_type (cat, cols[j], NULL);

      g_checksum_update (sum, (const guchar *) &t32, sizeof (guint32));
    }

    for (t = NCM_CATALOG_COL_TYPE_DOUBLE; t < NCM_CATALOG_COL_TYPE_LEN; t++)
    {
      for (i = 0; i < len; i++)
      {
        for (j = 0; j < ncols; j++)
        {
          if (ncm_catalog_get_col_type (cat, cols[j], NULL) != t)
            continue;

          switch (t)
          {
            case NCM_CATALOG_COL_TYPE_DOUBLE:
            {
              const gdouble v = ncm_catalog_get (cat, cols[j], i, NULL);

              g_checksum_update (sum, (const guchar *) &v, sizeof (gdouble));
              break;
            }
            case NCM_CATALOG_COL_TYPE_INT:
            {
              const gint64 v = ncm_catalog_get_int (cat, cols[j], i, NULL);

              g_checksum_update (sum, (const guchar *) &v, sizeof (gint64));
              break;
            }
            default:
            {
              const gboolean v = ncm_catalog_get_bool (cat, cols[j], i, NULL);

              g_checksum_update (sum, (const guchar *) &v, sizeof (gboolean));
              break;
            }
          }
        }
      }
    }

    g_assert_cmpstr (hash, ==, g_checksum_get_string (sum));

    g_checksum_free (sum);
    g_free (hash);
  }
}


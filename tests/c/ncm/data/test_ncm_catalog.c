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

typedef struct _TestNcmCatalog
{
  NcmCatalog *catalog;
} TestNcmCatalog;

static void _test_ncm_catalog_assert_column_not_found_error (GError **error, const gchar *col);

void test_ncm_catalog_new (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_free (TestNcmCatalog *test, gconstpointer pdata);

void test_ncm_catalog_ref (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_index (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_set_get (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_types (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_meta (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_serialize (TestNcmCatalog *test, gconstpointer pdata);

void test_ncm_catalog_invalid_get (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_invalid_set (TestNcmCatalog *test, gconstpointer pdata);
void test_ncm_catalog_invalid_col_type (TestNcmCatalog *test, gconstpointer pdata);

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

void
test_ncm_catalog_new (TestNcmCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[] = {"ra", "dec", "z", NULL};

  test->catalog = ncm_catalog_new (4, (GStrv) col_names);

  g_assert_true (NCM_IS_CATALOG (test->catalog));
  g_assert_cmpuint (ncm_catalog_len (test->catalog), ==, 4);
  g_assert_cmpuint (ncm_catalog_ncols (test->catalog), ==, 3);
}

void
test_ncm_catalog_free (TestNcmCatalog *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_catalog_free, test->catalog);
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

  /* Typed accessors round-trip through the double backing store. */
  ncm_catalog_set_int (typed, "id", 0, G_GINT64_CONSTANT (9007199254740991), NULL); /* 2^53 - 1 */
  ncm_catalog_set_bool (typed, "detected", 0, TRUE, NULL);
  ncm_catalog_set_bool (typed, "detected", 1, FALSE, NULL);

  g_assert_cmpint (ncm_catalog_get_int (typed, "id", 0, NULL), ==, G_GINT64_CONSTANT (9007199254740991));
  g_assert_true (ncm_catalog_get_bool (typed, "detected", 0, NULL));
  g_assert_false (ncm_catalog_get_bool (typed, "detected", 1, NULL));

  NCM_TEST_FREE (ncm_catalog_free, typed);
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
  const gchar *col_names[]            = {"id", "flag", NULL};
  const NcmCatalogColType col_types[] = {
    NCM_CATALOG_COL_TYPE_INT,
    NCM_CATALOG_COL_TYPE_BOOL,
  };
  NcmCatalog *typed = ncm_catalog_new_full (1, (GStrv) col_names, col_types, 2);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmVarDict *meta  = ncm_var_dict_new ();
  GVariant *var;
  NcmCatalog *dup;
  NcmVarDict *dup_meta;
  gdouble val;

  ncm_catalog_set_int (typed, "id", 0, 42, NULL);
  ncm_catalog_set_bool (typed, "flag", 0, TRUE, NULL);

  ncm_var_dict_set_double (meta, "cluster_ra", 30.4273);
  ncm_catalog_set_meta (typed, meta);

  var = ncm_serialize_to_variant (ser, G_OBJECT (typed));
  dup = NCM_CATALOG (ncm_serialize_from_variant (ser, var));

  g_assert_cmpuint (ncm_catalog_len (dup), ==, 1);
  g_assert_cmpuint (ncm_catalog_ncols (dup), ==, 2);
  g_assert_cmpuint (ncm_catalog_get_col_type (dup, "id", NULL), ==, NCM_CATALOG_COL_TYPE_INT);
  g_assert_cmpuint (ncm_catalog_get_col_type (dup, "flag", NULL), ==, NCM_CATALOG_COL_TYPE_BOOL);
  g_assert_cmpint (ncm_catalog_get_int (dup, "id", 0, NULL), ==, 42);
  g_assert_true (ncm_catalog_get_bool (dup, "flag", 0, NULL));

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


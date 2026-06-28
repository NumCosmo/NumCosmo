/***************************************************************************
 *            test_nc_halo_catalog.c
 *
 *  Sat Jun 13 16:30 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_halo_catalog.c
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

typedef struct _TestNcHaloCatalog
{
  NcHaloCatalog *hcat;
} TestNcHaloCatalog;

void test_nc_halo_catalog_new (TestNcHaloCatalog *test, gconstpointer pdata);
void test_nc_halo_catalog_free (TestNcHaloCatalog *test, gconstpointer pdata);

void test_nc_halo_catalog_ref (TestNcHaloCatalog *test, gconstpointer pdata);
void test_nc_halo_catalog_metadata (TestNcHaloCatalog *test, gconstpointer pdata);
void test_nc_halo_catalog_linkage (TestNcHaloCatalog *test, gconstpointer pdata);
void test_nc_halo_catalog_find_children (TestNcHaloCatalog *test, gconstpointer pdata);
void test_nc_halo_catalog_no_linkage (TestNcHaloCatalog *test, gconstpointer pdata);
void test_nc_halo_catalog_bad_column (TestNcHaloCatalog *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/halo_catalog/ref", TestNcHaloCatalog, NULL,
              &test_nc_halo_catalog_new,
              &test_nc_halo_catalog_ref,
              &test_nc_halo_catalog_free);

  g_test_add ("/nc/halo_catalog/metadata", TestNcHaloCatalog, NULL,
              &test_nc_halo_catalog_new,
              &test_nc_halo_catalog_metadata,
              &test_nc_halo_catalog_free);

  g_test_add ("/nc/halo_catalog/linkage", TestNcHaloCatalog, NULL,
              &test_nc_halo_catalog_new,
              &test_nc_halo_catalog_linkage,
              &test_nc_halo_catalog_free);

  g_test_add ("/nc/halo_catalog/find_children", TestNcHaloCatalog, NULL,
              &test_nc_halo_catalog_new,
              &test_nc_halo_catalog_find_children,
              &test_nc_halo_catalog_free);

  g_test_add ("/nc/halo_catalog/no_linkage", TestNcHaloCatalog, NULL,
              &test_nc_halo_catalog_new,
              &test_nc_halo_catalog_no_linkage,
              &test_nc_halo_catalog_free);

  g_test_add ("/nc/halo_catalog/bad_column", TestNcHaloCatalog, NULL,
              &test_nc_halo_catalog_new,
              &test_nc_halo_catalog_bad_column,
              &test_nc_halo_catalog_free);

  g_test_run ();

  return 0;
}

void
test_nc_halo_catalog_new (TestNcHaloCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[]            = {"cluster_id", "parent_id", "Mass", NULL};
  const NcmCatalogColType col_types[] = {
    NCM_CATALOG_COL_TYPE_INT,
    NCM_CATALOG_COL_TYPE_INT,
    NCM_CATALOG_COL_TYPE_DOUBLE,
  };
  const gint64 cluster_id[] = {100, 101, 102, 103};
  const gint64 parent_id[]  = {10, 10, 20, 0};
  guint i;

  test->hcat = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_CLUSTER, "cluster_id", "parent_id",
                                    4, (GStrv) col_names, col_types, 3);

  g_assert_true (NC_IS_HALO_CATALOG (test->hcat));
  g_assert_true (NCM_IS_CATALOG (test->hcat));

  for (i = 0; i < 4; i++)
  {
    ncm_catalog_set_int (NCM_CATALOG (test->hcat), "cluster_id", i, cluster_id[i], NULL);
    ncm_catalog_set_int (NCM_CATALOG (test->hcat), "parent_id", i, parent_id[i], NULL);
  }
}

void
test_nc_halo_catalog_free (TestNcHaloCatalog *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_halo_catalog_free, test->hcat);
}

void
test_nc_halo_catalog_ref (TestNcHaloCatalog *test, gconstpointer pdata)
{
  NcHaloCatalog *hcat_ref = nc_halo_catalog_ref (test->hcat);

  g_assert_true (hcat_ref == test->hcat);

  nc_halo_catalog_clear (&hcat_ref);
  g_assert_null (hcat_ref);

  g_assert_true (NC_IS_HALO_CATALOG (test->hcat));
}

void
test_nc_halo_catalog_metadata (TestNcHaloCatalog *test, gconstpointer pdata)
{
  g_assert_cmpuint (nc_halo_catalog_get_kind (test->hcat), ==, NC_HALO_CATALOG_KIND_CLUSTER);
  g_assert_cmpstr (nc_halo_catalog_peek_id_col (test->hcat), ==, "cluster_id");
  g_assert_cmpstr (nc_halo_catalog_peek_parent_id_col (test->hcat), ==, "parent_id");
  g_assert_cmpuint (ncm_catalog_len (NCM_CATALOG (test->hcat)), ==, 4);
}

void
test_nc_halo_catalog_linkage (TestNcHaloCatalog *test, gconstpointer pdata)
{
  GError *error = NULL;

  g_assert_cmpint (nc_halo_catalog_get_id (test->hcat, 0, &error), ==, 100);
  g_assert_no_error (error);
  g_assert_cmpint (nc_halo_catalog_get_parent_id (test->hcat, 0, &error), ==, 10);
  g_assert_no_error (error);
  g_assert_cmpint (nc_halo_catalog_get_id (test->hcat, 3, &error), ==, 103);
  g_assert_no_error (error);
  g_assert_cmpint (nc_halo_catalog_get_parent_id (test->hcat, 3, &error), ==, 0);
  g_assert_no_error (error);
}

void
test_nc_halo_catalog_find_children (TestNcHaloCatalog *test, gconstpointer pdata)
{
  GError *error    = NULL;
  GArray *children = nc_halo_catalog_find_children (test->hcat, 10, &error);

  g_assert_no_error (error);
  g_assert_cmpuint (children->len, ==, 2);
  g_assert_cmpuint (g_array_index (children, guint, 0), ==, 0);
  g_assert_cmpuint (g_array_index (children, guint, 1), ==, 1);
  g_array_unref (children);

  children = nc_halo_catalog_find_children (test->hcat, 20, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (children->len, ==, 1);
  g_assert_cmpuint (g_array_index (children, guint, 0), ==, 2);
  g_array_unref (children);

  children = nc_halo_catalog_find_children (test->hcat, 999, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (children->len, ==, 0);
  g_array_unref (children);
}

void
test_nc_halo_catalog_no_linkage (TestNcHaloCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[]            = {"Mass", NULL};
  const NcmCatalogColType col_types[] = {NCM_CATALOG_COL_TYPE_DOUBLE};
  NcHaloCatalog *hcat                 = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_HALO, NULL, NULL,
                                                            1, (GStrv) col_names, col_types, 1);
  GError *error = NULL;

  nc_halo_catalog_get_id (hcat, 0, &error);
  g_assert_error (error, NC_HALO_CATALOG_ERROR, NC_HALO_CATALOG_ERROR_NO_LINKAGE);
  g_clear_error (&error);

  nc_halo_catalog_get_parent_id (hcat, 0, &error);
  g_assert_error (error, NC_HALO_CATALOG_ERROR, NC_HALO_CATALOG_ERROR_NO_LINKAGE);
  g_clear_error (&error);

  g_assert_null (nc_halo_catalog_find_children (hcat, 0, &error));
  g_assert_error (error, NC_HALO_CATALOG_ERROR, NC_HALO_CATALOG_ERROR_NO_LINKAGE);
  g_clear_error (&error);

  NCM_TEST_FREE (nc_halo_catalog_free, hcat);
}

void
test_nc_halo_catalog_bad_column (TestNcHaloCatalog *test, gconstpointer pdata)
{
  const gchar *col_names[]            = {"Mass", NULL};
  const NcmCatalogColType col_types[] = {NCM_CATALOG_COL_TYPE_DOUBLE};
  /* id-col names a column that does not exist; the error surfaces from the base. */
  NcHaloCatalog *hcat = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_HALO, "missing", NULL,
                                             1, (GStrv) col_names, col_types, 1);
  GError *error = NULL;

  nc_halo_catalog_get_id (hcat, 0, &error);
  g_assert_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND);
  g_clear_error (&error);

  NCM_TEST_FREE (nc_halo_catalog_free, hcat);
}

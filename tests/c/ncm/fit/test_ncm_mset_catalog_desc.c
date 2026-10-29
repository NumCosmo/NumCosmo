/***************************************************************************
 *            test_ncm_mset_catalog_desc.c
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mset_catalog_desc.c
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
#include "build_cfg.h"
#include <numcosmo/numcosmo.h>

#define TEST_DESC_DIM 3
#define TEST_DESC_NCHAINS 4
#define TEST_DESC_RUN_TYPE "test-run"

typedef struct _TestNcmMSetCatalogDesc
{
  NcmModelMVND *model_mvnd;
  NcmMSet *mset;
  NcmMSetCatalogDesc *desc;
} TestNcmMSetCatalogDesc;

static const gchar *test_names[]   = {"m2lnL", "aux", NULL};
static const gchar *test_symbols[] = {"-2\\ln(L)", "a", NULL};

void test_ncm_mset_catalog_desc_new (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_free (TestNcmMSetCatalogDesc *test, gconstpointer pdata);

void test_ncm_mset_catalog_desc_columns (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_mset_copy (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_serialize (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_col_by_name (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_new_errors (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_incompatible (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_order (TestNcmMSetCatalogDesc *test, gconstpointer pdata);
void test_ncm_mset_catalog_desc_property_traps (void);
void test_ncm_mset_catalog_desc_property_subprocess (void);

static NcmMSet *_test_ncm_mset_catalog_desc_mset_new (NcmModelMVND *model_mvnd);
static void _test_ncm_mset_catalog_desc_assert_incompatible (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/mset/catalog/desc/columns", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_columns,
              &test_ncm_mset_catalog_desc_free);
  g_test_add ("/ncm/mset/catalog/desc/mset_copy", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_mset_copy,
              &test_ncm_mset_catalog_desc_free);
  g_test_add ("/ncm/mset/catalog/desc/serialize", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_serialize,
              &test_ncm_mset_catalog_desc_free);
  g_test_add ("/ncm/mset/catalog/desc/col_by_name", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_col_by_name,
              &test_ncm_mset_catalog_desc_free);
  g_test_add ("/ncm/mset/catalog/desc/new_errors", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_new_errors,
              &test_ncm_mset_catalog_desc_free);
  g_test_add ("/ncm/mset/catalog/desc/incompatible", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_incompatible,
              &test_ncm_mset_catalog_desc_free);
  g_test_add ("/ncm/mset/catalog/desc/order", TestNcmMSetCatalogDesc, NULL,
              &test_ncm_mset_catalog_desc_new,
              &test_ncm_mset_catalog_desc_order,
              &test_ncm_mset_catalog_desc_free);

  g_test_add_func ("/ncm/mset/catalog/desc/property/traps", &test_ncm_mset_catalog_desc_property_traps);
  g_test_add_func ("/ncm/mset/catalog/desc/property/subprocess", &test_ncm_mset_catalog_desc_property_subprocess);

  g_test_run ();
}

/* An mset with one TEST_DESC_DIM-dimensional NcmModelMVND, every parameter free, map prepared. */
static NcmMSet *
_test_ncm_mset_catalog_desc_mset_new (NcmModelMVND *model_mvnd)
{
  NcmMSet *mset = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  return mset;
}

static void
_test_ncm_mset_catalog_desc_assert_incompatible (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other)
{
  GError *error = NULL;

  g_assert_false (ncm_mset_catalog_desc_check_compatible (desc, other, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE);
  g_clear_error (&error);

  g_assert_false (ncm_mset_catalog_desc_check_compatible (other, desc, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE);
  g_clear_error (&error);
}

void
test_ncm_mset_catalog_desc_new (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  GError *error = NULL;

  test->model_mvnd = ncm_model_mvnd_new (TEST_DESC_DIM);
  test->mset       = _test_ncm_mset_catalog_desc_mset_new (test->model_mvnd);
  test->desc       = ncm_mset_catalog_desc_new (test->mset, TEST_DESC_NCHAINS, TRUE,
                                                test_names, test_symbols, 0,
                                                TEST_DESC_RUN_TYPE, &error);

  g_assert_no_error (error);
  g_assert_true (NCM_IS_MSET_CATALOG_DESC (test->desc));
}

void
test_ncm_mset_catalog_desc_free (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_mset_catalog_desc_free, test->desc);
  ncm_mset_clear (&test->mset);
  ncm_model_mvnd_clear (&test->model_mvnd);
}

void
test_ncm_mset_catalog_desc_columns (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  NcmMSetCatalogDesc *desc = test->desc;
  NcmMSet *mset            = ncm_mset_catalog_desc_peek_mset (desc);
  guint i;

  g_assert_cmpuint (ncm_mset_catalog_desc_nchains (desc), ==, TEST_DESC_NCHAINS);
  g_assert_true (ncm_mset_catalog_desc_weighted (desc));
  g_assert_cmpint (ncm_mset_catalog_desc_get_m2lnp_var (desc), ==, 0);
  g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (desc), ==, TEST_DESC_RUN_TYPE);

  /* The user's two additional values, then the weight column, then the free parameters. */
  g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (desc), ==, 3);
  g_assert_cmpuint (ncm_mset_catalog_desc_ncols (desc), ==, 3 + TEST_DESC_DIM);

  g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 0), ==, "m2lnL");
  g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 1), ==, "aux");
  g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 2), ==, NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME);
  g_assert_cmpstr (ncm_mset_catalog_desc_col_symb (desc, 0), ==, "-2\\ln(L)");
  g_assert_cmpstr (ncm_mset_catalog_desc_col_symb (desc, 1), ==, "a");
  g_assert_cmpstr (ncm_mset_catalog_desc_col_symb (desc, 2), ==, NCM_MSET_CATALOG_DESC_WEIGHT_SYMBOL);
  g_assert_cmpstr (ncm_mset_catalog_desc_col_full_name (desc, 2), ==, NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME);

  for (i = 0; i < TEST_DESC_DIM; i++)
  {
    g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 3 + i), ==, ncm_mset_fparam_name (mset, i));
    g_assert_cmpstr (ncm_mset_catalog_desc_col_full_name (desc, 3 + i), ==, ncm_mset_fparam_full_name (mset, i));
    g_assert_cmpstr (ncm_mset_catalog_desc_col_symb (desc, 3 + i), ==, ncm_mset_fparam_symbol (mset, i));
  }

  /* Unweighted: no weight column. */
  {
    GError *error               = NULL;
    NcmMSetCatalogDesc *desc_nw = ncm_mset_catalog_desc_new (test->mset, 1, FALSE, test_names, test_symbols,
                                                             -1, NULL, &error);

    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (desc_nw), ==, 2);
    g_assert_cmpuint (ncm_mset_catalog_desc_ncols (desc_nw), ==, 2 + TEST_DESC_DIM);
    g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (desc_nw), ==, NCM_MSET_CATALOG_RTYPE_UNDEFINED);

    ncm_mset_catalog_desc_free (desc_nw);
  }
}

void
test_ncm_mset_catalog_desc_mset_copy (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  NcmMSetCatalogDesc *desc = test->desc;
  NcmMSet *mset_copy       = ncm_mset_catalog_desc_peek_mset (desc);
  const gdouble value0     = ncm_mset_fparam_get (mset_copy, 0);

  g_assert_true (mset_copy != test->mset);
  g_assert_true (ncm_mset_catalog_desc_peek_mset (desc) == mset_copy);

  /* Changing the caller's mset reaches neither the copy's values nor its free parameters. */
  ncm_mset_fparam_set (test->mset, 0, value0 + 1.0);
  ncm_model_param_set_ftype (NCM_MODEL (test->model_mvnd), 0, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);

  g_assert_cmpuint (ncm_mset_fparams_len (test->mset), ==, TEST_DESC_DIM - 1);
  g_assert_cmpuint (ncm_mset_fparams_len (mset_copy), ==, TEST_DESC_DIM);
  g_assert_cmpfloat (ncm_mset_fparam_get (mset_copy, 0), ==, value0);
  g_assert_cmpuint (ncm_mset_catalog_desc_ncols (desc), ==, 3 + TEST_DESC_DIM);
}

void
test_ncm_mset_catalog_desc_serialize (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  NcmSerialize *ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmMSetCatalogDesc *desc = test->desc;
  NcmMSetCatalogDesc *dup  = NCM_MSET_CATALOG_DESC (ncm_serialize_dup_obj (ser, G_OBJECT (desc)));
  GError *error            = NULL;
  gchar **names            = NULL;
  guint i;

  /* The properties hold the user's values only, so the copy adds the weight column once. */
  g_object_get (dup, "nadd-val-names", &names, NULL);
  g_assert_cmpuint (g_strv_length (names), ==, 2);
  g_strfreev (names);

  g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (dup), ==, ncm_mset_catalog_desc_nadd_vals (desc));
  g_assert_cmpuint (ncm_mset_catalog_desc_ncols (dup), ==, ncm_mset_catalog_desc_ncols (desc));
  g_assert_cmpuint (ncm_mset_catalog_desc_nchains (dup), ==, TEST_DESC_NCHAINS);
  g_assert_true (ncm_mset_catalog_desc_weighted (dup));
  g_assert_cmpint (ncm_mset_catalog_desc_get_m2lnp_var (dup), ==, 0);
  g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (dup), ==, TEST_DESC_RUN_TYPE);

  for (i = 0; i < ncm_mset_catalog_desc_ncols (desc); i++)
  {
    g_assert_cmpstr (ncm_mset_catalog_desc_col_full_name (dup, i), ==, ncm_mset_catalog_desc_col_full_name (desc, i));
    g_assert_cmpstr (ncm_mset_catalog_desc_col_symb (dup, i), ==, ncm_mset_catalog_desc_col_symb (desc, i));
  }

  g_assert_true (ncm_mset_catalog_desc_check_compatible (desc, dup, &error));
  g_assert_no_error (error);
  g_assert_true (ncm_mset_catalog_desc_check_same_order (desc, dup, &error));
  g_assert_no_error (error);

  ncm_mset_catalog_desc_free (dup);
  ncm_serialize_free (ser);
}

void
test_ncm_mset_catalog_desc_col_by_name (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  NcmMSetCatalogDesc *desc = test->desc;
  NcmMSet *mset            = ncm_mset_catalog_desc_peek_mset (desc);
  guint col                = G_MAXUINT;

  g_assert_true (ncm_mset_catalog_desc_col_by_name (desc, "m2lnL", &col));
  g_assert_cmpuint (col, ==, 0);
  g_assert_true (ncm_mset_catalog_desc_col_by_name (desc, NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME, &col));
  g_assert_cmpuint (col, ==, 2);
  g_assert_true (ncm_mset_catalog_desc_col_by_name (desc, ncm_mset_fparam_full_name (mset, 1), &col));
  g_assert_cmpuint (col, ==, 3 + 1);
  g_assert_true (ncm_mset_catalog_desc_col_by_name (desc, "1", &col));
  g_assert_cmpuint (col, ==, 1);
  g_assert_false (ncm_mset_catalog_desc_col_by_name (desc, "no-such-column", &col));
}

void
test_ncm_mset_catalog_desc_new_errors (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  const gchar *one_symbol[] = {"-2\\ln(L)", NULL};
  NcmMSet *mset_nomap       = ncm_mset_new (NCM_MODEL (test->model_mvnd), NULL, NULL);
  GError *error             = NULL;

  g_assert_false (ncm_mset_fparam_map_valid (mset_nomap));

  g_assert_null (ncm_mset_catalog_desc_new (NULL, 1, FALSE, NULL, NULL, -1, NULL, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_MSET);
  g_clear_error (&error);

  g_assert_null (ncm_mset_catalog_desc_new (mset_nomap, 1, FALSE, NULL, NULL, -1, NULL, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_MSET);
  g_clear_error (&error);

  g_assert_null (ncm_mset_catalog_desc_new (test->mset, 1, FALSE, test_names, one_symbol, -1, NULL, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_ADD_VALS);
  g_clear_error (&error);

  /* m2lnp-var must be -1 or a user additional value; the weight column does not count. */
  g_assert_null (ncm_mset_catalog_desc_new (test->mset, 1, TRUE, test_names, test_symbols, 2, NULL, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_M2LNP_VAR);
  g_clear_error (&error);

  g_assert_null (ncm_mset_catalog_desc_new (test->mset, 1, FALSE, test_names, test_symbols, -2, NULL, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_M2LNP_VAR);
  g_clear_error (&error);

  ncm_mset_free (mset_nomap);
}

void
test_ncm_mset_catalog_desc_incompatible (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  const gchar *other_names[]   = {"m2lnL", "other", NULL};
  const gchar *other_symbols[] = {"-2\\ln(L)", "b", NULL};
  NcmMSetCatalogDesc *desc     = test->desc;
  NcmMSet *mset                = test->mset;
  GError *error                = NULL;
  NcmMSetCatalogDesc *other;

  /* Equal inputs are compatible. */
  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, TRUE, test_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset_catalog_desc_check_compatible (desc, other, &error));
  g_assert_no_error (error);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS + 1, TRUE, test_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, FALSE, test_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, TRUE, test_names, test_symbols, 1, TEST_DESC_RUN_TYPE, &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, TRUE, test_names, test_symbols, 0, "other-run", &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, TRUE, other_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, TRUE, test_names, other_symbols, 0, TEST_DESC_RUN_TYPE, &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (mset, TEST_DESC_NCHAINS, TRUE, NULL, NULL, -1, TEST_DESC_RUN_TYPE, &error);
  _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);
  ncm_mset_catalog_desc_free (other);

  /* Parameters 1 and 2 free against all three free, then against parameters 0 and 2 free. */
  {
    NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_DESC_DIM);
    NcmMSet *mset_free       = _test_ncm_mset_catalog_desc_mset_new (model_mvnd);
    NcmMSetCatalogDesc *other_02;

    ncm_model_param_set_ftype (NCM_MODEL (model_mvnd), 0, NCM_PARAM_TYPE_FIXED);
    ncm_mset_prepare_fparam_map (mset_free);

    other = ncm_mset_catalog_desc_new (mset_free, TEST_DESC_NCHAINS, TRUE, test_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
    g_assert_no_error (error);
    _test_ncm_mset_catalog_desc_assert_incompatible (desc, other);

    ncm_model_param_set_ftype (NCM_MODEL (model_mvnd), 0, NCM_PARAM_TYPE_FREE);
    ncm_model_param_set_ftype (NCM_MODEL (model_mvnd), 1, NCM_PARAM_TYPE_FIXED);
    ncm_mset_prepare_fparam_map (mset_free);

    other_02 = ncm_mset_catalog_desc_new (mset_free, TEST_DESC_NCHAINS, TRUE, test_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_desc_ncols (other_02), ==, ncm_mset_catalog_desc_ncols (other));
    _test_ncm_mset_catalog_desc_assert_incompatible (other, other_02);

    ncm_mset_catalog_desc_free (other_02);
    ncm_mset_catalog_desc_free (other);
    ncm_mset_free (mset_free);
    ncm_model_mvnd_free (model_mvnd);
  }
}

void
test_ncm_mset_catalog_desc_order (TestNcmMSetCatalogDesc *test, gconstpointer pdata)
{
  NcmMSetCatalogDesc *desc = test->desc;
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_DESC_DIM);
  NcmMSet *mset_rev        = _test_ncm_mset_catalog_desc_mset_new (model_mvnd);
  gchar **fmap             = ncm_mset_get_fmap (mset_rev);
  GError *error            = NULL;
  NcmMSetCatalogDesc *other;
  guint i;

  /* The same free parameters, in reverse order. */
  for (i = 0; i < TEST_DESC_DIM / 2; i++)
  {
    gchar *tmp = fmap[i];

    fmap[i]                     = fmap[TEST_DESC_DIM - 1 - i];
    fmap[TEST_DESC_DIM - 1 - i] = tmp;
  }

  ncm_mset_set_fmap (mset_rev, (const gchar * const *) fmap, FALSE, &error);
  g_assert_no_error (error);

  other = ncm_mset_catalog_desc_new (mset_rev, TEST_DESC_NCHAINS, TRUE, test_names, test_symbols, 0, TEST_DESC_RUN_TYPE, &error);
  g_assert_no_error (error);
  g_assert_cmpstr (ncm_mset_catalog_desc_col_full_name (other, 3), ==, ncm_mset_catalog_desc_col_full_name (desc, 3 + TEST_DESC_DIM - 1));

  g_assert_true (ncm_mset_catalog_desc_check_compatible (desc, other, &error));
  g_assert_no_error (error);

  g_assert_false (ncm_mset_catalog_desc_check_same_order (desc, other, &error));
  g_assert_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_ORDER);
  g_clear_error (&error);

  g_assert_true (ncm_mset_catalog_desc_check_same_order (desc, desc, &error));
  g_assert_no_error (error);

  ncm_mset_catalog_desc_free (other);
  g_strfreev (fmap);
  ncm_mset_free (mset_rev);
  ncm_model_mvnd_free (model_mvnd);
}

void
test_ncm_mset_catalog_desc_property_traps (void)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/desc/property/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the free-parameter map of the NcmMSet is not prepared*");
}

void
test_ncm_mset_catalog_desc_property_subprocess (void)
{
  /* Built from properties there is no error location, so an invalid input aborts. */
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_DESC_DIM);
  NcmMSet *mset_nomap      = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);

  g_object_new (NCM_TYPE_MSET_CATALOG_DESC, "mset", mset_nomap, NULL);

  g_assert_not_reached ();
}


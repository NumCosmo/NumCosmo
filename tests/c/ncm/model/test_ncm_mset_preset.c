/***************************************************************************
 *            test_ncm_mset_preset.c
 *
 *  Thu Nov 19 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mset_preset.c
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
#include <numcosmo/ncm/model/ncm_reparam_linear.h>

#include <math.h>
#include <string.h>
#include <glib.h>
#include <glib-object.h>
#include <glib/gstdio.h>

#define TEST_X1 "NcmModelRosenbrock:x1"
#define TEST_X2 "NcmModelRosenbrock:x2"

typedef struct _TestNcmMSetPreset
{
  NcmModel *model;
  NcmMSet *mset;
} TestNcmMSetPreset;

/* A preset file with comments, in the form of ncm_serialize_to_yaml(). */
static const gchar *test_yaml =
  "# Test reference, table 1.\n"
  "NcmMSetPreset:\n"
  "  desc: 'Test reference, table 1'\n"
  "  data-combination: 'A+B'\n"
  "  statistic: 'best fit'\n"
  "  value:\n"
  "    NcmModelRosenbrock:x1: 1.5   # column 2\n"
  "  scale:\n"
  "    NcmModelRosenbrock:x1: 0.25  # column 3\n"
  "  lower-bound:\n"
  "    NcmModelRosenbrock:x2: -3.0\n"
  "  fit:\n"
  "    NcmModelRosenbrock:x2: false\n";

static void test_ncm_mset_preset_new (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_free (TestNcmMSetPreset *test, gconstpointer pdata);

static void test_ncm_mset_preset_apply (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_serialize (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_from_yaml_comments (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_reparam (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_error_param_not_found (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_error_invalid_entry (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_from_name (TestNcmMSetPreset *test, gconstpointer pdata);
static void test_ncm_mset_preset_shipped (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/mset_preset/apply", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_apply,
              &test_ncm_mset_preset_free);

  g_test_add ("/ncm/mset_preset/serialize", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_serialize,
              &test_ncm_mset_preset_free);

  g_test_add ("/ncm/mset_preset/from_yaml_comments", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_from_yaml_comments,
              &test_ncm_mset_preset_free);

  g_test_add ("/ncm/mset_preset/reparam", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_reparam,
              &test_ncm_mset_preset_free);

  g_test_add ("/ncm/mset_preset/error/param_not_found", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_error_param_not_found,
              &test_ncm_mset_preset_free);

  g_test_add ("/ncm/mset_preset/error/invalid_entry", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_error_invalid_entry,
              &test_ncm_mset_preset_free);

  g_test_add ("/ncm/mset_preset/from_name", TestNcmMSetPreset, NULL,
              &test_ncm_mset_preset_new,
              &test_ncm_mset_preset_from_name,
              &test_ncm_mset_preset_free);

  g_test_add_func ("/ncm/mset_preset/shipped", &test_ncm_mset_preset_shipped);

  g_test_run ();
}

/* Asserts that the attributes of the two Rosenbrock parameters are their defaults. */
static void
_test_assert_unchanged (NcmModel *model)
{
  g_assert_cmpfloat (ncm_model_param_get (model, NCM_MODEL_ROSENBROCK_X1), ==, 0.0);
  g_assert_cmpfloat (ncm_model_param_get_scale (model, NCM_MODEL_ROSENBROCK_X1), ==, 1.0e-1);
  g_assert_cmpfloat (ncm_model_param_get_lower_bound (model, NCM_MODEL_ROSENBROCK_X1), ==, -200.0);
  g_assert_cmpfloat (ncm_model_param_get_upper_bound (model, NCM_MODEL_ROSENBROCK_X1), ==, 200.0);
  g_assert_cmpfloat (ncm_model_param_get_lower_bound (model, NCM_MODEL_ROSENBROCK_X2), ==, -400.0);
  g_assert_cmpint (ncm_model_param_get_ftype (model, NCM_MODEL_ROSENBROCK_X2), ==, NCM_PARAM_TYPE_FREE);
}

/* Asserts the attributes set by the preset of test_yaml. */
static void
_test_assert_test_yaml_applied (NcmModel *model)
{
  g_assert_cmpfloat (ncm_model_param_get (model, NCM_MODEL_ROSENBROCK_X1), ==, 1.5);
  g_assert_cmpfloat (ncm_model_param_get_scale (model, NCM_MODEL_ROSENBROCK_X1), ==, 0.25);
  g_assert_cmpfloat (ncm_model_param_get_lower_bound (model, NCM_MODEL_ROSENBROCK_X2), ==, -3.0);
  g_assert_cmpfloat (ncm_model_param_get_upper_bound (model, NCM_MODEL_ROSENBROCK_X2), ==, 800.0);
  g_assert_cmpint (ncm_model_param_get_ftype (model, NCM_MODEL_ROSENBROCK_X2), ==, NCM_PARAM_TYPE_FIXED);
}

/* Applies @preset, asserts @code in NCM_MSET_PRESET_ERROR and that @model is unchanged. */
static void
_test_assert_apply_error (NcmMSetPreset *preset, NcmMSet *mset, NcmModel *model, gint code)
{
  GError *error = NULL;

  ncm_mset_preset_apply (preset, mset, &error);
  g_assert_error (error, NCM_MSET_PRESET_ERROR, code);
  g_assert_true (g_str_has_prefix (error->message, "ncm_mset_preset_apply: "));
  g_clear_error (&error);

  _test_assert_unchanged (model);
}

static void
_test_write_file (const gchar *filename, const gchar *contents)
{
  GError *error = NULL;

  g_assert_true (g_file_set_contents (filename, contents, -1, &error));
  g_assert_no_error (error);
}

static void
test_ncm_mset_preset_new (TestNcmMSetPreset *test, gconstpointer pdata)
{
  test->model = NCM_MODEL (ncm_model_rosenbrock_new ());
  test->mset  = ncm_mset_new (test->model, NULL, NULL);

  /* Both parameters free, as their NcmSParam defaults. */
  ncm_model_params_set_default_ftype (test->model);
}

static void
test_ncm_mset_preset_free (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_mset_free, test->mset);
  NCM_TEST_FREE (ncm_model_free, test->model);
}

static void
test_ncm_mset_preset_apply (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NcmMSetPreset *preset = ncm_mset_preset_new ("Test reference");
  GError *error         = NULL;

  ncm_mset_preset_set_value (preset, TEST_X1, 2.0);
  ncm_mset_preset_set_scale (preset, TEST_X1, 0.5);
  ncm_mset_preset_set_lower_bound (preset, TEST_X1, -10.0);
  ncm_mset_preset_set_upper_bound (preset, TEST_X1, 10.0);
  ncm_mset_preset_set_abstol (preset, TEST_X1, 1.0e-3);
  ncm_mset_preset_set_fit (preset, TEST_X1, FALSE);
  ncm_mset_preset_set_scale (preset, TEST_X2, 3.0);

  g_assert_cmpstr (ncm_mset_preset_peek_desc (preset), ==, "Test reference");
  g_assert_null (ncm_mset_preset_peek_data_combination (preset));
  g_assert_null (ncm_mset_preset_peek_statistic (preset));

  ncm_mset_preset_apply (preset, test->mset, &error);
  g_assert_no_error (error);

  g_assert_cmpfloat (ncm_model_param_get (test->model, NCM_MODEL_ROSENBROCK_X1), ==, 2.0);
  g_assert_cmpfloat (ncm_model_param_get_scale (test->model, NCM_MODEL_ROSENBROCK_X1), ==, 0.5);
  g_assert_cmpfloat (ncm_model_param_get_lower_bound (test->model, NCM_MODEL_ROSENBROCK_X1), ==, -10.0);
  g_assert_cmpfloat (ncm_model_param_get_upper_bound (test->model, NCM_MODEL_ROSENBROCK_X1), ==, 10.0);
  g_assert_cmpfloat (ncm_model_param_get_abstol (test->model, NCM_MODEL_ROSENBROCK_X1), ==, 1.0e-3);
  g_assert_cmpint (ncm_model_param_get_ftype (test->model, NCM_MODEL_ROSENBROCK_X1), ==, NCM_PARAM_TYPE_FIXED);

  /* x2 appears only in scale: its other attributes are unchanged. */
  g_assert_cmpfloat (ncm_model_param_get_scale (test->model, NCM_MODEL_ROSENBROCK_X2), ==, 3.0);
  g_assert_cmpfloat (ncm_model_param_get (test->model, NCM_MODEL_ROSENBROCK_X2), ==, 0.0);
  g_assert_cmpfloat (ncm_model_param_get_lower_bound (test->model, NCM_MODEL_ROSENBROCK_X2), ==, -400.0);
  g_assert_cmpfloat (ncm_model_param_get_upper_bound (test->model, NCM_MODEL_ROSENBROCK_X2), ==, 800.0);
  g_assert_cmpint (ncm_model_param_get_ftype (test->model, NCM_MODEL_ROSENBROCK_X2), ==, NCM_PARAM_TYPE_FREE);

  NCM_TEST_FREE (ncm_mset_preset_free, preset);
}

static void
test_ncm_mset_preset_serialize (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmMSetPreset *preset = ncm_mset_preset_new ("Test reference");
  GError *error         = NULL;
  NcmMSetPreset *dup;
  gchar *yaml;

  g_object_set (preset, "data-combination", "A+B", "statistic", "best fit", NULL);
  ncm_mset_preset_set_value (preset, TEST_X1, 1.5);
  ncm_mset_preset_set_scale (preset, TEST_X1, 0.25);
  ncm_mset_preset_set_lower_bound (preset, TEST_X2, -3.0);
  ncm_mset_preset_set_fit (preset, TEST_X2, FALSE);

  yaml = ncm_serialize_to_yaml (ser, G_OBJECT (preset));
  dup  = NCM_MSET_PRESET (ncm_serialize_from_yaml (ser, yaml));

  g_assert_cmpstr (ncm_mset_preset_peek_desc (dup), ==, "Test reference");
  g_assert_cmpstr (ncm_mset_preset_peek_data_combination (dup), ==, "A+B");
  g_assert_cmpstr (ncm_mset_preset_peek_statistic (dup), ==, "best fit");

  ncm_mset_preset_apply (dup, test->mset, &error);
  g_assert_no_error (error);
  _test_assert_test_yaml_applied (test->model);

  g_free (yaml);
  ncm_serialize_free (ser);
  NCM_TEST_FREE (ncm_mset_preset_free, dup);
  NCM_TEST_FREE (ncm_mset_preset_free, preset);
}

static void
test_ncm_mset_preset_from_yaml_comments (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmMSetPreset *preset = NCM_MSET_PRESET (ncm_serialize_from_yaml (ser, test_yaml));
  GError *error         = NULL;

  g_assert_cmpstr (ncm_mset_preset_peek_desc (preset), ==, "Test reference, table 1");

  ncm_mset_preset_apply (preset, test->mset, &error);
  g_assert_no_error (error);
  _test_assert_test_yaml_applied (test->model);

  ncm_serialize_free (ser);
  NCM_TEST_FREE (ncm_mset_preset_free, preset);
}

static void
test_ncm_mset_preset_reparam (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NcmMSetPreset *preset = ncm_mset_preset_new ("Test reference");
  NcmModel *model       = NCM_MODEL (ncm_model_rosenbrock_new ());
  NcmMSet *mset         = ncm_mset_new (model, NULL, NULL);
  NcmMatrix *T          = ncm_matrix_new (2, 2);
  NcmVector *v          = ncm_vector_new (2);
  GError *error         = NULL;
  NcmReparamLinear *relin;

  /* The identity map, with x1 renamed to y1. */
  ncm_matrix_set_identity (T);
  ncm_vector_set_zero (v);
  relin = ncm_reparam_linear_new (2, T, v);
  ncm_reparam_set_compat_type (NCM_REPARAM (relin), NCM_TYPE_MODEL_ROSENBROCK);
  ncm_reparam_set_param_desc_full (NCM_REPARAM (relin), 0, "y1", "y_1", -10.0, 10.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  ncm_model_set_reparam (model, NCM_REPARAM (relin), &error);
  g_assert_no_error (error);

  ncm_mset_preset_set_value (preset, "NcmModelRosenbrock:y1", 2.0);

  /* The model has a reparametrization the preset does not name. */
  ncm_mset_preset_apply (preset, mset, &error);
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH);
  g_clear_error (&error);
  g_assert_cmpfloat (ncm_model_param_get (model, 0), ==, 0.0);

  /* The preset requires another reparametrization. */
  ncm_mset_preset_set_reparam (preset, "NcmModelRosenbrock", "NcmReparamOther");
  ncm_mset_preset_apply (preset, mset, &error);
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH);
  g_clear_error (&error);
  g_assert_cmpfloat (ncm_model_param_get (model, 0), ==, 0.0);

  /* The preset requires the reparametrization in force. */
  ncm_mset_preset_set_reparam (preset, "NcmModelRosenbrock", "NcmReparamLinear");
  ncm_mset_preset_apply (preset, mset, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (ncm_model_param_get (model, 0), ==, 2.0);

  /* The model without reparametrization fails the same requirement. */
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH);

  ncm_matrix_free (T);
  ncm_vector_free (v);
  ncm_reparam_free (NCM_REPARAM (relin));
  NCM_TEST_FREE (ncm_mset_free, mset);
  NCM_TEST_FREE (ncm_model_free, model);
  NCM_TEST_FREE (ncm_mset_preset_free, preset);
}

static void
test_ncm_mset_preset_error_param_not_found (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NcmMSetPreset *preset;

  /* A parameter the model does not have, after a valid entry. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_value (preset, TEST_X1, 2.0);
  ncm_mset_preset_set_value (preset, "NcmModelRosenbrock:x3", 2.0);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND);
  ncm_mset_preset_clear (&preset);

  /* A registered model that is not in the NcmMSet. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_scale (preset, "NcmModelFunnel:nu", 2.0);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND);
  ncm_mset_preset_clear (&preset);

  /* A reparametrization required for a model that is not in the NcmMSet. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_value (preset, TEST_X1, 2.0);
  ncm_mset_preset_set_reparam (preset, "NcmModelFunnel", "NcmReparamLinear");
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND);
  ncm_mset_preset_clear (&preset);

  g_assert_cmpstr (g_quark_to_string (NCM_MSET_PRESET_ERROR), ==, "ncm-mset-preset-error");
}

static void
test_ncm_mset_preset_error_invalid_entry (TestNcmMSetPreset *test, gconstpointer pdata)
{
  NcmMSetPreset *preset;
  NcmVarDict *vd;

  /* Lower bound not below the upper bound. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_lower_bound (preset, TEST_X1, 1.0);
  ncm_mset_preset_set_upper_bound (preset, TEST_X1, 1.0);
  ncm_mset_preset_set_value (preset, TEST_X1, 1.0);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);

  /* The current value outside the new bounds. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_lower_bound (preset, TEST_X1, 1.0);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);

  /* A non-positive scale. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_scale (preset, TEST_X1, 0.0);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);

  /* A negative absolute tolerance. */
  preset = ncm_mset_preset_new (NULL);
  ncm_mset_preset_set_abstol (preset, TEST_X1, -1.0);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);

  /*
   * Entries of the wrong type: an integer value, a double fit flag, an integer
   * reparametrization.
   */
  vd = ncm_var_dict_new ();
  ncm_var_dict_set_int (vd, TEST_X1, 1);
  preset = g_object_new (NCM_TYPE_MSET_PRESET, "value", vd, NULL);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);
  ncm_var_dict_unref (vd);

  vd = ncm_var_dict_new ();
  ncm_var_dict_set_double (vd, TEST_X2, 1.0);
  preset = g_object_new (NCM_TYPE_MSET_PRESET, "fit", vd, NULL);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);
  ncm_var_dict_unref (vd);

  vd = ncm_var_dict_new ();
  ncm_var_dict_set_int (vd, "NcmModelRosenbrock", 1);
  preset = g_object_new (NCM_TYPE_MSET_PRESET, "reparam", vd, NULL);
  _test_assert_apply_error (preset, test->mset, test->model, NCM_MSET_PRESET_ERROR_INVALID_ENTRY);
  ncm_mset_preset_clear (&preset);
  ncm_var_dict_unref (vd);
}

static void
test_ncm_mset_preset_from_name (TestNcmMSetPreset *test, gconstpointer pdata)
{
  const gchar *old_data_dir = g_getenv (NCM_CFG_DATA_DIR_ENV);
  gchar *saved_data_dir     = g_strdup (old_data_dir);
  GError *error             = NULL;
  gchar *tmp_dir            = g_dir_make_tmp ("test_ncm_mset_preset_XXXXXX", &error);
  gchar *preset_dir;
  gchar *preset_a, *preset_b, *other, *not_yaml;
  NcmMSetPreset *preset;
  GStrv names;

  g_assert_no_error (error);

  preset_dir = g_build_filename (tmp_dir, "data", NCM_MSET_PRESET_DATA_DIR, NULL);
  g_assert_cmpint (g_mkdir_with_parents (preset_dir, 0700), ==, 0);

  preset_a = g_build_filename (preset_dir, "test_a.yaml", NULL);
  preset_b = g_build_filename (preset_dir, "test_b.yaml", NULL);
  other    = g_build_filename (preset_dir, "other.yaml", NULL);
  not_yaml = g_build_filename (preset_dir, "readme.txt", NULL);

  _test_write_file (preset_a, test_yaml);
  _test_write_file (preset_b, test_yaml);
  _test_write_file (other, "NcmModelRosenbrock: {}\n");
  _test_write_file (not_yaml, "not a preset\n");

  g_setenv (NCM_CFG_DATA_DIR_ENV, tmp_dir, TRUE);

  names = ncm_mset_preset_list_names ();
  g_assert_cmpuint (g_strv_length (names), ==, 3);
  g_assert_cmpstr (names[0], ==, "other");
  g_assert_cmpstr (names[1], ==, "test_a");
  g_assert_cmpstr (names[2], ==, "test_b");
  g_strfreev (names);

  preset = ncm_mset_preset_new_from_name ("test_a", &error);
  g_assert_no_error (error);
  ncm_mset_preset_apply (preset, test->mset, &error);
  g_assert_no_error (error);
  _test_assert_test_yaml_applied (test->model);
  ncm_mset_preset_clear (&preset);

  g_assert_null (ncm_mset_preset_new_from_name ("missing", &error));
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND);
  g_clear_error (&error);

  g_assert_null (ncm_mset_preset_new_from_name ("../mset_presets/test_a", &error));
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND);
  g_clear_error (&error);

  g_assert_null (ncm_mset_preset_new_from_name ("", &error));
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND);
  g_clear_error (&error);

  g_assert_null (ncm_mset_preset_new_from_name ("other", &error));
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_FILE);
  g_assert_true (g_str_has_prefix (error->message, "ncm_mset_preset_new_from_name: "));
  g_clear_error (&error);

  g_assert_null (ncm_mset_preset_new_from_file (preset_dir, &error));
  g_assert_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND);
  g_clear_error (&error);

  if (saved_data_dir != NULL)
    g_setenv (NCM_CFG_DATA_DIR_ENV, saved_data_dir, TRUE);
  else
    g_unsetenv (NCM_CFG_DATA_DIR_ENV);

  g_assert_cmpint (g_remove (preset_a), ==, 0);
  g_assert_cmpint (g_remove (preset_b), ==, 0);
  g_assert_cmpint (g_remove (other), ==, 0);
  g_assert_cmpint (g_remove (not_yaml), ==, 0);
  g_assert_cmpint (g_rmdir (preset_dir), ==, 0);

  {
    gchar *data_dir = g_path_get_dirname (preset_dir);

    g_assert_cmpint (g_rmdir (data_dir), ==, 0);
    g_free (data_dir);
  }

  g_assert_cmpint (g_rmdir (tmp_dir), ==, 0);

  g_free (preset_a);
  g_free (preset_b);
  g_free (other);
  g_free (not_yaml);
  g_free (preset_dir);
  g_free (tmp_dir);
  g_free (saved_data_dir);
}

/*
 * Applies every distributed preset to the model it is written for: the CMB
 * parametrization of NcHICosmoDEXcdm with one massive neutrino, NcHIPrimPowerLaw and
 * NcHIReionCamb, plus the plik NcPlanckFI model for the nuisance presets. Each value
 * in the preset must then be the parameter value.
 */
static void
test_ncm_mset_preset_shipped (void)
{
  GStrv names = ncm_mset_preset_list_names ();
  guint i;

  g_assert_true (g_strv_contains ((const gchar * const *) names, "planck2018_tt_lowl_lowe"));
  g_assert_true (g_strv_contains ((const gchar * const *) names, "planck2018_ttteee_lowl_lowe_lensing_nuisance"));
  g_assert_true (g_strv_contains ((const gchar * const *) names, "desi_dr2_bao_cmb"));

  for (i = 0; names[i] != NULL; i++)
  {
    NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
    NcHIPrim *prim   = NC_HIPRIM (nc_hiprim_power_law_new ());
    NcHICosmo *cosmo = g_object_new (NC_TYPE_HICOSMO_DE_XCDM,
                                     "massnu-length", 1,
                                     "reion", reion,
                                     "prim", prim,
                                     NULL);
    NcmModel *fi  = NULL;
    GError *error = NULL;
    NcmMSet *mset;
    NcmMSetPreset *preset;
    NcmVarDict *values;
    GStrv keys;
    guint k;

    nc_hicosmo_de_cmb_params (NC_HICOSMO_DE (cosmo), &error);
    g_assert_no_error (error);

    if (g_str_has_suffix (names[i], "_nuisance"))
      fi = NCM_MODEL (g_object_new ((strstr (names[i], "ttteee") != NULL) ? NC_TYPE_PLANCK_FI_COR_TTTEEE : NC_TYPE_PLANCK_FI_COR_TT, NULL));

    mset = ncm_mset_new (NCM_MODEL (cosmo), &error, fi, NULL);
    g_assert_no_error (error);

    preset = ncm_mset_preset_new_from_name (names[i], &error);
    g_assert_no_error (error);

    ncm_mset_preset_apply (preset, mset, &error);
    g_assert_no_error (error);

    g_object_get (preset, "value", &values, NULL);
    keys = ncm_var_dict_keys (values);
    g_assert_cmpuint (g_strv_length (keys), >, 0);

    for (k = 0; keys[k] != NULL; k++)
    {
      NcmMSetPIndex *pi = ncm_mset_param_get_by_full_name (mset, keys[k], &error);
      gdouble value;

      g_assert_no_error (error);
      g_assert_nonnull (pi);
      g_assert_true (ncm_var_dict_get_double (values, keys[k], &value));
      g_assert_cmpfloat (ncm_mset_param_get (mset, pi->mid, pi->pid), ==, value);
      ncm_mset_pindex_free (pi);
    }

    g_free (keys);
    ncm_var_dict_unref (values);
    ncm_mset_preset_free (preset);
    ncm_mset_free (mset);

    if (fi != NULL)
      ncm_model_free (fi);

    nc_hicosmo_free (cosmo);
    nc_hiprim_free (prim);
    nc_hireion_free (reion);
  }

  g_strfreev (names);
}


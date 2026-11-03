/***************************************************************************
 *            test_ncm_model_builder.c
 *
 *  Mon Sep 28 16:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

void test_ncm_model_builder_create (void);
void test_ncm_model_builder_parent_params (void);
void test_ncm_model_builder_params_arrays (void);
void test_ncm_model_builder_recreate (void);
void test_ncm_model_builder_serialize (void);
void test_ncm_model_builder_serialize_chain (void);
void test_ncm_model_builder_peek_from_type (void);
void test_ncm_model_builder_create_errors (void);
void test_ncm_model_builder_errors (void);
void test_ncm_model_builder_clash_null_error_subprocess (void);
void test_ncm_model_builder_add_after_create_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/model_builder/create", &test_ncm_model_builder_create);
  g_test_add_func ("/ncm/model_builder/parent_params", &test_ncm_model_builder_parent_params);
  g_test_add_func ("/ncm/model_builder/params_arrays", &test_ncm_model_builder_params_arrays);
  g_test_add_func ("/ncm/model_builder/recreate", &test_ncm_model_builder_recreate);
  g_test_add_func ("/ncm/model_builder/serialize", &test_ncm_model_builder_serialize);
  g_test_add_func ("/ncm/model_builder/serialize_chain", &test_ncm_model_builder_serialize_chain);
  g_test_add_func ("/ncm/model_builder/peek_from_type", &test_ncm_model_builder_peek_from_type);
  g_test_add_func ("/ncm/model_builder/create_errors", &test_ncm_model_builder_create_errors);
  g_test_add_func ("/ncm/model_builder/errors", &test_ncm_model_builder_errors);
  g_test_add_func ("/ncm/model_builder/errors/clash_null_error/subprocess", &test_ncm_model_builder_clash_null_error_subprocess);
  g_test_add_func ("/ncm/model_builder/errors/add_after_create/subprocess", &test_ncm_model_builder_add_after_create_subprocess);

  g_test_run ();
}

/*
 * A builder with one scalar and one vector parameter; @z_default changes the
 * scalar parameter's default value and @extra adds a second scalar parameter.
 */
static NcmModelBuilder *
test_builder_new (GType ptype, const gchar *name, const gdouble z_default, const gboolean extra)
{
  NcmModelBuilder *mb = ncm_model_builder_new (ptype, name, "Test model");

  ncm_model_builder_add_sparam (mb, "z", "z", -1.0, 1.0, 0.1, 0.0, z_default, NCM_PARAM_TYPE_FIXED);
  ncm_model_builder_add_vparam (mb, 2, "v", "v", -1.0, 1.0, 0.1, 0.0, 0.25, NCM_PARAM_TYPE_FIXED);

  if (extra)
    ncm_model_builder_add_sparam (mb, "w", "w", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FIXED);

  return mb;
}

void
test_ncm_model_builder_create (void)
{
  NcmModelBuilder *mb = test_builder_new (NCM_TYPE_MODEL, "TestMBCreate", 0.5, FALSE);
  NcmModelBuilder *mb_ref;
  GType type = ncm_model_builder_create (mb, NULL);
  NcmModel *model;

  g_assert_true (g_type_is_a (type, NCM_TYPE_MODEL));
  g_assert_cmpstr (g_type_name (type), ==, "TestMBCreate");
  g_assert_true (ncm_model_builder_create (mb, NULL) == type);

  model = g_object_new (type, NULL);

  g_assert_cmpuint (ncm_model_len (model), ==, 3);
  g_assert_cmpstr (ncm_model_param_name (model, 0), ==, "z");
  g_assert_cmpstr (ncm_model_param_name (model, 1), ==, "v_0");
  g_assert_cmpstr (ncm_model_param_name (model, 2), ==, "v_1");
  g_assert_cmpfloat (ncm_model_param_get (model, 0), ==, 0.5);
  g_assert_cmpfloat (ncm_model_param_get (model, 1), ==, 0.25);

  /* Registered under its own name, as the parent NcmModel is not registered. */
  g_assert_cmpint (ncm_model_id (model), >=, 0);
  g_assert_cmpstr (ncm_mset_get_ns_by_id (ncm_model_id (model)), ==, "TestMBCreate");

  ncm_model_free (model);

  mb_ref = ncm_model_builder_ref (mb);
  g_assert_true (mb_ref == mb);
  ncm_model_builder_free (mb_ref);

  ncm_model_builder_clear (&mb);
  g_assert_null (mb);
  ncm_model_builder_clear (&mb);
}

void
test_ncm_model_builder_parent_params (void)
{
  NcmModelBuilder *mb = test_builder_new (NCM_TYPE_MODEL_ROSENBROCK, "TestMBRosenbrockExt", 0.5, FALSE);
  GType type          = ncm_model_builder_create (mb, NULL);
  NcmModel *model     = g_object_new (type, NULL);
  NcmModel *mrb       = NCM_MODEL (ncm_model_rosenbrock_new ());

  g_assert_true (g_type_is_a (type, NCM_TYPE_MODEL_ROSENBROCK));

  /* The new parameters come after the parent's. */
  g_assert_cmpuint (ncm_model_len (model), ==, 5);
  g_assert_cmpuint (ncm_model_sparam_len (model), ==, 3);
  g_assert_cmpuint (ncm_model_vparam_array_len (model), ==, 1);
  g_assert_cmpstr (ncm_model_param_name (model, 0), ==, ncm_model_param_name (mrb, 0));
  g_assert_cmpstr (ncm_model_param_name (model, 1), ==, ncm_model_param_name (mrb, 1));
  g_assert_cmpstr (ncm_model_param_name (model, 2), ==, "z");
  g_assert_cmpstr (ncm_model_param_name (model, 3), ==, "v_0");
  g_assert_cmpfloat (ncm_model_param_get (model, 0), ==, ncm_model_param_get (mrb, 0));
  g_assert_cmpfloat (ncm_model_param_get (model, 2), ==, 0.5);

  /* It shares the parent's model id, as C subclasses do. */
  g_assert_cmpint (ncm_model_id (model), ==, ncm_model_id (mrb));

  ncm_model_free (mrb);
  ncm_model_free (model);
  ncm_model_builder_free (mb);
}

void
test_ncm_model_builder_params_arrays (void)
{
  NcmModelBuilder *mb  = ncm_model_builder_new (NCM_TYPE_MODEL, "TestMBArrays", "Test model");
  NcmObjArray *sparams = ncm_obj_array_new ();
  NcmObjArray *vparams = ncm_obj_array_new ();
  NcmSParam *a         = ncm_sparam_new ("a", "a", -1.0, 1.0, 0.1, 0.0, 0.1, NCM_PARAM_TYPE_FIXED);
  NcmSParam *b         = ncm_sparam_new ("b", "b", -1.0, 1.0, 0.1, 0.0, 0.2, NCM_PARAM_TYPE_FREE);
  NcmVParam *u         = ncm_vparam_full_new (3, "u", "u", -1.0, 1.0, 0.1, 0.0, 0.3, NCM_PARAM_TYPE_FIXED);
  NcmObjArray *got;

  ncm_obj_array_add (sparams, G_OBJECT (a));
  ncm_obj_array_add (sparams, G_OBJECT (b));
  ncm_obj_array_add (vparams, G_OBJECT (u));

  ncm_model_builder_add_sparams (mb, sparams);
  ncm_model_builder_add_vparams (mb, vparams);

  got = ncm_model_builder_get_sparams (mb);
  g_assert_cmpuint (got->len, ==, 2);
  g_assert_true (ncm_obj_array_peek (got, 0) == G_OBJECT (a));
  g_assert_true (ncm_obj_array_peek (got, 1) == G_OBJECT (b));
  ncm_obj_array_unref (got);

  got = ncm_model_builder_get_vparams (mb);
  g_assert_cmpuint (got->len, ==, 1);
  g_assert_true (ncm_obj_array_peek (got, 0) == G_OBJECT (u));
  ncm_obj_array_unref (got);

  {
    NcmModel *model = g_object_new (ncm_model_builder_create (mb, NULL), NULL);

    g_assert_cmpuint (ncm_model_len (model), ==, 5);
    g_assert_cmpstr (ncm_model_param_name (model, 1), ==, "b");
    g_assert_cmpstr (ncm_model_param_name (model, 4), ==, "u_2");

    ncm_model_free (model);
  }

  ncm_sparam_free (a);
  ncm_sparam_free (b);
  ncm_vparam_free (u);
  ncm_obj_array_unref (sparams);
  ncm_obj_array_unref (vparams);
  ncm_model_builder_free (mb);
}

void
test_ncm_model_builder_recreate (void)
{
  NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBRecreate", 0.5, FALSE);
  NcmModelBuilder *mb2 = test_builder_new (NCM_TYPE_MODEL, "TestMBRecreate", 0.5, FALSE);

  /* A second builder with the same definition gets the existing type. */
  g_assert_true (ncm_model_builder_create (mb1, NULL) == ncm_model_builder_create (mb2, NULL));

  ncm_model_builder_free (mb1);
  ncm_model_builder_free (mb2);
}

void
test_ncm_model_builder_serialize (void)
{
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  NcmModelBuilder *mb = test_builder_new (NCM_TYPE_MODEL_ROSENBROCK, "TestMBSerialize", 0.5, FALSE);
  GType type          = ncm_model_builder_create (mb, NULL);
  gchar *mb_ser       = ncm_serialize_to_string (ser, G_OBJECT (mb), TRUE);
  NcmModelBuilder *mb_dup;
  gchar *parent_type_string;

  ncm_serialize_reset (ser, TRUE);
  mb_dup = NCM_MODEL_BUILDER (ncm_serialize_from_string (ser, mb_ser));

  g_object_get (mb_dup, "parent-type-string", &parent_type_string, NULL);
  g_assert_cmpstr (parent_type_string, ==, "NcmModelRosenbrock");
  g_free (parent_type_string);

  /* The rebuilt builder describes the same type. */
  g_assert_true (ncm_model_builder_create (mb_dup, NULL) == type);

  g_free (mb_ser);
  ncm_model_builder_free (mb_dup);
  ncm_model_builder_free (mb);
  ncm_serialize_free (ser);
}

/*
 * A builder whose parent is a builder type deserializes before that parent is
 * created; the types are then created in array order.
 */
void
test_ncm_model_builder_serialize_chain (void)
{
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmModelBuilder *mb_parent = test_builder_new (NCM_TYPE_MODEL, "TestMBChainParent", 0.5, FALSE);
  NcmModelBuilder *mb_child  = g_object_new (NCM_TYPE_MODEL_BUILDER,
                                             "parent-type-string", "TestMBChainParent",
                                             "name", "TestMBChainChild",
                                             "description", "Test model",
                                             NULL);
  NcmObjArray *oa = ncm_obj_array_new ();
  NcmObjArray *oa_dup;
  GVariant *var;
  GType ptype, ctype;
  NcmModel *model;

  ncm_model_builder_add_sparam (mb_child, "w", "w", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FIXED);
  ncm_obj_array_add (oa, G_OBJECT (mb_parent));
  ncm_obj_array_add (oa, G_OBJECT (mb_child));

  g_assert_true (g_type_from_name ("TestMBChainParent") == 0);

  var = ncm_serialize_array_to_variant (ser, oa);
  ncm_serialize_reset (ser, TRUE);
  oa_dup = ncm_serialize_array_from_variant (ser, var);

  ptype = ncm_model_builder_create (NCM_MODEL_BUILDER (ncm_obj_array_peek (oa_dup, 0)), NULL);
  ctype = ncm_model_builder_create (NCM_MODEL_BUILDER (ncm_obj_array_peek (oa_dup, 1)), NULL);

  g_assert_true (g_type_parent (ctype) == ptype);
  g_assert_true (ncm_model_builder_peek_from_type (ptype) == NCM_MODEL_BUILDER (ncm_obj_array_peek (oa_dup, 0)));
  g_assert_true (ncm_model_builder_peek_from_type (ctype) == NCM_MODEL_BUILDER (ncm_obj_array_peek (oa_dup, 1)));

  model = g_object_new (ctype, NULL);
  g_assert_cmpuint (ncm_model_len (model), ==, 4);
  g_assert_cmpstr (ncm_model_param_name (model, 1), ==, "w");

  ncm_model_free (model);
  g_variant_unref (var);
  ncm_obj_array_unref (oa_dup);
  ncm_obj_array_unref (oa);
  ncm_model_builder_free (mb_child);
  ncm_model_builder_free (mb_parent);
  ncm_serialize_free (ser);
}

void
test_ncm_model_builder_peek_from_type (void)
{
  NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBPeek", 0.5, FALSE);
  NcmModelBuilder *mb2 = test_builder_new (NCM_TYPE_MODEL, "TestMBPeek", 0.5, FALSE);
  GType type           = ncm_model_builder_create (mb1, NULL);
  gchar *name;

  g_assert_true (ncm_model_builder_peek_from_type (type) == mb1);

  /* The builder that registered the type stays, not a later equal one. */
  g_assert_true (ncm_model_builder_create (mb2, NULL) == type);
  g_assert_true (ncm_model_builder_peek_from_type (type) == mb1);

  g_assert_null (ncm_model_builder_peek_from_type (NCM_TYPE_MODEL));
  g_assert_null (ncm_model_builder_peek_from_type (NCM_TYPE_MODEL_ROSENBROCK));
  g_assert_null (ncm_model_builder_peek_from_type (NCM_TYPE_VECTOR));

  ncm_model_builder_free (mb2);
  ncm_model_builder_free (mb1);

  /* The type holds its builder after the caller's references are gone. */
  g_object_get (ncm_model_builder_peek_from_type (type), "name", &name, NULL);
  g_assert_cmpstr (name, ==, "TestMBPeek");
  g_free (name);
}

void
test_ncm_model_builder_create_errors (void)
{
  GError *error = NULL;

  {
    NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashParam", 0.5, FALSE);
    NcmModelBuilder *mb2 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashParam", 0.7, FALSE);

    ncm_model_builder_create (mb1, NULL);
    g_assert_true (ncm_model_builder_create (mb2, &error) == G_TYPE_INVALID);
    g_assert_error (error, NCM_MODEL_BUILDER_ERROR, NCM_MODEL_BUILDER_ERROR_NAME_CLASH);
    g_assert_nonnull (strstr (error->message, "scalar parameter 0 (`z' there, `z' here) differs"));
    g_clear_error (&error);

    ncm_model_builder_free (mb1);
    ncm_model_builder_free (mb2);
  }

  {
    NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashLen", 0.5, FALSE);
    NcmModelBuilder *mb2 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashLen", 0.5, TRUE);

    ncm_model_builder_create (mb1, NULL);
    g_assert_true (ncm_model_builder_create (mb2, &error) == G_TYPE_INVALID);
    g_assert_error (error, NCM_MODEL_BUILDER_ERROR, NCM_MODEL_BUILDER_ERROR_NAME_CLASH);
    g_assert_nonnull (strstr (error->message, "with 1 scalar and 1 vector parameter(s), but this builder has 2 and 1"));
    g_clear_error (&error);

    ncm_model_builder_free (mb1);
    ncm_model_builder_free (mb2);
  }

  {
    NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashParent", 0.5, FALSE);
    NcmModelBuilder *mb2 = test_builder_new (NCM_TYPE_MODEL_ROSENBROCK, "TestMBClashParent", 0.5, FALSE);

    ncm_model_builder_create (mb1, NULL);
    g_assert_true (ncm_model_builder_create (mb2, &error) == G_TYPE_INVALID);
    g_assert_error (error, NCM_MODEL_BUILDER_ERROR, NCM_MODEL_BUILDER_ERROR_NAME_CLASH);
    g_assert_nonnull (strstr (error->message, "with parent `NcmModel', but this builder has parent `NcmModelRosenbrock'"));
    g_clear_error (&error);

    ncm_model_builder_free (mb1);
    ncm_model_builder_free (mb2);
  }

  {
    NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashStack", 0.5, FALSE);
    NcmModelBuilder *mb2 = g_object_new (NCM_TYPE_MODEL_BUILDER,
                                         "parent-type-string", "NcmModel",
                                         "name", "TestMBClashStack",
                                         "description", "Test model",
                                         "stackable", TRUE,
                                         NULL);

    ncm_model_builder_add_sparam (mb2, "z", "z", -1.0, 1.0, 0.1, 0.0, 0.5, NCM_PARAM_TYPE_FIXED);
    ncm_model_builder_add_vparam (mb2, 2, "v", "v", -1.0, 1.0, 0.1, 0.0, 0.25, NCM_PARAM_TYPE_FIXED);

    ncm_model_builder_create (mb1, NULL);
    g_assert_true (ncm_model_builder_create (mb2, &error) == G_TYPE_INVALID);
    g_assert_error (error, NCM_MODEL_BUILDER_ERROR, NCM_MODEL_BUILDER_ERROR_NAME_CLASH);
    g_assert_nonnull (strstr (error->message, "with stackable FALSE, but this builder has stackable TRUE"));
    g_clear_error (&error);

    ncm_model_builder_free (mb1);
    ncm_model_builder_free (mb2);
  }

  {
    NcmModelBuilder *mb = g_object_new (NCM_TYPE_MODEL_BUILDER, "parent-type-string", "NoSuchModel", "name", "TestMBNoParent", NULL);

    g_assert_true (ncm_model_builder_create (mb, &error) == G_TYPE_INVALID);
    g_assert_error (error, NCM_MODEL_BUILDER_ERROR, NCM_MODEL_BUILDER_ERROR_PARENT_NOT_REGISTERED);
    g_assert_nonnull (strstr (error->message, "parent type `NoSuchModel' of `TestMBNoParent' is not registered"));
    g_clear_error (&error);

    ncm_model_builder_free (mb);
  }

  {
    NcmModelBuilder *mb = g_object_new (NCM_TYPE_MODEL_BUILDER, "parent-type-string", "NcmVector", "name", "TestMBVectorParent", NULL);

    g_assert_true (ncm_model_builder_create (mb, &error) == G_TYPE_INVALID);
    g_assert_error (error, NCM_MODEL_BUILDER_ERROR, NCM_MODEL_BUILDER_ERROR_PARENT_NOT_MODEL);
    g_assert_nonnull (strstr (error->message, "parent type `NcmVector' of `TestMBVectorParent' is not a NcmModel"));
    g_clear_error (&error);

    ncm_model_builder_free (mb);
  }
}

void
test_ncm_model_builder_errors (void)
{
  g_test_trap_subprocess ("/ncm/model_builder/errors/clash_null_error/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*a type named `TestMBClashNull' already exists, but its scalar parameter 0 "
                             "(`z' there, `z' here) differs in name, symbol, bounds, scale, tolerance, default value "
                             "or fit type*");

  g_test_trap_subprocess ("/ncm/model_builder/errors/add_after_create/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*model `TestMBAddAfter' was already created, cannot add the parameter `w'*");
}

/* With a NULL error location the clash aborts. */
void
test_ncm_model_builder_clash_null_error_subprocess (void)
{
  NcmModelBuilder *mb1 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashNull", 0.5, FALSE);
  NcmModelBuilder *mb2 = test_builder_new (NCM_TYPE_MODEL, "TestMBClashNull", 0.7, FALSE);

  ncm_model_builder_create (mb1, NULL);
  ncm_model_builder_create (mb2, NULL);

  ncm_model_builder_free (mb1);
  ncm_model_builder_free (mb2);
}

void
test_ncm_model_builder_add_after_create_subprocess (void)
{
  NcmModelBuilder *mb = test_builder_new (NCM_TYPE_MODEL, "TestMBAddAfter", 0.5, FALSE);

  ncm_model_builder_create (mb, NULL);
  ncm_model_builder_add_sparam (mb, "w", "w", -1.0, 1.0, 0.1, 0.0, 0.0, NCM_PARAM_TYPE_FIXED);

  ncm_model_builder_free (mb);
}


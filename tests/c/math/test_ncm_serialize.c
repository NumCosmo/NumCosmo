/***************************************************************************
 *            test_ncm_object_serialization.c
 *
 *  Wed June 20 21:56:09 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br> <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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


G_DECLARE_FINAL_TYPE (NcmObjectTest, ncm_object_test, NCM, OBJECT_TEST, GObject)

struct _NcmObjectTest
{
  GObject parent_instance;

  NcmDTuple2 *dt2;
  NcmDTuple3 *dt3;
  NcmVector *vector;
  NcmMatrix *matrix;
};

enum
{
  PROP_0,
  PROP_DT2,
  PROP_DT3,
  PROP_VECTOR,
  PROP_MATRIX,
  PROP_LEN,
};

G_DEFINE_TYPE (NcmObjectTest, ncm_object_test, G_TYPE_OBJECT)

static void
ncm_object_test_init (NcmObjectTest *self)
{
  self->dt2 = NULL;
  self->dt3 = NULL;
}

static void
ncm_object_test_finalize (GObject *object)
{
  NcmObjectTest *self = NCM_OBJECT_TEST (object);

  ncm_dtuple2_clear (&self->dt2);
  ncm_dtuple3_clear (&self->dt3);
  ncm_vector_clear (&self->vector);
  ncm_matrix_clear (&self->matrix);

  /* Chain up to the parent class */
  G_OBJECT_CLASS (ncm_object_test_parent_class)->finalize (object);
}

static void
ncm_object_test_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec   *pspec)
{
  NcmObjectTest *self = NCM_OBJECT_TEST (object);

  switch (prop_id)
  {
    case PROP_DT2:
      ncm_dtuple2_clear (&self->dt2);
      self->dt2 = g_value_dup_boxed (value);
      break;
    case PROP_DT3:
      ncm_dtuple3_clear (&self->dt3);
      self->dt3 = g_value_dup_boxed (value);
      break;
    case PROP_VECTOR:
      ncm_vector_clear (&self->vector);
      self->vector = g_value_dup_object (value);
      break;
    case PROP_MATRIX:
      ncm_matrix_clear (&self->matrix);
      self->matrix = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_object_test_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmObjectTest *self = NCM_OBJECT_TEST (object);

  switch (prop_id)
  {
    case PROP_DT2:
      g_value_set_boxed (value, self->dt2);
      break;
    case PROP_DT3:
      g_value_set_boxed (value, self->dt3);
      break;
    case PROP_VECTOR:
      g_value_set_object (value, self->vector);
      break;
    case PROP_MATRIX:
      g_value_set_object (value, self->matrix);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_object_test_class_init (NcmObjectTestClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize     = ncm_object_test_finalize;
  object_class->set_property = ncm_object_test_set_property;
  object_class->get_property = ncm_object_test_get_property;

  g_object_class_install_property (object_class,
                                   PROP_DT2,
                                   g_param_spec_boxed ("dt2",
                                                       NULL,
                                                       "dt2",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DT3,
                                   g_param_spec_boxed ("dt3",
                                                       NULL,
                                                       "dt3",
                                                       NCM_TYPE_DTUPLE3,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_VECTOR,
                                   g_param_spec_object ("vector",
                                                        NULL,
                                                        "vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MATRIX,
                                   g_param_spec_object ("matrix",
                                                        NULL,
                                                        "matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

typedef struct _TestNcmObjectTest
{
  NcmObjectTest *obj;
} TestNcmObjectTest;

typedef struct _TestNcmSerialize
{
  NcmSerialize *ser;
} TestNcmSerialize;

static void test_ncm_serialize_new (TestNcmSerialize *test, gconstpointer pdata);
void test_ncm_serialize_new_noclean_dup (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_free (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_global_from_string_plain (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_from_string_params (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_from_string_nest (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_from_string_plain_named (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_string_params_named (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_string_nest_named (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_string_nest_samename (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_to_file_from_file (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_to_binfile_from_binfile (TestNcmSerialize *test, gconstpointer pdata);

#ifdef HAVE_LIBFYAML
static void test_ncm_serialize_to_yaml_from_yaml (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_yaml_special_types (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_yaml_special_types_block_flow (TestNcmSerialize *test, gconstpointer pdata);

#endif /* HAVE_LIBFYAML */

static void test_ncm_serialize_reset_autosave_only (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_traps (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_invalid_from_string_syntax (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_invalid_from_string_nonexist (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_invalid_from_string_mnames (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_invalid_from_string_samename (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_invalid_from_string_wrongorder (TestNcmSerialize *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  ncm_cfg_register_obj (ncm_object_test_get_type ());

  g_test_add ("/ncm/serialize/global/from_string/plain", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_from_string_plain,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/global/from_string/params", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_from_string_params,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/global/from_string/nest", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_from_string_nest,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/plain/named", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_plain_named,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/params/named", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_params_named,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/nested/named", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_nest_named,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/nested/samename", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_nest_samename,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/plain/reset/autosave_only", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_reset_autosave_only,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/to_file/from_file", TestNcmSerialize, NULL,
              &test_ncm_serialize_new_noclean_dup,
              &test_ncm_serialize_to_file_from_file,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/to_binfile/from_binfile", TestNcmSerialize, NULL,
              &test_ncm_serialize_new_noclean_dup,
              &test_ncm_serialize_to_binfile_from_binfile,
              &test_ncm_serialize_free);
#ifdef HAVE_LIBFYAML
  g_test_add ("/ncm/serialize/to_yaml/from_yaml", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_to_yaml_from_yaml,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_yaml/special_types", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_yaml_special_types,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_yaml/special_types/block_flow", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_yaml_special_types_block_flow,
              &test_ncm_serialize_free);

#endif /* HAVE_LIBFYAML */

  g_test_add ("/ncm/serialize/traps", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_traps,
              &test_ncm_serialize_free);

  g_test_add ("/ncm/serialize/invalid/from_string/syntax/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_invalid_from_string_syntax,
              &test_ncm_serialize_free);

  g_test_add ("/ncm/serialize/invalid/from_string/nonexistent/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_invalid_from_string_nonexist,
              &test_ncm_serialize_free);

  g_test_add ("/ncm/serialize/invalid/from_string/multiplenames/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_invalid_from_string_mnames,
              &test_ncm_serialize_free);

  g_test_add ("/ncm/serialize/invalid/from_string/samename/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_invalid_from_string_samename,
              &test_ncm_serialize_free);

  g_test_add ("/ncm/serialize/invalid/from_string/wrongorder/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_invalid_from_string_wrongorder,
              &test_ncm_serialize_free);

  g_test_run ();
}

void
test_ncm_serialize_new (TestNcmSerialize *test, gconstpointer pdata)
{
  test->ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
}

void
test_ncm_serialize_new_noclean_dup (TestNcmSerialize *test, gconstpointer pdata)
{
  test->ser = ncm_serialize_new (0);
}

void
test_ncm_serialize_free (TestNcmSerialize *test, gconstpointer pdata)
{
  NcmSerialize *ser = test->ser;

  NCM_TEST_FREE (ncm_serialize_free, ser);
}

void
test_ncm_serialize_global_from_string_plain (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj       = ncm_serialize_global_from_string ("NcHICosmoLCDM");
  NcHICosmo *hic     = NC_HICOSMO (obj);
  gchar *obj_ser     = ncm_serialize_global_to_string (obj, TRUE);
  GObject *obj_new   = ncm_serialize_global_from_string (obj_ser);
  gchar *obj_new_ser = ncm_serialize_global_to_string (obj_new, TRUE);

  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);

  obj_ser     = ncm_serialize_global_to_string (obj, FALSE);
  obj_new_ser = ncm_serialize_global_to_string (obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  NCM_TEST_FREE (nc_hicosmo_free, hic);
}

void
test_ncm_serialize_global_from_string_params (TestNcmSerialize *test, gconstpointer pdata)
{
  NcmModel *m = NCM_MODEL (ncm_serialize_global_from_string ("NcHICosmoLCDM{'H0':<12.3>,'Omegac':<0.2>}"));

  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);

  NCM_TEST_FREE (ncm_model_free, m);

  m = NCM_MODEL (ncm_serialize_global_from_string ("('NcHICosmoLCDM', {'H0':<12.3>,'Omegac':<0.2>})"));
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);

  NCM_TEST_FREE (ncm_model_free, m);
}

void
test_ncm_serialize_global_from_string_nest (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_global_from_string ("NcPowspecMLTransfer{'transfer':<('NcTransferFuncEH',@a{sv} {})>}");

  NcPowspecML *ps_ml = NC_POWSPEC_ML (obj);

  gchar *obj_ser     = ncm_serialize_global_to_string (obj, TRUE);
  GObject *obj_new   = ncm_serialize_global_from_string (obj_ser);
  gchar *obj_new_ser = ncm_serialize_global_to_string (obj_new, TRUE);

  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);

  g_free (obj_ser);
  g_free (obj_new_ser);

  obj_ser     = ncm_serialize_global_to_string (obj, FALSE);
  obj_new_ser = ncm_serialize_global_to_string (obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  g_assert_true (NC_IS_TRANSFER_FUNC_EH (nc_powspec_ml_transfer_peek_tf (NC_POWSPEC_ML_TRANSFER (ps_ml))));

  NCM_TEST_FREE (nc_powspec_ml_free, ps_ml);
}

static void
test_ncm_serialize_from_string_plain_named (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj       = ncm_serialize_from_string (test->ser, "NcHICosmoLCDM[myname:is]");
  NcHICosmo *hic     = NC_HICOSMO (obj);
  gchar *obj_ser     = ncm_serialize_to_string (test->ser, obj, TRUE);
  GObject *obj_new   = ncm_serialize_from_string (test->ser, obj_ser);
  gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);

  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);

  obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);

  g_assert_cmpstr (obj_ser, ==, obj_new_ser);

  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  ncm_serialize_clear_instances (test->ser, FALSE);

  NCM_TEST_FREE (nc_hicosmo_free, hic);
}

static void
test_ncm_serialize_reset_autosave_only (TestNcmSerialize *test, gconstpointer pdata)
{
  /* Deserialized from string, since it is named it will be automatically added to the saved instances. */
  GObject *obj1 = ncm_serialize_from_string (test->ser, "NcHICosmoLCDM[myname:is]");
  /* Deserialized from string, since it is *not* named it will *not* be automatically added to the saved instances.*/
  GObject *obj2 = ncm_serialize_from_string (test->ser, "NcPowspecMLTransfer{'transfer':<('NcTransferFuncEH',@a{sv} {})>}");

  /* Serializing from object */
  gchar *obj_ser1 = ncm_serialize_to_string (test->ser, obj1, FALSE);
  /* Serializing from object, all objects will receive a name */
  gchar *obj_ser2 = ncm_serialize_to_string (test->ser, obj2, FALSE);

  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), ==, 1);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), >=, 2);

  {
    GObject *obj1_ref2 = ncm_serialize_from_string (test->ser, obj_ser1);
    GObject *obj2_ref2 = ncm_serialize_from_string (test->ser, obj_ser2);

    g_assert_true (obj1_ref2 == obj1);
    g_assert_true (obj2_ref2 != obj2);

    g_object_unref (obj1_ref2);
    g_object_unref (obj2_ref2);

    g_free (obj_ser1);
    g_free (obj_ser2);
  }

  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), >, 1);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), >=, 2);

  ncm_serialize_reset (test->ser, TRUE);

  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), ==, 1);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), ==, 0);

  ncm_serialize_reset (test->ser, FALSE);

  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), ==, 0);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), ==, 0);

  NCM_TEST_FREE (g_object_unref, obj1);
  NCM_TEST_FREE (g_object_unref, obj2);
}

void
test_ncm_serialize_from_string_params_named (TestNcmSerialize *test, gconstpointer pdata)
{
  NcmModel *m = NCM_MODEL (ncm_serialize_from_string (test->ser, "NcHICosmoLCDM[T0]{'H0':<12.3>,'Omegac':<0.2>}"));

  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);

  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (ncm_model_free, m);

  m = NCM_MODEL (ncm_serialize_from_string (test->ser, "('NcHICosmoLCDM[T1]', {'H0':<12.3>,'Omegac':<0.2>})"));
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);

  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (ncm_model_free, m);
}

void
test_ncm_serialize_from_string_nest_named (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser, "NcPowspecMLTransfer[T0]{'transfer':<('NcTransferFuncEH[T1]',@a{sv} {})>}");

  NcPowspecML *ps_ml = NC_POWSPEC_ML (obj);

  gchar *obj_ser     = ncm_serialize_to_string (test->ser, obj, TRUE);
  GObject *obj_new   = ncm_serialize_from_string (test->ser, obj_ser);
  gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);

  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);

  g_free (obj_ser);
  g_free (obj_new_ser);

  obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  g_assert_true (NC_IS_TRANSFER_FUNC_EH (nc_powspec_ml_transfer_peek_tf (NC_POWSPEC_ML_TRANSFER (ps_ml))));

  g_assert_cmphex (GPOINTER_TO_SIZE (ncm_serialize_peek_by_name (test->ser, "T0")), ==, GPOINTER_TO_SIZE (ps_ml));
  g_assert_cmphex (GPOINTER_TO_SIZE (ncm_serialize_peek_by_name (test->ser, "T1")), ==, GPOINTER_TO_SIZE (nc_powspec_ml_transfer_peek_tf (NC_POWSPEC_ML_TRANSFER (ps_ml))));

  g_assert_cmpstr (ncm_serialize_peek_name (test->ser, nc_powspec_ml_transfer_peek_tf (NC_POWSPEC_ML_TRANSFER (ps_ml))), ==, "T1");
  g_assert_cmpstr (ncm_serialize_peek_name (test->ser, ps_ml), ==, "T0");
  g_variant_unref (ncm_serialize_to_variant (test->ser, obj));

  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (nc_powspec_ml_free, ps_ml);
}

void
test_ncm_serialize_from_string_nest_samename (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser,
                                            "NcmSplineCubicNotaknot{'length':<6>, "
                                            "'x' : <('NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>})>, "
                                            "'y' : <('NcmVector[T0]', @a{sv} {})>}");

  NcmSpline *s = NCM_SPLINE (obj);
  gchar *obj_ser, *obj_new_ser;
  GObject *obj_new;

  ncm_serialize_clear_instances (test->ser, FALSE);

  obj_ser = ncm_serialize_to_string (test->ser, obj, TRUE);
  obj_new = ncm_serialize_from_string (test->ser, obj_ser);

  ncm_serialize_reset (test->ser, FALSE);

  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);

  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);

  g_free (obj_ser);
  g_free (obj_new_ser);

  ncm_serialize_reset (test->ser, FALSE);
  obj_ser = ncm_serialize_to_string (test->ser, obj, FALSE);

  ncm_serialize_reset (test->ser, FALSE);
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);

  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (ncm_spline_free, s);
}

static void
test_ncm_serialize_to_file_from_file (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj    = ncm_serialize_from_string (test->ser, "NcHICosmoDEXcdm{'w':<-2.0>}");
  NcHICosmo *hic  = NC_HICOSMO (obj);
  gchar *obj_ser  = ncm_serialize_to_string (test->ser, obj, TRUE);
  gchar *tmp_dir  = g_dir_make_tmp ("test_serialize_file_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/test_serialize_file.obj", tmp_dir);

  ncm_serialize_to_file (test->ser, obj, filename);

  {
    GObject *obj_new   = ncm_serialize_from_file (test->ser, filename);
    gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);

    g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
    g_assert_cmpstr (obj_ser, ==, obj_new_ser);
    g_free (obj_ser);
    g_free (obj_new_ser);

    obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
    obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);

    g_assert_cmpstr (obj_ser, ==, obj_new_ser);

    g_free (obj_ser);
    g_free (obj_new_ser);
    g_object_unref (obj_new);

    ncm_serialize_clear_instances (test->ser, FALSE);

    NCM_TEST_FREE (nc_hicosmo_free, hic);
  }

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);
}

static void
test_ncm_serialize_to_binfile_from_binfile (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj    = ncm_serialize_from_string (test->ser, "NcHICosmoDEXcdm{'w':<-2.0>}");
  NcHICosmo *hic  = NC_HICOSMO (obj);
  gchar *obj_ser  = ncm_serialize_to_string (test->ser, obj, TRUE);
  gchar *tmp_dir  = g_dir_make_tmp ("test_serialize_binfile_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/test_serialize_binfile.obj", tmp_dir);

  ncm_serialize_to_binfile (test->ser, obj, filename);

  {
    GObject *obj_new   = ncm_serialize_from_binfile (test->ser, filename);
    gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);

    g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
    g_assert_cmpstr (obj_ser, ==, obj_new_ser);
    g_free (obj_ser);
    g_free (obj_new_ser);

    obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
    obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);

    g_assert_cmpstr (obj_ser, ==, obj_new_ser);

    g_free (obj_ser);
    g_free (obj_new_ser);
    g_object_unref (obj_new);

    ncm_serialize_clear_instances (test->ser, FALSE);

    NCM_TEST_FREE (nc_hicosmo_free, hic);
  }

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);
}

#ifdef HAVE_LIBFYAML

static void
test_ncm_serialize_to_yaml_from_yaml (TestNcmSerialize *test, gconstpointer pdata)
{
  NcDistance *dist = nc_distance_new (3.0);
  NcmData *bao1    = nc_data_bao_create (dist, NC_DATA_BAO_A_EISENSTEIN2005);
  NcmData *bao2    = nc_data_bao_create (dist, NC_DATA_BAO_RDV_BOSS_QSO_ATA2017);
  GObject *obj     = G_OBJECT (ncm_dataset_new_list (bao1, bao2, NULL));
  gchar *yaml_str  = ncm_serialize_to_yaml (test->ser, obj);
  GObject *obj_new = ncm_serialize_from_yaml (test->ser, yaml_str);

  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));

  ncm_serialize_reset (test->ser, FALSE);
  {
    gchar *yaml_str_new = ncm_serialize_to_yaml (test->ser, obj_new);

    g_assert_cmpstr (yaml_str, ==, yaml_str_new);
    g_free (yaml_str_new);
  }

  nc_distance_free (dist);
  ncm_data_free (bao1);
  ncm_data_free (bao2);
  g_object_unref (obj);
  g_object_unref (obj_new);
  g_free (yaml_str);
}

static void
test_ncm_serialize_from_yaml_special_types (TestNcmSerialize *test, gconstpointer pdata)
{
  const gchar *yaml_str =
    "NcmObjectTest:\n"
    "  dt2: (3.4, -3.3)\n"
    "  dt3: (0.4, 4.3,1.0e44)\n"
    "  vector: [1.0, 2.0, 3.0]\n"
    "  matrix: [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]\n";
  GObject *obj_new = ncm_serialize_from_yaml (test->ser, yaml_str);

  g_assert_true (G_IS_OBJECT (obj_new));
  g_assert_true (NCM_IS_OBJECT_TEST (obj_new));

  {
    NcmDTuple2 *dt2;
    NcmDTuple3 *dt3;
    NcmVector *vector;
    NcmMatrix *matrix;

    g_object_get (obj_new, "dt2", &dt2,
                  "dt3", &dt3,
                  "vector", &vector,
                  "matrix", &matrix,
                  NULL);

    g_assert_cmpfloat (dt2->elements[0], ==, 3.4);
    g_assert_cmpfloat (dt2->elements[1], ==, -3.3);

    g_assert_cmpfloat (dt3->elements[0], ==, 0.4);
    g_assert_cmpfloat (dt3->elements[1], ==, 4.3);
    g_assert_cmpfloat (dt3->elements[2], ==, 1.0e44);

    g_assert_cmpfloat (ncm_vector_get (vector, 0), ==, 1.0);
    g_assert_cmpfloat (ncm_vector_get (vector, 1), ==, 2.0);
    g_assert_cmpfloat (ncm_vector_get (vector, 2), ==, 3.0);

    g_assert_cmpfloat (ncm_matrix_get (matrix, 0, 0), ==, 1.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 0, 1), ==, 2.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 0, 2), ==, 3.0);

    g_assert_cmpfloat (ncm_matrix_get (matrix, 1, 0), ==, 4.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 1, 1), ==, 5.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 1, 2), ==, 6.0);

    ncm_dtuple2_free (dt2);
    ncm_dtuple3_free (dt3);

    ncm_vector_free (vector);
    ncm_matrix_free (matrix);
  }
}

static void
test_ncm_serialize_from_yaml_special_types_block_flow (TestNcmSerialize *test, gconstpointer pdata)
{
  const gchar *yaml_str =
    "NcmObjectTest:\n"
    "  dt2: (3.4, -3.3)\n"
    "  dt3: (0.4, 4.3,1.0e44)\n"
    "  vector:\n"
    "   - 1.0\n"
    "   - 2.0\n"
    "   - 3.0\n"
    "  matrix:\n"
    "   - [1.0, 2.0, 3.0]\n"
    "   - [4.0, 5.0, 6.0]\n";
  GObject *obj_new = ncm_serialize_from_yaml (test->ser, yaml_str);

  g_assert_true (G_IS_OBJECT (obj_new));
  g_assert_true (NCM_IS_OBJECT_TEST (obj_new));

  {
    NcmDTuple2 *dt2;
    NcmDTuple3 *dt3;
    NcmVector *vector;
    NcmMatrix *matrix;

    g_object_get (obj_new, "dt2", &dt2,
                  "dt3", &dt3,
                  "vector", &vector,
                  "matrix", &matrix,
                  NULL);

    g_assert_cmpfloat (dt2->elements[0], ==, 3.4);
    g_assert_cmpfloat (dt2->elements[1], ==, -3.3);

    g_assert_cmpfloat (dt3->elements[0], ==, 0.4);
    g_assert_cmpfloat (dt3->elements[1], ==, 4.3);
    g_assert_cmpfloat (dt3->elements[2], ==, 1.0e44);

    g_assert_cmpfloat (ncm_vector_get (vector, 0), ==, 1.0);
    g_assert_cmpfloat (ncm_vector_get (vector, 1), ==, 2.0);
    g_assert_cmpfloat (ncm_vector_get (vector, 2), ==, 3.0);

    g_assert_cmpfloat (ncm_matrix_get (matrix, 0, 0), ==, 1.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 0, 1), ==, 2.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 0, 2), ==, 3.0);

    g_assert_cmpfloat (ncm_matrix_get (matrix, 1, 0), ==, 4.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 1, 1), ==, 5.0);
    g_assert_cmpfloat (ncm_matrix_get (matrix, 1, 2), ==, 6.0);

    ncm_dtuple2_free (dt2);
    ncm_dtuple3_free (dt3);

    ncm_vector_free (vector);
    ncm_matrix_free (matrix);
  }
}

#endif /* HAVE_LIBFYAML */

void
test_ncm_serialize_traps (TestNcmSerialize *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/syntax/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/nonexistent/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/multiplenames/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/samename/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/wrongorder/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_serialize_global_invalid_from_string_syntax (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_global_from_string ("NcPowspecMLTransfer Error{'transfer':<('NcTransferFuncEH',@a{sv} {})>}");

  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_global_invalid_from_string_nonexist (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_global_from_string ("NcPowspecMLErrorName   {'transfer':<('NcTransferFuncEH',@a{sv} {})>}");

  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_invalid_from_string_mnames (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser, "NcPowspecMLTransfer[T0]{'transfer':<('NcTransferFuncEH[T0]', @a{sv} {})>}");

  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_invalid_from_string_samename (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser,
                                            "NcmSplineCubicNotaknot{'length':<6>, "
                                            "'x' : <('NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>})>, "
                                            "'y' : <('NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>})>}");

  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_invalid_from_string_wrongorder (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser,
                                            "NcmSplineCubicNotaknot{'length':<6>, "
                                            "'x' : <('NcmVector[T0]', @a{sv} {})>, "
                                            "'y' : <('NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>})>}");

  NCM_TEST_FREE (g_object_unref, obj);
}


/***************************************************************************
 *            test_ncm_dtuple.c
 *
 *  Fri December 08 10:28:16 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_dtuple.c
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

typedef struct _TestNcmDTuple2
{
  NcmDTuple2 *dt2;
} TestNcmDTuple2;

void test_ncm_dtuple2_new (TestNcmDTuple2 *test, gconstpointer pdata);
void test_ncm_dtuple2_free (TestNcmDTuple2 *test, gconstpointer pdata);

void test_ncm_dtuple2_copy (TestNcmDTuple2 *test, gconstpointer pdata);
void test_ncm_dtuple2_static (TestNcmDTuple2 *test, gconstpointer pdata);
void test_ncm_dtuple2_serialize (TestNcmDTuple2 *test, gconstpointer pdata);

typedef struct _TestNcmDTuple3
{
  NcmDTuple3 *dt3;
} TestNcmDTuple3;

void test_ncm_dtuple3_new (TestNcmDTuple3 *test, gconstpointer pdata);
void test_ncm_dtuple3_free (TestNcmDTuple3 *test, gconstpointer pdata);

void test_ncm_dtuple3_copy (TestNcmDTuple3 *test, gconstpointer pdata);
void test_ncm_dtuple3_static (TestNcmDTuple3 *test, gconstpointer pdata);
void test_ncm_dtuple3_serialize (TestNcmDTuple3 *test, gconstpointer pdata);


G_DECLARE_FINAL_TYPE (NcmObjectWithTuples, ncm_object_with_tuples, NCM, OBJECT_WITH_TUPLES, GObject)

struct _NcmObjectWithTuples
{
  GObject parent_instance;

  NcmDTuple2 *dt2;
  NcmDTuple3 *dt3;
};

enum
{
  PROP_0,
  PROP_DT2,
  PROP_DT3,
  PROP_LEN,
};

G_DEFINE_TYPE (NcmObjectWithTuples, ncm_object_with_tuples, G_TYPE_OBJECT)

static void
ncm_object_with_tuples_init (NcmObjectWithTuples *self)
{
  self->dt2 = NULL;
  self->dt3 = NULL;
}

static void
ncm_object_with_tuples_finalize (GObject *object)
{
  NcmObjectWithTuples *self = NCM_OBJECT_WITH_TUPLES (object);

  ncm_dtuple2_clear (&self->dt2);
  ncm_dtuple3_clear (&self->dt3);

  /* Chain up to the parent class */
  G_OBJECT_CLASS (ncm_object_with_tuples_parent_class)->finalize (object);
}

static void
ncm_object_with_tuples_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec   *pspec)
{
  NcmObjectWithTuples *self = NCM_OBJECT_WITH_TUPLES (object);

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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_object_with_tuples_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmObjectWithTuples *self = NCM_OBJECT_WITH_TUPLES (object);

  switch (prop_id)
  {
    case PROP_DT2:
      g_value_set_boxed (value, self->dt2);
      break;
    case PROP_DT3:
      g_value_set_boxed (value, self->dt3);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_object_with_tuples_class_init (NcmObjectWithTuplesClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize     = ncm_object_with_tuples_finalize;
  object_class->set_property = ncm_object_with_tuples_set_property;
  object_class->get_property = ncm_object_with_tuples_get_property;

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
}

typedef struct _TestNcmObjectWithTuples
{
  NcmObjectWithTuples *obj;
} TestNcmObjectWithTuples;


void test_ncm_object_with_tuples_new (TestNcmObjectWithTuples *test, gconstpointer pdata);
void test_ncm_object_with_tuples_free (TestNcmObjectWithTuples *test, gconstpointer pdata);

void test_ncm_object_with_tuples_serialize (TestNcmObjectWithTuples *test, gconstpointer pdata);
void test_ncm_object_with_tuples_from_string (TestNcmObjectWithTuples *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/dtuple2/copy", TestNcmDTuple2, NULL,
              &test_ncm_dtuple2_new,
              &test_ncm_dtuple2_copy,
              &test_ncm_dtuple2_free);

  g_test_add ("/ncm/dtuple2/static", TestNcmDTuple2, NULL,
              &test_ncm_dtuple2_new,
              &test_ncm_dtuple2_static,
              &test_ncm_dtuple2_free);

  g_test_add ("/ncm/dtuple2/serialize", TestNcmDTuple2, NULL,
              &test_ncm_dtuple2_new,
              &test_ncm_dtuple2_serialize,
              &test_ncm_dtuple2_free);

  g_test_add ("/ncm/dtuple3/copy", TestNcmDTuple3, NULL,
              &test_ncm_dtuple3_new,
              &test_ncm_dtuple3_copy,
              &test_ncm_dtuple3_free);

  g_test_add ("/ncm/dtuple3/static", TestNcmDTuple3, NULL,
              &test_ncm_dtuple3_new,
              &test_ncm_dtuple3_static,
              &test_ncm_dtuple3_free);

  g_test_add ("/ncm/dtuple3/serialize", TestNcmDTuple3, NULL,
              &test_ncm_dtuple3_new,
              &test_ncm_dtuple3_serialize,
              &test_ncm_dtuple3_free);

  g_test_add ("/ncm/object_with_tuples/serialize", TestNcmObjectWithTuples, NULL,
              &test_ncm_object_with_tuples_new,
              &test_ncm_object_with_tuples_serialize,
              &test_ncm_object_with_tuples_free);

  g_test_add ("/ncm/object_with_tuples/from_string", TestNcmObjectWithTuples, NULL,
              &test_ncm_object_with_tuples_new,
              &test_ncm_object_with_tuples_from_string,
              &test_ncm_object_with_tuples_free);

  g_test_run ();

  return 0;
}

void
test_ncm_dtuple2_new (TestNcmDTuple2 *test, gconstpointer pdata)
{
  const gdouble elem0 = g_test_rand_double ();
  const gdouble elem1 = g_test_rand_double ();

  test->dt2 = ncm_dtuple2_new (elem0, elem1);

  g_assert_cmpfloat (test->dt2->elements[0], ==, elem0);
  g_assert_cmpfloat (test->dt2->elements[1], ==, elem1);
}

void
test_ncm_dtuple2_free (TestNcmDTuple2 *test, gconstpointer pdata)
{
  ncm_dtuple2_free (test->dt2);
}

void
test_ncm_dtuple2_copy (TestNcmDTuple2 *test, gconstpointer pdata)
{
  NcmDTuple2 *dt2_copy = ncm_dtuple2_copy (test->dt2);

  g_assert_cmpfloat (test->dt2->elements[0], ==, dt2_copy->elements[0]);
  g_assert_cmpfloat (test->dt2->elements[1], ==, dt2_copy->elements[1]);

  ncm_dtuple2_clear (&dt2_copy);
  g_assert_null (dt2_copy);
}

void
test_ncm_dtuple2_static (TestNcmDTuple2 *test, gconstpointer pdata)
{
  NcmDTuple2 dt2 = NCM_DTUPLE2_STATIC_INIT (test->dt2->elements[0], test->dt2->elements[1]);

  g_assert_cmpfloat (dt2.elements[0], ==, test->dt2->elements[0]);
  g_assert_cmpfloat (dt2.elements[1], ==, test->dt2->elements[1]);
}

void
test_ncm_dtuple2_serialize (TestNcmDTuple2 *test, gconstpointer pdata)
{
  GVariant *var   = ncm_dtuple2_serialize (test->dt2);
  NcmDTuple2 *dt2 = ncm_dtuple2_new_from_variant (var);

  g_assert_cmpfloat_with_epsilon (test->dt2->elements[0], dt2->elements[0], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->dt2->elements[1], dt2->elements[1], GSL_DBL_EPSILON);

  ncm_dtuple2_free (dt2);
  g_variant_unref (var);
}

void
test_ncm_dtuple3_new (TestNcmDTuple3 *test, gconstpointer pdata)
{
  const gdouble elem0 = g_test_rand_double ();
  const gdouble elem1 = g_test_rand_double ();
  const gdouble elem2 = g_test_rand_double ();

  test->dt3 = ncm_dtuple3_new (elem0, elem1, elem2);

  g_assert_cmpfloat (test->dt3->elements[0], ==, elem0);
  g_assert_cmpfloat (test->dt3->elements[1], ==, elem1);
  g_assert_cmpfloat (test->dt3->elements[2], ==, elem2);
}

void
test_ncm_dtuple3_free (TestNcmDTuple3 *test, gconstpointer pdata)
{
  ncm_dtuple3_free (test->dt3);
}

void
test_ncm_dtuple3_copy (TestNcmDTuple3 *test, gconstpointer pdata)
{
  NcmDTuple3 *dt3_copy = ncm_dtuple3_copy (test->dt3);

  g_assert_cmpfloat (test->dt3->elements[0], ==, dt3_copy->elements[0]);
  g_assert_cmpfloat (test->dt3->elements[1], ==, dt3_copy->elements[1]);
  g_assert_cmpfloat (test->dt3->elements[2], ==, dt3_copy->elements[2]);

  ncm_dtuple3_clear (&dt3_copy);
  g_assert_null (dt3_copy);
}

void
test_ncm_dtuple3_static (TestNcmDTuple3 *test, gconstpointer pdata)
{
  NcmDTuple3 dt3 = NCM_DTUPLE3_STATIC_INIT (test->dt3->elements[0], test->dt3->elements[1], test->dt3->elements[2]);

  g_assert_cmpfloat (dt3.elements[0], ==, test->dt3->elements[0]);
  g_assert_cmpfloat (dt3.elements[1], ==, test->dt3->elements[1]);
  g_assert_cmpfloat (dt3.elements[2], ==, test->dt3->elements[2]);
}

void
test_ncm_dtuple3_serialize (TestNcmDTuple3 *test, gconstpointer pdata)
{
  GVariant *var   = ncm_dtuple3_serialize (test->dt3);
  NcmDTuple3 *dt3 = ncm_dtuple3_new_from_variant (var);

  g_assert_cmpfloat_with_epsilon (test->dt3->elements[0], dt3->elements[0], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->dt3->elements[1], dt3->elements[1], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->dt3->elements[2], dt3->elements[2], GSL_DBL_EPSILON);

  ncm_dtuple3_free (dt3);
  g_variant_unref (var);
}

void
test_ncm_object_with_tuples_new (TestNcmObjectWithTuples *test, gconstpointer pdata)
{
  const gdouble elem0 = g_test_rand_double ();
  const gdouble elem1 = g_test_rand_double ();
  const gdouble elem2 = g_test_rand_double ();
  const gdouble elem3 = g_test_rand_double ();
  const gdouble elem4 = g_test_rand_double ();
  const gdouble elem5 = g_test_rand_double ();
  NcmDTuple2 *dt2     = ncm_dtuple2_new (elem0, elem1);
  NcmDTuple3 *dt3     = ncm_dtuple3_new (elem2, elem3, elem4);

  test->obj = g_object_new (ncm_object_with_tuples_get_type (),
                            "dt2", dt2,
                            "dt3", dt3,
                            NULL);

  g_assert_cmpfloat (test->obj->dt2->elements[0], ==, elem0);
  g_assert_cmpfloat (test->obj->dt2->elements[1], ==, elem1);
  g_assert_cmpfloat (test->obj->dt3->elements[0], ==, elem2);
  g_assert_cmpfloat (test->obj->dt3->elements[1], ==, elem3);
  g_assert_cmpfloat (test->obj->dt3->elements[2], ==, elem4);

  ncm_dtuple2_free (dt2);
  ncm_dtuple3_free (dt3);
}

void
test_ncm_object_with_tuples_free (TestNcmObjectWithTuples *test, gconstpointer pdata)
{
  NCM_TEST_FREE (g_object_unref, test->obj);
}

void
test_ncm_object_with_tuples_serialize (TestNcmObjectWithTuples *test, gconstpointer pdata)
{
  NcmSerialize *ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  GVariant *var            = ncm_serialize_to_variant (ser, G_OBJECT (test->obj));
  NcmObjectWithTuples *obj = NCM_OBJECT_WITH_TUPLES (ncm_serialize_from_variant (ser, var));

  g_assert_cmpfloat_with_epsilon (test->obj->dt2->elements[0], obj->dt2->elements[0], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt2->elements[1], obj->dt2->elements[1], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt3->elements[0], obj->dt3->elements[0], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt3->elements[1], obj->dt3->elements[1], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt3->elements[2], obj->dt3->elements[2], GSL_DBL_EPSILON);

  g_object_unref (obj);
  g_variant_unref (var);
  ncm_serialize_free (ser);
}

void
test_ncm_object_with_tuples_from_string (TestNcmObjectWithTuples *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *string     = g_strdup_printf ("('NcmObjectWithTuples', {'dt2': <(%24.16g,%24.16g)>, 'dt3': <(%24.16g,%24.16g,%24.16g)>})",
                                       test->obj->dt2->elements[0],
                                       test->obj->dt2->elements[1],
                                       test->obj->dt3->elements[0],
                                       test->obj->dt3->elements[1],
                                       test->obj->dt3->elements[2]);
  NcmObjectWithTuples *obj = NCM_OBJECT_WITH_TUPLES (ncm_serialize_from_string (ser, string));

  g_assert_cmpfloat_with_epsilon (test->obj->dt2->elements[0], obj->dt2->elements[0], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt2->elements[1], obj->dt2->elements[1], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt3->elements[0], obj->dt3->elements[0], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt3->elements[1], obj->dt3->elements[1], GSL_DBL_EPSILON);
  g_assert_cmpfloat_with_epsilon (test->obj->dt3->elements[2], obj->dt3->elements[2], GSL_DBL_EPSILON);

  g_object_unref (obj);
  g_free (string);
  ncm_serialize_free (ser);
}


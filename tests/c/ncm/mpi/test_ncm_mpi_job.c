/***************************************************************************
 *            test_ncm_mpi_job.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mpi_job.c
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
#include <string.h>

typedef void (*TestNcmMPIJobRunArray) (NcmMPIJob *mpi_job, GPtrArray *input_array, GPtrArray *ret_array);

/*
 * TestMPIJobShape: a job whose input and return messages have different lengths and
 * travel in the default pooled buffers of NcmMPIJob, which no library job uses. The
 * return is ret[k] = (k + 1) sum_j (j + 1) input[j]. With float-input the input travels
 * as MPI_FLOAT, so the return message is longer in bytes than return-len elements of
 * the input datatype.
 */
#define TEST_TYPE_MPI_JOB_SHAPE (test_mpi_job_shape_get_type ())
G_DECLARE_FINAL_TYPE (TestMPIJobShape, test_mpi_job_shape, TEST, MPI_JOB_SHAPE, NcmMPIJob)

struct _TestMPIJobShape
{
  NcmMPIJob parent_instance;
  guint input_len;
  guint return_len;
  gboolean float_input;
};

enum
{
  PROP_0,
  PROP_INPUT_LEN,
  PROP_RETURN_LEN,
  PROP_FLOAT_INPUT,
};

G_DEFINE_TYPE (TestMPIJobShape, test_mpi_job_shape, NCM_TYPE_MPI_JOB)

static void
test_mpi_job_shape_init (TestMPIJobShape *mjs)
{
  mjs->input_len   = 0;
  mjs->return_len  = 0;
  mjs->float_input = FALSE;
}

static void
_test_mpi_job_shape_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  TestMPIJobShape *mjs = TEST_MPI_JOB_SHAPE (object);

  switch (prop_id)
  {
    case PROP_INPUT_LEN:
      mjs->input_len = g_value_get_uint (value);
      break;
    case PROP_RETURN_LEN:
      mjs->return_len = g_value_get_uint (value);
      break;
    case PROP_FLOAT_INPUT:
      mjs->float_input = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_test_mpi_job_shape_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  TestMPIJobShape *mjs = TEST_MPI_JOB_SHAPE (object);

  switch (prop_id)
  {
    case PROP_INPUT_LEN:
      g_value_set_uint (value, mjs->input_len);
      break;
    case PROP_RETURN_LEN:
      g_value_set_uint (value, mjs->return_len);
      break;
    case PROP_FLOAT_INPUT:
      g_value_set_boolean (value, mjs->float_input);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static NcmMPIDatatype
_test_mpi_job_shape_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  TestMPIJobShape *mjs = TEST_MPI_JOB_SHAPE (mpi_job);

  len[0] = mjs->input_len;

  if (mjs->float_input)
  {
    size[0] = sizeof (gfloat) * mjs->input_len;

    return MPI_FLOAT;
  }

  size[0] = sizeof (gdouble) * mjs->input_len;

  return MPI_DOUBLE;
}

static NcmMPIDatatype
_test_mpi_job_shape_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  TestMPIJobShape *mjs = TEST_MPI_JOB_SHAPE (mpi_job);

  len[0]  = mjs->return_len;
  size[0] = sizeof (gdouble) * mjs->return_len;

  return MPI_DOUBLE;
}

static gpointer
_test_mpi_job_shape_create_input (NcmMPIJob *mpi_job)
{
  return ncm_vector_new (TEST_MPI_JOB_SHAPE (mpi_job)->input_len);
}

static gpointer
_test_mpi_job_shape_create_return (NcmMPIJob *mpi_job)
{
  return ncm_vector_new (TEST_MPI_JOB_SHAPE (mpi_job)->return_len);
}

static void
_test_mpi_job_shape_destroy_vector (NcmMPIJob *mpi_job, gpointer v)
{
  ncm_vector_free (v);
}

static gpointer
_test_mpi_job_shape_pack_input (NcmMPIJob *mpi_job, gpointer input)
{
  gpointer buf = ncm_mpi_job_get_input_buffer (mpi_job, input);
  guint j;

  if (TEST_MPI_JOB_SHAPE (mpi_job)->float_input)
    for (j = 0; j < ncm_vector_len (input); j++)
      ((gfloat *) buf)[j] = ncm_vector_get (input, j);

  else
    memcpy (buf, ncm_vector_data (input), sizeof (gdouble) * ncm_vector_len (input));

  return buf;
}

static void
_test_mpi_job_shape_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input)
{
  guint j;

  if (TEST_MPI_JOB_SHAPE (mpi_job)->float_input)
    for (j = 0; j < ncm_vector_len (input); j++)
      ncm_vector_set (input, j, ((gfloat *) buf)[j]);

  else
    memcpy (ncm_vector_data (input), buf, sizeof (gdouble) * ncm_vector_len (input));
}

static gpointer
_test_mpi_job_shape_pack_return (NcmMPIJob *mpi_job, gpointer ret)
{
  gdouble *buf = ncm_mpi_job_get_return_buffer (mpi_job, ret);

  memcpy (buf, ncm_vector_data (ret), sizeof (gdouble) * ncm_vector_len (ret));

  return buf;
}

static void
_test_mpi_job_shape_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret)
{
  memcpy (ncm_vector_data (ret), buf, sizeof (gdouble) * ncm_vector_len (ret));
}

static void
_test_mpi_job_shape_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret)
{
  gdouble s = 0.0;
  guint j, k;

  for (j = 0; j < ncm_vector_len (input); j++)
    s += (j + 1.0) * ncm_vector_get (input, j);

  for (k = 0; k < ncm_vector_len (ret); k++)
    ncm_vector_set (ret, k, (k + 1.0) * s);
}

static void
test_mpi_job_shape_class_init (TestMPIJobShapeClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcmMPIJobClass *mpi_job_class = NCM_MPI_JOB_CLASS (klass);

  object_class->set_property = &_test_mpi_job_shape_set_property;
  object_class->get_property = &_test_mpi_job_shape_get_property;

  g_object_class_install_property (object_class, PROP_INPUT_LEN,
                                   g_param_spec_uint ("input-len", NULL, "Input length", 1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));
  g_object_class_install_property (object_class, PROP_RETURN_LEN,
                                   g_param_spec_uint ("return-len", NULL, "Return length", 1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));
  g_object_class_install_property (object_class, PROP_FLOAT_INPUT,
                                   g_param_spec_boolean ("float-input", NULL, "Input as MPI_FLOAT", FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  mpi_job_class->input_datatype  = &_test_mpi_job_shape_input_datatype;
  mpi_job_class->return_datatype = &_test_mpi_job_shape_return_datatype;
  mpi_job_class->create_input    = &_test_mpi_job_shape_create_input;
  mpi_job_class->create_return   = &_test_mpi_job_shape_create_return;
  mpi_job_class->destroy_input   = &_test_mpi_job_shape_destroy_vector;
  mpi_job_class->destroy_return  = &_test_mpi_job_shape_destroy_vector;
  mpi_job_class->pack_input      = &_test_mpi_job_shape_pack_input;
  mpi_job_class->pack_return     = &_test_mpi_job_shape_pack_return;
  mpi_job_class->unpack_input    = &_test_mpi_job_shape_unpack_input;
  mpi_job_class->unpack_return   = &_test_mpi_job_shape_unpack_return;
  mpi_job_class->run             = &_test_mpi_job_shape_run;
}

/*
 * NcmMPIJobTest returns the entry of its vector at the index it receives. Each output
 * must be the entry at the index of its own input, whichever rank ran it.
 */
static void
_test_ncm_mpi_job_run_test_job (TestNcmMPIJobRunArray run_array, const guint len)
{
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, 20260930);
  NcmMPIJobTest *mjt = ncm_mpi_job_test_new ();
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GPtrArray *input_a = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  GPtrArray *ret_a   = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  NcmVector *vec     = NULL;
  guint i;

  ncm_mpi_job_test_set_rand_vector (mjt, 97, rng);
  g_object_get (mjt, "vector", &vec, NULL);

  for (i = 0; i < len; i++)
  {
    NcmVector *input = ncm_vector_new (1);

    ncm_vector_set (input, 0, (7 * i) % 97);
    g_ptr_array_add (input_a, input);
    g_ptr_array_add (ret_a, ncm_vector_new (1));
  }

  ncm_mpi_job_init_all_slaves (NCM_MPI_JOB (mjt), ser);
  run_array (NCM_MPI_JOB (mjt), input_a, ret_a);
  ncm_mpi_job_free_all_slaves (NCM_MPI_JOB (mjt));

  for (i = 0; i < len; i++)
    g_assert_cmpfloat (ncm_vector_get (g_ptr_array_index (ret_a, i), 0), ==, ncm_vector_get (vec, (7 * i) % 97));

  ncm_vector_free (vec);
  g_ptr_array_unref (input_a);
  g_ptr_array_unref (ret_a);
  ncm_serialize_free (ser);
  ncm_mpi_job_test_free (mjt);
  ncm_rng_free (rng);
}

/* The shape job on len inputs, input i with entries i + j / 2. */
static void
_test_ncm_mpi_job_run_shape_job_full (TestNcmMPIJobRunArray run_array, const guint input_len, const guint return_len, const guint len, gboolean float_input)
{
  NcmMPIJob *mpi_job = g_object_new (TEST_TYPE_MPI_JOB_SHAPE,
                                     "input-len", input_len,
                                     "return-len", return_len,
                                     "float-input", float_input,
                                     NULL);
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GPtrArray *input_a = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  GPtrArray *ret_a   = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  guint i, j, k;

  for (i = 0; i < len; i++)
  {
    NcmVector *input = ncm_vector_new (input_len);

    for (j = 0; j < input_len; j++)
      ncm_vector_set (input, j, i + 0.5 * j);

    g_ptr_array_add (input_a, input);
    g_ptr_array_add (ret_a, ncm_vector_new (return_len));
  }

  ncm_mpi_job_init_all_slaves (mpi_job, ser);
  run_array (mpi_job, input_a, ret_a);
  ncm_mpi_job_free_all_slaves (mpi_job);

  for (i = 0; i < len; i++)
  {
    gdouble s = 0.0;

    for (j = 0; j < input_len; j++)
      s += (j + 1.0) * (i + 0.5 * j);

    for (k = 0; k < return_len; k++)
      g_assert_cmpfloat (ncm_vector_get (g_ptr_array_index (ret_a, i), k), ==, (k + 1.0) * s);
  }

  g_ptr_array_unref (input_a);
  g_ptr_array_unref (ret_a);
  ncm_serialize_free (ser);
  ncm_mpi_job_free (mpi_job);
}

static void
_test_ncm_mpi_job_run_shape_job (TestNcmMPIJobRunArray run_array, const guint input_len, const guint return_len, const guint len)
{
  _test_ncm_mpi_job_run_shape_job_full (run_array, input_len, return_len, len, FALSE);
}

static void
test_ncm_mpi_job_run_array (void)
{
  _test_ncm_mpi_job_run_test_job (&ncm_mpi_job_run_array, 97);
}

static void
test_ncm_mpi_job_run_array_async (void)
{
  _test_ncm_mpi_job_run_test_job (&ncm_mpi_job_run_array_async, 97);
}

static void
test_ncm_mpi_job_shape_run_array (void)
{
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array, 3, 5, 53);
}

static void
test_ncm_mpi_job_shape_run_array_async (void)
{
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array_async, 3, 5, 53);
}

static void
test_ncm_mpi_job_float_input (void)
{
  /* Every value is exact in single precision. */
  _test_ncm_mpi_job_run_shape_job_full (&ncm_mpi_job_run_array, 3, 5, 53, TRUE);
  _test_ncm_mpi_job_run_shape_job_full (&ncm_mpi_job_run_array_async, 3, 5, 53, TRUE);
}

static void
test_ncm_mpi_job_empty (void)
{
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array, 3, 5, 0);
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array_async, 3, 5, 0);
}

static void
test_ncm_mpi_job_single_input (void)
{
  /* Fewer inputs than workers from three ranks up. */
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array, 3, 5, 1);
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array_async, 3, 5, 1);
}

static void
test_ncm_mpi_job_sequence (void)
{
  /* The workers serve one job after another, of either type and shape. */
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array, 3, 5, 20);
  _test_ncm_mpi_job_run_test_job (&ncm_mpi_job_run_array_async, 20);
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array_async, 7, 2, 20);
  _test_ncm_mpi_job_run_shape_job (&ncm_mpi_job_run_array, 1, 11, 20);
}

gint
main (gint argc, gchar *argv[])
{
  /* The workers enter their loop inside ncm_cfg_init and deserialize the jobs they are
   * sent, so the test job type must exist before it. */
  g_type_ensure (TEST_TYPE_MPI_JOB_SHAPE);

  ncm_cfg_init_full_ptr (&argc, &argv);
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/mpi/job/run_array", &test_ncm_mpi_job_run_array);
  g_test_add_func ("/ncm/mpi/job/run_array_async", &test_ncm_mpi_job_run_array_async);
  g_test_add_func ("/ncm/mpi/job/shape/run_array", &test_ncm_mpi_job_shape_run_array);
  g_test_add_func ("/ncm/mpi/job/shape/run_array_async", &test_ncm_mpi_job_shape_run_array_async);
  g_test_add_func ("/ncm/mpi/job/float_input", &test_ncm_mpi_job_float_input);
  g_test_add_func ("/ncm/mpi/job/empty", &test_ncm_mpi_job_empty);
  g_test_add_func ("/ncm/mpi/job/single_input", &test_ncm_mpi_job_single_input);
  g_test_add_func ("/ncm/mpi/job/sequence", &test_ncm_mpi_job_sequence);

  g_test_run ();
}


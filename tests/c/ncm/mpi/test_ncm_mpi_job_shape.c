/***************************************************************************
 *            test_ncm_mpi_job_shape.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mpi_job_shape.c
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
#include <string.h>

#include "test_ncm_mpi_job_shape.h"

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


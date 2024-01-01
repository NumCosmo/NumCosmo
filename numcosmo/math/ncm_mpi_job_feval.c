/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_mpi_job_feval.c
 *
 *  Fri February 21 11:16:32 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_feval.c
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/**
 * SECTION:ncm_mpi_job_feval
 * @title: NcmMPIJobFEval
 * @short_description: MPI job object for evaluating fit steps
 *
 * This object is a subclass of #NcmMPIJob and is designed to implement an MPI job
 * for evaluating fit steps. It is employed by #NcmFit to parallelize the evaluation
 * of the posterior function. The job entails computing the posterior function and,
 * if applicable, additional functions (e.g., derived quantities) at a specified
 * point within the parameter space.
 *
 * The MPI job is implemented as a function that takes a vector of parameters as input
 * and produces a vector of values as output. The first value represents the posterior
 * function's value, while the subsequent values correspond to those of additional functions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm_enum_types.h"
#include "math/ncm_mpi_job_feval.h"

#ifndef HAVE_MPI
#define MPI_DATATYPE_NULL (0)
#define MPI_DOUBLE (0)
#endif /* HAVE_MPI */

typedef struct _NcmMPIJobFEvalPrivate
{
  NcmFit *fit;
  NcmObjArray *func_oa;
  gint fparam_len;
  gint nadd_vals;
} NcmMPIJobFEvalPrivate;

enum
{
  PROP_0,
  PROP_FIT,
  PROP_FUNC_ARRAY,
  PROP_JOB_TYPE,
};

struct _NcmMPIJobFEval
{
  NcmMPIJob parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmMPIJobFEval, ncm_mpi_job_feval, NCM_TYPE_MPI_JOB)

static void
ncm_mpi_job_feval_init (NcmMPIJobFEval *mjfeval)
{
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  self->fit     = NULL;
  self->func_oa = NULL;

  self->fparam_len = 0;
  self->nadd_vals  = 0;
}

static void
_ncm_mpi_job_feval_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (object);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  g_return_if_fail (NCM_IS_MPI_JOB_FEVAL (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_assert (self->fit == NULL);
      self->fit = g_value_dup_object (value);
      g_assert (self->fit != NULL);
      break;
    case PROP_FUNC_ARRAY:
    {
      ncm_obj_array_clear (&self->func_oa);
      self->func_oa = g_value_dup_boxed (value);

      if (self->func_oa != NULL)
      {
        guint i;

        for (i = 0; i < self->func_oa->len; i++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, i));

          g_assert (NCM_IS_MSET_FUNC (func));
          g_assert (ncm_mset_func_is_scalar (func));
          g_assert (ncm_mset_func_is_const (func));
        }
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mpi_job_feval_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (object);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  g_return_if_fail (NCM_IS_MPI_JOB_FEVAL (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, self->fit);
      break;
    case PROP_FUNC_ARRAY:
      g_value_set_boxed (value, self->func_oa);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mpi_job_feval_constructed (GObject *object)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (object);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);
  NcmMSet *mset                      = ncm_fit_peek_mset (self->fit);

  self->fparam_len = ncm_mset_fparam_len (mset);
  self->nadd_vals  = 1 + ((self->func_oa != NULL) ? self->func_oa->len : 0);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mpi_job_feval_parent_class)->constructed (object);
}

static void
_ncm_mpi_job_feval_dispose (GObject *object)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (object);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  ncm_fit_clear (&self->fit);
  ncm_obj_array_clear (&self->func_oa);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mpi_job_feval_parent_class)->dispose (object);
}

static void
_ncm_mpi_job_feval_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mpi_job_feval_parent_class)->finalize (object);
}

static NcmMPIDatatype _ncm_mpi_job_feval_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);
static NcmMPIDatatype _ncm_mpi_job_feval_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);

static gpointer _ncm_mpi_job_feval_create_input (NcmMPIJob *mpi_job);
static gpointer _ncm_mpi_job_feval_create_return (NcmMPIJob *mpi_job);

static void _ncm_mpi_job_feval_destroy_input (NcmMPIJob *mpi_job, gpointer input);
static void _ncm_mpi_job_feval_destroy_return (NcmMPIJob *mpi_job, gpointer ret);

static gpointer _ncm_mpi_job_feval_get_input_buffer (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_feval_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_feval_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf);
static void _ncm_mpi_job_feval_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf);

static gpointer _ncm_mpi_job_feval_pack_input (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_feval_pack_return (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_feval_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input);
static void _ncm_mpi_job_feval_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret);

static void _ncm_mpi_job_feval_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret);

static void
ncm_mpi_job_feval_class_init (NcmMPIJobFEvalClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcmMPIJobClass *mpi_job_class = NCM_MPI_JOB_CLASS (klass);

  object_class->set_property = &_ncm_mpi_job_feval_set_property;
  object_class->get_property = &_ncm_mpi_job_feval_get_property;
  object_class->constructed  = &_ncm_mpi_job_feval_constructed;
  object_class->dispose      = &_ncm_mpi_job_feval_dispose;
  object_class->finalize     = &_ncm_mpi_job_feval_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FUNC_ARRAY,
                                   g_param_spec_boxed ("function-array",
                                                       NULL,
                                                       "Functions array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  mpi_job_class->input_datatype  = &_ncm_mpi_job_feval_input_datatype;
  mpi_job_class->return_datatype = &_ncm_mpi_job_feval_return_datatype;

  mpi_job_class->create_input  = &_ncm_mpi_job_feval_create_input;
  mpi_job_class->create_return = &_ncm_mpi_job_feval_create_return;

  mpi_job_class->destroy_input  = &_ncm_mpi_job_feval_destroy_input;
  mpi_job_class->destroy_return = &_ncm_mpi_job_feval_destroy_return;

  mpi_job_class->get_input_buffer  = &_ncm_mpi_job_feval_get_input_buffer;
  mpi_job_class->get_return_buffer = &_ncm_mpi_job_feval_get_return_buffer;

  mpi_job_class->destroy_input_buffer  = &_ncm_mpi_job_feval_destroy_input_buffer;
  mpi_job_class->destroy_return_buffer = &_ncm_mpi_job_feval_destroy_return_buffer;

  mpi_job_class->pack_input  = &_ncm_mpi_job_feval_pack_input;
  mpi_job_class->pack_return = &_ncm_mpi_job_feval_pack_return;

  mpi_job_class->unpack_input  = &_ncm_mpi_job_feval_unpack_input;
  mpi_job_class->unpack_return = &_ncm_mpi_job_feval_unpack_return;

  mpi_job_class->run = &_ncm_mpi_job_feval_run;
}

static NcmMPIDatatype
_ncm_mpi_job_feval_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (mpi_job);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  len[0]  = self->fparam_len;
  size[0] = sizeof (gdouble) * len[0];

  return MPI_DOUBLE;
}

static NcmMPIDatatype
_ncm_mpi_job_feval_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (mpi_job);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  len[0]  = self->nadd_vals;
  size[0] = sizeof (gdouble) * len[0];

  return MPI_DOUBLE;
}

static gpointer
_ncm_mpi_job_feval_create_input (NcmMPIJob *mpi_job)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (mpi_job);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  return ncm_vector_new (self->fparam_len);
}

static gpointer
_ncm_mpi_job_feval_create_return (NcmMPIJob *mpi_job)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (mpi_job);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);

  return ncm_vector_new (self->nadd_vals);
}

static void
_ncm_mpi_job_feval_destroy_input (NcmMPIJob *mpi_job, gpointer input)
{
  ncm_vector_free (input);
}

static void
_ncm_mpi_job_feval_destroy_return (NcmMPIJob *mpi_job, gpointer ret)
{
  ncm_vector_free (ret);
}

static gpointer
_ncm_mpi_job_feval_get_input_buffer (NcmMPIJob *mpi_job, gpointer input)
{
  return ncm_vector_data (input);
}

static gpointer
_ncm_mpi_job_feval_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret)
{
  return ncm_vector_data (ret);
}

static void
_ncm_mpi_job_feval_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf)
{
  g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (input)), ==, GPOINTER_TO_INT (buf));
}

static void
_ncm_mpi_job_feval_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf)
{
  g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (ret)), ==, GPOINTER_TO_INT (buf));
}

static gpointer
_ncm_mpi_job_feval_pack_input (NcmMPIJob *mpi_job, gpointer input)
{
  return ncm_vector_data (input);
}

static gpointer
_ncm_mpi_job_feval_pack_return (NcmMPIJob *mpi_job, gpointer ret)
{
  return ncm_vector_data (ret);
}

static void
_ncm_mpi_job_feval_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input)
{
  g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (input)), ==, GPOINTER_TO_INT (buf));
}

static void
_ncm_mpi_job_feval_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret)
{
  g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (ret)), ==, GPOINTER_TO_INT (buf));
}

static void
_ncm_mpi_job_feval_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret)
{
  NcmMPIJobFEval *mjfeval            = NCM_MPI_JOB_FEVAL (mpi_job);
  NcmMPIJobFEvalPrivate * const self = ncm_mpi_job_feval_get_instance_private (mjfeval);
  gdouble *m2lnL_star                = ncm_vector_ptr (ret, 0);
  NcmMSet *mset                      = ncm_fit_peek_mset (self->fit);

  ncm_fit_params_set_vector (self->fit, input);
  ncm_fit_m2lnL_val (self->fit, m2lnL_star);

  if (self->func_oa != NULL)
  {
    guint i;

    for (i = 0; i < self->func_oa->len; i++)
    {
      NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, i));
      const gdouble a_i = ncm_mset_func_eval0 (func, mset);

      ncm_vector_set (ret, i + 1, a_i);
    }
  }
}

/**
 * ncm_mpi_job_feval_new:
 * @fit: a #NcmFit
 * @func_oa: (array) (element-type NcmMSetFunc) (allow-none): a #NcmObjArray
 *
 * Creates a new #NcmMPIJobFEval object.
 *
 * Returns: a new #NcmMPIJobFEval.
 */
NcmMPIJobFEval *
ncm_mpi_job_feval_new (NcmFit *fit, NcmObjArray *func_oa)
{
  NcmMPIJobFEval *mjfeval = g_object_new (NCM_TYPE_MPI_JOB_FEVAL,
                                          "fit",            fit,
                                          "function-array", func_oa,
                                          NULL);

  return mjfeval;
}

/**
 * ncm_mpi_job_feval_ref:
 * @mjfeval: a #NcmMPIJobFEval
 *
 * Increase the reference of @mjfeval by one.
 *
 * Returns: (transfer full): @mjfeval.
 */
NcmMPIJobFEval *
ncm_mpi_job_feval_ref (NcmMPIJobFEval *mjfeval)
{
  return g_object_ref (mjfeval);
}

/**
 * ncm_mpi_job_feval_free:
 * @mjfeval: a #NcmMPIJobFEval
 *
 * Decrease the reference count of @mjfeval by one.
 *
 */
void
ncm_mpi_job_feval_free (NcmMPIJobFEval *mjfeval)
{
  g_object_unref (mjfeval);
}

/**
 * ncm_mpi_job_feval_clear:
 * @mjfeval: a #NcmMPIJobFEval
 *
 * Decrease the reference count of @mjfeval by one, and sets the pointer *@mjfeval to
 * NULL.
 *
 */
void
ncm_mpi_job_feval_clear (NcmMPIJobFEval **mjfeval)
{
  g_clear_object (mjfeval);
}


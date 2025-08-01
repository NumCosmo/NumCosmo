/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_mpi_job.c
 *
 *  Sat April 21 18:33:19 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mpi_job.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmMPIJob:
 *
 * Abstract class to implement MPI jobs.
 *
 * This abstract class simplifies the implementation of MPI jobs through a master/slave
 * model. Subclasses must implement virtual methods. The master dispatches jobs to the
 * slaves, awaits results, and slaves execute the job, sending results back to the
 * master.
 *
 * For example, NcmMPIJobMCMC is a specific subclass that implements an MCMC job. Each
 * slave creates an instance of #NcmLikelihood, evaluating the likelihood function for
 * samples received from the master.
 *
 * The master/slave model leverages MPI. The subclass is responsible for implementing
 * virtual methods for packing/unpacking input and return objects into/from MPI buffers.
 *
 * When NumCosmo is compiled with MPI support, the rank 0 process executes normally,
 * while other ranks wait for commands from the master. The method
 * ncm_mpi_job_init_all_slaves() sends a serialized version of itself to all slaves.
 * Slaves deserialize it and wait for commands.
 *
 * Each call to ncm_mpi_job_run_array() sends an array of inputs to the slaves in a
 * round-robin fashion. Slaves execute the job and return results to the master,
 * maintaining the input order. This method uses rank 0 to control the slaves, so the
 * master isn't involved in computations.
 *
 * Conversely, ncm_mpi_job_run_array_async() sends inputs to the slaves and creates a
 * lightweight thread to receive results while the master concurrently executes the job.
 *
 * After job completion, the master calls ncm_mpi_job_free_all_slaves() to release the
 * slaves. It's an error to use two or more instances of any subclass of #NcmMPIJob
 * simultaneously.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mpi_job.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_timer.h"
#include "math/ncm_util.h"

#ifndef HAVE_MPI
#define MPI_DATATYPE_NULL (0)
#define MPI_DOUBLE (0)
#endif /* HAVE_MPI */

typedef struct _NcmMPIJobPrivate
{
  guint placeholder;
  NcmMPIDatatype input_dtype;
  NcmMPIDatatype return_dtype;
  gint input_size;
  gint return_size;
  gint input_len;
  gint return_len;
  gint owned_slaves;
  NcmMemoryPool *input_buf_pool;
  NcmMemoryPool *return_buf_pool;
  GHashTable *input_buf_table;
  GHashTable *return_buf_table;
  NcmTimer *nt;
} NcmMPIJobPrivate;

enum
{
  PROP_0,
  PROP_PLACEHOLDER,
};

extern NcmMPIJobCtrl _mpi_ctrl;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmMPIJob, ncm_mpi_job, G_TYPE_OBJECT)

static gpointer _ncm_mpi_job_create_input_buffer (gpointer userdata);
static gpointer _ncm_mpi_job_create_return_buffer (gpointer userdata);
static void _ncm_mpi_job_destroy_buffer (gpointer p);

static void
ncm_mpi_job_init (NcmMPIJob *mpi_job)
{
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  self->placeholder      = 0;
  self->input_dtype      = MPI_DATATYPE_NULL;
  self->return_dtype     = MPI_DATATYPE_NULL;
  self->input_size       = 0;
  self->return_size      = 0;
  self->input_len        = 0;
  self->return_len       = 0;
  self->owned_slaves     = 0;
  self->input_buf_pool   = ncm_memory_pool_new (_ncm_mpi_job_create_input_buffer,  mpi_job, _ncm_mpi_job_destroy_buffer);
  self->return_buf_pool  = ncm_memory_pool_new (_ncm_mpi_job_create_return_buffer, mpi_job, _ncm_mpi_job_destroy_buffer);
  self->input_buf_table  = g_hash_table_new (g_direct_hash, g_direct_equal);
  self->return_buf_table = g_hash_table_new (g_direct_hash, g_direct_equal);
  self->nt               = ncm_timer_new ();
}

static gpointer
_ncm_mpi_job_create_input_buffer (gpointer userdata)
{
  NcmMPIJob *mpi_job            = NCM_MPI_JOB (userdata);
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  return g_malloc (self->input_size);
}

static gpointer
_ncm_mpi_job_create_return_buffer (gpointer userdata)
{
  NcmMPIJob *mpi_job            = NCM_MPI_JOB (userdata);
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  return g_malloc (self->return_size);
}

static void
_ncm_mpi_job_destroy_buffer (gpointer p)
{
  g_free (p);
}

static void
_ncm_mpi_job_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMPIJob *mpi_job            = NCM_MPI_JOB (object);
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  g_return_if_fail (NCM_IS_MPI_JOB (object));

  switch (prop_id)
  {
    case PROP_PLACEHOLDER:
      self->placeholder = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mpi_job_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMPIJob *mpi_job            = NCM_MPI_JOB (object);
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  g_return_if_fail (NCM_IS_MPI_JOB (object));

  switch (prop_id)
  {
    case PROP_PLACEHOLDER:
      g_value_set_uint (value, self->placeholder);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mpi_job_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_mpi_job_parent_class)->constructed (object);
  {
    NcmMPIJob *mpi_job            = NCM_MPI_JOB (object);
    NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

    self->input_dtype  = ncm_mpi_job_input_datatype  (mpi_job, &self->input_len,  &self->input_size);
    self->return_dtype = ncm_mpi_job_return_datatype (mpi_job, &self->return_len, &self->return_size);
  }
}

static void
_ncm_mpi_job_dispose (GObject *object)
{
  NcmMPIJob *mpi_job            = NCM_MPI_JOB (object);
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  g_clear_pointer (&self->input_buf_table, g_hash_table_unref);
  g_clear_pointer (&self->return_buf_table, g_hash_table_unref);

  ncm_timer_clear (&self->nt);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mpi_job_parent_class)->dispose (object);
}

static void
_ncm_mpi_job_finalize (GObject *object)
{
  NcmMPIJob *mpi_job            = NCM_MPI_JOB (object);
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

  ncm_mpi_job_free_all_slaves (mpi_job);

  ncm_memory_pool_free (self->input_buf_pool, TRUE);
  ncm_memory_pool_free (self->return_buf_pool, TRUE);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mpi_job_parent_class)->finalize (object);
}

static NcmMPIDatatype
_ncm_mpi_job_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  g_error ("_ncm_mpi_job_input_datatype: method not implemented.");

  return MPI_DATATYPE_NULL;
}

static NcmMPIDatatype
_ncm_mpi_job_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  g_error ("_ncm_mpi_job_return_datatype: method not implemented.");

  return MPI_DATATYPE_NULL;
}

static gpointer
_ncm_mpi_job_create_input (NcmMPIJob *mpi_job)
{
  g_error ("_ncm_mpi_job_create_input: method not implemented.");

  return NULL;
}

static gpointer
_ncm_mpi_job_create_return (NcmMPIJob *mpi_job)
{
  g_error ("_ncm_mpi_job_create_return: method not implemented.");

  return NULL;
}

static void
_ncm_mpi_job_destroy_input (NcmMPIJob *mpi_job, gpointer input)
{
  g_error ("_ncm_mpi_job_destroy_input: method not implemented.");
}

static void
_ncm_mpi_job_destroy_return (NcmMPIJob *mpi_job, gpointer ret)
{
  g_error ("_ncm_mpi_job_destroy_return: method not implemented.");
}

static gpointer _ncm_mpi_job_get_input_buffer (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf);
static void _ncm_mpi_job_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf);

static gpointer
_ncm_mpi_job_pack_input (NcmMPIJob *mpi_job, gpointer input)
{
  g_error ("_ncm_mpi_job_pack_input: method not implemented.");

  return NULL;
}

static gpointer
_ncm_mpi_job_pack_return (NcmMPIJob *mpi_job, gpointer ret)
{
  g_error ("_ncm_mpi_job_pack_return: method not implemented.");

  return NULL;
}

static void
_ncm_mpi_job_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input)
{
  g_error ("_ncm_mpi_job_unpack_input: method not implemented.");
}

static void
_ncm_mpi_job_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret)
{
  g_error ("_ncm_mpi_job_unpack_return: method not implemented.");
}

static void
_ncm_mpi_job_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret)
{
  g_error ("_ncm_mpi_job_run: method not implemented.");
}

static void
ncm_mpi_job_class_init (NcmMPIJobClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_mpi_job_set_property;
  object_class->get_property = &_ncm_mpi_job_get_property;
  object_class->constructed  = &_ncm_mpi_job_constructed;
  object_class->dispose      = &_ncm_mpi_job_dispose;
  object_class->finalize     = &_ncm_mpi_job_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PLACEHOLDER,
                                   g_param_spec_uint ("placeholder",
                                                      NULL,
                                                      "placeholder",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->work_init  = NULL;
  klass->work_clear = NULL;

  klass->input_datatype  = &_ncm_mpi_job_input_datatype;
  klass->return_datatype = &_ncm_mpi_job_return_datatype;

  klass->create_input  = &_ncm_mpi_job_create_input;
  klass->create_return = &_ncm_mpi_job_create_return;

  klass->destroy_input  = &_ncm_mpi_job_destroy_input;
  klass->destroy_return = &_ncm_mpi_job_destroy_return;

  klass->get_input_buffer  = &_ncm_mpi_job_get_input_buffer;
  klass->get_return_buffer = &_ncm_mpi_job_get_return_buffer;

  klass->destroy_input_buffer  = &_ncm_mpi_job_destroy_input_buffer;
  klass->destroy_return_buffer = &_ncm_mpi_job_destroy_return_buffer;

  klass->pack_input  = &_ncm_mpi_job_pack_input;
  klass->pack_return = &_ncm_mpi_job_pack_return;

  klass->unpack_input  = &_ncm_mpi_job_unpack_input;
  klass->unpack_return = &_ncm_mpi_job_unpack_return;

  klass->run = &_ncm_mpi_job_run;
}

static void
_ncm_mpi_job_work_clear (NcmMPIJob *mpi_job)
{
  /*NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);*/
}

static gpointer
_ncm_mpi_job_get_input_buffer (NcmMPIJob *mpi_job, gpointer input)
{
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
  gpointer *buf_ptr             = ncm_memory_pool_get (self->input_buf_pool);

  g_hash_table_insert (self->input_buf_table, *buf_ptr, buf_ptr);

  return *buf_ptr;
}

static gpointer
_ncm_mpi_job_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret)
{
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
  gpointer *buf_ptr             = ncm_memory_pool_get (self->return_buf_pool);

  g_hash_table_insert (self->return_buf_table, *buf_ptr, buf_ptr);

  return *buf_ptr;
}

static void
_ncm_mpi_job_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf)
{
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
  gpointer *buf_ptr             = g_hash_table_lookup (self->input_buf_table, buf);

  g_assert (buf_ptr != NULL);
  ncm_memory_pool_return (buf_ptr);
}

static void
_ncm_mpi_job_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf)
{
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
  gpointer *buf_ptr             = g_hash_table_lookup (self->return_buf_table, buf);

  g_assert (buf_ptr != NULL);
  ncm_memory_pool_return (buf_ptr);
}

/**
 * ncm_mpi_job_ref:
 * @mpi_job: a #NcmMPIJob
 *
 * Increase the reference of @mpi_job by one.
 *
 * Returns: (transfer full): @mpi_job.
 */
NcmMPIJob *
ncm_mpi_job_ref (NcmMPIJob *mpi_job)
{
  return g_object_ref (mpi_job);
}

/**
 * ncm_mpi_job_free:
 * @mpi_job: a #NcmMPIJob
 *
 * Decrease the reference count of @mpi_job by one.
 *
 */
void
ncm_mpi_job_free (NcmMPIJob *mpi_job)
{
  g_object_unref (mpi_job);
}

/**
 * ncm_mpi_job_clear:
 * @mpi_job: a #NcmMPIJob
 *
 * Decrease the reference count of @mpi_job by one, and sets the pointer *@mpi_job to
 * NULL.
 *
 */
void
ncm_mpi_job_clear (NcmMPIJob **mpi_job)
{
  g_clear_object (mpi_job);
}

/**
 * ncm_mpi_job_work_init: (virtual work_init)
 * @mpi_job: a #NcmMPIJob
 *
 * Method called after @mpi_job is initialized at the slave
 * and before start working.
 *
 */
void
ncm_mpi_job_work_init (NcmMPIJob *mpi_job)
{
  if (NCM_MPI_JOB_GET_CLASS (mpi_job)->work_init != NULL)
    NCM_MPI_JOB_GET_CLASS (mpi_job)->work_init (mpi_job);
}

/**
 * ncm_mpi_job_work_clear: (virtual work_clear)
 * @mpi_job: a #NcmMPIJob
 *
 * Method called during the working phase of @mpi_job and in the
 * end before object destruction. This method can be called multiple
 * times during the work phase.
 *
 */
void
ncm_mpi_job_work_clear (NcmMPIJob *mpi_job)
{
  if (NCM_MPI_JOB_GET_CLASS (mpi_job)->work_clear != NULL)
    NCM_MPI_JOB_GET_CLASS (mpi_job)->work_clear (mpi_job);

  _ncm_mpi_job_work_clear (mpi_job);
}

/**
 * ncm_mpi_job_input_datatype: (virtual input_datatype)
 * @mpi_job: a #NcmMPIJob
 * @len: (out): input length
 * @size: (out): input buffer size
 *
 * Computes the size and datatype of the input buffer.
 *
 * Returns: (transfer none): the input datatype.
 */
NcmMPIDatatype
ncm_mpi_job_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->input_datatype (mpi_job, len, size);
}

/**
 * ncm_mpi_job_return_datatype: (virtual return_datatype)
 * @mpi_job: a #NcmMPIJob
 * @len: (out): input length
 * @size: (out): input buffer size
 *
 * Computes the size and datatype of the return buffer.
 *
 * Returns: (transfer none): the return datatype.
 */
NcmMPIDatatype
ncm_mpi_job_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->return_datatype (mpi_job, len, size);
}

/**
 * ncm_mpi_job_create_input: (virtual create_input)
 * @mpi_job: a #NcmMPIJob
 *
 * Creates a new input object.
 *
 * Returns: (transfer none): the newly input object.
 */
gpointer
ncm_mpi_job_create_input (NcmMPIJob *mpi_job)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->create_input (mpi_job);
}

/**
 * ncm_mpi_job_create_return: (virtual create_return)
 * @mpi_job: a #NcmMPIJob
 *
 * Creates a new return object.
 *
 * Returns: (transfer none): the newly return object.
 */
gpointer
ncm_mpi_job_create_return (NcmMPIJob *mpi_job)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->create_return (mpi_job);
}

/**
 * ncm_mpi_job_destroy_input: (virtual destroy_input)
 * @mpi_job: a #NcmMPIJob
 * @input: an input object
 *
 * Destroy the @input object created with ncm_mpi_job_create_input().
 *
 */
void
ncm_mpi_job_destroy_input (NcmMPIJob *mpi_job, gpointer input)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->destroy_input (mpi_job, input);
}

/**
 * ncm_mpi_job_destroy_return: (virtual destroy_return)
 * @mpi_job: a #NcmMPIJob
 * @ret: a return object
 *
 * Destroy the @return object created with ncm_mpi_job_create_return().
 *
 */
void
ncm_mpi_job_destroy_return (NcmMPIJob *mpi_job, gpointer ret)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->destroy_return (mpi_job, ret);
}

/**
 * ncm_mpi_job_get_input_buffer: (virtual get_input_buffer)
 * @mpi_job: a #NcmMPIJob
 * @input: an input object
 *
 * Creates a buffer from @input compatible with ncm_mpi_job_input_datatype().
 *
 * Returns: (transfer none): the created buffer.
 */
gpointer
ncm_mpi_job_get_input_buffer (NcmMPIJob *mpi_job, gpointer input)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->get_input_buffer (mpi_job, input);
}

/**
 * ncm_mpi_job_get_return_buffer: (virtual get_return_buffer)
 * @mpi_job: a #NcmMPIJob
 * @ret: a return object
 *
 * Creates a buffer from @ret compatible with ncm_mpi_job_return_datatype().
 *
 * Returns: (transfer none): the created buffer.
 */
gpointer
ncm_mpi_job_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->get_return_buffer (mpi_job, ret);
}

/**
 * ncm_mpi_job_destroy_input_buffer: (virtual destroy_input_buffer)
 * @mpi_job: a #NcmMPIJob
 * @input: an input object
 * @buf: a input buffer
 *
 * Destroy @buf created with ncm_mpi_job_get_input_buffer()
 * or ncm_mpi_job_pack_input().
 *
 */
void
ncm_mpi_job_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->destroy_input_buffer (mpi_job, input, buf);
}

/**
 * ncm_mpi_job_destroy_return_buffer: (virtual destroy_return_buffer)
 * @mpi_job: a #NcmMPIJob
 * @ret: a return object
 * @buf: a return buffer
 *
 * Destroy @buf created with ncm_mpi_job_get_return_buffer()
 * or ncm_mpi_job_pack_return().
 *
 */
void
ncm_mpi_job_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->destroy_return_buffer (mpi_job, ret, buf);
}

/**
 * ncm_mpi_job_pack_input: (virtual pack_input)
 * @mpi_job: a #NcmMPIJob
 * @input: the input pointer
 *
 * Packs (when necessary) the input into the input buffer.
 *
 * Returns: (transfer none): the packed buffer.
 */
gpointer
ncm_mpi_job_pack_input (NcmMPIJob *mpi_job, gpointer input)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->pack_input (mpi_job, input);
}

/**
 * ncm_mpi_job_pack_return: (virtual pack_return)
 * @mpi_job: a #NcmMPIJob
 * @ret: the return pointer
 *
 * Packs (when necessary) the return into the return buffer @buf.
 *
 * Returns: (transfer none): the packed buffer.
 */
gpointer
ncm_mpi_job_pack_return (NcmMPIJob *mpi_job, gpointer ret)
{
  return NCM_MPI_JOB_GET_CLASS (mpi_job)->pack_return (mpi_job, ret);
}

/**
 * ncm_mpi_job_unpack_input: (virtual unpack_input)
 * @mpi_job: a #NcmMPIJob
 * @buf: the received buffer
 * @input: the unpacked buffer
 *
 * Unpacks (when necessary) the buffer @buf into the input pointer @input.
 *
 */
void
ncm_mpi_job_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->unpack_input (mpi_job, buf, input);
}

/**
 * ncm_mpi_job_unpack_return: (virtual unpack_return)
 * @mpi_job: a #NcmMPIJob
 * @buf: the received buffer
 * @ret: the unpacked buffer
 *
 * Unpacks (when necessary) the buffer @buf into the return pointer @return.
 *
 */
void
ncm_mpi_job_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->unpack_return (mpi_job, buf, ret);
}

/**
 * ncm_mpi_job_run: (virtual run)
 * @mpi_job: a #NcmMPIJob
 * @input: an input pointer
 * @ret: an return pointer
 *
 * Runs job @mpi_job using @input and returns in @ret.
 *
 */
void
ncm_mpi_job_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret)
{
  NCM_MPI_JOB_GET_CLASS (mpi_job)->run (mpi_job, input, ret);
}

/**
 * ncm_mpi_job_init_all_slaves:
 * @mpi_job: a #NcmMPIJob
 * @ser: a #NcmSerialize
 *
 * Initialize all available slaves with @mpi_job.
 *
 */
void
ncm_mpi_job_init_all_slaves (NcmMPIJob *mpi_job, NcmSerialize *ser)
{
#ifdef HAVE_MPI
  g_assert_cmpint (_mpi_ctrl.rank, ==, NCM_MPI_CTRL_MASTER_ID);

  if (_mpi_ctrl.size > 1)
  {
    NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
    GVariant *job_ser             = ncm_serialize_to_variant (ser, G_OBJECT (mpi_job));
    gconstpointer job_data        = g_variant_get_data (job_ser);
    gint length                   = g_variant_get_size (job_ser);
    MPI_Request *request          = g_new (MPI_Request, _mpi_ctrl.nslaves);
    MPI_Status *status            = g_new (MPI_Status, _mpi_ctrl.nslaves);
    gint *cmds                    = g_new (gint, _mpi_ctrl.nslaves);

    gint i;

    for (i = 0; i < _mpi_ctrl.nslaves; i++)
    {
      gint slave_id = i + 1;

      cmds[i] = NCM_MPI_CTRL_SLAVE_INIT;
      MPI_Isend (&cmds[i], 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &request[i]);
    }

    MPI_Waitall (_mpi_ctrl.nslaves, request, status);

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] All commands sent!\n", _mpi_ctrl.size, _mpi_ctrl.rank);

    for (i = 0; i < _mpi_ctrl.nslaves; i++)
    {
      gint slave_id = i + 1;

      MPI_Isend (job_data, length, MPI_BYTE, slave_id, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD, &request[i]);
    }

    MPI_Waitall (_mpi_ctrl.nslaves, request, status);

    g_free (request);
    g_free (status);
    g_free (cmds);

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] All objects sent!\n", _mpi_ctrl.size, _mpi_ctrl.rank);

    self->owned_slaves       = _mpi_ctrl.nslaves;
    _mpi_ctrl.working_slaves = _mpi_ctrl.nslaves;

    g_variant_unref (job_ser);
  }

#else
  g_error ("ncm_mpi_job_init_all_slaves: MPI unsupported.");
#endif /* HAVE_MPI */
}

/**
 * ncm_mpi_job_run_array:
 * @mpi_job: a #NcmMPIJob
 * @input_array: (array) (element-type GObject): an array of input pointers
 * @ret_array: (array) (element-type GObject): an array of (allocated) return pointers
 *
 * Send work to all slaves in a round-robin fashion. Both arrays @input_array and @ret_array
 * must have the same length and should be filled with the appropriated pointers.
 *
 */
void
ncm_mpi_job_run_array (NcmMPIJob *mpi_job, GPtrArray *input_array, GPtrArray *ret_array)
{
#ifdef HAVE_MPI
  g_assert_cmpint (_mpi_ctrl.rank, ==, NCM_MPI_CTRL_MASTER_ID);

  enum buf_type
  {
    input_type,
    ret_type,
    msg_type,
  };
  struct buf_desc
  {
    gpointer obj;
    gpointer buf;
    enum buf_type t;
  };

  if (_mpi_ctrl.size > 1)
  {
    NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
    const guint njobs             = input_array->len;
    const guint njobs_pa          = njobs - njobs / _mpi_ctrl.size;
    const guint prealloc          = (njobs < 100) ? njobs : 100;
    gint check_block              = 10;
    GArray *req_array             = g_array_sized_new (FALSE, TRUE, sizeof (MPI_Request), prealloc);
    GPtrArray *cmd_array          = g_ptr_array_new_with_free_func (g_free);
    GArray *buf_desc_a            = g_array_new (FALSE, TRUE, sizeof (struct buf_desc));
    guint i;

    g_assert_cmpuint (input_array->len, ==, ret_array->len);

    for (i = 0; i < njobs; i++)
    {
      gpointer input = g_ptr_array_index (input_array, i);
      gpointer ret   = g_ptr_array_index (ret_array, i);

      if (i < njobs_pa)
      {
        const gint slave_id = (i % _mpi_ctrl.nslaves) + 1;
        gint *cmd           = g_new (gint, 1);
        MPI_Request request;

        NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Sending job %d to slave %d.\n", _mpi_ctrl.size, _mpi_ctrl.rank, i, slave_id);

        *cmd = NCM_MPI_CTRL_SLAVE_WORK;
        g_ptr_array_add (cmd_array, cmd);

        {
          struct buf_desc bd = {NULL, NULL, msg_type};

          MPI_Isend (cmd, 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &request);
          g_array_append_val (req_array, request);
          g_array_append_val (buf_desc_a, bd);
        }

        {
          struct buf_desc bd = {input, ncm_mpi_job_pack_input (mpi_job, input), input_type};

          MPI_Isend (bd.buf, self->input_len, self->input_dtype, slave_id, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD, &request);
          g_array_append_val (req_array, request);
          g_array_append_val (buf_desc_a, bd);
        }

        {
          struct buf_desc bd = {ret, ncm_mpi_job_get_return_buffer (mpi_job, ret), ret_type};

          MPI_Irecv (bd.buf, self->return_len, self->input_dtype, slave_id, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, &request);
          g_array_append_val (req_array, request);
          g_array_append_val (buf_desc_a, bd);
        }
      }
      else
      {
        NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Running job %d on master.\n", _mpi_ctrl.size, _mpi_ctrl.rank, i);
        ncm_mpi_job_run (mpi_job, input, ret);
        check_block = 1;
      }

      if ((i > 0) && ((i % check_block) == 0))
      {
        gint j;

        NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Testing %d sends:\n", _mpi_ctrl.size, _mpi_ctrl.rank, req_array->len);

        for (j = (gint) req_array->len - 1; j >= 0; j--)
        {
          MPI_Status status;
          gint done = 0;

          MPI_Test (&g_array_index (req_array, MPI_Request, j), &done, &status);
          NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Send %d is %s!\n", _mpi_ctrl.size, _mpi_ctrl.rank, j, done ? "done" : "not done");

          if (done)
          {
            struct buf_desc bd = g_array_index (buf_desc_a, struct buf_desc, j);

            NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Send %d removing %d!\n", _mpi_ctrl.size, _mpi_ctrl.rank, j, bd.t);

            g_array_remove_index_fast (req_array,  j);
            g_array_remove_index_fast (buf_desc_a, j);

            switch (bd.t)
            {
              case input_type:
                ncm_mpi_job_destroy_input_buffer (mpi_job, bd.obj, bd.buf);
                break;
              case ret_type:
                ncm_mpi_job_unpack_return (mpi_job, bd.buf, bd.obj);
                ncm_mpi_job_destroy_return_buffer (mpi_job, bd.obj, bd.buf);
                break;
              case msg_type:
                break;
              default:
                g_assert_not_reached ();
                break;
            }

            NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Send %d removing %d, finished!\n", _mpi_ctrl.size, _mpi_ctrl.rank, j, bd.t);
          }
        }
      }
    }

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Waiting for the remaining return (%u)!\n", _mpi_ctrl.size, _mpi_ctrl.rank, req_array->len);
    MPI_Waitall (req_array->len, (MPI_Request *) req_array->data, MPI_STATUSES_IGNORE);

    for (i = 0; i < buf_desc_a->len; i++)
    {
      struct buf_desc bd = g_array_index (buf_desc_a, struct buf_desc, i);

      switch (bd.t)
      {
        case input_type:
          ncm_mpi_job_destroy_input_buffer (mpi_job, bd.obj, bd.buf);
          break;
        case ret_type:
          ncm_mpi_job_unpack_return (mpi_job, bd.buf, bd.obj);
          ncm_mpi_job_destroy_return_buffer (mpi_job, bd.obj, bd.buf);
          break;
        case msg_type:
          break;
        default:
          g_assert_not_reached ();
      }
    }

    g_array_unref (req_array);
    g_array_unref (buf_desc_a);
    g_ptr_array_unref (cmd_array);

    return;
  }
  else
  {
    const guint njobs = input_array->len;
    guint i;

    for (i = 0; i < njobs; i++)
    {
      gpointer input = g_ptr_array_index (input_array, i);
      gpointer ret   = g_ptr_array_index (ret_array, i);

      ncm_mpi_job_run (mpi_job, input, ret);
    }

    return;
  }

#else
  g_error ("ncm_mpi_job_run_array: MPI unsupported.");

  return;

#endif /* HAVE_MPI */
}

typedef struct _NcmMPIJobCtrlData
{
  NcmMPIJob *mpi_job;
  GAsyncQueue *jobs;
} NcmMPIJobCtrlData;

typedef struct _NcmMPIJobJob
{
  gpointer input;
  gpointer ret;
  gint id;
} NcmMPIJobJob;

enum buf_type
{
  input_type,
  ret_type,
  msg_type,
};

struct buf_desc
{
  gpointer obj;
  gpointer buf;
  enum buf_type t;
  gdouble t0;
};

#ifdef HAVE_MPI

static gpointer
_ncm_mpi_job_run_array_async_ctrl_thread (gpointer data)
{
  NcmMPIJobCtrlData *ctrl_data  = data;
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (ctrl_data->mpi_job);
  GAsyncQueue *jobs             = ctrl_data->jobs;
  const guint prealloc          = _mpi_ctrl.nslaves;
  GArray *send_req_array        = g_array_sized_new (FALSE, TRUE, sizeof (MPI_Request), prealloc);
  GArray *ret_req_array         = g_array_sized_new (FALSE, TRUE, sizeof (MPI_Request), prealloc);
  GPtrArray *cmd_array          = g_ptr_array_new_with_free_func (g_free);
  GArray *send_buf_desc_a       = g_array_new (FALSE, TRUE, sizeof (struct buf_desc));
  GArray *ret_buf_desc_a        = g_array_new (FALSE, TRUE, sizeof (struct buf_desc));
  gint i;

  for (i = 1; i <= _mpi_ctrl.nslaves; i++)
  {
    NcmMPIJobJob *j = g_async_queue_try_pop (jobs);

    if (j == NULL)
    {
      break;
    }
    else
    {
      gint *cmd = g_new (gint, 1);
      MPI_Request request;

#ifdef NCM_MPI_DEBUG
      const gdouble t0 = ncm_timer_elapsed (self->nt);
#endif /* NCM_MPI_DEBUG */

      NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Sending job %d to slave %d.\n", _mpi_ctrl.size, _mpi_ctrl.rank, j->id, i);

      *cmd = NCM_MPI_CTRL_SLAVE_WORK;
      g_ptr_array_add (cmd_array, cmd);

      {
        struct buf_desc bd = {NULL, NULL, msg_type, 0.0};

        MPI_Isend (cmd, 1, MPI_INT, i, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &request);
        g_array_append_val (send_req_array, request);
        g_array_append_val (send_buf_desc_a, bd);
      }

      {
        struct buf_desc bd = {j->input, ncm_mpi_job_pack_input (ctrl_data->mpi_job, j->input), input_type, 0.0};

        MPI_Isend (bd.buf, self->input_len, self->input_dtype, i, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD, &request);
        g_array_append_val (send_req_array, request);
        g_array_append_val (send_buf_desc_a, bd);
      }

      {
        struct buf_desc bd = {j->ret, ncm_mpi_job_get_return_buffer (ctrl_data->mpi_job, j->ret), ret_type, 0.0};
#ifdef NCM_MPI_DEBUG
        bd.t0 = t0;
#endif /* NCM_MPI_DEBUG */

        MPI_Irecv (bd.buf, self->return_len, self->input_dtype, i, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, &request);
        g_array_append_val (ret_req_array, request);
        g_array_append_val (ret_buf_desc_a, bd);
      }

      g_free (j);
    }
  }

  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Waiting for the first work messages to arrive (%u)!\n", _mpi_ctrl.size, _mpi_ctrl.rank, send_req_array->len);
  MPI_Waitall (send_req_array->len, (MPI_Request *) send_req_array->data, MPI_STATUSES_IGNORE);
  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] All first work messages received (%u)!\n", _mpi_ctrl.size, _mpi_ctrl.rank, send_req_array->len);

  for (i = 0; i < (gint) send_buf_desc_a->len; i++)
  {
    struct buf_desc bd = g_array_index (send_buf_desc_a, struct buf_desc, i);

    switch (bd.t)
    {
      case input_type:
        ncm_mpi_job_destroy_input_buffer (ctrl_data->mpi_job, bd.obj, bd.buf);
        break;
      case msg_type:
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

  g_clear_pointer (&send_req_array,  g_array_unref);
  g_clear_pointer (&send_buf_desc_a, g_array_unref);
  g_clear_pointer (&cmd_array,       g_ptr_array_unref);

  {
    NcmMPIJobJob *j;

    while ((j = g_async_queue_try_pop (jobs)))
    {
      struct buf_desc *bd = NULL;
      MPI_Status status;
      gint index, slave_id, cmd;
      gpointer input_buf;
      MPI_Request *request;

      while (TRUE)
      {
        gint flag = 0;

        MPI_Testany (ret_req_array->len, (MPI_Request *) ret_req_array->data, &index, &flag, &status);

        if (flag)
          break;
        else
          ncm_util_sleep_ms (10);
      }

      bd = &g_array_index (ret_buf_desc_a, struct buf_desc, index);

      ncm_mpi_job_unpack_return (ctrl_data->mpi_job, bd->buf, bd->obj);
      ncm_mpi_job_destroy_return_buffer (ctrl_data->mpi_job, bd->obj, bd->buf);

      request   = &g_array_index (ret_req_array, MPI_Request, index);
      slave_id  = status.MPI_SOURCE;
      cmd       = NCM_MPI_CTRL_SLAVE_WORK;
      input_buf = ncm_mpi_job_pack_input (ctrl_data->mpi_job, j->input);

      NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave %d has finished its job, sending another (%u)!\n", _mpi_ctrl.size, _mpi_ctrl.rank, slave_id, j->id);
      NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Job %d took %fs.\n", _mpi_ctrl.size, _mpi_ctrl.rank, j->id, ncm_timer_elapsed (self->nt) - bd->t0);

#ifdef NCM_MPI_DEBUG
      bd->t0 = ncm_timer_elapsed (self->nt);
#endif /* NCM_MPI_DEBUG */

      MPI_Send (&cmd, 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD);
      MPI_Send (input_buf, self->input_len, self->input_dtype, slave_id, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD);

      bd->obj = j->ret;
      bd->buf = ncm_mpi_job_get_return_buffer (ctrl_data->mpi_job, j->ret);
      bd->t   = ret_type;

      MPI_Irecv (bd->buf, self->return_len, self->input_dtype, slave_id, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, request);

      g_free (j);
    }
  }

  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Waiting for the slaves to finish last jobs (%u)!\n", _mpi_ctrl.size, _mpi_ctrl.rank, ret_req_array->len);
  MPI_Waitall (ret_req_array->len, (MPI_Request *) ret_req_array->data, MPI_STATUSES_IGNORE);

  for (i = 0; i < (gint) ret_buf_desc_a->len; i++)
  {
    struct buf_desc bd = g_array_index (ret_buf_desc_a, struct buf_desc, i);

    switch (bd.t)
    {
      case ret_type:
        ncm_mpi_job_unpack_return (ctrl_data->mpi_job, bd.buf, bd.obj);
        ncm_mpi_job_destroy_return_buffer (ctrl_data->mpi_job, bd.obj, bd.buf);
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

  g_clear_pointer (&ret_req_array,  g_array_unref);
  g_clear_pointer (&ret_buf_desc_a, g_array_unref);

  return NULL;
}

#endif /* HAVE_MPI */

/**
 * ncm_mpi_job_run_array_async:
 * @mpi_job: a #NcmMPIJob
 * @input_array: (array) (element-type GObject): an array of input pointers
 * @ret_array: (array) (element-type GObject): an array of (allocated) return pointers
 *
 * Send work to all slaves using an additional thread to control the slaves work.
 * The main execution thread runs jobs in parallel while it waits for the slaves to finish.
 * Both arrays @input_array and @ret_array must have the same length and should be
 * filled with the appropriated pointers.
 *
 */
void
ncm_mpi_job_run_array_async (NcmMPIJob *mpi_job, GPtrArray *input_array, GPtrArray *ret_array)
{
#ifdef HAVE_MPI
#ifdef NCM_MPI_DEBUG
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);
#endif /* NCM_MPI_DEBUG */
  g_assert_cmpint (_mpi_ctrl.rank, ==, NCM_MPI_CTRL_MASTER_ID);

  if (_mpi_ctrl.size > 1)
  {
    const guint njobs           = input_array->len;
    GAsyncQueue *jobs           = g_async_queue_new_full (g_free);
    NcmMPIJobCtrlData ctrl_data = {mpi_job, jobs};
    GThread *ctrl_thread;
    gint i;

    g_assert_cmpuint (input_array->len, ==, ret_array->len);

    g_async_queue_lock (jobs);

    for (i = 0; i < (gint) njobs; i++)
    {
      NcmMPIJobJob *j_i = g_new (NcmMPIJobJob, 1);

      j_i->input = g_ptr_array_index (input_array, i);
      j_i->ret   = g_ptr_array_index (ret_array, i);
      j_i->id    = i;

      g_async_queue_push_unlocked (jobs, j_i);
    }

    g_async_queue_unlock (jobs);

    ctrl_thread = g_thread_new ("ctrl_thread", &_ncm_mpi_job_run_array_async_ctrl_thread, &ctrl_data);

    {
      NcmMPIJobJob *j;

      while ((j = g_async_queue_try_pop (jobs)))
      {
#ifdef NCM_MPI_DEBUG
        const gdouble t0 = ncm_timer_elapsed (self->nt);
#endif /* NCM_MPI_DEBUG */
        NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Running job %d on master.\n", _mpi_ctrl.size, _mpi_ctrl.rank, j->id);
        ncm_mpi_job_run (mpi_job, j->input, j->ret);
        NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Job %d took %fs.\n", _mpi_ctrl.size, _mpi_ctrl.rank, j->id,
                                 ncm_timer_elapsed (self->nt) - t0);

        g_free (j);
      }
    }

    g_thread_join (ctrl_thread);
    g_async_queue_unref (jobs);
  }
  else
  {
    const guint njobs = input_array->len;
    gint i;

    for (i = 0; i < (gint) njobs; i++)
    {
      gpointer input = g_ptr_array_index (input_array, i);
      gpointer ret   = g_ptr_array_index (ret_array, i);

      ncm_mpi_job_run (mpi_job, input, ret);
    }

    return;
  }

#else
  g_error ("ncm_mpi_job_run_array: MPI unsupported.");

  return;

#endif /* HAVE_MPI */
}

/**
 * ncm_mpi_job_free_all_slaves:
 * @mpi_job: a #NcmMPIJob
 *
 * Frees all available slaves used by @mpi_job.
 *
 */
void
ncm_mpi_job_free_all_slaves (NcmMPIJob *mpi_job)
{
  NcmMPIJobPrivate * const self = ncm_mpi_job_get_instance_private (mpi_job);

#ifdef HAVE_MPI

  if (_mpi_ctrl.rank != NCM_MPI_CTRL_MASTER_ID)
    return;

  if (self->owned_slaves > 0)
  {
    MPI_Request *request = g_new (MPI_Request, _mpi_ctrl.nslaves);
    MPI_Status *status   = g_new (MPI_Status, _mpi_ctrl.nslaves);
    gint *cmds           = g_new (gint, _mpi_ctrl.nslaves);

    gint i;

    g_assert_cmpuint (self->owned_slaves, ==, _mpi_ctrl.nslaves);

    for (i = 0; i < _mpi_ctrl.nslaves; i++)
    {
      gint slave_id = i + 1;

      cmds[i] = NCM_MPI_CTRL_SLAVE_FREE;
      MPI_Isend (&cmds[i], 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &request[i]);
    }

    MPI_Waitall (_mpi_ctrl.nslaves, request, status);

    g_free (request);
    g_free (status);
    g_free (cmds);

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] All slaves free!\n", _mpi_ctrl.size, _mpi_ctrl.rank);

    _mpi_ctrl.working_slaves = 0;
    self->owned_slaves       = 0;
  }

#endif /* HAVE_MPI */
  g_assert_cmpint (self->owned_slaves, ==, 0);
}


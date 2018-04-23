/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job.c
 *
 *  Sat April 21 18:33:19 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_mpi_job
 * @title: NcmMPIJob
 * @short_description: Abstract class to implement MPI jobs
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mpi_job.h"

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_MPI
#include <mpi.h>
#endif /* HAVE_MPI */
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmMPIJobPrivate
{
	guint placeholder;
	guint ret_vec_size;
	guint in_vec_size;
	GPtrArray *vec_array;
};

enum
{
	PROP_0,
	PROP_PLACEHOLDER,
	PROP_RET_VEC_SIZE,
	PROP_IN_VEC_SIZE,
};

extern NcmMPIJobCtrl _mpi_ctrl;

G_DEFINE_ABSTRACT_TYPE (NcmMPIJob, ncm_mpi_job, G_TYPE_OBJECT);

static void
ncm_mpi_job_init (NcmMPIJob *mpi_job)
{
	NcmMPIJobPrivate * const self = mpi_job->priv = G_TYPE_INSTANCE_GET_PRIVATE (mpi_job, NCM_TYPE_MPI_JOB, NcmMPIJobPrivate);

	self->placeholder  = 0;
	self->ret_vec_size = 0;
	self->in_vec_size  = 0;
	self->vec_array    = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
}

static void
_ncm_mpi_job_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcmMPIJob *mpi_job = NCM_MPI_JOB (object);
	NcmMPIJobPrivate * const self = mpi_job->priv;
	g_return_if_fail (NCM_IS_MPI_JOB (object));

	switch (prop_id)
	{
		case PROP_PLACEHOLDER:
			self->placeholder = g_value_get_uint (value);
			break;
		case PROP_RET_VEC_SIZE:
			self->ret_vec_size = g_value_get_uint (value);
			break;
		case PROP_IN_VEC_SIZE:
			self->in_vec_size = g_value_get_uint (value);
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcmMPIJob *mpi_job = NCM_MPI_JOB (object);
	NcmMPIJobPrivate * const self = mpi_job->priv;
	g_return_if_fail (NCM_IS_MPI_JOB (object));

	switch (prop_id)
	{
		case PROP_PLACEHOLDER:
			g_value_set_uint (value, self->placeholder);
			break;
		case PROP_RET_VEC_SIZE:
			g_value_set_uint (value, self->ret_vec_size);
			break;
		case PROP_IN_VEC_SIZE:
			g_value_set_uint (value, self->in_vec_size);
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_dispose (GObject *object)
{
	NcmMPIJob *mpi_job = NCM_MPI_JOB (object);
	NcmMPIJobPrivate * const self = mpi_job->priv;

	g_clear_pointer (&self->vec_array, g_ptr_array_unref);
	
	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_parent_class)->dispose (object);
}

static void
_ncm_mpi_job_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_parent_class)->finalize (object);
}

static void _ncm_mpi_job_work_clear (NcmMPIJob *mpi_job);
static GObject *_ncm_mpi_job_run (NcmMPIJob *mpi_job, GObject *input) { g_error ("_ncm_mpi_job_run: method not implemented."); return NULL; }
static NcmVector *_ncm_mpi_job_run_vector (NcmMPIJob *mpi_job, NcmVector *input) { g_error ("_ncm_mpi_job_run_vector: method not implemented."); return NULL; }

static void
ncm_mpi_job_class_init (NcmMPIJobClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);

	g_type_class_add_private (klass, sizeof (NcmMPIJobPrivate));

	object_class->set_property = &_ncm_mpi_job_set_property;
	object_class->get_property = &_ncm_mpi_job_get_property;
	object_class->dispose      = &_ncm_mpi_job_dispose;
	object_class->finalize     = &_ncm_mpi_job_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_PLACEHOLDER,
	                                 g_param_spec_uint ("placeholder",
	                                                    NULL,
	                                                    "placeholder",
	                                                    0, G_MAXUINT, 0,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_RET_VEC_SIZE,
	                                 g_param_spec_uint ("ret-vec-size",
	                                                    NULL,
	                                                    "Return vector size",
	                                                    1, G_MAXUINT, 1,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_IN_VEC_SIZE,
	                                 g_param_spec_uint ("in-vec-size",
	                                                    NULL,
	                                                    "Input vector size",
	                                                    1, G_MAXUINT, 1,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	klass->work_init  = NULL;
	klass->work_clear = NULL;
	klass->run        = &_ncm_mpi_job_run;
	klass->run_vector = &_ncm_mpi_job_run_vector;
}

static void 
_ncm_mpi_job_work_clear (NcmMPIJob *mpi_job)
{
	NcmMPIJobPrivate * const self = mpi_job->priv;
	g_ptr_array_set_size (self->vec_array, 0);
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
 * ncm_mpi_job_set_vec_sizes:
 * @mpi_job: a #NcmMPIJob
 * @ret_size: return vector size
 * @in_size: input vector size
 *
 * Sets the size of the input/return #NcmVector. This method should be used by #NcmMPIJob
 * only, do not use it elsewhere!
 * 
 */
void 
ncm_mpi_job_set_vec_sizes (NcmMPIJob *mpi_job, const guint ret_size, const guint in_size)
{
	NcmMPIJobPrivate * const self = mpi_job->priv;
	self->ret_vec_size = ret_size;
	self->in_vec_size  = in_size;
}

/**
 * ncm_mpi_job_get_input_vec_size:
 * @mpi_job: a #NcmMPIJob
 *
 * Returns: the size of the input vector.
 */
guint 
ncm_mpi_job_get_input_vec_size (NcmMPIJob *mpi_job)
{
	NcmMPIJobPrivate * const self = mpi_job->priv;

	return self->in_vec_size;
}

/**
 * ncm_mpi_job_get_return_vec_size:
 * @mpi_job: a #NcmMPIJob
 *
 * Returns: the size of the return vector.
 */
guint 
ncm_mpi_job_get_return_vec_size (NcmMPIJob *mpi_job)
{
	NcmMPIJobPrivate * const self = mpi_job->priv;
	return self->ret_vec_size;
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
 * ncm_mpi_job_run:
 * @mpi_job: a #NcmMPIJob
 * @input: a #GObject
 * 
 * Runs job @mpi_job using @input.
 * 
 * Returns: (transfer none): the output #GObject.
 */ 
GObject *
ncm_mpi_job_run (NcmMPIJob *mpi_job, GObject *input)
{
	return NCM_MPI_JOB_GET_CLASS (mpi_job)->run (mpi_job, input);
}

/**
 * ncm_mpi_job_run_vector:
 * @mpi_job: a #NcmMPIJob
 * @input: a #NcmVector
 * 
 * Runs job @mpi_job using @input.
 * 
 * Returns: (transfer none): the output #NcmVector.
 */ 
NcmVector *
ncm_mpi_job_run_vector (NcmMPIJob *mpi_job, NcmVector *input)
{
	NcmMPIJobPrivate * const self = mpi_job->priv;

	NcmVector *vec     = NCM_MPI_JOB_GET_CLASS (mpi_job)->run_vector (mpi_job, input);
	NcmVector *vec_dup = ncm_vector_dup (vec);

	g_ptr_array_add (self->vec_array, ncm_vector_ref (vec_dup));

	return vec_dup;
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
	{
		GVariant *job_ser      = ncm_serialize_to_variant (ser, G_OBJECT (mpi_job));
		gconstpointer job_data = g_variant_get_data (job_ser);
		gsize length           = g_variant_get_size (job_ser);
		MPI_Request *request   = g_new (MPI_Request, _mpi_ctrl.nslaves);
		MPI_Status *status     = g_new (MPI_Status, _mpi_ctrl.nslaves);
		gint *cmds             = g_new (gint, _mpi_ctrl.nslaves);

		gint i;

		for (i = 0; i < _mpi_ctrl.nslaves; i++)
		{
			gint slave_id = i + 1;
			cmds[i] = NCM_MPI_CTRL_SLAVE_INIT;
			MPI_Isend (&cmds[i], 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &request[i]);	
		}
		MPI_Waitall (_mpi_ctrl.nslaves, request, status);

		/*printf ("#[%3d %3d] All commands sent!\n", _mpi_ctrl.size, _mpi_ctrl.rank);*/
		
		for (i = 0; i < _mpi_ctrl.nslaves; i++)
		{
			gint slave_id = i + 1;
			MPI_Isend (job_data, length, MPI_BYTE, slave_id, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD, &request[i]);	
		}
		MPI_Waitall (_mpi_ctrl.nslaves, request, status);

		g_free (request);
		g_free (status);
		g_free (cmds);

		/*printf ("#[%3d %3d] All objects sent!\n", _mpi_ctrl.size, _mpi_ctrl.rank);*/
	}
#else
	g_error ("ncm_mpi_job_init_all_slaves: MPI unsupported.");
#endif /* HAVE_MPI */
}

/**
 * ncm_mpi_job_run_vectors:
 * @mpi_job: a #NcmMPIJob
 * @input_array: (array) (element-type NcmVector): an array of #NcmVector
 * 
 * Send work to all slaves through vectors.
 * 
 * Returns: (transfer full) (array) (element-type NcmVector): an array of #NcmVector
 */
GPtrArray *
ncm_mpi_job_run_vectors (NcmMPIJob *mpi_job, GPtrArray *input_array)
{
#ifdef HAVE_MPI
	g_assert_cmpint (_mpi_ctrl.rank, ==, NCM_MPI_CTRL_MASTER_ID);
	{
		NcmMPIJobPrivate * const self = mpi_job->priv;
		const guint njobs      = input_array->len;
		const guint prealloc   = (njobs < 100) ? njobs : 100;
		const gint check_block = 10;
		GArray *req_array      = g_array_sized_new (FALSE, TRUE, sizeof (MPI_Request), prealloc);
		GPtrArray *cmd_array   = g_ptr_array_new_with_free_func (g_free);
		GPtrArray *ret_array   = g_ptr_array_new ();
		gint i;

		for (i = 0; i < njobs; i++)
		{
			NcmVector *input    = g_ptr_array_index (input_array, i);
			NcmVector *ret      = ncm_vector_new (self->ret_vec_size);
			const gint slave_id = (i % _mpi_ctrl.nslaves) + 1;
			gint *cmd           = g_new (gint, 1);
			MPI_Request request;
			
			*cmd = NCM_MPI_CTRL_SLAVE_WORK_VECTOR;

			g_ptr_array_add (ret_array, ret);
			g_ptr_array_add (cmd_array, cmd);

			MPI_Isend (cmd, 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &request);
			g_array_append_val (req_array, request);
			
			MPI_Isend (ncm_vector_data (input), ncm_vector_len (input), MPI_DOUBLE, slave_id, NCM_MPI_CTRL_TAG_WORK_INPUT,  MPI_COMM_WORLD, &request);
			g_array_append_val (req_array, request);

			MPI_Irecv (ncm_vector_data (ret),   self->ret_vec_size,     MPI_DOUBLE, slave_id, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, &request);
			g_array_append_val (req_array, request);

			if ((i > 0) && ((i % check_block) == 0))
			{
				gint j;
				/*printf ("#[%3d %3d] Testing %d sends:\n", _mpi_ctrl.size, _mpi_ctrl.rank, req_array->len);*/
				for (j = req_array->len - 1; j >= 0; j--)
				{
					MPI_Status status;
					gint done = 0;

					MPI_Test (&g_array_index (req_array, MPI_Request, j), &done, &status);
					/*printf ("#[%3d %3d] Send %d is %s!\n", _mpi_ctrl.size, _mpi_ctrl.rank, j, done ? "done" : "not done");*/
					
					if (done)
						g_array_remove_index_fast (req_array, j);
				}
			}
		}

		/*printf ("#[%3d %3d] Waiting for the remainers return (%u)!\n", _mpi_ctrl.size, _mpi_ctrl.rank, req_array->len);*/
		MPI_Waitall (req_array->len, (MPI_Request *) req_array->data, MPI_STATUSES_IGNORE);

		g_array_unref (req_array);
		g_ptr_array_unref (cmd_array);

		return ret_array;
	}
#else
	g_error ("ncm_mpi_job_run_vectors: MPI unsupported.");
	return NULL;
#endif /* HAVE_MPI */
}

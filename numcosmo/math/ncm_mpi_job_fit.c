/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job_fit.c
 *
 *  Mon April 23 17:04:31 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_fit.c
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
 * SECTION:ncm_mpi_job_fit
 * @title: NcmMPIJobFit
 * @short_description: MPI job object for running #NcmFit
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm_enum_types.h"
#include "math/ncm_mpi_job_fit.h"

#ifndef HAVE_MPI
#define MPI_DATATYPE_NULL (0)
#define MPI_DOUBLE (0)
#endif /* HAVE_MPI */

struct _NcmMPIJobFitPrivate
{
	NcmFit *fit;
	NcmObjArray *func_oa;
	guint fparam_len;
	guint nfuncs;
};

enum
{
	PROP_0,
	PROP_FIT,
	PROP_FUNC_ARRAY,
	PROP_JOB_TYPE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmMPIJobFit, ncm_mpi_job_fit, NCM_TYPE_MPI_JOB);

static void
ncm_mpi_job_fit_init (NcmMPIJobFit *mjfit)
{
	NcmMPIJobFitPrivate * const self = ncm_mpi_job_fit_get_instance_private (mjfit);

	self->fit            = NULL;
	self->func_oa        = NULL;
	self->fparam_len     = 0;
	self->nfuncs         = 0;
}

static void
_ncm_mpi_job_fit_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (object);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	g_return_if_fail (NCM_IS_MPI_JOB_FIT (object));

	switch (prop_id)
	{
		case PROP_FIT:
			g_assert (self->fit == NULL);
			self->fit        = g_value_dup_object (value);
			self->fparam_len = ncm_mset_fparam_len (self->fit->mset);
			break;
    case PROP_FUNC_ARRAY:
    {
			ncm_obj_array_clear (&self->func_oa);
      self->func_oa = g_value_dup_boxed (value);
			self->nfuncs  = 0;
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
				self->nfuncs = self->func_oa->len;
      }
      break;
    }
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_fit_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (object);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	g_return_if_fail (NCM_IS_MPI_JOB_FIT (object));

	switch (prop_id)
	{
		case PROP_FIT:
			g_value_set_object (value, self->fit);
			break;
    case PROP_FUNC_ARRAY:
      g_value_set_boxed (value, self->func_oa);
      break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_fit_dispose (GObject *object)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (object);
	NcmMPIJobFitPrivate * const self = mjfit->priv;

	ncm_fit_clear (&self->fit);
	ncm_obj_array_clear (&self->func_oa);
	
	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_fit_parent_class)->dispose (object);
}

static void
_ncm_mpi_job_fit_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_fit_parent_class)->finalize (object);
}

static MPI_Datatype _ncm_mpi_job_fit_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);
static MPI_Datatype _ncm_mpi_job_fit_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);

static gpointer _ncm_mpi_job_fit_create_input (NcmMPIJob *mpi_job);
static gpointer _ncm_mpi_job_fit_create_return (NcmMPIJob *mpi_job);

static void _ncm_mpi_job_fit_destroy_input (NcmMPIJob *mpi_job, gpointer input);
static void _ncm_mpi_job_fit_destroy_return (NcmMPIJob *mpi_job, gpointer ret);

static gpointer _ncm_mpi_job_fit_get_input_buffer (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_fit_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_fit_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf);
static void _ncm_mpi_job_fit_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf);

static gpointer _ncm_mpi_job_fit_pack_input (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_fit_pack_return (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_fit_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input);
static void _ncm_mpi_job_fit_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret);

static void _ncm_mpi_job_fit_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret);

static void
ncm_mpi_job_fit_class_init (NcmMPIJobFitClass *klass)
{
	GObjectClass* object_class    = G_OBJECT_CLASS (klass);
	NcmMPIJobClass *mpi_job_class = NCM_MPI_JOB_CLASS (klass);

	object_class->set_property = &_ncm_mpi_job_fit_set_property;
	object_class->get_property = &_ncm_mpi_job_fit_get_property;
	object_class->dispose      = &_ncm_mpi_job_fit_dispose;
	object_class->finalize     = &_ncm_mpi_job_fit_finalize;

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
	
	mpi_job_class->input_datatype        = &_ncm_mpi_job_fit_input_datatype;
	mpi_job_class->return_datatype       = &_ncm_mpi_job_fit_return_datatype;
	
	mpi_job_class->create_input          = &_ncm_mpi_job_fit_create_input;
	mpi_job_class->create_return         = &_ncm_mpi_job_fit_create_return;

	mpi_job_class->destroy_input         = &_ncm_mpi_job_fit_destroy_input;
	mpi_job_class->destroy_return        = &_ncm_mpi_job_fit_destroy_return;

	mpi_job_class->get_input_buffer      = &_ncm_mpi_job_fit_get_input_buffer;
	mpi_job_class->get_return_buffer     = &_ncm_mpi_job_fit_get_return_buffer;
	
	mpi_job_class->destroy_input_buffer  = &_ncm_mpi_job_fit_destroy_input_buffer;
	mpi_job_class->destroy_return_buffer = &_ncm_mpi_job_fit_destroy_return_buffer;
	
	mpi_job_class->pack_input            = &_ncm_mpi_job_fit_pack_input;
	mpi_job_class->pack_return           = &_ncm_mpi_job_fit_pack_return;
	
	mpi_job_class->unpack_input          = &_ncm_mpi_job_fit_unpack_input;
	mpi_job_class->unpack_return         = &_ncm_mpi_job_fit_unpack_return;

	mpi_job_class->run                   = &_ncm_mpi_job_fit_run;
}

static MPI_Datatype 
_ncm_mpi_job_fit_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (mpi_job);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	len[0]  = self->fparam_len;
	size[0] = sizeof (gdouble) * len[0];
	return MPI_DOUBLE;
}

static MPI_Datatype 
_ncm_mpi_job_fit_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (mpi_job);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	
	if (self->func_oa == NULL)
	{
		len[0]  = 1 + self->fparam_len;
		size[0] = sizeof (gdouble);
		return MPI_DOUBLE;
	}
	else
	{
		len[0]  = 1 + self->fparam_len + self->nfuncs;
		size[0] = sizeof (gdouble) * len[0];
		return MPI_DOUBLE;
	}
}

static gpointer 
_ncm_mpi_job_fit_create_input (NcmMPIJob *mpi_job)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (mpi_job);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	return ncm_vector_new (self->fparam_len);
}

static gpointer 
_ncm_mpi_job_fit_create_return (NcmMPIJob *mpi_job)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (mpi_job);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	return ncm_vector_new (1 + self->fparam_len + self->nfuncs);
}

static void 
_ncm_mpi_job_fit_destroy_input (NcmMPIJob *mpi_job, gpointer input)
{
	ncm_vector_free (input);
}

static void 
_ncm_mpi_job_fit_destroy_return (NcmMPIJob *mpi_job, gpointer ret)
{
	ncm_vector_free (ret);
}

static gpointer 
_ncm_mpi_job_fit_get_input_buffer (NcmMPIJob *mpi_job, gpointer input)
{
	return ncm_vector_data (input);
}

static gpointer 
_ncm_mpi_job_fit_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret)
{
	return ncm_vector_data (ret);
}

static void 
_ncm_mpi_job_fit_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (input)), ==, GPOINTER_TO_INT (buf));
}

static void 
_ncm_mpi_job_fit_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (ret)), ==, GPOINTER_TO_INT (buf));
}

static gpointer 
_ncm_mpi_job_fit_pack_input (NcmMPIJob *mpi_job, gpointer input)
{
	return ncm_vector_data (input);
}

static gpointer 
_ncm_mpi_job_fit_pack_return (NcmMPIJob *mpi_job, gpointer ret)
{
	return ncm_vector_data (ret);
}

static void 
_ncm_mpi_job_fit_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (input)), ==, GPOINTER_TO_INT (buf));
}

static void 
_ncm_mpi_job_fit_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (ret)), ==, GPOINTER_TO_INT (buf));
}

static void
_ncm_mpi_job_fit_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret)
{
	NcmMPIJobFit *mjt = NCM_MPI_JOB_FIT (mpi_job);
	NcmMPIJobFitPrivate * const self = mjt->priv;

	ncm_fit_params_set_vector (self->fit, input);

	ncm_fit_run (self->fit, NCM_FIT_RUN_MSGS_NONE);
	ncm_fit_m2lnL_val (self->fit, ncm_vector_ptr (ret, 0));
	ncm_mset_fparams_get_vector_offset (self->fit->mset, ret, 1);
	
	if (self->func_oa != NULL)
	{
		gint i;
		for (i = 0; i < self->func_oa->len; i++)
		{
			NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, i));
			const gdouble a_i = ncm_mset_func_eval0 (func, self->fit->mset);

			ncm_vector_set (ret, i + 1 + self->fparam_len, a_i);
		}
	}
}

/**
 * ncm_mpi_job_fit_new:
 * @fit: a #NcmFit
 * @func_oa: (array) (element-type NcmMSetFunc) (allow-none): a #NcmObjArray
 * 
 * Creates a new #NcmMPIJobFit object.
 * 
 * Returns: a new #NcmMPIJobFit.
 */
NcmMPIJobFit *
ncm_mpi_job_fit_new (NcmFit *fit, NcmObjArray *func_oa)
{
	NcmMPIJobFit *mjfit = g_object_new (NCM_TYPE_MPI_JOB_FIT,
	                                    "fit",            fit,
	                                    "function-array", func_oa,
	                                    NULL);
	return mjfit;
}

/**
 * ncm_mpi_job_fit_ref:
 * @mjfit: a #NcmMPIJobFit
 *
 * Increase the reference of @mjfit by one.
 *
 * Returns: (transfer full): @mjfit.
 */
NcmMPIJobFit *
ncm_mpi_job_fit_ref (NcmMPIJobFit *mjfit)
{
  return g_object_ref (mjfit);
}

/**
 * ncm_mpi_job_fit_free:
 * @mjfit: a #NcmMPIJobFit
 *
 * Decrease the reference count of @mjfit by one.
 *
 */
void
ncm_mpi_job_fit_free (NcmMPIJobFit *mjfit)
{
  g_object_unref (mjfit);
}

/**
 * ncm_mpi_job_fit_clear:
 * @mjfit: a #NcmMPIJobFit
 *
 * Decrease the reference count of @mjfit by one, and sets the pointer *@mjfit to
 * NULL.
 *
 */
void
ncm_mpi_job_fit_clear (NcmMPIJobFit **mjfit)
{
  g_clear_object (mjfit);
}

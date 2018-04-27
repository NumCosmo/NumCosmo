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

struct _NcmMPIJobFitPrivate
{
	NcmFit *fit;
  NcmMPIJobFitType job_type;
};

enum
{
	PROP_0,
	PROP_FIT,
	PROP_JOB_TYPE,
};

G_DEFINE_TYPE (NcmMPIJobFit, ncm_mpi_job_fit, NCM_TYPE_MPI_JOB);

static void
ncm_mpi_job_fit_init (NcmMPIJobFit *mjfit)
{
	NcmMPIJobFitPrivate * const self = mjfit->priv = G_TYPE_INSTANCE_GET_PRIVATE (mjfit, NCM_TYPE_MPI_JOB_FIT, NcmMPIJobFitPrivate);

	self->fit      = NULL;
	self->job_type = NCM_MPI_JOB_FIT_TYPE_LEN;

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
			self->fit = g_value_dup_object (value);
			break;
		case PROP_JOB_TYPE:
			g_assert_cmpint (self->job_type, ==, NCM_MPI_JOB_FIT_TYPE_LEN);
			self->job_type = g_value_get_enum (value);
			break;
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
		case PROP_JOB_TYPE:
			g_value_set_enum (value, self->job_type);
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

	g_type_class_add_private (klass, sizeof (NcmMPIJobFitPrivate));

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
	                                 PROP_JOB_TYPE,
	                                 g_param_spec_enum ("job-type",
	                                                    NULL,
	                                                    "MPI Job Fit type",
	                                                    NCM_TYPE_MPI_JOB_FIT_TYPE, NCM_MPI_JOB_FIT_TYPE_M2LNL_VAL,
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
	len[0]  = ncm_mset_fparam_len (self->fit->mset);
	size[0] = sizeof (gdouble) * len[0];
	return MPI_DOUBLE;
}

static MPI_Datatype 
_ncm_mpi_job_fit_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
	len[0]  = 1;
	size[0] = sizeof (gdouble);
	return MPI_DOUBLE;
}

static gpointer 
_ncm_mpi_job_fit_create_input (NcmMPIJob *mpi_job)
{
	NcmMPIJobFit *mjfit = NCM_MPI_JOB_FIT (mpi_job);
	NcmMPIJobFitPrivate * const self = mjfit->priv;
	return ncm_vector_new (ncm_mset_fparam_len (self->fit->mset));
}

static gpointer 
_ncm_mpi_job_fit_create_return (NcmMPIJob *mpi_job)
{
	return ncm_vector_new (1);
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
	
	switch (self->job_type)
	{
		case NCM_MPI_JOB_FIT_TYPE_M2LNL_VAL:
		{
			ncm_fit_params_set_vector (self->fit, input);
			ncm_fit_m2lnL_val (self->fit, ncm_vector_ptr (ret, 0));
			break;
		}
		default:
			g_assert_not_reached ();
			break;
	}
}

/**
 * ncm_mpi_job_fit_new:
 * @fit: a #NcmFit
 * @job_type: a #NcmMPIJobFitType
 * 
 * Creates a new #NcmMPIJobFit object.
 * 
 * Returns: a new #NcmMPIJobFit.
 */
NcmMPIJobFit *
ncm_mpi_job_fit_new (NcmFit *fit, NcmMPIJobFitType job_type)
{
	NcmMPIJobFit *mjfit = g_object_new (NCM_TYPE_MPI_JOB_FIT,
	                                    "fit",      fit,
	                                    "job-type", job_type,
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


/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job_mcmc.c
 *
 *  Fri April 27 16:32:14 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_mcmc.c
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
 * SECTION:ncm_mpi_job_mcmc
 * @title: NcmMPIJobMCMC
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
#include "math/ncm_mpi_job_mcmc.h"

struct _NcmMPIJobMCMCPrivate
{
	NcmFit *fit;
	NcmObjArray *func_oa;
};

enum
{
	PROP_0,
	PROP_FIT,
	PROP_FUNC_ARRAY,
	PROP_JOB_TYPE,
};

G_DEFINE_TYPE (NcmMPIJobMCMC, ncm_mpi_job_mcmc, NCM_TYPE_MPI_JOB);

static void
ncm_mpi_job_mcmc_init (NcmMPIJobMCMC *mjmcmc)
{
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv = G_TYPE_INSTANCE_GET_PRIVATE (mjmcmc, NCM_TYPE_MPI_JOB_MCMC, NcmMPIJobMCMCPrivate);

	self->fit      = NULL;
	self->func_oa  = NULL;

}

static void
_ncm_mpi_job_mcmc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (object);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;
	g_return_if_fail (NCM_IS_MPI_JOB_MCMC (object));

	switch (prop_id)
	{
		case PROP_FIT:
			g_assert (self->fit == NULL);
			self->fit = g_value_dup_object (value);
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
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_mcmc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (object);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;
	g_return_if_fail (NCM_IS_MPI_JOB_MCMC (object));

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
_ncm_mpi_job_mcmc_dispose (GObject *object)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (object);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;

	ncm_fit_clear (&self->fit);
	
	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_mcmc_parent_class)->dispose (object);
}

static void
_ncm_mpi_job_mcmc_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_mcmc_parent_class)->finalize (object);
}

static MPI_Datatype _ncm_mpi_job_mcmc_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);
static MPI_Datatype _ncm_mpi_job_mcmc_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size);

static gpointer _ncm_mpi_job_mcmc_create_input (NcmMPIJob *mpi_job);
static gpointer _ncm_mpi_job_mcmc_create_return (NcmMPIJob *mpi_job);

static void _ncm_mpi_job_mcmc_destroy_input (NcmMPIJob *mpi_job, gpointer input);
static void _ncm_mpi_job_mcmc_destroy_return (NcmMPIJob *mpi_job, gpointer ret);

static gpointer _ncm_mpi_job_mcmc_get_input_buffer (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_mcmc_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_mcmc_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf);
static void _ncm_mpi_job_mcmc_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf);

static gpointer _ncm_mpi_job_mcmc_pack_input (NcmMPIJob *mpi_job, gpointer input);
static gpointer _ncm_mpi_job_mcmc_pack_return (NcmMPIJob *mpi_job, gpointer ret);

static void _ncm_mpi_job_mcmc_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input);
static void _ncm_mpi_job_mcmc_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret);

static void _ncm_mpi_job_mcmc_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret);

static void
ncm_mpi_job_mcmc_class_init (NcmMPIJobMCMCClass *klass)
{
	GObjectClass* object_class    = G_OBJECT_CLASS (klass);
	NcmMPIJobClass *mpi_job_class = NCM_MPI_JOB_CLASS (klass);

	g_type_class_add_private (klass, sizeof (NcmMPIJobMCMCPrivate));

	object_class->set_property = &_ncm_mpi_job_mcmc_set_property;
	object_class->get_property = &_ncm_mpi_job_mcmc_get_property;
	object_class->dispose      = &_ncm_mpi_job_mcmc_dispose;
	object_class->finalize     = &_ncm_mpi_job_mcmc_finalize;

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
	
	mpi_job_class->input_datatype        = &_ncm_mpi_job_mcmc_input_datatype;
	mpi_job_class->return_datatype       = &_ncm_mpi_job_mcmc_return_datatype;
	
	mpi_job_class->create_input          = &_ncm_mpi_job_mcmc_create_input;
	mpi_job_class->create_return         = &_ncm_mpi_job_mcmc_create_return;

	mpi_job_class->destroy_input         = &_ncm_mpi_job_mcmc_destroy_input;
	mpi_job_class->destroy_return        = &_ncm_mpi_job_mcmc_destroy_return;

	mpi_job_class->get_input_buffer      = &_ncm_mpi_job_mcmc_get_input_buffer;
	mpi_job_class->get_return_buffer     = &_ncm_mpi_job_mcmc_get_return_buffer;
	
	mpi_job_class->destroy_input_buffer  = &_ncm_mpi_job_mcmc_destroy_input_buffer;
	mpi_job_class->destroy_return_buffer = &_ncm_mpi_job_mcmc_destroy_return_buffer;
	
	mpi_job_class->pack_input            = &_ncm_mpi_job_mcmc_pack_input;
	mpi_job_class->pack_return           = &_ncm_mpi_job_mcmc_pack_return;
	
	mpi_job_class->unpack_input          = &_ncm_mpi_job_mcmc_unpack_input;
	mpi_job_class->unpack_return         = &_ncm_mpi_job_mcmc_unpack_return;

	mpi_job_class->run                   = &_ncm_mpi_job_mcmc_run;
}

static MPI_Datatype 
_ncm_mpi_job_mcmc_input_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;
	len[0]  = ncm_mset_fparam_len (self->fit->mset);
	size[0] = sizeof (gdouble) * len[0];
	return MPI_DOUBLE;
}

static MPI_Datatype 
_ncm_mpi_job_mcmc_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;
	if (self->func_oa == NULL)
	{
		len[0]  = 1;
		size[0] = sizeof (gdouble);
		return MPI_DOUBLE;
	}
	else
	{
		len[0]  = 1 + self->func_oa->len;
		size[0] = sizeof (gdouble) * len[0];
		return MPI_DOUBLE;
	}
}

static gpointer 
_ncm_mpi_job_mcmc_create_input (NcmMPIJob *mpi_job)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;
	return ncm_vector_new (ncm_mset_fparam_len (self->fit->mset));
}

static gpointer 
_ncm_mpi_job_mcmc_create_return (NcmMPIJob *mpi_job)
{
	return ncm_vector_new (1);
}

static void 
_ncm_mpi_job_mcmc_destroy_input (NcmMPIJob *mpi_job, gpointer input)
{
	ncm_vector_free (input);
}

static void 
_ncm_mpi_job_mcmc_destroy_return (NcmMPIJob *mpi_job, gpointer ret)
{
	ncm_vector_free (ret);
}

static gpointer 
_ncm_mpi_job_mcmc_get_input_buffer (NcmMPIJob *mpi_job, gpointer input)
{
	return ncm_vector_data (input);
}

static gpointer 
_ncm_mpi_job_mcmc_get_return_buffer (NcmMPIJob *mpi_job, gpointer ret)
{
	return ncm_vector_data (ret);
}

static void 
_ncm_mpi_job_mcmc_destroy_input_buffer (NcmMPIJob *mpi_job, gpointer input, gpointer buf)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (input)), ==, GPOINTER_TO_INT (buf));
}

static void 
_ncm_mpi_job_mcmc_destroy_return_buffer (NcmMPIJob *mpi_job, gpointer ret, gpointer buf)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (ret)), ==, GPOINTER_TO_INT (buf));
}

static gpointer 
_ncm_mpi_job_mcmc_pack_input (NcmMPIJob *mpi_job, gpointer input)
{
	return ncm_vector_data (input);
}

static gpointer 
_ncm_mpi_job_mcmc_pack_return (NcmMPIJob *mpi_job, gpointer ret)
{
	return ncm_vector_data (ret);
}

static void 
_ncm_mpi_job_mcmc_unpack_input (NcmMPIJob *mpi_job, gpointer buf, gpointer input)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (input)), ==, GPOINTER_TO_INT (buf));
}

static void 
_ncm_mpi_job_mcmc_unpack_return (NcmMPIJob *mpi_job, gpointer buf, gpointer ret)
{
	g_assert_cmphex (GPOINTER_TO_INT (ncm_vector_data (ret)), ==, GPOINTER_TO_INT (buf));
}

static void
_ncm_mpi_job_mcmc_run (NcmMPIJob *mpi_job, gpointer input, gpointer ret)
{
	NcmMPIJobMCMC *mjt = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjt->priv;

	ncm_fit_params_set_vector (self->fit, input);
	ncm_fit_m2lnL_val (self->fit, ncm_vector_ptr (ret, 0));

	if (self->func_oa != NULL)
	{
		gint i;
		for (i = 0; i < self->func_oa->len; i++)
		{
			NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, i));
			const gdouble a_i = ncm_mset_func_eval0 (func, self->fit->mset);

			ncm_vector_set (ret, i + 1, a_i);
		}
	}
}

/**
 * ncm_mpi_job_mcmc_new:
 * @fit: a #NcmFit
 * @func_oa: (array) (element-type NcmMSetFunc) (allow-none): a #NcmObjArray
 * 
 * Creates a new #NcmMPIJobMCMC object.
 * 
 * Returns: a new #NcmMPIJobMCMC.
 */
NcmMPIJobMCMC *
ncm_mpi_job_mcmc_new (NcmFit *fit, NcmObjArray *func_oa)
{
	NcmMPIJobMCMC *mjmcmc = g_object_new (NCM_TYPE_MPI_JOB_MCMC,
	                                      "fit",            fit,
	                                      "function-array", func_oa,
	                                      NULL);
	return mjmcmc;
}

/**
 * ncm_mpi_job_mcmc_ref:
 * @mjmcmc: a #NcmMPIJobMCMC
 *
 * Increase the reference of @mjmcmc by one.
 *
 * Returns: (transfer full): @mjmcmc.
 */
NcmMPIJobMCMC *
ncm_mpi_job_mcmc_ref (NcmMPIJobMCMC *mjmcmc)
{
  return g_object_ref (mjmcmc);
}

/**
 * ncm_mpi_job_mcmc_free:
 * @mjmcmc: a #NcmMPIJobMCMC
 *
 * Decrease the reference count of @mjmcmc by one.
 *
 */
void
ncm_mpi_job_mcmc_free (NcmMPIJobMCMC *mjmcmc)
{
  g_object_unref (mjmcmc);
}

/**
 * ncm_mpi_job_mcmc_clear:
 * @mjmcmc: a #NcmMPIJobMCMC
 *
 * Decrease the reference count of @mjmcmc by one, and sets the pointer *@mjmcmc to
 * NULL.
 *
 */
void
ncm_mpi_job_mcmc_clear (NcmMPIJobMCMC **mjmcmc)
{
  g_clear_object (mjmcmc);
}


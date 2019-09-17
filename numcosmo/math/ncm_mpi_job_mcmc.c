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
#include "math/ncm_fit_esmcmc.h"

#ifndef HAVE_MPI
#define MPI_DATATYPE_NULL (0)
#define MPI_DOUBLE (0)
#endif /* HAVE_MPI */

struct _NcmMPIJobMCMCPrivate
{
	NcmFit *fit;
	NcmObjArray *func_oa;
	gint fparam_len;
	gint nadd_vals;
};

enum
{
	PROP_0,
	PROP_FIT,
	PROP_FUNC_ARRAY,
	PROP_JOB_TYPE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmMPIJobMCMC, ncm_mpi_job_mcmc, NCM_TYPE_MPI_JOB);

static void
ncm_mpi_job_mcmc_init (NcmMPIJobMCMC *mjmcmc)
{
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv = ncm_mpi_job_mcmc_get_instance_private (mjmcmc);

	self->fit        = NULL;
	self->func_oa    = NULL;

	self->fparam_len = 0;
	self->nadd_vals  = 0;
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
_ncm_mpi_job_mcmc_constructed (GObject *object)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (object);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;

	self->fparam_len = ncm_mset_fparam_len (self->fit->mset);
	self->nadd_vals  = 1 + ((self->func_oa != NULL) ? self->func_oa->len : 0);
	
	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_mcmc_parent_class)->constructed (object);
}

static void
_ncm_mpi_job_mcmc_dispose (GObject *object)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (object);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;

	ncm_fit_clear (&self->fit);
	ncm_obj_array_clear (&self->func_oa);

	
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

	object_class->set_property = &_ncm_mpi_job_mcmc_set_property;
	object_class->get_property = &_ncm_mpi_job_mcmc_get_property;
	object_class->constructed  = &_ncm_mpi_job_mcmc_constructed;
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

	len[0]  = self->fparam_len + NCM_FIT_ESMCMC_MPI_IN_LEN;
	size[0] = sizeof (gdouble) * len[0];

	return MPI_DOUBLE;
}

static MPI_Datatype 
_ncm_mpi_job_mcmc_return_datatype (NcmMPIJob *mpi_job, gint *len, gint *size)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;

	len[0]  = NCM_FIT_ESMCMC_MPI_OUT_LEN + self->nadd_vals;
	size[0] = sizeof (gdouble) * len[0];
	
	return MPI_DOUBLE;
}

static gpointer 
_ncm_mpi_job_mcmc_create_input (NcmMPIJob *mpi_job)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;

	return ncm_vector_new (self->fparam_len + NCM_FIT_ESMCMC_MPI_IN_LEN);
}

static gpointer 
_ncm_mpi_job_mcmc_create_return (NcmMPIJob *mpi_job)
{
	NcmMPIJobMCMC *mjmcmc = NCM_MPI_JOB_MCMC (mpi_job);
	NcmMPIJobMCMCPrivate * const self = mjmcmc->priv;
	
	return ncm_vector_new (NCM_FIT_ESMCMC_MPI_OUT_LEN + self->nadd_vals);
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
	const gdouble m2lnL_cur = ncm_vector_get (input, self->fparam_len + 0);
	const gdouble norm      = ncm_vector_get (input, self->fparam_len + 1);
	const gdouble jump      = ncm_vector_get (input, self->fparam_len + 2);
	gdouble *m2lnL_star     = ncm_vector_ptr (ret, NCM_FIT_ESMCMC_MPI_OUT_LEN);
	gdouble prob            = 0.0;
	gboolean accepted       = FALSE;
	
	ncm_fit_params_set_vector (self->fit, input);
	ncm_fit_m2lnL_val (self->fit, m2lnL_star);

	if (gsl_finite (m2lnL_star[0]))
	{
		if (jump >= 0.0)
		{
			prob     = exp ((m2lnL_cur - m2lnL_star[0]) * 0.5 + norm);
			prob     = MIN (prob, 1.0);
			accepted = (jump < prob);
		}
		else
			accepted = TRUE;
	}

	if (accepted)
	{
		ncm_vector_set (ret, 0, 1.0);
		/*printf ("# has oa?! %p\n", self->func_oa);*/
		if (self->func_oa != NULL)
		{
			gint i;
			for (i = 0; i < self->func_oa->len; i++)
			{
				NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, i));
				const gdouble a_i = ncm_mset_func_eval0 (func, self->fit->mset);
				/*printf ("### %d % 22.15g\n", i, a_i);*/
				ncm_vector_set (ret, i + 1 + NCM_FIT_ESMCMC_MPI_OUT_LEN, a_i);
			}
		}
	}
	else
		ncm_vector_set (ret, 0, 0.0);
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


/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_mpi_job_test.c
 *
 *  Sun April 22 14:47:49 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mpi_job_test.c
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
 * SECTION:ncm_mpi_job_test
 * @title: NcmMPIJobTest
 * @short_description: Test implementation of MPI job class
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mpi_job_test.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <unistd.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmMPIJobTestPrivate
{
	NcmVector *vec;
	GPtrArray *ret_array;
	NcmRNG *rng;
};

enum
{
	PROP_0,
	PROP_VECTOR
};

G_DEFINE_TYPE (NcmMPIJobTest, ncm_mpi_job_test, NCM_TYPE_MPI_JOB);

static void
ncm_mpi_job_test_init (NcmMPIJobTest *mjt)
{
	NcmMPIJobTestPrivate * const self = mjt->priv = G_TYPE_INSTANCE_GET_PRIVATE (mjt, NCM_TYPE_MPI_JOB_TEST, NcmMPIJobTestPrivate);

	self->vec       = NULL;
	self->ret_array = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
	self->rng       = ncm_rng_new (NULL);

	ncm_rng_set_random_seed (self->rng, FALSE);
}

static void
_ncm_mpi_job_test_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcmMPIJobTest *mjt = NCM_MPI_JOB_TEST (object);
	NcmMPIJobTestPrivate * const self = mjt->priv;
	g_return_if_fail (NCM_IS_MPI_JOB_TEST (object));

	switch (prop_id)
	{
		case PROP_VECTOR:
			ncm_vector_clear (&self->vec);
			self->vec = g_value_dup_object (value);
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_test_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcmMPIJobTest *mjt = NCM_MPI_JOB_TEST (object);
	NcmMPIJobTestPrivate * const self = mjt->priv;
	g_return_if_fail (NCM_IS_MPI_JOB_TEST (object));

	switch (prop_id)
	{
		case PROP_VECTOR:
			g_value_set_object (value, self->vec);
			ncm_mpi_job_set_vec_sizes (NCM_MPI_JOB (object), 1, 1);
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_mpi_job_test_dispose (GObject *object)
{
	NcmMPIJobTest *mjt = NCM_MPI_JOB_TEST (object);
	NcmMPIJobTestPrivate * const self = mjt->priv;
	
	ncm_vector_clear (&self->vec);
	ncm_rng_clear (&self->rng);
	g_clear_pointer (&self->ret_array, g_ptr_array_unref);

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_test_parent_class)->dispose (object);
}

static void
_ncm_mpi_job_test_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_mpi_job_test_parent_class)->finalize (object);
}

static GObject *_ncm_mpi_job_test_run (NcmMPIJob *mpi_job, GObject *input);
static NcmVector *_ncm_mpi_job_test_run_vector (NcmMPIJob *mpi_job, NcmVector *input);

static void
ncm_mpi_job_test_class_init (NcmMPIJobTestClass *klass)
{
	GObjectClass* object_class    = G_OBJECT_CLASS (klass);
	NcmMPIJobClass *mpi_job_class = NCM_MPI_JOB_CLASS (klass);

	g_type_class_add_private (klass, sizeof (NcmMPIJobTestPrivate));

	object_class->set_property = &_ncm_mpi_job_test_set_property;
	object_class->get_property = &_ncm_mpi_job_test_get_property;
	object_class->dispose      = &_ncm_mpi_job_test_dispose;
	object_class->finalize     = &_ncm_mpi_job_test_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_VECTOR,
	                                 g_param_spec_object ("vector",
	                                                      NULL,
	                                                      "vector",
	                                                      NCM_TYPE_VECTOR,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	mpi_job_class->run        = &_ncm_mpi_job_test_run;
	mpi_job_class->run_vector = &_ncm_mpi_job_test_run_vector;
}

static GObject *
_ncm_mpi_job_test_run (NcmMPIJob *mpi_job, GObject *input)
{
	return G_OBJECT (_ncm_mpi_job_test_run_vector (mpi_job, NCM_VECTOR (input)));
}

static NcmVector *
_ncm_mpi_job_test_run_vector (NcmMPIJob *mpi_job, NcmVector *input)
{
	g_assert_cmpuint  (ncm_vector_len (input), ==, 1);
	{
		NcmMPIJobTest *mjt = NCM_MPI_JOB_TEST (mpi_job);
		NcmMPIJobTestPrivate * const self = mjt->priv;
		
		NcmVector *ret = ncm_vector_new (1);
		guint index    = floor (ncm_vector_get (input, 0));

		g_assert_cmpuint (index, <, ncm_vector_len (self->vec));

		ncm_vector_set (ret, 0, ncm_vector_get (self->vec, index));
		g_ptr_array_add (self->ret_array, ret);

		/*printf ("# Received %.5u.\n", index);*/
		//sleep (1 + gsl_rng_uniform_int (self->rng->r, 4));
		/*printf ("# Received %.5u done!\n", index);*/

		return ret;
	}
}

/**
 * ncm_mpi_job_test_new:
 * 
 * Creates a new #NcmMPIJobTest object.
 * 
 * Returns: a new #NcmMPIJobTest.
 */
NcmMPIJobTest *
ncm_mpi_job_test_new (void)
{
  NcmMPIJobTest *mjt = g_object_new (NCM_TYPE_MPI_JOB_TEST,
                                     NULL);
  return mjt;
}

/**
 * ncm_mpi_job_test_ref:
 * @mjt: a #NcmMPIJobTest
 *
 * Increase the reference of @mjt by one.
 *
 * Returns: (transfer full): @mjt.
 */
NcmMPIJobTest *
ncm_mpi_job_test_ref (NcmMPIJobTest *mjt)
{
  return g_object_ref (mjt);
}

/**
 * ncm_mpi_job_test_free:
 * @mjt: a #NcmMPIJobTest
 *
 * Decrease the reference count of @mjt by one.
 *
 */
void
ncm_mpi_job_test_free (NcmMPIJobTest *mjt)
{
  g_object_unref (mjt);
}

/**
 * ncm_mpi_job_test_clear:
 * @mjt: a #NcmMPIJobTest
 *
 * Decrease the reference count of @mjt by one, and sets the pointer *@mjt to
 * NULL.
 *
 */
void
ncm_mpi_job_test_clear (NcmMPIJobTest **mjt)
{
  g_clear_object (mjt);
}

/**
 * ncm_mpi_job_test_set_rand_vector:
 * @mjt: a #NcmMPIJobTest
 * @len: vector length
 * @rng: a #NcmRNG
 *
 * Sets a random vector of length @len in @mjt.
 *
 */
void
ncm_mpi_job_test_set_rand_vector (NcmMPIJobTest *mjt, const guint len, NcmRNG *rng)
{
	NcmMPIJobTestPrivate * const self = mjt->priv;
	gint i;

	g_assert_cmpuint (len, >, 0);
	ncm_vector_clear (&self->vec);

	self->vec = ncm_vector_new (len);

	for (i = 0; i < len; i++)
	{
		const gdouble v_i = ncm_rng_gaussian_gen (rng, 0, 1.0);
		ncm_vector_set (self->vec, i, v_i);
	}
}

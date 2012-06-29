/***************************************************************************
 *            data_poisson.c
 *
 *  Sun Apr  4 21:57:39 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:data_poisson
 * @title: Poisson Distribution Data
 * @short_description: #NcData representing data from a Poisson distribution
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gmp.h>

/********************************************************************************************************/

static NcDataPoisson *
_poisson_data_alloc (guint np)
{
  NcDataPoisson *poisson = g_slice_new0 (NcDataPoisson);
	poisson->h = gsl_histogram_alloc (np);
  poisson->extra_data = NULL;

  return poisson;
}

/**********************************************************************************************/

static void
_poisson_data_copy (gpointer dest_ptr, gpointer src_ptr)
{
  NcDataPoisson *dest = (NcDataPoisson *) dest_ptr;
  NcDataPoisson *src = (NcDataPoisson *) src_ptr;

	gsl_histogram_memcpy (dest->h, src->h);

  if (src->extra_data && dest->extra_data)
  {
    if (NC_DATA_STRUCT_HAS_COPY(src->extra_data))
      NC_DATA_STRUCT_COPY (dest->extra_data, src->extra_data);
    else
      g_error ("NcDataStruct extra do not implement copy.");
  }

}

static void
_poisson_data_free (gpointer poisson_ptr)
{
  NcDataPoisson *poisson = (NcDataPoisson *) poisson_ptr;
	gsl_histogram_free (poisson->h);

  if (poisson->extra_data)
    nc_data_struct_free (poisson->extra_data);

  g_slice_free (NcDataPoisson, poisson);
}

static guint _poisson_get_length (gpointer poisson_ptr) { return ((NcDataPoisson *) poisson_ptr)->h->n; }

/********************************************************************************************************/

static void
_poisson_calc_m2lnL (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  NcmMSetFunc *int_func = NCM_MSET_FUNC (model);
  NcDataPoisson *poisson = (NcDataPoisson *) data;
  const gdouble x0 = poisson->h->range[0];
  gdouble lambdai = ncm_mset_func_eval1 (int_func, mset, x0);
  gdouble lambdaip1;
  gdouble lambda;
  gint i;
  *m2lnL = 0.0;

  for (i = 0; i < poisson->h->n; i++)
  {
    const gulong N_i = gsl_histogram_get (poisson->h, i);
    const gdouble x_i = poisson->h->range[i + 1];
    lambdaip1 = ncm_mset_func_eval1 (int_func, mset, x_i);
    lambda = (lambdaip1 - lambdai);
    *m2lnL += -2.0 * ( N_i * log(lambda) - lambda - gsl_sf_lngamma(N_i + 1.0));
    //printf ("Ni = % 20.15ld lambda = % 20.15g Ni*log(lambda) = % 20.15g lngamma(N_i + 1.0) = % 20.15g\n", N_i, lambda, N_i * log(lambda), gsl_sf_lngamma(N_i + 1.0));
    lambdai = lambdaip1;
  }

  return;
}

static void
_poisson_calc_leastsquares_f (NcmMSet *mset, gpointer model, gpointer data, NcmVector *v)
{
  NcmMSetFunc *int_func = NCM_MSET_FUNC (model);
  NcDataPoisson *poisson = (NcDataPoisson *) data;
  const gdouble x0 = poisson->h->range[0];
  gdouble lambdai = ncm_mset_func_eval1 (int_func, mset, x0);
  gdouble lambdaip1;
  gdouble lambda;
  gint i;

  for (i = 0; i < poisson->h->n; i++)
  {
    const gulong N_i = gsl_histogram_get (poisson->h, i);
    const gdouble x_i = poisson->h->range[i + 1];
    lambdaip1 = ncm_mset_func_eval1 (int_func, mset, x_i);
    lambda = (lambdaip1 - lambdai);

		ncm_vector_set (v, i, (lambda - N_i) / sqrt(lambda));
		lambdai = lambdaip1;
  }
}

/********************************************************************************************************/

static gpointer
_poisson_data_dup (gpointer poisson_ptr)
{
  NcDataPoisson *poisson = (NcDataPoisson *) poisson_ptr;
  NcDataPoisson *clone;
  clone = _poisson_data_alloc (poisson->h->n);
  _poisson_data_copy (clone, poisson);

  return clone;
}

/********************************************************************************************************/

static void
_poisson_data_end (gpointer poisson_ptr)
{
  NcDataPoisson *poisson = (NcDataPoisson *) poisson_ptr;

  if (poisson->extra_data)
  {
    if (NC_DATA_STRUCT_HAS_END (poisson->extra_data))
      NC_DATA_STRUCT_END (poisson->extra_data);
  }
}

/********************************************************************************************************/

static void
_poisson_data_begin (gpointer poisson_ptr)
{
  NcDataPoisson *poisson = (NcDataPoisson *) poisson_ptr;

  if (poisson->extra_data)
  {
    if (NC_DATA_STRUCT_HAS_BEGIN (poisson->extra_data))
      NC_DATA_STRUCT_BEGIN (poisson->extra_data);
  }
}

/********************************************************************************************************/

static void
_poisson_resample (NcmMSet *mset, gpointer model, gpointer data)
{
  NcmMSetFunc *int_func = NCM_MSET_FUNC (model);
  NcDataPoisson *poisson = (NcDataPoisson *) data;
  const gdouble x0 = poisson->h->range[0];
  gdouble lambdai = ncm_mset_func_eval1 (int_func, mset, x0);
  gdouble lambdaip1;
  gdouble lambda;
  gint i;
  gsl_rng *rng = ncm_get_rng ();

  for (i = 0; i < poisson->h->n; i++)
  {
    const gdouble x_i = poisson->h->range[i + 1];
    gulong N_i;

    lambdaip1 = ncm_mset_func_eval1 (int_func, mset, x_i);
    lambda = (lambdaip1 - lambdai);
    N_i = gsl_ran_poisson (rng, lambda);
		poisson->h->bin[i] = N_i;

    lambdai = lambdaip1;
  }
}

/********************************************************************************************************/

static NcDataStruct *
_nc_data_struct_poisson_new (void)
{
  NcDataStruct *dts = nc_data_struct_new ();
  dts->data       = NULL;
  dts->dup        = &_poisson_data_dup;
  dts->copy       = &_poisson_data_copy;
  dts->free       = &_poisson_data_free;
  dts->begin      = &_poisson_data_begin;
  dts->end        = &_poisson_data_end;
  dts->get_length = &_poisson_get_length;

  return dts;
}

/********************************************************************************************************/

/**
 * nc_data_poisson_new:
 * @poisson_type: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_poisson_new (NcDataPoissonType poisson_type)
{
  NcData *data = nc_data_new ();
  static gchar *name = "Poisson";

  data->name                  = name;
  data->type                  = poisson_type;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_init            = NULL;
  data->model_ref             = (NcDataRef) &ncm_mset_func_ref;
  data->model_free            = (NcDataFree) &ncm_mset_func_free;
  data->dts                   = _nc_data_struct_poisson_new ();
  data->resample              = &_poisson_resample;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = &_poisson_calc_leastsquares_f;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_poisson_calc_m2lnL;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

/********************************************************************************************************/

/**
 * nc_data_poisson_set_function:
 * @data: FIXME
 * @func: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_poisson_set_function (NcData *data, NcmMSetFunc *func)
{
  data->model = ncm_mset_func_ref (func);
}

/********************************************************************************************************/

/**
 * nc_data_poisson_set_prepare:
 * @data: FIXME
 * @prepare: (scope call): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_poisson_set_prepare (NcData *data, NcDataPrepare prepare)
{
  data->prepare = prepare;
}

/**
 * nc_data_poisson_set_resample:
 * @data: FIXME
 * @resample: (scope call): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_poisson_set_resample (NcData *data, NcDataResample resample)
{
  data->resample = resample;
}

/********************************************************************************************************/

/**
 * nc_data_poisson_init_from_vector: (skip)
 * @data: FIXME
 * @nodes: FIXME
 * @N: FIXME
 * @extra_data: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_poisson_init_from_vector (NcData *data, NcmVector *nodes, gsl_vector_ulong *N, NcDataStruct *extra_data)
{
  NcDataPoisson *poisson;
	gint i;

	poisson = _poisson_data_alloc (N->size);
	poisson->h->range[0] = ncm_vector_get (nodes, 0);
	for (i = 0; i < N->size; i++)
	{
		poisson->h->range[i + 1] = ncm_vector_get (nodes, i + 1);
		poisson->h->bin[i] = gsl_vector_ulong_get (N, i);
	}

  NC_DATA_STRUCT_DATA (data->dts) = poisson;
  poisson->extra_data = extra_data;

  nc_data_init (data);
}

/**
 * nc_data_poisson_init_from_histogram: (skip)
 * @data: FIXME
 * @h: FIXME
 * @steal: FIXME
 * @extra_data: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_poisson_init_from_histogram (NcData *data, gsl_histogram *h, gboolean steal, NcDataStruct *extra_data)
{
  NcDataPoisson *poisson;

  poisson = g_slice_new0 (NcDataPoisson);
	poisson->h = steal ? h : gsl_histogram_clone (h);

  NC_DATA_STRUCT_DATA (data->dts) = poisson;
  poisson->extra_data = extra_data;

  nc_data_init (data);
}

/**
 * nc_data_poisson_init_zero: (skip)
 * @data: FIXME
 * @nodes: FIXME
 * @extra_data: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_poisson_init_zero (NcData *data, NcmVector *nodes, NcDataStruct *extra_data)
{
  NcDataPoisson *poisson;
	gint i;

	poisson = _poisson_data_alloc (ncm_vector_len (nodes) - 1);
	poisson->h->range[0] = ncm_vector_get (nodes, 0);
	for (i = 0; i < poisson->h->n; i++)
		poisson->h->range[i + 1] = ncm_vector_get (nodes, i + 1);
	gsl_histogram_reset (poisson->h);

  NC_DATA_STRUCT_DATA (data->dts) = poisson;
  poisson->extra_data = extra_data;

  nc_data_init (data);
}

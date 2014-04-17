/***************************************************************************
 *            ncm_data_poisson.c
 *
 *  Sun Apr  4 21:57:39 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_data_poisson
 * @title: Poisson Data
 * @short_description: Poisson data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm_data_poisson.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

enum
{
  PROP_0,
  PROP_NPOINTS,
  PROP_MEAN,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataPoisson, ncm_data_poisson, NCM_TYPE_DATA);

static void
ncm_data_poisson_init (NcmDataPoisson *poisson)
{
  poisson->np       = 0;
  poisson->h        = NULL;
  poisson->log_Nfac = NULL;
}

static void
_ncm_data_poisson_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_poisson_parent_class)->constructed (object);
  {
  }
}

static void
ncm_data_poisson_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (object);
  g_return_if_fail (NCM_IS_DATA_POISSON (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      ncm_data_poisson_set_size (poisson, g_value_get_uint (value));
      break;
    case PROP_MEAN:
    {
      NcmVector *v = ncm_vector_new_data_static (poisson->h->bin, poisson->h->n, 1); 
      ncm_vector_set_from_variant (v, g_value_get_variant (value));
      ncm_vector_free (v);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_poisson_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (object);
  g_return_if_fail (NCM_IS_DATA_POISSON (object));

  switch (prop_id)
  {
    case PROP_NPOINTS:
      g_value_set_uint (value, poisson->np);
      break;
    case PROP_MEAN:
    {
      NcmVector *v = ncm_vector_new_data_static (poisson->h->bin, poisson->h->n, 1); 
      g_value_take_variant (value, ncm_vector_get_variant (v));
      ncm_vector_free (v);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_poisson_dispose (GObject *object)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (object);

  ncm_data_poisson_set_size (poisson, 0);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_poisson_parent_class)->dispose (object);
}

static void
ncm_data_poisson_finalize (GObject *object)
{
  /* NcmDataPoisson *poisson = NCM_DATA_POISSON (object); */
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_poisson_parent_class)->finalize (object);
}

static guint _ncm_data_poisson_get_length (NcmData *data);
static void _ncm_data_poisson_begin (NcmData *data);
static void _ncm_data_poisson_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_poisson_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_poisson_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint np);
static guint _ncm_data_poisson_get_size (NcmDataPoisson *poisson);

static void
ncm_data_poisson_class_init (NcmDataPoissonClass *klass)
{
  GObjectClass* object_class         = G_OBJECT_CLASS (klass);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_CLASS (klass);
  NcmDataClass *data_class           = NCM_DATA_CLASS (klass);

  object_class->constructed  = &_ncm_data_poisson_constructed;
  object_class->set_property = &ncm_data_poisson_set_property;
  object_class->get_property = &ncm_data_poisson_get_property;

  object_class->dispose      = &ncm_data_poisson_dispose;
  object_class->finalize     = &ncm_data_poisson_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NPOINTS,
                                   g_param_spec_uint ("n-points",
                                                      NULL,
                                                      "Data sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MEAN,
                                   g_param_spec_variant ("mean",
                                                         NULL,
                                                         "Data mean",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->get_length     = &_ncm_data_poisson_get_length;
  data_class->begin          = &_ncm_data_poisson_begin;

  data_class->resample       = &_ncm_data_poisson_resample;
  data_class->m2lnL_val      = &_ncm_data_poisson_m2lnL_val;
  data_class->leastsquares_f = &_ncm_data_poisson_leastsquares_f;

  poisson_class->mean_func   = NULL;
  poisson_class->set_size  = &_ncm_data_poisson_set_size;
  poisson_class->get_size  = &_ncm_data_poisson_get_size;
}


static guint 
_ncm_data_poisson_get_length (NcmData *data) 
{ 
  return NCM_DATA_POISSON (data)->np; 
}

static void
_ncm_data_poisson_begin (NcmData *data)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);
  guint i;

  for (i = 0; i < poisson->h->n; i++)
  {
    const gulong N_i = gsl_histogram_get (poisson->h, i);
    ncm_vector_fast_set (poisson->log_Nfac, i, lgamma (N_i + 1.0));
  }
}

static void
_ncm_data_poisson_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (data);
  guint i;

  ncm_rng_lock (rng);
  for (i = 0; i < poisson->h->n; i++)
  {
    const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);
    const gulong N_i = gsl_ran_poisson (rng->r, lambda_i);
    
		poisson->h->bin[i] = N_i;
  }
  ncm_rng_unlock (rng);
}

static void
_ncm_data_poisson_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (data);
  guint i;

  *m2lnL = 0.0;
  if (!ncm_data_bootstrap_enabled (data))
  {
    for (i = 0; i < poisson->h->n; i++)
    {
      const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);
      const gulong N_i = gsl_histogram_get (poisson->h, i);
      const gdouble log_Nfac_i = ncm_vector_fast_get (poisson->log_Nfac, i);
      *m2lnL += -2.0 * ( N_i * log (lambda_i) - lambda_i - log_Nfac_i);
    }
  }
  else
  {
    const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
    for (i = 0; i < bsize; i++)
    {
      guint k = ncm_bootstrap_get (data->bstrap, i);
      const gdouble lambda_k = poisson_class->mean_func (poisson, mset, k);
      const gulong N_k = gsl_histogram_get (poisson->h, k);
      const gdouble log_Nfac_k = ncm_vector_fast_get (poisson->log_Nfac, k);
      *m2lnL += -2.0 * ( N_k * log (lambda_k) - lambda_k - log_Nfac_k);
    }
  }
  return;
}

static void
_ncm_data_poisson_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (data);
  guint i;
  
  if (ncm_data_bootstrap_enabled (data))
    g_error ("NcmDataPoisson: does not support bootstrap with least squares");

  for (i = 0; i < poisson->h->n; i++)
  {
    const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);
    const gulong N_i = gsl_histogram_get (poisson->h, i);

		ncm_vector_set (v, i, (lambda_i - N_i) / sqrt (lambda_i));
  }
}

static void 
_ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint np)
{
  NcmData *data = NCM_DATA (poisson);
  if ((np == 0) || (np != poisson->np))
  {
    poisson->np = 0;
    ncm_vector_clear (&poisson->log_Nfac);
    if (poisson->h != NULL)
    {
      gsl_histogram_free (poisson->h);
      poisson->h = NULL;
    }
    data->init = FALSE;
  }
  if ((np != 0) && (np != poisson->np))
  {
    poisson->np       = np;
    poisson->log_Nfac = ncm_vector_new (poisson->np);
    poisson->h        = gsl_histogram_alloc (poisson->np);
    if (ncm_data_bootstrap_enabled (data))
    {
      ncm_bootstrap_set_fsize (data->bstrap, np);
      ncm_bootstrap_set_bsize (data->bstrap, np);
    }
    data->init = FALSE;
  }
}

static guint 
_ncm_data_poisson_get_size (NcmDataPoisson *poisson)
{
  return poisson->np;
}

/**
 * ncm_data_poisson_init_from_vector: (skip)
 * @data: FIXME
 * @nodes: FIXME
 * @N: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_data_poisson_init_from_vector (NcmData *data, NcmVector *nodes, gsl_vector_ulong *N)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);
	guint i;

  ncm_data_poisson_set_size (poisson, ncm_vector_len (nodes) - 1);
  
	poisson->h->range[0] = ncm_vector_get (nodes, 0);
	for (i = 0; i < N->size; i++)
	{
		poisson->h->range[i + 1] = ncm_vector_get (nodes, i + 1);
		poisson->h->bin[i] = gsl_vector_ulong_get (N, i);
	}

  ncm_data_set_init (data, TRUE);
}

/**
 * ncm_data_poisson_init_from_histogram: (skip)
 * @data: FIXME
 * @h: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_data_poisson_init_from_histogram (NcmData *data, gsl_histogram *h)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);

  ncm_data_poisson_set_size (poisson, h->n);
  gsl_histogram_memcpy (poisson->h, h);

  ncm_data_set_init (data, TRUE);
}

/**
 * ncm_data_poisson_init_zero: (skip)
 * @data: FIXME
 * @nodes: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
ncm_data_poisson_init_zero (NcmData *data, NcmVector *nodes)
{
  NcmDataPoisson *poisson = NCM_DATA_POISSON (data);
	guint i;

  ncm_data_poisson_set_size (poisson, ncm_vector_len (nodes) - 1);
    
	poisson->h->range[0] = ncm_vector_get (nodes, 0);
	for (i = 0; i < poisson->h->n; i++)
		poisson->h->range[i + 1] = ncm_vector_get (nodes, i + 1);
	gsl_histogram_reset (poisson->h);

  ncm_data_set_init (data, TRUE);
}

/**
 * ncm_data_poisson_set_size: (virtual set_size)
 * @poisson: a #NcmDataPoisson
 * @np: data size.
 *
 * Sets the data size to @np.
 * 
 */
void 
ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint np)
{
  NCM_DATA_POISSON_GET_CLASS (poisson)->set_size (poisson, np);
}

/**
 * ncm_data_poisson_get_size: (virtual get_size)
 * @poisson: a #NcmDataPoisson
 *
 * Gets the data size.
 * 
 * Returns: Data size.
 * 
 */
guint 
ncm_data_poisson_get_size (NcmDataPoisson *poisson)
{
  return NCM_DATA_POISSON_GET_CLASS (poisson)->get_size (poisson);
}

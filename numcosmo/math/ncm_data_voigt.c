/***************************************************************************
 *            ncm_data_voigt.c
 *
 *  Thu Mar  15 16:02:50 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:ncm_data_voigt
 * @title: NcmDataVoigt
 * @short_description: Abstract class for implementing voigt distributed data.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm_data_voigt.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_NBINS,
  PROP_MEANS,
  PROP_KNOTS,	 
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmDataVoigt, ncm_data_voigt, NCM_TYPE_DATA);

static void
ncm_data_voigt_init (NcmDataVoigt *voigt)
{
  voigt->nbins       = 0;
  voigt->h        = NULL;
	voigt->means    = NULL;
  voigt->log_Nfac = NULL;
}

static void
_ncm_data_voigt_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_voigt_parent_class)->constructed (object);
  {
  }
}

static void
ncm_data_voigt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (object);
  g_return_if_fail (NCM_IS_DATA_VOIGT (object));

  switch (prop_id)
  {
    case PROP_NBINS:
      ncm_data_voigt_set_size (voigt, g_value_get_uint (value));
      break;
    case PROP_MEANS:
      ncm_vector_memcpy (voigt->means, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_voigt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (object);
  g_return_if_fail (NCM_IS_DATA_VOIGT (object));

  switch (prop_id)
  {
    case PROP_NBINS:
      g_value_set_uint (value, voigt->nbins);
      break;
    case PROP_MEANS:
      g_value_set_object (value, voigt->means);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_voigt_dispose (GObject *object)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (object);

  ncm_data_voigt_set_size (voigt, 0);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_voigt_parent_class)->dispose (object);
}

static void
ncm_data_voigt_finalize (GObject *object)
{
  /* NcmDataVoigt *voigt = NCM_DATA_VOIGT (object); */
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_voigt_parent_class)->finalize (object);
}

static guint _ncm_data_voigt_get_length (NcmData *data);
static void _ncm_data_voigt_begin (NcmData *data);
static void _ncm_data_voigt_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_voigt_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_voigt_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_voigt_set_size (NcmDataVoigt *voigt, guint nbins);
static guint _ncm_data_voigt_get_size (NcmDataVoigt *voigt);

static void
ncm_data_voigt_class_init (NcmDataVoigtClass *klass)
{
  GObjectClass* object_class         = G_OBJECT_CLASS (klass);
  NcmDataVoigtClass *voigt_class = NCM_DATA_VOIGT_CLASS (klass);
  NcmDataClass *data_class           = NCM_DATA_CLASS (klass);

  object_class->constructed  = &_ncm_data_voigt_constructed;
  object_class->set_property = &ncm_data_voigt_set_property;
  object_class->get_property = &ncm_data_voigt_get_property;

  object_class->dispose      = &ncm_data_voigt_dispose;
  object_class->finalize     = &ncm_data_voigt_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NBINS,
                                   g_param_spec_uint ("n-bins",
                                                      NULL,
                                                      "Number of bins",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MEANS,
                                   g_param_spec_object ("mean",
                                                        NULL,
                                                        "Data mean",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->get_length     = &_ncm_data_voigt_get_length;
  data_class->begin          = &_ncm_data_voigt_begin;

  data_class->resample       = &_ncm_data_voigt_resample;
  data_class->m2lnL_val      = &_ncm_data_voigt_m2lnL_val;
  data_class->leastsquares_f = &_ncm_data_voigt_leastsquares_f;

  voigt_class->mean_func   = NULL;
  voigt_class->set_size  = &_ncm_data_voigt_set_size;
  voigt_class->get_size  = &_ncm_data_voigt_get_size;
}


static guint 
_ncm_data_voigt_get_length (NcmData *data) 
{ 
  return NCM_DATA_VOIGT (data)->nbins; 
}

static void
_ncm_data_voigt_begin (NcmData *data)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (data);
  guint i;

  for (i = 0; i < voigt->h->n; i++)
  {
    const gdouble N_i = gsl_histogram_get (voigt->h, i);
    ncm_vector_fast_set (voigt->log_Nfac, i, lgamma (N_i + 1.0));
  }
}

static void
_ncm_data_voigt_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (data);
  NcmDataVoigtClass *voigt_class = NCM_DATA_VOIGT_GET_CLASS (data);
  guint i;

  ncm_rng_lock (rng);
  for (i = 0; i < voigt->h->n; i++)
  {
    const gdouble lambda_i = voigt_class->mean_func (voigt, mset, i);
    const gdouble N_i      = gsl_ran_voigt (rng->r, lambda_i);
    
		voigt->h->bin[i] = N_i;
  }

	ncm_rng_unlock (rng);
}

static void
_ncm_data_voigt_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (data);
  NcmDataVoigtClass *voigt_class = NCM_DATA_VOIGT_GET_CLASS (data);
  guint i;

  *m2lnL = 0.0;
  if (!ncm_data_bootstrap_enabled (data))
  {
    for (i = 0; i < voigt->h->n; i++)
    {
      const gdouble lambda_i   = voigt_class->mean_func (voigt, mset, i);
      const gdouble N_i        = gsl_histogram_get (voigt->h, i);

			if (N_i > 0.0)
        *m2lnL += -2.0 * ( N_i * log (lambda_i / N_i) - lambda_i + N_i);
			else
				*m2lnL += -2.0 * ( - lambda_i);
    }
  }
  else
  {
    const guint bsize = ncm_bootstrap_get_bsize (data->bstrap);
    for (i = 0; i < bsize; i++)
    {
      guint k                = ncm_bootstrap_get (data->bstrap, i);
      const gdouble lambda_k = voigt_class->mean_func (voigt, mset, k);
      const gdouble N_k      = gsl_histogram_get (voigt->h, k);
			
			if (N_k > 0.0)
				*m2lnL += -2.0 * ( N_k * log (lambda_k / N_k) - lambda_k + N_k);
			else
				*m2lnL += -2.0 * ( - lambda_k);
    }
  }
  return;
}

static void
_ncm_data_voigt_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataVoigt *voigt = NCM_DATA_VOIGT (data);
  NcmDataVoigtClass *voigt_class = NCM_DATA_VOIGT_GET_CLASS (data);
  guint i;
  
  if (ncm_data_bootstrap_enabled (data))
    g_error ("NcmDataVoigt: does not support bootstrap with least squares");

  for (i = 0; i < voigt->h->n; i++)
  {
    const gdouble lambda_i = voigt_class->mean_func (voigt, mset, i);
    const gdouble N_i      = gsl_histogram_get (voigt->h, i);
		const gdouble m2lnL_i  = (N_i == 0.0) ? 2.0 * lambda_i : (-2.0 * ( N_i * log (lambda_i / N_i) - lambda_i + N_i));
		
		ncm_vector_set (v, i, sqrt (m2lnL_i));
  }
}

static void 
_ncm_data_voigt_set_size (NcmDataVoigt *voigt, guint nbins)
{
  NcmData *data = NCM_DATA (voigt);
	
  if (nbins != voigt->nbins)
  {
    voigt->nbins = 0;
    ncm_vector_clear (&voigt->log_Nfac);
    if (voigt->h != NULL)
    {
      gsl_histogram_free (voigt->h);
      voigt->h = NULL;
    }
		ncm_vector_clear (&voigt->means);
    data->init = FALSE;

		if (nbins > 0)
		{
			voigt->nbins    = nbins;
			voigt->log_Nfac = ncm_vector_new (voigt->nbins);
			voigt->h        = gsl_histogram_alloc (voigt->nbins);
			voigt->means    = ncm_vector_new_data_static (voigt->h->bin, voigt->h->n, 1); 

			if (ncm_data_bootstrap_enabled (data))
			{
				ncm_bootstrap_set_fsize (data->bstrap, nbins);
				ncm_bootstrap_set_bsize (data->bstrap, nbins);
			}
			data->init = FALSE;
		}		
  }
}

static guint 
_ncm_data_voigt_get_size (NcmDataVoigt *voigt)
{
  return voigt->nbins;
}

/**
 * ncm_data_voigt_init_from_vector:
 * @voigt: a #NcmDataVoigt
 * @nodes: FIXME
 * @N: FIXME
 *
 * FIXME
 *
 */
void
ncm_data_voigt_init_from_vector (NcmDataVoigt *voigt, NcmVector *nodes, NcmVector *N)
{
  guint i;

  ncm_data_voigt_set_size (voigt, ncm_vector_len (nodes) - 1);

	g_assert_cmpuint (ncm_vector_len (nodes), ==, ncm_vector_len (N) + 1);
  
	voigt->h->range[0] = ncm_vector_get (nodes, 0);
	
	for (i = 0; i < ncm_vector_len (N); i++)
	{
		voigt->h->range[i + 1] = ncm_vector_get (nodes, i + 1);
		voigt->h->bin[i]       = ncm_vector_get (N, i);
	}

  ncm_data_set_init (NCM_DATA (voigt), TRUE);
}

/**
 * ncm_data_voigt_init_from_binning:
 * @voigt: a #NcmDataVoigt
 * @nodes: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_data_voigt_init_from_binning (NcmDataVoigt *voigt, NcmVector *nodes, NcmVector *x)
{
  guint i;

  ncm_data_voigt_set_size (voigt, ncm_vector_len (nodes) - 1);

	if (ncm_vector_stride (nodes) == 1)
	{
		gsl_histogram_set_ranges (voigt->h, ncm_vector_data (nodes), ncm_vector_len (nodes));
	}
	else
	{
		for (i = 0; i < ncm_vector_len (nodes); i++)
			voigt->h->range[i] = ncm_vector_get (nodes, i);
	}

	for (i = 0; i < ncm_vector_len (x); i++)
		gsl_histogram_increment (voigt->h, ncm_vector_get (x, i));
	
  ncm_data_set_init (NCM_DATA (voigt), TRUE);
}


/**
 * ncm_data_voigt_init_from_histogram: (skip)
 * @voigt: a #NcmDataVoigt
 * @h: FIXME
 *
 * FIXME
 *
 */
void
ncm_data_voigt_init_from_histogram (NcmDataVoigt *voigt, gsl_histogram *h)
{
  ncm_data_voigt_set_size (voigt, h->n);
  gsl_histogram_memcpy (voigt->h, h);

  ncm_data_set_init (NCM_DATA (voigt), TRUE);
}

/**
 * ncm_data_voigt_init_zero: (skip)
 * @voigt: a #NcmDataVoigt
 * @nodes: a #NcmVector
 *
 * FIXME
 *
 */
void
ncm_data_voigt_init_zero (NcmDataVoigt *voigt, NcmVector *nodes)
{
	guint i;

  ncm_data_voigt_set_size (voigt, ncm_vector_len (nodes) - 1);
    
	voigt->h->range[0] = ncm_vector_get (nodes, 0);
	for (i = 0; i < voigt->h->n; i++)
		voigt->h->range[i + 1] = ncm_vector_get (nodes, i + 1);
	gsl_histogram_reset (voigt->h);

  ncm_data_set_init (NCM_DATA (voigt), TRUE);
}

/**
 * ncm_data_voigt_set_size: (virtual set_size)
 * @voigt: a #NcmDataVoigt
 * @nbins: number of bins.
 *
 * Sets the number of bins to @nbins.
 * 
 */
void 
ncm_data_voigt_set_size (NcmDataVoigt *voigt, guint nbins)
{
  NCM_DATA_VOIGT_GET_CLASS (voigt)->set_size (voigt, nbins);
}

/**
 * ncm_data_voigt_get_size: (virtual get_size)
 * @voigt: a #NcmDataVoigt
 *
 * Gets the data size.
 * 
 * Returns: Data size.
 * 
 */
guint 
ncm_data_voigt_get_size (NcmDataVoigt *voigt)
{
  return NCM_DATA_VOIGT_GET_CLASS (voigt)->get_size (voigt);
}

/**
 * ncm_data_voigt_get_sum:
 * @voigt: a #NcmDataVoigt
 *
 * Gets the sum of all bins.
 * 
 * Returns: Sum of all bins.
 * 
 */
gdouble 
ncm_data_voigt_get_sum (NcmDataVoigt *voigt)
{
	guint i;
	gdouble hsum = 0.0;

	for (i = 0; i < voigt->h->n; i++)
	{
		const gdouble N_i = gsl_histogram_get (voigt->h, i);
		hsum += N_i;
	}

	return hsum;
}

/**
 * ncm_data_voigt_get_hist_vals:
 * @voigt: a #NcmDataVoigt
 *
 * Gets the vector containing the bins values.
 * 
 * Returns: (transfer full): vector containing the bins values.
 */
NcmVector *
ncm_data_voigt_get_hist_vals (NcmDataVoigt *voigt)
{
	guint i;
	NcmVector *v = ncm_vector_new (voigt->h->n);

	for (i = 0; i < voigt->h->n; i++)
	{
		const gdouble N_i = gsl_histogram_get (voigt->h, i);
		ncm_vector_set (v, i, N_i);
	}

	return v;
}

/**
 * ncm_data_voigt_get_hist_means:
 * @voigt: a #NcmDataVoigt
 * @mset: a #NcmMSet
 *
 * Gets the vector containing the bins values.
 * 
 * Returns: (transfer full): vector containing the bins values.
 */
NcmVector *
ncm_data_voigt_get_hist_means (NcmDataVoigt *voigt, NcmMSet *mset)
{
  NcmDataVoigtClass *voigt_class = NCM_DATA_VOIGT_GET_CLASS (voigt);
	NcmVector *v = ncm_vector_new (voigt->h->n);
  guint i;

	ncm_data_prepare (NCM_DATA (voigt), mset);

	for (i = 0; i < voigt->h->n; i++)
	{
    const gdouble lambda_i = voigt_class->mean_func (voigt, mset, i);
		ncm_vector_set (v, i, lambda_i);
	}

	return v;
}

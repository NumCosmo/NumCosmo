/***************************************************************************
 *            ncm_data_poisson.c
 *
 *  Sun Apr  4 21:57:39 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmDataPoisson:
 *
 * Abstract class for binned Poisson data.
 *
 * The data are the counts $N_i$ in #NcmDataPoisson:n-bins bins, with edges
 * #NcmDataPoisson:bin-edges, each drawn from a Poisson distribution of mean
 * $\lambda_i$. Subclasses implement the mean_func virtual method, which returns
 * $\lambda_i$ for bin $i$ given the models in a #NcmMSet; ncm_data_poisson_get_bin_range()
 * gives the edges of the bin.
 *
 * $-2\ln L$ is the Poisson deviance, relative to the saturated model $\lambda_i = N_i$:
 * $$-2\ln L = -2\sum_i \left[N_i \ln(\lambda_i/N_i) - \lambda_i + N_i\right],$$
 * with the term $2\lambda_i$ for an empty bin. The least-squares vector holds the square
 * roots of the terms, and the Fisher matrix uses the variance $\lambda_i$.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/data/ncm_data_poisson.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_histogram.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_NBINS,
  PROP_MEANS,
  PROP_BIN_EDGES,
  PROP_SIZE,
};

typedef struct _NcmDataPoissonPrivate
{
  gsl_histogram *h;
  NcmVector *means;
  guint nbins;
} NcmDataPoissonPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmDataPoisson, ncm_data_poisson, NCM_TYPE_DATA)

static void
ncm_data_poisson_init (NcmDataPoisson *poisson)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  self->nbins = 0;
  self->h     = NULL;
  self->means = NULL;
}

static void _ncm_data_poisson_set_edges (NcmDataPoisson *poisson, NcmVector *edges);

static void
ncm_data_poisson_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (object);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  g_return_if_fail (NCM_IS_DATA_POISSON (object));

  switch (prop_id)
  {
    case PROP_NBINS:
      ncm_data_poisson_set_size (poisson, g_value_get_uint (value));
      break;
    case PROP_MEANS:
    {
      NcmVector *counts = g_value_get_object (value);

      if (counts == NULL)
        break;

      if (ncm_vector_len (counts) != self->nbins)
        g_error ("ncm_data_poisson_set_property: data `%s' has %u bins, but the counts have %u components.",
                 ncm_data_peek_desc (NCM_DATA (poisson)), self->nbins, ncm_vector_len (counts));

      ncm_vector_memcpy (self->means, counts);
      break;
    }
    case PROP_BIN_EDGES:
    {
      NcmVector *edges = g_value_get_object (value);

      if (edges != NULL)
        _ncm_data_poisson_set_edges (poisson, edges);

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_data_poisson_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (object);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  g_return_if_fail (NCM_IS_DATA_POISSON (object));

  switch (prop_id)
  {
    case PROP_NBINS:
      g_value_set_uint (value, self->nbins);
      break;
    case PROP_MEANS:
      g_value_set_object (value, self->means);
      break;
    case PROP_BIN_EDGES:
      g_value_take_object (value, (self->h != NULL) ? ncm_data_poisson_get_bin_edges (poisson) : NULL);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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

static guint _ncm_data_poisson_get_length (NcmData *data);
static void _ncm_data_poisson_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _ncm_data_poisson_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _ncm_data_poisson_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v);
static void _ncm_data_poisson_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu);
static void _ncm_data_poisson_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H);
static void _ncm_data_poisson_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f);
static void _ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint nbins);
static guint _ncm_data_poisson_get_size (NcmDataPoisson *poisson);

static void
ncm_data_poisson_class_init (NcmDataPoissonClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_CLASS (klass);
  NcmDataClass *data_class           = NCM_DATA_CLASS (klass);

  object_class->set_property = &ncm_data_poisson_set_property;
  object_class->get_property = &ncm_data_poisson_get_property;
  object_class->dispose      = &ncm_data_poisson_dispose;

  /**
   * NcmDataPoisson:n-bins:
   *
   * The number of bins; changing it reallocates the data, sets the edges to
   * $0, 1, \dots, n$ and marks the data not initialized.
   */
  g_object_class_install_property (object_class,
                                   PROP_NBINS,
                                   g_param_spec_uint ("n-bins",
                                                      NULL,
                                                      "Number of bins",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmDataPoisson:mean:
   *
   * The counts $N_i$ in each bin.
   */
  g_object_class_install_property (object_class,
                                   PROP_MEANS,
                                   g_param_spec_object ("mean",
                                                        NULL,
                                                        "Data mean",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmDataPoisson:bin-edges:
   *
   * The $n + 1$ bin edges, in increasing order.
   */
  g_object_class_install_property (object_class,
                                   PROP_BIN_EDGES,
                                   g_param_spec_object ("bin-edges",
                                                        NULL,
                                                        "Bin edges",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = TRUE;
  data_class->get_length = &_ncm_data_poisson_get_length;

  data_class->resample       = &_ncm_data_poisson_resample;
  data_class->m2lnL_val      = &_ncm_data_poisson_m2lnL_val;
  data_class->leastsquares_f = &_ncm_data_poisson_leastsquares_f;
  data_class->mean_vector    = &_ncm_data_poisson_mean_vector;
  data_class->inv_cov_UH     = &_ncm_data_poisson_inv_cov_UH;
  data_class->inv_cov_Uf     = &_ncm_data_poisson_inv_cov_Uf;

  poisson_class->mean_func = NULL;
  poisson_class->set_size  = &_ncm_data_poisson_set_size;
  poisson_class->get_size  = &_ncm_data_poisson_get_size;
}

static guint
_ncm_data_poisson_get_length (NcmData *data)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  return self->nbins;
}

static void
_ncm_data_poisson_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (data);
  guint i;

  ncm_rng_lock (rng);

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);
    const gdouble N_i      = ncm_rng_poisson_gen (rng, lambda_i);

    self->h->bin[i] = N_i;
  }

  ncm_rng_unlock (rng);
}

static void
_ncm_data_poisson_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (data);
  guint i;

  *m2lnL = 0.0;

  if (!ncm_data_bootstrap_enabled (data))
  {
    for (i = 0; i < self->h->n; i++)
    {
      const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);
      const gdouble N_i      = gsl_histogram_get (self->h, i);

      if (N_i > 0.0)
        *m2lnL += -2.0 * (N_i * log (lambda_i / N_i) - lambda_i + N_i);
      else
        *m2lnL += -2.0 * (-lambda_i);
    }
  }
  else
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
    const guint bsize    = ncm_bootstrap_get_bsize (bstrap);

    for (i = 0; i < bsize; i++)
    {
      guint k                = ncm_bootstrap_get (bstrap, i);
      const gdouble lambda_k = poisson_class->mean_func (poisson, mset, k);
      const gdouble N_k      = gsl_histogram_get (self->h, k);

      if (N_k > 0.0)
        *m2lnL += -2.0 * (N_k * log (lambda_k / N_k) - lambda_k + N_k);
      else
        *m2lnL += -2.0 * (-lambda_k);
    }
  }

  return;
}

static void
_ncm_data_poisson_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *v)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (data);
  guint i;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("_ncm_data_poisson_leastsquares_f: data `%s': bootstrap is not supported with least squares.",
             ncm_data_peek_desc (data));

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);
    const gdouble N_i      = gsl_histogram_get (self->h, i);
    const gdouble m2lnL_i  = (N_i == 0.0) ? 2.0 * lambda_i : (-2.0 * (N_i * log (lambda_i / N_i) - lambda_i + N_i));

    ncm_vector_set (v, i, sqrt (m2lnL_i));
  }
}

static void
_ncm_data_poisson_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (poisson);
  guint i;

  if (ncm_vector_len (mu) != self->h->n)
    g_error ("_ncm_data_poisson_mean_vector: data `%s' has %u bins, but the vector has %u components.",
             ncm_data_peek_desc (data), (guint) self->h->n, ncm_vector_len (mu));

  for (i = 0; i < self->h->n; i++)
  {
    ncm_vector_set (mu, i, poisson_class->mean_func (poisson, mset, i));
  }
}

static void
_ncm_data_poisson_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (poisson);
  guint i;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("_ncm_data_poisson_inv_cov_UH: data `%s': bootstrap is not supported with the Fisher matrix.",
             ncm_data_peek_desc (data));

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble mean_i = poisson_class->mean_func (poisson, mset, i);

    ncm_matrix_mul_col (H, i, 1.0 / sqrt (mean_i));
  }
}

static void
_ncm_data_poisson_inv_cov_Uf (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  NcmDataPoisson *poisson            = NCM_DATA_POISSON (data);
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (poisson);
  guint i;

  if (ncm_data_bootstrap_enabled (data))
    g_error ("_ncm_data_poisson_inv_cov_Uf: data `%s': bootstrap is not supported with the Fisher matrix.",
             ncm_data_peek_desc (data));

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble mean_i = poisson_class->mean_func (poisson, mset, i);

    ncm_vector_mulby (f, i, 1.0 / sqrt (mean_i));
  }
}

static void
_ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint nbins)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmData *data                      = NCM_DATA (poisson);

  if (nbins != self->nbins)
  {
    self->nbins = 0;

    if (self->h != NULL)
    {
      gsl_histogram_free (self->h);
      self->h = NULL;
    }

    ncm_vector_clear (&self->means);
    ncm_data_set_init (data, FALSE);

    if (nbins > 0)
    {
      NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

      self->nbins = nbins;
      self->h     = gsl_histogram_alloc (self->nbins);
      self->means = ncm_vector_new_data_static (self->h->bin, self->h->n, 1);

      /* gsl_histogram_alloc leaves the ranges undefined. */
      gsl_histogram_set_ranges_uniform (self->h, 0.0, self->nbins);
      gsl_histogram_reset (self->h);

      if (ncm_data_bootstrap_enabled (data))
      {
        ncm_bootstrap_set_fsize (bstrap, nbins);
        ncm_bootstrap_set_bsize (bstrap, nbins);
      }

      ncm_data_set_init (data, FALSE);
    }
  }
}

static guint
_ncm_data_poisson_get_size (NcmDataPoisson *poisson)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  return self->nbins;
}

/*
 * Sets the bin edges from @edges, which must have n-bins + 1 increasing components.
 */
static void
_ncm_data_poisson_set_edges (NcmDataPoisson *poisson, NcmVector *edges)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  const guint len                    = ncm_vector_len (edges);
  guint i;

  if (len != self->nbins + 1)
    g_error ("ncm_data_poisson: data `%s' has %u bins, but %u bin edges were given (expected %u).",
             ncm_data_peek_desc (NCM_DATA (poisson)), self->nbins, len, self->nbins + 1);

  for (i = 0; i + 1 < len; i++)
  {
    if (!(ncm_vector_get (edges, i) < ncm_vector_get (edges, i + 1)))
      g_error ("ncm_data_poisson: data `%s': the bin edges must be increasing, but edge %u is %g and edge %u is %g.",
               ncm_data_peek_desc (NCM_DATA (poisson)), i, ncm_vector_get (edges, i), i + 1, ncm_vector_get (edges, i + 1));
  }

  for (i = 0; i < len; i++)
    self->h->range[i] = ncm_vector_get (edges, i);
}

/* Resizes to ncm_vector_len (@nodes) - 1 bins, which must be at least one. */
static void
_ncm_data_poisson_resize_from_nodes (NcmDataPoisson *poisson, NcmVector *nodes, const gchar *func)
{
  if (ncm_vector_len (nodes) < 2)
    g_error ("%s: data `%s': at least two bin edges are needed, but %u were given.",
             func, ncm_data_peek_desc (NCM_DATA (poisson)), ncm_vector_len (nodes));

  ncm_data_poisson_set_size (poisson, ncm_vector_len (nodes) - 1);
}

/**
 * ncm_data_poisson_init_from_vector:
 * @poisson: a #NcmDataPoisson
 * @nodes: bin edges
 * @N: counts in each bin
 *
 * Sets the bins to those with edges @nodes and the counts to @N, which must have one
 * component less than @nodes, and marks the data initialized.
 *
 */
void
ncm_data_poisson_init_from_vector (NcmDataPoisson *poisson, NcmVector *nodes, NcmVector *N)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  guint i;

  if (ncm_vector_len (nodes) != ncm_vector_len (N) + 1)
    g_error ("ncm_data_poisson_init_from_vector: data `%s': %u bin edges for %u counts, expected %u.",
             ncm_data_peek_desc (NCM_DATA (poisson)), ncm_vector_len (nodes), ncm_vector_len (N), ncm_vector_len (N) + 1);

  _ncm_data_poisson_resize_from_nodes (poisson, nodes, "ncm_data_poisson_init_from_vector");
  _ncm_data_poisson_set_edges (poisson, nodes);

  for (i = 0; i < ncm_vector_len (N); i++)
    self->h->bin[i] = ncm_vector_get (N, i);

  ncm_data_set_init (NCM_DATA (poisson), TRUE);
}

/**
 * ncm_data_poisson_init_from_binning:
 * @poisson: a #NcmDataPoisson
 * @nodes: bin edges
 * @x: data to be binned
 *
 * Sets the bins to those with edges @nodes, counts the values in @x that fall in each
 * bin (values outside the edges are ignored) and marks the data initialized.
 *
 */
void
ncm_data_poisson_init_from_binning (NcmDataPoisson *poisson, NcmVector *nodes, NcmVector *x)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  guint i;

  _ncm_data_poisson_resize_from_nodes (poisson, nodes, "ncm_data_poisson_init_from_binning");
  _ncm_data_poisson_set_edges (poisson, nodes);
  gsl_histogram_reset (self->h);

  for (i = 0; i < ncm_vector_len (x); i++)
    gsl_histogram_increment (self->h, ncm_vector_get (x, i));

  ncm_data_set_init (NCM_DATA (poisson), TRUE);
}

/**
 * ncm_data_poisson_init_zero:
 * @poisson: a #NcmDataPoisson
 * @nodes: bin edges
 *
 * Sets the bins to those with edges @nodes, with zero counts, and marks the data
 * initialized.
 *
 */
void
ncm_data_poisson_init_zero (NcmDataPoisson *poisson, NcmVector *nodes)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  _ncm_data_poisson_resize_from_nodes (poisson, nodes, "ncm_data_poisson_init_zero");
  _ncm_data_poisson_set_edges (poisson, nodes);
  gsl_histogram_reset (self->h);

  ncm_data_set_init (NCM_DATA (poisson), TRUE);
}

/**
 * ncm_data_poisson_set_size: (virtual set_size)
 * @poisson: a #NcmDataPoisson
 * @nbins: number of bins.
 *
 * Sets the number of bins to @nbins.
 *
 */
void
ncm_data_poisson_set_size (NcmDataPoisson *poisson, guint nbins)
{
  NCM_DATA_POISSON_GET_CLASS (poisson)->set_size (poisson, nbins);
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

/**
 * ncm_data_poisson_get_sum:
 * @poisson: a #NcmDataPoisson
 *
 * Gets the sum of all bins.
 *
 * Returns: Sum of all bins.
 *
 */
gdouble
ncm_data_poisson_get_sum (NcmDataPoisson *poisson)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  gdouble hsum                       = 0.0;
  guint i;

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble N_i = gsl_histogram_get (self->h, i);

    hsum += N_i;
  }

  return hsum;
}

/**
 * ncm_data_poisson_get_hist_vals:
 * @poisson: a #NcmDataPoisson
 *
 * Gets the counts $N_i$ in each bin.
 *
 * Returns: (transfer full): a new vector with the counts.
 */
NcmVector *
ncm_data_poisson_get_hist_vals (NcmDataPoisson *poisson)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmVector *v                       = ncm_vector_new (self->h->n);
  guint i;

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble N_i = gsl_histogram_get (self->h, i);

    ncm_vector_set (v, i, N_i);
  }

  return v;
}

/**
 * ncm_data_poisson_get_hist_means:
 * @poisson: a #NcmDataPoisson
 * @mset: a #NcmMSet
 *
 * Computes the mean $\lambda_i$ of each bin given the models in @mset.
 *
 * Returns: (transfer full): the means $\lambda_i$.
 */
NcmVector *
ncm_data_poisson_get_hist_means (NcmDataPoisson *poisson, NcmMSet *mset)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_GET_CLASS (poisson);
  NcmVector *v                       = ncm_vector_new (self->h->n);
  guint i;

  ncm_data_prepare (NCM_DATA (poisson), mset);

  for (i = 0; i < self->h->n; i++)
  {
    const gdouble lambda_i = poisson_class->mean_func (poisson, mset, i);

    ncm_vector_set (v, i, lambda_i);
  }

  return v;
}

/**
 * ncm_data_poisson_get_bin_edges:
 * @poisson: a #NcmDataPoisson
 *
 * Gets the n-bins + 1 bin edges.
 *
 * Returns: (transfer full): a new vector with the bin edges.
 */
NcmVector *
ncm_data_poisson_get_bin_edges (NcmDataPoisson *poisson)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);
  NcmVector *edges;
  guint i;

  if (self->h == NULL)
    g_error ("ncm_data_poisson_get_bin_edges: data `%s' has no bins.", ncm_data_peek_desc (NCM_DATA (poisson)));

  edges = ncm_vector_new (self->nbins + 1);

  for (i = 0; i <= self->nbins; i++)
    ncm_vector_set (edges, i, self->h->range[i]);

  return edges;
}

/**
 * ncm_data_poisson_get_bin_range:
 * @poisson: a #NcmDataPoisson
 * @i: bin index
 * @lower: (out): lower edge of bin @i
 * @upper: (out): upper edge of bin @i
 *
 * Gets the edges of bin @i, for use in mean_func.
 *
 */
void
ncm_data_poisson_get_bin_range (NcmDataPoisson *poisson, const guint i, gdouble *lower, gdouble *upper)
{
  NcmDataPoissonPrivate * const self = ncm_data_poisson_get_instance_private (poisson);

  if (i >= self->nbins)
    g_error ("ncm_data_poisson_get_bin_range: data `%s' has %u bins, bin %u requested.",
             ncm_data_peek_desc (NCM_DATA (poisson)), self->nbins, i);

  *lower = self->h->range[i];
  *upper = self->h->range[i + 1];
}


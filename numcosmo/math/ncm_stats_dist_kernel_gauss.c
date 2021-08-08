/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel_gauss.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel_gauss.c
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
 * SECTION:ncm_stats_dist_kernel_gauss
 * @title: NcmStatsDistKernelGauss
 * @short_description: An N dimensional probability distributions using gaussian KDE
 *
 * An arbitrary N dimensional probability distribution using gaussian KDE.
 *
 *
 * This object provides the tools to perform a radial basis interpolation
 * in a multidimensional function, using a multivariate gaussian distribution.
 * such that we are able to sample the original functio from this new interpolation function.
 * For more informations about radial basis interpolation,
 * check #NcmStatsDist.
 * A brief description of the Multivariate gaussian function can be found below.
 * For more information, check [[The R Journal Vol. 5/2, December 2013](https://journal.r-project.org/archive/2013/RJ-2013-033/RJ-2013-033.pdf)]
 *
 * In this file, we use the Multivariate gaussian function as the radial basis function. The function has the stocastic representation given by
 *
 * \begin{align}
 * \boldsymbol{X}=\boldsymbol{\mu}+ A \boldsymbol{Z}
 * ,\end{align}
 * where $\boldsymbol{Z}$ is a p-dimensional random vector, $A$ is a $p \times p$ matrix and $\mu$ is the mean vector.
 *
 *$\boldsymbol{X}$ is fully determined by the covariance matrix $\Sigma = \boldsymbol{A}\boldsymbol{A}^t$ and the mean vector $\mu$.
 * Assuming that the covariance matrix is positive definite, $\boldsymbol{X}$ has the probability density
 *
 * \begin{align}
 * f_{X}(x)=\frac{1}{(2 \pi)^{d / 2} \sqrt{\operatorname{det} \Sigma}} \exp \left(-\frac{1}{2}(x-\mu)^{\top} \Sigma^{-1}(x-\mu)\right), \quad x \in \mathbb{R}^{d}
 *, \end{align}
 * where $p$ is the dimension and $x$ are the points to be evaluated.
 *
 * Once this object is initialized, we may use the methods in the #NcmStatsDistKernelClass to perform the interpolation
 * and to generate a sample from the interpolated polinomal function.
 *
 * The user must provide the following input values: $p$ - ncm_stats_dist_kernel_gauss_new(),
 * cv_type - ncm_stats_dist_kernel_gauss_new().
 **/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_kernel_gauss.h"
#include "math/ncm_stats_vec.h"
#include "math/ncm_c.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "math/ncm_stats_dist_kernel_private.h"

struct _NcmStatsDistKernelGaussPrivate
{
  GArray *t_array;
};

enum
{
  PROP_0,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistKernelGauss, ncm_stats_dist_kernel_gauss, NCM_TYPE_STATS_DIST_KERNEL);

static void
ncm_stats_dist_kernel_gauss_init (NcmStatsDistKernelGauss *sdkg)
{
  NcmStatsDistKernelGaussPrivate * const self = sdkg->priv = ncm_stats_dist_kernel_gauss_get_instance_private (sdkg);
  
  self->t_array = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void
_ncm_stats_dist_kernel_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcmStatsDistKernelGauss *sdkg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistKernelGaussPrivate * const self = sdkg->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL_GAUSS (object));
  
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_kernel_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcmStatsDistKernelGauss *sdkg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistKernelGaussPrivate * const self = sdkg->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL_GAUSS (object));
  
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_kernel_gauss_dispose (GObject *object)
{
  NcmStatsDistKernelGauss *sdkg               = NCM_STATS_DIST_KERNEL_GAUSS (object);
  NcmStatsDistKernelGaussPrivate * const self = sdkg->priv;
  
  g_clear_pointer (&self->t_array, g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kernel_gauss_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_kernel_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kernel_gauss_parent_class)->finalize (object);
}

static gdouble _ncm_stats_dist_kernel_gauss_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n);
static gdouble _ncm_stats_dist_kernel_gauss_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp);
static gdouble _ncm_stats_dist_kernel_gauss_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2);
static void _ncm_stats_dist_kernel_gauss_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku);
static void _ncm_stats_dist_kernel_gauss_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, gdouble *gamma, gdouble *lambda);
static void _ncm_stats_dist_kernel_gauss_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, gdouble *gamma, gdouble *lambda);
static void _ncm_stats_dist_kernel_gauss_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *x, NcmRNG *rng);

static void
ncm_stats_dist_kernel_gauss_class_init (NcmStatsDistKernelGaussClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmStatsDistKernelClass *sdk_class = NCM_STATS_DIST_KERNEL_CLASS (klass);
  
  object_class->set_property = &_ncm_stats_dist_kernel_gauss_set_property;
  object_class->get_property = &_ncm_stats_dist_kernel_gauss_get_property;
  object_class->dispose      = &_ncm_stats_dist_kernel_gauss_dispose;
  object_class->finalize     = &_ncm_stats_dist_kernel_gauss_finalize;
  
  sdk_class->get_rot_bandwidth      = &_ncm_stats_dist_kernel_gauss_get_rot_bandwidth;
  sdk_class->get_lnnorm             = &_ncm_stats_dist_kernel_gauss_get_lnnorm;
  sdk_class->eval_unnorm            = &_ncm_stats_dist_kernel_gauss_eval_unnorm;
  sdk_class->eval_unnorm_vec        = &_ncm_stats_dist_kernel_gauss_eval_unnorm_vec;
  sdk_class->eval_sum0_gamma_lambda = &_ncm_stats_dist_kernel_gauss_eval_sum0_gamma_lambda;
  sdk_class->eval_sum1_gamma_lambda = &_ncm_stats_dist_kernel_gauss_eval_sum1_gamma_lambda;
  sdk_class->sample                 = &_ncm_stats_dist_kernel_gauss_sample;
}

static gdouble
_ncm_stats_dist_kernel_gauss_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n)
{
  NcmStatsDistKernelPrivate * const pself = sdk->priv;

  return pow (4.0 / (n * (pself->d + 2.0)), 1.0 / (pself->d + 4.0));
}

static gdouble
_ncm_stats_dist_kernel_gauss_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  NcmStatsDistKernelPrivate * const pself = sdk->priv;
  
  return 0.5 * (pself->d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (cov_decomp));
}

static gdouble
_ncm_stats_dist_kernel_gauss_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2)
{
  return exp (-0.5 * chi2);
}

static void
_ncm_stats_dist_kernel_gauss_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku)
{
  /*NcmStatsDistKernelPrivate * const pself = sdk->priv;*/
  const guint n = ncm_vector_len (chi2);
  gint i;
  
  g_assert (ncm_vector_len (Ku) == n);
  
  if ((ncm_vector_stride (Ku) == 1) && (ncm_vector_stride (chi2) == 1))
  {
    for (i = 0; i < n; i++)
    {
      const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
      const gdouble Ku_i   = _ncm_stats_dist_kernel_gauss_eval_unnorm (sdk, chi2_i);
      
      ncm_vector_fast_set (Ku, i, Ku_i);
    }
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      const gdouble chi2_i = ncm_vector_get (chi2, i);
      const gdouble Ku_i   = _ncm_stats_dist_kernel_gauss_eval_unnorm (sdk, chi2_i);
      
      ncm_vector_set (Ku, i, Ku_i);
    }
  }
}

static void
_ncm_stats_dist_kernel_gauss_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, gdouble *gamma, gdouble *lambda)
{
  NcmStatsDistKernelGaussPrivate * const self = NCM_STATS_DIST_KERNEL_GAUSS (sdk)->priv;
  /*NcmStatsDistKernelPrivate * const pself  = sdk->priv;*/
  
  const guint n = ncm_vector_len (chi2);
  gdouble lnt_max = GSL_NEGINF;
  gint i, i_max = -1;
  
  g_assert (n == ncm_vector_len (weights));
  g_assert (n == ncm_vector_len (lnnorms));
  g_assert (1 == ncm_vector_stride (chi2));
  g_assert (1 == ncm_vector_stride (weights));
  g_assert (1 == ncm_vector_stride (lnnorms));
  
  g_array_set_size (self->t_array, n);
  
  for (i = 0; i < n; i++)
  {
    const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
    const gdouble w_i    = ncm_vector_fast_get (weights, i);
    const gdouble lnu_i  = ncm_vector_fast_get (lnnorms, i);
    
    const gdouble lnt_i = -0.5 * chi2_i - lnu_i + log (w_i);
    
    if (lnt_i > lnt_max)
    {
      i_max   = i;
      lnt_max = lnt_i;
    }
    
    g_array_index (self->t_array, gdouble, i) = lnt_i;
  }
  
  lambda[0] = 0.0;
  
  for (i = 0; i < i_max; i++)
    lambda[0] += exp (g_array_index (self->t_array, gdouble, i) - lnt_max);
  
  for (i = i_max + 1; i < n; i++)
    lambda[0] += exp (g_array_index (self->t_array, gdouble, i) - lnt_max);
  
  gamma[0] = lnt_max;
}

static void
_ncm_stats_dist_kernel_gauss_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, gdouble *gamma, gdouble *lambda)
{
  NcmStatsDistKernelGaussPrivate * const self = NCM_STATS_DIST_KERNEL_GAUSS (sdk)->priv;
  /*NcmStatsDistKernelPrivate * const pself  = sdk->priv;*/
  
  const guint n = ncm_vector_len (chi2);
  gdouble lnt_max = GSL_NEGINF;
  gint i, i_max = -1;
  
  g_assert (n == ncm_vector_len (weights));
  g_assert (1 == ncm_vector_stride (chi2));
  g_assert (1 == ncm_vector_stride (weights));
  
  g_array_set_size (self->t_array, n);
  
  for (i = 0; i < n; i++)
  {
    const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
    const gdouble w_i    = ncm_vector_fast_get (weights, i);
    
    const gdouble lnt_i = -0.5 * chi2_i + log (w_i);
    
    if (lnt_i > lnt_max)
    {
      i_max   = i;
      lnt_max = lnt_i;
    }
    
    g_array_index (self->t_array, gdouble, i) = lnt_i;
  }
  
  lambda[0] = 0.0;
  
  for (i = 0; i < i_max; i++)
    lambda[0] += exp (g_array_index (self->t_array, gdouble, i) - lnt_max);
  
  for (i = i_max + 1; i < n; i++)
    lambda[0] += exp (g_array_index (self->t_array, gdouble, i) - lnt_max);
  
  gamma[0] = lnt_max - lnnorm;
}

static void
_ncm_stats_dist_kernel_gauss_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistKernelPrivate * const pself = sdk->priv;
  gint i, ret;
  
  for (i = 0; i < pself->d; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    
    ncm_vector_set (x, i, u_i * href);
  }
  
  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (x));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_kernel_gauss_sample", ret);

  ncm_vector_add (x, mu);
}

/**
 * ncm_stats_dist_kernel_gauss_new:
 * @dim: sample space dimension
 *
 * Creates a new #NcmStatsDistKernelGauss object with sample dimension @dim.
 *
 * Returns: a new #NcmStatsDistKernelGauss.
 */
NcmStatsDistKernelGauss *
ncm_stats_dist_kernel_gauss_new (const guint dim)
{
  NcmStatsDistKernelGauss *sdkg = g_object_new (NCM_TYPE_STATS_DIST_KERNEL_GAUSS,
                                                "dimension", dim,
                                                NULL);
  
  return sdkg;
}

/**
 * ncm_stats_dist_kernel_gauss_ref:
 * @sdkg: a #NcmStatsDistKernelGauss
 *
 * Increase the reference of @stats_dist_nd_vbk_gauss by one.
 *
 * Returns: (transfer full): @stats_dist_nd_vbk_gauss.
 */
NcmStatsDistKernelGauss *
ncm_stats_dist_kernel_gauss_ref (NcmStatsDistKernelGauss *sdkg)
{
  return g_object_ref (sdkg);
}

/**
 * ncm_stats_dist_kernel_gauss_free:
 * @sdkg: a #NcmStatsDistKernelGauss
 *
 * Decrease the reference count of @stats_dist_nd_vbk_gauss by one.
 *
 */
void
ncm_stats_dist_kernel_gauss_free (NcmStatsDistKernelGauss *sdkg)
{
  g_object_unref (sdkg);
}

/**
 * ncm_stats_dist_kernel_gauss_clear:
 * @sdkg: a #NcmStatsDistKernelGauss
 *
 * Decrease the reference count of @stats_dist_nd_vbk_gauss by one, and sets the pointer *@stats_dist_nd_vbk_gauss to
 * NULL.
 *
 */
void
ncm_stats_dist_kernel_gauss_clear (NcmStatsDistKernelGauss **sdkg)
{
  g_clear_object (sdkg);
}


/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel.c
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
 * SECTION:ncm_stats_dist_kernel
 * @title: NcmStatsDistKernel
 * @short_description: An N dimensional probability distributions using gaussian KDE
 *
 * An arbitrary N dimensional probability distribution using gaussian KDE.
 *
 *
 * This object provides the tools to perform a radial basis interpolation
 * in a multidimensional function, using a multivariate gaussian distribution.
 * such that we are able to sample the original functio from this new interpolation function.
 * For more informations about radial basis interpolation,
 * check #NcmStatsDistKernel.
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
 * Once this object is initialized, we may use the methods in the #NcmStatsDistKernel class to perform the interpolation
 * and to generate a sample from the interpolated polynomial function.
 *
 * The user must provide the following input values: $p$ - ncm_stats_dist_kernel_new(),
 * cv_type - ncm_stats_dist_kernel_new().
 **/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_kernel.h"
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

enum
{
  PROP_0,
  PROP_DIM,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDistKernel, ncm_stats_dist_kernel, G_TYPE_OBJECT);

static void
ncm_stats_dist_kernel_init (NcmStatsDistKernel *sdk)
{
  NcmStatsDistKernelPrivate * const self = sdk->priv = ncm_stats_dist_kernel_get_instance_private (sdk);
  
  self->d = 0;
}

static void
_ncm_stats_dist_kernel_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernel *sdk = NCM_STATS_DIST_KERNEL (object);
  
  /*NcmStatsDistKernelPrivate * const self = sdk->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL (object));
  
  switch (prop_id)
  {
    case PROP_DIM:
      NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->set_dim (sdk, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_kernel_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernel *sdk = NCM_STATS_DIST_KERNEL (object);
  
  /*NcmStatsDistKernelPrivate * const self = sdk->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL (object));
  
  switch (prop_id)
  {
    case PROP_DIM:
      g_value_set_uint (value, ncm_stats_dist_kernel_get_dim (sdk));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_kernel_dispose (GObject *object)
{
  /*NcmStatsDistKernel *sdk = NCM_STATS_DIST_KERNEL (object);*/
  /*NcmStatsDistKernelPrivate * const self = sdk->priv;*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kernel_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_kernel_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kernel_parent_class)->finalize (object);
}

static void _ncm_stats_dist_kernel_set_dim (NcmStatsDistKernel *sdk, const guint dim);
static guint _ncm_stats_dist_kernel_get_dim (NcmStatsDistKernel *sdk);

static gdouble
_ncm_stats_dist_kernel_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n)
{
  g_error ("method get_rot_bandwidth not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
  
  return 0.0;
}

static gdouble
_ncm_stats_dist_kernel_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  g_error ("method get_lnnorm not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
  
  return 0.0;
}

static gdouble
_ncm_stats_dist_kernel_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2)
{
  g_error ("method eval_unnorm not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
  
  return 0.0;
}

static void
_ncm_stats_dist_kernel_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku)
{
  g_error ("method eval_unnorm_vec not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
}

static void
_ncm_stats_dist_kernel_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, gdouble *gamma, gdouble *lambda)
{
  g_error ("method eval_sum0_gamma_lambda not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
}

static void
_ncm_stats_dist_kernel_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, gdouble *gamma, gdouble *lambda)
{
  g_error ("method eval_sum1_gamma_lambda not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
}

static void
_ncm_stats_dist_kernel_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng)
{
  g_error ("method sample not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
}

static void
ncm_stats_dist_kernel_class_init (NcmStatsDistKernelClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmStatsDistKernelClass *sd_class = NCM_STATS_DIST_KERNEL_CLASS (klass);
  
  object_class->set_property = &_ncm_stats_dist_kernel_set_property;
  object_class->get_property = &_ncm_stats_dist_kernel_get_property;
  object_class->dispose      = &_ncm_stats_dist_kernel_dispose;
  object_class->finalize     = &_ncm_stats_dist_kernel_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dimension",
                                                      NULL,
                                                      "Kernel dimension",
                                                      2, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  sd_class->set_dim                = &_ncm_stats_dist_kernel_set_dim;
  sd_class->get_dim                = &_ncm_stats_dist_kernel_get_dim;
  sd_class->get_rot_bandwidth      = &_ncm_stats_dist_kernel_get_rot_bandwidth;
  sd_class->get_lnnorm             = &_ncm_stats_dist_kernel_get_lnnorm;
  sd_class->eval_unnorm            = &_ncm_stats_dist_kernel_eval_unnorm;
  sd_class->eval_unnorm_vec        = &_ncm_stats_dist_kernel_eval_unnorm_vec;
  sd_class->eval_sum0_gamma_lambda = &_ncm_stats_dist_kernel_eval_sum0_gamma_lambda;
  sd_class->eval_sum1_gamma_lambda = &_ncm_stats_dist_kernel_eval_sum1_gamma_lambda;
  sd_class->sample                 = &_ncm_stats_dist_kernel_sample;
}

static void
_ncm_stats_dist_kernel_set_dim (NcmStatsDistKernel *sdk, const guint dim)
{
  NcmStatsDistKernelPrivate * const self = sdk->priv;
  
  self->d = dim;
}

static guint
_ncm_stats_dist_kernel_get_dim (NcmStatsDistKernel *sdk)
{
  NcmStatsDistKernelPrivate * const self = sdk->priv;
  
  return self->d;
}

/**
 * ncm_stats_dist_kernel_ref:
 * @sdk: a #NcmStatsDistKernel
 *
 * Increase the reference of @sdk by one.
 *
 * Returns: (transfer full): @sdk.
 */
NcmStatsDistKernel *
ncm_stats_dist_kernel_ref (NcmStatsDistKernel *sdk)
{
  return g_object_ref (sdk);
}

/**
 * ncm_stats_dist_kernel_free:
 * @sdk: a #NcmStatsDistKernel
 *
 * Decrease the reference count of @sdk by one.
 *
 */
void
ncm_stats_dist_kernel_free (NcmStatsDistKernel *sdk)
{
  g_object_unref (sdk);
}

/**
 * ncm_stats_dist_kernel_clear:
 * @sdk: a #NcmStatsDistKernel
 *
 * Decrease the reference count of @stats_dist_nd_kde_gauss by one, and sets the pointer *@sdk to
 * NULL.
 *
 */
void
ncm_stats_dist_kernel_clear (NcmStatsDistKernel **sdk)
{
  g_clear_object (sdk);
}

/**
 * ncm_stats_dist_kernel_get_dim: (virtual get_dim)
 * @sdk: a #NcmStatsDistKernel
 *
 * Gets current kernel dimension.
 *
 * Returns: current kernel dimension.
 */
guint
ncm_stats_dist_kernel_get_dim (NcmStatsDistKernel *sdk)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->get_dim (sdk);
}

/**
 * ncm_stats_dist_kernel_get_rot_bandwidth: (virtual get_rot_bandwidth)
 * @sdk: a #NcmStatsDistKernel
 * @n: number of kernels
 *
 * Computes the rule-of-thumb bandwidth for a interpolation
 * using @n kernels.
 *
 * Returns: the rule-of-thumb bandwidth.
 */
gdouble
ncm_stats_dist_kernel_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->get_rot_bandwidth (sdk, n);
}

/**
 * ncm_stats_dist_kernel_get_lnnorm: (virtual get_lnnorm)
 * @sdk: a #NcmStatsDistKernel
 * @cov_decomp: Cholesky decomposition of the kernel covariance
 *
 * Computes the kernel normalization for a given covariance @cov_decomp.
 *
 * Returns: the kernel normalization logarithm.
 */
gdouble
ncm_stats_dist_kernel_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->get_lnnorm (sdk, cov_decomp);
}

/**
 * ncm_stats_dist_kernel_eval_unnorm: (virtual eval_unnorm)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: a double
 *
 * Computes the unnormalized kernel at $\chi^2=$@chi2.
 *
 * Returns: the unnormalized kernel at $\chi^2=$@chi2.
 */
gdouble
ncm_stats_dist_kernel_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_unnorm (sdk, chi2);
}

/**
 * ncm_stats_dist_kernel_eval_unnorm_vec: (virtual eval_unnorm_vec)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: a #NcmVector
 * @Ku: a #NcmVector
 *
 * Computes the unnormalized kernel at $\chi^2=$@chi2 for all elements of @chi2
 * and store the results at @Ku.
 *
 */
void
ncm_stats_dist_kernel_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku)
{
  NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_unnorm_vec (sdk, chi2, Ku);
}

/**
 * ncm_stats_dist_kernel_eval_sum0_gamma_lambda: (virtual eval_sum0_gamma_lambda)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: a #NcmVector
 * @weights: a #NcmVector
 * @lnnorms: a #NcmVector
 * @gamma: (out): $\gamma$
 * @lambda: (out): $\lambda$
 *
 * Computes the weighted sum of kernels at $\chi^2=$@chi2,
 * $$ e^\gamma (1+\lambda) = \sum_i w_i\bar{K} (\chi^2_i) / u_i,$$
 * where $\gamma = \ln(w_a\bar{K} (\chi^2_a) / u_a)$ and $a$ labels
 * is the largest term of the sum.
 *
 */
void
ncm_stats_dist_kernel_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, gdouble *gamma, gdouble *lambda)
{
  NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_sum0_gamma_lambda (sdk, chi2, weights, lnnorms, gamma, lambda);
}

/**
 * ncm_stats_dist_kernel_eval_sum1_gamma_lambda: (virtual eval_sum1_gamma_lambda)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: a #NcmVector
 * @weights: a #NcmVector
 * @lnnorm: a double
 * @gamma: (out): $\gamma$
 * @lambda: (out): $\lambda$
 *
 * Computes the weighted sum of kernels at $\chi^2=$@chi2,
 * $$ e^\gamma (1+\lambda) = \sum_i w_i\bar{K} (\chi^2_i) / u,$$
 * where $\gamma = \ln(w_a\bar{K} (\chi^2_a) / u)$ and $a$ labels
 * is the largest term of the sum.
 *
 */
void
ncm_stats_dist_kernel_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, gdouble *gamma, gdouble *lambda)
{
  NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_sum1_gamma_lambda (sdk, chi2, weights, lnnorm, gamma, lambda);
}

/**
 * ncm_stats_dist_kernel_sample: (virtual sample)
 * @sdk: a #NcmStatsDistKernel
 * @cov_decomp: Cholesky decomposition of the kernel covariance
 * @href: kernel bandwidth
 * @mu: kernel location vector
 * @y: output vector
 * @rng: a #NcmRNG
 *
 * Generates a random vector from the kernel distribution
 * using the covariance @cov_decomp, bandwidth @href and
 * location vector @mu. The result is stored in @y.
 *
 */
void
ncm_stats_dist_kernel_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->sample (sdk, cov_decomp, href, mu, y, rng);
}


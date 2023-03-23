/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * @short_description: An N-dimensional kernel used to compute the kernel density estimation function (KDE) in the #NcmStatsDist class.
 *
 * An N-dimensional kernel used to compute the kernel density estimation function (KDE) in the #NcmStatsDist class.
 *
 * This class provides the tools to generate a kernel function to be used in a kernel 
 * density estimation method. Below is a quick review of the kernel density estimation method 
 * and some properties of the kernel function, which are generalized for multidimensional problems. 
 * For further information, check [[Density Estimation for Statistics and Data Analysis, B.W. Silverman](https://www.routledge.com/Density-Estimation-for-Statistics-and-Data-Analysis/Silverman/p/book/9780412246203)].
 *
 * Starting with the uni-dimensional case, let $X_1,...,X_n$ be independent and identically 
 * distributed (iid) samples drawn from a distribution $f(x)$. The kernel density estimation of the function is 
 * \begin{align}
 * \tilde{f}(x) = \sum_{i=1}^{n}K\left(\frac{x-x_i}{h}\right)
 * ,\end{align}
 * where $K$ is the kernel function and $h$ is the bandwidth parameter. The kernel density 
 * estimator function must be close to the true density function $f(x)$, which can be tested 
 * by analyzing whether the estimator provides similar expected values as the function $f(x)$, 
 * that is, the function $\tilde{f}(x)$ must minimize the mean square error (MSE)
 * \begin{align}
 * \label{eqmse}
 * MSE_x(\tilde{f}) = E\left[\tilde{f}(x) - f(x)\right]^2
 * ,\end{align} 
 * where $E$ represents the expected value. This value depends on the choice of the kernel function, 
 * the data and the bandwidth. If the estimator $\tilde{f}(x)$ is close enough to the true function, 
 * it shall be used to generate samples that are distributed by $f(x)$.
 *
 * The kernel $K$ is a symmetric function that must satisfy
 * \begin{align}
 * &\int K(x)~dx = 1 
 * .\end{align}
 * Usually, the kernel function is a symmetric probability density function that is easy to sample from,
 * but it is totally under the user's control. Using simple kernels, such as the Gaussian kernel, makes 
 * the kernel density estimator method a better alternative to generate samples when the desired distribution is a complicated function.
 *
 * For the multidimensional case, given i.i.d d-dimensional sample points $X_1,.., X_n$ distributed by $f(x)$,
 * the multivariate kernel density estimator function $\tilde{f}(x)$ is given by
 * \begin{align}
 * \tilde{f}(x) = \frac{1}{h^d} \sum_{i=1}^n w_i K\left(\frac{x-x_i}{h}, \Sigma_i\right)
 * ,\end{align}
 * where $\Sigma_i$ is the covariance matrix of the $i$-th point (the kernels used in this library depend on the covariance matrix), 
 * $d$ is the dimension and $w_i$ is the weight attached to each kernel to find the minimal error in equation \eqref{eqmse}.
 *
 * The methods in this class define the type of kernel $K$, compute the bandwidth factor $h$, evaluate the kernel 
 * function at a given $d$-dimensional point $x$ or at a given vector of points $\vec{x}$, and, given the weights $w_i$, 
 * compute the kernel density estimation function $\tilde{f}(x)$.
 *
 * Besides the function ncm\_stats\_dist\_kernel\_get\_dim(), this class object only has virtual methods.
 * Therefore, to use this object, the user must initialize one of the child objects (#NcmStatsDistKernelGauss or #NcmStatsDistKernelST). 
 * Inside the child objects are the implemented functions, which must be defined for each specific type of kernel function. 
 * Check the childs documentations for more information. More information about how the algorithm should be implemented is described below:
 *		
 *		-This class is implemented in the #NcmStatsDist class, where the #NcmStatsDistKernel class shall define
 *		the type of kernel used in the interpolation function in #NcmStatsDist and how to compute values such as
 *		the weighted sum of the kernels, the bandwidth, and so on. Yet, the user may use these class objects
 *		to perform other kernel calculations, although some of the methods are not implemented outside the
 *		#NcmStatsDist class.
 *
 *		-This class does not possess the methods to compute the weights of each kernel. You may find this method in the
 *		#NcmStatsDist class. 
 *
 *		-Every child object of this class can be used either in the #NcmStatsDistKDE class or in the #NcmStatsDistVKDE class.
 *
 *
 *
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
 * Computes the weighted sum of kernels at $\chi^2=$@chi2 (the density estimator function),
 * $$ e^\gamma (1+\lambda) = \sum_i w_i\bar{K} (\chi^2_i) / u_i,$$
 * where $\gamma = \ln(w_a\bar{K} (\chi^2_a) / u_a)$ and $a$ labels
 * is the largest term of the sum. This function shall be used when
 * each kernel has a different normalization factor.
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
 * Computes the weighted sum of kernels at $\chi^2=$@chi2 (the density estimator function),
 * $$ e^\gamma (1+\lambda) = \sum_i w_i\bar{K} (\chi^2_i) / u,$$
 * where $\gamma = \ln(w_a\bar{K} (\chi^2_a) / u)$ and $a$ labels
 * is the largest term of the sum. This function shall be used when 
 * all the kernels have the same normalization factor.
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


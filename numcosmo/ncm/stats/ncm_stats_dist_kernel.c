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
 * NcmStatsDistKernel:
 *
 * Abstract kernel of the kernel mixture densities of #NcmStatsDist.
 *
 * A kernel is a symmetric density in $d$ dimensions with location $\mu$, scale matrix
 * $\Sigma$ and bandwidth $h$,
 * \begin{equation}
 * K(x) = \frac{\bar{K}(\chi^2)}{h^d\,u(\Sigma)}, \qquad
 * \chi^2 = \frac{(x - \mu)^T \Sigma^{-1} (x - \mu)}{h^2},
 * \end{equation}
 * where $\bar{K}$ is the unnormalized kernel (ncm_stats_dist_kernel_eval_unnorm()) and
 * $u(\Sigma)$ its normalization at $h = 1$ (ncm_stats_dist_kernel_get_lnnorm() returns
 * $\ln u$). The covariance of $K$ is $\kappa h^2 \Sigma$, with $\kappa$ given by
 * ncm_stats_dist_kernel_get_var_factor().
 *
 * Given points $x_1, \dots, x_n$ with weights $w_i$, #NcmStatsDist builds the mixture
 * \begin{equation}
 * \tilde{f}(x) = \sum_{i=1}^n w_i K_i(x),
 * \end{equation}
 * where $K_i$ has location $x_i$. #NcmStatsDist chooses the weights, the scale matrices
 * and the bandwidth, and adds the $d \ln h$ term to the normalization. The kernel
 * provides $\bar{K}$, $u$, $\kappa$, the rule-of-thumb bandwidth and samples.
 *
 * Apart from ncm_stats_dist_kernel_get_dim(), all methods are virtual. The
 * implementations are #NcmStatsDistKernelGauss and #NcmStatsDistKernelST.
 *
 * For background see [Density Estimation for Statistics and Data Analysis,
 * B.W. Silverman](https://www.routledge.com/Density-Estimation-for-Statistics-and-Data-Analysis/Silverman/p/book/9780412246203).
 *
 **/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist_kernel.h"
#include "ncm/stats/ncm_stats_vec.h"
#include "ncm/core/ncm_c.h"
#include "ncm/algebra/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#endif /* NUMCOSMO_GIR_SCAN */

#include "ncm/stats/ncm_stats_dist_kernel_private.h"

enum
{
  PROP_0,
  PROP_DIM,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDistKernel, ncm_stats_dist_kernel, G_TYPE_OBJECT)

static void
ncm_stats_dist_kernel_init (NcmStatsDistKernel *sdk)
{
  NcmStatsDistKernelPrivate * const self = ncm_stats_dist_kernel_get_instance_private (sdk);

  self->d = 0;
}

static void
_ncm_stats_dist_kernel_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernel *sdk = NCM_STATS_DIST_KERNEL (object);

  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL (object));

  switch (prop_id)
  {
    case PROP_DIM:
      NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->set_dim (sdk, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_kernel_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernel *sdk = NCM_STATS_DIST_KERNEL (object);

  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL (object));

  switch (prop_id)
  {
    case PROP_DIM:
      g_value_set_uint (value, ncm_stats_dist_kernel_get_dim (sdk));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
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
_ncm_stats_dist_kernel_get_var_factor (NcmStatsDistKernel *sdk)
{
  g_error ("method get_var_factor not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));

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
_ncm_stats_dist_kernel_eval_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *lnc, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  g_error ("method eval_gamma_lambda not implemented by %s.", G_OBJECT_TYPE_NAME (sdk));
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

  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dimension",
                                                      NULL,
                                                      "Kernel dimension",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd_class->set_dim           = &_ncm_stats_dist_kernel_set_dim;
  sd_class->get_dim           = &_ncm_stats_dist_kernel_get_dim;
  sd_class->get_rot_bandwidth = &_ncm_stats_dist_kernel_get_rot_bandwidth;
  sd_class->get_var_factor    = &_ncm_stats_dist_kernel_get_var_factor;
  sd_class->get_lnnorm        = &_ncm_stats_dist_kernel_get_lnnorm;
  sd_class->eval_unnorm       = &_ncm_stats_dist_kernel_eval_unnorm;
  sd_class->eval_unnorm_vec   = &_ncm_stats_dist_kernel_eval_unnorm_vec;
  sd_class->eval_gamma_lambda = &_ncm_stats_dist_kernel_eval_gamma_lambda;
  sd_class->sample            = &_ncm_stats_dist_kernel_sample;
}

static void
_ncm_stats_dist_kernel_set_dim (NcmStatsDistKernel *sdk, const guint dim)
{
  NcmStatsDistKernelPrivate * const self = ncm_stats_dist_kernel_get_instance_private (sdk);

  self->d = dim;
}

static guint
_ncm_stats_dist_kernel_get_dim (NcmStatsDistKernel *sdk)
{
  NcmStatsDistKernelPrivate * const self = ncm_stats_dist_kernel_get_instance_private (sdk);

  return self->d;
}

/**
 * ncm_stats_dist_kernel_ref:
 * @sdk: a #NcmStatsDistKernel
 *
 * Increases the reference count of @sdk by one.
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
 * Decreases the reference count of @sdk by one.
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
 * Decreases the reference count of *@sdk by one and sets *@sdk to NULL.
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
 * Returns: the kernel dimension $d$.
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
 * Computes the rule-of-thumb bandwidth $h$ for a mixture of @n kernels: the $h$
 * that minimizes the asymptotic mean integrated squared error when the estimated
 * density is the kernel itself with the scale matrix $\Sigma$ of the mixture. See
 * the implementations for the closed forms.
 *
 * Returns: the rule-of-thumb bandwidth.
 */
gdouble
ncm_stats_dist_kernel_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->get_rot_bandwidth (sdk, n);
}

/**
 * ncm_stats_dist_kernel_get_var_factor: (virtual get_var_factor)
 * @sdk: a #NcmStatsDistKernel
 *
 * Computes the factor $\kappa$ relating the kernel covariance to its scale matrix
 * $\Sigma$, that is $\mathrm{Cov} = \kappa \Sigma$. It is one for the Gaussian kernel
 * and $\nu / (\nu - 2)$ for the Student-t kernel with $\nu$ degrees of freedom,
 * which is infinite for $\nu \leq 2$ since such kernels have no covariance.
 *
 * Returns: the kernel variance factor $\kappa$.
 */
gdouble
ncm_stats_dist_kernel_get_var_factor (NcmStatsDistKernel *sdk)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->get_var_factor (sdk);
}

/**
 * ncm_stats_dist_kernel_get_lnnorm: (virtual get_lnnorm)
 * @sdk: a #NcmStatsDistKernel
 * @cov_decomp: upper-triangular Cholesky factor $U$ of the scale matrix, $\Sigma = U^T U$
 *
 * Computes $\ln u(\Sigma)$, the logarithm of the kernel normalization at $h = 1$.
 *
 * Returns: $\ln u(\Sigma)$.
 */
gdouble
ncm_stats_dist_kernel_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->get_lnnorm (sdk, cov_decomp);
}

/**
 * ncm_stats_dist_kernel_eval_unnorm: (virtual eval_unnorm)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: the scaled squared distance $\chi^2$
 *
 * Returns: the unnormalized kernel $\bar{K}(\chi^2)$ at $\chi^2 = $ @chi2.
 */
gdouble
ncm_stats_dist_kernel_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2)
{
  return NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_unnorm (sdk, chi2);
}

/**
 * ncm_stats_dist_kernel_eval_unnorm_vec: (virtual eval_unnorm_vec)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: a #NcmVector of $\chi^2$ values
 * @Ku: a #NcmVector of the same length
 *
 * Computes the unnormalized kernel $\bar{K}$ at every element of @chi2 and stores
 * the results in @Ku.
 *
 */
void
ncm_stats_dist_kernel_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku)
{
  NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_unnorm_vec (sdk, chi2, Ku);
}

/**
 * ncm_stats_dist_kernel_eval_gamma_lambda: (virtual eval_gamma_lambda)
 * @sdk: a #NcmStatsDistKernel
 * @chi2: a #NcmVector holding $\chi^2_i$, one entry per kernel
 * @lnc: a #NcmVector holding $\ln (w_i / u_i)$, one entry per kernel
 * @lnK: a #NcmVector that receives the logarithm of each term, $\ln (w_i \bar{K}(\chi^2_i) / u_i)$
 * @gamma: (out): $\gamma$
 * @lambda: (out): $\lambda$
 *
 * Computes the weighted sum of kernels (the mixture density at one point),
 * $$ e^\gamma (1+\lambda) = \sum_i w_i\bar{K} (\chi^2_i) / u_i,$$
 * where $\gamma = \ln(w_a\bar{K} (\chi^2_a) / u_a)$ and $a$ labels the largest term of
 * the sum. The three vectors must have the same length and unit stride.
 *
 * The weight and the normalization enter only through their ratio, so @lnc carries the
 * combination $\ln w_i - \ln u_i$ already formed. The caller builds it once per batch of
 * points instead of once per (point, kernel) pair, and a shared normalization is just the
 * same $u$ in every entry.
 *
 */
void
ncm_stats_dist_kernel_eval_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *lnc, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->eval_gamma_lambda (sdk, chi2, lnc, lnK, gamma, lambda);
}

/**
 * ncm_stats_dist_kernel_sample: (virtual sample)
 * @sdk: a #NcmStatsDistKernel
 * @cov_decomp: upper-triangular Cholesky factor $U$ of the scale matrix, $\Sigma = U^T U$
 * @href: kernel bandwidth $h$
 * @mu: kernel location $\mu$
 * @y: output vector
 * @rng: a #NcmRNG
 *
 * Draws a point from the kernel with location @mu, scale matrix $\Sigma$ and
 * bandwidth @href, and stores it in @y.
 *
 */
void
ncm_stats_dist_kernel_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng)
{
  NCM_STATS_DIST_KERNEL_GET_CLASS (sdk)->sample (sdk, cov_decomp, href, mu, y, rng);
}


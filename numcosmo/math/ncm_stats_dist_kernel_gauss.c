/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel_gauss.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel_gauss.c
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
 * SECTION:ncm_stats_dist_kernel_gauss
 * @title: NcmStatsDistKernelGauss
 * @short_description: An N-dimensional Gaussian kernel used to compute the kernel density estimation function (KDE) in the #NcmStatsDist class.
 * An N-dimensional Gaussian kernel used to compute the kernel density estimation function (KDE) in the #NcmStatsDist class.
 *
 * This object defines a multivariate Gaussian kernel to be used in the #NcmStatsDistKernel class. Also, this object implements
 * the virtual methods of the #NcmStatsDistKernel class. For more information, check the documentation of #NcmStatsDistKernel.
 * Below, there are some definitions of the multivariate Gaussian distribution. For more information, check  [[On Sampling from the Multivariate t Distribution, Marius Hofert](https://journal.r-project.org/archive/2013/RJ-2013-033/RJ-2013-033.pdf)].
 *
 * The multivariate Normal distribution has its stochastic representation as
 * \begin{align}
 * \textbf{X} &= \mu + A \textbf{Z}
 * ,\end{align}
 * where $\textbf{Z}= (Z_1,Z_2,...,Z_n)$ is a $n$-dimension random vector whose components are independent normal random variables.
 * $A$ is a $d \times n$ matrix and $\mu$ is a $d$-dimensional random vector that defines the mean of the distribution. The covariance
 * matrix is defined as $\Sigma =  AA^T$, such that the distribtuion of $\textbf{X}$ is uniquely defined by its covariance matrix and the
 * mean vector, that is,$ \textbf{X} \sim N(\mu, \Sigma)$.
 *
 * Assuming that $n=d$, the probability density function (pdf) of $\textbf{X}$ is
 * \begin{align}
 * \label{pdf}
 * f_{\textbf{X}(x)} = \frac{1}{(2\pi)^{\frac{d}{2}} \sqrt{det \Sigma}} \exp\left[-\frac{1}{2}(x-\mu)^T \Sigma^-1 (x-\mu)\right]
 * ,\end{align}
 * considering that the covariance matrix is positive definite and $x \in \mathbb{R^d}$. Also, the covariance matrix can be decomposed in its Cholesky
 * decomposition,
 * \begin{align}
 * \Sigma = LL^T
 * ,\end{align}
 * where $L$ is a triangular matrix with positive definite values. This decomposition can facilitate some computational calculations.
 *
 * This object uses the pdf given by equation \eqref{pdf} to define a Gaussian kernel, such that it can generate points distributed
 * by multivariate Gaussian distributions. The normal distribution is easy to sample from and therefore is commonly used as a kernel.
 *
 * The user must provide the following input value: @dim - ncm\_stats\_dist\_kernel\_gauss\_new(). Once this object is initialized,
 * the user can use the methods in the #NcmStatsDistKernel class with this object.
 *
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

typedef struct _NcmStatsDistKernelGaussPrivate
{
  gint place_holder;
} NcmStatsDistKernelGaussPrivate;

enum
{
  PROP_0,
};

struct _NcmStatsDistKernelGauss
{
  NcmStatsDistKernel parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistKernelGauss, ncm_stats_dist_kernel_gauss, NCM_TYPE_STATS_DIST_KERNEL)

static NcmStatsDistKernelPrivate *
ncm_stats_dist_kernel_get_instance_private (NcmStatsDistKernel * sdk)
{
  return g_type_instance_get_private ((GTypeInstance *) sdk, NCM_TYPE_STATS_DIST_KERNEL);
}

static void
ncm_stats_dist_kernel_gauss_init (NcmStatsDistKernelGauss *sdkg)
{
  NcmStatsDistKernelGaussPrivate * const self = ncm_stats_dist_kernel_gauss_get_instance_private (sdkg);

  self->place_holder = 0;
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
  /* NcmStatsDistKernelGauss *sdkg               = NCM_STATS_DIST_KERNEL_GAUSS (object); */
  /* NcmStatsDistKernelGaussPrivate * const self = sdkg->priv; */

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
static void _ncm_stats_dist_kernel_gauss_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, NcmVector *lnK, gdouble *gamma, gdouble *lambda);
static void _ncm_stats_dist_kernel_gauss_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, NcmVector *lnK, gdouble *gamma, gdouble *lambda);
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
  NcmStatsDistKernelPrivate * const pself = ncm_stats_dist_kernel_get_instance_private (sdk);

  return pow (4.0 / (n * (pself->d + 2.0)), 1.0 / (pself->d + 4.0));
}

static gdouble
_ncm_stats_dist_kernel_gauss_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  NcmStatsDistKernelPrivate * const pself = ncm_stats_dist_kernel_get_instance_private (sdk);

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
  guint i;

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
_ncm_stats_dist_kernel_gauss_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  /* NcmStatsDistKernelGaussPrivate * const self = NCM_STATS_DIST_KERNEL_GAUSS (sdk)->priv; */
  /*NcmStatsDistKernelPrivate * const pself  = sdk->priv;*/

  const guint n   = ncm_vector_len (chi2);
  gdouble lnt_max = GSL_NEGINF;
  guint i, i_max = 0;

  g_assert (n == ncm_vector_len (weights));
  g_assert (n == ncm_vector_len (lnnorms));
  g_assert (n == ncm_vector_len (lnK));
  g_assert (1 == ncm_vector_stride (chi2));
  g_assert (1 == ncm_vector_stride (weights));
  g_assert (1 == ncm_vector_stride (lnnorms));
  g_assert (1 == ncm_vector_stride (lnK));

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

    ncm_vector_fast_set (lnK, i, lnt_i);
  }

  lambda[0] = 0.0;

  for (i = 0; i < i_max; i++)
    lambda[0] += exp (ncm_vector_fast_get (lnK, i) - lnt_max);

  for (i = i_max + 1; i < n; i++)
    lambda[0] += exp (ncm_vector_fast_get (lnK, i) - lnt_max);

  gamma[0] = lnt_max;
}

static void
_ncm_stats_dist_kernel_gauss_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  /* NcmStatsDistKernelGaussPrivate * const self = NCM_STATS_DIST_KERNEL_GAUSS (sdk)->priv; */
  /* NcmStatsDistKernelPrivate * const pself  = sdk->priv; */

  const guint n   = ncm_vector_len (chi2);
  gdouble lnt_max = GSL_NEGINF;
  guint i, i_max = 0;

  g_assert (n == ncm_vector_len (weights));
  g_assert (n == ncm_vector_len (lnK));
  g_assert (1 == ncm_vector_stride (chi2));
  g_assert (1 == ncm_vector_stride (weights));
  g_assert (1 == ncm_vector_stride (lnK));

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

    ncm_vector_fast_set (lnK, i, lnt_i);
  }

  lambda[0] = 0.0;

  for (i = 0; i < i_max; i++)
    lambda[0] += exp (ncm_vector_fast_get (lnK, i) - lnt_max);

  for (i = i_max + 1; i < n; i++)
    lambda[0] += exp (ncm_vector_fast_get (lnK, i) - lnt_max);

  gamma[0] = lnt_max - lnnorm;
}

static void
_ncm_stats_dist_kernel_gauss_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistKernelPrivate * const pself = ncm_stats_dist_kernel_get_instance_private (sdk);
  gint ret;
  guint i;

  for (i = 0; i < pself->d; i++)
  {
    const gdouble u_i = ncm_rng_ugaussian_gen (rng);

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


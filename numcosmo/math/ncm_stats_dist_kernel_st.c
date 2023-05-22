/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel_st.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel_st.c
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
 * SECTION:ncm_stats_dist_kernel_st
 * @title: NcmStatsDistKernelST
 * @short_description: An N-dimensional Student's t kernel used to compute the kernel density estimation function (KDE) in the #NcmStatsDist class.
 *
 * An N-dimensional Student's t kernel used to compute the kernel density estimation function (KDE) in the #NcmStatsDist class.
 *
 * This object defines a multivariate Student's t kernel to be used in the #NcmStatsDistKernel class. Also, this object implements
 * the virtual methods of the #NcmStatsDistKernel class. For more information about the class, check the documentation of #NcmStatsDistKernel.
 * Below, there are some definitions of the multivariate Student t distribution.
 * For more information, check  [[On Sampling from the Multivariate t
 * Distribution, Marius Hofert](https://journal.r-project.org/archive/2013/RJ-2013-033/RJ-2013-033.pdf)].
 *
 * The multivariate t distribution with $\nu$ degrees of freedom has its stochastic representation as
 * \begin{align}
 * \label{st}
 * \textbf{X} &= \mu + \sqrt{W} A \textbf{Z}
 * ,\end{align}
 * where $\textbf{Z}=(Z_1,Z_2,...,Z_n)$ is a $n$-dimension random vector whose components are independent normal random variables.
 * $A$ is a $d \times n$ matrix, $\mu$ is a $d$-dimensional random vector that defines the mean of the distribution and
 * $W=\frac{\nu}{\chi^2}$, being $\chi^2$ a random variable following a chi-squared distribution with $\nu > 0$ degrees of freedom.
 * The covariance matrix is defined as $\Sigma =  AA^T$, such that the distribution of $\textbf{X}$ is uniquely defined by its
 * covariance matrix and the mean vector, that is, $\textbf{X} \sim t(\mu, \Sigma)$.
 *
 * Assuming that $n=d$, the probability density function (pdf) of $\textbf{X}$ is
 * \begin{align}
 * \label{pdfst}
 * f_{\textbf{X}(x)} = \frac{\Gamma\left(\frac{\nu + d}{2}\right)}{\Gamma\left(\frac{\nu}{2}\right)(2\pi)^{\frac{d}{2}} \sqrt{det \Sigma}} \left[1+\frac{(x-\mu)^T \Sigma^-1 (x-\mu)}{\nu}\right]^{-\frac{\nu + d}{2}}
 * ,\end{align}
 * considering that the covariance matrix is positive definite and $x \in \mathbb{R^d}$. Also, the covariance matrix can be
 * decomposed in its Cholesky decomposition,
 * \begin{align}
 * \Sigma = LL^t
 * ,\end{align}
 * where $L$ is a triangular matrix with positive definite values. This decomposition can facilitate some computational calculations.
 *
 * The $\sqrt{W}$ factor makes the multivariate t distribution more flexible than the multivariate Gaussian distribution,
 * especially on its tails. Therefore, for problems that require a smoother function, the multivariate t kernell shall be used.
 * Also, as seen in equation \eqref{st}, the Student's t distribtuion can be generated using normal random variables, which makes the
 * distribution easier to be generated. For the case $\nu \rightarrow \infty$, the multivariate t distribution becomes the Gaussian
 * distribution.
 *
 * This object uses the pdf given by equation \eqref{pdfst} to define a Student's t kernel, such that it can generate points distributed
 * by multivariate t distributions.
 *
 * The user must provide the following input value: @dim - ncm\_stats\_dist\_kernel\_st\_new(). Once this object is initialized,
 * the user can use the methods in the #NcmStatsDistKernel class with this object.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "gsl/gsl_sf_result.h"

#include "math/ncm_stats_dist_kernel_st.h"
#include "math/ncm_stats_vec.h"
#include "math/ncm_c.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_gamma.h>
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "math/ncm_stats_dist_kernel_private.h"

struct _NcmStatsDistKernelSTPrivate
{
  gdouble nu;
};

enum
{
  PROP_0,
  PROP_NU,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistKernelST, ncm_stats_dist_kernel_st, NCM_TYPE_STATS_DIST_KERNEL);

static void
ncm_stats_dist_kernel_st_init (NcmStatsDistKernelST *sdkst)
{
  NcmStatsDistKernelSTPrivate * const self = sdkst->priv = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  self->nu = 0.0;
}

static void
_ncm_stats_dist_kernel_st_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernelST *sdkst = NCM_STATS_DIST_KERNEL_ST (object);

  /*NcmStatsDistKernelSTPrivate * const self = sdkst->priv;*/
  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL_ST (object));

  switch (prop_id)
  {
    case PROP_NU:
      ncm_stats_dist_kernel_st_set_nu (sdkst, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_kernel_st_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernelST *sdkst = NCM_STATS_DIST_KERNEL_ST (object);

  /*NcmStatsDistKernelSTPrivate * const self = sdkst->priv;*/

  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL_ST (object));

  switch (prop_id)
  {
    case PROP_NU:
      g_value_set_double (value, ncm_stats_dist_kernel_st_get_nu (sdkst));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_kernel_st_dispose (GObject *object)
{
  /*NcmStatsDistKernelST *sdkst               = NCM_STATS_DIST_KERNEL_ST (object);*/
  /*NcmStatsDistKernelSTPrivate * const self = sdkst->priv;*/

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kernel_st_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_kernel_st_finalize (GObject *object)
{
  /* NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (object); */
  /* NcmStatsDistKernelSTPrivate * const self = sdkst->priv; */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kernel_st_parent_class)->finalize (object);
}

static gdouble _ncm_stats_dist_kernel_st_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n);
static gdouble _ncm_stats_dist_kernel_st_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp);
static gdouble _ncm_stats_dist_kernel_st_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2);
static void _ncm_stats_dist_kernel_st_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku);
static void _ncm_stats_dist_kernel_st_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, NcmVector *lnK, gdouble *gamma, gdouble *lambda);
static void _ncm_stats_dist_kernel_st_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, NcmVector *lnK, gdouble *gamma, gdouble *lambda);
static void _ncm_stats_dist_kernel_st_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng);

static void
ncm_stats_dist_kernel_st_class_init (NcmStatsDistKernelSTClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmStatsDistKernelClass *sdk_class = NCM_STATS_DIST_KERNEL_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_kernel_st_set_property;
  object_class->get_property = &_ncm_stats_dist_kernel_st_get_property;
  object_class->dispose      = &_ncm_stats_dist_kernel_st_dispose;
  object_class->finalize     = &_ncm_stats_dist_kernel_st_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NU,
                                   g_param_spec_double ("nu",
                                                        NULL,
                                                        "nu value of the function",
                                                        1.0, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sdk_class->get_rot_bandwidth      = &_ncm_stats_dist_kernel_st_get_rot_bandwidth;
  sdk_class->get_lnnorm             = &_ncm_stats_dist_kernel_st_get_lnnorm;
  sdk_class->eval_unnorm            = &_ncm_stats_dist_kernel_st_eval_unnorm;
  sdk_class->eval_unnorm_vec        = &_ncm_stats_dist_kernel_st_eval_unnorm_vec;
  sdk_class->eval_sum0_gamma_lambda = &_ncm_stats_dist_kernel_st_eval_sum0_gamma_lambda;
  sdk_class->eval_sum1_gamma_lambda = &_ncm_stats_dist_kernel_st_eval_sum1_gamma_lambda;
  sdk_class->sample                 = &_ncm_stats_dist_kernel_st_sample;
}

static gdouble
_ncm_stats_dist_kernel_st_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = sdkst->priv;

  const guint d    = ncm_stats_dist_kernel_get_dim (sdk);
  const gdouble nu = (self->nu >= 3.0) ? self->nu : 3.0;

  return pow (
    16.0 * gsl_pow_2 (nu - 2) * (1.0 + d + nu) * (3.0 + d + nu) /
    ((2.0 + d) * (d + nu) * (2.0 + d + nu) * (d + 2.0 * nu) * (2.0 + d + 2.0 * nu) * n),
    1.0 / (d + 4.0));
}

static gdouble
_ncm_stats_dist_kernel_st_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = sdkst->priv;

  const guint d             = ncm_stats_dist_kernel_get_dim (sdk);
  const gdouble lg_lnnorm   = lgamma (self->nu / 2.0) - lgamma ((self->nu + d) / 2.0);
  const gdouble chol_lnnorm = 0.5 * ncm_matrix_cholesky_lndet (cov_decomp);
  const gdouble nc_lnnorm   = (d / 2.0) * (ncm_c_lnpi () + log (self->nu));

  return lg_lnnorm + nc_lnnorm + chol_lnnorm;
}

static gdouble
_ncm_stats_dist_kernel_st_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2)
{
  NcmStatsDistKernelSTPrivate * const self = NCM_STATS_DIST_KERNEL_ST (sdk)->priv;
  NcmStatsDistKernelPrivate * const pself  = sdk->priv;

  return pow (1.0 + chi2 / self->nu, -0.5 * (self->nu + pself->d));
}

static void
_ncm_stats_dist_kernel_st_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku)
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
      const gdouble Ku_i   = _ncm_stats_dist_kernel_st_eval_unnorm (sdk, chi2_i);

      ncm_vector_fast_set (Ku, i, Ku_i);
    }
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      const gdouble chi2_i = ncm_vector_get (chi2, i);
      const gdouble Ku_i   = _ncm_stats_dist_kernel_st_eval_unnorm (sdk, chi2_i);

      ncm_vector_set (Ku, i, Ku_i);
    }
  }
}

static void
_ncm_stats_dist_kernel_st_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  NcmStatsDistKernelSTPrivate * const self = NCM_STATS_DIST_KERNEL_ST (sdk)->priv;
  NcmStatsDistKernelPrivate * const pself  = sdk->priv;

  const guint n = ncm_vector_len (chi2);
  const gdouble kappa = -0.5 * (self->nu + pself->d);
  gdouble lnt_max = GSL_NEGINF;
  gint i, i_max = -1;

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

    const gdouble lnt_i = kappa * log1p (chi2_i / self->nu) - lnu_i + log (w_i);

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
_ncm_stats_dist_kernel_st_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  NcmStatsDistKernelSTPrivate * const self = NCM_STATS_DIST_KERNEL_ST (sdk)->priv;
  NcmStatsDistKernelPrivate * const pself  = sdk->priv;

  const guint n = ncm_vector_len (chi2);
  const gdouble kappa = -0.5 * (self->nu + pself->d);
  gdouble lnt_max = GSL_NEGINF;
  gint i, i_max = -1;

  g_assert_cmpuint (n, ==, ncm_vector_len (weights));
  g_assert_cmpuint (n, ==, ncm_vector_len (lnK));
  g_assert_cmpuint (1, ==, ncm_vector_stride (chi2));
  g_assert_cmpuint (1, ==, ncm_vector_stride (weights));
  g_assert_cmpuint (1, ==, ncm_vector_stride (lnK));

  for (i = 0; i < n; i++)
  {
    const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
    const gdouble w_i    = ncm_vector_fast_get (weights, i);

    const gdouble lnt_i = kappa * log1p (chi2_i / self->nu) + log (w_i);

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
_ncm_stats_dist_kernel_st_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = sdkst->priv;
  NcmStatsDistKernelPrivate * const pself  = sdk->priv;
  gdouble chi_scale;
  gint i, ret;

  for (i = 0; i < pself->d; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);

    ncm_vector_set (x, i, u_i * href);
  }

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (x));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_kernel_st_sample", ret);

  chi_scale = sqrt (self->nu / gsl_ran_chisq (rng->r, self->nu));

  ncm_vector_scale (x, chi_scale);
  ncm_vector_add (x, mu);
}

/**
 * ncm_stats_dist_kernel_st_new:
 * @dim: sample space dimension
 * @nu: Student-t parameter $\nu$
 *
 * Creates a new #NcmStatsDistKernelST object with sample dimension @dim
 * and $\nu$ = @nu.
 *
 * Returns: a new #NcmStatsDistKernelST.
 */
NcmStatsDistKernelST *
ncm_stats_dist_kernel_st_new (const guint dim, const gdouble nu)
{
  NcmStatsDistKernelST *sdkst = g_object_new (NCM_TYPE_STATS_DIST_KERNEL_ST,
                                              "dimension", dim,
                                              "nu",        nu,
                                              NULL);

  return sdkst;
}

/**
 * ncm_stats_dist_kernel_st_ref:
 * @sdkst: a #NcmStatsDistKernelST
 *
 * Increase the reference of @stats_dist_kernel_st by one.
 *
 * Returns: (transfer full): @stats_dist_kernel_st.
 */
NcmStatsDistKernelST *
ncm_stats_dist_kernel_st_ref (NcmStatsDistKernelST *sdkst)
{
  return g_object_ref (sdkst);
}

/**
 * ncm_stats_dist_kernel_st_free:
 * @sdkst: a #NcmStatsDistKernelST
 *
 * Decrease the reference count of @stats_dist_kernel_st by one.
 *
 */
void
ncm_stats_dist_kernel_st_free (NcmStatsDistKernelST *sdkst)
{
  g_object_unref (sdkst);
}

/**
 * ncm_stats_dist_kernel_st_clear:
 * @sdkst: a #NcmStatsDistKernelST
 *
 * Decrease the reference count of @stats_dist_kernel_st by one, and sets the pointer *@stats_dist_kernel_st to
 * NULL.
 *
 */
void
ncm_stats_dist_kernel_st_clear (NcmStatsDistKernelST **sdkst)
{
  g_clear_object (sdkst);
}

/**
 * ncm_stats_dist_kernel_st_set_nu:
 * @sdkst: a #NcmStatsDistKernelST
 * @nu: the over-smooth factor
 *
 * Sets the over-smooth factor to @nu.
 *
 */
void
ncm_stats_dist_kernel_st_set_nu (NcmStatsDistKernelST *sdkst, const gdouble nu)
{
  NcmStatsDistKernelSTPrivate * const self = sdkst->priv;

  self->nu = nu;
}

/**
 * ncm_stats_dist_kernel_st_get_nu:
 * @sdkst: a #NcmStatsDistKernelST
 *
 * Returns: the over-smooth factor.
 */
gdouble
ncm_stats_dist_kernel_st_get_nu (NcmStatsDistKernelST *sdkst)
{
  NcmStatsDistKernelSTPrivate * const self = sdkst->priv;

  return self->nu;
}


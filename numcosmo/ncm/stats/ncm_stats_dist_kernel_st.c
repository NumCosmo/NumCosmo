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
 * NcmStatsDistKernelST:
 *
 * Multivariate Student's t kernel for #NcmStatsDist.
 *
 * The kernel of #NcmStatsDistKernel with $\nu$ degrees of freedom,
 * \begin{equation}
 * \bar{K}(\chi^2) = \left(1 + \frac{\chi^2}{\nu}\right)^{-(\nu + d)/2}, \qquad
 * u(\Sigma) = \frac{\Gamma(\nu/2)\,(\nu\pi)^{d/2}}{\Gamma\left((\nu + d)/2\right)}\sqrt{\det\Sigma},
 * \end{equation}
 * that is, the density of the multivariate t distribution with location $\mu$ and scale
 * matrix $h^2\Sigma$. For $\nu > 2$ its covariance is $\kappa h^2\Sigma$ with
 * $\kappa = \nu/(\nu - 2)$; for $\nu \leq 2$ it has none. A sample is
 * $x = \mu + h\sqrt{\nu/W}\,U^T z$, where $z$ holds $d$ independent standard normal
 * variables, $W$ is a chi-squared variable with $\nu$ degrees of freedom and
 * $\Sigma = U^T U$; see [On Sampling from the Multivariate t Distribution, Marius
 * Hofert](https://journal.r-project.org/archive/2013/RJ-2013-033/RJ-2013-033.pdf). As
 * $\nu \to \infty$ the kernel tends to #NcmStatsDistKernelGauss; for finite $\nu$ it
 * decays as a power of $\chi^2$.
 *
 * The rule-of-thumb bandwidth is
 * \begin{equation}
 * h = \left[\frac{16 (\nu - 2)^2 (1 + d + \nu)(3 + d + \nu)}
 * {(2 + d)(d + \nu)(2 + d + \nu)(d + 2\nu)(2 + d + 2\nu)\,n}\right]^{1/(d + 4)},
 * \end{equation}
 * which minimizes the asymptotic mean integrated squared error for $n$ points drawn from
 * a t density with $\nu$ degrees of freedom and scale matrix $\Sigma$. It is evaluated
 * at $\nu = 3$ when $\nu < 3$, and tends to the Gaussian rule as $\nu \to \infty$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist_kernel_st.h"
#include "ncm/stats/ncm_stats_vec.h"
#include "ncm/core/ncm_c.h"
#include "ncm/algebra/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

#include "ncm/stats/ncm_stats_dist_kernel_private.h"

typedef struct _NcmStatsDistKernelSTPrivate
{
  gdouble nu;
} NcmStatsDistKernelSTPrivate;

enum
{
  PROP_0,
  PROP_NU,
};

struct _NcmStatsDistKernelST
{
  NcmStatsDistKernel parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistKernelST, ncm_stats_dist_kernel_st, NCM_TYPE_STATS_DIST_KERNEL)

static NcmStatsDistKernelPrivate *
ncm_stats_dist_kernel_get_instance_private (NcmStatsDistKernel * sdk)
{
  return g_type_instance_get_private ((GTypeInstance *) sdk, NCM_TYPE_STATS_DIST_KERNEL);
}

static void
ncm_stats_dist_kernel_st_init (NcmStatsDistKernelST *sdkst)
{
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  self->nu = 0.0;
}

static void
_ncm_stats_dist_kernel_st_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernelST *sdkst = NCM_STATS_DIST_KERNEL_ST (object);

  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL_ST (object));

  switch (prop_id)
  {
    case PROP_NU:
      ncm_stats_dist_kernel_st_set_nu (sdkst, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_kernel_st_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKernelST *sdkst = NCM_STATS_DIST_KERNEL_ST (object);

  g_return_if_fail (NCM_IS_STATS_DIST_KERNEL_ST (object));

  switch (prop_id)
  {
    case PROP_NU:
      g_value_set_double (value, ncm_stats_dist_kernel_st_get_nu (sdkst));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _ncm_stats_dist_kernel_st_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n);
static gdouble _ncm_stats_dist_kernel_st_get_var_factor (NcmStatsDistKernel *sdk);
static gdouble _ncm_stats_dist_kernel_st_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp);
static gdouble _ncm_stats_dist_kernel_st_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2);
static void _ncm_stats_dist_kernel_st_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku);
static void _ncm_stats_dist_kernel_st_eval_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *lnc, NcmVector *lnK, gdouble *gamma, gdouble *lambda);
static void _ncm_stats_dist_kernel_st_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng);

static void
ncm_stats_dist_kernel_st_class_init (NcmStatsDistKernelSTClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmStatsDistKernelClass *sdk_class = NCM_STATS_DIST_KERNEL_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_kernel_st_set_property;
  object_class->get_property = &_ncm_stats_dist_kernel_st_get_property;

  g_object_class_install_property (object_class,
                                   PROP_NU,
                                   g_param_spec_double ("nu",
                                                        NULL,
                                                        "Degrees of freedom",
                                                        1.0, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sdk_class->get_rot_bandwidth = &_ncm_stats_dist_kernel_st_get_rot_bandwidth;
  sdk_class->get_var_factor    = &_ncm_stats_dist_kernel_st_get_var_factor;
  sdk_class->get_lnnorm        = &_ncm_stats_dist_kernel_st_get_lnnorm;
  sdk_class->eval_unnorm       = &_ncm_stats_dist_kernel_st_eval_unnorm;
  sdk_class->eval_unnorm_vec   = &_ncm_stats_dist_kernel_st_eval_unnorm_vec;
  sdk_class->eval_gamma_lambda = &_ncm_stats_dist_kernel_st_eval_gamma_lambda;
  sdk_class->sample            = &_ncm_stats_dist_kernel_st_sample;
}

static gdouble
_ncm_stats_dist_kernel_st_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  const guint d    = ncm_stats_dist_kernel_get_dim (sdk);
  const gdouble nu = (self->nu >= 3.0) ? self->nu : 3.0;

  return pow (
    16.0 * gsl_pow_2 (nu - 2) * (1.0 + d + nu) * (3.0 + d + nu) /
    ((2.0 + d) * (d + nu) * (2.0 + d + nu) * (d + 2.0 * nu) * (2.0 + d + 2.0 * nu) * n),
    1.0 / (d + 4.0));
}

static gdouble
_ncm_stats_dist_kernel_st_get_var_factor (NcmStatsDistKernel *sdk)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  /* The multivariate Student-t covariance is nu / (nu - 2) times its scale matrix;
   * for nu <= 2 the second moment does not exist. */
  if (self->nu <= 2.0)
    return GSL_POSINF;

  return self->nu / (self->nu - 2.0);
}

static gdouble
_ncm_stats_dist_kernel_st_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  const guint d             = ncm_stats_dist_kernel_get_dim (sdk);
  const gdouble lg_lnnorm   = lgamma (self->nu / 2.0) - lgamma ((self->nu + d) / 2.0);
  const gdouble chol_lnnorm = 0.5 * ncm_matrix_cholesky_lndet (cov_decomp);
  const gdouble nc_lnnorm   = (d / 2.0) * (ncm_c_lnpi () + log (self->nu));

  return lg_lnnorm + nc_lnnorm + chol_lnnorm;
}

static gdouble
_ncm_stats_dist_kernel_st_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);
  NcmStatsDistKernelPrivate * const pself  = ncm_stats_dist_kernel_get_instance_private (sdk);

  return pow (1.0 + chi2 / self->nu, -0.5 * (self->nu + pself->d));
}

static void
_ncm_stats_dist_kernel_st_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku)
{
  const guint n = ncm_vector_len (chi2);
  guint i;

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
_ncm_stats_dist_kernel_st_eval_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *lnc, NcmVector *lnK, gdouble *gamma, gdouble *lambda)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);
  NcmStatsDistKernelPrivate * const pself  = ncm_stats_dist_kernel_get_instance_private (sdk);

  const gdouble kappa = -0.5 * (self->nu + pself->d);
  const guint n       = ncm_vector_len (chi2);
  gdouble lnt_max     = GSL_NEGINF;
  guint i, i_max = 0;

  g_assert_cmpuint (n, ==, ncm_vector_len (lnc));
  g_assert_cmpuint (n, ==, ncm_vector_len (lnK));
  g_assert_cmpuint (1, ==, ncm_vector_stride (chi2));
  g_assert_cmpuint (1, ==, ncm_vector_stride (lnc));
  g_assert_cmpuint (1, ==, ncm_vector_stride (lnK));

  for (i = 0; i < n; i++)
  {
    const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
    const gdouble lnc_i  = ncm_vector_fast_get (lnc, i);

    const gdouble lnt_i = kappa * log1p (chi2_i / self->nu) + lnc_i;

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
_ncm_stats_dist_kernel_st_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistKernelST *sdkst              = NCM_STATS_DIST_KERNEL_ST (sdk);
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);
  NcmStatsDistKernelPrivate * const pself  = ncm_stats_dist_kernel_get_instance_private (sdk);
  gdouble chi_scale;
  guint i;

  for (i = 0; i < pself->d; i++)
  {
    const gdouble u_i = ncm_rng_ugaussian_gen (rng);

    ncm_vector_set (x, i, u_i * href);
  }

  /* x <- U^T x, the lower factor L = U^T applied to a standard normal draw. */
  ncm_matrix_dtrmv (cov_decomp, 'U', 'T', x);

  chi_scale = sqrt (self->nu / ncm_rng_chisq_gen (rng, self->nu));

  ncm_vector_scale (x, chi_scale);
  ncm_vector_add (x, mu);
}

/**
 * ncm_stats_dist_kernel_st_new:
 * @dim: sample space dimension
 * @nu: degrees of freedom $\nu$
 *
 * Creates a new #NcmStatsDistKernelST of dimension @dim with @nu degrees of freedom.
 *
 * Returns: (transfer full): a new #NcmStatsDistKernelST.
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
 * Increases the reference count of @sdkst by one.
 *
 * Returns: (transfer full): @sdkst.
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
 * Decreases the reference count of @sdkst by one.
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
 * Decreases the reference count of *@sdkst by one and sets *@sdkst to NULL.
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
 * @nu: degrees of freedom $\nu$
 *
 * Sets the degrees of freedom to @nu.
 *
 */
void
ncm_stats_dist_kernel_st_set_nu (NcmStatsDistKernelST *sdkst, const gdouble nu)
{
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  self->nu = nu;
}

/**
 * ncm_stats_dist_kernel_st_get_nu:
 * @sdkst: a #NcmStatsDistKernelST
 *
 * Returns: the degrees of freedom $\nu$.
 */
gdouble
ncm_stats_dist_kernel_st_get_nu (NcmStatsDistKernelST *sdkst)
{
  NcmStatsDistKernelSTPrivate * const self = ncm_stats_dist_kernel_st_get_instance_private (sdkst);

  return self->nu;
}


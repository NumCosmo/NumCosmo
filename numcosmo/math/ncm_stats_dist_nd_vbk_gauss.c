/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd_vbk_gauss.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_vbk_gauss.c
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
 * SECTION:ncm_stats_dist_nd_vbk_gauss
 * @title: NcmStatsDistNdVbkGauss
 * @short_description: An N dimensional probability distributions using gaussian KDE
 *
 * An arbitrary N dimensional probability distribution using gaussian KDE.
 *
 * 
 * This object provides the tools to perform a radial basis interpolation
 * in a multidimensional function, using a multivariate gaussian distribution.
 * such that we are able to sample the original functio from this new interpolation function.
 * For more informations about radial basis interpolation,   
 * check #NcmStatsDistNd.
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
 * Once this object is initialized, we may use the methods in the #NcmStatsDistNdVbkclass to perform the interpolation 
 * and to generate a sample from the interpolated polinomal function.    
 *
 * The user must provide the following input values: $p$ - ncm_stats_dist_nd_vbk_gauss_new(),
 * cv_type - ncm_stats_dist_nd_vbk_gauss_new().
 **/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_nd_vbk_gauss.h"
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

struct _NcmStatsDistNdVbkGaussPrivate
{
  GArray *eval_v;
};

enum
{
  PROP_0,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistNdVbkGauss, ncm_stats_dist_nd_vbk_gauss, NCM_TYPE_STATS_DIST_ND_VBK);

static void
ncm_stats_dist_nd_vbk_gauss_init (NcmStatsDistNdVbkGauss *dndg)
{
  NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv = ncm_stats_dist_nd_vbk_gauss_get_instance_private (dndg);

  self->eval_v = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void
_ncm_stats_dist_nd_vbk_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND_VBK_GAUSS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_vbk_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;*/

  g_return_if_fail (NCM_IS_STATS_DIST_ND_VBK_GAUSS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_vbk_gauss_dispose (GObject *object)
{
  NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_VBK_GAUSS (object);
  NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;

  g_clear_pointer (&self->eval_v, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_vbk_gauss_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_vbk_gauss_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_vbk_gauss_parent_class)->finalize (object);
}

static gdouble _ncm_stats_dist_nd_vbk_gauss_get_rot_bandwidth (NcmStatsDistNdVbk *dnd, const guint d, const gdouble n);
static gdouble _ncm_stats_dist_nd_vbk_gauss_get_kernel_lnnorm (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href);
static void _ncm_stats_dist_nd_vbk_gauss_prepare_IM (NcmStatsDistNdVbk *dnd, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, NcmMatrix *IM);

static gdouble _ncm_stats_dist_nd_vbk_gauss_eval (NcmStatsDistNdVbk *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href);
static gdouble _ncm_stats_dist_nd_vbk_gauss_eval_m2lnp (NcmStatsDistNdVbk *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href);
static void _ncm_stats_dist_nd_vbk_gauss_kernel_sample (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, const NcmVector *href, NcmRNG *rng);
static gdouble _ncm_stats_dist_nd_vbk_gauss_kernel_eval_m2lnp (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const NcmVector *href);

static void
ncm_stats_dist_nd_vbk_gauss_class_init (NcmStatsDistNdVbkGaussClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmStatsDistNdVbkClass *dnd_class  = NCM_STATS_DIST_ND_VBK_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_nd_vbk_gauss_set_property;
  object_class->get_property = &_ncm_stats_dist_nd_vbk_gauss_get_property;
  object_class->dispose      = &_ncm_stats_dist_nd_vbk_gauss_dispose;
  object_class->finalize     = &_ncm_stats_dist_nd_vbk_gauss_finalize;
  
  dnd_class->get_rot_bandwidth = &_ncm_stats_dist_nd_vbk_gauss_get_rot_bandwidth;
  dnd_class->get_kernel_lnnorm = &_ncm_stats_dist_nd_vbk_gauss_get_kernel_lnnorm;
  dnd_class->prepare_IM        = &_ncm_stats_dist_nd_vbk_gauss_prepare_IM;
  dnd_class->eval              = &_ncm_stats_dist_nd_vbk_gauss_eval;
  dnd_class->eval_m2lnp        = &_ncm_stats_dist_nd_vbk_gauss_eval_m2lnp;
  dnd_class->kernel_sample     = &_ncm_stats_dist_nd_vbk_gauss_kernel_sample;
  dnd_class->kernel_eval_m2lnp = &_ncm_stats_dist_nd_vbk_gauss_kernel_eval_m2lnp;
}

static gdouble
_ncm_stats_dist_nd_vbk_gauss_f (NcmStatsDistNdVbkGaussPrivate * const self, const gdouble chi2)
{
  return exp (- 0.5 * chi2);
}

static gdouble
_ncm_stats_dist_nd_vbk_gauss_get_rot_bandwidth (NcmStatsDistNdVbk *dnd, const guint d, const gdouble n)
{
  return pow (4.0 / (n * (d + 2.0)), 1.0 / (d + 4.0));
}

static gdouble
_ncm_stats_dist_nd_vbk_gauss_get_kernel_lnnorm (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href)
{
  register gdouble det_href = 1.0;
  gint i;
  for (i = 0; i < d; i++)
    det_href *= ncm_vector_fast_get (href, i);

  return 0.5 * (d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (cov_decomp)) + log (det_href);
}

static void 
_ncm_stats_dist_nd_vbk_gauss_prepare_IM (NcmStatsDistNdVbk *dnd, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, NcmMatrix *IM)
{
  NcmStatsDistNdVbkGaussPrivate * const self = NCM_STATS_DIST_ND_VBK_GAUSS (dnd)->priv;
  gint i;

  for (i = 0; i < n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (invUsample, i);
    gint j;

    ncm_matrix_set (IM, i, i, 1.0);

    for (j = i + 1; j < n; j++)
    {
      NcmVector *row_j = g_ptr_array_index (invUsample, j);
      gdouble chi2_ij = 0.0;
      gdouble p_ij;
      gint k;

      for (k = 0; k < d; k++)
      {
        chi2_ij += gsl_pow_2 ((ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (row_j, k)) / ncm_vector_fast_get (href, k));
      }

      p_ij = _ncm_stats_dist_nd_vbk_gauss_f (self, chi2_ij);
      
      ncm_matrix_set (IM, i, j, p_ij);
      ncm_matrix_set (IM, j, i, p_ij);
    }
  }
}

static gdouble 
_ncm_stats_dist_nd_vbk_gauss_eval (NcmStatsDistNdVbk *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href)
{
  NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_VBK_GAUSS (dnd);
  NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i;
  
  for (i = 0; i < n; i++)
  {
    NcmVector *row_i  = g_ptr_array_index (invUsample, i);
    const gdouble w_i = ncm_vector_get (weights, i);
    gdouble e_i, t, chi2_i = 0.0;
    gint k;

    for (k = 0; k < d; k++)
    {
      chi2_i += gsl_pow_2 ((ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (invUy, k)) / ncm_vector_fast_get (href, k));
    }

    e_i  = w_i * _ncm_stats_dist_nd_vbk_gauss_f (self, chi2_i);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
  }

  return s;
}

static gdouble
_ncm_stats_dist_nd_vbk_gauss_eval_m2lnp (NcmStatsDistNdVbk *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href)
{
  NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_VBK_GAUSS (dnd);
  NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i;

  for (i = 0; i < n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (invUsample, i);
    const gdouble w_i = ncm_vector_get (weights, i);
    gdouble e_i, t, chi2_i = 0.0;
    gint k;

    for (k = 0; k < d; k++)
    {
      chi2_i += gsl_pow_2 ((ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (invUy, k)) / ncm_vector_fast_get (href, k));
    }

    e_i  = w_i * _ncm_stats_dist_nd_vbk_gauss_f (self, chi2_i);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
  }

  return -2.0 * log (s);
}

static void
_ncm_stats_dist_nd_vbk_gauss_kernel_sample (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, const NcmVector *href, NcmRNG *rng)
{
  /*NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);*/
  /*NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;*/
  gint i, ret;

  for (i = 0; i < d; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    ncm_vector_set (y, i, u_i * ncm_vector_fast_get (href, i));
  }

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (y));
  NCM_TEST_GSL_RESULT ("ncm_stats_dist_nd_vbk_gauss_sample_mean_scale", ret);

  ncm_vector_add (y, mu);
}

static gdouble
_ncm_stats_dist_nd_vbk_gauss_kernel_eval_m2lnp (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const NcmVector *href)
{
  /*NcmStatsDistNdVbkGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);*/
  /*NcmStatsDistNdVbkGaussPrivate * const self = dndg->priv;*/
  gdouble m2lnp;
  gint ret;

  ncm_vector_memcpy (v, x);
  ncm_vector_sub (v, y);

  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (v));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_vbk_gauss_kernel_eval_m2lnp", ret);

  ncm_vector_div (v, href);

  ret = gsl_blas_ddot (ncm_vector_gsl (v), ncm_vector_gsl (v), &m2lnp);
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_vbk_gauss_kernel_eval_m2lnp", ret);

  return m2lnp;
}

/**
 * ncm_stats_dist_nd_vbk_gauss_new:
 * @dim: sample space dimension
 * @cv_type: a #NcmStatsDistNdCV
 * 
 * Creates a new #NcmStatsDistNdVbkGauss object with sample dimension @dim.
 * 
 * Returns: a new #NcmStatsDistNdVbkGauss.
 */
NcmStatsDistNdVbkGauss *
ncm_stats_dist_nd_vbk_gauss_new (const guint dim, const NcmStatsDistNdVbkCV cv_type)
{
  NcmStatsDistNdVbkGauss *dndg = g_object_new (NCM_TYPE_STATS_DIST_ND_VBK_GAUSS,
                                               "dimension", dim,
                                               "CV-type",   cv_type,
                                               NULL);
  return dndg;
}

/**
 * ncm_stats_dist_nd_vbk_gauss_ref:
 * @dndg: a #NcmStatsDistNdVbkGauss
 *
 * Increase the reference of @stats_dist_nd_vbk_gauss by one.
 *
 * Returns: (transfer full): @stats_dist_nd_vbk_gauss.
 */
NcmStatsDistNdVbkGauss *
ncm_stats_dist_nd_vbk_gauss_ref (NcmStatsDistNdVbkGauss *dndg)
{
  return g_object_ref (dndg);
}

/**
 * ncm_stats_dist_nd_vbk_gauss_free:
 * @dndg: a #NcmStatsDistNdVbkGauss
 *
 * Decrease the reference count of @stats_dist_nd_vbk_gauss by one.
 *
 */
void
ncm_stats_dist_nd_vbk_gauss_free (NcmStatsDistNdVbkGauss *dndg)
{
  g_object_unref (dndg);
}

/**
 * ncm_stats_dist_nd_vbk_gauss_clear:
 * @dndg: a #NcmStatsDistNdVbkGauss
 *
 * Decrease the reference count of @stats_dist_nd_vbk_gauss by one, and sets the pointer *@stats_dist_nd_vbk_gauss to
 * NULL.
 *
 */
void
ncm_stats_dist_nd_vbk_gauss_clear (NcmStatsDistNdVbkGauss **dndg)
{
  g_clear_object (dndg);
}


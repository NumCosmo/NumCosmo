/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd_kde_gauss.c
 *
 *  Wed November 07 17:41:47 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_kde_gauss.c
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
 * SECTION:ncm_stats_dist_nd_kde_gauss
 * @title: NcmStatsDistNdKDEGauss
 * @short_description: An N dimensional probability distributions using gaussian KDE
 *
 * An arbitrary N dimensional probability distribution using gaussian KDE.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_nd_kde_gauss.h"
#include "math/ncm_stats_vec.h"
#include "math/ncm_c.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include "gslextras/cqp/gsl_cqp.h"
#include "libqp.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsDistNdKDEGaussPrivate
{
  gint place_holder;
};

enum
{
  PROP_0,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistNdKDEGauss, ncm_stats_dist_nd_kde_gauss, NCM_TYPE_STATS_DIST_ND);

static void
ncm_stats_dist_nd_kde_gauss_init (NcmStatsDistNdKDEGauss *dndg)
{
  /*NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv = ncm_stats_dist_nd_kde_gauss_get_instance_private (dndg);*/

  dndg->priv = ncm_stats_dist_nd_kde_gauss_get_instance_private (dndg);
}

static void
_ncm_stats_dist_nd_kde_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_GAUSS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_kde_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;*/

  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_GAUSS (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_kde_gauss_dispose (GObject *object)
{
  /*NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);*/
  /*NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;*/

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_kde_gauss_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_kde_gauss_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_kde_gauss_parent_class)->finalize (object);
}

static gdouble _ncm_stats_dist_nd_kde_gauss_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n);
static gdouble _ncm_stats_dist_nd_kde_gauss_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const gdouble href);
static void _ncm_stats_dist_nd_kde_gauss_prepare_IM (NcmStatsDistNd *dnd, GPtrArray *invUsample, const gint d, const gint n, const gdouble href, const gdouble href2, NcmMatrix *IM);

static gdouble _ncm_stats_dist_nd_kde_gauss_eval (NcmStatsDistNd *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const gdouble href, const gdouble href2);
static void _ncm_stats_dist_nd_kde_gauss_kernel_sample (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, gdouble href, NcmRNG *rng);
static gdouble _ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const gdouble href, const gdouble href2);

static void
ncm_stats_dist_nd_kde_gauss_class_init (NcmStatsDistNdKDEGaussClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmStatsDistNdClass *dnd_class  = NCM_STATS_DIST_ND_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_nd_kde_gauss_set_property;
  object_class->get_property = &_ncm_stats_dist_nd_kde_gauss_get_property;
  object_class->dispose      = &_ncm_stats_dist_nd_kde_gauss_dispose;
  object_class->finalize     = &_ncm_stats_dist_nd_kde_gauss_finalize;
  
  dnd_class->get_rot_bandwidth = &_ncm_stats_dist_nd_kde_gauss_get_rot_bandwidth;
  dnd_class->get_kernel_lnnorm = &_ncm_stats_dist_nd_kde_gauss_get_kernel_lnnorm;
  dnd_class->prepare_IM        = &_ncm_stats_dist_nd_kde_gauss_prepare_IM;
  dnd_class->eval              = &_ncm_stats_dist_nd_kde_gauss_eval;
  dnd_class->kernel_sample     = &_ncm_stats_dist_nd_kde_gauss_kernel_sample;
  dnd_class->kernel_eval_m2lnp = &_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp;
}

static gdouble
_ncm_stats_dist_nd_kde_gauss_f (NcmStatsDistNdKDEGaussPrivate * const self, const gdouble chi2)
{
  return exp (- 0.5 * chi2);
}

static gdouble
_ncm_stats_dist_nd_kde_gauss_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n)
{
  return pow (4.0 / (n * (d + 2.0)), 1.0 / (d + 4.0));
}

static gdouble
_ncm_stats_dist_nd_kde_gauss_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const gdouble href)
{
  return 0.5 * (d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (cov_decomp)) + d * log (href);
}

static void 
_ncm_stats_dist_nd_kde_gauss_prepare_IM (NcmStatsDistNd *dnd, GPtrArray *invUsample, const gint d, const gint n, const gdouble href, const gdouble href2, NcmMatrix *IM)
{
  NcmStatsDistNdKDEGaussPrivate * const self = NCM_STATS_DIST_ND_KDE_GAUSS (dnd)->priv;
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
        chi2_ij += gsl_pow_2 (ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (row_j, k));
      }

      p_ij = _ncm_stats_dist_nd_kde_gauss_f (self, chi2_ij / href2);
      
      ncm_matrix_set (IM, i, j, p_ij);
      ncm_matrix_set (IM, j, i, p_ij);
    }
  }
}

static gdouble 
_ncm_stats_dist_nd_kde_gauss_eval (NcmStatsDistNd *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const gdouble href, const gdouble href2)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i;
  
  for (i = 0; i < n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (invUsample, i);
    gdouble e_i, t, chi2_i = 0.0;
    gint k;

    for (k = 0; k < d; k++)
    {
      chi2_i += gsl_pow_2 (ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (invUy, k));
    }

    e_i  = ncm_vector_get (weights, i) * _ncm_stats_dist_nd_kde_gauss_f (self, chi2_i / href2);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
  }

  return s;
}

static void
_ncm_stats_dist_nd_kde_gauss_kernel_sample (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, gdouble href, NcmRNG *rng)
{
  /*NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);*/
  /*NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;*/
  gint i, ret;

  for (i = 0; i < d; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    ncm_vector_set (y, i, u_i);
  }

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (y));
  NCM_TEST_GSL_RESULT ("ncm_stats_dist_nd_kde_gauss_sample_mean_scale", ret);

  ncm_vector_scale (y, href);
  ncm_vector_add (y, mu);
}

static gdouble
_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const gdouble href, const gdouble href2)
{
  /*NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);*/
  /*NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;*/
  gdouble m2lnp;
  gint ret;

  ncm_vector_memcpy (v, x);
  ncm_vector_sub (v, y);

  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (cov_decomp), ncm_vector_gsl (v));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp", ret);

  ret = gsl_blas_ddot (ncm_vector_gsl (v), ncm_vector_gsl (v), &m2lnp);
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp", ret);

  return m2lnp / href2;
}

/**
 * ncm_stats_dist_nd_kde_gauss_new:
 * @dim: sample space dimension
 * @cv_type: a #NcmStatsDistNdCV
 * 
 * Creates a new #NcmStatsDistNdKDEGauss object with sample dimension @dim.
 * 
 * Returns: a new #NcmStatsDistNdKDEGauss.
 */
NcmStatsDistNdKDEGauss *
ncm_stats_dist_nd_kde_gauss_new (const guint dim, const NcmStatsDistNdCV cv_type)
{
  NcmStatsDistNdKDEGauss *dndg = g_object_new (NCM_TYPE_STATS_DIST_ND_KDE_GAUSS,
                                               "dimension", dim,
                                               "CV-type",   cv_type,
                                               NULL);
  return dndg;
}

/**
 * ncm_stats_dist_nd_kde_gauss_ref:
 * @dndg: a #NcmStatsDistNdKDEGauss
 *
 * Increase the reference of @stats_dist_nd_kde_gauss by one.
 *
 * Returns: (transfer full): @stats_dist_nd_kde_gauss.
 */
NcmStatsDistNdKDEGauss *
ncm_stats_dist_nd_kde_gauss_ref (NcmStatsDistNdKDEGauss *dndg)
{
  return g_object_ref (dndg);
}

/**
 * ncm_stats_dist_nd_kde_gauss_free:
 * @dndg: a #NcmStatsDistNdKDEGauss
 *
 * Decrease the reference count of @stats_dist_nd_kde_gauss by one.
 *
 */
void
ncm_stats_dist_nd_kde_gauss_free (NcmStatsDistNdKDEGauss *dndg)
{
  g_object_unref (dndg);
}

/**
 * ncm_stats_dist_nd_kde_gauss_clear:
 * @dndg: a #NcmStatsDistNdKDEGauss
 *
 * Decrease the reference count of @stats_dist_nd_kde_gauss by one, and sets the pointer *@stats_dist_nd_kde_gauss to
 * NULL.
 *
 */
void
ncm_stats_dist_nd_kde_gauss_clear (NcmStatsDistNdKDEGauss **dndg)
{
  g_clear_object (dndg);
}


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
#include "gslextras/cqp/gsl_cqp.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsDistNdKDEGaussPrivate
{
  NcmStatsVec *sample;
  NcmMatrix *cov_decomp;
  NcmVector *weights;
  NcmVector *v;
  gdouble href;
  gdouble href2;
  gdouble lnnorm;
  gdouble min_m2lnp;
  guint n;
  guint d;
  GArray *sampling;
};

enum
{
  PROP_0,
  PROP_N
};

G_DEFINE_TYPE_WITH_CODE (NcmStatsDistNdKDEGauss, ncm_stats_dist_nd_kde_gauss, NCM_TYPE_STATS_DIST_ND, G_ADD_PRIVATE (NcmStatsDistNdKDEGauss));

static void
ncm_stats_dist_nd_kde_gauss_init (NcmStatsDistNdKDEGauss *dndg)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv = G_TYPE_INSTANCE_GET_PRIVATE (dndg, NCM_TYPE_STATS_DIST_ND_KDE_GAUSS, NcmStatsDistNdKDEGaussPrivate);

  self->sample     = NULL;
  self->cov_decomp = NULL;
  self->weights    = NULL;
  self->v          = NULL;
  self->href       = 0.0;
  self->href2      = 0.0;
  self->lnnorm     = 0.0;
  self->min_m2lnp  = 0.0;
  self->n          = 0;
  self->d          = 0;
  self->sampling   = g_array_new (FALSE, FALSE, sizeof (guint));
}

static void
_ncm_stats_dist_nd_kde_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_GAUSS (object));

  switch (prop_id)
  {
    case PROP_N:
      g_assert_not_reached ();
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_kde_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;

  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_GAUSS (object));

  switch (prop_id)
  {
    case PROP_N:
      g_value_set_uint (value, ncm_stats_vec_nitens (self->sample));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_kde_gauss_dispose (GObject *object)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;

  ncm_stats_vec_clear (&self->sample);
  ncm_matrix_clear (&self->cov_decomp);
  ncm_vector_clear (&self->weights);
  ncm_vector_clear (&self->v);

  g_clear_pointer (&self->sampling, g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_kde_gauss_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_kde_gauss_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_kde_gauss_parent_class)->finalize (object);
}

static void _ncm_stats_dist_nd_kde_gauss_set_dim (NcmStatsDistNd *dnd, const guint dim);
static void _ncm_stats_dist_nd_kde_gauss_prepare (NcmStatsDistNd *dnd);
static void _ncm_stats_dist_nd_kde_gauss_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp, gboolean normalized);
static gdouble _ncm_stats_dist_nd_kde_gauss_eval (NcmStatsDistNd *dnd, NcmVector *y);
static gdouble _ncm_stats_dist_nd_kde_gauss_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *y);
static void _ncm_stats_dist_nd_kde_gauss_sample (NcmStatsDistNd *dnd, NcmVector *y, NcmRNG *rng);

static void
ncm_stats_dist_nd_kde_gauss_class_init (NcmStatsDistNdKDEGaussClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmStatsDistNdClass *dnd_class  = NCM_STATS_DIST_ND_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_nd_kde_gauss_set_property;
  object_class->get_property = &_ncm_stats_dist_nd_kde_gauss_get_property;
  object_class->dispose      = &_ncm_stats_dist_nd_kde_gauss_dispose;
  object_class->finalize     = &_ncm_stats_dist_nd_kde_gauss_finalize;

  g_object_class_install_property (object_class,
                                   PROP_N,
                                   g_param_spec_uint ("N",
                                                      NULL,
                                                      "sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  dnd_class->set_dim        = &_ncm_stats_dist_nd_kde_gauss_set_dim;
  dnd_class->prepare        = &_ncm_stats_dist_nd_kde_gauss_prepare;
  dnd_class->prepare_interp = &_ncm_stats_dist_nd_kde_gauss_prepare_interp;
  dnd_class->eval           = &_ncm_stats_dist_nd_kde_gauss_eval;
  dnd_class->eval_m2lnp     = &_ncm_stats_dist_nd_kde_gauss_eval_m2lnp;
  dnd_class->sample         = &_ncm_stats_dist_nd_kde_gauss_sample;

}

static void 
_ncm_stats_dist_nd_kde_gauss_set_dim (NcmStatsDistNd *dnd, const guint dim)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;

  ncm_stats_vec_clear (&self->sample);

  self->sample     = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, TRUE);
  self->cov_decomp = ncm_matrix_new (dim, dim);
  self->v          = ncm_vector_new (dim);
}

static void 
_ncm_stats_dist_nd_kde_gauss_prepare (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gint i;

  self->n = ncm_stats_vec_nitens (self->sample);
  self->d = ncm_stats_vec_len (self->sample);
  
  ncm_matrix_memcpy (self->cov_decomp, ncm_stats_vec_peek_cov_matrix (self->sample, 0));
  ncm_vector_clear (&self->weights);

  self->weights = ncm_vector_new (self->n);
  
  if (ncm_matrix_cholesky_decomp (self->cov_decomp, 'U') != 0)
    ncm_matrix_nearPD (self->cov_decomp, 'U', TRUE);

  self->href   = pow (4.0 / (self->n * (self->d + 2.0)), 1.0 / (self->d + 4.0));
  self->href2  = self->href * self->href;
  self->lnnorm = 0.5 * (self->d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (self->cov_decomp)) + self->d * log (self->href);

  for (i = 0; i < self->n; i++)
    ncm_vector_set (self->weights, i, 1.0 / (1.0 * self->n));
}

void LowRankQP (gint *n, gint *m, gint *p, gint *method, gint *verbose, gint *niter, gdouble *Q, gdouble *c, gdouble *A, gdouble *b, gdouble *u, gdouble *alpha, gdouble *beta, gdouble *xi, gdouble *zeta);

static void 
_ncm_stats_dist_nd_kde_gauss_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp, gboolean normalized)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;

  _ncm_stats_dist_nd_kde_gauss_prepare (dnd);
  g_assert_cmpuint (ncm_vector_len (m2lnp), ==, self->n);

  {
    gint n            = self->n;
    gint nc           = 1;
    NcmMatrix *A      = ncm_matrix_new (nc, self->n);
    NcmMatrix *IM     = ncm_matrix_new (self->n, self->n);
    NcmVector *b      = ncm_vector_new (nc);
    NcmVector *beta   = ncm_vector_new (nc);
    NcmVector *d      = ncm_vector_new (self->n);
    NcmVector *alpha  = ncm_vector_new (self->n);
    NcmVector *xi     = ncm_vector_new (self->n);
    NcmVector *zeta   = ncm_vector_new (self->n);
    gint method       = 2;
    gint verbose      = 0;
    gint maxiter      = 200;
    gint i, ret;

    self->min_m2lnp = GSL_POSINF;
    
    if (self->n > 10000)
      g_warning ("_ncm_stats_dist_nd_kde_gauss_prepare_interp: very large system n = %u!", self->n);

    for (i = 0; i < self->n; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
      NcmVector *row_i      = ncm_stats_vec_peek_row (self->sample, i);
      gint j;

      self->min_m2lnp = MIN (self->min_m2lnp, m2lnp_i);

      ncm_matrix_set (IM, i, i, exp (- self->lnnorm));

      for (j = i + 1; j < self->n; j++)
      {
        NcmVector *row_j = ncm_stats_vec_peek_row (self->sample, j);
        gdouble m2lnp_ij;

        ncm_vector_memcpy (self->v, row_i);
        ncm_vector_sub (self->v, row_j);

        ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                              ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
        NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_prepare_interp", ret);

        ret = gsl_blas_ddot (ncm_vector_gsl (self->v), ncm_vector_gsl (self->v), &m2lnp_ij);
        NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_prepare_interp", ret);

        ncm_matrix_set (IM, i, j, exp (- 0.5 * m2lnp_ij / self->href2 - self->lnnorm));
        ncm_matrix_set (IM, j, i, exp (- 0.5 * m2lnp_ij / self->href2 - self->lnnorm));
      }
    }

    if (normalized)
      self->min_m2lnp = 0.0;

    /*printf ("min_m2lnp: % 22.15g\n", self->min_m2lnp);*/

    for (i = 0; i < self->n; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
      ncm_vector_set (self->weights, i, exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
    }
    
    /* Setting RHS */
    ncm_vector_scale (self->weights, -1.0);

    /* Setting constraint value */
    ncm_vector_set (b, 0, 1.0);

    /* Setting constraint matrix */
    ncm_matrix_add_constant (A, 1.0);

    /* Setting upper-bound */
    ncm_vector_set_all (d, 1.0);

    /* Setting output vectors to zero */
    ncm_vector_set_zero (alpha);
    ncm_vector_set_zero (beta);
    ncm_vector_set_zero (xi);
    ncm_vector_set_zero (zeta);

    if (normalized)
      nc = 1;
    else
      nc = 0;

    LowRankQP (&n, &n, &nc, &method, &verbose, &maxiter, 
               ncm_matrix_data (IM), ncm_vector_data (self->weights),
               ncm_matrix_data (A), ncm_vector_data (b), ncm_vector_data (d), 
               ncm_vector_data (alpha), ncm_vector_data (beta), ncm_vector_data (xi), ncm_vector_data (zeta));

    ncm_vector_memcpy (self->weights, alpha);

    if (FALSE)
    {
      gdouble W = 0.0;
      for (i = 0; i < self->n; i++)
      {
        W += ncm_vector_get (self->weights, i);
        printf("%9.6e ", ncm_vector_get (self->weights, i));
      }
      printf (" = % 22.15g\n", W);
    }
    
    ncm_matrix_free (A);
    ncm_matrix_free (IM);
    ncm_vector_free (b);
    ncm_vector_free (beta);
    ncm_vector_free (d);
    ncm_vector_free (alpha);
    ncm_vector_free (xi);
    ncm_vector_free (zeta);
  }
}

static gdouble 
_ncm_stats_dist_nd_kde_gauss_eval (NcmStatsDistNd *dnd, NcmVector *y)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i, ret;

  for (i = 0; i < self->n; i++)
  {
    NcmVector *row_i = ncm_stats_vec_peek_row (self->sample, i);
    gdouble e_i, t, m2lnp_i = 0.0;

    ncm_vector_memcpy (self->v, row_i);
    ncm_vector_sub (self->v, y);

    ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                          ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_eval", ret);

    ret = gsl_blas_ddot (ncm_vector_gsl (self->v), ncm_vector_gsl (self->v), &m2lnp_i);
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_eval", ret);

    e_i  = ncm_vector_get (self->weights, i) * exp (- 0.5 * m2lnp_i / self->href2 - self->lnnorm);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
  }

  return s * exp (-0.5 * self->min_m2lnp);
}

static gdouble 
_ncm_stats_dist_nd_kde_gauss_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *y)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i, ret;

  for (i = 0; i < self->n; i++)
  {
    NcmVector *row_i = ncm_stats_vec_peek_row (self->sample, i);
    gdouble e_i, t, m2lnp_i = 0.0;

    ncm_vector_memcpy (self->v, row_i);
    ncm_vector_sub (self->v, y);

    ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                          ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_eval", ret);

    ret = gsl_blas_ddot (ncm_vector_gsl (self->v), ncm_vector_gsl (self->v), &m2lnp_i);
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_eval", ret);

    e_i  = ncm_vector_get (self->weights, i) * exp (- 0.5 * m2lnp_i / self->href2 - self->lnnorm);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
  }

  /*printf ("% 22.15g\n", self->min_m2lnp);*/
  
  return -2.0 * log (s) + self->min_m2lnp;
}

static void 
_ncm_stats_dist_nd_kde_gauss_sample (NcmStatsDistNd *dnd, NcmVector *y, NcmRNG *rng)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  const guint n = ncm_vector_len (self->weights);
  gint i;
  
  g_array_set_size (self->sampling, ncm_vector_len (self->weights));
  gsl_ran_multinomial (rng->r, n, 1, ncm_vector_data (self->weights), (guint *)self->sampling->data);

  for (i = 0; i < n; i++)
  {
    if (g_array_index (self->sampling, guint, i) > 0)
    {
      NcmVector *y_i = ncm_stats_vec_peek_row (self->sample, i);
      gint j, ret;

      for (j = 0; j < self->d; j++)
      {
        const gdouble u_j = gsl_ran_ugaussian (rng->r);
        ncm_vector_set (y, j, u_j);
      }

      /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
      ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                            ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (y));
      NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_sample", ret);

      ncm_vector_scale (y, self->href);
      ncm_vector_add (y, y_i);

      break;
    }
  }
}

/**
 * ncm_stats_dist_nd_kde_gauss_new:
 * @dim: sample space dimension
 * 
 * Creates a new #NcmStatsDistNdKDEGauss object with sample dimension @dim.
 * 
 * Returns: a new #NcmStatsDistNdKDEGauss.
 */
NcmStatsDistNdKDEGauss *
ncm_stats_dist_nd_kde_gauss_new (const guint dim)
{
  NcmStatsDistNdKDEGauss *dndg = g_object_new (NCM_TYPE_STATS_DIST_ND_KDE_GAUSS,
                                               "dimension", dim,
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

/**
 * ncm_stats_dist_nd_kde_gauss_add_obs_weight:
 * @dndg: a #NcmStatsDistNdKDEGauss
 * @y: a #NcmVector
 * @w: weight
 *
 * Adds a new point @y to the sample with weight @w.
 * 
 */
void 
ncm_stats_dist_nd_kde_gauss_add_obs_weight (NcmStatsDistNdKDEGauss *dndg, NcmVector *y, const gdouble w)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  ncm_stats_vec_append_weight (self->sample, y, w, TRUE);
}

/**
 * ncm_stats_dist_nd_kde_gauss_add_obs:
 * @dndg: a #NcmStatsDistNdKDEGauss
 * @y: a #NcmVector
 *
 * Adds a new point @y to the sample with weight 1.0.
 * 
 */
void 
ncm_stats_dist_nd_kde_gauss_add_obs (NcmStatsDistNdKDEGauss *dndg, NcmVector *y)
{
  ncm_stats_dist_nd_kde_gauss_add_obs_weight (dndg, y, 1.0);
}


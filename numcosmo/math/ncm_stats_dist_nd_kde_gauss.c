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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsDistNdKDEGaussPrivate
{
  NcmStatsVec *sample;
  NcmMatrix *cov_decomp;
  NcmMatrix *log_cov;
  NcmMatrix *sample_matrix;
  GPtrArray *smatrix_rows;
  NcmVector *weights;
  NcmVector *v;
  gdouble over_smooth;
  gboolean LOOCV;
  gdouble href;
  gdouble href2;
  gdouble lnnorm;
  gdouble us_lnnorm;
  gdouble min_m2lnp;
  gdouble max_m2lnp;
  guint n;
  guint alloc_n;
  guint d;
  GArray *sampling;
  guint nearPD_maxiter;
  NcmMatrix *A;
  NcmMatrix *IM;
  NcmVector *b;
  NcmVector *beta;
  NcmVector *dv;
  NcmVector *alpha;
  NcmVector *xi;
  NcmVector *zeta;
  NcmLapackWS *lapack_ws;
  GArray *ipiv;
};

enum
{
  PROP_0,
  PROP_N,
  PROP_OVER_SMOOTH,
  PROP_NEARPD_MAXITER,
  PROP_LOOCV,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistNdKDEGauss, ncm_stats_dist_nd_kde_gauss, NCM_TYPE_STATS_DIST_ND);

static void
ncm_stats_dist_nd_kde_gauss_init (NcmStatsDistNdKDEGauss *dndg)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv = ncm_stats_dist_nd_kde_gauss_get_instance_private (dndg);

  self->sample           = NULL;
  self->cov_decomp       = NULL;
  self->log_cov          = NULL;
  self->sample_matrix    = NULL;
  self->smatrix_rows     = g_ptr_array_new ();
  self->weights          = NULL;
  self->v                = NULL;
  self->over_smooth      = 0.0;
  self->LOOCV            = FALSE;
  self->href             = 0.0;
  self->href2            = 0.0;
  self->lnnorm           = 0.0;
  self->us_lnnorm        = 0.0;
  self->min_m2lnp        = 0.0;
  self->max_m2lnp        = 0.0;
  self->n                = 0;
  self->alloc_n          = 0;
  self->d                = 0;
  self->sampling         = g_array_new (FALSE, FALSE, sizeof (guint));
  self->nearPD_maxiter   = 0;
  self->A                = NULL;
  self->IM               = NULL;
  self->b                = NULL;
  self->beta             = NULL;
  self->dv               = NULL;
  self->alpha            = NULL;
  self->xi               = NULL;
  self->zeta             = NULL;
  self->lapack_ws        = ncm_lapack_ws_new ();
  self->ipiv             = g_array_new (FALSE, FALSE, sizeof (guint));

  g_ptr_array_set_free_func (self->smatrix_rows, (GDestroyNotify) ncm_vector_free);
}

static void
_ncm_stats_dist_nd_kde_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (object);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND_KDE_GAUSS (object));

  switch (prop_id)
  {
    case PROP_N:
      g_assert_not_reached ();
      break;
    case PROP_OVER_SMOOTH:
      ncm_stats_dist_nd_kde_gauss_set_over_smooth (dndg, g_value_get_double (value));
      break;
    case PROP_NEARPD_MAXITER:
      self->nearPD_maxiter = g_value_get_uint (value);
      break;
    case PROP_LOOCV:
      ncm_stats_dist_nd_kde_gauss_set_LOOCV_bandwidth_adj (dndg, g_value_get_boolean (value));
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
    case PROP_OVER_SMOOTH:
      g_value_set_double (value, ncm_stats_dist_nd_kde_gauss_get_over_smooth (dndg));
      break;
    case PROP_NEARPD_MAXITER:
      g_value_set_uint (value, self->nearPD_maxiter);
      break;
    case PROP_LOOCV:
      g_value_set_boolean (value, ncm_stats_dist_nd_kde_gauss_get_LOOCV_bandwidth_adj (dndg));
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
  ncm_matrix_clear (&self->log_cov);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_vector_clear (&self->weights);
  ncm_vector_clear (&self->v);

  g_clear_pointer (&self->sampling, g_array_unref);
  g_clear_pointer (&self->smatrix_rows, g_ptr_array_unref);

  ncm_matrix_clear (&self->A);
  ncm_matrix_clear (&self->IM);
  ncm_vector_clear (&self->b);
  ncm_vector_clear (&self->beta);
  ncm_vector_clear (&self->dv);
  ncm_vector_clear (&self->alpha);
  ncm_vector_clear (&self->xi);
  ncm_vector_clear (&self->zeta);

  ncm_lapack_ws_clear (&self->lapack_ws);
  g_clear_pointer (&self->ipiv, g_array_unref);

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
static void _ncm_stats_dist_nd_kde_gauss_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp);
static gdouble _ncm_stats_dist_nd_kde_gauss_eval (NcmStatsDistNd *dnd, NcmVector *y);
static gdouble _ncm_stats_dist_nd_kde_gauss_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *y);
static void _ncm_stats_dist_nd_kde_gauss_sample (NcmStatsDistNd *dnd, NcmVector *y, NcmRNG *rng);
void _ncm_stats_dist_nd_kde_gauss_kernel_sample (NcmStatsDistNd *dnd, NcmVector *y, NcmVector *mu, gdouble scale, NcmRNG *rng);
gdouble _ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *y, gdouble scale);
static void _ncm_stats_dist_nd_kde_gauss_reset (NcmStatsDistNd *dnd);

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

  g_object_class_install_property (object_class,
                                   PROP_OVER_SMOOTH,
                                   g_param_spec_double ("over-smooth",
                                                        NULL,
                                                        "Oversmooth distribution",
                                                        1.0e-5, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_LOOCV,
                                   g_param_spec_boolean ("LOOCV",
                                                         NULL,
                                                         "Leave one out cross validation",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_NEARPD_MAXITER,
                                   g_param_spec_uint ("nearPD-maxiter",
                                                      NULL,
                                                      "Maximum number of iterations in the nearPD call",
                                                      1, G_MAXUINT, 200,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  dnd_class->set_dim           = &_ncm_stats_dist_nd_kde_gauss_set_dim;
  dnd_class->prepare           = &_ncm_stats_dist_nd_kde_gauss_prepare;
  dnd_class->prepare_interp    = &_ncm_stats_dist_nd_kde_gauss_prepare_interp;
  dnd_class->eval              = &_ncm_stats_dist_nd_kde_gauss_eval;
  dnd_class->eval_m2lnp        = &_ncm_stats_dist_nd_kde_gauss_eval_m2lnp;
  dnd_class->sample            = &_ncm_stats_dist_nd_kde_gauss_sample;
  dnd_class->kernel_sample     = &_ncm_stats_dist_nd_kde_gauss_kernel_sample;
  dnd_class->kernel_eval_m2lnp = &_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp;
  dnd_class->reset             = &_ncm_stats_dist_nd_kde_gauss_reset;
}

static void 
_ncm_stats_dist_nd_kde_gauss_set_dim (NcmStatsDistNd *dnd, const guint dim)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;

  ncm_stats_vec_clear (&self->sample);

  ncm_matrix_clear (&self->cov_decomp);
  ncm_matrix_clear (&self->log_cov);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_vector_clear (&self->v);

  self->sample     = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, TRUE);
  self->cov_decomp = ncm_matrix_new (dim, dim);
  self->log_cov    = ncm_matrix_new (dim, dim);
  self->v          = ncm_vector_new (dim);
}

static void 
_ncm_stats_dist_nd_kde_gauss_prepare_cov (NcmStatsDistNdKDEGaussPrivate * const self)
{
  gint i, ret;

  if (self->sample_matrix == NULL)
  {
    self->sample_matrix = ncm_matrix_new (self->n, self->d);
    g_ptr_array_set_size (self->smatrix_rows, 0);    
    for (i = 0; i < self->n; i++)
    {
      NcmVector *row_i = ncm_matrix_get_row (self->sample_matrix, i);
      g_ptr_array_add (self->smatrix_rows, row_i);
    }
  }
  else if ((self->n != ncm_matrix_nrows (self->sample_matrix)) || (self->d != ncm_matrix_ncols (self->sample_matrix)))
  {
    ncm_matrix_clear (&self->sample_matrix);
    self->sample_matrix = ncm_matrix_new (self->n, self->d);

    g_ptr_array_set_size (self->smatrix_rows, 0);    
    for (i = 0; i < self->n; i++)
    {
      NcmVector *row_i = ncm_matrix_get_row (self->sample_matrix, i);
      g_ptr_array_add (self->smatrix_rows, row_i);
    }
  }

  for (i = 0; i < self->n; i++)
    ncm_matrix_set_row (self->sample_matrix, i, ncm_stats_vec_peek_row (self->sample, i));

  ret = gsl_blas_dtrsm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, ncm_matrix_gsl (self->cov_decomp), ncm_matrix_gsl (self->sample_matrix));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_prepare", ret);  
}

static void 
_ncm_stats_dist_nd_kde_gauss_prepare (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gint i;

  self->n = ncm_stats_vec_nitens (self->sample);
  self->d = ncm_stats_vec_len (self->sample);

  g_assert_cmpuint (self->n, >, 1);
  
  ncm_matrix_memcpy (self->cov_decomp, ncm_stats_vec_peek_cov_matrix (self->sample, 0));
  if (ncm_matrix_cholesky_decomp (self->cov_decomp, 'U') != 0)
  {
    if (ncm_matrix_nearPD (self->cov_decomp, 'U', TRUE, self->nearPD_maxiter) != 0)
    {
      ncm_matrix_set_zero (self->cov_decomp);
      for (i = 0; i < self->d; i++)
      {
        ncm_matrix_set (self->cov_decomp, i, i, ncm_stats_vec_get_cov (self->sample, i, i));
      }
      /*ncm_matrix_log_vals (self->cov_decomp, "# COV: ", "% 22.15g");*/
      g_assert_cmpint (ncm_matrix_cholesky_decomp (self->cov_decomp, 'U'), ==, 0);
    }
  }

  self->href      = self->over_smooth * pow (4.0 / (self->n * (self->d + 2.0)), 1.0 / (self->d + 4.0));
  self->href2     = self->href * self->href;
  self->us_lnnorm = 0.5 * (self->d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (self->cov_decomp));
  self->lnnorm    = self->us_lnnorm + self->d * log (self->href);

  self->min_m2lnp = 0.0;
  self->max_m2lnp = 0.0;

  if (self->weights == NULL)
  {
    self->weights = ncm_vector_new (self->n);
  }
  else if (self->n != ncm_vector_len (self->weights))
  {
    ncm_vector_clear (&self->weights);
    self->weights = ncm_vector_new (self->n);
  }
  
  for (i = 0; i < self->n; i++)
  {
    const gdouble wn = 1.0 / (1.0 * self->n);
    ncm_vector_set (self->weights, i, wn);
  }

  _ncm_stats_dist_nd_kde_gauss_prepare_cov (self);
}

void LowRankQP (gint *n, gint *m, gint *p, gint *method, gint *verbose, gint *niter, gdouble *Q, gdouble *c, gdouble *A, gdouble *b, gdouble *u, gdouble *alpha, gdouble *beta, gdouble *xi, gdouble *zeta);

static void 
_ncm_stats_dist_nd_kde_gauss_prepare_IM (NcmStatsDistNdKDEGaussPrivate * const self)
{
  gint i;

  for (i = 0; i < self->n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (self->smatrix_rows, i);
    gint j;

    ncm_matrix_set (self->IM, i, i, 1.0);

    for (j = i + 1; j < self->n; j++)
    {
      NcmVector *row_j = g_ptr_array_index (self->smatrix_rows, j);
      gdouble m2lnp_ij = 0.0;
      gdouble p_ij;
      gint k;

      for (k = 0; k < self->d; k++)
      {
        m2lnp_ij += gsl_pow_2 (ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (row_j, k));
      }

      p_ij = exp (- 0.5 * m2lnp_ij / self->href2);
      
      ncm_matrix_set (self->IM, i, j, p_ij);
      ncm_matrix_set (self->IM, j, i, p_ij);
    }
  }
}

static gdouble 
_ncm_stats_dist_nd_kde_gauss_LOOCV_err2 (gdouble h, gpointer user_data)
{
  NcmStatsDistNdKDEGaussPrivate * const self = user_data;
  gdouble err       = 0.0;
	guint i;
	gint ret;

  self->href      = h;
  self->href2     = self->href * self->href;
  self->us_lnnorm = 0.5 * (self->d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (self->cov_decomp));
  self->lnnorm    = self->us_lnnorm + self->d * log (self->href);
  
  _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);
	ret = ncm_matrix_cholesky_decomp (self->IM, 'U');	

	if (ret != 0)
	{
    _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);

    g_array_set_size (self->ipiv, self->n);
    
    ret = ncm_lapack_dsytrf ('U', self->n, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), self->lapack_ws);
    if (ret != 0)
    {
      printf ("ERR DEC2 h % 22.15g err % 22.15e %d\n", h, 1.0e20, ret);
      return 1.0e20;
    }

    ncm_vector_memcpy (self->zeta, self->weights);
    ret = ncm_lapack_dsytrs ('U', self->n, 1, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), ncm_vector_data (self->zeta), ncm_vector_len (self->zeta));
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_lapack_dsytrs]: %d.", ret);
    
    ret = ncm_lapack_dsytri ('U', self->n, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), self->lapack_ws);
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_lapack_dsytrs]: %d.", ret);
	}
  else
  {
    ncm_vector_memcpy (self->zeta, self->weights);
    ret = ncm_matrix_cholesky_solve2 (self->IM, self->zeta, 'U');
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_matrix_cholesky_solve2]: %d.", ret);

    ret = ncm_matrix_cholesky_inverse (self->IM, 'U');
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_matrix_cholesky_inverse]: %d.", ret);
  }
  
	for (i = 0; i < self->n; i++)
	{
		err += gsl_pow_2 (ncm_vector_get (self->zeta, i) / (ncm_matrix_get (self->IM, i, i)));
	}

  return log (err);
} 

gdouble 
_ncm_stats_dist_nd_kde_gauss_LOOCV_err2_f (guint n, const gdouble *x, gdouble *grad, gpointer f_data)
{
  const gdouble fx = _ncm_stats_dist_nd_kde_gauss_LOOCV_err2 (x[0], f_data);
  /*printf ("# % 22.15g % 22.15g\n", x[0], fx);*/
  return fx;
}

gdouble 
_ncm_stats_dist_nd_kde_gauss_LOOCV_err2_full_f (guint n, const gdouble *x, gdouble *grad, gpointer f_data)
{
  NcmStatsDistNdKDEGaussPrivate * const self = f_data;
  gdouble err       = 0.0;
	guint i, j, k;
	gint ret;
  gdouble tr_log_cov = 0.0;

  k = 0;
  for (i = 0; i < self->d; i++)
  {
    for (j = i; j < self->d; j++)
    {
      ncm_matrix_set (self->log_cov, i, j, x[k]);
      if (i == j)
        tr_log_cov += x[k];
      k++;
    }
  }
  ncm_matrix_sym_exp_cholesky (self->log_cov, 'U', self->cov_decomp);
  _ncm_stats_dist_nd_kde_gauss_prepare_cov (self);
  
  self->href      = 1.0;
  self->href2     = self->href * self->href;
  self->us_lnnorm = 0.5 * (self->d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (self->cov_decomp));
  self->lnnorm    = self->us_lnnorm;
  
  _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);
	ret = ncm_matrix_cholesky_decomp (self->IM, 'U');	

	if (ret != 0)
	{
    _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);

    g_array_set_size (self->ipiv, self->n);
    
    ret = ncm_lapack_dsytrf ('U', self->n, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), self->lapack_ws);
    if (ret != 0)
    {
      /*printf ("ERR DEC2 h % 22.15g err % 22.15e %d\n", 0.0, 1.0e20, ret);*/
      return 1.0e20;
    }

    ncm_vector_memcpy (self->zeta, self->weights);
    ret = ncm_lapack_dsytrs ('U', self->n, 1, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), ncm_vector_data (self->zeta), ncm_vector_len (self->zeta));
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_lapack_dsytrs]: %d.", ret);
    
    ret = ncm_lapack_dsytri ('U', self->n, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), self->lapack_ws);
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_lapack_dsytrs]: %d.", ret);
	}
  else
  {
    ncm_vector_memcpy (self->zeta, self->weights);
    ret = ncm_matrix_cholesky_solve2 (self->IM, self->zeta, 'U');
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_matrix_cholesky_solve2]: %d.", ret);

    ret = ncm_matrix_cholesky_inverse (self->IM, 'U');
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_matrix_cholesky_inverse]: %d.", ret);
  }
  
	for (i = 0; i < self->n; i++)
	{
		err += gsl_pow_2 (ncm_vector_get (self->zeta, i) / (ncm_matrix_get (self->IM, i, i)));
	}

  if (FALSE)
  {
    NcmMatrix *cov = ncm_matrix_dup (self->cov_decomp);
    ncm_matrix_triang_to_sym (self->cov_decomp, 'U', TRUE, cov);

    k = 0;
    for (i = 0; i < self->d; i++)
    {
      gint j;
      for (j = 0; j < self->d; j++)
      {
        if (j >= i)
        {
          printf ("% 22.15g ", x[k]);
          k++;
        }
        else
          printf ("% 22.15g ", 0.0);
      }
      printf ("\n");
    }
    printf ("#--------------------------------------------------------------------------------------------------------------\n");
    
    for (i = 0; i < self->d; i++)
    {
      gint j;
      for (j = 0; j < self->d; j++)
      {
        printf ("% 22.15g ", ncm_matrix_get (cov, i, j));
      }
      printf ("\n");
    }
    printf ("#--------------------------------------------------------------------------------------------------------------\n");
    printf ("err: % 22.15g\n", err);
    
    ncm_matrix_free (cov);
  }
  
  return log (err);
}

void 
_ncm_stats_dist_nd_kde_gauss_LOOCV_err2_full_levmar_f (gdouble *x, gdouble *hx, gint m, gint n, gpointer f_data)
{
  NcmStatsDistNdKDEGaussPrivate * const self = f_data;
  gdouble err        = 0.0;
  gdouble tr_log_cov = 0.0;
	guint i, j, k;
	gint ret;

  k = 0;
  for (i = 0; i < self->d; i++)
  {
    tr_log_cov += x[k];
    for (j = i; j < self->d; j++)
    {
      ncm_matrix_set (self->log_cov, i, j, x[k]);
      k++;
    }
  }
  ncm_matrix_sym_exp_cholesky (self->log_cov, 'U', self->cov_decomp);
  _ncm_stats_dist_nd_kde_gauss_prepare_cov (self);
  
  self->href      = 1.0;
  self->href2     = self->href * self->href;
  self->us_lnnorm = 0.5 * (self->d * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (self->cov_decomp));
  self->lnnorm    = self->us_lnnorm;

  /*printf ("CMP DETS: % 22.15g % 22.15g\n", tr_log_cov, ncm_matrix_cholesky_lndet (self->cov_decomp));*/
  
  _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);
	ret = ncm_matrix_cholesky_decomp (self->IM, 'U');	

	if (ret != 0)
	{
    /*printf ("cholesky failed!\n");*/
    _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);

    g_array_set_size (self->ipiv, self->n);
    
    ret = ncm_lapack_dsytrf ('U', self->n, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), self->lapack_ws);
    if (ret != 0)
    {
      /*printf ("ERR DEC2 h % 22.15g err % 22.15e %d\n", 0.0, 1.0e20, ret);*/
    }

    ncm_vector_memcpy (self->zeta, self->weights);
    ret = ncm_lapack_dsytrs ('U', self->n, 1, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), ncm_vector_data (self->zeta), ncm_vector_len (self->zeta));
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_lapack_dsytrs]: %d.", ret);
    
    ret = ncm_lapack_dsytri ('U', self->n, ncm_matrix_data (self->IM), ncm_matrix_tda (self->IM), &g_array_index (self->ipiv, gint, 0), self->lapack_ws);
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_lapack_dsytrs]: %d.", ret);
	}
  else
  {
    ncm_vector_memcpy (self->zeta, self->weights);
    ret = ncm_matrix_cholesky_solve2 (self->IM, self->zeta, 'U');
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_matrix_cholesky_solve2]: %d.", ret);

    ret = ncm_matrix_cholesky_inverse (self->IM, 'U');
    if (ret != 0)
      g_error ("_ncm_stats_dist_nd_kde_gauss_LOOCV_err2[ncm_matrix_cholesky_inverse]: %d.", ret);
  }
  
	for (i = 0; i < self->n; i++)
	{
		hx[i] = ncm_vector_get (self->zeta, i) / ncm_matrix_get (self->IM, i, i);
    err += gsl_pow_2 (hx[i]);
	}

  if (FALSE)
  {
    NcmMatrix *cov = ncm_matrix_dup (self->cov_decomp);
    ncm_matrix_triang_to_sym (self->cov_decomp, 'U', TRUE, cov);

    k = 0;
    for (i = 0; i < self->d; i++)
    {
      gint j;
      for (j = 0; j < self->d; j++)
      {
        if (j >= i)
        {
          printf ("% 22.15g ", x[k]);
          k++;
        }
        else
          printf ("% 22.15g ", 0.0);
      }
      printf ("\n");
    }
    printf ("#--------------------------------------------------------------------------------------------------------------\n");
    
    for (i = 0; i < self->d; i++)
    {
      gint j;
      for (j = 0; j < self->d; j++)
      {
        printf ("% 22.15g ", ncm_matrix_get (cov, i, j));
      }
      printf ("\n");
    }
    printf ("#--------------------------------------------------------------------------------------------------------------\n");
    printf ("err: % 22.15g lnnorm % 22.15g\n", err, self->lnnorm);
    ncm_matrix_free (cov);
  }
}

static void
_ncm_stats_dist_nd_kde_gauss_calib_href (NcmStatsDistNdKDEGaussPrivate * const self, NcmVector *m2lnp)
{
  gint i;

  for (i = 0; i < self->n; i++)
  {
    const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
    ncm_vector_set (self->weights, i, exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
  }

  if (TRUE)
  {
    const gint n     = (self->d * (self->d + 1)) / 2;
    gdouble *dht     = g_new0 (gdouble, n);
    gdouble *ht      = g_new0 (gdouble, n);
    NcmMatrix *cov   = ncm_matrix_dup (ncm_stats_vec_peek_cov_matrix (self->sample, 0));
    NcmMatrix *ln_cm = ncm_matrix_dup (ncm_stats_vec_peek_cov_matrix (self->sample, 0));
    gint k = 0;
    
    ncm_matrix_scale (cov, gsl_pow_2 (self->over_smooth * pow (4.0 / (self->n * (self->d + 2.0)), 1.0 / (self->d + 4.0))));
    /*printf ("\n");*/
    /*ncm_matrix_log_vals (cov, "ORIG_COV: ", "% 22.15g");*/

    ncm_matrix_sym_posdef_log (cov, 'U', ln_cm);

    /*ncm_matrix_log_vals (ln_cm, "OLOG_COV: ", "% 22.15g");*/

    for (i = 0; i < self->d; i++)
    {
      gint j;
      for (j = i; j < self->d; j++)
      {
        ht[k]  = ncm_matrix_get (ln_cm, i, j);
        dht[k] = 1.0;
        k++;
      }
    }

    {
      gdouble info[LM_INFO_SZ];
      gdouble opts[LM_OPTS_SZ];
      gint ret;

      opts[0] = LM_INIT_MU; 
      opts[1] = 1.0e-15; 
      opts[2] = 1.0e-15;
      opts[3] = 1.0e-20;
      opts[4] = -LM_DIFF_DELTA;

      ret = dlevmar_dif (&_ncm_stats_dist_nd_kde_gauss_LOOCV_err2_full_levmar_f,
                         ht, NULL, n, self->n, 10000,   
                         opts, info, NULL, NULL, self);
      if (ret < 0)
        g_error ("_ncm_stats_dist_nd_kde_gauss_calib_href: could not find the bestfit covariance matrix.");

      g_free (dht);
      g_free (ht);
    }
  }
  else
  {
    gdouble ht      = self->href;
    gdouble hi      = ht / 1000.0;
    gdouble hf      = ht * 100.0;
    gdouble last_hi = hi;
    gdouble last_hf = hf;
    guint iter      = 0;
    gint status     = 0;
    gsl_function F;
    
    gsl_min_fminimizer *fmin = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);

    F.params   = self;
    F.function = &_ncm_stats_dist_nd_kde_gauss_LOOCV_err2;

    gsl_min_fminimizer_set (fmin, &F, ht, hi, hf);
    /*printf ("[%d] % 22.15e [% 22.15e % 22.15e]  \n", status, ht, hi, hf);*/
    do {
      iter++;
      status = gsl_min_fminimizer_iterate (fmin);
      if (status)
        g_error ("_ncm_stats_dist_nd_kde_gauss_prepare: Cannot find minimum (%s)", gsl_strerror (status));

      ht = gsl_min_fminimizer_x_minimum (fmin);
      hi = gsl_min_fminimizer_x_lower (fmin);
      hf = gsl_min_fminimizer_x_upper (fmin);

      status = gsl_min_test_interval (hi, hf, 0.0, 1.0e-2);

      if ((status == GSL_CONTINUE) && (hi == last_hi) && (hf == last_hf))
      {
        g_warning ("_ncm_stats_dist_nd_kde_gauss_prepare: minimization not improving, giving up...");
        break;
      }

      /*printf ("[%d] % 22.15e [% 22.15e % 22.15e]  \n", status, ht, hi, hf);*/

      last_hi = hi;
      last_hf = hf;

    } while (status == GSL_CONTINUE && iter < 10000);

    self->href = gsl_min_fminimizer_x_minimum (fmin);
    gsl_min_fminimizer_free (fmin);
  }
}

static void 
_ncm_stats_dist_nd_kde_gauss_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  const gdouble dbl_limit = 1.0;
  gint i;

  _ncm_stats_dist_nd_kde_gauss_prepare (dnd);
  g_assert_cmpuint (ncm_vector_len (m2lnp), ==, self->n);

  if (self->n != self->alloc_n)
  {
    ncm_matrix_clear (&self->A);
    ncm_matrix_clear (&self->IM);
    ncm_vector_clear (&self->b);
    ncm_vector_clear (&self->beta);
    ncm_vector_clear (&self->dv);
    ncm_vector_clear (&self->alpha);
    ncm_vector_clear (&self->xi);
    ncm_vector_clear (&self->zeta);

    self->A       = ncm_matrix_new (1, self->n);
    self->IM      = ncm_matrix_new (self->n, self->n);
    self->b       = ncm_vector_new (1);
    self->beta    = ncm_vector_new (1);
    self->dv      = ncm_vector_new (self->n);
    self->alpha   = ncm_vector_new (self->n);
    self->xi      = ncm_vector_new (self->n);
    self->zeta    = ncm_vector_new (self->n);

    self->alloc_n = self->n;
  }
  
  self->min_m2lnp = GSL_POSINF;
  self->max_m2lnp = GSL_NEGINF;
    
  for (i = 0; i < self->n; i++)
  {
    const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
    self->min_m2lnp = MIN (self->min_m2lnp, m2lnp_i);
    self->max_m2lnp = MAX (self->max_m2lnp, m2lnp_i);
  }

  if ((-0.5 * (self->max_m2lnp - self->min_m2lnp) < dbl_limit * GSL_LOG_DBL_EPSILON) && TRUE)
  {
    const guint fi = ncm_vector_get_min_index (m2lnp);

    g_assert_cmpuint (fi, <, ncm_vector_len (m2lnp));
    
    /*printf ("# Single point fi %u\n", fi);*/

    ncm_vector_set_zero (self->weights);
    ncm_vector_set (self->weights, fi, 1.0);

    self->href   = 1.0;
    self->href2  = self->href * self->href;
    self->lnnorm = self->us_lnnorm + self->d * log (self->href);

    return;
  }

	if (self->LOOCV)
	{
    _ncm_stats_dist_nd_kde_gauss_calib_href (self, m2lnp);
	}
  
  /*printf ("# Using INTERP? % 22.15g % 22.15g | % 22.15g % 22.15g\n", self->max_m2lnp, self->min_m2lnp, -0.5 * (self->max_m2lnp - self->min_m2lnp), GSL_LOG_DBL_EPSILON);*/
  /*if ((-0.5 * (self->max_m2lnp - self->min_m2lnp) >= 3.0 * GSL_LOG_DBL_EPSILON) && TRUE)*/
  {
    gint n       = self->n;
    gint nc      = 1;
    gint method  = 2;
    gint verbose = 0;
    gint maxiter = 400;

    if (self->n > 10000)
      g_warning ("_ncm_stats_dist_nd_kde_gauss_prepare_interp: very large system n = %u!", self->n);

    /*printf ("# Using INTERP!\n");*/
    _ncm_stats_dist_nd_kde_gauss_prepare_IM (self);

    /*printf ("min_m2lnp: % 22.15g\n", self->min_m2lnp);*/
    /*ncm_matrix_log_vals (IM, "# IM: ", "% 12.5g");*/

    for (i = 0; i < self->n; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
      ncm_vector_set (self->weights, i, exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
      /*printf ("RHS:    %4d % 22.15e % 22.15g\n", i, exp (-0.5 * (m2lnp_i - self->min_m2lnp)), m2lnp_i);*/
    }
    
    /* Setting RHS */
    ncm_vector_scale (self->weights, -1.0);

    /* Setting constraint value */
    ncm_vector_set (self->b, 0, 1.0);

    /* Setting constraint matrix */
    ncm_matrix_set_all (self->A, 1.0);

    /* Setting upper-bound */
    ncm_vector_set_all (self->dv, 1.0);

    /* Setting output vectors to zero */
    ncm_vector_set_zero (self->alpha);
    ncm_vector_set_zero (self->beta);
    ncm_vector_set_zero (self->xi);
    ncm_vector_set_zero (self->zeta);

    nc = 0; /* Normalization is not useful... */
    LowRankQP (&n, &n, &nc, &method, &verbose, &maxiter, 
               ncm_matrix_data (self->IM), ncm_vector_data (self->weights),
               ncm_matrix_data (self->A), ncm_vector_data (self->b), ncm_vector_data (self->dv), 
               ncm_vector_data (self->alpha), ncm_vector_data (self->beta), ncm_vector_data (self->xi), ncm_vector_data (self->zeta));

    ncm_vector_memcpy (self->weights, self->alpha);    
  }

  /*ncm_vector_set_all (self->weights, 1.0);*/
  ncm_vector_scale (self->weights, 1.0 / ncm_vector_sum_cpts (self->weights));
}

static gdouble 
_ncm_stats_dist_nd_kde_gauss_eval (NcmStatsDistNd *dnd, NcmVector *y)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i, ret;

  ncm_vector_memcpy (self->v, y);
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_eval", ret);
  
  for (i = 0; i < self->n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (self->smatrix_rows, i);
    gdouble e_i, t, m2lnp_i = 0.0;
    gint k;

    for (k = 0; k < self->d; k++)
    {
      m2lnp_i += gsl_pow_2 (ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (self->v, k));
    }

    e_i  = ncm_vector_get (self->weights, i) * exp (- 0.5 * m2lnp_i / self->href2);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
  }

  return s * exp (- self->lnnorm);
}

static gdouble 
_ncm_stats_dist_nd_kde_gauss_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *y)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gdouble s = 0.0;
  gdouble c = 0.0;
  gint i, ret;

  ncm_vector_memcpy (self->v, y);
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_eval_m2lnp", ret);

  for (i = 0; i < self->n; i++)
  {
    NcmVector *row_i = g_ptr_array_index (self->smatrix_rows, i);
    gdouble e_i, t, m2lnp_i = 0.0;
    gint k;

    for (k = 0; k < self->d; k++)
    {
      m2lnp_i += gsl_pow_2 (ncm_vector_fast_get (row_i, k) - ncm_vector_fast_get (self->v, k));
    }

    e_i  = ncm_vector_get (self->weights, i) * exp (- 0.5 * m2lnp_i / self->href2);
    t    = s + e_i;
    c   += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s    = t;
    /*printf ("% 22.15g % 22.15g % 22.15g\n", e_i, s, self->lnnorm);*/
  }

  return -2.0 * (log (s) - self->lnnorm);
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
      _ncm_stats_dist_nd_kde_gauss_kernel_sample (dnd, y, y_i, self->href, rng);
      break;
    }
  }
}

void 
_ncm_stats_dist_nd_kde_gauss_kernel_sample (NcmStatsDistNd *dnd, NcmVector *y, NcmVector *mu, gdouble scale, NcmRNG *rng)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gint i, ret;

  for (i = 0; i < self->d; i++)
  {
    const gdouble u_i = gsl_ran_ugaussian (rng->r);
    ncm_vector_set (y, i, u_i);
  }

  /* CblasLower, CblasNoTrans => CblasUpper, CblasTrans */
  ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (y));
  NCM_TEST_GSL_RESULT ("ncm_stats_dist_nd_kde_gauss_sample_mean_scale", ret);

  ncm_vector_scale (y, scale);
  ncm_vector_add (y, mu);
}

gdouble 
_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *y, const gdouble scale)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  gdouble m2lnp;
  gint ret;

  ncm_vector_memcpy (self->v, x);
  ncm_vector_sub (self->v, y);

  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp", ret);

  ret = gsl_blas_ddot (ncm_vector_gsl (self->v), ncm_vector_gsl (self->v), &m2lnp);
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_kde_gauss_kernel_eval_m2lnp", ret);

  return m2lnp / gsl_pow_2 (scale) + 2.0 * (self->us_lnnorm + self->d * log (scale));
}

static void 
_ncm_stats_dist_nd_kde_gauss_reset (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdKDEGauss *dndg = NCM_STATS_DIST_ND_KDE_GAUSS (dnd);
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;

  ncm_stats_vec_reset (self->sample, TRUE);
}

/**
 * ncm_stats_dist_nd_kde_gauss_new:
 * @dim: sample space dimension
 * @LOOCV: whether to use LOOCV to calibrate the bandwidth
 * 
 * Creates a new #NcmStatsDistNdKDEGauss object with sample dimension @dim.
 * 
 * Returns: a new #NcmStatsDistNdKDEGauss.
 */
NcmStatsDistNdKDEGauss *
ncm_stats_dist_nd_kde_gauss_new (const guint dim, const gboolean LOOCV)
{
  NcmStatsDistNdKDEGauss *dndg = g_object_new (NCM_TYPE_STATS_DIST_ND_KDE_GAUSS,
                                               "dimension", dim,
                                               "LOOCV",     LOOCV,
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
 * ncm_stats_dist_nd_kde_gauss_set_over_smooth:
 * @dndg: a #NcmStatsDistNdKDEGauss
 * @over_smooth: the over-smooth factor
 *
 * Sets the over-smooth factor to @over_smooth.
 * 
 */
void 
ncm_stats_dist_nd_kde_gauss_set_over_smooth (NcmStatsDistNdKDEGauss *dndg, const gdouble over_smooth)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  self->over_smooth = over_smooth;
}

/**
 * ncm_stats_dist_nd_kde_gauss_get_over_smooth:
 * @dndg: a #NcmStatsDistNdKDEGauss
 *
 * Returns: the over-smooth factor.
 */
gdouble 
ncm_stats_dist_nd_kde_gauss_get_over_smooth (NcmStatsDistNdKDEGauss *dndg)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  return self->over_smooth;
}

/**
 * ncm_stats_dist_nd_kde_gauss_set_LOOCV_bandwidth_adj:
 * @dndg: a #NcmStatsDistNdKDEGauss
 * @LOOCV: whether to use LOOCV to adjust the bandwidth
 *
 * Sets whether to use LOOCV to adjust the bandwidth.
 * 
 */
void 
ncm_stats_dist_nd_kde_gauss_set_LOOCV_bandwidth_adj (NcmStatsDistNdKDEGauss *dndg, gboolean LOOCV)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  self->LOOCV = LOOCV;
}

/**
 * ncm_stats_dist_nd_kde_gauss_get_LOOCV_bandwidth_adj:
 * @dndg: a #NcmStatsDistNdKDEGauss
 *
 * Returns: true if the objects uses LOOCV to adjust the bandwidth.
 */
gboolean 
ncm_stats_dist_nd_kde_gauss_get_LOOCV_bandwidth_adj (NcmStatsDistNdKDEGauss *dndg)
{
  NcmStatsDistNdKDEGaussPrivate * const self = dndg->priv;
  return self->LOOCV;
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


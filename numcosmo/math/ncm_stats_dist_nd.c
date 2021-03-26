/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd.c
 *
 *  Wed November 07 16:02:36 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd.c
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
 * SECTION:ncm_stats_dist_nd
 * @title: NcmStatsDistNd
 * @short_description: Abstract class for implementing N dimensional probability distributions
 *
 * Abstract class to reconstruct an arbitrary N dimensional probability distribution.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_nd.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include "gslextras/cqp/gsl_cqp.h"
#include "libqp.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsDistNdPrivate
{
  NcmStatsVec *sample;
  NcmMatrix *cov_decomp;
  NcmMatrix *log_cov;
  NcmMatrix *sample_matrix;
  GPtrArray *invUsample;
  NcmVector *weights;
  NcmVector *v;
  gdouble over_smooth;
  NcmStatsDistNdCV cv_type;
  gdouble href;
  gdouble href2;
  gdouble kernel_lnnorm;
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
};

enum
{
  PROP_0,
  PROP_DIM,
  PROP_SAMPLE_SIZE,
  PROP_OVER_SMOOTH,
  PROP_NEARPD_MAXITER,
  PROP_CV_TYPE,

};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDistNd, ncm_stats_dist_nd, G_TYPE_OBJECT);

static void
ncm_stats_dist_nd_init (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv = ncm_stats_dist_nd_get_instance_private (dnd);

  self->sample           = NULL;
  self->cov_decomp       = NULL;
  self->log_cov          = NULL;
  self->sample_matrix    = NULL;
  self->invUsample       = g_ptr_array_new ();
  self->weights          = NULL;
  self->v                = NULL;
  self->over_smooth      = 0.0;
  self->cv_type          = NCM_STATS_DIST_ND_CV_LEN;
  self->href             = 0.0;
  self->href2            = 0.0;
  self->kernel_lnnorm    = 0.0;
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

  g_ptr_array_set_free_func (self->invUsample, (GDestroyNotify) ncm_vector_free);
}

static void _ncm_stats_dist_nd_set_dim (NcmStatsDistNd *dnd, const guint dim);

static void
_ncm_stats_dist_nd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNd *dnd = NCM_STATS_DIST_ND (object);
  /*g_return_if_fail (NCM_IS_STATS_DIST_ND (object));*/

  switch (prop_id)
  {
    case PROP_DIM:
      _ncm_stats_dist_nd_set_dim (dnd, g_value_get_uint (value));
      break;
    case PROP_SAMPLE_SIZE:
      g_assert_not_reached ();
      break;
    case PROP_OVER_SMOOTH:
      ncm_stats_dist_nd_set_over_smooth (dnd, g_value_get_double (value));
      break;
    case PROP_NEARPD_MAXITER:
      ncm_stats_dist_nd_set_nearPD_maxiter (dnd, g_value_get_uint (value));
      break;
    case PROP_CV_TYPE:
      ncm_stats_dist_nd_set_cv_type (dnd, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNd *dnd = NCM_STATS_DIST_ND (object);
  NcmStatsDistNdPrivate * const self = dnd->priv;
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND (object));

  switch (prop_id)
  {
    case PROP_DIM:
      g_value_set_uint (value, self->d);
      break;
    case PROP_SAMPLE_SIZE:
      g_value_set_uint (value, ncm_stats_vec_nitens (self->sample));
      break;
    case PROP_OVER_SMOOTH:
      g_value_set_double (value, ncm_stats_dist_nd_get_over_smooth (dnd));
      break;
    case PROP_NEARPD_MAXITER:
      g_value_set_uint (value, self->nearPD_maxiter);
      break;
    case PROP_CV_TYPE:
      g_value_set_enum (value, ncm_stats_dist_nd_get_cv_type (dnd));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_dispose (GObject *object)
{
  NcmStatsDistNd *dnd = NCM_STATS_DIST_ND (object);
  NcmStatsDistNdPrivate * const self = dnd->priv;

  ncm_stats_vec_clear (&self->sample);
  ncm_matrix_clear (&self->cov_decomp);
  ncm_matrix_clear (&self->log_cov);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_vector_clear (&self->weights);
  ncm_vector_clear (&self->v);

  g_clear_pointer (&self->sampling, g_array_unref);
  g_clear_pointer (&self->invUsample, g_ptr_array_unref);

  ncm_matrix_clear (&self->A);
  ncm_matrix_clear (&self->IM);
  ncm_vector_clear (&self->b);
  ncm_vector_clear (&self->beta);
  ncm_vector_clear (&self->dv);
  ncm_vector_clear (&self->alpha);
  ncm_vector_clear (&self->xi);
  ncm_vector_clear (&self->zeta);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_parent_class)->finalize (object);
}

gdouble _ncm_stats_dist_nd_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n) { g_error ("method get_rot_bandwidth not implemented by %s.", G_OBJECT_TYPE_NAME (dnd)); return 0.0; }
gdouble _ncm_stats_dist_nd_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const gdouble href) { g_error ("method get_kernel_lnnorm not implemented by %s.", G_OBJECT_TYPE_NAME (dnd)); return 0.0; }
static void _ncm_stats_dist_nd_prepare_kernel_args (NcmStatsDistNd *dnd, NcmStatsVec *sample);
static void _ncm_stats_dist_nd_prepare (NcmStatsDistNd *dnd);
static void _ncm_stats_dist_nd_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp);
static void _ncm_stats_dist_nd_reset (NcmStatsDistNd *dnd);

static void
ncm_stats_dist_nd_class_init (NcmStatsDistNdClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  

  object_class->set_property = &_ncm_stats_dist_nd_set_property;
  object_class->get_property = &_ncm_stats_dist_nd_get_property;
  object_class->dispose      = &_ncm_stats_dist_nd_dispose;
  object_class->finalize     = &_ncm_stats_dist_nd_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dimension",
                                                      NULL,
                                                      "PDF dimension",
                                                      2, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SAMPLE_SIZE,
                                   g_param_spec_uint ("N",
                                                      NULL,
                                                      "sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_OVER_SMOOTH,
                                   g_param_spec_double ("over-smooth",
                                                        NULL,
                                                        "Over-smooth distribution",
                                                        1.0e-5, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_CV_TYPE,
                                   g_param_spec_enum ("CV-type",
                                                      NULL,
                                                      "Cross-validation method",
                                                      NCM_TYPE_STATS_DIST_ND_CV, NCM_STATS_DIST_ND_CV_NONE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NEARPD_MAXITER,
                                   g_param_spec_uint ("nearPD-maxiter",
                                                      NULL,
                                                      "Maximum number of iterations in the nearPD call",
                                                      1, G_MAXUINT, 200,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_dim             = &_ncm_stats_dist_nd_set_dim;
  klass->get_rot_bandwidth   = &_ncm_stats_dist_nd_get_rot_bandwidth;
  klass->get_kernel_lnnorm   = &_ncm_stats_dist_nd_get_kernel_lnnorm;
  klass->prepare_kernel_args = &_ncm_stats_dist_nd_prepare_kernel_args;
  klass->prepare             = &_ncm_stats_dist_nd_prepare;
  klass->prepare_interp      = &_ncm_stats_dist_nd_prepare_interp;
  klass->eval                = NULL;
  klass->kernel_sample       = NULL;
  klass->kernel_eval_m2lnp   = NULL;
  klass->reset               = &_ncm_stats_dist_nd_reset;
}

static void
_ncm_stats_dist_nd_set_dim (NcmStatsDistNd *dnd, const guint dim)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  ncm_stats_vec_clear (&self->sample);

  ncm_matrix_clear (&self->cov_decomp);
  ncm_matrix_clear (&self->log_cov);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_vector_clear (&self->v);

  self->d          = dim;
  self->sample     = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, TRUE);
  self->cov_decomp = ncm_matrix_new (dim, dim);
  self->log_cov    = ncm_matrix_new (dim, dim);
  self->v          = ncm_vector_new (dim);
}

static void
_ncm_stats_dist_nd_prepare_kernel_args (NcmStatsDistNd *dnd, NcmStatsVec *sample)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;
  gint i, ret;

  self->n = ncm_stats_vec_nitens (sample);
  g_assert_cmpuint (self->n, >, 1);

  ncm_matrix_memcpy (self->cov_decomp, ncm_stats_vec_peek_cov_matrix (sample, 0));
  if (ncm_matrix_cholesky_decomp (self->cov_decomp, 'U') != 0)
  {
    if (ncm_matrix_nearPD (self->cov_decomp, 'U', TRUE, self->nearPD_maxiter) != 0)
    {
      ncm_matrix_set_zero (self->cov_decomp);
      for (i = 0; i < self->d; i++)
      {
        ncm_matrix_set (self->cov_decomp, i, i, ncm_stats_vec_get_cov (sample, i, i));
      }
      g_assert_cmpint (ncm_matrix_cholesky_decomp (self->cov_decomp, 'U'), ==, 0);
    }
  }

  self->href          = self->over_smooth * ncm_stats_dist_nd_get_rot_bandwidth (dnd, self->d, self->n);
  self->href2         = self->href * self->href;
  self->kernel_lnnorm = ncm_stats_dist_nd_get_kernel_lnnorm (dnd, self->cov_decomp, self->d, self->n, self->href);

  if ((self->sample_matrix == NULL) || (self->n != ncm_matrix_nrows (self->sample_matrix)) || (self->d != ncm_matrix_ncols (self->sample_matrix)))
  {
    ncm_matrix_clear (&self->sample_matrix);
    self->sample_matrix = ncm_matrix_new (self->n, self->d);

    g_ptr_array_set_size (self->invUsample, 0);
    for (i = 0; i < self->n; i++)
    {
      NcmVector *row_i = ncm_matrix_get_row (self->sample_matrix, i);
      g_ptr_array_add (self->invUsample, row_i);
    }
  }

  for (i = 0; i < self->n; i++)
    ncm_matrix_set_row (self->sample_matrix, i, ncm_stats_vec_peek_row (self->sample, i));

  ret = gsl_blas_dtrsm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, ncm_matrix_gsl (self->cov_decomp), ncm_matrix_gsl (self->sample_matrix));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_prepare_kernel_args", ret);
}

static void
_ncm_stats_dist_nd_prepare (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd);
  NcmStatsDistNdPrivate * const self = dnd->priv;

  dnd_class->prepare_kernel_args (dnd, self->sample);

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

  ncm_vector_set_all (self->weights, 1.0 / (1.0 * self->n));
}

void LowRankQP (gint *n, gint *m, gint *p, gint *method, gint *verbose, gint *niter,
                gdouble *Q, gdouble *c, gdouble *A, gdouble *b, gdouble *u, gdouble *alpha,
                gdouble *beta, gdouble *xi, gdouble *zeta);

int nnls_c (double *a, const int *mda, const int *m, const int *n, double *b,
            double *x, double* rnorm, double* w, double* zz, int *index,
            int *mode);

static void
_ncm_stats_dist_nd_LowRankQP_solve_IMx_f (gint n, NcmMatrix *IM, NcmVector *x, NcmVector *f,
    NcmVector *b, NcmMatrix *A, NcmVector *ub, NcmVector *beta, NcmVector *xi, NcmVector *zeta)
{
  gint nc      = 0;
  gint method  = 2;
  gint verbose = 0;
  gint maxiter = 4000;

  /* Setting constraint value */
  ncm_vector_set (b, 0, 1.0);

  /* Setting constraint matrix */
  ncm_matrix_set_all (A, 1.0);

  /* Setting upper-bound */
  ncm_vector_set_all (ub, 1.0e20);

  /* Setting output vectors to zero */
  ncm_vector_set_zero (x);
  ncm_vector_set_zero (beta);
  ncm_vector_set_zero (xi);
  ncm_vector_set_zero (zeta);

  LowRankQP (&n, &n, &nc, &method, &verbose, &maxiter,
             ncm_matrix_data (IM), ncm_vector_data (f),
             ncm_matrix_data (A), ncm_vector_data (b), ncm_vector_data (ub),
             ncm_vector_data (x), ncm_vector_data (beta), ncm_vector_data (xi), ncm_vector_data (zeta));

}

static libqp_state_T
_ncm_stats_dist_nd_libqp_solve_IMx_f (gint n, NcmMatrix *IM, NcmVector *x, NcmVector *f,
    NcmVector *b, NcmMatrix *A, NcmVector *ub, NcmVector *beta, NcmVector *xi, NcmVector *zeta)
{
  GPtrArray *col = g_ptr_array_new ();
  uint32_t *II = g_new (uint32_t, n);
  gdouble bbb = 1.0;
  uint8_t S = 0;
  libqp_state_T res;
  gint i;

  for (i = 0; i < n; i++)
  {
    g_ptr_array_add (col, ncm_matrix_ptr (IM, i, 0));
  }

  ncm_vector_set_all (zeta, 1.0);
  ncm_vector_set_all (x, 0.0);
  ncm_vector_set_all (xi, GSL_POSINF);
  ncm_vector_set_zero (beta);

  for (i = 0; i < n; i++)
    II[i] = 1;

  res = libqp_splx_solver ((gdouble **)col->pdata,
                           ncm_vector_data (zeta),
                           ncm_vector_data (f),
                           &bbb,
                           II,
                           &S,
                           ncm_vector_data (x),
                           n, 10000000,
                           0.0, 1.0e-13, 1.0e-10, NULL);

  g_ptr_array_unref (col);
  g_free (II);

  return res;

}

static void
_ncm_stats_dist_nd_nnls_solve_IMx_f (gint n, NcmMatrix *IM, NcmVector *x, NcmVector *f,
    NcmVector *b, NcmMatrix *A, NcmVector *ub, NcmVector *beta, NcmVector *xi, NcmVector *zeta)
{
  gint nnn = n, mode = -31;
  gdouble rnorm;

  nnls_c (ncm_matrix_data (IM), &nnn, &nnn, &nnn, ncm_vector_data (f),
      ncm_vector_data (x), &rnorm, ncm_vector_data (zeta), ncm_vector_data (xi), (gint *)ncm_vector_data (beta),
      &mode);
}

enum {
  LOW_RANK_QP,
  LIBQP,
  NNLS,
};

static void
_ncm_stats_dist_nd_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd);
  NcmStatsDistNdPrivate * const self = dnd->priv;
  gint method = NNLS;
  gint i;

  ncm_stats_dist_nd_prepare (dnd);
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
    self->b       = ncm_vector_new (self->n);
    self->beta    = ncm_vector_new (self->n);
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

  switch (self->cv_type)
  {
    case NCM_STATS_DIST_ND_CV_NONE:
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  if (self->n > 10000)
    g_warning ("_ncm_stats_dist_nd_prepare_interp: very large system n = %u!", self->n);

  switch (method)
  {
    case LOW_RANK_QP:
      dnd_class->prepare_IM (dnd, self->invUsample, self->d, self->n, self->href, self->href2, self->IM);

      for (i = 0; i < self->n; i++)
      {
        const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
        ncm_vector_set (self->weights, i, - exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
      }

      _ncm_stats_dist_nd_LowRankQP_solve_IMx_f (self->n, self->IM, self->alpha, self->weights,
          self->b, self->A, self->dv, self->beta, self->xi, self->zeta);
      break;
    case LIBQP:
      dnd_class->prepare_IM (dnd, self->invUsample, self->d, self->n, self->href, self->href2, self->IM);

      for (i = 0; i < self->n; i++)
      {
        const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
        ncm_vector_set (self->weights, i, - exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
      }

      _ncm_stats_dist_nd_libqp_solve_IMx_f (self->n, self->IM, self->alpha, self->weights,
          self->b, self->A, self->dv, self->beta, self->xi, self->zeta);
      break;
    case NNLS:
      dnd_class->prepare_IM (dnd, self->invUsample, self->d, self->n, self->href, self->href2, self->IM);

      for (i = 0; i < self->n; i++)
      {
        const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);
        ncm_vector_set (self->weights, i, + exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
      }

      _ncm_stats_dist_nd_nnls_solve_IMx_f (self->n, self->IM, self->alpha, self->weights,
          self->b, self->A, self->dv, self->beta, self->xi, self->zeta);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  ncm_vector_memcpy (self->weights, self->alpha);
  /*ncm_vector_scale (self->weights, 1.0 / ncm_vector_sum_cpts (self->weights));*/
}

static void
_ncm_stats_dist_nd_reset (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  ncm_stats_vec_reset (self->sample, TRUE);
}

/**
 * ncm_stats_dist_nd_ref:
 * @dnd: a #NcmStatsDistNd
 *
 * Increases the reference count of @dnd.
 * 
 * Returns: (transfer full): @dnd.
 */
NcmStatsDistNd *
ncm_stats_dist_nd_ref (NcmStatsDistNd *dnd)
{
  return g_object_ref (dnd);
}

/**
 * ncm_stats_dist_nd_free:
 * @dnd: a #NcmStatsDistNd
 *
 * Decreases the reference count of @dnd.
 *
 */
void 
ncm_stats_dist_nd_free (NcmStatsDistNd *dnd)
{
  g_object_unref (dnd);
}

/**
 * ncm_stats_dist_nd_clear:
 * @dnd: a #NcmStatsDistNd
 *
 * Decreases the reference count of *@dnd and sets the pointer *@dnd to NULL.
 *
 */
void 
ncm_stats_dist_nd_clear (NcmStatsDistNd **dnd)
{
  g_clear_object (dnd);
}

/**
 * ncm_stats_dist_nd_get_dim:
 * @dnd: a #NcmStatsDistNd
 *
 * Returns: the dimension of the sample space.
 */
guint 
ncm_stats_dist_nd_get_dim (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  return self->d;
}

/**
 * ncm_stats_dist_nd_get_rot_bandwidth: (virtual get_rot_bandwidth)
 * @dnd: a #NcmStatsDistNd
 * @d: problem dimension
 * @n: sample size
 *
 * Returns: the rule-of-thumb bandwidth estimate.
 */
gdouble
ncm_stats_dist_nd_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd);

  return dnd_class->get_rot_bandwidth (dnd, d, n);
}

/**
 * ncm_stats_dist_nd_get_kernel_lnnorm: (virtual get_kernel_lnnorm)
 * @dnd: a #NcmStatsDistNd
 * @cov_decomp: Cholesky decomposition (U) of the covariance matrix
 * @d: problem dimension
 * @n: sample size
 * @href: interpolation bandwidth
 *
 * Returns: the log-normalization of a single kernel.
 */
gdouble
ncm_stats_dist_nd_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const gdouble href)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd);
  NcmStatsDistNdPrivate * const self = dnd->priv;

  return dnd_class->get_kernel_lnnorm (dnd, self->cov_decomp, d, n, self->href);
}

/**
 * ncm_stats_dist_nd_set_over_smooth:
 * @dnd: a #NcmStatsDistNd
 * @over_smooth: the over-smooth factor
 *
 * Sets the over-smooth factor to @over_smooth.
 *
 */
void
ncm_stats_dist_nd_set_over_smooth (NcmStatsDistNd *dnd, const gdouble over_smooth)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  self->over_smooth = over_smooth;
}

/**
 * ncm_stats_dist_nd_get_over_smooth:
 * @dnd: a #NcmStatsDistNd
 *
 * Returns: the over-smooth factor.
 */
gdouble
ncm_stats_dist_nd_get_over_smooth (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  return self->over_smooth;
}

/**
 * ncm_stats_dist_nd_set_nearPD_maxiter:
 * @dnd: a #NcmStatsDistNd
 * @maxiter: maximum number of iterations
 *
 * Sets the maximum number of iterations when finding the
 * nearest positive definite covariance matrix to @maxiter.
 *
 */
void
ncm_stats_dist_nd_set_nearPD_maxiter (NcmStatsDistNd *dnd, const guint maxiter)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  self->nearPD_maxiter = maxiter;
}

/**
 * ncm_stats_dist_nd_get_nearPD_maxiter:
 * @dnd: a #NcmStatsDistNd
 *
 * Returns: maximum number of iterations when finding the nearest positive definite covariance matrix.
 */
guint
ncm_stats_dist_nd_get_nearPD_maxiter (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  return self->nearPD_maxiter;
}

/**
 * ncm_stats_dist_nd_set_cv_type:
 * @dnd: a #NcmStatsDistNd
 * @cv_type: a #NcmStatsDistNdCV
 *
 * Sets the cross-validation method to @cv_type.
 *
 */
void
ncm_stats_dist_nd_set_cv_type (NcmStatsDistNd *dnd, const NcmStatsDistNdCV cv_type)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  self->cv_type = cv_type;
}

/**
 * ncm_stats_dist_nd_get_cv_type:
 * @dnd: a #NcmStatsDistNd
 *
 * Returns: current cross-validation method used.
 */
NcmStatsDistNdCV
ncm_stats_dist_nd_get_cv_type (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;

  return self->cv_type;
}

/**
 * ncm_stats_dist_nd_prepare: (virtual prepare)
 * @dnd: a #NcmStatsDistNd
 *
 * Prepares the object for calculations.
 */
void 
ncm_stats_dist_nd_prepare (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 

  if (dnd_class->prepare != NULL)
    dnd_class->prepare (dnd);
}

/**
 * ncm_stats_dist_nd_prepare_interp: (virtual prepare_interp)
 * @dnd: a #NcmStatsDistNd
 * @m2lnp: a #NcmVector containing the distribution values
 *
 * Prepares the object for calculations. Using the distribution values
 * at the sample points.
 * 
 */
void 
ncm_stats_dist_nd_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 

  g_assert (dnd_class->prepare_interp != NULL);
  dnd_class->prepare_interp (dnd, m2lnp);
}

/**
 * ncm_stats_dist_nd_eval:
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 *
 * Evaluate the distribution at $\vec{x}=$@x. If the distribution
 * was prepared using ncm_stats_dist_nd_prepare_interp(), the 
 * results will follow the interpolation and may not be properly 
 * normalized. In this case the method ncm_stats_dist_nd_eval_m2lnp()
 * should be used to avoid underflow.
 * 
 * Returns: $P(\vec{x})$.
 */
gdouble 
ncm_stats_dist_nd_eval (NcmStatsDistNd *dnd, NcmVector *x)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 
  NcmStatsDistNdPrivate * const self = dnd->priv;
  gint ret;

  ncm_vector_memcpy (self->v, x);
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("ncm_stats_dist_nd_eval", ret);

  return dnd_class->eval (dnd, self->weights, self->v, self->invUsample, self->d, self->n, self->href, self->href2) * exp (- self->kernel_lnnorm);
}

/**
 * ncm_stats_dist_nd_eval_m2lnp:
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 *
 * Evaluate the distribution at $\vec{x}=$@x. If the distribution
 * was prepared using ncm_stats_dist_nd_prepare_interp(), the 
 * results will follow the interpolation and may not be properly 
 * normalized.
 * 
 * Returns: $P(\vec{x})$.
 */
gdouble 
ncm_stats_dist_nd_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 
  NcmStatsDistNdPrivate * const self = dnd->priv;
  gint ret;

  ncm_vector_memcpy (self->v, x);
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                        ncm_matrix_gsl (self->cov_decomp), ncm_vector_gsl (self->v));
  NCM_TEST_GSL_RESULT ("ncm_stats_dist_nd_eval", ret);

  return -2.0 * (log (dnd_class->eval (dnd, self->weights, self->v, self->invUsample, self->d, self->n, self->href, self->href2)) - self->kernel_lnnorm);
}

/**
 * ncm_stats_dist_nd_sample:
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 * @rng: a #NcmRNG
 * 
 * Using the pseudo-random number generator @rng generates a 
 * point from the distribution and copy it to @x.
 * 
 */
void
ncm_stats_dist_nd_sample (NcmStatsDistNd *dnd, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistNdPrivate * const self = dnd->priv;
  const guint n = ncm_vector_len (self->weights);
  gint i;

  g_array_set_size (self->sampling, ncm_vector_len (self->weights));
  gsl_ran_multinomial (rng->r, n, 1, ncm_vector_data (self->weights), (guint *)self->sampling->data);

  for (i = 0; i < n; i++)
  {
    if (g_array_index (self->sampling, guint, i) > 0)
    {
      NcmVector *y_i = ncm_stats_vec_peek_row (self->sample, i);
      ncm_stats_dist_nd_kernel_sample (dnd, x, y_i, self->href, rng);
      break;
    }
  }
}

/**
 * ncm_stats_dist_nd_kernel_sample:
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 * @mu: a #NcmVector
 * @href: a double
 * @rng: a #NcmRNG
 * 
 * Using the pseudo-random number generator @rng generates a 
 * point from the distribution and copy it to @x.
 * 
 */
void
ncm_stats_dist_nd_kernel_sample (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *mu, const gdouble href, NcmRNG *rng)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 
  NcmStatsDistNdPrivate * const self = dnd->priv;

  dnd_class->kernel_sample (dnd, self->cov_decomp, self->d, x, mu, self->href, rng);
}

/**
 * ncm_stats_dist_nd_kernel_eval_m2lnp: (virtual kernel_eval_m2lnp)
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 * @y: a #NcmVector
 * @href: covariance href
 * 
 * Evaluates a single kernel at @x and @y and href @s, i.e., $K_s(x,y)$.
 * 
 * Returns: $K_s(x,y)$.
 */
gdouble
ncm_stats_dist_nd_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *y, const gdouble href)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd);
  NcmStatsDistNdPrivate * const self = dnd->priv;

  return dnd_class->kernel_eval_m2lnp (dnd, self->cov_decomp, self->d, x, y, self->v, href, href * href) + 2.0 * self->kernel_lnnorm;
}

/**
 * ncm_stats_dist_nd_add_obs_weight:
 * @dndg: a #NcmStatsDistNd
 * @y: a #NcmVector
 * @w: weight
 *
 * Adds a new point @y to the sample with weight @w.
 *
 */
void
ncm_stats_dist_nd_add_obs_weight (NcmStatsDistNd *dndg, NcmVector *y, const gdouble w)
{
  NcmStatsDistNdPrivate * const self = dndg->priv;
  ncm_stats_vec_append_weight (self->sample, y, w, TRUE);
}

/**
 * ncm_stats_dist_nd_add_obs:
 * @dndg: a #NcmStatsDistNd
 * @y: a #NcmVector
 *
 * Adds a new point @y to the sample with weight 1.0.
 *
 */
void
ncm_stats_dist_nd_add_obs (NcmStatsDistNd *dndg, NcmVector *y)
{
  ncm_stats_dist_nd_add_obs_weight (dndg, y, 1.0);
}

/**
 * ncm_stats_dist_nd_reset: (virtual reset)
 * @dnd: a #NcmStatsDistNd
 * 
 * Reset the object discarding all added points.
 * 
 */
void 
ncm_stats_dist_nd_reset (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd);
  dnd_class->reset (dnd);
}

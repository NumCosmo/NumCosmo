/***************************************************************************
 *            ncm_stats_dist_kde.c
 *
 *  Wed November 07 16:02:36 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kde.c
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
 * NcmStatsDistKDE:
 *
 * Fixed-bandwidth kernel estimator for #NcmStatsDist.
 *
 * Uses one covariance matrix $\Sigma$ for all sample points. After the covariance
 * matrix is computed,
 * the algorithm computes the Cholesky decomposition, that is \begin{align} \Sigma &=
 * AA^T ,\end{align} where $A$ is a triangular positive defined matrix and $A^T$ is its
 * transpose. The $A$ matrix is used in the least square squares calculation method that
 * is called in the #NcmStatsDist class.
 *
 *
 * The object also prepares the interpolation matrix to be implemented in the
 * least-squares problem, that is, given the relation
 *
 * $\left[\begin{array}{cccc}
 * \phi\left(\left\|\mathbf{x}_{1}-\mathbf{x}_{1}\right\|\right) & \phi\left(\left\|\mathbf{x}_{2}-\mathbf{x}_{1}\right\|\right) & \ldots & \phi\left(\left\|\mathbf{x}_{n}-\mathbf{x}_{1}\right\|\right) \newline
 * \phi\left(\left\|\mathbf{x}_{1}-\mathbf{x}_{2}\right\|\right) & \phi\left(\left\|\mathbf{x}_{2}-\mathbf{x}_{2}\right\|\right) & \ldots & \phi\left(\left\|\mathbf{x}_{n}-\mathbf{x}_{2}\right\|\right) \newline
 * \vdots & \vdots & & \vdots \newline
 * \phi\left(\left\|\mathbf{x}_{1}-\mathbf{x}_{n}\right\|\right) & \phi\left(\left\|\mathbf{x}_{2}-\mathbf{x}_{n}\right\|\right) & \ldots & \phi\left(\left\|\mathbf{x}_{n}-\mathbf{x}_{n}\right\|\right)
 * \end{array}\right]\left[\begin{array}{c}
 * \lambda_{1} \newline
 * \lambda_{2} \newline
 * \vdots \newline
 * \lambda_{n}
 * \end{array}\right]=\left[\begin{array}{c}
 * g_{1} \newline
 * g_{2} \newline
 * \vdots \newline
 * g_{n}
 * ,\end{array}\right]$
 *
 * which is explained in the #NcmStatsDist class, this object prepares the first matrix
 * for all the $n$ points in the sample, using the covariance matrix and the defined
 * kernel. The #NcmStatsDist class implements the solution for this relation and then
 * one can compute the distribution for a given vector $\vec{x}$ using a method of the
 * #NcmStatsDist class but that is implemented in this object.
 *
 *
 * The user must provide input the values: @sdk, @CV_type - ncm_stats_dist_kde_new(), @y
 * - ncm_stats_dist_add_obs(), @split_frac - ncm_stats_dist_set_split_frac(),
 * @over_smooth - ncm_stats_dist_set_over_smooth(), $v(x)$ -
 * ncm_stats_dist_prepare(). To see an example of how to use this object and the
 * main functions that are called within each function, check the flowchart at the end
 * of this documentation, where the order of the functions that should be called by the
 * user and some of the functions that the algorithm calls.
 *
 * ![kde_sketch](kde.png)
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist_kde.h"
#include "ncm/core/ncm_iset.h"
#include "ncm/algebra/ncm_nnls.h"
#include "ncm/algebra/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#endif /* NUMCOSMO_GIR_SCAN */

#include "ncm/stats/ncm_stats_dist_kernel_gauss.h"
#include "ncm/stats/ncm_stats_dist_kde_private.h"
#include "ncm/stats/ncm_stats_dist_private.h"

enum
{
  PROP_0,
  PROP_NEARPD_MAXITER,
  PROP_COV_TYPE,
  PROP_COV_FIXED,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistKDE, ncm_stats_dist_kde, NCM_TYPE_STATS_DIST)

static NcmStatsDistPrivate *
ncm_stats_dist_get_instance_private (NcmStatsDist * sd)
{
  return g_type_instance_get_private ((GTypeInstance *) sd, NCM_TYPE_STATS_DIST);
}

typedef struct _NcmStatsDistKDEEvalVars
{
  NcmVector *v;
  NcmVector *chi2;
  NcmVector *lnK;
  NcmVector *lnc;
} NcmStatsDistKDEEvalVars;

static gpointer
_ncm_stats_dist_kde_eval_vars_new (gpointer userdata)
{
  NcmStatsDist *sd                   = NCM_STATS_DIST (userdata);
  NcmStatsDistPrivate * const ppself = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistKDEEvalVars *ev        = g_new0 (NcmStatsDistKDEEvalVars, 1);

  ev->v    = ncm_vector_new (ppself->d);
  ev->chi2 = ncm_vector_new (ppself->n_kernels);
  ev->lnK  = ncm_vector_new (ppself->n_kernels);
  ev->lnc  = ncm_vector_new (ppself->n_kernels);

  return ev;
}

static void
_ncm_stats_dist_kde_eval_vars_free (gpointer userdata)
{
  NcmStatsDistKDEEvalVars *ev = (NcmStatsDistKDEEvalVars *) userdata;

  ncm_vector_free (ev->v);
  ncm_vector_free (ev->chi2);
  ncm_vector_free (ev->lnK);
  ncm_vector_free (ev->lnc);

  g_free (ev);
}

static void
ncm_stats_dist_kde_init (NcmStatsDistKDE *sdkde)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  self->sample            = NULL;
  self->cov_type          = NCM_STATS_DIST_KDE_COV_TYPE_LEN;
  self->cov               = NULL;
  self->cov_fixed         = NULL;
  self->cov_decomp        = NULL;
  self->cov_decomp0       = NULL;
  self->sample_matrix     = NULL;
  self->invUsample_matrix = NULL;
  self->invUsample_array  = g_ptr_array_new ();
  self->center_matrix     = NULL;
  self->invUcenter_matrix = NULL;
  self->invUcenter_array  = g_ptr_array_new ();
  self->kernel_lnnorm     = 0.0;
  self->nearPD_maxiter    = 0;

  self->mp_eval_vars = ncm_memory_pool_new (&_ncm_stats_dist_kde_eval_vars_new, sdkde,
                                            &_ncm_stats_dist_kde_eval_vars_free);
  self->mp_eval_vars_len = 0;


  g_ptr_array_set_free_func (self->invUsample_array, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->invUcenter_array, (GDestroyNotify) ncm_vector_free);
}

static void
_ncm_stats_dist_kde_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKDE *sdkde = NCM_STATS_DIST_KDE (object);

  /*g_return_if_fail (NCM_IS_STATS_DIST_KDE (object));*/

  switch (prop_id)
  {
    case PROP_NEARPD_MAXITER:
      ncm_stats_dist_kde_set_nearPD_maxiter (sdkde, g_value_get_uint (value));
      break;
    case PROP_COV_TYPE:
      ncm_stats_dist_kde_set_cov_type (sdkde, g_value_get_enum (value));
      break;
    case PROP_COV_FIXED:
      ncm_stats_dist_kde_set_cov_fixed (sdkde, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_kde_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistKDE *sdkde = NCM_STATS_DIST_KDE (object);

  /*NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);*/

  g_return_if_fail (NCM_IS_STATS_DIST_KDE (object));

  switch (prop_id)
  {
    case PROP_NEARPD_MAXITER:
      g_value_set_uint (value, ncm_stats_dist_kde_get_nearPD_maxiter (sdkde));
      break;
    case PROP_COV_TYPE:
      g_value_set_enum (value, ncm_stats_dist_kde_get_cov_type (sdkde));
      break;
    case PROP_COV_FIXED:
      g_value_set_object (value, ncm_stats_dist_kde_peek_cov_fixed (sdkde));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_kde_dispose (GObject *object)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (object);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  ncm_stats_vec_clear (&self->sample);
  ncm_matrix_clear (&self->cov);
  ncm_matrix_clear (&self->cov_fixed);
  ncm_matrix_clear (&self->cov_decomp);
  ncm_matrix_clear (&self->cov_decomp0);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_matrix_clear (&self->invUsample_matrix);
  ncm_matrix_clear (&self->center_matrix);
  ncm_matrix_clear (&self->invUcenter_matrix);

  g_clear_pointer (&self->invUsample_array, g_ptr_array_unref);
  g_clear_pointer (&self->invUcenter_array, g_ptr_array_unref);

  if (self->mp_eval_vars)
  {
    ncm_memory_pool_free (self->mp_eval_vars, TRUE);
    self->mp_eval_vars = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kde_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_kde_finalize (GObject *object)
{
  /*NcmStatsDistKDE *sdkde = NCM_STATS_DIST_KDE (object);*/
  /*NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);*/

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_kde_parent_class)->finalize (object);
}

static void _ncm_stats_dist_kde_set_dim (NcmStatsDist *sd, const guint dim);
static void _ncm_stats_dist_kde_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array);
static void _ncm_stats_dist_kde_prepare_kernels (NcmStatsDist *sd);
static void _ncm_stats_dist_kde_compute_IM (NcmStatsDist *sd, NcmMatrix *IM);
static gdouble _ncm_stats_dist_kde_amise (NcmStatsDist *sd);
static NcmMatrix *_ncm_stats_dist_kde_peek_cov_decomp (NcmStatsDist *sd, guint i);
static NcmMatrix *_ncm_stats_dist_kde_peek_full_cov_decomp (NcmStatsDist *sd);
static NcmMatrix *_ncm_stats_dist_kde_peek_full_cov (NcmStatsDist *sd);
static gdouble _ncm_stats_dist_kde_get_lnnorm (NcmStatsDist *sd, guint i);
static gdouble _ncm_stats_dist_kde_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
static gdouble _ncm_stats_dist_kde_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);

static void
ncm_stats_dist_kde_class_init (NcmStatsDistKDEClass *klass)
{
  GObjectClass *object_class  = G_OBJECT_CLASS (klass);
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_kde_set_property;
  object_class->get_property = &_ncm_stats_dist_kde_get_property;
  object_class->dispose      = &_ncm_stats_dist_kde_dispose;
  object_class->finalize     = &_ncm_stats_dist_kde_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NEARPD_MAXITER,
                                   g_param_spec_uint ("nearPD-maxiter",
                                                      NULL,
                                                      "Maximum number of iterations in the nearPD call",
                                                      1, G_MAXUINT, 200,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COV_TYPE,
                                   g_param_spec_enum ("cov-type",
                                                      NULL,
                                                      "Covariance type",
                                                      NCM_TYPE_STATS_DIST_KDE_COV_TYPE, NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COV_FIXED,
                                   g_param_spec_object ("cov-fixed",
                                                        NULL,
                                                        "Fixed covariance matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  sd_class->set_dim              = &_ncm_stats_dist_kde_set_dim;
  sd_class->prepare_shapes       = &_ncm_stats_dist_kde_prepare_shapes;
  sd_class->prepare_kernels      = &_ncm_stats_dist_kde_prepare_kernels;
  sd_class->compute_IM           = &_ncm_stats_dist_kde_compute_IM;
  sd_class->amise                = &_ncm_stats_dist_kde_amise;
  sd_class->peek_cov_decomp      = &_ncm_stats_dist_kde_peek_cov_decomp;
  sd_class->peek_full_cov_decomp = &_ncm_stats_dist_kde_peek_full_cov_decomp;
  sd_class->peek_full_cov        = &_ncm_stats_dist_kde_peek_full_cov;
  sd_class->get_lnnorm           = &_ncm_stats_dist_kde_get_lnnorm;
  sd_class->eval_weights         = &_ncm_stats_dist_kde_eval_weights;
  sd_class->eval_weights_m2lnp   = &_ncm_stats_dist_kde_eval_weights_m2lnp;
}

static void
_ncm_stats_dist_kde_set_dim (NcmStatsDist *sd, const guint dim)
{
  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_kde_parent_class)->set_dim (sd, dim);
  {
    NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
    NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

    ncm_stats_vec_clear (&self->sample);

    ncm_matrix_clear (&self->cov);
    ncm_matrix_clear (&self->cov_decomp);
    ncm_matrix_clear (&self->cov_decomp0);
    ncm_matrix_clear (&self->sample_matrix);
    ncm_matrix_clear (&self->invUsample_matrix);
    ncm_matrix_clear (&self->center_matrix);
    ncm_matrix_clear (&self->invUcenter_matrix);
    g_ptr_array_set_size (self->invUcenter_array, 0);

    self->sample      = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, TRUE);
    self->cov         = ncm_matrix_new (dim, dim);
    self->cov_decomp  = ncm_matrix_new (dim, dim);
    self->cov_decomp0 = ncm_matrix_new (dim, dim);
  }
}

static void
_ncm_stats_dist_kde_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);
  guint i;

  /*
   * Computing the covariance matrix considering the whole sample.
   */
  ncm_stats_vec_reset (self->sample, TRUE);

  for (i = 0; i < pself->n_kernels; i++)
  {
    NcmVector *theta_i = g_ptr_array_index (sample_array, i);

    ncm_stats_vec_append (self->sample, theta_i, FALSE);
  }

  /*
   * Computing the Cholesky decomposition of the total covariance.
   */
  /*ncm_matrix_log_vals (ncm_stats_vec_peek_cov_matrix (self->sample, 0), "COV: ", "%12.5g");*/

  switch (self->cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    {
      NcmMatrix *cov = ncm_stats_vec_peek_cov_matrix (self->sample, 0);

      _ncm_stats_dist_cholesky (self->cov_decomp, cov, self->nearPD_maxiter, "the sample covariance");
      ncm_matrix_memcpy (self->cov, cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      if (self->cov_fixed == NULL)
        g_error ("_ncm_stats_dist_kde_prepare_shapes: cov_type is FIXED but a fixed "
                 "covariance matrix was not provided, use ncm_stats_dist_kde_set_cov_fixed to set one.");

      ncm_matrix_memcpy (self->cov, self->cov_fixed);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
    {
      NcmMatrix *cov = ncm_stats_vec_compute_cov_robust_diag (self->sample);

      _ncm_stats_dist_cholesky (self->cov_decomp, cov, self->nearPD_maxiter, "the robust diagonal covariance");
      ncm_matrix_memcpy (self->cov, cov);
      ncm_matrix_free (cov);

      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    {
      NcmMatrix *cov = ncm_stats_vec_compute_cov_robust_ogk (self->sample);

      _ncm_stats_dist_cholesky (self->cov_decomp, cov, self->nearPD_maxiter, "the robust covariance");
      ncm_matrix_memcpy (self->cov, cov);
      ncm_matrix_free (cov);
    }
    break;
    default:
      g_assert_not_reached ();
      break;
  }

  if ((self->sample_matrix == NULL) ||
      (pself->n_obs != ncm_matrix_nrows (self->sample_matrix)) ||
      (pself->d != ncm_matrix_ncols (self->sample_matrix)))
  {
    ncm_matrix_clear (&self->sample_matrix);
    ncm_matrix_clear (&self->invUsample_matrix);

    self->sample_matrix     = ncm_matrix_new (pself->n_obs, pself->d);
    self->invUsample_matrix = ncm_matrix_new (pself->n_obs, pself->d);

    g_ptr_array_set_size (self->invUsample_array, 0);

    for (i = 0; i < pself->n_obs; i++)
    {
      NcmVector *row_i = ncm_matrix_get_row (self->invUsample_matrix, i);

      g_ptr_array_add (self->invUsample_array, row_i);
    }
  }

  for (i = 0; i < pself->n_obs; i++)
    ncm_matrix_set_row (self->sample_matrix, i, g_ptr_array_index (sample_array, i));

  ncm_matrix_memcpy (self->invUsample_matrix, self->sample_matrix);

  ncm_matrix_dtrsm (self->invUsample_matrix, 'R', 'U', 'N', 1.0, self->cov_decomp);

  if ((self->center_matrix == NULL) ||
      (pself->n_kernels != ncm_matrix_nrows (self->center_matrix)) ||
      (pself->d != ncm_matrix_ncols (self->center_matrix)))
  {
    ncm_matrix_clear (&self->center_matrix);
    ncm_matrix_clear (&self->invUcenter_matrix);
    self->center_matrix     = ncm_matrix_new (pself->n_kernels, pself->d);
    self->invUcenter_matrix = ncm_matrix_new (pself->n_kernels, pself->d);

    g_ptr_array_set_size (self->invUcenter_array, 0);

    for (i = 0; i < pself->n_kernels; i++)
    {
      NcmVector *row_i = ncm_matrix_get_row (self->invUcenter_matrix, i);

      g_ptr_array_add (self->invUcenter_array, row_i);
    }
  }

  /*
   * What center shrinkage needs: the factor of the sample covariance and the kernel
   * scale matrix; NcmStatsDistVKDE replaces the latter with the mean over its kernels.
   * The untransformed factor is kept, the applied one follows the bandwidth.
   */
  {
    NcmMatrix *C_decomp = pself->sample_decomp;
    NcmMatrix *mean_cov = pself->kernel_cov;

    _ncm_stats_dist_cholesky (C_decomp, ncm_stats_vec_peek_cov_matrix (self->sample, 0), self->nearPD_maxiter, "the sample covariance");
    ncm_matrix_memcpy (mean_cov, self->cov);

    ncm_matrix_memcpy (self->cov_decomp0, self->cov_decomp);
  }

  /*
   * Allocating the evaluation vector
   */
  if (self->mp_eval_vars_len != pself->n_kernels)
  {
    ncm_memory_pool_empty (self->mp_eval_vars, TRUE);
    self->mp_eval_vars_len = pself->n_kernels;
  }
}

static void
_ncm_stats_dist_kde_compute_IM (NcmStatsDist *sd, NcmMatrix *IM)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble href2                 = pself->href * pself->href;
  guint i;

  /*
   * Rows are the observation points, columns the kernel centers. With center shrinkage
   * the two differ and the whole block has to be computed. Without it the centers are
   * the sample points themselves, so the first n_kernels x n_kernels block is symmetric
   * and only half of it is worth computing.
   */
  if (pself->shrink.on)
  {
    for (i = 0; i < pself->n_obs; i++)
    {
      NcmVector *row_i = g_ptr_array_index (self->invUsample_array, i);
      guint j;

      for (j = 0; j < pself->n_kernels; j++)
      {
        NcmVector *center_j = g_ptr_array_index (self->invUcenter_array, j);
        gdouble chi2_ij     = ncm_vector_sqr_dist (row_i, center_j);

        ncm_matrix_set (IM, i, j, chi2_ij / href2);
      }
    }
  }
  else
  {
    for (i = 0; i < pself->n_kernels; i++)
    {
      NcmVector *row_i = g_ptr_array_index (self->invUsample_array, i);
      guint j;

      ncm_matrix_set (IM, i, i, 0.0);

      for (j = i + 1; j < pself->n_kernels; j++)
      {
        NcmVector *row_j = g_ptr_array_index (self->invUsample_array, j);
        gdouble chi2_ij  = ncm_vector_sqr_dist (row_i, row_j);

        chi2_ij = chi2_ij / href2;

        ncm_matrix_set (IM, i, j, chi2_ij);
        ncm_matrix_set (IM, j, i, chi2_ij);
      }
    }

    for (i = pself->n_kernels; i < pself->n_obs; i++)
    {
      NcmVector *row_i = g_ptr_array_index (self->invUsample_array, i);
      guint j;

      for (j = 0; j < pself->n_kernels; j++)
      {
        NcmVector *row_j = g_ptr_array_index (self->invUsample_array, j);
        gdouble chi2_ij  = ncm_vector_sqr_dist (row_i, row_j);

        ncm_matrix_set (IM, i, j, chi2_ij / href2);
      }
    }
  }

  {
    /* One view moved along the rows, rather than a vector allocated for each of them. */
    NcmVector *row_i = ncm_vector_new_data_static (ncm_matrix_ptr (IM, 0, 0), ncm_matrix_ncols (IM), 1);

    for (i = 0; i < pself->n_obs; i++)
    {
      ncm_vector_replace_data (row_i, ncm_matrix_ptr (IM, i, 0));
      ncm_stats_dist_kernel_eval_unnorm_vec (pself->kernel, row_i, row_i);
    }

    ncm_vector_free (row_i);
  }

  ncm_matrix_scale (IM, exp (-(self->kernel_lnnorm + pself->d * log (pself->href))));
}

/*
 * The AMISE objective for one Gaussian bandwidth in closed form. The integral of the
 * squared estimate is the mean of the interpolation matrix at sqrt(2) times the
 * bandwidth (the convolution of two Gaussians); the leave-one-out cross term is the
 * off-diagonal mean of the matrix at the bandwidth itself, which is also where the object
 * is left. Any other kernel takes the base class's Monte Carlo estimate. NcmStatsDistVKDE
 * inherits this; there the sqrt(2) rule is an approximation, since two kernels with
 * different covariances do not convolve to either at sqrt(2) h -- the behaviour the
 * closed form has always had.
 */
static gdouble
_ncm_stats_dist_kde_amise (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const pself = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistClass *sd_class       = NCM_STATS_DIST_GET_CLASS (sd);
  const gdouble href                = pself->href;
  const gdouble n                   = pself->n_kernels;
  gdouble amise                     = 0.0;
  guint i, j;

  if (!NCM_IS_STATS_DIST_KERNEL_GAUSS (pself->kernel))
    return _ncm_stats_dist_amise (sd);

  pself->href = sqrt (2.0) * href;
  sd_class->compute_IM (sd, pself->IM);

  for (i = 0; i < pself->n_kernels; i++)
    for (j = 0; j < pself->n_kernels; j++)
      amise += ncm_matrix_get (pself->IM, i, j) / (n * n);

  pself->href = href;
  sd_class->compute_IM (sd, pself->IM);

  for (i = 0; i < pself->n_kernels; i++)
    for (j = 0; j < pself->n_kernels; j++)
      if (j != i)
        amise -= 2.0 * ncm_matrix_get (pself->IM, i, j) / (n * (n - 1.0));

  if (pself->print_fit)
    ncm_message ("# over-smooth: % 22.15g, amise = % 22.15g\n", pself->over_smooth, amise);

  return amise;
}

static NcmMatrix *
_ncm_stats_dist_kde_peek_cov_decomp (NcmStatsDist *sd, guint i)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  return self->cov_decomp;
}

static NcmMatrix *
_ncm_stats_dist_kde_peek_full_cov_decomp (NcmStatsDist *sd)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  return self->cov_decomp;
}

static NcmMatrix *
_ncm_stats_dist_kde_peek_full_cov (NcmStatsDist *sd)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  return self->cov;
}

static gdouble
_ncm_stats_dist_kde_get_lnnorm (NcmStatsDist *sd, guint i)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);

  return self->kernel_lnnorm + pself->d * log (pself->href);
}

static gdouble
_ncm_stats_dist_kde_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble href2                 = pself->href * pself->href;
  NcmStatsDistKDEEvalVars **ev_ptr    = ncm_memory_pool_get (self->mp_eval_vars);
  NcmStatsDistKDEEvalVars *ev         = *ev_ptr;
  gdouble res;
  guint i;

  ncm_vector_memcpy (ev->v, x);
  ncm_matrix_dtrsv (self->cov_decomp, 'U', 'T', ev->v);

  for (i = 0; i < pself->n_kernels; i++)
  {
    NcmVector *row_i = g_ptr_array_index (self->invUcenter_array, i);
    gdouble chi2_i   = ncm_vector_sqr_dist (row_i, ev->v);

    chi2_i = chi2_i / href2;

    ncm_vector_fast_set (ev->chi2, i, chi2_i);
  }

  ncm_stats_dist_kernel_eval_unnorm_vec (pself->kernel, ev->chi2, ev->chi2);

  res = ncm_vector_dot (ev->chi2, weights) * exp (-(self->kernel_lnnorm + pself->d * log (pself->href)));

  ncm_memory_pool_return (ev_ptr);

  return res;
}

static gdouble
_ncm_stats_dist_kde_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble href2                 = pself->href * pself->href;
  NcmStatsDistKDEEvalVars **ev_ptr    = ncm_memory_pool_get (self->mp_eval_vars);
  NcmStatsDistKDEEvalVars *ev         = *ev_ptr;
  guint i;

  ncm_vector_memcpy (ev->v, x);
  ncm_matrix_dtrsv (self->cov_decomp, 'U', 'T', ev->v);

  for (i = 0; i < pself->n_kernels; i++)
  {
    NcmVector *row_i = g_ptr_array_index (self->invUcenter_array, i);
    gdouble chi2_i   = ncm_vector_sqr_dist (row_i, ev->v);

    chi2_i = chi2_i / href2;

    ncm_vector_fast_set (ev->chi2, i, chi2_i);
  }

  {
    gdouble gamma, lambda;

    for (i = 0; i < pself->n_kernels; i++)
      ncm_vector_fast_set (ev->lnc, i, log (ncm_vector_get (weights, i)) - self->kernel_lnnorm);

    ncm_stats_dist_kernel_eval_gamma_lambda (pself->kernel, ev->chi2, ev->lnc, ev->lnK, &gamma, &lambda);

    ncm_memory_pool_return (ev_ptr);

    return -2.0 * (gamma + log1p (lambda) - pself->d * log (pself->href));
  }
}

static void
_ncm_stats_dist_kde_prepare_kernels (NcmStatsDist *sd)
{
  NcmStatsDistKDE *sdkde              = NCM_STATS_DIST_KDE (sd);
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);
  guint i;

  /*
   * The applied factor follows the current center transform; the whitened sample was
   * built from the untransformed one in prepare_shapes() and has to follow too.
   */
  _ncm_stats_dist_refactor_decomp (sd, self->cov_decomp0, self->cov_decomp);
  ncm_matrix_memcpy (self->invUsample_matrix, self->sample_matrix);
  ncm_matrix_dtrsm (self->invUsample_matrix, 'R', 'U', 'N', 1.0, self->cov_decomp);
  self->kernel_lnnorm = ncm_stats_dist_kernel_get_lnnorm (pself->kernel, self->cov_decomp);

  for (i = 0; i < pself->n_kernels; i++)
    ncm_matrix_set_row (self->center_matrix, i, g_ptr_array_index (pself->center_array, i));

  ncm_matrix_memcpy (self->invUcenter_matrix, self->center_matrix);
  ncm_matrix_dtrsm (self->invUcenter_matrix, 'R', 'U', 'N', 1.0, self->cov_decomp);
}

/**
 * ncm_stats_dist_kde_new:
 * @sdk: a #NcmStatsDistKernel
 * @CV_type: a #NcmStatsDistCV
 *
 * Creates a new #NcmStatsDistKDE object using @sdk as
 * kernel and @CV_type as cross-validation method.
 *
 * Returns: (transfer full): the newly created #NcmStatsDistKDE object.
 */
NcmStatsDistKDE *
ncm_stats_dist_kde_new (NcmStatsDistKernel *sdk, NcmStatsDistCV CV_type)
{
  NcmStatsDistKDE *sdkde = g_object_new (NCM_TYPE_STATS_DIST_KDE,
                                         "kernel", sdk,
                                         "CV-type", CV_type,
                                         NULL);

  return sdkde;
}

/**
 * ncm_stats_dist_kde_ref:
 * @sdkde: a #NcmStatsDistKDE
 *
 * Increases the reference count of @sdkde.
 *
 * Returns: (transfer full): @sdkde.
 */
NcmStatsDistKDE *
ncm_stats_dist_kde_ref (NcmStatsDistKDE *sdkde)
{
  return g_object_ref (sdkde);
}

/**
 * ncm_stats_dist_kde_free:
 * @sdkde: a #NcmStatsDistKDE
 *
 * Decreases the reference count of @sdkde.
 *
 */
void
ncm_stats_dist_kde_free (NcmStatsDistKDE *sdkde)
{
  g_object_unref (sdkde);
}

/**
 * ncm_stats_dist_kde_clear:
 * @sdkde: a #NcmStatsDistKDE
 *
 * Decreases the reference count of *@sdkde and sets the pointer *@sdkde to NULL.
 *
 */
void
ncm_stats_dist_kde_clear (NcmStatsDistKDE **sdkde)
{
  g_clear_object (sdkde);
}

/**
 * ncm_stats_dist_kde_set_nearPD_maxiter:
 * @sdkde: a #NcmStatsDistKDE
 * @maxiter: maximum number of iterations
 *
 * Sets the maximum number of iterations when finding the
 * nearest positive definite covariance matrix to @maxiter. This function is implemented
 * as a property and used wherever the object factors a covariance.
 *
 */
void
ncm_stats_dist_kde_set_nearPD_maxiter (NcmStatsDistKDE *sdkde, const guint maxiter)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  self->nearPD_maxiter = maxiter;
}

/**
 * ncm_stats_dist_kde_get_nearPD_maxiter:
 * @sdkde: a #NcmStatsDistKDE
 *
 * Returns:an int nearPD_maxiter, the maximum number of iterations when finding the nearest positive definite covariance matrix.
 */
guint
ncm_stats_dist_kde_get_nearPD_maxiter (NcmStatsDistKDE *sdkde)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  return self->nearPD_maxiter;
}

/**
 * ncm_stats_dist_kde_set_cov_type:
 * @sdkde: a #NcmStatsDistKDE
 * @cov_type: covariance type
 *
 * Sets the covariance type to use in kernel interpolation.
 *
 */
void
ncm_stats_dist_kde_set_cov_type (NcmStatsDistKDE *sdkde, NcmStatsDistKDECovType cov_type)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  self->cov_type = cov_type;

  if ((self->cov_type == NCM_STATS_DIST_KDE_COV_TYPE_FIXED) && (self->cov_fixed != NULL))
    _ncm_stats_dist_cholesky (self->cov_decomp, self->cov_fixed, self->nearPD_maxiter, "the fixed covariance");
}

/**
 * ncm_stats_dist_kde_get_cov_type:
 * @sdkde: a #NcmStatsDistKDE
 *
 * Returns: the covariance type #NcmStatsDistKDECovType.
 */
NcmStatsDistKDECovType
ncm_stats_dist_kde_get_cov_type (NcmStatsDistKDE *sdkde)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  return self->cov_type;
}

/**
 * ncm_stats_dist_kde_set_cov_fixed:
 * @sdkde: a #NcmStatsDistKDE
 * @cov_fixed: the fixed covariance matrix #NcmMatrix
 *
 * Sets the covariance matrix to be used when #NcmStatsDistKDECovType is
 * set to #NCM_STATS_DIST_KDE_COV_TYPE_FIXED. A copy of the matrix
 * @cov_fixed is made and saved into the object.
 *
 */
void
ncm_stats_dist_kde_set_cov_fixed (NcmStatsDistKDE *sdkde, NcmMatrix *cov_fixed)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);
  NcmStatsDist *sd                    = NCM_STATS_DIST (sdkde);
  NcmStatsDistPrivate * const pself   = ncm_stats_dist_get_instance_private (sd);

  g_assert_cmpuint (ncm_matrix_ncols (cov_fixed), ==, pself->d);
  g_assert_cmpuint (ncm_matrix_nrows (cov_fixed), ==, pself->d);

  ncm_matrix_clear (&self->cov_fixed);

  self->cov_fixed = ncm_matrix_dup (cov_fixed);

  if (self->cov_type == NCM_STATS_DIST_KDE_COV_TYPE_FIXED)
    _ncm_stats_dist_cholesky (self->cov_decomp, self->cov_fixed, self->nearPD_maxiter, "the fixed covariance");
}

/**
 * ncm_stats_dist_kde_peek_cov_fixed:
 * @sdkde: a #NcmStatsDistKDE
 *
 * Gets the currently used fixed covariance matrix.
 *
 * Returns: (transfer none) (allow-none): the fixed covariance matrix
 */
NcmMatrix *
ncm_stats_dist_kde_peek_cov_fixed (NcmStatsDistKDE *sdkde)
{
  NcmStatsDistKDEPrivate * const self = ncm_stats_dist_kde_get_instance_private (sdkde);

  return self->cov_fixed;
}


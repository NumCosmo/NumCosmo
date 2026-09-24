/***************************************************************************
 *            ncm_stats_dist_vkde.c
 *
 *  Wed November 07 16:02:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_vkde.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmStatsDistVKDE
 *
 * Variable-bandwidth kernel estimator for #NcmStatsDist.
 *
 * Uses one covariance matrix per sample point, computed in
 * ncm_stats_dist_prepare_shapes() from the @local_frac nearest sample points.
 * The rest of the calculation follows #NcmStatsDist and #NcmStatsDistKDE, with
 * a different covariance matrix and normalization factor per kernel.
 *
 * With #NcmStatsDist:center-shrink enabled the shrinkage scale $s^2$ is the mean
 * of $\mathrm{tr}(C_i \Sigma^{-1}) / d$ over the local covariances $C_i$, so that the
 * mixture covariance matches the sample covariance $\Sigma$; the local covariances
 * themselves are still estimated around the original sample points.
 *
 * The caller must supply @sdk and @CV_type through ncm_stats_dist_vkde_new(),
 * @y through ncm_stats_dist_add_obs(), @split_frac through
 * ncm_stats_dist_set_split_frac(), @over_smooth through
 * ncm_stats_dist_set_over_smooth(), @local_frac through
 * ncm_stats_dist_vkde_set_local_frac(), and $v(x)$ through
 * ncm_stats_dist_prepare().
 *
 * The flowchart below gives the call order.
 *
 * ![vkde_sketch](vkde.png)
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist_vkde.h"
#include "ncm/core/ncm_iset.h"
#include "ncm/algebra/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "external/misc/kdtree.h"
#include "external/misc/rb_knn_list.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "ncm/stats/ncm_stats_dist_vkde_private.h"
#include "ncm/stats/ncm_stats_dist_kde_private.h"
#include "ncm/stats/ncm_stats_dist_private.h"

enum
{
  PROP_0,
  PROP_LOCAL_FRAC,
  PROP_USE_ROT_HREF,
  PROP_POINTS_PER_DIM,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistVKDE, ncm_stats_dist_vkde, NCM_TYPE_STATS_DIST_KDE)

static NcmStatsDistPrivate *
ncm_stats_dist_get_instance_private (NcmStatsDist * sd)
{
  return g_type_instance_get_private ((GTypeInstance *) sd, NCM_TYPE_STATS_DIST);
}

static NcmStatsDistKDEPrivate *
ncm_stats_dist_kde_get_instance_private (NcmStatsDistKDE *sd)
{
  return g_type_instance_get_private ((GTypeInstance *) sd, NCM_TYPE_STATS_DIST_KDE);
}

static gpointer
_ncm_stats_dist_vkde_stats_vec_new (gpointer userdata)
{
  NcmStatsDist *sd                   = NCM_STATS_DIST (userdata);
  NcmStatsDistPrivate * const ppself = ncm_stats_dist_get_instance_private (sd);
  NcmStatsVec *sample                = ncm_stats_vec_new (ppself->d, NCM_STATS_VEC_COV, TRUE);

  return sample;
}

typedef struct _NcmStatsDistVKDEEvalVars
{
  NcmVector *delta_x;
  NcmVector *chi2;
  NcmVector *lnK;
  NcmMatrix *chi2_M;
  NcmVector *chi2_row;
  NcmMatrix *X;
  NcmMatrix *delta_X;
  NcmVector *delta_q;
} NcmStatsDistVKDEEvalVars;

static gpointer
_ncm_stats_dist_vkde_eval_vars_new (gpointer userdata)
{
  NcmStatsDist *sd                   = NCM_STATS_DIST (userdata);
  NcmStatsDistPrivate * const ppself = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistVKDEEvalVars *ev       = g_new0 (NcmStatsDistVKDEEvalVars, 1);


  ev->delta_x = ncm_vector_new (ppself->d);
  ev->chi2    = ncm_vector_new (ppself->n_kernels);
  ev->lnK     = ncm_vector_new (ppself->n_kernels);

  /* Only the batched evaluator needs the tile, so it is built on first use: a thread
   * that never takes that path does not carry the n_kernels columns. */
  ev->chi2_M   = NULL;
  ev->chi2_row = NULL;
  ev->X        = NULL;
  ev->delta_X  = NULL;
  ev->delta_q  = NULL;

  return ev;
}

static void
_ncm_stats_dist_vkde_eval_vars_free (gpointer userdata)
{
  NcmStatsDistVKDEEvalVars *ev = (NcmStatsDistVKDEEvalVars *) userdata;

  ncm_vector_free (ev->delta_x);
  ncm_vector_free (ev->chi2);
  ncm_vector_free (ev->lnK);
  ncm_matrix_clear (&ev->chi2_M);
  ncm_vector_clear (&ev->chi2_row);
  ncm_matrix_clear (&ev->X);
  ncm_matrix_clear (&ev->delta_X);
  ncm_vector_clear (&ev->delta_q);

  g_free (ev);
}

static void
ncm_stats_dist_vkde_init (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  self->cov_array    = g_ptr_array_new ();
  self->cov_array0   = g_ptr_array_new ();
  self->lnnorms      = NULL;
  self->IM_delta     = NULL;
  self->IM_delta_row = NULL;

  self->local_frac     = 0.0;
  self->points_per_dim = 0.0;
  self->use_rot_href   = FALSE;

  self->mp_stats_vec = ncm_memory_pool_new (&_ncm_stats_dist_vkde_stats_vec_new, sdvkde,
                                            (GDestroyNotify) & ncm_stats_vec_free);

  self->mp_eval_vars = ncm_memory_pool_new (&_ncm_stats_dist_vkde_eval_vars_new, sdvkde,
                                            &_ncm_stats_dist_vkde_eval_vars_free);

  g_ptr_array_set_free_func (self->cov_array, (GDestroyNotify) ncm_matrix_free);
  g_ptr_array_set_free_func (self->cov_array0, (GDestroyNotify) ncm_matrix_free);
}

static void
_ncm_stats_dist_vkde_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistVKDE *sdvkde = NCM_STATS_DIST_VKDE (object);

  /*g_return_if_fail (NCM_IS_STATS_DIST (object));*/

  switch (prop_id)
  {
    case PROP_LOCAL_FRAC:
      ncm_stats_dist_vkde_set_local_frac (sdvkde, g_value_get_double (value));
      break;
    case PROP_USE_ROT_HREF:
      ncm_stats_dist_vkde_set_use_rot_href (sdvkde, g_value_get_boolean (value));
      break;
    case PROP_POINTS_PER_DIM:
      ncm_stats_dist_vkde_set_points_per_dim (sdvkde, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_vkde_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistVKDE *sdvkde = NCM_STATS_DIST_VKDE (object);

  /*NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);*/

  g_return_if_fail (NCM_IS_STATS_DIST_VKDE (object));

  switch (prop_id)
  {
    case PROP_LOCAL_FRAC:
      g_value_set_double (value, ncm_stats_dist_vkde_get_local_frac (sdvkde));
      break;
    case PROP_USE_ROT_HREF:
      g_value_set_boolean (value, ncm_stats_dist_vkde_get_use_rot_href (sdvkde));
      break;
    case PROP_POINTS_PER_DIM:
      g_value_set_double (value, ncm_stats_dist_vkde_get_points_per_dim (sdvkde));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_vkde_dispose (GObject *object)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (object);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  ncm_vector_clear (&self->lnnorms);
  ncm_matrix_clear (&self->IM_delta);
  ncm_vector_clear (&self->IM_delta_row);

  g_clear_pointer (&self->cov_array, g_ptr_array_unref);
  g_clear_pointer (&self->cov_array0, g_ptr_array_unref);

  if (self->mp_stats_vec != NULL)
  {
    ncm_memory_pool_free (self->mp_stats_vec, TRUE);
    self->mp_stats_vec = NULL;
  }

  if (self->mp_eval_vars != NULL)
  {
    ncm_memory_pool_free (self->mp_eval_vars, TRUE);
    self->mp_eval_vars = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_vkde_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_vkde_finalize (GObject *object)
{
  /* NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (object); */
  /* NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde); */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_vkde_parent_class)->finalize (object);
}

static void _ncm_stats_dist_vkde_set_dim (NcmStatsDist *sd, const guint dim);
static gdouble _ncm_stats_dist_vkde_bandwidth (NcmStatsDist *sd);
static void _ncm_stats_dist_vkde_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array);
static void _ncm_stats_dist_vkde_prepare_kernels (NcmStatsDist *sd);
static void _ncm_stats_dist_vkde_compute_IM (NcmStatsDist *sd, NcmMatrix *IM);
static NcmMatrix *_ncm_stats_dist_vkde_peek_cov_decomp (NcmStatsDist *sd, guint i);
static gdouble _ncm_stats_dist_vkde_get_lnnorm (NcmStatsDist *sd, guint i);
static gdouble _ncm_stats_dist_vkde_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
static gdouble _ncm_stats_dist_vkde_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
static void _ncm_stats_dist_vkde_eval_weights_m2lnp_vec (NcmStatsDist *sd, NcmVector *weights, GPtrArray *x_a, NcmVector *m2lnp);
static void _ncm_stats_dist_vkde_reset (NcmStatsDist *sd);

static void
ncm_stats_dist_vkde_class_init (NcmStatsDistVKDEClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcmStatsDistClass *base_class = NCM_STATS_DIST_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_vkde_set_property;
  object_class->get_property = &_ncm_stats_dist_vkde_get_property;
  object_class->dispose      = &_ncm_stats_dist_vkde_dispose;
  object_class->finalize     = &_ncm_stats_dist_vkde_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LOCAL_FRAC,
                                   g_param_spec_double ("local-frac",
                                                        NULL,
                                                        "Fraction to use in the local kernel covariance computation",
                                                        0.001, 1.0, 0.05,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_USE_ROT_HREF,
                                   g_param_spec_boolean ("use-rot-href",
                                                         NULL,
                                                         "Whether to use the href rule-of-thumb to compute the final bandwidth",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDistVKDE:points-per-dim:
   *
   * Number of nearest neighbors per dimension used for each local scale matrix,
   * $k = \min(n, \lceil c\, d \rceil)$ with $c$ this value, $d$ the dimension and $n$ the
   * sample size. It replaces #NcmStatsDistVKDE:local-frac when positive: the neighbor
   * count follows what a $d \times d$ covariance estimate needs, so a small sample turns
   * every local covariance into the global one (the KDE limit) while a large sample keeps
   * the kernels local. Zero (the default) keeps the fraction of the sample.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_POINTS_PER_DIM,
                                   g_param_spec_double ("points-per-dim",
                                                        NULL,
                                                        "Nearest neighbors per dimension for the local covariances (0: use local-frac)",
                                                        0.0, 1.0e6, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  base_class->set_dim                = &_ncm_stats_dist_vkde_set_dim;
  base_class->bandwidth              = &_ncm_stats_dist_vkde_bandwidth;
  base_class->prepare_shapes         = &_ncm_stats_dist_vkde_prepare_shapes;
  base_class->prepare_kernels        = &_ncm_stats_dist_vkde_prepare_kernels;
  base_class->compute_IM             = &_ncm_stats_dist_vkde_compute_IM;
  base_class->peek_cov_decomp        = &_ncm_stats_dist_vkde_peek_cov_decomp;
  base_class->get_lnnorm             = &_ncm_stats_dist_vkde_get_lnnorm;
  base_class->eval_weights           = &_ncm_stats_dist_vkde_eval_weights;
  base_class->eval_weights_m2lnp     = &_ncm_stats_dist_vkde_eval_weights_m2lnp;
  base_class->eval_weights_m2lnp_vec = &_ncm_stats_dist_vkde_eval_weights_m2lnp_vec;
  base_class->reset                  = &_ncm_stats_dist_vkde_reset;
}

static void
_ncm_stats_dist_vkde_set_dim (NcmStatsDist *sd, const guint dim)
{
  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->set_dim  (sd, dim);
  {
    NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
    NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

    g_ptr_array_set_size (self->cov_array, 0);
    g_ptr_array_set_size (self->cov_array0, 0);
    ncm_matrix_clear (&self->IM_delta);
    ncm_vector_clear (&self->IM_delta_row);
  }
}

static gdouble
_ncm_stats_dist_vkde_bandwidth (NcmStatsDist *sd)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);

  if (self->use_rot_href)
  {
    /* Chain up : start */
    const gdouble href_base = NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->bandwidth (sd);

    return href_base * ppself->n_obs / (1.0 * ncm_stats_dist_vkde_get_n_neighbors (sdvkde, ppself->n_obs));
  }
  else
  {
    return ppself->over_smooth;
  }
}

static void
_ncm_stats_dist_vkde_build_cov_array_kdtree (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistKDEPrivate * const pself = ncm_stats_dist_kde_get_instance_private (NCM_STATS_DIST_KDE (sd));
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);

  /*
   * Creates a near tree object and add all transformed vectors.
   */
  struct kdtree *tree = kdtree_init (ppself->d);
  guint i;

  g_assert_cmpint (ppself->n_obs, >, 2);

  for (i = 0; i < ppself->n_obs; i++)
  {
    NcmVector *invUtheta_i = g_ptr_array_index (pself->invUsample_array, i);

    /*
     * Inserting the transformed vector in the tree, saving also the index.
     */
    kdtree_insert (tree, ncm_vector_data (invUtheta_i));
  }

  kdtree_rebuild (tree);

  /*
   * Checking allocation of the norm vector and
   * covariance array.
   */
  if ((self->lnnorms == NULL) || (ncm_vector_len (self->lnnorms) != ppself->n_kernels))
  {
    ncm_vector_clear (&self->lnnorms);
    self->lnnorms = ncm_vector_new (ppself->n_kernels);

    ncm_memory_pool_empty (self->mp_eval_vars, TRUE);
  }

  /* Scratch for compute_IM: the sample, centred and whitened one kernel at a time. */
  if ((self->IM_delta == NULL) || (ncm_matrix_nrows (self->IM_delta) != ppself->n_obs))
  {
    ncm_matrix_clear (&self->IM_delta);
    ncm_vector_clear (&self->IM_delta_row);
    self->IM_delta     = ncm_matrix_new (ppself->n_obs, ppself->d);
    self->IM_delta_row = ncm_vector_new_data_static (ncm_matrix_ptr (self->IM_delta, 0, 0), ppself->d, 1);
  }

  /*
   * Checking allocation of the covariance array.
   */
  {
    guint cur_size = self->cov_array->len;

    g_ptr_array_set_size (self->cov_array, ppself->n_kernels);

    if (cur_size < ppself->n_kernels)
    {
      for (i = cur_size; i < ppself->n_kernels; i++)
      {
        g_ptr_array_index (self->cov_array, i) = ncm_matrix_new (ppself->d, ppself->d);
      }
    }
  }

  /*
   * Lets find the k nearest neighbors of each vector in the
   * sample and use them to define the local covariance at each
   * vector location.
   */
  {
    const size_t k = ncm_stats_dist_vkde_get_n_neighbors (sdvkde, ppself->n_obs);

    /* #pragma omp parallel for schedule(static) if (ppself->use_threads) */

    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmVector *invUtheta_i   = g_ptr_array_index (pself->invUsample_array, i);
      NcmStatsVec **sample_ptr = ncm_memory_pool_get (self->mp_stats_vec);
      NcmStatsVec *sample      = *sample_ptr;
      rb_knn_list_table_t *table;

      /*gint tid = omp_get_thread_num(); */
      /*printf("Hello world from omp thread %d\n", tid); */

      table = kdtree_knn_search (tree, ncm_vector_data (invUtheta_i), k);
      {
        rb_knn_list_traverser_t trav;
        knn_list_t *p;

        p = rb_knn_list_t_first (&trav, table);

        do {
          NcmVector *ni = g_ptr_array_index (sample_array, p->node->coord_index);

          ncm_stats_vec_append (sample, ni, FALSE);
        } while ((p = rb_knn_list_t_next (&trav)) != NULL);
      }
      rb_knn_list_destroy (table);

      /*
       * Saving the covariance for each vector.
       */
      {
        NcmMatrix *sample_cov = NULL;

        switch (pself->cov_type)
        {
          case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
          case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:

            sample_cov = ncm_matrix_ref (ncm_stats_vec_peek_cov_matrix (sample, 0));
            break;
          case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
            sample_cov = ncm_stats_vec_compute_cov_robust_diag (sample); /* */
            break;
          case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
            sample_cov = ncm_stats_vec_compute_cov_robust_ogk (sample); /* */
            break;
          default:
            g_assert_not_reached ();
            break;
        }

        {
          NcmMatrix *cov_decomp = g_ptr_array_index (self->cov_array, i);

          _ncm_stats_dist_cholesky (cov_decomp, sample_cov, pself->nearPD_maxiter, "a local kernel covariance");
          ncm_matrix_free (sample_cov);
        }
      }

      ncm_stats_vec_reset (sample, TRUE);
      ncm_memory_pool_return (sample_ptr);
    }

    /*
     * What center shrinkage needs from a variable-bandwidth estimator: the mean of the
     * kernel scale matrices. The factors are copied before the transform is applied,
     * since the applied ones follow the bandwidth.
     */
    {
      NcmMatrix *mean_cov = ppself->kernel_cov;
      guint j;

      {
        const guint cur_size = self->cov_array0->len;

        g_ptr_array_set_size (self->cov_array0, ppself->n_kernels);

        for (j = cur_size; j < ppself->n_kernels; j++)
          g_ptr_array_index (self->cov_array0, j) = ncm_matrix_new (ppself->d, ppself->d);
      }

      ncm_matrix_set_zero (mean_cov);

      for (j = 0; j < ppself->n_kernels; j++)
      {
        NcmMatrix *U_j = g_ptr_array_index (self->cov_array, j);

        /* The decomposition leaves the covariance's lower triangle under the factor;
         * cleared once so both copies are plain upper factors and dsyrk reads a clean
         * matrix. Only the diagonal and the upper triangle are read anywhere else. */
        ncm_matrix_zero_triangle (U_j, 'U');
        ncm_matrix_memcpy (g_ptr_array_index (self->cov_array0, j), U_j);
        ncm_matrix_dsyrk (mean_cov, 'U', 'T', 1.0 / (1.0 * ppself->n_kernels), U_j, 1.0);
      }

      ncm_matrix_copy_triangle (mean_cov, 'U');
    }
  }
  kdtree_destroy (tree);
}

static void
_ncm_stats_dist_vkde_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);

  {
    const guint k = ncm_stats_dist_vkde_get_n_neighbors (sdvkde, ppself->n_obs);

    /* A covariance from k <= d points is singular by construction. */
    if (k <= ppself->d)
      g_error ("_ncm_stats_dist_vkde_prepare_shapes: %u neighbors cannot define a local covariance in "
               "%u dimensions; raise local-frac (%g) or points-per-dim (%g).",
               k, ppself->d, self->local_frac, self->points_per_dim);
  }

  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->prepare_shapes (sd, sample_array);
  _ncm_stats_dist_vkde_build_cov_array_kdtree (sd, sample_array);
}

static void
_ncm_stats_dist_vkde_prepare_kernels (NcmStatsDist *sd)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistKernel *kernel           = ncm_stats_dist_peek_kernel (sd);
  guint i;

  /*
   * The applied factors follow the current center transform, from the untransformed
   * copies kept by prepare_shapes(); the normalizations follow the factors and the
   * kernel, whose nu the fit may have changed. The centers themselves are read from
   * NcmStatsDist's own center array.
   */
  g_assert_cmpuint (self->cov_array0->len, ==, ppself->n_kernels);

  for (i = 0; i < ppself->n_kernels; i++)
  {
    NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);

    _ncm_stats_dist_refactor_decomp (sd, g_ptr_array_index (self->cov_array0, i), cov_decomp_i);
    ncm_vector_set (self->lnnorms, i, ncm_stats_dist_kernel_get_lnnorm (kernel, cov_decomp_i));
  }
}

static void
_ncm_stats_dist_vkde_compute_IM (NcmStatsDist *sd, NcmMatrix *IM)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistKDEPrivate * const pself = ncm_stats_dist_kde_get_instance_private (NCM_STATS_DIST_KDE (sd));
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble href2                  = ppself->href * ppself->href;
  const gdouble one_href2              = 1.0 / href2;

  guint i;

  /* #pragma omp parallel if (ppself->use_threads) */
  {
    NcmMatrix *invUsample_matrix = self->IM_delta;
    NcmVector *theta_j           = self->IM_delta_row;

    /* #pragma omp for schedule(static) */
    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
      NcmVector *theta_i      = g_ptr_array_index (ppself->center_array, i);
      guint j;

      ncm_matrix_memcpy (invUsample_matrix, pself->sample_matrix);
      ncm_matrix_sub_row_vector (invUsample_matrix, theta_i);
      ncm_matrix_dtrsm (invUsample_matrix, 'R', 'U', 'N', 1.0, cov_decomp_i);

      for (j = 0; j < ppself->n_obs; j++)
      {
        ncm_vector_replace_data (theta_j, ncm_matrix_ptr (invUsample_matrix, j, 0));
        ncm_matrix_set (IM, j, i, ncm_vector_dot (theta_j, theta_j) * one_href2);
      }
    }
  }

  {
    const gdouble lnnorm_href = ppself->d * log (ppself->href);

    /* #pragma omp parallel for if (ppself->use_threads) */

    for (i = 0; i < ppself->n_obs; i++)
    {
      NcmVector *row_i = ncm_matrix_get_row (IM, i);

      ncm_stats_dist_kernel_eval_unnorm_vec (ppself->kernel, row_i, row_i);
      ncm_vector_free (row_i);
    }

    /* #pragma omp parallel for if (ppself->use_threads) */

    for (i = 0; i < ppself->n_kernels; i++)
    {
      const gdouble norm_i = exp (ncm_vector_fast_get (self->lnnorms, i) + lnnorm_href);

      ncm_matrix_mul_col (IM, i, 1.0 / norm_i);
    }
  }
}

static NcmMatrix *
_ncm_stats_dist_vkde_peek_cov_decomp (NcmStatsDist *sd, guint i)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  g_assert (i < self->cov_array->len);

  return g_ptr_array_index (self->cov_array, i);
}

static gdouble
_ncm_stats_dist_vkde_get_lnnorm (NcmStatsDist *sd, guint i)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);

  g_assert (i < self->cov_array->len);

  return ncm_vector_fast_get (self->lnnorms, i) + ppself->d * log (ppself->href);
}

static gdouble
_ncm_stats_dist_vkde_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble href2                  = ppself->href * ppself->href;
  const gdouble one_href2              = 1.0 / href2;
  NcmStatsDistVKDEEvalVars **ev_ptr    = ncm_memory_pool_get (self->mp_eval_vars);
  NcmStatsDistVKDEEvalVars *ev         = *ev_ptr;
  gdouble s                            = 0.0;
  guint i;

  {
    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
      NcmVector *theta_i      = g_ptr_array_index (ppself->center_array, i);

      ncm_vector_memcpy (ev->delta_x, x);
      ncm_vector_axpy (ev->delta_x, -1.0, theta_i);

      ncm_matrix_dtrsv (cov_decomp_i, 'U', 'T', ev->delta_x);

      {
        const gdouble chi2_i = ncm_vector_dot (ev->delta_x, ev->delta_x) * one_href2;

        ncm_vector_fast_set (ev->chi2, i, chi2_i);
      }
    }
  }

  ncm_stats_dist_kernel_eval_unnorm_vec (ppself->kernel, ev->chi2, ev->chi2);

  for (i = 0; i < ppself->n_kernels; i++)
  {
    const gdouble Ku_i = ncm_vector_fast_get (ev->chi2, i);
    const gdouble u_i  = exp (ncm_vector_fast_get (self->lnnorms, i));
    const gdouble w_i  = ncm_vector_fast_get (weights, i);

    s += w_i * (Ku_i / u_i);
  }

  ncm_memory_pool_return (ev_ptr);

  return s / pow (ppself->href, ppself->d);
}

static gdouble
_ncm_stats_dist_vkde_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble href2                  = ppself->href * ppself->href;
  const gdouble one_href2              = 1.0 / href2;
  NcmStatsDistVKDEEvalVars **ev_ptr    = ncm_memory_pool_get (self->mp_eval_vars);
  NcmStatsDistVKDEEvalVars *ev         = *ev_ptr;
  guint i;

  {
    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
      NcmVector *theta_i      = g_ptr_array_index (ppself->center_array, i);

      ncm_vector_memcpy (ev->delta_x, x);
      ncm_vector_axpy (ev->delta_x, -1.0, theta_i);

      ncm_matrix_dtrsv (cov_decomp_i, 'U', 'T', ev->delta_x);

      {
        const gdouble chi2_i = ncm_vector_dot (ev->delta_x, ev->delta_x) * one_href2;

        ncm_vector_fast_set (ev->chi2, i, chi2_i);
      }
    }
  }

  {
    gdouble gamma, lambda;

    ncm_stats_dist_kernel_eval_sum0_gamma_lambda (ppself->kernel, ev->chi2, weights, self->lnnorms, ev->lnK, &gamma, &lambda);

    ncm_memory_pool_return (ev_ptr);

    return -2.0 * (gamma + log1p (lambda) - ppself->d * log (ppself->href));
  }
}

/*
 * Largest number of points per sweep over the kernels. Each covariance factor is then
 * read once per tile instead of once per point, which is where the gain comes from: on a
 * 2600-kernel, d = 50 ensemble the speedup saturates by 256 points (4.2x against one
 * point at a time), and the chi2 tile is 5 MB rather than the 103 MB a full 5200-point
 * batch would need. The batch is cut into equal tiles no larger than this, so that no
 * tile is left nearly empty.
 */
#define _NCM_STATS_DIST_VKDE_EVAL_TILE (256)

static void
_ncm_stats_dist_vkde_eval_weights_m2lnp_vec (NcmStatsDist *sd, NcmVector *weights, GPtrArray *x_a, NcmVector *m2lnp)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);
  const gdouble one_href2              = 1.0 / (ppself->href * ppself->href);
  const gdouble lnnorm_href            = ppself->d * log (ppself->href);
  const guint np                       = x_a->len;
  const guint n_tiles                  = MAX (1u, (np + _NCM_STATS_DIST_VKDE_EVAL_TILE - 1) / _NCM_STATS_DIST_VKDE_EVAL_TILE);
  const guint nt                       = (np + n_tiles - 1) / n_tiles;

  NcmStatsDistVKDEEvalVars **ev_ptr = ncm_memory_pool_get (self->mp_eval_vars);
  NcmStatsDistVKDEEvalVars *ev      = *ev_ptr;
  guint p0;

  if ((ev->chi2_M == NULL) || (ncm_matrix_col_len (ev->chi2_M) < nt))
  {
    ncm_matrix_clear (&ev->chi2_M);
    ncm_vector_clear (&ev->chi2_row);

    ncm_matrix_clear (&ev->X);
    ncm_matrix_clear (&ev->delta_X);
    ncm_vector_clear (&ev->delta_q);
    ev->chi2_M   = ncm_matrix_new (nt, ppself->n_kernels);
    ev->chi2_row = ncm_vector_new_data_static (ncm_matrix_ptr (ev->chi2_M, 0, 0), ppself->n_kernels, 1);
    ev->X        = ncm_matrix_new (nt, ppself->d);
    ev->delta_X  = ncm_matrix_new (nt, ppself->d);
    ev->delta_q  = ncm_vector_new_data_static (ncm_matrix_ptr (ev->delta_X, 0, 0), ppself->d, 1);
  }

  for (p0 = 0; p0 < np; p0 += nt)
  {
    const guint ntp = MIN (nt, np - p0);
    guint p;

    /* The points of the tile, once; every kernel subtracts its own centre from them. The
     * tile has nt rows; the last one is shorter by less than n_tiles rows, and the solve
     * on its unused trailing rows is harmless. */
    for (p = 0; p < ntp; p++)
      ncm_matrix_set_row (ev->X, p, g_ptr_array_index (x_a, p0 + p));

    /* #pragma omp parallel if (ppself->use_threads) */
    {
      NcmMatrix *delta_X = ev->delta_X;
      NcmVector *delta_q = ev->delta_q;
      guint i;

      /* #pragma omp for schedule(static) */
      for (i = 0; i < ppself->n_kernels; i++)
      {
        NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
        NcmVector *theta_i      = g_ptr_array_index (ppself->center_array, i);
        guint q;

        ncm_matrix_memcpy (delta_X, ev->X);
        ncm_matrix_sub_row_vector (delta_X, theta_i);
        ncm_matrix_dtrsm (delta_X, 'R', 'U', 'N', 1.0, cov_decomp_i);

        for (q = 0; q < ntp; q++)
        {
          ncm_vector_replace_data (delta_q, ncm_matrix_ptr (delta_X, q, 0));
          ncm_matrix_set (ev->chi2_M, q, i, ncm_vector_dot (delta_q, delta_q) * one_href2);
        }
      }
    }

    for (p = 0; p < ntp; p++)
    {
      gdouble gamma, lambda;

      ncm_vector_replace_data (ev->chi2_row, ncm_matrix_ptr (ev->chi2_M, p, 0));

      ncm_stats_dist_kernel_eval_sum0_gamma_lambda (ppself->kernel, ev->chi2_row, weights, self->lnnorms, ev->lnK, &gamma, &lambda);
      ncm_vector_set (m2lnp, p0 + p, -2.0 * (gamma + log1p (lambda) - lnnorm_href));
    }
  }

  ncm_memory_pool_return (ev_ptr);
}

static void
_ncm_stats_dist_vkde_reset (NcmStatsDist *sd)
{
  /* Chain up : end */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->reset (sd);
}

/**
 * ncm_stats_dist_vkde_new:
 * @sdk: a #NcmStatsDistKernel
 * @CV_type: a #NcmStatsDistCV
 *
 * Creates a new #NcmStatsDistVKDE object using @sdk as
 * kernel and @CV_type as cross-validation method.
 *
 * Returns: (transfer full): the newly created #NcmStatsDistVKDE object.
 */
NcmStatsDistVKDE *
ncm_stats_dist_vkde_new (NcmStatsDistKernel *sdk, NcmStatsDistCV CV_type)
{
  NcmStatsDistVKDE *sdvkde = g_object_new (NCM_TYPE_STATS_DIST_VKDE,
                                           "kernel", sdk,
                                           "CV-type", CV_type,
                                           NULL);

  return sdvkde;
}

/**
 * ncm_stats_dist_vkde_ref:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Increases the reference count of @sdvkde.
 *
 * Returns: (transfer full): @sdvkde.
 */
NcmStatsDistVKDE *
ncm_stats_dist_vkde_ref (NcmStatsDistVKDE *sdvkde)
{
  return g_object_ref (sdvkde);
}

/**
 * ncm_stats_dist_vkde_free:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Decreases the reference count of @sdvkde.
 *
 */
void
ncm_stats_dist_vkde_free (NcmStatsDistVKDE *sdvkde)
{
  g_object_unref (sdvkde);
}

/**
 * ncm_stats_dist_vkde_clear:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Decreases the reference count of *@sdvkde and sets the pointer *@sdvkde to NULL.
 *
 */
void
ncm_stats_dist_vkde_clear (NcmStatsDistVKDE **sdvkde)
{
  g_clear_object (sdvkde);
}

/**
 * ncm_stats_dist_vkde_set_local_frac:
 * @sdvkde: a #NcmStatsDistVKDE
 * @local_frac: the over-smooth factor
 *
 * Sets local kernel fraction to @local_frac. This fraction
 * defines the amount of closest points from each sample point
 * that will be used to compute the covariance matrix of each point.
 *
 */
void
ncm_stats_dist_vkde_set_local_frac (NcmStatsDistVKDE *sdvkde, const gdouble local_frac)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  g_assert_cmpfloat (local_frac, >=, 0.001);
  g_assert_cmpfloat (local_frac, <=, 1.0);

  self->local_frac = local_frac;
}

/**
 * ncm_stats_dist_vkde_get_local_frac:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Returns: a double @local_frac, the local kernel fraction.
 */
gdouble
ncm_stats_dist_vkde_get_local_frac (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  return self->local_frac;
}

/**
 * ncm_stats_dist_vkde_set_use_rot_href:
 * @sdvkde: a #NcmStatsDistVKDE
 * @use_rot_href: whether to use the rule of thumb bandwidth
 *
 * Sets whether to use the rule of thumb bandwidth for the
 *
 */
void
ncm_stats_dist_vkde_set_use_rot_href (NcmStatsDistVKDE *sdvkde, const gboolean use_rot_href)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  self->use_rot_href = use_rot_href;
}

/**
 * ncm_stats_dist_vkde_get_use_rot_href:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Returns: whether to use the rule of thumb bandwidth.
 */
gboolean
ncm_stats_dist_vkde_get_use_rot_href (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  return self->use_rot_href;
}

/**
 * ncm_stats_dist_vkde_set_points_per_dim:
 * @sdvkde: a #NcmStatsDistVKDE
 * @points_per_dim: nearest neighbors per dimension, zero to use #NcmStatsDistVKDE:local-frac
 *
 * Sets #NcmStatsDistVKDE:points-per-dim. Takes effect at the next preparation.
 *
 */
void
ncm_stats_dist_vkde_set_points_per_dim (NcmStatsDistVKDE *sdvkde, const gdouble points_per_dim)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  g_assert_cmpfloat (points_per_dim, >=, 0.0);
  self->points_per_dim = points_per_dim;
}

/**
 * ncm_stats_dist_vkde_get_points_per_dim:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Returns: #NcmStatsDistVKDE:points-per-dim.
 */
gdouble
ncm_stats_dist_vkde_get_points_per_dim (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  return self->points_per_dim;
}

/**
 * ncm_stats_dist_vkde_get_n_neighbors:
 * @sdvkde: a #NcmStatsDistVKDE
 * @n_obs: sample size
 *
 * Number of nearest neighbors used for each local scale matrix with a sample of @n_obs
 * points: $\min(n_\mathrm{obs}, \lceil c\, d \rceil)$ when #NcmStatsDistVKDE:points-per-dim
 * $c$ is positive, otherwise #NcmStatsDistVKDE:local-frac times @n_obs; at least 2.
 *
 * Returns: the neighbor count.
 */
guint
ncm_stats_dist_vkde_get_n_neighbors (NcmStatsDistVKDE *sdvkde, const guint n_obs)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  const guint d                        = ncm_stats_dist_get_dim (NCM_STATS_DIST (sdvkde));
  gdouble k;

  if (self->points_per_dim > 0.0)
    k = GSL_MIN (1.0 * n_obs, ceil (self->points_per_dim * d));
  else
    k = self->local_frac * n_obs;

  return (guint) GSL_MAX (k, 2.0);
}


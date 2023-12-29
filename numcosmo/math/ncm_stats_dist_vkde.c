/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

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
 * SECTION:ncm_stats_dist_vkde
 * @title: NcmStatsDistVKDE
 * @short_description: Abstract class for implementing N-dimensional probability distributions with a variable density estimator kernel.
 *
 * Abstract object to reconstruct an arbitrary N-dimensional probability distribution.
 * This object provides the complementary tools to perform a radial basis interpolation
 * in a multidimensional function using the #NcmStatsDist class.
 *
 * This object sets the kernel $\phi$ to be used in the radial basis interpolation. This object also implements some
 * calculations needed in the #NcmStatsDist class, such as the covariance matrices of the whole sample points and its Cholesky decompositions,
 * the preparation of the interpolation matrix $IM$, the kernel normalization factors, and given a sample vector $\vec{x}$, the distribution
 * evaluated in these points. Some of these calculations are explained below.
 *
 * The #NcmStatsDistVKDE uses a different covariance matrix for each sample point. This feature is computed
 * in the ncm_stats_dist_prepare_kernel() function. In this algorithm, one should define the @local_frac parameter, that is,
 * the fraction of nearest sample points that will be used to compute each covariance matrix of each
 * sample point. This is done by calling the function ncm_stats_dist_vkde_set_local_frac().
 * The rest of the calculation follows the same procedure as the #NcmStatsDist and #NcmStatsDistKDE objects,
 * using now a different covariance matrix and normalization factor for each kernel. For more information about
 * how the #NcmStatsDist class works, check #NcmStatsDist and #NcmStatsDistKDE objects.
 *
 * The user must provide input the values: @sdk, @CV_type - ncm_stats_dist_vkde_new(), @y - ncm_stats_dist_add_obs(), @split_frac - ncm_stats_dist_set_split_frac(),
 * @over_smooth - ncm_stats_dist_set_over_smooth(), @local_Frac - ncm_stats_dist_vkde_set_local_frac(), $v(x)$ - ncm_stats_dist_prepare_interp().
 * To see an example of how to use this object and the main functions that are called within each function, check the fluxogram at the end of this documentation,
 * where the order of the functions that should be called by the user and some of the functions that the algorithm calls.
 *
 * ![vkde_sketch](vkde.png)
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_vkde.h"
#include "math/ncm_iset.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <omp.h>
#include "misc/kdtree.h"
#include "misc/rb_knn_list.h"
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "math/ncm_stats_dist_vkde_private.h"
#include "math/ncm_stats_dist_kde_private.h"
#include "math/ncm_stats_dist_private.h"

enum
{
  PROP_0,
  PROP_LOCAL_FRAC,
  PROP_USE_ROT_HREF,
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

  return ev;
}

static void
_ncm_stats_dist_vkde_eval_vars_free (gpointer userdata)
{
  NcmStatsDistVKDEEvalVars *ev = (NcmStatsDistVKDEEvalVars *) userdata;

  ncm_vector_free (ev->delta_x);
  ncm_vector_free (ev->chi2);
  ncm_vector_free (ev->lnK);

  g_free (ev);
}

static void
ncm_stats_dist_vkde_init (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  self->cov_array = g_ptr_array_new ();
  self->lnnorms   = NULL;

  self->local_frac   = 0.0;
  self->use_rot_href = FALSE;

  self->mp_stats_vec = ncm_memory_pool_new (&_ncm_stats_dist_vkde_stats_vec_new, sdvkde,
                                            (GDestroyNotify) & ncm_stats_vec_free);

  self->mp_eval_vars = ncm_memory_pool_new (&_ncm_stats_dist_vkde_eval_vars_new, sdvkde,
                                            &_ncm_stats_dist_vkde_eval_vars_free);

  g_ptr_array_set_free_func (self->cov_array, (GDestroyNotify) ncm_matrix_free);
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

  g_clear_pointer (&self->cov_array, g_ptr_array_unref);

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
static gdouble _ncm_stats_dist_vkde_get_href (NcmStatsDist *sd);
static void _ncm_stats_dist_vkde_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array);
static void _ncm_stats_dist_vkde_compute_IM (NcmStatsDist *sd, NcmMatrix *IM);
static NcmMatrix *_ncm_stats_dist_vkde_peek_cov_decomp (NcmStatsDist *sd, guint i);
static gdouble _ncm_stats_dist_vkde_get_lnnorm (NcmStatsDist *sd, guint i);
static gdouble _ncm_stats_dist_vkde_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
static gdouble _ncm_stats_dist_vkde_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
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

  base_class->set_dim            = &_ncm_stats_dist_vkde_set_dim;
  base_class->get_href           = &_ncm_stats_dist_vkde_get_href;
  base_class->prepare_kernel     = &_ncm_stats_dist_vkde_prepare_kernel;
  base_class->compute_IM         = &_ncm_stats_dist_vkde_compute_IM;
  base_class->peek_cov_decomp    = &_ncm_stats_dist_vkde_peek_cov_decomp;
  base_class->get_lnnorm         = &_ncm_stats_dist_vkde_get_lnnorm;
  base_class->eval_weights       = &_ncm_stats_dist_vkde_eval_weights;
  base_class->eval_weights_m2lnp = &_ncm_stats_dist_vkde_eval_weights_m2lnp;
  base_class->reset              = &_ncm_stats_dist_vkde_reset;
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
  }
}

static gdouble
_ncm_stats_dist_vkde_get_href (NcmStatsDist *sd)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);

  if (self->use_rot_href)
  {
    /* Chain up : start */
    const gdouble href_base = NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->get_href (sd);

    return href_base / self->local_frac;
  }
  else
  {
    NcmStatsDistPrivate * const ppself = ncm_stats_dist_get_instance_private (sd);

    return ppself->over_smooth;
  }
}

static void
_cholesky_decomp (NcmMatrix *cov_decomp, NcmMatrix *cov, const guint d, const guint maxiter)
{
  ncm_matrix_memcpy (cov_decomp, cov);

  if (ncm_matrix_cholesky_decomp (cov_decomp, 'U') != 0)
  {
    ncm_matrix_memcpy (cov_decomp, cov);

    if (ncm_matrix_nearPD (cov_decomp, 'U', TRUE, maxiter) != 0)
    {
      guint i;

      ncm_matrix_set_zero (cov_decomp);

      for (i = 0; i < d; i++)
      {
        ncm_matrix_set (cov_decomp, i, i, ncm_matrix_get (cov, i, i));
      }

      g_assert_cmpint (ncm_matrix_cholesky_decomp (cov_decomp, 'U'), ==, 0);
    }
  }
}

static void
_ncm_stats_dist_vkde_build_cov_array_kdtree (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistKDEPrivate * const pself = ncm_stats_dist_kde_get_instance_private (NCM_STATS_DIST_KDE (sd));
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistKernel *kernel           = ncm_stats_dist_peek_kernel (sd);

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
    const size_t k = GSL_MAX (self->local_frac * ppself->n_obs, 2);

    #pragma omp parallel for schedule(dynamic, 1) if (ppself->use_threads)

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
          gdouble lnnorm_i;

          _cholesky_decomp (cov_decomp, sample_cov, ppself->d, pself->nearPD_maxiter);
          ncm_matrix_free (sample_cov);

          lnnorm_i = ncm_stats_dist_kernel_get_lnnorm (kernel, cov_decomp);

          ncm_vector_set (self->lnnorms, i, lnnorm_i);
        }
      }

      ncm_stats_vec_reset (sample, TRUE);
      ncm_memory_pool_return (sample_ptr);
    }
  }
  kdtree_destroy (tree);
}

static void
_ncm_stats_dist_vkde_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  NcmStatsDistPrivate * const ppself   = ncm_stats_dist_get_instance_private (sd);

  if (self->local_frac * ppself->n_obs < 2)
    g_error ("Too few observations.\n"
             "\tThe number of observations is too small to use the local covariance method.\n"
             "\tThe local fraction = %f times number of observations %d is less than 2.",
             self->local_frac,
             ppself->n_obs);

  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->prepare_kernel (sd, sample_array);
  _ncm_stats_dist_vkde_build_cov_array_kdtree (sd, sample_array);
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

  #pragma omp parallel if (ppself->use_threads)
  {
    NcmMatrix *invUsample_matrix = ncm_matrix_new (ncm_matrix_col_len (pself->sample_matrix), ncm_matrix_row_len (pself->sample_matrix));

    #pragma omp for schedule(dynamic, 1)

    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
      NcmVector *theta_i      = g_ptr_array_index (ppself->sample_array, i);
      gint ret;
      guint j;

      /*ncm_matrix_memcpy (pself->invUsample_matrix, pself->sample_matrix); */
      ncm_matrix_memcpy (invUsample_matrix, pself->sample_matrix);

      /*#pragma omp parallel for if (ppself->use_threads) */

      for (j = 0; j < ppself->n_obs; j++)
      {
        /*NcmVector *theta_j = g_ptr_array_index (pself->invUsample_array, j); */
        gdouble *theta_j = ncm_matrix_ptr (invUsample_matrix, j, 0);

        /*ncm_vector_axpy (theta_j, -1.0, theta_i); */
        cblas_daxpy (ppself->d, -1.0, ncm_vector_const_data (theta_i), 1, theta_j, 1);
      }

      /*
       *  ret = gsl_blas_dtrsm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
       *                     1.0, ncm_matrix_gsl (cov_decomp_i),
       *                     ncm_matrix_gsl (pself->invUsample_matrix));
       */
      ret = gsl_blas_dtrsm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                            1.0, ncm_matrix_gsl (cov_decomp_i),
                            ncm_matrix_gsl (invUsample_matrix));
      NCM_TEST_GSL_RESULT ("_ncm_stats_dist_vkde_compute_IM", ret);

      /* #pragma omp parallel for if (ppself->use_threads) */

      for (j = 0; j < ppself->n_obs; j++)
      {
        /*NcmVector *theta_j = g_ptr_array_index (pself->invUsample_array, j); */
        gdouble *theta_j = ncm_matrix_ptr (invUsample_matrix, j, 0);
        gdouble chi2_ij;

        /*chi2_ij = ncm_vector_dot (theta_j, theta_j) * one_href2; */
        chi2_ij = cblas_ddot (ppself->d, theta_j, 1, theta_j, 1) * one_href2;

        ncm_matrix_set (IM, j, i, chi2_ij);
      }
    }

    ncm_matrix_free (invUsample_matrix);
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
  gint ret;
  guint i;

  {
    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
      NcmVector *theta_i      = g_ptr_array_index (ppself->sample_array, i);

      ncm_vector_memcpy (ev->delta_x, x);
      ncm_vector_axpy (ev->delta_x, -1.0, theta_i);

      ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (cov_decomp_i), ncm_vector_gsl (ev->delta_x));
      NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_vbk_studentt_eval", ret);

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
    const gdouble w_i  = ncm_vector_fast_get (ppself->weights, i);

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
  gint ret;
  guint i;

  {
    for (i = 0; i < ppself->n_kernels; i++)
    {
      NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
      NcmVector *theta_i      = g_ptr_array_index (ppself->sample_array, i);

      ncm_vector_memcpy (ev->delta_x, x);
      ncm_vector_axpy (ev->delta_x, -1.0, theta_i);

      ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (cov_decomp_i), ncm_vector_gsl (ev->delta_x));
      NCM_TEST_GSL_RESULT ("_ncm_stats_dist_vkde_eval_weights_m2lnp", ret);

      {
        const gdouble chi2_i = ncm_vector_dot (ev->delta_x, ev->delta_x) * one_href2;

        ncm_vector_fast_set (ev->chi2, i, chi2_i);
      }
    }
  }

  {
    gdouble gamma, lambda;

    ncm_stats_dist_kernel_eval_sum0_gamma_lambda (ppself->kernel, ev->chi2, ppself->weights, self->lnnorms, ev->lnK, &gamma, &lambda);

    ncm_memory_pool_return (ev_ptr);

    return -2.0 * (gamma + log1p (lambda) - ppself->d * log (ppself->href));
  }
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


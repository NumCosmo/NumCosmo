/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_vkde.c
 *
 *  Wed November 07 16:02:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_vkde.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * @title: NcmStatsDist
 * @short_description: Abstract class for implementing N dimensional probability distributions FIXME
 *
 * Abstract class to reconstruct an arbitrary N dimensional probability distribution FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_vkde.h"
#include "math/ncm_iset.h"
#include "math/ncm_nnls.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "misc/kdtree.h"
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "math/ncm_stats_dist_vkde_private.h"

enum
{
  PROP_0,
  PROP_LOCAL_FRAC,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDistVKDE, ncm_stats_dist_vkde, NCM_TYPE_STATS_DIST_KDE);

static void
ncm_stats_dist_vkde_init (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  
  self->sample       = NULL;
  self->cov_array    = g_ptr_array_new ();
  self->tmp_cov      = NULL;
  self->lnnorm_array = g_array_new (FALSE, FALSE, sizeof (gdouble));
  
  self->local_frac = 0.0;
  
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_vkde_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistVKDE *sdvkde = NCM_STATS_DIST_VKDE (object);
  
  /*NcmStatsDistVKDEPrivate * const self = sdvkde->priv;*/
  
  g_return_if_fail (NCM_IS_STATS_DIST_VKDE (object));
  
  switch (prop_id)
  {
    case PROP_LOCAL_FRAC:
      g_value_set_double (value, ncm_stats_dist_vkde_get_local_frac (sdvkde));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_vkde_dispose (GObject *object)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (object);
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  
  ncm_stats_vec_clear (&self->sample);
  ncm_matrix_clear (&self->tmp_cov);
  
  g_clear_pointer (&self->cov_array, g_ptr_array_unref);
  g_clear_pointer (&self->lnnorm_array, g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_vkde_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_vkde_finalize (GObject *object)
{
  /*NcmStatsDistVKDE *sdvkde                = NCM_STATS_DIST_VKDE (object);*/
  /*NcmStatsDistVKDEPrivate * const self = sdvkde->priv;*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_vkde_parent_class)->finalize (object);
}

static void _ncm_stats_dist_vkde_set_dim (NcmStatsDist *sd, const guint dim);
static void _ncm_stats_dist_vkde_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array);
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
  
  base_class->set_dim        = &_ncm_stats_dist_vkde_set_dim;
  base_class->prepare_kernel = &_ncm_stats_dist_vkde_prepare_kernel;
  base_class->reset          = &_ncm_stats_dist_vkde_reset;
}

static void
_ncm_stats_dist_vkde_set_dim (NcmStatsDist *sd, const guint dim)
{
  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->set_dim  (sd, dim);
  {
    NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
    NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
    
    ncm_stats_vec_clear (&self->sample);
    ncm_matrix_clear (&self->tmp_cov);
    
    g_ptr_array_set_size (self->cov_array, 0);
    
    self->sample  = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, TRUE);
    self->tmp_cov = ncm_matrix_new (dim, dim);
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
      gint i;
      
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
_ncm_stats_dist_vkde_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array)
{
  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->prepare_kernel (sd, sample_array);
  {
    NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
    NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
    const guint n                        = sample_array->len;
    const guint dim                      = ncm_stats_dist_get_dim (NCM_STATS_DIST (sdvkde));
    GPtrArray *invU_sample               = ncm_stats_dist_kde_peek_invU_sample (NCM_STATS_DIST_KDE (sd));
    gint i;
    
    /*
     * Creates a near tree object and add all transformed vector.
     */
    {
      struct kdtree *tree = kdtree_init (dim);
      
      for (i = 0; i < n; i++)
      {
        NcmVector *invUtheta_i = g_ptr_array_index (invU_sample, i);
        
        /*
         * Inserting the transformed vector in the tree, saving also the index.
         */
        kdtree_insert (tree, ncm_vector_data (invUtheta_i));
      }
      
      kdtree_rebuild (tree);
      
      /*
       * Lets find the k nearest neighbors of each vector in the
       * sample and use them to define the local covariance at each
       * vector location.
       */
      {
        const size_t k = self->local_frac * sample_array->len;
        
        g_ptr_array_set_size (self->cov_array, 0);
        
        for (i = 0; i < n; i++)
        {
          NcmVector *invUtheta_i = g_ptr_array_index (invU_sample, i);
          
          kdtree_knn_search_clean (tree);
          kdtree_knn_search (tree, ncm_vector_data (invUtheta_i), k);
          
          ncm_stats_vec_reset (self->sample, TRUE);
          {
            struct knn_list *p = tree->knn_list_head.next;
            
            while (p != &tree->knn_list_head)
            {
              NcmVector *ni = g_ptr_array_index (sample_array, p->node->coord_index);
              
              ncm_stats_vec_append (self->sample, ni, FALSE);
              p = p->next;
            }
          }
          
          /*
           * Saving the covariance for each vector.
           */
          {
            NcmMatrix *sample_cov = ncm_matrix_dup (ncm_stats_vec_peek_cov_matrix (self->sample, 0));
            
            g_ptr_array_add (self->cov_array, sample_cov);
          }
        }
        
        /*
         * Computing the bandwidth for each local kernel.
         */
        {
          /*self->href = self->over_smooth / self->local_frac * ncm_stats_dist_kernel_get_rot_bandwidth (kernel, n);*/
        }
      }
      kdtree_destroy (tree);
    }
    
    {
      NcmStatsDistKernel *kernel = ncm_stats_dist_peek_kernel (sd);
      const guint nearPD_maxiter = ncm_stats_dist_kde_get_nearPD_maxiter (NCM_STATS_DIST_KDE (sd));
      const guint href           = ncm_stats_dist_get_href (sd);
      
      g_array_set_size  (self->lnnorm_array, n);
      
      for (i = 0; i < n; i++)
      {
        NcmMatrix *cov_i = g_ptr_array_index (self->cov_array, i);
        gdouble lnnorm_i;
        
        ncm_matrix_memcpy (self->tmp_cov, cov_i);
        
        _cholesky_decomp (cov_i, self->tmp_cov, dim, nearPD_maxiter);
        
        lnnorm_i = ncm_stats_dist_kernel_get_lnnorm (kernel, cov_i, href);
        
        g_array_index (self->lnnorm_array, gdouble, i) = lnnorm_i;
      }
    }
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
 * Sets local kernel fraction to @local_frac.
 *
 */
void
ncm_stats_dist_vkde_set_local_frac (NcmStatsDistVKDE *sdvkde, const gdouble local_frac)
{
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  
  g_assert_cmpfloat (local_frac, >=, 0.001);
  g_assert_cmpfloat (local_frac, <=, 1.0);
  
  self->local_frac = local_frac;
}

/**
 * ncm_stats_dist_vkde_get_local_frac:
 * @sdvkde: a #NcmStatsDistVKDE
 *
 * Returns: the local kernel fraction.
 */
gdouble
ncm_stats_dist_vkde_get_local_frac (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  
  return self->local_frac;
}

/**
 * ncm_stats_dist_vkde_get_Ki:
 * @sdvkde: a #NcmStatsDistVKDE
 * @i: kernel index
 * @y_i: (out callee-allocates) (transfer full): kernel location
 * @cov_i: (out callee-allocates) (transfer full): kernel covariance U
 * @n_i: (out): kernel normalization
 * @w_i: (out): kernel weight
 *
 * Return all information about the @i-th kernel.
 *
 */
void
ncm_stats_dist_vkde_get_Ki (NcmStatsDistVKDE *sdvkde, const guint i, NcmVector **y_i, NcmMatrix **cov_i, gdouble *n_i, gdouble *w_i)
{
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  GPtrArray *sample_array              = ncm_stats_dist_peek_sample_array (NCM_STATS_DIST (sdvkde));
  NcmVector *weights                   = ncm_stats_dist_peek_weights (NCM_STATS_DIST (sdvkde));
  
  g_assert (i < ncm_stats_dist_get_sample_size (NCM_STATS_DIST (sdvkde)));
  
  y_i[0]   = ncm_vector_dup (g_ptr_array_index (sample_array, i));
  cov_i[0] = ncm_matrix_dup (g_ptr_array_index (self->cov_array, i));
  n_i[0]   = exp (g_array_index (self->lnnorm_array, gdouble, i));
  w_i[0]   = ncm_vector_get (weights, i);
  
  ncm_matrix_triang_to_sym (g_ptr_array_index (self->cov_array, i), 'U', TRUE, cov_i[0]);
}


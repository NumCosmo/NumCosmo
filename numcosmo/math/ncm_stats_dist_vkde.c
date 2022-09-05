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
#include "math/ncm_nnls.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "misc/kdtree.h"
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "math/ncm_stats_dist_vkde_private.h"
#include "math/ncm_stats_dist_kde_private.h"
#include "math/ncm_stats_dist_private.h"

enum
{
  PROP_0,
  PROP_LOCAL_FRAC,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmStatsDistVKDE, ncm_stats_dist_vkde, NCM_TYPE_STATS_DIST_KDE);

static void
ncm_stats_dist_vkde_init (NcmStatsDistVKDE *sdvkde)
{
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv = ncm_stats_dist_vkde_get_instance_private (sdvkde);
  
  self->sample    = NULL;
  self->cov_array = g_ptr_array_new ();
  self->tmp_cov   = NULL;
  self->delta_x   = NULL;
  self->lnnorms   = NULL;
  
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
  ncm_vector_clear (&self->delta_x);
  ncm_vector_clear (&self->lnnorms);
  
  g_clear_pointer (&self->cov_array, g_ptr_array_unref);
  
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
static gdouble _ncm_stats_dist_vkde_get_href (NcmStatsDist *sd);
static void _ncm_stats_dist_vkde_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array, NcmVector *m2lnp);
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
    NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
    
    ncm_stats_vec_clear (&self->sample);
    ncm_matrix_clear (&self->tmp_cov);
    ncm_vector_clear (&self->delta_x);
    
    g_ptr_array_set_size (self->cov_array, 0);
    
    self->sample  = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, FALSE);
    self->tmp_cov = ncm_matrix_new (dim, dim);
    self->delta_x = ncm_vector_new (dim);
  }
}

static gdouble
_ncm_stats_dist_vkde_get_href (NcmStatsDist *sd)
{
  /* Chain up : start */
  const gdouble href_base = NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->get_href (sd);
  {
    NcmStatsDistVKDEPrivate * const self = NCM_STATS_DIST_VKDE (sd)->priv;
    
    /*printf ("href % 22.15g % 22.15g % 22.15g\n", href_base / self->local_frac, href_base, self->local_frac);*/

    return href_base / self->local_frac;
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
_ncm_stats_dist_vkde_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array, NcmVector *m2lnp)
{
  /* Chain up : start */
  NCM_STATS_DIST_CLASS (ncm_stats_dist_vkde_parent_class)->prepare_kernel (sd, sample_array, m2lnp);
  {
    NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
    NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
    NcmStatsDistKDEPrivate * const pself = NCM_STATS_DIST_KDE (sd)->priv;
    NcmStatsDistPrivate * const ppself   = sd->priv;
    gint i;
    
    /*
     * Creates a near tree object and add all transformed vector.
     */
    if (m2lnp == NULL || FALSE)
    {
      struct kdtree *tree = kdtree_init (ppself->d);
      
      for (i = 0; i < ppself->n; i++)
      {
        NcmVector *invUtheta_i = g_ptr_array_index (pself->invUsample, i);
        
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
        
        for (i = 0; i < ppself->n; i++)
        {
          NcmVector *invUtheta_i = g_ptr_array_index (pself->invUsample, i);
          NcmVector *theta_i = g_ptr_array_index (sample_array, i);
          
          kdtree_knn_search_clean (tree);
          kdtree_knn_search (tree, ncm_vector_data (invUtheta_i), k);

          //printf ("points close to: \n");
          //ncm_vector_log_vals (g_ptr_array_index (sample_array, i), "CE : ", "% 12.5g", TRUE);
          //ncm_vector_log_vals (invUtheta_i, "CET: ", "% 12.5g", TRUE);
          
          ncm_stats_vec_reset (self->sample, TRUE);
          ncm_matrix_set_zero (self->tmp_cov);
          {
            struct knn_list *p = tree->knn_list_head.next;
            gdouble norma = 0.0;
            
            /*ncm_message ("#-----------------------------------------\n");*/
            while (p != &tree->knn_list_head)
            {
              NcmVector *ni = g_ptr_array_index (sample_array, p->node->coord_index);
              gint l;
              /*ncm_message ("dist: % 22.15g: ", sqrt (p->distance));*/

              //ncm_vector_log_vals (ni, "NI : ", "% 12.5g", TRUE);
              //ncm_vector_log_vals (g_ptr_array_index (pself->invUsample, p->node->coord_index), "NIT: ", "% 12.5g", TRUE);
              
              ncm_stats_vec_append (self->sample, ni, FALSE);

              for (l = 0; l < ppself->d; l++)
              {
                const gdouble dx_l = ncm_vector_get (ni, l) - ncm_vector_get (theta_i, l);
                gint m;
                /*ncm_message ("% 12.5g ", p->node->coord[l] - ncm_vector_get (invUtheta_i, l));*/

                ncm_matrix_addto (self->tmp_cov, l, l, dx_l * dx_l);
                for (m = l + 1; m < ppself->d; m++)
                {
                  const gdouble dx_m   = ncm_vector_get (ni, m) - ncm_vector_get (theta_i, m);
                  const gdouble cov_lm = dx_l * dx_m;

                  ncm_matrix_addto (self->tmp_cov, l, m, cov_lm);
                  ncm_matrix_addto (self->tmp_cov, m, l, cov_lm);
                }
              }
              /*ncm_vector_log_vals (ni, "V: ", "% 12.5g", TRUE);*/

              /*ncm_message ("\n");*/


              p = p->next;
              norma++;
            }
            ncm_matrix_scale (self->tmp_cov, 1.0 / norma);
          }
          
          /*
           * Saving the covariance for each vector.
           */
          {
            NcmMatrix *sample_cov = ncm_matrix_dup (ncm_stats_vec_peek_cov_matrix (self->sample, 0));
            
            //ncm_matrix_log_vals (sample_cov, "S1", "% 12.5e");
            //ncm_matrix_log_vals (self->tmp_cov, "S2", "% 12.5e");
            
            ncm_matrix_memcpy (sample_cov, self->tmp_cov);

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
    else /* Uses m2lnp as distances to determine the closest points*/
    {
      const size_t k = self->local_frac * sample_array->len;
      glong i;

      g_array_set_size (ppself->m2lnp_sort, ppself->n);
      gsl_sort_index (&g_array_index (ppself->m2lnp_sort, size_t, 0),
          ncm_vector_data (m2lnp),
          ncm_vector_stride (m2lnp),
          ncm_vector_len (m2lnp));

      g_ptr_array_set_size (self->cov_array, ppself->n);

      for (i = 0; i < ppself->n; i++)
      {
        const glong l        = g_array_index (ppself->m2lnp_sort, size_t, i);
        const gdouble m2lnp_l = ncm_vector_get (m2lnp, l);
        glong up = 1;
        glong down = 1;
        glong count = 0;

        ncm_stats_vec_reset (self->sample, TRUE);

        ncm_stats_vec_append (self->sample, g_ptr_array_index (sample_array, l), FALSE);
        /*printf ("# Preparing kernel %ld, m2lnp = % 22.15g\n", l, m2lnp_l);*/

        while (TRUE)
        {
          const glong i_up   = i + up;
          const glong i_down = i - down;

          if ((i_up < ppself->n) && (i_down >= 0))
          {
            const glong l_up        = g_array_index (ppself->m2lnp_sort, size_t, i_up);
            const glong l_down      = g_array_index (ppself->m2lnp_sort, size_t, i_down);
            const gdouble m2lnp_up   = ncm_vector_get (m2lnp, l_up);
            const gdouble m2lnp_down = ncm_vector_get (m2lnp, l_down);

/*
            printf ("# Testing up and down (%ld %ld) (%ld, %ld) (% 22.15g % 22.15g)\n",
                i_up, i_down, l_up, l_down,
                fabs (m2lnp_up - m2lnp_l), fabs (m2lnp_down - m2lnp_l));
*/

            if (fabs (m2lnp_up - m2lnp_l) < fabs (m2lnp_down - m2lnp_l))
            {
              NcmVector *n_l_up = g_ptr_array_index (sample_array, l_up);
              ncm_stats_vec_append (self->sample, n_l_up, FALSE);
              /*printf ("# Up smaller, adding it! (%ld %ld)\n", up, count);*/

              up++;
              count++;
            }
            else
            {
              NcmVector *n_l_down = g_ptr_array_index (sample_array, l_down);
              ncm_stats_vec_append (self->sample, n_l_down, FALSE);
              /*printf ("# Down smaller, adding it! (%ld %ld)\n", down, count);*/

              down++;
              count++;
            }
          }
          else if (i_down < 0)
          {
            const glong l_up = g_array_index (ppself->m2lnp_sort, size_t, i_up);
            NcmVector *n_l_up = g_ptr_array_index (sample_array, l_up);
            ncm_stats_vec_append (self->sample, n_l_up, FALSE);

            /*printf ("# Up smaller, adding it! (%ld %ld)\n", up, count);*/

            up++;
            count++;
          }
          else if (i_up >= ppself->n)
          {
            const glong l_down = g_array_index (ppself->m2lnp_sort, size_t, i_down);
            NcmVector *n_l_down = g_ptr_array_index (sample_array, l_down);
            ncm_stats_vec_append (self->sample, n_l_down, FALSE);

            /*printf ("# Down smaller, adding it! (%ld %ld)\n", down, count);*/
            down++;
            count++;
          }
          else
            g_assert_not_reached ();

          if (count >= k)
            break;
        }

        {
          NcmMatrix *sample_cov = ncm_matrix_dup (ncm_stats_vec_peek_cov_matrix (self->sample, 0));
          g_ptr_array_index (self->cov_array, l) = sample_cov;
        }

      }
    }
    
    {
      NcmStatsDistKernel *kernel = ncm_stats_dist_peek_kernel (sd);
      
      if ((self->lnnorms == NULL) || (ncm_vector_len (self->lnnorms) != ppself->n))
      {
        ncm_vector_clear (&self->lnnorms);
        self->lnnorms = ncm_vector_new (ppself->n);
      }
      
      for (i = 0; i < ppself->n; i++)
      {
        NcmMatrix *cov_i = g_ptr_array_index (self->cov_array, i);
        gdouble lnnorm_i;
        
        if (i < 10)
        {
          ncm_cfg_msg_sepa ();
          ncm_vector_log_vals (g_ptr_array_index (sample_array, i), "THETA: ", "% 12.5g", TRUE);
          ncm_matrix_log_vals (cov_i, "COV:   ", "% 12.5g");
        }

        ncm_matrix_memcpy (self->tmp_cov, cov_i);
        
        _cholesky_decomp (cov_i, self->tmp_cov, ppself->d, pself->nearPD_maxiter);
        
        lnnorm_i = ncm_stats_dist_kernel_get_lnnorm (kernel, cov_i);
        
        ncm_vector_set (self->lnnorms, i, lnnorm_i);
      }
    }
  }
}

static void
_ncm_stats_dist_vkde_compute_IM (NcmStatsDist *sd, NcmMatrix *IM)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  /*NcmStatsDistKDEPrivate * const pself = NCM_STATS_DIST_KDE (sd)->priv;*/
  NcmStatsDistPrivate * const ppself = sd->priv;
  const gdouble href2                = ppself->href * ppself->href;
  
  gint i, j, res;
  
  for (i = 0; i < ppself->n; i++)
  {
    NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
    NcmVector *theta_i      = g_ptr_array_index (ppself->sample_array, i);
    
    for (j = 0; j < ppself->n; j++)
    {
      NcmVector *theta_j = g_ptr_array_index (ppself->sample_array, j);
      gdouble chi2_ij;
      
      ncm_vector_memcpy (self->delta_x, theta_i);
      ncm_vector_sub (self->delta_x, theta_j);
      
      res = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (cov_decomp_i), ncm_vector_gsl (self->delta_x));
      NCM_TEST_GSL_RESULT ("_ncm_stats_dist_vkde_compute_IM", res);
      
      chi2_ij = ncm_vector_dot (self->delta_x, self->delta_x) / href2;
      /*p_ij = _ncm_stats_dist_nd_vbk_studentt_f (self, d, m2lnp_ij) / norm_i;*/
      
      ncm_matrix_set (IM, j, i, chi2_ij);
    }
  }
  
  {
    NcmVector *IMv            = ncm_matrix_as_vector (IM);
    const gdouble lnnorm_href = ppself->d * log (ppself->href);
    
    ncm_stats_dist_kernel_eval_unnorm_vec (ppself->kernel, IMv, IMv);
    
    for (i = 0; i < ppself->n; i++)
    {
      const gdouble norm_i = exp (ncm_vector_fast_get (self->lnnorms, i) + lnnorm_href);
      
      ncm_matrix_mul_col (IM, i, 1.0 / norm_i);
    }
    
    ncm_vector_free (IMv);
  }
}

static NcmMatrix *
_ncm_stats_dist_vkde_peek_cov_decomp (NcmStatsDist *sd, guint i)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  
  g_assert (i < self->cov_array->len);
  
  return g_ptr_array_index (self->cov_array, i);
}

static gdouble
_ncm_stats_dist_vkde_get_lnnorm (NcmStatsDist *sd, guint i)
{
  NcmStatsDistVKDE *sdvkde             = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  NcmStatsDistPrivate * const ppself   = sd->priv;
  
  g_assert (i < self->cov_array->len);
  
  return ncm_vector_fast_get (self->lnnorms, i) + ppself->d * log (ppself->href);
}

static gdouble
_ncm_stats_dist_vkde_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  NcmStatsDistVKDE *sdvkde = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  NcmStatsDistKDEPrivate * const pself = NCM_STATS_DIST_KDE (sd)->priv;
  NcmStatsDistPrivate * const ppself = sd->priv;
  const gdouble href2 = ppself->href * ppself->href;
  gdouble s = 0.0;
  gint i, ret;
  
  for (i = 0; i < ppself->n; i++)
  {
    NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
    NcmVector *theta_i      = g_ptr_array_index (ppself->sample_array, i);
    
    ncm_vector_memcpy (self->delta_x, x);
    ncm_vector_sub (self->delta_x, theta_i);
    
    ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (cov_decomp_i), ncm_vector_gsl (self->delta_x));
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_vbk_studentt_eval", ret);
    
    {
      const gdouble chi2_i = ncm_vector_dot (self->delta_x, self->delta_x) / href2;
      
      ncm_vector_fast_set (pself->chi2, i, chi2_i);
    }
  }
  
  ncm_stats_dist_kernel_eval_unnorm_vec (ppself->kernel, pself->chi2, pself->chi2);
  
  for (i = 0; i < ppself->n; i++)
  {
    const gdouble Ku_i = ncm_vector_fast_get (pself->chi2, i);
    const gdouble u_i  = exp (ncm_vector_fast_get (self->lnnorms, i));
    const gdouble w_i  = ncm_vector_fast_get (ppself->weights, i);
    
    s += w_i * (Ku_i / u_i);
  }
  
  return s / pow (ppself->href, ppself->d);
}

static gdouble
_ncm_stats_dist_vkde_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  NcmStatsDistVKDE *sdvkde = NCM_STATS_DIST_VKDE (sd);
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  NcmStatsDistKDEPrivate * const pself = NCM_STATS_DIST_KDE (sd)->priv;
  NcmStatsDistPrivate * const ppself = sd->priv;
  const gdouble href2 = ppself->href * ppself->href;
  gint i, ret;
  
  for (i = 0; i < ppself->n; i++)
  {
    NcmMatrix *cov_decomp_i = g_ptr_array_index (self->cov_array, i);
    NcmVector *theta_i      = g_ptr_array_index (ppself->sample_array, i);
    
    ncm_vector_memcpy (self->delta_x, x);
    ncm_vector_sub (self->delta_x, theta_i);
    
    ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (cov_decomp_i), ncm_vector_gsl (self->delta_x));
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_nd_vbk_studentt_eval", ret);
    
    {
      const gdouble chi2_i = ncm_vector_dot (self->delta_x, self->delta_x) / href2;
      
      ncm_vector_fast_set (pself->chi2, i, chi2_i);
    }
  }
  
  {
    gdouble gamma, lambda;
    
    ncm_stats_dist_kernel_eval_sum0_gamma_lambda (ppself->kernel, pself->chi2, ppself->weights, self->lnnorms, &gamma, &lambda);

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
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  
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
  NcmStatsDistVKDEPrivate * const self = sdvkde->priv;
  
  return self->local_frac;
}


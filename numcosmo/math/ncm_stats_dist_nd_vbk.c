/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd_vbk.c
 *
 *  Wed November 07 16:02:36 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_vbk.c
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
 * SECTION:ncm_stats_dist_nd_vbk
 * @title: NcmStatsDistNd
 * @short_description: Abstract class for implementing N dimensional probability distributions FIXME
 *
 * Abstract class to reconstruct an arbitrary N dimensional probability distribution FIXME
 *
 *  
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_stats_dist_nd_vbk.h"
#include "math/ncm_iset.h"
#include "math/ncm_nnls.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "misc/CVector.h"
#include "misc/CNearTree.h"
#include "levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmStatsDistNdVBKPrivate
{
  GPtrArray *sample_array;
  GPtrArray *cov_array;
  GArray *norm_array;
  GPtrArray *IM_array;
  
  NcmStatsVec *sample;
  NcmMatrix *cov_decomp;
  NcmMatrix *log_cov;
  NcmMatrix *href_cov;
  NcmMatrix *sample_matrix;
  GPtrArray *invUsample;
  NcmVector *weights;
  NcmVector *v;
  gdouble over_smooth;
  NcmStatsDistNdVBKCV cv_type;
  gdouble split_frac;
  NcmVector *href;
  gdouble kernel_lnnorm;
  gdouble min_m2lnp;
  gdouble max_m2lnp;
  NcmStatsVec *m2lnp_stats;
  gdouble rnorm;
  guint n;
  guint alloc_n;
  guint d;
  GArray *sampling;
  guint nearPD_maxiter;
  NcmNNLS *nnls;
  NcmMatrix *IM;
  NcmMatrix *sub_IM;
  NcmVector *sub_x;
  NcmVector *f;
  gsl_min_fminimizer *min;
  gsl_multimin_fminimizer *mmin;
  gsl_multimin_fminimizer *mmin_fitd;
  guint mmin_d;
  gdouble *levmar_workz;
  guint levmar_n;
};

enum
{
  PROP_0,
  PROP_DIM,
  PROP_SAMPLE_SIZE,
  PROP_OVER_SMOOTH,
  PROP_NEARPD_MAXITER,
  PROP_CV_TYPE,
  PROP_SPLIT_FRAC,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDistNdVBK, ncm_stats_dist_nd_vbk, G_TYPE_OBJECT);

static void
ncm_stats_dist_nd_vbk_init (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv = ncm_stats_dist_nd_vbk_get_instance_private (dnd);

  self->sample_array     = g_ptr_array_new ();
  self->norm_array    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->cov_array 	 = g_ptr_array_new ();
  self->IM_array	 = g_ptr_array_new ();
  
  self->sample           = NULL;
  self->cov_decomp       = NULL;
  self->log_cov          = NULL;
  self->href_cov         = NULL;
  self->sample_matrix    = NULL;
  self->invUsample       = g_ptr_array_new ();
  self->weights          = NULL;
  self->v                = NULL;
  self->over_smooth      = 0.0;
  self->cv_type          = NCM_STATS_DIST_ND_VBK_CV_LEN;
  self->split_frac       = 0.0;
  self->href             = NULL;
  self->kernel_lnnorm    = 0.0;
  self->min_m2lnp        = 0.0;
  self->max_m2lnp        = 0.0;
  self->m2lnp_stats      = ncm_stats_vec_new (1, NCM_STATS_VEC_MEAN, FALSE);
  self->rnorm            = 0.0;
  self->n                = 0;
  self->alloc_n          = 0;
  self->d                = 0;
  self->sampling         = g_array_new (FALSE, FALSE, sizeof (guint));
  self->nearPD_maxiter   = 0;
  self->nnls             = NULL;
  self->IM               = NULL;
  self->sub_IM           = NULL;
  self->sub_x            = NULL;
  self->f                = NULL;
  self->min              = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  self->mmin             = gsl_multimin_fminimizer_alloc ( gsl_multimin_fminimizer_nmsimplex2, 1);
  self->mmin_fitd        = NULL;
  self->mmin_d           = 0;
  self->levmar_workz     = NULL;
  self->levmar_n         = 0;

  ncm_stats_vec_enable_quantile (self->m2lnp_stats, 0.9999); /* 0.6827, 0.9545, 0.9973 */

  g_ptr_array_set_free_func (self->sample_array, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->cov_array, (GDestroyNotify) ncm_matrix_free);
  g_ptr_array_set_free_func (self->invUsample, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->IM_array, (GDestroyNotify) ncm_matrix_free);
}

static void _ncm_stats_dist_nd_vbk_set_dim (NcmStatsDistNdVBK *dnd, const guint dim);

static void
_ncm_stats_dist_nd_vbk_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNdVBK *dnd = NCM_STATS_DIST_ND_VBK (object);
  /*g_return_if_fail (NCM_IS_STATS_DIST_ND (object));*/

  switch (prop_id)
  {
    case PROP_DIM:
      _ncm_stats_dist_nd_vbk_set_dim (dnd, g_value_get_uint (value));
      break;
    case PROP_SAMPLE_SIZE:
      g_assert_not_reached ();
      break;
    case PROP_OVER_SMOOTH:
      ncm_stats_dist_nd_vbk_set_over_smooth (dnd, g_value_get_double (value));
      break;
    case PROP_NEARPD_MAXITER:
      ncm_stats_dist_nd_vbk_set_nearPD_maxiter (dnd, g_value_get_uint (value));
      break;
    case PROP_CV_TYPE:
      ncm_stats_dist_nd_vbk_set_cv_type (dnd, g_value_get_enum (value));
      break;
    case PROP_SPLIT_FRAC:
      ncm_stats_dist_nd_vbk_set_split_frac (dnd, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_vbk_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNdVBK *dnd = NCM_STATS_DIST_ND_VBK (object);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND_VBK (object));

  switch (prop_id)
  {
    case PROP_DIM:
      g_value_set_uint (value, self->d);
      break;
    case PROP_SAMPLE_SIZE:
      g_value_set_uint (value, ncm_stats_vec_nitens (self->sample));
      break;
    case PROP_OVER_SMOOTH:
      g_value_set_double (value, ncm_stats_dist_nd_vbk_get_over_smooth (dnd));
      break;
    case PROP_NEARPD_MAXITER:
      g_value_set_uint (value, self->nearPD_maxiter);
      break;
    case PROP_CV_TYPE:
      g_value_set_enum (value, ncm_stats_dist_nd_vbk_get_cv_type (dnd));
      break;
    case PROP_SPLIT_FRAC:
      g_value_set_double (value, ncm_stats_dist_nd_vbk_get_split_frac (dnd));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_vbk_dispose (GObject *object)
{
  NcmStatsDistNdVBK *dnd = NCM_STATS_DIST_ND_VBK (object);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  g_clear_pointer (&self->sample_array, g_ptr_array_unref);
  g_clear_pointer (&self->cov_array, g_ptr_array_unref);
  g_clear_pointer (&self->norm_array, g_array_unref);

  ncm_stats_vec_clear (&self->sample);
  ncm_matrix_clear (&self->cov_decomp);
  ncm_matrix_clear (&self->log_cov);
  ncm_matrix_clear (&self->href_cov);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_vector_clear (&self->href);
  ncm_vector_clear (&self->weights);
  ncm_vector_clear (&self->v);

  g_clear_pointer (&self->sampling, g_array_unref);
  g_clear_pointer (&self->invUsample, g_ptr_array_unref);
  g_clear_pointer (&self->IM_array, g_ptr_array_unref);

  ncm_stats_vec_clear (&self->m2lnp_stats);

  ncm_nnls_clear (&self->nnls);

  ncm_matrix_clear (&self->IM);
  ncm_matrix_clear (&self->sub_IM);
  ncm_vector_clear (&self->sub_x);
  ncm_vector_clear (&self->f);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_vbk_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_vbk_finalize (GObject *object)
{
  NcmStatsDistNdVBK *dnd = NCM_STATS_DIST_ND_VBK (object);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  g_clear_pointer (&self->min,  gsl_min_fminimizer_free);
  g_clear_pointer (&self->mmin, gsl_multimin_fminimizer_free);
  g_clear_pointer (&self->mmin_fitd, gsl_multimin_fminimizer_free);
  g_clear_pointer (&self->levmar_workz, g_free);

  self->mmin_d   = 0;
  self->levmar_n = 0;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_vbk_parent_class)->finalize (object);
}

gdouble _ncm_stats_dist_nd_vbk_get_rot_bandwidth (NcmStatsDistNdVBK *dnd, const guint d, const gdouble n) { g_error ("method get_rot_bandwidth not implemented by %s.", G_OBJECT_TYPE_NAME (dnd)); return 0.0; }
gdouble _ncm_stats_dist_nd_vbk_get_kernel_lnnorm (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href) { g_error ("method get_kernel_lnnorm not implemented by %s.", G_OBJECT_TYPE_NAME (dnd)); return 0.0; }
static void _ncm_stats_dist_nd_vbk_prepare_kernel_args (NcmStatsDistNdVBK *dnd, NcmStatsVec *sample);
static void _ncm_stats_dist_nd_vbk_prepare (NcmStatsDistNdVBK *dnd);
static void _ncm_stats_dist_nd_vbk_prepare_interp (NcmStatsDistNdVBK *dnd, NcmVector *m2lnp);
static void _ncm_stats_dist_nd_vbk_reset (NcmStatsDistNdVBK *dnd);

static void
ncm_stats_dist_nd_vbk_class_init (NcmStatsDistNdVBKClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  

  object_class->set_property = &_ncm_stats_dist_nd_vbk_set_property;
  object_class->get_property = &_ncm_stats_dist_nd_vbk_get_property;
  object_class->dispose      = &_ncm_stats_dist_nd_vbk_dispose;
  object_class->finalize     = &_ncm_stats_dist_nd_vbk_finalize;

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
                                                      NCM_TYPE_STATS_DIST_ND_VBKCV, NCM_STATS_DIST_ND_VBK_CV_NONE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SPLIT_FRAC,
                                   g_param_spec_double ("split-frac",
                                                        NULL,
                                                        "Fraction to use in the split cross-validation",
                                                        0.50, 0.95, 0.9,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NEARPD_MAXITER,
                                   g_param_spec_uint ("nearPD-maxiter",
                                                      NULL,
                                                      "Maximum number of iterations in the nearPD call",
                                                      1, G_MAXUINT, 200,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_dim             = &_ncm_stats_dist_nd_vbk_set_dim;
  klass->get_rot_bandwidth   = &_ncm_stats_dist_nd_vbk_get_rot_bandwidth;
  klass->get_kernel_lnnorm   = &_ncm_stats_dist_nd_vbk_get_kernel_lnnorm;
  klass->prepare_kernel_args = &_ncm_stats_dist_nd_vbk_prepare_kernel_args;
  klass->prepare             = &_ncm_stats_dist_nd_vbk_prepare;
  klass->prepare_interp      = &_ncm_stats_dist_nd_vbk_prepare_interp;
  klass->eval                = NULL;
  klass->eval_m2lnp          = NULL;
  klass->kernel_sample       = NULL;
  klass->kernel_eval_m2lnp   = NULL;
  klass->reset               = &_ncm_stats_dist_nd_vbk_reset;
}

static void
_ncm_stats_dist_nd_vbk_set_dim (NcmStatsDistNdVBK *dnd, const guint dim)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  ncm_stats_vec_clear (&self->sample);

  ncm_matrix_clear (&self->cov_decomp);
  ncm_matrix_clear (&self->log_cov);
  ncm_matrix_clear (&self->href_cov);
  ncm_matrix_clear (&self->sample_matrix);
  ncm_vector_clear (&self->href);
  ncm_vector_clear (&self->v);

  self->d          = dim;
  self->sample     = ncm_stats_vec_new (dim, NCM_STATS_VEC_COV, TRUE);
  self->href       = ncm_vector_new (dim);
  self->cov_decomp = ncm_matrix_new (dim, dim);
  self->log_cov    = ncm_matrix_new (dim, dim);
  self->href_cov   = ncm_matrix_new (dim, dim);
  self->v          = ncm_vector_new (dim);
}

static void
_ncm_stats_dist_nd_vbk_prepare_kernel_args (NcmStatsDistNdVBK *dnd, NcmStatsVec *sample)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;
  gint i, ret, j;
  
  {
    CNearTreeHandle treehandle;
    const size_t k        = 0.2 * self->sample_array->len;
    const gdouble dRadius = 1.0e10;
    CVectorHandle results_coords;
    CVectorHandle results_objs;
    NcmVector *scale_vec_mean = ncm_vector_new (self->d);
    NcmVector *scale_vec_sd = ncm_vector_new (self->d);
    NcmVector *scale_vec_ratio = ncm_vector_new (self->d);
    NcmStatsVec *theta_t;
        
    CVectorCreate (&results_objs,   sizeof (gpointer), k);
    CVectorCreate (&results_coords, sizeof (gpointer), k);

    CNearTreeCreate (&treehandle, self->d, CNEARTREE_TYPE_DOUBLE | CNEARTREE_NORM_L2);

    theta_t = ncm_stats_vec_new (self->sample_array->len, NCM_STATS_VEC_VAR ,FALSE);
    for (j=0; j < self->sample_array->len; j++)
    {
      NcmVector *theta_j = ncm_vector_dup(g_ptr_array_index (self->sample_array, j));
      for (i=0; i < self->d; i++)
     {
     ncm_stats_vec_set(theta_t, i, ncm_vector_get(theta_j, i));
     }
     ncm_stats_vec_update(theta_t);
    }
   
    for(i = 0; i < self->d; i++)
    {
    ncm_vector_set (scale_vec_mean, i, ncm_stats_vec_get_mean (theta_t, i));
    ncm_vector_set (scale_vec_sd, i, ncm_stats_vec_get_sd (theta_t, i));
    ncm_vector_set (scale_vec_ratio, i, ncm_stats_vec_get_mean (theta_t, i) / ncm_stats_vec_get_sd (theta_t, i));
    }        
   
    for (i = 0; i < self->sample_array->len; i++)
     {
      NcmVector *theta_i = ncm_vector_dup(g_ptr_array_index (self->sample_array, i));
      ncm_vector_div (theta_i, scale_vec_sd);
      ncm_vector_sub (theta_i, scale_vec_ratio);
      CNearTreeInsert (treehandle, ncm_vector_data (theta_i), NULL);
     }
    

    for (i = 0; i < self->sample_array->len; i++)
    {
      NcmVector *theta_i = ncm_vector_dup(g_ptr_array_index (self->sample_array, i)); 
      gdouble *vdata;
      gint res, j, l;

      res = CNearTreeFindKNearest (treehandle, k, dRadius,
          results_coords, results_objs, ncm_vector_data (theta_i), TRUE);

      ncm_stats_vec_reset (self->sample, TRUE);
      for (j = 0; j < k; j++)
      {
        CVectorGetElement (results_coords, &vdata, j);

        /*printf ("vector[%2d]:", j);*/
        for (l = 0; l < self->d; l++)
        {
          ncm_stats_vec_set (self->sample, l, (vdata[l] * ncm_vector_get(scale_vec_sd,l)) + ncm_vector_get(scale_vec_mean, l));
          
        /* printf (" % 22.15g", vdata[l]);*/
        }
        ncm_stats_vec_update (self->sample);
        /*printf ("\n");*/
      }
      {
      NcmMatrix *sample_cov = ncm_matrix_dup (ncm_stats_vec_peek_cov_matrix (self->sample, 0));
      g_ptr_array_add (self->cov_array, sample_cov);
      }
    }  /*printf ("res %d\n", res);*/
    CNearTreeFree (&treehandle); 
   }    
   
 self->n = (self->sample_array->len);
  g_assert_cmpuint (self->n, >, 1);

 {
 const gdouble href0 = self->over_smooth * ncm_stats_dist_nd_vbk_get_rot_bandwidth (dnd, self->d, 0.2 * self->n);
 ncm_vector_set_all (self->href, href0);
 }
  g_ptr_array_set_size (self->invUsample, 0);
  g_array_set_size  (self->norm_array, 0);
  for (i = 0; i < self->n; i++)
  {
  NcmMatrix *cov_i = g_ptr_array_index (self->cov_array, i);
  gdouble norm_i;
  ncm_matrix_cholesky_decomp(g_ptr_array_index(self->cov_array, i), 'U');
  norm_i = ncm_stats_dist_nd_vbk_get_kernel_lnnorm (dnd, cov_i, self->d, self->n, self->href);
  g_array_append_val (self->norm_array, norm_i);
  }

}

static void
_ncm_stats_dist_nd_vbk_prepare (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

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

typedef struct _NcmStatsDistNdVBKEval
{
  NcmStatsDistNdVBK *dnd;
  NcmStatsDistNdVBKPrivate * const self;
  NcmStatsDistNdVBKClass *dnd_class;
} NcmStatsDistNdVBKEval;

static gdouble
_ncm_stats_dist_nd_vbk_prepare_interp_fit_nnls (gdouble os, gpointer userdata)
{
  NcmStatsDistNdVBKEval *eval = userdata;
  const gdouble href0 = os * ncm_stats_dist_nd_vbk_get_rot_bandwidth (eval->dnd, eval->self->d, 0.2 * eval->self->n);
  gdouble rnorm;
  
  ncm_vector_set_all (eval->self->href, href0);
  
  eval->dnd_class->prepare_IM (eval->dnd, eval->self->cov_array, eval->self->d, eval->self->n, eval->self->href, eval->self->IM, eval->self->sample_array, eval->self->norm_array);
  rnorm = ncm_nnls_solve (eval->self->nnls, eval->self->sub_IM, eval->self->sub_x, eval->self->f);

  /*ncm_message ("RNORM: % 22.15g | href %.1e\n", rnorm, href0);*/

  return rnorm;
}

static void
_ncm_stats_dist_nd_vbk_prepare_interp (NcmStatsDistNdVBK *dnd, NcmVector *m2lnp)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;
  ncm_stats_dist_nd_vbk_prepare (dnd);
  g_assert_cmpuint ( ncm_vector_len (m2lnp), ==, self->n);
  {
    NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);
    NcmStatsDistNdVBKEval eval = {dnd, self, dnd_class};
    const gint nrows        = self->n;
    const gint ncols        = self->n;
    const gdouble dbl_limit = 1.0;
    gdouble qq;
    gint i;
    
    g_assert_cmpuint (nrows, >=, ncols);

    /*
     * Preparing allocations
     */
    if (self->n != self->alloc_n)
    {
      ncm_matrix_clear (&self->IM);
      ncm_vector_clear (&self->f);

      self->IM      = ncm_matrix_new (self->n, self->n);
      self->f       = ncm_vector_new (self->n);

      self->alloc_n = self->n;
    }

    if ((self->nnls == NULL) ||
        ((ncm_nnls_get_nrows (self->nnls) != nrows) || (ncm_nnls_get_ncols (self->nnls) != ncols)))
    {
      ncm_nnls_clear (&self->nnls);
      self->nnls = ncm_nnls_new (nrows, ncols);
      ncm_nnls_set_umethod (self->nnls, NCM_NNLS_UMETHOD_NORMAL);

      ncm_matrix_clear (&self->sub_IM);
      ncm_vector_clear (&self->sub_x);

      self->sub_IM = ncm_matrix_get_submatrix (self->IM, 0, 0, nrows, ncols);
      self->sub_x  = ncm_vector_get_subvector (self->weights, 0, ncols);
    }
    ncm_vector_set_zero (self->weights);
    /*
     * Evaluating the right-hand-side
     */

    self->min_m2lnp = GSL_POSINF;
    self->max_m2lnp = GSL_NEGINF;

    ncm_stats_vec_reset (self->m2lnp_stats, TRUE);
    for (i = 0; i < self->n; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);

      ncm_stats_vec_set (self->m2lnp_stats, 0, m2lnp_i);
      ncm_stats_vec_update (self->m2lnp_stats);
      self->min_m2lnp = MIN (self->min_m2lnp, m2lnp_i);
      self->max_m2lnp = MAX (self->max_m2lnp, m2lnp_i);
    }
    qq = ncm_stats_vec_get_quantile (self->m2lnp_stats, 0);

    /*printf ("% 22.15g % 22.15g % 22.15g\n", self->max_m2lnp, self->min_m2lnp, qq);*/

    if (-0.5 * (self->max_m2lnp - self->min_m2lnp) < dbl_limit * GSL_LOG_DBL_EPSILON)
    {
      const guint fi = ncm_vector_get_min_index (m2lnp);

      g_assert_cmpuint (fi, <, ncm_vector_len (m2lnp));

      ncm_vector_set_zero (self->weights);
      ncm_vector_set (self->weights, fi, 1.0);
      ncm_vector_set_all (self->href, 1.0);

      return;
    }

    for (i = 0; i < self->n; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);

      if (m2lnp_i > qq)
        ncm_vector_set (self->f, i, exp (-0.5 * (qq      - self->min_m2lnp)));
      else
        ncm_vector_set (self->f, i, exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
    }

    if (self->n > 10000)
      g_warning ("_ncm_stats_dist_nd_vbk_prepare_interp: very large system n = %u!", self->n);
    
    self->rnorm = _ncm_stats_dist_nd_vbk_prepare_interp_fit_nnls (self->over_smooth, &eval);
    
  }
}

static void
_ncm_stats_dist_nd_vbk_reset (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  g_ptr_array_set_size (self->sample_array, 0);
  g_ptr_array_set_size(self->cov_array, 0);
  g_array_set_size (self->norm_array, 0);
}

/**
 * ncm_stats_dist_nd_vbk_ref:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Increases the reference count of @dnd.
 * 
 * Returns: (transfer full): @dnd.
 */
NcmStatsDistNdVBK *
ncm_stats_dist_nd_vbk_ref (NcmStatsDistNdVBK *dnd)
{
  return g_object_ref (dnd);
}

/**
 * ncm_stats_dist_nd_vbk_free:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Decreases the reference count of @dnd.
 *
 */
void 
ncm_stats_dist_nd_vbk_free (NcmStatsDistNdVBK *dnd)
{
  g_object_unref (dnd);
}

/**
 * ncm_stats_dist_nd_vbk_clear:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Decreases the reference count of *@dnd and sets the pointer *@dnd to NULL.
 *
 */
void 
ncm_stats_dist_nd_vbk_clear (NcmStatsDistNdVBK **dnd)
{
  g_clear_object (dnd);
}

/**
 * ncm_stats_dist_nd_vbk_get_dim:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Returns: the dimension of the sample space.
 */
guint 
ncm_stats_dist_nd_vbk_get_dim (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return self->d;
}

/**
 * ncm_stats_dist_nd_vbk_get_rot_bandwidth: (virtual get_rot_bandwidth)
 * @dnd: a #NcmStatsDistNdVBK
 * @d: problem dimension
 * @n: sample size
 *
 * Returns: the rule-of-thumb bandwidth estimate.
 */
gdouble
ncm_stats_dist_nd_vbk_get_rot_bandwidth (NcmStatsDistNdVBK *dnd, const guint d, const gdouble n)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);

  return dnd_class->get_rot_bandwidth (dnd, d, n);
}

/**
 * ncm_stats_dist_nd_vbk_get_kernel_lnnorm: (virtual get_kernel_lnnorm)
 * @dnd: a #NcmStatsDistNdVBK
 * @cov_decomp:GPtrArray with Cholesky decomposition (U) of each covariance matrix
 * @d: problem dimension
 * @n: sample size
 * @href: interpolation bandwidth
 *
 * Returns: the log-normalization of a single kernel.
 */
gdouble
ncm_stats_dist_nd_vbk_get_kernel_lnnorm (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return dnd_class->get_kernel_lnnorm (dnd, cov_decomp, d, n, href);
}

/**
 * ncm_stats_dist_nd_vbk_set_over_smooth:
 * @dnd: a #NcmStatsDistNdVBK
 * @over_smooth: the over-smooth factor
 *
 * Sets the over-smooth factor to @over_smooth.
 *
 */
void
ncm_stats_dist_nd_vbk_set_over_smooth (NcmStatsDistNdVBK *dnd, const gdouble over_smooth)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  self->over_smooth = over_smooth;
}

/**
 * ncm_stats_dist_nd_vbk_get_over_smooth:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Returns: the over-smooth factor.
 */
gdouble
ncm_stats_dist_nd_vbk_get_over_smooth (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return self->over_smooth;
}

/**
 * ncm_stats_dist_nd_vbk_set_split_frac:
 * @dnd: a #NcmStatsDistNdVBK
 * @split_frac: the over-smooth factor
 *
 * Sets cross-correlation split fraction to @split_frac.
 *
 */
void
ncm_stats_dist_nd_vbk_set_split_frac (NcmStatsDistNdVBK *dnd, const gdouble split_frac)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  g_assert_cmpfloat (split_frac, >=, 0.5);
  g_assert_cmpfloat (split_frac, <=, 1.0);

  self->split_frac = split_frac;
}

/**
 * ncm_stats_dist_nd_vbk_get_split_frac:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Returns: the cross-correlation split fraction.
 */
gdouble
ncm_stats_dist_nd_vbk_get_split_frac (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return self->split_frac;
}

/**
 * ncm_stats_dist_nd_vbk_set_nearPD_maxiter:
 * @dnd: a #NcmStatsDistNdVBK
 * @maxiter: maximum number of iterations
 *
 * Sets the maximum number of iterations when finding the
 * nearest positive definite covariance matrix to @maxiter.
 *
 */
void
ncm_stats_dist_nd_vbk_set_nearPD_maxiter (NcmStatsDistNdVBK *dnd, const guint maxiter)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  self->nearPD_maxiter = maxiter;
}

/**
 * ncm_stats_dist_nd_vbk_get_nearPD_maxiter:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Returns: maximum number of iterations when finding the nearest positive definite covariance matrix.
 */
guint
ncm_stats_dist_nd_vbk_get_nearPD_maxiter (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return self->nearPD_maxiter;
}

/**
 * ncm_stats_dist_nd_vbk_set_cv_type:
 * @dnd: a #NcmStatsDistNdVBK
 * @cv_type: a #NcmStatsDistNdVBKCV
 *
 * Sets the cross-validation method to @cv_type.
 *
 */
void
ncm_stats_dist_nd_vbk_set_cv_type (NcmStatsDistNdVBK *dnd, const NcmStatsDistNdVBKCV cv_type)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  self->cv_type = cv_type;
}

/**
 * ncm_stats_dist_nd_vbk_get_cv_type:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Returns: current cross-validation method used.
 */
NcmStatsDistNdVBKCV
ncm_stats_dist_nd_vbk_get_cv_type (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return self->cv_type;
}

/**
 * ncm_stats_dist_nd_vbk_prepare: (virtual prepare)
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Prepares the object for calculations.
 */
void 
ncm_stats_dist_nd_vbk_prepare (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd); 

  if (dnd_class->prepare != NULL)
    dnd_class->prepare (dnd);
}

/**
 * ncm_stats_dist_nd_vbk_prepare_interp: (virtual prepare_interp)
 * @dnd: a #NcmStatsDistNdVBK
 * @m2lnp: a #NcmVector containing the distribution values
 *
 * Prepares the object for calculations. Using the distribution values
 * at the sample points.
 * 
 */
void 
ncm_stats_dist_nd_vbk_prepare_interp (NcmStatsDistNdVBK *dnd, NcmVector *m2lnp)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd); 

  g_assert (dnd_class->prepare_interp != NULL);
  dnd_class->prepare_interp (dnd, m2lnp);
}

/**
 * ncm_stats_dist_nd_vbk_eval:
 * @dnd: a #NcmStatsDistNdVBK
 * @x: a #NcmVector
 *
 * Evaluate the distribution at $\vec{x}=$@x. If the distribution
 * was prepared using ncm_stats_dist_nd_vbk_prepare_interp(), the 
 * results will follow the interpolation and may not be properly 
 * normalized. In this case the method ncm_stats_dist_nd_vbk_eval_m2lnp()
 * should be used to avoid underflow.
 * 
 * Returns: $P(\vec{x})$.
 */
gdouble 
ncm_stats_dist_nd_vbk_eval (NcmStatsDistNdVBK *dnd, NcmVector *x)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd); 
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return dnd_class->eval (dnd, self->weights, self->v, self->sample_array, self->d, self->n, self->href, self->cov_array, self->norm_array);
}

/**
 * ncm_stats_dist_nd_vbk_eval_m2lnp:
 * @dnd: a #NcmStatsDistNdVBK
 * @x: a #NcmVector
 *
 * Evaluate the distribution at $\vec{x}=$@x. If the distribution
 * was prepared using ncm_stats_dist_nd_vbk_prepare_interp(), the 
 * results will follow the interpolation and may not be properly 
 * normalized.
 * 
 * Returns: $P(\vec{x})$.
 */
gdouble 
ncm_stats_dist_nd_vbk_eval_m2lnp (NcmStatsDistNdVBK *dnd, NcmVector *x)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  ncm_vector_memcpy (self->v, x);
  
  if (dnd_class->eval_m2lnp != NULL)
    return dnd_class->eval_m2lnp (dnd, self->weights, self->v, self->invUsample, self->d, self->n, self->href, self->cov_array) + 2.0 * self->kernel_lnnorm;
  else
    return -2.0 * (log (dnd_class->eval (dnd, self->weights, self->v, self->sample_array, self->d, self->n, self->href, self->cov_array, self->norm_array))); /*- self->kernel_lnnorm);*/
}

/**
 * ncm_stats_dist_nd_vbk_sample:
 * @dnd: a #NcmStatsDistNdVBK
 * @x: a #NcmVector
 * @rng: a #NcmRNG
 * 
 * Using the pseudo-random number generator @rng generates a 
 * point from the distribution and copy it to @x.
 * 
 */
void
ncm_stats_dist_nd_vbk_sample (NcmStatsDistNdVBK *dnd, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;
  const guint n = ncm_vector_len (self->weights);
  gint i;

  g_array_set_size (self->sampling, ncm_vector_len (self->weights));
  gsl_ran_multinomial (rng->r, n, 1, ncm_vector_data (self->weights), (guint *)self->sampling->data);

  for (i = 0; i < n; i++)
  {
    if (g_array_index (self->sampling, guint, i) > 0)
    {
      NcmVector *y_i = g_ptr_array_index (self->sample_array, i);
      NcmMatrix *cov_decomp = g_ptr_array_index (self->cov_array,i);
      ncm_stats_dist_nd_vbk_kernel_sample (dnd, cov_decomp, x, y_i, self->href, rng);
      break;
    }
  }
}

/**
 * ncm_stats_dist_nd_vbk_get_rnorm:
 * @dnd: a #NcmStatsDistNdVBK
 *
 * Gets the value of the last $\chi^2$ fit obtained
 * when computing the interpolation through
 * ncm_stats_dist_nd_vbk_prepare_interp().
 *
 * Returns: the value of the $\chi^2$.
 */
gdouble
ncm_stats_dist_nd_vbk_get_rnorm (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;
  return self->rnorm * self->rnorm;
}

/**
 * ncm_stats_dist_nd_vbk_kernel_sample:
 * @dnd: a #NcmStatsDistNdVBK
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
ncm_stats_dist_nd_vbk_kernel_sample (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, NcmVector *x, NcmVector *mu, const NcmVector *href, NcmRNG *rng)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd); 
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  dnd_class->kernel_sample (dnd, cov_decomp, self->d, x, mu, href, rng);
}

/**
 * ncm_stats_dist_nd_vbk_kernel_eval_m2lnp: (virtual kernel_eval_m2lnp)
 * @dnd: a #NcmStatsDistNdVBK
 * @x: a #NcmVector
 * @y: a #NcmVector
 * @href: covariance href
 * 
 * Evaluates a single kernel at @x and @y and href @s, i.e., $K_s(x,y)$.
 * 
 * Returns: $K_s(x,y)$.
 */
gdouble
ncm_stats_dist_nd_vbk_kernel_eval_m2lnp (NcmStatsDistNdVBK *dnd, NcmVector *x, NcmVector *y, const NcmVector *href)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);
  NcmStatsDistNdVBKPrivate * const self = dnd->priv;

  return dnd_class->kernel_eval_m2lnp (dnd, self->cov_decomp, self->d, x, y, self->v, href) + 2.0 * self->kernel_lnnorm;
}

/**
 * ncm_stats_dist_nd_vbk_add_obs_weight:
 * @dndg: a #NcmStatsDistNdVBK
 * @y: a #NcmVector
 * @w: weight
 *
 * Adds a new point @y to the sample with weight @w.
 *
 */
void
ncm_stats_dist_nd_vbk_add_obs_weight (NcmStatsDistNdVBK *dndg, NcmVector *y, const gdouble w)
{
  NcmStatsDistNdVBKPrivate * const self = dndg->priv;

  g_ptr_array_add (self->sample_array, ncm_vector_dup (y));
}

/**
 * ncm_stats_dist_nd_vbk_add_obs:
 * @dndg: a #NcmStatsDistNdVBK
 * @y: a #NcmVector
 *
 * Adds a new point @y to the sample with weight 1.0.
 *
 */
void
ncm_stats_dist_nd_vbk_add_obs (NcmStatsDistNdVBK *dndg, NcmVector *y)
{
  ncm_stats_dist_nd_vbk_add_obs_weight (dndg, y, 1.0);
}

/**
 * ncm_stats_dist_nd_vbk_reset: (virtual reset)
 * @dnd: a #NcmStatsDistNdVBK
 * 
 * Reset the object discarding all added points.
 * 
 */
void 
ncm_stats_dist_nd_vbk_reset (NcmStatsDistNdVBK *dnd)
{
  NcmStatsDistNdVBKClass *dnd_class = NCM_STATS_DIST_ND_VBK_GET_CLASS (dnd);
  dnd_class->reset (dnd);
}

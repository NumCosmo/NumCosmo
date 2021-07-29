/***************************************************************************
 *            ncm_fit_esmcmc_walker_apes.c
 *
 *  Sat October 27 13:08:13 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_apes.c
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
 * SECTION:ncm_fit_esmcmc_walker_apes
 * @title: NcmFitESMCMCWalkerAPES
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker - aps move.
 *
 * Implementing aps move walker for #NcmFitESMCMC (affine invariant).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_apes.h"

#include "math/ncm_fit_esmcmc.h"
#include "math/ncm_stats_dist_vkde.h"
#include "math/ncm_stats_dist_kernel_st.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sort.h>
#endif /* NUMCOSMO_GIR_SCAN */


enum
{
  PROP_0,
  PROP_USE_INTERP,
};

struct _NcmFitESMCMCWalkerAPESPrivate
{
  guint size;
  guint size_2;
  guint nparams;
  NcmVector *m2lnp_star;
  NcmVector *m2lnp_cur;
  gchar *desc;
  NcmStatsDist *sd0;
  NcmStatsDist *sd1;
  GPtrArray *thetastar;
  NcmVector *m2lnL_s0;
  NcmVector *m2lnL_s1;
  gboolean use_interp;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFitESMCMCWalkerAPES, ncm_fit_esmcmc_walker_apes, NCM_TYPE_FIT_ESMCMC_WALKER);

static void
ncm_fit_esmcmc_walker_apes_init (NcmFitESMCMCWalkerAPES *aps)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv = ncm_fit_esmcmc_walker_apes_get_instance_private (aps);
  
  self->size       = 0;
  self->size_2     = 0;
  self->nparams    = 0;
  self->m2lnp_star = NULL;
  self->m2lnp_cur  = NULL;
  self->desc       = "APES-Move";
  self->sd0        = NULL;
  self->sd1        = NULL;
  self->thetastar  = g_ptr_array_new ();
  self->m2lnL_s0   = NULL;
  self->m2lnL_s1   = NULL;
  self->use_interp = FALSE;
  
  g_ptr_array_set_free_func (self->thetastar, (GDestroyNotify) ncm_vector_free);
}

static void
_ncm_fit_esmcmc_walker_apes_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerAPES *aps = NCM_FIT_ESMCMC_WALKER_APES (object);
  
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APES (object));
  
  switch (prop_id)
  {
    case PROP_USE_INTERP:
      ncm_fit_esmcmc_walker_apes_use_interp (aps, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_apes_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (object);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APES (object));
  
  switch (prop_id)
  {
    case PROP_USE_INTERP:
      g_value_set_boolean (value, self->use_interp);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_apes_dispose (GObject *object)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (object);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  ncm_vector_clear (&self->m2lnp_star);
  ncm_vector_clear (&self->m2lnp_cur);
  ncm_vector_clear (&self->m2lnL_s0);
  ncm_vector_clear (&self->m2lnL_s1);
  
  ncm_stats_dist_clear (&self->sd0);
  ncm_stats_dist_clear (&self->sd1);
  
  g_clear_pointer (&self->thetastar, g_ptr_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_apes_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_apes_finalize (GObject *object)
{
  /*NcmFitESMCMCWalkerAPES *aps = NCM_FIT_ESMCMC_WALKER_APES (object);*/
  /*NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_apes_parent_class)->finalize (object);
}

static void _ncm_fit_esmcmc_walker_apes_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_apes_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_apes_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_apes_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_apes_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_apes_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_apes_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_apes_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_apes_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_apes_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_apes_class_init (NcmFitESMCMCWalkerAPESClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);
  
  object_class->set_property = _ncm_fit_esmcmc_walker_apes_set_property;
  object_class->get_property = _ncm_fit_esmcmc_walker_apes_get_property;
  object_class->dispose      = _ncm_fit_esmcmc_walker_apes_dispose;
  object_class->finalize     = _ncm_fit_esmcmc_walker_apes_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_USE_INTERP,
                                   g_param_spec_boolean ("use-interp",
                                                         NULL,
                                                         "Whether to use interpolation to build the posterior approximation",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  walker_class->set_size    = &_ncm_fit_esmcmc_walker_apes_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_apes_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_apes_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_apes_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_apes_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_apes_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_apes_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_apes_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_apes_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_apes_desc;
}

static void
_ncm_fit_esmcmc_walker_apes_set_sys (NcmFitESMCMCWalker *walker, guint size, guint nparams)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (nparams, >, 0);
  
  if ((self->size != size) || (self->nparams != nparams))
  {
    gint i;
    
    ncm_stats_dist_clear (&self->sd0);
    ncm_stats_dist_clear (&self->sd1);
    ncm_vector_clear (&self->m2lnp_star);
    ncm_vector_clear (&self->m2lnp_cur);
    
    g_ptr_array_set_size (self->thetastar, 0);
    
    g_assert (size % 2 == 0);
    self->size    = size;
    self->size_2  = size / 2;
    self->nparams = nparams;
    
    self->m2lnp_star = ncm_vector_new (self->size);
    self->m2lnp_cur  = ncm_vector_new (self->size);
    self->m2lnL_s0   = ncm_vector_new (self->size_2);
    self->m2lnL_s1   = ncm_vector_new (self->size_2);
    
    if (TRUE)
    {
      NcmStatsDistKernel *sdk = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (self->nparams, 1.0));
      
      self->sd0 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (sdk, NCM_STATS_DIST_CV_NONE));
      self->sd1 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (sdk, NCM_STATS_DIST_CV_NONE));
    }
    
    ncm_stats_dist_set_over_smooth (self->sd0, 1.0);
    ncm_stats_dist_set_over_smooth (self->sd1, 1.0);
    ncm_stats_dist_set_split_frac (self->sd0, 1.0);
    ncm_stats_dist_set_split_frac (self->sd1, 1.0);
    
    for (i = 0; i < self->size; i++)
    {
      NcmVector *thetastar_i = ncm_vector_new (self->nparams);
      
      g_ptr_array_add (self->thetastar, thetastar_i);
    }
  }
}

static void
_ncm_fit_esmcmc_walker_apes_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  g_assert_cmpuint (size, >, 0);
  
  if (self->nparams != 0)
    _ncm_fit_esmcmc_walker_apes_set_sys (walker, size, self->nparams);
  else
    self->size = size;
}

static guint
_ncm_fit_esmcmc_walker_apes_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  return self->size;
}

static void
_ncm_fit_esmcmc_walker_apes_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != 0)
    _ncm_fit_esmcmc_walker_apes_set_sys (walker, self->size, nparams);
  else
    self->nparams = nparams;
}

static guint
_ncm_fit_esmcmc_walker_apes_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  return self->nparams;
}

static void
_ncm_fit_esmcmc_walker_apes_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  const gdouble T                            = 1.0;
  gint i;
  
  if (ki < self->size_2)
  {
    ncm_stats_dist_reset (self->sd0);
    
    for (i = self->size_2; i < self->size; i++)
    {
      /*NcmVector *theta_i = g_ptr_array_index (theta, i);*/
      ncm_vector_set (self->m2lnL_s0, i - self->size_2, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      /*ncm_stats_dist_add_obs (NCM_STATS_DIST_ND_VBK (self->dndg0), theta_i);*/
    }
    
    {
      const gint n = self->size - self->size_2;
      size_t p[n];
      
      gsl_sort_index (p, ncm_vector_data (self->m2lnL_s0), ncm_vector_stride (self->m2lnL_s0), n);
      
      for (i = self->size_2; i < self->size; i++)
      {
        const gint j       = self->size_2 + p[(i - self->size_2)];
        NcmVector *theta_j = g_ptr_array_index (theta, j);
        
        ncm_vector_set (self->m2lnL_s0, i - self->size_2, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));
        ncm_stats_dist_add_obs (self->sd0, theta_j);
        /*printf ("AddingA [%d %d %d] % 22.15g\n", i, n - 1 - (i-self->size_2), j, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));*/
      }
    }
    
    if (self->use_interp)
      ncm_stats_dist_prepare_interp (self->sd0, self->m2lnL_s0);
    else
      ncm_stats_dist_prepare (self->sd0);
    
    for (i = ki; i < self->size_2; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);
      
      /*do {*/
      ncm_stats_dist_sample (self->sd0, thetastar_i, rng);
      /*} while (!ncm_mset_fparam_validate_all (mset, thetastar_i));*/
    }
  }
  
  if (kf >= self->size_2)
  {
    ncm_stats_dist_reset (self->sd1);
    
    for (i = 0; i < self->size_2; i++)
    {
      /*NcmVector *theta_i = g_ptr_array_index (theta, i);*/
      ncm_vector_set (self->m2lnL_s1, i, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      /*ncm_stats_dist_add_obs (self->dndg1, theta_i);*/
    }
    
    {
      const gint n = self->size_2;
      size_t p[n];
      
      gsl_sort_index (p, ncm_vector_data (self->m2lnL_s1), ncm_vector_stride (self->m2lnL_s1), n);
      
      for (i = 0; i < self->size_2; i++)
      {
        const gint j       = p[i];
        NcmVector *theta_j = g_ptr_array_index (theta, j);
        
        ncm_vector_set (self->m2lnL_s1, i, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));
        ncm_stats_dist_add_obs (self->sd1, theta_j);
        /*printf ("AddingB [%d %d %d] % 22.15g\n", i, n - 1 - i, j, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));*/
      }
    }
    
    
    if (self->use_interp)
      ncm_stats_dist_prepare_interp (self->sd1, self->m2lnL_s1);
    else
      ncm_stats_dist_prepare (self->sd1);
    
    for (i = self->size_2; i < kf; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);
      
      /*do {*/
      ncm_stats_dist_sample (self->sd1, thetastar_i, rng);
      /*} while (!ncm_mset_fparam_validate_all (mset, thetastar_i));*/
    }
  }
}

static void
_ncm_fit_esmcmc_walker_apes_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  NcmVector *theta_k                         = g_ptr_array_index (theta, k);
  
  ncm_vector_memcpy (thetastar, g_ptr_array_index (self->thetastar, k));
  
  if (k < self->size_2)
  {
    const gdouble m2lnaps_star = ncm_stats_dist_eval_m2lnp (self->sd0, thetastar);
    const gdouble m2lnaps_cur  = ncm_stats_dist_eval_m2lnp (self->sd0, theta_k);
    
    ncm_vector_set (self->m2lnp_star, k, m2lnaps_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnaps_cur);
    /*printf ("%u % 22.15g | lnp_cur % 22.15g lnp_star % 22.15g\n", k, - 0.5 * (m2lnaps_cur - m2lnaps_star), - 0.5 * m2lnaps_cur, - 0.5 * m2lnaps_star);*/
  }
  
  if (k >= self->size_2)
  {
    const gdouble m2lnaps_star = ncm_stats_dist_eval_m2lnp (self->sd1, thetastar);
    const gdouble m2lnaps_cur  = ncm_stats_dist_eval_m2lnp (self->sd1, theta_k);
    
    ncm_vector_set (self->m2lnp_star, k, m2lnaps_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnaps_cur);
    /*printf ("%u % 22.15g | lnp_cur % 22.15g lnp_star % 22.15g\n", k, - 0.5 * (m2lnaps_cur - m2lnaps_star), - 0.5 * m2lnaps_cur, - 0.5 * m2lnaps_star);*/
  }
  
  /*ncm_vector_log_vals (theta_k,   "    THETA: ", "% 22.15g", TRUE);*/
  /*ncm_vector_log_vals (thetastar, "THETASTAR: ", "% 22.15g", TRUE);*/
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  const gdouble m2lnp_star                   = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur                    = ncm_vector_get (self->m2lnp_cur, k);
  
/*
 *  printf ("AAA m2lnL_star % 22.15g m2lnp_star % 22.15g m2lnL_cur % 22.15g m2lnp_cur % 22.15g L cur->star: %12.5g p star->cur: %12.5g | T %12.5g\n",
 *         m2lnL_star, m2lnp_star, m2lnL_cur, m2lnp_cur,
 *         exp (- 0.5 * (m2lnL_star - m2lnL_cur)), exp (- 0.5 * (m2lnp_cur - m2lnp_star)),
 *         MIN (exp (- 0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur))), 1.0));
 */
  
  return exp (-0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur)));
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  const gdouble m2lnp_star                   = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur                    = ncm_vector_get (self->m2lnp_cur, k);
  
  return -0.5 * (m2lnp_cur - m2lnp_star);
}

static void
_ncm_fit_esmcmc_walker_apes_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /*NcmFitESMCMCWalkerAPES *aps = NCM_FIT_ESMCMC_WALKER_APES (walker);*/
  /*NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;*/
  
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_apes_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *aps                = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  return self->desc;
}

/**
 * ncm_fit_esmcmc_walker_apes_new:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 *
 * Creates a new #NcmFitESMCMCWalkerAPES to be used
 * with @nwalkers.
 *
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerAPES.
 */
NcmFitESMCMCWalkerAPES *
ncm_fit_esmcmc_walker_apes_new (guint nwalkers, guint nparams)
{
  NcmFitESMCMCWalkerAPES *aps = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_APES,
                                              "size", nwalkers,
                                              "nparams", nparams,
                                              NULL);
  
  return aps;
}

/**
 * ncm_fit_esmcmc_walker_apes_use_interp:
 * @aps: a #NcmFitESMCMCWalkerAPES
 * @use_interp: whether to use interpolation of the posterior
 *
 * Sets whether to use interpolation of the posterior approximation (@use_interp == TRUE)
 * or kernel density estimate (@use_interp == FALSE).
 *
 */
void
ncm_fit_esmcmc_walker_apes_use_interp (NcmFitESMCMCWalkerAPES *aps, gboolean use_interp)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = aps->priv;
  
  self->use_interp = use_interp;
}


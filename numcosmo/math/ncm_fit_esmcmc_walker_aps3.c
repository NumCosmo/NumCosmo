/***************************************************************************
 *            ncm_fit_esmcmc_walker_aps.c
 *
 *  Sat October 27 13:08:13 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_aps.c
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
 * SECTION:ncm_fit_esmcmc_walker_aps
 * @title: NcmFitESMCMCWalkerAPS
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
#include "math/ncm_fit_esmcmc_walker_aps.h"

#include "math/ncm_fit_esmcmc.h"
#include "math/ncm_stats_dist_nd_vbk_studentt.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sort.h>
#endif /* NUMCOSMO_GIR_SCAN */


enum
{
  PROP_0,
  PROP_USE_INTERP,
};

struct _NcmFitESMCMCWalkerAPSPrivate
{
  guint size;
  guint size_2;
  guint nparams;
  NcmVector *m2lnp_star;
  NcmVector *m2lnp_cur;
  gchar *desc;
  NcmStatsDistNdVbk *dndg0;
  NcmStatsDistNdVbk *dndg1;
  GPtrArray *thetastar;
  NcmVector *m2lnL_s0;
  NcmVector *m2lnL_s1;
  gboolean use_interp;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFitESMCMCWalkerAPS, ncm_fit_esmcmc_walker_aps, NCM_TYPE_FIT_ESMCMC_WALKER);

static void
ncm_fit_esmcmc_walker_aps_init (NcmFitESMCMCWalkerAPS *aps)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv = ncm_fit_esmcmc_walker_aps_get_instance_private (aps);

  self->size            = 0;
  self->size_2          = 0;
  self->nparams         = 0;
  self->m2lnp_star      = NULL;
  self->m2lnp_cur       = NULL;
  self->desc            = "APS-Move";
  self->dndg0           = NULL;
  self->dndg1           = NULL;
  self->thetastar       = g_ptr_array_new ();
  self->m2lnL_s0        = NULL;
  self->m2lnL_s1        = NULL;
  self->use_interp      = FALSE;

  g_ptr_array_set_free_func (self->thetastar, (GDestroyNotify) ncm_vector_free);
}

static void
_ncm_fit_esmcmc_walker_aps_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APS (object));

  switch (prop_id)
  {
    case PROP_USE_INTERP:
      ncm_fit_esmcmc_walker_aps_use_interp (aps, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_aps_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APS (object));

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
_ncm_fit_esmcmc_walker_aps_dispose (GObject *object)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  ncm_vector_clear (&self->m2lnp_star);
  ncm_vector_clear (&self->m2lnp_cur);
  ncm_vector_clear (&self->m2lnL_s0);
  ncm_vector_clear (&self->m2lnL_s1);
  
  ncm_stats_dist_nd_vbk_clear (&self->dndg0);
  ncm_stats_dist_nd_vbk_clear (&self->dndg1);

  g_clear_pointer (&self->thetastar, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_aps_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_aps_finalize (GObject *object)
{
  /*NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (object);*/
  /*NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_aps_parent_class)->finalize (object);
}

static void _ncm_fit_esmcmc_walker_aps_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_aps_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_aps_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_aps_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_aps_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_aps_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_aps_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_aps_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_aps_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_aps_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_aps_class_init (NcmFitESMCMCWalkerAPSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = _ncm_fit_esmcmc_walker_aps_set_property;
  object_class->get_property = _ncm_fit_esmcmc_walker_aps_get_property;
  object_class->dispose      = _ncm_fit_esmcmc_walker_aps_dispose;
  object_class->finalize     = _ncm_fit_esmcmc_walker_aps_finalize;

  g_object_class_install_property (object_class,
                                   PROP_USE_INTERP,
                                   g_param_spec_boolean ("use-interp",
                                                         NULL,
                                                         "Whether to use interpolation to build the posterior approximation",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  walker_class->set_size    = &_ncm_fit_esmcmc_walker_aps_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_aps_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_aps_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_aps_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_aps_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_aps_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_aps_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_aps_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_aps_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_aps_desc;
}

static void 
_ncm_fit_esmcmc_walker_aps_set_sys (NcmFitESMCMCWalker *walker, guint size, guint nparams)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != size || self->nparams != nparams)
  {
    gint i;
    
    ncm_stats_dist_nd_vbk_clear (&self->dndg0);
    ncm_stats_dist_nd_vbk_clear (&self->dndg1);
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
      self->dndg0 = NCM_STATS_DIST_ND_VBK (ncm_stats_dist_nd_vbk_studentt_new (self->nparams, NCM_STATS_DIST_ND_VBK_CV_NONE, 1.0));
      self->dndg1 = NCM_STATS_DIST_ND_VBK (ncm_stats_dist_nd_vbk_studentt_new (self->nparams, NCM_STATS_DIST_ND_VBK_CV_NONE, 1.0));
    }
 /*   else
    {
      self->dndg0 = NCM_STATS_DIST_ND_VBK (ncm_stats_dist_nd_vbk_gauss_new (self->nparams, NCM_STATS_DIST_ND_CV_SPLIT));
      self->dndg1 = NCM_STATS_DIST_ND_VBK (ncm_stats_dist_nd_vbk_gauss_new (self->nparams, NCM_STATS_DIST_ND_CV_SPLIT));
    }*/

    ncm_stats_dist_nd_vbk_set_over_smooth (self->dndg0, 1.0);
    ncm_stats_dist_nd_vbk_set_over_smooth (self->dndg1, 1.0);
    ncm_stats_dist_nd_vbk_set_split_frac (self->dndg0, 1.0);
    ncm_stats_dist_nd_vbk_set_split_frac (self->dndg1, 1.0);
    
    for (i = 0; i < self->size; i++)
    {
      NcmVector *thetastar_i = ncm_vector_new (self->nparams);
      g_ptr_array_add (self->thetastar, thetastar_i);
    }
  }
}

static void 
_ncm_fit_esmcmc_walker_aps_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  g_assert_cmpuint (size, >, 0);
  
  if (self->nparams != 0)
    _ncm_fit_esmcmc_walker_aps_set_sys (walker, size, self->nparams);
  else
    self->size = size;  
}

static guint 
_ncm_fit_esmcmc_walker_aps_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  return self->size;
}

static void 
_ncm_fit_esmcmc_walker_aps_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  g_assert_cmpuint (nparams, >, 0);
  
  if (self->size != 0)
    _ncm_fit_esmcmc_walker_aps_set_sys (walker, self->size, nparams);
  else
    self->nparams = nparams;
}

static guint 
_ncm_fit_esmcmc_walker_aps_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;

  return self->nparams;
}

static void 
_ncm_fit_esmcmc_walker_aps_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  const gdouble T = 1.0;
  gint i;

  if (ki < self->size_2)
  {
    ncm_stats_dist_nd_vbk_reset (NCM_STATS_DIST_ND_VBK (self->dndg0));
    for (i = self->size_2; i < self->size; i++)
    {
      /*NcmVector *theta_i = g_ptr_array_index (theta, i);*/
      ncm_vector_set (self->m2lnL_s0, i - self->size_2, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      /*ncm_stats_dist_nd_vbk_add_obs (NCM_STATS_DIST_ND_VBK (self->dndg0), theta_i);*/
    }

    {
      const gint n = self->size - self->size_2;
      size_t p[n];
      gsl_sort_index (p, ncm_vector_data (self->m2lnL_s0), ncm_vector_stride (self->m2lnL_s0), n);

      for (i = self->size_2; i < self->size; i++)
      {
        const gint j = self->size_2 + p[n - 1 - (i-self->size_2)];
        NcmVector *theta_j = g_ptr_array_index (theta, j);
        ncm_vector_set (self->m2lnL_s0, i - self->size_2, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));
        ncm_stats_dist_nd_vbk_add_obs (self->dndg0, theta_j);
        /*printf ("AddingA [%d %d %d] % 22.15g\n", i, n - 1 - (i-self->size_2), j, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));*/
      }
    }

    if (self->use_interp)
      ncm_stats_dist_nd_vbk_prepare_interp (NCM_STATS_DIST_ND_VBK (self->dndg0), self->m2lnL_s0);
    else
      ncm_stats_dist_nd_vbk_prepare (NCM_STATS_DIST_ND_VBK (self->dndg0));

    for (i = ki; i < self->size_2; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);
      do {
        ncm_stats_dist_nd_vbk_sample (NCM_STATS_DIST_ND_VBK (self->dndg0), thetastar_i, rng);
      } while (!ncm_mset_fparam_validate_all (mset, thetastar_i));
    }
  }
  if (kf >= self->size_2)
  {
    ncm_stats_dist_nd_vbk_reset (NCM_STATS_DIST_ND_VBK (self->dndg1));
    
    for (i = 0; i < self->size_2; i++)
    {
      /*NcmVector *theta_i = g_ptr_array_index (theta, i);*/
      ncm_vector_set (self->m2lnL_s1, i, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      /*ncm_stats_dist_nd_vbk_add_obs (self->dndg1, theta_i);*/
    }

    {
      const gint n = self->size_2;
      size_t p[n];
      gsl_sort_index (p, ncm_vector_data (self->m2lnL_s1), ncm_vector_stride (self->m2lnL_s1), n);

      for (i = 0; i < self->size_2; i++)
      {
        const gint j = p[n - 1 - i];
        NcmVector *theta_j = g_ptr_array_index (theta, j);
        ncm_vector_set (self->m2lnL_s1, i, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));
        ncm_stats_dist_nd_vbk_add_obs (self->dndg1, theta_j);
        /*printf ("AddingB [%d %d %d] % 22.15g\n", i, n - 1 - i, j, ncm_vector_get (g_ptr_array_index (m2lnL, j), 0));*/
      }
    }

    
    if (self->use_interp)
      ncm_stats_dist_nd_vbk_prepare_interp (NCM_STATS_DIST_ND_VBK (self->dndg1), self->m2lnL_s1);
    else
      ncm_stats_dist_nd_vbk_prepare (NCM_STATS_DIST_ND_VBK (self->dndg1));
    
    for (i = self->size_2; i < kf; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);
      do {
        ncm_stats_dist_nd_vbk_sample (NCM_STATS_DIST_ND_VBK (self->dndg1), thetastar_i, rng);
      } while (!ncm_mset_fparam_validate_all (mset, thetastar_i));
    }
  }
}
static void
_ncm_fit_esmcmc_walker_aps_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  NcmVector *theta_k = g_ptr_array_index (theta, k);

  ncm_vector_memcpy (thetastar, g_ptr_array_index (self->thetastar, k));

  if (k < self->size_2)
  {
    const gdouble m2lnaps_star = ncm_stats_dist_nd_vbk_eval_m2lnp (NCM_STATS_DIST_ND_VBK (self->dndg0), thetastar);
    const gdouble m2lnaps_cur  = ncm_stats_dist_nd_vbk_eval_m2lnp (NCM_STATS_DIST_ND_VBK (self->dndg0), theta_k);
    
    ncm_vector_set (self->m2lnp_star, k, m2lnaps_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnaps_cur);
    /*printf ("%u % 22.15g | lnp_cur % 22.15g lnp_star % 22.15g\n", k, - 0.5 * (m2lnaps_cur - m2lnaps_star), - 0.5 * m2lnaps_cur, - 0.5 * m2lnaps_star);*/
  }

  if (k >= self->size_2)
  {
    const gdouble m2lnaps_star = ncm_stats_dist_nd_vbk_eval_m2lnp (NCM_STATS_DIST_ND_VBK (self->dndg1), thetastar);
    const gdouble m2lnaps_cur  = ncm_stats_dist_nd_vbk_eval_m2lnp (NCM_STATS_DIST_ND_VBK (self->dndg1), theta_k);
    
    ncm_vector_set (self->m2lnp_star, k, m2lnaps_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnaps_cur);
    /*printf ("%u % 22.15g | lnp_cur % 22.15g lnp_star % 22.15g\n", k, - 0.5 * (m2lnaps_cur - m2lnaps_star), - 0.5 * m2lnaps_cur, - 0.5 * m2lnaps_star);*/
  }
  
  /*ncm_vector_log_vals (theta_k,   "    THETA: ", "% 22.15g", TRUE);*/
  /*ncm_vector_log_vals (thetastar, "THETASTAR: ", "% 22.15g", TRUE);*/
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  const gdouble m2lnp_star = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur  = ncm_vector_get (self->m2lnp_cur, k);

/*
  printf ("AAA m2lnL_star % 22.15g m2lnp_star % 22.15g m2lnL_cur % 22.15g m2lnp_cur % 22.15g L cur->star: %12.5g p star->cur: %12.5g | T %12.5g\n",
          -0.5 * m2lnL_star, -0.5 * m2lnp_star, -0.5 * m2lnL_cur, -0.5 * m2lnp_cur,
          exp (- 0.5 * (m2lnL_star - m2lnL_cur)), exp (- 0.5 * (m2lnp_cur - m2lnp_star)),
          MIN (exp (- 0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur))), 1.0));
*/


  return exp (- 0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur)));
}

static gdouble 
_ncm_fit_esmcmc_walker_aps_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  const gdouble m2lnp_star = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur  = ncm_vector_get (self->m2lnp_cur, k);

  return -0.5 * (m2lnp_cur - m2lnp_star);
}

static void
_ncm_fit_esmcmc_walker_aps_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /*NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);*/
  /*NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;*/
  
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_aps_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPS *aps = NCM_FIT_ESMCMC_WALKER_APS (walker);  
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  return self->desc;
}

/**
 * ncm_fit_esmcmc_walker_aps_new:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 * 
 * Creates a new #NcmFitESMCMCWalkerAPS to be used
 * with @nwalkers.
 * 
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerAPS.
 */
NcmFitESMCMCWalkerAPS *
ncm_fit_esmcmc_walker_aps_new (guint nwalkers, guint nparams)
{
  NcmFitESMCMCWalkerAPS *aps = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_APS,
                                             "size", nwalkers,
                                             "nparams", nparams,
                                             NULL);
  
  return aps;
}

/**
 * ncm_fit_esmcmc_walker_aps_use_interp:
 * @aps: a #NcmFitESMCMCWalkerAPS
 * @use_interp: whether to use interpolation of the posterior
 * 
 * Sets whether to use interpolation of the posterior approximation (@use_interp == TRUE) 
 * or kernel density estimate (@use_interp == FALSE).
 * 
 */
void 
ncm_fit_esmcmc_walker_aps_use_interp (NcmFitESMCMCWalkerAPS *aps, gboolean use_interp)
{
  NcmFitESMCMCWalkerAPSPrivate * const self = aps->priv;
  
  self->use_interp = use_interp;
}

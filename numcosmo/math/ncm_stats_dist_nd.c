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

struct _NcmStatsDistNdPrivate
{
  guint dim;
};

enum
{
  PROP_0,
  PROP_DIM
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDistNd, ncm_stats_dist_nd, G_TYPE_OBJECT);

static void
ncm_stats_dist_nd_init (NcmStatsDistNd *dnd)
{
  NcmStatsDistNdPrivate * const self = dnd->priv = ncm_stats_dist_nd_get_instance_private (dnd);

  self->dim = 0;
}

static void
_ncm_stats_dist_nd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDistNd *dnd = NCM_STATS_DIST_ND (object);
  NcmStatsDistNdPrivate * const self = dnd->priv;
  
  g_return_if_fail (NCM_IS_STATS_DIST_ND (object));

  switch (prop_id)
  {
    case PROP_DIM:
      self->dim = g_value_get_uint (value);
      g_assert (NCM_STATS_DIST_ND_GET_CLASS (dnd)->set_dim != NULL);
      NCM_STATS_DIST_ND_GET_CLASS (dnd)->set_dim (dnd, g_value_get_uint (value));
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
      g_value_set_uint (value, self->dim);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_stats_dist_nd_dispose (GObject *object)
{
  /*NcmStatsDistNd *dnd = NCM_STATS_DIST_ND (object);*/
  /*NcmStatsDistNdPrivate * const self = dnd->priv;*/

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_nd_finalize (GObject *object)
{


  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_nd_parent_class)->finalize (object);
}

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
  klass->prepare           = NULL;
  klass->prepare_interp    = NULL;
  klass->set_dim           = NULL;
  klass->eval              = NULL;
  klass->eval_m2lnp        = NULL;
  klass->sample            = NULL;
  klass->kernel_sample     = NULL;
  klass->kernel_eval_m2lnp = NULL;
  klass->reset             = NULL;
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

  return self->dim;
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
 * ncm_stats_dist_nd_eval: (virtual eval)
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

  return dnd_class->eval (dnd, x);
}

/**
 * ncm_stats_dist_nd_eval_m2lnp: (virtual eval_m2lnp)
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

  return dnd_class->eval_m2lnp (dnd, x);
}

/**
 * ncm_stats_dist_nd_sample: (virtual sample)
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
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 
  dnd_class->sample (dnd, x, rng);
}

/**
 * ncm_stats_dist_nd_kernel_sample: (virtual kernel_sample)
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 * @mu: a #NcmVector
 * @scale: a double
 * @rng: a #NcmRNG
 * 
 * Using the pseudo-random number generator @rng generates a 
 * point from the distribution and copy it to @x.
 * 
 */
void
ncm_stats_dist_nd_kernel_sample (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *mu, const gdouble scale, NcmRNG *rng)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 
  dnd_class->kernel_sample (dnd, x, mu, scale, rng);
}

/**
 * ncm_stats_dist_nd_kernel_eval_m2lnp: (virtual kernel_eval_m2lnp)
 * @dnd: a #NcmStatsDistNd
 * @x: a #NcmVector
 * @y: a #NcmVector
 * @scale: covariance scale
 * 
 * Evaluates a single kernel at @x and @y and scale @s, i.e., $K_s(x,y)$.
 * 
 * Returns: $K_s(x,y)$.
 */
gdouble
ncm_stats_dist_nd_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *y, const gdouble scale)
{
  NcmStatsDistNdClass *dnd_class = NCM_STATS_DIST_ND_GET_CLASS (dnd); 
  return dnd_class->kernel_eval_m2lnp (dnd, x, y, scale);
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

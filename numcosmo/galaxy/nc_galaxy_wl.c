/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_likelihood.c
 *
 *  Mon May 08 16:12:03 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_likelihood.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION: nc_galaxy_wl_likelihood
 * @title: NcGalaxyWLLikelihood
 * @short_description: Class describing galaxy weak lensing distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy weak lensing distribution.
 * It is composed by three distributions: a shape distribution $P(s)$, a proxy redshift distribution $P(z_p)$, and a position distribution $P(z)P(r)$.
 * The shape distribution is defined by the abstract class #NcGalaxySDShape.
 * The proxy redshift distribution is defined by the abstract class #NcGalaxySDZProxy.
 * The position distribution is defined by the abstract class #NcGalaxySDPosition.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl_likelihood.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLLikelihoodPrivate
{
  NcGalaxySDShape *s_dist;
  NcGalaxySDZProxy *zp_dist;
  NcGalaxySDPosition *rz_dist;
};

enum
{
  PROP_0,
  PROP_S_DIST,
  PROP_ZP_DIST,
  PROP_RZ_DIST,
};

G_DEFINE_TYPE_WITH_PRIVATE(NcGalaxyWLLikelihood, nc_galaxy_wl_likelihood, G_TYPE_OBJECT);

static void
nc_galaxy_wl_likelihood_init (NcGalaxyWLLikelihood *gwl)
{
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv = nc_galaxy_wl_likelihood_get_instance_private (gwl);

  self->s_dist  = NULL;
  self->zp_dist = NULL;
  self->rz_dist = NULL;
}

static void
_nc_galaxy_wl_likelihood_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLLikelihood *gwl = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_S_DIST:
      self->s_dist = g_value_dup_object (value);
      break;
    case PROP_ZP_DIST:
      self->zp_dist = g_value_dup_object (value);
      break;
    case PROP_RZ_DIST:
      self->rz_dist = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}


static void
_nc_galaxy_wl_likelihood_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLLikelihood *gwl = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_S_DIST:
      g_value_set_object (value, self->s_dist);
      break;
    case PROP_ZP_DIST:
      g_value_set_object (value, self->zp_dist);
      break;
    case PROP_RZ_DIST:
      g_value_set_object (value, self->rz_dist);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_likelihood_dispose (GObject *object)
{
  NcGalaxyWLLikelihood *gwl = NC_GALAXY_WL_LIKELIHOOD (object);
  NcGalaxyWLLikelihoodPrivate * const self = gwl->priv;

  /* FIX ME */
  /* nc_galaxy_sd_shape_clear (&self->s_dist); */
  /* nc_galaxy_sd_z_proxy_clear (&self->zp_dist); */
  /* nc_galaxy_sd_position_clear (&self->rz_dist); */

  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_likelihood_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_wl_likelihood_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_likelihood_class_init (NcGalaxyWLLikelihoodClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_likelihood_set_property;
  object_class->get_property = &_nc_galaxy_wl_likelihood_get_property;
  object_class->dispose = &_nc_galaxy_wl_likelihood_dispose;
  object_class->finalize = &_nc_galaxy_wl_likelihood_finalize;

  /**
   * NcGalaxyWLLikelihood:s-dist:
   *
   * A #NcGalaxySDShape object.
   *
   */

  /* FIX ME */
  /* g_object_class_install_property (object_class,
                                   PROP_S_DIST.
                                   g_param_spec_object ("s-dist",
                                                        NULL,
                                                        "Galaxy sample shape distribution",
                                                        NC_TYPE_GALAXY_SD_SHAPE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); */

  /**
   * NcGalaxyWLLikelihood:zp-dist:
   *
   * A #NcGalaxySDZProxy object.
   *
   */

  /* FIX ME */
  /* g_object_class_install_property (object_class,
                                   PROP_ZP_DIST,
                                   g_param_spec_object ("zp-dist",
                                                        NULL,
                                                        "Galaxy sample proxy redshift distribution",
                                                        NC_TYPE_GALAXY_SD_Z_PROXY,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); */

  /**
   * NcGalaxyWLLikelihood:rz-dist:
   *
   * A #NcGalaxySDZPosition object.
   *
   */

  /* FIX ME */
  /* g_object_class_install_property (object_class,
                                   PROP_RZ_DIST,
                                   g_param_spec_object ("rz-dist",
                                                        NULL,
                                                        "Galaxy sample position distribution",
                                                        NC_TYPE_GALAXY_SD_POSITION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); */
}

/**
 * nc_galaxy_wl_likelihood_new:
 * @s_dist: a #NcGalaxySDShape
 * @zp_dist: a #NcGalaxySDZProxy
 * @rz_dist: a #NcGalaxySDPosition
 *
 * Creates a new galaxy weak lensing object.
 * Requires an instance of #NcGalaxySDShape, #NcGalaxySDZProxy, and #NcGalaxySDPosition.
 *
 * Returns: (transfer full): a new NcGalaxyWLLikelihood.
 */
NcGalaxyWLLikelihood *
nc_galaxy_wl_likelihood_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *rz_dist)
{
  NcGalaxyWLLikelihood *gwl = g_object_new (NC_TYPE_GALAXY_WL_LIKELIHOOD,
                                  "s-dist", s_dist,
                                  "zp-dist", zp_dist,
                                  "zp-dist", rz_dist,
                                  NULL);

  return gwl;
}

/**
 * nc_galaxy_wl_likelihood_ref:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Increase the reference of @gwl by one.
 *
 * Returns: (transfer full): @gwl.
 */
NcGalaxyWLLikelihood *
nc_galaxy_wl_likelihood_ref (NcGalaxyWLLikelihood *gwl)
{
  return g_object_ref (gwl);
}

/**
 * nc_galaxy_wl_likelihood_free:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Decrease the reference count of @gwl by one.
 *
 */
void
nc_galaxy_wl_likelihood_free (NcGalaxyWLLikelihood *gwl)
{
  g_object_unref (gwl);
}

/**
 * nc_galaxy_wl_likelihood_clear:
 * @gwl: a #NcGalaxyWLLikelihood
 *
 * Decrease the reference count of @gwl by one, and sets the pointer *@gwl to
 * NULL.
 *
 */
void
nc_galaxy_wl_likelihood_clear (NcGalaxyWLLikelihood **gwl)
{
  g_clear_object (gwl);
}
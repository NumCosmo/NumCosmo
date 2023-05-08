/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl.c
 *
 *  Mon May 08 16:12:03 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl.c
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
 * SECTION:nc_galaxy_wl
 * @title: NcGalaxyWL
 * @short_description: Class describing galaxy weak lensing distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy weak lensing distribution.
 * It is composed by three distributions: a shape distribution $P(s)$, a proxy redshift distribution $P(z_p)$, and a position distribution $P(z)P(r)$.
 * The shape distribution is defined by the abstract class #NcGSDShape.
 * The proxy redshift distribution is defined by the abstract class #NcGSDZProxy.
 * The position distribution is defined by the abstract class #NcGSDPosition.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl.h"
#include "galaxy/nc_gsd_shape.h"
#include "galaxy/nc_gsd_z_proxy.h"
#include "galaxy/nc_gsd_position.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLPrivate
{
  NcGSDShape *s_dist;
  NcGSDZProxy *zp_dist;
  NcGSDPosition *rz_dist;
};

enum
{
  PROP_0,
  PROP_S_DIST,
  PROP_ZP_DIST,
  PROP_RZ_DIST,
};

G_DEFINE_TYPE_WITH_PRIVATE(NcGalaxyWL, nc_galaxy_wl, G_TYPE_OBJECT);

static void
nc_galaxy_wl_init (NcGalaxyWL *gwl)
{
  NcGalaxyWLprivate * const self = gwl->priv = nc_galaxy_wl_get_instance_private (gwl);

  self->s_dist  = NULL;
  self->zp_dist = NULL;
  self->rz_dist = NULL;
}

static void
_nc_galaxy_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWL *gwl = NC_GALAXY_WL (object);
  NcGalaxyWLPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL (object));

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
_nc_galaxy_wl_get_property (GObject *object, guint prop_id, GValue, *value, GParamSpec *pspec)
{
  NcGalaxyWL *gwl = NC_GALAXY_WL (object);
  NcGalaxyWLPrivate * const self = gwl->priv;

  g_return_if_fail (NC_IS_GALAXY_WL (object));

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
    default;
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_dispose (GObject *object)
{
  NcGalaxyWL *gwl = NC_GALAXY_WL (object);
  NcGalaxyWLPrivate * const self = gwl->priv;

  /* FIX ME */
  /* nc_gsd_shape_clear (&self->s_dist); */
  /* nc_gsd_z_proxy_clear (&self->zp_dist); */
  /* nc_gsd_position_clear (&self->rz_dist); */

  G_OBJECT_CLASS (nc_galaxy_wl_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_wl_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_class_init (NcGalaxyWLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_set_property;
  object_class->get_property = &_nc_galaxy_wl_get_property;
  object_class->dispose = &_nc_galaxy_wl_dispose;
  object_class->finalize = &_nc_galaxy_wl_finalize;

  /**
   * NcGalaxyWL:s-dist:
   *
   * A #NcGSDShape object.
   *
   */

  /* FIX ME */
  /* g_object_class_install_property (object_class,
                                   PROP_S_DIST.
                                   g_param_spec_object ("s-dist",
                                                        NULL,
                                                        "Galaxy sample shape distribution",
                                                        NC_TYPE_GSD_SHAPE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); */

  /**
   * NcGalaxyWL:zp-dist:
   *
   * A #NcGSDZProxy object.
   *
   */

  /* FIX ME */
  /* g_object_class_install_property (object_class,
                                   PROP_ZP_DIST,
                                   g_param_spec_object ("zp-dist",
                                                        NULL,
                                                        "Galaxy sample proxy redshift distribution",
                                                        NC_TYPE_GSD_Z_PROXY,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); */

  /**
   * NcGalaxyWL:rz-dist:
   *
   * A #NcGSDZPosition object.
   *
   */

  /* FIX ME */
  /* g_object_class_install_property (object_class,
                                   PROP_RZ_DIST,
                                   g_param_spec_object ("rz-dist",
                                                        NULL,
                                                        "Galaxy sample position distribution",
                                                        NC_TYPE_GSD_POSITION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); */
}

/**
 * nc_galaxy_wl_new:
 * @s_dist: a #NcGSDShape
 * @zp_dist: a #NcGSDZProxy
 * @rz_dist: a #NcGSDPosition
 *
 * Creates a new galaxy weak lensing object.
 * Requires an instance of #NcGSDShape, #NcGSDZProxy, and #NcGSDPosition.
 *
 * Returns: (transfer full): a new NcGalaxyWL.
 */
NcGalaxyWL *
nc_galaxy_wl_new (NcGSDShape *s_dist, NcGSDZProxy *zp_dist, NcGSDPosition *rz_dist)
{
  NcGalaxyWL *gwl = g_object_new (NC_TYPE_GALAXY_WL,
                                  "s-dist", s_dist,
                                  "zp-dist", zp_dist,
                                  "zp-dist", rz_dist,
                                  NULL);

  return gwl;
}

/**
 * nc_galaxy_wl_ref:
 * @gwl: a #NcGalaxyWL
 *
 * Increase the reference of @gwl by one.
 *
 * Returns: (transfer full): @gwl.
 */
NcGalaxyWL *
nc_galaxy_wl_ref (NcGalaxyWL *gwl)
{
  return g_object_ref (gwl);
}

/**
 * nc_galaxy_wl_free:
 * @gwl: a #NcGalaxyWL
 *
 * Decrease the reference count of @gwl by one.
 *
 */
void
nc_galaxy_wl_free (NcGalaxyWL *gwl)
{
  g_object_unref (gwl);
}

/**
 * nc_galaxy_wl_clear:
 * @gwl: a #NcGalaxyWL
 *
 * Decrease the reference count of @gwl by one, and sets the pointer *@gwl to
 * NULL.
 *
 */
void
nc_galaxy_wl_clear (NcGalaxyWL **gwl)
{
  g_clear_object (gwl);
}
/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_redshift.c
 *
 *  Thu Aug 1 00:45:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION: nc_galaxy_redshift
 * @title: NcGalaxyRedshift
 * @short_description: Class describing galaxy observed redshift.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy observed redshift. It is used to store the
 * redshift of a galaxy and its associated error.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxyRedshiftPrivate
{
  gint placeholder;
} NcGalaxyRedshiftPrivate;

enum
{
  PROP_0,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshift, nc_galaxy_redshift, G_TYPE_OBJECT);

static void
nc_galaxy_redshift_init (NcGalaxyRedshift *gz)
{
}

static void
_nc_galaxy_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshift *gz                 = NC_GALAXY_REDSHIFT (object);
  NcGalaxyRedshiftPrivate * const self = nc_galaxy_redshift_get_instance_private (gz);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT (gz));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshift *gz = NC_GALAXY_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT (gz));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_parent_class)->finalize (object);
}

/*  LCOV_EXCL_START */
static gdouble
_nc_galaxy_redshift_gen (NcGalaxyRedshift *gz, NcmRNG *rng, NcmVector *data)
{
  g_error ("_nc_galaxy_redshift_gen: method not implemented");

  return 0.0;
}

static gdouble
_nc_galaxy_redshift_integ (NcGalaxyRedshift *gz, NcmVector *data)
{
  g_error ("_nc_galaxy_redshift_integ: method not implemented");

  return 0.0;
}

static gboolean
_nc_galaxy_redshift_get_header (NcGalaxyRedshift *gz)
{
  g_error ("_nc_galaxy_redshift_get_header: method not implemented");

  return FALSE;
}

/*  LCOV_EXCL_STOP */

static void
nc_galaxy_redshift_class_init (NcGalaxyRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_finalize;

  klass->gen        = &_nc_galaxy_redshift_gen;
  klass->integ      = &_nc_galaxy_redshift_integ;
  klass->get_header = &_nc_galaxy_redshift_get_header;
}

/**
 * nc_galaxy_redshift_ref:
 * @gz: a #NcGalaxyRedshift
 *
 * Increases the reference count of @gz by one.
 *
 * Returns: (transfer full): @gz
 */
NcGalaxyRedshift *
nc_galaxy_redshift_ref (NcGalaxyRedshift *gz)
{
  return g_object_ref (gz);
}

/**
 * nc_galaxy_redshift_free:
 * @gz: a #NcGalaxyRedshift
 *
 * Decreases the reference count of @gz by one.
 *
 */
void
nc_galaxy_redshift_free (NcGalaxyRedshift *gz)
{
  g_object_unref (gz);
}

/**
 * nc_galaxy_redshift_clear:
 * @gz: a #NcGalaxyRedshift
 *
 * Decreases the reference count of @gz by one, and sets the pointer *@gz to
 * NULL.
 *
 */
void
nc_galaxy_redshift_clear (NcGalaxyRedshift **gz)
{
  g_clear_object (gz);
}

/**
 * nc_galaxy_redshift_gen:
 * @gz: a #NcGalaxyRedshift
 * @rng: a #NcmRNG
 * @data: a #NcmVector
 *
 * Generates an observed redshift value from the distribution using @rng.
 *
 * Returns: the generated redshift value.
 */
gdouble
nc_galaxy_redshift_gen (NcGalaxyRedshift *gz, NcmRNG *rng, NcmVector *data)
{
  return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->gen (gz, rng, data);
}

/**
 * nc_galaxy_redshift_integ:
 * @gz: a #NcGalaxyRedshift
 * @data: a #NcmVector
 *
 * Computes the probability density of the observable $z_p$.
 *
 * Returns: the probability density at $z_p$, $P(z_p)$.
 */
gdouble
nc_galaxy_redshift_integ (NcGalaxyRedshift *gz, NcmVector *data)
{
  return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->integ (gz, data);
}

/**
 * nc_galaxy_redshift_get_header:
 * @gz: a #NcGalaxyRedshift
 *
 * Gets the header of the galaxy redshift data.
 *
 * Returns: the header of the galaxy redshift data.
 */
gboolean
nc_galaxy_redshift_get_header (NcGalaxyRedshift *gz)
{
  return NC_GALAXY_REDSHIFT_GET_CLASS (gz)->get_header (gz);
}


/***************************************************************************
 *            nc_galaxy_redshift_spec.c
 *
 *  Thu Aug 1 15:10:22 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_spec.c
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
 * SECTION:nc_galaxy_redshift_spec
 * @title: NcGalaxyRedshiftSpec
 * @short_description: Class describing spectroscopic redshift observations
 * @stability: Unstable
 *
 *
 * Class describing spectroscopic redshift observations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_redshift.h"
#include "galaxy/nc_galaxy_redshift_spec.h"
#include "galaxy/nc_galaxy_sd_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxyRedshiftSpecPrivate
{
  NcGalaxySDRedshift *sdz;
} NcGalaxyRedshiftSpecPrivate;

struct _NcGalaxyRedshiftSpec
{
  NcGalaxyRedshift parent_instance;
};

enum
{
  PROP_0,
  PROP_SDZ,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftSpec, nc_galaxy_redshift_spec, NC_TYPE_GALAXY_REDSHIFT);

static void
nc_galaxy_redshift_spec_init (NcGalaxyRedshiftSpec *gzspec)
{
  NcGalaxyRedshiftSpecPrivate * const self = nc_galaxy_redshift_spec_get_instance_private (gzspec);

  self->sdz = NULL;
}

static void
_nc_galaxy_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftSpec *gzspec             = NC_GALAXY_REDSHIFT_SPEC (object);
  NcGalaxyRedshiftSpecPrivate * const self = nc_galaxy_redshift_spec_get_instance_private (gzspec);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_SPEC (gzspec));

  switch (prop_id)
  {
    case PROP_SDZ:
      self->sdz = g_value_get_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftSpec *gzspec             = NC_GALAXY_REDSHIFT_SPEC (object);
  NcGalaxyRedshiftSpecPrivate * const self = nc_galaxy_redshift_spec_get_instance_private (gzspec);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_SPEC (gzspec));

  switch (prop_id)
  {
    case PROP_SDZ:
      g_value_set_object (value, self->sdz);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_redshift_spec_dispose (GObject *object)
{
  NcGalaxyRedshiftSpec *gzspec             = NC_GALAXY_REDSHIFT_SPEC (object);
  NcGalaxyRedshiftSpecPrivate * const self = nc_galaxy_redshift_spec_get_instance_private (gzspec);

  nc_galaxy_sd_redshift_clear (&self->sdz);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_redshift_spec_parent_class)->dispose (object);
}

static void
nc_galaxy_redshift_spec_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_redshift_spec_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_redshift_spec_gen (NcGalaxyRedshift *gz, NcmRNG *rng, NcmVector *data);
static gdouble _nc_galaxy_redshift_spec_integ (NcGalaxyRedshift *gz, gdouble z, NcmVector *data);
static GStrv _nc_galaxy_redshift_spec_get_header (NcGalaxyRedshift *gz);

static void
nc_galaxy_sd_redshift_spec_class_init (NcGalaxyRedshiftSpecClass *klass)
{
  NcGalaxyRedshiftClass *gz_class = NC_GALAXY_REDSHIFT_CLASS (klass);
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_redshift_spec_dispose;
  object_class->finalize = &nc_galaxy_redshift_spec_finalize;

  /**
   * NcGalaxyRedshiftSpec:sdz:
   *
   * The galaxy redshift sample distribution.
   *
   */
  g_object_class_install_property (object_class, PROP_SDZ,
                                   g_param_spec_object ("sdz",
                                                        "Galaxy redshift sample distribution",
                                                        "The galaxy redshift sample distribution",
                                                        NC_TYPE_GALAXY_SD_REDSHIFT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gz_class->gen        = &_nc_galaxy_redshift_spec_gen;
  gz_class->integ      = &_nc_galaxy_redshift_spec_integ;
  gz_class->get_header = &_nc_galaxy_redshift_spec_get_header;
}

static gdouble
_nc_galaxy_redshift_spec_gen (NcGalaxyRedshift *gz, NcmRNG *rng, NcmVector *data)
{
  NcGalaxyRedshiftSpec *gzspec             = NC_GALAXY_REDSHIFT_SPEC (gz);
  NcGalaxyRedshiftSpecPrivate * const self = nc_galaxy_redshift_spec_get_instance_private (gzspec);

  return nc_galaxy_sd_redshift_gen (self->sdz, rng);
}

static gdouble
_nc_galaxy_redshift_spec_integ (NcGalaxyRedshift *gz, gdouble z, NcmVector *data)
{
  NcGalaxyRedshiftSpec *gzspec             = NC_GALAXY_REDSHIFT_SPEC (gz);
  NcGalaxyRedshiftSpecPrivate * const self = nc_galaxy_redshift_spec_get_instance_private (gzspec);

  return nc_galaxy_sd_redshift_integ (self->sdz, z);
}

static GStrv
_nc_galaxy_redshift_spec_get_header (NcGalaxyRedshift *gz)
{
  GStrv header = g_strsplit ("z", " ", -1);

  return header;
}

/**
 * nc_galaxy_redshift_spec_new:
 * @sdz: The galaxy redshift sample distribution.
 *
 * Creates a new #NcGalaxyRedshiftSpec object.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftSpec object.
 */
NcGalaxyRedshiftSpec *
nc_galaxy_redshift_spec_new (NcGalaxySDRedshift *sdz)
{
  return g_object_new (NC_TYPE_GALAXY_REDSHIFT_SPEC, "sdz", sdz, NULL);
}

/**
 * nc_galaxy_redshift_spec_ref:
 * @gzspec: a #NcGalaxyRedshiftSpec
 *
 * Increases the reference count of @gzspec by one.
 *
 * Returns: (transfer full): @gzspec.
 */
NcGalaxyRedshiftSpec *
nc_galaxy_redshift_spec_ref (NcGalaxyRedshiftSpec *gzspec)
{
  return g_object_ref (gzspec);
}

/**
 * nc_galaxy_redshift_spec_free:
 * @gzspec: a #NcGalaxyRedshiftSpec
 *
 * Decreases the reference count of @gzspec by one.
 *
 */
void
nc_galaxy_redshift_spec_free (NcGalaxyRedshiftSpec *gzspec)
{
  g_object_unref (gzspec);
}

/**
 * nc_galaxy_redshift_spec_clear:
 * @gzspec: a #NcGalaxyRedshiftSpec
 *
 * Decreases the reference count of @gzspec by one, and sets the pointer *@gzspec to
 * NULL.
 *
 */
void
nc_galaxy_redshift_spec_clear (NcGalaxyRedshiftSpec **gzspec)
{
  g_clear_object (gzspec);
}


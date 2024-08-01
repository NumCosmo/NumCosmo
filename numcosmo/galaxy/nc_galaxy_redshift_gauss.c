/***************************************************************************
 *            nc_galaxy_redshift_gauss.c
 *
 *  Thu Aug 1 20:03:55 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_gauss.c
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
 * SECTION:nc_galaxy_redshift_gauss
 * @title: NcGalaxyRedshiftGauss
 * @short_description: Class describing photometric redshift observations with gaussian errors.
 * @stability: Unstable
 *
 *
 * Class describing photometric redshift observations with gaussian errors.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "gsl/gsl_randist.h"
#include "galaxy/nc_galaxy_redshift.h"
#include "galaxy/nc_galaxy_redshift_gauss.h"
#include "galaxy/nc_galaxy_sd_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxyRedshiftGaussPrivate
{
  NcGalaxySDRedshift *sdz;
} NcGalaxyRedshiftGaussPrivate;

struct _NcGalaxyRedshiftGauss
{
  NcGalaxyRedshift parent_instance;
};

enum
{
  PROP_0,
  PROP_SDZ,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftGauss, nc_galaxy_redshift_gauss, NC_TYPE_GALAXY_REDSHIFT);

static void
nc_galaxy_redshift_gauss_init (NcGalaxyRedshiftGauss *gzgauss)
{
  NcGalaxyRedshiftGaussPrivate * const self = nc_galaxy_redshift_gauss_get_instance_private (gzgauss);

  self->sdz = NULL;
}

static void
_nc_galaxy_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftGauss *gzgauss            = NC_GALAXY_REDSHIFT_GAUSS (object);
  NcGalaxyRedshiftGaussPrivate * const self = nc_galaxy_redshift_gauss_get_instance_private (gzgauss);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_GAUSS (gzgauss));

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
  NcGalaxyRedshiftGauss *gzgauss            = NC_GALAXY_REDSHIFT_GAUSS (object);
  NcGalaxyRedshiftGaussPrivate * const self = nc_galaxy_redshift_gauss_get_instance_private (gzgauss);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_GAUSS (gzgauss));

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
_nc_galaxy_redshift_gauss_dispose (GObject *object)
{
  NcGalaxyRedshiftGauss *gzgauss            = NC_GALAXY_REDSHIFT_GAUSS (object);
  NcGalaxyRedshiftGaussPrivate * const self = nc_galaxy_redshift_gauss_get_instance_private (gzgauss);

  nc_galaxy_sd_redshift_clear (&self->sdz);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_redshift_gauss_parent_class)->dispose (object);
}

static void
nc_galaxy_redshift_gauss_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_redshift_gauss_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_redshift_gauss_gen (NcGalaxyRedshift *gz, NcmRNG *rng, NcmVector *data);
static gdouble _nc_galaxy_redshift_gauss_integ (NcGalaxyRedshift *gz, gdouble z, NcmVector *data);
static GStrv _nc_galaxy_redshift_gauss_get_header (NcGalaxyRedshift *gz);

static void
nc_galaxy_sd_redshift_gauss_class_init (NcGalaxyRedshiftGaussClass *klass)
{
  NcGalaxyRedshiftClass *gz_class = NC_GALAXY_REDSHIFT_CLASS (klass);
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_redshift_gauss_dispose;
  object_class->finalize = &nc_galaxy_redshift_gauss_finalize;

  /**
   * NcGalaxyRedshiftGauss:sdz:
   *
   * The galaxy redshift sample distribution.
   *
   */
  g_object_class_install_property (object_class, PROP_SDZ,
                                   g_param_gauss_object ("sdz",
                                                         "Galaxy redshift sample distribution",
                                                         "The galaxy redshift sample distribution",
                                                         NC_TYPE_GALAXY_SD_REDSHIFT,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gz_class->gen        = &_nc_galaxy_redshift_gauss_gen;
  gz_class->integ      = &_nc_galaxy_redshift_gauss_integ;
  gz_class->get_header = &_nc_galaxy_redshift_gauss_get_header;
}

static gdouble
_nc_galaxy_redshift_gauss_gen (NcGalaxyRedshift *gz, NcmRNG *rng, NcmVector *data)
{
  NcGalaxyRedshiftGauss *gzgauss            = NC_GALAXY_REDSHIFT_GAUSS (gz);
  NcGalaxyRedshiftGaussPrivate * const self = nc_galaxy_redshift_gauss_get_instance_private (gzgauss);
  gdouble sigma                             = ncm_vector_get (data, 1);
  gdouble z_min                             = 0.0;
  gdouble z_max                             = 0.0;
  gdouble zp;
  gdouble z;

  if (!nc_galaxy_sd_redshift_get_lim (self->sdz, &z_min, &z_max))
    g_error ("Failed to get redshift limits.");

  do {
    z  = nc_galaxy_sd_redshift_gen (self->sdz, rng);
    zp = ncm_rng_gaussian_gen (rng, z, sigma * (1.0 + z));
  } while ((zp > z_max) || (zp < z_min));

  return zp;
}

static gdouble
_nc_galaxy_redshift_gauss_integ (NcGalaxyRedshift *gz, gdouble z, NcmVector *data)
{
  NcGalaxyRedshiftGauss *gzgauss            = NC_GALAXY_REDSHIFT_GAUSS (gz);
  NcGalaxyRedshiftGaussPrivate * const self = nc_galaxy_redshift_gauss_get_instance_private (gzgauss);
  gdouble zp                                = ncm_vector_get (data, 0);
  gdouble sigma                             = ncm_vector_get (data, 1);

  return nc_galaxy_sd_redshift_integ (self->sdz, z) * gsl_ran_gaussian_pdf (zp - z, sigma * (1.0 + z));
}

static GStrv
_nc_galaxy_redshift_gauss_get_header (NcGalaxyRedshift *gz)
{
  GStrv header = g_strsplit ("zp zp_sigma", " ", -1);

  return header;
}

/**
 * nc_galaxy_redshift_gauss_new:
 * @sdz: The galaxy redshift sample distribution.
 *
 * Creates a new #NcGalaxyRedshiftGauss object.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftGauss object.
 */
NcGalaxyRedshiftGauss *
nc_galaxy_redshift_gauss_new (NcGalaxySDRedshift *sdz)
{
  return g_object_new (NC_TYPE_GALAXY_REDSHIFT_GAUSS, "sdz", sdz, NULL);
}

/**
 * nc_galaxy_redshift_gauss_ref:
 * @gzgauss: a #NcGalaxyRedshiftGauss
 *
 * Increases the reference count of @gzgauss by one.
 *
 * Returns: (transfer full): @gzgauss.
 */
NcGalaxyRedshiftGauss *
nc_galaxy_redshift_gauss_ref (NcGalaxyRedshiftGauss *gzgauss)
{
  return g_object_ref (gzgauss);
}

/**
 * nc_galaxy_redshift_gauss_free:
 * @gzgauss: a #NcGalaxyRedshiftGauss
 *
 * Decreases the reference count of @gzgauss by one.
 *
 */
void
nc_galaxy_redshift_gauss_free (NcGalaxyRedshiftGauss *gzgauss)
{
  g_object_unref (gzgauss);
}

/**
 * nc_galaxy_redshift_gauss_clear:
 * @gzgauss: a #NcGalaxyRedshiftGauss
 *
 * Decreases the reference count of @gzgauss by one, and sets the pointer *@gzgauss to
 * NULL.
 *
 */
void
nc_galaxy_redshift_gauss_clear (NcGalaxyRedshiftGauss **gzgauss)
{
  g_clear_object (gzgauss);
}


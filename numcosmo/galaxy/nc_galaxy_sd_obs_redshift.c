/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_obs_redshift.c
 *
 *  Thu Aug 1 00:45:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift.c
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
 * SECTION: nc_galaxy_sd_obs_redshift
 * @title: NcGalaxySDObsRedshift
 * @short_description: Class describing galaxy sample observed redshift
 * distribution.
 * @stability: Unstable
 *
 *
 * This class describes galaxy sample observed redshift distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxySDObsRedshiftPrivate
{
  gint placeholder;
} NcGalaxySDObsRedshiftPrivate;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshift, nc_galaxy_sd_obs_redshift, NCM_TYPE_MODEL);

static void
nc_galaxy_sd_obs_redshift_init (NcGalaxySDObsRedshift *gsdor)
{
}

static void
_nc_galaxy_sd_obs_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshift *gsdor              = NC_GALAXY_SD_OBS_REDSHIFT (object);
  NcGalaxySDObsRedshiftPrivate * const self = nc_galaxy_sd_obs_redshift_get_instance_private (gsdor);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshift *gsdor = NC_GALAXY_SD_OBS_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_obs_redshift, NC_TYPE_GALAXY_SD_OBS_REDSHIFT)

/*  LCOV_EXCL_START */
static gdouble
_nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_obs_redshift_gen: method not implemented");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_obs_redshift_integ: method not implemented");

  return 0.0;
}

static GStrv
_nc_galaxy_sd_obs_redshift_get_header (NcGalaxySDObsRedshift *gsdor)
{
  g_error ("_nc_galaxy_sd_obs_redshift_get_header: method not implemented");

  return NULL;
}

/*  LCOV_EXCL_STOP */

static void
nc_galaxy_sd_obs_redshift_class_init (NcGalaxySDObsRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_get_property;
  object_class->finalize    = &_nc_galaxy_sd_obs_redshift_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample observed redshift distribution", "GalaxySDObsRedshift");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  ncm_model_class_check_params_info (model_class);

  klass->gen        = &_nc_galaxy_sd_obs_redshift_gen;
  klass->integ      = &_nc_galaxy_sd_obs_redshift_integ;
  klass->get_header = &_nc_galaxy_sd_obs_redshift_get_header;
}

/**
 * nc_galaxy_sd_obs_redshift_ref:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Increases the reference count of @gsdor by one.
 *
 * Returns: (transfer full): @gsdor
 */
NcGalaxySDObsRedshift *
nc_galaxy_sd_obs_redshift_ref (NcGalaxySDObsRedshift *gsdor)
{
  return g_object_ref (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_free:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Decreases the reference count of @gsdor by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_free (NcGalaxySDObsRedshift *gsdor)
{
  g_object_unref (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_clear:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Decreases the reference count of @gsdor by one, and sets the pointer *@gsdor to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_clear (NcGalaxySDObsRedshift **gsdor)
{
  g_clear_object (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_gen:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @rng: a #NcmRNG
 * @data: a #NcmVector
 *
 * Generates an observed redshift value from the distribution using @rng.
 *
 * Returns: the generated redshift value.
 */
gdouble
nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->gen (gsdor, rng, data);
}

/**
 * nc_galaxy_sd_obs_redshift_integ:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @data: a #NcmVector
 *
 * Computes the probability density of the observable $z_p$.
 *
 * Returns: the probability density at $z_p$, $P(z_p)$.
 */
gdouble
nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->integ (gsdor, z, data);
}

/**
 * nc_galaxy_sd_obs_redshift_get_header:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Gets the header of the galaxy redshift data.
 *
 * Returns: (transfer full): the header of the galaxy redshift data.
 */
GStrv
nc_galaxy_sd_obs_redshift_get_header (NcGalaxySDObsRedshift *gsdor)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->get_header (gsdor);
}


/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_spec.c
 *
 *  Thu Aug 1 15:10:22 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_spec.c
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
 * SECTION:nc_galaxy_sd_obs_redshift_spec
 * @title: NcGalaxySDObsRedshiftSpec
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

#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_spec.h"
#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxySDObsRedshiftSpecPrivate
{
  NcGalaxySDTrueRedshift *sdz;
} NcGalaxySDObsRedshiftSpecPrivate;

struct _NcGalaxySDObsRedshiftSpec
{
  NcGalaxySDObsRedshift parent_instance;
};

enum
{
  PROP_0,
  PROP_SDZ,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshiftSpec, nc_galaxy_sd_obs_redshift_spec, NC_TYPE_GALAXY_SD_OBS_REDSHIFT);

static void
nc_galaxy_sd_obs_redshift_spec_init (NcGalaxySDObsRedshiftSpec *gsdorspec)
{
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  self->sdz = NULL;
}

static void
_nc_galaxy_sd_obs_redshift_spec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (object);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdorspec));

  switch (prop_id)
  {
    case PROP_SDZ:
      self->sdz = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_spec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (object);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdorspec));

  switch (prop_id)
  {
    case PROP_SDZ:
      g_value_set_object (value, self->sdz);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_spec_dispose (GObject *object)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (object);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  nc_galaxy_sd_true_redshift_clear (&self->sdz);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_spec_parent_class)->dispose (object);
}

static void
nc_galaxy_sd_obs_redshift_spec_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_spec_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_obs_redshift_spec_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data);
static gdouble _nc_galaxy_sd_obs_redshift_spec_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data);
static GStrv _nc_galaxy_sd_obs_redshift_spec_get_header (NcGalaxySDObsRedshift *gsdor);

static void
nc_galaxy_sd_obs_redshift_spec_class_init (NcGalaxySDObsRedshiftSpecClass *klass)
{
  NcGalaxySDObsRedshiftClass *gsdor_class = NC_GALAXY_SD_OBS_REDSHIFT_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_spec_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_spec_get_property;
  object_class->dispose     = &_nc_galaxy_sd_obs_redshift_spec_dispose;
  object_class->finalize    = &nc_galaxy_sd_obs_redshift_spec_finalize;

  ncm_model_class_set_name_nick (model_class, "Spectroscopic Observed Redshift", "GalaxySDObsRedshiftSpec");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDObsRedshiftSpec:sdz:
   *
   * The galaxy redshift sample distribution.
   *
   */
  g_object_class_install_property (object_class, PROP_SDZ,
                                   g_param_spec_object ("sdz",
                                                        "Galaxy redshift sample distribution",
                                                        "The galaxy redshift sample distribution",
                                                        NC_TYPE_GALAXY_SD_TRUE_REDSHIFT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen        = &_nc_galaxy_sd_obs_redshift_spec_gen;
  gsdor_class->integ      = &_nc_galaxy_sd_obs_redshift_spec_integ;
  gsdor_class->get_header = &_nc_galaxy_sd_obs_redshift_spec_get_header;
}

static void
_nc_galaxy_sd_obs_redshift_spec_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdor);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  ncm_vector_set (data, 0, nc_galaxy_sd_true_redshift_gen (self->sdz, rng));
}

static gdouble
_nc_galaxy_sd_obs_redshift_spec_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdor);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  return nc_galaxy_sd_true_redshift_integ (self->sdz, z);
}

static GStrv
_nc_galaxy_sd_obs_redshift_spec_get_header (NcGalaxySDObsRedshift *gsdor)
{
  GStrv header = g_strsplit ("z", " ", -1);

  return header;
}

/**
 * nc_galaxy_sd_obs_redshift_spec_new:
 * @sdz: a #NcGalaxySDTrueRedshift.
 *
 * Creates a new #NcGalaxySDObsRedshiftSpec object.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftSpec object.
 */
NcGalaxySDObsRedshiftSpec *
nc_galaxy_sd_obs_redshift_spec_new (NcGalaxySDTrueRedshift *sdz)
{
  return g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_SPEC, "sdz", sdz, NULL);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_ref:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 *
 * Increases the reference count of @gsdorspec by one.
 *
 * Returns: (transfer full): @gsdorspec.
 */
NcGalaxySDObsRedshiftSpec *
nc_galaxy_sd_obs_redshift_spec_ref (NcGalaxySDObsRedshiftSpec *gsdorspec)
{
  return g_object_ref (gsdorspec);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_free:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 *
 * Decreases the reference count of @gsdorspec by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_spec_free (NcGalaxySDObsRedshiftSpec *gsdorspec)
{
  g_object_unref (gsdorspec);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_clear:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 *
 * Decreases the reference count of @gsdorspec by one, and sets the pointer *@gsdorspec to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_spec_clear (NcGalaxySDObsRedshiftSpec **gsdorspec)
{
  g_clear_object (gsdorspec);
}


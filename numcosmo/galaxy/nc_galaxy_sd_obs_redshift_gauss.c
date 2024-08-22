/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_gauss.c
 *
 *  Thu Aug 1 20:03:55 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_gauss.c
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
 * SECTION:nc_galaxy_sd_obs_redshift_gauss
 * @title: NcGalaxySDObsRedshiftGauss
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
#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_gauss.h"
#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxySDObsRedshiftGaussPrivate
{
  NcGalaxySDTrueRedshift *sdz;
} NcGalaxySDObsRedshiftGaussPrivate;

struct _NcGalaxySDObsRedshiftGauss
{
  NcGalaxySDObsRedshift parent_instance;
};

enum
{
  PROP_0,
  PROP_SDZ,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshiftGauss, nc_galaxy_sd_obs_redshift_gauss, NC_TYPE_GALAXY_SD_OBS_REDSHIFT);

static void
nc_galaxy_sd_obs_redshift_gauss_init (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->sdz = NULL;
}

static void
_nc_galaxy_sd_obs_redshift_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

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
_nc_galaxy_sd_obs_redshift_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

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
_nc_galaxy_sd_obs_redshift_gauss_dispose (GObject *object)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  nc_galaxy_sd_true_redshift_clear (&self->sdz);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->dispose (object);
}

static void
nc_galaxy_sd_obs_redshift_gauss_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data);
static gdouble _nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data);
static GStrv _nc_galaxy_sd_obs_redshift_gauss_get_header (NcGalaxySDObsRedshift *gsdor);

static void
nc_galaxy_sd_obs_redshift_gauss_class_init (NcGalaxySDObsRedshiftGaussClass *klass)
{
  NcGalaxySDObsRedshiftClass *gsdor_class = NC_GALAXY_SD_OBS_REDSHIFT_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_gauss_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_sd_obs_redshift_gauss_dispose;
  object_class->finalize    = &nc_galaxy_sd_obs_redshift_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian Observed Redshift", "GalaxySDObsRedshiftGauss");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDObsRedshiftGauss:sdz:
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

  /**
   * NcGalaxySDObsRedshiftGauss:sigma:
   *
   * The standard deviation of the gaussian distribution.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SIGMA, "\\sigma", "sigma",
                              0.0, 1.0, 1.0e-3, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_DEFAULT_SIGMA, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen        = &_nc_galaxy_sd_obs_redshift_gauss_gen;
  gsdor_class->integ      = &_nc_galaxy_sd_obs_redshift_gauss_integ;
  gsdor_class->get_header = &_nc_galaxy_sd_obs_redshift_gauss_get_header;
}

#define VECTOR (NCM_MODEL (gsdor))
#define SIGMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SIGMA))

static void
_nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshift *gsdor, NcmRNG *rng, NcmVector *data)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  gdouble sigma                                  = SIGMA;
  gdouble z_min                                  = 0.0;
  gdouble z_max                                  = 0.0;
  gdouble zp;
  gdouble z;

  do {
    z  = nc_galaxy_sd_true_redshift_gen (self->sdz, rng);
    zp = ncm_rng_gaussian_gen (rng, z, sigma * (1.0 + z));
  } while ((zp > z_max) || (zp < z_min));

  ncm_vector_set (data, 0, zp);
  ncm_vector_set (data, 1, sigma);
}

static gdouble
_nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor, gdouble z, NcmVector *data)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  gdouble zp                                     = ncm_vector_get (data, 0);
  gdouble sigma                                  = ncm_vector_get (data, 1);

  return nc_galaxy_sd_true_redshift_integ (self->sdz, z) * exp (-0.5 * (zp - z) * (zp - z) / (sigma * sigma * (1.0 + z) * (1.0 + z))) / (sqrt (2.0 * M_PI) * sigma * (1.0 + z));
}

static GStrv
_nc_galaxy_sd_obs_redshift_gauss_get_header (NcGalaxySDObsRedshift *gsdor)
{
  GStrv header = g_strsplit ("zp zp_sigma", " ", -1);

  return header;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_new:
 * @sdz: a #NcGalaxySDTrueRedshift.
 *
 * Creates a new #NcGalaxySDObsRedshiftGauss object.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftGauss object.
 */
NcGalaxySDObsRedshiftGauss *
nc_galaxy_sd_obs_redshift_gauss_new (NcGalaxySDTrueRedshift *sdz)
{
  return g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_GAUSS, "sdz", sdz, NULL);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_ref:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Increases the reference count of @gsdorgauss by one.
 *
 * Returns: (transfer full): @gsdorgauss.
 */
NcGalaxySDObsRedshiftGauss *
nc_galaxy_sd_obs_redshift_gauss_ref (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  return g_object_ref (gsdorgauss);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_free:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Decreases the reference count of @gsdorgauss by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_free (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  g_object_unref (gsdorgauss);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_clear:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Decreases the reference count of @gsdorgauss by one, and sets the pointer *@gsdorgauss to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_clear (NcGalaxySDObsRedshiftGauss **gsdorgauss)
{
  g_clear_object (gsdorgauss);
}


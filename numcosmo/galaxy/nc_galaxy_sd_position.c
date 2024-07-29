/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position.c
 *
 *  Sat May 20 17:52:48 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION: nc_galaxy_sd_position
 * @title: NcGalaxySDPosition
 * @short_description: Class describing galaxy sample position distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample position distributions.
 * It is composed by two distributions: an image position distribution $P(\theta)$ and a redshift distribution $P(z)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_dtuple.h"
#include <math.h>
#include <gsl/gsl_math.h>


typedef struct _NcGalaxySDPositionPrivate
{
  gint placeholder;
} NcGalaxySDPositionPrivate;

enum
{
  PROP_0,
  PROP_Z_LIM,
  PROP_THETA_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPosition, nc_galaxy_sd_position, NCM_TYPE_MODEL)

static void
nc_galaxy_sd_position_init (NcGalaxySDPosition *gsdp)
{
}

static void
_nc_galaxy_sd_position_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPosition *gsdp = NC_GALAXY_SD_POSITION (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
    {
      NcmDTuple2 *z_lim = g_value_get_boxed (value);

      if (z_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: z_lim is NULL.");

      nc_galaxy_sd_position_set_z_lim (gsdp, z_lim->elements[0], z_lim->elements[1]);
      break;
    }
    case PROP_THETA_LIM:
    {
      NcmDTuple2 *theta_lim = g_value_get_boxed (value);

      if (theta_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: theta_lim is NULL.");

      nc_galaxy_sd_position_set_theta_lim (gsdp, theta_lim->elements[0], theta_lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_position_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPosition *gsdp = NC_GALAXY_SD_POSITION (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
    {
      gdouble z_min, z_max;

      nc_galaxy_sd_position_get_z_lim (gsdp, &z_min, &z_max);

      g_value_take_boxed (value, ncm_dtuple2_new (z_min, z_max));
      break;
    }
    case PROP_THETA_LIM:
    {
      gdouble theta_min, theta_max;

      nc_galaxy_sd_position_get_theta_lim (gsdp, &theta_min, &theta_max);

      g_value_take_boxed (value, ncm_dtuple2_new (theta_min, theta_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_position_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_position_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_position, NC_TYPE_GALAXY_SD_POSITION)

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_sd_position_gen_theta (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen_theta: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_position_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen_z: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, const gdouble theta, const gdouble z)
{
  g_error ("_nc_galaxy_sd_position_integ: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_sd_position_set_z_lim (NcGalaxySDPosition *gsdp, const gdouble z_min, const gdouble z_max)
{
  g_error ("_nc_galaxy_sd_position_set_z_lim: method not implemented.");
}

static void
_nc_galaxy_sd_position_set_theta_lim (NcGalaxySDPosition *gsdp, const gdouble r_min, const gdouble r_max)
{
  g_error ("_nc_galaxy_sd_position_set_theta_lim: method not implemented.");
}

static void
_nc_galaxy_sd_position_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_sd_position_get_z_lim: method not implemented.");
}

static void
_nc_galaxy_sd_position_get_theta_lim (NcGalaxySDPosition *gsdp, gdouble *r_min, gdouble *r_max)
{
  g_error ("_nc_galaxy_sd_position_get_theta_lim: method not implemented.");
}

/* LCOV_LINE_STOP */

static void
nc_galaxy_sd_position_class_init (NcGalaxySDPositionClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_position_set_property;
  model_class->get_property = &_nc_galaxy_sd_position_get_property;
  object_class->finalize    = &_nc_galaxy_sd_position_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample position distribution", "NcGalaxySDPosition");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDPosition:z-lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LIM,
                                   g_param_spec_boxed ("z-lim",
                                                       NULL,
                                                       "Galaxy sample redshift distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPosition:theta-lim:
   *
   * Galaxy sample radius distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_THETA_LIM,
                                   g_param_spec_boxed ("theta-lim",
                                                       NULL,
                                                       "Galaxy sample radius distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class,
                              "NcGalaxySDPosition",
                              "Galaxy sample position distribution.",
                              NULL,
                              TRUE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  klass->gen_theta     = &_nc_galaxy_sd_position_gen_theta;
  klass->gen_z         = &_nc_galaxy_sd_position_gen_z;
  klass->integ         = &_nc_galaxy_sd_position_integ;
  klass->set_z_lim     = &_nc_galaxy_sd_position_set_z_lim;
  klass->set_theta_lim = &_nc_galaxy_sd_position_set_theta_lim;
  klass->get_z_lim     = &_nc_galaxy_sd_position_get_z_lim;
  klass->get_theta_lim = &_nc_galaxy_sd_position_get_theta_lim;
}

/**
 * nc_galaxy_sd_position_ref:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Increase the reference of @gsdp by one.
 *
 * Returns: (transfer full): @gsdp.
 */
NcGalaxySDPosition *
nc_galaxy_sd_position_ref (NcGalaxySDPosition *gsdp)
{
  return g_object_ref (gsdp);
}

/**
 * nc_galaxy_sd_position_free:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Decrease the reference count of @gsdp by one.
 *
 */
void
nc_galaxy_sd_position_free (NcGalaxySDPosition *gsdp)
{
  g_object_unref (gsdp);
}

/**
 * nc_galaxy_sd_position_clear:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Decrease the reference count of @gsdp by one, and sets the pointer *@gsdp to
 * NULL.
 *
 */
void
nc_galaxy_sd_position_clear (NcGalaxySDPosition **gsdp)
{
  g_clear_object (gsdp);
}

/**
 * nc_galaxy_sd_position_set_z_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @z_min: a #gdouble representing the minimum redshift
 * @z_max: a #gdouble representing the maximum redshift
 *
 * Sets the redshift limits of the distribution.
 *
 *
 */
void
nc_galaxy_sd_position_set_z_lim (NcGalaxySDPosition *gsdp, const gdouble z_min, const gdouble z_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_z_lim (gsdp, z_min, z_max);
}

/**
 * nc_galaxy_sd_position_get_z_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @z_min: (out): a #gdouble representing the minimum redshift
 * @z_max: (out): a #gdouble representing the maximum redshift
 *
 * Gets the redshift limits of the distribution.
 *
 *
 */
void
nc_galaxy_sd_position_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_z_lim (gsdp, z_min, z_max);
}

/**
 * nc_galaxy_sd_position_set_theta_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @theta_min: a #gdouble representing the minimum angular radial position
 * @theta_max: a #gdouble representing the maximum angular radial position
 *
 * Sets the radial position limits of the distribution.
 *
 *
 */
void
nc_galaxy_sd_position_set_theta_lim (NcGalaxySDPosition *gsdp, const gdouble theta_min, const gdouble theta_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_theta_lim (gsdp, theta_min, theta_max);
}

/**
 * nc_galaxy_sd_position_get_theta_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @theta_min: (out): a #gdouble representing the minimum angular radial position
 * @theta_max: (out): a #gdouble representing the maximum angular radial position
 *
 * Gets the radial position limits of the distribution.
 *
 *
 */
void
nc_galaxy_sd_position_get_theta_lim (NcGalaxySDPosition *gsdp, gdouble *theta_min, gdouble *theta_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_theta_lim (gsdp, theta_min, theta_max);
}

/**
 * nc_galaxy_sd_position_gen_theta: (virtual gen_theta)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 *
 * Generates a $\theta$ value from the distribution using @rng
 * and saves it in @theta.
 *
 * Returns: the generated $\theta$ value.
 */
gdouble
nc_galaxy_sd_position_gen_theta (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen_theta (gsdp, rng);
}

/**
 * nc_galaxy_sd_position_gen_z: (virtual gen_z)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 *
 * Generates a $z$ value from the distribution using @rng
 * and saves it in @z.
 *
 * Returns: the generated $z$ value.
 */
gdouble
nc_galaxy_sd_position_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen_z (gsdp, rng);
}

/**
 * nc_galaxy_sd_position_integ: (virtual integ)
 * @gsdp: a #NcGalaxySDPosition
 * @theta: a #gdouble representing the angular radial position
 * @z: a #gdouble representing the redshift
 *
 * Computes the probability density of the observables $rtheta and $z$ given the redshift.
 * The probability density is given by $P(z)P(\theta)$.
 *
 * Returns: the probability density at $(\theta, z)$, $P(z)P(\theta)$.
 */
gdouble
nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, const gdouble theta, const gdouble z)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->integ (gsdp, theta, z);
}


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
 * It is composed by two distributions: an image position distribution $P(r)$ and a redshift distribution $P(z)$.
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
  PROP_R_LIM,
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
    case PROP_R_LIM:
    {
      NcmDTuple2 *r_lim = g_value_get_boxed (value);

      if (r_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: r_lim is NULL.");

      nc_galaxy_sd_position_set_r_lim (gsdp, r_lim->elements[0], r_lim->elements[1]);
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
    case PROP_R_LIM:
    {
      gdouble r_min, r_max;

      nc_galaxy_sd_position_get_r_lim (gsdp, &r_min, &r_max);

      g_value_take_boxed (value, ncm_dtuple2_new (r_min, r_max));
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

static gdouble
_nc_galaxy_sd_position_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen_r: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_position_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen_z: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z)
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
_nc_galaxy_sd_position_set_r_lim (NcGalaxySDPosition *gsdp, const gdouble r_min, const gdouble r_max)
{
  g_error ("_nc_galaxy_sd_position_set_r_lim: method not implemented.");
}

static void
_nc_galaxy_sd_position_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_sd_position_get_z_lim: method not implemented.");
}

static void
_nc_galaxy_sd_position_get_r_lim (NcGalaxySDPosition *gsdp, gdouble *r_min, gdouble *r_max)
{
  g_error ("_nc_galaxy_sd_position_get_r_lim: method not implemented.");
}

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
   * NcGalaxySDPosition:r-lim:
   *
   * Galaxy sample radius distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_LIM,
                                   g_param_spec_boxed ("r-lim",
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

  klass->gen_r     = &_nc_galaxy_sd_position_gen_r;
  klass->gen_z     = &_nc_galaxy_sd_position_gen_z;
  klass->integ     = &_nc_galaxy_sd_position_integ;
  klass->set_z_lim = &_nc_galaxy_sd_position_set_z_lim;
  klass->set_r_lim = &_nc_galaxy_sd_position_set_r_lim;
  klass->get_z_lim = &_nc_galaxy_sd_position_get_z_lim;
  klass->get_r_lim = &_nc_galaxy_sd_position_get_r_lim;
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
 * nc_galaxy_sd_position_set_r_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @r_min: a #gdouble representing the minimum radial position
 * @r_max: a #gdouble representing the maximum radial position
 *
 * Sets the radial position limits of the distribution.
 *
 *
 */
void
nc_galaxy_sd_position_set_r_lim (NcGalaxySDPosition *gsdp, const gdouble r_min, const gdouble r_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_r_lim (gsdp, r_min, r_max);
}

/**
 * nc_galaxy_sd_position_get_r_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @r_min: (out): a #gdouble representing the minimum radial position
 * @r_max: (out): a #gdouble representing the maximum radial position
 *
 * Gets the radial position limits of the distribution.
 *
 *
 */
void
nc_galaxy_sd_position_get_r_lim (NcGalaxySDPosition *gsdp, gdouble *r_min, gdouble *r_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_r_lim (gsdp, r_min, r_max);
}

/**
 * nc_galaxy_sd_position_gen_r: (virtual gen_r)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 *
 * Generates a $r$ value from the distribution using @rng
 * and saves it in @r.
 *
 * Returns: the generated $r$ value.
 */
gdouble
nc_galaxy_sd_position_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen_r (gsdp, rng);
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
 * @r: a #gdouble representing the radial position
 * @z: a #gdouble representing the redshift
 *
 * Computes the probability density of the observables $r$ and $z$ given the redshift.
 * The probability density is given by $P(z)P(r)$.
 *
 * Returns: the probability density at $(r, z)$, $P(z)P(r)$.
 */
gdouble
nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->integ (gsdp, r, z);
}


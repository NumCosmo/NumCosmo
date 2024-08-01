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
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_dtuple.h"
#include "math.h"
#include "gsl/gsl_math.h"


typedef struct _NcGalaxySDPositionPrivate
{
  gint placeholder;
} NcGalaxySDPositionPrivate;

enum
{
  PROP_0,
  PROP_RA_LIM,
  PROP_DEC_LIM,
  PROP_HEADER,
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
    case PROP_RA_LIM:
    {
      NcmDTuple2 *ra_lim = g_value_get_boxed (value);

      if (ra_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: ra_lim is NULL.");

      nc_galaxy_sd_position_set_ra_lim (gsdp, ra_lim->elements[0], ra_lim->elements[1]);
      break;
    }
    case PROP_DEC_LIM:
    {
      NcmDTuple2 *dec_lim = g_value_get_boxed (value);

      if (dec_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: dec_lim is NULL.");

      nc_galaxy_sd_position_set_dec_lim (gsdp, dec_lim->elements[0], dec_lim->elements[1]);
      break;
    }
    case PROP_HEADER:
    {
      nc_galaxy_sd_position_set_header (gsdp, g_value_get_string (value));
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
    case PROP_RA_LIM:
    {
      gdouble ra_min, ra_max;

      nc_galaxy_sd_position_get_ra_lim (gsdp, &ra_min, &ra_max);

      g_value_take_boxed (value, ncm_dtuple2_new (ra_min, ra_max));
      break;
    }
    case PROP_DEC_LIM:
    {
      gdouble dec_min, dec_max;

      nc_galaxy_sd_position_get_dec_lim (gsdp, &dec_min, &dec_max);

      g_value_take_boxed (value, ncm_dtuple2_new (dec_min, dec_max));
      break;
    }
    case PROP_HEADER:
    {
      g_value_set_string (value, nc_galaxy_sd_position_get_header (gsdp));
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
_nc_galaxy_sd_position_gen_ra (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen_ra: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_position_gen_dec (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen_dec: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_position_integ: method not implemented.");

  return 0.0;
}

static gboolean
_nc_galaxy_sd_position_set_ra_lim (NcGalaxySDPosition *gsdp, const gdouble ra_min, const gdouble ra_max)
{
  g_error ("_nc_galaxy_sd_position_set_ra_lim: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_set_dec_lim (NcGalaxySDPosition *gsdp, const gdouble dec_min, const gdouble dec_max)
{
  g_error ("_nc_galaxy_sd_position_set_dec_lim: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_set_header (NcGalaxySDPosition *gsdp, const GStrv header)
{
  g_error ("_nc_galaxy_sd_position_set_header: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max)
{
  g_error ("_nc_galaxy_sd_position_get_ra_lim: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max)
{
  g_error ("_nc_galaxy_sd_position_get_dec_lim: method not implemented.");

  return FALSE;
}

GStrv
_nc_galaxy_sd_position_get_header (NcGalaxySDPosition *gsdp)
{
  g_error ("_nc_galaxy_sd_position_get_header: method not implemented.");

  return NULL;
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
   * NcGalaxySDPosition:ra-lim:
   *
   * Galaxy sample righ ascension distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RA_LIM,
                                   g_param_spec_boxed ("ra-lim",
                                                       NULL,
                                                       "Galaxy sample right ascension distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPosition:dec-lim:
   *
   * Galaxy sample declination distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEC_LIM,
                                   g_param_spec_boxed ("dec-lim",
                                                       NULL,
                                                       "Galaxy sample declination distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPosition:header:
   *
   * Galaxy sample position distribution header.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_HEADER,
                                   g_param_spec_string ("header",
                                                        NULL,
                                                        "Galaxy sample position distribution header",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class,
                              "NcGalaxySDPosition",
                              "Galaxy sample position distribution.",
                              NULL,
                              TRUE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  klass->gen_ra      = &_nc_galaxy_sd_position_gen_ra;
  klass->gen_dec     = &_nc_galaxy_sd_position_gen_dec;
  klass->integ       = &_nc_galaxy_sd_position_integ;
  klass->set_ra_lim  = &_nc_galaxy_sd_position_set_ra_lim;
  klass->set_dec_lim = &_nc_galaxy_sd_position_set_dec_lim;
  klass->set_header  = &_nc_galaxy_sd_position_set_header;
  klass->get_ra_lim  = &_nc_galaxy_sd_position_get_ra_lim;
  klass->get_dec_lim = &_nc_galaxy_sd_position_get_dec_lim;
  klass->get_header  = &_nc_galaxy_sd_position_get_header;
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
 * nc_galaxy_sd_position_set_ra_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @ra_min: the minimum right ascension
 * @ra_max: the maximum right ascension
 *
 * Sets the right ascension limits of the distribution.
 *
 * Returns: TRUE if the limits were set successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_set_ra_lim (NcGalaxySDPosition *gsdp, const gdouble ra_min, const gdouble ra_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_ra_lim (gsdp, ra_min, ra_max);
}

/**
 * nc_galaxy_sd_position_set_dec_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @dec_min: the minimum declination
 * @dec_max: the maximum declination
 *
 * Sets the declination limits of the distribution.
 *
 * Returns: TRUE if the limits were set successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_set_dec_lim (NcGalaxySDPosition *gsdp, const gdouble dec_min, const gdouble dec_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_dec_lim (gsdp, dec_min, dec_max);
}

/**
 * nc_galaxy_sd_position_set_header:
 * @gsdp: a #NcGalaxySDPosition
 * @header: a #GStrv
 *
 * Sets the header of the distribution.
 *
 * Returns: TRUE if the header was set successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_set_header (NcGalaxySDPosition *gsdp, const GStrv header)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_header (gsdp, header);
}

/**
 * nc_galaxy_sd_position_get_ra_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @ra_min: (out): the minimum right ascension
 * @ra_max: (out): the maximum right ascension
 *
 * Gets the right ascension limits of the distribution.
 *
 * Returns: TRUE if the limits were retrieved successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_ra_lim (gsdp, ra_min, ra_max);
}

/**
 * nc_galaxy_sd_position_get_dec_lim:
 * @gsdp: a #NcGalaxySDPosition
 * @dec_min: (out): the minimum declination
 * @dec_max: (out): the maximum declination
 *
 * Gets the declination limits of the distribution.
 *
 * Returns: TRUE if the limits were retrieved successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_dec_lim (gsdp, dec_min, dec_max);
}

/**
 * nc_galaxy_sd_position_get_header:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Gets the header of the distribution.
 *
 * Returns: the header of the distribution.
 */
GStrv
nc_galaxy_sd_position_get_header (NcGalaxySDPosition *gsdp)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_header (gsdp);
}

/**
 * nc_galaxy_sd_position_gen_ra: (virtual gen_ra)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 *
 * Generates a right ascension value from the distribution using @rng.
 *
 * Returns: the generated right ascension value.
 */
gdouble
nc_galaxy_sd_position_gen_ra (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen_ra (gsdp, rng);
}

/**
 * nc_galaxy_sd_position_gen_dec: (virtual gen_dec)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 *
 * Generates a declination value from the distribution using @rng
 *
 * Returns: the generated declination value.
 */
gdouble
nc_galaxy_sd_position_gen_dec (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen_dec (gsdp, rng);
}

/**
 * nc_galaxy_sd_position_integ: (virtual integ)
 * @gsdp: a #NcGalaxySDPosition
 * @data: a #NcmVector
 *
 * Computes the probability density of the right ascension and declination.
 *
 * Returns: the probability density at $(\text{ra}, \text{dec})$, $P(\text{ra})P(\text{dec})$.
 */
gdouble
nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, NcmVector *data)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->integ (gsdp, data);
}


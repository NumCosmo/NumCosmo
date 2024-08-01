/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_redshift.c
 *
 *  Wed Jul 31 20:52:43 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_redshift.c
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
 * SECTION: nc_galaxy_sd_redshift
 * @title: NcGalaxySDRedshift
 * @short_description: Class describing galaxy sample redshift distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample redshift distributions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_redshift.h"
#include "math/ncm_dtuple.h"

typedef struct _NcGalaxySDRedshiftPrivate
{
  gint placeholder;
} NcGalaxySDRedshiftPrivate;

enum
{
  PROP_0,
  PROP_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDRedshift, nc_galaxy_sd_redshift, NCM_TYPE_MODEL)

static void
nc_galaxy_sd_redshift_init (NcGalaxySDRedshift *gsdr)
{
}

static void
_nc_galaxy_sd_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDRedshift *gsdr = NC_GALAXY_SD_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_REDSHIFT (gsdr));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_sd_redshift_set_property: lim is NULL");

      nc_galaxy_sd_redshift_set_lim (gsdr, lim->elements[0], lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDRedshift *gsdr = NC_GALAXY_SD_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_REDSHIFT (gsdr));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      gdouble z_min, z_max;

      nc_galaxy_sd_redshift_get_lim (gsdr, &z_min, &z_max);

      g_value_set_boxed (value, ncm_dtuple2_new (z_min, z_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_redshift_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_redshift_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_redshift, NC_TYPE_GALAXY_SD_REDSHIFT)

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_sd_redshift_gen (NcGalaxySDRedshift *gsdr, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_redshift_gen_z: not implemented");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_redshift_integ (NcGalaxySDRedshift *gsdr, const gdouble z)
{
  g_error ("_nc_galaxy_sd_redshift_integ: not implemented");

  return 0.0;
}

static void
_nc_galaxy_sd_redshift_set_lim (NcGalaxySDRedshift *gsdr, const gdouble z_min, const gdouble z_max)
{
  g_error ("_nc_galaxy_sd_redshift_set_lim: not implemented");
}

static void
_nc_galaxy_sd_redshift_get_lim (NcGalaxySDRedshift *gsdr, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_sd_redshift_get_lim: not implemented");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_redshift_class_init (NcGalaxySDRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_redshift_set_property;
  model_class->get_property = &_nc_galaxy_sd_redshift_get_property;
  object_class->finalize    = &_nc_galaxy_sd_redshift_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample redshift distribution", "GalaxySDRedshift");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDRedshift:lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LIM,
                                   g_param_spec_boxed ("lim",
                                                       NULL,
                                                       "Galaxy sample redshift distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NCM_MODEL_CLASS_CHECK_PARAMS_INFO (NCM_MODEL_CLASS (klass));

  klass->gen     = &_nc_galaxy_sd_redshift_gen;
  klass->integ   = &_nc_galaxy_sd_redshift_integ;
  klass->set_lim = &_nc_galaxy_sd_redshift_set_lim;
  klass->get_lim = &_nc_galaxy_sd_redshift_get_lim;
}

/**
 * nc_galaxy_sd_redshift_ref:
 * @gsdr: a #NcGalaxySDRedshift
 *
 * Increases the reference count of @gsdr by one.
 *
 * Returns: (transfer full): @gsdr.
 */
NcGalaxySDRedshift *
nc_galaxy_sd_redshift_ref (NcGalaxySDRedshift *gsdr)
{
  return g_object_ref (gsdr);
}

/**
 * nc_galaxy_sd_redshift_free:
 * @gsdr: a #NcGalaxySDRedshift
 *
 * Decreases the reference count of @gsdr by one.
 *
 */
void
nc_galaxy_sd_redshift_free (NcGalaxySDRedshift *gsdr)
{
  g_object_unref (gsdr);
}

/**
 * nc_galaxy_sd_redshift_clear:
 * @gsdr: a #NcGalaxySDRedshift
 *
 * Decreases the reference count of @gsdr by one, and sets the
 * pointer @gsdr to NULL.
 *
 */
void
nc_galaxy_sd_redshift_clear (NcGalaxySDRedshift **gsdr)
{
  g_clear_object (gsdr);
}

/**
 * nc_galaxy_sd_redshift_set_lim:
 * @gsdr: a #NcGalaxySDRedshift
 * @z_min: a #gdouble representing minimum redshift
 * @z_max: a #gdouble representing maximum redshift
 *
 * Sets the redshift limits of the galaxy sample redshift distribution.
 *
 */
void
nc_galaxy_sd_redshift_set_lim (NcGalaxySDRedshift *gsdr, const gdouble z_min, const gdouble z_max)
{
  NC_GALAXY_SD_REDSHIFT_GET_CLASS (gsdr)->set_lim (gsdr, z_min, z_max);
}

/**
 * nc_galaxy_sd_redshift_get_lim:
 * @gsdr: a #NcGalaxySDRedshift
 * @z_min: a #gdouble pointer representing minimum redshift
 * @z_max: a #gdouble pointer representing maximum redshift
 *
 * Gets the redshift limits of the galaxy sample redshift distribution.
 *
 */
void
nc_galaxy_sd_redshift_get_lim (NcGalaxySDRedshift *gsdr, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_SD_REDSHIFT_GET_CLASS (gsdr)->get_lim (gsdr, z_min, z_max);
}

/**
 * nc_galaxy_sd_redshift_gen:
 * @gsdr: a #NcGalaxySDRedshift
 * @rng: a #NcmRNG
 *
 * Generates a redshift value from the galaxy sample redshift distribution.
 *
 * Returns: the generated redshift.
 */
gdouble
nc_galaxy_sd_redshift_gen (NcGalaxySDRedshift *gsdr, NcmRNG *rng)
{
  return NC_GALAXY_SD_REDSHIFT_GET_CLASS (gsdr)->gen (gsdr, rng);
}

/**
 * nc_galaxy_sd_redshift_integ:
 * @gsdr: a #NcGalaxySDRedshift
 * @z: a #gdouble representing the redshift
 *
 * Evaluates the galaxy sample redshift distribution at redshift @z.
 *
 * Returns: the probability density at $z$, $P(z)$.
 */
gdouble
nc_galaxy_sd_redshift_integ (NcGalaxySDRedshift *gsdr, const gdouble z)
{
  return NC_GALAXY_SD_REDSHIFT_GET_CLASS (gsdr)->integ (gsdr, z);
}


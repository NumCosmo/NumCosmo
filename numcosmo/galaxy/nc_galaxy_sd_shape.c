/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape.c
 *
 *  Sat May 21 20:43:32 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape.c
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
 * SECTION: nc_galaxy_sd_shape
 * @title: NcGalaxySDShape
 * @short_description: Class describing galaxy sample shape distribution.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample shape distribution.
 * It is composed by a distribution $P(s)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "nc_enum_types.h"
#include "nc_hicosmo.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"

typedef struct _NcGalaxySDShapePrivate
{
  gint placeholder;
} NcGalaxySDShapePrivate;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShape, nc_galaxy_sd_shape, NCM_TYPE_MODEL);

static void
nc_galaxy_sd_shape_init (NcGalaxySDShape *gsds)
{
}

static void
_nc_galaxy_sd_shape_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE (gsds));

  switch (property_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE (gsds));

  switch (property_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static gboolean
_nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmRNG *rng, const gdouble ra, const gdouble dec, const gdouble z, gdouble *e1, gdouble *e2, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_shape_gen: method not implemented.");

  return FALSE;
}

static gdouble
_nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_shape_integ: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_sd_shape_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_shape_integ_optzs_prep: method not implemented.");
}

static gdouble
_nc_galaxy_sd_shape_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  g_error ("_nc_galaxy_sd_shape_integ_optzs: method not implemented.");

  return 0.0;
}

static GStrv
_nc_galaxy_sd_shape_get_header (NcGalaxySDShape *gsds)
{
  g_error ("_nc_galaxy_sd_shape_get_header: method not implemented.");

  return NULL;
}

static gboolean
_nc_galaxy_sd_shape_set_coord (NcGalaxySDShape *gsds, NcGalaxyWLObsCoord coord)
{
  g_error ("_nc_galaxy_sd_shape_set_coord: method not implemented.");

  return FALSE;
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_shape_class_init (NcGalaxySDShapeClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_get_property;
  object_class->finalize    = &_nc_galaxy_sd_shape_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample shape distribution", "GalaxySDShape");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  ncm_model_class_check_params_info (model_class);

  klass->gen              = &_nc_galaxy_sd_shape_gen;
  klass->integ            = &_nc_galaxy_sd_shape_integ;
  klass->integ_optzs_prep = &_nc_galaxy_sd_shape_integ_optzs_prep;
  klass->integ_optzs      = &_nc_galaxy_sd_shape_integ_optzs;
  klass->get_header       = &_nc_galaxy_sd_shape_get_header;
  klass->set_coord        = &_nc_galaxy_sd_shape_set_coord;
}

/**
 * nc_galaxy_sd_shape_ref:
 * @gsds: a #NcGalaxySDShape
 *
 * Increases the reference count of @gsds by one.
 *
 * Returns: (transfer full): @gsds.
 */
NcGalaxySDShape *
nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds)
{
  return g_object_ref (gsds);
}

/**
 * nc_galaxy_sd_shape_free:
 * @gsds: a #NcGalaxySDShape
 *
 * Decreases the reference count of @gsds by one.
 *
 */
void
nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds)
{
  g_object_unref (gsds);
}

/**
 * nc_galaxy_sd_shape_clear:
 * @gsds: a #NcGalaxySDShape
 *
 * Decreases the reference count of @gsds by one, and sets the pointer *@gsds to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds)
{
  g_clear_object (gsds);
}

/**
 * nc_galaxy_sd_shape_gen: (virtual gen)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @rng: a #NcmRNG
 * @ra: the right ascension of the galaxy
 * @dec: the declination of the galaxy
 * @z: the redshift of the galaxy
 * @e1: (out): the generated first ellipticity component $e_1$
 * @e2: (out): the generated second ellipticity component $e_2$
 * @data: a #NcmVector
 *
 * Generates a shape value from the position using @rng.
 *
 */
gboolean
nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmRNG *rng, const gdouble ra, const gdouble dec, const gdouble z, gdouble *e1, gdouble *e2, NcmVector *data)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->gen (gsds, cosmo, dp, smd, rng, ra, dec, z, e1, e2, data);

  return FALSE;
}

/**
 * nc_galaxy_sd_shape_integ: (virtual integ)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @data: a #NcmVector
 *
 * Computes the probability density of the observable shape given the position.
 *
 * Returns: the probability density of observable shape, $P(s)$.
 */
gdouble
nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ (gsds, cosmo, dp, smd, data);
}

/**
 * nc_galaxy_sd_shape_integ_optzs_prep: (virtual integ_optzs_prep)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @data: a #NcmVector
 *
 * Prepares $gsds$ to compute the probability density of the observable shape given the position.
 *
 */
void
nc_galaxy_sd_shape_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ_optzs_prep (gsds, cosmo, dp, smd, data);
}

/**
 * nc_galaxy_sd_shape_integ_optzs: (virtual integ_optzs)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @data: a #NcmVector
 *
 * Computes the probability density of the observable shape given the position.
 *
 * Returns: the probability density of observable shape, $P(s)$.
 */
gdouble
nc_galaxy_sd_shape_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ_optzs (gsds, cosmo, dp, smd, data);
}

/**
 * nc_galaxy_sd_shape_get_header: (virtual get_header)
 * @gsds: a #NcGalaxySDShape
 *
 * Gets the header of the distribution.
 *
 * Returns: (transfer full): the header of the data.
 */
GStrv
nc_galaxy_sd_shape_get_header (NcGalaxySDShape *gsds)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->get_header (gsds);
}

/**
 * nc_galaxy_sd_shape_set_coord: (virtual set_coord)
 * @gsds: a #NcGalaxySDShape
 * @coord: a #NcGalaxyWLObsCoord
 *
 * Sets the coordinate of the galaxy.
 *
 * Returns: TRUE if the coordinate was set, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_shape_set_coord (NcGalaxySDShape *gsds, NcGalaxyWLObsCoord coord)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->set_coord (gsds, coord);
}


/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape_gauss.c
 *
 *  Mon June 5 14:56:41 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_gauss.c
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
 * SECTION:nc_galaxy_sd_shape_gauss
 * @title: NcGalaxySDShapeGauss
 * @short_description: Class describing galaxy sample shape gaussian distribution
 * @stability: Unstable
 *
 * This class describes a galaxy sample shape gaussian
 * probability distribution $P(s)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape_gauss.h"
#include "galaxy/nc_galaxy_sd_shape.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */


typedef struct _NcGalaxySDShapeGaussPrivate
{
  NcDistance *dist;
  NcGalaxyWLObsCoord coord;
  NcWLSurfaceMassDensityOptzs optzs;
  gdouble phi;
} NcGalaxySDShapeGaussPrivate;

struct _NcGalaxySDShapeGauss
{
  NcGalaxySDShape parent_instance;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_COORD,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShapeGauss, nc_galaxy_sd_shape_gauss, NC_TYPE_GALAXY_SD_SHAPE);

static void
nc_galaxy_sd_shape_gauss_init (NcGalaxySDShapeGauss *gsdsgauss)
{
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  self->dist = NULL;
  self->phi  = 0.0;
}

static void
_nc_galaxy_sd_shape_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_private (gsdsgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsgauss));

  switch (prop_id)
  {
    case PROP_DIST:
      self->dist = g_value_get_object (value);
      break;
    case PROP_COORD:
      self->coord = g_value_get_enum (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_private (gsdsgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsgauss));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
    case PROP_COORD:
      g_value_set_enum (value, self->coord);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_gauss_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_shape_gauss_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmRNG *rng, const gdouble ra, const gdouble dec, const gdouble z, gdouble *e1, gdouble *e2, NcmVector *data);
static gdouble _nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data);
static void _nc_galaxy_sd_shape_gauss_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data);
static gdouble _nc_galaxy_sd_shape_gauss_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data);
static GStrv _nc_galaxy_sd_shape_gauss_get_header (NcGalaxySDShape *gsds);
static gboolean _nc_galaxy_sd_shape_gauss_set_coord (NcGalaxySDShape *gsds, NcGalaxyWLObsCoord coord);

static void
nc_galaxy_sd_shape_gauss_class_init (NcGalaxySDShapeGaussClass *klass)
{
  NcGalaxySDShapeClass *sd_position_class = NC_GALAXY_SD_SHAPE_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  object_class->set_property = &_nc_galaxy_sd_shape_gauss_set_property;
  object_class->get_property = &_nc_galaxy_sd_shape_gauss_get_property;
  object_class->dispose      = &_nc_galaxy_sd_shape_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_sd_shape_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian galaxy shape distribution", "Gaussian shape");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_SHAPE_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDShapeGauss:sigma:
   *
   * The standard deviation of the gaussian distribution.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_SHAPE_GAUSS_SIGMA, "\\sigma", "sigma", 1e-8, 2.0, 1e-1, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_SIGMA, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  sd_position_class->gen              = &_nc_galaxy_sd_shape_gauss_gen;
  sd_position_class->integ            = &_nc_galaxy_sd_shape_gauss_integ;
  sd_position_class->integ_optzs_prep = &_nc_galaxy_sd_shape_gauss_integ_optzs_prep;
  sd_position_class->integ_optzs      = &_nc_galaxy_sd_shape_gauss_integ_optzs;
  sd_position_class->get_header       = &_nc_galaxy_sd_shape_gauss_get_header;
  sd_position_class->set_coord        = &_nc_galaxy_sd_shape_gauss_set_coord;
}

#define VECTOR (NCM_MODEL (gsds))
#define SIGMA  (ncm_model_get_param (VECTOR, NC_GALAXY_SD_SHAPE_GAUSS_SIGMA))

static void
_nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmRNG *rng, const gdouble ra, const gdouble dec, const gdouble z, gdouble *e1, gdouble *e2, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  gdouble ra_cl      = ncm_vector_get (data, 0);
  gdouble dec_cl     = ncm_vector_get (data, 1);
  gdouble z_cl       = ncm_vector_get (data, 2);
  gdouble theta      = 0.0;
  gdouble phi        = 0.0;
  gdouble et_s       = ncm_rng_gaussian_gen (rng, 0.0, SIGMA);
  gdouble ex_s       = ncm_rng_gaussian_gen (rng, 0.0, SIGMA);
  complex double e_s = et_s + I * ex_s;
  complex double e_o = e_s;
  complex double gt  = 0.0;
  gdouble r;

  ncm_util_polar_angles (ra_cl, dec_cl, ra, dec, &theta, &phi);

  if (self->coord == NC_GALAXY_WL_OBS_COORD_CELESTIAL)
    phi = M_PI - phi;

  r = ncm_util_projected_radius (theta, nc_distance_angular_diameter (self->dist, cosmo, z_cl) * nc_distance_hubble (self->dist, cosmo));

  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cl, z_cl);

    if (cabs (gt) > 1.0)
      e_o = (1.0 + gt * conj (e_s)) / (conj (e_s) + conj (gt));
    else
      e_o = (e_s + gt) / (1.0 + conj (gt) * e_s);
  }

  *e1 = -creal (e_o) * cos (2 * phi) + cimag (e_o) * sin (2 * phi);
  *e2 = -creal (e_o) * sin (2 * phi) - cimag (e_o) * cos (2 * phi);
}

static gdouble
_nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble ra_cl                            = ncm_vector_get (data, 0);
  gdouble dec_cl                           = ncm_vector_get (data, 1);
  gdouble z_cl                             = ncm_vector_get (data, 2);
  gdouble ra                               = ncm_vector_get (data, 3);
  gdouble dec                              = ncm_vector_get (data, 4);
  gdouble z                                = ncm_vector_get (data, 5);
  gdouble e1                               = ncm_vector_get (data, 6);
  gdouble e1_sigma                         = ncm_vector_get (data, 7);
  gdouble e2                               = ncm_vector_get (data, 8);
  gdouble e2_sigma                         = ncm_vector_get (data, 9);
  complex double e_o                       = e1 + I * e2;
  complex double e_s                       = e_o;
  complex double g                         = 0.0;
  gdouble theta                            = 0.0;
  gdouble phi                              = 0.0;
  gdouble gt;
  gdouble r;

  ncm_util_polar_angles (ra_cl, dec_cl, ra, dec, &theta, &phi);

  if (self->coord == NC_GALAXY_WL_OBS_COORD_CELESTIAL)
    phi = M_PI - phi;

  r = ncm_util_projected_radius (theta, nc_distance_angular_diameter (self->dist, cosmo, z_cl) * nc_distance_hubble (self->dist, cosmo));

  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cl, z_cl);

    g = -(creal (gt) * cos (2 * phi) - cimag (gt) * sin (2 * phi)) - I * (creal (gt) * sin (2 * phi) - cimag (gt) * cos (2 * phi));

    if (cabs (gt) > 1.0)
      e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
    else
      e_s = (e_o - g) / (1.0 - conj (g) * e_o);
  }

  return gsl_ran_gaussian_pdf (creal (e_s), e1_sigma) * gsl_ran_gaussian_pdf (cimag (e_s), e2_sigma);
}

static void
_nc_galaxy_sd_shape_gauss_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble ra_cl                            = ncm_vector_get (data, 0);
  gdouble dec_cl                           = ncm_vector_get (data, 1);
  gdouble z_cl                             = ncm_vector_get (data, 2);
  gdouble ra                               = ncm_vector_get (data, 3);
  gdouble dec                              = ncm_vector_get (data, 4);
  gdouble theta                            = 0.0;
  gdouble phi                              = 0.0;
  gdouble r;

  ncm_util_polar_angles (ra_cl, dec_cl, ra, dec, &theta, &phi);

  r = ncm_util_projected_radius (theta, nc_distance_angular_diameter (self->dist, cosmo, z_cl) * nc_distance_hubble (self->dist, cosmo));

  if (self->coord == NC_GALAXY_WL_OBS_COORD_CELESTIAL)
    phi = M_PI - phi;

  self->phi = phi;

  nc_wl_surface_mass_density_reduced_shear_optzs_prep (smd, dp, cosmo, r, z_cl, z_cl, &self->optzs);
}

static gdouble
_nc_galaxy_sd_shape_gauss_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble z_cl                             = ncm_vector_get (data, 2);
  gdouble z                                = ncm_vector_get (data, 5);
  gdouble e1                               = ncm_vector_get (data, 6);
  gdouble e1_sigma                         = ncm_vector_get (data, 7);
  gdouble e2                               = ncm_vector_get (data, 8);
  gdouble e2_sigma                         = ncm_vector_get (data, 9);
  complex double e_o                       = e1 + I * e2;
  complex double e_s                       = e_s;
  complex double g                         = 0.0;
  gdouble gt;
  gdouble r;

  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear_optzs (smd, dp, cosmo, z, z_cl, &self->optzs);

    g = -(creal (gt) * cos (2 * self->phi) - cimag (gt) * sin (2 * self->phi)) - I * (creal (gt) * sin (2 * self->phi) - cimag (gt) * cos (2 * self->phi));

    if (cabs (gt) > 1.0)
      e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
    else
      e_s = (e_o - g) / (1.0 - conj (g) * e_o);
  }

  return gsl_ran_gaussian_pdf (creal (e_s), e1_sigma) * gsl_ran_gaussian_pdf (cimag (e_s), e2_sigma);
}

static GStrv
_nc_galaxy_sd_shape_gauss_get_header (NcGalaxySDShape *gsds)
{
  GStrv header = g_strsplit ("ra_cl dec_cl z_cl ra dec zp e1 e1_sigma e2 e2_sigma", " ", -1);

  return header;
}

static gboolean
_nc_galaxy_sd_shape_gauss_set_coord (NcGalaxySDShape *gsds, NcGalaxyWLObsCoord coord)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  self->coord = coord;

  return TRUE;
}

/**
 * nc_galaxy_sd_shape_gauss_new:
 *
 * Creates a new #NcGalaxySDShapeGauss
 *
 * Returns: (transfer full): a new NcGalaxySDShapeGauss.
 */
NcGalaxySDShapeGauss *
nc_galaxy_sd_shape_gauss_new (NcDistance *dist)
{
  NcGalaxySDShapeGauss *gsdsgauss = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_GAUSS,
                                                  NULL);

  return gsdsgauss;
}

/**
 * nc_galaxy_sd_shape_gauss_ref:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Increase the reference of @gsdsgauss by one.
 *
 * Returns: (transfer full): @gsdsgauss.
 */
NcGalaxySDShapeGauss *
nc_galaxy_sd_shape_gauss_ref (NcGalaxySDShapeGauss *gsdsgauss)
{
  return g_object_ref (gsdsgauss);
}

/**
 * nc_galaxy_sd_shape_gauss_free:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Decrease the reference count of @gsdsgauss by one.
 *
 */
void
nc_galaxy_sd_shape_gauss_free (NcGalaxySDShapeGauss *gsdsgauss)
{
  g_object_unref (gsdsgauss);
}

/**
 * nc_galaxy_sd_shape_gauss_clear:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Decrease the reference count of @gsdsgauss by one, and sets the pointer *@gsdsgauss to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_gauss_clear (NcGalaxySDShapeGauss **gsdsgauss)
{
  g_clear_object (gsdsgauss);
}


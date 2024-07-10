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
 * @short_description: Class describing galaxy sample shape distribution with gaussian distribution.
 * @stability: Unstable
 *
 *
 * his class describes a galaxy sample shape gaussian
 * probability distribution $P(s)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "galaxy/nc_galaxy_sd_shape_gauss.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include <math.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */


struct _NcGalaxySDShapeGaussPrivate
{
  NcWLSurfaceMassDensityOptzs optzs;
  gdouble sigma;
};

enum
{
  PROP_0,
  PROP_SIGMA,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShapeGauss, nc_galaxy_sd_shape_gauss, NC_TYPE_GALAXY_SD_SHAPE);

static void
nc_galaxy_sd_shape_gauss_init (NcGalaxySDShapeGauss *gsdsgauss)
{
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  self->sigma = 0.0;
}

static void
_nc_galaxy_sd_shape_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss = NC_GALAXY_SD_SHAPE_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (object));

  switch (prop_id)
  {
    case PROP_SIGMA:
      nc_galaxy_sd_shape_gauss_set_sigma (gsdsgauss, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_shape_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss = NC_GALAXY_SD_SHAPE_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (object));

  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;

  switch (prop_id)
  {
    case PROP_SIGMA:
      g_value_set_double (value, self->sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_shape_gauss_dispose (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_shape_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, NcmRNG *rng, gdouble theta, gdouble z, gdouble *et, gdouble *ex);
static gdouble _nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble theta, const gdouble z, const gdouble et, const gdouble ex);
static void _nc_galaxy_sd_shape_gauss_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble theta);
static gdouble _nc_galaxy_sd_shape_gauss_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z_source, const gdouble et, const gdouble ex);

static void
nc_galaxy_sd_shape_gauss_class_init (NcGalaxySDShapeGaussClass *klass)
{
  NcGalaxySDShapeClass *sd_position_class = NC_GALAXY_SD_SHAPE_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_sd_shape_gauss_set_property;
  object_class->get_property = &_nc_galaxy_sd_shape_gauss_get_property;
  object_class->dispose      = &_nc_galaxy_sd_shape_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_sd_shape_gauss_finalize;

  /**
   * NcGalaxySDShapeGauss:sigma:
   *
   * Galaxy sample shape distribution standard deviation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA,
                                   g_param_spec_double ("sigma",
                                                        NULL,
                                                        "Galaxy sample shape standard deviation",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd_position_class->gen              = &_nc_galaxy_sd_shape_gauss_gen;
  sd_position_class->integ            = &_nc_galaxy_sd_shape_gauss_integ;
  sd_position_class->integ_optzs_prep = &_nc_galaxy_sd_shape_gauss_integ_optzs_prep;
  sd_position_class->integ_optzs      = &_nc_galaxy_sd_shape_gauss_integ_optzs;
}

static void
_nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, NcmRNG *rng, const gdouble theta, const gdouble z, gdouble *et, gdouble *ex)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;

  gdouble et_source       = ncm_rng_gaussian_gen (rng, 0, self->sigma);
  gdouble ex_source       = ncm_rng_gaussian_gen (rng, 0, self->sigma);
  complex double e_source = et_source + I * ex_source;
  complex double e_obs    = e_source;
  complex double gt       = 0.0;

  if (z > z_cluster)
  {
    gt    = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, theta, z, z_cluster, z_cluster);
    e_obs = (e_source + gt) / (1.0 + conj (gt) * e_source);
  }

  *et = creal (e_obs);
  *ex = cimag (e_obs);
}

static gdouble
_nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble theta, const gdouble z, const gdouble et, const gdouble ex)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;
  complex double e_obs                     = et + I * ex;
  complex double e_source                  = e_obs;
  complex double gt                        = 0.0;

  if (z > z_cluster)
  {
    gt = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, theta, z, z_cluster, z_cluster);

    if (cabs (gt) > 1.0)
      e_source = (1.0 - gt * conj (e_obs)) / (conj (e_obs) - conj (gt));
    else
      e_source = (e_obs - gt) / (1.0 - conj (gt) * e_obs);
  }

  return exp (-0.5 * gsl_pow_2 (cabs (e_source) / self->sigma)) / (2.0 * M_PI * self->sigma * self->sigma);
}

static void
_nc_galaxy_sd_shape_gauss_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble theta)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;

  nc_wl_surface_mass_density_reduced_shear_optzs_prep (smd, dp, cosmo, theta, z_cluster, z_cluster, &self->optzs);
}

static gdouble
_nc_galaxy_sd_shape_gauss_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z_source, const gdouble et, const gdouble ex)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;
  complex double e_obs                     = et + I * ex;
  complex double e_source                  = e_obs;
  complex double gt                        = 0.0;

  if (z_source > z_cluster)
  {
    gt = nc_wl_surface_mass_density_reduced_shear_optzs (smd, dp, cosmo, z_source, z_cluster, &self->optzs);

    if (cabs (gt) > 1.0)
      e_source = (1.0 - gt * conj (e_obs)) / (conj (e_obs) - conj (gt));
    else
      e_source = (e_obs - gt) / (1.0 - conj (gt) * e_obs);
  }

  return exp (-0.5 * gsl_pow_2 (cabs (e_source) / self->sigma)) / (2.0 * M_PI * self->sigma * self->sigma);
}

/**
 * nc_galaxy_sd_shape_gauss_new:
 *
 * Creates a new #NcGalaxySDShapeGauss
 *
 * Returns: (transfer full): a new NcGalaxySDShapeGauss.
 */
NcGalaxySDShapeGauss *
nc_galaxy_sd_shape_gauss_new ()
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

/**
 * nc_galaxy_sd_shape_gauss_set_sigma:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 * @sigma: the standard deviation
 *
 * Sets the standard deviation $\sigma$.
 */
void
nc_galaxy_sd_shape_gauss_set_sigma (NcGalaxySDShapeGauss *gsdsgauss, gdouble sigma)
{
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;

  self->sigma = sigma;
}

/**
 * nc_galaxy_sd_shape_gauss_get_sigma:
 * @gsdsgauss: a #NcGalaxySDShapeGauss
 *
 * Gets the standard deviation $\sigma$.
 *
 * Returns: the standard deviation $\sigma$.
 */
gdouble
nc_galaxy_sd_shape_gauss_get_sigma (NcGalaxySDShapeGauss *gsdsgauss)
{
  NcGalaxySDShapeGaussPrivate * const self = gsdsgauss->priv;

  return self->sigma;
}


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
#include "lss/nc_halo_position.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */


typedef struct _NcGalaxySDShapeGaussPrivate
{
  NcWLSurfaceMassDensityOptzs optzs;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_hp;
} NcGalaxySDShapeGaussPrivate;

struct _NcGalaxySDShapeGauss
{
  NcGalaxySDShape parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShapeGauss, nc_galaxy_sd_shape_gauss, NC_TYPE_GALAXY_SD_SHAPE);

static void
nc_galaxy_sd_shape_gauss_init (NcGalaxySDShapeGauss *gsdsgauss)
{
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  self->ctrl_cosmo = NULL;
  self->ctrl_hp    = NULL;
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsgauss));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE_GAUSS (gsdsgauss));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_shape_gauss_dispose (GObject *object)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (object);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  if (self->ctrl_cosmo)
    ncm_model_ctrl_clear (&self->ctrl_cosmo);

  if (self->ctrl_hp)
    ncm_model_ctrl_clear (&self->ctrl_hp);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_shape_gauss_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_shape_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, NcmRNG *rng, NcGalaxyWLObsCoord coord, NcmVector *data_p, NcmVector *data_z, NcmVector *data_s);
static gdouble _nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, gdouble z, NcmVector *data);
static void _nc_galaxy_sd_shape_gauss_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, NcmVector *data);
static gdouble _nc_galaxy_sd_shape_gauss_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, gdouble z, NcmVector *data);
static gboolean _nc_galaxy_sd_shape_gauss_prepare (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hp, NcGalaxyWLObsCoord coord, gboolean force, NcmObjArray *data, NcmObjArray *data_prep);
static GStrv _nc_galaxy_sd_shape_gauss_get_header (NcGalaxySDShape *gsds);
static guint _nc_galaxy_sd_shape_gauss_get_vec_size (NcGalaxySDShape *gsds);
static gboolean _nc_galaxy_sd_shape_gauss_set_models (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hp);

static void
nc_galaxy_sd_shape_gauss_class_init (NcGalaxySDShapeGaussClass *klass)
{
  NcGalaxySDShapeClass *sd_position_class = NC_GALAXY_SD_SHAPE_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_gauss_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_sd_shape_gauss_dispose;
  object_class->finalize    = &_nc_galaxy_sd_shape_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian galaxy shape distribution", "Gaussian shape");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_SHAPE_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDShapeGauss:e-rms:
   *
   * The intrinsic shape dispersion.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_E_RMS, "\\e-rms", "e-rms", 0.0, 1.0, 1.0e-1, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_E_RMS, NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxySDShapeGauss:e-sigma:
   *
   * The shape measurement error.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_SHAPE_GAUSS_E_SIGMA, "\\e-sigma", "e-sigma", 0.0, 2.0, 1.0e-1, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_E_SIGMA, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  sd_position_class->gen              = &_nc_galaxy_sd_shape_gauss_gen;
  sd_position_class->integ            = &_nc_galaxy_sd_shape_gauss_integ;
  sd_position_class->integ_optzs_prep = &_nc_galaxy_sd_shape_gauss_integ_optzs_prep;
  sd_position_class->integ_optzs      = &_nc_galaxy_sd_shape_gauss_integ_optzs;
  sd_position_class->prepare          = &_nc_galaxy_sd_shape_gauss_prepare;
  sd_position_class->get_header       = &_nc_galaxy_sd_shape_gauss_get_header;
  sd_position_class->get_vec_size     = &_nc_galaxy_sd_shape_gauss_get_vec_size;
  sd_position_class->set_models       = &_nc_galaxy_sd_shape_gauss_set_models;
}

#define VECTOR  (NCM_MODEL (gsds))
#define E_RMS   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_SHAPE_GAUSS_E_RMS))
#define E_SIGMA (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_SHAPE_GAUSS_E_SIGMA))

static void
_nc_galaxy_sd_shape_gauss_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, NcmRNG *rng, NcGalaxyWLObsCoord coord, NcmVector *data_p, NcmVector *data_z, NcmVector *data_s)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble ra_cl                            = ncm_model_param_get_by_name (NCM_MODEL (hp), "ra");
  gdouble dec_cl                           = ncm_model_param_get_by_name (NCM_MODEL (hp), "dec");
  gdouble z_cl                             = ncm_model_param_get_by_name (NCM_MODEL (hp), "z");
  gdouble ra                               = ncm_vector_get (data_p, 0);
  gdouble dec                              = ncm_vector_get (data_p, 1);
  gdouble z                                = ncm_vector_get (data_z, 0);
  gdouble theta                            = 0.0;
  gdouble phi                              = 0.0;
  gdouble et_s                             = ncm_rng_gaussian_gen (rng, 0.0, E_RMS);
  gdouble ex_s                             = ncm_rng_gaussian_gen (rng, 0.0, E_RMS);
  gdouble gt                               = 0.0;
  complex double e_s                       = et_s + I * ex_s;
  complex double e_o                       = e_s;
  gdouble e1, e2;
  gdouble r;

  nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);

  r = nc_halo_position_projected_radius (hp, cosmo, theta);

  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cl, z_cl);

    if (gt > 1.0)
      e_o = (1.0 + gt * conj (e_s)) / (conj (e_s) + gt);
    else
      e_o = (e_s + gt) / (1.0 + gt * e_s);
  }

  e1 = -creal (e_o) * cos (2 * phi) + cimag (e_o) * sin (2 * phi) + ncm_rng_gaussian_gen (rng, 0.0, E_SIGMA);
  e2 = -creal (e_o) * sin (2 * phi) - cimag (e_o) * cos (2 * phi) + ncm_rng_gaussian_gen (rng, 0.0, E_SIGMA);

  if (coord == NC_GALAXY_WL_OBS_COORD_CELESTIAL)
    e2 = -e2;

  ncm_vector_set (data_s, 0, ra);
  ncm_vector_set (data_s, 1, dec);
  ncm_vector_set (data_s, 2, e1);
  ncm_vector_set (data_s, 3, e2);
  ncm_vector_set (data_s, 4, E_RMS);
  ncm_vector_set (data_s, 5, E_SIGMA);
}

static gdouble
_nc_galaxy_sd_shape_gauss_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, gdouble z, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble z_cl                             = ncm_model_param_get_by_name (NCM_MODEL (hp), "z");
  gdouble r                                = ncm_vector_get (data, 0);
  gdouble phi                              = ncm_vector_get (data, 1);
  gdouble e1                               = ncm_vector_get (data, 2);
  gdouble e2                               = ncm_vector_get (data, 3);
  gdouble e_rms                            = ncm_vector_get (data, 4);
  gdouble e_sigma                          = ncm_vector_get (data, 5);
  gdouble gt                               = 0.0;
  complex double e_o                       = e1 + I * e2;
  complex double e_s                       = e_o;
  complex double g                         = 0.0;

  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, r, z, z_cl, z_cl);
    g  = -creal (gt) * cos (2 * phi) - I * creal (gt) * sin (2 * phi);

    if (gt > 1.0)
      e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
    else
      e_s = (e_o - g) / (1.0 - conj (g) * e_o);
  }

  return exp (-0.5 * (cabs (e_s) / (gsl_pow_2 (e_rms) + gsl_pow_2 (e_sigma)))) / sqrt (2.0 * M_PI) / (gsl_pow_2 (e_rms) + gsl_pow_2 (e_sigma));
}

static void
_nc_galaxy_sd_shape_gauss_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble z_cl                             = ncm_model_param_get_by_name (NCM_MODEL (hp), "z");
  gdouble r                                = ncm_vector_get (data, 0);

  nc_wl_surface_mass_density_reduced_shear_optzs_prep (smd, dp, cosmo, r, z_cl, z_cl, &self->optzs);
}

static gdouble
_nc_galaxy_sd_shape_gauss_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, gdouble z, NcmVector *data)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);
  gdouble z_cl                             = ncm_model_param_get_by_name (NCM_MODEL (hp), "z");
  gdouble phi                              = ncm_vector_get (data, 1);
  gdouble e1                               = ncm_vector_get (data, 2);
  gdouble e2                               = ncm_vector_get (data, 3);
  gdouble e_rms                            = ncm_vector_get (data, 4);
  gdouble e_sigma                          = ncm_vector_get (data, 5);
  gdouble gt                               = 0.0;
  complex double e_o                       = e1 + I * e2;
  complex double e_s                       = e_s;
  complex double g                         = 0.0;

  if (z > z_cl)
  {
    gt = nc_wl_surface_mass_density_reduced_shear_optzs (smd, dp, cosmo, z, z_cl, &self->optzs);
    g  = -creal (gt) * cos (2 * phi) - I * creal (gt) * sin (2 * phi);

    if (gt > 1.0)
      e_s = (1.0 - g * conj (e_o)) / (conj (e_o) - conj (g));
    else
      e_s = (e_o - g) / (1.0 - conj (g) * e_o);
  }

  return exp (-0.5 * (cabs (e_s) / (gsl_pow_2 (e_rms) + gsl_pow_2 (e_sigma)))) / sqrt (2.0 * M_PI) / (gsl_pow_2 (e_rms) + gsl_pow_2 (e_sigma));
}

static gboolean
_nc_galaxy_sd_shape_gauss_prepare (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hp, NcGalaxyWLObsCoord coord, gboolean force, NcmObjArray *data, NcmObjArray *data_prep)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  if (ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo)) || ncm_model_ctrl_update (self->ctrl_hp, NCM_MODEL (hp)) || force)
  {
    guint i;

    for (i = 0; i < ncm_obj_array_len (data); i++)
    {
      NcmVector *data_i      = NCM_VECTOR (ncm_obj_array_peek (data, i));
      NcmVector *data_prep_i = NCM_VECTOR (ncm_obj_array_peek (data_prep, i));
      gdouble ra_cl          = ncm_model_param_get_by_name (NCM_MODEL (hp), "ra");
      gdouble dec_cl         = ncm_model_param_get_by_name (NCM_MODEL (hp), "dec");
      gdouble z_cl           = ncm_model_param_get_by_name (NCM_MODEL (hp), "z");
      gdouble ra             = ncm_vector_get (data_i, 0);
      gdouble dec            = ncm_vector_get (data_i, 1);
      gdouble e1             = ncm_vector_get (data_i, 2);
      gdouble e2             = ncm_vector_get (data_i, 3);
      gdouble e_rms          = ncm_vector_get (data_i, 4);
      gdouble e_sigma        = ncm_vector_get (data_i, 5);
      gdouble theta          = 0.0;
      gdouble phi            = 0.0;
      gdouble r;

      nc_halo_position_polar_angles (hp, ra_cl, dec_cl, &theta, &phi);

      r = nc_halo_position_projected_radius (hp, cosmo, theta);

      if (coord == NC_GALAXY_WL_OBS_COORD_CELESTIAL)
        phi = M_PI - phi;

      ncm_vector_set (data_prep_i, 0, r);
      ncm_vector_set (data_prep_i, 1, phi);
      ncm_vector_set (data_prep_i, 2, e1);
      ncm_vector_set (data_prep_i, 3, e2);
      ncm_vector_set (data_prep_i, 4, e_rms);
      ncm_vector_set (data_prep_i, 5, e_sigma);
    }

    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static GStrv
_nc_galaxy_sd_shape_gauss_get_header (NcGalaxySDShape *gsds)
{
  GStrv header = g_strsplit ("ra dec e1 e2 e_rms e_sigma", " ", -1);

  return header;
}

static guint
_nc_galaxy_sd_shape_gauss_get_vec_size (NcGalaxySDShape *gsds)
{
  return 6;
}

static gboolean
_nc_galaxy_sd_shape_gauss_set_models (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hp)
{
  NcGalaxySDShapeGauss *gsdsgauss          = NC_GALAXY_SD_SHAPE_GAUSS (gsds);
  NcGalaxySDShapeGaussPrivate * const self = nc_galaxy_sd_shape_gauss_get_instance_private (gsdsgauss);

  self->ctrl_cosmo = ncm_model_ctrl_new (NCM_MODEL (cosmo));
  self->ctrl_hp    = ncm_model_ctrl_new (NCM_MODEL (hp));

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
nc_galaxy_sd_shape_gauss_new ()
{
  NcGalaxySDShapeGauss *gsdsgauss = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_GAUSS, NULL);

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


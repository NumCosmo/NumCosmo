/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_reduced_shear_gauss.c
 *
 *  Mon July 27 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_reduced_shear_gauss.c
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:nc_galaxy_wl_reduced_shear_gauss
 * @title: NcGalaxyWLReducedShearGauss
 * @short_description: Abstract class describing galaxy weak lensing reduced shear Gaussian distribution
 * @stability: Unstable
 *
 *
 * Class defining a galaxy weak lensing reduced shear normally distributed.
 * probability distribution $P_\mathrm{wl}(g)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_reduced_shear_gauss.h"
#include "nc_enum_types.h"

#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLReducedShearGaussPrivate
{
  NcmMatrix *obs;
  NcGalaxyWLReducedShearGaussPos pos;
  gdouble norma;
  gdouble r;
  gdouble twolnN;
  guint len;
  NcWLSurfaceMassDensityOptzs optzs;
};

enum
{
  PROP_0,
  PROP_POS,
  PROP_OBS,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLReducedShearGauss, nc_galaxy_wl_reduced_shear_gauss, NC_TYPE_GALAXY_WL_DIST);

static void
nc_galaxy_wl_reduced_shear_gauss_init (NcGalaxyWLReducedShearGauss *grsg)
{
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv = nc_galaxy_wl_reduced_shear_gauss_get_instance_private (grsg);
  
  self->obs    = NULL;
  self->pos    = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_LEN;
  self->norma  = 0.0;
  self->r      = 0.0;
  self->twolnN = 0.0;
  self->len    = 0;
}

static void
_nc_galaxy_wl_reduced_shear_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLReducedShearGauss *grsg = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS (object);
  
  g_return_if_fail (NC_IS_GALAXY_WL_REDUCED_SHEAR_GAUSS (object));
  
  switch (prop_id)
  {
    case PROP_POS:
      nc_galaxy_wl_reduced_shear_gauss_set_pos (grsg, g_value_get_enum (value));
      break;
    case PROP_OBS:
      nc_galaxy_wl_reduced_shear_gauss_set_obs (grsg, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_reduced_shear_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLReducedShearGauss *grsg = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS (object);
  
  g_return_if_fail (NC_IS_GALAXY_WL_REDUCED_SHEAR_GAUSS (object));
  
  switch (prop_id)
  {
    case PROP_POS:
      g_value_set_enum (value, nc_galaxy_wl_reduced_shear_gauss_get_pos (grsg));
      break;
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_reduced_shear_gauss_peek_obs (grsg));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_reduced_shear_gauss_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_reduced_shear_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_reduced_shear_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_reduced_shear_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_wl_reduced_shear_gauss_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i);
static gdouble _nc_galaxy_wl_reduced_shear_gauss_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
static gdouble _nc_galaxy_wl_reduced_shear_gauss_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
static guint _nc_galaxy_wl_reduced_shear_gauss_len (NcGalaxyWLDist *gwld);

static void
nc_galaxy_wl_reduced_shear_gauss_class_init (NcGalaxyWLReducedShearGaussClass *klass)
{
  NcGalaxyWLDistClass *wl_dist_class = NC_GALAXY_WL_DIST_CLASS (klass);
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_nc_galaxy_wl_reduced_shear_gauss_set_property;
  object_class->get_property = &_nc_galaxy_wl_reduced_shear_gauss_get_property;
  object_class->dispose      = &_nc_galaxy_wl_reduced_shear_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_wl_reduced_shear_gauss_finalize;
  
  /**
   * NcGalaxyWLReducedShearGauss:pos:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_POS,
                                   g_param_spec_enum ("pos",
                                                      NULL,
                                                      "Observable position type",
                                                      NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS, NC_GALAXY_WL_REDUCED_SHEAR_GAUSS_POS_R,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcGalaxyWLReducedShearGauss:obs:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Galaxy observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  wl_dist_class->m2lnP_prep = &_nc_galaxy_wl_reduced_shear_gauss_m2lnP_prep;
  wl_dist_class->m2lnP      = &_nc_galaxy_wl_reduced_shear_gauss_m2lnP;
  wl_dist_class->gen        = &_nc_galaxy_wl_reduced_shear_gauss_gen;
  wl_dist_class->len        = &_nc_galaxy_wl_reduced_shear_gauss_len;
}

static void
_nc_galaxy_wl_reduced_shear_gauss_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i)
{
  NcGalaxyWLReducedShearGauss *grsg               = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS (gwld);
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  const gdouble R                                 = ncm_matrix_get (self->obs, gal_i, 0);
  
  nc_wl_surface_mass_density_reduced_shear_optzs_prep (smd, dp, cosmo, R, z_cluster, z_cluster, &self->optzs);
}

static gdouble
_nc_galaxy_wl_reduced_shear_gauss_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  NcGalaxyWLReducedShearGauss *grsg               = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS (gwld);
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  const gdouble g_mean                            = ncm_matrix_get (self->obs, gal_i, 1);
  const gdouble sigma_g                           = ncm_matrix_get (self->obs, gal_i, 2);
  const gdouble g_th                              = nc_wl_surface_mass_density_reduced_shear_optzs (smd, dp, cosmo, z, z_cluster, &self->optzs); /*nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, R, z, z_cluster, z_cluster); */
  
  return gsl_pow_2 ((g_th - g_mean) / sigma_g); /* + self->twolnN; */
}

static gdouble
_nc_galaxy_wl_reduced_shear_gauss_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLReducedShearGauss *grsg               = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS (gwld);
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  const gdouble sigma_g                           = ncm_matrix_get (self->obs, 0, 2);
  
  return ncm_rng_gaussian_gen (rng, g_true, sigma_g);
}

static guint
_nc_galaxy_wl_reduced_shear_gauss_len (NcGalaxyWLDist *gwld)
{
  NcGalaxyWLReducedShearGauss *grsg               = NC_GALAXY_WL_REDUCED_SHEAR_GAUSS (gwld);
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  
  return self->len;
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_new:
 * @pos: a #NcGalaxyWLReducedShearGaussPos
 *
 * Creates a new #NcGalaxyWLReducedShearGauss using
 * @pos as the position type.
 *
 * Returns: (transfer full): a new NcGalaxyWLReducedShearGauss.
 */
NcGalaxyWLReducedShearGauss *
nc_galaxy_wl_reduced_shear_gauss_new (NcGalaxyWLReducedShearGaussPos pos)
{
  NcGalaxyWLReducedShearGauss *grsg = g_object_new (NC_TYPE_GALAXY_WL_REDUCED_SHEAR_GAUSS,
                                                    "pos", pos,
                                                    NULL);
  
  return grsg;
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_ref:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 *
 * Increase the reference of @grsg by one.
 *
 * Returns: (transfer full): @grsg.
 */
NcGalaxyWLReducedShearGauss *
nc_galaxy_wl_reduced_shear_gauss_ref (NcGalaxyWLReducedShearGauss *grsg)
{
  return g_object_ref (grsg);
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_free:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 *
 * Decrease the reference count of @grsg by one.
 *
 */
void
nc_galaxy_wl_reduced_shear_gauss_free (NcGalaxyWLReducedShearGauss *grsg)
{
  g_object_unref (grsg);
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_clear:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 *
 * Decrease the reference count of @grsg by one, and sets the pointer *@grsg to
 * NULL.
 *
 */
void
nc_galaxy_wl_reduced_shear_gauss_clear (NcGalaxyWLReducedShearGauss **grsg)
{
  g_clear_object (grsg);
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_set_pos:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 * @pos: a #NcGalaxyWLReducedShearGaussPos
 *
 * Sets the position observable type.
 */
void
nc_galaxy_wl_reduced_shear_gauss_set_pos (NcGalaxyWLReducedShearGauss *grsg, NcGalaxyWLReducedShearGaussPos pos)
{
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  
  self->pos = pos;
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_get_pos:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 *
 * Gets the position observable type.
 *
 * Returns: the position observable type.
 */
NcGalaxyWLReducedShearGaussPos
nc_galaxy_wl_reduced_shear_gauss_get_pos (NcGalaxyWLReducedShearGauss *grsg)
{
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  
  return self->pos;
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_set_obs:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 * @obs: a #NcmMatrix
 *
 * Sets the observables matrix @obs.
 */
void
nc_galaxy_wl_reduced_shear_gauss_set_obs (NcGalaxyWLReducedShearGauss *grsg, NcmMatrix *obs)
{
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  
  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);
  
  ncm_matrix_clear (&self->obs);
  
  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_reduced_shear_gauss_peek_obs:
 * @grsg: a #NcGalaxyWLReducedShearGauss
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_reduced_shear_gauss_peek_obs (NcGalaxyWLReducedShearGauss *grsg)
{
  NcGalaxyWLReducedShearGaussPrivate * const self = grsg->priv;
  
  return self->obs;
}


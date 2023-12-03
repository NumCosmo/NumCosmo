/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_ellipticity_gauss.c
 *
 *  Mon July 27 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_ellipticity_gauss.c
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_galaxy_wl_ellipticity_gauss
 * @title: NcGalaxyWLEllipticityGauss
 * @short_description: Abstract class describing galaxy weak lensing ellipticity Gaussian distribution
 * @stability: Unstable
 *
 *
 * Class defining a galaxy weak lensing ellipticity normally distributed.
 * probability distribution $P_\mathrm{wl}(g)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_ellipticity_gauss.h"
#include "nc_enum_types.h"

#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLEllipticityGaussPrivate
{
  NcmMatrix *obs;
  NcGalaxyWLEllipticityGaussPos pos;
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

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLEllipticityGauss, nc_galaxy_wl_ellipticity_gauss, NC_TYPE_GALAXY_WL_DIST)

static void
nc_galaxy_wl_ellipticity_gauss_init (NcGalaxyWLEllipticityGauss *gegauss)
{
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv = nc_galaxy_wl_ellipticity_gauss_get_instance_private (gegauss);

  self->obs    = NULL;
  self->pos    = NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_LEN;
  self->norma  = 0.0;
  self->r      = 0.0;
  self->twolnN = 0.0;
  self->len    = 0;
}

static void
_nc_galaxy_wl_ellipticity_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityGauss *gegauss = NC_GALAXY_WL_ELLIPTICITY_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_GAUSS (object));

  switch (prop_id)
  {
    case PROP_POS:
      nc_galaxy_wl_ellipticity_gauss_set_pos (gegauss, g_value_get_enum (value));
      break;
    case PROP_OBS:
      nc_galaxy_wl_ellipticity_gauss_set_obs (gegauss, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLEllipticityGauss *gegauss = NC_GALAXY_WL_ELLIPTICITY_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_WL_ELLIPTICITY_GAUSS (object));

  switch (prop_id)
  {
    case PROP_POS:
      g_value_set_enum (value, nc_galaxy_wl_ellipticity_gauss_get_pos (gegauss));
      break;
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_ellipticity_gauss_peek_obs (gegauss));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_ellipticity_gauss_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_ellipticity_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_ellipticity_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_wl_ellipticity_gauss_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i);
static gdouble _nc_galaxy_wl_ellipticity_gauss_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
static gdouble _nc_galaxy_wl_ellipticity_gauss_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
static guint _nc_galaxy_wl_ellipticity_gauss_len (NcGalaxyWLDist *gwld);

static void
nc_galaxy_wl_ellipticity_gauss_class_init (NcGalaxyWLEllipticityGaussClass *klass)
{
  NcGalaxyWLDistClass *wl_dist_class = NC_GALAXY_WL_DIST_CLASS (klass);
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_ellipticity_gauss_set_property;
  object_class->get_property = &_nc_galaxy_wl_ellipticity_gauss_get_property;
  object_class->dispose      = &_nc_galaxy_wl_ellipticity_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_wl_ellipticity_gauss_finalize;

  /**
   * NcGalaxyWLEllipticityGauss:pos:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_POS,
                                   g_param_spec_enum ("pos",
                                                      NULL,
                                                      "Observable position type",
                                                      NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS_POS, NC_GALAXY_WL_ELLIPTICITY_GAUSS_POS_R,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLEllipticityGauss:obs:
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

  wl_dist_class->m2lnP_prep = &_nc_galaxy_wl_ellipticity_gauss_m2lnP_prep;
  wl_dist_class->m2lnP      = &_nc_galaxy_wl_ellipticity_gauss_m2lnP;
  wl_dist_class->gen        = &_nc_galaxy_wl_ellipticity_gauss_gen;
  wl_dist_class->len        = &_nc_galaxy_wl_ellipticity_gauss_len;
}

static void
_nc_galaxy_wl_ellipticity_gauss_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i)
{
  NcGalaxyWLEllipticityGauss *gegauss            = NC_GALAXY_WL_ELLIPTICITY_GAUSS (gwld);
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;
  const gdouble R                                = ncm_matrix_get (self->obs, gal_i, 0);

  nc_wl_surface_mass_density_reduced_shear_optzs_prep (smd, dp, cosmo, R, z_cluster, z_cluster, &self->optzs);
}

static gdouble
_nc_galaxy_wl_ellipticity_gauss_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  NcGalaxyWLEllipticityGauss *gegauss            = NC_GALAXY_WL_ELLIPTICITY_GAUSS (gwld);
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;
  const gdouble g_mean                           = ncm_matrix_get (self->obs, gal_i, 1);
  const gdouble sigma_g                          = ncm_matrix_get (self->obs, gal_i, 2);
  const gdouble g_th                             = nc_wl_surface_mass_density_reduced_shear_optzs (smd, dp, cosmo, z, z_cluster, &self->optzs); /*nc_wl_surface_mass_density_ellipticity (smd, dp, cosmo, R, z, z_cluster, z_cluster); */

  return gsl_pow_2 ((g_th - g_mean) / sigma_g); /* + self->twolnN; */
}

static gdouble
_nc_galaxy_wl_ellipticity_gauss_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLEllipticityGauss *gegauss            = NC_GALAXY_WL_ELLIPTICITY_GAUSS (gwld);
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;
  const gdouble sigma_g                          = ncm_matrix_get (self->obs, 0, 2);

  return ncm_rng_gaussian_gen (rng, g_true, sigma_g);
}

static guint
_nc_galaxy_wl_ellipticity_gauss_len (NcGalaxyWLDist *gwld)
{
  NcGalaxyWLEllipticityGauss *gegauss            = NC_GALAXY_WL_ELLIPTICITY_GAUSS (gwld);
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;

  return self->len;
}

/**
 * nc_galaxy_wl_ellipticity_gauss_new:
 * @pos: a #NcGalaxyWLEllipticityGaussPos
 *
 * Creates a new #NcGalaxyWLEllipticityGauss using
 * @pos as the position type.
 *
 * Returns: (transfer full): a new NcGalaxyWLEllipticityGauss.
 */
NcGalaxyWLEllipticityGauss *
nc_galaxy_wl_ellipticity_gauss_new (NcGalaxyWLEllipticityGaussPos pos)
{
  NcGalaxyWLEllipticityGauss *gegauss = g_object_new (NC_TYPE_GALAXY_WL_ELLIPTICITY_GAUSS,
                                                      "pos", pos,
                                                      NULL);

  return gegauss;
}

/**
 * nc_galaxy_wl_ellipticity_gauss_ref:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 *
 * Increase the reference of @gegauss by one.
 *
 * Returns: (transfer full): @gegauss.
 */
NcGalaxyWLEllipticityGauss *
nc_galaxy_wl_ellipticity_gauss_ref (NcGalaxyWLEllipticityGauss *gegauss)
{
  return g_object_ref (gegauss);
}

/**
 * nc_galaxy_wl_ellipticity_gauss_free:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 *
 * Decrease the reference count of @gegauss by one.
 *
 */
void
nc_galaxy_wl_ellipticity_gauss_free (NcGalaxyWLEllipticityGauss *gegauss)
{
  g_object_unref (gegauss);
}

/**
 * nc_galaxy_wl_ellipticity_gauss_clear:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 *
 * Decrease the reference count of @gegauss by one, and sets the pointer *@gegauss to
 * NULL.
 *
 */
void
nc_galaxy_wl_ellipticity_gauss_clear (NcGalaxyWLEllipticityGauss **gegauss)
{
  g_clear_object (gegauss);
}

/**
 * nc_galaxy_wl_ellipticity_gauss_set_pos:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 * @pos: a #NcGalaxyWLEllipticityGaussPos
 *
 * Sets the position observable type.
 */
void
nc_galaxy_wl_ellipticity_gauss_set_pos (NcGalaxyWLEllipticityGauss *gegauss, NcGalaxyWLEllipticityGaussPos pos)
{
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;

  self->pos = pos;
}

/**
 * nc_galaxy_wl_ellipticity_gauss_get_pos:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 *
 * Gets the position observable type.
 *
 * Returns: the position observable type.
 */
NcGalaxyWLEllipticityGaussPos
nc_galaxy_wl_ellipticity_gauss_get_pos (NcGalaxyWLEllipticityGauss *gegauss)
{
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;

  return self->pos;
}

/**
 * nc_galaxy_wl_ellipticity_gauss_set_obs:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 * @obs: a #NcmMatrix
 *
 * Sets the observable matrix @obs.
 */
void
nc_galaxy_wl_ellipticity_gauss_set_obs (NcGalaxyWLEllipticityGauss *gegauss, NcmMatrix *obs)
{
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;

  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);

  ncm_matrix_clear (&self->obs);

  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_ellipticity_gauss_peek_obs:
 * @gegauss: a #NcGalaxyWLEllipticityGauss
 *
 * Gets the observable matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_ellipticity_gauss_peek_obs (NcGalaxyWLEllipticityGauss *gegauss)
{
  NcGalaxyWLEllipticityGaussPrivate * const self = gegauss->priv;

  return self->obs;
}


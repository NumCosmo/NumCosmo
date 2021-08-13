/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_proj.c
 *
 *  Sun Aug 02 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_proj.c
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
 * SECTION:nc_galaxy_wl_proj
 * @title: NcGalaxyWLProj
 * @short_description: Abstract class describing galaxy weak lensing reduced shear Gaussian distribution
 * @stability: Unstable
 *
 * Class defining a galaxy weak lensing projected distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_proj.h"
#include "nc_enum_types.h"

#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLProjPrivate
{
  NcmMatrix *obs;
  NcGalaxyWLProjPos pos;
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

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLProj, nc_galaxy_wl_proj, NC_TYPE_GALAXY_WL_DIST);

static void
nc_galaxy_wl_proj_init (NcGalaxyWLProj *gwlp)
{
  NcGalaxyWLProjPrivate * const self = gwlp->priv = nc_galaxy_wl_proj_get_instance_private (gwlp);
  
  self->obs    = NULL;
  self->pos    = NC_GALAXY_WL_PROJ_POS_LEN;
  self->norma  = 0.0;
  self->r      = 0.0;
  self->twolnN = 0.0;
  self->len    = 0;
}

static void
_nc_galaxy_wl_proj_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLProj *gwlp = NC_GALAXY_WL_PROJ (object);
  
  g_return_if_fail (NC_IS_GALAXY_WL_PROJ (object));
  
  switch (prop_id)
  {
    case PROP_POS:
      nc_galaxy_wl_proj_set_pos (gwlp, g_value_get_enum (value));
      break;
    case PROP_OBS:
      nc_galaxy_wl_proj_set_obs (gwlp, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_proj_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLProj *gwlp = NC_GALAXY_WL_PROJ (object);
  
  g_return_if_fail (NC_IS_GALAXY_WL_PROJ (object));
  
  switch (prop_id)
  {
    case PROP_POS:
      g_value_set_enum (value, nc_galaxy_wl_proj_get_pos (gwlp));
      break;
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_wl_proj_peek_obs (gwlp));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_proj_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_proj_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_proj_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_proj_parent_class)->finalize (object);
}

static void _nc_galaxy_wl_proj_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i);
static gdouble _nc_galaxy_wl_proj_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
static gdouble _nc_galaxy_wl_proj_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
static guint _nc_galaxy_wl_proj_len (NcGalaxyWLDist *gwld);

static void
nc_galaxy_wl_proj_class_init (NcGalaxyWLProjClass *klass)
{
  NcGalaxyWLDistClass *wl_dist_class = NC_GALAXY_WL_DIST_CLASS (klass);
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_nc_galaxy_wl_proj_set_property;
  object_class->get_property = &_nc_galaxy_wl_proj_get_property;
  object_class->dispose      = &_nc_galaxy_wl_proj_dispose;
  object_class->finalize     = &_nc_galaxy_wl_proj_finalize;
  
  /**
   * NcGalaxyWLProj:pos:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_POS,
                                   g_param_spec_enum ("pos",
                                                      NULL,
                                                      "Observable position type",
                                                      NC_TYPE_GALAXY_WL_PROJ_POS, NC_GALAXY_WL_PROJ_POS_R,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcGalaxyWLProj:obs:
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
  
  wl_dist_class->m2lnP_prep = &_nc_galaxy_wl_proj_m2lnP_prep;
  wl_dist_class->m2lnP      = &_nc_galaxy_wl_proj_m2lnP;
  wl_dist_class->gen        = &_nc_galaxy_wl_proj_gen;
  wl_dist_class->len        = &_nc_galaxy_wl_proj_len;
}

static void
_nc_galaxy_wl_proj_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i)
{
  NcGalaxyWLProj *gwlp               = NC_GALAXY_WL_PROJ (gwld);
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  const gdouble R                    = ncm_matrix_get (self->obs, gal_i, 0);
  
  nc_wl_surface_mass_density_reduced_shear_optzs_prep (smd, dp, cosmo, R, z_cluster, z_cluster, &self->optzs);
}

static gdouble
_nc_galaxy_wl_proj_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  NcGalaxyWLProj *gwlp               = NC_GALAXY_WL_PROJ (gwld);
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  const gdouble g_mean               = ncm_matrix_get (self->obs, gal_i, 1);
  const gdouble sigma_g              = ncm_matrix_get (self->obs, gal_i, 2);
  const gdouble g_th                 = nc_wl_surface_mass_density_reduced_shear_optzs (smd, dp, cosmo, z, z_cluster, &self->optzs); /*nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, R, z, z_cluster, z_cluster); */
  
  return gsl_pow_2 ((g_th - g_mean) / sigma_g); /* + self->twolnN; */
}

static gdouble
_nc_galaxy_wl_proj_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLProj *gwlp               = NC_GALAXY_WL_PROJ (gwld);
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  const gdouble sigma_g              = ncm_matrix_get (self->obs, 0, 2);
  
  return ncm_rng_gaussian_gen (rng, g_true, sigma_g);
}

static guint
_nc_galaxy_wl_proj_len (NcGalaxyWLDist *gwld)
{
  NcGalaxyWLProj *gwlp               = NC_GALAXY_WL_PROJ (gwld);
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  
  return self->len;
}

/**
 * nc_galaxy_wl_proj_new:
 * @pos: a #NcGalaxyWLProjPos
 *
 * Creates a new #NcGalaxyWLProj using
 * @pos as the position type.
 *
 * Returns: (transfer full): a new NcGalaxyWLProj.
 */
NcGalaxyWLProj *
nc_galaxy_wl_proj_new (NcGalaxyWLProjPos pos)
{
  NcGalaxyWLProj *gwlp = g_object_new (NC_TYPE_GALAXY_WL_PROJ,
                                       "pos", pos,
                                       NULL);
  
  return gwlp;
}

/**
 * nc_galaxy_wl_proj_ref:
 * @gwlp: a #NcGalaxyWLProj
 *
 * Increase the reference of @gwlp by one.
 *
 * Returns: (transfer full): @gwlp.
 */
NcGalaxyWLProj *
nc_galaxy_wl_proj_ref (NcGalaxyWLProj *gwlp)
{
  return g_object_ref (gwlp);
}

/**
 * nc_galaxy_wl_proj_free:
 * @gwlp: a #NcGalaxyWLProj
 *
 * Decrease the reference count of @gwlp by one.
 *
 */
void
nc_galaxy_wl_proj_free (NcGalaxyWLProj *gwlp)
{
  g_object_unref (gwlp);
}

/**
 * nc_galaxy_wl_proj_clear:
 * @gwlp: a #NcGalaxyWLProj
 *
 * Decrease the reference count of @gwlp by one, and sets the pointer *@gwlp to
 * NULL.
 *
 */
void
nc_galaxy_wl_proj_clear (NcGalaxyWLProj **gwlp)
{
  g_clear_object (gwlp);
}

/**
 * nc_galaxy_wl_proj_set_pos:
 * @gwlp: a #NcGalaxyWLProj
 * @pos: a #NcGalaxyWLProjPos
 *
 * Sets the position observable type.
 */
void
nc_galaxy_wl_proj_set_pos (NcGalaxyWLProj *gwlp, NcGalaxyWLProjPos pos)
{
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  
  self->pos = pos;
}

/**
 * nc_galaxy_wl_proj_get_pos:
 * @gwlp: a #NcGalaxyWLProj
 *
 * Gets the position observable type.
 *
 * Returns: the position observable type.
 */
NcGalaxyWLProjPos
nc_galaxy_wl_proj_get_pos (NcGalaxyWLProj *gwlp)
{
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  
  return self->pos;
}

/**
 * nc_galaxy_wl_proj_set_obs:
 * @gwlp: a #NcGalaxyWLProj
 * @obs: a #NcmMatrix
 *
 * Sets the observables matrix @obs.
 */
void
nc_galaxy_wl_proj_set_obs (NcGalaxyWLProj *gwlp, NcmMatrix *obs)
{
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  
  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 3);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);
  
  ncm_matrix_clear (&self->obs);
  
  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
}

/**
 * nc_galaxy_wl_proj_peek_obs:
 * @gwlp: a #NcGalaxyWLProj
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcmMatrix *
nc_galaxy_wl_proj_peek_obs (NcGalaxyWLProj *gwlp)
{
  NcGalaxyWLProjPrivate * const self = gwlp->priv;
  
  return self->obs;
}


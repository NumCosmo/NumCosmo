/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_proj.c
 *
 *  Thu Jul 30 11:12:53 2020
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
 *
 * Class defining a galaxy weak lensing reduced shear normally distributed.
 * probability distribution $P_\mathrm{wl}(g)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_proj.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLProjPrivate
{
  gdouble g_mean;
  gdouble sigma_g;
  gdouble R;
  gdouble norma;
  gdouble r;
  gdouble twolnN;
};

enum
{
  PROP_0,
  PROP_G_MEAN,
  PROP_SIGMA_G,
  PROP_R,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLProj, nc_galaxy_wl_proj, NC_TYPE_GALAXY_WL_DIST);

static void
nc_galaxy_wl_proj_init (NcGalaxyWLProj *gp)
{
  NcGalaxyWLProjPrivate * const self = gp->priv = nc_galaxy_wl_proj_get_instance_private (gp);
  
  self->g_mean  = 0.0;
  self->sigma_g = 0.0;
  self->R       = 0.0;
  self->norma   = 0.0;
  self->r       = 0.0;
  self->twolnN  = 0.0;
}

static void
_nc_galaxy_wl_proj_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLProj *gp = NC_GALAXY_WL_PROJ (object);

  g_return_if_fail (NC_IS_GALAXY_WL_PROJ (object));

  switch (prop_id)
  {
    case PROP_G_MEAN:
      nc_galaxy_wl_proj_set_g_obs (gp, g_value_get_double (value));
      break;
    case PROP_SIGMA_G:
      nc_galaxy_wl_proj_set_sigma_g (gp, g_value_get_double (value));
      break;
    case PROP_R:
      nc_galaxy_wl_proj_set_R (gp, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_proj_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLProj *gp = NC_GALAXY_WL_PROJ (object);

  g_return_if_fail (NC_IS_GALAXY_WL_PROJ (object));

  switch (prop_id)
  {
    case PROP_G_MEAN:
      g_value_set_double (value, nc_galaxy_wl_proj_get_g_obs (gp));
      break;
    case PROP_SIGMA_G:
      g_value_set_double (value, nc_galaxy_wl_proj_get_sigma_g (gp));
      break;
    case PROP_R:
      g_value_set_double (value, nc_galaxy_wl_proj_get_R (gp));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_proj_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_proj_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_wl_proj_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z);
static gdouble _nc_galaxy_wl_proj_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);

static void
nc_galaxy_wl_proj_class_init (NcGalaxyWLProjClass *klass)
{
  NcGalaxyWLDistClass *wl_dist_class = NC_GALAXY_WL_DIST_CLASS (klass);
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_nc_galaxy_wl_proj_set_property;
  object_class->get_property = &_nc_galaxy_wl_proj_get_property;
  object_class->finalize     = &_nc_galaxy_wl_proj_finalize;

  /**
   * NcGalaxyWLProj:g-mean:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_G_MEAN,
                                   g_param_spec_double ("g-mean",
                                                        NULL,
                                                        "Observed value for g",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcGalaxyWLProj:sigma-g:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_G,
                                   g_param_spec_double ("sigma-g",
                                                        NULL,
                                                        "Observed value for g standard deviation",
                                                        1.0e-200, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcGalaxyWLProj:R:
   *
   * FIXME [Mpc]
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R,
                                   g_param_spec_double ("R",
                                                        NULL,
                                                        "Distance with respect to cluster center [Mpc]",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  wl_dist_class->m2lnP = &_nc_galaxy_wl_proj_m2lnP;
  wl_dist_class->gen   = &_nc_galaxy_wl_proj_gen;
}

static gdouble
_nc_galaxy_wl_proj_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z)
{  
  NcGalaxyWLProj *gp = NC_GALAXY_WL_PROJ (gwld);
  NcGalaxyWLProjPrivate * const self = gp->priv;
  const gdouble g_th = nc_wl_surface_mass_density_reduced_shear (smd, dp, cosmo, self->R, z, z_cluster, z_cluster);

  return gsl_pow_2 (g_th - self->g_mean) * self->r;// + self->twolnN;
}

static gdouble
_nc_galaxy_wl_proj_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  NcGalaxyWLProj *gp = NC_GALAXY_WL_PROJ (gwld);
  NcGalaxyWLProjPrivate * const self = gp->priv;
  return ncm_rng_gaussian_gen (rng, g_true, self->sigma_g);
}

/**
 * nc_galaxy_wl_proj_new:
 * @g_obs: the observed value of $g$
 * @sigma_g: the observed value of $g$ standard deviation
 *
 * Returns: (transfer full): a new NcGalaxyWLProj.
 */
NcGalaxyWLProj *
nc_galaxy_wl_proj_new (const gdouble g_obs, const gdouble sigma_g)
{
  NcGalaxyWLProj *gp = g_object_new (NC_TYPE_GALAXY_WL_PROJ,
                                                    "g-mean", g_obs,
                                                    "sigma-g", sigma_g,
                                                    NULL);

  return gp;
}

/**
 * nc_galaxy_wl_proj_ref:
 * @gp: a #NcGalaxyWLProj
 *
 * Increase the reference of @gp by one.
 *
 * Returns: (transfer full): @gp.
 */
NcGalaxyWLProj *
nc_galaxy_wl_proj_ref (NcGalaxyWLProj *gp)
{
  return g_object_ref (gp);
}

/**
 * nc_galaxy_wl_proj_free:
 * @gp: a #NcGalaxyWLProj
 *
 * Decrease the reference count of @gp by one.
 *
 */
void
nc_galaxy_wl_proj_free (NcGalaxyWLProj *gp)
{
  g_object_unref (gp);
}

/**
 * nc_galaxy_wl_proj_clear:
 * @gp: a #NcGalaxyWLProj
 *
 * Decrease the reference count of @gp by one, and sets the pointer *@gp to
 * NULL.
 *
 */
void
nc_galaxy_wl_proj_clear (NcGalaxyWLProj **gp)
{
  g_clear_object (gp);
}

/**
 * nc_galaxy_wl_proj_set_g_obs:
 * @gp: a #NcGalaxyWLProj
 * @g_obs: the observed value of $g$
 * 
 * Sets the value of $g_\mathrm{obs}$.
 */
void
nc_galaxy_wl_proj_set_g_obs (NcGalaxyWLProj *gp, const gdouble g_obs)
{
  NcGalaxyWLProjPrivate * const self = gp->priv;
  
  self->g_mean = g_obs;
}

/**
 * nc_galaxy_wl_proj_set_sigma_g:
 * @gp: a #NcGalaxyWLProj
 * @sigma_g: the observed value of $g$ standard deviation
 * 
 * Sets the value of $\sigma_g$.
 * 
 */
void
nc_galaxy_wl_proj_set_sigma_g (NcGalaxyWLProj *gp, const gdouble sigma_g)
{
  NcGalaxyWLProjPrivate * const self = gp->priv;

  g_assert_cmpfloat (sigma_g, >, 0.0);
  
  self->sigma_g = sigma_g;
  self->norma   = sqrt (2.0 * M_PI) * sigma_g;
  self->r       = 1.0 / gsl_pow_2 (sigma_g);
  self->twolnN  = 2.0 * log (self->norma);
}

/**
 * nc_galaxy_wl_proj_set_R:
 * @gp: a #NcGalaxyWLProj
 * @R: the radial distance from the cluster center
 *
 * Sets the value of $R\,\mathrm{Mpc}$.
 *
 */
void
nc_galaxy_wl_proj_set_R (NcGalaxyWLProj *gp, const gdouble R)
{
  NcGalaxyWLProjPrivate * const self = gp->priv;

  g_assert_cmpfloat (R, >, 0.0);

  self->R = R;
}

/**
 * nc_galaxy_wl_proj_get_g_obs:
 * @gp: a #NcGalaxyWLProj
 * 
 * Gets the value of $g_\mathrm{obs}$.
 *
 * Returns: $g_\mathrm{obs}$.
 */
gdouble
nc_galaxy_wl_proj_get_g_obs (NcGalaxyWLProj *gp)
{
  NcGalaxyWLProjPrivate * const self = gp->priv;
  return self->g_mean;
}

/**
 * nc_galaxy_wl_proj_get_sigma_g:
 * @gp: a #NcGalaxyWLProj
 * 
 * Gets the value of $\sigma_g$.
 * 
 * Returns: $\sigma_g$.
 */
gdouble
nc_galaxy_wl_proj_get_sigma_g (NcGalaxyWLProj *gp)
{
  NcGalaxyWLProjPrivate * const self = gp->priv;

  return self->sigma_g;
}

/**
 * nc_galaxy_wl_proj_get_R:
 * @gp: a #NcGalaxyWLProj
 *
 * Gets the value of the radial distance from the cluster center $R\,\mathrm{Mpc}$.
 *
 * Returns: $R\,\mathrm{Mpc}$.
 */
gdouble
nc_galaxy_wl_proj_get_R (NcGalaxyWLProj *gp)
{
  NcGalaxyWLProjPrivate * const self = gp->priv;

  return self->R;
}


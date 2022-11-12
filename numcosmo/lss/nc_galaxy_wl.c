/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl.c
 *
 *  Mon July 27 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl.c
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
 * SECTION:nc_galaxy_wl
 * @title: NcGalaxyWL
 * @short_description: Object containing galaxy weak lensing data
 * @stability: Unstable
 *
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl.h"
#include "lss/nc_galaxy_wl_dist.h"
#include "lss/nc_galaxy_redshift.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxyWLPrivate
{
  NcGalaxyWLDist *wl_dist;
  NcGalaxyRedshift *gz_dist;
  guint len;
};

enum
{
  PROP_0,
  PROP_WL_DIST,
  PROP_GZ_DIST,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWL, nc_galaxy_wl, G_TYPE_OBJECT);

static void
nc_galaxy_wl_init (NcGalaxyWL *gwl)
{
  NcGalaxyWLPrivate * const self = gwl->priv = nc_galaxy_wl_get_instance_private (gwl);
  
  self->gz_dist = NULL;
  self->wl_dist = NULL;
  self->len     = 0;
}

static void
_nc_galaxy_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWL *gwl                = NC_GALAXY_WL (object);
  NcGalaxyWLPrivate * const self = gwl->priv;
  
  g_return_if_fail (NC_IS_GALAXY_WL (object));
  
  switch (prop_id)
  {
    case PROP_WL_DIST:
      self->wl_dist = g_value_dup_object (value);
      
      if ((self->gz_dist != NULL) && (self->wl_dist != NULL))
      {
        g_assert_cmpuint (nc_galaxy_wl_dist_len (self->wl_dist), ==, nc_galaxy_redshift_len (self->gz_dist));
        self->len = nc_galaxy_wl_dist_len (self->wl_dist);
      }
      
      break;
    case PROP_GZ_DIST:
      self->gz_dist = g_value_dup_object (value);
      
      if ((self->gz_dist != NULL) && (self->wl_dist != NULL))
      {
        g_assert_cmpuint (nc_galaxy_wl_dist_len (self->wl_dist), ==, nc_galaxy_redshift_len (self->gz_dist));
        self->len = nc_galaxy_wl_dist_len (self->wl_dist);
      }
      
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWL *gwl                = NC_GALAXY_WL (object);
  NcGalaxyWLPrivate * const self = gwl->priv;
  
  g_return_if_fail (NC_IS_GALAXY_WL (object));
  
  switch (prop_id)
  {
    case PROP_WL_DIST:
      g_value_set_object (value, self->wl_dist);
      break;
    case PROP_GZ_DIST:
      g_value_set_object (value, self->gz_dist);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_wl_dispose (GObject *object)
{
  NcGalaxyWL *gwl                = NC_GALAXY_WL (object);
  NcGalaxyWLPrivate * const self = gwl->priv;
  
  nc_galaxy_wl_dist_clear (&self->wl_dist);
  nc_galaxy_redshift_clear (&self->gz_dist);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_class_init (NcGalaxyWLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_nc_galaxy_wl_set_property;
  object_class->get_property = &_nc_galaxy_wl_get_property;
  object_class->dispose      = &_nc_galaxy_wl_dispose;
  object_class->finalize     = &_nc_galaxy_wl_finalize;
  
  /**
   * NcGalaxyWL:wl-dist:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_WL_DIST,
                                   g_param_spec_object ("wl-dist",
                                                        NULL,
                                                        "Weak Lensing distribution",
                                                        NC_TYPE_GALAXY_WL_DIST,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcGalaxyWL:gz-dist:
   *
   * FIXME
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_GZ_DIST,
                                   g_param_spec_object ("gz-dist",
                                                        NULL,
                                                        "Galaxy redshift distribution",
                                                        NC_TYPE_GALAXY_REDSHIFT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_galaxy_wl_new:
 * @wl_dist: a #NcGalaxyWLDist
 * @gz_dist: a #NcGalaxyRedshift
 *
 * Creates a new galaxy weak lensing object.
 *
 * Returns: (transfer full): a new NcGalaxyWL.
 */
NcGalaxyWL *
nc_galaxy_wl_new (NcGalaxyWLDist *wl_dist, NcGalaxyRedshift *gz_dist)
{
  NcGalaxyWL *gwl = g_object_new (NC_TYPE_GALAXY_WL,
                                  "wl-dist", wl_dist,
                                  "gz-dist", gz_dist,
                                  NULL);
  
  return gwl;
}

/**
 * nc_galaxy_wl_ref:
 * @gwl: a #NcGalaxyWL
 *
 * Increase the reference of @gwl by one.
 *
 * Returns: (transfer full): @gwl.
 */
NcGalaxyWL *
nc_galaxy_wl_ref (NcGalaxyWL *gwl)
{
  return g_object_ref (gwl);
}

/**
 * nc_galaxy_wl_free:
 * @gwl: a #NcGalaxyWL
 *
 * Decrease the reference count of @gwl by one.
 *
 */
void
nc_galaxy_wl_free (NcGalaxyWL *gwl)
{
  g_object_unref (gwl);
}

/**
 * nc_galaxy_wl_clear:
 * @gwl: a #NcGalaxyWL
 *
 * Decrease the reference count of @gwl by one, and sets the pointer *@gwl to
 * NULL.
 *
 */
void
nc_galaxy_wl_clear (NcGalaxyWL **gwl)
{
  g_clear_object (gwl);
}

typedef struct _NcGalaxyWLEval
{
  NcGalaxyWLDist *gwld;
  NcHICosmo *cosmo;
  NcHaloDensityProfile *dp;
  NcWLSurfaceMassDensity *smd;
  const gdouble z_cluster;
  gint gal_i;
} NcGalaxyWLEval;

static gdouble
_nc_galaxy_wl_Pz_integ (const gdouble z, gpointer userdata)
{
  NcGalaxyWLEval *gwleval = (NcGalaxyWLEval *) userdata;
  
  return nc_galaxy_wl_dist_m2lnP (gwleval->gwld, gwleval->cosmo, gwleval->dp, gwleval->smd, gwleval->z_cluster, gwleval->gal_i, z);
}

/**
 * nc_galaxy_wl_eval_m2lnP:
 * @gwl: a #NcGalaxyWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 *
 * Computes the galaxy probability given the theoretical modeling.
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_galaxy_wl_eval_m2lnP (NcGalaxyWL *gwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLPrivate * const self = gwl->priv;
  NcGalaxyWLEval gwleval         = {self->wl_dist, cosmo, dp, smd, z_cluster, 0};
  gdouble res                    = 0.0;
  gint gal_i;
  
  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    gwleval.gal_i = gal_i;
    nc_galaxy_wl_dist_m2lnP_prep (self->wl_dist, cosmo, dp, smd, z_cluster, gal_i);
    res += nc_galaxy_redshift_compute_mean_m2lnf (self->gz_dist, gal_i, &_nc_galaxy_wl_Pz_integ, &gwleval);
  }
  
  return res;
}

/**
 * nc_galaxy_wl_len:
 * @gwl: a #NcGalaxyWL
 *
 * Returns: the number of galaxies in @gwl.
 */
guint
nc_galaxy_wl_len (NcGalaxyWL *gwl)
{
  NcGalaxyWLPrivate * const self = gwl->priv;
  
  return self->len;
}


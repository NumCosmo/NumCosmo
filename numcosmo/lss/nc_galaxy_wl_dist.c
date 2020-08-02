/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_dist.c
 *
 *  Mon July 27 11:12:53 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_dist.c
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
 * SECTION:nc_galaxy_wl_dist
 * @title: NcGalaxyWLDist
 * @short_description: Abstract class describing galaxy wl_dists.
 *
 * Abstract class used to define a generic galaxy weak lensing observable 
 * probability distribution $P_\mathrm{wl}(g)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_wl_dist.h"

struct _NcGalaxyWLDistPrivate
{
  gint placeholder;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLDist, nc_galaxy_wl_dist, G_TYPE_OBJECT);

static void
nc_galaxy_wl_dist_init (NcGalaxyWLDist *gwld)
{
  gwld->priv = nc_galaxy_wl_dist_get_instance_private (gwld);
}

static void
_nc_galaxy_wl_dist_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_dist_parent_class)->finalize (object);
}

static gdouble
_nc_galaxy_wl_dist_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, guint gal_i, const gdouble z)
{
  g_error ("_nc_galaxy_wl_dist_m2lnP: method not implemented.");
  
  return 0.0;
}

static gdouble
_nc_galaxy_wl_dist_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  g_error ("_nc_galaxy_wl_dist_gen: method not implemented.");
  
  return 0.0;
}

static guint
_nc_galaxy_wl_dist_len (NcGalaxyWLDist *gwld)
{
  g_error ("_nc_galaxy_wl_dist_len: method not implemented.");

  return 0;
}

static void
nc_galaxy_wl_dist_class_init (NcGalaxyWLDistClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->finalize = &_nc_galaxy_wl_dist_finalize;
  
  klass->m2lnP           = &_nc_galaxy_wl_dist_m2lnP;
  klass->gen             = &_nc_galaxy_wl_dist_gen;
  klass->len             = &_nc_galaxy_wl_dist_len;
}

/**
 * nc_galaxy_wl_dist_ref:
 * @gwld: a #NcGalaxyWLDist
 *
 * Increase the reference of @gwld by one.
 *
 * Returns: (transfer full): @gwld.
 */
NcGalaxyWLDist *
nc_galaxy_wl_dist_ref (NcGalaxyWLDist *gwld)
{
  return g_object_ref (gwld);
}

/**
 * nc_galaxy_wl_dist_free:
 * @gwld: a #NcGalaxyWLDist
 *
 * Decrease the reference count of @gwld by one.
 *
 */
void
nc_galaxy_wl_dist_free (NcGalaxyWLDist *gwld)
{
  g_object_unref (gwld);
}

/**
 * nc_galaxy_wl_dist_clear:
 * @gwld: a #NcGalaxyWLDist
 *
 * Decrease the reference count of @gwld by one, and sets the pointer *@gwld to
 * NULL.
 *
 */
void
nc_galaxy_wl_dist_clear (NcGalaxyWLDist **gwld)
{
  g_clear_object (gwld);
}

/**
 * nc_galaxy_wl_dist_m2lnP_prep: (virtual m2lnP_prep)
 * @gwld: a #NcGalaxyWLDist
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 *
 * Prepare to compute nc_galaxy_wl_dist_m2lnP() at different redshifts.
 *
 */
/**
 * nc_galaxy_wl_dist_m2lnP: (virtual m2lnP)
 * @gwld: a #NcGalaxyWLDist
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 * @z: source redshift
 *
 * Returns: the probability density at @g, $-2\ln\left[P_\mathrm{wl}(g)\right]$.
 */
/**
 * nc_galaxy_wl_dist_gen: (virtual gen)
 * @gwld: a #NcGalaxyWLDist
 * @g_true: true value of the observable
 * @rng: a #NcmRNG
 *
 * Generates a g value from the distribution using @rng.
 *
 * Returns: the generated value $g$.
 */
/**
 * nc_galaxy_wl_dist_len: (virtual len)
 * @gwld: a #NcGalaxyWLDist
 *
 * Returns: the number of galaxies in the object.
 */

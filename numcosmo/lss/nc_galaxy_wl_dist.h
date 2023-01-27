/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_wl_dist.h
 *
 *  Mon July 27 11:13:56 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_wl_dist.h
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

#ifndef _NC_GALAXY_WL_DIST_H_
#define _NC_GALAXY_WL_DIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>
#include <numcosmo/lss/nc_galaxy_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_DIST             (nc_galaxy_wl_dist_get_type ())
#define NC_GALAXY_WL_DIST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_WL_DIST, NcGalaxyWLDist))
#define NC_GALAXY_WL_DIST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_WL_DIST, NcGalaxyWLDistClass))
#define NC_IS_GALAXY_WL_DIST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_WL_DIST))
#define NC_IS_GALAXY_WL_DIST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_WL_DIST))
#define NC_GALAXY_WL_DIST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_WL_DIST, NcGalaxyWLDistClass))

typedef struct _NcGalaxyWLDistClass NcGalaxyWLDistClass;
typedef struct _NcGalaxyWLDist NcGalaxyWLDist;
typedef struct _NcGalaxyWLDistPrivate NcGalaxyWLDistPrivate;

struct _NcGalaxyWLDistClass
{
  /*< private >*/
  GObjectClass parent_class;
  
  void (*m2lnP_initial_prep) (NcGalaxyWLDist *gwld, NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
  void (*m2lnP_prep) (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i);
  gdouble (*m2lnP) (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
  gdouble (*gen) (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
  guint (*len) (NcGalaxyWLDist *gwld);
};

struct _NcGalaxyWLDist
{
  /*< private >*/
  GObject parent_instance;
  NcGalaxyWLDistPrivate *priv;
};

GType nc_galaxy_wl_dist_get_type (void) G_GNUC_CONST;

NcGalaxyWLDist *nc_galaxy_wl_dist_ref (NcGalaxyWLDist *gwld);

void nc_galaxy_wl_dist_free (NcGalaxyWLDist *gwld);
void nc_galaxy_wl_dist_clear (NcGalaxyWLDist **gwld);

NCM_INLINE void nc_galaxy_wl_dist_m2lnP_initial_prep (NcGalaxyWLDist *gwld,NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
NCM_INLINE void nc_galaxy_wl_dist_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i);
NCM_INLINE gdouble nc_galaxy_wl_dist_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z);
NCM_INLINE gdouble nc_galaxy_wl_dist_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng);
NCM_INLINE guint nc_galaxy_wl_dist_len (NcGalaxyWLDist *gwld);

G_END_DECLS

#endif /* _NC_GALAXY_WL_DIST_H_ */

#ifndef _NC_GALAXY_WL_DIST_INLINE_H_
#define _NC_GALAXY_WL_DIST_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
nc_galaxy_wl_dist_m2lnP_initial_prep (NcGalaxyWLDist *gwld,NcGalaxyRedshift *gz, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster)
{
  NcGalaxyWLDistClass *gwld_class = NC_GALAXY_WL_DIST_GET_CLASS (gwld);
  
  if (gwld_class->m2lnP_initial_prep != NULL)
    return NC_GALAXY_WL_DIST_GET_CLASS (gwld)->m2lnP_initial_prep (gwld, gz, cosmo, dp, smd, z_cluster);
}

NCM_INLINE void
nc_galaxy_wl_dist_m2lnP_prep (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i)
{
  NcGalaxyWLDistClass *gwld_class = NC_GALAXY_WL_DIST_GET_CLASS (gwld);
  
  if (gwld_class->m2lnP_prep != NULL)
    return NC_GALAXY_WL_DIST_GET_CLASS (gwld)->m2lnP_prep (gwld, cosmo, dp, smd, z_cluster, gal_i);
}

NCM_INLINE gdouble
nc_galaxy_wl_dist_m2lnP (NcGalaxyWLDist *gwld, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const guint gal_i, const gdouble z)
{
  return NC_GALAXY_WL_DIST_GET_CLASS (gwld)->m2lnP (gwld, cosmo, dp, smd, z_cluster, gal_i, z);
}

NCM_INLINE gdouble
nc_galaxy_wl_dist_gen (NcGalaxyWLDist *gwld, const gdouble g_true, NcmRNG *rng)
{
  return NC_GALAXY_WL_DIST_GET_CLASS (gwld)->gen (gwld, g_true, rng);
}

NCM_INLINE guint
nc_galaxy_wl_dist_len (NcGalaxyWLDist *gwld)
{
  return NC_GALAXY_WL_DIST_GET_CLASS (gwld)->len (gwld);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_GALAXY_WL_DIST_INLINE_H_ */


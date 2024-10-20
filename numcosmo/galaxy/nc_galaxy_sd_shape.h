/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape.h
 *
 *  Sat May 21 22:00:06 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape.h
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
#ifndef _NC_GALAXY_SD_SHAPE_H_
#define _NC_GALAXY_SD_SHAPE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>


G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE             (nc_galaxy_sd_shape_get_type ())
#define NC_GALAXY_SD_SHAPE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SD_SHAPE, NcGalaxySDShape))
#define NC_GALAXY_SD_SHAPE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SD_SHAPE, NcGalaxySDShapeClass))
#define NC_IS_GALAXY_SD_SHAPE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SD_SHAPE))
#define NC_IS_GALAXY_SD_SHAPE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SD_SHAPE))
#define NC_GALAXY_SD_SHAPE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SD_SHAPE, NcGalaxySDShapeClass))

typedef struct _NcGalaxySDShapeClass NcGalaxySDShapeClass;
typedef struct _NcGalaxySDShape NcGalaxySDShape;
typedef struct _NcGalaxySDShapePrivate NcGalaxySDShapePrivate;

struct _NcGalaxySDShapeClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*gen) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, NcmRNG *rng, const gdouble r, const gdouble z, gdouble *et, gdouble *ex);
  gdouble (*integ) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble r, const gdouble z, const gdouble et, const gdouble ex);
  void (*integ_optzs_prep) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble R);
  gdouble (*integ_optzs) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z_source, const gdouble et, const gdouble ex);
};

struct _NcGalaxySDShape
{
  /*< private >*/
  GObject parent_instance;
  NcGalaxySDShapePrivate *priv;
};

GType nc_galaxy_sd_shape_get_type (void) G_GNUC_CONST;

NcGalaxySDShape *nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds);

void nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds);
void nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds);

void nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, NcmRNG *rng, gdouble r, gdouble z, gdouble *et, gdouble *ex);
gdouble nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble r, const gdouble z, const gdouble et, const gdouble ex);
void nc_galaxy_sd_shape_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble R);
gdouble nc_galaxy_sd_shape_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z_source, const gdouble et, const gdouble ex);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_H_ */


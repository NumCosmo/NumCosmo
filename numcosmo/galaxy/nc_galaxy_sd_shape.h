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
#include <numcosmo/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/nc_enum_types.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDShape, nc_galaxy_sd_shape, NC, GALAXY_SD_SHAPE, NcmModel)

struct _NcGalaxySDShapeClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gboolean (*gen) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, NcmRNG *rng, NcGalaxyWLObsCoord coord, const gdouble ra, const gdouble dec, const gdouble z, gdouble *e1, gdouble *e2);
  gdouble (*integ) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, gdouble z, NcmVector *data);
  void (*integ_optzs_prep) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, NcmVector *data);
  gdouble (*integ_optzs) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, gdouble z, NcmVector *data);
  gboolean (*prepare) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hc, NcGalaxyWLObsCoord coord, NcmObjArray *data, NcmObjArray *data_prep);
  GStrv (*get_header) (NcGalaxySDShape *gsds);
  guint (*get_vec_size) (NcGalaxySDShape *gsds);
  gboolean (*set_models) (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hc);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[10];
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_shape);

NcGalaxySDShape *nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds);

void nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds);
void nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds);

gboolean nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, NcmRNG *rng, NcGalaxyWLObsCoord coord, const gdouble ra, const gdouble dec, const gdouble z, gdouble *e1, gdouble *e2);
gdouble nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, gdouble z, NcmVector *data);
void nc_galaxy_sd_shape_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, NcmVector *data);
gdouble nc_galaxy_sd_shape_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, gdouble z, NcmVector *data);
gboolean nc_galaxy_sd_shape_prepare (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hc, NcGalaxyWLObsCoord coord, NcmObjArray *data, NcmObjArray *data_prep);

GStrv nc_galaxy_sd_shape_get_header (NcGalaxySDShape *gsds);
guint nc_galaxy_sd_shape_get_vec_size (NcGalaxySDShape *gsds);
gboolean nc_galaxy_sd_shape_set_models (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloPosition *hc);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_H_ */


/***************************************************************************
 *            nc_galaxy_halo_wl.h
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_halo_wl.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_HALO_WL_H_
#define _NC_GALAXY_HALO_WL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/nc/lss/wl/nc_wl_cluster_lens.h>
#include <numcosmo/nc/lss/wl/nc_wl_surface_mass_density.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile.h>
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_position_factor.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_HALO_WL      (nc_galaxy_halo_wl_get_type ())
#define NC_TYPE_GALAXY_HALO_WL_DATA (nc_galaxy_halo_wl_data_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyHaloWL, nc_galaxy_halo_wl, NC, GALAXY_HALO_WL, GObject)

typedef struct _NcGalaxyHaloWLData NcGalaxyHaloWLData;

GType nc_galaxy_halo_wl_data_get_type (void) G_GNUC_CONST;

NcGalaxyHaloWLData *nc_galaxy_halo_wl_data_new (void);
NcGalaxyHaloWLData *nc_galaxy_halo_wl_data_ref (NcGalaxyHaloWLData *data);
void nc_galaxy_halo_wl_data_unref (NcGalaxyHaloWLData *data);
gdouble nc_galaxy_halo_wl_data_get_radius (NcGalaxyHaloWLData *data);
void nc_galaxy_halo_wl_data_get_rotation (NcGalaxyHaloWLData *data, gdouble *rot_re, gdouble *rot_im);
void nc_galaxy_halo_wl_data_get_sigma (NcGalaxyHaloWLData *data, gdouble *sigma, gdouble *mean_sigma);

NcGalaxyHaloWL *nc_galaxy_halo_wl_new (NcWLClusterLens *lens);
NcGalaxyHaloWL *nc_galaxy_halo_wl_ref (NcGalaxyHaloWL *hwl);
void nc_galaxy_halo_wl_free (NcGalaxyHaloWL *hwl);
void nc_galaxy_halo_wl_clear (NcGalaxyHaloWL **hwl);

NcWLClusterLens *nc_galaxy_halo_wl_peek_lens (NcGalaxyHaloWL *hwl);

void nc_galaxy_halo_wl_prepare (NcGalaxyHaloWL *hwl, NcmMSet *mset);
guint64 nc_galaxy_halo_wl_get_radius_hash (NcGalaxyHaloWL *hwl);
guint64 nc_galaxy_halo_wl_get_profile_hash (NcGalaxyHaloWL *hwl);
void nc_galaxy_halo_wl_update_data_radius (NcGalaxyHaloWL *hwl, NcGalaxyHaloWLData *data, NcGalaxyPositionFactorData *pos_data);
void nc_galaxy_halo_wl_update_data_profile (NcGalaxyHaloWL *hwl, NcHaloDensityProfile *dp, NcGalaxyHaloWLData *data);

void nc_galaxy_halo_wl_prepare_bound (NcGalaxyHaloWL *hwl, NcmMSet *mset);
gboolean nc_galaxy_halo_wl_bound_ready (NcGalaxyHaloWL *hwl);
gboolean nc_galaxy_halo_wl_bound_holds (NcGalaxyHaloWL *hwl);
gdouble nc_galaxy_halo_wl_get_radius_min (NcGalaxyHaloWL *hwl, NcGalaxyHaloWLData *data);
void nc_galaxy_halo_wl_eval_bound (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble R_min, const gdouble z_s, gdouble *kappa, gdouble *g_t);
void nc_galaxy_halo_wl_set_shape_bound (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, NcGalaxyHaloWLData *data, const gdouble z_max, NcGalaxyShapeFactorData *s_data);

G_END_DECLS

#endif /* _NC_GALAXY_HALO_WL_H_ */


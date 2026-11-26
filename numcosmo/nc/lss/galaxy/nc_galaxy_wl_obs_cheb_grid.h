/***************************************************************************
 *            nc_galaxy_wl_obs_cheb_grid.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_cheb_grid.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_GALAXY_WL_OBS_CHEB_GRID_H_
#define _NC_GALAXY_WL_OBS_CHEB_GRID_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_OBS_CHEB_GRID (nc_galaxy_wl_obs_cheb_grid_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLObsChebGrid, nc_galaxy_wl_obs_cheb_grid, NC, GALAXY_WL_OBS_CHEB_GRID, GObject)

NcGalaxyWLObsChebGrid *nc_galaxy_wl_obs_cheb_grid_new (const guint len, const guint n_weights, const guint n_nodes, const gdouble z_cl, const gchar *z_type, const gchar *data_hash);
NcGalaxyWLObsChebGrid *nc_galaxy_wl_obs_cheb_grid_ref (NcGalaxyWLObsChebGrid *grid);
void nc_galaxy_wl_obs_cheb_grid_free (NcGalaxyWLObsChebGrid *grid);
void nc_galaxy_wl_obs_cheb_grid_clear (NcGalaxyWLObsChebGrid **grid);

guint nc_galaxy_wl_obs_cheb_grid_len (NcGalaxyWLObsChebGrid *grid);
guint nc_galaxy_wl_obs_cheb_grid_get_n_weights (NcGalaxyWLObsChebGrid *grid);
guint nc_galaxy_wl_obs_cheb_grid_get_n_nodes (NcGalaxyWLObsChebGrid *grid);
gdouble nc_galaxy_wl_obs_cheb_grid_get_z_cl (NcGalaxyWLObsChebGrid *grid);
const gchar *nc_galaxy_wl_obs_cheb_grid_peek_z_type (NcGalaxyWLObsChebGrid *grid);
const gchar *nc_galaxy_wl_obs_cheb_grid_peek_data_hash (NcGalaxyWLObsChebGrid *grid);

void nc_galaxy_wl_obs_cheb_grid_set_galaxy (NcGalaxyWLObsChebGrid *grid, const guint i, const guint n_total, GArray *weights, NcmVector *z_nodes, const gdouble norm, const gdouble bg_norm, const gdouble fg_norm);
guint nc_galaxy_wl_obs_cheb_grid_get_galaxy (NcGalaxyWLObsChebGrid *grid, const guint i, GArray **weights, NcmVector **z_nodes, gdouble *norm, gdouble *bg_norm, gdouble *fg_norm);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_CHEB_GRID_H_ */


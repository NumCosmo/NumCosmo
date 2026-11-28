/***************************************************************************
 *            nc_galaxy_wl_obs_binned.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_binned.h
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

#ifndef _NC_GALAXY_WL_OBS_BINNED_H_
#define _NC_GALAXY_WL_OBS_BINNED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_catalog.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_OBS_BINNED (nc_galaxy_wl_obs_binned_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLObsBinned, nc_galaxy_wl_obs_binned, NC, GALAXY_WL_OBS_BINNED, NcmCatalog)

/* Columns present in every binned profile */
#define NC_GALAXY_WL_OBS_BINNED_COL_R_LO "r_lo"
#define NC_GALAXY_WL_OBS_BINNED_COL_R_HI "r_hi"
#define NC_GALAXY_WL_OBS_BINNED_COL_RADIUS "radius"
#define NC_GALAXY_WL_OBS_BINNED_COL_VALUE "value"
#define NC_GALAXY_WL_OBS_BINNED_COL_OFFSET "offset"
#define NC_GALAXY_WL_OBS_BINNED_COL_N_GAL "n_gal"

/* Columns of the members catalog (one row per galaxy in a bin) */
#define NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_ROW "row"
#define NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_BIN "bin"
#define NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_RADIUS "radius"
#define NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_FACTOR "factor"

NcGalaxyWLObsBinned * nc_galaxy_wl_obs_binned_new (guint n_bins, GStrv extra_cols);
NcGalaxyWLObsBinned *nc_galaxy_wl_obs_binned_new_from_arrays (GArray *radius, GArray *value, NcmMatrix *cov, GStrv extra_cols);
NcGalaxyWLObsBinned *nc_galaxy_wl_obs_binned_ref (NcGalaxyWLObsBinned *binned);
void nc_galaxy_wl_obs_binned_free (NcGalaxyWLObsBinned *binned);
void nc_galaxy_wl_obs_binned_clear (NcGalaxyWLObsBinned **binned);

guint nc_galaxy_wl_obs_binned_len (NcGalaxyWLObsBinned *binned);

void nc_galaxy_wl_obs_binned_set (NcGalaxyWLObsBinned *binned, const gchar *col, const guint i, const gdouble val);
gdouble nc_galaxy_wl_obs_binned_get (NcGalaxyWLObsBinned *binned, const gchar *col, const guint i);
gboolean nc_galaxy_wl_obs_binned_has_column (NcGalaxyWLObsBinned *binned, const gchar *col);

void nc_galaxy_wl_obs_binned_set_column (NcGalaxyWLObsBinned *binned, const gchar *col, GArray *v);
NcmVector *nc_galaxy_wl_obs_binned_get_column (NcGalaxyWLObsBinned *binned, const gchar *col);

void nc_galaxy_wl_obs_binned_set_cov (NcGalaxyWLObsBinned *binned, NcmMatrix *cov);
NcmMatrix *nc_galaxy_wl_obs_binned_peek_cov (NcGalaxyWLObsBinned *binned);

void nc_galaxy_wl_obs_binned_set_z_cluster (NcGalaxyWLObsBinned *binned, const gdouble z_cluster);
gdouble nc_galaxy_wl_obs_binned_get_z_cluster (NcGalaxyWLObsBinned *binned);

void nc_galaxy_wl_obs_binned_set_estimator (NcGalaxyWLObsBinned *binned, const gchar *estimator);
const gchar *nc_galaxy_wl_obs_binned_peek_estimator (NcGalaxyWLObsBinned *binned);

void nc_galaxy_wl_obs_binned_set_members (NcGalaxyWLObsBinned *binned, NcGalaxyWLObs *source, NcmCatalog *members);
NcGalaxyWLObs *nc_galaxy_wl_obs_binned_peek_source (NcGalaxyWLObsBinned *binned);
NcmCatalog *nc_galaxy_wl_obs_binned_peek_members (NcGalaxyWLObsBinned *binned);
gboolean nc_galaxy_wl_obs_binned_has_members (NcGalaxyWLObsBinned *binned);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_BINNED_H_ */


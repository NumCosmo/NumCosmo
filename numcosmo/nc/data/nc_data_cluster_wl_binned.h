/***************************************************************************
 *            nc_data_cluster_wl_binned.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl_binned.h
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

#ifndef _NC_DATA_CLUSTER_WL_BINNED_H_
#define _NC_DATA_CLUSTER_WL_BINNED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/data/ncm_data.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs_binned.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_binned_estimator.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WL_BINNED (nc_data_cluster_wl_binned_get_type ())

G_DECLARE_FINAL_TYPE (NcDataClusterWLBinned, nc_data_cluster_wl_binned, NC, DATA_CLUSTER_WL_BINNED, NcmData)

NcDataClusterWLBinned *nc_data_cluster_wl_binned_new (NcGalaxyWLBinnedEstimator * est, NcGalaxyWLObsBinned * binned);
NcDataClusterWLBinned *nc_data_cluster_wl_binned_new_from_obs (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObs *obs, GArray *edges);
NcDataClusterWLBinned *nc_data_cluster_wl_binned_ref (NcDataClusterWLBinned *dcwlb);
void nc_data_cluster_wl_binned_free (NcDataClusterWLBinned *dcwlb);
void nc_data_cluster_wl_binned_clear (NcDataClusterWLBinned **dcwlb);

void nc_data_cluster_wl_binned_set_binned (NcDataClusterWLBinned *dcwlb, NcGalaxyWLObsBinned *binned);
NcGalaxyWLObsBinned *nc_data_cluster_wl_binned_peek_binned (NcDataClusterWLBinned *dcwlb);
NcGalaxyWLBinnedEstimator *nc_data_cluster_wl_binned_peek_estimator (NcDataClusterWLBinned *dcwlb);

void nc_data_cluster_wl_binned_set_exact (NcDataClusterWLBinned *dcwlb, gboolean exact);
gboolean nc_data_cluster_wl_binned_get_exact (NcDataClusterWLBinned *dcwlb);

void nc_data_cluster_wl_binned_eval_model (NcDataClusterWLBinned *dcwlb, NcmMSet *mset, NcmVector *model);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WL_BINNED_H_ */


/***************************************************************************
 *            nc_data_cluster_wl_factor.h
 *
 *  Sun Jul 5 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_wl_factor.h
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_DATA_CLUSTER_WL_FACTOR_H_
#define _NC_DATA_CLUSTER_WL_FACTOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/data/ncm_data.h>
#include <numcosmo/nc/data/nc_data_cluster_wl.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_position_factor.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_factor.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WL_FACTOR (nc_data_cluster_wl_factor_get_type ())

G_DECLARE_FINAL_TYPE (NcDataClusterWLFactor, nc_data_cluster_wl_factor, NC, DATA_CLUSTER_WL_FACTOR, NcmData)

NcDataClusterWLFactor *nc_data_cluster_wl_factor_new (NcGalaxyPositionFactor * position_factor, NcGalaxyRedshiftFactor * redshift_factor, NcGalaxyShapeFactor * shape_factor);
NcDataClusterWLFactor *nc_data_cluster_wl_factor_ref (NcDataClusterWLFactor *dcwlf);

void nc_data_cluster_wl_factor_free (NcDataClusterWLFactor *dcwlf);
void nc_data_cluster_wl_factor_clear (NcDataClusterWLFactor **dcwlf);

void nc_data_cluster_wl_factor_set_prec (NcDataClusterWLFactor *dcwlf, gdouble prec);
void nc_data_cluster_wl_factor_set_obs (NcDataClusterWLFactor *dcwlf, NcGalaxyWLObs *obs);
void nc_data_cluster_wl_factor_set_cut (NcDataClusterWLFactor *dcwlf, const gdouble r_min, const gdouble r_max);

void nc_data_cluster_wl_factor_set_integ_method (NcDataClusterWLFactor *dcwlf, NcDataClusterWLIntegMethod integ_method);
NcDataClusterWLIntegMethod nc_data_cluster_wl_factor_get_integ_method (NcDataClusterWLFactor *dcwlf);
void nc_data_cluster_wl_factor_set_n_nodes (NcDataClusterWLFactor *dcwlf, guint n_nodes);
guint nc_data_cluster_wl_factor_get_n_nodes (NcDataClusterWLFactor *dcwlf);
void nc_data_cluster_wl_factor_set_rule_n (NcDataClusterWLFactor *dcwlf, guint rule_n);
guint nc_data_cluster_wl_factor_get_rule_n (NcDataClusterWLFactor *dcwlf);
void nc_data_cluster_wl_factor_set_auto_nodes (NcDataClusterWLFactor *dcwlf, gboolean auto_nodes);
gboolean nc_data_cluster_wl_factor_get_auto_nodes (NcDataClusterWLFactor *dcwlf);
void nc_data_cluster_wl_factor_set_node_reltol (NcDataClusterWLFactor *dcwlf, gdouble node_reltol);
gdouble nc_data_cluster_wl_factor_get_node_reltol (NcDataClusterWLFactor *dcwlf);
void nc_data_cluster_wl_factor_set_max_total_nodes (NcDataClusterWLFactor *dcwlf, guint max_total_nodes);
guint nc_data_cluster_wl_factor_get_max_total_nodes (NcDataClusterWLFactor *dcwlf);

NcGalaxyWLObs *nc_data_cluster_wl_factor_peek_obs (NcDataClusterWLFactor *dcwlf);

void nc_data_cluster_wl_factor_set_resample_flag (NcDataClusterWLFactor *dcwlf, NcDataClusterWLResampleFlag resample_flag);
NcDataClusterWLResampleFlag nc_data_cluster_wl_factor_get_resample_flag (NcDataClusterWLFactor *dcwlf);

GPtrArray *nc_data_cluster_wl_factor_peek_data_array (NcDataClusterWLFactor *dcwlf);

void nc_data_cluster_wl_factor_eval_m2lnP_gal (NcDataClusterWLFactor *dcwlf, NcmMSet *mset, NcmVector *m2lnP_gal);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WL_FACTOR_H_ */


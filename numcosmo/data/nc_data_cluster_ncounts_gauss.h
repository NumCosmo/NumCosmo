/***************************************************************************
 *            nc_data_cluster_ncounts_gauss.h
 *
 *  Thu Apr 22 15:31:19 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_CLUSTER_NCOUNTS_GAUSS_H_
#define _NC_DATA_CLUSTER_NCOUNTS_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/lss/nc_cluster_abundance.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_NCOUNTS_GAUSS             (nc_data_cluster_ncounts_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcDataClusterNCountsGauss, nc_data_cluster_ncounts_gauss, NC, DATA_CLUSTER_NCOUNTS_GAUSS, NcmDataGaussCov)


NcDataClusterNCountsGauss *nc_data_cluster_ncounts_gauss_new (NcClusterAbundance *cad);
void nc_data_cluster_ncounts_gauss_free (NcDataClusterNCountsGauss *ncounts_gauss);


void nc_data_cluster_ncounts_gauss_set_z_obs (NcDataClusterNCountsGauss *ncounts_gauss, NcmVector *z_obs);
void nc_data_cluster_ncounts_gauss_set_z_obs_params (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *z_obs_params);
void nc_data_cluster_ncounts_gauss_set_lnM_obs (NcDataClusterNCountsGauss *ncounts_gauss, NcmVector *lnM_obs);
void nc_data_cluster_ncounts_gauss_set_lnM_obs_params (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *lnM_obs_params);
void nc_data_cluster_ncounts_gauss_set_has_ssc (NcDataClusterNCountsGauss *ncounts_gauss, gboolean on);
void nc_data_cluster_ncounts_gauss_set_s_matrix (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *s_matrix);
void nc_data_cluster_ncounts_gauss_set_resample_s_matrix (NcDataClusterNCountsGauss *ncounts_gauss, NcmMatrix *s_matrix);
void nc_data_cluster_ncounts_gauss_set_fix_cov (NcDataClusterNCountsGauss *ncounts_gauss, gboolean on);


NcmVector *nc_data_cluster_ncounts_gauss_get_z_obs (NcDataClusterNCountsGauss *ncounts_gauss);
NcmMatrix *nc_data_cluster_ncounts_gauss_get_z_obs_params (NcDataClusterNCountsGauss *ncounts_gauss);
NcmVector *nc_data_cluster_ncounts_gauss_get_lnM_obs (NcDataClusterNCountsGauss *ncounts_gauss);
NcmMatrix *nc_data_cluster_ncounts_gauss_get_lnM_obs_params (NcDataClusterNCountsGauss *ncounts_gauss);
NcmMatrix *nc_data_cluster_ncounts_gauss_get_s_matrix (NcDataClusterNCountsGauss *ncounts_gauss);
NcmMatrix *nc_data_cluster_ncounts_gauss_get_resample_s_matrix (NcDataClusterNCountsGauss *ncounts_gauss);
gboolean nc_data_cluster_ncounts_gauss_get_has_ssc (NcDataClusterNCountsGauss *ncounts_gauss);
gboolean nc_data_cluster_ncounts_gauss_get_fix_cov (NcDataClusterNCountsGauss *ncounts_gauss);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_NCOUNTS_GAUSS_H_ */


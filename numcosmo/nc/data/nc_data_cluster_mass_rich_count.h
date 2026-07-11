/***************************************************************************
 *            nc_data_cluster_mass_rich_count.h
 *
 *  Fri Jul 10 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_mass_rich_count.h
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

#ifndef _NC_DATA_CLUSTER_MASS_RICH_COUNT_H_
#define _NC_DATA_CLUSTER_MASS_RICH_COUNT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/ncm/data/ncm_data.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_MASS_RICH_COUNT (nc_data_cluster_mass_rich_count_get_type ())

G_DECLARE_FINAL_TYPE (NcDataClusterMassRichCount, nc_data_cluster_mass_rich_count, NC, DATA_CLUSTER_MASS_RICH_COUNT, NcmData);

NcDataClusterMassRichCount *nc_data_cluster_mass_rich_count_new (void);
NcDataClusterMassRichCount *nc_data_cluster_mass_rich_count_ref (NcDataClusterMassRichCount *dmrc);
void nc_data_cluster_mass_rich_count_free (NcDataClusterMassRichCount *dmrc);
void nc_data_cluster_mass_rich_count_clear (NcDataClusterMassRichCount **dmrc);
void nc_data_cluster_mass_rich_count_apply_cut (NcDataClusterMassRichCount *dmrc, guint N_min);

void nc_data_cluster_mass_rich_count_set_data (NcDataClusterMassRichCount *dmrc, NcmVector *lnM, NcmVector *z, NcmVector *N);
NcmVector *nc_data_cluster_mass_rich_count_peek_lnM (NcDataClusterMassRichCount *dmrc);
NcmVector *nc_data_cluster_mass_rich_count_peek_z (NcDataClusterMassRichCount *dmrc);
NcmVector *nc_data_cluster_mass_rich_count_peek_N (NcDataClusterMassRichCount *dmrc);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_MASS_RICH_COUNT_H_ */


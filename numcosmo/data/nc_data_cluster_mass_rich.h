/***************************************************************************
 *            nc_data_cluster_mass_rich.h
 *
 *  Tue Oct 24 22:22:00 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_mass_rich.c
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_CLUSTER_MASS_RICH_H_
#define _NC_DATA_CLUSTER_MASS_RICH_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_MASS_RICH (nc_data_cluster_mass_rich_get_type ())

G_DECLARE_FINAL_TYPE (NcDataClusterMassRich, nc_data_cluster_mass_rich, NC, DATA_CLUSTER_MASS_RICH, NcmData);

NcDataClusterMassRich *nc_data_cluster_mass_rich_new (void);
NcDataClusterMassRich *nc_data_cluster_mass_rich_ref (NcDataClusterMassRich *dmr);
void nc_data_cluster_mass_rich_free (NcDataClusterMassRich *dmr);
void nc_data_cluster_mass_rich_clear (NcDataClusterMassRich **dmr);

void nc_data_cluster_mass_rich_set_data (NcDataClusterMassRich *dmr, NcmVector *lnM, NcmVector *z, NcmVector *lnR);


G_END_DECLS

#endif /* _NC_DATA_CLUSTER_MASS_RICH_H_ */
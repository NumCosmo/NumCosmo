/***************************************************************************
 *           nc_data_cluster_wl.h
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 *  Tue Jun 15 16:24:17 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.h
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_CLUSTER_WL_H_
#define _NC_DATA_CLUSTER_WL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>
#include <numcosmo/galaxy/nc_galaxy_sd_obs_redshift.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WL (nc_data_cluster_wl_get_type ())

G_DECLARE_FINAL_TYPE (NcDataClusterWL, nc_data_cluster_wl, NC, DATA_CLUSTER_WL, NcmData);

typedef enum _NcDataClusterWLResampleFlag /*< flags,underscore_name=NC_DATA_CLUSTER_WL_RESAMPLE_FLAG >*/
{
  NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_POSITION = 1 << 0,
  NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_REDSHIFT = 1 << 1,
  NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_SHAPE    = 1 << 2,
  NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL      = (1 << 3) - 1,
} NcDataClusterWLResampleFlag;

typedef struct _NcDataClusterWLPrivate NcDataClusterWLPrivate;

NcDataClusterWL *nc_data_cluster_wl_new (void);
NcDataClusterWL *nc_data_cluster_wl_ref (NcDataClusterWL *dcwl);

void nc_data_cluster_wl_free (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_clear (NcDataClusterWL **dcwl);

void nc_data_cluster_wl_set_prec (NcDataClusterWL *dcwl, gdouble prec);
void nc_data_cluster_wl_set_obs (NcDataClusterWL *dcwl, NcGalaxyWLObs *obs);
void nc_data_cluster_wl_set_cut (NcDataClusterWL *dcwl, const gdouble r_min, const gdouble r_max);
NcGalaxyWLObs *nc_data_cluster_wl_peek_obs (NcDataClusterWL *dcwl);

void nc_data_cluster_wl_set_resample_flag (NcDataClusterWL *dcwl, NcDataClusterWLResampleFlag resample_flag);
NcDataClusterWLResampleFlag nc_data_cluster_wl_get_resample_flag (NcDataClusterWL *dcwl);

NcmObjArray *nc_data_cluster_wl_peek_data_array (NcDataClusterWL *dcwl);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WL_H_ */


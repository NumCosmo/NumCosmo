/***************************************************************************
 *           nc_data_cluster_wl.h
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 *  Tue Jun 15 16:24:17 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.h
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_DATA_CLUSTER_WL_H_
#define _NC_DATA_CLUSTER_WL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_stats_dist.h>
#include <numcosmo/math/ncm_stats_dist_kde.h>
#include <numcosmo/math/ncm_stats_dist_kernel_gauss.h>
#include <numcosmo/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>
#include <numcosmo/galaxy/nc_galaxy_sd_obs_redshift.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>


G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WL             (nc_data_cluster_wl_get_type ())
G_DECLARE_FINAL_TYPE (NcDataClusterWL, nc_data_cluster_wl, NC, DATA_CLUSTER_WL, NcmData);

typedef struct _NcDataClusterWLPrivate NcDataClusterWLPrivate;

struct _NcDataClusterWLClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

struct _NcDataClusterWL
{
  /*< private >*/
  NcmData parent_instance;
  NcDataClusterWLPrivate *priv;
};


NcDataClusterWL *nc_data_cluster_wl_new (NcGalaxySDShape *s_dist, NcGalaxySDObsRedshift *z_dist, NcGalaxySDPosition *p_dist);
NcDataClusterWL *nc_data_cluster_wl_ref (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_prepare_kde (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd);
gdouble nc_data_cluster_wl_kde_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *m2lnP_gal);
gdouble nc_data_cluster_wl_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, NcmVector *m2lnP_gal);
void nc_data_cluster_wl_free (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_clear (NcDataClusterWL **dcwl);
void nc_data_cluster_wl_set_use_kde (NcDataClusterWL *dcwl, gboolean kde);
void nc_data_cluster_wl_set_prec (NcDataClusterWL *dcwl, gdouble prec);
void nc_data_cluster_wl_set_ndata (NcDataClusterWL *dcwl, gdouble ndata);
void nc_data_cluster_wl_set_obs (NcDataClusterWL *dcwl, NcGalaxyWLObs *obs);
void nc_data_cluster_wl_set_cut (NcDataClusterWL *dcwl, const gdouble theta_min, const gdouble theta_max);
void nc_data_cluster_wl_gen_obs (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hc, guint nobs, NcmRNG *rng);
NcGalaxyWLObs *nc_data_cluster_wl_peek_obs (NcDataClusterWL *dcwl);
NcmStatsDist *nc_data_cluster_wl_peek_kde (NcDataClusterWL *dcwl);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WL_H_ */


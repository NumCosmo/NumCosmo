/***************************************************************************
 *           nc_data_cluster_wl.h
 *
 *  Tue Jun 15 16:24:17 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.h
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
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>
#include <numcosmo/galaxy/nc_galaxy_sd_z_proxy.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>


G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_WL             (nc_data_cluster_wl_get_type ())
#define NC_DATA_CLUSTER_WL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_WL, NcDataClusterWL))
#define NC_DATA_CLUSTER_WL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_WL, NcDataClusterWLClass))
#define NC_IS_DATA_CLUSTER_WL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_WL))
#define NC_IS_DATA_CLUSTER_WL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_WL))
#define NC_DATA_CLUSTER_WL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_WL, NcDataClusterWLClass))

typedef struct _NcDataClusterWLClass NcDataClusterWLClass;
typedef struct _NcDataClusterWL NcDataClusterWL;
typedef struct _NcDataClusterWLPrivate NcDataClusterWLPrivate;

struct _NcDataClusterWLClass
{
  /*< private >*/
  NcmDataClass parent_class;
};

/**
 * NcDataClusterWLObs:
 * @NC_DATA_CLUSTER_WL_ZCLUSTER: cluster redshift
 * @NC_DATA_CLUSTER_WL_GOBS: measured reduced shear
 * @NC_DATA_CLUSTER_WL_PZ: redshift distribution (photometric)
 *
 */
typedef enum _NcDataClusterWLObs
{
  NC_DATA_CLUSTER_WL_ZCLUSTER = 0,
  NC_DATA_CLUSTER_WL_GOBS,
  NC_DATA_CLUSTER_WL_PZ,
  /* < private > */
  NC_DATA_CLUSTER_WL_LEN, /*< skip >*/
} NcDataClusterWLObs;

struct _NcDataClusterWL
{
  /*< private >*/
  NcmData parent_instance;
  NcDataClusterWLPrivate *priv;
};

GType nc_data_cluster_wl_get_type (void) G_GNUC_CONST;

NcDataClusterWL *nc_data_cluster_wl_new (NcGalaxySDShape *s_dist, NcGalaxySDZProxy *zp_dist, NcGalaxySDPosition *rz_dist);
NcDataClusterWL *nc_data_cluster_wl_new_from_file (const gchar *filename);
NcDataClusterWL *nc_data_cluster_wl_ref (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_prepare_kde (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster);
gdouble nc_data_cluster_wl_kde_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *m2lnP_gal);
gdouble nc_data_cluster_wl_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcmVector *m2lnP_gal);
void nc_data_cluster_wl_free (NcDataClusterWL *dcwl);
void nc_data_cluster_wl_clear (NcDataClusterWL **dcwl);
void nc_data_cluster_wl_set_use_kde (NcDataClusterWL *dcwl, gboolean kde);
void nc_data_cluster_wl_set_prec (NcDataClusterWL *dcwl, gdouble prec);
void nc_data_cluster_wl_set_ndata (NcDataClusterWL *dcwl, gdouble ndata);
void nc_data_cluster_wl_set_obs (NcDataClusterWL *dcwl, NcmMatrix *obs);
void nc_data_cluster_wl_set_cut (NcDataClusterWL *dcwl, const gdouble r_min, const gdouble r_max);
void nc_data_cluster_wl_gen_obs (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, guint nobs, NcmRNG *rng);
NcmMatrix *nc_data_cluster_wl_peek_obs (NcDataClusterWL *dcwl);
NcmStatsDist *nc_data_cluster_wl_peek_kde (NcDataClusterWL *dcwl);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_WL_H_ */


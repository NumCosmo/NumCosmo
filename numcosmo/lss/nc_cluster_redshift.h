/***************************************************************************
 *            nc_cluster_redshift.h
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_CLUSTER_REDSHIFT_H_
#define _NC_CLUSTER_REDSHIFT_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_REDSHIFT             (nc_cluster_redshift_get_type ())
#define NC_CLUSTER_REDSHIFT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_REDSHIFT, NcClusterRedshift))
#define NC_CLUSTER_REDSHIFT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_REDSHIFT, NcClusterRedshiftClass))
#define NC_IS_CLUSTER_REDSHIFT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_REDSHIFT))
#define NC_IS_CLUSTER_REDSHIFT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_REDSHIFT))
#define NC_CLUSTER_REDSHIFT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_REDSHIFT, NcClusterRedshiftClass))

typedef struct _NcClusterRedshiftClass NcClusterRedshiftClass;
typedef struct _NcClusterRedshift NcClusterRedshift;

/**
 * NcClusterRedshiftImpl:
 * @NC_CLUSTER_REDSHIFT_P: FIXME
 * @NC_CLUSTER_REDSHIFT_INTP: FIXME
 * @NC_CLUSTER_REDSHIFT_RESAMPLE: FIXME
 * @NC_CLUSTER_REDSHIFT_P_LIMITS: FIXME
 * @NC_CLUSTER_REDSHIFT_N_LIMTS: FIXME
 * 
 */ 
typedef enum _NcClusterRedshiftImpl
{
  NC_CLUSTER_REDSHIFT_P        = 1 << 0,
  NC_CLUSTER_REDSHIFT_INTP     = 1 << 1,
  NC_CLUSTER_REDSHIFT_RESAMPLE = 1 << 2,
  NC_CLUSTER_REDSHIFT_P_LIMITS = 1 << 3,
  NC_CLUSTER_REDSHIFT_N_LIMTS  = 1 << 4,
} NcClusterRedshiftImpl;

#define NC_CLUSTER_REDSHIFT_IMPL_ALL (~0)

struct _NcClusterRedshift
{
  /*< private >*/
  GObject parent_instance;
};

struct _NcClusterRedshiftClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*P) (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params);
  gdouble (*intP) (NcClusterRedshift *clusterz, gdouble lnM, gdouble z);
  gboolean (*resample) (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params);
  void (*P_limits) (NcClusterRedshift *clusterz, gdouble *z_obs, gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
  void (*N_limits) (NcClusterRedshift *clusterz, gdouble *z_lower, gdouble *z_upper);
  guint (*obs_len) (NcClusterRedshift *clusterz);
  guint (*obs_params_len) (NcClusterRedshift *clusterz);
  NcClusterRedshiftImpl impl;
};

GType nc_cluster_redshift_get_type (void) G_GNUC_CONST;

NcClusterRedshift *nc_cluster_redshift_new_from_name (gchar *redshift_name);
NcClusterRedshift *nc_cluster_redshift_ref (NcClusterRedshift *clusterz);
void nc_cluster_redshift_free (NcClusterRedshift *clusterz);
void nc_cluster_redshift_clear (NcClusterRedshift **clusterz);

NcClusterRedshiftImpl nc_cluster_redshift_impl (NcClusterRedshift *clusterz);

guint nc_cluster_redshift_obs_len (NcClusterRedshift *clusterz);
guint nc_cluster_redshift_obs_params_len (NcClusterRedshift *clusterz);

gdouble nc_cluster_redshift_p (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params);
gdouble nc_cluster_redshift_intp (NcClusterRedshift *clusterz, gdouble lnM, gdouble z);
gboolean nc_cluster_redshift_resample (NcClusterRedshift *clusterz, gdouble lnM, gdouble z, gdouble *z_obs, gdouble *z_obs_params);
void nc_cluster_redshift_p_limits (NcClusterRedshift *clusterz, gdouble *z_obs, gdouble *z_obs_params, gdouble *z_lower, gdouble *z_upper);
void nc_cluster_redshift_n_limits (NcClusterRedshift *clusterz, gdouble *z_lower, gdouble *z_upper);

G_END_DECLS

#endif /* _NC_CLUSTER_REDSHIFT_H_ */

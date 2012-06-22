/***************************************************************************
 *            nc_cluster_mass.h
 *
 *  Thu June 21 23:27:30 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_CLUSTER_MASS_H_
#define _NC_CLUSTER_MASS_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS             (nc_cluster_mass_get_type ())
#define NC_CLUSTER_MASS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS, NcClusterMass))
#define NC_CLUSTER_MASS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS, NcClusterMassClass))
#define NC_IS_CLUSTER_MASS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS))
#define NC_IS_CLUSTER_MASS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS))
#define NC_CLUSTER_MASS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS, NcClusterMassClass))

typedef struct _NcClusterMassClass NcClusterMassClass;
typedef struct _NcClusterMass NcClusterMass;

struct _NcClusterMassClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*dist_eval) (NcClusterMass *clusterm, gdouble z, gdouble lnM, gdouble *lnM_obs, gdouble *lnM_obs_params);
  gboolean (*resample) (NcClusterMass *clusterm, gdouble z, gdouble lnM, gdouble *lnM_obs, gdouble *lnM_obs_params);
  void (*integ_limits) (NcClusterMass *clusterm, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
  void (*lnm_limits) (NcClusterMass *clusterm, gdouble *lnM_lower, gdouble *lnM_upper);
  guint (*obs_len) (NcClusterMass *clusterm);
  guint (*obs_params_len) (NcClusterMass *clusterm);
};

struct _NcClusterMass
{
  /*< private >*/
  GObject parent_instance;
};

GType nc_cluster_mass_get_type (void) G_GNUC_CONST;

NcClusterMass *nc_cluster_mass_new_from_name (gchar *mass_name);
NcClusterMass *nc_cluster_mass_ref (NcClusterMass *clusterm);
void nc_cluster_mass_free (NcClusterMass *clusterm);

guint nc_cluster_mass_obs_len (NcClusterMass *clusterm);
guint nc_cluster_mass_obs_params_len (NcClusterMass *clusterm);

gdouble nc_cluster_mass_dist_eval (NcClusterMass *clusterm, gdouble z, gdouble lnM, gdouble *lnM_obs, gdouble *lnM_obs_params);
gboolean nc_cluster_mass_resample (NcClusterMass *clusterm, gdouble z, gdouble lnM, gdouble *lnM_obs, gdouble *lnM_obs_params);
void nc_cluster_mass_integ_limits (NcClusterMass *clusterm, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
void nc_cluster_mass_lnm_limits (NcClusterMass *clusterm, gdouble *lnM_lower, gdouble *lnM_upper);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_H_ */

/***************************************************************************
 *            nc_cluster_mass.h
 *
 *  Thu June 21 23:27:30 2012
 *  Copyright  2012  Mariana Penna Lima, Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com> 
 * Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS            (nc_cluster_mass_get_type ())
#define NC_CLUSTER_MASS(obj)            (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS, NcClusterMass))
#define NC_CLUSTER_MASS_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS, NcClusterMassClass))
#define NC_IS_CLUSTER_MASS(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS))
#define NC_IS_CLUSTER_MASS_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS))
#define NC_CLUSTER_MASS_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS, NcClusterMassClass))

typedef struct _NcClusterMassClass NcClusterMassClass;
typedef struct _NcClusterMass NcClusterMass;

/**
 * NcClusterMassImpl:
 * @NC_CLUSTER_MASS_P: FIXME
 * @NC_CLUSTER_MASS_INTP: FIXME
 * @NC_CLUSTER_MASS_RESAMPLE: FIXME
 * @NC_CLUSTER_MASS_P_LIMITS: FIXME
 * @NC_CLUSTER_MASS_N_LIMITS: FIXME
 * 
 */ 
typedef enum _NcClusterMassImpl
{
  NC_CLUSTER_MASS_P        = 1 << 0,
  NC_CLUSTER_MASS_INTP     = 1 << 1,
  NC_CLUSTER_MASS_RESAMPLE = 1 << 2,
  NC_CLUSTER_MASS_P_LIMITS = 1 << 3,
  NC_CLUSTER_MASS_N_LIMITS = 1 << 4,
} NcClusterMassImpl;

#define NC_CLUSTER_MASS_IMPL_ALL (~0)

struct _NcClusterMassClass
{
  /*< private >*/
  NcmModelClass parent_class;
  gdouble (*P) (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params);
  gdouble (*intP) (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z);
  gboolean (*resample) (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params);
  void (*P_limits) (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
  void (*N_limits) (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper);
  guint (*obs_len) (NcClusterMass *clusterm);
  guint (*obs_params_len) (NcClusterMass *clusterm);
  NcClusterMassImpl impl;
};

struct _NcClusterMass
{
  /*< private >*/
  NcmModel parent_instance;
};

GType nc_cluster_mass_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_cluster_mass);

NcClusterMass *nc_cluster_mass_new_from_name (gchar *mass_name);
NcClusterMass *nc_cluster_mass_ref (NcClusterMass *clusterm);
void nc_cluster_mass_free (NcClusterMass *clusterm);
void nc_cluster_mass_clear (NcClusterMass **clusterm);

NcClusterMassImpl nc_cluster_mass_impl (NcClusterMass *clusterm);

guint nc_cluster_mass_obs_len (NcClusterMass *clusterm);
guint nc_cluster_mass_obs_params_len (NcClusterMass *clusterm);

gdouble nc_cluster_mass_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params);
gdouble nc_cluster_mass_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z);
gboolean nc_cluster_mass_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params);
void nc_cluster_mass_p_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
void nc_cluster_mass_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper);

void nc_cluster_mass_log_all_models ();


G_END_DECLS

#endif /* _NC_CLUSTER_MASS_H_ */

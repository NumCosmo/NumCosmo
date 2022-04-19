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
typedef struct _NcClusterMassPrivate NcClusterMassPrivate;

/**
 * NcClusterMassImpl:
 * @NC_CLUSTER_MASS_P: probability density function of the true-observable cluster masses
 * @NC_CLUSTER_MASS_INTP: probability distribution (integration over the observable mass(es))
 * @NC_CLUSTER_MASS_RESAMPLE: resample function to generate the cluster masses following
 * the underlying cluster mass distribution.
 * @NC_CLUSTER_MASS_P_LIMITS: function to set the lower and upper limits of the to compute
 * the integral of the cluster mass distribution.
 * @NC_CLUSTER_MASS_N_LIMITS: function to set the lower and upper thresholds of
 * the observable cluster mass to compute the normalization of the cluster mass distribution.
 *
 */
typedef enum _NcClusterMassImpl
{
  NC_CLUSTER_MASS_P = 0,
  NC_CLUSTER_MASS_INTP,
  NC_CLUSTER_MASS_RESAMPLE,
  NC_CLUSTER_MASS_P_LIMITS,
  NC_CLUSTER_MASS_N_LIMITS,
} NcClusterMassImpl;

#define NC_CLUSTER_MASS_IMPL_ALL NCM_MODEL_CLASS_IMPL_ALL

struct _NcClusterMassClass
{
  /*< private >*/
  NcmModelClass parent_class;
  
  gdouble (*P) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
  gdouble (*intP) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z);
  gdouble (*intP_bin) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *lnM_bin_lower, const gdouble *lnM_bin_upper, const gdouble *lnM_obs_params);
  gboolean (*resample) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
  void (*P_limits) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
  void (*P_bin_limits) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
  void (*N_limits) (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
  gdouble (*volume) (NcClusterMass *clusterm);
  void (*P_vec_z_lnMobs) (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, GArray *res);
  guint obs_len;
  guint obs_params_len;
};

struct _NcClusterMass
{
  /*< private >*/
  NcmModel parent_instance;
  NcClusterMassPrivate *priv;
};

GType nc_cluster_mass_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_cluster_mass);

guint nc_cluster_mass_class_obs_len (NcClusterMassClass *clusterm_class);
guint nc_cluster_mass_class_obs_params_len (NcClusterMassClass *clusterm_class);

NcClusterMass *nc_cluster_mass_new_from_name (gchar *mass_name);
NcClusterMass *nc_cluster_mass_ref (NcClusterMass *clusterm);
void nc_cluster_mass_free (NcClusterMass *clusterm);
void nc_cluster_mass_clear (NcClusterMass **clusterm);

gdouble nc_cluster_mass_p (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
gdouble nc_cluster_mass_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z);
gdouble nc_cluster_mass_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
gboolean nc_cluster_mass_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
void nc_cluster_mass_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
void nc_cluster_mass_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
void nc_cluster_mass_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
gdouble nc_cluster_mass_volume (NcClusterMass *clusterm);

void nc_cluster_mass_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, GArray *res);

guint nc_cluster_mass_obs_len (NcClusterMass *clusterm);
guint nc_cluster_mass_obs_params_len (NcClusterMass *clusterm);

void nc_cluster_mass_log_all_models (void);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_H_ */

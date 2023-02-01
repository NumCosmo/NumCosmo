/***************************************************************************
 *            nc_cor_cluster_cmb_lens_limber.h
 *
 *  Wed June 11 17:19:50 2014
 *  Copyright  2014  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_cor_cluster_cmb_lens_limber.h
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_COR_CLUSTER_CMB_LENS_LIMBER_H_
#define _NC_COR_CLUSTER_CMB_LENS_LIMBER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>
//#include <numcosmo/lss/nc_matter_var.h>
#include <numcosmo/lss/nc_growth_func.h>
//#include <numcosmo/lss/nc_mass_function.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_halo_bias.h>
#include <numcosmo/lss/nc_cluster_abundance.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER             (nc_cor_cluster_cmb_lens_limber_get_type ())
#define NC_COR_CLUSTER_CMB_LENS_LIMBER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER, NcCorClusterCmbLensLimber))
#define NC_COR_CLUSTER_CMB_LENS_LIMBER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER, NcCorClusterCmbLensLimberClass))
#define NC_IS_COR_CLUSTER_CMB_LENS_LIMBER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER))
#define NC_IS_COR_CLUSTER_CMB_LENS_LIMBER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER))
#define NC_COR_CLUSTER_CMB_LENS_LIMBER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER, NcCorClusterCmbLensLimberClass))

typedef struct _NcCorClusterCmbLensLimberClass NcCorClusterCmbLensLimberClass;
typedef struct _NcCorClusterCmbLensLimber NcCorClusterCmbLensLimber;

struct _NcCorClusterCmbLensLimberClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcCorClusterCmbLensLimber
{
  /*< private >*/
  GObject parent_instance;
  NcmSpline *oneh_int_mass_spline;
};

GType nc_cor_cluster_cmb_lens_limber_get_type (void) G_GNUC_CONST;

NcCorClusterCmbLensLimber *nc_cor_cluster_cmb_lens_limber_new (void);

gdouble nc_cor_cluster_cmb_lens_limber_oneh_int_mass (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcClusterMass *clusterm, NcHICosmo *cosmo, NcHaloDensityProfile *dp, gdouble k, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params);
gdouble nc_cor_cluster_cmb_lens_limber_oneh_term (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcHICosmo *cosmo, NcDistance *dist, NcHaloDensityProfile *dp, gint l, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *z_obs, gdouble *z_obs_params);
gdouble nc_cor_cluster_cmb_lens_limber_twoh_int_mass1 (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble z);
gdouble nc_cor_cluster_cmb_lens_limber_twoh_int_mm (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcHICosmo *cosmo, NcHaloDensityProfile *dp, gdouble k, gdouble z);
gdouble nc_cor_cluster_cmb_lens_limber_twoh_int_mass2 (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcClusterMass *clusterm, NcHICosmo *cosmo, NcHaloDensityProfile *dp, gdouble k, gdouble z);
gdouble nc_cor_cluster_cmb_lens_limber_twoh_term (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcHICosmo *cosmo, NcDistance *dist, NcHaloDensityProfile *dp, gint l, gdouble *z_obs, gdouble *z_obs_params);

G_END_DECLS

#endif /* _NC_COR_CLUSTER_CMB_LENS_LIMBER_H_ */


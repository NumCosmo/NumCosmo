/***************************************************************************
 *            nc_cluster_mass_benson_xray.h
 *
 *  Tue July 9 17:01:24 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_CLUSTER_MASS_BENSON_XRAY_H_
#define _NC_CLUSTER_MASS_BENSON_XRAY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_BENSON_XRAY             (nc_cluster_mass_benson_xray_get_type ())
#define NC_CLUSTER_MASS_BENSON_XRAY(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_BENSON_XRAY, NcClusterMassBensonXRay))
#define NC_CLUSTER_MASS_BENSON_XRAY_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_BENSON_XRAY, NcClusterMassBensonXRayClass))
#define NC_IS_CLUSTER_MASS_BENSON_XRAY(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_BENSON_XRAY))
#define NC_IS_CLUSTER_MASS_BENSON_XRAY_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_BENSON_XRAY))
#define NC_CLUSTER_MASS_BENSON_XRAY_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_BENSON_XRAY, NcClusterMassBensonXRayClass))

typedef struct _NcClusterMassBensonXRayClass NcClusterMassBensonXRayClass;
typedef struct _NcClusterMassBensonXRay NcClusterMassBensonXRay;

/**
 * NcClusterMassBensonXRayParams:
 * @NC_CLUSTER_MASS_BENSON_XRAY_A_SZ: normalization of the SZ mass-observable relation
 * @NC_CLUSTER_MASS_BENSON_XRAY_B_SZ: FIXME
 * @NC_CLUSTER_MASS_BENSON_XRAY_C_SZ: FIXME
 * @NC_CLUSTER_MASS_BENSON_XRAY_D_SZ: standard deviation of the mass-observable relation 
 * @NC_CLUSTER_MASS_BENSON_XRAY_A_X: normalization of the X-ray mass-observable relation
 * @NC_CLUSTER_MASS_BENSON_XRAY_B_X: FIXME
 * @NC_CLUSTER_MASS_BENSON_XRAY_C_X: FIXME
 * @NC_CLUSTER_MASS_BENSON_XRAY_D_X: standard deviation of the X-ray mass-observable relation  
 *
 * FIXME
 */
typedef enum _NcClusterMassBensonXRayParams
{
  NC_CLUSTER_MASS_BENSON_XRAY_A_SZ = 0,
  NC_CLUSTER_MASS_BENSON_XRAY_B_SZ,
  NC_CLUSTER_MASS_BENSON_XRAY_C_SZ,
  NC_CLUSTER_MASS_BENSON_XRAY_D_SZ,
  NC_CLUSTER_MASS_BENSON_XRAY_A_X,
  NC_CLUSTER_MASS_BENSON_XRAY_B_X,
  NC_CLUSTER_MASS_BENSON_XRAY_C_X,
  NC_CLUSTER_MASS_BENSON_XRAY_D_X, /*< private >*/
  NC_CLUSTER_MASS_BENSON_XRAY_SPARAM_LEN, /*< skip >*/
} NcClusterMassBensonXRayParams;

#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_A_SZ  (5.58)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_B_SZ  (1.32)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_C_SZ  (0.87)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_D_SZ  (0.24)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_A_X  (5.77)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_B_X  (0.57)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_C_X  (-0.4)
#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_D_X  (0.12)

#define NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL (0.0)


struct _NcClusterMassBensonXRayClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassBensonXRay
{
  /*< private >*/
  NcClusterMass parent_instance;
  gdouble signif_obs_min;
  gdouble signif_obs_max;  
  gdouble z0;
  gdouble M0; 
};

GType nc_cluster_mass_benson_xray_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_BENSON_XRAY_H_ */

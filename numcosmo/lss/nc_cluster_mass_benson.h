/***************************************************************************
 *            nc_cluster_mass_benson.h
 *
 *  Tue July 9 14:18:11 2012
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

#ifndef _NC_CLUSTER_MASS_BENSON_H_
#define _NC_CLUSTER_MASS_BENSON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_BENSON             (nc_cluster_mass_benson_get_type ())
#define NC_CLUSTER_MASS_BENSON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_BENSON, NcClusterMassBenson))
#define NC_CLUSTER_MASS_BENSON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_BENSON, NcClusterMassBensonClass))
#define NC_IS_CLUSTER_MASS_BENSON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_BENSON))
#define NC_IS_CLUSTER_MASS_BENSON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_BENSON))
#define NC_CLUSTER_MASS_BENSON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_BENSON, NcClusterMassBensonClass))

typedef struct _NcClusterMassBensonClass NcClusterMassBensonClass;
typedef struct _NcClusterMassBenson NcClusterMassBenson;

/**
 * NcClusterMassBensonSParams:
 * @NC_CLUSTER_MASS_BENSON_A_SZ: normalization of the mass-observable relation
 * @NC_CLUSTER_MASS_BENSON_B_SZ: FIXME
 * @NC_CLUSTER_MASS_BENSON_C_SZ: FIXME
 * @NC_CLUSTER_MASS_BENSON_D_SZ: standard deviation of the mass-observable relation 
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_BENSON_SPARAMS >*/
{
  NC_CLUSTER_MASS_BENSON_A_SZ = 0,
  NC_CLUSTER_MASS_BENSON_B_SZ,
  NC_CLUSTER_MASS_BENSON_C_SZ,
  NC_CLUSTER_MASS_BENSON_D_SZ, 
  /* < private > */
  NC_CLUSTER_MASS_BENSON_SPARAM_LEN, /*< skip >*/
} NcClusterMassBensonSParams;

#define NC_CLUSTER_MASS_BENSON_DEFAULT_A_SZ  (5.58)
#define NC_CLUSTER_MASS_BENSON_DEFAULT_B_SZ  (1.32)
#define NC_CLUSTER_MASS_BENSON_DEFAULT_C_SZ  (0.87)
#define NC_CLUSTER_MASS_BENSON_DEFAULT_D_SZ  (0.24)

#define NC_CLUSTER_MASS_BENSON_DEFAULT_PARAMS_ABSTOL (0.0)

#define NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND (1.0e13)
#define NC_CLUSTER_MASS_BENSON_XI_ZETA_DIST_CUT (2.0)

struct _NcClusterMassBensonClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassBenson
{
  /*< private >*/
  NcClusterMass parent_instance;
  gdouble signif_obs_min;
  gdouble signif_obs_max;  
  gdouble z0;
  gdouble M0;
};

GType nc_cluster_mass_benson_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_BENSON_H_ */

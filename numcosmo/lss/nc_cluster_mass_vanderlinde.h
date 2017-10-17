/***************************************************************************
 *            nc_cluster_mass_vanderlinde.h
 *
 *  Tue July 3 15:21:05 2012
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

#ifndef _NC_CLUSTER_MASS_VANDERLINDE_H_
#define _NC_CLUSTER_MASS_VANDERLINDE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_VANDERLINDE             (nc_cluster_mass_vanderlinde_get_type ())
#define NC_CLUSTER_MASS_VANDERLINDE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_VANDERLINDE, NcClusterMassVanderlinde))
#define NC_CLUSTER_MASS_VANDERLINDE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_VANDERLINDE, NcClusterMassVanderlindeClass))
#define NC_IS_CLUSTER_MASS_VANDERLINDE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_VANDERLINDE))
#define NC_IS_CLUSTER_MASS_VANDERLINDE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_VANDERLINDE))
#define NC_CLUSTER_MASS_VANDERLINDE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_VANDERLINDE, NcClusterMassVanderlindeClass))

typedef struct _NcClusterMassVanderlindeClass NcClusterMassVanderlindeClass;
typedef struct _NcClusterMassVanderlinde NcClusterMassVanderlinde;

/**
 * NcClusterMassVanderlindeSParams:
 * @NC_CLUSTER_MASS_VANDERLINDE_A_SZ: normalization of the mass-observable relation
 * @NC_CLUSTER_MASS_VANDERLINDE_B_SZ: FIXME
 * @NC_CLUSTER_MASS_VANDERLINDE_C_SZ: FIXME
 * @NC_CLUSTER_MASS_VANDERLINDE_D_SZ: standard deviation of the mass-observable relation 
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_VANDERLINDE_SPARAMS >*/
{
  NC_CLUSTER_MASS_VANDERLINDE_A_SZ = 0,
  NC_CLUSTER_MASS_VANDERLINDE_B_SZ,
  NC_CLUSTER_MASS_VANDERLINDE_C_SZ,
  NC_CLUSTER_MASS_VANDERLINDE_D_SZ, 
  /* < private > */
  NC_CLUSTER_MASS_VANDERLINDE_SPARAM_LEN, /*< skip >*/
} NcClusterMassVanderlindeSParams;

#define NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_A_SZ  (6.01)
#define NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_B_SZ  (1.31)
#define NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_C_SZ  (1.6)
#define NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_D_SZ  (0.21)

#define NC_CLUSTER_MASS_VANDERLINDE_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterMassVanderlindeClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassVanderlinde
{
  /*< private >*/
  NcClusterMass parent_instance; 
  gdouble signif_obs_min;
  gdouble signif_obs_max;  
  gdouble z0;
  gdouble M0;
};

GType nc_cluster_mass_vanderlinde_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_VANDERLINDE_H_ */

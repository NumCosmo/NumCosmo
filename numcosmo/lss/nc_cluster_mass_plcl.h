/***************************************************************************
 *            nc_cluster_mass_plcl.h
 *
 *  Sun Mar 1 21:58:11 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2015 <pennalima@gmail.com>
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

#ifndef _NC_CLUSTER_MASS_PLCL_H_
#define _NC_CLUSTER_MASS_PLCL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_PLCL             (nc_cluster_mass_plcl_get_type ())
#define NC_CLUSTER_MASS_PLCL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_PLCL, NcClusterMassPlCL))
#define NC_CLUSTER_MASS_PLCL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_PLCL, NcClusterMassPlCLClass))
#define NC_IS_CLUSTER_MASS_PLCL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_PLCL))
#define NC_IS_CLUSTER_MASS_PLCL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_PLCL))
#define NC_CLUSTER_MASS_PLCL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_PLCL, NcClusterMassPlCLClass))

typedef struct _NcClusterMassPlCLClass NcClusterMassPlCLClass;
typedef struct _NcClusterMassPlCL NcClusterMassPlCL;

/**
 * NcClusterMassPlclParams:
 * @NC_CLUSTER_MASS_PLCL_A_SZ: slope of the mass-SZ relation
 * @NC_CLUSTER_MASS_PLCL_B_SZ: FIXME
 * @NC_CLUSTER_MASS_PLCL_SD_SZ: standard deviation of the mass-SZ relation
 * @NC_CLUSTER_MASS_PLCL_A_L: slope of the mass-lensing relation
 * @NC_CLUSTER_MASS_PLCL_B_L: FIXME
 * @NC_CLUSTER_MASS_PLCL_SD_L: standard deviation of the mass-lensing relation 
 * @NC_CLUSTER_MASS_PLCL_COR: correlation coefficient between the SZ and lensing masses
 * @NC_CLUSTER_MASS_PLCL_MCUT: lower mass cut-off 
 *
 * FIXME
 */
typedef enum _NcClusterMassPlCLParams
{
  NC_CLUSTER_MASS_PLCL_A_SZ = 0,
  NC_CLUSTER_MASS_PLCL_B_SZ,
  NC_CLUSTER_MASS_PLCL_SD_SZ,
  NC_CLUSTER_MASS_PLCL_A_L,
  NC_CLUSTER_MASS_PLCL_B_L,
  NC_CLUSTER_MASS_PLCL_SD_L, 
  NC_CLUSTER_MASS_PLCL_COR, 
  NC_CLUSTER_MASS_PLCL_MCUT, /*< private >*/
  NC_CLUSTER_MASS_PLCL_SPARAM_LEN, /*< skip >*/
} NcClusterMassPlclParams;

#define NC_CLUSTER_MASS_PLCL_DEFAULT_A_SZ  (5.58)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_B_SZ  (1.32)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_SD_SZ (1.32)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_A_L  (0.87)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_B_L  (0.24)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_SD_L (1.32)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_COR  (1.32)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_MCUT  (1.32)

#define NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL (0.0)

//#define NC_CLUSTER_MASS_PLCL_M_LOWER_BOUND (1.0e13)

struct _NcClusterMassPlCL
{
  /*< private >*/
  NcClusterMass parent_instance;
  //gdouble sz_obs_min;
  //gdouble sz_obs_max;
  //gdouble lens_obs_min;
  //gdouble lens_obs_max;
  //gdouble M0;
};

struct _NcClusterMassPlCLClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

GType nc_cluster_mass_plcl_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_PLCL_H_ */

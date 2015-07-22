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
 * NcClusterMassPlCLParams:
 * @NC_CLUSTER_MASS_PLCL_A_SZ: slope of the mass-SZ relation
 * @NC_CLUSTER_MASS_PLCL_B_SZ: FIXME
 * @NC_CLUSTER_MASS_PLCL_SD_SZ: standard deviation of the mass-SZ relation
 * @NC_CLUSTER_MASS_PLCL_A_L: slope of the mass-lensing relation
 * @NC_CLUSTER_MASS_PLCL_B_L: FIXME
 * @NC_CLUSTER_MASS_PLCL_SD_L: standard deviation of the mass-lensing relation 
 * @NC_CLUSTER_MASS_PLCL_COR: correlation coefficient between the SZ and lensing masses 
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
  NC_CLUSTER_MASS_PLCL_COR, /*< private >*/
  NC_CLUSTER_MASS_PLCL_SPARAM_LEN, /*< skip >*/
} NcClusterMassPlCLParams;

#define NC_CLUSTER_MASS_PLCL_DEFAULT_A_SZ  (1.041)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_B_SZ  (0.17)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_SD_SZ (0.234)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_A_L  (0.9737)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_B_L  (0.0)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_SD_L (0.4)
#define NC_CLUSTER_MASS_PLCL_DEFAULT_COR  (0.488)

#define NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL (0.0)

#define NC_CLUSTER_MASS_PLCL_MPL (0)
#define NC_CLUSTER_MASS_PLCL_SD_PL (0)
#define NC_CLUSTER_MASS_PLCL_MCL (1)
#define NC_CLUSTER_MASS_PLCL_SD_CL (1)

struct _NcClusterMassPlCL
{
  /*< private >*/
  NcClusterMass parent_instance;
  gdouble M0;
  gdouble *workz;
};

struct _NcClusterMassPlCLClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

GType nc_cluster_mass_plcl_get_type (void) G_GNUC_CONST;

void nc_cluster_mass_plcl_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, NcClusterMassPlCL *mszl, gdouble lnM, const gdouble *Mobs, const gdouble *Mobs_params);
void nc_cluster_mass_plcl_levmar_f_new_variables (gdouble *p, gdouble *hx, gint m, gint n, NcClusterMassPlCL *mszl, gdouble lnM, const gdouble *Mobs, const gdouble *Mobs_params);
void nc_cluster_mass_plcl_levmar_J_new_variables (gdouble *p, gdouble *j, gint m, gint n, NcClusterMassPlCL *mszl, gdouble lnM, const gdouble *Mobs, const gdouble *Mobs_params);
void nc_cluster_mass_plcl_peak_new_variables (gdouble N, gdouble *lb, gdouble *ub, NcClusterMassPlCL *mszl, gdouble lnM, const gdouble *Mobs, const gdouble *Mobs_params);
gdouble nc_cluster_mass_plcl_pdf (NcClusterMass *clusterm, gdouble lnM500, gdouble lnMsz, gdouble lnMl, const gdouble *Mobs, const gdouble *Mobs_params);
gdouble nc_cluster_mass_plcl_pdf_new_variables (NcClusterMass *clusterm, gdouble lnM_M0, gdouble w1, gdouble w2, const gdouble *Mobs, const gdouble *Mobs_params);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_PLCL_H_ */

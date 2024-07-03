/***************************************************************************
 *            nc_cluster_mass_lnrich_ext.c
 *
 *  Tue Oct 31 16:15:11 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti and Cinthia Nunes de Lima
 *  <vitenti@uel.br>, <cinthia.n.lima@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti and Cinthia Nunes de Lima 2023
 * <vitenti@uel.br>, <cinthia.n.lima@uel.br>
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

#ifndef _NC_CLUSTER_MASS_LNRICH_EXT_H_
#define _NC_CLUSTER_MASS_LNRICH_EXT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_LNRICH_EXT             (nc_cluster_mass_lnrich_ext_get_type ())
#define NC_CLUSTER_MASS_LNRICH_EXT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_LNRICH_EXT, NcClusterMassLnrichExt))
#define NC_CLUSTER_MASS_LNRICH_EXT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_LNRICH_EXT, NcClusterMassLnrichExtClass))
#define NC_IS_CLUSTER_MASS_LNRICH_EXT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_LNRICH_EXT))
#define NC_IS_CLUSTER_MASS_LNRICH_EXT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_LNRICH_EXT))
#define NC_CLUSTER_MASS_LNRICH_EXT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_LNRICH_EXT, NcClusterMassLnrichExtClass))

typedef struct _NcClusterMassLnrichExtClass NcClusterMassLnrichExtClass;
typedef struct _NcClusterMassLnrichExt NcClusterMassLnrichExt;
typedef struct _NcClusterMassLnrichExtPrivate NcClusterMassLnrichExtPrivate;

/**
 * NcClusterMassLnrichExtSParams:
 * @NC_CLUSTER_MASS_LNRICH_EXT_MU: bias of the mean
 * @NC_CLUSTER_MASS_LNRICH_EXT_MU_M1: slope on the mean
 * @NC_CLUSTER_MASS_LNRICH_EXT_MU_Z1: redshift dependency on the mean
 * @NC_CLUSTER_MASS_LNRICH_EXT_MU_M2: quadratic slope on the mean
 * @NC_CLUSTER_MASS_LNRICH_EXT_MU_Z2: quadratic redshift dependency on the mean
 * @NC_CLUSTER_MASS_LNRICH_EXT_MU_MZ: cross term on the mean
 * @NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_0: bias of the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M1: slope on the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z1: redshift dependency standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M2: quadratic slope on the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z2: quadratic redshift dependency on the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_MZ: cross term on the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_LNRICH_EXT_A0: 
 * @NC_CLUSTER_MASS_LNRICH_EXT_CUT: cut in richness
 * @NC_CLUSTER_MASS_LNRICH_EXT_CUT_M1: slope on the cut in richness
 * @NC_CLUSTER_MASS_LNRICH_EXT_CUT_Z1: redshift dependency on the cut in richness
 * Parameters of the extended richness-mass relation.
 *
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_LNRICH_EXT_SPARAMS >*/
{
  NC_CLUSTER_MASS_LNRICH_EXT_MU,
  NC_CLUSTER_MASS_LNRICH_EXT_MU_M1,
  NC_CLUSTER_MASS_LNRICH_EXT_MU_Z1,
  NC_CLUSTER_MASS_LNRICH_EXT_MU_M2,
  NC_CLUSTER_MASS_LNRICH_EXT_MU_Z2,
  NC_CLUSTER_MASS_LNRICH_EXT_MU_MZ,
  NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_0,
  NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M1,
  NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z1,
  NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_M2,
  NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_Z2,
  NC_CLUSTER_MASS_LNRICH_EXT_SIGMA_MZ,
  NC_CLUSTER_MASS_LNRICH_EXT_A0,
  NC_CLUSTER_MASS_LNRICH_EXT_CUT,
  NC_CLUSTER_MASS_LNRICH_EXT_CUT_M1,
  NC_CLUSTER_MASS_LNRICH_EXT_CUT_Z1,
  /* < private > */
  NC_CLUSTER_MASS_LNRICH_EXT_SPARAM_LEN, /*< skip >*/
} NcClusterMassLnrichExtSParams;

#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU  (3.19)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_M1  (2.0 / M_LN10)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_Z1  (-0.7 / M_LN10)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_M2  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_Z2  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_MU_MZ  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_0  (0.33)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_M1  (-0.08 / M_LN10)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_Z1  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_M2  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_Z2  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_SIGMA_MZ  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_A0  (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_CUT       (6)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_PARAMS_ABSTOL (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_CUT_M1 (0.0)
#define NC_CLUSTER_MASS_LNRICH_EXT_DEFAULT_CUT_Z1 (0.0)

struct _NcClusterMassLnrichExtClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassLnrichExt
{
  /*< private >*/
  NcClusterMass parent_instance;
  NcClusterMassLnrichExtPrivate *priv;
};

GType nc_cluster_mass_lnrich_ext_get_type (void) G_GNUC_CONST;

gdouble nc_cluster_mass_lnrich_ext_get_mean_richness (NcClusterMassLnrichExt *lnrich_ext, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_lnrich_ext_get_std_richness (NcClusterMassLnrichExt *lnrich_ext, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_lnrich_ext_get_cut (NcClusterMassLnrichExt *lnrich_ext, gdouble lnM, gdouble z);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_LNRICH_EXT_H_ */


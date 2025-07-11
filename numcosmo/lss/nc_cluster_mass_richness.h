/***************************************************************************
 *            nc_cluster_mass_richness.h
 *
 *  Fri June 22 17:49:11 2012
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

#ifndef _NC_CLUSTER_MASS_RICHNESS_H_
#define _NC_CLUSTER_MASS_RICHNESS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_RICHNESS             (nc_cluster_mass_richness_get_type ())
#define NC_CLUSTER_MASS_RICHNESS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_RICHNESS, NcClusterMassRichness))
#define NC_CLUSTER_MASS_RICHNESS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_RICHNESS, NcClusterMassRichnessClass))
#define NC_IS_CLUSTER_MASS_RICHNESS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_RICHNESS))
#define NC_IS_CLUSTER_MASS_RICHNESS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_RICHNESS))
#define NC_CLUSTER_MASS_RICHNESS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_RICHNESS, NcClusterMassRichnessClass))

typedef struct _NcClusterMassRichnessClass NcClusterMassRichnessClass;
typedef struct _NcClusterMassRichness NcClusterMassRichness;
typedef struct _NcClusterMassRichnessPrivate NcClusterMassRichnessPrivate;

/**
 * NcClusterMassRichnessSParams:
 * @NC_CLUSTER_MASS_RICHNESS_MU: bias of the mean
 * @NC_CLUSTER_MASS_RICHNESS_MU_M: slope on the mean
 * @NC_CLUSTER_MASS_RICHNESS_MU_Z: redshift dependency on the mean
 * @NC_CLUSTER_MASS_RICHNESS_MU_MZ: mixed dependency on the mean
 * @NC_CLUSTER_MASS_RICHNESS_SIGMA_0: mass-richness relation variance
 * 
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_RICHNESS_SPARAMS >*/
{
  NC_CLUSTER_MASS_RICHNESS_MU,
  NC_CLUSTER_MASS_RICHNESS_MU_M,
  NC_CLUSTER_MASS_RICHNESS_MU_Z,
  NC_CLUSTER_MASS_RICHNESS_MU_MZ,
  NC_CLUSTER_MASS_RICHNESS_SIGMA_0,
  /* < private > */
  NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN, /*< skip >*/
} NcClusterMassRichnessSParams;

#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU  (3.19)
#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU_M  (2.0 / M_LN10)
#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU_Z  (-0.7 / M_LN10)
#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_MU_MZ  (0.0)
#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_SIGMA_0  (0.33)
#define NC_CLUSTER_MASS_RICHNESS_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterMassRichnessClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassRichness
{
  /*< private >*/
  NcClusterMass parent_instance;
  NcClusterMassRichnessPrivate *priv;
};

GType nc_cluster_mass_richness_get_type (void) G_GNUC_CONST;
#define NC_CLUSTER_MASS_RICHNESS_SIGMA_OBS (0)

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_RICHNESS_H_ */


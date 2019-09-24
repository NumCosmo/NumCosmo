/***************************************************************************
 *            nc_cluster_mass_ascaso.h
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Mariana Penna Lima and Begoña Ascaso
 *  <pennalima@gmail.com>, <bego.ascaso.work@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima and Begoña Ascaso 2017 <pennalima@gmail.com>, <bego.ascaso.work@gmail.com>
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

#ifndef _NC_CLUSTER_MASS_ASCASO_H_
#define _NC_CLUSTER_MASS_ASCASO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_ASCASO             (nc_cluster_mass_ascaso_get_type ())
#define NC_CLUSTER_MASS_ASCASO(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_ASCASO, NcClusterMassAscaso))
#define NC_CLUSTER_MASS_ASCASO_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_ASCASO, NcClusterMassAscasoClass))
#define NC_IS_CLUSTER_MASS_ASCASO(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_ASCASO))
#define NC_IS_CLUSTER_MASS_ASCASO_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_ASCASO))
#define NC_CLUSTER_MASS_ASCASO_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_ASCASO, NcClusterMassAscasoClass))

typedef struct _NcClusterMassAscasoClass NcClusterMassAscasoClass;
typedef struct _NcClusterMassAscaso NcClusterMassAscaso;

/**
 * NcClusterMassAscasoSParams:
 * @NC_CLUSTER_MASS_ASCASO_P0: bias 
 * @NC_CLUSTER_MASS_ASCASO_P1: slope 
 * @NC_CLUSTER_MASS_ASCASO_P2: redshift dependency
 * @NC_CLUSTER_MASS_ASCASO_SIGMA: standard deviation of the log-normal distribution
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_ASCASO_SPARAMS >*/
{
  NC_CLUSTER_MASS_ASCASO_P0 = 0,
  NC_CLUSTER_MASS_ASCASO_P1,
  NC_CLUSTER_MASS_ASCASO_P2,
  NC_CLUSTER_MASS_ASCASO_SIGMA,  
  /* < private > */
  NC_CLUSTER_MASS_ASCASO_SPARAM_LEN, /*< skip >*/
} NcClusterMassAscasoSParams;

#define NC_CLUSTER_MASS_ASCASO_DEFAULT_P0  (0.0)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_P1  (1.0)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_P2  (1.0)
#define NC_CLUSTER_MASS_ASCASO_DEFAULT_SIGMA  (0.04)

#define NC_CLUSTER_MASS_ASCASO_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterMassAscasoClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassAscaso
{
  /*< private >*/
  NcClusterMass parent_instance;
  gdouble M0;
  gdouble lnRichness_max;
  gdouble lnRichness_min;
};

GType nc_cluster_mass_ascaso_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_ASCASO_H_ */

/***************************************************************************
 *            nc_cluster_mass_lnnormal.h
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

#ifndef _NC_CLUSTER_MASS_LNNORMAL_H_
#define _NC_CLUSTER_MASS_LNNORMAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_LNNORMAL             (nc_cluster_mass_lnnormal_get_type ())
#define NC_CLUSTER_MASS_LNNORMAL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_LNNORMAL, NcClusterMassLnnormal))
#define NC_CLUSTER_MASS_LNNORMAL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_LNNORMAL, NcClusterMassLnnormalClass))
#define NC_IS_CLUSTER_MASS_LNNORMAL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_LNNORMAL))
#define NC_IS_CLUSTER_MASS_LNNORMAL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_LNNORMAL))
#define NC_CLUSTER_MASS_LNNORMAL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_LNNORMAL, NcClusterMassLnnormalClass))

typedef struct _NcClusterMassLnnormalClass NcClusterMassLnnormalClass;
typedef struct _NcClusterMassLnnormal NcClusterMassLnnormal;

/**
 * NcClusterMassLnnormalSParams:
 * @NC_CLUSTER_MASS_LNNORMAL_BIAS: FIXME
 * @NC_CLUSTER_MASS_LNNORMAL_SIGMA: standard deviation of the log-normal distribution
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_LNNORMAL_SPARAMS >*/
{
  NC_CLUSTER_MASS_LNNORMAL_BIAS = 0,
  NC_CLUSTER_MASS_LNNORMAL_SIGMA,  
  /* < private > */
  NC_CLUSTER_MASS_LNNORMAL_SPARAM_LEN, /*< skip >*/
} NcClusterMassLnnormalSParams;

#define NC_CLUSTER_MASS_LNNORMAL_DEFAULT_BIAS  (0.0)
#define NC_CLUSTER_MASS_LNNORMAL_DEFAULT_SIGMA  (0.04)

#define NC_CLUSTER_MASS_LNNORMAL_DEFAULT_PARAMS_ABSTOL (0.0)

struct _NcClusterMassLnnormalClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassLnnormal
{
  /*< private >*/
  NcClusterMass parent_instance;
  gdouble lnMobs_max;
  gdouble lnMobs_min;
};

GType nc_cluster_mass_lnnormal_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_LNNORMAL_H_ */

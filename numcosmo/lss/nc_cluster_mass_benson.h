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

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_BENSON             (nc_cluster_mass_benson_get_type ())
#define NC_CLUSTER_MASS_BENSON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_BENSON, NcClusterMassBenson))
#define NC_CLUSTER_MASS_BENSON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_BENSON, NcClusterMassBensonClass))
#define NC_IS_CLUSTER_MASS_BENSON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_BENSON))
#define NC_IS_CLUSTER_MASS_BENSON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_BENSON))
#define NC_CLUSTER_MASS_BENSON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_BENSON, NcClusterMassBensonClass))

typedef struct _NcClusterMassBensonClass NcClusterMassBensonClass;
typedef struct _NcClusterMassBenson NcClusterMassBenson;

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
  gdouble Asz;
  gdouble Bsz;
  gdouble Csz;
  gdouble Dsz;
};

GType nc_cluster_mass_benson_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_BENSON_H_ */

/***************************************************************************
 *            nc_cluster_mass_nodist.h
 *
 *  Fri June 22 13:43:03 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_CLUSTER_MASS_NODIST_H_
#define _NC_CLUSTER_MASS_NODIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_NODIST             (nc_cluster_mass_nodist_get_type ())
#define NC_CLUSTER_MASS_NODIST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_MASS_NODIST, NcClusterMassNodist))
#define NC_CLUSTER_MASS_NODIST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_MASS_NODIST, NcClusterMassNodistClass))
#define NC_IS_CLUSTER_MASS_NODIST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_MASS_NODIST))
#define NC_IS_CLUSTER_MASS_NODIST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_MASS_NODIST))
#define NC_CLUSTER_MASS_NODIST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_MASS_NODIST, NcClusterMassNodistClass))

typedef struct _NcClusterMassNodistClass NcClusterMassNodistClass;
typedef struct _NcClusterMassNodist NcClusterMassNodist;
typedef struct _NcClusterMassNodistPrivate NcClusterMassNodistPrivate;

struct _NcClusterMassNodistClass
{
  /*< private >*/
  NcClusterMassClass parent_class;
};

struct _NcClusterMassNodist
{
  /*< private >*/
  NcClusterMass parent_instance;
  NcClusterMassNodistPrivate *priv;
};

GType nc_cluster_mass_nodist_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_NODIST_H_ */


/***************************************************************************
 *            nc_data_cluster_poisson.h
 *
 *  Tue Apr  6 01:12:58 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_CLUSTER_POISSON_H_
#define _NC_DATA_CLUSTER_POISSON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_poisson.h>
#include <numcosmo/nc_data_cluster_ncount.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_POISSON             (nc_data_cluster_poisson_get_type ())
#define NC_DATA_CLUSTER_POISSON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_POISSON, NcDataClusterPoisson))
#define NC_DATA_CLUSTER_POISSON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_POISSON, NcDataClusterPoissonClass))
#define NC_IS_DATA_CLUSTER_POISSON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_POISSON))
#define NC_IS_DATA_CLUSTER_POISSON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_POISSON))
#define NC_DATA_CLUSTER_POISSON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_POISSON, NcDataClusterPoissonClass))

typedef struct _NcDataClusterPoissonClass NcDataClusterPoissonClass;
typedef struct _NcDataClusterPoisson NcDataClusterPoisson;

struct _NcDataClusterPoisson
{
  /*< private >*/
  NcmDataPoisson parent_instance;
  NcDataClusterNCount *ncount;
};

struct _NcDataClusterPoissonClass
{
  /*< private >*/
  NcmDataPoissonClass parent_class;
};

GType nc_data_cluster_poisson_get_type (void) G_GNUC_CONST;

NcmData *nc_data_cluster_poisson_new_cad (NcClusterAbundance *cad);
NcmData *nc_data_cluster_poisson_new (NcDataClusterNCount *ncount);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_POISSON_H_ */


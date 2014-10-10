/***************************************************************************
 *            nc_abc_cluster_ncount.h
 *
 *  Mon October 13 13:29:08 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_abc_cluster_ncount.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_ABC_CLUSTER_NCOUNT_H_
#define _NC_ABC_CLUSTER_NCOUNT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_abc.h>

#include <gsl/gsl_histogram2d.h>

G_BEGIN_DECLS

#define NC_TYPE_ABC_CLUSTER_NCOUNT             (nc_abc_cluster_ncount_get_type ())
#define NC_ABC_CLUSTER_NCOUNT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_ABC_CLUSTER_NCOUNT, NcABCClusterNCount))
#define NC_ABC_CLUSTER_NCOUNT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_ABC_CLUSTER_NCOUNT, NcABCClusterNCountClass))
#define NC_IS_ABC_CLUSTER_NCOUNT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_ABC_CLUSTER_NCOUNT))
#define NC_IS_ABC_CLUSTER_NCOUNT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_ABC_CLUSTER_NCOUNT))
#define NC_ABC_CLUSTER_NCOUNT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_ABC_CLUSTER_NCOUNT, NcABCClusterNCountClass))

typedef struct _NcABCClusterNCountClass NcABCClusterNCountClass;
typedef struct _NcABCClusterNCount NcABCClusterNCount;

struct _NcABCClusterNCountClass
{
  /*< private >*/
  NcmABCClass parent_class;

};

struct _NcABCClusterNCount
{
  /*< private >*/
  NcmABC parent_instance;
  gdouble epsilon;
  gsl_histogram2d *data_summary;
  NcmMatrix *covar;
};

GType nc_abc_cluster_ncount_get_type (void) G_GNUC_CONST;

NcABCClusterNCount *nc_abc_cluster_ncount_new (NcmMSet *mset, NcmMSetTransKern *prior, NcmDataset *dset);

G_END_DECLS

#endif /* _NC_ABC_CLUSTER_NCOUNT_H_ */


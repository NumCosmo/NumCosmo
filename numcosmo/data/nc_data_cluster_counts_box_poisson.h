/***************************************************************************
 *            nc_data_cluster_counts_box_poisson.h
 *
 *  Mon Feb 20 15:31:47 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_DATA_CLUSTER_COUNTS_BOX_POISSON_H_
#define _NC_DATA_CLUSTER_COUNTS_BOX_POISSON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_poisson.h>
#include <numcosmo/lss/nc_halo_mass_function.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON             (nc_data_cluster_counts_box_poisson_get_type ())
#define NC_DATA_CLUSTER_COUNTS_BOX_POISSON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON, NcDataClusterCountsBoxPoisson))
#define NC_DATA_CLUSTER_COUNTS_BOX_POISSON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON, NcDataClusterCountsBoxPoissonClass))
#define NC_IS_DATA_CLUSTER_COUNTS_BOX_POISSON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON))
#define NC_IS_DATA_CLUSTER_COUNTS_BOX_POISSON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON))
#define NC_DATA_CLUSTER_COUNTS_BOX_POISSON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON, NcDataClusterCountsBoxPoissonClass))

typedef struct _NcDataClusterCountsBoxPoissonClass NcDataClusterCountsBoxPoissonClass;
typedef struct _NcDataClusterCountsBoxPoisson NcDataClusterCountsBoxPoisson;

struct _NcDataClusterCountsBoxPoissonClass
{
  /*< private >*/
  NcmDataPoissonClass parent_class;
};

struct _NcDataClusterCountsBoxPoisson
{
  /*< private >*/
  NcmDataPoisson parent_instance;
  NcHaloMassFunction *mfp;
  NcmVector *mass_knots;
	gdouble redshift;
  gdouble volume;
	NcmSpline *dndlog10M;
};

GType nc_data_cluster_counts_box_poisson_get_type (void) G_GNUC_CONST;

NcDataClusterCountsBoxPoisson *nc_data_cluster_counts_box_poisson_new (NcHaloMassFunction *mfp);
void nc_data_cluster_counts_box_poisson_init_from_sampling (NcDataClusterCountsBoxPoisson *cpoisson, NcmMSet *mset, NcmVector *mass_knots, const gdouble volume, const gdouble redshift, NcmRNG *rng);

G_END_DECLS

#endif /* _NC_DATA_CLUSTER_COUNTS_BOX_POISSON_H_ */


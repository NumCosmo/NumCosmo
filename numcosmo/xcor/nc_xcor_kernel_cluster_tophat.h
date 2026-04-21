/***************************************************************************
 *            nc_xcor_kernel_cluster_tophat.h
 *
 *  Mon April 21 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_XCOR_KERNEL_CLUSTER_TOPHAT_H_
#define _NC_XCOR_KERNEL_CLUSTER_TOPHAT_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/xcor/nc_xcor_kernel_cluster.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_CLUSTER_TOPHAT (nc_xcor_kernel_cluster_tophat_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelClusterTophat, nc_xcor_kernel_cluster_tophat, NC, XCOR_KERNEL_CLUSTER_TOPHAT, NcXcorKernelCluster)

NcXcorKernelClusterTophat *nc_xcor_kernel_cluster_tophat_new (NcDistance *dist, NcmPowspec *ps, gdouble z_lower, gdouble z_upper);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_CLUSTER_TOPHAT_H_ */

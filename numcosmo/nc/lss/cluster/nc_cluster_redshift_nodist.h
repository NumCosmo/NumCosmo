/***************************************************************************
 *            nc_cluster_redshift_nodist.h
 *
 *  Fri June 22 13:45:04 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NC_CLUSTER_REDSHIFT_NODIST_H_
#define _NC_CLUSTER_REDSHIFT_NODIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_REDSHIFT_NODIST             (nc_cluster_redshift_nodist_get_type ())

G_DECLARE_FINAL_TYPE (NcClusterRedshiftNodist, nc_cluster_redshift_nodist, NC, CLUSTER_REDSHIFT_NODIST, NcClusterRedshift)


G_END_DECLS

#endif /* _NC_CLUSTER_REDSHIFT_NODIST_H_ */


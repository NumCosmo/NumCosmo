/***************************************************************************
 *            nc_cluster_richness_projection.h
 *
 *  Wed September 17 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_cluster_richness_projection.h
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

#ifndef _NC_CLUSTER_RICHNESS_PROJECTION_H_
#define _NC_CLUSTER_RICHNESS_PROJECTION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_RICHNESS_PROJECTION (nc_cluster_richness_projection_get_type ())

G_DECLARE_FINAL_TYPE (NcClusterRichnessProjection, nc_cluster_richness_projection, NC, CLUSTER_RICHNESS_PROJECTION, GObject)

NcClusterRichnessProjection *nc_cluster_richness_projection_new (void);
NcClusterRichnessProjection *nc_cluster_richness_projection_ref (NcClusterRichnessProjection *crp);

void nc_cluster_richness_projection_free (NcClusterRichnessProjection *crp);
void nc_cluster_richness_projection_clear (NcClusterRichnessProjection **crp);

void nc_cluster_richness_projection_set_lnlambda_range (NcClusterRichnessProjection *crp, gdouble lnlambda_min, gdouble lnlambda_max);
void nc_cluster_richness_projection_set_reltol (NcClusterRichnessProjection *crp, gdouble reltol);
gdouble nc_cluster_richness_projection_get_reltol (NcClusterRichnessProjection *crp);

void nc_cluster_richness_projection_prepare (NcClusterRichnessProjection *crp, gdouble mu, gdouble sigma, gdouble tau);

gdouble nc_cluster_richness_projection_eval (NcClusterRichnessProjection *crp, gdouble lnlambda);
gdouble nc_cluster_richness_projection_eval_lnlambda (NcClusterRichnessProjection *crp, gdouble lnlambda);
gdouble nc_cluster_richness_projection_eval_int (NcClusterRichnessProjection *crp, gdouble lnlambda_lo, gdouble lnlambda_hi);

#define NC_CLUSTER_RICHNESS_PROJECTION_DEFAULT_RELTOL (1.0e-11)

G_END_DECLS

#endif /* _NC_CLUSTER_RICHNESS_PROJECTION_H_ */


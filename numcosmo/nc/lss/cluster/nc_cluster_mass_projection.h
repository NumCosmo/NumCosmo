/***************************************************************************
 *            nc_cluster_mass_projection.h
 *
 *  Thu October 02 10:00:00 2026
 *  Copyright  2026  Cinthia N. Lima
 *  <cinthia.n.lima@uel.br>
 ****************************************************************************/
/*
 * nc_cluster_mass_projection.h
 * Copyright (C) 2026 Cinthia N. Lima <cinthia.n.lima@uel.br>
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

#ifndef _NC_CLUSTER_MASS_PROJECTION_H_
#define _NC_CLUSTER_MASS_PROJECTION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_richness.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_PROJECTION (nc_cluster_mass_projection_get_type ())

G_DECLARE_FINAL_TYPE (NcClusterMassProjection, nc_cluster_mass_projection, NC, CLUSTER_MASS_PROJECTION, NcClusterMassRichness)

/**
 * NcClusterMassProjectionSParams:
 * @NC_CLUSTER_MASS_PROJECTION_MU_P0: constant term (bias) in the mean log-richness
 * @NC_CLUSTER_MASS_PROJECTION_MU_P1: linear mass coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_PROJECTION_MU_P2: redshift evolution coefficient in the mean log-richness
 * @NC_CLUSTER_MASS_PROJECTION_SIGMA_P0: constant term (bias) in the standard deviation
 * @NC_CLUSTER_MASS_PROJECTION_SIGMA_P1: linear mass coefficient in the standard deviation
 * @NC_CLUSTER_MASS_PROJECTION_SIGMA_P2: redshift evolution coefficient in the standard deviation
 * @NC_CLUSTER_MASS_PROJECTION_F_PRJ: fraction of clusters affected by projection
 * @NC_CLUSTER_MASS_PROJECTION_TAU: rate of the exponential richness boost
 *
 * Parameters of the mass-richness relation with projection effects.
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_PROJECTION_SPARAMS,prefix=NC_CLUSTER_MASS_PROJECTION >*/
{
  NC_CLUSTER_MASS_PROJECTION_MU_P0 = NC_CLUSTER_MASS_RICHNESS_SPARAM_LEN,
  NC_CLUSTER_MASS_PROJECTION_MU_P1,
  NC_CLUSTER_MASS_PROJECTION_MU_P2,
  NC_CLUSTER_MASS_PROJECTION_SIGMA_P0,
  NC_CLUSTER_MASS_PROJECTION_SIGMA_P1,
  NC_CLUSTER_MASS_PROJECTION_SIGMA_P2,
  NC_CLUSTER_MASS_PROJECTION_F_PRJ,
  NC_CLUSTER_MASS_PROJECTION_TAU,
  /* < private > */
  NC_CLUSTER_MASS_PROJECTION_SPARAM_LEN, /*< skip >*/
} NcClusterMassProjectionSParams;

#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P0    (3.19)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P1    (2.0 / M_LN10)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P2    (-0.7 / M_LN10)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P0 (0.33)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P1 (-0.08 / M_LN10)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P2 (0.0)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_F_PRJ    (0.15)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_TAU      (0.1)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL (0.0)

#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_RELTOL (1.0e-7)

void nc_cluster_mass_projection_set_reltol (NcClusterMassProjection *mp, gdouble reltol);
gdouble nc_cluster_mass_projection_get_reltol (NcClusterMassProjection *mp);

gdouble nc_cluster_mass_projection_f_prj (NcClusterMassProjection *mp);
gdouble nc_cluster_mass_projection_tau (NcClusterMassProjection *mp);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_PROJECTION_H_ */

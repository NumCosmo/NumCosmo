/***************************************************************************
 * nc_cluster_mass_projection.h
 *
 * Thu Sep 10 18:25:11 2026
 * Copyright 2026 Cinthia Nunes de Lima and Henrique Lettieri Projection
 * <cinthia.nlima@gmail.com>, <henrique.lettieri@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Cinthia Nunes de Lima and Henrique Lettieri Projection 2026 <cinthia.nlima@gmail.com>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_CLUSTER_MASS_PROJECTION_H_
#define _NC_CLUSTER_MASS_PROJECTION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/spline/ncm_spline2d_bicubic.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_MASS_PROJECTION (nc_cluster_mass_projection_get_type ())
G_DECLARE_FINAL_TYPE (NcClusterMassProjection, nc_cluster_mass_projection, NC, CLUSTER_MASS_PROJECTION, NcClusterMass)


/**
 * NcClusterMassProjectionSParams:
 * @NC_CLUSTER_MASS_PROJECTION_MU_P0: bias of the mean
 * @NC_CLUSTER_MASS_PROJECTION_MU_P1: slope on the mean
 * @NC_CLUSTER_MASS_PROJECTION_MU_P2: redshift dependency on the mean
 * @NC_CLUSTER_MASS_PROJECTION_SIGMA_P0: bias of the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_PROJECTION_SIGMA_P1: slope on the standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_PROJECTION_SIGMA_P2: redshift dependency standard deviation of the log-normal distribution
 * @NC_CLUSTER_MASS_PROJECTION_CUT: cut in richness
 * @NC_CLUSTER_MASS_PROJECTION_TAU: exponential-tail scale of the projection contamination
 * @NC_CLUSTER_MASS_PROJECTION_F_PROJ: fraction of projected (contaminated) clusters
 *
 * Cluster mass distribution with projection effects, model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_MASS_PROJECTION_SPARAMS,prefix=NC_CLUSTER_MASS_PROJECTION >*/
{
  NC_CLUSTER_MASS_PROJECTION_MU_P0,
  NC_CLUSTER_MASS_PROJECTION_MU_P1,
  NC_CLUSTER_MASS_PROJECTION_MU_P2,
  NC_CLUSTER_MASS_PROJECTION_SIGMA_P0,
  NC_CLUSTER_MASS_PROJECTION_SIGMA_P1,
  NC_CLUSTER_MASS_PROJECTION_SIGMA_P2,
  NC_CLUSTER_MASS_PROJECTION_CUT,
  NC_CLUSTER_MASS_PROJECTION_TAU,
  NC_CLUSTER_MASS_PROJECTION_F_PROJ,
  /* < private > */
  NC_CLUSTER_MASS_PROJECTION_SPARAM_LEN, /*< skip >*/
} NcClusterMassProjectionSParams;

#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P0    (3.19)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P1    (2.0 / M_LN10)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_MU_P2    (-0.7 / M_LN10)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P0 (0.33)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P1 (-0.08 / M_LN10)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_SIGMA_P2 (0.0)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_CUT      (0.0)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_TAU      (1.0)
#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_F_PROJ   (0.2)

#define NC_CLUSTER_MASS_PROJECTION_DEFAULT_PARAMS_ABSTOL (0.0)

void nc_cluster_mass_projection_set_enable_rejection (NcClusterMassProjection *projection, gboolean on);
void nc_cluster_mass_projection_set_ipurity (NcClusterMassProjection *projection, NcmSpline2dBicubic *ipurity);
void nc_cluster_mass_projection_set_completeness (NcClusterMassProjection *projection, NcmSpline2dBicubic *completeness);
void nc_cluster_mass_projection_set_lnM_limits (NcClusterMassProjection *projection, NcmVector *lnM_limits);

gdouble nc_cluster_mass_projection_completeness (NcClusterMassProjection *projection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_projection_ipurity (NcClusterMassProjection *projection, gdouble lnM_obs, gdouble z);
gdouble nc_cluster_mass_projection_get_mean_richness (NcClusterMassProjection *projection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_projection_get_std_richness (NcClusterMassProjection *projection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_projection_get_cut (NcClusterMassProjection *projection);
gdouble nc_cluster_mass_projection_get_mean (NcClusterMassProjection *projection, gdouble lnM, gdouble z);
gdouble nc_cluster_mass_projection_get_std (NcClusterMassProjection *projection, gdouble lnM, gdouble z);
gboolean nc_cluster_mass_projection_get_enable_rejection (NcClusterMassProjection *projection);

G_END_DECLS

#endif /* _NC_CLUSTER_MASS_PROJECTION_H_ */


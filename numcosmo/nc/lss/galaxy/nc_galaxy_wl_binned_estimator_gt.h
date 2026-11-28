/***************************************************************************
 *            nc_galaxy_wl_binned_estimator_gt.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator_gt.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_WL_BINNED_ESTIMATOR_GT_H_
#define _NC_GALAXY_WL_BINNED_ESTIMATOR_GT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_binned_estimator.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR_GT (nc_galaxy_wl_binned_estimator_gt_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLBinnedEstimatorGt, nc_galaxy_wl_binned_estimator_gt, NC, GALAXY_WL_BINNED_ESTIMATOR_GT, NcGalaxyWLBinnedEstimator)

/**
 * NcGalaxyWLBinnedEstimatorGtError:
 * @NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_STD_ERR: the standard error of the weighted mean
 *   of the bin, from the scatter of the per-galaxy shears: the weighted variance of
 *   the bin times $\sum w_i^2 / (\sum w_i)^2$.
 * @NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_SHAPE_NOISE: the propagated shape noise of the
 *   weighted mean, from the per-galaxy intrinsic ellipticity rms and measurement noise.
 *
 * The error of the binned reduced tangential shear.
 */
typedef enum _NcGalaxyWLBinnedEstimatorGtError /*< enum,prefix=NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR >*/
{
  NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_STD_ERR,
  NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_SHAPE_NOISE,
  /* < private > */
  NC_GALAXY_WL_BINNED_ESTIMATOR_GT_ERROR_LEN, /*< skip >*/
} NcGalaxyWLBinnedEstimatorGtError;

#define NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S_MEAN "beta_s_mean"
#define NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_BETA_S2_MEAN "beta_s2_mean"
#define NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_RESPONSIVITY "responsivity"
#define NC_GALAXY_WL_BINNED_ESTIMATOR_GT_COL_M_MEAN "m_mean"

#define NC_GALAXY_WL_BINNED_ESTIMATOR_GT_DEFAULT_SIGMA_POP (0.26)

NcGalaxyWLBinnedEstimatorGt *nc_galaxy_wl_binned_estimator_gt_new (void);
NcGalaxyWLBinnedEstimatorGt *nc_galaxy_wl_binned_estimator_gt_ref (NcGalaxyWLBinnedEstimatorGt *est_gt);
void nc_galaxy_wl_binned_estimator_gt_free (NcGalaxyWLBinnedEstimatorGt *est_gt);
void nc_galaxy_wl_binned_estimator_gt_clear (NcGalaxyWLBinnedEstimatorGt **est_gt);

void nc_galaxy_wl_binned_estimator_gt_set_sigma_pop (NcGalaxyWLBinnedEstimatorGt *est_gt, const gdouble sigma_pop);
gdouble nc_galaxy_wl_binned_estimator_gt_get_sigma_pop (NcGalaxyWLBinnedEstimatorGt *est_gt);
void nc_galaxy_wl_binned_estimator_gt_set_error_type (NcGalaxyWLBinnedEstimatorGt *est_gt, NcGalaxyWLBinnedEstimatorGtError error_type);
NcGalaxyWLBinnedEstimatorGtError nc_galaxy_wl_binned_estimator_gt_get_error_type (NcGalaxyWLBinnedEstimatorGt *est_gt);
void nc_galaxy_wl_binned_estimator_gt_set_beta_weighted (NcGalaxyWLBinnedEstimatorGt *est_gt, gboolean beta_weighted);
gboolean nc_galaxy_wl_binned_estimator_gt_get_beta_weighted (NcGalaxyWLBinnedEstimatorGt *est_gt);

G_END_DECLS

#endif /* _NC_GALAXY_WL_BINNED_ESTIMATOR_GT_H_ */


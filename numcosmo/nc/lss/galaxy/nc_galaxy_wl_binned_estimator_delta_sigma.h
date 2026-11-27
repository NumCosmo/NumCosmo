/***************************************************************************
 *            nc_galaxy_wl_binned_estimator_delta_sigma.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator_delta_sigma.h
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

#ifndef _NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_H_
#define _NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_binned_estimator.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA (nc_galaxy_wl_binned_estimator_delta_sigma_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLBinnedEstimatorDeltaSigma, nc_galaxy_wl_binned_estimator_delta_sigma, NC, GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA, NcGalaxyWLBinnedEstimator)

#define NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_INV_SIGMA_CRIT "inv_sigma_crit"
#define NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_RESPONSIVITY "responsivity"
#define NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_COL_MULT_BIAS "mult_bias"

NcGalaxyWLBinnedEstimatorDeltaSigma * nc_galaxy_wl_binned_estimator_delta_sigma_new (void);
NcGalaxyWLBinnedEstimatorDeltaSigma *nc_galaxy_wl_binned_estimator_delta_sigma_ref (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds);
void nc_galaxy_wl_binned_estimator_delta_sigma_free (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds);
void nc_galaxy_wl_binned_estimator_delta_sigma_clear (NcGalaxyWLBinnedEstimatorDeltaSigma **est_ds);

void nc_galaxy_wl_binned_estimator_delta_sigma_set_subtract_additive_bias (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds, gboolean subtract);
gboolean nc_galaxy_wl_binned_estimator_delta_sigma_get_subtract_additive_bias (NcGalaxyWLBinnedEstimatorDeltaSigma *est_ds);

G_END_DECLS

#endif /* _NC_GALAXY_WL_BINNED_ESTIMATOR_DELTA_SIGMA_H_ */


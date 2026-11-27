/***************************************************************************
 *            nc_galaxy_wl_binned_estimator.h
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator.h
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

#ifndef _NC_GALAXY_WL_BINNED_ESTIMATOR_H_
#define _NC_GALAXY_WL_BINNED_ESTIMATOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs_binned.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR (nc_galaxy_wl_binned_estimator_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyWLBinnedEstimator, nc_galaxy_wl_binned_estimator, NC, GALAXY_WL_BINNED_ESTIMATOR, GObject)

typedef struct _NcGalaxyWLBinnedEstimatorSources NcGalaxyWLBinnedEstimatorSources;

struct _NcGalaxyWLBinnedEstimatorClass
{
  /*< private >*/
  GObjectClass parent_class;

  const gchar *name;
  GStrv (*extra_columns) (NcGalaxyWLBinnedEstimator *est);
  void (*build_bins) (NcGalaxyWLBinnedEstimator *est, const NcGalaxyWLBinnedEstimatorSources *src, NcGalaxyWLObsBinned *binned, GArray *factor);
  void (*eval_model) (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model);

  /* Padding to allow adding up to 8 more virtual functions without breaking ABI. */
  gpointer padding[8];
};

/* Optional columns of the per-galaxy catalog; absent ones are taken as zero */
#define NC_GALAXY_WL_BINNED_ESTIMATOR_COL_E_RMS "e_rms"

#define NC_GALAXY_WL_BINNED_ESTIMATOR_DEFAULT_PZ_INTERVALS (200)
#define NC_GALAXY_WL_BINNED_ESTIMATOR_DEFAULT_EXACT_INTERVALS (20)

NcGalaxyWLBinnedEstimator *nc_galaxy_wl_binned_estimator_ref (NcGalaxyWLBinnedEstimator *est);
void nc_galaxy_wl_binned_estimator_free (NcGalaxyWLBinnedEstimator *est);
void nc_galaxy_wl_binned_estimator_clear (NcGalaxyWLBinnedEstimator **est);

const gchar *nc_galaxy_wl_binned_estimator_get_name (NcGalaxyWLBinnedEstimator *est);

void nc_galaxy_wl_binned_estimator_set_e_rms (NcGalaxyWLBinnedEstimator *est, const gdouble e_rms);
gdouble nc_galaxy_wl_binned_estimator_get_e_rms (NcGalaxyWLBinnedEstimator *est);
void nc_galaxy_wl_binned_estimator_set_pz_intervals (NcGalaxyWLBinnedEstimator *est, const guint pz_intervals);
guint nc_galaxy_wl_binned_estimator_get_pz_intervals (NcGalaxyWLBinnedEstimator *est);
void nc_galaxy_wl_binned_estimator_set_exact_intervals (NcGalaxyWLBinnedEstimator *est, const guint exact_intervals);
guint nc_galaxy_wl_binned_estimator_get_exact_intervals (NcGalaxyWLBinnedEstimator *est);
void nc_galaxy_wl_binned_estimator_set_redshift_factor (NcGalaxyWLBinnedEstimator *est, NcGalaxyRedshiftFactor *redshift_factor);
NcGalaxyRedshiftFactor *nc_galaxy_wl_binned_estimator_peek_redshift_factor (NcGalaxyWLBinnedEstimator *est);

GList *nc_galaxy_wl_binned_estimator_required_columns (NcGalaxyWLBinnedEstimator *est);

NcGalaxyWLObsBinned *nc_galaxy_wl_binned_estimator_build (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObs *obs, GArray *edges);
void nc_galaxy_wl_binned_estimator_eval_model (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model);
void nc_galaxy_wl_binned_estimator_eval_model_exact (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model);
gdouble nc_galaxy_wl_binned_estimator_mean_inv_sigma_crit (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObs *obs, const guint i, gdouble *mean_inv_sigma_crit2);

G_END_DECLS

#endif /* _NC_GALAXY_WL_BINNED_ESTIMATOR_H_ */


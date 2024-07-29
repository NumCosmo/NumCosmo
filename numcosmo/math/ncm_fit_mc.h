/***************************************************************************
 *            ncm_fit_mc.h
 *
 *  Sat December 01 17:19:03 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_FIT_MC_H_
#define _NCM_FIT_MC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_mset_catalog.h>
#include <numcosmo/math/ncm_timer.h>
#include <numcosmo/math/ncm_memory_pool.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_MC (ncm_fit_mc_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitMC, ncm_fit_mc, NCM, FIT_MC, GObject)

/**
 * NcmFitMCResampleType:
 * @NCM_FIT_MC_RESAMPLE_FROM_MODEL: Montecarlo resampling from models
 * @NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX: Montecarlo bootstraping each #NcmData separately.
 * @NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX: Montecarlo bootstraping mixing all #NcmData in the bootstrap process.
 *
 * Montecarlo resample options
 *
 */
typedef enum _NcmFitMCResampleType
{
  NCM_FIT_MC_RESAMPLE_FROM_MODEL = 0,
  NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX,
  NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX,
  /* < private > */
  NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN, /*< skip >*/
} NcmFitMCResampleType;

typedef void (*NcmFitMCResample) (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng);

NcmFitMC *ncm_fit_mc_new (NcmFit *fit, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype);
void ncm_fit_mc_free (NcmFitMC *mc);
void ncm_fit_mc_clear (NcmFitMC **mc);

void ncm_fit_mc_set_data_file (NcmFitMC *mc, const gchar *filename);

void ncm_fit_mc_set_mtype (NcmFitMC *mc, NcmFitRunMsgs mtype);
void ncm_fit_mc_set_rtype (NcmFitMC *mc, NcmFitMCResampleType rtype);
void ncm_fit_mc_set_nthreads (NcmFitMC *mc, guint nthreads);
void ncm_fit_mc_keep_order (NcmFitMC *mc, gboolean keep_order);
void ncm_fit_mc_set_fiducial (NcmFitMC *mc, NcmMSet *fiduc);
void ncm_fit_mc_set_rng (NcmFitMC *mc, NcmRNG *rng);

gboolean ncm_fit_mc_is_running (NcmFitMC *mc);

void ncm_fit_mc_start_run (NcmFitMC *mc);
void ncm_fit_mc_end_run (NcmFitMC *mc);
void ncm_fit_mc_reset (NcmFitMC *mc);
void ncm_fit_mc_set_first_sample_id (NcmFitMC *mc, gint first_sample_id);
void ncm_fit_mc_run (NcmFitMC *mc, guint n);
void ncm_fit_mc_run_lre (NcmFitMC *mc, guint prerun, gdouble lre);
void ncm_fit_mc_mean_covar (NcmFitMC *mc);

NcmMSetCatalog *ncm_fit_mc_get_catalog (NcmFitMC *mc);
NcmMSetCatalog *ncm_fit_mc_peek_catalog (NcmFitMC *mc);

#define NCM_FIT_MC_MIN_SYNC_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_FIT_MC_H_ */


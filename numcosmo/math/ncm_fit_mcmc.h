/***************************************************************************
 *            ncm_fit_mcmc.h
 *
 *  Sun May 25 16:42:07 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_FIT_MCMC_H_
#define _NCM_FIT_MCMC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_mset_catalog.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>
#include <numcosmo/math/ncm_timer.h>
#include <numcosmo/math/ncm_memory_pool.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_MCMC (ncm_fit_mcmc_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitMCMC, ncm_fit_mcmc, NCM, FIT_MCMC, GObject)

NcmFitMCMC *ncm_fit_mcmc_new (NcmFit * fit, NcmMSetTransKern * tkern, NcmFitRunMsgs mtype);
void ncm_fit_mcmc_free (NcmFitMCMC *mcmc);
void ncm_fit_mcmc_clear (NcmFitMCMC **mcmc);

void ncm_fit_mcmc_set_data_file (NcmFitMCMC *mcmc, const gchar *filename);

void ncm_fit_mcmc_set_mtype (NcmFitMCMC *mcmc, NcmFitRunMsgs mtype);
void ncm_fit_mcmc_set_trans_kern (NcmFitMCMC *mcmc, NcmMSetTransKern *tkern);
void ncm_fit_mcmc_set_nthreads (NcmFitMCMC *mcmc, guint nthreads);
void ncm_fit_mcmc_set_fiducial (NcmFitMCMC *mcmc, NcmMSet *fiduc);
void ncm_fit_mcmc_set_rng (NcmFitMCMC *mcmc, NcmRNG *rng);

gdouble ncm_fit_mcmc_get_accept_ratio (NcmFitMCMC *mcmc);

void ncm_fit_mcmc_start_run (NcmFitMCMC *mcmc);
void ncm_fit_mcmc_end_run (NcmFitMCMC *mcmc);
void ncm_fit_mcmc_reset (NcmFitMCMC *mcmc);
void ncm_fit_mcmc_set_first_sample_id (NcmFitMCMC *mcmc, gint first_sample_id);
void ncm_fit_mcmc_run (NcmFitMCMC *mcmc, guint n);
void ncm_fit_mcmc_run_lre (NcmFitMCMC *mcmc, guint prerun, gdouble lre);
void ncm_fit_mcmc_mean_covar (NcmFitMCMC *mcmc);

NcmMSetCatalog *ncm_fit_mcmc_get_catalog (NcmFitMCMC *mcmc);

#define NCM_FIT_MCMC_MIN_SYNC_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_FIT_MCMC_H_ */


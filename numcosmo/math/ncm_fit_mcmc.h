/***************************************************************************
 *            ncm_fit_mcmc.h
 *
 *  Sun May 25 16:42:07 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
#include <numcosmo/math/memory_pool.h>
#include <gsl/gsl_histogram.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

G_BEGIN_DECLS

#define NCM_TYPE_FIT_MCMC             (ncm_fit_mcmc_get_type ())
#define NCM_FIT_MCMC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_MCMC, NcmFitMCMC))
#define NCM_FIT_MCMC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_MCMC, NcmFitMCMCClass))
#define NCM_IS_FIT_MCMC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_MCMC))
#define NCM_IS_FIT_MCMC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_MCMC))
#define NCM_FIT_MCMC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_MCMC, NcmFitMCMCClass))

typedef struct _NcmFitMCMCClass NcmFitMCMCClass;
typedef struct _NcmFitMCMC NcmFitMCMC;

struct _NcmFitMCMC
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmMSetCatalog *mcat;
  NcmFitRunMsgs mtype;
  NcmTimer *nt;
  NcmSerialize *ser;
  NcmMSetTransKern *tkern;
  NcmVector *theta;
  NcmVector *thetastar;
  guint nthreads;
  guint n;
  NcmMemoryPool *mp;
  gint write_index;
  gint cur_sample_id;
  guint naccepted;
  guint ntotal;
  gboolean started;
  GMutex dup_fit;
  GMutex resample_lock;
  GMutex update_lock;
  GCond write_cond;
};

struct _NcmFitMCMCClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_fit_mcmc_get_type (void) G_GNUC_CONST;

NcmFitMCMC *ncm_fit_mcmc_new (NcmFit *fit, NcmMSetTransKern *tkern, NcmFitRunMsgs mtype);
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

#define NCM_FIT_MCMC_MIN_FLUSH_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_FIT_MCMC_H_ */

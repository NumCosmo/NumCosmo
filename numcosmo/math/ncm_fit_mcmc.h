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
#include <numcosmo/math/ncm_fit_catalog.h>
#include <numcosmo/math/ncm_mc_sampler.h>
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
  NcmFitCatalog *fcat;
  NcmFitRunMsgs mtype;
  NcmRNG *rng;
  NcmTimer *nt;
  NcmSerialize *ser;
  NcmMCSampler *mcs;
  NcmVector *theta;
  NcmVector *thetastar;
  guint nthreads;
  guint n;
  NcmMemoryPool *mp;
  gint write_index;
  gint cur_sample_id;
  gint first_sample_id;
  gboolean started;
  GMutex *dup_fit;
  GMutex *resample_lock;
  GCond *write_cond;
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32))
  GMutex dup_fit_m;
  GMutex resample_lock_m;
  GCond write_cond_m;
#endif
};

struct _NcmFitMCMCClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_fit_mcmc_get_type (void) G_GNUC_CONST;

NcmFitMCMC *ncm_fit_mcmc_new (NcmFit *fit, NcmMCSampler *mcs, NcmFitRunMsgs mtype);
void ncm_fit_mcmc_free (NcmFitMCMC *mc);
void ncm_fit_mcmc_clear (NcmFitMCMC **mc);

void ncm_fit_mcmc_set_data_file (NcmFitMCMC *mc, const gchar *filename);

void ncm_fit_mcmc_set_mtype (NcmFitMCMC *mc, NcmFitRunMsgs mtype);
void ncm_fit_mcmc_set_sampler (NcmFitMCMC *mc, NcmMCSampler *mcs);
void ncm_fit_mcmc_set_nthreads (NcmFitMCMC *mc, guint nthreads);
void ncm_fit_mcmc_set_fiducial (NcmFitMCMC *mc, NcmMSet *fiduc);
void ncm_fit_mcmc_set_rng (NcmFitMCMC *mc, NcmRNG *rng);

void ncm_fit_mcmc_start_run (NcmFitMCMC *mc);
void ncm_fit_mcmc_end_run (NcmFitMCMC *mc);
void ncm_fit_mcmc_reset (NcmFitMCMC *mc);
void ncm_fit_mcmc_set_first_sample_id (NcmFitMCMC *mc, gint first_sample_id);
void ncm_fit_mcmc_run (NcmFitMCMC *mc, guint n);
void ncm_fit_mcmc_run_lre (NcmFitMCMC *mc, guint prerun, gdouble lre);
void ncm_fit_mcmc_mean_covar (NcmFitMCMC *mc);

#define NCM_FIT_MCMC_MIN_FLUSH_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_FIT_MCMC_H_ */

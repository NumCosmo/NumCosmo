/***************************************************************************
 *            ncm_fit_mc.h
 *
 *  Sat December 01 17:19:03 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_FIT_MC_H_
#define _NCM_FIT_MC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_mset_catalog.h>
#include <numcosmo/math/ncm_timer.h>
#include <numcosmo/math/memory_pool.h>
#include <gsl/gsl_histogram.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

G_BEGIN_DECLS

#define NCM_TYPE_FIT_MC             (ncm_fit_mc_get_type ())
#define NCM_FIT_MC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_MC, NcmFitMC))
#define NCM_FIT_MC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_MC, NcmFitMCClass))
#define NCM_IS_FIT_MC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_MC))
#define NCM_IS_FIT_MC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_MC))
#define NCM_FIT_MC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_MC, NcmFitMCClass))

typedef struct _NcmFitMCClass NcmFitMCClass;
typedef struct _NcmFitMC NcmFitMC;

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
  NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX, /*< private >*/
  NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN, /*< skip >*/
} NcmFitMCResampleType;

typedef void (*NcmFitMCResample) (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng);

struct _NcmFitMC
{
  /*< private >*/
  GObject parent_instance;
  NcmFitMCResample resample;
  NcmFit *fit;
  NcmMSet *fiduc;
  NcmMSetCatalog *mcat;
  NcmFitRunMsgs mtype;
  NcmFitMCResampleType rtype;
  NcmVector *bf;
  NcmTimer *nt;
  NcmSerialize *ser;
  guint nthreads;
  guint n;
  gboolean keep_order;
  NcmMemoryPool *mp;
  gint write_index;
  gint cur_sample_id;
  gint first_sample_id;
  gboolean started;
  GMutex dup_fit;
  GMutex resample_lock;
  GMutex update_lock;
  GCond write_cond;
};

struct _NcmFitMCClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_fit_mc_get_type (void) G_GNUC_CONST;

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

void ncm_fit_mc_start_run (NcmFitMC *mc);
void ncm_fit_mc_end_run (NcmFitMC *mc);
void ncm_fit_mc_reset (NcmFitMC *mc);
void ncm_fit_mc_set_first_sample_id (NcmFitMC *mc, gint first_sample_id);
void ncm_fit_mc_run (NcmFitMC *mc, guint n);
void ncm_fit_mc_run_lre (NcmFitMC *mc, guint prerun, gdouble lre);
void ncm_fit_mc_mean_covar (NcmFitMC *mc);

NcmMSetCatalog *ncm_fit_mc_get_catalog (NcmFitMC *mc);

#define NCM_FIT_MC_MIN_FLUSH_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_FIT_MC_H_ */

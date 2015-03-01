/***************************************************************************
 *            ncm_fit_esmcmc.h
 *
 *  Tue January 20 16:59:54 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti & Mariana Penna-Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
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

#ifndef _NCM_FIT_ESMCMC_H_
#define _NCM_FIT_ESMCMC_H_

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

#define NCM_TYPE_FIT_ESMCMC             (ncm_fit_esmcmc_get_type ())
#define NCM_FIT_ESMCMC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_ESMCMC, NcmFitESMCMC))
#define NCM_FIT_ESMCMC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_ESMCMC, NcmFitESMCMCClass))
#define NCM_IS_FIT_ESMCMC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_ESMCMC))
#define NCM_IS_FIT_ESMCMC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_ESMCMC))
#define NCM_FIT_ESMCMC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_ESMCMC, NcmFitESMCMCClass))

typedef struct _NcmFitESMCMCClass NcmFitESMCMCClass;
typedef struct _NcmFitESMCMC NcmFitESMCMC;

/**
 * NcmFitESMCMCMoveType:
 * @NCM_FIT_ESMCMC_MOVE_TYPE_STRETCH: Stretch move.
 * 
 * Ensemble Sampler Move Type. 
 * 
 */
typedef enum _NcmFitESMCMCMoveType
{
  NCM_FIT_ESMCMC_MOVE_TYPE_STRETCH,  /*< private >*/
  NCM_FIT_ESMCMC_MOVE_TYPE_LEN,      /*< skip >*/
} NcmFitESMCMCMoveType;

struct _NcmFitESMCMC
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  GPtrArray *walker_fits;
  NcmMSetTransKern *sampler;
  NcmMSetCatalog *mcat;
  NcmFitRunMsgs mtype;
  NcmTimer *nt;
  NcmSerialize *ser;
  NcmFitESMCMCMoveType mt;
  GPtrArray *theta_k;
  GPtrArray *theta_j;
  GPtrArray *thetastar;
  gdouble a;
  guint fparam_len;
  guint nthreads;
  guint n;
  gint nwalkers;
  gint write_index;
  gint cur_sample_id;
  guint ntotal;
  guint naccepted;
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

struct _NcmFitESMCMCClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_fit_esmcmc_get_type (void) G_GNUC_CONST;

NcmFitESMCMC *ncm_fit_esmcmc_new (NcmFit *fit, gint nwalkers, NcmMSetTransKern *sampler, NcmFitESMCMCMoveType mt, NcmFitRunMsgs mtype);
void ncm_fit_esmcmc_free (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_clear (NcmFitESMCMC **esmcmc);

void ncm_fit_esmcmc_set_data_file (NcmFitESMCMC *esmcmc, const gchar *filename);

void ncm_fit_esmcmc_set_sampler (NcmFitESMCMC *esmcmc, NcmMSetTransKern *sampler);
void ncm_fit_esmcmc_set_mtype (NcmFitESMCMC *esmcmc, NcmFitRunMsgs mtype);
void ncm_fit_esmcmc_set_move_type (NcmFitESMCMC *esmcmc, NcmFitESMCMCMoveType mt);
void ncm_fit_esmcmc_set_move_scale (NcmFitESMCMC *esmcmc, gdouble a);
void ncm_fit_esmcmc_set_nthreads (NcmFitESMCMC *esmcmc, guint nthreads);
void ncm_fit_esmcmc_set_rng (NcmFitESMCMC *esmcmc, NcmRNG *rng);

gdouble ncm_fit_esmcmc_get_accept_ratio (NcmFitESMCMC *esmcmc);

void ncm_fit_esmcmc_start_run (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_end_run (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_reset (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_run (NcmFitESMCMC *esmcmc, guint n);
void ncm_fit_esmcmc_run_lre (NcmFitESMCMC *esmcmc, guint prerun, gdouble lre);
void ncm_fit_esmcmc_mean_covar (NcmFitESMCMC *esmcmc);

NcmMSetCatalog *ncm_fit_esmcmc_get_catalog (NcmFitESMCMC *esmcmc);

#define NCM_FIT_ESMCMC_WALKER_ID "NcmFitESMCMC:Walker"
#define NCM_FIT_ESMCMC_MIN_FLUSH_INTERVAL (10.0)

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_H_ */

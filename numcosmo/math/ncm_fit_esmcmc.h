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
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_histogram.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC             (ncm_fit_esmcmc_get_type ())
#define NCM_FIT_ESMCMC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_ESMCMC, NcmFitESMCMC))
#define NCM_FIT_ESMCMC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_ESMCMC, NcmFitESMCMCClass))
#define NCM_IS_FIT_ESMCMC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_ESMCMC))
#define NCM_IS_FIT_ESMCMC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_ESMCMC))
#define NCM_FIT_ESMCMC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_ESMCMC, NcmFitESMCMCClass))

typedef struct _NcmFitESMCMCClass NcmFitESMCMCClass;
typedef struct _NcmFitESMCMC NcmFitESMCMC;
typedef struct _NcmFitESMCMCPrivate NcmFitESMCMCPrivate;

struct _NcmFitESMCMCClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmFitESMCMC
{
  /*< private >*/
  GObject parent_instance;
	NcmFitESMCMCPrivate *priv;
};

GType ncm_fit_esmcmc_get_type (void) G_GNUC_CONST;

NcmFitESMCMC *ncm_fit_esmcmc_new (NcmFit *fit, gint nwalkers, NcmMSetTransKern *sampler, NcmFitESMCMCWalker *walker, NcmFitRunMsgs mtype);
NcmFitESMCMC *ncm_fit_esmcmc_new_funcs_array (NcmFit *fit, gint nwalkers, NcmMSetTransKern *sampler, NcmFitESMCMCWalker *walker, NcmFitRunMsgs mtype, NcmObjArray *funcs_array);

NcmFitESMCMC *ncm_fit_esmcmc_ref (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_free (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_clear (NcmFitESMCMC **esmcmc);

void ncm_fit_esmcmc_set_data_file (NcmFitESMCMC *esmcmc, const gchar *filename);

void ncm_fit_esmcmc_set_sampler (NcmFitESMCMC *esmcmc, NcmMSetTransKern *sampler);
void ncm_fit_esmcmc_set_mtype (NcmFitESMCMC *esmcmc, NcmFitRunMsgs mtype);
void ncm_fit_esmcmc_set_nthreads (NcmFitESMCMC *esmcmc, guint nthreads);
void ncm_fit_esmcmc_use_mpi (NcmFitESMCMC *esmcmc, gboolean use_mpi);
void ncm_fit_esmcmc_set_rng (NcmFitESMCMC *esmcmc, NcmRNG *rng);
void ncm_fit_esmcmc_set_auto_trim (NcmFitESMCMC *esmcmc, gboolean enable);
void ncm_fit_esmcmc_set_auto_trim_div (NcmFitESMCMC *esmcmc, guint div);
void ncm_fit_esmcmc_set_min_runs (NcmFitESMCMC *esmcmc, guint min_runs);
void ncm_fit_esmcmc_set_max_runs_time (NcmFitESMCMC *esmcmc, gdouble max_runs_time);

gboolean ncm_fit_esmcmc_has_rng (NcmFitESMCMC *esmcmc);

gdouble ncm_fit_esmcmc_get_accept_ratio (NcmFitESMCMC *esmcmc);
gdouble ncm_fit_esmcmc_get_offboard_ratio (NcmFitESMCMC *esmcmc);
gdouble ncm_fit_esmcmc_get_accept_ratio_last_update (NcmFitESMCMC *esmcmc);
gdouble ncm_fit_esmcmc_get_offboard_ratio_last_update (NcmFitESMCMC *esmcmc);

void ncm_fit_esmcmc_start_run (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_end_run (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_reset (NcmFitESMCMC *esmcmc);
void ncm_fit_esmcmc_run (NcmFitESMCMC *esmcmc, guint n);
void ncm_fit_esmcmc_run_lre (NcmFitESMCMC *esmcmc, guint prerun, gdouble lre);
void ncm_fit_esmcmc_mean_covar (NcmFitESMCMC *esmcmc);

NcmSerialize *ncm_fit_esmcmc_peek_ser (NcmFitESMCMC *esmcmc);
NcmMSetCatalog *ncm_fit_esmcmc_get_catalog (NcmFitESMCMC *esmcmc);
NcmMSetCatalog *ncm_fit_esmcmc_peek_catalog (NcmFitESMCMC *esmcmc);

gboolean ncm_fit_esmcmc_validate (NcmFitESMCMC *esmcmc, gulong pi, gulong pf);

#define NCM_FIT_ESMCMC_MIN_SYNC_INTERVAL (10.0)
#define NCM_FIT_ESMCMC_M2LNL_ID (0)
#define NCM_FIT_ESMCMC_MPI_IN_LEN (3)
#define NCM_FIT_ESMCMC_MPI_OUT_LEN (1)

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_H_ */

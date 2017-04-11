/***************************************************************************
 *            ncm_mset_catalog.h
 *
 *  Tue February 18 10:49:59 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_MSET_CATALOG_H_
#define _NCM_MSET_CATALOG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_stats_vec.h>
#include <numcosmo/math/ncm_stats_dist1d_epdf.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector_complex.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

G_BEGIN_DECLS

#define NCM_TYPE_MSET_CATALOG             (ncm_mset_catalog_get_type ())
#define NCM_MSET_CATALOG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET_CATALOG, NcmMSetCatalog))
#define NCM_MSET_CATALOG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET_CATALOG, NcmMSetCatalogClass))
#define NCM_IS_MSET_CATALOG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET_CATALOG))
#define NCM_IS_MSET_CATALOG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET_CATALOG))
#define NCM_MSET_CATALOG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET_CATALOG, NcmMSetCatalogClass))

typedef struct _NcmMSetCatalogClass NcmMSetCatalogClass;
typedef struct _NcmMSetCatalog NcmMSetCatalog;

struct _NcmMSetCatalogClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmMSetCatalogSync:
 * @NCM_MSET_CATALOG_SYNC_DISABLE: Catalog will be synchronized only when closing the file or with an explicit call of ncm_mset_catalog_sync().
 * @NCM_MSET_CATALOG_SYNC_AUTO: Catalog will be synchronized in every catalog addition.
 * @NCM_MSET_CATALOG_SYNC_TIMED: Catalog will be synchronized with a minimum time interval between syncs.
 * 
 * Catalog sync modes. 
 * 
 */
typedef enum _NcmMSetCatalogSync
{
  NCM_MSET_CATALOG_SYNC_DISABLE,
  NCM_MSET_CATALOG_SYNC_AUTO,
  NCM_MSET_CATALOG_SYNC_TIMED, /*< private >*/
  NCM_MSET_CATALOG_SYNC_LEN,   /*< skip >*/
} NcmMSetCatalogSync;

/**
 * NcmMSetCatalogTrimType:
 * @NCM_MSET_CATALOG_TRIM_TYPE_ESS: trim the catalog using the maximum ess criterium.
 * @NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL: trim the catalog using the Heidelberger and Welch’s convergence diagnostic.
 * @NCM_MSET_CATALOG_TRIM_TYPE_ALL: trim the catalog using all tests above.
 * 
 * See ncm_mset_catalog_calc_max_ess_time() and ncm_mset_catalog_calc_heidel_diag().
 * 
 */
typedef enum _NcmMSetCatalogTrimType
{
  NCM_MSET_CATALOG_TRIM_TYPE_ESS    = 1 << 0,
  NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL = 1 << 1,
  NCM_MSET_CATALOG_TRIM_TYPE_ALL    = (1 << 2) - 1,
} NcmMSetCatalogTrimType;

/**
 * NcmMSetCatalogTauMethod:
 * @NCM_MSET_CATALOG_TAU_METHOD_ACOR: uses the autocorrelation to estimate $\tau$.
 * @NCM_MSET_CATALOG_TAU_METHOD_AR_MODEL: uses an autoregressive model fitting to estimate $\tau$.
 * 
 * Method used to estimate the autocorrelation time $\tau$.
 * 
 */
typedef enum _NcmMSetCatalogTauMethod
{
  NCM_MSET_CATALOG_TAU_METHOD_ACOR = 0,
  NCM_MSET_CATALOG_TAU_METHOD_AR_MODEL, /*< private >*/
  NCM_MSET_CATALOG_TAU_METHOD_LEN,      /*< skip >*/
} NcmMSetCatalogTauMethod;

struct _NcmMSetCatalog
{
  /*< private >*/
  GObject parent_instance;
  NcmMSet *mset;
  guint nadd_vals;
  GPtrArray *add_vals_names;
  GPtrArray *add_vals_symbs;
  NcmStatsVec *pstats;
  NcmMSetCatalogSync smode;
  gboolean readonly;
  NcmRNG *rng;
  gboolean weighted;
  gboolean first_flush;
  guint nchains;
  GPtrArray *chain_pstats;
  NcmStatsVec *mean_pstats;
  NcmStatsVec *e_stats;
  NcmStatsVec *e_mean_stats;
  GPtrArray *e_var_array;
  NcmVector *chain_means;
  NcmVector *chain_vars;
  NcmMatrix *chain_cov;
  NcmMatrix *chain_sM;
  gsl_eigen_nonsymm_workspace *chain_sM_ws;
  gsl_vector_complex *chain_sM_ev;
  NcmMSetCatalogTauMethod tau_method;
  NcmVector *tau;
  gchar *rng_inis;
  gchar *rng_stat;
  GTimer *sync_timer;
  gdouble sync_interval;
  gchar *file;
  gchar *mset_file;
  gchar *rtype_str;
  GArray *porder;
  NcmVector *quantile_ws;
  gint first_id;
  gint cur_id;
  gint file_first_id;
  gint file_cur_id;
  glong burnin;
#ifdef NUMCOSMO_HAVE_CFITSIO
  fitsfile *fptr;
#endif /* NUMCOSMO_HAVE_CFITSIO */
  NcmVector *params_max;
  NcmVector *params_min;
  glong pdf_i;
  gsl_histogram *h;
  gsl_histogram_pdf *h_pdf;
  gboolean constructed;
};

GType ncm_mset_catalog_get_type (void) G_GNUC_CONST;

NcmMSetCatalog *ncm_mset_catalog_new (NcmMSet *mset, guint nadd_vals, guint nchains, gboolean weighted, ...) G_GNUC_NULL_TERMINATED;
NcmMSetCatalog *ncm_mset_catalog_new_array (NcmMSet *mset, guint nadd_vals, guint nchains, gboolean weighted, gchar **names, gchar **symbols);

NcmMSetCatalog *ncm_mset_catalog_new_from_file (const gchar *filename, glong burnin);
NcmMSetCatalog *ncm_mset_catalog_new_from_file_ro (const gchar *filename, glong burnin);
NcmMSetCatalog *ncm_mset_catalog_ref (NcmMSetCatalog *mcat);
void ncm_mset_catalog_free (NcmMSetCatalog *mcat);
void ncm_mset_catalog_clear (NcmMSetCatalog **mcat);

void ncm_mset_catalog_set_file (NcmMSetCatalog *mcat, const gchar *filename);
void ncm_mset_catalog_set_sync_mode (NcmMSetCatalog *mcat, NcmMSetCatalogSync smode);
void ncm_mset_catalog_set_sync_interval (NcmMSetCatalog *mcat, gdouble interval);
void ncm_mset_catalog_set_first_id (NcmMSetCatalog *mcat, gint first_id);
void ncm_mset_catalog_set_run_type (NcmMSetCatalog *mcat, const gchar *rtype_str);
void ncm_mset_catalog_set_rng (NcmMSetCatalog *mcat, NcmRNG *rng);
void ncm_mset_catalog_sync (NcmMSetCatalog *mcat, gboolean check);
void ncm_mset_catalog_timed_sync (NcmMSetCatalog *mcat, gboolean check);
void ncm_mset_catalog_reset_stats (NcmMSetCatalog *mcat);
void ncm_mset_catalog_reset (NcmMSetCatalog *mcat);
void ncm_mset_catalog_erase_data (NcmMSetCatalog *mcat);

const gchar *ncm_mset_catalog_peek_filename (NcmMSetCatalog *mcat);
NcmRNG *ncm_mset_catalog_get_rng (NcmMSetCatalog *mcat);
NcmRNG *ncm_mset_catalog_peek_rng (NcmMSetCatalog *mcat);

gboolean ncm_mset_catalog_is_empty (NcmMSetCatalog *mcat);
gdouble ncm_mset_catalog_largest_error (NcmMSetCatalog *mcat);
guint ncm_mset_catalog_len (NcmMSetCatalog *mcat);
guint ncm_mset_catalog_max_time (NcmMSetCatalog *mcat);
guint ncm_mset_catalog_nchains (NcmMSetCatalog *mcat);
guint ncm_mset_catalog_get_row_from_time (NcmMSetCatalog *mcat, gint t);
gint ncm_mset_catalog_get_first_id (NcmMSetCatalog *mcat);
gint ncm_mset_catalog_get_cur_id (NcmMSetCatalog *mcat);

void ncm_mset_catalog_set_burnin (NcmMSetCatalog *mcat, glong burnin);
glong ncm_mset_catalog_get_burnin (NcmMSetCatalog *mcat);

void ncm_mset_catalog_set_tau_method (NcmMSetCatalog *mcat, NcmMSetCatalogTauMethod tau_method);
NcmMSetCatalogTauMethod ncm_mset_catalog_get_tau_method (NcmMSetCatalog *mcat);

void ncm_mset_catalog_add_from_mset (NcmMSetCatalog *mcat, NcmMSet *mset, ...) G_GNUC_NULL_TERMINATED;
void ncm_mset_catalog_add_from_mset_array (NcmMSetCatalog *mcat, NcmMSet *mset, gdouble *ax);
void ncm_mset_catalog_add_from_vector (NcmMSetCatalog *mcat, NcmVector *vals);
void ncm_mset_catalog_log_current_stats (NcmMSetCatalog *mcat);
void ncm_mset_catalog_log_current_chain_stats (NcmMSetCatalog *mcat);

NcmMSet *ncm_mset_catalog_get_mset (NcmMSetCatalog *mcat);
const gchar *ncm_mset_catalog_get_run_type (NcmMSetCatalog *mcat);

NcmVector *ncm_mset_catalog_peek_row (NcmMSetCatalog *mcat, guint i);
NcmVector *ncm_mset_catalog_peek_current_row (NcmMSetCatalog *mcat);
NcmVector *ncm_mset_catalog_peek_current_e_mean (NcmMSetCatalog *mcat);
NcmVector *ncm_mset_catalog_peek_current_e_var (NcmMSetCatalog *mcat);
NcmVector *ncm_mset_catalog_peek_e_mean_t (NcmMSetCatalog *mcat, guint t);
NcmVector *ncm_mset_catalog_peek_e_var_t (NcmMSetCatalog *mcat, guint t);

void ncm_mset_catalog_get_mean (NcmMSetCatalog *mcat, NcmVector  **mean);
void ncm_mset_catalog_get_covar (NcmMSetCatalog *mcat, NcmMatrix **cov);
void ncm_mset_catalog_get_full_covar (NcmMSetCatalog *mcat, NcmMatrix **cov);
void ncm_mset_catalog_log_full_covar (NcmMSetCatalog *mcat);

void ncm_mset_catalog_estimate_autocorrelation_tau (NcmMSetCatalog *mcat, gboolean force_single_chain);
NcmVector *ncm_mset_catalog_peek_autocorrelation_tau (NcmMSetCatalog *mcat);
gdouble ncm_mset_catalog_get_param_shrink_factor (NcmMSetCatalog *mcat, guint p);
gdouble ncm_mset_catalog_get_shrink_factor (NcmMSetCatalog *mcat);

void ncm_mset_catalog_param_pdf (NcmMSetCatalog *mcat, guint i);
gdouble ncm_mset_catalog_param_pdf_pvalue (NcmMSetCatalog *mcat, gdouble pvalue, gboolean both);

NcmMatrix *ncm_mset_catalog_calc_ci_direct (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmVector *x_v, GArray *p_val);
NcmMatrix *ncm_mset_catalog_calc_ci_interp (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmVector *x_v, GArray *p_val, guint nodes, NcmFitRunMsgs mtype);
NcmMatrix *ncm_mset_catalog_calc_pvalue (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmVector *x_v, GArray *lim, guint nodes, NcmFitRunMsgs mtype);
NcmStatsDist1d *ncm_mset_catalog_calc_distrib (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmFitRunMsgs mtype);
NcmStatsDist1d *ncm_mset_catalog_calc_param_distrib (NcmMSetCatalog *mcat, const NcmMSetPIndex *pi, NcmFitRunMsgs mtype);
NcmStatsDist1d *ncm_mset_catalog_calc_add_param_distrib (NcmMSetCatalog *mcat, guint add_param, NcmFitRunMsgs mtype);

void ncm_mset_catalog_calc_param_ensemble_evol (NcmMSetCatalog *mcat, const NcmMSetPIndex *pi, guint nsteps, NcmFitRunMsgs mtype, NcmVector **pval, NcmMatrix **t_evol);
void ncm_mset_catalog_calc_add_param_ensemble_evol (NcmMSetCatalog *mcat, guint add_param, guint nsteps, NcmFitRunMsgs mtype, NcmVector **pval, NcmMatrix **t_evol);

void ncm_mset_catalog_trim (NcmMSetCatalog *mcat, const guint tc);

guint ncm_mset_catalog_calc_max_ess_time (NcmMSetCatalog *mcat, const guint ntests, gdouble *max_ess, NcmFitRunMsgs mtype);
guint ncm_mset_catalog_calc_heidel_diag (NcmMSetCatalog *mcat, const guint ntests, const gdouble pvalue, NcmFitRunMsgs mtype);

void ncm_mset_catalog_trim_by_type (NcmMSetCatalog *mcat, guint ntests, NcmMSetCatalogTrimType trim_type, NcmFitRunMsgs mtype);

guint ncm_mset_catalog_max_ess_time_by_chain (NcmMSetCatalog *mcat, const guint ntests, gdouble *max_ess, NcmFitRunMsgs mtype);
guint ncm_mset_catalog_heidel_diag_by_chain (NcmMSetCatalog *mcat, const guint ntests, const gdouble pvalue, gdouble *wp_pvalue, NcmFitRunMsgs mtype);

#define NCM_MSET_CATALOG_EXTNAME "NcmMSetCatalog:DATA"
#define NCM_MSET_CATALOG_M2LNL_COLNAME "NcmFit:m2lnL"
#define NCM_MSET_CATALOG_M2LNL_SYMBOL "-2\\ln(L)"
#define NCM_MSET_CATALOG_FIRST_ID_LABEL "FIRST_ID"
#define NCM_MSET_CATALOG_RNG_ALGO_LABEL "RNG_ALGO"
#define NCM_MSET_CATALOG_RNG_SEED_LABEL "RNG_SEED"
#define NCM_MSET_CATALOG_RNG_STAT_LABEL "RNG_STAT"
#define NCM_MSET_CATALOG_RNG_INIS_LABEL "RNG_INIS"
#define NCM_MSET_CATALOG_NROWS_LABEL "NAXIS2"
#define NCM_MSET_CATALOG_RTYPE_LABEL "RTYPE"
#define NCM_MSET_CATALOG_NCHAINS_LABEL "NCHAINS"
#define NCM_MSET_CATALOG_NADDVAL_LABEL "NADDVAL"
#define NCM_MSET_CATALOG_WEIGHTED_LABEL "WEIGHTED"
#define NCM_MSET_CATALOG_RTYPE_BSTRAP_MEAN "bootstrap-mean"
#define NCM_MSET_CATALOG_RTYPE_UNDEFINED "undefined-run"
#define NCM_MSET_CATALOG_FSYMB_LABEL "FSYMB"
#define NCM_MSET_CATALOG_ASYMB_LABEL "ASYMB"
#define NCM_MSET_CATALOG_DIST_EST_SD_SCALE (1.0e-3)

G_END_DECLS

#endif /* _NCM_MSET_CATALOG_H_ */


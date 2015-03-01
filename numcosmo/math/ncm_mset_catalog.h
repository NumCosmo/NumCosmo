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

/**
 * NcmMSetCatalogFlush:
 * @NCM_MSET_CATALOG_FLUSH_DISABLE: Catalog will flush only when closing the file.
 * @NCM_MSET_CATALOG_FLUSH_AUTO: Catalog will flush in every catalog addition.
 * @NCM_MSET_CATALOG_FLUSH_TIMED: Catalog will flush with a minimum time interval between flushs.
 * 
 * Catalog flush modes. 
 * 
 */
typedef enum _NcmMSetCatalogFlush
{
  NCM_MSET_CATALOG_FLUSH_DISABLE,
  NCM_MSET_CATALOG_FLUSH_AUTO,
  NCM_MSET_CATALOG_FLUSH_TIMED, /*< private >*/
  NCM_MSET_CATALOG_FLUSH_LEN,   /*< skip >*/
} NcmMSetCatalogFlush;

struct _NcmMSetCatalog
{
  /*< private >*/
  GObject parent_instance;
  NcmMSet *mset;
  guint nadd_vals;
  GPtrArray *add_vals_names;
  NcmStatsVec *pstats;
  NcmMSetCatalogFlush fmode;
  NcmRNG *rng;
  gboolean weighted;
  gboolean first_flush;
  guint nchains;
  GPtrArray *chain_pstats;
  NcmStatsVec *mean_pstats;
  NcmVector *chain_means;
  NcmVector *chain_vars;
  NcmMatrix *chain_cov;
  NcmMatrix *chain_sM;
  gsl_eigen_nonsymm_workspace *chain_sM_ws;
  gsl_vector_complex *chain_sM_ev;
  NcmVector *tau;
  gchar *rng_inis;
  gchar *rng_stat;
  GTimer *flush_timer;
  gdouble flush_interval;
  gchar *file;
  gchar *rtype_str;
  GArray *porder;
  GArray *quantile_ws;
  gint first_id;
  gint cur_id;
  gint file_first_id;
  gint file_cur_id;
#ifdef NUMCOSMO_HAVE_CFITSIO
  fitsfile *fptr;
#endif /* NUMCOSMO_HAVE_CFITSIO */
  NcmVector *params_max;
  NcmVector *params_min;
  glong pdf_i;
  gsl_histogram *h;
  gsl_histogram_pdf *h_pdf;
};

struct _NcmMSetCatalogClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_mset_catalog_get_type (void) G_GNUC_CONST;

NcmMSetCatalog *ncm_mset_catalog_new (NcmMSet *mset, guint nadd_vals, guint nchains, gboolean weighted, ...);
NcmMSetCatalog *ncm_mset_catalog_ref (NcmMSetCatalog *mcat);
void ncm_mset_catalog_free (NcmMSetCatalog *mcat);
void ncm_mset_catalog_clear (NcmMSetCatalog **mcat);

void ncm_mset_catalog_set_add_val_name (NcmMSetCatalog *mcat, guint i, const gchar *name);
void ncm_mset_catalog_set_file (NcmMSetCatalog *mcat, const gchar *filename);
void ncm_mset_catalog_set_flush_mode (NcmMSetCatalog *mcat, NcmMSetCatalogFlush fmode);
void ncm_mset_catalog_set_flush_interval (NcmMSetCatalog *mcat, gdouble interval);
void ncm_mset_catalog_set_first_id (NcmMSetCatalog *mcat, gint first_id);
void ncm_mset_catalog_set_run_type (NcmMSetCatalog *mcat, const gchar *rtype_str);
void ncm_mset_catalog_set_rng (NcmMSetCatalog *mcat, NcmRNG *rng);
void ncm_mset_catalog_sync (NcmMSetCatalog *mcat, gboolean check);
void ncm_mset_catalog_reset (NcmMSetCatalog *mcat);
void ncm_mset_catalog_erase_data (NcmMSetCatalog *mcat);
const gchar *ncm_mset_catalog_peek_filename (NcmMSetCatalog *mcat);

gboolean ncm_mset_catalog_is_empty (NcmMSetCatalog *mcat);
gdouble ncm_mset_catalog_largest_error (NcmMSetCatalog *mcat);
guint ncm_mset_catalog_len (NcmMSetCatalog *mcat);

void ncm_mset_catalog_add_from_mset (NcmMSetCatalog *mcat, NcmMSet *mset, ...);
void ncm_mset_catalog_add_from_mset_array (NcmMSetCatalog *mcat, NcmMSet *mset, gdouble *ax);
void ncm_mset_catalog_add_from_vector (NcmMSetCatalog *mcat, NcmVector *vals);
void ncm_mset_catalog_log_current_stats (NcmMSetCatalog *mcat);
void ncm_mset_catalog_log_current_chain_stats (NcmMSetCatalog *mcat);

NcmVector *ncm_mset_catalog_peek_row (NcmMSetCatalog *mcat, guint i);
NcmVector *ncm_mset_catalog_peek_current_row (NcmMSetCatalog *mcat);

void ncm_mset_catalog_get_mean (NcmMSetCatalog *mcat, NcmVector  **mean);
void ncm_mset_catalog_get_covar (NcmMSetCatalog *mcat, NcmMatrix **cov);

void ncm_mset_catalog_estimate_autocorrelation_tau (NcmMSetCatalog *mcat);
NcmVector *ncm_mset_catalog_peek_autocorrelation_tau (NcmMSetCatalog *mcat);
gdouble ncm_mset_catalog_get_param_shrink_factor (NcmMSetCatalog *mcat, guint p);
gdouble ncm_mset_catalog_get_shrink_factor (NcmMSetCatalog *mcat);

void ncm_mset_catalog_param_pdf (NcmMSetCatalog *mcat, guint i);
gdouble ncm_mset_catalog_param_pdf_pvalue (NcmMSetCatalog *mcat, gdouble pval, gboolean both);

NcmMatrix *ncm_mset_catalog_calc_ci (NcmMSetCatalog *mcat, NcmMSetFunc *func, gdouble *x, GArray *p_val);

#define NCM_MSET_CATALOG_EXTNAME "NcmMSetCatalog:DATA"
#define NCM_MSET_CATALOG_M2LNL_COLNAME "NcmFit:m2lnL"
#define NCM_MSET_CATALOG_FIRST_ID_LABEL "FIRST_ID"
#define NCM_MSET_CATALOG_RNG_ALGO_LABEL "RNG_ALGO"
#define NCM_MSET_CATALOG_RNG_SEED_LABEL "RNG_SEED"
#define NCM_MSET_CATALOG_RNG_STAT_LABEL "RNG_STAT"
#define NCM_MSET_CATALOG_RNG_INIS_LABEL "RNG_INIS"
#define NCM_MSET_CATALOG_NROWS_LABEL "NAXIS2"
#define NCM_MSET_CATALOG_RTYPE_LABEL "RTYPE"
#define NCM_MSET_CATALOG_RTYPE_BSTRAP_MEAN "bootstrap-mean"
#define NCM_MSET_CATALOG_RTYPE_UNDEFINED "undefined-run"
#define NCM_MSET_CATALOG_FSYMB_LABEL "FSYMB"

G_END_DECLS

#endif /* _NCM_MSET_CATALOG_H_ */


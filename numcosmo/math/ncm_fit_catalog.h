/***************************************************************************
 *            ncm_fit_catalog.h
 *
 *  Tue February 18 10:49:59 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_catalog.h
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

#ifndef _NCM_FIT_CATALOG_H_
#define _NCM_FIT_CATALOG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_stats_vec.h>
#include <gsl/gsl_histogram.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

G_BEGIN_DECLS

#define NCM_TYPE_FIT_CATALOG             (ncm_fit_catalog_get_type ())
#define NCM_FIT_CATALOG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_CATALOG, NcmFitCatalog))
#define NCM_FIT_CATALOG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_CATALOG, NcmFitCatalogClass))
#define NCM_IS_FIT_CATALOG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_CATALOG))
#define NCM_IS_FIT_CATALOG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_CATALOG))
#define NCM_FIT_CATALOG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_CATALOG, NcmFitCatalogClass))

typedef struct _NcmFitCatalogClass NcmFitCatalogClass;
typedef struct _NcmFitCatalog NcmFitCatalog;

/**
 * NcmFitCatalogFlush:
 * @NCM_FIT_CATALOG_FLUSH_DISABLE: Catalog will flush only when closing the file.
 * @NCM_FIT_CATALOG_FLUSH_AUTO: Catalog will flush in every catalog addition.
 * @NCM_FIT_CATALOG_FLUSH_TIMED: Catalog will flush with a minimum time interval between flushs.
 * 
 * Catalog flush modes. 
 * 
 */
typedef enum _NcmFitCatalogFlush
{
  NCM_FIT_CATALOG_FLUSH_DISABLE,
  NCM_FIT_CATALOG_FLUSH_AUTO,
  NCM_FIT_CATALOG_FLUSH_TIMED, /*< private >*/
  NCM_FIT_CATALOG_FLUSH_LEN,   /*< skip >*/
} NcmFitCatalogFlush;

struct _NcmFitCatalog
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmStatsVec *pstats;
  NcmFitCatalogFlush fmode;
  GTimer *flush_timer;
  gdouble flush_interval;
  gchar *file;
  gchar *rtype_str;
  GArray *porder;
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

struct _NcmFitCatalogClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_fit_catalog_get_type (void) G_GNUC_CONST;

NcmFitCatalog *ncm_fit_catalog_new (NcmFit *fit);
void ncm_fit_catalog_free (NcmFitCatalog *fcat);
void ncm_fit_catalog_clear (NcmFitCatalog **fcat);

void ncm_fit_catalog_set_file (NcmFitCatalog *fcat, const gchar *filename);
void ncm_fit_catalog_set_flush_mode (NcmFitCatalog *fcat, NcmFitCatalogFlush fmode);
void ncm_fit_catalog_set_flush_interval (NcmFitCatalog *fcat, gdouble interval);
void ncm_fit_catalog_set_first_id (NcmFitCatalog *fcat, gint first_id);
void ncm_fit_catalog_set_run_type (NcmFitCatalog *fcat, const gchar *rtype_str);
void ncm_fit_catalog_sync (NcmFitCatalog *fcat, gboolean check);
void ncm_fit_catalog_reset (NcmFitCatalog *fcat);
void ncm_fit_catalog_erase_data (NcmFitCatalog *fcat);
const gchar *ncm_fit_catalog_peek_filename (NcmFitCatalog *fcat);

gboolean ncm_fit_catalog_is_empty (NcmFitCatalog *fcat);
gdouble ncm_fit_catalog_largest_error (NcmFitCatalog *fcat);
guint ncm_fit_catalog_len (NcmFitCatalog *fcat);
gboolean ncm_fit_catalog_get_prng (NcmFitCatalog *fcat, gchar **prng_algo, gulong *seed);

void ncm_fit_catalog_set_prng (NcmFitCatalog *fcat, NcmRNG *rng);
void ncm_fit_catalog_add_from_fit (NcmFitCatalog *fcat, NcmFit *fit);
void ncm_fit_catalog_add_from_vector (NcmFitCatalog *fcat, NcmVector *vals);
void ncm_fit_catalog_log_current_stats (NcmFitCatalog *fcat);

void ncm_fit_catalog_set_fit_mean_covar (NcmFitCatalog *fcat);
void ncm_fit_catalog_param_pdf (NcmFitCatalog *fcat, guint i);
gdouble ncm_fit_catalog_param_pdf_pvalue (NcmFitCatalog *fcat, gdouble pval, gboolean both);

#define NCM_FIT_CATALOG_EXTNAME "NcmFitCatalog:DATA"
#define NCM_FIT_CATALOG_M2LNL_COLNAME "NcmFit:m2lnL"
#define NCM_FIT_CATALOG_FIRST_ID_LABEL "FIRST_ID"
#define NCM_FIT_CATALOG_PRNG_ALGO_LABEL "RNG_ALGO"
#define NCM_FIT_CATALOG_PRNG_SEED_LABEL "RNG_SEED"
#define NCM_FIT_CATALOG_NROWS_LABEL "NAXIS2"
#define NCM_FIT_CATALOG_RTYPE_LABEL "RTYPE"
#define NCM_FIT_CATALOG_RTYPE_BSTRAP_MEAN "bootstrap-mean"
#define NCM_FIT_CATALOG_RTYPE_UNDEFINED "undefined-run"

G_END_DECLS

#endif /* _NCM_FIT_CATALOG_H_ */


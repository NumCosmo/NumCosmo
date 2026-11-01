/***************************************************************************
 *            ncm_mset_catalog.c
 *
 *  Tue February 18 10:49:26 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti & Mariana Penna Lima (January 2017)
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_mset_catalog.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br> & Mariana Penna Lima
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

/**
 * NcmMSetCatalog:
 *
 * Ordered table of free-parameter values of a #NcmMSet, kept in a
 * #NcmMSetCatalogStorage, with running statistics of its rows.
 *
 * Each row holds #NcmMSetCatalog:nadd-vals additional values ($-2\ln L$ and derived
 * functions, named by #NcmMSetCatalog:nadd-val-names; for a #NcmMSetCatalog:weighted
 * catalog the last one is the row weight) followed by the free parameters. With
 * #NcmMSetCatalog:nchains chains the rows are interleaved, the row with id $i$ belonging to chain
 * $i \bmod n_\mathrm{chains}$ (the walkers of an ensemble sampler), and the catalog
 * keeps statistics of every chain, of the chain means and of each ensemble. Row ids are
 * row indices; #NcmMSetCatalog:markovian-id marks where the Markov chain starts and
 * #NcmMSetCatalog:burnin drops leading rows when reading a file.
 *
 * The rows, the #NcmMSetCatalogDesc, the #NcmRNG of the run, the Markovian id and the
 * sampler segments live in the storage: a #NcmMSetCatalogStorageMemory for a catalog
 * created by ncm_mset_catalog_new(), a #NcmMSetCatalogStorageFits for one created by
 * ncm_mset_catalog_new_from_file() or given a file by ncm_mset_catalog_set_file(). The
 * catalog opens a recording on the storage at the first row added and closes it in
 * ncm_mset_catalog_reset(), ncm_mset_catalog_set_file() and when it is disposed; the
 * sampler strings (ncm_mset_catalog_set_sampler(), ncm_mset_catalog_set_sampler_options(),
 * ncm_mset_catalog_set_initial_sampler()) name the segments it opens. The rows are read
 * through ncm_mset_catalog_peek_row() and ncm_mset_catalog_get_row().
 *
 * #NcmFitMC stores the best fit of each realization, in the order of the realizations,
 * which reproduces a resampling from a seed; #NcmFitMCMC and #NcmFitESMCMC store their
 * chains. The running means and variances allow convergence to be followed row by row,
 * and the catalog provides autocorrelation times, effective sample sizes and other
 * convergence diagnostics.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_catalog.h"
#include "ncm/fit/ncm_mset_catalog_desc.h"
#include "ncm/fit/ncm_mset_catalog_storage.h"
#include "ncm/fit/ncm_mset_catalog_storage_memory.h"
#include "ncm/fit/ncm_mset_catalog_storage_fits.h"
#include "ncm/data/ncm_data_gauss_cov_mvnd.h"
#include "ncm/model/ncm_model_mvnd.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_func_eval.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_cdf.h>
#include <math.h>
#include <glib/gstdio.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmMSetCatalogPrivate
{
  NcmMSetCatalogStorage *storage;
  NcmMSetCatalogDesc *desc;
  NcmMSet *mset;
  NcmObjArray *functions;
  guint nadd_vals;
  gint m2lnp_var;
  NcmVector *bestfit_row;
  gdouble bestfit;
  gdouble post_lnnorm;
  gdouble post_lnnorm_sd;
  gboolean post_lnnorm_up;
  GArray *order_m2lnp;
  gboolean order_sorted;
  GPtrArray *add_vals_names;
  GPtrArray *add_vals_symbs;
  guint ncols;
  NcmStatsVec *pstats;
  NcmMSetCatalogSync smode;
  gboolean readonly;
  NcmRNG *rng;
  gboolean weighted;
  guint nchains;
  GPtrArray *chain_pstats;
  NcmStatsVec *mean_pstats;
  NcmStatsVec *e_stats;
  NcmStatsVec *e_mean_stats;
  GPtrArray *last_rows;
  guint read_block;
  GMutex row_lock;
  GPtrArray *e_mean_rows;
  NcmStatsRows *diag_rows;
  GPtrArray *e_var_array;
  gint naccepted;
  GArray *accept_ratio;
  NcmVector *chain_means;
  NcmVector *chain_vars;
  NcmMatrix *chain_cov;
  NcmMatrix *chain_sM;
  gsl_eigen_nonsymm_workspace *chain_sM_ws;
  gsl_vector_complex *chain_sM_ev;
  NcmStatsAcorr *acorr;
  NcmStatsAcorrMethod tau_method;
  NcmMSetCatalogPostNormMethod post_lnnorm_method;
  NcmVector *tau;
  GTimer *sync_timer;
  gdouble sync_interval;
  gchar *file;
  gchar *rtype_str;
  gchar *sampler_str;
  gchar *sampler_opts_str;
  gchar *init_sampler_str;
  gboolean init_recording;
  NcmVector *quantile_ws;
  gint cur_id;
  gint markovian_id;
  glong burnin;
  NcmVector *params_max;
  NcmVector *params_min;
  glong pdf_i;
  gsl_histogram *h;
  gsl_histogram_pdf *h_pdf;
  gboolean constructed;
} NcmMSetCatalogPrivate;

struct _NcmMSetCatalog
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSetCatalog, ncm_mset_catalog, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_MSET,
  PROP_NADD_VALS,
  PROP_ADD_VAL_NAMES,
  PROP_ADD_VAL_SYMBS,
  PROP_M2LNP_VAR,
  PROP_WEIGHTED,
  PROP_NCHAINS,
  PROP_BURNIN,
  PROP_MARKOVIAN_ID,
  PROP_TAU_METHOD,
  PROP_POST_LNNORM_METHOD,
  PROP_RNG,
  PROP_FILE,
  PROP_RUN_TYPE_STR,
  PROP_SAMPLER,
  PROP_SAMPLER_OPTS,
  PROP_INIT_SAMPLER,
  PROP_SYNC_MODE,
  PROP_SYNC_INTERVAL,
  PROP_READ_BLOCK,
  PROP_READONLY,
  PROP_STORAGE,
};

static gint
_ncm_mset_catalog_double_compare (gconstpointer a, gconstpointer b)
{
  const gdouble a_m2lnp = *(const gdouble *) a;
  const gdouble b_m2lnp = *(const gdouble *) b;

  return (a_m2lnp == b_m2lnp) ? 0.0 : ((a_m2lnp > b_m2lnp) ? +1 : -1);
}

static void
ncm_mset_catalog_init (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->storage        = NULL;
  self->desc           = NULL;
  self->mset           = NULL;
  self->functions      = NULL;
  self->nadd_vals      = 0;
  self->m2lnp_var      = -1;
  self->bestfit_row    = NULL;
  self->bestfit        = GSL_POSINF;
  self->post_lnnorm    = 0.0;
  self->post_lnnorm_sd = GSL_NAN;
  self->post_lnnorm_up = FALSE;
  self->order_m2lnp    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->order_sorted   = FALSE;
  self->add_vals_names = g_ptr_array_new_with_free_func (g_free);
  self->add_vals_symbs = g_ptr_array_new_with_free_func (g_free);
  self->ncols          = 0;
  self->pstats         = NULL;
  self->smode          = NCM_MSET_CATALOG_SYNC_LEN;
  self->readonly       = FALSE;
  self->rng            = NULL;
  self->weighted       = FALSE;
  self->nchains        = 0;
  self->chain_pstats   = g_ptr_array_new ();
  g_ptr_array_set_free_func (self->chain_pstats, (GDestroyNotify) & ncm_stats_vec_free);
  self->mean_pstats = NULL;
  self->e_var_array = g_ptr_array_new ();
  g_ptr_array_set_free_func (self->e_var_array, (GDestroyNotify) & ncm_vector_free);
  self->e_stats      = NULL;
  self->e_mean_stats = NULL;
  self->last_rows    = g_ptr_array_new ();
  self->read_block   = 0;
  g_mutex_init (&self->row_lock);
  self->e_mean_rows  = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vector_free);
  self->diag_rows    = NULL;
  self->naccepted    = 0;
  self->accept_ratio = NULL;
  self->chain_means  = NULL;
  self->chain_vars   = NULL;
  self->chain_cov    = NULL;
  self->chain_sM     = NULL;
  self->chain_sM_ws  = NULL;
  self->chain_sM_ev  = NULL;
  self->acorr        = NULL;
  self->tau_method   = NCM_STATS_ACORR_METHOD_MAX;

  self->post_lnnorm_method = NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX;
  self->tau                = NULL;

  self->sync_timer       = g_timer_new ();
  self->sync_interval    = 0.0;
  self->cur_id           = -1; /* the id of the last added row, the row index; -1 for no rows */
  self->markovian_id     = 0;  /* first row of the Markovian chain */
  self->burnin           = 0;  /* Number of elements to ignore when reading a catalog */
  self->file             = NULL;
  self->rtype_str        = NULL;
  self->sampler_str      = NULL;
  self->sampler_opts_str = NULL;
  self->init_sampler_str = NULL;
  self->init_recording   = FALSE;
  self->quantile_ws      = NULL;
  self->pdf_i            = -1;
  self->h                = NULL;
  self->h_pdf            = NULL;
  self->params_max       = NULL;
  self->params_min       = NULL;

  self->constructed = FALSE;
}

/* Releases the copies in last_rows and makes it n NULL entries. */
static void
_ncm_mset_catalog_reset_last_rows (NcmMSetCatalogPrivate *self, const guint n)
{
  guint i;

  for (i = 0; i < self->last_rows->len; i++)
    g_clear_pointer (&g_ptr_array_index (self->last_rows, i), ncm_vector_free);

  g_ptr_array_set_size (self->last_rows, n);
}

/*
 * Allocates the statistics for ncols columns and nchains chains: the running statistics
 * of every row, of each chain, of the chain means and of each ensemble, the
 * autocorrelation accumulator, the extrema of every column and the copies of the last
 * row of each chain that the acceptance count compares with.
 */
static void
_ncm_mset_catalog_alloc_stats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint free_params_len = ncm_mset_fparams_len (self->mset);
  const guint total           = self->ncols;
  guint i;

  _ncm_mset_catalog_reset_last_rows (self, self->nchains);

  self->pstats     = ncm_stats_vec_new (total, NCM_STATS_VEC_COV, FALSE);
  self->diag_rows  = ncm_stats_rows_new (total);
  self->params_max = ncm_vector_new (total);
  self->params_min = ncm_vector_new (total);

  ncm_vector_set_all (self->params_max, GSL_NEGINF);
  ncm_vector_set_all (self->params_min, GSL_POSINF);

  if (self->nchains > 1)
  {
    for (i = 0; i < self->nchains; i++)
    {
      NcmStatsVec *pstats = ncm_stats_vec_new (total, NCM_STATS_VEC_COV, FALSE);

      g_ptr_array_add (self->chain_pstats, pstats);
    }

    self->mean_pstats  = ncm_stats_vec_new (free_params_len, NCM_STATS_VEC_COV, FALSE);
    self->e_stats      = ncm_stats_vec_new (total, NCM_STATS_VEC_VAR, FALSE);
    self->e_mean_stats = ncm_stats_vec_new (total, NCM_STATS_VEC_VAR, FALSE);
    self->accept_ratio = g_array_new (FALSE, FALSE, sizeof (gdouble));

    self->chain_means = ncm_vector_new (self->nchains);
    self->chain_vars  = ncm_vector_new (self->nchains);
    self->chain_cov   = ncm_matrix_new (free_params_len, free_params_len);
    self->chain_sM    = ncm_matrix_new (free_params_len, free_params_len);
    self->chain_sM_ws = gsl_eigen_nonsymm_alloc (free_params_len);
    self->chain_sM_ev = gsl_vector_complex_alloc (free_params_len);
  }

  self->acorr = ncm_stats_acorr_new (total);
  ncm_stats_acorr_set_method (self->acorr, self->tau_method);

  self->tau = ncm_vector_new (total);
  ncm_vector_set_all (self->tau, 1.0);
}

/* Releases what _ncm_mset_catalog_alloc_stats() allocated, and the rows derived from it. */
static void
_ncm_mset_catalog_free_stats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  ncm_vector_clear (&self->bestfit_row);
  g_array_set_size (self->order_m2lnp, 0);
  self->order_sorted = FALSE;

  ncm_stats_vec_clear (&self->pstats);
  ncm_stats_rows_clear (&self->diag_rows);
  ncm_vector_clear (&self->params_max);
  ncm_vector_clear (&self->params_min);

  g_ptr_array_set_size (self->chain_pstats, 0);
  ncm_stats_vec_clear (&self->mean_pstats);
  ncm_stats_vec_clear (&self->e_stats);
  ncm_stats_vec_clear (&self->e_mean_stats);
  g_clear_pointer (&self->accept_ratio, g_array_unref);
  g_ptr_array_set_size (self->e_mean_rows, 0);
  g_ptr_array_set_size (self->e_var_array, 0);
  _ncm_mset_catalog_reset_last_rows (self, 0);
  self->naccepted = 0;

  ncm_vector_clear (&self->chain_means);
  ncm_vector_clear (&self->chain_vars);
  ncm_matrix_clear (&self->chain_cov);
  ncm_matrix_clear (&self->chain_sM);
  g_clear_pointer (&self->chain_sM_ws, gsl_eigen_nonsymm_free);
  g_clear_pointer (&self->chain_sM_ev, gsl_vector_complex_free);
  ncm_stats_acorr_clear (&self->acorr);
  ncm_vector_clear (&self->tau);
  ncm_vector_clear (&self->quantile_ws);
}

static void _ncm_mset_catalog_take_desc (NcmMSetCatalog *mcat, NcmMSetCatalogDesc *desc);
static NcmMSetCatalogStorageFits *_ncm_mset_catalog_peek_fits_storage (NcmMSetCatalog *mcat);
static void _ncm_mset_catalog_load_storage (NcmMSetCatalog *mcat);
static NcmMSetCatalogStorage *_ncm_mset_catalog_new_fits_storage (NcmMSetCatalog *mcat, const gchar *filename, gboolean read_only);
static NcmMSetCatalogStorage *_ncm_mset_catalog_new_memory_storage (NcmMSetCatalog *mcat);
static NcmMSetCatalogDesc *_ncm_mset_catalog_new_desc (NcmMSetCatalog *mcat, guint nchains, gint m2lnp_var, const gchar *rtype_str);

/*
 * A catalog given a storage takes its description, rows and run metadata. A catalog
 * without a model set is loaded from its file: the description, the model set, the rows
 * and the run metadata come from the file. A catalog with a model set builds its
 * description from its properties and stores it in a memory storage or, with a file, in a
 * FITS storage, binding it to the file's description when the file exists.
 */
static void
_ncm_mset_catalog_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->constructed (object);
  {
    NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
    NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

    self->constructed = TRUE;

    if (self->storage != NULL)
    {
      NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_peek_fits_storage (mcat);
      NcmMSetCatalogDesc *desc;
      GError *error = NULL;

      if ((self->mset != NULL) || (self->file != NULL))
        g_error ("_ncm_mset_catalog_constructed: a catalog built from a storage takes no NcmMSet and no file.");

      desc = ncm_mset_catalog_storage_peek_desc (self->storage, &error);

      if (error != NULL)
        g_error ("_ncm_mset_catalog_constructed: %s", error->message);

      if (desc == NULL)
        g_error ("_ncm_mset_catalog_constructed: the storage has no description.");

      g_object_get (self->storage, "smode", &self->smode, "sync-interval", &self->sync_interval, NULL);

      if (storage_fits != NULL)
      {
        g_object_get (storage_fits, "read-block", &self->read_block, NULL);
        self->readonly = ncm_mset_catalog_storage_fits_get_read_only (storage_fits);
      }

      _ncm_mset_catalog_take_desc (mcat, desc);
    }
    else if (self->mset == NULL)
    {
      NcmMSetCatalogDesc *desc;
      GError *error = NULL;

      if (self->file == NULL)
        g_error ("_ncm_mset_catalog_constructed: a catalog needs an NcmMSet or a file.");

      if (!g_file_test (self->file, G_FILE_TEST_EXISTS))
        g_error ("_ncm_mset_catalog_constructed: cannot create catalog file `%s' not found.",
                 self->file);

      self->storage = _ncm_mset_catalog_new_fits_storage (mcat, self->file, self->readonly);
      desc          = ncm_mset_catalog_storage_peek_desc (self->storage, &error);

      if (error != NULL)
        g_error ("_ncm_mset_catalog_constructed: %s", error->message);

      if (desc == NULL)
        g_error ("_ncm_mset_catalog_constructed: file `%s' holds no catalog.", self->file);

      _ncm_mset_catalog_take_desc (mcat, desc);
    }
    else
    {
      NcmMSetCatalogDesc *desc;
      GError *error = NULL;

      g_assert_cmpuint (self->add_vals_names->len, ==, self->add_vals_symbs->len);
      g_assert_cmpuint (self->add_vals_names->len, ==, self->nadd_vals);

      desc = _ncm_mset_catalog_new_desc (mcat, self->nchains, self->m2lnp_var, self->rtype_str);

      if (self->file != NULL)
        self->storage = _ncm_mset_catalog_new_fits_storage (mcat, self->file, self->readonly);
      else
        self->storage = _ncm_mset_catalog_new_memory_storage (mcat);

      ncm_mset_catalog_storage_set_desc (self->storage, desc, &error);

      if (error != NULL)
        g_error ("_ncm_mset_catalog_constructed: %s", error->message);

      _ncm_mset_catalog_take_desc (mcat, desc);
      ncm_mset_catalog_desc_free (desc);
    }

    g_clear_pointer (&self->file, g_free);

    _ncm_mset_catalog_load_storage (mcat);
  }
}

/*
 * Builds the description of the catalog from its model set and its user's additional
 * values (the names and symbols, which include the weight column once the catalog has a
 * description), with the given number of chains, -2ln(p) column and run type. Aborts on
 * an invalid input.
 */
static NcmMSetCatalogDesc *
_ncm_mset_catalog_new_desc (NcmMSetCatalog *mcat, guint nchains, gint m2lnp_var, const gchar *rtype_str)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint n_user          = self->add_vals_names->len - ((self->weighted && (self->desc != NULL)) ? 1 : 0);
  GPtrArray *names            = g_ptr_array_new ();
  GPtrArray *symbols          = g_ptr_array_new ();
  GError *error               = NULL;
  NcmMSetCatalogDesc *desc;
  guint i;

  for (i = 0; i < n_user; i++)
  {
    g_ptr_array_add (names, g_ptr_array_index (self->add_vals_names, i));
    g_ptr_array_add (symbols, g_ptr_array_index (self->add_vals_symbs, i));
  }

  g_ptr_array_add (names, NULL);
  g_ptr_array_add (symbols, NULL);

  desc = ncm_mset_catalog_desc_new (self->mset, nchains, self->weighted,
                                    (const gchar * const *) names->pdata,
                                    (const gchar * const *) symbols->pdata,
                                    m2lnp_var, rtype_str, &error);

  if (error != NULL)
    g_error ("_ncm_mset_catalog_new_desc: %s", error->message);

  g_ptr_array_unref (names);
  g_ptr_array_unref (symbols);

  return desc;
}

/*
 * Makes desc the description of the catalog: the model set when the catalog has none,
 * the number of columns, chains and additional values (the weight included), the
 * -2ln(p) column, the run type and the names and symbols of the additional values. The
 * statistics are allocated for it, replacing those of a previous description.
 */
static void
_ncm_mset_catalog_take_desc (NcmMSetCatalog *mcat, NcmMSetCatalogDesc *desc)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint i;

  ncm_mset_catalog_desc_clear (&self->desc);
  self->desc = ncm_mset_catalog_desc_ref (desc);

  if (self->mset == NULL)
    self->mset = ncm_mset_ref (ncm_mset_catalog_desc_peek_mset (desc));

  self->nadd_vals = ncm_mset_catalog_desc_nadd_vals (desc);
  self->ncols     = ncm_mset_catalog_desc_ncols (desc);
  self->nchains   = ncm_mset_catalog_desc_nchains (desc);
  self->weighted  = ncm_mset_catalog_desc_weighted (desc);
  self->m2lnp_var = ncm_mset_catalog_desc_get_m2lnp_var (desc);

  g_clear_pointer (&self->rtype_str, g_free);
  self->rtype_str = g_strdup (ncm_mset_catalog_desc_get_run_type (desc));

  g_ptr_array_set_size (self->add_vals_names, 0);
  g_ptr_array_set_size (self->add_vals_symbs, 0);

  for (i = 0; i < self->nadd_vals; i++)
  {
    g_ptr_array_add (self->add_vals_names, g_strdup (ncm_mset_catalog_desc_col_name (desc, i)));
    g_ptr_array_add (self->add_vals_symbs, g_strdup (ncm_mset_catalog_desc_col_symb (desc, i)));
  }

  g_assert_cmpuint (self->ncols, ==, ncm_mset_fparams_len (self->mset) + self->nadd_vals);

  _ncm_mset_catalog_free_stats (mcat);
  _ncm_mset_catalog_alloc_stats (mcat);
}

static void _ncm_mset_catalog_set_add_val_name_array (NcmMSetCatalog *mcat, gchar **names);
static void _ncm_mset_catalog_set_add_val_symbol_array (NcmMSetCatalog *mcat, gchar **symbols);

static void
_ncm_mset_catalog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_return_if_fail (NCM_IS_MSET_CATALOG (object));

  switch (prop_id)
  {
    case PROP_MSET:
      self->mset = g_value_dup_object (value);
      break;
    case PROP_NADD_VALS:
      self->nadd_vals = g_value_get_uint (value);
      break;
    case PROP_ADD_VAL_NAMES:
      _ncm_mset_catalog_set_add_val_name_array (mcat, g_value_get_boxed (value));
      break;
    case PROP_ADD_VAL_SYMBS:
      _ncm_mset_catalog_set_add_val_symbol_array (mcat, g_value_get_boxed (value));
      break;
    case PROP_M2LNP_VAR:
      ncm_mset_catalog_set_m2lnp_var (mcat, g_value_get_int (value));
      break;
    case PROP_WEIGHTED:
      self->weighted = g_value_get_boolean (value);
      break;
    case PROP_NCHAINS:
      self->nchains = g_value_get_uint (value);
      break;
    case PROP_BURNIN:
      ncm_mset_catalog_set_burnin (mcat, g_value_get_long (value));
      break;
    case PROP_MARKOVIAN_ID:
      ncm_mset_catalog_set_markovian_id (mcat, g_value_get_int (value));
      break;
    case PROP_TAU_METHOD:
      ncm_mset_catalog_set_tau_method (mcat, g_value_get_enum (value));
      break;
    case PROP_POST_LNNORM_METHOD:
      ncm_mset_catalog_set_post_lnnorm_method (mcat, g_value_get_enum (value));
      break;
    case PROP_RNG:
      ncm_mset_catalog_set_rng (mcat, g_value_get_object (value));
      break;
    case PROP_FILE:
      ncm_mset_catalog_set_file (mcat, g_value_get_string (value));
      break;
    case PROP_RUN_TYPE_STR:
      ncm_mset_catalog_set_run_type (mcat, g_value_get_string (value));
      break;
    case PROP_SAMPLER:
      ncm_mset_catalog_set_sampler (mcat, g_value_get_string (value));
      break;
    case PROP_SAMPLER_OPTS:
      ncm_mset_catalog_set_sampler_options (mcat, g_value_get_string (value));
      break;
    case PROP_INIT_SAMPLER:
      ncm_mset_catalog_set_initial_sampler (mcat, g_value_get_string (value));
      break;
    case PROP_SYNC_MODE:
      ncm_mset_catalog_set_sync_mode (mcat, g_value_get_enum (value));
      break;
    case PROP_SYNC_INTERVAL:
      ncm_mset_catalog_set_sync_interval (mcat, g_value_get_double (value));
      break;
    case PROP_READ_BLOCK:
      ncm_mset_catalog_set_read_block (mcat, g_value_get_uint (value));
      break;
    case PROP_READONLY:
      self->readonly = g_value_get_boolean (value);
      break;
    case PROP_STORAGE:
      self->storage = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_return_if_fail (NCM_IS_MSET_CATALOG (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, self->mset);
      break;
    case PROP_NADD_VALS:
      g_value_set_uint (value, self->nadd_vals);
      break;
    case PROP_ADD_VAL_NAMES:
    {
      gchar **names = g_new (gchar *, self->add_vals_names->len + 1);
      guint i;

      for (i = 0; i < self->add_vals_names->len; i++)
      {
        names[i] = g_strdup (g_ptr_array_index (self->add_vals_names, i));
      }

      names[i] = NULL;

      g_value_take_boxed (value, names);

      break;
    }
    case PROP_ADD_VAL_SYMBS:
    {
      gchar **symbs = g_new (gchar *, self->add_vals_symbs->len + 1);
      guint i;

      for (i = 0; i < self->add_vals_symbs->len; i++)
      {
        symbs[i] = g_strdup (g_ptr_array_index (self->add_vals_symbs, i));
      }

      symbs[i] = NULL;

      g_value_take_boxed (value, symbs);

      break;
    }
    case PROP_M2LNP_VAR:
      g_value_set_int (value, ncm_mset_catalog_get_m2lnp_var (mcat));
      break;
    case PROP_WEIGHTED:
      g_value_set_boolean (value, self->weighted);
      break;
    case PROP_NCHAINS:
      g_value_set_uint (value, self->nchains);
      break;
    case PROP_BURNIN:
      g_value_set_long (value, ncm_mset_catalog_get_burnin (mcat));
      break;
    case PROP_MARKOVIAN_ID:
      g_value_set_int (value, ncm_mset_catalog_get_markovian_id (mcat));
      break;
    case PROP_TAU_METHOD:
      g_value_set_enum (value, ncm_mset_catalog_get_tau_method (mcat));
      break;
    case PROP_POST_LNNORM_METHOD:
      g_value_set_enum (value, ncm_mset_catalog_get_post_lnnorm_method (mcat));
      break;
    case PROP_RNG:
      g_value_set_object (value, self->rng);
      break;
    case PROP_FILE:
      g_value_set_string (value, ncm_mset_catalog_peek_filename (mcat));
      break;
    case PROP_RUN_TYPE_STR:
      g_value_set_string (value, self->rtype_str);
      break;
    case PROP_SAMPLER:
      g_value_set_string (value, self->sampler_str);
      break;
    case PROP_SAMPLER_OPTS:
      g_value_set_string (value, self->sampler_opts_str);
      break;
    case PROP_INIT_SAMPLER:
      g_value_set_string (value, self->init_sampler_str);
      break;
    case PROP_SYNC_MODE:
      g_value_set_enum (value, self->smode);
      break;
    case PROP_SYNC_INTERVAL:
      g_value_set_double (value, self->sync_interval);
      break;
    case PROP_READ_BLOCK:
      g_value_set_uint (value, self->read_block);
      break;
    case PROP_READONLY:
      g_value_set_boolean (value, self->readonly);
      break;
    case PROP_STORAGE:
      g_value_set_object (value, self->storage);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_mset_catalog_end_recording (NcmMSetCatalog *mcat);

static void
_ncm_mset_catalog_dispose (GObject *object)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->storage != NULL)
    _ncm_mset_catalog_end_recording (mcat);

  _ncm_mset_catalog_free_stats (mcat);
  g_clear_pointer (&self->order_m2lnp, g_array_unref);

  ncm_mset_catalog_storage_clear (&self->storage);
  ncm_mset_catalog_desc_clear (&self->desc);
  ncm_mset_clear (&self->mset);
  ncm_obj_array_clear (&self->functions);
  ncm_rng_clear (&self->rng);

  g_clear_pointer (&self->chain_pstats, g_ptr_array_unref);
  g_clear_pointer (&self->e_mean_rows, g_ptr_array_unref);
  g_clear_pointer (&self->e_var_array, g_ptr_array_unref);
  g_clear_pointer (&self->last_rows, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->dispose (object);
}

static void
_ncm_mset_catalog_finalize (GObject *object)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->h != NULL)
    gsl_histogram_free (self->h);

  if (self->h_pdf != NULL)
    gsl_histogram_pdf_free (self->h_pdf);

  g_clear_pointer (&self->rtype_str, g_free);
  g_clear_pointer (&self->sampler_str, g_free);
  g_clear_pointer (&self->sampler_opts_str, g_free);
  g_clear_pointer (&self->init_sampler_str, g_free);

  g_timer_destroy (self->sync_timer);
  g_mutex_clear (&self->row_lock);

  g_ptr_array_unref (self->add_vals_names);
  g_ptr_array_unref (self->add_vals_symbs);

  g_clear_pointer (&self->file, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->finalize (object);
}

static void
ncm_mset_catalog_class_init (NcmMSetCatalogClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_mset_catalog_constructed;
  object_class->set_property = &_ncm_mset_catalog_set_property;
  object_class->get_property = &_ncm_mset_catalog_get_property;
  object_class->dispose      = &_ncm_mset_catalog_dispose;
  object_class->finalize     = &_ncm_mset_catalog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "NcmMSet object",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NADD_VALS,
                                   g_param_spec_uint ("nadd-vals",
                                                      NULL,
                                                      "Number of additional values",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ADD_VAL_NAMES,
                                   g_param_spec_boxed ("nadd-val-names",
                                                       NULL,
                                                       "Additional value names",
                                                       G_TYPE_STRV,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ADD_VAL_SYMBS,
                                   g_param_spec_boxed ("nadd-val-symbols",
                                                       NULL,
                                                       "Additional value symbols",
                                                       G_TYPE_STRV,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M2LNP_VAR,
                                   g_param_spec_int ("m2lnp-var",
                                                     NULL,
                                                     "Index of the variable representing m2lnp",
                                                     -1, G_MAXINT, -1,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_WEIGHTED,
                                   g_param_spec_boolean ("weighted",
                                                         NULL,
                                                         "Catalog with weighted rows",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NCHAINS,
                                   g_param_spec_uint ("nchains",
                                                      NULL,
                                                      "Number of different chains in the catalog",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_BURNIN,
                                   g_param_spec_long ("burnin",
                                                      NULL,
                                                      "Burn-in size",
                                                      0, G_MAXLONG, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:markovian-id:
   *
   * Id of the first row produced by a Markovian step after the last move that did not
   * satisfy detailed balance (the initial ensemble, an exploration phase of the walker).
   * Rows from this id on form a valid Markov chain whose initial state is the ensemble of
   * the #NcmMSetCatalog:nchains rows before it. A row index, like the current id; a file
   * without the key reads as 0 (every row Markovian). When a file is read with a burn-in
   * the stored value is shifted by it, like every other row id.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MARKOVIAN_ID,
                                   g_param_spec_int ("markovian-id",
                                                     NULL,
                                                     "First row id of the Markovian chain",
                                                     0, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_TAU_METHOD,
                                   g_param_spec_enum ("tau-method",
                                                      NULL,
                                                      "Method used to calculate the autocorrelation time",
                                                      NCM_TYPE_STATS_ACORR_METHOD, NCM_STATS_ACORR_METHOD_MAX,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_POST_LNNORM_METHOD,
                                   g_param_spec_enum ("post-lnnorm-method",
                                                      NULL,
                                                      "Method used to estimate the log evidence",
                                                      NCM_TYPE_MSET_CATALOG_POST_NORM_METHOD, NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RNG,
                                   g_param_spec_object ("rng",
                                                        NULL,
                                                        "Random number generator object",
                                                        NCM_TYPE_RNG,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:smode:
   *
   * When the storage writes the rows added, see #NcmMSetCatalogStorage:smode.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SYNC_MODE,
                                   g_param_spec_enum ("smode",
                                                      NULL,
                                                      "Catalog sync mode",
                                                      NCM_TYPE_MSET_CATALOG_SYNC, NCM_MSET_CATALOG_SYNC_AUTO,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:filename:
   *
   * The FITS file of the catalog, see ncm_mset_catalog_set_file(); %NULL for a catalog
   * in memory.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_FILE,
                                   g_param_spec_string ("filename",
                                                        NULL,
                                                        "Catalog filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RUN_TYPE_STR,
                                   g_param_spec_string ("run-type-string",
                                                        NULL,
                                                        "Run type string",
                                                        NCM_MSET_CATALOG_RTYPE_UNDEFINED,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:sampler:
   *
   * A description of the sampler that fills the rows: the name of the sampler segment
   * the catalog opens on its storage, see ncm_mset_catalog_set_sampler(). It is never
   * compared: a catalog may be continued with any sampler.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SAMPLER,
                                   g_param_spec_string ("sampler",
                                                        NULL,
                                                        "Description of the sampler that filled the rows",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:sampler-options:
   *
   * The sampler's tunable settings, as a colon separated list of `name=value' pairs: the
   * options of the sampler segment the catalog opens on its storage. Kept separate from
   * #NcmMSetCatalog:sampler, so that the structure and its tuning can be read apart:
   * two runs of the same sampler differ here and nowhere else.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SAMPLER_OPTS,
                                   g_param_spec_string ("sampler-options",
                                                        NULL,
                                                        "Tunable settings of the sampler that filled the rows",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:initial-sampler:
   *
   * A description of the sampler that draws the initial points: the name of the first
   * sampler segment of an empty storage, which covers the first ensemble
   * (#NcmMSetCatalog:nchains rows), see ncm_mset_catalog_set_initial_sampler().
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_INIT_SAMPLER,
                                   g_param_spec_string ("initial-sampler",
                                                        NULL,
                                                        "Description of the sampler that drew the initial points",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SYNC_INTERVAL,
                                   g_param_spec_double ("sync-interval",
                                                        NULL,
                                                        "Data sync interval",
                                                        0.0, 1.0e3, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:read-block:
   *
   * Number of consecutive rows read from the file at once,
   * #NcmMSetCatalogStorageFits:read-block of the FITS storage. The storage keeps in
   * memory the rows not yet written and the last two ensembles; see
   * ncm_mset_catalog_peek_row(). A value at least as large as the catalog holds all its
   * rows once they are read.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_READ_BLOCK,
                                   g_param_spec_uint ("read-block",
                                                      NULL,
                                                      "Number of rows read from the file at once",
                                                      1, G_MAXUINT, 1024,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:read-only:
   *
   * Whether the file given at construction is opened read-only. A read-only catalog
   * refuses every write to its storage; ncm_mset_catalog_set_file() gives it a writable
   * storage.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_READONLY,
                                   g_param_spec_boolean ("read-only",
                                                         NULL,
                                                         "Whether the FITS file is opened read-only",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalog:storage:
   *
   * The #NcmMSetCatalogStorage the catalog is built on, see
   * ncm_mset_catalog_new_from_storage(). Given, it must hold a description, and the
   * catalog takes no #NcmMSetCatalog:mset and no #NcmMSetCatalog:filename. The catalog's
   * sync mode, sync interval and read block are then the storage's.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_STORAGE,
                                   g_param_spec_object ("storage",
                                                        NULL,
                                                        "The storage of the catalog",
                                                        NCM_TYPE_MSET_CATALOG_STORAGE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_mset_catalog_new:
 * @mset: a #NcmMSet
 * @nadd_vals: number of additional values
 * @nchains: number of chains, at least one
 * @weighted: whether the rows are weighted
 * @...: a name and a symbol for each additional value
 *
 * Creates a #NcmMSetCatalog for the free parameters of @mset, which must not change
 * while the catalog is in use, on a #NcmMSetCatalogStorageMemory. With more than one
 * chain the statistics of each chain are kept separately.
 *
 * Returns: (transfer full): a new #NcmMSetCatalog
 */
NcmMSetCatalog *
ncm_mset_catalog_new (NcmMSet *mset, guint nadd_vals, guint nchains, gboolean weighted, ...)
{
  va_list ap;
  gchar **names   = g_new (gchar *, nadd_vals + 1);
  gchar **symbols = g_new (gchar *, nadd_vals + 1);
  guint i;

  va_start (ap, weighted);

  for (i = 0; i < nadd_vals; i++)
  {
    gchar *name   = va_arg (ap, gchar *);
    gchar *symbol = (name != NULL) ? va_arg (ap, gchar *) : NULL;

    if ((name == NULL) || (symbol == NULL))
      g_error ("ncm_mset_catalog_new: missing %u-th name and or symbol.", i + 1);

    names[i]   = name;
    symbols[i] = symbol;
  }

  va_end (ap);

  names[i]   = NULL;
  symbols[i] = NULL;

  {
    NcmMSetCatalog *mcat = g_object_new (NCM_TYPE_MSET_CATALOG,
                                         "mset",             mset,
                                         "nadd-vals",        nadd_vals,
                                         "nchains",          nchains,
                                         "weighted",         weighted,
                                         "nadd-val-names",   names,
                                         "nadd-val-symbols", symbols,
                                         NULL);

    g_free (names);
    g_free (symbols);

    return mcat;
  }
}

/**
 * ncm_mset_catalog_new_array:
 * @mset: a #NcmMSet
 * @nadd_vals: number of additional values
 * @nchains: number of chains, at least one
 * @weighted: whether the rows are weighted
 * @names: (array zero-terminated=1): names of the additional values
 * @symbols: (array zero-terminated=1): symbols of the additional values
 *
 * Creates a #NcmMSetCatalog as ncm_mset_catalog_new(), with the names and symbols as
 * arrays.
 *
 * Returns: (transfer full): a new #NcmMSetCatalog
 */
NcmMSetCatalog *
ncm_mset_catalog_new_array (NcmMSet *mset, guint nadd_vals, guint nchains, gboolean weighted, gchar **names, gchar **symbols)
{
  NcmMSetCatalog *mcat = g_object_new (NCM_TYPE_MSET_CATALOG,
                                       "mset",             mset,
                                       "nadd-vals",        nadd_vals,
                                       "nchains",          nchains,
                                       "weighted",         weighted,
                                       "nadd-val-names",   names,
                                       "nadd-val-symbols", symbols,
                                       NULL);

  return mcat;
}

/**
 * ncm_mset_catalog_new_from_file:
 * @filename: catalog FITS file name
 * @burnin: number of leading rows to drop
 *
 * Creates a #NcmMSetCatalog on a writable #NcmMSetCatalogStorageFits of @filename, which
 * must exist. The description, the model set, the rows and the run metadata are read
 * from the file (see #NcmMSetCatalogStorageFits for the files of the previous library).
 * @burnin is applied as in ncm_mset_catalog_set_burnin().
 *
 * Returns: (transfer full): a new #NcmMSetCatalog
 */
NcmMSetCatalog *
ncm_mset_catalog_new_from_file (const gchar *filename, glong burnin)
{
  NcmMSetCatalog *mcat = g_object_new (NCM_TYPE_MSET_CATALOG,
                                       "filename", filename,
                                       "burnin", burnin,
                                       NULL);

  return mcat;
}

/**
 * ncm_mset_catalog_new_from_file_ro:
 * @filename: catalog FITS file name
 * @burnin: number of leading rows to drop
 *
 * Creates a #NcmMSetCatalog as ncm_mset_catalog_new_from_file(), opening @filename
 * read-only.
 *
 * Returns: (transfer full): a new #NcmMSetCatalog
 */
NcmMSetCatalog *
ncm_mset_catalog_new_from_file_ro (const gchar *filename, glong burnin)
{
  NcmMSetCatalog *mcat = g_object_new (NCM_TYPE_MSET_CATALOG,
                                       "filename", filename,
                                       "read-only", TRUE,
                                       "burnin", burnin,
                                       NULL);

  return mcat;
}

/**
 * ncm_mset_catalog_new_from_storage:
 * @storage: a #NcmMSetCatalogStorage holding a description
 * @burnin: number of leading rows to drop
 *
 * Creates a #NcmMSetCatalog on @storage: its description, rows and run metadata are
 * the storage's. @burnin is applied as in ncm_mset_catalog_set_burnin(). A storage
 * without a description aborts.
 *
 * Returns: (transfer full): a new #NcmMSetCatalog
 */
NcmMSetCatalog *
ncm_mset_catalog_new_from_storage (NcmMSetCatalogStorage *storage, glong burnin)
{
  NcmMSetCatalog *mcat = g_object_new (NCM_TYPE_MSET_CATALOG,
                                       "storage", storage,
                                       "burnin", burnin,
                                       NULL);

  return mcat;
}

/**
 * ncm_mset_catalog_peek_storage:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer none): the #NcmMSetCatalogStorage of @mcat
 */
NcmMSetCatalogStorage *
ncm_mset_catalog_peek_storage (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->storage;
}

/* A read-only FITS storage of filename, with its description read. Aborts on error. */
static NcmMSetCatalogStorageFits *
_ncm_mset_catalog_open_ro (const gchar *filename, const gchar *func)
{
  GError *error                           = NULL;
  NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (filename, TRUE, &error);

  if (error != NULL)
    g_error ("%s: %s", func, error->message);

  if (ncm_mset_catalog_storage_peek_desc (NCM_MSET_CATALOG_STORAGE (storage_fits), NULL) == NULL)
    g_error ("%s: file `%s' holds no catalog.", func, filename);

  return storage_fits;
}

/**
 * ncm_mset_catalog_peek_info_from_file:
 * @filename: filename of the catalog fits
 * @nrows: (out): number of rows (samples) currently in the catalog
 * @nchains: (out): number of chains (walkers) in the catalog
 *
 * Reads the number of rows and of chains of @filename from its header, without loading
 * the catalog, e.g. to convert a burn-in from iterations to rows (see
 * ncm_mset_catalog_set_burnin()) before ncm_mset_catalog_new_from_file_ro().
 *
 */
void
ncm_mset_catalog_peek_info_from_file (const gchar *filename, glong *nrows, guint *nchains)
{
  NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_open_ro (filename, "ncm_mset_catalog_peek_info_from_file");
  NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);
  GError *error                           = NULL;

  *nrows   = ncm_mset_catalog_storage_nrows (storage, &error);
  *nchains = ncm_mset_catalog_desc_nchains (ncm_mset_catalog_storage_peek_desc (storage, NULL));

  if (error != NULL)
    g_error ("ncm_mset_catalog_peek_info_from_file: %s", error->message);

  ncm_mset_catalog_storage_fits_free (storage_fits);
}

/**
 * ncm_mset_catalog_peek_markovian_id_from_file:
 * @filename: catalog file name
 *
 * Reads #NcmMSetCatalog:markovian-id from the file header without loading the catalog,
 * as a row index. A file without the key returns 0 (every row Markovian).
 *
 * Returns: the index of the first row of the Markovian chain.
 */
gint
ncm_mset_catalog_peek_markovian_id_from_file (const gchar *filename)
{
  NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_open_ro (filename, "ncm_mset_catalog_peek_markovian_id_from_file");
  NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);
  GError *error                           = NULL;
  const guint markid                      = ncm_mset_catalog_storage_get_markovian_id (storage, &error);

  if (error != NULL)
    g_error ("ncm_mset_catalog_peek_markovian_id_from_file: %s", error->message);

  ncm_mset_catalog_storage_fits_free (storage_fits);

  return (gint) markid;
}

/**
 * ncm_mset_catalog_ref:
 * @mcat: a #NcmMSetCatalog
 *
 * Increases the reference count of @mcat atomically.
 *
 * Returns: (transfer full): @mcat.
 */
NcmMSetCatalog *
ncm_mset_catalog_ref (NcmMSetCatalog *mcat)
{
  return g_object_ref (mcat);
}

/**
 * ncm_mset_catalog_free:
 * @mcat: a #NcmMSetCatalog
 *
 * Saves @mcat to its file (ncm_mset_catalog_sync()) and decreases its reference count.
 *
 */
void
ncm_mset_catalog_free (NcmMSetCatalog *mcat)
{
  ncm_mset_catalog_sync (mcat, TRUE);
  g_object_unref (mcat);
}

/**
 * ncm_mset_catalog_clear:
 * @mcat: a #NcmMSetCatalog
 *
 * Saves *@mcat to its file, when not %NULL, decreases its reference count and sets it to
 * %NULL.
 *
 */
void
ncm_mset_catalog_clear (NcmMSetCatalog **mcat)
{
  if (*mcat != NULL)
    ncm_mset_catalog_sync (*mcat, TRUE);

  g_clear_object (mcat);
}

/*
 * A FITS storage of filename with the catalog's sync mode, sync interval and read block,
 * keeping the last two ensembles in memory; the file is opened when it exists. Aborts
 * when it cannot be read.
 */
static NcmMSetCatalogStorage *
_ncm_mset_catalog_new_fits_storage (NcmMSetCatalog *mcat, const gchar *filename, gboolean read_only)
{
  NcmMSetCatalogPrivate *self    = ncm_mset_catalog_get_instance_private (mcat);
  GError *error                  = NULL;
  NcmMSetCatalogStorage *storage = g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_FITS,
                                                 "filename", filename,
                                                 "read-only", read_only,
                                                 "read-block", self->read_block,
                                                 "keep-last", 2,
                                                 "smode", self->smode,
                                                 "sync-interval", self->sync_interval,
                                                 NULL);

  ncm_mset_catalog_storage_peek_desc (storage, &error);

  if (error != NULL)
    g_error ("ncm_mset_catalog: cannot open catalog file `%s': %s", filename, error->message);

  return storage;
}

/* A memory storage with the catalog's sync mode and sync interval. */
static NcmMSetCatalogStorage *
_ncm_mset_catalog_new_memory_storage (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_MEMORY,
                       "smode", self->smode,
                       "sync-interval", self->sync_interval,
                       NULL);
}

/* The storage as a FITS storage, NULL for any other storage. */
static NcmMSetCatalogStorageFits *
_ncm_mset_catalog_peek_fits_storage (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if ((self->storage != NULL) && NCM_IS_MSET_CATALOG_STORAGE_FITS (self->storage))
    return NCM_MSET_CATALOG_STORAGE_FITS (self->storage);
  else
    return NULL;
}

/* Aborts with the message of error, prefixed by func, when error is set. */
static void
_ncm_mset_catalog_check_error (GError **error, const gchar *func)
{
  if (*error != NULL)
    g_error ("%s: %s", func, (*error)->message);
}

static void _ncm_mset_catalog_post_update_stats (NcmMSetCatalog *mcat, NcmVector *x);

/*
 * Takes the rows and the run metadata of the storage: the burn-in, the statistics of the
 * rows after the burn-in, the Markovian id, the RNG and the sampler strings. The catalog
 * must hold no statistics.
 */
static void
_ncm_mset_catalog_load_storage (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;
  const guint nrows           = ncm_mset_catalog_storage_nrows (self->storage, &error);
  const gchar *filename       = ncm_mset_catalog_peek_filename (mcat);
  guint markid;
  guint i;

  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");

  g_assert_cmpuint (ncm_stats_vec_nitens (self->pstats), ==, 0);

  if (self->burnin % self->nchains != 0)
    g_error ("_ncm_mset_catalog_load_storage: burnin of %ld row(s) is not a multiple of "
             "the number of chains %u of catalog `%s'.",
             self->burnin, self->nchains, filename != NULL ? filename : "(memory)");

  if (nrows < (guint) self->burnin)
    g_error ("_ncm_mset_catalog_load_storage: burnin of %ld row(s) exceeds "
             "catalog `%s' (%u row(s), %u chain(s), i.e. %u iteration(s)).",
             self->burnin, filename != NULL ? filename : "(memory)", nrows, self->nchains, nrows / self->nchains);

  self->cur_id = -1;

  for (i = self->burnin; i < nrows; i++)
  {
    NcmVector *row = ncm_mset_catalog_storage_peek_row (self->storage, i, &error);

    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");
    _ncm_mset_catalog_post_update_stats (mcat, row);
  }

  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, nrows - self->burnin);

  /* Rows before the burn-in are dropped and the ids renumbered from 0, so the stored
   * Markovian id shifts by the same amount; it is never below 0. */
  markid = ncm_mset_catalog_storage_get_markovian_id (self->storage, &error);
  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");
  self->markovian_id = MAX (self->markovian_id, (gint) markid - (gint) self->burnin);

  {
    NcmRNG *rng = ncm_mset_catalog_storage_get_rng (self->storage, &error);

    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");

    if (rng != NULL)
    {
      /* The RNG of the file continues the run; one set before only has to be the same
       * generator (algorithm and seed), the state saved with the rows is what counts. */
      if ((self->rng != NULL) && (self->rng != rng))
        if ((strcmp (ncm_rng_get_algo (self->rng), ncm_rng_get_algo (rng)) != 0) ||
            (ncm_rng_get_seed (self->rng) != ncm_rng_get_seed (rng)))
          g_error ("_ncm_mset_catalog_load_storage: catalog `%s' holds an RNG (%s, seed %lu) and the "
                   "catalog has another (%s, seed %lu); do not set an RNG when resuming from a file.",
                   filename, ncm_rng_get_algo (rng), ncm_rng_get_seed (rng),
                   ncm_rng_get_algo (self->rng), ncm_rng_get_seed (self->rng));


      ncm_rng_clear (&self->rng);
      self->rng = rng;
    }
    else if (self->rng != NULL)
    {
      ncm_mset_catalog_storage_set_rng (self->storage, self->rng, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");
    }
  }

  {
    const guint n = ncm_mset_catalog_storage_n_samplers (self->storage, &error);

    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");

    if (n > 0)
    {
      gchar *sampler = NULL;
      gchar *options = NULL;

      ncm_mset_catalog_storage_get_sampler (self->storage, n - 1, NULL, &sampler, &options, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");

      if (self->sampler_str == NULL)
      {
        self->sampler_str = sampler;
        sampler           = NULL;
      }

      if ((self->sampler_opts_str == NULL) && (strlen (options) > 0))
      {
        self->sampler_opts_str = options;
        options                = NULL;
      }

      g_free (sampler);
      g_free (options);
    }

    if ((n > 1) && (self->init_sampler_str == NULL))
    {
      ncm_mset_catalog_storage_get_sampler (self->storage, 0, NULL, &self->init_sampler_str, NULL, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_load_storage");
    }
  }
}

/*
 * Writes the catalog's Markovian id to the storage when it is beyond the stored one; only
 * inside a recording.
 */
static void
_ncm_mset_catalog_push_markovian_id (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;
  const gint markid           = ncm_mset_catalog_get_markovian_id (mcat) + self->burnin;
  guint stored;

  if (!ncm_mset_catalog_storage_is_recording (self->storage))
    return;

  stored = ncm_mset_catalog_storage_get_markovian_id (self->storage, &error);
  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_push_markovian_id");

  if (markid > (gint) stored)
  {
    ncm_mset_catalog_storage_set_markovian_id (self->storage, markid, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_push_markovian_id");
  }
}

/*
 * Opens the recording the next row is added to, when none is open (a sampler may have
 * opened one on the storage): the initial sampler's segment on an empty storage with an
 * initial sampler, the sampler's segment otherwise (named by the run type when no sampler
 * was set). A pending Markovian id is written.
 */
static void
_ncm_mset_catalog_ensure_recording (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;
  guint nrows;

  if (ncm_mset_catalog_storage_is_recording (self->storage))
  {
    _ncm_mset_catalog_push_markovian_id (mcat);

    return;
  }

  nrows = ncm_mset_catalog_storage_nrows (self->storage, &error);
  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_ensure_recording");

  if ((nrows == 0) && (self->init_sampler_str != NULL))
  {
    ncm_mset_catalog_storage_start_recording (self->storage, self->init_sampler_str, NULL, &error);
    self->init_recording = TRUE;
  }
  else
  {
    ncm_mset_catalog_storage_start_recording (self->storage,
                                              (self->sampler_str != NULL) ? self->sampler_str : self->rtype_str,
                                              self->sampler_opts_str, &error);
    self->init_recording = FALSE;
  }

  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_ensure_recording");

  _ncm_mset_catalog_push_markovian_id (mcat);
}

/* Closes the open recording, if any, writing what it holds. */
static void
_ncm_mset_catalog_end_recording (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;

  if (!ncm_mset_catalog_storage_is_recording (self->storage))
    return;

  ncm_mset_catalog_storage_end_recording (self->storage, &error);
  self->init_recording = FALSE;
  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_end_recording");
}

/*
 * Appends row to the storage, inside the recording of _ncm_mset_catalog_ensure_recording(),
 * and feeds the statistics with the storage's copy of it, so that the references the
 * statistics keep are to the storage's vectors. The initial sampler's recording is closed
 * when the first ensemble is complete, so that the sampler's segment starts at the next
 * row.
 */
static void
_ncm_mset_catalog_add_row (NcmMSetCatalog *mcat, NcmVector *row)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;

  _ncm_mset_catalog_ensure_recording (mcat);

  ncm_mset_catalog_storage_add_row (self->storage, row, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_add_row");

  {
    const guint nrows = ncm_mset_catalog_storage_nrows (self->storage, &error);
    NcmVector *stored;

    _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_add_row");
    stored = ncm_mset_catalog_storage_peek_row (self->storage, nrows - 1, &error);
    _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_add_row");

    _ncm_mset_catalog_post_update_stats (mcat, stored);
  }

  if (self->init_recording && (ncm_mset_catalog_len (mcat) == self->nchains))
    _ncm_mset_catalog_end_recording (mcat);
}

/*
 * Copies the rows of src after the catalog's burn-in into dst, one recording per sampler
 * segment of src (one recording named by the sampler string, or the run type, when src
 * has no segment), and writes the Markovian id. dst must hold the description and no
 * rows.
 */
static void
_ncm_mset_catalog_copy_rows (NcmMSetCatalog *mcat, NcmMSetCatalogStorage *src, NcmMSetCatalogStorage *dst)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;
  const guint nrows           = ncm_mset_catalog_storage_nrows (src, &error);
  const guint nsegs           = ncm_mset_catalog_storage_n_samplers (src, &error);
  const gint markid           = ncm_mset_catalog_get_markovian_id (mcat);
  gboolean recording          = FALSE;
  guint i;

  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");

  if (nrows <= (guint) self->burnin)
    return;

  g_object_set (dst, "smode", NCM_MSET_CATALOG_SYNC_DISABLE, NULL);

  if (nsegs == 0)
  {
    ncm_mset_catalog_storage_start_recording (dst, (self->sampler_str != NULL) ? self->sampler_str : self->rtype_str,
                                              self->sampler_opts_str, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
    recording = TRUE;

    for (i = self->burnin; i < nrows; i++)
    {
      NcmVector *row = ncm_mset_catalog_storage_peek_row (src, i, &error);

      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
      ncm_mset_catalog_storage_add_row (dst, row, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
    }
  }
  else
  {
    guint s;

    for (s = 0; s < nsegs; s++)
    {
      gchar *sampler = NULL;
      gchar *options = NULL;
      guint start    = 0;
      guint end      = nrows;

      ncm_mset_catalog_storage_get_sampler (src, s, &start, &sampler, &options, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");

      if (s + 1 < nsegs)
      {
        ncm_mset_catalog_storage_get_sampler (src, s + 1, &end, NULL, NULL, &error);
        _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
      }

      start = MAX (start, (guint) self->burnin);

      if (start < end)
      {
        if (recording)
        {
          ncm_mset_catalog_storage_end_recording (dst, &error);
          _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
        }

        ncm_mset_catalog_storage_start_recording (dst, sampler, options, &error);
        _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
        recording = TRUE;

        for (i = start; i < end; i++)
        {
          NcmVector *row = ncm_mset_catalog_storage_peek_row (src, i, &error);

          _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
          ncm_mset_catalog_storage_add_row (dst, row, &error);
          _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
        }
      }

      g_free (sampler);
      g_free (options);
    }
  }

  g_assert (recording);

  if (markid > 0)
  {
    ncm_mset_catalog_storage_set_markovian_id (dst, markid, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");
  }

  ncm_mset_catalog_storage_end_recording (dst, &error);
  _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_copy_rows");

  g_object_set (dst, "smode", self->smode, NULL);
}

/*
 * Replaces the storage by new_storage, which holds no description: the description, the
 * RNG and the rows after the burn-in are copied into it, and the burn-in is then zero.
 * The previous storage is released; a FITS one is left as it was.
 */
static void
_ncm_mset_catalog_switch_storage (NcmMSetCatalog *mcat, NcmMSetCatalogStorage *new_storage)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;

  _ncm_mset_catalog_end_recording (mcat);

  ncm_mset_catalog_storage_set_desc (new_storage, self->desc, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_set_file");

  if (self->rng != NULL)
  {
    ncm_mset_catalog_storage_set_rng (new_storage, self->rng, &error);
    _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_set_file");
  }

  _ncm_mset_catalog_copy_rows (mcat, self->storage, new_storage);

  ncm_mset_catalog_storage_clear (&self->storage);
  self->storage = new_storage;
  self->burnin  = 0;
}

/*
 * Replaces the storage of an empty catalog by a new one of the same kind holding desc,
 * which becomes the catalog's description, and the catalog's RNG. A FITS file is
 * recreated. A read-only catalog aborts.
 */
static void
_ncm_mset_catalog_rebuild_empty (NcmMSetCatalog *mcat, NcmMSetCatalogDesc *desc, const gchar *func)
{
  NcmMSetCatalogPrivate *self             = ncm_mset_catalog_get_instance_private (mcat);
  NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_peek_fits_storage (mcat);
  NcmMSetCatalogStorage *new_storage;
  GError *error = NULL;

  g_assert (ncm_mset_catalog_is_empty (mcat));

  if (storage_fits != NULL)
  {
    gchar *filename = g_strdup (ncm_mset_catalog_storage_fits_peek_filename (storage_fits));

    if (ncm_mset_catalog_storage_fits_get_read_only (storage_fits))
      g_error ("%s: cannot change a read-only catalog (`%s').", func, filename);

    _ncm_mset_catalog_end_recording (mcat);
    ncm_mset_catalog_storage_clear (&self->storage);

    if (g_file_test (filename, G_FILE_TEST_EXISTS))
      g_unlink (filename);

    new_storage = _ncm_mset_catalog_new_fits_storage (mcat, filename, FALSE);
    g_free (filename);
  }
  else
  {
    _ncm_mset_catalog_end_recording (mcat);
    ncm_mset_catalog_storage_clear (&self->storage);
    new_storage = _ncm_mset_catalog_new_memory_storage (mcat);
  }

  ncm_mset_catalog_storage_set_desc (new_storage, desc, &error);
  _ncm_mset_catalog_check_error (&error, func);

  if (self->rng != NULL)
  {
    ncm_mset_catalog_storage_set_rng (new_storage, self->rng, &error);
    _ncm_mset_catalog_check_error (&error, func);
  }

  self->storage = new_storage;

  if (desc != self->desc)
    _ncm_mset_catalog_take_desc (mcat, desc);

  self->cur_id = -1;
}

static void
_ncm_mset_catalog_set_add_val_name_array (NcmMSetCatalog *mcat, gchar **names)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (names != NULL)
  {
    const guint len = g_strv_length (names);
    guint i;

    g_ptr_array_set_size (self->add_vals_names, len);

    for (i = 0; i < len; i++)
    {
      g_clear_pointer (&g_ptr_array_index (self->add_vals_names, i), g_free);
      g_ptr_array_index (self->add_vals_names, i) = g_strdup (names[i]);
    }
  }
  else
  {
    g_ptr_array_set_size (self->add_vals_names, 0);
  }
}

static void
_ncm_mset_catalog_set_add_val_symbol_array (NcmMSetCatalog *mcat, gchar **symbols)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (symbols != NULL)
  {
    const guint len = g_strv_length (symbols);
    guint i;

    g_ptr_array_set_size (self->add_vals_symbs, len);

    for (i = 0; i < len; i++)
    {
      g_clear_pointer (&g_ptr_array_index (self->add_vals_symbs, i), g_free);
      g_ptr_array_index (self->add_vals_symbs, i) = g_strdup (symbols[i]);
    }
  }
  else
  {
    g_ptr_array_set_size (self->add_vals_symbs, 0);
  }
}

/**
 * ncm_mset_catalog_set_file:
 * @mcat: a #NcmMSetCatalog
 * @filename: (nullable): a FITS file name, or %NULL
 *
 * Puts @mcat on a writable #NcmMSetCatalogStorageFits of @filename. A catalog with rows
 * writes them to @filename, which must then not hold rows; an empty catalog attaches to
 * @filename, created when it does not exist, and takes its rows (with the burn-in of
 * ncm_mset_catalog_set_burnin()) and run metadata. The description of @filename must be
 * compatible with the catalog's (ncm_mset_catalog_desc_is_compatible()), else it
 * aborts. With %NULL the catalog moves to a #NcmMSetCatalogStorageMemory holding its
 * rows; the file is left as it was. The open recording is closed first.
 *
 */
void
ncm_mset_catalog_set_file (NcmMSetCatalog *mcat, const gchar *filename)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const gchar *cur_filename;

  if (!self->constructed)
  {
    if (self->file != NULL)
      g_error ("ncm_mset_catalog_set_file: file already set before construction.");

    self->file = g_strdup (filename);

    return;
  }

  cur_filename = ncm_mset_catalog_peek_filename (mcat);

  if ((cur_filename != NULL) && (filename != NULL) && (strcmp (cur_filename, filename) == 0))
    return;

  if (filename == NULL)
  {
    if (cur_filename != NULL)
      _ncm_mset_catalog_switch_storage (mcat, _ncm_mset_catalog_new_memory_storage (mcat));

    return;
  }

  {
    NcmMSetCatalogStorage *new_storage = _ncm_mset_catalog_new_fits_storage (mcat, filename, FALSE);
    GError *error                      = NULL;

    if (!ncm_mset_catalog_is_empty (mcat))
    {
      const guint file_nrows = ncm_mset_catalog_storage_nrows (new_storage, &error);

      _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_set_file");

      if (file_nrows > 0)
        g_error ("ncm_mset_catalog_set_file: the catalog has %u row(s) and file `%s' has %u; "
                 "a catalog with rows is written to a file without rows.",
                 ncm_mset_catalog_len (mcat), filename, file_nrows);

      _ncm_mset_catalog_switch_storage (mcat, new_storage);
    }
    else
    {
      _ncm_mset_catalog_end_recording (mcat);

      ncm_mset_catalog_storage_set_desc (new_storage, self->desc, &error);
      _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_set_file");

      ncm_mset_catalog_storage_clear (&self->storage);
      self->storage = new_storage;

      ncm_mset_catalog_reset_stats (mcat);
      _ncm_mset_catalog_load_storage (mcat);
    }
  }
}

/**
 * ncm_mset_catalog_set_sync_mode:
 * @mcat: a #NcmMSetCatalog
 * @smode: sync mode
 *
 * Sets #NcmMSetCatalog:smode, the sync mode of the storage.
 *
 */
void
ncm_mset_catalog_set_sync_mode (NcmMSetCatalog *mcat, NcmMSetCatalogSync smode)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->smode = smode;

  if (self->storage != NULL)
    g_object_set (self->storage, "smode", smode, NULL);
}

/**
 * ncm_mset_catalog_set_sync_interval:
 * @mcat: a #NcmMSetCatalog
 * @interval: minimum time between syncs, in seconds
 *
 * Sets the minimum time between the syncs of ncm_mset_catalog_timed_sync() and of the
 * storage in %NCM_MSET_CATALOG_SYNC_TIMED mode.
 *
 */
void
ncm_mset_catalog_set_sync_interval (NcmMSetCatalog *mcat, gdouble interval)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->sync_interval = interval;

  if (self->storage != NULL)
    g_object_set (self->storage, "sync-interval", interval, NULL);
}

/**
 * ncm_mset_catalog_set_read_block:
 * @mcat: a #NcmMSetCatalog
 * @read_block: number of rows, at least one
 *
 * Sets #NcmMSetCatalog:read-block; nothing changes for a catalog in memory.
 *
 */
void
ncm_mset_catalog_set_read_block (NcmMSetCatalog *mcat, guint read_block)
{
  NcmMSetCatalogPrivate *self             = ncm_mset_catalog_get_instance_private (mcat);
  NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_peek_fits_storage (mcat);

  g_assert_cmpuint (read_block, >, 0);

  self->read_block = read_block;

  if (storage_fits != NULL)
    g_object_set (storage_fits, "read-block", read_block, NULL);
}

/**
 * ncm_mset_catalog_get_read_block:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets #NcmMSetCatalog:read-block.
 *
 * Returns: the number of rows read from the file at once.
 */
guint
ncm_mset_catalog_get_read_block (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->read_block;
}

/**
 * ncm_mset_catalog_set_run_type:
 * @mcat: a #NcmMSetCatalog
 * @rtype_str: description of the run
 *
 * Sets #NcmMSetCatalog:run-type-string, the run type of the catalog's description
 * (ncm_mset_catalog_desc_get_run_type()). On an empty catalog a new value rebuilds the
 * storage with the new description (a FITS file is recreated); on a non-empty one it
 * aborts.
 *
 */
void
ncm_mset_catalog_set_run_type (NcmMSetCatalog *mcat, const gchar *rtype_str)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert (rtype_str != NULL);

  if ((self->rtype_str != NULL) && (strcmp (self->rtype_str, rtype_str) == 0))
    return;

  if (!self->constructed)
  {
    g_clear_pointer (&self->rtype_str, g_free);
    self->rtype_str = g_strdup (rtype_str);

    return;
  }

  if (!ncm_mset_catalog_is_empty (mcat))
    g_error ("ncm_mset_catalog_set_run_type: cannot change run type string in a non-empty catalog, actual: `%s' new: `%s'.",
             self->rtype_str, rtype_str);

  {
    NcmMSetCatalogDesc *desc = _ncm_mset_catalog_new_desc (mcat, self->nchains, self->m2lnp_var, rtype_str);

    _ncm_mset_catalog_rebuild_empty (mcat, desc, "ncm_mset_catalog_set_run_type");
    ncm_mset_catalog_desc_free (desc);
  }
}

/**
 * ncm_mset_catalog_set_rng:
 * @mcat: a #NcmMSetCatalog
 * @rng: a #NcmRNG
 *
 * Makes @rng the random number generator of @mcat, owned by its storage
 * (ncm_mset_catalog_storage_set_rng()). On a non-empty catalog it aborts: the catalog
 * continues from the generator state saved with its rows (ncm_mset_catalog_peek_rng()),
 * and a new generator would repeat rows already in it. On an empty catalog that has a
 * generator the storage is rebuilt with @rng (a FITS file is recreated).
 */
void
ncm_mset_catalog_set_rng (NcmMSetCatalog *mcat, NcmRNG *rng)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;

  if (!ncm_mset_catalog_is_empty (mcat))
    g_error ("ncm_mset_catalog_set_rng: refusing to set RNG in a non-empty catalog (%u rows): "
             "this would discard the saved RNG state and repeat rows already computed. "
             "Do not pass an explicit RNG/seed when resuming; the catalog's saved state is used.",
             ncm_mset_catalog_len (mcat));

  if (self->rng == rng)
    return;

  if (!self->constructed)
  {
    ncm_rng_clear (&self->rng);
    self->rng = ncm_rng_ref (rng);

    return;
  }

  if (self->rng != NULL)
  {
    ncm_rng_clear (&self->rng);
    self->rng = ncm_rng_ref (rng);
    _ncm_mset_catalog_rebuild_empty (mcat, self->desc, "ncm_mset_catalog_set_rng");

    return;
  }

  self->rng = ncm_rng_ref (rng);
  ncm_mset_catalog_storage_set_rng (self->storage, rng, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_set_rng");
}

/**
 * ncm_mset_catalog_sync:
 * @mcat: a #NcmMSetCatalog
 * @check: unused
 *
 * Writes the rows not yet in the file, and the generator state, through
 * ncm_mset_catalog_storage_sync(); nothing is done outside a recording (no row was added
 * since the catalog was created, reset or given a file).
 *
 */
void
ncm_mset_catalog_sync (NcmMSetCatalog *mcat, gboolean check)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;

  if ((self->storage == NULL) || !ncm_mset_catalog_storage_is_recording (self->storage))
    return;

  ncm_mset_catalog_storage_sync (self->storage, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_sync");
}

/**
 * ncm_mset_catalog_timed_sync:
 * @mcat: a #NcmMSetCatalog
 * @check: unused
 *
 * Calls ncm_mset_catalog_sync() when the time since the last sync exceeds the interval
 * of ncm_mset_catalog_set_sync_interval().
 *
 */
void
ncm_mset_catalog_timed_sync (NcmMSetCatalog *mcat, gboolean check)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (g_timer_elapsed (self->sync_timer, NULL) > self->sync_interval)
  {
    g_timer_start (self->sync_timer);
    ncm_mset_catalog_sync (mcat, check);
  }
}

/**
 * ncm_mset_catalog_reset_stats:
 * @mcat: a #NcmMSetCatalog
 *
 * Resets the statistics of @mcat (means, variances, chain and ensemble statistics,
 * autocorrelation, best fit, per-ensemble arrays) and keeps its rows.
 *
 */
void
ncm_mset_catalog_reset_stats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  ncm_stats_vec_reset (self->pstats, FALSE);
  ncm_stats_acorr_reset (self->acorr);

  if (self->nchains > 1)
  {
    guint i;

    for (i = 0; i < self->nchains; i++)
    {
      NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, i);

      ncm_stats_vec_reset (pstats, FALSE);
    }

    ncm_stats_vec_reset (self->mean_pstats, FALSE);
    ncm_stats_vec_reset (self->e_stats, FALSE);
    ncm_stats_vec_reset (self->e_mean_stats, FALSE);
    g_ptr_array_set_size (self->e_var_array, 0);
    g_array_set_size (self->accept_ratio, 0);
    self->naccepted = 0;
  }

  g_ptr_array_set_size (self->e_mean_rows, 0);
  _ncm_mset_catalog_reset_last_rows (self, self->nchains);

  ncm_vector_set_all (self->params_max, GSL_NEGINF);
  ncm_vector_set_all (self->params_min, GSL_POSINF);

  ncm_vector_clear (&self->bestfit_row);

  self->bestfit        = GSL_POSINF;
  self->post_lnnorm    = 0.0;
  self->post_lnnorm_up = FALSE;

  g_array_set_size (self->order_m2lnp, 0);
  self->order_sorted = FALSE;
}

/**
 * ncm_mset_catalog_reset:
 * @mcat: a #NcmMSetCatalog
 *
 * Erases the rows of @mcat, with all the statistics, through
 * ncm_mset_catalog_storage_erase(): the sampler segments and the Markovian id go with
 * them; the description, the file and the random number generator are kept. A read-only
 * catalog aborts.
 *
 */
void
ncm_mset_catalog_reset (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;

  _ncm_mset_catalog_end_recording (mcat);

  ncm_mset_catalog_storage_erase (self->storage, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_reset");

  ncm_mset_catalog_reset_stats (mcat);

  self->burnin       = 0;
  self->cur_id       = -1;
  self->markovian_id = 0; /* a new chain: nothing non-Markovian yet */
}

/**
 * ncm_mset_catalog_erase_data:
 * @mcat: a #NcmMSetCatalog
 *
 * Erases the rows of @mcat, as ncm_mset_catalog_reset().
 *
 */
void
ncm_mset_catalog_erase_data (NcmMSetCatalog *mcat)
{
  ncm_mset_catalog_reset (mcat);
}

/**
 * ncm_mset_catalog_set_m2lnp_var:
 * @mcat: a #NcmMSetCatalog
 * @p: $-2\ln(p)$ parameter index
 *
 * Sets @p as the $-2\ln(p)$ parameter index, part of the catalog's description. After
 * construction a new value rebuilds the storage of an empty catalog, as
 * ncm_mset_catalog_set_run_type(), and aborts on a non-empty one.
 *
 */
void
ncm_mset_catalog_set_m2lnp_var (NcmMSetCatalog *mcat, const gint p)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (p == self->m2lnp_var)
    return;

  if (!self->constructed)
  {
    self->m2lnp_var = p;

    return;
  }

  if (!ncm_mset_catalog_is_empty (mcat))
    g_error ("ncm_mset_catalog_set_m2lnp_var: cannot change the -2ln(p) column in a non-empty catalog, actual: %d new: %d.",
             self->m2lnp_var, p);

  {
    NcmMSetCatalogDesc *desc = _ncm_mset_catalog_new_desc (mcat, self->nchains, p, self->rtype_str);

    _ncm_mset_catalog_rebuild_empty (mcat, desc, "ncm_mset_catalog_set_m2lnp_var");
    ncm_mset_catalog_desc_free (desc);
  }
}

/**
 * ncm_mset_catalog_get_m2lnp_var:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the $-2\ln(p)$ parameter index.
 *
 * Returns: the $-2\ln(p)$ parameter index, -1 if unset.
 */
gint
ncm_mset_catalog_get_m2lnp_var (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->m2lnp_var;
}

/**
 * ncm_mset_catalog_peek_filename:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the file of the #NcmMSetCatalogStorageFits of @mcat.
 *
 * Returns: (transfer none) (nullable): the file name, %NULL for a catalog in memory.
 */
const gchar *
ncm_mset_catalog_peek_filename (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_peek_fits_storage (mcat);

  if (storage_fits != NULL)
    return ncm_mset_catalog_storage_fits_peek_filename (storage_fits);
  else
    return NULL;
}

/**
 * ncm_mset_catalog_get_rng:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer full) (allow-none): the random number generator of @mcat, or %NULL
 */
NcmRNG *
ncm_mset_catalog_get_rng (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->rng != NULL)
    return ncm_rng_ref (self->rng);
  else
    return NULL;
}

/**
 * ncm_mset_catalog_peek_rng:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer none) (allow-none): the random number generator of @mcat, or %NULL
 */
NcmRNG *
ncm_mset_catalog_peek_rng (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->rng;
}

/**
 * ncm_mset_catalog_is_empty:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: TRUE when the catalog is empty.
 */
gboolean
ncm_mset_catalog_is_empty (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return (self->cur_id < 0);
}

/**
 * ncm_mset_catalog_largest_error:
 * @mcat: a #NcmMSetCatalog
 *
 * This function calculates the largest proportional error of the parameters included,
 * $\text{lre} = \sigma_{\hat{p}}/(|\hat{p}|\sqrt{n_\mathrm{eff}})$, where
 * $n_\mathrm{eff}$ is the effective sample size of that parameter,
 * ncm_mset_catalog_get_ess(), $\hat{p}$ is its estimated mean and $\sigma_{\hat{p}}$ its
 * standard deviation. The effective sample size already carries both the autocorrelation
 * in time and the correlation between chains at a given iteration, so no assumption is
 * made about either.
 *
 * A mean compatible with zero makes the relative error about one: with at least ten
 * samples, a parameter whose $\text{lre}$ lies in $[1, 2)$ contributes its absolute
 * error $\sigma_{\hat{p}}/\sqrt{n_\mathrm{eff}}$ instead.
 *
 * Returns: the largest proportional error $\text{lre}$.
 */
gdouble
ncm_mset_catalog_largest_error (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint free_params_len       = ncm_mset_fparams_len (self->mset);
  const gdouble n             = ncm_stats_vec_get_weight (self->pstats);
  const gdouble fpi           = self->nadd_vals;
  const gdouble fpf           = free_params_len + self->nadd_vals;
  gdouble lerror              = 0.0;
  guint i;

  for (i = fpi; i < fpf; i++)
  {
    const gdouble mu       = ncm_stats_vec_get_mean (self->pstats, i);
    const gdouble sd       = ncm_stats_vec_get_sd (self->pstats, i);
    const gdouble ess      = ncm_mset_catalog_get_ess (mcat, i);
    const gdouble sqrt_ess = sqrt (GSL_MAX (ess, 1.0));
    gdouble lerror_i       = fabs (sd / (mu * sqrt_ess));

    if ((n >= 10) && ((guint) lerror_i == 1))
      lerror_i = fabs (sd / sqrt_ess);

    lerror = GSL_MAX (lerror, lerror_i);
  }

  return lerror;
}

/**
 * ncm_mset_catalog_len:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the number of rows of @mcat
 */
guint
ncm_mset_catalog_len (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_stats_vec_nitens (self->pstats);
}

/**
 * ncm_mset_catalog_max_time:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the number of rows divided by the number of chains: the number of iterations
 *   (ensembles)
 */
guint
ncm_mset_catalog_max_time (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->nchains > 1)
    return ncm_stats_vec_nitens (self->e_mean_stats);
  else
    return ncm_stats_vec_nitens (self->pstats);
}

/**
 * ncm_mset_catalog_nchains:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the number of chains of @mcat
 */
guint
ncm_mset_catalog_nchains (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->nchains;
}

/**
 * ncm_mset_catalog_nadd_vals:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the number of additional values per row, the weight included
 */
guint
ncm_mset_catalog_nadd_vals (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->nadd_vals;
}

/**
 * ncm_mset_catalog_weighted:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: whether the rows of @mcat are weighted
 */
gboolean
ncm_mset_catalog_weighted (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->weighted;
}

/**
 * ncm_mset_catalog_get_row_from_time:
 * @mcat: a #NcmMSetCatalog
 * @t: row id
 *
 * Returns: the index in @mcat of the row with id @t, which is @t
 */
guint
ncm_mset_catalog_get_row_from_time (NcmMSetCatalog *mcat, gint t)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert_cmpint (t, >=, 0);
  g_assert_cmpuint ((guint) t, <, ncm_stats_vec_nitens (self->pstats));

  return (guint) t;
}

/**
 * ncm_mset_catalog_get_cur_id:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the id of the last row added, its row index; -1 when empty
 */
gint
ncm_mset_catalog_get_cur_id (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->cur_id;
}

/**
 * ncm_mset_catalog_ncols:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: total number of columns in the catalog.
 */
guint
ncm_mset_catalog_ncols (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->ncols;
}

/**
 * ncm_mset_catalog_col_name:
 * @mcat: a #NcmMSetCatalog
 * @i: column index
 *
 * Returns: (transfer none): the name of the @i-th column.
 */
const gchar *
ncm_mset_catalog_col_name (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_mset_catalog_desc_col_name (self->desc, i);
}

/**
 * ncm_mset_catalog_col_full_name:
 * @mcat: a #NcmMSetCatalog
 * @i: column index
 *
 * Returns: (transfer none): the full name of the @i-th column.
 */
const gchar *
ncm_mset_catalog_col_full_name (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_mset_catalog_desc_col_full_name (self->desc, i);
}

/**
 * ncm_mset_catalog_col_symb:
 * @mcat: a #NcmMSetCatalog
 * @i: column index
 *
 * Returns: (transfer none): the symbol of the @i-th column.
 */
const gchar *
ncm_mset_catalog_col_symb (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_mset_catalog_desc_col_symb (self->desc, i);
}

/**
 * ncm_mset_catalog_col_by_name:
 * @mcat: a #NcmMSetCatalog
 * @name: column name
 * @col_index: (out): column index
 *
 * Finds the column given by @name: a free-parameter name, an additional value name or a
 * column number (ncm_mset_catalog_desc_col_by_name()).
 *
 * Returns: whether @name was found
 */
gboolean
ncm_mset_catalog_col_by_name (NcmMSetCatalog *mcat, const gchar *name, guint *col_index)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_mset_catalog_desc_col_by_name (self->desc, name, col_index);
}

/**
 * ncm_mset_catalog_set_burnin:
 * @mcat: a #NcmMSetCatalog
 * @burnin: number of leading rows to drop
 *
 * Sets the number of leading rows of the file that are not read; the ids of the rows
 * read start after them. It applies only when rows are loaded from a file, so it must be
 * set before that, and is meant for analyzing a catalog. Loading aborts when @burnin is
 * not a multiple of the number of chains, since the rows read would then be assigned to
 * the wrong chains. The burn-in is zero again once the catalog leaves the file
 * (ncm_mset_catalog_set_file(), ncm_mset_catalog_reset()).
 *
 */
void
ncm_mset_catalog_set_burnin (NcmMSetCatalog *mcat, glong burnin)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (_ncm_mset_catalog_peek_fits_storage (mcat) != NULL)
    g_error ("ncm_mset_catalog_set_burnin: cannot set burnin with an already loaded catalog");

  self->burnin = burnin;
}

/**
 * ncm_mset_catalog_get_burnin:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the burn-in size, see ncm_mset_catalog_set_burnin().
 *
 * Returns: the burn-in.
 */
glong
ncm_mset_catalog_get_burnin (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->burnin;
}

/**
 * ncm_mset_catalog_set_markovian_id:
 * @mcat: a #NcmMSetCatalog
 * @markovian_id: the id of the first row of the Markovian chain
 *
 * Sets #NcmMSetCatalog:markovian-id. The id can only move forward and cannot pass the row
 * after the current one; the producer calls it after the initial ensemble and after every
 * iteration that used a non-Markovian move. It is written to the storage inside the
 * recording (at once when one is open, else when the next row opens one).
 *
 */
void
ncm_mset_catalog_set_markovian_id (NcmMSetCatalog *mcat, gint markovian_id)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const gint cur              = ncm_mset_catalog_get_markovian_id (mcat);

  if (markovian_id < cur)
    g_error ("ncm_mset_catalog_set_markovian_id: the Markovian id cannot move backwards (%d < %d).", markovian_id, cur);

  if (markovian_id > self->cur_id + 1)
    g_error ("ncm_mset_catalog_set_markovian_id: id %d is beyond the row after the current one (%d).", markovian_id, self->cur_id + 1);

  self->markovian_id = markovian_id;

  if (self->storage != NULL)
    _ncm_mset_catalog_push_markovian_id (mcat);
}

/**
 * ncm_mset_catalog_get_markovian_id:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: #NcmMSetCatalog:markovian-id, never below 0.
 */
gint
ncm_mset_catalog_get_markovian_id (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return MAX (self->markovian_id, 0);
}

/**
 * ncm_mset_catalog_get_markovian_burnin:
 * @mcat: a #NcmMSetCatalog
 *
 * Number of complete iterations (ensemble steps) before the Markovian chain starts, i.e.
 * the smallest burn-in that leaves only Markovian rows; 0 when every row is Markovian.
 *
 * Returns: the burn-in in iterations.
 */
guint
ncm_mset_catalog_get_markovian_burnin (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const gint rows             = ncm_mset_catalog_get_markovian_id (mcat);

  return (rows + self->nchains - 1) / self->nchains;
}

/**
 * ncm_mset_catalog_set_tau_method:
 * @mcat: a #NcmMSetCatalog
 * @tau_method: a #NcmStatsAcorrMethod
 *
 * Sets the estimator the autocorrelation time is computed with.
 *
 */
void
ncm_mset_catalog_set_tau_method (NcmMSetCatalog *mcat, NcmStatsAcorrMethod tau_method)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->tau_method = tau_method;

  /* tau-method is a G_PARAM_CONSTRUCT property, so it is set before constructed() builds
   * the accumulator; the value is kept here and applied there. */
  if (self->acorr != NULL)
    ncm_stats_acorr_set_method (self->acorr, tau_method);
}

/**
 * ncm_mset_catalog_get_tau_method:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the estimator the autocorrelation time is computed with.
 */
NcmStatsAcorrMethod
ncm_mset_catalog_get_tau_method (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->tau_method;
}

/**
 * ncm_mset_catalog_set_post_lnnorm_method:
 * @mcat: a #NcmMSetCatalog
 * @method: a #NcmMSetCatalogPostNormMethod
 *
 * Sets the estimator of ncm_mset_catalog_get_post_lnnorm(). A different estimator
 * discards the kept estimate.
 *
 */
void
ncm_mset_catalog_set_post_lnnorm_method (NcmMSetCatalog *mcat, NcmMSetCatalogPostNormMethod method)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (method != self->post_lnnorm_method)
  {
    self->post_lnnorm_method = method;
    self->post_lnnorm_up     = FALSE;
  }
}

/**
 * ncm_mset_catalog_get_post_lnnorm_method:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the estimator of ncm_mset_catalog_get_post_lnnorm().
 */
NcmMSetCatalogPostNormMethod
ncm_mset_catalog_get_post_lnnorm_method (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->post_lnnorm_method;
}

/*
 * The autocorrelation accumulator forms unweighted lagged sums. A weighted catalog is
 * therefore not fed (see _ncm_mset_catalog_post_update_stats) and every quantity derived
 * from the accumulator refuses rather than returning a number built from the wrong sums.
 */
static void
_ncm_mset_catalog_check_unweighted (NcmMSetCatalogPrivate *self, const gchar *func)
{
  if (self->weighted)
    g_error ("%s: the autocorrelation accumulator does not support weighted catalogs: it "
             "forms unweighted lagged sums, so the autocorrelation time and the effective "
             "sample size built from them would not describe this catalog.", func);
}

/**
 * ncm_mset_catalog_peek_acorr:
 * @mcat: a #NcmMSetCatalog
 *
 * The autocorrelation accumulator, updated as rows are added. For a catalog with more
 * than one chain it is fed the ensemble mean of each complete iteration, which is the
 * series whose spectral density at zero gives the error of the reported mean; for a
 * single chain it is fed every row.
 *
 * Returns: (transfer none): the #NcmStatsAcorr.
 */
NcmStatsAcorr *
ncm_mset_catalog_peek_acorr (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_check_unweighted (self, G_STRFUNC);

  return self->acorr;
}

/*
 * Feeds the statistics with the row x, the storage's vector of the row just added at
 * index len (the rows of the catalog before it): the extrema of every column, the
 * running statistics of all rows and of its chain, the acceptance count (x against the
 * row of the same chain one iteration before, a reference kept in last_rows), the
 * ensemble statistics when x completes an ensemble, the best fit and the -2ln(p)
 * order. The current id moves to x.
 */
static void
_ncm_mset_catalog_post_update_stats (NcmMSetCatalog *mcat, NcmVector *x)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint len             = self->ncols;
  const guint index           = ncm_stats_vec_nitens (self->pstats);
  const guint chain           = index % self->nchains;
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble p_i         = ncm_vector_get (x, i);
    const gdouble cur_max_p_i = ncm_vector_get (self->params_max, i);
    const gdouble cur_min_p_i = ncm_vector_get (self->params_min, i);

    ncm_vector_set (self->params_max, i, GSL_MAX (p_i, cur_max_p_i));
    ncm_vector_set (self->params_min, i, GSL_MIN (p_i, cur_min_p_i));
  }

  /* Counting the number of accepted steps */
  if (index >= self->nchains)
  {
    NcmVector *last_x = g_ptr_array_index (self->last_rows, chain);

    if (ncm_vector_cmp2 (last_x, x, 1.0e-14, 0.0) > 0)
      self->naccepted++;
  }

  g_clear_pointer (&g_ptr_array_index (self->last_rows, chain), ncm_vector_free);
  g_ptr_array_index (self->last_rows, chain) = ncm_vector_ref (x);

  if (self->weighted)
  {
    /* The accumulator is not fed here: it forms unweighted lagged sums, and everything
     * that reads it refuses for a weighted catalog. */
    if (self->nchains > 1)
    {
      guint chain_id      = (self->cur_id + 1) % self->nchains;
      NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, chain_id);

      ncm_stats_vec_append_weight (pstats,        x, ncm_vector_get (x, self->nadd_vals - 1), FALSE);
      ncm_stats_vec_append_weight (self->e_stats, x, ncm_vector_get (x, self->nadd_vals - 1), FALSE);
    }

    ncm_stats_vec_append_weight (self->pstats, x, ncm_vector_get (x, self->nadd_vals - 1), FALSE);
  }
  else
  {
    if (self->nchains > 1)
    {
      guint chain_id      = (self->cur_id + 1) % self->nchains;
      NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, chain_id);

      ncm_stats_vec_append (pstats,        x, FALSE);
      ncm_stats_vec_append (self->e_stats, x, FALSE);
    }
    else
    {
      ncm_stats_acorr_update (self->acorr, x);
    }

    ncm_stats_vec_append (self->pstats, x, FALSE);
  }

  self->cur_id++;

  if (self->nchains > 1)
  {
    if ((self->cur_id + 1) % self->nchains == 0)
    {
      NcmVector *e_mean = ncm_stats_vec_peek_mean (self->e_stats);
      const guint len   = ncm_vector_len (e_mean);
      NcmVector *e_var  = ncm_vector_new (len);
      guint i;

      for (i = 0; i < len; i++)
      {
        ncm_vector_set (e_var, i, ncm_stats_vec_get_var (self->e_stats, i));
      }

      g_assert_cmpuint (ncm_stats_vec_nitens (self->e_stats), ==, self->nchains);

      ncm_stats_vec_append (self->e_mean_stats, e_mean, FALSE);
      g_ptr_array_add (self->e_mean_rows, ncm_vector_dup (e_mean));
      ncm_stats_acorr_update (self->acorr, e_mean);
      g_ptr_array_add (self->e_var_array, e_var);

      if (ncm_stats_vec_nitens (self->pstats) > self->nchains)
      {
        const gdouble acc_ratio = self->naccepted / (gdouble) self->nchains;

        g_array_append_val (self->accept_ratio, acc_ratio);
      }

      ncm_stats_vec_reset (self->e_stats, FALSE);
      self->naccepted = 0;
    }
  }

  self->post_lnnorm_up = FALSE;

  if (self->m2lnp_var >= 0)
  {
    gdouble *m2lnp = ncm_vector_ptr (x, self->m2lnp_var);

    if (*m2lnp < self->bestfit)
    {
      ncm_vector_clear (&self->bestfit_row);
      self->bestfit_row = ncm_vector_ref (x);
      self->bestfit     = *m2lnp;

      ncm_mset_fparams_set_vector_offset (self->mset, x, self->nadd_vals);
    }

    g_array_append_val (self->order_m2lnp, *m2lnp);
    self->order_sorted = FALSE;
  }
}

/**
 * ncm_mset_catalog_add_from_mset:
 * @mcat: a #NcmMSetCatalog
 * @mset: a #NcmMSet
 * @...: additional values
 *
 * Adds a row with the additional values @... (the weight last, for a weighted catalog)
 * and the free parameters of @mset, which must have the free parameters of the catalog.
 *
 */
void
ncm_mset_catalog_add_from_mset (NcmMSetCatalog *mcat, NcmMSet *mset, ...)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmVector *row_i            = ncm_vector_new (self->ncols);
  va_list ap;
  guint i;

  g_assert_cmpuint (ncm_mset_fparams_len (mset) + self->nadd_vals, ==, self->ncols);

  va_start (ap, mset);

  for (i = 0; i < self->nadd_vals; i++)
  {
    gdouble val = va_arg (ap, gdouble);

    ncm_vector_set (row_i, i, val);
  }

  va_end (ap);

  ncm_mset_fparams_get_vector_offset (mset, row_i, self->nadd_vals);

  _ncm_mset_catalog_add_row (mcat, row_i);
  ncm_vector_free (row_i);
}

/**
 * ncm_mset_catalog_add_from_mset_array:
 * @mcat: a #NcmMSetCatalog
 * @mset: a #NcmMSet
 * @ax: (array) (element-type double): additional values array
 *
 * Adds a row with the additional values @ax (the weight last, for a weighted catalog)
 * and the free parameters of @mset, which must have the free parameters of the catalog.
 *
 */
void
ncm_mset_catalog_add_from_mset_array (NcmMSetCatalog *mcat, NcmMSet *mset, gdouble *ax)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmVector *row_i            = ncm_vector_new (self->ncols);
  guint i;

  g_assert_cmpuint (ncm_mset_fparams_len (mset) + self->nadd_vals, ==, self->ncols);

  for (i = 0; i < self->nadd_vals; i++)
    ncm_vector_set (row_i, i, ax[i]);

  ncm_mset_fparams_get_vector_offset (mset, row_i, self->nadd_vals);

  _ncm_mset_catalog_add_row (mcat, row_i);
  ncm_vector_free (row_i);
}

/**
 * ncm_mset_catalog_add_from_vector:
 * @mcat: a #NcmMSetCatalog
 * @vals: a whole row: the additional values followed by the free parameters
 *
 * Adds the row @vals to the storage (ncm_mset_catalog_storage_add_row()), opening a
 * recording when none is open, and updates the statistics.
 *
 */
void
ncm_mset_catalog_add_from_vector (NcmMSetCatalog *mcat, NcmVector *vals)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert_cmpuint (ncm_vector_len (vals), ==, self->ncols);

  _ncm_mset_catalog_add_row (mcat, vals);
}

/**
 * ncm_mset_catalog_add_from_vector_array:
 * @mcat: a #NcmMSetCatalog
 * @vals: the free parameters
 * @ax: (array) (element-type double): the additional values
 *
 * Adds a row with the additional values @ax and the free parameters @vals.
 *
 */
void
ncm_mset_catalog_add_from_vector_array (NcmMSetCatalog *mcat, NcmVector *vals, gdouble *ax)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmVector *row_i            = ncm_vector_new (self->ncols);
  guint i;

  g_assert_cmpuint (ncm_vector_len (vals) + self->nadd_vals, ==, self->ncols);

  for (i = 0; i < self->nadd_vals; i++)
    ncm_vector_set (row_i, i, ax[i]);

  ncm_vector_memcpy2 (row_i, vals, self->nadd_vals, 0, ncm_vector_len (vals));

  _ncm_mset_catalog_add_row (mcat, row_i);
  ncm_vector_free (row_i);
}

static gdouble
_fmeansd (gdouble v_i, guint i, gpointer user_data)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (user_data);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return sqrt (ncm_stats_vec_get_var (self->pstats, i) * ncm_vector_get (self->tau, i) / ncm_stats_vec_nitens (self->pstats));
}

static gdouble
_fsd (gdouble v_i, guint i, gpointer user_data)
{
  NcmStatsVec *pstats = NCM_STATS_VEC (user_data);

  return ncm_stats_vec_get_sd (pstats, i);
}

static gdouble
_fvar (gdouble v_i, guint i, gpointer user_data)
{
  NcmStatsVec *pstats = NCM_STATS_VEC (user_data);

  return ncm_stats_vec_get_var (pstats, i);
}

static gdouble
_ftau (gdouble v_i, guint i, gpointer user_data)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (user_data);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_vector_get (self->tau, i);
}

/**
 * ncm_mset_catalog_log_current_stats:
 * @mcat: a #NcmMSetCatalog
 *
 * Logs the current means and standard deviations of the catalog's parameters.
 *
 */
void
ncm_mset_catalog_log_current_stats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmVector *mean             = ncm_stats_vec_peek_mean (self->pstats);

  ncm_vector_log_vals (mean,      "# NcmMSetCatalog: Current mean:  ", "% -12.5g", TRUE);
  ncm_vector_log_vals_func (mean, "# NcmMSetCatalog: Current msd:   ", "% -12.5g", &_fmeansd, mcat);
  ncm_vector_log_vals_func (mean, "# NcmMSetCatalog: Current sd:    ", "% -12.5g", &_fsd, self->pstats);
  ncm_vector_log_vals_func (mean, "# NcmMSetCatalog: Current var:   ", "% -12.5g", &_fvar, self->pstats);
  ncm_vector_log_vals_func (mean, "# NcmMSetCatalog: Current tau:   ", "% -12.5g", &_ftau, mcat);
}

/**
 * ncm_mset_catalog_get_mset:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer full): the #NcmMSet of @mcat
 */
NcmMSet *
ncm_mset_catalog_get_mset (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_mset_ref (self->mset);
}

/**
 * ncm_mset_catalog_peek_mset:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer none): the #NcmMSet of @mcat
 */
NcmMSet *
ncm_mset_catalog_peek_mset (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->mset;
}

/**
 * ncm_mset_catalog_set_functions_array:
 * @mcat: a #NcmMSetCatalog
 * @functions: (nullable): a #NcmObjArray of #NcmMSetFunc computing @mcat's extra columns
 *
 * Sets the array of functions that compute the additional values of @mcat. It is kept
 * in memory only: the description records the names and symbols of the columns, and
 * the functions belong to the experiment that produced the catalog.
 *
 */
void
ncm_mset_catalog_set_functions_array (NcmMSetCatalog *mcat, NcmObjArray *functions)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  ncm_obj_array_clear (&self->functions);

  if (functions != NULL)
    self->functions = ncm_obj_array_ref (functions);
}

/**
 * ncm_mset_catalog_peek_functions_array:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the array set by ncm_mset_catalog_set_functions_array(); a catalog loaded from a
 * file has none.
 *
 * Returns: (transfer none) (nullable): the #NcmObjArray of #NcmMSetFunc, or %NULL if none.
 */
NcmObjArray *
ncm_mset_catalog_peek_functions_array (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->functions;
}

/**
 * ncm_mset_catalog_get_run_type:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the run type string from @mcat.
 *
 * Returns: (transfer none): the run type string.
 */
const gchar *
ncm_mset_catalog_get_run_type (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->rtype_str;
}

/*
 * Sets one of the sampler strings. A change while the sampler's recording is open closes
 * it, so that the next row opens a segment with the new value.
 */
static void
_ncm_mset_catalog_set_record (NcmMSetCatalog *mcat, gchar **field, const gchar *value)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if ((value == NULL) && (*field == NULL))
    return;

  if ((value != NULL) && (*field != NULL) && (strcmp (*field, value) == 0))
    return;

  g_clear_pointer (field, g_free);
  *field = g_strdup (value);

  if ((self->storage != NULL) && !self->init_recording)
    _ncm_mset_catalog_end_recording (mcat);
}

/**
 * ncm_mset_catalog_set_sampler:
 * @mcat: a #NcmMSetCatalog
 * @sampler: (nullable): a description of the sampler about to fill the rows
 *
 * Sets #NcmMSetCatalog:sampler, the name of the sampler segment the next rows are
 * recorded in. A change closes the open recording, so the next row starts a segment with
 * the new name; a sampler that adds no rows is not recorded. The run type names the
 * segment when no sampler is set.
 *
 */
void
ncm_mset_catalog_set_sampler (NcmMSetCatalog *mcat, const gchar *sampler)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_set_record (mcat, &self->sampler_str, sampler);
}

/**
 * ncm_mset_catalog_set_sampler_options:
 * @mcat: a #NcmMSetCatalog
 * @options: (nullable): the sampler's tunable settings
 *
 * Sets #NcmMSetCatalog:sampler-options, the options of the sampler segment the next
 * rows are recorded in, as ncm_mset_catalog_set_sampler().
 *
 */
void
ncm_mset_catalog_set_sampler_options (NcmMSetCatalog *mcat, const gchar *options)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_set_record (mcat, &self->sampler_opts_str, options);
}

/**
 * ncm_mset_catalog_set_initial_sampler:
 * @mcat: a #NcmMSetCatalog
 * @sampler: (nullable): a description of the sampler that draws the initial points
 *
 * Sets #NcmMSetCatalog:initial-sampler. On an empty storage the first row opens a
 * segment with this name, closed when the first ensemble is complete; after that the
 * value is a record only.
 *
 */
void
ncm_mset_catalog_set_initial_sampler (NcmMSetCatalog *mcat, const gchar *sampler)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if ((sampler == NULL) && (self->init_sampler_str == NULL))
    return;

  if ((sampler != NULL) && (self->init_sampler_str != NULL) && (strcmp (self->init_sampler_str, sampler) == 0))
    return;

  g_clear_pointer (&self->init_sampler_str, g_free);
  self->init_sampler_str = g_strdup (sampler);
}

/**
 * ncm_mset_catalog_get_sampler:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer none) (nullable): the sampler description, see #NcmMSetCatalog:sampler.
 */
const gchar *
ncm_mset_catalog_get_sampler (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->sampler_str;
}

/**
 * ncm_mset_catalog_get_sampler_options:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer none) (nullable): the sampler's tunable settings, see #NcmMSetCatalog:sampler-options.
 */
const gchar *
ncm_mset_catalog_get_sampler_options (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->sampler_opts_str;
}

/**
 * ncm_mset_catalog_get_initial_sampler:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer none) (nullable): the initial points sampler description, see #NcmMSetCatalog:initial-sampler.
 */
const gchar *
ncm_mset_catalog_get_initial_sampler (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->init_sampler_str;
}

/**
 * ncm_mset_catalog_log_current_chain_stats:
 * @mcat: a #NcmMSetCatalog
 *
 * Logs the current means and standard deviations of the catalog's parameters
 * for each chain.
 *
 */
void
ncm_mset_catalog_log_current_chain_stats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint fparams_len     = ncm_mset_fparams_len (self->mset);

  if (self->nchains > 1)
  {
    NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
    const gdouble shrink_factor = ncm_mset_catalog_get_shrink_factor (mcat);
    guint i;

    g_message ("# NcmMSetCatalog: Current skfac:");

    for (i = 0; i < fparams_len + self->nadd_vals; i++)
    {
      const gdouble shrink_factor_i = ncm_mset_catalog_get_param_shrink_factor (mcat, i);

      g_message (" % -12.5g", shrink_factor_i);
    }

    g_message ("\n");

    g_message ("# NcmMSetCatalog: Maximal Shrink factor = %-22.6g\n", shrink_factor);

    ncm_mset_catalog_calc_const_break (mcat, self->m2lnp_var, NCM_FIT_RUN_MSGS_SIMPLE);
  }
}

/**
 * ncm_mset_catalog_peek_pstats:
 * @mcat: a #NcmMSetCatalog
 *
 * Peeks the online statistics of all rows (mean, variance, covariance). It keeps no rows;
 * they are read with ncm_mset_catalog_peek_row().
 *
 * Returns: (transfer none): the #NcmStatsVec object.
 */
NcmStatsVec *
ncm_mset_catalog_peek_pstats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->pstats;
}

/**
 * ncm_mset_catalog_peek_e_mean_stats:
 * @mcat: a #NcmMSetCatalog
 *
 * Peeks the online statistics of the ensemble means. It keeps no rows; the series is read
 * with ncm_mset_catalog_peek_e_mean_t() or ncm_mset_catalog_get_diag_rows().
 *
 * Returns: (transfer none): the #NcmStatsVec object.
 */
NcmStatsVec *
ncm_mset_catalog_peek_e_mean_stats (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->e_mean_stats;
}

/**
 * ncm_mset_catalog_peek_chain_pstats:
 * @mcat: a #NcmMSetCatalog
 * @i: chain index
 *
 * Peeks the online statistics of chain @i. It keeps no rows.
 *
 * Returns: (transfer none): the #NcmStatsVec object.
 */
NcmStatsVec *
ncm_mset_catalog_peek_chain_pstats (NcmMSetCatalog *mcat, const guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return g_ptr_array_index (self->chain_pstats, i);
}

static GPtrArray *_ncm_mset_catalog_get_diag_series (NcmMSetCatalog *mcat);
static GPtrArray *_ncm_mset_catalog_collect_rows (NcmMSetCatalog *mcat, const guint first, const guint step);

/**
 * ncm_mset_catalog_get_diag_rows:
 * @mcat: a #NcmMSetCatalog
 *
 * Creates a #NcmStatsRows given the series the burn-in diagnostics of @mcat read: the
 * rows when @mcat has one chain, the ensemble means when it has several. It holds a
 * reference to each row present when it is created, so it stays valid when @mcat adds
 * rows or is reset; rows added later are not part of it.
 *
 * Returns: (transfer full): a new #NcmStatsRows.
 */
NcmStatsRows *
ncm_mset_catalog_get_diag_rows (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmStatsRows *rows          = ncm_stats_rows_new (self->ncols);
  GPtrArray *sample           = _ncm_mset_catalog_get_diag_series (mcat);

  ncm_stats_rows_set_rows (rows, sample);
  g_ptr_array_unref (sample);

  return rows;
}

/*
 * A new array holding a reference to each element of the series the burn-in diagnostics
 * read: the rows when the catalog has one chain, the ensemble means when it has several.
 */
static GPtrArray *
_ncm_mset_catalog_get_diag_series (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->nchains == 1)
  {
    return _ncm_mset_catalog_collect_rows (mcat, 0, 1);
  }
  else
  {
    GPtrArray *dup = g_ptr_array_new_full (self->e_mean_rows->len, (GDestroyNotify) & ncm_vector_free);
    guint i;

    for (i = 0; i < self->e_mean_rows->len; i++)
      g_ptr_array_add (dup, ncm_vector_ref (g_ptr_array_index (self->e_mean_rows, i)));

    return dup;
  }
}

/*
 * A new array holding a reference to the rows first, first + step, first + 2 step, ...,
 * up to the last row of the catalog; empty when first is past the last row.
 */
static GPtrArray *
_ncm_mset_catalog_collect_rows (NcmMSetCatalog *mcat, const guint first, const guint step)
{
  const guint len   = ncm_mset_catalog_len (mcat);
  const guint nrows = (len > first) ? (len - first + step - 1) / step : 0;
  GPtrArray *rows   = g_ptr_array_new_full (nrows, (GDestroyNotify) & ncm_vector_free);
  guint t;

  for (t = 0; t < nrows; t++)
    g_ptr_array_add (rows, ncm_vector_ref (ncm_mset_catalog_peek_row (mcat, first + t * step)));

  return rows;
}

/**
 * ncm_mset_catalog_peek_accept_ratio_array:
 * @mcat: a #NcmMSetCatalog
 *
 * Peeks the acceptance ratio array.
 *
 * Returns: (transfer none) (element-type gdouble): the acceptance ratio array.
 */
GArray *
ncm_mset_catalog_peek_accept_ratio_array (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->accept_ratio;
}

/**
 * ncm_mset_catalog_peek_row:
 * @mcat: a #NcmMSetCatalog
 * @i: the row index
 *
 * Gets the @i-th row from the storage (ncm_mset_catalog_storage_peek_row()): the row
 * @i + burn-in of the storage. The vector is readable until the next call on @mcat that
 * reads another row, adds a row, syncs, resets, changes the file, or until @mcat is
 * disposed; hold a reference (ncm_vector_ref()) to keep it longer. Not safe while
 * another thread reads rows of @mcat; use ncm_mset_catalog_get_row() there.
 *
 * Returns: (transfer none): the row with index @i or NULL if not available.
 */
NcmVector *
ncm_mset_catalog_peek_row (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;
  NcmVector *row;

  if (i >= ncm_mset_catalog_len (mcat))
    return NULL;

  row = ncm_mset_catalog_storage_peek_row (self->storage, i + self->burnin, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_peek_row");

  return row;
}

/**
 * ncm_mset_catalog_get_row:
 * @mcat: a #NcmMSetCatalog
 * @i: the row index
 *
 * Gets the @i-th row as ncm_mset_catalog_peek_row() does, holding a lock and returning a
 * new reference, so that several threads can read rows of @mcat at the same time. No
 * thread may add rows or change the file of @mcat meanwhile.
 *
 * Returns: (transfer full) (nullable): the row with index @i or NULL if not available.
 */
NcmVector *
ncm_mset_catalog_get_row (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GError *error               = NULL;
  NcmVector *row;

  if (i >= ncm_mset_catalog_len (mcat))
    return NULL;

  g_mutex_lock (&self->row_lock);
  row = ncm_mset_catalog_storage_get_row (self->storage, i + self->burnin, &error);
  g_mutex_unlock (&self->row_lock);

  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_get_row");

  return row;
}

/**
 * ncm_mset_catalog_peek_current_row:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the last added row.
 *
 * Returns: (transfer none): last added row or NULL if empty.
 */
NcmVector *
ncm_mset_catalog_peek_current_row (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (ncm_stats_vec_nitens (self->pstats) == 0)
    return NULL;
  else
    return ncm_mset_catalog_peek_row (mcat, ncm_stats_vec_nitens (self->pstats) - 1);
}

/**
 * ncm_mset_catalog_peek_current_e_mean:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the last ensemble mean.
 *
 * Returns: (transfer none): last added ensemble mean or NULL if empty.
 */
NcmVector *
ncm_mset_catalog_peek_current_e_mean (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (ncm_stats_vec_nitens (self->e_mean_stats) > 0)
    return g_ptr_array_index (self->e_mean_rows, ncm_stats_vec_nitens (self->e_mean_stats) - 1);
  else
    return NULL;
}

/**
 * ncm_mset_catalog_peek_current_e_var:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the last ensemble variance.
 *
 * Returns: (transfer none): last added ensemble variance or NULL if empty.
 */
NcmVector *
ncm_mset_catalog_peek_current_e_var (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->e_var_array->len > 0)
    return g_ptr_array_index (self->e_var_array, self->e_var_array->len - 1);
  else
    return NULL;
}

/**
 * ncm_mset_catalog_peek_e_mean_t:
 * @mcat: a #NcmMSetCatalog
 * @t: time
 *
 * Gets the mean of the @t-th ensemble.
 *
 * Returns: (transfer none): @t-th ensemble mean.
 */
NcmVector *
ncm_mset_catalog_peek_e_mean_t (NcmMSetCatalog *mcat, guint t)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->nchains > 1)
  {
    g_assert_cmpuint (t, <, ncm_stats_vec_nitens (self->e_mean_stats));

    return g_ptr_array_index (self->e_mean_rows, t);
  }
  else
  {
    return NULL;
  }
}

/**
 * ncm_mset_catalog_peek_e_var_t:
 * @mcat: a #NcmMSetCatalog
 * @t: MCMC time
 *
 * Gets the variance of the @t-th ensemble.
 *
 * Returns: (transfer none): @t-th ensemble variance.
 */
NcmVector *
ncm_mset_catalog_peek_e_var_t (NcmMSetCatalog *mcat, guint t)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert_cmpuint (t, <, self->e_var_array->len);

  return g_ptr_array_index (self->e_var_array, t);
}

#define NCM_MSET_CATALOG_RESCALE_COV (0.80)

static gdouble _ncm_mset_catalog_get_post_lnnorm_sum (NcmMSetCatalog *mcat, NcmVector *mean, NcmMatrix *cov, const gdouble lnNorma, const gdouble R2_cut, gdouble *lnnorm_sd);

static gdouble
_ncm_mset_catalog_get_post_lnnorm_ellipsoid (NcmMSetCatalog *mcat, gdouble *post_lnnorm_sd)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint fparams_len     = ncm_mset_fparams_len (self->mset);
  NcmMatrix *cov              = NULL;
  NcmVector *mean             = NULL;
  gdouble level               = 0.50;
  gdouble R2_cut              = gsl_cdf_chisq_Pinv (level, fparams_len);
  gdouble lnNorma             = 0.0;
  gdouble R_max               = 1.0e300;
  gdouble R2_max, post_lnnorm;
  gint ret;
  guint i;

  ncm_mset_catalog_get_covar (mcat, &cov);
  ncm_mset_catalog_get_mean (mcat, &mean);
  ncm_matrix_scale (cov, gsl_pow_2 (NCM_MSET_CATALOG_RESCALE_COV));

  for (i = 0; i < fparams_len; i++)
  {
    const gdouble mu_i    = ncm_vector_get (mean, i);
    const gdouble sigma_i = sqrt (ncm_matrix_get (cov, i, i));
    const gdouble l_i     = ncm_mset_fparam_get_lower_bound (self->mset, i);
    const gdouble u_i     = ncm_mset_fparam_get_upper_bound (self->mset, i);

    R_max = MIN (R_max, +(mu_i - l_i) / sigma_i);
    R_max = MIN (R_max, -(mu_i - u_i) / sigma_i);
  }

  R2_max = R_max * R_max;

  if (R2_cut > R2_max)
  {
    R2_cut = R2_max;
    level  = gsl_cdf_chisq_P (R2_cut, fparams_len);
  }

  ret = ncm_matrix_cholesky_decomp (cov, 'U');

  if (ret != 0)
  {
    g_warning ("ncm_mset_catalog_get_post_lnnorm[ncm_matrix_cholesky_decomp]: %d. Non-positive definite covariance, more points are necessary.", ret);

    ncm_vector_clear (&mean);
    ncm_matrix_clear (&cov);
    post_lnnorm_sd[0] = GSL_NAN;

    return 0.0;
  }

  lnNorma     = 0.5 * (fparams_len * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (cov)) + log (level);
  post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_sum (mcat, mean, cov, lnNorma, R2_cut, post_lnnorm_sd);

  ncm_vector_clear (&mean);
  ncm_matrix_clear (&cov);

  return post_lnnorm;
}

static gdouble
_ncm_mset_catalog_get_post_lnnorm_sum (NcmMSetCatalog *mcat, NcmVector *mean, NcmMatrix *cov, const gdouble lnNorma, const gdouble R2_cut, gdouble *lnnorm_sd)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint fparams_len     = ncm_mset_fparams_len (self->mset);
  const guint cat_len         = ncm_mset_catalog_len (mcat);
  NcmVector *v                = ncm_vector_new (fparams_len);
  NcmStatsVec *slnnorm        = ncm_stats_vec_new (2, NCM_STATS_VEC_VAR, FALSE);
  gdouble s                   = 0.0;
  gdouble c                   = 0.0;
  guint slice_min             = 100;
  guint slice_div             = 100000;
  guint slice_size            = cat_len;
  guint slice_res             = 0;
  guint nslices               = 1;
  gdouble mean_lnnorm;
  gint ret;
  guint i, j, w;

  while ((cat_len / slice_div) < slice_min)
  {
    slice_div = slice_div / 10;

    if (slice_div <= 1)
    {
      slice_size = cat_len;
      slice_res  = 0;
      nslices    = 1;
      g_warning ("_ncm_mset_catalog_get_post_lnnorm_sum: catalog too small to estimate error on the posterior norm.");
      break;
    }
    else
    {
      slice_size = cat_len / slice_div;
      slice_res  = cat_len % slice_div;
      nslices    = slice_div;
    }
  }

  w = 0;

  for (j = 0; j < nslices; j++)
  {
    const guint slice_size1 = (j == 0) ? (slice_size + slice_res) : slice_size;

    s = 0.0;
    c = 0.0;

    for (i = 0; i < slice_size1; i++)
    {
      NcmVector *row_i      = ncm_mset_catalog_peek_row (mcat, i + w);
      const gdouble m2lnL_i = ncm_vector_get (row_i, self->m2lnp_var);
      gdouble m2lnp_i       = 0.0;
      gdouble e_i, t;

      ncm_vector_memcpy2 (v, row_i, 0, self->nadd_vals, fparams_len);
      ncm_vector_sub (v, mean);

      ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                            ncm_matrix_gsl (cov), ncm_vector_gsl (v));
      NCM_TEST_GSL_RESULT ("ncm_mset_catalog_get_post_lnnorm", ret);

      ret = gsl_blas_ddot (ncm_vector_gsl (v), ncm_vector_gsl (v), &m2lnp_i);
      NCM_TEST_GSL_RESULT ("ncm_mset_catalog_get_post_lnnorm", ret);

      /* A row outside the cut contributes zero. */
      if (m2lnp_i > R2_cut)
        continue;

      e_i = exp (0.5 * ((m2lnL_i - self->bestfit) - m2lnp_i));
      t   = s + e_i;
      c  += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
      s   = t;
    }

    ncm_stats_vec_set (slnnorm, 0, (s + c) / (slice_size1));
    ncm_stats_vec_set (slnnorm, 1, -(log ((s + c) / (slice_size1)) - lnNorma + 0.5 * self->bestfit));

    ncm_stats_vec_update_weight (slnnorm, slice_size1);

    w += slice_size1;

/*
 *
 *   printf ("# SL LNNORM: % 22.15g : % 22.15g % 22.15g % 22.15g % 22.15g\n",
 *       - (log ((s + c) / slice_size1) - lnNorma + 0.5 * self->bestfit),
 *       - (log (ncm_stats_vec_get_mean (slnnorm, 0)) - lnNorma + 0.5 * self->bestfit),
 *       ncm_stats_vec_get_mean (slnnorm, 1),
 *       ncm_stats_vec_get_sd (slnnorm, 1),
 *       ncm_stats_vec_get_sd (slnnorm, 1) / sqrt (j + 1.0)
 *   );
 */
  }

  lnnorm_sd[0] = ncm_stats_vec_get_sd (slnnorm, 1) / sqrt (1.0 * nslices);
  mean_lnnorm  = -(log (ncm_stats_vec_get_mean (slnnorm, 0)) - lnNorma + 0.5 * self->bestfit);

  ncm_vector_free (v);
  ncm_stats_vec_clear (&slnnorm);

  return mean_lnnorm;
}

static gdouble
_ncm_mset_catalog_get_post_lnnorm_sum_bs (NcmMSetCatalog *mcat, NcmVector *mean, NcmMatrix *cov, const gdouble lnNorma, const gdouble reltol, gdouble *lnnorm_sd, NcmRNG *rng)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint max_iter        = 1000000;
  const guint fparams_len     = ncm_mset_fparams_len (self->mset);
  const guint cat_len         = ncm_mset_catalog_len (mcat);
  NcmVector *v                = ncm_vector_new (fparams_len);
  NcmVector *m2lnp_v          = ncm_vector_new (cat_len);
  NcmBootstrap *bs            = ncm_bootstrap_sized_new (cat_len);
  NcmStatsVec *slnnorm        = ncm_stats_vec_new (2, NCM_STATS_VEC_VAR, FALSE);
  gdouble s                   = 0.0;
  gdouble c                   = 0.0;
  gdouble mean_lnnorm;
  gint ret;
  guint i, j;

  for (i = 0; i < cat_len; i++)
  {
    NcmVector *row_i = ncm_mset_catalog_peek_row (mcat, i);
    gdouble m2lnp_i  = 0.0;

    ncm_vector_memcpy2 (v, row_i, 0, self->nadd_vals, fparams_len);
    ncm_vector_sub (v, mean);

    ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                          ncm_matrix_gsl (cov), ncm_vector_gsl (v));
    NCM_TEST_GSL_RESULT ("ncm_mset_catalog_get_post_lnnorm", ret);

    ret = gsl_blas_ddot (ncm_vector_gsl (v), ncm_vector_gsl (v), &m2lnp_i);
    NCM_TEST_GSL_RESULT ("ncm_mset_catalog_get_post_lnnorm", ret);

    ncm_vector_set (m2lnp_v, i, m2lnp_i);
  }

  for (j = 0; j < max_iter; j++)
  {
    GArray *bs_array;
    guint bs_size;

    ncm_bootstrap_resample (bs, rng);

    bs_array = ncm_bootstrap_get_sortncomp (bs);
    s        = 0.0;
    c        = 0.0;

    bs_size = bs_array->len / 2;

    for (i = 0; i < bs_size; i++)
    {
      const guint c_i       = g_array_index (bs_array, guint, 2 * i + 0);
      const gdouble f_i     = g_array_index (bs_array, guint, 2 * i + 1);
      const gdouble m2lnL_i = ncm_vector_fast_get (ncm_mset_catalog_peek_row (mcat, c_i), self->m2lnp_var);
      const gdouble m2lnp_i = ncm_vector_fast_get (m2lnp_v, c_i);
      const gdouble e_i     = f_i * exp (0.5 * ((m2lnL_i - self->bestfit) - m2lnp_i));
      const gdouble t       = s + e_i;

      c += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
      s  = t;
    }

    ncm_stats_vec_set (slnnorm, 0, (s + c) / cat_len);
    ncm_stats_vec_set (slnnorm, 1, -(log ((s + c) / cat_len) - lnNorma + 0.5 * self->bestfit));

    ncm_stats_vec_update (slnnorm);

    g_array_unref (bs_array);

    /* Stops when the mean over the resamples is known to reltol. */
    if (j > 10)
    {
      const gdouble norm_mean    = ncm_stats_vec_get_mean (slnnorm, 0);
      const gdouble norm_mean_sd = ncm_stats_vec_get_sd (slnnorm, 0) / sqrt (j + 1.0);

      if (norm_mean_sd / norm_mean < reltol)
        break;
    }
  }

  lnnorm_sd[0] = ncm_stats_vec_get_sd (slnnorm, 1);
  mean_lnnorm  = -(log (ncm_stats_vec_get_mean (slnnorm, 0)) - lnNorma + 0.5 * self->bestfit);

  ncm_vector_free (v);
  ncm_vector_free (m2lnp_v);
  ncm_bootstrap_clear (&bs);
  ncm_stats_vec_clear (&slnnorm);

  return mean_lnnorm;
}

static gdouble
_ncm_mset_catalog_get_post_lnnorm_hyperbox (NcmMSetCatalog *mcat, gboolean use_bs, gdouble *post_lnnorm_sd)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint fparams_len     = ncm_mset_fparams_len (self->mset);
  NcmMatrix *cov              = NULL;
  NcmVector *mean             = NULL;
  NcmRNG *rng                 = ncm_rng_new (NULL);
  gdouble lnNorma             = 0.0;
  NcmMatrix *cov_mvnd;
  gdouble ratio, post_lnnorm;
  gint ret;

  ncm_mset_catalog_get_covar (mcat, &cov);
  ncm_mset_catalog_get_mean (mcat, &mean);
  ncm_matrix_scale (cov, gsl_pow_2 (NCM_MSET_CATALOG_RESCALE_COV));

  /* Factored before the ratio below, which would abort on the same covariance. */
  cov_mvnd = ncm_matrix_dup (cov);
  ret      = ncm_matrix_cholesky_decomp (cov, 'U');

  if (ret != 0)
  {
    g_warning ("ncm_mset_catalog_get_post_lnnorm[ncm_matrix_cholesky_decomp]: %d. Non-positive definite covariance, more points are necessary.", ret);

    ncm_vector_clear (&mean);
    ncm_matrix_clear (&cov);
    ncm_matrix_clear (&cov_mvnd);
    ncm_rng_clear (&rng);
    post_lnnorm_sd[0] = GSL_NAN;

    return 0.0;
  }

  {
    NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new (fparams_len);
    NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (fparams_len);
    NcmMSet *mset_mvnd             = ncm_mset_new (model_mvnd, NULL, NULL);

    ncm_data_gauss_cov_mvnd_set_cov_mean (data_mvnd, mean, cov_mvnd);

    ncm_mset_param_set_all_ftype (mset_mvnd, NCM_PARAM_TYPE_FREE);
    ncm_mset_prepare_fparam_map (mset_mvnd);
    ncm_mset_fparams_set_vector (mset_mvnd, mean);

    ratio = ncm_data_gauss_cov_mvnd_est_ratio (data_mvnd, mset_mvnd, self->mset, (NcmDataGaussCovMVNDBound) ncm_mset_fparam_valid_bounds, NULL, NULL, 1.0e-3, rng);

    ncm_mset_clear (&mset_mvnd);
    ncm_model_mvnd_clear (&model_mvnd);
    ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  }

  ncm_matrix_clear (&cov_mvnd);

  lnNorma = 0.5 * (fparams_len * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (cov)) + log (ratio);

  if (use_bs)
    post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_sum_bs (mcat, mean, cov, lnNorma, 1.0e-2, post_lnnorm_sd, rng);
  else
    post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_sum (mcat, mean, cov, lnNorma, GSL_POSINF, post_lnnorm_sd);

  ncm_vector_clear (&mean);
  ncm_matrix_clear (&cov);
  ncm_rng_clear (&rng);

  return post_lnnorm;
}

/**
 * ncm_mset_catalog_get_post_lnnorm:
 * @mcat: a #NcmMSetCatalog
 * @post_lnnorm_sd: (out): standard deviation of the estimate
 *
 * Estimates $\ln Z$, $Z = \int e^{-m/2}\,\mathrm{d}\theta$ over the box of the
 * parameter bounds, $m$ the #NcmMSetCatalog:m2lnp-var column: the log evidence for a
 * flat prior of unit density on the box (for a normalized flat prior subtract the log of
 * the box volume). The estimate averages $g/e^{-m/2}$ over the rows, $g$ the Gaussian
 * with the rows' mean and $0.8^2$ times their covariance, normalized as
 * #NcmMSetCatalog:post-lnnorm-method sets (#NcmMSetCatalogPostNormMethod): to the box, by
 * a Monte Carlo estimate of its mass there, or to an ellipsoid inside it. The result is
 * kept until rows are added or the method changes.
 *
 * @post_lnnorm_sd is the spread of the estimate over slices of the rows, or over
 * bootstrap resamples; below 1000 rows the slices give none, and it is NaN with a
 * warning. Without #NcmMSetCatalog:m2lnp-var it warns and returns zero; a covariance
 * that is not positive definite warns and gives zero, with a NaN error.
 *
 * Returns: the estimate of $\ln Z$
 */
gdouble
ncm_mset_catalog_get_post_lnnorm (NcmMSetCatalog *mcat, gdouble *post_lnnorm_sd)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->m2lnp_var < 0)
  {
    g_warning ("ncm_mset_catalog_get_post_lnnorm: m2lnp variable not set, cannot calculate the posterior normalization.");

    return 0.0;
  }

  if (!self->post_lnnorm_up)
  {
    switch (self->post_lnnorm_method)
    {
      case NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX:
        self->post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_hyperbox (mcat, FALSE, &self->post_lnnorm_sd);
        break;
      case NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX_BS:
        self->post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_hyperbox (mcat, TRUE, &self->post_lnnorm_sd);
        break;
      case NCM_MSET_CATALOG_POST_LNNORM_METHOD_ELLIPSOID:
        self->post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_ellipsoid (mcat, &self->post_lnnorm_sd);
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    self->post_lnnorm_up = TRUE;
  }

  /* The cached value comes with its error. */
  post_lnnorm_sd[0] = self->post_lnnorm_sd;

  return self->post_lnnorm;
}

void
_ncm_mset_catalog_sort_by_m2lnp (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->m2lnp_var < 0)
  {
    g_error ("_ncm_mset_catalog_sort_by_m2lnp: m2lnp variable not set, cannot sort the catalog.");

    return;
  }

  if (!self->order_sorted)
  {
    g_array_sort (self->order_m2lnp, &_ncm_mset_catalog_double_compare);
    self->order_sorted = TRUE;
  }
}

/**
 * ncm_mset_catalog_get_post_lnvol:
 * @mcat: a #NcmMSetCatalog
 * @level: fraction of the posterior, in $(0, 1)$
 * @glnvol: (out) (nullable): log volume of the region for the Gaussian approximation
 *
 * Estimates the log volume of the highest-posterior region holding the fraction @level
 * of the rows: $Z$ (ncm_mset_catalog_get_post_lnnorm()) times the average over the rows
 * of $e^{m/2}$ on that region, $m$ the #NcmMSetCatalog:m2lnp-var column. @glnvol gets
 * the log volume of the same region for a Gaussian with the rows' covariance $C$, the
 * ellipsoid $\theta^T C^{-1}\theta \leq \chi^2_n(\text{@level})$. With fewer rows than
 * $1/\text{@level}$ it warns and returns zero.
 *
 * Returns: the estimate of the log volume
 */
gdouble
ncm_mset_catalog_get_post_lnvol (NcmMSetCatalog *mcat, const gdouble level, gdouble *glnvol)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gdouble post_lnnorm_sd;
  const gdouble lnnorm = ncm_mset_catalog_get_post_lnnorm (mcat, &post_lnnorm_sd);
  const guint cat_len  = ncm_mset_catalog_len (mcat);
  const gint nitens    = cat_len * level;
  gdouble s            = 0.0;
  gint i;

  g_assert_cmpfloat (level, >, 0.0);
  g_assert_cmpfloat (level, <, 1.0);

  if (nitens <= 0)
  {
    g_warning ("ncm_mset_catalog_get_post_lnvol: too few points in the catalog `%u', or level too close to zero %2f%%.\n",
               cat_len, level * 100.0);

    if (glnvol != NULL)
      glnvol[0] = GSL_NAN;

    return 0.0;
  }

  _ncm_mset_catalog_sort_by_m2lnp (mcat);

  for (i = 0; i < nitens; i++)
  {
    const gdouble m2lnp_i = g_array_index (self->order_m2lnp, gdouble, i);

    s += exp (0.5 * (m2lnp_i - self->bestfit));
  }

  if (glnvol != NULL)
  {
    NcmMatrix *cov             = NULL;
    gint signp                 = 0;
    const guint fparams_len    = ncm_mset_fparams_len (self->mset);
    const gdouble lnVnball     = fparams_len * ncm_c_lnpi () / 2.0 - lgamma_r (fparams_len / 2.0 + 1.0, &signp);
    const gdouble lnsigma_size = log (gsl_cdf_chisq_Pinv (level, fparams_len));

    ncm_mset_catalog_get_covar (mcat, &cov);

    ncm_matrix_cholesky_decomp (cov, 'U');

    glnvol[0] = lnVnball + 0.5 * ncm_matrix_cholesky_lndet (cov) + 0.5 * fparams_len * lnsigma_size;

    ncm_matrix_clear (&cov);
  }

  return lnnorm + 0.5 * self->bestfit + log (s / cat_len);
}

/**
 * ncm_mset_catalog_get_nth_m2lnL_percentile:
 * @mcat: a #NcmMSetCatalog
 * @p: fraction of the rows, in $(0, 1)$
 * @nth: (out) (allow-none): number of rows with a smaller value
 *
 * Finds the value of the #NcmMSetCatalog:m2lnp-var column at sorted position
 * $\lfloor p\,n\rfloor$: the fraction @p of the $n$ rows have a smaller value.
 *
 * Returns: the @p quantile of $-2\ln L$
 */
gdouble
ncm_mset_catalog_get_nth_m2lnL_percentile (NcmMSetCatalog *mcat, const gdouble p, guint *nth)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint cat_len         = ncm_mset_catalog_len (mcat);
  const gint nitens           = cat_len * p;
  gdouble m2lnL               = 0.0;

  g_assert_cmpfloat (p, >, 0.0);
  g_assert_cmpfloat (p, <, 1.0);

  if (nitens <= 0)
  {
    g_warning ("ncm_mset_catalog_get_nth_m2lnL_percentile: too few points in the catalog `%u', or percentile too close to zero %2f%%.\n",
               cat_len, p * 100.0);

    return 0.0;
  }

  _ncm_mset_catalog_sort_by_m2lnp (mcat);

  m2lnL = g_array_index (self->order_m2lnp, gdouble, nitens);

  if (nth != NULL)
    nth[0] = nitens;

  return m2lnL;
}

/**
 * ncm_mset_catalog_get_bestfit_m2lnL:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the smallest value of the #NcmMSetCatalog:m2lnp-var column so far
 */
gdouble
ncm_mset_catalog_get_bestfit_m2lnL (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->bestfit;
}

/**
 * ncm_mset_catalog_get_bestfit_row:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: (transfer full) (nullable): a copy of the row with the smallest
 *   #NcmMSetCatalog:m2lnp-var, or %NULL for an empty catalog
 */
NcmVector *
ncm_mset_catalog_get_bestfit_row (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->bestfit_row)
    return ncm_vector_dup (self->bestfit_row);
  else
    return NULL;
}

/**
 * ncm_mset_catalog_get_mean:
 * @mcat: a #NcmMSetCatalog
 * @mean: (inout) (allow-none) (transfer full): a #NcmVector, allocated when *@mean is %NULL
 *
 * Sets *@mean to the mean of the free parameters over the rows (weighted for a weighted
 * catalog).
 *
 */
void
ncm_mset_catalog_get_mean (NcmMSetCatalog *mcat, NcmVector **mean)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (*mean == NULL)
    *mean = ncm_vector_new (self->ncols - self->nadd_vals);

  ncm_stats_vec_get_mean_vector (self->pstats, *mean, self->nadd_vals);
}

/**
 * ncm_mset_catalog_get_covar:
 * @mcat: a #NcmMSetCatalog
 * @cov: (inout) (allow-none) (transfer full): a #NcmMatrix, allocated when *@cov is %NULL
 *
 * Sets *@cov to the covariance of the free parameters over the rows.
 *
 */
void
ncm_mset_catalog_get_covar (NcmMSetCatalog *mcat, NcmMatrix **cov)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (*cov == NULL)
    *cov = ncm_matrix_new (self->ncols - self->nadd_vals, self->ncols - self->nadd_vals);

  ncm_stats_vec_get_cov_matrix (self->pstats, *cov, self->nadd_vals);
}

/**
 * ncm_mset_catalog_get_full_covar:
 * @mcat: a #NcmMSetCatalog
 * @cov: (inout) (allow-none) (transfer full): a #NcmMatrix
 *
 * Sets *@cov to the covariance of all the columns (additional values and free
 * parameters) over the rows.
 *
 */
void
ncm_mset_catalog_get_full_covar (NcmMSetCatalog *mcat, NcmMatrix **cov)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (*cov == NULL)
    *cov = ncm_matrix_new (self->ncols, self->ncols);

  ncm_stats_vec_get_cov_matrix (self->pstats, *cov, 0);
}

/**
 * ncm_mset_catalog_log_full_covar:
 * @mcat: a #NcmMSetCatalog
 *
 * Logs the covariance of all the columns (additional values and free parameters).
 *
 */
void
ncm_mset_catalog_log_full_covar (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint params_len      = self->ncols;
  const gchar *box            = "---------------";
  guint name_size             = ncm_mset_max_fparam_name (self->mset);
  guint i, j;

  for (i = 0; i < self->add_vals_names->len; i++)
  {
    name_size = GSL_MAX (name_size, strlen (g_ptr_array_index (self->add_vals_names, i)));
  }

  ncm_cfg_msg_sepa ();
  g_message ("# NcmMSetCatalog full covariance matrix\n");
  g_message ("#                                                           ");

  for (i = 0; i < name_size; i++)
    g_message (" ");

  for (i = 0; i < params_len; i++)
    i ? g_message ("%s", box) : g_message ("-%s", box);

  if (i)
    g_message ("\n");

  for (i = 0; i < params_len; i++)
  {
    if (i < self->nadd_vals)
    {
      const gchar *pname = g_ptr_array_index (self->add_vals_names, i);

      g_message ("# %*s[%05d:%02d] = % -12.4g (% -12.4g) +/- % -12.4g |",
                 name_size, pname, 99999, i,
                 ncm_stats_vec_get_mean (self->pstats, i),
                 self->bestfit_row != NULL ? ncm_vector_get (self->bestfit_row, i) : GSL_NAN,
                 ncm_stats_vec_get_sd (self->pstats, i)
      );
    }
    else
    {
      const gint fpi          = i - self->nadd_vals;
      const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (self->mset, fpi);
      const gchar *pname      = ncm_mset_fparam_name (self->mset, fpi);

      g_message ("# %*s[%05d:%02d] = % -12.4g (% -12.4g) +/- % -12.4g |",
                 name_size, pname, pi->mid, pi->pid,
                 ncm_stats_vec_get_mean (self->pstats, i),
                 self->bestfit_row != NULL ? ncm_vector_get (self->bestfit_row, i) : GSL_NAN,
                 ncm_stats_vec_get_sd (self->pstats, i)
      );
    }

    for (j = 0; j < params_len; j++)
    {
      g_message (" % -12.4g |", ncm_stats_vec_get_cor (self->pstats, i, j));
    }

    g_message ("\n");
  }

  g_message ("#                                                           ");

  for (i = 0; i < name_size; i++)
    g_message (" ");

  for (i = 0; i < params_len; i++)
    i ? g_message ("%s", box) : g_message ("-%s", box);

  if (i)
    g_message ("\n");

  return;
}

/**
 * ncm_mset_catalog_estimate_autocorrelation_tau:
 * @mcat: a #NcmMSetCatalog
 * @force_single_chain: whether to force the catalog to be treated as a single chain
 *
 * Updates the internal estimates of the integrated autocorrelation time, in ensemble
 * iterations, from the accumulator the catalog keeps as rows are added.
 *
 * With @force_single_chain the interleaved rows are treated as one series instead, which
 * is only meaningful for a catalog whose rows are a single chain.
 *
 */
void
ncm_mset_catalog_estimate_autocorrelation_tau (NcmMSetCatalog *mcat, gboolean force_single_chain)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint total           = ncm_vector_len (self->tau);
  guint p;

  _ncm_mset_catalog_check_unweighted (self, G_STRFUNC);

  if (force_single_chain && (self->nchains > 1))
  {
    NcmStatsAcorr *rows = ncm_stats_acorr_new (total);
    const guint nitens  = ncm_stats_vec_nitens (self->pstats);
    guint i;

    for (i = 0; i < nitens; i++)
      ncm_stats_acorr_update (rows, ncm_mset_catalog_peek_row (mcat, i));

    for (p = 0; p < total; p++)
      ncm_vector_set (self->tau, p, ncm_stats_acorr_get_tau (rows, p));

    ncm_stats_acorr_free (rows);

    return;
  }

  for (p = 0; p < total; p++)
    ncm_vector_set (self->tau, p, ncm_stats_acorr_get_tau (self->acorr, p));
}

/**
 * ncm_mset_catalog_get_tau_diag:
 * @mcat: a #NcmMSetCatalog
 * @p: parameter id
 *
 * Conditions attached to the autocorrelation time of parameter @p. Anything other than
 * %NCM_STATS_ACORR_DIAG_OK means the estimate is not to be read as converged.
 *
 * Returns: the #NcmStatsAcorrDiag flags.
 */
NcmStatsAcorrDiag
ncm_mset_catalog_get_tau_diag (NcmMSetCatalog *mcat, guint p)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_check_unweighted (self, G_STRFUNC);

  return ncm_stats_acorr_get_diag (self->acorr, p);
}

/**
 * ncm_mset_catalog_log_tau_diag:
 * @mcat: a #NcmMSetCatalog
 *
 * Logs the parameters whose autocorrelation time carries a condition, and which one. A
 * catalog with nothing to report logs nothing.
 *
 * Returns: TRUE when at least one parameter is flagged.
 */
gboolean
ncm_mset_catalog_log_tau_diag (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint total           = ncm_vector_len (self->tau);
  gboolean any                = FALSE;
  guint p;

  _ncm_mset_catalog_check_unweighted (self, G_STRFUNC);

  for (p = 0; p < total; p++)
  {
    const NcmStatsAcorrDiag diag = ncm_stats_acorr_get_diag (self->acorr, p);

    if (diag != NCM_STATS_ACORR_DIAG_OK)
    {
      gchar *diag_str = ncm_stats_acorr_diag_to_string (diag);

      if (!any)
        g_message ("# NcmMSetCatalog: the autocorrelation time of the following "
                   "parameters is not to be read as converged:\n");

      g_message ("#   %-30s tau = % 12.4g  %s\n",
                 ncm_mset_catalog_col_full_name (mcat, p),
                 ncm_vector_get (self->tau, p), diag_str);

      g_free (diag_str);
      any = TRUE;
    }
  }

  return any;
}

/**
 * ncm_mset_catalog_tau_needs_more:
 * @mcat: a #NcmMSetCatalog
 * @required_niter: (out) (optional): iterations the reliability factor asks for
 *
 * Whether any free parameter carries %NCM_STATS_ACORR_DIAG_SHORT_CHAIN, that is, whether
 * the chain is shorter than #NcmStatsAcorr:reliability-factor times its own $\tau$. That
 * is the one condition sampling more clears on its own, and @required_niter is the length
 * at which it does: the reliability factor times the largest $\tau$ over the free
 * parameters.
 *
 * The other conditions are reported by ncm_mset_catalog_log_tau_diag() and do not ask for
 * more sampling, because more sampling does not remove them. A burn-in still in the
 * catalog keeps %NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT and
 * %NCM_STATS_ACORR_DIAG_DRIFT set however long the chain grows, since the early
 * iterations stay in the first half of it; what removes those is trimming the catalog. A
 * parameter that never moved keeps %NCM_STATS_ACORR_DIAG_ZERO_VARIANCE. Two estimators
 * differing is a statement about the estimate and not about the length of the chain.
 *
 * Two conditions are refused outright, because the effective sample size they lead to is
 * not a sample size: a free parameter whose ensemble mean has no variance
 * (%NCM_STATS_ACORR_DIAG_ZERO_VARIANCE: the walkers did not move), and a free parameter
 * whose $K_\mathrm{eff}$ exceeds twice the number of chains, which is what walkers frozen
 * at distinct positions look like (the rows vary, the ensemble mean does not). Both return
 * TRUE with @required_niter as above; #NcmFitESMCMC bounds the number of rounds.
 *
 * Returns: TRUE when more sampling is called for.
 */
gboolean
ncm_mset_catalog_tau_needs_more (NcmMSetCatalog *mcat, guint *required_niter)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint free_params_len = ncm_mset_fparams_len (self->mset);
  const guint fpi             = self->nadd_vals;
  const guint fpf             = free_params_len + self->nadd_vals;
  const gdouble factor        = ncm_stats_acorr_get_reliability_factor (self->acorr);
  gdouble max_tau             = 1.0;
  gboolean needs              = FALSE;
  guint p;

  _ncm_mset_catalog_check_unweighted (self, G_STRFUNC);

  for (p = fpi; p < fpf; p++)
  {
    const NcmStatsAcorrDiag diag = ncm_stats_acorr_get_diag (self->acorr, p);

    if (diag & (NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_ZERO_VARIANCE))
      needs = TRUE;

    /* K_eff is estimated from niter ensemble means, so it scatters around K by
     * ~sqrt (2 / niter); frozen walkers at distinct positions give K_eff >> K. */
    if ((self->nchains > 1) && (ncm_mset_catalog_get_keff (mcat, p) > 2.0 * self->nchains))
      needs = TRUE;

    max_tau = GSL_MAX (max_tau, ncm_stats_acorr_get_tau (self->acorr, p));
  }

  if (required_niter != NULL)
    required_niter[0] = (guint) ceil (factor * max_tau);

  return needs;
}

/**
 * ncm_mset_catalog_get_keff:
 * @mcat: a #NcmMSetCatalog
 * @p: parameter id
 *
 * Effective number of independent chains per iteration for parameter @p, the variance
 * over all rows divided by the variance of the ensemble mean. It is the number of chains
 * when they are independent at a given iteration, one when they move together, and above
 * the number of chains when the ensemble mean is steadier than independent chains would
 * make it, which is what walkers held in place at different positions look like.
 *
 * For a catalog with a single chain it is one.
 *
 * Returns: the effective number of chains.
 */
gdouble
ncm_mset_catalog_get_keff (NcmMSetCatalog *mcat, guint p)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gdouble var_e;

  _ncm_mset_catalog_check_unweighted (self, G_STRFUNC);

  if (self->nchains == 1)
    return 1.0;

  var_e = ncm_stats_acorr_get_var (self->acorr, p);

  if (var_e <= 0.0)
    return 1.0;

  return ncm_stats_vec_get_var (self->pstats, p) / var_e;
}

/**
 * ncm_mset_catalog_get_ess:
 * @mcat: a #NcmMSetCatalog
 * @p: parameter id
 *
 * Effective sample size of parameter @p over the whole catalog,
 * $K_\mathrm{eff}\,n_\mathrm{iter}/\tau$. No assumption is made about the chains being
 * independent within an iteration: whatever correlation they carry is already in
 * $K_\mathrm{eff}$, see ncm_mset_catalog_get_keff().
 *
 * Returns: the effective sample size.
 */
gdouble
ncm_mset_catalog_get_ess (NcmMSetCatalog *mcat, guint p)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return ncm_mset_catalog_get_keff (mcat, p) * ncm_stats_acorr_get_ess (self->acorr, p);
}

/**
 * ncm_mset_catalog_peek_autocorrelation_tau:
 * @mcat: a #NcmMSetCatalog
 *
 * Gets the last estimate of the autocorrelation tau calculated
 * in the last call of ncm_mset_catalog_estimate_autocorrelation_tau().
 *
 * Returns: (transfer none): the last estimate of the autocorrelation tau.
 */
NcmVector *
ncm_mset_catalog_peek_autocorrelation_tau (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->tau;
}

/**
 * ncm_mset_catalog_get_param_shrink_factor:
 * @mcat: a #NcmMSetCatalog
 * @p: parameter id
 *
 * The potential scale reduction factor of column @p [Gelman and Rubin
 * (1992)](https://doi.org/10.1214/ss/1177011136),
 * $$\hat{R} = \sqrt{\frac{n - 1}{n} + \frac{m + 1}{m}\frac{B}{nW}},$$
 * where $n$ is the length of one chain, $m$ the number of chains, $W$ the mean of the
 * variances of the chains and $B/n$ the variance of their means. A catalog with a single
 * chain returns one. It warns when the chains do not all have the same length.
 *
 * Returns: the shrink factor of @p.
 */
gdouble
ncm_mset_catalog_get_param_shrink_factor (NcmMSetCatalog *mcat, guint p)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint i;
  gdouble W, B_n, shrink_factor;
  guint n = ncm_stats_vec_nitens (self->pstats);

  if (self->nchains == 1)
    return 1.0;

  if (n % self->nchains != 0)
    g_warning ("ncm_mset_catalog_get_param_shrink_factor: not all chains have the same size [%u %u] %u.", n, self->nchains, (n % self->nchains));

  n = n / self->nchains;

  for (i = 0; i < self->nchains; i++)
  {
    NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, i);

    ncm_vector_set (self->chain_means, i, ncm_stats_vec_get_mean (pstats, p));
    ncm_vector_set (self->chain_vars, i, ncm_stats_vec_get_var (pstats, p));
  }

  W   = gsl_stats_mean (ncm_vector_ptr (self->chain_vars, 0), ncm_vector_stride (self->chain_vars), ncm_vector_len (self->chain_vars));
  B_n = gsl_stats_variance (ncm_vector_ptr (self->chain_means, 0), ncm_vector_stride (self->chain_means), ncm_vector_len (self->chain_means));

  shrink_factor = sqrt ((n - 1.0) /  (1.0 * n) + (self->nchains + 1.0) / (1.0 * self->nchains) * B_n / W);

  return shrink_factor;
}

/**
 * ncm_mset_catalog_get_shrink_factor:
 * @mcat: a #NcmMSetCatalog
 *
 * The multivariate potential scale reduction factor over the free parameters [Brooks and
 * Gelman (1998)](https://doi.org/10.1080/10618600.1998.10474787),
 * $$\hat{R}^p = \sqrt{\frac{n - 1}{n} + \frac{m + 1}{m} \lambda_1},$$
 * where $n$ is the length of one chain, $m$ the number of chains and $\lambda_1$ the
 * largest eigenvalue of $W^{-1}B/n$, with $W$ the mean of the covariance matrices of the
 * chains and $B/n$ the covariance of their means.
 *
 * A catalog with a single chain returns one. When the first entry of $W$ is not finite
 * or $W$ is not positive definite the factor is not computed and the return value is
 * $10^{10}$. It warns when the chains do not all have the same length.
 *
 * Returns: the shrink factor $\hat{R}^p$.
 */
gdouble
ncm_mset_catalog_get_shrink_factor (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gint ret;
  guint i;
  guint n                     = ncm_stats_vec_nitens (self->pstats);
  const guint free_params_len = self->ncols - self->nadd_vals;
  gdouble shrink_factor       = 1.0e10;

  if (self->nchains == 1)
    return 1.0;

  if (n % self->nchains != 0)
    g_warning ("ncm_mset_catalog_get_shrink_factor: not all chains have the same size [%u %u] %u.", n, self->nchains, (n % self->nchains));

  n = n / self->nchains;

  ncm_stats_vec_reset (self->mean_pstats, TRUE);
  ncm_matrix_set_zero (self->chain_cov);

  for (i = 0; i < self->nchains; i++)
  {
    NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, i);
    NcmMatrix *cov      = ncm_stats_vec_peek_cov_matrix (pstats, self->nadd_vals);
    guint p;

    for (p = 0; p < free_params_len; p++)
    {
      ncm_stats_vec_set (self->mean_pstats, p, ncm_stats_vec_get_mean (pstats, p + self->nadd_vals));
      /*printf ("chain %u param %u mean % 20.15g\n", i, p, ncm_stats_vec_get_mean (pstats, p + self->nadd_vals));*/
    }

    ncm_stats_vec_update (self->mean_pstats);

    ncm_matrix_add_mul (self->chain_cov, 1.0, cov);
  }

  ncm_matrix_scale (self->chain_cov, 1.0 / (1.0 * self->nchains));

  {
    NcmMatrix *cov = ncm_stats_vec_peek_cov_matrix (self->mean_pstats, 0);

/*
 *   ncm_matrix_log_vals (self->chain_cov, "# mean cov", "% 10.5g");
 *   ncm_matrix_log_vals (cov, "# cov mean", "% 10.5g");
 */
    if (gsl_finite (ncm_matrix_get (self->chain_cov, 0, 0)))
    {
      do {
        gdouble lev = 0.0;

        ret = ncm_matrix_cholesky_decomp (self->chain_cov, 'U');

        if (ret != 0)
          break;

        ret = ncm_matrix_cholesky_inverse (self->chain_cov, 'U');

        if (ret != 0)
          break;

        ncm_matrix_dsymm (self->chain_sM, 'U', 1.0, self->chain_cov, cov, 0.0);

        gsl_eigen_nonsymm_params (0, 1, self->chain_sM_ws);
        gsl_eigen_nonsymm (ncm_matrix_gsl (self->chain_sM), self->chain_sM_ev, self->chain_sM_ws);

        for (i = 0; i < free_params_len; i++)
        {
          lev = GSL_MAX (lev, GSL_VECTOR_REAL (self->chain_sM_ev, i));

          /*if (GSL_VECTOR_IMAG (self->chain_sM_ev, i) != 0.0) */
          /*  g_warning ("ncm_mset_catalog_get_shrink_factor: complex eigenvalue in SM matrix, unreliable shrink factor, try using more chains."); */
        }

        shrink_factor = sqrt ((n - 1.0) / (1.0 * n) + (self->nchains + 1.0) * lev / (self->nchains * 1.0));
      } while (FALSE);
    }
  }

  return shrink_factor;
}

/**
 * ncm_mset_catalog_param_pdf:
 * @mcat: a #NcmMSetCatalog
 * @i: column index
 *
 * Bins column @i of the rows into a histogram for ncm_mset_catalog_param_pdf_pvalue().
 * The histogram has one bin per ten rows, and at least ten, spanning the sampled range
 * of the column. Row weights are not used.
 *
 */
void
ncm_mset_catalog_param_pdf (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint n               = ncm_stats_vec_nitens (self->pstats);
  const guint nbins           = n / 10 >= 10 ? n / 10 : 10;
  const gdouble p_max         = ncm_vector_get (self->params_max, i);
  const gdouble p_min         = ncm_vector_get (self->params_min, i);
  guint k;

  self->pdf_i = i;

  if ((self->h != NULL) && (self->h->n != nbins))
  {
    gsl_histogram_free (self->h);
    self->h = NULL;
  }

  if (self->h == NULL)
    self->h = gsl_histogram_alloc (nbins);

  if ((self->h_pdf != NULL) && (self->h_pdf->n != nbins))
  {
    gsl_histogram_pdf_free (self->h_pdf);
    self->h_pdf = NULL;
  }

  if (self->h_pdf == NULL)
    self->h_pdf = gsl_histogram_pdf_alloc (nbins);

  gsl_histogram_set_ranges_uniform (self->h, p_min, p_max);
  /* The bins are half open, so the last edge moves just past the largest value. */
  self->h->range[nbins] = nextafter (p_max, GSL_POSINF);

  for (k = 0; k < ncm_stats_vec_nitens (self->pstats); k++)
  {
    NcmVector *row = ncm_mset_catalog_peek_row (mcat, k);

    gsl_histogram_increment (self->h, ncm_vector_get (row, i));
  }

  gsl_histogram_pdf_init (self->h_pdf, self->h);
}

/**
 * ncm_mset_catalog_param_pdf_pvalue:
 * @mcat: a #NcmMSetCatalog
 * @pvalue: value of the column binned by ncm_mset_catalog_param_pdf()
 * @both: unused
 *
 * The upper tail p-value of @pvalue in the histogram of the last call to
 * ncm_mset_catalog_param_pdf(): the fraction of rows at or above the lower edge of the
 * bin holding @pvalue. A value outside the sampled range warns and gives one below it,
 * zero above it.
 *
 * Returns: the p-value.
 */
gdouble
ncm_mset_catalog_param_pdf_pvalue (NcmMSetCatalog *mcat, gdouble pvalue, gboolean both)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert_cmpint (self->pdf_i, >=, 0);
  g_assert (self->h_pdf != NULL);

  {
    const gdouble p_max = ncm_vector_get (self->params_max, self->pdf_i);
    const gdouble p_min = ncm_vector_get (self->params_min, self->pdf_i);
    gsize i             = 0;

    NCM_UNUSED (both);

    if ((pvalue < p_min) || (pvalue > p_max))
    {
      const gdouble p = (pvalue < p_min) ? 1.0 : 0.0;

      g_warning ("ncm_mset_catalog_param_pdf_pvalue: value % 20.15g outside the sampled interval [% 20.15g % 20.15g]. Assuming p-value %g.",
                 pvalue, p_min, p_max, p);

      return p;
    }

    gsl_histogram_find (self->h, pvalue, &i);
    g_assert_cmpint (i, <, self->h_pdf->n);

    /* sum[i] is the fraction of rows below the lower edge of bin i. */
    return 1.0 - self->h_pdf->sum[i];
  }
}

/**
 * ncm_mset_catalog_calc_ci_direct:
 * @mcat: a #NcmMSetCatalog
 * @func: a #NcmMSetFunc of dimension one and one argument
 * @x_v: #NcmVector of arguments of @func
 * @p_val: (element-type double): probabilities of the confidence intervals, in $(0, 1)$
 *
 * The mean and the central confidence intervals of @func over the rows, at each
 * argument in @x_v. Row $k$ of the result belongs to the $k$-th argument; column 0 holds
 * the mean and columns $1 + 2j$ and $2 + 2j$ the lower and upper bounds of the interval
 * of probability $p_j$, the quantiles $(1 - p_j)/2$ and $(1 + p_j)/2$. For @p_val =
 * (0.6827, 0.9545) the columns are the mean, the $1\sigma$ bounds and the $2\sigma$
 * bounds.
 *
 * The quantiles are those of the sorted values, so the function holds the catalog length
 * times the length of @x_v in memory; ncm_mset_catalog_calc_ci_interp() does not. Row
 * weights are not used. The free parameters of the #NcmMSet are restored at the end.
 *
 * Returns: (transfer full): a #NcmMatrix with the means and the interval bounds.
 */
NcmMatrix *
ncm_mset_catalog_calc_ci_direct (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmVector *x_v, GArray *p_val)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint dim             = ncm_vector_len (x_v);

  g_assert_cmpuint (p_val->len, >, 0);
  {
    const guint nelem      = p_val->len * 2 + 1;
    NcmMatrix *res         = ncm_matrix_new (dim, nelem);
    NcmVector *save_params = ncm_vector_new (ncm_mset_fparams_len (self->mset));
    const guint cat_len    = ncm_mset_catalog_len (mcat);
    guint i, j;

    ncm_mset_fparams_get_vector (self->mset, save_params);

    ncm_vector_clear (&self->quantile_ws);
    self->quantile_ws = ncm_vector_new (cat_len * dim);

    for (i = 0; i < cat_len; i++)
    {
      NcmVector *row = ncm_mset_catalog_peek_row (mcat, i);
      NcmVector *qws = ncm_vector_get_subvector (self->quantile_ws, i * dim, dim);

      ncm_mset_fparams_set_vector_offset (self->mset, row, self->nadd_vals);

      ncm_mset_func_eval_vector (func, self->mset, x_v, qws);

      ncm_vector_free (qws);
    }

    for (i = 0; i < dim; i++)
    {
      gdouble *ret_i = ncm_vector_ptr (self->quantile_ws, i);

      gsl_sort (ret_i, dim, cat_len);
      ncm_matrix_set (res, i, 0, gsl_stats_mean (ret_i, dim, cat_len));
    }

    for (j = 0; j < p_val->len; j++)
    {
      const gdouble p = g_array_index (p_val, gdouble, j);

      g_assert_cmpfloat (p, >, 0.0);
      g_assert_cmpfloat (p, <, 1.0);

      for (i = 0; i < dim; i++)
      {
        gdouble *ret_i        = ncm_vector_ptr (self->quantile_ws, i);
        const gdouble lb_prob = (1.0 - p) / 2.0;
        const gdouble ub_prob = (1.0 + p) / 2.0;
        const gdouble lb      = gsl_stats_quantile_from_sorted_data (ret_i, dim, cat_len, lb_prob);
        const gdouble ub      = gsl_stats_quantile_from_sorted_data (ret_i, dim, cat_len, ub_prob);

        ncm_matrix_set (res, i, 1 + j * 2 + 0, lb);
        ncm_matrix_set (res, i, 1 + j * 2 + 1, ub);
      }
    }

    ncm_vector_clear (&self->quantile_ws);
    ncm_mset_fparams_set_vector (self->mset, save_params);
    ncm_vector_free (save_params);

    return res;
  }
}

/**
 * ncm_mset_catalog_calc_ci_interp:
 * @mcat: a #NcmMSetCatalog
 * @func: a #NcmMSetFunc of dimension one and one argument
 * @x_v: #NcmVector of arguments of @func
 * @p_val: (element-type double): probabilities of the confidence intervals, in $(0, 1)$
 * @nodes: unused
 * @mtype: #NcmFitRunMsgs log level
 *
 * The mean and the confidence intervals of @func over the rows, at each argument in
 * @x_v, with the quantiles taken from a #NcmStatsDist1dEPDF of the values at each
 * argument. Row $k$ of the result belongs to the $k$-th argument and has $1 + 4 n_p$
 * columns, $n_p$ the length of @p_val: column 0 holds the mean; columns $1 + 2j$ and
 * $2 + 2j$ the central interval of probability $p_j$, as in
 * ncm_mset_catalog_calc_ci_direct(); columns $1 + 2n_p + 2j$ and $2 + 2n_p + 2j$ the
 * one-sided bounds, the quantiles $p_j$ and $1 - p_j$.
 *
 * Row weights are not used. The free parameters of the #NcmMSet are restored at the end.
 *
 * Returns: (transfer full): a #NcmMatrix with the means and the interval bounds.
 */
NcmMatrix *
ncm_mset_catalog_calc_ci_interp (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmVector *x_v, GArray *p_val, guint nodes, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint dim             = ncm_vector_len (x_v);

  g_assert_cmpuint (p_val->len, >, 0);
  {
    const guint nelem      = p_val->len * 4 + 1;
    NcmMatrix *res         = ncm_matrix_new (dim, nelem);
    NcmVector *save_params = ncm_vector_new (ncm_mset_fparams_len (self->mset));
    const guint cat_len    = ncm_mset_catalog_len (mcat);
    const guint div        = cat_len > 100 ? cat_len / 100 : 1;
    GPtrArray *epdf_a      = g_ptr_array_sized_new (dim);
    guint i, j;

    ncm_mset_fparams_get_vector (self->mset, save_params);

    if (self->quantile_ws == NULL)
    {
      self->quantile_ws = ncm_vector_new (dim);
    }
    else if (ncm_vector_len (self->quantile_ws) != dim)
    {
      ncm_vector_clear (&self->quantile_ws);
      self->quantile_ws = ncm_vector_new (dim);
    }

    g_ptr_array_set_free_func (epdf_a, (GDestroyNotify) ncm_stats_dist1d_free);

    for (i = 0; i < dim; i++)
    {
      NcmStatsDist1dEPDF *epdf = ncm_stats_dist1d_epdf_new (NCM_MSET_CATALOG_DIST_EST_SD_SCALE);

      g_ptr_array_add (epdf_a, epdf);
    }

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_message ("# Calculating %u models in catalog: \n# - |", cat_len);

      for (i = 0; i < 100; i++)
        ncm_message ("-");

      ncm_message ("|\n# - |");
    }

    for (i = 0; i < cat_len; i++)
    {
      NcmVector *row = ncm_mset_catalog_peek_row (mcat, i);

      ncm_mset_fparams_set_vector_offset (self->mset, row, self->nadd_vals);

      ncm_mset_func_eval_vector (func, self->mset, x_v, self->quantile_ws);

      for (j = 0; j < dim; j++)
      {
        NcmStatsDist1dEPDF *epdf = g_ptr_array_index (epdf_a, j);

        ncm_stats_dist1d_epdf_add_obs (epdf, ncm_vector_get (self->quantile_ws, j));
      }

      if ((mtype > NCM_FIT_RUN_MSGS_NONE) && (i % div == 0))
        ncm_message ("=");
    }

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      if (i % div != 0)
        ncm_message ("=");

      ncm_message ("|\n");
      ncm_message ("# - |");

      for (i = 0; i < 100; i++)
        ncm_message ("-");

      ncm_message ("|\n");
    }

    for (i = 0; i < dim; i++)
    {
      NcmStatsDist1dEPDF *epdf = g_ptr_array_index (epdf_a, i);
      const gdouble mean       = ncm_stats_dist1d_epdf_get_obs_mean (epdf);

      ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf));

      ncm_matrix_set (res, i, 0, mean);

      for (j = 0; j < p_val->len; j++)
      {
        const gdouble p       = g_array_index (p_val, gdouble, j);
        const gdouble lb_prob = (1.0 - p) / 2.0;
        const gdouble ub_prob = (1.0 + p) / 2.0;
        const gdouble lb      = ncm_stats_dist1d_eval_inv_pdf (NCM_STATS_DIST1D (epdf), lb_prob);
        const gdouble ub      = ncm_stats_dist1d_eval_inv_pdf (NCM_STATS_DIST1D (epdf), ub_prob);

        g_assert_cmpfloat (p, >, 0.0);
        g_assert_cmpfloat (p, <, 1.0);

        ncm_matrix_set (res, i, 1 + j * 2 + 0, lb);
        ncm_matrix_set (res, i, 1 + j * 2 + 1, ub);
      }

      for (j = 0; j < p_val->len; j++)
      {
        const gdouble p  = g_array_index (p_val, gdouble, j);
        const gdouble lb = ncm_stats_dist1d_eval_inv_pdf (NCM_STATS_DIST1D (epdf), p);
        const gdouble ub = ncm_stats_dist1d_eval_inv_pdf_tail (NCM_STATS_DIST1D (epdf), p);

        g_assert_cmpfloat (p, >, 0.0);
        g_assert_cmpfloat (p, <, 1.0);

        ncm_matrix_set (res, i, 1 + 2 * p_val->len + j * 2 + 0, lb);
        ncm_matrix_set (res, i, 1 + 2 * p_val->len + j * 2 + 1, ub);
      }
    }

    g_ptr_array_unref (epdf_a);
    ncm_mset_fparams_set_vector (self->mset, save_params);
    ncm_vector_free (save_params);

    return res;
  }
}

/**
 * ncm_mset_catalog_calc_pvalue:
 * @mcat: a #NcmMSetCatalog
 * @func: a #NcmMSetFunc of dimension one and one argument
 * @x_v: #NcmVector of arguments of @func
 * @lim: (element-type double): limits of the p-values
 * @nodes: unused
 * @mtype: #NcmFitRunMsgs log level
 *
 * The lower tail probability of @func at each limit in @lim, at each argument in @x_v,
 * from a #NcmStatsDist1dEPDF of the values at each argument. Row $k$ of the result
 * belongs to the $k$-th argument and column $j$ to the $j$-th limit; a limit below the
 * range of the distribution gives zero and one above it gives one.
 *
 * Row weights are not used. The free parameters of the #NcmMSet are restored at the end.
 *
 * Returns: (transfer full): a #NcmMatrix with the p-values.
 */
NcmMatrix *
ncm_mset_catalog_calc_pvalue (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmVector *x_v, GArray *lim, guint nodes, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint dim             = ncm_vector_len (x_v);

  g_assert_cmpuint (lim->len, >, 0);
  {
    const guint nelem      = lim->len;
    NcmMatrix *res         = ncm_matrix_new (dim, nelem);
    NcmVector *save_params = ncm_vector_new (ncm_mset_fparams_len (self->mset));
    const guint cat_len    = ncm_mset_catalog_len (mcat);
    const guint div        = cat_len > 100 ? cat_len / 100 : 1;
    GPtrArray *epdf_a      = g_ptr_array_sized_new (dim);
    guint i, j;

    ncm_mset_fparams_get_vector (self->mset, save_params);

    if (self->quantile_ws == NULL)
    {
      self->quantile_ws = ncm_vector_new (dim);
    }
    else if (ncm_vector_len (self->quantile_ws) != dim)
    {
      ncm_vector_clear (&self->quantile_ws);
      self->quantile_ws = ncm_vector_new (dim);
    }

    g_ptr_array_set_free_func (epdf_a, (GDestroyNotify) ncm_stats_dist1d_free);

    for (i = 0; i < dim; i++)
    {
      NcmStatsDist1dEPDF *epdf = ncm_stats_dist1d_epdf_new (NCM_MSET_CATALOG_DIST_EST_SD_SCALE);

      g_ptr_array_add (epdf_a, epdf);
    }

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_message ("# Calculating %u models in catalog: \n# - |", cat_len);

      for (i = 0; i < 100; i++)
        ncm_message ("-");

      ncm_message ("|\n# - |");
    }

    for (i = 0; i < cat_len; i++)
    {
      NcmVector *row = ncm_mset_catalog_peek_row (mcat, i);

      ncm_mset_fparams_set_vector_offset (self->mset, row, self->nadd_vals);

      ncm_mset_func_eval_vector (func, self->mset, x_v, self->quantile_ws);

      for (j = 0; j < dim; j++)
      {
        NcmStatsDist1dEPDF *epdf = g_ptr_array_index (epdf_a, j);

        ncm_stats_dist1d_epdf_add_obs (epdf, ncm_vector_get (self->quantile_ws, j));
      }

      if ((mtype > NCM_FIT_RUN_MSGS_NONE) && (i % div == 0))
        ncm_message ("=");
    }

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      if (i % div != 0)
        ncm_message ("=");

      ncm_message ("|\n");
      ncm_message ("# - |");

      for (i = 0; i < 100; i++)
        ncm_message ("-");

      ncm_message ("|\n");
    }

    for (i = 0; i < dim; i++)
    {
      NcmStatsDist1dEPDF *epdf = g_ptr_array_index (epdf_a, i);
      NcmStatsDist1d *dist1d   = NCM_STATS_DIST1D (epdf);

      ncm_stats_dist1d_prepare (dist1d);

      for (j = 0; j < lim->len; j++)
      {
        const gdouble l = g_array_index (lim, gdouble, j);
        gdouble pvalue;

        if (l < ncm_stats_dist1d_get_xi (dist1d))
          pvalue = 0.0;

        else if (l > ncm_stats_dist1d_get_xf (dist1d))
          pvalue = 1.0;

        else
          pvalue = ncm_stats_dist1d_eval_pdf (dist1d, l);

        ncm_matrix_set (res, i, j, pvalue);
      }
    }

    g_ptr_array_unref (epdf_a);
    ncm_mset_fparams_set_vector (self->mset, save_params);
    ncm_vector_free (save_params);

    return res;
  }
}

/**
 * ncm_mset_catalog_calc_distrib:
 * @mcat: a #NcmMSetCatalog
 * @func: a #NcmMSetFunc of dimension one and no arguments
 * @mtype: #NcmFitRunMsgs log level
 *
 * The distribution of @func over the rows, a prepared #NcmStatsDist1dEPDF of its values.
 * Row weights are not used. The free parameters of the #NcmMSet are restored at the end.
 *
 * Returns: (transfer full): a #NcmStatsDist1d describing the distribution.
 */
NcmStatsDist1d *
ncm_mset_catalog_calc_distrib (NcmMSetCatalog *mcat, NcmMSetFunc *func, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint dim                   = ncm_mset_func_get_dim (func);

  g_assert_cmpuint (dim, ==, 1);
  {
    NcmStatsDist1dEPDF *epdf1d = ncm_stats_dist1d_epdf_new (NCM_MSET_CATALOG_DIST_EST_SD_SCALE);
    NcmVector *save_params     = ncm_vector_new (ncm_mset_fparams_len (self->mset));
    const guint cat_len        = ncm_mset_catalog_len (mcat);
    const guint div            = cat_len > 100 ? cat_len / 100 : 1;
    guint i;

    ncm_mset_fparams_get_vector (self->mset, save_params);

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_message ("# Calculating %u models in catalog: \n# - |", cat_len);

      for (i = 0; i < 100; i++)
        ncm_message ("-");

      ncm_message ("|\n# - |");
    }

    for (i = 0; i < cat_len; i++)
    {
      NcmVector *row = ncm_mset_catalog_peek_row (mcat, i);
      gdouble x;

      ncm_mset_fparams_set_vector_offset (self->mset, row, self->nadd_vals);
      x = ncm_mset_func_eval0 (func, self->mset);
      ncm_stats_dist1d_epdf_add_obs (epdf1d, x);

      if ((mtype > NCM_FIT_RUN_MSGS_NONE) && (i % div == 0))
        ncm_message ("=");
    }

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      if (i % div != 0)
        ncm_message ("=");

      ncm_message ("|\n");
      ncm_message ("# - |");

      for (i = 0; i < 100; i++)
        ncm_message ("-");

      ncm_message ("|\n");
    }

    ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf1d));

    ncm_mset_fparams_set_vector (self->mset, save_params);
    ncm_vector_free (save_params);

    return NCM_STATS_DIST1D (epdf1d);
  }
}

static NcmStatsDist1d *
_ncm_mset_catalog_calc_distrib (NcmMSetCatalog *mcat, guint vi, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmStatsDist1dEPDF *epdf1d  = ncm_stats_dist1d_epdf_new (NCM_MSET_CATALOG_DIST_EST_SD_SCALE);
  NcmVector *save_params      = ncm_vector_new (ncm_mset_fparams_len (self->mset));
  const guint cat_len         = ncm_mset_catalog_len (mcat);
  const guint div             = cat_len > 100 ? cat_len / 100 : 1;
  guint i;

  ncm_mset_fparams_get_vector (self->mset, save_params);

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_message ("# Calculating %u models in catalog: \n# - |", cat_len);

    for (i = 0; i < 100; i++)
      ncm_message ("-");

    ncm_message ("|\n# - |");
  }

  for (i = 0; i < cat_len; i++)
  {
    NcmVector *row = ncm_mset_catalog_peek_row (mcat, i);
    gdouble x      = ncm_vector_get (row, vi);

    ncm_stats_dist1d_epdf_add_obs (epdf1d, x);

    if ((mtype > NCM_FIT_RUN_MSGS_NONE) && (i % div == 0))
      ncm_message ("=");
  }

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (i % div != 0)
      ncm_message ("=");

    ncm_message ("|\n");
    ncm_message ("# - |");

    for (i = 0; i < 100; i++)
      ncm_message ("-");

    ncm_message ("|\n");
  }

  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf1d));

  ncm_mset_fparams_set_vector (self->mset, save_params);
  ncm_vector_free (save_params);

  return NCM_STATS_DIST1D (epdf1d);
}

/**
 * ncm_mset_catalog_calc_param_distrib:
 * @mcat: a #NcmMSetCatalog
 * @pi: a #NcmMSetPIndex of a free parameter
 * @mtype: #NcmFitRunMsgs log level
 *
 * The distribution of the free parameter @pi over the rows, as in
 * ncm_mset_catalog_calc_distrib().
 *
 * Returns: (transfer full): a #NcmStatsDist1d describing the distribution.
 */
NcmStatsDist1d *
ncm_mset_catalog_calc_param_distrib (NcmMSetCatalog *mcat, const NcmMSetPIndex *pi, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint fpi                   = ncm_mset_fparam_get_fpi (self->mset, pi->mid, pi->pid);

  g_assert (ncm_mset_param_get_ftype (self->mset, pi->mid, pi->pid, NULL) == NCM_PARAM_TYPE_FREE);

  return _ncm_mset_catalog_calc_distrib (mcat, fpi + self->nadd_vals, mtype);
}

/**
 * ncm_mset_catalog_calc_add_param_distrib:
 * @mcat: a #NcmMSetCatalog
 * @add_param: additional value index
 * @mtype: #NcmFitRunMsgs log level
 *
 * The distribution of the additional value @add_param over the rows, as in
 * ncm_mset_catalog_calc_distrib().
 *
 * Returns: (transfer full): a #NcmStatsDist1d describing the distribution.
 */
NcmStatsDist1d *
ncm_mset_catalog_calc_add_param_distrib (NcmMSetCatalog *mcat, guint add_param, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert_cmpuint (add_param, <, self->nadd_vals);

  return _ncm_mset_catalog_calc_distrib (mcat, add_param, mtype);
}

static void
_ncm_mset_catalog_calc_ensemble_evol (NcmMSetCatalog *mcat, guint vi, guint nsteps, NcmFitRunMsgs mtype, NcmVector **pval, NcmMatrix **t_evol)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmStatsDist1dEPDF *epdf1d  = ncm_stats_dist1d_epdf_new (NCM_MSET_CATALOG_DIST_EST_SD_SCALE);
  NcmStatsDist1d *sd1         = NCM_STATS_DIST1D (epdf1d);
  NcmVector *save_params      = ncm_vector_new (ncm_mset_fparams_len (self->mset));
  const guint max_t           = ncm_mset_catalog_max_time (mcat);
  NcmMatrix *res              = ncm_matrix_new (max_t, nsteps);
  NcmVector *pv               = ncm_vector_new (nsteps);
  const gdouble pmin          = ncm_vector_get (self->params_min, vi);
  const gdouble pmax          = ncm_vector_get (self->params_max, vi);
  const guint div             = max_t > 100 ? max_t / 100 : 1;

  guint i, t;

  *pval   = pv;
  *t_evol = res;

  ncm_mset_fparams_get_vector (self->mset, save_params);

  g_assert_cmpuint (self->nchains, >, 1);
  g_assert_cmpuint (nsteps, >, 1);

  for (i = 0; i < nsteps; i++)
  {
    const gdouble pv_i = pmin + (pmax - pmin) / (nsteps - 1.0) * i;

    ncm_vector_set (pv, i, pv_i);
  }

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_message ("# Calculating evolution to time %u: \n# - |", max_t);

    for (i = 0; i < 100; i++)
      ncm_message ("-");

    ncm_message ("|\n# - |");
  }

  for (t = 0; t < max_t; t++)
  {
    for (i = 0; i < self->nchains; i++)
    {
      NcmVector *row = ncm_mset_catalog_peek_row (mcat, t * self->nchains + i);
      gdouble x      = ncm_vector_get (row, vi);

      ncm_stats_dist1d_epdf_add_obs (epdf1d, x);
    }

    ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf1d));

    for (i = 0; i < nsteps; i++)
    {
      const gdouble pv_i = ncm_vector_get (pv, i);
      const gdouble p_i  = ncm_stats_dist1d_eval_p (sd1, pv_i);

      ncm_matrix_set (res, t, i, p_i);
    }

    ncm_stats_dist1d_epdf_reset (epdf1d);

    if ((mtype > NCM_FIT_RUN_MSGS_NONE) && (t % div == 0))
      ncm_message ("=");
  }

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (t % div != 0)
      ncm_message ("=");

    ncm_message ("|\n");
    ncm_message ("# - |");

    for (i = 0; i < 100; i++)
      ncm_message ("-");

    ncm_message ("|\n");
  }

  ncm_stats_dist1d_free (NCM_STATS_DIST1D (epdf1d));

  ncm_mset_fparams_set_vector (self->mset, save_params);
  ncm_vector_free (save_params);
}

/**
 * ncm_mset_catalog_calc_param_ensemble_evol:
 * @mcat: a #NcmMSetCatalog
 * @pi: a #NcmMSetPIndex of a free parameter
 * @nsteps: number of grid points, at least two
 * @mtype: #NcmFitRunMsgs log level
 * @pval: (out callee-allocates): the grid of parameter values
 * @t_evol: (out callee-allocates): the density at each iteration
 *
 * The distribution of the free parameter @pi across the ensemble at each iteration. @pval
 * is a uniform grid of @nsteps points over the sampled range of the parameter; row $t$ of
 * @t_evol holds, at those points, the normalized density of a #NcmStatsDist1dEPDF of the
 * $n_\mathrm{chains}$ values at iteration $t$. The catalog must have more than one chain.
 *
 */
void
ncm_mset_catalog_calc_param_ensemble_evol (NcmMSetCatalog *mcat, const NcmMSetPIndex *pi, guint nsteps, NcmFitRunMsgs mtype, NcmVector **pval, NcmMatrix **t_evol)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint fpi                   = ncm_mset_fparam_get_fpi (self->mset, pi->mid, pi->pid);

  g_assert (ncm_mset_param_get_ftype (self->mset, pi->mid, pi->pid, NULL) == NCM_PARAM_TYPE_FREE);

  _ncm_mset_catalog_calc_ensemble_evol (mcat, fpi + self->nadd_vals, nsteps, mtype, pval, t_evol);
}

/**
 * ncm_mset_catalog_calc_add_param_ensemble_evol:
 * @mcat: a #NcmMSetCatalog
 * @add_param: additional value index
 * @nsteps: number of grid points, at least two
 * @mtype: #NcmFitRunMsgs log level
 * @pval: (out callee-allocates): the grid of values
 * @t_evol: (out callee-allocates): the density at each iteration
 *
 * The distribution of the additional value @add_param across the ensemble at each
 * iteration, as in ncm_mset_catalog_calc_param_ensemble_evol().
 *
 */
void
ncm_mset_catalog_calc_add_param_ensemble_evol (NcmMSetCatalog *mcat, guint add_param, guint nsteps, NcmFitRunMsgs mtype, NcmVector **pval, NcmMatrix **t_evol)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert_cmpuint (add_param, <, self->nadd_vals);
  _ncm_mset_catalog_calc_ensemble_evol (mcat, add_param, nsteps, mtype, pval, t_evol);
}

/*
 * Replaces the rows of the catalog by rows, with the given Markovian id, written in one
 * recording named by the sampler string (the run type without one). A
 * FITS storage is rebuilt: the file is renamed to <file>.<n>.bak, with the first free n,
 * and a new file written under its name; a memory storage is erased. The description
 * and the RNG are kept, and the burn-in is zero after it.
 */
static void
_ncm_mset_catalog_rewrite (NcmMSetCatalog *mcat, GPtrArray *rows, const gint markovian_id)
{
  NcmMSetCatalogPrivate *self             = ncm_mset_catalog_get_instance_private (mcat);
  NcmMSetCatalogStorageFits *storage_fits = _ncm_mset_catalog_peek_fits_storage (mcat);
  GError *error                           = NULL;
  guint i;

  _ncm_mset_catalog_end_recording (mcat);

  if (storage_fits != NULL)
  {
    gchar *file = g_strdup (ncm_mset_catalog_storage_fits_peek_filename (storage_fits));
    guint bak_n = 0;

    ncm_mset_catalog_storage_clear (&self->storage);

    while (TRUE)
    {
      gchar *bak_name = g_strdup_printf ("%s.%d.bak", file, bak_n);

      if (!g_file_test (bak_name, G_FILE_TEST_EXISTS))
      {
        g_rename (file, bak_name);
        g_free (bak_name);

        break;
      }

      g_free (bak_name);
      bak_n++;
    }

    self->storage = _ncm_mset_catalog_new_fits_storage (mcat, file, FALSE);
    g_free (file);

    ncm_mset_catalog_storage_set_desc (self->storage, self->desc, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_rewrite");

    if (self->rng != NULL)
    {
      ncm_mset_catalog_storage_set_rng (self->storage, self->rng, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_rewrite");
    }
  }
  else
  {
    ncm_mset_catalog_storage_erase (self->storage, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_rewrite");
  }

  ncm_mset_catalog_reset_stats (mcat);

  self->burnin       = 0;
  self->cur_id       = -1;
  self->markovian_id = 0;

  if (rows->len > 0)
  {
    g_object_set (self->storage, "smode", NCM_MSET_CATALOG_SYNC_DISABLE, NULL);

    ncm_mset_catalog_storage_start_recording (self->storage,
                                              (self->sampler_str != NULL) ? self->sampler_str : self->rtype_str,
                                              self->sampler_opts_str, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_rewrite");

    for (i = 0; i < rows->len; i++)
    {
      NcmVector *row = g_ptr_array_index (rows, i);

      ncm_mset_catalog_storage_add_row (self->storage, row, &error);
      _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_rewrite");
      _ncm_mset_catalog_post_update_stats (mcat, row);
    }

    self->markovian_id = markovian_id;
    _ncm_mset_catalog_push_markovian_id (mcat);

    ncm_mset_catalog_storage_end_recording (self->storage, &error);
    _ncm_mset_catalog_check_error (&error, "_ncm_mset_catalog_rewrite");

    g_object_set (self->storage, "smode", self->smode, NULL);
  }
  else
  {
    self->markovian_id = markovian_id;
  }
}

/**
 * ncm_mset_catalog_trim:
 * @mcat: a #NcmMSetCatalog
 * @tc: number of iterations to drop
 * @thin: keep one iteration in every @thin
 *
 * Drops the first @tc iterations, the first $t_c n_\mathrm{chains}$ rows, and keeps one
 * iteration in every @thin of the rest, all chains of it. The markovian id follows the
 * rows it pointed to. The catalog
 * must hold complete ensembles and must not be read-only.
 *
 * When a file is attached it is renamed to `<file>.<n>.bak`, with the first free $n$,
 * and the trimmed catalog is written to a new file under the original name, in one
 * sampler segment. Nothing is done when @tc is zero and @thin is one.
 *
 */
void
ncm_mset_catalog_trim (NcmMSetCatalog *mcat, const guint tc, const guint thin)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  /* This rewrites the catalog in place: the file is renamed to <file>.<n>.bak and a fresh
   * one written, which goes around the read-only FITS handle opened by
   * ncm_mset_catalog_new_from_file_ro(). No analysis path needs that - the numcosmo CLI
   * resolves --burnin/--tail to a row count at read time and thins the analysis rather
   * than the file - so a read-only catalog reaching here is a programming error. The
   * superseded tools/mcat_analyze.c is the one caller that did it, silently destroying
   * the catalog it was asked only to read. */
  g_return_if_fail (!self->readonly);

  g_assert_cmpuint (thin, >=, 1);
  g_assert_cmpuint (tc,   >=, 0);

  if ((thin != 1) || (tc > 0))
  {
    GPtrArray *rows = _ncm_mset_catalog_collect_rows (mcat, 0, 1);
    GPtrArray *kept = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vector_free);
    guint t, tf;

    /* Iterations before the Markovian chain, in the old numbering; after the cut and the
     * thinning the same rows are the first ceil ((mark_iters - tc) / thin) kept iterations. */
    const gint mark_iters = (ncm_mset_catalog_get_markovian_id (mcat) + (gint) self->nchains - 1) / (gint) self->nchains;
    const gint kept_iters = (mark_iters > (gint) tc) ? (mark_iters - (gint) tc + (gint) thin - 1) / (gint) thin : 0;

    g_assert_cmpint (rows->len % self->nchains, ==, 0);

    tf = rows->len / self->nchains;

    for (t = tc; t < tf; t += thin)
    {
      guint i;

      for (i = 0; i < self->nchains; i++)
        g_ptr_array_add (kept, ncm_vector_ref (g_ptr_array_index (rows, t * self->nchains + i)));
    }

    _ncm_mset_catalog_rewrite (mcat, kept, MIN (kept_iters * (gint) self->nchains, (gint) kept->len));

    g_ptr_array_unref (kept);
    g_ptr_array_unref (rows);
  }
}

/**
 * ncm_mset_catalog_trim_p:
 * @mcat: a #NcmMSetCatalog
 * @p: fraction of the iterations to drop, in $(0, 1)$
 *
 * Drops the first fraction @p of the iterations with ncm_mset_catalog_trim().
 *
 */
void
ncm_mset_catalog_trim_p (NcmMSetCatalog *mcat, const gdouble p)
{
  g_assert_cmpfloat (p, >, 0.0);
  g_assert_cmpfloat (p, <, 1.0);

  ncm_mset_catalog_trim (mcat, ncm_mset_catalog_max_time (mcat) * p, 1);
}

/**
 * ncm_mset_catalog_trim_oob:
 * @mcat: a #NcmMSetCatalog
 * @out_file: (nullable): output filename, %NULL for a catalog in memory
 *
 * Removes the rows whose free parameters lie outside the bounds of the #NcmMSet of the
 * catalog. The catalog is renumbered as a single chain, with the markovian id after the
 * kept rows that preceded it, and written in one sampler segment
 * to @out_file, which must not hold rows; the storage used before is left as it was.
 *
 * Returns: the number of rows removed.
 */
guint
ncm_mset_catalog_trim_oob (NcmMSetCatalog *mcat, const gchar *out_file)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GPtrArray *rows             = _ncm_mset_catalog_collect_rows (mcat, 0, 1);
  const guint mark_rows       = ncm_mset_catalog_get_markovian_id (mcat);
  NcmMSetCatalogDesc *desc    = _ncm_mset_catalog_new_desc (mcat, 1, self->m2lnp_var, self->rtype_str);
  GError *error               = NULL;
  NcmMSetCatalogStorage *out;
  guint kept_before = 0;
  guint i, ndel;

  if (out_file != NULL)
    out = _ncm_mset_catalog_new_fits_storage (mcat, out_file, FALSE);
  else
    out = _ncm_mset_catalog_new_memory_storage (mcat);

  if (ncm_mset_catalog_storage_nrows (out, &error) > 0)
    g_error ("ncm_mset_catalog_trim_oob: file `%s' holds rows.", out_file);

  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_trim_oob");

  ncm_mset_catalog_storage_set_desc (out, desc, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_trim_oob");

  if (self->rng != NULL)
  {
    ncm_mset_catalog_storage_set_rng (out, self->rng, &error);
    _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_trim_oob");
  }

  _ncm_mset_catalog_end_recording (mcat);
  ncm_mset_catalog_storage_clear (&self->storage);
  self->storage = out;

  _ncm_mset_catalog_take_desc (mcat, desc);
  ncm_mset_catalog_desc_free (desc);

  self->burnin       = 0;
  self->cur_id       = -1;
  self->markovian_id = 0;

  g_object_set (self->storage, "smode", NCM_MSET_CATALOG_SYNC_DISABLE, NULL);

  ncm_mset_catalog_storage_start_recording (self->storage,
                                            (self->sampler_str != NULL) ? self->sampler_str : self->rtype_str,
                                            self->sampler_opts_str, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_trim_oob");

  ndel = 0;

  for (i = 0; i < rows->len; i++)
  {
    NcmVector *row_t = g_ptr_array_index (rows, i);

    if (ncm_mset_fparam_valid_bounds_offset (self->mset, row_t, self->nadd_vals))
    {
      ncm_mset_catalog_storage_add_row (self->storage, row_t, &error);
      _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_trim_oob");
      _ncm_mset_catalog_post_update_stats (mcat, row_t);

      if (i < mark_rows)
        kept_before++;
    }
    else
    {
      ndel++;
    }
  }

  /* Rows are renumbered: the Markovian chain starts after the kept rows that preceded it. */
  self->markovian_id = kept_before;
  _ncm_mset_catalog_push_markovian_id (mcat);

  ncm_mset_catalog_storage_end_recording (self->storage, &error);
  _ncm_mset_catalog_check_error (&error, "ncm_mset_catalog_trim_oob");

  g_object_set (self->storage, "smode", self->smode, NULL);

  g_ptr_array_unref (rows);

  return ndel;
}

/**
 * ncm_mset_catalog_remove_last_ensemble:
 * @mcat: a #NcmMSetCatalog
 *
 * Removes the last ensemble, the last $n_\mathrm{chains}$ rows, or the rows of the last
 * incomplete ensemble when there is one. The markovian id moves back to the row after
 * the new last one when it pointed past it. An empty catalog warns and is left as it is;
 * a read-only one is refused.
 *
 * When a file is attached it is renamed to `<file>.<n>.bak`, with the first free $n$,
 * and the catalog is written to a new file under the original name, in one sampler
 * segment.
 *
 */
void
ncm_mset_catalog_remove_last_ensemble (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  /* These three rewrite the catalog in place: the file is renamed to <file>.<n>.bak and a
   * fresh one written, which goes around the read-only FITS handle opened by
   * ncm_mset_catalog_new_from_file_ro(). Burn-in belongs at read time - the burnin argument
   * of the _ro constructor, or --burnin on catalog analyze and plot-corner - so a read-only
   * catalog reaching here is a programming error, not a request. */
  g_return_if_fail (!self->readonly);

  const gint last_t = ncm_mset_catalog_len (mcat);

  if (last_t > 0)
  {
    const guint nrows_to_remove = (last_t % self->nchains == 0) ? self->nchains : (last_t % self->nchains);
    const guint new_len         = last_t  - nrows_to_remove;
    const gint markovian_id     = ncm_mset_catalog_get_markovian_id (mcat);
    GPtrArray *rows             = _ncm_mset_catalog_collect_rows (mcat, 0, 1);

    g_assert_cmpuint (new_len, <, last_t);
    g_ptr_array_set_size (rows, new_len);

    /* The removed rows may have been the ones the chain was to start at: the id then sits
     * at the row after the new last one, i.e. the phase is still open. */
    _ncm_mset_catalog_rewrite (mcat, rows, MIN (markovian_id, (gint) new_len));

    g_ptr_array_unref (rows);
  }
  else
  {
    g_warning ("ncm_mset_catalog_remove_last_ensemble: called on an empty catalog.");
  }
}

/**
 * ncm_mset_catalog_calc_max_ess_time:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @max_ess: (out) (optional): the smallest effective sample size over the columns at $t_m$
 * @mtype: #NcmFitRunMsgs log level
 *
 * The starting time $t_m$ that maximizes the smallest effective sample size (ESS) over
 * the columns, see ncm_stats_rows_max_ess_time(). A catalog with one chain is read row
 * by row; one with several chains is read through the ensemble means, so $t_m$ is in
 * iterations. Fewer than ten rows (ensembles) give no estimate: $t_m$ is zero and
 * @max_ess zero.
 *
 * Returns: the starting time $t_m$.
 */
guint
ncm_mset_catalog_calc_max_ess_time (NcmMSetCatalog *mcat, const guint ntests, gdouble *max_ess, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GPtrArray *sample           = _ncm_mset_catalog_get_diag_series (mcat);
  const gint last_t           = sample->len;
  NcmVector *esss             = NULL;
  gdouble wp_ess              = 0.0;
  gint bindex                 = 0;
  guint wp_order              = 0;
  guint wp                    = 0;

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    ncm_message ("# NcmMSetCatalog: Calculating catalog effective sample size from chain %d => 0 using %u blocks:\n", last_t, ntests);
  }

  if (last_t < 10)
  {
    if (mtype > NCM_FIT_RUN_MSGS_NONE)
      ncm_message ("# NcmMSetCatalog: Catalog too small.\n");

    if (max_ess != NULL)
      *max_ess = 0.0;

    g_ptr_array_unref (sample);

    return 0;
  }

  ncm_stats_rows_set_rows (self->diag_rows, sample);
  esss = ncm_stats_rows_max_ess_time (self->diag_rows, ntests, &bindex, &wp, &wp_order, &wp_ess);
  ncm_stats_rows_reset (self->diag_rows);
  g_ptr_array_unref (sample);

  if (max_ess != NULL)
    *max_ess = wp_ess;

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    const gchar *name = NULL;

    if (wp < self->nadd_vals)
      name = g_ptr_array_index (self->add_vals_names, wp);
    else
      name = ncm_mset_fparam_full_name (self->mset, wp - self->nadd_vals);

    ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4u\n", bindex);

    if (self->nchains > 1)
    {
      ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u (%04u)\n", last_t, last_t * self->nchains);
      ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u (%04u)\n", last_t - bindex, (last_t - bindex) * self->nchains);
    }
    else
    {
      ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", last_t);
      ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", last_t - bindex);
    }

    ncm_message ("# NcmMSetCatalog: - worst parameter:          %s[%02u]\n", name, wp);
    ncm_message ("# NcmMSetCatalog: - worst parameter ess:      %-.2f\n", wp_ess);
    ncm_message ("# NcmMSetCatalog: - worst parameter ar order: %-4u\n", wp_order);
    ncm_vector_log_vals (esss, "# NcmMSetCatalog: - ess's:                    ", "%-.2f", TRUE);
  }

  ncm_vector_clear (&esss);

  return bindex;
}

static gdouble
_f_one_m_pval (gdouble v_i, guint i, gpointer user_data)
{
  return (1.0 - v_i) * 100.0;
}

/**
 * ncm_mset_catalog_calc_heidel_diag:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @pvalue: the required Schruben test p-value
 * @mtype: #NcmFitRunMsgs log level
 *
 * Applies the Heidelberger and Welch convergence diagnostic, see
 * ncm_stats_rows_heidel_diag(). A catalog with one chain is read row by row; one with
 * several chains is read through the ensemble means, so the time is in iterations. When
 * @pvalue is zero the required p-value is $1 - 0.95^{1/n}$, $n$ the number of free
 * parameters, so that all $n$ pass together with probability 0.95.
 *
 * Returns: the smallest starting time at which every column passes, or zero when none
 * does or the catalog has fewer than ten rows (ensembles).
 */
guint
ncm_mset_catalog_calc_heidel_diag (NcmMSetCatalog *mcat, const guint ntests, const gdouble pvalue, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GPtrArray *sample           = _ncm_mset_catalog_get_diag_series (mcat);
  const gdouble pvalue_lef    = (pvalue == 0.0) ? NCM_STATS_VEC_HEIDEL_PVAL_COR (0.05, ncm_mset_fparams_len (self->mset)) : pvalue;
  gint bindex                 = 0;
  guint wp = 0, wp_order = 0;
  gdouble wp_pvalue  = 0.0;
  const guint last_t = sample->len;
  NcmVector *pvals;

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    ncm_cfg_msg_sepa ();
    ncm_message ("# NcmMSetCatalog: Applying the Heidelberger and Welch's convergence diagnostic from chain %d => 0 using %u blocks:\n",
                 last_t, ntests);
  }

  if (last_t < 10)
  {
    if (mtype > NCM_FIT_RUN_MSGS_NONE)
      ncm_message ("# NcmMSetCatalog: Catalog too small.\n");

    g_ptr_array_unref (sample);

    return 0;
  }

  ncm_stats_rows_set_rows (self->diag_rows, sample);
  pvals = ncm_stats_rows_heidel_diag (self->diag_rows, ntests, pvalue_lef, &bindex, &wp, &wp_order, &wp_pvalue);
  ncm_stats_rows_reset (self->diag_rows);
  g_ptr_array_unref (sample);

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    const gchar *name = NULL;

    if (wp < self->nadd_vals)
      name = g_ptr_array_index (self->add_vals_names, wp);
    else
      name = ncm_mset_fparam_full_name (self->mset, wp - self->nadd_vals);

    if (bindex < 0)
    {
      if (ntests > 1)
        ncm_message ("# NcmMSetCatalog: - test failed at all points.\n");
      else
        ncm_message ("# NcmMSetCatalog: - test failed for the full catalog.\n");
    }
    else
    {
      ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4u\n", bindex);
    }

    ncm_message ("# NcmMSetCatalog: - worst parameter:          %s[%02u]\n", name, wp);
    ncm_message ("# NcmMSetCatalog: - worst parameter pvalue:   %5.2f%%\n", (1.0 - wp_pvalue) * 100.0);
    ncm_message ("# NcmMSetCatalog: - worst parameter ar order: %-4u\n", wp_order);
    ncm_message ("# NcmMSetCatalog: - target pvalue:            %5.2f%%\n", pvalue_lef * 100.0);
    ncm_vector_log_vals_func (pvals, "# NcmMSetCatalog: - pvalues:                  ", "%5.2f%%", &_f_one_m_pval, NULL);
  }

  ncm_vector_free (pvals);

  return (bindex >= 0) ? bindex : 0;
}

/**
 * ncm_mset_catalog_calc_const_break:
 * @mcat: a #NcmMSetCatalog
 * @p: column index
 * @mtype: #NcmFitRunMsgs log level
 *
 * The time $t_0$ from which column @p stays near its robust mean, see
 * ncm_stats_rows_estimate_const_break(). A catalog with one chain is read row by row;
 * one with several chains is read through the ensemble means, so $t_0$ is in iterations.
 *
 * Returns: $t_0$, rounded up.
 */
guint
ncm_mset_catalog_calc_const_break (NcmMSetCatalog *mcat, guint p, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GPtrArray *sample           = _ncm_mset_catalog_get_diag_series (mcat);
  guint tc, n;

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    if (self->nchains > 1)
      n = ncm_stats_vec_nitens (self->e_mean_stats);
    else
      n = ncm_stats_vec_nitens (self->pstats);

    ncm_cfg_msg_sepa ();
    ncm_message ("# NcmMSetCatalog: Computing the constant break point for parameter `%d', sample size `%u':\n", p, n);
  }

  ncm_stats_rows_set_rows (self->diag_rows, sample);
  tc = ceil (ncm_stats_rows_estimate_const_break (self->diag_rows, p));
  ncm_stats_rows_reset (self->diag_rows);
  g_ptr_array_unref (sample);

  if (mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_message ("# NcmMSetCatalog: Constant break point at `%d':\n", tc);

  return tc;
}

/**
 * ncm_mset_catalog_trim_by_type:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @trim_type: the criteria, a #NcmMSetCatalogTrimType
 * @mtype: #NcmFitRunMsgs log level
 *
 * Trims the catalog with ncm_mset_catalog_trim() at the largest of the times given by
 * the criteria in @trim_type: ncm_mset_catalog_calc_max_ess_time(),
 * ncm_mset_catalog_calc_heidel_diag() with the default p-value, and
 * ncm_mset_catalog_calc_const_break() on the $-2\ln(L)$ column.
 *
 */
void
ncm_mset_catalog_trim_by_type (NcmMSetCatalog *mcat, const guint ntests, NcmMSetCatalogTrimType trim_type, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint t_c                   = 0;

  if (trim_type & NCM_MSET_CATALOG_TRIM_TYPE_ESS)
  {
    gdouble max_ess_time = 0.0;
    const guint t_c_ess  = ncm_mset_catalog_calc_max_ess_time (mcat, ntests, &max_ess_time, mtype);

    t_c = GSL_MAX (t_c, t_c_ess);
  }

  if (trim_type & NCM_MSET_CATALOG_TRIM_TYPE_HEIDEL)
  {
    const guint t_c_heidel = ncm_mset_catalog_calc_heidel_diag (mcat, ntests, 0.0, mtype);

    t_c = GSL_MAX (t_c, t_c_heidel);
  }

  if (trim_type & NCM_MSET_CATALOG_TRIM_TYPE_CK)
  {
    guint t_c_cb = ncm_mset_catalog_calc_const_break (mcat, self->m2lnp_var, mtype);

    t_c = GSL_MAX (t_c, t_c_cb);
  }

  if ((t_c > 0) && (mtype >= NCM_FIT_RUN_MSGS_SIMPLE))
    g_message ("# NcmMSetCatalog: Trimming the catalog at: %u.\n", t_c);

  ncm_mset_catalog_trim (mcat, t_c, 1);
}

typedef struct _NcmMSetCatalogESSRes
{
  guint chain_n;
  gint best_cutoff;
  guint size;
  guint size_left;
  guint wp;
  guint wp_order;
  gdouble wp_ess;
} NcmMSetCatalogESSRes;

static gint
_ess_res_cmp (gconstpointer a, gconstpointer b)
{
  const NcmMSetCatalogESSRes *res_a = a;
  const NcmMSetCatalogESSRes *res_b = b;

  if (res_a->best_cutoff < 0)
  {
    if (res_b->best_cutoff < 0)
    {
      if (res_a->wp_ess > res_b->wp_ess)
        return -1;
      else if (res_a->wp_ess == res_b->wp_ess)
        return 0;
      else
        return 1;
    }
    else
    {
      return -1;
    }
  }
  else if (res_b->best_cutoff < 0)
  {
    return 1;
  }

  if (res_a->best_cutoff > res_b->best_cutoff)
    return -1;
  else if (res_a->best_cutoff == res_b->best_cutoff)
    return 0;
  else
    return 1;
}

static void _ncm_mset_catalog_set_chain_rows (NcmMSetCatalog *mcat, const guint i);

/**
 * ncm_mset_catalog_max_ess_time_by_chain:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @max_ess: (out) (optional): the smallest effective sample size of the worst chain
 * @mtype: #NcmFitRunMsgs log level
 *
 * Applies ncm_stats_rows_max_ess_time() to each chain on its own. The worst chain is the
 * one with the latest starting time, and its time and ESS are returned. A catalog with
 * one chain is passed to ncm_mset_catalog_calc_max_ess_time(). Fewer than ten iterations
 * give zero, and @max_ess zero.
 *
 * Returns: the latest of the per-chain starting times, in iterations.
 */
guint
ncm_mset_catalog_max_ess_time_by_chain (NcmMSetCatalog *mcat, const guint ntests, gdouble *max_ess, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->nchains == 1)
  {
    return ncm_mset_catalog_calc_max_ess_time (mcat, ntests, max_ess, mtype);
  }
  else
  {
    gint tb_index = -1;
    guint twp = 0, twp_order = 0;
    gdouble twp_ess = 0.0;
    guint i, ti = 0;
    GArray *res_a;

    if (ncm_mset_catalog_max_time (mcat) < 10)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_message ("# NcmMSetCatalog: Catalog too small.\n");

      if (max_ess != NULL)
        max_ess[0] = 0.0;

      return 0;
    }

    res_a = g_array_new (FALSE, FALSE, sizeof (NcmMSetCatalogESSRes));

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      ncm_message ("# NcmMSetCatalog: Computing max ESS for all %u chains in the catalog:\n", self->nchains);
    }

    for (i = 0; i < self->nchains; i++)
    {
      gint bindex = 0;
      guint wp = 0, wp_order = 0;
      gdouble wp_ess = 0.0;
      NcmVector *esss;
      guint size;

      _ncm_mset_catalog_set_chain_rows (mcat, i);
      esss = ncm_stats_rows_max_ess_time (self->diag_rows, ntests, &bindex, &wp, &wp_order, &wp_ess);
      size = ncm_stats_rows_get_nrows (self->diag_rows);

      if (mtype > NCM_FIT_RUN_MSGS_SIMPLE)
      {
        NcmMSetCatalogESSRes res_i;

        res_i.chain_n     = i;
        res_i.best_cutoff = bindex;
        res_i.size        = size;
        res_i.size_left   = size - bindex;
        res_i.wp          = wp;
        res_i.wp_order    = wp_order;
        res_i.wp_ess      = wp_ess;

        g_array_append_val (res_a, res_i);
      }

      if (bindex > tb_index)
      {
        tb_index  = bindex;
        twp       = wp;
        twp_order = wp_order;
        twp_ess   = wp_ess;
        ti        = i;
      }

      ncm_vector_free (esss);
    }

    ncm_stats_rows_reset (self->diag_rows);

    if (mtype > NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_array_sort (res_a, &_ess_res_cmp);

      ncm_message ("# NcmMSetCatalog: Chains results from worst to best:\n");

      for (i = 0; i < self->nchains; i++)
      {
        NcmMSetCatalogESSRes *res_i = &g_array_index (res_a, NcmMSetCatalogESSRes, i);

        ncm_cfg_msg_sepa ();
        ncm_message ("# NcmMSetCatalog: - chain index:              %-4u\n", res_i->chain_n);
        ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4u\n", res_i->best_cutoff);
        ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", res_i->size);
        ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", res_i->size_left);
        ncm_message ("# NcmMSetCatalog: - worst parameter:          %-4u\n", res_i->wp);
        ncm_message ("# NcmMSetCatalog: - worst parameter order:    %-4u\n", res_i->wp_order);
        ncm_message ("# NcmMSetCatalog: - worst parameter ESS:      %-.2f\n", res_i->wp_ess);
      }
    }

    g_array_unref (res_a);

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      guint size = ncm_stats_vec_nitens (g_ptr_array_index (self->chain_pstats, ti));

      ncm_cfg_msg_sepa ();
      ncm_message ("# NcmMSetCatalog: - Worst chain:\n");
      ncm_cfg_msg_sepa ();
      ncm_message ("# NcmMSetCatalog: - worst chain no:           %-4u\n", ti);
      ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4d\n", tb_index);
      ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", size);
      ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", size - tb_index);
      ncm_message ("# NcmMSetCatalog: - worst parameter:          %-4u\n", twp);
      ncm_message ("# NcmMSetCatalog: - worst parameter order:    %-4u\n", twp_order);
      ncm_message ("# NcmMSetCatalog: - worst parameter ESS:      %-.2f\n", twp_ess);
    }

    if (max_ess != NULL)
      max_ess[0] = twp_ess;

    return tb_index;
  }
}

/*
 * Gives self->diag_rows the rows of chain i, rows i, i + nchains, i + 2 nchains, ..., in
 * order; the caller resets self->diag_rows after its last chain. Row t belongs to chain
 * t mod nchains.
 */
static void
_ncm_mset_catalog_set_chain_rows (NcmMSetCatalog *mcat, const guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GPtrArray *chain_rows       = _ncm_mset_catalog_collect_rows (mcat, i, self->nchains);

  ncm_stats_rows_set_rows (self->diag_rows, chain_rows);
  g_ptr_array_unref (chain_rows);
}

/**
 * ncm_mset_catalog_heidel_diag_by_chain:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @pvalue: required p-value, the default of ncm_mset_catalog_calc_heidel_diag() when zero
 * @wp_pvalue: (out) (optional): the worst value of the worst chain
 * @mtype: #NcmFitRunMsgs log level
 *
 * Applies ncm_stats_rows_heidel_diag() to each chain on its own. The worst chain is one
 * that passes from no starting time, the one with the largest value among those, or
 * else the one with the latest starting time; @wp_pvalue is its largest Cramér-von Mises
 * cumulative distribution value, one minus its p-value. Fewer than ten iterations give
 * zero, and @wp_pvalue zero.
 *
 * Returns: the smallest time from which every chain passes, in iterations, or zero when
 * a chain passes from none.
 */
guint
ncm_mset_catalog_heidel_diag_by_chain (NcmMSetCatalog *mcat, const guint ntests, const gdouble pvalue, gdouble *wp_pvalue, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  {
    const gdouble pvalue_lef = (pvalue == 0.0) ? NCM_STATS_VEC_HEIDEL_PVAL_COR (0.05, ncm_mset_fparams_len (self->mset)) : pvalue;
    gint tb_index            = -1;
    guint twp = 0, twp_order = 0;
    gdouble twp_pvalue = 0.0;
    guint i, ti = 0;
    GArray *res_a;

    if (ncm_mset_catalog_max_time (mcat) < 10)
    {
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_message ("# NcmMSetCatalog: Catalog too small.\n");

      if (wp_pvalue != NULL)
        wp_pvalue[0] = 0.0;

      return 0;
    }

    res_a = g_array_new (FALSE, FALSE, sizeof (NcmMSetCatalogESSRes));

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      ncm_message ("# NcmMSetCatalog: Computing Heidelberger and Welch's convergence diagnostic for all %u chains in the catalog:\n", self->nchains);
    }

    for (i = 0; i < self->nchains; i++)
    {
      gint bindex = 0;
      guint wp = 0, wp_order = 0;
      gdouble lwp_pvalue = 0.0;
      NcmVector *pvals;
      guint size;

      _ncm_mset_catalog_set_chain_rows (mcat, i);
      pvals = ncm_stats_rows_heidel_diag (self->diag_rows, ntests, pvalue_lef, &bindex, &wp, &wp_order, &lwp_pvalue);
      size  = ncm_stats_rows_get_nrows (self->diag_rows);

      if (mtype > NCM_FIT_RUN_MSGS_SIMPLE)
      {
        NcmMSetCatalogESSRes res_i;

        res_i.chain_n     = i;
        res_i.best_cutoff = bindex;
        res_i.size        = size;
        res_i.size_left   = size - bindex;
        res_i.wp          = wp;
        res_i.wp_order    = wp_order;
        res_i.wp_ess      = lwp_pvalue;

        g_array_append_val (res_a, res_i);
      }

      if ((i == 0) || ((tb_index >= 0) && (bindex > tb_index)) || ((bindex == -1) && (tb_index >= 0)) || ((bindex == -1) && (lwp_pvalue > twp_pvalue)))
      {
        tb_index   = bindex;
        twp        = wp;
        twp_order  = wp_order;
        twp_pvalue = lwp_pvalue;
        ti         = i;
      }

      ncm_vector_free (pvals);
    }

    ncm_stats_rows_reset (self->diag_rows);

    if (mtype > NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_array_sort (res_a, &_ess_res_cmp);

      ncm_message ("# NcmMSetCatalog: Chains results from worst to best:\n");

      for (i = 0; i < self->nchains; i++)
      {
        NcmMSetCatalogESSRes *res_i = &g_array_index (res_a, NcmMSetCatalogESSRes, i);

        ncm_cfg_msg_sepa ();
        ncm_message ("# NcmMSetCatalog: - chain index:              %-4u\n", res_i->chain_n);

        if (res_i->best_cutoff == -1)
        {
          ncm_message ("# NcmMSetCatalog: - best cutoff time:         not satisfied!\n");
          ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", res_i->size);
          ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", 0);
        }
        else
        {
          ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4u\n", res_i->best_cutoff);
          ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", res_i->size);
          ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", res_i->size_left);
        }

        ncm_message ("# NcmMSetCatalog: - worst parameter:          %-4u\n", res_i->wp);
        ncm_message ("# NcmMSetCatalog: - worst parameter order:    %-4u\n", res_i->wp_order);
        ncm_message ("# NcmMSetCatalog: - worst parameter pvalue:   %6.2f%%\n", (1.0 - res_i->wp_ess) * 100.0);
      }
    }

    g_array_unref (res_a);

    if (mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      guint size = ncm_stats_vec_nitens ((self->nchains == 1) ? self->pstats : g_ptr_array_index (self->chain_pstats, ti));

      ncm_cfg_msg_sepa ();
      ncm_message ("# NcmMSetCatalog: - Worst chain:\n");
      ncm_cfg_msg_sepa ();
      ncm_message ("# NcmMSetCatalog: - worst chain index:        %-4u\n", ti);

      if (tb_index == -1)
      {
        ncm_message ("# NcmMSetCatalog: - best cutoff time:         not satisfied!\n");
        ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", size);
        ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", 0);
      }
      else
      {
        ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4u\n", tb_index);
        ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", size);
        ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", size - tb_index);
      }

      ncm_message ("# NcmMSetCatalog: - worst parameter:          %-4u\n", twp);
      ncm_message ("# NcmMSetCatalog: - worst parameter order:    %-4u\n", twp_order);
      ncm_message ("# NcmMSetCatalog: - worst parameter pvalue:   %6.2f%%\n", (1.0 - twp_pvalue) * 100.0);
    }

    if (wp_pvalue != NULL)
      wp_pvalue[0] = twp_pvalue;

    return (tb_index >= 0) ? tb_index : 0;
  }
}


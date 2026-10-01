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
 * Ordered table of free-parameter values of a #NcmMSet, optionally kept in sync with a
 * FITS file.
 *
 * Each row holds #NcmMSetCatalog:nadd-vals additional values ($-2\ln L$ and derived
 * functions, named by #NcmMSetCatalog:nadd-val-names; for a #NcmMSetCatalog:weighted
 * catalog the last one is the row weight) followed by the free parameters. With
 * #NcmMSetCatalog:nchains chains the rows are interleaved, the row with id $i$ belonging to chain
 * $i \bmod n_\mathrm{chains}$ (the walkers of an ensemble sampler), and the catalog
 * keeps statistics of every chain, of the chain means and of each ensemble. Row ids start
 * at the first id (ncm_mset_catalog_get_first_id()); #NcmMSetCatalog:markovian-id marks
 * where the Markov chain starts and #NcmMSetCatalog:burnin drops leading rows when
 * reading a file.
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
#include "ncm/data/ncm_data_gauss_cov_mvnd.h"
#include "ncm/model/ncm_model_mvnd.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_func_eval.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_cdf.h>
#include <math.h>
#ifdef HAVE_CFITSIO
#include <fitsio.h>
#endif /* HAVE_CFITSIO */
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmMSetCatalogPrivate
{
  NcmMSet *mset;
  NcmObjArray *functions;
  guint nadd_vals;
  gint m2lnp_var;
  NcmVector *bestfit_row;
  gdouble bestfit;
  gdouble post_lnnorm;
  gdouble post_lnnorm_sd;
  gboolean post_lnnorm_up;
  GPtrArray *order_cat;
  gboolean order_cat_sort;
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
  NcmVector *tau;
  gchar *rng_inis;
  gchar *rng_stat;
  GTimer *sync_timer;
  gdouble sync_interval;
  gchar *file;
  gchar *mset_file;
  gchar *rtype_str;
  gchar *sampler_str;
  gchar *sampler_opts_str;
  gchar *init_sampler_str;
  GArray *porder;
  NcmVector *quantile_ws;
  gint first_id;
  gint cur_id;
  gint file_first_id;
  gint file_cur_id;
  gint markovian_id;
  gint file_markovian_id;
  glong burnin;

#ifdef HAVE_CFITSIO
  fitsfile *fptr;

#endif /* HAVE_CFITSIO */
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
  PROP_RNG,
  PROP_FILE,
  PROP_RUN_TYPE_STR,
  PROP_SAMPLER,
  PROP_SAMPLER_OPTS,
  PROP_INIT_SAMPLER,
  PROP_SYNC_MODE,
  PROP_SYNC_INTERVAL,
  PROP_READONLY,
};

static gint
_ncm_mset_catalog_double_compare (gconstpointer a, gconstpointer b, gpointer data)
{
  const NcmVector **row_a = (const NcmVector **) a;
  const NcmVector **row_b = (const NcmVector **) b;
  const gint i            = GPOINTER_TO_INT (data);
  const gdouble a_m2lnp   = ncm_vector_get (*row_a, i);
  const gdouble b_m2lnp   = ncm_vector_get (*row_b, i);

  return (a_m2lnp == b_m2lnp) ? 0.0 : ((a_m2lnp > b_m2lnp) ? +1 : -1);
}

static void
ncm_mset_catalog_init (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->mset           = NULL;
  self->functions      = NULL;
  self->nadd_vals      = 0;
  self->m2lnp_var      = -1;
  self->bestfit_row    = NULL;
  self->bestfit        = GSL_POSINF;
  self->post_lnnorm    = 0.0;
  self->post_lnnorm_sd = GSL_NAN;
  self->post_lnnorm_up = FALSE;
  self->order_cat      = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vector_free);
  self->order_cat_sort = FALSE;
  self->add_vals_names = g_ptr_array_new_with_free_func (g_free);
  self->add_vals_symbs = g_ptr_array_new_with_free_func (g_free);
  self->pstats         = NULL;
  self->smode          = NCM_MSET_CATALOG_SYNC_LEN;
  self->readonly       = FALSE;
  self->rng            = NULL;
  self->weighted       = FALSE;
  self->first_flush    = FALSE;
  self->nchains        = 0;
  self->chain_pstats   = g_ptr_array_new ();
  g_ptr_array_set_free_func (self->chain_pstats, (GDestroyNotify) & ncm_stats_vec_free);
  self->mean_pstats = NULL;
  self->e_var_array = g_ptr_array_new ();
  g_ptr_array_set_free_func (self->e_var_array, (GDestroyNotify) & ncm_vector_free);
  self->e_stats      = NULL;
  self->e_mean_stats = NULL;
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
  self->tau          = NULL;

  self->rng_inis          = NULL;
  self->rng_stat          = NULL;
  self->sync_timer        = g_timer_new ();
  self->cur_id            = -1; /* Represents that there are no elements in the catalog, i.e., the id of the last added row. */
  self->first_id          = 0;  /* The element to be in the catalog will be the one with index == 0, cross catalog index */
  self->file_cur_id       = -1; /* Represents that no elements in the catalog file, i.e., the id of the last added row. */
  self->file_first_id     = 0;  /* The element to be in the catalog file will be the one with index == 0, cross catalog index */
  self->markovian_id      = 0;  /* First row of the Markovian chain, cross catalog index; effective value is max (markovian_id, first_id) */
  self->file_markovian_id = -1; /* Value last written to the file, -1: never written */
  self->burnin            = 0;  /* Number of elements to ignore when reading a catalog */
  self->file              = NULL;
  self->mset_file         = NULL;
  self->rtype_str         = NULL;
  self->sampler_str       = NULL;
  self->sampler_opts_str  = NULL;
  self->init_sampler_str  = NULL;
  self->porder            = g_array_new (FALSE, FALSE, sizeof (gint));
  self->quantile_ws       = NULL;
#ifdef HAVE_CFITSIO
  self->fptr = NULL;
#endif /* HAVE_CFITSIO */
  self->pdf_i      = -1;
  self->h          = NULL;
  self->h_pdf      = NULL;
  self->params_max = NULL;
  self->params_min = NULL;

  self->constructed = FALSE;
}

#ifdef HAVE_CFITSIO
static void _ncm_mset_catalog_open_create_file (NcmMSetCatalog *mcat, gboolean load_from_cat);
static void _ncm_mset_catalog_flush_file (NcmMSetCatalog *mcat);
static NcmMSet *_ncm_mset_catalog_load_mset (NcmMSetCatalog *mcat);

#endif /* HAVE_CFITSIO */

static void
_ncm_mset_catalog_constructed_alloc_chains (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint free_params_len = ncm_mset_fparams_len (self->mset);
  const guint total           = free_params_len + self->nadd_vals + (self->weighted ? 1 : 0);
  guint i;

  self->pstats     = ncm_stats_vec_new (total, NCM_STATS_VEC_COV, TRUE);
  self->params_max = ncm_vector_new (total);
  self->params_min = ncm_vector_new (total);

  ncm_vector_set_all (self->params_max, GSL_NEGINF);
  ncm_vector_set_all (self->params_min, GSL_POSINF);

  if (self->nchains > 1)
  {
    for (i = 0; i < self->nchains; i++)
    {
      NcmStatsVec *pstats = ncm_stats_vec_new (total, NCM_STATS_VEC_COV, TRUE);

      g_ptr_array_add (self->chain_pstats, pstats);
    }

    self->mean_pstats  = ncm_stats_vec_new (free_params_len, NCM_STATS_VEC_COV, FALSE);
    self->e_stats      = ncm_stats_vec_new (total, NCM_STATS_VEC_VAR, FALSE);
    self->e_mean_stats = ncm_stats_vec_new (total, NCM_STATS_VEC_VAR, TRUE);
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

static void
_ncm_mset_catalog_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->constructed (object);
  {
    NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
    NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

    g_assert_cmpuint (self->add_vals_names->len, ==, self->add_vals_symbs->len);
    g_assert_cmpuint (self->add_vals_names->len, ==, self->nadd_vals);

    self->constructed = TRUE;

    if (self->file != NULL)
    {
      gchar *file = self->file;

      self->file = NULL;

      ncm_mset_catalog_set_file (mcat, file);

      g_free (file);
    }

    if (self->mset == NULL)
    {
#ifdef HAVE_CFITSIO

      if (!g_file_test (self->file, G_FILE_TEST_EXISTS))
        g_error ("_ncm_mset_catalog_constructed: cannot create catalog file `%s' not found.",
                 self->file);

      self->mset = _ncm_mset_catalog_load_mset (mcat);

      _ncm_mset_catalog_open_create_file (mcat, TRUE);
      _ncm_mset_catalog_constructed_alloc_chains (mcat);

      ncm_mset_catalog_sync (mcat, TRUE);
#else
      g_error ("_ncm_mset_catalog_constructed: cannot create catalog without mset.");
#endif /* HAVE_CFITSIO */
    }
    else
    {
      const guint free_params_len = ncm_mset_fparams_len (self->mset);
      const guint total           = free_params_len + self->nadd_vals + (self->weighted ? 1 : 0);

      g_array_set_size (self->porder, total);
      _ncm_mset_catalog_constructed_alloc_chains (mcat);

      if (self->weighted)
      {
        g_ptr_array_add (self->add_vals_names, g_strdup ("NcmMSetCatalog:Row-weights"));
        g_ptr_array_add (self->add_vals_symbs, g_strdup ("W"));
        self->nadd_vals++;
      }

      if (self->file != NULL)
      {
#ifdef HAVE_CFITSIO
        _ncm_mset_catalog_open_create_file (mcat, FALSE);
        ncm_mset_catalog_sync (mcat, TRUE);
#else
        g_error ("_ncm_mset_catalog_constructed: cannot create catalog without cfitsio.");
#endif /* HAVE_CFITSIO */
      }
    }
  }
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
    case PROP_READONLY:
      self->readonly = g_value_get_boolean (value);
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
    case PROP_RNG:
      g_value_set_object (value, self->rng);
      break;
    case PROP_FILE:
      g_value_set_string (value, self->file);
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
    case PROP_READONLY:
      g_value_set_boolean (value, self->readonly);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_dispose (GObject *object)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  ncm_vector_clear (&self->bestfit_row);
  g_clear_pointer (&self->order_cat, g_ptr_array_unref);

  ncm_mset_clear (&self->mset);
  ncm_obj_array_clear (&self->functions);
  ncm_rng_clear (&self->rng);
  ncm_stats_vec_clear (&self->pstats);
  ncm_vector_clear (&self->params_max);
  ncm_vector_clear (&self->params_min);

  g_clear_pointer (&self->chain_pstats, g_ptr_array_unref);
  ncm_stats_vec_clear (&self->mean_pstats);
  ncm_stats_vec_clear (&self->e_stats);
  ncm_stats_vec_clear (&self->e_mean_stats);

  g_clear_pointer (&self->e_var_array, g_ptr_array_unref);
  g_clear_pointer (&self->accept_ratio, g_array_unref);

  ncm_vector_clear (&self->chain_means);
  ncm_vector_clear (&self->chain_vars);
  ncm_matrix_clear (&self->chain_cov);
  ncm_matrix_clear (&self->chain_sM);
  g_clear_pointer (&self->chain_sM_ws, gsl_eigen_nonsymm_free);
  g_clear_pointer (&self->chain_sM_ev, gsl_vector_complex_free);
  ncm_stats_acorr_clear (&self->acorr);
  ncm_vector_clear (&self->tau);
  ncm_vector_clear (&self->quantile_ws);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_parent_class)->dispose (object);
}

#ifdef HAVE_CFITSIO
static void _ncm_mset_catalog_close_file (NcmMSetCatalog *mcat);

#endif /* HAVE_CFITSIO */

static void
_ncm_mset_catalog_finalize (GObject *object)
{
  NcmMSetCatalog *mcat        = NCM_MSET_CATALOG (object);
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->h != NULL)
    gsl_histogram_free (self->h);

  if (self->h_pdf != NULL)
    gsl_histogram_pdf_free (self->h_pdf);

#ifdef HAVE_CFITSIO
  _ncm_mset_catalog_close_file (mcat);
#endif /* HAVE_CFITSIO */

  g_clear_pointer (&self->rtype_str, g_free);
  g_clear_pointer (&self->sampler_str, g_free);
  g_clear_pointer (&self->sampler_opts_str, g_free);
  g_clear_pointer (&self->init_sampler_str, g_free);

  g_array_unref (self->porder);
  g_timer_destroy (self->sync_timer);

  g_ptr_array_unref (self->add_vals_names);
  g_ptr_array_unref (self->add_vals_symbs);

  g_clear_pointer (&self->rng_inis, g_free);
  g_clear_pointer (&self->rng_stat, g_free);

  g_clear_pointer (&self->file, g_free);
  g_clear_pointer (&self->mset_file, g_free);

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
   * the #NcmMSetCatalog:nchains rows before it. Same numbering as the first id
   * (ncm_mset_catalog_get_first_id()) and the current id; the effective value is never below the first id, and a file
   * without the key reads as its first id (every row Markovian). When a file is read with
   * a burn-in the stored value is shifted by it, like every other row id.
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
                                   PROP_RNG,
                                   g_param_spec_object ("rng",
                                                        NULL,
                                                        "Random number generator object",
                                                        NCM_TYPE_RNG,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SYNC_MODE,
                                   g_param_spec_enum ("smode",
                                                      NULL,
                                                      "Catalog sync mode",
                                                      NCM_TYPE_MSET_CATALOG_SYNC, NCM_MSET_CATALOG_SYNC_AUTO,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
   * A description of the sampler that filled the rows, for the record only: it is
   * written to the file's SAMPLER header key and, at every change, to a HISTORY card,
   * and it is never compared. A catalog may be continued with any sampler.
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
   * The sampler's tunable settings, as a colon separated list of `name=value' pairs. Kept
   * like #NcmMSetCatalog:sampler (the SAMPOPT header key and a HISTORY card, never
   * compared) and separate from it, so that the structure and its tuning can be read
   * apart: two runs of the same sampler differ here and nowhere else.
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
   * A description of the sampler that drew the initial points, kept like
   * #NcmMSetCatalog:sampler: the INITSMP header key and a HISTORY card, never compared.
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
  g_object_class_install_property (object_class,
                                   PROP_READONLY,
                                   g_param_spec_boolean ("read-only",
                                                         NULL,
                                                         "Whether the FITS file is opened read-only",
                                                         FALSE,
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
 * while the catalog is in use. With more than one chain the statistics of each chain
 * are kept separately.
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
 * Creates a #NcmMSetCatalog from the catalog in @filename, synchronized with it. The
 * model set is read from the file, or for older files from the file with the same base
 * name and the .mset extension.
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
 * ncm_mset_catalog_peek_info_from_file:
 * @filename: filename of the catalog fits
 * @nrows: (out): number of rows (samples) currently in the catalog
 * @nchains: (out): number of chains (walkers) in the catalog
 * @first_id: (out): id of the first element in the catalog
 *
 * Reads the number of rows, of chains and the first id of @filename from its header,
 * without loading the catalog, e.g. to convert a burn-in from iterations to rows (see
 * ncm_mset_catalog_set_burnin()) before ncm_mset_catalog_new_from_file_ro().
 *
 */
void
ncm_mset_catalog_peek_info_from_file (const gchar *filename, glong *nrows, guint *nchains, gint *first_id)
{
#ifdef HAVE_CFITSIO
  fitsfile *fptr = NULL;
  gint status    = 0;
  gint nchains_i = 0;

  fits_open_file (&fptr, filename, READONLY, &status);
  NCM_FITS_ERROR (status);

  fits_movnam_hdu (fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
  NCM_FITS_ERROR (status);

  fits_read_key (fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, first_id, NULL, &status);
  NCM_FITS_ERROR (status);

  fits_read_key (fptr, TINT, NCM_MSET_CATALOG_NCHAINS_LABEL, &nchains_i, NULL, &status);
  NCM_FITS_ERROR (status);
  g_assert_cmpint (nchains_i, >, 0);
  *nchains = (guint) nchains_i;

  fits_get_num_rows (fptr, nrows, &status);
  NCM_FITS_ERROR (status);

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

#else
  g_error ("ncm_mset_catalog_peek_info_from_file: cannot read file without cfitsio.");

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_peek_markovian_id_from_file:
 * @filename: catalog file name
 *
 * Reads #NcmMSetCatalog:markovian-id from the file header without loading the catalog,
 * in the file's own row numbering. A file without the key returns its first id (every
 * row Markovian).
 *
 * Returns: the id of the first row of the Markovian chain.
 */
gint
ncm_mset_catalog_peek_markovian_id_from_file (const gchar *filename)
{
#ifdef HAVE_CFITSIO
  fitsfile *fptr = NULL;
  gint status    = 0;
  gint first_id  = 0;
  gint markid    = 0;

  fits_open_file (&fptr, filename, READONLY, &status);
  NCM_FITS_ERROR (status);

  fits_movnam_hdu (fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
  NCM_FITS_ERROR (status);

  fits_read_key (fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, &first_id, NULL, &status);
  NCM_FITS_ERROR (status);

  fits_read_key (fptr, TINT, NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL, &markid, NULL, &status);

  if (status == KEY_NO_EXIST)
  {
    markid = first_id;
    status = 0;
  }
  else
  {
    NCM_FITS_ERROR (status);
  }

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

  return MAX (markid, first_id);

#else
  g_error ("ncm_mset_catalog_peek_markovian_id_from_file: numcosmo built without cfitsio.");

  return 0;

#endif /* HAVE_CFITSIO */
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

#ifdef HAVE_CFITSIO

static void
_ncm_fits_update_key_str (fitsfile *fptr, gchar *keyname, gchar *value, gchar *comment, gboolean overwrite)
{
  gint status = 0;

  if (overwrite)
  {
    fits_update_key_str (fptr, keyname, value, comment, &status);
    NCM_FITS_ERROR (status);
  }
  else
  {
    gchar key_text[FLEN_VALUE];
    gchar comment_text[FLEN_COMMENT];

    fits_read_key_str (fptr, keyname, key_text, comment_text, &status);

    if (status == 0)
    {
      g_assert_cmpstr (key_text, ==, value);

      if (comment != NULL)
        g_assert_cmpstr (comment_text, ==, comment);
    }
    else if (status == KEY_NO_EXIST)
    {
      status = 0;
      fits_update_key_str (fptr, keyname, value, comment, &status);
      NCM_FITS_ERROR (status);
    }
    else
    {
      NCM_FITS_ERROR (status);
    }
  }
}

static void
_ncm_fits_update_key_longstr (fitsfile *fptr, gchar *keyname, gchar *value, gchar *comment, gboolean overwrite)
{
  gint status = 0;

  if (overwrite)
  {
    fits_update_key_longstr (fptr, keyname, value, comment, &status);
    NCM_FITS_ERROR (status);
  }
  else
  {
    gchar *key_text;
    gchar comment_text[FLEN_COMMENT];

    fits_read_key_longstr (fptr, keyname, &key_text, comment_text, &status);

    if (status == 0)
    {
      g_assert_cmpstr (key_text, ==, value);

      if (comment != NULL)
        g_assert_cmpstr (comment_text, ==, comment);
    }
    else if (status == KEY_NO_EXIST)
    {
      status = 0;
      fits_update_key_longstr (fptr, keyname, value, comment, &status);
      NCM_FITS_ERROR (status);
    }
    else
    {
      NCM_FITS_ERROR (status);
    }

    fits_free_memory (key_text, &status);
    NCM_FITS_ERROR (status);
  }
}

/*
 * Records one of the informative strings: the header key @label holds the current value,
 * a HISTORY card with the date and @what keeps every value the catalog has seen.
 */
static void
_ncm_mset_catalog_write_record (fitsfile *fptr, const gchar *label, const gchar *what, const gchar *value)
{
  GDateTime *now = g_date_time_new_now_local ();
  gchar *stamp   = g_date_time_format (now, "%Y-%m-%dT%H:%M:%S");
  gchar *history = g_strdup_printf ("%s %s: %s", stamp, what, value);
  gchar *comment = g_strdup_printf ("%s; informative only.", what);
  gint status    = 0;

  _ncm_fits_update_key_longstr (fptr, (gchar *) label, (gchar *) value, comment, TRUE);
  fits_write_history (fptr, history, &status);
  NCM_FITS_ERROR (status);

  g_free (comment);
  g_free (history);
  g_free (stamp);
  g_date_time_unref (now);
}

/*
 * On opening a file: @field takes the file's value when loading from it or when unset;
 * otherwise a different value in memory is recorded, never rejected.
 */
static void
_ncm_mset_catalog_load_record (fitsfile *fptr, const gchar *label, const gchar *what, gchar **field,
                               gboolean load_from_cat, gboolean readonly)
{
  gchar *value = NULL;
  gint status  = 0;

  fits_read_key_longstr (fptr, (gchar *) label, &value, NULL, &status);

  if (status == KEY_NO_EXIST)
    status = 0;

  NCM_FITS_ERROR (status);

  if (load_from_cat || (*field == NULL))
  {
    if (value != NULL)
    {
      g_clear_pointer (field, g_free);
      *field = g_strdup (value);
    }
  }
  else if (((value == NULL) || (strcmp (value, *field) != 0)) && !readonly)
  {
    _ncm_mset_catalog_write_record (fptr, label, what, *field);
  }

  if (value != NULL)
    fits_free_memory (value, &status);
}

static void
_ncm_fits_update_key_int (fitsfile *fptr, gchar *keyname, gint value, gchar *comment, gboolean overwrite)
{
  gint status  = 0;
  glong lvalue = value;

  if (overwrite)
  {
    fits_update_key_lng (fptr, keyname, lvalue, comment, &status);
    NCM_FITS_ERROR (status);
  }
  else
  {
    glong key_value;
    gchar comment_text[FLEN_COMMENT];

    fits_read_key_lng (fptr, keyname, &key_value, comment_text, &status);

    if (status == 0)
    {
      g_assert_cmpint (key_value, ==, value);

      if (comment != NULL)
        g_assert_cmpstr (comment_text, ==, comment);
    }
    else if (status == KEY_NO_EXIST)
    {
      status = 0;
      fits_update_key_lng (fptr, keyname, lvalue, comment, &status);
      NCM_FITS_ERROR (status);
    }
    else
    {
      NCM_FITS_ERROR (status);
    }
  }
}

static void
_ncm_fits_update_key_ulong (fitsfile *fptr, gchar *keyname, gulong value, gchar *comment, gboolean overwrite)
{
  gint status = 0;

  if (overwrite)
  {
    fits_update_key_lng (fptr, keyname, value, comment, &status);
    NCM_FITS_ERROR (status);
  }
  else
  {
    glong key_value;
    gchar comment_text[FLEN_COMMENT];

    fits_read_key_lng (fptr, keyname, &key_value, comment_text, &status);

    if (status == 0)
    {
      g_assert_cmpint (key_value, ==, value);

      if (comment != NULL)
        g_assert_cmpstr (comment_text, ==, comment);
    }
    else if (status == KEY_NO_EXIST)
    {
      status = 0;
      fits_update_key_lng (fptr, keyname, value, comment, &status);
      NCM_FITS_ERROR (status);
    }
    else
    {
      NCM_FITS_ERROR (status);
    }
  }
}

/*
 * Writes the Markovian id in the file numbering (the in-memory id plus the burn-in
 * dropped at load) when it differs from what the file holds. Returns TRUE if written.
 */
static gboolean
_ncm_mset_catalog_sync_markovian_id (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const gint file_markid      = ncm_mset_catalog_get_markovian_id (mcat) + self->burnin;

  if ((self->fptr == NULL) || self->readonly || (file_markid == self->file_markovian_id))
    return FALSE;

  _ncm_fits_update_key_int (self->fptr, NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL, file_markid, "First row id of the Markovian chain.", TRUE);
  self->file_markovian_id = file_markid;

  return TRUE;
}

static void
_ncm_mset_catalog_sync_rng (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gchar key_text[FLEN_VALUE];
  gint status = 0;

  fits_read_key_str (self->fptr, NCM_MSET_CATALOG_RNG_ALGO_LABEL,
                     key_text, NULL, &status);

  if (status == 0)
  {
    glong seed  = 0;
    gchar *inis = NULL;

    fits_read_key_lng (self->fptr, NCM_MSET_CATALOG_RNG_SEED_LABEL,
                       &seed, NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_key_longstr (self->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, &inis, NULL, &status);
    NCM_FITS_ERROR (status);

    if (self->rng != NULL)
    {
      const gchar *cat_algo = ncm_rng_get_algo (self->rng);

      g_assert_cmpstr (cat_algo, ==, key_text);
      g_assert_cmpstr (inis, ==, self->rng_inis);
    }
    else
    {
      NcmRNG *rng = ncm_rng_new (key_text);

      ncm_rng_set_state (rng, inis);
      ncm_mset_catalog_set_rng (mcat, rng);
      ncm_rng_free (rng);
    }

    fits_free_memory (inis, &status);
    NCM_FITS_ERROR (status);
  }
  else if (status == KEY_NO_EXIST)
  {
    status = 0;

    if (self->rng != NULL)
    {
      gulong seed = ncm_rng_get_seed (self->rng);

      _ncm_fits_update_key_str (self->fptr, NCM_MSET_CATALOG_RNG_ALGO_LABEL, (gchar *) ncm_rng_get_algo (self->rng), "RNG Algorithm name.", FALSE);
      _ncm_fits_update_key_longstr (self->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, self->rng_inis, NULL, FALSE);
      _ncm_fits_update_key_ulong (self->fptr, NCM_MSET_CATALOG_RNG_SEED_LABEL, seed, "RNG Algorithm seed.", FALSE);
    }
  }
  else
  {
    NCM_FITS_ERROR (status);
  }
}

/*
 * Writes @mset (and, if non-%NULL, @functions) into the primary HDU (HDU0)
 * of @fptr, which must have just been created by fits_create_file() and
 * hold no data yet. This must happen before any table extension is
 * created, since resizing the primary HDU later on would require shifting
 * every following HDU.
 *
 * The envelope is a #NcmVarDict with a "model-set" entry (always) and a
 * "functions" entry (only when @functions is given), so the catalog file
 * is self-sufficient: it carries both the model-set and, when applicable,
 * the #NcmMSetFunc array that computed its extra columns, with no need for
 * a companion experiment or functions file to make sense of it.
 */
static void
_ncm_mset_catalog_write_hdu0 (fitsfile *fptr, NcmMSet *mset, NcmObjArray *functions)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmVarDict *vd    = ncm_var_dict_new ();
  GVariant *var;
  glong naxes[1];
  gint status = 0;

  ncm_var_dict_set_object (vd, NCM_MSET_CATALOG_HDU0_MSET_KEY, ser, G_OBJECT (mset));

  if (functions != NULL)
    ncm_var_dict_set_object_array (vd, NCM_MSET_CATALOG_HDU0_FUNCTIONS_KEY, ser, functions);

  var      = ncm_serialize_var_dict_to_variant (ser, vd);
  naxes[0] = (glong) g_variant_get_size (var);

  fits_create_img (fptr, BYTE_IMG, 1, naxes, &status);
  NCM_FITS_ERROR (status);

  fits_write_img (fptr, TBYTE, 1, naxes[0], (gpointer) g_variant_get_data (var), &status);
  NCM_FITS_ERROR (status);

  fits_write_key_str (fptr, NCM_MSET_CATALOG_MSET_FORMAT_LABEL, NCM_MSET_CATALOG_MSET_FORMAT_VARDICT,
                      "Format of the data stored in this HDU.", &status);
  NCM_FITS_ERROR (status);

  g_variant_unref (var);
  ncm_var_dict_unref (vd);
  ncm_serialize_free (ser);
}

/*
 * Reads back the NcmMSet (and, if present, the functions array) embedded
 * by _ncm_mset_catalog_write_hdu0(), if any. Leaves both out parameters
 * %NULL for legacy files whose primary HDU is the default empty one (no
 * data written to HDU0). Understands both the current vardict format and
 * the older bare-mset-object format (pre-existing catalogs that only ever
 * embedded the model-set).
 */
static void
_ncm_mset_catalog_read_hdu0 (fitsfile *fptr, NcmMSet **mset, NcmObjArray **functions)
{
  gint status    = 0;
  gint bitpix    = 0;
  gint naxis     = 0;
  gint hdutype   = 0;
  glong naxes[1] = { 0 };

  *mset      = NULL;
  *functions = NULL;

  fits_movabs_hdu (fptr, 1, &hdutype, &status);
  NCM_FITS_ERROR (status);

  fits_get_img_param (fptr, 1, &bitpix, &naxis, naxes, &status);
  NCM_FITS_ERROR (status);

  if ((naxis == 1) && (naxes[0] > 0))
  {
    const glong len = naxes[0];
    guint8 *data    = g_malloc (len);
    gint any_null   = 0;
    gchar format[FLEN_VALUE];
    NcmSerialize *ser;

    fits_read_img (fptr, TBYTE, 1, len, NULL, data, &any_null, &status);
    NCM_FITS_ERROR (status);

    fits_read_key (fptr, TSTRING, NCM_MSET_CATALOG_MSET_FORMAT_LABEL, format, NULL, &status);

    if (status == KEY_NO_EXIST)
    {
      status = 0;
      g_strlcpy (format, NCM_MSET_CATALOG_MSET_FORMAT_OBJECT, sizeof (format));
    }

    NCM_FITS_ERROR (status);

    ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

    if (strcmp (format, NCM_MSET_CATALOG_MSET_FORMAT_VARDICT) == 0)
    {
      GVariant *var = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE),
                                               data, len, TRUE, g_free, data);
      NcmVarDict *vd    = ncm_serialize_var_dict_from_variant (ser, var);
      GObject *mset_obj = NULL;

      if (!ncm_var_dict_get_object (vd, NCM_MSET_CATALOG_HDU0_MSET_KEY, ser, &mset_obj))
        g_error ("_ncm_mset_catalog_read_hdu0: HDU0 vardict has no `%s' entry.",
                 NCM_MSET_CATALOG_HDU0_MSET_KEY);

      g_assert (NCM_IS_MSET (mset_obj));
      *mset = NCM_MSET (mset_obj);
      ncm_var_dict_get_object_array (vd, NCM_MSET_CATALOG_HDU0_FUNCTIONS_KEY, ser, functions);

      ncm_var_dict_unref (vd);
      g_variant_unref (var);
    }
    else if (strcmp (format, NCM_MSET_CATALOG_MSET_FORMAT_OBJECT) == 0)
    {
      GVariant *var = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE),
                                               data, len, TRUE, g_free, data);
      GObject *mset_obj = ncm_serialize_from_variant (ser, var);

      g_assert (NCM_IS_MSET (mset_obj));
      *mset = NCM_MSET (mset_obj);
      g_variant_unref (var);
    }
    else
    {
      g_error ("_ncm_mset_catalog_read_hdu0: unknown HDU0 format `%s'.", format);
    }

    ncm_serialize_free (ser);
  }
}

/*
 * Resolves the NcmMSet (and functions array, if any) for a catalog being
 * loaded from an existing file: the ones embedded in HDU0 if present,
 * otherwise the legacy `.mset' GKeyFile sidecar for the mset (read-only
 * support; new catalogs no longer write that file) and no functions.
 */
static NcmMSet *
_ncm_mset_catalog_load_mset (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  NcmMSet *mset               = NULL;
  NcmObjArray *functions      = NULL;
  fitsfile *fptr              = NULL;
  gint status                 = 0;

  fits_open_file (&fptr, self->file, READONLY, &status);
  NCM_FITS_ERROR (status);

  _ncm_mset_catalog_read_hdu0 (fptr, &mset, &functions);

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

  if (mset == NULL)
  {
    if ((self->mset_file == NULL) || !g_file_test (self->mset_file, G_FILE_TEST_EXISTS))
      g_error ("_ncm_mset_catalog_load_mset: catalog file `%s' has no embedded mset "
               "and no legacy mset file `%s' was found.",
               self->file, self->mset_file != NULL ? self->mset_file : "(null)");

    {
      NcmSerialize *ser = ncm_serialize_global ();

      mset = ncm_mset_load (self->mset_file, ser, NULL);
      ncm_serialize_free (ser);
    }
  }

  ncm_obj_array_clear (&self->functions);
  self->functions = functions;

  return mset;
}

static void
_ncm_mset_catalog_open_create_file (NcmMSetCatalog *mcat, gboolean load_from_cat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint fparam_len            = ncm_mset_fparam_len (self->mset);
  gchar key_text[FLEN_VALUE];
  gint status = 0;
  guint i;

  g_assert (self->file != NULL);
  g_assert (self->fptr == NULL);

  if (g_file_test (self->file, G_FILE_TEST_EXISTS))
  {
    gboolean weighted       = FALSE;
    gint nchains            = 0;
    gint nadd_vals          = 0;
    gboolean remap          = FALSE;
    GPtrArray *remap_remove = g_ptr_array_new ();
    glong nrows;

    g_ptr_array_set_free_func (remap_remove, g_free);

    if (self->readonly)
    {
      fits_open_file (&self->fptr, self->file, READONLY, &status);
      NCM_FITS_ERROR (status);
    }
    else
    {
      fits_open_file (&self->fptr, self->file, READWRITE, &status);
      NCM_FITS_ERROR (status);
    }

    fits_movnam_hdu (self->fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
    NCM_FITS_ERROR (status);

    fits_read_key (self->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL,
                   &self->file_first_id, NULL, &status);
    NCM_FITS_ERROR (status);

    {
      const gint m2lnp_var = self->m2lnp_var;

      fits_read_key (self->fptr, TINT, NCM_MSET_CATALOG_M2LNP_ID_LABEL,
                     &self->m2lnp_var, NULL, &status);

      if (status == KEY_NO_EXIST)
      {
        g_warning ("_ncm_mset_catalog_open_create_file: catalog does not contain `%s' key, using original value `%d'.",
                   NCM_MSET_CATALOG_M2LNP_ID_LABEL,
                   m2lnp_var);
        self->m2lnp_var = m2lnp_var;
      }
      else
      {
        NCM_FITS_ERROR (status);
      }

      status = 0;
    }

    {
      gint file_markid = self->file_first_id;

      fits_read_key (self->fptr, TINT, NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL, &file_markid, NULL, &status);

      if (status == KEY_NO_EXIST)
        file_markid = self->file_first_id;  /* written before the key existed: every row Markovian */
      else
        NCM_FITS_ERROR (status);

      status = 0;
      /* The file numbering; converted to the in-memory numbering once the burn-in is known below. */
      self->file_markovian_id = file_markid;
    }

    fits_read_key (self->fptr, TSTRING, NCM_MSET_CATALOG_RTYPE_LABEL,
                   key_text, NULL, &status);
    NCM_FITS_ERROR (status);

    if (load_from_cat)
      ncm_mset_catalog_set_run_type (mcat, key_text);
    else if (strcmp (self->rtype_str, key_text) != 0)
      g_error ("_ncm_mset_catalog_open_create_file: incompatible run type strings from catalog and file, catalog: `%s' file: `%s'.",
               self->rtype_str, key_text);

    /* Records, not constraints: a catalog started by one sampler may be continued by another. */
    _ncm_mset_catalog_load_record (self->fptr, NCM_MSET_CATALOG_INIT_SAMPLER_LABEL, "initial points sampler",
                                   &self->init_sampler_str, load_from_cat, self->readonly);
    _ncm_mset_catalog_load_record (self->fptr, NCM_MSET_CATALOG_SAMPLER_LABEL, "sampler",
                                   &self->sampler_str, load_from_cat, self->readonly);
    _ncm_mset_catalog_load_record (self->fptr, NCM_MSET_CATALOG_SAMPLER_OPTS_LABEL, "sampler options",
                                   &self->sampler_opts_str, load_from_cat, self->readonly);

    fits_read_key (self->fptr, TINT, NCM_MSET_CATALOG_NCHAINS_LABEL,
                   &nchains, NULL, &status);
    NCM_FITS_ERROR (status);
    g_assert_cmpint (nchains, >, 0);

    if (load_from_cat)
      self->nchains = nchains;
    else if (nchains != (gint) self->nchains)
      g_error ("_ncm_mset_catalog_open_create_file: catalog has %d chains and file contains %d.", self->nchains, nchains);

    fits_read_key (self->fptr, TINT, NCM_MSET_CATALOG_NADDVAL_LABEL,
                   &nadd_vals, NULL, &status);
    NCM_FITS_ERROR (status);

    if (load_from_cat)
      self->nadd_vals = nadd_vals;
    else if (nadd_vals != (gint) self->nadd_vals)
      g_error ("_ncm_mset_catalog_open_create_file: catalog has %d additional values and file contains %d.", self->nadd_vals, nadd_vals);

    fits_read_key (self->fptr, TLOGICAL, NCM_MSET_CATALOG_WEIGHTED_LABEL,
                   &weighted, NULL, &status);
    NCM_FITS_ERROR (status);

    if (load_from_cat)
      self->weighted = weighted ? TRUE : FALSE;
    else if ((weighted && !self->weighted) || (!weighted && self->weighted))
      g_error ("_ncm_mset_catalog_open_create_file: catalog %s weighted and file %s.",
               self->weighted ? "is" : "is not",
               weighted ? "is" : "is not");

    fits_get_num_rows (self->fptr, &nrows, &status);
    NCM_FITS_ERROR (status);

    if (nrows < self->burnin)
      g_error ("_ncm_mset_catalog_open_create_file: burnin of %ld row(s) exceeds "
               "catalog `%s' (%ld row(s), %u chain(s), i.e. %ld iteration(s)).",
               self->burnin, self->file, nrows, self->nchains, nrows / self->nchains);
    else
      nrows -= self->burnin;

    if (self->file_first_id != self->first_id)
    {
      if (nrows == 0)
      {
        if (self->file_first_id != 0)
          g_warning ("_ncm_mset_catalog_open_create_file: Empty data file with "NCM_MSET_CATALOG_FIRST_ID_LABEL " different from first_id: %d != %d. Setting to first_id.\n",
                     self->file_first_id, self->first_id);

        self->file_first_id = self->first_id;
      }
      else if (ncm_mset_catalog_is_empty (mcat))
      {
        if (self->first_id != 0)
          g_warning ("_ncm_mset_catalog_open_create_file: Empty memory catalog with first_id different from "NCM_MSET_CATALOG_FIRST_ID_LABEL ": %d != %d. Setting to "NCM_MSET_CATALOG_FIRST_ID_LABEL ".\n",
                     self->first_id, self->file_first_id);

        self->first_id = self->file_first_id;
        self->cur_id   = self->file_first_id - 1;
      }
    }

    self->file_cur_id = self->file_first_id + nrows - 1;

    /* Rows before the burn-in are dropped and the ids renumbered from file_first_id, so the
     * stored Markovian id shifts by the same amount; it is never below the first id. */
    self->markovian_id = MAX (self->markovian_id, (gint) (self->file_markovian_id - self->burnin));

    if (load_from_cat)
    {
      gchar colname[FLEN_VALUE];
      gint cindex = 0;
      guint total = fparam_len + self->nadd_vals + (self->weighted ? 1 : 0);

      g_array_set_size (self->porder, total);

      i = 0;

      while (fits_get_colname (self->fptr, CASESEN, "*", colname, &cindex, &status) == COL_NOT_UNIQUE)
      {
        gchar *d_colname = g_strdup (colname);

        g_assert_cmpint (i + 1, ==, cindex);

        status = 0;

        if (i < self->nadd_vals)
        {
          g_ptr_array_add (self->add_vals_names, d_colname);
          g_array_index (self->porder, gint, i) = cindex;

          {
            gchar symbol_s[FLEN_VALUE];
            gchar *symbol;
            gchar *asymbi = g_strdup_printf ("%s%d", NCM_MSET_CATALOG_ASYMB_LABEL, i + 1);

            fits_read_key (self->fptr, TSTRING, asymbi, &symbol_s, NULL, &status);

            if (status == KEY_NO_EXIST)
            {
              symbol = g_strdup ("no-symbol");
              status = 0;
            }
            else
            {
              symbol = g_strdup (symbol_s);
              NCM_FITS_ERROR (status);
            }

            g_ptr_array_add (self->add_vals_symbs, symbol);

            g_free (asymbi);
          }
        }
        else
        {
          NcmMSetPIndex *pi = ncm_mset_param_get_by_full_name (self->mset, d_colname, NULL);

          if (pi == NULL)
          {
            g_error ("_ncm_mset_catalog_open_create_file: cannot find parameter `%s' in mset file.", d_colname);
          }
          else
          {
            NcmParamType ftype = ncm_mset_param_get_ftype (self->mset, pi->mid, pi->pid, NULL);

            if (ftype != NCM_PARAM_TYPE_FREE)
            {
              g_warning ("_ncm_mset_catalog_open_create_file: parameter `%s' found but not free on the catalog, setting it to NCM_PARAM_TYPE_FREE.",
                         d_colname);
              ncm_mset_param_set_ftype (self->mset, pi->mid, pi->pid, NCM_PARAM_TYPE_FREE);
              remap = TRUE;
            }

            ncm_mset_pindex_free (pi);
          }

          g_free (d_colname);
        }

        status = COL_NOT_UNIQUE;
        i++;
      }

      status = 0;
    }
    else
    {
      for (i = 0; i < self->nadd_vals; i++)
      {
        const gchar *cname   = g_ptr_array_index (self->add_vals_names, i);
        const gchar *csymbol = g_ptr_array_index (self->add_vals_symbs, i);
        gchar *asymbi        = g_strdup_printf ("%s%d", NCM_MSET_CATALOG_ASYMB_LABEL, i + 1);
        gchar symbol_s[FLEN_VALUE];
        gint cindex = 0;

        if (fits_get_colnum (self->fptr, CASESEN, (gchar *) cname, &cindex, &status))
          g_error ("_ncm_mset_catalog_open_create_file: Additional column %s not found, invalid fits file.", cname);

        fits_read_key (self->fptr, TSTRING, asymbi, &symbol_s, NULL, &status);

        if (status == KEY_NO_EXIST)
          g_error ("_ncm_mset_catalog_open_create_file: symbol %s not found", asymbi);
        else
          g_assert_cmpstr (symbol_s, ==, csymbol);

        if (cindex != (gint) (i + 1))
          g_error ("_ncm_mset_catalog_open_create_file: Additional column %s is not the %d-th column [%d], invalid fits file.",
                   cname, i + 1, cindex);

        g_array_index (self->porder, gint, i) = cindex;
        g_free (asymbi);
      }
    }

    if (remap)
    {
      guint total;

      ncm_mset_prepare_fparam_map (self->mset);

      fparam_len = ncm_mset_fparam_len (self->mset);
      total      = fparam_len + self->nadd_vals + (self->weighted ? 1 : 0);

      g_array_set_size (self->porder, total);
      remap = FALSE;
    }

    for (i = 0; i < fparam_len; i++)
    {
      const gchar *fparam_fullname = ncm_mset_fparam_full_name (self->mset, i);

      if (fits_get_colnum (self->fptr, CASESEN, (gchar *) fparam_fullname, &g_array_index (self->porder, gint, i + self->nadd_vals), &status))
      { /* I don't like this too ^^^^^^^^^  */
        g_warning ("_ncm_mset_catalog_open_create_file: Parameter `%s' set free in mset but not found on the fits file, setting it to NCM_PARAM_TYPE_FIXED.", fparam_fullname);
        g_ptr_array_add (remap_remove, (gpointer) g_strdup (fparam_fullname));
        remap  = TRUE;
        status = 0;
      }
    }

    if (remap)
    {
      guint total;

      for (i = 0; i < remap_remove->len; i++)
      {
        NcmMSetPIndex *pi = ncm_mset_param_get_by_full_name (self->mset, g_ptr_array_index (remap_remove, i), NULL);

        if (pi == NULL)
        {
          g_error ("_ncm_mset_catalog_open_create_file: unknown error should never happen! Cannot find parameter `%s' in mset file.",
                   (gchar *) g_ptr_array_index (remap_remove, i));
        }
        else
        {
          ncm_mset_param_set_ftype (self->mset, pi->mid, pi->pid, NCM_PARAM_TYPE_FIXED);
          ncm_mset_pindex_free (pi);
        }
      }

      ncm_mset_prepare_fparam_map (self->mset);

      fparam_len = ncm_mset_fparam_len (self->mset);
      total      = fparam_len + self->nadd_vals + (self->weighted ? 1 : 0);

      g_array_set_size (self->porder, total);

      for (i = 0; i < fparam_len; i++)
      {
        const gchar *fparam_fullname = ncm_mset_fparam_full_name (self->mset, i);

        if (fits_get_colnum (self->fptr, CASESEN, (gchar *) fparam_fullname, &g_array_index (self->porder, gint, i + self->nadd_vals), &status)) /* I don't like this too ^^^^^^^^^  */
          g_error ("_ncm_mset_catalog_open_create_file: Parameter `%s' set free in mset but not found on the fits file, this should never happen!", fparam_fullname);
      }

      remap = FALSE;
    }

    g_ptr_array_unref (remap_remove);
  }
  else
  {
    GPtrArray *ttype_array = g_ptr_array_sized_new (10);
    GPtrArray *tform_array = g_ptr_array_sized_new (10);

    fits_create_file (&self->fptr, self->file, &status);
    NCM_FITS_ERROR (status);

    _ncm_mset_catalog_write_hdu0 (self->fptr, self->mset, self->functions);

    for (i = 0; i < self->nadd_vals; i++)
    {
      const gchar *cname = g_ptr_array_index (self->add_vals_names, i);

      g_ptr_array_add (ttype_array, (gchar *) cname);
      g_ptr_array_add (tform_array, "1D");
      g_array_index (self->porder, gint, i) = tform_array->len;
    }

    for (i = 0; i < fparam_len; i++)
    {
      const gchar *fparam_fullname = ncm_mset_fparam_full_name (self->mset, i);

      g_ptr_array_add (ttype_array, (gchar *) fparam_fullname);
      /* I don't like this too ^^^^^^^^^  */
      g_ptr_array_add (tform_array, "1D");
      g_array_index (self->porder, gint, i + self->nadd_vals) = tform_array->len;
    }

    /* append a new empty binary table onto the FITS file */
    fits_create_tbl (self->fptr, BINARY_TBL, 0, fparam_len + self->nadd_vals, (gchar **) ttype_array->pdata, (gchar **) tform_array->pdata,
                     NULL, NCM_MSET_CATALOG_EXTNAME, &status);
    NCM_FITS_ERROR (status);

    fits_update_key (self->fptr, TSTRING, NCM_MSET_CATALOG_RTYPE_LABEL, self->rtype_str, "Run type string.", &status);
    NCM_FITS_ERROR (status);

    if (self->init_sampler_str != NULL)
      _ncm_mset_catalog_write_record (self->fptr, NCM_MSET_CATALOG_INIT_SAMPLER_LABEL, "initial points sampler", self->init_sampler_str);

    if (self->sampler_str != NULL)
      _ncm_mset_catalog_write_record (self->fptr, NCM_MSET_CATALOG_SAMPLER_LABEL, "sampler", self->sampler_str);

    if (self->sampler_opts_str != NULL)
      _ncm_mset_catalog_write_record (self->fptr, NCM_MSET_CATALOG_SAMPLER_OPTS_LABEL, "sampler options", self->sampler_opts_str);

    fits_update_key (self->fptr, TINT, NCM_MSET_CATALOG_NCHAINS_LABEL, &self->nchains, "Number of chains.", &status);
    NCM_FITS_ERROR (status);

    fits_update_key (self->fptr, TINT, NCM_MSET_CATALOG_NADDVAL_LABEL, &self->nadd_vals, "Number of additional values.", &status);
    NCM_FITS_ERROR (status);

    fits_update_key (self->fptr, TLOGICAL, NCM_MSET_CATALOG_WEIGHTED_LABEL, &self->weighted, "Whether the catalog is weighted.", &status);
    NCM_FITS_ERROR (status);

    for (i = 0; i < self->nadd_vals; i++)
    {
      const gchar *aname = g_ptr_array_index (self->add_vals_names, i);
      const gchar *asymb = g_ptr_array_index (self->add_vals_symbs, i);

      gchar *asymbi     = g_strdup_printf ("%s%d", NCM_MSET_CATALOG_ASYMB_LABEL, i + 1);
      gchar *asymb_desc = g_strdup_printf ("Symbol for additional value %s[%d]",
                                           aname,
                                           i + 1);

      fits_update_key (self->fptr, TSTRING, asymbi, (gchar *) asymb, asymb_desc, &status);
      NCM_FITS_ERROR (status);

      g_free (asymbi);
      g_free (asymb_desc);
    }

    for (i = 0; i < fparam_len; i++)
    {
      gchar *fsymbi     = g_strdup_printf ("%s%d", NCM_MSET_CATALOG_FSYMB_LABEL, i + 1);
      gchar *fsymb_desc = g_strdup_printf ("Symbol for parameter %s[%d]",
                                           ncm_mset_fparam_name (self->mset, i),
                                           i + 1);
      const gchar *fsymb = ncm_mset_fparam_symbol (self->mset, i);

      fits_update_key (self->fptr, TSTRING, fsymbi, (gchar *) fsymb, fsymb_desc, &status);
      NCM_FITS_ERROR (status);

      g_free (fsymbi);
      g_free (fsymb_desc);
    }

    self->file_first_id = self->first_id;
    self->file_cur_id   = self->first_id - 1;

    g_ptr_array_unref (ttype_array);
    g_ptr_array_unref (tform_array);
  }

  _ncm_mset_catalog_sync_rng (mcat);

  _ncm_fits_update_key_int (self->fptr, NCM_MSET_CATALOG_FIRST_ID_LABEL, self->file_first_id, "Id of the first element.", !self->readonly);
  _ncm_fits_update_key_int (self->fptr, NCM_MSET_CATALOG_M2LNP_ID_LABEL, self->m2lnp_var,     "Id of the m2lnp variable.", !self->readonly);
  _ncm_mset_catalog_sync_markovian_id (mcat);

  if (!self->readonly)
  {
    fits_flush_file (self->fptr, &status);
    NCM_FITS_ERROR (status);
  }
}

#endif /* HAVE_CFITSIO */

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
 * @filename: a filename
 *
 * Sets the data filename to be used to sync/save data.
 *
 */
void
ncm_mset_catalog_set_file (NcmMSetCatalog *mcat, const gchar *filename)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

#ifdef HAVE_CFITSIO

  if (!self->constructed)
  {
    if (self->file != NULL)
      g_error ("ncm_mset_catalog_set_file: file already set before construction.");

    self->file = g_strdup (filename);
  }

  if ((self->file != NULL) && (filename != NULL) && (strcmp (self->file, filename) == 0))
    return;

  _ncm_mset_catalog_close_file (mcat);

  g_clear_pointer (&self->file, g_free);
  g_clear_pointer (&self->mset_file, g_free);

  if (filename == NULL)
    return;

  self->file = g_strdup (filename);
  {
    gchar *base_name = ncm_util_basename_fits (self->file);

    self->mset_file = g_strdup_printf ("%s.mset", base_name);
    g_free (base_name);
  }

  if (self->mset != NULL)
  {
    _ncm_mset_catalog_open_create_file (mcat, FALSE);
    ncm_mset_catalog_sync (mcat, TRUE);
  }

#else
  g_error ("ncm_mset_catalog_set_file: cannot set file without cfitsio.");
#endif /* HAVE_CFITSIO */

  self->first_flush = TRUE;
}

/**
 * ncm_mset_catalog_set_sync_mode:
 * @mcat: a #NcmMSetCatalog
 * @smode: sync mode
 *
 * Sets the sync mode to @smode.
 *
 */
void
ncm_mset_catalog_set_sync_mode (NcmMSetCatalog *mcat, NcmMSetCatalogSync smode)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->smode = smode;
}

/**
 * ncm_mset_catalog_set_sync_interval:
 * @mcat: a #NcmMSetCatalog
 * @interval: minimum time between syncs, in seconds
 *
 * Sets the minimum time between the syncs of ncm_mset_catalog_timed_sync().
 *
 */
void
ncm_mset_catalog_set_sync_interval (NcmMSetCatalog *mcat, gdouble interval)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->sync_interval = interval;
}

/**
 * ncm_mset_catalog_set_first_id:
 * @mcat: a #NcmMSetCatalog
 * @first_id: id of the first row
 *
 * Sets the id of the first row, e.g. the realization a #NcmFitMC catalog starts at. On a
 * non-empty catalog it aborts.
 *
 */
void
ncm_mset_catalog_set_first_id (NcmMSetCatalog *mcat, gint first_id)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (first_id == self->first_id)
    return;

  g_assert_cmpint (self->file_first_id, ==, self->first_id);
  g_assert_cmpint (self->file_cur_id,   ==, self->cur_id);

  if (!ncm_mset_catalog_is_empty (mcat))
    g_error ("ncm_mset_catalog_set_first_id: cannot modify first_id to %d in a non-empty catalog, catalog first id: %d, catalog current id: %d.",
             first_id, self->first_id, self->cur_id);

  self->first_id = first_id;
  self->cur_id   = first_id - 1;

  self->file_first_id = first_id;
  self->file_cur_id   = first_id - 1;
#ifdef HAVE_CFITSIO

  if (self->fptr != NULL)
  {
    _ncm_fits_update_key_int (self->fptr, NCM_MSET_CATALOG_FIRST_ID_LABEL, self->file_first_id, "Id of the first element.", !self->readonly);
    ncm_mset_catalog_sync (mcat, TRUE);
  }

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_set_run_type:
 * @mcat: a #NcmMSetCatalog
 * @rtype_str: description of the run
 *
 * Sets #NcmMSetCatalog:run-type-string, recorded in the file.
 *
 */
void
ncm_mset_catalog_set_run_type (NcmMSetCatalog *mcat, const gchar *rtype_str)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  g_assert (rtype_str != NULL);

  if (self->rtype_str != NULL)
  {
    if (strcmp (self->rtype_str, rtype_str) != 0)
    {
      if (!ncm_mset_catalog_is_empty (mcat))
        g_error ("ncm_mset_catalog_set_run_type: cannot change run type string in a non-empty catalog, actual: `%s' new: `%s'.",
                 self->rtype_str, rtype_str);
      else
        g_clear_pointer (&self->rtype_str, g_free);
    }
    else
    {
      return;
    }
  }

  self->rtype_str = g_strdup (rtype_str);
#ifdef HAVE_CFITSIO

  if (self->fptr != NULL)
    _ncm_fits_update_key_str (self->fptr, NCM_MSET_CATALOG_RTYPE_LABEL, self->rtype_str, NULL, !self->readonly);

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_set_rng:
 * @mcat: a #NcmMSetCatalog
 * @rng: a #NcmRNG
 *
 * Makes @rng the random number generator of @mcat. On a non-empty catalog it aborts:
 * the catalog continues from the generator state saved with its rows
 * (ncm_mset_catalog_peek_rng()), and a new generator would repeat rows already in it.
 */
void
ncm_mset_catalog_set_rng (NcmMSetCatalog *mcat, NcmRNG *rng)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (!ncm_mset_catalog_is_empty (mcat))
    g_error ("ncm_mset_catalog_set_rng: refusing to set RNG in a non-empty catalog (first id: %d, current id: %d): "
             "this would discard the saved RNG state and repeat rows already computed. "
             "Do not pass an explicit RNG/seed when resuming; the catalog's saved state is used.",
             self->first_id, self->cur_id);

  ncm_rng_clear (&self->rng);
  self->rng = ncm_rng_ref (rng);

  g_clear_pointer (&self->rng_inis, g_free);
  g_clear_pointer (&self->rng_stat, g_free);

  self->rng_inis = ncm_rng_get_state (rng);
  self->rng_stat = g_strdup (self->rng_inis);
#ifdef HAVE_CFITSIO

  if (self->fptr != NULL)
  {
    _ncm_mset_catalog_sync_rng (mcat);

    if (!self->readonly)
    {
      gint status = 0;

      fits_flush_file (self->fptr, &status);
      NCM_FITS_ERROR (status);
    }
  }

#endif /* HAVE_CFITSIO */
}

#ifdef HAVE_CFITSIO

static void
_ncm_mset_catalog_flush_file (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gint status                 = 0;

  /*gint64 nrows = self->file_cur_id - self->file_first_id + 1;*/

  /*printf ("# Flush: Updating to %ld nrows AXIS2!\n", nrows);*/
  /*fits_update_key (self->fptr, TLONGLONG, NCM_MSET_CATALOG_NROWS_LABEL, &nrows, NULL, &status);*/
  /*NCM_FITS_ERROR (status);*/

  if (G_UNLIKELY (self->first_flush))
  {
    fits_flush_file (self->fptr, &status);
    NCM_FITS_ERROR (status);
    self->first_flush = FALSE;
  }
  else
  {
    fits_flush_buffer (self->fptr, 0, &status);
    NCM_FITS_ERROR (status);
  }
}

static void
_ncm_mset_catalog_close_file (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gint status                 = 0;

  if (self->fptr != NULL)
  {
    ncm_mset_catalog_sync (mcat, FALSE);
    fits_close_file (self->fptr, &status);
    NCM_FITS_ERROR (status);
    self->fptr = NULL;

    g_clear_pointer (&self->file, g_free);
  }
}

static void
_ncm_mset_catalog_write_row (NcmMSetCatalog *mcat, NcmVector *row, guint row_index)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gint status                 = 0;
  guint i;

  /*printf ("Writing %u\n", row_index);*/

  for (i = 0; i < ncm_vector_len (row); i++)
  {
    fits_write_col_dbl (self->fptr, g_array_index (self->porder, gint, i), row_index + self->burnin,
                        1, 1, ncm_vector_ptr (row, i), &status);
    /*printf ("writing[%d]... %u %u == % 20.15g\n", g_array_index (self->porder, gint, i), i, row_index, ncm_vector_get (row, i));*/
    NCM_FITS_ERROR (status);
  }
}

static void
_ncm_mset_catalog_read_row (NcmMSetCatalog *mcat, NcmVector *row, guint row_index)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint i;
  gint status         = 0;
  const gdouble dnull = 0.0;

  /*printf ("Reading size: row %ld\n", row_index);*/

  for (i = 0; i < ncm_vector_len (row); i++)
  {
    fits_read_col_dbl (self->fptr, g_array_index (self->porder, gint, i), row_index + self->burnin,
                       1, 1, dnull, ncm_vector_ptr (row, i), NULL, &status);
    /*printf ("reading[%d]... %u %u == % 20.15g <%d>\n", g_array_index (self->porder, gint, i), i, row_index, ncm_vector_get (row, i), status);*/
    NCM_FITS_ERROR (status);
  }
}

#endif /* HAVE_CFITSIO */

static void _ncm_mset_catalog_post_update (NcmMSetCatalog *mcat, NcmVector *x);

/**
 * ncm_mset_catalog_sync:
 * @mcat: a #NcmMSetCatalog
 * @check: whether to check that the file and the memory agree
 *
 * Writes the rows not yet in the file, and the generator state, to the file of @mcat;
 * without a file it does nothing.
 *
 */
void
ncm_mset_catalog_sync (NcmMSetCatalog *mcat, gboolean check)
{
#ifdef HAVE_CFITSIO
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  gint status                 = 0;
  guint i;
  gboolean need_flush = FALSE;

  /*printf ("# Sync: start!\n");*/

  if (self->file == NULL)
    return;

  g_assert (self->fptr != NULL);

  /*printf ("# Sync: check %d\n", check);*/
  if (check)
  {
    gchar fptr_filename[FLEN_FILENAME];

    fits_file_name (self->fptr, fptr_filename, &status);
    NCM_FITS_ERROR (status);

    g_assert_cmpstr (fptr_filename, ==, self->file);

    if ((self->file_cur_id < self->first_id - 1) || (self->cur_id < self->file_first_id - 1))
      g_error ("ncm_mset_catalog_sync: file data & catalog mismatch, they do not intersect each other: file data [%d, %d] catalog [%d, %d]",
               self->file_first_id, self->file_cur_id,
               self->first_id, self->cur_id);
  }

  /*printf ("# Sync: file_first_id != self->first_id %d != %d\n", self->file_first_id, self->first_id);*/
  if (self->file_first_id != self->first_id)
  {
    if (self->file_first_id > self->first_id)
    {
      guint rows_to_add = self->file_first_id - self->first_id;

      fits_insert_rows (self->fptr, 0, rows_to_add, &status);
      NCM_FITS_ERROR (status);

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_stats_vec_peek_row (self->pstats, i);

        _ncm_mset_catalog_write_row (mcat, row, i + 1);
      }

      self->file_first_id = self->first_id;

      if (self->rng != NULL)
      {
        fits_update_key_longstr (self->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, self->rng_inis, NULL, &status);
        NCM_FITS_ERROR (status);
      }

      fits_update_key (self->fptr, TINT, NCM_MSET_CATALOG_FIRST_ID_LABEL, &self->file_first_id, "Id of the first element.", &status);
      NCM_FITS_ERROR (status);

      fits_update_key (self->fptr, TINT, NCM_MSET_CATALOG_M2LNP_ID_LABEL, &self->m2lnp_var,     "Id of the m2lnp variable.", &status);
      NCM_FITS_ERROR (status);

      need_flush = TRUE;
    }
    else if (self->file_first_id < self->first_id)
    {
      guint rows_to_add = self->first_id - self->file_first_id;
      GPtrArray *rows   = g_ptr_array_new ();
      gchar *inis       = NULL;

      g_ptr_array_set_size (rows, rows_to_add);

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_vector_dup (ncm_stats_vec_peek_x (self->pstats));

        _ncm_mset_catalog_read_row (mcat, row, i + 1);
        g_ptr_array_index (rows, i) = row;
      }

      ncm_stats_vec_prepend_data (self->pstats, rows, FALSE);

      if (self->nchains > 1)
      {
        for (i = 0; i < rows->len; i++)
        {
          NcmVector *x        = g_ptr_array_index (rows, i);
          guint chain_id      = (self->file_first_id + i) % self->nchains;
          NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, chain_id);

          ncm_stats_vec_prepend (pstats, x, FALSE);
        }
      }

      g_ptr_array_unref (rows);
      self->first_id = self->file_first_id;

      if (self->rng != NULL)
      {
        fits_read_key_longstr (self->fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, &inis, NULL, &status);
        NCM_FITS_ERROR (status);

        g_clear_pointer (&self->rng_inis, g_free);
        self->rng_inis = g_strdup (inis);

        fits_free_memory (inis, &status);
        NCM_FITS_ERROR (status);
      }
    }

    g_assert_cmpint (self->file_first_id, ==, self->first_id);
  }

  /*printf ("# Sync: self->file_cur_id != self->cur_id %d != %d\n", self->file_cur_id, self->cur_id);*/
  if (self->file_cur_id != self->cur_id)
  {
    if (self->file_cur_id < self->cur_id)
    {
      guint rows_to_add = self->cur_id - self->file_cur_id;
      guint offset      = self->file_cur_id + 1 - self->file_first_id;

      /*printf ("Adding %u rows after %u\n", rows_to_add, offset);*/
      fits_insert_rows (self->fptr, offset, rows_to_add, &status);
      NCM_FITS_ERROR (status);

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_stats_vec_peek_row (self->pstats, offset + i);

        _ncm_mset_catalog_write_row (mcat, row, offset + i + 1);
      }

      self->file_cur_id = self->cur_id;

      if (self->rng != NULL)
      {
        g_clear_pointer (&self->rng_stat, g_free);
        self->rng_stat = ncm_rng_get_state (self->rng);

        fits_update_key_longstr (self->fptr, NCM_MSET_CATALOG_RNG_STAT_LABEL, self->rng_stat, NULL, &status);
        NCM_FITS_ERROR (status);
      }

      need_flush = TRUE;
    }
    else if (self->file_cur_id > self->cur_id)
    {
      guint rows_to_add        = self->file_cur_id - self->cur_id;
      guint offset             = self->cur_id + 1 - self->first_id;
      gchar *stat              = NULL;
      NcmMSetCatalogSync smode = self->smode;

      self->smode = NCM_MSET_CATALOG_SYNC_DISABLE;

      for (i = 0; i < rows_to_add; i++)
      {
        NcmVector *row = ncm_vector_new (ncm_stats_vec_len (self->pstats));

        _ncm_mset_catalog_read_row (mcat, row, offset + i + 1);
        _ncm_mset_catalog_post_update (mcat, row);
        ncm_vector_free (row);
      }

      self->smode = smode;

      g_assert_cmpint (self->cur_id, ==, self->file_cur_id);

      if (self->rng != NULL)
      {
        fits_read_key_longstr (self->fptr, NCM_MSET_CATALOG_RNG_STAT_LABEL, &stat, NULL, &status);
        NCM_FITS_ERROR (status);

        g_clear_pointer (&self->rng_stat, g_free);
        self->rng_stat = g_strdup (stat);

        fits_free_memory (stat, &status);
        NCM_FITS_ERROR (status);

        ncm_rng_set_state (self->rng, self->rng_stat);
      }
    }
  }

  /*printf ("# Sync: status %d %d, %d %d\n", self->file_first_id, self->first_id, self->file_cur_id, self->cur_id);*/
  if (_ncm_mset_catalog_sync_markovian_id (mcat))
    need_flush = TRUE;

  /*printf ("# Sync: need flush %d\n", need_flush);*/
  if (need_flush)
    _ncm_mset_catalog_flush_file (mcat);

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_timed_sync:
 * @mcat: a #NcmMSetCatalog
 * @check: whether to check that the file and the memory agree
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

  ncm_vector_set_all (self->params_max, GSL_NEGINF);
  ncm_vector_set_all (self->params_min, GSL_POSINF);

  ncm_vector_clear (&self->bestfit_row);

  self->bestfit        = GSL_POSINF;
  self->post_lnnorm    = 0.0;
  self->post_lnnorm_up = FALSE;

  g_ptr_array_set_size (self->order_cat, 0);
  self->order_cat_sort = FALSE;
}

/**
 * ncm_mset_catalog_reset:
 * @mcat: a #NcmMSetCatalog
 *
 * Erases the rows of @mcat from memory and from its file, with all the statistics;
 * the configuration (model set, chains, file, random number generator) is kept.
 *
 */
void
ncm_mset_catalog_reset (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  ncm_mset_catalog_erase_data (mcat);

  ncm_stats_vec_reset (self->pstats, TRUE);
  ncm_stats_acorr_reset (self->acorr);

  if (self->nchains > 1)
  {
    guint i;

    for (i = 0; i < self->nchains; i++)
    {
      NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, i);

      ncm_stats_vec_reset (pstats, TRUE);
    }

    ncm_stats_vec_reset (self->mean_pstats, TRUE);
    ncm_stats_vec_reset (self->e_stats, TRUE);
    ncm_stats_vec_reset (self->e_mean_stats, TRUE);
    g_ptr_array_set_size (self->e_var_array, 0);
    g_array_set_size (self->accept_ratio, 0);
    self->naccepted = 0;
  }

  ncm_vector_set_all (self->params_max, GSL_NEGINF);
  ncm_vector_set_all (self->params_min, GSL_POSINF);

  ncm_vector_clear (&self->bestfit_row);

  self->bestfit        = GSL_POSINF;
  self->post_lnnorm    = 0.0;
  self->post_lnnorm_up = FALSE;

  g_ptr_array_set_size (self->order_cat, 0);
  self->order_cat_sort = FALSE;

  self->cur_id       = self->first_id - 1;
  self->markovian_id = self->first_id; /* a new chain: nothing non-Markovian yet */
#ifdef HAVE_CFITSIO
  self->file_cur_id = self->file_first_id - 1;
  _ncm_mset_catalog_close_file (mcat);

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_erase_data:
 * @mcat: a #NcmMSetCatalog
 *
 * Erases the rows of the file of @mcat; the rows in memory are kept.
 *
 */
void
ncm_mset_catalog_erase_data (NcmMSetCatalog *mcat)
{
#ifdef HAVE_CFITSIO
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (self->fptr != NULL)
  {
    gint status = 0;
    gint nrows  = self->file_cur_id - self->file_first_id + 1;

    if (nrows > 0)
    {
      fits_delete_rows (self->fptr, 1, nrows, &status);
      NCM_FITS_ERROR (status);

      self->file_cur_id = self->file_first_id - 1;
      _ncm_mset_catalog_flush_file (mcat);
    }
  }

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_set_m2lnp_var:
 * @mcat: a #NcmMSetCatalog
 * @p: $-2\ln(p)$ parameter index
 *
 * Sets @p as the $-2\ln(p)$ parameter index.
 *
 */
void
ncm_mset_catalog_set_m2lnp_var (NcmMSetCatalog *mcat, const gint p)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  self->m2lnp_var = p;
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
 * Gets the filename associated with @mcat.
 *
 * Returns: (transfer none): filename or NULL.
 */
const gchar *
ncm_mset_catalog_peek_filename (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->file;
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

  return (self->cur_id < self->first_id);
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
 * Returns: the index in @mcat of the row with id @t
 */
guint
ncm_mset_catalog_get_row_from_time (NcmMSetCatalog *mcat, gint t)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint row_n           = t - self->first_id;

  g_assert_cmpint (t, >=, self->first_id);
  g_assert_cmpuint (row_n, <, ncm_stats_vec_nitens (self->pstats));

  return row_n;
}

/**
 * ncm_mset_catalog_get_first_id:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the id of the first row in this catalog.
 */
gint
ncm_mset_catalog_get_first_id (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return self->first_id;
}

/**
 * ncm_mset_catalog_get_cur_id:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: the id of the last row added; the first id minus one when empty
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

  return ncm_mset_fparam_len (self->mset) + self->nadd_vals;
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

  if (i < self->nadd_vals)
    return g_ptr_array_index (self->add_vals_names, i);
  else
    return ncm_mset_fparam_name (self->mset, i - self->nadd_vals);
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

  if (i < self->nadd_vals)
    return g_ptr_array_index (self->add_vals_names, i);
  else
    return ncm_mset_fparam_full_name (self->mset, i - self->nadd_vals);
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

  if (i < self->nadd_vals)
    return g_ptr_array_index (self->add_vals_symbs, i);
  else
    return ncm_mset_fparam_symbol (self->mset, i - self->nadd_vals);
}

/**
 * ncm_mset_catalog_col_by_name:
 * @mcat: a #NcmMSetCatalog
 * @name: column name
 * @col_index: (out): column index
 *
 * Finds the column given by @name: a free-parameter name, an additional value name or a
 * column number.
 *
 * Returns: whether @name was found
 */
gboolean
ncm_mset_catalog_col_by_name (NcmMSetCatalog *mcat, const gchar *name, guint *col_index)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const NcmMSetPIndex *pi     = ncm_mset_fparam_get_pi_by_name (self->mset, name, NULL);

  if (pi == NULL)
  {
    if (g_ptr_array_find_with_equal_func (self->add_vals_names, name, g_str_equal, col_index))
    {
      return TRUE;
    }
    else
    {
      gchar *end_ptr  = NULL;
      glong add_param = strtol (name, &end_ptr, 10);

      if (name != end_ptr)
      {
        col_index[0] = add_param;

        return TRUE;
      }
      else
      {
        return FALSE;
      }
    }
  }
  else
  {
    col_index[0] = ncm_mset_fparam_get_fpi (self->mset, pi->mid, pi->pid) + ncm_mset_catalog_nadd_vals (mcat);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_mset_catalog_set_burnin:
 * @mcat: a #NcmMSetCatalog
 * @burnin: number of leading rows to drop
 *
 * Sets the number of leading rows of the file that are not read; the ids of the rows
 * read start after them. It applies only when rows are loaded from a file, so it must be
 * set before that, and is meant for analyzing a catalog.
 *
 */
void
ncm_mset_catalog_set_burnin (NcmMSetCatalog *mcat, glong burnin)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

#ifdef HAVE_CFITSIO

  if (self->fptr != NULL)
    g_error ("ncm_mset_catalog_set_burnin: cannot set burnin with an already loaded catalog");

#endif /* HAVE_CFITSIO */
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
 * iteration that used a non-Markovian move.
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
}

/**
 * ncm_mset_catalog_get_markovian_id:
 * @mcat: a #NcmMSetCatalog
 *
 * Returns: #NcmMSetCatalog:markovian-id, never below the first id.
 */
gint
ncm_mset_catalog_get_markovian_id (NcmMSetCatalog *mcat)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return MAX (self->markovian_id, self->first_id);
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
  const gint rows             = ncm_mset_catalog_get_markovian_id (mcat) - self->first_id;

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

/*
 * The autocorrelation accumulator forms unweighted lagged sums. A weighted catalog is
 * therefore not fed (see _ncm_mset_catalog_post_update) and every quantity derived from
 * the accumulator refuses rather than returning a number built from the wrong sums.
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

static void
_ncm_mset_catalog_post_update (NcmMSetCatalog *mcat, NcmVector *x)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint len             = ncm_stats_vec_len (self->pstats);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble p_i         = ncm_vector_get (x, i);
    const gdouble cur_max_p_i = ncm_vector_get (self->params_max, i);
    const gdouble cur_min_p_i = ncm_vector_get (self->params_min, i);

    ncm_vector_set (self->params_max, i, GSL_MAX (p_i, cur_max_p_i));
    ncm_vector_set (self->params_min, i, GSL_MIN (p_i, cur_min_p_i));
  }

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

  /* Counting the number of accepted steps */
  if (ncm_stats_vec_nitens (self->pstats) > self->nchains)
  {
    const gint current_index = ncm_stats_vec_nitens (self->pstats) - 1;
    const gint last_index    = current_index - self->nchains;
    NcmVector *last_x        = ncm_stats_vec_peek_row (self->pstats, last_index);

    if (ncm_vector_cmp2 (last_x, x, 1.0e-14, 0.0) > 0)
      self->naccepted++;
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

      ncm_stats_vec_append (self->e_mean_stats, e_mean, TRUE);
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

    g_ptr_array_add (self->order_cat, ncm_vector_ref (x));
    self->order_cat_sort = FALSE;
  }

  switch (self->smode)
  {
    case NCM_MSET_CATALOG_SYNC_DISABLE:
      break;
    case NCM_MSET_CATALOG_SYNC_AUTO:
      ncm_mset_catalog_sync (mcat, FALSE);
      break;
    case NCM_MSET_CATALOG_SYNC_TIMED:
      ncm_mset_catalog_timed_sync (mcat, FALSE);
      break;
    default:
      g_assert_not_reached ();
      break;
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
  NcmVector *row_i            = ncm_vector_new (ncm_stats_vec_len (self->pstats));
  va_list ap;
  guint i;

  va_start (ap, mset);

  for (i = 0; i < self->nadd_vals; i++)
  {
    gdouble val = va_arg (ap, gdouble);

    ncm_vector_set (row_i, i, val);
  }

  va_end (ap);

  ncm_mset_fparams_get_vector_offset (mset, row_i, self->nadd_vals);

  _ncm_mset_catalog_post_update (mcat, row_i);
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
  NcmVector *row_i            = ncm_vector_new (ncm_stats_vec_len (self->pstats));
  guint i;

  for (i = 0; i < self->nadd_vals; i++)
    ncm_vector_set (row_i, i, ax[i]);

  ncm_mset_fparams_get_vector_offset (mset, row_i, self->nadd_vals);

  _ncm_mset_catalog_post_update (mcat, row_i);
  ncm_vector_free (row_i);
}

/**
 * ncm_mset_catalog_add_from_vector:
 * @mcat: a #NcmMSetCatalog
 * @vals: a whole row: the additional values followed by the free parameters
 *
 * Adds the row @vals.
 *
 */
void
ncm_mset_catalog_add_from_vector (NcmMSetCatalog *mcat, NcmVector *vals)
{
  NcmVector *row_i = ncm_vector_dup (vals);

  _ncm_mset_catalog_post_update (mcat, row_i);
  ncm_vector_free (row_i);
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
  NcmVector *row_i            = ncm_vector_new (ncm_stats_vec_len (self->pstats));
  guint i;

  for (i = 0; i < self->nadd_vals; i++)
    ncm_vector_set (row_i, i, ax[i]);

  ncm_vector_memcpy2 (row_i, vals, self->nadd_vals, 0, ncm_vector_len (vals));

  _ncm_mset_catalog_post_update (mcat, row_i);
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
 * Sets the array of functions used to compute @mcat's extra columns, so it is
 * embedded alongside the model-set the next time @mcat's primary HDU is
 * written, making the catalog file self-sufficient (no companion functions
 * file needed to know what those columns mean). Must be called before
 * ncm_mset_catalog_set_file() creates a new catalog file.
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
 * Gets the array of functions used to compute @mcat's extra columns, either
 * set through ncm_mset_catalog_set_functions_array() or read back from the
 * catalog file's primary HDU.
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

static void
_ncm_mset_catalog_set_record (NcmMSetCatalog *mcat, gchar **field, const gchar *label, const gchar *what, const gchar *value)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  if (value == NULL)
  {
    g_clear_pointer (field, g_free);

    return;
  }

  if ((*field != NULL) && (strcmp (*field, value) == 0))
    return;

  g_clear_pointer (field, g_free);
  *field = g_strdup (value);
#ifdef HAVE_CFITSIO

  if ((self->fptr != NULL) && !self->readonly)
    _ncm_mset_catalog_write_record (self->fptr, label, what, *field);

#endif /* HAVE_CFITSIO */
}

/**
 * ncm_mset_catalog_set_sampler:
 * @mcat: a #NcmMSetCatalog
 * @sampler: (nullable): a description of the sampler about to fill the rows
 *
 * Records which sampler produces the rows, see #NcmMSetCatalog:sampler. It is written to
 * the file when there is one and it is never checked: unlike the run type, it does not
 * have to match what the file holds, and a catalog may be continued with another sampler.
 *
 */
void
ncm_mset_catalog_set_sampler (NcmMSetCatalog *mcat, const gchar *sampler)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_set_record (mcat, &self->sampler_str, NCM_MSET_CATALOG_SAMPLER_LABEL, "sampler", sampler);
}

/**
 * ncm_mset_catalog_set_sampler_options:
 * @mcat: a #NcmMSetCatalog
 * @options: (nullable): the sampler's tunable settings
 *
 * Records how the sampler is tuned, see #NcmMSetCatalog:sampler-options. Kept and written
 * like ncm_mset_catalog_set_sampler(), and never checked.
 *
 */
void
ncm_mset_catalog_set_sampler_options (NcmMSetCatalog *mcat, const gchar *options)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_set_record (mcat, &self->sampler_opts_str, NCM_MSET_CATALOG_SAMPLER_OPTS_LABEL, "sampler options", options);
}

/**
 * ncm_mset_catalog_set_initial_sampler:
 * @mcat: a #NcmMSetCatalog
 * @sampler: (nullable): a description of the sampler that drew the initial points
 *
 * Records which sampler drew the initial points, see #NcmMSetCatalog:initial-sampler.
 * Kept and written like ncm_mset_catalog_set_sampler(), and never checked.
 *
 */
void
ncm_mset_catalog_set_initial_sampler (NcmMSetCatalog *mcat, const gchar *sampler)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  _ncm_mset_catalog_set_record (mcat, &self->init_sampler_str, NCM_MSET_CATALOG_INIT_SAMPLER_LABEL, "initial points sampler", sampler);
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
 * Peeks the parameters #NcmStatsVec object.
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
 * Peeks the ensemble #NcmStatsVec object.
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
 * Peeks the chain #NcmStatsVec object.
 *
 * Returns: (transfer none): the #NcmStatsVec object.
 */
NcmStatsVec *
ncm_mset_catalog_peek_chain_pstats (NcmMSetCatalog *mcat, const guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);

  return g_ptr_array_index (self->chain_pstats, i);
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
 * Gets the @i-th row.
 *
 * Returns: (transfer none): the row with index @i or NULL if not available.
 */
NcmVector *
ncm_mset_catalog_peek_row (NcmMSetCatalog *mcat, guint i)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  guint nrows                 = ncm_stats_vec_nrows (self->pstats);

  if (i >= nrows)
    return NULL;
  else
    return ncm_stats_vec_peek_row (self->pstats, i);
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
    return ncm_stats_vec_peek_row (self->pstats, ncm_stats_vec_nitens (self->pstats) - 1);
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
    return ncm_stats_vec_peek_row (self->e_mean_stats, ncm_stats_vec_nitens (self->e_mean_stats) - 1);
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

    return ncm_stats_vec_peek_row (self->e_mean_stats, t);
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

static gdouble
_ncm_mset_catalog_get_post_lnnorm_elipsoid (NcmMSetCatalog *mcat, gdouble *post_lnnorm_sd)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  const guint fparams_len     = ncm_mset_fparams_len (self->mset);
  const guint cat_len         = ncm_mset_catalog_len (mcat);
  NcmMatrix *cov              = NULL;
  NcmVector *mean             = NULL;
  NcmVector *v                = ncm_vector_new (fparams_len);
  gdouble level               = 0.50;
  gdouble R2_cut              = gsl_cdf_chisq_Pinv (level, fparams_len);
  gdouble lnNorma             = 0.0;
  gdouble s                   = 0.0;
  gdouble c                   = 0.0;
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
    ncm_vector_free (v);

    return 0.0;
  }

  lnNorma = 0.5 * (fparams_len * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (cov)) + log (level);

  for (i = 0; i < cat_len; i++)
  {
    NcmVector *row_i      = ncm_mset_catalog_peek_row (mcat, i);
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

    if (m2lnp_i > R2_cut)
      continue;

    e_i = exp (0.5 * ((m2lnL_i - self->bestfit) - m2lnp_i));
    t   = s + e_i;
    c  += (s >= e_i) ? ((s - t) + e_i) : ((e_i - t) + s);
    s   = t;
  }

  post_lnnorm       = -(log ((s + c) / cat_len) - lnNorma + 0.5 * self->bestfit);
  post_lnnorm_sd[0] = GSL_NAN;

  ncm_vector_clear (&mean);
  ncm_matrix_clear (&cov);
  ncm_vector_free (v);

  return post_lnnorm;
}

static gdouble
_ncm_mset_catalog_get_post_lnnorm_sum (NcmMSetCatalog *mcat, NcmVector *mean, NcmMatrix *cov, const gdouble lnNorma, gdouble *lnnorm_sd)
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


    printf ("# BS LNNORM: % 22.15g : % 22.15g % 22.15g % 22.15g % 22.15g\n",
            -(log ((s + c) / cat_len) - lnNorma + 0.5 * self->bestfit),
            -(log (ncm_stats_vec_get_mean (slnnorm, 0)) - lnNorma + 0.5 * self->bestfit),
            ncm_stats_vec_get_mean (slnnorm, 1),
            ncm_stats_vec_get_sd (slnnorm, 1),
            ncm_stats_vec_get_sd (slnnorm, 1) / sqrt (j + 1.0)
    );


    g_array_unref (bs_array);

    if (j > 10)
    {
      const gdouble lnnorm_mean    = ncm_stats_vec_get_sd (slnnorm, 0);
      const gdouble lnnorm_mean_sd = ncm_stats_vec_get_sd (slnnorm, 0) / sqrt (j + 1.0);

      if (lnnorm_mean_sd / lnnorm_mean < reltol)
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
    post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_sum (mcat, mean, cov, lnNorma, post_lnnorm_sd);

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
 * with the rows' mean and $0.8^2$ times their covariance, normalized to the box by a
 * Monte Carlo estimate of its mass there; the result is kept until rows are added.
 * Without #NcmMSetCatalog:m2lnp-var it warns and returns zero.
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
    NcmMSetCatalogPostNormMethod method = NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX;

    switch (method)
    {
      case NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX:
        self->post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_hyperbox (mcat, FALSE, &self->post_lnnorm_sd);
        break;
      case NCM_MSET_CATALOG_POST_LNNORM_METHOD_HYPERBOX_BS:
        self->post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_hyperbox (mcat, TRUE, &self->post_lnnorm_sd);
        break;
      case NCM_MSET_CATALOG_POST_LNNORM_METHOD_ELLIPSOID:
        self->post_lnnorm = _ncm_mset_catalog_get_post_lnnorm_elipsoid (mcat, &self->post_lnnorm_sd);
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

  if (!self->order_cat_sort)
  {
    g_ptr_array_sort_with_data (self->order_cat, &_ncm_mset_catalog_double_compare, GINT_TO_POINTER (self->m2lnp_var));
    self->order_cat_sort = TRUE;
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
    const gdouble m2lnp_i = ncm_vector_get (g_ptr_array_index (self->order_cat, i), self->m2lnp_var);

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

  m2lnL = ncm_vector_get (g_ptr_array_index (self->order_cat, nitens), self->m2lnp_var);

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
    *mean = ncm_vector_new (ncm_stats_vec_len (self->pstats) - self->nadd_vals);

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
    *cov = ncm_matrix_new (ncm_stats_vec_len (self->pstats) - self->nadd_vals, ncm_stats_vec_len (self->pstats) - self->nadd_vals);

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
    *cov = ncm_matrix_new (ncm_stats_vec_len (self->pstats), ncm_stats_vec_len (self->pstats));

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
  const guint params_len      = ncm_mset_fparam_len (self->mset) + self->nadd_vals;
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
      ncm_stats_acorr_update (rows, ncm_stats_vec_peek_row (self->pstats, i));

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
  const guint free_params_len = ncm_stats_vec_len (self->pstats) - self->nadd_vals;
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
    NcmVector *row = ncm_stats_vec_peek_row (self->pstats, k);

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

/**
 * ncm_mset_catalog_trim:
 * @mcat: a #NcmMSetCatalog
 * @tc: number of iterations to drop
 * @thin: thinning factor, at least one
 *
 * Drops the first @tc iterations, the first $t_c n_\mathrm{chains}$ rows, and keeps one
 * iteration in every @thin of the rest, all chains of it. The first id moves forward by
 * $t_c n_\mathrm{chains}$ and the markovian id follows the rows it pointed to. The catalog
 * must hold complete ensembles and must not be read-only.
 *
 * When a file is attached it is renamed to `<file>.<n>.bak`, with the first free $n$,
 * and the trimmed catalog is written to a new file under the original name. Nothing is
 * done when @tc is zero and @thin is one.
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
    GPtrArray *rows = ncm_stats_vec_dup_saved_x (self->pstats);
    gchar *file     = g_strdup (ncm_mset_catalog_peek_filename (mcat));
    guint t, tf;

    if (file != NULL)
      ncm_mset_catalog_set_file (mcat, NULL);

    /* Iterations before the Markovian chain, in the old numbering; after the cut and the
     * thinning the same rows are the first ceil ((mark_iters - tc) / thin) kept iterations. */
    const gint mark_iters = (ncm_mset_catalog_get_markovian_id (mcat) - self->first_id + (gint) self->nchains - 1) / (gint) self->nchains;
    const gint kept_iters = (mark_iters > (gint) tc) ? (mark_iters - (gint) tc + (gint) thin - 1) / (gint) thin : 0;

    ncm_mset_catalog_reset (mcat);
    ncm_mset_catalog_set_first_id (mcat, self->first_id + tc * self->nchains);

    g_assert_cmpint (rows->len % self->nchains, ==, 0);

    tf = rows->len / self->nchains;

    for (t = tc; t < tf; t += thin)
    {
      guint i;

      for (i = 0; i < self->nchains; i++)
      {
        NcmVector *row_t = g_ptr_array_index (rows, t * self->nchains + i);

        _ncm_mset_catalog_post_update (mcat, row_t);
      }
    }

    self->markovian_id = MIN (self->first_id + kept_iters * (gint) self->nchains, self->cur_id + 1);

    if (file != NULL)
    {
      guint bak_n = 0;

      while (TRUE)
      {
        gchar *bak_name = g_strdup_printf ("%s.%d.bak", file, bak_n);

        if (!g_file_test (bak_name, G_FILE_TEST_EXISTS))
        {
          g_rename (file, bak_name);
          ncm_mset_catalog_set_file (mcat, file);
          g_free (bak_name);

          break;
        }

        g_free (bak_name);
        bak_n++;
      }
    }

    g_ptr_array_unref (rows);
    g_free (file);
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
 * @out_file: output filename
 *
 * Removes the rows whose free parameters lie outside the bounds of the #NcmMSet of the
 * catalog. The catalog is renumbered as a single chain from the same first id, with the
 * markovian id after the kept rows that preceded it, and written to @out_file; the file
 * attached before is left as it was.
 *
 * Returns: the number of rows removed.
 */
guint
ncm_mset_catalog_trim_oob (NcmMSetCatalog *mcat, const gchar *out_file)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
  GPtrArray *rows             = ncm_stats_vec_dup_saved_x (self->pstats);
  gchar *file                 = g_strdup (ncm_mset_catalog_peek_filename (mcat));
  guint i, ndel;

  if (file != NULL)
    ncm_mset_catalog_set_file (mcat, NULL);

  const guint mark_rows = ncm_mset_catalog_get_markovian_id (mcat) - self->first_id;
  guint kept_before     = 0;

  ncm_mset_catalog_reset (mcat);
  self->nchains = 1;
  ncm_mset_catalog_set_first_id (mcat, self->first_id);

  ndel = 0;

  for (i = 0; i < rows->len; i++)
  {
    NcmVector *row_t = g_ptr_array_index (rows, i);

    if (ncm_mset_fparam_valid_bounds_offset (self->mset, row_t, self->nadd_vals))
    {
      _ncm_mset_catalog_post_update (mcat, row_t);

      if (i < mark_rows)
        kept_before++;
    }
    else
    {
      ndel++;
    }
  }

  /* Rows are renumbered: the Markovian chain starts after the kept rows that preceded it. */
  self->markovian_id = self->first_id + kept_before;

  ncm_mset_catalog_set_file (mcat, out_file);
  ncm_mset_catalog_sync (mcat, TRUE);

  g_ptr_array_unref (rows);
  g_free (file);

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
 * and the catalog is written to a new file under the original name.
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

  const gint last_t = ncm_stats_vec_nrows (self->pstats);

  if (last_t > 0)
  {
    const guint nrows_to_remove = (last_t % self->nchains == 0) ? self->nchains : (last_t % self->nchains);
    const guint new_len         = last_t  - nrows_to_remove;

    g_assert_cmpuint (new_len, <, last_t);

    {
      GPtrArray *rows = ncm_stats_vec_dup_saved_x (self->pstats);
      gchar *file     = g_strdup (ncm_mset_catalog_peek_filename (mcat));
      guint t;

      const gint markovian_id = ncm_mset_catalog_get_markovian_id (mcat);

      ncm_mset_catalog_set_file (mcat, NULL);
      ncm_mset_catalog_reset (mcat);

      for (t = 0; t < new_len; t++)
      {
        NcmVector *row_t = g_ptr_array_index (rows, t);

        _ncm_mset_catalog_post_update (mcat, row_t);
      }

      /* The removed rows may have been the ones the chain was to start at: the id then sits
       * at the row after the new last one, i.e. the phase is still open. */
      self->markovian_id = MIN (markovian_id, self->cur_id + 1);

      if (file != NULL)
      {
        guint bak_n = 0;

        while (TRUE)
        {
          gchar *bak_name = g_strdup_printf ("%s.%d.bak", file, bak_n);

          if (!g_file_test (bak_name, G_FILE_TEST_EXISTS))
          {
            g_rename (file, bak_name);
            ncm_mset_catalog_set_file (mcat, file);
            g_free (bak_name);

            break;
          }

          g_free (bak_name);
          bak_n++;
        }
      }

      g_ptr_array_unref (rows);
      g_free (file);
    }
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
 * the columns, see ncm_stats_vec_max_ess_time(). A catalog with one chain is read row
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
  NcmStatsVec *pstats         = (self->nchains == 1) ? self->pstats : self->e_mean_stats;
  const gint last_t           = ncm_stats_vec_nrows (pstats);
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

    return 0;
  }

  esss = ncm_stats_vec_max_ess_time (pstats, ntests, &bindex, &wp, &wp_order, &wp_ess);

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
_fonempval (gdouble v_i, guint i, gpointer user_data)
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
 * ncm_stats_vec_heidel_diag(). A catalog with one chain is read row by row; one with
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
  NcmStatsVec *pstats         = (self->nchains == 1) ? self->pstats : self->e_mean_stats;
  const gdouble pvalue_lef    = (pvalue == 0.0) ? NCM_STATS_VEC_HEIDEL_PVAL_COR (0.05, ncm_mset_fparams_len (self->mset)) : pvalue;
  gint bindex                 = 0;
  guint wp = 0, wp_order = 0;
  gdouble wp_pvalue  = 0.0;
  const guint last_t = ncm_stats_vec_nrows (pstats);
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

    return 0;
  }

  pvals = ncm_stats_vec_heidel_diag (pstats, ntests, pvalue_lef, &bindex, &wp, &wp_order, &wp_pvalue);

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
    ncm_vector_log_vals_func (pvals, "# NcmMSetCatalog: - pvalues:                  ", "%5.2f%%", &_fonempval, NULL);
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
 * ncm_stats_vec_estimate_const_break(). A catalog with one chain is read row by row;
 * one with several chains is read through the ensemble means, so $t_0$ is in iterations.
 *
 * Returns: $t_0$, rounded up.
 */
guint
ncm_mset_catalog_calc_const_break (NcmMSetCatalog *mcat, guint p, NcmFitRunMsgs mtype)
{
  NcmMSetCatalogPrivate *self = ncm_mset_catalog_get_instance_private (mcat);
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

  if (self->nchains > 1)
    tc = ceil (ncm_stats_vec_estimate_const_break (self->e_mean_stats, p));
  else
    tc = ceil (ncm_stats_vec_estimate_const_break (self->pstats, p));

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

/**
 * ncm_mset_catalog_max_ess_time_by_chain:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @max_ess: (out) (optional): the smallest effective sample size of the worst chain
 * @mtype: #NcmFitRunMsgs log level
 *
 * Applies ncm_stats_vec_max_ess_time() to each chain on its own. The worst chain is the
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
    gint tbindex = -1;
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
      gdouble wp_ess      = 0.0;
      NcmStatsVec *pstats = g_ptr_array_index (self->chain_pstats, i);
      NcmVector *esss     = ncm_stats_vec_max_ess_time (pstats, ntests, &bindex, &wp, &wp_order, &wp_ess);
      guint size          = ncm_stats_vec_nitens (pstats);

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

      if (bindex > tbindex)
      {
        tbindex   = bindex;
        twp       = wp;
        twp_order = wp_order;
        twp_ess   = wp_ess;
        ti        = i;
      }

      ncm_vector_free (esss);
    }

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
      ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4d\n", tbindex);
      ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", size);
      ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", size - tbindex);
      ncm_message ("# NcmMSetCatalog: - worst parameter:          %-4u\n", twp);
      ncm_message ("# NcmMSetCatalog: - worst parameter order:    %-4u\n", twp_order);
      ncm_message ("# NcmMSetCatalog: - worst parameter ESS:      %-.2f\n", twp_ess);
    }

    if (max_ess != NULL)
      max_ess[0] = twp_ess;

    return tbindex;
  }
}

/**
 * ncm_mset_catalog_heidel_diag_by_chain:
 * @mcat: a #NcmMSetCatalog
 * @ntests: number of starting times tested, 10 when zero
 * @pvalue: required p-value, the default of ncm_mset_catalog_calc_heidel_diag() when zero
 * @wp_pvalue: (out) (optional): the worst value of the worst chain
 * @mtype: #NcmFitRunMsgs log level
 *
 * Applies ncm_stats_vec_heidel_diag() to each chain on its own. The worst chain is one
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
    gint tbindex             = -1;
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
      gdouble lwp_pvalue  = 0.0;
      NcmStatsVec *pstats = (self->nchains == 1) ? self->pstats : g_ptr_array_index (self->chain_pstats, i);
      NcmVector *pvals    = ncm_stats_vec_heidel_diag (pstats, ntests, pvalue_lef, &bindex, &wp, &wp_order, &lwp_pvalue);
      guint size          = ncm_stats_vec_nitens (pstats);

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

      if ((i == 0) || ((tbindex >= 0) && (bindex > tbindex)) || ((bindex == -1) && (tbindex >= 0)) || ((bindex == -1) && (lwp_pvalue > twp_pvalue)))
      {
        tbindex    = bindex;
        twp        = wp;
        twp_order  = wp_order;
        twp_pvalue = lwp_pvalue;
        ti         = i;
      }

      ncm_vector_free (pvals);
    }

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

      if (tbindex == -1)
      {
        ncm_message ("# NcmMSetCatalog: - best cutoff time:         not satisfied!\n");
        ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", size);
        ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", 0);
      }
      else
      {
        ncm_message ("# NcmMSetCatalog: - best cutoff time:         %-4u\n", tbindex);
        ncm_message ("# NcmMSetCatalog: - total number of points:   %-4u\n", size);
        ncm_message ("# NcmMSetCatalog: - number of points left:    %-4u\n", size - tbindex);
      }

      ncm_message ("# NcmMSetCatalog: - worst parameter:          %-4u\n", twp);
      ncm_message ("# NcmMSetCatalog: - worst parameter order:    %-4u\n", twp_order);
      ncm_message ("# NcmMSetCatalog: - worst parameter pvalue:   %6.2f%%\n", (1.0 - twp_pvalue) * 100.0);
    }

    if (wp_pvalue != NULL)
      wp_pvalue[0] = twp_pvalue;

    return (tbindex >= 0) ? tbindex : 0;
  }
}


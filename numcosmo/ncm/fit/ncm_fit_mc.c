/***************************************************************************
 *            ncm_fit_mc.c
 *
 *  Sat December 01 17:19:10 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmFitMC:
 *
 * Monte Carlo study of a best-fit estimator.
 *
 * Each realization resamples the data of the fit's #NcmDataset, runs the fit from the
 * fiducial parameters and stores $-2\ln L$, the best-fit free parameters and the values
 * of the optional functions (#NcmFitMC:function-array) as a row of a #NcmMSetCatalog.
 * The resampling is chosen by #NcmFitMCResampleType: from the fiducial model
 * (#NcmFitMC:fiducial, ncm_dataset_resample()), or a bootstrap of the data, of each
 * #NcmData separately or of all of them together. The catalog then holds the
 * distribution of the estimator, e.g. its mean and covariance (ncm_fit_mc_mean_covar()).
 *
 * A run is started with ncm_fit_mc_start_run(), extended with ncm_fit_mc_run() or
 * ncm_fit_mc_run_lre() and closed with ncm_fit_mc_end_run(). With a data file the
 * catalog is saved during the run and a restarted run continues from its last row.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_fit_mc.h"
#include "ncm/fit/ncm_mset_catalog_storage_fits.h"
#include "ncm/fit/ncm_mset_catalog_storage_memory.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"
#include "ncm/core/ncm_func_eval.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_statistics_double.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_FIT,
  PROP_RTYPE,
  PROP_FIDUC,
  PROP_MTYPE,
  PROP_USE_THREADS,
  PROP_KEEP_ORDER,
  PROP_STORAGE,
  PROP_FUNC_ARRAY,
};

struct _NcmFitMC
{
  /*< private >*/
  GObject parent_instance;
  NcmFitMCResample resample;
  NcmFit *fit;
  NcmMSet *fiduc;
  NcmMSetCatalog *mcat;
  NcmMSetCatalogStorage *storage;
  NcmRNG *rng;
  NcmFitRunMsgs mtype;
  NcmFitMCResampleType rtype;
  NcmVector *bf;
  NcmTimer *nt;
  NcmSerialize *ser;
  gboolean use_threads;
  guint n;
  NcmMemoryPool *mp;
  gint write_index;
  gint cur_sample_id;
  gboolean started;
  gboolean keep_order;
  NcmObjArray *func_oa;
  gchar *func_oa_file;
  guint nadd_vals;
  guint theta_len;
  gboolean constructed;
};

G_DEFINE_TYPE (NcmFitMC, ncm_fit_mc, G_TYPE_OBJECT)

static void
ncm_fit_mc_init (NcmFitMC *mc)
{
  mc->fit           = NULL;
  mc->fiduc         = NULL;
  mc->mcat          = NULL;
  mc->storage       = NULL;
  mc->rng           = NULL;
  mc->mtype         = NCM_FIT_RUN_MSGS_NONE;
  mc->rtype         = NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN;
  mc->bf            = NULL;
  mc->nt            = ncm_timer_new ();
  mc->ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  mc->use_threads   = FALSE;
  mc->n             = 0;
  mc->cur_sample_id = -1; /* Represents that no samples were calculated yet. */
  mc->write_index   = 0;
  mc->started       = FALSE;
  mc->keep_order    = FALSE;
  mc->func_oa       = NULL;
  mc->func_oa_file  = NULL;
  mc->nadd_vals     = 0;
  mc->theta_len     = 0;
  mc->constructed   = FALSE;
}

static void
ncm_fit_mc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMC *mc = NCM_FIT_MC (object);

  g_return_if_fail (NCM_IS_FIT_MC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_assert (mc->fit == NULL);
      mc->fit = g_value_dup_object (value);

      break;
    case PROP_RTYPE:
      mc->rtype = g_value_get_enum (value);

      break;
    case PROP_FIDUC:
      ncm_fit_mc_set_fiducial (mc, g_value_get_object (value));
      break;
    case PROP_MTYPE:
      ncm_fit_mc_set_mtype (mc, g_value_get_enum (value));
      break;
    case PROP_USE_THREADS:
      ncm_fit_mc_set_use_threads (mc, g_value_get_boolean (value));
      break;
    case PROP_KEEP_ORDER:
      ncm_fit_mc_keep_order (mc, g_value_get_boolean (value));
      break;
    case PROP_STORAGE:
      ncm_fit_mc_set_storage (mc, g_value_get_object (value));
      break;
    case PROP_FUNC_ARRAY:
    {
      ncm_obj_array_clear (&mc->func_oa);
      mc->func_oa = g_value_dup_boxed (value);

      if (mc->func_oa != NULL)
      {
        guint i;

        for (i = 0; i < mc->func_oa->len; i++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (mc->func_oa, i));

          g_assert (NCM_IS_MSET_FUNC (func));
          g_assert (ncm_mset_func_is_scalar (func));
          g_assert (ncm_mset_func_is_const (func));
        }
      }

      break;
    }

    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_fit_mc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitMC *mc = NCM_FIT_MC (object);

  g_return_if_fail (NCM_IS_FIT_MC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, mc->fit);
      break;
    case PROP_RTYPE:
      g_value_set_enum (value, mc->rtype);
      break;
    case PROP_FIDUC:
      g_value_set_object (value, mc->fiduc);
      break;
    case PROP_MTYPE:
      g_value_set_enum (value, mc->mtype);
      break;
    case PROP_USE_THREADS:
      g_value_set_boolean (value, mc->use_threads);
      break;
    case PROP_KEEP_ORDER:
      g_value_set_boolean (value, mc->keep_order);
      break;
    case PROP_STORAGE:
      g_value_set_object (value, mc->storage);
      break;
    case PROP_FUNC_ARRAY:
      g_value_set_boxed (value, mc->func_oa);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_fit_mc_apply_rtype (NcmFitMC *mc);

static void
ncm_fit_mc_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->constructed (object);
  {
    NcmFitMC *mc           = NCM_FIT_MC (object);
    NcmMSet *mset          = ncm_fit_peek_mset (mc->fit);
    const guint nfuncs     = (mc->func_oa != NULL) ? mc->func_oa->len : 0;
    const guint nadd_vals  = mc->nadd_vals = nfuncs + 1;
    const guint fparam_len = ncm_mset_fparam_len (mset);

    mc->theta_len = fparam_len + nadd_vals;

    _ncm_fit_mc_apply_rtype (mc);
    mc->constructed = TRUE;
  }
}

static void
ncm_fit_mc_dispose (GObject *object)
{
  NcmFitMC *mc = NCM_FIT_MC (object);

  ncm_fit_clear (&mc->fit);
  ncm_mset_clear (&mc->fiduc);
  ncm_vector_clear (&mc->bf);
  ncm_timer_clear (&mc->nt);
  ncm_serialize_clear (&mc->ser);
  ncm_mset_catalog_clear (&mc->mcat);
  ncm_rng_clear (&mc->rng);
  ncm_mset_catalog_storage_clear (&mc->storage);

  ncm_obj_array_clear (&mc->func_oa);
  g_clear_pointer (&mc->func_oa_file, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->dispose (object);
}

static void
ncm_fit_mc_finalize (GObject *object)
{
  /* NcmFitMC *mc = NCM_FIT_MC (object); */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->finalize (object);
}

static void
ncm_fit_mc_class_init (NcmFitMCClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_fit_mc_set_property;
  object_class->get_property = &ncm_fit_mc_get_property;
  object_class->constructed  = &ncm_fit_mc_constructed;
  object_class->dispose      = &ncm_fit_mc_dispose;
  object_class->finalize     = &ncm_fit_mc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RTYPE,
                                   g_param_spec_enum ("rtype",
                                                      NULL,
                                                      "Monte Carlo run type",
                                                      NCM_TYPE_FIT_MC_RESAMPLE_TYPE, NCM_FIT_MC_RESAMPLE_FROM_MODEL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FIDUC,
                                   g_param_spec_object ("fiducial",
                                                        NULL,
                                                        "Fiducial model to sample from",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MTYPE,
                                   g_param_spec_enum ("mtype",
                                                      NULL,
                                                      "Run messages type",
                                                      NCM_TYPE_FIT_RUN_MSGS, NCM_FIT_RUN_MSGS_SIMPLE,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_THREADS,
                                   g_param_spec_boolean ("use-threads",
                                                         NULL,
                                                         "Whether to use OpenMP threads (real thread count is controlled by OMP_NUM_THREADS)",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_KEEP_ORDER,
                                   g_param_spec_boolean ("keep-order",
                                                         NULL,
                                                         "Whether to keep the catalog in order of sampling under multi-threaded runs",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitMC:storage:
   *
   * The #NcmMSetCatalogStorage of the catalog that ncm_fit_mc_start_run() builds, see
   * ncm_fit_mc_set_storage(); without one the catalog is kept in a
   * #NcmMSetCatalogStorageMemory.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_STORAGE,
                                   g_param_spec_object ("storage",
                                                        NULL,
                                                        "Catalog storage",
                                                        NCM_TYPE_MSET_CATALOG_STORAGE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FUNC_ARRAY,
                                   g_param_spec_boxed ("function-array",
                                                       NULL,
                                                       "Functions array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_fit_mc_new:
 * @fit: a #NcmFit
 * @rtype: a #NcmFitMCResampleType
 * @mtype: a #NcmFitRunMsgs
 *
 * Creates a #NcmFitMC for @fit with resampling @rtype and messages @mtype. The
 * fiducial model is a copy of the model set of @fit.
 *
 * Returns: (transfer full): a new #NcmFitMC.
 */
NcmFitMC *
ncm_fit_mc_new (NcmFit *fit, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype)
{
  NcmFitMC *mc = g_object_new (NCM_TYPE_FIT_MC,
                               "fit", fit,
                               "rtype", rtype,
                               "mtype", mtype,
                               NULL);

  return mc;
}

/**
 * ncm_fit_mc_new_funcs_array:
 * @fit: a #NcmFit
 * @rtype: a #NcmFitMCResampleType
 * @mtype: a #NcmFitRunMsgs
 * @funcs_array: a #NcmObjArray of scalar constant #NcmMSetFunc
 *
 * Creates a #NcmFitMC as ncm_fit_mc_new(), with one more catalog column per function
 * of @funcs_array, evaluated at each best fit.
 *
 * Returns: (transfer full): a new #NcmFitMC.
 */
NcmFitMC *
ncm_fit_mc_new_funcs_array (NcmFit *fit, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype, NcmObjArray *funcs_array)
{
  NcmFitMC *mc = g_object_new (NCM_TYPE_FIT_MC,
                               "fit", fit,
                               "rtype", rtype,
                               "mtype", mtype,
                               "function-array", funcs_array,
                               NULL);

  return mc;
}

/**
 * ncm_fit_mc_free:
 * @mc: a #NcmFitMC
 *
 * Decreases the reference count of @mc.
 *
 */
void
ncm_fit_mc_free (NcmFitMC *mc)
{
  g_object_unref (mc);
}

/**
 * ncm_fit_mc_clear:
 * @mc: a #NcmFitMC
 *
 * Decreases the reference count of *@mc and sets it to %NULL.
 */
void
ncm_fit_mc_clear (NcmFitMC **mc)
{
  g_clear_object (mc);
}

static void ncm_fit_mc_intern_skip (NcmFitMC *mc, guint n);

/**
 * ncm_fit_mc_set_storage:
 * @mc: a #NcmFitMC
 * @storage: a #NcmMSetCatalogStorage
 *
 * Makes @storage the storage of the catalog that ncm_fit_mc_start_run() builds: the
 * catalog's description is set on it then, a #NcmMSetCatalogStorageFits being created
 * when its file does not exist and continued when it does. Once the catalog is built
 * the storage cannot change, so calling it after the first ncm_fit_mc_start_run()
 * aborts.
 *
 */
void
ncm_fit_mc_set_storage (NcmFitMC *mc, NcmMSetCatalogStorage *storage)
{
  g_assert_nonnull (storage);

  if (mc->mcat != NULL)
    g_error ("ncm_fit_mc_set_storage: the catalog is built by ncm_fit_mc_start_run() and its storage cannot change.");

  ncm_mset_catalog_storage_clear (&mc->storage);
  mc->storage = ncm_mset_catalog_storage_ref (storage);
}

/*
 * Builds the catalog at the first start_run: a description of the free parameters, the
 * -2ln(L) column and the function columns, on a FITS storage of the data file (created or
 * continued) or a memory storage, and gives it the RNG set before.
 */
static void
_ncm_fit_mc_build_catalog (NcmFitMC *mc, GError **error)
{
  NcmMSet *mset          = ncm_fit_peek_mset (mc->fit);
  const guint nfuncs     = (mc->func_oa != NULL) ? mc->func_oa->len : 0;
  const GEnumValue *eval = ncm_cfg_enum_get_value (NCM_TYPE_FIT_MC_RESAMPLE_TYPE, mc->rtype);
  NcmMSetCatalogDesc *desc;
  gchar **names;
  gchar **symbols;
  guint k;

  if (mc->mcat != NULL)
    return;

  names   = g_new (gchar *, nfuncs + 2);
  symbols = g_new (gchar *, nfuncs + 2);

  names[0]   = g_strdup (NCM_MSET_CATALOG_M2LNL_COLNAME);
  symbols[0] = g_strdup (NCM_MSET_CATALOG_M2LNL_SYMBOL);

  for (k = 0; k < nfuncs; k++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (mc->func_oa, k));

    g_assert (NCM_IS_MSET_FUNC (func));

    names[1 + k]   = g_strdup (ncm_mset_func_peek_uname (func));
    symbols[1 + k] = g_strdup (ncm_mset_func_peek_usymbol (func));
  }

  names[1 + k]   = NULL;
  symbols[1 + k] = NULL;

  desc = ncm_mset_catalog_desc_new (mset, 1, FALSE, (const gchar * const *) names, (const gchar * const *) symbols,
                                    0, NCM_MSET_CATALOG_RTYPE_MC, error);
  NCM_UTIL_ON_ERROR_RETURN (error, { g_strfreev (names);
                                     g_strfreev (symbols);
                            }, );

  if (mc->storage == NULL)
    mc->storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_memory_new ());

  ncm_mset_catalog_storage_set_desc (mc->storage, desc, error);
  NCM_UTIL_ON_ERROR_RETURN (error, { ncm_mset_catalog_desc_free (desc);
                                     g_strfreev (names);
                                     g_strfreev (symbols);
                            }, );

  mc->mcat = ncm_mset_catalog_new_from_storage (mc->storage, 0);
  ncm_mset_catalog_set_functions_array (mc->mcat, mc->func_oa);
  ncm_mset_catalog_set_sampler (mc->mcat, eval->value_nick);

  if (mc->rng != NULL)
  {
    NcmRNG *cat_rng = ncm_mset_catalog_peek_rng (mc->mcat);

    if (cat_rng == NULL)
      ncm_mset_catalog_set_rng (mc->mcat, mc->rng);
    else if ((strcmp (ncm_rng_get_algo (cat_rng), ncm_rng_get_algo (mc->rng)) != 0) ||
             (ncm_rng_get_seed (cat_rng) != ncm_rng_get_seed (mc->rng)))
      ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG,
                                  "ncm_fit_mc_start_run: catalog `%s' holds an RNG (%s, seed %lu) and another was set (%s, seed %lu); "
                                  "do not set an RNG when resuming from a file",
                                  ncm_mset_catalog_peek_filename (mc->mcat), ncm_rng_get_algo (cat_rng), ncm_rng_get_seed (cat_rng),
                                  ncm_rng_get_algo (mc->rng), ncm_rng_get_seed (mc->rng));

    ncm_rng_clear (&mc->rng);
  }

  ncm_mset_catalog_desc_free (desc);
  g_strfreev (names);
  g_strfreev (symbols);
}

/* Opens a recording of the catalog's storage named by the resampling type, closing one left open. */
static void
_ncm_fit_mc_start_recording (NcmFitMC *mc, GError **error)
{
  NcmMSetCatalogStorage *storage = ncm_mset_catalog_peek_storage (mc->mcat);
  const GEnumValue *eval         = ncm_cfg_enum_get_value (NCM_TYPE_FIT_MC_RESAMPLE_TYPE, mc->rtype);

  if (ncm_mset_catalog_storage_is_recording (storage))
  {
    ncm_mset_catalog_storage_end_recording (storage, error);
    NCM_UTIL_ON_ERROR_RETURN (error, , );
  }

  ncm_mset_catalog_storage_start_recording (storage, eval->value_nick, NULL, error);
}

/* Closes the recording of the catalog's storage, writing what it holds; none open is fine. */
static void
_ncm_fit_mc_end_recording (NcmFitMC *mc)
{
  NcmMSetCatalogStorage *storage = ncm_mset_catalog_peek_storage (mc->mcat);
  GError *error                  = NULL;

  if (!ncm_mset_catalog_storage_is_recording (storage))
    return;

  ncm_mset_catalog_storage_end_recording (storage, &error);

  if (error != NULL)
    g_error ("ncm_fit_mc_end_run: %s", error->message);
}

/**
 * ncm_fit_mc_set_mtype:
 * @mc: a #NcmFitMC
 * @mtype: a #NcmFitRunMsgs
 *
 * Sets the run messages type of @mc to @mtype.
 *
 */
void
ncm_fit_mc_set_mtype (NcmFitMC *mc, NcmFitRunMsgs mtype)
{
  mc->mtype = mtype;
}

static void _ncm_fit_mc_resample_bstrap (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng);

static gint
_ncm_fit_mc_resample (NcmFitMC *mc, NcmFit *fit)
{
  NcmLikelihood *lh = ncm_fit_peek_likelihood (fit);
  NcmDataset *dset  = ncm_likelihood_peek_dataset (lh);

  mc->resample (dset, mc->fiduc, ncm_mset_catalog_peek_rng (mc->mcat));
  mc->cur_sample_id++;

  return mc->cur_sample_id;
}

/*
 * Sets the resampling function of NcmFitMC:rtype and the bootstrap mode of the dataset.
 */
static void
_ncm_fit_mc_apply_rtype (NcmFitMC *mc)
{
  NcmLikelihood *lh                = ncm_fit_peek_likelihood (mc->fit);
  NcmDataset *dset                 = ncm_likelihood_peek_dataset (lh);
  const NcmFitMCResampleType rtype = mc->rtype;

  g_assert_cmpint (rtype, <, NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN);

  switch (rtype)
  {
    case NCM_FIT_MC_RESAMPLE_FROM_MODEL:
      mc->resample = &ncm_dataset_resample;
      ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_DISABLE);
      break;
    case NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX:
      mc->resample = &_ncm_fit_mc_resample_bstrap;
      ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_PARTIAL);
      break;
    case NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX:
      mc->resample = &_ncm_fit_mc_resample_bstrap;
      ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_TOTAL);
      break;
    case NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_fit_mc_set_use_threads:
 * @mc: a #NcmFitMC
 * @use_threads: whether to use OpenMP threads
 *
 * Sets whether @mc should run its Monte Carlo loop through the OpenMP
 * parallel code path. The real number of threads used is controlled by the
 * ambient OpenMP runtime state (the OMP_NUM_THREADS environment variable or
 * an explicit omp_set_num_threads() call), not by this property.
 *
 */
void
ncm_fit_mc_set_use_threads (NcmFitMC *mc, gboolean use_threads)
{
  mc->use_threads = use_threads;
}

/**
 * ncm_fit_mc_get_use_threads:
 * @mc: a #NcmFitMC
 *
 * Returns: whether @mc is set to use OpenMP threads.
 */
gboolean
ncm_fit_mc_get_use_threads (NcmFitMC *mc)
{
  return mc->use_threads;
}

/**
 * ncm_fit_mc_keep_order:
 * @mc: a #NcmFitMC
 * @keep_order: whether to keep the catalog in order of sampling
 *
 * When running with more than one thread, resample() calls race to claim
 * the next catalog row, so which physical realization lands in which row
 * is otherwise scheduling-dependent (not reproducible run to run, even for
 * the same seed). Setting @keep_order to %TRUE forces the resample step to
 * execute in strict loop-iteration order across all threads (only the
 * subsequent, expensive fit still runs in parallel), making the catalog's
 * row order deterministic and reproducible regardless of thread count.
 *
 */
void
ncm_fit_mc_keep_order (NcmFitMC *mc, gboolean keep_order)
{
  mc->keep_order = keep_order;
}

/**
 * ncm_fit_mc_set_fiducial:
 * @mc: a #NcmFitMC
 * @fiduc: (nullable): a #NcmMSet
 *
 * Makes @fiduc the fiducial model of @mc: the model the data are resampled from and
 * whose parameters start every fit. %NULL, or the model set of the fit itself, selects a
 * copy of the model set of the fit; any other @fiduc must have the same models and
 * parameters (ncm_mset_cmp()). After a run the model set of the fit holds the last best
 * fit.
 *
 */
void
ncm_fit_mc_set_fiducial (NcmFitMC *mc, NcmMSet *fiduc)
{
  NcmMSet *mset = ncm_fit_peek_mset (mc->fit);

  ncm_mset_clear (&mc->fiduc);

  if ((fiduc == NULL) || (fiduc == mset))
  {
    mc->fiduc = ncm_mset_dup (mset, mc->ser);
    ncm_serialize_reset (mc->ser, TRUE);
  }
  else
  {
    mc->fiduc = ncm_mset_ref (fiduc);
    g_assert (ncm_mset_cmp (mset, fiduc, FALSE));
  }
}

/**
 * ncm_fit_mc_set_rng:
 * @mc: a #NcmFitMC
 * @rng: a #NcmRNG
 *
 * Makes @rng the random number generator of the catalog of @mc, used by the
 * resampling, given to the catalog when ncm_fit_mc_start_run() builds it. Calling it
 * during a run aborts; a run started without one creates one with a random seed.
 *
 */
void
ncm_fit_mc_set_rng (NcmFitMC *mc, NcmRNG *rng)
{
  if (mc->started)
    g_error ("ncm_fit_mc_set_rng: Cannot change the RNG object during a run, call ncm_fit_mc_end_run() first.");

  if (mc->mcat != NULL)
  {
    ncm_mset_catalog_set_rng (mc->mcat, rng);
  }
  else
  {
    ncm_rng_clear (&mc->rng);
    mc->rng = ncm_rng_ref (rng);
  }
}

/**
 * ncm_fit_mc_is_running:
 * @mc: a #NcmFitMC
 *
 * Returns: whether ncm_fit_mc_start_run() was called without a later
 *   ncm_fit_mc_end_run()
 */
gboolean
ncm_fit_mc_is_running (NcmFitMC *mc)
{
  return mc->started;
}

static void _ncm_fit_mc_update_post (NcmFitMC *mc);

static void
_ncm_fit_mc_update_from_theta (NcmFitMC *mc, NcmVector *theta)
{
  ncm_mset_catalog_add_from_vector (mc->mcat, theta);

  _ncm_fit_mc_update_post (mc);
}

static void
_ncm_fit_mc_update_post (NcmFitMC *mc)
{
  const guint part = 5;
  const guint step = (mc->n / part) == 0 ? 1 : (mc->n / part);

  ncm_timer_task_increment (mc->nt);

  switch (mc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi          = ncm_timer_task_completed (mc->nt) % step;
      gboolean log_timeout = FALSE;

      if (ncm_timer_elapsed_since_last_log (mc->nt) > 60.0)
        log_timeout = TRUE;

      if (log_timeout || (stepi == 0) || ncm_timer_task_has_ended (mc->nt))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (mc->mcat);
        /* ncm_timer_task_accumulate (mc->nt, acc); */
        ncm_timer_task_log_elapsed (mc->nt);
        ncm_timer_task_log_mean_time (mc->nt);
        ncm_timer_task_log_time_left (mc->nt);
        ncm_timer_task_log_cur_datetime (mc->nt);
        ncm_timer_task_log_end_datetime (mc->nt);
      }

      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_mset_catalog_log_current_stats (mc->mcat);
      /* ncm_timer_task_increment (mc->nt); */
      ncm_timer_task_log_elapsed (mc->nt);
      ncm_timer_task_log_mean_time (mc->nt);
      ncm_timer_task_log_time_left (mc->nt);
      ncm_timer_task_log_cur_datetime (mc->nt);
      ncm_timer_task_log_end_datetime (mc->nt);
      break;
  }
}

static void
_ncm_fit_mc_resample_bstrap (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng)
{
  ncm_dataset_bootstrap_resample (dset, rng);
}

/**
 * ncm_fit_mc_start_run:
 * @mc: a #NcmFitMC
 * @error: a #GError
 *
 * Builds the catalog at the first call (ncm_fit_mc_set_storage()), failing with the
 * errors of #NcmMSetCatalogStorage when the file cannot be opened, its description is
 * incompatible or it holds an RNG of another algorithm or seed; a run also fails where
 * a recording cannot be opened on the storage.
 *
 * Starts a run: records the fiducial parameters, creates a random number generator if
 * the catalog has none, and skips the realizations the catalog already holds. It
 * computes no realization; ncm_fit_mc_run() or ncm_fit_mc_run_lre() do. Starting a
 * running @mc aborts.
 *
 */
void
ncm_fit_mc_start_run (NcmFitMC *mc, GError **error)
{
  NcmLikelihood *lh     = ncm_fit_peek_likelihood (mc->fit);
  NcmMSet *mset         = ncm_fit_peek_mset (mc->fit);
  const guint param_len = ncm_mset_total_len (mc->fiduc);
  NcmDataset *dset      = ncm_likelihood_peek_dataset (lh);
  gint mcat_cur_id;
  NcmRNG *mcat_rng;

  if (mc->started)
    g_error ("ncm_fit_mc_start_run: run already started, run ncm_fit_mc_end_run() first.");

  _ncm_fit_mc_build_catalog (mc, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mc_start_run: ");
  mcat_cur_id = ncm_mset_catalog_get_cur_id (mc->mcat);
  mcat_rng    = ncm_mset_catalog_peek_rng (mc->mcat);

  ncm_vector_clear (&mc->bf);
  mc->bf = ncm_vector_new (param_len);
  ncm_mset_param_get_vector (mc->fiduc, mc->bf);

  switch (mc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Starting Monte Carlo...\n");
      ncm_dataset_log_info (dset);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Fiducial model set:\n");
      ncm_mset_pretty_log (mc->fiduc);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Fitting model set:\n");
      ncm_mset_pretty_log (mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (mcat_rng == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);

    ncm_rng_set_random_seed (rng, FALSE);
    ncm_fit_mc_set_rng (mc, rng);

    if (mc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmFitMC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));

    ncm_rng_free (rng);
  }

  mc->started = TRUE;

  ncm_dataset_register_shared (dset, mc->ser);

  g_object_set (ncm_mset_catalog_peek_storage (mc->mcat),
                "smode", NCM_MSET_CATALOG_SYNC_TIMED,
                "sync-interval", NCM_FIT_MC_MIN_SYNC_INTERVAL,
                NULL);

  _ncm_fit_mc_start_recording (mc, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mc_start_run: ");

  if (mcat_cur_id > mc->cur_sample_id)
  {
    ncm_fit_mc_intern_skip (mc, mcat_cur_id - mc->cur_sample_id);
    g_assert_cmpint (mc->cur_sample_id, ==, mcat_cur_id);
  }
  else if (mcat_cur_id < mc->cur_sample_id)
  {
    g_error ("ncm_fit_mc_start_run: the catalog has fewer rows than the run [%d < %d].",
             mcat_cur_id, mc->cur_sample_id);
  }
}

/**
 * ncm_fit_mc_end_run:
 * @mc: a #NcmFitMC
 *
 * Ends the run: saves the catalog and disables the dataset bootstrap.
 *
 */
void
ncm_fit_mc_end_run (NcmFitMC *mc)
{
  NcmLikelihood *lh = ncm_fit_peek_likelihood (mc->fit);
  NcmDataset *dset  = ncm_likelihood_peek_dataset (lh);

  if (ncm_timer_task_is_running (mc->nt))
    ncm_timer_task_end (mc->nt);

  _ncm_fit_mc_end_recording (mc);
  ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_DISABLE);

  /* Releases the objects ncm_dataset_register_shared() kept for this run. */
  ncm_serialize_reset (mc->ser, FALSE);

  mc->started = FALSE;
}

/**
 * ncm_fit_mc_reset:
 * @mc: a #NcmFitMC
 *
 * Erases all realizations: the catalog is emptied and the next run starts from the
 * first realization.
 *
 */
void
ncm_fit_mc_reset (NcmFitMC *mc)
{
  mc->n             = 0;
  mc->cur_sample_id = -1;
  mc->write_index   = 0;
  mc->started       = FALSE;

  if (mc->mcat != NULL)
    ncm_mset_catalog_reset (mc->mcat);
}

static void
ncm_fit_mc_intern_skip (NcmFitMC *mc, guint n)
{
  if (n == 0)
    return;

  switch (mc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Skipping %u realizations, will start at %u-th realization.\n", n, mc->cur_sample_id + n + 1 + 1);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  mc->cur_sample_id += n;
  mc->write_index    = mc->cur_sample_id + 1;
}

static void _ncm_fit_mc_run_single (NcmFitMC *mc);
static void _ncm_fit_mc_mt_eval (glong i, glong f, gpointer data);

/**
 * ncm_fit_mc_run:
 * @mc: a #NcmFitMC
 * @n: total number of realizations
 *
 * Computes realizations until the run holds @n of them, counting those already done
 * and skipped. It requires a started run; with #NcmFitMC:use-threads the fits run in
 * OpenMP threads, each with its own copy of the fit.
 *
 */
void
ncm_fit_mc_run (NcmFitMC *mc, guint n)
{
  if (!mc->started)
    g_error ("ncm_fit_mc_run: run not started, run ncm_fit_mc_start_run() first.");

  if (n <= (guint) (mc->cur_sample_id + 1))
  {
    if (mc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Nothing to do, current Monte Carlo run is %d\n", mc->cur_sample_id + 1);
    }

    return;
  }

  mc->n = n - (mc->cur_sample_id + 1);

  switch (mc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      const GEnumValue *eval = ncm_cfg_enum_get_value (NCM_TYPE_FIT_MC_RESAMPLE_TYPE, mc->rtype);

      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Calculating [%06d] Monte Carlo fits [%s]\n", mc->n, eval->value_nick);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_timer_task_is_running (mc->nt))
  {
    ncm_timer_task_add_tasks (mc->nt, mc->n);
    ncm_timer_task_continue (mc->nt);
  }
  else
  {
    ncm_timer_task_start (mc->nt, mc->n);
    ncm_timer_set_name (mc->nt, "NcmFitMC");
  }

  if (mc->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (mc->nt);

  if (!mc->use_threads || (mc->n <= 1))
    _ncm_fit_mc_run_single (mc);
  else
    _ncm_fit_mc_mt_eval (0, mc->n, mc);

  ncm_timer_task_pause (mc->nt);
}

static void
_ncm_fit_mc_run_single (NcmFitMC *mc)
{
  NcmMSet *mset    = ncm_fit_peek_mset (mc->fit);
  NcmVector *theta = ncm_vector_new (mc->theta_len);
  guint i;

  for (i = 0; i < mc->n; i++)
  {
    ncm_mset_param_set_vector (mset, mc->bf);
    _ncm_fit_mc_resample (mc, mc->fit);
    ncm_fit_run (mc->fit, NCM_FIT_RUN_MSGS_NONE);

    {
      NcmFitState *fstate = ncm_fit_peek_state (mc->fit);

      ncm_vector_set (theta, 0, ncm_fit_state_get_m2lnL_curval (fstate));
      ncm_mset_fparams_get_vector_offset (mset, theta, mc->nadd_vals);
    }

    if (mc->func_oa != NULL)
    {
      glong k;

      for (k = 0; k < mc->func_oa->len; k++)
      {
        NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (mc->func_oa, k));
        const gdouble a_k = ncm_mset_func_eval0 (func, mset);

        ncm_vector_set (theta, k + 1, a_k);
      }
    }

    _ncm_fit_mc_update_from_theta (mc, theta);
    mc->write_index++;
  }

  ncm_vector_free (theta);
}

static void
_ncm_fit_mc_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitMC *mc                 = NCM_FIT_MC (data);
  const glong n                = f - i;
  const glong cur_pos          = mc->cur_sample_id + 1;
  const glong last_write_index = mc->write_index + n;
  GPtrArray *thetas            = g_ptr_array_sized_new (n);
  GArray *ready                = g_array_sized_new (FALSE, FALSE, sizeof (gboolean), n);
  glong li;

  g_ptr_array_set_free_func (thetas, (GDestroyNotify) ncm_vector_free);

  for (li = 0; li < n; li++)
  {
    NcmVector *theta = ncm_vector_new (mc->theta_len);

    g_ptr_array_add (thetas, theta);
    g_array_index (ready, gboolean, li) = FALSE;
  }

#pragma omp parallel shared (mc, thetas, ready)
  {
    /* Each thread gets a fit */
    NcmFit *fit              = NULL;
    NcmMSet *mset            = NULL;
    NcmFitState *fstate      = NULL;
    NcmVector *theta         = NULL;
    NcmObjArray *funcs_array = NULL;
    NcmVector *oa_vals       = NULL;
    glong j, k, sample_index;

    #pragma omp critical(dup_fit)
    {
      fit    = ncm_fit_dup (mc->fit, mc->ser);
      mset   = ncm_fit_peek_mset (fit);
      fstate = ncm_fit_peek_state (fit);

      if ((mc->func_oa != NULL) && (mc->func_oa->len > 0))
      {
        funcs_array = ncm_serialize_dup_array (mc->ser, mc->func_oa);
        oa_vals     = ncm_vector_new (mc->func_oa->len);
      }

      ncm_serialize_reset (mc->ser, TRUE);
    }

    #pragma omp for ordered

    for (j = 0; j < n; j++)
    {
      ncm_mset_param_set_vector (mset, mc->bf);

      /* keep_order is fixed for the run, so every iteration takes the same branch, as
       * "omp ordered" requires. With it the resamplings (and the sample indices) follow
       * the loop order; without it "critical" only serializes the shared RNG and
       * counter, and the indices follow the thread arrival order. */
      if (mc->keep_order)
      {
        #pragma omp ordered
        {
          sample_index = _ncm_fit_mc_resample (mc, fit);
        }
      }
      else
      {
        #pragma omp critical(resample_phase)
        {
          sample_index = _ncm_fit_mc_resample (mc, fit);
        }
      }

      theta = g_ptr_array_index (thetas, sample_index - cur_pos);

      /* This is the slow section that we want to parallelize over */
      ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

      if ((mc->func_oa != NULL) && (mc->func_oa->len > 0))
      {
        for (k = 0; k < mc->func_oa->len; k++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (funcs_array, k));
          const gdouble a_k = ncm_mset_func_eval0 (func, mset);

          ncm_vector_set (oa_vals, k, a_k);
        }
      }

      /* End of the slow section */

      #pragma omp critical(update_phase)
      {
        ncm_vector_set (theta, 0, ncm_fit_state_get_m2lnL_curval (fstate));
        ncm_mset_fparams_get_vector_offset (mset, theta, mc->nadd_vals);

        if ((mc->func_oa != NULL) && (mc->func_oa->len > 0))
        {
          for (k = 0; k < mc->func_oa->len; k++)
          {
            const gdouble a_k = ncm_vector_get (oa_vals, k);

            ncm_vector_set (theta, k + 1, a_k);
          }
        }

        g_array_index (ready, gboolean, sample_index - cur_pos) = TRUE;

        {
          const glong cur_write_index = mc->write_index;

          for (k = cur_write_index; k < last_write_index; k++)
          {
            if (g_array_index (ready, gboolean, k - cur_pos))
            {
              NcmVector *saving_theta = g_ptr_array_index (thetas, k - cur_pos);

              _ncm_fit_mc_update_from_theta (mc, saving_theta);
              mc->write_index++;
            }
            else
            {
              break;
            }
          }
        }
      }
    }

    ncm_fit_clear (&fit);
    ncm_obj_array_clear (&funcs_array);
  }

  {
    const glong cur_write_index = mc->write_index;

    for (li = cur_write_index; li < last_write_index; li++)
    {
      NcmVector *theta = g_ptr_array_index (thetas, li - cur_pos);

      g_assert (g_array_index (ready, gboolean, li - cur_pos));

      _ncm_fit_mc_update_from_theta (mc, theta);
      mc->write_index++;
    }
  }

  g_ptr_array_free (thetas, TRUE);
  g_array_free (ready, TRUE);
}

/**
 * ncm_fit_mc_run_lre:
 * @mc: a #NcmFitMC
 * @prerun: number of pre-runs
 * @lre: largest relative error
 *
 * Computes at least @prerun realizations (100 when zero), then more until the largest
 * relative error of the parameter means (ncm_mset_catalog_largest_error()) is below
 * @lre.
 *
 */
void
ncm_fit_mc_run_lre (NcmFitMC *mc, guint prerun, gdouble lre)
{
  gdouble lerror;
  const gdouble lre2 = lre * lre;

  g_assert_cmpfloat (lre, >, 0.0);
  /* g_assert_cmpfloat (lre, <, 1.0); */

  prerun = (prerun == 0) ? 100 : prerun;

  if (ncm_mset_catalog_len (mc->mcat) < prerun)
  {
    guint prerun_left = prerun - ncm_mset_catalog_len (mc->mcat);

    if (mc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitMC: Running first %u pre-runs...\n", prerun_left);

    ncm_fit_mc_run (mc, prerun);
  }

  lerror = ncm_mset_catalog_largest_error (mc->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n             = ncm_mset_catalog_len (mc->mcat);
    gdouble m             = n * lerror2 / lre2;
    guint runs            = ((m - n) > 1000.0) ? ceil ((m - n) * 0.25) : ceil (m - n);

    if (mc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmFitMC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmFitMC: Running more %u runs...\n", runs);
    }

    ncm_fit_mc_run (mc, mc->cur_sample_id + runs + 1);
    lerror = ncm_mset_catalog_largest_error (mc->mcat);
  }

  if (mc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    g_message ("# NcmFitMC: Largest relative error %e attained: %e\n", lre, lerror);
}

/**
 * ncm_fit_mc_mean_covar:
 * @mc: a #NcmFitMC
 *
 * Stores the mean and covariance of the best fits in the catalog as the parameters and
 * covariance of the #NcmFitState of the fit, and sets the model set of the catalog to
 * the mean.
 *
 */
void
ncm_fit_mc_mean_covar (NcmFitMC *mc)
{
  NcmMSet *mset       = ncm_mset_catalog_peek_mset (mc->mcat);
  NcmFitState *fstate = ncm_fit_peek_state (mc->fit);
  NcmVector *fparams  = ncm_fit_state_peek_fparams (fstate);
  NcmMatrix *covar    = ncm_fit_state_peek_covar (fstate);

  ncm_mset_catalog_get_mean (mc->mcat, &fparams);
  ncm_mset_catalog_get_covar (mc->mcat, &covar);
  ncm_mset_fparams_set_vector (mset, fparams);

  ncm_fit_state_set_has_covar (fstate, TRUE);
}

/**
 * ncm_fit_mc_get_catalog:
 * @mc: a #NcmFitMC
 *
 * Returns: (transfer full): the catalog of @mc
 */
NcmMSetCatalog *
ncm_fit_mc_get_catalog (NcmFitMC *mc)
{
  return ncm_mset_catalog_ref (ncm_fit_mc_peek_catalog (mc));
}

/**
 * ncm_fit_mc_peek_catalog:
 * @mc: a #NcmFitMC
 *
 * The catalog is built by the first ncm_fit_mc_start_run(); before it there is none and
 * this aborts.
 *
 * Returns: (transfer none): the catalog of @mc
 */
NcmMSetCatalog *
ncm_fit_mc_peek_catalog (NcmFitMC *mc)
{
  if (mc->mcat == NULL)
    g_error ("ncm_fit_mc_peek_catalog: the catalog is built by ncm_fit_mc_start_run(), call it first.");

  return mc->mcat;
}


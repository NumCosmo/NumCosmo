/***************************************************************************
 *            ncm_fit_mcbs.c
 *
 *  Tue February 11 13:54:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>>
 ****************************************************************************/
/*
 * ncm_fit_mcbs.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmFitMCBS:
 *
 * Monte Carlo study of the bootstrap mean.
 *
 * For each data realization resampled from a fiducial model it runs a bootstrap of that
 * realization (two #NcmFitMC) and stores the mean of the bootstrap best fits as a row
 * of a #NcmMSetCatalog. The catalog shows how the bootstrap mean is distributed over
 * data realizations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_fit_mcbs.h"
#include "ncm/fit/ncm_mset_catalog_storage_fits.h"
#include "ncm/fit/ncm_mset_catalog_storage_memory.h"
#include "ncm_enum_types.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"

#include <gio/gio.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_STORAGE,
};

struct _NcmFitMCBS
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmFitMC *mc_resample;
  NcmFitMC *mc_bstrap;
  NcmMSetCatalog *mcat;
  NcmMSetCatalogStorage *storage;
  gchar *base_name;
  NcmRNG *rng;
};


G_DEFINE_TYPE (NcmFitMCBS, ncm_fit_mcbs, G_TYPE_OBJECT)

static void
ncm_fit_mcbs_init (NcmFitMCBS *mcbs)
{
  mcbs->fit         = NULL;
  mcbs->mc_resample = NULL;
  mcbs->mc_bstrap   = NULL;
  mcbs->mcat        = NULL;
  mcbs->storage     = NULL;
  mcbs->base_name   = NULL;
  mcbs->rng         = NULL;
}

static void
ncm_fit_mcbs_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMCBS *mcbs = NCM_FIT_MCBS (object);

  g_return_if_fail (NCM_IS_FIT_MCBS (object));
  {
    switch (prop_id)
    {
      case PROP_FIT:
        mcbs->fit         = g_value_dup_object (value);
        mcbs->mc_resample = ncm_fit_mc_new (mcbs->fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
        break;
      case PROP_STORAGE:
        ncm_fit_mcbs_set_storage (mcbs, g_value_get_object (value));
        break;
      default:                                                      /* LCOV_EXCL_LINE */
        G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
        break;                                                      /* LCOV_EXCL_LINE */
    }
  }
}

static void
ncm_fit_mcbs_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitMCBS *mcbs = NCM_FIT_MCBS (object);

  g_return_if_fail (NCM_IS_FIT_MCBS (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, mcbs->fit);
      break;
    case PROP_STORAGE:
      g_value_set_object (value, mcbs->storage);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fit_mcbs_dispose (GObject *object)
{
  NcmFitMCBS *mcbs = NCM_FIT_MCBS (object);

  ncm_fit_clear (&mcbs->fit);

  ncm_fit_mc_clear (&mcbs->mc_resample);
  ncm_fit_mc_clear (&mcbs->mc_bstrap);

  ncm_mset_catalog_clear (&mcbs->mcat);
  ncm_mset_catalog_storage_clear (&mcbs->storage);
  ncm_rng_clear (&mcbs->rng);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mcbs_parent_class)->dispose (object);
}

static void
_ncm_fit_mcbs_finalize (GObject *object)
{
  NcmFitMCBS *mcbs = NCM_FIT_MCBS (object);

  g_clear_pointer (&mcbs->base_name, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mcbs_parent_class)->finalize (object);
}

static void
ncm_fit_mcbs_class_init (NcmFitMCBSClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_fit_mcbs_set_property;
  object_class->get_property = &ncm_fit_mcbs_get_property;
  object_class->dispose      = &_ncm_fit_mcbs_dispose;
  object_class->finalize     = &_ncm_fit_mcbs_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitMCBS:storage:
   *
   * The #NcmMSetCatalogStorage of the catalog that ncm_fit_mcbs_run() builds, see
   * ncm_fit_mcbs_set_storage(); without one the catalog is kept in a
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
}

/**
 * ncm_fit_mcbs_new:
 * @fit: a #NcmFit
 *
 * Creates a #NcmFitMCBS for @fit.
 *
 * Returns: (transfer full): a new #NcmFitMCBS.
 */
NcmFitMCBS *
ncm_fit_mcbs_new (NcmFit *fit)
{
  return g_object_new (NCM_TYPE_FIT_MCBS,
                       "fit", fit,
                       NULL);
}

/**
 * ncm_fit_mcbs_free:
 * @mcbs: a #NcmFitMCBS
 *
 * Decreases the reference count of @mcbs.
 *
 */
void
ncm_fit_mcbs_free (NcmFitMCBS *mcbs)
{
  g_object_unref (mcbs);
}

/**
 * ncm_fit_mcbs_clear:
 * @mcbs: a #NcmFitMCBS
 *
 * Decreases the reference count of *@mcbs and sets it to %NULL.
 *
 */
void
ncm_fit_mcbs_clear (NcmFitMCBS **mcbs)
{
  g_clear_object (mcbs);
}

/**
 * ncm_fit_mcbs_set_storage:
 * @mcbs: a #NcmFitMCBS
 * @storage: a #NcmMSetCatalogStorage
 *
 * Makes @storage the storage of the catalog that ncm_fit_mcbs_run() builds; the rows it
 * holds are erased by the run. With a #NcmMSetCatalogStorageFits each bootstrap catalog
 * is saved next to its file, as the file name without its .fits extension followed by
 * -bstrap-NNNNNN.fits, NNNNNN the realization index; with a memory storage the bootstrap
 * catalogs are kept in memory. Once the catalog is built its storage cannot change, so
 * calling it after the first run aborts.
 *
 */
void
ncm_fit_mcbs_set_storage (NcmFitMCBS *mcbs, NcmMSetCatalogStorage *storage)
{
  g_assert_nonnull (storage);

  if (mcbs->mcat != NULL)
    g_error ("ncm_fit_mcbs_set_storage: the catalog is built by ncm_fit_mcbs_run() and its storage cannot change.");

  ncm_mset_catalog_storage_clear (&mcbs->storage);
  mcbs->storage = ncm_mset_catalog_storage_ref (storage);
}

/*
 * Builds the catalog of the bootstrap means at the first run: a description of the free
 * parameters and the -2ln(L) column, on the storage given (its rows erased) or a memory
 * storage, with the RNG set before; the base name of the bootstrap files is taken from a
 * FITS storage's file name. The resampling NcmFitMC is seeded as the catalog's RNG here,
 * once: its catalog refuses an RNG once it has rows.
 */
static void
_ncm_fit_mcbs_build_catalog (NcmFitMCBS *mcbs, GError **error)
{
  NcmMSet *mset          = ncm_fit_peek_mset (mcbs->fit);
  const gchar *names[]   = {NCM_MSET_CATALOG_M2LNL_COLNAME, NULL};
  const gchar *symbols[] = {NCM_MSET_CATALOG_M2LNL_SYMBOL, NULL};
  NcmMSetCatalogDesc *desc;

  if (mcbs->mcat != NULL)
    return;

  desc = ncm_mset_catalog_desc_new (mset, 1, FALSE, names, symbols, 0, NCM_MSET_CATALOG_RTYPE_MCBS, error);
  NCM_UTIL_ON_ERROR_RETURN (error, , );

  if (mcbs->storage == NULL)
    mcbs->storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_memory_new ());

  ncm_mset_catalog_storage_set_desc (mcbs->storage, desc, error);
  NCM_UTIL_ON_ERROR_RETURN (error, ncm_mset_catalog_desc_free (desc), );

  mcbs->mcat = ncm_mset_catalog_new_from_storage (mcbs->storage);

  g_clear_pointer (&mcbs->base_name, g_free);

  if (NCM_IS_MSET_CATALOG_STORAGE_FITS (mcbs->storage))
    mcbs->base_name = ncm_util_basename_fits (ncm_mset_catalog_storage_fits_peek_filename (NCM_MSET_CATALOG_STORAGE_FITS (mcbs->storage)));

  if (!ncm_mset_catalog_is_empty (mcbs->mcat))
    ncm_mset_catalog_reset (mcbs->mcat);

  ncm_mset_catalog_set_sampler (mcbs->mcat, "from-model");

  if (mcbs->rng != NULL)
  {
    ncm_mset_catalog_set_rng (mcbs->mcat, mcbs->rng);
    ncm_rng_clear (&mcbs->rng);
  }

  {
    NcmRNG *mcat_rng = ncm_mset_catalog_peek_rng (mcbs->mcat);

    if (mcat_rng != NULL)
    {
      NcmRNG *rng = ncm_rng_seeded_new (ncm_rng_get_algo (mcat_rng), ncm_rng_get_seed (mcat_rng));

      ncm_fit_mc_set_rng (mcbs->mc_resample, rng);
      ncm_rng_free (rng);
    }
  }

  ncm_mset_catalog_desc_free (desc);
}

/**
 * ncm_fit_mcbs_set_rng:
 * @mcbs: a #NcmFitMCBS
 * @rng: a #NcmRNG
 *
 * Makes @rng the random number generator of the catalog of @mcbs; the realizations are
 * resampled with a generator seeded as @rng. Calling it during a run aborts.
 *
 */
void
ncm_fit_mcbs_set_rng (NcmFitMCBS *mcbs, NcmRNG *rng)
{
  if ((mcbs->mc_bstrap != NULL) && ncm_fit_mc_is_running (mcbs->mc_bstrap))
    g_error ("ncm_fit_mcbs_set_rng: Cannot change the RNG object during a run.");

  if (mcbs->mcat != NULL)
  {
    ncm_mset_catalog_set_rng (mcbs->mcat, rng);
  }
  else
  {
    ncm_rng_clear (&mcbs->rng);
    mcbs->rng = ncm_rng_ref (rng);
  }
}

/**
 * ncm_fit_mcbs_run:
 * @mcbs: a #NcmFitMCBS
 * @fiduc: (nullable): the fiducial #NcmMSet, see ncm_fit_mc_set_fiducial()
 * @nf: number of realizations
 * @nbstraps: number of bootstrap fits per realization
 * @rtype: bootstrap #NcmFitMCResampleType
 * @mtype: #NcmFitRunMsgs of the bootstrap runs
 * @bsmt: whether the bootstrap fits use OpenMP threads
 * @error: a #GError
 *
 * Builds the catalog at the first call (ncm_fit_mcbs_set_storage()), failing with the
 * errors of #NcmMSetCatalogStorage when the file cannot be opened or its description
 * is incompatible; a run also fails where a recording cannot be opened on a storage.
 *
 * Resamples @nf realizations from @fiduc and, for each, runs @nbstraps
 * bootstrap fits of type @rtype, adding their mean to the catalog of @mcbs. At the end
 * the mean and covariance of the catalog become the parameters and covariance of the
 * #NcmFitState of the fit. @rtype must be a bootstrap; #NCM_FIT_MC_RESAMPLE_FROM_MODEL
 * aborts.
 *
 * Multi-threaded bootstraps (@bsmt = %TRUE) may give incorrect results; use %FALSE.
 *
 */
void
ncm_fit_mcbs_run (NcmFitMCBS *mcbs, NcmMSet *fiduc, guint nf, guint nbstraps, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype, gboolean bsmt, GError **error)
{
  NcmFitState *fstate    = ncm_fit_peek_state (mcbs->fit);
  NcmLikelihood *lh      = ncm_fit_peek_likelihood (mcbs->fit);
  NcmDataset *dset       = ncm_likelihood_peek_dataset (lh);
  const GEnumValue *eval = ncm_cfg_enum_get_value (NCM_TYPE_FIT_MC_RESAMPLE_TYPE, rtype);
  NcmMSetCatalogStorage *storage;
  NcmDatasetBStrapType bstrap_type;
  guint i;

  if (rtype == NCM_FIT_MC_RESAMPLE_FROM_MODEL)
    g_error ("ncm_fit_mcbs_run: the internal run must be a bootstrap: NCM_FIT_MC_RESAMPLE_BOOTSTRAP_*.");

  bstrap_type = (rtype == NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX) ? NCM_DATASET_BSTRAP_PARTIAL : NCM_DATASET_BSTRAP_TOTAL;

  _ncm_fit_mcbs_build_catalog (mcbs, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mcbs_run: ");
  storage = ncm_mset_catalog_peek_storage (mcbs->mcat);

  ncm_fit_mc_set_mtype (mcbs->mc_resample, NCM_FIT_RUN_MSGS_SIMPLE);
  ncm_fit_mc_set_fiducial (mcbs->mc_resample, fiduc);

  ncm_fit_mc_start_run (mcbs->mc_resample, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mcbs_run: ");

  /* Without an RNG the first resampling run drew one; the catalog takes it. */
  if (ncm_mset_catalog_peek_rng (mcbs->mcat) == NULL)
    ncm_mset_catalog_set_rng (mcbs->mcat,
                              ncm_mset_catalog_peek_rng (ncm_fit_mc_peek_catalog (mcbs->mc_resample)));

  ncm_mset_catalog_storage_start_recording (storage, "from-model", eval->value_nick, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mcbs_run: ");

  /* One NcmFitMC per realization: its catalog is built by its start_run on the
   * realization's own file. */
  for (i = 0; i < nf; i++)
  {
    ncm_dataset_bootstrap_set (dset, NCM_DATASET_BSTRAP_DISABLE);
    ncm_fit_mc_run (mcbs->mc_resample, i + 1);
    ncm_dataset_bootstrap_set (dset, bstrap_type);

    ncm_fit_mc_clear (&mcbs->mc_bstrap);
    mcbs->mc_bstrap = ncm_fit_mc_new (mcbs->fit, rtype, mtype);
    ncm_fit_mc_set_use_threads (mcbs->mc_bstrap, bsmt);

    if (mcbs->base_name != NULL)
    {
      gchar *bstrap_str                       = g_strdup_printf ("%s-bstrap-%06d.fits", mcbs->base_name, i);
      NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (bstrap_str, FALSE, error);

      NCM_UTIL_ON_ERROR_FORWARD (error, g_free (bstrap_str), , "ncm_fit_mcbs_run: ");
      ncm_fit_mc_set_storage (mcbs->mc_bstrap, NCM_MSET_CATALOG_STORAGE (storage_fits));
      ncm_mset_catalog_storage_fits_free (storage_fits);
      g_free (bstrap_str);
    }

    ncm_fit_mc_start_run (mcbs->mc_bstrap, error);
    NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mcbs_run: ");
    ncm_fit_mc_run (mcbs->mc_bstrap, nbstraps);

    ncm_mset_catalog_add_from_vector (mcbs->mcat,
                                      ncm_stats_vec_peek_mean (ncm_mset_catalog_peek_pstats (ncm_fit_mc_peek_catalog (mcbs->mc_bstrap))));
    ncm_mset_catalog_log_current_stats (mcbs->mcat);

    ncm_fit_mc_end_run (mcbs->mc_bstrap);
  }

  ncm_mset_catalog_storage_end_recording (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_fit_mcbs_run: ");

  {
    NcmVector *fparams = ncm_fit_state_peek_fparams (fstate);
    NcmMatrix *covar   = ncm_fit_state_peek_covar (fstate);

    ncm_mset_catalog_get_mean (mcbs->mcat, &fparams);
    ncm_mset_catalog_get_covar (mcbs->mcat, &covar);

    ncm_fit_state_set_has_covar (fstate, TRUE);
  }
  ncm_fit_mc_end_run (mcbs->mc_resample);
}

/**
 * ncm_fit_mcbs_get_catalog:
 * @mcbs: a #NcmFitMCBS
 *
 * The catalog is built by the first ncm_fit_mcbs_run(); before it there is none and
 * this aborts.
 *
 * Returns: (transfer full): the catalog of @mcbs, one bootstrap mean per realization
 */
NcmMSetCatalog *
ncm_fit_mcbs_get_catalog (NcmFitMCBS *mcbs)
{
  if (mcbs->mcat == NULL)
    g_error ("ncm_fit_mcbs_get_catalog: the catalog is built by ncm_fit_mcbs_run(), call it first.");

  return ncm_mset_catalog_ref (mcbs->mcat);
}


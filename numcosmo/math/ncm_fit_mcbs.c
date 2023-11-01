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
 * SECTION:ncm_fit_mcbs
 * @title: NcmFitMCBS
 * @short_description: Monte Carlo and bootstrap analysis.
 * @include: numcosmo/math/ncm_fit_mcbs.h
 * @stability: Unstable
 *
 * This class implements the Monte Carlo and bootstrap analysis.
 * It performs a Monte Carlo analysis where for each sample of the
 * likelihood function, a bootstrap analysis is performed.
 * The results are stored in a #NcmMSetCatalog.
 *
 * The main objective of this class is to provide a way to estimate
 * the accuracy of the bootstrap analysis.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_mcbs.h"
#include "math/ncm_cfg.h"

#include <gio/gio.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_FILE,
};

struct _NcmFitMCBS
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmFitMC *mc_resample;
  NcmFitMC *mc_bstrap;
  NcmMSetCatalog *mcat;
  gchar *base_name;
};


G_DEFINE_TYPE (NcmFitMCBS, ncm_fit_mcbs, G_TYPE_OBJECT);

static void
ncm_fit_mcbs_init (NcmFitMCBS *mcbs)
{
  mcbs->fit         = NULL;
  mcbs->mc_resample = NULL;
  mcbs->mc_bstrap   = NULL;
  mcbs->mcat        = NULL;
  mcbs->base_name   = NULL;
}

static void
ncm_fit_mcbs_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMCBS *mcbs = NCM_FIT_MCBS (object);

  g_return_if_fail (NCM_IS_FIT_MCBS (object));
  {
    NcmMSet *mset = ncm_fit_peek_mset (mcbs->fit);

    switch (prop_id)
    {
      case PROP_FIT:
        mcbs->fit         = g_value_dup_object (value);
        mcbs->mc_resample = ncm_fit_mc_new (mcbs->fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
        mcbs->mc_bstrap   = ncm_fit_mc_new (mcbs->fit, NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, NCM_FIT_RUN_MSGS_NONE);
        mcbs->mcat        = ncm_mset_catalog_new (mset, 1, 1, FALSE,
                                                  NCM_MSET_CATALOG_M2LNL_COLNAME, NCM_MSET_CATALOG_M2LNL_SYMBOL,
                                                  NULL);
        ncm_mset_catalog_set_m2lnp_var (mcbs->mcat, 0);
        ncm_mset_catalog_set_run_type (mcbs->mcat, NCM_MSET_CATALOG_RTYPE_BSTRAP_MEAN);
        break;
      case PROP_FILE:
        ncm_fit_mcbs_set_filename (mcbs, g_value_get_string (value));
        break;
      default:
        G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
        break;
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
    case PROP_FILE:
      g_value_set_string (value, ncm_mset_catalog_peek_filename (mcbs->mcat));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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

  g_object_class_install_property (object_class,
                                   PROP_FILE,
                                   g_param_spec_string ("filename",
                                                        NULL,
                                                        "Data filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_fit_mcbs_new:
 * @fit: a #NcmFit
 *
 * Creates a new #NcmFitMCBS object.
 *
 * Returns: (transfer full): a new #NcmFitMCBS object.
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
 * Decreases the reference count of @mcbs by one.
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
 * If *@mcbs is not %NULL, decreases its reference count by one and
 * sets *@mcbs to %NULL.
 *
 */
void
ncm_fit_mcbs_clear (NcmFitMCBS **mcbs)
{
  g_clear_object (mcbs);
}

/**
 * ncm_fit_mcbs_set_filename:
 * @mcbs: a #NcmFitMCBS
 * @filename: a filename
 *
 * Sets the filename of the data to be used in the analysis.
 *
 */
void
ncm_fit_mcbs_set_filename (NcmFitMCBS *mcbs, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (mcbs->mcat);

  if ((cur_filename != NULL) && (strcmp (cur_filename, filename) == 0))
  {
    return;
  }
  else
  {
    GError *error          = NULL;
    GMatchInfo *match_info = NULL;
    GRegex *fits_ext       = g_regex_new ("(.*)\\.[fF][iI][tT][sS]$", 0, 0, &error);

    g_clear_pointer (&mcbs->base_name, g_free);

    if (g_regex_match (fits_ext, filename, 0, &match_info))
    {
      mcbs->base_name = g_match_info_fetch (match_info, 1);
      g_match_info_free (match_info);
    }
    else
    {
      mcbs->base_name = g_strdup (filename);
    }

    g_regex_unref (fits_ext);
    {
      gchar *resample_str = g_strdup_printf ("%s-resample.fits", mcbs->base_name);

      ncm_mset_catalog_set_file (mcbs->mcat, filename);

      ncm_mset_catalog_reset (mcbs->mcat);
      ncm_mset_catalog_erase_data (mcbs->mcat);

      g_free (resample_str);
    }
  }
}

/**
 * ncm_fit_mcbs_set_rng:
 * @mcbs: a #NcmFitMCBS
 * @rng: a #NcmRNG
 *
 * Sets the random number generator to be used in the analysis.
 *
 */
void
ncm_fit_mcbs_set_rng (NcmFitMCBS *mcbs, NcmRNG *rng)
{
  if (ncm_fit_mc_is_running (mcbs->mc_bstrap))
    g_error ("ncm_fit_mcbs_set_rng: Cannot change the RNG object during a run.");

  ncm_mset_catalog_set_rng (mcbs->mcat, rng);
}

/**
 * ncm_fit_mcbs_run:
 * @mcbs: a #NcmFitMCBS
 * @fiduc: a #NcmMSet
 * @ni: index of the first sample to be used in the analysis
 * @nf: index of the last sample to be used in the analysis
 * @nbstraps: number of bootstrap samples
 * @rtype: a #NcmFitMCResampleType
 * @mtype: a #NcmFitRunMsgs
 * @bsmt: number of threads to be used in the bootstrap analysis
 *
 * Runs the Monte Carlo and bootstrap analysis.
 * The results are stored in the catalog of @mcbs.
 * The catalog is cleared before the analysis.
 *
 * WARNING not working correctly with bsmt > 0 FIXME
 *
 */
void
ncm_fit_mcbs_run (NcmFitMCBS *mcbs, NcmMSet *fiduc, guint ni, guint nf, guint nbstraps, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype, guint bsmt)
{
  NcmRNG *mcat_rng     = ncm_mset_catalog_peek_rng (mcbs->mcat);
  NcmFitState *fstate  = ncm_fit_peek_state (mcbs->fit);
  NcmLikelihood *lh    = ncm_fit_peek_likelihood (mcbs->fit);
  gboolean cat_has_rng = FALSE;
  guint i;

  ncm_fit_mc_set_rtype (mcbs->mc_resample, NCM_FIT_MC_RESAMPLE_FROM_MODEL);
  ncm_fit_mc_set_mtype (mcbs->mc_resample, NCM_FIT_RUN_MSGS_SIMPLE);
  ncm_fit_mc_set_fiducial (mcbs->mc_resample, fiduc);

  if (mcat_rng != NULL)
  {
    NcmRNG *rng = ncm_rng_seeded_new (ncm_rng_get_algo (mcat_rng), ncm_rng_get_seed (mcat_rng));

    ncm_fit_mc_set_rng (mcbs->mc_resample, rng);
    ncm_rng_free (rng);
    cat_has_rng = TRUE;
  }

  ncm_fit_mc_start_run (mcbs->mc_resample);

  if (!cat_has_rng)
    ncm_mset_catalog_set_rng (mcbs->mcat,
                              ncm_mset_catalog_peek_rng (ncm_fit_mc_peek_catalog (mcbs->mc_resample)));

  if (ni > 0)
    ncm_fit_mc_set_first_sample_id (mcbs->mc_resample, ni);

  ncm_fit_mc_set_rtype (mcbs->mc_bstrap, rtype);
  ncm_fit_mc_set_mtype (mcbs->mc_bstrap, mtype);
  ncm_fit_mc_set_nthreads (mcbs->mc_bstrap, bsmt);

  if (rtype == NCM_FIT_MC_RESAMPLE_FROM_MODEL)
    g_error ("ncm_fit_mcbs_run: the internal run must be a bootstrap: NCM_FIT_MC_RESAMPLE_BOOTSTRAP_*.");

  for (i = ni; i < nf; i++)
  {
    ncm_dataset_bootstrap_set (lh->dset, NCM_DATASET_BSTRAP_DISABLE);
    /*ncm_fit_mc_set_first_sample_id (mcbs->mc_resample, i + 1); */
    ncm_fit_mc_run (mcbs->mc_resample, i + 1);
    ncm_dataset_bootstrap_set (lh->dset, NCM_DATASET_BSTRAP_TOTAL); /* FIXME */

    if (mcbs->base_name != NULL)
    {
      gchar *bstrap_str = g_strdup_printf ("%s-bstrap-%06d.fits", mcbs->base_name, i);

      ncm_fit_mc_set_data_file (mcbs->mc_bstrap, bstrap_str);
      g_free (bstrap_str);
    }

    ncm_fit_mc_start_run (mcbs->mc_bstrap);
    ncm_fit_mc_run (mcbs->mc_bstrap, nbstraps);

    ncm_mset_catalog_add_from_vector (mcbs->mcat,
                                      ncm_mset_catalog_peek_pstats (ncm_fit_mc_peek_catalog (mcbs->mc_bstrap))->mean);
    ncm_mset_catalog_log_current_stats (mcbs->mcat);

    ncm_fit_mc_end_run (mcbs->mc_bstrap);
    ncm_fit_mc_reset (mcbs->mc_bstrap);
  }

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
 * Gets the generated catalog of @mcbs.
 *
 * Returns: (transfer full): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_mcbs_get_catalog (NcmFitMCBS *mcbs)
{
  return ncm_mset_catalog_ref (mcbs->mcat);
}


/***************************************************************************
 *            ncm_fit_mcbs.c
 *
 *  Tue February 11 13:54:57 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>>
 ****************************************************************************/
/*
 * ncm_fit_mcbs.c
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

/**
 * SECTION:ncm_fit_mcbs
 * @title: Monte Carlo Bootstrap Analysis
 * @short_description: Object implementing Monte Carlo of Bootstrap analysis
 *
 * FIXME
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

G_DEFINE_TYPE (NcmFitMCBS, ncm_fit_mcbs, G_TYPE_OBJECT);

static void
ncm_fit_mcbs_init (NcmFitMCBS *mcbs)
{
  mcbs->fit = NULL;
  mcbs->mc_resample = NULL;
  mcbs->mc_bstrap = NULL;
  mcbs->fcat = NULL;
  mcbs->rng = NULL;
  mcbs->base_name = NULL;
}

static void
ncm_fit_mcbs_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMCBS *mcbs = NCM_FIT_MCBS (object);
  g_return_if_fail (NCM_IS_FIT_MCBS (object));

  switch (prop_id)
  {
    case PROP_FIT:
      mcbs->fit = g_value_dup_object (value);
      mcbs->mc_resample = ncm_fit_mc_new (mcbs->fit, NCM_FIT_MC_RESAMPLE_FROM_MODEL, NCM_FIT_RUN_MSGS_NONE);
      mcbs->mc_bstrap = ncm_fit_mc_new (mcbs->fit, NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX, NCM_FIT_RUN_MSGS_NONE);
      mcbs->fcat = ncm_fit_catalog_new (mcbs->fit);
      ncm_fit_catalog_set_run_type (mcbs->fcat, NCM_FIT_CATALOG_RTYPE_BSTRAP_MEAN);
      break;
    case PROP_FILE:
      ncm_fit_mcbs_set_filename (mcbs, g_value_get_string (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
      g_value_set_string (value, mcbs->fcat->file);
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

  ncm_fit_catalog_clear (&mcbs->fcat);

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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
                                   PROP_FIT,
                                   g_param_spec_string ("filename",
                                                        NULL,
                                                        "Data filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_fit_mcbs_new:
 * @fit: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
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
 * @mcbs: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcbs_free (NcmFitMCBS *mcbs)
{
  g_object_unref (mcbs);
}

/**
 * ncm_fit_mcbs_clear:
 * @mcbs: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcbs_clear (NcmFitMCBS **mcbs)
{
  g_clear_object (mcbs);
}

/**
 * ncm_fit_mcbs_set_filename:
 * @mcbs: FIXME
 * @filename: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcbs_set_filename (NcmFitMCBS *mcbs, const gchar *filename)
{
  const gchar *cur_filename = ncm_fit_catalog_peek_filename (mcbs->fcat); 

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;
  else
  {
    GError *error = NULL;
    GMatchInfo *match_info = NULL;
    GRegex *fits_ext = g_regex_new ("(.*)\\.[fF][iI][tT][sS]$", 0, 0, &error);

    g_clear_pointer (&mcbs->base_name, g_free);
    
    if (g_regex_match (fits_ext, filename, 0, &match_info))
    {
      mcbs->base_name = g_match_info_fetch (match_info, 1);
      g_match_info_free (match_info);
    }
    else
      mcbs->base_name = g_strdup (filename);
    g_regex_unref (fits_ext);
    {
      gchar *resample_str = g_strdup_printf ("%s-resample.fits", mcbs->base_name);

      ncm_fit_catalog_set_file (mcbs->fcat, filename);

      ncm_fit_catalog_reset (mcbs->fcat);
      ncm_fit_catalog_erase_data (mcbs->fcat);
      
      if (mcbs->rng != NULL)
      {
        gchar *algo = NULL;
        gulong seed;

        if (ncm_fit_catalog_get_prng (mcbs->fcat, &algo, &seed))
        {
          if (strcmp (algo, ncm_rng_get_algo (mcbs->rng)) != 0 || seed != ncm_rng_get_seed (mcbs->rng))
            g_error ("ncm_fit_mcbs_set_rng: catalog has PRNG with algorithm: `%s' and seed: %lu, and the montecarlo object has algorithm: `%s' and seed: %lu.",
                     algo, seed, ncm_rng_get_algo (mcbs->rng), ncm_rng_get_seed (mcbs->rng));
          g_free (algo);
        }
        else
          ncm_fit_catalog_set_prng (mcbs->fcat, mcbs->rng);
      }
      
      g_free (resample_str);
    }
  }
}

/**
 * ncm_fit_mcbs_set_rng:
 * @mcbs: FIXME
 * @rng: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_mcbs_set_rng (NcmFitMCBS *mcbs, NcmRNG *rng)
{
  gchar *algo = NULL;
  gulong seed = 0;
  
  ncm_rng_clear (&mcbs->rng);
  mcbs->rng = ncm_rng_ref (rng);

  if (ncm_fit_catalog_peek_filename (mcbs->fcat) != NULL)
  {
    if (ncm_fit_catalog_get_prng (mcbs->fcat, &algo, &seed))
    {
      if (strcmp (algo, ncm_rng_get_algo (rng)) != 0 || seed != ncm_rng_get_seed (rng))
        g_error ("ncm_fit_mcbs_set_rng: catalog has PRNG with algorithm: `%s' and seed: %lu, and the montecarlo object has algorithm: `%s' and seed: %lu.",
                 algo, seed, ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));
      g_free (algo);
    }
    else
      ncm_fit_catalog_set_prng (mcbs->fcat, rng);
  }
}

/**
 * ncm_fit_mcbs_run:
 * @mcbs: FIXME
 * @fiduc: FIXME
 * @ni: FIXME
 * @nf: FIXME
 * @nbstraps: FIXME
 * @rtype: FIXME
 * @mtype: FIXME
 * @bsmt: FIXME
 * 
 * FIXME
 * 
 * WARNING not working correctly with bsmt > 0 FIXME
 *
 */
void 
ncm_fit_mcbs_run (NcmFitMCBS *mcbs, NcmMSet *fiduc, guint ni, guint nf, guint nbstraps, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype, guint bsmt)
{
  guint i;
  gboolean cat_has_rng = FALSE;

  ncm_fit_mc_set_rtype (mcbs->mc_resample, NCM_FIT_MC_RESAMPLE_FROM_MODEL);
  ncm_fit_mc_set_mtype (mcbs->mc_resample, NCM_FIT_RUN_MSGS_SIMPLE);
  ncm_fit_mc_set_fiducial (mcbs->mc_resample, fiduc);

  if (mcbs->rng != NULL)
  {
    NcmRNG *rng = ncm_rng_seeded_new (ncm_rng_get_algo (mcbs->rng), ncm_rng_get_seed (mcbs->rng));
    ncm_fit_mc_set_rng (mcbs->mc_resample, rng);
    ncm_rng_free (rng);
    cat_has_rng = TRUE;
  }

  ncm_fit_mc_start_run (mcbs->mc_resample);

  if (!cat_has_rng)
    ncm_fit_catalog_set_prng (mcbs->fcat, mcbs->mc_resample->rng);
  
  if (ni > 0)
    ncm_fit_mc_set_first_sample_id (mcbs->mc_resample, ni);

  ncm_fit_mc_set_rtype (mcbs->mc_bstrap, rtype);
  ncm_fit_mc_set_mtype (mcbs->mc_bstrap, mtype);
  ncm_fit_mc_set_nthreads (mcbs->mc_bstrap, bsmt);
  
  if (rtype == NCM_FIT_MC_RESAMPLE_FROM_MODEL)
    g_error ("ncm_fit_mcbs_run: the internal run must be a bootstrap: NCM_FIT_MC_RESAMPLE_BOOTSTRAP_*.");
  
  for (i = ni; i < nf; i++)
  {
    ncm_dataset_bootstrap_set (mcbs->fit->lh->dset, NCM_DATASET_BSTRAP_DISABLE);
    //ncm_fit_mc_set_first_sample_id (mcbs->mc_resample, i + 1);
    ncm_fit_mc_run (mcbs->mc_resample, i + 1);
    ncm_dataset_bootstrap_set (mcbs->fit->lh->dset, rtype);
    
    if (mcbs->base_name != NULL)
    {
      gchar *bstrap_str = g_strdup_printf ("%s-bstrap-%06d.fits", mcbs->base_name, i);
      ncm_fit_mc_set_data_file (mcbs->mc_bstrap, bstrap_str);
      g_free (bstrap_str);
    }

    ncm_fit_mc_start_run (mcbs->mc_bstrap);
    ncm_fit_mc_run (mcbs->mc_bstrap, nbstraps);

    ncm_fit_catalog_add_from_vector (mcbs->fcat, mcbs->mc_bstrap->fcat->pstats->mean);
    ncm_fit_catalog_log_current_stats (mcbs->fcat);

    ncm_fit_mc_end_run (mcbs->mc_bstrap);
    ncm_fit_mc_reset (mcbs->mc_bstrap);
  }

  ncm_fit_catalog_set_fit_mean_covar (mcbs->fcat);
  ncm_fit_mc_end_run (mcbs->mc_resample);  
}

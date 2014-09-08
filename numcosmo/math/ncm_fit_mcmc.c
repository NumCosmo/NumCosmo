/***************************************************************************
 *            ncm_fit_mcmc.c
 *
 *  Sun May 25 16:51:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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

/**
 * SECTION:ncm_fit_mcmc
 * @title: Markov Chain Monte Carlo Analysis
 * @short_description: Object implementing Markov Chain Monte Carlo analysis
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_mcmc.h"
#include "math/ncm_cfg.h"
#include "math/ncm_func_eval.h"
#include "ncm_enum_types.h"

#include <gsl/gsl_statistics_double.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_SAMPLER,
  PROP_MTYPE,
  PROP_NTHREADS,
  PROP_DATA_FILE,
};

G_DEFINE_TYPE (NcmFitMCMC, ncm_fit_mcmc, G_TYPE_OBJECT);

static void
ncm_fit_mcmc_init (NcmFitMCMC *mcmc)
{
  mcmc->fit             = NULL;
  mcmc->mtype           = NCM_FIT_RUN_MSGS_NONE;
  mcmc->mcs             = NULL;
  mcmc->rng             = NULL;
  mcmc->nt              = ncm_timer_new ();
  mcmc->ser             = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  mcmc->theta           = NULL;
  mcmc->thetastar       = NULL;
  mcmc->nthreads        = 0;
  mcmc->n               = 0;
  mcmc->mp              = NULL;
  mcmc->cur_sample_id   = -1; /* Represents that no samples were calculated yet. */
  mcmc->write_index     = 0;
  mcmc->started         = FALSE;
#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  mcmc->dup_fit         = g_mutex_new ();
  mcmc->resample_lock   = g_mutex_new ();
  mcmc->write_cond      = g_cond_new ();
#else
  g_mutex_init (&mcmc->dup_fit_m);
  mcmc->dup_fit = &mcmc->dup_fit_m;

  g_mutex_init (&mcmc->resample_lock_m);
  mcmc->resample_lock = &mcmc->resample_lock_m;

  g_cond_init (&mcmc->write_cond_m);
  mcmc->write_cond = &mcmc->write_cond_m;
#endif
}

static void _ncm_fit_mcmc_set_fit_obj (NcmFitMCMC *mcmc, NcmFit *fit);

static void
ncm_fit_mcmc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (object);
  g_return_if_fail (NCM_IS_FIT_MCMC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      _ncm_fit_mcmc_set_fit_obj (mcmc, g_value_get_object (value));
      break;
    case PROP_SAMPLER:
      ncm_fit_mcmc_set_sampler (mcmc, g_value_get_object (value));
      break;      
    case PROP_MTYPE:
      ncm_fit_mcmc_set_mtype (mcmc, g_value_get_enum (value));
      break;
    case PROP_NTHREADS:
      ncm_fit_mcmc_set_nthreads (mcmc, g_value_get_uint (value));
      break;
    case PROP_DATA_FILE:
      ncm_fit_mcmc_set_data_file (mcmc, g_value_get_string (value));
      break;    
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_mcmc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (object);
  g_return_if_fail (NCM_IS_FIT_MCMC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, mcmc->fit);
      break;
    case PROP_SAMPLER:
      g_value_set_object (value, mcmc->mcs);
      break;      
    case PROP_MTYPE:
      g_value_set_enum (value, mcmc->mtype);
      break;
    case PROP_NTHREADS:
      g_value_set_uint (value, mcmc->nthreads);
      break;
    case PROP_DATA_FILE:
      g_value_set_string (value, ncm_fit_catalog_peek_filename (mcmc->fcat));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_mcmc_dispose (GObject *object)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (object);

  ncm_fit_clear (&mcmc->fit);
  ncm_timer_clear (&mcmc->nt);
  ncm_serialize_clear (&mcmc->ser);
  ncm_fit_catalog_clear (&mcmc->fcat);
  ncm_rng_clear (&mcmc->rng);
  ncm_vector_clear (&mcmc->theta);
  ncm_vector_clear (&mcmc->thetastar);

  if (mcmc->mp != NULL)
  {
    ncm_memory_pool_free (mcmc->mp, TRUE);
    mcmc->mp = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mcmc_parent_class)->dispose (object);
}

static void
ncm_fit_mcmc_finalize (GObject *object)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (object);

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  g_mutex_free (mcmc->dup_fit);
  g_mutex_free (mcmc->resample_lock);
  g_cond_free (mcmc->write_cond);
#else
  g_mutex_clear (mcmc->dup_fit);
  g_mutex_clear (mcmc->resample_lock);
  g_cond_clear (mcmc->write_cond);
#endif
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mcmc_parent_class)->finalize (object);
}

static void
ncm_fit_mcmc_class_init (NcmFitMCMCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_fit_mcmc_set_property;
  object_class->get_property = &ncm_fit_mcmc_get_property;
  object_class->dispose      = &ncm_fit_mcmc_dispose;
  object_class->finalize     = &ncm_fit_mcmc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SAMPLER,
                                   g_param_spec_object ("sampler",
                                                        NULL,
                                                        "Montecarlo sampler",
                                                        NCM_TYPE_MC_SAMPLER,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_MTYPE,
                                   g_param_spec_enum ("mtype",
                                                      NULL,
                                                      "Run messages type",
                                                      NCM_TYPE_FIT_RUN_MSGS, NCM_FIT_RUN_MSGS_SIMPLE,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NTHREADS,
                                   g_param_spec_uint ("nthreads",
                                                      NULL,
                                                      "Number of threads to run",
                                                      0, 100, 0,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void 
_ncm_fit_mcmc_set_fit_obj (NcmFitMCMC *mcmc, NcmFit *fit)
{
  g_assert (mcmc->fit == NULL);
  mcmc->fit = ncm_fit_ref (fit);
  mcmc->fcat = ncm_fit_catalog_new (fit);
}

/**
 * ncm_fit_mcmc_new:
 * @fit: a #NcmFit
 * @mcs: a #NcmMCSampler.
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitMCMC *
ncm_fit_mcmc_new (NcmFit *fit, NcmMCSampler *mcs, NcmFitRunMsgs mtype)
{
  NcmFitMCMC *mcmc = g_object_new (NCM_TYPE_FIT_MCMC, 
                                "fit", fit,
                                "sampler", mcs,
                                "mtype", mtype,
                                NULL);
  return mcmc;
}

/**
 * ncm_fit_mcmc_free:
 * @mc: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcmc_free (NcmFitMCMC *mcmc)
{
  g_object_unref (mcmc);
}

/**
 * ncm_fit_mcmc_clear:
 * @mc: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcmc_clear (NcmFitMCMC **mcmc)
{
  g_clear_object (mcmc);
}

/**
 * ncm_fit_mcmc_set_data_file:
 * @mc: FIXME
 * @filename: a filename.
 *
 * FIXME
 *
 */
void 
ncm_fit_mcmc_set_data_file (NcmFitMCMC *mcmc, const gchar *filename)
{
  const gchar *cur_filename = ncm_fit_catalog_peek_filename (mcmc->fcat); 
  if (mcmc->started && cur_filename != NULL)
    g_error ("ncm_fit_mcmc_set_data_file: Cannot change data file during a run, call ncm_fit_mcmc_end_run() first.");    

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;

  ncm_fit_catalog_set_file (mcmc->fcat, filename);

  if (mcmc->started)
  {
    gulong seed;
    gchar *algo = NULL;
    if (mcmc->fcat->cur_id > mcmc->cur_sample_id)
    {
//      ncm_fit_mcmc_intern_skip (mcmc, mcmc->fcat->cur_id - mcmc->cur_sample_id);
      g_assert_cmpint (mcmc->cur_sample_id, ==, mcmc->fcat->cur_id);
    }
    else if (mcmc->fcat->cur_id < mcmc->cur_sample_id)
      g_error ("ncm_fit_mcmc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
               mcmc->fcat->cur_id, mcmc->cur_sample_id);

    if (ncm_fit_catalog_get_prng (mcmc->fcat, &algo, &seed))
    {
      g_assert_cmpstr (algo, ==, ncm_rng_get_algo (mcmc->rng));
      g_assert_cmpuint (seed, ==, ncm_rng_get_seed (mcmc->rng));
      g_free (algo);
    }
    else
      ncm_fit_catalog_set_prng (mcmc->fcat, mcmc->rng);
  }
  else
  {
    gulong seed;
    gchar *algo = NULL;
    if (ncm_fit_catalog_get_prng (mcmc->fcat, &algo, &seed))
    {
      if (mcmc->rng == NULL)
      {
        NcmRNG *rng = ncm_rng_seeded_new (algo, seed);
        ncm_fit_mcmc_set_rng (mcmc, rng);
        ncm_rng_free (rng);
      }
      else
      {
        g_assert_cmpstr (algo, ==, ncm_rng_get_algo (mcmc->rng));
        g_assert_cmpuint (seed, ==, ncm_rng_get_seed (mcmc->rng));
      }
      g_free (algo);
    }
    else if (mcmc->rng != NULL)
      ncm_fit_catalog_set_prng (mcmc->fcat, mcmc->rng);
  }
}

/**
 * ncm_fit_mcmc_set_mtype:
 * @mc: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcmc_set_mtype (NcmFitMCMC *mcmc, NcmFitRunMsgs mtype)
{
  mcmc->mtype = mtype;
}

void
ncm_fit_mcmc_metropolis (NcmMSet *mset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  gint i;
  for (i = 0; i < ncm_vector_len (theta); i++)
  {
    gdouble scale = ncm_mset_fparam_get_scale (mset, i);
    ncm_vector_set (thetastar, i, ncm_vector_get (theta, i) + gsl_ran_gaussian (rng->r, scale * 10.0));
  }
}

/**
 * ncm_fit_mcmc_set_mcmctype:
 * @mc: FIXME
 * @mcs: a #NcmMCSampler.
 *
 * FIXME
 *
 */
void 
ncm_fit_mcmc_set_sampler (NcmFitMCMC *mcmc, NcmMCSampler *mcs)
{
  if (mcmc->started)
    g_error ("ncm_fit_mcmc_set_rtype: Cannot change resample type during a run, call ncm_fit_mcmc_end_run() first.");

  mcmc->mcs = ncm_mc_sampler_ref (mcs);

  ncm_fit_catalog_set_run_type (mcmc->fcat, ncm_mc_sampler_get_name (mcs));
  ncm_mc_sampler_set_mset (mcs, mcmc->fit->mset);
}

/**
 * ncm_fit_mcmc_set_nthreads:
 * @mc: FIXME
 * @nthreads: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mcmc_set_nthreads (NcmFitMCMC *mcmc, guint nthreads)
{
  mcmc->nthreads = nthreads;
}

/**
 * ncm_fit_mcmc_set_rng:
 * @mc: FIXME
 * @rng: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_set_rng (NcmFitMCMC *mcmc, NcmRNG *rng)
{
  gchar *algo = NULL;
  gulong seed = 0;
  
  if (mcmc->started)
    g_error ("ncm_fit_mcmc_set_rng: Cannot change the PRNG object during a run, call ncm_fit_mcmc_end_run() first.");

  ncm_rng_clear (&mcmc->rng);
  mcmc->rng = ncm_rng_ref (rng);

  if (ncm_fit_catalog_peek_filename (mcmc->fcat) != NULL)
  {
    if (ncm_fit_catalog_get_prng (mcmc->fcat, &algo, &seed))
    {
      if (strcmp (algo, ncm_rng_get_algo (rng)) != 0 || seed != ncm_rng_get_seed (rng))
        g_error ("ncm_fit_mcmc_set_rng: catalog has PRNG with algorithm: `%s' and seed: %lu, and the montecarlo object has algorithm: `%s' and seed: %lu.",
                 algo, seed, ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));
      g_free (algo);
    }
    else
      ncm_fit_catalog_set_prng (mcmc->fcat, rng);
  }
}

void
_ncm_fit_mcmc_update (NcmFitMCMC *mcmc, NcmFit *fit)
{
  const guint part = 5;
  const guint step = (mcmc->n / part) == 0 ? 1 : (mcmc->n / part);

  ncm_fit_catalog_add_from_fit (mcmc->fcat, fit);
  ncm_timer_task_increment (mcmc->nt);

  switch (mcmc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = mcmc->nt->task_pos % step;
      if ((stepi == 0) || (mcmc->nt->task_pos == mcmc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_fit_catalog_log_current_stats (mcmc->fcat);
        /* ncm_timer_task_accumulate (mcmc->nt, acc); */
        ncm_timer_task_log_elapsed (mcmc->nt);
        ncm_timer_task_log_mean_time (mcmc->nt);
        ncm_timer_task_log_time_left (mcmc->nt);
        ncm_timer_task_log_end_datetime (mcmc->nt);
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      fit->mtype = mcmc->mtype;
      ncm_fit_log_start (fit);
      ncm_fit_log_end (fit);
      ncm_fit_catalog_log_current_stats (mcmc->fcat);
      /* ncm_timer_task_increment (mcmc->nt); */
      ncm_timer_task_log_elapsed (mcmc->nt);
      ncm_timer_task_log_mean_time (mcmc->nt);
      ncm_timer_task_log_time_left (mcmc->nt);
      ncm_timer_task_log_end_datetime (mcmc->nt);
      break;      
  }

  if ((mcmc->fcat->fmode != NCM_FIT_MCMC_MIN_FLUSH_INTERVAL) &&
      (ncm_timer_task_mean_time (mcmc->nt) < NCM_FIT_MCMC_MIN_FLUSH_INTERVAL))
  {
    ncm_fit_catalog_set_flush_mode (mcmc->fcat, NCM_FIT_CATALOG_FLUSH_TIMED);
    ncm_fit_catalog_set_flush_interval (mcmc->fcat, NCM_FIT_MCMC_MIN_FLUSH_INTERVAL);
  }
}

void
_ncm_fit_mcmc_resample_bstrap (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng)
{
  ncm_dataset_bootstrap_resample (dset, rng);
}

/**
 * ncm_fit_mcmc_start_run:
 * @mc: a #NcmFitMCMC
 * 
 * FIXME
 * 
 */
void 
ncm_fit_mcmc_start_run (NcmFitMCMC *mcmc)
{
  if (mcmc->started)
    g_error ("ncm_fit_mcmc_start_run: run already started, run ncm_fit_mcmc_end_run() first.");

  switch (mcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Starting Markov Chain Montecarlo...\n");
      ncm_dataset_log_info (mcmc->fit->lh->dset);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Model set:\n");
      ncm_mset_pretty_log (mcmc->fit->mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (mcmc->rng == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);
    ncm_rng_set_random_seed (rng, FALSE);
    ncm_fit_mcmc_set_rng (mcmc, rng);
    if (mcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmFitMCMC: No PRNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));
    ncm_rng_free (rng);
  }
  
  mcmc->started = TRUE;

  {
    guint fparam_len = ncm_mset_fparam_len (mcmc->fit->mset);
    if (mcmc->theta != NULL)
    {
      ncm_vector_free (mcmc->theta);
      ncm_vector_free (mcmc->thetastar);
    }
    mcmc->theta = ncm_vector_new (fparam_len);
    mcmc->thetastar = ncm_vector_new (fparam_len);
  }
  
  ncm_fit_catalog_sync (mcmc->fcat, TRUE);
  if (mcmc->fcat->cur_id > mcmc->cur_sample_id)
  {
//    ncm_fit_mcmc_intern_skip (mcmc, mcmc->fcat->cur_id - mcmc->cur_sample_id);
    g_assert_cmpint (mcmc->cur_sample_id, ==, mcmc->fcat->cur_id);
  }
  else if (mcmc->fcat->cur_id < mcmc->cur_sample_id)
    g_error ("ncm_fit_mcmc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
             mcmc->fcat->cur_id, mcmc->cur_sample_id);
  
}

/**
 * ncm_fit_mcmc_end_run:
 * @mc: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fit_mcmc_end_run (NcmFitMCMC *mcmc)
{
  if (ncm_timer_task_is_running (mcmc->nt))
    ncm_timer_task_end (mcmc->nt);

  if (mcmc->mp != NULL)
  {
    ncm_memory_pool_free (mcmc->mp, TRUE);
    mcmc->mp = NULL;
  }

  ncm_fit_catalog_sync (mcmc->fcat, TRUE);
  
  mcmc->started = FALSE;
}

/**
 * ncm_fit_mcmc_reset:
 * @mc: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_fit_mcmc_reset (NcmFitMCMC *mcmc)
{
  mcmc->n               = 0;
  mcmc->cur_sample_id   = -1;
  mcmc->write_index     = 0;
  mcmc->started         = FALSE;  
  ncm_fit_catalog_reset (mcmc->fcat);
  ncm_rng_clear (&mcmc->rng);
}

/**
 * ncm_fit_mcmc_set_first_sample_id:
 * @mc: FIXME
 * @first_sample_id: FIXME
 * 
 * FIXME
 *
 */
void 
ncm_fit_mcmc_set_first_sample_id (NcmFitMCMC *mcmc, gint first_sample_id)
{
  if (mcmc->fcat->first_id == first_sample_id)
    return;

  if (!mcmc->started)
    g_error ("ncm_fit_mcmc_set_first_sample_id: run not started, run ncm_fit_mcmc_start_run() first.");

  if (first_sample_id <= mcmc->cur_sample_id)
    g_error ("ncm_fit_mcmc_set_first_sample_id: cannot move first sample id backwards to: %d, catalog first id: %d, current sample id: %d.",
             first_sample_id, mcmc->fcat->first_id, mcmc->cur_sample_id);

  ncm_fit_catalog_set_first_id (mcmc->fcat, first_sample_id);
}

static void _ncm_fit_mcmc_run_single (NcmFitMCMC *mcmc);
static void _ncm_fit_mcmc_run_mt (NcmFitMCMC *mcmc);

/**
 * ncm_fit_mcmc_run:
 * @mc: a #NcmFitMCMC
 * @n: total number of realizations to run
 * 
 * Runs the montecarlo until it reaches the @n-th realization. Note that
 * if the first_id is non-zero it will run @n - first_id realizations.
 *
 */
void 
ncm_fit_mcmc_run (NcmFitMCMC *mcmc, guint n)
{
  if (!mcmc->started)
    g_error ("ncm_fit_mcmc_run: run not started, run ncm_fit_mcmc_start_run() first.");

  if (n <= (mcmc->cur_sample_id + 1))
  {
    if (mcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Nothing to do, current Montecarlo run is %d\n", mcmc->cur_sample_id + 1);
    }
    return;
  }
  
  mcmc->n = n - (mcmc->cur_sample_id + 1);
  
  switch (mcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Calculating [%06d] Markov Chain Montecarlo runs [%s]\n", mcmc->n, ncm_mc_sampler_get_name (mcmc->mcs));
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }
  
  if (ncm_timer_task_is_running (mcmc->nt))
  {
    ncm_timer_task_add_tasks (mcmc->nt, mcmc->n);
    ncm_timer_task_continue (mcmc->nt);
  }
  else
  {
    ncm_timer_task_start (mcmc->nt, mcmc->n);
    ncm_timer_set_name (mcmc->nt, "NcmFitMCMC");
  }
  if (mcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (mcmc->nt);

  if (mcmc->nthreads <= 1)
    _ncm_fit_mcmc_run_single (mcmc);
  else
    _ncm_fit_mcmc_run_mt (mcmc);

  ncm_timer_task_pause (mcmc->nt);
}

static void 
_ncm_fit_mcmc_run_single (NcmFitMCMC *mcmc)
{
  guint i = 0;
  gdouble m2lnL_cur;
  ncm_fit_m2lnL_val (mcmc->fit, &m2lnL_cur);
  
  for (i = 0; i < mcmc->n; i++)
  {
    gdouble m2lnL_star, prob, jump = 0.0;
    
    ncm_mset_fparams_get_vector (mcmc->fit->mset, mcmc->theta);
    ncm_mc_sampler_generate (mcmc->mcs, mcmc->theta, mcmc->thetastar, mcmc->rng);
    ncm_mset_fparams_set_vector (mcmc->fit->mset, mcmc->thetastar);

    ncm_fit_m2lnL_val (mcmc->fit, &m2lnL_star);

/*    
    ncm_vector_log_vals (mcmc->theta, "# Theta  : ", "% 8.5g");
    ncm_vector_log_vals (mcmc->thetastar, "# Theta* : ", "% 8.5g");
*/    
    prob = GSL_MIN (exp ((m2lnL_cur - m2lnL_star) * 0.5), 1.0);
    ncm_fit_state_set_m2lnL_curval (mcmc->fit->fstate, m2lnL_star);
    /*printf ("# Prob %e\n", prob);*/    
    if (prob != 1.0)
    {
      jump = gsl_rng_uniform (mcmc->rng->r);
      if (jump > prob)
      {
        ncm_mset_fparams_set_vector (mcmc->fit->mset, mcmc->theta);
        ncm_fit_state_set_m2lnL_curval (mcmc->fit->fstate, m2lnL_cur);
      }
    }

    m2lnL_cur = ncm_fit_state_get_m2lnL_curval (mcmc->fit->fstate);
    
    _ncm_fit_mcmc_update (mcmc, mcmc->fit);
    mcmc->write_index++;
  }
}

static gpointer
_ncm_fit_mcmc_dup_fit (gpointer userdata)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (userdata);
  g_mutex_lock (mcmc->dup_fit);
  {
    NcmFit *fit = ncm_fit_dup (mcmc->fit, mcmc->ser);
    ncm_serialize_reset (mcmc->ser);
    g_mutex_unlock (mcmc->dup_fit);
    return fit;
  }
}

static void 
_ncm_fit_mcmc_mt_eval (glong i, glong f, gpointer data)
{
  _NCM_STATIC_MUTEX_DECL (update_lock);
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (data);
  NcmFit **fit_ptr = ncm_memory_pool_get (mcmc->mp);
  NcmFit *fit = *fit_ptr;
  guint j;

  for (j = i; j < f; j++)
  {
    gint sample_index;

//    ncm_mset_param_set_vector (fit->mset, mcmc->bf);

    g_mutex_lock (mcmc->resample_lock);
//    sample_index = _ncm_fit_mcmc_resample (mc, fit);
    sample_index = 0;
    g_mutex_unlock (mcmc->resample_lock);

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

    _NCM_MUTEX_LOCK (&update_lock);    
    while (mcmc->write_index != sample_index)
      g_cond_wait (mcmc->write_cond, &update_lock);

    _ncm_fit_mcmc_update (mcmc, fit);
    mcmc->write_index++;
    g_cond_broadcast (mcmc->write_cond);

    _NCM_MUTEX_UNLOCK (&update_lock);
  }

  ncm_memory_pool_return (fit_ptr);
}

static void
_ncm_fit_mcmc_run_mt (NcmFitMCMC *mcmc)
{
  const guint nthreads = mcmc->n > mcmc->nthreads ? mcmc->nthreads : (mcmc->n - 1);

  if (nthreads == 0)
  {
    _ncm_fit_mcmc_run_single (mcmc);
    return;
  }
  
  if (mcmc->mp != NULL)
    ncm_memory_pool_free (mcmc->mp, TRUE);
  mcmc->mp = ncm_memory_pool_new (&_ncm_fit_mcmc_dup_fit, mcmc, 
                                (GDestroyNotify) &ncm_fit_free);

  /*
   * The line below added de main fit object to the pool, but can cause
   * several race conditions as it is used to make the copies for the other
   * threads. So, no.
   */
  //ncm_memory_pool_add (mc->mp, mc->fit);

  g_assert_cmpuint (mcmc->nthreads, >, 1);

  ncm_func_eval_threaded_loop_full (&_ncm_fit_mcmc_mt_eval, 0, mcmc->n, mcmc);  
}

/**
 * ncm_fit_mcmc_run_lre:
 * @mc: FIXME
 * @prerun: FIXME
 * @lre: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_fit_mcmc_run_lre (NcmFitMCMC *mcmc, guint prerun, gdouble lre)
{
  gdouble lerror;
  const gdouble lre2 = lre * lre;

  g_assert_cmpfloat (lre, >, 0.0);
  /* g_assert_cmpfloat (lre, <, 1.0); */

  if (prerun == 0)
  {
    guint fparam_len = ncm_mset_fparam_len (mcmc->fit->mset);
    prerun = fparam_len * 100;
  }

  if (ncm_fit_catalog_len (mcmc->fcat) < prerun)
  {
    guint prerun_left = prerun - ncm_fit_catalog_len (mcmc->fcat);
    if (mcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitMCMC: Running first %u pre-runs...\n", prerun_left);
    ncm_fit_mcmc_run (mcmc, prerun);
  }

  lerror = ncm_fit_catalog_largest_error (mcmc->fcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n = ncm_fit_catalog_len (mcmc->fcat);
    gdouble m = n * lerror2 / lre2;
    guint runs = ((m - n) > 1000.0) ? ceil ((m - n) * 0.25) : ceil (m - n);

    if (mcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmFitMCMC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmFitMCMC: Running more %u runs...\n", runs);
    }
    ncm_fit_mcmc_run (mcmc, mcmc->cur_sample_id + runs + 1);
    lerror = ncm_fit_catalog_largest_error (mcmc->fcat);
  }

  if (mcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    g_message ("# NcmFitMCMC: Largest relative error %e attained: %e\n", lre, lerror);
}

/**
 * ncm_fit_mcmc_mean_covar:
 * @mc: FIXME
 *
 * FIXME
 */
void
ncm_fit_mcmc_mean_covar (NcmFitMCMC *mcmc)
{
  ncm_fit_catalog_set_fit_mean_covar (mcmc->fcat);
}

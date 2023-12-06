/***************************************************************************
 *            ncm_fit_mcmc.c
 *
 *  Sun May 25 16:51:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_fit_mcmc
 * @title: NcmFitMCMC
 * @short_description: Markov Chain Monte Carlo analysis.
 *
 * FIXME
 *
 * Metropolis–Hastings sampler.
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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_statistics_double.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_FIT,
  PROP_SAMPLER,
  PROP_MTYPE,
  PROP_NTHREADS,
  PROP_DATA_FILE,
};

struct _NcmFitMCMC
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmMSetCatalog *mcat;
  NcmFitRunMsgs mtype;
  NcmTimer *nt;
  NcmSerialize *ser;
  NcmMSetTransKern *tkern;
  NcmVector *theta;
  NcmVector *thetastar;
  guint nthreads;
  guint n;
  NcmMemoryPool *mp;
  gint write_index;
  gint cur_sample_id;
  guint naccepted;
  guint ntotal;
  gboolean started;
  GMutex dup_fit;
  GMutex resample_lock;
  GMutex update_lock;
  GCond write_cond;
};

G_DEFINE_TYPE (NcmFitMCMC, ncm_fit_mcmc, G_TYPE_OBJECT)

static void
ncm_fit_mcmc_init (NcmFitMCMC *mcmc)
{
  mcmc->fit           = NULL;
  mcmc->mtype         = NCM_FIT_RUN_MSGS_NONE;
  mcmc->tkern         = NULL;
  mcmc->nt            = ncm_timer_new ();
  mcmc->ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  mcmc->theta         = NULL;
  mcmc->thetastar     = NULL;
  mcmc->nthreads      = 0;
  mcmc->n             = 0;
  mcmc->mp            = NULL;
  mcmc->cur_sample_id = -1; /* Represents that no samples were calculated yet. */
  mcmc->naccepted     = 0;
  mcmc->ntotal        = 0;
  mcmc->write_index   = 0;
  mcmc->started       = FALSE;

  g_mutex_init (&mcmc->dup_fit);
  g_mutex_init (&mcmc->resample_lock);
  g_mutex_init (&mcmc->update_lock);
  g_cond_init (&mcmc->write_cond);
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
      ncm_fit_mcmc_set_trans_kern (mcmc, g_value_get_object (value));
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
      g_value_set_object (value, mcmc->tkern);
      break;
    case PROP_MTYPE:
      g_value_set_enum (value, mcmc->mtype);
      break;
    case PROP_NTHREADS:
      g_value_set_uint (value, mcmc->nthreads);
      break;
    case PROP_DATA_FILE:
      g_value_set_string (value, ncm_mset_catalog_peek_filename (mcmc->mcat));
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
  ncm_mset_trans_kern_clear (&mcmc->tkern);
  ncm_timer_clear (&mcmc->nt);
  ncm_serialize_clear (&mcmc->ser);
  ncm_mset_catalog_clear (&mcmc->mcat);
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

  g_mutex_clear (&mcmc->dup_fit);
  g_mutex_clear (&mcmc->resample_lock);
  g_mutex_clear (&mcmc->update_lock);
  g_cond_clear (&mcmc->write_cond);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mcmc_parent_class)->finalize (object);
}

static void
ncm_fit_mcmc_class_init (NcmFitMCMCClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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
                                                        "Metropolis–Hastings sampler",
                                                        NCM_TYPE_MSET_TRANS_KERN,
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
  NcmMSet *mset = ncm_fit_peek_mset (fit);

  g_assert (mcmc->fit == NULL);
  mcmc->fit  = ncm_fit_ref (fit);
  mcmc->mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE,
                                     NCM_MSET_CATALOG_M2LNL_COLNAME, NCM_MSET_CATALOG_M2LNL_SYMBOL,
                                     NULL);
  ncm_mset_catalog_set_m2lnp_var (mcmc->mcat, 0);
}

/**
 * ncm_fit_mcmc_new:
 * @fit: a #NcmFit
 * @tkern: a #NcmMSetTransKern.
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitMCMC *
ncm_fit_mcmc_new (NcmFit *fit, NcmMSetTransKern *tkern, NcmFitRunMsgs mtype)
{
  NcmFitMCMC *mcmc = g_object_new (NCM_TYPE_FIT_MCMC,
                                   "fit", fit,
                                   "sampler", tkern,
                                   "mtype", mtype,
                                   NULL);

  return mcmc;
}

/**
 * ncm_fit_mcmc_free:
 * @mcmc: a #NcmFitMCMC
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
 * @mcmc: a #NcmFitMCMC
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
 * @mcmc: a #NcmFitMCMC
 * @filename: a filename.
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_set_data_file (NcmFitMCMC *mcmc, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (mcmc->mcat);

  if (mcmc->started && (cur_filename != NULL))
    g_error ("ncm_fit_mcmc_set_data_file: Cannot change data file during a run, call ncm_fit_mcmc_end_run() first.");

  if ((cur_filename != NULL) && (strcmp (cur_filename, filename) == 0))
    return;

  ncm_mset_catalog_set_file (mcmc->mcat, filename);

  if (mcmc->started)
    g_assert_cmpint (mcmc->cur_sample_id, ==, ncm_mset_catalog_get_cur_id (mcmc->mcat));
}

/**
 * ncm_fit_mcmc_set_mtype:
 * @mcmc: a #NcmFitMCMC
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

/**
 * ncm_fit_mcmc_set_trans_kern:
 * @mcmc: a #NcmFitMCMC
 * @tkern: a #NcmMSetTransKern.
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_set_trans_kern (NcmFitMCMC *mcmc, NcmMSetTransKern *tkern)
{
  NcmMSet *mset = ncm_fit_peek_mset (mcmc->fit);

  if (mcmc->started)
    g_error ("ncm_fit_mcmc_set_rtype: Cannot change resample type during a run, call ncm_fit_mcmc_end_run() first.");

  mcmc->tkern = ncm_mset_trans_kern_ref (tkern);

  ncm_mset_catalog_set_run_type (mcmc->mcat, ncm_mset_trans_kern_get_name (tkern));
  ncm_mset_trans_kern_set_mset (tkern, mset);
}

/**
 * ncm_fit_mcmc_set_nthreads:
 * @mcmc: a #NcmFitMCMC
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
 * @mcmc: a #NcmFitMCMC
 * @rng: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_set_rng (NcmFitMCMC *mcmc, NcmRNG *rng)
{
  if (mcmc->started)
    g_error ("ncm_fit_mcmc_set_rng: Cannot change the RNG object during a run, call ncm_fit_mcmc_end_run() first.");

  ncm_mset_catalog_set_rng (mcmc->mcat, rng);
}

/**
 * ncm_fit_mcmc_get_accept_ratio:
 * @mcmc: a #NcmFitMCMC
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_fit_mcmc_get_accept_ratio (NcmFitMCMC *mcmc)
{
  return mcmc->naccepted * 1.0 / (mcmc->ntotal * 1.0);
}

void
_ncm_fit_mcmc_update (NcmFitMCMC *mcmc, NcmFit *fit)
{
  NcmMSet *mset       = ncm_fit_peek_mset (fit);
  NcmFitState *fstate = ncm_fit_peek_state (fit);
  const guint part    = 5;
  const guint step    = (mcmc->n / part) == 0 ? 1 : (mcmc->n / part);

  ncm_mset_catalog_add_from_mset (mcmc->mcat, mset, ncm_fit_state_get_m2lnL_curval (fstate), NULL);
  ncm_timer_task_increment (mcmc->nt);

  switch (mcmc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi          = mcmc->nt->task_pos % step;
      gboolean log_timeout = FALSE;

      if ((mcmc->nt->pos_time - mcmc->nt->last_log_time) > 60.0)
        log_timeout = TRUE;

      if (log_timeout || (stepi == 0) || (mcmc->nt->task_pos == mcmc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (mcmc->mcat);
        g_message ("# NcmFitMCMC:acceptance ratio %7.4f%%.\n", ncm_fit_mcmc_get_accept_ratio (mcmc) * 100.0);

        /* ncm_timer_task_accumulate (mcmc->nt, acc); */
        ncm_timer_task_log_elapsed (mcmc->nt);
        ncm_timer_task_log_mean_time (mcmc->nt);
        ncm_timer_task_log_time_left (mcmc->nt);
        ncm_timer_task_log_cur_datetime (mcmc->nt);
        ncm_timer_task_log_end_datetime (mcmc->nt);
      }

      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_fit_set_messages (fit, mcmc->mtype);
      ncm_fit_log_state (fit);
      ncm_mset_catalog_log_current_stats (mcmc->mcat);
      g_message ("# NcmFitMCMC:acceptance ratio %7.4f%%.\n", ncm_fit_mcmc_get_accept_ratio (mcmc) * 100.0);
      /* ncm_timer_task_increment (mcmc->nt); */
      ncm_timer_task_log_elapsed (mcmc->nt);
      ncm_timer_task_log_mean_time (mcmc->nt);
      ncm_timer_task_log_time_left (mcmc->nt);
      ncm_timer_task_log_cur_datetime (mcmc->nt);
      ncm_timer_task_log_end_datetime (mcmc->nt);
      break;
  }
}

static void ncm_fit_mcmc_intern_skip (NcmFitMCMC *mcmc, guint n);

/**
 * ncm_fit_mcmc_start_run:
 * @mcmc: a #NcmFitMCMC
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_start_run (NcmFitMCMC *mcmc)
{
  NcmMSet *mset       = ncm_fit_peek_mset (mcmc->fit);
  NcmLikelihood *lh   = ncm_fit_peek_likelihood (mcmc->fit);
  NcmFitState *fstate = ncm_fit_peek_state (mcmc->fit);
  NcmDataset *dset    = ncm_likelihood_peek_dataset (lh);
  const gint cur_id   = ncm_mset_catalog_get_cur_id (mcmc->mcat);

  if (mcmc->started)
    g_error ("ncm_fit_mcmc_start_run: run already started, run ncm_fit_mcmc_end_run() first.");

  switch (mcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Starting Markov Chain Monte Carlo...\n");
      ncm_dataset_log_info (dset);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Model set:\n");
      ncm_mset_pretty_log (mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_mset_catalog_peek_rng (mcmc->mcat) == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);

    ncm_rng_set_random_seed (rng, FALSE);
    ncm_fit_mcmc_set_rng (mcmc, rng);

    if (mcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmFitMCMC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));

    ncm_rng_free (rng);
  }

  mcmc->started = TRUE;

  {
    guint fparam_len = ncm_mset_fparam_len (mset);

    if (mcmc->theta != NULL)
    {
      ncm_vector_free (mcmc->theta);
      ncm_vector_free (mcmc->thetastar);
    }

    mcmc->theta     = ncm_vector_new (fparam_len);
    mcmc->thetastar = ncm_vector_new (fparam_len);
  }

  mcmc->naccepted = 0;
  mcmc->ntotal    = 0;

  ncm_mset_catalog_set_sync_mode (mcmc->mcat, NCM_MSET_CATALOG_SYNC_TIMED);
  ncm_mset_catalog_set_sync_interval (mcmc->mcat, NCM_FIT_MCMC_MIN_SYNC_INTERVAL);

  ncm_mset_catalog_sync (mcmc->mcat, TRUE);

  if (cur_id > mcmc->cur_sample_id)
  {
    ncm_fit_mcmc_intern_skip (mcmc, cur_id - mcmc->cur_sample_id);
    g_assert_cmpint (mcmc->cur_sample_id, ==, cur_id);
  }
  else if (cur_id < mcmc->cur_sample_id)
  {
    g_error ("ncm_fit_mcmc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].",
             cur_id, mcmc->cur_sample_id);
  }

  {
    NcmVector *cur_row = NULL;

    cur_row = ncm_mset_catalog_peek_current_row (mcmc->mcat);

    if (cur_row != NULL)
    {
      ncm_mset_fparams_set_vector_offset (mset, cur_row, 1);
      ncm_fit_state_set_m2lnL_curval (fstate, ncm_vector_get (cur_row, 0));
    }
    else
    {
      gdouble m2lnL = 0.0;

      ncm_fit_m2lnL_val (mcmc->fit, &m2lnL);
      ncm_fit_state_set_m2lnL_curval (fstate, m2lnL);
    }
  }
}

/**
 * ncm_fit_mcmc_end_run:
 * @mcmc: a #NcmFitMCMC
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

  ncm_mset_catalog_sync (mcmc->mcat, TRUE);

  mcmc->started = FALSE;
}

/**
 * ncm_fit_mcmc_reset:
 * @mcmc: a #NcmFitMCMC
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_reset (NcmFitMCMC *mcmc)
{
  mcmc->n             = 0;
  mcmc->cur_sample_id = -1;
  mcmc->write_index   = 0;
  mcmc->ntotal        = 0;
  mcmc->naccepted     = 0;
  mcmc->started       = FALSE;
  ncm_mset_catalog_reset (mcmc->mcat);
}

static void
ncm_fit_mcmc_intern_skip (NcmFitMCMC *mcmc, guint n)
{
  if (n == 0)
    return;

  switch (mcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Skipping %u tries, will start at %u-th try.\n", n, mcmc->cur_sample_id + n + 1 + 1);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  mcmc->cur_sample_id += n;
  mcmc->write_index    = mcmc->cur_sample_id + 1;
}

/**
 * ncm_fit_mcmc_set_first_sample_id:
 * @mcmc: a #NcmFitMCMC
 * @first_sample_id: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_set_first_sample_id (NcmFitMCMC *mcmc, gint first_sample_id)
{
  const gint first_id = ncm_mset_catalog_get_first_id (mcmc->mcat);

  if (first_id == first_sample_id)
    return;

  if (!mcmc->started)
    g_error ("ncm_fit_mcmc_set_first_sample_id: run not started, run ncm_fit_mcmc_start_run() first.");

  if (first_sample_id <= mcmc->cur_sample_id)
    g_error ("ncm_fit_mcmc_set_first_sample_id: cannot move first sample id backwards to: %d, catalog first id: %d, current sample id: %d.",
             first_sample_id, first_id, mcmc->cur_sample_id);

  ncm_mset_catalog_set_first_id (mcmc->mcat, first_sample_id);
}

static void _ncm_fit_mcmc_run_single (NcmFitMCMC *mcmc);
static void _ncm_fit_mcmc_run_mt (NcmFitMCMC *mcmc);

/**
 * ncm_fit_mcmc_run:
 * @mcmc: a #NcmFitMCMC
 * @n: total number of realizations to run
 *
 * Runs the Markov Chain Monte Carlo until it reaches the @n-th realization. Note that
 * if the first_id is non-zero it will run @n - first_id realizations.
 *
 */
void
ncm_fit_mcmc_run (NcmFitMCMC *mcmc, guint n)
{
  if (!mcmc->started)
    g_error ("ncm_fit_mcmc_run: run not started, run ncm_fit_mcmc_start_run() first.");

  if (n <= (guint) (mcmc->cur_sample_id + 1))
  {
    if (mcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMCMC: Nothing to do, current Monte Carlo run is %d\n", mcmc->cur_sample_id + 1);
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
      g_message ("# NcmFitMCMC: Calculating [%06d] Markov Chain Monte Carlo runs [%s]\n", mcmc->n, ncm_mset_trans_kern_get_name (mcmc->tkern));
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
  NcmFitState *fstate = ncm_fit_peek_state (mcmc->fit);
  NcmMSet *mset       = ncm_fit_peek_mset (mcmc->fit);
  NcmRNG *rng         = ncm_mset_catalog_peek_rng (mcmc->mcat);
  guint i             = 0;

  for (i = 0; i < mcmc->n; i++)
  {
    gdouble m2lnL_cur = ncm_fit_state_get_m2lnL_curval (fstate);
    gdouble m2lnL_star, prob, jump = 0.0;

    ncm_mset_fparams_get_vector (mset, mcmc->theta);
    ncm_mset_trans_kern_generate (mcmc->tkern, mcmc->theta, mcmc->thetastar, rng);
    ncm_mset_fparams_set_vector (mset, mcmc->thetastar);

    ncm_fit_m2lnL_val (mcmc->fit, &m2lnL_star);
    mcmc->ntotal++;
    mcmc->naccepted++;

/*
 *   ncm_vector_log_vals (mcmc->theta, "# Theta  : ", "% 8.5g");
 *   ncm_vector_log_vals (mcmc->thetastar, "# Theta* : ", "% 8.5g");
 */
    prob = GSL_MIN (exp ((m2lnL_cur - m2lnL_star) * 0.5), 1.0);
    ncm_fit_state_set_m2lnL_curval (fstate, m2lnL_star);

    /*printf ("# Prob %e [% 21.16g % 21.16g] % 21.16g\n", prob, m2lnL_cur, m2lnL_star, m2lnL_cur - m2lnL_star);*/

    if (prob != 1.0)
    {
      jump = gsl_rng_uniform (rng->r);

      if (jump > prob)
      {
        ncm_mset_fparams_set_vector (mset, mcmc->theta);
        ncm_fit_state_set_m2lnL_curval (fstate, m2lnL_cur);
        mcmc->naccepted--;
      }
    }

    _ncm_fit_mcmc_update (mcmc, mcmc->fit);
    mcmc->write_index++;
  }
}

static gpointer
_ncm_fit_mcmc_dup_fit (gpointer userdata)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (userdata);

  g_mutex_lock (&mcmc->dup_fit);
  {
    NcmFit *fit = ncm_fit_dup (mcmc->fit, mcmc->ser);

    ncm_serialize_reset (mcmc->ser, TRUE);
    g_mutex_unlock (&mcmc->dup_fit);

    return fit;
  }
}

static void
_ncm_fit_mcmc_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitMCMC *mcmc = NCM_FIT_MCMC (data);
  NcmFit **fit_ptr = ncm_memory_pool_get (mcmc->mp);
  NcmFit *fit      = *fit_ptr;
  guint j;

  g_assert_not_reached (); /* FIXME not implemented */

  for (j = i; j < f; j++)
  {
    gint sample_index;

/*    ncm_mset_param_set_vector (fit->mset, mcmc->bf); */

    g_mutex_lock (&mcmc->resample_lock);
/*    sample_index = _ncm_fit_mcmc_resample (mc, fit); */
    sample_index = 0;
    g_mutex_unlock (&mcmc->resample_lock);

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

    g_mutex_lock (&mcmc->update_lock);

    while (mcmc->write_index != sample_index)
      g_cond_wait (&mcmc->write_cond, &mcmc->update_lock);

    _ncm_fit_mcmc_update (mcmc, fit);
    mcmc->write_index++;
    g_cond_broadcast (&mcmc->write_cond);

    g_mutex_unlock (&mcmc->update_lock);
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
                                  (GDestroyNotify) & ncm_fit_free);

  /*
   * The line below added de main fit object to the pool, but can cause
   * several race conditions as it is used to make the copies for the other
   * threads. So, no.
   */
  /*ncm_memory_pool_add (mc->mp, mc->fit); */

  g_assert_cmpuint (mcmc->nthreads, >, 1);

  ncm_func_eval_threaded_loop_full (&_ncm_fit_mcmc_mt_eval, 0, mcmc->n, mcmc);
}

/**
 * ncm_fit_mcmc_run_lre:
 * @mcmc: a #NcmFitMCMC
 * @prerun: FIXME
 * @lre: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_mcmc_run_lre (NcmFitMCMC *mcmc, guint prerun, gdouble lre)
{
  NcmMSet *mset          = ncm_fit_peek_mset (mcmc->fit);
  const gdouble lre2     = lre * lre;
  const guint fparam_len = ncm_mset_fparam_len (mset);
  gdouble lerror;

  g_assert_cmpfloat (lre, >, 0.0);

  prerun = GSL_MAX (prerun, fparam_len * 1000);

  if (ncm_mset_catalog_len (mcmc->mcat) < prerun)
  {
    guint prerun_left = prerun - ncm_mset_catalog_len (mcmc->mcat);

    if (mcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitMCMC: Running first %u pre-runs...\n", prerun_left);

    ncm_fit_mcmc_run (mcmc, prerun);
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (mcmc->mcat, FALSE);
  lerror = ncm_mset_catalog_largest_error (mcmc->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n             = ncm_mset_catalog_len (mcmc->mcat);
    gdouble m             = n * lerror2 / lre2;
    guint runs            = ((m - n) > 1000.0) ? ceil ((m - n) * 1.0e-1) : ceil (m - n);

    if (mcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmFitMCMC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmFitMCMC: Running more %u runs...\n", runs);
    }

    ncm_fit_mcmc_run (mcmc, mcmc->cur_sample_id + runs + 1);
    ncm_mset_catalog_estimate_autocorrelation_tau (mcmc->mcat, FALSE);
    lerror = ncm_mset_catalog_largest_error (mcmc->mcat);
  }

  if (mcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    g_message ("# NcmFitMCMC: Largest relative error %e attained: %e\n", lre, lerror);
}

/**
 * ncm_fit_mcmc_mean_covar:
 * @mcmc: a #NcmFitMCMC
 *
 * Computes the mean and covariance of the generated catalog.
 *
 */
void
ncm_fit_mcmc_mean_covar (NcmFitMCMC *mcmc)
{
  NcmFitState *fstate = ncm_fit_peek_state (mcmc->fit);
  NcmVector *fparams  = ncm_fit_state_peek_fparams (fstate);
  NcmMatrix *covar    = ncm_fit_state_peek_covar (fstate);

  ncm_mset_catalog_get_mean (mcmc->mcat, &fparams);
  ncm_mset_catalog_get_covar (mcmc->mcat, &covar);
  ncm_mset_fparams_set_vector (ncm_mset_catalog_peek_mset (mcmc->mcat), fparams);

  ncm_fit_state_set_has_covar (fstate, TRUE);
}

/**
 * ncm_fit_mcmc_get_catalog:
 * @mcmc: a #NcmFitMCMC
 *
 * Gets the generated catalog of @mcmc.
 *
 * Returns: (transfer full): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_mcmc_get_catalog (NcmFitMCMC *mcmc)
{
  return ncm_mset_catalog_ref (mcmc->mcat);
}


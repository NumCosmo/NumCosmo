/***************************************************************************
 *            ncm_fit_mc.c
 *
 *  Sat December 01 17:19:10 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_fit_mc
 * @title: NcmFitMC
 * @short_description: Monte Carlo analysis.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_mc.h"
#include "math/ncm_cfg.h"
#include "math/ncm_func_eval.h"
#include "ncm_enum_types.h"

#include <gsl/gsl_statistics_double.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_RTYPE,
  PROP_FIDUC,
  PROP_MTYPE,
  PROP_NTHREADS,
  PROP_KEEP_ORDER,
  PROP_DATA_FILE,
};

G_DEFINE_TYPE (NcmFitMC, ncm_fit_mc, G_TYPE_OBJECT);

static void
ncm_fit_mc_init (NcmFitMC *mc)
{
  mc->fit             = NULL;
  mc->fiduc           = NULL;
  mc->mcat            = NULL;
  mc->mtype           = NCM_FIT_RUN_MSGS_NONE;
  mc->rtype           = NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN;
  mc->bf              = NULL;
  mc->nt              = ncm_timer_new ();
  mc->ser             = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  mc->nthreads        = 0;
  mc->n               = 0;
  mc->keep_order      = FALSE;
  mc->mp              = NULL;
  mc->cur_sample_id   = -1; /* Represents that no samples were calculated yet. */
  mc->write_index     = 0;
  mc->started         = FALSE;
#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  mc->dup_fit         = g_mutex_new ();
  mc->resample_lock   = g_mutex_new ();
  mc->write_cond      = g_cond_new ();
#else
  g_mutex_init (&mc->dup_fit_m);
  mc->dup_fit = &mc->dup_fit_m;

  g_mutex_init (&mc->resample_lock_m);
  mc->resample_lock = &mc->resample_lock_m;

  g_cond_init (&mc->write_cond_m);
  mc->write_cond = &mc->write_cond_m;
#endif
}

static void _ncm_fit_mc_set_fit_obj (NcmFitMC *mc, NcmFit *fit);

static void
ncm_fit_mc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitMC *mc = NCM_FIT_MC (object);
  g_return_if_fail (NCM_IS_FIT_MC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      _ncm_fit_mc_set_fit_obj (mc, g_value_get_object (value));
      break;
    case PROP_RTYPE:
      ncm_fit_mc_set_rtype (mc, g_value_get_enum (value));
      break;      
    case PROP_FIDUC:
      ncm_fit_mc_set_fiducial (mc, g_value_get_object (value));
      break;
    case PROP_MTYPE:
      ncm_fit_mc_set_mtype (mc, g_value_get_enum (value));
      break;
    case PROP_NTHREADS:
      ncm_fit_mc_set_nthreads (mc, g_value_get_uint (value));
      break;
    case PROP_KEEP_ORDER:
      ncm_fit_mc_keep_order (mc, g_value_get_boolean (value));
      break;
    case PROP_DATA_FILE:
      ncm_fit_mc_set_data_file (mc, g_value_get_string (value));
      break;    
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    case PROP_NTHREADS:
      g_value_set_uint (value, mc->nthreads);
      break;
    case PROP_KEEP_ORDER:
      g_value_set_boolean (value, mc->keep_order);
      break;
    case PROP_DATA_FILE:
      g_value_set_string (value, ncm_mset_catalog_peek_filename (mc->mcat));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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

  if (mc->mp != NULL)
  {
    ncm_memory_pool_free (mc->mp, TRUE);
    mc->mp = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->dispose (object);
}

static void
ncm_fit_mc_finalize (GObject *object)
{
  NcmFitMC *mc = NCM_FIT_MC (object);

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  g_mutex_free (mc->dup_fit);
  g_mutex_free (mc->resample_lock);
  g_cond_free (mc->write_cond);
#else
  g_mutex_clear (mc->dup_fit);
  g_mutex_clear (mc->resample_lock);
  g_cond_clear (mc->write_cond);
#endif
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->finalize (object);
}

static void
ncm_fit_mc_class_init (NcmFitMCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_fit_mc_set_property;
  object_class->get_property = &ncm_fit_mc_get_property;
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
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
                                   PROP_KEEP_ORDER,
                                   g_param_spec_boolean ("keep-order",
                                                         NULL,
                                                         "Whether keep the catalog in order of sampling",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NTHREADS,
                                   g_param_spec_uint ("nthreads",
                                                      NULL,
                                                      "Number of threads to run",
                                                      0, 100, 0,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void 
_ncm_fit_mc_set_fit_obj (NcmFitMC *mc, NcmFit *fit)
{
  g_assert (mc->fit == NULL);
  mc->fit = ncm_fit_ref (fit);
  mc->mcat = ncm_mset_catalog_new (fit->mset, 1, 1, FALSE, NCM_MSET_CATALOG_M2LNL_COLNAME);
}

/**
 * ncm_fit_mc_new:
 * @fit: FIXME
 * @rtype: FIXME
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitMC *
ncm_fit_mc_new (NcmFit *fit, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype)
{
  NcmFitMC * mc = g_object_new (NCM_TYPE_FIT_MC, 
                                "fit", fit,
                                "rtype", rtype,
                                "mtype", mtype,
                                NULL);
  return mc;
}

/**
 * ncm_fit_mc_free:
 * @mc: a #NcmFitMC
 *
 * FIXME
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
 * FIXME
 *
 */
void 
ncm_fit_mc_clear (NcmFitMC **mc)
{
  g_clear_object (mc);
}

static void ncm_fit_mc_intern_skip (NcmFitMC *mc, guint n);

/**
 * ncm_fit_mc_set_data_file:
 * @mc: a #NcmFitMC
 * @filename: a filename.
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_set_data_file (NcmFitMC *mc, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (mc->mcat); 
  if (mc->started && cur_filename != NULL)
    g_error ("ncm_fit_mc_set_data_file: Cannot change data file during a run, call ncm_fit_mc_end_run() first.");    

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;

  ncm_mset_catalog_set_file (mc->mcat, filename);

  if (mc->started)
  {
    if (mc->mcat->cur_id > mc->cur_sample_id)
    {
      ncm_fit_mc_intern_skip (mc, mc->mcat->cur_id - mc->cur_sample_id);
      g_assert_cmpint (mc->cur_sample_id, ==, mc->mcat->cur_id);
    }
    else if (mc->mcat->cur_id < mc->cur_sample_id)
    {
      g_error ("ncm_fit_mc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
               mc->mcat->cur_id, mc->cur_sample_id);
    }
  }
}

/**
 * ncm_fit_mc_set_mtype:
 * @mc: a #NcmFitMC
 * @mtype: FIXME
 *
 * FIXME
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
  mc->resample (fit->lh->dset, mc->fiduc, mc->mcat->rng);
  mc->cur_sample_id++;
  return mc->cur_sample_id;
}

/**
 * ncm_fit_mc_set_rtype:
 * @mc: a #NcmFitMC
 * @rtype: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_set_rtype (NcmFitMC *mc, NcmFitMCResampleType rtype)
{
  const GEnumValue *eval = ncm_cfg_enum_get_value (NCM_TYPE_FIT_MC_RESAMPLE_TYPE, rtype);

  if (mc->started)
    g_error ("ncm_fit_mc_set_rtype: Cannot change resample type during a run, call ncm_fit_mc_end_run() first.");

  mc->rtype = rtype;

  ncm_mset_catalog_set_run_type (mc->mcat, eval->value_nick);
  
  switch (rtype)
  {
    case NCM_FIT_MC_RESAMPLE_FROM_MODEL:
      mc->resample = &ncm_dataset_resample;
      ncm_dataset_bootstrap_set (mc->fit->lh->dset, NCM_DATASET_BSTRAP_DISABLE);
      break;
    case NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX:
      mc->resample = &_ncm_fit_mc_resample_bstrap;
      ncm_dataset_bootstrap_set (mc->fit->lh->dset, NCM_DATASET_BSTRAP_PARTIAL);
      break;
    case NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX:
      mc->resample = &_ncm_fit_mc_resample_bstrap;
      ncm_dataset_bootstrap_set (mc->fit->lh->dset, NCM_DATASET_BSTRAP_TOTAL);
      break;
    case NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_fit_mc_set_nthreads:
 * @mc: a #NcmFitMC
 * @nthreads: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_set_nthreads (NcmFitMC *mc, guint nthreads)
{
  mc->nthreads = nthreads;
}

/**
 * ncm_fit_mc_keep_order:
 * @mc: a #NcmFitMC
 * @keep_order: FIXME
 *
 * FIXME
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
 * @fiduc: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_mc_set_fiducial (NcmFitMC *mc, NcmMSet *fiduc)
{
  if (mc->fiduc != NULL)
    ncm_mset_clear (&mc->fiduc);

  if (fiduc == NULL || fiduc == mc->fit->mset)
  {
    mc->fiduc = ncm_mset_dup (mc->fit->mset, mc->ser);
    ncm_serialize_clear_instances (mc->ser);
  }
  else
  {
    mc->fiduc = ncm_mset_ref (fiduc);
    g_assert (ncm_mset_cmp (mc->fit->mset, fiduc, FALSE));
  }
}

/**
 * ncm_fit_mc_set_rng:
 * @mc: a #NcmFitMC
 * @rng: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_mc_set_rng (NcmFitMC *mc, NcmRNG *rng)
{  
  if (mc->started)
    g_error ("ncm_fit_mc_set_rng: Cannot change the RNG object during a run, call ncm_fit_mc_end_run() first.");

  ncm_mset_catalog_set_rng (mc->mcat, rng);
}

void
_ncm_fit_mc_update (NcmFitMC *mc, NcmFit *fit)
{
  const guint part = 5;
  const guint step = (mc->n / part) == 0 ? 1 : (mc->n / part);

  ncm_mset_catalog_add_from_mset (mc->mcat, fit->mset, ncm_fit_state_get_m2lnL_curval (fit->fstate));
  ncm_timer_task_increment (mc->nt);

  switch (mc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = mc->nt->task_pos % step;
      gboolean log_timeout = FALSE;
      if ((mc->nt->pos_time - mc->nt->last_log_time) > 60.0)
        log_timeout = TRUE;
      if (log_timeout || (stepi == 0) || (mc->nt->task_pos == mc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (mc->mcat);
        /* ncm_timer_task_accumulate (mc->nt, acc); */
        ncm_timer_task_log_elapsed (mc->nt);
        ncm_timer_task_log_mean_time (mc->nt);
        ncm_timer_task_log_time_left (mc->nt);
        ncm_timer_task_log_end_datetime (mc->nt);
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      fit->mtype = mc->mtype;
      ncm_fit_log_start (fit);
      ncm_fit_log_end (fit);
      ncm_mset_catalog_log_current_stats (mc->mcat);
      /* ncm_timer_task_increment (mc->nt); */
      ncm_timer_task_log_elapsed (mc->nt);
      ncm_timer_task_log_mean_time (mc->nt);
      ncm_timer_task_log_time_left (mc->nt);
      ncm_timer_task_log_end_datetime (mc->nt);
      break;      
  }

  if ((mc->mcat->fmode != NCM_FIT_MC_MIN_FLUSH_INTERVAL) &&
      (ncm_timer_task_mean_time (mc->nt) < NCM_FIT_MC_MIN_FLUSH_INTERVAL))
  {
    ncm_mset_catalog_set_flush_mode (mc->mcat, NCM_MSET_CATALOG_FLUSH_TIMED);
    ncm_mset_catalog_set_flush_interval (mc->mcat, NCM_FIT_MC_MIN_FLUSH_INTERVAL);
  }
}

void
_ncm_fit_mc_resample_bstrap (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng)
{
  ncm_dataset_bootstrap_resample (dset, rng);
}

/**
 * ncm_fit_mc_start_run:
 * @mc: a #NcmFitMC
 * 
 * FIXME
 * 
 */
void 
ncm_fit_mc_start_run (NcmFitMC *mc)
{
  guint param_len = ncm_mset_param_len (mc->fit->mset);

  if (mc->started)
    g_error ("ncm_fit_mc_start_run: run already started, run ncm_fit_mc_end_run() first.");

  ncm_vector_clear (&mc->bf);
  mc->bf = ncm_vector_new (param_len);
  ncm_mset_param_get_vector (mc->fit->mset, mc->bf);

  switch (mc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Starting Monte Carlo...\n");
      ncm_dataset_log_info (mc->fit->lh->dset);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Fiducial model set:\n");
      ncm_mset_pretty_log (mc->fiduc);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitMC: Fitting model set:\n");
      ncm_mset_pretty_log (mc->fit->mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (mc->mcat->rng == NULL)
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

  ncm_mset_catalog_sync (mc->mcat, TRUE);
  if (mc->mcat->cur_id > mc->cur_sample_id)
  {
    ncm_fit_mc_intern_skip (mc, mc->mcat->cur_id - mc->cur_sample_id);
    g_assert_cmpint (mc->cur_sample_id, ==, mc->mcat->cur_id);
  }
  else if (mc->mcat->cur_id < mc->cur_sample_id)
    g_error ("ncm_fit_mc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
             mc->mcat->cur_id, mc->cur_sample_id);
  
}

/**
 * ncm_fit_mc_end_run:
 * @mc: a #NcmFitMC
 * 
 * FIXME
 * 
 */
void
ncm_fit_mc_end_run (NcmFitMC *mc)
{
  if (ncm_timer_task_is_running (mc->nt))
    ncm_timer_task_end (mc->nt);

  if (mc->mp != NULL)
  {
    ncm_memory_pool_free (mc->mp, TRUE);
    mc->mp = NULL;
  }

  ncm_mset_catalog_sync (mc->mcat, TRUE);
  
  mc->started = FALSE;
}

/**
 * ncm_fit_mc_reset:
 * @mc: a #NcmFitMC
 * 
 * FIXME
 * 
 */
void 
ncm_fit_mc_reset (NcmFitMC *mc)
{
  mc->n               = 0;
  mc->cur_sample_id   = -1;
  mc->write_index     = 0;
  mc->started         = FALSE;  
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
  mc->write_index = mc->cur_sample_id + 1;
}


/**
 * ncm_fit_mc_set_first_sample_id:
 * @mc: a #NcmFitMC
 * @first_sample_id: FIXME
 * 
 * FIXME
 *
 */
void 
ncm_fit_mc_set_first_sample_id (NcmFitMC *mc, gint first_sample_id)
{
  if (mc->mcat->first_id == first_sample_id)
    return;

  if (!mc->started)
    g_error ("ncm_fit_mc_set_first_sample_id: run not started, run ncm_fit_mc_start_run() first.");

  if (first_sample_id <= mc->cur_sample_id)
    g_error ("ncm_fit_mc_set_first_sample_id: cannot move first sample id backwards to: %d, catalog first id: %d, current sample id: %d.",
             first_sample_id, mc->mcat->first_id, mc->cur_sample_id);

  ncm_mset_catalog_set_first_id (mc->mcat, first_sample_id);
  ncm_fit_mc_intern_skip (mc, first_sample_id - mc->cur_sample_id - 1);
}

static void _ncm_fit_mc_run_single (NcmFitMC *mc);
static void _ncm_fit_mc_run_mt (NcmFitMC *mc);

/**
 * ncm_fit_mc_run:
 * @mc: a #NcmFitMC
 * @n: total number of realizations to run
 * 
 * Runs the Monte Carlo until it reaches the @n-th realization. Note that
 * if the first_id is non-zero it will run @n - first_id realizations.
 *
 */
void 
ncm_fit_mc_run (NcmFitMC *mc, guint n)
{
  if (!mc->started)
    g_error ("ncm_fit_mc_run: run not started, run ncm_fit_mc_start_run() first.");

  if (n <= (mc->cur_sample_id + 1))
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

  if (mc->nthreads <= 1)
    _ncm_fit_mc_run_single (mc);
  else
    _ncm_fit_mc_run_mt (mc);

  ncm_timer_task_pause (mc->nt);
}

static void 
_ncm_fit_mc_run_single (NcmFitMC *mc)
{
  guint i;

  for (i = 0; i < mc->n; i++)
  {
    ncm_mset_param_set_vector (mc->fit->mset, mc->bf);
    _ncm_fit_mc_resample (mc, mc->fit);
    ncm_fit_run (mc->fit, NCM_FIT_RUN_MSGS_NONE);

    _ncm_fit_mc_update (mc, mc->fit);
    mc->write_index++;
  }
}

static gpointer
_ncm_fit_mc_dup_fit (gpointer userdata)
{
  NcmFitMC *mc = NCM_FIT_MC (userdata);
  g_mutex_lock (mc->dup_fit);
  {
    NcmFit *fit = ncm_fit_dup (mc->fit, mc->ser);
    ncm_serialize_clear_instances (mc->ser);
    g_mutex_unlock (mc->dup_fit);
    return fit;
  }
}

static void 
_ncm_fit_mc_mt_eval_keep_order (glong i, glong f, gpointer data)
{
  _NCM_STATIC_MUTEX_DECL (update_lock);
  NcmFitMC *mc = NCM_FIT_MC (data);
  NcmFit **fit_ptr = ncm_memory_pool_get (mc->mp);
  NcmFit *fit = *fit_ptr;
  guint j;

  for (j = i; j < f; j++)
  {
    gint sample_index;

    ncm_mset_param_set_vector (fit->mset, mc->bf);

    g_mutex_lock (mc->resample_lock);
    sample_index = _ncm_fit_mc_resample (mc, fit);
    g_mutex_unlock (mc->resample_lock);

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

    _NCM_MUTEX_LOCK (&update_lock);    
    while (mc->write_index != sample_index)
      g_cond_wait (mc->write_cond, &update_lock);

    _ncm_fit_mc_update (mc, fit);
    mc->write_index++;
    g_cond_broadcast (mc->write_cond);

    _NCM_MUTEX_UNLOCK (&update_lock);
  }

  ncm_memory_pool_return (fit_ptr);
}


static void 
_ncm_fit_mc_mt_eval (glong i, glong f, gpointer data)
{
  _NCM_STATIC_MUTEX_DECL (update_lock);
  NcmFitMC *mc = NCM_FIT_MC (data);
  NcmFit **fit_ptr = ncm_memory_pool_get (mc->mp);
  NcmFit *fit = *fit_ptr;
  guint j;

  for (j = i; j < f; j++)
  {
    ncm_mset_param_set_vector (fit->mset, mc->bf);

    g_mutex_lock (mc->resample_lock);
    _ncm_fit_mc_resample (mc, fit);
    g_mutex_unlock (mc->resample_lock);

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

    _NCM_MUTEX_LOCK (&update_lock);    
    _ncm_fit_mc_update (mc, fit);
    mc->write_index++;
    _NCM_MUTEX_UNLOCK (&update_lock);
  }

  ncm_memory_pool_return (fit_ptr);
}


static void
_ncm_fit_mc_run_mt (NcmFitMC *mc)
{
  const guint nthreads = mc->n > mc->nthreads ? mc->nthreads : (mc->n - 1);

  if (nthreads == 0)
  {
    _ncm_fit_mc_run_single (mc);
    return;
  }
  
  if (mc->mp != NULL)
    ncm_memory_pool_free (mc->mp, TRUE);
  mc->mp = ncm_memory_pool_new (&_ncm_fit_mc_dup_fit, mc, 
                                (GDestroyNotify) &ncm_fit_free);

  /*
   * The line below added de main fit object to the pool, but can cause
   * several race conditions as it is used to make the copies for the other
   * threads. So, no.
   */
  //ncm_memory_pool_add (mc->mp, mc->fit);

  g_assert_cmpuint (mc->nthreads, >, 1);

  if (mc->keep_order)
    ncm_func_eval_threaded_loop_full (&_ncm_fit_mc_mt_eval_keep_order, 0, mc->n, mc);
  else
    ncm_func_eval_threaded_loop_full (&_ncm_fit_mc_mt_eval, 0, mc->n, mc);  
}

/**
 * ncm_fit_mc_run_lre:
 * @mc: a #NcmFitMC
 * @prerun: FIXME
 * @lre: FIXME
 *
 * FIXME
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
    gdouble n = ncm_mset_catalog_len (mc->mcat);
    gdouble m = n * lerror2 / lre2;
    guint runs = ((m - n) > 1000.0) ? ceil ((m - n) * 0.25) : ceil (m - n);

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
 * FIXME
 */
void
ncm_fit_mc_mean_covar (NcmFitMC *mc)
{
  ncm_mset_catalog_get_mean (mc->mcat, &mc->fit->fstate->fparams);
  ncm_mset_catalog_get_covar (mc->mcat, &mc->fit->fstate->covar);
  ncm_mset_fparams_set_vector (mc->mcat->mset, mc->fit->fstate->fparams);
  mc->fit->fstate->has_covar = TRUE;
}

/**
 * ncm_fit_mc_get_catalog:
 * @mc: a #NcmFitMC
 *
 * Gets the generated catalog of @mc.
 * 
 * Returns: (transfer full): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_mc_get_catalog (NcmFitMC *mc)
{
  return ncm_mset_catalog_ref (mc->mcat);
}

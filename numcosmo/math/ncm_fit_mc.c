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
 * SECTION:ncm_fit_mc
 * @title: NcmFitMC
 * @short_description: Monte Carlo analysis.
 *
 * This object implements a Monte Carlo analysis. This object is initialized
 * with a #NcmFit object and a #NcmMSet object. The #NcmFit object is used to
 * calculate the likelihood of the #NcmMSet object. The #NcmMSet object is
 * used to sample the parameter space.
 *
 * The #NcmFitMC object will resample the likelihood using the input #NcmMSet
 * object as a fiducial model. The resampling can be done in three ways:
 * #NCM_FIT_MC_RESAMPLE_FROM_MODEL, #NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX and
 * #NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX. The first option will resample the
 * likelihood from the input #NcmMSet object. The other two options will
 * resample the likelihood from the original likelihood data using bootstrap
 * resampling. The difference between the last two options is that the
 * #NCM_FIT_MC_RESAMPLE_BOOTSTRAP_NOMIX will resample each #NcmData separately
 * while the #NCM_FIT_MC_RESAMPLE_BOOTSTRAP_MIX will resample all #NcmData
 * together.
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
  PROP_NTHREADS,
  PROP_KEEP_ORDER,
  PROP_DATA_FILE,
};

struct _NcmFitMC
{
  /*< private >*/
  GObject parent_instance;
  NcmFitMCResample resample;
  NcmFit *fit;
  NcmMSet *fiduc;
  NcmMSetCatalog *mcat;
  NcmFitRunMsgs mtype;
  NcmFitMCResampleType rtype;
  NcmVector *bf;
  NcmTimer *nt;
  NcmSerialize *ser;
  guint nthreads;
  guint n;
  gboolean keep_order;
  NcmMemoryPool *mp;
  gint write_index;
  gint cur_sample_id;
  gint first_sample_id;
  gboolean started;
  GMutex dup_fit;
  GMutex resample_lock;
  GMutex update_lock;
  GCond write_cond;
};

G_DEFINE_TYPE (NcmFitMC, ncm_fit_mc, G_TYPE_OBJECT)

static void
ncm_fit_mc_init (NcmFitMC *mc)
{
  mc->fit           = NULL;
  mc->fiduc         = NULL;
  mc->mcat          = NULL;
  mc->mtype         = NCM_FIT_RUN_MSGS_NONE;
  mc->rtype         = NCM_FIT_MC_RESAMPLE_BOOTSTRAP_LEN;
  mc->bf            = NULL;
  mc->nt            = ncm_timer_new ();
  mc->ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  mc->nthreads      = 0;
  mc->n             = 0;
  mc->keep_order    = FALSE;
  mc->mp            = NULL;
  mc->cur_sample_id = -1; /* Represents that no samples were calculated yet. */
  mc->write_index   = 0;
  mc->started       = FALSE;
  g_mutex_init (&mc->dup_fit);
  g_mutex_init (&mc->resample_lock);
  g_mutex_init (&mc->update_lock);
  g_cond_init (&mc->write_cond);
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

  g_mutex_clear (&mc->dup_fit);
  g_mutex_clear (&mc->resample_lock);
  g_mutex_clear (&mc->update_lock);
  g_cond_clear (&mc->write_cond);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_mc_parent_class)->finalize (object);
}

static void
ncm_fit_mc_class_init (NcmFitMCClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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
  NcmMSet *mset = ncm_fit_peek_mset (fit);

  g_assert (mc->fit == NULL);
  mc->fit  = ncm_fit_ref (fit);
  mc->mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE,
                                   NCM_MSET_CATALOG_M2LNL_COLNAME, NCM_MSET_CATALOG_M2LNL_SYMBOL,
                                   NULL);
  ncm_mset_catalog_set_m2lnp_var (mc->mcat, 0);
}

/**
 * ncm_fit_mc_new:
 * @fit: a #NcmFit
 * @rtype: a #NcmFitMCResampleType
 * @mtype: a #NcmFitRunMsgs
 *
 * Creates a new #NcmFitMC object with the fit object @fit, the resample type
 * @rtype and the run messages type @mtype.
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
 * ncm_fit_mc_free:
 * @mc: a #NcmFitMC
 *
 * Decrement the reference count atomically by one. If the reference count
 * reaches zero, all memory allocated by the object is released and the object
 * is freed.
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
 * When this function is invoked, it first checks whether *@mc is not %NULL. If *@mc is
 * not NULL, the function reduces the reference count of *@mc by one. If the reference
 * count reaches zero, all memory allocated by the object is released, the object is
 * freed, and *@mc is set to %NULL.
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
 * @filename: a filename
 *
 * Sets the data file to be used by the #NcmMSetCatalog object of @mc.
 *
 */
void
ncm_fit_mc_set_data_file (NcmFitMC *mc, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (mc->mcat);

  if (mc->started && (cur_filename != NULL))
    g_error ("ncm_fit_mc_set_data_file: Cannot change data file during a run, call ncm_fit_mc_end_run() first.");

  if ((cur_filename != NULL) && (strcmp (cur_filename, filename) == 0))
    return;

  ncm_mset_catalog_set_file (mc->mcat, filename);

  if (mc->started)
  {
    const gint mcat_cur_id = ncm_mset_catalog_get_cur_id (mc->mcat);

    if (mcat_cur_id > mc->cur_sample_id)
    {
      ncm_fit_mc_intern_skip (mc, mcat_cur_id - mc->cur_sample_id);
      g_assert_cmpint (mc->cur_sample_id, ==, mcat_cur_id);
    }
    else if (mcat_cur_id < mc->cur_sample_id)
    {
      g_error ("ncm_fit_mc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].",
               mcat_cur_id, mc->cur_sample_id);
    }
  }
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

/**
 * ncm_fit_mc_set_rtype:
 * @mc: a #NcmFitMC
 * @rtype: a #NcmFitMCResampleType
 *
 * Sets the resample type of @mc to @rtype.
 *
 */
void
ncm_fit_mc_set_rtype (NcmFitMC *mc, NcmFitMCResampleType rtype)
{
  const GEnumValue *eval = ncm_cfg_enum_get_value (NCM_TYPE_FIT_MC_RESAMPLE_TYPE, rtype);
  NcmLikelihood *lh      = ncm_fit_peek_likelihood (mc->fit);
  NcmDataset *dset       = ncm_likelihood_peek_dataset (lh);

  if (mc->started)
    g_error ("ncm_fit_mc_set_rtype: Cannot change resample type during a run, call ncm_fit_mc_end_run() first.");

  mc->rtype = rtype;

  ncm_mset_catalog_set_run_type (mc->mcat, eval->value_nick);

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
 * ncm_fit_mc_set_nthreads:
 * @mc: a #NcmFitMC
 * @nthreads: number of threads
 *
 * Sets the number of threads to be used by @mc to @nthreads.
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
 * @keep_order: whether keep the catalog in order of sampling
 *
 * Sets whether keep the catalog in order of sampling. When performing
 * parallel runs, the catalog is not kept in order of sampling. This
 * is done to avoid the overhead of locking the catalog. If you want
 * to keep the catalog in order of sampling, set this property to
 * %TRUE.
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
 * @fiduc: a #NcmMSet
 *
 * Sets the fiducial model of @mc to @fiduc. If @fiduc is %NULL, the fiducial
 * model is set to the model set of the fit object of @mc. If @fiduc is not
 * %NULL, the fiducial model is set to @fiduc. If @fiduc is not %NULL, it must
 * be equal to the model set of the fit object of @mc.
 *
 * Note that in the end of an analysis, the #NcmMSet in @mc will be equal to
 * the last sampled model set.
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
 * Sets the RNG object of @mc to @rng.
 *
 */
void
ncm_fit_mc_set_rng (NcmFitMC *mc, NcmRNG *rng)
{
  if (mc->started)
    g_error ("ncm_fit_mc_set_rng: Cannot change the RNG object during a run, call ncm_fit_mc_end_run() first.");

  ncm_mset_catalog_set_rng (mc->mcat, rng);
}

/**
 * ncm_fit_mc_is_running:
 * @mc: a #NcmFitMC
 *
 * Checks whether a run is running, that is whether ncm_fit_mc_start_run()
 * was called and ncm_fit_mc_end_run() was not called yet.
 *
 * Returns: %TRUE if a run is running, %FALSE otherwise.
 */
gboolean
ncm_fit_mc_is_running (NcmFitMC *mc)
{
  return mc->started;
}

void
_ncm_fit_mc_update (NcmFitMC *mc, NcmFit *fit)
{
  NcmMSet *mset       = ncm_fit_peek_mset (fit);
  NcmFitState *fstate = ncm_fit_peek_state (fit);

  const guint part = 5;
  const guint step = (mc->n / part) == 0 ? 1 : (mc->n / part);

  ncm_mset_catalog_add_from_mset (mc->mcat, mset, ncm_fit_state_get_m2lnL_curval (fstate), NULL);
  ncm_timer_task_increment (mc->nt);

  switch (mc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi          = mc->nt->task_pos % step;
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
        ncm_timer_task_log_cur_datetime (mc->nt);
        ncm_timer_task_log_end_datetime (mc->nt);
      }

      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_fit_set_messages (fit, mc->mtype);
      ncm_fit_log_state (fit);
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

void
_ncm_fit_mc_resample_bstrap (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng)
{
  ncm_dataset_bootstrap_resample (dset, rng);
}

/**
 * ncm_fit_mc_start_run:
 * @mc: a #NcmFitMC
 *
 * Starts a Monte Carlo run. This function will start a Monte Carlo run
 * using the #NcmFit object and the #NcmMSet object of @mc.
 *
 * This method will not compute the any likelihood. It will only start
 * the run. To compute samples, call ncm_fit_mc_run() or ncm_fit_mc_run_lre()
 * after this function.
 *
 * To finish the run, call ncm_fit_mc_end_run().
 *
 */
void
ncm_fit_mc_start_run (NcmFitMC *mc)
{
  NcmLikelihood *lh      = ncm_fit_peek_likelihood (mc->fit);
  NcmMSet *mset          = ncm_fit_peek_mset (mc->fit);
  const guint param_len  = ncm_mset_total_len (mc->fiduc);
  const gint mcat_cur_id = ncm_mset_catalog_get_cur_id (mc->mcat);
  NcmRNG *mcat_rng       = ncm_mset_catalog_peek_rng (mc->mcat);
  NcmDataset *dset       = ncm_likelihood_peek_dataset (lh);

  if (mc->started)
    g_error ("ncm_fit_mc_start_run: run already started, run ncm_fit_mc_end_run() first.");

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

  ncm_mset_catalog_set_sync_mode (mc->mcat, NCM_MSET_CATALOG_SYNC_TIMED);
  ncm_mset_catalog_set_sync_interval (mc->mcat, NCM_FIT_MC_MIN_SYNC_INTERVAL);

  ncm_mset_catalog_sync (mc->mcat, TRUE);

  if (mcat_cur_id > mc->cur_sample_id)
  {
    ncm_fit_mc_intern_skip (mc, mcat_cur_id - mc->cur_sample_id);
    g_assert_cmpint (mc->cur_sample_id, ==, mcat_cur_id);
  }
  else if (mcat_cur_id < mc->cur_sample_id)
  {
    g_error ("ncm_fit_mc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].",
             mcat_cur_id, mc->cur_sample_id);
  }
}

/**
 * ncm_fit_mc_end_run:
 * @mc: a #NcmFitMC
 *
 * Ends a Monte Carlo run. This function will end a Monte Carlo run
 * using the #NcmFit object and the #NcmMSet object of @mc.
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
 * Resets the Monte Carlo run. This function will reset the Monte Carlo run
 * erase all samples and reset the catalog.
 *
 */
void
ncm_fit_mc_reset (NcmFitMC *mc)
{
  mc->n             = 0;
  mc->cur_sample_id = -1;
  mc->write_index   = 0;
  mc->started       = FALSE;
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

/**
 * ncm_fit_mc_set_first_sample_id:
 * @mc: a #NcmFitMC
 * @first_sample_id: first sample id
 *
 * Sets the first sample id of the Monte Carlo run to @first_sample_id.
 * This function will skip all samples until it reaches the @first_sample_id.
 *
 */
void
ncm_fit_mc_set_first_sample_id (NcmFitMC *mc, gint first_sample_id)
{
  const gint mcat_first_id = ncm_mset_catalog_get_first_id (mc->mcat);

  if (mcat_first_id == first_sample_id)
    return;

  if (!mc->started)
    g_error ("ncm_fit_mc_set_first_sample_id: run not started, run ncm_fit_mc_start_run() first.");

  if (first_sample_id <= mc->cur_sample_id)
    g_error ("ncm_fit_mc_set_first_sample_id: cannot move first sample id backwards to: %d, catalog first id: %d, current sample id: %d.",
             first_sample_id, mcat_first_id, mc->cur_sample_id);

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

  if (mc->nthreads <= 1)
    _ncm_fit_mc_run_single (mc);
  else
    _ncm_fit_mc_run_mt (mc);

  ncm_timer_task_pause (mc->nt);
}

static void
_ncm_fit_mc_run_single (NcmFitMC *mc)
{
  NcmMSet *mset = ncm_fit_peek_mset (mc->fit);
  guint i;

  for (i = 0; i < mc->n; i++)
  {
    ncm_mset_param_set_vector (mset, mc->bf);
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

  g_mutex_lock (&mc->dup_fit);
  {
    NcmFit *fit = ncm_fit_dup (mc->fit, mc->ser);

    ncm_serialize_reset (mc->ser, TRUE);
    g_mutex_unlock (&mc->dup_fit);

    return fit;
  }
}

static void
_ncm_fit_mc_mt_eval_keep_order (glong i, glong f, gpointer data)
{
  NcmFitMC *mc     = NCM_FIT_MC (data);
  NcmMSet *mset    = ncm_fit_peek_mset (mc->fit);
  NcmFit **fit_ptr = ncm_memory_pool_get (mc->mp);
  NcmFit *fit      = *fit_ptr;
  guint j;

  for (j = i; j < f; j++)
  {
    gint sample_index;

    ncm_mset_param_set_vector (mset, mc->bf);

    g_mutex_lock (&mc->resample_lock);
    sample_index = _ncm_fit_mc_resample (mc, fit);
    g_mutex_unlock (&mc->resample_lock);

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

    g_mutex_lock (&mc->update_lock);

    while (mc->write_index != sample_index)
      g_cond_wait (&mc->write_cond, &mc->update_lock);

    _ncm_fit_mc_update (mc, fit);
    mc->write_index++;
    g_cond_broadcast (&mc->write_cond);

    g_mutex_unlock (&mc->update_lock);
  }

  ncm_memory_pool_return (fit_ptr);
}

static void
_ncm_fit_mc_mt_eval (glong i, glong f, gpointer data)
{
  G_LOCK_DEFINE_STATIC (update_lock);

  NcmFitMC *mc     = NCM_FIT_MC (data);
  NcmMSet *mset    = ncm_fit_peek_mset (mc->fit);
  NcmFit **fit_ptr = ncm_memory_pool_get (mc->mp);
  NcmFit *fit      = *fit_ptr;
  guint j;

  for (j = i; j < f; j++)
  {
    ncm_mset_param_set_vector (mset, mc->bf);

    g_mutex_lock (&mc->resample_lock);
    _ncm_fit_mc_resample (mc, fit);
    g_mutex_unlock (&mc->resample_lock);

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

    G_LOCK (update_lock);
    _ncm_fit_mc_update (mc, fit);
    mc->write_index++;
    G_UNLOCK (update_lock);
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
                                (GDestroyNotify) & ncm_fit_free);

  /*
   * The line below added de main fit object to the pool, but can cause
   * several race conditions as it is used to make the copies for the other
   * threads. So, no.
   */
  /*ncm_memory_pool_add (mc->mp, mc->fit); */

  g_assert_cmpuint (mc->nthreads, >, 1);

  if (mc->keep_order)
    ncm_func_eval_threaded_loop_full (&_ncm_fit_mc_mt_eval_keep_order, 0, mc->n, mc);
  else
    ncm_func_eval_threaded_loop_full (&_ncm_fit_mc_mt_eval, 0, mc->n, mc);
}

/**
 * ncm_fit_mc_run_lre:
 * @mc: a #NcmFitMC
 * @prerun: number of pre-runs
 * @lre: largest relative error
 *
 * Runs the Monte Carlo until the largest relative error considering the
 * erros on the parameter means is less than @lre.
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
 * Computes the mean and covariance of the Monte Carlo run.
 * The mean and covariance are stored in the #NcmFit object of @mc.
 * The mean is stored in the #NcmFitState object of the #NcmFit object
 * and in the #NcmMSetCatalog object of @mc.
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
 * Gets the generated catalog of @mc.
 *
 * Returns: (transfer full): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_mc_get_catalog (NcmFitMC *mc)
{
  return ncm_mset_catalog_ref (mc->mcat);
}

/**
 * ncm_fit_mc_peek_catalog:
 * @mc: a #NcmFitMC
 *
 * Peeks the generated catalog of @mc.
 *
 * Returns: (transfer none): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_mc_peek_catalog (NcmFitMC *mc)
{
  return mc->mcat;
}


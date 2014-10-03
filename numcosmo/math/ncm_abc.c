/***************************************************************************
 *            ncm_abc.c
 *
 *  Tue September 30 15:46:48 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_abc.c
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
 * SECTION:ncm_abc
 * @title: Monte Carlo ABC Analysis
 * @short_description: Object implementing abstract Approximate Bayesian Computation (ABC)
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_abc.h"
#include "math/ncm_func_eval.h"

enum
{
  PROP_0,
  PROP_MSET,
  PROP_PRIOR,
  PROP_TKERN,
  PROP_DATASET,
  PROP_LEN
};

G_DEFINE_ABSTRACT_TYPE (NcmABC, ncm_abc, G_TYPE_OBJECT);

static void
ncm_abc_init (NcmABC *abc)
{
  abc->mcat          = NULL;
  abc->dset          = NULL;
  abc->mp            = NULL;
  abc->tkern         = NULL;
  abc->prior         = NULL;
  abc->nt            = ncm_timer_new ();
  abc->ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  abc->mtype         = NCM_FIT_RUN_MSGS_NONE;
  abc->theta         = NULL;
  abc->thetastar     = NULL;  
  abc->started       = FALSE;
  abc->cur_sample_id = -1; /* Represents that no samples were calculated yet. */
  abc->nthreads      = 0;
  abc->n             = 0;
}

static void
_ncm_abc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmABC *abc = NCM_ABC (object);
  g_return_if_fail (NCM_IS_ABC (object));

  switch (prop_id)
  {
    case PROP_MSET:
    {
      NcmMSet *mset = g_value_get_object (value);
      g_assert (abc->mcat == NULL);
      if (mset == NULL)
        g_error ("ncm_abc_new: mset (NcmMSet) cannot be NULL."); 
      abc->mcat = ncm_mset_catalog_new (mset, 2, "NcmABC:Distance", "NcmABC:Weight");
      break;
    }
    case PROP_PRIOR:
    {
      NcmMSetTransKern *prior = g_value_dup_object (value); 
      g_assert (abc->prior == NULL);
      if (prior == NULL)
        g_error ("ncm_abc_new: prior (NcmMSetTransKern) cannot be NULL.");       
      abc->prior = prior;
      break;
    }
    case PROP_TKERN:
      abc->tkern = g_value_dup_object (value); 
      break;
    case PROP_DATASET:
    {
      NcmDataset *dset = g_value_dup_object (value); 
      g_assert (abc->dset == NULL);
      if (dset == NULL)
        g_error ("ncm_abc_new: dset (NcmDataset) cannot be NULL.");       
      abc->dset = dset;
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_abc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmABC *abc = NCM_ABC (object);
  g_return_if_fail (NCM_IS_ABC (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, abc->mcat->mset);
      break;
    case PROP_PRIOR:
      g_value_set_object (value, abc->prior);
      break;
    case PROP_TKERN:
      g_value_set_object (value, abc->tkern);
      break;
    case PROP_DATASET:
      g_value_set_object (value, abc->dset);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_abc_dispose (GObject *object)
{
  NcmABC *abc = NCM_ABC (object);

  ncm_mset_catalog_clear (&abc->mcat);
  ncm_mset_trans_kern_clear (&abc->prior);
  ncm_mset_trans_kern_clear (&abc->tkern);
  ncm_dataset_clear (&abc->dset);
  ncm_timer_clear (&abc->nt);
  ncm_serialize_clear (&abc->ser);

  if (abc->mp != NULL)
  {
    ncm_memory_pool_free (abc->mp, TRUE);
    abc->mp = NULL;
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_abc_parent_class)->dispose (object);
}

static void
_ncm_abc_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_abc_parent_class)->finalize (object);
}

static void
ncm_abc_class_init (NcmABCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_abc_set_property;
  object_class->get_property = &_ncm_abc_get_property;
  object_class->dispose      = &_ncm_abc_dispose;
  object_class->finalize     = &_ncm_abc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "Model Set",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_PRIOR,
                                   g_param_spec_object ("prior",
                                                        NULL,
                                                        "Prior Sampler",
                                                        NCM_TYPE_MSET_TRANS_KERN,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_TKERN,
                                   g_param_spec_object ("trans-kernel",
                                                        NULL,
                                                        "Transition Kernel",
                                                        NCM_TYPE_MSET_TRANS_KERN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DATASET,
                                   g_param_spec_object ("data-set",
                                                        NULL,
                                                        "Dataset",
                                                        NCM_TYPE_DATASET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_abc_free:
 * @abc: a #NcmABC
 *
 * FIXME
 *
 */
void 
ncm_abc_free (NcmABC *abc)
{
  g_object_unref (abc);
}

/**
 * ncm_abc_clear:
 * @abc: a #NcmABC
 *
 * FIXME *
 */
void 
ncm_abc_clear (NcmABC **abc)
{
  g_clear_object (abc);
}

/**
 * ncm_abc_data_summary: (virtual data_summary)
 * @abc: a #NcmABC.
 *
 * Calculates the data summary and stores internally.
 * 
 * Returns: if the summary calculation was successful.
 */
gboolean
ncm_abc_data_summary (NcmABC *abc)
{
  return NCM_ABC_GET_CLASS (abc)->data_summary (abc);
}

/**
 * ncm_abc_mock_distance: (virtual mock_distance)
 * @abc: a #NcmABC.
 * @dset: a #NcmDataset.
 * @theta: a #NcmVector.
 * @thetastar: a #NcmVector.
 * @rng: a #NcmRNG.
 *
 * Calculates the distance of the new point given by @thetastar 
 * given the old point @theta.
 * 
 * Returns: the distance to the new point @thetastar.
 */
gdouble
ncm_abc_mock_distance (NcmABC *abc, NcmDataset *dset, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  return NCM_ABC_GET_CLASS (abc)->mock_distance (abc, dset, theta, thetastar, rng);
}

/**
 * ncm_abc_distance_prob: (virtual distance_prob)
 * @abc: a #NcmABC.
 * @distance: the distance.
 *
 * Calculates the probability of the distance been accepted. 
 * 
 * Returns: the probability of accepting the @distance.
 */
gdouble
ncm_abc_distance_prob (NcmABC *abc, gdouble distance)
{
  return NCM_ABC_GET_CLASS (abc)->distance_prob (abc, distance);
}

/**
 * ncm_abc_get_desc: (virtual get_desc)
 * @abc: a #NcmABC.
 *
 * Gets the description of the current ABC implementation. 
 * 
 * Returns: (transfer none): the description of the ABC implementation.
 */
const gchar *
ncm_abc_get_desc (NcmABC *abc)
{
  return NCM_ABC_GET_CLASS (abc)->get_desc (abc);
}

/**
 * ncm_abc_log_info: (virtual log_info)
 * @abc: a #NcmABC.
 *
 * Gets the informations about the current run of ABC. 
 * 
 * Returns: (transfer none): the informations about the current run.
 */
const gchar *
ncm_abc_log_info (NcmABC *abc)
{
  return NCM_ABC_GET_CLASS (abc)->log_info (abc);
}

/**
 * ncm_abc_set_mtype:
 * @abc: a #NcmFitMC
 * @mtype: FIXME
 *
 * FIXME
 *
 */
void 
ncm_abc_set_mtype (NcmABC *abc, NcmFitRunMsgs mtype)
{
  abc->mtype = mtype;
}

/**
 * ncm_abc_set_data_file:
 * @abc: a #NcmABC
 * @filename: a filename.
 *
 * FIXME
 *
 */
void 
ncm_abc_set_data_file (NcmABC *abc, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (abc->mcat);
  
  if (abc->started && cur_filename != NULL)
    g_error ("ncm_abc_set_data_file: Cannot change data file during a run, call ncm_abc_end_run() first.");    

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;

  ncm_mset_catalog_set_file (abc->mcat, filename);
  
  if (abc->started)
    g_assert_cmpint (abc->cur_sample_id, ==, abc->mcat->cur_id);
}

/**
 * ncm_abc_set_nthreads:
 * @abc: a #NcmABC
 * @nthreads: FIXME
 *
 * FIXME
 *
 */
void 
ncm_abc_set_nthreads (NcmABC *abc, guint nthreads)
{
  abc->nthreads = nthreads;
}

/**
 * ncm_abc_set_rng:
 * @abc: a #NcmABC
 * @rng: FIXME
 *
 * FIXME
 *
 */
void
ncm_abc_set_rng (NcmABC *abc, NcmRNG *rng)
{
  if (abc->started)
    g_error ("ncm_abc_set_rng: Cannot change the RNG object during a run, call ncm_abc_end_run() first.");

  ncm_mset_catalog_set_rng (abc->mcat, rng);
}

void
_ncm_abc_update (NcmABC *abc, NcmMSet *mset, gdouble dist)
{
  const guint part = 5;
  const guint step = (abc->n / part) == 0 ? 1 : (abc->n / part);

  ncm_mset_catalog_add_from_mset (abc->mcat, mset, dist, 1.0);
  ncm_timer_task_increment (abc->nt);

  switch (abc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = abc->nt->task_pos % step;
      if ((stepi == 0) || (abc->nt->task_pos == abc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (abc->mcat);
        /* ncm_timer_task_accumulate (abc->nt, acc); */
        ncm_timer_task_log_elapsed (abc->nt);
        ncm_timer_task_log_mean_time (abc->nt);
        ncm_timer_task_log_time_left (abc->nt);
        ncm_timer_task_log_end_datetime (abc->nt);
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_mset_catalog_log_current_stats (abc->mcat);
      /* ncm_timer_task_increment (abc->nt); */
      ncm_timer_task_log_elapsed (abc->nt);
      ncm_timer_task_log_mean_time (abc->nt);
      ncm_timer_task_log_time_left (abc->nt);
      ncm_timer_task_log_end_datetime (abc->nt);
      break;      
  }

  if ((abc->mcat->fmode != NCM_ABC_MIN_FLUSH_INTERVAL) &&
      (ncm_timer_task_mean_time (abc->nt) < NCM_ABC_MIN_FLUSH_INTERVAL))
  {
    ncm_mset_catalog_set_flush_mode (abc->mcat, NCM_MSET_CATALOG_FLUSH_TIMED);
    ncm_mset_catalog_set_flush_interval (abc->mcat, NCM_ABC_MIN_FLUSH_INTERVAL);
  }
}

static void ncm_abc_intern_skip (NcmABC *abc, guint n);

/**
 * ncm_abc_start_run:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 */
void 
ncm_abc_start_run (NcmABC *abc)
{
  if (abc->started)
    g_error ("ncm_abc_start_run: run already started, run ncm_abc_end_run() first.");

  switch (abc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Starting ABC (%s)...\n", ncm_abc_get_desc (abc));
      ncm_message ("%s", ncm_abc_log_info (abc));
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Model set:\n");
      ncm_mset_pretty_log (abc->mcat->mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (abc->mcat->rng == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);
    ncm_rng_set_random_seed (rng, FALSE);
    ncm_abc_set_rng (abc, rng);
    if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmABC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));
    ncm_rng_free (rng);
  }

  abc->started = TRUE;

  {
    guint fparam_len = ncm_mset_fparam_len (abc->mcat->mset);
    if (abc->theta != NULL)
    {
      ncm_vector_free (abc->theta);
      ncm_vector_free (abc->thetastar);
    }
    abc->theta = ncm_vector_new (fparam_len);
    abc->thetastar = ncm_vector_new (fparam_len);
  }
  
  ncm_mset_catalog_sync (abc->mcat, TRUE);
  if (abc->mcat->cur_id > abc->cur_sample_id)
  {
    ncm_abc_intern_skip (abc, abc->mcat->cur_id - abc->cur_sample_id);
    g_assert_cmpint (abc->cur_sample_id, ==, abc->mcat->cur_id);
  }
  else if (abc->mcat->cur_id < abc->cur_sample_id)
    g_error ("ncm_abc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
             abc->mcat->cur_id, abc->cur_sample_id);
  
  {
    NcmVector *cur_row = NULL;
    
    cur_row = ncm_mset_catalog_peek_current_row (abc->mcat);
    if (cur_row != NULL)
    {
      ncm_mset_fparams_set_vector_offset (abc->mcat->mset, cur_row, 1);
    }
  }

  if (!ncm_abc_data_summary (abc))
    g_error ("ncm_abc_start_run: error calculating summary data.");
}

/**
 * ncm_abc_end_run:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 */
void
ncm_abc_end_run (NcmABC *abc)
{
  if (ncm_timer_task_is_running (abc->nt))
    ncm_timer_task_end (abc->nt);

  ncm_mset_catalog_sync (abc->mcat, TRUE);
  
  abc->started = FALSE;
}

/**
 * ncm_abc_reset:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 */
void 
ncm_abc_reset (NcmABC *abc)
{
  abc->n               = 0;
  abc->cur_sample_id   = -1;
  abc->started         = FALSE;  
  ncm_mset_catalog_reset (abc->mcat);
}

static void 
ncm_abc_intern_skip (NcmABC *abc, guint n)
{
  if (n == 0)
    return;

  switch (abc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Skipping %u tries, will start at %u-th try.\n", n, abc->cur_sample_id + n + 1 + 1);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  abc->cur_sample_id += n;
}

/**
 * ncm_abc_set_first_sample_id:
 * @abc: a #NcmABC
 * @first_sample_id: FIXME
 * 
 * FIXME
 *
 */
void 
ncm_abc_set_first_sample_id (NcmABC *abc, gint first_sample_id)
{
  if (abc->mcat->first_id == first_sample_id)
    return;

  if (!abc->started)
    g_error ("ncm_abc_set_first_sample_id: run not started, run ncm_abc_start_run () first.");

  if (first_sample_id <= abc->cur_sample_id)
    g_error ("ncm_abc_set_first_sample_id: cannot move first sample id backwards to: %d, catalog first id: %d, current sample id: %d.",
             first_sample_id, abc->mcat->first_id, abc->cur_sample_id);

  ncm_mset_catalog_set_first_id (abc->mcat, first_sample_id);
}

static void _ncm_abc_run_single (NcmABC *abc);
static void _ncm_abc_run_mt (NcmABC *abc);

/**
 * ncm_abc_run:
 * @abc: a #NcmABC
 * @n: total number of realizations to run
 * 
 * Runs the montecarlo until it reaches the @n-th realization. Note that
 * if the first_id is non-zero it will run @n - first_id realizations.
 *
 */
void 
ncm_abc_run (NcmABC *abc, guint n)
{
  if (!abc->started)
    g_error ("ncm_abc_run: run not started, run ncm_abc_start_run() first.");

  if (n <= (abc->cur_sample_id + 1))
  {
    if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Nothing to do, current ABC particle number is %d\n", abc->cur_sample_id + 1);
    }
    return;
  }
  
  abc->n = n - (abc->cur_sample_id + 1);
  
  switch (abc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Calculating [%06d] ABC particles [%s]\n", abc->n, ncm_abc_get_desc (abc));
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }
  
  if (ncm_timer_task_is_running (abc->nt))
  {
    ncm_timer_task_add_tasks (abc->nt, abc->n);
    ncm_timer_task_continue (abc->nt);
  }
  else
  {
    ncm_timer_task_start (abc->nt, abc->n);
    ncm_timer_set_name (abc->nt, "NcmABC");
  }
  if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (abc->nt);

  if (abc->nthreads <= 1)
    _ncm_abc_run_single (abc);
  else
    _ncm_abc_run_mt (abc);

  ncm_timer_task_pause (abc->nt);
}

static void 
_ncm_abc_run_single (NcmABC *abc)
{
  guint i = 0;
  
  for (i = 0; i < abc->n;)
  {
    gdouble dist = 0.0, prob = 0.0;
    ncm_mset_trans_kern_prior_sample (abc->prior, abc->thetastar, abc->mcat->rng);

    ncm_mset_fparams_set_vector (abc->mcat->mset, abc->thetastar);
    ncm_dataset_resample (abc->dset, abc->mcat->mset, abc->mcat->rng);
    
    dist = ncm_abc_mock_distance (abc, abc->dset, abc->theta, abc->thetastar, abc->mcat->rng);
    prob = ncm_abc_distance_prob (abc, dist);
    
    if (prob == 1.0 || (prob != 0.0 && gsl_rng_uniform (abc->mcat->rng->r) < prob))
    {
      _ncm_abc_update (abc, abc->mcat->mset, dist);
      i++;
    }
  }
}

typedef struct _NcmABCThread
{
  NcmMSet *mset;
  NcmDataset *dset;
  NcmVector *thetastar;
  NcmRNG *rng;
} NcmABCThread;

static gpointer
_ncm_abc_dup_thread (gpointer userdata)
{
  _NCM_STATIC_MUTEX_DECL (dup_thread);
  NcmABC *abc = NCM_ABC (userdata);
  NcmABCThread *abct = g_new (NcmABCThread, 1);
  
  _NCM_MUTEX_LOCK (&dup_thread);
  {
    abct->mset      = ncm_mset_dup (abc->mcat->mset, abc->ser);
    abct->dset      = ncm_dataset_dup (abc->dset, abc->ser);

    abct->thetastar = ncm_vector_dup (abc->thetastar);
    abct->rng       = ncm_rng_new (NULL);
    ncm_rng_set_seed (abct->rng, gsl_rng_get (abc->mcat->rng->r));

    ncm_serialize_reset (abc->ser);
    _NCM_MUTEX_UNLOCK (&dup_thread);
    return abct;
  }
}

static void
_ncm_abc_free_thread (gpointer data)
{
  NcmABCThread *abct = (NcmABCThread *) data;
  
  ncm_mset_clear (&abct->mset);
  ncm_dataset_clear (&abct->dset);
  ncm_vector_clear (&abct->thetastar);
  ncm_rng_clear (&abct->rng);
  g_free (abct);  
}

static void 
_ncm_abc_thread_eval (glong i, glong f, gpointer data)
{
  _NCM_STATIC_MUTEX_DECL (update_lock);
  NcmABC *abc = NCM_ABC (data);
  NcmABCThread **abct_ptr = ncm_memory_pool_get (abc->mp);
  NcmABCThread *abct = *abct_ptr;
  guint j;

  for (j = i; j < f;)
  {
    gdouble dist, prob;
    ncm_mset_trans_kern_prior_sample (abc->prior, abct->thetastar, abct->rng);
    ncm_mset_fparams_set_vector (abct->mset, abct->thetastar);
    ncm_dataset_resample (abct->dset, abct->mset, abct->rng);
    
    dist = ncm_abc_mock_distance (abc, abct->dset, abct->thetastar, abct->thetastar, abct->rng);
    prob = ncm_abc_distance_prob (abc, dist);
    
    if (prob == 1.0 || (prob != 0.0 && gsl_rng_uniform (abct->rng->r) < prob))
    {
      _NCM_MUTEX_LOCK (&update_lock);    
      _ncm_abc_update (abc, abct->mset, dist);
      abc->cur_sample_id++;
      j++;
      _NCM_MUTEX_UNLOCK (&update_lock);
    }
  }

  ncm_memory_pool_return (abct_ptr);
}

static void
_ncm_abc_run_mt (NcmABC *abc)
{
  const guint nthreads = abc->n > abc->nthreads ? abc->nthreads : (abc->n - 1);

  if (nthreads == 0)
  {
    _ncm_abc_run_single (abc);
    return;
  }
  
  if (abc->mp != NULL)
    ncm_memory_pool_free (abc->mp, TRUE);
  abc->mp = ncm_memory_pool_new (&_ncm_abc_dup_thread, abc, 
                                 (GDestroyNotify) &_ncm_abc_free_thread);
  
  g_assert_cmpuint (abc->nthreads, >, 1);

  ncm_func_eval_threaded_loop_full (&_ncm_abc_thread_eval, 0, abc->n, abc);
}

/**
 * ncm_abc_run_lre:
 * @abc: a #NcmABC
 * @prerun: FIXME
 * @lre: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_abc_run_lre (NcmABC *abc, guint prerun, gdouble lre)
{
  gdouble lerror;
  const gdouble lre2 = lre * lre;

  g_assert_cmpfloat (lre, >, 0.0);
  /* g_assert_cmpfloat (lre, <, 1.0); */

  if (prerun == 0)
  {
    guint fparam_len = ncm_mset_fparam_len (abc->mcat->mset);
    prerun = fparam_len * 100;
  }

  if (ncm_mset_catalog_len (abc->mcat) < prerun)
  {
    guint prerun_left = prerun - ncm_mset_catalog_len (abc->mcat);
    if (abc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmABC: Running first %u pre-runs...\n", prerun_left);
    ncm_abc_run (abc, prerun);
  }

  lerror = ncm_mset_catalog_largest_error (abc->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n = ncm_mset_catalog_len (abc->mcat);
    gdouble m = n * lerror2 / lre2;
    guint runs = ((m - n) > 1000.0) ? ceil ((m - n) * 0.25) : ceil (m - n);

    if (abc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmABC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmABC: Running more %u runs...\n", runs);
    }
    ncm_abc_run (abc, abc->cur_sample_id + runs + 1);
    lerror = ncm_mset_catalog_largest_error (abc->mcat);
  }

  if (abc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    g_message ("# NcmABC: Largest relative error %e attained: %e\n", lre, lerror);
}

/**
 * ncm_abc_mean_covar:
 * @abc: a #NcmABC
 * @fit: a #NcmFit
 *
 * FIXME
 */
void
ncm_abc_mean_covar (NcmABC *abc, NcmFit *fit)
{
  ncm_mset_catalog_get_mean (abc->mcat, &fit->fstate->fparams);
  ncm_mset_catalog_get_covar (abc->mcat, &fit->fstate->covar);
  ncm_mset_fparams_set_vector (abc->mcat->mset, fit->fstate->fparams);
  fit->fstate->has_covar = TRUE;
}

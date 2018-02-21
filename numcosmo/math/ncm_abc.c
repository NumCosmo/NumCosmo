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
 * @title: NcmABC
 * @short_description: Abstract class for Approximate Bayesian Computation (ABC).
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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_statistics_double.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_MSET,
  PROP_PRIOR,
  PROP_TKERN,
  PROP_DATASET,
  PROP_EPSILON,
  PROP_NPARTICLES,
  PROP_LEN
};

G_DEFINE_ABSTRACT_TYPE (NcmABC, ncm_abc, G_TYPE_OBJECT);

static void
ncm_abc_init (NcmABC *abc)
{
  abc->mcat          = NULL;
  abc->dset          = NULL;
  abc->dset_mock     = NULL;
  abc->mp            = NULL;
  abc->tkern         = NULL;
  abc->prior         = NULL;
  abc->nt            = ncm_timer_new ();
  abc->ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  abc->mtype         = NCM_FIT_RUN_MSGS_NONE;
  abc->theta         = NULL;
  abc->thetastar     = NULL;  
  abc->covar         = NULL;
  abc->weights       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  abc->weights_tm1   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  abc->dists         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  abc->dists_sorted  = FALSE;
  abc->epsilon       = 0.0;
  abc->depsilon      = 0.0;
  abc->wran          = NULL;
  abc->started       = FALSE;
  abc->started_up    = FALSE;
  abc->cur_sample_id = -1; /* Represents that no samples were calculated yet. */
  abc->nthreads      = 0;
  abc->nparticles    = 0;
  abc->n             = 0;
  abc->ntotal        = 0;
  abc->naccepted     = 0;
  abc->nupdates      = 0;
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
      abc->mcat = ncm_mset_catalog_new (mset, 1, 1, TRUE, "NcmABC:Distance", "d", NULL);
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
      ncm_abc_set_trans_kern (abc, g_value_get_object (value));
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
    case PROP_EPSILON:
      abc->epsilon  = g_value_get_double (value);
      abc->depsilon = abc->epsilon;
      break;
    case PROP_NPARTICLES:
      g_assert (!abc->started);
      g_assert_not_reached ();
      abc->nparticles = g_value_get_uint (value);
      break;
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
      g_value_set_object (value, ncm_mset_catalog_peek_mset (abc->mcat));
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
    case PROP_EPSILON:
      g_value_set_double (value, abc->epsilon);
      break;
    case PROP_NPARTICLES:
      g_value_set_uint (value, abc->nparticles);
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
  ncm_dataset_clear (&abc->dset_mock);
  ncm_timer_clear (&abc->nt);
  ncm_serialize_clear (&abc->ser);

  ncm_vector_clear (&abc->theta);
  ncm_vector_clear (&abc->thetastar);
  ncm_matrix_clear (&abc->covar);

  g_clear_pointer (&abc->weights, g_array_unref);
  g_clear_pointer (&abc->weights_tm1, g_array_unref);
  g_clear_pointer (&abc->dists, g_array_unref);
  g_clear_pointer (&abc->wran, gsl_ran_discrete_free);

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
  g_object_class_install_property (object_class,
                                   PROP_EPSILON,
                                   g_param_spec_double ("epsilon",
                                                        NULL,
                                                        "epsilon",
                                                        0.0, G_MAXDOUBLE, 1.0e20,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NPARTICLES,
                                   g_param_spec_uint ("nparticles",
                                                      NULL,
                                                      "Number of particles",
                                                      0, G_MAXUINT, 100,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * FIXME 
 */
void 
ncm_abc_clear (NcmABC **abc)
{
  g_clear_object (abc);
}

/**
 * ncm_abc_data_summary: (virtual data_summary)
 * @abc: a #NcmABC
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
 * @abc: a #NcmABC
 * @dset: a #NcmDataset
 * @theta: a #NcmVector
 * @thetastar: a #NcmVector
 * @rng: a #NcmRNG
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
 * @abc: a #NcmABC
 * @distance: the distance
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
 * ncm_abc_update_tkern: (virtual update_tkern)
 * @abc: a #NcmABC
 *
 * Updates the transition kernel present in @abc->tkern. 
 * 
 */
void
ncm_abc_update_tkern (NcmABC *abc)
{
  NCM_ABC_GET_CLASS (abc)->update_tkern (abc);
}


/**
 * ncm_abc_get_desc: (virtual get_desc)
 * @abc: a #NcmABC
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
 * @abc: a #NcmABC
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
 * @mtype: a #NcmFitRunMsgs
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
 * @filename: a filename
 *
 * FIXME
 *
 */
void 
ncm_abc_set_data_file (NcmABC *abc, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (abc->mcat);
  
  if (abc->cur_sample_id >= 0 && cur_filename != NULL)
    g_error ("ncm_abc_set_data_file: Cannot change data file during a run, call ncm_abc_end_run() first.");    

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;

  ncm_mset_catalog_set_file (abc->mcat, filename);
  
  if (abc->started)
    g_assert_cmpint (abc->cur_sample_id, ==, ncm_mset_catalog_get_cur_id (abc->mcat));
}

/**
 * ncm_abc_set_nthreads:
 * @abc: a #NcmABC
 * @nthreads: number of threads
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
 * @rng: a #NcmRNG
 *
 * FIXME
 *
 */
void
ncm_abc_set_rng (NcmABC *abc, NcmRNG *rng)
{
  if (abc->cur_sample_id >= 0)
    g_error ("ncm_abc_set_rng: Cannot change the RNG object in a non-empty catalog.");

  ncm_mset_catalog_set_rng (abc->mcat, rng);
}

/**
 * ncm_abc_set_trans_kern:
 * @abc: a #NcmABC
 * @tkern: a #NcmMSetTransKern
 * 
 * FIXME
 *
 */
void 
ncm_abc_set_trans_kern (NcmABC *abc, NcmMSetTransKern *tkern)
{
  if (abc->tkern != NULL)
    g_warning ("Transition kernel already set, replacing.");
  ncm_mset_trans_kern_clear (&abc->tkern);
  abc->tkern = ncm_mset_trans_kern_ref (tkern);
}

static gint 
_compare (gconstpointer a, gconstpointer b)
{
  if (*(gdouble*)a > *(gdouble*)b) 
    return 1;
  else if (*(gdouble*)a < *(gdouble*)b) 
    return -1;
  else 
    return 0;  
}

/**
 * ncm_abc_get_dist_quantile:
 * @abc: a #NcmABC
 * @p: FIXME
 * 
 * FIXME
 *
 */
gdouble 
ncm_abc_get_dist_quantile (NcmABC *abc, gdouble p)
{
  if (abc->started)
    g_error ("ncm_abc_get_dist_quantile: Cannot get quantiles during a run, call ncm_abc_end_run() first.");
  
  if (abc->dists->len < 1)
    g_error ("ncm_abc_get_dist_quantile: Cannot get quantiles, no particles calculated.");
  
  if (!abc->dists_sorted)
  {
    g_array_sort (abc->dists, &_compare);
    abc->dists_sorted = TRUE;
  }
  
  return gsl_stats_quantile_from_sorted_data ((gdouble *)abc->dists->data, 1, abc->dists->len, p);
}

/**
 * ncm_abc_get_accept_ratio:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_abc_get_accept_ratio (NcmABC *abc)
{
  return abc->naccepted * 1.0 / (abc->ntotal * 1.0);
}

/**
 * ncm_abc_update_epsilon:
 * @abc: a #NcmABC
 * @epsilon: new epsilon.
 * 
 * FIXME
 * 
 */
void 
ncm_abc_update_epsilon (NcmABC *abc, gdouble epsilon)
{
  if (epsilon >= abc->epsilon)
    g_warning ("ncm_abc_update_epsilon: increasing epsilon.");

  if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
  {
    guint i;
    g_message ("# NcmABC: ");
    for (i = 0; i < 8; i++)
    {
      gdouble p = (15.0 * i) / 100.0;
      p = p > 1.0 ? 1.0 : p;
      g_message ("[%02.0f%% %-8.2g] ", 100.0 * p, ncm_abc_get_dist_quantile (abc, p));
    }
    g_message ("\n");
    g_message ("# NcmABC: epsilon_t        = %g.\n", 
               abc->epsilon);
    g_message ("# NcmABC: epsilon_t+1      = %g.\n", 
               epsilon);
    g_message ("# NcmABC: depsilon/epsilon = %04.2f%%.\n", 100.0 * (abc->epsilon - epsilon) / abc->epsilon);
  }

  abc->depsilon = abc->epsilon - epsilon;
  abc->epsilon  = epsilon;
}

/**
 * ncm_abc_get_epsilon:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_abc_get_epsilon (NcmABC *abc)
{
  return abc->epsilon;
}

/**
 * ncm_abc_get_depsilon:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_abc_get_depsilon (NcmABC *abc)
{
  return abc->depsilon;
}

void
_ncm_abc_update (NcmABC *abc, NcmMSet *mset, gdouble dist, gdouble weight)
{
  const guint part = 5;
  const guint step = (abc->n / part) == 0 ? 1 : (abc->n / part);
  const guint pindex = abc->cur_sample_id % abc->nparticles;

  ncm_mset_catalog_add_from_mset (abc->mcat, mset, dist, weight, NULL);
  ncm_timer_task_increment (abc->nt);
  g_array_index (abc->weights, gdouble, pindex) = weight;
  g_array_index (abc->dists,   gdouble, pindex) = dist;
  abc->dists_sorted = FALSE;
  
  switch (abc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = (abc->cur_sample_id + 1) % step;
      gboolean log_timeout = FALSE;    
      if ((abc->nt->pos_time - abc->nt->last_log_time) > 60.0)
        log_timeout = TRUE;
    
      if (log_timeout || (stepi == 0) || (abc->nt->task_pos == abc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (abc->mcat);
        g_message ("# NcmABC:acceptance ratio %7.4f%%.\n", ncm_abc_get_accept_ratio (abc) * 100.0);
        /* ncm_timer_task_accumulate (abc->nt, acc); */
        ncm_timer_task_log_elapsed (abc->nt);
        ncm_timer_task_log_mean_time (abc->nt);
        ncm_timer_task_log_time_left (abc->nt);
        ncm_timer_task_log_end_datetime (abc->nt);
        ncm_cfg_logfile_flush_now ();
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_mset_catalog_log_current_stats (abc->mcat);
      g_message ("# NcmABC:acceptance ratio %7.4f%%.\n", ncm_abc_get_accept_ratio (abc) * 100.0);
      /* ncm_timer_task_increment (abc->nt); */
      ncm_timer_task_log_elapsed (abc->nt);
      ncm_timer_task_log_mean_time (abc->nt);
      ncm_timer_task_log_time_left (abc->nt);
      ncm_timer_task_log_end_datetime (abc->nt);
      ncm_cfg_logfile_flush_now ();
      break;      
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
  NcmMSet *mset     = ncm_mset_catalog_peek_mset (abc->mcat);
  const gint cur_id = ncm_mset_catalog_get_cur_id (abc->mcat);
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
      ncm_mset_pretty_log (mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_mset_catalog_peek_rng (abc->mcat) == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);
    
    ncm_rng_set_random_seed (rng, FALSE);
    ncm_abc_set_rng (abc, rng);

    if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmABC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));

    ncm_rng_free (rng);
  }

  abc->dists_sorted = FALSE;
  abc->started = TRUE;

  {
    guint fparam_len = ncm_mset_fparam_len (mset);
    if (abc->theta != NULL)
    {
      ncm_vector_free (abc->theta);
      ncm_vector_free (abc->thetastar);
      ncm_matrix_free (abc->covar);
    }
    abc->theta     = ncm_vector_new (fparam_len);
    abc->thetastar = ncm_vector_new (fparam_len);
    abc->covar     = ncm_matrix_new (fparam_len, fparam_len);
  }

  ncm_mset_catalog_set_sync_mode (abc->mcat, NCM_MSET_CATALOG_SYNC_TIMED);
  ncm_mset_catalog_set_sync_interval (abc->mcat, NCM_ABC_MIN_SYNC_INTERVAL);

  ncm_mset_catalog_sync (abc->mcat, TRUE);
  if (cur_id > abc->cur_sample_id)
  {
    ncm_abc_intern_skip (abc, cur_id - abc->cur_sample_id);
    g_assert_cmpint (abc->cur_sample_id, ==, cur_id);
  }
  else if (cur_id < abc->cur_sample_id)
    g_error ("ncm_abc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
             cur_id, abc->cur_sample_id);
  
  {
    NcmVector *cur_row = NULL;
    
    cur_row = ncm_mset_catalog_peek_current_row (abc->mcat);
    if (cur_row != NULL)
    {
      ncm_mset_fparams_set_vector_offset (mset, cur_row, 2);
    }
  }

  if (!ncm_abc_data_summary (abc))
    g_error ("ncm_abc_start_run: error calculating summary data.");

  ncm_dataset_clear (&abc->dset_mock);
  abc->dset_mock = ncm_dataset_dup (abc->dset, abc->ser);
  ncm_serialize_reset (abc->ser, TRUE);

  abc->ntotal = 0;
  abc->naccepted = 0;
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
  NcmMSet *mset = ncm_mset_catalog_peek_mset (abc->mcat);
  gdouble WT    = 0.0;
  guint i;

  g_assert (abc->started);
  
  if (ncm_timer_task_is_running (abc->nt))
    ncm_timer_task_end (abc->nt);

  ncm_mset_catalog_sync (abc->mcat, TRUE);

  for (i = 0; i < abc->nparticles; i++)
    WT += g_array_index (abc->weights, gdouble, i);
  
  for (i = 0; i < abc->nparticles; i++)
    g_array_index (abc->weights, gdouble, i) = g_array_index (abc->weights, gdouble, i) / WT;

  switch (abc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    case NCM_FIT_RUN_MSGS_FULL:
      g_message ("# NcmABC:Current covariance matrix:\n");
      ncm_mset_catalog_get_mean (abc->mcat, &abc->theta);
      ncm_mset_fparams_set_vector (mset, abc->theta);
      ncm_mset_catalog_get_covar (abc->mcat, &abc->covar);
      ncm_mset_fparams_log_covar (mset, abc->covar);
      break;
    default:
      break;
  }
  
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
  abc->ntotal          = 0;
  abc->naccepted       = 0;
  abc->nupdates        = 0;
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

static void _ncm_abc_run_single (NcmABC *abc);
static void _ncm_abc_run_mt (NcmABC *abc);

/**
 * ncm_abc_run:
 * @abc: a #NcmABC
 * @nparticles: total number of particles to generate
 * 
 * Generates particlea until @n particles are accepted.
 *
 */
void 
ncm_abc_run (NcmABC *abc, guint nparticles)
{
  if (!abc->started)
    g_error ("ncm_abc_run: run not started, run ncm_abc_start_run() first.");
  
  if (nparticles <= (abc->cur_sample_id + 1))
  {
    if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Nothing to do, current ABC particle number is already %d of %u\n", abc->cur_sample_id + 1, nparticles);
    }
    return;
  }
  
  abc->n = nparticles - (abc->cur_sample_id + 1);

  if (abc->n > 0 && abc->nupdates > 0)
    g_error ("ncm_abc_run: cannot generate new particles when time t = %u != 0.", abc->nupdates);
  
  abc->nparticles = nparticles;
  g_array_set_size (abc->weights, nparticles);
  g_array_set_size (abc->dists, nparticles);
  
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
  
  g_assert_cmpuint (abc->nparticles, ==, abc->cur_sample_id + 1);
}

static void 
_ncm_abc_run_single (NcmABC *abc)
{
  NcmMSet *mset = ncm_mset_catalog_peek_mset (abc->mcat);
  NcmRNG *rng   = ncm_mset_catalog_peek_rng (abc->mcat);

  guint i = 0;
  
  for (i = 0; i < abc->n;)
  {
    gdouble dist = 0.0, prob = 0.0;
    ncm_mset_trans_kern_prior_sample (abc->prior, abc->thetastar, rng);

    ncm_mset_fparams_set_vector (mset, abc->thetastar);
    ncm_dataset_resample (abc->dset_mock, mset, rng);
    
    dist = ncm_abc_mock_distance (abc, abc->dset_mock, abc->theta, abc->thetastar, rng);
    prob = ncm_abc_distance_prob (abc, dist);

    abc->ntotal++;
    
    if (prob == 1.0 || (prob != 0.0 && gsl_rng_uniform (rng->r) < prob))
    {
      abc->cur_sample_id++;
      abc->naccepted++;
      _ncm_abc_update (abc, mset, dist, 1.0);
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
  G_LOCK_DEFINE_STATIC (dup_thread);
  NcmABC *abc        = NCM_ABC (userdata);
  NcmABCThread *abct = g_new (NcmABCThread, 1);
  NcmMSet *mset      = ncm_mset_catalog_peek_mset (abc->mcat);
  NcmRNG *rng        = ncm_mset_catalog_peek_rng (abc->mcat);

  G_LOCK (dup_thread);
  {
    abct->mset      = ncm_mset_dup (mset, abc->ser);
    abct->dset      = ncm_dataset_dup (abc->dset, abc->ser);
    abct->thetastar = ncm_vector_dup (abc->thetastar);
    abct->rng       = ncm_rng_new (NULL);
    
    ncm_rng_set_seed (abct->rng, gsl_rng_get (rng->r));

    ncm_serialize_reset (abc->ser, TRUE);

    G_UNLOCK (dup_thread);
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
  G_LOCK_DEFINE_STATIC (update_lock);
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

    G_LOCK (update_lock);
    abc->ntotal++;
    G_UNLOCK (update_lock);
    if (prob == 1.0 || (prob != 0.0 && gsl_rng_uniform (abct->rng->r) < prob))
    {
      G_LOCK (update_lock);
      abc->cur_sample_id++;
      abc->naccepted++;
      _ncm_abc_update (abc, abct->mset, dist, 1.0);
      j++;
      G_UNLOCK (update_lock);
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
 * @lre: largest relative error
 *
 * FIXME
 * 
 */
void 
ncm_abc_run_lre (NcmABC *abc, guint prerun, gdouble lre)
{
  NcmMSet *mset = ncm_mset_catalog_peek_mset (abc->mcat);
 
  gdouble lerror;
  const gdouble lre2 = lre * lre;

  g_assert_cmpfloat (lre, >, 0.0);
  /* g_assert_cmpfloat (lre, <, 1.0); */

  if (prerun == 0)
  {
    guint fparam_len = ncm_mset_fparam_len (mset);
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
 * 
 */
void
ncm_abc_mean_covar (NcmABC *abc, NcmFit *fit)
{
  NcmMSet *mset = ncm_mset_catalog_peek_mset (abc->mcat);
  
  ncm_mset_catalog_get_mean (abc->mcat, &fit->fstate->fparams);
  ncm_mset_catalog_get_covar (abc->mcat, &fit->fstate->covar);
  ncm_mset_fparams_set_vector (mset, fit->fstate->fparams);
  fit->fstate->has_covar = TRUE;
}

/**
 * ncm_abc_start_update:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 */
void 
ncm_abc_start_update (NcmABC *abc)
{
  NcmMSet *mset     = ncm_mset_catalog_peek_mset (abc->mcat);
  const gint cur_id = ncm_mset_catalog_get_cur_id (abc->mcat);
  
  if (abc->started)
    g_error ("ncm_abc_start_update: particle generation not finished, run ncm_abc_end_run() first.");
  if (abc->started_up)
    g_error ("ncm_abc_start_update: run already started, run ncm_abc_end_update() first.");

  if (abc->tkern == NULL)
    g_error ("ncm_abc_start_update: no transition kernel defined.");

  switch (abc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Starting ABC %d-th particle update (%s)...\n", abc->nupdates, ncm_abc_get_desc (abc));
      ncm_message ("%s", ncm_abc_log_info (abc));
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Model set:\n");
      ncm_mset_pretty_log (mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Starting ABC %d-th particle update (%s)...\n", abc->nupdates, ncm_abc_get_desc (abc));
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_mset_catalog_peek_rng (abc->mcat) == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);
    
    ncm_rng_set_random_seed (rng, FALSE);
    ncm_abc_set_rng (abc, rng);

    if (abc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmABC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));
    
    ncm_rng_free (rng);
  }

  {
    guint fparam_len = ncm_mset_fparam_len (mset);
    if (abc->theta != NULL)
    {
      ncm_vector_free (abc->theta);
      ncm_vector_free (abc->thetastar);
      ncm_matrix_free (abc->covar);
    }
    abc->theta     = ncm_vector_new (fparam_len);
    abc->thetastar = ncm_vector_new (fparam_len);
    abc->covar     = ncm_matrix_new (fparam_len, fparam_len);
  }
  
  ncm_mset_catalog_sync (abc->mcat, TRUE);
  if (cur_id > abc->cur_sample_id)
  {
    ncm_abc_intern_skip (abc, cur_id - abc->cur_sample_id);
    g_assert_cmpint (abc->cur_sample_id, ==, cur_id);
  }
  else if (cur_id < abc->cur_sample_id)
    g_error ("ncm_abc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
             cur_id, abc->cur_sample_id);

  g_assert_cmpuint (abc->weights->len, ==, abc->nparticles);
  g_clear_pointer (&abc->wran, gsl_ran_discrete_free);
  g_array_set_size (abc->weights_tm1, abc->nparticles);

  ncm_abc_update_tkern (abc);
  
  memcpy (abc->weights_tm1->data, abc->weights->data, sizeof (gdouble) * abc->nparticles);
  abc->wran = gsl_ran_discrete_preproc (abc->nparticles, (gdouble *)abc->weights->data);

  ncm_mset_catalog_reset_stats (abc->mcat);
  
  if (!ncm_abc_data_summary (abc))
    g_error ("ncm_abc_start_run: error calculating summary data.");

  ncm_dataset_clear (&abc->dset_mock);
  abc->dset_mock = ncm_dataset_dup (abc->dset, abc->ser);
  ncm_serialize_reset (abc->ser, TRUE);

  abc->dists_sorted = FALSE;
  abc->started_up = TRUE;
  abc->ntotal = 0;
  abc->naccepted = 0;
}

/**
 * ncm_abc_end_update:
 * @abc: a #NcmABC
 * 
 * FIXME
 * 
 */
void
ncm_abc_end_update (NcmABC *abc)
{
  NcmMSet *mset = ncm_mset_catalog_peek_mset (abc->mcat);
  guint i;
  gdouble WT = 0.0;

  g_assert (abc->started_up);
  
  if (ncm_timer_task_is_running (abc->nt))
    ncm_timer_task_end (abc->nt);

  g_clear_pointer (&abc->wran, gsl_ran_discrete_free);

  for (i = 0; i < abc->nparticles; i++)
    WT += g_array_index (abc->weights, gdouble, i);

  for (i = 0; i < abc->nparticles; i++)
    g_array_index (abc->weights, gdouble, i) = g_array_index (abc->weights, gdouble, i) / WT;

  ncm_mset_catalog_sync (abc->mcat, TRUE);
  switch (abc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    case NCM_FIT_RUN_MSGS_FULL:
      g_message ("# NcmABC:Current covariance matrix:\n");
      ncm_mset_catalog_get_mean (abc->mcat, &abc->theta);
      ncm_mset_fparams_set_vector (mset, abc->theta);
      ncm_mset_catalog_get_covar (abc->mcat, &abc->covar);
      ncm_mset_fparams_log_covar (mset, abc->covar);
      break;
    default:
      break;
  }

  abc->started_up = FALSE;
}

static void _ncm_abc_update_single (NcmABC *abc);
static void _ncm_abc_update_mt (NcmABC *abc);

/**
 * ncm_abc_update:
 * @abc: a #NcmABC
 * 
 * Runs the Monte Carlo until it reaches the @n-th realization. Note that
 * if the first_id is non-zero it will run @n - first_id realizations.
 *
 */
void 
ncm_abc_update (NcmABC *abc)
{
  if (!abc->started_up)
    g_error ("ncm_abc_update: run not started, run ncm_abc_start_update() first.");

  switch (abc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      g_message ("# NcmABC: Calculating [%06d] ABC particles updates [%s]\n", abc->n, ncm_abc_get_desc (abc));
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
    _ncm_abc_update_single (abc);
  else
    _ncm_abc_update_mt (abc);

  abc->nupdates++;
  ncm_timer_task_pause (abc->nt);
}

static void 
_ncm_abc_update_single (NcmABC *abc)
{
  NcmMSet *mset = ncm_mset_catalog_peek_mset (abc->mcat);
  NcmRNG *rng   = ncm_mset_catalog_peek_rng (abc->mcat);
  guint i = 0;
  
  for (i = 0; i < abc->n;)
  {
    gdouble dist = 0.0, prob = 0.0;
    gsize np = gsl_ran_discrete (rng->r, abc->wran);

    NcmVector *row = ncm_mset_catalog_peek_row (abc->mcat, abc->nparticles * abc->nupdates + np);
    NcmVector *theta = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);
    
    ncm_mset_trans_kern_generate (abc->tkern, theta, abc->thetastar, rng);
    ncm_mset_fparams_set_vector (mset, abc->thetastar);
    ncm_dataset_resample (abc->dset_mock, mset, rng);
    
    dist = ncm_abc_mock_distance (abc, abc->dset_mock, abc->theta, abc->thetastar, rng);
    prob = ncm_abc_distance_prob (abc, dist);
    ncm_vector_free (theta);

    abc->ntotal++;
    
    if (prob == 1.0 || (prob != 0.0 && gsl_rng_uniform (rng->r) < prob))
    {
      gdouble new_weight = ncm_mset_trans_kern_prior_pdf (abc->prior, abc->thetastar);
      gdouble denom = 0.0;
      guint j;
      for (j = 0; j < abc->nparticles; j++)
      {        
        row    = ncm_mset_catalog_peek_row (abc->mcat, abc->nparticles * abc->nupdates + j);
        theta  = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);
        denom += g_array_index (abc->weights_tm1, gdouble, j) * ncm_mset_trans_kern_pdf (abc->tkern, theta, abc->thetastar);
        ncm_vector_free (theta);
      }
      new_weight = new_weight / denom; 

      abc->naccepted++;
      abc->cur_sample_id++;
      _ncm_abc_update (abc, mset, dist, new_weight);
      i++;
    }
  }
}

static void 
_ncm_abc_thread_update_eval (glong i, glong f, gpointer data)
{
  G_LOCK_DEFINE_STATIC (update_lock);
  NcmABC *abc = NCM_ABC (data);
  NcmABCThread **abct_ptr = ncm_memory_pool_get (abc->mp);
  NcmABCThread *abct = *abct_ptr;
  guint j;

  for (j = i; j < f;)
  {
    gdouble dist = 0.0, prob = 0.0;
    gsize np = gsl_ran_discrete (abct->rng->r, abc->wran);

    NcmVector *row = ncm_mset_catalog_peek_row (abc->mcat, abc->nparticles * abc->nupdates + np);
    NcmVector *theta = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);

    ncm_mset_trans_kern_generate (abc->tkern, theta, abct->thetastar, abct->rng);

    ncm_mset_fparams_set_vector (abct->mset, abct->thetastar);
    ncm_dataset_resample (abct->dset, abct->mset, abct->rng);

    dist = ncm_abc_mock_distance (abc, abct->dset, theta, abct->thetastar, abct->rng);
    prob = ncm_abc_distance_prob (abc, dist);
/*
    G_LOCK (update_lock);
    printf ("# choice %zu\n", np);
    ncm_vector_log_vals (theta,           "# v_choice = ", "% 20.15g");
    ncm_vector_log_vals (abct->thetastar, "# v_tkern  = ", "% 20.15g");
    printf ("# dist % 20.15g prob % 20.15g\n", dist, prob);
    G_UNLOCK (update_lock);
*/
    ncm_vector_free (theta);
    
    G_LOCK (update_lock);
    abc->ntotal++;
    G_UNLOCK (update_lock);
    
    if (prob == 1.0 || (prob != 0.0 && gsl_rng_uniform (abct->rng->r) < prob))
    {
      gdouble new_weight = ncm_mset_trans_kern_prior_pdf (abc->prior, abct->thetastar);
      gdouble denom = 0.0;
      guint k;
      for (k = 0; k < abc->nparticles; k++)
      {
        row    = ncm_mset_catalog_peek_row (abc->mcat, abc->nparticles * abc->nupdates + k);
        theta  = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);
        denom += g_array_index (abc->weights_tm1, gdouble, k) * ncm_mset_trans_kern_pdf (abc->tkern, theta, abct->thetastar);
        ncm_vector_free (theta);
      }
      new_weight = new_weight / denom; 

      G_LOCK (update_lock);
      abc->cur_sample_id++;
      abc->naccepted++;
      _ncm_abc_update (abc, abct->mset, dist, new_weight);
      j++;
      G_UNLOCK (update_lock);
    }
  }

  ncm_memory_pool_return (abct_ptr);
}

static void
_ncm_abc_update_mt (NcmABC *abc)
{
  const guint nthreads = abc->n > abc->nthreads ? abc->nthreads : (abc->n - 1);

  if (nthreads == 0)
  {
    _ncm_abc_update_single (abc);
    return;
  }
  
  if (abc->mp != NULL)
    ncm_memory_pool_free (abc->mp, TRUE);
  abc->mp = ncm_memory_pool_new (&_ncm_abc_dup_thread, abc, 
                                 (GDestroyNotify) &_ncm_abc_free_thread);
  
  g_assert_cmpuint (abc->nthreads, >, 1);

  ncm_func_eval_threaded_loop_full (&_ncm_abc_thread_update_eval, 0, abc->n, abc);
}


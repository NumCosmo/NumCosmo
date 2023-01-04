/***************************************************************************
 *            ncm_abc.c
 *
 *  Tue September 30 15:46:48 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_abc.c
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

typedef struct _NcmABCPrivate
{
  NcmMSetCatalog *mcat;
  NcmDataset *dset;
  NcmDataset *dset_mock;
  NcmMemoryPool *mp;
  NcmMSetTransKern *prior;
  NcmMSetTransKern *tkern;
  NcmTimer *nt;
  NcmSerialize *ser;
  NcmFitRunMsgs mtype;
  NcmVector *theta;
  NcmVector *thetastar;
  NcmMatrix *covar;
  GArray *weights;
  GArray *weights_tm1;
  GArray *pchoice;
  GArray *dists;
  gdouble epsilon;
  gdouble depsilon;
  gboolean dists_sorted;
  gsl_ran_discrete_t *wran;
  gboolean started;
  gboolean started_up;
  gint cur_sample_id;
  guint ntotal;
  guint naccepted;
  guint nthreads;
  guint nupdates;
  guint n;
  guint nparticles;
} NcmABCPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmABC, ncm_abc, G_TYPE_OBJECT);

static void
ncm_abc_init (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  self->mcat          = NULL;
  self->dset          = NULL;
  self->dset_mock     = NULL;
  self->mp            = NULL;
  self->tkern         = NULL;
  self->prior         = NULL;
  self->nt            = ncm_timer_new ();
  self->ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  self->mtype         = NCM_FIT_RUN_MSGS_NONE;
  self->theta         = NULL;
  self->thetastar     = NULL;
  self->covar         = NULL;
  self->weights       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->weights_tm1   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->dists         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->dists_sorted  = FALSE;
  self->epsilon       = 0.0;
  self->depsilon      = 0.0;
  self->wran          = NULL;
  self->started       = FALSE;
  self->started_up    = FALSE;
  self->cur_sample_id = -1; /* Represents that no samples were calculated yet. */
  self->nthreads      = 0;
  self->nparticles    = 0;
  self->n             = 0;
  self->ntotal        = 0;
  self->naccepted     = 0;
  self->nupdates      = 0;
}

static void
_ncm_abc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmABC *abc                = NCM_ABC (object);
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  g_return_if_fail (NCM_IS_ABC (object));

  switch (prop_id)
  {
    case PROP_MSET:
    {
      NcmMSet *mset = g_value_get_object (value);

      g_assert (self->mcat == NULL);

      if (mset == NULL)
        g_error ("ncm_abc_new: mset (NcmMSet) cannot be NULL.");

      self->mcat = ncm_mset_catalog_new (mset, 1, 1, TRUE, "NcmABC:Distance", "d", NULL);
      break;
    }
    case PROP_PRIOR:
    {
      NcmMSetTransKern *prior = g_value_dup_object (value);

      g_assert (self->prior == NULL);

      if (prior == NULL)
        g_error ("ncm_abc_new: prior (NcmMSetTransKern) cannot be NULL.");

      self->prior = prior;
      break;
    }
    case PROP_TKERN:
      ncm_abc_set_trans_kern (abc, g_value_get_object (value));
      break;
    case PROP_DATASET:
    {
      NcmDataset *dset = g_value_dup_object (value);

      g_assert (self->dset == NULL);

      if (dset == NULL)
        g_error ("ncm_abc_new: dset (NcmDataset) cannot be NULL.");

      self->dset = dset;
      break;
    }
    case PROP_EPSILON:
      self->epsilon  = g_value_get_double (value);
      self->depsilon = self->epsilon;
      break;
    case PROP_NPARTICLES:
      g_assert (!self->started);
      g_assert_not_reached ();
      self->nparticles = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_abc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmABC *abc                = NCM_ABC (object);
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  g_return_if_fail (NCM_IS_ABC (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, ncm_mset_catalog_peek_mset (self->mcat));
      break;
    case PROP_PRIOR:
      g_value_set_object (value, self->prior);
      break;
    case PROP_TKERN:
      g_value_set_object (value, self->tkern);
      break;
    case PROP_DATASET:
      g_value_set_object (value, self->dset);
      break;
    case PROP_EPSILON:
      g_value_set_double (value, self->epsilon);
      break;
    case PROP_NPARTICLES:
      g_value_set_uint (value, self->nparticles);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_abc_dispose (GObject *object)
{
  NcmABC *abc                = NCM_ABC (object);
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  ncm_mset_catalog_clear (&self->mcat);
  ncm_mset_trans_kern_clear (&self->prior);
  ncm_mset_trans_kern_clear (&self->tkern);
  ncm_dataset_clear (&self->dset);
  ncm_dataset_clear (&self->dset_mock);
  ncm_timer_clear (&self->nt);
  ncm_serialize_clear (&self->ser);

  ncm_vector_clear (&self->theta);
  ncm_vector_clear (&self->thetastar);
  ncm_matrix_clear (&self->covar);

  g_clear_pointer (&self->weights, g_array_unref);
  g_clear_pointer (&self->weights_tm1, g_array_unref);
  g_clear_pointer (&self->dists, g_array_unref);
  g_clear_pointer (&self->wran, gsl_ran_discrete_free);

  if (self->mp != NULL)
  {
    ncm_memory_pool_free (self->mp, TRUE);
    self->mp = NULL;
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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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
 * Updates the transition kernel present in @self->tkern.
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  self->mtype = mtype;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  const gchar *cur_filename  = ncm_mset_catalog_peek_filename (self->mcat);

  if ((self->cur_sample_id >= 0) && (cur_filename != NULL))
    g_error ("ncm_abc_set_data_file: Cannot change data file during a run, call ncm_abc_end_run() first.");

  if ((cur_filename != NULL) && (strcmp (cur_filename, filename) == 0))
    return;

  ncm_mset_catalog_set_file (self->mcat, filename);

  if (self->started)
    g_assert_cmpint (self->cur_sample_id, ==, ncm_mset_catalog_get_cur_id (self->mcat));
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  self->nthreads = nthreads;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (self->cur_sample_id >= 0)
    g_error ("ncm_abc_set_rng: Cannot change the RNG object in a non-empty catalog.");

  ncm_mset_catalog_set_rng (self->mcat, rng);
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (self->tkern != NULL)
    g_warning ("Transition kernel already set, replacing.");

  ncm_mset_trans_kern_clear (&self->tkern);
  self->tkern = ncm_mset_trans_kern_ref (tkern);
}

static gint
_compare (gconstpointer a, gconstpointer b)
{
  if (*(gdouble *) a > *(gdouble *) b)
    return 1;
  else if (*(gdouble *) a < *(gdouble *) b)
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (self->started)
    g_error ("ncm_abc_get_dist_quantile: Cannot get quantiles during a run, call ncm_abc_end_run() first.");

  if (self->dists->len < 1)
    g_error ("ncm_abc_get_dist_quantile: Cannot get quantiles, no particles calculated.");

  if (!self->dists_sorted)
  {
    g_array_sort (self->dists, &_compare);
    self->dists_sorted = TRUE;
  }

  return gsl_stats_quantile_from_sorted_data ((gdouble *) self->dists->data, 1, self->dists->len, p);
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->naccepted * 1.0 / (self->ntotal * 1.0);
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (epsilon >= self->epsilon)
    g_warning ("ncm_abc_update_epsilon: increasing epsilon.");

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
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
               self->epsilon);
    g_message ("# NcmABC: epsilon_t+1      = %g.\n",
               epsilon);
    g_message ("# NcmABC: depsilon/epsilon = %04.2f%%.\n", 100.0 * (self->epsilon - epsilon) / self->epsilon);
  }

  self->depsilon = self->epsilon - epsilon;
  self->epsilon  = epsilon;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->epsilon;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->depsilon;
}

/**
 * ncm_abc_get_mtype:
 * @abc: a #NcmABC
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitRunMsgs
ncm_abc_get_mtype (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->mtype;
}

/**
 * ncm_abc_peek_catalog:
 * @abc: a #NcmABC
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmMSetCatalog *
ncm_abc_peek_catalog (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->mcat;
}

/**
 * ncm_abc_peek_dataset:
 * @abc: a #NcmABC
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmDataset *
ncm_abc_peek_dataset (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->dset;
}

/**
 * ncm_abc_peek_covar:
 * @abc: a #NcmABC
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmMatrix *
ncm_abc_peek_covar (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->covar;
}

/**
 * ncm_abc_peek_trans_kern:
 * @abc: a #NcmABC
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmMSetTransKern *
ncm_abc_peek_trans_kern (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  return self->tkern;
}

void
_ncm_abc_update (NcmABC *abc, NcmMSet *mset, gdouble dist, gdouble weight)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  const guint part   = 5;
  const guint step   = (self->n / part) == 0 ? 1 : (self->n / part);
  const guint pindex = self->cur_sample_id % self->nparticles;

  ncm_mset_catalog_add_from_mset (self->mcat, mset, dist, weight, NULL);
  ncm_timer_task_increment (self->nt);
  g_array_index (self->weights, gdouble, pindex) = weight;
  g_array_index (self->dists,   gdouble, pindex) = dist;
  self->dists_sorted                             = FALSE;

  switch (self->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi          = (self->cur_sample_id + 1) % step;
      gboolean log_timeout = FALSE;

      if ((self->nt->pos_time - self->nt->last_log_time) > 60.0)
        log_timeout = TRUE;

      if (log_timeout || (stepi == 0) || (self->nt->task_pos == self->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (self->mcat);
        g_message ("# NcmABC:acceptance ratio %7.4f%%.\n", ncm_abc_get_accept_ratio (abc) * 100.0);
        /* ncm_timer_task_accumulate (self->nt, acc); */
        ncm_timer_task_log_elapsed (self->nt);
        ncm_timer_task_log_mean_time (self->nt);
        ncm_timer_task_log_time_left (self->nt);
        ncm_timer_task_log_end_datetime (self->nt);
        ncm_cfg_logfile_flush_now ();
      }

      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_mset_catalog_log_current_stats (self->mcat);
      g_message ("# NcmABC:acceptance ratio %7.4f%%.\n", ncm_abc_get_accept_ratio (abc) * 100.0);
      /* ncm_timer_task_increment (self->nt); */
      ncm_timer_task_log_elapsed (self->nt);
      ncm_timer_task_log_mean_time (self->nt);
      ncm_timer_task_log_time_left (self->nt);
      ncm_timer_task_log_end_datetime (self->nt);
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  NcmMSet *mset     = ncm_mset_catalog_peek_mset (self->mcat);
  const gint cur_id = ncm_mset_catalog_get_cur_id (self->mcat);

  if (self->started)
    g_error ("ncm_abc_start_run: run already started, run ncm_abc_end_run() first.");

  switch (self->mtype)
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

  if (ncm_mset_catalog_peek_rng (self->mcat) == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);

    ncm_rng_set_random_seed (rng, FALSE);
    ncm_abc_set_rng (abc, rng);

    if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmABC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));

    ncm_rng_free (rng);
  }

  self->dists_sorted = FALSE;
  self->started      = TRUE;

  {
    guint fparam_len = ncm_mset_fparam_len (mset);

    if (self->theta != NULL)
    {
      ncm_vector_free (self->theta);
      ncm_vector_free (self->thetastar);
      ncm_matrix_free (self->covar);
    }

    self->theta     = ncm_vector_new (fparam_len);
    self->thetastar = ncm_vector_new (fparam_len);
    self->covar     = ncm_matrix_new (fparam_len, fparam_len);
  }

  ncm_mset_catalog_set_sync_mode (self->mcat, NCM_MSET_CATALOG_SYNC_TIMED);
  ncm_mset_catalog_set_sync_interval (self->mcat, NCM_ABC_MIN_SYNC_INTERVAL);

  ncm_mset_catalog_sync (self->mcat, TRUE);

  if (cur_id > self->cur_sample_id)
  {
    ncm_abc_intern_skip (abc, cur_id - self->cur_sample_id);
    g_assert_cmpint (self->cur_sample_id, ==, cur_id);
  }
  else if (cur_id < self->cur_sample_id)
  {
    g_error ("ncm_abc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].",
             cur_id, self->cur_sample_id);
  }

  {
    NcmVector *cur_row = NULL;

    cur_row = ncm_mset_catalog_peek_current_row (self->mcat);

    if (cur_row != NULL)
      ncm_mset_fparams_set_vector_offset (mset, cur_row, 2);
  }

  if (!ncm_abc_data_summary (abc))
    g_error ("ncm_abc_start_run: error calculating summary data.");

  ncm_dataset_clear (&self->dset_mock);
  self->dset_mock = ncm_dataset_dup (self->dset, self->ser);
  ncm_serialize_reset (self->ser, TRUE);

  self->ntotal    = 0;
  self->naccepted = 0;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  NcmMSet *mset = ncm_mset_catalog_peek_mset (self->mcat);
  gdouble WT    = 0.0;
  guint i;

  g_assert (self->started);

  if (ncm_timer_task_is_running (self->nt))
    ncm_timer_task_end (self->nt);

  ncm_mset_catalog_sync (self->mcat, TRUE);

  for (i = 0; i < self->nparticles; i++)
    WT += g_array_index (self->weights, gdouble, i);

  for (i = 0; i < self->nparticles; i++)
    g_array_index (self->weights, gdouble, i) = g_array_index (self->weights, gdouble, i) / WT;

  switch (self->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    case NCM_FIT_RUN_MSGS_FULL:
      g_message ("# NcmABC:Current covariance matrix:\n");
      ncm_mset_catalog_get_mean (self->mcat, &self->theta);
      ncm_mset_fparams_set_vector (mset, self->theta);
      ncm_mset_catalog_get_covar (self->mcat, &self->covar);
      ncm_mset_fparams_log_covar (mset, self->covar);
      break;
    default:
      break;
  }

  self->started = FALSE;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  self->n             = 0;
  self->cur_sample_id = -1;
  self->started       = FALSE;
  self->ntotal        = 0;
  self->naccepted     = 0;
  self->nupdates      = 0;
  ncm_mset_catalog_reset (self->mcat);
}

static void
ncm_abc_intern_skip (NcmABC *abc, guint n)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (n == 0)
    return;

  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Skipping %u tries, will start at %u-th try.\n", n, self->cur_sample_id + n + 1 + 1);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  self->cur_sample_id += n;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (!self->started)
    g_error ("ncm_abc_run: run not started, run ncm_abc_start_run() first.");

  if (nparticles <= (self->cur_sample_id + 1))
  {
    if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Nothing to do, current ABC particle number is already %d of %u\n", self->cur_sample_id + 1, nparticles);
    }

    return;
  }

  self->n = nparticles - (self->cur_sample_id + 1);

  if ((self->n > 0) && (self->nupdates > 0))
    g_error ("ncm_abc_run: cannot generate new particles when time t = %u != 0.", self->nupdates);

  self->nparticles = nparticles;
  g_array_set_size (self->weights, nparticles);
  g_array_set_size (self->dists, nparticles);

  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Calculating [%06d] ABC particles [%s]\n", self->n, ncm_abc_get_desc (abc));
    }
    break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_timer_task_is_running (self->nt))
  {
    ncm_timer_task_add_tasks (self->nt, self->n);
    ncm_timer_task_continue (self->nt);
  }
  else
  {
    ncm_timer_task_start (self->nt, self->n);
    ncm_timer_set_name (self->nt, "NcmABC");
  }

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (self->nt);

  if (self->nthreads <= 1)
    _ncm_abc_run_single (abc);
  else
    _ncm_abc_run_mt (abc);

  ncm_timer_task_pause (self->nt);

  g_assert_cmpuint (self->nparticles, ==, self->cur_sample_id + 1);
}

static void
_ncm_abc_run_single (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  NcmMSet *mset = ncm_mset_catalog_peek_mset (self->mcat);
  NcmRNG *rng   = ncm_mset_catalog_peek_rng (self->mcat);

  guint i = 0;

  for (i = 0; i < self->n;)
  {
    gdouble dist = 0.0, prob = 0.0;

    ncm_mset_trans_kern_prior_sample (self->prior, self->thetastar, rng);

    ncm_mset_fparams_set_vector (mset, self->thetastar);
    ncm_dataset_resample (self->dset_mock, mset, rng);

    dist = ncm_abc_mock_distance (abc, self->dset_mock, self->theta, self->thetastar, rng);
    prob = ncm_abc_distance_prob (abc, dist);

    self->ntotal++;

    if ((prob == 1.0) || ((prob != 0.0) && (gsl_rng_uniform (rng->r) < prob)))
    {
      self->cur_sample_id++;
      self->naccepted++;
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

  NcmABC *abc                = NCM_ABC (userdata);
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  NcmABCThread *abct         = g_new (NcmABCThread, 1);
  NcmMSet *mset              = ncm_mset_catalog_peek_mset (self->mcat);
  NcmRNG *rng                = ncm_mset_catalog_peek_rng (self->mcat);

  G_LOCK (dup_thread);
  {
    abct->mset      = ncm_mset_dup (mset, self->ser);
    abct->dset      = ncm_dataset_dup (self->dset, self->ser);
    abct->thetastar = ncm_vector_dup (self->thetastar);
    abct->rng       = ncm_rng_new (NULL);

    ncm_rng_set_seed (abct->rng, gsl_rng_get (rng->r));

    ncm_serialize_reset (self->ser, TRUE);

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

  NcmABC *abc                = NCM_ABC (data);
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  NcmABCThread **abct_ptr    = ncm_memory_pool_get (self->mp);
  NcmABCThread *abct         = *abct_ptr;
  guint j;

  for (j = i; j < f;)
  {
    gdouble dist, prob;

    ncm_mset_trans_kern_prior_sample (self->prior, abct->thetastar, abct->rng);

    ncm_mset_fparams_set_vector (abct->mset, abct->thetastar);

    ncm_dataset_resample (abct->dset, abct->mset, abct->rng);
    dist = ncm_abc_mock_distance (abc, abct->dset, abct->thetastar, abct->thetastar, abct->rng);
    prob = ncm_abc_distance_prob (abc, dist);

    G_LOCK (update_lock);
    self->ntotal++;
    G_UNLOCK (update_lock);

    if ((prob == 1.0) || ((prob != 0.0) && (gsl_rng_uniform (abct->rng->r) < prob)))
    {
      G_LOCK (update_lock);
      self->cur_sample_id++;
      self->naccepted++;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  const guint nthreads = self->n > self->nthreads ? self->nthreads : (self->n - 1);

  if (nthreads == 0)
  {
    _ncm_abc_run_single (abc);

    return;
  }

  if (self->mp != NULL)
    ncm_memory_pool_free (self->mp, TRUE);

  self->mp = ncm_memory_pool_new (&_ncm_abc_dup_thread, abc,
                                  (GDestroyNotify) & _ncm_abc_free_thread);

  g_assert_cmpuint (self->nthreads, >, 1);

  ncm_func_eval_threaded_loop_full (&_ncm_abc_thread_eval, 0, self->n, abc);
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  NcmMSet *mset = ncm_mset_catalog_peek_mset (self->mcat);

  gdouble lerror;
  const gdouble lre2 = lre * lre;

  g_assert_cmpfloat (lre, >, 0.0);
  /* g_assert_cmpfloat (lre, <, 1.0); */

  if (prerun == 0)
  {
    guint fparam_len = ncm_mset_fparam_len (mset);

    prerun = fparam_len * 100;
  }

  if (ncm_mset_catalog_len (self->mcat) < prerun)
  {
    guint prerun_left = prerun - ncm_mset_catalog_len (self->mcat);

    if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmABC: Running first %u pre-runs...\n", prerun_left);

    ncm_abc_run (abc, prerun);
  }

  lerror = ncm_mset_catalog_largest_error (self->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n             = ncm_mset_catalog_len (self->mcat);
    gdouble m             = n * lerror2 / lre2;
    guint runs            = ((m - n) > 1000.0) ? ceil ((m - n) * 0.25) : ceil (m - n);

    if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmABC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmABC: Running more %u runs...\n", runs);
    }

    ncm_abc_run (abc, self->cur_sample_id + runs + 1);
    lerror = ncm_mset_catalog_largest_error (self->mcat);
  }

  if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  NcmMSet *mset = ncm_mset_catalog_peek_mset (self->mcat);

  ncm_mset_catalog_get_mean (self->mcat, &fit->fstate->fparams);
  ncm_mset_catalog_get_covar (self->mcat, &fit->fstate->covar);
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  NcmMSet *mset              = ncm_mset_catalog_peek_mset (self->mcat);
  const gint cur_id          = ncm_mset_catalog_get_cur_id (self->mcat);

  if (self->started)
    g_error ("ncm_abc_start_update: particle generation not finished, run ncm_abc_end_run() first.");

  if (self->started_up)
    g_error ("ncm_abc_start_update: run already started, run ncm_abc_end_update() first.");

  if (self->tkern == NULL)
    g_error ("ncm_abc_start_update: no transition kernel defined.");

  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Starting ABC %d-th particle update (%s)...\n", self->nupdates, ncm_abc_get_desc (abc));
      ncm_message ("%s", ncm_abc_log_info (abc));
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Model set:\n");
      ncm_mset_pretty_log (mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmABC: Starting ABC %d-th particle update (%s)...\n", self->nupdates, ncm_abc_get_desc (abc));
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_mset_catalog_peek_rng (self->mcat) == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);

    ncm_rng_set_random_seed (rng, FALSE);
    ncm_abc_set_rng (abc, rng);

    if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmABC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));

    ncm_rng_free (rng);
  }

  {
    guint fparam_len = ncm_mset_fparam_len (mset);

    if (self->theta != NULL)
    {
      ncm_vector_free (self->theta);
      ncm_vector_free (self->thetastar);
      ncm_matrix_free (self->covar);
    }

    self->theta     = ncm_vector_new (fparam_len);
    self->thetastar = ncm_vector_new (fparam_len);
    self->covar     = ncm_matrix_new (fparam_len, fparam_len);
  }

  ncm_mset_catalog_sync (self->mcat, TRUE);

  if (cur_id > self->cur_sample_id)
  {
    ncm_abc_intern_skip (abc, cur_id - self->cur_sample_id);
    g_assert_cmpint (self->cur_sample_id, ==, cur_id);
  }
  else if (cur_id < self->cur_sample_id)
  {
    g_error ("ncm_abc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].",
             cur_id, self->cur_sample_id);
  }

  g_assert_cmpuint (self->weights->len, ==, self->nparticles);
  g_clear_pointer (&self->wran, gsl_ran_discrete_free);
  g_array_set_size (self->weights_tm1, self->nparticles);

  ncm_abc_update_tkern (abc);

  memcpy (self->weights_tm1->data, self->weights->data, sizeof (gdouble) * self->nparticles);
  self->wran = gsl_ran_discrete_preproc (self->nparticles, (gdouble *) self->weights->data);

  ncm_mset_catalog_reset_stats (self->mcat);

  if (!ncm_abc_data_summary (abc))
    g_error ("ncm_abc_start_run: error calculating summary data.");

  ncm_dataset_clear (&self->dset_mock);
  self->dset_mock = ncm_dataset_dup (self->dset, self->ser);
  ncm_serialize_reset (self->ser, TRUE);

  self->dists_sorted = FALSE;
  self->started_up   = TRUE;
  self->ntotal       = 0;
  self->naccepted    = 0;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  NcmMSet *mset              = ncm_mset_catalog_peek_mset (self->mcat);
  guint i;
  gdouble WT = 0.0;

  g_assert (self->started_up);

  if (ncm_timer_task_is_running (self->nt))
    ncm_timer_task_end (self->nt);

  g_clear_pointer (&self->wran, gsl_ran_discrete_free);

  for (i = 0; i < self->nparticles; i++)
    WT += g_array_index (self->weights, gdouble, i);

  for (i = 0; i < self->nparticles; i++)
    g_array_index (self->weights, gdouble, i) = g_array_index (self->weights, gdouble, i) / WT;

  ncm_mset_catalog_sync (self->mcat, TRUE);

  switch (self->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    case NCM_FIT_RUN_MSGS_FULL:
      g_message ("# NcmABC:Current covariance matrix:\n");
      ncm_mset_catalog_get_mean (self->mcat, &self->theta);
      ncm_mset_fparams_set_vector (mset, self->theta);
      ncm_mset_catalog_get_covar (self->mcat, &self->covar);
      ncm_mset_fparams_log_covar (mset, self->covar);
      break;
    default:
      break;
  }

  self->started_up = FALSE;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  if (!self->started_up)
    g_error ("ncm_abc_update: run not started, run ncm_abc_start_update() first.");

  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      g_message ("# NcmABC: Calculating [%06d] ABC particles updates [%s]\n", self->n, ncm_abc_get_desc (abc));
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_timer_task_is_running (self->nt))
  {
    ncm_timer_task_add_tasks (self->nt, self->n);
    ncm_timer_task_continue (self->nt);
  }
  else
  {
    ncm_timer_task_start (self->nt, self->n);
    ncm_timer_set_name (self->nt, "NcmABC");
  }

  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (self->nt);

  if (self->nthreads <= 1)
    _ncm_abc_update_single (abc);
  else
    _ncm_abc_update_mt (abc);

  self->nupdates++;
  ncm_timer_task_pause (self->nt);
}

static void
_ncm_abc_update_single (NcmABC *abc)
{
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);

  NcmMSet *mset = ncm_mset_catalog_peek_mset (self->mcat);
  NcmRNG *rng   = ncm_mset_catalog_peek_rng (self->mcat);
  guint i       = 0;

  for (i = 0; i < self->n;)
  {
    gdouble dist = 0.0, prob = 0.0;
    gsize np = gsl_ran_discrete (rng->r, self->wran);

    NcmVector *row   = ncm_mset_catalog_peek_row (self->mcat, self->nparticles * self->nupdates + np);
    NcmVector *theta = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);

    ncm_mset_trans_kern_generate (self->tkern, theta, self->thetastar, rng);
    ncm_mset_fparams_set_vector (mset, self->thetastar);
    ncm_dataset_resample (self->dset_mock, mset, rng);

    dist = ncm_abc_mock_distance (abc, self->dset_mock, self->theta, self->thetastar, rng);
    prob = ncm_abc_distance_prob (abc, dist);
    ncm_vector_free (theta);

    self->ntotal++;

    if ((prob == 1.0) || ((prob != 0.0) && (gsl_rng_uniform (rng->r) < prob)))
    {
      gdouble new_weight = ncm_mset_trans_kern_prior_pdf (self->prior, self->thetastar);
      gdouble denom      = 0.0;
      guint j;

      for (j = 0; j < self->nparticles; j++)
      {
        row    = ncm_mset_catalog_peek_row (self->mcat, self->nparticles * self->nupdates + j);
        theta  = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);
        denom += g_array_index (self->weights_tm1, gdouble, j) * ncm_mset_trans_kern_pdf (self->tkern, theta, self->thetastar);
        ncm_vector_free (theta);
      }

      new_weight = new_weight / denom;

      self->naccepted++;
      self->cur_sample_id++;
      _ncm_abc_update (abc, mset, dist, new_weight);
      i++;
    }
  }
}

static void
_ncm_abc_thread_update_eval (glong i, glong f, gpointer data)
{
  G_LOCK_DEFINE_STATIC (update_lock);

  NcmABC *abc                = NCM_ABC (data);
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  NcmABCThread **abct_ptr    = ncm_memory_pool_get (self->mp);
  NcmABCThread *abct         = *abct_ptr;
  guint j;

  for (j = i; j < f;)
  {
    gdouble dist = 0.0, prob = 0.0;
    gsize np = gsl_ran_discrete (abct->rng->r, self->wran);

    NcmVector *row   = ncm_mset_catalog_peek_row (self->mcat, self->nparticles * self->nupdates + np);
    NcmVector *theta = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);

    ncm_mset_trans_kern_generate (self->tkern, theta, abct->thetastar, abct->rng);

    ncm_mset_fparams_set_vector (abct->mset, abct->thetastar);
    ncm_dataset_resample (abct->dset, abct->mset, abct->rng);

    dist = ncm_abc_mock_distance (abc, abct->dset, theta, abct->thetastar, abct->rng);
    prob = ncm_abc_distance_prob (abc, dist);

/*
 *   G_LOCK (update_lock);
 *   printf ("# choice %zu\n", np);
 *   ncm_vector_log_vals (theta,           "# v_choice = ", "% 20.15g");
 *   ncm_vector_log_vals (abct->thetastar, "# v_tkern  = ", "% 20.15g");
 *   printf ("# dist % 20.15g prob % 20.15g\n", dist, prob);
 *   G_UNLOCK (update_lock);
 */
    ncm_vector_free (theta);

    G_LOCK (update_lock);
    self->ntotal++;
    G_UNLOCK (update_lock);

    if ((prob == 1.0) || ((prob != 0.0) && (gsl_rng_uniform (abct->rng->r) < prob)))
    {
      gdouble new_weight = ncm_mset_trans_kern_prior_pdf (self->prior, abct->thetastar);
      gdouble denom      = 0.0;
      guint k;

      for (k = 0; k < self->nparticles; k++)
      {
        row    = ncm_mset_catalog_peek_row (self->mcat, self->nparticles * self->nupdates + k);
        theta  = ncm_vector_get_subvector (row, 2, ncm_vector_len (row) - 2);
        denom += g_array_index (self->weights_tm1, gdouble, k) * ncm_mset_trans_kern_pdf (self->tkern, theta, abct->thetastar);
        ncm_vector_free (theta);
      }

      new_weight = new_weight / denom;

      G_LOCK (update_lock);
      self->cur_sample_id++;
      self->naccepted++;
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
  NcmABCPrivate * const self = ncm_abc_get_instance_private (abc);
  const guint nthreads       = self->n > self->nthreads ? self->nthreads : (self->n - 1);

  if (nthreads == 0)
  {
    _ncm_abc_update_single (abc);

    return;
  }

  if (self->mp != NULL)
    ncm_memory_pool_free (self->mp, TRUE);

  self->mp = ncm_memory_pool_new (&_ncm_abc_dup_thread, abc,
                                  (GDestroyNotify) & _ncm_abc_free_thread);

  g_assert_cmpuint (self->nthreads, >, 1);

  ncm_func_eval_threaded_loop_full (&_ncm_abc_thread_update_eval, 0, self->n, abc);
}


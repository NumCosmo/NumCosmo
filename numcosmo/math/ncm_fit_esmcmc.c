/***************************************************************************
 *            ncm_fit_esmcmc.c
 *
 *  Tue January 20 16:59:36 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti & Mariana Penna-Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit_esmcmc
 * @title: NcmFitESMCMC
 * @short_description: Ensemble sampler Markov Chain Monte Carlo analysis.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc.h"
#include "math/ncm_cfg.h"
#include "math/ncm_func_eval.h"
#include "ncm_enum_types.h"

#include <gsl/gsl_statistics_double.h>

enum
{
  PROP_0,
  PROP_FIT,
  PROP_NWALKERS,
  PROP_SAMPLER,
  PROP_MOVE_TYPE,
  PROP_A,
  PROP_MTYPE,
  PROP_NTHREADS,
  PROP_DATA_FILE,
};

G_DEFINE_TYPE (NcmFitESMCMC, ncm_fit_esmcmc, G_TYPE_OBJECT);

static void
ncm_fit_esmcmc_init (NcmFitESMCMC *esmcmc)
{
  esmcmc->fit             = NULL;
  esmcmc->walker_fits     = g_ptr_array_new ();
  g_ptr_array_set_free_func (esmcmc->walker_fits, (GDestroyNotify) &ncm_fit_free);

  esmcmc->sampler         = NULL;
  esmcmc->mcat            = NULL;
  esmcmc->mtype           = NCM_FIT_RUN_MSGS_NONE;
  esmcmc->nt              = ncm_timer_new ();
  esmcmc->ser             = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  esmcmc->mt              = NCM_FIT_ESMCMC_MOVE_TYPE_LEN;
  esmcmc->a               = 0.0;
  esmcmc->fparam_len      = 0;

  esmcmc->theta_k         = g_ptr_array_new ();
  esmcmc->theta_j         = g_ptr_array_new ();
  esmcmc->thetastar       = g_ptr_array_new ();
  g_ptr_array_set_free_func (esmcmc->theta_k, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (esmcmc->theta_j, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (esmcmc->thetastar, (GDestroyNotify) &ncm_vector_free);
  
  esmcmc->nthreads        = 0;
  esmcmc->n               = 0;
  esmcmc->nwalkers        = 0;
  esmcmc->cur_sample_id   = -1; /* Represents that no samples were calculated yet. */
  esmcmc->ntotal          = 0;
  esmcmc->naccepted       = 0;
  esmcmc->write_index     = 0;
  esmcmc->started         = FALSE;

  g_mutex_init (&esmcmc->dup_fit);
  g_mutex_init (&esmcmc->resample_lock);
  g_mutex_init (&esmcmc->update_lock);
  g_cond_init (&esmcmc->write_cond);
}

static void
_ncm_fit_esmcmc_constructed (GObject *object)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
  guint k;
  g_assert_cmpint (esmcmc->nwalkers, >, 0);
  esmcmc->fparam_len = ncm_mset_fparam_len (esmcmc->fit->mset);

  esmcmc->mcat = ncm_mset_catalog_new (esmcmc->fit->mset, 1, esmcmc->nwalkers, FALSE, NCM_MSET_CATALOG_M2LNL_COLNAME);
  ncm_mset_catalog_set_run_type (esmcmc->mcat, "Ensemble Sampler MCMC -- Stretch Move");
  
  for (k = 0; k < esmcmc->nwalkers; k++)
  {
    NcmFit *fit = ncm_fit_dup (esmcmc->fit, esmcmc->ser);

    NcmVector *theta_k = ncm_vector_new (esmcmc->fparam_len);
    NcmVector *theta_j = ncm_vector_new (esmcmc->fparam_len);
    NcmVector *thetastar = ncm_vector_new (esmcmc->fparam_len);
     
    ncm_serialize_clear_instances (esmcmc->ser);

    g_ptr_array_add (esmcmc->theta_k, theta_k);
    g_ptr_array_add (esmcmc->theta_j, theta_j);
    g_ptr_array_add (esmcmc->thetastar, thetastar);
    
    g_ptr_array_add (esmcmc->walker_fits, fit);
  }
}

static void _ncm_fit_esmcmc_set_fit_obj (NcmFitESMCMC *esmcmc, NcmFit *fit);

static void
ncm_fit_esmcmc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      _ncm_fit_esmcmc_set_fit_obj (esmcmc, g_value_get_object (value));
      break;
    case PROP_NWALKERS:
      esmcmc->nwalkers = g_value_get_int (value);
      break;      
    case PROP_SAMPLER:
      ncm_fit_esmcmc_set_sampler (esmcmc, g_value_get_object (value));
      break;      
    case PROP_MOVE_TYPE:
      ncm_fit_esmcmc_set_move_type (esmcmc, g_value_get_enum (value));
      break;
    case PROP_A:
      ncm_fit_esmcmc_set_move_scale (esmcmc, g_value_get_double (value));
      break;
    case PROP_MTYPE:
      ncm_fit_esmcmc_set_mtype (esmcmc, g_value_get_enum (value));
      break;
    case PROP_NTHREADS:
      ncm_fit_esmcmc_set_nthreads (esmcmc, g_value_get_uint (value));
      break;
    case PROP_DATA_FILE:
      ncm_fit_esmcmc_set_data_file (esmcmc, g_value_get_string (value));
      break;    
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_esmcmc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, esmcmc->fit);
      break;
    case PROP_NWALKERS:
      g_value_set_int (value, esmcmc->nwalkers);
      break;
    case PROP_SAMPLER:
      g_value_set_object (value, esmcmc->sampler);
      break;      
    case PROP_MOVE_TYPE:
      g_value_set_enum (value, esmcmc->mt);
      break;
    case PROP_A:
      g_value_set_double (value, esmcmc->a);
      break;
    case PROP_MTYPE:
      g_value_set_enum (value, esmcmc->mtype);
      break;
    case PROP_NTHREADS:
      g_value_set_uint (value, esmcmc->nthreads);
      break;
    case PROP_DATA_FILE:
      g_value_set_string (value, ncm_mset_catalog_peek_filename (esmcmc->mcat));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_esmcmc_dispose (GObject *object)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);

  ncm_fit_clear (&esmcmc->fit);
  ncm_mset_trans_kern_clear (&esmcmc->sampler);
  ncm_timer_clear (&esmcmc->nt);
  ncm_serialize_clear (&esmcmc->ser);
  ncm_mset_catalog_clear (&esmcmc->mcat);

  g_clear_pointer (&esmcmc->theta_k, g_ptr_array_unref);
  g_clear_pointer (&esmcmc->theta_j, g_ptr_array_unref);
  g_clear_pointer (&esmcmc->thetastar, g_ptr_array_unref);
  
  g_clear_pointer (&esmcmc->walker_fits, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_parent_class)->dispose (object);
}

static void
ncm_fit_esmcmc_finalize (GObject *object)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);

  g_mutex_clear (&esmcmc->dup_fit);
  g_mutex_clear (&esmcmc->resample_lock);
  g_mutex_clear (&esmcmc->update_lock);
  g_cond_clear (&esmcmc->write_cond);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_parent_class)->finalize (object);
}

static void
ncm_fit_esmcmc_class_init (NcmFitESMCMCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_esmcmc_constructed;
  object_class->set_property = &ncm_fit_esmcmc_set_property;
  object_class->get_property = &ncm_fit_esmcmc_get_property;
  object_class->dispose      = &ncm_fit_esmcmc_dispose;
  object_class->finalize     = &ncm_fit_esmcmc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_object ("fit",
                                                        NULL,
                                                        "Fit object",
                                                        NCM_TYPE_FIT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NWALKERS,
                                   g_param_spec_int ("nwalkers",
                                                      NULL,
                                                      "Number of walkers",
                                                      1, G_MAXINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SAMPLER,
                                   g_param_spec_object ("sampler",
                                                        NULL,
                                                        "Initial points sampler",
                                                        NCM_TYPE_MSET_TRANS_KERN,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_MOVE_TYPE,
                                   g_param_spec_enum ("move-type",
                                                      NULL,
                                                      "Move type",
                                                      NCM_TYPE_FIT_ESMCMC_MOVE_TYPE, NCM_FIT_ESMCMC_MOVE_TYPE_STRETCH,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("move-scale",
                                                      NULL,
                                                      "Move scale (a)",
                                                      1.0, G_MAXDOUBLE, 2.0,
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
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void 
_ncm_fit_esmcmc_set_fit_obj (NcmFitESMCMC *esmcmc, NcmFit *fit)
{
  g_assert (esmcmc->fit == NULL);
  esmcmc->fit = ncm_fit_ref (fit);
}

/**
 * ncm_fit_esmcmc_new:
 * @fit: a #NcmFit.
 * @nwalkers: number of walkers.
 * @sampler: inital points sampler #NcmMSetTransKern.
 * @mt: move type from #NcmFitESMCMCMoveType.
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitESMCMC *
ncm_fit_esmcmc_new (NcmFit *fit, gint nwalkers, NcmMSetTransKern *sampler, NcmFitESMCMCMoveType mt, NcmFitRunMsgs mtype)
{
  NcmFitESMCMC *esmcmc = g_object_new (NCM_TYPE_FIT_ESMCMC, 
                                "fit", fit,
                                "nwalkers", nwalkers,
                                "sampler", sampler,
                                "move-type", mt,
                                "mtype", mtype,
                                NULL);
  return esmcmc;
}

/**
 * ncm_fit_esmcmc_free:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_free (NcmFitESMCMC *esmcmc)
{
  g_object_unref (esmcmc);
}

/**
 * ncm_fit_esmcmc_clear:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_clear (NcmFitESMCMC **esmcmc)
{
  g_clear_object (esmcmc);
}

/**
 * ncm_fit_esmcmc_set_data_file:
 * @esmcmc: a #NcmFitESMCMC
 * @filename: a filename.
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_set_data_file (NcmFitESMCMC *esmcmc, const gchar *filename)
{
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (esmcmc->mcat);
  
  if (esmcmc->started && cur_filename != NULL)
    g_error ("ncm_fit_esmcmc_set_data_file: Cannot change data file during a run, call ncm_fit_esmcmc_end_run() first.");    

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;

  ncm_mset_catalog_set_file (esmcmc->mcat, filename);
  
  if (esmcmc->started)
    g_assert_cmpint (esmcmc->cur_sample_id, ==, esmcmc->mcat->cur_id);
}

/**
 * ncm_fit_esmcmc_set_mtype:
 * @esmcmc: a #NcmFitESMCMC
 * @mtype: FIXME
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_set_mtype (NcmFitESMCMC *esmcmc, NcmFitRunMsgs mtype)
{
  esmcmc->mtype = mtype;
}

/**
 * ncm_fit_esmcmc_set_move_type:
 * @esmcmc: a #NcmFitESMCMC.
 * @mt: a #NcmFitESMCMCMoveType.
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_set_move_type (NcmFitESMCMC *esmcmc, NcmFitESMCMCMoveType mt)
{
  esmcmc->mt = mt;
}

/**
 * ncm_fit_esmcmc_set_move_scale:
 * @esmcmc: a #NcmFitESMCMC
 * @a: a double given the move scale (>1.0).
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_set_move_scale (NcmFitESMCMC *esmcmc, gdouble a)
{
  esmcmc->a = a;
}


/**
 * ncm_fit_esmcmc_set_trans_kern:
 * @esmcmc: a #NcmFitESMCMC
 * @tkern: a #NcmMSetTransKern.
 *
 * FIXME
 *
 */
void 
ncm_fit_esmcmc_set_sampler (NcmFitESMCMC *esmcmc, NcmMSetTransKern *tkern)
{
  ncm_mset_trans_kern_clear (&esmcmc->sampler);
  esmcmc->sampler = ncm_mset_trans_kern_ref (tkern);
}

/**
 * ncm_fit_esmcmc_set_nthreads:
 * @esmcmc: a #NcmFitESMCMC
 * @nthreads: numbers of simultaneous walkers updates.
 *
 * If @nthreads is larger than nwalkers / 2, it will be set to
 * nwalkers / 2.
 *
 */
void 
ncm_fit_esmcmc_set_nthreads (NcmFitESMCMC *esmcmc, guint nthreads)
{
  if (nthreads > 0)
  {
    if (esmcmc->nwalkers % 2 == 1)
      g_error ("ncm_fit_esmcmc_set_nthreads: cannot parallelize with an odd number of walkers [%u].", esmcmc->nwalkers);
    if (nthreads > esmcmc->nwalkers / 2)
      nthreads = esmcmc->nwalkers / 2;
  }
  esmcmc->nthreads = nthreads;
}

/**
 * ncm_fit_esmcmc_set_rng:
 * @esmcmc: a #NcmFitESMCMC
 * @rng: FIXME
 *
 * FIXME
 *
 */
void
ncm_fit_esmcmc_set_rng (NcmFitESMCMC *esmcmc, NcmRNG *rng)
{
  if (esmcmc->started)
    g_error ("ncm_fit_esmcmc_set_rng: Cannot change the RNG object during a run, call ncm_fit_esmcmc_end_run() first.");

  ncm_mset_catalog_set_rng (esmcmc->mcat, rng);
}

/**
 * ncm_fit_esmcmc_get_accept_ratio:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_fit_esmcmc_get_accept_ratio (NcmFitESMCMC *esmcmc)
{
  return esmcmc->naccepted * 1.0 / (esmcmc->ntotal * 1.0);
}

void
_ncm_fit_esmcmc_update (NcmFitESMCMC *esmcmc, NcmFit *fit, guint k)
{
  const guint part = 5;
  const guint step = esmcmc->nwalkers * ((esmcmc->n / part) == 0 ? 1 : (esmcmc->n / part));

  ncm_mset_catalog_add_from_mset (esmcmc->mcat, fit->mset, ncm_fit_state_get_m2lnL_curval (fit->fstate));
  esmcmc->cur_sample_id++;
  ncm_timer_task_increment (esmcmc->nt);

  switch (esmcmc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = (esmcmc->cur_sample_id + 1) % step;
      gboolean log_timeout = FALSE;
      if ((esmcmc->nt->pos_time - esmcmc->nt->last_log_time) > 60.0)
        log_timeout = TRUE && ((esmcmc->cur_sample_id + 1) % esmcmc->nwalkers == 0);
      if (log_timeout || (stepi == 0) || (esmcmc->nt->task_pos == esmcmc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (esmcmc->mcat);
        ncm_mset_catalog_log_current_chain_stats (esmcmc->mcat);
        g_message ("# NcmFitESMCMC:acceptance ratio %7.4f%%.\n", ncm_fit_esmcmc_get_accept_ratio (esmcmc) * 100.0);
        /* ncm_timer_task_accumulate (esmcmc->nt, acc); */
        ncm_timer_task_log_elapsed (esmcmc->nt);
        ncm_timer_task_log_mean_time (esmcmc->nt);
        ncm_timer_task_log_time_left (esmcmc->nt);
        ncm_timer_task_log_end_datetime (esmcmc->nt);
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    {
      if ((esmcmc->cur_sample_id + 1) % esmcmc->nwalkers == 0)
      {
        fit->mtype = esmcmc->mtype;
        ncm_fit_log_start (fit);
        ncm_fit_log_end (fit);
        ncm_mset_catalog_log_current_stats (esmcmc->mcat);
        ncm_mset_catalog_log_current_chain_stats (esmcmc->mcat);
        g_message ("# NcmFitESMCMC:acceptance ratio %7.4f%%.\n", ncm_fit_esmcmc_get_accept_ratio (esmcmc) * 100.0);
        /* ncm_timer_task_increment (esmcmc->nt); */
        ncm_timer_task_log_elapsed (esmcmc->nt);
        ncm_timer_task_log_mean_time (esmcmc->nt);
        ncm_timer_task_log_time_left (esmcmc->nt);
        ncm_timer_task_log_end_datetime (esmcmc->nt);
      }
      break;
    }
  }

  if ((esmcmc->mcat->fmode != NCM_MSET_CATALOG_FLUSH_TIMED) &&
      (ncm_timer_task_mean_time (esmcmc->nt) < NCM_FIT_ESMCMC_MIN_FLUSH_INTERVAL))
  {
    ncm_mset_catalog_set_flush_mode (esmcmc->mcat, NCM_MSET_CATALOG_FLUSH_TIMED);
    ncm_mset_catalog_set_flush_interval (esmcmc->mcat, NCM_FIT_ESMCMC_MIN_FLUSH_INTERVAL);
  }
}

static void ncm_fit_esmcmc_intern_skip (NcmFitESMCMC *esmcmc, guint n);

static void 
_ncm_fit_esmcmc_gen_init_points_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (data);
  glong k;
  
  for (k = i; k < f; k++)
  {
    gdouble m2lnL = 0.0;
    NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);
    NcmVector *thetastar = g_ptr_array_index (esmcmc->thetastar, k);

    while (TRUE)
    {
      g_mutex_lock (&esmcmc->resample_lock);
      ncm_mset_trans_kern_prior_sample (esmcmc->sampler, thetastar, esmcmc->mcat->rng);
      g_mutex_unlock (&esmcmc->resample_lock);

      ncm_mset_fparams_set_vector (fit_k->mset, thetastar);
      ncm_fit_m2lnL_val (fit_k, &m2lnL);

      if (gsl_finite (m2lnL))
        break;
    }
    
    ncm_fit_state_set_m2lnL_curval (fit_k->fstate, m2lnL);

    g_mutex_lock (&esmcmc->update_lock);
    while (esmcmc->write_index != k)
      g_cond_wait (&esmcmc->write_cond, &esmcmc->update_lock);

    esmcmc->ntotal++;
    esmcmc->naccepted++;
    _ncm_fit_esmcmc_update (esmcmc, fit_k, k);
    esmcmc->write_index++;
    
    g_cond_broadcast (&esmcmc->write_cond);
    g_mutex_unlock (&esmcmc->update_lock);
  }
}


static void 
_ncm_fit_esmcmc_gen_init_points (NcmFitESMCMC *esmcmc)
{
  if (esmcmc->cur_sample_id + 1 == esmcmc->nwalkers)
    return;
  else if (esmcmc->cur_sample_id + 1 > esmcmc->nwalkers)
    g_error ("_ncm_fit_esmcmc_gen_init_points: initial points already generated.");
  
  if (esmcmc->nthreads > 1)
  {
    ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_gen_init_points_mt_eval, esmcmc->cur_sample_id + 1, esmcmc->nwalkers, esmcmc);
  }
  else
  {
    gint k;

    for (k = esmcmc->cur_sample_id + 1; k < esmcmc->nwalkers; k++)
    {
      NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);
      NcmVector *thetastar = g_ptr_array_index (esmcmc->thetastar, k);
      gdouble m2lnL = 0.0;

      while (TRUE)
      {
        ncm_mset_trans_kern_prior_sample (esmcmc->sampler, thetastar, esmcmc->mcat->rng);
        ncm_mset_fparams_set_vector (fit_k->mset, thetastar);

        ncm_fit_m2lnL_val (fit_k, &m2lnL);
        if (gsl_finite (m2lnL))
          break;
      }
      ncm_fit_state_set_m2lnL_curval (fit_k->fstate, m2lnL);
      esmcmc->ntotal++;
      esmcmc->naccepted++;
      _ncm_fit_esmcmc_update (esmcmc, fit_k, k);
      esmcmc->write_index++;
    }
  }    
}

/**
 * ncm_fit_esmcmc_start_run:
 * @esmcmc: a #NcmFitESMCMC
 * 
 * FIXME
 * 
 */
void 
ncm_fit_esmcmc_start_run (NcmFitESMCMC *esmcmc)
{
  gboolean init_point_task = FALSE;
  if (esmcmc->started)
    g_error ("ncm_fit_esmcmc_start_run: run already started, run ncm_fit_esmcmc_end_run() first.");

  switch (esmcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Starting Ensamble Sampler Markov Chain Monte Carlo...\n");
      ncm_dataset_log_info (esmcmc->fit->lh->dset);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Model set:\n");
      ncm_mset_pretty_log (esmcmc->fit->mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (esmcmc->mcat->rng == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);
    ncm_rng_set_random_seed (rng, FALSE);
    ncm_fit_esmcmc_set_rng (esmcmc, rng);
    if (esmcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmFitESMCMC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));
    ncm_rng_free (rng);
  }

  esmcmc->started = TRUE;
    
  ncm_mset_catalog_sync (esmcmc->mcat, TRUE);
  esmcmc->ntotal = 0;
  esmcmc->naccepted = 0;

  if (esmcmc->mcat->first_id > 0)
    g_error ("ncm_fit_esmcmc_start_run: cannot use catalogs with first_id > 0.");

  if (esmcmc->mcat->cur_id > esmcmc->cur_sample_id)
  {
    ncm_fit_esmcmc_intern_skip (esmcmc, esmcmc->mcat->cur_id - esmcmc->cur_sample_id);
    g_assert_cmpint (esmcmc->cur_sample_id, ==, esmcmc->mcat->cur_id);
  }
  else if (esmcmc->mcat->cur_id < esmcmc->cur_sample_id)
    g_error ("ncm_fit_esmcmc_set_data_file: Unknown error cur_id < cur_sample_id [%d < %d].", 
             esmcmc->mcat->cur_id, esmcmc->cur_sample_id);

  if (esmcmc->cur_sample_id + 1 < esmcmc->nwalkers)
  {
    ncm_timer_task_start (esmcmc->nt, esmcmc->nwalkers - esmcmc->cur_sample_id - 1);
    ncm_timer_set_name (esmcmc->nt, "NcmFitESMCMC");
    init_point_task = TRUE;
  }

  if (esmcmc->cur_sample_id + 1 <= esmcmc->nwalkers)
  {
    gint k;
    for (k = 0; k <= esmcmc->cur_sample_id; k++)
    {
      NcmVector *cur_row = ncm_mset_catalog_peek_row (esmcmc->mcat, k);
      NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);
      g_assert (cur_row != NULL);
      ncm_mset_fparams_set_vector_offset (fit_k->mset, cur_row, NCM_FIT_ESMCMC_NADD_VALS);
      ncm_fit_state_set_m2lnL_curval (fit_k->fstate, ncm_vector_get (cur_row, NCM_FIT_ESMCMC_M2LNL_ID));
    }
    esmcmc->write_index = k;
    _ncm_fit_esmcmc_gen_init_points (esmcmc);
    esmcmc->ntotal    = 0;
    esmcmc->naccepted = 0;
  }
  else
  {
    guint t = (esmcmc->cur_sample_id + 1) / esmcmc->nwalkers;
    guint ki = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;
    guint tm1 = t - 1;
    gint k;

    for (k = 0; k < ki; k++)
    {
      NcmVector *cur_row = ncm_mset_catalog_peek_row (esmcmc->mcat, esmcmc->nwalkers * t + k);
      NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);
      g_assert (cur_row != NULL);
      ncm_mset_fparams_set_vector_offset (fit_k->mset, cur_row, NCM_FIT_ESMCMC_NADD_VALS);
      ncm_fit_state_set_m2lnL_curval (fit_k->fstate, ncm_vector_get (cur_row, NCM_FIT_ESMCMC_M2LNL_ID));
    }

    for (k = ki; k < esmcmc->nwalkers; k++)
    {
      NcmVector *cur_row = ncm_mset_catalog_peek_row (esmcmc->mcat, esmcmc->nwalkers * tm1 + k);
      NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);
      g_assert (cur_row != NULL);
      ncm_mset_fparams_set_vector_offset (fit_k->mset, cur_row, NCM_FIT_ESMCMC_NADD_VALS);
      ncm_fit_state_set_m2lnL_curval (fit_k->fstate, ncm_vector_get (cur_row, NCM_FIT_ESMCMC_M2LNL_ID));
    }
    esmcmc->write_index = k;
  }
  if (init_point_task)
    ncm_timer_task_pause (esmcmc->nt);
}

/**
 * ncm_fit_esmcmc_end_run:
 * @esmcmc: a #NcmFitESMCMC
 * 
 * FIXME
 * 
 */
void
ncm_fit_esmcmc_end_run (NcmFitESMCMC *esmcmc)
{
  if (ncm_timer_task_is_running (esmcmc->nt))
    ncm_timer_task_end (esmcmc->nt);

  ncm_mset_catalog_sync (esmcmc->mcat, TRUE);
  
  esmcmc->started = FALSE;
}

/**
 * ncm_fit_esmcmc_reset:
 * @esmcmc: a #NcmFitESMCMC
 * 
 * FIXME
 * 
 */
void 
ncm_fit_esmcmc_reset (NcmFitESMCMC *esmcmc)
{
  esmcmc->n               = 0;
  esmcmc->cur_sample_id   = -1;
  esmcmc->ntotal          = 0;
  esmcmc->naccepted       = 0;
  esmcmc->write_index     = 0;
  esmcmc->started         = FALSE;  
  ncm_mset_catalog_reset (esmcmc->mcat);
}

static void 
ncm_fit_esmcmc_intern_skip (NcmFitESMCMC *esmcmc, guint n)
{
  if (n == 0)
    return;

  switch (esmcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Skipping %u points (%f iterations), will start at %u-th point.\n", n, n * 1.0 / esmcmc->nwalkers, esmcmc->cur_sample_id + n + 1 + 1);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  esmcmc->cur_sample_id += n;
  esmcmc->write_index = esmcmc->cur_sample_id + 1;
}

static void _ncm_fit_esmcmc_run_single (NcmFitESMCMC *esmcmc);
static void _ncm_fit_esmcmc_run_mt (NcmFitESMCMC *esmcmc);

/**
 * ncm_fit_esmcmc_run:
 * @esmcmc: a #NcmFitESMCMC
 * @n: total number of realizations to run
 * 
 * Runs the Monte Carlo until it reaches the @n-th realization. Note that
 * if the first_id is non-zero it will run @n - first_id realizations.
 *
 */
void 
ncm_fit_esmcmc_run (NcmFitESMCMC *esmcmc, guint n)
{
  guint ti = (esmcmc->cur_sample_id + 1) / esmcmc->nwalkers;
  guint ki = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;
  
  if (!esmcmc->started)
    g_error ("ncm_fit_esmcmc_run: run not started, run ncm_fit_esmcmc_start_run() first.");

  if (n <= ti)
  {
    if (esmcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Nothing to do, current Monte Carlo run is %d\n", ti);
    }
    return;
  }
  
  esmcmc->n = n - ti;
  
  switch (esmcmc->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      switch (esmcmc->mt)
      {
        case NCM_FIT_ESMCMC_MOVE_TYPE_STRETCH:
        g_message ("# NcmFitESMCMC: Calculating [%06d] Ensemble Sampler Markov Chain Monte Carlo runs [Stretch Move]\n", esmcmc->n);
          break;
        default:
          g_assert_not_reached ();
      }
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_timer_task_is_running (esmcmc->nt))
  {
    ncm_timer_task_add_tasks (esmcmc->nt, esmcmc->n * esmcmc->nwalkers - ki);
    ncm_timer_task_continue (esmcmc->nt);
  }
  else
  {
    ncm_timer_task_start (esmcmc->nt, esmcmc->n * esmcmc->nwalkers - ki);
    ncm_timer_set_name (esmcmc->nt, "NcmFitESMCMC");
  }
  if (esmcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (esmcmc->nt);

  if (esmcmc->nthreads <= 1)
    _ncm_fit_esmcmc_run_single (esmcmc);
  else
    _ncm_fit_esmcmc_run_mt (esmcmc);

  ncm_timer_task_pause (esmcmc->nt);
}

static void
ncm_fit_esmcmc_move_try (NcmFitESMCMC *esmcmc, gdouble z, NcmVector *theta_k, NcmVector *theta_j, NcmVector *thetastar)
{
  const guint len = ncm_vector_len (theta_k);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble theta_j_i = ncm_vector_get (theta_j, i);
    const gdouble theta_k_i = ncm_vector_get (theta_k, i);
    
    ncm_vector_set (thetastar, i, theta_j_i + z * (theta_k_i - theta_j_i));
  }
}

static void 
_ncm_fit_esmcmc_run_single (NcmFitESMCMC *esmcmc)
{
  guint i = 0, k;
  guint ki = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;

  for (i = 0; i < esmcmc->n; i++)
  {
    k = ki;
    while (k < esmcmc->nwalkers)
    {
      NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);

      NcmVector *theta_k = g_ptr_array_index (esmcmc->theta_k, k);
      NcmVector *theta_j = g_ptr_array_index (esmcmc->theta_j, k);
      NcmVector *thetastar = g_ptr_array_index (esmcmc->thetastar, k);

      gdouble m2lnL_cur = ncm_fit_state_get_m2lnL_curval (fit_k->fstate);
      gdouble m2lnL_star, prob, jump, z, u;
      gulong j;

      jump = gsl_rng_uniform (esmcmc->mcat->rng->r);
      u    = gsl_rng_uniform (esmcmc->mcat->rng->r);
      z    = gsl_pow_2 (1.0 + (esmcmc->a - 1.0) * u) / esmcmc->a;
      j    = gsl_rng_uniform_int (esmcmc->mcat->rng->r, esmcmc->nwalkers - 1);

      if (j >= k)
        j++;

      ncm_mset_fparams_get_vector (fit_k->mset, theta_k);
      ncm_mset_fparams_get_vector (NCM_FIT (g_ptr_array_index (esmcmc->walker_fits, j))->mset, theta_j);

      ncm_fit_esmcmc_move_try (esmcmc, z, theta_k, theta_j, thetastar);
/*
      ncm_vector_log_vals (theta_k, "# theta_k  ", "% 20.15g");
      ncm_vector_log_vals (theta_j, "# theta_j  ", "% 20.15g");
      ncm_vector_log_vals (thetastar, "# thetastar", "% 20.15g");
*/
      if (!ncm_mset_fparam_valid_bounds (fit_k->mset, thetastar))
        continue;
      
      ncm_mset_fparams_set_vector (fit_k->mset, thetastar);

      if (!ncm_mset_params_valid_bounds (fit_k->mset))

      ncm_fit_m2lnL_val (fit_k, &m2lnL_star);
      esmcmc->ntotal++;
      esmcmc->naccepted++;
      /*
       ncm_vector_log_vals (esmcmc->theta, "# Theta  : ", "% 8.5g");
       ncm_vector_log_vals (esmcmc->thetastar, "# Theta* : ", "% 8.5g");
       */
      if (gsl_finite (m2lnL_star))
      {
        prob = pow (z, esmcmc->fparam_len - 1.0) * exp ((m2lnL_cur - m2lnL_star) * 0.5);
        prob = GSL_MIN (prob, 1.0);
        ncm_fit_state_set_m2lnL_curval (fit_k->fstate, m2lnL_star);
      }
      else
        prob = 0.0;

      /*printf ("# Prob %e [% 21.16g % 21.16g] % 21.16g\n", prob, m2lnL_cur, m2lnL_star, m2lnL_cur - m2lnL_star);*/    

      if (prob != 1.0)
      {
        if (jump > prob)
        {
          ncm_mset_fparams_set_vector (fit_k->mset, theta_k);
          ncm_fit_state_set_m2lnL_curval (fit_k->fstate, m2lnL_cur);
          esmcmc->naccepted--;
        }
      }

      _ncm_fit_esmcmc_update (esmcmc, fit_k, k);
      esmcmc->write_index++;
      k++;
    }
    ki = 0;
  }
}

static void 
_ncm_fit_esmcmc_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (data);
  const guint nwalkers_2  = esmcmc->nwalkers / 2;
  const guint subensemble = (i < nwalkers_2) ? nwalkers_2 : 0;
  guint k = i;

  while (k < f)
  {
    NcmFit *fit_k = g_ptr_array_index (esmcmc->walker_fits, k);
    NcmVector *theta_k = g_ptr_array_index (esmcmc->theta_k, k);
    NcmVector *theta_j = g_ptr_array_index (esmcmc->theta_j, k);
    NcmVector *thetastar = g_ptr_array_index (esmcmc->thetastar, k);

    gdouble m2lnL_cur = ncm_fit_state_get_m2lnL_curval (fit_k->fstate);
    gdouble m2lnL_star, prob, jump, z, u;
    gulong j;
    gboolean accepted = TRUE;

    g_mutex_lock (&esmcmc->resample_lock);
    jump = gsl_rng_uniform (esmcmc->mcat->rng->r);
    u    = gsl_rng_uniform (esmcmc->mcat->rng->r);
    j    = gsl_rng_uniform_int (esmcmc->mcat->rng->r, nwalkers_2);
    g_mutex_unlock (&esmcmc->resample_lock);

    z    = gsl_pow_2 (1.0 + (esmcmc->a - 1.0) * u) / esmcmc->a;
    j   += subensemble;

    ncm_mset_fparams_get_vector (fit_k->mset, theta_k);
    ncm_mset_fparams_get_vector (NCM_FIT (g_ptr_array_index (esmcmc->walker_fits, j))->mset, theta_j);

    ncm_fit_esmcmc_move_try (esmcmc, z, theta_k, theta_j, thetastar);
    if (!ncm_mset_fparam_valid_bounds (fit_k->mset, thetastar))
      continue;

    ncm_mset_fparams_set_vector (fit_k->mset, thetastar);

    ncm_fit_m2lnL_val (fit_k, &m2lnL_star);
    /*
     ncm_vector_log_vals (esmcmc->theta, "# Theta  : ", "% 8.5g");
     ncm_vector_log_vals (esmcmc->thetastar, "# Theta* : ", "% 8.5g");
     */
    if (gsl_finite (m2lnL_star))
    {
      prob = pow (z, esmcmc->fparam_len - 1.0) * exp ((m2lnL_cur - m2lnL_star) * 0.5);
      prob = GSL_MIN (prob, 1.0);
      ncm_fit_state_set_m2lnL_curval (fit_k->fstate, m2lnL_star);
    }
    else
      prob = 0.0;
    
    /*printf ("# Prob %e [% 21.16g % 21.16g] % 21.16g\n", prob, m2lnL_cur, m2lnL_star, m2lnL_cur - m2lnL_star);*/    

    if (prob != 1.0)
    {
      if (jump > prob)
      {
        ncm_mset_fparams_set_vector (fit_k->mset, theta_k);
        ncm_fit_state_set_m2lnL_curval (fit_k->fstate, m2lnL_cur);
        accepted = FALSE;
      }
    }
    
    g_mutex_lock (&esmcmc->update_lock);    
    while (esmcmc->write_index != k)
      g_cond_wait (&esmcmc->write_cond, &esmcmc->update_lock);

    esmcmc->ntotal++;
    if (accepted)
      esmcmc->naccepted++;

    _ncm_fit_esmcmc_update (esmcmc, fit_k, k);
    esmcmc->write_index++;
    g_cond_broadcast (&esmcmc->write_cond);

    g_mutex_unlock (&esmcmc->update_lock);
    k++;
  }
}

static void
_ncm_fit_esmcmc_run_mt (NcmFitESMCMC *esmcmc)
{
  guint i;

  if (esmcmc->nthreads <= 1)
  {
    _ncm_fit_esmcmc_run_single (esmcmc);
    return;
  }

  g_assert_cmpuint (esmcmc->nthreads, >, 1);

  {
    const guint nwalkers_2 = esmcmc->nwalkers / 2;
    guint ki = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;

    if (esmcmc->n > 0)
    {
      esmcmc->write_index = ki;
      if (ki < nwalkers_2)
      {
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, ki, nwalkers_2, esmcmc);
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, nwalkers_2, esmcmc->nwalkers, esmcmc);
      }
      else
      {
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, ki, esmcmc->nwalkers, esmcmc);
      }

      for (i = 1; i < esmcmc->n; i++)
      {
        esmcmc->write_index = 0;
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, 0, nwalkers_2, esmcmc);
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, nwalkers_2, esmcmc->nwalkers, esmcmc);
      }
    }
  }
}

/**
 * ncm_fit_esmcmc_run_lre:
 * @esmcmc: a #NcmFitESMCMC
 * @prerun: FIXME
 * @lre: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_fit_esmcmc_run_lre (NcmFitESMCMC *esmcmc, guint prerun, gdouble lre)
{
  gdouble lerror;
  const gdouble lre2 = lre * lre;
  const guint fparam_len = ncm_mset_fparam_len (esmcmc->fit->mset);
  const guint catlen = ncm_mset_catalog_len (esmcmc->mcat) / esmcmc->nwalkers;

  g_assert_cmpfloat (lre, >, 0.0);

  prerun = GSL_MAX (prerun, fparam_len * 100);

  if (catlen < prerun)
  {
    guint prerun_left = prerun - catlen;
    if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitESMCMC: Running first %u pre-runs...\n", prerun_left);
    ncm_fit_esmcmc_run (esmcmc, prerun);
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (esmcmc->mcat);
  lerror = ncm_mset_catalog_largest_error (esmcmc->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n = ncm_mset_catalog_len (esmcmc->mcat);
    gdouble m = n * lerror2 / lre2;
    guint runs = ((m - n) > 1000.0) ? ceil ((m - n) * 1.0e-1) : ceil (m - n);
    guint ti = (esmcmc->cur_sample_id + 1) / esmcmc->nwalkers;
    
    runs = runs / esmcmc->nwalkers + 1;

    if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmFitESMCMC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmFitESMCMC: Running more %u runs...\n", runs);
    }
    ncm_fit_esmcmc_run (esmcmc, ti + runs);
    ncm_mset_catalog_estimate_autocorrelation_tau (esmcmc->mcat);
    lerror = ncm_mset_catalog_largest_error (esmcmc->mcat);
  }

  if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    g_message ("# NcmFitESMCMC: Largest relative error %e attained: %e\n", lre, lerror);
}

/**
 * ncm_fit_esmcmc_mean_covar:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 */
void
ncm_fit_esmcmc_mean_covar (NcmFitESMCMC *esmcmc)
{
  ncm_mset_catalog_get_mean (esmcmc->mcat, &esmcmc->fit->fstate->fparams);
  ncm_mset_catalog_get_covar (esmcmc->mcat, &esmcmc->fit->fstate->covar);
  ncm_mset_fparams_set_vector (esmcmc->mcat->mset, esmcmc->fit->fstate->fparams);
  esmcmc->fit->fstate->has_covar = TRUE;
}

/**
 * ncm_fit_esmcmc_get_catalog:
 * @esmcmc: a #NcmFitESMCMC
 *
 * Gets the generated catalog of @esmcmc.
 * 
 * Returns: (transfer full): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_esmcmc_get_catalog (NcmFitESMCMC *esmcmc)
{
  return ncm_mset_catalog_ref (esmcmc->mcat);
}

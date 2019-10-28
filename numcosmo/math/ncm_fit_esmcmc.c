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
#include "math/ncm_timer.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_mpi_job_mcmc.h"
#include "math/ncm_c.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmFitESMCMCPrivate
{
  NcmFit *fit;
	NcmMPIJob *mj;
  NcmMemoryPool *walker_pool;
  NcmMSetTransKern *sampler;
  NcmMSetCatalog *mcat;
  NcmFitRunMsgs mtype;
  NcmTimer *nt;
  NcmSerialize *ser;
  NcmFitESMCMCWalker *walker;
  gboolean auto_trim;
  guint auto_trim_div;
  gdouble lre_step;
  NcmMSetCatalogTrimType trim_type;
  guint min_runs;
  gdouble max_runs_time;
  gdouble log_time_interval;
  guint interm_log;
  GPtrArray *full_theta;
  GPtrArray *full_thetastar;
  GPtrArray *full_thetastar_inout;
  GPtrArray *m2lnL;
  GPtrArray *theta;
  GPtrArray *thetastar;
  GPtrArray *thetastar_in;
  GPtrArray *thetastar_out;
  NcmVector *jumps;
  GArray *accepted;
  GArray *offboard;
  NcmObjArray *func_oa;
  gchar *func_oa_file;
  guint nadd_vals;
  guint fparam_len;
  guint nthreads;
  gboolean use_mpi;
  gboolean has_mpi;
  guint nslaves;
  guint n;
  gint nwalkers;
  gint cur_sample_id;
  guint ntotal;
  guint naccepted;
  guint noffboard;
  guint ntotal_lup;
  guint naccepted_lup;
  guint noffboard_lup;
  gboolean started;
  GMutex dup_fit;
  GMutex resample_lock;
  GMutex update_lock;
  GCond write_cond;
};

enum
{
  PROP_0,
  PROP_FIT,
  PROP_NWALKERS,
  PROP_SAMPLER,
  PROP_WALKER,
  PROP_LRE_STEP,
  PROP_AUTO_TRIM,
  PROP_AUTO_TRIM_DIV,
  PROP_TRIM_TYPE,
  PROP_MIN_RUNS,
  PROP_MAX_RUNS_TIME,
  PROP_LOG_TIME_INTERVAL,
  PROP_INTERM_LOG,
  PROP_MTYPE,
  PROP_NTHREADS,
  PROP_USE_MPI,
  PROP_DATA_FILE,
  PROP_FUNC_ARRAY,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFitESMCMC, ncm_fit_esmcmc, G_TYPE_OBJECT);

static gpointer _ncm_fit_esmcmc_worker_dup (gpointer userdata);
static void _ncm_fit_esmcmc_worker_free (gpointer p);

static void
ncm_fit_esmcmc_init (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv = ncm_fit_esmcmc_get_instance_private (esmcmc);
	
  self->fit               = NULL;
	self->mj                = NULL;
  self->walker_pool       = ncm_memory_pool_new (&_ncm_fit_esmcmc_worker_dup, esmcmc, 
                                                 &_ncm_fit_esmcmc_worker_free);

  self->sampler           = NULL;
  self->mcat              = NULL;
  self->mtype             = NCM_FIT_RUN_MSGS_NONE;
  self->nt                = ncm_timer_new ();
  self->ser               = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  self->walker            = NULL;
  self->lre_step          = 0.0;
  self->auto_trim         = FALSE;
  self->auto_trim_div     = 0;
  self->trim_type         = 0;
  self->min_runs          = 0;
  self->max_runs_time     = 0.0;
  self->log_time_interval = 0.0;
  self->nadd_vals         = 0;
  self->fparam_len        = 0;

  self->m2lnL                = g_ptr_array_new ();
  self->theta                = g_ptr_array_new ();
  self->thetastar            = g_ptr_array_new ();
  self->full_theta           = g_ptr_array_new ();
  self->full_thetastar       = g_ptr_array_new ();
  self->full_thetastar_inout = g_ptr_array_new ();
  self->thetastar_in         = g_ptr_array_new ();
  self->thetastar_out        = g_ptr_array_new ();

  g_ptr_array_set_free_func (self->m2lnL, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (self->theta, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (self->thetastar, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (self->thetastar_in, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (self->thetastar_out, (GDestroyNotify) &ncm_vector_free);

  g_ptr_array_set_free_func (self->full_theta, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (self->full_thetastar, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (self->full_thetastar_inout, (GDestroyNotify) &ncm_vector_free);

  self->jumps           = NULL;
  self->accepted        = g_array_new (TRUE, TRUE, sizeof (gboolean));
  self->offboard        = g_array_new (TRUE, TRUE, sizeof (gboolean));

  self->func_oa         = NULL;
  self->func_oa_file    = NULL;
  
  self->nthreads        = 0;
  self->use_mpi         = FALSE;
  self->has_mpi         = FALSE;
  self->nslaves         = 0;
  self->n               = 0;
  self->nwalkers        = 0;
  self->cur_sample_id   = -1; /* Represents that no samples were calculated yet, i.e., id of the last added point. */
  self->ntotal          = 0;
  self->naccepted       = 0;
  self->noffboard       = 0;
  self->ntotal_lup      = 0;
  self->naccepted_lup   = 0;
  self->noffboard_lup   = 0;
  self->started         = FALSE;

  g_mutex_init (&self->dup_fit);
  g_mutex_init (&self->resample_lock);
  g_mutex_init (&self->update_lock);
  g_cond_init (&self->write_cond);
}

static void
_ncm_fit_esmcmc_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_esmcmc_parent_class)->constructed (object);
  {
    NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
		NcmFitESMCMCPrivate * const self = esmcmc->priv;
    guint k;

    g_assert_cmpint (self->nwalkers, >, 0);
    self->fparam_len = ncm_mset_fparam_len (self->fit->mset);
    {
      const guint nfuncs    = (self->func_oa != NULL) ? self->func_oa->len : 0;
      const guint nadd_vals = self->nadd_vals = nfuncs + 1;
      const guint theta_len = self->fparam_len + nadd_vals;

      if (nfuncs > 0)
      {
        gchar **names   = g_new (gchar *, nadd_vals + 1);
        gchar **symbols = g_new (gchar *, nadd_vals + 1);

        names[0]   = g_strdup (NCM_MSET_CATALOG_M2LNL_COLNAME);
        symbols[0] = g_strdup (NCM_MSET_CATALOG_M2LNL_SYMBOL);

        for (k = 0; k < nfuncs; k++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, k));
          g_assert (NCM_IS_MSET_FUNC (func));

          names[1 + k]   = g_strdup (ncm_mset_func_peek_uname (func));
          symbols[1 + k] = g_strdup (ncm_mset_func_peek_usymbol (func));
        }
        names[1 + k]   = NULL;
        symbols[1 + k] = NULL;

        self->mcat = ncm_mset_catalog_new_array (self->fit->mset, nadd_vals, self->nwalkers, FALSE, 
                                                   names, symbols);

        g_strfreev (names);
        g_strfreev (symbols);
      }
      else
      {
				self->mcat = ncm_mset_catalog_new (self->fit->mset, nadd_vals, self->nwalkers, FALSE, 
				                                   NCM_MSET_CATALOG_M2LNL_COLNAME, NCM_MSET_CATALOG_M2LNL_SYMBOL, 
				                                   NULL);
			}
			
      ncm_mset_catalog_set_m2lnp_var (self->mcat, 0);
      ncm_mset_catalog_set_run_type (self->mcat, "Ensemble Sampler MCMC");

      for (k = 0; k < self->nwalkers; k++)
      {
				NcmVector *full_thetastar_inout_k = ncm_vector_new (theta_len + NCM_FIT_ESMCMC_MPI_IN_LEN + NCM_FIT_ESMCMC_MPI_OUT_LEN);
        NcmVector *full_theta_k           = ncm_vector_new (theta_len);
        NcmVector *full_thetastar_k       = ncm_vector_get_subvector (full_thetastar_inout_k, NCM_FIT_ESMCMC_MPI_OUT_LEN, theta_len);
        NcmVector *m2lnL_k                = ncm_vector_get_subvector (full_theta_k, NCM_FIT_ESMCMC_M2LNL_ID, 1);
        NcmVector *theta_k                = ncm_vector_get_subvector (full_theta_k, nadd_vals, self->fparam_len);
        NcmVector *thetastar_k            = ncm_vector_get_subvector (full_thetastar_k, nadd_vals, self->fparam_len);
        NcmVector *thetastar_in_k         = ncm_vector_get_subvector (full_thetastar_inout_k, nadd_vals + NCM_FIT_ESMCMC_MPI_OUT_LEN, self->fparam_len + NCM_FIT_ESMCMC_MPI_OUT_LEN);
        NcmVector *thetastar_out_k        = ncm_vector_get_subvector (full_thetastar_inout_k, 0, nadd_vals + NCM_FIT_ESMCMC_MPI_OUT_LEN);

        g_ptr_array_add (self->m2lnL,     m2lnL_k);
        g_ptr_array_add (self->theta,     theta_k);
        g_ptr_array_add (self->thetastar, thetastar_k);

				g_ptr_array_add (self->thetastar_in,  thetastar_in_k);
        g_ptr_array_add (self->thetastar_out, thetastar_out_k);

        g_ptr_array_add (self->full_theta,           full_theta_k);
        g_ptr_array_add (self->full_thetastar,       full_thetastar_k);
        g_ptr_array_add (self->full_thetastar_inout, full_thetastar_inout_k);
      }
    }

    self->jumps = ncm_vector_new (self->nwalkers);
    g_array_set_size (self->accepted, self->nwalkers);
    g_array_set_size (self->offboard, self->nwalkers);
    
    if (self->walker == NULL)
      self->walker = ncm_fit_esmcmc_walker_new_from_name ("NcmFitESMCMCWalkerStretch");

    ncm_fit_esmcmc_walker_set_size (self->walker, self->nwalkers);
    ncm_fit_esmcmc_walker_set_nparams (self->walker, self->fparam_len);

		g_assert (self->mj == NULL);
		self->mj = NCM_MPI_JOB (ncm_mpi_job_mcmc_new (self->fit, self->func_oa));
  }
}

static void _ncm_fit_esmcmc_set_fit_obj (NcmFitESMCMC *esmcmc, NcmFit *fit);

static void
_ncm_fit_esmcmc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  g_return_if_fail (NCM_IS_FIT_ESMCMC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      _ncm_fit_esmcmc_set_fit_obj (esmcmc, g_value_get_object (value));
      break;
    case PROP_NWALKERS:
      self->nwalkers = g_value_get_int (value);
      break;      
    case PROP_SAMPLER:
      ncm_fit_esmcmc_set_sampler (esmcmc, g_value_get_object (value));
      break;      
    case PROP_WALKER:
      self->walker = g_value_dup_object (value);
      break;
    case PROP_LRE_STEP:
      self->lre_step = g_value_get_double (value);
      break;
    case PROP_AUTO_TRIM:
      ncm_fit_esmcmc_set_auto_trim (esmcmc, g_value_get_boolean (value));
      break;
    case PROP_AUTO_TRIM_DIV:
      ncm_fit_esmcmc_set_auto_trim_div (esmcmc, g_value_get_uint (value));
      break;
    case PROP_TRIM_TYPE:
      ncm_fit_esmcmc_set_auto_trim_type (esmcmc, g_value_get_flags (value));
      break;
    case PROP_MIN_RUNS:
      self->min_runs = g_value_get_uint (value);
      break;
    case PROP_MAX_RUNS_TIME:
      self->max_runs_time = g_value_get_double (value);
      break;
    case PROP_LOG_TIME_INTERVAL:
      self->log_time_interval = g_value_get_double (value);
      break;
    case PROP_INTERM_LOG:
      self->interm_log = g_value_get_uint (value);
      break;
    case PROP_MTYPE:
      ncm_fit_esmcmc_set_mtype (esmcmc, g_value_get_enum (value));
      break;
    case PROP_NTHREADS:
      ncm_fit_esmcmc_set_nthreads (esmcmc, g_value_get_uint (value));
      break;
    case PROP_USE_MPI:
      ncm_fit_esmcmc_use_mpi (esmcmc, g_value_get_boolean (value));
      break;
    case PROP_DATA_FILE:
      ncm_fit_esmcmc_set_data_file (esmcmc, g_value_get_string (value));
      break;
    case PROP_FUNC_ARRAY:
    {
			ncm_obj_array_clear (&self->func_oa);
      self->func_oa = g_value_dup_boxed (value);
      if (self->func_oa != NULL)
      {
        guint i;
        for (i = 0; i < self->func_oa->len; i++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (self->func_oa, i));
          g_assert (NCM_IS_MSET_FUNC (func));
          g_assert (ncm_mset_func_is_scalar (func));
          g_assert (ncm_mset_func_is_const (func));
        }
      }
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  g_return_if_fail (NCM_IS_FIT_ESMCMC (object));

  switch (prop_id)
  {
    case PROP_FIT:
      g_value_set_object (value, self->fit);
      break;
    case PROP_NWALKERS:
      g_value_set_int (value, self->nwalkers);
      break;
    case PROP_SAMPLER:
      g_value_set_object (value, self->sampler);
      break;      
    case PROP_WALKER:
      g_value_set_object (value, self->walker);
      break;
    case PROP_LRE_STEP:
      g_value_set_double (value, self->lre_step);
      break;
    case PROP_AUTO_TRIM:
      g_value_set_boolean (value, self->auto_trim);
      break;
    case PROP_AUTO_TRIM_DIV:
      g_value_set_uint (value, self->auto_trim_div);
      break;
    case PROP_TRIM_TYPE:
      g_value_set_flags (value, self->trim_type);
      break;
    case PROP_MIN_RUNS:
      g_value_set_uint (value, self->min_runs);
      break;
    case PROP_MAX_RUNS_TIME:
      g_value_set_double (value, self->max_runs_time);
      break;
    case PROP_LOG_TIME_INTERVAL:
      g_value_set_double (value, self->log_time_interval);
      break;
    case PROP_INTERM_LOG:
      g_value_set_uint (value, self->interm_log);
      break;
    case PROP_MTYPE:
      g_value_set_enum (value, self->mtype);
      break;
    case PROP_NTHREADS:
      g_value_set_uint (value, self->nthreads);
      break;
    case PROP_USE_MPI:
      g_value_set_uint (value, self->use_mpi);
      break;
    case PROP_DATA_FILE:
      g_value_set_string (value, ncm_mset_catalog_peek_filename (self->mcat));
      break;
    case PROP_FUNC_ARRAY:
      g_value_set_boxed (value, self->func_oa);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_dispose (GObject *object)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;

  ncm_fit_clear (&self->fit);
	ncm_mpi_job_clear (&self->mj);
	
  ncm_mset_trans_kern_clear (&self->sampler);
  ncm_timer_clear (&self->nt);
  ncm_serialize_clear (&self->ser);
  ncm_mset_catalog_clear (&self->mcat);

  ncm_fit_esmcmc_walker_clear (&self->walker);

  ncm_vector_clear (&self->jumps);

  ncm_obj_array_clear (&self->func_oa);

  g_clear_pointer (&self->func_oa_file, g_free);
  
  g_clear_pointer (&self->m2lnL, g_ptr_array_unref);
  g_clear_pointer (&self->theta, g_ptr_array_unref);
  g_clear_pointer (&self->thetastar, g_ptr_array_unref);
  g_clear_pointer (&self->thetastar_in, g_ptr_array_unref);
  g_clear_pointer (&self->thetastar_out, g_ptr_array_unref);

  g_clear_pointer (&self->full_theta, g_ptr_array_unref);
  g_clear_pointer (&self->full_thetastar, g_ptr_array_unref);
  g_clear_pointer (&self->full_thetastar_inout, g_ptr_array_unref);

  g_clear_pointer (&self->accepted, g_array_unref);
  g_clear_pointer (&self->offboard, g_array_unref);

  if (self->walker_pool != NULL)
  {
    ncm_memory_pool_free (self->walker_pool, TRUE);
    self->walker_pool = NULL;
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_finalize (GObject *object)
{
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;

  g_mutex_clear (&self->dup_fit);
  g_mutex_clear (&self->resample_lock);
  g_mutex_clear (&self->update_lock);
  g_cond_clear (&self->write_cond);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_parent_class)->finalize (object);
}

static void
ncm_fit_esmcmc_class_init (NcmFitESMCMCClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_esmcmc_constructed;
  object_class->set_property = &_ncm_fit_esmcmc_set_property;
  object_class->get_property = &_ncm_fit_esmcmc_get_property;
  object_class->dispose      = &_ncm_fit_esmcmc_dispose;
  object_class->finalize     = &_ncm_fit_esmcmc_finalize;

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
                                   PROP_WALKER,
                                   g_param_spec_object ("walker",
                                                        NULL,
                                                        "Walker object",
                                                        NCM_TYPE_FIT_ESMCMC_WALKER,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LRE_STEP,
                                   g_param_spec_double ("lre-step",
                                                         NULL,
                                                         "Step size in the lre run",
                                                         1.0e-2, 1.0, 0.1,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_AUTO_TRIM,
                                   g_param_spec_boolean ("auto-trim",
                                                         NULL,
                                                         "Whether to automatically trim the catalog",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_AUTO_TRIM_DIV,
                                   g_param_spec_uint ("auto-trim-div",
                                                      NULL,
                                                      "Automatically trim divisor",
                                                      1, G_MAXUINT, 100,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRIM_TYPE,
                                   g_param_spec_flags ("trim-type",
                                                       NULL,
                                                       "Trimming tests to apply",
                                                       NCM_TYPE_MSET_CATALOG_TRIM_TYPE, NCM_MSET_CATALOG_TRIM_TYPE_ESS,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MIN_RUNS,
                                   g_param_spec_uint ("min-runs",
                                                      NULL,
                                                      "Minumum number of runs",
                                                      1, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MAX_RUNS_TIME,
                                   g_param_spec_double ("max-runs-time",
                                                        NULL,
                                                        "Maximum time between runs",
                                                        1.0, G_MAXDOUBLE, 2.0 * 60.0 * 60.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LOG_TIME_INTERVAL,
                                   g_param_spec_double ("log-time-interval",
                                                        NULL,
                                                        "Time interval between log",
                                                        1.0, G_MAXDOUBLE, 60.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INTERM_LOG,
                                   g_param_spec_uint ("intermediary-log",
                                                      NULL,
                                                      "Number of intermediary logs",
                                                      1, G_MAXUINT, 5,
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
	g_object_class_install_property (object_class,
	                                 PROP_USE_MPI,
	                                 g_param_spec_boolean ("use-mpi",
	                                                       NULL,
	                                                       "Use MPI instead of threads",
	                                                       TRUE,
	                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT  | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
                                   PROP_DATA_FILE,
                                   g_param_spec_string ("data-file",
                                                        NULL,
                                                        "Data filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FUNC_ARRAY,
                                   g_param_spec_boxed ("function-array",
                                                       NULL,
                                                       "Functions array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

typedef struct _NcmFitESMCMCWorker
{
  NcmFit *fit;
  NcmObjArray *funcs_array;
} NcmFitESMCMCWorker;

static gpointer
_ncm_fit_esmcmc_worker_dup (gpointer userdata)
{
  G_LOCK_DEFINE_STATIC (dup_thread);
  NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (userdata);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;

  G_LOCK (dup_thread);
  {
    NcmFitESMCMCWorker *fw = g_new (NcmFitESMCMCWorker, 1);

    fw->fit = ncm_fit_dup (self->fit, self->ser);
		
    if (self->func_oa != NULL)
      fw->funcs_array = ncm_obj_array_dup (self->func_oa, self->ser);    
    else
      fw->funcs_array = NULL;

    ncm_serialize_reset (self->ser, TRUE);
    
    G_UNLOCK (dup_thread);

    return fw;
  }
}

static void
_ncm_fit_esmcmc_worker_free (gpointer userdata)
{
  NcmFitESMCMCWorker *fw = (NcmFitESMCMCWorker *) userdata;

  ncm_fit_clear (&fw->fit);
  ncm_obj_array_clear (&fw->funcs_array);

  g_free (fw);
}

static void 
_ncm_fit_esmcmc_set_fit_obj (NcmFitESMCMC *esmcmc, NcmFit *fit)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  g_assert (self->fit == NULL);
  self->fit = ncm_fit_ref (fit);
}

/**
 * ncm_fit_esmcmc_new:
 * @fit: a #NcmFit
 * @nwalkers: number of walkers
 * @sampler: inital points sampler #NcmMSetTransKern
 * @walker: (allow-none): a #NcmFitESMCMCWalker
 * @mtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitESMCMC *
ncm_fit_esmcmc_new (NcmFit *fit, gint nwalkers, NcmMSetTransKern *sampler, NcmFitESMCMCWalker *walker, NcmFitRunMsgs mtype)
{
  NcmFitESMCMC *esmcmc = g_object_new (NCM_TYPE_FIT_ESMCMC, 
                                       "fit", fit,
                                       "nwalkers", nwalkers,
                                       "sampler", sampler,
                                       "walker", walker,
                                       "mtype", mtype,
                                       NULL);
  return esmcmc;
}

/**
 * ncm_fit_esmcmc_new_funcs_array:
 * @fit: a #NcmFit
 * @nwalkers: number of walkers
 * @sampler: inital points sampler #NcmMSetTransKern
 * @walker: (allow-none): a #NcmFitESMCMCWalker
 * @mtype: FIXME
 * @funcs_array: a #NcmObjArray of scalar functions to include in the catalog.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFitESMCMC *
ncm_fit_esmcmc_new_funcs_array (NcmFit *fit, gint nwalkers, NcmMSetTransKern *sampler, NcmFitESMCMCWalker *walker, NcmFitRunMsgs mtype, NcmObjArray *funcs_array)
{
  NcmFitESMCMC *esmcmc = g_object_new (NCM_TYPE_FIT_ESMCMC, 
                                       "fit", fit,
                                       "nwalkers",       nwalkers,
                                       "sampler",        sampler,
                                       "walker",         walker,
                                       "mtype",          mtype,
                                       "function-array", funcs_array,
                                       NULL);
  return esmcmc;
}


/**
 * ncm_fit_esmcmc_ref:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmFitESMCMC *
ncm_fit_esmcmc_ref (NcmFitESMCMC *esmcmc)
{
  return g_object_ref (esmcmc);
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const gchar *cur_filename = ncm_mset_catalog_peek_filename (self->mcat);
  
  if (self->started && cur_filename != NULL)
    g_error ("ncm_fit_esmcmc_set_data_file: Cannot change data file during a run, call ncm_fit_esmcmc_end_run() first.");    

  if (cur_filename != NULL && strcmp (cur_filename, filename) == 0)
    return;

  ncm_mset_catalog_set_file (self->mcat, filename);

  if (self->func_oa != NULL)
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    gchar *base_name      = ncm_util_basename_fits (filename);
    self->func_oa_file = g_strdup_printf ("%s.oa", base_name);
    g_free (base_name); 

    ncm_obj_array_save (self->func_oa, ser, self->func_oa_file, TRUE);

    ncm_serialize_free (ser);
  }
  
  if (self->started)
    g_assert_cmpint (self->cur_sample_id, ==, ncm_mset_catalog_get_cur_id (self->mcat));
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  self->mtype = mtype;
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  ncm_mset_trans_kern_clear (&self->sampler);
  self->sampler = ncm_mset_trans_kern_ref (tkern);
}

/**
 * ncm_fit_esmcmc_set_nthreads:
 * @esmcmc: a #NcmFitESMCMC
 * @nthreads: numbers of simultaneous walkers updates
 *
 * If @nthreads is larger than nwalkers / 2, it will be set to
 * nwalkers / 2.
 *
 */
void 
ncm_fit_esmcmc_set_nthreads (NcmFitESMCMC *esmcmc, guint nthreads)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  if (nthreads > 0)
  {
    if (self->nwalkers % 2 == 1)
      g_error ("ncm_fit_esmcmc_set_nthreads: cannot parallelize with an odd number of walkers [%u].", self->nwalkers);
    if (nthreads > self->nwalkers / 2)
      nthreads = self->nwalkers / 2;
  }
  self->nthreads = nthreads;
}

/**
 * ncm_fit_esmcmc_use_mpi:
 * @esmcmc: a #NcmFitESMCMC
 * @use_mpi: whether to prefer MPI
 *
 * If @use_mpi is TRUE then the paralelization will be acomplished
 * using MPI if any slaves are available. If no slaves are available
 * then it falls back to threads. 
 * 
 * Note that paralelization will only occur if the number of threads
 * set using ncm_fit_esmcmc_set_nthreads() is larger than one.
 * 
 */
void 
ncm_fit_esmcmc_use_mpi (NcmFitESMCMC *esmcmc, gboolean use_mpi)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
	const guint nslaves = ncm_cfg_mpi_nslaves ();
	self->use_mpi = use_mpi;
	self->nslaves = nslaves;
	self->has_mpi = use_mpi && (nslaves > 0);
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  if (self->started)
    g_error ("ncm_fit_esmcmc_set_rng: Cannot change the RNG object during a run, call ncm_fit_esmcmc_end_run() first.");

  ncm_mset_catalog_set_rng (self->mcat, rng);
}

/**
 * ncm_fit_esmcmc_set_auto_trim:
 * @esmcmc: a #NcmFitESMCMC
 * @enable: a boolean
 *
 * If @enable is TRUE turns on the auto-trimming when performing a
 * run_lre.
 *
 */
void 
ncm_fit_esmcmc_set_auto_trim (NcmFitESMCMC *esmcmc, gboolean enable)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  self->auto_trim = enable;
}

/**
 * ncm_fit_esmcmc_set_auto_trim_div:
 * @esmcmc: a #NcmFitESMCMC
 * @div: a unsigned integer
 * 
 * Sets the divisor for the auto trim tests.
 *
 */
void 
ncm_fit_esmcmc_set_auto_trim_div (NcmFitESMCMC *esmcmc, guint div)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  self->auto_trim_div = div;
}

/**
 * ncm_fit_esmcmc_set_auto_trim_type:
 * @esmcmc: a #NcmFitESMCMC
 * @ttype: a #NcmMSetCatalogTrimType
 * 
 * Sets the trim type.
 *
 */
void 
ncm_fit_esmcmc_set_auto_trim_type (NcmFitESMCMC *esmcmc, NcmMSetCatalogTrimType ttype)
{
  NcmFitESMCMCPrivate * const self = esmcmc->priv;
  
  g_assert (ttype & (NCM_MSET_CATALOG_TRIM_TYPE_ALL));

  self->trim_type = ttype;
}

/**
 * ncm_fit_esmcmc_set_min_runs:
 * @esmcmc: a #NcmFitESMCMC
 * @min_runs: a unsigned integer
 * 
 * Sets the minimum number of runs between tests.
 *
 */
void 
ncm_fit_esmcmc_set_min_runs (NcmFitESMCMC *esmcmc, guint min_runs)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  g_assert_cmpuint (min_runs, >, 0);
  self->min_runs = min_runs;
}

/**
 * ncm_fit_esmcmc_set_max_runs_time:
 * @esmcmc: a #NcmFitESMCMC
 * @max_runs_time: a unsigned integer
 * 
 * Sets the maximum time for the runs between tests.
 *
 */
void 
ncm_fit_esmcmc_set_max_runs_time (NcmFitESMCMC *esmcmc, gdouble max_runs_time)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  g_assert_cmpfloat (max_runs_time, >=, 1.0);
  self->max_runs_time = max_runs_time;
}

/**
 * ncm_fit_esmcmc_has_rng:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 *
 * Returns: whether there is a #NcmRNG set.
 */
gboolean
ncm_fit_esmcmc_has_rng (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  return (ncm_mset_catalog_peek_rng (self->mcat) != NULL);
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const gdouble accept_ratio = self->naccepted * 1.0 / (self->ntotal * 1.0);
  
  return accept_ratio;
}

/**
 * ncm_fit_esmcmc_get_offboard_ratio:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_fit_esmcmc_get_offboard_ratio (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const gdouble offboard_ratio = self->noffboard * 1.0 / (self->ntotal * 1.0);
  
  return offboard_ratio;
}

/**
 * ncm_fit_esmcmc_get_accept_ratio_last_update:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_fit_esmcmc_get_accept_ratio_last_update (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const gdouble accept_ratio = self->naccepted_lup * 1.0 / (self->ntotal_lup * 1.0);
  
  return accept_ratio;
}

/**
 * ncm_fit_esmcmc_get_offboard_ratio_last_update:
 * @esmcmc: a #NcmFitESMCMC
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_fit_esmcmc_get_offboard_ratio_last_update (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const gdouble offboard_ratio = self->noffboard_lup * 1.0 / (self->ntotal_lup * 1.0);
  
  return offboard_ratio;
}

void
_ncm_fit_esmcmc_update (NcmFitESMCMC *esmcmc, guint ki, guint kf)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const guint step = self->nwalkers * ((self->n / self->interm_log) == 0 ? 1 : (self->n / self->interm_log));
  guint k;

  g_assert_cmpuint (ki, <, kf);
  g_assert_cmpuint (kf, <=, self->nwalkers);

  self->ntotal_lup    = 0;
  self->naccepted_lup = 0;
  self->noffboard_lup = 0;

  for (k = ki; k < kf; k++)
  {
    NcmVector *full_theta_k = g_ptr_array_index (self->full_theta, k);

    ncm_mset_catalog_add_from_vector (self->mcat, full_theta_k);

    self->cur_sample_id++;
    ncm_timer_task_increment (self->nt);

    self->ntotal++;
    self->ntotal_lup++;
    if (g_array_index (self->accepted, gboolean, k))
    {
      self->naccepted++;
      self->naccepted_lup++;
      g_array_index (self->accepted, gboolean, k) = FALSE;
    }
    if (g_array_index (self->offboard, gboolean, k))
    {
      self->noffboard++;
      self->noffboard_lup++;
      g_array_index (self->offboard, gboolean, k) = FALSE;
    }
  }

  switch (self->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = (self->cur_sample_id + 1) % step;
      gboolean log_timeout = FALSE;

      if ((self->nt->pos_time - self->nt->last_log_time) > self->log_time_interval)
      {
        log_timeout = ((self->cur_sample_id + 1) % self->nwalkers == 0);
      }
    
      if (log_timeout || (stepi == 0) || (self->nt->task_pos == self->nt->task_len))
      {
        NcmVector *e_var = ncm_mset_catalog_peek_current_e_var (self->mcat);
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (self->mcat);
        ncm_mset_catalog_log_current_chain_stats (self->mcat);
        g_message ("# NcmFitESMCMC:acceptance ratio %7.4f%% (last update %7.4f%%), offboard ratio %7.4f%% (last update %7.4f%%).\n", 
                   ncm_fit_esmcmc_get_accept_ratio (esmcmc) * 100.0, ncm_fit_esmcmc_get_accept_ratio_last_update (esmcmc) * 100.0,
                   ncm_fit_esmcmc_get_offboard_ratio (esmcmc) * 100.0, ncm_fit_esmcmc_get_offboard_ratio_last_update (esmcmc) * 100.0);
        g_message ("# NcmFitESMCMC:last ensemble variance of -2ln(L): % 22.15g (2n = %d), min(-2ln(L)) = % 22.15g.\n", 
                    e_var != NULL ? ncm_vector_get (e_var, NCM_FIT_ESMCMC_M2LNL_ID) : GSL_POSINF,
                   2 * self->fparam_len,
                   ncm_mset_catalog_get_bestfit_m2lnL (self->mcat));
        /* ncm_timer_task_accumulate (self->nt, acc); */
        ncm_timer_task_log_elapsed (self->nt);
        ncm_timer_task_log_mean_time (self->nt);
        ncm_timer_task_log_time_left (self->nt);
        ncm_timer_task_log_cur_datetime (self->nt);
        ncm_timer_task_log_end_datetime (self->nt);
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    {
      if ((self->cur_sample_id + 1) % self->nwalkers == 0)
      {
        NcmVector *e_mean = ncm_mset_catalog_peek_current_e_mean (self->mcat);
        NcmVector *e_var  = ncm_mset_catalog_peek_current_e_var (self->mcat);
        self->fit->mtype = self->mtype;

        if (e_mean != NULL)
        {
          ncm_mset_fparams_set_vector_offset (self->fit->mset, e_mean, self->nadd_vals);
          ncm_fit_state_set_m2lnL_curval (self->fit->fstate, ncm_vector_get (e_mean, NCM_FIT_ESMCMC_M2LNL_ID));
        }

        ncm_fit_log_state (self->fit);
        ncm_mset_catalog_log_current_stats (self->mcat);
        ncm_mset_catalog_log_current_chain_stats (self->mcat);
        g_message ("# NcmFitESMCMC:acceptance ratio %7.4f%% (last update %7.4f%%), offboard ratio %7.4f%% (last update %7.4f%%).\n", 
                   ncm_fit_esmcmc_get_accept_ratio (esmcmc) * 100.0, ncm_fit_esmcmc_get_accept_ratio_last_update (esmcmc) * 100.0,
                   ncm_fit_esmcmc_get_offboard_ratio (esmcmc) * 100.0, ncm_fit_esmcmc_get_offboard_ratio_last_update (esmcmc) * 100.0);
        g_message ("# NcmFitESMCMC:last ensemble variance of -2ln(L): % 22.15g (2n = %d), min(-2ln(L)) = % 22.15g.\n", 
                    e_var != NULL ? ncm_vector_get (e_var, NCM_FIT_ESMCMC_M2LNL_ID) : GSL_POSINF,
                   2 * self->fparam_len,
                   ncm_mset_catalog_get_bestfit_m2lnL (self->mcat));
        /* ncm_timer_task_increment (self->nt); */
        ncm_timer_task_log_elapsed (self->nt);
        ncm_timer_task_log_mean_time (self->nt);
        ncm_timer_task_log_time_left (self->nt);
        ncm_timer_task_log_cur_datetime (self->nt);
        ncm_timer_task_log_end_datetime (self->nt);
      }
      break;
    }
  }
}

static void ncm_fit_esmcmc_intern_skip (NcmFitESMCMC *esmcmc, guint n);

static void 
_ncm_fit_esmcmc_gen_init_points_mpi (NcmFitESMCMC *esmcmc, const glong i, const glong f)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmRNG *rng   = ncm_mset_catalog_peek_rng (self->mcat);
	GPtrArray *thetastar_in_a  = g_ptr_array_new ();
	GPtrArray *thetastar_out_a = g_ptr_array_new ();
  glong k, j;

  for (k = i; k < f; k++)
  {
    NcmVector *thetastar_in_k  = g_ptr_array_index (self->thetastar_in, k);
    NcmVector *thetastar_out_k = g_ptr_array_index (self->thetastar_out, k);
    NcmVector *thetastar       = g_ptr_array_index (self->thetastar, k);
		
		ncm_mset_trans_kern_prior_sample (self->sampler, thetastar, rng);
		for (j = 0; j < NCM_FIT_ESMCMC_MPI_IN_LEN; j++)
			ncm_vector_set (thetastar_in_k, self->fparam_len + j, -1.0);

		/*ncm_vector_log_vals (g_ptr_array_index (self->full_thetastar_inout, k), "#  FULL IN: ", "% 22.15g", TRUE);*/
		
		g_ptr_array_add (thetastar_in_a,  thetastar_in_k);
		g_ptr_array_add (thetastar_out_a, thetastar_out_k);
	}

	ncm_mpi_job_run_array (self->mj, thetastar_in_a, thetastar_out_a);

	k = i;
  for (j = 0; j < thetastar_out_a->len; j++)
  {
		NcmVector *thetastar_out_k = g_ptr_array_index (thetastar_out_a, j);

		/*ncm_vector_log_vals (g_ptr_array_index (self->full_thetastar_inout, k), "# FULL OUT: ", "% 22.15g", TRUE);*/
		
		if (ncm_vector_get (thetastar_out_k, 0) != 0.0)
		{
			NcmVector *full_theta_k     = g_ptr_array_index (self->full_theta, k);
			NcmVector *full_thetastar_k = g_ptr_array_index (self->full_thetastar, k);

			/*printf ("# AHA %ld % 22.15g % 22.15g\n", k, ncm_vector_get (full_thetastar_k, 0), ncm_vector_get (full_thetastar_k, 1));*/
			
			ncm_vector_memcpy (full_theta_k, full_thetastar_k);

			g_array_index (self->accepted, gboolean, k) = TRUE;
			k++;
		}
	}

	g_ptr_array_unref (thetastar_in_a);
	g_ptr_array_unref (thetastar_out_a);

	if (k < f)
		_ncm_fit_esmcmc_gen_init_points_mpi (esmcmc, k, f);
}

static void 
_ncm_fit_esmcmc_gen_init_points_mt_eval (glong i, glong f, gpointer data)
{
	NcmFitESMCMC *esmcmc        = NCM_FIT_ESMCMC (data);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmFitESMCMCWorker **fk_ptr = ncm_memory_pool_get (self->walker_pool);
  NcmFit *fit_k               = fk_ptr[0]->fit;
  NcmRNG *rng                 = ncm_mset_catalog_peek_rng (self->mcat);
  glong k;

  for (k = i; k < f; k++)
  {
    NcmVector *full_theta_k = g_ptr_array_index (self->full_theta, k);
    NcmVector *theta_k      = g_ptr_array_index (self->theta, k);
    gdouble *m2lnL          = ncm_vector_ptr (full_theta_k, NCM_FIT_ESMCMC_M2LNL_ID);

    do
    {
      g_mutex_lock (&self->resample_lock);
      ncm_mset_trans_kern_prior_sample (self->sampler, theta_k, rng);
      g_mutex_unlock (&self->resample_lock);

      ncm_mset_fparams_set_vector (fit_k->mset, theta_k);
      ncm_fit_m2lnL_val (fit_k, m2lnL);

    } while (!gsl_finite (m2lnL[0]));

    if (fk_ptr[0]->funcs_array != NULL)
    {
      guint j;
      for (j = 0; j < fk_ptr[0]->funcs_array->len; j++)
      {
        NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (fk_ptr[0]->funcs_array, j));
        const gdouble a_j = ncm_mset_func_eval0 (func, fit_k->mset);

        ncm_vector_set (full_theta_k, j + 1, a_j);
      }
    }

    g_array_index (self->accepted, gboolean, k) = TRUE;
  }

  ncm_memory_pool_return (fk_ptr);
}

static void 
_ncm_fit_esmcmc_gen_init_points (NcmFitESMCMC *esmcmc, const gboolean use_mpi)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  if (self->cur_sample_id + 1 == self->nwalkers)
    return;
  else if (self->cur_sample_id + 1 > self->nwalkers)
    g_error ("_ncm_fit_esmcmc_gen_init_points: initial points already generated.");

	if (self->nthreads > 1)
  {
    ncm_mset_catalog_set_sync_mode (self->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
		if (use_mpi)
		{
			_ncm_fit_esmcmc_gen_init_points_mpi (esmcmc, self->cur_sample_id + 1, self->nwalkers);
		}
		else
			ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_gen_init_points_mt_eval, self->cur_sample_id + 1, self->nwalkers, esmcmc);
  }
  else
  {
    ncm_mset_catalog_set_sync_mode (self->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
    _ncm_fit_esmcmc_gen_init_points_mt_eval (self->cur_sample_id + 1, self->nwalkers, esmcmc);
  }

	_ncm_fit_esmcmc_update (esmcmc, self->cur_sample_id + 1, self->nwalkers); 
  ncm_mset_catalog_sync (self->mcat, FALSE);
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const gint mcat_cur_id = ncm_mset_catalog_get_cur_id (self->mcat);
  gboolean init_point_task = FALSE;
	gboolean read_from_cat   = FALSE;
	const gboolean use_mpi   = (self->has_mpi && (self->nthreads > 0)) ? TRUE : FALSE;
  
  if (self->started)
    g_error ("ncm_fit_esmcmc_start_run: run already started, run ncm_fit_esmcmc_end_run() first.");

  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Starting Ensamble Sampler Markov Chain Monte Carlo.\n");
      g_message ("#   Number of walkers: %.4d.\n", self->nwalkers);
      g_message ("#   Number of threads: %.4d.\n", self->nthreads);
			if (self->nthreads > 0)
				g_message ("#   Using MPI:         %s.\n", self->use_mpi ? ((self->nslaves > 0) ? "yes" : "no - use MPI enabled but no slaves available") : "no");
      ncm_dataset_log_info (self->fit->lh->dset);
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Model set:\n");
      ncm_mset_pretty_log (self->fit->mset);
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Starting Ensamble Sampler Markov Chain Monte Carlo.\n");
      g_message ("#   Number of walkers: %.4d.\n", self->nwalkers);
      g_message ("#   Number of threads: %.4d.\n", self->nthreads);
			if (self->nthreads > 0)
				g_message ("#   Using MPI:         %s.\n", self->use_mpi ? ((self->nslaves > 0) ? "yes" : "no - use MPI enabled but no slaves available") : "no");
      break;
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_mset_catalog_peek_rng (self->mcat) == NULL)
  {
    NcmRNG *rng = ncm_rng_new (NULL);
    
    ncm_rng_set_random_seed (rng, FALSE);
    ncm_fit_esmcmc_set_rng (esmcmc, rng);

    if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
      g_message ("# NcmFitESMCMC: No RNG was defined, using algorithm: `%s' and seed: %lu.\n",
                 ncm_rng_get_algo (rng), ncm_rng_get_seed (rng));

    ncm_rng_free (rng);
  }

  self->started = TRUE;

  ncm_mset_catalog_set_sync_mode (self->mcat, NCM_MSET_CATALOG_SYNC_TIMED);
  ncm_mset_catalog_set_sync_interval (self->mcat, NCM_FIT_ESMCMC_MIN_SYNC_INTERVAL);
  
  ncm_mset_catalog_sync (self->mcat, TRUE);

	if (use_mpi)
	{
		ncm_mpi_job_init_all_slaves (self->mj, self->ser);
		ncm_serialize_reset (self->ser, TRUE);
	}
	
  g_mutex_lock (&self->update_lock);
  self->ntotal        = 0;
  self->naccepted     = 0;
  self->noffboard     = 0;
  self->ntotal_lup    = 0;
  self->naccepted_lup = 0;
  self->noffboard_lup = 0;
  g_mutex_unlock (&self->update_lock);

  if (mcat_cur_id > self->cur_sample_id)
  {
    ncm_fit_esmcmc_intern_skip (esmcmc, mcat_cur_id - self->cur_sample_id);
    g_assert_cmpint (self->cur_sample_id, ==, mcat_cur_id);
  }
  else if (mcat_cur_id < self->cur_sample_id)
    g_error ("ncm_fit_esmcmc_start_run: Unknown error cur_id < cur_sample_id [%d < %d].", 
             mcat_cur_id, self->cur_sample_id);

  if (self->cur_sample_id + 1 < self->nwalkers)
  {
    ncm_timer_task_start (self->nt, self->nwalkers - self->cur_sample_id - 1);
    ncm_timer_set_name (self->nt, "NcmFitESMCMC");
    init_point_task = TRUE;
  }

  if (self->cur_sample_id + 1 < self->nwalkers)
  {
    gint k;

    if (ncm_mset_catalog_get_first_id (self->mcat) > 0)
      g_error ("ncm_fit_esmcmc_start_run: cannot use catalogs with first_id > 0 without a complete first ensemble.");
    
    for (k = 0; k <= self->cur_sample_id; k++)
    {
      NcmVector *cur_row      = ncm_mset_catalog_peek_row (self->mcat, k);
      NcmVector *full_theta_k = g_ptr_array_index (self->full_theta, k);
      g_assert (cur_row != NULL);

			read_from_cat = TRUE;
      ncm_vector_memcpy (full_theta_k, cur_row);
    }

		_ncm_fit_esmcmc_gen_init_points (esmcmc, use_mpi);

    g_mutex_lock (&self->update_lock);
    self->ntotal        = 0;
    self->naccepted     = 0;
    self->noffboard     = 0;
    self->ntotal_lup    = 0;
    self->naccepted_lup = 0;
    self->noffboard_lup = 0;
    g_mutex_unlock (&self->update_lock);
  }
  else
  {
    const guint t   = (self->cur_sample_id + 1) / self->nwalkers;
    const guint ki  = (self->cur_sample_id + 1) % self->nwalkers;
    const guint tm1 = t - 1;
    gint k;

    for (k = 0; k < ki; k++)
    {
      const guint row_n       = ncm_mset_catalog_get_row_from_time (self->mcat, self->nwalkers * t + k);
      NcmVector *cur_row      = ncm_mset_catalog_peek_row (self->mcat, row_n);
      NcmVector *full_theta_k = g_ptr_array_index (self->full_theta, k);

      g_assert (cur_row != NULL);

			read_from_cat = TRUE;
      ncm_vector_memcpy (full_theta_k, cur_row);
    }

    for (k = ki; k < self->nwalkers; k++)
    {
      const guint row_n       = ncm_mset_catalog_get_row_from_time (self->mcat, self->nwalkers * tm1 + k);
      NcmVector *cur_row      = ncm_mset_catalog_peek_row (self->mcat, row_n);
      NcmVector *full_theta_k = g_ptr_array_index (self->full_theta, k);

      g_assert (cur_row != NULL);

			read_from_cat = TRUE;
      ncm_vector_memcpy (full_theta_k, cur_row);
    }
  }

	if (read_from_cat)
	{
		const guint len = ncm_mset_catalog_len (self->mcat);
		if (!ncm_fit_esmcmc_validate (esmcmc, len - self->nwalkers, len))
		{
			if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
			{
				ncm_cfg_msg_sepa ();
				g_message ("# NcmFitESMCMC: Last ensemble failed in the m2lnL check, the catalog may be corrupted, removing last ensemble and retrying...\n");
			}
			ncm_fit_esmcmc_end_run (esmcmc);
			
			ncm_mset_catalog_remove_last_ensemble (self->mcat);

			self->cur_sample_id -= self->nwalkers;
			self->started        = FALSE;

			ncm_fit_esmcmc_start_run (esmcmc);
		}
	}

  if (init_point_task)
    ncm_timer_task_pause (self->nt);
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;

  if (!self->started)
    g_error ("ncm_fit_esmcmc_end_run: run not started, run ncm_fit_esmcmc_start_run() first.");

	ncm_mpi_job_free_all_slaves (self->mj);
	
  if (ncm_timer_task_is_running (self->nt))
    ncm_timer_task_end (self->nt);

  ncm_mset_catalog_sync (self->mcat, TRUE);
  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_mset_catalog_log_current_stats (self->mcat);
  
  self->started = FALSE;
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  self->n               = 0;
  self->cur_sample_id   = -1;
  self->ntotal          = 0;
  self->naccepted       = 0;
  self->noffboard       = 0;
  self->ntotal_lup      = 0;
  self->naccepted_lup   = 0;
  self->noffboard_lup   = 0;
  self->started         = FALSE;  
  ncm_mset_catalog_reset (self->mcat);
}

static void 
ncm_fit_esmcmc_intern_skip (NcmFitESMCMC *esmcmc, guint n)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  if (n == 0)
    return;

  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Skipping %u points (%f iterations), will start at %u-th point.\n", n, n * 1.0 / self->nwalkers, self->cur_sample_id + n + 1 + 1);
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  self->cur_sample_id += n;
}

static void _ncm_fit_esmcmc_run (NcmFitESMCMC *esmcmc);

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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  guint ti = (self->cur_sample_id + 1) / self->nwalkers;
  guint ki = (self->cur_sample_id + 1) % self->nwalkers;
  
  if (!self->started)
    g_error ("ncm_fit_esmcmc_run: run not started, run ncm_fit_esmcmc_start_run() first.");

  if (n <= ti)
  {
    if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Nothing to do, current Monte Carlo run is %d\n", ti);
    }
    return;
  }
  
  self->n = n - ti;
  
  switch (self->mtype)
  {
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      ncm_cfg_msg_sepa ();
      g_message ("# NcmFitESMCMC: Calculating [%06d] Ensemble Sampler Markov Chain Monte Carlo runs [%s]\n", 
                 self->n, ncm_fit_esmcmc_walker_desc (self->walker));
    }
    case NCM_FIT_RUN_MSGS_NONE:
      break;
  }

  if (ncm_timer_task_is_running (self->nt))
  {
    ncm_timer_task_add_tasks (self->nt, self->n * self->nwalkers - ki);
    ncm_timer_task_continue (self->nt);
  }
  else
  {
    ncm_timer_task_start (self->nt, self->n * self->nwalkers - ki);
    ncm_timer_set_name (self->nt, "NcmFitESMCMC");
  }
  if (self->mtype > NCM_FIT_RUN_MSGS_NONE)
    ncm_timer_task_log_start_datetime (self->nt);

  _ncm_fit_esmcmc_run (esmcmc);

  ncm_timer_task_pause (self->nt);
}

static void 
_ncm_fit_esmcmc_eval_mpi (NcmFitESMCMC *esmcmc, const glong i, const glong f)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
	GPtrArray *thetastar_in_a  = g_ptr_array_new ();
	GPtrArray *thetastar_out_a = g_ptr_array_new ();
  glong k;

  for (k = i; k < f; k++)
  {
		NcmVector *full_theta_k    = g_ptr_array_index (self->full_theta, k);
    NcmVector *thetastar_in_k  = g_ptr_array_index (self->thetastar_in, k);
    NcmVector *thetastar_out_k = g_ptr_array_index (self->thetastar_out, k);
    NcmVector *thetastar       = g_ptr_array_index (self->thetastar, k);

		ncm_fit_esmcmc_walker_step (self->walker, self->theta, self->m2lnL, thetastar, k);

		ncm_vector_set (thetastar_in_k, self->fparam_len + 0, ncm_vector_get (full_theta_k, NCM_FIT_ESMCMC_M2LNL_ID));
		ncm_vector_set (thetastar_in_k, self->fparam_len + 1, ncm_fit_esmcmc_walker_prob_norm (self->walker, self->theta, self->m2lnL, thetastar, k));
		ncm_vector_set (thetastar_in_k, self->fparam_len + 2, ncm_vector_get (self->jumps, k));

		/*ncm_vector_log_vals (g_ptr_array_index (self->full_thetastar_inout, k), "#  FULL IN: ", "% 22.15g", TRUE);*/

		/*if (ncm_mset_fparam_valid_bounds (self->fit->mset, thetastar))*/
    if (ncm_mset_fparam_validate_all (self->fit->mset, thetastar)) /* TEST! */
		{
			g_ptr_array_add (thetastar_in_a,  thetastar_in_k);
			g_ptr_array_add (thetastar_out_a, thetastar_out_k);
		}
		else
		{
			ncm_vector_set (thetastar_out_k, 0, 0.0);
		}
	}

	ncm_mpi_job_run_array (self->mj, thetastar_in_a, thetastar_out_a);

  for (k = i; k < f; k++)
  {
		NcmVector *thetastar_out_k = g_ptr_array_index (self->thetastar_out, k);

		/*ncm_vector_log_vals (g_ptr_array_index (self->full_thetastar_inout, k), "# FULL OUT: ", "% 22.15g", TRUE);*/
		
		if (ncm_vector_get (thetastar_out_k, 0) != 0.0)
		{
			NcmVector *full_theta_k     = g_ptr_array_index (self->full_theta, k);
			NcmVector *full_thetastar_k = g_ptr_array_index (self->full_thetastar, k);
			
			ncm_vector_memcpy (full_theta_k, full_thetastar_k);

			g_array_index (self->accepted, gboolean, k) = TRUE;
		}
	}

	g_ptr_array_unref (thetastar_in_a);
	g_ptr_array_unref (thetastar_out_a);
}

static void 
_ncm_fit_esmcmc_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc        = NCM_FIT_ESMCMC (data);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmFitESMCMCWorker **fk_ptr = ncm_memory_pool_get (self->walker_pool);
  NcmFit *fit_k               = fk_ptr[0]->fit;
  guint k = i;

  while (k < f)
  {
    NcmVector *full_thetastar = g_ptr_array_index (self->full_thetastar, k);
    NcmVector *full_theta_k   = g_ptr_array_index (self->full_theta, k);
    NcmVector *thetastar      = g_ptr_array_index (self->thetastar, k);

    gdouble *m2lnL_cur  = ncm_vector_ptr (full_theta_k, NCM_FIT_ESMCMC_M2LNL_ID);
    gdouble *m2lnL_star = ncm_vector_ptr (full_thetastar, NCM_FIT_ESMCMC_M2LNL_ID);
    gdouble jump        = ncm_vector_get (self->jumps, k);
    gdouble prob = 0.0;

    ncm_fit_esmcmc_walker_step (self->walker, self->theta, self->m2lnL, thetastar, k);

		/*if (ncm_mset_fparam_valid_bounds (fit_k->mset, thetastar))*/
    if (ncm_mset_fparam_validate_all (fit_k->mset, thetastar)) /* TEST! */
    {
      ncm_mset_fparams_set_vector (fit_k->mset, thetastar);
      ncm_fit_m2lnL_val (fit_k, m2lnL_star);

      if (gsl_finite (m2lnL_star[0]))
      {
        prob = ncm_fit_esmcmc_walker_prob (self->walker, self->theta, self->m2lnL, thetastar, k, m2lnL_cur[0], m2lnL_star[0]);
        prob = GSL_MIN (prob, 1.0);
      }
    }
    else
    {
      g_array_index (self->offboard, gboolean, k) = TRUE;
    }
    
    if (jump < prob)
    {
      if (fk_ptr[0]->funcs_array != NULL)
      {
        guint j;
        for (j = 0; j < fk_ptr[0]->funcs_array->len; j++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (fk_ptr[0]->funcs_array, j));
          const gdouble a_j = ncm_mset_func_eval0 (func, fit_k->mset);

          ncm_vector_set (full_thetastar, j + 1, a_j);
        }
      }

      ncm_vector_memcpy (full_theta_k, full_thetastar);
      g_array_index (self->accepted, gboolean, k) = TRUE;
    }
    k++;
  }

  ncm_memory_pool_return (fk_ptr);
}

static void
_ncm_fit_esmcmc_get_jumps (NcmFitESMCMC *esmcmc, guint ki, guint kf)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmRNG *rng = ncm_mset_catalog_peek_rng (self->mcat);
  guint k;
  
  for (k = ki; k < kf; k++)
  {
    const gdouble jump = gsl_rng_uniform (rng->r);
    ncm_vector_set (self->jumps, k, jump);
  }
}

static void 
_ncm_fit_esmcmc_run_serial (NcmFitESMCMC *esmcmc, const glong i, const glong f)
{
	_ncm_fit_esmcmc_mt_eval (i, f, esmcmc);
}

static void 
_ncm_fit_esmcmc_run_mt (NcmFitESMCMC *esmcmc, const glong i, const glong f)
{
	ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, i, f, esmcmc);
}

static void
_ncm_fit_esmcmc_run (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmRNG *rng      = ncm_mset_catalog_peek_rng (self->mcat);
  gboolean mthread = (self->nthreads > 1);
	void (*run) (NcmFitESMCMC *, const glong, const glong) = mthread ? (self->has_mpi ? _ncm_fit_esmcmc_eval_mpi : _ncm_fit_esmcmc_run_mt) : _ncm_fit_esmcmc_run_serial;
	const guint nwalkers_2 = self->nwalkers / 2;
	guint ki               = (self->cur_sample_id + 1) % self->nwalkers;
	guint i;

	ncm_mset_catalog_set_sync_mode (self->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
	if (self->n > 0)
	{
		_ncm_fit_esmcmc_get_jumps (esmcmc, ki, self->nwalkers);
		ncm_fit_esmcmc_walker_setup (self->walker, self->theta, self->m2lnL, ki, self->nwalkers, rng);

		if (ki < nwalkers_2)
		{
			run (esmcmc, ki, nwalkers_2);
			run (esmcmc, nwalkers_2, self->nwalkers);
		}
		else
		{
			run (esmcmc, ki, self->nwalkers);
		}

		ncm_fit_esmcmc_walker_clean (self->walker, ki, self->nwalkers);

		_ncm_fit_esmcmc_update (esmcmc, ki, self->nwalkers);
		ncm_mset_catalog_timed_sync (self->mcat, FALSE);

		for (i = 1; i < self->n; i++)
		{
			_ncm_fit_esmcmc_get_jumps (esmcmc, 0, self->nwalkers);
			ncm_fit_esmcmc_walker_setup (self->walker, self->theta, self->m2lnL, 0, self->nwalkers, rng);

			run (esmcmc, 0, nwalkers_2);
			run (esmcmc, nwalkers_2, self->nwalkers);

			ncm_fit_esmcmc_walker_clean (self->walker, 0, self->nwalkers);

			_ncm_fit_esmcmc_update (esmcmc, 0, self->nwalkers);
			ncm_mset_catalog_timed_sync (self->mcat, FALSE);
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  gdouble lerror, post_lnnorm_sd;
  const gdouble lre2 = lre * lre;
  const guint catlen = ncm_mset_catalog_len (self->mcat) / self->nwalkers;

  g_assert_cmpfloat (lre, >, 0.0);

  prerun = GSL_MAX (prerun, 100);

  if (catlen < prerun)
  {
    guint prerun_left = prerun - catlen;
    
    if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitESMCMC: Running first %u pre-runs...\n", prerun_left);

    ncm_fit_esmcmc_run (esmcmc, prerun);
    if (self->auto_trim)
    {
      ncm_mset_catalog_trim_by_type (self->mcat, self->auto_trim_div, self->trim_type, self->mtype);
    }
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (self->mcat, FALSE);
  lerror = ncm_mset_catalog_largest_error (self->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n             = ncm_mset_catalog_len (self->mcat);
    gdouble m             = n * lerror2 / lre2;
    guint runs            = ((m - n) > 1000.0) ? MIN (ceil ((m - n) * self->lre_step), 100000) : ceil (m - n);
    guint ti              = (self->cur_sample_id + 1) / self->nwalkers;

    runs = GSL_MIN (ncm_timer_task_estimate_by_time (self->nt, self->max_runs_time), runs);
    runs = GSL_MAX (runs / self->nwalkers + 1, self->min_runs);
    
    if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      gdouble glnvol;
      const gdouble lnevol = ncm_mset_catalog_get_post_lnvol (self->mcat, ncm_c_stats_1sigma (), &glnvol);
      g_message ("# NcmFitESMCMC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmFitESMCMC: ln (eVol) = % 22.15g; ln (gVol) = % 22.15g; lnNorm = % 22.15g\n", lnevol, glnvol, ncm_mset_catalog_get_post_lnnorm (self->mcat, &post_lnnorm_sd));
      g_message ("# NcmFitESMCMC: Running more %u runs...\n", runs);
    }

    ncm_fit_esmcmc_run (esmcmc, ti + runs);
    if (self->auto_trim)
    {
      ncm_mset_catalog_trim_by_type (self->mcat, self->auto_trim_div, self->trim_type, self->mtype);
    }

    ncm_mset_catalog_estimate_autocorrelation_tau (self->mcat, FALSE);
    lerror = ncm_mset_catalog_largest_error (self->mcat);
  }

  if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# NcmFitESMCMC: Largest relative error %e attained: %e\n", lre, lerror);
  }
}

/**
 * ncm_fit_esmcmc_run_burnin:
 * @esmcmc: a #NcmFitESMCMC
 * @prerun: FIXME
 * @ntimes: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_fit_esmcmc_run_burnin (NcmFitESMCMC *esmcmc, guint prerun, guint ntimes)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  const guint catlen = ncm_mset_catalog_len (self->mcat) / self->nwalkers;
  NcmVector *var;
  gdouble m2lnL_var;
  guint cb;
  guint ti;

  prerun = GSL_MAX (prerun, 20);

  if (catlen < prerun)
  {
    guint prerun_left = prerun - catlen;
    
    if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitESMCMC: Running first %u pre-runs...\n", prerun_left);

    ncm_fit_esmcmc_run (esmcmc, prerun);
    if (self->auto_trim)
    {
      ncm_mset_catalog_trim_by_type (self->mcat, self->auto_trim_div, self->trim_type, self->mtype);
    }
  }
  ncm_mset_catalog_estimate_autocorrelation_tau (self->mcat, FALSE);
  
  var       = ncm_mset_catalog_peek_current_e_var (self->mcat);
  m2lnL_var = ncm_vector_get (var, NCM_FIT_ESMCMC_M2LNL_ID);
  cb        = ncm_mset_catalog_calc_const_break (self->mcat, NCM_FIT_ESMCMC_M2LNL_ID, self->mtype);
  ti        = (self->cur_sample_id + 1) / self->nwalkers;

  while (ntimes * cb > ti)
  {
    const guint truns = ntimes * cb;
    const guint runs  = (truns > ti) ? truns : (ti + 10);

    if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmFitESMCMC: `%02d'-times the estimated constant break not attained: %6u %6u, var(m2lnL) = %22.15g.\n", ntimes, ti, runs, m2lnL_var);
      g_message ("# NcmFitESMCMC: Running more %u runs...\n", runs - ti);
    }

    ncm_fit_esmcmc_run (esmcmc, runs);
    if (self->auto_trim)
    {
      ncm_mset_catalog_trim_by_type (self->mcat, self->auto_trim_div, self->trim_type, self->mtype);
    }
    ncm_mset_catalog_estimate_autocorrelation_tau (self->mcat, FALSE);

    var       = ncm_mset_catalog_peek_current_e_var (self->mcat);
    m2lnL_var = ncm_vector_get (var, NCM_FIT_ESMCMC_M2LNL_ID);
    cb        = ncm_mset_catalog_calc_const_break (self->mcat, NCM_FIT_ESMCMC_M2LNL_ID, self->mtype);
    ti        = (self->cur_sample_id + 1) / self->nwalkers;
  }
  if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
  {
    g_message ("# NcmFitESMCMC: `%02d'-times the estimated constant break achieved, var(m2lnL) = %22.15g.\n", ntimes, m2lnL_var);
  }
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmMSet *mset = ncm_mset_catalog_peek_mset (self->mcat);
  
  ncm_mset_catalog_get_mean (self->mcat, &self->fit->fstate->fparams);
  ncm_mset_catalog_get_covar (self->mcat, &self->fit->fstate->covar);
  ncm_mset_fparams_set_vector (mset, self->fit->fstate->fparams);
  
  self->fit->fstate->has_covar = TRUE;
}

/**
 * ncm_fit_esmcmc_peek_ser:
 * @esmcmc: a #NcmFitESMCMC
 *
 * Peeks the internal #NcmSerialize object from @esmcmc.
 * 
 * Returns: (transfer none): the internal #NcmSerialize object.
 */
NcmSerialize *
ncm_fit_esmcmc_peek_ser (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  return self->ser;
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
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  return ncm_mset_catalog_ref (self->mcat);
}

/**
 * ncm_fit_esmcmc_peek_catalog:
 * @esmcmc: a #NcmFitESMCMC
 *
 * Gets the generated catalog of @esmcmc.
 * 
 * Returns: (transfer none): the generated catalog.
 */
NcmMSetCatalog *
ncm_fit_esmcmc_peek_catalog (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  return self->mcat;
}

/**
 * ncm_fit_esmcmc_peek_fit:
 * @esmcmc: a #NcmFitESMCMC
 *
 * Gets the #NcmFit object used by the sampler.
 * 
 * Returns: (transfer none): the #NcmFit object.
 */
NcmFit *
ncm_fit_esmcmc_peek_fit (NcmFitESMCMC *esmcmc)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  return self->fit;
}

static void 
_ncm_fit_esmcmc_validade_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc        = NCM_FIT_ESMCMC (data);
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  NcmFitESMCMCWorker **fk_ptr = ncm_memory_pool_get (self->walker_pool);
  NcmFit *fit_k               = fk_ptr[0]->fit;
  glong k;

  for (k = i; k < f; k++)
  {
    NcmVector *cur_row      = ncm_mset_catalog_peek_row (self->mcat, k);
    const gdouble row_m2lnL = ncm_vector_get (cur_row, NCM_FIT_ESMCMC_M2LNL_ID);
    gdouble m2lnL           = 0.0;
    gdouble diff;
    
    g_assert (cur_row != NULL);

    ncm_mset_fparams_set_vector_offset (fit_k->mset, cur_row, self->nadd_vals);

    ncm_fit_m2lnL_val (fit_k, &m2lnL);

    diff = fabs ((row_m2lnL - m2lnL) / row_m2lnL);
    
    if (diff > 1.0e-3)
    {
      if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      {
				g_mutex_lock (&self->resample_lock);
        ncm_message ("# NcmFitESMCMC: Catalog row %5lu: m2lnL = %20.15g, recalculated to % 20.15g, diff = %8.5e <====== FAILED.\n",
                     k, row_m2lnL, m2lnL, diff);
				ncm_message ("# NcmFitESMCMC: row %5lu values: ", k);
				ncm_vector_log_vals (cur_row, "", "%22.15g", TRUE);
				g_mutex_unlock (&self->resample_lock);
      }

			g_mutex_lock (&self->update_lock);
			g_array_index (self->accepted, gboolean, 0) = FALSE;
			g_mutex_unlock (&self->update_lock);
		}
    else if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
			g_mutex_lock (&self->resample_lock);
      ncm_message ("# NcmFitESMCMC: Catalog row %5lu: m2lnL = %20.15g, recalculated to % 20.15g, diff = %8.5e, SUCCEEDED!\n",
                   k, row_m2lnL, m2lnL, diff);
			if (self->mtype >= NCM_FIT_RUN_MSGS_FULL)
			{
				ncm_message ("# NcmFitESMCMC: row %5lu values:", k);
				ncm_vector_log_vals (cur_row, "", "%22.15g", TRUE);
			}
			g_mutex_unlock (&self->resample_lock);
    }
  }

  ncm_memory_pool_return (fk_ptr);
}

static void 
_ncm_fit_esmcmc_validate_mpi (NcmFitESMCMC *esmcmc, const glong i, const glong f)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
	GPtrArray *thetastar_in_a  = g_ptr_array_new ();
	GPtrArray *thetastar_out_a = g_ptr_array_new ();
	const glong len            = MIN (f - i, self->thetastar->len);
  glong j, k;

  for (j = 0; j < len; j++)
  {
		glong k                     = i + j;
		NcmVector *cur_row          = ncm_mset_catalog_peek_row (self->mcat, k);
    NcmVector *thetastar_in_j   = g_ptr_array_index (self->thetastar_in, j);
    NcmVector *thetastar_out_j  = g_ptr_array_index (self->thetastar_out, j);
    NcmVector *full_thetastar_j = g_ptr_array_index (self->full_thetastar, j);
		
		ncm_vector_memcpy (full_thetastar_j, cur_row);
		
		ncm_vector_set (thetastar_in_j, self->fparam_len + 0, ncm_vector_get (cur_row, NCM_FIT_ESMCMC_M2LNL_ID));
		ncm_vector_set (thetastar_in_j, self->fparam_len + 1, -1.0);
		ncm_vector_set (thetastar_in_j, self->fparam_len + 2, -1.0);
		
		g_ptr_array_add (thetastar_in_a,  thetastar_in_j);
		g_ptr_array_add (thetastar_out_a, thetastar_out_j);
	}

	ncm_mpi_job_run_array (self->mj, thetastar_in_a, thetastar_out_a);

	k = i;
  for (j = 0; j < thetastar_out_a->len; j++)
  {
		NcmVector *thetastar_out_j = g_ptr_array_index (thetastar_out_a, j);
		NcmVector *thetastar_in_j  = g_ptr_array_index (thetastar_in_a, j);
		const gdouble m2lnL_cur    = ncm_vector_get (thetastar_in_j, self->fparam_len + 0);
		const gdouble m2lnL        = ncm_vector_get (thetastar_out_j, NCM_FIT_ESMCMC_MPI_OUT_LEN + NCM_FIT_ESMCMC_M2LNL_ID);
    const gdouble diff         = fabs ((m2lnL_cur - m2lnL) / m2lnL_cur);
    
    if (diff > 1.0e-3)
    {
      if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      {
				NcmVector *cur_row = ncm_mset_catalog_peek_row (self->mcat, k);
        ncm_message ("# NcmFitESMCMC: Catalog row %5lu: m2lnL = %20.15g, recalculated to % 20.15g, diff = %8.5e <====== FAILED.\n",
                     k, m2lnL_cur, m2lnL, diff);
				ncm_message ("# NcmFitESMCMC: row %5lu values: ", k);
				ncm_vector_log_vals (cur_row, "", "%22.15g", TRUE);
      }
			g_array_index (self->accepted, gboolean, 0) = FALSE;
		}
    else if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      ncm_message ("# NcmFitESMCMC: Catalog row %5lu: m2lnL = %20.15g, recalculated to % 20.15g, diff = %8.5e, SUCCEEDED!\n",
                   k, m2lnL_cur, m2lnL, diff);
			if (self->mtype >= NCM_FIT_RUN_MSGS_FULL)
			{
				NcmVector *cur_row = ncm_mset_catalog_peek_row (self->mcat, k);
				ncm_message ("# NcmFitESMCMC: row %5lu values:", k);
				ncm_vector_log_vals (cur_row, "", "%22.15g", TRUE);
			}
    }
		k++;
	}

	g_ptr_array_unref (thetastar_in_a);
	g_ptr_array_unref (thetastar_out_a);

	if (k < f)
		_ncm_fit_esmcmc_validate_mpi (esmcmc, k, f);
}


/**
 * ncm_fit_esmcmc_validate:
 * @esmcmc: a #NcmFitESMCMC
 * @pi: initial position
 * @pf: final position
 *
 * Recalculates the value of $-2\ln(L)$ and compares
 * with the values found in the catalog. This function
 * is particularly useful to check if any problem occured
 * during a multithread evaluation of the likelihood.
 *
 * Choosing @pf == 0 performs the validation from  @pi to the end.
 * 
 * Returns: Whether the validation was TRUE or FALSE.
 */
gboolean
ncm_fit_esmcmc_validate (NcmFitESMCMC *esmcmc, gulong pi, gulong pf)
{
	NcmFitESMCMCPrivate * const self = esmcmc->priv;
  gulong len;
  
  ncm_mset_catalog_sync (self->mcat, TRUE);

  len = ncm_mset_catalog_len (self->mcat);

	if (pf == 0)
		pf = len;

  g_assert_cmpuint (pi, <, pf);
  g_assert_cmpuint (pf, <=, len);

  if (self->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# NcmFitESMCMC: validating catalog rows [%lu, %lu)\n", pi, pf);
		if (self->mtype >= NCM_FIT_RUN_MSGS_FULL)
			ncm_fit_log_info (self->fit);
	}

	g_array_index (self->accepted, gboolean, 0) = TRUE;
  if (self->nthreads > 1)
	{
		if (self->has_mpi)
			_ncm_fit_esmcmc_validate_mpi (esmcmc, pi, pf);
		else
			ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_validade_mt_eval, pi, pf, esmcmc);
	}
  else
    _ncm_fit_esmcmc_validade_mt_eval (pi, pf, esmcmc);

  return g_array_index (self->accepted, gboolean, 0);
}

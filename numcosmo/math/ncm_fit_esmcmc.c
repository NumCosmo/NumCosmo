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
  PROP_WALKER,
  PROP_AUTO_TRIM,
  PROP_AUTO_TRIM_DIV,
  PROP_TRIM_TYPE,
  PROP_MIN_RUNS,
  PROP_MAX_RUNS_TIME,
  PROP_MTYPE,
  PROP_NTHREADS,
  PROP_DATA_FILE,
  PROP_FUNCS_ARRAY,
};

G_DEFINE_TYPE (NcmFitESMCMC, ncm_fit_esmcmc, G_TYPE_OBJECT);

static gpointer _ncm_fit_esmcmc_worker_dup (gpointer userdata);
static void _ncm_fit_esmcmc_worker_free (gpointer p);

static void
ncm_fit_esmcmc_init (NcmFitESMCMC *esmcmc)
{
  esmcmc->fit             = NULL;
  esmcmc->walker_pool     = ncm_memory_pool_new (&_ncm_fit_esmcmc_worker_dup, esmcmc, 
                                                 &_ncm_fit_esmcmc_worker_free);

  esmcmc->sampler         = NULL;
  esmcmc->mcat            = NULL;
  esmcmc->mtype           = NCM_FIT_RUN_MSGS_NONE;
  esmcmc->nt              = ncm_timer_new ();
  esmcmc->ser             = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  esmcmc->walker          = NULL;
  esmcmc->auto_trim       = FALSE;
  esmcmc->auto_trim_div   = 0;
  esmcmc->trim_type       = 0;
  esmcmc->min_runs        = 0;
  esmcmc->max_runs_time   = 0.0;
  esmcmc->nadd_vals       = 0;
  esmcmc->fparam_len      = 0;

  esmcmc->theta           = g_ptr_array_new ();
  esmcmc->thetastar       = g_ptr_array_new ();
  esmcmc->full_theta      = g_ptr_array_new ();
  esmcmc->full_thetastar  = g_ptr_array_new ();

  g_ptr_array_set_free_func (esmcmc->theta, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (esmcmc->thetastar, (GDestroyNotify) &ncm_vector_free);

  g_ptr_array_set_free_func (esmcmc->full_theta, (GDestroyNotify) &ncm_vector_free);
  g_ptr_array_set_free_func (esmcmc->full_thetastar, (GDestroyNotify) &ncm_vector_free);

  esmcmc->jumps           = NULL;
  esmcmc->accepted        = g_array_new (TRUE, TRUE, sizeof (gboolean));
  esmcmc->offboard        = g_array_new (TRUE, TRUE, sizeof (gboolean));

  esmcmc->funcs_oa        = NULL;
  esmcmc->funcs_oa_file   = NULL;
  
  esmcmc->nthreads        = 0;
  esmcmc->n               = 0;
  esmcmc->nwalkers        = 0;
  esmcmc->cur_sample_id   = -1; /* Represents that no samples were calculated yet, i.e., id of the last added point. */
  esmcmc->ntotal          = 0;
  esmcmc->naccepted       = 0;
  esmcmc->noffboard       = 0;
  esmcmc->started         = FALSE;

  g_mutex_init (&esmcmc->dup_fit);
  g_mutex_init (&esmcmc->resample_lock);
  g_mutex_init (&esmcmc->update_lock);
  g_cond_init (&esmcmc->write_cond);
}

static void
_ncm_fit_esmcmc_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_esmcmc_parent_class)->constructed (object);
  {
    NcmFitESMCMC *esmcmc = NCM_FIT_ESMCMC (object);
    guint k;

    g_assert_cmpint (esmcmc->nwalkers, >, 0);
    esmcmc->fparam_len = ncm_mset_fparam_len (esmcmc->fit->mset);
    {
      const guint nfuncs    = (esmcmc->funcs_oa != NULL) ? esmcmc->funcs_oa->len : 0;
      const guint nadd_vals = esmcmc->nadd_vals = nfuncs + 1;
      const guint theta_len = esmcmc->fparam_len + nadd_vals;

      if (nfuncs > 0)
      {
        gchar **names   = g_new (gchar *, nadd_vals + 1);
        gchar **symbols = g_new (gchar *, nadd_vals + 1);

        names[0]   = g_strdup (NCM_MSET_CATALOG_M2LNL_COLNAME);
        symbols[0] = g_strdup (NCM_MSET_CATALOG_M2LNL_SYMBOL);

        for (k = 0; k < nfuncs; k++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (esmcmc->funcs_oa, k));
          g_assert (NCM_IS_MSET_FUNC (func));

          names[1 + k]   = g_strdup (ncm_mset_func_peek_uname (func));
          symbols[1 + k] = g_strdup (ncm_mset_func_peek_usymbol (func));
        }
        names[1 + k]   = NULL;
        symbols[1 + k] = NULL;

        esmcmc->mcat = ncm_mset_catalog_new_array (esmcmc->fit->mset, nadd_vals, esmcmc->nwalkers, FALSE, 
                                                   names, symbols);

        g_strfreev (names);
        g_strfreev (symbols);
      }
      else
      {
        esmcmc->mcat = ncm_mset_catalog_new (esmcmc->fit->mset, nadd_vals, esmcmc->nwalkers, FALSE, 
                                             NCM_MSET_CATALOG_M2LNL_COLNAME, NCM_MSET_CATALOG_M2LNL_SYMBOL, 
                                             NULL);
      }
      ncm_mset_catalog_set_run_type (esmcmc->mcat, "Ensemble Sampler MCMC");

      for (k = 0; k < esmcmc->nwalkers; k++)
      {
        NcmVector *full_theta_k     = ncm_vector_new (theta_len);
        NcmVector *full_thetastar_k = ncm_vector_new (theta_len);
        NcmVector *theta_k          = ncm_vector_get_subvector (full_theta_k, nadd_vals, esmcmc->fparam_len);
        NcmVector *thetastar_k      = ncm_vector_get_subvector (full_thetastar_k, nadd_vals, esmcmc->fparam_len);

        g_ptr_array_add (esmcmc->theta, theta_k);
        g_ptr_array_add (esmcmc->thetastar, thetastar_k);

        g_ptr_array_add (esmcmc->full_theta, full_theta_k);
        g_ptr_array_add (esmcmc->full_thetastar, full_thetastar_k);
      }
    }

    esmcmc->jumps = ncm_vector_new (esmcmc->nwalkers);
    g_array_set_size (esmcmc->accepted, esmcmc->nwalkers);
    g_array_set_size (esmcmc->offboard, esmcmc->nwalkers);
    
    if (esmcmc->walker == NULL)
      esmcmc->walker = ncm_fit_esmcmc_walker_new_from_name ("NcmFitESMCMCWalkerStretch");

    ncm_fit_esmcmc_walker_set_size (esmcmc->walker, esmcmc->nwalkers);
    ncm_fit_esmcmc_walker_set_nparams (esmcmc->walker, esmcmc->fparam_len);
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
    case PROP_WALKER:
      esmcmc->walker = g_value_dup_object (value);
      break;
    case PROP_AUTO_TRIM:
      esmcmc->auto_trim = g_value_get_boolean (value);
      break;
    case PROP_AUTO_TRIM_DIV:
      esmcmc->auto_trim_div = g_value_get_uint (value);
      break;
    case PROP_TRIM_TYPE:
      esmcmc->trim_type = g_value_get_flags (value);
      break;
    case PROP_MIN_RUNS:
      esmcmc->min_runs = g_value_get_uint (value);
      break;
    case PROP_MAX_RUNS_TIME:
      esmcmc->max_runs_time = g_value_get_double (value);
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
    case PROP_FUNCS_ARRAY:
    {
      esmcmc->funcs_oa = g_value_dup_boxed (value);
      if (esmcmc->funcs_oa != NULL)
      {
        guint i;
        for (i = 0; i < esmcmc->funcs_oa->len; i++)
        {
          NcmMSetFunc *func = NCM_MSET_FUNC (ncm_obj_array_peek (esmcmc->funcs_oa, i));
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
    case PROP_WALKER:
      g_value_set_object (value, esmcmc->walker);
      break;
    case PROP_AUTO_TRIM:
      g_value_set_boolean (value, esmcmc->auto_trim);
      break;
    case PROP_AUTO_TRIM_DIV:
      g_value_set_uint (value, esmcmc->auto_trim_div);
      break;
    case PROP_TRIM_TYPE:
      g_value_set_flags (value, esmcmc->trim_type);
      break;
    case PROP_MIN_RUNS:
      g_value_set_uint (value, esmcmc->min_runs);
      break;
    case PROP_MAX_RUNS_TIME:
      g_value_set_double (value, esmcmc->max_runs_time);
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
    case PROP_FUNCS_ARRAY:
      g_value_set_boxed (value, esmcmc->funcs_oa);
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

  ncm_fit_esmcmc_walker_clear (&esmcmc->walker);

  ncm_vector_clear (&esmcmc->jumps);

  ncm_obj_array_clear (&esmcmc->funcs_oa);

  g_clear_pointer (&esmcmc->funcs_oa_file, g_free);
  
  g_clear_pointer (&esmcmc->theta, g_ptr_array_unref);
  g_clear_pointer (&esmcmc->thetastar, g_ptr_array_unref);

  g_clear_pointer (&esmcmc->full_theta, g_ptr_array_unref);
  g_clear_pointer (&esmcmc->full_thetastar, g_ptr_array_unref);

  g_clear_pointer (&esmcmc->accepted, g_array_unref);
  g_clear_pointer (&esmcmc->offboard, g_array_unref);

  if (esmcmc->walker_pool != NULL)
  {
    ncm_memory_pool_free (esmcmc->walker_pool, TRUE);
    esmcmc->walker_pool = NULL;
  }
  
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
                                   PROP_WALKER,
                                   g_param_spec_object ("walker",
                                                        NULL,
                                                        "Walker object",
                                                        NCM_TYPE_FIT_ESMCMC_WALKER,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
                                   PROP_DATA_FILE,
                                   g_param_spec_string ("data-file",
                                                        NULL,
                                                        "Data filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FUNCS_ARRAY,
                                   g_param_spec_boxed ("functions-array",
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

  G_LOCK (dup_thread);
  {
    NcmFitESMCMCWorker *fw = g_new (NcmFitESMCMCWorker, 1);
    
    fw->fit = ncm_fit_dup (esmcmc->fit, esmcmc->ser);

    if (esmcmc->funcs_oa != NULL)
      fw->funcs_array = ncm_obj_array_dup (esmcmc->funcs_oa, esmcmc->ser);    
    else
      fw->funcs_array = NULL;

    ncm_serialize_reset (esmcmc->ser, TRUE);
    
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
  g_assert (esmcmc->fit == NULL);
  esmcmc->fit = ncm_fit_ref (fit);
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
                                       "nwalkers",        nwalkers,
                                       "sampler",         sampler,
                                       "walker",          walker,
                                       "mtype",           mtype,
                                       "functions-array", funcs_array,
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

  if (esmcmc->funcs_oa != NULL)
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    gchar *base_name      = ncm_util_basename_fits (filename);
    esmcmc->funcs_oa_file = g_strdup_printf ("%s.oa", base_name);
    g_free (base_name); 

    ncm_obj_array_save (esmcmc->funcs_oa, ser, esmcmc->funcs_oa_file, TRUE);

    ncm_serialize_free (ser);
  }
  
  if (esmcmc->started)
    g_assert_cmpint (esmcmc->cur_sample_id, ==, ncm_mset_catalog_get_cur_id (esmcmc->mcat));
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
  esmcmc->auto_trim = enable;
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
  esmcmc->auto_trim_div = div;
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
  g_assert_cmpuint (min_runs, >, 0);
  esmcmc->min_runs = min_runs;
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
  g_assert_cmpfloat (max_runs_time, >=, 1.0);
  esmcmc->max_runs_time = max_runs_time;
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
  return (ncm_mset_catalog_peek_rng (esmcmc->mcat) != NULL);
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
  gdouble accept_ratio;
  accept_ratio = esmcmc->naccepted * 1.0 / (esmcmc->ntotal * 1.0);
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
  gdouble offboard_ratio;
  offboard_ratio = esmcmc->noffboard * 1.0 / (esmcmc->ntotal * 1.0);
  return offboard_ratio;
}

void
_ncm_fit_esmcmc_update (NcmFitESMCMC *esmcmc, guint ki, guint kf)
{
  const guint part = 5;
  const guint step = esmcmc->nwalkers * ((esmcmc->n / part) == 0 ? 1 : (esmcmc->n / part));
  guint k;

  g_assert_cmpuint (ki, <, kf);
  g_assert_cmpuint (kf, <=, esmcmc->nwalkers);

  for (k = ki; k < kf; k++)
  {
    NcmVector *full_theta_k = g_ptr_array_index (esmcmc->full_theta, k);

    ncm_mset_catalog_add_from_vector (esmcmc->mcat, full_theta_k);

    esmcmc->cur_sample_id++;
    ncm_timer_task_increment (esmcmc->nt);

    esmcmc->ntotal++;
    if (g_array_index (esmcmc->accepted, gboolean, k))
    {
      esmcmc->naccepted++;
      g_array_index (esmcmc->accepted, gboolean, k) = FALSE;
    }
    if (g_array_index (esmcmc->offboard, gboolean, k))
    {
      esmcmc->noffboard++;
      g_array_index (esmcmc->offboard, gboolean, k) = FALSE;
    }
    
  }

  switch (esmcmc->mtype)
  {
    case NCM_FIT_RUN_MSGS_NONE:
      break;
    case NCM_FIT_RUN_MSGS_SIMPLE:
    {
      guint stepi = (esmcmc->cur_sample_id + 1) % step;
      gboolean log_timeout = FALSE;

      if ((esmcmc->nt->pos_time - esmcmc->nt->last_log_time) > 60.0)
      {
        log_timeout = ((esmcmc->cur_sample_id + 1) % esmcmc->nwalkers == 0);
      }
    
      if (log_timeout || (stepi == 0) || (esmcmc->nt->task_pos == esmcmc->nt->task_len))
      {
        /* guint acc = stepi == 0 ? step : stepi; */
        ncm_mset_catalog_log_current_stats (esmcmc->mcat);
        ncm_mset_catalog_log_current_chain_stats (esmcmc->mcat);
        g_message ("# NcmFitESMCMC:acceptance ratio %7.4f%%, offboard ratio %7.4f%%.\n", 
                   ncm_fit_esmcmc_get_accept_ratio (esmcmc) * 100.0,
                   ncm_fit_esmcmc_get_offboard_ratio (esmcmc) * 100.0);
        /* ncm_timer_task_accumulate (esmcmc->nt, acc); */
        ncm_timer_task_log_elapsed (esmcmc->nt);
        ncm_timer_task_log_mean_time (esmcmc->nt);
        ncm_timer_task_log_time_left (esmcmc->nt);
        ncm_timer_task_log_cur_datetime (esmcmc->nt);
        ncm_timer_task_log_end_datetime (esmcmc->nt);
      }
      break;
    }
    default:
    case NCM_FIT_RUN_MSGS_FULL:
    {
      if ((esmcmc->cur_sample_id + 1) % esmcmc->nwalkers == 0)
      {
        NcmVector *e_mean = ncm_mset_catalog_peek_current_e_mean (esmcmc->mcat);
        esmcmc->fit->mtype = esmcmc->mtype;

        if (e_mean != NULL)
        {
          ncm_mset_fparams_set_vector_offset (esmcmc->fit->mset, e_mean, esmcmc->nadd_vals);
          ncm_fit_state_set_m2lnL_curval (esmcmc->fit->fstate, ncm_vector_get (e_mean, NCM_FIT_ESMCMC_M2LNL_ID));
        }

        ncm_fit_log_state (esmcmc->fit);
        ncm_mset_catalog_log_current_stats (esmcmc->mcat);
        ncm_mset_catalog_log_current_chain_stats (esmcmc->mcat);
        g_message ("# NcmFitESMCMC:acceptance ratio %7.4f%%, offboard ratio %7.4f%%.\n", 
                   ncm_fit_esmcmc_get_accept_ratio (esmcmc) * 100.0,
                   ncm_fit_esmcmc_get_offboard_ratio (esmcmc) * 100.0);
        /* ncm_timer_task_increment (esmcmc->nt); */
        ncm_timer_task_log_elapsed (esmcmc->nt);
        ncm_timer_task_log_mean_time (esmcmc->nt);
        ncm_timer_task_log_time_left (esmcmc->nt);
        ncm_timer_task_log_cur_datetime (esmcmc->nt);
        ncm_timer_task_log_end_datetime (esmcmc->nt);
      }
      break;
    }
  }
}

static void ncm_fit_esmcmc_intern_skip (NcmFitESMCMC *esmcmc, guint n);

static void 
_ncm_fit_esmcmc_gen_init_points_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc        = NCM_FIT_ESMCMC (data);
  NcmFitESMCMCWorker **fk_ptr = ncm_memory_pool_get (esmcmc->walker_pool);
  NcmFit *fit_k               = fk_ptr[0]->fit;
  NcmRNG *rng                 = ncm_mset_catalog_peek_rng (esmcmc->mcat);
  glong k;

  for (k = i; k < f; k++)
  {
    NcmVector *full_theta_k = g_ptr_array_index (esmcmc->full_theta, k);
    NcmVector *theta_k      = g_ptr_array_index (esmcmc->theta, k);
    gdouble *m2lnL          = ncm_vector_ptr (full_theta_k, NCM_FIT_ESMCMC_M2LNL_ID);

    do
    {
      g_mutex_lock (&esmcmc->resample_lock);
      ncm_mset_trans_kern_prior_sample (esmcmc->sampler, theta_k, rng);
      g_mutex_unlock (&esmcmc->resample_lock);

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

    g_array_index (esmcmc->accepted, gboolean, k) = TRUE;
  }

  ncm_memory_pool_return (fk_ptr);
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
    ncm_mset_catalog_set_sync_mode (esmcmc->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
    ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_gen_init_points_mt_eval, esmcmc->cur_sample_id + 1, esmcmc->nwalkers, esmcmc);
  }
  else
  {
    ncm_mset_catalog_set_sync_mode (esmcmc->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
    _ncm_fit_esmcmc_gen_init_points_mt_eval (esmcmc->cur_sample_id + 1, esmcmc->nwalkers, esmcmc);
  }

  _ncm_fit_esmcmc_update (esmcmc, esmcmc->cur_sample_id + 1, esmcmc->nwalkers); 
  ncm_mset_catalog_sync (esmcmc->mcat, FALSE);
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
  const gint mcat_cur_id = ncm_mset_catalog_get_cur_id (esmcmc->mcat);
  gboolean init_point_task = FALSE;
	gboolean read_from_cat   = FALSE;
  
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

  if (ncm_mset_catalog_peek_rng (esmcmc->mcat) == NULL)
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

  ncm_mset_catalog_set_sync_mode (esmcmc->mcat, NCM_MSET_CATALOG_SYNC_TIMED);
  ncm_mset_catalog_set_sync_interval (esmcmc->mcat, NCM_FIT_ESMCMC_MIN_SYNC_INTERVAL);
  
  ncm_mset_catalog_sync (esmcmc->mcat, TRUE);

  g_mutex_lock (&esmcmc->update_lock);
  esmcmc->ntotal    = 0;
  esmcmc->naccepted = 0;
  esmcmc->noffboard = 0;
  g_mutex_unlock (&esmcmc->update_lock);

  if (mcat_cur_id > esmcmc->cur_sample_id)
  {
    ncm_fit_esmcmc_intern_skip (esmcmc, mcat_cur_id - esmcmc->cur_sample_id);
    g_assert_cmpint (esmcmc->cur_sample_id, ==, mcat_cur_id);
  }
  else if (mcat_cur_id < esmcmc->cur_sample_id)
    g_error ("ncm_fit_esmcmc_start_run: Unknown error cur_id < cur_sample_id [%d < %d].", 
             mcat_cur_id, esmcmc->cur_sample_id);

  if (esmcmc->cur_sample_id + 1 < esmcmc->nwalkers)
  {
    ncm_timer_task_start (esmcmc->nt, esmcmc->nwalkers - esmcmc->cur_sample_id - 1);
    ncm_timer_set_name (esmcmc->nt, "NcmFitESMCMC");
    init_point_task = TRUE;
  }

  if (esmcmc->cur_sample_id + 1 < esmcmc->nwalkers)
  {
    gint k;

    if (ncm_mset_catalog_get_first_id (esmcmc->mcat) > 0)
      g_error ("ncm_fit_esmcmc_start_run: cannot use catalogs with first_id > 0 without a complete first ensemble.");
    
    for (k = 0; k <= esmcmc->cur_sample_id; k++)
    {
      NcmVector *cur_row      = ncm_mset_catalog_peek_row (esmcmc->mcat, k);
      NcmVector *full_theta_k = g_ptr_array_index (esmcmc->full_theta, k);
      g_assert (cur_row != NULL);

			read_from_cat = TRUE;
      ncm_vector_memcpy (full_theta_k, cur_row);
    }
    _ncm_fit_esmcmc_gen_init_points (esmcmc);

    g_mutex_lock (&esmcmc->update_lock);
    esmcmc->ntotal    = 0;
    esmcmc->naccepted = 0;
    esmcmc->noffboard = 0;
    g_mutex_unlock (&esmcmc->update_lock);
  }
  else
  {
    const guint t   = (esmcmc->cur_sample_id + 1) / esmcmc->nwalkers;
    const guint ki  = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;
    const guint tm1 = t - 1;
    gint k;

    for (k = 0; k < ki; k++)
    {
      const guint row_n       = ncm_mset_catalog_get_row_from_time (esmcmc->mcat, esmcmc->nwalkers * t + k);
      NcmVector *cur_row      = ncm_mset_catalog_peek_row (esmcmc->mcat, row_n);
      NcmVector *full_theta_k = g_ptr_array_index (esmcmc->full_theta, k);

      g_assert (cur_row != NULL);

			read_from_cat = TRUE;
      ncm_vector_memcpy (full_theta_k, cur_row);
    }

    for (k = ki; k < esmcmc->nwalkers; k++)
    {
      const guint row_n       = ncm_mset_catalog_get_row_from_time (esmcmc->mcat, esmcmc->nwalkers * tm1 + k);
      NcmVector *cur_row      = ncm_mset_catalog_peek_row (esmcmc->mcat, row_n);
      NcmVector *full_theta_k = g_ptr_array_index (esmcmc->full_theta, k);

      g_assert (cur_row != NULL);

			read_from_cat = TRUE;
      ncm_vector_memcpy (full_theta_k, cur_row);
    }
  }

	if (read_from_cat)
	{
		const guint len = ncm_mset_catalog_len (esmcmc->mcat);
		if (!ncm_fit_esmcmc_validate (esmcmc, len - esmcmc->nwalkers, len))
		{
			if (esmcmc->mtype > NCM_FIT_RUN_MSGS_NONE)
			{
				ncm_cfg_msg_sepa ();
				g_message ("# NcmFitESMCMC: Last ensemble failed in the m2lnL check, the catalog may be corrupted, removing last ensemble and retrying...\n");
			}

			ncm_mset_catalog_remove_last_ensemble (esmcmc->mcat);

			esmcmc->cur_sample_id -= esmcmc->nwalkers;
			esmcmc->started        = FALSE;

			ncm_fit_esmcmc_start_run (esmcmc);
		}
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
  esmcmc->noffboard       = 0;
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
      g_message ("# NcmFitESMCMC: Calculating [%06d] Ensemble Sampler Markov Chain Monte Carlo runs [%s]\n", 
                 esmcmc->n, ncm_fit_esmcmc_walker_desc (esmcmc->walker));
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

  _ncm_fit_esmcmc_run (esmcmc);

  ncm_timer_task_pause (esmcmc->nt);
}

static void 
_ncm_fit_esmcmc_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc        = NCM_FIT_ESMCMC (data);
  NcmFitESMCMCWorker **fk_ptr = ncm_memory_pool_get (esmcmc->walker_pool);
  NcmFit *fit_k               = fk_ptr[0]->fit;
  guint k = i;

  while (k < f)
  {
    NcmVector *full_thetastar = g_ptr_array_index (esmcmc->full_thetastar, k);
    NcmVector *full_theta_k   = g_ptr_array_index (esmcmc->full_theta, k);
    NcmVector *thetastar      = g_ptr_array_index (esmcmc->thetastar, k);

    gdouble *m2lnL_cur  = ncm_vector_ptr (full_theta_k, NCM_FIT_ESMCMC_M2LNL_ID);
    gdouble *m2lnL_star = ncm_vector_ptr (full_thetastar, NCM_FIT_ESMCMC_M2LNL_ID);
    gdouble jump        = ncm_vector_get (esmcmc->jumps, k);
    gdouble prob = 0.0;

    ncm_fit_esmcmc_walker_step (esmcmc->walker, esmcmc->theta, thetastar, k);

    if (ncm_mset_fparam_valid_bounds (fit_k->mset, thetastar))
    {
      ncm_mset_fparams_set_vector (fit_k->mset, thetastar);
      ncm_fit_m2lnL_val (fit_k, m2lnL_star);

      if (gsl_finite (m2lnL_star[0]))
      {
        prob = ncm_fit_esmcmc_walker_prob (esmcmc->walker, esmcmc->theta, thetastar, k, m2lnL_cur[0], m2lnL_star[0]);
        prob = GSL_MIN (prob, 1.0);
      }
    }
    else
    {
      g_array_index (esmcmc->offboard, gboolean, k) = TRUE;
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
      g_array_index (esmcmc->accepted, gboolean, k) = TRUE;
    }
    k++;
  }

  ncm_memory_pool_return (fk_ptr);
}

static void
_ncm_fit_esmcmc_get_jumps (NcmFitESMCMC *esmcmc, guint ki, guint kf)
{
  NcmRNG *rng = ncm_mset_catalog_peek_rng (esmcmc->mcat);
  guint k;
  
  for (k = ki; k < kf; k++)
  {
    const gdouble jump = gsl_rng_uniform (rng->r);; 
    ncm_vector_set (esmcmc->jumps, k, jump);
  }
}


static void
_ncm_fit_esmcmc_run (NcmFitESMCMC *esmcmc)
{
  NcmRNG *rng      = ncm_mset_catalog_peek_rng (esmcmc->mcat);
  gboolean mthread = (esmcmc->nthreads > 1);
  guint i;

  if (mthread)
  {
    const guint nwalkers_2 = esmcmc->nwalkers / 2;
    guint ki = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;

    ncm_mset_catalog_set_sync_mode (esmcmc->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
    if (esmcmc->n > 0)
    {
      _ncm_fit_esmcmc_get_jumps (esmcmc, ki, esmcmc->nwalkers);
      ncm_fit_esmcmc_walker_setup (esmcmc->walker, esmcmc->theta, ki, esmcmc->nwalkers, rng);
      
      if (ki < nwalkers_2)
      {
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, ki, nwalkers_2, esmcmc);
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, nwalkers_2, esmcmc->nwalkers, esmcmc);
      }
      else
      {
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, ki, esmcmc->nwalkers, esmcmc);
      }

      ncm_fit_esmcmc_walker_clean (esmcmc->walker, ki, esmcmc->nwalkers);

      _ncm_fit_esmcmc_update (esmcmc, ki, esmcmc->nwalkers);
      ncm_mset_catalog_timed_sync (esmcmc->mcat, FALSE);

      for (i = 1; i < esmcmc->n; i++)
      {
        _ncm_fit_esmcmc_get_jumps (esmcmc, 0, esmcmc->nwalkers);
        ncm_fit_esmcmc_walker_setup (esmcmc->walker, esmcmc->theta, 0, esmcmc->nwalkers, rng);
        
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, 0, nwalkers_2, esmcmc);
        ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_mt_eval, nwalkers_2, esmcmc->nwalkers, esmcmc);

        ncm_fit_esmcmc_walker_clean (esmcmc->walker, 0, esmcmc->nwalkers);

        _ncm_fit_esmcmc_update (esmcmc, 0, esmcmc->nwalkers);
        ncm_mset_catalog_timed_sync (esmcmc->mcat, FALSE);
      }
    }
  }
  else
  {
    const guint nwalkers_2 = esmcmc->nwalkers / 2;
    guint ki = (esmcmc->cur_sample_id + 1) % esmcmc->nwalkers;

    ncm_mset_catalog_set_sync_mode (esmcmc->mcat, NCM_MSET_CATALOG_SYNC_DISABLE);
    if (esmcmc->n > 0)
    {
      _ncm_fit_esmcmc_get_jumps (esmcmc, ki, esmcmc->nwalkers);
      ncm_fit_esmcmc_walker_setup (esmcmc->walker, esmcmc->theta, ki, esmcmc->nwalkers, rng);
      
      if (ki < nwalkers_2)
      {
        _ncm_fit_esmcmc_mt_eval (ki, nwalkers_2, esmcmc);
        _ncm_fit_esmcmc_mt_eval (nwalkers_2, esmcmc->nwalkers, esmcmc);
      }
      else
      {
        _ncm_fit_esmcmc_mt_eval (ki, esmcmc->nwalkers, esmcmc);
      }
      
      ncm_fit_esmcmc_walker_clean (esmcmc->walker, ki, esmcmc->nwalkers);

      _ncm_fit_esmcmc_update (esmcmc, ki, esmcmc->nwalkers);
      ncm_mset_catalog_timed_sync (esmcmc->mcat, FALSE);

      for (i = 1; i < esmcmc->n; i++)
      {
        _ncm_fit_esmcmc_get_jumps (esmcmc, 0, esmcmc->nwalkers);
        ncm_fit_esmcmc_walker_setup (esmcmc->walker, esmcmc->theta, 0, esmcmc->nwalkers, rng);

        _ncm_fit_esmcmc_mt_eval (0, nwalkers_2, esmcmc);
        _ncm_fit_esmcmc_mt_eval (nwalkers_2, esmcmc->nwalkers, esmcmc);

        ncm_fit_esmcmc_walker_clean (esmcmc->walker, 0, esmcmc->nwalkers);
        
        _ncm_fit_esmcmc_update (esmcmc, 0, esmcmc->nwalkers);
        ncm_mset_catalog_timed_sync (esmcmc->mcat, FALSE);
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
  const guint catlen = ncm_mset_catalog_len (esmcmc->mcat) / esmcmc->nwalkers;

  g_assert_cmpfloat (lre, >, 0.0);

  prerun = GSL_MAX (prerun, 100);

  if (catlen < prerun)
  {
    guint prerun_left = prerun - catlen;
    if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      g_message ("# NcmFitESMCMC: Running first %u pre-runs...\n", prerun_left);
    ncm_fit_esmcmc_run (esmcmc, prerun);
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (esmcmc->mcat, FALSE);
  lerror = ncm_mset_catalog_largest_error (esmcmc->mcat);

  while (lerror > lre)
  {
    const gdouble lerror2 = lerror * lerror;
    gdouble n             = ncm_mset_catalog_len (esmcmc->mcat);
    gdouble m             = n * lerror2 / lre2;
    guint runs            = ((m - n) > 1000.0) ? ceil ((m - n) * 1.0e-1) : ceil (m - n);
    guint ti              = (esmcmc->cur_sample_id + 1) / esmcmc->nwalkers;

    runs = GSL_MIN (ncm_timer_task_estimate_by_time (esmcmc->nt, esmcmc->max_runs_time), runs);
    runs = GSL_MAX (runs / esmcmc->nwalkers + 1, esmcmc->min_runs);
    
    if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
      g_message ("# NcmFitESMCMC: Largest relative error %e not attained: %e\n", lre, lerror);
      g_message ("# NcmFitESMCMC: Running more %u runs...\n", runs);
    }

    ncm_fit_esmcmc_run (esmcmc, ti + runs);
    if (esmcmc->auto_trim)
    {
      ncm_mset_catalog_trim_by_type (esmcmc->mcat, esmcmc->auto_trim_div, esmcmc->trim_type, esmcmc->mtype);
    }

    ncm_mset_catalog_estimate_autocorrelation_tau (esmcmc->mcat, FALSE);
    lerror = ncm_mset_catalog_largest_error (esmcmc->mcat);
  }

  if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# NcmFitESMCMC: Largest relative error %e attained: %e\n", lre, lerror);
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
  ncm_mset_catalog_get_mean (esmcmc->mcat, &esmcmc->fit->fstate->fparams);
  ncm_mset_catalog_get_covar (esmcmc->mcat, &esmcmc->fit->fstate->covar);
  ncm_mset_fparams_set_vector (esmcmc->mcat->mset, esmcmc->fit->fstate->fparams);
  esmcmc->fit->fstate->has_covar = TRUE;
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
  return esmcmc->ser;
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

static void 
_ncm_fit_esmcmc_validade_mt_eval (glong i, glong f, gpointer data)
{
  NcmFitESMCMC *esmcmc        = NCM_FIT_ESMCMC (data);
  NcmFitESMCMCWorker **fk_ptr = ncm_memory_pool_get (esmcmc->walker_pool);
  NcmFit *fit_k               = fk_ptr[0]->fit;
  glong k;

  for (k = i; k < f; k++)
  {
    NcmVector *cur_row      = ncm_mset_catalog_peek_row (esmcmc->mcat, k);
    const gdouble row_m2lnL = ncm_vector_get (cur_row, NCM_FIT_ESMCMC_M2LNL_ID);
    gdouble m2lnL           = 0.0;
    gdouble diff;
    
    g_assert (cur_row != NULL);

    ncm_mset_fparams_set_vector_offset (fit_k->mset, cur_row, esmcmc->nadd_vals);

    ncm_fit_m2lnL_val (fit_k, &m2lnL);

    diff = fabs ((row_m2lnL - m2lnL) / row_m2lnL);
    
    if (diff > 1.0e-3)
    {
      if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
      {
				g_mutex_lock (&esmcmc->resample_lock);
        ncm_message ("# NcmFitESMCMC: Catalog row %5lu: m2lnL = %20.15g, recalculated to % 20.15g, diff = %8.5e <====== FAILED.\n",
                     k, row_m2lnL, m2lnL, diff);
				ncm_message ("# NcmFitESMCMC: row %5lu values: ", k);
				ncm_vector_log_vals (cur_row, "", "%22.15g", TRUE);
				g_mutex_unlock (&esmcmc->resample_lock);
      }

			g_mutex_lock (&esmcmc->update_lock);
			g_array_index (esmcmc->accepted, gboolean, 0) = FALSE;
			g_mutex_unlock (&esmcmc->update_lock);
		}
    else if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
    {
			g_mutex_lock (&esmcmc->resample_lock);
      ncm_message ("# NcmFitESMCMC: Catalog row %5lu: m2lnL = %20.15g, recalculated to % 20.15g, diff = %8.5e, SUCCEEDED!\n",
                   k, row_m2lnL, m2lnL, diff);
			if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_FULL)
			{
				ncm_message ("# NcmFitESMCMC: row %5lu values:", k);
				ncm_vector_log_vals (cur_row, "", "%22.15g", TRUE);
			}
			g_mutex_unlock (&esmcmc->resample_lock);
    }
  }

  ncm_memory_pool_return (fk_ptr);
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
  gulong len;
  
  ncm_mset_catalog_sync (esmcmc->mcat, TRUE);

  len = ncm_mset_catalog_len (esmcmc->mcat);

	if (pf == 0)
		pf = len;

  g_assert_cmpuint (pi, <, pf);
  g_assert_cmpuint (pf, <=, len);

  if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_SIMPLE)
  {
    ncm_cfg_msg_sepa ();
    g_message ("# NcmFitESMCMC: validating catalog rows [%lu, %lu)\n", pi, pf);
		if (esmcmc->mtype >= NCM_FIT_RUN_MSGS_FULL)
			ncm_fit_log_info (esmcmc->fit);
	}

	g_array_index (esmcmc->accepted, gboolean, 0) = TRUE;
  if (esmcmc->nthreads > 1)
    ncm_func_eval_threaded_loop_full (&_ncm_fit_esmcmc_validade_mt_eval, pi, pf, esmcmc);
  else
    _ncm_fit_esmcmc_validade_mt_eval (pi, pf, esmcmc);

  return g_array_index (esmcmc->accepted, gboolean, 0);
}

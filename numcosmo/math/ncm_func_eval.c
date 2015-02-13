/***************************************************************************
 *            ncm_func_eval.c
 *
 *  Sun Jul 25 19:50:02 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:ncm_func_eval
 * @title: NcmFuncEval
 * @short_description: A general purpose multi-threaded function evaluator.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_func_eval.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include <stdio.h>

typedef struct _NcmFuncEvalCtrl
{
  gint active_threads;
  GMutex *update;
  GCond *finish;
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32))
  GMutex update_m;
  GCond finish_c;
#endif
} NcmFuncEvalCtrl;

typedef struct _NcmFuncEvalLoopEval
{
  NcmFuncEvalLoop lfunc;
  glong i;
  glong f;
  gpointer data;
  NcmFuncEvalCtrl *ctrl;
} NcmFuncEvalLoopEval;

static GThreadPool *_function_thread_pool = NULL;

static void
func (gpointer data, gpointer empty)
{
  NcmFuncEvalLoopEval *arg = (NcmFuncEvalLoopEval *)data;
  NcmFuncEvalCtrl *ctrl = arg->ctrl;
  NCM_UNUSED (empty);
  arg->lfunc (arg->i, arg->f, arg->data);
  g_slice_free (NcmFuncEvalLoopEval, arg);

  g_mutex_lock (ctrl->update);

  ctrl->active_threads--;
  if (ctrl->active_threads == 0)
    g_cond_signal (ctrl->finish);

  g_mutex_unlock (ctrl->update);

  return;
}

/**
 * ncm_func_eval_get_pool: (skip)
 *
 * Allocate if its not yet allocated and return the
 * internal GThreadPool pool.
 *
 * Returns: the pointer to the internal GThreadPool pool
 */
GThreadPool *
ncm_func_eval_get_pool ()
{
  _NCM_STATIC_MUTEX_DECL (create_lock);
  GError *err = NULL;

  _NCM_MUTEX_LOCK (&create_lock);
  if (_function_thread_pool == NULL)
  {
    _function_thread_pool = g_thread_pool_new (func, NULL, NCM_THREAD_POOL_MAX, TRUE, &err);
    if (err != NULL)
      g_error ("ncm_func_eval_get_pool: %s", err->message);
    g_clear_error (&err);
  }
  _NCM_MUTEX_UNLOCK (&create_lock);

  return _function_thread_pool;
}

/**
 * ncm_func_eval_set_max_threads:
 * @mt: new max threads to be used in the pool, -1 means unlimited
 *
 * Set the new maximun number of threads to be used by the pool. Note that this
 * function is global changing this will affect every place which uses these
 * functions.
 *
 */
void
ncm_func_eval_set_max_threads (gint mt)
{
  GError *err = NULL;
  ncm_func_eval_get_pool ();
  g_thread_pool_set_max_threads (_function_thread_pool, mt, &err);
  if (err != NULL)
    g_error ("ncm_func_eval_set_max_threads: %s", err->message);
}

/**
 * ncm_func_eval_threaded_loop_nw:
 * @lfunc: (scope notified): #NcmFuncEvalLoop to be evaluated in threads
 * @i: initial index
 * @f: final index
 * @data: pointer to be passed to @fl
 * @nworkers: number of workers.
 *
 * Using the thread pool, evaluate @fl in each value of (@f-@i)/@nwork.
 *
 */
#if NCM_THREAD_POOL_MAX > 1
void
ncm_func_eval_threaded_loop_nw (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data, guint nworkers)
{
  NcmFuncEvalCtrl ctrl = {0, NULL, NULL, 
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32))
    {NULL},
    {NULL},
#endif
  };
  guint delta, res;

  ncm_func_eval_get_pool ();
#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  ctrl.update = g_mutex_new ();
  ctrl.finish = g_cond_new ();
#else
  g_mutex_init (&ctrl.update_m);
  g_cond_init (&ctrl.finish_c);
  ctrl.update = &ctrl.update_m;
  ctrl.finish = &ctrl.finish_c;
#endif

  g_assert_cmpuint (f, >, i);
  g_assert_cmpuint (f - i, >, nworkers);
  g_assert_cmpuint (nworkers, >, 0);

  delta = (f - i) / nworkers;
  res = (f - i) % nworkers;

  if (delta == 0)
  {
    lfunc (i, f, data);
  }
  else
  {
    GError *err = NULL;
    glong li = i;
    glong lf = i + delta + res;
    ctrl.active_threads = nworkers;

    do {
      NcmFuncEvalLoopEval *arg = g_slice_new (NcmFuncEvalLoopEval);
      arg->lfunc = lfunc;
      arg->i = li;
      arg->f = lf;
      arg->data = data;
      arg->ctrl = &ctrl;
      g_thread_pool_push (_function_thread_pool, arg, &err);
      li = lf;
      lf += delta;
    } while (--nworkers);
  }

  g_mutex_lock (ctrl.update);
  while (ctrl.active_threads != 0)
    g_cond_wait (ctrl.finish, ctrl.update);
  g_mutex_unlock (ctrl.update);

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  g_mutex_free (ctrl.update);
  g_cond_free (ctrl.finish);
#else
  g_mutex_clear (ctrl.update);
  g_cond_clear (ctrl.finish);
#endif
}
#else
void
ncm_func_eval_threaded_loop_nw (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data, guint nworkers)
{
  lfunc (i, f, data);
}
#endif

/**
 * ncm_func_eval_threaded_loop:
 * @lfunc: (scope notified): #NcmFuncEvalLoop to be evaluated in threads
 * @i: initial index
 * @f: final index
 * @data: pointer to be passed to @fl
 *
 * Using the thread pool, evaluate @fl in each value of (@f-@i)/nthreads
 *
 */
void
ncm_func_eval_threaded_loop (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data)
{
  ncm_func_eval_get_pool ();
  {
    guint nthreads = g_thread_pool_get_max_threads (_function_thread_pool);
    ncm_func_eval_threaded_loop_nw (lfunc, i, f, data, nthreads);
  }
}

/**
 * ncm_func_eval_threaded_loop_full:
 * @lfunc: (scope notified): #NcmFuncEvalLoop to be evaluated in threads
 * @i: initial index
 * @f: final index
 * @data: pointer to be passed to @fl
 *
 * Using the thread pool, evaluate @fl sending one worker per index.
 *
 */
#if NCM_THREAD_POOL_MAX > 1
void
ncm_func_eval_threaded_loop_full (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data)
{
  NcmFuncEvalCtrl ctrl = {0, NULL, NULL, 
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32))
    {NULL},
    {NULL},
#endif
  };

  ncm_func_eval_get_pool ();
#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  ctrl.update = g_mutex_new ();
  ctrl.finish = g_cond_new ();
#else
  g_mutex_init (&ctrl.update_m);
  g_cond_init (&ctrl.finish_c);
  ctrl.update = &ctrl.update_m;
  ctrl.finish = &ctrl.finish_c;
#endif

  g_assert_cmpuint (f, >, i);

  {
    GError *err = NULL;
    glong l;
    ctrl.active_threads = f - i;

    for (l = i; l < f; l++)
    {
      NcmFuncEvalLoopEval *arg = g_slice_new (NcmFuncEvalLoopEval);
      arg->lfunc = lfunc;
      arg->i = l;
      arg->f = l + 1;
      arg->data = data;
      arg->ctrl = &ctrl;
      g_thread_pool_push (_function_thread_pool, arg, &err);
    }
  }

  g_mutex_lock (ctrl.update);
  while (ctrl.active_threads != 0)
    g_cond_wait (ctrl.finish, ctrl.update);
  g_mutex_unlock (ctrl.update);

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
  g_mutex_free (ctrl.update);
  g_cond_free (ctrl.finish);
#else
  g_mutex_clear (ctrl.update);
  g_cond_clear (ctrl.finish);
#endif
}
#else
void
ncm_func_eval_threaded_loop_full (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data)
{
  lfunc (i, f, data);
}
#endif

void 
ncm_func_eval_log_pool_stats ()
{
  ncm_func_eval_get_pool ();
  g_message  ("# NcmThreadPool:Unused:      %d\n", g_thread_pool_get_num_unused_threads ());
  g_message  ("# NcmThreadPool:Max Unused:  %d\n", g_thread_pool_get_max_unused_threads ());
  g_message  ("# NcmThreadPool:Running:     %d\n", g_thread_pool_get_num_threads (_function_thread_pool));
  g_message  ("# NcmThreadPool:Unprocessed: %d\n", g_thread_pool_unprocessed (_function_thread_pool));
  g_message  ("# NcmThreadPool:Unused:      %d\n", g_thread_pool_get_max_threads (_function_thread_pool));  
}

/***************************************************************************
 *            ncm_func_eval.c
 *
 *  Sun Jul 25 19:50:02 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcmFuncEval:
 *
 * Parallel loops over an index range on a global thread pool.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_func_eval.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <stdio.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmFuncEvalCtrl
{
  gint active_threads;
  GMutex update;
  GCond finish;
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
  NcmFuncEvalLoopEval *arg = (NcmFuncEvalLoopEval *) data;
  NcmFuncEvalCtrl *ctrl    = arg->ctrl;

  NCM_UNUSED (empty);
  arg->lfunc (arg->i, arg->f, arg->data);
  g_slice_free (NcmFuncEvalLoopEval, arg);

  g_mutex_lock (&ctrl->update);

  ctrl->active_threads--;

  if (ctrl->active_threads == 0)
    g_cond_signal (&ctrl->finish);

  g_mutex_unlock (&ctrl->update);

  return;
}

/**
 * ncm_func_eval_get_pool: (skip)
 *
 * Returns the global thread pool, creating it on first call.
 *
 * Returns: the global #GThreadPool.
 */
GThreadPool *
ncm_func_eval_get_pool (void)
{
  G_LOCK_DEFINE_STATIC (create_lock);

  GError *err = NULL;

  G_LOCK (create_lock);

  if (_function_thread_pool == NULL)
  {
    _function_thread_pool = g_thread_pool_new (func, NULL, NCM_THREAD_POOL_MAX, FALSE, &err);

    if (err != NULL)
      g_error ("ncm_func_eval_get_pool: %s", err->message);

    g_clear_error (&err);
  }

  G_UNLOCK (create_lock);

  return _function_thread_pool;
}

/**
 * ncm_func_eval_set_max_threads:
 * @mt: maximum number of threads
 *
 * Sets the maximum number of threads of the global pool, used by every caller
 * of the ncm_func_eval functions. A value of -1 means unlimited.
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
 * @lfunc: (scope notified): loop function
 * @i: first index
 * @f: one past the last index
 * @data: user data
 * @nworkers: number of workers
 *
 * Splits [@i, @f) into @nworkers contiguous ranges, the first one also taking
 * the remainder, and calls @lfunc with @data on each range in the pool. Returns
 * after all calls finish. Uses at most @f - @i workers, and runs a single serial
 * call if the pool has zero threads. Requires @f > @i and @nworkers > 0.
 */
void
ncm_func_eval_threaded_loop_nw (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data, guint nworkers)
{
  NcmFuncEvalCtrl ctrl = {0};
  guint delta, res;

  ncm_func_eval_get_pool ();

  g_mutex_init (&ctrl.update);
  g_cond_init (&ctrl.finish);

  g_assert_cmpint (f, >, i);
  g_assert_cmpuint (nworkers, >, 0);

  nworkers = MIN (nworkers, (gulong) (f - i));
  delta    = (f - i) / nworkers;
  res      = (f - i) % nworkers;

  if ((g_thread_pool_get_max_threads (_function_thread_pool) == 0) || (delta == 0))
  {
    lfunc (i, f, data);
  }
  else
  {
    GError *err = NULL;
    glong li    = i;
    glong lf    = i + delta + res;

    ctrl.active_threads = nworkers;

    do {
      NcmFuncEvalLoopEval *arg = g_slice_new (NcmFuncEvalLoopEval);

      arg->lfunc = lfunc;
      arg->i     = li;
      arg->f     = lf;
      arg->data  = data;
      arg->ctrl  = &ctrl;
      g_thread_pool_push (_function_thread_pool, arg, &err);
      li  = lf;
      lf += delta;
    } while (--nworkers);
  }

  g_mutex_lock (&ctrl.update);

  while (ctrl.active_threads != 0)
    g_cond_wait (&ctrl.finish, &ctrl.update);

  g_mutex_unlock (&ctrl.update);

  g_mutex_clear (&ctrl.update);
  g_cond_clear (&ctrl.finish);
}

/**
 * ncm_func_eval_threaded_loop:
 * @lfunc: (scope notified): loop function
 * @i: first index
 * @f: one past the last index
 * @data: user data
 *
 * Calls ncm_func_eval_threaded_loop_nw() with one worker per pool thread, one per
 * processor for an unlimited pool, and serially for a pool of zero threads.
 */
void
ncm_func_eval_threaded_loop (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data)
{
  ncm_func_eval_get_pool ();
  {
    const gint max_threads = g_thread_pool_get_max_threads (_function_thread_pool);
    const guint nworkers   = (max_threads > 0) ? (guint) max_threads : ((max_threads < 0) ? g_get_num_processors () : 1);

    ncm_func_eval_threaded_loop_nw (lfunc, i, f, data, nworkers);
  }
}

/**
 * ncm_func_eval_threaded_loop_full:
 * @lfunc: (scope notified): loop function
 * @i: first index
 * @f: one past the last index
 * @data: user data
 *
 * Calls @lfunc with @data on each index of [@i, @f) as a separate pool task
 * and returns after all finish. Runs a single serial call if the pool has zero
 * threads.
 */
#if NCM_THREAD_POOL_MAX > 1

void
ncm_func_eval_threaded_loop_full (NcmFuncEvalLoop lfunc, glong i, glong f, gpointer data)
{
  NcmFuncEvalCtrl ctrl = {0};

  ncm_func_eval_get_pool ();
  g_mutex_init (&ctrl.update);
  g_cond_init (&ctrl.finish);

  g_assert_cmpuint (f, >, i);

  if (g_thread_pool_get_max_threads (_function_thread_pool) > 0)
  {
    GError *err = NULL;
    glong l;

    ctrl.active_threads = f - i;

    for (l = i; l < f; l++)
    {
      NcmFuncEvalLoopEval *arg = g_slice_new (NcmFuncEvalLoopEval);

      arg->lfunc = lfunc;
      arg->i     = l;
      arg->f     = l + 1;
      arg->data  = data;
      arg->ctrl  = &ctrl;
      g_thread_pool_push (_function_thread_pool, arg, &err);
    }
  }
  else
  {
    lfunc (i, f, data);
  }

  g_mutex_lock (&ctrl.update);

  while (ctrl.active_threads != 0)
    g_cond_wait (&ctrl.finish, &ctrl.update);

  g_mutex_unlock (&ctrl.update);

  g_mutex_clear (&ctrl.update);
  g_cond_clear (&ctrl.finish);
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


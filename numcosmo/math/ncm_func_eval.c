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
 * @title: Function Evaluator
 * @short_description: A general purpose multi-threaded function evaluator
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_func_eval.h"
#include "math/ncm_cfg.h"
#include <stdio.h>

typedef struct _NcmFunctionEvalCtrl
{
  gint active_threads;
  GMutex *update;
  GCond *finish;
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32))
  GMutex update_m;
  GCond finish_c;
#endif
} NcmFunctionEvalCtrl;

typedef struct _NcmLoopFuncEval
{
  NcmLoopFunc lfunc;
  glong i;
  glong f;
  gpointer data;
  NcmFunctionEvalCtrl *ctrl;
} NcmLoopFuncEval;

static GThreadPool *_function_thread_pool = NULL;

static void
func (gpointer data, gpointer empty)
{
  NcmLoopFuncEval *arg = (NcmLoopFuncEval *)data;
  NcmFunctionEvalCtrl *ctrl = arg->ctrl;

  arg->lfunc (arg->i, arg->f, arg->data);
  g_slice_free (NcmLoopFuncEval, arg);

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
	g_clear_error (&err);
  }
  _NCM_MUTEX_UNLOCK (&create_lock);

  return _function_thread_pool;
}

/**
 * ncm_func_eval_set_max_threads:
 * @mt: new max threads to be used in the pool, -1 means unlimited
 *
 * Set the new maximun number of threads to be used by the pool
 *
 */
void
ncm_func_eval_set_max_threads (gint mt)
{
  GError *err = NULL;
  ncm_func_eval_get_pool ();
  g_thread_pool_set_max_threads (_function_thread_pool, mt, &err);
}

/**
 * ncm_func_eval_threaded_loop:
 * @lfunc: (scope notified): #NcmLoopFunc to be evaluated in threads
 * @i: initial index
 * @f: final index
 * @data: pointer to be passed to @fl
 *
 * Using the thread pool, evaluate @fl in each value of (@f-@i)/nthreads
 *
 */
#if NCM_THREAD_POOL_MAX > 1
void
ncm_func_eval_threaded_loop (NcmLoopFunc lfunc, glong i, glong f, gpointer data)
{
  NcmFunctionEvalCtrl ctrl = {0, NULL, NULL};
  gint nthreads, delta, res;

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

  g_assert (f >= i);
  nthreads = g_thread_pool_get_max_threads (_function_thread_pool);
  delta = (f-i) / nthreads;
  res = (f-i) % nthreads;

  if (delta == 0)
  {
	lfunc (i, f, data);
  }
  else
  {
	GError *err = NULL;
	glong li = i;
	glong lf = delta + res;
	ctrl.active_threads = nthreads;

	do {
	  NcmLoopFuncEval *arg = g_slice_new (NcmLoopFuncEval);
	  arg->lfunc = lfunc;
	  arg->i = li;
	  arg->f = lf;
	  arg->data = data;
	  arg->ctrl = &ctrl;
	  g_thread_pool_push (_function_thread_pool, arg, &err);
	  li = lf;
	  lf += delta;
	} while (--nthreads);
  }

  g_mutex_lock (ctrl.update);
  while (ctrl.active_threads != 0)
	g_cond_wait (ctrl.finish, ctrl.update);
  g_mutex_unlock (ctrl.update);

  if (FALSE)
  {
	printf  ("Unused:      %d\n", g_thread_pool_get_num_unused_threads ());fflush (stdout);
	printf  ("Max Unused:  %d\n", g_thread_pool_get_max_unused_threads ());fflush (stdout);
	printf  ("Running:     %d\n", g_thread_pool_get_num_threads (_function_thread_pool));fflush (stdout);
	printf  ("Unprocessed: %d\n", g_thread_pool_unprocessed (_function_thread_pool));fflush (stdout);
	printf  ("Unused:      %d\n", g_thread_pool_get_max_threads (_function_thread_pool));fflush (stdout);
  }

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
ncm_func_eval_threaded_loop (NcmLoopFunc lfunc, glong i, glong f, gpointer data)
{
  lfunc (i, f, data);
}
#endif

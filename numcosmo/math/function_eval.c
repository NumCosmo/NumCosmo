/***************************************************************************
 *            function_eval.c
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
 * SECTION:function_eval
 * @title: Function Evaluator
 * @short_description: A general purpose multi-threaded function evaluator
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>

typedef struct _NcmFunctionEvalCtrl
{
  gint active_threads;
  gsl_function *F;
  GMutex *update;
  GCond *finish;
} NcmFunctionEvalCtrl;

typedef struct _NcmFunctionEvalArg
{
  gdouble *x;
	gdouble *val;
  NcmFunctionEvalCtrl *ctrl;
} NcmFunctionEvalArg;

static GThreadPool *_function_thread_pool = NULL;

static void
func (gpointer data, gpointer empty)
{
  NcmFunctionEvalArg *arg = (NcmFunctionEvalArg *)data;
	arg->val[0] = GSL_FN_EVAL (arg->ctrl->F, arg->x[0]);
	
  g_mutex_lock (arg->ctrl->update);
  arg->ctrl->active_threads--;
  g_slice_free (NcmFunctionEvalArg, arg);
  if (arg->ctrl->active_threads == 0)
    g_cond_signal (arg->ctrl->finish);
  g_mutex_unlock (arg->ctrl->update);
	
  return;
}

/**
 * ncm_function_eval_get_pool: (skip)
 * 
 * Allocate if its not yet allocated and return the
 * internal GThreadPool pool.
 * 
 * Returns: the pointer to the internal GThreadPool pool
 */
GThreadPool *
ncm_function_eval_get_pool ()
{
  static GStaticMutex create_lock = G_STATIC_MUTEX_INIT;
  GError *err = NULL;
  
  g_static_mutex_lock (&create_lock);
  if (_function_thread_pool == NULL)
  {
    _function_thread_pool = g_thread_pool_new (func, NULL, NC_THREAD_POOL_MAX, TRUE, &err);
    g_clear_error (&err);
  }
  g_static_mutex_unlock (&create_lock);

  return _function_thread_pool;
}

/**
 * ncm_function_eval_set_max_threads:
 * @mt: new max threads to be used in the pool, -1 means unlimited
 * 
 * Set the new maximun number of threads to be used by the pool
 * 
 */
void
ncm_function_eval_set_max_threads (gint mt)
{
  GError *err = NULL;
  ncm_function_eval_get_pool ();
  g_thread_pool_set_max_threads (_function_thread_pool, mt, &err);
}

/**
 * ncm_function_eval_threaded: (skip)
 * @F: gsl_function to be evaluated in threads
 * @x: array of values to evaluate the function
 * @val: array to store the values of F evaluated in x
 * @n: number of elements in x
 * @x_stride: the space between elements in the x array
 * @val_stride: the space between elements in the val array
 * 
 * Using the thread pool, evaluate the gsl_function in each value of 
 * the array x[x_stride * i] and stores the result in val[val_stride * i]
 * for i in [0, n-1]. Note that F must contain a reentrant function.
 * 
 */
void
ncm_function_eval_threaded (gsl_function *F, gdouble *x, gdouble *val, gulong n, guint x_stride, guint val_stride)
{
	NcmFunctionEvalCtrl ctrl = {n, F, NULL, NULL};
  gint i;

  ncm_function_eval_get_pool ();
  ctrl.update = g_mutex_new ();
  ctrl.finish = g_cond_new ();

#if NC_THREAD_POOL_MAX > 1
  {
    GError *err = NULL;
    for (i = 0; i < n; i++)
    {
      NcmFunctionEvalArg *arg = g_slice_new (NcmFunctionEvalArg);
      arg->x = &x[i * x_stride];
      arg->val = &val[i * x_stride];
      arg->ctrl = &ctrl;
      g_thread_pool_push (_function_thread_pool, arg, &err);
    }
  }
#else
  for (i = 0; i < n; i++)
    val[i * x_stride] = GSL_FN_EVAL (F, x[i * x_stride]);
#endif

  printf ("Waiting for the sun\n");fflush (stdout);
  g_mutex_lock (ctrl.update);
	while (ctrl.active_threads != 0)
    g_cond_wait (ctrl.finish, ctrl.update);
  g_mutex_unlock (ctrl.update);
  
  printf  ("Unused:      %d\n", g_thread_pool_get_num_unused_threads ());fflush (stdout);
  printf  ("Max Unused:  %d\n", g_thread_pool_get_max_unused_threads ());fflush (stdout);
  printf  ("Running:     %d\n", g_thread_pool_get_num_threads (_function_thread_pool));fflush (stdout);
  printf  ("Unprocessed: %d\n", g_thread_pool_unprocessed (_function_thread_pool));fflush (stdout);
  printf  ("Unused:      %d\n", g_thread_pool_get_max_threads (_function_thread_pool));fflush (stdout);

  g_mutex_free (ctrl.update);
  g_cond_free (ctrl.finish);
}

/**
 * ncm_function_eval_threaded_vec: (skip)
 * @F: gsl_function to be evaluated in threads
 * @x: gsl_vector of values to evaluate the function
 * @val: gsl_vector to store the values of F evaluated in x
 * 
 * Using the thread pool, evaluate the gsl_function in each value of 
 * the gsl_vector x and stores the result in the gsl_vector val
 * for the x->size values of x. Note that F must contain a reentrant function.
 * 
 */
void
ncm_function_eval_threaded_vec (gsl_function *F, gsl_vector *x, gsl_vector *val)
{
  g_assert (x->size >= val->size);
  ncm_function_eval_threaded (F, x->data, val->data, x->size, x->stride, val->stride);
}

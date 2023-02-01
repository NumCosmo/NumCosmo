/***************************************************************************
 *            ncm_func_eval.c
 *
 *  Fri April 10 16:25:22 2015
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2015 <vitenti@uel.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcmSparam
{
  guint ntests;
} TestNcmSparam;

void test_ncm_func_eval_new (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_func_eval_free (TestNcmSparam *test, gconstpointer pdata);

void test_ncm_func_eval_run (TestNcmSparam *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/ncm/func_eval/run", TestNcmSparam, NULL,
              &test_ncm_func_eval_new,
              &test_ncm_func_eval_run,
              &test_ncm_func_eval_free);
  
  g_test_run ();
}

void
test_ncm_func_eval_new (TestNcmSparam *test, gconstpointer pdata)
{
  test->ntests = 10000;
}

void
test_ncm_func_eval_free (TestNcmSparam *test, gconstpointer pdata)
{
}

void
test_ncm_func_eval_run_func (glong i, glong f, gpointer data)
{
  G_LOCK_DEFINE_STATIC (save_data);
  glong k;
  gdouble *res = (gdouble *) data;
  gdouble part = 0.0;
  
  for (k = 0; k < 1000; k++)
  {
    gdouble v;
    
    v     = cos (k + M_PI * 8.9);
    v     = sin (part);
    v     = exp (log (fabs (part)) * 0.9);
    part += v;
  }
  
  G_LOCK (save_data);
  *res += part;
  G_UNLOCK (save_data);
}

void
test_ncm_func_eval_run (TestNcmSparam *test, gconstpointer pdata)
{
  gdouble res = 0.0;
  
  ncm_func_eval_threaded_loop_full (test_ncm_func_eval_run_func, 0, test->ntests, &res);
}


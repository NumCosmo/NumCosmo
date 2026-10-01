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
void test_ncm_func_eval_threaded_loop_cover (void);

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
  g_test_add_func ("/ncm/func_eval/threaded_loop/cover", &test_ncm_func_eval_threaded_loop_cover);

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

static void
_test_ncm_func_eval_count (glong i, glong f, gpointer data)
{
  gint *count = data;
  glong k;

  for (k = i; k < f; k++)
    g_atomic_int_inc (&count[k]);
}

/* Every index is visited once, for any pool size, including unlimited and zero, and any length */
void
test_ncm_func_eval_threaded_loop_cover (void)
{
  const gint max_threads[] = {-1, 0, 1, 3, NCM_THREAD_POOL_MAX};
  const glong lengths[]    = {1, 2, 3, 7, 100};
  guint t, l;

  for (t = 0; t < G_N_ELEMENTS (max_threads); t++)
  {
    ncm_func_eval_set_max_threads (max_threads[t]);

    for (l = 0; l < G_N_ELEMENTS (lengths); l++)
    {
      gint *count = g_new0 (gint, lengths[l]);
      guint pass;

      for (pass = 0; pass < 3; pass++)
      {
        glong k;

        if (pass == 0)
          ncm_func_eval_threaded_loop (&_test_ncm_func_eval_count, 0, lengths[l], count);
        else if (pass == 1)
          ncm_func_eval_threaded_loop_nw (&_test_ncm_func_eval_count, 0, lengths[l], count, 8);
        else
          ncm_func_eval_threaded_loop_full (&_test_ncm_func_eval_count, 0, lengths[l], count);

        for (k = 0; k < lengths[l]; k++)
          g_assert_cmpint (count[k], ==, pass + 1);
      }

      g_free (count);
    }
  }

  ncm_func_eval_set_max_threads (NCM_THREAD_POOL_MAX);
}


/***************************************************************************
 *            test_ncm_timer.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_timer.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

static void
test_ncm_timer_elapsed (void)
{
  NcmTimer *nt = ncm_timer_new ();
  guint day, hour, min;
  gdouble sec, elapsed;

  ncm_timer_start (nt);
  g_usleep (2000);
  ncm_timer_stop (nt);

  /* A stopped timer does not advance */
  elapsed = ncm_timer_elapsed (nt);
  g_assert_cmpfloat (elapsed, >=, 2.0e-3);
  g_usleep (2000);
  g_assert_cmpfloat (ncm_timer_elapsed (nt), ==, elapsed);

  ncm_timer_elapsed_dhms (nt, &day, &hour, &min, &sec);
  g_assert_cmpuint (day, ==, 0);
  g_assert_cmpuint (hour, ==, 0);
  g_assert_cmpuint (min, ==, 0);
  g_assert_cmpfloat (sec, ==, elapsed);
  {
    /* From the measured time: a 2 ms sleep can last more than 10 ms on a busy machine */
    gchar *expected = g_strdup_printf ("00:00:" NCM_TIMER_SEC_FORMAT, elapsed);

    g_assert_cmpstr (ncm_timer_elapsed_dhms_str (nt), ==, expected);
    g_free (expected);
  }

  ncm_timer_continue (nt);
  g_usleep (2000);
  g_assert_cmpfloat (ncm_timer_elapsed (nt), >, elapsed);

  ncm_timer_free (nt);
}

static void
test_ncm_timer_task (void)
{
  NcmTimer *nt = ncm_timer_new ();
  guint task_len, task_pos;
  gchar *name;

  ncm_timer_set_name (nt, "test-task");
  g_object_get (nt, "name", &name, NULL);
  g_assert_cmpstr (name, ==, "test-task");
  g_free (name);

  g_assert_false (ncm_timer_task_is_running (nt));

  ncm_timer_task_start (nt, 5);
  g_assert_true (ncm_timer_task_is_running (nt));

  ncm_timer_task_increment (nt);
  ncm_timer_task_accumulate (nt, 2);
  g_assert_cmpuint (ncm_timer_task_completed (nt), ==, 3);
  g_assert_false (ncm_timer_task_has_ended (nt));

  g_object_get (nt, "task-len", &task_len, "task-pos", &task_pos, NULL);
  g_assert_cmpuint (task_len, ==, 5);
  g_assert_cmpuint (task_pos, ==, 3);

  g_assert_true (g_strrstr (ncm_timer_task_elapsed_str (nt), "test-task, completed: 3 of 5") != NULL);

  ncm_timer_task_add_tasks (nt, 2);
  ncm_timer_task_accumulate (nt, 4);
  g_assert_true (ncm_timer_task_has_ended (nt));
  g_assert_true (ncm_timer_task_end (nt));
  g_assert_false (ncm_timer_task_is_running (nt));

  /* Ending before all items are completed */
  ncm_timer_task_start (nt, 3);
  ncm_timer_task_increment (nt);
  g_assert_false (ncm_timer_task_end (nt));

  ncm_timer_free (nt);
}

static void
test_ncm_timer_task_estimates (void)
{
  NcmTimer *nt = ncm_timer_new ();
  gdouble mean_time;
  guint i;

  ncm_timer_task_start (nt, 10);

  for (i = 0; i < 4; i++)
  {
    g_usleep (2000);
    ncm_timer_task_increment (nt);
  }

  mean_time = ncm_timer_task_mean_time (nt);
  g_assert_cmpfloat (mean_time, >=, 2.0e-3);
  g_assert_cmpfloat (ncm_timer_task_time_left (nt), ==, 6.0 * mean_time);
  g_assert_cmpuint (ncm_timer_task_estimate_by_time (nt, 2.5 * mean_time), ==, 3);

  /* A paused task does not advance */
  ncm_timer_task_pause (nt);
  {
    const gdouble elapsed = ncm_timer_elapsed (nt);

    g_usleep (2000);
    g_assert_cmpfloat (ncm_timer_elapsed (nt), ==, elapsed);
  }
  ncm_timer_task_continue (nt);

  g_assert_true (g_strrstr (ncm_timer_task_mean_time_str (nt), "mean time: ") != NULL);
  g_assert_true (g_strrstr (ncm_timer_task_time_left_str (nt), "time left: ") != NULL);
  g_assert_true (g_strrstr (ncm_timer_task_start_datetime_str (nt), "started at: ") != NULL);
  g_assert_true (g_strrstr (ncm_timer_task_end_datetime_str (nt), "estimated to end at: ") != NULL);
  g_assert_true (g_strrstr (ncm_timer_task_cur_datetime_str (nt), "current time: ") != NULL);

  /* Logging records the time of the last completed item */
  ncm_timer_task_log_elapsed (nt);
  ncm_timer_task_log_mean_time (nt);
  ncm_timer_task_log_time_left (nt);
  ncm_timer_task_log_start_datetime (nt);
  ncm_timer_task_log_cur_datetime (nt);
  ncm_timer_task_log_end_datetime (nt);
  g_assert_cmpfloat (ncm_timer_elapsed_since_last_log (nt), >=, 0.0);
  g_assert_cmpfloat (ncm_timer_elapsed_since_last_log (nt), <=, ncm_timer_elapsed (nt));

  ncm_timer_free (nt);
}

static void
test_ncm_timer_task_start_zero_subprocess (void)
{
  NcmTimer *nt = ncm_timer_new ();

  ncm_timer_task_start (nt, 0);
}

static void
test_ncm_timer_start_during_task_subprocess (void)
{
  NcmTimer *nt = ncm_timer_new ();

  ncm_timer_task_start (nt, 2);
  ncm_timer_start (nt);
}

static void
test_ncm_timer_increment_past_end_subprocess (void)
{
  NcmTimer *nt = ncm_timer_new ();

  ncm_timer_task_start (nt, 1);
  ncm_timer_task_increment (nt);
  ncm_timer_task_increment (nt);
}

static void
test_ncm_timer_traps (void)
{
  g_test_trap_subprocess ("/ncm/timer/task_start_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/timer/start_during_task/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/timer/increment_past_end/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

/* Estimates print "unknown" until defined, never nan */
static void
test_ncm_timer_str_undefined (void)
{
  NcmTimer *nt = ncm_timer_new ();
  guint n;

  ncm_timer_task_start (nt, 10);

  for (n = 0; n < 3; n++)
  {
    const gchar *strs[3];
    guint k;

    strs[0] = g_strdup (ncm_timer_task_mean_time_str (nt));
    strs[1] = g_strdup (ncm_timer_task_time_left_str (nt));
    strs[2] = g_strdup (ncm_timer_task_end_datetime_str (nt));

    for (k = 0; k < 3; k++)
    {
      g_assert_null (g_strstr_len (strs[k], -1, "nan"));

      if (n == 0)
        g_assert_nonnull (g_strstr_len (strs[k], -1, "unknown"));
      else if (n == 2)
        g_assert_null (g_strstr_len (strs[k], -1, "unknown"));

      g_free ((gchar *) strs[k]);
    }

    ncm_timer_task_increment (nt);
  }

  ncm_timer_free (nt);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/timer/str_undefined", &test_ncm_timer_str_undefined);
  g_test_add_func ("/ncm/timer/elapsed", &test_ncm_timer_elapsed);
  g_test_add_func ("/ncm/timer/task", &test_ncm_timer_task);
  g_test_add_func ("/ncm/timer/task_estimates", &test_ncm_timer_task_estimates);
  g_test_add_func ("/ncm/timer/traps", &test_ncm_timer_traps);
  g_test_add_func ("/ncm/timer/task_start_zero/subprocess", &test_ncm_timer_task_start_zero_subprocess);
  g_test_add_func ("/ncm/timer/start_during_task/subprocess", &test_ncm_timer_start_during_task_subprocess);
  g_test_add_func ("/ncm/timer/increment_past_end/subprocess", &test_ncm_timer_increment_past_end_subprocess);

  g_test_run ();
}


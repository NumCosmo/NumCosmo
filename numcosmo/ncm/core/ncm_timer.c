/***************************************************************************
 *            ncm_timer.c
 *
 *  Thu August 01 11:22:51 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_timer.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmTimer:
 *
 * Wall-clock timer with progress tracking for a task of a known number of items.
 *
 * A task is started with ncm_timer_task_start() and advanced with
 * ncm_timer_task_increment() or ncm_timer_task_accumulate(). The time per item is
 * accumulated to give the mean time per item, the time left and the estimated end
 * time. Functions returning a string return a buffer owned by the timer, which the
 * next call to any of them overwrites.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_timer.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"
#include "ncm/stats/ncm_stats_vec.h"

enum
{
  PROP_0,
  PROP_NAME,
  PROP_TASK_LEN,
  PROP_TASK_POS,
  PROP_SIZE,
};

struct _NcmTimer
{
  /*< private >*/
  GObject parent_instance;
  GTimer *gt;
  gchar *name;
  guint task_len;
  guint task_pos;
  gdouble pos_time;
  gdouble last_log_time;
  NcmStatsVec *time_stats;
  GString *msg;
  GString *msg_tmp1;
  GString *msg_tmp2;
};

G_DEFINE_TYPE (NcmTimer, ncm_timer, G_TYPE_OBJECT)

#define _NCM_TIMER_MSG_PREALLOC_SIZE 100

static void
ncm_timer_init (NcmTimer *nt)
{
  nt->gt         = g_timer_new ();
  nt->msg        = g_string_sized_new (_NCM_TIMER_MSG_PREALLOC_SIZE);
  nt->msg_tmp1   = g_string_sized_new (_NCM_TIMER_MSG_PREALLOC_SIZE);
  nt->msg_tmp2   = g_string_sized_new (_NCM_TIMER_MSG_PREALLOC_SIZE);
  nt->name       = NULL;
  nt->task_len   = 0;
  nt->task_pos   = 0;
  nt->pos_time   = 0.0;
  nt->time_stats = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
}

static void
_ncm_timer_dispose (GObject *object)
{
  NcmTimer *nt = NCM_TIMER (object);

  ncm_g_string_clear (&nt->msg);
  ncm_g_string_clear (&nt->msg_tmp1);
  ncm_g_string_clear (&nt->msg_tmp2);
  ncm_stats_vec_clear (&nt->time_stats);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_timer_parent_class)->dispose (object);
}

static void
_ncm_timer_finalize (GObject *object)
{
  NcmTimer *nt = NCM_TIMER (object);

  if (nt->gt != NULL)
  {
    g_timer_destroy (nt->gt);
    nt->gt = NULL;
  }

  g_free (nt->name);
  nt->name = NULL;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_timer_parent_class)->finalize (object);
}

static void
_ncm_timer_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmTimer *nt = NCM_TIMER (object);

  g_return_if_fail (NCM_IS_TIMER (object));

  switch (prop_id)
  {
    case PROP_NAME:
      ncm_timer_set_name (nt, g_value_get_string (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_timer_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmTimer *nt = NCM_TIMER (object);

  g_return_if_fail (NCM_IS_TIMER (object));

  switch (prop_id)
  {
    case PROP_NAME:
      g_value_set_string (value, nt->name);
      break;
    case PROP_TASK_LEN:
      g_value_set_uint (value, nt->task_len);
      break;
    case PROP_TASK_POS:
      g_value_set_uint (value, nt->task_pos);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_timer_class_init (NcmTimerClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = &_ncm_timer_dispose;
  object_class->finalize     = &_ncm_timer_finalize;
  object_class->set_property = &_ncm_timer_set_property;
  object_class->get_property = &_ncm_timer_get_property;

  /**
   * NcmTimer:name:
   *
   * The name used in the task messages.
   */
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Timer's name",
                                                        "timer",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmTimer:task-len:
   *
   * The number of items of the running task, zero if no task is running.
   */
  g_object_class_install_property (object_class,
                                   PROP_TASK_LEN,
                                   g_param_spec_uint ("task-len",
                                                      NULL,
                                                      "Task length",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmTimer:task-pos:
   *
   * The number of completed items of the task, in [0, `NcmTimer:task-len`].
   */
  g_object_class_install_property (object_class,
                                   PROP_TASK_POS,
                                   g_param_spec_uint ("task-pos",
                                                      NULL,
                                                      "Task position",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_timer_new:
 *
 * Creates a new #NcmTimer.
 *
 * Returns: (transfer full): a new #NcmTimer.
 */
NcmTimer *
ncm_timer_new (void)
{
  NcmTimer *nt = g_object_new (NCM_TYPE_TIMER, NULL);

  return nt;
}

/**
 * ncm_timer_ref:
 * @nt: a #NcmTimer
 *
 * Increases the reference count of @nt by one.
 *
 * Returns: (transfer full): @nt.
 */
NcmTimer *
ncm_timer_ref (NcmTimer *nt)
{
  return g_object_ref (nt);
}

/**
 * ncm_timer_free:
 * @nt: a #NcmTimer
 *
 * Decreases the reference count of @nt by one.
 */
void
ncm_timer_free (NcmTimer *nt)
{
  g_object_unref (nt);
}

/**
 * ncm_timer_clear:
 * @nt: a #NcmTimer
 *
 * Decreases the reference count of *@nt by one and sets *@nt to %NULL.
 */
void
ncm_timer_clear (NcmTimer **nt)
{
  g_clear_object (nt);
}

/**
 * ncm_timer_set_name:
 * @nt: a #NcmTimer
 * @name: the name
 *
 * Sets the name of @nt, used in the task messages.
 */
void
ncm_timer_set_name (NcmTimer *nt, const gchar *name)
{
  if (nt->name != NULL)
    g_free (nt->name);

  nt->name = g_strdup (name);
}

/**
 * ncm_timer_elapsed:
 * @nt: a #NcmTimer
 *
 * Returns: the seconds elapsed since @nt was started, up to the moment it was
 * stopped if it is stopped.
 */
gdouble
ncm_timer_elapsed (NcmTimer *nt)
{
  return g_timer_elapsed (nt->gt, NULL);
}

static void
_ncm_timer_sec_to_dhms (gdouble t, guint *elap_day, guint *elap_hour, guint *elap_min, gdouble *elap_sec)
{
  *elap_sec  = t;
  *elap_min  = *elap_sec / 60;
  *elap_hour = *elap_min / 60;
  *elap_day  = *elap_hour / 24;

  *elap_sec  = fmod (*elap_sec, 60);
  *elap_min  = *elap_min % 60;
  *elap_hour = *elap_hour % 24;
}

static void
_ncm_timer_dhms_to_string (GString *s, guint elap_day, guint elap_hour, guint elap_min, gdouble elap_sec)
{
  switch (elap_day)
  {
    case 0:
      g_string_printf (s, "%02u:%02u:"NCM_TIMER_SEC_FORMAT, elap_hour, elap_min, elap_sec);
      break;
    case 1:
      g_string_printf (s, "1 day, %02u:%02u:"NCM_TIMER_SEC_FORMAT, elap_hour, elap_min, elap_sec);
      break;
    default:
      g_string_printf (s, "%02u days, %02u:%02u:"NCM_TIMER_SEC_FORMAT, elap_day, elap_hour, elap_min, elap_sec);
      break;
  }
}

/**
 * ncm_timer_elapsed_dhms:
 * @nt: a #NcmTimer
 * @elap_day: (out): days
 * @elap_hour: (out): hours
 * @elap_min: (out): minutes
 * @elap_sec: (out): seconds
 *
 * Splits ncm_timer_elapsed() into days, hours, minutes and seconds.
 */
void
ncm_timer_elapsed_dhms (NcmTimer *nt, guint *elap_day, guint *elap_hour, guint *elap_min, gdouble *elap_sec)
{
  gdouble elapsed = ncm_timer_elapsed (nt);

  _ncm_timer_sec_to_dhms (elapsed, elap_day, elap_hour, elap_min, elap_sec);
}

/**
 * ncm_timer_elapsed_dhms_str:
 * @nt: a #NcmTimer
 *
 * Formats ncm_timer_elapsed() as "[D days, ]HH:MM:SS.ssss".
 *
 * Returns: (transfer none): the formatted string.
 */
gchar *
ncm_timer_elapsed_dhms_str (NcmTimer *nt)
{
  guint elap_day, elap_hour, elap_min;
  gdouble elap_sec;

  ncm_timer_elapsed_dhms (nt, &elap_day, &elap_hour, &elap_min, &elap_sec);
  _ncm_timer_dhms_to_string (nt->msg, elap_day, elap_hour, elap_min, elap_sec);

  return nt->msg->str;
}

/**
 * ncm_timer_start:
 * @nt: a #NcmTimer
 *
 * Starts @nt, resetting the elapsed time. Aborts if a task is running.
 */
void
ncm_timer_start (NcmTimer *nt)
{
  if (nt->task_len != 0)
    g_error ("ncm_timer_start: cannot start timer during a task, call task_end first.");

  g_timer_start (nt->gt);
}

/**
 * ncm_timer_stop:
 * @nt: a #NcmTimer
 *
 * Stops @nt. Aborts if a task is running.
 */
void
ncm_timer_stop (NcmTimer *nt)
{
  if (nt->task_len != 0)
    g_error ("ncm_timer_stop: cannot end timer during a task, call task_end first.");

  g_timer_stop (nt->gt);
}

/**
 * ncm_timer_continue:
 * @nt: a #NcmTimer
 *
 * Resumes @nt after ncm_timer_stop(). Aborts if a task is running.
 */
void
ncm_timer_continue (NcmTimer *nt)
{
  if (nt->task_len != 0)
    g_error ("ncm_timer_continue: cannot continue timer during a task, call task_end first.");

  g_timer_continue (nt->gt);
}

/**
 * ncm_timer_task_start:
 * @nt: a #NcmTimer
 * @task_len: number of items
 *
 * Starts a task of @task_len items and restarts the timer. Aborts if a task is
 * running or if @task_len is zero.
 */
void
ncm_timer_task_start (NcmTimer *nt, guint task_len)
{
  if (nt->task_len != 0)
    g_error ("ncm_timer_task_start: cannot start a new task during a task, call task_end first.");

  if (task_len == 0)
    g_error ("ncm_timer_task_start: cannot start task with 0 items.");

  ncm_timer_start (nt);
  nt->task_len      = task_len;
  nt->task_pos      = 0;
  nt->pos_time      = g_timer_elapsed (nt->gt, NULL);
  nt->last_log_time = 0.0;
  ncm_stats_vec_reset (nt->time_stats, TRUE);
}

/**
 * ncm_timer_task_increment:
 * @nt: a #NcmTimer
 *
 * Marks one item as completed, recording the time since the previous item. Aborts
 * if no task is running or if the task is already complete.
 */
void
ncm_timer_task_increment (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  {
    const gdouble tot_elapsed = ncm_timer_elapsed (nt);
    const gdouble elap_sec    = tot_elapsed - nt->pos_time;

    nt->pos_time = tot_elapsed;
    nt->task_pos++;
    ncm_stats_vec_set (nt->time_stats, 0, elap_sec);
    ncm_stats_vec_update (nt->time_stats);

    if (nt->task_pos > nt->task_len)
      g_error ("ncm_timer_task_increment: incrementing past the end of the task.");
  }
}

/**
 * ncm_timer_task_accumulate:
 * @nt: a #NcmTimer
 * @nitens: number of items
 *
 * Marks @nitens items as completed, recording the time since the previous update
 * divided by @nitens as the time of each. Aborts if no task is running or if this
 * passes the end of the task.
 */
void
ncm_timer_task_accumulate (NcmTimer *nt, guint nitens)
{
  g_assert (nt->task_len != 0);
  {
    const gdouble tot_elapsed = ncm_timer_elapsed (nt);
    const gdouble elap_sec    = tot_elapsed - nt->pos_time;

    nt->pos_time  = tot_elapsed;
    nt->task_pos += nitens;

    ncm_stats_vec_set (nt->time_stats, 0, elap_sec / nitens);
    ncm_stats_vec_update_weight (nt->time_stats, nitens);

    if (nt->task_pos > nt->task_len)
      g_error ("ncm_timer_task_accumulate: incrementing past the end of the task.");
  }
}

/**
 * ncm_timer_task_completed:
 * @nt: a #NcmTimer
 *
 * Aborts if no task is running.
 *
 * Returns: the number of completed items.
 */
guint
ncm_timer_task_completed (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);

  return nt->task_pos;
}

/**
 * ncm_timer_task_estimate_by_time:
 * @nt: a #NcmTimer
 * @sec: time in seconds
 *
 * Returns: the number of items that take @sec seconds at the mean time per item,
 * rounded up.
 */
guint
ncm_timer_task_estimate_by_time (NcmTimer *nt, gdouble sec)
{
  return ceil (sec / ncm_timer_task_mean_time (nt));
}

/**
 * ncm_timer_task_pause:
 * @nt: a #NcmTimer
 *
 * Stops the timer during a task. Aborts if no task is running.
 */
void
ncm_timer_task_pause (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  g_timer_stop (nt->gt);
}

/**
 * ncm_timer_task_continue:
 * @nt: a #NcmTimer
 *
 * Resumes the timer after ncm_timer_task_pause(). Aborts if no task is running.
 */
void
ncm_timer_task_continue (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  g_timer_continue (nt->gt);
}

/**
 * ncm_timer_task_add_tasks:
 * @nt: a #NcmTimer
 * @ptasks: number of items
 *
 * Adds @ptasks items to the running task. Aborts if no task is running.
 */
void
ncm_timer_task_add_tasks (NcmTimer *nt, guint ptasks)
{
  g_assert (nt->task_len != 0);
  nt->task_len += ptasks;
}

/**
 * ncm_timer_task_is_running:
 * @nt: a #NcmTimer
 *
 * Returns: whether a task is running.
 */
gboolean
ncm_timer_task_is_running (NcmTimer *nt)
{
  return (nt->task_len != 0);
}

/**
 * ncm_timer_task_has_ended:
 * @nt: a #NcmTimer
 *
 * Aborts if no task is running.
 *
 * Returns: whether all items of the task are completed.
 */
gboolean
ncm_timer_task_has_ended (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  {
    gboolean ok = nt->task_len == nt->task_pos;

    return ok;
  }
}

/**
 * ncm_timer_task_end:
 * @nt: a #NcmTimer
 *
 * Ends the running task. Aborts if no task is running.
 *
 * Returns: whether all items of the task were completed.
 */
gboolean
ncm_timer_task_end (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  {
    gboolean ok = nt->task_len == nt->task_pos;

    nt->task_len = 0;

    return ok;
  }
}

/**
 * ncm_timer_elapsed_since_last_log:
 * @nt: a #NcmTimer
 *
 * The ncm_timer_task_log functions record the time of the last completed item as
 * the time of the log.
 *
 * Returns: the seconds elapsed since the last log.
 */
gdouble
ncm_timer_elapsed_since_last_log (NcmTimer *nt)
{
  return ncm_timer_elapsed (nt) - nt->last_log_time;
}

/**
 * ncm_timer_task_mean_time:
 * @nt: a #NcmTimer
 *
 * Returns: the mean time per completed item, in seconds.
 */
gdouble
ncm_timer_task_mean_time (NcmTimer *nt)
{
  return ncm_stats_vec_get_mean (nt->time_stats, 0);
}

/**
 * ncm_timer_task_time_left:
 * @nt: a #NcmTimer
 *
 * Returns: the mean time per item times the number of items left, in seconds.
 */
gdouble
ncm_timer_task_time_left (NcmTimer *nt)
{
  const gdouble mean_time = ncm_timer_task_mean_time (nt);
  const guint task_left   = nt->task_len - nt->task_pos;

  return mean_time * task_left;
}

/**
 * ncm_timer_task_elapsed_str:
 * @nt: a #NcmTimer
 *
 * Formats the number of completed items and the elapsed time. Aborts if no task is
 * running.
 *
 * Returns: (transfer none): the formatted string.
 */
const gchar *
ncm_timer_task_elapsed_str (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  {
    guint elap_day, elap_hour, elap_min;
    gdouble elap_sec;

    ncm_timer_elapsed_dhms (nt, &elap_day, &elap_hour, &elap_min, &elap_sec);
    _ncm_timer_dhms_to_string (nt->msg_tmp1, elap_day, elap_hour, elap_min, elap_sec);
    g_string_printf (nt->msg,
                     "# Task:%s, completed: %u of %u, elapsed time: %s",
                     nt->name, nt->task_pos, nt->task_len, nt->msg_tmp1->str);

    return nt->msg->str;
  }
}

/**
 * ncm_timer_task_mean_time_str:
 * @nt: a #NcmTimer
 *
 * Formats the mean time per item with its standard error. Aborts if no task is
 * running.
 *
 * Returns: (transfer none): the formatted string.
 */
const gchar *
ncm_timer_task_mean_time_str (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  {
    guint day, hour, min;
    gdouble sec;
    gdouble mean_time  = ncm_timer_task_mean_time (nt);
    gdouble sigma_time = ncm_stats_vec_get_sd (nt->time_stats, 0) / sqrt (nt->task_pos);

    _ncm_timer_sec_to_dhms (mean_time, &day, &hour, &min, &sec);
    _ncm_timer_dhms_to_string (nt->msg_tmp1, day, hour, min, sec);

    _ncm_timer_sec_to_dhms (sigma_time, &day, &hour, &min, &sec);
    _ncm_timer_dhms_to_string (nt->msg_tmp2, day, hour, min, sec);

    g_string_printf (nt->msg,
                     "# Task:%s, mean time: %s +/- %s",
                     nt->name, nt->msg_tmp1->str, nt->msg_tmp2->str);

    return nt->msg->str;
  }
}

/**
 * ncm_timer_task_time_left_str:
 * @nt: a #NcmTimer
 *
 * Formats ncm_timer_task_time_left() with its standard error. Aborts if no task is
 * running.
 *
 * Returns: (transfer none): the formatted string.
 */
const gchar *
ncm_timer_task_time_left_str (NcmTimer *nt)
{
  g_assert (nt->task_len != 0);
  {
    guint day, hour, min;
    gdouble sec;
    const gdouble mean_time       = ncm_stats_vec_get_mean (nt->time_stats, 0);
    const gdouble sigma_time      = ncm_stats_vec_get_sd (nt->time_stats, 0) / sqrt (nt->task_pos);
    const guint task_left         = nt->task_len - nt->task_pos;
    const gdouble mean_time_left  = mean_time * task_left;
    const gdouble sigma_time_left = sigma_time * task_left;

    _ncm_timer_sec_to_dhms (mean_time_left, &day, &hour, &min, &sec);
    _ncm_timer_dhms_to_string (nt->msg_tmp1, day, hour, min, sec);

    _ncm_timer_sec_to_dhms (sigma_time_left, &day, &hour, &min, &sec);
    _ncm_timer_dhms_to_string (nt->msg_tmp2, day, hour, min, sec);

    g_string_printf (nt->msg,
                     "# Task:%s, time left: %s +/- %s",
                     nt->name, nt->msg_tmp1->str, nt->msg_tmp2->str);

    return nt->msg->str;
  }
}

/**
 * ncm_timer_task_start_datetime_str:
 * @nt: a #NcmTimer
 *
 * Formats the local date and time now minus ncm_timer_elapsed(), which is the
 * start of the timer if it was never stopped.
 *
 * Returns: (transfer none): the formatted string.
 */
const gchar *
ncm_timer_task_start_datetime_str (NcmTimer *nt)
{
  GDateTime *dt_now   = g_date_time_new_now_local ();
  const gdouble elap  = ncm_timer_elapsed (nt);
  GDateTime *dt_start = g_date_time_add_seconds (dt_now, -elap);
  gchar *start_str    = g_date_time_format (dt_start, "%a %b %d %Y, %T");

  g_string_printf (nt->msg,
                   "# Task:%s, started at: %s",
                   nt->name, start_str);

  g_date_time_unref (dt_now);
  g_date_time_unref (dt_start);
  g_free (start_str);

  return nt->msg->str;
}

/**
 * ncm_timer_task_end_datetime_str:
 * @nt: a #NcmTimer
 *
 * Formats the local date and time at which the task is estimated to end, with the
 * standard error of the time left.
 *
 * Returns: (transfer none): the formatted string.
 */
const gchar *
ncm_timer_task_end_datetime_str (NcmTimer *nt)
{
  GDateTime *dt_now             = g_date_time_new_now_local ();
  const gdouble mean_time       = ncm_stats_vec_get_mean (nt->time_stats, 0);
  const gdouble sigma_time      = ncm_stats_vec_get_sd (nt->time_stats, 0) / sqrt (nt->task_pos);
  const guint task_left         = nt->task_len - nt->task_pos;
  const gdouble mean_time_left  = mean_time * task_left;
  const gdouble sigma_time_left = sigma_time * task_left;
  GDateTime *dt_end             = g_date_time_add_seconds (dt_now, mean_time_left);
  gchar *end_str                = g_date_time_format (dt_end, "%a %b %d %Y, %T");
  guint day, hour, min;
  gdouble sec;

  _ncm_timer_sec_to_dhms (sigma_time_left, &day, &hour, &min, &sec);
  _ncm_timer_dhms_to_string (nt->msg_tmp1, day, hour, min, sec);

  g_string_printf (nt->msg,
                   "# Task:%s, estimated to end at: %s +/- %s",
                   nt->name, end_str, nt->msg_tmp1->str);

  g_date_time_unref (dt_now);
  g_date_time_unref (dt_end);
  g_free (end_str);

  return nt->msg->str;
}

/**
 * ncm_timer_task_cur_datetime_str:
 * @nt: a #NcmTimer
 *
 * Formats the current local date and time.
 *
 * Returns: (transfer none): the formatted string.
 */
const gchar *
ncm_timer_task_cur_datetime_str (NcmTimer *nt)
{
  GDateTime *dt_now = g_date_time_new_now_local ();
  gchar *now_str    = g_date_time_format (dt_now, "%a %b %d %Y, %T");

  g_string_printf (nt->msg,
                   "# Task:%s, current time:        %s",
                   nt->name, now_str);

  g_date_time_unref (dt_now);
  g_free (now_str);

  return nt->msg->str;
}

/**
 * ncm_timer_task_log_elapsed:
 * @nt: a #NcmTimer
 *
 * Logs ncm_timer_task_elapsed_str().
 */
void
ncm_timer_task_log_elapsed (NcmTimer *nt)
{
  g_message ("%s\n", ncm_timer_task_elapsed_str (nt));
  nt->last_log_time = nt->pos_time;
}

/**
 * ncm_timer_task_log_mean_time:
 * @nt: a #NcmTimer
 *
 * Logs ncm_timer_task_mean_time_str().
 */
void
ncm_timer_task_log_mean_time (NcmTimer *nt)
{
  g_message ("%s\n", ncm_timer_task_mean_time_str (nt));
  nt->last_log_time = nt->pos_time;
}

/**
 * ncm_timer_task_log_time_left:
 * @nt: a #NcmTimer
 *
 * Logs ncm_timer_task_time_left_str().
 */
void
ncm_timer_task_log_time_left (NcmTimer *nt)
{
  g_message ("%s\n", ncm_timer_task_time_left_str (nt));
  nt->last_log_time = nt->pos_time;
}

/**
 * ncm_timer_task_log_start_datetime:
 * @nt: a #NcmTimer
 *
 * Logs ncm_timer_task_start_datetime_str().
 */
void
ncm_timer_task_log_start_datetime (NcmTimer *nt)
{
  g_message ("%s\n", ncm_timer_task_start_datetime_str (nt));
  nt->last_log_time = nt->pos_time;
}

/**
 * ncm_timer_task_log_cur_datetime:
 * @nt: a #NcmTimer
 *
 * Logs ncm_timer_task_cur_datetime_str().
 */
void
ncm_timer_task_log_cur_datetime (NcmTimer *nt)
{
  g_message ("%s\n", ncm_timer_task_cur_datetime_str (nt));
  nt->last_log_time = nt->pos_time;
}

/**
 * ncm_timer_task_log_end_datetime:
 * @nt: a #NcmTimer
 *
 * Logs ncm_timer_task_end_datetime_str().
 */
void
ncm_timer_task_log_end_datetime (NcmTimer *nt)
{
  g_message ("%s\n", ncm_timer_task_end_datetime_str (nt));
  nt->last_log_time = nt->pos_time;
}


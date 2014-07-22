/***************************************************************************
 *            ncm_timer.h
 *
 *  Thu August 01 11:23:06 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_timer.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_TIMER_H_
#define _NCM_TIMER_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_stats_vec.h>

G_BEGIN_DECLS

#define NCM_TYPE_TIMER             (ncm_timer_get_type ())
#define NCM_TIMER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_TIMER, NcmTimer))
#define NCM_TIMER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_TIMER, NcmTimerClass))
#define NCM_IS_TIMER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_TIMER))
#define NCM_IS_TIMER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_TIMER))
#define NCM_TIMER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_TIMER, NcmTimerClass))

typedef struct _NcmTimerClass NcmTimerClass;
typedef struct _NcmTimer NcmTimer;

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

struct _NcmTimerClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_timer_get_type (void) G_GNUC_CONST;

NcmTimer *ncm_timer_new (void);
NcmTimer *ncm_timer_ref (NcmTimer *nt);
void ncm_timer_free (NcmTimer *nt);
void ncm_timer_clear (NcmTimer **nt);
void ncm_timer_set_name (NcmTimer *nt, const gchar *name);

gdouble ncm_timer_elapsed (NcmTimer *nt);
void ncm_timer_elapsed_dhms (NcmTimer *nt, guint *elap_day, guint *elap_hour, guint *elap_min, gdouble *elap_sec);
gchar *ncm_timer_elapsed_dhms_str (NcmTimer *nt);
void ncm_timer_start (NcmTimer *nt);
void ncm_timer_stop (NcmTimer *nt);
void ncm_timer_continue (NcmTimer *nt);
void ncm_timer_task_start (NcmTimer *nt, guint task_len);
void ncm_timer_task_pause (NcmTimer *nt);
void ncm_timer_task_continue (NcmTimer *nt);
void ncm_timer_task_add_tasks (NcmTimer *nt, guint ptasks);
void ncm_timer_task_increment (NcmTimer *nt);
void ncm_timer_task_accumulate (NcmTimer *nt, guint nitens);
gboolean ncm_timer_task_is_running (NcmTimer *nt);
gboolean ncm_timer_task_end (NcmTimer *nt);

gdouble ncm_timer_task_mean_time (NcmTimer *nt);
gdouble ncm_timer_task_time_left (NcmTimer *nt);

gchar *ncm_timer_task_elapsed_str (NcmTimer *nt);
gchar *ncm_timer_task_mean_time_str (NcmTimer *nt);
gchar *ncm_timer_task_time_left_str (NcmTimer *nt);
gchar *ncm_timer_task_start_datetime_str (NcmTimer *nt);
gchar *ncm_timer_task_end_datetime_str (NcmTimer *nt);

void ncm_timer_task_log_elapsed (NcmTimer *nt);
void ncm_timer_task_log_mean_time (NcmTimer *nt);
void ncm_timer_task_log_time_left (NcmTimer *nt);
void ncm_timer_task_log_start_datetime (NcmTimer *nt);
void ncm_timer_task_log_end_datetime (NcmTimer *nt);

#define NCM_TIMER_SEC_FORMAT "%07.4f"

G_END_DECLS

#endif /* _NCM_TIMER_H_ */

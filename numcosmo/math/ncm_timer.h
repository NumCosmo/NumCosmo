/***************************************************************************
 *            ncm_timer.h
 *
 *  Thu August 01 11:23:06 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_timer.h
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

#ifndef _NCM_TIMER_H_
#define _NCM_TIMER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NCM_TYPE_TIMER (ncm_timer_get_type ())

G_DECLARE_FINAL_TYPE (NcmTimer, ncm_timer, NCM, TIMER, GObject)

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
guint ncm_timer_task_completed (NcmTimer *nt);
guint ncm_timer_task_estimate_by_time (NcmTimer *nt, gdouble sec);
gboolean ncm_timer_task_is_running (NcmTimer *nt);
gboolean ncm_timer_task_has_ended (NcmTimer *nt);
gboolean ncm_timer_task_end (NcmTimer *nt);

gdouble ncm_timer_elapsed_since_last_log (NcmTimer *nt);
gdouble ncm_timer_task_mean_time (NcmTimer *nt);
gdouble ncm_timer_task_time_left (NcmTimer *nt);

const gchar *ncm_timer_task_elapsed_str (NcmTimer *nt);
const gchar *ncm_timer_task_mean_time_str (NcmTimer *nt);
const gchar *ncm_timer_task_time_left_str (NcmTimer *nt);
const gchar *ncm_timer_task_start_datetime_str (NcmTimer *nt);
const gchar *ncm_timer_task_cur_datetime_str (NcmTimer *nt);
const gchar *ncm_timer_task_end_datetime_str (NcmTimer *nt);

void ncm_timer_task_log_elapsed (NcmTimer *nt);
void ncm_timer_task_log_mean_time (NcmTimer *nt);
void ncm_timer_task_log_time_left (NcmTimer *nt);
void ncm_timer_task_log_start_datetime (NcmTimer *nt);
void ncm_timer_task_log_cur_datetime (NcmTimer *nt);
void ncm_timer_task_log_end_datetime (NcmTimer *nt);

#define NCM_TIMER_SEC_FORMAT "%07.4f"

G_END_DECLS

#endif /* _NCM_TIMER_H_ */


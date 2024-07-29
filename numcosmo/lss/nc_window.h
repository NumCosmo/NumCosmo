/***************************************************************************
 *            nc_window.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_WINDOW_H_
#define _NC_WINDOW_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_WINDOW             (nc_window_get_type ())
#define NC_WINDOW(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_WINDOW, NcWindow))
#define NC_WINDOW_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_WINDOW, NcWindowClass))
#define NC_IS_WINDOW(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_WINDOW))
#define NC_IS_WINDOW_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_WINDOW))
#define NC_WINDOW_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_WINDOW, NcWindowClass))

typedef struct _NcWindowClass NcWindowClass;
typedef struct _NcWindow NcWindow;

struct _NcWindowClass
{
  /*< private > */
  GObjectClass parent_class;
  gdouble volume; /* Volume de uma janela de raio 1 */
  
  gdouble (*eval_fourier) (const NcWindow *wf, const gdouble k, const gdouble R);
  gdouble (*deriv_fourier) (const NcWindow *wf, const gdouble k, const gdouble R);
  gdouble (*eval_real) (const NcWindow *wf, const gdouble r, const gdouble R);
};

struct _NcWindow
{
  /*< private >*/
  GObject parent_instance;
};

GType nc_window_get_type (void) G_GNUC_CONST;

gdouble nc_window_volume (NcWindow *wf);
gdouble nc_window_eval_fourier (const NcWindow *wf, const gdouble k, const gdouble R);
gdouble nc_window_deriv_fourier (const NcWindow *wf, const gdouble k, const gdouble R);
gdouble nc_window_eval_realspace (const NcWindow *wf, const gdouble r, const gdouble R);
void nc_window_free (NcWindow *wf);
void nc_window_clear (NcWindow **wf);

G_END_DECLS

#endif /* _NC_WINDOW_H_ */


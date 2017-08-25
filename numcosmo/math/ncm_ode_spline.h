/***************************************************************************
 *            ncm_ode_spline.h
 *
 *  Wed Aug 13 21:13:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_ODE_SPLINE_H_
#define _NCM_ODE_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_model_ctrl.h>

G_BEGIN_DECLS

#define NCM_TYPE_ODE_SPLINE             (ncm_ode_spline_get_type ())
#define NCM_ODE_SPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_ODE_SPLINE, NcmOdeSpline))
#define NCM_ODE_SPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_ODE_SPLINE, NcmOdeSplineClass))
#define NCM_IS_ODE_SPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_ODE_SPLINE))
#define NCM_IS_ODE_SPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_ODE_SPLINE))
#define NCM_ODE_SPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_ODE_SPLINE, NcmOdeSplineClass))

typedef struct _NcmOdeSplineClass NcmOdeSplineClass;
typedef struct _NcmOdeSpline NcmOdeSpline;

typedef gdouble (*NcmOdeSplineDydx) (gdouble y, gdouble x, gpointer userdata);

struct _NcmOdeSplineClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmOdeSpline
{
  /*< private >*/
  GObject parent_instance;
  gpointer cvode;
  N_Vector y;
  GArray *y_array;
  GArray *x_array;
  gdouble xi;
  gdouble xf;
  gdouble yi;
	gdouble yf;
  gdouble reltol;
  gdouble abstol;
  NcmOdeSplineDydx dydx;
  NcmSpline *s;
  gboolean s_init;
  gboolean cvode_init;
  gboolean hnil;
  gboolean stop_hnil;
  NcmModelCtrl *ctrl;
};

GType ncm_ode_spline_get_type (void) G_GNUC_CONST;

NcmOdeSpline *ncm_ode_spline_new (NcmSpline *s, NcmOdeSplineDydx dydx);
NcmOdeSpline *ncm_ode_spline_new_full (NcmSpline *s, NcmOdeSplineDydx dydx, gdouble yi, gdouble xi, gdouble xf);
void ncm_ode_spline_prepare (NcmOdeSpline *os, gpointer userdata);
void ncm_ode_spline_free (NcmOdeSpline *os);
void ncm_ode_spline_clear (NcmOdeSpline **os);

void ncm_ode_spline_set_interval (NcmOdeSpline *os, gdouble yi, gdouble xi, gdouble xf);
void ncm_ode_spline_set_reltol (NcmOdeSpline *os, gdouble reltol);
void ncm_ode_spline_set_abstol (NcmOdeSpline *os, gdouble abstol);
void ncm_ode_spline_set_xi (NcmOdeSpline *os, gdouble xi);
void ncm_ode_spline_set_xf (NcmOdeSpline *os, gdouble xf);
void ncm_ode_spline_set_yi (NcmOdeSpline *os, gdouble yi);
void ncm_ode_spline_set_yf (NcmOdeSpline *os, gdouble yf);

#define NCM_ODE_SPLINE_DEFAULT_RELTOL (1.0e-13)
#define NCM_ODE_SPLINE_DEFAULT_ABSTOL (1.0e-80)
#define NCM_ODE_SPLINE_MIN_STEP (1.0e-10)

G_END_DECLS

#endif /* _NCM_ODE_SPLINE_H_ */

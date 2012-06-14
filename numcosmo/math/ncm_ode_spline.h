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

#ifndef _NCM_ODE_SPLINE_H
#define _NCM_ODE_SPLINE_H

#include <glib.h>

G_BEGIN_DECLS

typedef gdouble (*NcmOdeSplineDydx) (gdouble y, gdouble x, gpointer userdata);

typedef struct _NcmOdeSpline NcmOdeSpline;

/**
 * NcmOdeSpline:
 *
 * FIXME
 */
struct _NcmOdeSpline
{
  /*< private >*/
  gpointer cvode;
  N_Vector y;
  GArray *y_array;
  GArray *x_array;
  gdouble xi;
  gdouble xf;
  gdouble yi;
  NcmOdeSplineDydx dydx;
  NcmSpline *s;
  gboolean s_init;
  NcmModelCtrl *ctrl;
};

NcmOdeSpline *ncm_ode_spline_new (NcmSpline *s, NcmOdeSplineDydx dydx, gpointer userdata, gdouble yi, gdouble xi, gdouble xf);
void ncm_ode_spline_prepare (NcmOdeSpline *os, gpointer userdata);
void ncm_ode_spline_free (NcmOdeSpline *os);

G_END_DECLS

#endif /* _NCM_ODE_SPLINE_H */

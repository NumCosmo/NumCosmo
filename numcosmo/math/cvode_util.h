/***************************************************************************
 *            cvode_util.h
 *
 *  Wed Nov 12 15:37:44 2008
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

#ifndef _NC_CVODE_UTIL_H
#define _NC_CVODE_UTIL_H

#include <glib.h>
#include <glib-object.h>

#include <gsl/gsl_spline.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#endif

G_BEGIN_DECLS

typedef struct _NcmCVode
{
  gpointer cvode;
  guint n;
  gint spline_index;
  CVRhsFn linear_step;
  CVDlsDenseJacFn linear_J;
  gdouble ti;
  gdouble tf;
  gdouble t;
  N_Vector y0;
  N_Vector y;
  N_Vector abstol;
  gdouble sabstol;
  gdouble reltol;
  gsl_interp_accel *t_accel;
  gsl_spline **splines;
} NcmCVode;

N_Vector ncm_cvode_util_nvector_new (guint n);
gboolean ncm_cvode_util_check_flag (gpointer flagvalue, gchar *funcname, gint opt);
gboolean ncm_cvode_util_print_stats (gpointer cvode);

#define CVODE_CHECK(chk,name,val,ret) \
if (!ncm_cvode_util_check_flag(chk,name,val)) return ret

G_END_DECLS

#endif /* CVODE_UTIL_H */

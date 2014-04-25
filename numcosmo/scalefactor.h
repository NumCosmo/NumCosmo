/***************************************************************************
 *            scalefactor.h
 *
 *  Wed Nov 12 14:46:40 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_SCALEFACTOR_H
#define _NC_SCALEFACTOR_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#endif

G_BEGIN_DECLS

GType nc_scale_factor_get_type (void);

/**
 * NcScaleFactorTimeType:
 * @NC_TIME_COSMIC: Cosmic time
 * @NC_TIME_CONFORMAL: Conformal time
 *
 * FIXME
 */
typedef enum _NcScaleFactorTimeType
{
  NC_TIME_COSMIC = 0,
  NC_TIME_CONFORMAL
} NcScaleFactorTimeType;

typedef struct _NcScaleFactor NcScaleFactor;

/**
 * NcScaleFactor:
 * 
 * FIXME
 */
struct _NcScaleFactor
{
  /*< private >*/
  NcDistance *dist;
  NcScaleFactorTimeType ttype;
  NcmModelCtrl *ctrl;
  NcmSpline *a_t;
  NcmSpline *t_a;
  gdouble a0;
  gdouble zf;
  gdouble ti;
  gdouble tf;
  gboolean spline_init;
  gboolean cvode_malloc;
  gpointer cvode;
  gdouble reltol;
  gdouble abstol;
  N_Vector y;
  CVRhsFn dz_dt_f;
  CVDlsDenseJacFn dz_dt_J;
};

NcScaleFactor *nc_scale_factor_new (NcScaleFactorTimeType ttype, gdouble zf);
NcScaleFactor *nc_scale_factor_copy (NcScaleFactor *a);
void nc_scale_factor_free (NcScaleFactor *a);
void nc_scale_factor_prepare (NcScaleFactor *a, NcHICosmo *model);
void nc_scale_factor_prepare_if_needed (NcScaleFactor *a, NcHICosmo *model);

gdouble nc_scale_factor_z_t (NcScaleFactor *a, gdouble t);
gdouble nc_scale_factor_a_t (NcScaleFactor *a, gdouble t);

gdouble nc_scale_factor_t_z (NcScaleFactor *a, gdouble z);
gdouble nc_scale_factor_t_x (NcScaleFactor *a, gdouble x);

G_END_DECLS

#endif /* SCALEFACTOR_H */

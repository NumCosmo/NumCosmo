/***************************************************************************
 *            ncm_spline_func.h
 *
 *  Wed Aug 13 21:13:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_SPLINE_FUNC_H
#define _NCM_SPLINE_FUNC_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/spline/ncm_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

/**
 * NcmSplineFuncType:
 * @NCM_SPLINE_FUNCTION_SPLINE: adaptive, midpoints in $x$
 * @NCM_SPLINE_FUNCTION_SPLINE_LNKNOT: adaptive, midpoints in $\ln x$, for $x_i > 0$
 * @NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT: adaptive, midpoints in $\sinh^{-1} x$
 * @NCM_SPLINE_FUNC_GRID_LINEAR: fixed grid uniform in $x$
 * @NCM_SPLINE_FUNC_GRID_LOG: fixed grid uniform in $\ln x$, for $x_i > 0$
 *
 * The knot placement methods of #NcmSplineFunc.
 */
typedef enum _NcmSplineFuncType /*< prefix=NCM_SPLINE >*/
{
  NCM_SPLINE_FUNCTION_SPLINE          = 1,
  NCM_SPLINE_FUNCTION_SPLINE_LNKNOT   = 2,
  NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT = 3,
  NCM_SPLINE_FUNC_GRID_LINEAR         = 4,
  NCM_SPLINE_FUNC_GRID_LOG            = 5,
} NcmSplineFuncType;

typedef gdouble (*NcmSplineFuncF) (gdouble x, GObject *obj);

void ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error);
void ncm_spline_set_func_scale (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble scale, gint refine, gdouble refine_ns);
void ncm_spline_set_func1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error);
void ncm_spline_set_func_grid1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize nnodes);
void ncm_spline_set_func_grid (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize nnodes);

/**
 * NCM_SPLINE_FUNC_DEFAULT_MAX_NODES:
 *
 * Upper bound on the number of knots of ncm_spline_set_func_grid().
 */
#define NCM_SPLINE_FUNC_DEFAULT_MAX_NODES 10000000

/**
 * NCM_SPLINE_KNOT_DIFF_TOL:
 *
 * Smallest relative knot spacing of the adaptive #NcmSplineFunc methods.
 */
#define NCM_SPLINE_KNOT_DIFF_TOL (GSL_DBL_EPSILON * 1.0e2)

G_END_DECLS

#endif /* _NCM_SPLINE_FUNC_H */


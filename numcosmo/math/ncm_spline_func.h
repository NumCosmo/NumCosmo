/***************************************************************************
 *            ncm_spline_func.h
 *
 *  Wed Aug 13 21:13:59 2008
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

#ifndef _NCM_SPLINE_FUNC_H
#define _NCM_SPLINE_FUNC_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

/**
 * NcmSplineFuncType:
 * @NCM_SPLINE_FUNCTION_4POINTS: FIXME 
 * @NCM_SPLINE_FUNCTION_SPLINE: The knots are evenly distributed on a linear base at each step. The test points are place at $\overline{\mathbf{x}} = \frac{\mathbf{x}^{i+1} + \mathbf{x}^{i}}{2}$. 
 * @NCM_SPLINE_FUNCTION_SPLINE_LNKNOT: The knots are evenly distributed on a logarithm base at each step. The test points are place at $\overline{\mathbf{x}} = \mathrm{exp}\left( \frac{\ln \mathbf{x}^{i+1} + \ln \mathbf{x}^{i}}{2}  \right)$. This method is only applied for positive intervals and is indicated for functions that changes orders of magnitude across the interval.
 * @NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT: The knots are evenly distributed on a hyperbolic sine base at each step. The test points are place at $\overline{\mathbf{x}} = \sinh \left[ \frac{\sinh^{-1} \left( \mathbf{x}^{i+1} \right) + \sinh^{-1} \left( \mathbf{x}^{i} \right)}{2} \right]$. This method is indicated for functions that changes orders of magnitude across the interval.
 * @NCM_SPLINE_FUNC_GRID_LINEAR: The knots are evenly distributed on a linear base in the entire range [@xi, @xf].
 * @NCM_SPLINE_FUNC_GRID_LOG: The knots are evenly distributed on a natural logarithmic base in the entire range [@xi > 0, @xf > xi > 0]. 
 *
 * Enumeration to choose which of the functions to be applied when interpolating the input #gsl_function *@F, $f$, 
 * with the desired @rel_error in the range [@xi, @xf]. 
 * The interpolation knots, $\mathbf{x}$, are automatically defined internally by the functions. 
 * For more details see [description][numcosmo-NcmSplineFunc.description] above.
 * 
 */
typedef enum _NcmSplineFuncType
{
  NCM_SPLINE_FUNCTION_4POINTS,
  NCM_SPLINE_FUNCTION_SPLINE,
  NCM_SPLINE_FUNCTION_SPLINE_LNKNOT,
  NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT,
  NCM_SPLINE_FUNC_GRID_LINEAR,
  NCM_SPLINE_FUNC_GRID_LOG,
} NcmSplineFuncType;

typedef gdouble (*NcmSplineFuncF) (gdouble x, GObject *obj);

void ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error);
void ncm_spline_set_func_scale (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble scale, gboolean refine);
void ncm_spline_set_func1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error);
void ncm_spline_set_func_grid1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize nnodes);
void ncm_spline_set_func_grid (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize nnodes);

#define NCM_SPLINE_FUNC_DEFAULT_MAX_NODES 10000000
#define NCM_SPLINE_KNOT_DIFF_TOL (GSL_DBL_EPSILON * 1.0e2)

G_END_DECLS

#endif /* _NCM_SPLINE_FUNC_H */


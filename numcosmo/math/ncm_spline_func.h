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
 * @NCM_SPLINE_FUNCTION_2x2POINTS: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE_LNKNOT: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT: FIXME
 * 
 * FIXME
 */ 
typedef enum _NcmSplineFuncType
{
  NCM_SPLINE_FUNCTION_4POINTS,
  NCM_SPLINE_FUNCTION_2x2POINTS,
  NCM_SPLINE_FUNCTION_SPLINE,
  NCM_SPLINE_FUNCTION_SPLINE_LNKNOT,
  NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT,
} NcmSplineFuncType;

void ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error);

#define NCM_SPLINE_FUNC_DEFAULT_MAX_NODES 10000000
#define NCM_SPLINE_KNOT_DIFF_TOL (GSL_DBL_EPSILON * 1.0e2)

G_END_DECLS

#endif /* _NCM_SPLINE_FUNC_H */

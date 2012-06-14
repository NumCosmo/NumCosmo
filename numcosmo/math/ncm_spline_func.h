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

G_BEGIN_DECLS

/**
 * NcmSplineFuncType:
 * @NCM_SPLINE_FUNCTION_4POINTS: FIXME
 * @NCM_SPLINE_FUNCTION_2x2POINTS: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE_LNKNOT: FIXME
 * 
 * FIXME
 */ 
typedef enum _NcmSplineFuncType
{
  NCM_SPLINE_FUNCTION_4POINTS,
  NCM_SPLINE_FUNCTION_2x2POINTS,
  NCM_SPLINE_FUNCTION_SPLINE,
  NCM_SPLINE_FUNCTION_SPLINE_LNKNOT,
} NcmSplineFuncType;

void ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error);

G_END_DECLS

#endif /* _NCM_SPLINE_FUNC_H */

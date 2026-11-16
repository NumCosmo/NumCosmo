/***************************************************************************
 *            ncm_spline_bspline_private.h
 *
 *  Sun Nov 15 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * ncm_spline_bspline_private.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NCM_SPLINE_BSPLINE_PRIVATE_H_
#define _NCM_SPLINE_BSPLINE_PRIVATE_H_

#include <glib.h>
#include "ncm/spline/ncm_spline_bspline.h"

G_BEGIN_DECLS

/*
 * The internal layout of a prepared B-spline, for code of the library that copies it and
 * evaluates it without the object (see ncm_spline_bspline.c). Not part of the public API.
 */
void ncm_spline_bspline_peek_basis (NcmSplineBSpline *sbs, gsize *n_control, const gdouble **knots, const gdouble **c);

gdouble ncm_spline_bspline_eval_basis (const guint order, const gsize n_control, const gdouble *knots, const gdouble *c, const gdouble x);

G_END_DECLS

#endif /* _NCM_SPLINE_BSPLINE_PRIVATE_H_ */


/***************************************************************************
 *            ncm_spline_bspline.h
 *
 *  Wed Aug 20 10:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _NCM_SPLINE_BSPLINE_H_
#define _NCM_SPLINE_BSPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/spline/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_BSPLINE (ncm_spline_bspline_get_type ())

G_DECLARE_FINAL_TYPE (NcmSplineBSpline, ncm_spline_bspline, NCM, SPLINE_BSPLINE, NcmSpline)

NcmSplineBSpline *ncm_spline_bspline_new (guint order);
NcmSplineBSpline *ncm_spline_bspline_new_full (guint order, NcmVector *xv, NcmVector *yv, gboolean init);
NcmSplineBSpline *ncm_spline_bspline_new_tol (gdouble reltol, gdouble abstol);

void ncm_spline_bspline_set_order (NcmSplineBSpline *sbs, guint order);
guint ncm_spline_bspline_get_order (NcmSplineBSpline *sbs);
gdouble ncm_spline_bspline_get_achieved_error (NcmSplineBSpline *sbs);

/**
 * NCM_SPLINE_BSPLINE_DEFAULT_ORDER:
 *
 * Default B-spline order (degree 7). Degree 3 cannot reach machine precision at any
 * sample density; degree 9 and above gain nothing and lose conditioning.
 */
#define NCM_SPLINE_BSPLINE_DEFAULT_ORDER (8)

/**
 * NCM_SPLINE_BSPLINE_MAX_ORDER:
 *
 * Largest supported order. Above this the interpolation matrix conditioning costs more
 * than the extra order gains.
 */
#define NCM_SPLINE_BSPLINE_MAX_ORDER (10)

G_END_DECLS

#endif /* _NCM_SPLINE_BSPLINE_H_ */


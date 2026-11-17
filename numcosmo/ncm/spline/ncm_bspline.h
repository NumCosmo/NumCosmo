/***************************************************************************
 *            ncm_bspline.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_bspline.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NCM_BSPLINE_H_
#define _NCM_BSPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_BSPLINE (ncm_bspline_get_type ())

G_DECLARE_FINAL_TYPE (NcmBSpline, ncm_bspline, NCM, BSPLINE, GObject)

NcmBSpline *ncm_bspline_new (const guint order, NcmVector * knots, NcmVector * control_points);
NcmBSpline *ncm_bspline_ref (NcmBSpline *bs);
void ncm_bspline_free (NcmBSpline *bs);
void ncm_bspline_clear (NcmBSpline **bs);

guint ncm_bspline_get_order (NcmBSpline *bs);
NcmVector *ncm_bspline_peek_knots (NcmBSpline *bs);
NcmVector *ncm_bspline_peek_control_points (NcmBSpline *bs);
void ncm_bspline_get_range (NcmBSpline *bs, gdouble *x_lo, gdouble *x_hi);

gdouble ncm_bspline_eval (const NcmBSpline *bs, const gdouble x);

G_END_DECLS

#endif /* _NCM_BSPLINE_H_ */


/***************************************************************************
 *            ncm_spline_vec.h
 *
 *  Sat Mar 15 19:53:22 2026
 *  Copyright  2026
 *  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_SPLINE_VEC_H_
#define _NCM_SPLINE_VEC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_VEC (ncm_spline_vec_get_type ())

G_DECLARE_FINAL_TYPE (NcmSplineVec, ncm_spline_vec, NCM, SPLINE_VEC, GObject)

NcmSplineVec *ncm_spline_vec_new (const NcmSpline * s, NcmVector * xv, NcmMatrix * ym, const gboolean init);
NcmSplineVec *ncm_spline_vec_new_gpa (const NcmSpline *s, NcmVector *xv, GPtrArray *yv, const gboolean init);
NcmSplineVec *ncm_spline_vec_ref (NcmSplineVec *sv);

void ncm_spline_vec_set (NcmSplineVec *sv, NcmVector *xv, NcmMatrix *ym, gboolean init);
void ncm_spline_vec_set_gpa (NcmSplineVec *sv, NcmVector *xv, GPtrArray *yv, gboolean init);

void ncm_spline_vec_free (NcmSplineVec *sv);
void ncm_spline_vec_clear (NcmSplineVec **sv);

void ncm_spline_vec_prepare (NcmSplineVec *sv);
gboolean ncm_spline_vec_is_init (NcmSplineVec *sv);

guint ncm_spline_vec_get_len (NcmSplineVec *sv);
NcmSpline *ncm_spline_vec_peek_spline (NcmSplineVec *sv, guint i);

void ncm_spline_vec_eval (NcmSplineVec *sv, const gdouble x, NcmVector *res);
void ncm_spline_vec_deriv (NcmSplineVec *sv, const gdouble x, NcmVector *res);
void ncm_spline_vec_integ (NcmSplineVec *sv, const gdouble xi, const gdouble xf, NcmVector *res);

void ncm_spline_vec_eval_array (NcmSplineVec *sv, const gdouble x, GArray **res);
void ncm_spline_vec_deriv_array (NcmSplineVec *sv, const gdouble x, GArray **res);
void ncm_spline_vec_integ_array (NcmSplineVec *sv, const gdouble xi, const gdouble xf, GArray **res);

G_END_DECLS

#endif /* _NCM_SPLINE_VEC_H_ */


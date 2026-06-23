/***************************************************************************
 *            ncm_spline.h
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

#ifndef _NCM_SPLINE_H_
#define _NCM_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_spline.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE (ncm_spline_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmSpline, ncm_spline, NCM, SPLINE, GObject)

/**
 * NcmSplineCurvatureType:
 * @NCM_SPLINE_CURVATURE_D2: plain second derivative, $c(x) = f''(x)$
 * @NCM_SPLINE_CURVATURE_GEOMETRIC: geometric curvature of the curve $(x, f(x))$,
 *   $c(x) = f''(x) / \left(1 + f'(x)^2\right)^{3/2}$
 *
 * Choice of the curvature density $c(x)$ used by the $L_p$ curvature functionals
 * ncm_spline_curvature_lp_norm() and ncm_spline_curvature_max().
 *
 */
typedef enum _NcmSplineCurvatureType
{
  NCM_SPLINE_CURVATURE_D2,
  NCM_SPLINE_CURVATURE_GEOMETRIC,
} NcmSplineCurvatureType;

struct _NcmSplineClass
{
  /*< private >*/
  GObjectClass parent_class;
  const gchar *(*name) (NcmSpline *s);

  void (*reset) (NcmSpline *s);
  void (*prepare) (NcmSpline *s);
  void (*prepare_base) (NcmSpline *s);
  gsize (*min_size) (const NcmSpline *s);
  gdouble (*eval) (const NcmSpline *s, const gdouble x);
  gdouble (*deriv) (const NcmSpline *s, const gdouble x);
  gdouble (*deriv2) (const NcmSpline *s, const gdouble x);
  gdouble (*deriv_nmax) (const NcmSpline *s, const gdouble x);
  gdouble (*integ) (const NcmSpline *s, const gdouble xi, const gdouble xf);
  NcmSpline *(*copy_empty) (const NcmSpline *s);
  gdouble (*eval_idx) (const NcmSpline *s, const gdouble x, const gsize i);
  gdouble (*deriv_idx) (const NcmSpline *s, const gdouble x, const gsize i);
  gdouble (*integ_idx) (const NcmSpline *s, const gdouble xi, const gsize i, const gdouble xf, const gsize f);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[4];
};

NcmSpline *ncm_spline_copy_empty (const NcmSpline *s);
NcmSpline *ncm_spline_copy (const NcmSpline *s);

NcmSpline *ncm_spline_new (const NcmSpline *s, NcmVector *xv, NcmVector *yv, const gboolean init);
NcmSpline *ncm_spline_new_array (const NcmSpline *s, GArray *x, GArray *y, const gboolean init);
NcmSpline *ncm_spline_new_data (const NcmSpline *s, gdouble *x, gdouble *y, const gsize len, const gboolean init);
NcmSpline *ncm_spline_set (NcmSpline *s, NcmVector *xv, NcmVector *yv, gboolean init);
NcmSpline *ncm_spline_ref (NcmSpline *s);

void ncm_spline_acc (NcmSpline *s, gboolean enable);
gsl_interp_accel *ncm_spline_peek_acc (NcmSpline *s);

void ncm_spline_set_len (NcmSpline *s, guint len);
void ncm_spline_set_xv (NcmSpline *s, NcmVector *xv, gboolean init);
void ncm_spline_set_yv (NcmSpline *s, NcmVector *yv, gboolean init);
void ncm_spline_set_array (NcmSpline *s, GArray *x, GArray *y, gboolean init);
void ncm_spline_set_data_static (NcmSpline *s, gdouble *x, gdouble *y, gsize len, gboolean init);

guint ncm_spline_get_len (NcmSpline *s);
NcmVector *ncm_spline_get_xv (NcmSpline *s);
NcmVector *ncm_spline_get_yv (NcmSpline *s);
NcmVector *ncm_spline_peek_xv (NcmSpline *s);
NcmVector *ncm_spline_peek_yv (NcmSpline *s);

void ncm_spline_get_bounds (NcmSpline *s, gdouble *lb, gdouble *ub);
gboolean ncm_spline_is_init (NcmSpline *s);

void ncm_spline_free (NcmSpline *s);
void ncm_spline_clear (NcmSpline **s);

void ncm_spline_prepare (NcmSpline *s);
void ncm_spline_post_prepare (NcmSpline *s);
void ncm_spline_prepare_base (NcmSpline *s);
gdouble ncm_spline_eval (const NcmSpline *s, const gdouble x);
gdouble ncm_spline_eval_idx (const NcmSpline *s, const gdouble x, const gsize i);
gdouble ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x);
gdouble ncm_spline_eval_deriv_idx (const NcmSpline *s, const gdouble x, const gsize i);
gdouble ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x);
gdouble ncm_spline_eval_deriv_nmax (const NcmSpline *s, const gdouble x);
gdouble ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
gdouble ncm_spline_eval_integ_idx (const NcmSpline *s, const gdouble xi, const gsize i, const gdouble xf, const gsize f);
gboolean ncm_spline_is_empty (const NcmSpline *s);
gsize ncm_spline_min_size (const NcmSpline *s);
guint ncm_spline_get_index (const NcmSpline *s, const gdouble x);

gdouble ncm_spline_curvature_density (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble x);
gdouble ncm_spline_curvature_lp_norm (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble p, const gdouble xi, const gdouble xf);
gdouble ncm_spline_curvature_max (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble xi, const gdouble xf);

/* Utilities -- internal use */

gdouble _ncm_spline_util_integ_eval (const gdouble ai, const gdouble bi, const gdouble ci, const gdouble di, const gdouble xi, const gdouble a, const gdouble b);

G_END_DECLS

#endif /* _NCM_SPLINE_H_ */


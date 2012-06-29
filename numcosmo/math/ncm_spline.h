/***************************************************************************
 *            ncm_spline.h
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

#ifndef _NCM_SPLINE_H_
#define _NCM_SPLINE_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE             (ncm_spline_get_type ())
#define NCM_SPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE, NcmSpline))
#define NCM_SPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE, NcmSplineClass))
#define NCM_IS_SPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE))
#define NCM_IS_SPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE))
#define NCM_SPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE, NcmSplineClass))

typedef struct _NcmSplineClass NcmSplineClass;
typedef struct _NcmSpline NcmSpline;

struct _NcmSplineClass
{
  /*< private >*/
	GObjectClass parent_class;
	gchar *name;
	void (*reset) (NcmSpline *s);
	void (*prepare) (NcmSpline *s);
	void (*prepare_base) (NcmSpline *s);
	gsize (*min_size) (const NcmSpline *s);
	gdouble (*eval) (const NcmSpline *s, const gdouble x);
	gdouble (*deriv) (const NcmSpline *s, const gdouble x);
	gdouble (*deriv2) (const NcmSpline *s, const gdouble x);
	gdouble (*integ) (const NcmSpline *s, const gdouble xi, const gdouble xf);
	NcmSpline *(*copy_empty) (const NcmSpline *s);
};

struct _NcmSpline
{
  /*< private >*/
  GObject parent_instance;
  gsize len;
  NcmVector *xv;
  NcmVector *yv;
  gsl_interp_accel *acc;
  gboolean init;
  gboolean empty;
};

GType ncm_spline_get_type (void) G_GNUC_CONST;

NcmSpline *ncm_spline_copy_empty (const NcmSpline *s);
NcmSpline *ncm_spline_copy (const NcmSpline *s);

NcmSpline *ncm_spline_new (const NcmSpline *s, NcmVector *xv, NcmVector *yv, const gboolean init);
NcmSpline *ncm_spline_new_array (const NcmSpline *s, GArray *x, GArray *y, const gboolean init);
NcmSpline *ncm_spline_new_data (const NcmSpline *s, gdouble *x, gdouble *y, const gsize len, const gboolean init);
NcmSpline *ncm_spline_set (NcmSpline *s, NcmVector *xv, NcmVector *yv, gboolean init);
NcmSpline *ncm_spline_ref (NcmSpline *s);

void ncm_spline_set_xv (NcmSpline *s, NcmVector *xv, gboolean init);
void ncm_spline_set_yv (NcmSpline *s, NcmVector *yv, gboolean init);
void ncm_spline_set_array (NcmSpline *s, GArray *x, GArray *y, gboolean init);
void ncm_spline_set_data_static (NcmSpline *s, gdouble *x, gdouble *y, gsize len, gboolean init);

NcmVector *ncm_spline_get_xv (NcmSpline *s);
NcmVector *ncm_spline_get_yv (NcmSpline *s);

void ncm_spline_free (NcmSpline *s);

G_INLINE_FUNC void ncm_spline_prepare (NcmSpline *s);
G_INLINE_FUNC void ncm_spline_prepare_base (NcmSpline *s);
G_INLINE_FUNC gdouble ncm_spline_eval (const NcmSpline *s, const gdouble x);
G_INLINE_FUNC gdouble ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x);
G_INLINE_FUNC gdouble ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x);
G_INLINE_FUNC gdouble ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
G_INLINE_FUNC gboolean ncm_spline_is_empty (const NcmSpline *s);
G_INLINE_FUNC gsize ncm_spline_min_size (const NcmSpline *s);
G_INLINE_FUNC guint ncm_spline_get_index (const NcmSpline *s, const gdouble x);

/* Utilities -- internal use */

G_INLINE_FUNC gdouble _ncm_spline_util_integ_eval (const gdouble ai, const gdouble bi, const gdouble ci, const gdouble di, const gdouble xi, const gdouble a, const gdouble b);

G_END_DECLS

#endif /* _NCM_SPLINE_H_ */

#ifndef _NCM_SPLINE_INLINE_H_
#define _NCM_SPLINE_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC void
ncm_spline_prepare (NcmSpline *s)
{
  s->init = TRUE;
  NCM_SPLINE_GET_CLASS (s)->prepare (s);
}

G_INLINE_FUNC void
ncm_spline_prepare_base (NcmSpline *s)
{
	if (NCM_SPLINE_GET_CLASS (s)->prepare_base)
		NCM_SPLINE_GET_CLASS (s)->prepare_base (s);
}

G_INLINE_FUNC gdouble
ncm_spline_eval (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->eval (s, x);
}

G_INLINE_FUNC gdouble
ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->deriv (s, x);
}

G_INLINE_FUNC gdouble
ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->deriv2 (s, x);
}

G_INLINE_FUNC gdouble
ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  return NCM_SPLINE_GET_CLASS (s)->integ (s, x0, x1);
}

G_INLINE_FUNC gboolean
ncm_spline_is_empty (const NcmSpline *s)
{
	return s->empty;
}

G_INLINE_FUNC gsize
ncm_spline_min_size (const NcmSpline *s)
{
	return NCM_SPLINE_GET_CLASS (s)->min_size (s);
}

G_INLINE_FUNC guint
ncm_spline_get_index (const NcmSpline *s, const gdouble x)
{
  if (s->acc && ncm_vector_stride (s->xv) == 1)
  	return gsl_interp_accel_find (s->acc, ncm_vector_ptr (s->xv, 0), s->len, x);
	else
		return gsl_interp_bsearch (ncm_vector_ptr (s->xv, 0), x, 0, ncm_vector_len (s->xv) - 1);
}

/* Utilities -- internal use */

G_INLINE_FUNC gdouble
_ncm_spline_util_integ_eval (const gdouble ai, const gdouble bi, const gdouble ci, const gdouble di, const gdouble xi, const gdouble a, const gdouble b)
{
  const gdouble r1 = a - xi;
  const gdouble r2 = b - xi;
  const gdouble r12 = r1 + r2;
  const gdouble bterm = 0.5 * bi * r12;
  const gdouble cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
  const gdouble dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);

  return (b - a) * (ai + bterm + cterm + dterm);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SPLINE_INLINE_H_ */

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

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_spline.h>
#endif /* NUMCOSMO_GIR_SCAN */

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

void ncm_spline_acc (NcmSpline *s, gboolean enable);

void ncm_spline_set_len (NcmSpline *s, guint len);
void ncm_spline_set_xv (NcmSpline *s, NcmVector *xv, gboolean init);
void ncm_spline_set_yv (NcmSpline *s, NcmVector *yv, gboolean init);
void ncm_spline_set_array (NcmSpline *s, GArray *x, GArray *y, gboolean init);
void ncm_spline_set_data_static (NcmSpline *s, gdouble *x, gdouble *y, gsize len, gboolean init);

guint ncm_spline_get_len (NcmSpline *s);
NcmVector *ncm_spline_get_xv (NcmSpline *s);
NcmVector *ncm_spline_get_yv (NcmSpline *s);
void ncm_spline_get_bounds (NcmSpline *s, gdouble *lb, gdouble *ub);

void ncm_spline_free (NcmSpline *s);
void ncm_spline_clear (NcmSpline **s);

NCM_INLINE void ncm_spline_prepare (NcmSpline *s);
NCM_INLINE void ncm_spline_prepare_base (NcmSpline *s);
NCM_INLINE gdouble ncm_spline_eval (const NcmSpline *s, const gdouble x);
NCM_INLINE gdouble ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x);
NCM_INLINE gdouble ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x);
NCM_INLINE gdouble ncm_spline_eval_deriv_nmax (const NcmSpline *s, const gdouble x);
NCM_INLINE gdouble ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
NCM_INLINE gboolean ncm_spline_is_empty (const NcmSpline *s);
NCM_INLINE gsize ncm_spline_min_size (const NcmSpline *s);
NCM_INLINE guint ncm_spline_get_index (const NcmSpline *s, const gdouble x);

/* Utilities -- internal use */

NCM_INLINE gdouble _ncm_spline_util_integ_eval (const gdouble ai, const gdouble bi, const gdouble ci, const gdouble di, const gdouble xi, const gdouble a, const gdouble b);

G_END_DECLS

#endif /* _NCM_SPLINE_H_ */

#ifndef _NCM_SPLINE_INLINE_H_
#define _NCM_SPLINE_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
ncm_spline_prepare (NcmSpline *s)
{
  s->init = TRUE;
  NCM_SPLINE_GET_CLASS (s)->prepare (s);
}

NCM_INLINE void
ncm_spline_prepare_base (NcmSpline *s)
{
	if (NCM_SPLINE_GET_CLASS (s)->prepare_base)
		NCM_SPLINE_GET_CLASS (s)->prepare_base (s);
}

NCM_INLINE gdouble
ncm_spline_eval (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->eval (s, x);
}

NCM_INLINE gdouble
ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->deriv (s, x);
}

NCM_INLINE gdouble
ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->deriv2 (s, x);
}

NCM_INLINE gdouble
ncm_spline_eval_deriv_nmax (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS (s)->deriv_nmax (s, x);
}

NCM_INLINE gdouble
ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  return NCM_SPLINE_GET_CLASS (s)->integ (s, x0, x1);
}

NCM_INLINE gboolean
ncm_spline_is_empty (const NcmSpline *s)
{
	return s->empty;
}

NCM_INLINE gsize
ncm_spline_min_size (const NcmSpline *s)
{
	return NCM_SPLINE_GET_CLASS (s)->min_size (s);
}

NCM_INLINE gsize
_ncm_spline_bsearch_stride (const gdouble x_array[], const guint stride, const gdouble x, gsize index_lo, gsize index_hi)
{
  gsize ilo = index_lo;
  gsize ihi = index_hi;
	
  while (ihi > ilo + 1) 
	{
    gsize i = (ihi + ilo)/2;
    if (x_array[i * stride] > x)
      ihi = i;
    else
      ilo = i;
  }
  
  return ilo;
}

NCM_INLINE gsize
_ncm_spline_accel_find (gsl_interp_accel *a, const gdouble xa[], const guint stride, gsize len, gdouble x)
{
  gsize x_index = a->cache;
 
  if (x < xa[x_index * stride]) 
	{
    a->miss_count++;
    a->cache = _ncm_spline_bsearch_stride (xa, stride, x, 0, x_index);
  }
  else if (x >= xa[stride * (x_index + 1)]) 
	{
    a->miss_count++;
    a->cache = _ncm_spline_bsearch_stride (xa, stride, x, x_index, len - 1);
  }
  else 
	{
    a->hit_count++;
  }
  
  return a->cache;
}

NCM_INLINE guint
ncm_spline_get_index (const NcmSpline *s, const gdouble x)
{
	if (ncm_vector_stride (s->xv) == 1)
	{
		if (s->acc)
			return gsl_interp_accel_find (s->acc, ncm_vector_ptr (s->xv, 0), s->len, x);
		else
			return gsl_interp_bsearch (ncm_vector_ptr (s->xv, 0), x, 0, ncm_vector_len (s->xv) - 1);
	}
	else
	{
		if (s->acc)
			return _ncm_spline_accel_find (s->acc, ncm_vector_ptr (s->xv, 0), ncm_vector_stride (s->xv), s->len, x);
		else
			return _ncm_spline_bsearch_stride (ncm_vector_ptr (s->xv, 0), ncm_vector_stride (s->xv), x, 0, ncm_vector_len (s->xv) - 1);
	}
}

/* Utilities -- internal use */

NCM_INLINE gdouble
_ncm_spline_util_integ_eval (const gdouble ai, const gdouble bi, const gdouble ci, const gdouble di, const gdouble xi, const gdouble a, const gdouble b)
{
  const gdouble r1    = (a - xi);
  const gdouble r2    = (b - xi);
  const gdouble r12   = (r1 + r2);
  const gdouble bterm = 0.5 * bi * r12;
  const gdouble cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
  const gdouble dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);

  return (b - a) * (ai + bterm + cterm + dterm);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SPLINE_INLINE_H_ */

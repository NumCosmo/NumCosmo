/***************************************************************************
 *            ncm_spline2d.h
 *
 *  Sun Aug  1 17:17:20 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <sandro@isoftware.com.br>
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

#ifndef _NCM_SPLINE2D_H_
#define _NCM_SPLINE2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_spline_func.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE2D             (ncm_spline2d_get_type ())
#define NCM_SPLINE2D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE2D, NcmSpline2d))
#define NCM_SPLINE2D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE2D, NcmSpline2dClass))
#define NCM_IS_SPLINE2D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE2D))
#define NCM_IS_SPLINE2D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE2D))
#define NCM_SPLINE2D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE2D, NcmSpline2dClass))

typedef struct _NcmSpline2dClass NcmSpline2dClass;
typedef struct _NcmSpline2d NcmSpline2d;

struct _NcmSpline2dClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmSpline2d *(*copy_empty) (const NcmSpline2d *s2d);
  void (*reset) (NcmSpline2d *s2d);
  void (*prepare) (NcmSpline2d *s2d);
  gdouble (*eval) (NcmSpline2d *s2d, gdouble x, gdouble y);
  gdouble (*dzdx) (NcmSpline2d *s2d, gdouble x, gdouble y);
  gdouble (*dzdy) (NcmSpline2d *s2d, gdouble x, gdouble y);
  gdouble (*d2zdxy) (NcmSpline2d *s2d, gdouble x, gdouble y);
  gdouble (*d2zdx2) (NcmSpline2d *s2d, gdouble x, gdouble y);
  gdouble (*d2zdy2) (NcmSpline2d *s2d, gdouble x, gdouble y);
  gdouble (*int_dx) (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
  gdouble (*int_dy) (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
  gdouble (*int_dxdy) (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
  NcmSpline *(*int_dx_spline) (NcmSpline2d *s2d, gdouble xl, gdouble xu);
  NcmSpline *(*int_dy_spline) (NcmSpline2d *s2d, gdouble yl, gdouble yu);
};

struct _NcmSpline2d
{
  /*< private >*/
  GObject parent_instance;
  gboolean empty;
  gboolean init;
  gboolean to_init;
  NcmSpline *s;
  NcmVector *xv;
  NcmVector *yv;
  NcmMatrix *zm;
  gsl_interp_accel *acc_x;
  gsl_interp_accel *acc_y;
  gboolean use_acc;
  gboolean no_stride;
};

GType ncm_spline2d_get_type (void) G_GNUC_CONST;

void ncm_spline2d_set (NcmSpline2d *s2d, NcmVector *xv, NcmVector *yv, NcmMatrix *zm, gboolean init);
void ncm_spline2d_set_function (NcmSpline2d *s2d, NcmSplineFuncType ftype, gsl_function *Fx, gsl_function *Fy, gdouble xl, gdouble xu, gdouble yl, gdouble yu, gdouble rel_err);
void ncm_spline2d_prepare (NcmSpline2d *s2d);
guint ncm_spline2d_min_size (NcmSpline2d *s2d);

NcmSpline2d *ncm_spline2d_copy_empty (const NcmSpline2d *s2d);
NcmSpline2d *ncm_spline2d_copy (NcmSpline2d *s2d);

NcmSpline2d *ncm_spline2d_new (const NcmSpline2d *s2d, NcmVector *xv, NcmVector *yv, NcmMatrix *zm, gboolean init);

void ncm_spline2d_free (NcmSpline2d *s2d);
void ncm_spline2d_clear (NcmSpline2d **s2d);

void ncm_spline2d_use_acc (NcmSpline2d *s2d, gboolean use_acc);

NCM_INLINE gdouble ncm_spline2d_eval (NcmSpline2d *s2d, gdouble x, gdouble y);
gdouble ncm_spline2d_integ_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
gdouble ncm_spline2d_integ_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
NcmSpline *ncm_spline2d_integ_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu);
NcmSpline *ncm_spline2d_integ_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dx_spline_val (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
gdouble ncm_spline2d_integ_dy_spline_val (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dxdy_spline_x (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dxdy_spline_y (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);

NCM_INLINE gdouble ncm_spline2d_deriv_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y);
NCM_INLINE gdouble ncm_spline2d_deriv_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y);
NCM_INLINE gdouble ncm_spline2d_deriv_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y);
NCM_INLINE gdouble ncm_spline2d_deriv_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y);
NCM_INLINE gdouble ncm_spline2d_deriv_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y);

NCM_INLINE gdouble ncm_spline2dim_integ_total (NcmSpline2d *s2d);

G_END_DECLS

#endif /* _NCM_SPLINE2D_H_ */

#ifndef _NCM_SPLINE2D_INLINE_H_
#define _NCM_SPLINE2D_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
ncm_spline2d_eval (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  return NCM_SPLINE2D_GET_CLASS (s2d)->eval (s2d, x, y);
}

NCM_INLINE gdouble
ncm_spline2dim_integ_total (NcmSpline2d *s2d)
{
	return ncm_spline2d_integ_dxdy (s2d,
	                                ncm_vector_get (s2d->xv, 0),
	                                ncm_vector_get (s2d->xv, ncm_vector_len (s2d->xv) - 1),
	                                ncm_vector_get (s2d->yv, 0),
	                                ncm_vector_get (s2d->yv, ncm_vector_len (s2d->yv) - 1)
	                                );
}

NCM_INLINE gdouble 
ncm_spline2d_deriv_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->dzdx (s2d, x, y);
}

NCM_INLINE gdouble 
ncm_spline2d_deriv_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->dzdy (s2d, x, y);
}

NCM_INLINE gdouble 
ncm_spline2d_deriv_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->d2zdxy (s2d, x, y);
}

NCM_INLINE gdouble 
ncm_spline2d_deriv_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->d2zdx2 (s2d, x, y);
}

NCM_INLINE gdouble 
ncm_spline2d_deriv_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->d2zdy2 (s2d, x, y);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SPLINE2D_INLINE_H_ */

/***************************************************************************
 *            ncm_spline2d.h
 *
 *  Sun Aug  1 17:17:20 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <vitenti@uel.br>
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

#define NCM_TYPE_SPLINE2D (ncm_spline2d_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmSpline2d, ncm_spline2d, NCM, SPLINE2D, GObject)

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
  void (*eval_vec_y) (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[3];
};

void ncm_spline2d_set (NcmSpline2d *s2d, NcmVector *xv, NcmVector *yv, NcmMatrix *zm, gboolean init);
void ncm_spline2d_set_function (NcmSpline2d *s2d, NcmSplineFuncType ftype, gsl_function *Fx, gsl_function *Fy, gdouble xl, gdouble xu, gdouble yl, gdouble yu, gdouble rel_err);
void ncm_spline2d_prepare (NcmSpline2d *s2d);
guint ncm_spline2d_min_size (NcmSpline2d *s2d);

NcmSpline2d *ncm_spline2d_copy_empty (const NcmSpline2d *s2d);
NcmSpline2d *ncm_spline2d_copy (NcmSpline2d *s2d);

NcmSpline2d *ncm_spline2d_new (const NcmSpline2d *s2d, NcmVector *xv, NcmVector *yv, NcmMatrix *zm, gboolean init);
NcmSpline2d *ncm_spline2d_ref (NcmSpline2d *s2d);
void ncm_spline2d_free (NcmSpline2d *s2d);
void ncm_spline2d_clear (NcmSpline2d **s2d);

void ncm_spline2d_use_acc (NcmSpline2d *s2d, gboolean use_acc);
void ncm_spline2d_set_init (NcmSpline2d *s2d, gboolean init);

NcmSpline *ncm_spline2d_peek_spline (NcmSpline2d *s2d);
NcmVector *ncm_spline2d_peek_xv (NcmSpline2d *s2d);
NcmVector *ncm_spline2d_peek_yv (NcmSpline2d *s2d);
NcmMatrix *ncm_spline2d_peek_zm (NcmSpline2d *s2d);
gsl_interp_accel *ncm_spline2d_peek_acc_x (NcmSpline2d *s2d);
gsl_interp_accel *ncm_spline2d_peek_acc_y (NcmSpline2d *s2d);

gboolean ncm_spline2d_is_init (NcmSpline2d *s2d);
gboolean ncm_spline2d_has_no_stride (NcmSpline2d *s2d);
gboolean ncm_spline2d_using_acc (NcmSpline2d *s2d);

gdouble ncm_spline2d_eval (NcmSpline2d *s2d, gdouble x, gdouble y);

gdouble ncm_spline2d_integ_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
gdouble ncm_spline2d_integ_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
NcmSpline *ncm_spline2d_integ_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu);
NcmSpline *ncm_spline2d_integ_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu);

gdouble ncm_spline2d_integ_dx_spline_val (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
gdouble ncm_spline2d_integ_dy_spline_val (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dxdy_spline_x (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
gdouble ncm_spline2d_integ_dxdy_spline_y (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);

gdouble ncm_spline2d_deriv_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y);
gdouble ncm_spline2d_deriv_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y);
gdouble ncm_spline2d_deriv_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y);
gdouble ncm_spline2d_deriv_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y);
gdouble ncm_spline2d_deriv_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y);

gdouble ncm_spline2dim_integ_total (NcmSpline2d *s2d);

void ncm_spline2d_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res);

G_END_DECLS

#endif /* _NCM_SPLINE2D_H_ */


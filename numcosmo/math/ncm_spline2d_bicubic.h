/***************************************************************************
 *            ncm_spline2d_bicubic.h
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

#ifndef _NCM_SPLINE2D_BICUBIC_H_
#define _NCM_SPLINE2D_BICUBIC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_spline_cubic.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE2D_BICUBIC             (ncm_spline2d_bicubic_get_type ())
#define NCM_SPLINE2D_BICUBIC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE2D_BICUBIC, NcmSpline2dBicubic))
#define NCM_SPLINE2D_BICUBIC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE2D_BICUBIC, NcmSpline2dBicubicClass))
#define NCM_IS_SPLINE2D_BICUBIC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE2D_BICUBIC))
#define NCM_IS_SPLINE2D_BICUBIC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE2D_BICUBIC))
#define NCM_SPLINE2D_BICUBIC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE2D_BICUBIC, NcmSpline2dBicubicClass))

typedef struct _NcmSpline2dBicubicClass NcmSpline2dBicubicClass;
typedef struct _NcmSpline2dBicubic NcmSpline2dBicubic;
typedef struct __NcmSpline2dBicubicOptimizeInt _NcmSpline2dBicubicOptimizeInt;

struct _NcmSpline2dBicubicClass
{
  /*< private >*/
  NcmSpline2dClass parent_class;
};

/**
 * NcmSpline2dBicubicCoeffs:
 *
 * FIXME
 */
typedef struct _NcmSpline2dBicubicCoeffs NcmSpline2dBicubicCoeffs;

struct _NcmSpline2dBicubicCoeffs
{
  /*< private >*/
  gdouble ij[4][4];
};

struct __NcmSpline2dBicubicOptimizeInt
{
  /*< private >*/
  gdouble l;
  gdouble u;
  gboolean init;
  NcmSpline *s;
};

struct _NcmSpline2dBicubic
{
  /*< private >*/
  NcmSpline2d parent_instance;
  NcmSpline *z_x;
  NcmSpline *dzdy_x;
  NcmSpline *z_y;
  NcmSpline2dBicubicCoeffs *bicoeff;
  _NcmSpline2dBicubicOptimizeInt optimize_dx;
  _NcmSpline2dBicubicOptimizeInt optimize_dy;
};

GType ncm_spline2d_bicubic_get_type (void) G_GNUC_CONST;

NcmSpline2d *ncm_spline2d_bicubic_new (NcmSpline *s);
NcmSpline2d *ncm_spline2d_bicubic_notaknot_new (void);

#define NCM_SPLINE2D_BICUBIC_00 (0)
#define NCM_SPLINE2D_BICUBIC_10 (1)
#define NCM_SPLINE2D_BICUBIC_01 (2)
#define NCM_SPLINE2D_BICUBIC_11 (3)
#define NCM_SPLINE2D_BICUBIC_F   (0)
#define NCM_SPLINE2D_BICUBIC_FX  (1)
#define NCM_SPLINE2D_BICUBIC_FY  (2)
#define NCM_SPLINE2D_BICUBIC_FXY (3)
#define NCM_SPLINE2D_BICUBIC_COEFF_INDEX(s2d,i,j) ((s2d)->z_x->len * (i) + (j))
#define NCM_SPLINE2D_BICUBIC_COEFF(s2d,i,j) ((s2d)->bicoeff[NCM_SPLINE2D_BICUBIC_COEFF_INDEX(s2d,i,j)].ij)
#define NCM_SPLINE2D_BICUBIC_STRUCT(s2d,i,j) ((s2d)->bicoeff[NCM_SPLINE2D_BICUBIC_COEFF_INDEX(s2d,i,j)])

/* Utilities -- internal use */

NCM_INLINE gdouble ncm_spline2d_bicubic_eval_poly (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y);
NCM_INLINE void ncm_spline2d_bicubic_fij_to_aij (NcmSpline2dBicubicCoeffs *sf, const gdouble dx, const gdouble dy, NcmSpline2dBicubicCoeffs *sa);
NCM_INLINE gdouble ncm_spline2d_bicubic_bi (NcmSplineCubic *sc, NcmVector *xv, NcmVector *yv, gsize i);
NCM_INLINE void ncm_spline2d_bicubic_bi_bip1 (NcmSplineCubic *sc, NcmVector *xv, NcmVector *yv, gsize i, gdouble *b_i, gdouble *b_ip1);
NCM_INLINE void ncm_spline2d_bicubic_integ_dx_coeffs (NcmSpline2dBicubicCoeffs *aij, gdouble dy, gdouble *coeffs);
NCM_INLINE void ncm_spline2d_bicubic_integ_dy_coeffs (NcmSpline2dBicubicCoeffs *aij, gdouble dx, gdouble *coeffs);
NCM_INLINE gdouble ncm_spline2d_bicubic_integ_eval2d (NcmSpline2dBicubicCoeffs *aij, const gdouble x0, const gdouble xl, const gdouble xu, const gdouble y0, const gdouble yl, const gdouble yu);

G_END_DECLS

#endif /* _NCM_SPLINE2D_BICUBIC_H_ */

#ifndef _NCM_SPLINE2D_BICUBIC_INLINE_H_
#define _NCM_SPLINE2D_BICUBIC_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
ncm_spline2d_bicubic_eval_poly (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y)
{
  const gdouble a0 = (sa->ij[0][0] + x * (sa->ij[1][0] + x * (sa->ij[2][0] + x * sa->ij[3][0])));
  const gdouble a1 = (sa->ij[0][1] + x * (sa->ij[1][1] + x * (sa->ij[2][1] + x * sa->ij[3][1])));
  const gdouble a2 = (sa->ij[0][2] + x * (sa->ij[1][2] + x * (sa->ij[2][2] + x * sa->ij[3][2])));
  const gdouble a3 = (sa->ij[0][3] + x * (sa->ij[1][3] + x * (sa->ij[2][3] + x * sa->ij[3][3])));
  return a0 + y * (a1 + y * (a2 + y * a3));
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_eval_poly_dzdx (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y)
{
  const gdouble a0 = (sa->ij[1][0] + x * (2.0 * sa->ij[2][0] + x * 3.0 * sa->ij[3][0]));
  const gdouble a1 = (sa->ij[1][1] + x * (2.0 * sa->ij[2][1] + x * 3.0 * sa->ij[3][1]));
  const gdouble a2 = (sa->ij[1][2] + x * (2.0 * sa->ij[2][2] + x * 3.0 * sa->ij[3][2]));
  const gdouble a3 = (sa->ij[1][3] + x * (2.0 * sa->ij[2][3] + x * 3.0 * sa->ij[3][3]));
  return a0 + y * (a1 + y * (a2 + y * a3));
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_eval_poly_dzdy (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y)
{
  const gdouble a1 = (sa->ij[0][1] + x * (sa->ij[1][1] + x * (sa->ij[2][1] + x * sa->ij[3][1])));
  const gdouble a2 = (sa->ij[0][2] + x * (sa->ij[1][2] + x * (sa->ij[2][2] + x * sa->ij[3][2])));
  const gdouble a3 = (sa->ij[0][3] + x * (sa->ij[1][3] + x * (sa->ij[2][3] + x * sa->ij[3][3])));
  return a1 + y * (2.0 * a2 + y * 3.0 * a3);
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_eval_poly_d2zdxy (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y)
{
  const gdouble a1 = (sa->ij[1][1] + x * (2.0 * sa->ij[2][1] + x * 3.0 * sa->ij[3][1]));
  const gdouble a2 = (sa->ij[1][2] + x * (2.0 * sa->ij[2][2] + x * 3.0 * sa->ij[3][2]));
  const gdouble a3 = (sa->ij[1][3] + x * (2.0 * sa->ij[2][3] + x * 3.0 * sa->ij[3][3]));
  return a1 + y * (2.0 * a2 + y * 3.0 * a3);
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_eval_poly_d2zdx2 (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y)
{
  const gdouble a0 = (2.0 * sa->ij[2][0] + x * 6.0 * sa->ij[3][0]);
  const gdouble a1 = (2.0 * sa->ij[2][1] + x * 6.0 * sa->ij[3][1]);
  const gdouble a2 = (2.0 * sa->ij[2][2] + x * 6.0 * sa->ij[3][2]);
  const gdouble a3 = (2.0 * sa->ij[2][3] + x * 6.0 * sa->ij[3][3]);
  return a0 + y * (a1 + y * (a2 + y * a3));
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_eval_poly_d2zdy2 (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y)
{
  const gdouble a2 = (sa->ij[0][2] + x * (sa->ij[1][2] + x * (sa->ij[2][2] + x * sa->ij[3][2])));
  const gdouble a3 = (sa->ij[0][3] + x * (sa->ij[1][3] + x * (sa->ij[2][3] + x * sa->ij[3][3])));
  return 2.0 * a2 + y * 6.0 * a3;
}

NCM_INLINE void
ncm_spline2d_bicubic_fij_to_aij (NcmSpline2dBicubicCoeffs *sf, const gdouble dx, const gdouble dy, NcmSpline2dBicubicCoeffs *sa)
{
  gdouble (*a)[4] = sa->ij;
  gdouble (*fij)[4] = sf->ij;

  const gdouble dx2 = dx  * dx;
  const gdouble dx3 = dx2 * dx;
  const gdouble dy2 = dy  * dy;
  const gdouble dy3 = dy2 * dy;

  const gdouble f_00 = fij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F];
  const gdouble f_x0 = fij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_F];
  const gdouble f_0y = fij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F];
  const gdouble f_xy = fij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_F];

  const gdouble fx_00 = fij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX];
  const gdouble fx_x0 = fij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FX];
  const gdouble fx_0y = fij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX];
  const gdouble fx_xy = fij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FX];

  const gdouble fy_00 = fij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY];
  const gdouble fy_x0 = fij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FY];
  const gdouble fy_0y = fij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY];
  const gdouble fy_xy = fij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FY];

  const gdouble fxy_00 = fij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY];
  const gdouble fxy_x0 = fij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FXY];
  const gdouble fxy_0y = fij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY];
  const gdouble fxy_xy = fij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FXY];

  a[0][0] = f_00;
  a[1][0] = fx_00;
  a[2][0] = (3.0 * (-a[0][0] + f_x0) - (2.0 * a[1][0] + fx_x0) * dx) / dx2;
  a[3][0] = (2.0 * ( a[0][0] - f_x0) + (a[1][0] + fx_x0) * dx) / dx3;

  a[0][1] = fy_00;
  a[1][1] = fxy_00;
  a[2][1] = (3.0 * (fy_x0 - a[0][1]) - (2.0 * a[1][1] + fxy_x0) * dx) / dx2;
  a[3][1] = (2.0 * (a[0][1] - fy_x0) + (a[1][1] + fxy_x0) * dx) / dx3;

  a[0][2] = (3.0 * (-a[0][0] + f_0y) - (2.0 * a[0][1] + fy_0y) * dy) / dy2;
  a[1][2] = (3.0 * (fx_0y - a[1][0]) - (2.0 * a[1][1] + fxy_0y) * dy) / dy2;
  a[2][2] = (9.0 * (f_00 - f_x0 - f_0y + f_xy) +
             3.0 * (2.0 * fx_00 + fx_x0 - 2.0 * fx_0y - fx_xy) * dx +
             3.0 * (2.0 * fy_00 - 2.0 * fy_x0 + fy_0y - fy_xy) * dy +
             (4.0 * fxy_00 + 2.0 * fxy_x0 + 2.0 * fxy_0y + fxy_xy) * dx * dy) / (dx2 * dy2);
  a[3][2] = (6.0 * (-f_00 + f_x0 + f_0y - f_xy) +
             3.0 * (-fx_00 - fx_x0 + fx_0y + fx_xy) * dx +
             2.0 * (-2.0 * fy_00 + 2.0 * fy_x0 - fy_0y + fy_xy) * dy -
             (2.0 * fxy_00 + 2.0 * fxy_x0 + fxy_0y + fxy_xy) * dx * dy) / (dx3 * dy2);

  a[0][3] = (2.0 * (a[0][0] - f_0y) + (a[0][1] + fy_0y) * dy) / dy3;
  a[1][3] = (2.0 * (a[1][0] - fx_0y) + (a[1][1] + fxy_0y) * dy) / dy3;
  a[2][3] = (6.0 * (-f_00 + f_x0 + f_0y - f_xy) +
             2.0 * (-2.0 * fx_00 - fx_x0 + 2.0 * fx_0y + fx_xy) * dx +
             3.0 * (-fy_00 + fy_x0 - fy_0y + fy_xy) * dy -
             (2.0 * fxy_00 + fxy_x0 + 2.0 * fxy_0y + fxy_xy) * dx * dy) / (dx2 * dy3);
  a[3][3] = (4.0 * (f_00 - f_x0 - f_0y + f_xy) +
             2.0 * (fx_00 + fx_x0 - fx_0y - fx_xy) * dx +
             2.0 * (fy_00 - fy_x0 + fy_0y - fy_xy) * dy +
             (fxy_00 + fxy_x0 + fxy_0y + fxy_xy) * dx * dy) / (dx3 * dy3);
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_bi (NcmSplineCubic *sc, NcmVector *xv, NcmVector *yv, gsize i)
{
  const gdouble dx = ncm_vector_get (xv, i + 1) - ncm_vector_get (xv, i);
  const gdouble dy = ncm_vector_get (yv, i + 1) - ncm_vector_get (yv, i);
  const gdouble c_ip1 = ncm_vector_get (sc->c, i + 1);
  const gdouble c_i = ncm_vector_get (sc->c, i);
  const gdouble dy_dx = (dy / dx);
  const gdouble b_i = dy_dx - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  return b_i;
}

NCM_INLINE void
ncm_spline2d_bicubic_bi_bip1 (NcmSplineCubic *sc, NcmVector *xv, NcmVector *yv, gsize i, gdouble *b_i, gdouble *b_ip1)
{
  const gdouble dx = ncm_vector_get (xv, i + 1) - ncm_vector_get (xv, i);
  const gdouble dy = ncm_vector_get (yv, i + 1) - ncm_vector_get (yv, i);
  const gdouble c_ip1 = ncm_vector_get (sc->c, i + 1);
  const gdouble c_i = ncm_vector_get (sc->c, i);
  const gdouble dy_dx = (dy / dx);
  *b_i = dy_dx - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *b_ip1 = 3.0 * dy_dx - 2.0 * (*b_i) - c_i * dx;
}

NCM_INLINE void
ncm_spline2d_bicubic_integ_dx_coeffs (NcmSpline2dBicubicCoeffs *aij, gdouble dy, gdouble *coeffs)
{
  gdouble dy2 = dy * dy;
  gdouble dy3 = dy2 * dy;

  coeffs[0] = aij->ij[0][0] + aij->ij[0][1] * dy + aij->ij[0][2] * dy2 + aij->ij[0][3] * dy3;
  coeffs[1] = aij->ij[1][0] + aij->ij[1][1] * dy + aij->ij[1][2] * dy2 + aij->ij[1][3] * dy3;
  coeffs[2] = aij->ij[2][0] + aij->ij[2][1] * dy + aij->ij[2][2] * dy2 + aij->ij[2][3] * dy3;
  coeffs[3] = aij->ij[3][0] + aij->ij[3][1] * dy + aij->ij[3][2] * dy2 + aij->ij[3][3] * dy3;
}

NCM_INLINE void
ncm_spline2d_bicubic_integ_dy_coeffs (NcmSpline2dBicubicCoeffs *aij, gdouble dx, gdouble *coeffs)
{
  gdouble dx2 = dx * dx;
  gdouble dx3 = dx2 * dx;

  coeffs[0] = aij->ij[0][0] + aij->ij[1][0] * dx + aij->ij[2][0] * dx2 + aij->ij[3][0] * dx3;
  coeffs[1] = aij->ij[0][1] + aij->ij[1][1] * dx + aij->ij[2][1] * dx2 + aij->ij[3][1] * dx3;
  coeffs[2] = aij->ij[0][2] + aij->ij[1][2] * dx + aij->ij[2][2] * dx2 + aij->ij[3][2] * dx3;
  coeffs[3] = aij->ij[0][3] + aij->ij[1][3] * dx + aij->ij[2][3] * dx2 + aij->ij[3][3] * dx3;
}

NCM_INLINE gdouble
ncm_spline2d_bicubic_integ_eval2d (NcmSpline2dBicubicCoeffs *aij, const gdouble x0, const gdouble xl, const gdouble xu, const gdouble y0, const gdouble yl, const gdouble yu)
{
  const gdouble dxl  = xl- x0;
  const gdouble dxl2 = dxl * dxl;
  const gdouble dxl3 = dxl2 * dxl;
  const gdouble dxl4 = dxl3 * dxl;

  const gdouble dxu  = xu- x0;
  const gdouble dxu2 = dxu * dxu;
  const gdouble dxu3 = dxu2 * dxu;
  const gdouble dxu4 = dxu3 * dxu;

  const gdouble dyl  = yl- y0;
  const gdouble dyl2 = dyl * dyl;
  const gdouble dyl3 = dyl2 * dyl;
  const gdouble dyl4 = dyl3 * dyl;

  const gdouble dyu  = yu- y0;
  const gdouble dyu2 = dyu * dyu;
  const gdouble dyu3 = dyu2 * dyu;
  const gdouble dyu4 = dyu3 * dyu;

  const gdouble dxu_dxl   =  dxu - dxl;
  const gdouble dxu2_dxl2 =  0.5 * (dxu2 - dxl2);
  const gdouble dxu3_dxl3 =  (dxu3 - dxl3) / 3.0;
  const gdouble dxu4_dxl4 =  0.25 * (dxu4 - dxl4);

  gdouble result = (aij->ij[0][0] * dxu_dxl + aij->ij[1][0] * dxu2_dxl2 + aij->ij[2][0] * dxu3_dxl3 + aij->ij[3][0] * dxu4_dxl4) * (dyu - dyl)
	+ (aij->ij[0][1] * dxu_dxl + aij->ij[1][1] * dxu2_dxl2 + aij->ij[2][1] * dxu3_dxl3 + aij->ij[3][1] * dxu4_dxl4) * 0.5 * (dyu2 - dyl2)
	+ (aij->ij[0][2] * dxu_dxl + aij->ij[1][2] * dxu2_dxl2 + aij->ij[2][2] * dxu3_dxl3 + aij->ij[3][2] * dxu4_dxl4) * (dyu3 - dyl3) / 3.0
	+ (aij->ij[0][3] * dxu_dxl + aij->ij[1][3] * dxu2_dxl2 + aij->ij[2][3] * dxu3_dxl3 + aij->ij[3][3] * dxu4_dxl4) * 0.25 * (dyu4 - dyl4);

  return result;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SPLINE2D_BICUBIC_INLINE_H_ */

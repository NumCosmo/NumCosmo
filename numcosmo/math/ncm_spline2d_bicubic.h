/***************************************************************************
 *            ncm_spline2d_bicubic.h
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
  guint z_x_len;
  NcmSpline2dBicubicCoeffs *bicoeff;
  _NcmSpline2dBicubicOptimizeInt optimize_dx;
  _NcmSpline2dBicubicOptimizeInt optimize_dy;
};

GType ncm_spline2d_bicubic_get_type (void) G_GNUC_CONST;

NcmSpline2d *ncm_spline2d_bicubic_new (NcmSpline *s);
NcmSpline2d *ncm_spline2d_bicubic_notaknot_new (void);

/* Utilities -- internal use */

gdouble ncm_spline2d_bicubic_eval_poly (const NcmSpline2dBicubicCoeffs *sa, const gdouble x, const gdouble y);
void ncm_spline2d_bicubic_fij_to_aij (NcmSpline2dBicubicCoeffs *sf, const gdouble dx, const gdouble dy, NcmSpline2dBicubicCoeffs *sa);
gdouble ncm_spline2d_bicubic_bi (NcmSplineCubic *sc, NcmVector *xv, NcmVector *yv, gsize i);
void ncm_spline2d_bicubic_bi_bip1 (NcmSplineCubic *sc, NcmVector *xv, NcmVector *yv, gsize i, gdouble *b_i, gdouble *b_ip1);
void ncm_spline2d_bicubic_integ_dx_coeffs (NcmSpline2dBicubicCoeffs *aij, gdouble dy, gdouble *coeffs);
void ncm_spline2d_bicubic_integ_dy_coeffs (NcmSpline2dBicubicCoeffs *aij, gdouble dx, gdouble *coeffs);
gdouble ncm_spline2d_bicubic_integ_eval2d (NcmSpline2dBicubicCoeffs *aij, const gdouble x0, const gdouble xl, const gdouble xu, const gdouble y0, const gdouble yl, const gdouble yu);

#define NCM_SPLINE2D_BICUBIC_00 (0)
#define NCM_SPLINE2D_BICUBIC_10 (1)
#define NCM_SPLINE2D_BICUBIC_01 (2)
#define NCM_SPLINE2D_BICUBIC_11 (3)
#define NCM_SPLINE2D_BICUBIC_F   (0)
#define NCM_SPLINE2D_BICUBIC_FX  (1)
#define NCM_SPLINE2D_BICUBIC_FY  (2)
#define NCM_SPLINE2D_BICUBIC_FXY (3)

G_END_DECLS

#endif /* _NCM_SPLINE2D_BICUBIC_H_ */


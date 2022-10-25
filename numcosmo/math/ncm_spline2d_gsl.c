/***************************************************************************
 *            ncm_spline2d_gsl.c
 *
 *  Sun Aug  1 17:17:08 2010
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

/**
 * SECTION:ncm_spline2d_gsl
 * @title: NcmSpline2dGsl
 * @short_description: Implements spline from spline method using The GNU Scientific Library (GSL) as base splines.
 * @stability: Stable
 * @include: numcosmo/math/ncm_spline2d_gsl.h
 *
 * This object implements bidimensional splines with the method
 * given by the #NcmSplineGsl class.
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline2d_gsl.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcmSpline2dGsl, ncm_spline2d_gsl, NCM_TYPE_SPLINE2D);

static void
ncm_spline2d_gsl_init (NcmSpline2dGsl *s2dgsl)
{
  s2dgsl->zdiff       = NULL;
  s2dgsl->vertv       = NULL;
  s2dgsl->vertintv    = NULL;
  s2dgsl->s_hor       = NULL;
  s2dgsl->s_dzdy      = NULL;
  s2dgsl->s_ver       = NULL;
  s2dgsl->s_ver_integ = NULL;
  s2dgsl->s_hor_len   = 0;
}

static void
_ncm_spline2d_gsl_clear (NcmSpline2dGsl *s2dgsl)
{
  guint i;

  ncm_matrix_clear (&s2dgsl->zdiff);
  ncm_vector_clear (&s2dgsl->vertv);
  ncm_vector_clear (&s2dgsl->vertintv);

  for (i = 0; i < s2dgsl->s_hor_len; i++)
    ncm_spline_clear (&s2dgsl->s_hor[i]);

  for (i = 0; i < s2dgsl->s_hor_len; i++)
    ncm_spline_clear (&s2dgsl->s_dzdy[i]);

  ncm_spline_clear (&s2dgsl->s_ver);
  ncm_spline_clear (&s2dgsl->s_ver_integ);
}

static void
_ncm_spline2d_gsl_dispose (GObject *object)
{
  NcmSpline2d *s2d       = NCM_SPLINE2D (object);
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (object);

  _ncm_spline2d_gsl_clear (s2dgsl);
  s2d->init = FALSE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_gsl_parent_class)->dispose (object);
}

static void
_ncm_spline2d_gsl_finalize (GObject *object)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (object);

  if (s2dgsl->s_hor != NULL)
  {
    g_free (s2dgsl->s_hor);
    s2dgsl->s_hor = NULL;
  }
  
  if (s2dgsl->s_dzdy != NULL)
  {
    g_free (s2dgsl->s_dzdy);
    s2dgsl->s_dzdy = NULL;
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_gsl_parent_class)->finalize (object);
}

NcmSpline2d *_ncm_spline2d_gsl_copy_empty (const NcmSpline2d *s2d);
static void _ncm_spline2d_gsl_reset (NcmSpline2d *s2d);
static void _ncm_spline2d_gsl_prepare (NcmSpline2d *s2d);
static gdouble _ncm_spline2d_gsl_eval (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_gsl_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_gsl_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_gsl_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_gsl_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_gsl_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_gsl_int_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
static gdouble _ncm_spline2d_gsl_int_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
static gdouble _ncm_spline2d_gsl_int_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
static NcmSpline *_ncm_spline2d_gsl_int_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu);
static NcmSpline *_ncm_spline2d_gsl_int_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu);

static void
ncm_spline2d_gsl_class_init (NcmSpline2dGslClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmSpline2dClass *parent_class = NCM_SPLINE2D_CLASS (klass);
  
  object_class->dispose  = &_ncm_spline2d_gsl_dispose;
  object_class->finalize = &_ncm_spline2d_gsl_finalize;

  parent_class->copy_empty    = &_ncm_spline2d_gsl_copy_empty;
  parent_class->reset         = &_ncm_spline2d_gsl_reset;
  parent_class->prepare       = &_ncm_spline2d_gsl_prepare;
  parent_class->eval          = &_ncm_spline2d_gsl_eval;
  parent_class->dzdx          = &_ncm_spline2d_gsl_dzdx;
  parent_class->dzdy          = &_ncm_spline2d_gsl_dzdy;
  parent_class->d2zdxy        = &_ncm_spline2d_gsl_d2zdxy;
  parent_class->d2zdx2        = &_ncm_spline2d_gsl_d2zdx2;
  parent_class->d2zdy2        = &_ncm_spline2d_gsl_d2zdy2;
  parent_class->int_dx        = &_ncm_spline2d_gsl_int_dx;
  parent_class->int_dy        = &_ncm_spline2d_gsl_int_dy;
  parent_class->int_dxdy      = &_ncm_spline2d_gsl_int_dxdy;
  parent_class->int_dx_spline = &_ncm_spline2d_gsl_int_dx_spline;
  parent_class->int_dy_spline = &_ncm_spline2d_gsl_int_dy_spline;
}

NcmSpline2d *
_ncm_spline2d_gsl_copy_empty (const NcmSpline2d *s2d)
{
  return ncm_spline2d_gsl_new (s2d->s);
}

static void
_ncm_spline2d_gsl_alloc (NcmSpline2dGsl *s2dgsl)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (s2dgsl);
  guint i;
  
  s2dgsl->s_hor_len = ncm_matrix_col_len (s2d->zm);
  
  s2dgsl->zdiff    = ncm_matrix_new (s2dgsl->s_hor_len, ncm_matrix_row_len (s2d->zm));
  s2dgsl->vertv    = ncm_vector_new (s2dgsl->s_hor_len);
  s2dgsl->vertintv = ncm_vector_new (s2dgsl->s_hor_len);
  
  s2dgsl->s_hor = g_new0 (NcmSpline *, s2dgsl->s_hor_len);
  
  for (i = 0; i < s2dgsl->s_hor_len; i++)
  {
    NcmVector *zm_row_i = ncm_matrix_get_row (s2d->zm, i);
    
    s2dgsl->s_hor[i] = ncm_spline_new (s2d->s, s2d->xv, zm_row_i, s2d->init);
    ncm_vector_free (zm_row_i);
  }
  
  s2dgsl->s_ver       = ncm_spline_new (s2d->s, s2d->yv, s2dgsl->vertv, FALSE);
  s2dgsl->s_ver_integ = ncm_spline_new (s2d->s, s2d->yv, s2dgsl->vertintv, FALSE);
  
  s2dgsl->s_dzdy = g_new0 (NcmSpline *, s2dgsl->s_hor_len);
  
  for (i = 0; i < s2dgsl->s_hor_len; i++)
  {
    NcmVector *zdiff_row_i = ncm_matrix_get_row (s2dgsl->zdiff, i);
    
    s2dgsl->s_dzdy[i] = ncm_spline_new (s2d->s, s2d->xv, zdiff_row_i, s2d->init);
    ncm_vector_free (zdiff_row_i);
  }
}

static void
_ncm_spline2d_gsl_free (NcmSpline2dGsl *s2dgsl)
{
  _ncm_spline2d_gsl_clear (s2dgsl);
  
  g_free (s2dgsl->s_hor);
  s2dgsl->s_hor = NULL;
  g_free (s2dgsl->s_dzdy);
  s2dgsl->s_dzdy = NULL;
}

static void
_ncm_spline2d_gsl_reset (NcmSpline2d *s2d)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  if (s2d->init)
  {
    if ((ncm_matrix_nrows (s2d->zm) != ncm_vector_len (s2dgsl->vertv)) ||
        (ncm_matrix_ncols (s2d->zm) != ncm_vector_len (s2dgsl->s_hor[0]->xv)))
    {
      _ncm_spline2d_gsl_free (s2dgsl);
      _ncm_spline2d_gsl_alloc (s2dgsl);
    }
  }
  else
  {
    _ncm_spline2d_gsl_alloc (s2dgsl);
  }
}

static void
_ncm_spline2d_gsl_prepare (NcmSpline2d *s2d)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  guint i, j;
  
  for (i = 0; i < s2dgsl->s_hor_len; i++)
    ncm_spline_prepare (s2dgsl->s_hor[i]);
  
  for (j = 0; j < ncm_matrix_row_len (s2d->zm); j++)
  {
    for (i = 0; i < s2dgsl->s_hor_len; i++)
      ncm_vector_set (s2dgsl->vertv, i, ncm_matrix_get (s2d->zm, i, j));
    
    ncm_spline_prepare (s2dgsl->s_ver);
    
    for (i = 0; i < s2dgsl->s_hor_len; i++)
    {
      const gdouble yi       = ncm_vector_get (s2d->yv, i);
      const gdouble diff_val = ncm_spline_eval_deriv (s2dgsl->s_ver, yi);
      
      ncm_matrix_set (s2dgsl->zdiff, i, j, diff_val);
    }
  }
  
  for (i = 0; i < s2dgsl->s_hor_len; i++)
    ncm_spline_prepare (s2dgsl->s_dzdy[i]);
  
  s2d->init = TRUE;
}

static void
_ncm_bicubic_calc_aij (NcmSpline2dGsl *s2dgsl, gsize i, gsize j, const gdouble x0, const gdouble x1, const gdouble y0, const gdouble y1, NcmSpline2dBicubicCoeffs *a)
{
  const gdouble dx = x1 - x0;
  const gdouble dy = y1 - y0;
  NcmSpline2dBicubicCoeffs f;
  
  NCM_UNUSED (j);
  
  f.ij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F] = ncm_spline_eval (s2dgsl->s_hor[i], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_F] = ncm_spline_eval (s2dgsl->s_hor[i], x1);
  f.ij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F] = ncm_spline_eval (s2dgsl->s_hor[i + 1], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_F] = ncm_spline_eval (s2dgsl->s_hor[i + 1], x1);
  
  f.ij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = ncm_spline_eval (s2dgsl->s_dzdy[i], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FY] = ncm_spline_eval (s2dgsl->s_dzdy[i], x1);
  f.ij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = ncm_spline_eval (s2dgsl->s_dzdy[i + 1], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FY] = ncm_spline_eval (s2dgsl->s_dzdy[i + 1], x1);
  
  f.ij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX] = ncm_spline_eval_deriv (s2dgsl->s_hor[i], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FX] = ncm_spline_eval_deriv (s2dgsl->s_hor[i], x1);
  f.ij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX] = ncm_spline_eval_deriv (s2dgsl->s_hor[i + 1], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FX] = ncm_spline_eval_deriv (s2dgsl->s_hor[i + 1], x1);
  
  f.ij[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY] = ncm_spline_eval_deriv (s2dgsl->s_dzdy[i], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FXY] = ncm_spline_eval_deriv (s2dgsl->s_dzdy[i], x1);
  f.ij[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY] = ncm_spline_eval_deriv (s2dgsl->s_dzdy[i + 1], x0);
  f.ij[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FXY] = ncm_spline_eval_deriv (s2dgsl->s_dzdy[i + 1], x1);
  
  ncm_spline2d_bicubic_fij_to_aij (&f, dx, dy, a);
}

static gdouble
_ncm_spline2d_gsl_eval (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j          = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i          = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs a;
  
  _ncm_bicubic_calc_aij (s2dgsl, i, j, x0, x1, y0, y1, &a);
  
  
  return ncm_spline2d_bicubic_eval_poly (&a, x - x0, y - y0);
}

static gdouble
_ncm_spline2d_gsl_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j          = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i          = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs a;
  gdouble result;
  
  const gdouble dx  = x - x0;
  const gdouble dx2 = dx * dx;
  
  const gdouble dy  = y - y0;
  const gdouble dy2 = dy * dy;
  const gdouble dy3 = dy2 * dy;
  
  _ncm_bicubic_calc_aij (s2dgsl, i, j, x0, x1, y0, y1, &a);
  
  result = a.ij[1][0] + a.ij[1][1] * dy + a.ij[1][2] * dy2 + a.ij[1][3] * dy3
           + 2.0 * (a.ij[2][0] + a.ij[2][1] * dy + a.ij[2][2] * dy2 + a.ij[2][3] * dy3) * dx
           + 3.0 * (a.ij[3][0] + a.ij[3][1] * dy + a.ij[3][2] * dy2 + a.ij[3][3] * dy3) * dx2;
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j          = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i          = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs a;
  gdouble result;
  
  const gdouble dx  = x - x0;
  const gdouble dx2 = dx * dx;
  const gdouble dx3 = dx2 * dx;
  
  const gdouble dy  = y - y0;
  const gdouble dy2 = dy * dy;
  
  _ncm_bicubic_calc_aij (s2dgsl, i, j, x0, x1, y0, y1, &a);
  
  result = a.ij[0][1] + a.ij[1][1] * dx + a.ij[2][1] * dx2 + a.ij[3][1] * dx3
           + 2.0 * (a.ij[0][2] + a.ij[1][2] * dx + a.ij[2][2] * dx2 + a.ij[3][2] * dx3) * dy
           + 3.0 * (a.ij[0][3] + a.ij[1][3] * dx + a.ij[2][3] * dx2 + a.ij[3][3] * dx3) * dy2;
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j          = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i          = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs a;
  gdouble result;
  
  const gdouble dx = x - x0;
  
  const gdouble dy  = y - y0;
  const gdouble dy2 = dy * dy;
  const gdouble dy3 = dy2 * dy;
  
  _ncm_bicubic_calc_aij (s2dgsl, i, j, x0, x1, y0, y1, &a);
  
  result = 2.0 * (a.ij[2][0] + a.ij[2][1] * dy + a.ij[2][2] * dy2 + a.ij[2][3] * dy3)
           + 6.0 * (a.ij[3][0] + a.ij[3][1] * dy + a.ij[3][2] * dy2 + a.ij[3][3] * dy3) * dx;
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j          = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i          = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs a;
  gdouble result;
  
  const gdouble dx  = x - x0;
  const gdouble dx2 = dx * dx;
  const gdouble dx3 = dx2 * dx;
  
  const gdouble dy = y - y0;
  
  _ncm_bicubic_calc_aij (s2dgsl, i, j, x0, x1, y0, y1, &a);
  
  result = 2.0 * (a.ij[0][2] + a.ij[1][2] * dx + a.ij[2][2] * dx2 + a.ij[3][2] * dx3)
           + 6.0 * (a.ij[0][3] + a.ij[1][3] * dx + a.ij[2][3] * dx2 + a.ij[3][3] * dx3) * dy;
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j          = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i          = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs a;
  gdouble result;
  
  const gdouble dx  = x - x0;
  const gdouble dx2 = dx * dx;
  
  const gdouble dy  = y - y0;
  const gdouble dy2 = dy * dy;
  
  _ncm_bicubic_calc_aij (s2dgsl, i, j, x0, x1, y0, y1, &a);
  
  result = a.ij[1][1] + 2.0 * (a.ij[2][1] * dx + a.ij[1][2] * dy) + 4.0 * a.ij[2][2] * dx * dy
           + 3.0 * (a.ij[3][1] * dx2 + a.ij[1][3] * dy2) + 9.0 * a.ij[3][3] * dx2 * dy2
           + 6.0 * (a.ij[3][2] * dx2 * dy + a.ij[2][3] * dx * dy2);
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_int_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize jl = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xl, 0, ncm_vector_len (s2d->xv) - 1);
  gsize ju = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xu, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  gdouble x0, x1, result;
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  const gdouble y1 = ncm_vector_get (s2d->yv, i + 1);
  NcmSpline2dBicubicCoeffs aij;
  guint k;
  gdouble coeffs[4];
  
  g_assert (jl <= ju);
  
  if (jl == ju)
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    x1 = ncm_vector_get (s2d->xv, jl + 1);
    _ncm_bicubic_calc_aij (s2dgsl, i, jl, x0, x1, y0, y1, &aij);
    ncm_spline2d_bicubic_integ_dx_coeffs (&aij, y - y0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], x0, xl, xu);
  }
  else
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    x1 = ncm_vector_get (s2d->xv, jl + 1);
    _ncm_bicubic_calc_aij (s2dgsl, i, jl, x0, x1, y0, y1, &aij);
    ncm_spline2d_bicubic_integ_dx_coeffs (&aij, y - y0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], x0, xl, x1);
    
    for (k = jl + 1; k < ju; k++)
    {
      x0 = ncm_vector_get (s2d->xv, k);
      x1 = ncm_vector_get (s2d->xv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, i, k, x0, x1, y0, y1, &aij);
      ncm_spline2d_bicubic_integ_dx_coeffs (&aij, y - y0, coeffs);
      
      {
        const gdouble dx = x1 - x0;
        
        result += dx * (coeffs[0] + dx * (0.5 * coeffs[1] + dx * (coeffs[2] / 3.0 + 0.25 * coeffs[3] * dx)));
      }
    }
    
    k = ju;
    {
      x0 = ncm_vector_get (s2d->xv, k);
      x1 = ncm_vector_get (s2d->xv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, i, k, x0, x1, y0, y1, &aij);
      ncm_spline2d_bicubic_integ_dx_coeffs (&aij, y - y0, coeffs);
      result += _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], x0, x0, xu);
    }
  }
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_int_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize j = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize il = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yl, 0, ncm_vector_len (s2d->yv) - 1);
  gsize iu = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yu, 0, ncm_vector_len (s2d->yv) - 1);
  gdouble y0, y1, result;
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  const gdouble x1 = ncm_vector_get (s2d->xv, j + 1);
  NcmSpline2dBicubicCoeffs aij;
  guint k;
  gdouble coeffs[4];
  
  g_assert (il <= iu);
  
  if (il == iu)
  {
    y0 = ncm_vector_get (s2d->yv, il);
    y1 = ncm_vector_get (s2d->yv, il + 1);
    _ncm_bicubic_calc_aij (s2dgsl, il, j, x0, x1, y0, y1, &aij);
    ncm_spline2d_bicubic_integ_dy_coeffs (&aij, x - x0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], y0, yl, yu);
  }
  else
  {
    y0 = ncm_vector_get (s2d->yv, il);
    y1 = ncm_vector_get (s2d->yv, il + 1);
    _ncm_bicubic_calc_aij (s2dgsl, il, j, x0, x1, y0, y1, &aij);
    ncm_spline2d_bicubic_integ_dy_coeffs (&aij, x - x0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], y0, yl, y1);
    
    for (k = il + 1; k < iu; k++)
    {
      y0 = ncm_vector_get (s2d->yv, k);
      y1 = ncm_vector_get (s2d->yv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, k, j, x0, x1, y0, y1, &aij);
      ncm_spline2d_bicubic_integ_dy_coeffs (&aij, x - x0, coeffs);
      
      {
        const gdouble dy = y1 - y0;
        
        result += dy * (coeffs[0] + dy * (0.5 * coeffs[1] + dy * (coeffs[2] / 3.0 + 0.25 * coeffs[3] * dy)));
      }
    }
    
    k = iu;
    {
      y0 = ncm_vector_get (s2d->yv, k);
      y1 = ncm_vector_get (s2d->yv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, k, j, x0, x1, y0, y1, &aij);
      ncm_spline2d_bicubic_integ_dy_coeffs (&aij, x - x0, coeffs);
      result += _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], y0, y0, yu);
    }
  }
  
  return result;
}

static gdouble
_ncm_spline2d_gsl_int_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  NcmSpline2dGsl *s2dgsl = NCM_SPLINE2D_GSL (s2d);
  
  gsize jl = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xl, 0, ncm_vector_len (s2d->xv) - 1);
  gsize ju = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xu, 0, ncm_vector_len (s2d->xv) - 1);
  gsize il = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yl, 0, ncm_vector_len (s2d->yv) - 1);
  gsize iu = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yu, 0, ncm_vector_len (s2d->yv) - 1);
  gdouble x0, x1, y0, y1, result;
  NcmSpline2dBicubicCoeffs aij;
  guint k, m;
  
  g_assert (jl <= ju || il <= iu);
  
  if ((jl == ju) && (il == iu))
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    x1 = ncm_vector_get (s2d->xv, jl + 1);
    y0 = ncm_vector_get (s2d->yv, il);
    y1 = ncm_vector_get (s2d->yv, il + 1);
    _ncm_bicubic_calc_aij (s2dgsl, il, jl, x0, x1, y0, y1, &aij);
    result = ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, xu, y0, yl, yu);
  }
  else if (jl == ju)
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    x1 = ncm_vector_get (s2d->xv, jl + 1);
    y0 = ncm_vector_get (s2d->yv, il);
    y1 = ncm_vector_get (s2d->yv, il + 1);
    _ncm_bicubic_calc_aij (s2dgsl, il, jl, x0, x1, y0, y1, &aij);
    result = ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, xu, y0, yl, y1);
    
    for (k = il + 1; k < iu; k++)
    {
      y0 = ncm_vector_get (s2d->yv, k);
      y1 = ncm_vector_get (s2d->yv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, k, jl, x0, x1, y0, y1, &aij);
      result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, xu, y0, y0, y1);
    }
    
    k = iu;
    {
      y0 = ncm_vector_get (s2d->yv, k);
      y1 = ncm_vector_get (s2d->yv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, k, jl, x0, x1, y0, y1, &aij);
      result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, xu, y0, y0, yu);
    }
  }
  else if (il == iu)
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    x1 = ncm_vector_get (s2d->xv, jl + 1);
    y0 = ncm_vector_get (s2d->yv, il);
    y1 = ncm_vector_get (s2d->yv, il + 1);
    _ncm_bicubic_calc_aij (s2dgsl, il, jl, x0, x1, y0, y1, &aij);
    result = ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, x1, y0, yl, yu);
    
    for (k = jl + 1; k < ju; k++)
    {
      x0 = ncm_vector_get (s2d->xv, k);
      x1 = ncm_vector_get (s2d->xv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, il, k, x0, x1, y0, y1, &aij);
      result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, x1, y0, yl, yu);
    }
    
    k = ju;
    {
      x0 = ncm_vector_get (s2d->xv, k);
      x1 = ncm_vector_get (s2d->xv, k + 1);
      _ncm_bicubic_calc_aij (s2dgsl, il, k, x0, x1, y0, y1, &aij);
      result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, xu, y0, yl, yu);
    }
  }
  else
  {
    m = jl;
    {
      x0 = ncm_vector_get (s2d->xv, jl);
      x1 = ncm_vector_get (s2d->xv, jl + 1);
      y0 = ncm_vector_get (s2d->yv, il);
      y1 = ncm_vector_get (s2d->yv, il + 1);
      _ncm_bicubic_calc_aij (s2dgsl, il, m, x0, x1, y0, y1, &aij);
      result = ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, x1, y0, yl, y1);
      
      for (k = il + 1; k < iu; k++)
      {
        y0 = ncm_vector_get (s2d->yv, k);
        y1 = ncm_vector_get (s2d->yv, k + 1);
        _ncm_bicubic_calc_aij (s2dgsl, k, m, x0, x1, y0, y1, &aij);
        result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, x1, y0, y0, y1);
      }
      
      k = iu;
      {
        y0 = ncm_vector_get (s2d->yv, k);
        y1 = ncm_vector_get (s2d->yv, k + 1);
        _ncm_bicubic_calc_aij (s2dgsl, k, m, x0, x1, y0, y1, &aij);
        result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, xl, x1, y0, y0, yu);
      }
    }
    
    for (m = jl + 1; m < ju; m++)
    {
      x0 = ncm_vector_get (s2d->xv, m);
      x1 = ncm_vector_get (s2d->xv, m + 1);
      y0 = ncm_vector_get (s2d->yv, il);
      y1 = ncm_vector_get (s2d->yv, il + 1);
      _ncm_bicubic_calc_aij (s2dgsl, il, m, x0, x1, y0, y1, &aij);
      result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, x1, y0, yl, y1);
      
      for (k = il + 1; k < iu; k++)
      {
        y0 = ncm_vector_get (s2d->yv, k);
        y1 = ncm_vector_get (s2d->yv, k + 1);
        _ncm_bicubic_calc_aij (s2dgsl, k, m, x0, x1, y0, y1, &aij);
        result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, x1, y0, y0, y1);
      }
      
      k = iu;
      {
        y0 = ncm_vector_get (s2d->yv, k);
        y1 = ncm_vector_get (s2d->yv, k + 1);
        _ncm_bicubic_calc_aij (s2dgsl, k, m, x0, x1, y0, y1, &aij);
        result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, x1, y0, y0, yu);
      }
    }
    
    m = ju;
    {
      x0 = ncm_vector_get (s2d->xv, m);
      x1 = ncm_vector_get (s2d->xv, m + 1);
      y0 = ncm_vector_get (s2d->yv, il);
      y1 = ncm_vector_get (s2d->yv, il + 1);
      _ncm_bicubic_calc_aij (s2dgsl, il, m, x0, x1, y0, y1, &aij);
      result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, xu, y0, yl, y1);
      
      for (k = il + 1; k < iu; k++)
      {
        y0 = ncm_vector_get (s2d->yv, k);
        y1 = ncm_vector_get (s2d->yv, k + 1);
        _ncm_bicubic_calc_aij (s2dgsl, k, m, x0, x1, y0, y1, &aij);
        result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, xu, y0, y0, y1);
      }
      
      k = iu;
      {
        y0 = ncm_vector_get (s2d->yv, k);
        y1 = ncm_vector_get (s2d->yv, k + 1);
        _ncm_bicubic_calc_aij (s2dgsl, k, m, x0, x1, y0, y1, &aij);
        result += ncm_spline2d_bicubic_integ_eval2d (&aij, x0, x0, xu, y0, y0, yu);
      }
    }
  }
  
  return result;
}

/* LCOV_EXCL_START */

static NcmSpline *
_ncm_spline2d_gsl_int_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu)
{
  g_assert_not_reached ();
  return NULL;
}

static NcmSpline *
_ncm_spline2d_gsl_int_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu)
{
  g_assert_not_reached ();
  return NULL;
}

/* LCOV_EXCL_STOP */

/**
 * ncm_spline2d_gsl_new:
 * @s: a #NcmSplineGsl derived #NcmSpline
 *
 * This function initializes a #NcmSpline2d of
 * [GSL](https://www.gnu.org/software/gsl/) type given in @s.
 *
 * Returns: A new #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_gsl_new (NcmSpline *s)
{
  NcmSpline2d *s2d;

  g_assert (NCM_IS_SPLINE_GSL (s));

  s2d = g_object_new (NCM_TYPE_SPLINE2D_GSL, "spline", s, NULL);

  return s2d;
}

/**
 * ncm_spline2d_gsl_natural_new:
 *
 * This function initializes a #NcmSpline2d of [GSL](https://www.gnu.org/software/gsl/)
 * type [gsl_interp_cspline](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_cspline).
 *
 * Returns: A new #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_gsl_natural_new ()
{
  NcmSpline *s     = ncm_spline_gsl_new (gsl_interp_cspline);
  NcmSpline2d *s2d = ncm_spline2d_gsl_new (s);

  ncm_spline_free (s);

  return s2d;
}


/***************************************************************************
 *            ncm_spline2d_bicubic.c
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
 * SECTION:ncm_spline2d_bicubic
 * @title: NcmSpline2dBicubic
 * @short_description: Bidimensional bicubic spline
 * @stability: Stable
 * @include: numcosmo/math/ncm_spline2d_bicubic.h
 *
 * This class implements the functions which use a bicubic interpolation method.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcmSpline2dBicubic, ncm_spline2d_bicubic, NCM_TYPE_SPLINE2D);

static void
ncm_spline2d_bicubic_init (NcmSpline2dBicubic *object)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (object);

  s2dbc->z_x     = NULL;
  s2dbc->dzdy_x  = NULL;
  s2dbc->z_y     = NULL;
  s2dbc->bicoeff = NULL;

  s2dbc->optimize_dx.l    = GSL_NAN;
  s2dbc->optimize_dx.u    = GSL_NAN;
  s2dbc->optimize_dx.init = FALSE;
  s2dbc->optimize_dx.s    = NULL;

  s2dbc->optimize_dy.l    = GSL_NAN;
  s2dbc->optimize_dy.u    = GSL_NAN;
  s2dbc->optimize_dy.init = FALSE;
  s2dbc->optimize_dy.s    = NULL;
}

static void _ncm_spline2d_bicubic_clear (NcmSpline2dBicubic *s2dbc);

static void
_ncm_spline2d_bicubic_dispose (GObject *object)
{
  NcmSpline2d *s2d          = NCM_SPLINE2D (object);
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (object);

  _ncm_spline2d_bicubic_clear (s2dbc);
  s2d->init = FALSE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_bicubic_parent_class)->dispose (object);
}

static void
_ncm_spline2d_bicubic_finalize (GObject *object)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (object);

  g_clear_pointer (&s2dbc->bicoeff, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_bicubic_parent_class)->finalize (object);
}

static NcmSpline2d *_ncm_spline2d_bicubic_copy_empty (const NcmSpline2d *s2d);
static void _ncm_spline2d_bicubic_alloc (NcmSpline2dBicubic *s2dbc);
static void _ncm_spline2d_bicubic_reset (NcmSpline2d *s2d);
static void _ncm_spline2d_bicubic_prepare (NcmSpline2d *s2d);
static gdouble _ncm_spline2d_bicubic_eval (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_bicubic_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_bicubic_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_bicubic_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_bicubic_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_bicubic_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y);
static gdouble _ncm_spline2d_bicubic_int_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y);
static gdouble _ncm_spline2d_bicubic_int_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu);
static gdouble _ncm_spline2d_bicubic_int_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu);
static NcmSpline *_ncm_spline2d_bicubic_int_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu);
static NcmSpline *_ncm_spline2d_bicubic_int_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu);
static void _ncm_spline2d_bicubic_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res);

static void
ncm_spline2d_bicubic_class_init (NcmSpline2dBicubicClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmSpline2dClass *parent_class = NCM_SPLINE2D_CLASS (klass);

  object_class->dispose  = &_ncm_spline2d_bicubic_dispose;
  object_class->finalize = &_ncm_spline2d_bicubic_finalize;

  parent_class->copy_empty    = &_ncm_spline2d_bicubic_copy_empty;
  parent_class->reset         = &_ncm_spline2d_bicubic_reset;
  parent_class->prepare       = &_ncm_spline2d_bicubic_prepare;
  parent_class->eval          = &_ncm_spline2d_bicubic_eval;
  parent_class->dzdx          = &_ncm_spline2d_bicubic_dzdx;
  parent_class->dzdy          = &_ncm_spline2d_bicubic_dzdy;
  parent_class->d2zdxy        = &_ncm_spline2d_bicubic_d2zdxy;
  parent_class->d2zdx2        = &_ncm_spline2d_bicubic_d2zdx2;
  parent_class->d2zdy2        = &_ncm_spline2d_bicubic_d2zdy2;
  parent_class->int_dx        = &_ncm_spline2d_bicubic_int_dx;
  parent_class->int_dy        = &_ncm_spline2d_bicubic_int_dy;
  parent_class->int_dxdy      = &_ncm_spline2d_bicubic_int_dxdy;
  parent_class->int_dx_spline = &_ncm_spline2d_bicubic_int_dx_spline;
  parent_class->int_dy_spline = &_ncm_spline2d_bicubic_int_dy_spline;
  parent_class->eval_vec_y    = &_ncm_spline2d_bicubic_eval_vec_y;
}

/**
 * ncm_spline2d_bicubic_new:
 * @s: a #NcmSplineCubic derived #NcmSpline
 *
 * This function initializes a #NcmSpline2d of bicubic type given @s.
 *
 * Returns: A new #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_bicubic_new (NcmSpline *s)
{
  NcmSpline2d *s2d;

  g_assert (NCM_IS_SPLINE_CUBIC (s));

  s2d = g_object_new (NCM_TYPE_SPLINE2D_BICUBIC, "spline", s, NULL);

  return s2d;
}

/**
 * ncm_spline2d_bicubic_notaknot_new:
 *
 * This function initializes a #NcmSpline2d of bi-cubic not-a-knot type
 * (See #NcmSplineCubicNotaknot).
 *
 * Returns: A new #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_bicubic_notaknot_new ()
{
  NcmSpline *s     = ncm_spline_cubic_notaknot_new ();
  NcmSpline2d *s2d = ncm_spline2d_bicubic_new (s);

  ncm_spline_free (s);

  return s2d;
}

static NcmSpline2d *
_ncm_spline2d_bicubic_copy_empty (const NcmSpline2d *s2d)
{
  return ncm_spline2d_bicubic_new (s2d->s);
}

static void
_ncm_spline2d_bicubic_alloc (NcmSpline2dBicubic *s2dbc)
{
  NcmSpline2d *s2d        = NCM_SPLINE2D (s2dbc);
  NcmVector *zm_first_row = ncm_matrix_get_row (s2d->zm, 0);
  NcmVector *zm_first_col = ncm_matrix_get_col (s2d->zm, 0);
  NcmVector *dx_yv        = ncm_vector_new (ncm_vector_len (s2d->yv));
  NcmVector *dy_yv        = ncm_vector_new (ncm_vector_len (s2d->xv));

  s2dbc->z_x = ncm_spline_new (s2d->s, s2d->xv, zm_first_row, s2d->init);

  s2dbc->z_y    = ncm_spline_new (s2d->s, s2d->yv, zm_first_col, s2d->init);
  s2dbc->dzdy_x = ncm_spline_new (s2d->s, s2d->xv, zm_first_row, s2d->init);

  s2dbc->bicoeff = g_new0 (NcmSpline2dBicubicCoeffs, (ncm_matrix_col_len (s2d->zm) - 1) * ncm_matrix_row_len (s2d->zm));

  s2dbc->optimize_dx.s = ncm_spline_set (ncm_spline_cubic_notaknot_new (),
                                         s2d->yv, dx_yv, FALSE);
  s2dbc->optimize_dx.init = FALSE;

  s2dbc->optimize_dy.s = ncm_spline_set (ncm_spline_cubic_notaknot_new (),
                                         s2d->xv, dy_yv, FALSE);
  s2dbc->optimize_dy.init = FALSE;

  ncm_vector_free (zm_first_row);
  ncm_vector_free (zm_first_col);
  ncm_vector_free (dx_yv);
  ncm_vector_free (dy_yv);
}

static void
_ncm_spline2d_bicubic_clear (NcmSpline2dBicubic *s2dbc)
{
  s2dbc->optimize_dx.l    = GSL_NAN;
  s2dbc->optimize_dx.u    = GSL_NAN;
  s2dbc->optimize_dx.init = FALSE;

  s2dbc->optimize_dy.l    = GSL_NAN;
  s2dbc->optimize_dy.u    = GSL_NAN;
  s2dbc->optimize_dy.init = FALSE;

  ncm_spline_clear (&s2dbc->z_x);
  ncm_spline_clear (&s2dbc->z_y);
  ncm_spline_clear (&s2dbc->dzdy_x);
  ncm_spline_clear (&s2dbc->optimize_dx.s);
  ncm_spline_clear (&s2dbc->optimize_dy.s);
}

static void
_ncm_spline2d_bicubic_free (NcmSpline2dBicubic *s2dbc)
{
  _ncm_spline2d_bicubic_clear (s2dbc);

  if (s2dbc->bicoeff != NULL)
  {
    g_free (s2dbc->bicoeff);
    s2dbc->bicoeff = NULL;
  }
}

static void
_ncm_spline2d_bicubic_reset (NcmSpline2d *s2d)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (s2d->init)
  {
    if ((ncm_matrix_nrows (s2d->zm) != ncm_vector_len (s2dbc->z_y->xv)) ||
        (ncm_matrix_ncols (s2d->zm) != ncm_vector_len (s2dbc->z_x->xv)))
    {
      _ncm_spline2d_bicubic_free (s2dbc);
      _ncm_spline2d_bicubic_alloc (s2dbc);
    }
  }
  else
  {
    _ncm_spline2d_bicubic_alloc (s2dbc);
  }
}

static void
_ncm_spline2d_bicubic_prepare (NcmSpline2d *s2d)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  gsize i, j;

  s2dbc->optimize_dx.init = FALSE;
  s2dbc->optimize_dy.init = FALSE;

  g_assert (s2d->no_stride);

  /* First calculate z and dzdy in all knots and save in NCM_SPLINE2D_BICUBIC_COEFF (s2d, i, j) */
  /* First column */

  j = 0;
  {
    NcmSplineCubic *sc  = NCM_SPLINE_CUBIC (s2dbc->z_y);
    NcmVector *yv       = ncm_spline_get_yv (s2dbc->z_y);
    NcmVector *zm_col_j = ncm_matrix_get_col (s2d->zm, j);

    ncm_spline_set_yv (s2dbc->z_y, zm_col_j, FALSE);
    ncm_spline_prepare_base (s2dbc->z_y);
    ncm_vector_free (zm_col_j);

    i = 0;
    {
      const gdouble b_i = ncm_spline2d_bicubic_bi (sc, s2dbc->z_y->xv, s2dbc->z_y->yv, i);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = b_i;
    }

    for (i = 1; i < ncm_matrix_col_len (s2d->zm) - 2; i++)
    {
      const gdouble b_i = ncm_spline2d_bicubic_bi (sc, s2dbc->z_y->xv, s2dbc->z_y->yv, i);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = b_i;
    }

    i = ncm_matrix_col_len (s2d->zm) - 2;
    {
      gdouble b_i, b_ip1;

      ncm_spline2d_bicubic_bi_bip1 (sc, s2dbc->z_y->xv, s2dbc->z_y->yv, i, &b_i, &b_ip1);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i + 1);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = b_ip1;
    }

    ncm_spline_set_yv (s2dbc->z_y, yv, FALSE);
    ncm_vector_free (yv);
  }

  /* Second to the one before the one before last */
  for (j = 1; j < ncm_matrix_row_len (s2d->zm); j++)
  {
    NcmSplineCubic *sc  = NCM_SPLINE_CUBIC (s2dbc->z_y);
    NcmVector *zm_col_j = ncm_matrix_get_col (s2d->zm, j);

    ncm_spline_set_yv (s2dbc->z_y, zm_col_j, FALSE);
    ncm_spline_prepare_base (s2dbc->z_y);
    ncm_vector_free (zm_col_j);

    i = 0;
    {
      const gdouble b_i = ncm_spline2d_bicubic_bi (sc, s2dbc->z_y->xv, s2dbc->z_y->yv, i);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FY] = b_i;
    }

    for (i = 1; i < ncm_matrix_col_len (s2d->zm) - 2; i++)
    {
      const gdouble b_i = ncm_spline2d_bicubic_bi (sc, s2dbc->z_y->xv, s2dbc->z_y->yv, i);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FY] = b_i;
    }

    i = ncm_matrix_col_len (s2d->zm) - 2;
    {
      gdouble b_i, b_ip1;

      ncm_spline2d_bicubic_bi_bip1 (sc, s2dbc->z_y->xv, s2dbc->z_y->yv, i, &b_i, &b_ip1);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i + 1);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY] = b_ip1;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FY] = b_i;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_F]  = ncm_vector_get (s2dbc->z_y->yv, i + 1);
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FY] = b_ip1;
    }
  }

  /* Second calculate dzdx and d2zdxdy in all knots and save in NCM_SPLINE2D_BICUBIC_COEFF (s2d, i, j) */
  /* First row */
  i = 0;
  {
    gdouble *d          = &NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, 0)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY];
    NcmVector *dzdyv    = ncm_vector_new_data_static (d, ncm_matrix_row_len (s2d->zm), 16);
    NcmSplineCubic *sc  = NCM_SPLINE_CUBIC (s2dbc->z_x);
    NcmSplineCubic *dsc = NCM_SPLINE_CUBIC (s2dbc->dzdy_x);
    NcmVector *zm_row_i = ncm_matrix_get_row (s2d->zm, i);

    ncm_spline_set_yv (s2dbc->z_x, zm_row_i, FALSE);
    ncm_spline_set_yv (s2dbc->dzdy_x, dzdyv, FALSE);
    ncm_spline_prepare_base (s2dbc->z_x);
    ncm_spline_prepare_base (s2dbc->dzdy_x);
    ncm_vector_free (zm_row_i);
    ncm_vector_free (dzdyv);

    j = 0;
    {
      const gdouble b_j  = ncm_spline2d_bicubic_bi (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j);
      const gdouble db_j = ncm_spline2d_bicubic_bi (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX]  = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }

    for (j = 1; j < ncm_matrix_row_len (s2d->zm) - 1; j++)
    {
      const gdouble b_j  = ncm_spline2d_bicubic_bi (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j);
      const gdouble db_j = ncm_spline2d_bicubic_bi (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }

    j = ncm_matrix_row_len (s2d->zm) - 1;
    {
      gdouble b_jm1, b_j;
      gdouble db_jm1, db_j;

      ncm_spline2d_bicubic_bi_bip1 (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j - 1, &b_jm1, &b_j);
      ncm_spline2d_bicubic_bi_bip1 (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j - 1, &db_jm1, &db_j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }
  }

  /* Second row to the one before last */
  for (i = 1; i < ncm_matrix_col_len (s2d->zm) - 1; i++)
  {
    gdouble *d          = &NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, 0)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FY];
    NcmVector *dzdyv    = ncm_vector_new_data_static (d, ncm_matrix_row_len (s2d->zm), 16);
    NcmSplineCubic *sc  = NCM_SPLINE_CUBIC (s2dbc->z_x);
    NcmSplineCubic *dsc = NCM_SPLINE_CUBIC (s2dbc->dzdy_x);
    NcmVector *zm_row_i = ncm_matrix_get_row (s2d->zm, i);

    ncm_spline_set_yv (s2dbc->z_x, zm_row_i, FALSE);
    ncm_spline_set_yv (s2dbc->dzdy_x, dzdyv, FALSE);
    ncm_spline_prepare_base (s2dbc->z_x);
    ncm_spline_prepare_base (s2dbc->dzdy_x);
    ncm_vector_free (zm_row_i);
    ncm_vector_free (dzdyv);

    j = 0;
    {
      const gdouble b_j  = ncm_spline2d_bicubic_bi (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j);
      const gdouble db_j = ncm_spline2d_bicubic_bi (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX]  = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY] = db_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX]  = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }

    for (j = 1; j < ncm_matrix_row_len (s2d->zm) - 1; j++)
    {
      const gdouble b_j  = ncm_spline2d_bicubic_bi (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j);
      const gdouble db_j = ncm_spline2d_bicubic_bi (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FXY] = db_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }

    j = ncm_matrix_row_len (s2d->zm) - 1;
    {
      gdouble b_jm1, b_j;
      gdouble db_jm1, db_j;

      ncm_spline2d_bicubic_bi_bip1 (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j - 1, &b_jm1, &b_j);
      ncm_spline2d_bicubic_bi_bip1 (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j - 1, &db_jm1, &db_j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j)[NCM_SPLINE2D_BICUBIC_00][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i,     j - 1)[NCM_SPLINE2D_BICUBIC_10][NCM_SPLINE2D_BICUBIC_FXY] = db_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }
  }

  /* The last row */
  i = ncm_matrix_col_len (s2d->zm) - 1;
  {
    gdouble *d          = &NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, 0)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FY];
    NcmVector *dzdyv    = ncm_vector_new_data_static (d, ncm_matrix_row_len (s2d->zm), 16);
    NcmSplineCubic *sc  = NCM_SPLINE_CUBIC (s2dbc->z_x);
    NcmSplineCubic *dsc = NCM_SPLINE_CUBIC (s2dbc->dzdy_x);
    NcmVector *zm_row_i = ncm_matrix_get_row (s2d->zm, i);

    ncm_spline_set_yv (s2dbc->z_x, zm_row_i, FALSE);
    ncm_spline_set_yv (s2dbc->dzdy_x, dzdyv, FALSE);
    ncm_spline_prepare_base (s2dbc->z_x);
    ncm_spline_prepare_base (s2dbc->dzdy_x);
    ncm_vector_free (zm_row_i);
    ncm_vector_free (dzdyv);

    j = 0;
    {
      const gdouble b_j  = ncm_spline2d_bicubic_bi (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j);
      const gdouble db_j = ncm_spline2d_bicubic_bi (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX]  = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }

    for (j = 1; j < ncm_matrix_row_len (s2d->zm) - 1; j++)
    {
      const gdouble b_j  = ncm_spline2d_bicubic_bi (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j);
      const gdouble db_j = ncm_spline2d_bicubic_bi (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }

    j = ncm_matrix_row_len (s2d->zm) - 1;
    {
      gdouble b_jm1, b_j;
      gdouble db_jm1, db_j;

      ncm_spline2d_bicubic_bi_bip1 (sc, s2dbc->z_x->xv, s2dbc->z_x->yv, j - 1, &b_jm1, &b_j);
      ncm_spline2d_bicubic_bi_bip1 (dsc, s2dbc->dzdy_x->xv, s2dbc->dzdy_x->yv, j - 1, &db_jm1, &db_j);

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FX]     = b_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FX] = b_j;

      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j)[NCM_SPLINE2D_BICUBIC_01][NCM_SPLINE2D_BICUBIC_FXY]     = db_j;
      NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i - 1, j - 1)[NCM_SPLINE2D_BICUBIC_11][NCM_SPLINE2D_BICUBIC_FXY] = db_j;
    }
  }

  for (j = 0; j < ncm_matrix_row_len (s2d->zm) - 1; j++)
  {
    const gdouble dx = ncm_vector_get (s2d->xv, j + 1) - ncm_vector_get (s2d->xv, j);

    for (i = 0; i < ncm_matrix_col_len (s2d->zm) - 1; i++)
    {
      const gdouble dy = ncm_vector_get (s2d->yv, i + 1) - ncm_vector_get (s2d->yv, i);

      /*printf ("Coeffs of % 20.15g %20.15g | [%zd %zd]\n", ncm_vector_get (s2d->xv, j), ncm_vector_get (s2d->yv, i), i, j); */
      ncm_spline2d_bicubic_fij_to_aij (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), dx, dy, &NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j));
    }
  }

  s2d->init = TRUE;
}

static gdouble
_ncm_spline2d_bicubic_eval (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (!s2d->use_acc)
  {
    const gsize j    = gsl_interp_bsearch (s2d->x_data, x, 0, s2d->x_interv);
    const gsize i    = gsl_interp_bsearch (s2d->y_data, y, 0, s2d->y_interv);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
  else
  {
    const gsize j    = gsl_interp_accel_find (s2d->acc_x, s2d->x_data, s2d->x_interv + 1, x);
    const gsize i    = gsl_interp_accel_find (s2d->acc_y, s2d->y_data, s2d->y_interv + 1, y);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
}

static gdouble
_ncm_spline2d_bicubic_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (!s2d->use_acc)
  {
    const gsize j    = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
    const gsize i    = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
    const gdouble x0 = ncm_vector_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_dzdx (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
  else
  {
    const gsize j    = gsl_interp_accel_find (s2d->acc_x, ncm_vector_ptr (s2d->xv, 0), ncm_vector_len (s2d->xv), x);
    const gsize i    = gsl_interp_accel_find (s2d->acc_y, ncm_vector_ptr (s2d->yv, 0), ncm_vector_len (s2d->yv), y);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_dzdx (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
}

static gdouble
_ncm_spline2d_bicubic_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (!s2d->use_acc)
  {
    const gsize j    = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
    const gsize i    = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
    const gdouble x0 = ncm_vector_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_dzdy (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
  else
  {
    const gsize j    = gsl_interp_accel_find (s2d->acc_x, ncm_vector_ptr (s2d->xv, 0), ncm_vector_len (s2d->xv), x);
    const gsize i    = gsl_interp_accel_find (s2d->acc_y, ncm_vector_ptr (s2d->yv, 0), ncm_vector_len (s2d->yv), y);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_dzdy (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
}

static gdouble
_ncm_spline2d_bicubic_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (!s2d->use_acc)
  {
    const gsize j    = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
    const gsize i    = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
    const gdouble x0 = ncm_vector_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_d2zdx2 (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
  else
  {
    const gsize j    = gsl_interp_accel_find (s2d->acc_x, ncm_vector_ptr (s2d->xv, 0), ncm_vector_len (s2d->xv), x);
    const gsize i    = gsl_interp_accel_find (s2d->acc_y, ncm_vector_ptr (s2d->yv, 0), ncm_vector_len (s2d->yv), y);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_d2zdx2 (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
}

static gdouble
_ncm_spline2d_bicubic_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (!s2d->use_acc)
  {
    const gsize j    = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
    const gsize i    = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
    const gdouble x0 = ncm_vector_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_d2zdy2 (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
  else
  {
    const gsize j    = gsl_interp_accel_find (s2d->acc_x, ncm_vector_ptr (s2d->xv, 0), ncm_vector_len (s2d->xv), x);
    const gsize i    = gsl_interp_accel_find (s2d->acc_y, ncm_vector_ptr (s2d->yv, 0), ncm_vector_len (s2d->yv), y);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_d2zdy2 (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
}

static gdouble
_ncm_spline2d_bicubic_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  if (!s2d->use_acc)
  {
    const gsize j    = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
    const gsize i    = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
    const gdouble x0 = ncm_vector_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_d2zdxy (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
  else
  {
    const gsize j    = gsl_interp_accel_find (s2d->acc_x, ncm_vector_ptr (s2d->xv, 0), ncm_vector_len (s2d->xv), x);
    const gsize i    = gsl_interp_accel_find (s2d->acc_y, ncm_vector_ptr (s2d->yv, 0), ncm_vector_len (s2d->yv), y);
    const gdouble x0 = ncm_vector_fast_get (s2d->xv, j);
    const gdouble y0 = ncm_vector_fast_get (s2d->yv, i);

    return ncm_spline2d_bicubic_eval_poly_d2zdxy (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j), x - x0, y - y0);
  }
}

static gdouble
_ncm_spline2d_bicubic_int_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  gsize jl = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xl, 0, ncm_vector_len (s2d->xv) - 1);
  gsize ju = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xu, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), y, 0, ncm_vector_len (s2d->yv) - 1);
  gdouble x0, x1, result;
  const gdouble y0 = ncm_vector_get (s2d->yv, i);
  guint k;
  gdouble coeffs[4];

  g_assert (jl <= ju);

  if (jl == ju)
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    /* x1 = ncm_vector_get (s2d->xv, jl + 1); */
    ncm_spline2d_bicubic_integ_dx_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, jl), y - y0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], x0, xl, xu);
  }
  else
  {
    x0 = ncm_vector_get (s2d->xv, jl);
    x1 = ncm_vector_get (s2d->xv, jl + 1);
    ncm_spline2d_bicubic_integ_dx_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, jl), y - y0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], x0, xl, x1);

    for (k = jl + 1; k < ju; k++)
    {
      x0 = ncm_vector_get (s2d->xv, k);
      x1 = ncm_vector_get (s2d->xv, k + 1);
      ncm_spline2d_bicubic_integ_dx_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, k), y - y0, coeffs);

      {
        const gdouble dx = x1 - x0;

        result += dx * (coeffs[0] + dx * (0.5 * coeffs[1] + dx * (coeffs[2] / 3.0 + 0.25 * coeffs[3] * dx)));
      }
    }

    k = ju;
    {
      x0 = ncm_vector_get (s2d->xv, ju);
      /* x1 = ncm_vector_get (s2d->xv, ju + 1); */
      ncm_spline2d_bicubic_integ_dx_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, k), y - y0, coeffs);
      result += _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], x0, x0, xu);
    }
  }

  return result;
}

static NcmSpline *
_ncm_spline2d_bicubic_int_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  gsize jl = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xl, 0, ncm_vector_len (s2d->xv) - 1);
  gsize ju = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xu, 0, ncm_vector_len (s2d->xv) - 1);
  gsize i, j;
  gdouble x0, x1;
  NcmSplineCubic *sc;
  guint y_len_m1 = ncm_vector_len (s2d->yv) - 1;

  if ((s2dbc->optimize_dx.init && (s2dbc->optimize_dx.l == xl) && (s2dbc->optimize_dx.u == xu)))
    return s2dbc->optimize_dx.s;

  sc = NCM_SPLINE_CUBIC (s2dbc->optimize_dx.s);

#define _NC_AIJ NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)
#define _NCM_INTEGRAL_A (s2dbc->optimize_dx.s->yv)
#define _NCM_INTEGRAL_B (sc->b)
#define _NCM_INTEGRAL_C (sc->c)
#define _NCM_INTEGRAL_D (sc->d)

  g_assert (jl <= ju);

  s2dbc->optimize_dx.l = xl;
  s2dbc->optimize_dx.u = xu;

  if (jl == ju)
  {
    j  = jl;
    x0 = ncm_vector_get (s2d->xv, jl);
    /* x1 = ncm_vector_get (s2d->xv, jl + 1); */

    for (i = 0; i < y_len_m1; i++)
    {
      ncm_vector_set (_NCM_INTEGRAL_A, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][0], _NC_AIJ[1][0], _NC_AIJ[2][0], _NC_AIJ[3][0], x0, xl, xu));
      ncm_vector_set (_NCM_INTEGRAL_B, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][1], _NC_AIJ[1][1], _NC_AIJ[2][1], _NC_AIJ[3][1], x0, xl, xu));
      ncm_vector_set (_NCM_INTEGRAL_C, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][2], _NC_AIJ[1][2], _NC_AIJ[2][2], _NC_AIJ[3][2], x0, xl, xu));
      ncm_vector_set (_NCM_INTEGRAL_D, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][3], _NC_AIJ[1][3], _NC_AIJ[2][3], _NC_AIJ[3][3], x0, xl, xu));
    }
  }
  else
  {
    for (i = 0; i < y_len_m1; i++)
    {
      j  = jl;
      x0 = ncm_vector_get (s2d->xv, jl);
      x1 = ncm_vector_get (s2d->xv, jl + 1);

      ncm_vector_set (_NCM_INTEGRAL_A, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][0], _NC_AIJ[1][0], _NC_AIJ[2][0], _NC_AIJ[3][0], x0, xl, x1));
      ncm_vector_set (_NCM_INTEGRAL_B, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][1], _NC_AIJ[1][1], _NC_AIJ[2][1], _NC_AIJ[3][1], x0, xl, x1));
      ncm_vector_set (_NCM_INTEGRAL_C, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][2], _NC_AIJ[1][2], _NC_AIJ[2][2], _NC_AIJ[3][2], x0, xl, x1));
      ncm_vector_set (_NCM_INTEGRAL_D, i,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][3], _NC_AIJ[1][3], _NC_AIJ[2][3], _NC_AIJ[3][3], x0, xl, x1));

      for (j = jl + 1; j < ju; j++)
      {
        x0 = ncm_vector_get (s2d->xv, j);
        x1 = ncm_vector_get (s2d->xv, j + 1);
        {
          const gdouble dx = x1 - x0;

          ncm_vector_addto (_NCM_INTEGRAL_A, i,
                            dx * (_NC_AIJ[0][0] + dx * (0.5 * _NC_AIJ[1][0] + dx * (_NC_AIJ[2][0] / 3.0 + 0.25 * _NC_AIJ[3][0] * dx))));
          ncm_vector_addto (_NCM_INTEGRAL_B, i,
                            dx * (_NC_AIJ[0][1] + dx * (0.5 * _NC_AIJ[1][1] + dx * (_NC_AIJ[2][1] / 3.0 + 0.25 * _NC_AIJ[3][1] * dx))));
          ncm_vector_addto (_NCM_INTEGRAL_C, i,
                            dx * (_NC_AIJ[0][2] + dx * (0.5 * _NC_AIJ[1][2] + dx * (_NC_AIJ[2][2] / 3.0 + 0.25 * _NC_AIJ[3][2] * dx))));
          ncm_vector_addto (_NCM_INTEGRAL_D, i,
                            dx * (_NC_AIJ[0][3] + dx * (0.5 * _NC_AIJ[1][3] + dx * (_NC_AIJ[2][3] / 3.0 + 0.25 * _NC_AIJ[3][3] * dx))));
        }
      }

      j = ju;
      {
        x0 = ncm_vector_get (s2d->xv, ju);
        /* x1 = ncm_vector_get (s2d->xv, ju + 1); */
        ncm_vector_addto (_NCM_INTEGRAL_A, i,
                          _ncm_spline_util_integ_eval (_NC_AIJ[0][0], _NC_AIJ[1][0], _NC_AIJ[2][0], _NC_AIJ[3][0], x0, x0, xu));
        ncm_vector_addto (_NCM_INTEGRAL_B, i,
                          _ncm_spline_util_integ_eval (_NC_AIJ[0][1], _NC_AIJ[1][1], _NC_AIJ[2][1], _NC_AIJ[3][1], x0, x0, xu));
        ncm_vector_addto (_NCM_INTEGRAL_C, i,
                          _ncm_spline_util_integ_eval (_NC_AIJ[0][2], _NC_AIJ[1][2], _NC_AIJ[2][2], _NC_AIJ[3][2], x0, x0, xu));
        ncm_vector_addto (_NCM_INTEGRAL_D, i,
                          _ncm_spline_util_integ_eval (_NC_AIJ[0][3], _NC_AIJ[1][3], _NC_AIJ[2][3], _NC_AIJ[3][3], x0, x0, xu));
      }
    }
  }

  s2dbc->optimize_dx.init = TRUE;
#undef _NC_AIJ
#undef _NCM_INTEGRAL_A
#undef _NCM_INTEGRAL_B
#undef _NCM_INTEGRAL_C
#undef _NCM_INTEGRAL_D

  return s2dbc->optimize_dx.s;
}

static gdouble
_ncm_spline2d_bicubic_int_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  gsize j = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), x, 0, ncm_vector_len (s2d->xv) - 1);
  gsize il = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yl, 0, ncm_vector_len (s2d->yv) - 1);
  gsize iu = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yu, 0, ncm_vector_len (s2d->yv) - 1);
  gdouble y0, y1, result;
  const gdouble x0 = ncm_vector_get (s2d->xv, j);
  guint k;
  gdouble coeffs[4];

  g_assert (il <= iu);

  if (il == iu)
  {
    y0 = ncm_vector_get (s2d->yv, il);
    /* y1 = ncm_vector_get (s2d->yv, il + 1); */
    ncm_spline2d_bicubic_integ_dy_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, j), x - x0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], y0, yl, yu);
  }
  else
  {
    y0 = ncm_vector_get (s2d->yv, il);
    y1 = ncm_vector_get (s2d->yv, il + 1);
    ncm_spline2d_bicubic_integ_dy_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, j), x - x0, coeffs);
    result = _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], y0, yl, y1);

    for (k = il + 1; k < iu; k++)
    {
      y0 = ncm_vector_get (s2d->yv, k);
      y1 = ncm_vector_get (s2d->yv, k + 1);
      ncm_spline2d_bicubic_integ_dy_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, j), x - x0, coeffs);

      {
        const gdouble dy = y1 - y0;

        result += dy * (coeffs[0] + dy * (0.5 * coeffs[1] + dy * (coeffs[2] / 3.0 + 0.25 * coeffs[3] * dy)));
      }
    }

    k = iu;
    {
      y0 = ncm_vector_get (s2d->yv, k);
      /* y1 = ncm_vector_get (s2d->yv, k + 1); */
      ncm_spline2d_bicubic_integ_dy_coeffs (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, j), x - x0, coeffs);
      result += _ncm_spline_util_integ_eval (coeffs[0], coeffs[1], coeffs[2], coeffs[3], y0, y0, yu);
    }
  }

  return result;
}

static NcmSpline *
_ncm_spline2d_bicubic_int_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  gsize il = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yl, 0, ncm_vector_len (s2d->yv) - 1);
  gsize iu = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yu, 0, ncm_vector_len (s2d->yv) - 1);
  gsize i, j;
  gdouble y0, y1;
  NcmSplineCubic *sc;
  guint x_len_m1 = ncm_vector_len (s2d->xv) - 1;

  if ((s2dbc->optimize_dy.init && (s2dbc->optimize_dy.l == yl) && (s2dbc->optimize_dy.u == yu)))
    return s2dbc->optimize_dy.s;

  sc = NCM_SPLINE_CUBIC (s2dbc->optimize_dy.s);

#define _NC_AIJ NCM_SPLINE2D_BICUBIC_COEFF (s2dbc, i, j)
#define _NCM_INTEGRAL_A (s2dbc->optimize_dy.s->yv)
#define _NCM_INTEGRAL_B (sc->b)
#define _NCM_INTEGRAL_C (sc->c)
#define _NCM_INTEGRAL_D (sc->d)

  g_assert (il <= iu);

  s2dbc->optimize_dy.l = yl;
  s2dbc->optimize_dy.u = yu;

  if (il == iu)
  {
    i  = il;
    y0 = ncm_vector_get (s2d->yv, il);
    /* y1 = ncm_vector_get (s2d->yv, il + 1); */

    for (j = 0; j < x_len_m1; j++)
    {
      ncm_vector_set (_NCM_INTEGRAL_A, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][0], _NC_AIJ[0][1], _NC_AIJ[0][2], _NC_AIJ[0][3], y0, yl, yu));
      ncm_vector_set (_NCM_INTEGRAL_B, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[1][0], _NC_AIJ[1][1], _NC_AIJ[1][2], _NC_AIJ[1][3], y0, yl, yu));
      ncm_vector_set (_NCM_INTEGRAL_C, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[2][0], _NC_AIJ[2][1], _NC_AIJ[2][2], _NC_AIJ[2][3], y0, yl, yu));
      ncm_vector_set (_NCM_INTEGRAL_D, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[3][0], _NC_AIJ[3][1], _NC_AIJ[3][2], _NC_AIJ[3][3], y0, yl, yu));
    }
  }
  else
  {
    for (j = 0; j < x_len_m1; j++)
    {
      i  = il;
      y0 = ncm_vector_get (s2d->yv, il);
      y1 = ncm_vector_get (s2d->yv, il + 1);

      ncm_vector_set (_NCM_INTEGRAL_A, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[0][0], _NC_AIJ[0][1], _NC_AIJ[0][2], _NC_AIJ[0][3], y0, yl, y1));
      ncm_vector_set (_NCM_INTEGRAL_B, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[1][0], _NC_AIJ[1][1], _NC_AIJ[1][2], _NC_AIJ[1][3], y0, yl, y1));
      ncm_vector_set (_NCM_INTEGRAL_C, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[2][0], _NC_AIJ[2][1], _NC_AIJ[2][2], _NC_AIJ[2][3], y0, yl, y1));
      ncm_vector_set (_NCM_INTEGRAL_D, j,
                      _ncm_spline_util_integ_eval (_NC_AIJ[3][0], _NC_AIJ[3][1], _NC_AIJ[3][2], _NC_AIJ[3][3], y0, yl, y1));

      for (i = il + 1; i < iu; i++)
      {
        y0 = ncm_vector_get (s2d->yv, i);
        y1 = ncm_vector_get (s2d->yv, i + 1);

        {
          const gdouble dy = y1 - y0;

          ncm_vector_addto (_NCM_INTEGRAL_A, j,
                            dy * (_NC_AIJ[0][0] + dy * (0.5 * _NC_AIJ[0][1] + dy * (_NC_AIJ[0][2] / 3.0 + 0.25 * _NC_AIJ[0][3] * dy))));
          ncm_vector_addto (_NCM_INTEGRAL_B, j,
                            dy * (_NC_AIJ[1][0] + dy * (0.5 * _NC_AIJ[1][1] + dy * (_NC_AIJ[1][2] / 3.0 + 0.25 * _NC_AIJ[1][3] * dy))));
          ncm_vector_addto (_NCM_INTEGRAL_C, j,
                            dy * (_NC_AIJ[2][0] + dy * (0.5 * _NC_AIJ[2][1] + dy * (_NC_AIJ[2][2] / 3.0 + 0.25 * _NC_AIJ[2][3] * dy))));
          ncm_vector_addto (_NCM_INTEGRAL_D, j,
                            dy * (_NC_AIJ[3][0] + dy * (0.5 * _NC_AIJ[3][1] + dy * (_NC_AIJ[3][2] / 3.0 + 0.25 * _NC_AIJ[3][3] * dy))));
        }
      }

      i = iu;
      {
        y0 = ncm_vector_get (s2d->yv, iu);
        /* y1 = ncm_vector_get (s2d->yv, iu + 1); */
        ncm_vector_addto (_NCM_INTEGRAL_A, j,
                          _ncm_spline_util_integ_eval (_NC_AIJ[0][0], _NC_AIJ[0][1], _NC_AIJ[0][2], _NC_AIJ[0][3], y0, y0, yu));
        ncm_vector_addto (_NCM_INTEGRAL_B, j,
                          _ncm_spline_util_integ_eval (_NC_AIJ[1][0], _NC_AIJ[1][1], _NC_AIJ[1][2], _NC_AIJ[1][3], y0, y0, yu));
        ncm_vector_addto (_NCM_INTEGRAL_C, j,
                          _ncm_spline_util_integ_eval (_NC_AIJ[2][0], _NC_AIJ[2][1], _NC_AIJ[2][2], _NC_AIJ[2][3], y0, y0, yu));
        ncm_vector_addto (_NCM_INTEGRAL_D, j,
                          _ncm_spline_util_integ_eval (_NC_AIJ[3][0], _NC_AIJ[3][1], _NC_AIJ[3][2], _NC_AIJ[3][3], y0, y0, yu));
      }
    }
  }

  s2dbc->optimize_dy.init = TRUE;
#undef _NC_AIJ
#undef _NCM_INTEGRAL_A
#undef _NCM_INTEGRAL_B
#undef _NCM_INTEGRAL_C
#undef _NCM_INTEGRAL_D

  return s2dbc->optimize_dy.s;
}

static gdouble
_ncm_spline2d_bicubic_int_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);

  gsize jl = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xl, 0, ncm_vector_len (s2d->xv) - 1);
  gsize ju = gsl_interp_bsearch (ncm_vector_ptr (s2d->xv, 0), xu, 0, ncm_vector_len (s2d->xv) - 1);
  gsize il = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yl, 0, ncm_vector_len (s2d->yv) - 1);
  gsize iu = gsl_interp_bsearch (ncm_vector_ptr (s2d->yv, 0), yu, 0, ncm_vector_len (s2d->yv) - 1);
  gdouble x0, x1, y0, y1, result;
  guint k, m;

  g_assert (jl <= ju || il <= iu);

  if ((jl == ju) && (il == iu))
  {
    x0     = ncm_vector_get (s2d->xv, jl);
    y0     = ncm_vector_get (s2d->yv, il);
    result = ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, jl), x0, xl, xu, y0, yl, yu);
  }
  else if (jl == ju)
  {
    x0     = ncm_vector_get (s2d->xv, jl);
    y0     = ncm_vector_get (s2d->yv, il);
    y1     = ncm_vector_get (s2d->yv, il + 1);
    result = ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, jl), x0, xl, xu, y0, yl, y1);

    for (k = il + 1; k < iu; k++)
    {
      y0      = ncm_vector_get (s2d->yv, k);
      y1      = ncm_vector_get (s2d->yv, k + 1);
      result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, jl), x0, xl, xu, y0, y0, y1);
    }

    k = iu;
    {
      y0      = ncm_vector_get (s2d->yv, k);
      result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, jl), x0, xl, xu, y0, y0, yu);
    }
  }
  else if (il == iu)
  {
    x0     = ncm_vector_get (s2d->xv, jl);
    x1     = ncm_vector_get (s2d->xv, jl + 1);
    y0     = ncm_vector_get (s2d->yv, il);
    result = ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, jl), x0, xl, x1, y0, yl, yu);

    for (k = jl + 1; k < ju; k++)
    {
      x0      = ncm_vector_get (s2d->xv, k);
      x1      = ncm_vector_get (s2d->xv, k + 1);
      result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, k), x0, x0, x1, y0, yl, yu);
    }

    k = ju;
    {
      x0      = ncm_vector_get (s2d->xv, k);
      result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, k), x0, x0, xu, y0, yl, yu);
    }
  }
  else
  {
    m = jl;
    {
      x0     = ncm_vector_get (s2d->xv, jl);
      x1     = ncm_vector_get (s2d->xv, jl + 1);
      y0     = ncm_vector_get (s2d->yv, il);
      y1     = ncm_vector_get (s2d->yv, il + 1);
      result = ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, m), x0, xl, x1, y0, yl, y1);

      for (k = il + 1; k < iu; k++)
      {
        y0      = ncm_vector_get (s2d->yv, k);
        y1      = ncm_vector_get (s2d->yv, k + 1);
        result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, m), x0, xl, x1, y0, y0, y1);
      }

      k = iu;
      {
        y0      = ncm_vector_get (s2d->yv, k);
        result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, m), x0, xl, x1, y0, y0, yu);
      }
    }

    for (m = jl + 1; m < ju; m++)
    {
      x0      = ncm_vector_get (s2d->xv, m);
      x1      = ncm_vector_get (s2d->xv, m + 1);
      y0      = ncm_vector_get (s2d->yv, il);
      y1      = ncm_vector_get (s2d->yv, il + 1);
      result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, m), x0, x0, x1, y0, yl, y1);

      for (k = il + 1; k < iu; k++)
      {
        y0      = ncm_vector_get (s2d->yv, k);
        y1      = ncm_vector_get (s2d->yv, k + 1);
        result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, m), x0, x0, x1, y0, y0, y1);
      }

      k = iu;
      {
        y0      = ncm_vector_get (s2d->yv, k);
        result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, m), x0, x0, x1, y0, y0, yu);
      }
    }

    m = ju;
    {
      x0      = ncm_vector_get (s2d->xv, m);
      y0      = ncm_vector_get (s2d->yv, il);
      y1      = ncm_vector_get (s2d->yv, il + 1);
      result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, il, m), x0, x0, xu, y0, yl, y1);

      for (k = il + 1; k < iu; k++)
      {
        y0      = ncm_vector_get (s2d->yv, k);
        y1      = ncm_vector_get (s2d->yv, k + 1);
        result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, m), x0, x0, xu, y0, y0, y1);
      }

      k = iu;
      {
        y0      = ncm_vector_get (s2d->yv, k);
        result += ncm_spline2d_bicubic_integ_eval2d (&NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, k, m), x0, x0, xu, y0, y0, yu);
      }
    }
  }

  return result;
}

static void
_ncm_spline2d_bicubic_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res)
{
  NcmSpline2dBicubic *s2dbc = NCM_SPLINE2D_BICUBIC (s2d);
  const gsize j             = gsl_interp_bsearch (s2d->x_data, x, 0, s2d->x_interv);
  const gdouble x0          = ncm_vector_fast_get (s2d->xv, j);
  const gdouble dx          = x - x0;
  const guint len           = ncm_vector_len (y);
  NcmSpline2dBicubicCoeffs *sa;
  gint k, i;
  size_t l;
  gdouble y0, y_l, a0, a1, a2, a3;

  g_assert_cmpuint (len, ==, order->len);
  g_assert_cmpuint (len, ==, res->len);

  l   = g_array_index (order, size_t, 0);
  y_l = ncm_vector_get (y, l);
  i   = gsl_interp_bsearch (s2d->y_data, y_l, 0, s2d->y_interv);
  y0  = ncm_vector_fast_get (s2d->yv, i);
  sa  = &NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j);

  a0 = (sa->ij[0][0] + dx * (sa->ij[1][0] + dx * (sa->ij[2][0] + dx * sa->ij[3][0])));
  a1 = (sa->ij[0][1] + dx * (sa->ij[1][1] + dx * (sa->ij[2][1] + dx * sa->ij[3][1])));
  a2 = (sa->ij[0][2] + dx * (sa->ij[1][2] + dx * (sa->ij[2][2] + dx * sa->ij[3][2])));
  a3 = (sa->ij[0][3] + dx * (sa->ij[1][3] + dx * (sa->ij[2][3] + dx * sa->ij[3][3])));

  for (k = 0; k < len; k++)
  {
    l   = g_array_index (order, size_t, k);
    y_l = ncm_vector_get (y, l);

    if (y_l > ncm_vector_fast_get (s2d->yv, i + 1))
    {
      do {
        i++;
      } while (y_l > ncm_vector_fast_get (s2d->yv, i + 1));

      y0 = ncm_vector_fast_get (s2d->yv, i);
      sa = &NCM_SPLINE2D_BICUBIC_STRUCT (s2dbc, i, j);

      a0 = (sa->ij[0][0] + dx * (sa->ij[1][0] + dx * (sa->ij[2][0] + dx * sa->ij[3][0])));
      a1 = (sa->ij[0][1] + dx * (sa->ij[1][1] + dx * (sa->ij[2][1] + dx * sa->ij[3][1])));
      a2 = (sa->ij[0][2] + dx * (sa->ij[1][2] + dx * (sa->ij[2][2] + dx * sa->ij[3][2])));
      a3 = (sa->ij[0][3] + dx * (sa->ij[1][3] + dx * (sa->ij[2][3] + dx * sa->ij[3][3])));
    }

    {
      const gdouble dy    = y_l - y0;
      const gdouble res_l = a0 + dy * (a1 + dy * (a2 + dy * a3));

      g_array_index (res, gdouble, l) = res_l;
    }
  }
}


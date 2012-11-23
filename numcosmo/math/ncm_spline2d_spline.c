/***************************************************************************
 *            ncm_spline2d_spline.c
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
 * SECTION:ncm_spline2d_spline
 * @title: Two Dimensional Spline from Spline
 * @short_description: Abstract class for two dimensional splines.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline2d_spline.h"

G_DEFINE_TYPE (NcmSpline2dSpline, ncm_spline2d_spline, NCM_TYPE_SPLINE2D);

/**
 * ncm_spline2d_spline_new:
 * @s: a #NcmSpline.
 *
 * This function initializes a #NcmSpline2d
 * FIXME
 *
 * Returns: A new #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_spline_new (NcmSpline *s)
{
  NcmSpline2d *s2d;

  g_assert (NCM_IS_SPLINE (s));

  s2d = g_object_new (NCM_TYPE_SPLINE2D_SPLINE, "spline", s, NULL);

  return s2d;
}

NcmSpline2d *
_ncm_spline2d_spline_copy_empty (const NcmSpline2d *s2d)
{
  return ncm_spline2d_spline_new (s2d->s);
}

static void
_ncm_spline2d_spline_alloc (NcmSpline2d *s2d)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;

  s2ds->first_prepare = TRUE;

  s2ds->vertv = ncm_vector_new_sunk (NCM_MATRIX_NROWS (s2d->zm));
  s2ds->vertintv = ncm_vector_new_sunk (NCM_MATRIX_NROWS (s2d->zm));

  s2ds->s_hor = g_slice_alloc (sizeof (NcmSpline *) * NCM_MATRIX_NROWS (s2d->zm));
  for (i = 0; i < NCM_MATRIX_NROWS (s2d->zm); i++)
    s2ds->s_hor[i] = ncm_spline_new (s2d->s, s2d->xv, ncm_matrix_get_row (s2d->zm, i), FALSE);

  s2ds->s_ver = ncm_spline_new (s2d->s, s2d->yv, s2ds->vertv, FALSE);
  s2ds->s_ver_integ = ncm_spline_new (s2d->s, s2d->yv, s2ds->vertintv, FALSE);
}

static void
_ncm_spline2d_spline_clear (NcmSpline2d *s2d)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;

  for (i = 0; i < ncm_vector_len (s2ds->vertv); i++)
    ncm_spline_clear (&s2ds->s_hor[i]);

  ncm_spline_clear (&s2ds->s_ver);
  ncm_spline_clear (&s2ds->s_ver_integ);
}

static void
_ncm_spline2d_spline_free (NcmSpline2d *s2d)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);

  _ncm_spline2d_spline_clear (s2d);
  
  g_slice_free1 (sizeof (NcmSpline *) * NCM_MATRIX_COL_LEN (s2d->zm), s2ds->s_hor);
}

static void
_ncm_spline2d_spline_reset (NcmSpline2d *s2d)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  if (s2d->init)
  {
    if ((NCM_MATRIX_NROWS (s2d->zm) != ncm_vector_len (s2ds->vertv)) ||
        (NCM_MATRIX_NCOLS (s2d->zm) != ncm_vector_len (s2ds->s_hor[0]->yv)))
    {
      _ncm_spline2d_spline_free (s2d);
      _ncm_spline2d_spline_alloc (s2d);
    }
  }
  else
    _ncm_spline2d_spline_alloc (s2d);
}


static void
_ncm_spline2d_spline_prepare (NcmSpline2d *s2d)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;
  for (i = 0; i < NCM_MATRIX_COL_LEN (s2d->zm); i++)
	ncm_spline_prepare (s2ds->s_hor[i]);

  s2d->init = TRUE;
  s2ds->first_prepare = TRUE;
  s2ds->first_prepare_integ = TRUE;
}

static gdouble
_ncm_spline2d_spline_eval (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;
  if (!s2d->init)
	ncm_spline2d_prepare (s2d);
  if (s2ds->last_x != x || s2ds->first_prepare)
  {
	for (i = 0; i < NCM_MATRIX_COL_LEN (s2d->zm); i++)
	  ncm_vector_set (s2ds->vertv, i, ncm_spline_eval (s2ds->s_hor[i], x));

	ncm_spline_prepare (s2ds->s_ver);
	s2ds->last_x = x;
	s2ds->first_prepare = FALSE;
  }

  return ncm_spline_eval (s2ds->s_ver, y);
}

static gdouble _ncm_spline2d_spline_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y) { g_error ("spsp does not implement dzdx"); return 0.0; }
static gdouble _ncm_spline2d_spline_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y) { g_error ("spsp does not implement dzdy"); return 0.0; }
static gdouble _ncm_spline2d_spline_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y) { g_error ("spsp does not implement d2zdx2"); return 0.0; }
static gdouble _ncm_spline2d_spline_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y) { g_error ("spsp does not implement d2zdy2"); return 0.0; }
static gdouble _ncm_spline2d_spline_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y) { g_error ("spsp does not implement d2zdxy"); return 0.0; }

static gdouble
_ncm_spline2d_spline_int_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;
  if (!s2d->init)
	ncm_spline2d_prepare (s2d);

  if (s2ds->last_xl != xl || s2ds->last_xu != xu || s2ds->first_prepare_integ)
  {
	for (i = 0; i < NCM_MATRIX_COL_LEN (s2d->zm); i++)
	  ncm_vector_set (s2ds->vertintv, i, ncm_spline_eval_integ (s2ds->s_hor[i], xl, xu));
	ncm_spline_prepare (s2ds->s_ver_integ);
	s2ds->last_xl = xl;
	s2ds->last_xu = xu;
	s2ds->first_prepare_integ = FALSE;
  }

  return ncm_spline_eval (s2ds->s_ver_integ, y);
}

static gdouble
_ncm_spline2d_spline_int_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;
  if (!s2d->init)
	ncm_spline2d_prepare (s2d);
  if (s2ds->last_x != x || s2ds->first_prepare)
  {
	for (i = 0; i < NCM_MATRIX_COL_LEN (s2d->zm); i++)
	  ncm_vector_set (s2ds->vertv, i, ncm_spline_eval (s2ds->s_hor[i], x));
	ncm_spline_prepare (s2ds->s_ver);
	s2ds->last_x = x;
	s2ds->first_prepare = FALSE;
  }

  return ncm_spline_eval_integ (s2ds->s_ver, yl, yu);
}

static gdouble
_ncm_spline2d_spline_int_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);
  guint i;
  if (!s2d->init)
	ncm_spline2d_prepare (s2d);
  if (s2ds->last_xl != xl || s2ds->last_xu != xu || s2ds->first_prepare_integ)
  {
	for (i = 0; i < NCM_MATRIX_COL_LEN (s2d->zm); i++)
	  ncm_vector_set (s2ds->vertintv, i, ncm_spline_eval_integ (s2ds->s_hor[i], xl, xu));
	ncm_spline_prepare (s2ds->s_ver_integ);
	s2ds->last_xl = xl;
	s2ds->last_xu = xu;
	s2ds->first_prepare_integ = FALSE;
  }

  return ncm_spline_eval_integ (s2ds->s_ver_integ, yl, yu);
}

static NcmSpline *
_ncm_spline2d_spline_int_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu)
{
  g_assert_not_reached ();
}

static NcmSpline *
_ncm_spline2d_spline_int_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu)
{
  g_assert_not_reached ();
}

static void
ncm_spline2d_spline_init (NcmSpline2dSpline *s2ds)
{
  s2ds->first_prepare = FALSE;
  s2ds->first_prepare_integ = FALSE;
  s2ds->last_x = GSL_NAN;
  s2ds->last_xl = GSL_NAN;
  s2ds->last_xu = GSL_NAN;
  s2ds->last_yl = GSL_NAN;
  s2ds->last_yu = GSL_NAN;
  s2ds->vertv = NULL;
  s2ds->vertintv = NULL;
  s2ds->s_hor = NULL;
  s2ds->s_ver = NULL;
  s2ds->s_ver_integ = NULL;
}

static void
ncm_spline2d_spline_dispose (GObject *object)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (object);

  _ncm_spline2d_spline_clear (s2d);
  s2d->init = FALSE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_spline_parent_class)->dispose (object);
}

static void
ncm_spline2d_spline_finalize (GObject *object)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (object);
  NcmSpline2dSpline *s2ds = NCM_SPLINE2D_SPLINE (s2d);

  g_slice_free1 (sizeof (NcmSpline *) * NCM_MATRIX_COL_LEN (s2d->zm), s2ds->s_hor);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_spline_parent_class)->finalize (object);
}

static void
ncm_spline2d_spline_class_init (NcmSpline2dSplineClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmSpline2dClass* parent_class = NCM_SPLINE2D_CLASS (klass);

  parent_class->copy_empty = &_ncm_spline2d_spline_copy_empty;
  parent_class->reset = &_ncm_spline2d_spline_reset;
  parent_class->prepare = &_ncm_spline2d_spline_prepare;
  parent_class->eval = &_ncm_spline2d_spline_eval;
  parent_class->dzdx = &_ncm_spline2d_spline_dzdx;
  parent_class->dzdy = &_ncm_spline2d_spline_dzdy;
  parent_class->d2zdxy = &_ncm_spline2d_spline_d2zdxy;
  parent_class->d2zdx2 = &_ncm_spline2d_spline_d2zdx2;
  parent_class->d2zdy2 = &_ncm_spline2d_spline_d2zdy2;
  parent_class->int_dx = &_ncm_spline2d_spline_int_dx;
  parent_class->int_dy = &_ncm_spline2d_spline_int_dy;
  parent_class->int_dxdy = &_ncm_spline2d_spline_int_dxdy;
  parent_class->int_dx_spline = &_ncm_spline2d_spline_int_dx_spline;
  parent_class->int_dy_spline = &_ncm_spline2d_spline_int_dy_spline;

  object_class->dispose = ncm_spline2d_spline_dispose;
  object_class->finalize = ncm_spline2d_spline_finalize;
}

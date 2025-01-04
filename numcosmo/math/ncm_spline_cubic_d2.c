/***************************************************************************
 *            ncm_spline_cubic_d2.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
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

/**
 * NcmSplineCubicD2:
 *
 * Cubic spline implementation given second derivatives.
 *
 * This object implements the necessary functions to compute a cubic spline with where
 * the user provides the second derivatives of the function.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_d2.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_lapack.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmSplineCubicD2
{
  /*< private >*/
  NcmSplineCubic parent_instance;
  NcmVector *d2;
  gsize len;
};

G_DEFINE_TYPE (NcmSplineCubicD2, ncm_spline_cubic_d2, NCM_TYPE_SPLINE_CUBIC)

static void
ncm_spline_cubic_d2_init (NcmSplineCubicD2 *s)
{
  s->d2  = NULL;
  s->len = 0;
}

static void
_ncm_spline_cubic_d2_dispose (GObject *object)
{
  NcmSplineCubicD2 *s = NCM_SPLINE_CUBIC_D2 (object);

  ncm_vector_clear (&s->d2);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_cubic_d2_parent_class)->dispose (object);
}

static void
_ncm_spline_cubic_d2_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_cubic_d2_parent_class)->finalize (object);
}

static const gchar *_ncm_spline_cubic_d2_name (NcmSpline *s);
static void _ncm_spline_cubic_d2_reset (NcmSpline *s);
static void _ncm_spline_cubic_d2_prepare_base (NcmSpline *s);
static void _ncm_spline_cubic_d2_prepare (NcmSpline *s);
static gsize _ncm_spline_cubic_d2_min_size (const NcmSpline *s);
static NcmSpline *_ncm_spline_cubic_d2_copy_empty (const NcmSpline *s);

static void
ncm_spline_cubic_d2_class_init (NcmSplineCubicD2Class *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmSplineClass *s_class    = NCM_SPLINE_CLASS (klass);

  object_class->dispose  = &_ncm_spline_cubic_d2_dispose;
  object_class->finalize = &_ncm_spline_cubic_d2_finalize;

  s_class->reset        = &_ncm_spline_cubic_d2_reset;
  s_class->name         = &_ncm_spline_cubic_d2_name;
  s_class->prepare      = &_ncm_spline_cubic_d2_prepare;
  s_class->prepare_base = &_ncm_spline_cubic_d2_prepare_base;
  s_class->min_size     = &_ncm_spline_cubic_d2_min_size;
  s_class->copy_empty   = &_ncm_spline_cubic_d2_copy_empty;
}

static void
_ncm_spline_cubic_d2_reset (NcmSpline *s)
{
  /* Chain up : start */
  NCM_SPLINE_CLASS (ncm_spline_cubic_d2_parent_class)->reset (s);
  {
    NcmSplineCubicD2 *scd2 = NCM_SPLINE_CUBIC_D2 (s);
    const guint sc_len     = ncm_spline_get_len (s);

    if (sc_len != scd2->len)
    {
      ncm_vector_clear (&scd2->d2);

      scd2->d2  = ncm_vector_new (sc_len);
      scd2->len = sc_len;
    }
  }
}

static const gchar *
_ncm_spline_cubic_d2_name (NcmSpline *s)
{
  return "NcmSplineCubicD2";
}

#define NUMCOSMO_CHECK_SPLINE_NODES 1

static void
_ncm_spline_cubic_d2_prepare_base (NcmSpline *s)
{
  NcmSplineCubic *sc     = NCM_SPLINE_CUBIC (s);
  NcmSplineCubicD2 *scd2 = NCM_SPLINE_CUBIC_D2 (s);
  NcmVector *c_vec       = ncm_spline_cubic_peek_c_vec (sc);

  ncm_vector_memcpy (c_vec, scd2->d2);
  ncm_vector_scale (c_vec, 0.5);

  return;
}

static void
_ncm_spline_cubic_d2_prepare (NcmSpline *s)
{
  NcmSplineCubic *sc = NCM_SPLINE_CUBIC (s);
  const size_t size  = ncm_spline_get_len (s);
  const size_t n     = size - 1;

  size_t i;

  _ncm_spline_cubic_d2_prepare_base (s);
  {
    NcmVector *xv    = ncm_spline_peek_xv (s);
    NcmVector *yv    = ncm_spline_peek_yv (s);
    NcmVector *b_vec = ncm_spline_cubic_peek_b_vec (sc);
    NcmVector *c_vec = ncm_spline_cubic_peek_c_vec (sc);
    NcmVector *d_vec = ncm_spline_cubic_peek_d_vec (sc);

    for (i = 0; i < n; i++)
    {
      const gdouble dx    = ncm_vector_get (xv, i + 1) - ncm_vector_get (xv, i);
      const gdouble dy    = ncm_vector_get (yv, i + 1) - ncm_vector_get (yv, i);
      const gdouble c_ip1 = ncm_vector_fast_get (c_vec, i + 1);
      const gdouble c_i   = ncm_vector_fast_get (c_vec, i);

      ncm_vector_fast_set (b_vec, i, (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0);
      ncm_vector_fast_set (d_vec, i, (c_ip1 - c_i) / (3.0 * dx));
    }
  }

  return;
}

static NcmSpline *
_ncm_spline_cubic_d2_copy_empty (const NcmSpline *s)
{
  NCM_UNUSED (s);

  return g_object_new (NCM_TYPE_SPLINE_CUBIC_D2, NULL);
}

static gsize
_ncm_spline_cubic_d2_min_size (const NcmSpline *s)
{
  NCM_UNUSED (s);

  return 2;
}

/**
 * ncm_spline_cubic_d2_new:
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of the values of the function to be interpolated, computed at @xv
 * @d2yv: #NcmVector of the values of the second derivative of the function to be interpolated, computed at @xv
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it
 *
 * This function returns a new #NcmSpline setting all its members.
 * It makes a copy of the @d2yv vector and saves it internally.
 *
 * Returns: a new #NcmSpline.
 */
NcmSplineCubicD2 *
ncm_spline_cubic_d2_new (NcmVector *xv, NcmVector *yv, NcmVector *d2yv, gboolean init)
{
  NcmSplineCubicD2 *scd2 = g_object_new (NCM_TYPE_SPLINE_CUBIC_D2, NULL);

  g_assert_cmpuint (ncm_vector_len (yv), ==, ncm_vector_len (d2yv));

  ncm_spline_set (NCM_SPLINE (scd2), xv, yv, FALSE);

  ncm_vector_memcpy (scd2->d2, d2yv);

  if (init)
    ncm_spline_prepare (NCM_SPLINE (scd2));

  return scd2;
}


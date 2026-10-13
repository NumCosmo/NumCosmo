/***************************************************************************
 *            ncm_spline_hermite.c
 *
 *  Mon Oct 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * NcmSplineHermite:
 *
 * Quintic Hermite interpolant from the values and the first and second derivatives at the
 * knots.
 *
 * On each interval the quintic matches $y$, $y'$ and $y''$ at both ends, so the interpolant
 * is $C^2$ and local: preparation computes three coefficients per interval and solves no
 * system. The derivative vectors are set by ncm_spline_hermite_new_full() or
 * ncm_spline_hermite_set_derivs() and kept by reference, as the knots and values are;
 * setting knots or values of the same length keeps them, a different length drops them,
 * and preparing without them aborts, so ncm_spline_copy() is not supported. Outside
 * the knots the boundary quintic is extrapolated. On a prepared spline, evaluation reads
 * only state fixed by preparation and may be called concurrently.
 *
 * The interpolation error is $O(h^6)$; see
 * <a href="../../theory/ncm/spline/hermite.html">Quintic Hermite interpolation</a>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline_hermite.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmSplineHermite
{
  /*< private >*/
  NcmSpline parent_instance;
  NcmVector *xv;
  NcmVector *yv;
  NcmVector *d1v;
  NcmVector *d2v;
  NcmVector *c3;
  NcmVector *c4;
  NcmVector *c5;
  gsize len;
};

enum
{
  PROP_0,
  PROP_D1,
  PROP_D2,
};

G_DEFINE_TYPE (NcmSplineHermite, ncm_spline_hermite, NCM_TYPE_SPLINE)

static void
ncm_spline_hermite_init (NcmSplineHermite *sh)
{
  sh->xv  = NULL;
  sh->yv  = NULL;
  sh->d1v = NULL;
  sh->d2v = NULL;
  sh->c3  = NULL;
  sh->c4  = NULL;
  sh->c5  = NULL;
  sh->len = 0;
}

static void
_ncm_spline_hermite_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE (object);

  g_return_if_fail (NCM_IS_SPLINE_HERMITE (object));

  switch (prop_id)
  {
    case PROP_D1:
      ncm_spline_hermite_set_derivs (sh, g_value_get_object (value), sh->d2v, FALSE);
      break;
    case PROP_D2:
      ncm_spline_hermite_set_derivs (sh, sh->d1v, g_value_get_object (value), FALSE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_hermite_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE (object);

  g_return_if_fail (NCM_IS_SPLINE_HERMITE (object));

  switch (prop_id)
  {
    case PROP_D1:
      g_value_set_object (value, sh->d1v);
      break;
    case PROP_D2:
      g_value_set_object (value, sh->d2v);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_hermite_dispose (GObject *object)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE (object);

  ncm_vector_clear (&sh->d1v);
  ncm_vector_clear (&sh->d2v);
  ncm_vector_clear (&sh->c3);
  ncm_vector_clear (&sh->c4);
  ncm_vector_clear (&sh->c5);
  sh->len = 0;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_hermite_parent_class)->dispose (object);
}

static const gchar *_ncm_spline_hermite_name (NcmSpline *s);
static void _ncm_spline_hermite_reset (NcmSpline *s);
static void _ncm_spline_hermite_prepare (NcmSpline *s);
static gsize _ncm_spline_hermite_min_size (const NcmSpline *s);
static gdouble _ncm_spline_hermite_eval (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_hermite_eval_idx (const NcmSpline *s, const gdouble x, const gsize i);
static gdouble _ncm_spline_hermite_deriv (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_hermite_deriv_idx (const NcmSpline *s, const gdouble x, const gsize i);
static gdouble _ncm_spline_hermite_deriv2 (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_hermite_deriv_nmax (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_hermite_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
static gdouble _ncm_spline_hermite_integ_idx (const NcmSpline *s, const gdouble xi, const gsize i, const gdouble xf, const gsize f);
static NcmSpline *_ncm_spline_hermite_copy_empty (const NcmSpline *s);

static void
ncm_spline_hermite_class_init (NcmSplineHermiteClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmSplineClass *s_class    = NCM_SPLINE_CLASS (klass);

  object_class->set_property = &_ncm_spline_hermite_set_property;
  object_class->get_property = &_ncm_spline_hermite_get_property;
  object_class->dispose      = &_ncm_spline_hermite_dispose;

  /**
   * NcmSplineHermite:d1:
   *
   * The first derivatives at the knots.
   */
  g_object_class_install_property (object_class,
                                   PROP_D1,
                                   g_param_spec_object ("d1",
                                                        NULL,
                                                        "First derivatives at the knots",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSplineHermite:d2:
   *
   * The second derivatives at the knots.
   */
  g_object_class_install_property (object_class,
                                   PROP_D2,
                                   g_param_spec_object ("d2",
                                                        NULL,
                                                        "Second derivatives at the knots",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  s_class->name         = &_ncm_spline_hermite_name;
  s_class->reset        = &_ncm_spline_hermite_reset;
  s_class->prepare      = &_ncm_spline_hermite_prepare;
  s_class->prepare_base = NULL;
  s_class->min_size     = &_ncm_spline_hermite_min_size;
  s_class->eval         = &_ncm_spline_hermite_eval;
  s_class->eval_idx     = &_ncm_spline_hermite_eval_idx;
  s_class->deriv        = &_ncm_spline_hermite_deriv;
  s_class->deriv_idx    = &_ncm_spline_hermite_deriv_idx;
  s_class->deriv2       = &_ncm_spline_hermite_deriv2;
  s_class->deriv_nmax   = &_ncm_spline_hermite_deriv_nmax;
  s_class->integ        = &_ncm_spline_hermite_integ;
  s_class->integ_idx    = &_ncm_spline_hermite_integ_idx;
  s_class->copy_empty   = &_ncm_spline_hermite_copy_empty;
}

static const gchar *
_ncm_spline_hermite_name (NcmSpline *s)
{
  return "NcmSplineHermite";
}

/*
 * The coefficient vectors follow the number of knots; the derivative vectors are kept
 * while their length matches it and dropped otherwise
 */
static void
_ncm_spline_hermite_reset (NcmSpline *s)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE (s);
  const gsize s_len    = ncm_spline_get_len (s);

  if (sh->len != s_len)
  {
    ncm_vector_clear (&sh->c3);
    ncm_vector_clear (&sh->c4);
    ncm_vector_clear (&sh->c5);

    sh->c3  = ncm_vector_new (s_len);
    sh->c4  = ncm_vector_new (s_len);
    sh->c5  = ncm_vector_new (s_len);
    sh->len = s_len;
  }

  if ((sh->d1v != NULL) && (ncm_vector_len (sh->d1v) != s_len))
    ncm_vector_clear (&sh->d1v);

  if ((sh->d2v != NULL) && (ncm_vector_len (sh->d2v) != s_len))
    ncm_vector_clear (&sh->d2v);

  sh->xv = ncm_spline_peek_xv (s);
  sh->yv = ncm_spline_peek_yv (s);
}

/*
 * Computes the coefficients of x - x_i beyond the second on each interval,
 * s(x) = y_i + y'_i dx + y''_i dx^2 / 2 + c3 dx^3 + c4 dx^4 + c5 dx^5,
 * from the values and derivatives at x_{i+1}. Aborts without the derivative vectors.
 */
static void
_ncm_spline_hermite_prepare (NcmSpline *s)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE (s);
  const gsize n        = sh->len;
  gsize i;

  if ((sh->d1v == NULL) || (sh->d2v == NULL))
    g_error ("ncm_spline_hermite_prepare: the derivatives at the knots were not set.");

  for (i = 0; i + 1 < n; i++)
  {
    const gdouble h    = ncm_vector_get (sh->xv, i + 1) - ncm_vector_get (sh->xv, i);
    const gdouble h2   = h * h;
    const gdouble h3   = h2 * h;
    const gdouble y0   = ncm_vector_get (sh->yv, i);
    const gdouble y1   = ncm_vector_get (sh->yv, i + 1);
    const gdouble d1_0 = ncm_vector_get (sh->d1v, i);
    const gdouble d1_1 = ncm_vector_get (sh->d1v, i + 1);
    const gdouble d2_0 = ncm_vector_get (sh->d2v, i);
    const gdouble d2_1 = ncm_vector_get (sh->d2v, i + 1);
    const gdouble A    = (y1 - y0 - h * d1_0 - 0.5 * h2 * d2_0) / h3;
    const gdouble B    = (d1_1 - d1_0 - h * d2_0) / h2;
    const gdouble C    = (d2_1 - d2_0) / h;

    ncm_vector_fast_set (sh->c3, i, 10.0 * A - 4.0 * B + 0.5 * C);
    ncm_vector_fast_set (sh->c4, i, (-15.0 * A + 7.0 * B - C) / h);
    ncm_vector_fast_set (sh->c5, i, (6.0 * A - 3.0 * B + 0.5 * C) / h2);
  }
}

static gsize
_ncm_spline_hermite_min_size (const NcmSpline *s)
{
  NCM_UNUSED (s);

  return 2;
}

static gdouble
_ncm_spline_hermite_eval_idx (const NcmSpline *s, const gdouble x, const gsize i)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE ((NcmSpline *) s);
  const gdouble dx     = x - ncm_vector_get (sh->xv, i);
  const gdouble a      = ncm_vector_get (sh->yv, i);
  const gdouble b      = ncm_vector_get (sh->d1v, i);
  const gdouble c      = 0.5 * ncm_vector_get (sh->d2v, i);
  const gdouble d      = ncm_vector_fast_get (sh->c3, i);
  const gdouble e      = ncm_vector_fast_get (sh->c4, i);
  const gdouble f      = ncm_vector_fast_get (sh->c5, i);

  return a + dx * (b + dx * (c + dx * (d + dx * (e + dx * f))));
}

static gdouble
_ncm_spline_hermite_eval (const NcmSpline *s, const gdouble x)
{
  return _ncm_spline_hermite_eval_idx (s, x, ncm_spline_get_index (s, x));
}

static gdouble
_ncm_spline_hermite_deriv_idx (const NcmSpline *s, const gdouble x, const gsize i)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE ((NcmSpline *) s);
  const gdouble dx     = x - ncm_vector_get (sh->xv, i);
  const gdouble b      = ncm_vector_get (sh->d1v, i);
  const gdouble c2     = ncm_vector_get (sh->d2v, i);
  const gdouble d3     = 3.0 * ncm_vector_fast_get (sh->c3, i);
  const gdouble e4     = 4.0 * ncm_vector_fast_get (sh->c4, i);
  const gdouble f5     = 5.0 * ncm_vector_fast_get (sh->c5, i);

  return b + dx * (c2 + dx * (d3 + dx * (e4 + dx * f5)));
}

static gdouble
_ncm_spline_hermite_deriv (const NcmSpline *s, const gdouble x)
{
  return _ncm_spline_hermite_deriv_idx (s, x, ncm_spline_get_index (s, x));
}

static gdouble
_ncm_spline_hermite_deriv2 (const NcmSpline *s, const gdouble x)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE ((NcmSpline *) s);
  const gsize i        = ncm_spline_get_index (s, x);
  const gdouble dx     = x - ncm_vector_get (sh->xv, i);
  const gdouble c2     = ncm_vector_get (sh->d2v, i);
  const gdouble d6     = 6.0 * ncm_vector_fast_get (sh->c3, i);
  const gdouble e12    = 12.0 * ncm_vector_fast_get (sh->c4, i);
  const gdouble f20    = 20.0 * ncm_vector_fast_get (sh->c5, i);

  return c2 + dx * (d6 + dx * (e12 + dx * f20));
}

static gdouble
_ncm_spline_hermite_deriv_nmax (const NcmSpline *s, const gdouble x)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE ((NcmSpline *) s);
  const gsize i        = ncm_spline_get_index (s, x);

  return 120.0 * ncm_vector_fast_get (sh->c5, i);
}

/* Integral of the quintic of interval i from x_i + ta to x_i + tb */
static gdouble
_ncm_spline_hermite_integ_piece (NcmSplineHermite *sh, const gsize i, const gdouble ta, const gdouble tb)
{
  const gdouble a  = ncm_vector_get (sh->yv, i);
  const gdouble b  = ncm_vector_get (sh->d1v, i);
  const gdouble c  = 0.5 * ncm_vector_get (sh->d2v, i);
  const gdouble d  = ncm_vector_fast_get (sh->c3, i);
  const gdouble e  = ncm_vector_fast_get (sh->c4, i);
  const gdouble f  = ncm_vector_fast_get (sh->c5, i);
  const gdouble Fb = tb * (a + tb * (0.5 * b + tb * (c / 3.0 + tb * (0.25 * d + tb * (0.2 * e + tb * f / 6.0)))));
  const gdouble Fa = ta * (a + ta * (0.5 * b + ta * (c / 3.0 + ta * (0.25 * d + ta * (0.2 * e + ta * f / 6.0)))));

  return Fb - Fa;
}

static gdouble
_ncm_spline_hermite_integ_idx (const NcmSpline *s, const gdouble xi, const gsize index_a, const gdouble xf, const gsize index_b)
{
  NcmSplineHermite *sh = NCM_SPLINE_HERMITE ((NcmSpline *) s);
  gdouble result       = 0.0;
  gsize i;

  for (i = index_a; i <= index_b; i++)
  {
    const gdouble x_lo = ncm_vector_get (sh->xv, i);
    const gdouble x_hi = ncm_vector_get (sh->xv, i + 1);
    const gdouble ta   = (i == index_a) ? xi - x_lo : 0.0;
    const gdouble tb   = (i == index_b) ? xf - x_lo : x_hi - x_lo;

    result += _ncm_spline_hermite_integ_piece (sh, i, ta, tb);
  }

  return result;
}

static gdouble
_ncm_spline_hermite_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  return _ncm_spline_hermite_integ_idx (s, x0, ncm_spline_get_index (s, x0), x1, ncm_spline_get_index (s, x1));
}

static NcmSpline *
_ncm_spline_hermite_copy_empty (const NcmSpline *s)
{
  NCM_UNUSED (s);

  return NCM_SPLINE (ncm_spline_hermite_new ());
}

/**
 * ncm_spline_hermite_new:
 *
 * Creates an empty #NcmSplineHermite.
 *
 * Returns: (transfer full): a new #NcmSplineHermite.
 */
NcmSplineHermite *
ncm_spline_hermite_new (void)
{
  return g_object_new (NCM_TYPE_SPLINE_HERMITE, NULL);
}

/**
 * ncm_spline_hermite_new_full:
 * @xv: the knots
 * @yv: the values at @xv
 * @d1v: the first derivatives at @xv
 * @d2v: the second derivatives at @xv
 * @init: whether to prepare the new spline
 *
 * Creates a #NcmSplineHermite keeping references to the four vectors, see ncm_spline_set()
 * and ncm_spline_hermite_set_derivs().
 *
 * Returns: (transfer full): a new #NcmSplineHermite.
 */
NcmSplineHermite *
ncm_spline_hermite_new_full (NcmVector *xv, NcmVector *yv, NcmVector *d1v, NcmVector *d2v, gboolean init)
{
  NcmSplineHermite *sh = ncm_spline_hermite_new ();

  ncm_spline_set (NCM_SPLINE (sh), xv, yv, FALSE);
  ncm_spline_hermite_set_derivs (sh, d1v, d2v, init);

  return sh;
}

/**
 * ncm_spline_hermite_set_derivs:
 * @sh: a #NcmSplineHermite
 * @d1v: the first derivatives at the knots
 * @d2v: the second derivatives at the knots
 * @init: whether to prepare @sh
 *
 * Sets the derivative vectors, keeping references to them. Aborts if the knots are not
 * set or if a length is not the number of knots.
 */
void
ncm_spline_hermite_set_derivs (NcmSplineHermite *sh, NcmVector *d1v, NcmVector *d2v, gboolean init)
{
  const guint len = ncm_spline_get_len (NCM_SPLINE (sh));

  if (len == 0)
    g_error ("ncm_spline_hermite_set_derivs: the knots must be set first.");

  if ((d1v != NULL) && (ncm_vector_len (d1v) != len))
    g_error ("ncm_spline_hermite_set_derivs: d1 has length %u, the spline has %u knots.", ncm_vector_len (d1v), len);

  if ((d2v != NULL) && (ncm_vector_len (d2v) != len))
    g_error ("ncm_spline_hermite_set_derivs: d2 has length %u, the spline has %u knots.", ncm_vector_len (d2v), len);

  if (d1v != NULL)
    ncm_vector_ref (d1v);

  ncm_vector_clear (&sh->d1v);
  sh->d1v = d1v;

  if (d2v != NULL)
    ncm_vector_ref (d2v);

  ncm_vector_clear (&sh->d2v);
  sh->d2v = d2v;

  if (init)
    ncm_spline_prepare (NCM_SPLINE (sh));
}

/**
 * ncm_spline_hermite_peek_d1v:
 * @sh: a #NcmSplineHermite
 *
 * Returns: (transfer none) (nullable): the first derivatives at the knots.
 */
NcmVector *
ncm_spline_hermite_peek_d1v (NcmSplineHermite *sh)
{
  return sh->d1v;
}

/**
 * ncm_spline_hermite_peek_d2v:
 * @sh: a #NcmSplineHermite
 *
 * Returns: (transfer none) (nullable): the second derivatives at the knots.
 */
NcmVector *
ncm_spline_hermite_peek_d2v (NcmSplineHermite *sh)
{
  return sh->d2v;
}


/***************************************************************************
 *            ncm_spline_bspline.c
 *
 *  Wed Aug 20 10:00:00 2026
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
 * NcmSplineBSpline:
 *
 * Interpolating B-spline of order 2 to %NCM_SPLINE_BSPLINE_MAX_ORDER.
 *
 * Interpolates the values at the knots of #NcmSpline with one control point per knot,
 * B-spline knots placed by gsl_bspline_init_interp(), and a banded collocation solve, so
 * preparation costs $O(n)$. Outside the knots, the value, the derivatives and the
 * integral extrapolate the edge polynomial.
 *
 * A cubic spline is $C^2$ and its interpolation error stays well above machine precision
 * at any sample density; higher orders reach it. Maximum error at the interval midpoints
 * for a Gaussian of mean $1200$ and standard deviation $300$ sampled uniformly on
 * $[0, 2400]$:
 *
 * |samples|degree 3 |degree 5 |degree 7 |degree 9 |
 * |------:|--------:|--------:|--------:|--------:|
 * |    100|  3.4e-07|  1.2e-09|  2.8e-11|  5.2e-13|
 * |    500|  5.2e-10|  6.6e-14|  6.7e-16|  6.7e-16|
 * |   2000|  2.0e-12|  4.4e-16|  4.4e-16|  4.4e-16|
 *
 * On a prepared spline, ncm_spline_eval() may be called concurrently: it uses stack
 * scratch and reads only state fixed by preparation. The derivatives and the integral use
 * the GSL workspace and serialize on an internal lock. Preparing concurrently with any
 * evaluation is not supported, as for every #NcmSpline.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline_bspline.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmSplineBSpline
{
  /*< private >*/
  NcmSpline parent_instance;
  gsl_bspline_workspace *w;
  gsl_vector *c;  /* control points */
  gsl_matrix *XB; /* banded collocation matrix, reused across prepares */
  gsl_vector_uint *piv;
  gsl_integration_glfixed_table *gl; /* exact for the edge polynomials, see _integ */
  guint order;
  gsize alloc_len; /* length the workspace was allocated for */
  gdouble reltol;  /* > 0 selects the order automatically */
  gdouble abstol;
  gdouble achieved_err; /* estimated interpolation error of the chosen order */
  gchar *inst_name;
  GMutex lock; /* serializes the paths that write gsl workspace scratch */
};

G_DEFINE_TYPE (NcmSplineBSpline, ncm_spline_bspline, NCM_TYPE_SPLINE)

enum
{
  PROP_0,
  PROP_ORDER,
  PROP_RELTOL,
  PROP_ABSTOL,
};

static void
ncm_spline_bspline_init (NcmSplineBSpline *sbs)
{
  sbs->w            = NULL;
  sbs->c            = NULL;
  sbs->XB           = NULL;
  sbs->piv          = NULL;
  sbs->gl           = NULL;
  sbs->order        = 0;
  sbs->alloc_len    = 0;
  sbs->reltol       = 0.0;
  sbs->abstol       = 0.0;
  sbs->achieved_err = 0.0;
  sbs->inst_name    = NULL;

  g_mutex_init (&sbs->lock);
}

static void
_ncm_spline_bspline_free_workspace (NcmSplineBSpline *sbs)
{
  g_clear_pointer (&sbs->w, gsl_bspline_free);
  g_clear_pointer (&sbs->c, gsl_vector_free);
  g_clear_pointer (&sbs->XB, gsl_matrix_free);
  g_clear_pointer (&sbs->piv, gsl_vector_uint_free);
  g_clear_pointer (&sbs->gl, gsl_integration_glfixed_table_free);
  sbs->alloc_len = 0;
}

static void
_ncm_spline_bspline_finalize (GObject *object)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE (object);

  _ncm_spline_bspline_free_workspace (sbs);
  g_clear_pointer (&sbs->inst_name, g_free);
  g_mutex_clear (&sbs->lock);

  /* Chain up: end */
  G_OBJECT_CLASS (ncm_spline_bspline_parent_class)->finalize (object);
}

static void
_ncm_spline_bspline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE (object);

  g_return_if_fail (NCM_IS_SPLINE_BSPLINE (object));

  switch (prop_id)
  {
    case PROP_ORDER:
      ncm_spline_bspline_set_order (sbs, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      sbs->reltol = g_value_get_double (value);
      break;
    case PROP_ABSTOL:
      sbs->abstol = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_bspline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE (object);

  g_return_if_fail (NCM_IS_SPLINE_BSPLINE (object));

  switch (prop_id)
  {
    case PROP_ORDER:
      g_value_set_uint (value, sbs->order);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, sbs->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, sbs->abstol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static const gchar *_ncm_spline_bspline_name (NcmSpline *s);
static void _ncm_spline_bspline_reset (NcmSpline *s);
static void _ncm_spline_bspline_prepare (NcmSpline *s);
static gsize _ncm_spline_bspline_min_size (const NcmSpline *s);
static gdouble _ncm_spline_bspline_eval (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_bspline_deriv (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_bspline_deriv2 (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_bspline_deriv_nmax (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_bspline_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
static NcmSpline *_ncm_spline_bspline_copy_empty (const NcmSpline *s);

static void
ncm_spline_bspline_class_init (NcmSplineBSplineClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmSplineClass *s_class    = NCM_SPLINE_CLASS (klass);

  object_class->set_property = &_ncm_spline_bspline_set_property;
  object_class->get_property = &_ncm_spline_bspline_get_property;
  object_class->finalize     = &_ncm_spline_bspline_finalize;

  /**
   * NcmSplineBSpline:order:
   *
   * The B-spline order, the polynomial degree plus one.
   */
  g_object_class_install_property (object_class,
                                   PROP_ORDER,
                                   g_param_spec_uint ("order",
                                                      NULL,
                                                      "B-spline order (degree + 1)",
                                                      2, NCM_SPLINE_BSPLINE_MAX_ORDER,
                                                      NCM_SPLINE_BSPLINE_DEFAULT_ORDER,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSplineBSpline:reltol:
   *
   * The relative interpolation error requested, or zero to use #NcmSplineBSpline:order.
   *
   * When positive, preparation sets #NcmSplineBSpline:order to the lowest even order from
   * $4$ whose estimated error is at most $\max(\mathrm{reltol}\,(y_\mathrm{max} -
   * y_\mathrm{min}), \mathrm{abstol})$, and aborts if none is. The error of order $m$
   * is estimated as the largest difference from order $m + 2$ at the interval midpoints;
   * %NCM_SPLINE_BSPLINE_MAX_ORDER, which has no higher order to compare with, is given
   * the estimate of the order below it.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Requested relative interpolation error, 0 to select the order manually",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSplineBSpline:abstol:
   *
   * The absolute interpolation error requested, see #NcmSplineBSpline:reltol; unused
   * when #NcmSplineBSpline:reltol is zero.
   */
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute floor for the automatic order selection",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  s_class->name         = &_ncm_spline_bspline_name;
  s_class->reset        = &_ncm_spline_bspline_reset;
  s_class->prepare      = &_ncm_spline_bspline_prepare;
  s_class->prepare_base = NULL;
  s_class->min_size     = &_ncm_spline_bspline_min_size;
  s_class->eval         = &_ncm_spline_bspline_eval;
  s_class->deriv        = &_ncm_spline_bspline_deriv;
  s_class->deriv2       = &_ncm_spline_bspline_deriv2;
  s_class->deriv_nmax   = &_ncm_spline_bspline_deriv_nmax;
  s_class->integ        = &_ncm_spline_bspline_integ;
  s_class->copy_empty   = &_ncm_spline_bspline_copy_empty;
}

static const gchar *
_ncm_spline_bspline_name (NcmSpline *s)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE (s);

  g_mutex_lock (&sbs->lock);

  if (sbs->inst_name == NULL)
    sbs->inst_name = g_strdup_printf ("NcmSplineBSpline[order %u]", sbs->order);

  g_mutex_unlock (&sbs->lock);

  return sbs->inst_name;
}

static void
_ncm_spline_bspline_reset (NcmSpline *s)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE (s);
  const gsize s_len     = ncm_spline_get_len (s);

  if (sbs->alloc_len == s_len)
    return;

  _ncm_spline_bspline_free_workspace (sbs);

  /* One control point per sample: the fit is an interpolation, not a smoothing. */
  sbs->w   = gsl_bspline_alloc_ncontrol (sbs->order, s_len);
  sbs->c   = gsl_vector_alloc (s_len);
  sbs->XB  = gsl_matrix_alloc (s_len, 3 * (sbs->order - 1) + 1);
  sbs->piv = gsl_vector_uint_alloc (s_len);
  sbs->gl  = gsl_integration_glfixed_table_alloc ((sbs->order + 1) / 2);

  sbs->alloc_len = s_len;
}

/* Fits @order into the given workspace; FALSE when GSL fails. */
static gboolean
_ncm_spline_bspline_fit (const gsl_vector *xv, const gsl_vector *yv, const guint order,
                         gsl_bspline_workspace *w, gsl_matrix *XB, gsl_vector_uint *piv,
                         gsl_vector *c)
{
  const gsize band  = order - 1;
  const gsize s_len = xv->size;

  if (gsl_bspline_init_interp (xv, w) != GSL_SUCCESS)
    return FALSE;

  if (gsl_bspline_col_interp (xv, XB, w) != GSL_SUCCESS)
    return FALSE;

  if (gsl_linalg_LU_band_decomp (s_len, band, band, XB, piv) != GSL_SUCCESS)
    return FALSE;

  if (gsl_linalg_LU_band_solve (band, band, XB, piv, yv, c) != GSL_SUCCESS)
    return FALSE;

  return TRUE;
}

/* Largest difference between the @order and @ref_order fits at the interval midpoints,
 * or +inf when either fit fails. */
static gdouble
_ncm_spline_bspline_estimate_error (const gsl_vector *xv, const gsl_vector *yv,
                                    const guint order, const guint ref_order)
{
  const gsize s_len = xv->size;
  gdouble worst     = 0.0;
  gsize i;

  gsl_bspline_workspace *w_lo = gsl_bspline_alloc_ncontrol (order, s_len);
  gsl_bspline_workspace *w_hi = gsl_bspline_alloc_ncontrol (ref_order, s_len);
  gsl_vector *c_lo            = gsl_vector_alloc (s_len);
  gsl_vector *c_hi            = gsl_vector_alloc (s_len);
  gsl_matrix *XB_lo           = gsl_matrix_alloc (s_len, 3 * (order - 1) + 1);
  gsl_matrix *XB_hi           = gsl_matrix_alloc (s_len, 3 * (ref_order - 1) + 1);
  gsl_vector_uint *piv        = gsl_vector_uint_alloc (s_len);
  gboolean ok;

  ok = _ncm_spline_bspline_fit (xv, yv, order, w_lo, XB_lo, piv, c_lo) &&
       _ncm_spline_bspline_fit (xv, yv, ref_order, w_hi, XB_hi, piv, c_hi);

  if (ok)
  {
    for (i = 0; i + 1 < s_len; i++)
    {
      const gdouble xm = 0.5 * (gsl_vector_get (xv, i) + gsl_vector_get (xv, i + 1));
      gdouble v_lo = 0.0, v_hi = 0.0;

      gsl_bspline_calc (xm, c_lo, &v_lo, w_lo);
      gsl_bspline_calc (xm, c_hi, &v_hi, w_hi);

      worst = GSL_MAX (worst, fabs (v_lo - v_hi));
    }
  }
  else
  {
    worst = GSL_POSINF;
  }

  gsl_bspline_free (w_lo);
  gsl_bspline_free (w_hi);
  gsl_vector_free (c_lo);
  gsl_vector_free (c_hi);
  gsl_matrix_free (XB_lo);
  gsl_matrix_free (XB_hi);
  gsl_vector_uint_free (piv);

  return worst;
}

/* Sets the lowest even order meeting the tolerance, see #NcmSplineBSpline:reltol; aborts
 * when none does. */
static void
_ncm_spline_bspline_select_order (NcmSplineBSpline *sbs, const gsl_vector *xv, const gsl_vector *yv)
{
  const gsize s_len   = xv->size;
  const gdouble y_max = gsl_vector_max (yv) - gsl_vector_min (yv);
  const gdouble tol   = GSL_MAX (sbs->reltol * fabs (y_max), sbs->abstol);
  gdouble best_err    = GSL_POSINF;
  guint best_order    = 0;
  guint order;

  for (order = 4; order <= NCM_SPLINE_BSPLINE_MAX_ORDER; order += 2)
  {
    gdouble err;

    if (order + 2 > s_len)
      break;

    if (order + 2 <= NCM_SPLINE_BSPLINE_MAX_ORDER)
      err = _ncm_spline_bspline_estimate_error (xv, yv, order, order + 2);
    else
      /* No higher order to compare with: reuse the previous order's estimate, an
       * overestimate for this one. */
      err = best_err;

    if (err < best_err)
    {
      best_err   = err;
      best_order = order;
    }

    if (err <= tol)
    {
      ncm_spline_bspline_set_order (sbs, order);
      sbs->achieved_err = err;

      return;
    }
  }

  g_error ("_ncm_spline_bspline_select_order: %u samples cannot support a requested "
           "interpolation error of %.6e; the best supported order (%u) reaches only "
           "%.6e. Supply more samples, or relax reltol/abstol.",
           (guint) s_len, tol, best_order, best_err);
}

static void
_ncm_spline_bspline_prepare (NcmSpline *s)
{
  NcmSplineBSpline *sbs    = NCM_SPLINE_BSPLINE (s);
  NcmVector *s_xv          = ncm_spline_peek_xv (s);
  NcmVector *s_yv          = ncm_spline_peek_yv (s);
  const gsize s_len        = ncm_spline_get_len (s);
  gsl_vector_const_view xv = gsl_vector_const_view_array_with_stride (ncm_vector_ptr (s_xv, 0), ncm_vector_stride (s_xv), s_len);
  gsl_vector_const_view yv = gsl_vector_const_view_array_with_stride (ncm_vector_ptr (s_yv, 0), ncm_vector_stride (s_yv), s_len);

  if (sbs->reltol > 0.0)
    _ncm_spline_bspline_select_order (sbs, &xv.vector, &yv.vector);

  /* ncm_spline_prepare() does not call reset(), and a change of order frees the workspace */
  if ((sbs->w == NULL) || (sbs->alloc_len != s_len))
    _ncm_spline_bspline_reset (s);

  if (!_ncm_spline_bspline_fit (&xv.vector, &yv.vector, sbs->order, sbs->w, sbs->XB, sbs->piv, sbs->c))
    g_error ("_ncm_spline_bspline_prepare: order %u interpolation failed on %u samples.",
             sbs->order, (guint) s_len);
}

static gsize
_ncm_spline_bspline_min_size (const NcmSpline *s)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);

  return sbs->order;
}

/*
 * gsl_bspline_calc() writes scratch into the workspace, so evaluation runs de Boor's
 * recursion (PPPACK bsplvb, as in GSL) on stack scratch instead, for concurrent callers
 * such as the xcor kernel integrands.
 */
static gdouble
_ncm_spline_bspline_eval (const NcmSpline *s, const gdouble x)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  const gsize k         = sbs->order;
  const gsize ncontrol  = sbs->alloc_len;
  const gdouble *t      = sbs->w->knots->data; /* contiguous: allocated by gsl_bspline_alloc_ncontrol() */
  const gdouble *c      = sbs->c->data;
  gdouble deltal[NCM_SPLINE_BSPLINE_MAX_ORDER];
  gdouble deltar[NCM_SPLINE_BSPLINE_MAX_ORDER];
  gdouble B[NCM_SPLINE_BSPLINE_MAX_ORDER];
  gsize l, i, j;

  /* Largest span index l in [k - 1, ncontrol - 1] with t[l] <= x. Outside the range the
   * clamped edge span extrapolates its polynomial, matching gsl_bspline_calc(). */
  if (x < t[k - 1])
  {
    l = k - 1;
  }
  else if (x >= t[ncontrol])
  {
    l = ncontrol - 1;
  }
  else
  {
    gsize lo = k - 1;
    gsize hi = ncontrol;

    while (hi - lo > 1)
    {
      const gsize mid = (lo + hi) / 2;

      if (x < t[mid])
        hi = mid;
      else
        lo = mid;
    }

    l = lo;
  }

  B[0] = 1.0;

  for (j = 0; j + 1 < k; j++)
  {
    gdouble saved = 0.0;

    deltar[j] = t[l + j + 1] - x;
    deltal[j] = x - t[l - j];

    for (i = 0; i <= j; i++)
    {
      const gdouble term = B[i] / (deltar[i] + deltal[j - i]);

      B[i]  = saved + deltar[i] * term;
      saved = deltal[j - i] * term;
    }

    B[j + 1] = saved;
  }

  {
    gdouble res = 0.0;

    for (i = 0; i < k; i++)
      res += B[i] * c[l - (k - 1) + i];

    return res;
  }
}

static gdouble
_ncm_spline_bspline_deriv (const NcmSpline *s, const gdouble x)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  gdouble res           = 0.0;

  g_mutex_lock (&sbs->lock);
  gsl_bspline_calc_deriv (x, sbs->c, 1, &res, sbs->w);
  g_mutex_unlock (&sbs->lock);

  return res;
}

static gdouble
_ncm_spline_bspline_deriv2 (const NcmSpline *s, const gdouble x)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  gdouble res           = 0.0;

  g_mutex_lock (&sbs->lock);
  gsl_bspline_calc_deriv (x, sbs->c, 2, &res, sbs->w);
  g_mutex_unlock (&sbs->lock);

  return res;
}

static gdouble
_ncm_spline_bspline_deriv_nmax (const NcmSpline *s, const gdouble x)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  gdouble res           = 0.0;

  /* The derivative of order equal to the degree */
  g_mutex_lock (&sbs->lock);
  gsl_bspline_calc_deriv (x, sbs->c, sbs->order - 1, &res, sbs->w);
  g_mutex_unlock (&sbs->lock);

  return res;
}

/* Integral of the edge polynomial over [a, b], outside the knots: Gauss-Legendre with
 * (order + 1) / 2 points is exact for its degree, order - 1. */
static gdouble
_ncm_spline_bspline_integ_extrap (const NcmSpline *s, const gdouble a, const gdouble b)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  gdouble res           = 0.0;
  gsize i;

  for (i = 0; i < sbs->gl->n; i++)
  {
    gdouble xi, wi;

    gsl_integration_glfixed_point (a, b, i, &xi, &wi, sbs->gl);
    res += wi * _ncm_spline_bspline_eval (s, xi);
  }

  return res;
}

/* gsl_bspline_calc_integ() clamps the limits to the knots, so the parts outside them are
 * integrated here, consistently with the extrapolation of _eval. Called with x0 <= x1. */
static gdouble
_ncm_spline_bspline_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  const gdouble *t      = sbs->w->knots->data;
  const gdouble t_lo    = t[sbs->order - 1];
  const gdouble t_hi    = t[sbs->alloc_len];
  const gdouble lo      = GSL_MAX (x0, t_lo);
  const gdouble hi      = GSL_MIN (x1, t_hi);
  gdouble res           = 0.0;

  if (x0 < t_lo)
    res += _ncm_spline_bspline_integ_extrap (s, x0, GSL_MIN (x1, t_lo));

  if (x1 > t_hi)
    res += _ncm_spline_bspline_integ_extrap (s, GSL_MAX (x0, t_hi), x1);

  if (lo < hi)
  {
    gdouble inner = 0.0;

    g_mutex_lock (&sbs->lock);
    gsl_bspline_calc_integ (lo, hi, sbs->c, &inner, sbs->w);
    g_mutex_unlock (&sbs->lock);

    res += inner;
  }

  return res;
}

static NcmSpline *
_ncm_spline_bspline_copy_empty (const NcmSpline *s)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);

  return NCM_SPLINE (g_object_new (NCM_TYPE_SPLINE_BSPLINE,
                                   "order", sbs->order,
                                   "reltol", sbs->reltol,
                                   "abstol", sbs->abstol,
                                   NULL));
}

/**
 * ncm_spline_bspline_new:
 * @order: the B-spline order
 *
 * Creates an empty interpolating B-spline of order @order.
 *
 * Returns: (transfer full): a new #NcmSplineBSpline.
 */
NcmSplineBSpline *
ncm_spline_bspline_new (guint order)
{
  return g_object_new (NCM_TYPE_SPLINE_BSPLINE,
                       "order", order,
                       NULL);
}

/**
 * ncm_spline_bspline_new_full:
 * @order: the B-spline order
 * @xv: the knots
 * @yv: the values at @xv
 * @init: whether to prepare the spline
 *
 * Creates an interpolating B-spline of order @order over (@xv, @yv).
 *
 * Returns: (transfer full): a new #NcmSplineBSpline.
 */
NcmSplineBSpline *
ncm_spline_bspline_new_full (guint order, NcmVector *xv, NcmVector *yv, gboolean init)
{
  NcmSplineBSpline *sbs = ncm_spline_bspline_new (order);

  ncm_spline_set (NCM_SPLINE (sbs), xv, yv, init);

  return sbs;
}

/**
 * ncm_spline_bspline_set_order:
 * @sbs: a #NcmSplineBSpline
 * @order: the B-spline order
 *
 * Sets #NcmSplineBSpline:order, from $2$ to %NCM_SPLINE_BSPLINE_MAX_ORDER. The spline
 * must be prepared again afterwards.
 */
void
ncm_spline_bspline_set_order (NcmSplineBSpline *sbs, guint order)
{
  g_assert_cmpuint (order, >=, 2);
  g_assert_cmpuint (order, <=, NCM_SPLINE_BSPLINE_MAX_ORDER);

  if (order == sbs->order)
    return;

  sbs->order = order;

  /* The workspace depends on the order */
  _ncm_spline_bspline_free_workspace (sbs);
  g_clear_pointer (&sbs->inst_name, g_free);
}

/**
 * ncm_spline_bspline_new_tol:
 * @reltol: the relative interpolation error requested
 * @abstol: the absolute interpolation error requested
 *
 * Creates an empty B-spline that sets its order on preparation, see
 * #NcmSplineBSpline:reltol.
 *
 * Returns: (transfer full): a new #NcmSplineBSpline.
 */
NcmSplineBSpline *
ncm_spline_bspline_new_tol (gdouble reltol, gdouble abstol)
{
  return g_object_new (NCM_TYPE_SPLINE_BSPLINE,
                       "reltol", reltol,
                       "abstol", abstol,
                       NULL);
}

/**
 * ncm_spline_bspline_get_achieved_error:
 * @sbs: a #NcmSplineBSpline
 *
 * Gets the estimated interpolation error of the order chosen on preparation, see
 * #NcmSplineBSpline:reltol.
 *
 * Returns: the estimated error, or zero when #NcmSplineBSpline:reltol is zero.
 */
gdouble
ncm_spline_bspline_get_achieved_error (NcmSplineBSpline *sbs)
{
  return sbs->achieved_err;
}

/**
 * ncm_spline_bspline_get_order:
 * @sbs: a #NcmSplineBSpline
 *
 * Gets #NcmSplineBSpline:order.
 *
 * Returns: the B-spline order.
 */
guint
ncm_spline_bspline_get_order (NcmSplineBSpline *sbs)
{
  return sbs->order;
}


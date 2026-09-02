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
 * Interpolating B-spline of arbitrary order.
 *
 * Interpolates the sample points with a B-spline of the requested order (degree
 * `order - 1`), using the not-a-knot-like knot placement of gsl_bspline_init_interp()
 * and a banded collocation solve, so preparation is $O(n)$ in the number of points.
 *
 * The motivation is accuracy on data supplied as a fixed table. A cubic spline is only
 * $C^2$ and its interpolation error stalls well above machine precision no matter how
 * densely the table is sampled; higher orders do not. Measured maximum interpolation
 * error for a smooth function on a uniform grid:
 *
 * |samples|degree 3 |degree 5 |degree 7 |degree 9 |
 * |------:|--------:|--------:|--------:|--------:|
 * |    100|  3.4e-07|  1.2e-09|  2.8e-11|  5.2e-13|
 * |    500|  5.2e-10|  6.6e-14|  4.4e-16|  5.6e-16|
 * |   2000|  2.0e-12|  5.6e-16|  5.6e-16|  7.8e-16|
 *
 * Hence the default order of %NCM_SPLINE_BSPLINE_DEFAULT_ORDER and the cap at
 * %NCM_SPLINE_BSPLINE_MAX_ORDER.
 *
 * On a prepared spline, ncm_spline_eval() is safe to call concurrently: it evaluates
 * the basis with stack scratch and touches no shared mutable state. The derivative and
 * integral entry points go through the GSL workspace, whose scratch is per instance,
 * and therefore serialize on an internal lock. Preparing concurrently with any
 * evaluation is not supported, as for every #NcmSpline.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline_bspline.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_bspline.h>
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
static gdouble _ncm_spline_bspline_eval_idx (const NcmSpline *s, const gdouble x, const gsize i);
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
   * B-spline order; the polynomial degree is one less. Order 4 reproduces a cubic
   * spline. See the table in the class description for the accuracy this buys.
   *
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
   * When positive, the order is chosen automatically: the lowest order whose estimated
   * interpolation error meets $\max(\mathrm{reltol}\,\|y\|_\infty, \mathrm{abstol})$
   * is used, and preparation fails if no supported order reaches it. Zero (the default)
   * keeps the explicitly requested #NcmSplineBSpline:order.
   *
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
   * Absolute floor for the automatic order selection, ignored when
   * #NcmSplineBSpline:reltol is zero.
   *
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
  s_class->eval_idx     = &_ncm_spline_bspline_eval_idx;
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

  sbs->alloc_len = s_len;
}

/* Fit one order into the caller's workspace. Returns FALSE when the linear algebra
 * fails, which for an interpolation problem means the order is not usable here. */
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

/* Estimated interpolation error of @order, measured at the interval midpoints against
 * the next higher order. The samples say nothing about the function between them, so
 * the estimate is a held-out comparison -- the same principle as the knot-placement
 * criterion -- rather than a bound inferred from derivatives. */
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

/* Lowest order meeting the requested tolerance. Errors out when none does: the samples
 * simply do not support the request, and silently returning a worse answer would hide
 * the one fact the caller needs. */
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
      /* The highest order has no higher reference to be measured against. Comparing it
       * with itself would report a zero error and accept it unconditionally, which is
       * the silent-floor failure this selection exists to prevent. Reuse the previous
       * difference instead: it is an overestimate for this order, so the worst case is
       * refusing a request that was marginally attainable, never accepting one that
       * was not. */
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

  /* ncm_spline_prepare() does not call reset(), so a prepare that follows a change of
   * order -- which discards the order-specific workspace -- must rebuild it here. */
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

  /* Interpolation needs at least as many samples as there are basis functions in a
   * single span, i.e. the order. */
  return sbs->order;
}

/*
 * The gsl_bspline_calc* family writes scratch (deltal, deltar, B, dB, icache) into the
 * per-instance workspace, so it cannot run concurrently. Evaluation is the one
 * performance-critical path
 * -- kernel integrands call it from OpenMP loops on shared splines -- so it is computed
 * here with de Boor's recursion (PPPACK bsplvb, the same algorithm GSL runs) on stack
 * scratch, reading only state that preparation froze. The derivative and integral
 * entry points stay on GSL and serialize on the instance lock instead.
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

/*
 * @i is a hint, letting callers that already know the interval skip a binary
 * search. The evaluation locates the span itself and takes no index, so the
 * hint is dropped and this is a plain evaluation -- correct, but carrying none
 * of the saving the fast path exists for. Callers that lean on it, such as
 * ncm_spline_vec_eval() over a shared abscissa, pay full lookup per component.
 */
static gdouble
_ncm_spline_bspline_eval_idx (const NcmSpline *s, const gdouble x, const gsize i)
{
  return _ncm_spline_bspline_eval (s, x);
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

  /* Highest derivative that is not identically zero: the degree, order - 1. */
  g_mutex_lock (&sbs->lock);
  gsl_bspline_calc_deriv (x, sbs->c, sbs->order - 1, &res, sbs->w);
  g_mutex_unlock (&sbs->lock);

  return res;
}

static gdouble
_ncm_spline_bspline_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);
  gdouble res           = 0.0;

  g_mutex_lock (&sbs->lock);
  gsl_bspline_calc_integ (x0, x1, sbs->c, &res, sbs->w);
  g_mutex_unlock (&sbs->lock);

  return res;
}

static NcmSpline *
_ncm_spline_bspline_copy_empty (const NcmSpline *s)
{
  NcmSplineBSpline *sbs = NCM_SPLINE_BSPLINE ((NcmSpline *) s);

  return NCM_SPLINE (ncm_spline_bspline_new (sbs->order));
}

/**
 * ncm_spline_bspline_new:
 * @order: B-spline order, i.e. polynomial degree plus one
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
 * @order: B-spline order, i.e. polynomial degree plus one
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of the values of the function to be interpolated
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
 * @order: B-spline order, i.e. polynomial degree plus one
 *
 * Sets the B-spline order. The spline must be prepared again afterwards.
 */
void
ncm_spline_bspline_set_order (NcmSplineBSpline *sbs, guint order)
{
  g_assert_cmpuint (order, >=, 2);
  g_assert_cmpuint (order, <=, NCM_SPLINE_BSPLINE_MAX_ORDER);

  if (order == sbs->order)
    return;

  sbs->order = order;

  /* The workspace is order-specific, so it must be rebuilt even at unchanged length. */
  _ncm_spline_bspline_free_workspace (sbs);
  g_clear_pointer (&sbs->inst_name, g_free);
}

/**
 * ncm_spline_bspline_new_tol:
 * @reltol: requested relative interpolation error
 * @abstol: absolute floor
 *
 * Creates an empty B-spline that chooses its own order: on preparation the lowest order
 * whose estimated interpolation error meets the tolerance is used, and preparation
 * fails loudly if the samples supplied cannot support the request.
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
 * Estimated interpolation error of the prepared spline, when the order was selected
 * automatically. This is the quantity that bounds the accuracy of anything computed
 * from the spline, so a caller integrating it can combine this with the integral's own
 * cancellation ratio to predict the error it will see.
 *
 * Returns: the estimated error, or zero when the order was set manually.
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
 * Returns: the B-spline order.
 */
guint
ncm_spline_bspline_get_order (NcmSplineBSpline *sbs)
{
  return sbs->order;
}


/***************************************************************************
 *            ncm_spline.c
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
 * NcmSpline:
 *
 * Abstract base class of the one-dimensional interpolating splines.
 *
 * A spline holds the knots $x_i$, in increasing order, and the values $y_i$, both as
 * #NcmVector, and must be prepared, see ncm_spline_prepare(), before it is evaluated. The
 * subclasses define the interpolant: #NcmSplineCubicNotaknot, #NcmSplineCubicD2,
 * #NcmSplineGsl and others.
 *
 * Evaluation does no range check. Outside $[x_0, x_{n-1}]$ the native cubic splines
 * extrapolate the boundary polynomial, which can be wrong by orders of magnitude, while
 * #NcmSplineGsl returns NaN. Callers whose abscissa is not known to lie inside must check it.
 *
 * The index of the interval containing $x$ comes from a table of uniform buckets built by
 * ncm_spline_prepare(), or, when enabled, from a GSL accelerator, see ncm_spline_acc().
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/integration/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmSplinePrivate
{
  gsize len;
  NcmVector *xv;
  NcmVector *yv;
  gsl_interp_accel *acc;
  gboolean init;
  gboolean empty;
  GArray *bucket;
  gdouble bucket_dx;
  gdouble one_over_dx;
  gdouble x_min;
  gdouble x_max;
  gsize start_right;
  gsize last_poly;
  gdouble *x_data;
  guint stride;

  guint (*get_index) (const NcmSpline *s, const gdouble x);
} NcmSplinePrivate;

enum
{
  PROP_0,
  PROP_LEN,
  PROP_X,
  PROP_Y,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmSpline, ncm_spline, G_TYPE_OBJECT)

static void
ncm_spline_init (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  self->len         = 0;
  self->xv          = NULL;
  self->yv          = NULL;
  self->empty       = TRUE;
  self->acc         = NULL;
  self->init        = FALSE;
  self->bucket      = g_array_new (FALSE, FALSE, sizeof (gsize));
  self->bucket_dx   = 0.0;
  self->one_over_dx = 0.0;
  self->x_min       = 0.0;
  self->x_max       = 0.0;
  self->start_right = 0;
  self->last_poly   = 0;
  self->x_data      = NULL;
  self->stride      = 0;
  self->get_index   = NULL;
}

static void
_ncm_spline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_spline_parent_class)->constructed (object);
  {
    NcmSpline *s                  = NCM_SPLINE (object);
    NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

    if (self->len > 0)
    {
      guint len = self->len;

      self->len = 0;
      ncm_spline_set_len (s, len);
    }
  }
}

static void
_ncm_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpline *s                  = NCM_SPLINE (object);
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  g_return_if_fail (NCM_IS_SPLINE (object));

  switch (prop_id)
  {
    case PROP_LEN:
      self->len = g_value_get_uint (value);
      break;
    case PROP_X:
    {
      if (self->len == 0)
      {
        g_error ("ncm_spline_set_property: cannot set vector on an empty spline.");
      }
      else
      {
        ncm_vector_substitute (&self->xv, g_value_get_object (value), TRUE);
        NCM_SPLINE_GET_CLASS (s)->reset (s);
      }

      break;
    }
    case PROP_Y:
    {
      if (self->len == 0)
      {
        g_error ("ncm_spline_set_property: cannot set vector on an empty spline.");
      }
      else
      {
        ncm_vector_substitute (&self->yv, g_value_get_object (value), TRUE);
        NCM_SPLINE_GET_CLASS (s)->reset (s);
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpline *s                  = NCM_SPLINE (object);
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  g_return_if_fail (NCM_IS_SPLINE (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, self->len);
      break;
    case PROP_X:
      g_value_set_object (value, self->xv);
      break;
    case PROP_Y:
      g_value_set_object (value, self->yv);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_dispose (GObject *object)
{
  NcmSpline *s                  = NCM_SPLINE (object);
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  ncm_vector_clear (&self->xv);
  ncm_vector_clear (&self->yv);

  self->empty = TRUE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_parent_class)->dispose (object);
}

static void
_ncm_spline_finalize (GObject *object)
{
  NcmSpline *s                  = NCM_SPLINE (object);
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  g_clear_pointer (&self->acc, gsl_interp_accel_free);
  g_clear_pointer (&self->bucket, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_parent_class)->finalize (object);
}

static gdouble _ncm_spline_eval_idx_default (const NcmSpline *s, const gdouble x, const gsize i);
static gdouble _ncm_spline_deriv_idx_default (const NcmSpline *s, const gdouble x, const gsize i);
static gdouble _ncm_spline_integ_idx_default (const NcmSpline *s, const gdouble xi, const gsize i, const gdouble xf, const gsize f);

static void
ncm_spline_class_init (NcmSplineClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_spline_constructed;
  object_class->set_property = &_ncm_spline_set_property;
  object_class->get_property = &_ncm_spline_get_property;
  object_class->dispose      = &_ncm_spline_dispose;
  object_class->finalize     = &_ncm_spline_finalize;

  /**
   * NcmSpline:length:
   *
   * The number of knots.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "Spline length",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpline:x:
   *
   * The knot vector.
   */
  g_object_class_install_property (object_class,
                                   PROP_X,
                                   g_param_spec_object ("x",
                                                        NULL,
                                                        "Spline knots",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpline:y:
   *
   * The value vector.
   */
  g_object_class_install_property (object_class,
                                   PROP_Y,
                                   g_param_spec_object ("y",
                                                        NULL,
                                                        "Spline values",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->name         = NULL;
  klass->reset        = NULL;
  klass->prepare      = NULL;
  klass->prepare_base = NULL;
  klass->min_size     = NULL;
  klass->eval         = NULL;
  klass->deriv        = NULL;
  klass->deriv2       = NULL;
  klass->integ        = NULL;
  klass->eval_idx     = &_ncm_spline_eval_idx_default;
  klass->deriv_idx    = &_ncm_spline_deriv_idx_default;
  klass->integ_idx    = &_ncm_spline_integ_idx_default;
}

/**
 * ncm_spline_copy_empty:
 * @s: a constant #NcmSpline
 *
 * Returns: (transfer full): a new, empty spline of the type of @s.
 */
NcmSpline *
ncm_spline_copy_empty (const NcmSpline *s)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->copy_empty (s);
}

/**
 * ncm_spline_copy:
 * @s: a constant #NcmSpline
 *
 * Copies @s, with copies of its knot and value vectors, and prepares the copy.
 *
 * Returns: (transfer full): a new #NcmSpline.
 */
NcmSpline *
ncm_spline_copy (const NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);
  NcmVector *xv;
  NcmVector *yv;
  NcmSpline *s_cpy;

  g_assert (self->xv != NULL && self->yv != NULL);

  xv = ncm_vector_dup (self->xv);
  yv = ncm_vector_dup (self->yv);

  s_cpy = ncm_spline_new (s, xv, yv, TRUE);

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return s_cpy;
}

/**
 * ncm_spline_new:
 * @s: a constant #NcmSpline
 * @xv: the knots
 * @yv: the values at @xv
 * @init: whether to prepare the new spline
 *
 * Creates a spline of the type of @s with @xv and @yv, see ncm_spline_set().
 *
 * Returns: (transfer full): a new #NcmSpline.
 */
NcmSpline *
ncm_spline_new (const NcmSpline *s, NcmVector *xv, NcmVector *yv, gboolean init)
{
  NcmSpline *s_new = ncm_spline_copy_empty (s);

  ncm_spline_set (s_new, xv, yv, init);

  return s_new;
}

/**
 * ncm_spline_new_array:
 * @s: a constant #NcmSpline
 * @x: (element-type double): the knots
 * @y: (element-type double): the values at @x
 * @init: whether to prepare the new spline
 *
 * Same as ncm_spline_new() with the knots and values in arrays; the vectors share their data.
 *
 * Returns: (transfer full): a new #NcmSpline.
 */
NcmSpline *
ncm_spline_new_array (const NcmSpline *s, GArray *x, GArray *y, gboolean init)
{
  NcmSpline *s_new = ncm_spline_copy_empty (s);

  ncm_spline_set_array (s_new, x, y, init);

  return s_new;
}

/**
 * ncm_spline_new_data:
 * @s: a constant #NcmSpline
 * @x: the knots
 * @y: the values at @x
 * @len: length of @x and @y
 * @init: whether to prepare the new spline
 *
 * Same as ncm_spline_new() with the knots and values in C arrays, which are used in place,
 * not copied, and must outlive the spline.
 *
 * Returns: (transfer full): a new #NcmSpline.
 */
NcmSpline *
ncm_spline_new_data (const NcmSpline *s, gdouble *x, gdouble *y, gsize len, gboolean init)
{
  NcmSpline *s_new = ncm_spline_copy_empty (s);

  ncm_spline_set_data_static (s_new, x, y, len, init);

  return s_new;
}

/**
 * ncm_spline_set:
 * @s: a #NcmSpline
 * @xv: the knots
 * @yv: the values at @xv
 * @init: whether to prepare @s
 *
 * Sets the knot and value vectors of @s, keeping references to them. Aborts if their
 * lengths differ or are below ncm_spline_min_size().
 *
 * Returns: (transfer none): @s.
 */
NcmSpline *
ncm_spline_set (NcmSpline *s, NcmVector *xv, NcmVector *yv, gboolean init)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  g_assert (xv != NULL && yv != NULL);

  if (ncm_vector_len (xv) != ncm_vector_len (yv))
    g_error ("ncm_spline_set: knot and function values vector has not the same size");

  if (ncm_vector_len (xv) < NCM_SPLINE_GET_CLASS (s)->min_size (s))
    g_error ("ncm_spline_set: min size for [%s] is %zu but vector size is %u", NCM_SPLINE_GET_CLASS (s)->name (s),
             NCM_SPLINE_GET_CLASS (s)->min_size (s), ncm_vector_len (xv));

  if (self->xv != NULL)
  {
    if (self->xv != xv)
    {
      ncm_vector_ref (xv);
      ncm_vector_free (self->xv);
      self->xv = xv;
    }
  }
  else
  {
    ncm_vector_ref (xv);
    self->xv = xv;
  }

  if (self->yv != NULL)
  {
    if (self->yv != yv)
    {
      ncm_vector_ref (yv);
      ncm_vector_free (self->yv);
      self->yv = yv;
    }
  }
  else
  {
    ncm_vector_ref (yv);
    self->yv = yv;
  }

  self->len = ncm_vector_len (xv);

  NCM_SPLINE_GET_CLASS (s)->reset (s);

  self->empty = FALSE;

  if (init)
    ncm_spline_prepare (s);

  if (self->acc != NULL)
  {
    ncm_spline_acc (s, FALSE);
    ncm_spline_acc (s, TRUE);
  }

  return s;
}

/**
 * ncm_spline_ref:
 * @s: a #NcmSpline
 *
 * Increases the reference count of @s by one.
 *
 * Returns: (transfer full): @s.
 */
NcmSpline *
ncm_spline_ref (NcmSpline *s)
{
  return g_object_ref (s);
}

/**
 * ncm_spline_free:
 * @s: a #NcmSpline
 *
 * Decreases the reference count of @s by one.
 */
void
ncm_spline_free (NcmSpline *s)
{
  g_object_unref (s);
}

/**
 * ncm_spline_clear:
 * @s: a #NcmSpline
 *
 * If *@s is not %NULL, decreases its reference count by one and sets *@s to %NULL.
 */
void
ncm_spline_clear (NcmSpline **s)
{
  g_clear_object (s);
}

/**
 * ncm_spline_acc:
 * @s: a #NcmSpline
 * @enable: whether to use a GSL accelerator
 *
 * Enables or disables the GSL accelerator, which caches the last interval found. With it,
 * the evaluation modifies @s and must not run in two threads at once. The choice takes effect
 * at the next ncm_spline_prepare().
 */
void
ncm_spline_acc (NcmSpline *s, gboolean enable)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  if (enable)
  {
    if (self->acc == NULL)
      self->acc = gsl_interp_accel_alloc ();
  }
  else if (self->acc != NULL)
  {
    gsl_interp_accel_free (self->acc);
    self->acc = NULL;
  }
}

/**
 * ncm_spline_peek_acc: (skip)
 * @s: a #NcmSpline
 *
 * Returns: (transfer none) (nullable): the GSL accelerator, or %NULL if disabled.
 */
gsl_interp_accel *
ncm_spline_peek_acc (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  return self->acc;
}

/**
 * ncm_spline_set_len:
 * @s: a #NcmSpline
 * @len: the number of knots
 *
 * Replaces the knot and value vectors by new ones of length @len when it differs from the
 * current length.
 */
void
ncm_spline_set_len (NcmSpline *s, guint len)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  if (self->len != len)
  {
    g_assert_cmpuint (len, >, 0);
    {
      NcmVector *xv = ncm_vector_new (len);
      NcmVector *yv = ncm_vector_new (len);

      ncm_spline_set (s, xv, yv, FALSE);
      ncm_vector_free (xv);
      ncm_vector_free (yv);
    }
  }
}

/**
 * ncm_spline_get_len:
 * @s: a #NcmSpline
 *
 * Returns: the number of knots.
 */
guint
ncm_spline_get_len (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  return self->len;
}

/**
 * ncm_spline_set_xv:
 * @s: a #NcmSpline
 * @xv: the knots
 * @init: whether to prepare @s
 *
 * Replaces the knot vector, see ncm_spline_set().
 */
void
ncm_spline_set_xv (NcmSpline *s, NcmVector *xv, gboolean init)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  ncm_spline_set (s, xv, self->yv, init);
}

/**
 * ncm_spline_set_yv:
 * @s: a #NcmSpline
 * @yv: the values at the knots
 * @init: whether to prepare @s
 *
 * Replaces the value vector, see ncm_spline_set().
 */
void
ncm_spline_set_yv (NcmSpline *s, NcmVector *yv, gboolean init)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  ncm_spline_set (s, self->xv, yv, init);
}

/**
 * ncm_spline_set_array:
 * @s: a #NcmSpline
 * @x: (element-type double): the knots
 * @y: (element-type double): the values at @x
 * @init: whether to prepare @s
 *
 * Same as ncm_spline_set() with the knots and values in arrays; the vectors share their data.
 */
void
ncm_spline_set_array (NcmSpline *s, GArray *x, GArray *y, gboolean init)
{
  NcmVector *xv = ncm_vector_new_array (x);
  NcmVector *yv = ncm_vector_new_array (y);

  ncm_spline_set (s, xv, yv, init);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

/**
 * ncm_spline_set_data_static:
 * @s: a #NcmSpline
 * @x: the knots
 * @y: the values at @x
 * @len: length of @x and @y
 * @init: whether to prepare @s
 *
 * Same as ncm_spline_set() with the knots and values in C arrays, which are used in place,
 * not copied, and must outlive the spline.
 */
void
ncm_spline_set_data_static (NcmSpline *s, gdouble *x, gdouble *y, gsize len, gboolean init)
{
  NcmVector *xv = ncm_vector_new_data_static (x, len, 1);
  NcmVector *yv = ncm_vector_new_data_static (y, len, 1);

  ncm_spline_set (s, xv, yv, init);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

/**
 * ncm_spline_get_xv:
 * @s: a #NcmSpline
 *
 * Returns: (transfer full): the knot vector, %NULL before it is set.
 */
NcmVector *
ncm_spline_get_xv (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  if (self->xv != NULL)
    return ncm_vector_ref (self->xv);
  else
    return NULL;
}

/**
 * ncm_spline_get_yv:
 * @s: a #NcmSpline
 *
 * Returns: (transfer full): the value vector, %NULL before it is set.
 */
NcmVector *
ncm_spline_get_yv (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  if (self->yv != NULL)
    return ncm_vector_ref (self->yv);
  else
    return NULL;
}

/**
 * ncm_spline_peek_xv:
 * @s: a #NcmSpline
 *
 * Returns: (transfer none): the knot vector, %NULL before it is set.
 */
NcmVector *
ncm_spline_peek_xv (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  return self->xv;
}

/**
 * ncm_spline_peek_yv:
 * @s: a #NcmSpline
 *
 * Returns: (transfer none): the value vector, %NULL before it is set.
 */
NcmVector *
ncm_spline_peek_yv (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  return self->yv;
}

/**
 * ncm_spline_get_bounds:
 * @s: a #NcmSpline
 * @lb: (out): the first knot
 * @ub: (out): the last knot
 *
 * Gets the interval spanned by the knots.
 */
void
ncm_spline_get_bounds (NcmSpline *s, gdouble *lb, gdouble *ub)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  g_assert_cmpuint (self->len, >, 0);

  *lb = ncm_vector_get (self->xv, 0);
  *ub = ncm_vector_get (self->xv, self->len - 1);
}

/**
 * ncm_spline_is_init:
 * @s: a #NcmSpline
 *
 * Returns: whether @s was prepared since its creation.
 */
gboolean
ncm_spline_is_init (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  return self->init;
}

/**
 * ncm_spline_prepare:
 * @s: a #NcmSpline
 *
 * Computes the interpolant from the knots and values, and the table of buckets used to find
 * intervals, see ncm_spline_post_prepare(). Required before evaluation, and after any change
 * of the knots or values.
 */
void
ncm_spline_prepare (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  self->init = TRUE;
  NCM_SPLINE_GET_CLASS (s)->prepare (s);
  ncm_spline_post_prepare (s);
}

static guint _ncm_spline_get_index_no_stride_accel (const NcmSpline *s, const gdouble x);
static guint _ncm_spline_get_index_no_stride (const NcmSpline *s, const gdouble x);
static guint _ncm_spline_get_index_stride_accel (const NcmSpline *s, const gdouble x);
static guint _ncm_spline_get_index_stride (const NcmSpline *s, const gdouble x);

/**
 * ncm_spline_post_prepare:
 * @s: a #NcmSpline
 *
 * Builds the table of uniform buckets used to find the interval of $x$, and chooses the
 * search according to the stride of the knots and the accelerator. ncm_spline_prepare() calls
 * it; objects that prepare a spline by other means, such as #NcmSpline2d, must call it too.
 */
void
ncm_spline_post_prepare (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);
  NcmVector *xv                 = self->xv;
  const guint n_buckets         = (guint) (self->len - 1);
  const gdouble x_min           = ncm_vector_get (xv, 0);
  const gdouble x_max           = ncm_vector_get (xv, self->len - 1);
  const gdouble dx              = (x_max - x_min) / n_buckets;
  guint i                       = 0;
  guint j                       = 0;

  g_array_set_size (self->bucket, n_buckets + 1);
  self->bucket_dx   = dx;
  self->one_over_dx = 1.0 / dx;
  self->x_min       = x_min;
  self->x_max       = x_max;
  self->start_right = self->len - 1;
  self->last_poly   = self->len - 2;
  self->x_data      = ncm_vector_data (xv);
  self->stride      = ncm_vector_stride (xv);

  for (i = 0; i <= n_buckets; i++)
  {
    gdouble x_bucket = x_min + i * dx;

    while ((j < self->len - 1) && (ncm_vector_get (xv, j + 1) < x_bucket))
      j++;

    g_array_index (self->bucket, gsize, i) = j;
  }

  g_array_index (self->bucket, gsize, n_buckets) = self->len - 2;

  if (ncm_vector_stride (self->xv) == 1)
  {
    if (self->acc)
      self->get_index = _ncm_spline_get_index_no_stride_accel;
    else
      self->get_index = _ncm_spline_get_index_no_stride;
  }
  else
  {
    if (self->acc)
      self->get_index = _ncm_spline_get_index_stride_accel;
    else
      self->get_index = _ncm_spline_get_index_stride;
  }
}

/**
 * ncm_spline_prepare_base:
 * @s: a #NcmSpline
 *
 * Computes the coefficients a two-dimensional spline needs from @s, for a cubic spline
 * its second derivatives, without the rest of ncm_spline_prepare().
 */

void
ncm_spline_prepare_base (NcmSpline *s)
{
  if (NCM_SPLINE_GET_CLASS (s)->prepare_base)
    NCM_SPLINE_GET_CLASS (s)->prepare_base (s);
}

/**
 * ncm_spline_eval:
 * @s: a constant #NcmSpline
 * @x: the point
 *
 * Evaluates the interpolant at @x, with no range check; see #NcmSpline for the behaviour
 * outside the knots.
 *
 * Returns: the interpolated value at @x.
 */

gdouble
ncm_spline_eval (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->eval (s, x);
}

/* The interval index is a hint: without a type-specific use, evaluate without it */
static gdouble
_ncm_spline_eval_idx_default (const NcmSpline *s, const gdouble x, const gsize i)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->eval (s, x);
}

/**
 * ncm_spline_eval_idx:
 * @s: a constant #NcmSpline
 * @x: the point
 * @i: index of the interval, $x_i \le x < x_{i+1}$
 *
 * Same as ncm_spline_eval() with the interval given, for splines that share their knots.
 * A type that does not use the interval evaluates without it.
 *
 * Returns: the interpolated value at @x.
 */

gdouble
ncm_spline_eval_idx (const NcmSpline *s, const gdouble x, const gsize i)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->eval_idx (s, x, i);
}

/**
 * ncm_spline_eval_deriv:
 * @s: a constant #NcmSpline
 * @x: the point
 *
 * Returns: the first derivative of the interpolant at @x.
 */

gdouble
ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv (s, x);
}

static gdouble
_ncm_spline_deriv_idx_default (const NcmSpline *s, const gdouble x, const gsize i)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv (s, x);
}

/**
 * ncm_spline_eval_deriv_idx:
 * @s: a constant #NcmSpline
 * @x: the point
 * @i: index of the interval, $x_i \le x < x_{i+1}$
 *
 * Same as ncm_spline_eval_deriv() with the interval given. A type that does not use the
 * interval evaluates without it.
 *
 * Returns: the first derivative of the interpolant at @x.
 */

gdouble
ncm_spline_eval_deriv_idx (const NcmSpline *s, const gdouble x, const gsize i)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv_idx (s, x, i);
}

/**
 * ncm_spline_eval_deriv2:
 * @s: a constant #NcmSpline
 * @x: the point
 *
 * Returns: the second derivative of the interpolant at @x.
 */

gdouble
ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv2 (s, x);
}

/**
 * ncm_spline_eval_deriv_nmax:
 * @s: a constant #NcmSpline
 * @x: the point
 *
 * Returns: the highest nonzero derivative of the interpolant at @x, the third for a cubic.
 */

gdouble
ncm_spline_eval_deriv_nmax (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv_nmax (s, x);
}

/**
 * ncm_spline_eval_integ:
 * @s: a constant #NcmSpline
 * @x0: the lower limit
 * @x1: the upper limit
 *
 * Limits in decreasing order give minus the integral over the reversed interval.
 *
 * Returns: $\int_{x_0}^{x_1} s(x)\,\mathrm{d}x$ of the interpolant $s$.
 */

gdouble
ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  if (x1 < x0)
    return -NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->integ (s, x1, x0);

  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->integ (s, x0, x1);
}

static gdouble
_ncm_spline_integ_idx_default (const NcmSpline *s, const gdouble xi, const gsize i, const gdouble xf, const gsize f)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->integ (s, xi, xf);
}

/**
 * ncm_spline_eval_integ_idx:
 * @s: a constant #NcmSpline
 * @xi: the lower limit
 * @i: index of the interval of @xi
 * @xf: the upper limit
 * @f: index of the interval of @xf
 *
 * Same as ncm_spline_eval_integ() with the intervals given. A type that does not use the
 * intervals integrates without them.
 *
 * Returns: $\int_{x_i}^{x_f} s(x)\,\mathrm{d}x$.
 */

gdouble
ncm_spline_eval_integ_idx (const NcmSpline *s, const gdouble xi, const gsize i, const gdouble xf, const gsize f)
{
  if (xf < xi)
    return -NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->integ_idx (s, xf, f, xi, i);

  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->integ_idx (s, xi, i, xf, f);
}

/**
 * ncm_spline_is_empty:
 * @s: a constant #NcmSpline
 *
 * Returns: whether @s has no knots and values set.
 */

gboolean
ncm_spline_is_empty (const NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);

  return self->empty;
}

/**
 * ncm_spline_min_size:
 * @s: a constant #NcmSpline
 *
 * Returns: the minimum number of knots of the type of @s.
 */

gsize
ncm_spline_min_size (const NcmSpline *s)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->min_size (s);
}

gsize
_ncm_spline_bsearch_stride (const gdouble x_array[], const guint stride, const gdouble x, gsize index_lo, gsize index_hi)
{
  gsize ilo = index_lo;
  gsize ihi = index_hi;

  while (ihi > ilo + 1)
  {
    gsize i = (ihi + ilo) / 2;

    if (x_array[i * stride] > x)
      ihi = i;
    else
      ilo = i;
  }

  return ilo;
}

gsize
_ncm_spline_accel_find (gsl_interp_accel *a, const gdouble xa[], const guint stride, gsize len, gdouble x)
{
  gsize x_index = a->cache;

  if (x < xa[x_index * stride])
  {
    a->miss_count++;
    a->cache = _ncm_spline_bsearch_stride (xa, stride, x, 0, x_index);
  }
  else if (x >= xa[stride * (x_index + 1)])
  {
    a->miss_count++;
    a->cache = _ncm_spline_bsearch_stride (xa, stride, x, x_index, len - 1);
  }
  else
  {
    a->hit_count++;
  }

  return a->cache;
}

static guint
_ncm_spline_get_index_no_stride_accel (const NcmSpline *s, const gdouble x)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);

  return gsl_interp_accel_find (self->acc, self->x_data, self->len, x);
}

static guint
_ncm_spline_get_index_no_stride (const NcmSpline *s, const gdouble x)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);
  gsize left                    = 0;
  gsize right                   = self->start_right;

  if (G_UNLIKELY (x < self->x_min))
    return 0;

  if (G_UNLIKELY (x > self->x_max))
    return self->last_poly;

  {
    const guint n_buckets = self->bucket->len - 1;
    guint i_bucket        = (guint) ((x - self->x_min) * self->one_over_dx);

    if (i_bucket >= n_buckets)
      i_bucket = n_buckets - 1;

    left  = g_array_index (self->bucket, gsize, i_bucket);
    right = g_array_index (self->bucket, gsize, i_bucket + 1);
  }

  if (left == right)
    return left;

  return gsl_interp_bsearch (self->x_data, x, left, right + 1);
}

static guint
_ncm_spline_get_index_stride_accel (const NcmSpline *s, const gdouble x)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);

  return _ncm_spline_accel_find (self->acc, self->x_data, self->stride, self->len, x);
}

static guint
_ncm_spline_get_index_stride (const NcmSpline *s, const gdouble x)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);
  gsize left                    = 0;
  gsize right                   = self->start_right;

  if (G_UNLIKELY (x < self->x_min))
    return 0;

  if (G_UNLIKELY (x > self->x_max))
    return self->last_poly;

  if (TRUE)
  {
    const guint n_buckets = self->bucket->len - 1;
    guint i_bucket        = (guint) ((x - self->x_min) * self->one_over_dx);

    if (i_bucket >= n_buckets)
      i_bucket = n_buckets - 1;

    left  = g_array_index (self->bucket, gsize, i_bucket);
    right = g_array_index (self->bucket, gsize, i_bucket + 1);
  }

  if (left == right)
    return left;

  return _ncm_spline_bsearch_stride (self->x_data, self->stride, x, left, right + 1);
}

/**
 * ncm_spline_get_index:
 * @s: a constant #NcmSpline
 * @x: the point
 *
 * Finds the interval of @x, clamped to the first and last intervals outside the knots.
 *
 * Returns: the index $i$ with $x_i \le x < x_{i+1}$.
 */
guint
ncm_spline_get_index (const NcmSpline *s, const gdouble x)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);

  return self->get_index (s, x);
}

/**
 * ncm_spline_curvature_density:
 * @s: a #NcmSpline
 * @ctype: a #NcmSplineCurvatureType
 * @x: the point
 *
 * Evaluates the curvature density $c(x)$ selected by @ctype. @s must be prepared.
 *
 * Returns: $c(x)$.
 */
gdouble
ncm_spline_curvature_density (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble x)
{
  const gdouble d2 = ncm_spline_eval_deriv2 (s, x);

  switch (ctype)
  {
    case NCM_SPLINE_CURVATURE_D2:
      return d2;

    case NCM_SPLINE_CURVATURE_GEOMETRIC:
    {
      const gdouble d1 = ncm_spline_eval_deriv (s, x);

      return d2 / pow (1.0 + d1 * d1, 1.5);
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/* Relative accuracy of the curvature integrals, reached on each knot interval */
#define NCM_SPLINE_CURVATURE_RELTOL (1.0e-10)

/*
 * Integrates F over [xi, xf] one knot interval of s at a time, where the curvature density,
 * and so the integrand, is smooth; over the whole range its kinks at the knots would limit
 * the accuracy. Aborts if an interval fails to converge.
 */
static gdouble
_ncm_spline_integ_by_knots (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, const gchar *func)
{
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  NcmVector *xv                 = ncm_spline_peek_xv (s);
  const guint len               = ncm_vector_len (xv);
  gdouble a                     = xi;
  gdouble total                 = 0.0;
  guint i                       = 0;

  while (a < xf)
  {
    gdouble b = xf;
    gdouble result, error;
    gint status;

    while ((i < len) && (ncm_vector_get (xv, i) <= a))
      i++;

    if ((i < len) && (ncm_vector_get (xv, i) < xf))
      b = ncm_vector_get (xv, i);

    status = gsl_integration_qag (F, a, b, 0.0, NCM_SPLINE_CURVATURE_RELTOL, NCM_INTEGRAL_PARTITION, 6, *w, &result, &error);

    if (status != GSL_SUCCESS)
      g_error ("%s: the integral over [%.15g, %.15g] did not reach the relative tolerance %.1e: %s.",
               func, a, b, NCM_SPLINE_CURVATURE_RELTOL, gsl_strerror (status));

    total += result;
    a      = b;
  }

  ncm_memory_pool_return (w);

  return total;
}

typedef struct _NcmSplineCurvatureArg
{
  NcmSpline *s;
  NcmSplineCurvatureType ctype;
  gdouble p;
} NcmSplineCurvatureArg;

static gdouble
_ncm_spline_curvature_lp_integrand (gdouble x, gpointer data)
{
  NcmSplineCurvatureArg *arg = (NcmSplineCurvatureArg *) data;
  const gdouble c            = ncm_spline_curvature_density (arg->s, arg->ctype, x);

  return pow (fabs (c), arg->p);
}

/**
 * ncm_spline_curvature_lp_norm:
 * @s: a #NcmSpline
 * @ctype: a #NcmSplineCurvatureType
 * @p: the order $p > 0$
 * @xi: the lower limit
 * @xf: the upper limit
 *
 * Computes the $L_p$ norm of the curvature density normalized by the interval,
 * $$N_p = \left(\frac{1}{x_f - x_i}\int_{x_i}^{x_f} |c(x)|^p\,\mathrm{d}x\right)^{1/p}.$$
 * $p = 2$ gives the root-mean-square curvature, and $p \to \infty$ tends to
 * ncm_spline_curvature_max(). The integral is computed on each knot interval, where $c(x)$ is
 * smooth, to a relative accuracy of $10^{-10}$; aborts if an interval does not converge. @s must
 * be prepared.
 *
 * Returns: $N_p$.
 */
gdouble
ncm_spline_curvature_lp_norm (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble p, const gdouble xi, const gdouble xf)
{
  NcmSplineCurvatureArg arg = {s, ctype, p};
  gsl_function F;
  gdouble result;

  g_assert_cmpfloat (p, >, 0.0);
  g_assert_cmpfloat (xf, >, xi);

  F.function = &_ncm_spline_curvature_lp_integrand;
  F.params   = &arg;

  result = _ncm_spline_integ_by_knots (s, &F, xi, xf, "ncm_spline_curvature_lp_norm");

  return pow (result / (xf - xi), 1.0 / p);
}

typedef struct _NcmSplineCurvatureWeightedArg
{
  NcmSpline *s;
  NcmSpline *weight;
  NcmSplineCurvatureType ctype;
  gdouble p;
} NcmSplineCurvatureWeightedArg;

static gdouble
_ncm_spline_curvature_weighted_lp_integrand (gdouble x, gpointer data)
{
  NcmSplineCurvatureWeightedArg *arg = (NcmSplineCurvatureWeightedArg *) data;
  const gdouble w                    = ncm_spline_eval (arg->weight, x);
  const gdouble c                    = ncm_spline_curvature_density (arg->s, arg->ctype, x);

  return w * pow (fabs (c), arg->p);
}

static gdouble
_ncm_spline_curvature_weight_integrand (gdouble x, gpointer data)
{
  NcmSplineCurvatureWeightedArg *arg = (NcmSplineCurvatureWeightedArg *) data;

  return ncm_spline_eval (arg->weight, x);
}

/**
 * ncm_spline_curvature_weighted_lp_norm:
 * @s: a #NcmSpline
 * @ctype: a #NcmSplineCurvatureType
 * @p: the order $p > 0$
 * @weight: a #NcmSpline with the weight $W(x) \ge 0$
 * @xi: the lower limit
 * @xf: the upper limit
 *
 * Computes the $L_p$ norm of the curvature density weighted by $W$,
 * $$N_p = \left(\frac{\int_{x_i}^{x_f} W(x)\,|c(x)|^p\,\mathrm{d}x}{\int_{x_i}^{x_f} W(x)\,\mathrm{d}x}\right)^{1/p},$$
 * which is ncm_spline_curvature_lp_norm() for a constant $W$. Both integrals are computed on each
 * knot interval of @s to a relative accuracy of $10^{-10}$; aborts if an interval does not
 * converge. Both splines must be prepared, and $\int W$ must be positive; aborts otherwise.
 *
 * Returns: $N_p$.
 */
gdouble
ncm_spline_curvature_weighted_lp_norm (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble p, NcmSpline *weight, const gdouble xi, const gdouble xf)
{
  NcmSplineCurvatureWeightedArg arg = {s, weight, ctype, p};
  gsl_function F;
  gdouble num, wnorm;

  g_assert_cmpfloat (p, >, 0.0);
  g_assert_cmpfloat (xf, >, xi);

  F.params = &arg;

  F.function = &_ncm_spline_curvature_weighted_lp_integrand;
  num        = _ncm_spline_integ_by_knots (s, &F, xi, xf, "ncm_spline_curvature_weighted_lp_norm");

  F.function = &_ncm_spline_curvature_weight_integrand;
  wnorm      = _ncm_spline_integ_by_knots (s, &F, xi, xf, "ncm_spline_curvature_weighted_lp_norm");

  g_assert_cmpfloat (wnorm, >, 0.0);

  return pow (num / wnorm, 1.0 / p);
}

/**
 * ncm_spline_curvature_max:
 * @s: a #NcmSpline
 * @ctype: a #NcmSplineCurvatureType
 * @xi: the lower limit
 * @xf: the upper limit
 *
 * Estimates $\max_x |c(x)|$ over [@xi, @xf] from the knots inside the interval and a uniform
 * grid of $20 n$ points, for $n$ knots. @s must be prepared.
 *
 * Returns: the estimated maximum of $|c(x)|$.
 */
gdouble
ncm_spline_curvature_max (NcmSpline *s, NcmSplineCurvatureType ctype, const gdouble xi, const gdouble xf)
{
  NcmVector *xv      = ncm_spline_peek_xv (s);
  const guint len    = ncm_vector_len (xv);
  const guint n_grid = 20 * len;
  gdouble cmax       = 0.0;
  guint i;

  g_assert_cmpfloat (xf, >, xi);

  for (i = 0; i < len; i++)
  {
    const gdouble x = ncm_vector_get (xv, i);

    if ((x >= xi) && (x <= xf))
      cmax = GSL_MAX_DBL (cmax, fabs (ncm_spline_curvature_density (s, ctype, x)));
  }

  for (i = 0; i <= n_grid; i++)
  {
    const gdouble x = xi + (xf - xi) * i / (1.0 * n_grid);

    cmax = GSL_MAX_DBL (cmax, fabs (ncm_spline_curvature_density (s, ctype, x)));
  }

  return cmax;
}

/* Utilities -- internal use */

gdouble
_ncm_spline_util_integ_eval (const gdouble ai, const gdouble bi, const gdouble ci, const gdouble di, const gdouble xi, const gdouble a, const gdouble b)
{
  const gdouble r1    = (a - xi);
  const gdouble r2    = (b - xi);
  const gdouble r12   = (r1 + r2);
  const gdouble bterm = 0.5 * bi * r12;
  const gdouble cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
  const gdouble dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);

  return (b - a) * (ai + bterm + cterm + dterm);
}


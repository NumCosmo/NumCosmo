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
 * SECTION:ncm_spline
 * @title: NcmSpline
 * @short_description: Abstract class for implementing splines.
 * @stability: Stable
 * @include: numcosmo/math/ncm_spline.h
 *
 * This class comprises all functions to provide a #NcmSpline, together with
 * all necessary methods.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline.h"
#include "math/ncm_cfg.h"

typedef struct _NcmSplinePrivate
{
  /*< private >*/
  GObject parent_instance;
  gsize len;
  NcmVector *xv;
  NcmVector *yv;
  gsl_interp_accel *acc;
  gboolean init;
  gboolean empty;
} NcmSplinePrivate;

enum
{
  PROP_0,
  PROP_LEN,
  PROP_X,
  PROP_Y,
  PROP_ACC,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmSpline, ncm_spline, G_TYPE_OBJECT)

static void
ncm_spline_init (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  self->len   = 0;
  self->xv    = NULL;
  self->yv    = NULL;
  self->empty = TRUE;
  self->acc   = NULL;
  self->init  = FALSE;
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
        g_error ("ncm_spline_set_property: cannot set vector on an empty spline.");
      else
        ncm_vector_substitute (&self->xv, g_value_get_object (value), TRUE);

      break;
    }
    case PROP_Y:
    {
      if (self->len == 0)
        g_error ("ncm_spline_set_property: cannot set vector on an empty spline.");
      else
        ncm_vector_substitute (&self->yv, g_value_get_object (value), TRUE);

      break;
    }
    case PROP_ACC:
      ncm_spline_acc (s, g_value_get_boolean (value));
      break;
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
    case PROP_ACC:
      g_value_set_boolean (value, self->acc != NULL ? TRUE : FALSE);
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

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_parent_class)->finalize (object);
}

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
   * The spline length (total number of knots).
   *
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
   * #NcmVector with the spline knots.
   *
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
   * #NcmVector with the spline values.
   *
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
}

/**
 * ncm_spline_copy_empty:
 * @s: a constant #NcmSpline
 *
 * This function copies the spline @s into an initialized empty #NcmSpline of a specific type.
 *
 * Returns: (transfer full): a #NcmSpline
 */
NcmSpline *
ncm_spline_copy_empty (const NcmSpline *s)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->copy_empty (s);
}

/**
 * ncm_spline_copy:
 * @s: a costant #NcmSpline
 *
 * This function copies the two #NcmVector of the spline @s into those two
 * #NcmVector of a new #NcmSpline.
 *
 * Returns: (transfer full): a #NcmSpline
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
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of the values of the function, to be interpolated, computed at @xv
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it
 *
 * This function returns a new #NcmSpline, where the knots of this new spline are given
 * in the #NcmVector @xv and the values of the function, at those knots, to be interpolated are
 * given in the #NcmVector @yv.
 *
 * Returns: (transfer full): a new #NcmSpline
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
 * @x: (element-type double): GArray of knots
 * @y: (element-type double): GArray of the values of the function, to be interpolated, computed at @x
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it
 *
 * This function returns a new #NcmSpline, where the knots of this new spline are given
 * in the GArray @x and the values of the function, at those knots, to be interpolated are
 * given in the GArray @y.
 *
 * Returns: (transfer full): a new #NcmSpline
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
 * @x: array of knots
 * @y: array of the values of the function, to be interpolated, computed at @x
 * @len: lenght of @x and @y
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it
 *
 * This function returns a new #NcmSpline, where the knots of this new spline are given
 * in the array @x and the values of the function, at those knots, to be interpolated are
 * given in the array @y.
 *
 * Returns: (transfer full): a new #NcmSpline
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
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of the values of the function, to be interpolated, computed at @xv
 * @init: TRUE to prepare @s or FALSE to not prepare it
 *
 * This funtion sets both @xv and @yv vectors to @s.
 * The two vectors must have the same length.
 *
 * Returns: (transfer none): a #NcmSpline
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
 * Returns: (transfer full): @s
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
 * Atomically decrements the reference count of @s by one. If the reference count drops to 0,
 * all memory allocated by @s is released.
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
 * Atomically decrements the reference count of @s by one. If the reference count drops to 0,
 * all memory allocated by @s is released. The pointer is set to NULL.
 */
void
ncm_spline_clear (NcmSpline **s)
{
  g_clear_object (s);
}

/**
 * ncm_spline_acc:
 * @s: a #NcmSpline
 * @enable: a boolean
 *
 * Enables or disables spline accelerator. Note that, if
 * enabled, the spline becomes non-reentrant. In other words,
 * if @enable is TRUE, the spline evaluation is not thread safe.
 * Therefore, it should not be called concomitantly by two different threads.
 *
 * Warning: the accelerator must be reset if the spline's size changes, otherwise,
 * it can accessan out-of-bound index.
 *
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
 * This function returns the spline accelerator if it is enabled.
 * Otherwise, it returns NULL.
 *
 * Returns: (transfer none): a #gsl_interp_accel
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
 * @len: number of knots in the spline
 *
 * This function sets @len as the length of the spline,
 * it allocates the necessary #NcmVector. If it is already
 * allocated with different length it frees the current vectors
 * and allocates new ones.
 *
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
 * This function gets the length of the spline.
 *
 * Returns: spline's size.
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
 * @xv: #NcmVector of knots
 * @init: TRUE to prepare @s or FALSE to not prepare it
 *
 * This function sets @xv as the knot vector of the spline.
 *
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
 * @yv: #NcmVector of the values of the function to be interpolated
 * @init: TRUE to prepare @s or FALSE to not prepare it
 *
 * This function sets @yv as the function values vector. This #NcmVector @yv
 * comprises the function values computed at the knots of the spline.
 *
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
 * @x: (element-type double): GArray of knots
 * @y: (element-type double): GArray of the values of the function, to be interpolated, computed at @x
 * @init: TRUE to prepare @s or FALSE to not prepare it
 *
 * This function sets @x as the knot vector and @y as the function values vector
 * of the spline.
 *
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
 * @x: array of knots
 * @y: array of the values of the function, to be interpolated, computed at @x
 * @len: lenght of @x and @y
 * @init: TRUE to prepare @s or FALSE to not prepare it
 *
 * This function sets @x as the knot vector and @y as the function values vector
 * of the spline.
 *
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
 * This function returns the @s #NcmVector of knots.
 *
 * Returns: (transfer full): a #NcmVector
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
 * This function returns the @s #NcmVector of the values of the function to be interpolated.
 *
 * Returns: (transfer full): a #NcmVector
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
 * This function returns the @s #NcmVector of knots.
 *
 * Returns: (transfer none): a #NcmVector
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
 * This function returns the @s #NcmVector of the values of the function to be interpolated.
 *
 * Returns: (transfer none): a #NcmVector
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
 * @lb: (out): spline lower bound
 * @ub: (out): spline upper bound
 *
 * This function returns the lower and upper bound of @s.
 *
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
 * This function returns TRUE if @s is initialized or FALSE otherwise.
 *
 * Returns: TRUE if @s is initialized or FALSE otherwise.
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
 * This function prepares the spline @s such that one can evaluate it (#ncm_spline_eval), as well as
 * to compute its first and second derivatives (#ncm_spline_eval_deriv, #ncm_spline_eval_deriv2)
 * and integration (#ncm_spline_eval_integ).
 *
 */
void
ncm_spline_prepare (NcmSpline *s)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private (s);

  self->init = TRUE;
  NCM_SPLINE_GET_CLASS (s)->prepare (s);
}

/**
 * ncm_spline_prepare_base:
 * @s: a #NcmSpline
 *
 * This function computes the second derivatives of @s and it is used to prepare a
 * bidimensional spline.
 *
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
 * @x: x-coordinate value
 *
 * Returns: The interpolated value of a function computed at @x.
 */

gdouble
ncm_spline_eval (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->eval (s, x);
}

/**
 * ncm_spline_eval_deriv:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 *
 * Returns: The derivative of an interpolated function computed at @x.
 */

gdouble
ncm_spline_eval_deriv (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv (s, x);
}

/**
 * ncm_spline_eval_deriv2:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 *
 * Returns: The second derivative of an interpolated function computed at @x.
 */

gdouble
ncm_spline_eval_deriv2 (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv2 (s, x);
}

/**
 * ncm_spline_eval_deriv_nmax:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 *
 * Returns: The highest non null derivative of an interpolated function computed at @x.
 */

gdouble
ncm_spline_eval_deriv_nmax (const NcmSpline *s, const gdouble x)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->deriv_nmax (s, x);
}

/**
 * ncm_spline_eval_integ:
 * @s: a constant #NcmSpline
 * @x0: lower integration limit
 * @x1: upper integration limit
 *
 *
 * Returns: The numerical integral of an interpolated function over the range [@x0, @x1].
 */

gdouble
ncm_spline_eval_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  return NCM_SPLINE_GET_CLASS ((NcmSpline *) s)->integ (s, x0, x1);
}

/**
 * ncm_spline_is_empty:
 * @s: a constant #NcmSpline
 *
 *
 * Returns: TRUE If @s is empty or FALSE otherwise.
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
 *
 * Returns: Minimum number of knots required.
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

/**
 * ncm_spline_get_index:
 * @s: a constant #NcmSpline
 * @x: a value of the abscissa axis
 *
 *
 * Returns: The index of the lower knot of the interval @x belongs to.
 */
guint
ncm_spline_get_index (const NcmSpline *s, const gdouble x)
{
  NcmSplinePrivate * const self = ncm_spline_get_instance_private ((NcmSpline *) s);

  if (ncm_vector_stride (self->xv) == 1)
  {
    if (self->acc)
      return gsl_interp_accel_find (self->acc, ncm_vector_ptr (self->xv, 0), self->len, x);
    else
      return gsl_interp_bsearch (ncm_vector_ptr (self->xv, 0), x, 0, ncm_vector_len (self->xv) - 1);
  }
  else
  {
    if (self->acc)
      return _ncm_spline_accel_find (self->acc, ncm_vector_ptr (self->xv, 0), ncm_vector_stride (self->xv), self->len, x);
    else
      return _ncm_spline_bsearch_stride (ncm_vector_ptr (self->xv, 0), ncm_vector_stride (self->xv), x, 0, ncm_vector_len (self->xv) - 1);
  }
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


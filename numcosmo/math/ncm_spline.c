/***************************************************************************
 *            ncm_spline.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

enum
{
  PROP_0,
  PROP_LEN,
  PROP_X,
  PROP_Y,
  PROP_ACC,
};

G_DEFINE_ABSTRACT_TYPE (NcmSpline, ncm_spline, G_TYPE_OBJECT);

static void
ncm_spline_init (NcmSpline *s)
{
  s->len   = 0;
  s->xv    = NULL;
  s->yv    = NULL;
  s->empty = TRUE;
  s->acc   = NULL;
}

static void
_ncm_spline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_spline_parent_class)->constructed (object);
  {
    NcmSpline *s = NCM_SPLINE (object);
    if (s->len > 0)
    {
      guint len = s->len;
      s->len = 0;
      ncm_spline_set_len (s, len);
    }
  }
}

static void
_ncm_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpline *s = NCM_SPLINE (object);
  g_return_if_fail (NCM_IS_SPLINE (object));

  switch (prop_id)
  {
    case PROP_LEN:
      s->len = g_value_get_uint (value);
      break;
    case PROP_X:
    {
      if (s->len == 0)
        g_error ("ncm_spline_set_property: cannot set vector on an empty spline.");
      else
        ncm_vector_substitute (&s->xv, g_value_get_object (value), TRUE);
      break;
    }
    case PROP_Y:
    {
      if (s->len == 0)
        g_error ("ncm_spline_set_property: cannot set vector on an empty spline.");
      else
        ncm_vector_substitute (&s->yv, g_value_get_object (value), TRUE);
      break;
    }
    case PROP_ACC:
      ncm_spline_acc (s, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpline *s = NCM_SPLINE (object);
  g_return_if_fail (NCM_IS_SPLINE (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, s->len);
      break;
    case PROP_X:
      g_value_set_object (value, s->xv);
      break;
    case PROP_Y:
      g_value_set_object (value, s->yv);
      break;
    case PROP_ACC:
      g_value_set_boolean (value, s->acc != NULL ? TRUE : FALSE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_dispose (GObject *object)
{
	NcmSpline *s = NCM_SPLINE (object);

  ncm_vector_clear (&s->xv);
  ncm_vector_clear (&s->yv);
  
	s->empty = TRUE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_parent_class)->dispose (object);
}

static void
_ncm_spline_finalize (GObject *object)
{
  NcmSpline *s = NCM_SPLINE (object);

  g_clear_pointer (&s->acc, gsl_interp_accel_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_parent_class)->finalize (object);
}

static void
ncm_spline_class_init (NcmSplineClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_spline_constructed;
  object_class->set_property = &_ncm_spline_set_property;
  object_class->get_property = &_ncm_spline_get_property;
  object_class->dispose      = &_ncm_spline_dispose;
  object_class->finalize     = &_ncm_spline_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "Spline length",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_X,
                                   g_param_spec_object ("x",
                                                        NULL,
                                                        "Spline knots",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
	return NCM_SPLINE_GET_CLASS (s)->copy_empty (s);
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
  g_assert (s->xv != NULL && s->yv != NULL);
	return ncm_spline_new (s, ncm_vector_dup (s->xv), ncm_vector_dup (s->yv), TRUE);
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
	g_assert (xv != NULL && yv != NULL);
	if (ncm_vector_len (xv) != ncm_vector_len (yv))
		g_error ("ncm_spline_set: knot and function values vector has not the same size");
	if (ncm_vector_len (xv) < NCM_SPLINE_GET_CLASS (s)->min_size (s))
		g_error ("ncm_spline_set: min size for [%s] is %zu but vector size is %u", NCM_SPLINE_GET_CLASS (s)->name (s),
		         NCM_SPLINE_GET_CLASS (s)->min_size (s), ncm_vector_len (xv));

	if (s->xv != NULL)
	{
		if (s->xv != xv)
		{
			ncm_vector_free (s->xv);
			s->xv = xv;
			ncm_vector_ref (xv);
		}
	}
	else
	{
		s->xv = xv;
		ncm_vector_ref (xv);
	}

	if (s->yv != NULL)
	{
		if (s->yv != yv)
		{
			ncm_vector_free (s->yv);
			s->yv = yv;
			ncm_vector_ref (yv);
		}
	}
	else
	{
		s->yv = yv;
		ncm_vector_ref (yv);
	}

	s->len = ncm_vector_len (xv);

	NCM_SPLINE_GET_CLASS (s)->reset (s);

	s->empty = FALSE;

	if (init)
		ncm_spline_prepare (s);

  if (s->acc != NULL)
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
  if (enable)
  {
    if (s->acc == NULL)
      s->acc = gsl_interp_accel_alloc ();
  }
  else if (s->acc != NULL)
  {
    gsl_interp_accel_free (s->acc);
    s->acc = NULL;
  }
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
  if (s->len != len)
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
  return s->len;
}

/**
 * ncm_spline_set_xv:
 * @s: a #NcmSpline.
 * @xv: #NcmVector of knots.
 * @init: TRUE to prepare @s or FALSE to not prepare it.
 *
 * This function sets @xv as the knot vector of the spline.
 *
 */
void
ncm_spline_set_xv (NcmSpline *s, NcmVector *xv, gboolean init)
{
  ncm_spline_set (s, xv, s->yv, init);
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
  ncm_spline_set (s, s->xv, yv, init);
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
	if (s->xv != NULL)
    return ncm_vector_ref (s->xv);
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
  if (s->yv != NULL)
    return ncm_vector_ref (s->yv);
  else
    return NULL;
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
  g_assert_cmpuint (s->len, >, 0);

  *lb = ncm_vector_get (s->xv, 0);
  *ub = ncm_vector_get (s->xv, s->len - 1);
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
/**
 * ncm_spline_prepare_base:
 * @s: a #NcmSpline
 *
 * This function computes the second derivatives of @s and it is used to prepare a
 * bidimensional spline.
 *
 */
/**
 * ncm_spline_eval:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 * Returns: The interpolated value of a function computed at @x.
 */
/**
 * ncm_spline_eval_deriv:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 *
 * Returns: The derivative of an interpolated function computed at @x.
 */
/**
 * ncm_spline_eval_deriv2:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 *
 * Returns: The second derivative of an interpolated function computed at @x.
 */
/**
 * ncm_spline_eval_deriv_nmax:
 * @s: a constant #NcmSpline
 * @x: x-coordinate value
 *
 *
 * Returns: The highest non null derivative of an interpolated function computed at @x.
 */
/**
 * ncm_spline_eval_integ:
 * @s: a constant #NcmSpline
 * @x0: lower integration limit
 * @x1: upper integration limit
 *
 *
 * Returns: The numerical integral of an interpolated function over the range [@x0, @x1].
 */
/**
 * ncm_spline_min_size:
 * @s: a constant #NcmSpline
 *
 *
 * Returns: Minimum number of knots required.
 */
/**
 * ncm_spline_get_index:
 * @s: a constant #NcmSpline
 * @x: a value of the abscissa axis
 *
 *
 * Returns: The index of the lower knot of the interval @x belongs to.
 */

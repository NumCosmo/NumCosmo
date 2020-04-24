/***************************************************************************
 *            ncm_spline2d.c
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
 * SECTION:ncm_spline2d
 * @title: NcmSpline2d
 * @short_description: Abstract class for implementing bidimensional splines.
 * @stability: Stable
 * @include: numcosmo/math/ncm_spline2d.h
 *
 * This class comprises all functions to provide a #NcmSpline2d, get its properties
 * and evaluate it given an interpolation method.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline2d.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_SPLINE,
  PROP_XV,
  PROP_YV,
  PROP_ZM,
  PROP_INIT,
  PROP_USE_ACC,
};

G_DEFINE_ABSTRACT_TYPE (NcmSpline2d, ncm_spline2d, G_TYPE_OBJECT);

static void
ncm_spline2d_init (NcmSpline2d *s2d)
{
  s2d->xv        = NULL;
  s2d->yv        = NULL;
  s2d->zm        = NULL;
  s2d->s         = NULL;
  s2d->empty     = TRUE;
  s2d->init      = FALSE;
  s2d->to_init   = FALSE;
  s2d->acc_x     = gsl_interp_accel_alloc ();
  s2d->acc_y     = gsl_interp_accel_alloc ();
  s2d->use_acc   = FALSE;
  s2d->no_stride = FALSE;
}

static void _ncm_spline2d_makeup (NcmSpline2d *s2d);

static void
_ncm_spline2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (object);
  g_return_if_fail (NCM_IS_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      s2d->s = g_value_dup_object (value);
      break;
    case PROP_XV:
      s2d->xv = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_YV:
      s2d->yv = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_ZM:
      s2d->zm = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_INIT:
    {
      s2d->to_init = g_value_get_boolean (value);
      if (s2d->to_init && (s2d->xv != NULL) && (s2d->yv != NULL) && (s2d->zm != NULL))
      {
        ncm_spline2d_prepare (s2d);
      }
      break;
    }
    case PROP_USE_ACC:
      ncm_spline2d_use_acc (s2d, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (object);
  g_return_if_fail (NCM_IS_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      g_value_set_object (value, s2d->s);
      break;
    case PROP_XV:
      g_value_set_object (value, s2d->xv);
      break;
    case PROP_YV:
      g_value_set_object (value, s2d->yv);
      break;
    case PROP_ZM:
      g_value_set_object (value, s2d->zm);
      break;
    case PROP_INIT:
      g_value_set_boolean (value, s2d->init);
      break;
    case PROP_USE_ACC:
      g_value_set_boolean (value, s2d->use_acc);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_spline2d_dispose (GObject *object)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (object);

  ncm_vector_clear (&s2d->xv);
  ncm_vector_clear (&s2d->yv);
  ncm_matrix_clear (&s2d->zm);
  ncm_spline_clear (&s2d->s);

  s2d->empty = TRUE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_parent_class)->dispose (object);
}

static void
ncm_spline2d_finalize (GObject *object)
{
  NcmSpline2d *s2d = NCM_SPLINE2D (object);

  g_clear_pointer (&s2d->acc_x, gsl_interp_accel_free);
  g_clear_pointer (&s2d->acc_y, gsl_interp_accel_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_parent_class)->finalize (object);
}

static void
ncm_spline2d_class_init (NcmSpline2dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_spline2d_set_property;
  object_class->get_property = &_ncm_spline2d_get_property;
  object_class->dispose      = &ncm_spline2d_dispose;
  object_class->finalize     = &ncm_spline2d_finalize;

  klass->copy_empty    = NULL;
  klass->reset         = NULL;
  klass->prepare       = NULL;
  klass->eval          = NULL;
  klass->dzdx          = NULL;
  klass->dzdy          = NULL;
  klass->d2zdxy        = NULL;
  klass->d2zdx2        = NULL;
  klass->d2zdy2        = NULL;
  klass->int_dx        = NULL;
  klass->int_dy        = NULL;
  klass->int_dxdy      = NULL;
  klass->int_dx_spline = NULL;
  klass->int_dy_spline = NULL;

  /**
   * NcmSpline2d:spline:
   *
   * #NcmSpline object used internally.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE,
                                   g_param_spec_object ("spline",
                                                        NULL,
                                                        "Spline",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmSpline2d:xv:
   *
   * #NcmVector x-knots.
   */
  g_object_class_install_property (object_class,
                                   PROP_XV,
                                   g_param_spec_object ("x-vector",
                                                        NULL,
                                                        "x vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmSpline2d:yv:
   *
   * #NcmVector y-knots.
   */
  g_object_class_install_property (object_class,
                                   PROP_YV,
                                   g_param_spec_object ("y-vector",
                                                        NULL,
                                                        "y vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmSpline2d:zm:
   *
   * #NcmMatrix z-values.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZM,
                                   g_param_spec_object ("z-matrix",
                                                        NULL,
                                                        "z matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmSpline2d:init:
   *
   * boolean whether to prepare the NcmSpline2d.
   */
  g_object_class_install_property (object_class,
                                   PROP_INIT,
                                   g_param_spec_boolean ("init",
                                                         NULL,
                                                         "init",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcmSpline2d:use_acc:
   *
   * boolean whether to use acc. 
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_ACC,
                                   g_param_spec_boolean ("use-acc",
                                                         NULL,
                                                         "Use accelerated bsearch",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
_ncm_spline2d_makeup (NcmSpline2d *s2d)
{
  if ((s2d->xv != NULL) && (s2d->yv != NULL) && (s2d->zm != NULL))
  {
    g_assert_cmpuint (ncm_vector_len (s2d->xv), ==, ncm_matrix_row_len (s2d->zm));
    g_assert_cmpuint (ncm_vector_len (s2d->yv), ==, ncm_matrix_col_len (s2d->zm));
    g_assert_cmpuint (ncm_vector_len (s2d->xv), >=, ncm_spline2d_min_size (s2d));
    g_assert_cmpuint (ncm_vector_len (s2d->yv), >=, ncm_spline2d_min_size (s2d));

    s2d->empty = FALSE;

    if ((ncm_vector_stride (s2d->xv) == 1) && (ncm_vector_stride (s2d->yv) == 1))
      s2d->no_stride = TRUE;
    else
      s2d->no_stride = FALSE;

    if (s2d->use_acc && !s2d->no_stride)
    {
      g_warning ("_ncm_spline2d_makeup: use-acc true but strided knots vectors, disabling use-add.");
      s2d->use_acc = FALSE;
    }

    NCM_SPLINE2D_GET_CLASS (s2d)->reset (s2d);

    if (s2d->to_init)
      ncm_spline2d_prepare (s2d);
  }
  return;
}

/**
 * ncm_spline2d_set:
 * @s2d: a #NcmSpline2d
 * @xv: a #NcmVector of knots
 * @yv: a #NcmVector of knots
 * @zm: a #NcmMatrix of the values of the function, to be interpolated, computed at @xv and @yv
 * @init: TRUE to prepare the #NcmSpline2d or FALSE to not prepare it
 *
 * This funtion sets @xv and @yv vectors and @zm matrix to @s2d.
 *
 */
void
ncm_spline2d_set (NcmSpline2d *s2d, NcmVector *xv, NcmVector *yv, NcmMatrix *zm, gboolean init)
{
  g_assert ((xv != NULL) && (yv != NULL) && (zm != NULL));

  s2d->xv = ncm_vector_ref (xv);
  s2d->yv = ncm_vector_ref (yv);
  s2d->zm = ncm_matrix_ref (zm);

  s2d->to_init = init;
  _ncm_spline2d_makeup (s2d);
  return;
}

/**
 * ncm_spline2d_copy_empty:
 * @s2d: a #NcmSpline2d
 *
 * This function copies the bidimensional spline @s2d into an initialized
 * empty #NcmSpline2d of a specific type.
 *
 * Returns: (transfer full): a #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_copy_empty (const NcmSpline2d *s2d)
{
  return NCM_SPLINE2D_GET_CLASS (s2d)->copy_empty (s2d);
}

/**
 * ncm_spline2d_copy:
 * @s2d: a #NcmSpline2d
 *
 * This function copies the two #NcmVector and the #NcmMatrix of the bidimensional
 * spline @s2d into those two #NcmVector and #NcmMatrix of a new #NcmSpline2d.
 *
 * Returns: (transfer full): A #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_copy (NcmSpline2d *s2d)
{
  NcmSpline2d *new_s2d = ncm_spline2d_copy_empty (s2d);
  if (!s2d->empty)
  {
    ncm_spline2d_set (new_s2d,
                      ncm_vector_dup (s2d->xv),
                      ncm_vector_dup (s2d->yv),
                      ncm_matrix_dup (s2d->zm),
                      s2d->init
                      );
  }

  return new_s2d;
}


/**
 * ncm_spline2d_new:
 * @s2d: a constant #NcmSpline2d
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of knots
 * @zm: #NcmMatrix of the values of the function, to be interpolated, computed at @xv and @yv
 * @init: TRUE to prepare the new #NcmSpline2d or FALSE to not prepare it
 *
 * This function returns a new #NcmSpline2d, where the knots of this new spline are given
 * in the #NcmVector @xv and @yv. The values of the function, at those knots, to be interpolated are
 * given in the #NcmMatrix @zm.
 *
 * Returns: (transfer full): A new #NcmSpline2d.
 */
NcmSpline2d *
ncm_spline2d_new (const NcmSpline2d *s2d, NcmVector *xv, NcmVector *yv, NcmMatrix *zm, gboolean init)
{
  NcmSpline2d *s2d_new = ncm_spline2d_copy_empty (s2d);
  ncm_spline2d_set (s2d_new, xv, yv, zm, init);
  return s2d_new;
}

/**
 * ncm_spline2d_min_size:
 * @s2d: a #NcmSpline2d
 *
 * Returns: The size of the #NcmSpline member of @s2d.
 */
guint
ncm_spline2d_min_size (NcmSpline2d *s2d)
{
  return ncm_spline_min_size (s2d->s);
}

/**
 * ncm_spline2d_prepare:
 * @s2d: a #NcmSpline2d
 *
 * This function prepares the bidimensional spline @s2d such that one can evaluate it (#ncm_spline2d_eval),
 * as well as to compute its integration in x, y or both directions.
 */
void ncm_spline2d_prepare (NcmSpline2d *s2d)
{
  NCM_SPLINE2D_GET_CLASS (s2d)->prepare (s2d);
}

/**
 * ncm_spline2d_free:
 * @s2d: a #NcmSpline2d
 *
 * Atomically decrements the reference count of @s2d by one. If the reference count drops to 0,
 * all memory allocated by @s2d is released.
 */
void
ncm_spline2d_free (NcmSpline2d *s2d)
{
  g_object_unref (s2d);
}

/**
 * ncm_spline2d_clear:
 * @s2d: a #NcmSpline2d
 *
 * Atomically decrements the reference count of @s2d by one. If the reference count drops to 0,
 * all memory allocated by @s2d is released. Set pointer to NULL.
 */
void
ncm_spline2d_clear (NcmSpline2d **s2d)
{
  g_clear_object (s2d);
}

/**
 * ncm_spline2d_use_acc:
 * @s2d: a #NcmSpline2d
 * @use_acc: a boolean
 * 
 * Whether to use accelerated bsearch to find the
 * right knots. When enabled evaluation functions
 * are not reentrant.
 * 
 */
void 
ncm_spline2d_use_acc (NcmSpline2d *s2d, gboolean use_acc)
{
  s2d->use_acc = use_acc;
  if ((s2d->xv != NULL) && (s2d->yv != NULL) && (s2d->zm != NULL))
  {
    if (s2d->use_acc && !s2d->no_stride)
    {
      g_warning ("ncm_spline2d_use_acc: use-acc true but strided knots vectors, disabling use-add.");
      s2d->use_acc = FALSE;
    }
  }
}

/**
 * ncm_spline2d_integ_dx: (virtual int_dx)
 * @s2d: a #NcmSpline2d
 * @xl: lower limit of integration
 * @xu: upper limit of integration
 * @y: y-coordinate value
 *
 * This function computes the integration in x over the interval [@xl, @xu] and
 * at @y.
 *
 * Returns: The numerical integral in x of an interpolated function over the range [@xl, @xu] and at @y.
 */
gdouble
ncm_spline2d_integ_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->int_dx (s2d, xl, xu, y);
}

/**
 * ncm_spline2d_integ_dy: (virtual int_dy)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @yl: lower limit of integration
 * @yu: upper limit of integration
 *
 * This function computes the integration in y over the interval [@yl, @yu] and
 * at @x.
 *
 * Returns: The numerical integral in y of an interpolated function over the range [@yl, @yu] and at @x.
 */
gdouble
ncm_spline2d_integ_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->int_dy (s2d, x, yl, yu);
}

/**
 * ncm_spline2d_integ_dxdy: (virtual int_dxdy)
 * @s2d: a #NcmSpline2d
 * @xl: lower limit of integration in the x-direction
 * @xu: upper limit of integration in the x-direction
 * @yl: lower limit of integration in the y-direction
 * @yu: upper limit of integration in the y-direction
 *
 * This function computes the integration in both x and y directions over the intervals
 * [@xl, @xu] and [@yl, @yu].
 *
 * Returns: The numerical integral in x and y of an interpolated function over the ranges [@xl, @xu] and [@yl, @yu].
 */
gdouble
ncm_spline2d_integ_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return NCM_SPLINE2D_GET_CLASS (s2d)->int_dxdy (s2d, xl, xu, yl, yu);
}

/**
 * ncm_spline2d_integ_dx_spline:
 * @s2d: a #NcmSpline2d
 * @xl: lower limit of integration x
 * @xu: upper limit of integration x
 *
 * This function computes the integral in x of the bidimensional interpolated function
 * over the range [@xl, @xu] resulting in a one dimensional function.
 *
 * Returns: (transfer full): A #NcmSpline.
 */
NcmSpline *
ncm_spline2d_integ_dx_spline (NcmSpline2d *s2d, gdouble xl, gdouble xu)
{
  g_assert (!s2d->empty);
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return ncm_spline_copy (NCM_SPLINE2D_GET_CLASS (s2d)->int_dx_spline (s2d, xl, xu));
}

/**
 * ncm_spline2d_integ_dy_spline:
 * @s2d: a #NcmSpline2d
 * @yl: lower limit of integration
 * @yu: upper limit of integration
 *
 * This function computes the integral in y of the bidimensional interpolated function
 * over the range [@yl, @yu] resulting in a one dimensional function.
 *
 * Returns: (transfer full): A #NcmSpline.
 */
NcmSpline *
ncm_spline2d_integ_dy_spline (NcmSpline2d *s2d, gdouble yl, gdouble yu)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return ncm_spline_copy (NCM_SPLINE2D_GET_CLASS (s2d)->int_dy_spline (s2d, yl, yu));
}

/**
 * ncm_spline2d_integ_dx_spline_val:
 * @s2d: a #NcmSpline2d
 * @xl: lower limit of integration
 * @xu: upper limit of integration
 * @y: y-coordinate value
 *
 * This function calls #ncm_spline2d_integ_dx_spline and evaluates the resulting
 * #NcmSpline at @y.
 *
 * Returns: The value of @s2d integrated in x over the range [@xl, @xu] and computed at @y.
 */
gdouble
ncm_spline2d_integ_dx_spline_val (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return ncm_spline_eval (NCM_SPLINE2D_GET_CLASS (s2d)->int_dx_spline (s2d, xl, xu), y);
}

/**
 * ncm_spline2d_integ_dy_spline_val:
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @yl: lower limit of integration
 * @yu: upper limit of integration
 *
 * This function calls #ncm_spline2d_integ_dy_spline and evaluates the resulting
 * #NcmSpline at @x.
 *
 * Returns: The value of @s2d integrated in y over the range [@yl, @yu] and computed at @x.
 */
gdouble
ncm_spline2d_integ_dy_spline_val (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return ncm_spline_eval (NCM_SPLINE2D_GET_CLASS (s2d)->int_dy_spline (s2d, yl, yu), x);
}

/**
 * ncm_spline2d_integ_dxdy_spline_x:
 * @s2d: a #NcmSpline2d
 * @xl: lower limit of integration in the x-direction
 * @xu: upper limit of integration in the x-direction
 * @yl: lower limit of integration in the y-direction
 * @yu: upper limit of integration in the y-direction
 *
 * This function calls #ncm_spline2d_integ_dx_spline and integrates the resulting
 * #NcmSpline over the interval [@yl, @yu].
 *
 * Returns: The value of @s2d integrated in x and y over the ranges [@xl, @xu] and [@yl, @yu], respectively.
 */
gdouble
ncm_spline2d_integ_dxdy_spline_x (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return ncm_spline_eval_integ (NCM_SPLINE2D_GET_CLASS (s2d)->int_dx_spline (s2d, xl, xu), yl, yu);
}

/**
 * ncm_spline2d_integ_dxdy_spline_y:
 * @s2d: a #NcmSpline2d
 * @xl: lower limit of integration in the x-direction
 * @xu: upper limit of integration in the x-direction
 * @yl: lower limit of integration in the y-direction
 * @yu: upper limit of integration in the y-direction
 *
 * This function calls #ncm_spline2d_integ_dy_spline and integrates the resulting
 * #NcmSpline over the interval [@xl, @xu].
 *
 * Returns: The value of @s2d integrated in x and y over the ranges [@xl, @xu] and [@yl, @yu], respectively.
 */
gdouble
ncm_spline2d_integ_dxdy_spline_y (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  if (!s2d->init)
    ncm_spline2d_prepare (s2d);
  return ncm_spline_eval_integ (NCM_SPLINE2D_GET_CLASS (s2d)->int_dy_spline (s2d, yl, yu), xl, xu);
}

/**
 * ncm_spline2d_eval:
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 *
 *
 * Returns: The interpolated value of a function computed at the point (@x, @y).
 */
/**
 * ncm_spline2d_deriv_dzdx: (virtual dzdx)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 * 
 * Returns: The interpolated derivative $\mathrm{d}z/\mathrm{d}x$ computed at the point (@x, @y).
 */
/**
 * ncm_spline2d_deriv_dzdy: (virtual dzdy)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 * 
 * Returns: The interpolated derivative $\mathrm{d}z/\mathrm{d}y$ computed at the point (@x, @y).
 */
/**
 * ncm_spline2d_deriv_d2zdxy: (virtual d2zdxy)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 * 
 * Returns: The interpolated derivative $\mathrm{d}^2z/\mathrm{d}x\mathrm{d}y$ computed at the point (@x, @y).
 */
/**
 * ncm_spline2d_deriv_d2zdx2: (virtual d2zdx2)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 * 
 * Returns: The interpolated derivative $\mathrm{d}^2z/\mathrm{d}x^2$ computed at the point (@x, @y).
 */
/**
 * ncm_spline2d_deriv_d2zdy2: (virtual d2zdy2)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 * 
 * Returns: The interpolated derivative $\mathrm{d}^2z/\mathrm{d}y^2$ computed at the point (@x, @y).
 */

/**
 * ncm_spline2dim_integ_total:
 * @s2d: a #NcmSpline2d
 * 
 * Returns: The numerical integral in both x and y directions of an interpolated function
 * over the entire valid ranges of x and y coordinates.
 */

/*******************************************************************************
 * Autoknots
 *******************************************************************************/

typedef struct __NcFunction2D_args
{
  gpointer data;
  gdouble x;
  gdouble y;
} _NcFunction2D_args;

/**
 * ncm_spline2d_set_function: (skip)
 * @s2d: a #NcmSpline2d
 * @ftype: a #NcmSplineFuncType
 * @Fx: function of x variable to be approximated by spline functions
 * @Fy: function of y variable to be approximated by spline functions
 * @xl: lower knot of x-coordinate
 * @xu: upper knot of x-coordinate
 * @yl: lower knot of y-coordinate
 * @yu: upper knot of y-coordinate
 * @rel_err: relative error between the function to be interpolated and the spline result
 * 
 * This function automatically determines the knots of @s2d in the intervals [@xl, @xu] and
 * [@yl, @yu] given a @ftype and @rel_error.
 * 
 * The functions @Fx and @Fy are the bidimensional function given at specific values of y and x, respectively.
 * These x and y values must be in the the intervals [@xl, @xu] and [@yl, @yu].
 */
void
ncm_spline2d_set_function (NcmSpline2d *s2d, NcmSplineFuncType ftype, gsl_function *Fx, gsl_function *Fy, gdouble xl, gdouble xu, gdouble yl, gdouble yu, gdouble rel_err)
{
  NcmSpline *s_x = ncm_spline_copy_empty (s2d->s);
  NcmSpline *s_y = ncm_spline_copy_empty (s2d->s);

  ncm_spline_set_func (s_x, ftype, Fx, xl, xu, 0, rel_err);
  ncm_spline_set_func (s_y, ftype, Fy, yl, yu, 0, rel_err);

/*
  printf ("x % 22.15g % 22.15g %u | y % 22.15g % 22.15g %u\n", xl, xu, ncm_vector_len (s_y->xv), yl, yu, ncm_vector_len (s_x->xv));
  ncm_vector_log_vals (s_y->xv, "XV", "% 22.15g", TRUE);
  ncm_vector_log_vals (s_x->xv, "XV", "% 22.15g", TRUE);
*/  
  {
    NcmMatrix *s_z = ncm_matrix_new (ncm_vector_len (s_y->xv), ncm_vector_len (s_x->xv));
    ncm_spline2d_set (s2d, s_x->xv, s_y->xv, s_z, FALSE);
    ncm_matrix_free (s_z);
  }
  
  ncm_spline_free (s_x);
  ncm_spline_free (s_y);

  return;
}


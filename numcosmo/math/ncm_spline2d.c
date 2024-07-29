/***************************************************************************
 *            ncm_spline2d.c
 *
 *  Sun Aug  1 17:17:08 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <vitenti@uel.br>
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

typedef struct _NcmSpline2dPrivate
{
  /*< private >*/
  GObject parent_instance;
  gboolean empty;
  gboolean init;
  gboolean to_init;
  NcmSpline *s;
  NcmVector *xv;
  NcmVector *yv;
  guint x_interv;
  guint y_interv;
  gdouble *x_data;
  gdouble *y_data;
  NcmMatrix *zm;
  gsl_interp_accel *acc_x;
  gsl_interp_accel *acc_y;
  gboolean use_acc;
  gboolean no_stride;
} NcmSpline2dPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmSpline2d, ncm_spline2d, G_TYPE_OBJECT)

static void
ncm_spline2d_init (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  self->xv        = NULL;
  self->yv        = NULL;
  self->x_interv  = 0;
  self->y_interv  = 0;
  self->x_data    = NULL;
  self->y_data    = NULL;
  self->zm        = NULL;
  self->s         = NULL;
  self->empty     = TRUE;
  self->init      = FALSE;
  self->to_init   = FALSE;
  self->acc_x     = gsl_interp_accel_alloc ();
  self->acc_y     = gsl_interp_accel_alloc ();
  self->use_acc   = FALSE;
  self->no_stride = FALSE;
}

static void _ncm_spline2d_makeup (NcmSpline2d *s2d);

static void
_ncm_spline2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpline2d *s2d                = NCM_SPLINE2D (object);
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  g_return_if_fail (NCM_IS_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      self->s = g_value_dup_object (value);
      break;
    case PROP_XV:
      self->xv = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_YV:
      self->yv = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_ZM:
      self->zm = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_INIT:
    {
      self->to_init = g_value_get_boolean (value);

      if (self->to_init && (self->xv != NULL) && (self->yv != NULL) && (self->zm != NULL))
        ncm_spline2d_prepare (s2d);

      break;
    }
    case PROP_USE_ACC:
      ncm_spline2d_use_acc (s2d, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpline2d *s2d                = NCM_SPLINE2D (object);
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  g_return_if_fail (NCM_IS_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      g_value_set_object (value, self->s);
      break;
    case PROP_XV:
      g_value_set_object (value, self->xv);
      break;
    case PROP_YV:
      g_value_set_object (value, self->yv);
      break;
    case PROP_ZM:
      g_value_set_object (value, self->zm);
      break;
    case PROP_INIT:
      g_value_set_boolean (value, self->init);
      break;
    case PROP_USE_ACC:
      g_value_set_boolean (value, self->use_acc);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_spline2d_dispose (GObject *object)
{
  NcmSpline2d *s2d                = NCM_SPLINE2D (object);
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  ncm_vector_clear (&self->xv);
  ncm_vector_clear (&self->yv);
  ncm_matrix_clear (&self->zm);
  ncm_spline_clear (&self->s);

  self->empty = TRUE;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_parent_class)->dispose (object);
}

static void
ncm_spline2d_finalize (GObject *object)
{
  NcmSpline2d *s2d                = NCM_SPLINE2D (object);
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  g_clear_pointer (&self->acc_x, gsl_interp_accel_free);
  g_clear_pointer (&self->acc_y, gsl_interp_accel_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline2d_parent_class)->finalize (object);
}

static void
ncm_spline2d_class_init (NcmSpline2dClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

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
  klass->eval_vec_y    = NULL;

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if ((self->xv != NULL) && (self->yv != NULL) && (self->zm != NULL))
  {
    g_assert_cmpuint (ncm_vector_len (self->xv), ==, ncm_matrix_row_len (self->zm));
    g_assert_cmpuint (ncm_vector_len (self->yv), ==, ncm_matrix_col_len (self->zm));
    g_assert_cmpuint (ncm_vector_len (self->xv), >=, ncm_spline2d_min_size (s2d));
    g_assert_cmpuint (ncm_vector_len (self->yv), >=, ncm_spline2d_min_size (s2d));

    self->empty    = FALSE;
    self->x_interv = ncm_vector_len (self->xv) - 1;
    self->y_interv = ncm_vector_len (self->yv) - 1;
    self->x_data   = ncm_vector_data (self->xv);
    self->y_data   = ncm_vector_data (self->yv);

    if ((ncm_vector_stride (self->xv) == 1) && (ncm_vector_stride (self->yv) == 1))
      self->no_stride = TRUE;
    else
      self->no_stride = FALSE;

    if (self->use_acc && !self->no_stride)
    {
      g_warning ("_ncm_spline2d_makeup: use-acc true but strided knots vectors, disabling use-add.");
      self->use_acc = FALSE;
    }

    NCM_SPLINE2D_GET_CLASS (s2d)->reset (s2d);

    if (self->to_init)
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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);
  NcmVector *old_xv               = self->xv;
  NcmVector *old_yv               = self->yv;
  NcmMatrix *old_zm               = self->zm;

  g_assert ((xv != NULL) && (yv != NULL) && (zm != NULL));

  self->xv = ncm_vector_ref (xv);
  self->yv = ncm_vector_ref (yv);
  self->zm = ncm_matrix_ref (zm);

  ncm_vector_clear (&old_xv);
  ncm_vector_clear (&old_yv);
  ncm_matrix_clear (&old_zm);

  self->to_init = init;
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
  return NCM_SPLINE2D_GET_CLASS ((NcmSpline2d *) s2d)->copy_empty (s2d);
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
  NcmSpline2d *new_s2d            = ncm_spline2d_copy_empty (s2d);
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->empty)
  {
    NcmVector *xv = ncm_vector_dup (self->xv);
    NcmVector *yv = ncm_vector_dup (self->yv);
    NcmMatrix *zm = ncm_matrix_dup (self->zm);

    ncm_spline2d_set (new_s2d, xv, yv, zm, self->init);

    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_matrix_free (zm);
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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return ncm_spline_min_size (self->s);
}

/**
 * ncm_spline2d_prepare:
 * @s2d: a #NcmSpline2d
 *
 * This function prepares the bi-dimensional spline @s2d such that one can evaluate it (#ncm_spline2d_eval),
 * as well as to compute its integration in x, y or both directions.
 */
void
ncm_spline2d_prepare (NcmSpline2d *s2d)
{
  NCM_SPLINE2D_GET_CLASS (s2d)->prepare (s2d);
}

/**
 * ncm_spline2d_ref:
 * @s2d: a #NcmSpline2d
 *
 * Atomically increases the reference count of @s2d by one.
 *
 * Returns: (transfer full): the same object @s2d.
 */
NcmSpline2d *
ncm_spline2d_ref (NcmSpline2d *s2d)
{
  return g_object_ref (s2d);
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
 * ncm_spline2d_set_init:
 * @s2d: a #NcmSpline2d
 * @init: a boolean
 *
 * Whether to mark the #NcmSpline2d as initialized.
 * This method is intended for internal use only.
 *
 */
void
ncm_spline2d_set_init (NcmSpline2d *s2d, gboolean init)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  self->init = init;
}

/**
 * ncm_spline2d_peek_spline:
 * @s2d: a #NcmSpline2d
 *
 * Get the #NcmSpline of the #NcmSpline2d.
 * This method is intended for internal use only.
 *
 * Returns: (transfer none): The #NcmSpline of the #NcmSpline2d.
 */
NcmSpline *
ncm_spline2d_peek_spline (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->s;
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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  self->use_acc = use_acc;

  if ((self->xv != NULL) && (self->yv != NULL) && (self->zm != NULL))
  {
    if (self->use_acc && !self->no_stride)
    {
      g_warning ("ncm_spline2d_use_acc: use-acc true but strided knots vectors, disabling use-add.");
      self->use_acc = FALSE;
    }
  }
}

/**
 * ncm_spline2d_peek_xv:
 * @s2d: a #NcmSpline2d
 *
 * Get the #NcmVector of knots in the x-direction.
 * This method is intended for internal use only.
 *
 * Returns: (transfer none): The #NcmVector of knots in the x-direction.
 */
NcmVector *
ncm_spline2d_peek_xv (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->xv;
}

/**
 * ncm_spline2d_peek_yv:
 * @s2d: a #NcmSpline2d
 *
 * Get the #NcmVector of knots in the y-direction.
 * This method is intended for internal use only.
 *
 * Returns: (transfer none): The #NcmVector of knots in the y-direction.
 */
NcmVector *
ncm_spline2d_peek_yv (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->yv;
}

/**
 * ncm_spline2d_peek_zm:
 * @s2d: a #NcmSpline2d
 *
 * Get the #NcmMatrix of the values of the function, to be interpolated, computed at the knots.
 * This method is intended for internal use only.
 *
 * Returns: (transfer none): The #NcmMatrix of the values of the function, to be interpolated, computed at the knots.
 */
NcmMatrix *
ncm_spline2d_peek_zm (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->zm;
}

/**
 * ncm_spline2d_peek_acc_x: (skip)
 * @s2d: a #NcmSpline2d
 *
 * Get the #gsl_interp_accel of the #NcmSpline2d in the x-direction.
 * This method is intended for internal use only.
 *
 * Returns: (transfer none): The #gsl_interp_accel of the #NcmSpline2d in the x-direction.
 */
gsl_interp_accel *
ncm_spline2d_peek_acc_x (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->acc_x;
}

/**
 * ncm_spline2d_peek_acc_y: (skip)
 * @s2d: a #NcmSpline2d
 *
 * Get the #gsl_interp_accel of the #NcmSpline2d in the y-direction.
 * This method is intended for internal use only.
 *
 * Returns: (transfer none): The #gsl_interp_accel of the #NcmSpline2d in the y-direction.
 */
gsl_interp_accel *
ncm_spline2d_peek_acc_y (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->acc_y;
}

/**
 * ncm_spline2d_is_init:
 * @s2d: a #NcmSpline2d
 *
 * Whether the #NcmSpline2d is initialized.
 *
 * Returns: TRUE if the #NcmSpline2d is initialized.
 */
gboolean
ncm_spline2d_is_init (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->init;
}

/**
 * ncm_spline2d_has_no_stride:
 * @s2d: a #NcmSpline2d
 *
 * Whether the #NcmSpline2d has stride 1 in both knots vectors.
 *
 * Returns: TRUE if the #NcmSpline2d has stride 1 in both knots vectors.
 */
gboolean
ncm_spline2d_has_no_stride (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->no_stride;
}

/**
 * ncm_spline2d_using_acc:
 * @s2d: a #NcmSpline2d
 *
 * Whether the #NcmSpline2d is using accelerated bsearch to find the
 * right knots.
 *
 * Returns: TRUE if the #NcmSpline2d is using accelerated bsearch to find the right knots.
 */
gboolean
ncm_spline2d_using_acc (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return self->use_acc;
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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  g_assert (!self->empty);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

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

gdouble
ncm_spline2d_eval (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  return NCM_SPLINE2D_GET_CLASS (s2d)->eval (s2d, x, y);
}

/**
 * ncm_spline2d_deriv_dzdx: (virtual dzdx)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 *
 * Returns: The interpolated derivative $\mathrm{d}z/\mathrm{d}x$ computed at the point (@x, @y).
 */
gdouble
ncm_spline2d_deriv_dzdx (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  return NCM_SPLINE2D_GET_CLASS (s2d)->dzdx (s2d, x, y);
}

/**
 * ncm_spline2d_deriv_dzdy: (virtual dzdy)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 *
 * Returns: The interpolated derivative $\mathrm{d}z/\mathrm{d}y$ computed at the point (@x, @y).
 */

gdouble
ncm_spline2d_deriv_dzdy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  return NCM_SPLINE2D_GET_CLASS (s2d)->dzdy (s2d, x, y);
}

/**
 * ncm_spline2d_deriv_d2zdxy: (virtual d2zdxy)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 *
 * Returns: The interpolated derivative $\mathrm{d}^2z/\mathrm{d}x\mathrm{d}y$ computed at the point (@x, @y).
 */

gdouble
ncm_spline2d_deriv_d2zdxy (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  return NCM_SPLINE2D_GET_CLASS (s2d)->d2zdxy (s2d, x, y);
}

/**
 * ncm_spline2d_deriv_d2zdx2: (virtual d2zdx2)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 *
 * Returns: The interpolated derivative $\mathrm{d}^2z/\mathrm{d}x^2$ computed at the point (@x, @y).
 */

gdouble
ncm_spline2d_deriv_d2zdx2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  return NCM_SPLINE2D_GET_CLASS (s2d)->d2zdx2 (s2d, x, y);
}

/**
 * ncm_spline2d_deriv_d2zdy2: (virtual d2zdy2)
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: y-coordinate value
 *
 * Returns: The interpolated derivative $\mathrm{d}^2z/\mathrm{d}y^2$ computed at the point (@x, @y).
 */

gdouble
ncm_spline2d_deriv_d2zdy2 (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  return NCM_SPLINE2D_GET_CLASS (s2d)->d2zdy2 (s2d, x, y);
}

/**
 * ncm_spline2dim_integ_total:
 * @s2d: a #NcmSpline2d
 *
 * Returns: The numerical integral in both x and y directions of an interpolated function
 * over the entire valid ranges of x and y coordinates.
 */
gdouble
ncm_spline2dim_integ_total (NcmSpline2d *s2d)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  return ncm_spline2d_integ_dxdy (s2d,
                                  ncm_vector_get (self->xv, 0),
                                  ncm_vector_get (self->xv, ncm_vector_len (self->xv) - 1),
                                  ncm_vector_get (self->yv, 0),
                                  ncm_vector_get (self->yv, ncm_vector_len (self->yv) - 1)
                                 );
}

/**
 * ncm_spline2d_eval_vec_y:
 * @s2d: a #NcmSpline2d
 * @x: x-coordinate value
 * @y: a #NcmVector
 * @order: (element-type size_t) (allow-none): an array of containing the order of the indices of @y
 * @res: (element-type gdouble): an array of the same size as @y to store the interpolated values
 *
 * Computes the interpolated values of a function computed at the point (@x, @y) for
 * each element of @y. The order of the indices of @y is given by @order.
 *
 */
void
ncm_spline2d_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res)
{
  NCM_SPLINE2D_GET_CLASS (s2d)->eval_vec_y (s2d, x, y, order, res);
}

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
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);
  NcmSpline *s_x                  = ncm_spline_copy_empty (self->s);
  NcmSpline *s_y                  = ncm_spline_copy_empty (self->s);

  ncm_spline_set_func (s_x, ftype, Fx, xl, xu, 0, rel_err);
  ncm_spline_set_func (s_y, ftype, Fy, yl, yu, 0, rel_err);
  {
    NcmVector *s_x_xv = ncm_spline_peek_xv (s_x);
    NcmVector *s_y_xv = ncm_spline_peek_xv (s_y);
    NcmMatrix *s_z    = ncm_matrix_new (ncm_vector_len (s_y_xv), ncm_vector_len (s_x_xv));

    ncm_spline2d_set (s2d, s_x_xv, s_y_xv, s_z, FALSE);
    ncm_matrix_free (s_z);
  }

  ncm_spline_free (s_x);
  ncm_spline_free (s_y);

  return;
}


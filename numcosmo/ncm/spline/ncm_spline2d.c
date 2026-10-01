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
 * NcmSpline2d:
 *
 * Abstract class for splines of two variables on a rectangular grid.
 *
 * The knots are #NcmSpline2d:x-vector and #NcmSpline2d:y-vector, and the value
 * $z(x_j, y_i)$ is the element in row $i$ and column $j$ of #NcmSpline2d:z-matrix. The
 * vectors and the matrix are referenced, not copied. #NcmSpline2d:spline sets the
 * interpolation along each direction and the minimum number of knots.
 *
 * Evaluation, the derivatives and the integrals prepare the spline when it is not
 * prepared. The integrals accept limits in either order, reversed limits changing the
 * sign; the functions returning a spline, and those built on them, require increasing
 * limits. With
 * #NcmSpline2d:use-acc, the knot searches use GSL accelerators held by the object, so
 * evaluation is not reentrant.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline2d.h"
#include "ncm/core/ncm_cfg.h"

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
      ncm_spline_clear (&self->s);
      self->s = g_value_dup_object (value);
      break;
    case PROP_XV:
      ncm_vector_clear (&self->xv);
      self->xv = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_YV:
      ncm_vector_clear (&self->yv);
      self->yv = g_value_dup_object (value);
      _ncm_spline2d_makeup (s2d);
      break;
    case PROP_ZM:
      ncm_matrix_clear (&self->zm);
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

static void _ncm_spline2d_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res);

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
  klass->eval_vec_y    = &_ncm_spline2d_eval_vec_y;

  /**
   * NcmSpline2d:spline:
   *
   * The spline type used along each direction.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE,
                                   g_param_spec_object ("spline",
                                                        NULL,
                                                        "Spline",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpline2d:x-vector:
   *
   * The knots in $x$.
   */
  g_object_class_install_property (object_class,
                                   PROP_XV,
                                   g_param_spec_object ("x-vector",
                                                        NULL,
                                                        "x vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpline2d:y-vector:
   *
   * The knots in $y$.
   */
  g_object_class_install_property (object_class,
                                   PROP_YV,
                                   g_param_spec_object ("y-vector",
                                                        NULL,
                                                        "y vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpline2d:z-matrix:
   *
   * The values, one row per knot in $y$ and one column per knot in $x$.
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
   * Whether the spline is prepared; setting it to %TRUE prepares the spline once the
   * knots and the values are set.
   */
  g_object_class_install_property (object_class,
                                   PROP_INIT,
                                   g_param_spec_boolean ("init",
                                                         NULL,
                                                         "init",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSpline2d:use-acc:
   *
   * Whether the knot searches use GSL accelerators, see ncm_spline2d_use_acc().
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_ACC,
                                   g_param_spec_boolean ("use-acc",
                                                         NULL,
                                                         "Use accelerated bsearch",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* Evaluates each element in turn; implementations can use @order to walk the knots once */
static void
_ncm_spline2d_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res)
{
  const guint len = ncm_vector_len (y);
  guint l;

  g_assert_cmpuint (len, ==, res->len);

  for (l = 0; l < len; l++)
    g_array_index (res, gdouble, l) = ncm_spline2d_eval (s2d, x, ncm_vector_get (y, l));
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
      g_warning ("_ncm_spline2d_makeup: use-acc true but strided knots vectors, disabling use-acc.");
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
 * @xv: the knots in $x$
 * @yv: the knots in $y$
 * @zm: the values, one row per element of @yv and one column per element of @xv
 * @init: whether to prepare the spline
 *
 * Sets the knots and the values of @s2d, see #NcmSpline2d; aborts when the dimensions do
 * not match or a direction has fewer knots than ncm_spline2d_min_size().
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
 * Creates an empty spline of the type and configuration of @s2d.
 *
 * Returns: (transfer full): a new #NcmSpline2d.
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
 * Creates a spline of the type of @s2d on copies of its knots and values, prepared when
 * @s2d is.
 *
 * Returns: (transfer full): a new #NcmSpline2d.
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
 * @xv: the knots in $x$
 * @yv: the knots in $y$
 * @zm: the values, one row per element of @yv and one column per element of @xv
 * @init: whether to prepare the spline
 *
 * Creates a spline of the type of @s2d and sets it, see ncm_spline2d_set().
 *
 * Returns: (transfer full): a new #NcmSpline2d.
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
 * Returns: the minimum number of knots in each direction, that of #NcmSpline2d:spline.
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
 * Prepares @s2d for evaluation.
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
 * Increases the reference count of @s2d by one.
 *
 * Returns: (transfer full): @s2d.
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
 * Decreases the reference count of @s2d by one.
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
 * If *@s2d is not %NULL, decreases its reference count by one and sets *@s2d to %NULL.
 */
void
ncm_spline2d_clear (NcmSpline2d **s2d)
{
  g_clear_object (s2d);
}

/**
 * ncm_spline2d_set_init:
 * @s2d: a #NcmSpline2d
 * @init: whether the spline is prepared
 *
 * Marks @s2d as prepared or not; for the implementations of prepare.
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
 * Gets #NcmSpline2d:spline.
 *
 * Returns: (transfer none): the spline type used along each direction.
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
 * @use_acc: whether to use GSL accelerators
 *
 * Sets #NcmSpline2d:use-acc. The accelerators require knot vectors with unit stride; for
 * strided knots it warns and leaves it off. With it on, evaluation is not reentrant.
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
      g_warning ("ncm_spline2d_use_acc: use-acc true but strided knots vectors, disabling use-acc.");
      self->use_acc = FALSE;
    }
  }
}

/**
 * ncm_spline2d_peek_xv:
 * @s2d: a #NcmSpline2d
 *
 * Gets #NcmSpline2d:x-vector.
 *
 * Returns: (transfer none): the knots in $x$.
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
 * Gets #NcmSpline2d:y-vector.
 *
 * Returns: (transfer none): the knots in $y$.
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
 * Gets #NcmSpline2d:z-matrix.
 *
 * Returns: (transfer none): the values at the knots.
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
 * Gets the accelerator of the knot search in $x$, for the implementations.
 *
 * Returns: (transfer none): the #gsl_interp_accel in $x$.
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
 * Gets the accelerator of the knot search in $y$, for the implementations.
 *
 * Returns: (transfer none): the #gsl_interp_accel in $y$.
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
 * Gets whether @s2d is prepared.
 *
 * Returns: %TRUE if @s2d is prepared.
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
 * Gets whether both knot vectors have unit stride.
 *
 * Returns: %TRUE if both knot vectors have unit stride.
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
 * Gets whether the knot searches use GSL accelerators, see ncm_spline2d_use_acc().
 *
 * Returns: %TRUE if the accelerators are in use.
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
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 * @y: the point in $y$
 *
 * Returns: $\int_{x_l}^{x_u} z(x, y)\,\mathrm{d}x$.
 */
gdouble
ncm_spline2d_integ_dx (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  if (xu < xl)
    return -NCM_SPLINE2D_GET_CLASS (s2d)->int_dx (s2d, xu, xl, y);

  return NCM_SPLINE2D_GET_CLASS (s2d)->int_dx (s2d, xl, xu, y);
}

/**
 * ncm_spline2d_integ_dy: (virtual int_dy)
 * @s2d: a #NcmSpline2d
 * @x: the point in $x$
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 *
 * Returns: $\int_{y_l}^{y_u} z(x, y)\,\mathrm{d}y$.
 */
gdouble
ncm_spline2d_integ_dy (NcmSpline2d *s2d, gdouble x, gdouble yl, gdouble yu)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  if (yu < yl)
    return -NCM_SPLINE2D_GET_CLASS (s2d)->int_dy (s2d, x, yu, yl);

  return NCM_SPLINE2D_GET_CLASS (s2d)->int_dy (s2d, x, yl, yu);
}

/**
 * ncm_spline2d_integ_dxdy: (virtual int_dxdy)
 * @s2d: a #NcmSpline2d
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 *
 * Returns: $\int_{y_l}^{y_u}\int_{x_l}^{x_u} z(x, y)\,\mathrm{d}x\,\mathrm{d}y$.
 */
gdouble
ncm_spline2d_integ_dxdy (NcmSpline2d *s2d, gdouble xl, gdouble xu, gdouble yl, gdouble yu)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);  /* LCOV_EXCL_LINE */

  {
    const gdouble sign = ((xu < xl) != (yu < yl)) ? -1.0 : 1.0;

    return sign * NCM_SPLINE2D_GET_CLASS (s2d)->int_dxdy (s2d, GSL_MIN (xl, xu), GSL_MAX (xl, xu), GSL_MIN (yl, yu), GSL_MAX (yl, yu));
  }
}

/**
 * ncm_spline2d_integ_dx_spline:
 * @s2d: a #NcmSpline2d
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 *
 * Returns: (transfer full): a spline in $y$ of $\int_{x_l}^{x_u} z(x, y)\,\mathrm{d}x$.
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
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 *
 * Returns: (transfer full): a spline in $x$ of $\int_{y_l}^{y_u} z(x, y)\,\mathrm{d}y$.
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
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 * @y: the point in $y$
 *
 * Evaluates the spline of ncm_spline2d_integ_dx_spline() at @y.
 *
 * Returns: $\int_{x_l}^{x_u} z(x, y)\,\mathrm{d}x$ from that spline.
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
 * @x: the point in $x$
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 *
 * Evaluates the spline of ncm_spline2d_integ_dy_spline() at @x.
 *
 * Returns: $\int_{y_l}^{y_u} z(x, y)\,\mathrm{d}y$ from that spline.
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
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 *
 * Integrates the spline of ncm_spline2d_integ_dx_spline() over [@yl, @yu].
 *
 * Returns: the double integral of $z$ from that spline.
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
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 *
 * Integrates the spline of ncm_spline2d_integ_dy_spline() over [@xl, @xu].
 *
 * Returns: the double integral of $z$ from that spline.
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
 * @x: the point in $x$
 * @y: the point in $y$
 *
 * Returns: $z(x, y)$.
 */

gdouble
ncm_spline2d_eval (NcmSpline2d *s2d, gdouble x, gdouble y)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);

  return NCM_SPLINE2D_GET_CLASS (s2d)->eval (s2d, x, y);
}

/**
 * ncm_spline2d_deriv_dzdx: (virtual dzdx)
 * @s2d: a #NcmSpline2d
 * @x: the point in $x$
 * @y: the point in $y$
 *
 * Returns: $\partial z/\partial x$ at (@x, @y).
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
 * @x: the point in $x$
 * @y: the point in $y$
 *
 * Returns: $\partial z/\partial y$ at (@x, @y).
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
 * @x: the point in $x$
 * @y: the point in $y$
 *
 * Returns: $\partial^2 z/\partial x\partial y$ at (@x, @y).
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
 * @x: the point in $x$
 * @y: the point in $y$
 *
 * Returns: $\partial^2 z/\partial x^2$ at (@x, @y).
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
 * @x: the point in $x$
 * @y: the point in $y$
 *
 * Returns: $\partial^2 z/\partial y^2$ at (@x, @y).
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
 * Returns: the integral of $z$ over the whole grid.
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
 * @x: the point in $x$
 * @y: the points in $y$
 * @order: (element-type size_t): the indices of @y in increasing order of its elements
 * @res: (element-type gdouble): the output array, of the length of @y
 *
 * Computes $z(x, y_l)$ into element $l$ of @res for every element $y_l$ of @y.
 * #NcmSpline2dBicubic visits them in the order given by @order, locating the cells in a
 * single pass; the other types evaluate each element separately and ignore @order.
 */
void
ncm_spline2d_eval_vec_y (NcmSpline2d *s2d, gdouble x, const NcmVector *y, GArray *order, GArray *res)
{
  NcmSpline2dPrivate * const self = ncm_spline2d_get_instance_private (s2d);

  if (!self->init)
    ncm_spline2d_prepare (s2d);

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
 * @Fx: the function along $x$, at a fixed $y$
 * @Fy: the function along $y$, at a fixed $x$
 * @xl: the lower limit in $x$
 * @xu: the upper limit in $x$
 * @yl: the lower limit in $y$
 * @yu: the upper limit in $y$
 * @rel_err: the relative tolerance
 *
 * Places the knots in $x$ with ncm_spline_set_func() on @Fx and those in $y$ on @Fy, and
 * sets @s2d to that grid with an unprepared value matrix filled with NaN: the caller
 * fills #NcmSpline2d:z-matrix and then calls ncm_spline2d_prepare().
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

    /* NaN until the caller fills it */
    ncm_matrix_set_all (s_z, GSL_NAN);
    ncm_spline2d_set (s2d, s_x_xv, s_y_xv, s_z, FALSE);
    ncm_matrix_free (s_z);
  }

  ncm_spline_free (s_x);
  ncm_spline_free (s_y);

  return;
}


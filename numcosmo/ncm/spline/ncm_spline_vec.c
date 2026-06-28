/***************************************************************************
 *            ncm_spline_vec.c
 *
 *  Sat Mar 15 19:53:22 2026
 *  Copyright  2026
 *  Sandro Dias Pinto Vitenti
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
 * NcmSplineVec:
 *
 * A vector-valued spline function $\vec{F}(x): \mathbb{R} \to \mathbb{R}^n$.
 *
 * This class represents a vector-valued function where each component is
 * interpolated by a spline sharing the same x-vector. The key optimization
 * is that a single binary search (via ncm_spline_get_index()) is performed,
 * followed by direct evaluation of all components using the index-based
 * methods (*_idx()).
 *
 * The object can be constructed from:
 * - An x-vector and a matrix where each row is a component y-vector
 * - An x-vector and a GPtrArray of y-vectors
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline_vec.h"
#include "ncm/core/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_SPLINE,
  PROP_LEN,
};

struct _NcmSplineVec
{
  /*< private >*/
  GObject parent_instance;
  NcmSpline *base_spline;
  GPtrArray *spline_array;
  guint len;
  gboolean init;
};

G_DEFINE_TYPE (NcmSplineVec, ncm_spline_vec, G_TYPE_OBJECT)

static void
ncm_spline_vec_init (NcmSplineVec *sv)
{
  sv->base_spline  = NULL;
  sv->spline_array = g_ptr_array_new ();
  sv->len          = 0;
  sv->init         = FALSE;

  g_ptr_array_set_free_func (sv->spline_array, (GDestroyNotify) ncm_spline_free);
}

static void
_ncm_spline_vec_dispose (GObject *object)
{
  NcmSplineVec *sv = NCM_SPLINE_VEC (object);

  ncm_spline_clear (&sv->base_spline);
  g_clear_pointer (&sv->spline_array, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_vec_parent_class)->dispose (object);
}

static void
_ncm_spline_vec_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_vec_parent_class)->finalize (object);
}

static void
_ncm_spline_vec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineVec *sv = NCM_SPLINE_VEC (object);

  g_return_if_fail (NCM_IS_SPLINE_VEC (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      ncm_spline_clear (&sv->base_spline);
      sv->base_spline = g_value_dup_object (value);
      break;
    case PROP_LEN:
      sv->len = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_vec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineVec *sv = NCM_SPLINE_VEC (object);

  g_return_if_fail (NCM_IS_SPLINE_VEC (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      g_value_set_object (value, sv->base_spline);
      break;
    case PROP_LEN:
      g_value_set_uint (value, sv->len);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_spline_vec_class_init (NcmSplineVecClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = &_ncm_spline_vec_dispose;
  object_class->finalize     = &_ncm_spline_vec_finalize;
  object_class->set_property = &_ncm_spline_vec_set_property;
  object_class->get_property = &_ncm_spline_vec_get_property;

  /**
   * NcmSplineVec:spline:
   *
   * The base spline type used to create component splines.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE,
                                   g_param_spec_object ("spline",
                                                        NULL,
                                                        "Base spline",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSplineVec:len:
   *
   * The number of component splines (dimension of the vector function).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      NULL,
                                                      "Number of components",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spline_vec_new:
 * @s: a #NcmSpline
 * @xv: a #NcmVector containing the x-coordinates
 * @ym: a #NcmMatrix where each row is a y-vector for a component
 * @init: whether to initialize the splines immediately
 *
 * Creates a new #NcmSplineVec from a matrix of y-values. Each row of @ym
 * corresponds to one component of the vector function. All components share
 * the same x-vector @xv.
 *
 * Returns: (transfer full): a new #NcmSplineVec
 */
NcmSplineVec *
ncm_spline_vec_new (const NcmSpline *s, NcmVector *xv, NcmMatrix *ym, const gboolean init)
{
  NcmSplineVec *sv = g_object_new (NCM_TYPE_SPLINE_VEC,
                                   "spline", s,
                                   NULL);

  ncm_spline_vec_set (sv, xv, ym, init);

  return sv;
}

/**
 * ncm_spline_vec_new_gpa:
 * @s: a #NcmSpline
 * @xv: a #NcmVector containing the x-coordinates
 * @yv: (element-type NcmVector): a #GPtrArray of #NcmVector containing the y-coordinates for each component
 * @init: whether to initialize the splines immediately
 *
 * Creates a new #NcmSplineVec from a #GPtrArray of y-vectors. Each element
 * of @yv corresponds to one component of the vector function. All components
 * share the same x-vector @xv.
 *
 * Returns: (transfer full): a new #NcmSplineVec
 */
NcmSplineVec *
ncm_spline_vec_new_gpa (const NcmSpline *s, NcmVector *xv, GPtrArray *yv, const gboolean init)
{
  NcmSplineVec *sv = g_object_new (NCM_TYPE_SPLINE_VEC,
                                   "spline", s,
                                   NULL);

  ncm_spline_vec_set_gpa (sv, xv, yv, init);

  return sv;
}

/**
 * ncm_spline_vec_ref:
 * @sv: a #NcmSplineVec
 *
 * Increases the reference count of @sv atomically.
 *
 * Returns: (transfer full): @sv
 */
NcmSplineVec *
ncm_spline_vec_ref (NcmSplineVec *sv)
{
  return g_object_ref (sv);
}

/**
 * ncm_spline_vec_free:
 * @sv: a #NcmSplineVec
 *
 * Atomically decrements the reference count of @sv by one. If the reference
 * count drops to 0, all memory allocated by @sv is released.
 */
void
ncm_spline_vec_free (NcmSplineVec *sv)
{
  g_object_unref (sv);
}

/**
 * ncm_spline_vec_clear:
 * @sv: a #NcmSplineVec
 *
 * Atomically decrements the reference count of @sv by one. If the reference
 * count drops to 0, all memory allocated by @sv is released. Sets the pointer
 * to NULL.
 */
void
ncm_spline_vec_clear (NcmSplineVec **sv)
{
  g_clear_object (sv);
}

/**
 * ncm_spline_vec_set:
 * @sv: a #NcmSplineVec
 * @xv: a #NcmVector containing the x-coordinates
 * @ym: a #NcmMatrix where each row is a y-vector for a component
 * @init: whether to initialize the splines immediately
 *
 * Sets the data for @sv from a matrix. Each row of @ym corresponds to one
 * component of the vector function.
 */
void
ncm_spline_vec_set (NcmSplineVec *sv, NcmVector *xv, NcmMatrix *ym, gboolean init)
{
  const guint nrows = ncm_matrix_nrows (ym);
  const guint ncols = ncm_matrix_ncols (ym);
  guint i;

  g_assert_cmpuint (ncm_vector_len (xv), ==, ncols);

  /* Clear existing splines */
  g_ptr_array_set_size (sv->spline_array, 0);
  sv->len  = nrows;
  sv->init = FALSE;

  /* Create a spline for each row */
  for (i = 0; i < nrows; i++)
  {
    NcmSpline *s_i  = ncm_spline_copy_empty (sv->base_spline);
    NcmVector *yv_i = ncm_matrix_get_row (ym, i);

    ncm_spline_set (s_i, xv, yv_i, FALSE);
    ncm_vector_free (yv_i);

    g_ptr_array_add (sv->spline_array, s_i);
  }

  if (init)
    ncm_spline_vec_prepare (sv);
}

/**
 * ncm_spline_vec_set_gpa:
 * @sv: a #NcmSplineVec
 * @xv: a #NcmVector containing the x-coordinates
 * @yv: (element-type NcmVector): a #GPtrArray of #NcmVector containing the y-coordinates for each component
 * @init: whether to initialize the splines immediately
 *
 * Sets the data for @sv from a #GPtrArray of y-vectors. Each element of @yv
 * corresponds to one component of the vector function.
 */
void
ncm_spline_vec_set_gpa (NcmSplineVec *sv, NcmVector *xv, GPtrArray *yv, gboolean init)
{
  const guint len = yv->len;
  guint i;

  /* Clear existing splines */
  g_ptr_array_set_size (sv->spline_array, 0);
  sv->len  = len;
  sv->init = FALSE;

  /* Create a spline for each y-vector */
  for (i = 0; i < len; i++)
  {
    NcmSpline *s_i  = ncm_spline_copy_empty (sv->base_spline);
    NcmVector *yv_i = g_ptr_array_index (yv, i);

    g_assert_cmpuint (ncm_vector_len (xv), ==, ncm_vector_len (yv_i));

    ncm_spline_set (s_i, xv, yv_i, FALSE);

    g_ptr_array_add (sv->spline_array, s_i);
  }

  if (init)
    ncm_spline_vec_prepare (sv);
}

/**
 * ncm_spline_vec_prepare:
 * @sv: a #NcmSplineVec
 *
 * Prepares all component splines for evaluation.
 */
void
ncm_spline_vec_prepare (NcmSplineVec *sv)
{
  guint i;

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i = g_ptr_array_index (sv->spline_array, i);

    ncm_spline_prepare (s_i);
  }

  sv->init = TRUE;
}

/**
 * ncm_spline_vec_is_init:
 * @sv: a #NcmSplineVec
 *
 * Returns: whether the spline vector is initialized
 */
gboolean
ncm_spline_vec_is_init (NcmSplineVec *sv)
{
  return sv->init;
}

/**
 * ncm_spline_vec_get_len:
 * @sv: a #NcmSplineVec
 *
 * Returns: the number of component splines
 */
guint
ncm_spline_vec_get_len (NcmSplineVec *sv)
{
  return sv->len;
}

/**
 * ncm_spline_vec_get_nknots:
 * @sv: a #NcmSplineVec
 *
 * Gets the number of knot points (x-values) in the spline. All component
 * splines share the same x-vector, so this returns the length of that vector.
 *
 * Returns: the number of knot points
 */
guint
ncm_spline_vec_get_nknots (NcmSplineVec *sv)
{
  g_assert (sv->len > 0);
  g_assert (sv->init);

  return ncm_spline_get_len (g_ptr_array_index (sv->spline_array, 0));
}

/**
 * ncm_spline_vec_peek_spline:
 * @sv: a #NcmSplineVec
 * @i: component index
 *
 * Gets the @i-th component spline without increasing its reference count.
 *
 * Returns: (transfer none): the @i-th component spline
 */
NcmSpline *
ncm_spline_vec_peek_spline (NcmSplineVec *sv, guint i)
{
  g_assert_cmpuint (i, <, sv->len);

  return g_ptr_array_index (sv->spline_array, i);
}

/**
 * ncm_spline_vec_eval:
 * @sv: a #NcmSplineVec
 * @x: x-coordinate
 * @res: a #NcmVector to store the result
 *
 * Evaluates the vector function at @x. The key optimization is that a single
 * binary search is performed via ncm_spline_get_index(), then all components
 * are evaluated using ncm_spline_eval_idx().
 *
 * The vector @res must have length equal to the number of components.
 */
void
ncm_spline_vec_eval (NcmSplineVec *sv, const gdouble x, NcmVector *res)
{
  const gsize idx = ncm_spline_get_index (g_ptr_array_index (sv->spline_array, 0), x);
  guint i;

  g_assert_cmpuint (ncm_vector_len (res), ==, sv->len);
  g_assert (sv->init);

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i    = g_ptr_array_index (sv->spline_array, i);
    const gdouble y_i = ncm_spline_eval_idx (s_i, x, idx);

    ncm_vector_set (res, i, y_i);
  }
}

/**
 * ncm_spline_vec_deriv:
 * @sv: a #NcmSplineVec
 * @x: x-coordinate
 * @res: a #NcmVector to store the result
 *
 * Evaluates the derivative of the vector function at @x. Uses a single
 * binary search followed by index-based derivative evaluation for all
 * components.
 *
 * The vector @res must have length equal to the number of components.
 */
void
ncm_spline_vec_deriv (NcmSplineVec *sv, const gdouble x, NcmVector *res)
{
  const gsize idx = ncm_spline_get_index (g_ptr_array_index (sv->spline_array, 0), x);
  guint i;

  g_assert_cmpuint (ncm_vector_len (res), ==, sv->len);
  g_assert (sv->init);

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i     = g_ptr_array_index (sv->spline_array, i);
    const gdouble dy_i = ncm_spline_eval_deriv_idx (s_i, x, idx);

    ncm_vector_set (res, i, dy_i);
  }
}

/**
 * ncm_spline_vec_integ:
 * @sv: a #NcmSplineVec
 * @xi: initial x-coordinate
 * @xf: final x-coordinate
 * @res: a #NcmVector to store the result
 *
 * Evaluates the integral of the vector function from @xi to @xf. Uses index
 * lookups for both endpoints followed by index-based integration for all
 * components.
 *
 * The vector @res must have length equal to the number of components.
 */
void
ncm_spline_vec_integ (NcmSplineVec *sv, const gdouble xi, const gdouble xf, NcmVector *res)
{
  NcmSpline *s_0    = g_ptr_array_index (sv->spline_array, 0);
  const gsize idx_i = ncm_spline_get_index (s_0, xi);
  const gsize idx_f = ncm_spline_get_index (s_0, xf);
  guint i;

  g_assert_cmpuint (ncm_vector_len (res), ==, sv->len);
  g_assert (sv->init);

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i        = g_ptr_array_index (sv->spline_array, i);
    const gdouble integ_i = ncm_spline_eval_integ_idx (s_i, xi, idx_i, xf, idx_f);

    ncm_vector_set (res, i, integ_i);
  }
}

/**
 * ncm_spline_vec_eval_array:
 * @sv: a #NcmSplineVec
 * @x: x-coordinate
 * @res: (out callee-allocates) (element-type gdouble): result array
 *
 * Evaluates the vector function at @x and stores the result in @res.
 * If *@res is %NULL, a new #GArray is created. Otherwise, the existing array
 * is reused.
 */
void
ncm_spline_vec_eval_array (NcmSplineVec *sv, const gdouble x, GArray **res)
{
  const gsize idx = ncm_spline_get_index (g_ptr_array_index (sv->spline_array, 0), x);
  guint i;

  g_assert (sv->init);

  if (*res == NULL)
    *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), sv->len);

  g_array_set_size (*res, sv->len);

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i    = g_ptr_array_index (sv->spline_array, i);
    const gdouble y_i = ncm_spline_eval_idx (s_i, x, idx);

    g_array_index (*res, gdouble, i) = y_i;
  }
}

/**
 * ncm_spline_vec_deriv_array:
 * @sv: a #NcmSplineVec
 * @x: x-coordinate
 * @res: (out callee-allocates) (element-type gdouble): result array
 *
 * Evaluates the derivative of the vector function at @x and stores the result
 * in @res. If *@res is %NULL, a new #GArray is created. Otherwise, the
 * existing array is reused.
 */
void
ncm_spline_vec_deriv_array (NcmSplineVec *sv, const gdouble x, GArray **res)
{
  const gsize idx = ncm_spline_get_index (g_ptr_array_index (sv->spline_array, 0), x);
  guint i;

  g_assert (sv->init);

  if (*res == NULL)
    *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), sv->len);

  g_array_set_size (*res, sv->len);

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i     = g_ptr_array_index (sv->spline_array, i);
    const gdouble dy_i = ncm_spline_eval_deriv_idx (s_i, x, idx);

    g_array_index (*res, gdouble, i) = dy_i;
  }
}

/**
 * ncm_spline_vec_integ_array:
 * @sv: a #NcmSplineVec
 * @xi: initial x-coordinate
 * @xf: final x-coordinate
 * @res: (out callee-allocates) (element-type gdouble): result array
 *
 * Evaluates the integral of the vector function from @xi to @xf and stores
 * the result in @res. If *@res is %NULL, a new #GArray is created.
 * Otherwise, the existing array is reused.
 */
void
ncm_spline_vec_integ_array (NcmSplineVec *sv, const gdouble xi, const gdouble xf, GArray **res)
{
  NcmSpline *s_0    = g_ptr_array_index (sv->spline_array, 0);
  const gsize idx_i = ncm_spline_get_index (s_0, xi);
  const gsize idx_f = ncm_spline_get_index (s_0, xf);
  guint i;

  g_assert (sv->init);

  if (*res == NULL)
    *res = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), sv->len);

  g_array_set_size (*res, sv->len);

  for (i = 0; i < sv->len; i++)
  {
    NcmSpline *s_i        = g_ptr_array_index (sv->spline_array, i);
    const gdouble integ_i = ncm_spline_eval_integ_idx (s_i, xi, idx_i, xf, idx_f);

    g_array_index (*res, gdouble, i) = integ_i;
  }
}


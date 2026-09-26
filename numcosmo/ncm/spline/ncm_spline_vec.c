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
 * Vector-valued function $\vec{F}(x)$ whose components are splines on the same knots.
 *
 * Each component is an empty copy of #NcmSplineVec:spline, see ncm_spline_copy_empty(),
 * holding the shared knot vector and its own values; both vectors are referenced, not
 * copied. The evaluation functions locate the interval once, with ncm_spline_get_index()
 * on the first component, and pass it to ncm_spline_eval_idx(),
 * ncm_spline_eval_deriv_idx() or ncm_spline_eval_integ_idx(). A #NcmSplineVec has at
 * least one component.
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
   * The spline copied for each component.
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
   * The number of components.
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
 * @s: the spline copied for each component
 * @xv: the knots
 * @ym: the values at @xv, one row per component
 * @init: whether to prepare the splines
 *
 * Creates a #NcmSplineVec with one component per row of @ym, see ncm_spline_vec_set().
 *
 * Returns: (transfer full): a new #NcmSplineVec.
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
 * @s: the spline copied for each component
 * @xv: the knots
 * @yv: (element-type NcmVector): the values at @xv, one vector per component
 * @init: whether to prepare the splines
 *
 * Creates a #NcmSplineVec with one component per element of @yv, see
 * ncm_spline_vec_set_gpa().
 *
 * Returns: (transfer full): a new #NcmSplineVec.
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
 * Increases the reference count of @sv by one.
 *
 * Returns: (transfer full): @sv.
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
 * Decreases the reference count of @sv by one.
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
 * If *@sv is not %NULL, decreases its reference count by one and sets *@sv to %NULL.
 */
void
ncm_spline_vec_clear (NcmSplineVec **sv)
{
  g_clear_object (sv);
}

/**
 * ncm_spline_vec_set:
 * @sv: a #NcmSplineVec
 * @xv: the knots
 * @ym: the values at @xv, one row per component
 * @init: whether to prepare the splines
 *
 * Replaces the components of @sv by one per row of @ym, all on the knots @xv. The
 * components reference @xv and views of the rows of @ym. The number of columns of @ym
 * must equal the length of @xv, and @ym must have at least one row.
 */
void
ncm_spline_vec_set (NcmSplineVec *sv, NcmVector *xv, NcmMatrix *ym, gboolean init)
{
  const guint nrows = ncm_matrix_nrows (ym);
  const guint ncols = ncm_matrix_ncols (ym);
  guint i;

  g_assert_cmpuint (ncm_vector_len (xv), ==, ncols);

  if (nrows == 0)
    g_error ("ncm_spline_vec_set: the matrix has no rows, so there are no components.");

  g_ptr_array_set_size (sv->spline_array, 0);
  sv->len  = nrows;
  sv->init = FALSE;

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
 * @xv: the knots
 * @yv: (element-type NcmVector): the values at @xv, one vector per component
 * @init: whether to prepare the splines
 *
 * Replaces the components of @sv by one per element of @yv, all on the knots @xv. The
 * components reference @xv and the elements of @yv, each of which must have the length
 * of @xv; @yv must not be empty.
 */
void
ncm_spline_vec_set_gpa (NcmSplineVec *sv, NcmVector *xv, GPtrArray *yv, gboolean init)
{
  const guint len = yv->len;
  guint i;

  if (len == 0)
    g_error ("ncm_spline_vec_set_gpa: the array is empty, so there are no components.");

  g_ptr_array_set_size (sv->spline_array, 0);
  sv->len  = len;
  sv->init = FALSE;

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
 * Prepares every component.
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
 * Gets whether @sv was prepared since its components were last set.
 *
 * Returns: %TRUE if @sv is prepared.
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
 * Gets #NcmSplineVec:len.
 *
 * Returns: the number of components.
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
 * Gets the number of knots. @sv must be prepared and have at least one component.
 *
 * Returns: the number of knots.
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
 * @i: the component index
 *
 * Gets the component @i.
 *
 * Returns: (transfer none): the spline of component @i.
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
 * @x: the point
 * @res: the output vector, of length #NcmSplineVec:len
 *
 * Computes $\vec{F}(x)$ into @res.
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
 * @x: the point
 * @res: the output vector, of length #NcmSplineVec:len
 *
 * Computes $\mathrm{d}\vec{F}/\mathrm{d}x$ at @x into @res.
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
 * @xi: the lower limit
 * @xf: the upper limit
 * @res: the output vector, of length #NcmSplineVec:len
 *
 * Computes $\int_{x_i}^{x_f}\vec{F}(x)\,\mathrm{d}x$ into @res.
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
 * @x: the point
 * @res: (out callee-allocates) (element-type gdouble): the output array
 *
 * Computes $\vec{F}(x)$ into *@res, resized to #NcmSplineVec:len; a new #GArray is
 * created when *@res is %NULL.
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
 * @x: the point
 * @res: (out callee-allocates) (element-type gdouble): the output array
 *
 * Computes $\mathrm{d}\vec{F}/\mathrm{d}x$ at @x into *@res, resized to
 * #NcmSplineVec:len; a new #GArray is created when *@res is %NULL.
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
 * @xi: the lower limit
 * @xf: the upper limit
 * @res: (out callee-allocates) (element-type gdouble): the output array
 *
 * Computes $\int_{x_i}^{x_f}\vec{F}(x)\,\mathrm{d}x$ into *@res, resized to
 * #NcmSplineVec:len; a new #GArray is created when *@res is %NULL.
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


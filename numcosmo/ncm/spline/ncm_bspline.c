/***************************************************************************
 *            ncm_bspline.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_bspline.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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
 * NcmBSpline:
 *
 * A B-spline given by its order, knots and control points, read-only.
 *
 * The object is the function
 * $$
 * s(x) = \sum_{i=0}^{n-1} c^i B_{i,k}(x),
 * $$
 * where $k$ is #NcmBSpline:order, $c^i$ are the $n$ #NcmBSpline:control-points and
 * $B_{i,k}$ are the B-spline basis functions of order $k$ on the $n + k$ non-decreasing
 * #NcmBSpline:knots $t_j$. It holds no samples and fits nothing:
 * ncm_spline_bspline_dup_bspline() builds one from a prepared #NcmSplineBSpline, and
 * both evaluate with the same code. Outside $[t_{k-1}, t_n]$, the range given by
 * ncm_bspline_get_range(), the value extrapolates the edge polynomial.
 *
 * All three properties are set at construction only. The two vectors must be contiguous
 * (stride one). They are referenced, not copied, and must not be modified afterwards.
 * ncm_bspline_eval() uses stack scratch and writes nothing in the object, so it may be
 * called concurrently.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_bspline.h"
#include "ncm/spline/ncm_bspline_private.h"
#include "ncm/spline/ncm_spline_bspline.h"

struct _NcmBSpline
{
  /*< private >*/
  GObject parent_instance;
  guint order;
  NcmVector *knots;
  NcmVector *control_points;
  gsize n_control;  /* length of control_points */
  const gdouble *t; /* data of knots */
  const gdouble *c; /* data of control_points */
};

G_DEFINE_TYPE (NcmBSpline, ncm_bspline, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_ORDER,
  PROP_KNOTS,
  PROP_CONTROL_POINTS,
};

static void
ncm_bspline_init (NcmBSpline *bs)
{
  bs->order          = 0;
  bs->knots          = NULL;
  bs->control_points = NULL;
  bs->n_control      = 0;
  bs->t              = NULL;
  bs->c              = NULL;
}

static void
_ncm_bspline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmBSpline *bs = NCM_BSPLINE (object);

  g_return_if_fail (NCM_IS_BSPLINE (object));

  switch (prop_id)
  {
    case PROP_ORDER:
      bs->order = g_value_get_uint (value);
      break;
    case PROP_KNOTS:
      bs->knots = g_value_dup_object (value);
      break;
    case PROP_CONTROL_POINTS:
      bs->control_points = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_bspline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmBSpline *bs = NCM_BSPLINE (object);

  g_return_if_fail (NCM_IS_BSPLINE (object));

  switch (prop_id)
  {
    case PROP_ORDER:
      g_value_set_uint (value, bs->order);
      break;
    case PROP_KNOTS:
      g_value_set_object (value, bs->knots);
      break;
    case PROP_CONTROL_POINTS:
      g_value_set_object (value, bs->control_points);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Checks the three properties against each other and stores what ncm_bspline_eval()
 * reads. Aborts when a vector is missing or not contiguous, when the lengths do not
 * satisfy len (knots) = len (control_points) + order with at least order control
 * points, or when the knots decrease.
 */
static void
_ncm_bspline_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (ncm_bspline_parent_class)->constructed (object);
  {
    NcmBSpline *bs = NCM_BSPLINE (object);
    gsize n_knots, i;

    if ((bs->knots == NULL) || (bs->control_points == NULL))
      g_error ("NcmBSpline: both the knots and the control points are required.");

    if ((ncm_vector_stride (bs->knots) != 1) || (ncm_vector_stride (bs->control_points) != 1))
      g_error ("NcmBSpline: the knots and the control points must have stride one.");

    n_knots       = ncm_vector_len (bs->knots);
    bs->n_control = ncm_vector_len (bs->control_points);

    if ((n_knots != bs->n_control + bs->order) || (bs->n_control < bs->order))
      g_error ("NcmBSpline: order %u needs at least %u control points and as many knots as control points plus %u, got %u control points and %u knots.",
               bs->order, bs->order, bs->order, (guint) bs->n_control, (guint) n_knots);

    bs->t = ncm_vector_data (bs->knots);
    bs->c = ncm_vector_data (bs->control_points);

    for (i = 1; i < n_knots; i++)
    {
      if (bs->t[i] < bs->t[i - 1])
        g_error ("NcmBSpline: the knots must not decrease, knot %u is % 22.15g after % 22.15g.",
                 (guint) i, bs->t[i], bs->t[i - 1]);
    }
  }
}

static void
_ncm_bspline_dispose (GObject *object)
{
  NcmBSpline *bs = NCM_BSPLINE (object);

  ncm_vector_clear (&bs->knots);
  ncm_vector_clear (&bs->control_points);

  /* Chain up: end */
  G_OBJECT_CLASS (ncm_bspline_parent_class)->dispose (object);
}

static void
ncm_bspline_class_init (NcmBSplineClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_bspline_set_property;
  object_class->get_property = &_ncm_bspline_get_property;
  object_class->constructed  = &_ncm_bspline_constructed;
  object_class->dispose      = &_ncm_bspline_dispose;

  /**
   * NcmBSpline:order:
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
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmBSpline:knots:
   *
   * The knots, non-decreasing: as many as the control points plus #NcmBSpline:order.
   */
  g_object_class_install_property (object_class,
                                   PROP_KNOTS,
                                   g_param_spec_object ("knots",
                                                        NULL,
                                                        "B-spline knots",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmBSpline:control-points:
   *
   * The control points, the coefficients of the basis functions: at least
   * #NcmBSpline:order of them.
   */
  g_object_class_install_property (object_class,
                                   PROP_CONTROL_POINTS,
                                   g_param_spec_object ("control-points",
                                                        NULL,
                                                        "B-spline control points",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_bspline_new:
 * @order: the B-spline order
 * @knots: the knots
 * @control_points: the control points
 *
 * Creates the B-spline of order @order with @knots and @control_points, which it
 * references. Aborts when a vector is not contiguous, when the length of @knots is not
 * that of @control_points plus @order, when there are fewer than @order control points,
 * or when the knots decrease.
 *
 * Returns: (transfer full): a new #NcmBSpline.
 */
NcmBSpline *
ncm_bspline_new (const guint order, NcmVector *knots, NcmVector *control_points)
{
  return g_object_new (NCM_TYPE_BSPLINE,
                       "order", order,
                       "knots", knots,
                       "control-points", control_points,
                       NULL);
}

/**
 * ncm_bspline_ref:
 * @bs: a #NcmBSpline
 *
 * Increases the reference count of @bs by one.
 *
 * Returns: (transfer full): @bs.
 */
NcmBSpline *
ncm_bspline_ref (NcmBSpline *bs)
{
  return g_object_ref (bs);
}

/**
 * ncm_bspline_free:
 * @bs: a #NcmBSpline
 *
 * Decreases the reference count of @bs by one.
 */
void
ncm_bspline_free (NcmBSpline *bs)
{
  g_object_unref (bs);
}

/**
 * ncm_bspline_clear:
 * @bs: a #NcmBSpline
 *
 * Decreases the reference count of *@bs by one and sets *@bs to %NULL; does nothing when
 * *@bs is %NULL.
 */
void
ncm_bspline_clear (NcmBSpline **bs)
{
  g_clear_object (bs);
}

/**
 * ncm_bspline_get_order:
 * @bs: a #NcmBSpline
 *
 * Gets #NcmBSpline:order.
 *
 * Returns: the B-spline order.
 */
guint
ncm_bspline_get_order (NcmBSpline *bs)
{
  return bs->order;
}

/**
 * ncm_bspline_peek_knots:
 * @bs: a #NcmBSpline
 *
 * Gets #NcmBSpline:knots. The vector must not be modified.
 *
 * Returns: (transfer none): the knots.
 */
NcmVector *
ncm_bspline_peek_knots (NcmBSpline *bs)
{
  return bs->knots;
}

/**
 * ncm_bspline_peek_control_points:
 * @bs: a #NcmBSpline
 *
 * Gets #NcmBSpline:control-points. The vector must not be modified.
 *
 * Returns: (transfer none): the control points.
 */
NcmVector *
ncm_bspline_peek_control_points (NcmBSpline *bs)
{
  return bs->control_points;
}

/**
 * ncm_bspline_get_range:
 * @bs: a #NcmBSpline
 * @x_lo: (out): the lower end of the range
 * @x_hi: (out): the upper end of the range
 *
 * Gets the range of @bs: from the knot of index #NcmBSpline:order minus one to the knot
 * of index equal to the number of control points. Outside it ncm_bspline_eval()
 * extrapolates the edge polynomial.
 */
void
ncm_bspline_get_range (NcmBSpline *bs, gdouble *x_lo, gdouble *x_hi)
{
  *x_lo = bs->t[bs->order - 1];
  *x_hi = bs->t[bs->n_control];
}

/**
 * ncm_bspline_eval:
 * @bs: a #NcmBSpline
 * @x: a point
 *
 * Evaluates @bs at @x. No argument is checked. It may be called concurrently on one
 * object.
 *
 * Returns: the value of the B-spline at @x.
 */
gdouble
ncm_bspline_eval (const NcmBSpline *bs, const gdouble x)
{
  return _ncm_bspline_eval (bs->order, bs->n_control, bs->t, bs->c, x);
}


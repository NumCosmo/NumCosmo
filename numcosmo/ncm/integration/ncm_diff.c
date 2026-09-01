/***************************************************************************
 *            ncm_diff.c
 *
 *  Fri July 21 12:59:36 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_diff.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmDiff:
 *
 * Numerical differentiation object.
 *
 * Computes first and second derivatives of user functions by finite
 * differences plus Richardson extrapolation.
 *
 * The extrapolation uses precomputed tables. For a step sequence
 * $g_k$ (powers of #NcmDiff:richardson-step), order $n$ holds the nodes
 * $h_i = 1/g_i$, $i = 0, \dots, n+1$, and the Lagrange weights
 * $\lambda_i = \prod_{j \neq i} 1/(1 - g_j/g_i)$ that extrapolate a
 * polynomial in $h$ (forward) or $h^2$ (central) to $h \to 0$. The
 * extrapolated derivative at order $n$ is $\sum_i \lambda_i D(h_0 h_i)$,
 * where $D(h)$ is the finite-difference quotient at step $h$.
 *
 * Each result carries an error estimate combining the truncation error
 * (difference between consecutive extrapolation orders, padded by
 * #NcmDiff:terr-pad) and the propagated round-off of the difference
 * quotients (padded by #NcmDiff:round-off-pad). Orders are increased
 * while any component still improves; the best value per component is
 * returned.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_diff.h"
#include "ncm/core/ncm_cfg.h"

typedef struct _NcmDiffPrivate
{
  guint maxorder;
  gdouble rs;
  gdouble terr_pad;
  gdouble roff_pad;
  gdouble ini_h;
  gboolean dual_series;
  GPtrArray *central_tables;
  GPtrArray *forward_tables;
  GPtrArray *backward_tables;
} NcmDiffPrivate;

/* Extrapolation nodes h_i = 1/g_i and Lagrange weights lambda_i for one order. */
typedef struct _NcmDiffTable
{
  NcmVector *h;
  NcmVector *lambda;
} NcmDiffTable;

static NcmDiffTable *_ncm_diff_table_new (const guint n);
static void _ncm_diff_table_free (gpointer dtable_ptr);

enum
{
  PROP_0,
  PROP_MAXORDER,
  PROP_RS,
  PROP_ROFF_PAD,
  PROP_TERR_PAD,
  PROP_INI_H,
  PROP_DUAL_SERIES,
  PROP_SIZE,
};

struct _NcmDiff
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmDiff, ncm_diff, G_TYPE_OBJECT)

static void
ncm_diff_init (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  self->maxorder    = 0;
  self->rs          = 0.0;
  self->terr_pad    = 0.0;
  self->roff_pad    = 0.0;
  self->ini_h       = 0.0;
  self->dual_series = FALSE;

  self->central_tables  = g_ptr_array_new ();
  self->forward_tables  = g_ptr_array_new ();
  self->backward_tables = g_ptr_array_new ();

  g_ptr_array_set_free_func (self->central_tables,  &_ncm_diff_table_free);
  g_ptr_array_set_free_func (self->forward_tables,  &_ncm_diff_table_free);
  g_ptr_array_set_free_func (self->backward_tables, &_ncm_diff_table_free);
}

static void
_ncm_diff_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDiff *diff = NCM_DIFF (object);

  g_return_if_fail (NCM_IS_DIFF (object));

  switch (prop_id)
  {
    case PROP_MAXORDER:
      ncm_diff_set_max_order (diff, g_value_get_uint (value));
      break;
    case PROP_RS:
      ncm_diff_set_richardson_step (diff, g_value_get_double (value));
      break;
    case PROP_ROFF_PAD:
      ncm_diff_set_round_off_pad (diff, g_value_get_double (value));
      break;
    case PROP_TERR_PAD:
      ncm_diff_set_trunc_error_pad (diff, g_value_get_double (value));
      break;
    case PROP_INI_H:
      ncm_diff_set_ini_h (diff, g_value_get_double (value));
      break;
    case PROP_DUAL_SERIES:
      ncm_diff_set_dual_series (diff, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_diff_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDiff *diff = NCM_DIFF (object);

  g_return_if_fail (NCM_IS_DIFF (object));

  switch (prop_id)
  {
    case PROP_MAXORDER:
      g_value_set_uint (value, ncm_diff_get_max_order (diff));
      break;
    case PROP_RS:
      g_value_set_double (value, ncm_diff_get_richardson_step (diff));
      break;
    case PROP_ROFF_PAD:
      g_value_set_double (value, ncm_diff_get_round_off_pad (diff));
      break;
    case PROP_TERR_PAD:
      g_value_set_double (value, ncm_diff_get_trunc_error_pad (diff));
      break;
    case PROP_INI_H:
      g_value_set_double (value, ncm_diff_get_ini_h (diff));
      break;
    case PROP_DUAL_SERIES:
      g_value_set_boolean (value, ncm_diff_get_dual_series (diff));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_diff_build_diff_tables (NcmDiff *diff);

static void
_ncm_diff_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_diff_parent_class)->constructed (object);
  {
    NcmDiff *diff = NCM_DIFF (object);


    _ncm_diff_build_diff_tables (diff);
  }
}

static void
_ncm_diff_dispose (GObject *object)
{
  NcmDiff *diff               = NCM_DIFF (object);
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_clear_pointer (&self->central_tables,  g_ptr_array_unref);
  g_clear_pointer (&self->forward_tables,  g_ptr_array_unref);
  g_clear_pointer (&self->backward_tables, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_diff_parent_class)->dispose (object);
}

static void
_ncm_diff_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_diff_parent_class)->finalize (object);
}

static void
ncm_diff_class_init (NcmDiffClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);


  object_class->set_property = &_ncm_diff_set_property;
  object_class->get_property = &_ncm_diff_get_property;
  object_class->constructed  = &_ncm_diff_constructed;
  object_class->dispose      = &_ncm_diff_dispose;
  object_class->finalize     = &_ncm_diff_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MAXORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum order",
                                                      1, G_MAXUINT, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RS,
                                   g_param_spec_double ("richardson-step",
                                                        NULL,
                                                        "Richardson extrapolation step",
                                                        1.1, G_MAXDOUBLE, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ROFF_PAD,
                                   g_param_spec_double ("round-off-pad",
                                                        NULL,
                                                        "Round off padding",
                                                        1.01, G_MAXDOUBLE, 3.0e4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TERR_PAD,
                                   g_param_spec_double ("terr-pad",
                                                        NULL,
                                                        "Truncation error padding",
                                                        1.1, G_MAXDOUBLE, 3.0e4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INI_H,
                                   g_param_spec_double ("ini-h",
                                                        NULL,
                                                        "Initial h",
                                                        GSL_DBL_EPSILON, G_MAXDOUBLE, pow (GSL_DBL_EPSILON, 1.0 / 8.0),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DUAL_SERIES,
                                   g_param_spec_boolean ("dual-series",
                                                         NULL,
                                                         "Use two parallel extrapolation series",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static NcmDiffTable *
_ncm_diff_table_new (const guint n)
{
  NcmDiffTable *dtable = g_new (NcmDiffTable, 1);


  dtable->h      = ncm_vector_new (n);
  dtable->lambda = ncm_vector_new (n);

  return dtable;
}

static void
_ncm_diff_table_free (gpointer dtable_ptr)
{
  NcmDiffTable *dtable = (NcmDiffTable *) dtable_ptr;


  ncm_vector_free (dtable->h);
  ncm_vector_free (dtable->lambda);
  g_free (dtable);
}

static void
_ncm_diff_build_diff_table (NcmDiff *diff, GPtrArray *tables, const guint maxorder, gdouble (*g) (const guint, gpointer), gpointer user_data)
{
  if (tables->len == maxorder)
  {
    return;
  }
  else if (tables->len > maxorder)
  {
    g_ptr_array_set_size (tables, maxorder);
  }
  else
  {
    guint k;


    for (k = tables->len; k < maxorder; k++)
    {
      const guint n = k + 2;


      if (k == 0)
      {
        const gdouble g0     = g (0, user_data);
        const gdouble g1     = g (1, user_data);
        NcmDiffTable *dtable = _ncm_diff_table_new (n);


        ncm_vector_set (dtable->h, 0, 1.0 / g0);
        ncm_vector_set (dtable->h, 1, 1.0 / g1);

        ncm_vector_set (dtable->lambda, 0, g0 / (g0 - g1));
        ncm_vector_set (dtable->lambda, 1, g1 / (g1 - g0));

        g_ptr_array_add (tables, dtable);
      }
      else
      {
        const guint l         = n - 1;
        const gdouble gl      = g (l, user_data);
        NcmDiffTable *ldtable = g_ptr_array_index (tables, k - 1);
        NcmDiffTable *dtable  = _ncm_diff_table_new (n);
        gdouble lambda_l      = 1.0;
        guint i;


        for (i = 0; i < l; i++)
        {
          const gdouble gi       = g (i, user_data);
          const gdouble hi       = ncm_vector_get (ldtable->h, i);
          const gdouble lambda_i = ncm_vector_get (ldtable->lambda, i) / (1.0 - gl / gi);


          ncm_vector_set (dtable->h,      i, hi);
          ncm_vector_set (dtable->lambda, i, lambda_i);

          lambda_l *= 1.0 / (1.0 - gi / gl);
        }

        ncm_vector_set (dtable->h,      l, 1.0 / gl);
        ncm_vector_set (dtable->lambda, l, lambda_l);

        g_ptr_array_add (tables, dtable);
      }
    }
  }
}

static gdouble
_ncm_diff_central_g (const guint k, gpointer user_data)
{
  NcmDiff *diff               = NCM_DIFF (user_data);
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  const gdouble g             = pow (self->rs, 2 * k);


  return +g;
}

static gdouble
_ncm_diff_forward_g (const guint k, gpointer user_data)
{
  NcmDiff *diff               = NCM_DIFF (user_data);
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  const gdouble g             = pow (self->rs, k);


  return +g;
}

static gdouble
_ncm_diff_backward_g (const guint k, gpointer user_data)
{
  NcmDiff *diff               = NCM_DIFF (user_data);
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  const gdouble g             = pow (self->rs, k);


  return -g;
}

static void
_ncm_diff_build_diff_tables (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  _ncm_diff_build_diff_table (diff, self->central_tables,  self->maxorder, _ncm_diff_central_g,  diff);
  _ncm_diff_build_diff_table (diff, self->forward_tables,  self->maxorder, _ncm_diff_forward_g,  diff);
  _ncm_diff_build_diff_table (diff, self->backward_tables, self->maxorder, _ncm_diff_backward_g, diff);
}

/**
 * ncm_diff_new:
 *
 * Creates a new #NcmDiff object.
 *
 * Returns: a new #NcmDiff.
 */
NcmDiff *
ncm_diff_new (void)
{
  NcmDiff *diff = g_object_new (NCM_TYPE_DIFF,
                                NULL);


  return diff;
}

/**
 * ncm_diff_ref:
 * @diff: a #NcmDiff
 *
 * Increase the reference of @diff by one.
 *
 * Returns: (transfer full): @diff.
 */
NcmDiff *
ncm_diff_ref (NcmDiff *diff)
{
  return g_object_ref (diff);
}

/**
 * ncm_diff_free:
 * @diff: a #NcmDiff
 *
 * Decrease the reference count of @diff by one.
 *
 */
void
ncm_diff_free (NcmDiff *diff)
{
  g_object_unref (diff);
}

/**
 * ncm_diff_clear:
 * @diff: a #NcmDiff
 *
 * Decrease the reference count of @diff by one, and sets the pointer *@diff to
 * NULL.
 *
 */
void
ncm_diff_clear (NcmDiff **diff)
{
  g_clear_object (diff);
}

/**
 * ncm_diff_get_max_order:
 * @diff: a #NcmDiff
 *
 * Gets the maximum order used when calculating the derivatives.
 *
 * Returns: the maximum order.
 */
guint
ncm_diff_get_max_order (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->maxorder;
}

/**
 * ncm_diff_get_richardson_step:
 * @diff: a #NcmDiff
 *
 * Gets the current Richardson step used in the tables.
 *
 * Returns: the maximum order.
 */
gdouble
ncm_diff_get_richardson_step (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->rs;
}

/**
 * ncm_diff_get_round_off_pad:
 * @diff: a #NcmDiff
 *
 * Gets the current round-off padding used in calculations.
 *
 * Returns: the round-off padding.
 */
gdouble
ncm_diff_get_round_off_pad (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->roff_pad;
}

/**
 * ncm_diff_get_trunc_error_pad:
 * @diff: a #NcmDiff
 *
 * Gets the current truncation error padding used in calculations.
 *
 * Returns: the truncation error padding.
 */
gdouble
ncm_diff_get_trunc_error_pad (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->terr_pad;
}

/**
 * ncm_diff_get_ini_h:
 * @diff: a #NcmDiff
 *
 * Gets the current initial step used in calculations.
 *
 * Returns: the initial step.
 */
gdouble
ncm_diff_get_ini_h (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->ini_h;
}

/**
 * ncm_diff_set_max_order:
 * @diff: a #NcmDiff
 * @maxorder: the new maximum order
 *
 * Sets the maximum order used when calculating the derivatives to @maxorder.
 *
 */
void
ncm_diff_set_max_order (NcmDiff *diff, const guint maxorder)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpuint (maxorder, >, 0);

  if (maxorder != self->maxorder)
  {
    self->maxorder = maxorder;

    if (self->central_tables->len > 0)
      _ncm_diff_build_diff_tables (diff);
  }
}

/**
 * ncm_diff_set_richardson_step:
 * @diff: a #NcmDiff
 * @rs: the new Richardson step
 *
 * Sets the Richardson step used in the tables.
 *
 */
void
ncm_diff_set_richardson_step (NcmDiff *diff, const gdouble rs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (rs, >, 1.1);

  if (rs != self->rs)
  {
    self->rs = rs;

    if (self->central_tables->len > 0)
      _ncm_diff_build_diff_tables (diff);
  }
}

/**
 * ncm_diff_set_round_off_pad:
 * @diff: a #NcmDiff
 * @roff_pad: the new round-off padding
 *
 * Sets the round-off padding used in the calculations.
 *
 */
void
ncm_diff_set_round_off_pad (NcmDiff *diff, const gdouble roff_pad)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (roff_pad, >, 1.01);
  self->roff_pad = roff_pad;
}

/**
 * ncm_diff_set_trunc_error_pad:
 * @diff: a #NcmDiff
 * @terr_pad: the new truncation error padding
 *
 * Sets the truncation error padding used in the calculations.
 *
 */
void
ncm_diff_set_trunc_error_pad (NcmDiff *diff, const gdouble terr_pad)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (terr_pad, >, 1.01);
  self->terr_pad = terr_pad;
}

/**
 * ncm_diff_set_ini_h:
 * @diff: a #NcmDiff
 * @ini_h: the new initial step
 *
 * Sets the initial step used in the calculations.
 *
 */
void
ncm_diff_set_ini_h (NcmDiff *diff, const gdouble ini_h)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (ini_h, >, GSL_DBL_EPSILON);
  self->ini_h = ini_h;
}

/**
 * ncm_diff_set_dual_series:
 * @diff: a #NcmDiff
 * @dual_series: whether to use two parallel extrapolation series
 *
 * Enables or disables the dual-series scheme. When enabled, every derivative
 * runs two Richardson extrapolation series in parallel, started from the
 * initial steps $h_0$ and $h_0/\sqrt{r_s}$, where $r_s$ is
 * #NcmDiff:richardson-step. The disagreement between the two series at the
 * same order replaces the difference between consecutive orders as the
 * truncation error estimate, which is tighter and does not lag one order
 * behind. The returned derivative comes from the smaller-step series.
 *
 * The scheme roughly doubles the number of function evaluations.
 *
 */
void
ncm_diff_set_dual_series (NcmDiff *diff, const gboolean dual_series)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  self->dual_series = dual_series;
}

/**
 * ncm_diff_get_dual_series:
 * @diff: a #NcmDiff
 *
 * Gets whether the dual-series scheme is enabled, see ncm_diff_set_dual_series().
 *
 * Returns: whether the dual-series scheme is enabled.
 */
gboolean
ncm_diff_get_dual_series (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->dual_series;
}

/**
 * ncm_diff_log_central_tables:
 * @diff: a #NcmDiff
 *
 * Logs all central tables.
 *
 */
void
ncm_diff_log_central_tables (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  guint i;


  for (i = 0; i < self->central_tables->len; i++)
  {
    NcmDiffTable *dtable = g_ptr_array_index (self->central_tables, i);


    ncm_message ("# NcmDiff[central]  order: %u\n", i + 1);
    ncm_vector_log_vals (dtable->h,      "# NcmDiff[central]  h:      ", "% 22.15g", TRUE);
    ncm_vector_log_vals (dtable->lambda, "# NcmDiff[central]  lambda: ", "% 22.15g", TRUE);
  }
}

/**
 * ncm_diff_log_forward_tables:
 * @diff: a #NcmDiff
 *
 * Logs all central tables.
 *
 */
void
ncm_diff_log_forward_tables (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  guint i;


  for (i = 0; i < self->forward_tables->len; i++)
  {
    NcmDiffTable *dtable = g_ptr_array_index (self->forward_tables, i);


    ncm_message ("# NcmDiff[forward]  order: %u\n", i + 1);
    ncm_vector_log_vals (dtable->h,      "# NcmDiff[forward]  h:      ", "% 22.15g", TRUE);
    ncm_vector_log_vals (dtable->lambda, "# NcmDiff[forward]  lambda: ", "% 22.15g", TRUE);
  }
}

/**
 * ncm_diff_log_backward_tables:
 * @diff: a #NcmDiff
 *
 * Logs all central tables.
 *
 */
void
ncm_diff_log_backward_tables (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  guint i;


  for (i = 0; i < self->backward_tables->len; i++)
  {
    NcmDiffTable *dtable = g_ptr_array_index (self->backward_tables, i);


    ncm_message ("# NcmDiff[backward] order: %u\n", i + 1);
    ncm_vector_log_vals (dtable->h,      "# NcmDiff[backward] h:      ", "% 22.15g", TRUE);
    ncm_vector_log_vals (dtable->lambda, "# NcmDiff[backward] lambda: ", "% 22.15g", TRUE);
  }
}

typedef void (*NcmDiffStepAlgo) (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff);
typedef void (*NcmDiffHessianStepAlgo) (NcmDiff *diff, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hx, const guint b, const gdouble y, const gdouble hy, NcmVector *x_v, const gdouble fval, gdouble *df, gdouble *roff);

static void
_ncm_diff_rf_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  /*ncm_vector_log_vals (x_v,  "x_v  ", "% 22.15g", TRUE);*/
  ncm_vector_addto (x_v, a, h);

  /*ncm_vector_log_vals (x_v,  "xh_v ", "% 22.15g", TRUE);*/
  /*ncm_vector_log_vals (f_v,  "f_v  ", "% 22.15g", TRUE);*/
  f (x_v, df, user_data);

  /*ncm_vector_log_vals (yh_v, "yh_v ", "% 22.15g", TRUE);*/

  ncm_vector_memcpy (roff, df);
  ncm_vector_sub_round_off (roff, f_v);

  ncm_vector_sub   (df, f_v);
  ncm_vector_scale (df, 1.0 / h);

  ncm_vector_mul (roff, df);

  NCM_UNUSED (yh1_v);
  NCM_UNUSED (yh2_v);

  /*
   *  printf ("# t = %u\n", t);
   *  ncm_vector_log_vals (df,   "df   ", "% 22.15g", TRUE);
   *  ncm_vector_log_vals (roff, "roff ", "% 22.15g", TRUE);
   */
}

static void
_ncm_diff_rc_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  ncm_vector_memcpy (roff, df);
  ncm_vector_sub_round_off (roff, yh1_v);

  ncm_vector_sub   (df, yh1_v);
  ncm_vector_scale (df, 0.5 / h);

  ncm_vector_mul (roff, df);

  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rc_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  ncm_vector_add (df, yh1_v);
  ncm_vector_scale (df, 0.5);

  ncm_vector_memcpy (roff, df);
  ncm_vector_sub_round_off (roff, f_v);

  ncm_vector_sub   (df, f_v);
  ncm_vector_scale (df, 2.0 / (h * h));

  ncm_vector_mul (roff, df);

  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rf_Hessian_step (NcmDiff *diff, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hx, const guint b, const gdouble y, const gdouble hy, NcmVector *x_v, const gdouble fval, gdouble *df, gdouble *roff)
{
  gdouble f_hx, f_hy, f_hxhy;


  ncm_vector_addto (x_v, a, hx);
  f_hx = f (x_v, user_data);

  ncm_vector_set (x_v, a, x);
  ncm_vector_addto (x_v, b, hy);
  f_hy = f (x_v, user_data);

  ncm_vector_addto (x_v, a, hx);
  f_hxhy = f (x_v, user_data);

  df[0] = ((fval + f_hxhy) - (f_hx + f_hy)) / (hx * hy);

  {
    const gdouble s1 = fval + f_hxhy;
    const gdouble s2 = f_hx + f_hy;

    const gdouble d1 = s1 - s2;


    if (G_UNLIKELY (d1 == 0.0))
    {
      roff[0] = 1.0;
    }
    else
    {
      const gdouble abs_s1  = fabs (s1);
      const gdouble abs_s2  = fabs (s2);
      const gdouble max_s12 = GSL_MAX (abs_s1, abs_s2);


      roff[0] = fabs (max_s12 * GSL_DBL_EPSILON / d1);
    }
  }
}

#define NCM_DIFF_ERR_PAD (1.0e0)
#define NCM_DIFF_NTRY_CONV (3)

/*
 * Per-component convergence tracker, shared by the vector and Hessian
 * drivers. At each extrapolation order it receives the current and previous
 * extrapolated values with their round-off estimates, updates the best value
 * and error seen so far, and reports whether this component still wants a
 * higher order.
 */
typedef struct _NcmDiffConv
{
  gdouble df_best;
  gdouble err_best;
  gdouble err_last_max;
  guchar not_conv;
  guchar cstarted;
} NcmDiffConv;

static void
_ncm_diff_conv_init (NcmDiffConv *cs)
{
  cs->df_best      = 0.0;
  cs->err_best     = GSL_POSINF;
  cs->err_last_max = 0.0;
  cs->not_conv     = NCM_DIFF_NTRY_CONV;
  cs->cstarted     = 0;
}

/* Relative difference used for the convergence test, with the same zero
 * handling as ncm_vector_cmp(). */
static gdouble
_ncm_diff_rel_diff (const gdouble x1, const gdouble x2)
{
  if (G_UNLIKELY (x1 == 0.0))
    return (x2 == 0.0) ? 0.0 : fabs (x2);
  else if (G_UNLIKELY (x2 == 0.0))
    return fabs (x1);
  else
    return fabs ((x1 - x2) / GSL_MIN (fabs (x1), fabs (x2)));
}

static gboolean
_ncm_diff_conv_update (NcmDiffConv *cs, const gdouble terr_pad, const gdouble roff_pad,
                       const gdouble df_curr, const gdouble df_last,
                       const gdouble roff_curr, const gdouble roff_last)
{
  const gdouble err_trunc    = fabs (df_curr - df_last) * terr_pad;
  const gdouble err_err      = _ncm_diff_rel_diff (df_curr, df_last);
  const gdouble Eroff_last   = fabs (roff_last) * roff_pad;
  const gdouble Eroff_curr   = fabs (roff_curr) * roff_pad;
  const gdouble err_curr_max = GSL_MAX (err_trunc, GSL_MAX (Eroff_last, Eroff_curr));
  gdouble err_curr_best      = cs->err_best;
  gboolean improve           = FALSE;

  /*
   * Estimates fluctuate in the beginning. Thus we only start checking for
   * convergence after they agree at 1.0e-3.
   */
  if (cs->not_conv && (err_err < 1.0e-3))
  {
    cs->not_conv--;

    if (!cs->not_conv)
      cs->cstarted = 1;
  }
  else if (err_err > 1.0e-3)
  {
    cs->not_conv = NCM_DIFF_NTRY_CONV;
  }

  /*
   * If the current maximum error is smaller than the best error estimate
   * improve it again. It also sets the best error estimate to the average of
   * the last two to avoid fake convergence due to fluctuations on error
   * estimates.
   */
  if (cs->cstarted || ((err_curr_max < cs->err_best) && !cs->not_conv))
  {
    cs->df_best  = df_curr;
    cs->err_best = 0.5 * (err_curr_max + cs->err_last_max);

    err_curr_best = cs->err_best;
    improve       = TRUE;
  }
  else
  {
    const gdouble rel_error      = fabs (err_curr_max / df_curr);
    const gdouble best_rel_error = fabs (cs->err_best / cs->df_best);

    if (rel_error < best_rel_error)
    {
      cs->df_best  = df_curr;
      cs->err_best = err_curr_max;

      err_curr_best = err_curr_max;
    }
  }

  if (cs->not_conv || (Eroff_curr < err_curr_best))
    improve = TRUE;

  cs->cstarted     = 0;
  cs->err_last_max = err_curr_max;

  return improve;
}

/*
 * Per-component tracker for the dual-series scheme. The two series A and B
 * share the extrapolation tables but start from different initial steps, so
 * their disagreement at the same order estimates the truncation error
 * directly. The error also includes the propagated round-off of both series.
 * A component asks to stop after NCM_DIFF_DUAL_MAX_WORSE orders without
 * improvement.
 */
#define NCM_DIFF_DUAL_ERR_PAD (10.0)
#define NCM_DIFF_DUAL_MAX_WORSE (2)
#define NCM_DIFF_DUAL_MIN_ORDER (3)

typedef struct _NcmDiffDualConv
{
  gdouble df_best;
  gdouble err_best;
  guchar n_worse;
} NcmDiffDualConv;

static void
_ncm_diff_dual_conv_init (NcmDiffDualConv *cs)
{
  cs->df_best  = 0.0;
  cs->err_best = GSL_POSINF;
  cs->n_worse  = 0;
}

static gboolean
_ncm_diff_dual_conv_update (NcmDiffDualConv *cs, const gdouble df_A, const gdouble df_B,
                            const gdouble roff_A, const gdouble roff_B)
{
  const gdouble cross = fabs (df_A - df_B);
  const gdouble roff  = GSL_MAX (fabs (roff_A), fabs (roff_B));
  const gdouble err   = GSL_MAX (cross, roff) * NCM_DIFF_DUAL_ERR_PAD;

  if (err < cs->err_best)
  {
    /* B (smaller steps) has the smaller truncation error; when round-off
     * dominates the disagreement, the ladder with less round-off wins. */
    if (fabs (roff_B) <= cross)
      cs->df_best = df_B;
    else
      cs->df_best = (fabs (roff_A) < fabs (roff_B)) ? df_A : df_B;

    cs->err_best = err;
    cs->n_worse  = 0;

    return TRUE;
  }

  /*
   * While the round-off is still far below the best error the series has not
   * reached its floor: the current stagnation is pre-asymptotic (e.g. steps
   * still larger than the scale of variation of f), so keep refining.
   */
  if (roff * NCM_DIFF_DUAL_ERR_PAD < cs->err_best)
  {
    cs->n_worse = 0;

    return TRUE;
  }

  cs->n_worse++;

  return cs->n_worse < NCM_DIFF_DUAL_MAX_WORSE;
}

/* Evaluates one finite-difference step and appends the difference quotient
 * and its round-off estimate to dfs/roffs. */
static void
_ncm_diff_eval_step (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffFuncNtoM f, gpointer user_data,
                     const guint a, const gdouble x, const gdouble h0, const gdouble ht, const guint po,
                     NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v,
                     const guint dim, GPtrArray *dfs, GPtrArray *roffs)
{
  const gdouble ho      = h0 * ((po == 0) ? ht : sqrt (ht));
  volatile gdouble temp = x + ho;
  const gdouble h       = temp - x;
  NcmVector *df_t       = ncm_vector_new (dim);
  NcmVector *roff_t     = ncm_vector_new (dim);

  step_algo (diff, f, user_data, a, x, h, x_v, f_v, yh1_v, yh2_v, df_t, roff_t);

  g_ptr_array_add (dfs,   df_t);
  g_ptr_array_add (roffs, roff_t);

  ncm_vector_set (x_v, a, x);
}

/* Richardson extrapolation at the order of dtable: df = sum_i lambda_i dfs[i],
 * with the round-off estimates combined in quadrature. */
static void
_ncm_diff_accum_extrap (NcmDiffTable *dtable, const guint nt, GPtrArray *dfs, GPtrArray *roffs,
                        NcmVector *df_curr, NcmVector *roff_curr)
{
  guint i;

  ncm_vector_memcpy (df_curr, g_ptr_array_index (dfs, 0));
  ncm_vector_scale (df_curr, ncm_vector_get (dtable->lambda, 0));

  ncm_vector_memcpy (roff_curr, g_ptr_array_index (roffs, 0));
  ncm_vector_scale (roff_curr, ncm_vector_get (dtable->lambda, 0));

  for (i = 1; i < nt; i++)
  {
    ncm_vector_axpy (df_curr, ncm_vector_get (dtable->lambda, i), g_ptr_array_index (dfs, i));
    ncm_vector_hypot (roff_curr, ncm_vector_get (dtable->lambda, i), g_ptr_array_index (roffs, i));
  }
}

static void
_ncm_diff_accum_extrap_s (NcmDiffTable *dtable, const guint nt, GArray *dfs, GArray *roffs,
                          gdouble *df, gdouble *roff)
{
  gdouble df_curr   = 0.0;
  gdouble roff_curr = 0.0;
  guint i;

  for (i = 0; i < nt; i++)
  {
    const gdouble lambda_i = ncm_vector_get (dtable->lambda, i);

    df_curr  += lambda_i * g_array_index (dfs, gdouble, i);
    roff_curr = hypot (roff_curr, lambda_i * g_array_index (roffs, gdouble, i));
  }

  df[0]   = df_curr;
  roff[0] = roff_curr;
}

static GArray *
_ncm_diff_by_step_algo_single (NcmDiff *diff, NcmDiffStepAlgo step_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GPtrArray *tables           = (po == 0) ? self->forward_tables : self->central_tables;
  GPtrArray *dfs              = g_ptr_array_new ();
  GPtrArray *roffs            = g_ptr_array_new ();
  NcmVector *x_v              = NULL;
  NcmVector *f_v              = NULL;
  NcmVector *yh1_v            = NULL;
  NcmVector *yh2_v            = NULL;
  NcmVector *df_last          = NULL;
  NcmVector *df_curr          = NULL;
  NcmVector *roff_last        = NULL;
  NcmVector *roff_curr        = NULL;
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *conv                = g_array_new (FALSE, FALSE, sizeof (NcmDiffConv));
  const guint nvar            = x_a->len;
  NcmMatrix *Eerr_m           = NULL;
  NcmMatrix *df_m;
  guint a;


  g_array_set_size (df, dim * nvar);
  df_m = ncm_matrix_new_array (df, dim);

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, dim * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, dim);
  }

  g_ptr_array_set_free_func (dfs,   (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (roffs, (GDestroyNotify) ncm_vector_free);

  g_array_set_size (conv, dim);

  x_v   = ncm_vector_new_array (x_a);
  f_v   = ncm_vector_new (dim);
  yh1_v = ncm_vector_new (dim);
  yh2_v = ncm_vector_new (dim);

  df_last = ncm_vector_new (dim);
  df_curr = ncm_vector_new (dim);

  roff_last = ncm_vector_new (dim);
  roff_curr = ncm_vector_new (dim);

  f (x_v, f_v, user_data);

  for (a = 0; a < nvar; a++)
  {
    const gdouble x       = g_array_index (x_a, gdouble, a);
    const gdouble scale   = (x == 0.0) ? 1.0 : fabs (x);
    const gdouble h0      = self->ini_h * scale;
    NcmDiffTable *ldtable = NULL;
    guint order_index;
    guint i;
    guint t = 0;


    g_ptr_array_set_size (dfs, 0);
    g_ptr_array_set_size (roffs, 0);

    for (i = 0; i < dim; i++)
      _ncm_diff_conv_init (&g_array_index (conv, NcmDiffConv, i));

    for (order_index = 0; order_index < self->maxorder; order_index++)
    {
      const guint nt       = order_index + 2;
      NcmDiffTable *dtable = g_ptr_array_index (tables, order_index);


      for ( ; t < nt; t++)
        _ncm_diff_eval_step (diff, step_algo, f, user_data, a, x, h0, ncm_vector_get (dtable->h, t), po,
                             x_v, f_v, yh1_v, yh2_v, dim, dfs, roffs);

      if (ldtable == NULL)
      {
        ncm_vector_memcpy (df_last, g_ptr_array_index (dfs, 0));
        ncm_vector_memcpy (roff_last, g_ptr_array_index (roffs, 0));

        for (i = 0; i < dim; i++)
          g_array_index (conv, NcmDiffConv, i).df_best = ncm_vector_get (df_last, i);
      }

      _ncm_diff_accum_extrap (dtable, nt, dfs, roffs, df_curr, roff_curr);

      {
        gboolean improve = FALSE;

        for (i = 0; i < dim; i++)
        {
          NcmDiffConv *cs = &g_array_index (conv, NcmDiffConv, i);

          if (_ncm_diff_conv_update (cs, self->terr_pad, self->roff_pad,
                                     ncm_vector_get (df_curr, i), ncm_vector_get (df_last, i),
                                     ncm_vector_get (roff_curr, i), ncm_vector_get (roff_last, i)))
            improve = TRUE;
        }

        if (!improve)
          break;
      }

      ncm_vector_memcpy (df_last, df_curr);
      ncm_vector_memcpy (roff_last, roff_curr);
      ldtable = dtable;
    }

    for (i = 0; i < dim; i++)
    {
      const NcmDiffConv *cs = &g_array_index (conv, NcmDiffConv, i);

      ncm_matrix_set (df_m, a, i, cs->df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, i, cs->err_best);
    }
  }

  if (Eerr_m != NULL)
    ncm_matrix_scale (Eerr_m, NCM_DIFF_ERR_PAD);

  {
    g_array_unref (conv);

    g_ptr_array_unref (dfs);
    g_ptr_array_unref (roffs);

    ncm_vector_clear (&x_v);
    ncm_vector_clear (&f_v);
    ncm_vector_clear (&yh1_v);
    ncm_vector_clear (&yh2_v);

    ncm_vector_clear (&df_last);
    ncm_vector_clear (&df_curr);

    ncm_vector_clear (&roff_last);
    ncm_vector_clear (&roff_curr);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static GArray *
_ncm_diff_by_step_algo_dual (NcmDiff *diff, NcmDiffStepAlgo step_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GPtrArray *tables           = (po == 0) ? self->forward_tables : self->central_tables;
  GPtrArray *dfs_A            = g_ptr_array_new ();
  GPtrArray *roffs_A          = g_ptr_array_new ();
  GPtrArray *dfs_B            = g_ptr_array_new ();
  GPtrArray *roffs_B          = g_ptr_array_new ();
  NcmVector *x_v              = NULL;
  NcmVector *f_v              = NULL;
  NcmVector *yh1_v            = NULL;
  NcmVector *yh2_v            = NULL;
  NcmVector *df_A             = NULL;
  NcmVector *df_B             = NULL;
  NcmVector *roff_A           = NULL;
  NcmVector *roff_B           = NULL;
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *conv                = g_array_new (FALSE, FALSE, sizeof (NcmDiffDualConv));
  const guint nvar            = x_a->len;
  NcmMatrix *Eerr_m           = NULL;
  NcmMatrix *df_m;
  guint a;


  g_array_set_size (df, dim * nvar);
  df_m = ncm_matrix_new_array (df, dim);

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, dim * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, dim);
  }

  g_ptr_array_set_free_func (dfs_A,   (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (roffs_A, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (dfs_B,   (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (roffs_B, (GDestroyNotify) ncm_vector_free);

  g_array_set_size (conv, dim);

  x_v   = ncm_vector_new_array (x_a);
  f_v   = ncm_vector_new (dim);
  yh1_v = ncm_vector_new (dim);
  yh2_v = ncm_vector_new (dim);

  df_A = ncm_vector_new (dim);
  df_B = ncm_vector_new (dim);

  roff_A = ncm_vector_new (dim);
  roff_B = ncm_vector_new (dim);

  f (x_v, f_v, user_data);

  for (a = 0; a < nvar; a++)
  {
    const gdouble x     = g_array_index (x_a, gdouble, a);
    const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);
    const gdouble h0_A  = self->ini_h * scale;
    const gdouble h0_B  = h0_A / sqrt (self->rs);
    guint order_index;
    guint i;
    guint t = 0;


    g_ptr_array_set_size (dfs_A, 0);
    g_ptr_array_set_size (roffs_A, 0);
    g_ptr_array_set_size (dfs_B, 0);
    g_ptr_array_set_size (roffs_B, 0);

    for (i = 0; i < dim; i++)
      _ncm_diff_dual_conv_init (&g_array_index (conv, NcmDiffDualConv, i));

    for (order_index = 0; order_index < self->maxorder; order_index++)
    {
      const guint nt       = order_index + 2;
      NcmDiffTable *dtable = g_ptr_array_index (tables, order_index);


      for ( ; t < nt; t++)
      {
        const gdouble ht = ncm_vector_get (dtable->h, t);

        _ncm_diff_eval_step (diff, step_algo, f, user_data, a, x, h0_A, ht, po,
                             x_v, f_v, yh1_v, yh2_v, dim, dfs_A, roffs_A);
        _ncm_diff_eval_step (diff, step_algo, f, user_data, a, x, h0_B, ht, po,
                             x_v, f_v, yh1_v, yh2_v, dim, dfs_B, roffs_B);
      }

      _ncm_diff_accum_extrap (dtable, nt, dfs_A, roffs_A, df_A, roff_A);
      _ncm_diff_accum_extrap (dtable, nt, dfs_B, roffs_B, df_B, roff_B);

      {
        gboolean improve = FALSE;

        for (i = 0; i < dim; i++)
        {
          NcmDiffDualConv *cs = &g_array_index (conv, NcmDiffDualConv, i);

          if (_ncm_diff_dual_conv_update (cs,
                                          ncm_vector_get (df_A, i), ncm_vector_get (df_B, i),
                                          ncm_vector_get (roff_A, i), ncm_vector_get (roff_B, i)))
            improve = TRUE;
        }

        if ((order_index + 1 >= NCM_DIFF_DUAL_MIN_ORDER) && !improve)
          break;
      }
    }

    for (i = 0; i < dim; i++)
    {
      const NcmDiffDualConv *cs = &g_array_index (conv, NcmDiffDualConv, i);

      ncm_matrix_set (df_m, a, i, cs->df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, i, cs->err_best);
    }
  }

  {
    g_array_unref (conv);

    g_ptr_array_unref (dfs_A);
    g_ptr_array_unref (roffs_A);
    g_ptr_array_unref (dfs_B);
    g_ptr_array_unref (roffs_B);

    ncm_vector_clear (&x_v);
    ncm_vector_clear (&f_v);
    ncm_vector_clear (&yh1_v);
    ncm_vector_clear (&yh2_v);

    ncm_vector_clear (&df_A);
    ncm_vector_clear (&df_B);

    ncm_vector_clear (&roff_A);
    ncm_vector_clear (&roff_B);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static GArray *
ncm_diff_by_step_algo (NcmDiff *diff, NcmDiffStepAlgo step_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  if (self->dual_series)
    return _ncm_diff_by_step_algo_dual (diff, step_algo, po, x_a, dim, f, user_data, Eerr);

  return _ncm_diff_by_step_algo_single (diff, step_algo, po, x_a, dim, f, user_data, Eerr);
}

/* Evaluates one Hessian cross-term step and appends the difference quotient
 * and its round-off estimate to dfs/roffs. */
static void
_ncm_diff_hessian_eval_step (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, NcmDiffFuncNto1 f, gpointer user_data,
                             const guint a, const gdouble x, const gdouble hx0,
                             const guint b, const gdouble y, const gdouble hy0,
                             const gdouble ht, const guint po, NcmVector *x_v, const gdouble fval,
                             GArray *dfs, GArray *roffs)
{
  const gdouble hto    = (po == 0) ? ht : sqrt (ht);
  const gdouble hxo    = hx0 * hto;
  const gdouble hyo    = hy0 * hto;
  volatile gdouble t_x = x + hxo;
  const gdouble hx     = t_x - x;
  volatile gdouble t_y = y + hyo;
  const gdouble hy     = t_y - y;
  gdouble df_t         = 0.0;
  gdouble roff_t       = 0.0;

  Hstep_algo (diff, f, user_data, a, x, hx, b, y, hy, x_v, fval, &df_t, &roff_t);

  g_array_append_val (dfs,   df_t);
  g_array_append_val (roffs, roff_t);

  ncm_vector_set (x_v, a, x);
  ncm_vector_set (x_v, b, y);
}

static GArray *
_ncm_diff_Hessian_by_step_algo_single (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  GPtrArray *tables = NULL;
  GArray *dfs       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *roffs     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v    = NULL;
  GArray *df        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar  = x_a->len;
  gdouble fval      = 0.0;
  NcmMatrix *Eerr_m = NULL;
  NcmMatrix *df_m;
  guint a;


  g_array_set_size (df, nvar * nvar);
  df_m = ncm_matrix_new_array (df, nvar);

  if (po == 0)
    tables = self->forward_tables;
  else
    tables = self->central_tables;

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, nvar * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, nvar);
  }

  x_v = ncm_vector_new_array (x_a);

  fval = f (x_v, user_data);

  for (a = 0; a < x_a->len; a++)
  {
    guint b;


    for (b = a + 1; b < x_a->len; b++)
    {
      const gdouble x       = g_array_index (x_a, gdouble, a);
      const gdouble y       = g_array_index (x_a, gdouble, b);
      const gdouble scale_x = (x == 0.0) ? 1.0 : fabs (x);
      const gdouble scale_y = (y == 0.0) ? 1.0 : fabs (y);
      const gdouble hx0     = self->ini_h * scale_x;
      const gdouble hy0     = self->ini_h * scale_y;
      NcmDiffTable *ldtable = NULL;
      gdouble df_last       = 0.0;
      gdouble df_curr       = 0.0;
      gdouble roff_last     = 0.0;
      gdouble roff_curr     = 0.0;
      guint t               = 0;
      NcmDiffConv cs;
      guint order_index;


      _ncm_diff_conv_init (&cs);

      g_array_set_size (dfs, 0);
      g_array_set_size (roffs, 0);

      for (order_index = 0; order_index < self->maxorder; order_index++)
      {
        const guint nt       = order_index + 2;
        NcmDiffTable *dtable = g_ptr_array_index (tables, order_index);


        for ( ; t < nt; t++)
          _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data, a, x, hx0, b, y, hy0,
                                       ncm_vector_get (dtable->h, t), po, x_v, fval, dfs, roffs);

        if (ldtable == NULL)
        {
          df_last   = g_array_index (dfs, gdouble, 0);
          roff_last = fabs (g_array_index (roffs, gdouble, 0));
        }

        _ncm_diff_accum_extrap_s (dtable, nt, dfs, roffs, &df_curr, &roff_curr);

        if (ldtable == NULL)
          cs.df_best = df_curr;

        if (!_ncm_diff_conv_update (&cs, self->terr_pad, self->roff_pad,
                                    df_curr, df_last, roff_curr, roff_last))
          break;

        df_last   = df_curr;
        roff_last = roff_curr;
        ldtable   = dtable;
      }

      ncm_matrix_set (df_m, a, b, cs.df_best);
      ncm_matrix_set (df_m, b, a, cs.df_best);

      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, cs.err_best);
        ncm_matrix_set (Eerr_m, b, a, cs.err_best);
      }
    }
  }

  if (Eerr_m != NULL)
    ncm_matrix_scale (Eerr_m, NCM_DIFF_ERR_PAD);

  {
    g_array_unref (dfs);
    g_array_unref (roffs);

    ncm_vector_clear (&x_v);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static GArray *
_ncm_diff_Hessian_by_step_algo_dual (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  GPtrArray *tables = (po == 0) ? self->forward_tables : self->central_tables;
  GArray *dfs_A     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *roffs_A   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *dfs_B     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *roffs_B   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v    = NULL;
  GArray *df        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar  = x_a->len;
  gdouble fval      = 0.0;
  NcmMatrix *Eerr_m = NULL;
  NcmMatrix *df_m;
  guint a;


  g_array_set_size (df, nvar * nvar);
  df_m = ncm_matrix_new_array (df, nvar);

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, nvar * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, nvar);
  }

  x_v = ncm_vector_new_array (x_a);

  fval = f (x_v, user_data);

  for (a = 0; a < nvar; a++)
  {
    guint b;


    for (b = a + 1; b < nvar; b++)
    {
      const gdouble x       = g_array_index (x_a, gdouble, a);
      const gdouble y       = g_array_index (x_a, gdouble, b);
      const gdouble scale_x = (x == 0.0) ? 1.0 : fabs (x);
      const gdouble scale_y = (y == 0.0) ? 1.0 : fabs (y);
      const gdouble hx0_A   = self->ini_h * scale_x;
      const gdouble hy0_A   = self->ini_h * scale_y;
      const gdouble hx0_B   = hx0_A / sqrt (self->rs);
      const gdouble hy0_B   = hy0_A / sqrt (self->rs);
      gdouble df_A          = 0.0;
      gdouble df_B          = 0.0;
      gdouble roff_A        = 0.0;
      gdouble roff_B        = 0.0;
      guint t               = 0;
      NcmDiffDualConv cs;
      guint order_index;


      _ncm_diff_dual_conv_init (&cs);

      g_array_set_size (dfs_A, 0);
      g_array_set_size (roffs_A, 0);
      g_array_set_size (dfs_B, 0);
      g_array_set_size (roffs_B, 0);

      for (order_index = 0; order_index < self->maxorder; order_index++)
      {
        const guint nt       = order_index + 2;
        NcmDiffTable *dtable = g_ptr_array_index (tables, order_index);


        for ( ; t < nt; t++)
        {
          const gdouble ht = ncm_vector_get (dtable->h, t);

          _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data, a, x, hx0_A, b, y, hy0_A,
                                       ht, po, x_v, fval, dfs_A, roffs_A);
          _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data, a, x, hx0_B, b, y, hy0_B,
                                       ht, po, x_v, fval, dfs_B, roffs_B);
        }

        _ncm_diff_accum_extrap_s (dtable, nt, dfs_A, roffs_A, &df_A, &roff_A);
        _ncm_diff_accum_extrap_s (dtable, nt, dfs_B, roffs_B, &df_B, &roff_B);

        if (!_ncm_diff_dual_conv_update (&cs, df_A, df_B, roff_A, roff_B) &&
            (order_index + 1 >= NCM_DIFF_DUAL_MIN_ORDER))
          break;
      }

      ncm_matrix_set (df_m, a, b, cs.df_best);
      ncm_matrix_set (df_m, b, a, cs.df_best);

      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, cs.err_best);
        ncm_matrix_set (Eerr_m, b, a, cs.err_best);
      }
    }
  }

  {
    g_array_unref (dfs_A);
    g_array_unref (roffs_A);
    g_array_unref (dfs_B);
    g_array_unref (roffs_B);

    ncm_vector_clear (&x_v);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static GArray *
ncm_diff_Hessian_by_step_algo (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  if (self->dual_series)
    return _ncm_diff_Hessian_by_step_algo_dual (diff, Hstep_algo, po, x_a, f, user_data, Eerr);

  return _ncm_diff_Hessian_by_step_algo_single (diff, Hstep_algo, po, x_a, f, user_data, Eerr);
}

/**
 * ncm_diff_rf_d1_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N\to \mathbb{R}^M$,
 * where $N = $ length of @x_a and $M = $ @dim.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rf_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, dim, f, user_data, Eerr);
}

/**
 * ncm_diff_rc_d1_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N\to \mathbb{R}^M$,
 * where $N = $ length of @x_a and $M = $ @dim.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, dim, f, user_data, Eerr);
}

/**
 * ncm_diff_rc_d2_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N\to \mathbb{R}^M$,
 * where $N = $ length of @x_a and $M = $ @dim.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d2_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, dim, f, user_data, Eerr);
}

typedef struct _NcmDiffFuncParams
{
  NcmDiffFunc1toM f_1_to_M;
  NcmDiffFuncNto1 f_N_to_1;
  NcmDiffFunc1to1 f_1_to_1;
  gpointer user_data;
} NcmDiffFuncParams;

static void
_ncm_diff_trans_1_to_M (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDiffFuncParams *fp = (NcmDiffFuncParams *) user_data;


  fp->f_1_to_M (ncm_vector_get (x, 0), y, fp->user_data);
}

static void
_ncm_diff_trans_N_to_1 (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDiffFuncParams *fp = (NcmDiffFuncParams *) user_data;


  ncm_vector_set (y, 0, fp->f_N_to_1 (x, fp->user_data));
}

static void
_ncm_diff_trans_1_to_1 (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDiffFuncParams *fp = (NcmDiffFuncParams *) user_data;


  ncm_vector_set (y, 0, fp->f_1_to_1 (ncm_vector_get (x, 0), fp->user_data));
}

/**
 * ncm_diff_rf_d1_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rf_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;


  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a =  ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);

  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_rc_d1_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;


  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);
  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_rc_d2_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d2_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;


  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);
  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_rf_d1_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rf_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};


  return ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_rc_d1_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};


  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_rc_d2_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_rc_d2_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};


  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_rf_Hessian_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the Hessian of @f $\partial_i\partial_j f$ using the forward method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R}^N \to \mathbb{R}$,
 * where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The Hessian of @f at @x_a.
 */
GArray *
ncm_diff_rf_Hessian_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};
  GArray *dEerr        = NULL;

  GArray *diag = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, &dEerr);
  GArray *res  = ncm_diff_Hessian_by_step_algo (diff, _ncm_diff_rf_Hessian_step, 0, x_a, f, user_data, Eerr);

  guint i;


  g_assert_cmpuint (diag->len * diag->len, ==, res->len);

  for (i = 0; i < diag->len; i++)
  {
    g_array_index (res, gdouble, i * diag->len + i) = g_array_index (diag, gdouble, i);

    if (Eerr != NULL)
      g_array_index (*Eerr, gdouble, i * diag->len + i) = g_array_index (dEerr, gdouble, i);
  }

  g_array_unref (dEerr);
  g_array_unref (diag);

  return res;
}

/**
 * ncm_diff_rf_d1_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the forward method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R} \to \mathbb{R}$.
 *
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_rf_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr         = NULL;
  GArray *df_a;
  gdouble df;


  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, 0, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);

  g_array_unref (Eerr);

  return df;
}

/**
 * ncm_diff_rc_d1_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R} \to \mathbb{R}$.
 *
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_rc_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr         = NULL;
  GArray *df_a;
  gdouble df;


  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);

  g_array_unref (Eerr);

  return df;
}

/**
 * ncm_diff_rc_d2_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 *
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the central method plus
 * Richardson extrapolation. The function $f$ is considered as a $f:\mathbb{R} \to \mathbb{R}$.
 *
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_rc_d2_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr         = NULL;
  GArray *df_a;
  gdouble df;


  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);

  g_array_unref (Eerr);

  return df;
}


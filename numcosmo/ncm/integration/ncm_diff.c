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
 * quotients (padded by #NcmDiff:round-off-pad). Two orders agree when
 * they differ by less than a thousandth of the largest difference quotient
 * of the row, so a zero derivative converges like any other, and only
 * while the round-off of the row is below that level. Orders are increased
 * while any component still improves; the best value per component is
 * returned.
 *
 * The initial step is #NcmDiff:ini-h times the coordinate (times one at a
 * zero coordinate). When it is far below the scale on which $f$ varies the
 * first quotients are cancellation noise, and the step is moved up: first
 * to the step a zero coordinate gets, then by factors of
 * #NcmDiff:richardson-step, up to a few factors above it. No point of the
 * domain is treated as an edge, so near one the caller chooses the forward
 * scheme, which steps to one side only, or a smaller #NcmDiff:ini-h.
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
  gdouble spectral_window;
  NcmVector *lb;
  NcmVector *ub;
  gboolean domain_warnings;
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
  PROP_SPECTRAL_WINDOW,
  PROP_DOMAIN_WARNINGS,
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

  self->maxorder        = 0;
  self->rs              = 0.0;
  self->terr_pad        = 0.0;
  self->roff_pad        = 0.0;
  self->ini_h           = 0.0;
  self->dual_series     = FALSE;
  self->spectral_window = 0.0;
  self->lb              = NULL;
  self->ub              = NULL;
  self->domain_warnings = TRUE;

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
    case PROP_SPECTRAL_WINDOW:
      ncm_diff_set_spectral_window (diff, g_value_get_double (value));
      break;
    case PROP_DOMAIN_WARNINGS:
      ncm_diff_set_domain_warnings (diff, g_value_get_boolean (value));
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
    case PROP_SPECTRAL_WINDOW:
      g_value_set_double (value, ncm_diff_get_spectral_window (diff));
      break;
    case PROP_DOMAIN_WARNINGS:
      g_value_set_boolean (value, ncm_diff_get_domain_warnings (diff));
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

  ncm_vector_clear (&self->lb);
  ncm_vector_clear (&self->ub);

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
  g_object_class_install_property (object_class,
                                   PROP_SPECTRAL_WINDOW,
                                   g_param_spec_double ("spectral-window",
                                                        NULL,
                                                        "Initial spectral window half-width in units of the variable scale",
                                                        GSL_DBL_EPSILON, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DOMAIN_WARNINGS,
                                   g_param_spec_boolean ("domain-warnings",
                                                         NULL,
                                                         "Warn when a central difference falls back to a one-sided one at an edge of the domain",
                                                         TRUE,
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
 * ncm_diff_set_spectral_window:
 * @diff: a #NcmDiff
 * @spectral_window: the new initial spectral window half-width
 *
 * Sets the initial half-width of the window used by the spectral methods
 * (ncm_diff_sc_d1_N_to_M() and related), in units of the variable scale.
 * The window search starts from this value and expands or shrinks it to
 * match the scale of variation of the function.
 *
 */
void
ncm_diff_set_spectral_window (NcmDiff *diff, const gdouble spectral_window)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (spectral_window, >, 0.0);
  self->spectral_window = spectral_window;
}

/**
 * ncm_diff_set_domain:
 * @diff: a #NcmDiff
 * @lb: (nullable): lower bounds of the coordinates, $-\infty$ where there is none
 * @ub: (nullable): upper bounds of the coordinates, $+\infty$ where there is none
 *
 * Sets the domain of the functions to differentiate: no point is evaluated
 * outside $[\mathrm{lb}, \mathrm{ub}]$. Coordinates are matched by index; the
 * scalar functions use entry zero. When a central difference cannot keep its
 * points inside the domain it falls back to a forward or backward one toward
 * the open side, with a warning unless #NcmDiff:domain-warnings is off.
 * Passing %NULL for both clears the domain.
 */
void
ncm_diff_set_domain (NcmDiff *diff, NcmVector *lb, NcmVector *ub)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  if ((lb != NULL) && (ub != NULL))
    g_assert_cmpuint (ncm_vector_len (lb), ==, ncm_vector_len (ub));

  ncm_vector_clear (&self->lb);
  ncm_vector_clear (&self->ub);

  self->lb = (lb != NULL) ? ncm_vector_ref (lb) : NULL;
  self->ub = (ub != NULL) ? ncm_vector_ref (ub) : NULL;
}

/**
 * ncm_diff_clear_domain:
 * @diff: a #NcmDiff
 *
 * Removes the domain set by ncm_diff_set_domain().
 */
void
ncm_diff_clear_domain (NcmDiff *diff)
{
  ncm_diff_set_domain (diff, NULL, NULL);
}

/**
 * ncm_diff_set_domain_warnings:
 * @diff: a #NcmDiff
 * @domain_warnings: whether to warn on a fallback to a one-sided difference
 *
 * Sets #NcmDiff:domain-warnings.
 */
void
ncm_diff_set_domain_warnings (NcmDiff *diff, const gboolean domain_warnings)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  self->domain_warnings = domain_warnings;
}

/**
 * ncm_diff_get_domain_warnings:
 * @diff: a #NcmDiff
 *
 * Returns: whether a fallback to a one-sided difference at an edge of the domain is warned about.
 */
gboolean
ncm_diff_get_domain_warnings (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->domain_warnings;
}

/**
 * ncm_diff_get_spectral_window:
 * @diff: a #NcmDiff
 *
 * Gets the initial spectral window half-width, see ncm_diff_set_spectral_window().
 *
 * Returns: the initial spectral window half-width.
 */
gdouble
ncm_diff_get_spectral_window (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->spectral_window;
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

/*
 * Difference quotient (f1 - f2) scale, in place of f1, and its absolute
 * round-off error eps max (|f1|, |f2|) |scale|. The round-off does not depend
 * on the difference, so a quotient that rounds to exactly zero still carries
 * the error of the values it came from.
 */
static void
_ncm_diff_step_quotient (NcmVector *f1, const NcmVector *f2, const gdouble scale, NcmVector *roff)
{
  const guint len = ncm_vector_len (f1);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble f1_i = ncm_vector_get (f1, i);
    const gdouble f2_i = ncm_vector_get (f2, i);

    ncm_vector_set (f1,   i, (f1_i - f2_i) * scale);
    ncm_vector_set (roff, i, GSL_MAX (fabs (f1_i), fabs (f2_i)) * GSL_DBL_EPSILON * fabs (scale));
  }
}

static void
_ncm_diff_rf_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  _ncm_diff_step_quotient (df, f_v, 1.0 / h, roff);

  NCM_UNUSED (yh1_v);
  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rc_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  _ncm_diff_step_quotient (df, yh1_v, 0.5 / h, roff);

  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rc_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  const guint len     = ncm_vector_len (df);
  const gdouble scale = 2.0 / (h * h);
  guint i;

  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  /*
   * The average of f (x + h) and f (x - h) can cancel exactly (an odd f at
   * zero); the round-off is that of the three values entering the quotient.
   */
  for (i = 0; i < len; i++)
  {
    const gdouble fp_i  = ncm_vector_get (df, i);
    const gdouble fm_i  = ncm_vector_get (yh1_v, i);
    const gdouble f0_i  = ncm_vector_get (f_v, i);
    const gdouble avg_i = (fp_i + fm_i) * 0.5;
    const gdouble mag_i = GSL_MAX (GSL_MAX (fabs (fp_i), fabs (fm_i)), fabs (f0_i));

    ncm_vector_set (df,   i, (avg_i - f0_i) * scale);
    ncm_vector_set (roff, i, mag_i * GSL_DBL_EPSILON * scale);
  }

  NCM_UNUSED (yh2_v);
}

/*
 * Forward second difference (f (x + 2h) - 2 f (x + h) + f (x)) / h^2, the
 * one-sided scheme rc_d2 falls back to at an edge of the domain; a negative
 * h makes it backward.
 */
static void
_ncm_diff_rf_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df, NcmVector *roff)
{
  const guint len     = ncm_vector_len (df);
  const gdouble scale = 1.0 / (h * h);
  guint i;

  ncm_vector_addto (x_v, a, h);

  f (x_v, df, user_data);

  ncm_vector_addto (x_v, a, h);

  f (x_v, yh1_v, user_data);

  for (i = 0; i < len; i++)
  {
    const gdouble f1_i  = ncm_vector_get (df, i);
    const gdouble f2_i  = ncm_vector_get (yh1_v, i);
    const gdouble f0_i  = ncm_vector_get (f_v, i);
    const gdouble mag_i = GSL_MAX (GSL_MAX (fabs (f1_i), fabs (f2_i)), fabs (f0_i));

    ncm_vector_set (df,   i, ((f2_i - f1_i) - (f1_i - f0_i)) * scale);
    ncm_vector_set (roff, i, 4.0 * mag_i * GSL_DBL_EPSILON * scale);
  }

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
    /* Absolute round-off error of df[0], as in the other step algorithms. */
    const gdouble max_s12 = GSL_MAX (fabs (fval + f_hxhy), fabs (f_hx + f_hy));

    roff[0] = max_s12 * GSL_DBL_EPSILON / fabs (hx * hy);
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
  gboolean informative;
} NcmDiffConv;

static void
_ncm_diff_conv_init (NcmDiffConv *cs)
{
  cs->df_best      = 0.0;
  cs->err_best     = GSL_POSINF;
  cs->err_last_max = 0.0;
  cs->not_conv     = NCM_DIFF_NTRY_CONV;
  cs->cstarted     = 0;
  cs->informative  = FALSE;
}

/*
 * The agreement of two orders is measured against scale, the largest
 * difference quotient of the row: for a nonzero derivative that is the
 * derivative itself, for a zero one the size of the terms being cancelled.
 */
static gboolean
_ncm_diff_conv_update (NcmDiffConv *cs, const gdouble terr_pad, const gdouble roff_pad,
                       const gdouble df_curr, const gdouble df_last,
                       const gdouble roff_curr, const gdouble roff_last, const gdouble scale)
{
  const gdouble err_trunc    = fabs (df_curr - df_last) * terr_pad;
  const gdouble err_err      = (scale > 0.0) ? fabs (df_curr - df_last) / scale : 0.0;
  const gdouble Eroff_last   = fabs (roff_last) * roff_pad;
  const gdouble Eroff_curr   = fabs (roff_curr) * roff_pad;
  const gdouble err_curr_max = GSL_MAX (err_trunc, GSL_MAX (Eroff_last, Eroff_curr));
  gdouble err_curr_best      = cs->err_best;
  gboolean improve           = FALSE;
  gboolean informative, agree;

  /*
   * A row is informative while its round-off is below 1.0e-3 of the scale;
   * rows above it (steps far below the scale f carries, quotients rounded to
   * zero) can only agree by accident. Two orders agree when they differ by
   * less than 1.0e-3 of the scale on an informative row. When every quotient
   * is zero further steps can only add round-off, so the rows agree and the
   * ladder stops with zero and the round-off as its error; such a row is
   * informative only when f is exactly zero on the samples, with no
   * round-off at all.
   */
  informative     = (scale > 0.0) ? (Eroff_curr < 1.0e-3 * scale) : (Eroff_curr == 0.0);
  agree           = (scale > 0.0) ? ((err_err < 1.0e-3) && informative) : TRUE;
  cs->informative = informative;

  /*
   * Estimates fluctuate in the beginning. Thus we only start checking for
   * convergence after they agree NCM_DIFF_NTRY_CONV times in a row.
   */
  if (cs->not_conv && agree)
  {
    cs->not_conv--;

    if (!cs->not_conv)
      cs->cstarted = 1;
  }
  else if (!agree)
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
  else if (err_curr_max < cs->err_best)
  {
    /* Not converged: keep the row with the smallest total error. */
    cs->df_best  = df_curr;
    cs->err_best = err_curr_max;

    err_curr_best = err_curr_max;
  }

  /*
   * The round-off of the rows grows with the order and bounds their total
   * error from below, so once it passes the best error no later row can rank
   * better. A ladder that has not converged goes on while its rows are
   * informative, since they may yet agree and replace a best value taken from
   * a false plateau (steps aliased with an oscillation of f).
   */
  if ((Eroff_curr < err_curr_best) || (cs->not_conv && informative))
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

/*
 * One Richardson ladder for one scalar component: the difference quotient
 * and round-off estimate of each step h0 rs^-t, the extrapolated rows and the
 * convergence control. The vector drivers keep one ladder per component and
 * append every evaluated step to all of them; a converged ladder is no longer
 * updated.
 */
typedef struct _NcmDiffLadder
{
  GPtrArray *tables;
  guint order;
  GArray *dfs;
  GArray *roffs;
  gdouble df_curr;
  gdouble df_last;
  gdouble roff_curr;
  gdouble roff_last;
  gdouble dmax;
  gdouble max_term;
  gdouble rho;
  NcmDiffConv conv;
  gboolean converged;
  gboolean first_informative;
} NcmDiffLadder;

static void
_ncm_diff_ladder_init (NcmDiffLadder *ladder, GPtrArray *tables)
{
  ladder->tables = tables;
  ladder->dfs    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ladder->roffs  = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void
_ncm_diff_ladder_clear (NcmDiffLadder *ladder)
{
  g_clear_pointer (&ladder->dfs, g_array_unref);
  g_clear_pointer (&ladder->roffs, g_array_unref);
}

/* Restarts the ladder with no steps. */
static void
_ncm_diff_ladder_reset (NcmDiffLadder *ladder)
{
  ladder->order     = 0;
  ladder->df_curr   = 0.0;
  ladder->df_last   = 0.0;
  ladder->roff_curr = 0.0;
  ladder->roff_last = 0.0;
  ladder->dmax      = 0.0;
  ladder->max_term  = 0.0;
  ladder->rho       = 0.0;
  ladder->converged = FALSE;

  ladder->first_informative = FALSE;

  g_array_set_size (ladder->dfs, 0);
  g_array_set_size (ladder->roffs, 0);

  _ncm_diff_conv_init (&ladder->conv);
}

/* Restarts the rows and the control, keeping the steps. */
static void
_ncm_diff_ladder_restart (NcmDiffLadder *ladder)
{
  ladder->order     = 0;
  ladder->df_curr   = 0.0;
  ladder->df_last   = 0.0;
  ladder->roff_curr = 0.0;
  ladder->roff_last = 0.0;
  ladder->dmax      = 0.0;
  ladder->max_term  = 0.0;
  ladder->rho       = 0.0;
  ladder->converged = FALSE;

  ladder->first_informative = FALSE;

  _ncm_diff_conv_init (&ladder->conv);
}

/* Extrapolation order k uses the first k + 2 steps. */
static gboolean
_ncm_diff_ladder_needs_step (NcmDiffLadder *ladder)
{
  return !ladder->converged && (ladder->dfs->len < ladder->order + 2);
}

static void
_ncm_diff_ladder_add_step (NcmDiffLadder *ladder, const gdouble df, const gdouble roff)
{
  g_array_append_val (ladder->dfs, df);
  g_array_append_val (ladder->roffs, roff);
}

/* A step larger than every stored one becomes the new top of the ladder. */
static void
_ncm_diff_ladder_prepend_step (NcmDiffLadder *ladder, const gdouble df, const gdouble roff)
{
  g_array_prepend_val (ladder->dfs, df);
  g_array_prepend_val (ladder->roffs, roff);
}

/*
 * Richardson extrapolation at the order of dtable: df = sum_t lambda_t dfs[t],
 * with the round-off estimates combined in quadrature. Also measures the row:
 * the largest quotient dmax = max_t |dfs[t]|, the scale of what is being
 * extrapolated; the largest term max_term = max_t |lambda_t dfs[t]|; and
 * rho = |df| / max_term, the fraction of the terms surviving the sum, near
 * 1 / max |lambda| for a nonzero derivative and far below it when the terms
 * cancel (a zero derivative), zero when every quotient is zero.
 */
static void
_ncm_diff_ladder_accum (NcmDiffLadder *ladder, NcmDiffTable *dtable, const guint nt)
{
  gdouble dmax     = 0.0;
  gdouble max_term = 0.0;
  guint t;

  ladder->df_curr   = 0.0;
  ladder->roff_curr = 0.0;

  for (t = 0; t < nt; t++)
  {
    const gdouble lambda_t = ncm_vector_get (dtable->lambda, t);
    const gdouble df_t     = g_array_index (ladder->dfs, gdouble, t);
    const gdouble roff_t   = g_array_index (ladder->roffs, gdouble, t);

    /* Fused, as the BLAS daxpy of the former vector accumulation. */
    ladder->df_curr   = fma (lambda_t, df_t, ladder->df_curr);
    ladder->roff_curr = hypot (ladder->roff_curr, lambda_t * roff_t);

    dmax     = GSL_MAX (dmax, fabs (df_t));
    max_term = GSL_MAX (max_term, fabs (lambda_t * df_t));
  }

  ladder->dmax     = dmax;
  ladder->max_term = max_term;
  ladder->rho      = (max_term > 0.0) ? fabs (ladder->df_curr) / max_term : 0.0;
}

/*
 * Extrapolates at the current order, updates the control unless the ladder
 * has converged, and moves to the next order. The rows of a converged ladder
 * are still computed for the controls that read them.
 */
static void
_ncm_diff_ladder_extrapolate (NcmDiffLadder *ladder, const gdouble terr_pad, const gdouble roff_pad)
{
  NcmDiffTable *dtable;

  /* No table for this order, or not enough steps for its row (a ladder that
   * converged at the last stored step of a replay). */
  if ((ladder->order == ladder->tables->len) || (ladder->order + 2 > ladder->dfs->len))
    return;

  dtable = g_ptr_array_index (ladder->tables, ladder->order);

  if (ladder->order == 0)
  {
    ladder->df_last      = g_array_index (ladder->dfs, gdouble, 0);
    ladder->roff_last    = g_array_index (ladder->roffs, gdouble, 0);
    ladder->conv.df_best = ladder->df_last;
  }

  _ncm_diff_ladder_accum (ladder, dtable, ladder->order + 2);

  if (!ladder->converged && !_ncm_diff_conv_update (&ladder->conv, terr_pad, roff_pad,
                                                    ladder->df_curr, ladder->df_last,
                                                    ladder->roff_curr, ladder->roff_last, ladder->dmax))
    ladder->converged = TRUE;

  if (ladder->order == 0)
    ladder->first_informative = ladder->conv.informative;

  ladder->df_last   = ladder->df_curr;
  ladder->roff_last = ladder->roff_curr;
  ladder->order++;

  if (ladder->order == ladder->tables->len)
    ladder->converged = TRUE;
}

/*
 * Replays the rows and the control over the stored steps after the top has
 * changed, up to convergence or the last stored step; the same as a fresh
 * run from the new top, without evaluating f again.
 */
static void
_ncm_diff_ladder_replay (NcmDiffLadder *ladder, const gdouble terr_pad, const gdouble roff_pad)
{
  _ncm_diff_ladder_restart (ladder);

  while (!ladder->converged && (ladder->order + 2 <= ladder->dfs->len))
    _ncm_diff_ladder_extrapolate (ladder, terr_pad, roff_pad);
}

/*
 * Step of index k, h0 rs^-k, through the formula of the tables so that
 * k >= 0 reproduces their steps; k < 0 are the steps above the initial one.
 */
static gdouble
_ncm_diff_step_h (NcmDiffPrivate *self, const guint po, const gdouble h0, const gint k)
{
  if (po == 0)
    return h0 * (1.0 / pow (self->rs, k));
  else
    return h0 * sqrt (1.0 / pow (self->rs, 2 * k));
}

/* Evaluates the step ho along coordinate a into df_t and roff_t, restoring x_v. */
static void
_ncm_diff_eval_step (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffFuncNtoM f, gpointer user_data,
                     const guint a, const gdouble x, const gdouble ho,
                     NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *df_t, NcmVector *roff_t)
{
  volatile gdouble temp = x + ho;
  const gdouble h       = temp - x;

  step_algo (diff, f, user_data, a, x, h, x_v, f_v, yh1_v, yh2_v, df_t, roff_t);
  ncm_vector_set (x_v, a, x);
}

static gboolean
_ncm_diff_vector_finite (NcmVector *v)
{
  const guint len = ncm_vector_len (v);
  guint i;

  for (i = 0; i < len; i++)
    if (!gsl_finite (ncm_vector_get (v, i)))
      return FALSE;

  return TRUE;
}

#define NCM_DIFF_MOVE_UP_GUARD (3)

/*
 * Largest step a move up may take along a coordinate: the guard above the
 * step a zero coordinate gets. Zero is not treated as an edge of the domain:
 * no point is unless the caller sets a domain, the initial step already
 * crosses an edge that lies within ini_h |x| of x, and near one the caller
 * chooses the forward scheme or a smaller ini_h.
 */
static gdouble
_ncm_diff_move_up_cap (NcmDiffPrivate *self, const gdouble x)
{
  const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);

  return self->ini_h * GSL_MAX (1.0, scale) * pow (self->rs, NCM_DIFF_MOVE_UP_GUARD);
}

/* Distances from x to the edges of the domain of coordinate a, infinite without a domain. */
static void
_ncm_diff_room (NcmDiffPrivate *self, const guint a, const gdouble x, gdouble *lo, gdouble *hi)
{
  lo[0] = GSL_POSINF;
  hi[0] = GSL_POSINF;

  if ((self->lb != NULL) && (a < ncm_vector_len (self->lb)))
    lo[0] = x - ncm_vector_get (self->lb, a);

  if ((self->ub != NULL) && (a < ncm_vector_len (self->ub)))
    hi[0] = ncm_vector_get (self->ub, a) - x;

  if ((lo[0] < 0.0) || (hi[0] < 0.0))
    g_error ("NcmDiff: coordinate %u at % .15g is outside its domain [% .15g, % .15g].", a, x,
             (self->lb != NULL) ? ncm_vector_get (self->lb, a) : GSL_NEGINF,
             (self->ub != NULL) ? ncm_vector_get (self->ub, a) : GSL_POSINF);
}

/*
 * The scheme in use along one coordinate: the step function and tables, the
 * sign of the steps (a backward difference is a forward one with negative
 * steps), the magnitude of the initial step and the room the domain leaves
 * on each side. A central scheme that cannot keep its points inside the
 * domain falls back to the one-sided scheme toward the farther edge.
 */
typedef struct _NcmDiffScheme
{
  NcmDiffStepAlgo algo;
  NcmDiffStepAlgo onesided;
  GPtrArray *tables;
  guint po;
  gboolean central;
  gdouble sign;
  gdouble h0;
  gdouble lo;
  gdouble hi;
} NcmDiffScheme;

/* Room the scheme has: the nearer edge for a central one, the edge ahead for a one-sided one. */
static gdouble
_ncm_diff_scheme_room (const NcmDiffScheme *sch)
{
  if (sch->central)
    return GSL_MIN (sch->lo, sch->hi);

  return (sch->sign > 0.0) ? sch->hi : sch->lo;
}

/* Largest step the scheme may take: the move-up guard, and a step that keeps the points strictly inside the domain. */
static gdouble
_ncm_diff_scheme_cap (NcmDiffPrivate *self, const NcmDiffScheme *sch, const gdouble x)
{
  return GSL_MIN (_ncm_diff_move_up_cap (self, x), _ncm_diff_scheme_room (sch) * (1.0 - 1.0 / self->rs));
}

static void
_ncm_diff_scheme_warn (NcmDiffPrivate *self, const NcmDiffScheme *sch, const guint a, const gdouble x, const gchar *from)
{
  if (self->domain_warnings)
    g_warning ("NcmDiff: coordinate %u at % .15g is %.3e from the edge of its domain, "
               "falling back from a %s difference to a %s one.",
               a, x, _ncm_diff_scheme_room (sch), from, (sch->sign > 0.0) ? "forward" : "backward");
}

/*
 * Switches to the one-sided scheme toward the farther edge, or flips the
 * side of a one-sided one, when that gives more room than the present
 * scheme has. Returns whether the scheme changed.
 */
static gboolean
_ncm_diff_scheme_fallback (NcmDiffPrivate *self, NcmDiffScheme *sch, const guint a, const gdouble x)
{
  const gdouble room = _ncm_diff_scheme_room (sch);
  const gdouble sign = (sch->hi >= sch->lo) ? 1.0 : -1.0;
  const gdouble far  = GSL_MAX (sch->lo, sch->hi);

  if (far <= room)
    return FALSE;

  if (sch->central)
  {
    sch->central = FALSE;
    sch->algo    = sch->onesided;
    sch->tables  = self->forward_tables;
    sch->po      = 0;
    sch->sign    = sign;
    _ncm_diff_scheme_warn (self, sch, a, x, "central");
  }
  else
  {
    sch->sign = sign;
    _ncm_diff_scheme_warn (self, sch, a, x, (sign > 0.0) ? "backward" : "forward");
  }

  return TRUE;
}

/*
 * Starts coordinate a with the scheme the caller asked for and the initial
 * step ini_h |x|, within the room of the domain; without room for it on the
 * asked side the scheme falls back at once.
 */
static void
_ncm_diff_scheme_init (NcmDiffPrivate *self, NcmDiffScheme *sch, NcmDiffStepAlgo algo, NcmDiffStepAlgo onesided, const guint po,
                       const guint a, const gdouble x)
{
  const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);

  sch->algo     = algo;
  sch->onesided = onesided;
  sch->tables   = (po == 0) ? self->forward_tables : self->central_tables;
  sch->po       = po;
  sch->central  = (po != 0);
  sch->sign     = 1.0;
  sch->h0       = self->ini_h * scale;

  _ncm_diff_room (self, a, x, &sch->lo, &sch->hi);

  if (sch->h0 > _ncm_diff_scheme_cap (self, sch, x))
  {
    _ncm_diff_scheme_fallback (self, sch, a, x);
    sch->h0 = GSL_MIN (sch->h0, _ncm_diff_scheme_cap (self, sch, x));
  }
}

static GArray *
_ncm_diff_by_step_algo_single (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GArray *ladders             = g_array_new (FALSE, FALSE, sizeof (NcmDiffLadder));
  NcmVector *x_v              = NULL;
  NcmVector *f_v              = NULL;
  NcmVector *yh1_v            = NULL;
  NcmVector *yh2_v            = NULL;
  NcmVector *df_t             = NULL;
  NcmVector *roff_t           = NULL;
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar            = x_a->len;
  NcmMatrix *Eerr_m           = NULL;
  NcmMatrix *df_m;
  guint a, i;


  g_array_set_size (df, dim * nvar);
  df_m = ncm_matrix_new_array (df, dim);

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, dim * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, dim);
  }

  g_array_set_size (ladders, dim);

  for (i = 0; i < dim; i++)
    _ncm_diff_ladder_init (&g_array_index (ladders, NcmDiffLadder, i), self->forward_tables);

  x_v    = ncm_vector_new_array (x_a);
  f_v    = ncm_vector_new (dim);
  yh1_v  = ncm_vector_new (dim);
  yh2_v  = ncm_vector_new (dim);
  df_t   = ncm_vector_new (dim);
  roff_t = ncm_vector_new (dim);

  f (x_v, f_v, user_data);

  for (a = 0; a < nvar; a++)
  {
    const gdouble x     = g_array_index (x_a, gdouble, a);
    const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);
    gint k_top          = 0;
    gboolean may_move   = TRUE;
    gboolean running    = TRUE;
    NcmDiffScheme sch;


    _ncm_diff_scheme_init (self, &sch, step_algo, onesided_algo, po, a, x);

    for (i = 0; i < dim; i++)
    {
      NcmDiffLadder *ladder = &g_array_index (ladders, NcmDiffLadder, i);

      ladder->tables = sch.tables;
      _ncm_diff_ladder_reset (ladder);
    }

    while (running)
    {
      gboolean needs_step    = FALSE;
      gboolean uninformative = FALSE;

      for (i = 0; i < dim; i++)
        needs_step = needs_step || _ncm_diff_ladder_needs_step (&g_array_index (ladders, NcmDiffLadder, i));

      if (needs_step)
      {
        const guint nsteps = g_array_index (ladders, NcmDiffLadder, 0).dfs->len;

        _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * _ncm_diff_step_h (self, sch.po, sch.h0, k_top + nsteps),
                             x_v, f_v, yh1_v, yh2_v, df_t, roff_t);

        for (i = 0; i < dim; i++)
          _ncm_diff_ladder_add_step (&g_array_index (ladders, NcmDiffLadder, i), ncm_vector_get (df_t, i), ncm_vector_get (roff_t, i));

        continue;
      }

      for (i = 0; i < dim; i++)
      {
        NcmDiffLadder *ladder = &g_array_index (ladders, NcmDiffLadder, i);

        if (!ladder->converged)
          _ncm_diff_ladder_extrapolate (ladder, self->terr_pad, self->roff_pad);
      }

      /*
       * A first row that is not informative means the initial step is far
       * below the scale of f and the quotients are cancellation noise; the
       * only cure is a larger step. The first move jumps to the step a zero
       * coordinate gets and restarts the ladders there, since the steps
       * below are no longer consecutive; further moves go up by the
       * Richardson factor, within the guard, and the controls replay over
       * the steps from the new top. A step that leaves the domain of f (a
       * non-finite value) is rejected. When the domain set by the caller
       * leaves no room for the move, the scheme falls back to a one-sided
       * one and restarts; a domain with no room on either side is an error.
       */
      for (i = 0; i < dim; i++)
        uninformative = uninformative || !g_array_index (ladders, NcmDiffLadder, i).first_informative;

      if (may_move && uninformative)
      {
        const gboolean jump = (k_top == 0) && (scale < 1.0);
        const gdouble cap   = _ncm_diff_scheme_cap (self, &sch, x);
        const gint k_new    = jump ? (gint) ceil (log (sch.h0 / GSL_MIN (self->ini_h, cap)) / log (self->rs)) : k_top - 1;
        const gdouble h_new = _ncm_diff_step_h (self, sch.po, sch.h0, k_new);

        if ((k_new >= k_top) || (h_new > cap))
        {
          if (_ncm_diff_scheme_fallback (self, &sch, a, x))
          {
            sch.h0 = GSL_MIN (self->ini_h * GSL_MAX (1.0, scale), _ncm_diff_scheme_cap (self, &sch, x));
            k_top  = 0;

            for (i = 0; i < dim; i++)
            {
              NcmDiffLadder *ladder = &g_array_index (ladders, NcmDiffLadder, i);

              ladder->tables = sch.tables;
              _ncm_diff_ladder_reset (ladder);
            }

            continue;
          }

          may_move = FALSE;

          if (_ncm_diff_scheme_room (&sch) < self->ini_h * scale)
            g_error ("NcmDiff: coordinate %u at % .15g: its domain leaves room for steps of at most %.3e, "
                     "below the initial step %.3e, and the quotients are round-off.",
                     a, x, _ncm_diff_scheme_room (&sch), self->ini_h * scale);
        }
        else
        {
          _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * h_new, x_v, f_v, yh1_v, yh2_v, df_t, roff_t);

          if (!_ncm_diff_vector_finite (df_t))
          {
            may_move = FALSE;
          }
          else
          {
            for (i = 0; i < dim; i++)
            {
              NcmDiffLadder *ladder = &g_array_index (ladders, NcmDiffLadder, i);

              if (jump)
              {
                _ncm_diff_ladder_reset (ladder);
                _ncm_diff_ladder_add_step (ladder, ncm_vector_get (df_t, i), ncm_vector_get (roff_t, i));
              }
              else
              {
                _ncm_diff_ladder_prepend_step (ladder, ncm_vector_get (df_t, i), ncm_vector_get (roff_t, i));
                _ncm_diff_ladder_replay (ladder, self->terr_pad, self->roff_pad);
              }
            }

            k_top = k_new;

            continue;
          }
        }
      }

      running = FALSE;

      for (i = 0; i < dim; i++)
        running = running || !g_array_index (ladders, NcmDiffLadder, i).converged;
    }

    for (i = 0; i < dim; i++)
    {
      const NcmDiffLadder *ladder = &g_array_index (ladders, NcmDiffLadder, i);

      ncm_matrix_set (df_m, a, i, ladder->conv.df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, i, ladder->conv.err_best);
    }
  }

  if (Eerr_m != NULL)
    ncm_matrix_scale (Eerr_m, NCM_DIFF_ERR_PAD);

  {
    for (i = 0; i < dim; i++)
      _ncm_diff_ladder_clear (&g_array_index (ladders, NcmDiffLadder, i));

    g_array_unref (ladders);

    ncm_vector_clear (&x_v);
    ncm_vector_clear (&f_v);
    ncm_vector_clear (&yh1_v);
    ncm_vector_clear (&yh2_v);
    ncm_vector_clear (&df_t);
    ncm_vector_clear (&roff_t);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

/*
 * Two ladders sharing the step factors, A from h0 and B from h0 / sqrt (rs),
 * each with its own control, and a cross control that decides from their
 * disagreement at equal order.
 */
typedef struct _NcmDiffDual
{
  NcmDiffLadder A;
  NcmDiffLadder B;
  NcmDiffDualConv cross;
  gboolean converged;
} NcmDiffDual;

static void
_ncm_diff_dual_init (NcmDiffDual *dual, GPtrArray *tables)
{
  _ncm_diff_ladder_init (&dual->A, tables);
  _ncm_diff_ladder_init (&dual->B, tables);
}

static void
_ncm_diff_dual_clear (NcmDiffDual *dual)
{
  _ncm_diff_ladder_clear (&dual->A);
  _ncm_diff_ladder_clear (&dual->B);
}

static void
_ncm_diff_dual_reset (NcmDiffDual *dual)
{
  _ncm_diff_ladder_reset (&dual->A);
  _ncm_diff_ladder_reset (&dual->B);
  _ncm_diff_dual_conv_init (&dual->cross);
  dual->converged = FALSE;
}

static gboolean
_ncm_diff_dual_needs_step (NcmDiffDual *dual)
{
  return !dual->converged && (dual->A.dfs->len < dual->A.order + 2);
}

static void _ncm_diff_dual_extrapolate (NcmDiffDual *dual, const gdouble terr_pad, const gdouble roff_pad);

/* Replays both ladders and the cross control over the stored steps. */
static void
_ncm_diff_dual_replay (NcmDiffDual *dual, const gdouble terr_pad, const gdouble roff_pad)
{
  _ncm_diff_ladder_restart (&dual->A);
  _ncm_diff_ladder_restart (&dual->B);
  _ncm_diff_dual_conv_init (&dual->cross);
  dual->converged = FALSE;

  while (!dual->converged && (dual->A.order + 2 <= dual->A.dfs->len))
    _ncm_diff_dual_extrapolate (dual, terr_pad, roff_pad);
}

/* The first row of either ladder is not informative. */
static gboolean
_ncm_diff_dual_first_uninformative (NcmDiffDual *dual)
{
  return !dual->A.first_informative || !dual->B.first_informative;
}

/*
 * Extrapolates both ladders at the current order and lets the cross control
 * decide; the first NCM_DIFF_DUAL_MIN_ORDER orders are always taken.
 */
static void
_ncm_diff_dual_extrapolate (NcmDiffDual *dual, const gdouble terr_pad, const gdouble roff_pad)
{
  gboolean improve;

  _ncm_diff_ladder_extrapolate (&dual->A, terr_pad, roff_pad);
  _ncm_diff_ladder_extrapolate (&dual->B, terr_pad, roff_pad);

  improve = _ncm_diff_dual_conv_update (&dual->cross, dual->A.df_curr, dual->B.df_curr,
                                        dual->A.roff_curr, dual->B.roff_curr);

  if (((dual->A.order >= NCM_DIFF_DUAL_MIN_ORDER) && !improve) || (dual->A.order == dual->A.tables->len))
    dual->converged = TRUE;
}

static GArray *
_ncm_diff_by_step_algo_dual (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GArray *duals               = g_array_new (FALSE, FALSE, sizeof (NcmDiffDual));
  NcmVector *x_v              = NULL;
  NcmVector *f_v              = NULL;
  NcmVector *yh1_v            = NULL;
  NcmVector *yh2_v            = NULL;
  NcmVector *df_t[2]          = {NULL, NULL};
  NcmVector *roff_t[2]        = {NULL, NULL};
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar            = x_a->len;
  NcmMatrix *Eerr_m           = NULL;
  NcmMatrix *df_m;
  guint a, i, s;


  g_array_set_size (df, dim * nvar);
  df_m = ncm_matrix_new_array (df, dim);

  if (Eerr != NULL)
  {
    *Eerr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    g_array_set_size (*Eerr, dim * nvar);
    Eerr_m = ncm_matrix_new_array (*Eerr, dim);
  }

  g_array_set_size (duals, dim);

  for (i = 0; i < dim; i++)
    _ncm_diff_dual_init (&g_array_index (duals, NcmDiffDual, i), self->forward_tables);

  x_v   = ncm_vector_new_array (x_a);
  f_v   = ncm_vector_new (dim);
  yh1_v = ncm_vector_new (dim);
  yh2_v = ncm_vector_new (dim);

  for (s = 0; s < 2; s++)
  {
    df_t[s]   = ncm_vector_new (dim);
    roff_t[s] = ncm_vector_new (dim);
  }

  f (x_v, f_v, user_data);

  for (a = 0; a < nvar; a++)
  {
    const gdouble x     = g_array_index (x_a, gdouble, a);
    const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);
    gint k_top          = 0;
    gboolean may_move   = TRUE;
    gboolean may_jump   = (scale < 1.0);
    gboolean running    = TRUE;
    NcmDiffScheme sch;
    gdouble h0[2];


    _ncm_diff_scheme_init (self, &sch, step_algo, onesided_algo, po, a, x);
    h0[0] = sch.h0;
    h0[1] = sch.h0 / sqrt (self->rs);

    for (i = 0; i < dim; i++)
    {
      NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

      dual->A.tables = sch.tables;
      dual->B.tables = sch.tables;
      _ncm_diff_dual_reset (dual);
    }

    while (running)
    {
      gboolean needs_step    = FALSE;
      gboolean uninformative = FALSE;

      for (i = 0; i < dim; i++)
        needs_step = needs_step || _ncm_diff_dual_needs_step (&g_array_index (duals, NcmDiffDual, i));

      if (needs_step)
      {
        const gint k = k_top + g_array_index (duals, NcmDiffDual, 0).A.dfs->len;

        for (s = 0; s < 2; s++)
        {
          _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * _ncm_diff_step_h (self, sch.po, h0[s], k),
                               x_v, f_v, yh1_v, yh2_v, df_t[s], roff_t[s]);

          for (i = 0; i < dim; i++)
          {
            NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

            _ncm_diff_ladder_add_step ((s == 0) ? &dual->A : &dual->B, ncm_vector_get (df_t[s], i), ncm_vector_get (roff_t[s], i));
          }
        }

        continue;
      }

      for (i = 0; i < dim; i++)
      {
        NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

        if (!dual->converged)
          _ncm_diff_dual_extrapolate (dual, self->terr_pad, self->roff_pad);
      }

      /* Moving up, as in the single driver, both ladders together. */
      for (i = 0; i < dim; i++)
        uninformative = uninformative || _ncm_diff_dual_first_uninformative (&g_array_index (duals, NcmDiffDual, i));

      if (may_move && uninformative)
      {
        const gdouble cap       = _ncm_diff_scheme_cap (self, &sch, x);
        const gdouble base      = GSL_MIN (self->ini_h, cap);
        const gdouble h0_new[2] = {
          may_jump ? base : h0[0],
          may_jump ? base / sqrt (self->rs) : h0[1]
        };
        const gint k_new = may_jump ? 0 : k_top - 1;

        if ((may_jump && (base <= h0[0])) || (_ncm_diff_step_h (self, sch.po, h0_new[0], k_new) > cap))
        {
          if (_ncm_diff_scheme_fallback (self, &sch, a, x))
          {
            sch.h0 = GSL_MIN (self->ini_h * GSL_MAX (1.0, scale), _ncm_diff_scheme_cap (self, &sch, x));
            h0[0]  = sch.h0;
            h0[1]  = sch.h0 / sqrt (self->rs);
            k_top  = 0;

            for (i = 0; i < dim; i++)
            {
              NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

              dual->A.tables = sch.tables;
              dual->B.tables = sch.tables;
              _ncm_diff_dual_reset (dual);
            }

            continue;
          }

          may_move = FALSE;

          if (_ncm_diff_scheme_room (&sch) < self->ini_h * scale)
            g_error ("NcmDiff: coordinate %u at % .15g: its domain leaves room for steps of at most %.3e, "
                     "below the initial step %.3e, and the quotients are round-off.",
                     a, x, _ncm_diff_scheme_room (&sch), self->ini_h * scale);
        }
        else
        {
          gboolean finite = TRUE;

          for (s = 0; s < 2; s++)
          {
            _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * _ncm_diff_step_h (self, sch.po, h0_new[s], k_new),
                                 x_v, f_v, yh1_v, yh2_v, df_t[s], roff_t[s]);
            finite = finite && _ncm_diff_vector_finite (df_t[s]);
          }

          if (!finite)
          {
            may_move = FALSE;
          }
          else
          {
            for (i = 0; i < dim; i++)
            {
              NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

              if (may_jump)
                _ncm_diff_dual_reset (dual);

              for (s = 0; s < 2; s++)
              {
                NcmDiffLadder *ladder = (s == 0) ? &dual->A : &dual->B;

                if (may_jump)
                  _ncm_diff_ladder_add_step (ladder, ncm_vector_get (df_t[s], i), ncm_vector_get (roff_t[s], i));
                else
                  _ncm_diff_ladder_prepend_step (ladder, ncm_vector_get (df_t[s], i), ncm_vector_get (roff_t[s], i));
              }

              if (!may_jump)
                _ncm_diff_dual_replay (dual, self->terr_pad, self->roff_pad);
            }

            if (may_jump)
            {
              h0[0]    = h0_new[0];
              h0[1]    = h0_new[1];
              may_jump = FALSE;
            }
            else
            {
              k_top = k_new;
            }

            continue;
          }
        }
      }

      running = FALSE;

      for (i = 0; i < dim; i++)
        running = running || !g_array_index (duals, NcmDiffDual, i).converged;
    }

    for (i = 0; i < dim; i++)
    {
      const NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

      ncm_matrix_set (df_m, a, i, dual->cross.df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, i, dual->cross.err_best);
    }
  }

  {
    for (i = 0; i < dim; i++)
      _ncm_diff_dual_clear (&g_array_index (duals, NcmDiffDual, i));

    g_array_unref (duals);

    ncm_vector_clear (&x_v);
    ncm_vector_clear (&f_v);
    ncm_vector_clear (&yh1_v);
    ncm_vector_clear (&yh2_v);

    for (s = 0; s < 2; s++)
    {
      ncm_vector_clear (&df_t[s]);
      ncm_vector_clear (&roff_t[s]);
    }

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static GArray *
ncm_diff_by_step_algo (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  if (self->dual_series)
    return _ncm_diff_by_step_algo_dual (diff, step_algo, onesided_algo, po, x_a, dim, f, user_data, Eerr);

  return _ncm_diff_by_step_algo_single (diff, step_algo, onesided_algo, po, x_a, dim, f, user_data, Eerr);
}

/* Evaluates one Hessian cross-term step with offsets hxo and hyo: the
 * difference quotient df_t and its round-off estimate roff_t. */
static void
_ncm_diff_hessian_eval_step (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, NcmDiffFuncNto1 f, gpointer user_data,
                             const guint a, const gdouble x, const gdouble hxo,
                             const guint b, const gdouble y, const gdouble hyo,
                             NcmVector *x_v, const gdouble fval,
                             gdouble *df_t, gdouble *roff_t)
{
  volatile gdouble t_x = x + hxo;
  const gdouble hx     = t_x - x;
  volatile gdouble t_y = y + hyo;
  const gdouble hy     = t_y - y;

  Hstep_algo (diff, f, user_data, a, x, hx, b, y, hy, x_v, fval, df_t, roff_t);

  ncm_vector_set (x_v, a, x);
  ncm_vector_set (x_v, b, y);
}

/*
 * Moving up along the two axes of a cross term: the jump sets each axis
 * to the smaller of the step a zero coordinate gets and its cap, a further
 * move goes one factor up on both. Returns whether the move fits the caps;
 * when it does not, an axis without room falls back to the other side,
 * which restarts the pair, and an axis with room on neither side is an
 * error when the room is below its initial step.
 */
static gboolean
_ncm_diff_hessian_move_fits (NcmDiffPrivate *self, NcmDiffScheme *sx, NcmDiffScheme *sy, const guint a, const gdouble x, const guint b, const gdouble y,
                             const gboolean jump, const gint k_new, gdouble *hx0_new, gdouble *hy0_new, gboolean *restart)
{
  const gdouble scale_x = (x == 0.0) ? 1.0 : fabs (x);
  const gdouble scale_y = (y == 0.0) ? 1.0 : fabs (y);
  const gdouble cx      = _ncm_diff_scheme_cap (self, sx, x);
  const gdouble cy      = _ncm_diff_scheme_cap (self, sy, y);
  gboolean blocked_x, blocked_y;

  hx0_new[0] = jump ? GSL_MIN (self->ini_h * GSL_MAX (1.0, scale_x), cx) : sx->h0;
  hy0_new[0] = jump ? GSL_MIN (self->ini_h * GSL_MAX (1.0, scale_y), cy) : sy->h0;
  restart[0] = FALSE;

  if (jump)
  {
    blocked_x = (hx0_new[0] <= sx->h0);
    blocked_y = (hy0_new[0] <= sy->h0);

    if (!(blocked_x && blocked_y))
      return TRUE;
  }
  else
  {
    blocked_x = (_ncm_diff_step_h (self, 0, sx->h0, k_new) > cx);
    blocked_y = (_ncm_diff_step_h (self, 0, sy->h0, k_new) > cy);

    if (!blocked_x && !blocked_y)
      return TRUE;
  }

  if (blocked_x && _ncm_diff_scheme_fallback (self, sx, a, x))
    restart[0] = TRUE;

  if (blocked_y && _ncm_diff_scheme_fallback (self, sy, b, y))
    restart[0] = TRUE;

  if (!restart[0])
  {
    if (blocked_x && (_ncm_diff_scheme_room (sx) < self->ini_h * scale_x))
      g_error ("NcmDiff: coordinate %u at % .15g: its domain leaves room for steps of at most %.3e, "
               "below the initial step %.3e, and the quotients are round-off.",
               a, x, _ncm_diff_scheme_room (sx), self->ini_h * scale_x);

    if (blocked_y && (_ncm_diff_scheme_room (sy) < self->ini_h * scale_y))
      g_error ("NcmDiff: coordinate %u at % .15g: its domain leaves room for steps of at most %.3e, "
               "below the initial step %.3e, and the quotients are round-off.",
               b, y, _ncm_diff_scheme_room (sy), self->ini_h * scale_y);
  }

  return FALSE;
}

/* Initial steps of the two axes of a cross term, within their caps, after a fallback. */
static void
_ncm_diff_hessian_restart_bases (NcmDiffPrivate *self, NcmDiffScheme *sx, NcmDiffScheme *sy, const gdouble x, const gdouble y)
{
  const gdouble scale_x = (x == 0.0) ? 1.0 : fabs (x);
  const gdouble scale_y = (y == 0.0) ? 1.0 : fabs (y);

  sx->h0 = GSL_MIN (self->ini_h * GSL_MAX (1.0, scale_x), _ncm_diff_scheme_cap (self, sx, x));
  sy->h0 = GSL_MIN (self->ini_h * GSL_MAX (1.0, scale_y), _ncm_diff_scheme_cap (self, sy, y));
}

static GArray *
_ncm_diff_Hessian_by_step_algo_single (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GPtrArray *tables           = (po == 0) ? self->forward_tables : self->central_tables;
  NcmVector *x_v              = NULL;
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar            = x_a->len;
  gdouble fval                = 0.0;
  NcmMatrix *Eerr_m           = NULL;
  NcmDiffLadder ladder;
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

  _ncm_diff_ladder_init (&ladder, tables);

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
      gint k_top            = 0;
      gboolean may_move     = TRUE;
      gboolean may_jump     = (scale_x < 1.0) || (scale_y < 1.0);
      NcmDiffScheme sx, sy;


      _ncm_diff_scheme_init (self, &sx, NULL, NULL, 0, a, x);
      _ncm_diff_scheme_init (self, &sy, NULL, NULL, 0, b, y);
      _ncm_diff_ladder_reset (&ladder);

      while (!ladder.converged)
      {
        while (_ncm_diff_ladder_needs_step (&ladder))
        {
          const gint k = k_top + ladder.dfs->len;
          gdouble df_t, roff_t;

          _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                       a, x, sx.sign * _ncm_diff_step_h (self, 0, sx.h0, k),
                                       b, y, sy.sign * _ncm_diff_step_h (self, 0, sy.h0, k),
                                       x_v, fval, &df_t, &roff_t);
          _ncm_diff_ladder_add_step (&ladder, df_t, roff_t);
        }

        _ncm_diff_ladder_extrapolate (&ladder, self->terr_pad, self->roff_pad);

        /* Moving up, as in the vector driver: the two axes share the factor. */
        if (may_move && !ladder.first_informative)
        {
          const gint k_new = may_jump ? 0 : k_top - 1;
          gdouble hx0_new, hy0_new;
          gboolean restart;

          if (!_ncm_diff_hessian_move_fits (self, &sx, &sy, a, x, b, y, may_jump, k_new, &hx0_new, &hy0_new, &restart))
          {
            if (restart)
            {
              _ncm_diff_hessian_restart_bases (self, &sx, &sy, x, y);
              _ncm_diff_ladder_reset (&ladder);
              k_top    = 0;
              may_jump = FALSE;
              continue;
            }

            may_move = FALSE;
          }
          else
          {
            gdouble df_t, roff_t;

            _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                         a, x, sx.sign * _ncm_diff_step_h (self, 0, hx0_new, k_new),
                                         b, y, sy.sign * _ncm_diff_step_h (self, 0, hy0_new, k_new),
                                         x_v, fval, &df_t, &roff_t);

            if (!gsl_finite (df_t))
            {
              may_move = FALSE;
            }
            else if (may_jump)
            {
              _ncm_diff_ladder_reset (&ladder);
              _ncm_diff_ladder_add_step (&ladder, df_t, roff_t);
              sx.h0    = hx0_new;
              sy.h0    = hy0_new;
              may_jump = FALSE;
            }
            else
            {
              _ncm_diff_ladder_prepend_step (&ladder, df_t, roff_t);
              _ncm_diff_ladder_replay (&ladder, self->terr_pad, self->roff_pad);
              k_top = k_new;
            }
          }
        }
      }

      ncm_matrix_set (df_m, a, b, ladder.conv.df_best);
      ncm_matrix_set (df_m, b, a, ladder.conv.df_best);

      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, ladder.conv.err_best);
        ncm_matrix_set (Eerr_m, b, a, ladder.conv.err_best);
      }
    }
  }

  if (Eerr_m != NULL)
    ncm_matrix_scale (Eerr_m, NCM_DIFF_ERR_PAD);

  {
    _ncm_diff_ladder_clear (&ladder);

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
  GPtrArray *tables           = (po == 0) ? self->forward_tables : self->central_tables;
  NcmVector *x_v              = NULL;
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar            = x_a->len;
  gdouble fval                = 0.0;
  NcmMatrix *Eerr_m           = NULL;
  NcmDiffDual dual;
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

  _ncm_diff_dual_init (&dual, tables);

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
      const gdouble srs     = sqrt (self->rs);
      gint k_top            = 0;
      gboolean may_move     = TRUE;
      gboolean may_jump     = (scale_x < 1.0) || (scale_y < 1.0);
      NcmDiffScheme sx, sy;


      _ncm_diff_scheme_init (self, &sx, NULL, NULL, 0, a, x);
      _ncm_diff_scheme_init (self, &sy, NULL, NULL, 0, b, y);
      _ncm_diff_dual_reset (&dual);

      while (!dual.converged)
      {
        while (_ncm_diff_dual_needs_step (&dual))
        {
          const gint k = k_top + dual.A.dfs->len;
          guint s;

          for (s = 0; s < 2; s++)
          {
            const gdouble fac = (s == 0) ? 1.0 : 1.0 / srs;
            gdouble df_t, roff_t;

            _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                         a, x, sx.sign * _ncm_diff_step_h (self, 0, sx.h0 * fac, k),
                                         b, y, sy.sign * _ncm_diff_step_h (self, 0, sy.h0 * fac, k),
                                         x_v, fval, &df_t, &roff_t);
            _ncm_diff_ladder_add_step ((s == 0) ? &dual.A : &dual.B, df_t, roff_t);
          }
        }

        _ncm_diff_dual_extrapolate (&dual, self->terr_pad, self->roff_pad);

        /* Moving up, as in the single Hessian driver, both ladders together. */
        if (may_move && _ncm_diff_dual_first_uninformative (&dual))
        {
          const gint k_new = may_jump ? 0 : k_top - 1;
          gdouble hx0_new, hy0_new;
          gboolean restart;

          if (!_ncm_diff_hessian_move_fits (self, &sx, &sy, a, x, b, y, may_jump, k_new, &hx0_new, &hy0_new, &restart))
          {
            if (restart)
            {
              _ncm_diff_hessian_restart_bases (self, &sx, &sy, x, y);
              _ncm_diff_dual_reset (&dual);
              k_top    = 0;
              may_jump = FALSE;
              continue;
            }

            may_move = FALSE;
          }
          else
          {
            gdouble df_t[2], roff_t[2];
            gboolean finite = TRUE;
            guint s;

            for (s = 0; s < 2; s++)
            {
              const gdouble fac = (s == 0) ? 1.0 : 1.0 / srs;

              _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                           a, x, sx.sign * _ncm_diff_step_h (self, 0, hx0_new * fac, k_new),
                                           b, y, sy.sign * _ncm_diff_step_h (self, 0, hy0_new * fac, k_new),
                                           x_v, fval, &df_t[s], &roff_t[s]);
              finite = finite && gsl_finite (df_t[s]);
            }

            if (!finite)
            {
              may_move = FALSE;
            }
            else if (may_jump)
            {
              _ncm_diff_dual_reset (&dual);
              _ncm_diff_ladder_add_step (&dual.A, df_t[0], roff_t[0]);
              _ncm_diff_ladder_add_step (&dual.B, df_t[1], roff_t[1]);
              sx.h0    = hx0_new;
              sy.h0    = hy0_new;
              may_jump = FALSE;
            }
            else
            {
              _ncm_diff_ladder_prepend_step (&dual.A, df_t[0], roff_t[0]);
              _ncm_diff_ladder_prepend_step (&dual.B, df_t[1], roff_t[1]);
              _ncm_diff_dual_replay (&dual, self->terr_pad, self->roff_pad);
              k_top = k_new;
            }
          }
        }
      }

      ncm_matrix_set (df_m, a, b, dual.cross.df_best);
      ncm_matrix_set (df_m, b, a, dual.cross.df_best);

      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, dual.cross.err_best);
        ncm_matrix_set (Eerr_m, b, a, dual.cross.err_best);
      }
    }
  }

  {
    _ncm_diff_dual_clear (&dual);

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

/*
 * Spectral (Chebyshev) derivatives.
 *
 * The derivative along each variable is computed from a Chebyshev fit of the
 * function on a window [x - R, x + R]. The window half-width R is found by
 * probing a dyadic ladder of candidates with a fixed-order fit and scoring
 * each by the estimated derivative error (series tail plus propagated
 * round-off): shrinking resolves sharp features, expanding lowers the
 * round-off amplification of flat ones. The accepted window is then refined
 * on nested Chebyshev-Lobatto grids.
 *
 * Convergence is judged on the derivative itself between refinement levels,
 * never on the coefficients: a relative test on the coefficient norm would be
 * dominated by components the derivative does not use (e.g. a large constant
 * offset), declaring convergence while the derivative-carrying coefficients
 * are still unresolved. The round-off floor eps * sum |c_k|, propagated
 * through the coefficient differentiation, carries such offsets explicitly.
 */

#define NCM_DIFF_SC_PROBE_N (17) /* fixed order of the window probes */
#define NCM_DIFF_SC_MAX_N (65)   /* refinement grids: 17, 33, 65 nodes */
#define NCM_DIFF_SC_ERR_PAD (10.0)
#define NCM_DIFF_SC_MAX_SHRINK (45)
#define NCM_DIFF_SC_MAX_EXPAND (8)
#define NCM_DIFF_SC_BAD_Q (1.0e-3) /* q_fit above this marks an unresolved window */
#define NCM_DIFF_SC_PATIENCE (2)   /* non-improving steps allowed after a resolved window */
#define NCM_DIFF_SC_IMPROVE (0.9)  /* required score reduction factor */

typedef struct _NcmDiffSCData
{
  NcmDiffFuncNtoM f;
  gpointer user_data;
  NcmVector *x_v;
  guint a;
} NcmDiffSCData;

/*
 * Evaluates the function at the Chebyshev-Lobatto nodes of [x - R, x + R]
 * along variable data->a, storing component c at fvals[c][i]. With N_old == 0
 * all N nodes are evaluated; with N_new == 2 * N_old - 1 the old values are
 * spread to the even indices and only the new odd nodes are evaluated.
 * Returns FALSE when any value is non-finite (fvals is then partial).
 */
static gboolean
_ncm_diff_sc_eval_nodes (NcmDiffSCData *data, const gdouble x, const gdouble R,
                         const guint dim, NcmVector *y_v, NcmMatrix *fvals,
                         const guint N_old, const guint N_new)
{
  gboolean finite = TRUE;
  guint i, c;

  if (N_old > 0)
  {
    g_assert_cmpuint (N_new, ==, 2 * N_old - 1);

    for (i = N_old - 1; ; i--)
    {
      for (c = 0; c < dim; c++)
        ncm_matrix_set (fvals, c, 2 * i, ncm_matrix_get (fvals, c, i));

      if (i == 0)
        break;
    }
  }

  for (i = 0; i < N_new; i++)
  {
    const gboolean new_node = (N_old == 0) || ((i % 2) == 1);

    if (new_node)
    {
      const gdouble t  = cos (M_PI * i / (N_new - 1.0));
      const gdouble xi = x + R * t;

      ncm_vector_set (data->x_v, data->a, xi);
      data->f (data->x_v, y_v, data->user_data);

      for (c = 0; c < dim; c++)
      {
        const gdouble yc = ncm_vector_get (y_v, c);

        ncm_matrix_set (fvals, c, i, yc);

        if (!gsl_finite (yc))
          finite = FALSE;
      }
    }
  }

  ncm_vector_set (data->x_v, data->a, x);

  return finite;
}

/*
 * Direct DCT-I of the node values: coeffs[c][k] such that component c is
 * sum_k coeffs[c][k] T_k(t) on the window. O(N^2) per component, negligible
 * against the function evaluations for the N used here.
 */
static void
_ncm_diff_sc_dct (NcmMatrix *fvals, const guint dim, const guint N, NcmMatrix *coeffs)
{
  const guint two_Nm1 = 2 * (N - 1);
  gdouble *cosm       = g_new (gdouble, two_Nm1);
  guint i, k, c;

  for (i = 0; i < two_Nm1; i++)
    cosm[i] = cos ((M_PI * i) / (N - 1.0));

  for (c = 0; c < dim; c++)
  {
    const gdouble f_first = ncm_matrix_get (fvals, c, 0);
    const gdouble f_last  = ncm_matrix_get (fvals, c, N - 1);

    for (k = 0; k < N; k++)
    {
      gdouble s = 0.5 * (f_first + (((k % 2) == 0) ? f_last : -f_last));

      for (i = 1; i < N - 1; i++)
        s += ncm_matrix_get (fvals, c, i) * cosm[(k * i) % two_Nm1];

      ncm_matrix_set (coeffs, c, k, s * (((k == 0) || (k == N - 1)) ? 1.0 : 2.0) / (N - 1.0));
    }
  }

  g_free (cosm);
}

/*
 * Derivative of a Chebyshev series in coefficient space: given
 * f(t) = sum_{k=0}^{n-1} c_k T_k(t), fills b with the n - 1 coefficients of
 * f'(t) in the same convention. With non-negative input the output is
 * non-negative, so the same recurrence propagates error magnitudes.
 */
static void
_ncm_diff_sc_cheb_deriv (const gdouble *c, const guint n, gdouble *b)
{
  guint j;

  b[n - 2] = 2.0 * (n - 1.0) * c[n - 1];

  for (j = n - 2; j >= 1; j--)
    b[j - 1] = ((j + 1 <= n - 2) ? b[j + 1] : 0.0) + 2.0 * j * c[j];

  b[0] *= 0.5;
}

/* Chebyshev series value at the window center, t = 0. */
static gdouble
_ncm_diff_sc_eval0 (const gdouble *c, const guint n)
{
  gdouble s = 0.0;
  guint k;

  for (k = 0; k < n; k += 2)
    s += (((k % 4) == 0) ? c[k] : -c[k]);

  return s;
}

/*
 * Derivative of @order of one fitted component at the window center, with an
 * unpadded error estimate in two parts:
 *
 * - tail: the top quarter of the fitted coefficients, propagated in absolute
 *   value through the differentiation (truncation). Measured on the original
 *   coefficients, never on the differentiated ones: differentiation makes a
 *   slowly decaying series bottom-heavy, so its own tail can look converged
 *   while the fit has not resolved the function at all.
 * - roff: the coefficient round-off eps * sum |c_k|, propagated the same way
 *   (the amplification of the k-th coefficient grows as k^2 per order).
 *
 * Also returns q_fit, the fraction of the coefficient mass in the top
 * quarter: a scale-free measure of how resolved the fit is.
 */
static void
_ncm_diff_sc_deriv_est (const gdouble *c, const guint len, const gdouble R, const guint order,
                        gdouble *deriv, gdouble *tail, gdouble *roff, gdouble *q_fit)
{
  const gdouble Rinv     = 1.0 / R;
  const guint tail_start = (3 * len) / 4;
  gdouble *b             = g_new (gdouble, len);
  gdouble *nb            = g_new (gdouble, len);
  gdouble *tb            = g_new (gdouble, len);
  gdouble *tmp           = g_new (gdouble, len);
  gdouble sum_abs_c      = 0.0;
  gdouble sum_abs_tail   = 0.0;
  gdouble scale          = 1.0;
  guint n                = len;
  guint k, j;

  g_assert_cmpuint (len, >, order + 1);

  for (k = 0; k < len; k++)
  {
    const gdouble abs_ck = fabs (c[k]);

    b[k]       = c[k];
    tb[k]      = (k >= tail_start) ? abs_ck : 0.0;
    sum_abs_c += abs_ck;

    if (k >= tail_start)
      sum_abs_tail += abs_ck;
  }

  q_fit[0] = sum_abs_tail / (sum_abs_c + GSL_DBL_MIN);

  {
    const gdouble nu = GSL_DBL_EPSILON * sum_abs_c;

    for (k = 0; k < len; k++)
      nb[k] = nu;
  }

  for (j = 0; j < order; j++)
  {
    _ncm_diff_sc_cheb_deriv (b, n, tmp);
    memcpy (b, tmp, sizeof (gdouble) * (n - 1));

    _ncm_diff_sc_cheb_deriv (nb, n, tmp);
    memcpy (nb, tmp, sizeof (gdouble) * (n - 1));

    _ncm_diff_sc_cheb_deriv (tb, n, tmp);
    memcpy (tb, tmp, sizeof (gdouble) * (n - 1));

    n--;
    scale *= Rinv;
  }

  deriv[0] = _ncm_diff_sc_eval0 (b, n) * scale;

  {
    gdouble tail_sum = 0.0;
    gdouble roff_sum = 0.0;

    for (k = 0; k < n; k++)
    {
      tail_sum += tb[k];
      roff_sum += nb[k];
    }

    tail[0] = tail_sum * scale;
    roff[0] = roff_sum * scale;
  }

  g_free (b);
  g_free (nb);
  g_free (tb);
  g_free (tmp);
}

/*
 * Probes the window [x - R, x + R]: node evaluation plus fixed-order fit.
 * Produces two figures of merit over the components:
 *
 * - score: the estimated derivative errors, each in units of the fixed
 *   per-component weight w[c], summed. Window-independent normalization, so
 *   scores of different windows compare directly. With w_set FALSE the
 *   weights are first filled from this probe.
 * - q: the worst q_fit over components, the fraction of coefficient mass in
 *   the top quarter of the series. A resolved window sits at round-off
 *   (q ~ 1e-12); an unresolved one at q >~ 1e-3.
 *
 * Returns FALSE (score GSL_POSINF) when the function is non-finite on the window.
 */
static gboolean
_ncm_diff_sc_probe (NcmDiffSCData *data, const gdouble x, const gdouble R, const guint dim,
                    const guint order, NcmVector *y_v, NcmMatrix *fvals, NcmMatrix *coeffs,
                    gdouble *w, gboolean w_set, gdouble *score, gdouble *q)
{
  gdouble score_sum = 0.0;
  gdouble q_max     = 0.0;
  guint c;

  score[0] = GSL_POSINF;
  q[0]     = 1.0;

  if (!_ncm_diff_sc_eval_nodes (data, x, R, dim, y_v, fvals, 0, NCM_DIFF_SC_PROBE_N))
    return FALSE;

  _ncm_diff_sc_dct (fvals, dim, NCM_DIFF_SC_PROBE_N, coeffs);

  for (c = 0; c < dim; c++)
  {
    gdouble deriv, tail, roff, q_fit, err;

    _ncm_diff_sc_deriv_est (ncm_matrix_ptr (coeffs, c, 0), NCM_DIFF_SC_PROBE_N, R, order, &deriv, &tail, &roff, &q_fit);

    err = tail + roff;

    if (!w_set)
      w[c] = fabs (deriv) + err + GSL_DBL_MIN;

    score_sum += err / w[c];
    q_max      = GSL_MAX (q_max, q_fit);
  }

  score[0] = score_sum;
  q[0]     = q_max;

  return TRUE;
}

/*
 * Window search: hill descent on the probe score over a dyadic ladder of
 * half-widths, shrinking first and expanding only when shrinking never beat
 * the initial window. While no resolved window (q below NCM_DIFF_SC_BAD_Q)
 * has been seen the shrinking never gives up: for a feature much narrower
 * than the window the score can even move away from the eventual minimum
 * until the window reaches the feature scale. Expansion only makes sense for
 * a resolved fit (it lowers the round-off amplification), so it requires one.
 * Stores the node values of the best window in fvals_best and returns its
 * half-width, or 0.0 when no finite window was found.
 */
static gdouble
_ncm_diff_sc_scan_window (NcmDiffSCData *data, const gdouble x, const gdouble R0, const guint dim,
                          const guint order, NcmVector *y_v, NcmMatrix *fvals, NcmMatrix *coeffs,
                          NcmMatrix *fvals_best, gdouble *w)
{
  gdouble R_best     = 0.0;
  gdouble score_best = GSL_POSINF;
  gdouble q_best     = 1.0;
  gboolean w_set     = FALSE;
  gdouble R          = R0;
  guint worse        = 0;
  guint iter;
  gdouble score, q;

  if (_ncm_diff_sc_probe (data, x, R0, dim, order, y_v, fvals, coeffs, w, w_set, &score, &q))
  {
    w_set      = TRUE;
    score_best = score;
    q_best     = q;
    R_best     = R0;
    ncm_matrix_memcpy (fvals_best, fvals);
  }

  for (iter = 0; iter < NCM_DIFF_SC_MAX_SHRINK; iter++)
  {
    gboolean improved = FALSE;

    if (q_best < 1.0e-14)
      break;

    R *= 0.5;

    if (_ncm_diff_sc_probe (data, x, R, dim, order, y_v, fvals, coeffs, w, w_set, &score, &q))
      w_set = TRUE;

    /* A first resolved window always wins over an unresolved best: the score
     * of an unresolved fit only measures how wrong it knows itself to be. */
    if ((q_best > NCM_DIFF_SC_BAD_Q) && (q <= NCM_DIFF_SC_BAD_Q))
      improved = TRUE;
    else if ((score < NCM_DIFF_SC_IMPROVE * score_best) && (q <= GSL_MAX (q_best, NCM_DIFF_SC_BAD_Q)))
      improved = TRUE;

    if (improved)
    {
      score_best = score;
      q_best     = q;
      R_best     = R;
      worse      = 0;
      ncm_matrix_memcpy (fvals_best, fvals);
    }
    else if (q_best <= NCM_DIFF_SC_BAD_Q)
    {
      worse++;

      if (worse >= NCM_DIFF_SC_PATIENCE)
        break;
    }
  }

  if ((R_best == R0) && (q_best <= NCM_DIFF_SC_BAD_Q))
  {
    R     = R0;
    worse = 0;

    for (iter = 0; iter < NCM_DIFF_SC_MAX_EXPAND; iter++)
    {
      if (q_best < 1.0e-14)
        break;

      R *= 2.0;

      _ncm_diff_sc_probe (data, x, R, dim, order, y_v, fvals, coeffs, w, w_set, &score, &q);

      if ((score < NCM_DIFF_SC_IMPROVE * score_best) && (q <= NCM_DIFF_SC_BAD_Q))
      {
        score_best = score;
        q_best     = q;
        R_best     = R;
        worse      = 0;
        ncm_matrix_memcpy (fvals_best, fvals);
      }
      else
      {
        worse++;

        if (worse >= NCM_DIFF_SC_PATIENCE)
          break;
      }
    }
  }

  return R_best;
}

static GArray *
_ncm_diff_sc_dn (NcmDiff *diff, const guint order, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  const guint nvar            = x_a->len;
  GArray *df                  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v              = ncm_vector_new_array (x_a);
  NcmVector *y_v              = ncm_vector_new (dim);
  NcmMatrix *fvals            = ncm_matrix_new (dim, NCM_DIFF_SC_MAX_N);
  NcmMatrix *coeffs           = ncm_matrix_new (dim, NCM_DIFF_SC_MAX_N);
  NcmMatrix *fvals_probe      = ncm_matrix_new (dim, NCM_DIFF_SC_PROBE_N);
  NcmMatrix *coeffs_probe     = ncm_matrix_new (dim, NCM_DIFF_SC_PROBE_N);
  NcmMatrix *fvals_best       = ncm_matrix_new (dim, NCM_DIFF_SC_PROBE_N);
  GArray *d_prev              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *w_a                 = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *conv                = g_array_new (FALSE, FALSE, sizeof (NcmDiffDualConv));
  NcmDiffSCData data          = {f, user_data, x_v, 0};
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

  g_array_set_size (d_prev, dim);
  g_array_set_size (w_a, dim);
  g_array_set_size (conv, dim);

  for (a = 0; a < nvar; a++)
  {
    const gdouble x     = g_array_index (x_a, gdouble, a);
    const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);
    const gdouble R0    = self->spectral_window * scale;
    gdouble R;
    guint N_cur = 0;
    guint c;

    data.a = a;

    R = _ncm_diff_sc_scan_window (&data, x, R0, dim, order, y_v, fvals_probe, coeffs_probe, fvals_best,
                                  &g_array_index (w_a, gdouble, 0));

    if (R == 0.0)
      g_error ("ncm_diff_sc: no window around x[%u] = % 22.15g with finite function values.", a, x);

    for (c = 0; c < dim; c++)
      _ncm_diff_dual_conv_init (&g_array_index (conv, NcmDiffDualConv, c));

    /* Refinement on nested grids, starting from the best probe. */
    {
      guint N;

      for (c = 0; c < dim; c++)
      {
        guint i;

        for (i = 0; i < NCM_DIFF_SC_PROBE_N; i++)
          ncm_matrix_set (fvals, c, i, ncm_matrix_get (fvals_best, c, i));
      }

      N_cur = NCM_DIFF_SC_PROBE_N;

      for (N = NCM_DIFF_SC_PROBE_N; N <= NCM_DIFF_SC_MAX_N; N = 2 * N - 1)
      {
        gboolean improve  = FALSE;
        gboolean resolved = TRUE;

        if (N > N_cur)
        {
          if (!_ncm_diff_sc_eval_nodes (&data, x, R, dim, y_v, fvals, N_cur, N))
            break;

          N_cur = N;
        }

        _ncm_diff_sc_dct (fvals, dim, N, coeffs);

        for (c = 0; c < dim; c++)
        {
          NcmDiffDualConv *cs = &g_array_index (conv, NcmDiffDualConv, c);
          gdouble d_c, tail_c, roff_c, q_fit_c;

          _ncm_diff_sc_deriv_est (ncm_matrix_ptr (coeffs, c, 0), N, R, order, &d_c, &tail_c, &roff_c, &q_fit_c);

          if (tail_c > roff_c)
            resolved = FALSE;

          if (N == NCM_DIFF_SC_PROBE_N)
          {
            /* No cross-level check yet: record the value, leave the error
             * unknown. It only reaches the caller when the refinement cannot
             * run (non-finite values between the probe nodes). */
            cs->df_best = d_c;
            improve     = TRUE;
          }
          else
          {
            const gdouble err_c = GSL_MAX (fabs (d_c - g_array_index (d_prev, gdouble, c)), tail_c + roff_c);

            if (err_c < cs->err_best)
            {
              cs->df_best  = d_c;
              cs->err_best = err_c;
              improve      = TRUE;
            }
          }

          g_array_index (d_prev, gdouble, c) = d_c;
        }

        /* Once every component's series tail is below its round-off, a finer
         * grid can only add noise. Never stop before one cross-level check. */
        if ((N > NCM_DIFF_SC_PROBE_N) && (!improve || resolved))
          break;
      }
    }

    for (c = 0; c < dim; c++)
    {
      const NcmDiffDualConv *cs = &g_array_index (conv, NcmDiffDualConv, c);

      ncm_matrix_set (df_m, a, c, cs->df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, c, cs->err_best * NCM_DIFF_SC_ERR_PAD);
    }
  }

  {
    g_array_unref (d_prev);
    g_array_unref (w_a);
    g_array_unref (conv);

    ncm_vector_clear (&x_v);
    ncm_vector_clear (&y_v);

    ncm_matrix_clear (&fvals);
    ncm_matrix_clear (&coeffs);
    ncm_matrix_clear (&fvals_probe);
    ncm_matrix_clear (&coeffs_probe);
    ncm_matrix_clear (&fvals_best);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
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
  return ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, _ncm_diff_rf_d1_step, 0, x_a, dim, f, user_data, Eerr);
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
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, _ncm_diff_rf_d1_step, 1, x_a, dim, f, user_data, Eerr);
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
  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, _ncm_diff_rf_d2_step, 1, x_a, dim, f, user_data, Eerr);
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

  df_a =  ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, _ncm_diff_rf_d1_step, 0, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);

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

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, _ncm_diff_rf_d1_step, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);
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

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, _ncm_diff_rf_d2_step, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);
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


  return ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, _ncm_diff_rf_d1_step, 0, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
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


  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, _ncm_diff_rf_d1_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
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


  return ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, _ncm_diff_rf_d2_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
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

  GArray *diag = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, _ncm_diff_rf_d2_step, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, &dEerr);
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

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rf_d1_step, _ncm_diff_rf_d1_step, 0, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

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

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d1_step, _ncm_diff_rf_d1_step, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

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

  df_a = ncm_diff_by_step_algo (diff, _ncm_diff_rc_d2_step, _ncm_diff_rf_d2_step, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);

  g_array_unref (Eerr);

  return df;
}

/**
 * ncm_diff_sc_d1_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f: $\partial_i f$ using the spectral method:
 * for each variable the function is expanded in Chebyshev polynomials on a window
 * whose half-width is searched for automatically (starting from
 * #NcmDiff:spectral-window times the variable scale), and the expansion is
 * differentiated analytically at the window center. The function $f$ is considered
 * as a $f:\mathbb{R}^N\to \mathbb{R}^M$, where $N = $ length of @x_a and $M = $ @dim.
 *
 * Compared to the finite-difference methods, the spectral method samples the
 * function on a wide window instead of a shrinking neighborhood, which lowers the
 * round-off amplification and handles functions whose scale of variation differs
 * from the magnitude of the variable. It uses more function evaluations.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_sc_d1_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return _ncm_diff_sc_dn (diff, 1, x_a, dim, f, user_data, Eerr);
}

/**
 * ncm_diff_sc_d2_N_to_M:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the second derivative of @f: $\partial_i^2 f$ using the spectral
 * method, see ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R}^N\to \mathbb{R}^M$, where $N = $ length of @x_a and $M = $ @dim.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_sc_d2_N_to_M (NcmDiff *diff, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  return _ncm_diff_sc_dn (diff, 2, x_a, dim, f, user_data, Eerr);
}

/**
 * ncm_diff_sc_d1_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the first derivative of @f using the spectral method, see
 * ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R}\to \mathbb{R}^M$, where $M = $ @dim.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x.
 */
GArray *
ncm_diff_sc_d1_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;

  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = _ncm_diff_sc_dn (diff, 1, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);

  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_sc_d2_1_to_M:
 * @diff: a #NcmDiff
 * @x: function argument
 * @dim: dimension of @f
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the second derivative of @f using the spectral method, see
 * ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R}\to \mathbb{R}^M$, where $M = $ @dim.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x.
 */
GArray *
ncm_diff_sc_d2_1_to_M (NcmDiff *diff, const gdouble x, const guint dim, NcmDiffFunc1toM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {f, NULL, NULL, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *df_a;

  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = _ncm_diff_sc_dn (diff, 2, x_a, dim, &_ncm_diff_trans_1_to_M, &fp, Eerr);

  g_array_unref (x_a);

  return df_a;
}

/**
 * ncm_diff_sc_d1_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the gradient of @f using the spectral method, see
 * ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R}^N \to \mathbb{R}$, where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_sc_d1_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};

  return _ncm_diff_sc_dn (diff, 1, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_sc_d2_N_to_1:
 * @diff: a #NcmDiff
 * @x_a: (array) (element-type double) (in): function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @Eerr: (array) (element-type double) (out) (transfer full): estimated errors
 *
 * Calculates the second derivatives $\partial_i^2 f$ using the spectral method,
 * see ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R}^N \to \mathbb{R}$, where $N = $ length of @x_a.
 *
 * Returns: (transfer full) (array) (element-type double): The derivative of @f at @x_a.
 */
GArray *
ncm_diff_sc_d2_N_to_1 (NcmDiff *diff, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffFuncParams fp = {NULL, f, NULL, user_data};

  return _ncm_diff_sc_dn (diff, 2, x_a, 1, &_ncm_diff_trans_N_to_1, &fp, Eerr);
}

/**
 * ncm_diff_sc_d1_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 *
 * Calculates the first derivative of @f using the spectral method, see
 * ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R} \to \mathbb{R}$.
 *
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_sc_d1_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr         = NULL;
  GArray *df_a;
  gdouble df;

  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = _ncm_diff_sc_dn (diff, 1, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);

  g_array_unref (Eerr);

  return df;
}

/**
 * ncm_diff_sc_d2_1_to_1:
 * @diff: a #NcmDiff
 * @x: function argument
 * @f: (scope call): function to differentiate
 * @user_data: (nullable): function user data
 * @err: (out) (nullable): estimated error
 *
 * Calculates the second derivative of @f using the spectral method, see
 * ncm_diff_sc_d1_N_to_M(). The function $f$ is considered as a
 * $f:\mathbb{R} \to \mathbb{R}$.
 *
 * Returns: The derivative of @f at @x.
 */
gdouble
ncm_diff_sc_d2_1_to_1 (NcmDiff *diff, const gdouble x, NcmDiffFunc1to1 f, gpointer user_data, gdouble *err)
{
  NcmDiffFuncParams fp = {NULL, NULL, f, user_data};
  GArray *x_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *Eerr         = NULL;
  GArray *df_a;
  gdouble df;

  g_array_set_size (x_a, 1);
  g_array_index (x_a, gdouble, 0) = x;

  df_a = _ncm_diff_sc_dn (diff, 2, x_a, 1, &_ncm_diff_trans_1_to_1, &fp, &Eerr);

  df = g_array_index (df_a, gdouble, 0);

  g_array_unref (x_a);
  g_array_unref (df_a);

  if (err != NULL)
    *err = g_array_index (Eerr, gdouble, 0);

  g_array_unref (Eerr);

  return df;
}


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
 * Computes first and second derivatives and Hessians of user functions by
 * finite differences plus Richardson extrapolation, and first and second
 * derivatives by Chebyshev fits (the ncm_diff_sc_ functions). The methods,
 * their error estimates and the step control are derived in
 * <a href="../../theory/ncm/integration/diff.html">Numerical Differentiation</a>.
 *
 * The extrapolation uses precomputed tables. For a step sequence
 * $g_k$ (powers of #NcmDiff:richardson-step), order $n$ holds the nodes
 * $h_i = 1/g_i$, $i = 0, \dots, n+1$, and the Lagrange weights
 * $\lambda_i = \prod_{j \neq i} 1/(1 - g_j/g_i)$ that extrapolate a
 * polynomial in $h$ (forward, $g_k = r_s^k$) or $h^2$ (central,
 * $g_k = r_s^{2k}$) to $h \to 0$. The extrapolated derivative at order $n$
 * is $\sum_i \lambda_i D(h_0 h_i)$ (forward) or
 * $\sum_i \lambda_i D(h_0 \sqrt{h_i})$ (central), where $D(h)$ is the
 * finite-difference quotient at step $h$.
 *
 * By default every derivative runs two Richardson ladders whose steps are not
 * commensurate, see #NcmDiff:dual-series: their disagreement at the same order
 * measures the truncation error and the scatter of the values of $f$, and a
 * result is accepted only where they agree, which rejects a plateau that one
 * ladder alone would take for convergence when it aliases an oscillation of
 * $f$. With #NcmDiff:dual-series off a single ladder is used: it costs about
 * two thirds of the evaluations, and relies on heuristics to reject such
 * plateaus.
 *
 * A single ladder carries an error estimate combining the truncation error
 * (difference between consecutive extrapolation orders, times
 * #NcmDiff:trunc-change-ratio) and the cancellation scale of the difference
 * quotients, the size of the error of the values of $f$ that survives their
 * subtraction (for values of relative precision #NcmDiff:func-precision and
 * absolute precision #NcmDiff:func-abs-precision). Two orders agree when
 * they differ by less than a thousandth of the largest difference quotient
 * of the row and the cancellation scale of the row is below that level; a zero
 * derivative then converges like any other. Orders are increased while any
 * component still improves; the best value per component is returned.
 *
 * The initial step is #NcmDiff:ini-h times $|x|$, or #NcmDiff:ini-h at
 * $x = 0$. When the cancellation scale of the first row exceeds a thousandth
 * of its largest difference quotient the step is
 * increased: to #NcmDiff:ini-h first when $|x| < 1$, then by factors of
 * #NcmDiff:richardson-step, up to $r_s^3$ #NcmDiff:ini-h $\max(1, |x|)$.
 * With a domain set by ncm_diff_set_domain() no point is evaluated outside
 * it; without one every point is evaluated, including points across an
 * edge of the domain of $f$.
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
  gdouble trunc_ratio;
  gdouble func_prec;
  gdouble canc_pad;
  gdouble func_abs_prec;
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
  PROP_FUNC_PRECISION,
  PROP_TRUNC_CHANGE_RATIO,
  PROP_INI_H,
  PROP_DUAL_SERIES,
  PROP_SPECTRAL_WINDOW,
  PROP_DOMAIN_WARNINGS,
  PROP_FUNC_ABS_PRECISION,
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
  self->trunc_ratio     = 0.0;
  self->func_prec       = 0.0;
  self->canc_pad        = 0.0;
  self->func_abs_prec   = 0.0;
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
    case PROP_FUNC_PRECISION:
      ncm_diff_set_func_precision (diff, g_value_get_double (value));
      break;
    case PROP_TRUNC_CHANGE_RATIO:
      ncm_diff_set_trunc_change_ratio (diff, g_value_get_double (value));
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
    case PROP_FUNC_ABS_PRECISION:
      ncm_diff_set_func_abs_precision (diff, g_value_get_double (value));
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
    case PROP_FUNC_PRECISION:
      g_value_set_double (value, ncm_diff_get_func_precision (diff));
      break;
    case PROP_TRUNC_CHANGE_RATIO:
      g_value_set_double (value, ncm_diff_get_trunc_change_ratio (diff));
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
    case PROP_FUNC_ABS_PRECISION:
      g_value_set_double (value, ncm_diff_get_func_abs_precision (diff));
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
                                                      "Maximum extrapolation order",
                                                      1, G_MAXUINT, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RS,
                                   g_param_spec_double ("richardson-step",
                                                        NULL,
                                                        "Ratio between consecutive steps",
                                                        1.1, G_MAXDOUBLE, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FUNC_PRECISION,
                                   g_param_spec_double ("func-precision",
                                                        NULL,
                                                        "Relative precision of the values of f",
                                                        0.5 * GSL_DBL_EPSILON, 1.0, 1.5e4 * GSL_DBL_EPSILON,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TRUNC_CHANGE_RATIO,
                                   g_param_spec_double ("trunc-change-ratio",
                                                        NULL,
                                                        "Largest ratio of a row's truncation error to its change from the previous row",
                                                        1.1, G_MAXDOUBLE, 3.0e4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INI_H,
                                   g_param_spec_double ("ini-h",
                                                        NULL,
                                                        "Initial step relative to |x|",
                                                        GSL_DBL_EPSILON, G_MAXDOUBLE, pow (GSL_DBL_EPSILON, 1.0 / 8.0),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DUAL_SERIES,
                                   g_param_spec_boolean ("dual-series",
                                                         NULL,
                                                         "Use two extrapolation series",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SPECTRAL_WINDOW,
                                   g_param_spec_double ("spectral-window",
                                                        NULL,
                                                        "Largest spectral window half-width in units of max (1, |x|)",
                                                        GSL_DBL_EPSILON, G_MAXDOUBLE, 256.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DOMAIN_WARNINGS,
                                   g_param_spec_boolean ("domain-warnings",
                                                         NULL,
                                                         "Warn when a central difference falls back to a one-sided one at an edge of the domain",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FUNC_ABS_PRECISION,
                                   g_param_spec_double ("func-abs-precision",
                                                        NULL,
                                                        "Absolute precision of the values of f",
                                                        0.0, G_MAXDOUBLE, 0.0,
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
 * Increases the reference count of @diff by one.
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
 * Decreases the reference count of @diff by one.
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
 * Decreases the reference count of @diff by one and sets *@diff to %NULL.
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
 * Gets #NcmDiff:max-order.
 *
 * Returns: the maximum extrapolation order.
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
 * Gets #NcmDiff:richardson-step.
 *
 * Returns: the ratio $r_s$ between consecutive steps.
 */
gdouble
ncm_diff_get_richardson_step (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->rs;
}

/**
 * ncm_diff_get_func_precision:
 * @diff: a #NcmDiff
 *
 * Gets #NcmDiff:func-precision.
 *
 * Returns: the relative precision of the values of $f$.
 */
gdouble
ncm_diff_get_func_precision (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->func_prec;
}

/**
 * ncm_diff_get_trunc_change_ratio:
 * @diff: a #NcmDiff
 *
 * Gets #NcmDiff:trunc-change-ratio.
 *
 * Returns: the largest ratio of a row's truncation error to its change from the previous row.
 */
gdouble
ncm_diff_get_trunc_change_ratio (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->trunc_ratio;
}

/**
 * ncm_diff_get_ini_h:
 * @diff: a #NcmDiff
 *
 * Gets #NcmDiff:ini-h.
 *
 * Returns: the initial step relative to $|x|$.
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
 * @maxorder: maximum extrapolation order
 *
 * Sets #NcmDiff:max-order; order $n$ uses $n + 2$ steps. Requires
 * @maxorder $\geq 1$.
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
 * @rs: ratio between consecutive steps
 *
 * Sets #NcmDiff:richardson-step, the ratio $r_s$ between consecutive steps
 * of the ladders, and rebuilds the tables. Requires @rs $\geq 1.1$.
 */
void
ncm_diff_set_richardson_step (NcmDiff *diff, const gdouble rs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (rs, >=, 1.1);

  if (rs != self->rs)
  {
    self->rs = rs;

    if (self->central_tables->len > 0)
      _ncm_diff_build_diff_tables (diff);
  }
}

/**
 * ncm_diff_set_func_precision:
 * @diff: a #NcmDiff
 * @func_prec: relative precision of the values of $f$
 *
 * Sets #NcmDiff:func-precision, a bound $s$ on the relative error of each
 * value of $f$ around a smooth function: $\epsilon/2$ for a correctly
 * rounded value, the tolerance of the algorithm for a value computed by one.
 * The Richardson methods scale the cancellation scale of each difference
 * quotient, computed for values correct to $\epsilon/2$, by $2s/\epsilon$.
 * The spectral methods do not use it: their coefficient error assumes
 * values correct to $\epsilon/2$, and a larger relative scatter of the
 * values shows in the coefficient tail of the fit.
 * The default, $1.5 \times 10^4 \epsilon \approx 3.3 \times 10^{-12}$, is the
 * precision assumed for a function nothing is stated about; it covers, for
 * example, $\sin(w x)$ at $|w x|$ up to $10^4$, whose rounded argument
 * moves the values by up to $10^{-12}$. For a function accurate to a few
 * ulps the estimates are then about $10^4$ times its actual error. A
 * smaller stated precision gives tighter estimates and also lowers the
 * scale that keeps the ladder from stopping early, see
 * #NcmDiff:trunc-change-ratio. Requires $\epsilon/2 \leq$ @func_prec $\leq 1$.
 */
void
ncm_diff_set_func_precision (NcmDiff *diff, const gdouble func_prec)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (func_prec, >=, 0.5 * GSL_DBL_EPSILON);
  g_assert_cmpfloat (func_prec, <=, 1.0);
  self->func_prec = func_prec;
  self->canc_pad  = 2.0 * func_prec / GSL_DBL_EPSILON;
}

/**
 * ncm_diff_set_trunc_change_ratio:
 * @diff: a #NcmDiff
 * @trunc_ratio: largest ratio of a row's truncation error to its change
 *
 * Sets #NcmDiff:trunc-change-ratio, $P_{\rm t}$, the largest ratio assumed
 * between the truncation error of a Richardson row and its change from the
 * previous row. Once the rows are asymptotic the change is the truncation
 * error of the previous row and already exceeds that of the current one;
 * rows that agree before that, early or on a plateau aliased with an
 * oscillation of $f$, can have a change below their error. $P_{\rm t}$
 * times the change enters the error estimate of each row, and through the
 * best error the rule that keeps the ladder going past an agreement, which
 * is what lets it leave an aliased plateau. It applies to the single
 * ladder only: with #NcmDiff:dual-series the disagreement of the two
 * ladders at the same order replaces it in the error estimate, and their
 * agreement gate in the protection against aliased plateaus. Requires
 * @trunc_ratio $\geq 1.1$.
 */
void
ncm_diff_set_trunc_change_ratio (NcmDiff *diff, const gdouble trunc_ratio)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (trunc_ratio, >=, 1.1);
  self->trunc_ratio = trunc_ratio;
}

/**
 * ncm_diff_set_ini_h:
 * @diff: a #NcmDiff
 * @ini_h: initial step relative to $|x|$
 *
 * Sets #NcmDiff:ini-h: the first step is @ini_h $|x|$, or @ini_h at
 * $x = 0$. Requires @ini_h $\geq$ %GSL_DBL_EPSILON.
 */
void
ncm_diff_set_ini_h (NcmDiff *diff, const gdouble ini_h)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (ini_h, >=, GSL_DBL_EPSILON);
  self->ini_h = ini_h;
}

/**
 * ncm_diff_set_dual_series:
 * @diff: a #NcmDiff
 * @dual_series: whether to use two extrapolation series
 *
 * Enables or disables the dual-series scheme, on by default. When enabled,
 * every derivative runs two Richardson extrapolation series, started from
 * the initial steps $h_0$ and $h_0/\sqrt{r_s}$, where $r_s$ is
 * #NcmDiff:richardson-step. The difference between the two series at the
 * same order is the truncation error estimate, in place of the difference
 * between consecutive orders, and a row is accepted only where they agree:
 * a single series can converge on an oscillation of $f$ it aliases, the
 * other does not alias it the same way. Disabled, a single series is used,
 * with about two thirds of the function evaluations, and it rejects aliased
 * plateaus by heuristics, the lead share of a row and
 * #NcmDiff:trunc-change-ratio.
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
 * @spectral_window: the largest spectral window half-width
 *
 * Sets the largest half-width of the window used by the spectral methods
 * (ncm_diff_sc_d1_N_to_M() and related), in units of $\max(1, |x|)$: no
 * point farther than @spectral_window $\max(1, |x|)$ from $x$ is evaluated.
 * The search starts from #NcmDiff:ini-h times $|x|$ (or #NcmDiff:ini-h at
 * $x = 0$) and grows the window by factors of 4 until its highest
 * coefficients exceed the rounding error of the values; a 17-node fit then expands or shrinks
 * it to reduce the estimated derivative error. With a domain set the window
 * also stays inside it.
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
 * @lb: (nullable): lower bounds of the coordinates
 * @ub: (nullable): upper bounds of the coordinates
 *
 * Sets the domain of the functions to differentiate: no point is evaluated
 * outside $[\mathrm{lb}, \mathrm{ub}]$. An entry $-\infty$ in @lb or
 * $+\infty$ in @ub is no bound. Coordinates are matched by index; the scalar
 * functions use entry zero. A central difference without room for its steps
 * falls back to a forward or backward one toward the farther edge, with a
 * warning unless #NcmDiff:domain-warnings is %FALSE; a domain with room for
 * no useful step on either side is an error. %NULL for both clears the
 * domain.
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
 * Returns: whether a fallback to a one-sided difference is warned about.
 */
gboolean
ncm_diff_get_domain_warnings (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->domain_warnings;
}

/**
 * ncm_diff_set_func_abs_precision:
 * @diff: a #NcmDiff
 * @func_abs_prec: absolute precision of the values of $f$
 *
 * Sets #NcmDiff:func-abs-precision, a bound on the absolute error of each
 * value of $f$ around a smooth function, for values whose error does not
 * scale with their size: a function that subtracts close terms inside, such
 * as $\ln(1 + x)$ at small $x$, or one computed by an algorithm run with an
 * absolute tolerance. The cancellation scale of each difference quotient
 * includes the error these values leave in it, added to the one from
 * #NcmDiff:func-precision. The default, zero, adds nothing. Requires
 * @func_abs_prec $\geq 0$.
 */
void
ncm_diff_set_func_abs_precision (NcmDiff *diff, const gdouble func_abs_prec)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  g_assert_cmpfloat (func_abs_prec, >=, 0.0);
  g_assert (gsl_finite (func_abs_prec));
  self->func_abs_prec = func_abs_prec;
}

/**
 * ncm_diff_get_func_abs_precision:
 * @diff: a #NcmDiff
 *
 * Gets #NcmDiff:func-abs-precision.
 *
 * Returns: the absolute precision of the values of $f$.
 */
gdouble
ncm_diff_get_func_abs_precision (NcmDiff *diff)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  return self->func_abs_prec;
}

/**
 * ncm_diff_get_spectral_window:
 * @diff: a #NcmDiff
 *
 * Gets the largest spectral window half-width, see ncm_diff_set_spectral_window().
 *
 * Returns: the largest spectral window half-width, in units of $\max(1, |x|)$.
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
 * Logs the steps and weights of every order of the central tables.
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
 * Logs the steps and weights of every order of the forward tables.
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
 * Logs the steps and weights of every order of the backward tables.
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

typedef void (*NcmDiffStepAlgo) (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot, NcmVector *canc, gdouble *canc_abs);
typedef void (*NcmDiffHessianStepAlgo) (NcmDiff *diff, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hx, const guint b, const gdouble y, const gdouble hy, NcmVector *x_v, const gdouble fval, gdouble *quot, gdouble *canc, gdouble *canc_abs);

#define NCM_DIFF_ERR_PAD (1.0e0)
#define NCM_DIFF_NTRY_CONV (3)

/*
 * Largest share of a row its largest term may carry for the row to agree or
 * to replace a best value from such rows, see _ncm_diff_ladder_accum().
 */
#define NCM_DIFF_LEAD_SHARE_MAX (0.99)

/*
 * Control of one component in the dual-series scheme. Series A and B share
 * the tables and start from different steps: their difference at equal
 * order estimates the truncation error, and the cancellation scale of both is added.
 * It stops after NCM_DIFF_DUAL_MAX_WORSE orders without improvement.
 */
#define NCM_DIFF_DUAL_ERR_PAD (10.0)
#define NCM_DIFF_DUAL_MAX_WORSE (2)
#define NCM_DIFF_DUAL_MIN_ORDER (3)

#define NCM_DIFF_MOVE_UP_GUARD (3)

/*
 * Control of one component. At each order it receives the current and
 * previous rows with their cancellation scales, updates the best value and error, and
 * reports whether a higher order is needed.
 */
typedef struct _NcmDiffControl
{
  gdouble df_best;
  gdouble err_best;
  gdouble err_last_max;
  guchar agreements_left;
  guchar just_converged;
  gboolean informative;
  gboolean best_eligible;
} NcmDiffControl;

static void
_ncm_diff_control_init (NcmDiffControl *cs)
{
  cs->df_best         = 0.0;
  cs->err_best        = GSL_POSINF;
  cs->err_last_max    = 0.0;
  cs->agreements_left = NCM_DIFF_NTRY_CONV;
  cs->just_converged  = 0;
  cs->informative     = FALSE;
  cs->best_eligible   = FALSE;
}

typedef struct _NcmDiffCrossControl
{
  gdouble df_best;
  gdouble err_best;
  gdouble rel_err_best;
  gdouble cross_best;
  gdouble row_lo;
  gdouble row_hi;
  gdouble elig_lo;
  gdouble elig_hi;
  GArray *elig_val;
  GArray *elig_err;
  guchar n_worse;
  gboolean best_eligible;
} NcmDiffCrossControl;

static void
_ncm_diff_cross_control_init (NcmDiffCrossControl *cs)
{
  cs->df_best       = 0.0;
  cs->err_best      = GSL_POSINF;
  cs->rel_err_best  = GSL_POSINF;
  cs->cross_best    = GSL_POSINF;
  cs->row_lo        = GSL_POSINF;
  cs->row_hi        = GSL_NEGINF;
  cs->elig_lo       = GSL_POSINF;
  cs->elig_hi       = GSL_NEGINF;
  cs->n_worse       = 0;
  cs->best_eligible = FALSE;
  cs->elig_val      = NULL;
  cs->elig_err      = NULL;
}

/*
 * The error of the best value. Two ladders can agree on a row while both are
 * off, when the leading truncation term nearly vanishes at x, so the error
 * covers the spread around the best of the agreeing rows that follow it
 * (earlier ones agree worse, by the ranking, and pass the gate only because
 * it is relative to the largest quotient, large near a zero derivative); without any
 * agreeing row nothing establishes which row is right, and it covers the
 * spread of every row seen. The spreads enter only the reported error: the
 * rows are ranked and the ladders stopped by the errors of the rows alone.
 *
 * Every agreeing row k also bounds the best value through its own error,
 * |best - d| <= |best - R_k| + err_k, so with an agreeing best the error is
 * the smallest of these bounds, at most that of the best row, and at least
 * the spread of the rows that follow it.
 */
static gdouble
_ncm_diff_cross_control_err (const NcmDiffCrossControl *cs)
{
  const gdouble lo     = cs->best_eligible ? cs->elig_lo : cs->row_lo;
  const gdouble hi     = cs->best_eligible ? cs->elig_hi : cs->row_hi;
  const gdouble spread = NCM_DIFF_DUAL_ERR_PAD * GSL_MAX (hi - cs->df_best, cs->df_best - lo);
  gdouble bound        = cs->err_best;

  if (!cs->best_eligible)
    return GSL_MAX (cs->err_best, spread);

  if (cs->elig_val != NULL)
  {
    guint k;

    for (k = 0; k < cs->elig_val->len; k++)
      bound = GSL_MIN (bound, g_array_index (cs->elig_err, gdouble, k) + fabs (cs->df_best - g_array_index (cs->elig_val, gdouble, k)));
  }

  return GSL_MAX (bound, spread);
}

/*
 * One Richardson ladder for one scalar component: the difference quotient
 * and cancellation scale of each step h0 rs^-t, the extrapolated rows and the
 * convergence control. The vector drivers keep one ladder per component and
 * append every evaluated step to all of them; a converged ladder is no longer
 * updated.
 */
typedef struct _NcmDiffLadder
{
  GPtrArray *tables;
  guint order;
  GArray *quots;
  GArray *cancs;
  GArray *cancs_abs;
  gdouble row;
  gdouble row_prev;
  gdouble row_canc;
  gdouble row_canc_prev;
  gdouble row_canc_abs;
  gdouble row_canc_abs_prev;
  gdouble quot_scale;
  gdouble max_term;
  gdouble rho;
  gdouble lead_share;
  NcmDiffControl control;
  gboolean converged;
  gboolean first_informative;
} NcmDiffLadder;

static void
_ncm_diff_ladder_init (NcmDiffLadder *ladder, GPtrArray *tables)
{
  ladder->tables    = tables;
  ladder->quots     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ladder->cancs     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ladder->cancs_abs = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void
_ncm_diff_ladder_clear (NcmDiffLadder *ladder)
{
  g_clear_pointer (&ladder->quots, g_array_unref);
  g_clear_pointer (&ladder->cancs, g_array_unref);
  g_clear_pointer (&ladder->cancs_abs, g_array_unref);
}

/* Empties the ladder: no steps, no rows, the control reset. */
static void
_ncm_diff_ladder_reset (NcmDiffLadder *ladder)
{
  ladder->order         = 0;
  ladder->row           = 0.0;
  ladder->row_prev      = 0.0;
  ladder->row_canc      = 0.0;
  ladder->row_canc_prev = 0.0;

  ladder->row_canc_abs      = 0.0;
  ladder->row_canc_abs_prev = 0.0;
  ladder->quot_scale        = 0.0;
  ladder->max_term          = 0.0;
  ladder->rho               = 0.0;
  ladder->lead_share        = 0.0;
  ladder->converged         = FALSE;

  ladder->first_informative = FALSE;

  g_array_set_size (ladder->quots, 0);
  g_array_set_size (ladder->cancs, 0);
  g_array_set_size (ladder->cancs_abs, 0);

  _ncm_diff_control_init (&ladder->control);
}

/* Restarts the rows and the control, keeping the steps. */
static void
_ncm_diff_ladder_restart (NcmDiffLadder *ladder)
{
  ladder->order         = 0;
  ladder->row           = 0.0;
  ladder->row_prev      = 0.0;
  ladder->row_canc      = 0.0;
  ladder->row_canc_prev = 0.0;

  ladder->row_canc_abs      = 0.0;
  ladder->row_canc_abs_prev = 0.0;
  ladder->quot_scale        = 0.0;
  ladder->max_term          = 0.0;
  ladder->rho               = 0.0;
  ladder->lead_share        = 0.0;
  ladder->converged         = FALSE;

  ladder->first_informative = FALSE;

  _ncm_diff_control_init (&ladder->control);
}

static void
_ncm_diff_ladder_add_step (NcmDiffLadder *ladder, const gdouble quot, const gdouble canc, const gdouble canc_abs)
{
  g_array_append_val (ladder->quots, quot);
  g_array_append_val (ladder->cancs, canc);
  g_array_append_val (ladder->cancs_abs, canc_abs);
}

/* A step larger than every stored one becomes the new top of the ladder. */
static void
_ncm_diff_ladder_prepend_step (NcmDiffLadder *ladder, const gdouble quot, const gdouble canc, const gdouble canc_abs)
{
  g_array_prepend_val (ladder->quots, quot);
  g_array_prepend_val (ladder->cancs, canc);
  g_array_prepend_val (ladder->cancs_abs, canc_abs);
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
  NcmDiffCrossControl cross;
  GArray *elig_val;
  GArray *elig_err;
  gboolean converged;
} NcmDiffDual;

/* Restarts the cross control on the arrays of agreeing rows the dual owns, emptied. */
static void
_ncm_diff_dual_cross_restart (NcmDiffDual *dual)
{
  _ncm_diff_cross_control_init (&dual->cross);
  g_array_set_size (dual->elig_val, 0);
  g_array_set_size (dual->elig_err, 0);
  dual->cross.elig_val = dual->elig_val;
  dual->cross.elig_err = dual->elig_err;
}

static void
_ncm_diff_dual_init (NcmDiffDual *dual, GPtrArray *tables)
{
  _ncm_diff_ladder_init (&dual->A, tables);
  _ncm_diff_ladder_init (&dual->B, tables);
  dual->elig_val = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dual->elig_err = g_array_new (FALSE, FALSE, sizeof (gdouble));
  _ncm_diff_dual_cross_restart (dual);
}

static void
_ncm_diff_dual_clear (NcmDiffDual *dual)
{
  _ncm_diff_ladder_clear (&dual->A);
  _ncm_diff_ladder_clear (&dual->B);
  g_clear_pointer (&dual->elig_val, g_array_unref);
  g_clear_pointer (&dual->elig_err, g_array_unref);
}

static void
_ncm_diff_dual_reset (NcmDiffDual *dual)
{
  _ncm_diff_ladder_reset (&dual->A);
  _ncm_diff_ladder_reset (&dual->B);
  _ncm_diff_dual_cross_restart (dual);
  dual->converged = FALSE;
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

static GArray *_ncm_diff_by_step_algo_single (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr);
static GArray *_ncm_diff_by_step_algo_dual (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr);

static GArray *
ncm_diff_by_step_algo (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  if (self->dual_series)
    return _ncm_diff_by_step_algo_dual (diff, step_algo, onesided_algo, po, x_a, dim, f, user_data, Eerr);

  return _ncm_diff_by_step_algo_single (diff, step_algo, onesided_algo, po, x_a, dim, f, user_data, Eerr);
}

static void _ncm_diff_scheme_init (NcmDiffPrivate *self, NcmDiffScheme *sch, NcmDiffStepAlgo algo, NcmDiffStepAlgo onesided, const guint po, const guint a, const gdouble x);
static gboolean _ncm_diff_ladder_needs_step (NcmDiffLadder *ladder);
static gdouble _ncm_diff_step_h (NcmDiffPrivate *self, const guint po, const gdouble h0, const gint k);
static void _ncm_diff_eval_step (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble ho, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot_t, NcmVector *canc_t, gdouble *canc_abs_t);
static void _ncm_diff_ladder_extrapolate (NcmDiffLadder *ladder, const gdouble trunc_ratio, const gdouble canc_pad);
static gdouble _ncm_diff_scheme_cap (NcmDiffPrivate *self, const NcmDiffScheme *sch, const gdouble x);
static gboolean _ncm_diff_scheme_fallback (NcmDiffPrivate *self, NcmDiffScheme *sch, const guint a, const gdouble x);
static gdouble _ncm_diff_scheme_room (const NcmDiffScheme *sch);
static gboolean _ncm_diff_vector_finite (NcmVector *v);
static void _ncm_diff_ladder_replay (NcmDiffLadder *ladder, const gdouble trunc_ratio, const gdouble canc_pad);

static GArray *
_ncm_diff_by_step_algo_single (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GArray *ladders             = g_array_new (FALSE, FALSE, sizeof (NcmDiffLadder));
  NcmVector *x_v              = NULL;
  NcmVector *f_v              = NULL;
  NcmVector *yh1_v            = NULL;
  NcmVector *yh2_v            = NULL;
  NcmVector *quot_t           = NULL;
  NcmVector *canc_t           = NULL;
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
  quot_t = ncm_vector_new (dim);
  canc_t = ncm_vector_new (dim);

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
        const guint nsteps = g_array_index (ladders, NcmDiffLadder, 0).quots->len;

        gdouble canc_abs_t;

        _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * _ncm_diff_step_h (self, sch.po, sch.h0, k_top + nsteps),
                             x_v, f_v, yh1_v, yh2_v, quot_t, canc_t, &canc_abs_t);

        for (i = 0; i < dim; i++)
          _ncm_diff_ladder_add_step (&g_array_index (ladders, NcmDiffLadder, i), ncm_vector_get (quot_t, i), ncm_vector_get (canc_t, i), canc_abs_t);

        continue;
      }

      for (i = 0; i < dim; i++)
      {
        NcmDiffLadder *ladder = &g_array_index (ladders, NcmDiffLadder, i);

        if (!ladder->converged)
          _ncm_diff_ladder_extrapolate (ladder, self->trunc_ratio, self->canc_pad);
      }

      /*
       * A first row that is not informative means the initial step is far
       * below the scale f varies on and the quotients are dominated by cancellation. The
       * first move sets the step to the one a zero coordinate gets and
       * restarts the ladders there, since the stored steps are not
       * consecutive with it; each further move prepends a step one
       * Richardson factor larger, within the guard, and replays the
       * controls. A move whose evaluation is not finite is rejected. When
       * the domain leaves no room for the move the scheme falls back to a
       * one-sided one and restarts; a domain without room on either side is
       * an error.
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
                     "below the initial step %.3e, and the quotients are dominated by cancellation.",
                     a, x, _ncm_diff_scheme_room (&sch), self->ini_h * scale);
        }
        else
        {
          gdouble canc_abs_t;

          _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * h_new, x_v, f_v, yh1_v, yh2_v, quot_t, canc_t, &canc_abs_t);

          if (!_ncm_diff_vector_finite (quot_t))
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
                _ncm_diff_ladder_add_step (ladder, ncm_vector_get (quot_t, i), ncm_vector_get (canc_t, i), canc_abs_t);
              }
              else
              {
                _ncm_diff_ladder_prepend_step (ladder, ncm_vector_get (quot_t, i), ncm_vector_get (canc_t, i), canc_abs_t);
                _ncm_diff_ladder_replay (ladder, self->trunc_ratio, self->canc_pad);
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

      ncm_matrix_set (df_m, a, i, ladder->control.df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, i, ladder->control.err_best);
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
    ncm_vector_clear (&quot_t);
    ncm_vector_clear (&canc_t);

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static void _ncm_diff_room (NcmDiffPrivate *self, const guint a, const gdouble x, gdouble *lo, gdouble *hi);
static void _ncm_diff_rf_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot, NcmVector *canc, gdouble *canc_abs);

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

static gdouble _ncm_diff_move_up_cap (NcmDiffPrivate *self, const gdouble x);

/* Largest step of the scheme: the smaller of the move-up guard and the room, with a margin that keeps the points strictly inside the domain. */
static gdouble
_ncm_diff_scheme_cap (NcmDiffPrivate *self, const NcmDiffScheme *sch, const gdouble x)
{
  const gdouble margin = 1.0 - 1.0 / self->rs;
  /* The forward second difference evaluates x + 2h. */
  const gdouble fraction = (sch->algo == _ncm_diff_rf_d2_step) ? GSL_MIN (margin, 0.5) : margin;

  return GSL_MIN (_ncm_diff_move_up_cap (self, x), _ncm_diff_scheme_room (sch) * fraction);
}

/*
 * Largest step a move up may take along a coordinate,
 * NCM_DIFF_MOVE_UP_GUARD Richardson factors above the step a zero
 * coordinate gets. Edges are enforced only when the caller sets a domain,
 * see ncm_diff_set_domain().
 */
static gdouble
_ncm_diff_move_up_cap (NcmDiffPrivate *self, const gdouble x)
{
  const gdouble scale = (x == 0.0) ? 1.0 : fabs (x);

  return self->ini_h * GSL_MAX (1.0, scale) * pow (self->rs, NCM_DIFF_MOVE_UP_GUARD);
}

/* Room of the scheme: the distance to the nearer edge for a central one, to the edge its steps point to for a one-sided one. */
static gdouble
_ncm_diff_scheme_room (const NcmDiffScheme *sch)
{
  if (sch->central)
    return GSL_MIN (sch->lo, sch->hi);

  return (sch->sign > 0.0) ? sch->hi : sch->lo;
}

static void _ncm_diff_scheme_warn (NcmDiffPrivate *self, const NcmDiffScheme *sch, const guint a, const gdouble x, const gchar *from);

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

static void
_ncm_diff_scheme_warn (NcmDiffPrivate *self, const NcmDiffScheme *sch, const guint a, const gdouble x, const gchar *from)
{
  if (self->domain_warnings)
    g_warning ("NcmDiff: coordinate %u at % .15g is %.3e from the nearer edge of its domain, "
               "falling back from a %s difference to a %s one.",
               a, x, GSL_MIN (sch->lo, sch->hi), from, (sch->sign > 0.0) ? "forward" : "backward");
}

/* Extrapolation order k uses the first k + 2 steps. */
static gboolean
_ncm_diff_ladder_needs_step (NcmDiffLadder *ladder)
{
  return !ladder->converged && (ladder->quots->len < ladder->order + 2);
}

/*
 * Step of index k, h0 rs^-k, through the formula of the tables so that
 * k >= 0 reproduces their steps; k < 0 are the steps above the initial one.
 */
static gdouble
_ncm_diff_step_h (NcmDiffPrivate *self, const guint po, const gdouble h0, const gint k)
{
  /* Above the initial step rs^(-2k) overflows for large -k; h0 rs^(-k) is the same step. */
  if (k < 0)
    return h0 * pow (self->rs, -k);
  else if (po == 0)
    return h0 * (1.0 / pow (self->rs, k));
  else
    return h0 * sqrt (1.0 / pow (self->rs, 2 * k));
}

/* Evaluates the step ho along coordinate a into quot_t and canc_t, restoring x_v. */
static void
_ncm_diff_eval_step (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffFuncNtoM f, gpointer user_data,
                     const guint a, const gdouble x, const gdouble ho,
                     NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot_t, NcmVector *canc_t, gdouble *canc_abs_t)
{
  volatile gdouble temp = x + ho;
  const gdouble h       = temp - x;

  step_algo (diff, f, user_data, a, x, h, x_v, f_v, yh1_v, yh2_v, quot_t, canc_t, canc_abs_t);
  ncm_vector_set (x_v, a, x);
}

static void _ncm_diff_ladder_accum (NcmDiffLadder *ladder, NcmDiffTable *dtable, const guint nt);
static gboolean _ncm_diff_control_update (NcmDiffControl *cs, const gdouble trunc_ratio, const gdouble canc_pad, const gdouble row, const gdouble row_prev, const gdouble row_canc, const gdouble row_canc_prev, const gdouble row_canc_abs, const gdouble row_canc_abs_prev, const gdouble quot_scale, const gdouble lead_share);

/*
 * Extrapolates at the current order, updates the control unless the ladder
 * has converged, and moves to the next order. The rows of a converged ladder
 * are still computed for the controls that read them.
 */
static void
_ncm_diff_ladder_extrapolate (NcmDiffLadder *ladder, const gdouble trunc_ratio, const gdouble canc_pad)
{
  NcmDiffTable *dtable;

  /* No table for this order, or not enough steps for its row: a ladder that
   * converged at the last row of a replay. */
  if ((ladder->order == ladder->tables->len) || (ladder->order + 2 > ladder->quots->len))
    return;

  dtable = g_ptr_array_index (ladder->tables, ladder->order);

  if (ladder->order == 0)
  {
    ladder->row_prev          = g_array_index (ladder->quots, gdouble, 0);
    ladder->row_canc_prev     = g_array_index (ladder->cancs, gdouble, 0);
    ladder->row_canc_abs_prev = g_array_index (ladder->cancs_abs, gdouble, 0);
    ladder->control.df_best   = ladder->row_prev;
  }

  _ncm_diff_ladder_accum (ladder, dtable, ladder->order + 2);

  if (!ladder->converged && !_ncm_diff_control_update (&ladder->control, trunc_ratio, canc_pad,
                                                       ladder->row, ladder->row_prev,
                                                       ladder->row_canc, ladder->row_canc_prev,
                                                       ladder->row_canc_abs, ladder->row_canc_abs_prev,
                                                       ladder->quot_scale, ladder->lead_share))
    ladder->converged = TRUE;

  if (ladder->order == 0)
    ladder->first_informative = ladder->control.informative;

  ladder->row_prev          = ladder->row;
  ladder->row_canc_prev     = ladder->row_canc;
  ladder->row_canc_abs_prev = ladder->row_canc_abs;
  ladder->order++;

  if (ladder->order == ladder->tables->len)
    ladder->converged = TRUE;
}

/*
 * Richardson extrapolation at the order of dtable: df = sum_t lambda_t quots[t],
 * with the cancellation scales combined in quadrature. Also measures the row:
 * the largest quotient quot_scale = max_t |quots[t]|, the scale of what is being
 * extrapolated; the largest term max_term = max_t |lambda_t quots[t]|; and
 * rho = |df| / max_term, near 1 / max |lambda| for a nonzero derivative, far
 * below it when the terms cancel (a zero derivative), and zero when every
 * quotient is zero. Since sum_t lambda_t = 1, quotients that agree give
 * rho = 1 / max |lambda|, while rho = 1 means the largest term is the whole
 * row and the others vanished or cancelled, which a smooth error series does
 * not do; lead_share = (rho - 1 / max |lambda|) / (1 - 1 / max |lambda|)
 * places the row between the two, 0 and 1.
 */
static void
_ncm_diff_ladder_accum (NcmDiffLadder *ladder, NcmDiffTable *dtable, const guint nt)
{
  gdouble quot_scale = 0.0;
  gdouble max_term   = 0.0;
  gdouble lambda_max = 0.0;
  guint t;

  ladder->row          = 0.0;
  ladder->row_canc     = 0.0;
  ladder->row_canc_abs = 0.0;

  for (t = 0; t < nt; t++)
  {
    const gdouble lambda_t = ncm_vector_get (dtable->lambda, t);
    const gdouble quot_t   = g_array_index (ladder->quots, gdouble, t);
    const gdouble canc_t   = g_array_index (ladder->cancs, gdouble, t);
    const gdouble cabs_t   = g_array_index (ladder->cancs_abs, gdouble, t);

    /* fma: one rounding per term. */
    ladder->row          = fma (lambda_t, quot_t, ladder->row);
    ladder->row_canc     = hypot (ladder->row_canc, lambda_t * canc_t);
    ladder->row_canc_abs = hypot (ladder->row_canc_abs, lambda_t * cabs_t);

    quot_scale = GSL_MAX (quot_scale, fabs (quot_t));
    max_term   = GSL_MAX (max_term, fabs (lambda_t * quot_t));
    lambda_max = GSL_MAX (lambda_max, fabs (lambda_t));
  }

  ladder->quot_scale = quot_scale;
  ladder->max_term   = max_term;
  ladder->rho        = (max_term > 0.0) ? fabs (ladder->row) / max_term : 0.0;
  ladder->lead_share = (lambda_max > 1.0) ? (ladder->rho - 1.0 / lambda_max) / (1.0 - 1.0 / lambda_max) : 0.0;
}

/*
 * The agreement of two orders is measured against scale, the largest
 * difference quotient of the row: for a nonzero derivative that is the
 * derivative itself, for a zero one the size of the terms being cancelled.
 */
static gboolean
_ncm_diff_control_update (NcmDiffControl *cs, const gdouble trunc_ratio, const gdouble canc_pad,
                          const gdouble row, const gdouble row_prev,
                          const gdouble row_canc, const gdouble row_canc_prev,
                          const gdouble row_canc_abs, const gdouble row_canc_abs_prev,
                          const gdouble quot_scale, const gdouble lead_share)
{
  const gdouble trunc_est  = fabs (row - row_prev) * trunc_ratio;
  const gdouble rel_change = (quot_scale > 0.0) ? fabs (row - row_prev) / quot_scale : 0.0;

  /* The padded cancellation scale of the rounding plus the unpadded one of
   * the stated absolute precision: both come from the same values. */
  const gdouble canc_est_prev = fabs (row_canc_prev) * canc_pad + row_canc_abs_prev;
  const gdouble canc_est      = fabs (row_canc) * canc_pad + row_canc_abs;
  const gdouble err_curr_max  = GSL_MAX (trunc_est, GSL_MAX (canc_est_prev, canc_est));
  gdouble err_curr_best       = cs->err_best;
  gboolean improve            = FALSE;
  gboolean informative, eligible, agree;

  /*
   * A row is informative while its padded cancellation scale is below 1.0e-3
   * of quot_scale, the largest quotient of the row; rows above it (steps far
   * below the scale f varies on, quotients rounded to zero) can only agree by
   * accident. Two orders agree when they differ by less than 1.0e-3 of
   * quot_scale on an informative row. When every quotient is zero further
   * steps can only add cancellation, so the rows agree and the ladder stops
   * with zero and the padded cancellation scale as its error; such a row is
   * informative only when f is exactly zero on the samples, with a zero
   * cancellation scale. A row carried by its largest term (lead_share above
   * NCM_DIFF_LEAD_SHARE_MAX) is not a combination of quotients that agree,
   * whatever its cancellation scale says: such a row is not eligible to
   * agree, and never replaces a best value taken from eligible rows.
   */
  informative     = (quot_scale > 0.0) ? (canc_est < 1.0e-3 * quot_scale) : (canc_est == 0.0);
  eligible        = (lead_share < NCM_DIFF_LEAD_SHARE_MAX);
  agree           = (quot_scale > 0.0) ? ((rel_change < 1.0e-3) && informative && eligible) : TRUE;
  cs->informative = informative;

  /* The first rows fluctuate: convergence is checked only after
   * NCM_DIFF_NTRY_CONV consecutive agreeing rows. */
  if (cs->agreements_left && agree)
  {
    cs->agreements_left--;

    if (!cs->agreements_left)
      cs->just_converged = 1;
  }
  else if (!agree)
  {
    cs->agreements_left = NCM_DIFF_NTRY_CONV;
  }

  /* A converged row replaces the best value; its error is the mean of the
   * last two row errors, since a single row error fluctuates. */
  if (cs->just_converged || ((err_curr_max < cs->err_best) && !cs->agreements_left))
  {
    cs->df_best       = row;
    cs->err_best      = 0.5 * (err_curr_max + cs->err_last_max);
    cs->best_eligible = TRUE;

    err_curr_best = cs->err_best;
    improve       = TRUE;
  }
  else if ((err_curr_max < cs->err_best) && (eligible || !cs->best_eligible))
  {
    /* Not converged: keep the row with the smallest total error, eligible
     * rows first. */
    cs->df_best       = row;
    cs->err_best      = err_curr_max;
    cs->best_eligible = eligible;

    err_curr_best = err_curr_max;
  }

  /*
   * The cancellation scale of the rows grows with the order and bounds their total
   * error from below, so once it passes the best error no later row can rank
   * better. A ladder that has not converged goes on while its rows are
   * informative, since they may yet agree and replace a best value taken
   * from rows that agreed because the steps were aliased with an oscillation
   * of f.
   */
  if ((canc_est < err_curr_best) || (cs->agreements_left && informative))
    improve = TRUE;

  cs->just_converged = 0;
  cs->err_last_max   = err_curr_max;

  return improve;
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

/*
 * Replays the rows and the control over the stored steps after the top has
 * changed, up to convergence or the last stored step; the same as a fresh
 * run from the new top, without evaluating f again.
 */
static void
_ncm_diff_ladder_replay (NcmDiffLadder *ladder, const gdouble trunc_ratio, const gdouble canc_pad)
{
  _ncm_diff_ladder_restart (ladder);

  while (!ladder->converged && (ladder->order + 2 <= ladder->quots->len))
    _ncm_diff_ladder_extrapolate (ladder, trunc_ratio, canc_pad);
}

static gboolean _ncm_diff_dual_needs_step (NcmDiffDual *dual);
static void _ncm_diff_dual_extrapolate (NcmDiffDual *dual, const gdouble trunc_ratio, const gdouble canc_pad);
static gboolean _ncm_diff_dual_first_uninformative (NcmDiffDual *dual);
static void _ncm_diff_dual_replay (NcmDiffDual *dual, const gdouble trunc_ratio, const gdouble canc_pad);

static GArray *
_ncm_diff_by_step_algo_dual (NcmDiff *diff, NcmDiffStepAlgo step_algo, NcmDiffStepAlgo onesided_algo, guint po, GArray *x_a, const guint dim, NcmDiffFuncNtoM f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  GArray *duals               = g_array_new (FALSE, FALSE, sizeof (NcmDiffDual));
  NcmVector *x_v              = NULL;
  NcmVector *f_v              = NULL;
  NcmVector *yh1_v            = NULL;
  NcmVector *yh2_v            = NULL;
  NcmVector *quot_t[2]        = {NULL, NULL};
  NcmVector *canc_t[2]        = {NULL, NULL};
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
    quot_t[s] = ncm_vector_new (dim);
    canc_t[s] = ncm_vector_new (dim);
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
        const gint k = k_top + g_array_index (duals, NcmDiffDual, 0).A.quots->len;

        for (s = 0; s < 2; s++)
        {
          gdouble canc_abs_t;

          _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * _ncm_diff_step_h (self, sch.po, h0[s], k),
                               x_v, f_v, yh1_v, yh2_v, quot_t[s], canc_t[s], &canc_abs_t);

          for (i = 0; i < dim; i++)
          {
            NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

            _ncm_diff_ladder_add_step ((s == 0) ? &dual->A : &dual->B, ncm_vector_get (quot_t[s], i), ncm_vector_get (canc_t[s], i), canc_abs_t);
          }
        }

        continue;
      }

      for (i = 0; i < dim; i++)
      {
        NcmDiffDual *dual = &g_array_index (duals, NcmDiffDual, i);

        if (!dual->converged)
          _ncm_diff_dual_extrapolate (dual, self->trunc_ratio, self->canc_pad);
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
                     "below the initial step %.3e, and the quotients are dominated by cancellation.",
                     a, x, _ncm_diff_scheme_room (&sch), self->ini_h * scale);
        }
        else
        {
          gboolean finite = TRUE;
          gdouble canc_abs_t[2];

          for (s = 0; s < 2; s++)
          {
            _ncm_diff_eval_step (diff, sch.algo, f, user_data, a, x, sch.sign * _ncm_diff_step_h (self, sch.po, h0_new[s], k_new),
                                 x_v, f_v, yh1_v, yh2_v, quot_t[s], canc_t[s], &canc_abs_t[s]);
            finite = finite && _ncm_diff_vector_finite (quot_t[s]);
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
                  _ncm_diff_ladder_add_step (ladder, ncm_vector_get (quot_t[s], i), ncm_vector_get (canc_t[s], i), canc_abs_t[s]);
                else
                  _ncm_diff_ladder_prepend_step (ladder, ncm_vector_get (quot_t[s], i), ncm_vector_get (canc_t[s], i), canc_abs_t[s]);
              }

              if (!may_jump)
                _ncm_diff_dual_replay (dual, self->trunc_ratio, self->canc_pad);
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
        ncm_matrix_set (Eerr_m, a, i, _ncm_diff_cross_control_err (&dual->cross));
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
      ncm_vector_clear (&quot_t[s]);
      ncm_vector_clear (&canc_t[s]);
    }

    ncm_matrix_clear (&df_m);
    ncm_matrix_clear (&Eerr_m);

    return df;
  }
}

static gboolean
_ncm_diff_dual_needs_step (NcmDiffDual *dual)
{
  return !dual->converged && (dual->A.quots->len < dual->A.order + 2);
}

static gboolean _ncm_diff_cross_control_update (NcmDiffCrossControl *cs, const gdouble row_A, const gdouble row_B, const gdouble canc_A, const gdouble canc_B, const gdouble quot_scale);

/*
 * Extrapolates both ladders at the current order and lets the cross control
 * decide; the first NCM_DIFF_DUAL_MIN_ORDER orders are always taken.
 */
static void
_ncm_diff_dual_extrapolate (NcmDiffDual *dual, const gdouble trunc_ratio, const gdouble canc_pad)
{
  gboolean improve;

  _ncm_diff_ladder_extrapolate (&dual->A, trunc_ratio, canc_pad);
  _ncm_diff_ladder_extrapolate (&dual->B, trunc_ratio, canc_pad);

  improve = _ncm_diff_cross_control_update (&dual->cross, dual->A.row, dual->B.row,
                                            canc_pad * dual->A.row_canc + dual->A.row_canc_abs,
                                            canc_pad * dual->B.row_canc + dual->B.row_canc_abs,
                                            GSL_MAX (dual->A.quot_scale, dual->B.quot_scale));

  if (((dual->A.order >= NCM_DIFF_DUAL_MIN_ORDER) && !improve) || (dual->A.order == dual->A.tables->len))
    dual->converged = TRUE;
}

static gboolean
_ncm_diff_cross_control_update (NcmDiffCrossControl *cs, const gdouble row_A, const gdouble row_B,
                                const gdouble canc_A, const gdouble canc_B, const gdouble quot_scale)
{
  const gdouble cross = fabs (row_A - row_B);
  const gdouble canc  = GSL_MAX (fabs (canc_A), fabs (canc_B));
  const gdouble err   = GSL_MAX (cross, canc) * NCM_DIFF_DUAL_ERR_PAD;

  /*
   * The two ladders agree when they differ by less than 1.0e-3 of
   * quot_scale, the largest quotient of either row, with both cancellation
   * scales below that level, the gate of a single ladder. Only agreeing
   * rows are eligible: two ladders whose steps are not commensurate alias
   * an oscillation of f differently, so an aliased plateau of one ladder is
   * not one of the other. Before any agreement the row with the smallest
   * error relative to quot_scale is kept, since rows from steps far apart
   * can differ in size by orders of magnitude. Agreeing rows have the size
   * of the derivative and are ranked by their disagreement alone: it
   * measures the scatter of the values of f together with the truncation
   * error, while the cancellation scale assumes the stated precision, which
   * for a function more precise than stated would stop the ladders where its
   * noise has not yet appeared. The cancellation scale stays in the gate and
   * in the error reported. A disagreeing row never replaces an agreeing
   * one, and the first agreeing row replaces a disagreeing one whatever
   * their errors.
   */
  const gboolean eligible = (quot_scale > 0.0) ? ((cross < 1.0e-3 * quot_scale) && (canc < 1.0e-3 * quot_scale)) : (canc == 0.0);
  /* Rows of quotients that are all exactly zero rank among themselves by their errors. */
  const gdouble rel_err = err / GSL_MAX (quot_scale, GSL_DBL_MIN);
  const gboolean better = cs->best_eligible ? (eligible && (cross < cs->cross_best)) : (eligible || (rel_err < cs->rel_err_best));

  /* B (smaller steps) has the smaller truncation error; when cancellation
   * dominates the disagreement, the ladder with the smaller cancellation scale is
   * taken. */
  const gdouble row_sel = (fabs (canc_B) <= cross) ? row_B : ((fabs (canc_A) < fabs (canc_B)) ? row_A : row_B);

  cs->row_lo = GSL_MIN (cs->row_lo, row_sel);
  cs->row_hi = GSL_MAX (cs->row_hi, row_sel);

  if (eligible)
  {
    cs->elig_lo = GSL_MIN (cs->elig_lo, row_sel);
    cs->elig_hi = GSL_MAX (cs->elig_hi, row_sel);

    if (cs->elig_val != NULL)
    {
      g_array_append_val (cs->elig_val, row_sel);
      g_array_append_val (cs->elig_err, err);
    }
  }

  if (better)
  {
    cs->df_best = row_sel;
    cs->elig_lo = row_sel;
    cs->elig_hi = row_sel;

    cs->err_best      = err;
    cs->rel_err_best  = rel_err;
    cs->cross_best    = cross;
    cs->n_worse       = 0;
    cs->best_eligible = eligible;

    return TRUE;
  }

  /* While the cancellation scale is far below the best error the stagnation is
   * pre-asymptotic (steps larger than the scale f varies on): continue. Before
   * the ladders first agree the best error comes from rows that disagree and
   * says nothing about how far they are from agreeing, so the ladders continue
   * while the cancellation scale is below the agreement level. */
  if ((canc * NCM_DIFF_DUAL_ERR_PAD < cs->err_best) || (!cs->best_eligible && (canc < 1.0e-3 * quot_scale)))
  {
    cs->n_worse = 0;

    return TRUE;
  }

  cs->n_worse++;

  return cs->n_worse < NCM_DIFF_DUAL_MAX_WORSE;
}

/* The first row of either ladder is not informative. */
static gboolean
_ncm_diff_dual_first_uninformative (NcmDiffDual *dual)
{
  return !dual->A.first_informative || !dual->B.first_informative;
}

/* Replays both ladders and the cross control over the stored steps. */
static void
_ncm_diff_dual_replay (NcmDiffDual *dual, const gdouble trunc_ratio, const gdouble canc_pad)
{
  _ncm_diff_ladder_restart (&dual->A);
  _ncm_diff_ladder_restart (&dual->B);
  _ncm_diff_dual_cross_restart (dual);
  dual->converged = FALSE;

  while (!dual->converged && (dual->A.order + 2 <= dual->A.quots->len))
    _ncm_diff_dual_extrapolate (dual, trunc_ratio, canc_pad);
}

static GArray *_ncm_diff_Hessian_by_step_algo_single (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr);
static GArray *_ncm_diff_Hessian_by_step_algo_dual (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr);

static GArray *
ncm_diff_Hessian_by_step_algo (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, guint po, GArray *x_a, NcmDiffFuncNto1 f, gpointer user_data, GArray **Eerr)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  if (self->dual_series)
    return _ncm_diff_Hessian_by_step_algo_dual (diff, Hstep_algo, po, x_a, f, user_data, Eerr);

  return _ncm_diff_Hessian_by_step_algo_single (diff, Hstep_algo, po, x_a, f, user_data, Eerr);
}

static void _ncm_diff_hessian_eval_step (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hxo, const guint b, const gdouble y, const gdouble hyo, NcmVector *x_v, const gdouble fval, gdouble *quot_t, gdouble *canc_t, gdouble *canc_abs_t);
static gboolean _ncm_diff_hessian_move_fits (NcmDiffPrivate *self, NcmDiffScheme *sx, NcmDiffScheme *sy, const guint a, const gdouble x, const guint b, const gdouble y, const gboolean jump, const gint k_new, gdouble *hx0_new, gdouble *hy0_new, gboolean *restart);
static void _ncm_diff_hessian_restart_bases (NcmDiffPrivate *self, NcmDiffScheme *sx, NcmDiffScheme *sy, const gdouble x, const gdouble y);

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
          const gint k = k_top + ladder.quots->len;
          gdouble quot_t, canc_t, canc_abs_t;

          _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                       a, x, sx.sign * _ncm_diff_step_h (self, 0, sx.h0, k),
                                       b, y, sy.sign * _ncm_diff_step_h (self, 0, sy.h0, k),
                                       x_v, fval, &quot_t, &canc_t, &canc_abs_t);
          _ncm_diff_ladder_add_step (&ladder, quot_t, canc_t, canc_abs_t);
        }

        _ncm_diff_ladder_extrapolate (&ladder, self->trunc_ratio, self->canc_pad);

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
            gdouble quot_t, canc_t, canc_abs_t;

            _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                         a, x, sx.sign * _ncm_diff_step_h (self, 0, hx0_new, k_new),
                                         b, y, sy.sign * _ncm_diff_step_h (self, 0, hy0_new, k_new),
                                         x_v, fval, &quot_t, &canc_t, &canc_abs_t);

            if (!gsl_finite (quot_t))
            {
              may_move = FALSE;
            }
            else if (may_jump)
            {
              _ncm_diff_ladder_reset (&ladder);
              _ncm_diff_ladder_add_step (&ladder, quot_t, canc_t, canc_abs_t);
              sx.h0    = hx0_new;
              sy.h0    = hy0_new;
              may_jump = FALSE;
            }
            else
            {
              _ncm_diff_ladder_prepend_step (&ladder, quot_t, canc_t, canc_abs_t);
              _ncm_diff_ladder_replay (&ladder, self->trunc_ratio, self->canc_pad);
              k_top = k_new;
            }
          }
        }
      }

      ncm_matrix_set (df_m, a, b, ladder.control.df_best);
      ncm_matrix_set (df_m, b, a, ladder.control.df_best);

      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, ladder.control.err_best);
        ncm_matrix_set (Eerr_m, b, a, ladder.control.err_best);
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

/* Evaluates one Hessian cross-term step with offsets hxo and hyo: the
 * difference quotient quot_t and its cancellation scale canc_t. */
static void
_ncm_diff_hessian_eval_step (NcmDiff *diff, NcmDiffHessianStepAlgo Hstep_algo, NcmDiffFuncNto1 f, gpointer user_data,
                             const guint a, const gdouble x, const gdouble hxo,
                             const guint b, const gdouble y, const gdouble hyo,
                             NcmVector *x_v, const gdouble fval,
                             gdouble *quot_t, gdouble *canc_t, gdouble *canc_abs_t)
{
  volatile gdouble t_x = x + hxo;
  const gdouble hx     = t_x - x;
  volatile gdouble t_y = y + hyo;
  const gdouble hy     = t_y - y;

  Hstep_algo (diff, f, user_data, a, x, hx, b, y, hy, x_v, fval, quot_t, canc_t, canc_abs_t);

  ncm_vector_set (x_v, a, x);
  ncm_vector_set (x_v, b, y);
}

/*
 * Moving up along the two axes of a cross term: the first move sets each
 * axis to the smaller of the step a zero coordinate gets and its cap, a
 * further move goes one factor up on both. Returns whether the move fits the
 * caps. When it does not, an axis without room falls back to the other side
 * and the pair restarts; an axis with room on neither side below its initial
 * step is an error.
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
               "below the initial step %.3e, and the quotients are dominated by cancellation.",
               a, x, _ncm_diff_scheme_room (sx), self->ini_h * scale_x);

    if (blocked_y && (_ncm_diff_scheme_room (sy) < self->ini_h * scale_y))
      g_error ("NcmDiff: coordinate %u at % .15g: its domain leaves room for steps of at most %.3e, "
               "below the initial step %.3e, and the quotients are dominated by cancellation.",
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
          const gint k = k_top + dual.A.quots->len;
          guint s;

          for (s = 0; s < 2; s++)
          {
            const gdouble fac = (s == 0) ? 1.0 : 1.0 / srs;
            gdouble quot_t, canc_t, canc_abs_t;

            _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                         a, x, sx.sign * _ncm_diff_step_h (self, 0, sx.h0 * fac, k),
                                         b, y, sy.sign * _ncm_diff_step_h (self, 0, sy.h0 * fac, k),
                                         x_v, fval, &quot_t, &canc_t, &canc_abs_t);
            _ncm_diff_ladder_add_step ((s == 0) ? &dual.A : &dual.B, quot_t, canc_t, canc_abs_t);
          }
        }

        _ncm_diff_dual_extrapolate (&dual, self->trunc_ratio, self->canc_pad);

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
            gdouble quot_t[2], canc_t[2], canc_abs_t[2];
            gboolean finite = TRUE;
            guint s;

            for (s = 0; s < 2; s++)
            {
              const gdouble fac = (s == 0) ? 1.0 : 1.0 / srs;

              _ncm_diff_hessian_eval_step (diff, Hstep_algo, f, user_data,
                                           a, x, sx.sign * _ncm_diff_step_h (self, 0, hx0_new * fac, k_new),
                                           b, y, sy.sign * _ncm_diff_step_h (self, 0, hy0_new * fac, k_new),
                                           x_v, fval, &quot_t[s], &canc_t[s], &canc_abs_t[s]);
              finite = finite && gsl_finite (quot_t[s]);
            }

            if (!finite)
            {
              may_move = FALSE;
            }
            else if (may_jump)
            {
              _ncm_diff_dual_reset (&dual);
              _ncm_diff_ladder_add_step (&dual.A, quot_t[0], canc_t[0], canc_abs_t[0]);
              _ncm_diff_ladder_add_step (&dual.B, quot_t[1], canc_t[1], canc_abs_t[1]);
              sx.h0    = hx0_new;
              sy.h0    = hy0_new;
              may_jump = FALSE;
            }
            else
            {
              _ncm_diff_ladder_prepend_step (&dual.A, quot_t[0], canc_t[0], canc_abs_t[0]);
              _ncm_diff_ladder_prepend_step (&dual.B, quot_t[1], canc_t[1], canc_abs_t[1]);
              _ncm_diff_dual_replay (&dual, self->trunc_ratio, self->canc_pad);
              k_top = k_new;
            }
          }
        }
      }

      ncm_matrix_set (df_m, a, b, dual.cross.df_best);
      ncm_matrix_set (df_m, b, a, dual.cross.df_best);

      if (Eerr_m != NULL)
      {
        ncm_matrix_set (Eerr_m, a, b, _ncm_diff_cross_control_err (&dual.cross));
        ncm_matrix_set (Eerr_m, b, a, _ncm_diff_cross_control_err (&dual.cross));
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

static void _ncm_diff_step_quotient (NcmVector *f1, const NcmVector *f2, const gdouble scale, NcmVector *canc);

static void
_ncm_diff_rf_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot, NcmVector *canc, gdouble *canc_abs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  ncm_vector_addto (x_v, a, h);

  f (x_v, quot, user_data);

  _ncm_diff_step_quotient (quot, f_v, 1.0 / h, canc);

  /* Two values, each within func_abs_prec. */
  canc_abs[0] = 2.0 * self->func_abs_prec / fabs (h);

  NCM_UNUSED (yh1_v);
  NCM_UNUSED (yh2_v);
}

/*
 * Difference quotient (f1 - f2) scale, in place of f1, and its cancellation
 * scale eps max (|f1|, |f2|) |scale|. Subtracting close values loses digits:
 * when each value carries at most eps/2 of its size and they are within a
 * factor of 2, the difference is exact but keeps their errors, at most
 * eps max (|f1|, |f2|), so its relative accuracy is
 * eps max (|f1|, |f2|) / |f1 - f2|. Multiplied by the quotient
 * |f1 - f2| |scale| this gives the absolute error, in which |f1 - f2|
 * cancels; without it in the denominator the scale stays finite when the
 * difference is exactly zero. See "Cancellation and propagated evaluation
 * error" in docs/theory/ncm/integration/diff.qmd.
 */
static void
_ncm_diff_step_quotient (NcmVector *f1, const NcmVector *f2, const gdouble scale, NcmVector *canc)
{
  const guint len = ncm_vector_len (f1);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble f1_i = ncm_vector_get (f1, i);
    const gdouble f2_i = ncm_vector_get (f2, i);

    ncm_vector_set (f1,   i, (f1_i - f2_i) * scale);
    ncm_vector_set (canc, i, GSL_MAX (fabs (f1_i), fabs (f2_i)) * GSL_DBL_EPSILON * fabs (scale));
  }
}

static void
_ncm_diff_rc_d1_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot, NcmVector *canc, gdouble *canc_abs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);

  ncm_vector_addto (x_v, a, h);

  f (x_v, quot, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  _ncm_diff_step_quotient (quot, yh1_v, 0.5 / h, canc);

  /* Two values, each within func_abs_prec, over 2 h. */
  canc_abs[0] = self->func_abs_prec / fabs (h);

  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rc_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot, NcmVector *canc, gdouble *canc_abs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  const guint len             = ncm_vector_len (quot);
  const gdouble scale         = 2.0 / (h * h);
  guint i;

  ncm_vector_addto (x_v, a, h);

  f (x_v, quot, user_data);

  ncm_vector_set (x_v, a, x - h);

  f (x_v, yh1_v, user_data);

  /*
   * Cancellation error of the quotient: each of the three values carries
   * at most eps/2 of its size and the sum f (x + h) + f (x - h) one more
   * rounding, at most 3 eps max |f| in the numerator. It does not depend on
   * the difference, which can cancel exactly (an odd f at zero).
   */
  for (i = 0; i < len; i++)
  {
    const gdouble fp_i  = ncm_vector_get (quot, i);
    const gdouble fm_i  = ncm_vector_get (yh1_v, i);
    const gdouble f0_i  = ncm_vector_get (f_v, i);
    const gdouble avg_i = (fp_i + fm_i) * 0.5;
    const gdouble mag_i = GSL_MAX (GSL_MAX (fabs (fp_i), fabs (fm_i)), fabs (f0_i));

    ncm_vector_set (quot,   i, (avg_i - f0_i) * scale);
    ncm_vector_set (canc, i, 1.5 * mag_i * GSL_DBL_EPSILON * scale);
  }

  /* f (x + h) + f (x - h) - 2 f (x): at most 4 func_abs_prec, the quotient
   * scales (f (x + h) + f (x - h)) / 2 - f (x) by scale. */
  canc_abs[0] = 2.0 * self->func_abs_prec * scale;

  NCM_UNUSED (yh2_v);
}

/*
 * Forward second difference (f (x + 2h) - 2 f (x + h) + f (x)) / h^2, the
 * one-sided scheme rc_d2 falls back to at an edge of the domain; a negative
 * h makes it backward.
 */
static void
_ncm_diff_rf_d2_step (NcmDiff *diff, NcmDiffFuncNtoM f, gpointer user_data, const guint a, const gdouble x, const gdouble h, NcmVector *x_v, NcmVector *f_v, NcmVector *yh1_v, NcmVector *yh2_v, NcmVector *quot, NcmVector *canc, gdouble *canc_abs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  const guint len             = ncm_vector_len (quot);
  const gdouble scale         = 1.0 / (h * h);
  guint i;

  ncm_vector_addto (x_v, a, h);

  f (x_v, quot, user_data);

  ncm_vector_addto (x_v, a, h);

  f (x_v, yh1_v, user_data);

  for (i = 0; i < len; i++)
  {
    const gdouble f1_i  = ncm_vector_get (quot, i);
    const gdouble f2_i  = ncm_vector_get (yh1_v, i);
    const gdouble f0_i  = ncm_vector_get (f_v, i);
    const gdouble mag_i = GSL_MAX (GSL_MAX (fabs (f1_i), fabs (f2_i)), fabs (f0_i));

    ncm_vector_set (quot,   i, ((f2_i - f1_i) - (f1_i - f0_i)) * scale);
    ncm_vector_set (canc, i, 4.0 * mag_i * GSL_DBL_EPSILON * scale);
  }

  /* f (x + 2h) - 2 f (x + h) + f (x): at most 4 func_abs_prec. */
  canc_abs[0] = 4.0 * self->func_abs_prec * scale;

  NCM_UNUSED (yh2_v);
}

static void
_ncm_diff_rf_Hessian_step (NcmDiff *diff, NcmDiffFuncNto1 f, gpointer user_data, const guint a, const gdouble x, const gdouble hx, const guint b, const gdouble y, const gdouble hy, NcmVector *x_v, const gdouble fval, gdouble *quot, gdouble *canc, gdouble *canc_abs)
{
  NcmDiffPrivate * const self = ncm_diff_get_instance_private (diff);
  gdouble f_hx, f_hy, f_hxhy;


  ncm_vector_addto (x_v, a, hx);
  f_hx = f (x_v, user_data);

  ncm_vector_set (x_v, a, x);
  ncm_vector_addto (x_v, b, hy);
  f_hy = f (x_v, user_data);

  ncm_vector_addto (x_v, a, hx);
  f_hxhy = f (x_v, user_data);

  quot[0] = ((fval + f_hxhy) - (f_hx + f_hy)) / (hx * hy);

  {
    /* Cancellation error of quot[0]: the four values carry at most eps/2 of
     * their size each and the two sums one rounding each, at most
     * 4 eps max |f| in the numerator. */
    const gdouble max_f = GSL_MAX (GSL_MAX (fabs (fval), fabs (f_hxhy)), GSL_MAX (fabs (f_hx), fabs (f_hy)));

    canc[0] = 4.0 * max_f * GSL_DBL_EPSILON / fabs (hx * hy);
  }

  /* Four values, each within func_abs_prec. */
  canc_abs[0] = 4.0 * self->func_abs_prec / fabs (hx * hy);
}

/*
 * Spectral (Chebyshev) derivatives.
 *
 * The derivative along each variable is computed from a Chebyshev fit of the
 * function on a window [x - R, x + R], R at most spectral_window max (1, |x|).
 * Degree-four probes first grow R from ini_h |x| until their highest
 * coefficients exceed the rounding error of the values. The window
 * half-width R is then found by probing a dyadic ladder with a fixed-order
 * fit and scoring each by the estimated derivative error (series tail plus propagated
 * rounding error): shrinking resolves sharp features, expanding lowers the
 * amplification of the rounding error on flat ones. The accepted window is then refined
 * on nested Chebyshev-Lobatto grids.
 *
 * Convergence is judged on the derivative itself between refinement levels,
 * never on the coefficients: a relative test on the coefficient norm would be
 * dominated by components the derivative does not use (e.g. a large constant
 * offset), declaring convergence while the derivative-carrying coefficients
 * are still unresolved. The rounding error eps * sum |c_k|, propagated
 * through the coefficient differentiation, carries such offsets explicitly.
 */

#define NCM_DIFF_SC_PROBE_N (17) /* fixed order of the window probes */
#define NCM_DIFF_SC_MAX_N (65)   /* refinement grids: 17, 33, 65 nodes */
#define NCM_DIFF_SC_ERR_PAD (10.0)
#define NCM_DIFF_SC_MAX_SHRINK (45)
#define NCM_DIFF_SC_MAX_EXPAND (8)
#define NCM_DIFF_SC_BAD_Q (1.0e-3)        /* q_fit above this marks an unresolved window */
#define NCM_DIFF_SC_PATIENCE (2)          /* non-improving steps allowed after a resolved window */
#define NCM_DIFF_SC_IMPROVE (0.9)         /* required score reduction factor */
#define NCM_DIFF_SC_MAX_SCALE_STEPS (540) /* from GSL_DBL_MIN to the default cap 256 is 515 factors of 4 */

typedef struct _NcmDiffSCData
{
  NcmDiffFuncNtoM f;
  gpointer user_data;
  NcmVector *x_v;
  guint a;
  gdouble coeff_abs;
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
      const gdouble t  = (2 * i == N_new - 1) ? 0.0 : cos (M_PI * i / (N_new - 1.0));
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
    /* The center value is subtracted before the transform, so the rounding
     * of its products with the cosines does not enter every coefficient; it
     * is added back to c_0, where the error model keeps its rounding error. */
    const gdouble offset  = ncm_matrix_get (fvals, c, (N - 1) / 2);
    const gdouble f_first = ncm_matrix_get (fvals, c, 0) - offset;
    const gdouble f_last  = ncm_matrix_get (fvals, c, N - 1) - offset;

    for (k = 0; k < N; k++)
    {
      gdouble s = 0.5 * (f_first + (((k % 2) == 0) ? f_last : -f_last));

      for (i = 1; i < N - 1; i++)
        s += (ncm_matrix_get (fvals, c, i) - offset) * cosm[(k * i) % two_Nm1];

      ncm_matrix_set (coeffs, c, k, s * (((k == 0) || (k == N - 1)) ? 1.0 : 2.0) / (N - 1.0) + ((k == 0) ? offset : 0.0));
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
 *   coefficients, never on the differentiated ones: differentiation
 *   concentrates a slowly decaying series in its low coefficients, so its
 *   own tail is small while the fit has not resolved the function.
 * - coeff_round: the error nu = eps * sum |c_k| + coeff_abs of each
 *   coefficient, propagated the same way (the amplification of the k-th
 *   coefficient grows as k^2 per order). A value carrying eps/2 of its size
 *   moves each coefficient by at most eps max |f|, one carrying the absolute
 *   precision a by at most 2 a = coeff_abs.
 *
 * Also returns q_fit, the top-quarter coefficient mass above noise divided
 * by the nonconstant coefficient mass, with the offset retained in the noise.
 */
static void
_ncm_diff_sc_deriv_est (const gdouble *c, const guint len, const gdouble R, const guint order, const gdouble coeff_abs,
                        gdouble *deriv, gdouble *tail, gdouble *coeff_round, gdouble *q_fit)
{
  const gdouble Rinv     = 1.0 / R;
  const guint tail_start = (3 * len) / 4;
  gdouble *b             = g_new (gdouble, len);
  gdouble *nb            = g_new (gdouble, len);
  gdouble *tb            = g_new (gdouble, len);
  gdouble *tmp           = g_new (gdouble, len);
  gdouble sum_abs_c      = 0.0;
  gdouble sum_variation  = 0.0;
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

    if (k > 0)
      sum_variation += abs_ck;

    if (k >= tail_start)
      sum_abs_tail += abs_ck;
  }

  /* The denominator excludes c_0, so a constant offset cannot make an
   * unresolved tail small; only the coefficient noise is subtracted here,
   * the derivative error below keeps it. */
  q_fit[0] = GSL_MAX (0.0, sum_abs_tail - 8.0 * (len - tail_start) * (GSL_DBL_EPSILON * sum_abs_c + coeff_abs))
             / (sum_variation + GSL_DBL_MIN);

  {
    const gdouble nu = GSL_DBL_EPSILON * sum_abs_c + coeff_abs;

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
    gdouble tail_sum        = 0.0;
    gdouble coeff_round_sum = 0.0;

    for (k = 0; k < n; k++)
    {
      tail_sum        += tb[k];
      coeff_round_sum += nb[k];
    }

    tail[0]        = tail_sum * scale;
    coeff_round[0] = coeff_round_sum * scale;
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
 * - q: the worst q_fit over components, the fraction of nonconstant mass in
 *   the top quarter above its noise floor. A noise-only tail has q = 0;
 *   this alone says nothing about the accuracy of the derivative.
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
    gdouble deriv, tail, coeff_round, q_fit, err;

    _ncm_diff_sc_deriv_est (ncm_matrix_ptr (coeffs, c, 0), NCM_DIFF_SC_PROBE_N, R, order, data->coeff_abs, &deriv, &tail, &coeff_round, &q_fit);

    err = tail + coeff_round;

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
 * Outward search for the scale of f. Each probe is a degree-four fit on 5
 * Chebyshev-Lobatto nodes (3 nodes, then their 2 nested additions). The
 * half-width grows by 4 until c_3 and c_4 exceed both 16 times the
 * coefficient error nu = eps sum |c_k| + 2 a and 1e-3 of the nonconstant coefficients, or it
 * reaches cap. When the first probe is at the rounding level (nonconstant
 * coefficients below 1e3 nu in some component) the half-width is set to
 * R_jump at once, the probe a zero coordinate starts with. The result is
 * only the starting half-width of _ncm_diff_sc_scan_window(); a zero
 * derivative is never required to exceed the rounding error.
 */
static gdouble
_ncm_diff_sc_find_scale (NcmDiffSCData *data, const gdouble x, const gdouble R0, const gdouble R_jump, const gdouble cap,
                         const guint dim, NcmVector *y_v, NcmMatrix *fvals, NcmMatrix *coeffs)
{
  gdouble R       = GSL_MIN (R0, cap);
  gboolean jumped = FALSE;
  guint iter;

  for (iter = 0; iter < NCM_DIFF_SC_MAX_SCALE_STEPS; iter++)
  {
    gboolean feature = FALSE;
    guint c;

    if (!_ncm_diff_sc_eval_nodes (data, x, R, dim, y_v, fvals, 0, 3)
        || !_ncm_diff_sc_eval_nodes (data, x, R, dim, y_v, fvals, 3, 5))
      break;

    _ncm_diff_sc_dct (fvals, dim, 5, coeffs);

    {
      gboolean round_off = FALSE;

      for (c = 0; c < dim; c++)
      {
        gdouble variation  = 0.0;
        const gdouble tail = fabs (ncm_matrix_get (coeffs, c, 3)) + fabs (ncm_matrix_get (coeffs, c, 4));
        gdouble noise;
        guint k;

        for (k = 1; k < 5; k++)
          variation += fabs (ncm_matrix_get (coeffs, c, k));

        noise     = GSL_DBL_EPSILON * (fabs (ncm_matrix_get (coeffs, c, 0)) + variation) + data->coeff_abs;
        feature   = feature || (tail > GSL_MAX (16.0 * noise, NCM_DIFF_SC_BAD_Q * variation));
        round_off = round_off || (variation <= 1.0e3 * noise);
      }

      if (feature || (R >= cap))
        break;

      if (!jumped && round_off && (R < R_jump))
      {
        jumped = TRUE;
        R      = GSL_MIN (R_jump, cap);
        continue;
      }
    }

    jumped = TRUE;
    R      = GSL_MIN (4.0 * R, cap);
  }

  return R;
}

/*
 * Window search: hill descent on the probe score over a dyadic ladder of
 * half-widths, shrinking first and expanding only when no shrinking improved
 * on the initial window. While no resolved window (q below NCM_DIFF_SC_BAD_Q)
 * has been seen the shrinking continues: for a feature much narrower than
 * the window the score can move away from the eventual minimum until the
 * window reaches the feature scale. Expansion requires a resolved fit, since
 * it only lowers the amplification of the rounding error.
 * Stores the node values of the best window in fvals_best and returns its
 * half-width, or 0.0 when no finite window was found.
 */
static gdouble
_ncm_diff_sc_scan_window (NcmDiffSCData *data, const gdouble x, const gdouble R0, const gdouble cap, const guint dim,
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

    /* A resolved window replaces an unresolved best regardless of score: the
     * score of an unresolved fit measures its tail, not its distance to f. */
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
      if (R >= cap)
        break;

      R = GSL_MIN (2.0 * R, cap);

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
  GArray *conv                = g_array_new (FALSE, FALSE, sizeof (NcmDiffCrossControl));
  NcmDiffSCData data          = {f, user_data, x_v, 0, 2.0 * self->func_abs_prec};
  NcmMatrix *Eerr_m           = NULL;
  NcmMatrix *df_m;
  gboolean fallback = FALSE;
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
    const gdouble x      = g_array_index (x_a, gdouble, a);
    const gdouble scale  = (x == 0.0) ? 1.0 : fabs (x);
    const gdouble R0     = GSL_MAX (self->ini_h * scale, GSL_DBL_MIN);
    const gdouble R_jump = self->ini_h * GSL_MAX (1.0, scale);
    gdouble lo, hi, cap, seed;
    gdouble R;
    guint N_cur = 0;
    guint c;

    data.a = a;

    _ncm_diff_room (self, a, x, &lo, &hi);
    cap = GSL_MIN (GSL_MIN (lo, hi) * 0.5, self->spectral_window * GSL_MAX (1.0, scale));
    cap = GSL_MIN (cap, GSL_DBL_MAX * 0.25);

    if (!(cap > 0.0))
    {
      fallback = TRUE;
      break;
    }

    seed = _ncm_diff_sc_find_scale (&data, x, R0, R_jump, cap, dim, y_v, fvals_probe, coeffs_probe);
    R    = _ncm_diff_sc_scan_window (&data, x, seed, cap, dim, order, y_v, fvals_probe, coeffs_probe, fvals_best,
                                     &g_array_index (w_a, gdouble, 0));

    if (R == 0.0)
      g_error ("ncm_diff_sc: no window around x[%u] = % 22.15g with finite function values.", a, x);

    for (c = 0; c < dim; c++)
      _ncm_diff_cross_control_init (&g_array_index (conv, NcmDiffCrossControl, c));

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
          NcmDiffCrossControl *cs = &g_array_index (conv, NcmDiffCrossControl, c);
          gdouble d_c, tail_c, coeff_round_c, q_fit_c;

          _ncm_diff_sc_deriv_est (ncm_matrix_ptr (coeffs, c, 0), N, R, order, data.coeff_abs, &d_c, &tail_c, &coeff_round_c, &q_fit_c);

          if (tail_c > coeff_round_c)
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
            const gdouble err_c = GSL_MAX (fabs (d_c - g_array_index (d_prev, gdouble, c)), tail_c + coeff_round_c);

            if (err_c < cs->err_best)
            {
              cs->df_best  = d_c;
              cs->err_best = err_c;
              improve      = TRUE;
            }
          }

          g_array_index (d_prev, gdouble, c) = d_c;
        }

        /* Once every component's series tail is below its rounding error, a finer
         * grid can only add noise. Never stop before one cross-level check. */
        if ((N > NCM_DIFF_SC_PROBE_N) && (!improve || resolved))
          break;
      }
    }

    for (c = 0; c < dim; c++)
    {
      const NcmDiffCrossControl *cs = &g_array_index (conv, NcmDiffCrossControl, c);

      /* Near an edge no symmetric window inside the domain may give an
       * informative derivative; the call is then redone with the Richardson
       * driver, which can use a one-sided scheme. */
      if ((gsl_finite (lo) || gsl_finite (hi)) && (R >= cap * 0.5)
          && (cs->err_best > 1.0e-3 * fabs (cs->df_best)))
        fallback = TRUE;

      ncm_matrix_set (df_m, a, c, cs->df_best);

      if (Eerr_m != NULL)
        ncm_matrix_set (Eerr_m, a, c, cs->err_best * NCM_DIFF_SC_ERR_PAD);
    }

    if (fallback)
      break;
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

    if (fallback)
    {
      if (self->domain_warnings)
        g_warning ("NcmDiff: coordinate %u at % .15g has no symmetric spectral window inside its domain "
                   "giving an informative derivative, falling back to the Richardson central difference.",
                   a, g_array_index (x_a, gdouble, a));

      g_array_unref (df);

      if (Eerr != NULL)
        g_clear_pointer (Eerr, g_array_unref);

      return (order == 1) ? ncm_diff_rc_d1_N_to_M (diff, x_a, dim, f, user_data, Eerr)
                         : ncm_diff_rc_d2_N_to_M (diff, x_a, dim, f, user_data, Eerr);
    }

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
 * Calculates the first derivatives $\partial_i f_j$ of $f:\mathbb{R}^N\to \mathbb{R}^M$ at
 * @x_a using the forward method plus Richardson extrapolation, where $N$ is
 * the length of @x_a and $M = $ @dim. Element $i M + j$ of the result and of
 * @Eerr is the value and the estimated absolute error of $\partial_i f_j$.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the first derivatives $\partial_i f_j$ of $f:\mathbb{R}^N\to \mathbb{R}^M$ at
 * @x_a using the central method plus Richardson extrapolation, where $N$ is
 * the length of @x_a and $M = $ @dim. Element $i M + j$ of the result and of
 * @Eerr is the value and the estimated absolute error of $\partial_i f_j$.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the second derivatives $\partial_i^2 f_j$ of $f:\mathbb{R}^N\to \mathbb{R}^M$ at
 * @x_a using the central method plus Richardson extrapolation, where $N$ is
 * the length of @x_a and $M = $ @dim. Element $i M + j$ of the result and of
 * @Eerr is the value and the estimated absolute error of $\partial_i^2 f_j$.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the first derivatives of $f:\mathbb{R} \to \mathbb{R}^M$ at @x, $M = $
 * @dim, using the forward method plus Richardson extrapolation. Element $j$
 * of @Eerr is the estimated absolute error of component $j$.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x.
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
 * Calculates the first derivatives of $f:\mathbb{R} \to \mathbb{R}^M$ at @x, $M = $
 * @dim, using the central method plus Richardson extrapolation. Element $j$
 * of @Eerr is the estimated absolute error of component $j$.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x.
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
 * Calculates the second derivatives of $f:\mathbb{R} \to \mathbb{R}^M$ at @x, $M = $
 * @dim, using the central method plus Richardson extrapolation. Element $j$
 * of @Eerr is the estimated absolute error of component $j$.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x.
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
 * Calculates the first derivatives $\partial_i f$ of $f:\mathbb{R}^N \to \mathbb{R}$ at
 * @x_a using the forward method plus Richardson extrapolation, where $N$ is
 * the length of @x_a. Element $i$ of @Eerr is the estimated absolute error.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the first derivatives $\partial_i f$ of $f:\mathbb{R}^N \to \mathbb{R}$ at
 * @x_a using the central method plus Richardson extrapolation, where $N$ is
 * the length of @x_a. Element $i$ of @Eerr is the estimated absolute error.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the second derivatives $\partial_i^2 f$ of $f:\mathbb{R}^N \to \mathbb{R}$ at
 * @x_a using the central method plus Richardson extrapolation, where $N$ is
 * the length of @x_a. Element $i$ of @Eerr is the estimated absolute error.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the Hessian $\partial_i\partial_j f$ of $f:\mathbb{R}^N \to \mathbb{R}$ at
 * @x_a, where $N$ is the length of @x_a, plus Richardson extrapolation: the
 * diagonal by the central second difference, the entries $i \neq j$ by the
 * forward mixed difference with the steps of both coordinates. Element
 * $i N + j$ of the result and of @Eerr is the value and the estimated
 * absolute error of $\partial_i\partial_j f$.
 *
 * Returns: (transfer full) (array) (element-type double): the Hessian of @f at @x_a.
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
 * Calculates $f'(x)$ of $f:\mathbb{R} \to \mathbb{R}$ using the forward method plus
 * Richardson extrapolation; @err receives its estimated absolute error.
 *
 * Returns: the derivative of @f at @x.
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
 * Calculates $f'(x)$ of $f:\mathbb{R} \to \mathbb{R}$ using the central method plus
 * Richardson extrapolation; @err receives its estimated absolute error.
 *
 * Returns: the derivative of @f at @x.
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
 * Calculates $f''(x)$ of $f:\mathbb{R} \to \mathbb{R}$ using the central method plus
 * Richardson extrapolation; @err receives its estimated absolute error.
 *
 * Returns: the derivative of @f at @x.
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
 * Calculates the first derivatives $\partial_i f_j$ of
 * $f:\mathbb{R}^N\to \mathbb{R}^M$ at @x_a, where $N$ is the length of @x_a
 * and $M = $ @dim, using the spectral method: for each variable $f$ is
 * expanded in Chebyshev polynomials on a window $[x - R, x + R]$, and the
 * expansion is differentiated at the window center. The search for $R$
 * starts from #NcmDiff:ini-h $|x|$ (#NcmDiff:ini-h at $x = 0$), and $R$ is
 * at most #NcmDiff:spectral-window $\max(1, |x|)$. Element $i M + j$ of the
 * result and of @Eerr is the value and the estimated absolute error of
 * $\partial_i f_j$.
 *
 * The spectral method samples $f$ on a window matched to its scale of
 * variation instead of a shrinking neighborhood of $x$, with a smaller
 * amplification of the rounding error and more function evaluations than the
 * finite-difference methods.
 * With a domain set, the window stays inside it. When no symmetric window
 * inside the domain gives an informative derivative, the whole call is
 * redone with the Richardson central method, which can fall back to
 * one-sided differences, with a warning unless #NcmDiff:domain-warnings is
 * %FALSE.
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the second derivatives $\partial_i^2 f_j$ of
 * $f:\mathbb{R}^N\to \mathbb{R}^M$ at @x_a using the spectral method, with the
 * layout of ncm_diff_sc_d1_N_to_M().
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the first derivatives of $f:\mathbb{R}\to \mathbb{R}^M$ at @x,
 * $M = $ @dim, using the spectral method, see ncm_diff_sc_d1_N_to_M().
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x.
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
 * Calculates the second derivatives of $f:\mathbb{R}\to \mathbb{R}^M$ at @x,
 * $M = $ @dim, using the spectral method, see ncm_diff_sc_d1_N_to_M().
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x.
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
 * Calculates the gradient of $f:\mathbb{R}^N \to \mathbb{R}$ at @x_a using the
 * spectral method, see ncm_diff_sc_d1_N_to_M().
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates the second derivatives $\partial_i^2 f$ of
 * $f:\mathbb{R}^N \to \mathbb{R}$ at @x_a using the spectral method, see
 * ncm_diff_sc_d1_N_to_M().
 *
 * Returns: (transfer full) (array) (element-type double): the derivatives of @f at @x_a.
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
 * Calculates $f'(x)$ of $f:\mathbb{R} \to \mathbb{R}$ using the spectral method,
 * see ncm_diff_sc_d1_N_to_M(); @err receives its estimated absolute error.
 *
 * Returns: the derivative of @f at @x.
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
 * Calculates $f''(x)$ of $f:\mathbb{R} \to \mathbb{R}$ using the spectral method,
 * see ncm_diff_sc_d1_N_to_M(); @err receives its estimated absolute error.
 *
 * Returns: the derivative of @f at @x.
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


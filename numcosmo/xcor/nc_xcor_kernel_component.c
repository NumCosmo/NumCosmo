/***************************************************************************
 *            nc_xcor_kernel_component.c
 *
 *  Wed February 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_kernel_component.c
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
 * NcXcorKernelComponent:
 *
 * Abstract base class for kernel components in cross-correlation calculations.
 *
 * This class provides a framework for defining physical components of
 * cross-correlation kernels. Each component represents a distinct physical
 * contribution, such as galaxy number counts, magnification bias, or ISW effect.
 * Components can be combined to form multi-component kernels.
 *
 * Subclasses must implement:
 * - `eval_kernel`: evaluates K(k, xi) for the component
 * - `eval_prefactor`: evaluates any k and $\ell$-dependent prefactor
 *
 * Optionally, subclasses can implement:
 * - `get_limits`: returns valid integration ranges for xi and k
 *
 * The class provides automatic kernel analysis functionality that studies the behavior
 * of K*xi(k, y/k) to optimize integration strategies.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "xcor/nc_xcor_kernel_component.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcXcorKernelComponentPrivate
{
  NcmSpline *k_max_spline;       /* k_max(y) - k value that maximizes K*xi(k, y/k) */
  NcmSpline *K_max_spline;       /* K*xi_max(y) - maximum value of K*xi(k, y/k) */
  NcmSpline *k_epsilon_spline;   /* k_epsilon(y) - k where K*xi drops to epsilon*K*xi_max */
  gsl_min_fminimizer *minimizer; /* GSL minimizer for finding k_max */
  gsl_root_fsolver *root_solver; /* GSL root solver for finding k_epsilon */
  gdouble epsilon;
  gdouble sqrt_epsilon;
  guint ny;       /* Number of y points for analysis */
  guint max_iter; /* Maximum iterations for GSL solvers */
  gdouble tol;    /* Tolerance for GSL solvers */
} NcXcorKernelComponentPrivate;

enum
{
  PROP_0,
  PROP_EPSILON,
  PROP_NY,
  PROP_MAX_ITER,
  PROP_TOL,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorKernelComponent, nc_xcor_kernel_component, G_TYPE_OBJECT)

static void
nc_xcor_kernel_component_init (NcXcorKernelComponent *comp)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  self->k_max_spline     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->K_max_spline     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->k_epsilon_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  self->minimizer        = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  self->root_solver      = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  self->epsilon          = 0.0;
  self->sqrt_epsilon     = 0.0;
  self->ny               = 0.0;
  self->max_iter         = 0.0;
  self->tol              = 0.0;
}

static void
nc_xcor_kernel_component_dispose (GObject *object)
{
  NcXcorKernelComponent *comp        = NC_XCOR_KERNEL_COMPONENT (object);
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  ncm_spline_clear (&self->k_max_spline);
  ncm_spline_clear (&self->K_max_spline);
  ncm_spline_clear (&self->k_epsilon_spline);

  /* Chain up */
  G_OBJECT_CLASS (nc_xcor_kernel_component_parent_class)->dispose (object);
}

static void
nc_xcor_kernel_component_finalize (GObject *object)
{
  NcXcorKernelComponent *comp        = NC_XCOR_KERNEL_COMPONENT (object);
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  /* Free GSL solvers */
  if (self->minimizer != NULL)
    gsl_min_fminimizer_free (self->minimizer);

  if (self->root_solver != NULL)
    gsl_root_fsolver_free (self->root_solver);

  /* Chain up */
  G_OBJECT_CLASS (nc_xcor_kernel_component_parent_class)->finalize (object);
}

static void
nc_xcor_kernel_component_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelComponent *comp = NC_XCOR_KERNEL_COMPONENT (object);

  switch (prop_id)
  {
    case PROP_EPSILON:
      nc_xcor_kernel_component_set_epsilon (comp, g_value_get_double (value));
      break;
    case PROP_NY:
      nc_xcor_kernel_component_set_ny (comp, g_value_get_uint (value));
      break;
    case PROP_MAX_ITER:
      nc_xcor_kernel_component_set_max_iter (comp, g_value_get_uint (value));
      break;
    case PROP_TOL:
      nc_xcor_kernel_component_set_tol (comp, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_xcor_kernel_component_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelComponent *comp = NC_XCOR_KERNEL_COMPONENT (object);

  switch (prop_id)
  {
    case PROP_EPSILON:
      g_value_set_double (value, nc_xcor_kernel_component_get_epsilon (comp));
      break;
    case PROP_NY:
      g_value_set_uint (value, nc_xcor_kernel_component_get_ny (comp));
      break;
    case PROP_MAX_ITER:
      g_value_set_uint (value, nc_xcor_kernel_component_get_max_iter (comp));
      break;
    case PROP_TOL:
      g_value_set_double (value, nc_xcor_kernel_component_get_tol (comp));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_xcor_kernel_component_class_init (NcXcorKernelComponentClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = nc_xcor_kernel_component_dispose;
  object_class->finalize     = nc_xcor_kernel_component_finalize;
  object_class->set_property = nc_xcor_kernel_component_set_property;
  object_class->get_property = nc_xcor_kernel_component_get_property;

  klass->eval_kernel    = NULL;
  klass->eval_prefactor = NULL;
  klass->get_limits     = NULL;

  /**
   * NcXcorKernelComponent:epsilon:
   *
   * The epsilon value for kernel analysis, determining where K*xi(k, y/k)
   * drops to epsilon * K*xi_max.
   */
  g_object_class_install_property (object_class,
                                   PROP_EPSILON,
                                   g_param_spec_double ("epsilon",
                                                        NULL,
                                                        "Epsilon value for kernel analysis",
                                                        0.0, 1.0, NC_XCOR_KERNEL_COMPONENT_DEFAULT_EPSILON,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelComponent:ny:
   *
   * Number of y points for kernel analysis.
   */
  g_object_class_install_property (object_class,
                                   PROP_NY,
                                   g_param_spec_uint ("ny",
                                                      NULL,
                                                      "Number of y points",
                                                      1, G_MAXUINT, 600,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelComponent:max-iter:
   *
   * Maximum iterations for GSL solvers.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ITER,
                                   g_param_spec_uint ("max-iter",
                                                      NULL,
                                                      "Maximum iterations for GSL solvers",
                                                      1, G_MAXUINT, 100,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelComponent:tol:
   *
   * Tolerance for GSL solvers.
   */
  g_object_class_install_property (object_class,
                                   PROP_TOL,
                                   g_param_spec_double ("tol",
                                                        NULL,
                                                        "Tolerance for GSL solvers",
                                                        0.0, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_kernel_component_ref:
 * @comp: a #NcXcorKernelComponent
 *
 * Increases the reference count of @comp by one.
 *
 * Returns: (transfer full): @comp
 */
NcXcorKernelComponent *
nc_xcor_kernel_component_ref (NcXcorKernelComponent *comp)
{
  return g_object_ref (comp);
}

/**
 * nc_xcor_kernel_component_free:
 * @comp: a #NcXcorKernelComponent
 *
 * Decreases the reference count of @comp by one. If the reference count
 * reaches zero, the object is freed.
 */
void
nc_xcor_kernel_component_free (NcXcorKernelComponent *comp)
{
  g_object_unref (comp);
}

/**
 * nc_xcor_kernel_component_clear:
 * @comp: a #NcXcorKernelComponent
 *
 * Decreases the reference count of *@comp by one and sets *@comp to NULL.
 */
void
nc_xcor_kernel_component_clear (NcXcorKernelComponent **comp)
{
  g_clear_object (comp);
}

/**
 * nc_xcor_kernel_component_set_epsilon:
 * @comp: a #NcXcorKernelComponent
 * @epsilon: the epsilon value for kernel analysis
 *
 * Sets the epsilon value used in kernel analysis to determine where K*xi(k, y/k) drops
 * to epsilon * K*xi_max.
 */
void
nc_xcor_kernel_component_set_epsilon (NcXcorKernelComponent *comp, gdouble epsilon)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  self->epsilon      = epsilon;
  self->sqrt_epsilon = sqrt (epsilon);
}

/**
 * nc_xcor_kernel_component_get_epsilon:
 * @comp: a #NcXcorKernelComponent
 *
 * Gets the epsilon value used in kernel analysis.
 *
 * Returns: the epsilon value
 */
gdouble
nc_xcor_kernel_component_get_epsilon (NcXcorKernelComponent *comp)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  return self->epsilon;
}

/**
 * nc_xcor_kernel_component_set_ny:
 * @comp: a #NcXcorKernelComponent
 * @ny: number of y points for analysis
 *
 * Sets the number of y points to use in kernel analysis.
 */
void
nc_xcor_kernel_component_set_ny (NcXcorKernelComponent *comp, guint ny)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  self->ny = ny;
}

/**
 * nc_xcor_kernel_component_get_ny:
 * @comp: a #NcXcorKernelComponent
 *
 * Gets the number of y points used in kernel analysis.
 *
 * Returns: the number of y points
 */
guint
nc_xcor_kernel_component_get_ny (NcXcorKernelComponent *comp)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  return self->ny;
}

/**
 * nc_xcor_kernel_component_set_max_iter:
 * @comp: a #NcXcorKernelComponent
 * @max_iter: maximum iterations for GSL solvers
 *
 * Sets the maximum number of iterations for GSL minimizer and root finder.
 */
void
nc_xcor_kernel_component_set_max_iter (NcXcorKernelComponent *comp, guint max_iter)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  self->max_iter = max_iter;
}

/**
 * nc_xcor_kernel_component_get_max_iter:
 * @comp: a #NcXcorKernelComponent
 *
 * Gets the maximum number of iterations for GSL solvers.
 *
 * Returns: the maximum iterations
 */
guint
nc_xcor_kernel_component_get_max_iter (NcXcorKernelComponent *comp)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  return self->max_iter;
}

/**
 * nc_xcor_kernel_component_set_tol:
 * @comp: a #NcXcorKernelComponent
 * @tol: tolerance for GSL solvers
 *
 * Sets the tolerance for GSL minimizer and root finder.
 */
void
nc_xcor_kernel_component_set_tol (NcXcorKernelComponent *comp, gdouble tol)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  self->tol = tol;
}

/**
 * nc_xcor_kernel_component_get_tol:
 * @comp: a #NcXcorKernelComponent
 *
 * Gets the tolerance for GSL solvers.
 *
 * Returns: the tolerance value
 */
gdouble
nc_xcor_kernel_component_get_tol (NcXcorKernelComponent *comp)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  return self->tol;
}

/**
 * nc_xcor_kernel_component_eval_k_max:
 * @comp: a #NcXcorKernelComponent
 * @y: the y value (y = k * xi)
 *
 * Evaluates k_max at the given y value from kernel analysis, where k_max
 * is the value of k that maximizes K*xi(k, y/k) for this y.
 *
 * Returns: the k_max value at y
 */
gdouble
nc_xcor_kernel_component_eval_k_max (NcXcorKernelComponent *comp, gdouble y)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  g_assert (self->k_max_spline != NULL);

  return ncm_spline_eval (self->k_max_spline, y);
}

/**
 * nc_xcor_kernel_component_eval_K_max:
 * @comp: a #NcXcorKernelComponent
 * @y: the y value (y = k * xi)
 *
 * Evaluates the maximum value of K*xi(k, y/k) at the given y value from kernel analysis.
 * This is the value of K*xi at k = k_max(y).
 *
 * Returns: the K*xi_max value at y
 */
gdouble
nc_xcor_kernel_component_eval_K_max (NcXcorKernelComponent *comp, gdouble y)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  g_assert (self->K_max_spline != NULL);

  return ncm_spline_eval (self->K_max_spline, y);
}

/**
 * nc_xcor_kernel_component_eval_k_epsilon:
 * @comp: a #NcXcorKernelComponent
 * @y: the y value (y = k * xi)
 *
 * Evaluates k_epsilon at the given y value from kernel analysis, where k_epsilon
 * is the value of k (beyond k_max) where K*xi(k, y/k) drops to epsilon times K*xi_max.
 *
 * Returns: the k_epsilon value at y
 */
gdouble
nc_xcor_kernel_component_eval_k_epsilon (NcXcorKernelComponent *comp, gdouble y)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);

  g_assert (self->k_epsilon_spline != NULL);

  return ncm_spline_eval (self->k_epsilon_spline, y);
}

/* Structure for GSL function evaluation */
typedef struct _NcXcorKernelAnalysisData
{
  NcXcorKernelComponent *comp;
  NcHICosmo *cosmo;
  gdouble y;
  gdouble xi_min;
  gdouble xi_max;
  gdouble Kxi_threshold;
} NcXcorKernelAnalysisData;

/* GSL function: Returns -K*xi(k, y/k) for minimization (we want maximum, so negate) */
static gdouble
_nc_xcor_kernel_component_minus_Kxi (gdouble k, void *params)
{
  NcXcorKernelAnalysisData *data = (NcXcorKernelAnalysisData *) params;
  const gdouble xi               = data->y / k;
  const gdouble K_val            = nc_xcor_kernel_component_eval_kernel (data->comp, data->cosmo, xi, k);
  const gdouble Kxi              = xi * K_val;

  return -fabs (Kxi);
}

static gdouble
_nc_xcor_kernel_component_Kxi_minus_threshold (gdouble k, void *params)
{
  NcXcorKernelAnalysisData *data = (NcXcorKernelAnalysisData *) params;
  const gdouble xi               = data->y / k;
  const gdouble K_val            = nc_xcor_kernel_component_eval_kernel (data->comp, data->cosmo, xi, k);
  const gdouble Kxi              = xi * K_val;

  return log (fabs (Kxi / data->Kxi_threshold));
}

static void
_nc_xcor_kernel_component_find_k_max (NcXcorKernelComponent    *comp,
                                      NcXcorKernelAnalysisData *data,
                                      gdouble                  k_valid_min,
                                      gdouble                  k_valid_max,
                                      gdouble                  k_guess,
                                      gdouble                  *k_at_max,
                                      gdouble                  *K_max)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);
  gdouble k_lower                    = k_valid_min;
  gdouble k_upper                    = k_valid_max;
  gdouble k_init                     = k_guess;
  guint iter                         = 0;
  gdouble lk                         = 0.0;
  gdouble lk_lower                   = 0.0;
  gdouble lk_upper                   = 0.0;
  gsl_function F_min;
  gint status;

  F_min.function = &_nc_xcor_kernel_component_minus_Kxi;
  F_min.params   = data;

  if (k_lower >= k_upper)
  {
    if (fabs (k_lower / k_upper - 1.0) < GSL_DBL_EPSILON)
    {
      *k_at_max = k_lower;
      *K_max    = -_nc_xcor_kernel_component_minus_Kxi (*k_at_max, data);

      return;
    }
    else
    {
      g_error ("Invalid k range for finding k_max: [% 22.15g, % 22.15g]", k_lower, k_upper);

      return;
    }
  }

  if ((k_init < k_lower) || (k_init > k_upper))
    k_init = sqrt (k_lower * k_upper);

  /*
   * First we go through k with 50 steps in log space to find a good starting point for
   * the minimizer.
   */
  {
    const gdouble k_log_min  = log10 (k_lower);
    const gdouble k_log_max  = log10 (k_upper);
    const guint nsteps       = 50;
    const gdouble k_log_step = (k_log_max - k_log_min) / (nsteps - 1.0);
    gdouble k_best           = 0.0;
    gdouble K_best           = -G_MAXDOUBLE;
    guint best_idx           = 0;

    for (guint i = 0; i < nsteps; i++)
    {
      const gdouble k_test   = pow (10.0, k_log_min + i * k_log_step);
      const gdouble Kxi_test = -_nc_xcor_kernel_component_minus_Kxi (k_test, data);

      if (Kxi_test > K_best)
      {
        K_best   = Kxi_test;
        k_best   = k_test;
        best_idx = i;
      }
    }

    if (k_best > 0.0)
      k_init = k_best;

    if (best_idx == 0)
    {
      *k_at_max = k_lower;
      *K_max    = K_best;

      return;
    }
    else if (best_idx == nsteps - 1)
    {
      *k_at_max = k_upper;
      *K_max    = K_best;

      return;
    }
    else
    {
      k_lower = pow (10.0, k_log_min + (best_idx - 1) * k_log_step);
      k_upper = pow (10.0, k_log_min + (best_idx + 1) * k_log_step);
    }
  }

  gsl_min_fminimizer_set (self->minimizer, &F_min, k_init, k_lower, k_upper);

  do {
    iter++;
    status    = gsl_min_fminimizer_iterate (self->minimizer);
    *k_at_max = gsl_min_fminimizer_x_minimum (self->minimizer);
    k_lower   = gsl_min_fminimizer_x_lower (self->minimizer);
    k_upper   = gsl_min_fminimizer_x_upper (self->minimizer);
    status    = gsl_min_test_interval (k_lower, k_upper, 0.0, self->tol);

    if ((*k_at_max == lk) && (k_lower == lk_lower) && (k_upper == lk_upper))
    {
      const gdouble f_lower = gsl_min_fminimizer_f_lower (self->minimizer);
      const gdouble f_upper = gsl_min_fminimizer_f_upper (self->minimizer);

      if (fabs (2.0 * (f_upper - f_lower) / (f_upper + f_lower)) < self->tol)
      {
        status = GSL_SUCCESS;
      }
      else
      {
        g_error ("_nc_xcor_kernel_component_find_k_max: minimizer stuck at k = % 22.15g, y = % 22.15g",
                 *k_at_max, data->y);
        break;
      }
    }

    lk       = *k_at_max;
    lk_lower = k_lower;
    lk_upper = k_upper;
  } while ((status == GSL_CONTINUE) && (iter < self->max_iter));

  *K_max = -gsl_min_fminimizer_f_minimum (self->minimizer);
}

static gdouble
_nc_xcor_kernel_component_find_k_epsilon_high (NcXcorKernelComponent    *comp,
                                               NcXcorKernelAnalysisData *data,
                                               gdouble                  k_at_max,
                                               gdouble                  k_valid_max)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);
  gdouble k_low                      = k_at_max;
  gdouble k_high                     = k_valid_max;
  const gdouble f_low                = _nc_xcor_kernel_component_Kxi_minus_threshold (k_low, data);
  const gdouble f_high               = _nc_xcor_kernel_component_Kxi_minus_threshold (k_high, data);
  gsl_function F_root;


  F_root.function = &_nc_xcor_kernel_component_Kxi_minus_threshold;
  F_root.params   = data;

  if ((f_low > 0.0) && (f_high < 0.0))
  {
    guint iter = 0;
    gdouble k_epsilon;
    gint status;

    gsl_root_fsolver_set (self->root_solver, &F_root, k_low, k_high);

    do {
      iter++;
      status    = gsl_root_fsolver_iterate (self->root_solver);
      k_epsilon = gsl_root_fsolver_root (self->root_solver);
      k_low     = gsl_root_fsolver_x_lower (self->root_solver);
      k_high    = gsl_root_fsolver_x_upper (self->root_solver);
      status    = gsl_root_test_interval (k_low, k_high, 0.0, self->tol);
    } while (status == GSL_CONTINUE && iter < self->max_iter);

    return k_epsilon;
  }
  else
  {
    return k_high;
  }
}

/**
 * nc_xcor_kernel_component_prepare:
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 *
 * Prepares the kernel component by analyzing its behavior over the valid ranges. This
 * method calls get_limits to obtain the integration ranges, then studies K*xi(k, y/k)
 * to compute k_max(y), K_max(y), and k_epsilon(y) using GSL Brent minimizer and root
 * finder with warm starts.
 */
void
nc_xcor_kernel_component_prepare (NcXcorKernelComponent *comp, NcHICosmo *cosmo)
{
  NcXcorKernelComponentPrivate *self = nc_xcor_kernel_component_get_instance_private (comp);
  NcXcorKernelComponentClass *klass  = NC_XCOR_KERNEL_COMPONENT_GET_CLASS (comp);
  gdouble xi_min = 0.0, xi_max = 0.0, k_min = 0.0, k_max = 0.0;

  klass->get_limits (comp, cosmo, &xi_min, &xi_max, &k_min, &k_max);

  {
    NcmVector *yv                 = ncm_vector_new (self->ny);
    NcmVector *k_max_v            = ncm_vector_new (self->ny);
    NcmVector *K_max_v            = ncm_vector_new (self->ny);
    NcmVector *k_epsilon_v        = ncm_vector_new (self->ny);
    const gdouble y_min           = GSL_MAX (k_min * xi_min, 0.5);
    const gdouble y_max           = GSL_MIN (k_max * xi_max, 1000.0);
    NcXcorKernelAnalysisData data = {
      .comp          = comp,
      .cosmo         = cosmo,
      .y             = 0.0,
      .xi_min        = xi_min,
      .xi_max        = xi_max,
      .Kxi_threshold = 0.0
    };
    gdouble k_guess = 0.0;
    guint i;

    for (i = 0; i < self->ny; i++)
    {
      const gdouble log_y         = log (y_min) + (log (y_max) - log (y_min)) * i / (self->ny - 1.0);
      const gdouble y             = exp (log_y);
      const gdouble k_from_xi_min = y / xi_max;
      const gdouble k_from_xi_max = y / xi_min;
      const gdouble k_valid_min   = GSL_MAX (k_min, k_from_xi_min);
      const gdouble k_valid_max   = GSL_MIN (k_max, k_from_xi_max);
      gdouble k_at_max, K_max;

      data.y = y;
      ncm_vector_set (yv, i, y);

      if (k_valid_min >= k_valid_max)
      {
        g_warning ("# Skipping y = % 22.15g: no valid k range [% 22.15g, % 22.15g]\n", y, k_valid_min, k_valid_max);
        ncm_vector_set (k_max_v, i, 0.5 * (k_valid_min + k_valid_max));
        ncm_vector_set (K_max_v, i, 0.0);
        ncm_vector_set (k_epsilon_v, i, 0.5 * (k_valid_min + k_valid_max));

        continue;
      }

      _nc_xcor_kernel_component_find_k_max (comp, &data, k_valid_min, k_valid_max,
                                            k_guess, &k_at_max, &K_max);
      k_guess            = k_at_max;
      data.Kxi_threshold = self->sqrt_epsilon * K_max;

      {
        const gdouble k_epsilon = _nc_xcor_kernel_component_find_k_epsilon_high (comp, &data, k_at_max, k_valid_max);

        ncm_vector_set (k_max_v, i, k_at_max);
        ncm_vector_set (K_max_v, i, K_max);
        ncm_vector_set (k_epsilon_v, i, k_epsilon);
      }
    }

    ncm_spline_set (self->k_max_spline, yv, k_max_v, TRUE);
    ncm_spline_set (self->K_max_spline, yv, K_max_v, TRUE);
    ncm_spline_set (self->k_epsilon_spline, yv, k_epsilon_v, TRUE);

    ncm_vector_free (yv);
    ncm_vector_free (k_max_v);
    ncm_vector_free (K_max_v);
    ncm_vector_free (k_epsilon_v);
  }
}

/**
 * nc_xcor_kernel_component_eval_kernel: (virtual eval_kernel)
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 * @xi: comoving distance
 * @k: wave number
 *
 * Evaluates the kernel function K(k, xi) for this component.
 *
 * Returns: the value of K(k, xi)
 */
gdouble
nc_xcor_kernel_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k)
{
  NcXcorKernelComponentClass *klass = NC_XCOR_KERNEL_COMPONENT_GET_CLASS (comp);

  return klass->eval_kernel (comp, cosmo, xi, k);
}

/**
 * nc_xcor_kernel_component_eval_prefactor: (virtual eval_prefactor)
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 * @k: wave number
 * @l: multipole
 *
 * Evaluates the prefactor that may depend on k and ell.
 *
 * Returns: the prefactor value
 */
gdouble
nc_xcor_kernel_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  NcXcorKernelComponentClass *klass = NC_XCOR_KERNEL_COMPONENT_GET_CLASS (comp);

  g_assert (klass->eval_prefactor != NULL);

  return klass->eval_prefactor (comp, cosmo, k, l);
}

/**
 * nc_xcor_kernel_component_get_limits: (virtual get_limits)
 * @comp: a #NcXcorKernelComponent
 * @cosmo: a #NcHICosmo
 * @xi_min: (out): minimum comoving distance
 * @xi_max: (out): maximum comoving distance
 * @k_min: (out): minimum wave number
 * @k_max: (out): maximum wave number
 *
 * Gets the valid integration ranges for this component.
 */
void
nc_xcor_kernel_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max)
{
  NcXcorKernelComponentClass *klass = NC_XCOR_KERNEL_COMPONENT_GET_CLASS (comp);

  g_assert (klass->get_limits != NULL);

  klass->get_limits (comp, cosmo, xi_min, xi_max, k_min, k_max);
}


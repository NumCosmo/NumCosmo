/***************************************************************************
 *            nc_cluster_richness_projection.c
 *
 *  Wed September 17 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_cluster_richness_projection.c
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
 * NcClusterRichnessProjection:
 *
 * Projection contribution to the observed richness distribution.
 *
 * Computes the density of $\lambda = \lambda_\mathrm{true} + \Delta$, where
 * $\lambda_\mathrm{true}$ is log-normal with parameters $(\mu, \sigma)$ and $\Delta$
 * is an independent exponential of rate $\tau$,
 * \begin{equation}
 * T(\lambda) = \tau \int_0^\lambda f_\mathrm{LN}(t \mid \mu, \sigma)
 *              \, e^{-\tau (\lambda - t)} \, \mathrm{d}t .
 * \end{equation}
 * It is a normalized density on its own and carries no mixture weight; the caller
 * combines it with the unprojected term.
 *
 * The convolution is evaluated as an ODE in $x = \ln\lambda$ using a
 * #NcmOdeSpline, so a single prepare() serves every $\lambda$, and its integral
 * follows from a closed-form identity rather than a second solve. See
 * <a href="../../theory/cluster_richness_projection.html">Richness projection</a>
 * for the derivation.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/cluster/nc_cluster_richness_projection.h"
#include "ncm/core/ncm_c.h"
#include "ncm/spline/ncm_ode_spline.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Start of the integration, in units of sigma below mu. The log-normal CDF there,
 * 1.0e-23, bounds both T and the error of the initial value below. */
#define _NC_CLUSTER_RICHNESS_PROJECTION_NSIGMA (10.0)

/* Minimum knot density of the solver output. T is f_LN smoothed by a kernel of
 * width 1/tau, so in x it never varies faster than sigma. */
#define _NC_CLUSTER_RICHNESS_PROJECTION_KNOTS_PER_SIGMA (32.0)

typedef struct _NcClusterRichnessProjectionPrivate
{
  gdouble lnlambda_min;
  gdouble lnlambda_max;
  gdouble reltol;
  gdouble mu;
  gdouble sigma;
  gdouble tau;
  gdouble xi;
  NcmOdeSpline *T_ode;
  NcmSpline *T;
  gboolean prepared;
} NcClusterRichnessProjectionPrivate;

struct _NcClusterRichnessProjection
{
  GObject parent_instance;
};

enum
{
  PROP_0,
  PROP_LNLAMBDA_MIN,
  PROP_LNLAMBDA_MAX,
  PROP_RELTOL,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterRichnessProjection, nc_cluster_richness_projection, G_TYPE_OBJECT)

/*
 * dT/dx = tau * (g(x) - e^x T), with x = ln(lambda) and
 * g(x) = lambda * f_LN(lambda) the log-normal density in ln(lambda).
 */
static gdouble
_nc_cluster_richness_projection_dTdx (gdouble T, gdouble x, gpointer userdata)
{
  NcClusterRichnessProjectionPrivate * const self = (NcClusterRichnessProjectionPrivate *) userdata;
  const gdouble u                                 = (x - self->mu) / self->sigma;
  const gdouble g                                 = exp (-0.5 * u * u) / (ncm_c_sqrt_2pi () * self->sigma);

  return self->tau * (g - exp (x) * T);
}

static void
nc_cluster_richness_projection_init (NcClusterRichnessProjection *crp)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);
  NcmSpline *s                                    = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  self->lnlambda_min = 0.0;
  self->lnlambda_max = 0.0;
  self->reltol       = NC_CLUSTER_RICHNESS_PROJECTION_DEFAULT_RELTOL;
  self->mu           = 0.0;
  self->sigma        = 0.0;
  self->tau          = 0.0;
  self->xi           = 0.0;
  self->T_ode        = ncm_ode_spline_new (s, &_nc_cluster_richness_projection_dTdx);
  self->T            = NULL;
  self->prepared     = FALSE;

  ncm_spline_free (s);
}

static void
_nc_cluster_richness_projection_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterRichnessProjection *crp                = NC_CLUSTER_RICHNESS_PROJECTION (object);
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  g_return_if_fail (NC_IS_CLUSTER_RICHNESS_PROJECTION (object));

  switch (prop_id)
  {
    case PROP_LNLAMBDA_MIN:
      nc_cluster_richness_projection_set_lnlambda_range (crp, g_value_get_double (value), self->lnlambda_max);
      break;
    case PROP_LNLAMBDA_MAX:
      nc_cluster_richness_projection_set_lnlambda_range (crp, self->lnlambda_min, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      nc_cluster_richness_projection_set_reltol (crp, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_cluster_richness_projection_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterRichnessProjection *crp                = NC_CLUSTER_RICHNESS_PROJECTION (object);
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  g_return_if_fail (NC_IS_CLUSTER_RICHNESS_PROJECTION (object));

  switch (prop_id)
  {
    case PROP_LNLAMBDA_MIN:
      g_value_set_double (value, self->lnlambda_min);
      break;
    case PROP_LNLAMBDA_MAX:
      g_value_set_double (value, self->lnlambda_max);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_cluster_richness_projection_dispose (GObject *object)
{
  NcClusterRichnessProjection *crp                = NC_CLUSTER_RICHNESS_PROJECTION (object);
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  ncm_ode_spline_clear (&self->T_ode);
  self->T        = NULL;
  self->prepared = FALSE;

  G_OBJECT_CLASS (nc_cluster_richness_projection_parent_class)->dispose (object);
}

static void
_nc_cluster_richness_projection_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_cluster_richness_projection_parent_class)->finalize (object);
}

static void
nc_cluster_richness_projection_class_init (NcClusterRichnessProjectionClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_cluster_richness_projection_set_property;
  object_class->get_property = &_nc_cluster_richness_projection_get_property;
  object_class->dispose      = &_nc_cluster_richness_projection_dispose;
  object_class->finalize     = &_nc_cluster_richness_projection_finalize;

  /**
   * NcClusterRichnessProjection:lnlambda-min:
   *
   * Lower end $\ln\lambda_\mathrm{min}$ of the range over which the density is
   * evaluated.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNLAMBDA_MIN,
                                   g_param_spec_double ("lnlambda-min",
                                                        NULL,
                                                        "Minimum ln(lambda)",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterRichnessProjection:lnlambda-max:
   *
   * Upper end $\ln\lambda_\mathrm{max}$ of the range over which the density is
   * evaluated.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNLAMBDA_MAX,
                                   g_param_spec_double ("lnlambda-max",
                                                        NULL,
                                                        "Maximum ln(lambda)",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterRichnessProjection:reltol:
   *
   * Relative tolerance of the ODE solver.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        GSL_DBL_EPSILON, 1.0, NC_CLUSTER_RICHNESS_PROJECTION_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_cluster_richness_projection_new:
 *
 * Creates a new #NcClusterRichnessProjection. The range must be set with
 * nc_cluster_richness_projection_set_lnlambda_range() before preparing.
 *
 * Returns: (transfer full): a new #NcClusterRichnessProjection.
 */
NcClusterRichnessProjection *
nc_cluster_richness_projection_new (void)
{
  return g_object_new (NC_TYPE_CLUSTER_RICHNESS_PROJECTION, NULL);
}

/**
 * nc_cluster_richness_projection_ref:
 * @crp: a #NcClusterRichnessProjection
 *
 * Increases the reference count of @crp by one.
 *
 * Returns: (transfer full): @crp.
 */
NcClusterRichnessProjection *
nc_cluster_richness_projection_ref (NcClusterRichnessProjection *crp)
{
  return g_object_ref (crp);
}

/**
 * nc_cluster_richness_projection_free:
 * @crp: a #NcClusterRichnessProjection
 *
 * Decreases the reference count of @crp by one.
 *
 */
void
nc_cluster_richness_projection_free (NcClusterRichnessProjection *crp)
{
  g_object_unref (crp);
}

/**
 * nc_cluster_richness_projection_clear:
 * @crp: a #NcClusterRichnessProjection
 *
 * If *@crp is not %NULL, decreases its reference count by one and sets *@crp to
 * %NULL.
 *
 */
void
nc_cluster_richness_projection_clear (NcClusterRichnessProjection **crp)
{
  g_clear_object (crp);
}

/**
 * nc_cluster_richness_projection_set_lnlambda_range:
 * @crp: a #NcClusterRichnessProjection
 * @lnlambda_min: minimum $\ln\lambda$
 * @lnlambda_max: maximum $\ln\lambda$
 *
 * Sets the range over which the density is evaluated. Calling it invalidates any
 * previous nc_cluster_richness_projection_prepare().
 *
 */
void
nc_cluster_richness_projection_set_lnlambda_range (NcClusterRichnessProjection *crp, gdouble lnlambda_min, gdouble lnlambda_max)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  g_assert_cmpfloat (lnlambda_min, <, lnlambda_max);
  /* exp (lnlambda_max) enters the ODE right-hand side and must not overflow. */
  g_assert_cmpfloat (lnlambda_max, <, 0.5 * GSL_LOG_DBL_MAX);

  self->lnlambda_min = lnlambda_min;
  self->lnlambda_max = lnlambda_max;
  self->prepared     = FALSE;
}

/**
 * nc_cluster_richness_projection_set_reltol:
 * @crp: a #NcClusterRichnessProjection
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance of the ODE solver. Calling it invalidates any
 * previous nc_cluster_richness_projection_prepare().
 *
 */
void
nc_cluster_richness_projection_set_reltol (NcClusterRichnessProjection *crp, gdouble reltol)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  g_assert_cmpfloat (reltol, >=, GSL_DBL_EPSILON);
  g_assert_cmpfloat (reltol, <, 1.0);

  self->reltol   = reltol;
  self->prepared = FALSE;
}

/**
 * nc_cluster_richness_projection_get_reltol:
 * @crp: a #NcClusterRichnessProjection
 *
 * Returns: the relative tolerance of the ODE solver.
 */
gdouble
nc_cluster_richness_projection_get_reltol (NcClusterRichnessProjection *crp)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  return self->reltol;
}

/**
 * nc_cluster_richness_projection_prepare:
 * @crp: a #NcClusterRichnessProjection
 * @mu: log-normal location $\mu$
 * @sigma: log-normal scale $\sigma$
 * @tau: exponential rate $\tau$
 *
 * Solves $\mathrm{d}T/\mathrm{d}x = \tau (g(x) - e^x T)$, with $x = \ln\lambda$ and $g$
 * the log-normal density in $\ln\lambda$, over the range set by
 * nc_cluster_richness_projection_set_lnlambda_range(). One call serves every
 * $\lambda$ in that range.
 *
 */
void
nc_cluster_richness_projection_prepare (NcClusterRichnessProjection *crp, gdouble mu, gdouble sigma, gdouble tau)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);
  gdouble xi, yi, F_LN_i, f_LN_i;
  guint min_sub;

  g_assert_cmpfloat (sigma, >, 0.0);
  g_assert_cmpfloat (tau, >, 0.0);
  g_assert_cmpfloat (self->lnlambda_min, <, self->lnlambda_max);

  self->mu    = mu;
  self->sigma = sigma;
  self->tau   = tau;

  xi = mu - _NC_CLUSTER_RICHNESS_PROJECTION_NSIGMA * sigma;

  g_assert_cmpfloat (xi, <, self->lnlambda_max);

  self->xi = xi;

  min_sub = (guint) ceil ((self->lnlambda_max - xi) * _NC_CLUSTER_RICHNESS_PROJECTION_KNOTS_PER_SIGMA / sigma);

  /* Initial value at xi. T is bounded by tau * F_LN and by f_LN, each attained in
   * one of the two limits of tau, so the smaller is correct to the 1.0e-23 scale
   * of both. A zero initial value leaves the solver with no scale and it fails
   * the error test on the first step. */
  F_LN_i = 0.5 * erfc (_NC_CLUSTER_RICHNESS_PROJECTION_NSIGMA / M_SQRT2);
  f_LN_i = exp (-0.5 * gsl_pow_2 (_NC_CLUSTER_RICHNESS_PROJECTION_NSIGMA))
           / (ncm_c_sqrt_2pi () * sigma * exp (xi));
  yi = MIN (tau * F_LN_i, f_LN_i);

  ncm_ode_spline_set_interval (self->T_ode, yi, xi, self->lnlambda_max);
  ncm_ode_spline_set_reltol (self->T_ode, self->reltol);
  ncm_ode_spline_set_abstol (self->T_ode, self->reltol * yi);
  ncm_ode_spline_set_min_subdivisions (self->T_ode, MAX (min_sub, 16));

  /* The step derived from the initial slope is too small to advance x. */
  ncm_ode_spline_set_ini_step (self->T_ode, (self->lnlambda_max - xi) / MAX (min_sub, 16));
  ncm_ode_spline_prepare (self->T_ode, self);

  self->T        = ncm_ode_spline_peek_spline (self->T_ode);
  self->prepared = TRUE;
}

/**
 * nc_cluster_richness_projection_eval:
 * @crp: a #NcClusterRichnessProjection
 * @lnlambda: $\ln\lambda$
 *
 * Evaluates the density with respect to $\mathrm{d}\lambda$. Below $\mu - 10\sigma$,
 * where $T$ is under $10^{-23}$ of its peak, zero is returned.
 *
 * Returns: $T(\lambda)$.
 */
gdouble
nc_cluster_richness_projection_eval (NcClusterRichnessProjection *crp, gdouble lnlambda)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  g_assert (self->prepared);
  g_assert_cmpfloat (lnlambda, >=, self->lnlambda_min);
  g_assert_cmpfloat (lnlambda, <=, self->lnlambda_max);

  if (lnlambda < self->xi)
    return 0.0;

  return ncm_spline_eval (self->T, lnlambda);
}

/**
 * nc_cluster_richness_projection_eval_lnlambda:
 * @crp: a #NcClusterRichnessProjection
 * @lnlambda: $\ln\lambda$
 *
 * Evaluates the density with respect to $\mathrm{d}\ln\lambda$, that is
 * $\lambda T(\lambda)$.
 *
 * Returns: $\lambda T(\lambda)$.
 */
gdouble
nc_cluster_richness_projection_eval_lnlambda (NcClusterRichnessProjection *crp, gdouble lnlambda)
{
  return exp (lnlambda) * nc_cluster_richness_projection_eval (crp, lnlambda);
}

/**
 * nc_cluster_richness_projection_eval_int:
 * @crp: a #NcClusterRichnessProjection
 * @lnlambda_lo: lower end $\ln\lambda_\mathrm{lo}$
 * @lnlambda_hi: upper end $\ln\lambda_\mathrm{hi}$
 *
 * Integrates the density over $[\lambda_\mathrm{lo}, \lambda_\mathrm{hi}]$ using
 * $\int_0^\lambda T = F_\mathrm{LN}(\lambda) - T(\lambda)/\tau$, a consequence of the
 * ODE, so no quadrature is involved. The two terms cancel to leading order for
 * $\tau\lambda \ll 1$, costing about $\epsilon / (\tau\lambda)$ in relative accuracy.
 *
 * Returns: $\int_{\lambda_\mathrm{lo}}^{\lambda_\mathrm{hi}} T(\lambda) \, \mathrm{d}\lambda$.
 */
gdouble
nc_cluster_richness_projection_eval_int (NcClusterRichnessProjection *crp, gdouble lnlambda_lo, gdouble lnlambda_hi)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);
  const gdouble a_lo                              = (self->mu - lnlambda_lo) / (M_SQRT2 * self->sigma);
  const gdouble a_hi                              = (self->mu - lnlambda_hi) / (M_SQRT2 * self->sigma);
  const gdouble T_lo                              = nc_cluster_richness_projection_eval (crp, lnlambda_lo);
  const gdouble T_hi                              = nc_cluster_richness_projection_eval (crp, lnlambda_hi);
  gdouble dF;

  g_assert_cmpfloat (lnlambda_lo, <=, lnlambda_hi);

  /* F_LN(hi) - F_LN(lo), taking the erfc argument on the side where erfc is small
   * so the difference never subtracts two saturated values. */
  if (a_hi >= 0.0)
    dF = 0.5 * (erfc (a_hi) - erfc (a_lo));
  else
    dF = 0.5 * (erfc (-a_lo) - erfc (-a_hi));

  return dF - (T_hi - T_lo) / self->tau;
}


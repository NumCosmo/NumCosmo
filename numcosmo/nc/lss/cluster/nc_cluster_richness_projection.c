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
 * Writing the integral over the standardized $\ln t$, that is $u = (\ln t -
 * \mu)/\sigma$, puts it in the form
 * \begin{equation}
 * T(\lambda) = \frac{\tau}{\sqrt{2\pi}} \int_{-\infty}^{X} e^{\psi(u)} \mathrm{d}u,
 * \qquad \psi(u) = -\frac{u^2}{2} - s \left(1 - e^{-\sigma (X - u)}\right),
 * \end{equation}
 * with $X = (\ln\lambda - \mu)/\sigma$ and $s = \tau\lambda$. The integrand is
 * positive and bounded by the standard normal, and is evaluated by a composite
 * Gauss-Legendre rule: a single call returns one $\lambda$, so nothing is computed
 * that the caller does not consume. Its integral follows from a closed-form
 * identity rather than a second quadrature. See
 * <a href="../../theory/nc/lss/cluster/cluster_richness_projection.html">Richness projection</a>
 * for the derivation.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/cluster/nc_cluster_richness_projection.h"
#include "ncm/core/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Lower limit of the integration, in units of sigma below mu. The log-normal CDF
 * there, 1.0e-23, bounds both T and the error of the truncation. */
#define _NC_CLUSTER_RICHNESS_PROJECTION_NSIGMA (10.0)

/* Widest panel of the composite rule, in units of sigma. One unit is the scale on
 * which the log-normal bulk varies, so this resolves it with room to spare. */
#define _NC_CLUSTER_RICHNESS_PROJECTION_STEP_CAP (2.0)

/* Geometric grading of the panels toward u = X, where the exponential kernel
 * confines the integrand to a layer of width 1 / (tau * lambda * sigma). */
#define _NC_CLUSTER_RICHNESS_PROJECTION_RATIO (4.0)

/* Smallest panel allowed, so that a degenerate layer cannot stall the march. */
#define _NC_CLUSTER_RICHNESS_PROJECTION_MIN_STEP (1.0e-14)

typedef struct _NcClusterRichnessProjectionPrivate
{
  gdouble reltol;
  gdouble log_drop;
  gdouble mu;
  gdouble sigma;
  gdouble tau;
  gsl_integration_glfixed_table *gl;
  gboolean prepared;
} NcClusterRichnessProjectionPrivate;

struct _NcClusterRichnessProjection
{
  GObject parent_instance;
};

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcClusterRichnessProjection, nc_cluster_richness_projection, G_TYPE_OBJECT)

/*
 * psi(u) = -u^2/2 - s (1 - e^{-sigma (X - u)}), the log of the integrand up to the
 * normal prefactor. Written with expm1 so that the kernel keeps its relative
 * accuracy inside the layer, where sigma (X - u) underflows the cancellation.
 */
static inline gdouble
_nc_cluster_richness_projection_psi (const gdouble u, const gdouble X, const gdouble s, const gdouble sigma)
{
  return -0.5 * u * u + s * expm1 (-sigma * (X - u));
}

/* psi'(u). Used only to tell whether the march has passed the peak. */
static inline gdouble
_nc_cluster_richness_projection_dpsi (const gdouble u, const gdouble X, const gdouble s, const gdouble sigma)
{
  return -u + s * sigma * exp (-sigma * (X - u));
}

static void
nc_cluster_richness_projection_init (NcClusterRichnessProjection *crp)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  self->reltol   = 0.0;
  self->log_drop = 0.0;
  self->mu       = 0.0;
  self->sigma    = 0.0;
  self->tau      = 0.0;
  self->gl       = NULL;
  self->prepared = FALSE;

  nc_cluster_richness_projection_set_reltol (crp, NC_CLUSTER_RICHNESS_PROJECTION_DEFAULT_RELTOL);
}

static void
_nc_cluster_richness_projection_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterRichnessProjection *crp = NC_CLUSTER_RICHNESS_PROJECTION (object);

  g_return_if_fail (NC_IS_CLUSTER_RICHNESS_PROJECTION (object));

  switch (prop_id)
  {
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

  g_clear_pointer (&self->gl, gsl_integration_glfixed_table_free);
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
   * NcClusterRichnessProjection:reltol:
   *
   * Relative accuracy asked of the quadrature. It selects the order of the
   * Gauss-Legendre rule and the depth at which the integrand is treated as
   * negligible; the cost grows with the logarithm of it, not with its inverse.
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
 * Creates a new #NcClusterRichnessProjection.
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
 * nc_cluster_richness_projection_set_reltol:
 * @crp: a #NcClusterRichnessProjection
 * @reltol: relative tolerance
 *
 * Sets the relative accuracy asked of the quadrature.
 *
 */
void
nc_cluster_richness_projection_set_reltol (NcClusterRichnessProjection *crp, gdouble reltol)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);
  gsize nodes;

  g_assert_cmpfloat (reltol, >=, GSL_DBL_EPSILON);
  g_assert_cmpfloat (reltol, <, 1.0);

  if (reltol == self->reltol)
    return;

  /* Measured against a refined rule over the whole parameter box: 1.3e-6 at
   * eight nodes, 1.4e-11 at twelve, and the double-precision floor of the
   * composite rule, 1.1e-11, from sixteen on. */
  if (reltol > 1.0e-5)
    nodes = 8;
  else if (reltol > 1.0e-10)
    nodes = 12;
  else
    nodes = 16;

  /* Depth, in e-folds below the peak, at which the integrand stops contributing. */
  self->log_drop = CLAMP (-log (reltol) + 25.0, 40.0, 80.0);
  self->reltol   = reltol;

  g_clear_pointer (&self->gl, gsl_integration_glfixed_table_free);
  self->gl = gsl_integration_glfixed_table_alloc (nodes);
}

/**
 * nc_cluster_richness_projection_get_reltol:
 * @crp: a #NcClusterRichnessProjection
 *
 * Returns: the relative accuracy asked of the quadrature.
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
 * Sets the parameters the next evaluations refer to. The work is done per
 * $\lambda$ by nc_cluster_richness_projection_eval(), so this call is $O(1)$ and
 * carries no range of its own.
 *
 */
void
nc_cluster_richness_projection_prepare (NcClusterRichnessProjection *crp, gdouble mu, gdouble sigma, gdouble tau)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);

  g_assert_cmpfloat (sigma, >, 0.0);
  g_assert_cmpfloat (tau, >, 0.0);

  self->mu       = mu;
  self->sigma    = sigma;
  self->tau      = tau;
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
 * The integrand $e^{\psi}$ carries two scales: the log-normal bulk, of width one
 * in $u$, and a layer of width $1/(s\sigma)$ at the upper limit, left by the
 * exponential kernel. Since $\psi'' = -1 + s\sigma^2 e^{-\sigma(X-u)}$ grows with
 * $u$, the integrand is concave below $u_\mathrm{inf} = X - \ln(s\sigma^2)/\sigma$
 * and convex above it, so the two scales are the only ones there are. Panels
 * march down from $X$, geometrically graded so that the first one matches the
 * layer, and the march stops once it is past the peak of the concave part and the
 * integrand has dropped below the accuracy asked for.
 *
 * Returns: $T(\lambda)$.
 */
gdouble
nc_cluster_richness_projection_eval (NcClusterRichnessProjection *crp, gdouble lnlambda)
{
  NcClusterRichnessProjectionPrivate * const self = nc_cluster_richness_projection_get_instance_private (crp);
  const gdouble sigma                             = self->sigma;
  const gdouble X                                 = (lnlambda - self->mu) / sigma;
  const gdouble s                                 = self->tau * exp (lnlambda);
  const gdouble nsigma                            = _NC_CLUSTER_RICHNESS_PROJECTION_NSIGMA;
  const gdouble L                                 = X + nsigma;
  const gsize nodes                               = self->gl->n;
  gdouble u_inf, step, v, tot, psi_max;

  g_assert (self->prepared);

  if (X <= -nsigma)
    return 0.0;

  /* tau * lambda beyond the double range: nothing is added, T is the log-normal. */
  if (!gsl_finite (s))
    return exp (-0.5 * X * X - lnlambda) / (ncm_c_sqrt_2pi () * sigma);

  {
    const gdouble ssig2 = s * sigma * sigma;

    u_inf = (ssig2 > 1.0) ? X - log (ssig2) / sigma : X;
  }

  step = MIN (_NC_CLUSTER_RICHNESS_PROJECTION_STEP_CAP,
              MAX (1.0 / (s * sigma), _NC_CLUSTER_RICHNESS_PROJECTION_MIN_STEP));

  v       = 0.0;
  tot     = 0.0;
  psi_max = -G_MAXDOUBLE;

  while (v < L)
  {
    const gdouble b = X - v;
    gdouble a, psi_a, psi_b;

    step  = MIN (step, L - v);
    a     = b - step;
    psi_a = _nc_cluster_richness_projection_psi (a, X, s, sigma);
    psi_b = _nc_cluster_richness_projection_psi (b, X, s, sigma);

    /* Where psi is convex the panel lies below its endpoints, so a panel whose
     * endpoints are already negligible can be skipped outright. */
    if (!((a >= u_inf) && (MAX (psi_a, psi_b) < psi_max - self->log_drop)))
    {
      gsize i;

      for (i = 0; i < nodes; i++)
      {
        gdouble u_i, w_i, psi_i;

        gsl_integration_glfixed_point (a, b, i, &u_i, &w_i, self->gl);
        psi_i    = _nc_cluster_richness_projection_psi (u_i, X, s, sigma);
        tot     += w_i * exp (psi_i);
        psi_max  = MAX (psi_max, psi_i);
      }
    }

    psi_max = MAX (psi_max, MAX (psi_a, psi_b));
    v      += step;

    /* Concave from here down and past the peak: psi only decreases further left. */
    if ((a < u_inf) &&
        (_nc_cluster_richness_projection_dpsi (a, X, s, sigma) > 0.0) &&
        (psi_a < psi_max - self->log_drop))
      break;

    step = MIN (_NC_CLUSTER_RICHNESS_PROJECTION_STEP_CAP,
                _NC_CLUSTER_RICHNESS_PROJECTION_RATIO * step);
  }

  return self->tau * tot / ncm_c_sqrt_2pi ();
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
 * $\int_0^\lambda T = F_\mathrm{LN}(\lambda) - T(\lambda)/\tau$, a property of $T$
 * itself, so no quadrature of the cumulative is involved. The two terms cancel to
 * leading order for $\tau\lambda \ll 1$, costing about $\epsilon / (\tau\lambda)$
 * in relative accuracy.
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

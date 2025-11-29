/***************************************************************************
 *            ncm_pln1d.c
 *
 *  Fri Nov 28 11:24:36 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_pln1d.h
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmPLN1D:
 *
 * A Poisson–Lognormal 1D integrator using mode finding (GSL Lambert-W), shifted
 * Gauss–Hermite, or Laplace fallback.
 *
 * Poisson–Lognormal 1D Distribution
 *
 * Definitions
 * ===========
 *
 * The Poisson component for an integer count R ≥ 0 and rate λ > 0 is
 * $$
 *   P_{\rm P}(R \mid \lambda) = e^{-\lambda} \frac{\lambda^{R}}{R!}.
 * $$
 * The lognormal prior for λ is defined by
 * $$
 *   P_{\rm LN}(\lambda \mid \mu, \sigma) = \frac{1}{\lambda \, \sigma \sqrt{2\pi}}
 *          \exp\!\left[-\frac{(\ln\lambda - \mu)^2}{2\sigma^2}\right],
 * $$
 * where $\mu$ is the mean of $\log\lambda$ and $\sigma$ is the standard deviation of
 * $\log\lambda$.
 *
 * The Poisson–Lognormal likelihood is the integral
 * $$
 *   P(R \mid \mu,\sigma) = \int_{0}^{\infty} P_{\rm P}(R\mid\lambda)\,
 *             P_{\rm LN}(\lambda\mid\mu,\sigma)\; \mathrm{d}\lambda.
 * $$
 *
 * Change of Variables
 * ===================
 *
 * We rewrite $\lambda$ using the substitution
 * $$
 *   \lambda = e^{\sigma x},
 * $$
 * where  $x \in (-\infty,\infty)$ becomes the integration variable.
 *
 * The Jacobian is
 * $$
 *   \mathrm{d}\lambda = \sigma\, e^{\sigma x}\, \mathrm{d}x.
 * $$
 * Applying this change of variables and simplifying yields the integrand
 * $$
 *   I(x; R, \mu, \sigma) = -\frac{R^2\sigma^2}{2} + u R \sigma
 *                        - \left(x - u\right)^2 - e^{\sigma x},
 * $$
 * where $u = \mu/\sigma + R\sigma$. Thus,
 * $$
 *   P(R\mid\mu,\sigma)= \int_{-\infty}^{\infty}
 *   \frac{\exp(I(x;R,\mu,\sigma))}{\sqrt{2\pi}\, R!}\, \mathrm{d}x.
 * $$
 *
 * Integrand Mode Finding
 * ======================
 *
 * To find the mode of the integrand, we set the derivative of I with respect to x to zero:
 * $$
 *   \frac{\mathrm{d}}{\mathrm{d}x}I(x;R,\mu,\sigma) = 0.
 * $$
 * Solving for x gives
 * $$
 *   x = u - \frac{W_0(\sigma^2 \exp(\sigma u))}{\sigma},
 * $$
 * Here, $W_0$ is the principal branch of the Lambert $W$ function.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_pln1d.h"
#include "math/ncm_c.h"
#include "lintegrate/logadd.c"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_lambert.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmPLN1DPrivate
{
  guint gh_order;
  gsl_integration_fixed_workspace *gh_workspace;
  gdouble *nodes;
  gdouble *weights;
  gdouble *ln_weights;
} NcmPLN1DPrivate;


enum
{
  PROP_0,
  PROP_GH_ORDER,
  PROP_LEN,
};

struct _NcmPLN1D
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmPLN1D, ncm_pln1d, G_TYPE_OBJECT);

static void
ncm_pln1d_init (NcmPLN1D *pln1d)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln1d);

  self->gh_order     = 0;
  self->gh_workspace = NULL;
  self->nodes        = NULL;
  self->weights      = NULL;
  self->ln_weights   = NULL;
}

static void
_ncm_pln1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPLN1D *pln1d = NCM_PLN1D (object);

  g_return_if_fail (NCM_IS_PLN1D (object));

  switch (prop_id)
  {
    case PROP_GH_ORDER:
      ncm_pln1d_set_order (pln1d, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_pln1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPLN1D *pln1d = NCM_PLN1D (object);

  g_return_if_fail (NCM_IS_PLN1D (object));

  switch (prop_id)
  {
    case PROP_GH_ORDER:
      g_value_set_uint (value, ncm_pln1d_get_order (pln1d));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_pln1d_dispose (GObject *object)
{
  NcmPLN1D *pln1d              = NCM_PLN1D (object);
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln1d);

  if (self->gh_workspace != NULL)
  {
    gsl_integration_fixed_free (self->gh_workspace);
    self->gh_workspace = NULL;
    self->gh_order     = 0;
    self->nodes        = NULL;
    self->weights      = NULL;
  }


  /* Chain up : end */
  G_OBJECT_CLASS (ncm_pln1d_parent_class)->dispose (object);
}

static void
_ncm_pln1d_finalize (GObject *object)
{
  NcmPLN1D *pln1d              = NCM_PLN1D (object);
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln1d);

  if (self->ln_weights != NULL)
  {
    g_free (self->ln_weights);
    self->ln_weights = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_pln1d_parent_class)->finalize (object);
}

static void
ncm_pln1d_class_init (NcmPLN1DClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_pln1d_set_property;
  object_class->get_property = &_ncm_pln1d_get_property;
  object_class->dispose      = &_ncm_pln1d_dispose;
  object_class->finalize     = &_ncm_pln1d_finalize;

  /**
   * NcmPLN1D:gh-order:
   *
   * The Gauss-Hermite order to be used in the integration.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_GH_ORDER,
                                   g_param_spec_uint ("gh-order",
                                                      "Gauss-Hermite order",
                                                      "Order of the Gauss-Hermite quadrature to be used in the integration",
                                                      0, G_MAXUINT, 60,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_pln1d_new:
 *
 * Creates a new #NcmPLN1D object.
 *
 * Returns: a new #NcmPLN1D.
 */
NcmPLN1D *
ncm_pln1d_new (guint gh_order)
{
  NcmPLN1D *pln1d = g_object_new (NCM_TYPE_PLN1D,
                                  "gh-order", gh_order,
                                  NULL);

  return pln1d;
}

/**
 * ncm_pln1d_ref:
 * @pln1d: a #NcmPLN1D
 *
 * Increase the reference of @pln1d by one.
 *
 * Returns: (transfer full): @pln1d.
 */
NcmPLN1D *
ncm_pln1d_ref (NcmPLN1D *pln1d)
{
  return g_object_ref (pln1d);
}

/**
 * ncm_pln1d_free:
 * @pln1d: a #NcmPLN1D
 *
 * Decrease the reference count of @pln1d by one.
 *
 */
void
ncm_pln1d_free (NcmPLN1D *pln1d)
{
  g_object_unref (pln1d);
}

/**
 * ncm_pln1d_clear:
 * @pln1d: a #NcmPLN1D
 *
 * Decrease the reference count of @pln1d by one, and sets the pointer *@pln1d to NULL.
 *
 */
void
ncm_pln1d_clear (NcmPLN1D **pln1d)
{
  g_clear_object (pln1d);
}

/* Configuration */
void
ncm_pln1d_set_order (NcmPLN1D *pln, guint gh_order)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln);

  g_assert_cmpuint (gh_order, >, 0);

  if (self->gh_order == gh_order)
    return;

  if (self->gh_workspace != NULL)
  {
    gsl_integration_fixed_free (self->gh_workspace);
    self->gh_workspace = NULL;
  }

  self->gh_workspace = gsl_integration_fixed_alloc (gsl_integration_fixed_hermite,
                                                    gh_order, 0.0, 0.5, 0.0, 0.0);
  self->nodes   = gsl_integration_fixed_nodes (self->gh_workspace);
  self->weights = gsl_integration_fixed_weights (self->gh_workspace);

  if (self->ln_weights != NULL)
    g_free (self->ln_weights);

  self->ln_weights = g_new (gdouble, gh_order);

  for (guint i = 0; i < gh_order; i++)
    self->ln_weights[i] = log (self->weights[i]);

  self->gh_order = gh_order;
}

guint
ncm_pln1d_get_order (NcmPLN1D *pln)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln);

  return self->gh_order;
}

/**
 * ncm_pln1d_mode:
 * @R: Poisson rate parameter
 * @mu: log-normal location parameter (log-space mean)
 * @sigma: log-normal scale parameter (log-space standard deviation)
 *
 * Compute the mode of the Poisson–Lognormal integrand.
 *
 * Returns: the mode of the integrand.
 */
double
ncm_pln1d_mode (gdouble R, gdouble mu, gdouble sigma)
{
  const gdouble s2 = sigma * sigma;
  const gdouble u  = mu / sigma + R * sigma;
  const gdouble y  = s2 * exp (sigma * u);
  const gdouble W  = gsl_sf_lambert_W0 (y);

  return u - W / sigma;
}

static inline gdouble
_ncm_pln1d_log_integrand0 (gdouble x, gdouble u, gdouble sigma)
{
  return u * x - exp (sigma * x);
}

static gdouble
_ncm_pln1d_eval_gh_lnp (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln);
  const guint n                = self->gh_order;
  const gdouble z              = ncm_pln1d_mode (R, mu, sigma);
  const gdouble u              = mu / sigma + R * sigma;
  const gdouble ln_scale       = -0.5 * (gsl_pow_2 (mu / sigma) + z * z) - lgamma (R + 1.0);
  const gdouble ln_2pi         = ncm_c_ln2pi ();
  gdouble ln_sum               = -INFINITY;
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble y    = self->nodes[i];
    const gdouble ln_w = self->ln_weights[i];
    const gdouble logf = _ncm_pln1d_log_integrand0 (y + z, u, sigma) - y * z;

    ln_sum = logaddexp (ln_sum, ln_w + logf);
  }

  return ln_sum + ln_scale - 0.5 * ln_2pi;
}

static gdouble
_ncm_pln1d_eval_laplace_lnp (gdouble R, gdouble mu, gdouble sigma)
{
  const gdouble z       = ncm_pln1d_mode (R, mu, sigma);
  const gdouble u       = mu / sigma + R * sigma;
  const gdouble umz     = u - z;
  const gdouble R_sigma = R * sigma;
  const gdouble lnf     = -0.5 * (umz * umz + R_sigma * R_sigma) - umz / sigma + R_sigma * u;

  return lnf - lgamma (R + 1.0) - 0.5 * log1p (sigma * umz);
}

/**
 * ncm_pln1d_eval_range_sum_lnp:
 * @pln: a #NcmPLN1D
 * @R_min: minimum Poisson rate parameter (inclusive)
 * @R_max: maximum Poisson rate parameter (inclusive)
 * @mu: log-normal mean (location parameter)
 * @sigma: log-normal standard deviation (scale parameter)
 *
 * Evaluate the cumulative Poisson–Lognormal probability by summing over a range of R values.
 * This computes: ∑_{R=R_min}^{R_max} P(R|μ,σ)
 *
 * This is more efficient than calling ncm_pln1d_eval_p repeatedly
 * as it reuses the shifted Gauss-Hermite nodes for all R values.
 *
 * The shift point z is computed using the central R value to ensure
 * good accuracy across the range.
 *
 * Returns: the logarithm of the cumulative probability.
 */
gdouble
ncm_pln1d_eval_range_sum_lnp (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln);
  const guint n                = self->gh_order;
  const gdouble ln_2pi         = ncm_c_ln2pi ();
  gdouble ln_sum               = -INFINITY;

  g_assert_cmpuint (R_min, <=, R_max);

  if (sigma > 1.0e-4)
  {
    /* Use Gauss-Hermite quadrature with a single shift point for the range */
    /* Compute z using the central R value for better accuracy across range */
    const gdouble R_central = 0.5 * (R_min + R_max);
    const gdouble z         = ncm_pln1d_mode (R_central, mu, sigma);
    const gdouble u_base    = mu / sigma;

    for (guint j = R_min; j <= R_max; j++)
    {
      const gdouble R        = (gdouble) j;
      const gdouble u        = u_base + R * sigma;
      const gdouble ln_scale = -0.5 * (gsl_pow_2 (u_base) + z * z) - lgamma (R + 1.0);
      gdouble ln_R_sum       = -INFINITY;

      for (guint i = 0; i < n; i++)
      {
        const gdouble y    = self->nodes[i];
        const gdouble ln_w = self->ln_weights[i];
        const gdouble logf = _ncm_pln1d_log_integrand0 (y + z, u, sigma) - y * z;

        ln_R_sum = logaddexp (ln_R_sum, ln_w + logf);
      }

      {
        const gdouble p_R = ln_R_sum + ln_scale - 0.5 * ln_2pi;

        ln_sum = logaddexp (ln_sum, p_R);
      }
    }
  }
  else
  {
    /* Use Laplace approximation for small sigma */
    for (guint j = R_min; j <= R_max; j++)
    {
      const gdouble R   = (gdouble) j;
      const gdouble p_R = _ncm_pln1d_eval_laplace_lnp (R, mu, sigma);

      ln_sum = logaddexp (ln_sum, p_R);
    }
  }

  return ln_sum;
}

/**
 * ncm_pln1d_eval_range_sum:
 * @pln: a #NcmPLN1D
 * @R_min: minimum Poisson rate parameter (inclusive)
 * @R_max: maximum Poisson rate parameter (inclusive)
 * @mu: log-normal mean (location parameter)
 * @sigma: log-normal standard deviation (scale parameter)
 *
 * Evaluate the cumulative Poisson–Lognormal probability by summing over a range of R values.
 * This computes: ∑_{R=R_min}^{R_max} P(R|μ,σ)
 *
 * This is the non-log version of ncm_pln1d_eval_range_sum_lnp.
 *
 * Returns: the cumulative probability.
 */
gdouble
ncm_pln1d_eval_range_sum (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma)
{
  return exp (ncm_pln1d_eval_range_sum_lnp (pln, R_min, R_max, mu, sigma));
}

/**
 * ncm_pln1d_eval_lnp:
 * @pln: a #NcmPLN1D
 * @R: Poisson rate parameter
 * @mu: log-normal mean (location parameter)
 * @sigma: log-normal standard deviation (scale parameter)
 *
 * Evaluate the Poisson–Lognormal integral.
 *
 * Returns: the logarithm of the integral.
 */
gdouble
ncm_pln1d_eval_lnp (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma)
{
  if (sigma > 1.0e-4)
    return _ncm_pln1d_eval_gh_lnp (pln, R, mu, sigma);

  return _ncm_pln1d_eval_laplace_lnp (R, mu, sigma);
}

/**
 * ncm_pln1d_eval_p:
 * @pln: a #NcmPLN1D
 * @R: Poisson rate parameter
 * @mu: log-normal mean (location parameter)
 * @sigma: log-normal standard deviation (scale parameter)
 *
 * Evaluate the Poisson–Lognormal integral.
 *
 * Returns: the value of the integral.
 */
gdouble
ncm_pln1d_eval_p (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma)
{
  return exp (ncm_pln1d_eval_lnp (pln, R, mu, sigma));
}


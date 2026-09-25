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
 * Poisson-Lognormal probability of an integer count.
 *
 * Evaluates
 * $$
 *   P(R \mid \mu, \sigma) = \int_0^\infty e^{-\lambda} \frac{\lambda^R}{R!}
 *   \frac{1}{\lambda \sigma \sqrt{2\pi}}
 *   \exp\left[-\frac{(\ln\lambda - \mu)^2}{2\sigma^2}\right] \mathrm{d}\lambda,
 * $$
 * where $\mu$ and $\sigma$ are the mean and standard deviation of $\ln\lambda$. For
 * $\sigma > 10^{-4}$ it uses Gauss-Hermite quadrature with nodes centered on the mode of
 * the integrand and scaled to the width of its peak, otherwise the Laplace approximation
 * about the mode.
 *
 * See <a href="../../theory/poisson_lognormal.html">Poisson-Lognormal Distribution</a>
 * for the derivation.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_pln1d.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_util.h"
#include "external/lintegrate/logadd.c"

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
   * Number of Gauss-Hermite nodes.
   */
  g_object_class_install_property (object_class,
                                   PROP_GH_ORDER,
                                   g_param_spec_uint ("gh-order",
                                                      "Gauss-Hermite order",
                                                      "Order of the Gauss-Hermite quadrature to be used in the integration",
                                                      1, G_MAXUINT, 60,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_pln1d_new:
 * @gh_order: number of Gauss-Hermite nodes
 *
 * Creates a new #NcmPLN1D.
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
 * Increases the reference count of @pln1d by one.
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
 * Decreases the reference count of @pln1d by one.
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
 * Decreases the reference count of *@pln1d by one and sets *@pln1d to %NULL.
 */
void
ncm_pln1d_clear (NcmPLN1D **pln1d)
{
  g_clear_object (pln1d);
}

/**
 * ncm_pln1d_set_order:
 * @pln: a #NcmPLN1D
 * @gh_order: number of Gauss-Hermite nodes
 *
 * Sets the number of Gauss-Hermite nodes, which must be positive.
 */
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

  /* Logarithm of the weights times e^{y^2/2}, the Gauss-Hermite weight function */
  for (guint i = 0; i < gh_order; i++)
    self->ln_weights[i] = log (self->weights[i]) + 0.5 * gsl_pow_2 (self->nodes[i]);

  self->gh_order = gh_order;
}

/**
 * ncm_pln1d_get_order:
 * @pln: a #NcmPLN1D
 *
 * Returns: the number of Gauss-Hermite nodes.
 */
guint
ncm_pln1d_get_order (NcmPLN1D *pln)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln);

  return self->gh_order;
}

/**
 * ncm_pln1d_mode:
 * @R: count
 * @mu: mean of $\ln\lambda$
 * @sigma: standard deviation of $\ln\lambda$
 *
 * Computes the mode of the Poisson-Lognormal integrand in the variable
 * $x = \ln(\lambda)/\sigma$.
 *
 * Returns: the mode $z$.
 */
double
ncm_pln1d_mode (gdouble R, gdouble mu, gdouble sigma)
{
  const gdouble u    = mu / sigma + R * sigma;
  const gdouble ln_y = 2.0 * log (sigma) + sigma * u;
  const gdouble W    = ncm_util_lambert_W0_ln (ln_y);

  return u - W / sigma;
}

static gdouble
_ncm_pln1d_eval_gh_lnp (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma)
{
  NcmPLN1DPrivate * const self = ncm_pln1d_get_instance_private (pln);
  const guint n                = self->gh_order;
  const gdouble z              = ncm_pln1d_mode (R, mu, sigma);
  const gdouble u              = mu / sigma + R * sigma;
  const gdouble zmu            = z - u;
  const gdouble scale          = 1.0 / sqrt (1.0 - sigma * zmu);
  const gdouble ln_const       = log (scale) + R * mu + 0.5 * gsl_pow_2 (R * sigma) - lgamma (R + 1.0) - 0.5 * ncm_c_ln2pi ();
  gdouble ln_sum               = -INFINITY;
  guint i;

  /* Nodes at x = z + scale y, with scale the width of the peak of I(x) */
  for (i = 0; i < n; i++)
  {
    const gdouble y    = self->nodes[i];
    const gdouble t    = zmu + scale * y;
    const gdouble logf = self->ln_weights[i] - 0.5 * t * t - exp (sigma * (z + scale * y));

    ln_sum = logaddexp (ln_sum, logf);
  }

  return ln_sum + ln_const;
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
 * @R_min: smallest count
 * @R_max: largest count
 * @mu: mean of $\ln\lambda$
 * @sigma: standard deviation of $\ln\lambda$
 *
 * Computes $\ln\sum_{R=R_\mathrm{min}}^{R_\mathrm{max}} P(R \mid \mu, \sigma)$, with
 * $R_\mathrm{min} \leq R_\mathrm{max}$.
 *
 * Returns: the logarithm of the summed probability.
 */
gdouble
ncm_pln1d_eval_range_sum_lnp (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma)
{
  gdouble ln_sum = -INFINITY;
  guint R;

  g_assert_cmpuint (R_min, <=, R_max);

  for (R = R_min; R <= R_max; R++)
    ln_sum = logaddexp (ln_sum, ncm_pln1d_eval_lnp (pln, R, mu, sigma));

  return ln_sum;
}

/**
 * ncm_pln1d_eval_range_sum:
 * @pln: a #NcmPLN1D
 * @R_min: smallest count
 * @R_max: largest count
 * @mu: mean of $\ln\lambda$
 * @sigma: standard deviation of $\ln\lambda$
 *
 * Exponential of ncm_pln1d_eval_range_sum_lnp().
 *
 * Returns: the summed probability.
 */
gdouble
ncm_pln1d_eval_range_sum (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma)
{
  return exp (ncm_pln1d_eval_range_sum_lnp (pln, R_min, R_max, mu, sigma));
}

/**
 * ncm_pln1d_eval_lnp:
 * @pln: a #NcmPLN1D
 * @R: count
 * @mu: mean of $\ln\lambda$
 * @sigma: standard deviation of $\ln\lambda$
 *
 * Computes $\ln P(R \mid \mu, \sigma)$.
 *
 * Returns: the logarithm of the probability.
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
 * @R: count
 * @mu: mean of $\ln\lambda$
 * @sigma: standard deviation of $\ln\lambda$
 *
 * Exponential of ncm_pln1d_eval_lnp().
 *
 * Returns: the probability.
 */
gdouble
ncm_pln1d_eval_p (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma)
{
  return exp (ncm_pln1d_eval_lnp (pln, R, mu, sigma));
}


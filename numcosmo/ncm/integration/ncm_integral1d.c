/***************************************************************************
 *            ncm_integral1d.c
 *
 *  Sat February 20 14:29:30 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_integral1d.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmIntegral1d:
 *
 * Abstract class for one-dimensional integrals of an integrand $F$.
 *
 * Subclasses provide the integrand, see #NcmIntegral1dPtr. Every integral uses the
 * adaptive Gauss-Kronrod quadrature of GSL, gsl_integration_qag(), with
 * #NcmIntegral1d:rule, at most #NcmIntegral1d:partition subintervals and the tolerances
 * #NcmIntegral1d:reltol and #NcmIntegral1d:abstol; a GSL failure aborts. The integrals
 * over infinite ranges with a Gaussian or exponential weight change to a variable
 * $\alpha$ on a finite interval, stated in each function, and pass $\alpha$ to the
 * integrand as its argument $w$; the finite-range integrals pass $w = 1$. The GSL
 * workspace belongs to the object, so evaluation is not reentrant.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integral1d.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "external/lintegrate/lintegrate.h"
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmIntegral1dPrivate
{
  guint partition;
  gdouble reltol;
  gdouble abstol;
  guint rule;
  gsl_integration_workspace *ws;
  gsl_integration_cquad_workspace *cquad_ws;
} NcmIntegral1dPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmIntegral1d, ncm_integral1d, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_INTEGRAND,
  PROP_PARTITION,
  PROP_RULE,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_SIZE,
};

static void
ncm_integral1d_init (NcmIntegral1d *int1d)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  self->partition = 0;
  self->rule      = 0;
  self->reltol    = 0.0;
  self->abstol    = 0.0;
  self->ws        = NULL;
  self->cquad_ws  = NULL;
}

static void
ncm_integral1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegral1d *int1d = NCM_INTEGRAL1D (object);

  g_return_if_fail (NCM_IS_INTEGRAL1D (object));

  switch (prop_id)
  {
    case PROP_PARTITION:
      ncm_integral1d_set_partition (int1d, g_value_get_uint (value));
      break;
    case PROP_RULE:
      ncm_integral1d_set_rule (int1d, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      ncm_integral1d_set_reltol (int1d, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_integral1d_set_abstol (int1d, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_integral1d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmIntegral1d *int1d = NCM_INTEGRAL1D (object);

  g_return_if_fail (NCM_IS_INTEGRAL1D (object));

  switch (prop_id)
  {
    case PROP_PARTITION:
      g_value_set_uint (value, ncm_integral1d_get_partition (int1d));
      break;
    case PROP_RULE:
      g_value_set_uint (value, ncm_integral1d_get_rule (int1d));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_integral1d_get_reltol (int1d));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_integral1d_get_abstol (int1d));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_integral1d_finalize (GObject *object)
{
  NcmIntegral1d *int1d              = NCM_INTEGRAL1D (object);
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  g_clear_pointer (&self->ws,       gsl_integration_workspace_free);
  g_clear_pointer (&self->cquad_ws, gsl_integration_cquad_workspace_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integral1d_parent_class)->finalize (object);
}

static void
ncm_integral1d_class_init (NcmIntegral1dClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_integral1d_set_property;
  object_class->get_property = &ncm_integral1d_get_property;
  object_class->finalize     = &ncm_integral1d_finalize;

  /**
   * NcmIntegral1d:partition:
   *
   * The maximum number of subintervals.
   */
  g_object_class_install_property (object_class,
                                   PROP_PARTITION,
                                   g_param_spec_uint ("partition",
                                                      NULL,
                                                      "Integral maximum partition",
                                                      10, G_MAXUINT32, NCM_INTEGRAL1D_DEFAULT_PARTITION,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegral1d:rule:
   *
   * The Gauss-Kronrod rule, 1 to 6 for 15 to 61 points (GSL_INTEG_GAUSS15 to
   * GSL_INTEG_GAUSS61).
   */
  g_object_class_install_property (object_class,
                                   PROP_RULE,
                                   g_param_spec_uint ("rule",
                                                      NULL,
                                                      "Integration rule",
                                                      1, 6, NCM_INTEGRAL1D_DEFAULT_ALG,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegral1d:reltol:
   *
   * The relative tolerance.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Integral relative tolerance",
                                                        0.0, 1.0, NCM_INTEGRAL1D_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegral1d:abstol:
   *
   * The absolute tolerance.
   */
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Integral absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_INTEGRAL1D_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_integral1d_ref:
 * @int1d: a #NcmIntegral1d
 *
 * Increases the reference count of @int1d by one.
 *
 * Returns: (transfer full): @int1d.
 */
NcmIntegral1d *
ncm_integral1d_ref (NcmIntegral1d *int1d)
{
  return g_object_ref (int1d);
}

/**
 * ncm_integral1d_free:
 * @int1d: a #NcmIntegral1d
 *
 * Decreases the reference count of @int1d by one.
 */
void
ncm_integral1d_free (NcmIntegral1d *int1d)
{
  g_object_unref (int1d);
}

/**
 * ncm_integral1d_clear:
 * @int1d: a #NcmIntegral1d
 *
 * If *@int1d is not %NULL, decreases its reference count by one and sets *@int1d to %NULL.
 */
void
ncm_integral1d_clear (NcmIntegral1d **int1d)
{
  g_clear_object (int1d);
}

/**
 * ncm_integral1d_set_partition:
 * @int1d: a #NcmIntegral1d
 * @partition: the maximum number of subintervals
 *
 * Sets #NcmIntegral1d:partition, at least 10, reallocating the workspace.
 */
void
ncm_integral1d_set_partition (NcmIntegral1d *int1d, guint partition)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  g_assert_cmpuint (partition, >=, 10);

  if (self->partition != partition)
  {
    g_clear_pointer (&self->ws, gsl_integration_workspace_free);

    self->ws        = gsl_integration_workspace_alloc (partition);
    self->partition = partition;
  }
}

/**
 * ncm_integral1d_set_rule:
 * @int1d: a #NcmIntegral1d
 * @rule: the Gauss-Kronrod rule
 *
 * Sets #NcmIntegral1d:rule.
 */
void
ncm_integral1d_set_rule (NcmIntegral1d *int1d, guint rule)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  self->rule = rule;
}

/**
 * ncm_integral1d_set_reltol:
 * @int1d: a #NcmIntegral1d
 * @reltol: the relative tolerance
 *
 * Sets #NcmIntegral1d:reltol.
 */
void
ncm_integral1d_set_reltol (NcmIntegral1d *int1d, gdouble reltol)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  self->reltol = reltol;
}

/**
 * ncm_integral1d_set_abstol:
 * @int1d: a #NcmIntegral1d
 * @abstol: the absolute tolerance
 *
 * Sets #NcmIntegral1d:abstol.
 */
void
ncm_integral1d_set_abstol (NcmIntegral1d *int1d, gdouble abstol)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  self->abstol = abstol;
}

/**
 * ncm_integral1d_get_partition:
 * @int1d: a #NcmIntegral1d
 *
 * Gets #NcmIntegral1d:partition.
 *
 * Returns: the maximum number of subintervals.
 */
guint
ncm_integral1d_get_partition (NcmIntegral1d *int1d)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  return self->partition;
}

/**
 * ncm_integral1d_get_rule:
 * @int1d: a #NcmIntegral1d
 *
 * Gets #NcmIntegral1d:rule.
 *
 * Returns: the Gauss-Kronrod rule.
 */
guint
ncm_integral1d_get_rule (NcmIntegral1d *int1d)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  return self->rule;
}

/**
 * ncm_integral1d_get_reltol:
 * @int1d: a #NcmIntegral1d
 *
 * Gets #NcmIntegral1d:reltol.
 *
 * Returns: the relative tolerance.
 */
gdouble
ncm_integral1d_get_reltol (NcmIntegral1d *int1d)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  return self->reltol;
}

/**
 * ncm_integral1d_get_abstol:
 * @int1d: a #NcmIntegral1d
 *
 * Gets #NcmIntegral1d:abstol.
 *
 * Returns: the absolute tolerance.
 */
gdouble
ncm_integral1d_get_abstol (NcmIntegral1d *int1d)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);

  return self->abstol;
}

/**
 * ncm_integral1d_integrand:
 * @int1d: a #NcmIntegral1d
 * @x: the point
 * @w: the variable $\alpha$ of the change of variables, or 1
 *
 * Returns: the integrand $F(x)$, see #NcmIntegral1d.
 */
gdouble
ncm_integral1d_integrand (NcmIntegral1d *int1d, const gdouble x, const gdouble w)
{
  return NCM_INTEGRAL1D_GET_CLASS (int1d)->integrand (int1d, x, w);
}

typedef struct _NcIntegral1dHermite
{
  NcmIntegral1d *int1d;
  gdouble mu;
  gdouble r;
} NcIntegral1dHermite;

static gdouble
_ncm_integral1d_eval_p (const gdouble x, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;

  return ncm_integral1d_integrand (int1d_H->int1d, x, 1.0);
}

static gdouble
_ncm_integral1d_eval_gauss_hermite_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = gsl_cdf_ugaussian_Qinv (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble
_ncm_integral1d_eval_gauss_hermite1_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = sqrt (-2.0 * log (alpha));

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble
_ncm_integral1d_eval_gauss_hermite_r_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = gsl_cdf_ugaussian_Qinv (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x / int1d_H->r, alpha);
}

static gdouble
_ncm_integral1d_eval_gauss_hermite1_r_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = sqrt (-2.0 * log (alpha)) / int1d_H->r;

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble
_ncm_integral1d_eval_gauss_hermite (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = gsl_cdf_ugaussian_Qinv (alpha);

  return (ncm_integral1d_integrand (int1d_H->int1d, x, alpha) + ncm_integral1d_integrand (int1d_H->int1d, -x, alpha));
}

static gdouble
_ncm_integral1d_eval_gauss_hermite_mur (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = gsl_cdf_ugaussian_Qinv (alpha);
  const gdouble y              = x / int1d_H->r;

  return (ncm_integral1d_integrand (int1d_H->int1d,  int1d_H->mu + y, alpha) + ncm_integral1d_integrand (int1d_H->int1d, int1d_H->mu - y, alpha));
}

static gdouble
_ncm_integral1d_eval_gauss_laguerre (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = -log (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble
_ncm_integral1d_eval_gauss_laguerre_r (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x              = -log (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x / int1d_H->r, alpha);
}

/**
 * ncm_integral1d_eval:
 * @int1d: a #NcmIntegral1d
 * @xi: the lower limit $x_i$
 * @xf: the upper limit $x_f$
 * @err: (out): the error estimate
 *
 * Returns: $\int_{x_i}^{x_f}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval (NcmIntegral1d *int1d, const gdouble xi, const gdouble xf, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, 0.0};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_p;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, xi, xf, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval: %s.", gsl_strerror (ret));

  return result;
}

/**
 * ncm_integral1d_eval_lnint:
 * @int1d: a #NcmIntegral1d
 * @xi: the lower limit $x_i$
 * @xf: the upper limit $x_f$
 * @err: (out): the error estimate
 *
 * Integrates with the integrand taken as the logarithm of the function to integrate,
 * summing in log space so that $e^{F}$ may be outside the double range.
 *
 * Returns: $\ln\int_{x_i}^{x_f}e^{F(x)}\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_lnint (NcmIntegral1d *int1d, const gdouble xi, const gdouble xf, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, 0.0};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_p;
  F.params   = &int1d_H;

  ret = lintegration_qag (&F, xi, xf, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_lnint: %s.", gsl_strerror (ret));

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_p:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = Q(x)$, the upper tail probability of the standard normal
 * distribution, over $(0, 1/2]$.
 *
 * Returns: $\int_0^\infty e^{-x^2/2}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_hermite_p (NcmIntegral1d *int1d, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, 0.0};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_hermite_p;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 0.5, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite_p: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result;
  err[0] = ncm_c_sqrt_2pi () * err[0];

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = Q(|x|)$, the upper tail probability of the standard normal
 * distribution, over $(0, 1/2]$, evaluating $F$ at $\pm x$.
 *
 * Returns: $\int_{-\infty}^\infty e^{-x^2/2}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_hermite (NcmIntegral1d *int1d, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, 0.0};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_hermite;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 0.5, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result;
  err[0] = ncm_c_sqrt_2pi () * err[0];

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_r_p:
 * @int1d: a #NcmIntegral1d
 * @r: the inverse Gaussian width $r > 0$
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = Q(rx)$ over $(0, 1/2]$, see
 * ncm_integral1d_eval_gauss_hermite_p().
 *
 * Returns: $\int_0^\infty e^{-r^2x^2/2}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_hermite_r_p (NcmIntegral1d *int1d, const gdouble r, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, r};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);

  F.function = &_ncm_integral1d_eval_gauss_hermite_r_p;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 0.5, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite_r_p: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result / r;
  err[0] = ncm_c_sqrt_2pi () * err[0] / r;

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_mur:
 * @int1d: a #NcmIntegral1d
 * @r: the inverse Gaussian width $r > 0$
 * @mu: the Gaussian mean $\mu$
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = Q(r|x - \mu|)$ over $(0, 1/2]$, evaluating $F$ at
 * $\mu \pm |x - \mu|$.
 *
 * Returns: $\int_{-\infty}^\infty e^{-r^2(x - \mu)^2/2}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_hermite_mur (NcmIntegral1d *int1d, const gdouble r, const gdouble mu, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, mu, r};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);

  F.function = &_ncm_integral1d_eval_gauss_hermite_mur;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 0.5, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite_mur: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result / r;
  err[0] = ncm_c_sqrt_2pi () * err[0] / r;

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite1_p:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = e^{-x^2/2}$ over $(0, 1]$.
 *
 * Returns: $\int_0^\infty x e^{-x^2/2}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_hermite1_p (NcmIntegral1d *int1d, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, 0.0};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_hermite1_p;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 1.0, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite1_p: %s.", gsl_strerror (ret));

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite1_r_p:
 * @int1d: a #NcmIntegral1d
 * @r: the inverse Gaussian width $r > 0$
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = e^{-r^2x^2/2}$ over $(0, 1]$.
 *
 * Returns: $\int_0^\infty x e^{-r^2x^2/2}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_hermite1_r_p (NcmIntegral1d *int1d, const gdouble r, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, r};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);

  F.function = &_ncm_integral1d_eval_gauss_hermite1_r_p;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 1.0, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite1_r_p: %s.", gsl_strerror (ret));

  result = result / (r * r);
  err[0] = err[0] / (r * r);

  return result;
}

/**
 * ncm_integral1d_eval_gauss_laguerre:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = e^{-x}$ over $(0, 1]$.
 *
 * Returns: $\int_0^\infty e^{-x}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_laguerre (NcmIntegral1d *int1d, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, 0.0};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_laguerre;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 1.0, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_laguerre: %s.", gsl_strerror (ret));

  return result;
}

/**
 * ncm_integral1d_eval_gauss_laguerre_r:
 * @int1d: a #NcmIntegral1d
 * @r: the rate $r > 0$
 * @err: (out): the error estimate
 *
 * Integrates in $\alpha = e^{-rx}$ over $(0, 1]$.
 *
 * Returns: $\int_0^\infty e^{-rx}F(x)\,\mathrm{d}x$.
 */
gdouble
ncm_integral1d_eval_gauss_laguerre_r (NcmIntegral1d *int1d, const gdouble r, gdouble *err)
{
  NcmIntegral1dPrivate * const self = ncm_integral1d_get_instance_private (int1d);
  NcIntegral1dHermite int1d_H       = {int1d, 0.0, r};
  gdouble result                    = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);

  F.function = &_ncm_integral1d_eval_gauss_laguerre_r;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 1.0, self->abstol, self->reltol, self->partition, self->rule, self->ws, &result, err);

  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_laguerre_r: %s.", gsl_strerror (ret));

  result = result / r;
  err[0] = err[0] / r;

  return result;
}


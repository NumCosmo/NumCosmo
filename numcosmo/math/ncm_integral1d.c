/***************************************************************************
 *            ncm_integral1d.c
 *
 *  Sat February 20 14:29:30 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_integral1d.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_integral1d
 * @title: NcmIntegral1d
 * @short_description: One dimensional integration object.
 * @stability: Stable
 * @include: numcosmo/math/ncm_integral1d.h
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_integral1d.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmIntegral1dPrivate
{
  guint partition;
  gdouble reltol;
  gdouble abstol;
  guint rule;
  gsl_integration_workspace *ws;
  gsl_integration_cquad_workspace *cquad_ws;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmIntegral1d, ncm_integral1d, G_TYPE_OBJECT);

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
  int1d->priv            = ncm_integral1d_get_instance_private (int1d);
  int1d->priv->partition = 0;
  int1d->priv->rule      = 0;
  int1d->priv->reltol    = 0.0;
  int1d->priv->abstol    = 0.0;
  int1d->priv->ws        = NULL;
  int1d->priv->cquad_ws  = NULL;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_integral1d_finalize (GObject *object)
{
  NcmIntegral1d *int1d = NCM_INTEGRAL1D (object);

  g_clear_pointer (&int1d->priv->ws,       gsl_integration_workspace_free);
  g_clear_pointer (&int1d->priv->cquad_ws, gsl_integration_cquad_workspace_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integral1d_parent_class)->finalize (object);
}

static void
ncm_integral1d_class_init (NcmIntegral1dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_integral1d_set_property;
  object_class->get_property = &ncm_integral1d_get_property;
  object_class->finalize     = &ncm_integral1d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PARTITION,
                                   g_param_spec_uint ("partition",
                                                      NULL,
                                                      "Integral maximum partititon",
                                                      10, G_MAXUINT32, NCM_INTEGRAL1D_DEFAULT_PARTITION,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_RULE,
                                   g_param_spec_uint ("rule",
                                                      NULL,
                                                      "Integration rule",
                                                      1, 6, NCM_INTEGRAL1D_DEFAULT_ALG,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Integral relative tolerance",
                                                        0.0, 1.0, NCM_INTEGRAL1D_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Integral absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_INTEGRAL1D_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * 
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
 * If *@int1d is different from NULL, decreases the reference 
 * count of *@int1d by one and sets *@int1d to NULL.
 * 
 */
void 
ncm_integral1d_clear (NcmIntegral1d **int1d)
{
  g_clear_object (int1d);
}

/**
 * ncm_integral1d_set_partition:
 * @int1d: a #NcmIntegral1d
 * @partition: max number of subintervals
 * 
 * Sets the max number of subintervals to @partition.
 * 
 */
void 
ncm_integral1d_set_partition (NcmIntegral1d *int1d, guint partition)
{
  g_assert_cmpuint (partition, >=, 10);
  if (int1d->priv->partition != partition)
  {
    g_clear_pointer (&int1d->priv->ws, gsl_integration_workspace_free);

    int1d->priv->ws        = gsl_integration_workspace_alloc (partition);
    int1d->priv->partition = partition;
  }
}

/**
 * ncm_integral1d_set_rule:
 * @int1d: a #NcmIntegral1d
 * @rule: Gauss-Kronrod rule
 * 
 * Sets the Gauss-Kronrod @rule to use.
 * 
 */
void 
ncm_integral1d_set_rule (NcmIntegral1d *int1d, guint rule)
{
  int1d->priv->rule = rule;
}

/**
 * ncm_integral1d_set_reltol:
 * @int1d: a #NcmIntegral1d
 * @reltol: relative tolerance
 * 
 * Sets the relative tolerance @reltol to use.
 * 
 */
void 
ncm_integral1d_set_reltol (NcmIntegral1d *int1d, gdouble reltol)
{
  int1d->priv->reltol = reltol;
}

/**
 * ncm_integral1d_set_abstol:
 * @int1d: a #NcmIntegral1d
 * @abstol: absolute tolerance
 * 
 * Sets the absolute tolerance @reltol to use.
 * 
 */
void 
ncm_integral1d_set_abstol (NcmIntegral1d *int1d, gdouble abstol)
{
  int1d->priv->abstol = abstol;
}

/**
 * ncm_integral1d_get_partition:
 * @int1d: a #NcmIntegral1d
 * 
 * Returns: the maximum number of subdivisions used.
 */
guint 
ncm_integral1d_get_partition (NcmIntegral1d *int1d)
{
  return int1d->priv->partition;
}

/**
 * ncm_integral1d_get_rule:
 * @int1d: a #NcmIntegral1d
 * 
 * Returns: the Gauss-Kronrod rule used.
 */
guint 
ncm_integral1d_get_rule (NcmIntegral1d *int1d)
{
  return int1d->priv->rule;
}

/**
 * ncm_integral1d_get_reltol:
 * @int1d: a #NcmIntegral1d
 * 
 * Returns: the relative tolerance used.
 */
gdouble 
ncm_integral1d_get_reltol (NcmIntegral1d *int1d)
{
  return int1d->priv->reltol;
}

/**
 * ncm_integral1d_get_abstol:
 * @int1d: a #NcmIntegral1d
 * 
 * Returns: the absolute tolerance used.
 */
gdouble 
ncm_integral1d_get_abstol (NcmIntegral1d *int1d)
{
  return int1d->priv->abstol;
}

typedef struct _NcIntegral1dHermite
{
  NcmIntegral1d *int1d;
  gdouble mu;
  gdouble r;
  guint neval;
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
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite1_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = sqrt (- 2.0 * log (alpha));

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite_r_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x / int1d_H->r, alpha);
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite1_r_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = sqrt (- 2.0 * log (alpha) / int1d_H->r);

  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);

  return (ncm_integral1d_integrand (int1d_H->int1d, x, alpha) + ncm_integral1d_integrand (int1d_H->int1d, -x, alpha));
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite_mur (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);
  const gdouble y = x / int1d_H->r;

  return (ncm_integral1d_integrand (int1d_H->int1d,  int1d_H->mu + y, alpha) + ncm_integral1d_integrand (int1d_H->int1d, int1d_H->mu - y, alpha));
}

static gdouble 
_ncm_integral1d_eval_gauss_laguerre (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = - log (alpha);

  int1d_H->neval++;
  
  return ncm_integral1d_integrand (int1d_H->int1d, x, alpha);
}

static gdouble 
_ncm_integral1d_eval_gauss_laguerre_r (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = - log (alpha);

  return ncm_integral1d_integrand (int1d_H->int1d, x / int1d_H->r, alpha);
}

/**
 * ncm_integral1d_eval:
 * @int1d: a #NcmIntegral1d
 * @xi: inferior integration limit $x_i$
 * @xf: superior integration limit $x_f$
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $I_F(x_i, x_f) = \int_{x_i}^{x_f}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $I_F(x_i, x_f)$.
 */
gdouble 
ncm_integral1d_eval (NcmIntegral1d *int1d, const gdouble xi, const gdouble xf, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, 0.0};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;
  
  F.function = &_ncm_integral1d_eval_p;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, xi, xf, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval: %s.", gsl_strerror (ret));
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_p:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H^p_F = \int_{0}^{\infty}e^{-x^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H^p_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite_p (NcmIntegral1d *int1d, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, 0.0};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_hermite_p;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite_p: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H_F = \int_{-\infty}^{\infty}e^{-x^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite (NcmIntegral1d *int1d, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, 0.0};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_hermite;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_r_p:
 * @int1d: a #NcmIntegral1d
 * @r: Gaussian scale $r$
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H^p_F = \int_{0}^{\infty}e^{-x^2r^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H^p_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite_r_p (NcmIntegral1d *int1d, const gdouble r, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, r};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);
  
  F.function = &_ncm_integral1d_eval_gauss_hermite_r_p;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite_r_p: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result / r;
  err[0] = ncm_c_sqrt_2pi () * err[0] / r;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_mur:
 * @int1d: a #NcmIntegral1d
 * @r: Gaussian scale $r$
 * @mu: Gaussian mean $\mu$
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H_F = \int_{-\infty}^{\infty}e^{-(x-\mu)^2r^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite_mur (NcmIntegral1d *int1d, const gdouble r, const gdouble mu, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, mu, r};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);

  F.function = &_ncm_integral1d_eval_gauss_hermite_mur;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result / r;
  err[0] = ncm_c_sqrt_2pi () * err[0] / r;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite1_p:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H^p_F = \int_{0}^{\infty}xe^{-x^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H^p_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite1_p (NcmIntegral1d *int1d, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, 0.0};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_hermite1_p;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 1.0, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite1_p: %s.", gsl_strerror (ret));

  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite1_r_p:
 * @int1d: a #NcmIntegral1d
 * @r: Gaussian scale $r$
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H^p_F = \int_{0}^{\infty}xe^{-x^2r^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H^p_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite1_r_p (NcmIntegral1d *int1d, const gdouble r, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, r};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);
  
  F.function = &_ncm_integral1d_eval_gauss_hermite1_r_p;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 1.0, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite1_r_p: %s.", gsl_strerror (ret));

  result = result / r;
  err[0] = err[0] / r;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_laguerre:
 * @int1d: a #NcmIntegral1d
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $L_F = \int_{0}^{\infty}e^{-x}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $L_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_laguerre (NcmIntegral1d *int1d, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, 0.0, 0};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  F.function = &_ncm_integral1d_eval_gauss_laguerre;
  F.params   = &int1d_H;

  ret = gsl_integration_qag (&F, 0.0, 1.0, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_laguerre: %s.", gsl_strerror (ret));
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_laguerre_r:
 * @int1d: a #NcmIntegral1d
 * @r: exponential scale $r$
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $L_F = \int_{0}^{\infty}e^{-xr}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $L_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_laguerre_r (NcmIntegral1d *int1d, const gdouble r, gdouble *err)
{
  NcIntegral1dHermite int1d_H = {int1d, 0.0, r};
  gdouble result = 0.0;
  gsl_function F;
  gint ret;

  g_assert_cmpfloat (r, >, 0.0);

  F.function = &_ncm_integral1d_eval_gauss_laguerre_r;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 1.0, int1d->priv->abstol, int1d->priv->reltol, int1d->priv->partition, int1d->priv->rule, int1d->priv->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_laguerre_r: %s.", gsl_strerror (ret));

  result = result / r;
  err[0] = err[0] / r;
  
  return result;
}

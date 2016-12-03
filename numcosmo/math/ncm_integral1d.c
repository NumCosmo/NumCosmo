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
#include <gsl/gsl_cdf.h>

G_DEFINE_TYPE (NcmIntegral1d, ncm_integral1d, G_TYPE_OBJECT);

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
  int1d->F         = NULL;
  int1d->partition = 0;
  int1d->rule      = 0;
  int1d->reltol    = 0.0;
  int1d->abstol    = 0.0;
  int1d->ws        = NULL;
  int1d->cquad_ws  = NULL;
}

static void
ncm_integral1d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegral1d *int1d = NCM_INTEGRAL1D (object);
  g_return_if_fail (NCM_IS_INTEGRAL1D (object));

  switch (prop_id)
  {
    case PROP_INTEGRAND:
      int1d->F = g_value_get_pointer (value);
      break;
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
    case PROP_INTEGRAND:
      g_value_set_pointer (value, int1d->F);
      break;
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

  g_clear_pointer (&int1d->ws, gsl_integration_workspace_free);
  g_clear_pointer (&int1d->cquad_ws, gsl_integration_cquad_workspace_free);
  
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
                                   PROP_INTEGRAND,
                                   g_param_spec_pointer ("integrand",
                                                         NULL,
                                                         "Integrand function pointer",
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * ncm_integral1d_new:
 * @F: (scope notified): a #NcmIntegral1dF
 * 
 * Creates a new #NcmIntegral1d object for the integrand @F.
 * 
 * Returns: (transfer full): the new #NcmIntegral1d object. 
 */
NcmIntegral1d *
ncm_integral1d_new (NcmIntegral1dF F)
{
  NcmIntegral1d *int1d = g_object_new (NCM_TYPE_INTEGRAL1D,
                                       "integrand", F,
                                       NULL);
  return int1d;
}

/**
 * ncm_integral1d_new_full:
 * @F: (scope notified): a #NcmIntegral1dF
 * @reltol: the relative tolerance
 * @abstol: the absolute tolerance
 * @partition: the maximum subdivisions
 * @rule: integration rule to use in each subinterval
 * 
 * Creates a new #NcmIntegral1d object for the integrand @F.
 * 
 * Returns: (transfer full): the new #NcmIntegral1d object. 
 */
NcmIntegral1d *
ncm_integral1d_new_full (NcmIntegral1dF F, gdouble reltol, gdouble abstol, guint partition, guint rule)
{
  NcmIntegral1d *int1d = g_object_new (NCM_TYPE_INTEGRAL1D,
                                       "integrand", F,
                                       "reltol",    reltol,
                                       "abstol",    abstol,
                                       "partition", partition,
                                       "rule",      rule,
                                       NULL);
  return int1d;
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
  if (int1d->partition != partition)
  {
    g_clear_pointer (&int1d->ws, gsl_integration_workspace_free);
    int1d->ws = gsl_integration_workspace_alloc (partition);
    int1d->partition = partition;
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
  int1d->rule = rule;
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
  int1d->reltol = reltol;
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
  int1d->abstol = abstol;
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
  return int1d->partition;
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
  return int1d->rule;
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
  return int1d->reltol;
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
  return int1d->abstol;
}

/**
 * ncm_integral1d_eval:
 * @int1d: a #NcmIntegral1d
 * @xi: inferior integration limit $x_i$
 * @xf: superior integration limit $x_f$
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $I_F(x_i, x_f) = \int_{x_i}^{x_f}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $I_F(x_i, x_f)$.
 */
gdouble 
ncm_integral1d_eval (NcmIntegral1d *int1d, gdouble xi, gdouble xf, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;

  F.function = int1d->F;
  F.params   = userdata;

  ret = gsl_integration_qag (&F, xi, xf, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval: %s.", gsl_strerror (ret));
  
  return result;
}

typedef struct _NcIntegral1dHermite
{
  gdouble mu;
  gdouble r;
  NcmIntegral1dF F;
  gpointer userdata;
} NcIntegral1dHermite;

static gdouble 
_ncm_integral1d_eval_gauss_hermite_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);

  return int1d_H->F (x, int1d_H->userdata);
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite_r_p (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);

  return int1d_H->F (x / int1d_H->r, int1d_H->userdata);
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);

  return (int1d_H->F (x, int1d_H->userdata) + int1d_H->F (-x, int1d_H->userdata));
}

static gdouble 
_ncm_integral1d_eval_gauss_hermite_mur (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = gsl_cdf_ugaussian_Qinv (alpha);
  const gdouble y = x / int1d_H->r;

  return (int1d_H->F (int1d_H->mu + y, int1d_H->userdata) + int1d_H->F (int1d_H->mu - y, int1d_H->userdata));
}

/**
 * ncm_integral1d_eval_gauss_hermite_p:
 * @int1d: a #NcmIntegral1d
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H^p_F = \int_{0}^{\infty}e^{-x^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H^p_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite_p (NcmIntegral1d *int1d, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;
  NcIntegral1dHermite int1d_H;

  int1d_H.F        = int1d->F;
  int1d_H.userdata = userdata;

  F.function = &_ncm_integral1d_eval_gauss_hermite_p;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite_p: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite:
 * @int1d: a #NcmIntegral1d
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H_F = \int_{-\infty}^{\infty}e^{-x^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite (NcmIntegral1d *int1d, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;
  NcIntegral1dHermite int1d_H;

  int1d_H.F        = int1d->F;
  int1d_H.userdata = userdata;

  F.function = &_ncm_integral1d_eval_gauss_hermite;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_hermite: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result;
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_hermite_r_p:
 * @int1d: a #NcmIntegral1d
 * @r: Gaussian scale $r$
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H^p_F = \int_{0}^{\infty}e^{-x^2r^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H^p_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite_r_p (NcmIntegral1d *int1d, gdouble r, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;
  NcIntegral1dHermite int1d_H;

  g_assert_cmpfloat (r, >, 0.0);
  
  int1d_H.F        = int1d->F;
  int1d_H.userdata = userdata;
  int1d_H.r        = r;

  F.function = &_ncm_integral1d_eval_gauss_hermite_r_p;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
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
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $H_F = \int_{-\infty}^{\infty}e^{-(x-\mu)^2r^2/2}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $H_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_hermite_mur (NcmIntegral1d *int1d, gdouble r, gdouble mu, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;
  NcIntegral1dHermite int1d_H;

  g_assert_cmpfloat (r, >, 0.0);

  int1d_H.F        = int1d->F;
  int1d_H.userdata = userdata;
  int1d_H.r        = r;
  int1d_H.mu       = mu;

  F.function = &_ncm_integral1d_eval_gauss_hermite_mur;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 0.5, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval: %s.", gsl_strerror (ret));

  result = ncm_c_sqrt_2pi () * result / r;
  err[0] = ncm_c_sqrt_2pi () * err[0] / r;
  
  return result;
}

static gdouble 
_ncm_integral1d_eval_gauss_laguerre (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = - log (alpha);

  return int1d_H->F (x, int1d_H->userdata);
}

static gdouble 
_ncm_integral1d_eval_gauss_laguerre_r (gdouble alpha, gpointer userdata)
{
  NcIntegral1dHermite *int1d_H = (NcIntegral1dHermite *) userdata;
  const gdouble x = - log (alpha);

  return int1d_H->F (x / int1d_H->r, int1d_H->userdata);
}

/**
 * ncm_integral1d_eval_gauss_laguerre:
 * @int1d: a #NcmIntegral1d
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $L_F = \int_{0}^{\infty}e^{-x}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $L_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_laguerre (NcmIntegral1d *int1d, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;
  NcIntegral1dHermite int1d_H;

  int1d_H.F        = int1d->F;
  int1d_H.userdata = userdata;

  F.function = &_ncm_integral1d_eval_gauss_laguerre;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 1.0, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_laguerre: %s.", gsl_strerror (ret));
  
  return result;
}

/**
 * ncm_integral1d_eval_gauss_laguerre_r:
 * @int1d: a #NcmIntegral1d
 * @r: exponential scale $r$
 * @userdata: (allow-none): pointer to be passed to the integrand function
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $L_F = \int_{0}^{\infty}e^{-xr}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $L_F$.
 */
gdouble 
ncm_integral1d_eval_gauss_laguerre_r (NcmIntegral1d *int1d, gdouble r, gpointer userdata, gdouble *err)
{
  gsl_function F;
  gdouble result = 0.0;
  gint ret;
  NcIntegral1dHermite int1d_H;

  g_assert_cmpfloat (r, >, 0.0);

  int1d_H.F        = int1d->F;
  int1d_H.userdata = userdata;
  int1d_H.r        = r;

  F.function = &_ncm_integral1d_eval_gauss_laguerre_r;
  F.params   = &int1d_H;
  
  ret = gsl_integration_qag (&F, 0.0, 1.0, int1d->abstol, int1d->reltol, int1d->partition, int1d->rule, int1d->ws, &result, err);
  if (ret != GSL_SUCCESS)
    g_error ("ncm_integral1d_eval_gauss_laguerre_r: %s.", gsl_strerror (ret));

  result = result / r;
  err[0] = err[0] / r;
  
  return result;
}

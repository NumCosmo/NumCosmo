/***************************************************************************
 *            ncm_integrate.c
 *
 *  Wed Aug 13 20:35:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti & Mariana Penna Lima 2012 <vitenti@uel.br>, <pennalima@gmail.com>
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
 * NcmIntegral:
 *
 * Numerical integration helpers.
 *
 * This module provides functions to perform numerical integration. It uses GSL library
 * to perform the integration.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_integration.h>
#include <cuba.h>
#endif /* NUMCOSMO_GIR_SCAN */

static gpointer
_integral_ws_alloc (gpointer userdata)
{
  NCM_UNUSED (userdata);

  return gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);
}

static void
_integral_ws_free (gpointer p)
{
  gsl_integration_workspace_free ((gsl_integration_workspace *) p);
}

/**
 * ncm_integral_get_workspace: (skip)
 *
 * This function provides a workspace to be used by numerical integration
 * functions of GSL. It keeps a internal pool of workspaces and allocate a
 * new one if the function is called and the pool is empty. It is designed
 * to be used in a multi-thread environment. The workspace must be unlocked
 * in order to return to the pool. This must be done using the #ncm_memory_pool_return.
 *
 * Returns: a pointer to #gsl_integration_workspace structure.
 */
gsl_integration_workspace **
ncm_integral_get_workspace ()
{
  G_LOCK_DEFINE_STATIC (create_lock);

  static NcmMemoryPool *mp = NULL;

  G_LOCK (create_lock);

  if (mp == NULL)
    mp = ncm_memory_pool_new (_integral_ws_alloc, NULL, _integral_ws_free);

  G_UNLOCK (create_lock);

  return ncm_memory_pool_get (mp);
}

/**
 * ncm_integral_locked_a_b: (skip)
 * @F: a gsl_function wich is the integrand.
 * @a: lower integration limit.
 * @b: upper integration limit.
 * @abstol: absolute tolerance.
 * @reltol: relative tolerance.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function uses a workspace from the pool and gsl_integration_qag function to perform
 * the numerical integration in the [a, b] interval.
 *
 * Returns: the error code returned by gsl_integration_qag.
 */
gint
ncm_integral_locked_a_b (gsl_function *F, gdouble a, gdouble b, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error)
{
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gint error_code               = gsl_integration_qag (F, a, b, abstol, reltol, NCM_INTEGRAL_PARTITION, 6, *w, result, error);

  ncm_memory_pool_return (w);

  if ((error_code != GSL_SUCCESS) && (error_code != GSL_EROUND))
    *result = GSL_POSINF;

  return error_code;
}

/**
 * ncm_integral_locked_a_inf: (skip)
 * @F: a gsl_function which is the integrand.
 * @a: lower integration limit.
 * @abstol: absolute tolerance.
 * @reltol: relative tolerance.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function uses a workspace from the pool and gsl_integration_qagiu function to perform
 * the numerical integration in the $[a, \infty]$ interval.
 *
 * Returns: the error code returned by gsl_integration_qagiu.
 */
gint
ncm_integral_locked_a_inf (gsl_function *F, gdouble a, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error)
{
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gint error_code               = gsl_integration_qagiu (F, a, abstol, reltol, NCM_INTEGRAL_PARTITION, *w, result, error);

  ncm_memory_pool_return (w);

  if ((error_code != GSL_SUCCESS) && (error_code != GSL_EROUND))
  {
    g_error ("ncm_integral_locked_a_inf: %s", gsl_strerror (error_code));
    *result = GSL_POSINF;
  }

  return error_code;
}

/**
 * ncm_integral_cached_0_x: (skip)
 * @cache: a pointer to #NcmFunctionCache.
 * @F: a gsl_function wich is the integrand.
 * @x: upper integration limit.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function searches for the nearest x_near value previously chosen as the upper integration limit
 * and perform the integration at [x_near, x] interval. This result is summed to that obtained at [0, x_near]
 * and then it is saved in the cache.
 *
 * Returns: the error code returned by gsl_integration_qag.
 */
gint
ncm_integral_cached_0_x (NcmFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error)
{
  gdouble x_found     = 0.0;
  NcmVector *p_result = NULL;
  gint error_code     = GSL_SUCCESS;

  if (ncm_function_cache_get_near (cache, x, &x_found, &p_result, NCM_FUNCTION_CACHE_SEARCH_BOTH))
  {
    if (x == x_found)
    {
      *result = ncm_vector_get (p_result, 0);
    }
    else
    {
      error_code = ncm_integral_locked_a_b (F, x_found, x, 0.0, NCM_INTEGRAL_ERROR, result, error);
      *result   += ncm_vector_get (p_result, 0);
      ncm_function_cache_insert (cache, x, *result);
    }

    ncm_vector_clear (&p_result);
  }
  else
  {
    error_code = ncm_integral_locked_a_b (F, 0.0, x, 0.0, NCM_INTEGRAL_ERROR, result, error);
    ncm_function_cache_insert (cache, x, *result);
  }

  return error_code;
}

/**
 * ncm_integral_cached_x_inf: (skip)
 * @cache: a pointer to #NcmFunctionCache.
 * @F: a gsl_function wich is the integrand.
 * @x: lower integration limit.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function searchs for the nearest x_near value previously chosed as the lower integration limit
 * and perform the integration at $[x, x_{near}]$ interval. This result is summed to that
 * obtained at $[x_{near}, \infty]$ and then it is saved in the cache.
 *
 * Returns: the error code returned by gsl_integration_qagiu.
 */
gint
ncm_integral_cached_x_inf (NcmFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error)
{
  gdouble x_found     = 0.0;
  NcmVector *p_result = NULL;
  gint error_code     = GSL_SUCCESS;

  if (ncm_function_cache_get_near (cache, x, &x_found, &p_result, NCM_FUNCTION_CACHE_SEARCH_BOTH))
  {
    if (x == x_found)
    {
      *result = ncm_vector_get (p_result, 0);
    }
    else
    {
      error_code = ncm_integral_locked_a_b (F, x, x_found, ncm_function_cache_get_abstol (cache), ncm_function_cache_get_reltol (cache), result, error);
      *result   += ncm_vector_get (p_result, 0);

      ncm_function_cache_insert (cache, x, *result);
    }

    ncm_vector_clear (&p_result);
  }
  else
  {
    error_code = ncm_integral_locked_a_inf (F, x, ncm_function_cache_get_abstol (cache), ncm_function_cache_get_reltol (cache), result, error);
    ncm_function_cache_insert (cache, x, *result);
  }

  return error_code;
}

typedef struct _iCLIntegrand2dim
{
  NcmIntegrand2dim *integ;
  gdouble xi;
  gdouble xf;
  gdouble yi;
  gdouble yf;
  gint ldxgiven;
  NcmIntegralPeakfinder p;
} iCLIntegrand2dim;

static gint
_integrand_2dim (const gint *ndim, const gdouble x[], const gint *ncomp, gdouble f[], gpointer userdata)
{
  iCLIntegrand2dim *iinteg = (iCLIntegrand2dim *) userdata;

  NCM_UNUSED (ndim);
  NCM_UNUSED (ncomp);
  f[0] = iinteg->integ->f ((iinteg->xf - iinteg->xi) * x[0] + iinteg->xi, (iinteg->yf - iinteg->yi) * x[1] + iinteg->yi, iinteg->integ->userdata);

  return 0;
}

static void
_peakfinder_2dim (const gint *ndim, const gdouble b[], gint *n, gdouble x[], void *userdata)
{
  iCLIntegrand2dim *iinteg = (iCLIntegrand2dim *) userdata;
  gint i;
  gdouble newb[4];

  newb[0] = iinteg->xi + (iinteg->xf - iinteg->xi) * b[0];
  newb[1] = iinteg->xi + (iinteg->xf - iinteg->xi) * b[1];
  newb[2] = iinteg->yi + (iinteg->yf - iinteg->yi) * b[2];
  newb[3] = iinteg->yi + (iinteg->yf - iinteg->yi) * b[3];

  /*printf ("bounds: %.5g %.5g %.5g %.5g\n", b[0], b[1], b[2], b[3]); */
  /*printf ("new bounds: %.5g %.5g %.5g %.5g\n", newb[0], newb[1], newb[2], newb[3]); */

  iinteg->p (ndim, newb, n, x, iinteg->integ->userdata);

  /*printf ("real minimo: %.15g, %.15g\n", x[0], x[1]); */
  for (i = 0; i < *n; i++)
  {
    x[i * iinteg->ldxgiven + 0] = (x[i * iinteg->ldxgiven + 0] - iinteg->xi) / (iinteg->xf - iinteg->xi);
    x[i * iinteg->ldxgiven + 1] = (x[i * iinteg->ldxgiven + 1] - iinteg->yi) / (iinteg->yf - iinteg->yi);
  }

  /*printf ("minimo: %.15g, %.15g\n", x[0], x[1]); */
}

/**
 * ncm_integrate_2dim:
 * @integ: a pointer to #NcmIntegrand2dim.
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] using the Cuhre algorithm from the Cuba
 * library.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_2dim (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error)
{
  gboolean ret            = FALSE;
  const gint mineval      = 1;
  const gint maxeval      = 10000000;
  const gint key          = 11; /* 13 points rule */
  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, 0, NULL};
  gint nregions, neval, fail;
  gdouble prob;

#ifdef HAVE_LIBCUBA_3_1
  Cuhre (2, 1, &_integrand_2dim, &iinteg, epsrel, epsabs, 0, mineval, maxeval, key, NULL, &nregions, &neval, &fail, result, error, &prob);
#elif defined (HAVE_LIBCUBA_3_3)
  Cuhre (2, 1, &_integrand_2dim, &iinteg, 1, epsrel, epsabs, 0, mineval, maxeval, key, NULL, &nregions, &neval, &fail, result, error, &prob);
#elif defined (HAVE_LIBCUBA_4_0)
  Cuhre (2, 1, &_integrand_2dim, &iinteg, 1, epsrel, epsabs, 0, mineval, maxeval, key, NULL, NULL, &nregions, &neval, &fail, result, error, &prob);
#else
  Cuhre (2, 1, &_integrand_2dim, &iinteg, epsrel, epsabs, 0, mineval, maxeval, key, &nregions, &neval, &fail, result, error, &prob);
#endif /* HAVE_LIBCUBA_3_1 */

  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim: number of evaluations %d >= maximum number of evaluations %d (nregions %d, fail %d, result % 22.15g, error % 22.15g).\n",
               neval, maxeval, nregions, fail, *result, *error);

  *result *= (xf - xi) * (yf - yi);
  *error  *= (xf - xi) * (yf - yi);

  ret = (fail == 0);

  return ret;
}

typedef struct _iCLIntegrand3dim
{
  NcmIntegrand3dim *integ;
  gdouble xi;
  gdouble xf;
  gdouble yi;
  gdouble yf;
  gdouble zi;
  gdouble zf;
  gint ldxgiven;
} iCLIntegrand3dim;

static gint
_integrand_3dim (const gint *ndim, const gdouble x[], const gint *ncomp, gdouble f[], gpointer userdata)
{
  iCLIntegrand3dim *iinteg = (iCLIntegrand3dim *) userdata;

  NCM_UNUSED (ndim);
  NCM_UNUSED (ncomp);
  f[0] = iinteg->integ->f ((iinteg->xf - iinteg->xi) * x[0] + iinteg->xi, (iinteg->yf - iinteg->yi) * x[1] + iinteg->yi, (iinteg->zf - iinteg->zi) * x[2] + iinteg->zi, iinteg->integ->userdata);

  return 0;
}

/**
 * ncm_integrate_3dim:
 * @integ: a pointer to #NcmIntegrand3dim.
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @zi: gbouble which is the lower integration limit of variable z.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @zf: gbouble which is the upper integration limit of variable z.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] and [zi, zf] using the Cuhre algorithm from the Cuba
 * library.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_3dim (NcmIntegrand3dim *integ, gdouble xi, gdouble yi, gdouble zi, gdouble xf, gdouble yf, gdouble zf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error)
{
  gboolean ret            = FALSE;
  const gint mineval      = 1;
  const gint maxeval      = 10000000;
  const gint key          = 11; /* 11 points rule */
  iCLIntegrand3dim iinteg = {integ, xi, xf, yi, yf, zi, zf, 0};
  gint nregions, neval, fail;
  gdouble prob;

#ifdef HAVE_LIBCUBA_3_1
  Cuhre (3, 1, &_integrand_3dim, &iinteg, epsrel, epsabs, 0, mineval, maxeval, key, NULL, &nregions, &neval, &fail, result, error, &prob);
#elif defined (HAVE_LIBCUBA_3_3)
  Cuhre (3, 1, &_integrand_3dim, &iinteg, 1, epsrel, epsabs, 0, mineval, maxeval, key, NULL, &nregions, &neval, &fail, result, error, &prob);
#elif defined (HAVE_LIBCUBA_4_0)
  Cuhre (3, 1, &_integrand_3dim, &iinteg, 1, epsrel, epsabs, 0, mineval, maxeval, key, NULL, NULL, &nregions, &neval, &fail, result, error, &prob);
#else
  Cuhre (3, 1, &_integrand_3dim, &iinteg, epsrel, epsabs, 0, mineval, maxeval, key, &nregions, &neval, &fail, result, error, &prob);
#endif /* HAVE_LIBCUBA_3_1 */

  if (neval >= maxeval)
    g_warning ("ncm_integrate_3dim: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi) * (zf - zi);
  *error  *= (xf - xi) * (yf - yi) * (zf - zi);

  ret = (fail == 0);

  return ret;
}

/**
 * ncm_integrate_2dim_divonne:
 * @integ: a pointer to #NcmIntegrand2dim
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @ngiven: number of peaks
 * @ldxgiven: the leading dimension of xgiven, i.e. the offset between one
 * point and the next in memory (ref. libcuba documentation)
 * @xgiven: list of points where the integrand might have peaks (ref. libcuba documentation)
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] using the Divonne algorithm from the Cuba
 * library.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_2dim_divonne (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, const gint ngiven, const gint ldxgiven, gdouble xgiven[], gdouble *result, gdouble *error)
{
  gboolean ret               = FALSE;
  const gint nvec            = 1;
  const gint seed            = 0;
  const gint mineval         = 1;
  const gint maxeval         = 10000000;
  const gint key1            = 13; /* 13 points rule */
  const gint key2            = 13;
  const gint key3            = 1;
  const gint maxpass         = 10;
  const gdouble border       = 0.0;
  const gdouble maxchisq     = 0.10;
  const gdouble mindeviation = 0.25;
  const gint nextra          = 0;
  peakfinder_t peakfinder    = NULL;
  gint i;

  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, ldxgiven, NULL};
  gint nregions, neval, fail;
  gdouble prob;

  /*printf ("xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]); */
  for (i = 0; i < ngiven; i++)
  {
    xgiven[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
    /*printf ("conv xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]); */
  }

#ifdef HAVE_LIBCUBA_4_0
  Divonne (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border,
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra, peakfinder, NULL, NULL, &nregions, &neval, &fail,
           result, error, &prob);

  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim_divonne: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi);
  *error  *= (xf - xi) * (yf - yi);

  ret = (fail == 0);

  return ret;

#else
  g_error ("ncm_integrate_2dim_divonne: Needs libcuba > 4.0.");

  return FALSE;

#endif /*HAVE_LIBCUBA_4_0 */
}

/**
 * ncm_integrate_2dim_divonne_peakfinder:
 * @integ: a pointer to #NcmIntegrand2dim
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @ngiven: number of peaks
 * @ldxgiven: the leading dimension of xgiven, i.e. the offset between one
 * point and the next in memory (ref. libcuba documentation)
 * @xgiven: list of points where the integrand might have peaks (ref. libcuba documentation)
 * @nextra: the maximum number of extra points the peakfinder will return.
 * @peakfinder: (scope call): a #NcmIntegralPeakfinder
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] using the Divonne algorithm from the Cuba
 * library. It uses a peakfinder to find the peaks of the integrand and improve
 * the integration of concentrated functions.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_2dim_divonne_peakfinder (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, const gint ngiven, const gint ldxgiven, gdouble xgiven[], const gint nextra, NcmIntegralPeakfinder peakfinder, gdouble *result, gdouble *error)
{
  gboolean ret              = FALSE;
  const gint nvec           = 1;
  const gint seed           = 0;
  const gint mineval        = 1;
  const gint maxeval        = 100000000;
  const gint key1           = 13; /* 13 points rule */
  const gint key2           = 13;
  const gint key3           = 1;
  const int maxpass         = 10;
  const double border       = 0.0;
  const double maxchisq     = 0.10;
  const double mindeviation = 0.25;

  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, ldxgiven, peakfinder};
  gint nregions, neval, fail;
  gdouble prob;
  gint i;

  /*printf ("xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]); */
  for (i = 0; i < ngiven; i++)
  {
    xgiven[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
  }

  /*printf ("xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]); */
/*printf ("chamando divonne\n"); */
#ifdef HAVE_LIBCUBA_4_0
  Divonne (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border,
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra, _peakfinder_2dim, NULL, NULL, &nregions, &neval, &fail,
           result, error, &prob);
#else
  g_error ("ncm_integrate_2dim_divonne: Needs libcuba > 4.0.");
#endif /*HAVE_LIBCUBA_4_0 */

  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim_divonne_peakfinder: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi);
  *error  *= (xf - xi) * (yf - yi);

  ret = (fail == 0);

  return ret;
}

/**
 * ncm_integrate_2dim_vegas:
 * @integ: a pointer to #NcmIntegrand2dim
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @nstart: number of samples to start the first round of integration with.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] using the Vegas algorithm from the Cuba
 * library.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_2dim_vegas (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, const gint nstart, gdouble *result, gdouble *error)
{
  gboolean ret         = FALSE;
  const gint nvec      = 1;
  const gint seed      = 0;
  const gint mineval   = 1;
  const gint maxeval   = 10000;
  const gint nincrease = 500;
  const int nbatch     = 1000;
  const int gridno     = 0;

  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, 0, NULL};
  gint neval, fail;
  gdouble prob;

#ifdef HAVE_LIBCUBA_4_0
  Vegas (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval,
         nstart, nincrease, nbatch, gridno, NULL, NULL, &neval, &fail, result, error, &prob);
#else
  g_error ("ncm_integrate_2dim_vegas: Needs libcuba > 4.0.");
#endif /*HAVE_LIBCUBA_4_0 */

  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim_vegas: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi);
  *error  *= (xf - xi) * (yf - yi);

  ret = (fail == 0);

  return ret;
}

/**
 * ncm_integrate_3dim_divonne:
 * @integ: a pointer to #NcmIntegrand3dim
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @zi: gbouble which is the lower integration limit of variable z.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @zf: gbouble which is the upper integration limit of variable z.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @ngiven: number of peaks
 * @ldxgiven: the leading dimension of xgiven, i.e. the offset between one
 * point and the next in memory (ref. libcuba documentation)
 * @xgiven: list of points where the integrand might have peaks (ref. libcuba documentation)
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] and [zi, zf] using the Divonne algorithm from the Cuba
 * library.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_3dim_divonne (NcmIntegrand3dim *integ, gdouble xi, gdouble yi, gdouble zi, gdouble xf, gdouble yf, gdouble zf, gdouble epsrel, gdouble epsabs, const gint ngiven, const gint ldxgiven, gdouble xgiven[], gdouble *result, gdouble *error)
{
  gboolean ret              = FALSE;
  const gint nvec           = 1;
  const gint seed           = 0;
  const gint mineval        = 1; /*1000000; */
  const gint maxeval        = G_MAXINT;
  const gint key1           = 11; /* 11 points rule */
  const gint key2           = 11;
  const gint key3           = 1;
  const int maxpass         = 1;
  const double border       = 0.0;
  const double maxchisq     = 0.10;
  const double mindeviation = 0.25;
  const int nextra          = 0;
  peakfinder_t peakfinder   = NULL;

  iCLIntegrand3dim iinteg = {integ, xi, xf, yi, yf, zi, zf, ldxgiven};
  gint nregions, neval, fail, i;
  gdouble prob;

  for (i = 0; i < ngiven; i++)
  {
    xgiven[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
    xgiven[i * ldxgiven + 2] = (xgiven[i * ldxgiven + 2] - zi) / (zf - zi);
  }

#ifdef HAVE_LIBCUBA_4_0
  Divonne (3, 1, &_integrand_3dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border,
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra, peakfinder, NULL, NULL, &nregions, &neval, &fail,
           result, error, &prob);
#else
  g_error ("ncm_integrate_3dim_divonne: Needs libcuba > 4.0.");
#endif /*HAVE_LIBCUBA_4_0 */

  if (neval >= maxeval)
    g_warning ("ncm_integrate_3dim_divonne: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi) * (zf - zi);
  *error  *= (xf - xi) * (yf - yi) * (zf - zi);

  ret = (fail == 0);

  return ret;
}

/**
 * ncm_integrate_3dim_vegas:
 * @integ: a pointer to #NcmIntegrand3dim
 * @xi: gbouble which is the lower integration limit of variable x.
 * @yi: gbouble which is the lower integration limit of variable y.
 * @zi: gbouble which is the lower integration limit of variable z.
 * @xf: gbouble which is the upper integration limit of variable x.
 * @yf: gbouble which is the upper integration limit of variable y.
 * @zf: gbouble which is the upper integration limit of variable z.
 * @epsrel: relative error
 * @epsabs: absolute error
 * @nstart: number of samples to start the first round of integration with.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function computes the integral of the function @integ->f over the
 * interval [xi, xf] and [yi, yf] and [zi, zf] using the Vegas algorithm from the Cuba
 * library.
 *
 * Returns: whether the integration was successful.
 */
gboolean
ncm_integrate_3dim_vegas (NcmIntegrand3dim *integ, gdouble xi, gdouble yi, gdouble zi, gdouble xf, gdouble yf, gdouble zf, gdouble epsrel, gdouble epsabs, const gint nstart, gdouble *result, gdouble *error)
{
  gboolean ret         = FALSE;
  const gint nvec      = 1;
  const gint seed      = 0;
  const gint mineval   = 1;
  const gint maxeval   = G_MAXINT;
  const gint nincrease = 500;
  const int nbatch     = 1000;
  const int gridno     = 0;

  iCLIntegrand3dim iinteg = {integ, xi, xf, yi, yf, zi, zf, 0};
  gint neval, fail;
  gdouble prob;

#ifdef HAVE_LIBCUBA_4_0
  Vegas (3, 1, &_integrand_3dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval,
         nstart, nincrease, nbatch, gridno, NULL, NULL, &neval, &fail, result, error, &prob);
#else
  g_error ("ncm_integrate_3dim_vegas: Needs libcuba > 4.0.");
#endif /*HAVE_LIBCUBA_4_0 */

  if (neval >= maxeval)
    g_warning ("ncm_integrate_3dim_vegas: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi) * (zf - zi);
  *error  *= (xf - xi) * (yf - yi) * (zf - zi);

  ret = (fail == 0);

  return ret;
}

/**
 * ncm_integral_fixed_new: (skip)
 * @n_nodes: number of nodes in the full interval.
 * @rule_n: order of the Gauss-Legendre integration rule to be applied in each interval.
 * @xl: the interval lower limit.
 * @xu: the interval upper limit.
 *
 * This function prepares the #NcmIntegralFixed with a grid
 * with n_nodes - 1 intervals beteween xl and xu. In each interval it uses
 * a fixed order (rule_n) Gauss-Legendre integration rule to determine the
 * interval inner points. This results in a grid with (n_nodes - 1) * rule_n points.
 *
 * Returns: a pointer to the newly created #NcmIntegralFixed structure.
 */
NcmIntegralFixed *
ncm_integral_fixed_new (gulong n_nodes, gulong rule_n, gdouble xl, gdouble xu)
{
  NcmIntegralFixed *intf = g_slice_new (NcmIntegralFixed);

  intf->n_nodes   = n_nodes;
  intf->rule_n    = rule_n;
  intf->int_nodes = g_slice_alloc (sizeof (gdouble) * n_nodes * rule_n);
  intf->xl        = xl;
  intf->xu        = xu;

  intf->glt = gsl_integration_glfixed_table_alloc (rule_n);

  return intf;
}

/**
 * ncm_integral_fixed_free:
 * @intf: a pointer to #NcmIntegralFixed.
 *
 * This function frees the memory associated to #NcmIntegralFixed.
 *
 */
void
ncm_integral_fixed_free (NcmIntegralFixed *intf)
{
  g_slice_free1 (sizeof (gdouble) * intf->n_nodes * intf->rule_n, intf->int_nodes);
  gsl_integration_glfixed_table_free (intf->glt);
  g_slice_free (NcmIntegralFixed, intf);
}

/**
 * ncm_integral_fixed_calc_nodes: (skip)
 * @intf: a pointer to #NcmIntegralFixed.
 * @F: a pointer to a gsl_function.
 *
 * This function calculates the nodes of the #NcmIntegralFixed.
 * It uses the Gauss-Legendre integration rule to determine the
 * interval inner points.
 *
 */
void
ncm_integral_fixed_calc_nodes (NcmIntegralFixed *intf, gsl_function *F)
{
  const gulong r2         = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x   = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gulong i, j, k = 0;

  if (odd_rule)
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;

      for (j = 1; j < r2 + 1; j++)
        intf->int_nodes[k++] = GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->glt->w[j];

      intf->int_nodes[k++] = GSL_FN_EVAL (F, x1px0_2) * intf->glt->w[0];

      for (j = 1; j < r2 + 1; j++)
        intf->int_nodes[k++] = GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->glt->w[j];
    }
  }
  else
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;

      for (j = 0; j < r2; j++)
        intf->int_nodes[k++] = GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->glt->w[j];

      for (j = 0; j < r2; j++)
        intf->int_nodes[k++] = GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->glt->w[j];
    }
  }
}

/**
 * ncm_integral_fixed_nodes_eval:
 * @intf: a pointer to #NcmIntegralFixed.
 *
 * This function evaluates the integral of the function @integ->f over the
 * interval [xi, xf] using the nodes calculated by #ncm_integral_fixed_calc_nodes.
 *
 * Returns: the integral of the function @integ->f over the interval [xi, xf].
 */
gdouble
ncm_integral_fixed_nodes_eval (NcmIntegralFixed *intf)
{
  glong i;
  glong maxi            = (intf->n_nodes - 1) * intf->rule_n;
  gdouble res           = 0.0;
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);

  for (i = 0; i < maxi; i++)
    res += intf->int_nodes[i];

  return res * delta_x * 0.5;
}

/**
 * ncm_integral_fixed_integ_mult: (skip)
 * @intf: a pointer to #NcmIntegralFixed.
 * @F: a pointer to gsl_function.
 *
 * This function evaluates the integral of the function @integ->f over the
 * interval [xi, xf] using the nodes calculated by #ncm_integral_fixed_calc_nodes.
 * It uses the Gauss-Legendre integration rule to determine the
 * interval inner points. This function multiplies the integrand by the
 * function @F.
 *
 * Returns: the integral of the function @integ->f times @F over the interval [xi, xf].
 */
gdouble
ncm_integral_fixed_integ_mult (NcmIntegralFixed *intf, gsl_function *F)
{
  const gulong r2         = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x   = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gdouble res             = 0.0;
  gulong i, j, k = 0;

  if (odd_rule)
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;

      for (j = 1; j < r2 + 1; j++)
        res += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      res += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];

      for (j = 1; j < r2 + 1; j++)
        res += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
    }
  }
  else
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;

      for (j = 0; j < r2; j++)
        res += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      for (j = 0; j < r2; j++)
        res += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
    }
  }

  return res * delta_x * 0.5;
}

/**
 * ncm_integral_fixed_integ_posdef_mult: (skip)
 * @intf: a pointer to #NcmIntegralFixed
 * @F: a pointer to gsl_function
 * @max: maximum value of the integration interval
 * @reltol: relative tolerance
 *
 * This function computes the integral of the function @integ->f over the
 * interval starting at @max and going to the left using the nodes calculated
 * by #ncm_integral_fixed_calc_nodes. It uses the Gauss-Legendre integration
 * rule to determine the interval inner points. This function multiplies the
 * integrand by the function @F. It stops when the relative error is less than
 * @reltol.
 *
 * Returns: the integral of the function @integ->f times @F over the interval [max, xf].
 */
gdouble
ncm_integral_fixed_integ_posdef_mult (NcmIntegralFixed *intf, gsl_function *F, gdouble max, gdouble reltol)
{
  const gulong r2         = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x   = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  const glong mnode       = max / delta_x;
  gdouble res             = 0.0;
  glong i, j, k = 0;

  g_assert (mnode < (glong) intf->n_nodes);

  if (odd_rule)
  {
    for (i = mnode; i < (glong) (intf->n_nodes - 1); i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part          = 0.0;

      k = i * intf->rule_n;

      for (j = 1; j < (glong) (r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      part += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];

      for (j = 1; j < (glong) (r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      res += part;

      if (fabs (part / res) < reltol)
        break;
    }

    for (i = mnode - 1; i >= 0; i--)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part          = 0.0;

      k = i * intf->rule_n;

      for (j = 1; j < (glong) (r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      part += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];

      for (j = 1; j < (glong) (r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      res += part;

      if (fabs (part / res) < reltol)
        break;
    }
  }
  else
  {
    for (i = mnode; i < (glong) (intf->n_nodes - 1); i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part          = 0.0;

      k = i * intf->rule_n;

      for (j = 0; j < (glong) r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      for (j = 0; j < (glong) r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      res += part;

      if (fabs (part / res) < reltol)
        break;
    }

    for (i = mnode - 1; i >= 0; i--)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part          = 0.0;

      k = i * intf->rule_n;

      for (j = 0; j < (glong) r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      for (j = 0; j < (glong) r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];

      res += part;

      if (fabs (part / res) < reltol)
        break;
    }
  }

  return res * delta_x * 0.5;
}

/**
 * ncm_integral_fixed_integ_vec_mult:
 * @intf: a pointer to #NcmIntegralFixed whose nodes have been populated by
 *   #ncm_integral_fixed_calc_nodes.
 * @f_at_nodes: a #NcmVector of length `(n_nodes - 1) * rule_n` holding the
 *   integrand values evaluated at the abscissae returned by
 *   #ncm_integral_fixed_get_nodes.
 *
 * Computes the integral $\int F(x) G(x) \mathrm{d}x$ as the dot product of the
 * stored quadrature-weighted nodes (which already contain $w_k F(x_k)$) with
 * @f_at_nodes (the externally evaluated $G(x_k)$), scaled by `delta_x / 2`.
 *
 * Returns: the value of the integral.
 */
gdouble
ncm_integral_fixed_integ_vec_mult (NcmIntegralFixed *intf, const NcmVector *f_at_nodes)
{
  const gulong total_n  = (intf->n_nodes - 1) * intf->rule_n;
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gdouble res           = 0.0;
  gulong i;

  g_assert_nonnull (f_at_nodes);
  g_assert_cmpuint (ncm_vector_len (f_at_nodes), ==, total_n);

  for (i = 0; i < total_n; i++)
    res += intf->int_nodes[i] * ncm_vector_get (f_at_nodes, i);

  return res * delta_x * 0.5;
}

/**
 * ncm_integral_fixed_get_nodes:
 * @intf: a pointer to #NcmIntegralFixed.
 * @nodes (out): a #NcmVector of length `(n_nodes - 1) * rule_n` to receive the node abscissae.
 *
 * Fills @nodes with the canonical node abscissae used by
 * #ncm_integral_fixed_calc_nodes, in the same iteration order.
 */
void
ncm_integral_fixed_get_nodes (NcmIntegralFixed *intf, NcmVector *nodes)
{
  const gulong r2         = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x   = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  const gulong total_n    = (intf->n_nodes - 1) * intf->rule_n;
  gulong i, j, k = 0;

  g_assert_nonnull (nodes);
  g_assert_cmpuint (ncm_vector_len (nodes), ==, total_n);

  if (odd_rule)
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;

      for (j = 1; j < r2 + 1; j++)
        ncm_vector_set (nodes, k++, x1px0_2 - x1mx0_2 * intf->glt->x[j]);

      ncm_vector_set (nodes, k++, x1px0_2);

      for (j = 1; j < r2 + 1; j++)
        ncm_vector_set (nodes, k++, x1px0_2 + x1mx0_2 * intf->glt->x[j]);
    }
  }
  else
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0      = intf->xl + delta_x * i;
      const gdouble x1      = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;

      for (j = 0; j < r2; j++)
        ncm_vector_set (nodes, k++, x1px0_2 - x1mx0_2 * intf->glt->x[j]);

      for (j = 0; j < r2; j++)
        ncm_vector_set (nodes, k++, x1px0_2 + x1mx0_2 * intf->glt->x[j]);
    }
  }
}

/* Builds the fixed rule (@n_nodes, @rule_n) over [@xl, @xu], bakes @F, and tests
 * whether its INT F*G matches @I_ref and its INT F matches @guard_mass, both to
 * @reltol. The F*G relative error is returned in @err_I_out for fallback ranking. */
static gboolean
_ncm_integral_fixed_calib_try (gsl_function *F, gsl_function *G, gdouble xl, gdouble xu, gulong n_nodes, gulong rule_n, gdouble reltol, gdouble I_ref, gdouble denom_I, gdouble guard_mass, gdouble *err_I_out)
{
  NcmIntegralFixed *trial = ncm_integral_fixed_new (n_nodes, rule_n, xl, xu);
  const gdouble denom_m   = (guard_mass != 0.0) ? fabs (guard_mass) : 1.0;
  gdouble I_trial, mass, err_I;
  gboolean ok;

  ncm_integral_fixed_calc_nodes (trial, F);
  I_trial = ncm_integral_fixed_integ_mult (trial, G);
  mass    = ncm_integral_fixed_nodes_eval (trial);
  ncm_integral_fixed_free (trial);

  err_I = fabs (I_trial - I_ref) / denom_I;
  ok    = (err_I < reltol) && (fabs (mass - guard_mass) / denom_m < reltol);

  if (err_I_out != NULL)
    *err_I_out = err_I;

  return ok;
}

/**
 * ncm_integral_fixed_calibrate: (skip)
 * @F: a pointer to a gsl_function, the integration weight (baked into the nodes).
 * @G: a pointer to a gsl_function, the slowly-varying factor probed at the nodes.
 * @xl: the interval lower limit.
 * @xu: the interval upper limit.
 * @reltol: target relative tolerance.
 * @exact_F_integ: the exact value of $\int F \mathrm{d}x$ over [@xl, @xu] for the
 *   missed-mass guard, or %GSL_NAN to fall back to the reference's own
 *   $\int F \mathrm{d}x$.
 * @max_total_nodes: ceiling on the total node count $(n_\mathrm{nodes}-1)\,r$.
 * @n_nodes_out: (out): receives the selected number of nodes.
 * @rule_n_out: (out): receives the selected Gauss-Legendre rule order.
 *
 * Searches the $(n_\mathrm{nodes}, \mathrm{rule}_n)$ configuration space for the
 * fixed Gauss-Legendre rule of minimal *total* node count whose estimate of
 * $\int F(x) G(x) \mathrm{d}x$ over [@xl, @xu] matches a high-resolution
 * reference (n_nodes = 1000, rule_n = 7) to within @reltol. The selected
 * configuration must additionally reproduce $\int F \mathrm{d}x$ (via
 * #ncm_integral_fixed_nodes_eval) to within @reltol — this catches an @F feature
 * narrower than the reference panel width, where the $F\,G$ test alone could
 * false-pass because both reference and trial straddle it. The guard baseline is
 * @exact_F_integ when finite, otherwise the reference's own $\int F$.
 *
 * For each candidate @rule_n the threshold $n_\mathrm{nodes}$ is bracketed by
 * geometric growth and then pinned exactly by bisection (assuming convergence is
 * monotone in $n_\mathrm{nodes}$), so the returned configuration is the true
 * minimal-total rule rather than a geometric overshoot.
 *
 * If no configuration within @max_total_nodes meets the tolerance, the
 * configuration with the smallest $F\,G$ error seen is used and a warning is
 * emitted.
 *
 * Returns: (transfer full): a newly allocated #NcmIntegralFixed at the selected
 * configuration, with @F already baked into its nodes (see
 * #ncm_integral_fixed_calc_nodes).
 */
NcmIntegralFixed *
ncm_integral_fixed_calibrate (gsl_function *F, gsl_function *G, gdouble xl, gdouble xu, gdouble reltol, gdouble exact_F_integ, gulong max_total_nodes, guint *n_nodes_out, guint *rule_n_out)
{
  /* Reference resolution scaled to the target tolerance. A fixed Gauss-Legendre
   * rule resolves a piecewise-smooth weight at ~4th order, so panels growing as
   * reltol^(-1/4) keep the reference comfortably more accurate than @reltol while
   * staying just fine enough to exceed the selected configs (which themselves
   * grow as the tolerance tightens) - avoiding the cost of a fixed
   * ultra-high-resolution grid for loose tolerances. Clamped to a sane range. */
  const gulong ref_n_nodes      = (gulong) CLAMP ((glong) (10.0 * pow (reltol, -0.3) + 0.5), 96, 2048);
  const gulong ref_rule_n       = 7;
  const guint rule_candidates[] = { 3, 5, 7 };
  gdouble I_ref, denom_I, guard_mass;
  gulong best_n_nodes = 0, best_rule_n = 0, best_total = 0;
  gulong fb_n_nodes = 0, fb_rule_n = 0;
  gdouble fb_err     = GSL_POSINF;
  gboolean converged = FALSE;
  guint ri;

  g_assert (xu > xl);
  g_assert_cmpfloat (reltol, >, 0.0);

  /* High-resolution reference: INT F*G and (for the missed-mass guard) INT F.
   * The guard baseline is the caller's exact INT F when finite, else the
   * reference's own INT F - which still flags any trial that drops mass the
   * reference resolves, without requiring an (abort-prone) adaptive integral. */
  {
    NcmIntegralFixed *ref = ncm_integral_fixed_new (ref_n_nodes, ref_rule_n, xl, xu);

    ncm_integral_fixed_calc_nodes (ref, F);
    I_ref      = ncm_integral_fixed_integ_mult (ref, G);
    guard_mass = gsl_finite (exact_F_integ) ? exact_F_integ : ncm_integral_fixed_nodes_eval (ref);
    ncm_integral_fixed_free (ref);
  }

  denom_I = (I_ref != 0.0) ? fabs (I_ref) : 1.0;

  for (ri = 0; ri < G_N_ELEMENTS (rule_candidates); ri++)
  {
    const gulong rule_n = rule_candidates[ri];

    /* Largest n_nodes worth testing for this rule: its total must not exceed the
     * cap, nor (once a candidate exists) the best total already found - so only a
     * strictly smaller total can win. (n_ceil - 1) * rule_n <= cap_total. */
    const gulong cap_total = (best_total != 0) ? MIN (best_total - 1, max_total_nodes) : max_total_nodes;
    const gulong n_ceil    = cap_total / rule_n + 1;
    gulong n_nodes         = 2;
    gulong last_fail       = 0; /* largest n_nodes seen to fail (0 = none) */
    gulong hi              = 0; /* a passing n_nodes that brackets the threshold (0 = none) */
    gdouble err_I;

    if (n_ceil < 2)
      continue;  /* even the coarsest rule of this order cannot beat the best */

    /* Geometric phase: grow n_nodes (x1.5) only to BRACKET the pass threshold,
     * clamped to n_ceil. Convergence is assumed monotone in n_nodes. */
    while (TRUE)
    {
      const gboolean at_ceiling = (n_nodes >= n_ceil);
      const gulong n_try        = at_ceiling ? n_ceil : n_nodes;
      const gboolean ok         = _ncm_integral_fixed_calib_try (F, G, xl, xu, n_try, rule_n, reltol, I_ref, denom_I, guard_mass, &err_I);

      if (err_I < fb_err)
      {
        fb_err     = err_I;
        fb_n_nodes = n_try;
        fb_rule_n  = rule_n;
      }

      if (ok)
      {
        hi = n_try;
        break;
      }

      last_fail = n_try;

      if (at_ceiling)
        break;  /* ceiling fails => no passing config for this rule in range */

      {
        const gulong next = (gulong) ceil (n_nodes * 1.5);

        n_nodes = (next > n_nodes) ? next : (n_nodes + 1);
      }
    }

    if (hi == 0)
      continue;

    /* Bisection: exact smallest passing n_nodes in (last_fail, hi]. */
    {
      gulong lo = (last_fail > 0) ? last_fail : 1;

      while (hi - lo > 1)
      {
        const gulong mid = lo + (hi - lo) / 2;

        if (_ncm_integral_fixed_calib_try (F, G, xl, xu, mid, rule_n, reltol, I_ref, denom_I, guard_mass, &err_I))
          hi = mid;
        else
          lo = mid;

        if (err_I < fb_err)
        {
          fb_err     = err_I;
          fb_n_nodes = mid;
          fb_rule_n  = rule_n;
        }
      }
    }

    {
      const gulong total = (hi - 1) * rule_n;

      if ((best_total == 0) || (total < best_total))
      {
        best_n_nodes = hi;
        best_rule_n  = rule_n;
        best_total   = total;
        converged    = TRUE;
      }
    }
  }

  if (!converged)
  {
    best_n_nodes = fb_n_nodes;
    best_rule_n  = fb_rule_n;
    g_warning ("ncm_integral_fixed_calibrate: did not reach reltol %g within %lu total nodes "
               "(best rel error %g at n_nodes=%lu, rule_n=%lu).",
               reltol, max_total_nodes, fb_err, best_n_nodes, best_rule_n);
  }

  if (n_nodes_out != NULL)
    *n_nodes_out = (guint) best_n_nodes;

  if (rule_n_out != NULL)
    *rule_n_out = (guint) best_rule_n;

  {
    NcmIntegralFixed *intf = ncm_integral_fixed_new (best_n_nodes, best_rule_n, xl, xu);

    ncm_integral_fixed_calc_nodes (intf, F);

    return intf;
  }
}

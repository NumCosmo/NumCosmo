/***************************************************************************
 *            integral.c
 *
 *  Wed Aug 13 20:35:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti & Mariana Penna Lima 2012 <sandro@lapsandro>, <pennalima@gmail.com>
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
 * SECTION:integral
 * @title: NcmIntegral
 * @short_description: Numerical integration helpers.
 * @stability: Stable
 * @include: numcosmo/math/integral.h
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_util.h"

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
  gsl_integration_workspace_free ((gsl_integration_workspace *)p);
}

/**
 * ncm_integral_get_workspace: (skip)
 *
 * This function provides a workspace to be used by numerical integration
 * functions of GSL. It keeps a internal pool of workspaces and allocate a
 * new one if the function is called and the pool is empty. It is designed
 * to be used in a multithread enviroment. The workspace must be unlocked
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
  gsl_integration_workspace **w = ncm_integral_get_workspace();
  gint error_code = gsl_integration_qag (F, a, b, abstol, reltol, NCM_INTEGRAL_PARTITION, 6, *w, result, error);

  ncm_memory_pool_return (w);
  if (error_code != GSL_SUCCESS && error_code != GSL_EROUND)
    *result = GSL_POSINF;
  return error_code;
}

/**
 * ncm_integral_locked_a_inf: (skip)
 * @F: a gsl_function wich is the integrand.
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
  
  if (error_code != GSL_SUCCESS && error_code != GSL_EROUND)
  {
    g_warning ("ncm_integral_locked_a_inf: %s", gsl_strerror (error_code));
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
 * This function searchs for the nearest x_near value previously chosed as the upper integration limit
 * and perform the integration at [x_near, x] interval. This result is summed to that obtained at [0, x_near]
 * and then it is saved in the cache.
 *
 * Returns: the error code returned by gsl_integration_qag.
 */
gint
ncm_integral_cached_0_x (NcmFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error)
{
  gdouble x_found = 0.0;
  NcmVector *p_result;
  gint error_code = GSL_SUCCESS;

//printf ("[%p]SEARCH! -> %g\n", g_thread_self (), x);
  if (ncm_function_cache_get_near (cache, x, &x_found, &p_result, NCM_FUNCTION_CACHE_SEARCH_BOTH))
  {
//printf ("[%p]Found out %g %g [%p]\n", g_thread_self (), x_found, ncm_vector_get (p_result, 0), p_result);
    if (x == x_found)
      *result = ncm_vector_get (p_result, 0);
    else
    {
      error_code = ncm_integral_locked_a_b (F, x_found, x, 0.0, NCM_INTEGRAL_ERROR, result, error);
      *result += ncm_vector_get (p_result, 0);
      ncm_function_cache_insert (cache, x, *result);
    }
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
  gdouble x_found = 0.0;
  NcmVector *p_result;
  gint error_code = GSL_SUCCESS;

  if (ncm_function_cache_get_near (cache, x, &x_found, &p_result, NCM_FUNCTION_CACHE_SEARCH_BOTH))
  {
    if (x == x_found)
    {
      *result = ncm_vector_get (p_result, 0);
    }
    else
    {
      error_code = ncm_integral_locked_a_b (F, x, x_found, ncm_function_cache_get_abstol (cache), ncm_function_cache_get_reltol (cache), result, error);
      *result += ncm_vector_get(p_result, 0);

      ncm_function_cache_insert (cache, x, *result);
    }
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

  //printf ("bounds: %.5g %.5g %.5g %.5g\n", b[0], b[1], b[2], b[3]);
  //printf ("new bounds: %.5g %.5g %.5g %.5g\n", newb[0], newb[1], newb[2], newb[3]);

  iinteg->p (ndim, newb, n, x, iinteg->integ->userdata);

  //printf ("real minimo: %.15g, %.15g\n", x[0], x[1]);
  for (i = 0; i < *n; i++)
  {
    x[i * iinteg->ldxgiven + 0] = (x[i * iinteg->ldxgiven + 0] - iinteg->xi) / (iinteg->xf - iinteg->xi);
    x[i * iinteg->ldxgiven + 1] = (x[i * iinteg->ldxgiven + 1] - iinteg->yi) / (iinteg->yf - iinteg->yi);
  }
  //printf ("minimo: %.15g, %.15g\n", x[0], x[1]);
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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_2dim (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
	const gint mineval = 1;
	const gint maxeval = 10000000;
	const gint key = 13; /* 13 points rule */
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
    g_warning ("ncm_integrate_2dim: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);
  
	*result *= (xf - xi) * (yf - yi);
	*error *= (xf - xi) * (yf - yi);

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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_3dim (NcmIntegrand3dim *integ, gdouble xi, gdouble yi, gdouble zi, gdouble xf, gdouble yf, gdouble zf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
	const gint mineval = 1;
	const gint maxeval = 10000000;
	const gint key = 11; /* 11 points rule */
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
	*error *= (xf - xi) * (yf - yi) * (zf - zi);

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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_2dim_divonne (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, const gint ngiven, const gint ldxgiven, gdouble xgiven[], gdouble *result, gdouble *error)
{
  gboolean ret               = FALSE;
  const gint nvec            = 1;
  const gint seed            = 0;
	const gint mineval         = 1;
	const gint maxeval         = 10000000;
	const gint key1            = 13; // 13 points rule 
  const gint key2            = 13;
  const gint key3            = 1;
  const gint maxpass         = 10;
  const gdouble border       = 0.0;
  const gdouble maxchisq     = 0.10;
  const gdouble mindeviation = 0.25;
  const gint nextra = 0;
  peakfinder_t peakfinder = NULL;
  guint i;
  
  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, ldxgiven, NULL};
	gint nregions, neval, fail;
	gdouble prob;

  //printf ("xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]);
  for (i = 0; i < ngiven; i++)
  {
    xgiven[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
    //printf ("conv xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]);
  }
  
#ifdef HAVE_LIBCUBA_4_0
	Divonne (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border, 
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra, peakfinder, NULL, NULL, &nregions, &neval, &fail, 
           result, error, &prob);  
  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim_divonne: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);
	*result *= (xf - xi) * (yf - yi);
	*error *= (xf - xi) * (yf - yi);

	ret = (fail == 0);
	return ret;
#else
  g_error ("ncm_integrate_2dim_divonne: Needs libcuba > 4.0.");
  return FALSE;
#endif //HAVE_LIBCUBA_4_0      
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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_2dim_divonne_peakfinder (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, const gint ngiven, const gint ldxgiven, gdouble xgiven[], const gint nextra, NcmIntegralPeakfinder peakfinder, gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
  const gint nvec = 1;
  const gint seed = 0;
	const gint mineval = 1;
	const gint maxeval = 100000000;
	const gint key1 = 13; // 13 points rule 
  const gint key2 = 13;
  const gint key3 = 1;
  const int maxpass = 10;
  const double border = 0.0;
  const double maxchisq = 0.10;
  const double mindeviation = 0.25;
  
  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, ldxgiven, peakfinder};
	gint nregions, neval, fail;
	gdouble prob;
  guint i;

  //printf ("xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]);
  for (i = 0; i < ngiven; i++)
  {
    xgiven[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
  }
  //printf ("xgiven: %.20g %.20g\n", xgiven[0], xgiven[1]);
//printf ("chamando divonne\n");
#ifdef HAVE_LIBCUBA_4_0
	Divonne (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border, 
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra, _peakfinder_2dim, NULL, NULL, &nregions, &neval, &fail, 
           result, error, &prob);  
#else
  g_error ("ncm_integrate_2dim_divonne: Needs libcuba > 4.0.");
#endif //HAVE_LIBCUBA_4_0  

  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim_divonne_peakfinder: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);
  
	*result *= (xf - xi) * (yf - yi);
	*error *= (xf - xi) * (yf - yi);

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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_2dim_vegas (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, const gint nstart, gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
  const gint nvec = 1;
  const gint seed = 0;
	const gint mineval = 1;
	const gint maxeval = 10000;
	const gint nincrease = 500;
  const int nbatch = 1000;
  const int gridno = 0;
  
  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, 0, NULL};
	gint neval, fail;
	gdouble prob;
  
#ifdef HAVE_LIBCUBA_4_0
	Vegas (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, 
         nstart, nincrease, nbatch, gridno, NULL, NULL, &neval, &fail, result, error, &prob);  
#else
  g_error ("ncm_integrate_2dim_vegas: Needs libcuba > 4.0.");
#endif //HAVE_LIBCUBA_4_0  

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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_3dim_divonne (NcmIntegrand3dim *integ, gdouble xi, gdouble yi, gdouble zi, gdouble xf, gdouble yf, gdouble zf, gdouble epsrel, gdouble epsabs, const gint ngiven, const gint ldxgiven, gdouble xgiven[], gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
  const gint nvec = 1;
  const gint seed = 0;
	const gint mineval = 1; //1000000;
	const gint maxeval = G_MAXINT;
	const gint key1 = 11; // 11 points rule 
  const gint key2 = 11;
  const gint key3 = 1;
  const int maxpass = 1;
  const double border = 0.0;
  const double maxchisq = 0.10;
  const double mindeviation = 0.25;
  const int nextra = 0;
  peakfinder_t peakfinder = NULL;
  
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
#endif //HAVE_LIBCUBA_4_0  

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
 * This function FIXME
 *
 * Returns: a gboolean
 */
gboolean
ncm_integrate_3dim_vegas (NcmIntegrand3dim *integ, gdouble xi, gdouble yi, gdouble zi, gdouble xf, gdouble yf, gdouble zf, gdouble epsrel, gdouble epsabs, const gint nstart, gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
  const gint nvec = 1;
  const gint seed = 0;
	const gint mineval = 1;
	const gint maxeval = G_MAXINT;
	const gint nincrease = 500;
  const int nbatch = 1000;
  const int gridno = 0;
  
  iCLIntegrand3dim iinteg = {integ, xi, xf, yi, yf, zi, zf, 0};
	gint neval, fail;
	gdouble prob;
  
#ifdef HAVE_LIBCUBA_4_0
	Vegas (3, 1, &_integrand_3dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, 
         nstart, nincrease, nbatch, gridno, NULL, NULL, &neval, &fail, result, error, &prob);  
#else
  g_error ("ncm_integrate_3dim_vegas: Needs libcuba > 4.0.");
#endif //HAVE_LIBCUBA_4_0  

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
  intf->n_nodes = n_nodes;
  intf->rule_n = rule_n;
  intf->int_nodes = g_slice_alloc (sizeof(gdouble) * n_nodes * rule_n);
  intf->xl = xl;
  intf->xu = xu;

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
  g_slice_free1 (sizeof(gdouble) * intf->n_nodes * intf->rule_n, intf->int_nodes);
  gsl_integration_glfixed_table_free (intf->glt);
  g_slice_free (NcmIntegralFixed, intf);
}

/**
 * ncm_integral_fixed_calc_nodes: (skip)
 * @intf: a pointer to #NcmIntegralFixed.
 * @F: a pointer to a gsl_function.
 *
 * This function FIXME
 *
*/
void
ncm_integral_fixed_calc_nodes (NcmIntegralFixed *intf, gsl_function *F)
{
  const gulong r2 = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gulong i, j, k = 0;

  if (odd_rule)
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
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
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
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
 * This function
 *
 * Returns: FIXME
*/
gdouble
ncm_integral_fixed_nodes_eval (NcmIntegralFixed *intf)
{
  glong i;
  glong maxi = (intf->n_nodes - 1) * intf->rule_n;
  gdouble res = 0.0;
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
 * This function
 *
 * Returns: FIXME
*/
gdouble
ncm_integral_fixed_integ_mult (NcmIntegralFixed *intf, gsl_function *F)
{
  const gulong r2 = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gdouble res = 0.0;
  gulong i, j, k = 0;

  if (odd_rule)
  {
    for (i = 0; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
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
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
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
 * @intf: a pointer to #NcmIntegralFixed.
 * @F: a pointer to gsl_function.
 * @max: FIXME
 * @reltol: FIXME
 *
 * This function
 *
 * Returns: FIXME
*/
gdouble
ncm_integral_fixed_integ_posdef_mult (NcmIntegralFixed *intf, gsl_function *F, gdouble max, gdouble reltol)
{
  const gulong r2 = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  const glong mnode = max / delta_x;
  gdouble res = 0.0;
  glong i, j, k = 0;

  g_assert (mnode < (glong) intf->n_nodes);

  if (odd_rule)
  {
    for (i = mnode; i < (glong) (intf->n_nodes - 1); i++)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part = 0.0;
      k = i * intf->rule_n;

      for (j = 1; j < (glong) (r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      part += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];
      for (j = 1; j < (glong) (r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      res += part;

      if (fabs(part / res) < reltol)
        break;
    }

    for (i = mnode - 1; i >= 0; i--)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part = 0.0;
      k = i * intf->rule_n;

      for (j = 1; j < (glong)(r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      part += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];
      for (j = 1; j < (glong)(r2 + 1); j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      res += part;
      if (fabs(part / res) < reltol)
        break;
    }
  }
  else
  {
    for (i = mnode; i < (glong) (intf->n_nodes - 1); i++)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part = 0.0;
      k = i * intf->rule_n;

      for (j = 0; j < (glong)r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      for (j = 0; j < (glong)r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      res += part;

      if (fabs(part / res) < reltol)
        break;
    }

    for (i = mnode - 1; i >= 0; i--)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part = 0.0;
      k = i * intf->rule_n;

      for (j = 0; j < (glong)r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      for (j = 0; j < (glong)r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      res += part;

      if (fabs(part / res) < reltol)
        break;
    }
  }

  return res * delta_x * 0.5;
}

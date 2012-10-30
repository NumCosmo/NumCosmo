/***************************************************************************
 *            integral.c
 *
 *  Wed Aug 13 20:35:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * @title: Numerical Integration
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/integral.h"
#include "math/memory_pool.h"

#include <gsl/gsl_integration.h>
#ifdef HAVE_LIBCUBA
#include <cuba.h>  
#endif /* HAVE_LIBCUBA */


static gpointer
_integral_ws_alloc (void)
{
  return gsl_integration_workspace_alloc (NC_INT_PARTITION);
}

static void
_integral_ws_free (gpointer p)
{
  gsl_integration_workspace_free ((gsl_integration_workspace *)p);
}

/**
 * nc_integral_get_workspace: (skip)
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
nc_integral_get_workspace ()
{
  static GStaticMutex create_lock = G_STATIC_MUTEX_INIT;
  static NcmMemoryPool *mp = NULL;

  g_static_mutex_lock (&create_lock);
  if (mp == NULL)
    mp = ncm_memory_pool_new (_integral_ws_alloc, _integral_ws_free);
  g_static_mutex_unlock (&create_lock);

  return ncm_memory_pool_get (mp);
}

/**
 * nc_integral_locked_a_b: (skip)
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
nc_integral_locked_a_b (gsl_function *F, gdouble a, gdouble b, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error)
{
  gsl_integration_workspace **w = nc_integral_get_workspace();
  gint error_code = gsl_integration_qag (F, a, b, abstol, reltol, NC_INT_PARTITION, 6, *w, result, error);

  ncm_memory_pool_return (w);
  if (error_code != GSL_SUCCESS && error_code != GSL_EROUND)
    *result = GSL_POSINF;
  return error_code;
}

/**
 * nc_integral_locked_a_inf: (skip)
 * @F: a gsl_function wich is the integrand.
 * @a: lower integration limit.
 * @abstol: absolute tolerance.
 * @reltol: relative tolerance.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function uses a workspace from the pool and gsl_integration_qagiu function to perform
 * the numerical integration in the \f$ [a, \infty] \f$ interval.
 *
 * Returns: the error code returned by gsl_integration_qagiu.
 */
gint
nc_integral_locked_a_inf (gsl_function *F, gdouble a, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error)
{
  gsl_integration_workspace **w = nc_integral_get_workspace ();
  gint error_code = gsl_integration_qagiu (F, a, abstol, reltol, NC_INT_PARTITION, *w, result, error);
  ncm_memory_pool_return (w);
  if (error_code != GSL_SUCCESS && error_code != GSL_EROUND)
  {
    g_warning ("nc_integral_locked_a_inf: %s", gsl_strerror (error_code));
    *result = GSL_POSINF;
  }
  return error_code;
}

/**
 * nc_integral_cached_0_x: (skip)
 * @cache: a pointer to #NcFunctionCache.
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
nc_integral_cached_0_x (NcFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error)
{
  gdouble x_found = 0.0;
//  gdouble p_result = 0.0;
  gsl_vector *p_result;
  gint error_code = GSL_SUCCESS;

//printf ("[%p]SEARCH! -> %g\n", g_thread_self (), x);
  if (nc_function_cache_get_near (cache, x, &x_found, &p_result, NC_FUNCTION_CACHE_SEARCH_BOTH))
  {
//printf ("[%p]Found out %g %g [%p]\n", g_thread_self (), x_found, gsl_vector_get (p_result, 0), p_result);
    if (x == x_found)
      *result = gsl_vector_get (p_result, 0);
    else
    {
      error_code = nc_integral_locked_a_b (F, x_found, x, 0.0, NC_INT_ERROR, result, error);
      *result += gsl_vector_get (p_result, 0);
      nc_function_cache_insert (cache, x, *result);
    }
  }
  else
  {
    error_code = nc_integral_locked_a_b (F, 0.0, x, 0.0, NC_INT_ERROR, result, error);
    nc_function_cache_insert (cache, x, *result);
  }
  return error_code;
}

/**
 * nc_integral_cached_x_inf: (skip)
 * @cache: a pointer to #NcFunctionCache.
 * @F: a gsl_function wich is the integrand.
 * @x: lower integration limit.
 * @result: a pointer to a gdouble in which the function stores the result.
 * @error: a pointer to a gdouble in which the function stores the estimated error.
 *
 * This function searchs for the nearest x_near value previously chosed as the lower integration limit
 * and perform the integration at \f$ [x, x_{near}] \f$ interval. This result is summed to that
 * obtained at \f$ [x_{near}, \infty] \f$ and then it is saved in the cache.
 *
 * Returns: the error code returned by gsl_integration_qagiu.
 */
gint
nc_integral_cached_x_inf (NcFunctionCache *cache, gsl_function *F, gdouble x, gdouble *result, gdouble *error)
{
  gdouble x_found = 0.0;
//  gdouble p_result = 0.0;
  gsl_vector *p_result;
  gint error_code = GSL_SUCCESS;

  if (nc_function_cache_get_near (cache, x, &x_found, &p_result, NC_FUNCTION_CACHE_SEARCH_BOTH))
  {
    if (x == x_found)
      *result = gsl_vector_get(p_result, 0);
    else
    {
      error_code = nc_integral_locked_a_b (F, x, x_found, cache->abstol, cache->reltol, result, error);
      *result += gsl_vector_get(p_result, 0);
      nc_function_cache_insert (cache, x, *result);
    }
  }
  else
  {
    error_code = nc_integral_locked_a_inf (F, x, cache->abstol, cache->reltol, result, error);
    nc_function_cache_insert (cache, x, *result);
  }
  return error_code;
}

typedef struct _iCLIntegrand2dim
{
	NcIntegrand2dim *integ;
	gdouble xi;
	gdouble xf;
	gdouble yi;
	gdouble yf;
} iCLIntegrand2dim;

#ifdef HAVE_LIBCUBA
static gint
_integrand_2dim (const gint *ndim, const gdouble x[], const gint *ncomp, gdouble f[], gpointer userdata)
{
	iCLIntegrand2dim *iinteg = (iCLIntegrand2dim *) userdata;
	f[0] = iinteg->integ->f ((iinteg->xf - iinteg->xi) * x[0] + iinteg->xi, (iinteg->yf - iinteg->yi) * x[1] + iinteg->yi, iinteg->integ->userdata);
	return 0;
}
#endif /* HAVE_LIBCUBA */

/**
 * ncm_integrate_2dim:
 * @integ: a pointer to #NcIntegrand2dim.
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
ncm_integrate_2dim (NcIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error)
{
  gboolean ret = FALSE;
#ifdef HAVE_LIBCUBA
	const gint mineval = 1;
	const gint maxeval = 10000;
	const gint key = 13; /* 13 points rule */
  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf};
	gint nregions, neval, fail;
	gdouble prob;

	Cuhre (2, 1, &_integrand_2dim, &iinteg, epsrel, epsabs, 0, mineval, maxeval,
	      key, &nregions, &neval, &fail,
	      result, error, &prob);

	*result *= (xf - xi) * (yf - yi);
	*error *= (xf - xi) * (yf - yi);

	ret = (fail == 0);
#else
	g_assert_not_reached ();
#endif
	return ret;
}

/**
 * nc_integral_fixed_new: (skip)
 * @n_nodes: number of nodes in the full interval.
 * @rule_n: order of the Gauss-Legendre integration rule to be applied in each interval.
 * @xl: the interval lower limit.
 * @xu: the interval upper limit.
 *
 * This function prepares the #NcIntegralFixed with a grid
 * with n_nodes - 1 intervals beteween xl and xu. In each interval it uses
 * a fixed order (rule_n) Gauss-Legendre integration rule to determine the
 * interval inner points. This results in a grid with (n_nodes - 1) * rule_n points.
 *
 * Returns: a pointer to the newly created #NcIntegralFixed structure.
*/
NcIntegralFixed *
nc_integral_fixed_new (gulong n_nodes, gulong rule_n, gdouble xl, gdouble xu)
{
  NcIntegralFixed *intf = g_slice_new (NcIntegralFixed);
  intf->n_nodes = n_nodes;
  intf->rule_n = rule_n;
  intf->int_nodes = g_slice_alloc (sizeof(gdouble) * n_nodes * rule_n);
  intf->xl = xl;
  intf->xu = xu;

#ifdef HAVE_GSL_GLF
  intf->glt = gsl_integration_glfixed_table_alloc (rule_n);
#else
  g_error ("NcIntegralFixed: Need gsl version > 1.4\n");
#endif /* HAVE_GSL_GLF */

  return intf;
}

/**
 * nc_integral_fixed_free:
 * @intf: a pointer to #NcIntegralFixed.
 *
 * This function frees the memory associated to #NcIntegralFixed.
 *
*/
void
nc_integral_fixed_free (NcIntegralFixed *intf)
{
  g_slice_free1 (sizeof(gdouble) * intf->n_nodes * intf->rule_n, intf->int_nodes);
  gsl_integration_glfixed_table_free (intf->glt);
  g_slice_free (NcIntegralFixed, intf);
}

/**
 * nc_integral_fixed_calc_nodes: (skip)
 * @intf: a pointer to #NcIntegralFixed.
 * @F: a pointer to a gsl_function.
 *
 * This function FIXME
 *
*/
void
nc_integral_fixed_calc_nodes (NcIntegralFixed *intf, gsl_function *F)
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
 * nc_integral_fixed_nodes_eval:
 * @intf: a pointer to #NcIntegralFixed.
 *
 * This function
 *
 * Returns: FIXME
*/
gdouble
nc_integral_fixed_nodes_eval (NcIntegralFixed *intf)
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
 * nc_integral_fixed_integ_mult: (skip)
 * @intf: a pointer to #NcIntegralFixed.
 * @F: a pointer to gsl_function.
 *
 * This function
 *
 * Returns: FIXME
*/
gdouble
nc_integral_fixed_integ_mult (NcIntegralFixed *intf, gsl_function *F)
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
 * nc_integral_fixed_integ_posdef_mult: (skip)
 * @intf: a pointer to #NcIntegralFixed.
 * @F: a pointer to gsl_function.
 * @max: FIXME
 * @reltol: FIXME
 *
 * This function
 *
 * Returns: FIXME
*/
gdouble
nc_integral_fixed_integ_posdef_mult (NcIntegralFixed *intf, gsl_function *F, gdouble max, gdouble reltol)
{
  const gulong r2 = intf->rule_n / 2;
  const gboolean odd_rule = intf->rule_n & 1;
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  const glong mnode = max / delta_x;
  gdouble res = 0.0;
  glong i, j, k = 0;

  g_assert (mnode < intf->n_nodes);

  if (odd_rule)
  {
    for (i = mnode; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part = 0.0;
      k = i * intf->rule_n;

      for (j = 1; j < r2 + 1; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      part += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];
      for (j = 1; j < r2 + 1; j++)
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

      for (j = 1; j < r2 + 1; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      part += GSL_FN_EVAL (F, x1px0_2) * intf->int_nodes[k++];
      for (j = 1; j < r2 + 1; j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      res += part;
      if (fabs(part / res) < reltol)
        break;
    }
  }
  else
  {
    for (i = mnode; i < intf->n_nodes - 1; i++)
    {
      const gdouble x0 = intf->xl + delta_x * i;
      const gdouble x1 = x0 + delta_x;
      const gdouble x1px0_2 = (x1 + x0) / 2.0;
      const gdouble x1mx0_2 = (x1 - x0) / 2.0;
      gdouble part = 0.0;
      k = i * intf->rule_n;

      for (j = 0; j < r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      for (j = 0; j < r2; j++)
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

      for (j = 0; j < r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 - x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      for (j = 0; j < r2; j++)
        part += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->glt->x[j]) * intf->int_nodes[k++];
      res += part;

      if (fabs(part / res) < reltol)
        break;
    }
  }

  return res * delta_x * 0.5;
}

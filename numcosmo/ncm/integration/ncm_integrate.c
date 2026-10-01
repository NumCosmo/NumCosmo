/***************************************************************************
 *            ncm_integrate.c
 *
 *  Wed Aug 13 20:35:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br>, <pennalima@gmail.com>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti & Mariana Penna Lima 2012 <vitenti@uel.br>, <pennalima@gmail.com>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * - ncm_integral_locked_a_b() and ncm_integral_locked_a_inf(): adaptive quadrature of
 *   GSL with a workspace taken from a thread-safe pool, see ncm_integral_get_workspace();
 * - ncm_integral_cached_0_x() and ncm_integral_cached_x_inf(): integrals from $0$ or to
 *   $\infty$ that reuse the value cached in a #NcmFunctionCache at the nearest limit;
 * - ncm_integrate_2dim(), ncm_integrate_2dim_divonne() and ncm_integrate_3dim_divonne():
 *   integrals over rectangles and boxes with the Cuhre and Divonne algorithms of the Cuba
 *   library, which report failure by their return value;
 * - #NcmIntegralFixed: Gauss-Legendre rules on uniform panels, with a weight evaluated
 *   once at the nodes and reused for several integrands.
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
 * Takes a GSL integration workspace of %NCM_INTEGRAL_PARTITION subintervals from a
 * thread-safe pool, allocating one when the pool is empty. It must be given back with
 * ncm_memory_pool_return().
 *
 * Returns: a pointer to the workspace pointer.
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
 * @F: the integrand
 * @a: the lower limit
 * @b: the upper limit
 * @abstol: the absolute tolerance
 * @reltol: the relative tolerance
 * @result: (out): the integral
 * @error: (out): the error estimate
 *
 * Integrates @F over [@a, @b] with gsl_integration_qag(), the 61-point rule and a pooled
 * workspace; a GSL failure other than %GSL_EROUND aborts.
 *
 * Returns: the GSL status.
 */
gint
ncm_integral_locked_a_b (gsl_function *F, gdouble a, gdouble b, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error)
{
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gint error_code               = gsl_integration_qag (F, a, b, abstol, reltol, NCM_INTEGRAL_PARTITION, 6, *w, result, error);

  ncm_memory_pool_return (w);

  if ((error_code != GSL_SUCCESS) && (error_code != GSL_EROUND))
    g_error ("ncm_integral_locked_a_b: %s", gsl_strerror (error_code));

  return error_code;
}

/**
 * ncm_integral_locked_a_inf: (skip)
 * @F: the integrand
 * @a: the lower limit
 * @abstol: the absolute tolerance
 * @reltol: the relative tolerance
 * @result: (out): the integral
 * @error: (out): the error estimate
 *
 * Integrates @F over $[a, \infty)$ with gsl_integration_qagiu() and a pooled workspace;
 * a GSL failure other than %GSL_EROUND aborts.
 *
 * Returns: the GSL status.
 */
gint
ncm_integral_locked_a_inf (gsl_function *F, gdouble a, gdouble abstol, gdouble reltol, gdouble *result, gdouble *error)
{
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gint error_code               = gsl_integration_qagiu (F, a, abstol, reltol, NCM_INTEGRAL_PARTITION, *w, result, error);

  ncm_memory_pool_return (w);

  if ((error_code != GSL_SUCCESS) && (error_code != GSL_EROUND))
    g_error ("ncm_integral_locked_a_inf: %s", gsl_strerror (error_code));

  return error_code;
}

/**
 * ncm_integral_cached_0_x: (skip)
 * @cache: the cache of previous integrals
 * @F: the integrand
 * @x: the upper limit
 * @result: (out): the integral
 * @error: (out): the error estimate of the last segment integrated
 *
 * Computes $\int_0^x F$ as the cached integral up to the nearest cached limit $x_n$ plus
 * $\int_{x_n}^x F$, with the tolerances of @cache, and caches the result at @x.
 *
 * Returns: the GSL status.
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
      error_code = ncm_integral_locked_a_b (F, x_found, x, ncm_function_cache_get_abstol (cache), ncm_function_cache_get_reltol (cache), result, error);
      *result   += ncm_vector_get (p_result, 0);
      ncm_function_cache_insert (cache, x, *result);
    }

    ncm_vector_clear (&p_result);
  }
  else
  {
    error_code = ncm_integral_locked_a_b (F, 0.0, x, ncm_function_cache_get_abstol (cache), ncm_function_cache_get_reltol (cache), result, error);
    ncm_function_cache_insert (cache, x, *result);
  }

  return error_code;
}

/**
 * ncm_integral_cached_x_inf: (skip)
 * @cache: the cache of previous integrals
 * @F: the integrand
 * @x: the lower limit
 * @result: (out): the integral
 * @error: (out): the error estimate of the last segment integrated
 *
 * Computes $\int_x^\infty F$ as the cached integral from the nearest cached limit $x_n$
 * plus $\int_x^{x_n} F$, with the tolerances of @cache, and caches the result at @x.
 *
 * Returns: the GSL status.
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

/**
 * ncm_integrate_2dim:
 * @integ: the integrand
 * @xi: the lower limit in $x$
 * @yi: the lower limit in $y$
 * @xf: the upper limit in $x$
 * @yf: the upper limit in $y$
 * @epsrel: the relative tolerance
 * @epsabs: the absolute tolerance
 * @result: (out): the integral
 * @error: (out): the error estimate
 *
 * Integrates over [@xi, @xf] by [@yi, @yf] with the Cuhre algorithm of Cuba, at most
 * $10^7$ evaluations, warning when they are exhausted.
 *
 * Returns: whether Cuba reached the tolerance.
 */
gboolean
ncm_integrate_2dim (NcmIntegrand2dim *integ, gdouble xi, gdouble yi, gdouble xf, gdouble yf, gdouble epsrel, gdouble epsabs, gdouble *result, gdouble *error)
{
  gboolean ret            = FALSE;
  const gint mineval      = 1;
  const gint maxeval      = 10000000;
  const gint key          = 11; /* 13 points rule */
  iCLIntegrand2dim iinteg = {integ, xi, xf, yi, yf, 0};
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
 * ncm_integrate_2dim_divonne:
 * @integ: the integrand
 * @xi: the lower limit in $x$
 * @yi: the lower limit in $y$
 * @xf: the upper limit in $x$
 * @yf: the upper limit in $y$
 * @epsrel: the relative tolerance
 * @epsabs: the absolute tolerance
 * @ngiven: the number of points in @xgiven
 * @ldxgiven: the offset between consecutive points in @xgiven
 * @xgiven: the points where the integrand may peak
 * @result: (out): the integral
 * @error: (out): the error estimate
 *
 * Integrates over [@xi, @xf] by [@yi, @yf] with the Divonne algorithm of Cuba, at most
 * $10^7$ evaluations; aborts when Cuba is older than 4.0.
 *
 * Returns: whether Cuba reached the tolerance.
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
  iCLIntegrand2dim iinteg    = {integ, xi, xf, yi, yf, ldxgiven};
  gdouble *xgiven_u;
  gint nregions, neval, fail, i;
  gdouble prob;

  /* The points rescaled to the unit square, leaving the caller's array unchanged */
  xgiven_u = g_memdup2 (xgiven, sizeof (gdouble) * ngiven * ldxgiven);

  for (i = 0; i < ngiven; i++)
  {
    xgiven_u[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven_u[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
  }

#ifdef HAVE_LIBCUBA_4_0
  Divonne (2, 1, &_integrand_2dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border,
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven_u, nextra, peakfinder, NULL, NULL, &nregions, &neval, &fail,
           result, error, &prob);

  g_free (xgiven_u);

  if (neval >= maxeval)
    g_warning ("ncm_integrate_2dim_divonne: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi);
  *error  *= (xf - xi) * (yf - yi);

  ret = (fail == 0);

  return ret;

#else
  g_free (xgiven_u);
  g_error ("ncm_integrate_2dim_divonne: Needs libcuba > 4.0.");

  return FALSE;

#endif /*HAVE_LIBCUBA_4_0 */
}

/**
 * ncm_integrate_3dim_divonne:
 * @integ: the integrand
 * @xi: the lower limit in $x$
 * @yi: the lower limit in $y$
 * @zi: the lower limit in $z$
 * @xf: the upper limit in $x$
 * @yf: the upper limit in $y$
 * @zf: the upper limit in $z$
 * @epsrel: the relative tolerance
 * @epsabs: the absolute tolerance
 * @ngiven: the number of points in @xgiven
 * @ldxgiven: the offset between consecutive points in @xgiven
 * @xgiven: the points where the integrand may peak
 * @result: (out): the integral
 * @error: (out): the error estimate
 *
 * Integrates over the box [@xi, @xf] by [@yi, @yf] by [@zi, @zf] with the Divonne
 * algorithm of Cuba, with no evaluation limit; aborts when Cuba is older than 4.0.
 *
 * Returns: whether Cuba reached the tolerance.
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
  gdouble *xgiven_u;
  gint nregions, neval, fail, i;
  gdouble prob;

  /* The points rescaled to the unit box, leaving the caller's array unchanged */
  xgiven_u = g_memdup2 (xgiven, sizeof (gdouble) * ngiven * ldxgiven);

  for (i = 0; i < ngiven; i++)
  {
    xgiven_u[i * ldxgiven + 0] = (xgiven[i * ldxgiven + 0] - xi) / (xf - xi);
    xgiven_u[i * ldxgiven + 1] = (xgiven[i * ldxgiven + 1] - yi) / (yf - yi);
    xgiven_u[i * ldxgiven + 2] = (xgiven[i * ldxgiven + 2] - zi) / (zf - zi);
  }

#ifdef HAVE_LIBCUBA_4_0
  Divonne (3, 1, &_integrand_3dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border,
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven_u, nextra, peakfinder, NULL, NULL, &nregions, &neval, &fail,
           result, error, &prob);
  g_free (xgiven_u);
#else
  g_free (xgiven_u);
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
 * ncm_integral_fixed_new: (skip)
 * @n_nodes: the number of panel edges, at least 2
 * @rule_n: the number of Gauss-Legendre points per panel
 * @xl: the lower limit
 * @xu: the upper limit
 *
 * Creates a #NcmIntegralFixed with @n_nodes - 1 equal panels on [@xl, @xu] and an
 * @rule_n-point Gauss-Legendre rule in each, $(n_\mathrm{nodes} - 1)\,r$ nodes in all.
 *
 * Returns: a new #NcmIntegralFixed, to be freed with ncm_integral_fixed_free().
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
 * @intf: a #NcmIntegralFixed
 *
 * Frees @intf.
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
 * @intf: a #NcmIntegralFixed
 * @F: the weight
 *
 * Stores the weight @F times the Gauss-Legendre weight at every node, for the integrals
 * of @intf.
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
 * @intf: a #NcmIntegralFixed
 *
 * Returns: the integral of the weight $F$ stored by ncm_integral_fixed_calc_nodes().
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
 * @intf: a #NcmIntegralFixed
 * @F: the function multiplying the weight
 *
 * Returns: $\int F_w(x)\,F(x)\,\mathrm{d}x$, where $F_w$ is the weight stored by
 * ncm_integral_fixed_calc_nodes().
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
 * ncm_integral_fixed_integ_vec_mult:
 * @intf: a #NcmIntegralFixed with the weight stored by ncm_integral_fixed_calc_nodes()
 * @f_at_nodes: the values of $G$ at the nodes of ncm_integral_fixed_get_nodes()
 *
 * Same as ncm_integral_fixed_integ_mult() with $G$ given at the nodes.
 *
 * Returns: $\int F_w(x)\,G(x)\,\mathrm{d}x$.
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
 * @intf: a #NcmIntegralFixed
 * @nodes: the output vector, of length $(n_\mathrm{nodes} - 1)\,r$
 *
 * Computes the nodes into @nodes, in the order of the weights stored by
 * ncm_integral_fixed_calc_nodes().
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
 * @F: the weight, stored at the nodes
 * @G: the function multiplying the weight
 * @xl: the lower limit
 * @xu: the upper limit
 * @reltol: the relative tolerance
 * @exact_F_integ: the exact $\int F$ over [@xl, @xu], or %GSL_NAN
 * @max_total_nodes: the maximum total number of nodes $(n_\mathrm{nodes} - 1)\,r$
 * @n_nodes_out: (out): the number of panel edges chosen
 * @rule_n_out: (out): the number of Gauss-Legendre points per panel chosen
 * @relerr_out: (out) (nullable): the relative error in $\int F\,G$ reached
 *
 * Chooses the #NcmIntegralFixed with the fewest total nodes, among rules of 3, 5 and 7
 * points, whose $\int F\,G$ agrees with a reference to @reltol and whose $\int F$ agrees
 * to @reltol with @exact_F_integ, or with the reference's own $\int F$ when
 * @exact_F_integ is not finite. The second test catches features of @F narrower than a
 * panel, which both integrals of $F\,G$ could miss. The reference uses 7-point panels,
 * $10\,\mathrm{reltol}^{-0.3}$ of them, from 96 to 2048. For each rule the panel count
 * is bracketed by growing it geometrically and then found by bisection, assuming the
 * error decreases with it.
 *
 * When no configuration within @max_total_nodes meets @reltol, the one with the smallest
 * error in $\int F\,G$ is used and, if @relerr_out is %NULL, a warning is emitted.
 * @relerr_out receives the error of the configuration returned, relative to the
 * reference.
 *
 * Returns: (transfer full): a new #NcmIntegralFixed with @F stored, see
 * ncm_integral_fixed_calc_nodes().
 */
NcmIntegralFixed *
ncm_integral_fixed_calibrate (gsl_function *F, gsl_function *G, gdouble xl, gdouble xu, gdouble reltol, gdouble exact_F_integ, gulong max_total_nodes, guint *n_nodes_out, guint *rule_n_out, gdouble *relerr_out)
{
  /* Reference panels growing as reltol^(-0.3), faster than the ~4th-order convergence
   * of a piecewise-smooth weight, so the reference stays finer than any configuration
   * that passes */
  const gulong ref_n_nodes      = (gulong) CLAMP ((glong) (10.0 * pow (reltol, -0.3) + 0.5), 96, 2048);
  const gulong ref_rule_n       = 7;
  const guint rule_candidates[] = { 3, 5, 7 };
  gdouble I_ref, denom_I, guard_mass;
  gulong best_n_nodes = 0, best_rule_n = 0, best_total = 0;
  gdouble best_err = GSL_POSINF;
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
    gdouble hi_err         = GSL_POSINF;
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
        hi     = n_try;
        hi_err = err_I;
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
        {
          hi     = mid;
          hi_err = err_I;
        }
        else
        {
          lo = mid;
        }

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
        best_err     = hi_err;
        converged    = TRUE;
      }
    }
  }

  if (!converged)
  {
    best_n_nodes = fb_n_nodes;
    best_rule_n  = fb_rule_n;

    /* Silent when the caller asked for the achieved error: it is then reporting
     * on our behalf, and this function is called in a per-element loop where one
     * warning each is thousands of duplicates. */
    if (relerr_out == NULL)
      g_warning ("ncm_integral_fixed_calibrate: did not reach reltol %g within %lu total nodes "
                 "(best rel error %g at n_nodes=%lu, rule_n=%lu).",
                 reltol, max_total_nodes, fb_err, best_n_nodes, best_rule_n);
  }

  if (relerr_out != NULL)
    *relerr_out = converged ? best_err : fb_err;

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


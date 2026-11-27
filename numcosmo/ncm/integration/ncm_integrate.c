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

  Cuhre (2, 1, &_integrand_2dim, &iinteg, 1, epsrel, epsabs, 0, mineval, maxeval, key, NULL, NULL, &nregions, &neval, &fail, result, error, &prob);

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
 * $10^7$ evaluations.
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
 * algorithm of Cuba, with no evaluation limit.
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

  Divonne (3, 1, &_integrand_3dim, &iinteg, nvec, epsrel, epsabs, 0, seed, mineval, maxeval, key1, key2, key3, maxpass, border,
           maxchisq, mindeviation, ngiven, ldxgiven, xgiven_u, nextra, peakfinder, NULL, NULL, &nregions, &neval, &fail,
           result, error, &prob);
  g_free (xgiven_u);

  if (neval >= maxeval)
    g_warning ("ncm_integrate_3dim_divonne: number of evaluations %d >= maximum number of evaluations %d.\n", neval, maxeval);

  *result *= (xf - xi) * (yf - yi) * (zf - zi);
  *error  *= (xf - xi) * (yf - yi) * (zf - zi);

  ret = (fail == 0);

  return ret;
}

/* The rules of ncm_integral_gauss_legendre_rule(): at index n, the n nodes then the n weights */
static GPtrArray *_ncm_integral_gl_rules = NULL;

G_LOCK_DEFINE_STATIC (_ncm_integral_gl_rules);

static gint
_ncm_integral_gl_cmp (gconstpointer a, gconstpointer b)
{
  const gdouble xa = ((const gdouble *) a)[0];
  const gdouble xb = ((const gdouble *) b)[0];

  return (xa > xb) - (xa < xb);
}

/**
 * ncm_integral_gauss_legendre_rule: (skip)
 * @n: the number of points, at least 1
 * @x: (out) (transfer none): the nodes on $[-1, 1]$, in ascending order
 * @w: (out) (transfer none): the weights
 *
 * Gets the @n-point Gauss-Legendre rule on $[-1, 1]$, exact for polynomials of degree
 * $2n - 1$. The rule is the table of GSL when GSL has one for @n (gsl_integration_glfixed),
 * and otherwise the one of gsl_integration_fixed_legendre, which computes it from the
 * eigenvalues of the Jacobi matrix: the rules gsl_integration_glfixed_table_alloc()
 * computes for the other sizes carry errors near $10^{-11}$ from about 48 points. GSL has
 * tables for n = 2 to 20, 32, 64, 96, 100, 128, 256, 512 and 1024. Measured on the integrals
 * of $x^k$, $k \leq 2n - 1$, for every n up to 257: at most $1.1 \times 10^{-16}$ with a
 * table and $6.1 \times 10^{-15}$ without. Each rule is computed once and kept until the
 * process ends; the call is thread safe.
 */
void
ncm_integral_gauss_legendre_rule (const guint n, const gdouble **x, const gdouble **w)
{
  gdouble *rule;

  if (n == 0)
    g_error ("ncm_integral_gauss_legendre_rule: a rule needs at least one point.");

  G_LOCK (_ncm_integral_gl_rules);

  if (_ncm_integral_gl_rules == NULL)
    _ncm_integral_gl_rules = g_ptr_array_new ();

  if (_ncm_integral_gl_rules->len <= n)
    g_ptr_array_set_size (_ncm_integral_gl_rules, n + 1);

  rule = g_ptr_array_index (_ncm_integral_gl_rules, n);

  if (rule == NULL)
  {
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc (n);
    gdouble *pairs                   = g_new (gdouble, 2 * n);
    guint i;

    if (t->precomputed)
    {
      for (i = 0; i < n; i++)
        gsl_integration_glfixed_point (-1.0, 1.0, i, &pairs[2 * i], &pairs[2 * i + 1], t);
    }
    else
    {
      gsl_integration_fixed_workspace *ws = gsl_integration_fixed_alloc (gsl_integration_fixed_legendre, n, -1.0, 1.0, 0.0, 0.0);
      const gdouble *fx                   = gsl_integration_fixed_nodes (ws);
      const gdouble *fw                   = gsl_integration_fixed_weights (ws);

      for (i = 0; i < n; i++)
      {
        pairs[2 * i]     = fx[i];
        pairs[2 * i + 1] = fw[i];
      }

      gsl_integration_fixed_free (ws);
    }

    gsl_integration_glfixed_table_free (t);

    qsort (pairs, n, 2 * sizeof (gdouble), &_ncm_integral_gl_cmp);

    rule = g_new (gdouble, 2 * n);

    for (i = 0; i < n; i++)
    {
      rule[i]     = pairs[2 * i];
      rule[n + i] = pairs[2 * i + 1];
    }

    g_free (pairs);
    g_ptr_array_index (_ncm_integral_gl_rules, n) = rule;
  }

  G_UNLOCK (_ncm_integral_gl_rules);

  *x = rule;
  *w = rule + n;
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

  ncm_integral_gauss_legendre_rule (rule_n, &intf->gl_x, &intf->gl_w);

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
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gulong i, j, k = 0;

  for (i = 0; i < intf->n_nodes - 1; i++)
  {
    const gdouble x0      = intf->xl + delta_x * i;
    const gdouble x1      = x0 + delta_x;
    const gdouble x1px0_2 = (x1 + x0) / 2.0;
    const gdouble x1mx0_2 = (x1 - x0) / 2.0;

    for (j = 0; j < intf->rule_n; j++)
      intf->int_nodes[k++] = GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->gl_x[j]) * intf->gl_w[j];
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
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  gdouble res           = 0.0;
  gulong i, j, k = 0;

  for (i = 0; i < intf->n_nodes - 1; i++)
  {
    const gdouble x0      = intf->xl + delta_x * i;
    const gdouble x1      = x0 + delta_x;
    const gdouble x1px0_2 = (x1 + x0) / 2.0;
    const gdouble x1mx0_2 = (x1 - x0) / 2.0;

    for (j = 0; j < intf->rule_n; j++)
      res += GSL_FN_EVAL (F, x1px0_2 + x1mx0_2 * intf->gl_x[j]) * intf->int_nodes[k++];
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
 * ncm_integral_fixed_calc_nodes(): panel after panel, ascending in each.
 */
void
ncm_integral_fixed_get_nodes (NcmIntegralFixed *intf, NcmVector *nodes)
{
  const gdouble delta_x = (intf->xu - intf->xl) / (intf->n_nodes - 1.0);
  const gulong total_n  = (intf->n_nodes - 1) * intf->rule_n;
  gulong i, j, k = 0;

  g_assert_nonnull (nodes);
  g_assert_cmpuint (ncm_vector_len (nodes), ==, total_n);

  for (i = 0; i < intf->n_nodes - 1; i++)
  {
    const gdouble x0      = intf->xl + delta_x * i;
    const gdouble x1      = x0 + delta_x;
    const gdouble x1px0_2 = (x1 + x0) / 2.0;
    const gdouble x1mx0_2 = (x1 - x0) / 2.0;

    for (j = 0; j < intf->rule_n; j++)
      ncm_vector_set (nodes, k++, x1px0_2 + x1mx0_2 * intf->gl_x[j]);
  }
}


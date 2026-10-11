/***************************************************************************
 *            nc_xcor_kernel_integrand.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Thu October 08 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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

/*
 * NcXcorKernelIntegrand: the closure of W_l(k) of an ell block, as the k
 * integrals see it. The boxed type and its public API, and the two private
 * data layouts behind it, spline and Chebyshev, with the accessors each sets.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/powspec/ncm_powspec.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/specfunc/ncm_sbessel_ode_solver.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "ncm/stats/ncm_function_sample_set.h"
#include "nc/background/nc_distance.h"
#include "nc/xcor/nc_xcor_kernel.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/algebra/ncm_spectral.h"
#include "ncm/core/ncm_memory_pool.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor.h"
#include "nc_enum_types.h"
#include "nc/xcor/kernel/nc_xcor_kernel_private.h"

G_DEFINE_BOXED_TYPE (NcXcorKernelIntegrand, nc_xcor_kernel_integrand, nc_xcor_kernel_integrand_ref, nc_xcor_kernel_integrand_unref)

static NcmVector *
_spline_integrand_get_knots (gpointer data)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  /* Every component of the spline vec shares the same knots, so the first
   * component's knots are the knots of the whole block. */
  return ncm_spline_peek_xv (ncm_spline_vec_peek_spline (sid->spline_vec, 0));
}

static void
_spline_integrand_eval (gpointer data, gdouble k, gdouble *W)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;
  guint i;

  ncm_spline_vec_eval (sid->spline_vec, k, sid->eval_result);

  for (i = 0; i < sid->len; i++)
  {
    W[i] = ncm_vector_get (sid->eval_result, i);
  }
}

/*
 * Evaluates the components [offset, offset + len) at k, writing each at its own
 * index of W. All components share the knots, so the interval index found on
 * the first spline is used for every one of them.
 */
static void
_spline_integrand_eval_comps (gpointer data, gdouble k, guint offset, guint len, gdouble *W)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;
  NcmSpline *spline_0      = ncm_spline_vec_peek_spline (sid->spline_vec, 0);
  const gsize idx          = ncm_spline_get_index (spline_0, k);
  guint i;

  for (i = 0; i < len; i++)
  {
    NcmSpline *spline_i = ncm_spline_vec_peek_spline (sid->spline_vec, offset + i);

    W[offset + i] = ncm_spline_eval_idx (spline_i, k, idx);
  }
}

static void
_spline_integrand_get_range (gpointer data, gdouble *kmin, gdouble *kmax)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  *kmin = sid->k_min;
  *kmax = sid->k_max;
}

/*
 * Returns the k range on which component i is nonzero. Under the Limber
 * approximation this is the band [nu / chi_max, nu / chi_min] of multipole i
 * intersected with the block's domain, and W_i steps to zero at its edges; the
 * outer integral uses the range as its limits so that the steps fall on them.
 * Otherwise it is the block's domain.
 */
static void
_spline_integrand_get_range_comp (gpointer data, guint i, gdouble *kmin, gdouble *kmax)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  *kmin = sid->k_min_comp[i];
  *kmax = sid->k_max_comp[i];
}

static void
_spline_integrand_data_free (gpointer data)
{
  SplineIntegrandData *sid = (SplineIntegrandData *) data;

  nc_hicosmo_clear (&sid->cosmo);
  ncm_spline_vec_clear (&sid->spline_vec);
  ncm_vector_clear (&sid->eval_result);

  g_free (sid->k_min_comp);
  g_free (sid->k_max_comp);
  g_free (data);
}

/* Chebyshev closure of one ell block, built by both the Limber and the non-Limber paths. */

static gpointer
_nc_xcor_spectral_alloc (gpointer userdata)
{
  return ncm_spectral_new ();
}

static void
_nc_xcor_spectral_free (gpointer p)
{
  ncm_spectral_free (NCM_SPECTRAL (p));
}

/*
 * Returns a NcmSpectral from a process-wide pool, for the exclusive use of the
 * caller until it is given back with ncm_memory_pool_return(). A NcmSpectral
 * holds buffers for the values at the nodes, coefficient buffers and FFTW plans, so it cannot be
 * used by two threads at once, and closures are built and restricted
 * concurrently. The pool keeps the plans, whose creation costs more than the
 * expansion that uses them.
 */
NcmSpectral **
_nc_xcor_kernel_spectral_get (void)
{
  G_LOCK_DEFINE_STATIC (create_lock);

  static NcmMemoryPool *mp = NULL;

  G_LOCK (create_lock);

  if (mp == NULL)
    mp = ncm_memory_pool_new (_nc_xcor_spectral_alloc, NULL, _nc_xcor_spectral_free);

  G_UNLOCK (create_lock);

  return (NcmSpectral **) ncm_memory_pool_get (mp);
}

/*
 * Evaluates the components [offset, offset + len) at k. Outside the panels' span
 * the closure is zero: the domain was expanded until the block was negligible there.
 *
 * With per-multipole scales each component lives in its own x = k / scale[i] and is
 * nonzero on its own support [x_min_comp, x_max_comp], whose ends are panel edges.
 * Whether k is inside is decided in k, against that support times the scale, and x is
 * then clamped to the support: k / scale can round past an end that k reaches exactly.
 * The upper end is read from inside the support, on the panel that ends there, where
 * the panel that starts there may already be past the component's support.
 */
static void
_cheb_integrand_eval_comps (gpointer data, gdouble k, guint offset, guint len, gdouble *W)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;
  guint i;

  if (!cid->scaled)
  {
    gdouble x_min, x_max;

    ncm_spectral_panels_get_range (cid->sp, &x_min, &x_max);

    for (i = 0; i < len; i++)
    {
      if ((k < x_min) || (k > x_max))
        W[offset + i] = 0.0;
      else
        W[offset + i] = ncm_spectral_panels_eval_row (cid->sp, offset + i, k);
    }
  }
  else
  {
    for (i = 0; i < len; i++)
    {
      const guint r = offset + i;

      if ((k < cid->k_min_comp[r]) || (k > cid->k_max_comp[r]))
      {
        W[r] = 0.0;
      }
      else
      {
        const gdouble x = GSL_MIN (GSL_MAX (k / cid->scale[r], cid->x_min_comp[r]), cid->x_max_comp[r]);
        guint p         = ncm_spectral_panels_find_panel (cid->sp, x);

        if ((x == cid->x_max_comp[r]) && (p > 0))
        {
          gdouble a, b;
          guint N;

          ncm_spectral_panels_peek_panel (cid->sp, p, &a, &b, &N);

          if (a == x)
            p--;
        }

        W[r] = ncm_spectral_panels_eval_row_on (cid->sp, p, r, x);
      }
    }
  }
}

static void
_cheb_integrand_eval (gpointer data, gdouble k, gdouble *W)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  _cheb_integrand_eval_comps (data, k, 0, cid->len, W);
}

static gdouble
_cheb_integrand_get_scale (gpointer data, guint i)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  return cid->scale[i];
}

static void
_cheb_integrand_get_range (gpointer data, gdouble *kmin, gdouble *kmax)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  *kmin = cid->k_min;
  *kmax = cid->k_max;
}

static void
_cheb_integrand_get_range_comp (gpointer data, guint i, gdouble *kmin, gdouble *kmax)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  *kmin = cid->k_min_comp[i];
  *kmax = cid->k_max_comp[i];
}

static void
_cheb_integrand_data_free (gpointer data)
{
  ChebIntegrandData *cid = (ChebIntegrandData *) data;

  nc_hicosmo_clear (&cid->cosmo);
  ncm_spectral_panels_clear (&cid->sp);

  g_free (cid->scale);
  g_free (cid->k_min_comp);
  g_free (cid->k_max_comp);
  g_free (cid->x_min_comp);
  g_free (cid->x_max_comp);
  g_free (data);
}

/*
 * Creates the integrand of a Chebyshev closure from its data, taking ownership
 * of @cid, with every Chebyshev accessor set.
 */
NcXcorKernelIntegrand *
_nc_xcor_kernel_cheb_integrand_new (guint n_l, ChebIntegrandData *cid)
{
  NcXcorKernelIntegrand *integrand = nc_xcor_kernel_integrand_new (n_l,
                                                                   _cheb_integrand_eval,
                                                                   _cheb_integrand_get_range,
                                                                   cid,
                                                                   _cheb_integrand_data_free);

  nc_xcor_kernel_integrand_set_get_range_comp (integrand, _cheb_integrand_get_range_comp);
  nc_xcor_kernel_integrand_set_eval_comps (integrand, _cheb_integrand_eval_comps);
  nc_xcor_kernel_integrand_set_get_scale (integrand, _cheb_integrand_get_scale);
  nc_xcor_kernel_integrand_set_panels (integrand, cid->sp);

  return integrand;
}

/*
 * Creates the integrand of a spline closure from its data, taking ownership of
 * @sid, with every spline accessor set.
 */
NcXcorKernelIntegrand *
_nc_xcor_kernel_spline_integrand_new (guint n_l, SplineIntegrandData *sid)
{
  NcXcorKernelIntegrand *integrand = nc_xcor_kernel_integrand_new (n_l,
                                                                   _spline_integrand_eval,
                                                                   _spline_integrand_get_range,
                                                                   sid,
                                                                   _spline_integrand_data_free);

  nc_xcor_kernel_integrand_set_get_knots (integrand, _spline_integrand_get_knots);
  nc_xcor_kernel_integrand_set_get_range_comp (integrand, _spline_integrand_get_range_comp);
  nc_xcor_kernel_integrand_set_eval_comps (integrand, _spline_integrand_eval_comps);

  return integrand;
}

/**
 * nc_xcor_kernel_integrand_new:
 * @len: number of components in the integrand
 * @eval: (scope async): function to evaluate the integrand
 * @get_range: (scope async): function to get the k range
 * @data: (nullable): user data to pass to @eval and @get_range
 * @data_free: (nullable): function to free @data
 *
 * Creates a new #NcXcorKernelIntegrand with reference count of 1.
 *
 * Returns: (transfer full): a new #NcXcorKernelIntegrand
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_integrand_new (guint len, void (*eval) (gpointer, gdouble, gdouble *), void (*get_range) (gpointer, gdouble *, gdouble *), gpointer data, GDestroyNotify data_free)
{
  NcXcorKernelIntegrand *integrand = g_new (NcXcorKernelIntegrand, 1);

  integrand->refcount       = 1;
  integrand->len            = len;
  integrand->eval_func      = eval;
  integrand->get_range_func = get_range;
  integrand->data           = data;
  integrand->data_free      = data_free;
  integrand->get_knots_func = NULL;

  integrand->get_range_comp_func = NULL;
  integrand->eval_comps_func     = NULL;
  integrand->get_scale_func      = NULL;

  integrand->closure_error = NULL;
  integrand->panels        = NULL;
  integrand->reltol        = 0.0;
  integrand->peak_epsilon  = 0.0;

  return integrand;
}

/**
 * nc_xcor_kernel_integrand_set_get_knots: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_knots: (scope async): function returning @integrand's knots
 *
 * Declares @integrand as spline-backed, by installing the accessor returning
 * the knots its components are represented on. Left unset by
 * nc_xcor_kernel_integrand_new(), so integrands that are not spline-backed
 * report no knots.
 *
 */
void
nc_xcor_kernel_integrand_set_get_knots (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetKnots get_knots)
{
  integrand->get_knots_func = get_knots;
}

/**
 * nc_xcor_kernel_integrand_set_get_scale: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_scale: (scope async): function reporting the scale of one component
 *
 * Installs the accessor reporting the scale $s_i$ of each component. Left
 * unset by nc_xcor_kernel_integrand_new(), in which case every scale is 1.
 */
void
nc_xcor_kernel_integrand_set_get_scale (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetScale get_scale)
{
  integrand->get_scale_func = get_scale;
}

/**
 * nc_xcor_kernel_integrand_get_scale:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: component index
 *
 * The panels of a spectral integrand are in the variable $x = k / s_i$, with
 * $s_i$ the scale of component @i: $s_i = 1$ for a closure built in $k$, and
 * $s_i = \nu_i = \ell_i + 1/2$ for a Limber closure built in $u = k / \nu$, whose
 * panels are then common to the block. The panels of
 * nc_xcor_kernel_integrand_peek_panels() are in $x$.
 *
 * Returns: the scale of component @i
 */
gdouble
nc_xcor_kernel_integrand_get_scale (NcXcorKernelIntegrand *integrand, guint i)
{
  if (integrand->get_scale_func == NULL)
    return 1.0;

  return integrand->get_scale_func (integrand->data, i);
}

/**
 * nc_xcor_kernel_integrand_set_get_range_comp: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @get_range_comp: (scope async): function returning one component's k range
 *
 * Installs the accessor returning the k range a single component of @integrand
 * is supported on. Left unset by nc_xcor_kernel_integrand_new(), in which case
 * every component reports the whole range.
 *
 */
void
nc_xcor_kernel_integrand_set_get_range_comp (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetRangeComp get_range_comp)
{
  integrand->get_range_comp_func = get_range_comp;
}

/**
 * nc_xcor_kernel_integrand_set_eval_comps: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @eval_comps: (scope async): function evaluating a run of components
 *
 * Installs the accessor evaluating a contiguous run of @integrand's
 * components, for callers that integrate the run on its own. Left unset by
 * nc_xcor_kernel_integrand_new(), in which case a run is served by evaluating
 * every component.
 *
 */
void
nc_xcor_kernel_integrand_set_eval_comps (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandEvalComps eval_comps)
{
  integrand->eval_comps_func = eval_comps;
}

/**
 * nc_xcor_kernel_integrand_peek_knots:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the knots @integrand's components are represented on, shared by every
 * component (multipole) it carries, or %NULL when @integrand is not
 * spline-backed.
 *
 * These knots are what makes the outer $k$ integral exactly integrable, and
 * are why %NC_XCOR_METHOD_KERNEL_EXACT needs no tolerance. Each component is a
 * cubic spline in $k$, so on any interval over which both members of a pair
 * are a single cubic piece, the product $k^2 W_i(k) W_j(k)$ entering $C_\ell$
 * is a polynomial of degree $8$ and a $5$-node Gauss-Legendre rule integrates
 * it exactly.
 *
 * Two closures are built independently and so do not share knots: the
 * intervals with that property are the panels of the *common refinement* of
 * the two knot sets. Merging two sorted knot vectors is all the coupling the
 * argument needs -- building the closures jointly on one shared knot set
 * would also work but costs about twice as much to produce, and forces every
 * pair to carry every kernel's knots.
 *
 * Returns: (transfer none) (nullable): the knot vector, or %NULL.
 */

/**
 * nc_xcor_kernel_integrand_set_tolerances:
 * @integrand: a #NcXcorKernelIntegrand
 * @reltol: relative tolerance of the acceptance test
 * @peak_epsilon: absolute tolerance of the acceptance test, as a fraction of the closure's peak
 *
 * Records the two tolerances the closure was built with, see
 * #NcXcorKernel:reltol and #NcXcorKernel:peak-epsilon. nc_xcor_compute_full()
 * uses them in its error estimate of $C_\ell$ where no interpolation error was
 * recorded, see nc_xcor_kernel_integrand_set_closure_error().
 */
void
nc_xcor_kernel_integrand_set_tolerances (NcXcorKernelIntegrand *integrand, gdouble reltol, gdouble peak_epsilon)
{
  g_return_if_fail (integrand != NULL);
  g_return_if_fail (reltol >= 0.0);
  g_return_if_fail (peak_epsilon >= 0.0);

  integrand->reltol       = reltol;
  integrand->peak_epsilon = peak_epsilon;
}

/**
 * nc_xcor_kernel_integrand_get_reltol:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Returns: the relative tolerance of the acceptance test, or 0.0 when exact or
 * unknown. See nc_xcor_kernel_integrand_set_tolerances().
 */
gdouble
nc_xcor_kernel_integrand_get_reltol (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, 0.0);

  return integrand->reltol;
}

/**
 * nc_xcor_kernel_integrand_get_peak_epsilon:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Returns: the absolute tolerance of the acceptance test as a fraction of the
 * closure's peak, or 0.0 when there was none. See
 * nc_xcor_kernel_integrand_set_tolerances().
 */
gdouble
nc_xcor_kernel_integrand_get_peak_epsilon (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, 0.0);

  return integrand->peak_epsilon;
}

/**
 * nc_xcor_kernel_integrand_set_closure_error:
 * @integrand: a #NcXcorKernelIntegrand
 * @closure_error: (nullable): the interpolation error, or %NULL
 *
 * Records the interpolation error the closure achieved, one row per interval
 * on which it is a single polynomial (a knot interval or a panel) and one
 * column per component, as defined under #NcXcorKernel:track-closure-error.
 * nc_xcor_compute_full() uses it in its error estimate of $C_\ell$, and uses
 * the tolerances of nc_xcor_kernel_integrand_set_tolerances() on intervals
 * whose entry is NaN.
 */
void
nc_xcor_kernel_integrand_set_closure_error (NcXcorKernelIntegrand *integrand, NcmMatrix *closure_error)
{
  g_return_if_fail (integrand != NULL);

  ncm_matrix_clear (&integrand->closure_error);

  if (closure_error != NULL)
    integrand->closure_error = ncm_matrix_ref (closure_error);
}

/**
 * nc_xcor_kernel_integrand_peek_closure_error:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the interpolation error, or %NULL when the closure was built with
 * #NcXcorKernel:track-closure-error off. See
 * nc_xcor_kernel_integrand_set_closure_error().
 *
 * Returns: (transfer none) (nullable): the interpolation error matrix, or %NULL
 */
NcmMatrix *
nc_xcor_kernel_integrand_peek_closure_error (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, NULL);

  return integrand->closure_error;
}

/**
 * nc_xcor_kernel_integrand_set_panels:
 * @integrand: a #NcXcorKernelIntegrand
 * @sp: (nullable): the panels of the closure, or %NULL
 *
 * Sets the panels of a Chebyshev closure, one row per component, in the variable
 * $k / s_i$ of nc_xcor_kernel_integrand_get_scale().
 */
void
nc_xcor_kernel_integrand_set_panels (NcXcorKernelIntegrand *integrand, NcmSpectralPanels *sp)
{
  g_return_if_fail (integrand != NULL);

  ncm_spectral_panels_clear (&integrand->panels);

  if (sp != NULL)
    integrand->panels = ncm_spectral_panels_ref (sp);
}

/**
 * nc_xcor_kernel_integrand_peek_panels:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the panels of a Chebyshev closure, see nc_xcor_kernel_integrand_set_panels().
 *
 * Returns: (transfer none) (nullable): the panels, or %NULL when the closure has none
 */
NcmSpectralPanels *
nc_xcor_kernel_integrand_peek_panels (NcXcorKernelIntegrand *integrand)
{
  g_return_val_if_fail (integrand != NULL, NULL);

  return integrand->panels;
}

NcmVector *
nc_xcor_kernel_integrand_peek_knots (NcXcorKernelIntegrand *integrand)
{
  if (integrand->get_knots_func == NULL)
    return NULL;

  return integrand->get_knots_func (integrand->data);
}

/**
 * nc_xcor_kernel_integrand_ref:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Increases the reference count of @integrand by one atomically.
 *
 * Returns: (transfer full): @integrand
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_integrand_ref (NcXcorKernelIntegrand *integrand)
{
  g_atomic_int_inc (&integrand->refcount);

  return integrand;
}

/**
 * nc_xcor_kernel_integrand_unref:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Decreases the reference count of @integrand by one atomically.
 * When the reference count reaches zero, frees @integrand and its
 * associated data using the free function provided at creation time
 * (if any).
 */
void
nc_xcor_kernel_integrand_unref (NcXcorKernelIntegrand *integrand)
{
  if (g_atomic_int_dec_and_test (&integrand->refcount))
  {
    if (integrand->data_free != NULL)
      integrand->data_free (integrand->data);

    ncm_matrix_clear (&integrand->closure_error);
    ncm_spectral_panels_clear (&integrand->panels);

    g_free (integrand);
  }
}

/**
 * nc_xcor_kernel_integrand_clear:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * If *@integrand is not %NULL, decreases its reference count and
 * sets the pointer to %NULL.
 */
void
nc_xcor_kernel_integrand_clear (NcXcorKernelIntegrand **integrand)
{
  if (*integrand != NULL)
  {
    nc_xcor_kernel_integrand_unref (*integrand);
    *integrand = NULL;
  }
}

/**
 * nc_xcor_kernel_integrand_get_range:
 * @integrand: a #NcXcorKernelIntegrand
 * @k_min: (out): minimum k value
 * @k_max: (out): maximum k value
 *
 * Gets the valid k range for this integrand.
 */
/**
 * nc_xcor_kernel_integrand_get_range_comp:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: component index
 * @k_min: (out): minimum k value
 * @k_max: (out): maximum k value
 *
 * Gets the k range component @i is supported on, which can be a part of the
 * range nc_xcor_kernel_integrand_get_range() reports for the whole integrand:
 * a block of multipoles shares one domain, and under the Limber approximation
 * each of them vanishes outside its own band within it. Integrating a
 * component over its own range keeps that band edge on an integration limit
 * instead of leaving a step inside the interval.
 *
 * Falls back to the whole range for integrands that do not distinguish their
 * components.
 */
/**
 * nc_xcor_kernel_integrand_eval_comps: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 * @offset: index of the first component to evaluate
 * @len: number of components to evaluate
 * @W: (array) (out caller-allocates): full-length array to store results in
 *
 * Evaluates components [@offset, @offset + @len) at wavenumber @k, writing
 * them at their own indices in @W. Integrands that can only evaluate every
 * component at once do so, filling the whole of @W; either way the entries
 * the caller asked for are valid.
 */
/**
 * nc_xcor_kernel_integrand_eval: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 * @W: (array) (out caller-allocates): array of length @len to store results
 *
 * Evaluates the integrand at wavenumber @k, storing @len results in @W.
 */
/**
 * nc_xcor_kernel_integrand_get_len:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Gets the number of components in the integrand.
 *
 * Returns: the number of components
 */
/**
 * nc_xcor_kernel_integrand_eval_array:
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 *
 * Evaluates the integrand at wavenumber @k and returns the results
 * in a newly allocated #GArray. This is a convenience wrapper around
 * nc_xcor_kernel_integrand_eval() that handles array allocation.
 *
 * Returns: (transfer full) (element-type gdouble): a #GArray containing @len #gdouble values
 */


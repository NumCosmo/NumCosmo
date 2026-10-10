/***************************************************************************
 *            nc_xcor_kernel_closure.c
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
 * Construction of the closures of W_l(k): the component states shared by the
 * Limber and non-Limber paths, the tolerance checks, the seeds of the k domain,
 * the panelling and acceptance of the Chebyshev closure, and the spline and
 * Chebyshev builders.
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

void
_nc_xcor_kernel_component_state_init (ComponentState *state, NcXcorKernelComponent *comp, guint comp_idx, NcHICosmo *cosmo, guint n_l)
{
  state->comp                 = nc_xcor_kernel_component_ref (comp);
  state->comp_idx             = comp_idx;
  state->last_k_left          = G_MAXDOUBLE;
  state->last_k_right         = 0.0;
  state->left_boundary_found  = 0;
  state->right_boundary_found = 0;
  state->params.comp          = comp;
  state->params.cosmo         = cosmo;

  nc_xcor_kernel_component_get_limits (comp, cosmo,
                                       &state->chi_min, &state->chi_max,
                                       &state->k_min_hard, &state->k_max_hard);
}

/* Releases the component references taken by _nc_xcor_kernel_component_state_init(). */
void
_nc_xcor_kernel_component_states_clear (ComponentStates *comp_states)
{
  guint i;

  for (i = 0; i < comp_states->n_comp; i++)
    nc_xcor_kernel_component_clear (&comp_states->states[i].comp);

  ncm_vector_clear (&comp_states->k_vec);
  ncm_vector_clear (&comp_states->kf_vec);
}

GPtrArray *
_nc_xcor_kernel_validate_component_list (NcXcorKernel *xclk, guint n_l)
{
  NcXcorKernelClass *klass = NC_XCOR_KERNEL_GET_CLASS (xclk);
  GPtrArray *comp_list     = klass->get_component_list (xclk);

  if ((comp_list == NULL) || (comp_list->len == 0))
  {
    if (comp_list != NULL)
      g_ptr_array_unref (comp_list);

    g_error ("_nc_xcor_kernel_validate_component_list: kernel %s returned empty component list",
             G_OBJECT_TYPE_NAME (xclk));

    return NULL;
  }

  /* Hard errors rather than g_assert(): both bound writes into fixed-size
   * stack arrays (ComponentStates::last_values_*, kernel_out[][]), and asserts
   * compile out under -Dnumcosmo_assert=false, which would turn an
   * out-of-range block into a stack overflow instead of a clean abort. */
  if (n_l > MAX_ELL_BLOCK)
    g_error ("_nc_xcor_kernel_validate_component_list: kernel %s asked for %u multipoles "
             "in a single block, but at most %d fit (NC_XCOR_KERNEL_MAX_ELL_BLOCK). "
             "Split the range into blocks.",
             G_OBJECT_TYPE_NAME (xclk), n_l, MAX_ELL_BLOCK);

  if (comp_list->len > MAX_COMP_BLOCK)
    g_error ("_nc_xcor_kernel_validate_component_list: kernel %s has %u components, "
             "but at most %d fit in a single block.",
             G_OBJECT_TYPE_NAME (xclk), comp_list->len, MAX_COMP_BLOCK);

  return comp_list; /* Caller must unref */
}

/*
 * The refinement criterion is reltol * ||f||_2 + a * ||f||_2^max, a *sum*, so
 * the larger of the two terms sets what refinement stops at and the smaller one
 * is inert. Setting one far tighter than the other therefore buys nothing while
 * still being paid for in knots, and nothing else reports it.
 *
 * The two terms are not compared directly -- one is scaled to the block's norm
 * and the other to its peak, a factor of order sqrt(n_l) apart -- so the
 * threshold is deliberately loose at two orders, triggering only where the
 * imbalance cannot be anything else. Once per kernel: closures are built per
 * ell block, and the tolerances do not change between them.
 */
static void
_nc_xcor_kernel_check_tolerance_balance (NcXcorKernel *xclk)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  const gdouble ratio       = self->reltol / self->peak_epsilon;

  if (self->tolerance_balance_warned)
    return;

  if ((ratio > 1.0e2) || (ratio < 1.0e-2))
  {
    const gboolean reltol_inert = (ratio < 1.0e-2);

    self->tolerance_balance_warned = TRUE;

    g_warning ("_nc_xcor_kernel_check_tolerance_balance: %s has reltol %.3e and "
               "peak-epsilon %.3e, %.0f orders apart. The refinement criterion adds "
               "the two, so the looser one decides where refinement stops and %s is "
               "inert -- tightening it alone cannot improve the result, and it is "
               "still paid for in spline knots. Move them together.",
               G_OBJECT_TYPE_NAME (xclk), self->reltol, self->peak_epsilon,
               fabs (log10 (ratio)), reltol_inert ? "reltol" : "peak-epsilon");
  }
}

/*
 * Checks that the closure tolerances are not below the relative tolerance of
 * the Levin integrator, and aborts with g_error() when they are: a closure
 * cannot be more accurate than the values of W_l(k) it interpolates. The
 * looser of reltol and peak-epsilon is compared, since refinement stops at
 * the looser of the two. Other integrators are not checked.
 */
void
_nc_xcor_kernel_check_integrator_tolerance (NcXcorKernel *xclk, NcmSBesselIntegrator *sbi)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  const gdouble fit_tol     = GSL_MAX (self->reltol, self->peak_epsilon);
  gdouble integrator_reltol;

  if ((sbi == NULL) || !NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbi))
    return;

  {
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

    integrator_reltol = GSL_MAX (ncm_sbessel_integrator_levin_get_reltol (sbilv),
                                 ncm_sbessel_integrator_levin_get_cheb_reltol (sbilv));
  }

  if (fit_tol < integrator_reltol)
    g_error ("_nc_xcor_kernel_check_integrator_tolerance: kernel %s builds its closure of W_l(k) "
             "to %.17g but the integrator computes W_l(k) only to %.17g (the looser of its "
             "reltol and cheb-reltol). A closure cannot resolve what the values do not carry. "
             "Loosen NcXcorKernel:reltol and NcXcorKernel:peak-epsilon to at least %.17g, "
             "or construct the integrator with tighter tolerances.",
             G_OBJECT_TYPE_NAME (xclk), fit_tol, integrator_reltol, integrator_reltol);
}

/*
 * NcmSpectralFBatch takes user data first; the compute functions here take k
 * first, as NcmFunctionSampleSetFunc does. One adapter rather than changing
 * either.
 */
typedef struct _ChebCompute
{
  void (*compute_func) (const gdouble, NcmVector *, gpointer);

  ComponentStates *comp_states;
  NcmFunctionSampleSet *expansion; /* the values of W_l(k) computed during the domain expansion */
  GPtrArray *ends;                 /* per split depth, the 3 by n_l W_l at a, the midpoint and b */
  GPtrArray *ends_rows;            /* per split depth, the three row views of its ends matrix */
} ChebCompute;

static void
_cheb_compute_call (gpointer user_data, gdouble k, NcmVector *y)
{
  ChebCompute *compute = (ChebCompute *) user_data;

  compute->compute_func (k, y, compute->comp_states);
}

static gint
_nc_xcor_kernel_cmp_gdouble (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

/*
 * Collects into cuts the k strictly inside (k_min, k_max) at which the block
 * W_l(k) is discontinuous, sorted ascending with exact duplicates removed:
 * the truncation boundary of every component whose boundary was found, under
 * Limber the band edges nu / chi_max and nu / chi_min of every component for
 * every multipole of the block, and under Limber in u, for every component with
 * a Bessel derivative and every multipole, u = (1 + 1 / nu) / chi_max, where its
 * j_{l+1} term stops (see _component_states_compute_limber_u()). Each cut becomes a panel edge, with the
 * component or multipole off on the outer side, so that every panel is a
 * smooth function; a polynomial interpolating across a jump does not converge
 * below the size of the jump. The boundaries are fixed before this runs: only
 * the domain expansion moves them, and it has finished.
 *
 */
static void
_component_states_collect_cuts (ComponentStates *comp_states, gdouble k_min, gdouble k_max, GArray *cuts)
{
  guint ci;

  g_array_set_size (cuts, 0);

  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    const ComponentState *state = &comp_states->states[ci];

    if ((state->left_boundary_found >= comp_states->adaptive_boundary_tries) &&
        (state->last_k_left > k_min) && (state->last_k_left < k_max))
      g_array_append_val (cuts, state->last_k_left);

    if ((state->right_boundary_found >= comp_states->adaptive_boundary_tries) &&
        (state->last_k_right > k_min) && (state->last_k_right < k_max))
      g_array_append_val (cuts, state->last_k_right);

    /* Under Limber the window steps to zero at every multipole's band edge,
     * see _component_states_compute_limber_u(). */
    if (comp_states->is_limber)
    {
      guint j;

      for (j = 0; j < comp_states->n_l; j++)
      {
        if ((state->k_min_limber_ell[j] > k_min) && (state->k_min_limber_ell[j] < k_max))
          g_array_append_val (cuts, state->k_min_limber_ell[j]);

        if ((state->k_max_limber_ell[j] > k_min) && (state->k_max_limber_ell[j] < k_max))
          g_array_append_val (cuts, state->k_max_limber_ell[j]);
      }
    }

    /* The j_{l+1} term of a derivative component is dropped once its peak
     * (nu + 1) / k passes chi_max, a step at u = (1 + 1 / nu) / chi_max inside the band. */
    if (comp_states->is_limber && comp_states->in_u && (nc_xcor_kernel_component_get_bessel_deriv (state->comp) > 0))
    {
      guint j;

      for (j = 0; j < comp_states->n_l; j++)
      {
        const gdouble nu     = comp_states->lmin + j + 0.5;
        const gdouble u_step = (1.0 + 1.0 / nu) / state->chi_max;

        if ((u_step > k_min) && (u_step < k_max))
          g_array_append_val (cuts, u_step);
      }
    }
  }

  g_array_sort (cuts, _nc_xcor_kernel_cmp_gdouble);

  /* Two components stopping at the same k share the double. */
  {
    guint w = 0;

    for (ci = 0; ci < cuts->len; ci++)
      if ((w == 0) || (g_array_index (cuts, gdouble, ci) > g_array_index (cuts, gdouble, w - 1)))
        g_array_index (cuts, gdouble, w++) = g_array_index (cuts, gdouble, ci);

    g_array_set_size (cuts, w);
  }
}

/*
 * Fills edges with the initial panel edges of the Chebyshev closure on
 * [k_min, k_max], ascending, from k_min to k_max: the cuts of
 * _component_states_collect_cuts() and the geometric grid exp (j / panels_per_efold),
 * j integer, inside the range. The grid is anchored at 1 rather than at k_min
 * so that two closures share their grid edges wherever their ranges overlap,
 * and the exact integration then restricts neither of them on those cells. A
 * grid point closer than a quarter of the grid spacing, in ln k, to an edge
 * already present is dropped. With panels_per_efold zero the edges are the
 * cuts alone.
 */
static void
_nc_xcor_kernel_cheb_panel_edges (ComponentStates *comp_states, gdouble k_min, gdouble k_max,
                                  gdouble panels_per_efold, GArray *edges)
{
  GArray *cuts = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint i;

  _component_states_collect_cuts (comp_states, k_min, k_max, cuts);

  g_array_set_size (edges, 0);
  g_array_append_val (edges, k_min);
  g_array_append_vals (edges, cuts->data, cuts->len);
  g_array_append_val (edges, k_max);

  if (panels_per_efold > 0.0)
  {
    const gdouble h    = 1.0 / panels_per_efold;
    const glong j_min  = (glong) floor (log (k_min) / h) + 1;
    const glong j_max  = (glong) ceil (log (k_max) / h) - 1;
    const guint n_cuts = edges->len;
    glong j;

    for (j = j_min; j <= j_max; j++)
    {
      const gdouble k = exp (j * h);
      gboolean near   = FALSE;

      for (i = 0; i < n_cuts; i++)
      {
        if (fabs (log (k / g_array_index (edges, gdouble, i))) < 0.25 * h)
        {
          near = TRUE;
          break;
        }
      }

      if (!near)
        g_array_append_val (edges, k);
    }

    g_array_sort (edges, _nc_xcor_kernel_cmp_gdouble);
  }

  g_array_unref (cuts);
}

/* Default panel-order cap for Chebyshev closures. */
#define NC_XCOR_KERNEL_CHEB_PANEL_K_CAP (5)
#define NC_XCOR_KERNEL_CHEB_MIN_PANEL_FRAC (1.0e-6)

/*
 * Whether an accepted panel on (a, b) reproduces the values of W_l(k) computed
 * during the domain expansion at k strictly inside it.
 *
 * The doubling test compares two expansions of what the nodes saw, and a
 * feature narrower than the first grid's spacing is seen by neither: the
 * nodes of a 9-point grid on a domain five decades wide start past a peak
 * that sits at a thousandth of it, both levels read a small smooth function,
 * and the test passes. Measured on NcXcorKernelAnalyticMulti with
 * adaptive-epsilon 1e-11: the closure's peak came out at 4.9e-5 against 16.2,
 * with no diagnostic. The k of the domain expansion are independent of the
 * nodes, lie on every component's scale by construction, and W_l(k) is
 * already computed there.
 *
 * Values at the edges are excluded: a cut edge is two-valued there. The
 * factor of ten separates a missed feature from the acceptance test's own
 * slack, which is an l2 statement over the coefficients, not a pointwise one.
 */
static gboolean
_nc_xcor_kernel_cheb_panel_matches_expansion (const ChebPanel *panel, NcmFunctionSampleSet *samples,
                                              guint n_l, gdouble reltol, gdouble abstol)
{
  NcmFunctionSampleSetIter *iter = NULL;
  gboolean ok                    = TRUE;

  ncm_function_sample_set_iter_begin (samples, &iter);

  for ( ; ok && ncm_function_sample_set_iter_is_valid (iter); ncm_function_sample_set_iter_next (iter))
  {
    const gdouble x = ncm_function_sample_set_iter_get_x (iter);

    if ((x > panel->a) && (x < panel->b))
    {
      NcmVector *y    = ncm_function_sample_set_iter_get_y (iter);
      const gdouble s = ncm_spectral_x_to_s (panel->a, panel->b, x);
      guint c;

      for (c = 0; c < n_l; c++)
      {
        const gdouble yc    = ncm_vector_get (y, c);
        const gdouble value = _nc_xcor_kernel_cheb_panel_eval_one (panel, c, s);

        if (fabs (value - yc) > 10.0 * (reltol * fabs (yc) + abstol))
        {
          ok = FALSE;
          break;
        }
      }
    }
  }

  ncm_function_sample_set_iter_free (iter);

  return ok;
}

/*
 * Expands on [a, b], bisecting where the capped order does not converge or
 * where the accepted panel misses a value of W_l(k) computed during the domain
 * expansion.
 * Panels are appended in ascending order, so the result is contiguous.
 */
static void
_nc_xcor_kernel_cheb_split (NcmSpectral *spectral, ChebCompute *compute, guint n_l,
                            gdouble a, gdouble b, gdouble reltol, gdouble abstol,
                            guint level_min, guint k_cap, guint depth,
                            NcmVector *W_a, NcmVector *W_b, GArray *panels)
{
  NcmMatrix *coeffs = NULL;
  NcmMatrix *ends;
  guint k_ord;

  /* One ends matrix per depth, made the first time the depth is reached: a split hands
   * rows 0, 1 to its lower half and rows 1, 2 to its upper half, and the lower half's
   * own splits write one depth further down, so the rows survive for the upper half. */
  if (compute->ends->len <= depth)
  {
    NcmMatrix *m = ncm_matrix_new (3, n_l);
    guint r;

    g_ptr_array_add (compute->ends, m);

    for (r = 0; r < 3; r++)
      g_ptr_array_add (compute->ends_rows, ncm_matrix_get_row (m, r));
  }

  ends = g_ptr_array_index (compute->ends, depth);

  compute->comp_states->panel_mid = 0.5 * (a + b);

  k_ord = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (
    spectral, _cheb_compute_call, n_l, a, b, level_min,
    k_cap, reltol, abstol, FALSE, W_a, W_b, ends, &coeffs, compute);

  if (k_ord > 0)
  {
    ChebPanel panel = { a, b, coeffs, (1u << k_ord) + 1u };

    if (_nc_xcor_kernel_cheb_panel_matches_expansion (&panel, compute->expansion, n_l, reltol, abstol))
    {
      g_array_append_val (panels, panel);

      return;
    }

    ncm_matrix_clear (&coeffs);
  }

  {
    const gdouble mid = 0.5 * (a + b);

    /* A panel that will not converge however far it is split is not a
     * resolution problem; refusing to bisect past a fraction of the domain
     * turns a hang into a diagnosable expansion. */
    if ((b - a) < NC_XCOR_KERNEL_CHEB_MIN_PANEL_FRAC * b)
    {
      const guint k_forced = ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (
        spectral, _cheb_compute_call, n_l, a, b, level_min,
        k_cap, reltol, abstol, TRUE, W_a, W_b, NULL, &coeffs, compute);
      ChebPanel panel = { a, b, coeffs, (1u << k_forced) + 1u };

      g_array_append_val (panels, panel);

      return;
    }

    {
      NcmVector *W_lo  = g_ptr_array_index (compute->ends_rows, 3 * depth + 0);
      NcmVector *W_mid = g_ptr_array_index (compute->ends_rows, 3 * depth + 1);
      NcmVector *W_hi  = g_ptr_array_index (compute->ends_rows, 3 * depth + 2);

      _nc_xcor_kernel_cheb_split (spectral, compute, n_l, a, mid, reltol, abstol, level_min, k_cap, depth + 1, W_lo, W_mid, panels);
      _nc_xcor_kernel_cheb_split (spectral, compute, n_l, mid, b, reltol, abstol, level_min, k_cap, depth + 1, W_mid, W_hi, panels);
    }
  }
}

static gboolean
_is_new_k (gdouble k, gdouble *arr, guint n)
{
  guint m;

  for (m = 0; m < n; m++)
    if (gsl_fcmp (k, arr[m], 1e-3) == 0)
      return FALSE;

  return TRUE;
}

static void
_component_states_compute_k_seeds (ComponentStates *comp_states, GArray *k_seeds)
{
  gdouble log_k_center_sum = 0.0;
  gdouble k_comp_scales[MAX_COMP_BLOCK];
  gdouble k_min_soft = G_MAXDOUBLE;
  gdouble k_max_soft = 0.0;
  gdouble k_center;
  guint i, j;

  /* Clear the array and prepare for new values */
  g_array_set_size (k_seeds, 0);

  /* Compute soft limits and component scales */
  for (i = 0; i < comp_states->n_comp; i++)
  {
    ComponentState *state = &comp_states->states[i];
    gdouble ln_k_scale    = 0.0;
    gdouble n_k           = 0.0;

    for (j = 0; j < comp_states->n_l; j++)
    {
      const gdouble nu       = comp_states->lmin + j + 0.5;
      const gdouble k_max_ij = nc_xcor_kernel_component_eval_k_max (state->comp, nu) / (comp_states->in_u ? nu : 1.0);
      gdouble k_upper_ij     = k_max_ij * 1.01;
      gdouble k_lower_ij     = k_max_ij * 0.99;

      if ((k_lower_ij > comp_states->k_min_hard) && (k_upper_ij < comp_states->k_max_hard))
      {
        ln_k_scale += log (k_max_ij);
        n_k        += 1.0;
        k_max_soft  = GSL_MAX (k_max_soft, k_upper_ij);
        k_min_soft  = GSL_MIN (k_min_soft, k_lower_ij);
      }
    }

    if (n_k > 0.0)
    {
      log_k_center_sum += ln_k_scale / n_k;
      k_comp_scales[i]  = exp (ln_k_scale / n_k);
    }
    else
    {
      /* If no valid k_max was found for this component, use the geometric mean of hard limits */
      k_comp_scales[i]  = sqrt (comp_states->k_min_hard * comp_states->k_max_hard);
      log_k_center_sum += log (k_comp_scales[i]);
      k_max_soft        = GSL_MAX (k_max_soft, k_comp_scales[i] * (1.0 + 1.0e-5));
      k_min_soft        = GSL_MIN (k_min_soft, k_comp_scales[i] * (1.0 - 1.0e-5));
    }
  }

  k_center = exp (log_k_center_sum / comp_states->n_comp);

  if ((k_center < k_min_soft) || (k_center > k_max_soft))
    k_center = (k_min_soft + k_max_soft) / 2.0;

  g_assert_cmpfloat (k_min_soft, <, k_max_soft);
  g_assert_cmpfloat (comp_states->k_min_hard, <=, k_min_soft);
  g_assert_cmpfloat (k_max_soft, <=, comp_states->k_max_hard);
  g_assert_cmpfloat (k_min_soft, <, k_center);
  g_assert_cmpfloat (k_center, <, k_max_soft);

  /* Add k_min_soft */
  g_array_append_val (k_seeds, k_min_soft);

  /* Add k_center if unique */
  if (_is_new_k (k_center, (gdouble *) k_seeds->data, k_seeds->len))
    g_array_append_val (k_seeds, k_center);

  /* Add unique component scales */
  for (i = 0; i < comp_states->n_comp; i++)
  {
    if (_is_new_k (k_comp_scales[i], (gdouble *) k_seeds->data, k_seeds->len))
      g_array_append_val (k_seeds, k_comp_scales[i]);
  }

  /* Add k_max_soft */
  if (_is_new_k (k_max_soft, (gdouble *) k_seeds->data, k_seeds->len))
    g_array_append_val (k_seeds, k_max_soft);
}

/*
 * Builds the closure as a Chebyshev series rather than a refined spline.
 *
 * The domain is found exactly as the spline path finds it -- the seeds and
 * ncm_function_sample_set_expand_domain() are shared, since where W is
 * negligible is a property of the kernel and not of the representation. What
 * changes is everything after: instead of bisecting until the acceptance test is
 * met, the whole ell block is expanded on one Chebyshev-Lobatto grid per
 * panel, doubling the order until every multipole's coefficients converge.
 *
 * The panels start from the component boundaries the expansion found, since
 * the block W_l(k) jumps there (see _component_states_collect_cuts), and a
 * converged panel is checked against the values computed during the domain
 * expansion before it is kept (see _nc_xcor_kernel_cheb_panel_matches_expansion).
 */
NcXcorKernelIntegrand *
_nc_xcor_kernel_build_cheb_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax,
                                      ComponentStates *comp_states,
                                      void (*compute_func) (const gdouble, NcmVector *, gpointer),
                                      const gdouble reltol, const gdouble peak_epsilon)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  ChebIntegrandData *cid    = g_new0 (ChebIntegrandData, 1);
  const guint n_l           = lmax - lmin + 1;

  {
    NcmFunctionSampleSet *fss = ncm_function_sample_set_new (n_l);
    GArray *k_seeds           = g_array_new (FALSE, FALSE, sizeof (gdouble));
    ChebCompute compute       = {
      compute_func, comp_states, fss,
      g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_matrix_free),
      g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free)
    };
    NcmMatrix *closure_error = NULL;
    gdouble abstol;
    guint i;

    _component_states_compute_k_seeds (comp_states, k_seeds);

    for (i = 0; i < k_seeds->len; i++)
    {
      const gdouble k_seed = g_array_index (k_seeds, gdouble, i);

      ncm_function_sample_set_add_old_func (fss, k_seed, compute_func, comp_states);
    }

    g_array_unref (k_seeds);

    ncm_function_sample_set_expand_domain (
      fss,
      compute_func,
      comp_states->k_min_hard,
      comp_states->k_max_hard,
      self->expansion_factor,
      comp_states->epsilon,
      self->max_border_expansions,
      comp_states->adaptive_boundary_tries,
      comp_states
    );

    cid->k_min = ncm_function_sample_set_get_x_min (fss);
    cid->k_max = ncm_function_sample_set_get_x_max (fss);

    /* Same meaning the spline path gives it: a floor scaled to the smallest of
     * the block's peaks, so a sub-dominant multipole is not held to a
     * tolerance relative to its neighbours. */
    abstol = ncm_function_sample_set_get_absmaxF_min (fss) * peak_epsilon;

    cid->panels = g_array_new (FALSE, FALSE, sizeof (ChebPanel));
    cid->edges  = g_array_new (FALSE, FALSE, sizeof (gdouble));

    /* Expand each initial panel, bisecting where it does not converge. In
     * panel mode the compute function decides at the panel midpoint which
     * components and multipoles are on, so every discontinuity has to be an
     * initial edge, which _nc_xcor_kernel_cheb_panel_edges() guarantees. */
    {
      NcmSpectral **spectral = _nc_xcor_kernel_spectral_get ();
      GArray *edges0         = g_array_new (FALSE, FALSE, sizeof (gdouble));
      const guint k_cap      = (self->panel_order_cap == 0) ?
                               NC_XCOR_KERNEL_CHEB_PANEL_K_CAP : self->panel_order_cap;

      if (self->panel_level_min >= k_cap)
        g_error ("_nc_xcor_kernel_build_cheb_integrand: kernel %s has panel-level-min %u, "
                 "not below its panel-order-cap %u.",
                 G_OBJECT_TYPE_NAME (xclk), self->panel_level_min, k_cap);

      /* The grid serves the Limber closure in u, where the window does not
       * oscillate. A non-Limber W_l(k) oscillates with the period pi / chi_max,
       * and a grid panel whose width lets it alias on the first levels of the
       * ladder is accepted wrongly, which the bisection from the whole domain
       * has not shown on the certified windows. */
      _nc_xcor_kernel_cheb_panel_edges (comp_states, cid->k_min, cid->k_max,
                                        comp_states->in_u ? self->panels_per_efold : 0.0, edges0);

      comp_states->panel_mode = TRUE;

      /* An initial edge can be a cut, where W_l(k) jumps, so an initial panel starts
       * with no known ends; the halves of a split share them with their parent. */
      for (i = 0; i + 1 < edges0->len; i++)
        _nc_xcor_kernel_cheb_split (*spectral, &compute, n_l,
                                    g_array_index (edges0, gdouble, i), g_array_index (edges0, gdouble, i + 1),
                                    reltol, abstol, self->panel_level_min, k_cap, 0, NULL, NULL, cid->panels);

      comp_states->panel_mode = FALSE;

      g_array_unref (edges0);
      ncm_memory_pool_return (spectral);
    }

    ncm_function_sample_set_clear (&fss);
    g_ptr_array_unref (compute.ends_rows);
    g_ptr_array_unref (compute.ends);

    {
      const gdouble first = g_array_index (cid->panels, ChebPanel, 0).a;

      g_array_append_val (cid->edges, first);

      for (i = 0; i < cid->panels->len; i++)
        g_array_append_val (cid->edges, g_array_index (cid->panels, ChebPanel, i).b);
    }

    /* Record the interpolation error of each panel, one row per panel and one
     * column per multipole, as the sum of the absolute coefficients above the
     * order of the previous level; see NcXcorKernel:track-closure-error. */
    if (self->track_closure_error)
    {
      closure_error = ncm_matrix_new (cid->panels->len, n_l);

      for (i = 0; i < cid->panels->len; i++)
      {
        const ChebPanel *panel = &g_array_index (cid->panels, ChebPanel, i);
        const guint n_prev     = (panel->N + 1) / 2;
        guint c, j;

        for (c = 0; c < n_l; c++)
        {
          gdouble tail = 0.0;

          for (j = n_prev; j < panel->N; j++)
            tail += fabs (ncm_matrix_get (panel->coeffs, c, j));

          ncm_matrix_set (closure_error, i, c, tail);
        }
      }
    }

    cid->k_min_comp = g_new (gdouble, n_l);
    cid->k_max_comp = g_new (gdouble, n_l);
    cid->scale      = g_new (gdouble, n_l);
    cid->scaled     = comp_states->in_u;

    for (i = 0; i < n_l; i++)
      cid->scale[i] = comp_states->in_u ? (lmin + i + 0.5) : 1.0;

    /* The sampled range is in x; each multipole's k range is its scale times
     * the part of it inside the band, and the block's k range is their union. */
    for (i = 0; i < n_l; i++)
    {
      const gdouble x_min = cid->k_min;
      const gdouble x_max = cid->k_max;
      gdouble k_min_i     = x_min;
      gdouble k_max_i     = x_max;

      if (comp_states->is_limber)
      {
        gdouble band_min = G_MAXDOUBLE;
        gdouble band_max = 0.0;
        guint ci;

        for (ci = 0; ci < comp_states->n_comp; ci++)
        {
          band_min = GSL_MIN (band_min, comp_states->states[ci].k_min_limber_ell[i]);
          band_max = GSL_MAX (band_max, comp_states->states[ci].k_max_limber_ell[i]);
        }

        k_min_i = GSL_MAX (k_min_i, band_min);
        k_max_i = GSL_MIN (k_max_i, band_max);

        if (k_min_i >= k_max_i)
          k_min_i = k_max_i = x_min;
      }

      cid->k_min_comp[i] = k_min_i * cid->scale[i];
      cid->k_max_comp[i] = k_max_i * cid->scale[i];
    }

    if (cid->scaled)
    {
      cid->k_min = cid->k_min * cid->scale[0];
      cid->k_max = cid->k_max * cid->scale[n_l - 1];
    }

    cid->lmin   = lmin;
    cid->len    = n_l;
    cid->RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
    cid->cosmo  = nc_hicosmo_ref (cosmo);

    {
      NcXcorKernelIntegrand *integrand = _nc_xcor_kernel_cheb_integrand_new (n_l, cid);

      nc_xcor_kernel_integrand_set_tolerances (integrand, reltol, peak_epsilon);
      nc_xcor_kernel_integrand_set_closure_error (integrand, closure_error);
      ncm_matrix_clear (&closure_error);

      return integrand;
    }
  }
}

NcXcorKernelIntegrand *
_nc_xcor_kernel_build_spline_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax,
                                        ComponentStates *comp_states,
                                        void (*compute_func) (const gdouble, NcmVector *, gpointer),
                                        const gdouble reltol, const gdouble peak_epsilon)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  SplineIntegrandData *sid  = g_new0 (SplineIntegrandData, 1);
  const guint n_l           = lmax - lmin + 1;

  {
    NcmFunctionSampleSet *fss = ncm_function_sample_set_new (n_l);
    NcmSpline *spline         = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    GArray *k_seeds           = g_array_new (FALSE, FALSE, sizeof (gdouble));
    NcmMatrix *closure_error  = NULL;
    guint i;

    _nc_xcor_kernel_check_tolerance_balance (xclk);
    ncm_function_sample_set_set_track_residual (fss, self->track_closure_error);

    /* Compute the k seeds, the first k at which W is evaluated. Local to the call, not kernel
     * state: one kernel may be evaluated concurrently for different ell
     * blocks, each with its own integrator. */
    _component_states_compute_k_seeds (comp_states, k_seeds);

    /* Evaluate W at every seed */
    for (i = 0; i < k_seeds->len; i++)
    {
      const gdouble k_seed = g_array_index (k_seeds, gdouble, i);

      ncm_function_sample_set_add_old_func (fss, k_seed, compute_func, comp_states);
    }

    g_array_unref (k_seeds);

    {
      /* Domain expansion */
      ncm_function_sample_set_expand_domain (
        fss,
        compute_func,
        comp_states->k_min_hard,
        comp_states->k_max_hard,
        self->expansion_factor,
        comp_states->epsilon,
        self->max_border_expansions,
        comp_states->adaptive_boundary_tries,
        comp_states
      );
    }

    ncm_function_sample_set_mark_all_old (fss);
    ncm_function_sample_set_reset_interval_ok (fss);
    {
      const gdouble max_absF_total = ncm_function_sample_set_get_absmaxF_min (fss);

      ncm_function_sample_set_adaptive_midpoint (
        fss, compute_func,
        reltol, max_absF_total * peak_epsilon, self->max_iter, 1,
        spline, comp_states
      );
    }

    closure_error   = ncm_function_sample_set_get_residuals (fss);
    sid->spline_vec = ncm_function_sample_set_to_spline_vec (fss, spline);
    sid->k_min      = ncm_function_sample_set_get_x_min (fss);
    sid->k_max      = ncm_function_sample_set_get_x_max (fss);
    sid->k_min_comp = g_new (gdouble, n_l);
    sid->k_max_comp = g_new (gdouble, n_l);

    /* Per-multipole support within the block's shared domain. Only the Limber
     * branch confines a multipole to a band of its own; outside it the window
     * is zero, so the band edge falling inside the shared domain is a step.
     * The band is taken over all components, since the window is their sum. */
    for (i = 0; i < n_l; i++)
    {
      gdouble k_min_i = sid->k_min;
      gdouble k_max_i = sid->k_max;

      if (comp_states->is_limber)
      {
        gdouble band_min = G_MAXDOUBLE;
        gdouble band_max = 0.0;
        guint ci;

        for (ci = 0; ci < comp_states->n_comp; ci++)
        {
          band_min = GSL_MIN (band_min, comp_states->states[ci].k_min_limber_ell[i]);
          band_max = GSL_MAX (band_max, comp_states->states[ci].k_max_limber_ell[i]);
        }

        k_min_i = GSL_MAX (k_min_i, band_min);
        k_max_i = GSL_MIN (k_max_i, band_max);

        /* A band disjoint from the closure's domain leaves nothing to integrate;
         * report the empty range as the domain's lower edge rather than an
         * inverted one. */
        if (k_min_i >= k_max_i)
          k_min_i = k_max_i = sid->k_min;
      }

      sid->k_min_comp[i] = k_min_i;
      sid->k_max_comp[i] = k_max_i;
    }

    sid->lmin        = lmin;
    sid->len         = n_l;
    sid->RH_Mpc      = nc_hicosmo_RH_Mpc (cosmo);
    sid->cosmo       = nc_hicosmo_ref (cosmo);
    sid->eval_result = ncm_vector_new (n_l);

    ncm_function_sample_set_clear (&fss);
    ncm_spline_free (spline);

    {
      NcXcorKernelIntegrand *integrand = _nc_xcor_kernel_spline_integrand_new (n_l, sid);

      nc_xcor_kernel_integrand_set_tolerances (integrand, reltol, peak_epsilon);
      nc_xcor_kernel_integrand_set_closure_error (integrand, closure_error);
      ncm_matrix_clear (&closure_error);

      return integrand;
    }
  }
}


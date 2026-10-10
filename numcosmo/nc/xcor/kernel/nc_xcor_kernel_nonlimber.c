/***************************************************************************
 *            nc_xcor_kernel_nonlimber.c
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
 * The non-Limber path of NcXcorKernel: the component states of a block, the
 * evaluation of W_l(k) through the spherical Bessel integrator, and the
 * non-Limber integrand builder.
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

static gdouble
_nc_xcor_kernel_component_kernel_integ (gpointer params, gdouble chi, gdouble k)
{
  const ComponentParams *nlcp = (const ComponentParams *) params;
  const gdouble kernel        = nc_xcor_kernel_component_eval_kernel (nlcp->comp, nlcp->cosmo, chi, k);

  return kernel / (k * sqrt (k));
}

static ComponentStates
_component_states_init_non_limber (NcXcorKernel *xclk, gint lmin, guint n_l,
                                   GPtrArray *comp_list, NcHICosmo *cosmo,
                                   NcmSBesselIntegrator *sbi)
{
  NcXcorKernelPrivate *self   = _nc_xcor_kernel_get_private (xclk);
  const guint n_comp          = comp_list->len;
  ComponentStates comp_states = {
    .xclk                    = xclk,
    .sbi                     = sbi,
    .k_min_hard              = 0.0,
    .k_max_hard              = G_MAXDOUBLE,
    .l2_norm                 = 0.0,
    .n_comp                  = n_comp,
    .lmin                    = lmin,
    .n_l                     = n_l,
    .epsilon                 = self->adaptive_epsilon,
    .adaptive_boundary_tries = self->adaptive_boundary_tries,
    .is_limber               = FALSE,
    .in_u                    = FALSE
  };
  guint i;

  /* Initialize each component state */
  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];

    _nc_xcor_kernel_component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

    /* Compute global hard limits as intersection of all component limits */
    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard);

    state->group = i;
  }

  g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

  /*
   * Group the components for the k-range truncation. Components whose chi supports
   * touch or overlap are tested together on the norm of their sum and truncated at the
   * same k; grouping is transitive and implemented by relabelling. Testing such
   * components separately fails when the window is continuous across the shared edge:
   * the endpoint terms of the two radial integrals at that edge cancel in the sum, and
   * cutting one component alone leaves the other's term uncancelled over its whole
   * remaining range.
   */
  for (i = 0; i < n_comp; i++)
  {
    guint j;

    for (j = 0; j < i; j++)
    {
      ComponentState *si = &comp_states.states[i];
      ComponentState *sj = &comp_states.states[j];
      const gdouble tol  = 1.0e-9 * GSL_MAX (si->chi_max, sj->chi_max);

      if ((si->chi_min <= sj->chi_max + tol) && (sj->chi_min <= si->chi_max + tol))
      {
        const guint gi = si->group;
        const guint gj = sj->group;
        guint m;

        for (m = 0; m < n_comp; m++)
          if (comp_states.states[m].group == gj)
            comp_states.states[m].group = gi;
      }
    }
  }

  return comp_states;
}

static void
_component_states_compute_non_limber (const gdouble k, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  gdouble kernel_out[MAX_COMP_BLOCK][MAX_ELL_BLOCK];
  gdouble group_sum[MAX_COMP_BLOCK][MAX_ELL_BLOCK];
  gboolean integrated[MAX_COMP_BLOCK];
  gdouble l2_norm = 0.0;
  guint ci, i;

  memset (group_sum, 0, sizeof (group_sum));

  /* Compute kernel for each component */
  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    NcmVector *integ_result              = ncm_vector_new_data_static (kernel_out[ci], comp_states->n_l, 1);
    ComponentState *state                = &comp_states->states[ci];
    const gboolean right_boundary_found  = state->right_boundary_found >= comp_states->adaptive_boundary_tries;
    const gboolean left_boundary_found   = state->left_boundary_found >= comp_states->adaptive_boundary_tries;
    const gdouble k_side                 = comp_states->panel_mode ? comp_states->panel_mid : k;
    const gboolean within_left_boundary  = !left_boundary_found || (k_side >= state->last_k_left);
    const gboolean within_right_boundary = !right_boundary_found || (k_side <= state->last_k_right);

    integrated[ci] = within_left_boundary && within_right_boundary;

    if (integrated[ci])
    {
      /* Exact integration within boundaries */
      ncm_sbessel_integrator_integrate_deriv (
        comp_states->sbi, _nc_xcor_kernel_component_kernel_integ,
        state->chi_min, state->chi_max, k,
        nc_xcor_kernel_component_get_bessel_deriv (state->comp),
        integ_result, &state->params
      );

      for (i = 0; i < comp_states->n_l; i++)
      {
        const gdouble prefactor = nc_xcor_kernel_component_eval_prefactor (
          state->comp, state->params.cosmo, k, comp_states->lmin + i
        );

        kernel_out[ci][i]          *= prefactor * k * sqrt (k);
        group_sum[state->group][i] += kernel_out[ci][i];
      }
    }
    else
    {
      /* Exponential tail extrapolation beyond boundaries.
       *
       * At DECAY_RATE this falls by e^-100 within a *relative* 1e-8 of the
       * boundary, so it adds nothing to the integral. Its purpose is
       * continuity: the block W_l(k) has no jump at the boundary, so one
       * cubic spline can span it. The first derivative still jumps, by the
       * component's value at the boundary.
       *
       * A Chebyshev panel needs none of that: the boundary is one of its
       * edges, so the whole panel is outside and the component is zero on it,
       * including at the shared edge. */
      if (comp_states->panel_mode)
      {
        for (i = 0; i < comp_states->n_l; i++)
          kernel_out[ci][i] = 0.0;
      }
      else if (!within_right_boundary)
      {
        for (i = 0; i < comp_states->n_l; i++)
        {
          const gdouble val              = state->last_values_right[i];
          const gdouble delta_k          = k - state->last_k_right;
          const gdouble decay_rate       = DECAY_RATE;
          const gdouble val_extrapolated = val * exp (-decay_rate * delta_k / state->last_k_right);

          kernel_out[ci][i] = val_extrapolated;
        }
      }
      else
      {
        for (i = 0; i < comp_states->n_l; i++)
        {
          const gdouble val              = state->last_values_left[i];
          const gdouble delta_k          = state->last_k_left - k;
          const gdouble decay_rate       = DECAY_RATE;
          const gdouble val_extrapolated = val * exp (-decay_rate * delta_k / state->last_k_left);

          kernel_out[ci][i] = val_extrapolated;
        }
      }
    }

    ncm_vector_clear (&integ_result);
  }

  /*
   * Boundary tracking. For each component still integrated, record the
   * outermost k on each side and count the consecutive outermost k at which
   * the norm of its group's sum is below epsilon times the block's running
   * maximum; adaptive-boundary-tries such k mark the boundary on that side.
   * Skipped in panel mode: the boundaries are fixed by then and the panel
   * nodes do not move outward.
   */
  if (!comp_states->panel_mode)
  {
    for (ci = 0; ci < comp_states->n_comp; ci++)
    {
      ComponentState *state  = &comp_states->states[ci];
      gdouble group_l2_norm2 = 0.0;
      gboolean below_epsilon;

      if (!integrated[ci])
        continue;

      for (i = 0; i < comp_states->n_l; i++)
        group_l2_norm2 += gsl_pow_2 (group_sum[state->group][i]);

      below_epsilon = group_l2_norm2 < gsl_pow_2 (comp_states->epsilon * comp_states->l2_norm);

      if (k > state->last_k_right)
      {
        state->last_k_right = k;

        for (i = 0; i < comp_states->n_l; i++)
          state->last_values_right[i] = kernel_out[ci][i];

        if (below_epsilon)
          state->right_boundary_found++;
        else
          state->right_boundary_found = 0;
      }
      else if (k < state->last_k_left)
      {
        state->last_k_left = k;

        for (i = 0; i < comp_states->n_l; i++)
          state->last_values_left[i] = kernel_out[ci][i];

        if (below_epsilon)
          state->left_boundary_found++;
        else
          state->left_boundary_found = 0;
      }
    }
  }

  /* Sum contributions from all components and compute total L2 norm */
  g_assert_cmpuint (ncm_vector_len (y), ==, comp_states->n_l);

  for (i = 0; i < comp_states->n_l; i++)
  {
    gdouble sum = 0.0;

    for (ci = 0; ci < comp_states->n_comp; ci++)
      sum += kernel_out[ci][i];

    ncm_vector_set (y, i, sum);
    l2_norm += sum * sum;
  }

  l2_norm = sqrt (l2_norm);

  /* Update reference L2 norm for convergence testing */
  if (l2_norm > comp_states->l2_norm)
    comp_states->l2_norm = l2_norm;
}

NcXcorKernelIntegrand *
_nc_xcor_kernel_build_non_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi, NcXcorKernelClosure closure_type)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  const guint n_l           = lmax - lmin + 1;
  GPtrArray *comp_list      = _nc_xcor_kernel_validate_component_list (xclk, n_l);

  if (comp_list == NULL)
    return NULL;

  if (sbi == NULL)
  {
    g_ptr_array_unref (comp_list);
    g_error ("_nc_xcor_kernel_build_non_limber_integrand: no integrator for kernel %s. "
             "Either set the 'integrator' property or pass one to "
             "nc_xcor_kernel_get_eval_vectorized_full().",
             G_OBJECT_TYPE_NAME (xclk));

    return NULL;
  }

  _nc_xcor_kernel_check_integrator_tolerance (xclk, sbi);
  ncm_sbessel_integrator_set_ell_range (sbi, lmin, lmax);

  /* Initialize with standard (non-Limber) limits */
  {
    ComponentStates comp_states      = _component_states_init_non_limber (xclk, lmin, n_l, comp_list, cosmo, sbi);
    NcXcorKernelIntegrand *integrand = NULL;

    g_ptr_array_unref (comp_list);

    if (closure_type == NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV)
      integrand = _nc_xcor_kernel_build_cheb_integrand (xclk, cosmo, lmin, lmax,
                                                        &comp_states,
                                                        _component_states_compute_non_limber,
                                                        self->reltol, self->peak_epsilon);
    else
      integrand = _nc_xcor_kernel_build_spline_integrand (xclk, cosmo, lmin, lmax,
                                                          &comp_states,
                                                          _component_states_compute_non_limber,
                                                          self->reltol, self->peak_epsilon);

    _nc_xcor_kernel_component_states_clear (&comp_states);

    return integrand;
  }
}


/***************************************************************************
 *            nc_xcor_kernel_limber.c
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
 * The Limber path of NcXcorKernel: the Limber window of a component, the
 * component states of a block sampled in k (spline closure) or in u = k / nu
 * (Chebyshev closure), their evaluation, and the Limber integrand builder.
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

/*
 * Returns the Limber approximation of one component's contribution to W_l(k),
 * prefactor(k, l) int K(chi, k) j_l^(d)(k chi) dchi with
 * int K(chi, k) j_l(k chi) dchi ~ sqrt(pi / (2 nu)) K(nu / k, k) / k, nu = l + 1/2.
 * For d = 1, 2 the derivative is written as a combination of j_l and j_{l+1}
 * and each term is approximated at its own peak, nu / k and (nu + 1) / k. The
 * second peak lies above the first by 1 / k and is set to zero when it exceeds
 * chi_max, outside the component's support.
 */
static gdouble
_component_limber_eval (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi_max, gdouble k, gint l)
{
  const guint deriv       = nc_xcor_kernel_component_get_bessel_deriv (comp);
  const gdouble nu        = l + 0.5;
  const gdouble prefactor = nc_xcor_kernel_component_eval_prefactor (comp, cosmo, k, l);
  const gdouble peak_l    = sqrt (M_PI / (2.0 * nu)) * nc_xcor_kernel_component_eval_kernel (comp, cosmo, nu / k, k);
  gdouble val;

  if (deriv == 0)
  {
    val = peak_l;
  }
  else
  {
    const gdouble nup      = nu + 1.0;
    const gdouble chi_p    = nup / k;
    const gdouble peak_lp1 = (chi_p <= chi_max) ?
                             sqrt (M_PI / (2.0 * nup)) * nc_xcor_kernel_component_eval_kernel (comp, cosmo, chi_p, k) :
                             0.0;

    if (deriv == 1) /* j_l' = (l/x) j_l - j_{l+1} */
      val = (l / nu) * peak_l - peak_lp1;
    else /* j_l'' = (l (l-1)/x^2 - 1) j_l + (2/x) j_{l+1} */
      val = -(2.0 * l + 0.25) / (nu * nu) * peak_l + (2.0 / nup) * peak_lp1;
  }

  return prefactor * val / k;
}

/*
 * Initializes the component states of a Limber block sampled in k, the one the
 * spline closure uses. For each component and multipole it stores the band
 * [nu / chi_max, nu / chi_min], the k where chi = nu / k lies in the component's
 * support, and the window at the two band edges, which
 * _component_states_compute_limber() continues outside the band. The hard
 * limits are the intersection of the components' k ranges, narrowed to the
 * union of the bands.
 */
static ComponentStates
_component_states_init_limber (NcXcorKernel *xclk, gint lmin, guint n_l,
                               GPtrArray *comp_list, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self   = _nc_xcor_kernel_get_private (xclk);
  const guint n_comp          = comp_list->len;
  ComponentStates comp_states = {
    .xclk                    = xclk,
    .k_min_hard              = 0.0,
    .k_max_hard              = G_MAXDOUBLE,
    .l2_norm                 = 0.0,
    .n_comp                  = n_comp,
    .lmin                    = lmin,
    .n_l                     = n_l,
    .epsilon                 = self->adaptive_epsilon,
    .adaptive_boundary_tries = self->adaptive_boundary_tries,
    .is_limber               = TRUE,
    .in_u                    = FALSE
  };
  gdouble k_min_union = G_MAXDOUBLE;
  gdouble k_max_union = 0.0;
  guint i, j;

  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];

    _nc_xcor_kernel_component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard);

    /* The band of each multipole, and the union of all bands */
    for (j = 0; j < n_l; j++)
    {
      const gdouble nu_j = lmin + j + 0.5;

      state->k_min_limber_ell[j] = nu_j / state->chi_max;
      state->k_max_limber_ell[j] = nu_j / state->chi_min;

      k_min_union = GSL_MIN (k_min_union, state->k_min_limber_ell[j]);
      k_max_union = GSL_MAX (k_max_union, state->k_max_limber_ell[j]);
    }
  }

  comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, k_min_union);
  comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, k_max_union);

  g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

  /* The window at each band edge, continued outside the band by _component_states_compute_limber() */
  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];

    for (j = 0; j < n_l; j++)
    {
      const gint l_j = lmin + j;

      state->last_values_left[j]  = _component_limber_eval (state->comp, cosmo, state->chi_max, state->k_min_limber_ell[j], l_j);
      state->last_values_right[j] = _component_limber_eval (state->comp, cosmo, state->chi_max, state->k_max_limber_ell[j], l_j);
    }
  }

  return comp_states;
}

/*
 * Initializes the component states of a Limber block sampled in u = k / nu.
 * Every multipole has the band [1 / chi_max, 1 / chi_min] of its component, so
 * the per-multipole band arrays hold the same two values for every multipole.
 * The hard limits are the intersection over the multipoles of
 * [k_min / nu, k_max / nu], in u.
 */
static ComponentStates
_component_states_init_limber_u (NcXcorKernel *xclk, gint lmin, guint n_l,
                                 GPtrArray *comp_list, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self   = _nc_xcor_kernel_get_private (xclk);
  const guint n_comp          = comp_list->len;
  const gdouble nu_min        = lmin + 0.5;
  const gdouble nu_max        = lmin + n_l - 1 + 0.5;
  ComponentStates comp_states = {
    .xclk                    = xclk,
    .k_min_hard              = 0.0,
    .k_max_hard              = G_MAXDOUBLE,
    .l2_norm                 = 0.0,
    .n_comp                  = n_comp,
    .lmin                    = lmin,
    .n_l                     = n_l,
    .epsilon                 = self->adaptive_epsilon,
    .adaptive_boundary_tries = self->adaptive_boundary_tries,
    .is_limber               = TRUE,
    .in_u                    = TRUE,
    .k_vec                   = ncm_vector_new (n_l),
    .kf_vec                  = ncm_vector_new (n_l)
  };
  gdouble u_min_union = G_MAXDOUBLE;
  gdouble u_max_union = 0.0;
  guint i, j;

  for (i = 0; i < n_comp; i++)
  {
    ComponentState *state = &comp_states.states[i];
    gdouble u_band_min, u_band_max;

    _nc_xcor_kernel_component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard / nu_min);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard / nu_max);

    u_band_min = 1.0 / state->chi_max;
    u_band_max = 1.0 / state->chi_min;

    for (j = 0; j < n_l; j++)
    {
      state->k_min_limber_ell[j] = u_band_min;
      state->k_max_limber_ell[j] = u_band_max;
    }

    u_min_union = GSL_MIN (u_min_union, u_band_min);
    u_max_union = GSL_MAX (u_max_union, u_band_max);
  }

  comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, u_min_union);
  comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, u_max_union);

  g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

  return comp_states;
}

/*
 * Evaluates the Limber block at u = k / nu: the point chi = 1 / u is shared by
 * every multipole, so each component's window is one evaluation and only the
 * wave number factor at k = nu u is per multipole. Outside its band a
 * component is zero. A derivative component adds the j_{l+1} term at the point
 * chi_p = (1 + 1 / nu) / u, one more window evaluation per multipole, dropped
 * when that point is beyond chi_max. In panel mode both choices are made at the
 * panel midpoint, so a panel ending on one of these steps takes the limit from
 * its own side.
 */
static void
_component_states_compute_limber_u (const gdouble u, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  NcDistance *dist             = nc_xcor_kernel_peek_dist (comp_states->xclk);
  NcHICosmo *cosmo             = comp_states->states[0].params.cosmo;
  const gdouble u_side         = comp_states->panel_mode ? comp_states->panel_mid : u;
  const gdouble chi            = 1.0 / u;
  NcXcorKinetic xck;
  gdouble l2_norm = 0.0;
  guint ci, i;

  gdouble W[MAX_ELL_BLOCK] = { 0.0 };

  xck.chi_z = chi;
  xck.z     = nc_distance_inv_comoving (dist, cosmo, chi);
  xck.E_z   = nc_hicosmo_E (cosmo, xck.z);

  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    ComponentState *state       = &comp_states->states[ci];
    const gboolean within_range = (u_side >= state->k_min_limber_ell[0]) && (u_side <= state->k_max_limber_ell[0]);
    const guint deriv           = nc_xcor_kernel_component_get_bessel_deriv (state->comp);
    gdouble window;

    if (!within_range)
      continue;

    window = nc_xcor_kernel_component_eval_window (state->comp, cosmo, &xck);

    for (i = 0; i < comp_states->n_l; i++)
      ncm_vector_set (comp_states->k_vec, i, (comp_states->lmin + i + 0.5) * u);

    nc_xcor_kernel_component_eval_kfactor_vec (state->comp, cosmo, &xck, comp_states->k_vec, comp_states->kf_vec);

    for (i = 0; i < comp_states->n_l; i++)
    {
      const gint l            = comp_states->lmin + i;
      const gdouble nu        = l + 0.5;
      const gdouble k         = nu * u;
      const gdouble prefactor = nc_xcor_kernel_component_eval_prefactor (state->comp, cosmo, k, l);
      const gdouble peak_l    = sqrt (M_PI / (2.0 * nu)) * window * ncm_vector_get (comp_states->kf_vec, i);
      gdouble val;

      if (deriv == 0)
      {
        val = peak_l;
      }
      else
      {
        const gdouble nup   = nu + 1.0;
        const gdouble chi_p = nup / k;
        gdouble peak_lp1    = 0.0;

        if (nup / (nu * u_side) <= state->chi_max)
        {
          NcXcorKinetic xck_p;

          xck_p.chi_z = chi_p;
          xck_p.z     = nc_distance_inv_comoving (dist, cosmo, chi_p);
          xck_p.E_z   = nc_hicosmo_E (cosmo, xck_p.z);

          peak_lp1 = sqrt (M_PI / (2.0 * nup)) *
                     nc_xcor_kernel_component_eval_window (state->comp, cosmo, &xck_p) *
                     nc_xcor_kernel_component_eval_kfactor (state->comp, cosmo, &xck_p, k);
        }

        if (deriv == 1)
          val = (l / nu) * peak_l - peak_lp1;
        else
          val = -(2.0 * l + 0.25) / (nu * nu) * peak_l + (2.0 / nup) * peak_lp1;
      }

      W[i] += prefactor * val / k;
    }
  }

  for (i = 0; i < comp_states->n_l; i++)
  {
    ncm_vector_set (y, i, W[i]);
    l2_norm += W[i] * W[i];
  }

  l2_norm = sqrt (l2_norm);

  if (l2_norm > comp_states->l2_norm)
    comp_states->l2_norm = l2_norm;
}

static void
_component_states_compute_limber (const gdouble k, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  gdouble kernel_out[MAX_COMP_BLOCK][MAX_ELL_BLOCK];
  gdouble l2_norm = 0.0;
  guint ci, i;

  /* Compute kernel for each component using Limber approximation */
  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    ComponentState *state = &comp_states->states[ci];

    /* A multipole's window is supported on its band alone: chi = nu / k must
     * fall inside the component's support. Outside the band the step to zero
     * is replaced by a Gaussian fall-off from the edge value, of relative width
     * 1e-8 (DECAY_RATE), so that the midpoint refinement of the spline closure
     * stops at a passing test inside the fall-off instead of bisecting to
     * machine precision. */
    for (i = 0; i < comp_states->n_l; i++)
    {
      const gint l                = comp_states->lmin + i;
      const gboolean within_range = (k >= state->k_min_limber_ell[i]) && (k <= state->k_max_limber_ell[i]);

      if (within_range)
      {
        kernel_out[ci][i] = _component_limber_eval (state->comp, state->params.cosmo, state->chi_max, k, l);
      }
      else if (k < state->k_min_limber_ell[i])
      {
        const gdouble delta_k = state->k_min_limber_ell[i] - k;

        kernel_out[ci][i] = state->last_values_left[i] * exp (-gsl_pow_2 (DECAY_RATE * delta_k / state->k_min_limber_ell[i]));
      }
      else
      {
        const gdouble delta_k = k - state->k_max_limber_ell[i];

        kernel_out[ci][i] = state->last_values_right[i] * exp (-gsl_pow_2 (DECAY_RATE * delta_k / state->k_max_limber_ell[i]));
      }
    }
  }

  /* Sum contributions from all components and compute total L2 norm */
  for (i = 0; i < comp_states->n_l; i++)
  {
    gdouble sum = 0.0;

    for (ci = 0; ci < comp_states->n_comp; ci++)
      sum += kernel_out[ci][i];

    ncm_vector_set (y, i, sum);
    l2_norm += sum * sum;
  }

  l2_norm = sqrt (l2_norm);

  /* Update reference L2 norm for overall convergence testing */
  if (l2_norm > comp_states->l2_norm)
    comp_states->l2_norm = l2_norm;
}

NcXcorKernelIntegrand *
_nc_xcor_kernel_build_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcXcorKernelClosure closure_type)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  const guint n_l           = lmax - lmin + 1;
  GPtrArray *comp_list      = _nc_xcor_kernel_validate_component_list (xclk, n_l);

  if (comp_list == NULL)
    return NULL;

  /* Under Limber a multipole's window is supported on the band
   * [nu / chi_max, nu / chi_min] of each component and is zero outside it. The
   * Chebyshev closure is built in u = k / nu, where the band is the same for
   * every multipole of the block: the band edges are its only cuts, and the
   * point chi = 1 / u is shared by the block (_component_states_compute_limber_u()).
   * The spline closure is built in k and confines each multipole to its band
   * at integration time (_spline_integrand_get_range_comp()). */
  if (closure_type == NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV)
  {
    ComponentStates comp_states      = _component_states_init_limber_u (xclk, lmin, n_l, comp_list, cosmo);
    NcXcorKernelIntegrand *integrand = NULL;

    g_ptr_array_unref (comp_list);

    integrand = _nc_xcor_kernel_build_cheb_integrand (xclk, cosmo, lmin, lmax,
                                                      &comp_states,
                                                      _component_states_compute_limber_u,
                                                      self->reltol, self->peak_epsilon);

    _nc_xcor_kernel_component_states_clear (&comp_states);

    return integrand;
  }
  else
  {
    ComponentStates comp_states      = _component_states_init_limber (xclk, lmin, n_l, comp_list, cosmo);
    NcXcorKernelIntegrand *integrand = NULL;

    g_ptr_array_unref (comp_list);

    integrand = _nc_xcor_kernel_build_spline_integrand (xclk, cosmo, lmin, lmax,
                                                        &comp_states,
                                                        _component_states_compute_limber,
                                                        self->reltol, self->peak_epsilon);

    _nc_xcor_kernel_component_states_clear (&comp_states);

    return integrand;
  }
}


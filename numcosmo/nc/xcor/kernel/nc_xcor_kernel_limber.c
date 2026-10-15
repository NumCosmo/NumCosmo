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
 * The Limber path of NcXcorKernel: the component states of a block sampled in
 * u = k / nu, their evaluation, and the Limber integrand builder.
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
 * Initializes the component states of a Limber block sampled in u = k / nu.
 * Every multipole has the band [1 / chi_max, 1 / chi_min] of its component, the
 * same in u, but the component's hard k limits are [k_min / nu, k_max / nu] in u,
 * one range per multipole. The per-multipole band arrays hold the band within that
 * range, where the multipole is nonzero, and the block's hard limits are the union
 * over the multipoles, [k_min / nu_max, k_max / nu_min], narrowed to the union of
 * the bands. A component with a Bessel derivative also has a j_{l+1} term, sampled in
 * x = k / (nu + 1), where its band is again [1 / chi_max, 1 / chi_min] and its hard
 * limits are [k_min / (nu + 1), k_max / (nu + 1)]; the block then has two terms.
 */
static ComponentStates
_component_states_init_limber_u (NcXcorKernel *xclk, gint lmin, guint n_l,
                                 GPtrArray *comp_list, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  const guint n_comp        = comp_list->len;
  const gdouble nu_min      = lmin + 0.5;
  const gdouble nu_max      = lmin + n_l - 1 + 0.5;
  gboolean has_deriv        = FALSE;
  guint i, j;

  for (i = 0; i < n_comp; i++)
    has_deriv = has_deriv || (nc_xcor_kernel_component_get_bessel_deriv (g_ptr_array_index (comp_list, i)) > 0);

  {
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
      .n_terms                 = has_deriv ? 2 : 1,
      .k_vec                   = ncm_vector_new (n_l),
      .kf_vec                  = ncm_vector_new (n_l),
      .k_vec_p1                = has_deriv ? ncm_vector_new (n_l) : NULL,
      .kf_vec_p1               = has_deriv ? ncm_vector_new (n_l) : NULL,
      .y_terms                 = has_deriv ? ncm_vector_new (2 * n_l) : NULL
    };
    gdouble u_min_union = G_MAXDOUBLE;
    gdouble u_max_union = 0.0;

    for (i = 0; i < n_comp; i++)
    {
      ComponentState *state = &comp_states.states[i];
      const gboolean deriv  = (nc_xcor_kernel_component_get_bessel_deriv (g_ptr_array_index (comp_list, i)) > 0);
      gdouble u_band_min, u_band_max;

      _nc_xcor_kernel_component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

      comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard / (deriv ? nu_max + 1.0 : nu_max));
      comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard / nu_min);

      u_band_min = 1.0 / state->chi_max;
      u_band_max = 1.0 / state->chi_min;

      for (j = 0; j < n_l; j++)
      {
        const gdouble nu_j = lmin + j + 0.5;

        state->k_min_limber_ell[j] = GSL_MAX (u_band_min, state->k_min_hard / nu_j);
        state->k_max_limber_ell[j] = GSL_MIN (u_band_max, state->k_max_hard / nu_j);

        if (deriv)
        {
          state->k_min_limber_ell_p1[j] = GSL_MAX (u_band_min, state->k_min_hard / (nu_j + 1.0));
          state->k_max_limber_ell_p1[j] = GSL_MIN (u_band_max, state->k_max_hard / (nu_j + 1.0));
        }
      }

      u_min_union = GSL_MIN (u_min_union, u_band_min);
      u_max_union = GSL_MAX (u_max_union, u_band_max);
    }

    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, u_min_union);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, u_max_union);

    g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

    return comp_states;
  }
}

/*
 * Evaluates the Limber block at x: the point chi = 1 / x is shared by every multipole
 * and by both terms, so each component's window is one evaluation and only the wave
 * number factors are per multipole. Rows [0, n_l) are the j_l terms in u = k / nu, at
 * k = nu x; with two terms, rows [n_l, 2 n_l) are the j_{l+1} terms of the Bessel
 * derivative components in x = k / (nu + 1), at k = (nu + 1) x. Under Limber
 *
 *   int K(chi, k) j_l'(k chi) dchi = (l / nu) P_l(k) - P_{l+1}(k),
 *   int K(chi, k) j_l''(k chi) dchi = -(2 l + 1/4) / nu^2 P_l(k) + 2 / (nu + 1) P_{l+1}(k),
 *
 * with P_m(k) = sqrt(pi / (2 (m + 1/2))) K(chi_m, k) / k the j_m peak term, K the radial
 * kernel, the window factor times the wave number factor, at chi_m = (m + 1/2) / k. Outside its band a term is zero.
 * In panel mode the choice is made at the panel midpoint, so a panel ending on a band
 * edge takes the limit from its own side.
 */
static void
_component_states_compute_limber_u (const gdouble u, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  NcDistance *dist             = nc_xcor_kernel_peek_dist (comp_states->xclk);
  NcHICosmo *cosmo             = comp_states->states[0].params.cosmo;
  const gdouble u_side         = comp_states->panel_mode ? comp_states->interval_mid : u;
  const guint n_l              = comp_states->n_l;
  NcXcorKinetic xck;
  guint ci, i;

  gdouble W[MAX_ELL_BLOCK]    = { 0.0 };
  gdouble W_p1[MAX_ELL_BLOCK] = { 0.0 };

  xck.chi_z = 1.0 / u;
  xck.z     = nc_distance_inv_comoving (dist, cosmo, xck.chi_z);
  xck.E_z   = nc_hicosmo_E (cosmo, xck.z);

  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    ComponentState *state = &comp_states->states[ci];
    const guint deriv     = nc_xcor_kernel_component_get_bessel_deriv (state->comp);
    gboolean active[MAX_ELL_BLOCK];
    gboolean active_p1[MAX_ELL_BLOCK];
    gint first_active    = -1;
    gint first_active_p1 = -1;
    gdouble window;

    /* Each multipole is nonzero within its own range, which differs between multipoles
     * where the component's hard k limits bind */
    for (i = 0; i < n_l; i++)
    {
      active[i] = (u_side >= state->k_min_limber_ell[i]) && (u_side <= state->k_max_limber_ell[i]);

      if (active[i] && (first_active < 0))
        first_active = i;

      active_p1[i] = (deriv > 0) && (u_side >= state->k_min_limber_ell_p1[i]) && (u_side <= state->k_max_limber_ell_p1[i]);

      if (active_p1[i] && (first_active_p1 < 0))
        first_active_p1 = i;
    }

    if ((first_active < 0) && (first_active_p1 < 0))
      continue;

    window = nc_xcor_kernel_component_eval_window (state->comp, cosmo, &xck);

    if (first_active >= 0)
    {
      /* An inactive multipole takes the k of an active one, inside the component's k
       * range; its value is not used */
      for (i = 0; i < n_l; i++)
        ncm_vector_set (comp_states->k_vec, i, (comp_states->lmin + (active[i] ? i : (guint) first_active) + 0.5) * u);

      nc_xcor_kernel_component_eval_kfactor_vec (state->comp, cosmo, &xck, comp_states->k_vec, comp_states->kf_vec);

      for (i = 0; i < n_l; i++)
      {
        const gint l         = comp_states->lmin + i;
        const gdouble nu     = l + 0.5;
        const gdouble k      = nu * u;
        const gdouble peak_l = sqrt (M_PI / (2.0 * nu)) * window * ncm_vector_get (comp_states->kf_vec, i);
        const gdouble coef   = (deriv == 0) ? 1.0 : ((deriv == 1) ? l / nu : -(2.0 * l + 0.25) / (nu * nu));

        if (active[i])
          W[i] += nc_xcor_kernel_component_eval_prefactor (state->comp, cosmo, k, l) * coef * peak_l / k;
      }
    }

    if (first_active_p1 >= 0)
    {
      for (i = 0; i < n_l; i++)
        ncm_vector_set (comp_states->k_vec_p1, i, (comp_states->lmin + (active_p1[i] ? i : (guint) first_active_p1) + 1.5) * u);

      nc_xcor_kernel_component_eval_kfactor_vec (state->comp, cosmo, &xck, comp_states->k_vec_p1, comp_states->kf_vec_p1);

      for (i = 0; i < n_l; i++)
      {
        const gint l           = comp_states->lmin + i;
        const gdouble nup      = l + 1.5;
        const gdouble k        = nup * u;
        const gdouble peak_lp1 = sqrt (M_PI / (2.0 * nup)) * window * ncm_vector_get (comp_states->kf_vec_p1, i);
        const gdouble coef     = (deriv == 1) ? -1.0 : 2.0 / nup;

        if (active_p1[i])
          W_p1[i] += nc_xcor_kernel_component_eval_prefactor (state->comp, cosmo, k, l) * coef * peak_lp1 / k;
      }
    }
  }

  for (i = 0; i < n_l; i++)
    ncm_vector_set (y, i, W[i]);

  if (comp_states->n_terms == 2)
    for (i = 0; i < n_l; i++)
      ncm_vector_set (y, n_l + i, W_p1[i]);
}

/*
 * The magnitude of each multipole at x for the domain search, |j_l term| + |j_{l+1}
 * term| with two terms: the domain ends where both are negligible.
 */
static void
_component_states_expand_limber_u (const gdouble u, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  const guint n_l              = comp_states->n_l;
  guint i;

  if (comp_states->n_terms == 1)
  {
    _component_states_compute_limber_u (u, y, user_data);

    return;
  }

  _component_states_compute_limber_u (u, comp_states->y_terms, user_data);

  for (i = 0; i < n_l; i++)
    ncm_vector_set (y, i, fabs (ncm_vector_get (comp_states->y_terms, i)) + fabs (ncm_vector_get (comp_states->y_terms, n_l + i)));
}

NcXcorKernelIntegrand *
_nc_xcor_kernel_build_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, gdouble reltol, gdouble peak_epsilon)
{
  const guint n_l      = lmax - lmin + 1;
  GPtrArray *comp_list = _nc_xcor_kernel_validate_component_list (xclk, n_l);

  if (comp_list == NULL)
    return NULL;

  /* Under Limber a multipole's window is supported on the band
   * [nu / chi_max, nu / chi_min] of each component and is zero outside it. The
   * closure is built in u = k / nu, where the band is the same for every
   * multipole of the block: the band edges are its only cuts, and the point
   * chi = 1 / u is shared by the block (_component_states_compute_limber_u()). The
   * j_{l+1} terms of Bessel derivative components are rows of their own in
   * x = k / (nu + 1), on the same band and the same point. */
  {
    ComponentStates comp_states      = _component_states_init_limber_u (xclk, lmin, n_l, comp_list, cosmo);
    NcXcorKernelIntegrand *integrand = NULL;

    g_ptr_array_unref (comp_list);

    integrand = _nc_xcor_kernel_build_cheb_integrand (xclk, cosmo, lmin, lmax,
                                                      &comp_states,
                                                      _component_states_expand_limber_u,
                                                      _component_states_compute_limber_u,
                                                      reltol, peak_epsilon);

    _nc_xcor_kernel_component_states_clear (&comp_states);

    return integrand;
  }
}


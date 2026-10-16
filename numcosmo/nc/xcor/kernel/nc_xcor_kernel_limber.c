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
 * The coefficient of j_{l + t} in the Bessel weight of order deriv, zero when the weight has
 * no such term or l + t < 0. The derivatives of j_l are exact combinations of neighbouring
 * orders, DLMF 10.51.1:
 *
 *   j_l' = (l j_{l-1} - (l + 1) j_{l+1}) / (2 l + 1),
 *   j_l'' = a j_{l-2} - b j_l + c j_{l+2},
 *
 * with a = l (l - 1) / ((2 l - 1) (2 l + 1)), b = (2 l^2 + 2 l - 1) / ((2 l - 1) (2 l + 3))
 * and c = (l + 1) (l + 2) / ((2 l + 1) (2 l + 3)). At l = 0 and 1 the coefficients of the
 * negative orders vanish.
 */
static gdouble
_limber_order_coef (const guint deriv, const gint l, const gint t)
{
  const gdouble dl = l;

  if (l + t < 0)
    return 0.0;

  switch (deriv)
  {
    case 0:
      return (t == 0) ? 1.0 : 0.0;

    case 1:

      switch (t)
      {
        case -1:
          return dl / (2.0 * dl + 1.0);

        case 1:
          return -(dl + 1.0) / (2.0 * dl + 1.0);

        default:
          return 0.0;
      }

    case 2:

      switch (t)
      {
        case -2:
          return dl * (dl - 1.0) / ((2.0 * dl - 1.0) * (2.0 * dl + 1.0));

        case 0:
          return -(2.0 * dl * dl + 2.0 * dl - 1.0) / ((2.0 * dl - 1.0) * (2.0 * dl + 3.0));

        case 2:
          return (dl + 1.0) * (dl + 2.0) / ((2.0 * dl + 1.0) * (2.0 * dl + 3.0));

        default:
          return 0.0;
      }

    default:
      g_error ("_limber_order_coef: Bessel derivative order %u is not supported under Limber.", deriv);

      return 0.0;
  }
}

/*
 * Initializes the component states of a Limber block. The Bessel weight of a component is
 * a combination of the orders l + t (_limber_order_coef()), and each order offset t used by
 * a component of the block is a closure term, sampled in x = k / (nu + t). Its point is
 * chi = 1 / x for every multipole and every term, and its band is
 * [1 / chi_max, 1 / chi_min] in x for every multipole, but the component's hard k limits
 * are [k_min / (nu + t), k_max / (nu + t)] in x, one range per multipole and term. The
 * band arrays hold the band within that range, where the multipole's term is nonzero, and
 * are empty where the term does not apply. The block's hard limits are the intersection of
 * the components' hard limits over their terms, narrowed to the union of the bands.
 */
static ComponentStates
_component_states_init_limber_u (NcXcorKernel *xclk, gint lmin, guint n_l,
                                 GPtrArray *comp_list, NcHICosmo *cosmo)
{
  NcXcorKernelPrivate *self = _nc_xcor_kernel_get_private (xclk);
  const guint n_comp        = comp_list->len;
  guint used                = 0; /* bit t + 2 for each order offset t used in the block */
  guint n_terms             = 0;
  gint offsets[MAX_LIMBER_TERMS];
  gint t;
  guint i, j;

  for (i = 0; i < n_comp; i++)
  {
    const guint deriv = nc_xcor_kernel_component_get_bessel_deriv (g_ptr_array_index (comp_list, i));

    for (t = -2; t <= 2; t++)
      for (j = 0; j < n_l; j++)
        if (_limber_order_coef (deriv, lmin + j, t) != 0.0)
          used |= 1u << (t + 2);
  }

  for (t = -2; t <= 2; t++)
    if (used & (1u << (t + 2)))
      offsets[n_terms++] = t;

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
      .n_terms                 = n_terms,
      .k_vec                   = ncm_vector_new (n_l),
      .kf_vec                  = ncm_vector_new (n_l),
      .y_terms                 = (n_terms > 1) ? ncm_vector_new (n_terms * n_l) : NULL
    };
    gdouble u_min_union = G_MAXDOUBLE;
    gdouble u_max_union = 0.0;
    guint tau;

    for (tau = 0; tau < n_terms; tau++)
      comp_states.offsets[tau] = offsets[tau];

    for (i = 0; i < n_comp; i++)
    {
      ComponentState *state = &comp_states.states[i];
      const guint deriv     = nc_xcor_kernel_component_get_bessel_deriv (g_ptr_array_index (comp_list, i));
      gdouble nu_lo         = G_MAXDOUBLE;
      gdouble nu_hi         = 0.0;

      _nc_xcor_kernel_component_state_init (state, g_ptr_array_index (comp_list, i), i, cosmo, n_l);

      for (tau = 0; tau < n_terms; tau++)
      {
        for (j = 0; j < n_l; j++)
        {
          const gint l       = lmin + j;
          const gdouble nu_m = l + offsets[tau] + 0.5;

          if (_limber_order_coef (deriv, l, offsets[tau]) != 0.0)
          {
            state->k_min_limber_ell[tau][j] = GSL_MAX (1.0 / state->chi_max, state->k_min_hard / nu_m);
            state->k_max_limber_ell[tau][j] = GSL_MIN (1.0 / state->chi_min, state->k_max_hard / nu_m);

            nu_lo = GSL_MIN (nu_lo, nu_m);
            nu_hi = GSL_MAX (nu_hi, nu_m);
          }
          else
          {
            /* An empty band: the term is off everywhere */
            state->k_min_limber_ell[tau][j] = G_MAXDOUBLE;
            state->k_max_limber_ell[tau][j] = 0.0;
          }
        }
      }

      g_assert_cmpfloat (nu_lo, <=, nu_hi);

      comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, state->k_min_hard / nu_hi);
      comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, state->k_max_hard / nu_lo);

      u_min_union = GSL_MIN (u_min_union, 1.0 / state->chi_max);
      u_max_union = GSL_MAX (u_max_union, 1.0 / state->chi_min);
    }

    comp_states.k_min_hard = GSL_MAX (comp_states.k_min_hard, u_min_union);
    comp_states.k_max_hard = GSL_MIN (comp_states.k_max_hard, u_max_union);

    g_assert_cmpfloat (comp_states.k_min_hard, <, comp_states.k_max_hard);

    return comp_states;
  }
}

/*
 * Evaluates the Limber block at x: the point chi = 1 / x is shared by every multipole and
 * by every term, so each component's window is one evaluation and only the wave number
 * factors are per multipole and term. Row tau n_l + i is the j_{l + t} term of multipole i,
 * t = offsets[tau], at k = (nu + t) x:
 *
 *   c_{l,t} p_l(k) L_{l+t}(k),  L_m(k) = sqrt(pi / (2 (m + 1/2))) K(chi_m, k) / k,
 *
 * with c_{l,t} the coefficient of j_{l+t} in the component's Bessel weight
 * (_limber_order_coef()), p_l the prefactor, K the radial kernel, the window factor times
 * the wave number factor, at chi_m = (m + 1/2) / k, and L_m the Limber form of
 * int K j_m(k chi) dchi. Outside its band a term is zero. In panel mode the choice is made
 * at the panel midpoint, so a panel ending on a band edge takes the limit from its own side.
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
  guint ci, i, tau;

  gdouble W[MAX_LIMBER_TERMS][MAX_ELL_BLOCK] = {
    {
      0.0
    }
  };

  xck.chi_z = 1.0 / u;
  xck.z     = nc_distance_inv_comoving (dist, cosmo, xck.chi_z);
  xck.E_z   = nc_hicosmo_E (cosmo, xck.z);

  for (ci = 0; ci < comp_states->n_comp; ci++)
  {
    ComponentState *state = &comp_states->states[ci];
    const guint deriv     = nc_xcor_kernel_component_get_bessel_deriv (state->comp);
    gboolean have_window  = FALSE;
    gdouble window        = 0.0;

    for (tau = 0; tau < comp_states->n_terms; tau++)
    {
      const gint t = comp_states->offsets[tau];
      gboolean active[MAX_ELL_BLOCK];
      gint first_active = -1;

      /* Each multipole is nonzero within its own range, which differs between multipoles
       * where the component's hard k limits bind */
      for (i = 0; i < n_l; i++)
      {
        active[i] = (u_side >= state->k_min_limber_ell[tau][i]) && (u_side <= state->k_max_limber_ell[tau][i]);

        if (active[i] && (first_active < 0))
          first_active = i;
      }

      if (first_active < 0)
        continue;

      if (!have_window)
      {
        window      = nc_xcor_kernel_component_eval_window (state->comp, cosmo, &xck);
        have_window = TRUE;
      }

      /* An inactive multipole takes the k of an active one, inside the component's k
       * range; its value is not used */
      for (i = 0; i < n_l; i++)
        ncm_vector_set (comp_states->k_vec, i, ((gint) (comp_states->lmin + (active[i] ? i : (guint) first_active)) + t + 0.5) * u);

      nc_xcor_kernel_component_eval_kfactor_vec (state->comp, cosmo, &xck, comp_states->k_vec, comp_states->kf_vec);

      for (i = 0; i < n_l; i++)
      {
        if (active[i])
        {
          const gint l       = comp_states->lmin + i;
          const gdouble nu_m = l + t + 0.5;
          const gdouble k    = nu_m * u;
          const gdouble peak = sqrt (M_PI / (2.0 * nu_m)) * window * ncm_vector_get (comp_states->kf_vec, i);
          const gdouble coef = _limber_order_coef (deriv, l, t);

          W[tau][i] += nc_xcor_kernel_component_eval_prefactor (state->comp, cosmo, k, l) * coef * peak / k;
        }
      }
    }
  }

  for (tau = 0; tau < comp_states->n_terms; tau++)
    for (i = 0; i < n_l; i++)
      ncm_vector_set (y, tau * n_l + i, W[tau][i]);
}

/*
 * The magnitude of each multipole at x for the domain search, the sum of the magnitudes
 * of its terms: the domain ends where all are negligible.
 */
static void
_component_states_expand_limber_u (const gdouble u, NcmVector *y, gpointer user_data)
{
  ComponentStates *comp_states = (ComponentStates *) user_data;
  const guint n_l              = comp_states->n_l;
  guint i, tau;

  if (comp_states->n_terms == 1)
  {
    _component_states_compute_limber_u (u, y, user_data);

    return;
  }

  _component_states_compute_limber_u (u, comp_states->y_terms, user_data);

  for (i = 0; i < n_l; i++)
  {
    gdouble sum = 0.0;

    for (tau = 0; tau < comp_states->n_terms; tau++)
      sum += fabs (ncm_vector_get (comp_states->y_terms, tau * n_l + i));

    ncm_vector_set (y, i, sum);
  }
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
   * Bessel derivative weights are combinations of the orders l + t, each a closure term
   * of its own in x = k / (nu + t), on the same band and the same point. */
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


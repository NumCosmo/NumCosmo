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
 * the panelling and acceptance of the Chebyshev closure, and its builder.
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
  ncm_vector_clear (&comp_states->k_vec_p1);
  ncm_vector_clear (&comp_states->kf_vec_p1);
  ncm_vector_clear (&comp_states->y_terms);
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
   * stack arrays (ComponentState::k_min_limber_ell, kernel_out[][]), and asserts
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
 * Checks that the closure tolerances are not below the relative tolerance of
 * the Levin integrator, and aborts with g_error() when they are: a closure
 * cannot be more accurate than the values of W_l(k) it interpolates. The
 * looser of reltol and peak-epsilon is compared, since refinement stops at
 * the looser of the two. Other integrators are not checked.
 */
void
_nc_xcor_kernel_check_integrator_tolerance (NcXcorKernel *xclk, NcmSBesselIntegrator *sbi, gdouble reltol, gdouble peak_epsilon)
{
  const gdouble fit_tol = GSL_MAX (reltol, peak_epsilon);
  gdouble integrator_reltol;

  if ((sbi == NULL) || !NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbi))
    return;

  {
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

    integrator_reltol = ncm_sbessel_integrator_levin_get_reltol (sbilv);
  }

  if (fit_tol < integrator_reltol)
    g_error ("_nc_xcor_kernel_check_integrator_tolerance: kernel %s builds its closure of W_l(k) "
             "to %.17g but the integrator computes W_l(k) only to %.17g (its reltol). A "
             "closure cannot resolve what the values do not carry. "
             "Loosen NcXcorKernel:reltol and NcXcorKernel:peak-epsilon to at least %.17g, "
             "or construct the integrator with a tighter reltol.",
             G_OBJECT_TYPE_NAME (xclk), fit_tol, integrator_reltol, integrator_reltol);
}

/*
 * What the panel builder hands back to the compute functions: they take k first, as
 * NcmFunctionSampleSetFunc does, and decide which components are on at a point of the
 * continuous interval the panel lies in, so a panel ending on a cut takes the limit from
 * its own side.
 */
typedef struct _ChebCompute
{
  void (*compute_func) (const gdouble, NcmVector *, gpointer);

  ComponentStates *comp_states;
  GArray *interval_mid; /* gdouble, the midpoint of each continuous interval */
} ChebCompute;

static void
_cheb_compute_call (gpointer user_data, guint i, gdouble k, NcmVector *fval)
{
  ChebCompute *compute = (ChebCompute *) user_data;

  compute->comp_states->interval_mid = g_array_index (compute->interval_mid, gdouble, i);
  compute->compute_func (k, fval, compute->comp_states);
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
 * the truncation boundary of every component whose boundary was found, and under
 * Limber the band edges of every component for every multipole of the block, of its
 * j_{l+1} term too when it has a Bessel derivative (see
 * _component_states_compute_limber_u()). Each cut becomes a panel edge, with the
 * component or multipole off on the outer side, so that every panel is a
 * smooth function; a polynomial interpolating across a jump does not converge
 * below the size of the jump. The boundaries are fixed before this runs: only
 * the domain expansion moves them, and it has finished.
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

        if (nc_xcor_kernel_component_get_bessel_deriv (state->comp) > 0)
        {
          if ((state->k_min_limber_ell_p1[j] > k_min) && (state->k_min_limber_ell_p1[j] < k_max))
            g_array_append_val (cuts, state->k_min_limber_ell_p1[j]);

          if ((state->k_max_limber_ell_p1[j] > k_min) && (state->k_max_limber_ell_p1[j] < k_max))
            g_array_append_val (cuts, state->k_max_limber_ell_p1[j]);
        }
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

  /* A window that does not vanish at its near edge chi_min adds a feature of its own
   * where that edge reaches the turning point of j_l, k = nu / chi_min, falling as a power
   * of k beyond it. Between the main peak and that feature the block can stay below the
   * truncation threshold for several points in a row, where the domain expansion stops,
   * so the feature is seeded, at the block's largest nu. Not under Limber, where a band
   * edge is a step the panels cut at. */
  if (!comp_states->is_limber)
  {
    const gdouble nu_max = comp_states->lmin + comp_states->n_l - 1 + 0.5;

    for (i = 0; i < comp_states->n_comp; i++)
    {
      const gdouble chi_min = comp_states->states[i].chi_min;

      if (chi_min > 0.0)
      {
        const gdouble k_edge = nu_max / chi_min;

        if ((k_edge > comp_states->k_min_hard) && (k_edge < comp_states->k_max_hard) &&
            _is_new_k (k_edge, (gdouble *) k_seeds->data, k_seeds->len))
          g_array_append_val (k_seeds, k_edge);
      }
    }
  }
}

/*
 * Builds the closure as a Chebyshev series.
 *
 * The domain is found from the seeds by ncm_function_sample_set_expand_domain(),
 * which stops where W is negligible. On it the whole ell block is expanded on one
 * Chebyshev-Lobatto grid per panel, doubling the order until every multipole's
 * coefficients converge.
 *
 * The panels start from the component boundaries the expansion found, since
 * the block W_l(k) jumps there (see _component_states_collect_cuts), and from
 * the grid of #NcmSpectralPanelsBuilder:min-panels-per-efold, so that a feature narrower than
 * the domain is not missed by a first grid spread over all of it.
 */
NcXcorKernelIntegrand *
_nc_xcor_kernel_build_cheb_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax,
                                      ComponentStates *comp_states,
                                      void (*expand_func) (const gdouble, NcmVector *, gpointer),
                                      void (*compute_func) (const gdouble, NcmVector *, gpointer),
                                      const gdouble reltol, const gdouble peak_epsilon)
{
  NcXcorKernelPrivate *self        = _nc_xcor_kernel_get_private (xclk);
  NcXcorKernelIntegrand *integrand = _nc_xcor_kernel_integrand_new (lmax - lmin + 1, comp_states->n_terms);
  const guint n_l                  = lmax - lmin + 1;

  {
    NcmFunctionSampleSet *fss = ncm_function_sample_set_new (n_l);
    GArray *k_seeds           = g_array_new (FALSE, FALSE, sizeof (gdouble));
    ChebCompute compute       = {
      compute_func, comp_states, g_array_new (FALSE, FALSE, sizeof (gdouble))
    };
    gdouble abstol;
    guint i;

    /* The seeds are local to the call, not kernel state: one kernel may be
     * evaluated concurrently for different ell blocks, each with its own
     * integrator. */
    _component_states_compute_k_seeds (comp_states, k_seeds);

    for (i = 0; i < k_seeds->len; i++)
    {
      const gdouble k_seed = g_array_index (k_seeds, gdouble, i);

      ncm_function_sample_set_add_old_func (fss, k_seed, expand_func, comp_states);
    }

    g_array_unref (k_seeds);

    ncm_function_sample_set_expand_domain (
      fss,
      expand_func,
      comp_states->k_min_hard,
      comp_states->k_max_hard,
      self->expansion_factor,
      comp_states->epsilon,
      self->max_border_expansions,
      comp_states->adaptive_boundary_tries,
      comp_states
    );

    integrand->k_min = ncm_function_sample_set_get_x_min (fss);
    integrand->k_max = ncm_function_sample_set_get_x_max (fss);

    /* A floor scaled to the smallest of the block's peaks, so a sub-dominant
     * multipole is not held to a tolerance relative to its neighbours. */
    abstol = ncm_function_sample_set_get_absmaxF_min (fss) * peak_epsilon;

    /* Expand on the continuous intervals between the cuts, from the builder's grid,
     * bisecting where a panel does not converge. In panel mode the compute function
     * decides at a point of the panel's interval which components and multipoles are
     * on, so every discontinuity has to be a cut, which
     * _component_states_collect_cuts() guarantees. */
    {
      GArray *cuts = g_array_new (FALSE, FALSE, sizeof (gdouble));
      NcmVector *breakpoints;

      _component_states_collect_cuts (comp_states, integrand->k_min, integrand->k_max, cuts);

      for (i = 0; i <= cuts->len; i++)
      {
        const gdouble lo  = (i == 0) ? integrand->k_min : g_array_index (cuts, gdouble, i - 1);
        const gdouble hi  = (i == cuts->len) ? integrand->k_max : g_array_index (cuts, gdouble, i);
        const gdouble mid = 0.5 * (lo + hi);

        g_array_append_val (compute.interval_mid, mid);
      }

      breakpoints = (cuts->len > 0) ? ncm_vector_new_array (cuts) : NULL;

      comp_states->panel_mode = TRUE;
      ncm_spectral_panels_builder_build (self->builder, _cheb_compute_call, &compute,
                                         integrand->k_min, integrand->k_max, breakpoints, NULL, NULL,
                                         reltol, abstol, integrand->sp);
      comp_states->panel_mode = FALSE;

      ncm_vector_clear (&breakpoints);
      g_array_unref (cuts);
    }

    ncm_function_sample_set_clear (&fss);
    g_array_unref (compute.interval_mid);

    integrand->scaled = comp_states->in_u;

    /* The sampled range is in x; each panel row's k range is its scale times
     * the part of it inside the band, and the block's k range is their union. Row
     * i + n_l is the j_{l+1} term of multipole i, of the Bessel derivative components,
     * with scale nu + 1. */
    for (i = 0; i < comp_states->n_terms * n_l; i++)
    {
      NcXcorKernelIntegrandRow *row = &g_array_index (integrand->rows, NcXcorKernelIntegrandRow, i);
      const gboolean p1             = (i >= n_l);
      const guint i_l               = p1 ? i - n_l : i;
      const gdouble x_min           = integrand->k_min;
      const gdouble x_max           = integrand->k_max;
      gdouble k_min_i               = x_min;
      gdouble k_max_i               = x_max;

      if (comp_states->is_limber)
      {
        gdouble band_min = G_MAXDOUBLE;
        gdouble band_max = 0.0;
        guint ci;

        for (ci = 0; ci < comp_states->n_comp; ci++)
        {
          const ComponentState *state = &comp_states->states[ci];

          if (!p1)
          {
            band_min = GSL_MIN (band_min, state->k_min_limber_ell[i_l]);
            band_max = GSL_MAX (band_max, state->k_max_limber_ell[i_l]);
          }
          else if (nc_xcor_kernel_component_get_bessel_deriv (state->comp) > 0)
          {
            band_min = GSL_MIN (band_min, state->k_min_limber_ell_p1[i_l]);
            band_max = GSL_MAX (band_max, state->k_max_limber_ell_p1[i_l]);
          }
        }

        k_min_i = GSL_MAX (k_min_i, band_min);
        k_max_i = GSL_MIN (k_max_i, band_max);

        /* A band disjoint from the closure's domain leaves nothing to integrate;
         * report the empty range as the domain's lower edge rather than an
         * inverted one. */
        if (k_min_i >= k_max_i)
          k_min_i = k_max_i = x_min;
      }

      const gdouble scale = comp_states->in_u ? (lmin + i_l + (p1 ? 1.5 : 0.5)) : 1.0;

      ncm_vector_set (integrand->scales, i, scale);
      row->x_min = k_min_i;
      row->x_max = k_max_i;
      row->k_min = k_min_i * scale;
      row->k_max = k_max_i * scale;
    }

    if (integrand->scaled)
    {
      integrand->k_min *= ncm_vector_get (integrand->scales, 0);
      integrand->k_max *= ncm_vector_get (integrand->scales, comp_states->n_terms * n_l - 1);
    }

    return integrand;
  }
}


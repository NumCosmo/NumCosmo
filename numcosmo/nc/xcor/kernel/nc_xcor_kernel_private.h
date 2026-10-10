/***************************************************************************
 *            nc_xcor_kernel_private.h
 *
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
 * Private to the NcXcorKernel sources, nc_xcor_kernel.c and the files of this
 * directory: the private data of the object, the data layouts of the two
 * closures, the component states that build them, and the functions one source
 * file calls in another.
 */

#ifndef _NC_XCOR_KERNEL_PRIVATE_H_
#define _NC_XCOR_KERNEL_PRIVATE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/spline/ncm_spline_vec.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator.h>
#include <numcosmo/ncm/algebra/ncm_spectral.h>
#include <numcosmo/ncm/stats/ncm_function_sample_set.h>
#include <numcosmo/ncm/model/ncm_model_ctrl.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_component.h>

G_BEGIN_DECLS

typedef struct _NcXcorKernelPrivate
{
  /*< private >*/
  NcmModel parent_instance;
  NcDistance *dist;
  NcmPowspec *ps;
  NcmSBesselIntegrator *sbi;
  guint lmax;
  gint l_limber;
  gdouble adaptive_epsilon;
  guint adaptive_boundary_tries;
  gdouble reltol;
  gdouble peak_epsilon;
  guint max_border_expansions;
  guint max_iter;
  gdouble expansion_factor;
  guint panel_order_cap;
  gdouble panels_per_efold;
  guint panel_level_min;
  gboolean track_closure_error;
  gboolean tolerance_balance_warned;
  gboolean constructed;
  NcmModelCtrl *cosmo_ctrl;
  guint64 prepared_pkey;
  gboolean outdated;
} NcXcorKernelPrivate;

/* Spline closure of one ell block, built by both the Limber and the non-Limber paths. */
typedef struct _SplineIntegrandData
{
  NcHICosmo *cosmo;
  gdouble RH_Mpc;
  gint lmin;
  guint len;
  NcmSplineVec *spline_vec;
  NcmVector *eval_result;
  gdouble k_min;
  gdouble k_max;
  gdouble *k_min_comp;
  gdouble *k_max_comp;
} SplineIntegrandData;

/* One panel [a, b] of the Chebyshev closure, with N coefficients per multipole. */
typedef struct _ChebPanel
{
  gdouble a;
  gdouble b;
  NcmMatrix *coeffs; /* len x N, one row per multipole */
  guint N;
} ChebPanel;

/* Chebyshev closure of one ell block: contiguous panels of coefficients over the
 * domain [k_min, k_max], with the per-component ranges of SplineIntegrandData. */
typedef struct _ChebIntegrandData
{
  NcHICosmo *cosmo;
  gdouble RH_Mpc;
  gint lmin;
  guint len;
  GArray *panels;  /* ChebPanel, ascending and contiguous, in x = k / scale[i] */
  GArray *edges;   /* gdouble, panels->len + 1 entries, for the lookup */
  gdouble *scale;  /* per multipole; 1 when the closure is in k */
  gboolean scaled; /* whether any scale differs from 1 */
  gdouble k_min;   /* the k range of the whole block */
  gdouble k_max;
  gdouble *k_min_comp;
  gdouble *k_max_comp;
} ChebIntegrandData;

/*
 * Evaluates the Chebyshev series of component comp at s in [-1, 1] by the
 * Clenshaw recurrence, in O(N) operations.
 */
static inline gdouble
_nc_xcor_kernel_cheb_panel_eval_one (const ChebPanel *panel, guint comp, const gdouble s)
{
  const gdouble two_s = 2.0 * s;
  gdouble b_1         = 0.0;
  gdouble b_2         = 0.0;
  gint n;

  for (n = (gint) panel->N - 1; n >= 1; n--)
  {
    const gdouble b_0 = two_s * b_1 - b_2 + ncm_matrix_get (panel->coeffs, comp, n);

    b_2 = b_1;
    b_1 = b_0;
  }

  return s * b_1 - b_2 + ncm_matrix_get (panel->coeffs, comp, 0);
}

typedef struct _ComponentParams
{
  NcXcorKernelComponent *comp;
  NcHICosmo *cosmo;
} ComponentParams;

#define MAX_ELL_BLOCK NC_XCOR_KERNEL_MAX_ELL_BLOCK
#define MAX_COMP_BLOCK 6

typedef struct _ComponentState
{
  NcXcorKernelComponent *comp;
  guint comp_idx;
  gdouble chi_min;
  gdouble chi_max;
  gdouble k_min_hard;
  gdouble k_max_hard;
  gdouble last_k_left;
  gdouble last_k_right;
  gdouble last_values_left[MAX_ELL_BLOCK];
  gdouble last_values_right[MAX_ELL_BLOCK];
  guint left_boundary_found;
  guint right_boundary_found;
  guint group;                             /* Components truncated together, see _component_states_init_non_limber */
  gdouble k_min_limber_ell[MAX_ELL_BLOCK]; /* Per-ell minimum k for Limber */
  gdouble k_max_limber_ell[MAX_ELL_BLOCK]; /* Per-ell maximum k for Limber */
  ComponentParams params;
} ComponentState;

typedef struct _ComponentStates
{
  ComponentState states[MAX_COMP_BLOCK];
  NcXcorKernel *xclk;
  NcmSBesselIntegrator *sbi; /* Integrator for this call; not owned */
  gdouble k_min_hard;
  gdouble k_max_hard;
  gdouble l2_norm;
  const guint n_comp;
  const guint lmin;
  const guint n_l;
  const gdouble epsilon;
  const guint adaptive_boundary_tries;
  const gboolean is_limber; /* k_min_limber_ell/k_max_limber_ell are set only then */
  const gboolean in_u;      /* the sampled variable is u = k / nu: k_min_hard, k_max_hard,
                             * the band arrays and every sample are in u */

  /* Chebyshev path only. While a panel is being expanded a component is on or
   * off for the whole panel, decided at panel_mid against its boundaries. */
  gboolean panel_mode;
  gdouble panel_mid;

  /* Limber in u only: the block's wave numbers and wave number factors at one node */
  NcmVector *k_vec;
  NcmVector *kf_vec;
} ComponentStates;

#define DECAY_RATE 1.0e10

NcXcorKernelPrivate *_nc_xcor_kernel_get_private (NcXcorKernel *xclk);

NcmSpectral **_nc_xcor_kernel_spectral_get (void);
NcXcorKernelIntegrand *_nc_xcor_kernel_cheb_integrand_new (guint n_l, ChebIntegrandData *cid);
NcXcorKernelIntegrand *_nc_xcor_kernel_spline_integrand_new (guint n_l, SplineIntegrandData *sid);

void _nc_xcor_kernel_component_state_init (ComponentState *state, NcXcorKernelComponent *comp, guint comp_idx, NcHICosmo *cosmo, guint n_l);
void _nc_xcor_kernel_component_states_clear (ComponentStates *comp_states);
GPtrArray *_nc_xcor_kernel_validate_component_list (NcXcorKernel *xclk, guint n_l);
void _nc_xcor_kernel_check_integrator_tolerance (NcXcorKernel *xclk, NcmSBesselIntegrator *sbi);
NcXcorKernelIntegrand *_nc_xcor_kernel_build_cheb_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, ComponentStates *comp_states, void (*compute_func) (const gdouble, NcmVector *, gpointer), gdouble reltol, gdouble peak_epsilon);
NcXcorKernelIntegrand *_nc_xcor_kernel_build_spline_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, ComponentStates *comp_states, void (*compute_func) (const gdouble, NcmVector *, gpointer), gdouble reltol, gdouble peak_epsilon);

NcXcorKernelIntegrand *_nc_xcor_kernel_build_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcXcorKernelClosure closure_type);

NcXcorKernelIntegrand *_nc_xcor_kernel_build_non_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi, NcXcorKernelClosure closure_type);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_PRIVATE_H_ */


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
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator.h>
#include <numcosmo/ncm/algebra/ncm_spectral.h>
#include <numcosmo/ncm/algebra/ncm_spectral_panels.h>
#include <numcosmo/ncm/algebra/ncm_spectral_panels_builder.h>
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
  gdouble expansion_factor;
  NcmSpectralPanelsBuilder *builder; /* the panels of the Chebyshev closures */
  gboolean constructed;
  NcmModelCtrl *cosmo_ctrl;
  guint64 prepared_pkey;
  gboolean outdated;
} NcXcorKernelPrivate;

/* One panel row of a closure: its support in k and in its variable x = k / s, s its scale. */
typedef struct _NcXcorKernelIntegrandRow
{
  gdouble k_min;
  gdouble k_max;
  gdouble x_min; /* the ends of the support in x, panel edges: k_min and k_max over the scale */
  gdouble x_max;
} NcXcorKernelIntegrandRow;

/* The Chebyshev closure of one ell block: panels of coefficients over the block's
 * domain, n_terms groups of len rows; multipole i is the sum of the panel rows
 * i + t len over the terms t. */
struct _NcXcorKernelIntegrand
{
  gint refcount;
  guint len;             /* the multipoles of the block */
  guint n_terms;         /* the groups of len panel rows summed into each multipole */
  NcmSpectralPanels *sp; /* n_terms len rows, row r in x = k / scale_r */
  GArray *rows;          /* NcXcorKernelIntegrandRow, one per panel row */
  NcmVector *scales;     /* per panel row, 1 when the closure is in k, nu = l + 1/2 when it is in u = k / nu */
  gboolean scaled;       /* whether any scale differs from 1 */
  gdouble k_min;         /* the k range of the whole block */
  gdouble k_max;
  GPtrArray *windows; /* NcmSpectralPanels per multipole, its window as one row in k, the sum of its terms; NULL with one term */
};

/* The fraction of the fit's absolute floor below which the sum of a multipole's terms drops
 * trailing coefficients: at most that fraction times their number is added to its window */
#define NC_XCOR_KERNEL_WINDOW_SUM_FRACTION (1.0e-2)

/* The factor of the first fit's tolerances for a block of several terms: the fit meets
 * tolerances relative to the windows when their smallest peak is at least this fraction
 * of the terms' smallest peak, and is redone otherwise */
#define NC_XCOR_KERNEL_WINDOW_FIT_FACTOR (0.5)


typedef struct _ComponentParams
{
  NcXcorKernelComponent *comp;
  NcHICosmo *cosmo;
} ComponentParams;

#define MAX_ELL_BLOCK NC_XCOR_KERNEL_MAX_ELL_BLOCK
#define MAX_COMP_BLOCK 6
#define MAX_LIMBER_TERMS 5 /* the order offsets -2, ..., 2 of the Bessel derivative weights */


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
  guint left_boundary_found;
  guint right_boundary_found;
  guint group;                                               /* Components truncated together, see _component_states_init_non_limber */
  gdouble k_min_limber_ell[MAX_LIMBER_TERMS][MAX_ELL_BLOCK]; /* Limber: per term and multipole, the band in the term's x, empty when the term does not apply */
  gdouble k_max_limber_ell[MAX_LIMBER_TERMS][MAX_ELL_BLOCK];
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
  const gboolean is_limber;       /* k_min_limber_ell/k_max_limber_ell are set only then */
  const gboolean in_u;            /* the sampled variable is u = k / nu: k_min_hard, k_max_hard,
                                   * the band arrays and every sample are in u */
  const guint n_terms;            /* Limber: the order offsets of the block, the rows
                                   * [tau n_l, (tau + 1) n_l) holding the j_{l + t_tau} terms
                                   * in x = k / (nu + t_tau); 1 otherwise */
  gint offsets[MAX_LIMBER_TERMS]; /* t_tau, ascending; {0} outside Limber */

  /* While a panel is being expanded a component is on or off for the whole
   * panel, decided at interval_mid, the midpoint of the continuous interval
   * holding the panel, against its boundaries. */
  gboolean panel_mode;
  gdouble interval_mid;

  /* Limber only: the block's wave numbers and wave number factors of one term at one node;
   * with several terms the rows of all of them at one node */
  NcmVector *k_vec;
  NcmVector *kf_vec;
  NcmVector *y_terms;
} ComponentStates;

NcXcorKernelPrivate *_nc_xcor_kernel_get_private (NcXcorKernel *xclk);

NcXcorKernelIntegrand *_nc_xcor_kernel_integrand_new (guint len, guint n_terms);
void _nc_xcor_kernel_integrand_sum_terms (NcXcorKernelIntegrand *integrand, NcmSpectral *spectral, gdouble abstol);
gdouble _nc_xcor_kernel_integrand_window_peak_min (NcXcorKernelIntegrand *integrand);

void _nc_xcor_kernel_component_state_init (ComponentState *state, NcXcorKernelComponent *comp, guint comp_idx, NcHICosmo *cosmo, guint n_l);
void _nc_xcor_kernel_component_states_clear (ComponentStates *comp_states);
GPtrArray *_nc_xcor_kernel_validate_component_list (NcXcorKernel *xclk, guint n_l);
void _nc_xcor_kernel_check_integrator_tolerance (NcXcorKernel *xclk, NcmSBesselIntegrator *sbi, gdouble reltol, gdouble peak_epsilon);
NcXcorKernelIntegrand *_nc_xcor_kernel_build_cheb_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, ComponentStates *comp_states, void (*expand_func) (const gdouble, NcmVector *, gpointer), void (*compute_func) (const gdouble, NcmVector *, gpointer), gdouble reltol, gdouble peak_epsilon);

NcXcorKernelIntegrand *_nc_xcor_kernel_build_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, gdouble reltol, gdouble peak_epsilon);

NcXcorKernelIntegrand *_nc_xcor_kernel_build_non_limber_integrand (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi, gdouble reltol, gdouble peak_epsilon);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_PRIVATE_H_ */


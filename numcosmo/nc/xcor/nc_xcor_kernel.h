/***************************************************************************
 *            nc_xcor_kernel.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
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

#ifndef _NC_XCOR_KERNEL_H_
#define _NC_XCOR_KERNEL_H_

#include <glib-object.h>
#include <glib.h>

#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/powspec/nc_growth_func.h>
#include <numcosmo/nc/powspec/nc_transfer_func.h>
#include <numcosmo/ncm/core/ncm_c.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/powspec/ncm_powspec.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/background/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL (nc_xcor_kernel_get_type ())
#define NC_TYPE_XCOR_KERNEL_INTEGRAND (nc_xcor_kernel_integrand_get_type ())

/**
 * NC_XCOR_KERNEL_MAX_ELL_BLOCK:
 *
 * Hard cap on the number of multipoles in a single get_eval_vectorized()
 * call: #NcXcorKernel's internal per-block state uses fixed-size stack
 * arrays sized by this constant (not just the Levin integrator's own,
 * larger ell_cache_max). Exceeding it is a fatal, non-catchable g_error.
 * Public so callers planning ℓ-block tilings (e.g. #NcXcorSolver) can
 * respect it without duplicating the number.
 */
#define NC_XCOR_KERNEL_MAX_ELL_BLOCK 64

/**
 * NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL:
 *
 * Smallest #NcXcorKernel:scaled-abstol worth asking for. The tolerance is a
 * fraction of the peak of $W_i(k)$, but the quantity integrated to form
 * $C_\ell$ is $k^2 W_i W_j$, so it enters *squared*: this floor is $10^{-12}$
 * on the integrand, already past what the outer $k$ integral carries. Below it
 * nc_xcor_kernel_set_scaled_abstol() warns.
 */
#define NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL (1.0e-6)

G_DECLARE_DERIVABLE_TYPE (NcXcorKernel, nc_xcor_kernel, NC, XCOR_KERNEL, NcmModel);

typedef struct _NcXcorKinetic NcXcorKinetic;
typedef struct _NcXcorKernelIntegrand NcXcorKernelIntegrand;

/**
 * NcXcorKernelIntegrandEval:
 * @data: user data
 * @k: wavenumber
 * @W: (array) (out): output array to fill with integrand values
 *
 * Function type for evaluating kernel integrands.
 */
typedef void (*NcXcorKernelIntegrandEval) (gpointer data, gdouble k, gdouble *W);

/**
 * NcXcorKernelIntegrandGetRange:
 * @data: user data
 * @k_min: (out): minimum wavenumber
 * @k_max: (out): maximum wavenumber
 *
 * Function type for getting the valid k range.
 */
typedef void (*NcXcorKernelIntegrandGetRange) (gpointer data, gdouble *k_min, gdouble *k_max);

/**
 * NcXcorKernelIntegrandGetRangeComp:
 * @data: user data
 * @i: component index
 * @k_min: (out): minimum wavenumber
 * @k_max: (out): maximum wavenumber
 *
 * Function type for getting the valid k range of a single component. A block
 * of multipoles shares one k-domain, but each multipole may be supported on
 * only part of it -- under the Limber approximation a multipole's window
 * vanishes outside $[\nu/\xi_\mathrm{max}, \nu/\xi_\mathrm{min}]$, and the
 * edge of that band is a step in the shared domain. See
 * nc_xcor_kernel_integrand_get_range_comp().
 */
typedef void (*NcXcorKernelIntegrandGetRangeComp) (gpointer data, guint i, gdouble *k_min, gdouble *k_max);

/**
 * NcXcorKernelIntegrandEvalComps:
 * @data: user data
 * @k: wavenumber
 * @offset: index of the first component to evaluate
 * @len: number of components to evaluate
 * @W: (array): full-length buffer, of which only
 *   [@offset, @offset + @len) need be filled
 *
 * Function type for evaluating a contiguous run of components, leaving the
 * rest of @W untouched. See nc_xcor_kernel_integrand_eval_comps().
 */
typedef void (*NcXcorKernelIntegrandEvalComps) (gpointer data, gdouble k, guint offset, guint len, gdouble *W);

/**
 * NcXcorKernelIntegrandGetKnots:
 * @data: user data
 *
 * Function type for getting the knots the integrand is represented on, when it
 * is spline-backed. See nc_xcor_kernel_integrand_peek_knots().
 *
 * Returns: (transfer none): the knot vector.
 */
typedef NcmVector *(*NcXcorKernelIntegrandGetKnots) (gpointer data);

/**
 * NcXcorKernelIntegrand:
 * @refcount: atomic reference count
 * @len: number of components in the integrand
 * @eval_func: function to evaluate the integrand at @k, filling @W[@len]
 * @get_range_func: function to get the valid k range for this integrand
 * @get_knots_func: function to get the integrand's knots, or %NULL when it is
 *   not spline-backed
 * @get_range_comp_func: function to get the valid k range of one component, or
 *   %NULL when every component covers the whole range
 * @eval_comps_func: function to evaluate a run of components, or %NULL when
 *   only the whole vector can be evaluated at once
 * @reltol: the relative half of the fit criterion this integrand was built to,
 *   or 0.0 when it is exact or unknown
 * @scaled_abstol: the floor of that criterion, as a fraction of the fitted
 *   function's own peak, or 0.0 when there was none
 * @data: user data passed to @eval_func, @get_range_func and @get_knots_func
 * @data_free: function to free @data, or %NULL if no cleanup needed
 *
 * A reference-counted closure for computing kernel integrands.
 *
 * **One integrand must be evaluated by one thread at a time.** A spline-backed
 * integrand keeps a scratch vector for the result of each evaluation, so
 * concurrent nc_xcor_kernel_integrand_eval() calls on the *same* integrand
 * would race on it.
 *
 * #NcXcorSolver satisfies this by construction rather than by convention: an
 * integrand is built for one (kernel, ell-block) pair, and the ell block is
 * the unit of parallelism -- one integrator per block, blocks distributed
 * across the OpenMP team, kernels shared and read-only throughout. No two
 * threads ever hold the same integrand. Anything that made kernel *pairs* the
 * unit of parallelism instead would share integrands across threads and would
 * need per-thread evaluation scratch.
 * The @eval_func function should fill @len values in the @W array
 * for the given wavenumber @k.
 */
struct _NcXcorKernelIntegrand
{
  /*< private >*/
  gint refcount;
  /*< public >*/
  guint len;
  NcXcorKernelIntegrandEval eval_func;
  NcXcorKernelIntegrandGetRange get_range_func;
  gpointer data;
  GDestroyNotify data_free;
  NcXcorKernelIntegrandGetKnots get_knots_func;
  NcXcorKernelIntegrandGetRangeComp get_range_comp_func;
  NcXcorKernelIntegrandEvalComps eval_comps_func;
  gdouble reltol;
  gdouble scaled_abstol;
};

struct _NcXcorKernelClass
{
  /*< private >*/
  NcmModelClass parent_class;
  /* Original XcorKernel interface */
  void (*get_z_range) (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
  gdouble (*eval_limber_z) (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
  gdouble (*eval_limber_z_prefactor) (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
  void (*prepare) (NcXcorKernel *xclk, NcHICosmo *cosmo);
  void (*add_noise) (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
  guint (*obs_len) (NcXcorKernel *xclk);
  guint (*obs_params_len) (NcXcorKernel *xclk);
  /* End of original XcorKernel interface */
  GPtrArray *(*get_component_list) (NcXcorKernel *xclk);
};

/**
 * NcXcorKernelImpl:
 * @NC_XCOR_KERNEL_IMPL_EVAL_RADIAL_WEIGHT: implementation flag for kernel evaluation method
 * @NC_XCOR_KERNEL_IMPL_PREPARE: implementation flag for kernel preparation method
 * @NC_XCOR_KERNEL_IMPL_ADD_NOISE: implementation flag for noise addition method
 *
 */
typedef enum _NcXcorKernelImpl
{
  NC_XCOR_KERNEL_IMPL_EVAL_RADIAL_WEIGHT = 0,
  NC_XCOR_KERNEL_IMPL_PREPARE,
  NC_XCOR_KERNEL_IMPL_ADD_NOISE,
  /* < private > */
} NcXcorKernelImpl;

#define NC_XCOR_KERNEL_IMPL_ALL NCM_MODEL_CLASS_IMPL_ALL

/**
 * NcXcorKinetic:
 * @xi_z: comoving distance $\xi(z)$ at redshift $z$
 * @E_z: normalized Hubble function $E(z) = H(z)/H_0$ at redshift $z$
 *
 * A boxed type for the kinetic quantities necessary to compute the kernels.
 *
 */
struct _NcXcorKinetic
{
  gdouble xi_z;
  gdouble E_z;
};

NCM_MSET_MODEL_DECLARE_ID (nc_xcor_kernel);

NcXcorKernel *nc_xcor_kernel_ref (NcXcorKernel *xclk);
void nc_xcor_kernel_free (NcXcorKernel *xclk);
void nc_xcor_kernel_clear (NcXcorKernel **xclk);

NcXcorKinetic *nc_xcor_kinetic_copy (NcXcorKinetic *xck);
void nc_xcor_kinetic_free (NcXcorKinetic *xck);

guint nc_xcor_kernel_obs_len (NcXcorKernel *xclk);
guint nc_xcor_kernel_obs_params_len (NcXcorKernel *xclk);

guint nc_xcor_kernel_get_lmax (NcXcorKernel *xclk);
void nc_xcor_kernel_set_lmax (NcXcorKernel *xclk, guint lmax);

gint nc_xcor_kernel_get_l_limber (NcXcorKernel *xclk);
void nc_xcor_kernel_set_l_limber (NcXcorKernel *xclk, gint l_limber);

gdouble nc_xcor_kernel_get_adaptive_epsilon (NcXcorKernel *xclk);
void nc_xcor_kernel_set_adaptive_epsilon (NcXcorKernel *xclk, gdouble adaptive_epsilon);

guint nc_xcor_kernel_get_adaptive_boundary_tries (NcXcorKernel *xclk);
void nc_xcor_kernel_set_adaptive_boundary_tries (NcXcorKernel *xclk, guint adaptive_boundary_tries);

gdouble nc_xcor_kernel_get_reltol (NcXcorKernel *xclk);
void nc_xcor_kernel_set_reltol (NcXcorKernel *xclk, gdouble reltol);

gdouble nc_xcor_kernel_get_scaled_abstol (NcXcorKernel *xclk);
void nc_xcor_kernel_set_scaled_abstol (NcXcorKernel *xclk, gdouble scaled_abstol);

guint nc_xcor_kernel_get_max_border_expansions (NcXcorKernel *xclk);
void nc_xcor_kernel_set_max_border_expansions (NcXcorKernel *xclk, guint max_border_expansions);

guint nc_xcor_kernel_get_max_iter (NcXcorKernel *xclk);
void nc_xcor_kernel_set_max_iter (NcXcorKernel *xclk, guint max_iter);

gdouble nc_xcor_kernel_get_expansion_factor (NcXcorKernel *xclk);
void nc_xcor_kernel_set_expansion_factor (NcXcorKernel *xclk, gdouble expansion_factor);

NcDistance *nc_xcor_kernel_peek_dist (NcXcorKernel *xclk);
NcmPowspec *nc_xcor_kernel_peek_powspec (NcXcorKernel *xclk);
NcmSBesselIntegrator *nc_xcor_kernel_peek_integrator (NcXcorKernel *xclk);

void nc_xcor_kernel_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
void nc_xcor_kernel_get_k_range (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);
NcXcorKernelIntegrand *nc_xcor_kernel_get_eval (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
NcXcorKernelIntegrand *nc_xcor_kernel_get_eval_vectorized (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax);
NcXcorKernelIntegrand *nc_xcor_kernel_get_eval_vectorized_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gint lmin, gint lmax, NcmSBesselIntegrator *sbi);

gdouble nc_xcor_kernel_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
gdouble nc_xcor_kernel_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
gdouble nc_xcor_kernel_eval_limber_z_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l);

void nc_xcor_kernel_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
void nc_xcor_kernel_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);

GPtrArray *nc_xcor_kernel_get_component_list (NcXcorKernel *xclk);

void nc_xcor_kernel_log_all_models (void);

GType nc_xcor_kernel_integrand_get_type (void) G_GNUC_CONST;

NcXcorKernelIntegrand *nc_xcor_kernel_integrand_new (guint len, NcXcorKernelIntegrandEval eval, NcXcorKernelIntegrandGetRange get_range, gpointer data, GDestroyNotify data_free);
void nc_xcor_kernel_integrand_set_get_knots (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetKnots get_knots);
void nc_xcor_kernel_integrand_set_get_range_comp (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandGetRangeComp get_range_comp);
void nc_xcor_kernel_integrand_set_eval_comps (NcXcorKernelIntegrand *integrand, NcXcorKernelIntegrandEvalComps eval_comps);
NcmVector *nc_xcor_kernel_integrand_peek_knots (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_set_tolerances (NcXcorKernelIntegrand *integrand, gdouble reltol, gdouble scaled_abstol);
gdouble nc_xcor_kernel_integrand_get_reltol (NcXcorKernelIntegrand *integrand);
gdouble nc_xcor_kernel_integrand_get_scaled_abstol (NcXcorKernelIntegrand *integrand);
NcXcorKernelIntegrand *nc_xcor_kernel_integrand_ref (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_unref (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_clear (NcXcorKernelIntegrand **integrand);

NCM_INLINE guint nc_xcor_kernel_integrand_get_len (NcXcorKernelIntegrand *integrand);
NCM_INLINE void nc_xcor_kernel_integrand_eval (NcXcorKernelIntegrand *integrand, gdouble k, gdouble *W);
NCM_INLINE void nc_xcor_kernel_integrand_get_range (NcXcorKernelIntegrand *integrand, gdouble *k_min, gdouble *k_max);
NCM_INLINE void nc_xcor_kernel_integrand_get_range_comp (NcXcorKernelIntegrand *integrand, guint i, gdouble *k_min, gdouble *k_max);
NCM_INLINE void nc_xcor_kernel_integrand_eval_comps (NcXcorKernelIntegrand *integrand, gdouble k, guint offset, guint len, gdouble *W);
NCM_INLINE GArray *nc_xcor_kernel_integrand_eval_array (NcXcorKernelIntegrand *integrand, gdouble k);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_H_ */


#ifndef _NC_XCOR_KERNEL_INLINE_H_
#define _NC_XCOR_KERNEL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__
#ifndef NUMCOSMO_GIR_SCAN

G_BEGIN_DECLS

NCM_INLINE guint
nc_xcor_kernel_integrand_get_len (NcXcorKernelIntegrand *integrand)
{
  return integrand->len;
}

NCM_INLINE void
nc_xcor_kernel_integrand_eval (NcXcorKernelIntegrand *integrand, gdouble k, gdouble *W)
{
  integrand->eval_func (integrand->data, k, W);
}

NCM_INLINE void
nc_xcor_kernel_integrand_get_range (NcXcorKernelIntegrand *integrand, gdouble *k_min, gdouble *k_max)
{
  integrand->get_range_func (integrand->data, k_min, k_max);
}

NCM_INLINE void
nc_xcor_kernel_integrand_get_range_comp (NcXcorKernelIntegrand *integrand, guint i, gdouble *k_min, gdouble *k_max)
{
  if (integrand->get_range_comp_func != NULL)
    integrand->get_range_comp_func (integrand->data, i, k_min, k_max);
  else
    integrand->get_range_func (integrand->data, k_min, k_max);
}

/* Fills only W[offset, offset + len) when the integrand can evaluate a run of
 * components on its own, and the whole of W when it cannot -- either way those
 * entries are what the caller reads. */
NCM_INLINE void
nc_xcor_kernel_integrand_eval_comps (NcXcorKernelIntegrand *integrand, gdouble k, guint offset, guint len, gdouble *W)
{
  if (integrand->eval_comps_func != NULL)
    integrand->eval_comps_func (integrand->data, k, offset, len, W);
  else
    integrand->eval_func (integrand->data, k, W);
}

NCM_INLINE GArray *
nc_xcor_kernel_integrand_eval_array (NcXcorKernelIntegrand *integrand, gdouble k)
{
  GArray *arr = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), integrand->len);

  g_array_set_size (arr, integrand->len);
  integrand->eval_func (integrand->data, k, (gdouble *) arr->data);

  return arr;
}

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */
#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_XCOR_KERNEL_INLINE_H_ */


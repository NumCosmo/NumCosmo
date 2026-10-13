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

#include <numcosmo/ncm/algebra/ncm_spectral_panels.h>
#include <numcosmo/ncm/algebra/ncm_spectral_panels_builder.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/powspec/nc_growth_func.h>
#include <numcosmo/nc/powspec/nc_transfer_func.h>
#include <numcosmo/ncm/core/ncm_c.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/powspec/ncm_powspec.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/background/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL (nc_xcor_kernel_get_type ())
#define NC_TYPE_XCOR_KERNEL_INTEGRAND (nc_xcor_kernel_integrand_get_type ())

/**
 * NC_XCOR_KERNEL_MAX_ELL_BLOCK:
 *
 * Hard cap on the number of multipoles in a single get_eval_vectorized() call:
 * #NcXcorKernel's internal per-block state uses fixed-size stack arrays sized by this
 * constant (not just the Levin integrator's own, larger ell_cache_max). Exceeding it
 * is a fatal g_error.
 */
#define NC_XCOR_KERNEL_MAX_ELL_BLOCK 64

/**
 * NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON:
 *
 * Smallest #NcXcorKernel:peak-epsilon worth asking for. The tolerance is a fraction of
 * the peak of $W_i(k)$, but the quantity integrated to form $C_\ell$ is $k^2 W_i W_j$,
 * so it enters squared: this floor is $10^{-12}$ on the integrand, already past what
 * the outer $k$ integral carries. Below it nc_xcor_kernel_set_peak_epsilon() warns.
 */
#define NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON (1.0e-6)

G_DECLARE_DERIVABLE_TYPE (NcXcorKernel, nc_xcor_kernel, NC, XCOR_KERNEL, NcmModel);

typedef struct _NcXcorKinetic NcXcorKinetic;
typedef struct _NcXcorKernelIntegrand NcXcorKernelIntegrand;

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
typedef enum _NcXcorKernelImpl /*< prefix=NC_XCOR_KERNEL_IMPL >*/
{
  NC_XCOR_KERNEL_IMPL_EVAL_RADIAL_WEIGHT = 0,
  NC_XCOR_KERNEL_IMPL_PREPARE,
  NC_XCOR_KERNEL_IMPL_ADD_NOISE,
  /* < private > */
} NcXcorKernelImpl;

#define NC_XCOR_KERNEL_IMPL_ALL NCM_MODEL_CLASS_IMPL_ALL

/**
 * NcXcorKinetic:
 * @chi_z: comoving distance $\chi(z)$ at redshift $z$, in units of the Hubble radius
 * @E_z: normalized Hubble function $E(z) = H(z)/H_0$ at redshift $z$
 * @z: the redshift
 *
 * The point along the line of sight at which a kernel is evaluated: the
 * quantities shared by every kernel component and every multipole at that point.
 */
struct _NcXcorKinetic
{
  gdouble chi_z;
  gdouble E_z;
  gdouble z;
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

gdouble nc_xcor_kernel_get_peak_epsilon (NcXcorKernel *xclk);
void nc_xcor_kernel_set_peak_epsilon (NcXcorKernel *xclk, gdouble peak_epsilon);

guint nc_xcor_kernel_get_max_border_expansions (NcXcorKernel *xclk);
void nc_xcor_kernel_set_max_border_expansions (NcXcorKernel *xclk, guint max_border_expansions);

gdouble nc_xcor_kernel_get_expansion_factor (NcXcorKernel *xclk);
void nc_xcor_kernel_set_expansion_factor (NcXcorKernel *xclk, gdouble expansion_factor);
NcmSpectralPanelsBuilder *nc_xcor_kernel_peek_panels_builder (NcXcorKernel *xclk);
void nc_xcor_kernel_set_panels_builder (NcXcorKernel *xclk, NcmSpectralPanelsBuilder *builder);

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
void nc_xcor_kernel_prepare_if_needed (NcXcorKernel *xclk, NcHICosmo *cosmo);
void nc_xcor_kernel_mark_outdated (NcXcorKernel *xclk);
void nc_xcor_kernel_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);

GPtrArray *nc_xcor_kernel_get_component_list (NcXcorKernel *xclk);

void nc_xcor_kernel_log_all_models (void);

GType nc_xcor_kernel_integrand_get_type (void) G_GNUC_CONST;

NcXcorKernelIntegrand *nc_xcor_kernel_integrand_ref (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_unref (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_clear (NcXcorKernelIntegrand **integrand);

guint nc_xcor_kernel_integrand_get_len (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_get_range (NcXcorKernelIntegrand *integrand, gdouble *k_min, gdouble *k_max);
void nc_xcor_kernel_integrand_get_range_comp (NcXcorKernelIntegrand *integrand, guint i, gdouble *k_min, gdouble *k_max);
gdouble nc_xcor_kernel_integrand_get_scale (NcXcorKernelIntegrand *integrand, guint i);
NcmVector *nc_xcor_kernel_integrand_peek_scales (NcXcorKernelIntegrand *integrand);
NcmSpectralPanels *nc_xcor_kernel_integrand_peek_panels (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_eval (NcXcorKernelIntegrand *integrand, gdouble k, gdouble *W);
void nc_xcor_kernel_integrand_eval_comps (NcXcorKernelIntegrand *integrand, gdouble k, guint offset, guint len, gdouble *W);
GArray *nc_xcor_kernel_integrand_eval_array (NcXcorKernelIntegrand *integrand, gdouble k);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_H_ */


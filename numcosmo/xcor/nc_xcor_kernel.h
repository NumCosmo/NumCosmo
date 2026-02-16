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
#include <numcosmo/lss/nc_growth_func.h>
#include <numcosmo/lss/nc_transfer_func.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_sbessel_integrator.h>
#include <numcosmo/math/ncm_util.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL (nc_xcor_kernel_get_type ())
#define NC_TYPE_XCOR_KERNEL_INTEGRAND (nc_xcor_kernel_integrand_get_type ())

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
 * NcXcorKernelIntegrand:
 * @refcount: atomic reference count
 * @len: number of components in the integrand
 * @eval_func: function to evaluate the integrand at @k, filling @W[@len]
 * @get_range_func: function to get the valid k range for this integrand
 * @data: user data passed to @eval_func and @get_range_func
 * @data_free: function to free @data, or %NULL if no cleanup needed
 *
 * A reference-counted closure for computing kernel integrands.
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
  /* New XcorKernel interface - build closure and closure interval*/
  void (*get_k_range) (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);
  NcXcorKernelIntegrand *(*get_eval) (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
  /* Proposed new interface for Components*/
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

NcDistance *nc_xcor_kernel_peek_dist (NcXcorKernel *xclk);
NcmPowspec *nc_xcor_kernel_peek_powspec (NcXcorKernel *xclk);
NcmSBesselIntegrator *nc_xcor_kernel_peek_integrator (NcXcorKernel *xclk);

void nc_xcor_kernel_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
void nc_xcor_kernel_get_k_range (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);
NcXcorKernelIntegrand *nc_xcor_kernel_get_eval (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);

gdouble nc_xcor_kernel_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
gdouble nc_xcor_kernel_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
gdouble nc_xcor_kernel_eval_limber_z_full (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l);

void nc_xcor_kernel_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
void nc_xcor_kernel_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);

void nc_xcor_kernel_log_all_models (void);

GType nc_xcor_kernel_integrand_get_type (void) G_GNUC_CONST;

NcXcorKernelIntegrand *nc_xcor_kernel_integrand_new (guint len, NcXcorKernelIntegrandEval eval, NcXcorKernelIntegrandGetRange get_range, gpointer data, GDestroyNotify data_free);
NcXcorKernelIntegrand *nc_xcor_kernel_integrand_ref (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_unref (NcXcorKernelIntegrand *integrand);
void nc_xcor_kernel_integrand_clear (NcXcorKernelIntegrand **integrand);

NCM_INLINE guint nc_xcor_kernel_integrand_get_len (NcXcorKernelIntegrand *integrand);
NCM_INLINE void nc_xcor_kernel_integrand_eval (NcXcorKernelIntegrand *integrand, gdouble k, gdouble *W);
NCM_INLINE void nc_xcor_kernel_integrand_get_range (NcXcorKernelIntegrand *integrand, gdouble *k_min, gdouble *k_max);
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


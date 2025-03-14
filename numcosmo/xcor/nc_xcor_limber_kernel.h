/***************************************************************************
 *            nc_xcor_limber_kernel.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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

#ifndef _NC_XCOR_LIMBER_KERNEL_H_
#define _NC_XCOR_LIMBER_KERNEL_H_

#include <glib-object.h>
#include <glib.h>

#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_growth_func.h>
#include <numcosmo/lss/nc_transfer_func.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER_KERNEL (nc_xcor_limber_kernel_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcXcorLimberKernel, nc_xcor_limber_kernel, NC, XCOR_LIMBER_KERNEL, NcmModel);

typedef struct _NcXcorKinetic NcXcorKinetic;

struct _NcXcorLimberKernelClass
{
  /*< private >*/
  NcmModelClass parent_class;

  gdouble (*eval) (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
  void (*prepare) (NcXcorLimberKernel *xclk, NcHICosmo *cosmo);
  void (*add_noise) (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
  guint (*obs_len) (NcXcorLimberKernel *xclk);
  guint (*obs_params_len) (NcXcorLimberKernel *xclk);
};


/**
 * NcXcorLimberKernelImpl:
 * @NC_XCOR_LIMBER_KERNEL_IMPL_EVAL: FIXME
 * @NC_XCOR_LIMBER_KERNEL_IMPL_PREPARE: FIXME
 * @NC_XCOR_LIMBER_KERNEL_IMPL_ADD_NOISE: FIXME
 *
 */
typedef enum _NcXcorLimberKernelImpl
{
  NC_XCOR_LIMBER_KERNEL_IMPL_EVAL = 0,
  NC_XCOR_LIMBER_KERNEL_IMPL_PREPARE,
  NC_XCOR_LIMBER_KERNEL_IMPL_ADD_NOISE,
  /* < private > */
} NcXcorLimberKernelImpl;

#define NC_XCOR_LIMBER_KERNEL_IMPL_ALL NCM_MODEL_CLASS_IMPL_ALL

/**
 * NcXcorKinetic:
 * @xi_z: FIXME
 * @E_z: FIXME
 *
 * A boxed type for the kinetic quantities necessary to compute the kernels.
 *
 */
struct _NcXcorKinetic
{
  gdouble xi_z;
  gdouble E_z;
};

NCM_MSET_MODEL_DECLARE_ID (nc_xcor_limber_kernel);

NcXcorLimberKernel *nc_xcor_limber_kernel_ref (NcXcorLimberKernel *xclk);
void nc_xcor_limber_kernel_free (NcXcorLimberKernel *xclk);
void nc_xcor_limber_kernel_clear (NcXcorLimberKernel **xclk);

NcXcorKinetic *nc_xcor_kinetic_copy (NcXcorKinetic *xck);
void nc_xcor_kinetic_free (NcXcorKinetic *xck);

guint nc_xcor_limber_kernel_obs_len (NcXcorLimberKernel *xclk);
guint nc_xcor_limber_kernel_obs_params_len (NcXcorLimberKernel *xclk);

void nc_xcor_limber_kernel_set_z_range (NcXcorLimberKernel *xclk, gdouble zmin, gdouble zmax, gdouble zmid);
void nc_xcor_limber_kernel_get_z_range (NcXcorLimberKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);

void nc_xcor_limber_kernel_set_const_factor (NcXcorLimberKernel *xclk, gdouble cf);
gdouble nc_xcor_limber_kernel_get_const_factor (NcXcorLimberKernel *xclk);

gdouble nc_xcor_limber_kernel_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
gdouble nc_xcor_limber_kernel_eval_full (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, NcDistance *dist, gint l);
void nc_xcor_limber_kernel_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo);
void nc_xcor_limber_kernel_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);

void nc_xcor_limber_kernel_log_all_models (void);

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_H_ */


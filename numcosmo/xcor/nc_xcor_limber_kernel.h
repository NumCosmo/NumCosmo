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
#define NC_XCOR_LIMBER_KERNEL(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR_LIMBER_KERNEL, NcXcorLimberKernel))
#define NC_XCOR_LIMBER_KERNEL_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR_LIMBER_KERNEL, NcXcorLimberKernelClass))
#define NC_IS_XCOR_LIMBER_KERNEL(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR_LIMBER_KERNEL))
#define NC_IS_XCOR_LIMBER_KERNEL_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR_LIMBER_KERNEL))
#define NC_XCOR_LIMBER_KERNEL_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR_LIMBER_KERNEL, NcXcorLimberKernelClass))

typedef struct _NcXcorLimberKernelClass NcXcorLimberKernelClass;
typedef struct _NcXcorLimberKernel NcXcorLimberKernel;

/**
 * NcXcorLimberKernelImpl:
 * @NC_XCOR_LIMBER_KERNEL_IMPL_EVAL: FIXME
 * @NC_XCOR_LIMBER_KERNEL_IMPL_PREPARE: FIXME
 * @NC_XCOR_LIMBER_KERNEL_IMPL_ADD_NOISE: FIXME
 *
 */
typedef enum _NcXcorLimberKernelImpl
{
	NC_XCOR_LIMBER_KERNEL_IMPL_EVAL      = 1 << 0,
	NC_XCOR_LIMBER_KERNEL_IMPL_PREPARE   = 1 << 1,
	NC_XCOR_LIMBER_KERNEL_IMPL_ADD_NOISE = 1 << 2, /*< private >*/
} NcXcorLimberKernelImpl;

#define NC_XCOR_LIMBER_KERNEL_IMPL_ALL (~0)

struct _NcXcorLimberKernelClass
{
	/*< private >*/
	NcmModelClass parent_class;
	gdouble (*eval) (NcXcorLimberKernel* xclk, NcHICosmo* cosmo, gdouble z, gint l);//, gdouble geo_z[]);
	void (*prepare) (NcXcorLimberKernel* xclk, NcHICosmo* cosmo);
	void (*add_noise) (NcXcorLimberKernel* xclk, NcmVector* vp1, NcmVector* vp2, guint lmin);
	guint (*obs_len) (NcXcorLimberKernel* xclk);
	guint (*obs_params_len) (NcXcorLimberKernel* xclk);
	NcXcorLimberKernelImpl impl;
};

struct _NcXcorLimberKernel
{
	/*< private >*/
	NcmModel parent_instance;
	gdouble cons_factor;
	gdouble zmin, zmax;
};

GType nc_xcor_limber_kernel_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_xcor_limber_kernel);

NcXcorLimberKernel* nc_xcor_limber_kernel_new_from_name (gchar* xcor_name);
NcXcorLimberKernel* nc_xcor_limber_kernel_ref (NcXcorLimberKernel* xclk);
void nc_xcor_limber_kernel_free (NcXcorLimberKernel* xclk);
void nc_xcor_limber_kernel_clear (NcXcorLimberKernel** xclk);
// void nc_xcor_limber_kernel_constructed (NcXcorLimberKernel** xclk);

NcXcorLimberKernelImpl nc_xcor_limber_kernel_impl (NcXcorLimberKernel* xclk);

guint nc_xcor_limber_kernel_obs_len (NcXcorLimberKernel* xclk);
guint nc_xcor_limber_kernel_obs_params_len (NcXcorLimberKernel* xclk);

gdouble nc_xcor_limber_kernel_eval (NcXcorLimberKernel* xclk, NcHICosmo* cosmo, gdouble z, gint l);//, gdouble geo_z[]);
void nc_xcor_limber_kernel_prepare (NcXcorLimberKernel* xclk, NcHICosmo* cosmo);
// gdouble nc_xcor_limber_kernel_noise_spec (NcXcorLimberKernel* xclk, guint l);
void nc_xcor_limber_kernel_add_noise (NcXcorLimberKernel* xclk, NcmVector* vp1, NcmVector* vp2, guint lmin);

void nc_xcor_limber_kernel_log_all_models (void);


G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_H_ */

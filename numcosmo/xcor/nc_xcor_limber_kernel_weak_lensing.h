/***************************************************************************
 *            nc_xcor_limber_kernel_weak_lensing.h
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

#ifndef _NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_H_
#define _NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER_KERNEL_WEAK_LENSING (nc_xcor_limber_kernel_weak_lensing_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorLimberKernelWeakLensing, nc_xcor_limber_kernel_weak_lensing, NC, XCOR_LIMBER_KERNEL_WEAK_LENSING, NcXcorLimberKernel)

/**
 * NcXcorLimberKernelWeakLensingSParams:
 * @NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_SPARAM_LEN: FIXME
 *
 */
typedef enum _NcXcorLimberKernelWeakLensingSParams
{
  NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_SPARAM_LEN,
} NcXcorLimberKernelWeakLensingSParams;

/**
 * NcXcorLimberKernelWeakLensingVParams:
 * @NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_VPARAM_LEN: FIXME
 *
 */
typedef enum _NcXcorLimberKernelWeakLensingVParams
{
  NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_VPARAM_LEN,
} NcXcorLimberKernelWeakLensingVParams;

/* #define NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_BIAS_DEFAULT_LEN (1) */
/* #define NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_DEFAULT_BIAS (1.0) */

/* #define NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_DEFAULT_MAG_BIAS (0.4) */
/* #define NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_DEFAULT_NOISE_BIAS (0.0) */

/* #define NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN (200) */

#define NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_DEFAULT_PARAMS_ABSTOL (0.0)

NcXcorLimberKernelWeakLensing *nc_xcor_limber_kernel_weak_lensing_new (gdouble zmin, gdouble zmax, NcmSpline *dn_dz, gdouble nbar, gdouble intr_shear, NcDistance *dist);

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_H_ */


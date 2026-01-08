/***************************************************************************
 *            nc_xcor_kernel_weak_lensing.h
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

#ifndef _NC_XCOR_KERNEL_WEAK_LENSING_H_
#define _NC_XCOR_KERNEL_WEAK_LENSING_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_kernel.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_WEAK_LENSING (nc_xcor_kernel_weak_lensing_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelWeakLensing, nc_xcor_kernel_weak_lensing, NC, XCOR_KERNEL_WEAK_LENSING, NcXcorKernel)

/**
 * NcXcorKernelWeakLensingSParams:
 * @NC_XCOR_KERNEL_WEAK_LENSING_SPARAM_LEN: number of scalar parameters
 *
 */
typedef enum _NcXcorKernelWeakLensingSParams
{
  NC_XCOR_KERNEL_WEAK_LENSING_SPARAM_LEN,
} NcXcorKernelWeakLensingSParams;

/**
 * NcXcorKernelWeakLensingVParams:
 * @NC_XCOR_KERNEL_WEAK_LENSING_VPARAM_LEN: number of vector parameters
 *
 */
typedef enum _NcXcorKernelWeakLensingVParams
{
  NC_XCOR_KERNEL_WEAK_LENSING_VPARAM_LEN,
} NcXcorKernelWeakLensingVParams;

/* #define NC_XCOR_KERNEL_WEAK_LENSING_BIAS_DEFAULT_LEN (1) */
/* #define NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_BIAS (1.0) */

/* #define NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_MAG_BIAS (0.4) */
/* #define NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_NOISE_BIAS (0.0) */

/* #define NC_XCOR_KERNEL_WEAK_LENSING_G_FUNC_LEN (200) */

#define NC_XCOR_KERNEL_WEAK_LENSING_DEFAULT_PARAMS_ABSTOL (0.0)

NcXcorKernelWeakLensing *nc_xcor_kernel_weak_lensing_new (NcDistance *dist, NcmPowspec *ps, NcmSpline *dn_dz, gdouble nbar, gdouble intr_shear);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_WEAK_LENSING_H_ */


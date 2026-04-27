/***************************************************************************
 *            nc_xcor_kernel_CMB_lensing.h
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

#ifndef _NC_XCOR_KERNEL_CMB_LENSING_H_
#define _NC_XCOR_KERNEL_CMB_LENSING_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_kernel.h>
#include <numcosmo/nc_recomb.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_CMB_LENSING (nc_xcor_kernel_cmb_lensing_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelCMBLensing, nc_xcor_kernel_cmb_lensing, NC, XCOR_KERNEL_CMB_LENSING, NcXcorKernel)

/**
 * NcXcorKernelCMBLensingSParams:
 * @NC_XCOR_KERNEL_CMB_LENSING_SPARAM_LEN: number of scalar parameters
 *
 * Scalar parameters for CMB lensing kernel (currently none defined).
 */
typedef enum /*< enum,underscore_name=NC_XCOR_KERNEL_CMB_LENSING_SPARAMS >*/
{
  NC_XCOR_KERNEL_CMB_LENSING_SPARAM_LEN,
} NcXcorKernelCMBLensingSParams;

#define NC_XCOR_KERNEL_CMB_LENSING_DEFAULT_PARAMS_ABSTOL (0.0)

NcXcorKernelCMBLensing *nc_xcor_kernel_cmb_lensing_new (NcDistance *dist, NcmPowspec *ps, NcRecomb *recomb, NcmVector *Nl);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_CMB_LENSING_H_ */


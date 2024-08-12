/***************************************************************************
 *            nc_xcor_limber_kernel_cmb_isw.h
 *
 *  Tue Sept 28 17:17:26 2021
 *  Copyright  2021  Mariana Penna-Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna-Lima  <pennalima@gmail.com>
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

#ifndef _NC_XCOR_LIMBER_KERNEL_CMB_ISW_H_
#define _NC_XCOR_LIMBER_KERNEL_CMB_ISW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>
#include <numcosmo/nc_recomb.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW (nc_xcor_limber_kernel_cmb_isw_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorLimberKernelCMBISW, nc_xcor_limber_kernel_cmb_isw, NC, XCOR_LIMBER_KERNEL_CMB_ISW, NcXcorLimberKernel);

NcXcorLimberKernelCMBISW *nc_xcor_limber_kernel_cmb_isw_new (NcDistance *dist, NcmPowspec *ps, NcRecomb *recomb, NcmVector *Nl);

#define NC_XCOR_LIMBER_KERNEL_CMB_ISW_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_CMB_ISW_H_ */


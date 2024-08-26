/***************************************************************************
 *            nc_xcor_limber_kernel_gal.h
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

#ifndef _NC_XCOR_LIMBER_KERNEL_GAL_H_
#define _NC_XCOR_LIMBER_KERNEL_GAL_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER_KERNEL_GAL (nc_xcor_limber_kernel_gal_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorLimberKernelGal, nc_xcor_limber_kernel_gal, NC, XCOR_LIMBER_KERNEL_GAL, NcXcorLimberKernel)

/**
 * NcXcorLimberKernelGalSParams:
 * @NC_XCOR_LIMBER_KERNEL_GAL_MAG_BIAS: FIXME
 * @NC_XCOR_LIMBER_KERNEL_GAL_NOISE_BIAS: FIXME
 *
 */
typedef enum _NcXcorLimberKernelGalSParams
{
  NC_XCOR_LIMBER_KERNEL_GAL_MAG_BIAS = 0,
  NC_XCOR_LIMBER_KERNEL_GAL_NOISE_BIAS,
  /* < private > */
  NC_XCOR_LIMBER_KERNEL_GAL_SPARAM_LEN, /*< skip >*/
} NcXcorLimberKernelGalSParams;

/**
 * NcXcorLimberKernelGalVParams:
 * @NC_XCOR_LIMBER_KERNEL_GAL_BIAS: FIXME
 *
 */
typedef enum _NcXcorLimberKernelGalVParams
{
  NC_XCOR_LIMBER_KERNEL_GAL_BIAS,
  /* < private > */
  NC_XCOR_LIMBER_KERNEL_GAL_VPARAM_LEN, /*< skip >*/
} NcXcorLimberKernelGalVParams;


#define NC_XCOR_LIMBER_KERNEL_GAL_BIAS_DEFAULT_LEN (1)
#define NC_XCOR_LIMBER_KERNEL_GAL_DEFAULT_BIAS (1.0)

#define NC_XCOR_LIMBER_KERNEL_GAL_DEFAULT_MAG_BIAS (0.4)
#define NC_XCOR_LIMBER_KERNEL_GAL_DEFAULT_NOISE_BIAS (0.0)

#define NC_XCOR_LIMBER_KERNEL_GAL_G_FUNC_LEN (200)

#define NC_XCOR_LIMBER_KERNEL_GAL_DEFAULT_PARAMS_ABSTOL (0.0)

NcXcorLimberKernelGal *nc_xcor_limber_kernel_gal_new (gdouble zmin, gdouble zmax, gsize np, gdouble nbarm1, NcmSpline *dn_dz, NcDistance *dist, gboolean domagbias);

void nc_xcor_limber_kernel_gal_set_fast_update (NcXcorLimberKernelGal *xclk, gboolean fast_update);
gboolean nc_xcor_limber_kernel_gal_get_fast_update (NcXcorLimberKernelGal *xclk);

void nc_xcor_limber_kernel_gal_set_bias_old (NcXcorLimberKernelGal *xclk, gdouble bias_old, gdouble noise_bias_old);
void nc_xcor_limber_kernel_gal_get_bias (NcXcorLimberKernelGal *xclk, gdouble *bias, gdouble *bias_old, gdouble *noise_bias_old);

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_GAL_H_ */


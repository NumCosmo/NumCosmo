/***************************************************************************
 *            nc_xcor_kernel_gal.h
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

#ifndef _NC_XCOR_KERNEL_GAL_H_
#define _NC_XCOR_KERNEL_GAL_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_kernel.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_GAL (nc_xcor_kernel_gal_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelGal, nc_xcor_kernel_gal, NC, XCOR_KERNEL_GAL, NcXcorKernel)

/**
 * NcXcorKernelGalSParams:
 * @NC_XCOR_KERNEL_GAL_MAG_BIAS: magnification bias parameter $s$
 * @NC_XCOR_KERNEL_GAL_NOISE_BIAS: noise bias parameter
 *
 */
typedef enum _NcXcorKernelGalSParams
{
  NC_XCOR_KERNEL_GAL_MAG_BIAS = 0,
  NC_XCOR_KERNEL_GAL_NOISE_BIAS,
  /* < private > */
  NC_XCOR_KERNEL_GAL_SPARAM_LEN, /*< skip >*/
} NcXcorKernelGalSParams;

/**
 * NcXcorKernelGalVParams:
 * @NC_XCOR_KERNEL_GAL_BIAS: large-scale clustering bias $b(z)$
 *
 */
typedef enum _NcXcorKernelGalVParams
{
  NC_XCOR_KERNEL_GAL_BIAS,
  /* < private > */
  NC_XCOR_KERNEL_GAL_VPARAM_LEN, /*< skip >*/
} NcXcorKernelGalVParams;


#define NC_XCOR_KERNEL_GAL_BIAS_DEFAULT_LEN (1)
#define NC_XCOR_KERNEL_GAL_DEFAULT_BIAS (1.0)

#define NC_XCOR_KERNEL_GAL_DEFAULT_MAG_BIAS (0.4)
#define NC_XCOR_KERNEL_GAL_DEFAULT_NOISE_BIAS (0.0)

#define NC_XCOR_KERNEL_GAL_G_FUNC_LEN (200)

#define NC_XCOR_KERNEL_GAL_DEFAULT_PARAMS_ABSTOL (0.0)

NcXcorKernelGal *nc_xcor_kernel_gal_new (gdouble zmin, gdouble zmax, gsize np, gdouble nbarm1, NcmSpline *dn_dz, NcDistance *dist, gboolean domagbias);

void nc_xcor_kernel_gal_set_fast_update (NcXcorKernelGal *xclk, gboolean fast_update);
gboolean nc_xcor_kernel_gal_get_fast_update (NcXcorKernelGal *xclk);

void nc_xcor_kernel_gal_set_bias_old (NcXcorKernelGal *xclk, gdouble bias_old, gdouble noise_bias_old);
void nc_xcor_kernel_gal_get_bias (NcXcorKernelGal *xclk, gdouble *bias, gdouble *bias_old, gdouble *noise_bias_old);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_GAL_H_ */


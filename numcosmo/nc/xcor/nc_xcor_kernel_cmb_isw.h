/***************************************************************************
 *            nc_xcor_kernel_cmb_isw.h
 *
 *  Tue Sept 28 17:17:26 2021
 *  Copyright  2021  Mariana Penna-Lima
 *  <pennalima@gmail.com>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna-Lima  <pennalima@gmail.com>
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

#ifndef _NC_XCOR_KERNEL_CMB_ISW_H_
#define _NC_XCOR_KERNEL_CMB_ISW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/spline/ncm_spline.h>
#include <numcosmo/ncm/spline/ncm_spline_cubic_notaknot.h>
#include <numcosmo/ncm/spline/ncm_spline2d.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>
#include <numcosmo/nc/recomb/nc_recomb.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_CMB_ISW (nc_xcor_kernel_cmb_isw_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelCMBISW, nc_xcor_kernel_cmb_isw, NC, XCOR_KERNEL_CMB_ISW, NcXcorKernel);

/**
 * NcXcorKernelCMBISWSource:
 * @NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN: every CMB photon last scatters at the decoupling redshift of the #NcDistance
 * @NC_XCOR_KERNEL_CMB_ISW_SOURCE_VISIBILITY: the sources follow the visibility function of the #NcXcorKernelCMBISW:recomb object over the last-scattering shell, normalized to unit integral there
 * @NC_XCOR_KERNEL_CMB_ISW_SOURCE_VISIBILITY_REIONIZATION: the sources follow the full visibility function, the reionization bump included
 *
 * Where the CMB photons are placed along the line of sight. A photon that last
 * scattered at $\chi'$ integrates the decay of the potential over $\chi < \chi'$
 * only, so the kernel at $\chi$ is the thin-screen one times the fraction of
 * photons that last scatter beyond $\chi$: one below the shell, zero beyond it,
 * and lower by the rescattered fraction between the reionization bump and the
 * shell when reionization is included.
 */
typedef enum _NcXcorKernelCMBISWSource /*< prefix=NC_XCOR_KERNEL_CMB_ISW_SOURCE >*/
{
  NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN = 0,
  NC_XCOR_KERNEL_CMB_ISW_SOURCE_VISIBILITY,
  NC_XCOR_KERNEL_CMB_ISW_SOURCE_VISIBILITY_REIONIZATION,
} NcXcorKernelCMBISWSource;

NcXcorKernelCMBISW *nc_xcor_kernel_cmb_isw_new (NcDistance *dist, NcmPowspec *ps, NcRecomb *recomb, NcmVector *Nl);

void nc_xcor_kernel_cmb_isw_set_source (NcXcorKernelCMBISW *xcisw, NcXcorKernelCMBISWSource source);
NcXcorKernelCMBISWSource nc_xcor_kernel_cmb_isw_get_source (NcXcorKernelCMBISW *xcisw);

gdouble nc_xcor_kernel_cmb_isw_eval_k_max (NcXcorKernelCMBISW *xcisw, gdouble x);
gdouble nc_xcor_kernel_cmb_isw_eval_KL_max (NcXcorKernelCMBISW *xcisw, gdouble x);
gdouble nc_xcor_kernel_cmb_isw_eval_k_epsilon (NcXcorKernelCMBISW *xcisw, gdouble x);
void nc_xcor_kernel_cmb_isw_set_epsilon (NcXcorKernelCMBISW *xcisw, gdouble epsilon);
gdouble nc_xcor_kernel_cmb_isw_get_epsilon (NcXcorKernelCMBISW *xcisw);

#define NC_XCOR_KERNEL_CMB_ISW_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_CMB_ISW_H_ */


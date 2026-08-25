/***************************************************************************
 *            nc_xcor_kernel_analytic_power_exp.h
 *
 *  Thu August 21 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_XCOR_KERNEL_ANALYTIC_POWER_EXP_H_
#define _NC_XCOR_KERNEL_ANALYTIC_POWER_EXP_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/xcor/tests/nc_xcor_kernel_analytic.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_ANALYTIC_POWER_EXP (nc_xcor_kernel_analytic_power_exp_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelAnalyticPowerExp, nc_xcor_kernel_analytic_power_exp, NC, XCOR_KERNEL_ANALYTIC_POWER_EXP, NcXcorKernelAnalytic)

NcXcorKernelAnalyticPowerExp *nc_xcor_kernel_analytic_power_exp_new (NcDistance * dist, NcmPowspec * ps, gdouble chi_scale, gdouble alpha, gdouble beta, gdouble chi_lower, gdouble chi_upper);
NcXcorKernelAnalyticPowerExp *nc_xcor_kernel_analytic_power_exp_new_full (NcDistance *dist, NcmPowspec *ps, gdouble chi_scale, gdouble alpha, gdouble beta, gdouble chi_lower, gdouble chi_upper, NcmSBesselIntegrator *sbi);

gdouble nc_xcor_kernel_analytic_power_exp_get_chi_scale (NcXcorKernelAnalyticPowerExp *xckap);
gdouble nc_xcor_kernel_analytic_power_exp_get_alpha (NcXcorKernelAnalyticPowerExp *xckap);
gdouble nc_xcor_kernel_analytic_power_exp_get_beta (NcXcorKernelAnalyticPowerExp *xckap);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_ANALYTIC_POWER_EXP_H_ */


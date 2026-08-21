/***************************************************************************
 *            nc_xcor_kernel_analytic_multi.h
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

#ifndef _NC_XCOR_KERNEL_ANALYTIC_MULTI_H_
#define _NC_XCOR_KERNEL_ANALYTIC_MULTI_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_analytic.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_ANALYTIC_MULTI (nc_xcor_kernel_analytic_multi_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelAnalyticMulti, nc_xcor_kernel_analytic_multi, NC, XCOR_KERNEL_ANALYTIC_MULTI, NcXcorKernelAnalytic)

NcXcorKernelAnalyticMulti *nc_xcor_kernel_analytic_multi_new (NcDistance * dist, NcmPowspec * ps, NcmVector * chi_mean, NcmVector * chi_sigma, NcmVector * weight, gdouble n_sigma);
NcXcorKernelAnalyticMulti *nc_xcor_kernel_analytic_multi_new_full (NcDistance *dist, NcmPowspec *ps, NcmVector *chi_mean, NcmVector *chi_sigma, NcmVector *weight, gdouble n_sigma, NcmSBesselIntegrator *sbi);

guint nc_xcor_kernel_analytic_multi_get_n_bumps (NcXcorKernelAnalyticMulti *xckam);
gdouble nc_xcor_kernel_analytic_multi_get_n_sigma (NcXcorKernelAnalyticMulti *xckam);
NcmVector *nc_xcor_kernel_analytic_multi_peek_chi_mean (NcXcorKernelAnalyticMulti *xckam);
NcmVector *nc_xcor_kernel_analytic_multi_peek_chi_sigma (NcXcorKernelAnalyticMulti *xckam);
NcmVector *nc_xcor_kernel_analytic_multi_peek_weight (NcXcorKernelAnalyticMulti *xckam);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_ANALYTIC_MULTI_H_ */


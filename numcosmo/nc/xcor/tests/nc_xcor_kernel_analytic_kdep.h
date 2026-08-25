/***************************************************************************
 *            nc_xcor_kernel_analytic_kdep.h
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

#ifndef _NC_XCOR_KERNEL_ANALYTIC_KDEP_H_
#define _NC_XCOR_KERNEL_ANALYTIC_KDEP_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_ANALYTIC_KDEP (nc_xcor_kernel_analytic_kdep_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcXcorKernelAnalyticKDep, nc_xcor_kernel_analytic_kdep, NC, XCOR_KERNEL_ANALYTIC_KDEP, GObject)

struct _NcXcorKernelAnalyticKDepClass
{
  /*< private >*/
  GObjectClass parent_class;

  gdouble (*eval) (NcXcorKernelAnalyticKDep *kdep, gdouble chi, gdouble k);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[17];
};

NcXcorKernelAnalyticKDep *nc_xcor_kernel_analytic_kdep_ref (NcXcorKernelAnalyticKDep *kdep);
void nc_xcor_kernel_analytic_kdep_free (NcXcorKernelAnalyticKDep *kdep);
void nc_xcor_kernel_analytic_kdep_clear (NcXcorKernelAnalyticKDep **kdep);

gdouble nc_xcor_kernel_analytic_kdep_eval (NcXcorKernelAnalyticKDep *kdep, gdouble chi, gdouble k);

#define NC_TYPE_XCOR_KERNEL_ANALYTIC_KDEP_GROWTH (nc_xcor_kernel_analytic_kdep_growth_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelAnalyticKDepGrowth, nc_xcor_kernel_analytic_kdep_growth, NC, XCOR_KERNEL_ANALYTIC_KDEP_GROWTH, NcXcorKernelAnalyticKDep)

NcXcorKernelAnalyticKDepGrowth *nc_xcor_kernel_analytic_kdep_growth_new (gdouble amplitude, gdouble k_transition, gdouble chi_ref);

gdouble nc_xcor_kernel_analytic_kdep_growth_get_amplitude (NcXcorKernelAnalyticKDepGrowth *kdepg);
gdouble nc_xcor_kernel_analytic_kdep_growth_get_k_transition (NcXcorKernelAnalyticKDepGrowth *kdepg);
gdouble nc_xcor_kernel_analytic_kdep_growth_get_chi_ref (NcXcorKernelAnalyticKDepGrowth *kdepg);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_ANALYTIC_KDEP_H_ */


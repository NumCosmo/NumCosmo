/***************************************************************************
 *            nc_xcor_kernel_radial_kdep.h
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

#ifndef _NC_XCOR_KERNEL_RADIAL_KDEP_H_
#define _NC_XCOR_KERNEL_RADIAL_KDEP_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_RADIAL_KDEP (nc_xcor_kernel_radial_kdep_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcXcorKernelRadialKDep, nc_xcor_kernel_radial_kdep, NC, XCOR_KERNEL_RADIAL_KDEP, GObject)

struct _NcXcorKernelRadialKDepClass
{
  /*< private >*/
  GObjectClass parent_class;

  gdouble (*eval) (NcXcorKernelRadialKDep *kdep, gdouble chi, gdouble k);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[17];
};

NcXcorKernelRadialKDep *nc_xcor_kernel_radial_kdep_ref (NcXcorKernelRadialKDep *kdep);
void nc_xcor_kernel_radial_kdep_free (NcXcorKernelRadialKDep *kdep);
void nc_xcor_kernel_radial_kdep_clear (NcXcorKernelRadialKDep **kdep);

gdouble nc_xcor_kernel_radial_kdep_eval (NcXcorKernelRadialKDep *kdep, gdouble chi, gdouble k);

#define NC_TYPE_XCOR_KERNEL_RADIAL_KDEP_GROWTH (nc_xcor_kernel_radial_kdep_growth_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelRadialKDepGrowth, nc_xcor_kernel_radial_kdep_growth, NC, XCOR_KERNEL_RADIAL_KDEP_GROWTH, NcXcorKernelRadialKDep)

NcXcorKernelRadialKDepGrowth *nc_xcor_kernel_radial_kdep_growth_new (gdouble amplitude, gdouble k_transition, gdouble chi_ref);

gdouble nc_xcor_kernel_radial_kdep_growth_get_amplitude (NcXcorKernelRadialKDepGrowth *kdepg);
gdouble nc_xcor_kernel_radial_kdep_growth_get_k_transition (NcXcorKernelRadialKDepGrowth *kdepg);
gdouble nc_xcor_kernel_radial_kdep_growth_get_chi_ref (NcXcorKernelRadialKDepGrowth *kdepg);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_RADIAL_KDEP_H_ */


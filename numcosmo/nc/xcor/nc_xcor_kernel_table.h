/***************************************************************************
 *            nc_xcor_kernel_table.h
 *
 *  Sat August 30 2026
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_kernel_table.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_XCOR_KERNEL_TABLE_H_
#define _NC_XCOR_KERNEL_TABLE_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_radial.h>
#include <numcosmo/nc/xcor/nc_xcor_component_table.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_TABLE (nc_xcor_kernel_table_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorKernelTable, nc_xcor_kernel_table, NC, XCOR_KERNEL_TABLE, NcXcorKernelRadial)

/**
 * NcXcorKernelTableKind:
 * @NC_XCOR_KERNEL_TABLE_KIND_DENSITY: the window multiplies $j_\ell(k\chi)$ directly
 * @NC_XCOR_KERNEL_TABLE_KIND_SHEAR: the window carries the shear factors
 *
 * What the tabulated window is a window *of*, which fixes the factors applied
 * on top of it. %NC_XCOR_KERNEL_TABLE_KIND_SHEAR reproduces CCL's
 * `der_bessel = -1`, `der_angles = 2`: the kernel gains $1/(k\chi)^2$ and the
 * multipole gains $\sqrt{(\ell+2)(\ell+1)\ell(\ell-1)}$.
 *
 */
NcXcorKernelTable *nc_xcor_kernel_table_new (NcDistance * dist, NcmPowspec * ps, NcmVector * chi, NcmVector * W);
NcXcorKernelTable *nc_xcor_kernel_table_new_full (NcDistance *dist, NcmPowspec *ps, NcmVector *chi, NcmVector *W, NcXcorKernelTableKind kind, guint order, gboolean normalize, NcmSBesselIntegrator *sbi);
NcXcorKernelTable *nc_xcor_kernel_table_new_from_components (NcDistance *dist, NcmPowspec *ps, NcmObjArray *components, NcmSBesselIntegrator *sbi);

guint nc_xcor_kernel_table_get_n_components (NcXcorKernelTable *xckt);
NcXcorComponentTable *nc_xcor_kernel_table_peek_component (NcXcorKernelTable *xckt, guint i);

NcXcorKernelTableKind nc_xcor_kernel_table_get_kind (NcXcorKernelTable *xckt);
guint nc_xcor_kernel_table_get_order (NcXcorKernelTable *xckt);
gboolean nc_xcor_kernel_table_get_normalize (NcXcorKernelTable *xckt);
gdouble nc_xcor_kernel_table_get_norm (NcXcorKernelTable *xckt);
NcmSpline *nc_xcor_kernel_table_peek_spline (NcXcorKernelTable *xckt);
NcmVector *nc_xcor_kernel_table_peek_knots (NcXcorKernelTable *xckt);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_TABLE_H_ */


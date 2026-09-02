/***************************************************************************
 *            nc_xcor_component_table.h
 *
 *  Wed September 02 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_component_table.h
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

#ifndef _NC_XCOR_COMPONENT_TABLE_H_
#define _NC_XCOR_COMPONENT_TABLE_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/spline/ncm_spline.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_COMPONENT_TABLE (nc_xcor_component_table_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorComponentTable, nc_xcor_component_table, NC, XCOR_COMPONENT_TABLE, GObject)

/**
 * NcXcorKernelTableKind:
 * @NC_XCOR_KERNEL_TABLE_KIND_DENSITY: a density window, weighted by $j_\ell(k\chi)$
 * @NC_XCOR_KERNEL_TABLE_KIND_SHEAR: a lensing window read as shear, weighted by $j_\ell(k\chi)/(k\chi)^2$ with the spin-2 prefactor $\sqrt{(\ell+2)(\ell+1)\ell(\ell-1)}$
 * @NC_XCOR_KERNEL_TABLE_KIND_CONVERGENCE: a lensing window read as convergence, weighted by $j_\ell(k\chi)/(k\chi)^2$ with the prefactor $\ell(\ell+1)$
 * @NC_XCOR_KERNEL_TABLE_KIND_RSD: a redshift-space distortion window, weighted by $j_\ell''(k\chi)$
 *
 * What a tabulated window is a window of, which fixes the Bessel weight, the
 * $(\chi, k)$ factor and the $\ell$ prefactor applied on top of it. The four
 * cover the components of CCL's tracers: density (`der_bessel = 0`), shear and
 * intrinsic alignments (`der_bessel = -1`, `der_angles = 2`), CMB lensing and
 * magnification (`der_bessel = -1`, `der_angles = 1`) and RSD
 * (`der_bessel = 2`).
 */
typedef enum _NcXcorKernelTableKind /*< prefix=NC_XCOR_KERNEL_TABLE_KIND >*/
{
  NC_XCOR_KERNEL_TABLE_KIND_DENSITY = 0,
  NC_XCOR_KERNEL_TABLE_KIND_SHEAR,
  NC_XCOR_KERNEL_TABLE_KIND_CONVERGENCE,
  NC_XCOR_KERNEL_TABLE_KIND_RSD,
  /* < private > */
  NC_XCOR_KERNEL_TABLE_KIND_LEN, /*< skip >*/
} NcXcorKernelTableKind;

/**
 * NC_XCOR_COMPONENT_TABLE_CHI_FLOOR:
 *
 * Lower end of the support, in Mpc, of a component whose kind carries
 * $1/(k\chi)^2$ when its table reaches closer to the origin than this.
 */
#define NC_XCOR_COMPONENT_TABLE_CHI_FLOOR (1.0e-2)

NcXcorComponentTable *nc_xcor_component_table_new (NcmVector *chi, NcmVector *W);
NcXcorComponentTable *nc_xcor_component_table_new_full (NcmVector *chi, NcmVector *W, NcXcorKernelTableKind kind, guint order, gboolean normalize);
NcXcorComponentTable *nc_xcor_component_table_ref (NcXcorComponentTable *xcct);
void nc_xcor_component_table_free (NcXcorComponentTable *xcct);
void nc_xcor_component_table_clear (NcXcorComponentTable **xcct);

NcXcorKernelTableKind nc_xcor_component_table_get_kind (NcXcorComponentTable *xcct);
guint nc_xcor_component_table_get_order (NcXcorComponentTable *xcct);
gboolean nc_xcor_component_table_get_normalize (NcXcorComponentTable *xcct);
gdouble nc_xcor_component_table_get_norm (NcXcorComponentTable *xcct);
NcmSpline *nc_xcor_component_table_peek_spline (NcXcorComponentTable *xcct);
NcmVector *nc_xcor_component_table_peek_knots (NcXcorComponentTable *xcct);

void nc_xcor_component_table_get_support (NcXcorComponentTable *xcct, gdouble *chi_min, gdouble *chi_max);
gdouble nc_xcor_component_table_eval_W (NcXcorComponentTable *xcct, gdouble chi);
gdouble nc_xcor_component_table_eval_kernel_factor (NcXcorComponentTable *xcct, gdouble chi, gdouble k);
gdouble nc_xcor_component_table_eval_prefactor (NcXcorComponentTable *xcct, gint l);
guint nc_xcor_component_table_get_bessel_deriv (NcXcorComponentTable *xcct);

G_END_DECLS

#endif /* _NC_XCOR_COMPONENT_TABLE_H_ */


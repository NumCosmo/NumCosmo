/***************************************************************************
 *            nc_xcor_kernel_radial.h
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

#ifndef _NC_XCOR_KERNEL_RADIAL_H_
#define _NC_XCOR_KERNEL_RADIAL_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_radial_kdep.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_KERNEL_RADIAL (nc_xcor_kernel_radial_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcXcorKernelRadial, nc_xcor_kernel_radial, NC, XCOR_KERNEL_RADIAL, NcXcorKernel)

/**
 * NC_XCOR_KERNEL_RADIAL_MAX_COMPS:
 *
 * Largest number of components a #NcXcorKernelRadial may declare. Each becomes
 * one #NcXcorKernelComponent, and #NcXcorKernel integrates at most six of those
 * in a single multipole block.
 */
#define NC_XCOR_KERNEL_RADIAL_MAX_COMPS 6

struct _NcXcorKernelRadialClass
{
  /*< private >*/
  NcXcorKernelClass parent_class;

  guint (*get_n_comps) (NcXcorKernelRadial *xcka);
  gdouble (*eval_W_comp) (NcXcorKernelRadial *xcka, guint comp, gdouble chi);
  void (*get_comp_support) (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);
  gdouble (*eval_kernel_factor) (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k);
  gdouble (*eval_prefactor) (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[13];
};

guint nc_xcor_kernel_radial_get_n_comps (NcXcorKernelRadial *xcka);
gdouble nc_xcor_kernel_radial_eval_W_comp (NcXcorKernelRadial *xcka, guint comp, gdouble chi);
void nc_xcor_kernel_radial_get_comp_support (NcXcorKernelRadial *xcka, guint comp, gdouble *chi_min, gdouble *chi_max);

NcXcorKernelRadialKDep *nc_xcor_kernel_radial_peek_kdep (NcXcorKernelRadial *xcka);

gdouble nc_xcor_kernel_radial_eval_kernel_factor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gdouble chi, gdouble k);
gdouble nc_xcor_kernel_radial_eval_prefactor (NcXcorKernelRadial *xcka, NcHICosmo *cosmo, gint l);

gdouble nc_xcor_kernel_radial_eval_W (NcXcorKernelRadial *xcka, gdouble chi);
void nc_xcor_kernel_radial_get_support (NcXcorKernelRadial *xcka, gdouble *chi_min, gdouble *chi_max);

G_END_DECLS

#endif /* _NC_XCOR_KERNEL_RADIAL_H_ */


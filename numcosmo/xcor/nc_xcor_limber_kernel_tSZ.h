/***************************************************************************
 *            nc_xcor_limber_kernel_tSZ.h
 *
 *  Tue January 10 12:00:00 2022
 *  Copyright  2022  Arthur de Souza Molina
 *  <arthur.souza.molina@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2022 Arthur de Souza Molina <arthur.souza.molina@uel.br>
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

#ifndef _NC_XCOR_LIMBER_KERNEL_TSZ_H_
#define _NC_XCOR_LIMBER_KERNEL_TSZ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>
#include <numcosmo/math/ncm_c.h>
G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER_KERNEL_TSZ (nc_xcor_limber_kernel_tsz_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorLimberKerneltSZ, nc_xcor_limber_kernel_tsz, NC, XCOR_LIMBER_KERNEL_TSZ, NcXcorLimberKernel)


/**
 * NcXcorLimberKerneltSZSParams:
 * @NC_XCOR_LIMBER_KERNEL_TSZ_SPARAM_LEN: Number of parameters.
 *
 * Enum values for the tSZ kernel parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_XCOR_LIMBER_KERNEL_TSZ_SPARAMS >*/
{
  NC_XCOR_LIMBER_KERNEL_TSZ_SPARAM_LEN
} NcXcorLimberKerneltSZSParams;

NcXcorLimberKerneltSZ *nc_xcor_limber_kernel_tsz_new (gdouble zmax);

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_TSZ_H_ */


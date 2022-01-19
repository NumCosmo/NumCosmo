/***************************************************************************
 *            nc_xcor_limber_kernel_tSZ.h
 *
 *  Tue January 10 12:00:00 2021
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
#define NC_XCOR_LIMBER_KERNEL_TSZ(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR_LIMBER_KERNEL_TSZ, NcXcorLimberKerneltSZ))
#define NC_XCOR_LIMBER_KERNEL_TSZ_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR_LIMBER_KERNEL_TSZ, NcXcorLimberKerneltSZClass))
#define NC_IS_XCOR_LIMBER_KERNEL_TSZ(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR_LIMBER_KERNEL_TSZ))
#define NC_IS_XCOR_LIMBER_KERNEL_TSZ_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR_LIMBER_KERNEL_TSZ))
#define NC_XCOR_LIMBER_KERNEL_TSZ_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR_LIMBER_KERNEL_TSZ, NcXcorLimberKerneltSZClass))

typedef struct _NcXcorLimberKerneltSZClass NcXcorLimberKerneltSZClass;
typedef struct _NcXcorLimberKerneltSZ NcXcorLimberKerneltSZ;

/**
 * NcXcorLimberKerneltSZSParams:
 * @NC_XCOR_LIMBER_KERNEL_TSZ_SPARAM_LEN: FIXME
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_XCOR_LIMBER_KERNEL_TSZ_SPARAMS >*/
{
  NC_XCOR_LIMBER_KERNEL_TSZ_SPARAM_LEN

} NcXcorLimberKerneltSZSParams;

struct _NcXcorLimberKerneltSZ
{
  /*< private >*/
  NcXcorLimberKernel parent_instance;

  gdouble Zmax;

  gdouble Nchi;

};

struct _NcXcorLimberKerneltSZClass
{
  /*< private >*/
  NcXcorLimberKernelClass parent_class;
};

#define NC_XCOR_LIMBER_KERNEL_TSZ_PREFAC_PARAMS (4.01710079e-06)

#define NC_XCOR_LIMBER_KERNEL_TSZ_DEFAULT_Zmax (6.0)

#define NC_XCOR_LIMBER_KERNEL_TSZ_DEFAULT_Nchi (1024)

GType nc_xcor_limber_kernel_tsz_get_type (void) G_GNUC_CONST;

NcXcorLimberKerneltSZ *nc_xcor_limber_kernel_tsz_new ( gdouble Zmax , gdouble Nchi);

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_TSZ_H_ */


/***************************************************************************
 *            nc_xcor_limber_kernel_cmb_isw.h
 *
 *  Tue Sept 28 17:17:26 2021
 *  Copyright  2021  Mariana Penna-Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna-Lima  <pennalima@gmail.com>
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

#ifndef _NC_XCOR_LIMBER_KERNEL_CMB_ISW_H_
#define _NC_XCOR_LIMBER_KERNEL_CMB_ISW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>
#include <numcosmo/nc_recomb.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW (nc_xcor_limber_kernel_cmb_isw_get_type ())
#define NC_XCOR_LIMBER_KERNEL_CMB_ISW(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW, NcXcorLimberKernelCMBISW))
#define NC_XCOR_LIMBER_KERNEL_CMB_ISW_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW, NcXcorLimberKernelCMBISWClass))
#define NC_IS_XCOR_LIMBER_KERNEL_CMB_ISW(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW))
#define NC_IS_XCOR_LIMBER_KERNEL_CMB_ISW_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW))
#define NC_XCOR_LIMBER_KERNEL_CMB_ISW_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW, NcXcorLimberKernelCMBISWClass))

typedef struct _NcXcorLimberKernelCMBISWClass NcXcorLimberKernelCMBISWClass;
typedef struct _NcXcorLimberKernelCMBISW NcXcorLimberKernelCMBISW;
typedef struct _NcXcorLimberKernelCMBISWPrivate NcXcorLimberKernelCMBISWPrivate;

/**
 * NcXcorLimberKernelCMBISWSParams:
 * @NC_XCOR_LIMBER_KERNEL_CMB_ISW_SPARAM_LEN: FIXME
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_XCOR_LIMBER_KERNEL_CMB_ISW_SPARAMS >*/
{
  NC_XCOR_LIMBER_KERNEL_CMB_ISW_SPARAM_LEN,
} NcXcorLimberKernelCMBISWSParams;

struct _NcXcorLimberKernelCMBISWClass
{
  /*< private >*/
  NcXcorLimberKernelClass parent_class;
};

struct _NcXcorLimberKernelCMBISW
{
  /*< private >*/
  NcXcorLimberKernel parent_instance;
  NcXcorLimberKernelCMBISWPrivate *priv;
};

#define NC_XCOR_LIMBER_KERNEL_CMB_ISW_DEFAULT_PARAMS_ABSTOL (0.0)

GType nc_xcor_limber_kernel_cmb_isw_get_type (void) G_GNUC_CONST;

NcXcorLimberKernelCMBISW *nc_xcor_limber_kernel_cmb_isw_new (NcDistance *dist, NcRecomb *recomb, NcmVector *Nl);

G_END_DECLS

#endif /* _NC_XCOR_LIMBER_KERNEL_CMB_ISW_H_ */


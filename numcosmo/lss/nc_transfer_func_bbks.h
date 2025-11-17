/***************************************************************************
 *            nc_transfer_func_bbks.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_TRANSFER_FUNC_BBKS_H_
#define _NC_TRANSFER_FUNC_BBKS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_transfer_func.h>

G_BEGIN_DECLS

#define NC_TYPE_TRANSFER_FUNC_BBKS (nc_transfer_func_bbks_get_type ())

G_DECLARE_FINAL_TYPE (NcTransferFuncBBKS, nc_transfer_func_bbks, NC, TRANSFER_FUNC_BBKS, NcTransferFunc)

/**
 * NcTransferFuncBBKSType:
 * @NC_TRANSFER_FUNC_BBKS_TYPE_NOBARYONS: BBKS without Baryions correction.
 * @NC_TRANSFER_FUNC_BBKS_TYPE_BARYONS: BBKS with Baryions correction.
 * @NC_TRANSFER_FUNC_BBKS_TYPE_CCL: BBKS with Baryions correction without T_gamma term (as implemented in CCL). Check why they modify it like this, is it an typo?
 *
 * BBKS transfer function variant.
 *
 */
typedef enum _NcTransferFuncBBKSType
{
  NC_TRANSFER_FUNC_BBKS_TYPE_NOBARYONS,
  NC_TRANSFER_FUNC_BBKS_TYPE_BARYONS,
  NC_TRANSFER_FUNC_BBKS_TYPE_CCL,
} NcTransferFuncBBKSType;

NcTransferFunc *nc_transfer_func_bbks_new (void);

void nc_transfer_func_bbks_set_type (NcTransferFuncBBKS *tf_bbks, NcTransferFuncBBKSType bbks_type);

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_BBKS_H_ */



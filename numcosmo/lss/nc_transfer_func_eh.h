/***************************************************************************
 *            nc_transfer_func_eh.h
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

#ifndef _NC_TRANSFER_FUNC_EH_H_
#define _NC_TRANSFER_FUNC_EH_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_transfer_func.h>

G_BEGIN_DECLS

#define NC_TYPE_TRANSFER_FUNC_EH (nc_transfer_func_eh_get_type ())

G_DECLARE_FINAL_TYPE (NcTransferFuncEH, nc_transfer_func_eh, NC, TRANSFER_FUNC_EH, NcTransferFunc)

NcTransferFunc *nc_transfer_func_eh_new (void);

void nc_transfer_func_eh_set_CCL_comp (NcTransferFuncEH *tf_eh, gboolean CCL_comp);

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_EH_H_ */


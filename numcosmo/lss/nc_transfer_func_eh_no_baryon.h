/***************************************************************************
 *            nc_transfer_func_eh_no_baryon.h
 *
 *  Mon Nov 03 17:46:27 2025
 *  Copyright  2025  Mariana Penna-Lima <pennalima@unb.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna-Lima 2025 <pennalima@unb.br>
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

#ifndef _NC_TRANSFER_FUNC_EH_NO_BARYON_H_
#define _NC_TRANSFER_FUNC_EH_NO_BARYON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_transfer_func.h>

G_BEGIN_DECLS

#define NC_TYPE_TRANSFER_FUNC_EH_NO_BARYON (nc_transfer_func_eh_no_baryon_get_type ())

G_DECLARE_FINAL_TYPE (NcTransferFuncEHNoBaryon, nc_transfer_func_eh_no_baryon, NC, TRANSFER_FUNC_EH_NO_BARYON, NcTransferFunc)

NcTransferFunc *nc_transfer_func_eh_no_baryon_new (void);

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_EH_NO_BARYON_H_ */


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

#define NC_TYPE_TRANSFER_FUNC_EH             (nc_transfer_func_eh_get_type ())
#define NC_TRANSFER_FUNC_EH(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_TRANSFER_FUNC_EH, NcTransferFuncEH))
#define NC_TRANSFER_FUNC_EH_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_TRANSFER_FUNC_EH, NcTransferFuncEHClass))
#define NC_IS_TRANSFER_FUNC_EH(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_TRANSFER_FUNC_EH))
#define NC_IS_TRANSFER_FUNC_EH_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_TRANSFER_FUNC_EH))
#define NC_TRANSFER_FUNC_EH_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_TRANSFER_FUNC_EH, NcTransferFuncEHClass))

typedef struct _NcTransferFuncEHClass NcTransferFuncEHClass;
typedef struct _NcTransferFuncEH NcTransferFuncEH;
typedef struct _NcTransferFuncEHPrivate NcTransferFuncEHPrivate;

struct _NcTransferFuncEHClass
{
  /*< private >*/
  NcTransferFuncClass parent_class;
};

struct _NcTransferFuncEH
{
  /*< private >*/
  NcTransferFunc parent_instance;
  NcTransferFuncEHPrivate *priv;
};

GType nc_transfer_func_eh_get_type (void) G_GNUC_CONST;

NcTransferFunc *nc_transfer_func_eh_new (void);

void nc_transfer_func_eh_set_CCL_comp (NcTransferFuncEH *tf_eh, gboolean CCL_comp);

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_EH_H_ */


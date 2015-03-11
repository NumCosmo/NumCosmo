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

#define NC_TYPE_TRANSFER_FUNC_BBKS             (nc_transfer_func_bbks_get_type ())
#define NC_TRANSFER_FUNC_BBKS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_TRANSFER_FUNC_BBKS, NcTransferFuncBBKS))
#define NC_TRANSFER_FUNC_BBKS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_TRANSFER_FUNC_BBKS, NcTransferFuncBBKSClass))
#define NC_IS_TRANSFER_FUNC_BBKS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_TRANSFER_FUNC_BBKS))
#define NC_IS_TRANSFER_FUNC_BBKS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_TRANSFER_FUNC_BBKS))
#define NC_TRANSFER_FUNC_BBKS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_TRANSFER_FUNC_BBKS, NcTransferFuncBBKSClass))

typedef struct _NcTransferFuncBBKSClass NcTransferFuncBBKSClass;
typedef struct _NcTransferFuncBBKS NcTransferFuncBBKS;

struct _NcTransferFuncBBKS
{
  /*< private >*/
  NcTransferFunc parent_instance;
  gdouble c1;
  gdouble c2;
  gdouble c3;
  gdouble c4;
  gdouble c5_wm; /* c5_wm = c5/wm */
  gdouble h;  
};

struct _NcTransferFuncBBKSClass
{
  /*< private >*/
  NcTransferFuncClass parent_class;
};

GType nc_transfer_func_bbks_get_type (void) G_GNUC_CONST;

NcTransferFunc *nc_transfer_func_bbks_new (void);

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_BBKS_H_ */

/***************************************************************************
 *            nc_transfer_func_camb.h
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

#ifndef _NC_TRANSFER_FUNC_CAMB_H_
#define _NC_TRANSFER_FUNC_CAMB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/lss/nc_transfer_func.h>

G_BEGIN_DECLS

#define NC_TYPE_TRANSFER_FUNC_CAMB             (nc_transfer_func_camb_get_type ())
#define NC_TRANSFER_FUNC_CAMB(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_TRANSFER_FUNC_CAMB, NcTransferFuncCAMB))
#define NC_TRANSFER_FUNC_CAMB_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_TRANSFER_FUNC_CAMB, NcTransferFuncCAMBClass))
#define NC_IS_TRANSFER_FUNC_CAMB(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_TRANSFER_FUNC_CAMB))
#define NC_IS_TRANSFER_FUNC_CAMB_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_TRANSFER_FUNC_CAMB))
#define NC_TRANSFER_FUNC_CAMB_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_TRANSFER_FUNC_CAMB, NcTransferFuncCAMBClass))

typedef struct _NcTransferFuncCAMBClass NcTransferFuncCAMBClass;
typedef struct _NcTransferFuncCAMB NcTransferFuncCAMB;

struct _NcTransferFuncCAMB
{
  /*< private >*/
  NcTransferFunc parent_instance;
  NcmSpline *T_spline;
  gboolean init;
};

struct _NcTransferFuncCAMBClass
{
  /*< private >*/
  NcTransferFuncClass parent_class;
};

GType nc_transfer_func_camb_get_type (void) G_GNUC_CONST;

NcTransferFunc *nc_transfer_func_camb_new ();

extern gchar *camb_filename;

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_CAMB_H_ */

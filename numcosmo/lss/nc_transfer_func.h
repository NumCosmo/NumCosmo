/***************************************************************************
 *            nc_transfer_func.h
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

#ifndef _NC_TRANSFER_FUNC_H_
#define _NC_TRANSFER_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_hireion.h>
#include <numcosmo/math/ncm_model_ctrl.h>

G_BEGIN_DECLS

#define NC_TYPE_TRANSFER_FUNC             (nc_transfer_func_get_type ())
#define NC_TRANSFER_FUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_TRANSFER_FUNC, NcTransferFunc))
#define NC_TRANSFER_FUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_TRANSFER_FUNC, NcTransferFuncClass))
#define NC_IS_TRANSFER_FUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_TRANSFER_FUNC))
#define NC_IS_TRANSFER_FUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_TRANSFER_FUNC))
#define NC_TRANSFER_FUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_TRANSFER_FUNC, NcTransferFuncClass))

typedef struct _NcTransferFuncClass NcTransferFuncClass;
typedef struct _NcTransferFunc NcTransferFunc;

struct _NcTransferFuncClass
{
  /*< private > */
  GObjectClass parent_class;
  
  gpointer (*alloc) (void);
  void (*prepare) (NcTransferFunc *tf, NcHICosmo *cosmo);
  gdouble (*calc) (NcTransferFunc *tf, gdouble k);
};

struct _NcTransferFunc
{
  /*< private > */
  GObject parent_instance;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_reion;
};

GType nc_transfer_func_get_type (void) G_GNUC_CONST;

NcTransferFunc *nc_transfer_func_new_from_name (gchar *transfer_name);
NcTransferFunc *nc_transfer_func_ref (NcTransferFunc *tf);

void nc_transfer_func_free (NcTransferFunc *tf);
void nc_transfer_func_clear (NcTransferFunc **tf);

void nc_transfer_func_prepare (NcTransferFunc *tf, NcHICosmo *cosmo);
void nc_transfer_func_prepare_if_needed (NcTransferFunc *tf, NcHICosmo *cosmo);

gdouble nc_transfer_func_eval (NcTransferFunc *tf, NcHICosmo *cosmo, gdouble kh);

G_END_DECLS

#endif /* _NC_TRANSFER_FUNC_H_ */


/***************************************************************************
 *            nc_powspec_ml_transfer.h
 *
 *  Thu March 17 14:57:27 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_powspec_ml_transfer.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_POWSPEC_ML_TRANSFER_H_
#define _NC_POWSPEC_ML_TRANSFER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_powspec_ml.h>
#include <numcosmo/lss/nc_transfer_func.h>
#include <numcosmo/lss/nc_growth_func.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_ML_TRANSFER             (nc_powspec_ml_transfer_get_type ())
#define NC_POWSPEC_ML_TRANSFER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_POWSPEC_ML_TRANSFER, NcPowspecMLTransfer))
#define NC_POWSPEC_ML_TRANSFER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_POWSPEC_ML_TRANSFER, NcPowspecMLTransferClass))
#define NC_IS_POWSPEC_ML_TRANSFER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_POWSPEC_ML_TRANSFER))
#define NC_IS_POWSPEC_ML_TRANSFER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_POWSPEC_ML_TRANSFER))
#define NC_POWSPEC_ML_TRANSFER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_POWSPEC_ML_TRANSFER, NcPowspecMLTransferClass))

typedef struct _NcPowspecMLTransferClass NcPowspecMLTransferClass;
typedef struct _NcPowspecMLTransfer NcPowspecMLTransfer;

struct _NcPowspecMLTransferClass
{
  /*< private > */
  NcPowspecMLClass parent_class;
};

struct _NcPowspecMLTransfer
{
  /*< private > */
  NcPowspecML parent_instance;
  NcTransferFunc *tf;
  NcGrowthFunc *gf;
  gdouble Pm_k2Pzeta;
};

GType nc_powspec_ml_transfer_get_type (void) G_GNUC_CONST;

NcPowspecMLTransfer *nc_powspec_ml_transfer_new (NcTransferFunc *tf);

void nc_powspec_ml_transfer_set_tf (NcPowspecMLTransfer *ps_mlt, NcTransferFunc *tf);
void nc_powspec_ml_transfer_set_gf (NcPowspecMLTransfer *ps_mlt, NcGrowthFunc *gf);

NcTransferFunc *nc_powspec_ml_transfer_peek_tf (NcPowspecMLTransfer *ps_mlt);
NcGrowthFunc *nc_powspec_ml_transfer_peek_gf (NcPowspecMLTransfer *ps_mlt);

G_END_DECLS

#endif /* _NC_POWSPEC_ML_TRANSFER_H_ */

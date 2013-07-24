/***************************************************************************
 *            nc_data_bao_dvdv.h
 *
 *  Thu Apr 22 15:31:19 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_BAO_DVDV_H_
#define _NC_DATA_BAO_DVDV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_data_bao.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_BAO_DVDV             (nc_data_bao_dvdv_get_type ())
#define NC_DATA_BAO_DVDV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_BAO_DVDV, NcDataBaoDVDV))
#define NC_DATA_BAO_DVDV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_BAO_DVDV, NcDataBaoDVDVClass))
#define NC_IS_DATA_BAO_DVDV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_BAO_DVDV))
#define NC_IS_DATA_BAO_DVDV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_BAO_DVDV))
#define NC_DATA_BAO_DVDV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_BAO_DVDV, NcDataBaoDVDVClass))

typedef struct _NcDataBaoDVDVClass NcDataBaoDVDVClass;
typedef struct _NcDataBaoDVDV NcDataBaoDVDV;

struct _NcDataBaoDVDV
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcDistance *dist;
  NcDataBaoId id;
};

struct _NcDataBaoDVDVClass
{
  /*< private >*/
  NcmDataGaussDiagClass parent_class;
};

GType nc_data_bao_dvdv_get_type (void) G_GNUC_CONST;

NcmData *nc_data_bao_dvdv_new (NcDistance *dist, NcDataBaoId id);

void nc_data_bao_dvdv_set_sample (NcDataBaoDVDV *bao_dvdv, NcDataBaoId id);
NcDataBaoId nc_data_bao_dvdv_get_sample (NcDataBaoDVDV *bao_dvdv);

G_END_DECLS

#endif /* _NC_DATA_BAO_DVDV_H_ */


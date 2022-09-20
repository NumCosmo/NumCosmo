/***************************************************************************
 *            nc_data_bao_dtr_dhr.h
 *
 *  Mon Ago 15 14:38:40 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_bao_dtr_dhr.h
 * Copyright (C) 2022 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_BAO_DTR_DHR_H_
#define _NC_DATA_BAO_DTR_DHR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/data/nc_data_bao.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_BAO_DTR_DHR             (nc_data_bao_dtr_dhr_get_type ())
#define NC_DATA_BAO_DTR_DHR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_BAO_DTR_DHR, NcDataBaoDtrDHr))
#define NC_DATA_BAO_DTR_DHR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_BAO_DTR_DHR, NcDataBaoDtrDHrClass))
#define NC_IS_DATA_BAO_DTR_DHR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_BAO_DTR_DHR))
#define NC_IS_DATA_BAO_DTR_DHR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_BAO_DTR_DHR))
#define NC_DATA_BAO_DTR_DHR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_BAO_DTR_DHR, NcDataBaoDtrDHrClass))

typedef struct _NcDataBaoDtrDHrClass NcDataBaoDtrDHrClass;
typedef struct _NcDataBaoDtrDHr NcDataBaoDtrDHr;

struct _NcDataBaoDtrDHrClass
{
  /*< private >*/
  NcmDataGaussCovClass parent_class;
};

struct _NcDataBaoDtrDHr
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcDistance *dist;
  NcmVector *x;
};

GType nc_data_bao_dtr_dhr_get_type (void) G_GNUC_CONST;

NcDataBaoDtrDHr *nc_data_bao_dtr_dhr_new_from_file (const gchar *filename);
NcDataBaoDtrDHr *nc_data_bao_dtr_dhr_new_from_id (NcDistance *dist, NcDataBaoId id);
void nc_data_bao_dtr_dhr_set_dist (NcDataBaoDtrDHr *dhdt, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_DTR_DHR_H_ */

/***************************************************************************
 *            nc_data_bao_dtr_dhr.h
 *
 *  Mon Ago 15 14:38:40 2022
 *  Copyright  2022  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_bao_dtr_dhr.h
 * Copyright (C) 2022 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#define NC_TYPE_DATA_BAO_DTR_DHR (nc_data_bao_dtr_dhr_get_type ())

G_DECLARE_FINAL_TYPE (NcDataBaoDtrDHr, nc_data_bao_dtr_dhr, NC, DATA_BAO_DTR_DHR, NcmDataGaussCov)

NcDataBaoDtrDHr *nc_data_bao_dtr_dhr_new_from_file (const gchar *filename);
NcDataBaoDtrDHr *nc_data_bao_dtr_dhr_new_from_id (NcDistance *dist, NcDataBaoId id);
void nc_data_bao_dtr_dhr_set_dist (NcDataBaoDtrDHr *dhdt, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_DTR_DHR_H_ */


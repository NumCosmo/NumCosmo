/***************************************************************************
 *            nc_data_bao_dv.h
 *
 *  Thu Apr 22 15:31:19 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_BAO_DV_H_
#define _NC_DATA_BAO_DV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/data/nc_data_bao.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_BAO_DV             (nc_data_bao_dv_get_type ())

G_DECLARE_FINAL_TYPE (NcDataBaoDV, nc_data_bao_dv, NC, DATA_BAO_DV, NcmDataGaussDiag)

NcDataBaoDV *nc_data_bao_dv_new_from_file (const gchar *filename);
NcDataBaoDV *nc_data_bao_dv_new_from_id (NcDistance *dist, NcDataBaoId id);
void nc_data_bao_dv_set_dist (NcDataBaoDV *bao_dv, NcDistance *dist);

G_END_DECLS

#endif /* _NC_DATA_BAO_DV_H_ */


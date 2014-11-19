/***************************************************************************
 *            nc_data_bao.h
 *
 *  Thu November 22 20:41:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

#ifndef _NC_DATA_BAO_H_
#define _NC_DATA_BAO_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataBaoId:
 * @NC_DATA_BAO_A_EISENSTEIN2005: FIXME
 * @NC_DATA_BAO_DV_EISENSTEIN2005: FIXME
 * @NC_DATA_BAO_DVDV_PERCIVAL2007: FIXME
 * @NC_DATA_BAO_DVDV_PERCIVAL2010: FIXME
 * @NC_DATA_BAO_RDV_PERCIVAL2007: FIXME
 * @NC_DATA_BAO_RDV_PERCIVAL2010: FIXME
 * @NC_DATA_BAO_RDV_BEUTLER2011: FIXME
 * @NC_DATA_BAO_RDV_PADMANABHAN2012: FIXME
 * @NC_DATA_BAO_RDV_ANDERSON2012: FIXME
 * @NC_DATA_BAO_RDV_BLAKE2012: FIXME
 *
 * FIXME
 */
typedef enum _NcDataBaoId
{
  NC_DATA_BAO_A_EISENSTEIN2005 = 0,
  NC_DATA_BAO_DV_EISENSTEIN2005,
  NC_DATA_BAO_DVDV_PERCIVAL2007,
  NC_DATA_BAO_DVDV_PERCIVAL2010,
  NC_DATA_BAO_RDV_PERCIVAL2007,  
  NC_DATA_BAO_RDV_PERCIVAL2010,
  NC_DATA_BAO_RDV_BEUTLER2011,
  NC_DATA_BAO_RDV_PADMANABHAN2012,
  NC_DATA_BAO_RDV_ANDERSON2012,
  NC_DATA_BAO_RDV_BLAKE2012,    /*< private >*/
  NC_DATA_BAO_NSAMPLES,         /*< skip >*/
} NcDataBaoId;

NcmData *nc_data_bao_create (NcDistance *dist, NcDataBaoId id);

G_END_DECLS

#endif /* _NC_DATA_BAO_H_ */

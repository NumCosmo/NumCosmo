/***************************************************************************
 *            nc_data_bao.c
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

/**
 * SECTION:nc_data_bao
 * @title: BAO Data
 * @short_description: Object representing BAO data
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_bao.h"

#include "nc_data_bao_a.h"
#include "nc_data_bao_dv.h"
#include "nc_data_bao_rdv.h"
#include "nc_data_bao_dvdv.h"

/**
 * nc_data_bao_new:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_bao_new (NcDistance *dist, NcDataBaoId id)
{
  switch (id)
  {
    case NC_DATA_BAO_A_EISENSTEIN:
      return nc_data_bao_a_new (dist, id);
      break;      
    case NC_DATA_BAO_DV_EISENSTEIN:
      return nc_data_bao_dv_new (dist, id);
      break;
    case NC_DATA_BAO_DVDV_PERCIVAL:
      return nc_data_bao_dvdv_new (dist, id);
      break;
    case NC_DATA_BAO_RDV_PERCIVAL:
      return nc_data_bao_rdv_new (dist, id);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}


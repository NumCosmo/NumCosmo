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
 * @title: NcDataBao
 * @short_description: Helper functions for instantiating BAO data.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao.h"

#include "data/nc_data_bao_a.h"
#include "data/nc_data_bao_dv.h"
#include "data/nc_data_bao_rdv.h"
#include "data/nc_data_bao_dvdv.h"
#include "data/nc_data_bao_empirical_fit.h"

/**
 * nc_data_bao_create:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * FIXME
 *
 * Returns: (transfer full): a #NcmData
 */
NcmData *
nc_data_bao_create (NcDistance *dist, NcDataBaoId id)
{
  /* FIXME the switch is wrong but works since we have only two options. */
  switch (id)
  {
    case NC_DATA_BAO_A_EISENSTEIN2005:
      return nc_data_bao_a_new (dist, id);
      break;  
    case NC_DATA_BAO_DV_EISENSTEIN2005:
      return nc_data_bao_dv_new (dist, id);
      break;
    case NC_DATA_BAO_DVDV_PERCIVAL2007:
    case NC_DATA_BAO_DVDV_PERCIVAL2010:
      return nc_data_bao_dvdv_new (dist, id);
      break;
    case NC_DATA_BAO_RDV_PERCIVAL2007:
    case NC_DATA_BAO_RDV_PERCIVAL2010:
    case NC_DATA_BAO_RDV_BEUTLER2011:
    case NC_DATA_BAO_RDV_PADMANABHAN2012:
    case NC_DATA_BAO_RDV_ANDERSON2012:
    case NC_DATA_BAO_RDV_BLAKE2012:
    case NC_DATA_BAO_RDV_KAZIN2014:  
      return nc_data_bao_rdv_new (dist, id);
      break;
    case NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015:  
      return NCM_DATA (nc_data_bao_empirical_fit_new_from_id (dist, id));
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}


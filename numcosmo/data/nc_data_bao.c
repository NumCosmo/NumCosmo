/***************************************************************************
 *            nc_data_bao.c
 *
 *  Thu November 22 20:41:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

/**
 * NcDataBao:
 *
 * Helper functions for instantiating BAO data.
 *
 * A set of factory functions to instantiate BAO data.
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
#include "data/nc_data_bao_empirical_fit_2d.h"
#include "data/nc_data_bao_dhr_dar.h"
#include "data/nc_data_bao_dmr_hr.h"
#include "data/nc_data_bao_dtr_dhr.h"
#include "data/nc_data_bao_dvr_dtdh.h"

/**
 * nc_data_bao_create:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * Creates a #NcmData for the given @id.
 *
 * Returns: (transfer full): a #NcmData
 */
NcmData *
nc_data_bao_create (NcDistance *dist, NcDataBaoId id)
{
  NcmData *data;

  switch (id)
  {
    case NC_DATA_BAO_A_EISENSTEIN2005:
      data = NCM_DATA (nc_data_bao_a_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_DV_EISENSTEIN2005:
      data = NCM_DATA (nc_data_bao_dv_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_DVDV_PERCIVAL2007:
    case NC_DATA_BAO_DVDV_PERCIVAL2010:
      data = NCM_DATA (nc_data_bao_dvdv_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_RDV_PERCIVAL2007:
    case NC_DATA_BAO_RDV_PERCIVAL2010:
    case NC_DATA_BAO_RDV_BEUTLER2011:
    case NC_DATA_BAO_RDV_PADMANABHAN2012:
    case NC_DATA_BAO_RDV_ANDERSON2012:
    case NC_DATA_BAO_RDV_BLAKE2012:
    case NC_DATA_BAO_RDV_KAZIN2014:
    case NC_DATA_BAO_RDV_BOSS_QSO_ATA2017:
    case NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024:
    case NC_DATA_BAO_RDV_DESI_DR2_BGS_2025:
      data = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015:
    case NC_DATA_BAO_EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021:
      data = NCM_DATA (nc_data_bao_empirical_fit_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_EMPIRICAL_FIT_2D_BAUTISTA2017:
    case NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021:
    case NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021:
      data = NCM_DATA (nc_data_bao_empirical_fit_2d_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015:
    case NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015_LYAF_AUTO_CROSS:
      data = NCM_DATA (nc_data_bao_dhr_dar_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_DMR_HR_SDSS_DR12_2016:
      data = NCM_DATA (nc_data_bao_dmr_hr_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE:
    case NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021:
    case NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021:
    case NC_DATA_BAO_DTR_DHR_DESI_DR1_LYA_2025:
      data = NCM_DATA (nc_data_bao_dtr_dhr_new_from_id (dist, id));
      break;
    case NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024:
    case NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025:
      data = NCM_DATA (nc_data_bao_dvr_dtdh_new_from_id (dist, id));
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  g_assert (NCM_IS_DATA (data));

  return data;
}


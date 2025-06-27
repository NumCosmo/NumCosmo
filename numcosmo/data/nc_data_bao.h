/***************************************************************************
 *            nc_data_bao.h
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
 * @NC_DATA_BAO_A_EISENSTEIN2005: [Eisenstein et al. (2005)][XEisenstein2005]
 * @NC_DATA_BAO_DV_EISENSTEIN2005: [Eisenstein et al. (2005)][XEisenstein2005]
 * @NC_DATA_BAO_DVDV_PERCIVAL2007: [Percival et al. (2007)][XPercival2007]
 * @NC_DATA_BAO_DVDV_PERCIVAL2010: [Percival et al. (2010)][XPercival2010]
 * @NC_DATA_BAO_RDV_PERCIVAL2007: [Percival et al. (2007)][XPercival2007]
 * @NC_DATA_BAO_RDV_PERCIVAL2010: [Percival et al. (2010)][XPercival2010]
 * @NC_DATA_BAO_RDV_BEUTLER2011: [Beutler et al. (2011)][XBeutler2011]
 * @NC_DATA_BAO_RDV_PADMANABHAN2012: [Padmanabhan et al. (2012)][XPadmanabhan2012]
 * @NC_DATA_BAO_RDV_ANDERSON2012: [Anderson et al. (2012)][XAnderson2012]
 * @NC_DATA_BAO_RDV_BLAKE2012: [Blake et al. (2011)][XBlake2011]
 * @NC_DATA_BAO_RDV_KAZIN2014: [Kazin et al. (2014)][XKazin2014]
 * @NC_DATA_BAO_RDV_BOSS_QSO_ATA2017: [Ata et al. (2017)][XAta2017]
 * @NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024: Adame et al., arXiv:2404.03000, table 18
 * @NC_DATA_BAO_RDV_DESI_DR2_BGS_2025: Abdul Karim et al., arXiv:2503.14738, table IV
 * @NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015: [Ross et al. (2015)][XRoss2015]
 * @NC_DATA_BAO_EMPIRICAL_FIT_2D_BAUTISTA2017: [Bautista et al. (2017)][XBautista2017]
 * @NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015: [Delubac et al. (2015)][XDelubac2015]
 * @NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015_LYAF_AUTO_CROSS: [Aubourg et al. (2014)][XAubourg2014]
 * @NC_DATA_BAO_DMR_HR_SDSS_DR12_2016: [Alam et al. (2016)][XAlam2016]
 * @NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE: [Alam et al. (2016)][XAlam2016]
 * @NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021: [Alam et al. (2016)][XAlam2016]
 * @NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021: [Alam et al. (2016)][XAlam2016]
 * @NC_DATA_BAO_DTR_DHR_DESI_DR1_LYA_2025: Adame at el., JCAP 01 (2025) 124, arXiv:2404.03001
 * @NC_DATA_BAO_EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021: [Alam et al. (2021)][XAlam2021]
 * @NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021: [Alam et al. (2021)][XAlam2021]
 * @NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021: [Alam et al. (2021)][XAlam2021]
 * @NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024: Adame at el., arXiv:2404.03000, table 18
 * @NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025: Abdul Karim et al., arXiv:2503.14738, table IV
 *
 * FIXME
 *
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
  NC_DATA_BAO_RDV_BLAKE2012,
  NC_DATA_BAO_RDV_KAZIN2014,
  NC_DATA_BAO_RDV_BOSS_QSO_ATA2017,
  NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024,
  NC_DATA_BAO_RDV_DESI_DR2_BGS_2025,
  NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015,
  NC_DATA_BAO_EMPIRICAL_FIT_2D_BAUTISTA2017,
  NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015,
  NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015_LYAF_AUTO_CROSS,
  NC_DATA_BAO_DMR_HR_SDSS_DR12_2016,
  NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE,
  NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021,
  NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021,
  NC_DATA_BAO_DTR_DHR_DESI_DR1_LYA_2025,
  NC_DATA_BAO_EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021,
  NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021,
  NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021,
  NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024,
  NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025,
  /* < private > */
  NC_DATA_BAO_NSAMPLES, /*< skip >*/
} NcDataBaoId;

#define NC_DATA_BAO_RDV_FIRST NC_DATA_BAO_RDV_PERCIVAL2007
#define NC_DATA_BAO_RDV_LAST NC_DATA_BAO_RDV_DESI_DR2_BGS_2025
#define NC_DATA_BAO_RDV_LEN (NC_DATA_BAO_RDV_LAST - NC_DATA_BAO_RDV_FIRST + 1)

NcmData *nc_data_bao_create (NcDistance *dist, NcDataBaoId id);

G_END_DECLS

#endif /* _NC_DATA_BAO_H_ */


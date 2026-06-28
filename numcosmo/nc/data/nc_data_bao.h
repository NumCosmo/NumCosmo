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
#include <numcosmo/ncm/data/ncm_data.h>
#include <numcosmo/nc/background/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataBaoId:
 * @NC_DATA_BAO_A_EISENSTEIN2005: [Eisenstein et al. (2005)](https://arxiv.org/abs/astro-ph/0501171)
 * @NC_DATA_BAO_DV_EISENSTEIN2005: [Eisenstein et al. (2005)](https://arxiv.org/abs/astro-ph/0501171)
 * @NC_DATA_BAO_DVDV_PERCIVAL2007: [Percival et al. (2007)](https://arxiv.org/abs/0705.3323)
 * @NC_DATA_BAO_DVDV_PERCIVAL2010: [Percival et al. (2010)](https://arxiv.org/abs/0907.1660)
 * @NC_DATA_BAO_RDV_PERCIVAL2007: [Percival et al. (2007)](https://arxiv.org/abs/0705.3323)
 * @NC_DATA_BAO_RDV_PERCIVAL2010: [Percival et al. (2010)](https://arxiv.org/abs/0907.1660)
 * @NC_DATA_BAO_RDV_BEUTLER2011: [Beutler et al. (2011)](https://arxiv.org/abs/1106.3366)
 * @NC_DATA_BAO_RDV_PADMANABHAN2012: [Padmanabhan et al. (2012)](https://arxiv.org/abs/1202.0090)
 * @NC_DATA_BAO_RDV_ANDERSON2012: [Anderson et al. (2012)](https://arxiv.org/abs/1203.6594)
 * @NC_DATA_BAO_RDV_BLAKE2012: [Blake et al. (2011)](https://arxiv.org/abs/1108.2635)
 * @NC_DATA_BAO_RDV_KAZIN2014: [Kazin et al. (2014)](https://arxiv.org/abs/1401.0358)
 * @NC_DATA_BAO_RDV_BOSS_QSO_ATA2017: [Ata et al. (2017)](https://arxiv.org/abs/1705.06373)
 * @NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024: Adame et al., arXiv:2404.03000, table 18
 * @NC_DATA_BAO_RDV_DESI_DR2_BGS_2025: Abdul Karim et al., arXiv:2503.14738, table IV
 * @NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015: [Ross et al. (2015)](https://arxiv.org/abs/1409.3242)
 * @NC_DATA_BAO_EMPIRICAL_FIT_2D_BAUTISTA2017: [Bautista et al. (2017)](https://arxiv.org/abs/1702.00176)
 * @NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015: [Delubac et al. (2015)](https://arxiv.org/abs/1404.1801)
 * @NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015_LYAF_AUTO_CROSS: [Aubourg et al. (2014)](https://arxiv.org/abs/1411.1074)
 * @NC_DATA_BAO_DMR_HR_SDSS_DR12_2016: [Alam et al. (2016)](https://arxiv.org/abs/1607.03155)
 * @NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE: [Alam et al. (2016)](https://arxiv.org/abs/1607.03155)
 * @NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021: [Alam et al. (2016)](https://arxiv.org/abs/1607.03155)
 * @NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021: [Alam et al. (2016)](https://arxiv.org/abs/1607.03155)
 * @NC_DATA_BAO_DTR_DHR_DESI_DR1_LYA_2025: Adame at el., JCAP 01 (2025) 124, arXiv:2404.03001
 * @NC_DATA_BAO_EMPIRICAL_FIT_1D_SDSS_DR16_ELG_2021: [Alam et al. (2021)](https://arxiv.org/abs/2007.08991)
 * @NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021: [Alam et al. (2021)](https://arxiv.org/abs/2007.08991)
 * @NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021: [Alam et al. (2021)](https://arxiv.org/abs/2007.08991)
 * @NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024: Adame at el., arXiv:2404.03000, table 18
 * @NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025: Abdul Karim et al., arXiv:2503.14738, table IV
 *
 * Enumeration of available BAO (Baryon Acoustic Oscillations) datasets.
 *
 * These enumerate different Baryon Acoustic Oscillation observations from various
 * surveys and releases. Each identifier corresponds to a specific measurement type
 * and publication referenced in the parameter documentation above.
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


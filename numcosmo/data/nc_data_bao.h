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
 * @NC_DATA_BAO_A_EISENSTEIN2005: <link linkend="XEisenstein2005">Eisenstein et al. (2005)</link>
 * @NC_DATA_BAO_DV_EISENSTEIN2005: <link linkend="XEisenstein2005">Eisenstein et al. (2005)</link>
 * @NC_DATA_BAO_DVDV_PERCIVAL2007: <link linkend="XPercival2007">Percival et al. (2007)</link>
 * @NC_DATA_BAO_DVDV_PERCIVAL2010: <link linkend="XPercival2010">Percival et al. (2010)</link>
 * @NC_DATA_BAO_RDV_PERCIVAL2007: <link linkend="XPercival2007">Percival et al. (2007)</link>
 * @NC_DATA_BAO_RDV_PERCIVAL2010: <link linkend="XPercival2010">Percival et al. (2010)</link>
 * @NC_DATA_BAO_RDV_BEUTLER2011: <link linkend="XBeutler2011">Beutler et al. (2011)</link>
 * @NC_DATA_BAO_RDV_PADMANABHAN2012: <link linkend="XPadmanabhan2012">Padmanabhan et al. (2012)</link>
 * @NC_DATA_BAO_RDV_ANDERSON2012: <link linkend="XAnderson2012">Anderson et al. (2012)</link>
 * @NC_DATA_BAO_RDV_BLAKE2012: <link linkend="XBlake2011">Blake et al. (2011)</link>
 * @NC_DATA_BAO_RDV_KAZIN2014: <link linkend="XKazin2014">Kazin et al. (2014)</link>
 * @NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015: <link linkend="XRoss2014">Ross et al. (2015)</link>
 * @NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015: <link linkend="XDelubac2015">Delubac et al. (2015)</link>
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
  NC_DATA_BAO_RDV_BLAKE2012,
  NC_DATA_BAO_RDV_KAZIN2014,  
  NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015, 
  NC_DATA_BAO_DHR_DAR_SDSS_DR11_2015, /*< private >*/
  NC_DATA_BAO_NSAMPLES,               /*< skip >*/
} NcDataBaoId;

#define NC_DATA_BAO_RDV_FIRST NC_DATA_BAO_RDV_PERCIVAL2007
#define NC_DATA_BAO_RDV_LAST NC_DATA_BAO_RDV_KAZIN2014
#define NC_DATA_BAO_RDV_LEN (NC_DATA_BAO_RDV_LAST - NC_DATA_BAO_RDV_FIRST + 1)

NcmData *nc_data_bao_create (NcDistance *dist, NcDataBaoId id);

G_END_DECLS

#endif /* _NC_DATA_BAO_H_ */

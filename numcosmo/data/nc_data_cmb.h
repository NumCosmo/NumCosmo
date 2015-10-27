/***************************************************************************
 *            nc_data_cmb.h
 *
 *  Thu November 22 21:19:23 2012
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

#ifndef _NC_DATA_CMB_H_
#define _NC_DATA_CMB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataCMBDataType:
 * @NC_DATA_CMB_TYPE_TT: Temperature - Temperature correlation data
 * @NC_DATA_CMB_TYPE_EE: Mode E - Mode E correlation data
 * @NC_DATA_CMB_TYPE_BB: Mode B - Mode B correlation data
 * @NC_DATA_CMB_TYPE_TE: Temperature - Mode E correlation data
 * @NC_DATA_CMB_TYPE_TB: Temperature - Mode B correlation data
 * @NC_DATA_CMB_TYPE_EB: Mode E - Mode B correlation data
 * @NC_DATA_CMB_TYPE_PHIPHI: $\phi$ - $\phi$ correlation data
 * @NC_DATA_CMB_TYPE_ALL: All types above
 *
 * FIXME
 *
 */
typedef enum _NcDataCMBDataType
{
  NC_DATA_CMB_TYPE_TT     = 1 << 0,
  NC_DATA_CMB_TYPE_EE     = 1 << 1,
  NC_DATA_CMB_TYPE_BB     = 1 << 2,
  NC_DATA_CMB_TYPE_TE     = 1 << 3,
  NC_DATA_CMB_TYPE_TB     = 1 << 4,
  NC_DATA_CMB_TYPE_EB     = 1 << 5,
  NC_DATA_CMB_TYPE_PHIPHI = 1 << 6,
  NC_DATA_CMB_TYPE_ALL    = (1 << 7) - 1, /*< private >*/
  NC_DATA_CMB_TYPE_LEN,                   /*< skip >*/
} NcDataCMBDataType;

/**
 * NcDataCMBId:
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP3: FIXME
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP5: FIXME
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP7: FIXME
 * @NC_DATA_CMB_DIST_PRIORS_WMAP5: FIXME
 * @NC_DATA_CMB_DIST_PRIORS_WMAP7: FIXME
 * @NC_DATA_CMB_DIST_PRIORS_WMAP9: FIXME
 *
 * FIXME
 */
typedef enum _NcDataCMBId
{
  NC_DATA_CMB_SHIFT_PARAM_WMAP3 = 0,
  NC_DATA_CMB_SHIFT_PARAM_WMAP5,
  NC_DATA_CMB_SHIFT_PARAM_WMAP7,
  NC_DATA_CMB_DIST_PRIORS_WMAP5,
  NC_DATA_CMB_DIST_PRIORS_WMAP7,
  NC_DATA_CMB_DIST_PRIORS_WMAP9,     /*< private >*/
  NC_DATA_CMB_NSAMPLES,  /*< skip >*/
} NcDataCMBId;

NcmData *nc_data_cmb_create (NcDistance *dist, NcDataCMBId id);

G_END_DECLS

#endif /* _NC_DATA_CMB_H_ */


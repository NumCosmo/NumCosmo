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
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataCMBId:
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP3: FIXME
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP5: FIXME
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP7: FIXME
 * @NC_DATA_CMB_DIST_PRIORS_WMAP5: FIXME
 * @NC_DATA_CMB_DIST_PRIORS_WMAP7: FIXME
 *
 * FIXME
 */
typedef enum _NcDataCMBId
{
  NC_DATA_CMB_SHIFT_PARAM_WMAP3 = 0,
  NC_DATA_CMB_SHIFT_PARAM_WMAP5,
  NC_DATA_CMB_SHIFT_PARAM_WMAP7,
  NC_DATA_CMB_DIST_PRIORS_WMAP5,
  NC_DATA_CMB_DIST_PRIORS_WMAP7,     /*< private >*/
  NC_DATA_CMB_NSAMPLES,  /*< skip >*/
} NcDataCMBId;

NcmData *nc_data_cmb_create (NcDistance *dist, NcDataCMBId id);

G_END_DECLS

#endif /* _NC_DATA_CMB_H_ */


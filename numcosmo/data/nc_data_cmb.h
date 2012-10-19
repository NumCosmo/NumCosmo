/***************************************************************************
 *            nc_data_cmb.h
 *
 *  Thu Apr 22 15:56:59 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NC_DATA_CMB_H
#define _NC_DATA_CMB_H

#include <glib.h>

G_BEGIN_DECLS

/**
 * NcDataCMBId:
 * @NC_DATA_CMB_SHIFT_PARAMETER_WMAP3: FIXME
 * @NC_DATA_CMB_SHIFT_PARAMETER_WMAP5: FIXME
 * @NC_DATA_CMB_DISTANCE_PRIORS_WMAP5: FIXME
 *
 * FIXME
 */
typedef enum _NcDataCMBId
{
  NC_DATA_CMB_SHIFT_PARAMETER_WMAP3 = 0,
  NC_DATA_CMB_SHIFT_PARAMETER_WMAP5,
  NC_DATA_CMB_DISTANCE_PRIORS_WMAP5, /*< private >*/
  NC_DATA_CMB_NSAMPLES,              /*< skip >*/
} NcDataCMBId;

NcData *nc_data_cmb (NcDistance *dist, NcDataCMBId cmb_id);

G_END_DECLS

#endif /* _NC_DATA_CMB_H */

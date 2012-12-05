/***************************************************************************
 *            nc_data_snia.h
 *
 *  Mon December 10 00:20:48 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_SNIA_H_
#define _NC_DATA_SNIA_H_

#include <glib.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_snia_dist_cov.h>

G_BEGIN_DECLS

/**
 * NcDataSNIATypeId:
 * @NC_DATA_SNIA_TYPE_SIMPLE: FIXME
 * @NC_DATA_SNIA_TYPE_COV: FIXME
 *
 * FIXME
 */
typedef enum _NcDataSNIAType
{
  NC_DATA_SNIA_TYPE_SIMPLE = 0,
  NC_DATA_SNIA_TYPE_COV,        /*< private >*/
  NC_DATA_SNIA_TYPE_LEN,        /*< skip >*/
} NcDataSNIAType;

/**
 * NcDataSNIACatId:
 * @NC_DATA_SNIA_CAT_SNLS3_SYS_STAT: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcDataSNIACatId
{
  NC_DATA_SNIA_CAT_SNLS3_SYS_STAT = 0, 
  NC_DATA_SNIA_CAT_SNLS3_STAT_ONLY,    /*< private >*/
  NC_DATA_SNIA_CAT_LEN,                /*< skip >*/
} NcDataSNIACatId;

void nc_data_snia_load_cat (NcSNIADistCov *dcov, NcDataSNIACatId id);
gchar *nc_data_snia_get_fits (const gchar *filename, gboolean check_size);

G_END_DECLS

#endif /* _NC_DATA_SNIA_H_ */



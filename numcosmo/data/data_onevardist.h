/***************************************************************************
 *            data_onevardist.h
 *
 *  Thu Apr 15 11:17:25 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_DATA_ONEVARDIST_H
#define _NC_DATA_ONEVARDIST_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/data/data.h>

G_BEGIN_DECLS

typedef struct _NcDataOneVarDist NcDataOneVarDist;
typedef struct _NcDataModelOneVarDist NcDataModelOneVarDist;

/**
 * NcDataOneVarDist:
 *
 * FIXME
 */
struct _NcDataOneVarDist
{
  /*< private >*/
  guint np;
  NcmVector *x;
  NcDataStruct *extra_data;
};

/**
 * NcDataModelOneVarDist:
 *
 * FIXME
 */
struct _NcDataModelOneVarDist
{
  /*< private >*/
  NcmMSetFunc *dist;
  NcmMSetFunc *inv_pdf;
};

NcData *nc_data_onevardist_new (NcmMSetFunc *dist, NcmMSetFunc *inv_pdf);
void nc_data_onevardist_set_prepare (NcData *data, NcDataPrepare prepare);
void nc_data_onevardist_init_from_vector (NcData *data, NcmVector *x, NcDataStruct *extra_data);

G_END_DECLS

#endif /* NC_DATA_ONEVARDIST_H */

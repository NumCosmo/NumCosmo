/***************************************************************************
 *            nc_data_hubble.h
 *
 *  Thu Apr 22 14:35:37 2010
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

#ifndef _NC_DATA_HUBBLE_H
#define _NC_DATA_HUBBLE_H

#include <glib.h>
#include <glib-object.h>

#include <numcosmo/data/data.h>

G_BEGIN_DECLS

/**
 * NcDataHubbleId:
 * @NC_DATA_HUBBLE_VERDE: FIXME
 * @NC_DATA_HUBBLE_CABRE: FIXME
 *
 * FIXME
 */
typedef enum _NcDataHubbleId
{
  NC_DATA_HUBBLE_VERDE = 0,
  NC_DATA_HUBBLE_CABRE,    /*< private >*/
  NC_DATA_HUBBLE_NSAMPLES, /*< skip >*/
} NcDataHubbleId;

NcData *nc_data_hubble (NcDataHubbleId H_id);

G_END_DECLS

#endif /* _NC_DATA_HUBBLE_FUNCTION_H */

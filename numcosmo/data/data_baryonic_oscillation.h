/***************************************************************************
 *            data_baryonic_oscillation.h
 *
 *  Thu Apr 22 15:31:19 2010
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

#ifndef _NC_DATA_BARYONIC_OSCILLATION_H
#define _NC_DATA_BARYONIC_OSCILLATION_H

#include <glib.h>
#include <gsl/gsl_eigen.h>

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

G_BEGIN_DECLS

/**
 * NcDataBaoSampleId:
 * @NC_DATA_BAO_A_SAMPLE_EISENSTEIN: FIXME
 * @NC_DATA_BAO_Dv_SAMPLE_EISENSTEIN: FIXME
 * @NC_DATA_BAO_R_DV_SAMPLE_PERCIVAL: FIXME
 * @NC_DATA_BAO_DV_DV_SAMPLE_PERCIVAL: FIXME
 * 
 * FIXME
 */
typedef enum _NcDataBaoSampleId
{
  NC_DATA_BAO_A_SAMPLE_EISENSTEIN = 0,
  NC_DATA_BAO_Dv_SAMPLE_EISENSTEIN,
  NC_DATA_BAO_R_DV_SAMPLE_PERCIVAL,
  NC_DATA_BAO_DV_DV_SAMPLE_PERCIVAL
} NcDataBaoSampleId;

#define NC_DATA_BAO_NSAMPLES (NC_DATA_BAO_DV_DV_SAMPLE_PERCIVAL + 1)

NcData *nc_data_bao (NcDistance *dist, NcDataBaoSampleId bao_id);

G_END_DECLS

#endif /* _NC_DATA_BARYONIC_OSCILLATION_H */

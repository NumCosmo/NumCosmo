/***************************************************************************
 *            data_distance_modulus.h
 *
 *  Thu Apr 22 10:37:39 2010
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

/**
 * SECTION:data_distance_modulus
 * @title: Distance Modulus Data
 * @short_description: Data samples of distance modulus
 *
 * FIXME
 * 
 */

#ifndef _NC_DATA_DISTANCE_MODULUS_H
#define _NC_DATA_DISTANCE_MODULUS_H

#include <glib.h>
#include <gsl/gsl_eigen.h>

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

G_BEGIN_DECLS

/**
 * NcDataDistanceModulusSNIaSampleId:
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_GOLD_157: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_GOLD_182: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_GOLD_182_FULL: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_ESSENCE: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_LEGACY: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_UNION: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_CfA3: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_UNION2: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_UNION2_1: FIXME
 * @NC_DATA_DISTANCE_MODULUS_SNIA_NSAMPLES: FIXME
 */  
typedef enum _NcDataDistanceModulusSNIaSampleId
{
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_GOLD_157 = 0,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_GOLD_182,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_GOLD_182_FULL,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_ESSENCE,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_LEGACY,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_UNION,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_CfA3,
  NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_UNION2,
	NC_DATA_DISTANCE_MODULUS_SNIA_SAMPLE_UNION2_1,
  NC_DATA_DISTANCE_MODULUS_SNIA_NSAMPLES
} NcDataDistanceModulusSNIaSampleId;

NcData *nc_data_distance_modulus_snia (NcDistance *dist, NcDataDistanceModulusSNIaSampleId snia_id);

G_END_DECLS

#endif /* _NC_DATA_DISTANCE_MODULUS_H */

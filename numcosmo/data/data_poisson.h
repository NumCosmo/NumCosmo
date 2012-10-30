/***************************************************************************
 *            data_poisson.h
 *
 *  Sun Apr  4 21:57:50 2010
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

#ifndef _NC_DATA_POISSON_H
#define _NC_DATA_POISSON_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/data/data.h>
#include <gsl/gsl_histogram.h>

G_BEGIN_DECLS

typedef struct _NcDataPoisson NcDataPoisson;

/**
 * NcDataPoisson:
 *
 * FIXME
 */
struct _NcDataPoisson
{
  /* < private > */
  gsl_histogram *h;
  NcDataStruct *extra_data;
};

/**
 * NcDataPoissonType:
 * @NC_DATA_POISSON_INT: FIXME
 *
 * FIXME
 */
typedef enum _NcDataPoissonType
{
  NC_DATA_POISSON_INT,
} NcDataPoissonType;

NcData *nc_data_poisson_new (NcDataPoissonType poisson_type);

void nc_data_poisson_set_function (NcData *data, NcmMSetFunc *func);
void nc_data_poisson_set_prepare (NcData *data, NcDataPrepare prepare);
void nc_data_poisson_set_resample (NcData *data, NcDataResample resample);
void nc_data_poisson_init_from_vector (NcData *data, NcmVector *nodes, gsl_vector_ulong *N, NcDataStruct *extra_data);
void nc_data_poisson_init_from_histogram (NcData *data, gsl_histogram *h, gboolean steal, NcDataStruct *extra_data);
void nc_data_poisson_init_zero (NcData *data, NcmVector *nodes, NcDataStruct *extra_data);

G_END_DECLS

#endif /* _NC_DATA_POISSON_H */

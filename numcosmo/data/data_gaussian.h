/***************************************************************************
 *            data_gaussian.h
 *
 *  Fri Mar 19 14:57:54 2010
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

#ifndef _NC_DATA_GAUSSIAN_H
#define _NC_DATA_GAUSSIAN_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/data/data.h>
#include <gsl/gsl_eigen.h>

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

G_BEGIN_DECLS

typedef struct _NcDataGaussian NcDataGaussian;

/**
 * NcDataGaussian:
 *
 * FIXME
 */
struct _NcDataGaussian
{
  /*< private >*/
  gdouble wt;
  guint np;
  gsl_eigen_symmv_workspace *vw;
  NcmVector *x;
  NcmVector *y;
  NcmVector *weight;
  NcmVector *sigma;
  NcmVector *vp;
  NcmMatrix *inv_cov;
  NcmMatrix *orto;
  NcDataStruct *extra_data;
};

/**
 * NcDataGaussianType:
 * @NC_DATA_GAUSSIAN_SIGMA: FIXME
 * @NC_DATA_GAUSSIAN_COV: FIXME
 * @NC_DATA_GAUSSIAN_X_SIGMA: FIXME
 * @NC_DATA_GAUSSIAN_X_COV: FIXME
 * @NC_DATA_GAUSSIAN_X_WMEAN: FIXME
 *
 * FIXME
 */
typedef enum _NcDataGaussianType
{
  NC_DATA_GAUSSIAN_SIGMA,
  NC_DATA_GAUSSIAN_COV,
  NC_DATA_GAUSSIAN_X_SIGMA,
  NC_DATA_GAUSSIAN_X_COV,
  NC_DATA_GAUSSIAN_X_WMEAN,
} NcDataGaussianType;


NcData *nc_data_gaussian_new (NcDataGaussianType gauss_type);

void nc_data_gaussian_set_func_array (NcData *data, NcmMSetFunc **func, guint n);
void nc_data_gaussian_set_func (NcData *data, NcmMSetFunc *func);

void nc_data_gaussian_init_from_matrix (NcData *data, NcmMatrix *cm, NcDataStruct *extra_data);
#ifdef NUMCOSMO_HAVE_SQLITE3
void nc_data_gaussian_init_from_query (NcData *data, sqlite3 *db, gchar *query, NcDataStruct *extra_data);
#endif

G_END_DECLS

#endif /* NC_DATA_GAUSSIAN_H */

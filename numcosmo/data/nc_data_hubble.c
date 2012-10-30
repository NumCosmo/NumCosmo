/***************************************************************************
 *            nc_data_hubble.c
 *
 *  Thu Apr 22 14:34:54 2010
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
 * SECTION:nc_data_hubble
 * @title: Hubble Function Data
 * @short_description: Object representing Hubble Function data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_hubble.h"
#include "data/data_gaussian.h"
#include "math/ncm_mset_func.h"
#include "math/ncm_cfg.h"
#include "nc_hicosmo.h"

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

static gchar *_nc_data_hubble_function_query[] =
{
  "L. Verde sample", "SELECT z,p,s FROM kinematics WHERE param='H' ORDER BY z",
  "Cabre sample", "SELECT z,p,s FROM kinematics WHERE param='H_CABRE' ORDER BY z"
};

/**
 * nc_data_hubble:
 * @H_id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_hubble (NcDataHubbleId H_id)
{
  NcData *H_data = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_SIGMA);
  NcmMSetFunc *func = nc_hicosmo_func1_new (&nc_hicosmo_H);
  g_assert (H_id < NC_DATA_HUBBLE_NSAMPLES);

  nc_data_gaussian_set_func (H_data, func);
  ncm_mset_func_free (func);

#ifdef NUMCOSMO_HAVE_SQLITE3
  sqlite3 *db = ncm_cfg_get_default_sqlite3 ();
  H_data->name = _nc_data_hubble_function_query[H_id * 2];
  nc_data_gaussian_init_from_query (H_data, db, _nc_data_hubble_function_query[H_id * 2 + 1], NULL);
#else
  g_error (PACKAGE_NAME" compiled without support for sqlite3, SN Ia data not avaliable.");
#endif

  return H_data;
}

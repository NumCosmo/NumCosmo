/***************************************************************************
 *            nc_data_distance_mu.c
 *
 *  Thu Apr 22 10:37:22 2010
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
 * SECTION:nc_data_distance_mu
 * @title: Distance Modulus Data
 * @short_description: Data samples of distance modulus
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_distance_mu.h"
#include "data/data_gaussian.h"
#include "math/ncm_cfg.h"

#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

static gchar *_nc_data_snia_query[] =
{
  "Gold sample 157", "SELECT z,mu-0.32 AS muc,s FROM supernovae WHERE quality='Gold-2004' ORDER BY z",
  "Gold sample 182 - removed low redshift", "SELECT z,mu-0.32 AS muc,s FROM supernovae WHERE quality='Gold' AND z >= 0.0233 ORDER BY z",
  "Gold sample 182", "SELECT z,mu-0.32 AS muc,s FROM supernovae WHERE quality='Gold' ORDER BY z",
  "ESSENCE sample",  "SELECT z,mu,s FROM supernovae WHERE quality='ESSENCE' ORDER BY z",
  "Legacy sample",   "SELECT z,mu+19.308,s FROM supernovae WHERE quality='LEGACY' ORDER BY z",
  "Union sample",    "SELECT z,mu,s FROM supernovae WHERE quality='UNION' ORDER BY z",
  "CfA3 sample",     "SELECT z,mu,s FROM supernovae WHERE quality='CfA3' ORDER BY z",
  "Union2 sample",   "SELECT z,mu,s FROM supernovae WHERE quality='Union2' ORDER BY z",
  "Union2.1 sample", "SELECT z,mu,s FROM supernovae WHERE quality='Union2.1' ORDER BY z",
};

/**
 * nc_data_distance_modulus_snia:
 * @dist: a #NcDistance
 * @snia_id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_distance_mu_snia (NcDistance *dist, NcDataDistanceMuSNIaId snia_id)
{
  NcData *snia = nc_data_gaussian_new (NC_DATA_GAUSSIAN_X_WMEAN);
  NcmMSetFunc *dist_mod = nc_distance_func1_new (dist, &nc_distance_modulo);
  g_assert (snia_id < NC_DATA_DISTANCE_MU_SNIA_NSAMPLES);

  nc_data_gaussian_set_func (snia, dist_mod);
  ncm_mset_func_free (dist_mod);

#ifdef NUMCOSMO_HAVE_SQLITE3
  sqlite3 *db = ncm_cfg_get_default_sqlite3 ();
  snia->name = _nc_data_snia_query[snia_id * 2];
  nc_data_gaussian_init_from_query (snia, db, _nc_data_snia_query[snia_id * 2 + 1], NULL);
#else
  g_error (PACKAGE_NAME" compiled without support for sqlite3, SN Ia data not avaliable.");
#endif

  return snia;
}

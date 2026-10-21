/***************************************************************************
 *            nc_data_hubble_cov.h
 *
 *  Wed Oct 21 13:26:04 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_hubble_cov.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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

#ifndef _NC_DATA_HUBBLE_COV_H_
#define _NC_DATA_HUBBLE_COV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_data_gauss_cov.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_HUBBLE_COV (nc_data_hubble_cov_get_type ())

G_DECLARE_FINAL_TYPE (NcDataHubbleCov, nc_data_hubble_cov, NC, DATA_HUBBLE_COV, NcmDataGaussCov)

/**
 * NcDataHubbleCovId:
 * @NC_DATA_HUBBLE_COV_LOUBSER2025: [Loubser (2025)](https://arxiv.org/abs/2511.02730), DESI DR1 cosmic chronometers
 * @NC_DATA_HUBBLE_COV_MORESCO2020: D4000n cosmic chronometers of Moresco et al. (2012, 2015, 2016) with the covariance of [Moresco et al. (2020)](https://arxiv.org/abs/2003.07362)
 *
 * Identifiers for Hubble parameter datasets with a full covariance.
 */
typedef enum _NcDataHubbleCovId /*< prefix=NC_DATA_HUBBLE_COV >*/
{
  NC_DATA_HUBBLE_COV_LOUBSER2025 = 0,
  NC_DATA_HUBBLE_COV_MORESCO2020,
  /* < private > */
  NC_DATA_HUBBLE_COV_NSAMPLES, /*< skip >*/
} NcDataHubbleCovId;

NcDataHubbleCov *nc_data_hubble_cov_new_empty (void);
NcDataHubbleCov *nc_data_hubble_cov_new_from_file (const gchar *filename);
NcDataHubbleCov *nc_data_hubble_cov_new_from_id (NcDataHubbleCovId id);

G_END_DECLS

#endif /* _NC_DATA_HUBBLE_COV_H_ */


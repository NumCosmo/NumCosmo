/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_gauss.h
 *
 *  Thu Aug 1 20:02:19 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_gauss.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_H_
#define _NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/galaxy/nc_galaxy_sd_obs_redshift.h>
#include <numcosmo/galaxy/nc_galaxy_sd_true_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_OBS_REDSHIFT_GAUSS (nc_galaxy_sd_obs_redshift_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDObsRedshiftGauss, nc_galaxy_sd_obs_redshift_gauss, NC, GALAXY_SD_OBS_REDSHIFT_GAUSS, NcGalaxySDObsRedshift)

/**
 * NcGalaxySDObsRedshiftGaussParams:
 * @NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SIGMA: Standard deviation of the gaussian distribution
 *
 * Photometric redshift observation with gaussian errors model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_PARAMS >*/
{
  NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SIGMA = 0,
  /* < private > */
  NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SPARAM_LEN, /*< skip >*/
} NcGalaxySDObsRedshiftGaussParams;

#define NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_DEFAULT_SIGMA  (0.05)

#define NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDObsRedshiftGauss *nc_galaxy_sd_obs_redshift_gauss_new (NcGalaxySDTrueRedshift *sdz);
NcGalaxySDObsRedshiftGauss *nc_galaxy_sd_obs_redshift_gauss_ref (NcGalaxySDObsRedshiftGauss *gsdorgauss);

void nc_galaxy_sd_obs_redshift_gauss_free (NcGalaxySDObsRedshiftGauss *gsdorgauss);
void nc_galaxy_sd_obs_redshift_gauss_clear (NcGalaxySDObsRedshiftGauss **gsdorgauss);

G_END_DECLS

#endif /* _NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_H_ */

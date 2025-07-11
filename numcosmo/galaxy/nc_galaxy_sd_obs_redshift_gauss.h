/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_gauss.h
 *
 *  Thu Aug 1 20:02:19 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_gauss.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

NcGalaxySDObsRedshiftGauss *nc_galaxy_sd_obs_redshift_gauss_new (NcGalaxySDTrueRedshift * sdz, const gdouble zp_min, const gdouble zp_max);
NcGalaxySDObsRedshiftGauss *nc_galaxy_sd_obs_redshift_gauss_ref (NcGalaxySDObsRedshiftGauss *gsdorgauss);

void nc_galaxy_sd_obs_redshift_gauss_free (NcGalaxySDObsRedshiftGauss *gsdorgauss);
void nc_galaxy_sd_obs_redshift_gauss_clear (NcGalaxySDObsRedshiftGauss **gsdorgauss);

void nc_galaxy_sd_obs_redshift_gauss_set_lim (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gdouble zp_min, const gdouble zp_max);
void nc_galaxy_sd_obs_redshift_gauss_get_lim (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble *zp_min, gdouble *zp_max);

void nc_galaxy_sd_obs_redshift_gauss_set_use_true_z (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gboolean use_true_z);
gboolean nc_galaxy_sd_obs_redshift_gauss_get_use_true_z (NcGalaxySDObsRedshiftGauss *gsdorgauss);

void nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, const gdouble sigma0, NcmRNG *rng);
gboolean nc_galaxy_sd_obs_redshift_gauss_gen1 (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, const gdouble sigma0, NcmRNG *rng);

void nc_galaxy_sd_obs_redshift_gauss_data_set (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, const gdouble zp, const gdouble sigma0, const gdouble sigma_z);
void nc_galaxy_sd_obs_redshift_gauss_data_get (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, gdouble *zp, gdouble *sigma0, gdouble *sigma_z);

#define NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP "zp"
#define NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA "sigma_z"
#define NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0 "sigma_0"


G_END_DECLS

#endif /* _NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_H_ */


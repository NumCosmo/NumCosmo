/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_pz.h
 *
 *  Mon Nov 25 20:35:34 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_pz.h
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

#ifndef _NC_GALAXY_SD_OBS_REDSHIFT_PZ_H_
#define _NC_GALAXY_SD_OBS_REDSHIFT_PZ_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/galaxy/nc_galaxy_sd_obs_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_OBS_REDSHIFT_PZ (nc_galaxy_sd_obs_redshift_pz_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDObsRedshiftPz, nc_galaxy_sd_obs_redshift_pz, NC, GALAXY_SD_OBS_REDSHIFT_PZ, NcGalaxySDObsRedshift)

NcGalaxySDObsRedshiftPz *nc_galaxy_sd_obs_redshift_pz_new ();
NcGalaxySDObsRedshiftPz *nc_galaxy_sd_obs_redshift_pz_ref (NcGalaxySDObsRedshiftPz *gsdorpz);

void nc_galaxy_sd_obs_redshift_pz_free (NcGalaxySDObsRedshiftPz *gsdorpz);
void nc_galaxy_sd_obs_redshift_pz_clear (NcGalaxySDObsRedshiftPz **gsdorpz);

void nc_galaxy_sd_obs_redshift_pz_data_set (NcGalaxySDObsRedshiftPz *gsdorpz, NcGalaxySDObsRedshiftData *data, NcmSpline *spline);
void nc_galaxy_sd_obs_redshift_pz_data_get (NcGalaxySDObsRedshiftPz *gsdorpz, NcGalaxySDObsRedshiftData *data, NcmSpline **spline);

G_END_DECLS

#endif /* _NC_GALAXY_SD_OBS_REDSHIFT_PZ_H_ */


/***************************************************************************
 *            nc_galaxy_sd_position_srd_y1.h
 *
 *  Tue June 22 15:02:20 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_srd_y1.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_SD_POSITION_SRD_Y1_H_
#define _NC_GALAXY_SD_POSITION_SRD_Y1_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION_SRD_Y1 (nc_galaxy_sd_position_srd_y1_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDPositionSRDY1, nc_galaxy_sd_position_srd_y1, NC, GALAXY_SD_POSITION_SRD_Y1, NcGalaxySDPosition)

/**
 * NcGalaxySDPositionSRDY1SParams:
 * @NC_GALAXY_SD_POSITION_SRD_Y1_ALPHA: redshift exponential slope
 * @NC_GALAXY_SD_POSITION_SRD_Y1_BETA: redshift power law slope
 * @NC_GALAXY_SD_POSITION_SRD_Y1_Z0: Pivot redshift
 *
 * LSST SRD galaxy redshift distribution model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_POSITION_SPARAMS >*/
{
  NC_GALAXY_SD_POSITION_SRD_Y1_ALPHA = 0,
  NC_GALAXY_SD_POSITION_SRD_Y1_BETA,
  NC_GALAXY_SD_POSITION_SRD_Y1_Z0,
  /* < private > */
  NC_GALAXY_SD_POSITION_SRD_Y1_SPARAM_LEN, /*< skip >*/
} NcGalaxySDPositionSParams;

#define NC_GALAXY_SD_POSITION_SRD_Y1_DEFAULT_ALPHA  (0.78)
#define NC_GALAXY_SD_POSITION_SRD_Y1_DEFAULT_BETA   (2.00)
#define NC_GALAXY_SD_POSITION_SRD_Y1_DEFAULT_Z0     (0.13)

#define NC_GALAXY_SD_POSITION_SRD_Y1_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDPositionSRDY1 *nc_galaxy_sd_position_srd_y1_new ();
NcGalaxySDPositionSRDY1 *nc_galaxy_sd_position_srd_y1_ref (NcGalaxySDPositionSRDY1 *gsdpsrdy1);

void nc_galaxy_sd_position_srd_y1_free (NcGalaxySDPositionSRDY1 *gsdpsrdy1);
void nc_galaxy_sd_position_srd_y1_clear (NcGalaxySDPositionSRDY1 **gsdpsrdy1);

void nc_galaxy_sd_position_srd_y1_set_z_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1, NcmVector *lim);
NcmVector *nc_galaxy_sd_position_srd_y1_peek_z_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1);
void nc_galaxy_sd_position_srd_y1_set_r_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1, NcmVector *lim);
NcmVector *nc_galaxy_sd_position_srd_y1_peek_r_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1);

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_SRD_Y1_H_ */


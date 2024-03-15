/***************************************************************************
 *            nc_galaxy_sd_position_lsst_srd.h
 *
 *  Tue June 22 15:02:20 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_lsst_srd.h
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

#ifndef _NC_GALAXY_SD_POSITION_LSST_SRD_H_
#define _NC_GALAXY_SD_POSITION_LSST_SRD_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION_LSST_SRD (nc_galaxy_sd_position_lsst_srd_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDPositionLSSTSRD, nc_galaxy_sd_position_lsst_srd, NC, GALAXY_SD_POSITION_LSST_SRD, NcGalaxySDPosition)

/**
 * NcGalaxySDPositionLSSTSRDSParams:
 * @NC_GALAXY_SD_POSITION_LSST_SRD_ALPHA: redshift exponential slope
 * @NC_GALAXY_SD_POSITION_LSST_SRD_BETA: redshift power law slope
 * @NC_GALAXY_SD_POSITION_LSST_SRD_Z0: Pivot redshift
 *
 * LSST SRD galaxy redshift distribution model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_POSITION_SPARAMS >*/
{
  NC_GALAXY_SD_POSITION_LSST_SRD_ALPHA = 0,
  NC_GALAXY_SD_POSITION_LSST_SRD_BETA,
  NC_GALAXY_SD_POSITION_LSST_SRD_Z0,
  /* < private > */
  NC_GALAXY_SD_POSITION_LSST_SRD_SPARAM_LEN, /*< skip >*/
} NcGalaxySDPositionSParams;

#define NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_ALPHA  (0.78)
#define NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_BETA   (2.00)
#define NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_Z0     (0.13)

#define NC_GALAXY_SD_POSITION_LSST_SRD_Y10_ALPHA  (0.68)
#define NC_GALAXY_SD_POSITION_LSST_SRD_Y10_BETA   (2.00)
#define NC_GALAXY_SD_POSITION_LSST_SRD_Y10_Z0     (0.11)

#define NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDPositionLSSTSRD *nc_galaxy_sd_position_lsst_srd_new (const gdouble z_min, const gdouble z_max, const gdouble r_min, const gdouble r_max);
NcGalaxySDPositionLSSTSRD *nc_galaxy_sd_position_lsst_srd_new_y10 (const gdouble z_min, const gdouble z_max, const gdouble r_min, const gdouble r_max);
NcGalaxySDPositionLSSTSRD *nc_galaxy_sd_position_lsst_srd_ref (NcGalaxySDPositionLSSTSRD *gsdplsst);

void nc_galaxy_sd_position_lsst_srd_free (NcGalaxySDPositionLSSTSRD *gsdplsst);
void nc_galaxy_sd_position_lsst_srd_clear (NcGalaxySDPositionLSSTSRD **gsdplsst);

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_LSST_SRD_H_ */


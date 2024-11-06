/***************************************************************************
 *            nc_galaxy_sd_true_redshift_lsst_srd.h
 *
 *  Wed Jul 31 21:35:59 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_true_redshift_lssttroscopic.h
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

#ifndef _NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_H
#define _NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/galaxy/nc_galaxy_sd_true_redshift.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (nc_galaxy_sd_true_redshift_lsst_srd_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDTrueRedshiftLSSTSRD, nc_galaxy_sd_true_redshift_lsst_srd, NC, GALAXY_SD_TRUE_REDSHIFT_LSST_SRD, NcGalaxySDTrueRedshift)

/**
 * NcGalaxySDTrueRedshiftLSSTSRDSParams:
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_ALPHA: Alpha parameter
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_BETA: Beta parameter
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Z0: Z0 parameter
 *
 * LSST SRD galaxy redshift distribution model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_PARAMS >*/
{
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_ALPHA = 0,
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_BETA,
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Z0,
  /* < private > */
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_SPARAM_LEN, /*< skip >*/
} NcGalaxySDTrueRedshiftLSSTSRDSParams;

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_ALPHA  (0.78)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_BETA   (2.00)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z0     (0.13)

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_ALPHA  (0.68)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_BETA   (2.00)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_Z0     (0.11)

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_PARAMS_ABSTOL (0.0)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z_LOW         (0.0)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z_HIGH        (20.0)

NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new_y10 (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_ref (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst);

void nc_galaxy_sd_true_redshift_lsst_srd_free (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst);
void nc_galaxy_sd_true_redshift_lsst_srd_clear (NcGalaxySDTrueRedshiftLSSTSRD **gsdtrlsst);

G_END_DECLS

#endif /* _NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_H */


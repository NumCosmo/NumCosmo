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
 * nc_galaxy_sd_true_redshift_lsst_srd.h
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
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_true_redshift.h>
#include <numcosmo/ncm/core/ncm_rng.h>

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

/**
 * NcGalaxySDTrueRedshiftLSSTSRDType:
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE: Year 1 source parametrization
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS: Year 1 lens parametrization
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE: Year 10 source parametrization
 * @NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS: Year 10 lens parametrization
 *
 * LSST SRD galaxy redshift distribution types.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_TYPE >*/
{
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE = 0,
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS,
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE,
  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS,
} NcGalaxySDTrueRedshiftLSSTSRDType;

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE_ALPHA  (0.78)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE_BETA   (2.00)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE_Z0     (0.13)

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS_ALPHA  (0.94)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS_BETA   (2.00)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS_Z0     (0.26)

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE_ALPHA  (0.68)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE_BETA   (2.00)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE_Z0     (0.11)

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS_ALPHA  (0.90)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS_BETA   (2.00)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS_Z0     (0.28)

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_ALPHA  NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE_ALPHA
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_BETA   NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE_BETA
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z0     NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE_Z0

#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_PARAMS_ABSTOL (0.0)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z_LOW         (0.0)
#define NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z_HIGH        (20.0)

NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new_y1_source (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new_y1_lens (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new_y10_source (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new_y10_lens (void);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_new_from_type (NcGalaxySDTrueRedshiftLSSTSRDType type);
NcGalaxySDTrueRedshiftLSSTSRD *nc_galaxy_sd_true_redshift_lsst_srd_ref (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst);

void nc_galaxy_sd_true_redshift_lsst_srd_free (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst);
void nc_galaxy_sd_true_redshift_lsst_srd_clear (NcGalaxySDTrueRedshiftLSSTSRD **gsdtrlsst);

G_END_DECLS

#endif /* _NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_H */


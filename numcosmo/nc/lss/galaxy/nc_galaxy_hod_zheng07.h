/***************************************************************************
 *            nc_galaxy_hod_zheng07.h
 *
 *  Sun Jun 14 12:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_hod_zheng07.h
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

#ifndef _NC_GALAXY_HOD_ZHENG07_H
#define _NC_GALAXY_HOD_ZHENG07_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_hod.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_HOD_ZHENG07 (nc_galaxy_hod_zheng07_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyHODZheng07, nc_galaxy_hod_zheng07, NC, GALAXY_HOD_ZHENG07, NcGalaxyHOD)

/**
 * NcGalaxyHODZheng07SParams:
 * @NC_GALAXY_HOD_ZHENG07_LOG_MMIN: $\log_{10}$ of the central cutoff mass.
 * @NC_GALAXY_HOD_ZHENG07_SIGMA_LOG_M: width of the central transition.
 * @NC_GALAXY_HOD_ZHENG07_LOG_M0: $\log_{10}$ of the satellite cutoff mass.
 * @NC_GALAXY_HOD_ZHENG07_LOG_M1: $\log_{10}$ of the satellite normalization mass.
 * @NC_GALAXY_HOD_ZHENG07_ALPHA: satellite power-law slope.
 *
 * Zheng et al. (2007) HOD parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_HOD_ZHENG07_PARAMS >*/
{
  NC_GALAXY_HOD_ZHENG07_LOG_MMIN = 0,
  NC_GALAXY_HOD_ZHENG07_SIGMA_LOG_M,
  NC_GALAXY_HOD_ZHENG07_LOG_M0,
  NC_GALAXY_HOD_ZHENG07_LOG_M1,
  NC_GALAXY_HOD_ZHENG07_ALPHA,
  /* < private > */
  NC_GALAXY_HOD_ZHENG07_SPARAM_LEN, /*< skip >*/
} NcGalaxyHODZheng07SParams;

#define NC_GALAXY_HOD_ZHENG07_DEFAULT_LOG_MMIN (12.72)
#define NC_GALAXY_HOD_ZHENG07_DEFAULT_SIGMA_LOG_M (0.26)
#define NC_GALAXY_HOD_ZHENG07_DEFAULT_LOG_M0 (12.7)
#define NC_GALAXY_HOD_ZHENG07_DEFAULT_LOG_M1 (13.93)
#define NC_GALAXY_HOD_ZHENG07_DEFAULT_ALPHA (1.15)

#define NC_GALAXY_HOD_ZHENG07_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxyHODZheng07 *nc_galaxy_hod_zheng07_new (void);
NcGalaxyHODZheng07 *nc_galaxy_hod_zheng07_ref (NcGalaxyHODZheng07 *zheng07);

void nc_galaxy_hod_zheng07_free (NcGalaxyHODZheng07 *zheng07);
void nc_galaxy_hod_zheng07_clear (NcGalaxyHODZheng07 **zheng07);

G_END_DECLS

#endif /* _NC_GALAXY_HOD_ZHENG07_H */


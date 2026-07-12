/***************************************************************************
 *            nc_galaxy_redshift_obs_gauss.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs_gauss.h
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NC_GALAXY_REDSHIFT_OBS_GAUSS_H_
#define _NC_GALAXY_REDSHIFT_OBS_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_obs.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_OBS_GAUSS (nc_galaxy_redshift_obs_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyRedshiftObsGauss, nc_galaxy_redshift_obs_gauss, NC, GALAXY_REDSHIFT_OBS_GAUSS, NcGalaxyRedshiftObs)

NcGalaxyRedshiftObsGauss *nc_galaxy_redshift_obs_gauss_new (void);
NcGalaxyRedshiftObsGauss *nc_galaxy_redshift_obs_gauss_ref (NcGalaxyRedshiftObsGauss *gsdreg);

void nc_galaxy_redshift_obs_gauss_free (NcGalaxyRedshiftObsGauss *gsdreg);
void nc_galaxy_redshift_obs_gauss_clear (NcGalaxyRedshiftObsGauss **gsdreg);

void nc_galaxy_redshift_obs_gauss_data_set (NcGalaxyRedshiftObsGauss *gsdreg, NcGalaxyRedshiftObsData *data, const gdouble zp, const gdouble sigma0);
void nc_galaxy_redshift_obs_gauss_data_get (NcGalaxyRedshiftObsGauss *gsdreg, NcGalaxyRedshiftObsData *data, gdouble *zp, gdouble *sigma0);

#define NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_ZP "zp"
#define NC_GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0 "sigma0"

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_OBS_GAUSS_H_ */

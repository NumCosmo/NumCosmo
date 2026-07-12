/***************************************************************************
 *            nc_galaxy_redshift_obs_sel_gauss.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs_sel_gauss.h
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

#ifndef _NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_H_
#define _NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_obs_sel.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_OBS_SEL_GAUSS (nc_galaxy_redshift_obs_sel_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyRedshiftObsSelGauss, nc_galaxy_redshift_obs_sel_gauss, NC, GALAXY_REDSHIFT_OBS_SEL_GAUSS, NcGalaxyRedshiftObsSel)

/**
 * NcGalaxyRedshiftObsSelGaussSParams:
 * @NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SIGMA0: the population photo-z scatter, sigma_z = sigma0 (1 + z)
 *
 * Gaussian population photo-z observable model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SPARAMS >*/
{
  NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SIGMA0 = 0,
  /* < private > */
  NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_SPARAM_LEN, /*< skip >*/
} NcGalaxyRedshiftObsSelGaussSParams;

#define NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_DEFAULT_SIGMA0 (0.05)
#define NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxyRedshiftObsSelGauss *nc_galaxy_redshift_obs_sel_gauss_new (void);
NcGalaxyRedshiftObsSelGauss *nc_galaxy_redshift_obs_sel_gauss_ref (NcGalaxyRedshiftObsSelGauss *gsdropg);

void nc_galaxy_redshift_obs_sel_gauss_free (NcGalaxyRedshiftObsSelGauss *gsdropg);
void nc_galaxy_redshift_obs_sel_gauss_clear (NcGalaxyRedshiftObsSelGauss **gsdropg);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_OBS_SEL_GAUSS_H_ */

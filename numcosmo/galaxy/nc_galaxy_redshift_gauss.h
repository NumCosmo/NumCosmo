/***************************************************************************
 *            nc_galaxy_redshift_gauss.h
 *
 *  Thu Aug 1 20:02:19 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_gauss.h
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

#ifndef _NC_GALAXY_REDSHIFT_GAUSS_H_
#define _NC_GALAXY_REDSHIFT_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/galaxy/nc_galaxy_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_GAUSS (nc_galaxy_redshift_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyRedshiftGauss, nc_galaxy_redshift_gauss, NC, GALAXY_REDSHIFT_GAUSS, NcGalaxyRedshift)

NcGalaxyRedshiftGauss *nc_galaxy_redshift_gauss_new ();
NcGalaxyRedshiftGauss *nc_galaxy_redshift_gauss_ref (NcGalaxyRedshiftGauss *gzgauss);

void nc_galaxy_redshift_gauss_free (NcGalaxyRedshiftGauss *gzgauss);
void nc_galaxy_redshift_gauss_clear (NcGalaxyRedshiftGauss **gzgauss);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_GAUSS_H_ */


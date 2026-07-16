/***************************************************************************
 *            nc_galaxy_redshift_factor_spline.h
 *
 *  Sun Jul 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor_spline.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_H_
#define _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_factor.h>
#include <numcosmo/ncm/spline/ncm_spline.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_FACTOR_SPLINE (nc_galaxy_redshift_factor_spline_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyRedshiftFactorSpline, nc_galaxy_redshift_factor_spline, NC, GALAXY_REDSHIFT_FACTOR_SPLINE, NcGalaxyRedshiftFactor)

NcGalaxyRedshiftFactorSpline *nc_galaxy_redshift_factor_spline_new (void);
NcGalaxyRedshiftFactorSpline *nc_galaxy_redshift_factor_spline_ref (NcGalaxyRedshiftFactorSpline *gsdrs);

void nc_galaxy_redshift_factor_spline_free (NcGalaxyRedshiftFactorSpline *gsdrs);
void nc_galaxy_redshift_factor_spline_clear (NcGalaxyRedshiftFactorSpline **gsdrs);

void nc_galaxy_redshift_factor_spline_data_set (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data, NcmSpline *pz);
NcmSpline *nc_galaxy_redshift_factor_spline_data_peek (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_FACTOR_SPLINE_H_ */


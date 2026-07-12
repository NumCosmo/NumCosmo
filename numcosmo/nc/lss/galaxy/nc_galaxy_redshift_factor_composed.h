/***************************************************************************
 *            nc_galaxy_redshift_factor_composed.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor_composed.h
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

#ifndef _NC_GALAXY_REDSHIFT_FACTOR_COMPOSED_H_
#define _NC_GALAXY_REDSHIFT_FACTOR_COMPOSED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_factor.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_pop.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_obs.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_FACTOR_COMPOSED (nc_galaxy_redshift_factor_composed_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyRedshiftFactorComposed, nc_galaxy_redshift_factor_composed, NC, GALAXY_REDSHIFT_FACTOR_COMPOSED, NcGalaxyRedshiftFactor)

NcGalaxyRedshiftFactorComposed *nc_galaxy_redshift_factor_composed_new (const gdouble zp_min, const gdouble zp_max);
NcGalaxyRedshiftFactorComposed *nc_galaxy_redshift_factor_composed_ref (NcGalaxyRedshiftFactorComposed *gsdrc);

void nc_galaxy_redshift_factor_composed_free (NcGalaxyRedshiftFactorComposed *gsdrc);
void nc_galaxy_redshift_factor_composed_clear (NcGalaxyRedshiftFactorComposed **gsdrc);

void nc_galaxy_redshift_factor_composed_set_zp_lim (NcGalaxyRedshiftFactorComposed *gsdrc, const gdouble zp_min, const gdouble zp_max);
void nc_galaxy_redshift_factor_composed_get_zp_lim (NcGalaxyRedshiftFactorComposed *gsdrc, gdouble *zp_min, gdouble *zp_max);

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_FACTOR_COMPOSED_H_ */


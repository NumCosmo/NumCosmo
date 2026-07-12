/***************************************************************************
 *            nc_galaxy_position_factor_flat.h
 *
 *  Wed Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_position_factor_flat.h
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
#ifndef _NC_GALAXY_POSITION_FACTOR_FLAT_H_
#define _NC_GALAXY_POSITION_FACTOR_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_position_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_POSITION_FACTOR_FLAT (nc_galaxy_position_factor_flat_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyPositionFactorFlat, nc_galaxy_position_factor_flat, NC, GALAXY_POSITION_FACTOR_FLAT, NcGalaxyPositionFactor)

NcGalaxyPositionFactorFlat *nc_galaxy_position_factor_flat_new (const gdouble ra_min, const gdouble ra_max, const gdouble dec_min, const gdouble dec_max);
NcGalaxyPositionFactorFlat *nc_galaxy_position_factor_flat_ref (NcGalaxyPositionFactorFlat *gspfflat);

void nc_galaxy_position_factor_flat_free (NcGalaxyPositionFactorFlat *gspfflat);
void nc_galaxy_position_factor_flat_clear (NcGalaxyPositionFactorFlat **gspfflat);

void nc_galaxy_position_factor_flat_set_ra_lim (NcGalaxyPositionFactorFlat *gspfflat, const gdouble ra_min, const gdouble ra_max);
void nc_galaxy_position_factor_flat_get_ra_lim (NcGalaxyPositionFactorFlat *gspfflat, gdouble *ra_min, gdouble *ra_max);
void nc_galaxy_position_factor_flat_set_dec_lim (NcGalaxyPositionFactorFlat *gspfflat, const gdouble dec_min, const gdouble dec_max);
void nc_galaxy_position_factor_flat_get_dec_lim (NcGalaxyPositionFactorFlat *gspfflat, gdouble *dec_min, gdouble *dec_max);

G_END_DECLS

#endif /* _NC_GALAXY_POSITION_FACTOR_FLAT_H_ */

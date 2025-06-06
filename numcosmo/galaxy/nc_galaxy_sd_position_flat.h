/***************************************************************************
 *            nc_galaxy_sd_position_flat.h
 *
 *  Wed March 1 12:53:13 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_flat.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_SD_POSITION_FLAT_H_
#define _NC_GALAXY_SD_POSITION_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_POSITION_FLAT (nc_galaxy_sd_position_flat_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDPositionFlat, nc_galaxy_sd_position_flat, NC, GALAXY_SD_POSITION_FLAT, NcGalaxySDPosition)

NcGalaxySDPositionFlat *nc_galaxy_sd_position_flat_new (const gdouble ra_min, const gdouble ra_max, const gdouble dec_min, const gdouble dec_max);
NcGalaxySDPositionFlat *nc_galaxy_sd_position_flat_ref (NcGalaxySDPositionFlat *gsdpflat);

void nc_galaxy_sd_position_flat_free (NcGalaxySDPositionFlat *gsdpflat);
void nc_galaxy_sd_position_flat_clear (NcGalaxySDPositionFlat **gsdpflat);

void nc_galaxy_sd_position_flat_gen (NcGalaxySDPositionFlat *gsdpflat, NcmMSet *mset, NcGalaxySDPositionData *data, NcmRNG *rng);

G_END_DECLS

#endif /* _NC_GALAXY_SD_POSITION_FLAT_H_ */


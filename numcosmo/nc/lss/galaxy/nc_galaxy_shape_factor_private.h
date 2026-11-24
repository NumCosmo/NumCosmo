/***************************************************************************
 *            nc_galaxy_shape_factor_private.h
 *
 *  Tue Nov 24 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_private.h
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

/*
 * The critical densities at the redshift nodes of a galaxy as an immutable block, so that
 * the copies of a dataset made for a sampler share it. Not public API.
 */
#ifndef _NC_GALAXY_SHAPE_FACTOR_PRIVATE_H_
#define _NC_GALAXY_SHAPE_FACTOR_PRIVATE_H_

#ifndef NUMCOSMO_GIR_SCAN

#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

GBytes *_nc_galaxy_shape_factor_data_peek_crit (NcGalaxyShapeFactorData *data);

void _nc_galaxy_shape_factor_data_share_crit (NcGalaxyShapeFactorData *data, GBytes *crit, const guint n_nodes);

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _NC_GALAXY_SHAPE_FACTOR_PRIVATE_H_ */


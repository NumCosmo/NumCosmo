/***************************************************************************
 *            nc_galaxy_shape_factor_tilted_series.h
 *
 *  Thu Sep 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_tilted_series.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_H_
#define _NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_TILTED_SERIES (nc_galaxy_shape_factor_tilted_series_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorTiltedSeries, nc_galaxy_shape_factor_tilted_series, NC, GALAXY_SHAPE_FACTOR_TILTED_SERIES, NcGalaxyShapeFactor)

NcGalaxyShapeFactorTiltedSeries *nc_galaxy_shape_factor_tilted_series_new (NcGalaxyWLObsEllipConv ellip_conv, guint trunc_order);
NcGalaxyShapeFactorTiltedSeries *nc_galaxy_shape_factor_tilted_series_ref (NcGalaxyShapeFactorTiltedSeries *gsfts);

void nc_galaxy_shape_factor_tilted_series_free (NcGalaxyShapeFactorTiltedSeries *gsfts);
void nc_galaxy_shape_factor_tilted_series_clear (NcGalaxyShapeFactorTiltedSeries **gsfts);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_TILTED_SERIES_H_ */

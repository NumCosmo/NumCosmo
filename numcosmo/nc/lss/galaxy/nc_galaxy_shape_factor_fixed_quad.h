/***************************************************************************
 *            nc_galaxy_shape_factor_fixed_quad.h
 *
 *  Thu Jul 9 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_fixed_quad.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_H_
#define _NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_FIXED_QUAD (nc_galaxy_shape_factor_fixed_quad_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorFixedQuad, nc_galaxy_shape_factor_fixed_quad, NC, GALAXY_SHAPE_FACTOR_FIXED_QUAD, NcGalaxyShapeFactor)

NcGalaxyShapeFactorFixedQuad *nc_galaxy_shape_factor_fixed_quad_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorFixedQuad *nc_galaxy_shape_factor_fixed_quad_ref (NcGalaxyShapeFactorFixedQuad *gsffq);

void nc_galaxy_shape_factor_fixed_quad_free (NcGalaxyShapeFactorFixedQuad *gsffq);
void nc_galaxy_shape_factor_fixed_quad_clear (NcGalaxyShapeFactorFixedQuad **gsffq);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_FIXED_QUAD_H_ */


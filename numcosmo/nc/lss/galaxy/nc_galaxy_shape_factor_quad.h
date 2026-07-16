/***************************************************************************
 *            nc_galaxy_shape_factor_quad.h
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_quad.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_QUAD_H_
#define _NC_GALAXY_SHAPE_FACTOR_QUAD_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_QUAD (nc_galaxy_shape_factor_quad_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorQuad, nc_galaxy_shape_factor_quad, NC, GALAXY_SHAPE_FACTOR_QUAD, NcGalaxyShapeFactor)

NcGalaxyShapeFactorQuad *nc_galaxy_shape_factor_quad_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorQuad *nc_galaxy_shape_factor_quad_ref (NcGalaxyShapeFactorQuad *gsfq);

void nc_galaxy_shape_factor_quad_free (NcGalaxyShapeFactorQuad *gsfq);
void nc_galaxy_shape_factor_quad_clear (NcGalaxyShapeFactorQuad **gsfq);

void nc_galaxy_shape_factor_quad_set_bound (NcGalaxyShapeFactorQuad *gsfq, const gdouble bound);
gdouble nc_galaxy_shape_factor_quad_get_bound (NcGalaxyShapeFactorQuad *gsfq);

void nc_galaxy_shape_factor_quad_set_reltol (NcGalaxyShapeFactorQuad *gsfq, const gdouble reltol);
gdouble nc_galaxy_shape_factor_quad_get_reltol (NcGalaxyShapeFactorQuad *gsfq);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_QUAD_H_ */


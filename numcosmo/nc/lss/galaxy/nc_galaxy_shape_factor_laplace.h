/***************************************************************************
 *            nc_galaxy_shape_factor_laplace.h
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_laplace.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_LAPLACE_H_
#define _NC_GALAXY_SHAPE_FACTOR_LAPLACE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_LAPLACE (nc_galaxy_shape_factor_laplace_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorLaplace, nc_galaxy_shape_factor_laplace, NC, GALAXY_SHAPE_FACTOR_LAPLACE, NcGalaxyShapeFactor)

NcGalaxyShapeFactorLaplace *nc_galaxy_shape_factor_laplace_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorLaplace *nc_galaxy_shape_factor_laplace_ref (NcGalaxyShapeFactorLaplace *gsfl);

void nc_galaxy_shape_factor_laplace_free (NcGalaxyShapeFactorLaplace *gsfl);
void nc_galaxy_shape_factor_laplace_clear (NcGalaxyShapeFactorLaplace **gsfl);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_LAPLACE_H_ */


/***************************************************************************
 *            nc_galaxy_shape_factor_moments_tilt.h
 *
 *  Mon Sep 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moments_tilt.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_H_
#define _NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_MOMENTS_TILT (nc_galaxy_shape_factor_moments_tilt_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorMomentsTilt, nc_galaxy_shape_factor_moments_tilt, NC, GALAXY_SHAPE_FACTOR_MOMENTS_TILT, NcGalaxyShapeFactor)

NcGalaxyShapeFactorMomentsTilt *nc_galaxy_shape_factor_moments_tilt_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorMomentsTilt *nc_galaxy_shape_factor_moments_tilt_ref (NcGalaxyShapeFactorMomentsTilt *gsfmt);

void nc_galaxy_shape_factor_moments_tilt_free (NcGalaxyShapeFactorMomentsTilt *gsfmt);
void nc_galaxy_shape_factor_moments_tilt_clear (NcGalaxyShapeFactorMomentsTilt **gsfmt);

guint nc_galaxy_shape_factor_moments_tilt_get_solve_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt);
void nc_galaxy_shape_factor_moments_tilt_reset_solve_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt);
guint nc_galaxy_shape_factor_moments_tilt_get_table_build_count (NcGalaxyShapeFactorMomentsTilt *gsfmt);
guint nc_galaxy_shape_factor_moments_tilt_get_range_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt);
void nc_galaxy_shape_factor_moments_tilt_reset_range_error_count (NcGalaxyShapeFactorMomentsTilt *gsfmt);
void nc_galaxy_shape_factor_moments_tilt_reset_table_build_count (NcGalaxyShapeFactorMomentsTilt *gsfmt);

void nc_galaxy_shape_factor_moments_tilt_exact_moments (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble ghat, gdouble *mu, gdouble *Ex2, gdouble *Ey2);
void nc_galaxy_shape_factor_moments_tilt_eval_tilt (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, gdouble *lambda_1, gdouble *lambda_2, gdouble *lambda_3, gdouble *W);
void nc_galaxy_shape_factor_moments_tilt_peek_layout (NcGalaxyShapeFactorMomentsTilt *gsfmt, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, guint *n_panels, gdouble *top, GArray **degrees);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_MOMENTS_TILT_H_ */


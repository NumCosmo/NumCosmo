/***************************************************************************
 *            nc_galaxy_shape_factor_moments_gauss.h
 *
 *  Mon Sep 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moments_gauss.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_H_
#define _NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS (nc_galaxy_shape_factor_moments_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorMomentsGauss, nc_galaxy_shape_factor_moments_gauss, NC, GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS, NcGalaxyShapeFactor)

NcGalaxyShapeFactorMomentsGauss *nc_galaxy_shape_factor_moments_gauss_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorMomentsGauss *nc_galaxy_shape_factor_moments_gauss_ref (NcGalaxyShapeFactorMomentsGauss *gsfmg);

void nc_galaxy_shape_factor_moments_gauss_free (NcGalaxyShapeFactorMomentsGauss *gsfmg);
void nc_galaxy_shape_factor_moments_gauss_clear (NcGalaxyShapeFactorMomentsGauss **gsfmg);

guint nc_galaxy_shape_factor_moments_gauss_get_table_build_count (NcGalaxyShapeFactorMomentsGauss *gsfmg);
void nc_galaxy_shape_factor_moments_gauss_reset_table_build_count (NcGalaxyShapeFactorMomentsGauss *gsfmg);

void nc_galaxy_shape_factor_moments_gauss_exact_moments (NcGalaxyShapeFactorMomentsGauss *gsfmg, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble ghat, gdouble *mu, gdouble *var_x, gdouble *e_y2);
void nc_galaxy_shape_factor_moments_gauss_eval_moments (NcGalaxyShapeFactorMomentsGauss *gsfmg, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, gdouble *mu, gdouble *var_x, gdouble *e_y2);
void nc_galaxy_shape_factor_moments_gauss_peek_layout (NcGalaxyShapeFactorMomentsGauss *gsfmg, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, guint *n_panels, gdouble *top, GArray **degrees);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_MOMENTS_GAUSS_H_ */


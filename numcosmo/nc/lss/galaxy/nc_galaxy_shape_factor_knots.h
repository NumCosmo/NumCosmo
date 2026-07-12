/***************************************************************************
 *            nc_galaxy_shape_factor_knots.h
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_knots.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_KNOTS_H_
#define _NC_GALAXY_SHAPE_FACTOR_KNOTS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_KNOTS (nc_galaxy_shape_factor_knots_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorKnots, nc_galaxy_shape_factor_knots, NC, GALAXY_SHAPE_FACTOR_KNOTS, NcGalaxyShapeFactor)

/**
 * NcGalaxyShapeFactorKnotsMethod:
 * @NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN: Tensor-product Gauss-Legendre nodes on $[-B,B]^2$, fixed for every galaxy.
 * @NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN_ADAPTIVE: Same tensor-product rule, but re-centered on the joint (population x noise) mode and rescaled to its local curvature at every evaluation, capped at $B$ (see the class documentation).
 *
 * Fixed-node scheme used to lay out the $(u,v)$ knots and weights of the
 * plane-substituted integral evaluated by #NcGalaxyShapeFactorKnots (see the
 * class documentation for the substitution itself). Add new enum values here
 * together with a matching case in `_nc_galaxy_shape_factor_knots_regenerate()`
 * to try alternative schemes; everything else (re-centering, the integrand)
 * is shared.
 */
typedef enum _NcGalaxyShapeFactorKnotsMethod
{
  NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN,
  NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_GL_CARTESIAN_ADAPTIVE,
  /*< private >*/
  NC_GALAXY_SHAPE_FACTOR_KNOTS_METHOD_LEN
} NcGalaxyShapeFactorKnotsMethod;

NcGalaxyShapeFactorKnots *nc_galaxy_shape_factor_knots_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorKnots *nc_galaxy_shape_factor_knots_ref (NcGalaxyShapeFactorKnots *gsfk);

void nc_galaxy_shape_factor_knots_free (NcGalaxyShapeFactorKnots *gsfk);
void nc_galaxy_shape_factor_knots_clear (NcGalaxyShapeFactorKnots **gsfk);

void nc_galaxy_shape_factor_knots_set_bound (NcGalaxyShapeFactorKnots *gsfk, const gdouble bound);
gdouble nc_galaxy_shape_factor_knots_get_bound (NcGalaxyShapeFactorKnots *gsfk);

void nc_galaxy_shape_factor_knots_set_n (NcGalaxyShapeFactorKnots *gsfk, const guint n);
guint nc_galaxy_shape_factor_knots_get_n (NcGalaxyShapeFactorKnots *gsfk);

void nc_galaxy_shape_factor_knots_set_method (NcGalaxyShapeFactorKnots *gsfk, const NcGalaxyShapeFactorKnotsMethod method);
NcGalaxyShapeFactorKnotsMethod nc_galaxy_shape_factor_knots_get_method (NcGalaxyShapeFactorKnots *gsfk);

NcmVector *nc_galaxy_shape_factor_knots_peek_nodes (NcGalaxyShapeFactorKnots *gsfk);
NcmVector *nc_galaxy_shape_factor_knots_peek_weights (NcGalaxyShapeFactorKnots *gsfk);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_KNOTS_H_ */


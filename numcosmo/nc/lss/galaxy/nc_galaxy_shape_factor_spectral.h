/***************************************************************************
 *            nc_galaxy_shape_factor_spectral.h
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_spectral.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_SPECTRAL_H_
#define _NC_GALAXY_SHAPE_FACTOR_SPECTRAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_SPECTRAL (nc_galaxy_shape_factor_spectral_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorSpectral, nc_galaxy_shape_factor_spectral, NC, GALAXY_SHAPE_FACTOR_SPECTRAL, NcGalaxyShapeFactor)

NcGalaxyShapeFactorSpectral *nc_galaxy_shape_factor_spectral_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorSpectral *nc_galaxy_shape_factor_spectral_ref (NcGalaxyShapeFactorSpectral *gsfs);

void nc_galaxy_shape_factor_spectral_free (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_clear (NcGalaxyShapeFactorSpectral **gsfs);

void nc_galaxy_shape_factor_spectral_set_reltol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble reltol);
gdouble nc_galaxy_shape_factor_spectral_get_reltol (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_use_interp (NcGalaxyShapeFactorSpectral *gsfs, const gboolean use_interp);
gboolean nc_galaxy_shape_factor_spectral_get_use_interp (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs, const gdouble g_max);
gdouble nc_galaxy_shape_factor_spectral_get_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol);
gdouble nc_galaxy_shape_factor_spectral_get_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs, const guint n);
guint nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs);

gdouble nc_galaxy_shape_factor_spectral_eval_ln_direct (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise);
NcmMatrix *nc_galaxy_shape_factor_spectral_eval_panels (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, gdouble *ln_marginal);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_SPECTRAL_H_ */


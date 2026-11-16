/***************************************************************************
 *            nc_galaxy_shape_factor_spectral.h
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_spectral.h
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_SPECTRAL_H_
#define _NC_GALAXY_SHAPE_FACTOR_SPECTRAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_factor.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR_SPECTRAL (nc_galaxy_shape_factor_spectral_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapeFactorSpectral, nc_galaxy_shape_factor_spectral, NC, GALAXY_SHAPE_FACTOR_SPECTRAL, NcGalaxyShapeFactor)

/**
 * NcGalaxyShapeFactorSpectralFallback:
 * @NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION: the population is outside the box of the expansion of the galaxy, or its expansion did not converge
 * @NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_RANGE: the shear is beyond the range the interpolant of the galaxy covers
 * @NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_CAP: the shear is above #NcGalaxyShapeFactorSpectral:interp-g-max
 *
 * Why $\ln P$ of a galaxy was computed without its interpolant, with
 * #NcGalaxyShapeFactorSpectral:use-interp set (see
 * nc_galaxy_shape_factor_spectral_get_fallback_stats()).
 */
typedef enum _NcGalaxyShapeFactorSpectralFallback /*< prefix=NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK >*/
{
  NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION,
  NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_RANGE,
  NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_CAP,
  /*< private >*/
  NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_LEN, /*< skip >*/
} NcGalaxyShapeFactorSpectralFallback;

NcGalaxyShapeFactorSpectral *nc_galaxy_shape_factor_spectral_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxyShapeFactorSpectral *nc_galaxy_shape_factor_spectral_ref (NcGalaxyShapeFactorSpectral *gsfs);

void nc_galaxy_shape_factor_spectral_free (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_clear (NcGalaxyShapeFactorSpectral **gsfs);

void nc_galaxy_shape_factor_spectral_set_reltol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble reltol);
gdouble nc_galaxy_shape_factor_spectral_get_reltol (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_obs_source (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyWLObs *obs);

void nc_galaxy_shape_factor_spectral_set_use_interp (NcGalaxyShapeFactorSpectral *gsfs, const gboolean use_interp);
gboolean nc_galaxy_shape_factor_spectral_get_use_interp (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs, const gdouble g_max);
gdouble nc_galaxy_shape_factor_spectral_get_interp_g_max (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol);
gdouble nc_galaxy_shape_factor_spectral_get_interp_abstol (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs, const guint n);
guint nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_line_abstol (NcGalaxyShapeFactorSpectral *gsfs, const gdouble abstol);
gdouble nc_galaxy_shape_factor_spectral_get_line_abstol (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_restrict_range (NcGalaxyShapeFactorSpectral *gsfs, const gboolean restrict_range);
gboolean nc_galaxy_shape_factor_spectral_get_restrict_range (NcGalaxyShapeFactorSpectral *gsfs);
void nc_galaxy_shape_factor_spectral_set_pop_order (NcGalaxyShapeFactorSpectral *gsfs, const guint pop_order);
guint nc_galaxy_shape_factor_spectral_get_pop_order (NcGalaxyShapeFactorSpectral *gsfs);
guint nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (NcGalaxyShapeFactorSpectral *gsfs);
guint nc_galaxy_shape_factor_spectral_get_n_loaded_short (NcGalaxyShapeFactorSpectral *gsfs);

void nc_galaxy_shape_factor_spectral_get_fallback_stats (NcGalaxyShapeFactorSpectralFallback reason, guint64 *lines, guint64 *shears, gdouble *seconds);
void nc_galaxy_shape_factor_spectral_reset_fallback_stats (void);

gdouble nc_galaxy_shape_factor_spectral_eval_ln_direct (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise);
NcmMatrix *nc_galaxy_shape_factor_spectral_eval_panels (NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, gdouble *ln_marginal);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_SPECTRAL_H_ */


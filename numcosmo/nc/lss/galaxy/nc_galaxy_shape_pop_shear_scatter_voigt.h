/***************************************************************************
 *            nc_galaxy_shape_pop_shear_scatter_voigt.h
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_shear_scatter_voigt.h
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

#ifndef _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_H_
#define _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop_shear_scatter.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT (nc_galaxy_shape_pop_shear_scatter_voigt_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapePopShearScatterVoigt, nc_galaxy_shape_pop_shear_scatter_voigt, NC, GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT, NcGalaxyShapePopShearScatter)

/**
 * NcGalaxyShapePopShearScatterVoigtParams:
 * @NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_SIGMA: the width of the Gaussian core in reduced-shear units.
 * @NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_GAMMA: the half width at half maximum of the Lorentzian wings in reduced-shear units.
 *
 * Voigt shear-scatter parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_PARAMS,prefix=NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT >*/
{
  NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_SIGMA = 0,
  NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_GAMMA,
  /* < private > */
  NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_SPARAM_LEN, /*< skip >*/
} NcGalaxyShapePopShearScatterVoigtParams;

#define NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_SIGMA (0.25)
#define NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_GAMMA (0.01)
#define NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxyShapePopShearScatterVoigt *nc_galaxy_shape_pop_shear_scatter_voigt_new (void);
NcGalaxyShapePopShearScatterVoigt *nc_galaxy_shape_pop_shear_scatter_voigt_ref (NcGalaxyShapePopShearScatterVoigt *gspssv);

void nc_galaxy_shape_pop_shear_scatter_voigt_free (NcGalaxyShapePopShearScatterVoigt *gspssv);
void nc_galaxy_shape_pop_shear_scatter_voigt_clear (NcGalaxyShapePopShearScatterVoigt **gspssv);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_H_ */


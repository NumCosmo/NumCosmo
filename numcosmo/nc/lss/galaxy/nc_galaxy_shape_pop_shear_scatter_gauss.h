/***************************************************************************
 *            nc_galaxy_shape_pop_shear_scatter_gauss.h
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_shear_scatter_gauss.h
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

#ifndef _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_H_
#define _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop_shear_scatter.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (nc_galaxy_shape_pop_shear_scatter_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyShapePopShearScatterGauss, nc_galaxy_shape_pop_shear_scatter_gauss, NC, GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS, NcGalaxyShapePopShearScatter)

/**
 * NcGalaxyShapePopShearScatterGaussParams:
 * @NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_SIGMA: the Gaussian width in reduced-shear units.
 *
 * Gaussian shear-scatter parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_PARAMS,prefix=NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS >*/
{
  NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_SIGMA = 0,
  /* < private > */
  NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_SPARAM_LEN, /*< skip >*/
} NcGalaxyShapePopShearScatterGaussParams;

#define NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_DEFAULT_SIGMA (0.25)
#define NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_DEFAULT_PARAMS_ABSTOL (0.0)

#define NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_COL_E_RMS "e_rms"

NcGalaxyShapePopShearScatterGauss *nc_galaxy_shape_pop_shear_scatter_gauss_new (gboolean sigma_from_catalog);
NcGalaxyShapePopShearScatterGauss *nc_galaxy_shape_pop_shear_scatter_gauss_ref (NcGalaxyShapePopShearScatterGauss *gspssg);

void nc_galaxy_shape_pop_shear_scatter_gauss_free (NcGalaxyShapePopShearScatterGauss *gspssg);
void nc_galaxy_shape_pop_shear_scatter_gauss_clear (NcGalaxyShapePopShearScatterGauss **gspssg);

gboolean nc_galaxy_shape_pop_shear_scatter_gauss_get_sigma_from_catalog (NcGalaxyShapePopShearScatterGauss *gspssg);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_H_ */


/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape_gauss.h
 *
 *  Mon June 5 15:22:39 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_gauss.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NC_GALAXY_SD_SHAPE_GAUSS_H_
#define _NC_GALAXY_SD_SHAPE_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_center.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_GAUSS (nc_galaxy_sd_shape_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDShapeGauss, nc_galaxy_sd_shape_gauss, NC, GALAXY_SD_SHAPE_GAUSS, NcGalaxySDShape)

/**
 * NcGalaxySDShapeGaussParams:
 * @NC_GALAXY_SD_SHAPE_GAUSS_SIGMA: Standard deviation
 *
 * Gaussian galaxy shape distribution model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_SHAPE_GAUSS_PARAMS >*/
{
  NC_GALAXY_SD_SHAPE_GAUSS_SIGMA = 0,
  /* < private > */
  NC_GALAXY_SD_SHAPE_GAUSS_SPARAM_LEN, /*< skip >*/
} NcGalaxySDShapeGaussParams;

#define NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_SIGMA  (0.3)

#define NC_GALAXY_SD_SHAPE_GAUSS_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDShapeGauss *nc_galaxy_sd_shape_gauss_new (NcDistance *dist);
NcGalaxySDShapeGauss *nc_galaxy_sd_shape_gauss_ref (NcGalaxySDShapeGauss *gsdsgauss);

void nc_galaxy_sd_shape_gauss_free (NcGalaxySDShapeGauss *gsdsgauss);
void nc_galaxy_sd_shape_gauss_clear (NcGalaxySDShapeGauss **gsdsgauss);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_GAUSS_H_ */


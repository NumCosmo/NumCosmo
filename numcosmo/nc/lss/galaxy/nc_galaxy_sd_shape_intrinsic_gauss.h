/***************************************************************************
 *            nc_galaxy_sd_shape_intrinsic_gauss.h
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_intrinsic_gauss.h
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

#ifndef _NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_H_
#define _NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_shape_intrinsic.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC_GAUSS (nc_galaxy_sd_shape_intrinsic_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDShapeIntrinsicGauss, nc_galaxy_sd_shape_intrinsic_gauss, NC, GALAXY_SD_SHAPE_INTRINSIC_GAUSS, NcGalaxySDShapeIntrinsic)

/**
 * NcGalaxySDShapeIntrinsicGaussParams:
 * @NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_SIGMA: width of the truncated Gaussian.
 *
 * Truncated-Gaussian intrinsic ellipticity model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_PARAMS >*/
{
  NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_SIGMA = 0,
  /* < private > */
  NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_SPARAM_LEN, /*< skip >*/
} NcGalaxySDShapeIntrinsicGaussParams;

#define NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_DEFAULT_SIGMA (0.3)
#define NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDShapeIntrinsicGauss *nc_galaxy_sd_shape_intrinsic_gauss_new (void);
NcGalaxySDShapeIntrinsicGauss *nc_galaxy_sd_shape_intrinsic_gauss_ref (NcGalaxySDShapeIntrinsicGauss *gsig);

void nc_galaxy_sd_shape_intrinsic_gauss_free (NcGalaxySDShapeIntrinsicGauss *gsig);
void nc_galaxy_sd_shape_intrinsic_gauss_clear (NcGalaxySDShapeIntrinsicGauss **gsig);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_H_ */

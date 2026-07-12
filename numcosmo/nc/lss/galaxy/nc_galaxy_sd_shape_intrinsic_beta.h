/***************************************************************************
 *            nc_galaxy_sd_shape_intrinsic_beta.h
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_intrinsic_beta.h
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

#ifndef _NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_H_
#define _NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_shape_intrinsic.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC_BETA (nc_galaxy_sd_shape_intrinsic_beta_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDShapeIntrinsicBeta, nc_galaxy_sd_shape_intrinsic_beta, NC, GALAXY_SD_SHAPE_INTRINSIC_BETA, NcGalaxySDShapeIntrinsic)

/**
 * NcGalaxySDShapeIntrinsicBetaParams:
 * @NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_MU: mean of $x = |\chi_I|^2$ (controls typical ellipticity).
 * @NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_NU: concentration of the Beta distribution.
 *
 * Beta intrinsic ellipticity model parameters, with $\alpha = \mu\nu$ and
 * $\beta = (1-\mu)\nu$.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_PARAMS >*/
{
  NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_MU = 0,
  NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_NU,
  /* < private > */
  NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_SPARAM_LEN, /*< skip >*/
} NcGalaxySDShapeIntrinsicBetaParams;

#define NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_MU (0.18)
#define NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_NU (5.0)
#define NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDShapeIntrinsicBeta *nc_galaxy_sd_shape_intrinsic_beta_new (void);
NcGalaxySDShapeIntrinsicBeta *nc_galaxy_sd_shape_intrinsic_beta_ref (NcGalaxySDShapeIntrinsicBeta *gsib);

void nc_galaxy_sd_shape_intrinsic_beta_free (NcGalaxySDShapeIntrinsicBeta *gsib);
void nc_galaxy_sd_shape_intrinsic_beta_clear (NcGalaxySDShapeIntrinsicBeta **gsib);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_H_ */

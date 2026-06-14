/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape_hsm_gauss_global.h
 *
 *  Mon June 5 15:22:39 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_hsm_gauss_global.h
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_H_
#define _NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (nc_galaxy_sd_shape_hsm_gauss_global_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDShapeHSMGaussGlobal, nc_galaxy_sd_shape_hsm_gauss_global, NC, GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL, NcGalaxySDShape)

/**
 * NcGalaxySDShapeHSMGaussGlobalParams:
 * @NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SIGMA: Standard deviation of the ellipticity distribution.
 *
 * Gaussian galaxy shape distribution model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_PARAMS >*/
{
  NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SIGMA = 0,
  /* < private > */
  NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_SPARAM_LEN, /*< skip >*/
} NcGalaxySDShapeHSMGaussGlobalParams;

#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_DEFAULT_SIGMA   (0.3)
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_DEFAULT_PARAMS_ABSTOL (0.0)

NcGalaxySDShapeHSMGaussGlobal *nc_galaxy_sd_shape_hsm_gauss_global_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxySDShapeHSMGaussGlobal *nc_galaxy_sd_shape_hsm_gauss_global_ref (NcGalaxySDShapeHSMGaussGlobal *gsdsgauss);

void nc_galaxy_sd_shape_hsm_gauss_global_free (NcGalaxySDShapeHSMGaussGlobal *gsdsgauss);
void nc_galaxy_sd_shape_hsm_gauss_global_clear (NcGalaxySDShapeHSMGaussGlobal **gsdsgauss);

void nc_galaxy_sd_shape_hsm_gauss_global_gen (NcGalaxySDShapeHSMGaussGlobal *gsdsgauss, NcmMSet *mset, NcGalaxySDShapeData *data, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcGalaxyWLObsCoord coord, NcmRNG *rng);

void nc_galaxy_sd_shape_hsm_gauss_global_data_set (NcGalaxySDShapeHSMGaussGlobal *gsdsgauss, NcGalaxySDShapeData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m);
void nc_galaxy_sd_shape_hsm_gauss_global_data_get (NcGalaxySDShapeHSMGaussGlobal *gsdsgauss, NcGalaxySDShapeData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m);

#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_EPSILON_OBS_1 "epsilon_obs_1"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_EPSILON_OBS_2 "epsilon_obs_2"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_STD_NOISE "std_noise"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_C1 "c1"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_C2 "c2"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_COL_M "m"

gdouble nc_galaxy_sd_shape_hsm_gauss_global_sigma_from_std_shape (const gdouble std_shape);
gdouble nc_galaxy_sd_shape_hsm_gauss_global_std_shape_from_sigma (const gdouble sigma);

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL_H_ */


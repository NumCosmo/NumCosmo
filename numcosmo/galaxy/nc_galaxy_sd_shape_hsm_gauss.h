
/***************************************************************************
 *            nc_galaxy_sd_shape_hsm_gauss.h
 *
 *  Sun Jan 5 12:53:27 2025
 *  Copyright  2025  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_hsm_gauss.h
 * Copyright (C) 2025 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_GALAXY_SD_SHAPE_HSM_GAUSS_H_
#define _NC_GALAXY_SD_SHAPE_HSM_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_HSM_GAUSS (nc_galaxy_sd_shape_hsm_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDShapeHSMGauss, nc_galaxy_sd_shape_hsm_gauss, NC, GALAXY_SD_SHAPE_HSM_GAUSS, NcGalaxySDShape)

NcGalaxySDShapeHSMGauss *nc_galaxy_sd_shape_hsm_gauss_new (NcGalaxyWLObsEllipConv ellip_conv);
NcGalaxySDShapeHSMGauss *nc_galaxy_sd_shape_hsm_gauss_ref (NcGalaxySDShapeHSMGauss *gsdshsc);

void nc_galaxy_sd_shape_hsm_gauss_free (NcGalaxySDShapeHSMGauss *gsdshsc);
void nc_galaxy_sd_shape_hsm_gauss_clear (NcGalaxySDShapeHSMGauss **gsdshsc);

void nc_galaxy_sd_shape_hsm_gauss_gen (NcGalaxySDShapeHSMGauss *gsdshsc, NcmMSet *mset, NcGalaxySDShapeData *data, const gdouble std_shape, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcGalaxyWLObsCoord coord, NcmRNG *rng);

void nc_galaxy_sd_shape_hsm_gauss_data_set (NcGalaxySDShapeHSMGauss *gsdshsc, NcGalaxySDShapeData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_shape, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m);
void nc_galaxy_sd_shape_hsm_gauss_data_get (NcGalaxySDShapeHSMGauss *gsdshsc, NcGalaxySDShapeData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_shape, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m);

#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_1 "epsilon_obs_1"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_EPSILON_OBS_2 "epsilon_obs_2"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_SHAPE "std_shape"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_STD_NOISE "std_noise"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C1 "c1"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_C2 "c2"
#define NC_GALAXY_SD_SHAPE_HSM_GAUSS_COL_M "m"

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_HSM_GAUSS_H_ */


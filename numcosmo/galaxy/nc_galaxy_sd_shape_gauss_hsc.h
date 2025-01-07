/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape_gauss_hsc.h
 *
 *  Sun Jan 5 12:53:27 2025
 *  Copyright  2025  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_gauss_hsc.h
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

#ifndef _NC_GALAXY_SD_SHAPE_GAUSS_HSC_H_
#define _NC_GALAXY_SD_SHAPE_GAUSS_HSC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/galaxy/nc_galaxy_sd_shape.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE_GAUSS_HSC (nc_galaxy_sd_shape_gauss_hsc_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxySDShapeGaussHSC, nc_galaxy_sd_shape_gauss_hsc, NC, GALAXY_SD_SHAPE_GAUSS_HSC, NcGalaxySDShape)

/**
 * NcGalaxySDShapeGaussHSCParams:
 *
 * Gaussian galaxy shape distribution model parameters.
 *
 */
typedef enum /*< enum,underscore_name=NC_GALAXY_SD_SHAPE_GAUSS_HSC_PARAMS >*/
{
  /* < private > */
  NC_GALAXY_SD_SHAPE_GAUSS_HSC_SPARAM_LEN, /*< skip >*/
} NcGalaxySDShapeGaussHSCParams;


NcGalaxySDShapeGaussHSC *nc_galaxy_sd_shape_gauss_hsc_new ();
NcGalaxySDShapeGaussHSC *nc_galaxy_sd_shape_gauss_hsc_ref (NcGalaxySDShapeGaussHSC *gsdshsc);

void nc_galaxy_sd_shape_gauss_hsc_free (NcGalaxySDShapeGaussHSC *gsdshsc);
void nc_galaxy_sd_shape_gauss_hsc_clear (NcGalaxySDShapeGaussHSC **gsdshsc);

void nc_galaxy_sd_shape_gauss_hsc_gen (NcGalaxySDShapeGaussHSC *gsdshsc, NcmMSet *mset, NcGalaxySDShapeData *data, const gdouble sigma_int, const gdouble sigma_obs, NcGalaxyWLObsCoord coord, NcmRNG *rng);

void nc_galaxy_sd_shape_gauss_hsc_data_set (NcGalaxySDShapeGaussHSC *gsdshsc, NcGalaxySDShapeData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble sigma_int, const gdouble sigma_obs);
void nc_galaxy_sd_shape_gauss_hsc_data_get (NcGalaxySDShapeGaussHSC *gsdshsc, NcGalaxySDShapeData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *sigma_int, gdouble *sigma_obs);

#define NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_1 "epsilon_obs_1"
#define NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_EPSILON_OBS_2 "epsilon_obs_2"
#define NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_SIGMA_INT "sigma_int"
#define NC_GALAXY_SD_SHAPE_GAUSS_HSC_COL_SIGMA_OBS "sigma_obs"

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_GAUSS_HSC_H_ */


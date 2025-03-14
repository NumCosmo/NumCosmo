/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape.h
 *
 *  Sat May 21 22:00:06 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape.h
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
#ifndef _NC_GALAXY_SD_SHAPE_H_
#define _NC_GALAXY_SD_SHAPE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/galaxy/nc_galaxy_sd_position.h>
#include <numcosmo/lss/nc_halo_position.h>
#include <numcosmo/lss/nc_halo_density_profile.h>
#include <numcosmo/lss/nc_wl_surface_mass_density.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/nc_enum_types.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_get_type ())
#define NC_TYPE_GALAXY_SD_SHAPE_DATA (nc_galaxy_sd_shape_data_get_type ())
#define NC_TYPE_GALAXY_SD_SHAPE_INTEGRAND (nc_galaxy_sd_shape_integrand_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxySDShape, nc_galaxy_sd_shape, NC, GALAXY_SD_SHAPE, NcmModel)
typedef struct _NcGalaxySDShapeData NcGalaxySDShapeData;

NCM_UTIL_DECLARE_CALLBACK (NcGalaxySDShapeIntegrand,
                           NC_GALAXY_SD_SHAPE_INTEGRAND,
                           nc_galaxy_sd_shape_integrand,
                           gdouble,
                           NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxySDShapeData * data))


#define NC_GALAXY_SD_SHAPE_DATA(obj) ((NcGalaxySDShapeData *) (obj))

struct _NcGalaxySDShapeClass
{
  /*< private >*/
  NcmModelClass parent_class;

  void (*gen) (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng);
  NcGalaxySDShapeIntegrand *(*integ) (NcGalaxySDShape *gsds);
  gboolean (*prepare_data_array) (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array);
  void (*data_init) (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[10];
};

struct _NcGalaxySDShapeData
{
  NcGalaxySDPositionData *sdpos_data;
  NcGalaxyWLObsCoord coord;
  gdouble epsilon_int_1;
  gdouble epsilon_int_2;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxySDShapeData *data, GList *columns);
  gdouble (*ldata_get_radius) (NcGalaxySDShapeData *data);
  gatomicrefcount ref_count;
};

NCM_MSET_MODEL_DECLARE_ID (nc_galaxy_sd_shape);

GType nc_galaxy_sd_shape_data_get_type (void) G_GNUC_CONST;

NcGalaxySDShapeData *nc_galaxy_sd_shape_data_ref (NcGalaxySDShapeData *data);
void nc_galaxy_sd_shape_data_unref (NcGalaxySDShapeData *data);

void nc_galaxy_sd_shape_data_read_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_sd_shape_data_write_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_sd_shape_data_required_columns (NcGalaxySDShapeData *data);
gdouble nc_galaxy_sd_shape_data_get_radius (NcGalaxySDShapeData *data);

NcGalaxySDShape *nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds);

void nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds);
void nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds);

void nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng);
NcGalaxySDShapeIntegrand *nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds);
gboolean nc_galaxy_sd_shape_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array);

NcGalaxySDShapeData *nc_galaxy_sd_shape_data_new (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data);

#define NC_GALAXY_SD_SHAPE_COL_COORD "coord"
#define NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1 "epsilon_int_1"
#define NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2 "epsilon_int_2"

G_END_DECLS

#endif /* _NC_GALAXY_SD_SHAPE_H_ */


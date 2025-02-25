/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape.c
 *
 *  Sat May 21 20:43:32 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape.c
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

/**
 * NcGalaxySDShape:
 *
 * Class describing galaxy sample shape distribution.
 *
 * This class describes a galaxy sample shape distribution. It is composed by a
 * distribution $P(s)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_wl_obs.h"
#include "galaxy/nc_galaxy_sd_shape.h"
#include "nc_enum_types.h"
#include "nc_hicosmo.h"
#include "lss/nc_halo_position.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"

typedef struct _NcGalaxySDShapePrivate
{
  gint placeholder;
} NcGalaxySDShapePrivate;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShape, nc_galaxy_sd_shape, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxySDShapeData, nc_galaxy_sd_shape_data, nc_galaxy_sd_shape_data_ref, nc_galaxy_sd_shape_data_unref);
NCM_UTIL_DEFINE_CALLBACK (NcGalaxySDShapeIntegrand,
                          NC_GALAXY_SD_SHAPE_INTEGRAND,
                          nc_galaxy_sd_shape_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxySDShapeData * data),
                          NCM_UTIL_CALLBACK_ARGS (z, data))

static void
nc_galaxy_sd_shape_init (NcGalaxySDShape *gsds)
{
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE (gsds));

  switch (property_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_shape_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (object);

  g_return_if_fail (NC_IS_GALAXY_SD_SHAPE (gsds));

  switch (property_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_shape_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_shape, NC_TYPE_GALAXY_SD_SHAPE);

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_shape_gen: method not implemented.");
}

static NcGalaxySDShapeIntegrand *
_nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds)
{
  g_error ("_nc_galaxy_sd_shape_integ: method not implemented.");

  return NULL;
}

static gboolean
_nc_galaxy_sd_shape_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array)
{
  g_error ("_nc_galaxy_sd_shape_prepare: method not implemented.");

  return FALSE;
}

static void
_nc_galaxy_sd_shape_data_init (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data, NcGalaxySDShapeData *data)
{
  g_error ("_nc_galaxy_sd_shape_data_new: method not implemented.");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_shape_class_init (NcGalaxySDShapeClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_shape_set_property;
  model_class->get_property = &_nc_galaxy_sd_shape_get_property;
  object_class->finalize    = &_nc_galaxy_sd_shape_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample shape distribution", "GalaxySDShape");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxySDShape", "Galaxy sample shape distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->gen                = &_nc_galaxy_sd_shape_gen;
  klass->integ              = &_nc_galaxy_sd_shape_integ;
  klass->prepare_data_array = &_nc_galaxy_sd_shape_prepare_data_array;
  klass->data_init          = &_nc_galaxy_sd_shape_data_init;
}

/**
 * nc_galaxy_sd_shape_data_ref:
 * @data: a #NcGalaxySDShapeData
 *
 * Increases the reference count of @data by one.
 *
 */
NcGalaxySDShapeData *
nc_galaxy_sd_shape_data_ref (NcGalaxySDShapeData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_sd_shape_data_unref:
 * @data: a #NcGalaxySDShapeData
 *
 * Decreases the reference count of @data by one. If the reference count reaches 0, the
 * data is freed.
 *
 */
void
nc_galaxy_sd_shape_data_unref (NcGalaxySDShapeData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    nc_galaxy_sd_position_data_unref (data->sdpos_data);
    g_free (data);
  }
}

/**
 * nc_galaxy_sd_shape_data_read_row:
 * @data: a #NcGalaxySDShapeData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the row @i from the galaxy shape data.
 *
 */
void
nc_galaxy_sd_shape_data_read_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_sd_position_data_read_row (data->sdpos_data, obs, i);
  {
    data->coord         = nc_galaxy_wl_obs_get_coord (obs);
    data->epsilon_int_1 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i);
    data->epsilon_int_2 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i);

    data->ldata_read_row (data, obs, i);
  }
}

/**
 * nc_galaxy_sd_shape_data_write_row:
 * @data: a #NcGalaxySDShapeData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the row @i to the galaxy shape data.
 *
 */
void
nc_galaxy_sd_shape_data_write_row (NcGalaxySDShapeData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_sd_position_data_write_row (data->sdpos_data, obs, i);
  {
    nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1, i, data->epsilon_int_1);
    nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2, i, data->epsilon_int_2);

    data->ldata_write_row (data, obs, i);
  }
}

/**
 * nc_galaxy_sd_shape_data_required_columns:
 * @data: a #NcGalaxySDShapeData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy shape data.
 */
GList *
nc_galaxy_sd_shape_data_required_columns (NcGalaxySDShapeData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_SHAPE_COL_EPSILON_INT_2));
  data->ldata_required_columns (data, columns);

  {
    GList *sdpos_columns = nc_galaxy_sd_position_data_required_columns (data->sdpos_data);

    columns = g_list_concat (columns, sdpos_columns);
  }

  return columns;
}

/**
 * nc_galaxy_sd_shape_data_get_radius:
 * @data: a #NcGalaxySDShapeData
 *
 * Returns: the radius of the galaxy shape data.
 */
gdouble
nc_galaxy_sd_shape_data_get_radius (NcGalaxySDShapeData *data)
{
  return data->ldata_get_radius (data);
}

/**
 * nc_galaxy_sd_shape_integrand_new:
 * @func: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxySDShapeIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new galaxy shape integrand.
 *
 * Returns: (transfer full): a new NcGalaxySDShapeIntegrand object.
 */
/**
 * nc_galaxy_sd_shape_integrand_copy:
 * @callback_obj: a NcGalaxySDShapeIntegrand
 *
 * Copies the integrand for the galaxy shape data.
 *
 * Returns: (transfer full): a copy of @callback_obj
 */
/**
 * nc_galaxy_sd_shape_integrand_free:
 * @callback_obj: a NcGalaxySDShapeIntegrand
 *
 * Frees the integrand for the galaxy shape data.
 *
 */
/**
 * nc_galaxy_sd_shape_integrand_prepare:
 * @callback_obj: a NcGalaxySDShapeIntegrand
 * @mset: a #NcmMSet
 *
 * Prepares the integrand for the galaxy shape data.
 *
 */

/**
 * nc_galaxy_sd_shape_ref:
 * @gsds: a #NcGalaxySDShape
 *
 * Increases the reference count of @gsds by one.
 *
 * Returns: (transfer full): @gsds.
 */
NcGalaxySDShape *
nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds)
{
  return g_object_ref (gsds);
}

/**
 * nc_galaxy_sd_shape_free:
 * @gsds: a #NcGalaxySDShape
 *
 * Decreases the reference count of @gsds by one.
 *
 */
void
nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds)
{
  g_object_unref (gsds);
}

/**
 * nc_galaxy_sd_shape_clear:
 * @gsds: a #NcGalaxySDShape
 *
 * Decreases the reference count of @gsds by one, and sets the pointer *@gsds to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds)
{
  g_clear_object (gsds);
}

/**
 * nc_galaxy_sd_shape_gen: (virtual gen)
 * @gsds: a #NcGalaxySDShape
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDShapeData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy shape. The #NcGalaxySDShapeData object @data must be
 * initialized before calling this method.
 *
 */
void
nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcmMSet *mset, NcGalaxySDShapeData *data, NcmRNG *rng)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->gen (gsds, mset, data, rng);
}

/**
 * nc_galaxy_sd_shape_integ: (virtual integ)
 * @gsds: a #NcGalaxySDShape
 *
 * Creates a new galaxy shape integrand.
 *
 * Returns: (transfer full): a new NcGalaxySDShapeIntegrand object.
 */
NcGalaxySDShapeIntegrand *
nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ (gsds);
}

/**
 * nc_galaxy_sd_shape_prepare_data_array: (virtual prepare_data_array)
 * @gsds: a #NcGalaxySDShape
 * @mset: a #NcmMSet
 * @data_array: (element-type NcGalaxySDShapeData): a #GPtrArray of #NcGalaxySDShapeData
 *
 * Prepares the matrix to compute the probability density of the observaple shape.
 *
 * Returns: TRUE if the matrix was prepared, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_shape_prepare_data_array (NcGalaxySDShape *gsds, NcmMSet *mset, GPtrArray *data_array)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->prepare_data_array (gsds, mset, data_array);
}

/**
 * nc_galaxy_sd_shape_data_new:
 * @gsds: a #NcGalaxySDShape
 * @sdpos_data: a #NcGalaxySDPositionData
 *
 * Creates a new galaxy shape data.
 *
 * Returns: (transfer full): a new #NcGalaxySDShapeData object.
 */
NcGalaxySDShapeData *
nc_galaxy_sd_shape_data_new (NcGalaxySDShape *gsds, NcGalaxySDPositionData *sdpos_data)
{
  NcGalaxySDShapeData *data = g_new0 (NcGalaxySDShapeData, 1);

  data->sdpos_data             = nc_galaxy_sd_position_data_ref (sdpos_data);
  data->coord                  = NC_GALAXY_WL_OBS_COORD_CELESTIAL;
  data->epsilon_int_1          = 0.0;
  data->epsilon_int_2          = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->data_init (gsds, sdpos_data, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}


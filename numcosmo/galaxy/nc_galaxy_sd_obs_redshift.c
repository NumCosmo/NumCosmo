/***************************************************************************
 *            nc_galaxy_sd_obs_redshift.c
 *
 *  Thu Aug 1 00:45:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION: nc_galaxy_sd_obs_redshift
 * @title: NcGalaxySDObsRedshift
 * @short_description: Class describing galaxy sample observed redshift distribution.
 * @stability: Unstable
 *
 *
 * This class describes galaxy sample observed redshift distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxySDObsRedshiftPrivate
{
  gint placeholder;
} NcGalaxySDObsRedshiftPrivate;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshift, nc_galaxy_sd_obs_redshift, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxySDObsRedshiftData, nc_galaxy_sd_obs_redshift_data, nc_galaxy_sd_obs_redshift_data_ref, nc_galaxy_sd_obs_redshift_data_unref);
NCM_UTIL_DEFINE_CALLBACK (NcGalaxySDObsRedshiftIntegrand,
                          NC_GALAXY_SD_OBS_REDSHIFT_INTEGRAND,
                          nc_galaxy_sd_obs_redshift_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxySDObsRedshiftData * data),
                          NCM_UTIL_CALLBACK_ARGS (z, data))

static void
nc_galaxy_sd_obs_redshift_init (NcGalaxySDObsRedshift *gsdor)
{
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_obs_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshift *gsdor = NC_GALAXY_SD_OBS_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_obs_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshift *gsdor = NC_GALAXY_SD_OBS_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_obs_redshift_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_obs_redshift, NC_TYPE_GALAXY_SD_OBS_REDSHIFT)

/*  LCOV_EXCL_START */
static void
_nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_obs_redshift_gen: method not implemented");
}

static void
_nc_galaxy_sd_obs_redshift_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  g_error ("_nc_galaxy_sd_obs_redshift_prepare: method not implemented");
}

static NcGalaxySDObsRedshiftIntegrand *
_nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor)
{
  g_error ("_nc_galaxy_sd_obs_redshift_integ: method not implemented");

  return NULL;
}

static void
_nc_galaxy_sd_obs_redshift_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  g_error ("_nc_galaxy_sd_obs_redshift_data_new: method not implemented");
}

/*  LCOV_EXCL_STOP */

static void
nc_galaxy_sd_obs_redshift_class_init (NcGalaxySDObsRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_get_property;
  object_class->finalize    = &_nc_galaxy_sd_obs_redshift_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample observed redshift distribution", "GalaxySDObsRedshift");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxySDObsRedshift", "Galaxy sample observed redshift distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->gen       = &_nc_galaxy_sd_obs_redshift_gen;
  klass->prepare   = &_nc_galaxy_sd_obs_redshift_prepare;
  klass->integ     = &_nc_galaxy_sd_obs_redshift_integ;
  klass->data_init = &_nc_galaxy_sd_obs_redshift_data_init;
}

/**
 * nc_galaxy_sd_obs_redshift_data_ref:
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Increments the reference count of @data by one.
 *
 * Returns: (transfer full): a copy of @data
 */
NcGalaxySDObsRedshiftData *
nc_galaxy_sd_obs_redshift_data_ref (NcGalaxySDObsRedshiftData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_sd_obs_redshift_data_unref:
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Decreases the reference count of @data by one. If the reference count reaches 0, the
 * data is freed.
 */
void
nc_galaxy_sd_obs_redshift_data_unref (NcGalaxySDObsRedshiftData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_sd_obs_redshift_data_read_row:
 * @data: a #NcGalaxySDObsRedshiftData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the galaxy redshift data from the observation.
 *
 */
void
nc_galaxy_sd_obs_redshift_data_read_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->z = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i);
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_sd_obs_redshift_data_write_row:
 * @data: a #NcGalaxySDObsRedshiftData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the galaxy redshift data to the observation.
 *
 */
void
nc_galaxy_sd_obs_redshift_data_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i, data->z);
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_sd_obs_redshift_data_required_columns:
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy redshift data.
 */
GList *
nc_galaxy_sd_obs_redshift_data_required_columns (NcGalaxySDObsRedshiftData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_OBS_REDSHIFT_COL_Z));
  data->ldata_required_columns (data, columns);

  return columns;
}

/**
 * nc_galaxy_sd_obs_redshift_integrand_new:
 * @func: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new integrand for the galaxy redshift data.
 * The integrand is a function that takes the redshift @z and the galaxy redshift data @data as arguments.
 * The function should return the integrand value at @z.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftIntegrand object.
 */
/**
 * nc_galaxy_sd_obs_redshift_integrand_copy:
 * @callback_obj: a #NcGalaxySDObsRedshiftIntegrand
 *
 * Copies the integrand for the galaxy redshift data.
 *
 * Returns: (transfer full): a copy of @callback_obj
 */
/**
 * nc_galaxy_sd_obs_redshift_integrand_free:
 * @callback_obj: a #NcGalaxySDObsRedshiftIntegrand
 *
 * Frees the integrand for the galaxy redshift data.
 *
 */
/**
 * nc_galaxy_sd_obs_redshift_integrand_prepare:
 * @callback_obj: a #NcGalaxySDObsRedshiftIntegrand
 * @mset: a #NcmMSet
 *
 * Prepares the integrand for the galaxy redshift data.
 *
 */

/**
 * nc_galaxy_sd_obs_redshift_ref:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Increases the reference count of @gsdor by one.
 *
 * Returns: (transfer full): @gsdor
 */
NcGalaxySDObsRedshift *
nc_galaxy_sd_obs_redshift_ref (NcGalaxySDObsRedshift *gsdor)
{
  return g_object_ref (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_free:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Decreases the reference count of @gsdor by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_free (NcGalaxySDObsRedshift *gsdor)
{
  g_object_unref (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_clear:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Decreases the reference count of @gsdor by one, and sets the pointer *@gsdor to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_clear (NcGalaxySDObsRedshift **gsdor)
{
  g_clear_object (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_gen:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @data: a #NcGalaxySDObsRedshiftData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy redshift data. The #NcGalaxySDObsRedshiftData object @data must be
 * initialized before calling this method.
 *
 */
void
nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->gen (gsdor, data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_prepare:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Prepares the galaxy redshift data for generation.
 *
 */
void
nc_galaxy_sd_obs_redshift_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->prepare (gsdor, data);
}


/**
 * nc_galaxy_sd_obs_redshift_integ:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Prepares the integrand for the galaxy redshift data.
 *
 */
NcGalaxySDObsRedshiftIntegrand *
nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->integ (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_data_new:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Creates a new galaxy redshift data.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftData object.
 */
NcGalaxySDObsRedshiftData *
nc_galaxy_sd_obs_redshift_data_new (NcGalaxySDObsRedshift *gsdor)
{
  NcGalaxySDObsRedshiftData *data = g_new0 (NcGalaxySDObsRedshiftData, 1);

  data->z                      = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->data_init (gsdor, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}


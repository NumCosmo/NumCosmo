/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position.c
 *
 *  Sat May 20 17:52:48 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position.c
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
 * NcGalaxySDPosition:
 *
 * Class describing galaxy sample position distributions.
 *
 * This class describes a galaxy sample position distributions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_model.h"
#include "math/ncm_mset.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"


typedef struct _NcGalaxySDPositionPrivate
{
  gint placeholder;
} NcGalaxySDPositionPrivate;

enum
{
  PROP_0,
  PROP_RA_LIM,
  PROP_DEC_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPosition, nc_galaxy_sd_position, NCM_TYPE_MODEL)
G_DEFINE_BOXED_TYPE (NcGalaxySDPositionData, nc_galaxy_sd_position_data, nc_galaxy_sd_position_data_ref, nc_galaxy_sd_position_data_unref); /* LCOV_EXCL_LINE */
NCM_UTIL_DEFINE_CALLBACK (NcGalaxySDPositionIntegrand,
                          NC_GALAXY_SD_POSITION_INTEGRAND,
                          nc_galaxy_sd_position_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (NcGalaxySDPositionData * data),
                          NCM_UTIL_CALLBACK_ARGS (data))

static void
nc_galaxy_sd_position_init (NcGalaxySDPosition *gsdp)
{
}

static void
_nc_galaxy_sd_position_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPosition *gsdp = NC_GALAXY_SD_POSITION (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION (object));

  switch (prop_id)
  {
    case PROP_RA_LIM:
    {
      NcmDTuple2 *ra_lim = g_value_get_boxed (value);

      if (ra_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: ra_lim is NULL.");

      nc_galaxy_sd_position_set_ra_lim (gsdp, ra_lim->elements[0], ra_lim->elements[1]);
      break;
    }
    case PROP_DEC_LIM:
    {
      NcmDTuple2 *dec_lim = g_value_get_boxed (value);

      if (dec_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: dec_lim is NULL.");

      nc_galaxy_sd_position_set_dec_lim (gsdp, dec_lim->elements[0], dec_lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_position_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPosition *gsdp = NC_GALAXY_SD_POSITION (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION (object));

  switch (prop_id)
  {
    case PROP_RA_LIM:
    {
      gdouble ra_min, ra_max;

      nc_galaxy_sd_position_get_ra_lim (gsdp, &ra_min, &ra_max);

      g_value_take_boxed (value, ncm_dtuple2_new (ra_min, ra_max));
      break;
    }
    case PROP_DEC_LIM:
    {
      gdouble dec_min, dec_max;

      nc_galaxy_sd_position_get_dec_lim (gsdp, &dec_min, &dec_max);

      g_value_take_boxed (value, ncm_dtuple2_new (dec_min, dec_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_position_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_position_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_position, NC_TYPE_GALAXY_SD_POSITION)

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_position_gen (NcGalaxySDPosition *gsdp, NcGalaxySDPositionData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_position_gen: method not implemented.");
}

static NcGalaxySDPositionIntegrand *
_nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp)
{
  g_error ("_nc_galaxy_sd_position_integ: method not implemented.");

  return NULL;
}

static gboolean
_nc_galaxy_sd_position_set_ra_lim (NcGalaxySDPosition *gsdp, const gdouble ra_min, const gdouble ra_max)
{
  g_error ("_nc_galaxy_sd_position_set_ra_lim: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_set_dec_lim (NcGalaxySDPosition *gsdp, const gdouble dec_min, const gdouble dec_max)
{
  g_error ("_nc_galaxy_sd_position_set_dec_lim: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max)
{
  g_error ("_nc_galaxy_sd_position_get_ra_lim: method not implemented.");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_position_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max)
{
  g_error ("_nc_galaxy_sd_position_get_dec_lim: method not implemented.");

  return FALSE;
}

static void
_nc_galaxy_sd_position_data_init (NcGalaxySDPosition *gsdp, NcGalaxySDObsRedshiftData *sdz_data, NcGalaxySDPositionData *data)
{
  g_error ("_nc_galaxy_sd_position_data_new: method not implemented.");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_position_class_init (NcGalaxySDPositionClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_position_set_property;
  model_class->get_property = &_nc_galaxy_sd_position_get_property;
  object_class->finalize    = &_nc_galaxy_sd_position_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample position distribution", "NcGalaxySDPosition");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDPosition:ra-lim:
   *
   * Galaxy sample righ ascension distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RA_LIM,
                                   g_param_spec_boxed ("ra-lim",
                                                       NULL,
                                                       "Galaxy sample right ascension distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPosition:dec-lim:
   *
   * Galaxy sample declination distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEC_LIM,
                                   g_param_spec_boxed ("dec-lim",
                                                       NULL,
                                                       "Galaxy sample declination distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class,
                              "NcGalaxySDPosition",
                              "Galaxy sample position distribution.",
                              NULL,
                              TRUE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  klass->gen         = &_nc_galaxy_sd_position_gen;
  klass->integ       = &_nc_galaxy_sd_position_integ;
  klass->set_ra_lim  = &_nc_galaxy_sd_position_set_ra_lim;
  klass->set_dec_lim = &_nc_galaxy_sd_position_set_dec_lim;
  klass->get_ra_lim  = &_nc_galaxy_sd_position_get_ra_lim;
  klass->get_dec_lim = &_nc_galaxy_sd_position_get_dec_lim;
  klass->data_init   = &_nc_galaxy_sd_position_data_init;
}

/**
 * nc_galaxy_sd_position_data_ref:
 * @data: a #NcGalaxySDPositionData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): a copy of @data
 */
NcGalaxySDPositionData *
nc_galaxy_sd_position_data_ref (NcGalaxySDPositionData *data)
{
  g_return_val_if_fail (data != NULL, NULL);

  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_sd_position_data_unref:
 * @data: a #NcGalaxySDPositionData
 *
 * Decreases the reference count of @data by one. If the reference count reaches 0, the
 * data is freed.
 *
 */
void
nc_galaxy_sd_position_data_unref (NcGalaxySDPositionData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    nc_galaxy_sd_obs_redshift_data_unref (data->sdz_data);
    g_free (data);
  }
}

/**
 * nc_galaxy_sd_position_data_read_row:
 * @data: a #NcGalaxySDPositionData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the row @i from the galaxy position data.
 *
 */
void
nc_galaxy_sd_position_data_read_row (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_sd_obs_redshift_data_read_row (data->sdz_data, obs, i);
  {
    data->ra  = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_RA, i);
    data->dec = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_POSITION_COL_DEC, i);
    data->ldata_read_row (data, obs, i);
  }
}

/**
 * nc_galaxy_sd_position_data_write_row:
 * @data: a #NcGalaxySDPositionData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the galaxy position data to the observation.
 *
 */
void
nc_galaxy_sd_position_data_write_row (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_sd_obs_redshift_data_write_row (data->sdz_data, obs, i);
  {
    nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_POSITION_COL_RA, i, data->ra);
    nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_POSITION_COL_DEC, i, data->dec);
    data->ldata_write_row (data, obs, i);
  }
}

/**
 * nc_galaxy_sd_position_data_required_columns:
 * @data: a #NcGalaxySDPositionData
 *
 * Gets the required columns for the galaxy position data.
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy position data.
 */
GList *
nc_galaxy_sd_position_data_required_columns (NcGalaxySDPositionData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_POSITION_COL_RA));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_POSITION_COL_DEC));
  data->ldata_required_columns (data, columns);

  {
    GList *sdz_columns = nc_galaxy_sd_obs_redshift_data_required_columns (data->sdz_data);

    columns = g_list_concat (columns, sdz_columns);
  }

  return columns;
}

/**
 * nc_galaxy_sd_position_integrand_new:
 * @func: (scope async) (closure callback_data): a #NcGalaxySDPositionIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxySDPositionIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxySDPositionIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxySDPositionIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new galaxy position integrand.
 *
 */
/**
 * nc_galaxy_sd_position_integrand_copy:
 * @callback_obj: a NcGalaxySDPositionIntegrand
 *
 * Copies the galaxy position integrand.
 *
 * Returns: (transfer full): a copy of @callback_obj
 */
/**
 * nc_galaxy_sd_position_integrand_free:
 * @callback_obj: a NcGalaxySDPositionIntegrand
 *
 * Frees the galaxy position integrand.
 *
 */
/**
 * nc_galaxy_sd_position_integrand_prepare:
 * @callback_obj: a NcGalaxySDPositionIntegrand
 * @mset: a #NcmMSet
 *
 * Prepares the galaxy position integrand.
 *
 */

/**
 * nc_galaxy_sd_position_ref:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Increases the reference count of @gsdp by one.
 *
 * Returns: (transfer full): @gsdp.
 */
NcGalaxySDPosition *
nc_galaxy_sd_position_ref (NcGalaxySDPosition *gsdp)
{
  return g_object_ref (gsdp);
}

/**
 * nc_galaxy_sd_position_free:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Decreases the reference count of @gsdp by one.
 *
 */
void
nc_galaxy_sd_position_free (NcGalaxySDPosition *gsdp)
{
  g_object_unref (gsdp);
}

/**
 * nc_galaxy_sd_position_clear:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Decreases the reference count of @gsdp by one, and sets the pointer *@gsdp to
 * NULL.
 *
 */
void
nc_galaxy_sd_position_clear (NcGalaxySDPosition **gsdp)
{
  g_clear_object (gsdp);
}

/**
 * nc_galaxy_sd_position_set_ra_lim: (virtual set_ra_lim)
 * @gsdp: a #NcGalaxySDPosition
 * @ra_min: the minimum right ascension
 * @ra_max: the maximum right ascension
 *
 * Sets the right ascension limits of the distribution.
 *
 * Returns: TRUE if the limits were set successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_set_ra_lim (NcGalaxySDPosition *gsdp, const gdouble ra_min, const gdouble ra_max)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_ra_lim (gsdp, ra_min, ra_max);
}

/**
 * nc_galaxy_sd_position_set_dec_lim: (virtual set_dec_lim)
 * @gsdp: a #NcGalaxySDPosition
 * @dec_min: the minimum declination
 * @dec_max: the maximum declination
 *
 * Sets the declination limits of the distribution.
 *
 * Returns: TRUE if the limits were set successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_set_dec_lim (NcGalaxySDPosition *gsdp, const gdouble dec_min, const gdouble dec_max)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->set_dec_lim (gsdp, dec_min, dec_max);
}

/**
 * nc_galaxy_sd_position_get_ra_lim: (virtual get_ra_lim)
 * @gsdp: a #NcGalaxySDPosition
 * @ra_min: (out): the minimum right ascension
 * @ra_max: (out): the maximum right ascension
 *
 * Gets the right ascension limits of the distribution.
 *
 * Returns: TRUE if the limits were retrieved successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_ra_lim (gsdp, ra_min, ra_max);
}

/**
 * nc_galaxy_sd_position_get_dec_lim: (virtual get_dec_lim)
 * @gsdp: a #NcGalaxySDPosition
 * @dec_min: (out): the minimum declination
 * @dec_max: (out): the maximum declination
 *
 * Gets the declination limits of the distribution.
 *
 * Returns: TRUE if the limits were retrieved successfully, FALSE otherwise.
 */
gboolean
nc_galaxy_sd_position_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->get_dec_lim (gsdp, dec_min, dec_max);
}

/**
 * nc_galaxy_sd_position_gen: (virtual gen)
 * @gsdp: a #NcGalaxySDPosition
 * @data: a #NcGalaxySDPositionData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy position. The #NcGalaxySDPositionData object @data must be
 * initialized before calling this method.
 *
 */
void
nc_galaxy_sd_position_gen (NcGalaxySDPosition *gsdp, NcGalaxySDPositionData *data, NcmRNG *rng)
{
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->gen (gsdp, data, rng);
}

/**
 * nc_galaxy_sd_position_integ: (virtual integ)
 * @gsdp: a #NcGalaxySDPosition
 *
 * Prepares the integrand for the galaxy position distribution.
 *
 * Returns: (transfer full): a new NcGalaxySDPositionIntegrand object.
 */
NcGalaxySDPositionIntegrand *
nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp)
{
  return NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->integ (gsdp);
}

/**
 * nc_galaxy_sd_position_data_new:
 * @gsdp: a #NcGalaxySDPosition
 * @sdz_data: a #NcGalaxySDObsRedshiftData
 *
 * Creates a new galaxy position data.
 *
 * Returns: (transfer full): a new #NcGalaxySDPositionData object.
 */
NcGalaxySDPositionData *
nc_galaxy_sd_position_data_new (NcGalaxySDPosition *gsdp, NcGalaxySDObsRedshiftData *sdz_data)
{
  NcGalaxySDPositionData *data = g_new0 (NcGalaxySDPositionData, 1);

  data->sdz_data               = nc_galaxy_sd_obs_redshift_data_ref (sdz_data);
  data->ra                     = 0.0;
  data->dec                    = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SD_POSITION_GET_CLASS (gsdp)->data_init (gsdp, sdz_data, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}


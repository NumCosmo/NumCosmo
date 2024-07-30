/***************************************************************************
 *            nc_galaxy_wl_obs.c
 *
 *  Tue Jul 16 06:43:45 2024
 *  Copyright  2024 Caio Lima de Oliveira
 *  <caiooliveiracode@proton.me>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiooliveiracode@proton.me>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION: nc_galaxy_wl_obs
 * @title: NcGalaxyWLObs
 * @short_description: Galaxy weak lensing observation data.
 *
 * A class to store galaxy weak lensing observation data and information
 * about its coordinate system.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "stdarg.h"
#include "galaxy/nc_galaxy_wl_obs.h"
#include "math/ncm_matrix.h"
#include "math/ncm_obj_array.h"
#include "math/ncm_spline.h"
#include "math/ncm_vector.h"
#include "nc_enum_types.h"

struct _NcGalaxyWLObsPrivate
{
  NcmObjArray *data;
  NcmVarDict *header;
  NcGalaxyWLObsCoord coord;
  gdouble len;
};

enum
{
  PROP_0,
  PROP_DATA,
  PROP_HEADER,
  PROP_COORD,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLObs, nc_galaxy_wl_obs, G_TYPE_OBJECT)

static void
nc_galaxy_wl_obs_init (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = obs->priv = nc_galaxy_wl_obs_get_instance_private (obs);

  self->data   = NULL;
  self->header = NULL;
  self->coord  = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;
  self->len    = 0;
}

static void
nc_galaxy_wl_obs_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs         = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate * const self = obs->priv;

  switch (prop_id)
  {
    case PROP_DATA:
      g_value_set_object (value, self->data);
      break;
    case PROP_HEADER:
      g_value_set_object (value, nc_galaxy_wl_obs_peek_header (obs));
      break;
    case PROP_COORD:
      g_value_set_enum (value, nc_galaxy_wl_obs_get_coord (obs));
      break;
    case PROP_LEN:
      g_value_set_double (value, nc_galaxy_wl_obs_len (obs));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_galaxy_wl_obs_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs         = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate *priv = obs->priv;

  switch (prop_id)
  {
    case PROP_COORD:
      nc_galaxy_wl_obs_set_coord (obs, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_galaxy_wl_obs_dispose (GObject *object)
{
  NcGalaxyWLObs *obs                = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate * const self = obs->priv;

  ncm_var_dict_clear (&self->header);

  if (self->data)
    ncm_obj_array_clear (&self->data);

  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->dispose (object);
}

/**
 * nc_galaxy_wl_obs_finalize:
 * @object: a #GObject.
 *
 * Finalizes the #NcGalaxyWLObs object.
 *
 */
static void
nc_galaxy_wl_obs_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_obs_class_init (NcGalaxyWLObsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->get_property = nc_galaxy_wl_obs_get_property;
  object_class->set_property = nc_galaxy_wl_obs_set_property;
  object_class->dispose      = nc_galaxy_wl_obs_dispose;
  object_class->finalize     = nc_galaxy_wl_obs_finalize;

  /**
   * NcGalaxyWLObs:data:
   *
   * Weak lensing observation data.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DATA,
                                   g_param_spec_boxed ("data",
                                                       "Data",
                                                       "Weak lensing observation data",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));

  /**
   * NcGalaxyWLObs:header:
   *
   * Data columns headers.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_HEADER,
                                   g_param_spec_boxed ("header",
                                                       "Header",
                                                       "Data columns headers",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));

  /**
   * NcGalaxyWLObs:coord:
   *
   * Coordinate system used to store the data.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_COORD,
                                   g_param_spec_enum ("coord",
                                                      "Coordinate system",
                                                      "Coordinate system used to store the data",
                                                      NC_TYPE_GALAXY_WL_OBS_COORD,
                                                      NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));
}

/**
 * nc_galaxy_wl_obs_new:
 * @coord: the coordinate system used to store the data.
 * @nrows: the number of data rows.
 * @...: the data columns names, must end with NULL.
 *
 * Creates a new #NcGalaxyWLObs object.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new (NcGalaxyWLObsCoord coord, gdouble nrows, ...)
{
  NcmVarDict *header = ncm_var_dict_new ();
  NcmObjArray *data  = ncm_obj_array_sized_new (nrows);
  gboolean has_pz    = FALSE;
  gint i             = 0;
  NcGalaxyWLObs *obs;
  va_list args;
  gint j;

  va_start (args, nrows);

  while (TRUE)
  {
    const gchar *key = va_arg (args, const gchar *);

    if (!key)
      break;

    if (g_strcmp0 (key, "pz") == 0)
    {
      has_pz = TRUE;
    }
    else
    {
      ncm_var_dict_set_int (header, key, i);
      i++;
    }
  }

  va_end (args);

  for (j = 0; j < nrows; j++)
  {
    NcmObjArray *data_i = ncm_obj_array_sized_new (2);
    NcmVector *vec      = ncm_vector_new (i + 1);

    ncm_obj_array_set (data_i, 0, G_OBJECT (vec));
    ncm_obj_array_set (data, j, G_OBJECT (data_i));
  }

  obs = g_object_new (NC_TYPE_GALAXY_WL_OBS,
                      "data", data,
                      "header", header,
                      "coord", coord,
                      NULL);

  return obs;
}

/**
 * nc_galaxy_wl_obs_set:
 * @obs: a #NcGalaxyWLObs object.
 * @col: column name.
 * @i: row index.
 * @val: value to be set.
 *
 * Sets a value in the observation data.
 *
 */
void
nc_galaxy_wl_obs_set (NcGalaxyWLObs *obs, const gchar *col, const guint i, gdouble val)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;
  NcmObjArray *data_i               = NCM_OBJ_ARRAY (ncm_obj_array_peek (self->data, i));
  NcmVarDict *header                = self->header;
  NcmVector *vec;
  gint j;

  if (g_strcmp0 (col, "pz") == 0)
    g_error ("nc_galaxy_wl_obs_set: call nc_galaxy_wl_obs_set_pz to set P(z) splines.");

  vec = NCM_VECTOR (ncm_obj_array_peek (data_i, 0));

  ncm_var_dict_get_int (header, col, &j);
  ncm_vector_set (vec, j, val);
}

/**
 * nc_galaxy_wl_obs_set_pz:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 * @pz: a #NcmSpline.
 *
 * Sets a P(z) spline in the observation data.
 *
 */
void
nc_galaxy_wl_obs_set_pz (NcGalaxyWLObs *obs, const guint i, NcmSpline *pz)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;
  NcmObjArray *data_i               = NCM_OBJ_ARRAY (ncm_obj_array_peek (self->data, i));

  ncm_obj_array_set (data_i, 1, G_OBJECT (pz));
}

/**
 * nc_galaxy_wl_obs_get:
 * @obs: a #NcGalaxyWLObs object.
 * @col: column name.
 * @i: row index.
 *
 * Gets a value from the observation data.
 *
 */
gdouble
nc_galaxy_wl_obs_get (NcGalaxyWLObs *obs, const gchar *col, const guint i)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;
  NcmObjArray *data_i               = NCM_OBJ_ARRAY (ncm_obj_array_peek (self->data, i));
  NcmVarDict *header                = self->header;
  NcmVector *vec;
  gint j;

  if (g_strcmp0 (col, "pz") == 0)
    g_error ("nc_galaxy_wl_obs_get: call nc_galaxy_wl_obs_get_pz to get P(z) splines.");

  ncm_var_dict_get_int (header, col, &j);

  vec = NCM_VECTOR (ncm_obj_array_peek (data_i, 0));

  return ncm_vector_get (vec, j);
}

/**
 * nc_galaxy_wl_obs_peek_pz:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 *
 * Gets a P(z) spline from the observation data.
 *
 * Returns: (transfer full): the P(z) spline.
 *
 */
NcmSpline *
nc_galaxy_wl_obs_peek_pz (NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;
  NcmObjArray *data_i               = NCM_OBJ_ARRAY (ncm_obj_array_peek (self->data, i));

  return NCM_SPLINE (ncm_obj_array_peek (data_i, 1));
}

/**
 * nc_galaxy_wl_obs_peek_header:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the data columns headers.
 *
 * Returns: (transfer none): the data columns headers.
 *
 */
NcmVarDict *
nc_galaxy_wl_obs_peek_header (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;

  return self->header;
}

/**
 * nc_galaxy_wl_obs_set_coord:
 * @obs: a #NcGalaxyWLObs object.
 * @coord: the coordinate system.
 *
 * Sets the coordinate system used to store the data.
 *
 */
void
nc_galaxy_wl_obs_set_coord (NcGalaxyWLObs *obs, NcGalaxyWLObsCoord coord)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;

  self->coord = coord;
}

/**
 * nc_galaxy_wl_obs_get_coord:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the coordinate system used to store the data.
 *
 * Returns: the coordinate system.
 *
 */
NcGalaxyWLObsCoord
nc_galaxy_wl_obs_get_coord (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;

  return self->coord;
}

/**
 * nc_galaxy_wl_obs_len:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the number of rows of the observation data.
 *
 */
gdouble
nc_galaxy_wl_obs_len (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = obs->priv;

  return ncm_obj_array_len (self->data);
}

/**
 * nc_galaxy_wl_obs_ref:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Increases the reference count of the #NcGalaxyWLObs object.
 *
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_ref (NcGalaxyWLObs *obs)
{
  return g_object_ref (obs);
}

/**
 * nc_galaxy_wl_obs_free:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Frees the #NcGalaxyWLObs object.
 *
 */
void
nc_galaxy_wl_obs_free (NcGalaxyWLObs *obs)
{
  g_object_unref (obs);
}

/**
 * nc_galaxy_wl_obs_clear:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Clears the #NcGalaxyWLObs object.
 *
 */
void
nc_galaxy_wl_obs_clear (NcGalaxyWLObs **obs)
{
  g_clear_object (obs);
}


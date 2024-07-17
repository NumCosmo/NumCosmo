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

#include "galaxy/nc_galaxy_wl_obs.h"
#include "math/ncm_matrix.h"
#include "nc_enum_types.h"

enum
{
  PROP_0,
  PROP_DATA,
  PROP_COORD,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyWLObs, nc_galaxy_wl_obs, G_TYPE_OBJECT)

static void
nc_galaxy_wl_obs_init (NcGalaxyWLObs *obs)
{
  obs->data  = NULL;
  obs->coord = NC_GALAXY_WL_OBS_COORD_PIXEL;
  obs->len   = 0;
}

static void
nc_galaxy_wl_obs_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs         = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate *priv = obs->priv;

  switch (prop_id)
  {
    case PROP_DATA:
      g_value_set_object (value, nc_galaxy_wl_obs_get_data (obs));
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
    case PROP_DATA:
      nc_galaxy_wl_obs_set_data (obs, g_value_get_object (value));
      break;
    case PROP_COORD:
      nc_galaxy_wl_obs_set_coord (obs, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
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
                                   g_param_spec_object ("data",
                                                        "Data",
                                                        "Weak lensing observation data",
                                                        NCM_TYPE_MATRIX,
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
                                                      NC_GALAXY_WL_OBS_COORD_PIXEL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));
}

/**
 * nc_galaxy_wl_obs_new:
 * @nrows: number of rows of the observation data.
 *
 * Creates a new #NcGalaxyWLObs object.
 *
 * Returns: a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new (const guint nrows, NcGalaxyWLObsCoord coord)
{
  NcGalaxyWLObs *obs = g_object_new (NC_TYPE_GALAXY_WL_OBS,
                                     "data", ncm_matrix_new (nrows, 5),
                                     "coord", coord,
                                     NULL);

  return obs;
}

/**
 * nc_galaxy_wl_obs_set:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 * @j: column index.
 * @val: value to be set.
 *
 * Sets a value in the observation data.
 *
 */
void
nc_galaxy_wl_obs_set (NcGalaxyWLObs *obs, const guint i, const guint j, gdouble val)
{
  ncm_matrix_set (obs->data, i, j, val);
}

/**
 * nc_galaxy_wl_obs_get:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 * @j: column index.
 *
 * Gets a value from the observation data.
 *
 */
gdouble
nc_galaxy_wl_obs_get (NcGalaxyWLObs *obs, const guint i, const guint j)
{
  return ncm_matrix_get (obs->data, i, j);
}

/**
 * nc_galaxy_wl_obs_get_data:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the observation data.
 *
 * Returns: the observation data.
 *
 */
NcmMatrix *
nc_galaxy_wl_obs_get_data (NcGalaxyWLObs *obs)
{
  return obs->data;
}

/**
 * nc_galaxy_wl_obs_set_data:
 * @obs: a #NcGalaxyWLObs object.
 * @data: the observation data.
 *
 * Sets the observation data.
 *
 */
void
nc_galaxy_wl_obs_set_data (NcGalaxyWLObs *obs, NcmMatrix *data)
{
  g_assert_cmpuint (ncm_matrix_ncols (data), ==, 5);

  if (obs->data)
    g_object_unref (obs->data);

  obs->data = g_object_ref (data);
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
  return obs->coord;
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
  obs->coord = coord;
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
  return ncm_matrix_nrows (obs->data);
}

/**
 * nc_galaxy_wl_obs_dispose:
 * @object: a #GObject.
 *
 * Disposes the #NcGalaxyWLObs object.
 *
 */
static void
nc_galaxy_wl_obs_dispose (GObject *object)
{
  NcGalaxyWLObs *obs = NC_GALAXY_WL_OBS (object);

  if (obs->data)
    g_object_unref (obs->data);

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
  NcGalaxyWLObs *obs = NC_GALAXY_WL_OBS (object);

  if (obs->data)
    g_object_unref (obs->data);

  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->finalize (object);
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


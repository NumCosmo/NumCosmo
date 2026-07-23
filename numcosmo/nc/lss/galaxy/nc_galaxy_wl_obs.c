/***************************************************************************
 *            nc_galaxy_wl_obs.c
 *
 *  Tue Jul 16 06:43:45 2024
 *  Copyright  2024 Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyWLObs:
 *
 * Galaxy weak lensing observation data.
 *
 * A class to store galaxy weak lensing observation data and information about its
 * coordinate system.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "stdarg.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/core/ncm_obj_array.h"
#include "ncm/spline/ncm_spline.h"
#include "ncm/algebra/ncm_vector.h"
#include "nc_enum_types.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"


struct _NcGalaxyWLObs
{
  /*< private >*/
  NcmCatalog parent_instance;
};

struct _NcGalaxyWLObsPrivate
{
  NcmObjDictInt *pz;
  NcWLEllipticityFrame coord;
  NcGalaxyWLObsEllipConv ellip_conv;
};

enum
{
  PROP_PZ = 1000,
  PROP_COORD,
  PROP_ELLIP_CONV,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLObs, nc_galaxy_wl_obs, NCM_TYPE_CATALOG)

static void
nc_galaxy_wl_obs_init (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->pz         = ncm_obj_dict_int_new ();
  self->coord      = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;
  self->ellip_conv = 0;
}

static void
_nc_galaxy_wl_obs_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs                = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  switch (prop_id)
  {
    case PROP_PZ:
      g_value_set_boxed (value, self->pz);
      break;
    case PROP_COORD:
      g_value_set_enum (value, nc_galaxy_wl_obs_get_coord (obs));
      break;
    case PROP_ELLIP_CONV:
      g_value_set_enum (value, nc_galaxy_wl_obs_get_ellip_conv (obs));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _nc_galaxy_wl_obs_take_pz (NcGalaxyWLObsPrivate *self, NcmObjDictInt *pz);
static void _nc_galaxy_wl_obs_set_ellip_conv (NcGalaxyWLObs *obs, NcGalaxyWLObsEllipConv type);

static void
_nc_galaxy_wl_obs_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs         = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate *self = nc_galaxy_wl_obs_get_instance_private (obs);

  switch (prop_id)
  {
    case PROP_PZ:
      _nc_galaxy_wl_obs_take_pz (self, g_value_dup_boxed (value));
      break;
    case PROP_COORD:
      nc_galaxy_wl_obs_set_coord (obs, g_value_get_enum (value));
      break;
    case PROP_ELLIP_CONV:
      _nc_galaxy_wl_obs_set_ellip_conv (obs, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_obs_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->constructed (object);
}

static void
_nc_galaxy_wl_obs_dispose (GObject *object)
{
  NcGalaxyWLObs *obs                = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  ncm_obj_dict_int_clear (&self->pz);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->dispose (object);
}

/**
 * _nc_galaxy_wl_obs_finalize:
 * @object: a #GObject.
 *
 * Finalizes the #NcGalaxyWLObs object.
 *
 */
static void
_nc_galaxy_wl_obs_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_obs_class_init (NcGalaxyWLObsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->get_property = &_nc_galaxy_wl_obs_get_property;
  object_class->set_property = &_nc_galaxy_wl_obs_set_property;
  object_class->constructed  = &_nc_galaxy_wl_obs_constructed;
  object_class->dispose      = &_nc_galaxy_wl_obs_dispose;
  object_class->finalize     = &_nc_galaxy_wl_obs_finalize;

  /**
   * NcGalaxyWLObs:pz:
   *
   * P(z) splines.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ,
                                   g_param_spec_boxed ("pz",
                                                       "P(z)",
                                                       "P(z) splines",
                                                       NCM_TYPE_OBJ_DICT_INT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

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
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

/**
 * NcGalaxyWLObs:ellip-conv:
 *
 * Ellipticity convention.
 *
 */
  g_object_class_install_property (object_class,
                                   PROP_ELLIP_CONV,
                                   g_param_spec_enum ("ellip-conv",
                                                      "Ellipticity convention",
                                                      "Weak lensing observables ellipticity convention",
                                                      NC_TYPE_GALAXY_WL_OBS_ELLIP_CONV,
                                                      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));
}

void
_nc_galaxy_wl_obs_take_pz (NcGalaxyWLObsPrivate *self, NcmObjDictInt *pz)
{
  guint i;

  ncm_obj_dict_int_clear (&self->pz);
  self->pz = pz;

  if (self->pz != NULL)
  {
    GArray *keys  = ncm_obj_dict_int_keys (self->pz);
    const guint n = keys->len;

    for (i = 0; i < n; i++)
    {
      const gint key = g_array_index (keys, gint, i);
      GObject *obj   = ncm_obj_dict_int_peek (self->pz, key);

      g_assert (NCM_IS_SPLINE (obj));
      ncm_spline_prepare (NCM_SPLINE (obj));
    }

    g_array_unref (keys);
  }
  else
  {
    self->pz = ncm_obj_dict_int_new ();
  }
}

/**
 * nc_galaxy_wl_obs_new:
 * @ellip_conv: the ellipticity convention
 * @coord: the coordinate system used to store the data
 * @nrows: the number of data rows
 * @col_names: a #GStrv
 *
 * Creates a new #NcGalaxyWLObs object.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new (NcGalaxyWLObsEllipConv ellip_conv, NcWLEllipticityFrame coord, guint nrows, GStrv col_names)
{
  return g_object_new (NC_TYPE_GALAXY_WL_OBS,
                       "ellip-conv", ellip_conv,
                       "columns", col_names,
                       "coord", coord,
                       "len", nrows,
                       NULL);
}

/* Filenames of the curated Subaru HSC-SSP PDR1 catalogs in the NumCosmo
 * datafile-release-v1.0.0 GitHub release, indexed by #NcGalaxyWLObsCatalogId. */
static const gchar *_nc_galaxy_wl_obs_catalog_files[] = {
  "wl_obs_HWL16a-002.gvar",
  "wl_obs_HWL16a-007.gvar",
  "wl_obs_HWL16a-060.gvar",
  "wl_obs_HWL16a-064.gvar",
  "wl_obs_HWL16a-094.gvar",
};

/**
 * nc_galaxy_wl_obs_catalog_id_get_filename:
 * @id: a #NcGalaxyWLObsCatalogId
 *
 * Downloads (if not already cached) and returns the local filename of the
 * catalog identified by @id. The file is fetched from the NumCosmo
 * datafile-release-v1.0.0 GitHub release into the NumCosmo data directory
 * (see ncm_cfg_get_fullpath_base()) the first time it is requested; later
 * calls reuse the cached copy.
 *
 * Returns: (transfer full): Full path for the catalog file.
 */
gchar *
nc_galaxy_wl_obs_catalog_id_get_filename (NcGalaxyWLObsCatalogId id)
{
  g_return_val_if_fail (id < NC_GALAXY_WL_OBS_CATALOG_LEN, NULL);

  {
    const gchar *filename = _nc_galaxy_wl_obs_catalog_files[id];
    gchar *full_filename  = ncm_cfg_get_fullpath (filename);

    if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
    {
      gchar *url_str   = g_strdup_printf ("https://github.com/NumCosmo/NumCosmo/releases/download/datafile-release-v1.0.0/%s", filename);
      const gchar *dir = ncm_cfg_get_fullpath_base ();
      gchar *cmd[]     = {"wget", "--tries=3", "--timeout=30", "-O", full_filename, url_str, NULL };
      GError *error    = NULL;

      ncm_message ("# Downloading file [%s]...\n", url_str);

      if (!g_spawn_sync (dir, cmd, NULL,
                         G_SPAWN_SEARCH_PATH | G_SPAWN_STDOUT_TO_DEV_NULL | G_SPAWN_STDERR_TO_DEV_NULL,
                         NULL, NULL, NULL, NULL, NULL, &error))
        g_error ("nc_galaxy_wl_obs_catalog_id_get_filename: cannot download file: %s. Error: %s. "
                 "Please download the file manually from %s and extract it to %s.",
                 filename, error->message,
                 url_str, dir);

      g_free (url_str);
    }

    return full_filename;
  }
}

/**
 * nc_galaxy_wl_obs_new_from_catalog_id: (constructor)
 * @id: a #NcGalaxyWLObsCatalogId
 *
 * Downloads (if necessary) and loads the catalog identified by @id.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new_from_catalog_id (NcGalaxyWLObsCatalogId id)
{
  gchar *full_filename = nc_galaxy_wl_obs_catalog_id_get_filename (id);
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GObject *obj         = ncm_serialize_from_binfile (ser, full_filename);

  g_free (full_filename);
  ncm_serialize_free (ser);

  if (!NC_IS_GALAXY_WL_OBS (obj))
    g_error ("nc_galaxy_wl_obs_new_from_catalog_id: file for catalog id %d does not contain a NcGalaxyWLObs.", id);

  return NC_GALAXY_WL_OBS (obj);
}

/**
 * nc_galaxy_wl_obs_ref:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Increases the reference count of the #NcGalaxyWLObs object.
 *
 * Returns: (transfer full): the #NcGalaxyWLObs object.
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

/**
 * nc_galaxy_wl_obs_get_index:
 * @obs: a #NcGalaxyWLObs object
 * @col: column name
 * @i: a pointer to store the column index
 *
 * Gets the column index from the column name.
 *
 * Returns: %TRUE if the column was found, %FALSE otherwise.
 */
gboolean
nc_galaxy_wl_obs_get_index (NcGalaxyWLObs *obs, const gchar *col, guint *i)
{
  return ncm_catalog_get_index (NCM_CATALOG (obs), col, i);
}

/**
 * nc_galaxy_wl_obs_set:
 * @obs: a #NcGalaxyWLObs object
 * @col: column name
 * @i: row index
 * @val: value to be set
 * @error: a #GError for error reporting
 *
 * Sets a value in the observation data.
 *
 */
void
nc_galaxy_wl_obs_set (NcGalaxyWLObs *obs, const gchar *col, const guint i, gdouble val, GError **error)
{
  ncm_catalog_set (NCM_CATALOG (obs), col, i, val, error);
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
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  ncm_obj_dict_int_add (self->pz, i, G_OBJECT (pz));
}

/**
 * nc_galaxy_wl_obs_get:
 * @obs: a #NcGalaxyWLObs object.
 * @col: column name.
 * @i: row index.
 * @error: a #GError for error reporting
 *
 * Gets a value from the observation data.
 *
 * Returns: the requested value.
 */
gdouble
nc_galaxy_wl_obs_get (NcGalaxyWLObs *obs, const gchar *col, const guint i, GError **error)
{
  return ncm_catalog_get (NCM_CATALOG (obs), col, i, error);
}

/**
 * nc_galaxy_wl_obs_peek_data:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the observation data matrix.
 *
 * Returns: (transfer none): the observation data matrix.
 */
NcmMatrix *
nc_galaxy_wl_obs_peek_data (NcGalaxyWLObs *obs)
{
  return ncm_catalog_peek_data (NCM_CATALOG (obs));
}

/**
 * nc_galaxy_wl_obs_peek_pz:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 *
 * Gets a P(z) spline from the observation data.
 *
 * Returns: (transfer none): the P(z) spline.
 *
 */
NcmSpline *
nc_galaxy_wl_obs_peek_pz (NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return NCM_SPLINE (ncm_obj_dict_int_peek (self->pz, i));
}

/**
 * nc_galaxy_wl_obs_peek_columns:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the data columns names.
 *
 * Returns: (transfer none): the data columns names.
 */
GStrv
nc_galaxy_wl_obs_peek_columns (NcGalaxyWLObs *obs)
{
  return ncm_catalog_peek_columns (NCM_CATALOG (obs));
}

/**
 * nc_galaxy_wl_obs_set_coord:
 * @obs: a #NcGalaxyWLObs object.
 * @coord: the ellipticity handedness frame #NcWLEllipticityFrame.
 *
 * Sets the frame in which the stored ellipticity components are expressed. Sky
 * positions are always RA/Dec; this only selects the ellipticity basis (see
 * #NcWLEllipticityFrame).
 *
 */
void
nc_galaxy_wl_obs_set_coord (NcGalaxyWLObs *obs, NcWLEllipticityFrame coord)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->coord = coord;
}

static void
_nc_galaxy_wl_obs_set_ellip_conv (NcGalaxyWLObs *obs, NcGalaxyWLObsEllipConv type)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->ellip_conv = type;
}

/**
 * nc_galaxy_wl_obs_get_ellip_conv:
 * @obs: a #NcGalaxyWLObs object
 *
 * Gets the ellipse type used to store the data.
 *
 * Returns: the ellipse type #NcGalaxyWLObsEllipConv
 */
NcGalaxyWLObsEllipConv
nc_galaxy_wl_obs_get_ellip_conv (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return self->ellip_conv;
}

/**
 * nc_galaxy_wl_obs_get_coord:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the frame in which the stored ellipticity components are expressed (see
 * #NcWLEllipticityFrame).
 *
 * Returns: the ellipticity handedness frame #NcWLEllipticityFrame.
 *
 */
NcWLEllipticityFrame
nc_galaxy_wl_obs_get_coord (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return self->coord;
}

/**
 * nc_galaxy_wl_obs_len:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the number of rows of the observation data.
 *
 */
guint
nc_galaxy_wl_obs_len (NcGalaxyWLObs *obs)
{
  return ncm_catalog_len (NCM_CATALOG (obs));
}


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
 * A catalog may also carry the interpolants of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of its
 * galaxies (#NcGalaxyWLObs:shape-interp) with the record of what built them
 * (#NcGalaxyWLObs:shape-interp-source). They describe the rows as they were when the record
 * was set, by the key of #NcmCatalog (ncm_catalog_key_is_current()): a later change of the
 * rows makes nc_galaxy_wl_obs_peek_shape_interp() and
 * nc_galaxy_wl_obs_peek_shape_interp_source() return %NULL and the two properties return
 * empty dictionaries, so a serialized catalog never carries interpolants of other rows.
 *
 * The interpolants and their record may instead be those of a catalog file
 * (#NcGalaxyWLObs:shape-interp-file), as for the catalog returned by
 * nc_galaxy_wl_obs_new_from_file(). Such a catalog serializes the name of the file in place
 * of the interpolants and their record, and a deserialized one reads them from the file at
 * the first call of nc_galaxy_wl_obs_peek_shape_interp() or
 * nc_galaxy_wl_obs_peek_shape_interp_source() or the first read of one of the three
 * properties. They are taken only when the columns and the rows of the file equal those of
 * the catalog; otherwise the catalog has no interpolants. Setting an interpolant or a record
 * unsets the file.
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
#include "ncm/algebra/ncm_spectral_tensor.h"
#include "nc_enum_types.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/data/nc_data_download_priv.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"

/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-galaxy-wl-obs-error, nc_galaxy_wl_obs_error) 
/* *INDENT-ON* */

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
  NcmObjDictInt *shape_interp;
  NcmObjDictStr *shape_interp_source;
  guint64 shape_interp_key;           /* catalog key of the rows the interpolants describe */
  gchar *shape_interp_file;           /* catalog file that holds the interpolants and their record */
  gboolean shape_interp_file_pending; /* the file was set by the property and is not read yet */
};

enum
{
  PROP_PZ = 1000,
  PROP_COORD,
  PROP_ELLIP_CONV,
  PROP_SHAPE_INTERP,
  PROP_SHAPE_INTERP_SOURCE,
  PROP_SHAPE_INTERP_FILE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLObs, nc_galaxy_wl_obs, NCM_TYPE_CATALOG)

static void
nc_galaxy_wl_obs_init (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->pz                  = ncm_obj_dict_int_new ();
  self->coord               = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;
  self->ellip_conv          = 0;
  self->shape_interp        = ncm_obj_dict_int_new ();
  self->shape_interp_source = NULL;
  self->shape_interp_key    = NCM_CATALOG_KEY_PENDING;
  self->shape_interp_file   = NULL;

  self->shape_interp_file_pending = FALSE;
}

static gboolean _nc_galaxy_wl_obs_shape_interp_current (NcGalaxyWLObs *obs);
static void _nc_galaxy_wl_obs_unset_shape_interp_file (NcGalaxyWLObsPrivate *self);

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
    case PROP_SHAPE_INTERP:

      if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
        g_value_take_boxed (value, ncm_obj_dict_int_new ());
      else if (self->shape_interp_file == NULL)
        g_value_set_boxed (value, self->shape_interp);

      break;
    case PROP_SHAPE_INTERP_SOURCE:

      if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
        g_value_take_boxed (value, ncm_obj_dict_str_new ());
      else if (self->shape_interp_file == NULL)
        g_value_set_boxed (value, self->shape_interp_source);

      break;
    case PROP_SHAPE_INTERP_FILE:
      g_value_set_string (value, nc_galaxy_wl_obs_peek_shape_interp_file (obs));
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
    case PROP_SHAPE_INTERP:
    {
      NcmObjDictInt *shape_interp = g_value_dup_boxed (value);

      _nc_galaxy_wl_obs_unset_shape_interp_file (self);
      ncm_obj_dict_int_clear (&self->shape_interp);
      self->shape_interp = (shape_interp != NULL) ? shape_interp : ncm_obj_dict_int_new ();
      break;
    }
    case PROP_SHAPE_INTERP_SOURCE:
    {
      NcmObjDictStr *source = g_value_dup_boxed (value);

      _nc_galaxy_wl_obs_unset_shape_interp_file (self);
      ncm_obj_dict_str_clear (&self->shape_interp_source);

      /* An empty record is no record: it is what a stale catalog serializes */
      if ((source != NULL) && (ncm_obj_dict_str_len (source) > 0))
        self->shape_interp_source = source;
      else if (source != NULL)
        ncm_obj_dict_str_unref (source);

      self->shape_interp_key = NCM_CATALOG_KEY_PENDING;
      break;
    }
    case PROP_SHAPE_INTERP_FILE:
    {
      const gchar *filename = g_value_get_string (value);

      _nc_galaxy_wl_obs_unset_shape_interp_file (self);

      /* The file is read at the first check: the rows it is compared with may not be set yet */
      if ((filename != NULL) && (filename[0] != '\0'))
      {
        self->shape_interp_file         = g_strdup (filename);
        self->shape_interp_file_pending = TRUE;
      }

      break;
    }
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
  ncm_obj_dict_int_clear (&self->shape_interp);
  ncm_obj_dict_str_clear (&self->shape_interp_source);
  g_clear_pointer (&self->shape_interp_file, g_free);

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
   * The redshift densities of the galaxies: a #NcGalaxyPz for each row that has one,
   * keyed by the row.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ,
                                   g_param_spec_boxed ("pz",
                                                       "P(z)",
                                                       "Redshift density of each galaxy",
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

  /**
   * NcGalaxyWLObs:shape-interp:
   *
   * The interpolants of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of the galaxies, a
   * #NcmSpectralTensor for each row that has one, keyed by the row; empty when the
   * rows changed after #NcGalaxyWLObs:shape-interp-source was set, and %NULL when the
   * interpolants are those of #NcGalaxyWLObs:shape-interp-file.
   */
  g_object_class_install_property (object_class,
                                   PROP_SHAPE_INTERP,
                                   g_param_spec_boxed ("shape-interp",
                                                       "Shape interpolants",
                                                       "Interpolants of ln P in g of the galaxies",
                                                       NCM_TYPE_OBJ_DICT_INT,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObs:shape-interp-source:
   *
   * The record of what built #NcGalaxyWLObs:shape-interp, objects keyed by their role (see
   * nc_galaxy_shape_factor_spectral_set_obs_source()); empty when there is none or when the
   * rows changed after it was set, and %NULL when the record is that of
   * #NcGalaxyWLObs:shape-interp-file.
   */
  g_object_class_install_property (object_class,
                                   PROP_SHAPE_INTERP_SOURCE,
                                   g_param_spec_boxed ("shape-interp-source",
                                                       "Shape interpolants source",
                                                       "Record of what built the interpolants",
                                                       NCM_TYPE_OBJ_DICT_STR,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObs:shape-interp-file:
   *
   * The catalog file that holds the interpolants and their record, a serialized
   * #NcGalaxyWLObs with the columns and the rows of this one; %NULL when there is none or
   * when the rows changed after it was read. While it is set, the catalog serializes it in
   * place of #NcGalaxyWLObs:shape-interp and #NcGalaxyWLObs:shape-interp-source. Setting it
   * replaces the interpolants and the record by those of the file, read at the first use.
   * A file that cannot be read, that holds another type of object or that has no
   * interpolants aborts with an error.
   */
  g_object_class_install_property (object_class,
                                   PROP_SHAPE_INTERP_FILE,
                                   g_param_spec_string ("shape-interp-file",
                                                        "Shape interpolants file",
                                                        "Catalog file that holds the interpolants and their record",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

static void _nc_galaxy_wl_obs_read_shape_interp_file (NcGalaxyWLObs *obs);

/*
 * Whether the interpolants and their record describe the current rows. A file set by the
 * property is read here first.
 */
static gboolean
_nc_galaxy_wl_obs_shape_interp_current (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if (self->shape_interp_file_pending)
    _nc_galaxy_wl_obs_read_shape_interp_file (obs);

  return (self->shape_interp_source != NULL) && ncm_catalog_key_is_current (NCM_CATALOG (obs), &self->shape_interp_key);
}

static gboolean _nc_galaxy_wl_obs_same_rows (NcGalaxyWLObs *obs, NcGalaxyWLObs *other);

/*
 * Reads the catalog of shape_interp_file and takes its interpolants and its record, shared
 * with no other object, as those of the current rows of obs. They are taken only when the
 * columns and the rows of the file equal those of obs; otherwise obs is left without
 * interpolants, record and file.
 */
static void
_nc_galaxy_wl_obs_read_shape_interp_file (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);
  NcmSerialize *ser                 = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  GObject *file_obj                 = ncm_serialize_from_binfile (ser, self->shape_interp_file);

  self->shape_interp_file_pending = FALSE;

  if (!NC_IS_GALAXY_WL_OBS (file_obj))
    g_error ("nc_galaxy_wl_obs: the file %s of shape-interp-file holds a %s, not a NcGalaxyWLObs.",
             self->shape_interp_file, G_OBJECT_TYPE_NAME (file_obj));

  {
    NcGalaxyWLObs *file_obs                = NC_GALAXY_WL_OBS (file_obj);
    NcGalaxyWLObsPrivate * const file_self = nc_galaxy_wl_obs_get_instance_private (file_obs);

    if (!_nc_galaxy_wl_obs_shape_interp_current (file_obs))
      g_error ("nc_galaxy_wl_obs: the catalog of the file %s of shape-interp-file has no interpolants.",
               self->shape_interp_file);

    ncm_obj_dict_int_clear (&self->shape_interp);
    ncm_obj_dict_str_clear (&self->shape_interp_source);

    if (_nc_galaxy_wl_obs_same_rows (obs, file_obs))
    {
      self->shape_interp        = ncm_obj_dict_int_ref (file_self->shape_interp);
      self->shape_interp_source = ncm_obj_dict_str_ref (file_self->shape_interp_source);
      self->shape_interp_key    = ncm_catalog_get_key (NCM_CATALOG (obs));
    }
    else
    {
      self->shape_interp = ncm_obj_dict_int_new ();
      g_clear_pointer (&self->shape_interp_file, g_free);
    }

    nc_galaxy_wl_obs_clear (&file_obs);
  }

  ncm_serialize_free (ser);
}

/* Whether obs and other have the same column names and the same values in every row */
static gboolean
_nc_galaxy_wl_obs_same_rows (NcGalaxyWLObs *obs, NcGalaxyWLObs *other)
{
  NcmMatrix *data       = ncm_catalog_peek_data (NCM_CATALOG (obs));
  NcmMatrix *other_data = ncm_catalog_peek_data (NCM_CATALOG (other));
  const guint nrows     = ncm_matrix_nrows (data);
  const guint ncols     = ncm_matrix_ncols (data);
  guint i;

  if ((nrows != ncm_matrix_nrows (other_data)) || (ncols != ncm_matrix_ncols (other_data)))
    return FALSE;

  if (!g_strv_equal ((const gchar * const *) ncm_catalog_peek_columns (NCM_CATALOG (obs)),
                     (const gchar * const *) ncm_catalog_peek_columns (NCM_CATALOG (other))))
    return FALSE;

  for (i = 0; i < nrows; i++)
  {
    if (memcmp (ncm_matrix_ptr (data, i, 0), ncm_matrix_ptr (other_data, i, 0), sizeof (gdouble) * ncols) != 0)
      return FALSE;
  }

  return TRUE;
}

static void
_nc_galaxy_wl_obs_unset_shape_interp_file (NcGalaxyWLObsPrivate *self)
{
  g_clear_pointer (&self->shape_interp_file, g_free);
  self->shape_interp_file_pending = FALSE;
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

      /* Catalogs written by NumCosmo 0.x before NcGalaxyPz: this check is dropped in 1.0 */
      if (NCM_IS_SPLINE (obj))
        g_error ("nc_galaxy_wl_obs: the catalog holds the p(z) of row %d as a %s, the format written before NcGalaxyPz. "
                 "Update the file with tools/nc_galaxy_pz_update_catalog.py. A catalog downloaded by NumCosmo "
                 "(wl_obs_*.gvar in the NumCosmo data directory) can instead be removed, so that it is downloaded again.",
                 key, G_OBJECT_TYPE_NAME (obj));

      g_assert (NC_IS_GALAXY_PZ (obj));
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
 * datafile-release-v1.0.0 GitHub release into the NumCosmo user data directory
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
      gchar *url_str = g_strdup_printf ("https://github.com/NumCosmo/NumCosmo/releases/download/datafile-release-v1.0.0/%s", filename);
      gchar *lockdir = NULL;

      /* Serialize with any other process fetching the same catalog, and let it
       * publish by rename: writing straight to full_filename let a concurrent
       * reader open a half-written catalog, and an interrupted fetch left a
       * truncated one that every later run took for complete. */
      if (_nc_data_download_lock (full_filename, full_filename, 900, &lockdir))
      {
        _nc_data_download_file (url_str, full_filename, filename, lockdir);
        _nc_data_download_unlock (lockdir);
      }

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
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  GObject *obj         = ncm_serialize_from_binfile (ser, full_filename);

  g_free (full_filename);
  ncm_serialize_free (ser);

  if (!NC_IS_GALAXY_WL_OBS (obj))
    g_error ("nc_galaxy_wl_obs_new_from_catalog_id: file for catalog id %d does not contain a NcGalaxyWLObs.", id);

  return NC_GALAXY_WL_OBS (obj);
}

/**
 * nc_galaxy_wl_obs_new_from_file: (constructor)
 * @filename: the name of a file with a serialized #NcGalaxyWLObs
 * @error: a #GError
 *
 * Reads the catalog of the binary file @filename (ncm_serialize_from_binfile()). When the
 * file holds interpolants, the catalog returned has the absolute name of the file as its
 * #NcGalaxyWLObs:shape-interp-file: serializing it writes that name in place of the
 * interpolants and their record. A file that holds another type of object sets
 * %NC_GALAXY_WL_OBS_ERROR_WRONG_TYPE; a file that cannot be read aborts with an error.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new_from_file (const gchar *filename, GError **error)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  GObject *obj      = ncm_serialize_from_binfile (ser, filename);

  ncm_serialize_free (ser);

  if (!NC_IS_GALAXY_WL_OBS (obj))
  {
    ncm_util_set_or_call_error (error, NC_GALAXY_WL_OBS_ERROR, NC_GALAXY_WL_OBS_ERROR_WRONG_TYPE,
                                "nc_galaxy_wl_obs_new_from_file: the file %s holds a %s, not a NcGalaxyWLObs.",
                                filename, G_OBJECT_TYPE_NAME (obj));
    g_clear_object (&obj);

    return NULL;
  }

  {
    NcGalaxyWLObs *obs                = NC_GALAXY_WL_OBS (obj);
    NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

    /* A catalog that refers to another file keeps that file: @filename has no interpolants */
    if (_nc_galaxy_wl_obs_shape_interp_current (obs) && (self->shape_interp_file == NULL))
      self->shape_interp_file = g_canonicalize_filename (filename, NULL);

    return obs;
  }
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
 * @pz: a #NcGalaxyPz.
 *
 * Sets the redshift density of row @i to @pz, which is referenced. A #NcGalaxyPz does
 * not change after construction, so one object may be set in several rows and catalogs.
 *
 */
void
nc_galaxy_wl_obs_set_pz (NcGalaxyWLObs *obs, const guint i, NcGalaxyPz *pz)
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
 * Gets the redshift density of row @i.
 *
 * Returns: (transfer none) (nullable): the #NcGalaxyPz of row @i, or %NULL when the row
 * has none.
 *
 */
NcGalaxyPz *
nc_galaxy_wl_obs_peek_pz (NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return (NcGalaxyPz *) ncm_obj_dict_int_peek (self->pz, i);
}

/**
 * nc_galaxy_wl_obs_set_shape_interp:
 * @obs: a #NcGalaxyWLObs
 * @i: row index
 * @st: the interpolant of row @i
 *
 * Sets the interpolant of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of row @i, referencing
 * @st. It is returned by nc_galaxy_wl_obs_peek_shape_interp() only after
 * nc_galaxy_wl_obs_set_shape_interp_source() and while the rows are unchanged. The call
 * unsets #NcGalaxyWLObs:shape-interp-file.
 */
void
nc_galaxy_wl_obs_set_shape_interp (NcGalaxyWLObs *obs, const guint i, NcmSpectralTensor *st)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  _nc_galaxy_wl_obs_unset_shape_interp_file (self);
  ncm_obj_dict_int_set (self->shape_interp, i, G_OBJECT (st));
}

/**
 * nc_galaxy_wl_obs_peek_shape_interp:
 * @obs: a #NcGalaxyWLObs
 * @i: row index
 *
 * An object of another type in #NcGalaxyWLObs:shape-interp aborts with an error.
 *
 * Returns: (transfer none) (nullable): the interpolant of row @i, or %NULL when the row has
 * none, there is no record or the rows changed after the record was set.
 */
NcmSpectralTensor *
nc_galaxy_wl_obs_peek_shape_interp (NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);
  GObject *interp;

  if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
    return NULL;

  interp = ncm_obj_dict_int_peek (self->shape_interp, i);

  if ((interp != NULL) && !NCM_IS_SPECTRAL_TENSOR (interp))
    g_error ("nc_galaxy_wl_obs_peek_shape_interp: the interpolant of row %u is a %s, not a NcmSpectralTensor.",
             i, G_OBJECT_TYPE_NAME (interp));

  return (NcmSpectralTensor *) interp;
}

/**
 * nc_galaxy_wl_obs_set_shape_interp_source:
 * @obs: a #NcGalaxyWLObs
 * @source: the record of what built the interpolants
 *
 * Sets the record of the interpolants set with nc_galaxy_wl_obs_set_shape_interp(),
 * referencing @source, and records ncm_catalog_get_key(): the interpolants describe the rows
 * as they are now. Set it after the last change of the rows. The call unsets
 * #NcGalaxyWLObs:shape-interp-file.
 */
void
nc_galaxy_wl_obs_set_shape_interp_source (NcGalaxyWLObs *obs, NcmObjDictStr *source)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  ncm_obj_dict_str_ref (source);
  _nc_galaxy_wl_obs_unset_shape_interp_file (self);
  ncm_obj_dict_str_clear (&self->shape_interp_source);
  self->shape_interp_source = source;
  self->shape_interp_key    = ncm_catalog_get_key (NCM_CATALOG (obs));
}

/**
 * nc_galaxy_wl_obs_peek_shape_interp_source:
 * @obs: a #NcGalaxyWLObs
 *
 * Returns: (transfer none) (nullable): the record of the interpolants, or %NULL when there is
 * none or the rows changed after it was set.
 */
NcmObjDictStr *
nc_galaxy_wl_obs_peek_shape_interp_source (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
    return NULL;

  return self->shape_interp_source;
}

/**
 * nc_galaxy_wl_obs_peek_shape_interp_file:
 * @obs: a #NcGalaxyWLObs
 *
 * Returns: (transfer none) (nullable): the name of the catalog file that holds the
 * interpolants and their record (#NcGalaxyWLObs:shape-interp-file), or %NULL when there is
 * none or the rows changed after it was read.
 */
const gchar *
nc_galaxy_wl_obs_peek_shape_interp_file (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
    return NULL;

  return self->shape_interp_file;
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


/***************************************************************************
 *            nc_galaxy_wl_obs_shape_interp.c
 *
 *  Wed November 25 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_shape_interp.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxyWLObsShapeInterp:
 *
 * The interpolants of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of the galaxies of a
 * #NcGalaxyWLObs, with the record of what built them.
 *
 * The object holds a #NcmSpectralTensor for each row that has an interpolant, keyed by
 * the row, the record (objects keyed by their role, see
 * nc_galaxy_shape_factor_spectral_set_obs_source()) and the data hash of the rows the
 * interpolants describe (ncm_catalog_get_data_hash()).
 *
 * The interpolants are either in memory or in a file, the binary serialization of an
 * object of this type that holds them in memory
 * (nc_galaxy_wl_obs_shape_interp_new_from_file()). Serializing an object in memory writes
 * the interpolants; serializing one backed by a file writes the name of the file, the
 * record and the data hash. The file is read at the first call of
 * nc_galaxy_wl_obs_shape_interp_peek_interp(), and its data hash must equal the one of the
 * object. A change of the interpolants, of the record or of the data hash of an object
 * backed by a file first reads the file; the object then holds the interpolants in memory.
 *
 * The name of the file is absolute in memory. A relative name given at construction is
 * resolved against the directory of ncm_serialize_context_peek_current(), which is the
 * directory of the file being read during deserialization. Written into a file, the name is
 * the bare file name when that file is in the same directory; writing it into a file of
 * another directory aborts with an error. Serialized without a file, as by
 * ncm_serialize_dup_obj(), the name is absolute.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_obs_shape_interp.h"
#include "ncm/core/ncm_serialize.h"

struct _NcGalaxyWLObsShapeInterp
{
  /*< private >*/
  GObject parent_instance;
  NcmObjDictInt *interp; /* row -> NcmSpectralTensor, NULL until a file is read */
  NcmObjDictStr *source;
  gchar *data_hash;
  gchar *filename; /* absolute, NULL in memory */
  gint loaded;     /* whether the file was read, read atomically */
  GMutex lock;     /* serializes the read of the file */
};

enum
{
  PROP_0,
  PROP_INTERP,
  PROP_SOURCE,
  PROP_DATA_HASH,
  PROP_FILENAME,
};

G_DEFINE_TYPE (NcGalaxyWLObsShapeInterp, nc_galaxy_wl_obs_shape_interp, G_TYPE_OBJECT)

static void
nc_galaxy_wl_obs_shape_interp_init (NcGalaxyWLObsShapeInterp *si)
{
  si->interp    = NULL;
  si->source    = NULL;
  si->data_hash = NULL;
  si->filename  = NULL;
  si->loaded    = FALSE;
  g_mutex_init (&si->lock);
}

static void
_nc_galaxy_wl_obs_shape_interp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObsShapeInterp *si = NC_GALAXY_WL_OBS_SHAPE_INTERP (object);

  switch (prop_id)
  {
    case PROP_INTERP:
      ncm_obj_dict_int_clear (&si->interp);
      si->interp = g_value_dup_boxed (value);
      break;
    case PROP_SOURCE:
      ncm_obj_dict_str_clear (&si->source);
      si->source = g_value_dup_boxed (value);
      break;
    case PROP_DATA_HASH:
      g_clear_pointer (&si->data_hash, g_free);
      si->data_hash = g_value_dup_string (value);
      break;
    case PROP_FILENAME:
    {
      const gchar *filename = g_value_get_string (value);

      g_clear_pointer (&si->filename, g_free);

      /* A relative name is relative to the file being read, or to the base directory */
      if ((filename != NULL) && (filename[0] != '\0'))
        si->filename = g_canonicalize_filename (filename, ncm_serialize_context_peek_dir (ncm_serialize_context_peek_current ()));

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gchar *_nc_galaxy_wl_obs_shape_interp_written_name (NcGalaxyWLObsShapeInterp *si);

static void
_nc_galaxy_wl_obs_shape_interp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObsShapeInterp *si = NC_GALAXY_WL_OBS_SHAPE_INTERP (object);

  switch (prop_id)
  {
    case PROP_INTERP:

      /* Backed by a file, the object is written without its interpolants */
      if (si->filename == NULL)
        g_value_set_boxed (value, si->interp);

      break;
    case PROP_SOURCE:
      g_value_set_boxed (value, si->source);
      break;
    case PROP_DATA_HASH:
      g_value_set_string (value, si->data_hash);
      break;
    case PROP_FILENAME:

      if (si->filename != NULL)
        g_value_take_string (value, _nc_galaxy_wl_obs_shape_interp_written_name (si));

      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_obs_shape_interp_constructed (GObject *object)
{
  NcGalaxyWLObsShapeInterp *si = NC_GALAXY_WL_OBS_SHAPE_INTERP (object);

  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_shape_interp_parent_class)->constructed (object);

  if (si->filename != NULL)
  {
    if (si->interp != NULL)
      g_error ("nc_galaxy_wl_obs_shape_interp: an object backed by the file %s also holds interpolants.", si->filename);

    if (si->data_hash == NULL)
      g_error ("nc_galaxy_wl_obs_shape_interp: an object backed by the file %s has no data hash.", si->filename);
  }
  else if (si->interp == NULL)
  {
    si->interp = ncm_obj_dict_int_new ();
  }
}

static void
_nc_galaxy_wl_obs_shape_interp_dispose (GObject *object)
{
  NcGalaxyWLObsShapeInterp *si = NC_GALAXY_WL_OBS_SHAPE_INTERP (object);

  ncm_obj_dict_int_clear (&si->interp);
  ncm_obj_dict_str_clear (&si->source);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_shape_interp_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_obs_shape_interp_finalize (GObject *object)
{
  NcGalaxyWLObsShapeInterp *si = NC_GALAXY_WL_OBS_SHAPE_INTERP (object);

  g_clear_pointer (&si->data_hash, g_free);
  g_clear_pointer (&si->filename, g_free);
  g_mutex_clear (&si->lock);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_shape_interp_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_obs_shape_interp_class_init (NcGalaxyWLObsShapeInterpClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_obs_shape_interp_set_property;
  object_class->get_property = &_nc_galaxy_wl_obs_shape_interp_get_property;
  object_class->constructed  = &_nc_galaxy_wl_obs_shape_interp_constructed;
  object_class->dispose      = &_nc_galaxy_wl_obs_shape_interp_dispose;
  object_class->finalize     = &_nc_galaxy_wl_obs_shape_interp_finalize;

  /**
   * NcGalaxyWLObsShapeInterp:interp:
   *
   * The interpolants, a #NcmSpectralTensor for each row that has one, keyed by the row;
   * %NULL for an object backed by a file.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERP,
                                   g_param_spec_boxed ("interp",
                                                       "Interpolants",
                                                       "Interpolants of ln P in g of the galaxies",
                                                       NCM_TYPE_OBJ_DICT_INT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsShapeInterp:source:
   *
   * The record of what built the interpolants, objects keyed by their role (see
   * nc_galaxy_shape_factor_spectral_set_obs_source()).
   */
  g_object_class_install_property (object_class,
                                   PROP_SOURCE,
                                   g_param_spec_boxed ("source",
                                                       "Source",
                                                       "Record of what built the interpolants",
                                                       NCM_TYPE_OBJ_DICT_STR,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsShapeInterp:data-hash:
   *
   * The data hash (ncm_catalog_get_data_hash()) of the rows the interpolants describe.
   */
  g_object_class_install_property (object_class,
                                   PROP_DATA_HASH,
                                   g_param_spec_string ("data-hash",
                                                        "Data hash",
                                                        "Data hash of the rows the interpolants describe",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsShapeInterp:filename:
   *
   * The file that holds the interpolants, %NULL when they are in memory. A relative name
   * is resolved against the directory of ncm_serialize_context_peek_current(). Read, the
   * property gives the bare file name inside a file function of #NcmSerialize writing to
   * the same directory, aborts inside one writing to another directory, and gives the
   * absolute name otherwise.
   */
  g_object_class_install_property (object_class,
                                   PROP_FILENAME,
                                   g_param_spec_string ("filename",
                                                        "File name",
                                                        "File that holds the interpolants",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));
}

/*
 * The name of the file of si as written: the bare file name when the current
 * NcmSerializeContext is a file in the same directory, the absolute name when it is the base
 * context. A file of another directory aborts.
 */
static gchar *
_nc_galaxy_wl_obs_shape_interp_written_name (NcGalaxyWLObsShapeInterp *si)
{
  NcmSerializeContext *ctx = ncm_serialize_context_peek_current ();
  gchar *dir;
  gchar *name;

  if (!ncm_serialize_context_has_file (ctx))
    return g_strdup (si->filename);

  dir = g_path_get_dirname (si->filename);

  if (!g_str_equal (dir, ncm_serialize_context_peek_dir (ctx)))
    g_error ("nc_galaxy_wl_obs_shape_interp: the interpolants are in %s, which is not in %s, the "
             "directory of the file being written.", si->filename, ncm_serialize_context_peek_dir (ctx));

  name = g_path_get_basename (si->filename);
  g_free (dir);

  return name;
}

/*
 * Reads the interpolants from the file of si, once: the file must hold an object of this
 * type in memory with the data hash of si.
 */
static void
_nc_galaxy_wl_obs_shape_interp_load (NcGalaxyWLObsShapeInterp *si)
{
  if ((si->filename == NULL) || g_atomic_int_get (&si->loaded))
    return;

  g_mutex_lock (&si->lock);

  if (!g_atomic_int_get (&si->loaded))
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
    GObject *obj      = ncm_serialize_from_binfile (ser, si->filename);
    NcGalaxyWLObsShapeInterp *file_si;

    if (!NC_IS_GALAXY_WL_OBS_SHAPE_INTERP (obj))
      g_error ("nc_galaxy_wl_obs_shape_interp: the file %s holds a %s, not a NcGalaxyWLObsShapeInterp.",
               si->filename, G_OBJECT_TYPE_NAME (obj));

    file_si = NC_GALAXY_WL_OBS_SHAPE_INTERP (obj);

    if (file_si->filename != NULL)
      g_error ("nc_galaxy_wl_obs_shape_interp: the file %s refers to the file %s instead of holding the interpolants.",
               si->filename, file_si->filename);

    if (g_strcmp0 (file_si->data_hash, si->data_hash) != 0)
      g_error ("nc_galaxy_wl_obs_shape_interp: the file %s holds the interpolants of the data hash %s, not %s.",
               si->filename, file_si->data_hash, si->data_hash);

    si->interp = ncm_obj_dict_int_ref (file_si->interp);

    nc_galaxy_wl_obs_shape_interp_free (file_si);
    ncm_serialize_free (ser);

    g_atomic_int_set (&si->loaded, TRUE);
  }

  g_mutex_unlock (&si->lock);
}

/* Reads the file of si, when it has one, and keeps the interpolants in memory only */
static void
_nc_galaxy_wl_obs_shape_interp_to_memory (NcGalaxyWLObsShapeInterp *si)
{
  _nc_galaxy_wl_obs_shape_interp_load (si);
  g_clear_pointer (&si->filename, g_free);
}

/**
 * nc_galaxy_wl_obs_shape_interp_new:
 *
 * Creates an object with no interpolants in memory, no record and no data hash.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObsShapeInterp.
 */
NcGalaxyWLObsShapeInterp *
nc_galaxy_wl_obs_shape_interp_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP, NULL);
}

/**
 * nc_galaxy_wl_obs_shape_interp_new_from_file:
 * @filename: the file that holds the interpolants
 * @data_hash: the data hash of the rows they describe
 * @source: the record of what built them
 *
 * Creates an object backed by the file @filename, the binary serialization of a
 * #NcGalaxyWLObsShapeInterp in memory with the data hash @data_hash. The file is read at
 * the first call of nc_galaxy_wl_obs_shape_interp_peek_interp(). A relative @filename is
 * resolved against the directory of ncm_serialize_context_peek_current().
 *
 * Returns: (transfer full): a new #NcGalaxyWLObsShapeInterp.
 */
NcGalaxyWLObsShapeInterp *
nc_galaxy_wl_obs_shape_interp_new_from_file (const gchar *filename, const gchar *data_hash, NcmObjDictStr *source)
{
  return g_object_new (NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP,
                       "filename", filename,
                       "data-hash", data_hash,
                       "source", source,
                       NULL);
}

/**
 * nc_galaxy_wl_obs_shape_interp_ref:
 * @si: a #NcGalaxyWLObsShapeInterp
 *
 * Increases the reference count of @si.
 *
 * Returns: (transfer full): @si.
 */
NcGalaxyWLObsShapeInterp *
nc_galaxy_wl_obs_shape_interp_ref (NcGalaxyWLObsShapeInterp *si)
{
  return g_object_ref (si);
}

/**
 * nc_galaxy_wl_obs_shape_interp_free:
 * @si: a #NcGalaxyWLObsShapeInterp
 *
 * Decreases the reference count of @si.
 */
void
nc_galaxy_wl_obs_shape_interp_free (NcGalaxyWLObsShapeInterp *si)
{
  g_object_unref (si);
}

/**
 * nc_galaxy_wl_obs_shape_interp_clear:
 * @si: a #NcGalaxyWLObsShapeInterp
 *
 * Decreases the reference count of *@si when it is not %NULL, and sets *@si to %NULL.
 */
void
nc_galaxy_wl_obs_shape_interp_clear (NcGalaxyWLObsShapeInterp **si)
{
  g_clear_object (si);
}

/**
 * nc_galaxy_wl_obs_shape_interp_set_interp:
 * @si: a #NcGalaxyWLObsShapeInterp
 * @i: row index
 * @st: the interpolant of row @i
 *
 * Sets the interpolant of row @i, referencing @st. An object backed by a file reads it
 * first and then holds the interpolants in memory.
 */
void
nc_galaxy_wl_obs_shape_interp_set_interp (NcGalaxyWLObsShapeInterp *si, const guint i, NcmSpectralTensor *st)
{
  _nc_galaxy_wl_obs_shape_interp_to_memory (si);
  ncm_obj_dict_int_set (si->interp, i, G_OBJECT (st));
}

/**
 * nc_galaxy_wl_obs_shape_interp_peek_interp:
 * @si: a #NcGalaxyWLObsShapeInterp
 * @i: row index
 *
 * Gets the interpolant of row @i, reading the file of @si at the first call. A file that
 * cannot be read, that holds another type of object, that refers to another file or whose
 * data hash differs from the one of @si aborts with an error, and so does an object of
 * another type stored for row @i.
 *
 * Returns: (transfer none) (nullable): the interpolant of row @i, or %NULL when the row has
 * none.
 */
NcmSpectralTensor *
nc_galaxy_wl_obs_shape_interp_peek_interp (NcGalaxyWLObsShapeInterp *si, const guint i)
{
  GObject *interp;

  _nc_galaxy_wl_obs_shape_interp_load (si);

  interp = ncm_obj_dict_int_peek (si->interp, i);

  if ((interp != NULL) && !NCM_IS_SPECTRAL_TENSOR (interp))
    g_error ("nc_galaxy_wl_obs_shape_interp_peek_interp: the interpolant of row %u is a %s, not a NcmSpectralTensor.",
             i, G_OBJECT_TYPE_NAME (interp));

  return (NcmSpectralTensor *) interp;
}

/**
 * nc_galaxy_wl_obs_shape_interp_set_source:
 * @si: a #NcGalaxyWLObsShapeInterp
 * @source: the record of what built the interpolants
 *
 * Sets the record, referencing @source. An object backed by a file reads it first and then
 * holds the interpolants in memory.
 */
void
nc_galaxy_wl_obs_shape_interp_set_source (NcGalaxyWLObsShapeInterp *si, NcmObjDictStr *source)
{
  _nc_galaxy_wl_obs_shape_interp_to_memory (si);
  ncm_obj_dict_str_ref (source);
  ncm_obj_dict_str_clear (&si->source);
  si->source = source;
}

/**
 * nc_galaxy_wl_obs_shape_interp_peek_source:
 * @si: a #NcGalaxyWLObsShapeInterp
 *
 * Gets the record of what built the interpolants, without reading the file of @si.
 *
 * Returns: (transfer none) (nullable): the record, or %NULL when there is none.
 */
NcmObjDictStr *
nc_galaxy_wl_obs_shape_interp_peek_source (NcGalaxyWLObsShapeInterp *si)
{
  return si->source;
}

/**
 * nc_galaxy_wl_obs_shape_interp_set_data_hash:
 * @si: a #NcGalaxyWLObsShapeInterp
 * @data_hash: the data hash of the rows the interpolants describe
 *
 * Sets the data hash (ncm_catalog_get_data_hash()). An object backed by a file reads it
 * first and then holds the interpolants in memory.
 */
void
nc_galaxy_wl_obs_shape_interp_set_data_hash (NcGalaxyWLObsShapeInterp *si, const gchar *data_hash)
{
  _nc_galaxy_wl_obs_shape_interp_to_memory (si);
  g_free (si->data_hash);
  si->data_hash = g_strdup (data_hash);
}

/**
 * nc_galaxy_wl_obs_shape_interp_peek_data_hash:
 * @si: a #NcGalaxyWLObsShapeInterp
 *
 * Returns: (transfer none) (nullable): the data hash of the rows the interpolants describe,
 * or %NULL when it was not set.
 */
const gchar *
nc_galaxy_wl_obs_shape_interp_peek_data_hash (NcGalaxyWLObsShapeInterp *si)
{
  return si->data_hash;
}

/**
 * nc_galaxy_wl_obs_shape_interp_peek_filename:
 * @si: a #NcGalaxyWLObsShapeInterp
 *
 * Returns: (transfer none) (nullable): the absolute name of the file that holds the
 * interpolants, or %NULL when they are in memory.
 */
const gchar *
nc_galaxy_wl_obs_shape_interp_peek_filename (NcGalaxyWLObsShapeInterp *si)
{
  return si->filename;
}


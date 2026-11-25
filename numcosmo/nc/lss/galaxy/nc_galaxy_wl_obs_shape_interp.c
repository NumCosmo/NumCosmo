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
 * An object created by nc_galaxy_wl_obs_shape_interp_create_empty() writes its file:
 * nc_galaxy_wl_obs_shape_interp_set_interp() appends each interpolant to a temporary file,
 * and nc_galaxy_wl_obs_shape_interp_close() writes the file from it, after which the object
 * is backed by the file. The interpolants are never held in memory together.
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

#ifndef NUMCOSMO_GIR_SCAN
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <glib/gstdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyWLObsShapeInterp
{
  /*< private >*/
  GObject parent_instance;
  NcmObjDictInt *interp; /* row -> NcmSpectralTensor, NULL until a file is read */
  NcmObjDictStr *source;
  gchar *data_hash;
  gchar *filename;     /* absolute, NULL in memory */
  gint loaded;         /* whether the file was read, read atomically */
  FILE *records;       /* open between create_empty and close, NULL otherwise */
  gchar *records_name; /* filename with the suffix .records while writing */
  gboolean write_failed;
  GMutex lock; /* serializes the read of the file and every set_interp */
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
  si->interp       = NULL;
  si->source       = NULL;
  si->data_hash    = NULL;
  si->filename     = NULL;
  si->loaded       = FALSE;
  si->records      = NULL;
  si->records_name = NULL;
  si->write_failed = FALSE;
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

      if (si->records != NULL)
        g_error ("nc_galaxy_wl_obs_shape_interp: the file %s is being written; "
                 "nc_galaxy_wl_obs_shape_interp_close() must run first.", si->filename);

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

  /* An object never closed writes no file */
  if (si->records != NULL)
  {
    fclose (si->records);
    g_remove (si->records_name);
  }

  g_clear_pointer (&si->records_name, g_free);
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

  if (si->records != NULL)
    g_error ("nc_galaxy_wl_obs_shape_interp: the file %s is being written; "
             "nc_galaxy_wl_obs_shape_interp_close() must run first.", si->filename);

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

/*
 * Reads the file of si, when it has one, and keeps the interpolants in memory only. An object
 * being written stays in that state.
 */
static void
_nc_galaxy_wl_obs_shape_interp_to_memory (NcGalaxyWLObsShapeInterp *si)
{
  if (si->records != NULL)
    return;

  _nc_galaxy_wl_obs_shape_interp_load (si);

  g_mutex_lock (&si->lock);
  g_clear_pointer (&si->filename, g_free);
  g_mutex_unlock (&si->lock);
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
 * nc_galaxy_wl_obs_shape_interp_create_empty:
 * @filename: the file to write
 *
 * Creates an object that writes its interpolants to @filename, for interpolants that should
 * not be held in memory together. Each nc_galaxy_wl_obs_shape_interp_set_interp() appends
 * the interpolant as one record (its row, its byte length and its #GVariant) to a temporary
 * file, @filename with the suffix `.records`; an existing temporary file is overwritten.
 * The record of what built the interpolants and the data hash are set with
 * nc_galaxy_wl_obs_shape_interp_set_source() and
 * nc_galaxy_wl_obs_shape_interp_set_data_hash(), and
 * nc_galaxy_wl_obs_shape_interp_close() writes @filename. Until then, reading an
 * interpolant or serializing the object aborts with an error. An object freed before
 * nc_galaxy_wl_obs_shape_interp_close() removes the temporary file and writes nothing. A
 * relative @filename is resolved against the directory of
 * ncm_serialize_context_peek_current(). A temporary file that cannot be created aborts
 * with an error.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObsShapeInterp.
 */
NcGalaxyWLObsShapeInterp *
nc_galaxy_wl_obs_shape_interp_create_empty (const gchar *filename)
{
  NcGalaxyWLObsShapeInterp *si = nc_galaxy_wl_obs_shape_interp_new ();

  ncm_obj_dict_int_clear (&si->interp);
  si->filename     = g_canonicalize_filename (filename, ncm_serialize_context_peek_dir (ncm_serialize_context_peek_current ()));
  si->records_name = g_strconcat (si->filename, ".records", NULL);
  si->records      = g_fopen (si->records_name, "wb");

  if (si->records == NULL)
    g_error ("nc_galaxy_wl_obs_shape_interp_create_empty: cannot create %s: %s.", si->records_name, g_strerror (errno));

  return si;
}

typedef struct _NcGalaxyWLObsShapeInterpRecord
{
  guint64 row;
  guint64 seq;
  gsize off;
  gsize len;
} NcGalaxyWLObsShapeInterpRecord;

static gint
_nc_galaxy_wl_obs_shape_interp_record_cmp (gconstpointer a, gconstpointer b)
{
  const NcGalaxyWLObsShapeInterpRecord *ra = a;
  const NcGalaxyWLObsShapeInterpRecord *rb = b;

  if (ra->row != rb->row)
    return (ra->row < rb->row) ? -1 : 1;

  return (ra->seq < rb->seq) ? -1 : ((ra->seq > rb->seq) ? 1 : 0);
}

/**
 * nc_galaxy_wl_obs_shape_interp_close:
 * @si: a #NcGalaxyWLObsShapeInterp created by nc_galaxy_wl_obs_shape_interp_create_empty()
 *
 * Writes the file of @si, the binary serialization of a #NcGalaxyWLObsShapeInterp in memory
 * with the interpolants of the records, the record of what built them and the data hash of
 * @si, and removes the temporary file. For a row set more than once, the last interpolant
 * is written. @si is then backed by the file (nc_galaxy_wl_obs_shape_interp_new_from_file()).
 * The records are mapped and become the entries of the interpolants, in row order, without
 * being parsed, and the file is written through a second mapping, as @filename with the
 * suffix `.part` renamed to @filename: neither the interpolants nor the bytes of the file are
 * allocated, and an existing @filename is replaced, never truncated. Every
 * nc_galaxy_wl_obs_shape_interp_set_interp() must have returned. An object not created by
 * nc_galaxy_wl_obs_shape_interp_create_empty() or already closed, an object without a data
 * hash, a records file that does not hold whole records and a file that cannot be written
 * abort with an error.
 */
void
nc_galaxy_wl_obs_shape_interp_close (NcGalaxyWLObsShapeInterp *si)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GArray *recs      = g_array_new (FALSE, FALSE, sizeof (NcGalaxyWLObsShapeInterpRecord));
  GError *error     = NULL;
  GVariantBuilder dict, props;
  GMappedFile *map;
  GVariant *final;
  const guint8 *base;
  gchar *part_name;
  GBytes *all;
  gsize size = 0;
  gsize off  = 0;
  guint64 seq;
  guint i;

  if (si->records == NULL)
    g_error ("nc_galaxy_wl_obs_shape_interp_close: the object is not being written.");

  if (si->data_hash == NULL)
    g_error ("nc_galaxy_wl_obs_shape_interp_close: the object writing %s has no data hash.", si->filename);

  if ((fclose (si->records) != 0) || si->write_failed)
    g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot write %s.", si->records_name);

  si->records = NULL;
  map         = g_mapped_file_new (si->records_name, FALSE, &error);

  if (map == NULL)
    g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot read %s: %s.", si->records_name, error->message);

  all  = g_mapped_file_get_bytes (map);
  base = g_bytes_get_data (all, &size);

  for (seq = 0; off < size; seq++)
  {
    NcGalaxyWLObsShapeInterpRecord rec;
    guint64 head[2];

    if (size - off < sizeof (head))
      g_error ("nc_galaxy_wl_obs_shape_interp_close: %s ends inside a record.", si->records_name);

    memcpy (head, base + off, sizeof (head));
    off += sizeof (head);

    if ((head[0] > G_MAXINT32) || (head[1] > size - off))
      g_error ("nc_galaxy_wl_obs_shape_interp_close: %s has a record of row %" G_GUINT64_FORMAT " and %" G_GUINT64_FORMAT " bytes, "
               "outside the rows or the file.", si->records_name, head[0], head[1]);

    rec.row = head[0];
    rec.seq = seq;
    rec.off = off;
    rec.len = head[1];
    g_array_append_val (recs, rec);

    off += head[1] + (8 - head[1] % 8) % 8;
  }

  /* Row order; of the records of one row, the last written is kept */
  g_array_sort (recs, &_nc_galaxy_wl_obs_shape_interp_record_cmp);
  g_variant_builder_init (&dict, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_INT_TYPE));

  for (i = 0; i < recs->len; i++)
  {
    const NcGalaxyWLObsShapeInterpRecord *rec = &g_array_index (recs, NcGalaxyWLObsShapeInterpRecord, i);

    if ((i + 1 == recs->len) || (g_array_index (recs, NcGalaxyWLObsShapeInterpRecord, i + 1).row != rec->row))
    {
      GBytes *bytes   = g_bytes_new_from_bytes (all, rec->off, rec->len);
      GVariant *child = g_variant_new_from_bytes (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), bytes, TRUE);

      g_variant_builder_add_value (&dict, g_variant_new_dict_entry (g_variant_new_int32 ((gint32) rec->row), child));
      g_bytes_unref (bytes);
    }
  }

  g_variant_builder_init (&props, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE));
  g_variant_builder_add (&props, "{sv}", "interp", g_variant_builder_end (&dict));

  if (si->source != NULL)
  {
    GVariant *source_var = ncm_serialize_dict_str_to_variant (ser, si->source);

    g_variant_builder_add (&props, "{sv}", "source", source_var);
    g_variant_unref (source_var);
  }

  g_variant_builder_add (&props, "{sv}", "data-hash", g_variant_new_string (si->data_hash));
  final = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_FORMAT, G_OBJECT_TYPE_NAME (si), g_variant_builder_end (&props)));

  /*
   * The variant is serialized into a file mapped at its final size: the records are copied
   * from one mapped file to the other (g_variant_get_data() would allocate all the bytes).
   */
  part_name = g_strconcat (si->filename, ".part", NULL);
  {
    const gsize out_size = g_variant_get_size (final);
    const gint fd        = g_open (part_name, O_RDWR | O_CREAT | O_TRUNC, 0666);
    gpointer out;

    if (fd < 0)
      g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot create %s: %s.", part_name, g_strerror (errno));

    if (ftruncate (fd, (off_t) out_size) != 0)
      g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot set the size of %s: %s.", part_name, g_strerror (errno));

    out = mmap (NULL, out_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

    if (out == MAP_FAILED)
      g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot map %s: %s.", part_name, g_strerror (errno));

    g_variant_store (final, out);

    if ((msync (out, out_size, MS_SYNC) != 0) || (munmap (out, out_size) != 0) || (close (fd) != 0))
      g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot write %s: %s.", part_name, g_strerror (errno));
  }

  /* An existing file is replaced, never truncated: another object may have it mapped */
  if (g_rename (part_name, si->filename) != 0)
    g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot rename %s to %s: %s.", part_name, si->filename, g_strerror (errno));

  g_variant_unref (final);
  g_bytes_unref (all);
  g_mapped_file_unref (map);

  if (g_remove (si->records_name) != 0)
    g_error ("nc_galaxy_wl_obs_shape_interp_close: cannot remove %s: %s.", si->records_name, g_strerror (errno));

  g_clear_pointer (&si->records_name, g_free);
  g_array_unref (recs);
  g_free (part_name);
  ncm_serialize_free (ser);
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
 * Sets the interpolant of row @i, replacing the one it had. An object in memory references
 * @st. An object backed by a file reads it first and then holds the interpolants in memory.
 * An object created by nc_galaxy_wl_obs_shape_interp_create_empty() appends the
 * serialization of @st to its records file and does not keep @st. Calls on the same object
 * may run in parallel.
 */
void
nc_galaxy_wl_obs_shape_interp_set_interp (NcGalaxyWLObsShapeInterp *si, const guint i, NcmSpectralTensor *st)
{
  if (si->records != NULL)
  {
    /* Serialized outside the lock: only the write is one thread at a time */
    NcmSerialize *ser       = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    GVariant *var           = ncm_serialize_to_variant (ser, G_OBJECT (st));
    const guint64 head[2]   = {(guint64) i, (guint64) g_variant_get_size (var)};
    const gsize pad         = (8 - head[1] % 8) % 8;
    const gchar zeros[8]    = {0};
    const gconstpointer buf = g_variant_get_data (var);

    g_mutex_lock (&si->lock);

    if ((fwrite (head, sizeof (guint64), 2, si->records) != 2) || (fwrite (buf, 1, head[1], si->records) != head[1]) ||
        (fwrite (zeros, 1, pad, si->records) != pad))
      si->write_failed = TRUE;

    g_mutex_unlock (&si->lock);

    g_variant_unref (var);
    ncm_serialize_free (ser);
  }
  else
  {
    _nc_galaxy_wl_obs_shape_interp_to_memory (si);

    g_mutex_lock (&si->lock);
    ncm_obj_dict_int_set (si->interp, i, G_OBJECT (st));
    g_mutex_unlock (&si->lock);
  }
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


/***************************************************************************
 *            ncm_mset_catalog_storage_fits.c
 *
 *  Sat Oct 31 14:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_storage_fits.c
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
 * NcmMSetCatalogStorageFits:
 *
 * An #NcmMSetCatalogStorage backed by a FITS file, the file format of #NcmMSetCatalog.
 *
 * The file holds the serialized #NcmMSet in HDU0 and a binary table named
 * %NCM_MSET_CATALOG_EXTNAME with one column per catalog column: the additional values (the
 * weight last, when the catalog is weighted) and then the free parameters. The table header
 * holds the description (`RTYPE`, `NCHAINS`, `NADDVAL`, `WEIGHTED`, `M2LNP_ID`, `ASYMB<n>`,
 * `FSYMB<n>`) and the run metadata (`RNG_ALGO`, `RNG_SEED`, `RNG_INIS`, `RNG_STAT`,
 * `MARKID`, and the sampler segments as `NSEGS`, `SEGS<n>`, `SEGN<n>`, `SEGO<n>`; `SAMPLER`
 * and `SAMPOPT` repeat the last segment for the previous library).
 *
 * The file is opened at the first call that needs it. A storage on a file that does not
 * exist is uninitialized; ncm_mset_catalog_storage_set_desc() creates the file. A
 * read-only storage (#NcmMSetCatalogStorageFits:read-only) fails every write with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_READ_ONLY. Rows appended inside a recording are written to
 * the file at each sync; the last #NcmMSetCatalogStorageFits:keep-last rows per chain stay
 * in memory, and ncm_mset_catalog_storage_get_row() reads older rows from the file in blocks
 * of #NcmMSetCatalogStorageFits:read-block rows.
 *
 * Files written by the previous library have no `NSEGS` key. They are read through a
 * translation: the run type values "Ensemble Sampler MCMC", "bootstrap-mean", the
 * #NcmFitMC resampling nicks and the transition kernel names become the class names
 * %NCM_MSET_CATALOG_RTYPE_ESMCMC, %NCM_MSET_CATALOG_RTYPE_MCBS, %NCM_MSET_CATALOG_RTYPE_MC
 * and %NCM_MSET_CATALOG_RTYPE_MCMC, and the sampler segments are rebuilt from `INITSMP`,
 * `SAMPLER` and `SAMPOPT`. A read-only storage keeps the translation in memory; a writable
 * one writes it to the file at the first recording. `FIRST_ID` is read, so that `MARKID`,
 * written in the file's numbering, is converted to a row index; new files have `FIRST_ID`
 * 0. A file without `RNG_ALGO` and with an RNG state is read with the default algorithm of
 * ncm_rng_new(). Without CFITSIO every method fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NO_CFITSIO.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_catalog_storage_fits.h"
#include "ncm/fit/ncm_mset_catalog.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#ifdef HAVE_CFITSIO
#include <fitsio.h>
#endif /* HAVE_CFITSIO */
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmMSetCatalogStorageFitsSegment
{
  guint start;
  gchar *sampler;
  gchar *options;
} NcmMSetCatalogStorageFitsSegment;

struct _NcmMSetCatalogStorageFits
{
  NcmMSetCatalogStorage parent_instance;
  gchar *filename;
  gboolean read_only;
  guint keep_last;
  guint read_block;
  gpointer fptr;
  gboolean opened;
  gboolean legacy;
  NcmMSetCatalogDesc *desc;
  GArray *porder;
  gint first_id;
  GPtrArray *tail;
  guint tail_first;
  guint n_pending;
  GPtrArray *block;
  guint block_first;
  GArray *col_ws;
  gchar *rng_algo;
  gulong rng_seed;
  gchar *rng_inis;
  gchar *rng_stat;
  guint markovian_id;
  GArray *segments;
};

enum
{
  PROP_0,
  PROP_FILENAME,
  PROP_READ_ONLY,
  PROP_KEEP_LAST,
  PROP_READ_BLOCK,
};

G_DEFINE_TYPE (NcmMSetCatalogStorageFits, ncm_mset_catalog_storage_fits, NCM_TYPE_MSET_CATALOG_STORAGE)

static void
_ncm_mset_catalog_storage_fits_segment_clear (gpointer data)
{
  NcmMSetCatalogStorageFitsSegment *segment = data;

  g_clear_pointer (&segment->sampler, g_free);
  g_clear_pointer (&segment->options, g_free);
}

static void
ncm_mset_catalog_storage_fits_init (NcmMSetCatalogStorageFits *storage_fits)
{
  storage_fits->filename     = NULL;
  storage_fits->read_only    = FALSE;
  storage_fits->keep_last    = 0;
  storage_fits->read_block   = 0;
  storage_fits->fptr         = NULL;
  storage_fits->opened       = FALSE;
  storage_fits->legacy       = FALSE;
  storage_fits->desc         = NULL;
  storage_fits->porder       = g_array_new (FALSE, FALSE, sizeof (gint));
  storage_fits->first_id     = 0;
  storage_fits->tail         = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vector_free);
  storage_fits->tail_first   = 0;
  storage_fits->n_pending    = 0;
  storage_fits->block        = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vector_free);
  storage_fits->block_first  = 0;
  storage_fits->col_ws       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  storage_fits->rng_algo     = NULL;
  storage_fits->rng_seed     = 0;
  storage_fits->rng_inis     = NULL;
  storage_fits->rng_stat     = NULL;
  storage_fits->markovian_id = 0;
  storage_fits->segments     = g_array_new (FALSE, TRUE, sizeof (NcmMSetCatalogStorageFitsSegment));

  g_array_set_clear_func (storage_fits->segments, &_ncm_mset_catalog_storage_fits_segment_clear);
}

static void
_ncm_mset_catalog_storage_fits_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (object);

  g_return_if_fail (NCM_IS_MSET_CATALOG_STORAGE_FITS (object));

  switch (prop_id)
  {
    case PROP_FILENAME:
      storage_fits->filename = g_value_dup_string (value);
      break;
    case PROP_READ_ONLY:
      storage_fits->read_only = g_value_get_boolean (value);
      break;
    case PROP_KEEP_LAST:
      storage_fits->keep_last = g_value_get_uint (value);
      break;
    case PROP_READ_BLOCK:
      storage_fits->read_block = g_value_get_uint (value);
      g_ptr_array_set_size (storage_fits->block, 0);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_storage_fits_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (object);

  g_return_if_fail (NCM_IS_MSET_CATALOG_STORAGE_FITS (object));

  switch (prop_id)
  {
    case PROP_FILENAME:
      g_value_set_string (value, storage_fits->filename);
      break;
    case PROP_READ_ONLY:
      g_value_set_boolean (value, storage_fits->read_only);
      break;
    case PROP_KEEP_LAST:
      g_value_set_uint (value, storage_fits->keep_last);
      break;
    case PROP_READ_BLOCK:
      g_value_set_uint (value, storage_fits->read_block);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_mset_catalog_storage_fits_close (NcmMSetCatalogStorageFits *storage_fits);

static void
_ncm_mset_catalog_storage_fits_dispose (GObject *object)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (object);

  _ncm_mset_catalog_storage_fits_close (storage_fits);

  ncm_mset_catalog_desc_clear (&storage_fits->desc);
  g_clear_pointer (&storage_fits->tail, g_ptr_array_unref);
  g_clear_pointer (&storage_fits->block, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_fits_parent_class)->dispose (object);
}

static void
_ncm_mset_catalog_storage_fits_finalize (GObject *object)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (object);

  g_clear_pointer (&storage_fits->filename, g_free);
  g_clear_pointer (&storage_fits->porder, g_array_unref);
  g_clear_pointer (&storage_fits->col_ws, g_array_unref);
  g_clear_pointer (&storage_fits->rng_algo, g_free);
  g_clear_pointer (&storage_fits->rng_inis, g_free);
  g_clear_pointer (&storage_fits->rng_stat, g_free);
  g_clear_pointer (&storage_fits->segments, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_fits_parent_class)->finalize (object);
}

static NcmMSetCatalogDesc *_ncm_mset_catalog_storage_fits_peek_desc (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_fits_store_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);
static void _ncm_mset_catalog_storage_fits_bind_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);
static guint _ncm_mset_catalog_storage_fits_nrows (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_fits_append (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error);
static NcmVector *_ncm_mset_catalog_storage_fits_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error);
static void _ncm_mset_catalog_storage_fits_sync (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_fits_get_rng_init (NcmMSetCatalogStorage *storage, gchar **algo, gulong *seed, gchar **state, GError **error);
static void _ncm_mset_catalog_storage_fits_set_rng_init (NcmMSetCatalogStorage *storage, const gchar *algo, gulong seed, const gchar *state, GError **error);
static gchar *_ncm_mset_catalog_storage_fits_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_fits_set_rng_stat (NcmMSetCatalogStorage *storage, const gchar *state, GError **error);
static guint _ncm_mset_catalog_storage_fits_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_fits_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error);
static guint _ncm_mset_catalog_storage_fits_n_samplers (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_fits_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error);
static void _ncm_mset_catalog_storage_fits_set_sampler (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options, GError **error);

static void
ncm_mset_catalog_storage_fits_class_init (NcmMSetCatalogStorageFitsClass *klass)
{
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcmMSetCatalogStorageClass *storage_class = NCM_MSET_CATALOG_STORAGE_CLASS (klass);

  object_class->set_property = &_ncm_mset_catalog_storage_fits_set_property;
  object_class->get_property = &_ncm_mset_catalog_storage_fits_get_property;
  object_class->dispose      = &_ncm_mset_catalog_storage_fits_dispose;
  object_class->finalize     = &_ncm_mset_catalog_storage_fits_finalize;

  /**
   * NcmMSetCatalogStorageFits:filename:
   *
   * The FITS file. It is opened at the first call that needs it and created by
   * ncm_mset_catalog_storage_set_desc() when it does not exist.
   */
  g_object_class_install_property (object_class,
                                   PROP_FILENAME,
                                   g_param_spec_string ("filename",
                                                        NULL,
                                                        "Catalog file name",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageFits:read-only:
   *
   * Whether the file is opened for reading only; every write then fails with
   * %NCM_MSET_CATALOG_STORAGE_ERROR_READ_ONLY, and a file that does not exist is an error.
   */
  g_object_class_install_property (object_class,
                                   PROP_READ_ONLY,
                                   g_param_spec_boolean ("read-only",
                                                         NULL,
                                                         "Whether the file is read-only",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageFits:keep-last:
   *
   * Rows per chain kept in memory after a sync; older rows are read from the file.
   */
  g_object_class_install_property (object_class,
                                   PROP_KEEP_LAST,
                                   g_param_spec_uint ("keep-last",
                                                      NULL,
                                                      "Rows per chain kept in memory",
                                                      1, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageFits:read-block:
   *
   * Number of rows read from the file at once by ncm_mset_catalog_storage_get_row().
   */
  g_object_class_install_property (object_class,
                                   PROP_READ_BLOCK,
                                   g_param_spec_uint ("read-block",
                                                      NULL,
                                                      "Number of rows read from the file at once",
                                                      1, G_MAXUINT, 1024,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  storage_class->peek_desc        = &_ncm_mset_catalog_storage_fits_peek_desc;
  storage_class->store_desc       = &_ncm_mset_catalog_storage_fits_store_desc;
  storage_class->bind_desc        = &_ncm_mset_catalog_storage_fits_bind_desc;
  storage_class->nrows            = &_ncm_mset_catalog_storage_fits_nrows;
  storage_class->append           = &_ncm_mset_catalog_storage_fits_append;
  storage_class->get_row          = &_ncm_mset_catalog_storage_fits_get_row;
  storage_class->sync             = &_ncm_mset_catalog_storage_fits_sync;
  storage_class->get_rng_init     = &_ncm_mset_catalog_storage_fits_get_rng_init;
  storage_class->set_rng_init     = &_ncm_mset_catalog_storage_fits_set_rng_init;
  storage_class->get_rng_stat     = &_ncm_mset_catalog_storage_fits_get_rng_stat;
  storage_class->set_rng_stat     = &_ncm_mset_catalog_storage_fits_set_rng_stat;
  storage_class->get_markovian_id = &_ncm_mset_catalog_storage_fits_get_markovian_id;
  storage_class->set_markovian_id = &_ncm_mset_catalog_storage_fits_set_markovian_id;
  storage_class->n_samplers       = &_ncm_mset_catalog_storage_fits_n_samplers;
  storage_class->get_sampler      = &_ncm_mset_catalog_storage_fits_get_sampler;
  storage_class->set_sampler      = &_ncm_mset_catalog_storage_fits_set_sampler;
}

#ifdef HAVE_CFITSIO

#define FPTR(storage_fits) ((fitsfile *) (storage_fits)->fptr)

/* Sets error (BACKEND) with the CFITSIO message when status is not zero. */
static void
_ncm_mset_catalog_storage_fits_check (gint status, const gchar *where, GError **error)
{
  if (status != 0)
  {
    gchar msg[FLEN_STATUS];

    fits_get_errstatus (status, msg);
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "%s: FITS: %s", where, msg);
  }
}

/* Sets error (READ_ONLY) when the storage is read-only. */
static void
_ncm_mset_catalog_storage_fits_require_writable (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  if (storage_fits->read_only)
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_READ_ONLY,
                 "the file `%s' is read-only", storage_fits->filename);
}

/*
 * Header key readers. A missing key leaves the output as given and sets missing (when not
 * NULL); any other failure sets error.
 */

static void
_ncm_mset_catalog_storage_fits_read_int (fitsfile *fptr, const gchar *key, gint *value, gboolean *missing, GError **error)
{
  gint status = 0;

  if (missing != NULL)
    *missing = FALSE;

  fits_read_key (fptr, TINT, (gchar *) key, value, NULL, &status);

  if (status == KEY_NO_EXIST)
  {
    if (missing != NULL)
      *missing = TRUE;

    return;
  }

  _ncm_mset_catalog_storage_fits_check (status, key, error);
}

static void
_ncm_mset_catalog_storage_fits_read_ulong (fitsfile *fptr, const gchar *key, gulong *value, gboolean *missing, GError **error)
{
  gint status = 0;
  glong lvalue;

  if (missing != NULL)
    *missing = FALSE;

  fits_read_key_lng (fptr, (gchar *) key, &lvalue, NULL, &status);

  if (status == KEY_NO_EXIST)
  {
    if (missing != NULL)
      *missing = TRUE;

    return;
  }

  _ncm_mset_catalog_storage_fits_check (status, key, error);

  if (status == 0)
    *value = (gulong) lvalue;
}

static void
_ncm_mset_catalog_storage_fits_read_bool (fitsfile *fptr, const gchar *key, gboolean *value, gboolean *missing, GError **error)
{
  gint status = 0;
  gint ivalue = 0;

  if (missing != NULL)
    *missing = FALSE;

  fits_read_key (fptr, TLOGICAL, (gchar *) key, &ivalue, NULL, &status);

  if (status == KEY_NO_EXIST)
  {
    if (missing != NULL)
      *missing = TRUE;

    return;
  }

  _ncm_mset_catalog_storage_fits_check (status, key, error);

  if (status == 0)
    *value = ivalue ? TRUE : FALSE;
}

/* Returns a newly allocated string, NULL when the key is missing. */
static gchar *
_ncm_mset_catalog_storage_fits_read_str (fitsfile *fptr, const gchar *key, GError **error)
{
  gint status = 0;
  gchar value[FLEN_VALUE];

  fits_read_key_str (fptr, (gchar *) key, value, NULL, &status);

  if (status == KEY_NO_EXIST)
    return NULL;

  _ncm_mset_catalog_storage_fits_check (status, key, error);

  return (status == 0) ? g_strdup (value) : NULL;
}

/* Returns a newly allocated long string, NULL when the key is missing. */
static gchar *
_ncm_mset_catalog_storage_fits_read_longstr (fitsfile *fptr, const gchar *key, GError **error)
{
  gint status  = 0;
  gchar *value = NULL;
  gchar *res   = NULL;

  fits_read_key_longstr (fptr, (gchar *) key, &value, NULL, &status);

  if (status == KEY_NO_EXIST)
    return NULL;

  _ncm_mset_catalog_storage_fits_check (status, key, error);

  if (status == 0)
  {
    res = g_strdup (value);
    fits_free_memory (value, &status);
  }

  return res;
}

/* Header key writers; each sets error (BACKEND) on failure. */

static void
_ncm_mset_catalog_storage_fits_write_int (fitsfile *fptr, const gchar *key, gint value, const gchar *comment, GError **error)
{
  gint status = 0;

  fits_update_key (fptr, TINT, (gchar *) key, &value, (gchar *) comment, &status);
  _ncm_mset_catalog_storage_fits_check (status, key, error);
}

static void
_ncm_mset_catalog_storage_fits_write_ulong (fitsfile *fptr, const gchar *key, gulong value, const gchar *comment, GError **error)
{
  gint status  = 0;
  glong lvalue = (glong) value;

  fits_update_key (fptr, TLONG, (gchar *) key, &lvalue, (gchar *) comment, &status);
  _ncm_mset_catalog_storage_fits_check (status, key, error);
}

static void
_ncm_mset_catalog_storage_fits_write_bool (fitsfile *fptr, const gchar *key, gboolean value, const gchar *comment, GError **error)
{
  gint status = 0;
  gint ivalue = value ? 1 : 0;

  fits_update_key (fptr, TLOGICAL, (gchar *) key, &ivalue, (gchar *) comment, &status);
  _ncm_mset_catalog_storage_fits_check (status, key, error);
}

static void
_ncm_mset_catalog_storage_fits_write_str (fitsfile *fptr, const gchar *key, const gchar *value, const gchar *comment, GError **error)
{
  gint status = 0;

  fits_update_key (fptr, TSTRING, (gchar *) key, (gchar *) value, (gchar *) comment, &status);
  _ncm_mset_catalog_storage_fits_check (status, key, error);
}

static void
_ncm_mset_catalog_storage_fits_write_longstr (fitsfile *fptr, const gchar *key, const gchar *value, const gchar *comment, GError **error)
{
  gint status = 0;

  fits_update_key_longstr (fptr, (gchar *) key, (gchar *) value, (gchar *) comment, &status);
  _ncm_mset_catalog_storage_fits_check (status, key, error);
}

static void
_ncm_mset_catalog_storage_fits_flush (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  gint status = 0;

  fits_flush_file (FPTR (storage_fits), &status);
  _ncm_mset_catalog_storage_fits_check (status, "flush", error);
}

/*
 * Reads the mset of a file: the HDU0 vardict (`model-set' entry) or bare object, or the
 * `.mset' sidecar of files with an empty HDU0. Sets error (BACKEND) when none is found.
 */
static NcmMSet *
_ncm_mset_catalog_storage_fits_read_mset (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  fitsfile *fptr = FPTR (storage_fits);
  NcmMSet *mset  = NULL;
  gint status    = 0;
  gint hdutype   = 0;
  gint bitpix    = 0;
  gint naxis     = 0;
  glong naxes[1] = { 0 };

  fits_movabs_hdu (fptr, 1, &hdutype, &status);
  _ncm_mset_catalog_storage_fits_check (status, "HDU0", error);

  if (status == 0)
  {
    fits_get_img_param (fptr, 1, &bitpix, &naxis, naxes, &status);
    _ncm_mset_catalog_storage_fits_check (status, "HDU0", error);
  }

  if ((status == 0) && (naxis == 1) && (naxes[0] > 0))
  {
    const glong len = naxes[0];
    guint8 *data    = g_malloc (len);
    gint any_null   = 0;
    gchar *format;

    fits_read_img (fptr, TBYTE, 1, len, NULL, data, &any_null, &status);
    _ncm_mset_catalog_storage_fits_check (status, "HDU0", error);

    format = _ncm_mset_catalog_storage_fits_read_str (fptr, NCM_MSET_CATALOG_MSET_FORMAT_LABEL, error);

    if (format == NULL)
      format = g_strdup (NCM_MSET_CATALOG_MSET_FORMAT_OBJECT);

    if (*error == NULL)
    {
      NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

      if (strcmp (format, NCM_MSET_CATALOG_MSET_FORMAT_VARDICT) == 0)
      {
        GVariant *var = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE),
                                                 data, len, TRUE, g_free, data);
        NcmVarDict *vd    = ncm_serialize_var_dict_from_variant (ser, var);
        GObject *mset_obj = NULL;

        data = NULL;

        if (ncm_var_dict_get_object (vd, NCM_MSET_CATALOG_HDU0_MSET_KEY, ser, &mset_obj) && NCM_IS_MSET (mset_obj))
          mset = NCM_MSET (mset_obj);
        else
          g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                       "HDU0 of `%s' has no `%s' entry", storage_fits->filename, NCM_MSET_CATALOG_HDU0_MSET_KEY);

        ncm_var_dict_unref (vd);
        g_variant_unref (var);
      }
      else if (strcmp (format, NCM_MSET_CATALOG_MSET_FORMAT_OBJECT) == 0)
      {
        GVariant *var = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE),
                                                 data, len, TRUE, g_free, data);
        GObject *mset_obj = ncm_serialize_from_variant (ser, var);

        data = NULL;

        if (NCM_IS_MSET (mset_obj))
          mset = NCM_MSET (mset_obj);
        else
          g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                       "HDU0 of `%s' holds no NcmMSet", storage_fits->filename);

        g_variant_unref (var);
      }
      else
      {
        g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                     "unknown HDU0 format `%s' in `%s'", format, storage_fits->filename);
      }

      ncm_serialize_free (ser);
    }

    g_free (data);
    g_free (format);
  }
  else if (status == 0)
  {
    gchar *base_name = ncm_util_basename_fits (storage_fits->filename);
    gchar *mset_file = g_strdup_printf ("%s.mset", base_name);

    if (g_file_test (mset_file, G_FILE_TEST_EXISTS))
    {
      NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

      mset = ncm_mset_load (mset_file, ser, error);
      ncm_serialize_free (ser);
    }
    else
    {
      g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                   "`%s' has no embedded mset and no `%s' file", storage_fits->filename, mset_file);
    }

    g_free (mset_file);
    g_free (base_name);
  }

  if (*error == NULL)
  {
    fits_movnam_hdu (fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
    _ncm_mset_catalog_storage_fits_check (status, NCM_MSET_CATALOG_EXTNAME, error);
  }

  if ((*error != NULL) && (mset != NULL))
    ncm_mset_clear (&mset);

  return mset;
}

/* Appends a segment to the segment record. */
static void
_ncm_mset_catalog_storage_fits_add_segment (NcmMSetCatalogStorageFits *storage_fits, guint start, const gchar *sampler, const gchar *options)
{
  NcmMSetCatalogStorageFitsSegment segment = { start, g_strdup (sampler), g_strdup ((options != NULL) ? options : "") };

  g_array_append_val (storage_fits->segments, segment);
}

/*
 * The translation of files written by the previous library (no NSEGS key): the run type
 * values of each producer become the class names, and the sampler record becomes the
 * segments. Returns the run type of the Desc (newly allocated) and fills the segments.
 */
static gchar *
_ncm_mset_catalog_storage_fits_translate_legacy (NcmMSetCatalogStorageFits *storage_fits, const gchar *rtype, guint nchains, guint nrows, GError **error)
{
  static const struct
  {
    const gchar *old_rtype;
    const gchar *rtype;
  } table[] = {
    { "Ensemble Sampler MCMC", NCM_MSET_CATALOG_RTYPE_ESMCMC },
    { "bootstrap-mean", NCM_MSET_CATALOG_RTYPE_MCBS },
    { "from-model", NCM_MSET_CATALOG_RTYPE_MC },
    { "bootstrap-nomix", NCM_MSET_CATALOG_RTYPE_MC },
    { "bootstrap-mix", NCM_MSET_CATALOG_RTYPE_MC },
    { "Catalog Sampler", NCM_MSET_CATALOG_RTYPE_MCMC },
    { "Multivariate Flat Sampler", NCM_MSET_CATALOG_RTYPE_MCMC },
    { "Multivariate Gaussian Sampler", NCM_MSET_CATALOG_RTYPE_MCMC },
  };

  fitsfile *fptr    = FPTR (storage_fits);
  const gchar *rnew = NULL;
  gchar *sampler    = _ncm_mset_catalog_storage_fits_read_longstr (fptr, NCM_MSET_CATALOG_SAMPLER_LABEL, error);
  gchar *sampopt    = _ncm_mset_catalog_storage_fits_read_longstr (fptr, NCM_MSET_CATALOG_SAMPLER_OPTS_LABEL, error);
  gchar *initsmp    = _ncm_mset_catalog_storage_fits_read_longstr (fptr, NCM_MSET_CATALOG_INIT_SAMPLER_LABEL, error);
  gchar *res;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (table); i++)
  {
    if (strcmp (rtype, table[i].old_rtype) == 0)
    {
      rnew = table[i].rtype;
      break;
    }
  }

  if (*error == NULL)
  {
    if ((rnew != NULL) && (strcmp (rnew, NCM_MSET_CATALOG_RTYPE_ESMCMC) == 0))
    {
      if (initsmp != NULL)
        _ncm_mset_catalog_storage_fits_add_segment (storage_fits, 0, initsmp, NULL);

      if (sampler != NULL)
        _ncm_mset_catalog_storage_fits_add_segment (storage_fits, (initsmp != NULL) ? MIN (nchains, nrows) : 0, sampler, sampopt);
    }
    else if (rnew != NULL)
    {
      _ncm_mset_catalog_storage_fits_add_segment (storage_fits, 0, rtype, NULL);
    }
    else if (sampler != NULL)
    {
      _ncm_mset_catalog_storage_fits_add_segment (storage_fits, 0, sampler, sampopt);
    }
  }

  res = g_strdup ((rnew != NULL) ? rnew : rtype);

  g_free (sampler);
  g_free (sampopt);
  g_free (initsmp);

  return res;
}

/* Fills porder with the file column (1-based) of each column of desc. */
static void
_ncm_mset_catalog_storage_fits_set_porder (NcmMSetCatalogStorageFits *storage_fits, NcmMSetCatalogDesc *desc, GError **error)
{
  fitsfile *fptr    = FPTR (storage_fits);
  const guint ncols = ncm_mset_catalog_desc_ncols (desc);
  guint i;

  g_array_set_size (storage_fits->porder, ncols);

  for (i = 0; (i < ncols) && (*error == NULL); i++)
  {
    const gchar *name = ncm_mset_catalog_desc_col_full_name (desc, i);
    gint status       = 0;
    gint colnum       = 0;

    fits_get_colnum (fptr, CASESEN, (gchar *) name, &colnum, &status);

    if (status == COL_NOT_FOUND)
      g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                   "`%s' has no column `%s'", storage_fits->filename, name);
    else
      _ncm_mset_catalog_storage_fits_check (status, name, error);

    g_array_index (storage_fits->porder, gint, i) = colnum;
  }
}

/*
 * Opens an existing file and reads its description and run metadata. The table columns
 * are the additional values (NADDVAL of them, the weight last when WEIGHTED) and then the
 * free parameters of the HDU0 mset, each found by its full name.
 */
static void
_ncm_mset_catalog_storage_fits_open_existing (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  fitsfile *fptr    = NULL;
  gint status       = 0;
  gint nchains      = 0;
  gint nadd_vals    = 0;
  gint m2lnp_var    = -1;
  gint markid       = 0;
  gint nsegs        = 0;
  gint tfields      = 0;
  glong nrows       = 0;
  gboolean weighted = FALSE;
  gboolean missing  = FALSE;
  gboolean legacy   = FALSE;
  NcmMSet *mset     = NULL;
  gchar *rtype      = NULL;
  gchar *desc_rtype = NULL;
  GPtrArray *names;
  GPtrArray *symbs;
  guint i;

  fits_open_file (&fptr, storage_fits->filename, storage_fits->read_only ? READONLY : READWRITE, &status);
  _ncm_mset_catalog_storage_fits_check (status, storage_fits->filename, error);

  if (*error != NULL)
    return;

  storage_fits->fptr = fptr;
  names              = g_ptr_array_new_with_free_func (g_free);
  symbs              = g_ptr_array_new_with_free_func (g_free);

  mset = _ncm_mset_catalog_storage_fits_read_mset (storage_fits, error);

  if (*error == NULL)
  {
    fits_get_num_rows (fptr, &nrows, &status);
    _ncm_mset_catalog_storage_fits_check (status, "NAXIS2", error);
  }

  if (*error == NULL)
  {
    fits_get_num_cols (fptr, &tfields, &status);
    _ncm_mset_catalog_storage_fits_check (status, "TFIELDS", error);
  }

  if (*error == NULL)
    rtype = _ncm_mset_catalog_storage_fits_read_str (fptr, NCM_MSET_CATALOG_RTYPE_LABEL, error);

  if ((*error == NULL) && (rtype == NULL))
    rtype = g_strdup (NCM_MSET_CATALOG_RTYPE_UNDEFINED);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_int (fptr, NCM_MSET_CATALOG_NCHAINS_LABEL, &nchains, &missing, error);

  if ((*error == NULL) && (missing || (nchains <= 0)))
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "`%s' has no valid `%s' key", storage_fits->filename, NCM_MSET_CATALOG_NCHAINS_LABEL);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_int (fptr, NCM_MSET_CATALOG_NADDVAL_LABEL, &nadd_vals, &missing, error);

  if ((*error == NULL) && (missing || (nadd_vals < 0)))
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "`%s' has no valid `%s' key", storage_fits->filename, NCM_MSET_CATALOG_NADDVAL_LABEL);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_bool (fptr, NCM_MSET_CATALOG_WEIGHTED_LABEL, &weighted, NULL, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_int (fptr, NCM_MSET_CATALOG_M2LNP_ID_LABEL, &m2lnp_var, NULL, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_int (fptr, NCM_MSET_CATALOG_FIRST_ID_LABEL, &storage_fits->first_id, NULL, error);

  if (*error == NULL)
  {
    markid = storage_fits->first_id;
    _ncm_mset_catalog_storage_fits_read_int (fptr, NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL, &markid, NULL, error);
  }

  if ((*error == NULL) && ((guint) (ncm_mset_fparams_len (mset) + nadd_vals) != (guint) tfields))
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "`%s' has %d columns, its mset has %u free parameters and %d additional values",
                 storage_fits->filename, tfields, ncm_mset_fparams_len (mset), nadd_vals);

  /* Additional-value names and symbols: the first nadd_vals columns and ASYMB<n>. */
  for (i = 0; (i < (guint) nadd_vals) && (*error == NULL); i++)
  {
    gchar colname[FLEN_VALUE];
    gchar *asymbi = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_ASYMB_LABEL, i + 1);
    gchar *symbol;
    gchar *ttype = g_strdup_printf ("TTYPE%u", i + 1);

    fits_read_key_str (fptr, ttype, colname, NULL, &status);
    _ncm_mset_catalog_storage_fits_check (status, ttype, error);

    symbol = _ncm_mset_catalog_storage_fits_read_str (fptr, asymbi, error);

    if (*error == NULL)
    {
      g_ptr_array_add (names, g_strdup (colname));
      g_ptr_array_add (symbs, (symbol != NULL) ? symbol : g_strdup ("no-symbol"));
    }
    else
    {
      g_free (symbol);
    }

    g_free (ttype);
    g_free (asymbi);
  }

  if ((*error == NULL) && weighted)
  {
    if ((names->len == 0) || (strcmp (g_ptr_array_index (names, names->len - 1), NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME) != 0))
      g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                   "`%s' is weighted and its last additional value is not `%s'",
                   storage_fits->filename, NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME);
    else
      g_ptr_array_set_size (names, names->len - 1);

    g_ptr_array_set_size (symbs, names->len);
  }

  /* Sampler segments: NSEGS and SEGS<n>, SEGN<n>, SEGO<n>, or the translation of a file
   * without NSEGS. */
  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_int (fptr, NCM_MSET_CATALOG_STORAGE_FITS_NSEGS_LABEL, &nsegs, &legacy, error);

  if ((*error == NULL) && legacy)
  {
    desc_rtype = _ncm_mset_catalog_storage_fits_translate_legacy (storage_fits, rtype, nchains, nrows, error);
  }
  else if (*error == NULL)
  {
    desc_rtype = g_strdup (rtype);

    for (i = 0; (i < (guint) nsegs) && (*error == NULL); i++)
    {
      gchar *key_s   = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_STORAGE_FITS_SEG_START_LABEL, i);
      gchar *key_n   = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_STORAGE_FITS_SEG_SAMPLER_LABEL, i);
      gchar *key_o   = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_STORAGE_FITS_SEG_OPTIONS_LABEL, i);
      gint start     = 0;
      gchar *sampler = NULL;
      gchar *options = NULL;

      _ncm_mset_catalog_storage_fits_read_int (fptr, key_s, &start, &missing, error);

      if ((*error == NULL) && missing)
        g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                     "`%s' has %d segments and no `%s' key", storage_fits->filename, nsegs, key_s);

      if (*error == NULL)
        sampler = _ncm_mset_catalog_storage_fits_read_longstr (fptr, key_n, error);

      if (*error == NULL)
        options = _ncm_mset_catalog_storage_fits_read_longstr (fptr, key_o, error);

      if (*error == NULL)
        _ncm_mset_catalog_storage_fits_add_segment (storage_fits, start, (sampler != NULL) ? sampler : "", options);

      g_free (sampler);
      g_free (options);
      g_free (key_s);
      g_free (key_n);
      g_free (key_o);
    }
  }

  /* RNG record. */
  if (*error == NULL)
    storage_fits->rng_inis = _ncm_mset_catalog_storage_fits_read_longstr (fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, error);

  if (*error == NULL)
    storage_fits->rng_stat = _ncm_mset_catalog_storage_fits_read_longstr (fptr, NCM_MSET_CATALOG_RNG_STAT_LABEL, error);

  if (*error == NULL)
    storage_fits->rng_algo = _ncm_mset_catalog_storage_fits_read_str (fptr, NCM_MSET_CATALOG_RNG_ALGO_LABEL, error);

  if ((*error == NULL) && (storage_fits->rng_algo == NULL) && (storage_fits->rng_inis != NULL))
  {
    NcmRNG *rng = ncm_rng_new (NULL);

    storage_fits->rng_algo = g_strdup (ncm_rng_get_algo (rng));
    ncm_rng_free (rng);
  }

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_read_ulong (fptr, NCM_MSET_CATALOG_RNG_SEED_LABEL, &storage_fits->rng_seed, NULL, error);

  if ((*error == NULL) && (storage_fits->rng_algo != NULL) && (storage_fits->rng_inis == NULL))
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "`%s' has `%s' and no `%s' key", storage_fits->filename, NCM_MSET_CATALOG_RNG_ALGO_LABEL, NCM_MSET_CATALOG_RNG_INIS_LABEL);

  if (*error == NULL)
  {
    g_ptr_array_add (names, NULL);
    g_ptr_array_add (symbs, NULL);

    storage_fits->desc = ncm_mset_catalog_desc_new (mset, nchains, weighted,
                                                    (const gchar * const *) names->pdata,
                                                    (const gchar * const *) symbs->pdata,
                                                    m2lnp_var, desc_rtype, error);
  }

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_set_porder (storage_fits, storage_fits->desc, error);

  if (*error == NULL)
  {
    storage_fits->legacy       = legacy;
    storage_fits->tail_first   = (guint) nrows;
    storage_fits->n_pending    = 0;
    storage_fits->markovian_id = (markid > storage_fits->first_id) ? (guint) (markid - storage_fits->first_id) : 0;

    if (storage_fits->markovian_id > (guint) nrows)
      storage_fits->markovian_id = (guint) nrows;
  }

  ncm_mset_clear (&mset);
  g_free (rtype);
  g_free (desc_rtype);
  g_ptr_array_unref (names);
  g_ptr_array_unref (symbs);
}

/*
 * Opens the file at the first call that needs it. A storage whose file does not exist stays
 * uninitialized until store_desc creates it; read-only, the missing file is an error.
 */
static void
_ncm_mset_catalog_storage_fits_ensure_open (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  if (storage_fits->opened)
    return;

  if (storage_fits->filename == NULL)
  {
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "the storage has no file name");

    return;
  }

  if (g_file_test (storage_fits->filename, G_FILE_TEST_EXISTS))
  {
    _ncm_mset_catalog_storage_fits_open_existing (storage_fits, error);

    if (*error != NULL)
    {
      gint status = 0;

      if (storage_fits->fptr != NULL)
        fits_close_file (FPTR (storage_fits), &status);

      storage_fits->fptr = NULL;
      ncm_mset_catalog_desc_clear (&storage_fits->desc);
      g_array_set_size (storage_fits->segments, 0);
      g_clear_pointer (&storage_fits->rng_algo, g_free);
      g_clear_pointer (&storage_fits->rng_inis, g_free);
      g_clear_pointer (&storage_fits->rng_stat, g_free);

      return;
    }
  }
  else if (storage_fits->read_only)
  {
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
                 "the read-only file `%s' does not exist", storage_fits->filename);

    return;
  }

  storage_fits->opened = TRUE;
}

/* Writes the n_pending rows at the end of the tail to the file and drops from the tail the
 * rows beyond keep-last per chain. */
static void
_ncm_mset_catalog_storage_fits_write_pending (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  fitsfile *fptr    = FPTR (storage_fits);
  const guint ncols = storage_fits->porder->len;
  const guint n     = storage_fits->n_pending;
  const guint first = storage_fits->tail->len - n;
  const guint nrows = storage_fits->tail_first + first;
  gint status       = 0;
  gdouble *col;
  guint c, k;

  if (n == 0)
    return;

  g_array_set_size (storage_fits->col_ws, n);
  col = &g_array_index (storage_fits->col_ws, gdouble, 0);

  for (c = 0; (c < ncols) && (status == 0); c++)
  {
    for (k = 0; k < n; k++)
      col[k] = ncm_vector_get (g_ptr_array_index (storage_fits->tail, first + k), c);

    fits_write_col_dbl (fptr, g_array_index (storage_fits->porder, gint, c), nrows + 1, 1, n, col, &status);
  }

  _ncm_mset_catalog_storage_fits_check (status, "rows", error);

  if (*error == NULL)
  {
    const guint keep = storage_fits->keep_last * ncm_mset_catalog_desc_nchains (storage_fits->desc);

    storage_fits->n_pending = 0;

    if (storage_fits->tail->len > keep)
    {
      const guint drop = storage_fits->tail->len - keep;

      g_ptr_array_remove_range (storage_fits->tail, 0, drop);
      storage_fits->tail_first += drop;
    }
  }
}

/* Replaces the block by the n rows of the file starting at row first (0-based). */
static void
_ncm_mset_catalog_storage_fits_read_block (NcmMSetCatalogStorageFits *storage_fits, guint first, guint n, GError **error)
{
  fitsfile *fptr      = FPTR (storage_fits);
  const guint ncols   = storage_fits->porder->len;
  const gdouble dnull = 0.0;
  gint status         = 0;
  gdouble *col;
  guint c, k;

  g_ptr_array_set_size (storage_fits->block, 0);
  g_array_set_size (storage_fits->col_ws, n);
  col = &g_array_index (storage_fits->col_ws, gdouble, 0);

  for (k = 0; k < n; k++)
    g_ptr_array_add (storage_fits->block, ncm_vector_new (ncols));

  for (c = 0; (c < ncols) && (status == 0); c++)
  {
    fits_read_col_dbl (fptr, g_array_index (storage_fits->porder, gint, c), first + 1, 1, n, dnull, col, NULL, &status);

    for (k = 0; k < n; k++)
      ncm_vector_fast_set (g_ptr_array_index (storage_fits->block, k), c, col[k]);
  }

  _ncm_mset_catalog_storage_fits_check (status, "rows", error);

  if (*error == NULL)
    storage_fits->block_first = first;
  else
    g_ptr_array_set_size (storage_fits->block, 0);
}

/* Writes the segment record and the SAMPLER/SAMPOPT keys of the last segment. */
static void
_ncm_mset_catalog_storage_fits_write_segments (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  fitsfile *fptr = FPTR (storage_fits);
  const guint n  = storage_fits->segments->len;
  guint i;

  _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_STORAGE_FITS_NSEGS_LABEL, n, "Number of sampler segments.", error);

  for (i = 0; (i < n) && (*error == NULL); i++)
  {
    NcmMSetCatalogStorageFitsSegment *segment = &g_array_index (storage_fits->segments, NcmMSetCatalogStorageFitsSegment, i);
    gchar *key_s                              = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_STORAGE_FITS_SEG_START_LABEL, i);
    gchar *key_n                              = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_STORAGE_FITS_SEG_SAMPLER_LABEL, i);
    gchar *key_o                              = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_STORAGE_FITS_SEG_OPTIONS_LABEL, i);

    _ncm_mset_catalog_storage_fits_write_int (fptr, key_s, segment->start, "First row of the segment.", error);

    if (*error == NULL)
      _ncm_mset_catalog_storage_fits_write_longstr (fptr, key_n, segment->sampler, "Sampler of the segment.", error);

    if (*error == NULL)
      _ncm_mset_catalog_storage_fits_write_longstr (fptr, key_o, segment->options, "Sampler options of the segment.", error);

    g_free (key_s);
    g_free (key_n);
    g_free (key_o);
  }

  if ((*error == NULL) && (n > 0))
  {
    NcmMSetCatalogStorageFitsSegment *last = &g_array_index (storage_fits->segments, NcmMSetCatalogStorageFitsSegment, n - 1);

    _ncm_mset_catalog_storage_fits_write_longstr (fptr, NCM_MSET_CATALOG_SAMPLER_LABEL, last->sampler, "sampler; informative only.", error);

    if (*error == NULL)
      _ncm_mset_catalog_storage_fits_write_longstr (fptr, NCM_MSET_CATALOG_SAMPLER_OPTS_LABEL, last->options, "sampler options; informative only.", error);
  }
}

/* Writes the translated run type and the segment record of a file without NSEGS. */
static void
_ncm_mset_catalog_storage_fits_migrate_legacy (NcmMSetCatalogStorageFits *storage_fits, GError **error)
{
  if (!storage_fits->legacy)
    return;

  _ncm_mset_catalog_storage_fits_write_str (FPTR (storage_fits), NCM_MSET_CATALOG_RTYPE_LABEL,
                                            ncm_mset_catalog_desc_get_run_type (storage_fits->desc), "Run type string.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_segments (storage_fits, error);

  if (*error == NULL)
    storage_fits->legacy = FALSE;
}

/* Creates the file for desc: HDU0 with the mset, the table and the description keys. */
static void
_ncm_mset_catalog_storage_fits_create (NcmMSetCatalogStorageFits *storage_fits, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSet *mset          = ncm_mset_catalog_desc_peek_mset (desc);
  const guint nadd       = ncm_mset_catalog_desc_nadd_vals (desc);
  const guint ncols      = ncm_mset_catalog_desc_ncols (desc);
  const guint fparam_len = ncols - nadd;
  GPtrArray *ttype       = g_ptr_array_new ();
  GPtrArray *tform       = g_ptr_array_new ();
  fitsfile *fptr         = NULL;
  gint status            = 0;
  guint i;

  fits_create_file (&fptr, storage_fits->filename, &status);
  _ncm_mset_catalog_storage_fits_check (status, storage_fits->filename, error);

  if (*error != NULL)
  {
    g_ptr_array_unref (ttype);
    g_ptr_array_unref (tform);

    return;
  }

  storage_fits->fptr = fptr;

  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    NcmVarDict *vd    = ncm_var_dict_new ();
    GVariant *var;
    glong naxes[1];

    ncm_var_dict_set_object (vd, NCM_MSET_CATALOG_HDU0_MSET_KEY, ser, G_OBJECT (mset));
    var      = ncm_serialize_var_dict_to_variant (ser, vd);
    naxes[0] = (glong) g_variant_get_size (var);

    fits_create_img (fptr, BYTE_IMG, 1, naxes, &status);

    if (status == 0)
      fits_write_img (fptr, TBYTE, 1, naxes[0], (gpointer) g_variant_get_data (var), &status);

    if (status == 0)
      fits_write_key_str (fptr, NCM_MSET_CATALOG_MSET_FORMAT_LABEL, NCM_MSET_CATALOG_MSET_FORMAT_VARDICT,
                          "Format of the data stored in this HDU.", &status);

    _ncm_mset_catalog_storage_fits_check (status, "HDU0", error);

    g_variant_unref (var);
    ncm_var_dict_unref (vd);
    ncm_serialize_free (ser);
  }

  for (i = 0; i < ncols; i++)
  {
    g_ptr_array_add (ttype, (gchar *) ncm_mset_catalog_desc_col_full_name (desc, i));
    g_ptr_array_add (tform, "1D");
  }

  if (*error == NULL)
  {
    fits_create_tbl (fptr, BINARY_TBL, 0, ncols, (gchar **) ttype->pdata, (gchar **) tform->pdata,
                     NULL, NCM_MSET_CATALOG_EXTNAME, &status);
    _ncm_mset_catalog_storage_fits_check (status, NCM_MSET_CATALOG_EXTNAME, error);
  }

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_str (fptr, NCM_MSET_CATALOG_RTYPE_LABEL, ncm_mset_catalog_desc_get_run_type (desc), "Run type string.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_NCHAINS_LABEL, ncm_mset_catalog_desc_nchains (desc), "Number of chains.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_NADDVAL_LABEL, nadd, "Number of additional values.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_bool (fptr, NCM_MSET_CATALOG_WEIGHTED_LABEL, ncm_mset_catalog_desc_weighted (desc), "Whether the catalog is weighted.", error);

  for (i = 0; (i < nadd) && (*error == NULL); i++)
  {
    gchar *asymbi     = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_ASYMB_LABEL, i + 1);
    gchar *asymb_desc = g_strdup_printf ("Symbol for additional value %s[%u]", ncm_mset_catalog_desc_col_name (desc, i), i + 1);

    _ncm_mset_catalog_storage_fits_write_str (fptr, asymbi, ncm_mset_catalog_desc_col_symb (desc, i), asymb_desc, error);

    g_free (asymbi);
    g_free (asymb_desc);
  }

  for (i = 0; (i < fparam_len) && (*error == NULL); i++)
  {
    gchar *fsymbi     = g_strdup_printf ("%s%u", NCM_MSET_CATALOG_FSYMB_LABEL, i + 1);
    gchar *fsymb_desc = g_strdup_printf ("Symbol for parameter %s[%u]", ncm_mset_fparam_name (mset, i), i + 1);

    _ncm_mset_catalog_storage_fits_write_str (fptr, fsymbi, ncm_mset_fparam_symbol (mset, i), fsymb_desc, error);

    g_free (fsymbi);
    g_free (fsymb_desc);
  }

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_FIRST_ID_LABEL, 0, "Id of the first element.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_M2LNP_ID_LABEL, ncm_mset_catalog_desc_get_m2lnp_var (desc), "Id of the m2lnp variable.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL, 0, "First row id of the Markovian chain.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (fptr, NCM_MSET_CATALOG_STORAGE_FITS_NSEGS_LABEL, 0, "Number of sampler segments.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_flush (storage_fits, error);

  g_ptr_array_unref (ttype);
  g_ptr_array_unref (tform);
}

/* Writes the pending rows of a writable storage and closes the file. */
static void
_ncm_mset_catalog_storage_fits_close (NcmMSetCatalogStorageFits *storage_fits)
{
  if (storage_fits->fptr != NULL)
  {
    GError *error = NULL;
    gint status   = 0;

    if (!storage_fits->read_only)
      _ncm_mset_catalog_storage_fits_write_pending (storage_fits, &error);

    if ((error == NULL) && !storage_fits->read_only)
      _ncm_mset_catalog_storage_fits_flush (storage_fits, &error);

    if (error != NULL)
    {
      g_warning ("ncm_mset_catalog_storage_fits: closing `%s': %s", storage_fits->filename, error->message);
      g_clear_error (&error);
    }

    fits_close_file (FPTR (storage_fits), &status);
    storage_fits->fptr = NULL;
  }
}

/*
 * The virtual functions. Each opens the file first; a getter on a storage without a file
 * answers the empty state.
 */

static NcmMSetCatalogDesc *
_ncm_mset_catalog_storage_fits_peek_desc (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  return storage_fits->desc;
}

static void
_ncm_mset_catalog_storage_fits_store_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  g_assert_null (storage_fits->desc);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_require_writable (storage_fits, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_create (storage_fits, desc, error);

  if (*error == NULL)
  {
    storage_fits->desc = ncm_mset_catalog_desc_ref (desc);
    _ncm_mset_catalog_storage_fits_set_porder (storage_fits, desc, error);
  }
}

/* Binds desc: the column map follows its order and the rows in memory are permuted. */
static void
_ncm_mset_catalog_storage_fits_bind_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);
  NcmMSet *mset_stored                    = ncm_mset_catalog_desc_peek_mset (storage_fits->desc);
  NcmMSet *mset_bound                     = ncm_mset_catalog_desc_peek_mset (desc);
  const guint nadd                        = ncm_mset_catalog_desc_nadd_vals (desc);
  const guint nfree                       = ncm_mset_fparams_len (mset_bound);
  GArray *perm                            = g_array_sized_new (FALSE, FALSE, sizeof (guint), nfree);
  guint i, j;

  for (i = 0; i < nfree; i++)
  {
    const gchar *name = ncm_mset_fparam_full_name (mset_bound, i);

    for (j = 0; j < nfree; j++)
    {
      if (strcmp (name, ncm_mset_fparam_full_name (mset_stored, j)) == 0)
        break;
    }

    g_assert_cmpuint (j, <, nfree);
    g_array_append_val (perm, j);
  }

  for (i = 0; i < storage_fits->tail->len; i++)
  {
    NcmVector *row     = g_ptr_array_index (storage_fits->tail, i);
    NcmVector *new_row = ncm_vector_dup (row);

    for (j = 0; j < nfree; j++)
      ncm_vector_set (new_row, nadd + j, ncm_vector_get (row, nadd + g_array_index (perm, guint, j)));

    ncm_vector_free (row);
    g_ptr_array_index (storage_fits->tail, i) = new_row;
  }

  g_array_unref (perm);
  g_ptr_array_set_size (storage_fits->block, 0);

  _ncm_mset_catalog_storage_fits_set_porder (storage_fits, desc, error);

  if (*error == NULL)
  {
    ncm_mset_catalog_desc_clear (&storage_fits->desc);
    storage_fits->desc = ncm_mset_catalog_desc_ref (desc);
  }
}

static guint
_ncm_mset_catalog_storage_fits_nrows (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  return storage_fits->tail_first + storage_fits->tail->len;
}

static void
_ncm_mset_catalog_storage_fits_append (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_require_writable (storage_fits, error);

  if (*error == NULL)
  {
    g_ptr_array_add (storage_fits->tail, ncm_vector_dup (row));
    storage_fits->n_pending++;
  }
}

static NcmVector *
_ncm_mset_catalog_storage_fits_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  if (*error != NULL)
    return NULL;

  if (i >= storage_fits->tail_first)
    return ncm_vector_ref (g_ptr_array_index (storage_fits->tail, i - storage_fits->tail_first));

  if ((storage_fits->block->len == 0) || (i < storage_fits->block_first) || (i >= storage_fits->block_first + storage_fits->block->len))
  {
    const guint first = i - (i % storage_fits->read_block);

    _ncm_mset_catalog_storage_fits_read_block (storage_fits, first, MIN (storage_fits->read_block, storage_fits->tail_first - first), error);

    if (*error != NULL)
      return NULL;
  }

  return ncm_vector_ref (g_ptr_array_index (storage_fits->block, i - storage_fits->block_first));
}

static void
_ncm_mset_catalog_storage_fits_sync (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  if (storage_fits->read_only || (storage_fits->fptr == NULL))
    return;

  _ncm_mset_catalog_storage_fits_write_pending (storage_fits, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_flush (storage_fits, error);
}

static void
_ncm_mset_catalog_storage_fits_get_rng_init (NcmMSetCatalogStorage *storage, gchar **algo, gulong *seed, gchar **state, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  if (algo != NULL)
    *algo = g_strdup (storage_fits->rng_algo);

  if (seed != NULL)
    *seed = storage_fits->rng_seed;

  if (state != NULL)
    *state = g_strdup (storage_fits->rng_inis);
}

static void
_ncm_mset_catalog_storage_fits_set_rng_init (NcmMSetCatalogStorage *storage, const gchar *algo, gulong seed, const gchar *state, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);
  fitsfile *fptr                          = FPTR (storage_fits);

  _ncm_mset_catalog_storage_fits_require_writable (storage_fits, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_str (fptr, NCM_MSET_CATALOG_RNG_ALGO_LABEL, algo, "RNG Algorithm name.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_ulong (fptr, NCM_MSET_CATALOG_RNG_SEED_LABEL, seed, "RNG Algorithm seed.", error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_longstr (fptr, NCM_MSET_CATALOG_RNG_INIS_LABEL, state, NULL, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_flush (storage_fits, error);

  if (*error == NULL)
  {
    g_free (storage_fits->rng_algo);
    g_free (storage_fits->rng_inis);
    storage_fits->rng_algo = g_strdup (algo);
    storage_fits->rng_seed = seed;
    storage_fits->rng_inis = g_strdup (state);
  }
}

static gchar *
_ncm_mset_catalog_storage_fits_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  return g_strdup (storage_fits->rng_stat);
}

static void
_ncm_mset_catalog_storage_fits_set_rng_stat (NcmMSetCatalogStorage *storage, const gchar *state, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_require_writable (storage_fits, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_longstr (FPTR (storage_fits), NCM_MSET_CATALOG_RNG_STAT_LABEL, state, NULL, error);

  if (*error == NULL)
  {
    g_free (storage_fits->rng_stat);
    storage_fits->rng_stat = g_strdup (state);
  }
}

static guint
_ncm_mset_catalog_storage_fits_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  return storage_fits->markovian_id;
}

static void
_ncm_mset_catalog_storage_fits_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_require_writable (storage_fits, error);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_write_int (FPTR (storage_fits), NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL,
                                              (gint) markovian_id + storage_fits->first_id, "First row id of the Markovian chain.", error);

  if (*error == NULL)
    storage_fits->markovian_id = markovian_id;
}

static guint
_ncm_mset_catalog_storage_fits_n_samplers (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  return storage_fits->segments->len;
}

static void
_ncm_mset_catalog_storage_fits_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);
  NcmMSetCatalogStorageFitsSegment *segment;

  _ncm_mset_catalog_storage_fits_ensure_open (storage_fits, error);

  if (*error != NULL)
    return;

  if (i >= storage_fits->segments->len)
  {
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                 "sampler segment %u, the storage has %u", i, storage_fits->segments->len);

    return;
  }

  segment = &g_array_index (storage_fits->segments, NcmMSetCatalogStorageFitsSegment, i);

  if (start != NULL)
    *start = segment->start;

  if (sampler != NULL)
    *sampler = g_strdup (segment->sampler);

  if (options != NULL)
    *options = g_strdup (segment->options);
}

static void
_ncm_mset_catalog_storage_fits_set_sampler (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = NCM_MSET_CATALOG_STORAGE_FITS (storage);

  _ncm_mset_catalog_storage_fits_require_writable (storage_fits, error);

  if ((*error == NULL) && (i > storage_fits->segments->len))
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                 "sampler segment %u, the storage has %u", i, storage_fits->segments->len);

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_migrate_legacy (storage_fits, error);

  if (*error == NULL)
  {
    if (i == storage_fits->segments->len)
    {
      _ncm_mset_catalog_storage_fits_add_segment (storage_fits, start, sampler, options);
    }
    else
    {
      NcmMSetCatalogStorageFitsSegment *segment = &g_array_index (storage_fits->segments, NcmMSetCatalogStorageFitsSegment, i);

      _ncm_mset_catalog_storage_fits_segment_clear (segment);
      segment->start   = start;
      segment->sampler = g_strdup (sampler);
      segment->options = g_strdup ((options != NULL) ? options : "");
    }

    _ncm_mset_catalog_storage_fits_write_segments (storage_fits, error);
  }

  if (*error == NULL)
    _ncm_mset_catalog_storage_fits_flush (storage_fits, error);
}

#else /* HAVE_CFITSIO */

static void
_ncm_mset_catalog_storage_fits_no_cfitsio (GError **error)
{
  g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NO_CFITSIO,
               "NumCosmo was built without CFITSIO");
}

static void
_ncm_mset_catalog_storage_fits_close (NcmMSetCatalogStorageFits *storage_fits)
{
}

static NcmMSetCatalogDesc *
_ncm_mset_catalog_storage_fits_peek_desc (NcmMSetCatalogStorage *storage, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);

  return NULL;
}

static void
_ncm_mset_catalog_storage_fits_store_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static void
_ncm_mset_catalog_storage_fits_bind_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static guint
_ncm_mset_catalog_storage_fits_nrows (NcmMSetCatalogStorage *storage, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);

  return 0;
}

static void
_ncm_mset_catalog_storage_fits_append (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static NcmVector *
_ncm_mset_catalog_storage_fits_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);

  return NULL;
}

static void
_ncm_mset_catalog_storage_fits_sync (NcmMSetCatalogStorage *storage, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static void
_ncm_mset_catalog_storage_fits_get_rng_init (NcmMSetCatalogStorage *storage, gchar **algo, gulong *seed, gchar **state, GError **error)
{
  if (algo != NULL)
    *algo = NULL;

  if (seed != NULL)
    *seed = 0;

  if (state != NULL)
    *state = NULL;

  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static void
_ncm_mset_catalog_storage_fits_set_rng_init (NcmMSetCatalogStorage *storage, const gchar *algo, gulong seed, const gchar *state, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static gchar *
_ncm_mset_catalog_storage_fits_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);

  return NULL;
}

static void
_ncm_mset_catalog_storage_fits_set_rng_stat (NcmMSetCatalogStorage *storage, const gchar *state, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static guint
_ncm_mset_catalog_storage_fits_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);

  return 0;
}

static void
_ncm_mset_catalog_storage_fits_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static guint
_ncm_mset_catalog_storage_fits_n_samplers (NcmMSetCatalogStorage *storage, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);

  return 0;
}

static void
_ncm_mset_catalog_storage_fits_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

static void
_ncm_mset_catalog_storage_fits_set_sampler (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options, GError **error)
{
  _ncm_mset_catalog_storage_fits_no_cfitsio (error);
}

#endif /* HAVE_CFITSIO */

/**
 * ncm_mset_catalog_storage_fits_new:
 * @filename: the FITS file
 * @read_only: whether the file is opened for reading only
 * @error: a #GError
 *
 * Creates a storage on @filename and opens the file. A file that does not exist leaves the
 * storage uninitialized, to be created by ncm_mset_catalog_storage_set_desc(); read-only it
 * is an error. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND when the file cannot be
 * read and with %NCM_MSET_CATALOG_STORAGE_ERROR_NO_CFITSIO without CFITSIO.
 *
 * Returns: (transfer full) (nullable): a new #NcmMSetCatalogStorageFits, %NULL on error
 */
NcmMSetCatalogStorageFits *
ncm_mset_catalog_storage_fits_new (const gchar *filename, gboolean read_only, GError **error)
{
  NcmMSetCatalogStorageFits *storage_fits = g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_FITS,
                                                          "filename", filename,
                                                          "read-only", read_only,
                                                          NULL);
  GError *local_error = NULL;

  ncm_mset_catalog_storage_peek_desc (NCM_MSET_CATALOG_STORAGE (storage_fits), &local_error);

  if (local_error != NULL)
  {
    ncm_mset_catalog_storage_fits_clear (&storage_fits);
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_fits_new");

    return NULL;
  }

  return storage_fits;
}

/**
 * ncm_mset_catalog_storage_fits_ref:
 * @storage_fits: a #NcmMSetCatalogStorageFits
 *
 * Increases the reference count of @storage_fits by one.
 *
 * Returns: (transfer full): @storage_fits
 */
NcmMSetCatalogStorageFits *
ncm_mset_catalog_storage_fits_ref (NcmMSetCatalogStorageFits *storage_fits)
{
  return g_object_ref (storage_fits);
}

/**
 * ncm_mset_catalog_storage_fits_free:
 * @storage_fits: a #NcmMSetCatalogStorageFits
 *
 * Decreases the reference count of @storage_fits by one.
 */
void
ncm_mset_catalog_storage_fits_free (NcmMSetCatalogStorageFits *storage_fits)
{
  g_object_unref (storage_fits);
}

/**
 * ncm_mset_catalog_storage_fits_clear:
 * @storage_fits: a #NcmMSetCatalogStorageFits
 *
 * Decreases the reference count of *@storage_fits by one and sets *@storage_fits to %NULL.
 */
void
ncm_mset_catalog_storage_fits_clear (NcmMSetCatalogStorageFits **storage_fits)
{
  g_clear_object (storage_fits);
}

/**
 * ncm_mset_catalog_storage_fits_peek_filename:
 * @storage_fits: a #NcmMSetCatalogStorageFits
 *
 * Returns: (transfer none): the file name of @storage_fits
 */
const gchar *
ncm_mset_catalog_storage_fits_peek_filename (NcmMSetCatalogStorageFits *storage_fits)
{
  return storage_fits->filename;
}

/**
 * ncm_mset_catalog_storage_fits_get_read_only:
 * @storage_fits: a #NcmMSetCatalogStorageFits
 *
 * Returns: whether @storage_fits is read-only
 */
gboolean
ncm_mset_catalog_storage_fits_get_read_only (NcmMSetCatalogStorageFits *storage_fits)
{
  return storage_fits->read_only;
}


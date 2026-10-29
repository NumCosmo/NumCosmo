/***************************************************************************
 *            ncm_mset_catalog_storage.c
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_storage.c
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
 * NcmMSetCatalogStorage:
 *
 * Abstract storage of the rows, the #NcmMSetCatalogDesc and the run metadata of an
 * #NcmMSetCatalog.
 *
 * A storage without a description is uninitialized. ncm_mset_catalog_storage_set_desc()
 * stores a description on an uninitialized storage and binds a compatible one
 * (ncm_mset_catalog_desc_check_compatible()) on an initialized storage; rows are then read
 * and written in the order of the bound description.
 *
 * Rows and run metadata are written only inside a recording, opened by
 * ncm_mset_catalog_storage_start_recording() and closed by
 * ncm_mset_catalog_storage_end_recording(). One recording is open at a time. A recording
 * opens a sampler segment at the current number of rows and takes the run's #NcmRNG: on an
 * empty storage its state becomes the initial RNG state, on a non-empty one the stored RNG
 * state is restored into it. The stored RNG state is read from the #NcmRNG at each sync,
 * which happens only right after a row is added (as #NcmMSetCatalogStorage:smode requires)
 * and when the recording ends, so it is the state after the last written row; the caller
 * ends the recording before drawing from the #NcmRNG again.
 *
 * The run metadata are the initial and the current RNG state, the Markovian id (the first
 * row of the Markovian chain, never decreasing and never beyond the rows), the initial
 * sampler and the function array (each written once), and the sampler segments. The
 * backend is their only holder.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_catalog_storage.h"
#include "ncm/fit/ncm_mset_catalog.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-mset-catalog-storage-error, ncm_mset_catalog_storage_error)
/* *INDENT-ON* */

typedef struct _NcmMSetCatalogStoragePrivate
{
  gboolean recording;
  NcmRNG *rng;
  NcmMSetCatalogSync smode;
  gdouble sync_interval;
  GTimer *sync_timer;
} NcmMSetCatalogStoragePrivate;

enum
{
  PROP_0,
  PROP_SYNC_MODE,
  PROP_SYNC_INTERVAL,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmMSetCatalogStorage, ncm_mset_catalog_storage, G_TYPE_OBJECT)

static void
ncm_mset_catalog_storage_init (NcmMSetCatalogStorage *storage)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  self->recording     = FALSE;
  self->rng           = NULL;
  self->smode         = NCM_MSET_CATALOG_SYNC_DISABLE;
  self->sync_interval = 0.0;
  self->sync_timer    = g_timer_new ();
}

static void
_ncm_mset_catalog_storage_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogStorage *storage            = NCM_MSET_CATALOG_STORAGE (object);
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  g_return_if_fail (NCM_IS_MSET_CATALOG_STORAGE (object));

  switch (prop_id)
  {
    case PROP_SYNC_MODE:
      self->smode = g_value_get_enum (value);
      break;
    case PROP_SYNC_INTERVAL:
      self->sync_interval = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_storage_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogStorage *storage            = NCM_MSET_CATALOG_STORAGE (object);
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  g_return_if_fail (NCM_IS_MSET_CATALOG_STORAGE (object));

  switch (prop_id)
  {
    case PROP_SYNC_MODE:
      g_value_set_enum (value, self->smode);
      break;
    case PROP_SYNC_INTERVAL:
      g_value_set_double (value, self->sync_interval);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_storage_dispose (GObject *object)
{
  NcmMSetCatalogStorage *storage            = NCM_MSET_CATALOG_STORAGE (object);
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  /* The writer ends the recording; a storage released with one open keeps what the last
   * sync wrote. */
  self->recording = FALSE;
  ncm_rng_clear (&self->rng);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_parent_class)->dispose (object);
}

static void
_ncm_mset_catalog_storage_finalize (GObject *object)
{
  NcmMSetCatalogStorage *storage            = NCM_MSET_CATALOG_STORAGE (object);
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  g_clear_pointer (&self->sync_timer, g_timer_destroy);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_parent_class)->finalize (object);
}

static void
ncm_mset_catalog_storage_class_init (NcmMSetCatalogStorageClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_mset_catalog_storage_set_property;
  object_class->get_property = &_ncm_mset_catalog_storage_get_property;
  object_class->dispose      = &_ncm_mset_catalog_storage_dispose;
  object_class->finalize     = &_ncm_mset_catalog_storage_finalize;

  /**
   * NcmMSetCatalogStorage:smode:
   *
   * When ncm_mset_catalog_storage_add_row() syncs: after every row, after the
   * #NcmMSetCatalogStorage:sync-interval has elapsed since the last sync, or never (only
   * ncm_mset_catalog_storage_end_recording() syncs). A backend with nothing to write (the
   * memory storage) never syncs inside ncm_mset_catalog_storage_add_row().
   */
  g_object_class_install_property (object_class,
                                   PROP_SYNC_MODE,
                                   g_param_spec_enum ("smode",
                                                      NULL,
                                                      "Storage sync mode",
                                                      NCM_TYPE_MSET_CATALOG_SYNC, NCM_MSET_CATALOG_SYNC_AUTO,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorage:sync-interval:
   *
   * Minimum time between syncs, in seconds, for %NCM_MSET_CATALOG_SYNC_TIMED.
   */
  g_object_class_install_property (object_class,
                                   PROP_SYNC_INTERVAL,
                                   g_param_spec_double ("sync-interval",
                                                        NULL,
                                                        "Data sync interval",
                                                        0.0, 1.0e3, 10.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->peek_desc  = NULL;
  klass->store_desc = NULL;
  klass->bind_desc  = NULL;
  klass->nrows      = NULL;
  klass->append     = NULL;
  klass->get_row    = NULL;
  klass->sync       = NULL;

  klass->get_rng_inis        = NULL;
  klass->set_rng_inis        = NULL;
  klass->get_rng_stat        = NULL;
  klass->set_rng_stat        = NULL;
  klass->get_markovian_id    = NULL;
  klass->set_markovian_id    = NULL;
  klass->get_initial_sampler = NULL;
  klass->set_initial_sampler = NULL;
  klass->get_functions_array = NULL;
  klass->set_functions_array = NULL;
  klass->n_samplers          = NULL;
  klass->get_sampler         = NULL;
  klass->set_sampler         = NULL;
}

/**
 * ncm_mset_catalog_storage_ref:
 * @storage: a #NcmMSetCatalogStorage
 *
 * Increases the reference count of @storage by one.
 *
 * Returns: (transfer full): @storage
 */
NcmMSetCatalogStorage *
ncm_mset_catalog_storage_ref (NcmMSetCatalogStorage *storage)
{
  return g_object_ref (storage);
}

/**
 * ncm_mset_catalog_storage_free:
 * @storage: a #NcmMSetCatalogStorage
 *
 * Decreases the reference count of @storage by one.
 */
void
ncm_mset_catalog_storage_free (NcmMSetCatalogStorage *storage)
{
  g_object_unref (storage);
}

/**
 * ncm_mset_catalog_storage_clear:
 * @storage: a #NcmMSetCatalogStorage
 *
 * If *@storage is not %NULL, decreases its reference count and sets *@storage to %NULL.
 */
void
ncm_mset_catalog_storage_clear (NcmMSetCatalogStorage **storage)
{
  g_clear_object (storage);
}

/**
 * ncm_mset_catalog_storage_peek_desc:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns the description of @storage, %NULL when @storage is uninitialized. Sets @error
 * when the backend fails.
 *
 * Returns: (transfer none) (nullable): the #NcmMSetCatalogDesc of @storage
 */
NcmMSetCatalogDesc *
ncm_mset_catalog_storage_peek_desc (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->peek_desc (storage, error);
}

/* Returns the description of storage, or sets error (UNINITIALIZED or the backend's). */
static NcmMSetCatalogDesc *
_ncm_mset_catalog_storage_require_desc (NcmMSetCatalogStorage *storage, const gchar *func, GError **error)
{
  GError *local_error      = NULL;
  NcmMSetCatalogDesc *desc = ncm_mset_catalog_storage_peek_desc (storage, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "%s", func);

    return NULL;
  }

  if (desc == NULL)
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED,
                                "%s: the storage has no description.", func);

  return desc;
}

/**
 * ncm_mset_catalog_storage_set_desc:
 * @storage: a #NcmMSetCatalogStorage
 * @desc: a #NcmMSetCatalogDesc
 * @error: a #GError
 *
 * Stores @desc on an uninitialized @storage. On an initialized @storage it binds @desc when
 * it is compatible with the stored description (ncm_mset_catalog_desc_check_compatible()),
 * so rows are read and written in the order of @desc, and otherwise sets @error with the
 * incompatibility. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING while a recording is
 * open.
 *
 * Returns: whether @desc is now the description of @storage
 */
gboolean
ncm_mset_catalog_storage_set_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;
  NcmMSetCatalogDesc *stored;

  if (self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                                "ncm_mset_catalog_storage_set_desc: a recording is open.");

    return FALSE;
  }

  stored = ncm_mset_catalog_storage_peek_desc (storage, &local_error);

  if (local_error == NULL)
  {
    if (stored == NULL)
      klass->store_desc (storage, desc, &local_error);
    else if ((stored != desc) && ncm_mset_catalog_desc_check_compatible (stored, desc, &local_error))
      klass->bind_desc (storage, desc, &local_error);
  }

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_desc");

    return FALSE;
  }

  return TRUE;
}


static void _ncm_mset_catalog_storage_open_segment (NcmMSetCatalogStorage *storage, guint nrows, const gchar *sampler, const gchar *options, GError **error);
static void _ncm_mset_catalog_storage_take_rng (NcmMSetCatalogStorage *storage, guint nrows, NcmRNG *rng, GError **error);

/**
 * ncm_mset_catalog_storage_start_recording:
 * @storage: a #NcmMSetCatalogStorage
 * @sampler: the sampler description
 * @options: (nullable): the sampler options
 * @rng: (nullable): the #NcmRNG of the run
 * @error: a #GError
 *
 * Opens a recording, the only state in which rows and run metadata are written. It opens a
 * sampler segment starting at the current number of rows; a last segment that starts at the
 * same row, and so has no rows, is replaced. When @rng is given and @storage is empty, the
 * state of @rng becomes the initial and the current RNG state; when @storage has rows, the
 * stored RNG state, if any, is set on @rng. @storage keeps a reference to @rng until the
 * recording ends.
 *
 * Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED without a description and with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING when a recording is already open.
 *
 * Returns: whether the recording was opened
 */
gboolean
ncm_mset_catalog_storage_start_recording (NcmMSetCatalogStorage *storage, const gchar *sampler, const gchar *options, NcmRNG *rng, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  GError *local_error                       = NULL;
  guint nrows                               = 0;

  g_assert (sampler != NULL);

  if (self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                                "ncm_mset_catalog_storage_start_recording: a recording is already open.");

    return FALSE;
  }

  if (_ncm_mset_catalog_storage_require_desc (storage, "ncm_mset_catalog_storage_start_recording", error) == NULL)
    return FALSE;

  nrows = ncm_mset_catalog_storage_nrows (storage, &local_error);

  if (local_error == NULL)
    _ncm_mset_catalog_storage_open_segment (storage, nrows, sampler, options, &local_error);

  if ((local_error == NULL) && (rng != NULL))
    _ncm_mset_catalog_storage_take_rng (storage, nrows, rng, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_start_recording");

    return FALSE;
  }

  self->recording = TRUE;

  if (rng != NULL)
    self->rng = ncm_rng_ref (rng);

  g_timer_start (self->sync_timer);

  return TRUE;
}

/* Writes the segment (nrows, sampler, options): over the last segment when it starts at
 * nrows, after it otherwise. */
static void
_ncm_mset_catalog_storage_open_segment (NcmMSetCatalogStorage *storage, guint nrows, const gchar *sampler, const gchar *options, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  const guint n                     = klass->n_samplers (storage, &local_error);
  guint i                           = n;

  if ((local_error == NULL) && (n > 0))
  {
    guint last_start = 0;

    klass->get_sampler (storage, n - 1, &last_start, NULL, NULL, &local_error);

    if (last_start == nrows)
      i = n - 1;
  }

  if (local_error == NULL)
    klass->set_sampler (storage, i, nrows, sampler, (options != NULL) ? options : "", &local_error);

  ncm_util_forward_or_call_error (error, local_error, "_ncm_mset_catalog_storage_open_segment");
}

/* On an empty storage writes the state of rng as the initial and the current RNG state; on
 * a storage with rows sets the stored RNG state, if any, on rng. */
static void
_ncm_mset_catalog_storage_take_rng (NcmMSetCatalogStorage *storage, guint nrows, NcmRNG *rng, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;

  if (nrows == 0)
  {
    gchar *state = ncm_rng_get_state (rng);

    klass->set_rng_inis (storage, state, &local_error);

    if (local_error == NULL)
      klass->set_rng_stat (storage, state, &local_error);

    g_free (state);
  }
  else
  {
    gchar *state = klass->get_rng_stat (storage, &local_error);

    if (state != NULL)
      ncm_rng_set_state (rng, state);

    g_free (state);
  }

  ncm_util_forward_or_call_error (error, local_error, "_ncm_mset_catalog_storage_take_rng");
}

/* Writes the state of the recording's RNG, then lets the backend write what it holds. */
static void
_ncm_mset_catalog_storage_sync (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;

  if (self->rng != NULL)
  {
    gchar *state = ncm_rng_get_state (self->rng);

    klass->set_rng_stat (storage, state, &local_error);
    g_free (state);
  }

  if ((local_error == NULL) && (klass->sync != NULL))
    klass->sync (storage, &local_error);

  ncm_util_forward_or_call_error (error, local_error, "_ncm_mset_catalog_storage_sync");
}

/* Returns FALSE and sets error outside a recording. */
static gboolean
_ncm_mset_catalog_storage_require_recording (NcmMSetCatalogStorage *storage, const gchar *func, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  if (!self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING,
                                "%s: no recording is open.", func);

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_storage_add_row:
 * @storage: a #NcmMSetCatalogStorage
 * @row: a #NcmVector
 * @error: a #GError
 *
 * Appends a copy of @row, in the order of the bound description, then syncs as
 * #NcmMSetCatalogStorage:smode requires. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a recording and with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when the length of @row is not the number of columns.
 *
 * Returns: whether the row was added
 */
gboolean
ncm_mset_catalog_storage_add_row (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;
  NcmMSetCatalogDesc *desc;

  if (!_ncm_mset_catalog_storage_require_recording (storage, "ncm_mset_catalog_storage_add_row", error))
    return FALSE;

  desc = _ncm_mset_catalog_storage_require_desc (storage, "ncm_mset_catalog_storage_add_row", error);

  if (desc == NULL)
    return FALSE;

  if (ncm_vector_len (row) != ncm_mset_catalog_desc_ncols (desc))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_add_row: row of length %u, the description has %u columns.",
                                ncm_vector_len (row), ncm_mset_catalog_desc_ncols (desc));

    return FALSE;
  }

  klass->append (storage, row, &local_error);

  if ((local_error == NULL) && (klass->sync != NULL))
  {
    switch (self->smode)
    {
      case NCM_MSET_CATALOG_SYNC_AUTO:
        _ncm_mset_catalog_storage_sync (storage, &local_error);
        break;
      case NCM_MSET_CATALOG_SYNC_TIMED:

        if (g_timer_elapsed (self->sync_timer, NULL) > self->sync_interval)
        {
          g_timer_start (self->sync_timer);
          _ncm_mset_catalog_storage_sync (storage, &local_error);
        }

        break;
      default:
        break;
    }
  }

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_add_row");

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_storage_end_recording:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Syncs, writing the state of the recording's #NcmRNG, and closes the recording. The
 * recording is closed even when the sync fails. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a recording.
 *
 * Returns: whether the sync succeeded
 */
gboolean
ncm_mset_catalog_storage_end_recording (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  GError *local_error                       = NULL;

  if (!_ncm_mset_catalog_storage_require_recording (storage, "ncm_mset_catalog_storage_end_recording", error))
    return FALSE;

  _ncm_mset_catalog_storage_sync (storage, &local_error);

  self->recording = FALSE;
  ncm_rng_clear (&self->rng);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_end_recording");

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_storage_is_recording:
 * @storage: a #NcmMSetCatalogStorage
 *
 * Returns: whether a recording is open
 */
gboolean
ncm_mset_catalog_storage_is_recording (NcmMSetCatalogStorage *storage)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  return self->recording;
}

/**
 * ncm_mset_catalog_storage_nrows:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: the number of rows, zero for an uninitialized @storage
 */
guint
ncm_mset_catalog_storage_nrows (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->nrows (storage, error);
}

/**
 * ncm_mset_catalog_storage_get_row:
 * @storage: a #NcmMSetCatalogStorage
 * @i: row index
 * @error: a #GError
 *
 * Returns row @i in the order of the bound description. The returned vector must not be
 * modified. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when @i is not below the number
 * of rows.
 *
 * Returns: (transfer full) (nullable): row @i, or %NULL on error
 */
NcmVector *
ncm_mset_catalog_storage_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error)
{
  GError *local_error = NULL;
  const guint nrows   = ncm_mset_catalog_storage_nrows (storage, &local_error);
  NcmVector *row      = NULL;

  if ((local_error == NULL) && (i >= nrows))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_get_row: row %u, the storage has %u rows.", i, nrows);

    return NULL;
  }

  if (local_error == NULL)
    row = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_row (storage, i, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_get_row");
    ncm_vector_clear (&row);
  }

  return row;
}

/**
 * ncm_mset_catalog_storage_set_markovian_id:
 * @storage: a #NcmMSetCatalogStorage
 * @markovian_id: id of the first row of the Markovian chain
 * @error: a #GError
 *
 * Sets the Markovian id. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a
 * recording and with %NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID when @markovian_id is
 * below the stored one or above the number of rows.
 *
 * Returns: whether the Markovian id was set
 */
gboolean
ncm_mset_catalog_storage_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  guint current                     = 0;
  guint nrows                       = 0;

  if (!_ncm_mset_catalog_storage_require_recording (storage, "ncm_mset_catalog_storage_set_markovian_id", error))
    return FALSE;

  current = klass->get_markovian_id (storage, &local_error);

  if (local_error == NULL)
    nrows = ncm_mset_catalog_storage_nrows (storage, &local_error);

  if ((local_error == NULL) && ((markovian_id < current) || (markovian_id > nrows)))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID,
                                "ncm_mset_catalog_storage_set_markovian_id: id %u is below the stored id %u or above the %u rows.",
                                markovian_id, current, nrows);

    return FALSE;
  }

  if (local_error == NULL)
    klass->set_markovian_id (storage, markovian_id, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_markovian_id");

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_storage_set_initial_sampler:
 * @storage: a #NcmMSetCatalogStorage
 * @sampler: the description of the sampler of the first ensemble
 * @error: a #GError
 *
 * Sets the initial sampler, written once: setting the stored value again succeeds, setting
 * another fails with %NCM_MSET_CATALOG_STORAGE_ERROR_WRITE_ONCE. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a recording.
 *
 * Returns: whether @sampler is the stored initial sampler
 */
gboolean
ncm_mset_catalog_storage_set_initial_sampler (NcmMSetCatalogStorage *storage, const gchar *sampler, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  gchar *current;

  if (!_ncm_mset_catalog_storage_require_recording (storage, "ncm_mset_catalog_storage_set_initial_sampler", error))
    return FALSE;

  current = klass->get_initial_sampler (storage, &local_error);

  if (local_error == NULL)
  {
    if (current == NULL)
      klass->set_initial_sampler (storage, sampler, &local_error);
    else if (strcmp (current, sampler) != 0)
      g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_WRITE_ONCE,
                   "the initial sampler is already `%s'", current);
  }

  g_free (current);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_initial_sampler");

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_storage_set_functions_array:
 * @storage: a #NcmMSetCatalogStorage
 * @functions: a #NcmObjArray of #NcmMSetFunc
 * @error: a #GError
 *
 * Sets the function array whose values are additional values, written once and compared by
 * serialization: setting an equal array again succeeds, setting another fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_WRITE_ONCE. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a recording.
 *
 * Returns: whether @functions equals the stored function array
 */
gboolean
ncm_mset_catalog_storage_set_functions_array (NcmMSetCatalogStorage *storage, NcmObjArray *functions, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  NcmObjArray *current;

  if (!_ncm_mset_catalog_storage_require_recording (storage, "ncm_mset_catalog_storage_set_functions_array", error))
    return FALSE;

  current = klass->get_functions_array (storage, &local_error);

  if (local_error == NULL)
  {
    if (current == NULL)
    {
      klass->set_functions_array (storage, functions, &local_error);
    }
    else
    {
      NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
      GVariant *var       = ncm_serialize_array_to_variant (ser, functions);
      GVariant *var_cur   = ncm_serialize_array_to_variant (ser, current);
      const gboolean same = g_variant_equal (var, var_cur);

      if (!same)
        g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_WRITE_ONCE,
                     "the function array is already set to another one");

      g_variant_unref (var);
      g_variant_unref (var_cur);
      ncm_serialize_free (ser);
    }
  }

  g_clear_pointer (&current, ncm_obj_array_unref);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_functions_array");

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_storage_get_markovian_id:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: the Markovian id, zero when none was set
 */
guint
ncm_mset_catalog_storage_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_markovian_id (storage, error);
}

/**
 * ncm_mset_catalog_storage_get_initial_sampler:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: (transfer full) (nullable): the initial sampler, %NULL when unset
 */
gchar *
ncm_mset_catalog_storage_get_initial_sampler (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_initial_sampler (storage, error);
}

/**
 * ncm_mset_catalog_storage_get_functions_array:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: (transfer full) (nullable): the #NcmObjArray of #NcmMSetFunc, %NULL when unset
 */
NcmObjArray *
ncm_mset_catalog_storage_get_functions_array (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_functions_array (storage, error);
}

/**
 * ncm_mset_catalog_storage_get_rng_inis:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: (transfer full) (nullable): the RNG state before the first row, %NULL when unset
 */
gchar *
ncm_mset_catalog_storage_get_rng_inis (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_rng_inis (storage, error);
}

/**
 * ncm_mset_catalog_storage_get_rng_stat:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: (transfer full) (nullable): the RNG state after the last synced row, %NULL when unset
 */
gchar *
ncm_mset_catalog_storage_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_rng_stat (storage, error);
}

/**
 * ncm_mset_catalog_storage_n_samplers:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns: the number of sampler segments
 */
guint
ncm_mset_catalog_storage_n_samplers (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->n_samplers (storage, error);
}

/**
 * ncm_mset_catalog_storage_get_sampler:
 * @storage: a #NcmMSetCatalogStorage
 * @i: segment index
 * @start: (out) (optional): the first row of the segment
 * @sampler: (out) (optional) (transfer full): the sampler description
 * @options: (out) (optional) (transfer full): the sampler options
 * @error: a #GError
 *
 * Returns sampler segment @i. Segment @i covers the rows from @start up to the start of
 * segment @i + 1, the last one up to the number of rows. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when @i is not below the number of segments.
 *
 * Returns: whether segment @i exists
 */
gboolean
ncm_mset_catalog_storage_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  const guint n                     = klass->n_samplers (storage, &local_error);

  if ((local_error == NULL) && (i >= n))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_get_sampler: segment %u, the storage has %u segments.", i, n);

    return FALSE;
  }

  if (local_error == NULL)
    klass->get_sampler (storage, i, start, sampler, options, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_get_sampler");

    return FALSE;
  }

  return TRUE;
}

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
 * A storage without a description is uninitialized.
 * ncm_mset_catalog_storage_set_desc() stores a description on an uninitialized storage
 * and binds a compatible one (ncm_mset_catalog_desc_is_compatible()) on an initialized
 * storage; rows are then read and written in the order of the bound description.
 *
 * The storage owns the #NcmRNG of the run. ncm_mset_catalog_storage_set_rng() gives
 * one to a storage that has neither rows nor an RNG, and records its algorithm, seed
 * and state as the initial RNG state. ncm_mset_catalog_storage_get_rng() returns it;
 * on a loaded storage it is rebuilt from the stored algorithm, seed and current state.
 * A storage with rows and no RNG never gets one.
 *
 * Rows, sampler segments and the Markovian id are written only inside a recording,
 * opened by ncm_mset_catalog_storage_start_recording() and closed by
 * ncm_mset_catalog_storage_end_recording(). One recording is open at a time. A
 * recording opens a sampler segment at the current number of rows. The current RNG
 * state is read from the storage's #NcmRNG at each sync, which happens only right
 * after a row is added (as #NcmMSetCatalogStorage:smode requires) and when the
 * recording ends. The caller should draw from the #NcmRNG only inside a recording, so
 * the stored initial state is the state before the first row and the stored current
 * state is the state after the last written row.
 *
 * The run metadata are the RNG algorithm, seed and initial state (written once by
 * ncm_mset_catalog_storage_set_rng()), the current RNG state, the Markovian id (the
 * first row of the Markovian chain, never decreasing and never beyond the rows) and
 * the sampler segments. The backend is their only holder. The #NcmMSetFunc array
 * computing the additional values is not stored: the description records the names and
 * symbols of its columns, and the functions belong to the experiment that also holds
 * the likelihood, so that objects shared between them are deserialized once.
 *
 * A backend answers the getters on an uninitialized storage without error: peek_desc
 * %NULL, nrows 0, get_rng_init with a %NULL algorithm, get_markovian_id 0 and
 * n_samplers 0. A file-backed storage whose file exists opens it at the first of these
 * calls. The Desc is stored before any run metadata: the setters require one.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_catalog_storage.h"
#include "ncm/fit/ncm_mset_catalog.h"
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
  self->smode         = NCM_MSET_CATALOG_SYNC_LEN;
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

  /*
   * The writer ends the recording; a storage released with one open keeps what the
   * last sync wrote.
   */
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
   * #NcmMSetCatalogStorage:sync-interval has elapsed since the last sync, or never
   * (only ncm_mset_catalog_storage_end_recording() syncs). A backend with nothing to
   * write (the memory storage) never syncs inside ncm_mset_catalog_storage_add_row().
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

  klass->get_rng_init     = NULL;
  klass->set_rng_init     = NULL;
  klass->get_rng_stat     = NULL;
  klass->set_rng_stat     = NULL;
  klass->get_markovian_id = NULL;
  klass->set_markovian_id = NULL;
  klass->n_samplers       = NULL;
  klass->get_sampler      = NULL;
  klass->set_sampler      = NULL;
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
 * If *@storage is not %NULL, decreases its reference count and sets *@storage to
 * %NULL.
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
 * Returns the description of @storage, %NULL when @storage is uninitialized. Sets
 * @error when the backend fails.
 *
 * Returns: (transfer none) (nullable): the #NcmMSetCatalogDesc of @storage
 */
NcmMSetCatalogDesc *
ncm_mset_catalog_storage_peek_desc (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->peek_desc (storage, error);
}

/*
 * Returns the description of storage, or sets error (UNINITIALIZED or the backend's).
 */
static NcmMSetCatalogDesc *
_ncm_mset_catalog_storage_require_desc (NcmMSetCatalogStorage *storage, GError **error)
{
  GError *local_error      = NULL;
  NcmMSetCatalogDesc *desc = ncm_mset_catalog_storage_peek_desc (storage, &local_error);

  if ((local_error == NULL) && (desc == NULL))
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED,
                 "the storage has no description");

  if (local_error != NULL)
    g_propagate_error (error, local_error);

  return desc;
}

/**
 * ncm_mset_catalog_storage_set_desc:
 * @storage: a #NcmMSetCatalogStorage
 * @desc: a #NcmMSetCatalogDesc
 * @error: a #GError
 *
 * Stores @desc on an uninitialized @storage. On an initialized @storage it binds @desc
 * when it is compatible with the stored description
 * (ncm_mset_catalog_desc_is_compatible()), so rows are read and written in the order
 * of @desc, and otherwise fails with %NCM_MSET_CATALOG_STORAGE_ERROR_INCOMPATIBLE,
 * naming the difference. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING while a
 * recording is open.
 */
void
ncm_mset_catalog_storage_set_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;
  NcmMSetCatalogDesc *stored                = NULL;

  if (self->recording)
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                 "a recording is open");
  else
    stored = ncm_mset_catalog_storage_peek_desc (storage, &local_error);

  if (local_error == NULL)
  {
    if (stored == NULL)
    {
      klass->store_desc (storage, desc, &local_error);
    }
    else
    {
      if (stored != desc)
      {
        gchar *reason = NULL;

        if (!ncm_mset_catalog_desc_is_compatible (stored, desc, &reason))
          g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_INCOMPATIBLE,
                       "the description is not compatible with the stored one: %s", reason);

        g_free (reason);
      }

      if ((local_error == NULL) && (stored != desc))
        klass->bind_desc (storage, desc, &local_error);
    }
  }

  ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_desc");
}

/**
 * ncm_mset_catalog_storage_set_rng:
 * @storage: a #NcmMSetCatalogStorage
 * @rng: the #NcmRNG of the run
 * @error: a #GError
 *
 * Makes @rng the RNG of @storage and records its algorithm, seed and state as the
 * initial RNG state, which is also the current state until the first sync. @storage
 * keeps a reference to @rng. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED
 * without a description, with %NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING while a
 * recording is open, and with %NCM_MSET_CATALOG_STORAGE_ERROR_RNG when @storage has
 * rows or an RNG.
 */
void
ncm_mset_catalog_storage_set_rng (NcmMSetCatalogStorage *storage, NcmRNG *rng, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;
  gchar *stored_algo                        = NULL;
  guint nrows                               = 0;

  if (self->recording)
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                 "a recording is open");
  else
    _ncm_mset_catalog_storage_require_desc (storage, &local_error);

  if (local_error == NULL)
    nrows = ncm_mset_catalog_storage_nrows (storage, &local_error);

  if (local_error == NULL)
    klass->get_rng_init (storage, &stored_algo, NULL, NULL, &local_error);

  if ((local_error == NULL) && ((nrows > 0) || (stored_algo != NULL)))
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG,
                 "the storage has %u rows and %s RNG", nrows, (stored_algo != NULL) ? "an" : "no");

  if (local_error == NULL)
  {
    gchar *state = ncm_rng_get_state (rng);

    klass->set_rng_init (storage, ncm_rng_get_algo (rng), ncm_rng_get_seed (rng), state, &local_error);

    if (local_error == NULL)
      klass->set_rng_stat (storage, state, &local_error);

    g_free (state);
  }

  if (local_error == NULL)
  {
    ncm_rng_clear (&self->rng);
    self->rng = ncm_rng_ref (rng);
  }

  g_free (stored_algo);

  ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_rng");
}

/*
 * Builds the storage's RNG from the stored algorithm, seed and current state, when the
 * storage has one and it was not built yet.
 */
static void
_ncm_mset_catalog_storage_load_rng (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;
  gchar *algo                               = NULL;
  gchar *inis                               = NULL;
  gchar *stat                               = NULL;
  gulong seed                               = 0;

  if (self->rng != NULL)
    return;

  klass->get_rng_init (storage, &algo, &seed, &inis, &local_error);

  if ((local_error == NULL) && (algo != NULL))
    stat = klass->get_rng_stat (storage, &local_error);

  if ((local_error == NULL) && (algo != NULL))
  {
    self->rng = ncm_rng_seeded_new (algo, seed);
    ncm_rng_set_state (self->rng, (stat != NULL) ? stat : inis);
  }

  g_free (algo);
  g_free (inis);
  g_free (stat);

  if (local_error != NULL)
    g_propagate_error (error, local_error);
}

/**
 * ncm_mset_catalog_storage_get_rng:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Returns the RNG of @storage, the same object at every call. On a loaded storage it
 * is built from the stored algorithm and seed and set to the stored current state; its
 * state is written at each sync of a recording.
 *
 * Returns: (transfer full) (nullable): the #NcmRNG of @storage, %NULL when it has none
 */
NcmRNG *
ncm_mset_catalog_storage_get_rng (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  GError *local_error                       = NULL;

  _ncm_mset_catalog_storage_load_rng (storage, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_get_rng");

    return NULL;
  }

  return (self->rng != NULL) ? ncm_rng_ref (self->rng) : NULL;
}

static void _ncm_mset_catalog_storage_open_segment (NcmMSetCatalogStorage *storage, guint nrows, const gchar *sampler, const gchar *options, GError **error);

/**
 * ncm_mset_catalog_storage_start_recording:
 * @storage: a #NcmMSetCatalogStorage
 * @sampler: the sampler description
 * @options: (nullable): the sampler options
 * @error: a #GError
 *
 * Opens a recording, the only state in which rows, sampler segments and the Markovian
 * id are written. It opens a sampler segment for @sampler and @options starting at the
 * current number of rows; a last segment that starts at the same row, and so has no
 * rows, is replaced. A stored RNG is built (ncm_mset_catalog_storage_get_rng()) so
 * that the syncs of the recording write its state.
 *
 * Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED without a description and
 * with %NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING when a recording is already open; the
 * recording is then not opened.
 */
void
ncm_mset_catalog_storage_start_recording (NcmMSetCatalogStorage *storage, const gchar *sampler, const gchar *options, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  GError *local_error                       = NULL;
  guint nrows                               = 0;

  g_assert (sampler != NULL);

  if (self->recording)
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                 "a recording is already open");
  else
    _ncm_mset_catalog_storage_require_desc (storage, &local_error);

  if (local_error == NULL)
    nrows = ncm_mset_catalog_storage_nrows (storage, &local_error);

  if (local_error == NULL)
    _ncm_mset_catalog_storage_load_rng (storage, &local_error);

  if (local_error == NULL)
    _ncm_mset_catalog_storage_open_segment (storage, nrows, sampler, options, &local_error);

  if (local_error != NULL)
  {
    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_start_recording");

    return;
  }

  self->recording = TRUE;

  g_timer_start (self->sync_timer);
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

  if (local_error != NULL)
    g_propagate_error (error, local_error);
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

  if (local_error != NULL)
    g_propagate_error (error, local_error);
}

/* Sets error (NOT_RECORDING) outside a recording. */
static void
_ncm_mset_catalog_storage_require_recording (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  if (!self->recording)
    g_set_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING,
                 "no recording is open");
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
 */
void
ncm_mset_catalog_storage_add_row (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error                       = NULL;
  NcmMSetCatalogDesc *desc                  = NULL;

  _ncm_mset_catalog_storage_require_recording (storage, &local_error);

  if (local_error == NULL)
    desc = _ncm_mset_catalog_storage_require_desc (storage, &local_error);

  if ((local_error == NULL) && (ncm_vector_len (row) != ncm_mset_catalog_desc_ncols (desc)))
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                 "row of length %u, the description has %u columns",
                 ncm_vector_len (row), ncm_mset_catalog_desc_ncols (desc));

  if (local_error == NULL)
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

  ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_add_row");
}

/**
 * ncm_mset_catalog_storage_end_recording:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Syncs, writing the state of the storage's #NcmRNG, and closes the recording. The recording
 * is closed even when the sync fails. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a recording.
 */
void
ncm_mset_catalog_storage_end_recording (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  GError *local_error                       = NULL;

  _ncm_mset_catalog_storage_require_recording (storage, &local_error);

  if (local_error == NULL)
  {
    _ncm_mset_catalog_storage_sync (storage, &local_error);

    self->recording = FALSE;
  }

  ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_end_recording");
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
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                 "row %u, the storage has %u rows", i, nrows);

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
 */
void
ncm_mset_catalog_storage_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  guint current                     = 0;
  guint nrows                       = 0;

  _ncm_mset_catalog_storage_require_recording (storage, &local_error);

  if (local_error == NULL)
    current = klass->get_markovian_id (storage, &local_error);

  if (local_error == NULL)
    nrows = ncm_mset_catalog_storage_nrows (storage, &local_error);

  if ((local_error == NULL) && ((markovian_id < current) || (markovian_id > nrows)))
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID,
                 "id %u is below the stored id %u or above the %u rows",
                 markovian_id, current, nrows);

  if (local_error == NULL)
    klass->set_markovian_id (storage, markovian_id, &local_error);

  ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_set_markovian_id");
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
 * ncm_mset_catalog_storage_get_rng_init:
 * @storage: a #NcmMSetCatalogStorage
 * @algo: (out) (optional) (nullable) (transfer full): the RNG algorithm, %NULL when no RNG is stored
 * @seed: (out) (optional): the RNG seed
 * @state: (out) (optional) (nullable) (transfer full): the RNG state before the first row
 * @error: a #GError
 *
 * Returns the RNG record written by ncm_mset_catalog_storage_set_rng().
 */
void
ncm_mset_catalog_storage_get_rng_init (NcmMSetCatalogStorage *storage, gchar **algo, gulong *seed, gchar **state, GError **error)
{
  NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage)->get_rng_init (storage, algo, seed, state, error);
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
 * %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when @i is not below the number of segments; the
 * outputs are then not set.
 */
void
ncm_mset_catalog_storage_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  GError *local_error               = NULL;
  const guint n                     = klass->n_samplers (storage, &local_error);

  if ((local_error == NULL) && (i >= n))
    g_set_error (&local_error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                 "segment %u, the storage has %u segments", i, n);

  if (local_error == NULL)
    klass->get_sampler (storage, i, start, sampler, options, &local_error);

  ncm_util_forward_or_call_error (error, local_error, "ncm_mset_catalog_storage_get_sampler");
}


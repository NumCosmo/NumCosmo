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
 * ncm_mset_catalog_storage_erase() removes the rows, the sampler segments and the
 * Markovian id outside a recording; the description and the RNG stay.
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
  NcmMSetCatalogDesc *desc = ncm_mset_catalog_storage_peek_desc (storage, error);

  NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

  if (desc == NULL)
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED,
                                "the storage has no description");

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
  NcmMSetCatalogDesc *stored;

  if (self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                                "ncm_mset_catalog_storage_set_desc: a recording is open");

    return;
  }

  stored = ncm_mset_catalog_storage_peek_desc (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_desc: ");

  if (stored == NULL)
  {
    klass->store_desc (storage, desc, error);
    NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_desc: ");

    return;
  }

  if (stored == desc)
    return;

  {
    gchar *reason = NULL;

    if (!ncm_mset_catalog_desc_is_compatible (stored, desc, &reason))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_INCOMPATIBLE,
                                  "ncm_mset_catalog_storage_set_desc: the description is not compatible with the stored one: %s",
                                  reason);
      g_free (reason);

      return;
    }
  }

  klass->bind_desc (storage, desc, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_desc: ");
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
  gchar *stored_algo                        = NULL;
  gchar *state;
  guint nrows;

  if (self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                                "ncm_mset_catalog_storage_set_rng: a recording is open");

    return;
  }

  _ncm_mset_catalog_storage_require_desc (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_rng: ");

  nrows = ncm_mset_catalog_storage_nrows (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_rng: ");

  klass->get_rng_init (storage, &stored_algo, NULL, NULL, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, g_free (stored_algo), , "ncm_mset_catalog_storage_set_rng: ");

  if ((nrows > 0) || (stored_algo != NULL))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG,
                                "ncm_mset_catalog_storage_set_rng: the storage has %u rows and %s RNG",
                                nrows, (stored_algo != NULL) ? "an" : "no");
    g_free (stored_algo);

    return;
  }

  state = ncm_rng_get_state (rng);

  klass->set_rng_init (storage, ncm_rng_get_algo (rng), ncm_rng_get_seed (rng), state, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, g_free (state), , "ncm_mset_catalog_storage_set_rng: ");

  klass->set_rng_stat (storage, state, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, g_free (state), , "ncm_mset_catalog_storage_set_rng: ");

  g_free (state);

  ncm_rng_clear (&self->rng);
  self->rng = ncm_rng_ref (rng);
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
  gchar *algo                               = NULL;
  gchar *inis                               = NULL;
  gchar *stat                               = NULL;
  gulong seed                               = 0;

  if (self->rng != NULL)
    return;

  klass->get_rng_init (storage, &algo, &seed, &inis, error);
  NCM_UTIL_ON_ERROR_RETURN (error, { g_free (algo);
                                     g_free (inis);
                            }, );

  if (algo == NULL)
  {
    g_free (inis);

    return;
  }

  stat = klass->get_rng_stat (storage, error);
  NCM_UTIL_ON_ERROR_RETURN (error, { g_free (algo);
                                     g_free (inis);
                                     g_free (stat);
                            }, );

  self->rng = ncm_rng_seeded_new (algo, seed);
  ncm_rng_set_state (self->rng, (stat != NULL) ? stat : inis);

  g_free (algo);
  g_free (inis);
  g_free (stat);
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

  _ncm_mset_catalog_storage_load_rng (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_catalog_storage_get_rng: ");

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
  guint nrows;

  g_assert (sampler != NULL);

  if (self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                                "ncm_mset_catalog_storage_start_recording: a recording is already open");

    return;
  }

  _ncm_mset_catalog_storage_require_desc (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_start_recording: ");

  nrows = ncm_mset_catalog_storage_nrows (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_start_recording: ");

  _ncm_mset_catalog_storage_load_rng (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_start_recording: ");

  _ncm_mset_catalog_storage_open_segment (storage, nrows, sampler, options, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_start_recording: ");

  self->recording = TRUE;

  g_timer_start (self->sync_timer);
}

/*
 * Writes the segment (nrows, sampler, options): over the last segment when it starts
 * at nrows, after it otherwise.
 */
static void
_ncm_mset_catalog_storage_open_segment (NcmMSetCatalogStorage *storage, guint nrows, const gchar *sampler, const gchar *options, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  const guint n                     = klass->n_samplers (storage, error);
  guint i                           = n;

  NCM_UTIL_ON_ERROR_RETURN (error, , );

  if (n > 0)
  {
    guint last_start = 0;

    klass->get_sampler (storage, n - 1, &last_start, NULL, NULL, error);
    NCM_UTIL_ON_ERROR_RETURN (error, , );

    if (last_start == nrows)
      i = n - 1;
  }

  klass->set_sampler (storage, i, nrows, sampler, (options != NULL) ? options : "", error);
}

/*
 * Writes the state of the recording's RNG, then lets the backend write what it holds.
 */
static void
_ncm_mset_catalog_storage_sync (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);

  if (self->rng != NULL)
  {
    gchar *state = ncm_rng_get_state (self->rng);

    klass->set_rng_stat (storage, state, error);
    g_free (state);
    NCM_UTIL_ON_ERROR_RETURN (error, , );
  }

  if (klass->sync != NULL)
    klass->sync (storage, error);
}

/* Sets error (NOT_RECORDING) outside a recording. */
static void
_ncm_mset_catalog_storage_require_recording (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  if (!self->recording)
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING,
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
 * %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when the length of @row is not the number of
 * columns.
 */
void
ncm_mset_catalog_storage_add_row (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  NcmMSetCatalogDesc *desc;

  _ncm_mset_catalog_storage_require_recording (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_add_row: ");

  desc = _ncm_mset_catalog_storage_require_desc (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_add_row: ");

  if (ncm_vector_len (row) != ncm_mset_catalog_desc_ncols (desc))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_add_row: row of length %u, the description has %u columns",
                                ncm_vector_len (row), ncm_mset_catalog_desc_ncols (desc));

    return;
  }

  klass->append (storage, row, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_add_row: ");

  if (klass->sync != NULL)
  {
    switch (self->smode)
    {
      case NCM_MSET_CATALOG_SYNC_AUTO:
        _ncm_mset_catalog_storage_sync (storage, error);
        break;
      case NCM_MSET_CATALOG_SYNC_TIMED:

        if (g_timer_elapsed (self->sync_timer, NULL) > self->sync_interval)
        {
          g_timer_start (self->sync_timer);
          _ncm_mset_catalog_storage_sync (storage, error);
        }

        break;
      default:
        break;
    }

    NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_add_row: ");
  }
}

/**
 * ncm_mset_catalog_storage_end_recording:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Syncs, writing the state of the storage's #NcmRNG, and closes the recording. The
 * recording is closed even when the sync fails. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING outside a recording. This is where a
 * file-backed storage reports a failure to write the rows of the recording; a storage
 * dropped inside a recording writes them in dispose, which can only warn. Any holder of
 * the storage may call it.
 */
void
ncm_mset_catalog_storage_end_recording (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);

  _ncm_mset_catalog_storage_require_recording (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_end_recording: ");

  _ncm_mset_catalog_storage_sync (storage, error);
  self->recording = FALSE;
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_end_recording: ");
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
 * Returns row @i in the order of the bound description. The returned vector must not
 * be modified. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when @i is not below the
 * number of rows.
 *
 * Returns: (transfer full) (nullable): row @i, or %NULL on error
 */
NcmVector *
ncm_mset_catalog_storage_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  const guint nrows                 = klass->nrows (storage, error);
  NcmVector *row;

  NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_catalog_storage_get_row: ");

  if (i >= nrows)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_get_row: row %u, the storage has %u rows", i, nrows);

    return NULL;
  }

  row = klass->get_row (storage, i, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, ncm_vector_clear (&row), NULL, "ncm_mset_catalog_storage_get_row: ");

  return row;
}

/**
 * ncm_mset_catalog_storage_peek_row:
 * @storage: a #NcmMSetCatalogStorage
 * @i: row index
 * @error: a #GError
 *
 * Returns row @i in the order of the bound description without a new reference. The
 * vector must not be modified, and it is readable until the next call on @storage that
 * reads another row (ncm_mset_catalog_storage_peek_row() or
 * ncm_mset_catalog_storage_get_row() with another index), adds a row, syncs, erases,
 * sets the description, or until @storage is disposed. Concurrent readers use
 * ncm_mset_catalog_storage_get_row(). Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_ROW
 * when @i is not below the number of rows.
 *
 * Returns: (transfer none) (nullable): row @i, or %NULL on error
 */
NcmVector *
ncm_mset_catalog_storage_peek_row (NcmMSetCatalogStorage *storage, guint i, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  const guint nrows                 = klass->nrows (storage, error);
  NcmVector *row;

  NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_catalog_storage_peek_row: ");

  if (i >= nrows)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_peek_row: row %u, the storage has %u rows", i, nrows);

    return NULL;
  }

  row = klass->peek_row (storage, i, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_catalog_storage_peek_row: ");

  return row;
}

/**
 * ncm_mset_catalog_storage_erase:
 * @storage: a #NcmMSetCatalogStorage
 * @error: a #GError
 *
 * Removes every row, every sampler segment and the Markovian id. The description and
 * the RNG stay: the stored initial RNG state remains the state before the first row
 * ever added, and the current RNG state is left as written by the last sync. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED without a description and with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING inside a recording.
 */
void
ncm_mset_catalog_storage_erase (NcmMSetCatalogStorage *storage, GError **error)
{
  NcmMSetCatalogStoragePrivate * const self = ncm_mset_catalog_storage_get_instance_private (storage);
  NcmMSetCatalogStorageClass *klass         = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);

  _ncm_mset_catalog_storage_require_desc (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_erase: ");

  if (self->recording)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
                                "ncm_mset_catalog_storage_erase: a recording is open");

    return;
  }

  klass->erase (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_erase: ");
}

/**
 * ncm_mset_catalog_storage_set_markovian_id:
 * @storage: a #NcmMSetCatalogStorage
 * @markovian_id: id of the first row of the Markovian chain
 * @error: a #GError
 *
 * Sets the Markovian id. Fails with %NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING
 * outside a recording and with %NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID when
 * @markovian_id is below the stored one or above the number of rows.
 */
void
ncm_mset_catalog_storage_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  guint current;
  guint nrows;

  _ncm_mset_catalog_storage_require_recording (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_markovian_id: ");

  current = klass->get_markovian_id (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_markovian_id: ");

  nrows = klass->nrows (storage, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_markovian_id: ");

  if ((markovian_id < current) || (markovian_id > nrows))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID,
                                "ncm_mset_catalog_storage_set_markovian_id: id %u is below the stored id %u or above the %u rows",
                                markovian_id, current, nrows);

    return;
  }

  klass->set_markovian_id (storage, markovian_id, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_set_markovian_id: ");
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
 * Returns sampler segment @i. Segment @i covers the rows from @start up to the start
 * of segment @i + 1, the last one up to the number of rows. Fails with
 * %NCM_MSET_CATALOG_STORAGE_ERROR_ROW when @i is not below the number of segments; the
 * outputs are then not set.
 */
void
ncm_mset_catalog_storage_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error)
{
  NcmMSetCatalogStorageClass *klass = NCM_MSET_CATALOG_STORAGE_GET_CLASS (storage);
  const guint n                     = klass->n_samplers (storage, error);

  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_get_sampler: ");

  if (i >= n)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
                                "ncm_mset_catalog_storage_get_sampler: segment %u, the storage has %u segments", i, n);

    return;
  }

  klass->get_sampler (storage, i, start, sampler, options, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, , , "ncm_mset_catalog_storage_get_sampler: ");
}


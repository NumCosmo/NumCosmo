/***************************************************************************
 *            ncm_mset_catalog_storage.h
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_storage.h
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

#ifndef _NCM_MSET_CATALOG_STORAGE_H_
#define _NCM_MSET_CATALOG_STORAGE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/fit/ncm_mset_catalog_desc.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmMSetCatalogStorage, ncm_mset_catalog_storage, NCM, MSET_CATALOG_STORAGE, GObject)

/**
 * NcmMSetCatalogStorageError:
 * @NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED: the storage has no #NcmMSetCatalogDesc
 * @NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING: a recording is open
 * @NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING: a write outside a recording
 * @NCM_MSET_CATALOG_STORAGE_ERROR_ROW: a row or sampler-segment index out of range, or a row
 *   of the wrong length
 * @NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID: a Markovian id below the stored one or beyond
 *   the rows
 * @NCM_MSET_CATALOG_STORAGE_ERROR_WRITE_ONCE: a write-once item already holds another value
 * @NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND: the backend failed
 *
 * Error codes returned by the #NcmMSetCatalogStorage class.
 *
 */
typedef enum _NcmMSetCatalogStorageError /*< prefix=NCM_MSET_CATALOG_STORAGE_ERROR >*/
{
  NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED,
  NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING,
  NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING,
  NCM_MSET_CATALOG_STORAGE_ERROR_ROW,
  NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID,
  NCM_MSET_CATALOG_STORAGE_ERROR_WRITE_ONCE,
  NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND,
} NcmMSetCatalogStorageError;

GQuark ncm_mset_catalog_storage_error_quark (void);

#define NCM_MSET_CATALOG_STORAGE_ERROR (ncm_mset_catalog_storage_error_quark ())

struct _NcmMSetCatalogStorageClass
{
  /*< private >*/
  GObjectClass parent_class;

  NcmMSetCatalogDesc *(*peek_desc) (NcmMSetCatalogStorage *storage, GError **error);
  void (*store_desc) (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);
  void (*bind_desc) (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);
  guint (*nrows) (NcmMSetCatalogStorage *storage, GError **error);
  void (*append) (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error);
  NcmVector *(*get_row) (NcmMSetCatalogStorage *storage, guint i, GError **error);
  void (*sync) (NcmMSetCatalogStorage *storage, GError **error);
  gchar *(*get_rng_inis) (NcmMSetCatalogStorage *storage, GError **error);
  void (*set_rng_inis) (NcmMSetCatalogStorage *storage, const gchar *state, GError **error);
  gchar *(*get_rng_stat) (NcmMSetCatalogStorage *storage, GError **error);
  void (*set_rng_stat) (NcmMSetCatalogStorage *storage, const gchar *state, GError **error);
  guint (*get_markovian_id) (NcmMSetCatalogStorage *storage, GError **error);
  void (*set_markovian_id) (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error);
  gchar *(*get_initial_sampler) (NcmMSetCatalogStorage *storage, GError **error);
  void (*set_initial_sampler) (NcmMSetCatalogStorage *storage, const gchar *sampler, GError **error);
  NcmObjArray *(*get_functions_array) (NcmMSetCatalogStorage *storage, GError **error);
  void (*set_functions_array) (NcmMSetCatalogStorage *storage, NcmObjArray *functions, GError **error);
  guint (*n_samplers) (NcmMSetCatalogStorage *storage, GError **error);
  void (*get_sampler) (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error);
  void (*set_sampler) (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options, GError **error);
};

NcmMSetCatalogStorage *ncm_mset_catalog_storage_ref (NcmMSetCatalogStorage *storage);
void ncm_mset_catalog_storage_free (NcmMSetCatalogStorage *storage);
void ncm_mset_catalog_storage_clear (NcmMSetCatalogStorage **storage);

NcmMSetCatalogDesc *ncm_mset_catalog_storage_peek_desc (NcmMSetCatalogStorage *storage, GError **error);
gboolean ncm_mset_catalog_storage_set_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);

gboolean ncm_mset_catalog_storage_start_recording (NcmMSetCatalogStorage *storage, const gchar *sampler, const gchar *options, NcmRNG *rng, GError **error);
gboolean ncm_mset_catalog_storage_add_row (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error);
gboolean ncm_mset_catalog_storage_end_recording (NcmMSetCatalogStorage *storage, GError **error);
gboolean ncm_mset_catalog_storage_is_recording (NcmMSetCatalogStorage *storage);

guint ncm_mset_catalog_storage_nrows (NcmMSetCatalogStorage *storage, GError **error);
NcmVector *ncm_mset_catalog_storage_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error);

gboolean ncm_mset_catalog_storage_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error);
gboolean ncm_mset_catalog_storage_set_initial_sampler (NcmMSetCatalogStorage *storage, const gchar *sampler, GError **error);
gboolean ncm_mset_catalog_storage_set_functions_array (NcmMSetCatalogStorage *storage, NcmObjArray *functions, GError **error);

guint ncm_mset_catalog_storage_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error);
gchar *ncm_mset_catalog_storage_get_initial_sampler (NcmMSetCatalogStorage *storage, GError **error);
NcmObjArray *ncm_mset_catalog_storage_get_functions_array (NcmMSetCatalogStorage *storage, GError **error);
gchar *ncm_mset_catalog_storage_get_rng_inis (NcmMSetCatalogStorage *storage, GError **error);
gchar *ncm_mset_catalog_storage_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error);
guint ncm_mset_catalog_storage_n_samplers (NcmMSetCatalogStorage *storage, GError **error);
gboolean ncm_mset_catalog_storage_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error);

G_END_DECLS

#endif /* _NCM_MSET_CATALOG_STORAGE_H_ */

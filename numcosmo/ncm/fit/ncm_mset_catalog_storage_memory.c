/***************************************************************************
 *            ncm_mset_catalog_storage_memory.c
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_storage_memory.c
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
 * NcmMSetCatalogStorageMemory:
 *
 * An #NcmMSetCatalogStorage that keeps the description, every row and the run metadata
 * in memory.
 *
 * All of them are properties, so serialization writes the rows: it is slow for a long
 * catalog, which should use a file-backed storage. Rows are held in the order of the
 * bound description; binding a description whose free parameters are in another order
 * replaces every row by a permuted copy. The memory storage has nothing to write at a
 * sync.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_catalog_storage_memory.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NCM_MSET_CATALOG_STORAGE_MEMORY_SAMPLERS_TYPE "a(uss)"

typedef struct _NcmMSetCatalogStorageMemorySegment
{
  guint start;
  gchar *sampler;
  gchar *options;
} NcmMSetCatalogStorageMemorySegment;

struct _NcmMSetCatalogStorageMemory
{
  NcmMSetCatalogStorage parent_instance;
  NcmMSetCatalogDesc *desc;
  GPtrArray *rows;
  gchar *rng_algo;
  gulong rng_seed;
  gchar *rng_inis;
  gchar *rng_stat;
  guint markovian_id;
  GArray *samplers;
};

enum
{
  PROP_0,
  PROP_DESC,
  PROP_ROWS,
  PROP_RNG_ALGO,
  PROP_RNG_SEED,
  PROP_RNG_INIS,
  PROP_RNG_STAT,
  PROP_MARKOVIAN_ID,
  PROP_SAMPLERS,
};

G_DEFINE_TYPE (NcmMSetCatalogStorageMemory, ncm_mset_catalog_storage_memory, NCM_TYPE_MSET_CATALOG_STORAGE)

static void
_ncm_mset_catalog_storage_memory_segment_clear (gpointer data)
{
  NcmMSetCatalogStorageMemorySegment *segment = data;

  g_clear_pointer (&segment->sampler, g_free);
  g_clear_pointer (&segment->options, g_free);
}

static void
ncm_mset_catalog_storage_memory_init (NcmMSetCatalogStorageMemory *storage_mem)
{
  storage_mem->desc         = NULL;
  storage_mem->rows         = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  storage_mem->rng_algo     = NULL;
  storage_mem->rng_seed     = 0;
  storage_mem->rng_inis     = NULL;
  storage_mem->rng_stat     = NULL;
  storage_mem->markovian_id = 0;
  storage_mem->samplers     = g_array_new (FALSE, FALSE, sizeof (NcmMSetCatalogStorageMemorySegment));

  g_array_set_clear_func (storage_mem->samplers, _ncm_mset_catalog_storage_memory_segment_clear);
}

static void
_ncm_mset_catalog_storage_memory_set_rows (NcmMSetCatalogStorageMemory *storage_mem, NcmObjArray *rows)
{
  g_ptr_array_set_size (storage_mem->rows, 0);

  if (rows != NULL)
  {
    guint i;

    for (i = 0; i < rows->len; i++)
      g_ptr_array_add (storage_mem->rows, ncm_vector_ref (NCM_VECTOR (ncm_obj_array_peek (rows, i))));
  }
}

static NcmObjArray *
_ncm_mset_catalog_storage_memory_get_rows (NcmMSetCatalogStorageMemory *storage_mem)
{
  NcmObjArray *rows = ncm_obj_array_new ();
  guint i;

  for (i = 0; i < storage_mem->rows->len; i++)
    ncm_obj_array_add (rows, G_OBJECT (g_ptr_array_index (storage_mem->rows, i)));

  return rows;
}

static void
_ncm_mset_catalog_storage_memory_set_samplers (NcmMSetCatalogStorageMemory *storage_mem, GVariant *samplers)
{
  g_array_set_size (storage_mem->samplers, 0);

  if (samplers != NULL)
  {
    const guint n = g_variant_n_children (samplers);
    guint i;

    for (i = 0; i < n; i++)
    {
      NcmMSetCatalogStorageMemorySegment segment;

      g_variant_get_child (samplers, i, "(uss)", &segment.start, &segment.sampler, &segment.options);
      g_array_append_val (storage_mem->samplers, segment);
    }
  }
}

static GVariant *
_ncm_mset_catalog_storage_memory_get_samplers (NcmMSetCatalogStorageMemory *storage_mem)
{
  GVariantBuilder builder;
  guint i;

  g_variant_builder_init (&builder, G_VARIANT_TYPE (NCM_MSET_CATALOG_STORAGE_MEMORY_SAMPLERS_TYPE));

  for (i = 0; i < storage_mem->samplers->len; i++)
  {
    NcmMSetCatalogStorageMemorySegment *segment = &g_array_index (storage_mem->samplers, NcmMSetCatalogStorageMemorySegment, i);

    g_variant_builder_add (&builder, "(uss)", segment->start, segment->sampler, segment->options);
  }

  return g_variant_ref_sink (g_variant_builder_end (&builder));
}

static void
_ncm_mset_catalog_storage_memory_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (object);

  g_return_if_fail (NCM_IS_MSET_CATALOG_STORAGE_MEMORY (object));

  switch (prop_id)
  {
    case PROP_DESC:
      storage_mem->desc = g_value_dup_object (value);
      break;
    case PROP_ROWS:
      _ncm_mset_catalog_storage_memory_set_rows (storage_mem, g_value_get_boxed (value));
      break;
    case PROP_RNG_ALGO:
      storage_mem->rng_algo = g_value_dup_string (value);
      break;
    case PROP_RNG_SEED:
      storage_mem->rng_seed = g_value_get_ulong (value);
      break;
    case PROP_RNG_INIS:
      storage_mem->rng_inis = g_value_dup_string (value);
      break;
    case PROP_RNG_STAT:
      storage_mem->rng_stat = g_value_dup_string (value);
      break;
    case PROP_MARKOVIAN_ID:
      storage_mem->markovian_id = g_value_get_uint (value);
      break;
    case PROP_SAMPLERS:
      _ncm_mset_catalog_storage_memory_set_samplers (storage_mem, g_value_get_variant (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_storage_memory_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (object);

  g_return_if_fail (NCM_IS_MSET_CATALOG_STORAGE_MEMORY (object));

  switch (prop_id)
  {
    case PROP_DESC:
      g_value_set_object (value, storage_mem->desc);
      break;
    case PROP_ROWS:
      g_value_take_boxed (value, _ncm_mset_catalog_storage_memory_get_rows (storage_mem));
      break;
    case PROP_RNG_ALGO:
      g_value_set_string (value, storage_mem->rng_algo);
      break;
    case PROP_RNG_SEED:
      g_value_set_ulong (value, storage_mem->rng_seed);
      break;
    case PROP_RNG_INIS:
      g_value_set_string (value, storage_mem->rng_inis);
      break;
    case PROP_RNG_STAT:
      g_value_set_string (value, storage_mem->rng_stat);
      break;
    case PROP_MARKOVIAN_ID:
      g_value_set_uint (value, storage_mem->markovian_id);
      break;
    case PROP_SAMPLERS:
      g_value_take_variant (value, _ncm_mset_catalog_storage_memory_get_samplers (storage_mem));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_storage_memory_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_memory_parent_class)->constructed (object);
  {
    NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (object);
    const guint nrows                        = storage_mem->rows->len;
    guint i;

    /* Built from properties there is no error location: inconsistent properties abort. */
    if ((nrows > 0) && (storage_mem->desc == NULL))
      g_error ("_ncm_mset_catalog_storage_memory_constructed: %u rows without a description.", nrows);

    for (i = 0; i < nrows; i++)
    {
      NcmVector *row = g_ptr_array_index (storage_mem->rows, i);

      if (ncm_vector_len (row) != ncm_mset_catalog_desc_ncols (storage_mem->desc))
        g_error ("_ncm_mset_catalog_storage_memory_constructed: row %u has length %u, the description has %u columns.",
                 i, ncm_vector_len (row), ncm_mset_catalog_desc_ncols (storage_mem->desc));
    }

    if ((storage_mem->rng_algo == NULL) && ((storage_mem->rng_inis != NULL) || (storage_mem->rng_stat != NULL)))
      g_error ("_ncm_mset_catalog_storage_memory_constructed: an RNG state without an RNG algorithm.");

    if ((storage_mem->rng_algo != NULL) && (storage_mem->rng_inis == NULL))
      g_error ("_ncm_mset_catalog_storage_memory_constructed: an RNG algorithm without an initial RNG state.");

    if (storage_mem->markovian_id > nrows)
      g_error ("_ncm_mset_catalog_storage_memory_constructed: Markovian id %u above the %u rows.",
               storage_mem->markovian_id, nrows);

    for (i = 0; i < storage_mem->samplers->len; i++)
    {
      const guint start = g_array_index (storage_mem->samplers, NcmMSetCatalogStorageMemorySegment, i).start;

      if ((start > nrows) ||
          ((i > 0) && (start <= g_array_index (storage_mem->samplers, NcmMSetCatalogStorageMemorySegment, i - 1).start)))
        g_error ("_ncm_mset_catalog_storage_memory_constructed: sampler segment %u starts at %u, not after the "
                 "previous segment and within the %u rows.", i, start, nrows);
    }
  }
}

static void
_ncm_mset_catalog_storage_memory_dispose (GObject *object)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (object);

  ncm_mset_catalog_desc_clear (&storage_mem->desc);
  g_ptr_array_set_size (storage_mem->rows, 0);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_memory_parent_class)->dispose (object);
}

static void
_ncm_mset_catalog_storage_memory_finalize (GObject *object)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (object);

  g_clear_pointer (&storage_mem->rows, g_ptr_array_unref);
  g_clear_pointer (&storage_mem->samplers, g_array_unref);
  g_clear_pointer (&storage_mem->rng_algo, g_free);
  g_clear_pointer (&storage_mem->rng_inis, g_free);
  g_clear_pointer (&storage_mem->rng_stat, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_storage_memory_parent_class)->finalize (object);
}

static NcmMSetCatalogDesc *_ncm_mset_catalog_storage_memory_peek_desc (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_memory_store_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);
static void _ncm_mset_catalog_storage_memory_bind_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error);
static guint _ncm_mset_catalog_storage_memory_nrows (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_memory_append (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error);
static NcmVector *_ncm_mset_catalog_storage_memory_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error);
static void _ncm_mset_catalog_storage_memory_get_rng_init (NcmMSetCatalogStorage *storage, gchar **algo, gulong *seed, gchar **state, GError **error);
static void _ncm_mset_catalog_storage_memory_set_rng_init (NcmMSetCatalogStorage *storage, const gchar *algo, gulong seed, const gchar *state, GError **error);
static gchar *_ncm_mset_catalog_storage_memory_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_memory_set_rng_stat (NcmMSetCatalogStorage *storage, const gchar *state, GError **error);
static guint _ncm_mset_catalog_storage_memory_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_memory_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error);
static guint _ncm_mset_catalog_storage_memory_n_samplers (NcmMSetCatalogStorage *storage, GError **error);
static void _ncm_mset_catalog_storage_memory_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error);
static void _ncm_mset_catalog_storage_memory_set_sampler (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options, GError **error);

static void
ncm_mset_catalog_storage_memory_class_init (NcmMSetCatalogStorageMemoryClass *klass)
{
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcmMSetCatalogStorageClass *storage_class = NCM_MSET_CATALOG_STORAGE_CLASS (klass);

  object_class->set_property = &_ncm_mset_catalog_storage_memory_set_property;
  object_class->get_property = &_ncm_mset_catalog_storage_memory_get_property;
  object_class->constructed  = &_ncm_mset_catalog_storage_memory_constructed;
  object_class->dispose      = &_ncm_mset_catalog_storage_memory_dispose;
  object_class->finalize     = &_ncm_mset_catalog_storage_memory_finalize;

  /**
   * NcmMSetCatalogStorageMemory:desc:
   *
   * The stored #NcmMSetCatalogDesc, %NULL while the storage is uninitialized.
   */
  g_object_class_install_property (object_class,
                                   PROP_DESC,
                                   g_param_spec_object ("desc",
                                                        NULL,
                                                        "Catalog description",
                                                        NCM_TYPE_MSET_CATALOG_DESC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:rows:
   *
   * The rows, each a #NcmVector in the order of #NcmMSetCatalogStorageMemory:desc.
   */
  g_object_class_install_property (object_class,
                                   PROP_ROWS,
                                   g_param_spec_boxed ("rows",
                                                       NULL,
                                                       "Catalog rows",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:rng-algo:
   *
   * The algorithm of the run's #NcmRNG, %NULL when no RNG was set.
   */
  g_object_class_install_property (object_class,
                                   PROP_RNG_ALGO,
                                   g_param_spec_string ("rng-algo",
                                                        NULL,
                                                        "RNG algorithm",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:rng-seed:
   *
   * The seed of the run's #NcmRNG.
   */
  g_object_class_install_property (object_class,
                                   PROP_RNG_SEED,
                                   g_param_spec_ulong ("rng-seed",
                                                       NULL,
                                                       "RNG seed",
                                                       0, G_MAXULONG, 0,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:rng-inis:
   *
   * The state of the run's #NcmRNG before the first row (ncm_rng_get_state()), %NULL
   * when no RNG was set. It is set together with
   * #NcmMSetCatalogStorageMemory:rng-algo.
   */
  g_object_class_install_property (object_class,
                                   PROP_RNG_INIS,
                                   g_param_spec_string ("rng-inis",
                                                        NULL,
                                                        "RNG state before the first row",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:rng-stat:
   *
   * The state of the run's #NcmRNG after the last row, written at the end of each
   * recording; %NULL when no RNG was set.
   */
  g_object_class_install_property (object_class,
                                   PROP_RNG_STAT,
                                   g_param_spec_string ("rng-stat",
                                                        NULL,
                                                        "RNG state after the last synced row",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:markovian-id:
   *
   * The first row of the Markovian chain, at most the number of rows.
   */
  g_object_class_install_property (object_class,
                                   PROP_MARKOVIAN_ID,
                                   g_param_spec_uint ("markovian-id",
                                                      NULL,
                                                      "First row id of the Markovian chain",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogStorageMemory:samplers:
   *
   * The sampler segments, each (first row, sampler, options), in order of first row.
   */
  g_object_class_install_property (object_class,
                                   PROP_SAMPLERS,
                                   g_param_spec_variant ("samplers",
                                                         NULL,
                                                         "Sampler segments",
                                                         G_VARIANT_TYPE (NCM_MSET_CATALOG_STORAGE_MEMORY_SAMPLERS_TYPE),
                                                         NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  storage_class->peek_desc        = &_ncm_mset_catalog_storage_memory_peek_desc;
  storage_class->store_desc       = &_ncm_mset_catalog_storage_memory_store_desc;
  storage_class->bind_desc        = &_ncm_mset_catalog_storage_memory_bind_desc;
  storage_class->nrows            = &_ncm_mset_catalog_storage_memory_nrows;
  storage_class->append           = &_ncm_mset_catalog_storage_memory_append;
  storage_class->get_row          = &_ncm_mset_catalog_storage_memory_get_row;
  storage_class->get_rng_init     = &_ncm_mset_catalog_storage_memory_get_rng_init;
  storage_class->set_rng_init     = &_ncm_mset_catalog_storage_memory_set_rng_init;
  storage_class->get_rng_stat     = &_ncm_mset_catalog_storage_memory_get_rng_stat;
  storage_class->set_rng_stat     = &_ncm_mset_catalog_storage_memory_set_rng_stat;
  storage_class->get_markovian_id = &_ncm_mset_catalog_storage_memory_get_markovian_id;
  storage_class->set_markovian_id = &_ncm_mset_catalog_storage_memory_set_markovian_id;
  storage_class->n_samplers       = &_ncm_mset_catalog_storage_memory_n_samplers;
  storage_class->get_sampler      = &_ncm_mset_catalog_storage_memory_get_sampler;
  storage_class->set_sampler      = &_ncm_mset_catalog_storage_memory_set_sampler;
}

static NcmMSetCatalogDesc *
_ncm_mset_catalog_storage_memory_peek_desc (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->desc;
}

static void
_ncm_mset_catalog_storage_memory_store_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (storage);

  g_assert_null (storage_mem->desc);
  storage_mem->desc = ncm_mset_catalog_desc_ref (desc);
}

/*
 * Replaces each row by a copy with the free parameters in the order of desc, then
 * binds desc. The two descriptions are compatible, so they have the same additional
 * values and the same free parameters.
 */
static void
_ncm_mset_catalog_storage_memory_bind_desc (NcmMSetCatalogStorage *storage, NcmMSetCatalogDesc *desc, GError **error)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (storage);
  NcmMSet *mset_stored                     = ncm_mset_catalog_desc_peek_mset (storage_mem->desc);
  NcmMSet *mset_bound                      = ncm_mset_catalog_desc_peek_mset (desc);
  const guint nadd                         = ncm_mset_catalog_desc_nadd_vals (desc);
  const guint nfree                        = ncm_mset_fparams_len (mset_bound);
  GArray *perm                             = g_array_sized_new (FALSE, FALSE, sizeof (guint), nfree);
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

  for (i = 0; i < storage_mem->rows->len; i++)
  {
    NcmVector *row     = g_ptr_array_index (storage_mem->rows, i);
    NcmVector *new_row = ncm_vector_dup (row);

    for (j = 0; j < nfree; j++)
      ncm_vector_set (new_row, nadd + j, ncm_vector_get (row, nadd + g_array_index (perm, guint, j)));

    ncm_vector_free (row);
    g_ptr_array_index (storage_mem->rows, i) = new_row;
  }

  g_array_unref (perm);

  ncm_mset_catalog_desc_clear (&storage_mem->desc);
  storage_mem->desc = ncm_mset_catalog_desc_ref (desc);
}

static guint
_ncm_mset_catalog_storage_memory_nrows (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->rows->len;
}

static void
_ncm_mset_catalog_storage_memory_append (NcmMSetCatalogStorage *storage, NcmVector *row, GError **error)
{
  g_ptr_array_add (NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->rows, ncm_vector_dup (row));
}

static NcmVector *
_ncm_mset_catalog_storage_memory_get_row (NcmMSetCatalogStorage *storage, guint i, GError **error)
{
  return ncm_vector_ref (g_ptr_array_index (NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->rows, i));
}

static void
_ncm_mset_catalog_storage_memory_get_rng_init (NcmMSetCatalogStorage *storage, gchar **algo, gulong *seed, gchar **state, GError **error)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (storage);

  if (algo != NULL)
    *algo = g_strdup (storage_mem->rng_algo);

  if (seed != NULL)
    *seed = storage_mem->rng_seed;

  if (state != NULL)
    *state = g_strdup (storage_mem->rng_inis);
}

static void
_ncm_mset_catalog_storage_memory_set_rng_init (NcmMSetCatalogStorage *storage, const gchar *algo, gulong seed, const gchar *state, GError **error)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (storage);

  g_free (storage_mem->rng_algo);
  g_free (storage_mem->rng_inis);

  storage_mem->rng_algo = g_strdup (algo);
  storage_mem->rng_seed = seed;
  storage_mem->rng_inis = g_strdup (state);
}

static gchar *
_ncm_mset_catalog_storage_memory_get_rng_stat (NcmMSetCatalogStorage *storage, GError **error)
{
  return g_strdup (NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->rng_stat);
}

static void
_ncm_mset_catalog_storage_memory_set_rng_stat (NcmMSetCatalogStorage *storage, const gchar *state, GError **error)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (storage);

  g_free (storage_mem->rng_stat);
  storage_mem->rng_stat = g_strdup (state);
}

static guint
_ncm_mset_catalog_storage_memory_get_markovian_id (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->markovian_id;
}

static void
_ncm_mset_catalog_storage_memory_set_markovian_id (NcmMSetCatalogStorage *storage, guint markovian_id, GError **error)
{
  NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->markovian_id = markovian_id;
}

static guint
_ncm_mset_catalog_storage_memory_n_samplers (NcmMSetCatalogStorage *storage, GError **error)
{
  return NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->samplers->len;
}

static void
_ncm_mset_catalog_storage_memory_get_sampler (NcmMSetCatalogStorage *storage, guint i, guint *start, gchar **sampler, gchar **options, GError **error)
{
  NcmMSetCatalogStorageMemorySegment *segment = &g_array_index (NCM_MSET_CATALOG_STORAGE_MEMORY (storage)->samplers, NcmMSetCatalogStorageMemorySegment, i);

  if (start != NULL)
    *start = segment->start;

  if (sampler != NULL)
    *sampler = g_strdup (segment->sampler);

  if (options != NULL)
    *options = g_strdup (segment->options);
}

static void
_ncm_mset_catalog_storage_memory_set_sampler (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options, GError **error)
{
  NcmMSetCatalogStorageMemory *storage_mem = NCM_MSET_CATALOG_STORAGE_MEMORY (storage);
  NcmMSetCatalogStorageMemorySegment segment;

  g_assert_cmpuint (i, <=, storage_mem->samplers->len);

  segment.start   = start;
  segment.sampler = g_strdup (sampler);
  segment.options = g_strdup (options);

  if (i == storage_mem->samplers->len)
  {
    g_array_append_val (storage_mem->samplers, segment);
  }
  else
  {
    _ncm_mset_catalog_storage_memory_segment_clear (&g_array_index (storage_mem->samplers, NcmMSetCatalogStorageMemorySegment, i));
    g_array_index (storage_mem->samplers, NcmMSetCatalogStorageMemorySegment, i) = segment;
  }
}

/**
 * ncm_mset_catalog_storage_memory_new:
 *
 * Creates an uninitialized memory storage; ncm_mset_catalog_storage_set_desc() gives
 * it a description.
 *
 * Returns: (transfer full): a new #NcmMSetCatalogStorageMemory
 */
NcmMSetCatalogStorageMemory *
ncm_mset_catalog_storage_memory_new (void)
{
  return g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_MEMORY,
                       NULL);
}

/**
 * ncm_mset_catalog_storage_memory_ref:
 * @storage_mem: a #NcmMSetCatalogStorageMemory
 *
 * Increases the reference count of @storage_mem by one.
 *
 * Returns: (transfer full): @storage_mem
 */
NcmMSetCatalogStorageMemory *
ncm_mset_catalog_storage_memory_ref (NcmMSetCatalogStorageMemory *storage_mem)
{
  return g_object_ref (storage_mem);
}

/**
 * ncm_mset_catalog_storage_memory_free:
 * @storage_mem: a #NcmMSetCatalogStorageMemory
 *
 * Decreases the reference count of @storage_mem by one.
 */
void
ncm_mset_catalog_storage_memory_free (NcmMSetCatalogStorageMemory *storage_mem)
{
  g_object_unref (storage_mem);
}

/**
 * ncm_mset_catalog_storage_memory_clear:
 * @storage_mem: a #NcmMSetCatalogStorageMemory
 *
 * If *@storage_mem is not %NULL, decreases its reference count and sets *@storage_mem
 * to %NULL.
 */
void
ncm_mset_catalog_storage_memory_clear (NcmMSetCatalogStorageMemory **storage_mem)
{
  g_clear_object (storage_mem);
}


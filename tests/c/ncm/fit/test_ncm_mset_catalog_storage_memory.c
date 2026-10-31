/***************************************************************************
 *            test_ncm_mset_catalog_storage_memory.c
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mset_catalog_storage_memory.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include <numcosmo/numcosmo.h>

#define TEST_STORAGE_DIM 3
#define TEST_STORAGE_NCHAINS 2
#define TEST_STORAGE_RUN_TYPE "test-run"
#define TEST_STORAGE_NADD 1

typedef struct _TestNcmMSetCatalogStorageMemory
{
  NcmModelMVND *model_mvnd;
  NcmMSet *mset;
  NcmMSetCatalogDesc *desc;
  NcmMSetCatalogStorage *storage;
} TestNcmMSetCatalogStorageMemory;

static const gchar *test_names[]   = {"m2lnL", NULL};
static const gchar *test_symbols[] = {"-2\\ln(L)", NULL};

void test_ncm_mset_catalog_storage_memory_new (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_free (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);

void test_ncm_mset_catalog_storage_memory_uninitialized (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_set_desc (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_recording (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_samplers (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_rng (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_rng_after_rows (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_markovian_id (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_bind_order (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_serialize (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata);
void test_ncm_mset_catalog_storage_memory_property_traps (void);
void test_ncm_mset_catalog_storage_memory_property_subprocess (void);

static NcmVector *_test_ncm_mset_catalog_storage_memory_row (guint i);
static void _test_ncm_mset_catalog_storage_memory_add_rows (NcmMSetCatalogStorage *storage, guint first, guint n);
static void _test_ncm_mset_catalog_storage_memory_assert_row (NcmMSetCatalogStorage *storage, guint i);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/mset/catalog/storage/memory/uninitialized", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_uninitialized,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/set_desc", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_set_desc,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/recording", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_recording,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/samplers", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_samplers,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/rng", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_rng,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/rng_after_rows", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_rng_after_rows,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/markovian_id", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_markovian_id,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/bind_order", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_bind_order,
              &test_ncm_mset_catalog_storage_memory_free);
  g_test_add ("/ncm/mset/catalog/storage/memory/serialize", TestNcmMSetCatalogStorageMemory, NULL,
              &test_ncm_mset_catalog_storage_memory_new,
              &test_ncm_mset_catalog_storage_memory_serialize,
              &test_ncm_mset_catalog_storage_memory_free);

  g_test_add_func ("/ncm/mset/catalog/storage/memory/property/traps", &test_ncm_mset_catalog_storage_memory_property_traps);
  g_test_add_func ("/ncm/mset/catalog/storage/memory/property/subprocess", &test_ncm_mset_catalog_storage_memory_property_subprocess);

  g_test_run ();
}

/* Row i: column c holds 100 i + c, so every value identifies its row and column. */
static NcmVector *
_test_ncm_mset_catalog_storage_memory_row (guint i)
{
  NcmVector *row = ncm_vector_new (TEST_STORAGE_NADD + TEST_STORAGE_DIM);
  guint c;

  for (c = 0; c < ncm_vector_len (row); c++)
    ncm_vector_set (row, c, 100.0 * i + c);

  return row;
}

static void
_test_ncm_mset_catalog_storage_memory_add_rows (NcmMSetCatalogStorage *storage, guint first, guint n)
{
  guint i;

  for (i = first; i < first + n; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_memory_row (i);
    GError *error  = NULL;

    ncm_mset_catalog_storage_add_row (storage, row, &error);
    g_assert_no_error (error);

    ncm_vector_free (row);
  }
}

static void
_test_ncm_mset_catalog_storage_memory_assert_row (NcmMSetCatalogStorage *storage, guint i)
{
  GError *error     = NULL;
  NcmVector *row    = ncm_mset_catalog_storage_get_row (storage, i, &error);
  NcmVector *expect = _test_ncm_mset_catalog_storage_memory_row (i);
  guint c;

  g_assert_no_error (error);

  for (c = 0; c < ncm_vector_len (expect); c++)
    g_assert_cmpfloat (ncm_vector_get (row, c), ==, ncm_vector_get (expect, c));

  ncm_vector_free (row);
  ncm_vector_free (expect);
}

void
test_ncm_mset_catalog_storage_memory_new (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  GError *error = NULL;

  test->model_mvnd = ncm_model_mvnd_new (TEST_STORAGE_DIM);
  test->mset       = ncm_mset_new (NCM_MODEL (test->model_mvnd), NULL, NULL);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  test->desc = ncm_mset_catalog_desc_new (test->mset, TEST_STORAGE_NCHAINS, FALSE, test_names, test_symbols,
                                          0, TEST_STORAGE_RUN_TYPE, &error);
  g_assert_no_error (error);

  test->storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_memory_new ());
}

void
test_ncm_mset_catalog_storage_memory_free (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_mset_catalog_storage_free, test->storage);
  ncm_mset_catalog_desc_clear (&test->desc);
  ncm_mset_clear (&test->mset);
  ncm_model_mvnd_clear (&test->model_mvnd);
}

void
test_ncm_mset_catalog_storage_memory_uninitialized (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NcmVector *row = _test_ncm_mset_catalog_storage_memory_row (0);
  GError *error  = NULL;

  g_assert_null (ncm_mset_catalog_storage_peek_desc (test->storage, &error));
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (test->storage, &error), ==, 0);
  g_assert_false (ncm_mset_catalog_storage_is_recording (test->storage));

  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", NULL, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED);
  g_clear_error (&error);

  ncm_mset_catalog_storage_add_row (test->storage, row, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING);
  g_clear_error (&error);

  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING);
  g_clear_error (&error);

  {
    gchar *algo  = NULL;
    gchar *state = NULL;

    ncm_mset_catalog_storage_get_rng_init (test->storage, &algo, NULL, &state, &error);
    g_assert_null (algo);
    g_assert_null (state);
  }
  g_assert_null (ncm_mset_catalog_storage_get_rng (test->storage, &error));
  g_assert_no_error (error);

  {
    NcmRNG *rng = ncm_rng_seeded_new (NULL, 1234);

    ncm_mset_catalog_storage_set_rng (test->storage, rng, &error);
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED);
    g_clear_error (&error);
    ncm_rng_free (rng);
  }

  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (test->storage, &error), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (test->storage, &error), ==, 0);
  g_assert_no_error (error);

  ncm_vector_free (row);
}

void
test_ncm_mset_catalog_storage_memory_set_desc (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  GError *error = NULL;
  NcmMSetCatalogDesc *other;

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset_catalog_storage_peek_desc (test->storage, &error) == test->desc);

  /* A compatible description is bound and returned from then on. */
  other = ncm_mset_catalog_desc_new (test->mset, TEST_STORAGE_NCHAINS, FALSE, test_names, test_symbols,
                                     0, TEST_STORAGE_RUN_TYPE, &error);
  ncm_mset_catalog_storage_set_desc (test->storage, other, &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset_catalog_storage_peek_desc (test->storage, &error) == other);
  ncm_mset_catalog_desc_free (other);

  other = ncm_mset_catalog_desc_new (test->mset, TEST_STORAGE_NCHAINS + 1, FALSE, test_names, test_symbols,
                                     0, TEST_STORAGE_RUN_TYPE, &error);
  ncm_mset_catalog_storage_set_desc (test->storage, other, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_INCOMPATIBLE);
  g_clear_error (&error);

  /* No description change while a recording is open. */
  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", NULL, &error);
  g_assert_no_error (error);
  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING);
  g_clear_error (&error);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_desc_free (other);
}

void
test_ncm_mset_catalog_storage_memory_recording (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NcmVector *short_row = ncm_vector_new (TEST_STORAGE_DIM);
  GError *error        = NULL;
  guint i;

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", "opts", &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset_catalog_storage_is_recording (test->storage));

  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", NULL, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING);
  g_clear_error (&error);

  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 0, 5);

  ncm_mset_catalog_storage_add_row (test->storage, short_row, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW);
  g_clear_error (&error);

  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_no_error (error);
  g_assert_false (ncm_mset_catalog_storage_is_recording (test->storage));

  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (test->storage, &error), ==, 5);

  for (i = 0; i < 5; i++)
    _test_ncm_mset_catalog_storage_memory_assert_row (test->storage, i);

  g_assert_null (ncm_mset_catalog_storage_get_row (test->storage, 5, &error));
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW);
  g_clear_error (&error);

  ncm_vector_free (short_row);
}

void
test_ncm_mset_catalog_storage_memory_samplers (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  GError *error = NULL;
  gchar *sampler, *options;
  guint start;

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);

  ncm_mset_catalog_storage_start_recording (test->storage, "A", "a-opts", &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 0, 3);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);

  /* A recording without rows leaves a segment that the next recording replaces. */
  ncm_mset_catalog_storage_start_recording (test->storage, "B", NULL, &error);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (test->storage, &error), ==, 2);

  ncm_mset_catalog_storage_start_recording (test->storage, "C", "c-opts", &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 3, 2);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_no_error (error);

  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (test->storage, &error), ==, 2);

  ncm_mset_catalog_storage_get_sampler (test->storage, 0, &start, &sampler, &options, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (start, ==, 0);
  g_assert_cmpstr (sampler, ==, "A");
  g_assert_cmpstr (options, ==, "a-opts");
  g_free (sampler);
  g_free (options);

  ncm_mset_catalog_storage_get_sampler (test->storage, 1, &start, &sampler, &options, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (start, ==, 3);
  g_assert_cmpstr (sampler, ==, "C");
  g_assert_cmpstr (options, ==, "c-opts");
  g_free (sampler);
  g_free (options);

  ncm_mset_catalog_storage_get_sampler (test->storage, 2, NULL, NULL, NULL, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW);
  g_clear_error (&error);
}

void
test_ncm_mset_catalog_storage_memory_rng (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NcmRNG *rng   = ncm_rng_seeded_new ("mt19937", 1234);
  GError *error = NULL;
  gchar *state0 = ncm_rng_get_state (rng);
  NcmRNG *stored;
  gchar *state1, *state;

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);

  /* The RNG's state becomes the initial and the current state; the storage keeps the RNG. */
  ncm_mset_catalog_storage_set_rng (test->storage, rng, &error);
  g_assert_no_error (error);

  {
    gchar *algo = NULL;
    gulong seed = 0;

    ncm_mset_catalog_storage_get_rng_init (test->storage, &algo, &seed, &state, &error);
    g_assert_no_error (error);
    g_assert_cmpstr (algo, ==, "mt19937");
    g_assert_cmpuint (seed, ==, 1234);
    g_assert_cmpstr (state, ==, state0);
    g_free (algo);
    g_free (state);
  }

  stored = ncm_mset_catalog_storage_get_rng (test->storage, &error);
  g_assert_no_error (error);
  g_assert_true (stored == rng);
  ncm_rng_free (stored);

  /* A second RNG is refused, with or without rows. */
  ncm_mset_catalog_storage_set_rng (test->storage, rng, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG);
  g_clear_error (&error);

  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", NULL, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_set_rng (test->storage, rng, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING);
  g_clear_error (&error);

  ncm_rng_uniform_gen (rng, 0.0, 1.0);
  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 0, 2);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_no_error (error);

  /* The stored state is the state when the recording ended. */
  state1 = ncm_rng_get_state (rng);
  state  = ncm_mset_catalog_storage_get_rng_stat (test->storage, &error);
  g_assert_cmpstr (state, ==, state1);
  g_assert_cmpstr (state1, !=, state0);
  g_free (state);

  /* A storage rebuilt from its properties builds the RNG from the stored algorithm, seed and
   * current state. */
  {
    NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcmMSetCatalogStorage *dup = NCM_MSET_CATALOG_STORAGE (ncm_serialize_dup_obj (ser, G_OBJECT (test->storage)));
    NcmRNG *rng_dup            = ncm_mset_catalog_storage_get_rng (dup, &error);
    gchar *state_dup;

    g_assert_no_error (error);
    g_assert_nonnull (rng_dup);
    g_assert_cmpstr (ncm_rng_get_algo (rng_dup), ==, "mt19937");
    g_assert_cmpuint (ncm_rng_get_seed (rng_dup), ==, 1234);

    state_dup = ncm_rng_get_state (rng_dup);
    g_assert_cmpstr (state_dup, ==, state1);
    g_free (state_dup);

    ncm_rng_free (rng_dup);
    ncm_mset_catalog_storage_free (dup);
    ncm_serialize_free (ser);
  }

  g_free (state0);
  g_free (state1);
  ncm_rng_free (rng);
}

/* A storage with rows and no RNG never gets one. */
void
test_ncm_mset_catalog_storage_memory_rng_after_rows (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, 1234);
  GError *error = NULL;

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);
  ncm_mset_catalog_storage_start_recording (test->storage, "join", NULL, &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 0, 2);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_set_rng (test->storage, rng, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG);
  g_clear_error (&error);

  g_assert_null (ncm_mset_catalog_storage_get_rng (test->storage, &error));
  g_assert_no_error (error);

  ncm_rng_free (rng);
}

void
test_ncm_mset_catalog_storage_memory_markovian_id (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  GError *error = NULL;

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);

  ncm_mset_catalog_storage_set_markovian_id (test->storage, 0, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_NOT_RECORDING);
  g_clear_error (&error);

  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", NULL, &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 0, 4);

  ncm_mset_catalog_storage_set_markovian_id (test->storage, 2, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (test->storage, &error), ==, 2);

  /* The Markovian id may equal the number of rows: the chain starts at the next row. */
  ncm_mset_catalog_storage_set_markovian_id (test->storage, 4, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_set_markovian_id (test->storage, 3, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID);
  g_clear_error (&error);

  ncm_mset_catalog_storage_set_markovian_id (test->storage, 5, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_MARKOVIAN_ID);
  g_clear_error (&error);

  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (test->storage, &error), ==, 4);
}

void
test_ncm_mset_catalog_storage_memory_bind_order (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_STORAGE_DIM);
  NcmMSet *mset_rev        = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  GError *error            = NULL;
  NcmMSetCatalogDesc *desc_rev;
  gchar **fmap;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset_rev, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset_rev);

  fmap = ncm_mset_get_fmap (mset_rev);

  for (i = 0; i < TEST_STORAGE_DIM / 2; i++)
  {
    gchar *tmp = fmap[i];

    fmap[i]                        = fmap[TEST_STORAGE_DIM - 1 - i];
    fmap[TEST_STORAGE_DIM - 1 - i] = tmp;
  }

  ncm_mset_set_fmap (mset_rev, (const gchar * const *) fmap, FALSE, &error);
  g_assert_no_error (error);

  desc_rev = ncm_mset_catalog_desc_new (mset_rev, TEST_STORAGE_NCHAINS, FALSE, test_names, test_symbols,
                                        0, TEST_STORAGE_RUN_TYPE, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_set_desc (test->storage, test->desc, &error);
  ncm_mset_catalog_storage_start_recording (test->storage, "sampler", NULL, &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (test->storage, 0, 3);
  ncm_mset_catalog_storage_end_recording (test->storage, &error);
  g_assert_no_error (error);

  /* Binding the reversed description reverses the parameter columns of every row. */
  ncm_mset_catalog_storage_set_desc (test->storage, desc_rev, &error);
  g_assert_no_error (error);

  for (i = 0; i < 3; i++)
  {
    NcmVector *row    = ncm_mset_catalog_storage_get_row (test->storage, i, &error);
    NcmVector *expect = _test_ncm_mset_catalog_storage_memory_row (i);

    g_assert_cmpfloat (ncm_vector_get (row, 0), ==, ncm_vector_get (expect, 0));

    for (j = 0; j < TEST_STORAGE_DIM; j++)
      g_assert_cmpfloat (ncm_vector_get (row, TEST_STORAGE_NADD + j), ==,
                         ncm_vector_get (expect, TEST_STORAGE_NADD + TEST_STORAGE_DIM - 1 - j));

    ncm_vector_free (row);
    ncm_vector_free (expect);
  }

  ncm_mset_catalog_desc_free (desc_rev);
  g_strfreev (fmap);
  ncm_mset_free (mset_rev);
  ncm_model_mvnd_free (model_mvnd);
}

void
test_ncm_mset_catalog_storage_memory_serialize (TestNcmMSetCatalogStorageMemory *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 4321);
  GError *error     = NULL;
  NcmMSetCatalogStorage *storage, *dup;
  gchar *str, *str_dup, *sampler;
  guint i, start;

  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_memory_new ());

  ncm_mset_catalog_storage_set_desc (storage, test->desc, &error);
  ncm_mset_catalog_storage_set_rng (storage, rng, &error);
  ncm_mset_catalog_storage_start_recording (storage, "A", "a-opts", &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (storage, 0, 4);
  ncm_mset_catalog_storage_set_markovian_id (storage, 2, &error);
  ncm_mset_catalog_storage_end_recording (storage, &error);
  ncm_rng_uniform_gen (rng, 0.0, 1.0);
  ncm_mset_catalog_storage_start_recording (storage, "B", "b-opts", &error);
  _test_ncm_mset_catalog_storage_memory_add_rows (storage, 4, 2);
  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);

  dup = NCM_MSET_CATALOG_STORAGE (ncm_serialize_dup_obj (ser, G_OBJECT (storage)));

  g_assert_true (ncm_mset_catalog_desc_has_same_order (ncm_mset_catalog_storage_peek_desc (dup, &error), test->desc));
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (dup, &error), ==, 6);

  for (i = 0; i < 6; i++)
    _test_ncm_mset_catalog_storage_memory_assert_row (dup, i);

  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (dup, &error), ==, 2);

  ncm_mset_catalog_storage_get_rng_init (storage, NULL, NULL, &str, &error);
  ncm_mset_catalog_storage_get_rng_init (dup, NULL, NULL, &str_dup, &error);
  g_assert_cmpstr (str, ==, str_dup);
  g_free (str);
  g_free (str_dup);

  str     = ncm_mset_catalog_storage_get_rng_stat (storage, &error);
  str_dup = ncm_mset_catalog_storage_get_rng_stat (dup, &error);
  g_assert_cmpstr (str, ==, str_dup);
  g_free (str);
  g_free (str_dup);

  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (dup, &error), ==, 2);
  ncm_mset_catalog_storage_get_sampler (dup, 1, &start, &sampler, NULL, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (start, ==, 4);
  g_assert_cmpstr (sampler, ==, "B");
  g_free (sampler);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_free (dup);
  ncm_mset_catalog_storage_free (storage);
  ncm_rng_free (rng);
  ncm_serialize_free (ser);
}

void
test_ncm_mset_catalog_storage_memory_property_traps (void)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/storage/memory/property/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*1 rows without a description*");
}

void
test_ncm_mset_catalog_storage_memory_property_subprocess (void)
{
  /* Built from properties there is no error location, so inconsistent properties abort. */
  NcmObjArray *rows = ncm_obj_array_new ();
  NcmVector *row    = ncm_vector_new (2);

  ncm_obj_array_add (rows, G_OBJECT (row));

  g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_MEMORY, "rows", rows, NULL);

  g_assert_not_reached ();
}


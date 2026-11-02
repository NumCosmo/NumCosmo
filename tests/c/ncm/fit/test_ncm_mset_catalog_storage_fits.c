/***************************************************************************
 *            test_ncm_mset_catalog_storage_fits.c
 *
 *  Sat Oct 31 15:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mset_catalog_storage_fits.c
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

/*
 * Reads the fixtures of data/truth_tables/catalog/ through NcmMSetCatalogStorageFits and
 * checks the description, the run metadata and the rows (against the previous reader),
 * and writes, reopens and migrates files in a temporary directory.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include <numcosmo/numcosmo.h>

#ifdef HAVE_CFITSIO
#include <fitsio.h>
#endif /* HAVE_CFITSIO */

#define TEST_FITS_ESMCMC_NROWS 18
#define TEST_FITS_ESMCMC_NWALKERS 6
#define TEST_FITS_RNG_SEED 20261031
#define TEST_FITS_FIRST_ID 6

typedef struct _TestNcmMSetCatalogStorageFits
{
  gchar *tmpdir;
  gchar *filename;
  NcmMSetCatalogDesc *desc;
} TestNcmMSetCatalogStorageFits;

static void test_ncm_mset_catalog_storage_fits_new (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_free (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);

static void test_ncm_mset_catalog_storage_fits_fixture_esmcmc (void);
static void test_ncm_mset_catalog_storage_fits_fixture_esmcmc_variants (void);
static void test_ncm_mset_catalog_storage_fits_fixture_producers (void);
static void test_ncm_mset_catalog_storage_fits_fixture_weighted (void);
static void test_ncm_mset_catalog_storage_fits_fixture_rows (void);
static void test_ncm_mset_catalog_storage_fits_fixture_peek_row (void);
static void test_ncm_mset_catalog_storage_fits_fixture_read_only (void);
static void test_ncm_mset_catalog_storage_fits_errors (void);
static void test_ncm_mset_catalog_storage_fits_write_reopen (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_bind_order (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_migrate (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_erase (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_erase_legacy (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_serialize (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);
static void test_ncm_mset_catalog_storage_fits_float_column (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata);

static gchar *_test_ncm_mset_catalog_storage_fits_fixture_path (const gchar *name);
static NcmMSetCatalogStorage *_test_ncm_mset_catalog_storage_fits_open_fixture (const gchar *name);
static void _test_ncm_mset_catalog_storage_fits_check_segment (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options);
static NcmVector *_test_ncm_mset_catalog_storage_fits_row (guint i, guint ncols);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

#ifdef HAVE_CFITSIO
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/esmcmc", &test_ncm_mset_catalog_storage_fits_fixture_esmcmc);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/esmcmc_variants", &test_ncm_mset_catalog_storage_fits_fixture_esmcmc_variants);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/producers", &test_ncm_mset_catalog_storage_fits_fixture_producers);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/weighted", &test_ncm_mset_catalog_storage_fits_fixture_weighted);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/rows", &test_ncm_mset_catalog_storage_fits_fixture_rows);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/peek_row", &test_ncm_mset_catalog_storage_fits_fixture_peek_row);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/fixture/read_only", &test_ncm_mset_catalog_storage_fits_fixture_read_only);
  g_test_add_func ("/ncm/mset/catalog/storage/fits/errors", &test_ncm_mset_catalog_storage_fits_errors);
  g_test_add ("/ncm/mset/catalog/storage/fits/write_reopen", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_write_reopen,
              &test_ncm_mset_catalog_storage_fits_free);
  g_test_add ("/ncm/mset/catalog/storage/fits/bind_order", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_bind_order,
              &test_ncm_mset_catalog_storage_fits_free);
  g_test_add ("/ncm/mset/catalog/storage/fits/migrate", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_migrate,
              &test_ncm_mset_catalog_storage_fits_free);
  g_test_add ("/ncm/mset/catalog/storage/fits/erase", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_erase,
              &test_ncm_mset_catalog_storage_fits_free);
  g_test_add ("/ncm/mset/catalog/storage/fits/erase_legacy", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_erase_legacy,
              &test_ncm_mset_catalog_storage_fits_free);
  g_test_add ("/ncm/mset/catalog/storage/fits/serialize", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_serialize,
              &test_ncm_mset_catalog_storage_fits_free);
  g_test_add ("/ncm/mset/catalog/storage/fits/float_column", TestNcmMSetCatalogStorageFits, NULL,
              &test_ncm_mset_catalog_storage_fits_new,
              &test_ncm_mset_catalog_storage_fits_float_column,
              &test_ncm_mset_catalog_storage_fits_free);
#endif /* HAVE_CFITSIO */

  g_test_run ();
}

/*
 * A temporary directory, a file name in it and a Desc with two chains and two free
 * parameters of an NcmModelMVND.
 */
static void
test_ncm_mset_catalog_storage_fits_new (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  NcmModelMVND *model    = ncm_model_mvnd_new (2);
  NcmMSet *mset          = ncm_mset_new (NCM_MODEL (model), NULL, NULL);
  const gchar *names[]   = { "m2lnL", NULL };
  const gchar *symbols[] = { "-2\\ln(L)", NULL };
  GError *error          = NULL;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  test->tmpdir   = g_dir_make_tmp ("test_ncm_mset_catalog_storage_fits_XXXXXX", &error);
  test->filename = g_build_filename (test->tmpdir, "catalog.fits", NULL);
  test->desc     = ncm_mset_catalog_desc_new (mset, 2, FALSE, names, symbols, 0, NCM_MSET_CATALOG_RTYPE_ESMCMC, &error);
  g_assert_no_error (error);

  ncm_mset_free (mset);
  ncm_model_mvnd_free (model);
}

static void
test_ncm_mset_catalog_storage_fits_free (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  GDir *dir = g_dir_open (test->tmpdir, 0, NULL);
  const gchar *name;

  while ((name = g_dir_read_name (dir)) != NULL)
  {
    gchar *path = g_build_filename (test->tmpdir, name, NULL);

    g_unlink (path);
    g_free (path);
  }

  g_dir_close (dir);
  g_rmdir (test->tmpdir);

  ncm_mset_catalog_desc_free (test->desc);
  g_free (test->filename);
  g_free (test->tmpdir);
}

static gchar *
_test_ncm_mset_catalog_storage_fits_fixture_path (const gchar *name)
{
  gchar *rel_path = g_strdup_printf ("truth_tables/catalog/%s", name);
  gchar *path     = ncm_cfg_get_data_filename (rel_path, TRUE);

  g_free (rel_path);

  return path;
}

/* A catalog on a read-only FITS storage of filename. */
static NcmMSetCatalog *
_test_ncm_mset_catalog_storage_fits_open_catalog (const gchar *filename)
{
  NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (filename, TRUE, NULL);
  NcmMSetCatalog *mcat                    = ncm_mset_catalog_new_from_storage (NCM_MSET_CATALOG_STORAGE (storage_fits), 0, 1);

  ncm_mset_catalog_storage_fits_free (storage_fits);

  return mcat;
}

/* Opens data/truth_tables/catalog/name read-only. */
static NcmMSetCatalogStorage *
_test_ncm_mset_catalog_storage_fits_open_fixture (const gchar *name)
{
  gchar *path                             = _test_ncm_mset_catalog_storage_fits_fixture_path (name);
  GError *error                           = NULL;
  NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (path, TRUE, &error);

  g_assert_no_error (error);
  g_assert_nonnull (storage_fits);
  g_assert_true (ncm_mset_catalog_storage_get_read_only (NCM_MSET_CATALOG_STORAGE (storage_fits)));
  g_assert_cmpstr (ncm_mset_catalog_storage_fits_peek_filename (storage_fits), ==, path);

  g_free (path);

  return NCM_MSET_CATALOG_STORAGE (storage_fits);
}

static void
_test_ncm_mset_catalog_storage_fits_check_segment (NcmMSetCatalogStorage *storage, guint i, guint start, const gchar *sampler, const gchar *options)
{
  GError *error    = NULL;
  guint start_i    = 0;
  gchar *sampler_i = NULL;
  gchar *options_i = NULL;

  ncm_mset_catalog_storage_get_sampler (storage, i, &start_i, &sampler_i, &options_i, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (start_i, ==, start);
  g_assert_cmpstr (sampler_i, ==, sampler);
  g_assert_cmpstr (options_i, ==, options);

  g_free (sampler_i);
  g_free (options_i);
}

/* Row i: column c holds 100 i + c. */
static NcmVector *
_test_ncm_mset_catalog_storage_fits_row (guint i, guint ncols)
{
  NcmVector *row = ncm_vector_new (ncols);
  guint c;

  for (c = 0; c < ncols; c++)
    ncm_vector_set (row, c, 100.0 * i + c);

  return row;
}

static void
test_ncm_mset_catalog_storage_fits_fixture_esmcmc (void)
{
  NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc.fits");
  GError *error                  = NULL;
  NcmMSetCatalogDesc *desc       = ncm_mset_catalog_storage_peek_desc (storage, &error);
  gchar *algo                    = NULL;
  gchar *inis                    = NULL;
  gchar *stat                    = NULL;
  gulong seed                    = 0;
  NcmRNG *rng;

  g_assert_no_error (error);
  g_assert_nonnull (desc);
  g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (desc), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);
  g_assert_cmpuint (ncm_mset_catalog_desc_nchains (desc), ==, TEST_FITS_ESMCMC_NWALKERS);
  g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (desc), ==, 1);
  g_assert_cmpuint (ncm_mset_catalog_desc_ncols (desc), ==, 3);
  g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 0), ==, NCM_MSET_CATALOG_M2LNL_COLNAME);
  g_assert_cmpstr (ncm_mset_catalog_desc_col_symb (desc, 0), ==, NCM_MSET_CATALOG_M2LNL_SYMBOL);
  g_assert_cmpstr (ncm_mset_catalog_desc_col_full_name (desc, 1), ==, "NcmModelMVND:mu_0");
  g_assert_cmpstr (ncm_mset_catalog_desc_col_full_name (desc, 2), ==, "NcmModelMVND:mu_1");
  g_assert_false (ncm_mset_catalog_desc_weighted (desc));
  g_assert_cmpint (ncm_mset_catalog_desc_get_m2lnp_var (desc), ==, 0);

  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, TEST_FITS_ESMCMC_NROWS);
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, TEST_FITS_ESMCMC_NWALKERS);

  /* The translation of INITSMP, SAMPLER and SAMPOPT. */
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 2);
  _test_ncm_mset_catalog_storage_fits_check_segment (storage, 0, 0, "Multivariate Gaussian Sampler", "");
  _test_ncm_mset_catalog_storage_fits_check_segment (storage, 1, TEST_FITS_ESMCMC_NWALKERS, "Stretch-Move", "multi-stretch=no");

  ncm_mset_catalog_storage_get_rng_init (storage, &algo, &seed, &inis, &error);
  g_assert_no_error (error);
  g_assert_cmpstr (algo, ==, "mt19937");
  g_assert_cmpuint (seed, ==, TEST_FITS_RNG_SEED);
  g_assert_nonnull (inis);

  stat = ncm_mset_catalog_storage_get_rng_stat (storage, &error);
  g_assert_nonnull (stat);
  g_assert_cmpstr (stat, !=, inis);

  rng = ncm_mset_catalog_storage_get_rng (storage, &error);
  g_assert_no_error (error);
  g_assert_nonnull (rng);
  g_assert_cmpstr (ncm_rng_get_algo (rng), ==, "mt19937");
  g_assert_cmpuint (ncm_rng_get_seed (rng), ==, TEST_FITS_RNG_SEED);
  {
    gchar *rng_state = ncm_rng_get_state (rng);

    g_assert_cmpstr (rng_state, ==, stat);
    g_free (rng_state);
  }

  g_assert_false (ncm_mset_catalog_storage_is_recording (storage));

  ncm_rng_free (rng);
  g_free (algo);
  g_free (inis);
  g_free (stat);
  ncm_mset_catalog_storage_free (storage);
}

static void
test_ncm_mset_catalog_storage_fits_fixture_esmcmc_variants (void)
{
  GError *error = NULL;

  {
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_no_rng.fits");
    gchar *algo                    = NULL;

    ncm_mset_catalog_storage_get_rng_init (storage, &algo, NULL, NULL, &error);
    g_assert_no_error (error);
    g_assert_null (algo);
    g_assert_null (ncm_mset_catalog_storage_get_rng (storage, &error));
    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, TEST_FITS_ESMCMC_NROWS);

    ncm_mset_catalog_storage_free (storage);
  }

  {
    /* RNG state without RNG_ALGO: the default algorithm. */
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_no_rng_algo.fits");
    NcmRNG *rng_default            = ncm_rng_new (NULL);
    gchar *algo                    = NULL;
    gulong seed                    = 0;
    NcmRNG *rng;

    ncm_mset_catalog_storage_get_rng_init (storage, &algo, &seed, NULL, &error);
    g_assert_no_error (error);
    g_assert_cmpstr (algo, ==, ncm_rng_get_algo (rng_default));
    g_assert_cmpuint (seed, ==, TEST_FITS_RNG_SEED);

    rng = ncm_mset_catalog_storage_get_rng (storage, &error);
    g_assert_no_error (error);
    g_assert_nonnull (rng);

    ncm_rng_free (rng);
    ncm_rng_free (rng_default);
    g_free (algo);
    ncm_mset_catalog_storage_free (storage);
  }

  {
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_no_m2lnp_id.fits");
    NcmMSetCatalogDesc *desc       = ncm_mset_catalog_storage_peek_desc (storage, &error);

    g_assert_no_error (error);
    g_assert_cmpint (ncm_mset_catalog_desc_get_m2lnp_var (desc), ==, -1);

    ncm_mset_catalog_storage_free (storage);
  }

  {
    /* MARKID is in the file numbering, from FIRST_ID. */
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_first_id.fits");

    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, TEST_FITS_ESMCMC_NROWS);
    g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, TEST_FITS_ESMCMC_NWALKERS);
    g_assert_no_error (error);

    ncm_mset_catalog_storage_free (storage);
  }

  {
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_no_markid.fits");

    g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 0);
    g_assert_no_error (error);

    ncm_mset_catalog_storage_free (storage);
  }

  {
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_header_only.fits");

    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 0);
    g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 0);
    g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 2);
    _test_ncm_mset_catalog_storage_fits_check_segment (storage, 1, 0, "Stretch-Move", "multi-stretch=no");
    g_assert_no_error (error);

    ncm_mset_catalog_storage_free (storage);
  }
}

static void
test_ncm_mset_catalog_storage_fits_fixture_producers (void)
{
  const struct
  {
    const gchar *file;
    const gchar *rtype;
    const gchar *segment;
    guint nrows;
  } cases[] = {
    { "mc.fits", NCM_MSET_CATALOG_RTYPE_MC, "from-model", 10 },
    { "mcbs.fits", NCM_MSET_CATALOG_RTYPE_MCBS, "bootstrap-mean", 5 },
    { "mcmc.fits", NCM_MSET_CATALOG_RTYPE_MCMC, "Multivariate Gaussian Sampler", 10 },
  };

  GError *error = NULL;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture (cases[i].file);
    NcmMSetCatalogDesc *desc       = ncm_mset_catalog_storage_peek_desc (storage, &error);
    NcmRNG *rng;

    g_assert_no_error (error);
    g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (desc), ==, cases[i].rtype);
    g_assert_cmpuint (ncm_mset_catalog_desc_nchains (desc), ==, 1);
    g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (desc), ==, 1);
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, cases[i].nrows);
    g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 1);
    _test_ncm_mset_catalog_storage_fits_check_segment (storage, 0, 0, cases[i].segment, "");

    rng = ncm_mset_catalog_storage_get_rng (storage, &error);
    g_assert_no_error (error);
    g_assert_nonnull (rng);

    ncm_rng_free (rng);
    ncm_mset_catalog_storage_free (storage);
  }
}

static void
test_ncm_mset_catalog_storage_fits_fixture_weighted (void)
{
  GError *error = NULL;

  {
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("weighted.fits");
    NcmMSetCatalogDesc *desc       = ncm_mset_catalog_storage_peek_desc (storage, &error);
    NcmVector *row;

    g_assert_no_error (error);
    g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (desc), ==, "Fixture");
    g_assert_true (ncm_mset_catalog_desc_weighted (desc));
    g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (desc), ==, 2);
    g_assert_cmpuint (ncm_mset_catalog_desc_ncols (desc), ==, 4);
    g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 0), ==, "m2lnL");
    g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 1), ==, NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME);
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 10);
    g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 0);

    /* Row i has weight 1 / (i + 1). */
    row = ncm_mset_catalog_storage_get_row (storage, 3, &error);
    g_assert_no_error (error);
    g_assert_cmpfloat (ncm_vector_get (row, 1), ==, 0.25);
    ncm_vector_free (row);

    ncm_mset_catalog_storage_free (storage);
  }

  {
    /* The functions key in HDU0 is ignored. */
    NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("functions.fits");
    NcmMSetCatalogDesc *desc       = ncm_mset_catalog_storage_peek_desc (storage, &error);

    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_desc_nadd_vals (desc), ==, 2);
    g_assert_cmpstr (ncm_mset_catalog_desc_col_name (desc, 1), ==, "H0");
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 10);

    ncm_mset_catalog_storage_free (storage);
  }
}

/*
 * Every row of esmcmc.fits equals the row the previous reader gives, through the block
 * reads of a one-row block and of the default block.
 */
static void
test_ncm_mset_catalog_storage_fits_fixture_rows (void)
{
  gchar *path               = _test_ncm_mset_catalog_storage_fits_fixture_path ("esmcmc.fits");
  NcmMSetCatalog *mcat      = _test_ncm_mset_catalog_storage_fits_open_catalog (path);
  const guint read_blocks[] = { 1, 5, 1024 };
  GError *error             = NULL;
  guint b, i, c;

  for (b = 0; b < G_N_ELEMENTS (read_blocks); b++)
  {
    NcmMSetCatalogStorageFits *storage_fits = g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_FITS,
                                                            "filename", path,
                                                            "read-only", TRUE,
                                                            "read-block", read_blocks[b],
                                                            NULL);
    NcmMSetCatalogStorage *storage = NCM_MSET_CATALOG_STORAGE (storage_fits);

    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, ncm_mset_catalog_len (mcat));

    for (i = 0; i < TEST_FITS_ESMCMC_NROWS; i++)
    {
      const guint j     = (i * 7) % TEST_FITS_ESMCMC_NROWS;
      NcmVector *row    = ncm_mset_catalog_storage_get_row (storage, j, &error);
      NcmVector *mc_row = ncm_mset_catalog_peek_row (mcat, j);

      g_assert_no_error (error);
      g_assert_cmpuint (ncm_vector_len (row), ==, 3);

      for (c = 0; c < 3; c++)
        g_assert_cmpfloat (ncm_vector_get (row, c), ==, ncm_vector_get (mc_row, c));

      ncm_vector_free (row);
    }

    g_assert_null (ncm_mset_catalog_storage_get_row (storage, TEST_FITS_ESMCMC_NROWS, &error));
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW);
    g_clear_error (&error);

    ncm_mset_catalog_storage_free (storage);
  }

  ncm_mset_catalog_free (mcat);
  g_free (path);
}

/*
 * peek_row on the fixture for several block sizes: the values of the previous reader,
 * the same vector as get_row, and the pointer kept while the same row is peeked again.
 */
static void
test_ncm_mset_catalog_storage_fits_fixture_peek_row (void)
{
  gchar *path               = _test_ncm_mset_catalog_storage_fits_fixture_path ("esmcmc.fits");
  NcmMSetCatalog *mcat      = _test_ncm_mset_catalog_storage_fits_open_catalog (path);
  const guint read_blocks[] = { 1, 5, 1024 };
  GError *error             = NULL;
  guint b, i, c;

  for (b = 0; b < G_N_ELEMENTS (read_blocks); b++)
  {
    NcmMSetCatalogStorageFits *storage_fits = g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_FITS,
                                                            "filename", path,
                                                            "read-only", TRUE,
                                                            "read-block", read_blocks[b],
                                                            NULL);
    NcmMSetCatalogStorage *storage = NCM_MSET_CATALOG_STORAGE (storage_fits);

    for (i = 0; i < TEST_FITS_ESMCMC_NROWS; i++)
    {
      const guint j     = (i * 7) % TEST_FITS_ESMCMC_NROWS;
      NcmVector *row    = ncm_mset_catalog_storage_peek_row (storage, j, &error);
      NcmVector *mc_row = ncm_mset_catalog_peek_row (mcat, j);
      NcmVector *ref;

      g_assert_no_error (error);
      g_assert_nonnull (row);

      for (c = 0; c < 3; c++)
        g_assert_cmpfloat (ncm_vector_get (row, c), ==, ncm_vector_get (mc_row, c));

      g_assert_true (row == ncm_mset_catalog_storage_peek_row (storage, j, &error));

      ref = ncm_mset_catalog_storage_get_row (storage, j, &error);
      g_assert_true (ref == row);
      ncm_vector_free (ref);
    }

    g_assert_null (ncm_mset_catalog_storage_peek_row (storage, TEST_FITS_ESMCMC_NROWS, &error));
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW);
    g_clear_error (&error);

    ncm_mset_catalog_storage_free (storage);
  }

  ncm_mset_catalog_free (mcat);
  g_free (path);
}

static void
test_ncm_mset_catalog_storage_fits_fixture_read_only (void)
{
  NcmMSetCatalogStorage *storage = _test_ncm_mset_catalog_storage_fits_open_fixture ("esmcmc_no_rng.fits");
  GError *error                  = NULL;

  ncm_mset_catalog_storage_start_recording (storage, "A", NULL, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_READ_ONLY);
  g_clear_error (&error);
  g_assert_false (ncm_mset_catalog_storage_is_recording (storage));

  {
    NcmRNG *rng = ncm_rng_seeded_new (NULL, 1);

    ncm_mset_catalog_storage_set_rng (storage, rng, &error);
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_READ_ONLY);
    g_clear_error (&error);
    ncm_rng_free (rng);
  }

  ncm_mset_catalog_storage_free (storage);
}

static void
test_ncm_mset_catalog_storage_fits_errors (void)
{
  GError *error = NULL;
  NcmMSetCatalogStorageFits *storage_fits;

  storage_fits = ncm_mset_catalog_storage_fits_new ("/nonexistent/dir/catalog.fits", TRUE, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND);
  g_assert_null (storage_fits);
  g_clear_error (&error);

  {
    /* Not a catalog: the fixture of another subsystem. */
    gchar *path = ncm_cfg_get_data_filename ("truth_tables/enums/enum_nicks.json", TRUE);

    storage_fits = ncm_mset_catalog_storage_fits_new (path, TRUE, &error);
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND);
    g_assert_null (storage_fits);
    g_clear_error (&error);
    g_free (path);
  }

  {
    NcmMSetCatalogStorage *storage = NCM_MSET_CATALOG_STORAGE (g_object_new (NCM_TYPE_MSET_CATALOG_STORAGE_FITS, NULL));

    g_assert_null (ncm_mset_catalog_storage_peek_desc (storage, &error));
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_BACKEND);
    g_clear_error (&error);

    ncm_mset_catalog_storage_free (storage);
  }
}

/*
 * A new file: Desc, RNG, two recordings with syncs after every row, reopened writable
 * and read by the previous reader.
 */
static void
test_ncm_mset_catalog_storage_fits_write_reopen (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  GError *error                           = NULL;
  NcmRNG *rng                             = ncm_rng_seeded_new (NULL, 4321);
  NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error);
  NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);
  gchar *rng_state;
  guint i;

  g_assert_no_error (error);
  g_assert_false (g_file_test (test->filename, G_FILE_TEST_EXISTS));
  g_assert_null (ncm_mset_catalog_storage_peek_desc (storage, &error));
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 0);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_set_desc (storage, test->desc, &error);
  g_assert_no_error (error);
  g_assert_true (g_file_test (test->filename, G_FILE_TEST_EXISTS));
  g_assert_true (ncm_mset_catalog_storage_peek_desc (storage, &error) == test->desc);

  ncm_mset_catalog_storage_set_rng (storage, rng, &error);
  g_assert_no_error (error);

  g_object_set (storage, "smode", NCM_MSET_CATALOG_SYNC_AUTO, NULL);

  ncm_mset_catalog_storage_start_recording (storage, "A", "a-opts", &error);
  g_assert_no_error (error);

  for (i = 0; i < 8; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, 3);

    ncm_rng_uniform_gen (rng, 0.0, 1.0);
    ncm_mset_catalog_storage_add_row (storage, row, &error);
    g_assert_no_error (error);
    ncm_vector_free (row);
  }

  /* With keep-last 2 and two chains, rows before the last four come from the file. */
  for (i = 0; i < 8; i++)
  {
    NcmVector *row = ncm_mset_catalog_storage_get_row (storage, i, &error);

    g_assert_no_error (error);
    g_assert_cmpfloat (ncm_vector_get (row, 0), ==, 100.0 * i);
    g_assert_cmpfloat (ncm_vector_get (row, 2), ==, 100.0 * i + 2.0);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_set_markovian_id (storage, 2, &error);
  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_start_recording (storage, "B", NULL, &error);

  for (i = 8; i < 12; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, 3);

    ncm_rng_uniform_gen (rng, 0.0, 1.0);
    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);

  rng_state = ncm_rng_get_state (rng);
  ncm_mset_catalog_storage_free (storage);

  /* Reopened writable. */
  storage_fits = ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error);
  storage      = NCM_MSET_CATALOG_STORAGE (storage_fits);
  g_assert_no_error (error);

  {
    NcmMSetCatalogDesc *desc = ncm_mset_catalog_storage_peek_desc (storage, &error);
    gchar *algo              = NULL;
    gchar *stat              = NULL;
    gulong seed              = 0;
    NcmRNG *rng_file;

    g_assert_nonnull (desc);
    g_assert_true (ncm_mset_catalog_desc_is_compatible (desc, test->desc, NULL));
    g_assert_true (ncm_mset_catalog_desc_has_same_order (desc, test->desc));
    g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (desc), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);

    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 12);
    g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 2);
    g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 2);
    _test_ncm_mset_catalog_storage_fits_check_segment (storage, 0, 0, "A", "a-opts");
    _test_ncm_mset_catalog_storage_fits_check_segment (storage, 1, 8, "B", "");

    for (i = 0; i < 12; i++)
    {
      NcmVector *row = ncm_mset_catalog_storage_get_row (storage, i, &error);

      g_assert_no_error (error);
      g_assert_cmpfloat (ncm_vector_get (row, 1), ==, 100.0 * i + 1.0);
      ncm_vector_free (row);
    }

    ncm_mset_catalog_storage_get_rng_init (storage, &algo, &seed, NULL, &error);
    g_assert_cmpstr (algo, ==, ncm_rng_get_algo (rng));
    g_assert_cmpuint (seed, ==, 4321);

    stat = ncm_mset_catalog_storage_get_rng_stat (storage, &error);
    g_assert_cmpstr (stat, ==, rng_state);

    rng_file = ncm_mset_catalog_storage_get_rng (storage, &error);
    g_assert_no_error (error);
    g_assert_cmpfloat (ncm_rng_uniform_gen (rng_file, 0.0, 1.0), ==, ncm_rng_uniform_gen (rng, 0.0, 1.0));

    /* A second RNG is refused; a third recording appends. */
    ncm_mset_catalog_storage_set_rng (storage, rng, &error);
    g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RNG);
    g_clear_error (&error);

    ncm_mset_catalog_storage_start_recording (storage, "C", NULL, &error);
    {
      NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (12, 3);

      ncm_mset_catalog_storage_add_row (storage, row, &error);
      ncm_vector_free (row);
    }
    ncm_mset_catalog_storage_end_recording (storage, &error);
    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 13);
    g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 3);
    _test_ncm_mset_catalog_storage_fits_check_segment (storage, 2, 12, "C", "");

    ncm_rng_free (rng_file);
    g_free (algo);
    g_free (stat);
  }

  ncm_mset_catalog_storage_free (storage);

  /* The previous reader opens the new file. */
  {
    NcmMSetCatalog *mcat = _test_ncm_mset_catalog_storage_fits_open_catalog (test->filename);

    g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 13);
    g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, 2);
    g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);
    g_assert_cmpstr (ncm_mset_catalog_get_sampler (mcat), ==, "C");
    g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, 2);
    g_assert_cmpfloat (ncm_vector_get (ncm_mset_catalog_peek_row (mcat, 12), 0), ==, 1200.0);
    g_assert_nonnull (ncm_mset_catalog_peek_rng (mcat));

    ncm_mset_catalog_free (mcat);
  }

  g_free (rng_state);
  ncm_rng_free (rng);
}

/*
 * A Desc with the free parameters in the other order bound to a file: the rows come
 * out in the bound order, in memory and from the file.
 */
static void
test_ncm_mset_catalog_storage_fits_bind_order (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  GError *error                           = NULL;
  NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error);
  NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);
  NcmSerialize *ser                       = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmMSet *mset_rev                       = ncm_mset_dup (ncm_mset_catalog_desc_peek_mset (test->desc), ser);
  const gchar *fmap_rev[]                 = { "NcmModelMVND:mu_1", "NcmModelMVND:mu_0", NULL };
  const gchar *names[]                    = { "m2lnL", NULL };
  const gchar *symbols[]                  = { "-2\\ln(L)", NULL };
  NcmMSetCatalogDesc *desc_rev;
  guint i;

  ncm_mset_set_fmap (mset_rev, fmap_rev, FALSE, &error);
  g_assert_no_error (error);
  desc_rev = ncm_mset_catalog_desc_new (mset_rev, 2, FALSE, names, symbols, 0, NCM_MSET_CATALOG_RTYPE_ESMCMC, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_set_desc (storage, test->desc, &error);
  ncm_mset_catalog_storage_start_recording (storage, "A", NULL, &error);

  for (i = 0; i < 6; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, 3);

    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);

  /* Rows 0 and 1 are in the file only, rows 2 to 5 in memory. */
  ncm_mset_catalog_storage_set_desc (storage, desc_rev, &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset_catalog_storage_peek_desc (storage, &error) == desc_rev);

  for (i = 0; i < 6; i++)
  {
    NcmVector *row = ncm_mset_catalog_storage_get_row (storage, i, &error);

    g_assert_no_error (error);
    g_assert_cmpfloat (ncm_vector_get (row, 0), ==, 100.0 * i);
    g_assert_cmpfloat (ncm_vector_get (row, 1), ==, 100.0 * i + 2.0);
    g_assert_cmpfloat (ncm_vector_get (row, 2), ==, 100.0 * i + 1.0);
    ncm_vector_free (row);
  }

  /* A row added in the bound order is stored in the file order. */
  ncm_mset_catalog_storage_start_recording (storage, "B", NULL, &error);
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (6, 3);

    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }
  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);
  ncm_mset_catalog_storage_free (storage);

  {
    NcmMSetCatalogStorage *reopened = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, TRUE, &error));
    NcmVector *row                  = ncm_mset_catalog_storage_get_row (reopened, 6, &error);

    g_assert_no_error (error);
    g_assert_true (ncm_mset_catalog_desc_has_same_order (ncm_mset_catalog_storage_peek_desc (reopened, &error), test->desc));
    g_assert_no_error (error);
    g_assert_cmpfloat (ncm_vector_get (row, 1), ==, 602.0);
    g_assert_cmpfloat (ncm_vector_get (row, 2), ==, 601.0);

    ncm_vector_free (row);
    ncm_mset_catalog_storage_free (reopened);
  }

  ncm_mset_catalog_desc_free (desc_rev);
  ncm_mset_free (mset_rev);
  ncm_serialize_free (ser);
}

/*
 * A copy of esmcmc.fits opened writable: the first recording writes the translated run
 * type and the segment record; the file keeps the previous reader working.
 */
static void
test_ncm_mset_catalog_storage_fits_migrate (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  gchar *src    = _test_ncm_mset_catalog_storage_fits_fixture_path ("esmcmc.fits");
  GError *error = NULL;
  gchar *contents;
  gsize length;
  NcmMSetCatalogStorage *storage;
  NcmRNG *rng;
  guint i;

  g_file_get_contents (src, &contents, &length, &error);
  g_file_set_contents (test->filename, contents, length, &error);
  g_assert_no_error (error);
  g_free (contents);
  g_free (src);

  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error));
  g_assert_no_error (error);

  rng = ncm_mset_catalog_storage_get_rng (storage, &error);
  g_assert_nonnull (rng);

  ncm_mset_catalog_storage_start_recording (storage, "Walker-2", "w2", &error);
  g_assert_no_error (error);

  for (i = 0; i < TEST_FITS_ESMCMC_NWALKERS; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (TEST_FITS_ESMCMC_NROWS + i, 3);

    ncm_rng_uniform_gen (rng, 0.0, 1.0);
    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);
  ncm_rng_free (rng);
  ncm_mset_catalog_storage_free (storage);

  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, TRUE, &error));
  g_assert_no_error (error);
  g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (ncm_mset_catalog_storage_peek_desc (storage, &error)), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, TEST_FITS_ESMCMC_NROWS + TEST_FITS_ESMCMC_NWALKERS);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 3);
  _test_ncm_mset_catalog_storage_fits_check_segment (storage, 0, 0, "Multivariate Gaussian Sampler", "");
  _test_ncm_mset_catalog_storage_fits_check_segment (storage, 1, TEST_FITS_ESMCMC_NWALKERS, "Stretch-Move", "multi-stretch=no");
  _test_ncm_mset_catalog_storage_fits_check_segment (storage, 2, TEST_FITS_ESMCMC_NROWS, "Walker-2", "w2");
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, TEST_FITS_ESMCMC_NWALKERS);
  g_assert_no_error (error);
  ncm_mset_catalog_storage_free (storage);

  {
    NcmMSetCatalog *mcat = _test_ncm_mset_catalog_storage_fits_open_catalog (test->filename);

    g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FITS_ESMCMC_NROWS + TEST_FITS_ESMCMC_NWALKERS);
    g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);
    g_assert_cmpstr (ncm_mset_catalog_get_sampler (mcat), ==, "Walker-2");
    g_assert_cmpstr (ncm_mset_catalog_get_sampler_options (mcat), ==, "w2");
    g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, TEST_FITS_ESMCMC_NWALKERS);

    ncm_mset_catalog_free (mcat);
  }
}

/*
 * A new file with rows, two segments and a Markovian id, erased: empty in memory and
 * after reopening, the description and the RNG kept, a new recording from row 0, the
 * file read by the previous reader; erase refused inside a recording and read-only.
 */
static void
test_ncm_mset_catalog_storage_fits_erase (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  GError *error                  = NULL;
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, 99);
  NcmMSetCatalogStorage *storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error));
  guint i;

  g_assert_no_error (error);

  ncm_mset_catalog_storage_erase (storage, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_UNINITIALIZED);
  g_clear_error (&error);

  ncm_mset_catalog_storage_set_desc (storage, test->desc, &error);
  ncm_mset_catalog_storage_set_rng (storage, rng, &error);
  g_assert_no_error (error);
  g_object_set (storage, "smode", NCM_MSET_CATALOG_SYNC_AUTO, NULL);

  ncm_mset_catalog_storage_start_recording (storage, "A", "a-opts", &error);

  for (i = 0; i < 8; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, 3);

    ncm_rng_uniform_gen (rng, 0.0, 1.0);
    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_end_recording (storage, &error);
  ncm_mset_catalog_storage_start_recording (storage, "B", NULL, &error);
  ncm_mset_catalog_storage_set_markovian_id (storage, 4, &error);
  g_assert_no_error (error);

  ncm_mset_catalog_storage_erase (storage, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_RECORDING);
  g_clear_error (&error);

  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 2);

  ncm_mset_catalog_storage_erase (storage, &error);
  g_assert_no_error (error);

  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 0);
  g_assert_true (ncm_mset_catalog_storage_peek_desc (storage, &error) == test->desc);
  g_assert_null (ncm_mset_catalog_storage_peek_row (storage, 0, &error));
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_ROW);
  g_clear_error (&error);

  ncm_mset_catalog_storage_free (storage);

  /* Reopened writable: still empty, the RNG kept; a recording from row 0. */
  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error));
  g_assert_no_error (error);

  {
    gchar *algo = NULL;
    gulong seed = 0;

    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 0);
    g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 0);
    g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 0);
    g_assert_true (ncm_mset_catalog_desc_is_compatible (ncm_mset_catalog_storage_peek_desc (storage, &error), test->desc, NULL));

    ncm_mset_catalog_storage_get_rng_init (storage, &algo, &seed, NULL, &error);
    g_assert_no_error (error);
    g_assert_cmpstr (algo, ==, ncm_rng_get_algo (rng));
    g_assert_cmpuint (seed, ==, 99);
    g_free (algo);
  }

  ncm_mset_catalog_storage_start_recording (storage, "C", "c-opts", &error);

  for (i = 0; i < 3; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, 3);

    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);
  ncm_mset_catalog_storage_free (storage);

  /* Reopened read-only. */
  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, TRUE, &error));
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 3);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 1);
  _test_ncm_mset_catalog_storage_fits_check_segment (storage, 0, 0, "C", "c-opts");
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 0);

  for (i = 0; i < 3; i++)
  {
    NcmVector *row = ncm_mset_catalog_storage_peek_row (storage, i, &error);

    g_assert_no_error (error);
    g_assert_cmpfloat (ncm_vector_get (row, 0), ==, 100.0 * i);
  }

  ncm_mset_catalog_storage_erase (storage, &error);
  g_assert_error (error, NCM_MSET_CATALOG_STORAGE_ERROR, NCM_MSET_CATALOG_STORAGE_ERROR_READ_ONLY);
  g_clear_error (&error);
  ncm_mset_catalog_storage_free (storage);

  {
    NcmMSetCatalog *mcat = _test_ncm_mset_catalog_storage_fits_open_catalog (test->filename);

    g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 3);
    g_assert_cmpstr (ncm_mset_catalog_get_sampler (mcat), ==, "C");
    g_assert_cmpstr (ncm_mset_catalog_get_sampler_options (mcat), ==, "c-opts");
    g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, 0);
    g_assert_cmpfloat (ncm_vector_get (ncm_mset_catalog_peek_row (mcat, 2), 0), ==, 200.0);

    ncm_mset_catalog_free (mcat);
  }

  ncm_rng_free (rng);
}

/*
 * A copy of the ESMCMC fixture erased: no rows, no segments, the translated run type
 * written, readable by the previous reader.
 */
static void
test_ncm_mset_catalog_storage_fits_erase_legacy (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  gchar *src    = _test_ncm_mset_catalog_storage_fits_fixture_path ("esmcmc.fits");
  GError *error = NULL;
  gchar *contents;
  gsize length;
  NcmMSetCatalogStorage *storage;

  g_file_get_contents (src, &contents, &length, &error);
  g_file_set_contents (test->filename, contents, length, &error);
  g_assert_no_error (error);
  g_free (contents);
  g_free (src);

  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error));
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, TEST_FITS_ESMCMC_NROWS);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 2);

  ncm_mset_catalog_storage_erase (storage, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 0);
  ncm_mset_catalog_storage_free (storage);

  storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, TRUE, &error));
  g_assert_no_error (error);
  g_assert_cmpstr (ncm_mset_catalog_desc_get_run_type (ncm_mset_catalog_storage_peek_desc (storage, &error)), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (storage, &error), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_storage_get_markovian_id (storage, &error), ==, 0);
  {
    gchar *stat = ncm_mset_catalog_storage_get_rng_stat (storage, &error);

    g_assert_no_error (error);
    g_assert_nonnull (stat);
    g_free (stat);
  }
  ncm_mset_catalog_storage_free (storage);

  {
    NcmMSetCatalog *mcat = _test_ncm_mset_catalog_storage_fits_open_catalog (test->filename);

    g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 0);
    g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, NCM_MSET_CATALOG_RTYPE_ESMCMC);
    g_assert_null (ncm_mset_catalog_get_sampler (mcat));

    ncm_mset_catalog_free (mcat);
  }
}

/* A serialized storage names the same file and reopens it. */
static void
test_ncm_mset_catalog_storage_fits_serialize (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  GError *error                  = NULL;
  NcmSerialize *ser              = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmMSetCatalogStorage *storage = NCM_MSET_CATALOG_STORAGE (ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error));
  NcmMSetCatalogStorage *dup;
  guint i;

  ncm_mset_catalog_storage_set_desc (storage, test->desc, &error);
  ncm_mset_catalog_storage_start_recording (storage, "A", NULL, &error);

  for (i = 0; i < 4; i++)
  {
    NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, 3);

    ncm_mset_catalog_storage_add_row (storage, row, &error);
    ncm_vector_free (row);
  }

  ncm_mset_catalog_storage_end_recording (storage, &error);
  g_assert_no_error (error);

  dup = NCM_MSET_CATALOG_STORAGE (ncm_serialize_dup_obj (ser, G_OBJECT (storage)));
  ncm_mset_catalog_storage_free (storage);

  g_assert_cmpstr (ncm_mset_catalog_storage_fits_peek_filename (NCM_MSET_CATALOG_STORAGE_FITS (dup)), ==, test->filename);
  g_assert_cmpuint (ncm_mset_catalog_storage_nrows (dup, &error), ==, 4);
  g_assert_no_error (error);
  g_assert_cmpuint (ncm_mset_catalog_storage_n_samplers (dup, &error), ==, 1);

  ncm_mset_catalog_storage_free (dup);
  ncm_serialize_free (ser);
}

#ifdef HAVE_CFITSIO

/*
 * A table whose last column is a single float: the rows are read and written one column
 * at a time instead of as raw bytes. Four rows of integer values are written as doubles,
 * the column is replaced by a float one holding the same values, and the storage then
 * reads them back and appends two rows through the per-column path.
 */
static void
test_ncm_mset_catalog_storage_fits_float_column (TestNcmMSetCatalogStorageFits *test, gconstpointer pdata)
{
  const guint ncols = ncm_mset_catalog_desc_ncols (test->desc);
  GError *error     = NULL;
  guint i, c;

  {
    NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error);
    NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);

    g_assert_no_error (error);
    ncm_mset_catalog_storage_set_desc (storage, test->desc, &error);
    g_assert_no_error (error);
    ncm_mset_catalog_storage_start_recording (storage, "doubles", NULL, &error);
    g_assert_no_error (error);

    for (i = 0; i < 4; i++)
    {
      NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, ncols);

      ncm_mset_catalog_storage_add_row (storage, row, &error);
      g_assert_no_error (error);
      ncm_vector_free (row);
    }

    ncm_mset_catalog_storage_end_recording (storage, &error);
    g_assert_no_error (error);
    ncm_mset_catalog_storage_fits_free (storage_fits);
  }

  /* The last column becomes a float column with the same name and values. */
  {
    fitsfile *fptr = NULL;
    gint status    = 0;
    gchar name[FLEN_VALUE];
    gchar *ttype = g_strdup_printf ("TTYPE%u", ncols);
    gdouble values[4];

    fits_open_file (&fptr, test->filename, READWRITE, &status);
    fits_movnam_hdu (fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
    fits_read_key_str (fptr, ttype, name, NULL, &status);
    fits_read_col_dbl (fptr, ncols, 1, 1, 4, 0.0, values, NULL, &status);
    fits_delete_col (fptr, ncols, &status);
    fits_insert_col (fptr, ncols, name, "1E", &status);
    fits_write_col_dbl (fptr, ncols, 1, 1, 4, values, &status);
    fits_close_file (fptr, &status);
    g_assert_cmpint (status, ==, 0);
    g_free (ttype);
  }

  {
    NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (test->filename, FALSE, &error);
    NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);

    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 4);
    g_assert_no_error (error);

    for (i = 0; i < 4; i++)
    {
      NcmVector *row = ncm_mset_catalog_storage_peek_row (storage, i, &error);

      g_assert_no_error (error);

      for (c = 0; c < ncols; c++)
        g_assert_cmpfloat (ncm_vector_get (row, c), ==, 100.0 * i + c);
    }

    ncm_mset_catalog_storage_start_recording (storage, "floats", NULL, &error);
    g_assert_no_error (error);

    for (i = 4; i < 6; i++)
    {
      NcmVector *row = _test_ncm_mset_catalog_storage_fits_row (i, ncols);

      ncm_mset_catalog_storage_add_row (storage, row, &error);
      g_assert_no_error (error);
      ncm_vector_free (row);
    }

    ncm_mset_catalog_storage_end_recording (storage, &error);
    g_assert_no_error (error);
    ncm_mset_catalog_storage_fits_free (storage_fits);
  }

  {
    NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (test->filename, TRUE, &error);
    NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);

    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_catalog_storage_nrows (storage, &error), ==, 6);

    for (i = 0; i < 6; i++)
    {
      NcmVector *row = ncm_mset_catalog_storage_peek_row (storage, i, &error);

      g_assert_no_error (error);

      for (c = 0; c < ncols; c++)
        g_assert_cmpfloat (ncm_vector_get (row, c), ==, 100.0 * i + c);
    }

    ncm_mset_catalog_storage_fits_free (storage_fits);
  }
}

#endif /* HAVE_CFITSIO */


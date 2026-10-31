/***************************************************************************
 *            test_ncm_mset_catalog_fixtures.c
 *
 *  Sat Oct 31 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_mset_catalog_fixtures.c
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
 * Reads the catalog files of data/truth_tables/catalog/, written by
 * tests/tools/make_catalog_fixtures.py, with the read-only constructor and checks the
 * state each file carries: the run type and sampler strings of each producer, the
 * header variants (no RNG keys, no RNG_ALGO, no M2LNP_ID, FIRST_ID other than 0, no
 * MARKID, no rows), weighted rows and a `functions' key in HDU0.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include <numcosmo/numcosmo.h>

#define TEST_FIXTURES_ESMCMC_NROWS 18
#define TEST_FIXTURES_ESMCMC_NWALKERS 6
#define TEST_FIXTURES_FIRST_ID 6

static void test_ncm_mset_catalog_fixtures_esmcmc (void);
static void test_ncm_mset_catalog_fixtures_esmcmc_no_rng (void);
static void test_ncm_mset_catalog_fixtures_esmcmc_no_rng_algo (void);
static void test_ncm_mset_catalog_fixtures_esmcmc_no_m2lnp_id (void);
static void test_ncm_mset_catalog_fixtures_esmcmc_first_id (void);
static void test_ncm_mset_catalog_fixtures_esmcmc_no_markid (void);
static void test_ncm_mset_catalog_fixtures_esmcmc_header_only (void);
static void test_ncm_mset_catalog_fixtures_mc (void);
static void test_ncm_mset_catalog_fixtures_mcbs (void);
static void test_ncm_mset_catalog_fixtures_mcmc (void);
static void test_ncm_mset_catalog_fixtures_weighted (void);
static void test_ncm_mset_catalog_fixtures_functions (void);

static NcmMSetCatalog *_test_ncm_mset_catalog_fixtures_open (const gchar *name);
static void _test_ncm_mset_catalog_fixtures_check_esmcmc (NcmMSetCatalog *mcat);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

#ifdef HAVE_CFITSIO
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc", &test_ncm_mset_catalog_fixtures_esmcmc);
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc/no_rng", &test_ncm_mset_catalog_fixtures_esmcmc_no_rng);
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc/no_rng_algo", &test_ncm_mset_catalog_fixtures_esmcmc_no_rng_algo);
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc/no_m2lnp_id", &test_ncm_mset_catalog_fixtures_esmcmc_no_m2lnp_id);
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc/first_id", &test_ncm_mset_catalog_fixtures_esmcmc_first_id);
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc/no_markid", &test_ncm_mset_catalog_fixtures_esmcmc_no_markid);
  g_test_add_func ("/ncm/mset/catalog/fixtures/esmcmc/header_only", &test_ncm_mset_catalog_fixtures_esmcmc_header_only);
  g_test_add_func ("/ncm/mset/catalog/fixtures/mc", &test_ncm_mset_catalog_fixtures_mc);
  g_test_add_func ("/ncm/mset/catalog/fixtures/mcbs", &test_ncm_mset_catalog_fixtures_mcbs);
  g_test_add_func ("/ncm/mset/catalog/fixtures/mcmc", &test_ncm_mset_catalog_fixtures_mcmc);
  g_test_add_func ("/ncm/mset/catalog/fixtures/weighted", &test_ncm_mset_catalog_fixtures_weighted);
  g_test_add_func ("/ncm/mset/catalog/fixtures/functions", &test_ncm_mset_catalog_fixtures_functions);
#endif /* HAVE_CFITSIO */

  g_test_run ();
}

/* Opens data/truth_tables/catalog/name read-only, without burn-in. */
static NcmMSetCatalog *
_test_ncm_mset_catalog_fixtures_open (const gchar *name)
{
  gchar *rel_path      = g_strdup_printf ("truth_tables/catalog/%s", name);
  gchar *path          = ncm_cfg_get_data_filename (rel_path, TRUE);
  NcmMSetCatalog *mcat = ncm_mset_catalog_new_from_file_ro (path, 0);

  g_free (path);
  g_free (rel_path);

  return mcat;
}

/* The description and run strings shared by esmcmc.fits and its header variants. */
static void
_test_ncm_mset_catalog_fixtures_check_esmcmc (NcmMSetCatalog *mcat)
{
  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, "Ensemble Sampler MCMC");
  g_assert_cmpstr (ncm_mset_catalog_get_sampler (mcat), ==, "Stretch-Move");
  g_assert_cmpstr (ncm_mset_catalog_get_sampler_options (mcat), ==, "multi-stretch=no");
  g_assert_cmpstr (ncm_mset_catalog_get_initial_sampler (mcat), ==, "Multivariate Gaussian Sampler");
  g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, TEST_FIXTURES_ESMCMC_NWALKERS);
  g_assert_cmpuint (ncm_mset_catalog_nadd_vals (mcat), ==, 1);
  g_assert_cmpuint (ncm_mset_catalog_ncols (mcat), ==, 3);
  g_assert_cmpstr (ncm_mset_catalog_col_name (mcat, 0), ==, NCM_MSET_CATALOG_M2LNL_COLNAME);
  g_assert_cmpstr (ncm_mset_catalog_col_full_name (mcat, 1), ==, "NcmModelMVND:mu_0");
  g_assert_false (ncm_mset_catalog_weighted (mcat));
  g_assert_cmpuint (ncm_mset_fparams_len (ncm_mset_catalog_peek_mset (mcat)), ==, 2);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc.fits");
  NcmRNG *rng          = ncm_mset_catalog_peek_rng (mcat);

  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FIXTURES_ESMCMC_NROWS);
  g_assert_cmpint (ncm_mset_catalog_get_first_id (mcat), ==, 0);
  g_assert_cmpint (ncm_mset_catalog_get_m2lnp_var (mcat), ==, 0);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, TEST_FIXTURES_ESMCMC_NWALKERS);

  /*
   * The current library rebuilds the RNG from RNG_ALGO and RNG_INIS and does not
   * restore RNG_SEED.
   */
  g_assert_nonnull (rng);
  g_assert_cmpstr (ncm_rng_get_algo (rng), ==, "mt19937");

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc_no_rng (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc_no_rng.fits");

  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FIXTURES_ESMCMC_NROWS);
  g_assert_null (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc_no_rng_algo (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc_no_rng_algo.fits");

  /* The current library reads the RNG record only when RNG_ALGO is present. */
  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FIXTURES_ESMCMC_NROWS);
  g_assert_null (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc_no_m2lnp_id (void)
{
  NcmMSetCatalog *mcat;

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*does not contain `M2LNP_ID' key*");
  mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc_no_m2lnp_id.fits");
  g_test_assert_expected_messages ();

  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FIXTURES_ESMCMC_NROWS);
  g_assert_cmpint (ncm_mset_catalog_get_m2lnp_var (mcat), ==, -1);

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc_first_id (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc_first_id.fits");

  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FIXTURES_ESMCMC_NROWS);
  g_assert_cmpint (ncm_mset_catalog_get_first_id (mcat), ==, TEST_FIXTURES_FIRST_ID);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, TEST_FIXTURES_FIRST_ID + TEST_FIXTURES_ESMCMC_NWALKERS);

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc_no_markid (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc_no_markid.fits");

  /* Written before MARKID existed: every row is Markovian. */
  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, TEST_FIXTURES_ESMCMC_NROWS);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, 0);

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_esmcmc_header_only (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("esmcmc_header_only.fits");

  _test_ncm_mset_catalog_fixtures_check_esmcmc (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 0);
  g_assert_nonnull (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_mc (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("mc.fits");

  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, "from-model");
  g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, 1);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 10);
  g_assert_cmpuint (ncm_mset_catalog_nadd_vals (mcat), ==, 1);
  g_assert_cmpint (ncm_mset_catalog_get_m2lnp_var (mcat), ==, 0);
  g_assert_nonnull (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_mcbs (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("mcbs.fits");

  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, NCM_MSET_CATALOG_RTYPE_BSTRAP_MEAN);
  g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, 1);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 5);
  g_assert_cmpuint (ncm_mset_catalog_nadd_vals (mcat), ==, 1);
  g_assert_nonnull (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_mcmc (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("mcmc.fits");

  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, "Multivariate Gaussian Sampler");
  g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, 1);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 10);
  g_assert_cmpuint (ncm_mset_catalog_nadd_vals (mcat), ==, 1);
  g_assert_cmpint (ncm_mset_catalog_get_m2lnp_var (mcat), ==, 0);
  g_assert_nonnull (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_weighted (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_fixtures_open ("weighted.fits");

  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, "Fixture");
  g_assert_true (ncm_mset_catalog_weighted (mcat));
  g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, 1);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 10);
  g_assert_cmpuint (ncm_mset_catalog_nadd_vals (mcat), ==, 2);
  g_assert_cmpuint (ncm_mset_catalog_ncols (mcat), ==, 4);
  g_assert_cmpstr (ncm_mset_catalog_col_name (mcat, 0), ==, "m2lnL");
  g_assert_cmpstr (ncm_mset_catalog_col_name (mcat, 1), ==, NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME);
  g_assert_null (ncm_mset_catalog_peek_rng (mcat));

  ncm_mset_catalog_free (mcat);
}

static void
test_ncm_mset_catalog_fixtures_functions (void)
{
  NcmMSetCatalog *mcat   = _test_ncm_mset_catalog_fixtures_open ("functions.fits");
  NcmObjArray *functions = ncm_mset_catalog_peek_functions_array (mcat);

  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, "Fixture");
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 10);
  g_assert_cmpuint (ncm_mset_catalog_nadd_vals (mcat), ==, 2);
  g_assert_cmpstr (ncm_mset_catalog_col_name (mcat, 1), ==, "H0");

  g_assert_nonnull (functions);
  g_assert_cmpuint (ncm_obj_array_len (functions), ==, 1);
  g_assert_true (NCM_IS_MSET_FUNC (ncm_obj_array_peek (functions, 0)));

  ncm_mset_catalog_free (mcat);
}


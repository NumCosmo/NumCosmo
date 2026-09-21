/***************************************************************************
 *            test_ncm_mset_catalog.c
 *
 *  Wed February 07 10:36:55 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <vitenti@uel.br>
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
#ifdef HAVE_CFITSIO
#include <fitsio.h>
#endif /* HAVE_CFITSIO */

typedef struct _TestNcmMSetCatalog
{
  guint dim;
  NcmRNG *rng;
  NcmDataGaussCovMVND *data_mvnd;
  NcmMSetCatalog *mcat;
  guint ntests;
} TestNcmMSetCatalog;

void test_ncm_mset_catalog_new (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_new_2_chains (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_traps (TestNcmMSetCatalog *test, gconstpointer pdata);

void test_ncm_mset_catalog_free (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_mean (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_cov (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_norma (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_norma_bound (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_norma_unif (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_vol (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_bestfit (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_percentile (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_autocorrelation (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_tau_diagnostics (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_tau_frozen_walkers (void);
void test_ncm_mset_catalog_trim_oob_markovian (void);
void test_ncm_mset_catalog_weighted_tau_traps (void);
void test_ncm_mset_catalog_weighted_tau_subprocess (void);
void test_ncm_mset_catalog_accept_ratio_array (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_calc_param_ensemble_evol (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_calc_add_param_ensemble_evol (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_calc_add_param_ensemble_evol_short (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_invalid_run (TestNcmMSetCatalog *test, gconstpointer pdata);

#ifdef HAVE_CFITSIO
void test_ncm_mset_catalog_file_hdu0_roundtrip (void);
void test_ncm_mset_catalog_file_functions_array_roundtrip (void);
void test_ncm_mset_catalog_file_hdu0_legacy_object_format (void);
void test_ncm_mset_catalog_file_legacy_fallback (void);
void test_ncm_mset_catalog_file_missing_mset_traps (void);
void test_ncm_mset_catalog_file_missing_mset_subprocess (void);
void test_ncm_mset_catalog_file_peek_info (void);
void test_ncm_mset_catalog_file_multichain (void);
void test_ncm_mset_catalog_file_burnin_exceeds_traps (void);
void test_ncm_mset_catalog_file_burnin_exceeds_subprocess (void);
void test_ncm_mset_catalog_file_markovian_id_roundtrip (void);
void test_ncm_mset_catalog_file_markovian_id_missing_key (void);
void test_ncm_mset_catalog_file_markovian_id_backwards_traps (void);
void test_ncm_mset_catalog_file_markovian_id_backwards_subprocess (void);
void test_ncm_mset_catalog_file_markovian_id_trim (void);
void test_ncm_mset_catalog_file_markovian_id_frozen_needs_more (void);

#endif /* HAVE_CFITSIO */

typedef struct _TestNcmMSetCatalogTests
{
  const gchar *name;

  void (*func) (TestNcmMSetCatalog *test, gconstpointer pdata);
} TestNcmMSetCatalogTests;

TestNcmMSetCatalogTests fixtures[] =
{
  {"one_chain", test_ncm_mset_catalog_new},
  {"multiple_chain", test_ncm_mset_catalog_new_2_chains},
  {NULL, NULL}
};


TestNcmMSetCatalogTests tests[] =
{
  {"mean", test_ncm_mset_catalog_mean},
  {"cov", test_ncm_mset_catalog_cov},
  {"norma", test_ncm_mset_catalog_norma},
  {"norma/bound", test_ncm_mset_catalog_norma_bound},
  {"norma/unif", test_ncm_mset_catalog_norma_unif},
  {"vol", test_ncm_mset_catalog_vol},
  {"bestfit", test_ncm_mset_catalog_bestfit},
  {"percentile", test_ncm_mset_catalog_percentile},
  {"autocorrelation", test_ncm_mset_catalog_autocorrelation},
  {"tau_diagnostics", test_ncm_mset_catalog_tau_diagnostics},
  {"accept_ratio_array", test_ncm_mset_catalog_accept_ratio_array},
  {"calc_param_ensemble_evol", test_ncm_mset_catalog_calc_param_ensemble_evol},
  {"calc_add_param_ensemble_evol", test_ncm_mset_catalog_calc_add_param_ensemble_evol},
  {"calc_add_param_ensemble_evol/short", test_ncm_mset_catalog_calc_add_param_ensemble_evol_short},
  {NULL, NULL}
};

gint
main (gint argc, gchar *argv[])
{
  gint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  for (i = 0; fixtures[i].name != NULL; i++)
  {
    for (j = 0; tests[j].name != NULL; j++)
    {
      gchar *path = g_strdup_printf ("/ncm/mset/catalog/%s/%s",
                                     fixtures[i].name, tests[j].name);

      g_test_add (path, TestNcmMSetCatalog, NULL,
                  fixtures[i].func,
                  tests[j].func,
                  &test_ncm_mset_catalog_free);

      g_free (path);
    }
  }

  g_test_add ("/ncm/mset/catalog/traps", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_traps,
              &test_ncm_mset_catalog_free);

  g_test_add ("/ncm/mset/catalog/invalid/run/subprocess", TestNcmMSetCatalog, NULL,
              &test_ncm_mset_catalog_new,
              &test_ncm_mset_catalog_invalid_run,
              &test_ncm_mset_catalog_free);

#ifdef HAVE_CFITSIO
  g_test_add_func ("/ncm/mset/catalog/file/hdu0_roundtrip", &test_ncm_mset_catalog_file_hdu0_roundtrip);
  g_test_add_func ("/ncm/mset/catalog/file/functions_array_roundtrip", &test_ncm_mset_catalog_file_functions_array_roundtrip);
  g_test_add_func ("/ncm/mset/catalog/file/hdu0_legacy_object_format", &test_ncm_mset_catalog_file_hdu0_legacy_object_format);
  g_test_add_func ("/ncm/mset/catalog/file/legacy_fallback", &test_ncm_mset_catalog_file_legacy_fallback);
  g_test_add_func ("/ncm/mset/catalog/file/missing_mset/traps", &test_ncm_mset_catalog_file_missing_mset_traps);
  g_test_add_func ("/ncm/mset/catalog/file/missing_mset/subprocess", &test_ncm_mset_catalog_file_missing_mset_subprocess);
  g_test_add_func ("/ncm/mset/catalog/weighted/tau/traps", &test_ncm_mset_catalog_weighted_tau_traps);
  g_test_add_func ("/ncm/mset/catalog/weighted/tau/subprocess", &test_ncm_mset_catalog_weighted_tau_subprocess);
  g_test_add_func ("/ncm/mset/catalog/tau/frozen_keff", &test_ncm_mset_catalog_tau_frozen_walkers);
  g_test_add_func ("/ncm/mset/catalog/trim_oob/markovian", &test_ncm_mset_catalog_trim_oob_markovian);

  g_test_add_func ("/ncm/mset/catalog/file/peek_info", &test_ncm_mset_catalog_file_peek_info);
  g_test_add_func ("/ncm/mset/catalog/file/multichain", &test_ncm_mset_catalog_file_multichain);
  g_test_add_func ("/ncm/mset/catalog/file/burnin_exceeds/traps", &test_ncm_mset_catalog_file_burnin_exceeds_traps);
  g_test_add_func ("/ncm/mset/catalog/file/burnin_exceeds/subprocess", &test_ncm_mset_catalog_file_burnin_exceeds_subprocess);
  g_test_add_func ("/ncm/mset/catalog/file/markovian_id/roundtrip", &test_ncm_mset_catalog_file_markovian_id_roundtrip);
  g_test_add_func ("/ncm/mset/catalog/file/markovian_id/missing_key", &test_ncm_mset_catalog_file_markovian_id_missing_key);
  g_test_add_func ("/ncm/mset/catalog/file/markovian_id/backwards/traps", &test_ncm_mset_catalog_file_markovian_id_backwards_traps);
  g_test_add_func ("/ncm/mset/catalog/file/markovian_id/backwards/subprocess", &test_ncm_mset_catalog_file_markovian_id_backwards_subprocess);
  g_test_add_func ("/ncm/mset/catalog/file/markovian_id/trim", &test_ncm_mset_catalog_file_markovian_id_trim);
  g_test_add_func ("/ncm/mset/catalog/file/markovian_id/frozen_needs_more", &test_ncm_mset_catalog_file_markovian_id_frozen_needs_more);
#endif /* HAVE_CFITSIO */

  g_test_run ();
}

#define NTESTS_MIN 5000
#define NTESTS_MAX 10000

void
test_ncm_mset_catalog_new (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  const guint dim                = test->dim = g_test_rand_int_range (2, 5);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (dim, 5.0e-3, 1.0e-2, 1.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmVector *y                   = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmMSetCatalog *mcat;

  /*ncm_model_param_set_lower_bound (NCM_MODEL (model_mvnd), 0, 0.0);*/

  ncm_mset_param_set_vector (mset, y);
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE,
                               "m2lnL", "-2\\ln(L)",
                               NULL);

  ncm_mset_catalog_set_m2lnp_var (mcat, 0);

  test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd);
  test->mcat      = ncm_mset_catalog_ref (mcat);
  test->rng       = rng;

  g_assert_true (NCM_IS_MSET_CATALOG (test->mcat));

  ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_mset_clear (&mset);
  ncm_mset_catalog_clear (&mcat);
}

void
test_ncm_mset_catalog_new_2_chains (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  const gint dim                 = test->dim = g_test_rand_int_range (2, 5);
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (dim, 5.0e-3, 1.0e-2, 1.0, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmVector *y                   = ncm_data_gauss_cov_peek_mean (NCM_DATA_GAUSS_COV (data_mvnd));
  NcmMSetCatalog *mcat;

  /*ncm_model_param_set_lower_bound (NCM_MODEL (model_mvnd), 0, 0.0);*/

  ncm_mset_param_set_vector (mset, y);
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, 20, FALSE,
                               "m2lnL", "-2\\ln(L)",
                               NULL);

  ncm_mset_catalog_set_m2lnp_var (mcat, 0);

  test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd);
  test->mcat      = ncm_mset_catalog_ref (mcat);
  test->rng       = rng;

  g_assert_true (NCM_IS_MSET_CATALOG (test->mcat));

  ncm_data_gauss_cov_mvnd_clear (&data_mvnd);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_mset_clear (&mset);
  ncm_mset_catalog_clear (&mcat);
}

void
test_ncm_mset_catalog_free (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_mset_catalog_free, test->mcat);
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (test->data_mvnd));
  NCM_TEST_FREE (ncm_rng_free, test->rng);
}

void
test_ncm_mset_catalog_mean (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_data_resample (data, mset, test->rng);
    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  if (FALSE)
  {
    ncm_cfg_msg_sepa ();
    ncm_mset_catalog_estimate_autocorrelation_tau (test->mcat, FALSE);
    ncm_mset_catalog_log_current_stats (test->mcat);
    ncm_mset_params_log_vals (mset);
  }

  {
    NcmVector *mean = NULL;
    NcmVector *p    = ncm_vector_new (test->dim);

    ncm_mset_catalog_get_mean (test->mcat, &mean);
    ncm_mset_param_get_vector (mset, p);

    ncm_vector_cmp (p, mean);

    g_assert_cmpfloat (ncm_vector_get_max (p), <, 1.0e-2);

    ncm_vector_clear (&mean);
    ncm_vector_clear (&p);
  }
}

void
test_ncm_mset_catalog_cov (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  NcmMatrix *data_cov  = ncm_data_gauss_cov_peek_cov (cov);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_data_resample (data, mset, test->rng);
    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  {
    NcmMatrix *cat_cov = NULL;

    ncm_mset_catalog_get_covar (test->mcat, &cat_cov);

    g_assert_cmpfloat (ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0), <, 5.0e-1);

    ncm_matrix_cov2cor (data_cov, data_cov);
    ncm_matrix_cov2cor (cat_cov, cat_cov);

    g_assert_cmpfloat (ncm_matrix_cmp (cat_cov, data_cov, 1.0), <, 1.0e-1);

    if (FALSE)
    {
      ncm_cfg_msg_sepa ();
      ncm_matrix_log_vals (cat_cov,  "# CAT  COV: ", "% 12.5g");
      ncm_matrix_log_vals (data_cov, "# DATA COV: ", "% 12.5g");

      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cat_cov, data_cov, 0.0));
      printf ("# WDIFFD  : % 22.15e\n", ncm_matrix_cmp_diag (cat_cov, data_cov, 0.0));

      ncm_matrix_sub (cat_cov, data_cov);
      ncm_matrix_div_elements (cat_cov, data_cov);

      ncm_matrix_log_vals (cat_cov,  "# CMP     : ", "% 12.5e");
    }
  }
}

void
test_ncm_mset_catalog_norma (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX) * 10;
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  gdouble ratio;
  gdouble lnnorm_sd;
  gulong N;
  gulong Nin;
  guint i;

  N     = Nin = 0;
  ratio = 0.0;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;
    gulong Ni;

    ncm_data_gauss_cov_mvnd_gen (test->data_mvnd, mset, mset, (NcmDataGaussCovMVNDBound) ncm_mset_fparam_valid_bounds, test->rng, &Ni);
    N += Ni;

    ncm_data_m2lnL_val (data, mset, &m2lnL);
    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  while (TRUE)
  {
    gdouble err_rel;
    gulong Ni;

    ncm_data_gauss_cov_mvnd_gen (test->data_mvnd, mset, mset, (NcmDataGaussCovMVNDBound) ncm_mset_fparam_valid_bounds, test->rng, &Ni);
    N += Ni;
    i++;

    ratio   = i * 1.0 / (1.0 * N);
    err_rel = sqrt ((1.0 - ratio) / (N * ratio));

    if (err_rel < 1.0e-3)
      break;
  }

  ncm_mset_catalog_get_post_lnnorm (test->mcat, &lnnorm_sd);

  ncm_assert_cmpdouble_e (ncm_mset_catalog_get_post_lnnorm (test->mcat, &lnnorm_sd), ==, log (ratio), 0.2, 1.0e-3);
}

void
test_ncm_mset_catalog_norma_bound (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX) * 100;
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  gdouble ratio;
  gdouble lnnorm_sd;
  gulong N;
  gulong Nin;
  gulong i;

  for (i = 0; i < test->dim; i++)
  {
    ncm_model_param_set_upper_bound (ncm_mset_peek (mset, ncm_model_mvnd_id ()), i, ncm_model_param_get (ncm_mset_peek (mset, ncm_model_mvnd_id ()), i) * 1.01);
  }

  N     = Nin = 0;
  ratio = 0.0;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;
    gulong Ni;

    ncm_data_gauss_cov_mvnd_gen (test->data_mvnd, mset, mset, (NcmDataGaussCovMVNDBound) ncm_mset_fparam_valid_bounds, test->rng, &Ni);
    N += Ni;

    ncm_data_m2lnL_val (data, mset, &m2lnL);
    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  ratio = ncm_data_gauss_cov_mvnd_est_ratio (test->data_mvnd, mset, mset, (NcmDataGaussCovMVNDBound) ncm_mset_fparam_valid_bounds, &N, &i, 1.0e-3, test->rng);

  /*printf ("<% 22.15g % 22.15g % 12.5e %d>\n", ncm_mset_catalog_get_post_lnnorm (test->mcat), log (ratio), fabs (expm1 (ncm_mset_catalog_get_post_lnnorm (test->mcat) - log (ratio))), test->dim);*/

  ncm_assert_cmpdouble_e (ncm_mset_catalog_get_post_lnnorm (test->mcat, &lnnorm_sd), ==, log (ratio), 0.2, 1.0e-3);
}

void
test_ncm_mset_catalog_norma_unif (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX) * 10;
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  gdouble norma        = 1.0;
  gdouble lnnorm_sd;
  gulong i;

  for (i = 0; i < test->dim; i++)
  {
    ncm_model_param_set_upper_bound (ncm_mset_peek (mset, ncm_model_mvnd_id ()), i, ncm_model_param_get (ncm_mset_peek (mset, ncm_model_mvnd_id ()), i) * 1.05);
    {
      const gdouble ub = ncm_model_param_get_upper_bound (ncm_mset_peek (mset, ncm_model_mvnd_id ()), i);
      const gdouble lb = ncm_model_param_get_lower_bound (ncm_mset_peek (mset, ncm_model_mvnd_id ()), i);

      norma *= (ub - lb);
    }
  }

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;
    guint j;

    for (j = 0; j < test->dim; j++)
    {
      const gdouble ub = ncm_model_param_get_upper_bound (ncm_mset_peek (mset, ncm_model_mvnd_id ()), j);
      const gdouble lb = ncm_model_param_get_lower_bound (ncm_mset_peek (mset, ncm_model_mvnd_id ()), j);

      ncm_vector_set (y, j, ncm_rng_uniform_gen (test->rng, lb, ub));
    }

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  /*printf ("<% 22.15g % 22.15g % 12.5e %d>\n", ncm_mset_catalog_get_post_lnnorm (test->mcat), log (norma), fabs (expm1 (ncm_mset_catalog_get_post_lnnorm (test->mcat) - log (norma))), test->dim);*/

  ncm_assert_cmpdouble_e (ncm_mset_catalog_get_post_lnnorm (test->mcat, &lnnorm_sd), ==, log (norma), 0.2, 0.0);
}

void
test_ncm_mset_catalog_vol (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  gdouble glnvol;
  gdouble lnnorm_sd;
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  {
    const gdouble lnevol = ncm_mset_catalog_get_post_lnvol (test->mcat, 0.6827, &glnvol);

    ncm_assert_cmpdouble_e (lnevol, ==, glnvol, 0.2, 0.0);
  }

  if (FALSE)
  {
    printf ("# DIM %d LNNORMA = % 22.15g\n", test->dim, ncm_mset_catalog_get_post_lnnorm (test->mcat, &lnnorm_sd));
    printf ("# DIM %d VOL1SIG = % 22.15g ", test->dim, ncm_mset_catalog_get_post_lnvol (test->mcat, 0.6827, &glnvol));
    printf ("% 22.15g\n", glnvol);
  }
}

void
test_ncm_mset_catalog_bestfit (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  gdouble m2lnL_min    = GSL_POSINF;
  NcmVector *min_row   = ncm_vector_new (test->dim);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    if (m2lnL < m2lnL_min)
    {
      ncm_vector_memcpy (min_row, y);
      m2lnL_min = m2lnL;
    }

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  {
    NcmVector *bestfit          = ncm_mset_catalog_get_bestfit_row (test->mcat);
    NcmVector *bestfit_row      = ncm_vector_get_subvector (bestfit, 1, test->dim);
    const gdouble m2lnL_bestfit = ncm_mset_catalog_get_bestfit_m2lnL (test->mcat);

    g_assert_cmpfloat_with_epsilon (m2lnL_bestfit, m2lnL_min, 1.0e-15);
    g_assert_cmpfloat (ncm_vector_get (bestfit, 0), ==, m2lnL_bestfit);

    ncm_vector_cmp (min_row, bestfit_row);

    g_assert_cmpfloat (ncm_vector_sum_cpts (min_row), ==, 0.0);

    ncm_vector_free (bestfit);
    ncm_vector_free (bestfit_row);
  }

  ncm_vector_free (min_row);
}

static gint
_cmp_double (gconstpointer a, gconstpointer b)
{
  const gdouble *da = (const gdouble *) a;
  const gdouble *db = (const gdouble *) b;

  if (*da < *db)
    return -1;
  else if (*da > *db)
    return 1;
  else
    return 0;
}

void
test_ncm_mset_catalog_percentile (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  gdouble m2lnL_min    = GSL_POSINF;
  NcmVector *min_row   = ncm_vector_new (test->dim);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  GArray *m2lnL_array  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint nth            = 0;
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    if (m2lnL < m2lnL_min)
    {
      ncm_vector_memcpy (min_row, y);
      m2lnL_min = m2lnL;
    }

    g_array_append_val (m2lnL_array, m2lnL);
    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  g_array_sort (m2lnL_array, _cmp_double);

  {
    NcmVector *bestfit          = ncm_mset_catalog_get_bestfit_row (test->mcat);
    NcmVector *bestfit_row      = ncm_vector_get_subvector (bestfit, 1, test->dim);
    const gdouble m2lnL_bestfit = ncm_mset_catalog_get_bestfit_m2lnL (test->mcat);

    g_assert_cmpfloat_with_epsilon (m2lnL_bestfit, m2lnL_min, 1.0e-15);
    g_assert_cmpfloat (ncm_vector_get (bestfit, 0), ==, m2lnL_bestfit);
    g_assert_cmpfloat (ncm_vector_get (bestfit, 0), ==, g_array_index (m2lnL_array, gdouble, 0));

    ncm_vector_cmp (min_row, bestfit_row);

    g_assert_cmpfloat (ncm_vector_sum_cpts (min_row), ==, 0.0);

    ncm_vector_free (bestfit);
    ncm_vector_free (bestfit_row);
  }

  for (i = 0; i < 99; i++)
  {
    const gdouble p        = (i + 1.0) / 100.0;
    const gdouble m2lnL_p1 = ncm_mset_catalog_get_nth_m2lnL_percentile (test->mcat, p, &nth);
    const gdouble m2lnL_p2 = g_array_index (m2lnL_array, gdouble, (guint) (p * (gdouble) nt));

    g_assert_cmpuint (nth, ==, (guint) (p * (gdouble) nt));

    g_assert_cmpfloat_with_epsilon (m2lnL_p1, m2lnL_p2, 1.0e-15);
  }

  ncm_vector_free (min_row);
  g_array_unref (m2lnL_array);
}

void
test_ncm_mset_catalog_autocorrelation (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (test->mcat, FALSE);
  {
    NcmVector *tau = ncm_mset_catalog_peek_autocorrelation_tau (test->mcat);

    g_assert_true (ncm_vector_is_finite (tau));
  }

  {
    GArray *accept_ratio_array = ncm_mset_catalog_peek_accept_ratio_array (test->mcat);

    if (ncm_mset_catalog_nchains (test->mcat) > 1)
      g_assert_true (accept_ratio_array != NULL);
    else
      g_assert_true (accept_ratio_array == NULL);
  }
}

/* The autocorrelation accumulator the catalog keeps as rows are added, and the quantities
 * built from it: the effective number of chains per iteration, the effective sample size
 * that follows from it, and the conditions attached to each estimate. */

/* The autocorrelation accumulator forms unweighted lagged sums, so a weighted catalog has
 * to refuse rather than report a tau built from them. */
static NcmMSetCatalog *
_test_ncm_mset_catalog_weighted_new (void)
{
  NcmRNG *rng              = ncm_rng_seeded_new (NULL, 314159);
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (2);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  const gchar *names[]     = {"m2lnL", "w", NULL};
  const gchar *symbols[]   = {"-2\\ln(L)", "w", NULL};
  NcmMSetCatalog *mcat;
  NcmVector *row;
  guint i, p;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new_array (mset, 2, 1, TRUE, (gchar **) names, (gchar **) symbols);
  row  = ncm_vector_new (ncm_mset_catalog_ncols (mcat));

  for (i = 0; i < 100; i++)
  {
    ncm_vector_set (row, 0, ncm_rng_ugaussian_gen (rng));
    ncm_vector_set (row, 1, 1.0 + 0.1 * fabs (ncm_rng_ugaussian_gen (rng)));

    for (p = 2; p < ncm_mset_catalog_ncols (mcat); p++)
      ncm_vector_set (row, p, ncm_rng_ugaussian_gen (rng));

    ncm_mset_catalog_add_from_vector (mcat, row);
  }

  ncm_vector_free (row);
  ncm_mset_clear (&mset);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_rng_free (rng);

  return mcat;
}

void
test_ncm_mset_catalog_weighted_tau_traps (void)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/weighted/tau/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*does not support weighted catalogs*");
}

void
test_ncm_mset_catalog_weighted_tau_subprocess (void)
{
  NcmMSetCatalog *mcat = _test_ncm_mset_catalog_weighted_new ();

  /* The catalog itself works; only the autocorrelation quantities refuse. */
  g_assert_true (ncm_mset_catalog_weighted (mcat));
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 100);

  ncm_mset_catalog_estimate_autocorrelation_tau (mcat, FALSE);

  ncm_mset_catalog_free (mcat);
}

void
test_ncm_mset_catalog_tau_diagnostics (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nchains  = ncm_mset_catalog_nchains (test->mcat);
  const guint nt       = nchains * g_test_rand_int_range (NTESTS_MIN / nchains, NTESTS_MAX / nchains);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  NcmStatsAcorr *acorr = NULL;
  const guint total    = ncm_mset_catalog_ncols (test->mcat);
  gboolean any_flagged = FALSE;
  guint i, p;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (test->mcat, FALSE);

  /* The accumulator holds one entry per complete iteration, not one per row. */
  acorr = ncm_mset_catalog_peek_acorr (test->mcat);
  g_assert_true (NCM_IS_STATS_ACORR (acorr));
  g_assert_cmpuint (ncm_stats_acorr_len (acorr), ==, total);
  g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, 0), ==, nt / nchains);

  for (p = 0; p < total; p++)
  {
    const gdouble keff = ncm_mset_catalog_get_keff (test->mcat, p);
    const gdouble ess  = ncm_mset_catalog_get_ess (test->mcat, p);

    /* The catalog reports what the accumulator holds. */
    g_assert_cmpuint (ncm_mset_catalog_get_tau_diag (test->mcat, p), ==,
                      ncm_stats_acorr_get_diag (acorr, p));
    ncm_assert_cmpdouble_e (ncm_vector_get (ncm_mset_catalog_peek_autocorrelation_tau (test->mcat), p),
                            ==, ncm_stats_acorr_get_tau (acorr, p), 1.0e-14, 0.0);

    /* A single chain is one chain per iteration by definition; more than one is bounded
     * below by nothing but positivity, since walkers may be correlated or held apart. */
    if (nchains == 1)
      ncm_assert_cmpdouble_e (keff, ==, 1.0, 1.0e-14, 0.0);
    else
      g_assert_cmpfloat (keff, >, 0.0);

    g_assert_true (gsl_finite (keff));
    ncm_assert_cmpdouble_e (ess, ==, keff * ncm_stats_acorr_get_ess (acorr, p), 1.0e-12, 0.0);
    g_assert_cmpfloat (ess, >, 0.0);

    if (ncm_mset_catalog_get_tau_diag (test->mcat, p) != NCM_STATS_ACORR_DIAG_OK)
      any_flagged = TRUE;
  }

  g_assert_cmpint (ncm_mset_catalog_log_tau_diag (test->mcat), ==, any_flagged);

  /* The sampling gate asks for more only while a free parameter is shorter than the
   * reliability factor times its own tau, and the length it asks for is that product. */
  {
    guint req_niter    = 0;
    gboolean needs     = ncm_mset_catalog_tau_needs_more (test->mcat, &req_niter);
    gboolean any_short = FALSE;
    gdouble max_tau    = 1.0;
    const guint fpi    = ncm_mset_catalog_nadd_vals (test->mcat);
    const guint fpf    = fpi + ncm_mset_fparams_len (ncm_mset_catalog_peek_mset (test->mcat));

    for (p = fpi; p < fpf; p++)
    {
      /* The gate opens on any of the three conditions the catalog refuses to call a
       * sample size, not on the chain being short alone. */
      if (ncm_mset_catalog_get_tau_diag (test->mcat, p) &
          (NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_ZERO_VARIANCE))
        any_short = TRUE;

      if ((nchains > 1) && (ncm_mset_catalog_get_keff (test->mcat, p) > 2.0 * nchains))
        any_short = TRUE;

      max_tau = GSL_MAX (max_tau, ncm_stats_acorr_get_tau (acorr, p));
    }

    g_assert_cmpint (needs, ==, any_short);
    g_assert_cmpuint (req_niter, ==, (guint) ceil (ncm_stats_acorr_get_reliability_factor (acorr) * max_tau));

    /* A factor of one is met by any chain of more than one iteration, so the gate closes. */
    ncm_stats_acorr_set_reliability_factor (acorr, 1.0);
    g_assert_false (ncm_mset_catalog_tau_needs_more (test->mcat, NULL));

    /* Asking for far more than the chain holds opens it again. */
    ncm_stats_acorr_set_reliability_factor (acorr, 1.0e6);
    g_assert_true (ncm_mset_catalog_tau_needs_more (test->mcat, &req_niter));
    g_assert_cmpuint (req_niter, >, ncm_stats_acorr_nitens (acorr, 0));

    ncm_stats_acorr_set_reliability_factor (acorr, NCM_STATS_ACORR_DEFAULT_RELIABILITY_FACTOR);
  }

  /* The error the sampler runs against is built from the effective sample size. */
  {
    const gdouble lerror = ncm_mset_catalog_largest_error (test->mcat);

    g_assert_true (gsl_finite (lerror));
    g_assert_cmpfloat (lerror, >, 0.0);
  }

  /* Every estimator is accepted, reported back, and leaves the accumulated data alone. */
  {
    const NcmStatsAcorrMethod methods[] = {
      NCM_STATS_ACORR_METHOD_AR,
      NCM_STATS_ACORR_METHOD_GEYER,
      NCM_STATS_ACORR_METHOD_SOKAL,
      NCM_STATS_ACORR_METHOD_MAX
    };

    for (i = 0; i < G_N_ELEMENTS (methods); i++)
    {
      ncm_mset_catalog_set_tau_method (test->mcat, methods[i]);
      g_assert_cmpuint (ncm_mset_catalog_get_tau_method (test->mcat), ==, methods[i]);

      ncm_mset_catalog_estimate_autocorrelation_tau (test->mcat, FALSE);
      g_assert_true (ncm_vector_is_finite (ncm_mset_catalog_peek_autocorrelation_tau (test->mcat)));
      g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, 0), ==, nt / nchains);
    }
  }

  /* Treating the interleaved rows as one series is a different measurement of a different
   * thing, and it has to produce a number. */
  ncm_mset_catalog_estimate_autocorrelation_tau (test->mcat, TRUE);
  g_assert_true (ncm_vector_is_finite (ncm_mset_catalog_peek_autocorrelation_tau (test->mcat)));
}

void
test_ncm_mset_catalog_accept_ratio_array (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  guint i;

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  {
    GArray *accept_ratio_array = ncm_mset_catalog_peek_accept_ratio_array (test->mcat);

    if (ncm_mset_catalog_nchains (test->mcat) > 1)
    {
      const gint max_time = ncm_mset_catalog_max_time (test->mcat);

      g_assert_true (accept_ratio_array != NULL);

      g_assert_cmpuint (accept_ratio_array->len + 1, ==, max_time);
    }
    else
    {
      g_assert_true (accept_ratio_array == NULL);
    }
  }
}

void
test_ncm_mset_catalog_calc_param_ensemble_evol (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  guint i;

  if (ncm_mset_catalog_nchains (test->mcat) == 1)
  {
    g_test_skip ("Single chain");

    return;
  }

  for (i = 0; i < nt; i++)
  {
    gdouble m2lnL = 0.0;

    ncm_data_resample (data, mset, test->rng);
    ncm_data_m2lnL_val (data, mset, &m2lnL);

    ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
  }

  {
    NcmVector *pval             = NULL;
    NcmMatrix *t_evol           = NULL;
    const NcmMSetPIndex *pindex = ncm_mset_fparam_get_pi (mset, 0);
    guint i;

    ncm_mset_catalog_calc_param_ensemble_evol (test->mcat, pindex, 100, NCM_FIT_RUN_MSGS_NONE, &pval, &t_evol);

    g_assert_nonnull (pval);
    g_assert_nonnull (t_evol);
    g_assert_cmpuint (ncm_vector_len (pval), ==, ncm_matrix_ncols (t_evol));

    for (i = 0; i < ncm_vector_len (pval); i++)
    {
      const gdouble pval_i = ncm_vector_get (pval, i);

      g_assert_cmpfloat (pval_i, >=, 0.0);
      g_assert_cmpfloat (pval_i, <=, 3.0);
    }

    for (i = 0; i < ncm_matrix_nrows (t_evol); i++)
    {
      guint j;

      for (j = 0; j < ncm_matrix_ncols (t_evol); j++)
      {
        const gdouble tval_ij = ncm_matrix_get (t_evol, i, j);

        g_assert_cmpfloat (tval_ij, >=, 0.0);
      }
    }

    ncm_vector_clear (&pval);
    ncm_matrix_clear (&t_evol);
  }
}

void
test_ncm_mset_catalog_calc_add_param_ensemble_evol (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  const guint nt       = g_test_rand_int_range (NTESTS_MIN, NTESTS_MAX);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  guint i;

  if (ncm_mset_catalog_nchains (test->mcat) == 1)
  {
    g_test_skip ("Single chain");

    return;
  }

  if (g_test_subprocess ())
  {
    for (i = 0; i < nt; i++)
    {
      gdouble m2lnL = 0.0;

      ncm_data_resample (data, mset, test->rng);
      ncm_data_m2lnL_val (data, mset, &m2lnL);

      ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
    }

    {
      NcmVector *pval   = NULL;
      NcmMatrix *t_evol = NULL;
      guint i;

      ncm_mset_catalog_calc_add_param_ensemble_evol (test->mcat, 0, 100, NCM_FIT_RUN_MSGS_FULL, &pval, &t_evol);

      g_assert_nonnull (pval);
      g_assert_nonnull (t_evol);
      g_assert_cmpuint (ncm_vector_len (pval), ==, ncm_matrix_ncols (t_evol));

      for (i = 0; i < ncm_vector_len (pval); i++)
      {
        const gdouble pval_i = ncm_vector_get (pval, i);

        g_assert_true (gsl_finite (pval_i));
      }

      for (i = 0; i < ncm_matrix_nrows (t_evol); i++)
      {
        guint j;

        for (j = 0; j < ncm_matrix_ncols (t_evol); j++)
        {
          const gdouble tval_ij = ncm_matrix_get (t_evol, i, j);

          g_assert_cmpfloat (tval_ij, >=, 0.0);
        }
      }

      ncm_vector_clear (&pval);
      ncm_matrix_clear (&t_evol);
    }

    return;
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*Calculating evolution to time*");
}

void
test_ncm_mset_catalog_calc_add_param_ensemble_evol_short (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  NcmData *data        = NCM_DATA (test->data_mvnd);
  NcmDataGaussCov *cov = NCM_DATA_GAUSS_COV (test->data_mvnd);
  NcmMSet *mset        = ncm_mset_catalog_peek_mset (test->mcat);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (cov);
  guint i;

  if (ncm_mset_catalog_nchains (test->mcat) == 1)
  {
    g_test_skip ("Single chain");

    return;
  }

  if (g_test_subprocess ())
  {
    for (i = 0; i < ncm_mset_catalog_nchains (test->mcat) * 5; i++)
    {
      gdouble m2lnL = 0.0;

      ncm_data_resample (data, mset, test->rng);
      ncm_data_m2lnL_val (data, mset, &m2lnL);

      ncm_mset_catalog_add_from_vector_array (test->mcat, y, &m2lnL);
    }

    {
      NcmVector *pval   = NULL;
      NcmMatrix *t_evol = NULL;
      guint i;

      ncm_mset_catalog_calc_add_param_ensemble_evol (test->mcat, 0, 100, NCM_FIT_RUN_MSGS_FULL, &pval, &t_evol);

      g_assert_nonnull (pval);
      g_assert_nonnull (t_evol);
      g_assert_cmpuint (ncm_vector_len (pval), ==, ncm_matrix_ncols (t_evol));

      for (i = 0; i < ncm_vector_len (pval); i++)
      {
        const gdouble pval_i = ncm_vector_get (pval, i);

        g_assert_true (gsl_finite (pval_i));
      }

      for (i = 0; i < ncm_matrix_nrows (t_evol); i++)
      {
        guint j;

        for (j = 0; j < ncm_matrix_ncols (t_evol); j++)
        {
          const gdouble tval_ij = ncm_matrix_get (t_evol, i, j);

          g_assert_cmpfloat (tval_ij, >=, 0.0);
        }
      }

      ncm_vector_clear (&pval);
      ncm_matrix_clear (&t_evol);
    }

    return;
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*Calculating evolution to time*");
}

void
test_ncm_mset_catalog_traps (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/invalid/run/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_mset_catalog_invalid_run (TestNcmMSetCatalog *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

#ifdef HAVE_CFITSIO

/*
 * Builds a small free-parameter mset, writes a new-style catalog file for it
 * (mset embedded in HDU0), and returns both. The catalog is not left open.
 */
static NcmMSet *
_test_ncm_mset_catalog_new_file (const gchar *filename)
{
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (3);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmMSetCatalog *mcat;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE, "m2lnL", "-2\\ln(L)", NULL);
  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  ncm_mset_catalog_set_run_type (mcat, "test-run");
  ncm_mset_catalog_set_file (mcat, filename);

  {
    NcmVector *x  = ncm_vector_new (ncm_mset_fparams_len (mset));
    gdouble ax[1] = { 12.3 };

    ncm_mset_fparams_get_vector (mset, x);
    ncm_mset_catalog_add_from_vector_array (mcat, x, ax);
    ncm_mset_catalog_sync (mcat, TRUE);

    ncm_vector_clear (&x);
  }

  ncm_mset_catalog_clear (&mcat);
  ncm_model_mvnd_clear (&model_mvnd);

  return mset;
}

/*
 * Rewrites HDU0 in the pre-vardict format (a bare serialized #NcmMSet
 * object, format label "gvariant", no functions entry), simulating a
 * catalog file written before the vardict envelope was introduced, to
 * check that it is still read correctly.
 */
static void
_test_ncm_mset_catalog_rewrite_hdu0_legacy_object (const gchar *filename, NcmMSet *mset)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  GVariant *var     = ncm_serialize_to_variant (ser, G_OBJECT (mset));
  glong naxes[1]    = { (glong) g_variant_get_size (var) };
  fitsfile *fptr    = NULL;
  gint hdutype      = 0;
  gint status       = 0;

  fits_open_file (&fptr, filename, READWRITE, &status);
  g_assert_cmpint (status, ==, 0);

  fits_movabs_hdu (fptr, 1, &hdutype, &status);
  g_assert_cmpint (status, ==, 0);

  fits_resize_img (fptr, BYTE_IMG, 1, naxes, &status);
  g_assert_cmpint (status, ==, 0);

  fits_write_img (fptr, TBYTE, 1, naxes[0], (gpointer) g_variant_get_data (var), &status);
  g_assert_cmpint (status, ==, 0);

  fits_update_key_str (fptr, NCM_MSET_CATALOG_MSET_FORMAT_LABEL, NCM_MSET_CATALOG_MSET_FORMAT_OBJECT,
                       "Format of the data stored in this HDU.", &status);
  g_assert_cmpint (status, ==, 0);

  fits_close_file (fptr, &status);
  g_assert_cmpint (status, ==, 0);

  g_variant_unref (var);
  ncm_serialize_free (ser);
}

/*
 * Strips the mset embedded by _test_ncm_mset_catalog_new_file() back to an
 * empty primary HDU, simulating a catalog file written before HDU0 embedding
 * was introduced.
 */
static void
_test_ncm_mset_catalog_strip_hdu0 (const gchar *filename)
{
  fitsfile *fptr = NULL;
  gint hdutype   = 0;
  gint status    = 0;

  fits_open_file (&fptr, filename, READWRITE, &status);
  g_assert_cmpint (status, ==, 0);

  fits_movabs_hdu (fptr, 1, &hdutype, &status);
  g_assert_cmpint (status, ==, 0);

  fits_resize_img (fptr, 8, 0, NULL, &status);
  g_assert_cmpint (status, ==, 0);

  fits_close_file (fptr, &status);
  g_assert_cmpint (status, ==, 0);
}

void
test_ncm_mset_catalog_file_hdu0_roundtrip (void)
{
  gchar *tmp_dir      = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_hdu0_XXXXXX", NULL);
  gchar *filename     = g_strdup_printf ("%s/cat.fits", tmp_dir);
  gchar *mset_sidecar = g_strdup_printf ("%s/cat.mset", tmp_dir);
  NcmMSet *mset       = _test_ncm_mset_catalog_new_file (filename);
  NcmMSetCatalog *mcat2;
  NcmMSet *mset2;

  g_assert_false (g_file_test (mset_sidecar, G_FILE_TEST_EXISTS));

  mcat2 = ncm_mset_catalog_new_from_file_ro (filename, 0);
  mset2 = ncm_mset_catalog_peek_mset (mcat2);

  g_assert_true (ncm_mset_cmp (mset, mset2, TRUE));

  ncm_mset_catalog_clear (&mcat2);
  ncm_mset_clear (&mset);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (mset_sidecar);
  g_free (tmp_dir);
}

void
test_ncm_mset_catalog_file_functions_array_roundtrip (void)
{
  gchar *tmp_dir           = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_functions_XXXXXX", NULL);
  gchar *filename          = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (3);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmObjArray *functions   = ncm_obj_array_new ();
  NcmVector *v             = ncm_vector_new (2);
  NcmMSetCatalog *mcat;
  NcmMSetCatalog *mcat2;
  NcmObjArray *functions2;

  ncm_vector_set_all (v, 7.0);
  ncm_obj_array_add (functions, G_OBJECT (v));

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE, "m2lnL", "-2\\ln(L)", NULL);
  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  ncm_mset_catalog_set_run_type (mcat, "test-run");
  ncm_mset_catalog_set_functions_array (mcat, functions);
  ncm_mset_catalog_set_file (mcat, filename);

  g_assert_true (ncm_mset_catalog_peek_functions_array (mcat) == functions);

  {
    NcmVector *x  = ncm_vector_new (ncm_mset_fparams_len (mset));
    gdouble ax[1] = { 12.3 };

    ncm_mset_fparams_get_vector (mset, x);
    ncm_mset_catalog_add_from_vector_array (mcat, x, ax);
    ncm_mset_catalog_sync (mcat, TRUE);

    ncm_vector_clear (&x);
  }

  ncm_mset_catalog_clear (&mcat);

  mcat2      = ncm_mset_catalog_new_from_file_ro (filename, 0);
  functions2 = ncm_mset_catalog_peek_functions_array (mcat2);

  g_assert_nonnull (functions2);
  g_assert_cmpuint (ncm_obj_array_len (functions2), ==, 1);
  g_assert_true (NCM_IS_VECTOR (ncm_obj_array_peek (functions2, 0)));
  g_assert_cmpfloat (ncm_vector_get (NCM_VECTOR (ncm_obj_array_peek (functions2, 0)), 0), ==, 7.0);
  g_assert_true (ncm_obj_array_peek (functions2, 0) != G_OBJECT (v));

  ncm_mset_catalog_clear (&mcat2);
  ncm_obj_array_unref (functions);
  ncm_vector_free (v);
  ncm_mset_clear (&mset);
  ncm_model_mvnd_clear (&model_mvnd);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);
}

void
test_ncm_mset_catalog_file_hdu0_legacy_object_format (void)
{
  gchar *tmp_dir  = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_legacy_obj_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmMSet *mset   = _test_ncm_mset_catalog_new_file (filename);
  NcmMSetCatalog *mcat2;
  NcmMSet *mset2;

  _test_ncm_mset_catalog_rewrite_hdu0_legacy_object (filename, mset);

  mcat2 = ncm_mset_catalog_new_from_file_ro (filename, 0);
  mset2 = ncm_mset_catalog_peek_mset (mcat2);

  g_assert_true (ncm_mset_cmp (mset, mset2, TRUE));
  g_assert_null (ncm_mset_catalog_peek_functions_array (mcat2));

  ncm_mset_catalog_clear (&mcat2);
  ncm_mset_clear (&mset);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);
}

void
test_ncm_mset_catalog_file_legacy_fallback (void)
{
  gchar *tmp_dir      = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_legacy_XXXXXX", NULL);
  gchar *filename     = g_strdup_printf ("%s/cat.fits", tmp_dir);
  gchar *mset_sidecar = g_strdup_printf ("%s/cat.mset", tmp_dir);
  NcmMSet *mset       = _test_ncm_mset_catalog_new_file (filename);
  NcmMSetCatalog *mcat2;
  NcmMSet *mset2;
  NcmSerialize *ser;

  _test_ncm_mset_catalog_strip_hdu0 (filename);

  ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  ncm_mset_save (mset, ser, mset_sidecar, TRUE, NULL);
  ncm_serialize_clear (&ser);

  mcat2 = ncm_mset_catalog_new_from_file_ro (filename, 0);
  mset2 = ncm_mset_catalog_peek_mset (mcat2);

  g_assert_true (ncm_mset_cmp (mset, mset2, TRUE));

  ncm_mset_catalog_clear (&mcat2);
  ncm_mset_clear (&mset);

  g_unlink (mset_sidecar);
  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (mset_sidecar);
  g_free (tmp_dir);
}

void
test_ncm_mset_catalog_file_missing_mset_traps (void)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/file/missing_mset/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_mset_catalog_file_missing_mset_subprocess (void)
{
  gchar *tmp_dir  = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_missing_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmMSet *mset   = _test_ncm_mset_catalog_new_file (filename);

  _test_ncm_mset_catalog_strip_hdu0 (filename);
  ncm_mset_clear (&mset);

  /* Neither an embedded mset (just stripped) nor a `.mset' sidecar (never
   * written) exists: this must abort with a fatal error. */
  ncm_mset_catalog_new_from_file_ro (filename, 0);

  g_assert_not_reached ();
}

void
test_ncm_mset_catalog_file_peek_info (void)
{
  gchar *tmp_dir  = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_peek_info_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmMSet *mset   = _test_ncm_mset_catalog_new_file (filename);
  glong nrows     = -1;
  guint nchains   = 0;
  gint first_id   = -1;

  ncm_mset_catalog_peek_info_from_file (filename, &nrows, &nchains, &first_id);

  g_assert_cmpint (nrows, ==, 1);
  g_assert_cmpuint (nchains, ==, 1);
  g_assert_cmpint (first_id, ==, 0);

  ncm_mset_clear (&mset);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);
}

void
test_ncm_mset_catalog_file_burnin_exceeds_traps (void)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/file/burnin_exceeds/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exceeds catalog*");
}

void
test_ncm_mset_catalog_file_burnin_exceeds_subprocess (void)
{
  gchar *tmp_dir  = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_burnin_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmMSet *mset   = _test_ncm_mset_catalog_new_file (filename);

  ncm_mset_clear (&mset);

  /* The file has a single row: any burnin > 1 must abort with a clear,
   * unit-labeled error message (see _ncm_mset_catalog_open_create_file). */
  ncm_mset_catalog_new_from_file_ro (filename, 2);

  g_assert_not_reached ();
}

#endif /* HAVE_CFITSIO */


/*
 * A multi-chain catalog written to a file and read back.
 *
 * Rows are inserted directly, so nothing has to sample: what is under test is the
 * bookkeeping around them. Two things need more than the single-chain single-row
 * catalog the other file tests use -- the Gelman-Rubin shrink factor returns 1 outright
 * for one chain, and the RNG state is only read back when a file already carries it.
 */
#define TEST_CAT_NCHAINS 4
#define TEST_CAT_NROWS_PER_CHAIN 25
#define TEST_CAT_DIM 3

void
test_ncm_mset_catalog_file_multichain (void)
{
  gchar *tmp_dir           = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_multichain_XXXXXX", NULL);
  gchar *filename          = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_CAT_DIM);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);

  /* Two generators on purpose: the catalog's is the sampler's, and its recorded state
   * must not move, so the synthetic rows are drawn from a separate one. */
  NcmRNG *cat_rng  = ncm_rng_seeded_new (NULL, 987654321);
  NcmRNG *data_rng = ncm_rng_seeded_new (NULL, 13579);
  NcmMSetCatalog *mcat;
  gchar *state_at_set;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, TEST_CAT_NCHAINS, FALSE, "m2lnL", "-2\\ln(L)", NULL);

  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  ncm_mset_catalog_set_run_type (mcat, "multichain-run");
  ncm_mset_catalog_set_rng (mcat, cat_rng);
  ncm_mset_catalog_set_file (mcat, filename);

  /* The state the catalog recorded, which is the one the file carries: the generator
   * itself moves on as the rows below are drawn. */
  state_at_set = ncm_rng_get_state (cat_rng);

  {
    NcmVector *x = ncm_vector_new (ncm_mset_fparams_len (mset));

    for (i = 0; i < TEST_CAT_NROWS_PER_CHAIN * TEST_CAT_NCHAINS; i++)
    {
      gdouble ax[1] = { 1.0 + 0.01 * i };

      /* Chains that differ a little, so the between-chain covariance is not degenerate
       * and the shrink factor has something to measure. */
      for (j = 0; j < TEST_CAT_DIM; j++)
        ncm_vector_set (x, j, ncm_rng_gaussian_gen (data_rng, 0.1 * (i % TEST_CAT_NCHAINS), 1.0));

      ncm_mset_catalog_add_from_vector_array (mcat, x, ax);
    }

    ncm_vector_clear (&x);
  }

  ncm_mset_catalog_sync (mcat, TRUE);

  /* Every column answers to its own name, and the lookup is the inverse of the listing. */
  g_assert_cmpuint (ncm_mset_catalog_ncols (mcat), >, 0);

  for (i = 0; i < ncm_mset_catalog_ncols (mcat); i++)
  {
    const gchar *full = ncm_mset_catalog_col_full_name (mcat, i);
    guint back        = G_MAXUINT;

    g_assert_nonnull (full);
    g_assert_true (ncm_mset_catalog_col_by_name (mcat, full, &back));
    g_assert_cmpuint (back, ==, i);
  }

  {
    guint missing = G_MAXUINT;

    g_assert_false (ncm_mset_catalog_col_by_name (mcat, "no-such-column", &missing));
  }

  g_assert_cmpstr (ncm_mset_catalog_get_run_type (mcat), ==, "multichain-run");
  g_assert_nonnull (ncm_mset_catalog_peek_pstats (mcat));
  g_assert_nonnull (ncm_mset_catalog_peek_e_mean_stats (mcat));
  g_assert_cmpuint (ncm_mset_catalog_max_time (mcat), ==, TEST_CAT_NROWS_PER_CHAIN);

  /* Gelman-Rubin: at least 1 by construction, and finite for chains this similar. */
  {
    const gdouble shrink = ncm_mset_catalog_get_shrink_factor (mcat);

    g_assert_true (gsl_finite (shrink));
    g_assert_cmpfloat (shrink, >=, 1.0);
  }

  ncm_mset_catalog_clear (&mcat);

  /* The reopen path -- a catalog built with no generator adopting the one the file
   * carries, which is what lets a run resume rather than re-draw -- is not exercised
   * here: getting the recorded state to survive the write/reopen sequence needs
   * knowledge of the ordering this test could not establish from outside.
   */

  g_free (state_at_set);
  ncm_rng_free (cat_rng);
  ncm_rng_free (data_rng);
  ncm_mset_clear (&mset);
  ncm_model_mvnd_clear (&model_mvnd);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);
}

#ifdef HAVE_CFITSIO

/*
 * Writes TEST_CAT_NCHAINS chains x nrows_per_chain rows into a new catalog file, setting
 * the Markovian id when that row is reached; returns the model set (the catalog is closed).
 */
static NcmMSet *
_test_ncm_mset_catalog_new_markovian_file (const gchar *filename, const guint nrows_per_chain, const gint markovian_id)
{
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_CAT_DIM);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmRNG *rng              = ncm_rng_seeded_new (NULL, 2468);
  NcmMSetCatalog *mcat;
  NcmVector *x;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, TEST_CAT_NCHAINS, FALSE, "m2lnL", "-2\\ln(L)", NULL);
  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  ncm_mset_catalog_set_run_type (mcat, "markovian-run");
  ncm_mset_catalog_set_file (mcat, filename);

  /* Fresh catalog: every row is Markovian until told otherwise. */
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, 0);
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat), ==, 0);

  x = ncm_vector_new (ncm_mset_fparams_len (mset));

  for (i = 0; i < nrows_per_chain * TEST_CAT_NCHAINS; i++)
  {
    gdouble ax[1] = { 1.0 + 0.01 * i };

    for (j = 0; j < TEST_CAT_DIM; j++)
      ncm_vector_set (x, j, ncm_rng_gaussian_gen (rng, 0.0, 1.0));

    ncm_mset_catalog_add_from_vector_array (mcat, x, ax);

    if ((gint) i == markovian_id - 1)
      ncm_mset_catalog_set_markovian_id (mcat, markovian_id);
  }

  ncm_mset_catalog_sync (mcat, TRUE);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, markovian_id);

  ncm_vector_clear (&x);
  ncm_rng_free (rng);
  ncm_mset_catalog_clear (&mcat);
  ncm_model_mvnd_clear (&model_mvnd);

  return mset;
}

/* The id survives the file round trip, the file-level peek sees it, the conversion to
 * iterations is exact, and reading with a burn-in shifts it like every other row id. */
void
test_ncm_mset_catalog_file_markovian_id_roundtrip (void)
{
  gchar *tmp_dir          = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_markid_XXXXXX", NULL);
  gchar *filename         = g_strdup_printf ("%s/cat.fits", tmp_dir);
  const gint markovian_id = 3 * TEST_CAT_NCHAINS; /* the fourth ensemble starts the chain */
  NcmMSet *mset           = _test_ncm_mset_catalog_new_markovian_file (filename, 6, markovian_id);
  NcmMSetCatalog *mcat_ro = ncm_mset_catalog_new_from_file_ro (filename, 0);
  NcmMSetCatalog *mcat_b  = ncm_mset_catalog_new_from_file_ro (filename, 2 * TEST_CAT_NCHAINS);
  NcmMSetCatalog *mcat_b2 = ncm_mset_catalog_new_from_file_ro (filename, 5 * TEST_CAT_NCHAINS);
  gint prop_id            = -1;

  g_assert_cmpint (ncm_mset_catalog_peek_markovian_id_from_file (filename), ==, markovian_id);

  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat_ro), ==, markovian_id);
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat_ro), ==, 3);
  g_object_get (G_OBJECT (mcat_ro), "markovian-id", &prop_id, NULL);
  g_assert_cmpint (prop_id, ==, markovian_id);

  /* Two ensembles dropped at load: the chain now starts one ensemble in. */
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat_b), ==, markovian_id - 2 * TEST_CAT_NCHAINS);
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat_b), ==, 1);

  /* Burn-in past the id: everything loaded is Markovian, the id is the first row. */
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat_b2), ==, ncm_mset_catalog_get_first_id (mcat_b2));
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat_b2), ==, 0);

  ncm_mset_catalog_clear (&mcat_ro);
  ncm_mset_catalog_clear (&mcat_b);
  ncm_mset_catalog_clear (&mcat_b2);
  ncm_mset_clear (&mset);
  g_unlink (filename);
  g_rmdir (tmp_dir);
  g_free (filename);
  g_free (tmp_dir);
}

/* A file written before the key existed reads as fully Markovian. */
void
test_ncm_mset_catalog_file_markovian_id_missing_key (void)
{
  gchar *tmp_dir  = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_markid_old_XXXXXX", NULL);
  gchar *filename = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmMSet *mset   = _test_ncm_mset_catalog_new_markovian_file (filename, 4, 2 * TEST_CAT_NCHAINS);
  fitsfile *fptr  = NULL;
  gint status     = 0;
  NcmMSetCatalog *mcat;

  fits_open_file (&fptr, filename, READWRITE, &status);
  g_assert_cmpint (status, ==, 0);
  fits_movnam_hdu (fptr, BINARY_TBL, NCM_MSET_CATALOG_EXTNAME, 0, &status);
  g_assert_cmpint (status, ==, 0);
  fits_delete_key (fptr, NCM_MSET_CATALOG_MARKOVIAN_ID_LABEL, &status);
  g_assert_cmpint (status, ==, 0);
  fits_close_file (fptr, &status);
  g_assert_cmpint (status, ==, 0);

  g_assert_cmpint (ncm_mset_catalog_peek_markovian_id_from_file (filename), ==, 0);

  mcat = ncm_mset_catalog_new_from_file_ro (filename, 0);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, ncm_mset_catalog_get_first_id (mcat));
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat), ==, 0);

  ncm_mset_catalog_clear (&mcat);
  ncm_mset_clear (&mset);
  g_unlink (filename);
  g_rmdir (tmp_dir);
  g_free (filename);
  g_free (tmp_dir);
}

void
test_ncm_mset_catalog_file_markovian_id_backwards_traps (void)
{
  g_test_trap_subprocess ("/ncm/mset/catalog/file/markovian_id/backwards/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*cannot move backwards*");
}

void
test_ncm_mset_catalog_file_markovian_id_backwards_subprocess (void)
{
  gchar *tmp_dir       = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_markid_back_XXXXXX", NULL);
  gchar *filename      = g_strdup_printf ("%s/cat.fits", tmp_dir);
  NcmMSet *mset        = _test_ncm_mset_catalog_new_markovian_file (filename, 4, 2 * TEST_CAT_NCHAINS);
  NcmMSetCatalog *mcat = ncm_mset_catalog_new_from_file (filename, 0);

  ncm_mset_catalog_set_markovian_id (mcat, TEST_CAT_NCHAINS);

  ncm_mset_catalog_clear (&mcat);
  ncm_mset_clear (&mset);
  g_free (filename);
  g_free (tmp_dir);
}

/* trim, thinning and remove_last_ensemble keep the Markovian id on the same rows. */
void
test_ncm_mset_catalog_file_markovian_id_trim (void)
{
  gchar *tmp_dir          = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_markid_trim_XXXXXX", NULL);
  gchar *filename         = g_strdup_printf ("%s/cat.fits", tmp_dir);
  const gint markovian_id = 4 * TEST_CAT_NCHAINS; /* iterations 0-3 non-Markovian, 4-9 Markovian */
  NcmMSet *mset           = _test_ncm_mset_catalog_new_markovian_file (filename, 10, markovian_id);
  NcmMSetCatalog *mcat    = ncm_mset_catalog_new_from_file (filename, 0);

  /* Cut 2 of the 4 non-Markovian iterations: 2 remain before the chain. */
  ncm_mset_catalog_trim (mcat, 2, 1);
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat), ==, 2);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, ncm_mset_catalog_get_first_id (mcat) + 2 * TEST_CAT_NCHAINS);

  /* Thin by 2 without a cut: the 2 remaining non-Markovian iterations become 1. */
  ncm_mset_catalog_trim (mcat, 0, 2);
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat), ==, 1);

  /* Cut past the chain start: every row left is Markovian. */
  ncm_mset_catalog_trim (mcat, 3, 1);
  g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat), ==, 0);
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, ncm_mset_catalog_get_first_id (mcat));

  ncm_mset_catalog_clear (&mcat);
  ncm_mset_clear (&mset);

  /* Interrupted phase: the id sits after the last row; removing the last ensemble moves it
   * back to the row after the new last one, so the phase stays open. A fresh file: the
   * trims above rewrote the first one in place (the original is kept as .bak). */
  {
    gchar *filename2 = g_strdup_printf ("%s/cat2.fits", tmp_dir);

    mset = _test_ncm_mset_catalog_new_markovian_file (filename2, 10, markovian_id);
    mcat = ncm_mset_catalog_new_from_file (filename2, 0);
    ncm_mset_catalog_set_markovian_id (mcat, ncm_mset_catalog_get_cur_id (mcat) + 1);
    ncm_mset_catalog_remove_last_ensemble (mcat);
    g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==, ncm_mset_catalog_get_cur_id (mcat) + 1);
    g_assert_cmpuint (ncm_mset_catalog_get_markovian_burnin (mcat), ==, 9);

    ncm_mset_catalog_clear (&mcat);
    ncm_mset_clear (&mset);
    g_unlink (filename2);
    g_free (filename2);
  }

  {
    /* The trims and the removal renamed the originals to .N.bak; remove them too. */
    GDir *dir         = g_dir_open (tmp_dir, 0, NULL);
    const gchar *name = NULL;

    while ((name = g_dir_read_name (dir)) != NULL)
    {
      gchar *path = g_build_filename (tmp_dir, name, NULL);

      g_unlink (path);
      g_free (path);
    }

    g_dir_close (dir);
  }

  g_rmdir (tmp_dir);
  g_free (filename);
  g_free (tmp_dir);
}

/* Walkers frozen at distinct positions: the rows vary, the ensemble mean does not. The
 * effective sample size is then not a sample size, and the gate must ask for more. */
void
test_ncm_mset_catalog_file_markovian_id_frozen_needs_more (void)
{
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_CAT_DIM);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmMSetCatalog *mcat;
  NcmVector *x;
  guint required = 0;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, TEST_CAT_NCHAINS, FALSE, "m2lnL", "-2\\ln(L)", NULL);
  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  ncm_mset_catalog_set_run_type (mcat, "frozen-run");

  x = ncm_vector_new (ncm_mset_fparams_len (mset));

  for (i = 0; i < 200 * TEST_CAT_NCHAINS; i++)
  {
    gdouble ax[1] = { 1.0 + (i % TEST_CAT_NCHAINS) };

    for (j = 0; j < TEST_CAT_DIM; j++)
      ncm_vector_set (x, j, 1.0 * (i % TEST_CAT_NCHAINS) + 0.1 * j);

    ncm_mset_catalog_add_from_vector_array (mcat, x, ax);
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (mcat, FALSE);
  g_assert_true (ncm_mset_catalog_tau_needs_more (mcat, &required));

  ncm_vector_clear (&x);
  ncm_mset_catalog_clear (&mcat);
  ncm_mset_clear (&mset);
  ncm_model_mvnd_clear (&model_mvnd);
}

#endif /* HAVE_CFITSIO */

void
test_ncm_mset_catalog_tau_frozen_walkers (void)
{
  NcmRNG *rng              = ncm_rng_seeded_new (NULL, 271828);
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (2);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  const gchar *names[]     = {"m2lnL", NULL};
  const gchar *symbols[]   = {"-2\\ln(L)", NULL};
  const guint nchains      = 4;
  const guint nitens       = 200;
  NcmMSetCatalog *mcat;
  NcmVector *row;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new_array (mset, 1, nchains, FALSE, (gchar **) names, (gchar **) symbols);
  row  = ncm_vector_new (ncm_mset_catalog_ncols (mcat));

  /*
   * Walkers that all but stopped moving. Unlike the perfectly frozen catalog of
   * markovian_id/frozen_needs_more, each one still jitters, so the variance of the
   * ensemble mean is small rather than zero and the state shows up as a K_eff far above
   * the number of chains instead of as a refused variance. The last parameter never moves
   * at all, which is what the report has to name.
   */
  for (i = 0; i < nitens; i++)
  {
    for (j = 0; j < nchains; j++)
    {
      ncm_vector_set (row, 0, ncm_rng_ugaussian_gen (rng));
      ncm_vector_set (row, 1, 1.0 * j + 1.0e-8 * ncm_rng_ugaussian_gen (rng));
      ncm_vector_set (row, 2, 1.0);

      ncm_mset_catalog_add_from_vector (mcat, row);
    }
  }

  ncm_mset_catalog_estimate_autocorrelation_tau (mcat, FALSE);

  /* The frozen walkers are seen as a K_eff far above the number of chains: the rows vary,
   * the ensemble mean does not. */
  g_assert_cmpfloat (ncm_mset_catalog_get_keff (mcat, 1), >, 2.0 * nchains);

  /* The parameter that never moved carries the zero-variance condition. */
  g_assert_cmpuint (ncm_mset_catalog_get_tau_diag (mcat, 2) & NCM_STATS_ACORR_DIAG_ZERO_VARIANCE, !=, 0);

  g_assert_true (ncm_mset_catalog_log_tau_diag (mcat));
  g_assert_true (ncm_mset_catalog_tau_needs_more (mcat, NULL));

  /* Resetting drops everything accumulated from the rows. */
  ncm_mset_catalog_reset_stats (mcat);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, 0);

  ncm_vector_free (row);
  ncm_mset_catalog_free (mcat);
  ncm_mset_clear (&mset);
  ncm_model_mvnd_clear (&model_mvnd);
  ncm_rng_free (rng);
}

void
test_ncm_mset_catalog_trim_oob_markovian (void)
{
  gchar *tmp_dir           = g_dir_make_tmp ("tmp_test_ncm_mset_catalog_trim_oob_XXXXXX", NULL);
  gchar *out_file          = g_strdup_printf ("%s/cat_oob.fits", tmp_dir);
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_CAT_DIM);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  const guint nitens       = 10;
  const guint mark_itens   = 4;
  NcmMSetCatalog *mcat;
  NcmVector *x;
  guint i, j, ndel;
  guint oob = 0, oob_before = 0;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, TEST_CAT_NCHAINS, FALSE, "m2lnL", "-2\\ln(L)", NULL);
  ncm_mset_catalog_set_run_type (mcat, "trim-oob-run");

  x = ncm_vector_new (ncm_mset_fparams_len (mset));

  /*
   * Every fifth row sits outside the parameter bounds, so rows are dropped from both
   * sides of the Markovian id. Trimming renumbers what is left, and the id has to follow
   * the rows it marks rather than stay at its old number.
   */
  for (i = 0; i < nitens * TEST_CAT_NCHAINS; i++)
  {
    const gboolean row_oob = (i % 5 == 0);
    gdouble ax[1]          = { 1.0 * i };

    for (j = 0; j < TEST_CAT_DIM; j++)
      ncm_vector_set (x, j, row_oob ? 100.0 : 0.1 * j);

    ncm_mset_catalog_add_from_vector_array (mcat, x, ax);

    if (row_oob)
    {
      oob++;

      if (i < mark_itens * TEST_CAT_NCHAINS)
        oob_before++;
    }
  }

  ncm_mset_catalog_set_markovian_id (mcat, ncm_mset_catalog_get_first_id (mcat) + mark_itens * TEST_CAT_NCHAINS);

  ndel = ncm_mset_catalog_trim_oob (mcat, out_file);

  g_assert_cmpuint (ndel, ==, oob);
  g_assert_cmpuint (ncm_mset_catalog_len (mcat), ==, nitens * TEST_CAT_NCHAINS - oob);

  /* Out-of-bounds rows before the chain start are gone, so the start moves back by that many. */
  g_assert_cmpint (ncm_mset_catalog_get_markovian_id (mcat), ==,
                   ncm_mset_catalog_get_first_id (mcat) + mark_itens * TEST_CAT_NCHAINS - oob_before);

  /* Trimming always leaves a single chain. */
  g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, 1);

  ncm_vector_clear (&x);
  ncm_mset_catalog_clear (&mcat);
  ncm_mset_clear (&mset);
  ncm_model_mvnd_clear (&model_mvnd);

  {
    GDir *dir         = g_dir_open (tmp_dir, 0, NULL);
    const gchar *name = NULL;

    while ((name = g_dir_read_name (dir)) != NULL)
    {
      gchar *path = g_build_filename (tmp_dir, name, NULL);

      g_unlink (path);
      g_free (path);
    }

    g_dir_close (dir);
  }

  g_rmdir (tmp_dir);
  g_free (out_file);
  g_free (tmp_dir);
}


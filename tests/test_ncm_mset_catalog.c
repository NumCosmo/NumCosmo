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
#include <numcosmo/numcosmo.h>

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
void test_ncm_mset_catalog_accept_ratio_array (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_calc_param_ensemble_evol (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_calc_add_param_ensemble_evol (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_calc_add_param_ensemble_evol_short (TestNcmMSetCatalog *test, gconstpointer pdata);
void test_ncm_mset_catalog_invalid_run (TestNcmMSetCatalog *test, gconstpointer pdata);

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


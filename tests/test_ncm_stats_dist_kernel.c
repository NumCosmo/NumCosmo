/***************************************************************************
 *            test_ncm_stats_dist.c
 *
 *  Wed November 07 17:57:28 2018
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

typedef enum _NcmStatsDistKernelType
{
  NCM_STATS_DIST_KERNEL_GAUSS,
  NCM_STATS_DIST_KERNEL_ST3,
} NcmStatsDistKernelType;

typedef struct _TestNcmStatsDistKernel
{
  NcmStatsDistKernel *kernel;
  NcmMSet *mset;
  NcmStatsDistKernelType kernel_type;
  gdouble nu;
  guint dim;
  guint nfail;
} TestNcmStatsDistKernel;


static void test_ncm_stats_dist_kernel_new_st (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_new_gauss (TestNcmStatsDistKernel *test, gconstpointer pdata);

static void test_ncm_stats_dist_kernel_dim (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_bandwidth (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_norm (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_sum (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_sample (TestNcmStatsDistKernel *test, gconstpointer pdata);

static void test_ncm_stats_dist_kernel_free (TestNcmStatsDistKernel *test, gconstpointer pdata);

static void test_ncm_stats_dist_kernel_traps (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_invalid_stub (TestNcmStatsDistKernel *test, gconstpointer pdata);

typedef struct _TestNcmStatsDistKernelFunc
{
  const gchar *name;

  void (*test_func) (TestNcmStatsDistKernel *test, gconstpointer pdata);
} TestNcmStatsDistKernelFunc;

#define TEST_NCM_STATS_DIST_KERNEL_CONSTRUCTORS_LEN 2
#define TEST_NCM_STATS_DIST_KERNEL_TESTS_LEN 5

static TestNcmStatsDistKernelFunc constructors[TEST_NCM_STATS_DIST_KERNEL_CONSTRUCTORS_LEN] = {
  {"gauss", &test_ncm_stats_dist_kernel_new_gauss},
  {"st",  &test_ncm_stats_dist_kernel_new_st}
};

static TestNcmStatsDistKernelFunc tests[TEST_NCM_STATS_DIST_KERNEL_TESTS_LEN] = {
  {"dim",   &test_ncm_stats_dist_kernel_dim},
  {"band",        &test_ncm_stats_dist_kernel_bandwidth},
  {"gauss/norm",        &test_ncm_stats_dist_kernel_norm},
  {"sum",         &test_ncm_stats_dist_kernel_sum},
  {"sample",      &test_ncm_stats_dist_kernel_sample},
};

gint
main (gint argc, gchar *argv[])
{
  gint i, j;

  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < TEST_NCM_STATS_DIST_KERNEL_CONSTRUCTORS_LEN; i++)
  {
    for (j = 0; j < TEST_NCM_STATS_DIST_KERNEL_TESTS_LEN; j++)
    {
      gchar *test_name = g_strdup_printf ("/ncm/stats/dist/kernel/%s/%s", constructors[i].name, tests[j].name);

      g_test_add (test_name, TestNcmStatsDistKernel, NULL, constructors[i].test_func, tests[j].test_func, &test_ncm_stats_dist_kernel_free);
      g_free (test_name);
    }
  }

  g_test_add ("/ncm/stats/dist/kernel/gauss/traps", TestNcmStatsDistKernel, NULL,
              &test_ncm_stats_dist_kernel_new_gauss,
              &test_ncm_stats_dist_kernel_traps,
              &test_ncm_stats_dist_kernel_free);

#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/stats/dist/kernel/gauss/invalid/stub/subprocess", TestNcmStatsDistKernel, NULL,
              &test_ncm_stats_dist_kernel_new_gauss,
              &test_ncm_stats_dist_kernel_invalid_stub,
              &test_ncm_stats_dist_kernel_free);
#endif


  g_test_run ();
}

static void
test_ncm_stats_dist_kernel_new_gauss (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  const guint dim                    = g_test_rand_int_range (2, 4);
  NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);

  test->dim         = dim;
  test->kernel      = NCM_STATS_DIST_KERNEL (sdk_gauss);
  test->kernel_type = NCM_STATS_DIST_KERNEL_GAUSS;
  test->nfail       = 0;

  ncm_stats_dist_kernel_gauss_ref (sdk_gauss);
  ncm_stats_dist_kernel_gauss_free (sdk_gauss);
  {
    NcmStatsDistKernelGauss *sdk_gauss0 = ncm_stats_dist_kernel_gauss_ref (sdk_gauss);

    ncm_stats_dist_kernel_gauss_clear (&sdk_gauss0);
    g_assert_true (sdk_gauss0 == NULL);
  }
}

static void
test_ncm_stats_dist_kernel_new_st (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  const gdouble nu             = g_test_rand_double_range (3.0, 5.0);
  const guint dim              = g_test_rand_int_range (2, 4);
  NcmStatsDistKernelST *sdk_st = ncm_stats_dist_kernel_st_new (dim, nu);

  test->dim         = dim;
  test->nu          = nu;
  test->kernel_type = NCM_STATS_DIST_KERNEL_ST3;
  test->kernel      = NCM_STATS_DIST_KERNEL (sdk_st);
  test->nfail       = 0;

  ncm_stats_dist_kernel_st_ref (sdk_st);
  ncm_stats_dist_kernel_st_free (sdk_st);
  {
    NcmStatsDistKernelST *sdk_st0 = ncm_stats_dist_kernel_st_ref (sdk_st);

    ncm_stats_dist_kernel_st_clear (&sdk_st0);
    g_assert_true (sdk_st0 == NULL);
  }
}

static void
test_ncm_stats_dist_kernel_dim (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  g_assert_true (test->dim == ncm_stats_dist_kernel_get_dim (test->kernel));
}

static void
test_ncm_stats_dist_kernel_bandwidth (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  const gdouble n = g_test_rand_double_range (1.0, 200.0);
  const gdouble h = ncm_stats_dist_kernel_get_rot_bandwidth (test->kernel, n);
  gdouble h_test  = 0.0;

  switch (test->kernel_type)
  {
    case NCM_STATS_DIST_KERNEL_GAUSS:
      h_test = pow (4.0 / (n * (test->dim + 2.0)), 1.0 / (test->dim + 4.0));
      break;
    case NCM_STATS_DIST_KERNEL_ST3:
      h_test = pow (16.0 * gsl_pow_2 (test->nu - 2) * (1.0 + test->dim + test->nu) * (3.0 + test->dim + test->nu) / ((2.0 + test->dim) * (test->dim + test->nu) * (2.0 + test->dim + test->nu) * (test->dim + 2.0 * test->nu) * (2.0 + test->dim + 2.0 * test->nu) * n), 1.0 / (test->dim + 4.0));
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  ncm_assert_cmpdouble_e (h_test, ==, h, 1.0e-14, 0.0);
}

#define TESTMULT 200

static void
test_ncm_stats_dist_kernel_norm (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint np                 = TESTMULT * test->dim;
  NcmVector *m2lnp_v             = ncm_vector_new (np);
  const guint ntests             = 100 * g_test_rand_int_range (1, 5);
  guint i;
  NcmMatrix *cov    = NULL;
  gdouble lndet_cov = 0.0;

  ncm_data_gauss_cov_use_norma (NCM_DATA_GAUSS_COV (data_mvnd), FALSE);
  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < np; i++)
  {
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  cov       = ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd));
  lndet_cov = ncm_matrix_cholesky_lndet (cov);


  switch (test->kernel_type)
  {
    case NCM_STATS_DIST_KERNEL_GAUSS:
    {
      gdouble norm_test    = 0.5 * (test->dim * ncm_c_ln2pi () + lndet_cov);
      const gdouble lnnorm = ncm_stats_dist_kernel_get_lnnorm (test->kernel, cov);

      ncm_assert_cmpdouble_e (norm_test, ==, lnnorm, 1.0e-14, 0.0);
      break;
    }
    case NCM_STATS_DIST_KERNEL_ST3:
    {
      const guint d             = test->dim;
      const gdouble lg_lnnorm   = lgamma (test->nu / 2.0) - lgamma ((test->nu + d) / 2.0);
      const gdouble chol_lnnorm = 0.5 * lndet_cov;
      const gdouble nc_lnnorm   = (d / 2.0) * (ncm_c_lnpi () + log (test->nu));
      const gdouble lnnorm      = ncm_stats_dist_kernel_get_lnnorm (test->kernel, cov);

      ncm_assert_cmpdouble_e (lg_lnnorm + chol_lnnorm + nc_lnnorm, ==, lnnorm, 1.0e-14, 0.0);

      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  {
    gdouble *data         = g_new (gdouble, 2 * ntests);
    NcmVector *chi2_vec   = ncm_vector_new_full (data, ntests, 2, data, g_free);
    NcmVector *kernel_vec = ncm_vector_new (ntests);

    for (i = 0; i < ntests; i++)
    {
      const gdouble chi2_i = g_test_rand_double_range (0.0, 1.0e2);

      ncm_vector_set (chi2_vec, i, chi2_i);
    }

    ncm_stats_dist_kernel_eval_unnorm_vec (test->kernel, chi2_vec, kernel_vec);

    switch (test->kernel_type)
    {
      case NCM_STATS_DIST_KERNEL_GAUSS:

        for (i = 0; i < ntests; i++)
        {
          gdouble eval_test = exp (-0.5 * ncm_vector_get (chi2_vec, i));

          ncm_assert_cmpdouble_e (ncm_vector_get (kernel_vec, i), ==, ncm_stats_dist_kernel_eval_unnorm (test->kernel, ncm_vector_get (chi2_vec, i)), 1.0e-15, 0.0);
          ncm_assert_cmpdouble_e (eval_test, ==, ncm_stats_dist_kernel_eval_unnorm (test->kernel, ncm_vector_get (chi2_vec, i)), 1.0e-15, 0.0);
        }

        break;
      case NCM_STATS_DIST_KERNEL_ST3:

        for (i = 0; i < ntests; i++)
        {
          gdouble eval_test = pow (1 + ncm_vector_get (chi2_vec, i) / test->nu, -0.5 * (test->nu + test->dim));

          ncm_assert_cmpdouble_e (ncm_vector_get (kernel_vec, i), ==, ncm_stats_dist_kernel_eval_unnorm (test->kernel, ncm_vector_get (chi2_vec, i)), 1.0e-15, 0.0);
          ncm_assert_cmpdouble_e (eval_test, ==, ncm_stats_dist_kernel_eval_unnorm (test->kernel, ncm_vector_get (chi2_vec, i)), 1.0e-15, 0.0);
        }

        break;
      default:
        g_assert_not_reached ();
        break;
    }

    ncm_vector_free (chi2_vec);
    ncm_vector_free (kernel_vec);
  }
  ncm_model_mvnd_free (model_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
  ncm_matrix_free (cov);
}

static void
test_ncm_stats_dist_kernel_sum (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  guint i        = 0;
  const guint n  = g_test_rand_int_range (5, 100);
  gdouble lnnorm = g_test_rand_double_range (1.0, 200.0);
  gdouble lambda0, gamma0, gamma1, lambda1;
  gdouble lambda_test0   = 0.0;
  gdouble lambda_test1   = 0.0;
  gdouble lnt_i0         = 0.0;
  gdouble lnt_i1         = 0.0;
  NcmVector *weights     = ncm_vector_new (n);
  NcmVector *chi2        = ncm_vector_new (n);
  NcmVector *lnnorms_vec = ncm_vector_new (n);
  NcmVector *lnK         = ncm_vector_new (n);
  GArray *t_array0       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *t_array1       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble kappa    = -0.5 * (test->nu + test->dim);
  gdouble lnt_max0       = GSL_NEGINF;
  gdouble lnt_max1       = GSL_NEGINF;
  guint i_max0           = 0;
  guint i_max1           = 0;

  for (i = 0; i < n; i++)
  {
    const gdouble j = g_test_rand_double_range (1.0, 200.0);
    const gdouble k = g_test_rand_double_range (1.0, 200.0);
    const gdouble l = g_test_rand_double_range (1.0, 200.0);

    ncm_vector_set (weights, i, j);
    ncm_vector_set (lnnorms_vec, i, k);
    ncm_vector_set (chi2, i, l);
  }


  for (i = 0; i < n; i++)
  {
    const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
    const gdouble w_i    = ncm_vector_fast_get (weights, i);
    const gdouble lnu_i  = ncm_vector_fast_get (lnnorms_vec, i);

    switch (test->kernel_type)
    {
      case NCM_STATS_DIST_KERNEL_GAUSS:
        lnt_i0 = -0.5 * chi2_i - lnu_i + log (w_i);
        lnt_i1 = -0.5 * chi2_i + log (w_i);
        break;
      case NCM_STATS_DIST_KERNEL_ST3:
        lnt_i0 = kappa * log1p (chi2_i / test->nu) - lnu_i + log (w_i);
        lnt_i1 = kappa * log1p (chi2_i / test->nu) + log (w_i);
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    if (lnt_i0 > lnt_max0)
    {
      lnt_max0 = lnt_i0;
      i_max0   = i;
    }

    if (lnt_i1 > lnt_max1)
    {
      lnt_max1 = lnt_i1;
      i_max1   = i;
    }

    g_array_insert_val (t_array0, i, lnt_i0);
    g_array_insert_val (t_array1, i, lnt_i1);
  }

  ncm_stats_dist_kernel_eval_sum0_gamma_lambda (test->kernel, chi2, weights, lnnorms_vec, lnK, &gamma0, &lambda0);
  ncm_stats_dist_kernel_eval_sum1_gamma_lambda (test->kernel, chi2, weights, lnnorm, lnK, &gamma1, &lambda1);

  for (i = 0; i < i_max0; i++)
  {
    lambda_test0 += exp (g_array_index (t_array0, gdouble, i) - lnt_max0);
    ncm_assert_cmpdouble_e (g_array_index (t_array0, gdouble, i), <, gamma0, 1.0e-15, 0.0);
  }

  for (i = 0; i < i_max1; i++)
  {
    lambda_test1 += exp (g_array_index (t_array1, gdouble, i) - lnt_max1);
    ncm_assert_cmpdouble_e (g_array_index (t_array1, gdouble, i) - lnnorm, <, gamma1, 1.0e-15, 0.0);
  }

  for (i = i_max0 + 1; i < n; i++)
  {
    lambda_test0 += exp (g_array_index (t_array0, gdouble, i) - lnt_max0);
    ncm_assert_cmpdouble_e (g_array_index (t_array0, gdouble, i), <, gamma0, 1.0e-15, 0.0);
  }

  for (i = i_max1 + 1; i < n; i++)
  {
    lambda_test1 += exp (g_array_index (t_array1, gdouble, i) - lnt_max1);
    ncm_assert_cmpdouble_e (g_array_index (t_array1, gdouble, i) - lnnorm, <, gamma1, 1.0e-15, 0.0);
  }

  ncm_assert_cmpdouble_e (lnt_max0, ==, gamma0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (lnt_max1 - lnnorm, ==, gamma1, 1.0e-15, 0.0);

  ncm_assert_cmpdouble_e (lambda_test0, ==, lambda0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (lambda_test1, ==, lambda1, 1.0e-15, 0.0);

  ncm_assert_cmpdouble_e (lambda_test0 + exp (lnt_max0), ==, lambda0 + exp (gamma0), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (lambda_test1 + exp (lnt_max1 - lnnorm), ==, lambda1 + exp (gamma1), 1.0e-15, 0.0);

  ncm_vector_free (weights);
  ncm_vector_free (chi2);
  ncm_vector_free (lnnorms_vec);
  g_clear_pointer (&t_array0, g_array_unref);
  g_clear_pointer (&t_array1, g_array_unref);
}

static void
test_ncm_stats_dist_kernel_sample (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMatrix *cov_decomp   = ncm_matrix_new (test->dim, test->dim);
  gdouble href            = g_test_rand_double_range (1.0, 200.0);
  gdouble dif             = 0.0;
  NcmVector *mu           = ncm_vector_new (test->dim);
  NcmStatsVec *test_stats = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_VAR, FALSE);
  const guint ntests      = 300 * g_test_rand_int_range (1, 5);
  guint i, j;

  for (i = 0; i < test->dim; i++)
  {
    const gdouble val1 = g_test_rand_double_range (1.0, 200.0);

    ncm_vector_set (mu, i, val1);

    for (j = 0; j < test->dim; j++)
    {
      const gdouble val2 = g_test_rand_double_range (1.0, 200.0);

      ncm_matrix_set (cov_decomp, i, j, val2);
    }
  }

  ncm_matrix_cholesky_decomp (cov_decomp, 'U');

  while (test->nfail < 20)
  {
    for (i = 0; i < ntests; i++)
    {
      NcmVector *x_test = ncm_vector_new (test->dim);

      ncm_stats_dist_kernel_sample (test->kernel, cov_decomp, href, mu, x_test, rng);

      ncm_stats_vec_append (test_stats, x_test, FALSE);
    }

    for (i = 0; i < test->dim; i++)
    {
      gdouble dif_i;

      dif_i = fabs ((ncm_stats_vec_get_mean (test_stats, i) - ncm_vector_get (mu, i)) / ncm_vector_get (mu, i));
      dif   = dif + dif_i;
    }

    if (dif >= (test->dim * 0.2))
    {
      dif         = 0.0;
      test->nfail = test->nfail + 1;
    }
    else
    {
      test->nfail = 20;
    }
  }

  g_assert_cmpfloat (dif, <, (test->dim * 0.2));



  ncm_vector_free (mu);
  ncm_matrix_free (cov_decomp);
  ncm_rng_free (rng);
  ncm_stats_vec_free (test_stats);
}

static void
test_ncm_stats_dist_kernel_free (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_dist_kernel_free, test->kernel);
}

static void
test_ncm_stats_dist_kernel_traps (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  #if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/stats/dist/kernel/gauss/invalid/stub/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  #endif
}

static void
test_ncm_stats_dist_kernel_invalid_stub (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}


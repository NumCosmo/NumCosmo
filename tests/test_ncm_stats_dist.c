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
#include <gsl/gsl_sf_trig.h>

typedef struct _TestNcmStatsDist
{
  NcmStatsDistKernel *kernel;
  NcmStatsDist *sd;
  NcmMSet *mset;
  guint dim;
  guint nfail;
  guint np;
  guint ntests;
} TestNcmStatsDist;

static void test_ncm_stats_dist_new_kde_gauss (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_new_kde_studentt (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_new_vkde_gauss (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_new_vkde_studentt (TestNcmStatsDist *test, gconstpointer pdata);

static void test_ncm_stats_dist_dens_est (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp_sampling (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp_cv_split (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp_cv_split_nofit (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_sampling (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_serialize (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_get_kernel_info (TestNcmStatsDist *test, gconstpointer pdata);

static void test_ncm_stats_dist_free (TestNcmStatsDist *test, gconstpointer pdata);

static void test_ncm_stats_dist_traps (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_invalid_stub (TestNcmStatsDist *test, gconstpointer pdata);

typedef struct _TestNcmStatsDistFunc
{
  const gchar *name;

  void (*test_func) (TestNcmStatsDist *test, gconstpointer pdata);
} TestNcmStatsDistFunc;

#define TEST_NCM_STATS_DIST_CONSTRUCTORS_LEN 4
#define TEST_NCM_STATS_DIST_TESTS_LEN 8

static TestNcmStatsDistFunc constructors[TEST_NCM_STATS_DIST_CONSTRUCTORS_LEN] = {
  {"kde/gauss",           &test_ncm_stats_dist_new_kde_gauss, },
  {"kde/studentt",        &test_ncm_stats_dist_new_kde_studentt},
  {"vkde/gauss",          &test_ncm_stats_dist_new_vkde_gauss},
  {"vkde/studentt",       &test_ncm_stats_dist_new_vkde_studentt},
};

static TestNcmStatsDistFunc tests[TEST_NCM_STATS_DIST_TESTS_LEN] = {
  {"gauss/dens/est",                   &test_ncm_stats_dist_dens_est},
  {"gauss/dens/interp",                &test_ncm_stats_dist_dens_interp},
  {"gauss/dens/interp/sampling",       &test_ncm_stats_dist_dens_interp_sampling},
  {"gauss/dens/interp/cv_split",       &test_ncm_stats_dist_dens_interp_cv_split},
  {"gauss/dens/interp/cv_split_nofit", &test_ncm_stats_dist_dens_interp_cv_split_nofit},
  {"gauss/sampling",                   &test_ncm_stats_dist_sampling},
  {"gauss/serialize",                  &test_ncm_stats_dist_serialize},
  {"gauss/get_kernel_info",            &test_ncm_stats_dist_get_kernel_info},
};

gint
main (gint argc, gchar *argv[])
{
  GEnumClass *enum_class = g_type_class_ref (NCM_TYPE_STATS_DIST_KDE_COV_TYPE);
  gint i, j, k;

  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < TEST_NCM_STATS_DIST_CONSTRUCTORS_LEN; i++)
  {
    GEnumValue *ev;

    k = 0;

    while ((ev = g_enum_get_value (enum_class, k++)) != NULL)
    {
      for (j = 0; j < TEST_NCM_STATS_DIST_TESTS_LEN; j++)
      {
        gchar *test_name = g_strdup_printf ("/ncm/stats/dist/%s/%s/%s", constructors[i].name, tests[j].name, ev->value_nick);

        g_test_add (test_name,
                    TestNcmStatsDist,
                    GINT_TO_POINTER (ev->value),
                    constructors[i].test_func,
                    tests[j].test_func,
                    &test_ncm_stats_dist_free);

        g_free (test_name);
      }
    }
  }

  g_test_add ("/ncm/stats/dist/nd/kde/gauss/traps", TestNcmStatsDist, NULL,
              &test_ncm_stats_dist_new_kde_gauss,
              &test_ncm_stats_dist_traps,
              &test_ncm_stats_dist_free);
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/stats/dist/nd/kde/gauss/invalid/stub/subprocess", TestNcmStatsDist, NULL,
              &test_ncm_stats_dist_new_kde_gauss,
              &test_ncm_stats_dist_invalid_stub,
              &test_ncm_stats_dist_free);
#endif

  g_test_run ();
}

#define TESTMULT 200
#define NTESTS 2000

static void
test_ncm_stats_dist_new_kde_gauss (TestNcmStatsDist *test, gconstpointer pdata)
{
  const guint dim                    = g_test_rand_int_range (2, 4);
  NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);
  NcmStatsDistKDE *sdkde             = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (sdk_gauss), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type    = GPOINTER_TO_INT (pdata);

  test->dim    = dim;
  test->kernel = NCM_STATS_DIST_KERNEL (sdk_gauss);
  test->sd     = NCM_STATS_DIST (sdkde);
  test->nfail  = 0;
  test->np     = TESTMULT * test->dim;
  test->ntests = NTESTS;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->np *= 4;
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_stats_dist_kernel_gauss_ref (sdk_gauss);
  ncm_stats_dist_kernel_gauss_free (sdk_gauss);
  {
    NcmStatsDistKernelGauss *sdk_gauss0 = ncm_stats_dist_kernel_gauss_ref (sdk_gauss);

    ncm_stats_dist_kernel_gauss_clear (&sdk_gauss0);
    g_assert_true (sdk_gauss0 == NULL);
  }

  ncm_stats_dist_kde_ref (sdkde);
  ncm_stats_dist_kde_free (sdkde);
  {
    NcmStatsDistKDE *sdkde0 = ncm_stats_dist_kde_ref (sdkde);

    ncm_stats_dist_kde_clear (&sdkde0);
    g_assert_true (sdkde0 == 0);
  }

  ncm_stats_dist_kde_set_cov_type (sdkde, cov_type);
}

static void
test_ncm_stats_dist_new_kde_studentt (TestNcmStatsDist *test, gconstpointer pdata)
{
  const gdouble nu                = g_test_rand_double_range (3.0, 5.0);
  const guint dim                 = g_test_rand_int_range (2, 4);
  NcmStatsDistKernelST *sdk_st    = ncm_stats_dist_kernel_st_new (dim, nu);
  NcmStatsDistKDE *sdkde          = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (sdk_st), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);

  test->dim    = dim;
  test->kernel = NCM_STATS_DIST_KERNEL (sdk_st);
  test->sd     = NCM_STATS_DIST (sdkde);
  test->nfail  = 0;
  test->np     = TESTMULT * test->dim;
  test->ntests = NTESTS;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->np *= 4;
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_stats_dist_kernel_st_ref (sdk_st);
  ncm_stats_dist_kernel_st_free (sdk_st);
  {
    NcmStatsDistKernelST *sdk_st0 = ncm_stats_dist_kernel_st_ref (sdk_st);

    ncm_stats_dist_kernel_st_clear (&sdk_st0);
    g_assert_true (sdk_st0 == NULL);
  }

  ncm_stats_dist_kde_ref (sdkde);
  ncm_stats_dist_kde_free (sdkde);
  {
    NcmStatsDistKDE *sdkde0 = ncm_stats_dist_kde_ref (sdkde);

    ncm_stats_dist_kde_clear (&sdkde0);
    g_assert_true (sdkde0 == 0);
  }

  ncm_stats_dist_kde_set_cov_type (sdkde, cov_type);
}

static void
test_ncm_stats_dist_new_vkde_gauss (TestNcmStatsDist *test, gconstpointer pdata)
{
  const guint dim                    = g_test_rand_int_range (2, 4);
  NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);
  NcmStatsDistVKDE *sdvkde           = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (sdk_gauss), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type    = GPOINTER_TO_INT (pdata);

  test->dim    = dim;
  test->kernel = NCM_STATS_DIST_KERNEL (sdk_gauss);
  test->sd     = NCM_STATS_DIST (sdvkde);
  test->nfail  = 0;
  test->np     = TESTMULT * test->dim;
  test->ntests = NTESTS;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->np *= 4;
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_stats_dist_kernel_gauss_ref (sdk_gauss);
  ncm_stats_dist_kernel_gauss_free (sdk_gauss);
  {
    NcmStatsDistKernelGauss *sdk_gauss0 = ncm_stats_dist_kernel_gauss_ref (sdk_gauss);

    ncm_stats_dist_kernel_gauss_clear (&sdk_gauss0);
    g_assert_true (sdk_gauss0 == NULL);
  }

  ncm_stats_dist_vkde_ref (sdvkde);
  ncm_stats_dist_vkde_free (sdvkde);
  {
    NcmStatsDistVKDE *sdvkde0 = ncm_stats_dist_vkde_ref (sdvkde);

    ncm_stats_dist_vkde_clear (&sdvkde0);
    g_assert_true (sdvkde0 == 0);
  }

  ncm_stats_dist_set_over_smooth (test->sd, 0.2);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (sdvkde), cov_type);
}

static void
test_ncm_stats_dist_new_vkde_studentt (TestNcmStatsDist *test, gconstpointer pdata)
{
  const gdouble nu                = g_test_rand_double_range (3.0, 5.0);
  const guint dim                 = g_test_rand_int_range (2, 4);
  NcmStatsDistKernelST *sdk_st    = ncm_stats_dist_kernel_st_new (dim, nu);
  NcmStatsDistVKDE *sdvkde        = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (sdk_st), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);

  test->dim    = dim;
  test->kernel = NCM_STATS_DIST_KERNEL (sdk_st);
  test->sd     = NCM_STATS_DIST (sdvkde);
  test->nfail  = 0;
  test->np     = TESTMULT * test->dim;
  test->ntests = NTESTS;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->np *= 4;
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_stats_dist_kernel_st_ref (sdk_st);
  ncm_stats_dist_kernel_st_free (sdk_st);
  {
    NcmStatsDistKernelST *sdk_st0 = ncm_stats_dist_kernel_st_ref (sdk_st);

    ncm_stats_dist_kernel_st_clear (&sdk_st0);
    g_assert_true (sdk_st0 == NULL);
  }

  ncm_stats_dist_vkde_ref (sdvkde);
  ncm_stats_dist_vkde_free (sdvkde);
  {
    NcmStatsDistVKDE *sdvkde0 = ncm_stats_dist_vkde_ref (sdvkde);

    ncm_stats_dist_vkde_clear (&sdvkde0);
    g_assert_true (sdvkde0 == 0);
  }

  ncm_stats_dist_set_over_smooth (test->sd, 0.2);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (sdvkde), cov_type);
}

static void
test_ncm_stats_dist_free (TestNcmStatsDist *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_dist_free, test->sd);
  NCM_TEST_FREE (ncm_stats_dist_kernel_free, test->kernel);
}

static void
test_ncm_stats_dist_cmp_dist (TestNcmStatsDist *test, NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, NcmRNG *rng)
{
  NcmStatsVec *err_stats = ncm_stats_vec_new (2, NCM_STATS_VEC_VAR, FALSE);
  gulong N               = 0;
  gint i;

  for (i = 0; i < test->ntests; i++)
  {
    NcmVector *y;
    gdouble m2lnL, m2lnp_s, alpha;

    y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    m2lnp_s = ncm_stats_dist_eval_m2lnp (test->sd, y);

    /* Measuring the Kullback-Leibler divergence */
    alpha = 0.25 * (m2lnp_s - m2lnL);
    ncm_stats_vec_set (err_stats, 0, log1p (tanh (alpha)));

    ncm_stats_dist_sample (test->sd, y, rng);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    m2lnp_s = ncm_stats_dist_eval_m2lnp (test->sd, y);

    alpha = 0.25 * (m2lnL - m2lnp_s);
    ncm_stats_vec_set (err_stats, 1, log1p (tanh (alpha)));

    ncm_stats_vec_update (err_stats);
  }

  /* Measuring the Jensen-Shannon divergence */
  {
    gdouble KL_P_M  = ncm_stats_vec_get_mean (err_stats, 0);
    gdouble KL_PS_M = ncm_stats_vec_get_mean (err_stats, 1);
    gdouble JS_P_PS = 0.5 * (KL_P_M + KL_PS_M);

    g_assert_cmpfloat (JS_P_PS, <, 0.5);
    /*printf ("KL % 22.15g % 22.15g JS % 22.15g\n", KL_P_M, KL_PS_M, JS_P_PS); */
  }

  ncm_stats_vec_free (err_stats);
}

static void
test_ncm_stats_dist_dens_est (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  gulong N                        = 0;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);

    ncm_stats_dist_add_obs (test->sd, y);
  }

  ncm_stats_dist_prepare (test->sd);

  test_ncm_stats_dist_cmp_dist (test, data_mvnd, mset, rng);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_dens_interp (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  gulong N                        = 0;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));
  ncm_data_gauss_cov_use_norma (NCM_DATA_GAUSS_COV (data_mvnd), TRUE);

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_prepare_interp (test->sd, m2lnp_v);

  test_ncm_stats_dist_cmp_dist (test, data_mvnd, mset, rng);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_dens_interp_sampling (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  const guint ntests = 10000000;
  gulong N = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v = ncm_vector_new (test->np);
  NcmVector *weights, *cum;
  gdouble sum_weights;
  guint i, n;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_prepare_interp (test->sd, m2lnp_v);

  n   = ncm_stats_dist_get_sample_size (test->sd);
  cum = ncm_vector_new (n);

  ncm_vector_set_zero (cum);

  for (i = 0; i < ntests; i++)
  {
    const gint k_i = ncm_stats_dist_kernel_choose (test->sd, rng);

    g_assert_true (k_i < n);
    ncm_vector_addto (cum, k_i, 1.0);
  }

  ncm_vector_scale (cum, 1.0 / ntests);

  weights     = ncm_stats_dist_peek_weights (test->sd);
  sum_weights = ncm_vector_sum_cpts (weights);

  for (i = 0; i < test->np; i++)
  {
    const gdouble c_i = ncm_vector_get (cum, i);
    const gdouble w_i = ncm_vector_get (weights, i) / sum_weights;

    if (w_i == 0.0)
      g_assert_true (c_i == 0.0);
    else
      ncm_assert_cmpdouble_e (c_i, ==, w_i, 1.0e-1, 1.0e-4);
  }

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (cum);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_dens_interp_cv_split (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_SPLIT);
  ncm_stats_dist_prepare_interp (test->sd, m2lnp_v);

  test_ncm_stats_dist_cmp_dist (test, data_mvnd, mset, rng);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_dens_interp_cv_split_nofit (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  gdouble dm2lnL_mean             = 0.0;
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  guint ntests_fail;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_SPLIT_NOFIT);
  ncm_stats_dist_prepare_interp (test->sd, m2lnp_v);

/*
 *  sample_array = ncm_stats_dist_peek_sample_array (test->sd);
 *  for (i = 0; i < sample_array->len; i++)
 *  {
 *   NcmVector *y = g_ptr_array_index (sample_array, i);
 *   const gdouble m2lnp = ncm_stats_dist_eval_m2lnp (test->sd, y);
 *
 *   printf ("#A % 22.15g % 22.15g % 22.15g\n", m2lnp, ncm_vector_get (m2lnp_v, i), ncm_vector_get (m2lnp_v, i) - m2lnp);
 *  }
 */

  for (i = 0; i < test->ntests; i++)
  {
    NcmVector *y    = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnp_s = ncm_stats_dist_eval_m2lnp (test->sd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    dm2lnL_mean += (m2lnp_s - m2lnL);
  }

  dm2lnL_mean = dm2lnL_mean / test->ntests;

  ntests_fail = 0;

  for (i = 0; i < test->ntests; i++)
  {
    NcmVector *y    = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnp_s = ncm_stats_dist_eval_m2lnp (test->sd, y);
    gdouble m2lnL;

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);

    /*printf ("test % 22.15g\n", fabs (m2lnp_s / (m2lnL + dm2lnL_mean) - 1.0));*/
    if (fabs (m2lnp_s / (m2lnL + dm2lnL_mean) - 1.0) > 0.2)
      ntests_fail++;
  }

  /*printf ("%u %u %u %f\n", test->dim, ntests, ntests_fail, ntests_fail * 1.0 / ntests);*/

  g_assert_cmpfloat (ntests_fail * 1.0 / test->ntests, <, 0.5);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_sampling (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmVector *y                    = ncm_vector_new (test->dim);
  NcmStatsVec *test_stats         = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_COV, FALSE);
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    /*ncm_vector_log_vals (y, "Y: ", "% 12.5g", TRUE);*/

    ncm_stats_dist_add_obs (test->sd, y);

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
  }

  ncm_stats_dist_prepare (test->sd);

  for (i = 0; i < test->ntests; i++)
  {
    ncm_stats_dist_sample (test->sd, y, rng);
    ncm_stats_vec_append (test_stats, y, FALSE);
  }

  {
    NcmMatrix *cov_est = ncm_stats_vec_peek_cov_matrix (test_stats, 0);
    gdouble chi2;
    gint ret;

    if (
      (test->nfail < 10) &&
      ((ncm_matrix_cmp (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0) >= 0.5) ||
       (ncm_matrix_cmp_diag (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0) >= 0.5))
       )
    {
      guint nfail = test->nfail;

      test_ncm_stats_dist_free (test, pdata);
      test_ncm_stats_dist_new_kde_gauss (test, pdata);

      test->nfail = nfail + 1;
      test_ncm_stats_dist_sampling (test, pdata);
    }
    else
    {
      g_assert_cmpfloat (ncm_matrix_cmp (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0), <, 0.5);
      g_assert_cmpfloat (ncm_matrix_cmp_diag (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0), <, 0.5);
    }

    if (FALSE)
    {
      ncm_cfg_msg_sepa ();

      printf ("# WDIFF   : % 22.15e\n", ncm_matrix_cmp (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0));
      printf ("# WDIFFD  : % 22.15e\n", ncm_matrix_cmp_diag (cov_est, NCM_DATA_GAUSS_COV (data_mvnd)->cov, 0.0));


      ncm_vector_memcpy (y, NCM_DATA_GAUSS_COV (data_mvnd)->y);
      ncm_vector_sub (y, ncm_stats_vec_peek_mean (test_stats));

      ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                            ncm_matrix_gsl (NCM_DATA_GAUSS_COV (data_mvnd)->LLT), ncm_vector_gsl (y));
      NCM_TEST_GSL_RESULT ("test_ncm_stats_dist_gauss_sampling", ret);

      ret = gsl_blas_ddot (ncm_vector_gsl (y), ncm_vector_gsl (y), &chi2);
      NCM_TEST_GSL_RESULT ("test_ncm_stats_dist_gauss_sampling", ret);

      printf ("# chi2 == % 22.15g ~ %u +/- % 22.15g\n", chi2, test->dim, sqrt (2.0 * test->dim));

      ncm_vector_log_vals (NCM_DATA_GAUSS_COV (data_mvnd)->y,    "MEAN:     ", "% 12.5g", TRUE);
      ncm_vector_log_vals (ncm_stats_vec_peek_mean (test_stats), "MEAN_EST: ", "% 12.5g", TRUE);
      ncm_vector_log_vals (y,                                    "MEAN_CMP: ", "% 12.5g", TRUE);

      ncm_matrix_log_vals (NCM_DATA_GAUSS_COV (data_mvnd)->cov,  "COV:     ", "% 12.5g");
      ncm_matrix_log_vals (cov_est,                              "COV_EST: ", "% 12.5g");
    }
  }

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (y);
  ncm_stats_vec_free (test_stats);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_serialize (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  NcmStatsVec *cmp_stats          = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  gulong N                        = 0;
  NcmSerialize *ser;
  gchar *sd_ser;
  NcmStatsDist *sd_dup;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  sd_ser = ncm_serialize_to_string (ser, G_OBJECT (test->sd), TRUE);
  sd_dup = NCM_STATS_DIST (ncm_serialize_from_string (ser, sd_ser));

  ncm_data_gauss_cov_use_norma (NCM_DATA_GAUSS_COV (data_mvnd), FALSE);
  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);
    ncm_stats_dist_add_obs (sd_dup, y);

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_prepare_interp (test->sd, m2lnp_v);
  ncm_stats_dist_prepare_interp (sd_dup, m2lnp_v);

  for (i = 0; i < test->ntests; i++)
  {
    NcmVector *y     = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnp_s0 = ncm_stats_dist_eval_m2lnp (test->sd, y);
    gdouble m2lnp_s1 = ncm_stats_dist_eval_m2lnp (sd_dup, y);

    ncm_assert_cmpdouble_e (m2lnp_s0, ==, m2lnp_s1, 1.0e-14, 0.0);
  }

  g_free (sd_ser);
  ncm_stats_dist_free (sd_dup);
  ncm_serialize_free (ser);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
  ncm_stats_vec_clear (&cmp_stats);
}

static void
test_ncm_stats_dist_get_kernel_info (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-1, 1.0, -2.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), gcov->cov);
      break;
    }
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      break;
    default:
      g_assert_not_reached ();
  }

  ncm_data_gauss_cov_use_norma (NCM_DATA_GAUSS_COV (data_mvnd), FALSE);
  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);

    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  ncm_stats_dist_prepare (test->sd);
  ncm_stats_dist_prepare_kernel (test->sd, ncm_stats_dist_peek_sample_array (test->sd));
  ncm_stats_dist_prepare_interp (test->sd, m2lnp_v);

  {
    const gdouble rnorm        = ncm_stats_dist_get_rnorm (test->sd);
    const gdouble href         = ncm_stats_dist_get_href (test->sd);
    const guint dim            = ncm_stats_dist_get_dim (test->sd);
    const guint n              = ncm_stats_dist_get_sample_size (test->sd);
    GPtrArray *sample_array    = ncm_stats_dist_peek_sample_array (test->sd);
    NcmVector *weights         = ncm_stats_dist_peek_weights (test->sd);
    NcmStatsDistKernel *kernel = ncm_stats_dist_get_kernel (test->sd);


    g_assert_true (gsl_finite (rnorm));
    g_assert_true (sample_array->len == n);

    for (i = 0; i < n; i++)
    {
      NcmMatrix *cov_decomp = ncm_stats_dist_peek_cov_decomp (test->sd, i);
      const gdouble lnnorm_i = ncm_stats_dist_get_lnnorm (test->sd, i);
      const gdouble lnnorm0_i = ncm_stats_dist_kernel_get_lnnorm (kernel, cov_decomp);
      NcmVector *y_i = NULL;
      NcmMatrix *cov_i = NULL;
      gdouble w_i, n_i;

      ncm_stats_dist_get_Ki (test->sd, i, &y_i, &cov_i, &n_i, &w_i);

      ncm_assert_cmpdouble_e (n_i, ==, exp (lnnorm_i), 1.0e-14, 0.0);
      ncm_assert_cmpdouble_e (lnnorm_i, ==, lnnorm0_i + dim * log (href), 1.0e-14, 0.0);
      g_assert_true (w_i == ncm_vector_get (weights, i));

      {
        NcmVector *y_i_dup = ncm_vector_dup (y_i);

        ncm_vector_cmp (y_i_dup, g_ptr_array_index (sample_array, i));
        g_assert_true (ncm_vector_get_max (y_i_dup) < 1.0e-15);
        ncm_vector_free (y_i_dup);
      }

      ncm_vector_free (y_i);
      ncm_matrix_free (cov_i);
    }

    {
      gdouble *data         = g_new (gdouble, 2 * test->ntests);
      NcmVector *chi2_vec   = ncm_vector_new_full (data, test->ntests, 2, data, g_free);
      NcmVector *kernel_vec = ncm_vector_new (test->ntests);

      for (i = 0; i < test->ntests; i++)
      {
        const gdouble chi2_i = g_test_rand_double_range (0.0, 1.0e2);

        ncm_vector_set (chi2_vec, i, chi2_i);
      }

      ncm_stats_dist_kernel_eval_unnorm_vec (kernel, chi2_vec, kernel_vec);

      for (i = 0; i < test->ntests; i++)
      {
        ncm_assert_cmpdouble_e (ncm_vector_get (kernel_vec, i), ==, ncm_stats_dist_kernel_eval_unnorm (kernel, ncm_vector_get (chi2_vec, i)),
                                1.0e-15, 0.0);
      }

      ncm_vector_free (chi2_vec);
      ncm_vector_free (kernel_vec);
    }

    ncm_stats_dist_kernel_free (kernel);
  }

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_traps (TestNcmStatsDist *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/stats/dist/nd/kde/gauss/invalid/stub/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

static void
test_ncm_stats_dist_invalid_stub (TestNcmStatsDist *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}


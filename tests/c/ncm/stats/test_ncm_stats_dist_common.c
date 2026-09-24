/***************************************************************************
 *            test_ncm_stats_dist_common.c
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include "test_ncm_stats_dist_common.h"

typedef struct _TestNcmStatsDist
{
  NcmStatsDistKernel *kernel;
  NcmStatsDist *sd;
  NcmMSet *mset;
  guint dim;
  guint nfail;
  guint np;
  guint ntests;
  gdouble corr_level;
  gboolean divergence;
} TestNcmStatsDist;

/* Chosen once by main(): g_test_add()'s tdata already carries the covariance type, and
 * every test in one executable runs in the same mode. */
static TestNcmStatsDistMode _test_mode = TEST_NCM_STATS_DIST_MECHANICS;

static void test_ncm_stats_dist_new_kde_gauss (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_new_kde_studentt (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_new_vkde_gauss (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_new_vkde_studentt (TestNcmStatsDist *test, gconstpointer pdata);

static void test_ncm_stats_dist_sanity (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_prepare_too_few (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_est (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp_sampling (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp_cv_split_m2lnp (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_dens_interp_cv_loo (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_sampling (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_serialize (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_get_kernel_info (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_center_shrink (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_defensive (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_eval_vec (TestNcmStatsDist *test, gconstpointer pdata);
void test_ncm_stats_dist_cv_auto_kernel (TestNcmStatsDist *test, gconstpointer pdata);

static void test_ncm_stats_dist_free (TestNcmStatsDist *test, gconstpointer pdata);

static void test_ncm_stats_dist_traps (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_invalid_stub (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_invalid_center_shrink (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_vkde_points_per_dim (void);
static void test_ncm_stats_dist_cv_objectives (void);
static void test_ncm_stats_dist_loo_batch (void);
static void test_ncm_stats_dist_invalid_cv_accept (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_invalid_auto_kernel_gauss (TestNcmStatsDist *test, gconstpointer pdata);
static void test_ncm_stats_dist_split_underflowing_m2lnp (void);
static void test_ncm_stats_dist_split_drop_far_points (void);
static void test_ncm_stats_dist_print_fit (void);
static void test_ncm_stats_dist_kde_cov_fixed_nearPD (void);

typedef struct _TestNcmStatsDistFunc
{
  const gchar *name;

  void (*test_func) (TestNcmStatsDist *test, gconstpointer pdata);
} TestNcmStatsDistFunc;

#define TEST_NCM_STATS_DIST_CONSTRUCTORS_LEN 4
#define TEST_NCM_STATS_DIST_TESTS_LEN 14

static TestNcmStatsDistFunc constructors[TEST_NCM_STATS_DIST_CONSTRUCTORS_LEN] = {
  {"kde/gauss",           &test_ncm_stats_dist_new_kde_gauss, },
  {"kde/studentt",        &test_ncm_stats_dist_new_kde_studentt},
  {"vkde/gauss",          &test_ncm_stats_dist_new_vkde_gauss},
  {"vkde/studentt",       &test_ncm_stats_dist_new_vkde_studentt},
};

static TestNcmStatsDistFunc tests[TEST_NCM_STATS_DIST_TESTS_LEN] = {
  {"sanity",                           &test_ncm_stats_dist_sanity},
  {"prepare_too_few",                  &test_ncm_stats_dist_prepare_too_few},
  {"gauss/dens/est",                   &test_ncm_stats_dist_dens_est},
  {"gauss/dens/interp",                &test_ncm_stats_dist_dens_interp},
  {"gauss/dens/interp/sampling",       &test_ncm_stats_dist_dens_interp_sampling},
  {"gauss/dens/interp/cv_split_m2lnp", &test_ncm_stats_dist_dens_interp_cv_split_m2lnp},
  {"gauss/dens/interp/cv_loo",         &test_ncm_stats_dist_dens_interp_cv_loo},
  {"gauss/sampling",                   &test_ncm_stats_dist_sampling},
  {"gauss/serialize",                  &test_ncm_stats_dist_serialize},
  {"gauss/get_kernel_info",            &test_ncm_stats_dist_get_kernel_info},
  {"gauss/center_shrink",              &test_ncm_stats_dist_center_shrink},
  {"gauss/defensive",                  &test_ncm_stats_dist_defensive},
  {"eval_vec",                         &test_ncm_stats_dist_eval_vec},
  {"gauss/cv_auto_kernel",             &test_ncm_stats_dist_cv_auto_kernel},
};

/* The checks that assert the estimator reproduces its own distribution. In divergence
 * mode only these run; the others carry no statistical claim and would only repeat what
 * the mechanics run already covered. */
static gboolean
_test_ncm_stats_dist_is_divergence_check (const gchar *name)
{
  return g_str_has_prefix (name, "gauss/dens/") ||
         (g_strcmp0 (name, "gauss/sampling") == 0) ||
         (g_strcmp0 (name, "gauss/center_shrink") == 0);
}

gint
test_ncm_stats_dist_main (gint argc, gchar *argv[], TestNcmStatsDistMode mode)
{
  GEnumClass *enum_class = g_type_class_ref (NCM_TYPE_STATS_DIST_KDE_COV_TYPE);
  gint i, j, k;

  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  _test_mode = mode;

  for (i = 0; i < TEST_NCM_STATS_DIST_CONSTRUCTORS_LEN; i++)
  {
    GEnumValue *ev;

    k = 0;

    while ((ev = g_enum_get_value (enum_class, k++)) != NULL)
    {
      for (j = 0; j < TEST_NCM_STATS_DIST_TESTS_LEN; j++)
      {
        gchar *test_name;

        if ((mode == TEST_NCM_STATS_DIST_DIVERGENCE) && !_test_ncm_stats_dist_is_divergence_check (tests[j].name))
          continue;

        test_name = g_strdup_printf ("/ncm/stats/dist/%s/%s/%s", constructors[i].name, tests[j].name, ev->value_nick);

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

  if (mode == TEST_NCM_STATS_DIST_MECHANICS)
  {
    g_test_add ("/ncm/stats/dist/nd/kde/gauss/traps", TestNcmStatsDist, NULL,
                &test_ncm_stats_dist_new_kde_gauss,
                &test_ncm_stats_dist_traps,
                &test_ncm_stats_dist_free);

    g_test_add ("/ncm/stats/dist/nd/kde/gauss/invalid/stub/subprocess", TestNcmStatsDist, NULL,
                &test_ncm_stats_dist_new_kde_gauss,
                &test_ncm_stats_dist_invalid_stub,
                &test_ncm_stats_dist_free);

    g_test_add ("/ncm/stats/dist/nd/vkde/cauchy/invalid/center_shrink/subprocess", TestNcmStatsDist, NULL,
                &test_ncm_stats_dist_new_kde_gauss,
                &test_ncm_stats_dist_invalid_center_shrink,
                &test_ncm_stats_dist_free);

    g_test_add_func ("/ncm/stats/dist/nd/vkde/gauss/split/underflowing_m2lnp",
                     &test_ncm_stats_dist_split_underflowing_m2lnp);
    g_test_add_func ("/ncm/stats/dist/nd/kde/gauss/split/drop_far_points",
                     &test_ncm_stats_dist_split_drop_far_points);
    g_test_add_func ("/ncm/stats/dist/nd/print_fit",
                     &test_ncm_stats_dist_print_fit);
    g_test_add_func ("/ncm/stats/dist/nd/kde/gauss/cov_fixed/nearPD",
                     &test_ncm_stats_dist_kde_cov_fixed_nearPD);
    g_test_add_func ("/ncm/stats/dist/nd/vkde/gauss/points_per_dim",
                     &test_ncm_stats_dist_vkde_points_per_dim);
    g_test_add_func ("/ncm/stats/dist/nd/kde/gauss/cv_objectives",
                     &test_ncm_stats_dist_cv_objectives);
    g_test_add_func ("/ncm/stats/dist/nd/vkde/st/loo_batch",
                     &test_ncm_stats_dist_loo_batch);
    g_test_add ("/ncm/stats/dist/nd/kde/gauss/invalid/cv_accept_without_m2lnL/subprocess", TestNcmStatsDist, NULL,
                &test_ncm_stats_dist_new_kde_gauss,
                &test_ncm_stats_dist_invalid_cv_accept,
                &test_ncm_stats_dist_free);
    g_test_add ("/ncm/stats/dist/nd/kde/gauss/invalid/auto_kernel_gauss/subprocess", TestNcmStatsDist, NULL,
                &test_ncm_stats_dist_new_kde_gauss,
                &test_ncm_stats_dist_invalid_auto_kernel_gauss,
                &test_ncm_stats_dist_free);
  }

  return g_test_run ();
}

#define TESTMULT 200
#define NTESTS 500
#define RELTOL 0.5

/* Mechanics mode sizing. The interpolation fits dominate and their cost grows with the
 * sample, so this is what makes the difference; the guard the library applies is
 * n_obs > dim, which these clear by a wide margin. */
#define TESTMULT_MECHANICS 80
#define NTESTS_MECHANICS 20
#define SAMPLING_NTESTS_MECHANICS 10000

static void
test_ncm_stats_dist_new_kde_gauss (TestNcmStatsDist *test, gconstpointer pdata)
{
  const guint dim                    = g_test_rand_int_range (1, 4);
  NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);
  NcmStatsDistKDE *sdkde             = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (sdk_gauss), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type    = GPOINTER_TO_INT (pdata);

  test->dim        = dim;
  test->kernel     = NCM_STATS_DIST_KERNEL (sdk_gauss);
  test->sd         = NCM_STATS_DIST (sdkde);
  test->nfail      = 0;
  test->divergence = (_test_mode == TEST_NCM_STATS_DIST_DIVERGENCE);
  test->np         = (test->divergence ? TESTMULT : TESTMULT_MECHANICS) * test->dim;
  test->ntests     = test->divergence ? NTESTS : NTESTS_MECHANICS;
  test->corr_level = 1.0;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->corr_level = 100.0;
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
  const guint dim                 = g_test_rand_int_range (1, 4);
  NcmStatsDistKernelST *sdk_st    = ncm_stats_dist_kernel_st_new (dim, nu);
  NcmStatsDistKDE *sdkde          = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (sdk_st), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);

  test->dim        = dim;
  test->kernel     = NCM_STATS_DIST_KERNEL (sdk_st);
  test->sd         = NCM_STATS_DIST (sdkde);
  test->nfail      = 0;
  test->divergence = (_test_mode == TEST_NCM_STATS_DIST_DIVERGENCE);
  test->np         = (test->divergence ? TESTMULT : TESTMULT_MECHANICS) * test->dim;
  test->ntests     = test->divergence ? NTESTS : NTESTS_MECHANICS;
  test->corr_level = 1.0;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->corr_level = 100.0;
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
  const guint dim                    = g_test_rand_int_range (1, 4);
  NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);
  NcmStatsDistVKDE *sdvkde           = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (sdk_gauss), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type    = GPOINTER_TO_INT (pdata);

  test->dim        = dim;
  test->kernel     = NCM_STATS_DIST_KERNEL (sdk_gauss);
  test->sd         = NCM_STATS_DIST (sdvkde);
  test->nfail      = 0;
  test->divergence = (_test_mode == TEST_NCM_STATS_DIST_DIVERGENCE);
  test->np         = (test->divergence ? TESTMULT : TESTMULT_MECHANICS) * test->dim;
  test->ntests     = test->divergence ? NTESTS : NTESTS_MECHANICS;
  test->corr_level = 1.0;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->corr_level = 100.0;
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

  ncm_stats_dist_set_over_smooth (test->sd, 1.2);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (sdvkde), cov_type);
}

static void
test_ncm_stats_dist_new_vkde_studentt (TestNcmStatsDist *test, gconstpointer pdata)
{
  const gdouble nu                = g_test_rand_double_range (3.0, 5.0);
  const guint dim                 = g_test_rand_int_range (1, 4);
  NcmStatsDistKernelST *sdk_st    = ncm_stats_dist_kernel_st_new (dim, nu);
  NcmStatsDistVKDE *sdvkde        = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (sdk_st), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);

  test->dim        = dim;
  test->kernel     = NCM_STATS_DIST_KERNEL (sdk_st);
  test->sd         = NCM_STATS_DIST (sdvkde);
  test->nfail      = 0;
  test->divergence = (_test_mode == TEST_NCM_STATS_DIST_DIVERGENCE);
  test->np         = (test->divergence ? TESTMULT : TESTMULT_MECHANICS) * test->dim;
  test->ntests     = test->divergence ? NTESTS : NTESTS_MECHANICS;
  test->corr_level = 1.0;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
      break;
    case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
      test->corr_level = 100.0;
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

  ncm_stats_dist_set_over_smooth (test->sd, 1.2);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (sdvkde), cov_type);
}

/* Centre shrinkage: with it off the kernel centres are the sample points; with it on
 * they are mu + a (x_i - mu) with a = 1 / sqrt (1 + h^2 s^2), s = 1 for the fixed
 * bandwidth estimator, and the mixture covariance matches the sample covariance. */
static void
test_ncm_stats_dist_center_shrink (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmStatsVec *sample_stats       = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_COV, FALSE);
  NcmStatsVec *test_stats         = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_COV, FALSE);
  NcmVector *y                    = ncm_vector_new (test->dim);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  gulong N                        = 0;
  guint i, k;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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
    NcmVector *y_i = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y_i);
    ncm_stats_vec_append (sample_stats, y_i, FALSE);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  /* Off: the centres are the sample points and the factor is one. */
  ncm_stats_dist_set_center_shrink (test->sd, FALSE);
  g_assert_false (ncm_stats_dist_get_center_shrink (test->sd));

  ncm_stats_dist_prepare (test->sd, NULL);

  g_assert_cmpfloat (ncm_stats_dist_get_center_shrink_factor (test->sd), ==, 1.0);

  {
    GPtrArray *sample_array = ncm_stats_dist_peek_sample_array (test->sd);
    GPtrArray *center_array = ncm_stats_dist_peek_center_array (test->sd);

    g_assert_cmpuint (center_array->len, ==, ncm_stats_dist_get_n_kernels (test->sd));

    for (i = 0; i < center_array->len; i++)
    {
      NcmVector *x_i = g_ptr_array_index (sample_array, i);
      NcmVector *c_i = g_ptr_array_index (center_array, i);

      for (k = 0; k < test->dim; k++)
        g_assert_cmpfloat (ncm_vector_get (c_i, k), ==, ncm_vector_get (x_i, k));
    }
  }

  /* On: c_i = mu + a (x_i - mu). */
  ncm_stats_dist_set_center_shrink (test->sd, TRUE);
  g_assert_true (ncm_stats_dist_get_center_shrink (test->sd));

  ncm_stats_dist_prepare (test->sd, NULL);

  {
    GPtrArray *sample_array = ncm_stats_dist_peek_sample_array (test->sd);
    GPtrArray *center_array = ncm_stats_dist_peek_center_array (test->sd);
    NcmVector *mean         = ncm_stats_vec_peek_mean (sample_stats);
    const gdouble a         = ncm_stats_dist_get_center_shrink_factor (test->sd);
    const gdouble href      = ncm_stats_dist_get_href (test->sd);
    const gdouble kappa     = ncm_stats_dist_kernel_get_var_factor (test->kernel);

    g_assert_cmpfloat (a, >, 0.0);
    g_assert_cmpfloat (a, <, 1.0);
    g_assert_true (gsl_finite (kappa));

    /*
     * Center shrinkage contracts the centers by the matrix A and the kernel scale
     * matrices by Ahat = A / a, with det Ahat = 1, so that the mixture covariance equals
     * the sample covariance C exactly:
     *
     *   A C A^T + kappa href^2 <Sigma'> = C,
     *
     * where href = a h is the applied bandwidth and <Sigma'> the mean of the applied
     * kernel scale matrices U_i'^T U_i'. Everything on the left is public. For the fixed
     * bandwidth estimator with the sample covariance the transform is isotropic, A = a I.
     */
    {
      const guint d         = test->dim;
      const guint n_kernels = ncm_stats_dist_get_n_kernels (test->sd);
      NcmMatrix *A          = ncm_stats_dist_peek_center_shrink_matrix (test->sd);
      NcmMatrix *C          = ncm_stats_vec_peek_cov_matrix (sample_stats, 0);
      NcmMatrix *mean_cov   = ncm_matrix_new (d, d);
      NcmMatrix *B          = ncm_matrix_new (d, d);
      NcmMatrix *AC         = ncm_matrix_new (d, d);
      NcmMatrix *lhs        = ncm_matrix_new (d, d);
      gdouble max_C         = 0.0;
      gdouble max_dev       = 0.0;
      guint p, q, j;

      g_assert_nonnull (A);
      g_assert_cmpuint (ncm_matrix_nrows (A), ==, d);
      g_assert_cmpuint (ncm_matrix_ncols (A), ==, d);

      /* <Sigma'> from the applied factors; Cholesky leaves the strict lower triangle untouched. */
      gsl_matrix_set_zero (ncm_matrix_gsl (mean_cov));

      for (j = 0; j < n_kernels; j++)
      {
        ncm_matrix_memcpy (B, ncm_stats_dist_peek_cov_decomp (test->sd, j));

        for (p = 1; p < d; p++)
          for (q = 0; q < p; q++)
            ncm_matrix_set (B, p, q, 0.0);

        ncm_matrix_dgemm (mean_cov, 'T', 'N', 1.0 / (1.0 * n_kernels), B, B, 1.0);
      }

      /* lhs = A C A^T + kappa href^2 <Sigma'> */
      ncm_matrix_dgemm (AC, 'N', 'N', 1.0, A, C, 0.0);
      ncm_matrix_dgemm (lhs, 'N', 'T', 1.0, AC, A, 0.0);
      ncm_matrix_add_mul (lhs, kappa * href * href, mean_cov);

      for (p = 0; p < d; p++)
      {
        for (q = 0; q < d; q++)
        {
          max_C   = GSL_MAX (max_C, fabs (ncm_matrix_get (C, p, q)));
          max_dev = GSL_MAX (max_dev, fabs (ncm_matrix_get (lhs, p, q) - ncm_matrix_get (C, p, q)));
        }
      }

      g_assert_cmpfloat (max_dev, <, 1.0e-8 * max_C);

      /* det(A)^(1/d) is the reported scalar; for KDE with the sample covariance A = a I. */
      if (!NCM_IS_STATS_DIST_VKDE (test->sd) && (cov_type == NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE))
      {
        ncm_assert_cmpdouble_e (a, ==, 1.0 / sqrt (1.0 + kappa * href * href / (a * a)), 1.0e-10, 0.0);

        for (p = 0; p < d; p++)
          for (q = 0; q < d; q++)
            ncm_assert_cmpdouble_e (ncm_matrix_get (A, p, q), ==, (p == q) ? a : 0.0, 1.0e-10, 1.0e-12);
      }

      /* The centers are mu + A (x_i - mu). */
      g_assert_cmpuint (center_array->len, ==, n_kernels);

      for (i = 0; i < n_kernels; i++)
      {
        NcmVector *x_i = g_ptr_array_index (sample_array, i);
        NcmVector *c_i = g_ptr_array_index (center_array, i);
        NcmVector *dx  = ncm_vector_dup (x_i);
        gint ret;

        ncm_vector_sub (dx, mean);
        ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (A), ncm_vector_gsl (dx), 0.0, ncm_vector_gsl (y));
        g_assert_cmpint (ret, ==, 0);
        ncm_vector_add (y, mean);

        for (k = 0; k < d; k++)
          ncm_assert_cmpdouble_e (ncm_vector_get (c_i, k), ==, ncm_vector_get (y, k), 1.0e-10, 1.0e-10);

        ncm_vector_free (dx);
      }

      ncm_matrix_free (mean_cov);
      ncm_matrix_free (B);
      ncm_matrix_free (AC);
      ncm_matrix_free (lhs);
    }

    /* The density must still be finite and positive at the sample points. */
    for (i = 0; i < center_array->len; i++)
    {
      const gdouble p_i = ncm_stats_dist_eval (test->sd, g_ptr_array_index (sample_array, i));

      g_assert_true (gsl_finite (p_i));
      g_assert_cmpfloat (p_i, >, 0.0);
    }
  }

  /* Mixture covariance equals the sample covariance for a Gaussian kernel. That is a
   * statistical claim, so it belongs to the divergence lane. */
  if (test->divergence && NCM_IS_STATS_DIST_KERNEL_GAUSS (test->kernel))
  {
    NcmMatrix *cov_sample = ncm_stats_vec_peek_cov_matrix (sample_stats, 0);
    NcmMatrix *cov_est;

    for (i = 0; i < test->ntests; i++)
    {
      ncm_stats_dist_sample (test->sd, y, rng);
      ncm_stats_vec_append (test_stats, y, FALSE);
    }

    cov_est = ncm_stats_vec_peek_cov_matrix (test_stats, 0);

    g_assert_cmpfloat (ncm_matrix_cmp (cov_est, cov_sample, 1.0), <, 0.5);
    g_assert_cmpfloat (ncm_matrix_cmp_diag (cov_est, cov_sample, 1.0), <, 0.5);
  }

  /* Interpolated weights with shrunken centres. The rows of the interpolation matrix are
   * the sample points and the columns the centres, which no longer coincide: the whole
   * block has to be built instead of half of it. */
  {
    GPtrArray *sample_array = ncm_stats_dist_peek_sample_array (test->sd);

    ncm_stats_dist_prepare (test->sd, m2lnp_v);

    for (i = 0; i < sample_array->len; i++)
    {
      const gdouble m2lnp_i = ncm_stats_dist_eval_m2lnp (test->sd, g_ptr_array_index (sample_array, i));

      g_assert_true (gsl_finite (m2lnp_i));
    }
  }

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (y);
  ncm_vector_free (m2lnp_v);
  ncm_stats_vec_free (sample_stats);
  ncm_stats_vec_free (test_stats);
  ncm_mset_free (mset);
}

/* Defensive component: eps = 0 leaves eval, eval_m2lnp and sample untouched; eps > 0
 * gives q = (1 - eps) p + eps K with K the wide Student-t, eval and eval_m2lnp agree,
 * and the density stays positive far from the sample where the kernels alone vanish. */

/*
 * The batched evaluator against the scalar one, over more points than the VKDE puts in a
 * tile, so that a partial last tile is exercised too. The two are not required to agree
 * bit for bit: a tile shares one triangular solve across its points, which associates the
 * arithmetic differently from one solve per point. Run once with the defensive mixture
 * off and once with it on, since the mixture is applied by the wrapper rather than by
 * the subclass and so takes its own path in the batched version.
 */
static void
test_ncm_stats_dist_eval_vec (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  GPtrArray *x_a                 = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  const guint neval              = 300;
  NcmVector *m2lnp               = ncm_vector_new (neval);
  gulong N                       = 0;
  guint i;

  if (GPOINTER_TO_INT (pdata) == NCM_STATS_DIST_KDE_COV_TYPE_FIXED)
    ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd)));

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
    ncm_stats_dist_add_obs (test->sd, ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N));

  for (i = 0; i < neval; i++)
    g_ptr_array_add (x_a, ncm_vector_dup (ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N)));

  for (i = 0; i < 2; i++)
  {
    guint j;

    ncm_stats_dist_set_defensive_frac (test->sd, (i == 0) ? 0.0 : 0.05);
    ncm_stats_dist_prepare (test->sd, NULL);

    ncm_stats_dist_eval_m2lnp_vec (test->sd, x_a, m2lnp);

    for (j = 0; j < neval; j++)
    {
      const gdouble m2lnp_j = ncm_stats_dist_eval_m2lnp (test->sd, g_ptr_array_index (x_a, j));

      ncm_assert_cmpdouble_e (ncm_vector_get (m2lnp, j), ==, m2lnp_j, 1.0e-11, 0.0);
    }
  }

  ncm_vector_free (m2lnp);
  g_ptr_array_unref (x_a);
  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_defensive (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmStatsVec *sample_stats      = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_COV, FALSE);
  NcmVector *y                   = ncm_vector_new (test->dim);
  NcmVector *far                 = ncm_vector_new (test->dim);
  const gdouble eps              = 0.05;
  const gdouble scale            = 3.0;
  const gdouble nu               = 4.0;
  gulong N                       = 0;
  gdouble p0_far, p_far, m2lnp_far;
  guint i, k;

  if (GPOINTER_TO_INT (pdata) == NCM_STATS_DIST_KDE_COV_TYPE_FIXED)
    ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd)));

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y_i = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);

    ncm_stats_dist_add_obs (test->sd, y_i);
    ncm_stats_vec_append (sample_stats, y_i, FALSE);
  }

  /* Defaults: off. */
  g_assert_cmpfloat (ncm_stats_dist_get_defensive_frac (test->sd), ==, 0.0);
  g_assert_cmpfloat (ncm_stats_dist_get_defensive_scale (test->sd), ==, 4.0);
  g_assert_cmpfloat (ncm_stats_dist_get_defensive_nu (test->sd), ==, 3.0);

  ncm_stats_dist_prepare (test->sd, NULL);

  /* A point 40 sample standard deviations away along every axis. */
  {
    NcmVector *mean = ncm_stats_vec_peek_mean (sample_stats);

    for (k = 0; k < test->dim; k++)
      ncm_vector_set (far, k, ncm_vector_get (mean, k) + 40.0 * ncm_stats_vec_get_sd (sample_stats, k));
  }

  p0_far = ncm_stats_dist_eval (test->sd, far);

  /* eps = 0 must be the plain mixture, both evaluators. */
  {
    NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (test->sd);
    NcmVector *w                = ncm_stats_dist_peek_weights (test->sd);

    for (i = 0; i < 10; i++)
    {
      NcmVector *x_i = g_ptr_array_index (ncm_stats_dist_peek_sample_array (test->sd), i);

      g_assert_cmpfloat (ncm_stats_dist_eval (test->sd, x_i), ==, sd_class->eval_weights (test->sd, w, x_i));
      g_assert_cmpfloat (ncm_stats_dist_eval_m2lnp (test->sd, x_i), ==, sd_class->eval_weights_m2lnp (test->sd, w, x_i));
    }
  }

  ncm_stats_dist_set_defensive_frac (test->sd, eps);
  ncm_stats_dist_set_defensive_scale (test->sd, scale);
  ncm_stats_dist_set_defensive_nu (test->sd, nu);
  g_assert_cmpfloat (ncm_stats_dist_get_defensive_frac (test->sd), ==, eps);
  g_assert_cmpfloat (ncm_stats_dist_get_defensive_scale (test->sd), ==, scale);
  g_assert_cmpfloat (ncm_stats_dist_get_defensive_nu (test->sd), ==, nu);

  ncm_stats_dist_prepare (test->sd, NULL);

  /* q = (1 - eps) p + eps K, against an independent evaluation of K. */
  {
    NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (test->sd);
    NcmVector *w                = ncm_stats_dist_peek_weights (test->sd);
    NcmVector *mean             = ncm_stats_vec_peek_mean (sample_stats);
    NcmMatrix *C                = ncm_stats_vec_peek_cov_matrix (sample_stats, 0);
    NcmMatrix *U                = ncm_matrix_dup (C);
    NcmStatsDistKernelST *kst   = ncm_stats_dist_kernel_st_new (test->dim, nu);
    NcmVector *dx               = ncm_vector_new (test->dim);
    gdouble lnnorm;

    gsl_matrix_scale (ncm_matrix_gsl (U), scale);
    g_assert_cmpint (ncm_matrix_cholesky_decomp (U, 'U'), ==, 0);
    lnnorm = ncm_stats_dist_kernel_get_lnnorm (NCM_STATS_DIST_KERNEL (kst), U);

    for (i = 0; i < 10; i++)
    {
      NcmVector *x_i  = (i < 9) ? g_ptr_array_index (ncm_stats_dist_peek_sample_array (test->sd), i) : far;
      const gdouble p = sd_class->eval_weights (test->sd, w, x_i);
      gdouble chi2, K, q, m2lnq;

      ncm_vector_memcpy (dx, x_i);
      ncm_vector_sub (dx, mean);
      g_assert_cmpint (gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (U), ncm_vector_gsl (dx)), ==, 0);
      chi2  = ncm_vector_dot (dx, dx);
      K     = ncm_stats_dist_kernel_eval_unnorm (NCM_STATS_DIST_KERNEL (kst), chi2) / exp (lnnorm);
      q     = ncm_stats_dist_eval (test->sd, x_i);
      m2lnq = ncm_stats_dist_eval_m2lnp (test->sd, x_i);

      ncm_assert_cmpdouble_e (q, ==, (1.0 - eps) * p + eps * K, 1.0e-10, 0.0);
      ncm_assert_cmpdouble_e (m2lnq, ==, -2.0 * log (q), 1.0e-10, 0.0);
    }

    ncm_matrix_free (U);
    ncm_vector_free (dx);
    ncm_stats_dist_kernel_st_free (kst);
  }

  /* Far from the sample the wide component dominates and the density is positive. */
  p_far     = ncm_stats_dist_eval (test->sd, far);
  m2lnp_far = ncm_stats_dist_eval_m2lnp (test->sd, far);
  g_assert_true (gsl_finite (m2lnp_far));
  g_assert_cmpfloat (p_far, >, 0.0);
  g_assert_cmpfloat (p_far, >, p0_far);

  /* Draws: roughly eps of them come from the wide component, seen as a heavier tail in
   * the Mahalanobis distance; every draw is finite and the mean stays near the sample mean. */
  {
    NcmStatsVec *draw_stats = ncm_stats_vec_new (test->dim, NCM_STATS_VEC_COV, FALSE);
    NcmVector *mean         = ncm_stats_vec_peek_mean (sample_stats);
    const guint ndraws      = 2000;

    for (i = 0; i < ndraws; i++)
    {
      ncm_stats_dist_sample (test->sd, y, rng);

      for (k = 0; k < test->dim; k++)
        g_assert_true (gsl_finite (ncm_vector_get (y, k)));

      ncm_stats_vec_append (draw_stats, y, FALSE);
    }

    for (k = 0; k < test->dim; k++)
    {
      const gdouble sd_k = ncm_stats_vec_get_sd (sample_stats, k);

      g_assert_cmpfloat (fabs (ncm_stats_vec_get_mean (draw_stats, k) - ncm_vector_get (mean, k)), <, 0.5 * sd_k);
    }

    ncm_stats_vec_free (draw_stats);
  }

  /* A second preparation reuses the wide kernel: only its tail index changes, and the
   * density must follow. */
  {
    const gdouble p_nu = ncm_stats_dist_eval (test->sd, far);

    ncm_stats_dist_set_defensive_nu (test->sd, 1.0);
    ncm_stats_dist_prepare (test->sd, NULL);
    g_assert_cmpfloat (ncm_stats_dist_get_defensive_nu (test->sd), ==, 1.0);
    g_assert_cmpfloat (ncm_stats_dist_eval (test->sd, far), >, p_nu);

    ncm_stats_dist_set_defensive_nu (test->sd, nu);
  }

  /* Back to zero: the plain mixture again. */
  ncm_stats_dist_set_defensive_frac (test->sd, 0.0);
  ncm_stats_dist_prepare (test->sd, NULL);
  g_assert_cmpfloat (ncm_stats_dist_eval (test->sd, far), ==, p0_far);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (y);
  ncm_vector_free (far);
  ncm_stats_vec_free (sample_stats);
  ncm_mset_free (mset);
}

/* points-per-dim: the neighbor count is min (n, ceil (c d)). With c d >= n every local
 * covariance is the whole sample's, so all kernels share one scale matrix and it equals the
 * plain KDE's (the KDE limit); with c d < n the kernels differ. */
static void
test_ncm_stats_dist_vkde_points_per_dim (void)
{
  const guint d           = 3;
  const guint n           = 240;
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, 20260916);
  NcmStatsDistVKDE *vkde  = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistVKDE *vkde2 = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)), NCM_STATS_DIST_CV_NONE);
  NcmStatsDistKDE *kde    = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)), NCM_STATS_DIST_CV_NONE);
  NcmVector *y            = ncm_vector_new (d);
  guint i, j;

  g_assert_cmpfloat (ncm_stats_dist_vkde_get_points_per_dim (vkde), ==, 0.0);
  /* neighbor count: fraction when c = 0, min (n, ceil (c d)) otherwise */
  ncm_stats_dist_vkde_set_local_frac (vkde, 0.25);
  g_assert_cmpuint (ncm_stats_dist_vkde_get_n_neighbors (vkde, n), ==, 60);
  ncm_stats_dist_vkde_set_points_per_dim (vkde, 10.0);
  g_assert_cmpuint (ncm_stats_dist_vkde_get_n_neighbors (vkde, n), ==, 30);
  ncm_stats_dist_vkde_set_points_per_dim (vkde, 1000.0);
  g_assert_cmpuint (ncm_stats_dist_vkde_get_n_neighbors (vkde, n), ==, n);
  {
    gdouble c;

    g_object_get (G_OBJECT (vkde), "points-per-dim", &c, NULL);
    g_assert_cmpfloat (c, ==, 1000.0);
  }

  ncm_stats_dist_vkde_set_points_per_dim (vkde2, 10.0);

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < d; j++)
      ncm_vector_set (y, j, ncm_rng_gaussian_gen (rng, 0.0, 1.0 + j));

    ncm_stats_dist_add_obs (NCM_STATS_DIST (vkde), y);
    ncm_stats_dist_add_obs (NCM_STATS_DIST (vkde2), y);
    ncm_stats_dist_add_obs (NCM_STATS_DIST (kde), y);
  }

  ncm_stats_dist_prepare (NCM_STATS_DIST (vkde), NULL);
  ncm_stats_dist_prepare (NCM_STATS_DIST (vkde2), NULL);
  ncm_stats_dist_prepare (NCM_STATS_DIST (kde), NULL);

  /* KDE limit: every kernel's factor equals the KDE's global one. */
  {
    NcmMatrix *U_kde  = ncm_stats_dist_peek_cov_decomp (NCM_STATS_DIST (kde), 0);
    gboolean any_diff = FALSE;

    for (i = 0; i < ncm_stats_dist_get_n_kernels (NCM_STATS_DIST (vkde)); i++)
    {
      NcmMatrix *U_i = ncm_stats_dist_peek_cov_decomp (NCM_STATS_DIST (vkde), i);
      guint p, q;

      for (p = 0; p < d; p++)
        for (q = p; q < d; q++)
          ncm_assert_cmpdouble_e (ncm_matrix_get (U_i, p, q), ==, ncm_matrix_get (U_kde, p, q), 1.0e-10, 1.0e-12);
    }

    /* c d < n: local covariances differ between kernels. */
    for (i = 1; i < ncm_stats_dist_get_n_kernels (NCM_STATS_DIST (vkde2)) && !any_diff; i++)
      any_diff = ncm_matrix_cmp (ncm_stats_dist_peek_cov_decomp (NCM_STATS_DIST (vkde2), i), ncm_stats_dist_peek_cov_decomp (NCM_STATS_DIST (vkde2), 0), 0.0) > 1.0e-6;

    g_assert_true (any_diff);
  }

  ncm_vector_free (y);
  ncm_stats_dist_vkde_free (vkde);
  ncm_stats_dist_vkde_free (vkde2);
  ncm_stats_dist_kde_free (kde);
  ncm_rng_free (rng);
}

/* The batched leave-one-out sweep must reproduce the point-by-point one it replaces: drop
 * one kernel's weight, evaluate that point, restore. The batch reaches the same numbers by
 * blanking one entry of the kernel sum instead, over a sweep shared by every point. */
static void
test_ncm_stats_dist_loo_batch (void)
{
  const guint d    = 4;
  const guint np   = 250;
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, 20260924);
  NcmVector *m2lnL = ncm_vector_new (np);
  GPtrArray *x_a   = g_ptr_array_new ();
  guint t;

  for (t = 0; t < 2; t++)
  {
    NcmStatsDistKernel *kernel = (t == 0) ?
                                 NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (d, 3.0)) :
                                 NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d));
    NcmStatsDist *sd            = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
    NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);
    NcmVector *weights, *loo_weights, *ref, *batch;
    guint i, j, n;

    for (i = 0; i < np; i++)
    {
      NcmVector *y   = ncm_vector_new (d);
      gdouble chi2_i = 0.0;

      for (j = 0; j < d; j++)
      {
        const gdouble y_j = ncm_rng_gaussian_gen (rng, 0.0, 1.0);

        ncm_vector_set (y, j, y_j);
        chi2_i += y_j * y_j;
      }

      ncm_vector_set (m2lnL, i, chi2_i);
      ncm_stats_dist_add_obs (sd, y);
      ncm_vector_free (y);
    }

    ncm_stats_dist_prepare (sd, m2lnL);

    n           = ncm_stats_dist_get_n_kernels (sd);
    weights     = ncm_stats_dist_peek_weights (sd);
    loo_weights = ncm_vector_dup (weights);
    ref         = ncm_vector_new (n);
    batch       = ncm_vector_new (n);

    g_ptr_array_set_size (x_a, 0);

    for (i = 0; i < n; i++)
      g_ptr_array_add (x_a, g_ptr_array_index (ncm_stats_dist_peek_sample_array (sd), i));

    for (i = 0; i < n; i++)
    {
      const gdouble w_i = ncm_vector_get (weights, i);

      ncm_vector_set (loo_weights, i, 0.0);
      ncm_vector_set (ref, i, sd_class->eval_weights_m2lnp (sd, loo_weights, g_ptr_array_index (x_a, i)));
      ncm_vector_set (loo_weights, i, w_i);
    }

    sd_class->eval_weights_m2lnp_loo (sd, weights, x_a, batch);

    for (i = 0; i < n; i++)
    {
      ncm_assert_cmpdouble_e (ncm_vector_get (batch, i), ==, ncm_vector_get (ref, i), 1.0e-12, 0.0);
      g_assert_true (gsl_finite (ncm_vector_get (batch, i)));
    }

    ncm_vector_free (batch);
    ncm_vector_free (ref);
    ncm_vector_free (loo_weights);
    ncm_stats_dist_free (sd);
    ncm_stats_dist_kernel_free (kernel);
  }

  g_ptr_array_unref (x_a);
  ncm_vector_free (m2lnL);
  ncm_rng_free (rng);
}

/* The three bandwidth objectives on the same Gaussian sample: out-of-sample -2lnq
 * (SPLIT_M2LNP), out-of-sample acceptance (SPLIT_ACCEPT) and leave-one-out likelihood
 * (LOO_M2LNP) must each return a finite bandwidth inside the search range, agree within a
 * factor of three, and leave a normalized, positive density at the sample points. */
static void
test_ncm_stats_dist_cv_objectives (void)
{
  const guint d              = 3;
  const guint n              = 300;
  const NcmStatsDistCV cv[3] = {NCM_STATS_DIST_CV_SPLIT_M2LNP, NCM_STATS_DIST_CV_SPLIT_ACCEPT, NCM_STATS_DIST_CV_LOO_M2LNP};
  NcmRNG *rng                = ncm_rng_seeded_new (NULL, 20260916);
  NcmVector *x               = ncm_vector_new (d);
  NcmVector *m2lnL           = ncm_vector_new (n);
  GPtrArray *sample          = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  gdouble h[3];
  guint i, j, c;

  for (i = 0; i < n; i++)
  {
    NcmVector *y = ncm_vector_new (d);
    gdouble chi2 = 0.0;

    for (j = 0; j < d; j++)
    {
      const gdouble z = ncm_rng_ugaussian_gen (rng);

      ncm_vector_set (y, j, (1.0 + j) * z);
      chi2 += z * z;
    }

    ncm_vector_set (m2lnL, i, chi2);
    g_ptr_array_add (sample, y);
  }

  for (c = 0; c < 3; c++)
  {
    NcmStatsDist *sd = NCM_STATS_DIST (ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)), cv[c]));

    ncm_stats_dist_set_split_frac (sd, 0.8);
    ncm_stats_dist_set_over_smooth (sd, 1.0);
    ncm_stats_dist_set_uniform_weights (sd, TRUE);

    for (i = 0; i < n; i++)
      ncm_stats_dist_add_obs (sd, g_ptr_array_index (sample, i));

    ncm_stats_dist_prepare (sd, m2lnL);

    /* uniform-weights: no NNLS, every kernel weight is 1 / n_kernels */
    {
      NcmVector *w = ncm_stats_dist_peek_weights (sd);

      for (i = 0; i < ncm_stats_dist_get_n_kernels (sd); i++)
        ncm_assert_cmpdouble_e (ncm_vector_get (w, i), ==, 1.0 / ncm_stats_dist_get_n_kernels (sd), 1.0e-12, 0.0);
    }

    h[c] = ncm_stats_dist_get_over_smooth (sd);
    g_assert_true (gsl_finite (h[c]));
    g_assert_cmpfloat (h[c], >, 0.05);
    g_assert_cmpfloat (h[c], <, 20.0);

    for (i = 0; i < 20; i++)
    {
      const gdouble p = ncm_stats_dist_eval (sd, g_ptr_array_index (sample, i));

      g_assert_true (gsl_finite (p));
      g_assert_cmpfloat (p, >, 0.0);
    }

    ncm_stats_dist_free (sd);
  }

  g_assert_cmpfloat (h[1] / h[0], <, 3.0);
  g_assert_cmpfloat (h[1] / h[0], >, 1.0 / 3.0);
  g_assert_cmpfloat (h[2] / h[0], <, 3.0);
  g_assert_cmpfloat (h[2] / h[0], >, 1.0 / 3.0);

  g_ptr_array_unref (sample);
  ncm_vector_free (x);
  ncm_vector_free (m2lnL);
  ncm_rng_free (rng);
}

/* Auto-kernel tunes a Student-t kernel in place: a Gaussian one must be refused at prepare. */
static void
test_ncm_stats_dist_invalid_auto_kernel_gauss (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmStatsDist *sd = NCM_STATS_DIST (ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (2)), NCM_STATS_DIST_CV_SPLIT_M2LNP));
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, 5);
  guint i;

  ncm_stats_dist_set_auto_kernel (sd, TRUE);

  for (i = 0; i < 100; i++)
  {
    NcmVector *y = ncm_vector_new (2);

    ncm_vector_set (y, 0, ncm_rng_ugaussian_gen (rng));
    ncm_vector_set (y, 1, ncm_rng_ugaussian_gen (rng));
    ncm_stats_dist_add_obs (sd, y);
    ncm_vector_free (y);
  }

  /* Runs only as the subprocess of test_ncm_stats_dist_traps (); this aborts. */
  ncm_stats_dist_prepare (sd, NULL);
}

/* SPLIT_ACCEPT needs the sample's -2lnL: prepare () without it must be refused. */
static void
test_ncm_stats_dist_invalid_cv_accept (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmStatsDist *sd = NCM_STATS_DIST (ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (2)), NCM_STATS_DIST_CV_SPLIT_ACCEPT));
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, 5);
  guint i;

  for (i = 0; i < 100; i++)
  {
    NcmVector *y = ncm_vector_new (2);

    ncm_vector_set (y, 0, ncm_rng_ugaussian_gen (rng));
    ncm_vector_set (y, 1, ncm_rng_ugaussian_gen (rng));
    ncm_stats_dist_add_obs (sd, y);
    ncm_vector_free (y);
  }

  /* Runs only as the subprocess of test_ncm_stats_dist_traps (); this aborts. */
  ncm_stats_dist_prepare (sd, NULL);
}

static void
test_ncm_stats_dist_free (TestNcmStatsDist *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_dist_free, test->sd);
  NCM_TEST_FREE (ncm_stats_dist_kernel_free, test->kernel);
}

static void
test_ncm_stats_dist_sanity (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmStatsDist *sd = test->sd;

  ncm_stats_dist_ref (test->sd);
  ncm_stats_dist_free (test->sd);

  ncm_stats_dist_ref (test->sd);
  ncm_stats_dist_clear (&test->sd);

  g_assert_true (test->sd == NULL);

  test->sd = sd;

  ncm_stats_dist_set_use_threads (test->sd, TRUE);
  g_assert_true (ncm_stats_dist_get_use_threads (test->sd));

  ncm_stats_dist_set_use_threads (test->sd, FALSE);
  g_assert_true (!ncm_stats_dist_get_use_threads (test->sd));

  ncm_stats_dist_set_over_smooth (test->sd, 1.2);
  g_assert_cmpfloat (ncm_stats_dist_get_over_smooth (test->sd), ==, 1.2);

  ncm_stats_dist_set_split_frac (test->sd, 0.2);
  g_assert_cmpfloat (ncm_stats_dist_get_split_frac (test->sd), ==, 0.2);

  ncm_stats_dist_set_print_fit (test->sd, TRUE);
  g_assert_true (ncm_stats_dist_get_print_fit (test->sd));

  ncm_stats_dist_set_print_fit (test->sd, FALSE);
  g_assert_true (!ncm_stats_dist_get_print_fit (test->sd));

  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_NONE);
  g_assert_cmpint (ncm_stats_dist_get_cv_type (test->sd), ==, NCM_STATS_DIST_CV_NONE);

  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_SPLIT_ACCEPT);
  g_assert_cmpint (ncm_stats_dist_get_cv_type (test->sd), ==, NCM_STATS_DIST_CV_SPLIT_ACCEPT);

  {
    const guint dim                    = ncm_stats_dist_get_dim (test->sd);
    NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);
    NcmStatsDistKernel *kernel         = ncm_stats_dist_get_kernel (test->sd);

    g_assert_true (kernel == ncm_stats_dist_peek_kernel (test->sd));

    ncm_stats_dist_kernel_free (kernel);

    ncm_stats_dist_set_kernel (test->sd, NCM_STATS_DIST_KERNEL (sdk_gauss));
    g_assert_true (NCM_STATS_DIST_KERNEL (sdk_gauss) == ncm_stats_dist_get_kernel (test->sd));
  }
}

static void
test_ncm_stats_dist_prepare_too_few (TestNcmStatsDist *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
    NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
    NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (test->dim);
    NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
    gulong N                       = 0;

    ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

    {
      NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);

      ncm_stats_dist_add_obs (test->sd, y);
    }


    ncm_stats_dist_prepare (test->sd, NULL);
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the sample is too small*");
}

static void
test_ncm_stats_dist_cmp_dist (TestNcmStatsDist *test, NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, NcmRNG *rng)
{
  NcmStatsVec *err_stats = ncm_stats_vec_new (2, NCM_STATS_VEC_VAR, FALSE);
  gulong N               = 0;
  guint i;

  /* Mechanics: the estimator has to evaluate and sample without producing nonsense. That
   * it reproduces the distribution it was built from needs the full sample, and is the
   * divergence mode's claim. */
  if (!test->divergence)
  {
    for (i = 0; i < test->ntests; i++)
    {
      NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);

      g_assert_true (gsl_finite (ncm_stats_dist_eval_m2lnp (test->sd, y)));

      ncm_stats_dist_sample (test->sd, y, rng);
      g_assert_true (gsl_finite (ncm_stats_dist_eval_m2lnp (test->sd, y)));
    }

    ncm_stats_vec_free (err_stats);

    return;
  }

  for (i = 0; i < test->ntests; i++)
  {
    NcmVector *y;
    gdouble m2lnL, m2lnp_s, alpha0, alpha1;

    /* Measuring the Kullback-Leibler divergence */

    y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    m2lnp_s = ncm_stats_dist_eval_m2lnp (test->sd, y);
    alpha0  = 0.25 * (m2lnp_s - m2lnL);

    ncm_stats_dist_sample (test->sd, y, rng);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    m2lnp_s = ncm_stats_dist_eval_m2lnp (test->sd, y);
    alpha1  = 0.25 * (m2lnL - m2lnp_s);

    ncm_stats_vec_set (err_stats, 0, log1p (tanh (alpha0)));
    ncm_stats_vec_set (err_stats, 1, log1p (tanh (alpha1)));
    ncm_stats_vec_update (err_stats);
  }

  /* Measuring the Jensen-Shannon divergence */
  {
    gdouble KL_P_M  = ncm_stats_vec_get_mean (err_stats, 0);
    gdouble KL_PS_M = ncm_stats_vec_get_mean (err_stats, 1);
    gdouble JS_P_PS = 0.5 * (KL_P_M + KL_PS_M);

    /* printf ("JS(P||P_s) = % 22.15g\n", JS_P_PS); */
    g_assert_cmpfloat (JS_P_PS, <, RELTOL);
  }

  ncm_stats_vec_free (err_stats);
}

static void
test_ncm_stats_dist_dens_est (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  gulong N                        = 0;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_prepare (test->sd, NULL);

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
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  gulong N                        = 0;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_prepare (test->sd, m2lnp_v);

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
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  const guint ntests              = test->divergence ? 10000000 : SAMPLING_NTESTS_MECHANICS;
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  NcmVector *weights, *cum;
  gdouble sum_weights;
  guint i, n;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_prepare (test->sd, m2lnp_v);

  n   = ncm_stats_dist_get_sample_size (test->sd);
  cum = ncm_vector_new (n);

  ncm_vector_set_zero (cum);

  for (i = 0; i < ntests; i++)
  {
    const guint k_i = ncm_stats_dist_kernel_choose (test->sd, rng);

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

    /* A kernel of zero weight must never be drawn, whatever the draw count. Matching the
     * frequencies to 10% is what needs ten million draws, so it is asserted only where
     * that many are taken. */
    if (w_i == 0.0)
      g_assert_true (c_i == 0.0);
    else if (test->divergence)
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
test_ncm_stats_dist_dens_interp_cv_split_m2lnp (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_SPLIT_M2LNP);
  ncm_stats_dist_prepare (test->sd, m2lnp_v);

  test_ncm_stats_dist_cmp_dist (test, data_mvnd, mset, rng);

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}

static void
test_ncm_stats_dist_dens_interp_cv_loo (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  gulong N                        = 0;
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  /* ncm_stats_dist_set_print_fit (test->sd, TRUE); */
  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_LOO);
  ncm_stats_dist_prepare (test->sd, m2lnp_v);

  /* The closed-form objective evaluates the interpolation matrix at sqrt(2) h on the way;
   * the bandwidth the object keeps must be h itself. */
  if (!ncm_stats_dist_get_center_shrink (test->sd))
    ncm_assert_cmpdouble_e (ncm_stats_dist_get_href (test->sd), ==,
                            NCM_STATS_DIST_GET_CLASS (test->sd)->bandwidth (test->sd), 1.0e-12, 0.0);

  test_ncm_stats_dist_cmp_dist (test, data_mvnd, mset, rng);

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
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
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

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_prepare (test->sd, NULL);

  for (i = 0; i < test->ntests; i++)
  {
    ncm_stats_dist_sample (test->sd, y, rng);
    ncm_stats_vec_append (test_stats, y, FALSE);
  }

  {
    NcmMatrix *cov_est    = ncm_stats_vec_peek_cov_matrix (test_stats, 0);
    NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
    NcmMatrix *cov        = ncm_data_gauss_cov_peek_cov (gcov);

    /* Mechanics: the sampler has to produce a usable covariance. That it recovers the
     * one it was built from is the divergence mode's claim -- and the retry below
     * rebuilds the fixture, so it is not something to run in a lane that has to be
     * deterministic. */
    if (!test->divergence)
    {
      g_assert_true (gsl_finite (ncm_matrix_get (cov_est, 0, 0)));
    }
    else if (
      (test->nfail < 10) &&
      ((ncm_matrix_cmp (cov_est, cov, 0.0) >= 0.5) ||
       (ncm_matrix_cmp_diag (cov_est, cov, 0.0) >= 0.5))
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
      g_assert_cmpfloat (ncm_matrix_cmp (cov_est, cov, 1.0), <, 0.5);
      g_assert_cmpfloat (ncm_matrix_cmp_diag (cov_est, cov, 1.0), <, 0.5);
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
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
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

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_prepare (test->sd, m2lnp_v);
  ncm_stats_dist_prepare (sd_dup, m2lnp_v);

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
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  gulong N                        = 0;
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  guint i;

  switch (cov_type)
  {
    case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
    {
      NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);

      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (gcov));
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

  ncm_stats_dist_prepare (test->sd, NULL);
  ncm_stats_dist_prepare_shapes (test->sd, ncm_stats_dist_peek_sample_array (test->sd));
  ncm_stats_dist_prepare (test->sd, m2lnp_v);

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

    /* The whole-sample covariance and its Cholesky factor: U^T U = C over the upper
     * triangle, the only part of U the decomposition writes. */
    {
      NcmMatrix *full_cov = ncm_stats_dist_peek_full_cov (test->sd);
      NcmMatrix *U        = ncm_stats_dist_peek_full_cov_decomp (test->sd);
      guint a, b, k;

      g_assert_cmpuint (ncm_matrix_nrows (full_cov), ==, dim);
      g_assert_cmpuint (ncm_matrix_ncols (full_cov), ==, dim);
      g_assert_cmpuint (ncm_matrix_nrows (U), ==, dim);
      g_assert_cmpuint (ncm_matrix_ncols (U), ==, dim);

      for (a = 0; a < dim; a++)
      {
        for (b = a; b < dim; b++)
        {
          gdouble UtU_ab = 0.0;

          for (k = 0; k <= a; k++)
            UtU_ab += ncm_matrix_get (U, k, a) * ncm_matrix_get (U, k, b);

          ncm_assert_cmpdouble_e (UtU_ab, ==, ncm_matrix_get (full_cov, a, b), 1.0e-10, 1.0e-14);
        }
      }
    }

    for (i = 0; i < n; i++)
    {
      NcmMatrix *cov_decomp   = ncm_stats_dist_peek_cov_decomp (test->sd, i);
      const gdouble lnnorm_i  = ncm_stats_dist_get_lnnorm (test->sd, i);
      const gdouble lnnorm0_i = ncm_stats_dist_kernel_get_lnnorm (kernel, cov_decomp);
      NcmVector *y_i          = NULL;
      NcmMatrix *cov_i        = NULL;
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
  g_test_trap_subprocess ("/ncm/stats/dist/nd/kde/gauss/invalid/stub/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats/dist/nd/vkde/cauchy/invalid/center_shrink/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats/dist/nd/kde/gauss/invalid/cv_accept_without_m2lnL/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NCM_STATS_DIST_CV_SPLIT_ACCEPT needs the sample*");

  g_test_trap_subprocess ("/ncm/stats/dist/nd/kde/gauss/invalid/auto_kernel_gauss/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*needs an NcmStatsDistKernelST*");
}

static void
test_ncm_stats_dist_invalid_stub (TestNcmStatsDist *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

/* Centre shrinkage needs a kernel with a finite covariance; the Cauchy kernel
 * (Student-t with nu = 1) has none and must be refused. */
static void
test_ncm_stats_dist_invalid_center_shrink (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmStatsDistKernelST *sdk_st = ncm_stats_dist_kernel_st_new (2, 1.0);
  NcmStatsDistVKDE *sdvkde     = ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (sdk_st), NCM_STATS_DIST_CV_NONE);
  NcmStatsDist *sd             = NCM_STATS_DIST (sdvkde);
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, 123);
  guint i;

  g_assert_cmpint (gsl_isinf (ncm_stats_dist_kernel_get_var_factor (NCM_STATS_DIST_KERNEL (sdk_st))), ==, 1);

  ncm_stats_dist_set_center_shrink (sd, TRUE);
  ncm_stats_dist_vkde_set_local_frac (sdvkde, 0.5);

  for (i = 0; i < 100; i++)
  {
    NcmVector *y = ncm_vector_new (2);

    ncm_vector_set (y, 0, ncm_rng_ugaussian_gen (rng));
    ncm_vector_set (y, 1, ncm_rng_ugaussian_gen (rng));
    ncm_stats_dist_add_obs (sd, y);
    ncm_vector_free (y);
  }

  ncm_stats_dist_prepare (sd, NULL);

  g_assert_not_reached ();
}

static void
test_ncm_stats_dist_split_underflowing_m2lnp (void)
{
  /*
   * A cross-validation that splits the sample leaves n_kernels < n_obs. When the
   * posterior values span more than the double range, prepare() takes its cut
   * path, which used to sort n_obs indices into an array holding n_kernels of them and
   * corrupt the heap. Everything in that path is indexed by kernel, and the branch that
   * rebuilds the sample needs its own count over the observations.
   */
  const guint d                   = 4;
  const guint n_obs               = 2000;
  NcmStatsDistKernelGauss *kernel = ncm_stats_dist_kernel_gauss_new (d);
  NcmStatsDist *sd                = NCM_STATS_DIST (ncm_stats_dist_vkde_new (NCM_STATS_DIST_KERNEL (kernel),
                                                                             NCM_STATS_DIST_CV_SPLIT_M2LNP));
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, 42);
  NcmVector *m2lnp = ncm_vector_new (n_obs);
  NcmVector *x     = ncm_vector_new (d);
  guint i, k;

  ncm_stats_dist_set_split_frac (sd, 0.1);
  ncm_stats_dist_set_over_smooth (sd, 1.0);
  ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (sd), 0.4);

  for (i = 0; i < n_obs; i++)
  {
    gdouble m2lnp_i = 0.0;

    for (k = 0; k < d; k++)
    {
      const gdouble x_k = ncm_rng_gaussian_gen (rng, 0.0, 1.0);

      ncm_vector_set (x, k, x_k);
      m2lnp_i += x_k * x_k;
    }

    ncm_stats_dist_add_obs (sd, x);

    /*
     * A spread far beyond -2 * 2 * ln (DBL_EPSILON) ~ 144, alternating so that the kernel
     * block itself spans the range and the cut path is entered. Half the kernels sit
     * inside it, so the path taken is the one that reweights and returns; the sort that
     * overflowed runs before either branch.
     */
    ncm_vector_set (m2lnp, i, ((i % 2) == 0) ? m2lnp_i : m2lnp_i + 1.0e3);
  }

  ncm_stats_dist_prepare (sd, m2lnp);

  g_assert_cmpuint (ncm_stats_dist_get_n_kernels (sd), >, 0);
  g_assert_cmpuint (ncm_vector_len (ncm_stats_dist_peek_weights (sd)), ==,
                    ncm_stats_dist_get_n_kernels (sd));

  ncm_vector_free (x);
  ncm_vector_free (m2lnp);
  ncm_rng_free (rng);
  ncm_stats_dist_free (sd);
  ncm_stats_dist_kernel_gauss_free (kernel);
}

static void
test_ncm_stats_dist_split_drop_far_points (void)
{
  /*
   * The other branch of the split: at least half the kernel centres are within
   * NCM_STATS_DIST_M2LNL_RANGE of the best point, so the observations beyond it leave the
   * sample and everything is built from what remains. Here the far points are the last
   * tenth, outside the kernel block, and the sample must shrink to the kept nine tenths.
   */
  const guint d      = 3;
  const guint n      = 400;
  const guint n_far  = 40;
  const guint n_keep = n - n_far;
  NcmStatsDist *sd   = NCM_STATS_DIST (ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)),
                                                               NCM_STATS_DIST_CV_SPLIT_M2LNP));
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 20260923);
  NcmVector *m2lnL  = ncm_vector_new (n);
  GPtrArray *sample = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  guint i, j;

  ncm_stats_dist_set_split_frac (sd, 0.5);
  ncm_stats_dist_set_over_smooth (sd, 1.0);
  ncm_stats_dist_set_uniform_weights (sd, FALSE);

  for (i = 0; i < n; i++)
  {
    NcmVector *y = ncm_vector_new (d);
    gdouble chi2 = 0.0;

    for (j = 0; j < d; j++)
    {
      const gdouble z = ncm_rng_ugaussian_gen (rng);

      ncm_vector_set (y, j, z);
      chi2 += z * z;
    }

    ncm_vector_set (m2lnL, i, (i < n_keep) ? chi2 : chi2 + 1.0e3);
    g_ptr_array_add (sample, y);
    ncm_stats_dist_add_obs (sd, y);
  }

  ncm_stats_dist_prepare (sd, m2lnL);

  g_assert_cmpuint (ncm_stats_dist_get_sample_size (sd), ==, n_keep);
  g_assert_cmpuint (ncm_stats_dist_get_n_kernels (sd), ==, (guint) ceil (0.5 * n_keep));
  g_assert_cmpuint (ncm_vector_len (ncm_stats_dist_peek_weights (sd)), ==, ncm_stats_dist_get_n_kernels (sd));

  {
    GPtrArray *kept = ncm_stats_dist_peek_sample_array (sd);

    g_assert_cmpuint (kept->len, ==, n_keep);

    for (i = 0; i < n_keep; i++)
      for (j = 0; j < d; j++)
        g_assert_cmpfloat (ncm_vector_get (g_ptr_array_index (kept, i), j), ==, ncm_vector_get (g_ptr_array_index (sample, i), j));
  }

  for (i = 0; i < n; i += 20)
  {
    const gdouble p = ncm_stats_dist_eval (sd, g_ptr_array_index (sample, i));

    g_assert_true (gsl_finite (p));
    g_assert_cmpfloat (p, >=, 0.0);
  }

  g_ptr_array_unref (sample);
  ncm_vector_free (m2lnL);
  ncm_rng_free (rng);
  ncm_stats_dist_free (sd);
}

static void
_test_ncm_stats_dist_count_fit_message (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  guint *count = user_data;

  if (g_str_has_prefix (message, "# over-smooth:"))
    (*count)++;
}

static void
test_ncm_stats_dist_print_fit (void)
{
  /*
   * print-fit reports every objective evaluation as a "# over-smooth:" message, for
   * each cross-validation objective, both estimators and both AMISE routes (the closed
   * form of the Gaussian kernel, the Monte Carlo estimate of any other). The messages are
   * counted through a handler instead of being printed.
   */
  const guint d              = 2;
  const guint n              = 120;
  const NcmStatsDistCV cv[4] = {NCM_STATS_DIST_CV_SPLIT_M2LNP, NCM_STATS_DIST_CV_SPLIT_ACCEPT, NCM_STATS_DIST_CV_LOO_M2LNP, NCM_STATS_DIST_CV_LOO};
  NcmRNG *rng                = ncm_rng_seeded_new (NULL, 20260924);
  NcmVector *m2lnL           = ncm_vector_new (n);
  GPtrArray *sample          = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);
  guint i, j, c, m, k;

  for (i = 0; i < n; i++)
  {
    NcmVector *y = ncm_vector_new (d);
    gdouble chi2 = 0.0;

    for (j = 0; j < d; j++)
    {
      const gdouble z = ncm_rng_ugaussian_gen (rng);

      ncm_vector_set (y, j, z);
      chi2 += z * z;
    }

    ncm_vector_set (m2lnL, i, chi2);
    g_ptr_array_add (sample, y);
  }

  for (m = 0; m < 2; m++)
  {
    for (k = 0; k < 2; k++)
    {
      for (c = 0; c < 4; c++)
      {
        NcmStatsDistKernel *kernel = (k == 0) ? NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)) : NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (d, 3.0));
        NcmStatsDist *sd           = (m == 0) ? NCM_STATS_DIST (ncm_stats_dist_kde_new (kernel, cv[c])) : NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, cv[c]));
        guint count                = 0;
        guint handler_id;

        ncm_stats_dist_set_split_frac (sd, 0.5);
        ncm_stats_dist_set_over_smooth (sd, 1.0);
        ncm_stats_dist_set_print_fit (sd, TRUE);
        g_assert_true (ncm_stats_dist_get_print_fit (sd));

        for (i = 0; i < n; i++)
          ncm_stats_dist_add_obs (sd, g_ptr_array_index (sample, i));

        handler_id = g_log_set_handler ("NUMCOSMO", G_LOG_LEVEL_MESSAGE, &_test_ncm_stats_dist_count_fit_message, &count);
        ncm_stats_dist_prepare (sd, m2lnL);
        g_log_remove_handler ("NUMCOSMO", handler_id);

        g_assert_cmpuint (count, >, 1);
        g_assert_true (gsl_finite (ncm_stats_dist_eval (sd, g_ptr_array_index (sample, 0))));

        ncm_stats_dist_free (sd);
        ncm_stats_dist_kernel_free (kernel);
      }
    }
  }

  g_ptr_array_unref (sample);
  ncm_vector_free (m2lnL);
  ncm_rng_free (rng);
}

static void
test_ncm_stats_dist_kde_cov_fixed_nearPD (void)
{
  /*
   * A fixed covariance that is not positive definite is repaired by nearPD with a
   * warning, both when it is set and when the covariance type is switched to FIXED with
   * one already in place. The repaired factor is usable: the density evaluates.
   */
  const guint d          = 2;
  const guint n          = 100;
  NcmStatsDistKDE *sdkde = ncm_stats_dist_kde_new (NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (d)), NCM_STATS_DIST_CV_NONE);
  NcmStatsDist *sd       = NCM_STATS_DIST (sdkde);
  NcmMatrix *cov         = ncm_matrix_new (d, d);
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, 20260925);
  NcmVector *x           = ncm_vector_new (d);
  guint i, j;

  ncm_matrix_set (cov, 0, 0, 1.0);
  ncm_matrix_set (cov, 1, 1, 1.0);
  ncm_matrix_set (cov, 0, 1, 1.0 + 1.0e-3);
  ncm_matrix_set (cov, 1, 0, 1.0 + 1.0e-3);

  ncm_stats_dist_kde_set_cov_type (sdkde, NCM_STATS_DIST_KDE_COV_TYPE_FIXED);
  g_assert_cmpint (ncm_stats_dist_kde_get_cov_type (sdkde), ==, NCM_STATS_DIST_KDE_COV_TYPE_FIXED);

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*the fixed covariance was not positive definite*");
  ncm_stats_dist_kde_set_cov_fixed (sdkde, cov);
  g_test_assert_expected_messages ();

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*the fixed covariance was not positive definite*");
  ncm_stats_dist_kde_set_cov_type (sdkde, NCM_STATS_DIST_KDE_COV_TYPE_FIXED);
  g_test_assert_expected_messages ();

  {
    NcmMatrix *cov_fixed = ncm_stats_dist_kde_peek_cov_fixed (sdkde);

    g_assert_true (cov_fixed != cov);
    g_assert_cmpfloat (ncm_matrix_get (cov_fixed, 0, 1), ==, 1.0 + 1.0e-3);
  }

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < d; j++)
      ncm_vector_set (x, j, ncm_rng_ugaussian_gen (rng));

    ncm_stats_dist_add_obs (sd, x);
  }

  ncm_stats_dist_prepare (sd, NULL);

  {
    const gdouble p = ncm_stats_dist_eval (sd, x);

    g_assert_true (gsl_finite (p));
    g_assert_cmpfloat (p, >, 0.0);
  }

  ncm_vector_free (x);
  ncm_rng_free (rng);
  ncm_matrix_free (cov);
  ncm_stats_dist_free (sd);
}

void
test_ncm_stats_dist_cv_auto_kernel (TestNcmStatsDist *test, gconstpointer pdata)
{
  NcmRNG *rng                     = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmDataGaussCovMVND *data_mvnd  = ncm_data_gauss_cov_mvnd_new_full (test->dim, 1.0e-2, 5.0e-2, test->corr_level, 1.0, 2.0, rng);
  NcmModelMVND *model_mvnd        = ncm_model_mvnd_new (test->dim);
  NcmMSet *mset                   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmStatsDistKDECovType cov_type = GPOINTER_TO_INT (pdata);
  NcmVector *m2lnp_v              = ncm_vector_new (test->np);
  gulong N                        = 0;
  guint i;

  if (cov_type == NCM_STATS_DIST_KDE_COV_TYPE_FIXED)
    ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (test->sd), ncm_data_gauss_cov_peek_cov (NCM_DATA_GAUSS_COV (data_mvnd)));

  ncm_mset_param_set_vector (mset, ncm_data_gauss_cov_mvnd_peek_mean (data_mvnd));

  for (i = 0; i < test->np; i++)
  {
    NcmVector *y = ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, NULL, NULL, rng, &N);
    gdouble m2lnL;

    ncm_stats_dist_add_obs (test->sd, y);
    ncm_data_m2lnL_val (NCM_DATA (data_mvnd), mset, &m2lnL);
    ncm_vector_set (m2lnp_v, i, m2lnL);
  }

  /* With auto-kernel the bandwidth and the kernel tail index are fitted together, one
   * objective per cross-validation. Both objectives must leave a usable estimator and
   * report what they chose; print-fit is where that report comes out. The fit tunes the
   * Student-t kernel in place, so it is only run for those constructors. */
  ncm_stats_dist_set_auto_kernel (test->sd, TRUE);
  ncm_stats_dist_set_print_fit (test->sd, TRUE);
  g_assert_true (ncm_stats_dist_get_auto_kernel (test->sd));
  g_assert_true (ncm_stats_dist_get_print_fit (test->sd));

  if (NCM_IS_STATS_DIST_KERNEL_ST (test->kernel))
  {
    const NcmStatsDistCV cv_type[] = {NCM_STATS_DIST_CV_SPLIT_ACCEPT, NCM_STATS_DIST_CV_LOO_M2LNP};
    guint j;

    for (j = 0; j < G_N_ELEMENTS (cv_type); j++)
    {
      ncm_stats_dist_set_cv_type (test->sd, cv_type[j]);
      ncm_stats_dist_prepare (test->sd, m2lnp_v);

      g_assert_cmpfloat (ncm_stats_dist_get_over_smooth (test->sd), >, 0.0);
      g_assert_true (gsl_finite (ncm_stats_dist_get_href (test->sd)));

      for (i = 0; i < ncm_stats_dist_get_n_kernels (test->sd); i++)
      {
        NcmVector *x_i = g_ptr_array_index (ncm_stats_dist_peek_sample_array (test->sd), i);

        g_assert_true (gsl_finite (ncm_stats_dist_eval_m2lnp (test->sd, x_i)));
      }
    }
  }

  ncm_stats_dist_set_print_fit (test->sd, FALSE);

  /* Leave-one-out estimates its objective by Monte Carlo over antithetic pairs, which is
   * the one place the wide defensive component is drawn from outside ncm_stats_dist_sample ().
   */
  ncm_stats_dist_set_auto_kernel (test->sd, FALSE);
  ncm_stats_dist_set_defensive_frac (test->sd, 0.05);
  ncm_stats_dist_set_cv_type (test->sd, NCM_STATS_DIST_CV_LOO);
  ncm_stats_dist_prepare (test->sd, m2lnp_v);

  for (i = 0; i < ncm_stats_dist_get_n_kernels (test->sd); i++)
  {
    NcmVector *x_i = g_ptr_array_index (ncm_stats_dist_peek_sample_array (test->sd), i);

    g_assert_true (gsl_finite (ncm_stats_dist_eval_m2lnp (test->sd, x_i)));
  }

  ncm_stats_dist_set_defensive_frac (test->sd, 0.0);

  /* The rule-of-thumb bandwidth is derived for an estimator that sees the whole sample;
   * a local one sees only its neighborhood, so asking for it rescales the bandwidth by
   * the ratio of the two counts. */
  if (NCM_IS_STATS_DIST_VKDE (test->sd))
  {
    NcmStatsDistVKDE *sdvkde = NCM_STATS_DIST_VKDE (test->sd);

    g_assert_false (ncm_stats_dist_vkde_get_use_rot_href (sdvkde));

    ncm_stats_dist_vkde_set_use_rot_href (sdvkde, TRUE);
    g_assert_true (ncm_stats_dist_vkde_get_use_rot_href (sdvkde));

    ncm_stats_dist_prepare (test->sd, m2lnp_v);

    g_assert_true (gsl_finite (ncm_stats_dist_get_href (test->sd)));
    g_assert_cmpfloat (ncm_stats_dist_get_href (test->sd), >, 0.0);

    ncm_stats_dist_vkde_set_use_rot_href (sdvkde, FALSE);
  }

  ncm_model_mvnd_free (model_mvnd);
  ncm_data_gauss_cov_mvnd_free (data_mvnd);
  ncm_rng_free (rng);
  ncm_vector_free (m2lnp_v);
  ncm_mset_free (mset);
}


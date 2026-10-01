/***************************************************************************
 *            test_ncm_stats_dist_kernel.c
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

typedef enum _NcmStatsDistKernelType
{
  NCM_STATS_DIST_KERNEL_TYPE_GAUSS,
  NCM_STATS_DIST_KERNEL_TYPE_ST3,
} NcmStatsDistKernelType;

typedef struct _TestNcmStatsDistKernel
{
  NcmStatsDistKernel *kernel;
  NcmStatsDistKernelType kernel_type;
  gdouble nu;
  guint dim;
} TestNcmStatsDistKernel;

/* A subclass that overrides nothing, to reach the aborting defaults. */
#define TEST_TYPE_NCM_STATS_DIST_KERNEL_BARE (test_ncm_stats_dist_kernel_bare_get_type ())

G_DECLARE_FINAL_TYPE (TestNcmStatsDistKernelBare, test_ncm_stats_dist_kernel_bare, TEST, NCM_STATS_DIST_KERNEL_BARE, NcmStatsDistKernel)

struct _TestNcmStatsDistKernelBare
{
  NcmStatsDistKernel parent_instance;
};

G_DEFINE_TYPE (TestNcmStatsDistKernelBare, test_ncm_stats_dist_kernel_bare, NCM_TYPE_STATS_DIST_KERNEL)

static void
test_ncm_stats_dist_kernel_bare_init (TestNcmStatsDistKernelBare *bare)
{
}

static void
test_ncm_stats_dist_kernel_bare_class_init (TestNcmStatsDistKernelBareClass *klass)
{
}

static void test_ncm_stats_dist_kernel_new_st (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_new_gauss (TestNcmStatsDistKernel *test, gconstpointer pdata);

static void test_ncm_stats_dist_kernel_dim (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_bandwidth (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_norm (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_sum (TestNcmStatsDistKernel *test, gconstpointer pdata);
static void test_ncm_stats_dist_kernel_sample (TestNcmStatsDistKernel *test, gconstpointer pdata);

static void test_ncm_stats_dist_kernel_free (TestNcmStatsDistKernel *test, gconstpointer pdata);

static void test_ncm_stats_dist_kernel_st_nu_default (void);
static void test_ncm_stats_dist_kernel_traps (void);
static void test_ncm_stats_dist_kernel_unimplemented_subprocess (void);

typedef struct _TestNcmStatsDistKernelFunc
{
  const gchar *name;

  void (*test_func) (TestNcmStatsDistKernel *test, gconstpointer pdata);
} TestNcmStatsDistKernelFunc;

#define TEST_NCM_STATS_DIST_KERNEL_CONSTRUCTORS_LEN 2
#define TEST_NCM_STATS_DIST_KERNEL_TESTS_LEN 5

static TestNcmStatsDistKernelFunc constructors[TEST_NCM_STATS_DIST_KERNEL_CONSTRUCTORS_LEN] = {
  {"gauss", &test_ncm_stats_dist_kernel_new_gauss},
  {"st",    &test_ncm_stats_dist_kernel_new_st}
};

static TestNcmStatsDistKernelFunc tests[TEST_NCM_STATS_DIST_KERNEL_TESTS_LEN] = {
  {"dim",    &test_ncm_stats_dist_kernel_dim},
  {"band",   &test_ncm_stats_dist_kernel_bandwidth},
  {"norm",   &test_ncm_stats_dist_kernel_norm},
  {"sum",    &test_ncm_stats_dist_kernel_sum},
  {"sample", &test_ncm_stats_dist_kernel_sample},
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

  g_test_add_func ("/ncm/stats/dist/kernel/st/nu_default", &test_ncm_stats_dist_kernel_st_nu_default);
  g_test_add_func ("/ncm/stats/dist/kernel/traps", &test_ncm_stats_dist_kernel_traps);
  g_test_add_func ("/ncm/stats/dist/kernel/unimplemented/subprocess", &test_ncm_stats_dist_kernel_unimplemented_subprocess);

  g_test_run ();
}

static void
test_ncm_stats_dist_kernel_new_gauss (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  const guint dim                    = g_test_rand_int_range (1, 4);
  NcmStatsDistKernelGauss *sdk_gauss = ncm_stats_dist_kernel_gauss_new (dim);

  test->dim         = dim;
  test->kernel      = NCM_STATS_DIST_KERNEL (sdk_gauss);
  test->kernel_type = NCM_STATS_DIST_KERNEL_TYPE_GAUSS;

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
  const guint dim              = g_test_rand_int_range (1, 4);
  NcmStatsDistKernelST *sdk_st = ncm_stats_dist_kernel_st_new (dim, nu);

  test->dim         = dim;
  test->nu          = nu;
  test->kernel_type = NCM_STATS_DIST_KERNEL_TYPE_ST3;
  test->kernel      = NCM_STATS_DIST_KERNEL (sdk_st);

  g_assert_cmpfloat (ncm_stats_dist_kernel_st_get_nu (sdk_st), ==, nu);

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

/* Radial integrals of the unnormalized kernel g(s) = Kbar(s), s = r^2, at unit scale.
 * The derivatives are the closed forms of the two kernel shapes. */
typedef enum _TestRadial
{
  TEST_RADIAL_NORM,
  TEST_RADIAL_MOM2,
  TEST_RADIAL_RK,
  TEST_RADIAL_RLAP,
} TestRadial;

typedef struct _TestRadialArg
{
  TestNcmStatsDistKernel *test;
  TestRadial which;
} TestRadialArg;

static gdouble
_test_radial_integrand (gdouble r, gpointer userdata)
{
  TestRadialArg *arg           = userdata;
  TestNcmStatsDistKernel *test = arg->test;
  const gdouble d              = test->dim;
  const gdouble s              = r * r;
  const gdouble g              = ncm_stats_dist_kernel_eval_unnorm (test->kernel, s);
  const gdouble rdm1           = (test->dim == 1) ? 1.0 : gsl_pow_uint (r, test->dim - 1);
  gdouble dg, d2g;

  switch (test->kernel_type)
  {
    case NCM_STATS_DIST_KERNEL_TYPE_GAUSS:
      dg  = -0.5 * g;
      d2g = 0.25 * g;
      break;
    case NCM_STATS_DIST_KERNEL_TYPE_ST3:
    {
      const gdouble a = 0.5 * (test->nu + d);
      const gdouble q = 1.0 + s / test->nu;

      dg  = -a / test->nu * pow (q, -a - 1.0);
      d2g = a * (a + 1.0) / gsl_pow_2 (test->nu) * pow (q, -a - 2.0);
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  switch (arg->which)
  {
    case TEST_RADIAL_NORM:
      return rdm1 * g;

    case TEST_RADIAL_MOM2:
      return rdm1 * s * g;

    case TEST_RADIAL_RK:
      return rdm1 * g * g;

    case TEST_RADIAL_RLAP:
      /* Laplacian of g(r^2): 2 d g' + 4 s g'' */
      return rdm1 * gsl_pow_2 (2.0 * d * dg + 4.0 * s * d2g);

    default:
      g_assert_not_reached ();

      return 0.0;
  }
}

/* S_d int_0^inf r^(d-1) (...) dr, with S_d the area of the unit sphere */
static gdouble
_test_radial (TestNcmStatsDistKernel *test, TestRadial which)
{
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (1000);
  TestRadialArg arg             = {test, which};
  gsl_function F;
  gdouble res, err;

  F.function = &_test_radial_integrand;
  F.params   = &arg;

  gsl_integration_qagiu (&F, 0.0, 0.0, 1.0e-12, 1000, ws, &res, &err);
  gsl_integration_workspace_free (ws);

  return 2.0 * pow (M_PI, 0.5 * test->dim) / tgamma (0.5 * test->dim) * res;
}

static void
test_ncm_stats_dist_kernel_bandwidth (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  const gdouble n    = g_test_rand_double_range (1.0, 1.0e5);
  const gdouble h    = ncm_stats_dist_kernel_get_rot_bandwidth (test->kernel, n);
  const gdouble d    = test->dim;
  const gdouble I0   = _test_radial (test, TEST_RADIAL_NORM);
  const gdouble mu2  = _test_radial (test, TEST_RADIAL_MOM2) / (d * I0);
  const gdouble RK   = _test_radial (test, TEST_RADIAL_RK) / gsl_pow_2 (I0);
  const gdouble RLAP = _test_radial (test, TEST_RADIAL_RLAP) / gsl_pow_2 (I0);

  /* AMISE-optimal bandwidth when the estimated density is the kernel itself:
   * h^(d+4) = d R(K) / (mu_2^2 n R(Laplacian f)). */
  const gdouble h_amise = pow (d * RK / (gsl_pow_2 (mu2) * n * RLAP), 1.0 / (d + 4.0));

  ncm_assert_cmpdouble_e (ncm_stats_dist_kernel_get_var_factor (test->kernel), ==, mu2, 1.0e-10, 0.0);
  ncm_assert_cmpdouble_e (h, ==, h_amise, 1.0e-10, 0.0);

  if (test->kernel_type == NCM_STATS_DIST_KERNEL_TYPE_GAUSS)
  {
    ncm_assert_cmpdouble_e (h, ==, pow (4.0 / ((d + 2.0) * n), 1.0 / (d + 4.0)), 1.0e-15, 0.0);
  }
  else
  {
    NcmStatsDistKernelST *sdk_st = NCM_STATS_DIST_KERNEL_ST (test->kernel);
    gdouble h3;

    /* Below nu = 3 the rule is evaluated at nu = 3; nu <= 2 has no covariance. */
    ncm_stats_dist_kernel_st_set_nu (sdk_st, 3.0);
    h3 = ncm_stats_dist_kernel_get_rot_bandwidth (test->kernel, n);
    ncm_stats_dist_kernel_st_set_nu (sdk_st, 2.5);
    g_assert_cmpfloat (ncm_stats_dist_kernel_get_rot_bandwidth (test->kernel, n), ==, h3);
    ncm_stats_dist_kernel_st_set_nu (sdk_st, 2.0);
    g_assert_cmpfloat (ncm_stats_dist_kernel_get_var_factor (test->kernel), ==, GSL_POSINF);
    ncm_stats_dist_kernel_st_set_nu (sdk_st, test->nu);
  }
}

/* A random upper-triangular Cholesky factor, Sigma = U^T U. */
static NcmMatrix *
_test_random_cov_decomp (guint dim)
{
  NcmMatrix *U = ncm_matrix_new (dim, dim);
  guint i, j;

  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < dim; j++)
    {
      if (j < i)
        ncm_matrix_set (U, i, j, 0.0);
      else if (j == i)
        ncm_matrix_set (U, i, j, g_test_rand_double_range (0.5, 2.0));
      else
        ncm_matrix_set (U, i, j, g_test_rand_double_range (-1.0, 1.0));
    }
  }

  return U;
}

static void
test_ncm_stats_dist_kernel_norm (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  NcmMatrix *U        = _test_random_cov_decomp (test->dim);
  const guint ntests  = 100 * g_test_rand_int_range (1, 5);
  const gdouble I0    = _test_radial (test, TEST_RADIAL_NORM);
  gdouble lndet_Sigma = 0.0;
  guint i;

  for (i = 0; i < test->dim; i++)
    lndet_Sigma += 2.0 * log (ncm_matrix_get (U, i, i));

  /* int Kbar((x - mu)^T Sigma^-1 (x - mu)) dx = sqrt (det Sigma) I0 = u(Sigma) */
  ncm_assert_cmpdouble_e (ncm_stats_dist_kernel_get_lnnorm (test->kernel, U), ==, 0.5 * lndet_Sigma + log (I0), 1.0e-12, 1.0e-12);

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
      case NCM_STATS_DIST_KERNEL_TYPE_GAUSS:

        for (i = 0; i < ntests; i++)
        {
          gdouble eval_test = exp (-0.5 * ncm_vector_get (chi2_vec, i));

          ncm_assert_cmpdouble_e (ncm_vector_get (kernel_vec, i), ==, ncm_stats_dist_kernel_eval_unnorm (test->kernel, ncm_vector_get (chi2_vec, i)), 1.0e-15, 0.0);
          ncm_assert_cmpdouble_e (eval_test, ==, ncm_stats_dist_kernel_eval_unnorm (test->kernel, ncm_vector_get (chi2_vec, i)), 1.0e-15, 0.0);
        }

        break;
      case NCM_STATS_DIST_KERNEL_TYPE_ST3:

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

  ncm_matrix_free (U);
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
  NcmVector *lnc0        = ncm_vector_new (n);
  NcmVector *lnc1        = ncm_vector_new (n);
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


  /* Arm 0 gives every kernel its own normalization, arm 1 the one shared constant. Both
   * reach the kernel sum as log (w_i) minus that normalization, so the reference forms
   * the same combination in the same order and the comparison stays exact to 1e-15. */
  for (i = 0; i < n; i++)
  {
    const gdouble lnw_i = log (ncm_vector_fast_get (weights, i));

    ncm_vector_set (lnc0, i, lnw_i - ncm_vector_fast_get (lnnorms_vec, i));
    ncm_vector_set (lnc1, i, lnw_i - lnnorm);
  }

  for (i = 0; i < n; i++)
  {
    const gdouble chi2_i = ncm_vector_fast_get (chi2, i);
    const gdouble lnc_i0 = ncm_vector_get (lnc0, i);
    const gdouble lnc_i1 = ncm_vector_get (lnc1, i);

    switch (test->kernel_type)
    {
      case NCM_STATS_DIST_KERNEL_TYPE_GAUSS:
        lnt_i0 = -0.5 * chi2_i + lnc_i0;
        lnt_i1 = -0.5 * chi2_i + lnc_i1;
        break;
      case NCM_STATS_DIST_KERNEL_TYPE_ST3:
        lnt_i0 = kappa * log1p (chi2_i / test->nu) + lnc_i0;
        lnt_i1 = kappa * log1p (chi2_i / test->nu) + lnc_i1;
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

  ncm_stats_dist_kernel_eval_gamma_lambda (test->kernel, chi2, lnc0, lnK, &gamma0, &lambda0);
  ncm_stats_dist_kernel_eval_gamma_lambda (test->kernel, chi2, lnc1, lnK, &gamma1, &lambda1);

  for (i = 0; i < i_max0; i++)
  {
    lambda_test0 += exp (g_array_index (t_array0, gdouble, i) - lnt_max0);
    ncm_assert_cmpdouble_e (g_array_index (t_array0, gdouble, i), <, gamma0, 1.0e-15, 0.0);
  }

  for (i = 0; i < i_max1; i++)
  {
    lambda_test1 += exp (g_array_index (t_array1, gdouble, i) - lnt_max1);
    ncm_assert_cmpdouble_e (g_array_index (t_array1, gdouble, i), <, gamma1, 1.0e-15, 0.0);
  }

  for (i = i_max0 + 1; i < n; i++)
  {
    lambda_test0 += exp (g_array_index (t_array0, gdouble, i) - lnt_max0);
    ncm_assert_cmpdouble_e (g_array_index (t_array0, gdouble, i), <, gamma0, 1.0e-15, 0.0);
  }

  for (i = i_max1 + 1; i < n; i++)
  {
    lambda_test1 += exp (g_array_index (t_array1, gdouble, i) - lnt_max1);
    ncm_assert_cmpdouble_e (g_array_index (t_array1, gdouble, i), <, gamma1, 1.0e-15, 0.0);
  }

  ncm_assert_cmpdouble_e (lnt_max0, ==, gamma0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (lnt_max1, ==, gamma1, 1.0e-15, 0.0);

  ncm_assert_cmpdouble_e (lambda_test0, ==, lambda0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (lambda_test1, ==, lambda1, 1.0e-15, 0.0);

  /* The documented identity: e^gamma (1 + lambda) is the plain sum of the terms. */
  {
    gdouble sum0 = 0.0;
    gdouble sum1 = 0.0;

    for (i = 0; i < n; i++)
    {
      sum0 += exp (g_array_index (t_array0, gdouble, i));
      sum1 += exp (g_array_index (t_array1, gdouble, i));
    }

    ncm_assert_cmpdouble_e (exp (gamma0) * (1.0 + lambda0), ==, sum0, 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (exp (gamma1) * (1.0 + lambda1), ==, sum1, 1.0e-13, 0.0);
  }

  ncm_vector_free (weights);
  ncm_vector_free (chi2);
  ncm_vector_free (lnnorms_vec);
  ncm_vector_free (lnK);
  ncm_vector_free (lnc0);
  ncm_vector_free (lnc1);
  g_clear_pointer (&t_array0, g_array_unref);
  g_clear_pointer (&t_array1, g_array_unref);
}

static void
test_ncm_stats_dist_kernel_sample (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  const guint d         = test->dim;
  const guint nsamples  = 20000;
  const gdouble probs[] = {0.05, 0.25, 0.5, 0.75, 0.95};
  const guint nprobs    = G_N_ELEMENTS (probs);
  const gdouble href    = g_test_rand_double_range (0.1, 10.0);
  const gdouble kappa   = ncm_stats_dist_kernel_get_var_factor (test->kernel);
  NcmRNG *rng           = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  NcmMatrix *U          = _test_random_cov_decomp (d);
  NcmVector *mu         = ncm_vector_new (d);
  NcmVector *y          = ncm_vector_new (d);
  gsl_matrix *Sigma     = gsl_matrix_alloc (d, d);
  gsl_vector *v         = gsl_vector_alloc (d);
  gsl_vector *w         = gsl_vector_alloc (d);
  gdouble *mean         = g_new0 (gdouble, d);
  guint *count          = g_new0 (guint, nprobs);
  guint i, j, k;

  for (i = 0; i < d; i++)
    ncm_vector_set (mu, i, g_test_rand_double_range (-10.0, 10.0));

  /* Sigma = U^T U built and factored by GSL, independently of NcmMatrix. */
  for (i = 0; i < d; i++)
  {
    for (j = 0; j < d; j++)
    {
      gdouble S_ij = 0.0;

      for (k = 0; k < d; k++)
        S_ij += ncm_matrix_get (U, k, i) * ncm_matrix_get (U, k, j);

      gsl_matrix_set (Sigma, i, j, S_ij);
    }
  }

  gsl_linalg_cholesky_decomp1 (Sigma);

  for (i = 0; i < nsamples; i++)
  {
    gdouble chi2 = 0.0;

    ncm_stats_dist_kernel_sample (test->kernel, U, href, mu, y, rng);

    for (j = 0; j < d; j++)
    {
      const gdouble v_j = ncm_vector_get (y, j) - ncm_vector_get (mu, j);

      gsl_vector_set (v, j, v_j);
      mean[j] += v_j;
    }

    gsl_linalg_cholesky_solve (Sigma, v, w);

    for (j = 0; j < d; j++)
      chi2 += gsl_vector_get (v, j) * gsl_vector_get (w, j);

    chi2 /= href * href;

    /* chi2 follows chi-squared(d) for the Gaussian kernel, d F(d, nu) for Student-t */
    for (k = 0; k < nprobs; k++)
    {
      const gdouble q = (test->kernel_type == NCM_STATS_DIST_KERNEL_TYPE_GAUSS) ?
                        gsl_cdf_chisq_Pinv (probs[k], d) :
                        d *gsl_cdf_fdist_Pinv (probs[k], d, test->nu);

      count[k] += (chi2 <= q);
    }
  }

  for (k = 0; k < nprobs; k++)
  {
    const gdouble sd = sqrt (probs[k] * (1.0 - probs[k]) / nsamples);

    g_assert_cmpfloat (fabs (count[k] / (gdouble) nsamples - probs[k]), <, 5.0 * sd);
  }

  /* Each coordinate has variance kappa h^2 Sigma_jj about mu. */
  for (j = 0; j < d; j++)
  {
    gdouble Sigma_jj = 0.0;

    for (k = 0; k < d; k++)
      Sigma_jj += gsl_pow_2 (ncm_matrix_get (U, k, j));

    g_assert_cmpfloat (fabs (mean[j] / nsamples), <, 5.0 * href * sqrt (kappa * Sigma_jj / nsamples));
  }

  g_free (mean);
  g_free (count);
  gsl_matrix_free (Sigma);
  gsl_vector_free (v);
  gsl_vector_free (w);
  ncm_vector_free (mu);
  ncm_vector_free (y);
  ncm_matrix_free (U);
  ncm_rng_free (rng);
}

static void
test_ncm_stats_dist_kernel_free (TestNcmStatsDistKernel *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_dist_kernel_free, test->kernel);
}

static void
test_ncm_stats_dist_kernel_st_nu_default (void)
{
  NcmStatsDistKernelST *sdk_st = g_object_new (NCM_TYPE_STATS_DIST_KERNEL_ST, "dimension", 2, NULL);

  g_assert_cmpfloat (ncm_stats_dist_kernel_st_get_nu (sdk_st), ==, 3.0);
  ncm_assert_cmpdouble_e (ncm_stats_dist_kernel_get_var_factor (NCM_STATS_DIST_KERNEL (sdk_st)), ==, 3.0, 1.0e-15, 0.0);

  ncm_stats_dist_kernel_st_free (sdk_st);
}

static const gchar *unimplemented[] = {
  "get_rot_bandwidth", "get_var_factor", "get_lnnorm", "eval_unnorm",
  "eval_unnorm_vec", "eval_gamma_lambda", "sample"
};

static void
test_ncm_stats_dist_kernel_traps (void)
{
  guint k;

  for (k = 0; k < G_N_ELEMENTS (unimplemented); k++)
  {
    gchar *pattern = g_strdup_printf ("*method %s not implemented by TestNcmStatsDistKernelBare*", unimplemented[k]);

    g_setenv ("TEST_NCM_STATS_DIST_KERNEL_METHOD", unimplemented[k], TRUE);
    g_test_trap_subprocess ("/ncm/stats/dist/kernel/unimplemented/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (pattern);
    g_free (pattern);
  }

  g_unsetenv ("TEST_NCM_STATS_DIST_KERNEL_METHOD");
}

static void
test_ncm_stats_dist_kernel_unimplemented_subprocess (void)
{
  NcmStatsDistKernel *sdk = g_object_new (TEST_TYPE_NCM_STATS_DIST_KERNEL_BARE, "dimension", 2, NULL);
  const gchar *which      = g_getenv ("TEST_NCM_STATS_DIST_KERNEL_METHOD");
  gdouble gamma, lambda;

  g_assert_cmpuint (ncm_stats_dist_kernel_get_dim (sdk), ==, 2);

  if (g_str_equal (which, "get_rot_bandwidth"))
    ncm_stats_dist_kernel_get_rot_bandwidth (sdk, 10.0);
  else if (g_str_equal (which, "get_var_factor"))
    ncm_stats_dist_kernel_get_var_factor (sdk);
  else if (g_str_equal (which, "get_lnnorm"))
    ncm_stats_dist_kernel_get_lnnorm (sdk, NULL);
  else if (g_str_equal (which, "eval_unnorm"))
    ncm_stats_dist_kernel_eval_unnorm (sdk, 1.0);
  else if (g_str_equal (which, "eval_unnorm_vec"))
    ncm_stats_dist_kernel_eval_unnorm_vec (sdk, NULL, NULL);
  else if (g_str_equal (which, "eval_gamma_lambda"))
    ncm_stats_dist_kernel_eval_gamma_lambda (sdk, NULL, NULL, NULL, &gamma, &lambda);
  else
    ncm_stats_dist_kernel_sample (sdk, NULL, 1.0, NULL, NULL, NULL);
}


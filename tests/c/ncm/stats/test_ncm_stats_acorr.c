/***************************************************************************
 *            test_ncm_stats_acorr.c
 *
 *  Tue Sep 16 09:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

typedef struct _TestNcmStatsAcorr
{
  NcmStatsAcorr *acorr;
} TestNcmStatsAcorr;

void test_ncm_stats_acorr_new (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_free (TestNcmStatsAcorr *test, gconstpointer pdata);

void test_ncm_stats_acorr_basic (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_properties (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_feed_equivalence (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_levels (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_acf (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_mean_error (TestNcmStatsAcorr *test, gconstpointer pdata);

void test_ncm_stats_acorr_est_ar1_analytic (void);
void test_ncm_stats_acorr_est_ar2_analytic (void);
void test_ncm_stats_acorr_est_white_analytic (void);
void test_ncm_stats_acorr_est_degenerate (void);
void test_ncm_stats_acorr_est_not_positive_definite (void);
void test_ncm_stats_acorr_est_ar_criteria (void);
void test_ncm_stats_acorr_ar_fit_coefficients (void);

void test_ncm_stats_acorr_exact_vs_fft (void);
void test_ncm_stats_acorr_acov_fft_all_lags (void);
void test_ncm_stats_acorr_recenter (void);

void test_ncm_stats_acorr_tau_ar1_sampled (void);
void test_ncm_stats_acorr_tau_cascade (void);

void test_ncm_stats_acorr_diag_white (void);
void test_ncm_stats_acorr_diag_ramp (void);
void test_ncm_stats_acorr_diag_drift (void);
void test_ncm_stats_acorr_diag_truncated (void);
void test_ncm_stats_acorr_diag_disagreement (void);
void test_ncm_stats_acorr_diag_zero_variance (void);
void test_ncm_stats_acorr_diag_variance_shift (void);
void test_ncm_stats_acorr_diag_empty (void);
void test_ncm_stats_acorr_diag_to_string (void);

void test_ncm_stats_acorr_traps (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_invalid_var (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_invalid_level (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_invalid_method (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_invalid_ar_crit (TestNcmStatsAcorr *test, gconstpointer pdata);
void test_ncm_stats_acorr_invalid_factor (TestNcmStatsAcorr *test, gconstpointer pdata);

#define TEST_ACORR_PREC (1.0e-9)

/* An AR(1) process has C_k = rho^k and an integrated autocorrelation time of
 * (1 + rho) / (1 - rho). Both are exact, so the estimators are checked against a number
 * and not against a sample. */
static NcmVector *
_test_acorr_ar1_acov (gdouble rho, guint nlag)
{
  NcmVector *acov = ncm_vector_new (nlag + 1);
  guint k;

  for (k = 0; k <= nlag; k++)
    ncm_vector_set (acov, k, pow (rho, k));

  return acov;
}

static gdouble
_test_acorr_ar1_tau (gdouble rho)
{
  return (1.0 + rho) / (1.0 - rho);
}

/* AR(2): x_t = phi1 x_{t-1} + phi2 x_{t-2} + e_t with unit innovation variance. Its
 * autocovariances follow the Yule-Walker recursion from
 * C_0 = (1 - phi2) / [(1 + phi2) (1 - phi1 - phi2) (1 + phi1 - phi2)] and C_1 / C_0 =
 * phi1 / (1 - phi2), and the sum of all of them is S(0) = (1 - phi1 - phi2)^{-2}, so
 * tau = S(0) / C_0 is closed form as well. */
static NcmVector *
_test_acorr_ar2_acov (gdouble phi1, gdouble phi2, guint nlag, gdouble *tau)
{
  NcmVector *acov  = ncm_vector_new (nlag + 1);
  const gdouble r1 = phi1 / (1.0 - phi2);
  const gdouble c0 = (1.0 - phi2) / ((1.0 + phi2) * (1.0 - phi1 - phi2) * (1.0 + phi1 - phi2));
  guint k;

  ncm_vector_set (acov, 0, c0);
  ncm_vector_set (acov, 1, c0 * r1);

  for (k = 2; k <= nlag; k++)
    ncm_vector_set (acov, k, phi1 * ncm_vector_get (acov, k - 1) + phi2 * ncm_vector_get (acov, k - 2));

  tau[0] = 1.0 / (gsl_pow_2 (1.0 - phi1 - phi2) * c0);

  return acov;
}

static void
_test_acorr_fill_ar1 (NcmVector *x, gdouble rho, NcmRNG *rng)
{
  const guint n     = ncm_vector_len (x);
  const gdouble sig = sqrt (1.0 - rho * rho);
  guint i;

  ncm_vector_set (x, 0, ncm_rng_ugaussian_gen (rng));

  for (i = 1; i < n; i++)
    ncm_vector_set (x, i, rho * ncm_vector_get (x, i - 1) + sig * ncm_rng_ugaussian_gen (rng));
}

static NcmStatsAcorr *
_test_acorr_from_vector (NcmVector *x, guint max_lag, guint max_levels, NcmStatsAcorrMethod method)
{
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, max_lag, max_levels, method);

  ncm_stats_acorr_set_series (acorr, 0, x);

  return acorr;
}

void
test_ncm_stats_acorr_new (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  test->acorr = ncm_stats_acorr_new (3);

  g_assert_true (NCM_IS_STATS_ACORR (test->acorr));
}

void
test_ncm_stats_acorr_free (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_stats_acorr_free, test->acorr);
}

/* Construction, defaults and the state of an accumulator that has been fed and reset. */
void
test_ncm_stats_acorr_basic (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NcmStatsAcorr *acorr = test->acorr;
  NcmVector *x         = ncm_vector_new (3);
  guint i, p;

  g_assert_cmpuint (ncm_stats_acorr_len (acorr), ==, 3);
  g_assert_cmpuint (ncm_stats_acorr_get_max_lag (acorr), ==, NCM_STATS_ACORR_DEFAULT_MAX_LAG);
  g_assert_cmpuint (ncm_stats_acorr_get_max_levels (acorr), ==, NCM_STATS_ACORR_DEFAULT_MAX_LEVELS);
  g_assert_cmpuint (ncm_stats_acorr_get_method (acorr), ==, NCM_STATS_ACORR_METHOD_MAX);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_reliability_factor (acorr), ==,
                          NCM_STATS_ACORR_DEFAULT_RELIABILITY_FACTOR, TEST_ACORR_PREC, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_drift_threshold (acorr), ==,
                          NCM_STATS_ACORR_DEFAULT_DRIFT_THRESHOLD, TEST_ACORR_PREC, 0.0);

  for (p = 0; p < 3; p++)
    g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, p), ==, 0);

  for (i = 0; i < 100; i++)
  {
    for (p = 0; p < 3; p++)
      ncm_vector_set (x, p, 1.0 * i + 10.0 * p);

    ncm_stats_acorr_update (acorr, x);
  }

  for (p = 0; p < 3; p++)
  {
    g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, p), ==, 100);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (acorr, p), ==, 49.5 + 10.0 * p, 1.0e-12, 0.0);
  }

  ncm_stats_acorr_reset (acorr);

  for (p = 0; p < 3; p++)
  {
    g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, p), ==, 0);
    g_assert_cmpuint (ncm_stats_acorr_nlevels (acorr, p), ==, 1);
  }

  /* Every method is accepted and reported back, and switching does not touch the data. */
  {
    const NcmStatsAcorrMethod methods[] = {
      NCM_STATS_ACORR_METHOD_AR,
      NCM_STATS_ACORR_METHOD_GEYER,
      NCM_STATS_ACORR_METHOD_SOKAL,
      NCM_STATS_ACORR_METHOD_MAX
    };

    for (i = 0; i < 100; i++)
    {
      for (p = 0; p < 3; p++)
        ncm_vector_set (x, p, sin (0.1 * i) + 0.01 * p);

      ncm_stats_acorr_update (acorr, x);
    }

    for (i = 0; i < G_N_ELEMENTS (methods); i++)
    {
      ncm_stats_acorr_set_method (acorr, methods[i]);
      g_assert_cmpuint (ncm_stats_acorr_get_method (acorr), ==, methods[i]);
      g_assert_true (gsl_finite (ncm_stats_acorr_get_tau (acorr, 0)));
      g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, 0), ==, 100);
    }
  }

  ncm_vector_free (x);
}

/* The properties are the configuration, and they round-trip through GObject. */
void
test_ncm_stats_acorr_properties (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (2, 64, 5, NCM_STATS_ACORR_METHOD_GEYER);
  guint len, max_lag, max_levels;
  NcmStatsAcorrMethod method;
  gdouble factor, threshold;

  g_object_get (acorr,
                "len", &len,
                "max-lag", &max_lag,
                "max-levels", &max_levels,
                "method", &method,
                "reliability-factor", &factor,
                "drift-threshold", &threshold,
                NULL);

  g_assert_cmpuint (len, ==, 2);
  g_assert_cmpuint (max_lag, ==, 64);
  g_assert_cmpuint (max_levels, ==, 5);
  g_assert_cmpuint (method, ==, NCM_STATS_ACORR_METHOD_GEYER);

  g_object_set (acorr,
                "method", NCM_STATS_ACORR_METHOD_SOKAL,
                "reliability-factor", 10.0,
                "drift-threshold", 2.5,
                NULL);

  g_assert_cmpuint (ncm_stats_acorr_get_method (acorr), ==, NCM_STATS_ACORR_METHOD_SOKAL);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_reliability_factor (acorr), ==, 10.0, TEST_ACORR_PREC, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_drift_threshold (acorr), ==, 2.5, TEST_ACORR_PREC, 0.0);

  /* The order criterion is a property like the others, and changing it leaves the data. */
  g_object_set (acorr, "ar-criterion", NCM_STATS_ACORR_AR_CRIT_FPE, NULL);
  g_assert_cmpuint (ncm_stats_acorr_get_ar_criterion (acorr), ==, NCM_STATS_ACORR_AR_CRIT_FPE);
  {
    NcmStatsAcorrARCrit crit;

    g_object_get (acorr, "ar-criterion", &crit, NULL);
    g_assert_cmpuint (crit, ==, NCM_STATS_ACORR_AR_CRIT_FPE);
  }
  ncm_stats_acorr_set_ar_criterion (acorr, NCM_STATS_ACORR_AR_CRIT_AICC);
  g_assert_cmpuint (ncm_stats_acorr_get_ar_criterion (acorr), ==, NCM_STATS_ACORR_AR_CRIT_AICC);

  ncm_stats_acorr_set_reliability_factor (acorr, 25.0);
  ncm_stats_acorr_set_drift_threshold (acorr, 4.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_reliability_factor (acorr), ==, 25.0, TEST_ACORR_PREC, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_drift_threshold (acorr), ==, 4.0, TEST_ACORR_PREC, 0.0);

  ncm_stats_acorr_free (acorr);
}

/* The four ways of feeding the same numbers give the same accumulator. */
void
test_ncm_stats_acorr_feed_equivalence (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, 987654321);
  const guint n           = 2000;
  NcmMatrix *m            = ncm_matrix_new (n, 2);
  NcmVector *c0           = ncm_vector_new (n);
  NcmVector *c1           = ncm_vector_new (n);
  NcmStatsAcorr *a_update = ncm_stats_acorr_new (2);
  NcmStatsAcorr *a_matrix = ncm_stats_acorr_new (2);
  NcmStatsAcorr *a_series = ncm_stats_acorr_new (2);
  guint i, p;

  for (i = 0; i < n; i++)
  {
    ncm_matrix_set (m, i, 0, ncm_rng_ugaussian_gen (rng));
    ncm_matrix_set (m, i, 1, 3.0 + 2.0 * ncm_rng_ugaussian_gen (rng));
    ncm_vector_set (c0, i, ncm_matrix_get (m, i, 0));
    ncm_vector_set (c1, i, ncm_matrix_get (m, i, 1));
  }

  for (i = 0; i < n; i++)
  {
    NcmVector *row = ncm_matrix_get_row (m, i);

    ncm_stats_acorr_update (a_update, row);
    ncm_vector_free (row);
  }

  ncm_stats_acorr_set_series_matrix (a_matrix, m);
  ncm_stats_acorr_set_series (a_series, 0, c0);
  ncm_stats_acorr_set_series (a_series, 1, c1);

  for (p = 0; p < 2; p++)
  {
    const gdouble tau = ncm_stats_acorr_get_tau (a_update, p);

    g_assert_cmpuint (ncm_stats_acorr_nitens (a_matrix, p), ==, n);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (a_matrix, p), ==, ncm_stats_acorr_get_mean (a_update, p), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var (a_matrix, p), ==, ncm_stats_acorr_get_var (a_update, p), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (a_matrix, p), ==, tau, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (a_series, p), ==, tau, 1.0e-12, 0.0);
  }

  /* update_var feeds one series alone. */
  {
    NcmStatsAcorr *a_one = ncm_stats_acorr_new (1);

    for (i = 0; i < n; i++)
      ncm_stats_acorr_update_var (a_one, 0, ncm_vector_get (c0, i));

    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (a_one, 0), ==, ncm_stats_acorr_get_tau (a_update, 0), 1.0e-12, 0.0);
    ncm_stats_acorr_free (a_one);
  }

  ncm_stats_acorr_free (a_update);
  ncm_stats_acorr_free (a_matrix);
  ncm_stats_acorr_free (a_series);
  ncm_matrix_free (m);
  ncm_vector_free (c0);
  ncm_vector_free (c1);
  ncm_rng_free (rng);
}

/* Level j holds the means of 2^j consecutive samples: their count, and the mean they
 * preserve. */
void
test_ncm_stats_acorr_levels (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, 16, 4, NCM_STATS_ACORR_METHOD_AR);
  const guint n        = 1000;
  guint i, j;

  for (i = 0; i < n; i++)
    ncm_stats_acorr_update_var (acorr, 0, 1.0 * i);

  g_assert_cmpuint (ncm_stats_acorr_nlevels (acorr, 0), ==, 4);

  for (j = 0; j < ncm_stats_acorr_nlevels (acorr, 0); j++)
    g_assert_cmpuint (ncm_stats_acorr_level_nitens (acorr, 0, j), ==, n >> j);

  /* Level 1 of 0, 1, 2, ... is 0.5, 2.5, 4.5, ..., so its autocovariance at lag 0 is the
   * variance of that arithmetic sequence. */
  {
    NcmVector *acov1  = ncm_stats_acorr_get_acov (acorr, 0, 1);
    const guint n1    = n / 2;
    const gdouble var = (n1 * n1 - 1.0) / 3.0;

    ncm_assert_cmpdouble_e (ncm_vector_get (acov1, 0), ==, var, 1.0e-12, 0.0);
    ncm_vector_free (acov1);
  }

  ncm_stats_acorr_free (acorr);
}

/* The autocorrelation function is the autocovariance scaled by its first element. */
void
test_ncm_stats_acorr_acf (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, 13579);
  NcmVector *x         = ncm_vector_new (500);
  NcmStatsAcorr *acorr = NULL;
  NcmVector *acov, *acf;
  guint k;

  _test_acorr_fill_ar1 (x, 0.5, rng);
  acorr = _test_acorr_from_vector (x, 64, 4, NCM_STATS_ACORR_METHOD_MAX);

  acov = ncm_stats_acorr_get_acov (acorr, 0, 0);
  acf  = ncm_stats_acorr_get_acf (acorr, 0, 0);

  g_assert_cmpuint (ncm_vector_len (acf), ==, ncm_vector_len (acov));
  ncm_assert_cmpdouble_e (ncm_vector_get (acf, 0), ==, 1.0, 1.0e-14, 0.0);

  for (k = 1; k < ncm_vector_len (acf); k++)
    ncm_assert_cmpdouble_e (ncm_vector_get (acf, k), ==,
                            ncm_vector_get (acov, k) / ncm_vector_get (acov, 0), 1.0e-13, 0.0);

  /* A series with no variance has an autocorrelation function left unscaled. */
  {
    NcmStatsAcorr *flat = ncm_stats_acorr_new_full (1, 16, 2, NCM_STATS_ACORR_METHOD_MAX);
    NcmVector *flat_acf;

    for (k = 0; k < 50; k++)
      ncm_stats_acorr_update_var (flat, 0, 2.0);

    flat_acf = ncm_stats_acorr_get_acf (flat, 0, 0);
    ncm_assert_cmpdouble_e (ncm_vector_get (flat_acf, 0), ==, 0.0, 1.0e-14, 1.0e-14);

    ncm_vector_free (flat_acf);
    ncm_stats_acorr_free (flat);
  }

  ncm_vector_free (acov);
  ncm_vector_free (acf);
  ncm_vector_free (x);
  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* The error of the mean is built from the long-run variance and nothing else. */
void
test_ncm_stats_acorr_mean_error (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  NcmRNG *rng          = ncm_rng_pool_get ("test_ncm_stats_acorr");
  const guint n        = 20000;
  NcmVector *x         = ncm_vector_new (n);
  NcmStatsAcorr *acorr = NULL;
  gdouble tau, var, spec0;
  guint i;

  for (i = 0; i < n; i++)
    ncm_vector_set (x, i, ncm_rng_ugaussian_gen (rng));

  acorr = _test_acorr_from_vector (x, 512, 24, NCM_STATS_ACORR_METHOD_MAX);

  tau   = ncm_stats_acorr_get_tau (acorr, 0);
  var   = ncm_stats_acorr_get_var (acorr, 0);
  spec0 = ncm_stats_acorr_get_spec0 (acorr, 0);

  ncm_assert_cmpdouble_e (spec0, ==, tau * var, 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var_mean (acorr, 0), ==, spec0 / n, 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_sd_mean (acorr, 0), ==, sqrt (spec0 / n), 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_ess (acorr, 0), ==, n / tau, 1.0e-13, 0.0);

  /* Independent samples: the error of the mean is the usual one to within the scatter of
   * the estimate of tau. */
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_sd_mean (acorr, 0), ==, sqrt (var / n), 0.15, 0.0);

  ncm_vector_free (x);
  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* On the autocovariances of an AR(1) process all three estimators are exact. */
void
test_ncm_stats_acorr_est_ar1_analytic (void)
{
  const gdouble rhos[] = { 0.0, 0.3, 0.5, 0.9, 0.99, -0.5 };
  const guint64 n      = 10000000;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (rhos); i++)
  {
    const gdouble rho   = rhos[i];
    const gdouble exact = _test_acorr_ar1_tau (rho);
    NcmVector *acov     = _test_acorr_ar1_acov (rho, 2000);
    guint order = 0, window = 0;
    gdouble tau_ar, tau_geyer, tau_sokal;

    tau_ar    = ncm_stats_acorr_tau_ar (acov, n, NCM_STATS_ACORR_AR_CRIT_AICC, &order);
    tau_geyer = ncm_stats_acorr_tau_geyer (acov, &window);
    tau_sokal = ncm_stats_acorr_tau_sokal (acov, NCM_STATS_ACORR_SOKAL_C, NULL);

    ncm_assert_cmpdouble_e (tau_ar,    ==, exact, 1.0e-5, 0.0);
    ncm_assert_cmpdouble_e (tau_geyer, ==, exact, 1.0e-5, 0.0);

    /* Sokal's window rule assumes a non-negative autocorrelation function and reports 1
     * when it is not, which is the documented behaviour of that estimator. */
    if (rho >= 0.0)
      ncm_assert_cmpdouble_e (tau_sokal, ==, exact, 1.0e-3, 0.0);
    else
      ncm_assert_cmpdouble_e (tau_sokal, ==, 1.0, 1.0e-14, 0.0);

    /* AICc recovers the order of the process: one, or none when there is no correlation. */
    g_assert_cmpuint (order, ==, (rho == 0.0) ? 0 : 1);
    g_assert_cmpuint (window, >, 0);

    ncm_vector_free (acov);
  }
}

/* Same for an AR(2), where the order to be recovered is two. */
void
test_ncm_stats_acorr_est_ar2_analytic (void)
{
  gdouble exact           = 0.0;
  NcmVector *acov         = _test_acorr_ar2_acov (0.6, 0.2, 2000, &exact);
  guint order             = 0;
  const gdouble tau_ar    = ncm_stats_acorr_tau_ar (acov, 10000000, NCM_STATS_ACORR_AR_CRIT_AICC, &order);
  const gdouble tau_geyer = ncm_stats_acorr_tau_geyer (acov, NULL);

  ncm_assert_cmpdouble_e (tau_ar,    ==, exact, 1.0e-5, 0.0);
  ncm_assert_cmpdouble_e (tau_geyer, ==, exact, 1.0e-5, 0.0);
  g_assert_cmpuint (order, ==, 2);

  ncm_vector_free (acov);
}

/* Uncorrelated samples: tau is one, the auto-regressive order is zero. */
void
test_ncm_stats_acorr_est_white_analytic (void)
{
  NcmVector *acov = ncm_vector_new (100);
  guint order = 0, window = 0;

  ncm_vector_set_all (acov, 0.0);
  ncm_vector_set (acov, 0, 2.5);

  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_ar (acov, 100000, NCM_STATS_ACORR_AR_CRIT_AICC, &order), ==, 1.0, 1.0e-6, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_geyer (acov, &window), ==, 1.0, 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_sokal (acov, NCM_STATS_ACORR_SOKAL_C, NULL), ==, 1.0, 1.0e-12, 0.0);
  g_assert_cmpuint (order, ==, 0);

  ncm_vector_free (acov);
}

/* A sequence with no variance, and one whose first pair sum is not positive, are the two
 * inputs on which an integrated time is not defined. Both give one. */
void
test_ncm_stats_acorr_est_degenerate (void)
{
  NcmVector *zero = ncm_vector_new (10);
  NcmVector *neg  = ncm_vector_new (10);
  guint order = 1, window = 1;

  ncm_vector_set_all (zero, 0.0);

  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_ar (zero, 1000, NCM_STATS_ACORR_AR_CRIT_AICC, &order), ==, 1.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_geyer (zero, &window), ==, 1.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_sokal (zero, NCM_STATS_ACORR_SOKAL_C, NULL), ==, 1.0, 1.0e-14, 0.0);
  g_assert_cmpuint (order, ==, 0);
  g_assert_cmpuint (window, ==, 0);

  /* C_0 + C_1 <= 0: the initial sequence is empty. */
  ncm_vector_set_all (neg, 0.0);
  ncm_vector_set (neg, 0, 1.0);
  ncm_vector_set (neg, 1, -1.5);

  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_geyer (neg, &window), ==, 1.0, 1.0e-14, 0.0);
  g_assert_cmpuint (window, ==, 0);

  /* A short series is below what an auto-regressive fit is attempted on. */
  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_ar (neg, 4, NCM_STATS_ACORR_AR_CRIT_AICC, NULL), ==, 1.0, 1.0e-14, 0.0);

  /* Sokal never reaching its window sums the whole sequence. */
  {
    NcmVector *slow = _test_acorr_ar1_acov (0.99, 8);
    guint w         = 0;

    ncm_stats_acorr_tau_sokal (slow, NCM_STATS_ACORR_SOKAL_C, &w);
    g_assert_cmpuint (w, ==, 8);
    ncm_vector_free (slow);
  }

  ncm_vector_free (zero);
  ncm_vector_free (neg);
}

/* A sequence that is not a valid autocovariance stops the recursion instead of producing
 * a number from it. */
void
test_ncm_stats_acorr_est_not_positive_definite (void)
{
  NcmVector *acov = ncm_vector_new (10);
  guint order     = 5;

  ncm_vector_set_all (acov, 0.0);
  ncm_vector_set (acov, 0, 1.0);
  ncm_vector_set (acov, 1, 2.0);

  ncm_assert_cmpdouble_e (ncm_stats_acorr_tau_ar (acov, 10000, NCM_STATS_ACORR_AR_CRIT_AICC, &order), ==, 1.0, 1.0e-6, 0.0);
  g_assert_cmpuint (order, ==, 0);

  ncm_vector_free (acov);
}

/* Every order-selection rule on the autocovariances of a process whose order they should
 * agree on, plus the rule that selects nothing and takes the largest order offered. */
void
test_ncm_stats_acorr_est_ar_criteria (void)
{
  const NcmStatsAcorrARCrit crits[] = {
    NCM_STATS_ACORR_AR_CRIT_FPE,
    NCM_STATS_ACORR_AR_CRIT_AIC,
    NCM_STATS_ACORR_AR_CRIT_AICC
  };
  const gdouble rho   = 0.8;
  const gdouble exact = _test_acorr_ar1_tau (rho);
  NcmVector *acov     = _test_acorr_ar1_acov (rho, 2000);
  guint i;

  for (i = 0; i < G_N_ELEMENTS (crits); i++)
  {
    guint order       = 0;
    const gdouble tau = ncm_stats_acorr_tau_ar (acov, 10000000, crits[i], &order);

    /* An AR(1) is order one under every rule that selects, and exactly reproduced. */
    g_assert_cmpuint (order, ==, 1);
    ncm_assert_cmpdouble_e (tau, ==, exact, 1.0e-5, 0.0);
  }

  /* Selecting nothing takes the largest order the search offers, floor (10 log10 n), and
   * still reproduces tau: the extra coefficients of an exactly fitted AR(1) are zero. */
  {
    guint order       = 0;
    const gdouble tau = ncm_stats_acorr_tau_ar (acov, 10000000, NCM_STATS_ACORR_AR_CRIT_NONE, &order);

    g_assert_cmpuint (order, ==, (guint) floor (10.0 * log10 (10000000.0)));
    ncm_assert_cmpdouble_e (tau, ==, exact, 1.0e-5, 0.0);
  }

  ncm_vector_free (acov);
}

/* The fit itself: the coefficients of an AR(1) are (rho), and its reflection coefficients
 * are rho at lag one and zero after. */
void
test_ncm_stats_acorr_ar_fit_coefficients (void)
{
  const gdouble rho = 0.8;
  NcmVector *acov   = _test_acorr_ar1_acov (rho, 2000);
  NcmVector *phi = NULL, *pacf = NULL;
  gdouble ivar = 0.0;
  guint order  = 0;
  guint k;

  g_assert_true (ncm_stats_acorr_ar_fit (acov, 10000000, NCM_STATS_ACORR_AR_CRIT_AICC,
                                         &phi, &pacf, &ivar, &order));

  g_assert_cmpuint (order, ==, 1);
  g_assert_cmpuint (ncm_vector_len (phi), ==, 1);
  ncm_assert_cmpdouble_e (ncm_vector_get (phi, 0), ==, rho, 1.0e-12, 0.0);

  /* The innovation variance of an AR(1) with unit variance is 1 - rho^2. */
  ncm_assert_cmpdouble_e (ivar, ==, 1.0 - rho * rho, 1.0e-5, 0.0);

  /* The partial autocorrelations cut off after the order of the process. */
  g_assert_cmpuint (ncm_vector_len (pacf), >, 1);
  ncm_assert_cmpdouble_e (ncm_vector_get (pacf, 0), ==, rho, 1.0e-12, 0.0);

  for (k = 1; k < ncm_vector_len (pacf); k++)
    ncm_assert_cmpdouble_e (ncm_vector_get (pacf, k), ==, 0.0, 1.0e-12, 1.0e-10);

  ncm_vector_free (phi);
  ncm_vector_free (pacf);

  /* A sequence with no variance is fitted at order zero, and the outputs say so. */
  {
    NcmVector *zero = ncm_vector_new (10);

    ncm_vector_set_all (zero, 0.0);
    phi   = NULL;
    pacf  = NULL;
    order = 1;

    g_assert_false (ncm_stats_acorr_ar_fit (zero, 1000, NCM_STATS_ACORR_AR_CRIT_AICC,
                                            &phi, &pacf, &ivar, &order));
    g_assert_null (phi);
    g_assert_null (pacf);
    g_assert_cmpuint (order, ==, 0);

    ncm_vector_free (zero);
  }

  ncm_vector_free (acov);
}

/* The accumulated autocovariances and the ones computed in one pass over the whole
 * series by Fourier transform are the same numbers. */
void
test_ncm_stats_acorr_exact_vs_fft (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, 24680);
  const guint n        = 4096;
  NcmVector *x         = ncm_vector_new (n);
  NcmStatsAcorr *acorr = NULL;
  NcmVector *acov_acc, *acov_fft;
  guint k;

  _test_acorr_fill_ar1 (x, 0.7, rng);

  acorr    = _test_acorr_from_vector (x, 128, 8, NCM_STATS_ACORR_METHOD_MAX);
  acov_acc = ncm_stats_acorr_get_acov (acorr, 0, 0);
  acov_fft = ncm_stats_acorr_acov_fft (x, 128);

  g_assert_cmpuint (ncm_vector_len (acov_fft), ==, ncm_vector_len (acov_acc));

  for (k = 0; k < ncm_vector_len (acov_acc); k++)
    ncm_assert_cmpdouble_e (ncm_vector_get (acov_acc, k), ==, ncm_vector_get (acov_fft, k), 1.0e-10, 1.0e-13);

  ncm_vector_free (acov_acc);
  ncm_vector_free (acov_fft);
  ncm_vector_free (x);
  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* Asking for every lag gives one per sample less one. */
void
test_ncm_stats_acorr_acov_fft_all_lags (void)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, 111);
  const guint n = 512;
  NcmVector *x  = ncm_vector_new (n);
  NcmVector *acov;
  guint i;

  for (i = 0; i < n; i++)
    ncm_vector_set (x, i, ncm_rng_ugaussian_gen (rng));

  acov = ncm_stats_acorr_acov_fft (x, 0);
  g_assert_cmpuint (ncm_vector_len (acov), ==, n);
  g_assert_cmpfloat (ncm_vector_get (acov, 0), >, 0.0);
  ncm_vector_free (acov);

  /* A lag budget past the end of the series is the same as asking for all of them. */
  acov = ncm_stats_acorr_acov_fft (x, 10 * n);
  g_assert_cmpuint (ncm_vector_len (acov), ==, n);
  ncm_vector_free (acov);

  ncm_vector_free (x);
  ncm_rng_free (rng);
}

/* A large constant added to every sample cancels out: the accumulator moves its own
 * origin instead of subtracting two large numbers. */
void
test_ncm_stats_acorr_recenter (void)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 5150);
  const guint n     = 5000;
  NcmVector *x      = ncm_vector_new (n);
  NcmVector *x_off  = ncm_vector_new (n);
  const gdouble off = 1.0e8;
  NcmStatsAcorr *plain, *shifted;
  guint i;

  _test_acorr_fill_ar1 (x, 0.6, rng);

  for (i = 0; i < n; i++)
    ncm_vector_set (x_off, i, ncm_vector_get (x, i) + off);

  plain   = _test_acorr_from_vector (x,     256, 8, NCM_STATS_ACORR_METHOD_MAX);
  shifted = _test_acorr_from_vector (x_off, 256, 8, NCM_STATS_ACORR_METHOD_MAX);

  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (shifted, 0), ==,
                          ncm_stats_acorr_get_mean (plain, 0) + off, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var (shifted, 0), ==,
                          ncm_stats_acorr_get_var (plain, 0), 1.0e-9, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (shifted, 0), ==,
                          ncm_stats_acorr_get_tau (plain, 0), 1.0e-6, 0.0);

  /* A series walking away from where it started is what moves the origin, and it stays
   * exact while it does. */
  {
    NcmStatsAcorr *walking = ncm_stats_acorr_new_full (1, 32, 6, NCM_STATS_ACORR_METHOD_AR);
    gdouble sum = 0.0, sum2 = 0.0;

    for (i = 0; i < 2000; i++)
    {
      const gdouble v = 1.0 * i;

      ncm_stats_acorr_update_var (walking, 0, v);
      sum  += v;
      sum2 += v * v;
    }

    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (walking, 0), ==, sum / 2000.0, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var (walking, 0), ==,
                            sum2 / 2000.0 - gsl_pow_2 (sum / 2000.0), 1.0e-9, 0.0);
    ncm_stats_acorr_free (walking);
  }

  /* A chain entering from far outside its own scatter -- the first sample is the origin
   * the accumulator starts from, and the run then sits thousands of standard deviations
   * away from it. Moving the origin is exact, so the result is the one the same samples
   * give without the excursion. */
  {
    NcmStatsAcorr *entering = ncm_stats_acorr_new_full (1, 256, 8, NCM_STATS_ACORR_METHOD_MAX);
    NcmStatsAcorr *settled  = ncm_stats_acorr_new_full (1, 256, 8, NCM_STATS_ACORR_METHOD_MAX);

    ncm_stats_acorr_update_var (entering, 0, -1.0e4);

    for (i = 0; i < n; i++)
    {
      ncm_stats_acorr_update_var (entering, 0, ncm_vector_get (x, i));
      ncm_stats_acorr_update_var (settled,  0, ncm_vector_get (x, i));
    }

    /* The mean of the two differs by the one sample that is far away, and by nothing
     * else: the accumulated sums carry no loss from the excursion. */
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (entering, 0), ==,
                            (ncm_stats_acorr_get_mean (settled, 0) * n - 1.0e4) / (n + 1.0), 1.0e-12, 0.0);
    g_assert_true (gsl_finite (ncm_stats_acorr_get_tau (entering, 0)));
    g_assert_cmpfloat (ncm_stats_acorr_get_var (entering, 0), >, 0.0);

    ncm_stats_acorr_free (entering);
    ncm_stats_acorr_free (settled);
  }

  ncm_vector_free (x);
  ncm_vector_free (x_off);
  ncm_stats_acorr_free (plain);
  ncm_stats_acorr_free (shifted);
  ncm_rng_free (rng);
}

/* Sampled AR(1): the estimate has to land on the value the process is built with. */
void
test_ncm_stats_acorr_tau_ar1_sampled (void)
{
  NcmRNG *rng         = ncm_rng_pool_get ("test_ncm_stats_acorr");
  const gdouble rho   = 0.8;
  const gdouble exact = _test_acorr_ar1_tau (rho);
  const guint n       = 50000;
  NcmVector *x        = ncm_vector_new (n);
  NcmStatsAcorr *acorr;

  _test_acorr_fill_ar1 (x, rho, rng);
  acorr = _test_acorr_from_vector (x, 512, 24, NCM_STATS_ACORR_METHOD_MAX);

  /* The scatter of the estimate is of order sqrt (2 c tau / n); the bound is that with
   * room for the upward bias the conservative estimator carries by design. */
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, exact, 0.15, 0.0);
  g_assert_cmpuint (ncm_stats_acorr_get_level (acorr, 0), ==, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0), ==, NCM_STATS_ACORR_DIAG_OK);

  /* Each estimator alone, on the level the selection settled on. */
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_AR), ==, exact, 0.15, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_GEYER), ==, exact, 0.15, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_SOKAL), ==, exact, 0.15, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_MAX), ==, exact, 0.15, 0.0);

  g_assert_cmpuint (ncm_stats_acorr_get_ar_order (acorr, 0), >, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_window (acorr, 0), >, 0);

  /* The fit behind that order is available, and is the one the order came from. */
  {
    NcmVector *phi = NULL, *pacf = NULL;
    gdouble ivar = 0.0;
    guint order  = 0;

    g_assert_true (ncm_stats_acorr_get_ar_fit (acorr, 0, &phi, &pacf, &ivar, &order));
    g_assert_cmpuint (order, ==, ncm_stats_acorr_get_ar_order (acorr, 0));
    g_assert_cmpuint (ncm_vector_len (phi), ==, order);
    g_assert_cmpuint (ncm_vector_len (pacf), >=, order);
    g_assert_cmpfloat (ivar, >, 0.0);

    /* Coefficients of a stationary fit, and reflection coefficients inside the unit disc. */
    ncm_assert_cmpdouble_e (ncm_vector_get (phi, 0), ==, rho, 0.2, 0.0);

    for (order = 0; order < ncm_vector_len (pacf); order++)
      g_assert_cmpfloat (fabs (ncm_vector_get (pacf, order)), <, 1.0);

    ncm_vector_free (phi);
    ncm_vector_free (pacf);
  }

  ncm_vector_free (x);
  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* A correlation longer than the lag budget is resolved on a coarser level, not cut off
 * at the budget. */
void
test_ncm_stats_acorr_tau_cascade (void)
{
  NcmRNG *rng         = ncm_rng_pool_get ("test_ncm_stats_acorr");
  const gdouble rho   = 0.99;
  const gdouble exact = _test_acorr_ar1_tau (rho);
  const guint n       = 200000;
  const guint max_lag = 32;
  NcmVector *x        = ncm_vector_new (n);
  NcmStatsAcorr *acorr;

  _test_acorr_fill_ar1 (x, rho, rng);
  acorr = _test_acorr_from_vector (x, max_lag, 24, NCM_STATS_ACORR_METHOD_MAX);

  /* tau is six times the whole lag budget of a level. */
  g_assert_cmpfloat (exact, >, 6.0 * max_lag);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, exact, 0.25, 0.0);
  g_assert_cmpuint (ncm_stats_acorr_get_level (acorr, 0), >, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED, ==, 0);

  ncm_vector_free (x);
  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* Uncorrelated samples raise nothing. */
void
test_ncm_stats_acorr_diag_white (void)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, 424242);
  const guint n = 10000;
  NcmVector *x  = ncm_vector_new (n);
  NcmStatsAcorr *acorr;
  guint i;

  for (i = 0; i < n; i++)
    ncm_vector_set (x, i, ncm_rng_ugaussian_gen (rng));

  acorr = _test_acorr_from_vector (x, 512, 24, NCM_STATS_ACORR_METHOD_MAX);

  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0), ==, NCM_STATS_ACORR_DIAG_OK);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, 1.0, 0.25, 0.0);
  g_assert_cmpfloat (fabs (ncm_stats_acorr_get_drift_z (acorr, 0)), <, 3.0);

  ncm_vector_free (x);
  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* A series that only moves in one direction has an autocorrelation time of the order of
 * its own length: a handful of independent samples, reported as such. This is the
 * failure a fixed maximum lag hides, since it caps tau at the lag budget however long
 * the correlation is. */
void
test_ncm_stats_acorr_diag_ramp (void)
{
  const guint n        = 10000;
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, 512, 24, NCM_STATS_ACORR_METHOD_MAX);
  guint i;

  for (i = 0; i < n; i++)
    ncm_stats_acorr_update_var (acorr, 0, 1.0e-3 * i);

  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_SHORT_CHAIN, !=, 0);
  g_assert_cmpfloat (ncm_stats_acorr_get_tau (acorr, 0), >, 0.1 * n);
  g_assert_cmpfloat (ncm_stats_acorr_get_ess (acorr, 0), <, 10.0);
  g_assert_cmpuint (ncm_stats_acorr_get_level (acorr, 0), >, 0);

  ncm_stats_acorr_free (acorr);
}

/* A trend too small to show up in tau still moves the mean between the halves of the
 * series, and that is what the drift score is for. */
void
test_ncm_stats_acorr_diag_drift (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, 777);
  const guint n        = 10000;
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, 512, 24, NCM_STATS_ACORR_METHOD_MAX);
  guint i;

  for (i = 0; i < n; i++)
    ncm_stats_acorr_update_var (acorr, 0, ncm_rng_ugaussian_gen (rng) + 0.2 * i / (1.0 * n));

  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_DRIFT, !=, 0);
  g_assert_cmpfloat (fabs (ncm_stats_acorr_get_drift_z (acorr, 0)), >, 3.0);

  /* The correlation is still short: the trend is a fifth of a standard deviation. */
  g_assert_cmpfloat (ncm_stats_acorr_get_tau (acorr, 0), <, 5.0);

  /* Raising the threshold above the score clears the condition. */
  ncm_stats_acorr_set_drift_threshold (acorr, 1.0e3);
  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_DRIFT, ==, 0);

  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* With block averaging switched off there is no level able to resolve a long
 * correlation, and the estimate says so. */
void
test_ncm_stats_acorr_diag_truncated (void)
{
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, 31337);
  const guint n = 50000;
  NcmVector *x  = ncm_vector_new (n);
  NcmStatsAcorr *one_level, *cascade;

  _test_acorr_fill_ar1 (x, 0.999, rng);

  one_level = _test_acorr_from_vector (x, 512, 1,  NCM_STATS_ACORR_METHOD_MAX);
  cascade   = _test_acorr_from_vector (x, 512, 24, NCM_STATS_ACORR_METHOD_MAX);

  g_assert_cmpuint (ncm_stats_acorr_get_diag (one_level, 0) & NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED, !=, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_level (one_level, 0), ==, 0);
  g_assert_cmpuint (ncm_stats_acorr_nlevels (one_level, 0), ==, 1);

  g_assert_cmpuint (ncm_stats_acorr_get_diag (cascade, 0) & NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED, ==, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_level (cascade, 0), >, 0);

  ncm_vector_free (x);
  ncm_stats_acorr_free (one_level);
  ncm_stats_acorr_free (cascade);
  ncm_rng_free (rng);
}

/* An alternating series is where the two estimators of the largest-of-two method part
 * company: the auto-regressive fit reports the variance reduction, the initial sequence
 * is empty after the first pair. */
void
test_ncm_stats_acorr_diag_disagreement (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, 2718);
  const guint n        = 10000;
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, 512, 24, NCM_STATS_ACORR_METHOD_MAX);
  gdouble tau_ar, tau_geyer;
  guint i;

  for (i = 0; i < n; i++)
    ncm_stats_acorr_update_var (acorr, 0, ((i % 2 == 0) ? 1.0 : -1.0) + 0.1 * ncm_rng_ugaussian_gen (rng));

  tau_ar    = ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_AR);
  tau_geyer = ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_GEYER);

  g_assert_cmpfloat (tau_ar, <, 0.5 * tau_geyer);
  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT, !=, 0);

  /* The largest of the two is the one reported. */
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, GSL_MAX (tau_ar, tau_geyer), 1.0e-12, 0.0);

  /* A single estimator makes no claim about the other one. */
  ncm_stats_acorr_set_method (acorr, NCM_STATS_ACORR_METHOD_GEYER);
  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT, ==, 0);

  ncm_stats_acorr_free (acorr);
  ncm_rng_free (rng);
}

/* A series that never moves. */
void
test_ncm_stats_acorr_diag_zero_variance (void)
{
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, 16, 4, NCM_STATS_ACORR_METHOD_MAX);
  guint i;

  for (i = 0; i < 200; i++)
    ncm_stats_acorr_update_var (acorr, 0, 3.25);

  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_ZERO_VARIANCE, !=, 0);

  /* Never moved: perfectly correlated, so tau is reported at its cap and the effective
   * sample size is one, never the tau = 1 of a perfectly mixing chain. */
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, 200.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_ess (acorr, 0), ==, 1.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var (acorr, 0), ==, 0.0, 1.0e-14, 1.0e-14);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (acorr, 0), ==, 3.25, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_spec0 (acorr, 0), ==, 0.0, 1.0e-14, 1.0e-14);

  /* A series with no variance has no correlation time to report under any estimator. */
  {
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_GEYER), ==, 1.0, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau_method (acorr, 0, NCM_STATS_ACORR_METHOD_SOKAL), ==, 1.0, 1.0e-14, 0.0);
  }

  /* Too short to say anything about, as well. */
  {
    NcmStatsAcorr *few = ncm_stats_acorr_new_full (1, 16, 4, NCM_STATS_ACORR_METHOD_MAX);

    for (i = 0; i < 10; i++)
      ncm_stats_acorr_update_var (few, 0, 1.0);

    g_assert_cmpuint (ncm_stats_acorr_get_diag (few, 0) & NCM_STATS_ACORR_DIAG_SHORT_CHAIN, !=, 0);
    ncm_stats_acorr_free (few);
  }

  /* Varying, and still too short for any level to hold enough blocks for a drift score. */
  {
    NcmStatsAcorr *few = ncm_stats_acorr_new_full (1, 16, 4, NCM_STATS_ACORR_METHOD_MAX);

    for (i = 0; i < 10; i++)
      ncm_stats_acorr_update_var (few, 0, 1.0 * (i % 3));

    g_assert_cmpfloat (ncm_stats_acorr_get_var (few, 0), >, 0.0);
    ncm_assert_cmpdouble_e (ncm_stats_acorr_get_drift_z (few, 0), ==, 0.0, 1.0e-14, 1.0e-14);
    g_assert_cmpuint (ncm_stats_acorr_get_diag (few, 0) & NCM_STATS_ACORR_DIAG_SHORT_CHAIN, !=, 0);
    ncm_stats_acorr_free (few);
  }

  ncm_stats_acorr_free (acorr);
}

/* A chain whose burn-in has not been removed: its first part lives on a scale the rest
 * never returns to. The autocorrelation time of such a series is close to one, since the
 * excursion is most of its variance, so nothing in tau says the series is not one series.
 * The ratio of the two half variances does. */
void
test_ncm_stats_acorr_diag_variance_shift (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, 90210);
  const guint n        = 10000;
  NcmStatsAcorr *acorr = ncm_stats_acorr_new_full (1, 512, 24, NCM_STATS_ACORR_METHOD_MAX);
  NcmStatsAcorr *clean = ncm_stats_acorr_new_full (1, 512, 24, NCM_STATS_ACORR_METHOD_MAX);
  guint i;

  for (i = 0; i < n / 20; i++)
    ncm_stats_acorr_update_var (acorr, 0, 1.0e4 * ncm_rng_ugaussian_gen (rng));

  for (i = 0; i < n; i++)
  {
    const gdouble v = ncm_rng_ugaussian_gen (rng);

    ncm_stats_acorr_update_var (acorr, 0, v);
    ncm_stats_acorr_update_var (clean, 0, v);
  }

  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0) & NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT, !=, 0);
  g_assert_cmpfloat (ncm_stats_acorr_get_var_ratio (acorr, 0), <,
                     1.0 / NCM_STATS_ACORR_VARIANCE_SHIFT_FACTOR);

  /* tau alone says nothing: the excursion is most of the variance of the series. */
  g_assert_cmpfloat (ncm_stats_acorr_get_tau (acorr, 0), <, 10.0);

  /* The same samples without the excursion raise nothing. */
  g_assert_cmpuint (ncm_stats_acorr_get_diag (clean, 0), ==, NCM_STATS_ACORR_DIAG_OK);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var_ratio (clean, 0), ==, 1.0, 0.9, 0.0);

  ncm_stats_acorr_free (acorr);
  ncm_stats_acorr_free (clean);
  ncm_rng_free (rng);
}

/* Nothing accumulated, and a single sample. */
void
test_ncm_stats_acorr_diag_empty (void)
{
  NcmStatsAcorr *acorr           = ncm_stats_acorr_new (1);
  const NcmStatsAcorrDiag expect = NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_ZERO_VARIANCE;

  g_assert_cmpuint (ncm_stats_acorr_nitens (acorr, 0), ==, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0), ==, expect);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, 1.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (acorr, 0), ==, 0.0, 1.0e-14, 1.0e-14);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_ess (acorr, 0), ==, 0.0, 1.0e-14, 1.0e-14);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_var_mean (acorr, 0), ==, 0.0, 1.0e-14, 1.0e-14);

  ncm_stats_acorr_update_var (acorr, 0, 7.5);

  g_assert_cmpuint (ncm_stats_acorr_get_diag (acorr, 0), ==, expect);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_mean (acorr, 0), ==, 7.5, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_tau (acorr, 0), ==, 1.0, 1.0e-14, 0.0);
  g_assert_cmpuint (ncm_stats_acorr_get_level (acorr, 0), ==, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_window (acorr, 0), ==, 0);
  g_assert_cmpuint (ncm_stats_acorr_get_ar_order (acorr, 0), ==, 0);
  ncm_assert_cmpdouble_e (ncm_stats_acorr_get_drift_z (acorr, 0), ==, 0.0, 1.0e-14, 1.0e-14);

  ncm_stats_acorr_free (acorr);
}

/* Every condition has a name, and no condition has one. */
void
test_ncm_stats_acorr_diag_to_string (void)
{
  const struct
  {
    NcmStatsAcorrDiag diag;
    const gchar *str;
  } cases[] = {
    { NCM_STATS_ACORR_DIAG_OK,                  "ok"                   },
    { NCM_STATS_ACORR_DIAG_SHORT_CHAIN,         "short-chain"          },
    { NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED,    "window-truncated"     },
    { NCM_STATS_ACORR_DIAG_DRIFT,               "drift"                },
    { NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT, "method-disagreement"  },
    { NCM_STATS_ACORR_DIAG_ZERO_VARIANCE,       "zero-variance"        },
    { NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT,      "variance-shift"       },
    {
      NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_ZERO_VARIANCE,
      "short-chain,zero-variance"
    },
    {
      NCM_STATS_ACORR_DIAG_SHORT_CHAIN | NCM_STATS_ACORR_DIAG_WINDOW_TRUNCATED |
      NCM_STATS_ACORR_DIAG_DRIFT | NCM_STATS_ACORR_DIAG_METHOD_DISAGREEMENT |
      NCM_STATS_ACORR_DIAG_ZERO_VARIANCE | NCM_STATS_ACORR_DIAG_VARIANCE_SHIFT,
      "short-chain,window-truncated,drift,method-disagreement,zero-variance,variance-shift"
    },
  };

  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    gchar *str = ncm_stats_acorr_diag_to_string (cases[i].diag);

    g_assert_cmpstr (str, ==, cases[i].str);
    g_free (str);
  }
}

void
test_ncm_stats_acorr_traps (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/stats_acorr/invalid/var/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats_acorr/invalid/level/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats_acorr/invalid/method/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats_acorr/invalid/factor/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/stats_acorr/invalid/ar_crit/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_stats_acorr_invalid_var (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  ncm_stats_acorr_get_tau (test->acorr, 3);
}

void
test_ncm_stats_acorr_invalid_level (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  guint i;

  for (i = 0; i < 20; i++)
    ncm_stats_acorr_update_var (test->acorr, 0, 1.0 * i);

  ncm_vector_free (ncm_stats_acorr_get_acov (test->acorr, 0, 50));
}

void
test_ncm_stats_acorr_invalid_method (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  ncm_stats_acorr_set_method (test->acorr, NCM_STATS_ACORR_METHOD_LEN);
}

void
test_ncm_stats_acorr_invalid_ar_crit (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  ncm_stats_acorr_set_ar_criterion (test->acorr, NCM_STATS_ACORR_AR_CRIT_LEN);
}

void
test_ncm_stats_acorr_invalid_factor (TestNcmStatsAcorr *test, gconstpointer pdata)
{
  ncm_stats_acorr_set_reliability_factor (test->acorr, 0.5);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/stats_acorr/basic", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_basic,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/properties", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_properties,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/feed_equivalence", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_feed_equivalence,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/levels", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_levels,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/acf", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_acf,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/mean_error", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_mean_error,
              &test_ncm_stats_acorr_free);

  g_test_add_func ("/ncm/stats_acorr/estimator/ar1_analytic", &test_ncm_stats_acorr_est_ar1_analytic);
  g_test_add_func ("/ncm/stats_acorr/estimator/ar2_analytic", &test_ncm_stats_acorr_est_ar2_analytic);
  g_test_add_func ("/ncm/stats_acorr/estimator/white_analytic", &test_ncm_stats_acorr_est_white_analytic);
  g_test_add_func ("/ncm/stats_acorr/estimator/degenerate", &test_ncm_stats_acorr_est_degenerate);
  g_test_add_func ("/ncm/stats_acorr/estimator/not_positive_definite", &test_ncm_stats_acorr_est_not_positive_definite);
  g_test_add_func ("/ncm/stats_acorr/estimator/ar_criteria", &test_ncm_stats_acorr_est_ar_criteria);
  g_test_add_func ("/ncm/stats_acorr/ar_fit/coefficients", &test_ncm_stats_acorr_ar_fit_coefficients);

  g_test_add_func ("/ncm/stats_acorr/exact_vs_fft", &test_ncm_stats_acorr_exact_vs_fft);
  g_test_add_func ("/ncm/stats_acorr/acov_fft/all_lags", &test_ncm_stats_acorr_acov_fft_all_lags);
  g_test_add_func ("/ncm/stats_acorr/recenter", &test_ncm_stats_acorr_recenter);

  g_test_add_func ("/ncm/stats_acorr/tau/ar1_sampled", &test_ncm_stats_acorr_tau_ar1_sampled);
  g_test_add_func ("/ncm/stats_acorr/tau/cascade", &test_ncm_stats_acorr_tau_cascade);

  g_test_add_func ("/ncm/stats_acorr/diag/white", &test_ncm_stats_acorr_diag_white);
  g_test_add_func ("/ncm/stats_acorr/diag/ramp", &test_ncm_stats_acorr_diag_ramp);
  g_test_add_func ("/ncm/stats_acorr/diag/drift", &test_ncm_stats_acorr_diag_drift);
  g_test_add_func ("/ncm/stats_acorr/diag/truncated", &test_ncm_stats_acorr_diag_truncated);
  g_test_add_func ("/ncm/stats_acorr/diag/disagreement", &test_ncm_stats_acorr_diag_disagreement);
  g_test_add_func ("/ncm/stats_acorr/diag/zero_variance", &test_ncm_stats_acorr_diag_zero_variance);
  g_test_add_func ("/ncm/stats_acorr/diag/variance_shift", &test_ncm_stats_acorr_diag_variance_shift);
  g_test_add_func ("/ncm/stats_acorr/diag/empty", &test_ncm_stats_acorr_diag_empty);
  g_test_add_func ("/ncm/stats_acorr/diag/to_string", &test_ncm_stats_acorr_diag_to_string);

  g_test_add ("/ncm/stats_acorr/traps", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_traps,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/invalid/var/subprocess", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_invalid_var,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/invalid/level/subprocess", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_invalid_level,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/invalid/method/subprocess", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_invalid_method,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/invalid/factor/subprocess", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_invalid_factor,
              &test_ncm_stats_acorr_free);
  g_test_add ("/ncm/stats_acorr/invalid/ar_crit/subprocess", TestNcmStatsAcorr, NULL,
              &test_ncm_stats_acorr_new,
              &test_ncm_stats_acorr_invalid_ar_crit,
              &test_ncm_stats_acorr_free);

  g_test_run ();

  return 0;
}


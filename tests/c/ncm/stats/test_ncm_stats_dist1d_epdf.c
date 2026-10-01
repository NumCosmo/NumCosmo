/***************************************************************************
 *            test_ncm_stats_dist1d_epdf.c
 *
 *  Sat February 11 15:23:48 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br> & <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2017 <vitenti@uel.br>
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

typedef struct _TestNcmStatsDist1dEPDF
{
  NcmStatsDist1dEPDF *sd1;
  gdouble scale;
  gint max_tries;
} TestNcmStatsDist1dEPDF;

static void test_ncm_stats_dist1d_epdf_new (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_new_rot (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_new_small (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_new_large (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_gauss (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_gauss_reset (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_beta (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_isampling (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_serialize (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_gen (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_free (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);

static void test_ncm_stats_dist1d_epdf_mode_boundary (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);

static void test_ncm_stats_dist1d_epdf_traps (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_bw_regression (void);
static void test_ncm_stats_dist1d_epdf_bw_gauss_amise (void);
static void test_ncm_stats_dist1d_epdf_bounds (void);
static void test_ncm_stats_dist1d_epdf_empty_subprocess (void);
static void test_ncm_stats_dist1d_epdf_empty (void);
static void test_ncm_stats_dist1d_epdf_invalid_neg_weight (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_invalid_infinite_obs (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/stats_dist1d/epdf/bw_regression", &test_ncm_stats_dist1d_epdf_bw_regression);
  g_test_add_func ("/ncm/stats_dist1d/epdf/bw_gauss_amise", &test_ncm_stats_dist1d_epdf_bw_gauss_amise);
  g_test_add_func ("/ncm/stats_dist1d/epdf/bounds", &test_ncm_stats_dist1d_epdf_bounds);
  g_test_add_func ("/ncm/stats_dist1d/epdf/empty", &test_ncm_stats_dist1d_epdf_empty);
  g_test_add_func ("/ncm/stats_dist1d/epdf/empty/subprocess", &test_ncm_stats_dist1d_epdf_empty_subprocess);

  g_test_add ("/ncm/stats_dist1d/epdf/gauss", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_gauss,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/rot/gauss", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new_rot,
              &test_ncm_stats_dist1d_epdf_gauss,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/gauss/reset", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_gauss_reset,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/small/gauss", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new_small,
              &test_ncm_stats_dist1d_epdf_gauss,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/large/gauss", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new_large,
              &test_ncm_stats_dist1d_epdf_gauss,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/beta", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_beta,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/small/beta", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new_small,
              &test_ncm_stats_dist1d_epdf_beta,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/large/beta", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new_large,
              &test_ncm_stats_dist1d_epdf_beta,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/isampling", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_isampling,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/serialize", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_serialize,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/gen", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_gen,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/add/neg_weight/subprocess", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_invalid_neg_weight,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/add/infinite_obs/subprocess", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_invalid_infinite_obs,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/traps", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_traps,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_add ("/ncm/stats_dist1d/epdf/mode/boundary", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_mode_boundary,
              &test_ncm_stats_dist1d_epdf_free);

  g_test_run ();
}

#define TEST_NCM_STATS_DIST1D_EPDF_RELTOL (5.0e-1)

static void
test_ncm_stats_dist1d_epdf_gauss (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (test->sd1);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest   = 100000;
  const gdouble mu    = g_test_rand_double_range (-100.0, 100.0) * test->scale;
  const gdouble sigma = pow (10.0, g_test_rand_double_range (-2.0, 3.0)) * test->scale;
  const gdouble xl    = mu - 3.0 * sigma;
  const gdouble xu    = mu + 3.0 * sigma;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_gaussian_gen (rng, mu, sigma);

    ncm_stats_dist1d_epdf_add_obs (test->sd1, x);
  }

  {
    guint nobs  = 0;
    guint nobs0 = 0;

    g_object_get (test->sd1, "n-obs", &nobs, NULL);
    g_assert_cmpint (nobs, >, 0);

    ncm_stats_dist1d_epdf_add_obs_weight (test->sd1, 0.0, 0.0);
    g_object_get (test->sd1, "n-obs", &nobs0, NULL);
    g_assert_cmpint (nobs, ==, nobs0);
  }

  ncm_stats_dist1d_prepare (sd1);

  for (i = 0; i < ntest; i++)
  {
    const gdouble xi = xl + (xu - xl) / (1.0 * ntest - 1.0) * i;

    const gdouble ep_i = ncm_stats_dist1d_eval_p (sd1, xi);
    const gdouble ap_i = gsl_ran_gaussian_pdf (xi - mu, sigma);

    const gdouble epdf_i = ncm_stats_dist1d_eval_pdf (sd1, xi);
    const gdouble apdf_i = gsl_cdf_gaussian_P (xi - mu, sigma);

    /*printf ("% 22.15g % 22.15g % 22.15g %e % 22.15g % 22.15g %e\n", xi, ep_i, ap_i, fabs (ep_i / ap_i - 1.0), epdf_i, apdf_i, fabs (epdf_i / apdf_i - 1.0));*/
    ncm_assert_cmpdouble_e (epdf_i, ==, apdf_i, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
    ncm_assert_cmpdouble_e (ep_i,   ==, ap_i,   TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
  }

  {
    const gdouble xi          = ncm_stats_dist1d_get_xi (NCM_STATS_DIST1D (test->sd1));
    const gdouble xf          = ncm_stats_dist1d_get_xf (NCM_STATS_DIST1D (test->sd1));
    const gdouble cdf_minf_xl = gsl_cdf_gaussian_P (xi - mu, sigma);
    const gdouble cdf_minf_xu = gsl_cdf_gaussian_P (xf - mu, sigma);
    const gdouble isize       = xf - xi;

    g_assert_cmpfloat (isize, >, 0.0);

    for (i = 0; i < ntest; i++)
    {
      const gdouble ui       = 0.0 + 1.0 / (1.0 * ntest - 1.0) * i;
      const gdouble ui_gauss = cdf_minf_xl + (cdf_minf_xu - cdf_minf_xl) * ui;

      const gdouble einv_cdf_i = ncm_stats_dist1d_eval_inv_pdf (sd1, ui);
      const gdouble ainv_cdf_i = mu + gsl_cdf_gaussian_Pinv (ui_gauss, sigma);

      const gdouble einv_cdft_i = ncm_stats_dist1d_eval_inv_pdf_tail (sd1, ui);
      const gdouble ainv_cdft_i = mu + gsl_cdf_gaussian_Qinv (ui_gauss, sigma);

      /*printf ("[% 22.15g % 22.15g] (% 22.15g % 22.15g) <% 22.15g>: % 22.15g % 22.15g | % 22.15g % 22.15g\n", mu, sigma, ui, ui_gauss, isize * TEST_NCM_STATS_DIST1D_EPDF_RELTOL, einv_cdf_i, ainv_cdf_i, einv_cdft_i, ainv_cdft_i);*/

      ncm_assert_cmpdouble_e (einv_cdf_i,  ==, ainv_cdf_i,  TEST_NCM_STATS_DIST1D_EPDF_RELTOL, isize * TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
      ncm_assert_cmpdouble_e (einv_cdft_i, ==, ainv_cdft_i, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, isize * TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    }
  }

  {
    gdouble mode = ncm_stats_dist1d_eval_mode (sd1);

    if ((ncm_cmp (mode, mu, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0) != 0) && test->max_tries--)
    {
      NCM_TEST_FREE (ncm_rng_free, rng);
      ncm_stats_dist1d_epdf_reset (test->sd1);
      test_ncm_stats_dist1d_epdf_gauss (test, pdata);
    }
    else
    {
      ncm_assert_cmpdouble_e (mode, ==, mu, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
      NCM_TEST_FREE (ncm_rng_free, rng);
    }
  }

  {
    gdouble norma = 0.0;

    g_object_get (sd1, "norma", &norma, NULL);
    ncm_assert_cmpdouble_e (norma, !=, 0.0, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
    ncm_assert_cmpdouble_e (norma, ==, ncm_stats_dist1d_eval_norma (sd1), 1.0e-14, 0.0);
  }
}

static void
test_ncm_stats_dist1d_epdf_mode_boundary (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (test->sd1);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest   = 3000;
  guint i;

  /*
   * An exponential sample has a monotonically decreasing density: its true
   * mode sits at the domain edge, not strictly interior. This used to abort
   * the whole process (GSL_EINVAL from gsl_min_fminimizer_set's bracketing
   * precondition, wrapped in a fatal g_error); eval_mode() must now detect
   * that case and return the edge value instead.
   */
  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_exponential_gen (rng, 1.0);

    ncm_stats_dist1d_epdf_add_obs (test->sd1, x);
  }

  ncm_stats_dist1d_prepare (sd1);

  {
    const gdouble xi   = ncm_stats_dist1d_get_xi (sd1);
    const gdouble mode = ncm_stats_dist1d_eval_mode (sd1);

    ncm_assert_cmpdouble_e (mode, ==, xi, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
  }

  NCM_TEST_FREE (ncm_rng_free, rng);
}

static void
test_ncm_stats_dist1d_epdf_gauss_reset (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (test->sd1);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest   = 100000;
  const gdouble mu    = g_test_rand_double_range (-100.0, 100.0) * test->scale;
  const gdouble sigma = pow (10.0, g_test_rand_double_range (-2.0, 3.0)) * test->scale;
  const gdouble xl    = mu - 3.0 * sigma;
  const gdouble xu    = mu + 3.0 * sigma;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_uniform_gen (rng, 0.0, 1.0);

    ncm_stats_dist1d_epdf_add_obs (test->sd1, x);
  }

  ncm_stats_dist1d_prepare (sd1);

  ncm_stats_dist1d_epdf_reset (NCM_STATS_DIST1D_EPDF (sd1));

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_gaussian_gen (rng, mu, sigma);

    ncm_stats_dist1d_epdf_add_obs (test->sd1, x);
  }

  ncm_stats_dist1d_prepare (sd1);

  for (i = 0; i < ntest; i++)
  {
    const gdouble xi = xl + (xu - xl) / (1.0 * ntest - 1.0) * i;

    const gdouble ep_i = ncm_stats_dist1d_eval_p (sd1, xi);
    const gdouble ap_i = gsl_ran_gaussian_pdf (xi - mu, sigma);

    const gdouble epdf_i = ncm_stats_dist1d_eval_pdf (sd1, xi);
    const gdouble apdf_i = gsl_cdf_gaussian_P (xi - mu, sigma);

    /*printf ("% 22.15g % 22.15g % 22.15g %e % 22.15g % 22.15g %e\n", xi, ep_i, ap_i, fabs (ep_i / ap_i - 1.0), epdf_i, apdf_i, fabs (epdf_i / apdf_i - 1.0));*/
    ncm_assert_cmpdouble_e (epdf_i, ==, apdf_i, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
    ncm_assert_cmpdouble_e (ep_i,   ==, ap_i,   TEST_NCM_STATS_DIST1D_EPDF_RELTOL, 0.0);
  }

  {
    const gdouble xi          = ncm_stats_dist1d_get_xi (NCM_STATS_DIST1D (test->sd1));
    const gdouble xf          = ncm_stats_dist1d_get_xf (NCM_STATS_DIST1D (test->sd1));
    const gdouble cdf_minf_xl = gsl_cdf_gaussian_P (xi - mu, sigma);
    const gdouble cdf_minf_xu = gsl_cdf_gaussian_P (xf - mu, sigma);
    const gdouble isize       = xf - xi;

    g_assert_cmpfloat (isize, >, 0.0);

    for (i = 0; i < ntest; i++)
    {
      const gdouble ui       = 0.0 + 1.0 / (1.0 * ntest - 1.0) * i;
      const gdouble ui_gauss = cdf_minf_xl + (cdf_minf_xu - cdf_minf_xl) * ui;

      const gdouble einv_cdf_i = ncm_stats_dist1d_eval_inv_pdf (sd1, ui);
      const gdouble ainv_cdf_i = mu + gsl_cdf_gaussian_Pinv (ui_gauss, sigma);

      const gdouble einv_cdft_i = ncm_stats_dist1d_eval_inv_pdf_tail (sd1, ui);
      const gdouble ainv_cdft_i = mu + gsl_cdf_gaussian_Qinv (ui_gauss, sigma);

      /*printf ("[% 22.15g % 22.15g] (% 22.15g % 22.15g) <% 22.15g>: % 22.15g % 22.15g | % 22.15g % 22.15g\n", mu, sigma, ui, ui_gauss, isize * TEST_NCM_STATS_DIST1D_EPDF_RELTOL, einv_cdf_i, ainv_cdf_i, einv_cdft_i, ainv_cdft_i);*/

      ncm_assert_cmpdouble_e (einv_cdf_i,  ==, ainv_cdf_i,  TEST_NCM_STATS_DIST1D_EPDF_RELTOL, isize * TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
      ncm_assert_cmpdouble_e (einv_cdft_i, ==, ainv_cdft_i, TEST_NCM_STATS_DIST1D_EPDF_RELTOL, isize * TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    }
  }


  NCM_TEST_FREE (ncm_rng_free, rng);
}

static void
test_ncm_stats_dist1d_epdf_beta_meta (TestNcmStatsDist1dEPDF *test, gconstpointer pdata, const gdouble a, const gdouble b, const gdouble xl, const gdouble xu)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (test->sd1);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest   = 100000;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_beta_gen (rng, a, b) * test->scale;

    ncm_stats_dist1d_epdf_add_obs (test->sd1, x);
  }

  ncm_stats_dist1d_prepare (sd1);

  for (i = 0; i < ntest; i++)
  {
    const gdouble xi = xl + (xu - xl) / (1.0 * ntest - 1.0) * i;

    const gdouble ep_i = ncm_stats_dist1d_eval_p (sd1, xi * test->scale) * test->scale;
    const gdouble ap_i = gsl_ran_beta_pdf (xi, a, b);

    const gdouble epdf_i = ncm_stats_dist1d_eval_pdf (sd1, xi * test->scale);
    const gdouble apdf_i = gsl_cdf_beta_P (xi, a, b);

    /*printf ("% 22.15g % 22.15g % 22.15g %e % 22.15g % 22.15g %e\n", xi, ep_i, ap_i, fabs (ep_i / ap_i - 1.0), epdf_i, apdf_i, fabs (epdf_i / apdf_i - 1.0));*/

    if (apdf_i == 0.0)
      g_assert_cmpfloat (fabs (epdf_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else if (epdf_i == 0.0)
      g_assert_cmpfloat (fabs (apdf_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else
      g_assert_cmpfloat (fabs (epdf_i / apdf_i - 1.0), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);

    if (ap_i == 0.0)
      g_assert_cmpfloat (fabs (ep_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else if (ep_i == 0.0)
      g_assert_cmpfloat (fabs (ap_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else
      g_assert_cmpfloat (fabs (ep_i / ap_i - 1.0), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
  }

  NCM_TEST_FREE (ncm_rng_free, rng);
}

static void
test_ncm_stats_dist1d_epdf_beta (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  test_ncm_stats_dist1d_epdf_beta_meta (test, pdata, 0.5, 0.5, 0.2, 0.8);
  ncm_stats_dist1d_epdf_reset (test->sd1);

  test_ncm_stats_dist1d_epdf_beta_meta (test, pdata, 2.0, 5.0, 0.2, 0.8);
  ncm_stats_dist1d_epdf_reset (test->sd1);

  test_ncm_stats_dist1d_epdf_beta_meta (test, pdata, 2.0, 2.0, 0.1, 0.9);
  ncm_stats_dist1d_epdf_reset (test->sd1);
}

static void
test_ncm_stats_dist1d_epdf_isampling (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (test->sd1);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest   = 100000;
  const gdouble a     = 2.0;
  const gdouble b     = 5.0;
  const gdouble xl    = 0.1;
  const gdouble xu    = 0.9;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_uniform_gen (rng, 0.0, 1.0);

    ncm_stats_dist1d_epdf_add_obs_weight (test->sd1, x, gsl_ran_beta_pdf (x, a, b));
  }

  ncm_stats_dist1d_prepare (sd1);

  for (i = 0; i < ntest; i++)
  {
    const gdouble xi = xl + (xu - xl) / (1.0 * ntest - 1.0) * i;

    const gdouble ep_i = ncm_stats_dist1d_eval_p (sd1, xi);
    const gdouble ap_i = gsl_ran_beta_pdf (xi, a, b);

    const gdouble epdf_i = ncm_stats_dist1d_eval_pdf (sd1, xi);
    const gdouble apdf_i = gsl_cdf_beta_P (xi, a, b);

    /*printf ("% 22.15g % 22.15g % 22.15g %e % 22.15g % 22.15g %e\n", xi, ep_i, ap_i, fabs (ep_i / ap_i - 1.0), epdf_i, apdf_i, fabs (epdf_i / apdf_i - 1.0));*/

    if (apdf_i == 0.0)
      g_assert_cmpfloat (fabs (epdf_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else if (epdf_i == 0.0)
      g_assert_cmpfloat (fabs (apdf_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else
      g_assert_cmpfloat (fabs (epdf_i / apdf_i - 1.0), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);

    if (ap_i == 0.0)
      g_assert_cmpfloat (fabs (ep_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else if (ep_i == 0.0)
      g_assert_cmpfloat (fabs (ap_i), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
    else
      g_assert_cmpfloat (fabs (ep_i / ap_i - 1.0), <, TEST_NCM_STATS_DIST1D_EPDF_RELTOL);
  }

  NCM_TEST_FREE (ncm_rng_free, rng);
}

static void
test_ncm_stats_dist1d_epdf_serialize (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1d *sd1     = NCM_STATS_DIST1D (test->sd1);
  NcmSerialize *ser       = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmStatsDist1d *sd1_dup = NCM_STATS_DIST1D (ncm_serialize_dup_obj (ser, G_OBJECT (sd1)));
  NcmRNG *rng             = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest       = 100000;
  const gdouble mu        = g_test_rand_double_range (-100.0, 100.0) * test->scale;
  const gdouble sigma     = pow (10.0, g_test_rand_double_range (-2.0, 3.0)) * test->scale;
  const gdouble xl        = mu - 3.0 * sigma;
  const gdouble xu        = mu + 3.0 * sigma;
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_gaussian_gen (rng, mu, sigma);

    ncm_stats_dist1d_epdf_add_obs (NCM_STATS_DIST1D_EPDF (sd1), x);
    ncm_stats_dist1d_epdf_add_obs (NCM_STATS_DIST1D_EPDF (sd1_dup), x);
  }

  ncm_stats_dist1d_prepare (sd1);
  ncm_stats_dist1d_prepare (sd1_dup);

  for (i = 0; i < ntest; i++)
  {
    const gdouble xi = xl + (xu - xl) / (1.0 * ntest - 1.0) * i;

    const gdouble ep_i = ncm_stats_dist1d_eval_p (sd1, xi);
    const gdouble ap_i = ncm_stats_dist1d_eval_p (sd1_dup, xi);

    const gdouble epdf_i = ncm_stats_dist1d_eval_pdf (sd1, xi);
    const gdouble apdf_i = ncm_stats_dist1d_eval_pdf (sd1_dup, xi);

    const gdouble em2lnp_i = ncm_stats_dist1d_eval_m2lnp (sd1, xi);
    const gdouble am2lnp_i = ncm_stats_dist1d_eval_m2lnp (sd1_dup, xi);

    ncm_assert_cmpdouble_e (epdf_i,   ==, apdf_i,   1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ep_i,     ==, ap_i,     1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (em2lnp_i, ==, am2lnp_i, 1.0e-14, 0.0);
  }

  ncm_rng_free (rng);
  ncm_stats_dist1d_free (sd1_dup);
  ncm_serialize_free (ser);
}

static void
test_ncm_stats_dist1d_epdf_gen (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (test->sd1);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint ntest   = 100000;
  const gdouble mu    = g_test_rand_double_range (-100.0, 100.0) * test->scale;
  const gdouble sigma = pow (10.0, g_test_rand_double_range (-2.0, 3.0)) * test->scale;
  NcmStatsVec *svec   = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  guint i;

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_rng_gaussian_gen (rng, mu, sigma);

    ncm_stats_dist1d_epdf_add_obs (NCM_STATS_DIST1D_EPDF (sd1), x);
  }

  ncm_stats_dist1d_prepare (sd1);

  for (i = 0; i < ntest; i++)
  {
    const gdouble x = ncm_stats_dist1d_gen (sd1, rng);

    ncm_stats_vec_set (svec, 0, x);
    ncm_stats_vec_update_weight (svec, 1.0);
  }

  if ((ncm_cmp (ncm_stats_vec_get_mean (svec, 0), mu, 1.0e-2, 0.0) != 0) && test->max_tries--)
  {
    ncm_rng_free (rng);
    ncm_stats_vec_free (svec);
    ncm_stats_dist1d_epdf_reset (test->sd1);
    test_ncm_stats_dist1d_epdf_gen (test, pdata);
  }
  else
  {
    ncm_assert_cmpdouble_e (ncm_stats_vec_get_mean (svec, 0), ==, mu, 1.0e-2, 0.0);
    ncm_rng_free (rng);
    ncm_stats_vec_free (svec);
  }
}

static void
test_ncm_stats_dist1d_epdf_new (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  test->sd1       = ncm_stats_dist1d_epdf_new_full (50000, NCM_STATS_DIST1D_EPDF_BW_AUTO, 0.1, 1.0e-2);
  test->scale     = 1.0;
  test->max_tries = 5;

  g_assert_true (NCM_IS_STATS_DIST1D (test->sd1));
  g_assert_true (NCM_IS_STATS_DIST1D_EPDF (test->sd1));
}

static void
test_ncm_stats_dist1d_epdf_new_rot (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  test->sd1       = ncm_stats_dist1d_epdf_new_full (50000, NCM_STATS_DIST1D_EPDF_BW_RoT, 0.1, 1.0e-2);
  test->scale     = 1.0;
  test->max_tries = 5;

  g_assert_true (NCM_IS_STATS_DIST1D (test->sd1));
  g_assert_true (NCM_IS_STATS_DIST1D_EPDF (test->sd1));
}

static void
test_ncm_stats_dist1d_epdf_new_small (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  test->sd1       = ncm_stats_dist1d_epdf_new (1.0e-2);
  test->scale     = 1.0e-10;
  test->max_tries = 5;

  g_assert_true (NCM_IS_STATS_DIST1D (test->sd1));
  g_assert_true (NCM_IS_STATS_DIST1D_EPDF (test->sd1));
}

static void
test_ncm_stats_dist1d_epdf_new_large (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  test->sd1       = ncm_stats_dist1d_epdf_new (1.0e-2);
  test->scale     = 1.0e-10;
  test->max_tries = 5;

  g_assert_true (NCM_IS_STATS_DIST1D (test->sd1));
  g_assert_true (NCM_IS_STATS_DIST1D_EPDF (test->sd1));
}

static void
test_ncm_stats_dist1d_epdf_free (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1dEPDF *sd1  = test->sd1;
  NcmStatsDist1d *sd1_base = NCM_STATS_DIST1D (sd1);

  ncm_stats_dist1d_ref (sd1_base);
  ncm_stats_dist1d_clear (&sd1_base);
  g_assert_true (sd1_base == NULL);

  sd1_base = NCM_STATS_DIST1D (sd1);

  NCM_TEST_FREE (ncm_stats_dist1d_free, sd1_base);
}

static void
test_ncm_stats_dist1d_epdf_traps (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/stats_dist1d/epdf/add/neg_weight/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*invalid observation*");

  g_test_trap_subprocess ("/ncm/stats_dist1d/epdf/add/infinite_obs/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*invalid observation*");
}

static void
test_ncm_stats_dist1d_epdf_invalid_neg_weight (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  ncm_stats_dist1d_epdf_add_obs_weight (test->sd1, 0.0, -1.0);
}

static void
test_ncm_stats_dist1d_epdf_invalid_infinite_obs (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  ncm_stats_dist1d_epdf_add_obs_weight (test->sd1, GSL_POSINF, 1.0);
}

/* Draws n values of one of the regression shapes */
static void
_test_epdf_fill (NcmStatsDist1dEPDF *epdf, guint shape, guint n, NcmRNG *rng)
{
  guint i;

  for (i = 0; i < n; i++)
  {
    switch (shape)
    {
      case 0: /* weighted unit Gaussian */
      {
        /* Drawn in sequence: C leaves the order of argument evaluation unspecified */
        const gdouble w = ncm_rng_uniform_gen (rng, 0.5, 2.0);
        const gdouble x = ncm_rng_gaussian_gen (rng, 0.0, 1.0);

        ncm_stats_dist1d_epdf_add_obs_weight (epdf, x, w);
        break;
      }
      case 1: /* claw of Marron and Wand */
      {
        const gdouble u = ncm_rng_uniform01_gen (rng);
        const guint k   = (u < 0.5) ? 0 : 1 + GSL_MIN ((guint) ((u - 0.5) * 10.0), 4);

        ncm_stats_dist1d_epdf_add_obs (epdf, (k == 0) ? ncm_rng_gaussian_gen (rng, 0.0, 1.0) : ncm_rng_gaussian_gen (rng, 0.5 * (k - 1) - 1.0, 0.1));
        break;
      }
      case 2: /* unit exponential */
        ncm_stats_dist1d_epdf_add_obs (epdf, ncm_rng_exponential_gen (rng, 1.0));
        break;
      case 3: /* uniform on [2, 5] */
        ncm_stats_dist1d_epdf_add_obs (epdf, ncm_rng_uniform_gen (rng, 2.0, 5.0));
        break;
      default: /* arcsine, Beta (1/2, 1/2) */
        ncm_stats_dist1d_epdf_add_obs (epdf, ncm_rng_beta_gen (rng, 0.5, 0.5));
        break;
    }
  }
}

/*
 * Regression values of the automatic bandwidth, 20000 observations per shape, merging at
 * sd-min-scale 1e-3 while they are added (max-obs 1000). Compared to 1e-12 rather than bit
 * for bit: FFTW may pick another algorithm on another machine.
 */
static void
test_ncm_stats_dist1d_epdf_bw_regression (void)
{
  const gdouble h_ref[5] = {0.16568181535100246, 0.14716458563415075, 0.12453004358786544, 0.11489972618980605, 0.02824853152261754};
  guint shape;

  for (shape = 0; shape < 5; shape++)
  {
    NcmRNG *rng              = ncm_rng_seeded_new (NULL, 1234);
    NcmStatsDist1dEPDF *epdf = ncm_stats_dist1d_epdf_new_full (1000, NCM_STATS_DIST1D_EPDF_BW_AUTO, 0.1, 1.0e-3);

    _test_epdf_fill (epdf, shape, 20000, rng);
    ncm_stats_dist1d_set_compute_cdf (NCM_STATS_DIST1D (epdf), FALSE);
    ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf));

    g_test_message ("shape %u: h = %.17g", shape, ncm_stats_dist1d_get_current_h (NCM_STATS_DIST1D (epdf)));
    ncm_assert_cmpdouble_e (ncm_stats_dist1d_get_current_h (NCM_STATS_DIST1D (epdf)), ==, h_ref[shape], 1.0e-12, 0.0);

    ncm_stats_dist1d_epdf_free (epdf);
    ncm_rng_free (rng);
  }
}

/* Unit Gaussian, N = 1e5: the automatic bandwidth against (4 / (3 N))^(1/5), measured 1.007 to 1.015 */
static void
test_ncm_stats_dist1d_epdf_bw_gauss_amise (void)
{
  const guint N         = 100000;
  const gdouble h_amise = pow (4.0 / (3.0 * N), 0.2);
  guint seed;

  for (seed = 1; seed <= 3; seed++)
  {
    NcmRNG *rng              = ncm_rng_seeded_new (NULL, seed);
    NcmStatsDist1dEPDF *epdf = ncm_stats_dist1d_epdf_new_full (2 * N, NCM_STATS_DIST1D_EPDF_BW_AUTO, 0.1, 1.0e-6);
    guint i;

    for (i = 0; i < N; i++)
      ncm_stats_dist1d_epdf_add_obs (epdf, ncm_rng_gaussian_gen (rng, 0.0, 1.0));

    ncm_stats_dist1d_set_compute_cdf (NCM_STATS_DIST1D (epdf), FALSE);
    ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf));
    g_test_message ("seed %u: h / h_AMISE = %.4f", seed, ncm_stats_dist1d_get_current_h (NCM_STATS_DIST1D (epdf)) / h_amise);
    ncm_assert_cmpdouble_e (ncm_stats_dist1d_get_current_h (NCM_STATS_DIST1D (epdf)) / h_amise, ==, 1.0, 0.0, 0.05);

    ncm_stats_dist1d_epdf_free (epdf);
    ncm_rng_free (rng);
  }
}

/* Moving a bound after a prepare changes the support and the automatic bandwidth */
static void
test_ncm_stats_dist1d_epdf_bounds (void)
{
  NcmRNG *rng              = ncm_rng_seeded_new (NULL, 1);
  NcmStatsDist1dEPDF *epdf = ncm_stats_dist1d_epdf_new (1.0e-3);
  NcmStatsDist1d *sd1      = NCM_STATS_DIST1D (epdf);
  gdouble h0;

  _test_epdf_fill (epdf, 3, 5000, rng);
  ncm_stats_dist1d_prepare (sd1);
  h0 = ncm_stats_dist1d_get_current_h (sd1);

  ncm_stats_dist1d_epdf_set_min (epdf, 0.0);
  ncm_stats_dist1d_epdf_set_max (epdf, 7.0);
  ncm_stats_dist1d_prepare (sd1);

  g_assert_cmpfloat (ncm_stats_dist1d_get_xi (sd1), ==, 0.0);
  g_assert_cmpfloat (ncm_stats_dist1d_get_xf (sd1), ==, 7.0);
  g_assert_cmpfloat (ncm_stats_dist1d_get_current_h (sd1), !=, h0);

  ncm_stats_dist1d_epdf_reset (epdf);
  g_assert_cmpuint (ncm_stats_dist1d_epdf_get_bw_type (epdf), ==, NCM_STATS_DIST1D_EPDF_BW_AUTO);

  ncm_stats_dist1d_epdf_free (epdf);
  ncm_rng_free (rng);
}

static void
test_ncm_stats_dist1d_epdf_empty_subprocess (void)
{
  NcmStatsDist1dEPDF *epdf = ncm_stats_dist1d_epdf_new (1.0e-3);

  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (epdf));
}

static void
test_ncm_stats_dist1d_epdf_empty (void)
{
  g_test_trap_subprocess ("/ncm/stats_dist1d/epdf/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*no observations*");
}


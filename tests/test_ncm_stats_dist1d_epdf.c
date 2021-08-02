/***************************************************************************
 *            test_ncm_stats_dist1d_epdf.c
 *
 *  Sat February 11 15:23:48 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br> & <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2017 <sandro@isoftware.com.br>
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
} TestNcmStatsDist1dEPDF;

static void test_ncm_stats_dist1d_epdf_new (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_gauss (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_beta (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_isampling (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_free (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);

static void test_ncm_stats_dist1d_epdf_traps (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_invalid_neg_weight (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);
static void test_ncm_stats_dist1d_epdf_invalid_infinite_obs (TestNcmStatsDist1dEPDF *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_set_nonfatal_assertions ();
  
  g_test_add ("/ncm/stats_dist1d/epdf/gauss", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_gauss,
              &test_ncm_stats_dist1d_epdf_free);
  
  g_test_add ("/ncm/stats_dist1d/epdf/beta", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_beta,
              &test_ncm_stats_dist1d_epdf_free);
  
  g_test_add ("/ncm/stats_dist1d/epdf/isampling", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_isampling,
              &test_ncm_stats_dist1d_epdf_free);
  
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/stats_dist1d/epdf/add/neg_weight/subprocess", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_invalid_neg_weight,
              &test_ncm_stats_dist1d_epdf_free);
  
  g_test_add ("/ncm/stats_dist1d/epdf/add/infinite_obs/subprocess", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_invalid_infinite_obs,
              &test_ncm_stats_dist1d_epdf_free);
#endif
  
  g_test_add ("/ncm/stats_dist1d/epdf/traps", TestNcmStatsDist1dEPDF, NULL,
              &test_ncm_stats_dist1d_epdf_new,
              &test_ncm_stats_dist1d_epdf_traps,
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
  const gdouble mu    = g_test_rand_double_range (-100.0, 100.0);
  const gdouble sigma = pow (10.0, g_test_rand_double_range (-2.0, 3.0));
  const gdouble xl    = mu - 3.0 * sigma;
  const gdouble xu    = mu + 3.0 * sigma;
  guint i;
  
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
    const gdouble x = ncm_rng_beta_gen (rng, a, b);
    
    ncm_stats_dist1d_epdf_add_obs (test->sd1, x);
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
test_ncm_stats_dist1d_epdf_new (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  test->sd1 = ncm_stats_dist1d_epdf_new (1.0e-2);
  
  g_assert_true (NCM_IS_STATS_DIST1D (test->sd1));
  g_assert_true (NCM_IS_STATS_DIST1D_EPDF (test->sd1));
}

static void
test_ncm_stats_dist1d_epdf_free (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
  NcmStatsDist1dEPDF *sd1 = test->sd1;
  
  NCM_TEST_FREE (ncm_stats_dist1d_free, NCM_STATS_DIST1D (sd1));
}

static void
test_ncm_stats_dist1d_epdf_traps (TestNcmStatsDist1dEPDF *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/stats_dist1d/epdf/add/neg_weight/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/stats_dist1d/epdf/add/infinite_obs/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
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


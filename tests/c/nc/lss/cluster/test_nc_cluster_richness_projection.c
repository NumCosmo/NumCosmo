/***************************************************************************
 *            test_nc_cluster_richness_projection.c
 *
 *  Wed September 17 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_integration.h>

#define TEST_MU (log (30.0))
#define TEST_SIGMA (0.33)
#define TEST_LNL_MIN (log (1.0))
#define TEST_LNL_MAX (log (300.0))

typedef struct _TestNcClusterRichnessProjection
{
  NcClusterRichnessProjection *crp;
  gdouble mu;
  gdouble sigma;
  gdouble tau;
} TestNcClusterRichnessProjection;

void test_nc_cluster_richness_projection_new (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_free (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_properties (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_convolution (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_original_form (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_eval_lnlambda (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_eval_int (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_norma (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_lnnormal_limit (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_invalid_range (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_invalid_reltol (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_invalid_sigma (TestNcClusterRichnessProjection *test, gconstpointer pdata);
void test_nc_cluster_richness_projection_out_of_range (TestNcClusterRichnessProjection *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/cluster_richness_projection/properties", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_properties,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/convolution", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_convolution,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/original_form", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_original_form,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/eval_lnlambda", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_eval_lnlambda,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/eval_int", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_eval_int,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/norma", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_norma,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/lnnormal_limit", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_lnnormal_limit,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/invalid/range", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_invalid_range,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/invalid/reltol", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_invalid_reltol,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/invalid/sigma", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_invalid_sigma,
              &test_nc_cluster_richness_projection_free);
  g_test_add ("/nc/cluster_richness_projection/invalid/out_of_range", TestNcClusterRichnessProjection, NULL,
              &test_nc_cluster_richness_projection_new,
              &test_nc_cluster_richness_projection_out_of_range,
              &test_nc_cluster_richness_projection_free);

  g_test_run ();

  return 0;
}

void
test_nc_cluster_richness_projection_new (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  test->crp   = nc_cluster_richness_projection_new ();
  test->mu    = TEST_MU;
  test->sigma = TEST_SIGMA;
  test->tau   = 0.2;

  g_assert_true (NC_IS_CLUSTER_RICHNESS_PROJECTION (test->crp));

  nc_cluster_richness_projection_set_lnlambda_range (test->crp, TEST_LNL_MIN, TEST_LNL_MAX);
}

void
test_nc_cluster_richness_projection_free (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_cluster_richness_projection_free, test->crp);
}

/* Reference quadratures, deliberately independent of the class internals. */

typedef struct _TestRefParams
{
  gdouble mu;
  gdouble sigma;
  gdouble tau;
  gdouble lnlambda;
} TestRefParams;

/*
 * Integrand of T(lambda) = tau int_0^lambda f_LN(t) e^{-tau (lambda - t)} dt, in
 * u = ln t so that the log-normal lower tail is resolved.
 */
static gdouble
_test_conv_integrand (gdouble u, gpointer p)
{
  TestRefParams *pars = (TestRefParams *) p;
  const gdouble x     = (u - pars->mu) / pars->sigma;
  const gdouble lnl   = pars->lnlambda;

  return exp (-0.5 * x * x) / (ncm_c_sqrt_2pi () * pars->sigma)
         * exp (-pars->tau * (exp (lnl) - exp (u)));
}

static gdouble
_test_conv_quad (TestRefParams *pars)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (2000);
  gsl_function F;
  gdouble res, err;

  F.function = &_test_conv_integrand;
  F.params   = pars;

  gsl_integration_qag (&F, MIN (pars->mu, pars->lnlambda) - 12.0 * pars->sigma, pars->lnlambda,
                       0.0, 1.0e-13, 2000, GSL_INTEG_GAUSS61, w, &res, &err);
  gsl_integration_workspace_free (w);

  return pars->tau * res;
}

/* Integrand of the expression as originally written, in y. */
static gdouble
_test_original_integrand (gdouble y, gpointer p)
{
  TestRefParams *pars = (TestRefParams *) p;

  return exp (-y * y + pars->tau * exp (pars->mu + M_SQRT2 * pars->sigma * y));
}

static gdouble
_test_original_quad (TestRefParams *pars)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (2000);
  const gdouble L              = (pars->lnlambda - pars->mu) / (M_SQRT2 * pars->sigma);
  gsl_function F;
  gdouble res, err;

  F.function = &_test_original_integrand;
  F.params   = pars;

  gsl_integration_qag (&F, MIN (-12.0, L - 12.0), L, 0.0, 1.0e-13, 2000, GSL_INTEG_GAUSS61, w, &res, &err);
  gsl_integration_workspace_free (w);

  return pars->tau / sqrt (M_PI) * exp (-pars->tau * exp (pars->lnlambda)) * res;
}

/* Integrand of int T dlambda, in u = ln lambda. */
static gdouble
_test_int_integrand (gdouble u, gpointer p)
{
  NcClusterRichnessProjection *crp = NC_CLUSTER_RICHNESS_PROJECTION (p);

  return nc_cluster_richness_projection_eval_lnlambda (crp, u);
}

static gdouble
_test_int_quad (NcClusterRichnessProjection *crp, gdouble lnl_lo, gdouble lnl_hi)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (2000);
  gsl_function F;
  gdouble res, err;

  F.function = &_test_int_integrand;
  F.params   = crp;

  gsl_integration_qag (&F, lnl_lo, lnl_hi, 0.0, 1.0e-12, 2000, GSL_INTEG_GAUSS61, w, &res, &err);
  gsl_integration_workspace_free (w);

  return res;
}

void
test_nc_cluster_richness_projection_properties (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  NcClusterRichnessProjection *crp2;
  gdouble lnl_min, lnl_max, reltol;

  g_object_get (test->crp,
                "lnlambda-min", &lnl_min,
                "lnlambda-max", &lnl_max,
                "reltol", &reltol,
                NULL);

  ncm_assert_cmpdouble (lnl_min, ==, TEST_LNL_MIN);
  ncm_assert_cmpdouble (lnl_max, ==, TEST_LNL_MAX);
  ncm_assert_cmpdouble (reltol, ==, NC_CLUSTER_RICHNESS_PROJECTION_DEFAULT_RELTOL);
  ncm_assert_cmpdouble (nc_cluster_richness_projection_get_reltol (test->crp), ==, reltol);

  g_object_set (test->crp,
                "lnlambda-min", log (2.0),
                "lnlambda-max", log (400.0),
                "reltol", 1.0e-9,
                NULL);

  g_object_get (test->crp,
                "lnlambda-min", &lnl_min,
                "lnlambda-max", &lnl_max,
                "reltol", &reltol,
                NULL);

  ncm_assert_cmpdouble (lnl_min, ==, log (2.0));
  ncm_assert_cmpdouble (lnl_max, ==, log (400.0));
  ncm_assert_cmpdouble (reltol, ==, 1.0e-9);

  nc_cluster_richness_projection_set_reltol (test->crp, 1.0e-10);
  ncm_assert_cmpdouble (nc_cluster_richness_projection_get_reltol (test->crp), ==, 1.0e-10);

  crp2 = nc_cluster_richness_projection_ref (test->crp);
  g_assert_true (crp2 == test->crp);
  nc_cluster_richness_projection_clear (&crp2);
  g_assert_true (crp2 == NULL);

  g_assert_true (NC_IS_CLUSTER_RICHNESS_PROJECTION (test->crp));
}

void
test_nc_cluster_richness_projection_convolution (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  const gdouble tau_a[3]    = { 0.05, 0.2, 1.0 };
  const gdouble sigma_a[2]  = { 0.2, 0.5 };
  const gdouble lambda_a[5] = { 10.0, 20.0, 30.0, 50.0, 100.0 };
  guint i, j, k;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 2; j++)
    {
      nc_cluster_richness_projection_prepare (test->crp, test->mu, sigma_a[j], tau_a[i]);

      for (k = 0; k < 5; k++)
      {
        TestRefParams pars = { test->mu, sigma_a[j], tau_a[i], log (lambda_a[k]) };
        const gdouble T    = nc_cluster_richness_projection_eval (test->crp, pars.lnlambda);
        const gdouble ref  = _test_conv_quad (&pars);

        ncm_assert_cmpdouble_e (T, ==, ref, 1.0e-8, 0.0);
      }
    }
  }
}

void
test_nc_cluster_richness_projection_original_form (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  /* tau * lambda_max is kept moderate so that the original integrand, which grows
   * as e^{tau lambda}, stays inside the double range. */
  const gdouble tau         = 0.05;
  const gdouble lambda_a[4] = { 10.0, 20.0, 50.0, 100.0 };
  guint k;

  nc_cluster_richness_projection_prepare (test->crp, test->mu, test->sigma, tau);

  for (k = 0; k < 4; k++)
  {
    TestRefParams pars = { test->mu, test->sigma, tau, log (lambda_a[k]) };
    const gdouble T    = nc_cluster_richness_projection_eval (test->crp, pars.lnlambda);
    const gdouble ref  = _test_original_quad (&pars);

    ncm_assert_cmpdouble_e (T, ==, ref, 1.0e-8, 0.0);
  }
}

void
test_nc_cluster_richness_projection_eval_lnlambda (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  const gdouble lambda_a[3] = { 5.0, 40.0, 200.0 };
  guint k;

  nc_cluster_richness_projection_prepare (test->crp, test->mu, test->sigma, test->tau);

  for (k = 0; k < 3; k++)
  {
    const gdouble lnl = log (lambda_a[k]);

    ncm_assert_cmpdouble_e (nc_cluster_richness_projection_eval_lnlambda (test->crp, lnl), ==,
                            lambda_a[k] * nc_cluster_richness_projection_eval (test->crp, lnl),
                            1.0e-14, 0.0);
  }
}

void
test_nc_cluster_richness_projection_eval_int (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  const gdouble tau_a[3]    = { 0.05, 0.2, 1.0 };
  const gdouble bin_a[4][2] = {
    { 10.0, 20.0 }, { 20.0, 40.0 }, { 40.0, 80.0 }, { 80.0, 300.0 }
  };
  guint i, k;

  for (i = 0; i < 3; i++)
  {
    nc_cluster_richness_projection_prepare (test->crp, test->mu, test->sigma, tau_a[i]);

    for (k = 0; k < 4; k++)
    {
      const gdouble lnl_lo = log (bin_a[k][0]);
      const gdouble lnl_hi = log (bin_a[k][1]);
      const gdouble I_bin  = nc_cluster_richness_projection_eval_int (test->crp, lnl_lo, lnl_hi);
      const gdouble ref    = _test_int_quad (test->crp, lnl_lo, lnl_hi);

      ncm_assert_cmpdouble_e (I_bin, ==, ref, 1.0e-8, 1.0e-20);
    }
  }

  /* A degenerate bin integrates to zero. */
  ncm_assert_cmpdouble_e (nc_cluster_richness_projection_eval_int (test->crp, log (10.0), log (10.0)),
                          ==, 0.0, 0.0, 1.0e-16);
}

void
test_nc_cluster_richness_projection_norma (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  /* T is a normalized density: over a range wide enough to hold both the
   * log-normal and the exponential tail its integral is one. */
  const gdouble tau_a[2] = { 0.5, 5.0 };
  guint i;

  nc_cluster_richness_projection_set_lnlambda_range (test->crp, log (1.0e-6), log (1.0e4));

  for (i = 0; i < 2; i++)
  {
    nc_cluster_richness_projection_prepare (test->crp, test->mu, test->sigma, tau_a[i]);

    ncm_assert_cmpdouble_e (nc_cluster_richness_projection_eval_int (test->crp, log (1.0e-6), log (1.0e4)),
                            ==, 1.0, 1.0e-9, 0.0);
  }
}

void
test_nc_cluster_richness_projection_lnnormal_limit (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  /* As tau grows the added richness vanishes and T collapses onto f_LN. */
  const gdouble tau         = 1.0e4;
  const gdouble lambda_a[3] = { 20.0, 30.0, 45.0 };
  guint k;

  nc_cluster_richness_projection_prepare (test->crp, test->mu, test->sigma, tau);

  for (k = 0; k < 3; k++)
  {
    const gdouble lnl  = log (lambda_a[k]);
    const gdouble x    = (lnl - test->mu) / test->sigma;
    const gdouble f_LN = exp (-0.5 * x * x) / (ncm_c_sqrt_2pi () * test->sigma * lambda_a[k]);

    ncm_assert_cmpdouble_e (nc_cluster_richness_projection_eval (test->crp, lnl), ==, f_LN, 1.0e-3, 0.0);
  }
}

void
test_nc_cluster_richness_projection_invalid_range (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  NCM_TEST_FAIL (nc_cluster_richness_projection_set_lnlambda_range (test->crp, log (10.0), log (1.0)));
}

void
test_nc_cluster_richness_projection_invalid_reltol (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  NCM_TEST_FAIL (nc_cluster_richness_projection_set_reltol (test->crp, 2.0));
}

void
test_nc_cluster_richness_projection_invalid_sigma (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  NCM_TEST_FAIL (nc_cluster_richness_projection_prepare (test->crp, test->mu, -1.0, test->tau));
}

void
test_nc_cluster_richness_projection_out_of_range (TestNcClusterRichnessProjection *test, gconstpointer pdata)
{
  nc_cluster_richness_projection_prepare (test->crp, test->mu, test->sigma, test->tau);
  NCM_TEST_FAIL (nc_cluster_richness_projection_eval (test->crp, TEST_LNL_MAX + 1.0));
}


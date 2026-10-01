/***************************************************************************
 *            test_ncm_pln1d.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_pln1d.c
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

typedef struct _TestPLN1DArgs
{
  gdouble R;
  gdouble mu;
  gdouble sigma;
} TestPLN1DArgs;

/* Integrand of P(R | mu, sigma) in x = ln(lambda) */
static gdouble
_test_pln1d_integrand (gdouble x, gpointer userdata)
{
  TestPLN1DArgs *args = (TestPLN1DArgs *) userdata;
  const gdouble dx    = (x - args->mu) / args->sigma;

  return exp (args->R * x - exp (x) - lgamma (args->R + 1.0) - 0.5 * dx * dx) / (args->sigma * sqrt (2.0 * M_PI));
}

static gdouble
_test_pln1d_reference (gdouble R, gdouble mu, gdouble sigma)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
  TestPLN1DArgs args           = {R, mu, sigma};
  gsl_function F               = {&_test_pln1d_integrand, &args};
  gdouble result, abserr;

  gsl_integration_qagi (&F, 0.0, 1.0e-12, 1000, w, &result, &abserr);
  gsl_integration_workspace_free (w);

  return result;
}

static void
test_ncm_pln1d_order (void)
{
  NcmPLN1D *pln = g_object_new (NCM_TYPE_PLN1D, NULL);
  guint gh_order;

  g_assert_cmpuint (ncm_pln1d_get_order (pln), ==, 60);

  {
    GParamSpec *pspec = g_object_class_find_property (G_OBJECT_GET_CLASS (pln), "gh-order");

    g_assert_cmpuint (G_PARAM_SPEC_UINT (pspec)->minimum, ==, 1);
  }

  ncm_pln1d_free (pln);

  pln = ncm_pln1d_new (40);
  g_assert_cmpuint (ncm_pln1d_get_order (pln), ==, 40);

  ncm_pln1d_set_order (pln, 80);
  g_assert_cmpuint (ncm_pln1d_get_order (pln), ==, 80);

  g_object_get (pln, "gh-order", &gh_order, NULL);
  g_assert_cmpuint (gh_order, ==, 80);

  ncm_pln1d_clear (&pln);
  g_assert_null (pln);
}

static void
test_ncm_pln1d_mode (void)
{
  const gdouble R[]     = {0.0, 3.0, 20.0, 50.0, 1000.0};
  const gdouble mu[]    = {-1.0, 0.5, 3.0};
  const gdouble sigma[] = {1.0e-3, 0.3, 1.5};
  guint i, j, k;

  /*
   * At the mode, sigma e^{sigma z} = u - z with u = mu / sigma + R sigma. The subtraction
   * z = u - W / sigma loses a factor sigma u of relative precision.
   */
  for (i = 0; i < G_N_ELEMENTS (R); i++)
    for (j = 0; j < G_N_ELEMENTS (mu); j++)
      for (k = 0; k < G_N_ELEMENTS (sigma); k++)
      {
        const gdouble u = mu[j] / sigma[k] + R[i] * sigma[k];
        const gdouble z = ncm_pln1d_mode (R[i], mu[j], sigma[k]);

        ncm_assert_cmpdouble_e (sigma[k] * exp (sigma[k] * z), ==, u - z, 1.0e-15 * (1.0 + sigma[k] * fabs (u)), 1.0e-13 * (1.0 + fabs (u) + fabs (z)));
      }
}

static void
test_ncm_pln1d_eval_quadrature (void)
{
  NcmPLN1D *pln = ncm_pln1d_new (60);

  /* {R, mu, sigma} */
  const gdouble points[][3] = {
    {0.0, 1.0, 0.5},
    {1.0, 0.0, 0.1},
    {5.0, 1.5, 0.3},
    {5.0, 3.0, 0.8},
    {20.0, 3.0, 0.3},
    {20.0, 3.0, 0.8},
    {50.0, 3.0, 1.5},
    {2.0, 0.0, 0.5},
    {5.0, 1.0, 0.8},
    {1.0, -1.0, 1.0},
    {10.0, 2.0, 0.3},
    {0.1, -0.5, 0.8},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (points); i++)
  {
    const gdouble R     = points[i][0];
    const gdouble mu    = points[i][1];
    const gdouble sigma = points[i][2];
    const gdouble ref   = _test_pln1d_reference (R, mu, sigma);

    ncm_assert_cmpdouble_e (ncm_pln1d_eval_p (pln, R, mu, sigma), ==, ref, 1.0e-10, 0.0);
    ncm_assert_cmpdouble_e (ncm_pln1d_eval_lnp (pln, R, mu, sigma), ==, log (ref), 0.0, 1.0e-10);
  }

  ncm_pln1d_free (pln);
}

static void
test_ncm_pln1d_eval_laplace (void)
{
  NcmPLN1D *pln       = ncm_pln1d_new (60);
  const gdouble sigma = 1.0e-5;
  const gdouble R[]   = {0.0, 5.0, 20.0};
  const gdouble mu    = 1.5;
  guint i;

  /* Laplace approximation: tends to the Poisson probability with lambda = e^mu */
  for (i = 0; i < G_N_ELEMENTS (R); i++)
  {
    const gdouble lambda   = exp (mu);
    const gdouble lnp_pois = R[i] * mu - lambda - lgamma (R[i] + 1.0);

    ncm_assert_cmpdouble_e (ncm_pln1d_eval_lnp (pln, R[i], mu, sigma), ==, lnp_pois, 0.0, 1.0e-5);
  }

  ncm_pln1d_free (pln);
}

static gdouble
_test_pln1d_sum_eval_p (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma)
{
  gdouble sum = 0.0;
  guint R;

  for (R = R_min; R <= R_max; R++)
    sum += ncm_pln1d_eval_p (pln, R, mu, sigma);

  return sum;
}

static void
test_ncm_pln1d_eval_normalization (void)
{
  NcmPLN1D *pln = ncm_pln1d_new (60);

  /* The probabilities of all counts sum to one; the lognormal tail above R = 999 is 3e-12 */
  ncm_assert_cmpdouble_e (_test_pln1d_sum_eval_p (pln, 0, 999, 0.0, 1.0), ==, 1.0, 0.0, 1.0e-10);
  ncm_assert_cmpdouble_e (ncm_pln1d_eval_range_sum (pln, 0, 999, 0.0, 1.0), ==, 1.0, 0.0, 1.0e-10);

  ncm_pln1d_free (pln);
}

static void
test_ncm_pln1d_eval_range_sum (void)
{
  NcmPLN1D *pln     = ncm_pln1d_new (60);
  const gdouble sum = _test_pln1d_sum_eval_p (pln, 0, 10, 1.5, 0.3);

  ncm_assert_cmpdouble_e (ncm_pln1d_eval_range_sum (pln, 0, 10, 1.5, 0.3), ==, sum, 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (ncm_pln1d_eval_range_sum_lnp (pln, 0, 10, 1.5, 0.3), ==, log (sum), 1.0e-12, 0.0);

  ncm_assert_cmpdouble_e (ncm_pln1d_eval_range_sum (pln, 0, 20, 0.0, 1.0), ==, _test_pln1d_sum_eval_p (pln, 0, 20, 0.0, 1.0), 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (ncm_pln1d_eval_range_sum (pln, 5, 25, 3.0, 0.5), ==, _test_pln1d_sum_eval_p (pln, 5, 25, 3.0, 0.5), 1.0e-13, 0.0);

  /* Counts with sigma^2 e^{sigma u} beyond the double range */
  g_assert_true (gsl_finite (ncm_pln1d_eval_range_sum_lnp (pln, 300, 2000, 3.0, 1.5)));

  /* A single count reduces to ncm_pln1d_eval_p() */
  ncm_assert_cmpdouble_e (ncm_pln1d_eval_range_sum (pln, 4, 4, 1.5, 0.3), ==, ncm_pln1d_eval_p (pln, 4.0, 1.5, 0.3), 1.0e-14, 0.0);

  ncm_pln1d_free (pln);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/pln1d/order", &test_ncm_pln1d_order);
  g_test_add_func ("/ncm/pln1d/mode", &test_ncm_pln1d_mode);
  g_test_add_func ("/ncm/pln1d/eval/quadrature", &test_ncm_pln1d_eval_quadrature);
  g_test_add_func ("/ncm/pln1d/eval/laplace", &test_ncm_pln1d_eval_laplace);
  g_test_add_func ("/ncm/pln1d/eval/normalization", &test_ncm_pln1d_eval_normalization);
  g_test_add_func ("/ncm/pln1d/eval/range_sum", &test_ncm_pln1d_eval_range_sum);

  g_test_run ();
}


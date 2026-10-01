/***************************************************************************
 *            test_ncm_prior.c
 *
 *  Wed September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_prior.c
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

typedef struct _TestNcmPrior
{
  NcmModel *model;
  NcmMSet *mset;
} TestNcmPrior;

static void
test_ncm_prior_new (TestNcmPrior *test, gconstpointer pdata)
{
  test->model = NCM_MODEL (ncm_model_mvnd_new (1));
  test->mset  = ncm_mset_new (test->model, NULL, NULL);
}

static void
test_ncm_prior_free (TestNcmPrior *test, gconstpointer pdata)
{
  ncm_mset_free (test->mset);
  ncm_model_free (test->model);
}

static gdouble
_test_prior_at (TestNcmPrior *test, NcmPrior *prior, const gdouble x)
{
  ncm_model_orig_param_set (test->model, 0, x);

  return ncm_mset_func_eval0 (NCM_MSET_FUNC (prior), test->mset);
}

static void
test_ncm_prior_gauss (TestNcmPrior *test, gconstpointer pdata)
{
  const gdouble mu      = 0.3;
  const gdouble sigma   = 0.05;
  NcmPriorGaussParam *p = ncm_prior_gauss_param_new (test->model, 0, mu, sigma);
  guint i;

  g_assert_false (ncm_prior_is_m2lnL (NCM_PRIOR (p)));

  /* f = (x - mu) / sigma, so that -2 ln P = f^2 = (x - mu)^2 / sigma^2. */
  for (i = 0; i < 11; i++)
  {
    const gdouble x = -0.2 + 0.1 * i;

    ncm_assert_cmpdouble_e (_test_prior_at (test, NCM_PRIOR (p), x), ==, (x - mu) / sigma, 1.0e-14, 1.0e-14);
  }

  ncm_prior_gauss_param_free (p);
}

static void
test_ncm_prior_flat (TestNcmPrior *test, gconstpointer pdata)
{
  const gdouble x0     = -1.0;
  const gdouble x1     = 2.0;
  const gdouble s      = 1.0e-2;
  NcmPriorFlatParam *p = ncm_prior_flat_param_new (test->model, 0, x0, x1, s);
  const gdouble h0     = ncm_prior_flat_get_h0 (NCM_PRIOR_FLAT (p));
  guint i;

  g_assert_false (ncm_prior_is_m2lnL (NCM_PRIOR (p)));

  /* -2 ln P = f^2 is e^h0 at a limit, one half a width inside, e^-h0 a width inside.
   * ln f^2 moves by 2 h0 / s = 4000 per unit x, so the rounding of x - x_i sets the
   * precision: 4000 ulp (2) ~ 2e-12, hence 1e-11. */
  ncm_assert_cmpdouble_e (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x0)), ==, exp (h0), 1.0e-11, 0.0);
  ncm_assert_cmpdouble_e (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x0 + 0.5 * s)), ==, 1.0, 1.0e-11, 0.0);
  ncm_assert_cmpdouble_e (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x0 + s)), ==, exp (-h0), 1.0e-11, 0.0);
  ncm_assert_cmpdouble_e (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x1)), ==, exp (h0), 1.0e-11, 0.0);
  ncm_assert_cmpdouble_e (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x1 - 0.5 * s)), ==, 1.0, 1.0e-11, 0.0);
  ncm_assert_cmpdouble_e (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x1 - s)), ==, exp (-h0), 1.0e-11, 0.0);

  /* Inside [x0 + s, x1 - s] it adds less than e^-h0. */
  for (i = 0; i <= 100; i++)
  {
    const gdouble x = (x0 + s) + (x1 - x0 - 2.0 * s) * i / 100.0;

    g_assert_cmpfloat (gsl_pow_2 (_test_prior_at (test, NCM_PRIOR (p), x)), <=, exp (-h0) * (1.0 + 1.0e-12));
  }

  ncm_prior_flat_param_free (p);
}

static void
test_ncm_prior_gauss_mean_subprocess (void)
{
  TestNcmPrior test;
  NcmPriorGauss *p;

  test_ncm_prior_new (&test, NULL);
  p = g_object_new (NCM_TYPE_PRIOR_GAUSS, NULL);
  ncm_mset_func_eval0 (NCM_MSET_FUNC (p), test.mset);
}

static void
test_ncm_prior_flat_mean_subprocess (void)
{
  TestNcmPrior test;
  NcmPriorFlat *p;

  test_ncm_prior_new (&test, NULL);
  p = g_object_new (NCM_TYPE_PRIOR_FLAT, NULL);
  ncm_mset_func_eval0 (NCM_MSET_FUNC (p), test.mset);
}

static void
test_ncm_prior_traps (void)
{
  g_test_trap_subprocess ("/ncm/prior/gauss/mean/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*method mean not implemented by NcmPriorGauss*");

  g_test_trap_subprocess ("/ncm/prior/flat/mean/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*method mean not implemented by NcmPriorFlat*");
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/prior/gauss", TestNcmPrior, NULL, &test_ncm_prior_new, &test_ncm_prior_gauss, &test_ncm_prior_free);
  g_test_add ("/ncm/prior/flat", TestNcmPrior, NULL, &test_ncm_prior_new, &test_ncm_prior_flat, &test_ncm_prior_free);
  g_test_add_func ("/ncm/prior/traps", &test_ncm_prior_traps);
  g_test_add_func ("/ncm/prior/gauss/mean/subprocess", &test_ncm_prior_gauss_mean_subprocess);
  g_test_add_func ("/ncm/prior/flat/mean/subprocess", &test_ncm_prior_flat_mean_subprocess);

  g_test_run ();
}


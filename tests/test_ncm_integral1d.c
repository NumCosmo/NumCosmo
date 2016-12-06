/***************************************************************************
 *            test_ncm_integral1d.c
 *
 *  Sun February 21 22:44:16 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <sandro@isoftware.com.br>
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

typedef struct _TestNcmIntegral1d
{
  NcmIntegral1d *int1d;
  gdouble reltol;
} TestNcmIntegral1d;

void test_ncm_integral1d_new_sinx (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_new_x5_2_sinx (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_free (TestNcmIntegral1d *test, gconstpointer pdata);

void test_ncm_integral1d_sinx_eval (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_sinx_hermite_p (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_sinx_hermite (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_sinx_laguerre (TestNcmIntegral1d *test, gconstpointer pdata);

void test_ncm_integral1d_x5_2_sinx_eval (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_x5_2_sinx_hermite_p (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_x5_2_sinx_hermite (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_x5_2_sinx_laguerre (TestNcmIntegral1d *test, gconstpointer pdata);

void test_ncm_integral1d_traps (TestNcmIntegral1d *test, gconstpointer pdata);
void test_ncm_integral1d_invalid_test (TestNcmIntegral1d *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/ncm/integral1d/sinx/eval", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_sinx, 
              &test_ncm_integral1d_sinx_eval, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/sinx/hermite_p", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_sinx, 
              &test_ncm_integral1d_sinx_hermite_p, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/sinx/hermite", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_sinx, 
              &test_ncm_integral1d_sinx_hermite, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/sinx/laguerre", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_sinx, 
              &test_ncm_integral1d_sinx_laguerre, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/x5_2_sinx/eval", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_x5_2_sinx, 
              &test_ncm_integral1d_x5_2_sinx_eval, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/x5_2_sinx/hermite_p", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_x5_2_sinx, 
              &test_ncm_integral1d_x5_2_sinx_hermite_p, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/x5_2_sinx/hermite", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_x5_2_sinx, 
              &test_ncm_integral1d_x5_2_sinx_hermite, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/x5_2_sinx/laguerre", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_x5_2_sinx, 
              &test_ncm_integral1d_x5_2_sinx_laguerre, 
              &test_ncm_integral1d_free);

  g_test_add ("/ncm/integral1d/traps", TestNcmIntegral1d, NULL,
              &test_ncm_integral1d_new_sinx,
              &test_ncm_integral1d_traps,
              &test_ncm_integral1d_free);
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_add ("/ncm/integral1d/invalid/test/subprocess", TestNcmIntegral1d, NULL, 
              &test_ncm_integral1d_new_sinx, 
              &test_ncm_integral1d_invalid_test, 
              &test_ncm_integral1d_free);
#endif
  g_test_run ();
}

static gdouble 
test_sin (gpointer userdata, const gdouble x, const gdouble w)
{
  return sin (x);
}

static gdouble 
test_x5_2_sin (gpointer userdata, const gdouble x, const gdouble w)
{
  return x * x * sqrt (fabs (x)) * sin (x);
}

void
test_ncm_integral1d_new_sinx (TestNcmIntegral1d *test, gconstpointer pdata)
{
  test->int1d  = NCM_INTEGRAL1D (ncm_integral1d_ptr_new (&test_sin, NULL));
  test->reltol = ncm_integral1d_get_reltol (test->int1d);

  g_assert_cmpfloat (test->reltol, ==, NCM_INTEGRAL1D_DEFAULT_RELTOL);

  g_assert_cmpfloat (ncm_integral1d_get_abstol (test->int1d), ==, NCM_INTEGRAL1D_DEFAULT_ABSTOL);
  g_assert_cmpuint (ncm_integral1d_get_rule (test->int1d), ==, NCM_INTEGRAL1D_DEFAULT_ALG);
  g_assert_cmpuint (ncm_integral1d_get_partition (test->int1d), ==, NCM_INTEGRAL1D_DEFAULT_PARTITION);

  g_assert (NCM_IS_INTEGRAL1D (test->int1d));
}

void
test_ncm_integral1d_new_x5_2_sinx (TestNcmIntegral1d *test, gconstpointer pdata)
{
  test->int1d  = NCM_INTEGRAL1D (ncm_integral1d_ptr_new (&test_x5_2_sin, NULL));

  test->reltol = ncm_integral1d_get_reltol (test->int1d);
  g_assert_cmpfloat (test->reltol, ==, NCM_INTEGRAL1D_DEFAULT_RELTOL);

  g_assert_cmpfloat (ncm_integral1d_get_abstol (test->int1d), ==, NCM_INTEGRAL1D_DEFAULT_ABSTOL);
  g_assert_cmpuint (ncm_integral1d_get_rule (test->int1d), ==, NCM_INTEGRAL1D_DEFAULT_ALG);
  g_assert_cmpuint (ncm_integral1d_get_partition (test->int1d), ==, NCM_INTEGRAL1D_DEFAULT_PARTITION);
  
  g_assert (NCM_IS_INTEGRAL1D (test->int1d));
}


void _set_destroyed (gpointer b) { gboolean *destroyed = b; *destroyed = TRUE; }

void
test_ncm_integral1d_free (TestNcmIntegral1d *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_integral1d_free, test->int1d);
}

#define NCM_INTEGRAL1D_TESTCMP_ABS g_assert_cmpfloat (fabs (result), <=, prec)
#define NCM_INTEGRAL1D_TESTCMP(d) g_assert_cmpfloat (fabs ((result - (d)) / result), <=, prec)

void
test_ncm_integral1d_sinx_eval (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, prec);
  result = ncm_integral1d_eval (test->int1d, 0.0, 2.0 * ncm_c_pi (), &err);
  NCM_INTEGRAL1D_TESTCMP_ABS;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval (test->int1d, 0.0, 3.0 * ncm_c_pi (), &err);
  NCM_INTEGRAL1D_TESTCMP (2.0);
}

void
test_ncm_integral1d_sinx_hermite_p (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);

  result = ncm_integral1d_eval_gauss_hermite_p (test->int1d, &err);
  NCM_INTEGRAL1D_TESTCMP (0.724778459007076331818227967606L);

  result = ncm_integral1d_eval_gauss_hermite_r_p (test->int1d, 3.0 / 4.0, &err);
  NCM_INTEGRAL1D_TESTCMP (1.01985127282864194182872307604L);
}

void
test_ncm_integral1d_sinx_hermite (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval_gauss_hermite_mur (test->int1d, 3.0 / 4.0, 1.0, &err);
  NCM_INTEGRAL1D_TESTCMP (1.15618751869439848207550580739L);
}

void
test_ncm_integral1d_sinx_laguerre (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval_gauss_laguerre (test->int1d, &err);
  NCM_INTEGRAL1D_TESTCMP (0.5L);
  g_assert_cmpfloat (fabs (result - (0.5L)), <=, prec);

  result = ncm_integral1d_eval_gauss_laguerre_r (test->int1d, 8.0, &err);
  NCM_INTEGRAL1D_TESTCMP (0.0153846153846153846153846153846L);
}

void
test_ncm_integral1d_x5_2_sinx_eval (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval (test->int1d, 0.0, 2.0 * ncm_c_pi (), &err);
  NCM_INTEGRAL1D_TESTCMP (-91.8526176365063540430190011327L);

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval (test->int1d, 0.0, 3.0 * ncm_c_pi (), &err);
  NCM_INTEGRAL1D_TESTCMP (258.801799728660160318221782452L);
}

void
test_ncm_integral1d_x5_2_sinx_hermite_p (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);

  result = ncm_integral1d_eval_gauss_hermite_p (test->int1d, &err);
  NCM_INTEGRAL1D_TESTCMP (1.21497294652280121151075656545L);

  result = ncm_integral1d_eval_gauss_hermite_r_p (test->int1d, 3.0 / 4.0, &err);
  NCM_INTEGRAL1D_TESTCMP (2.15674842608261378569273430710L);
}

void
test_ncm_integral1d_x5_2_sinx_hermite (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval_gauss_hermite_mur (test->int1d, 3.0 / 4.0, 1.0, &err);
  NCM_INTEGRAL1D_TESTCMP (2.21858907956118443081664071147L);  
}

void
test_ncm_integral1d_x5_2_sinx_laguerre (TestNcmIntegral1d *test, gconstpointer pdata)
{
  gdouble err;
  gdouble result;
  const gdouble prec = 1.0e-11;

  ncm_integral1d_set_reltol (test->int1d, prec);
  ncm_integral1d_set_abstol (test->int1d, 0.0);
  result = ncm_integral1d_eval_gauss_laguerre (test->int1d, &err);
  NCM_INTEGRAL1D_TESTCMP (0.378105832435112890110047422239L);

  result = ncm_integral1d_eval_gauss_laguerre_r (test->int1d, 8.0, &err);
  NCM_INTEGRAL1D_TESTCMP (0.000941693635000236963348513154712L);
}

void
test_ncm_integral1d_traps (TestNcmIntegral1d *test, gconstpointer pdata)
{
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_trap_subprocess ("/ncm/integral1d/invalid/test/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_integral1d_invalid_test (TestNcmIntegral1d *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}

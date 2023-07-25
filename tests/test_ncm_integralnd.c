/***************************************************************************
 *            test_ncm_integral_nd.c
 *
 *  Fri July 21 09:44:16 2023
 *  Copyright  2023  Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Eduardo José Barroso 2023 <eduardo.jsbarroso@uel.br>
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

typedef struct _TestNcmIntegralND
{
  NcmIntegralND *intnd;
  gdouble reltol;
} TestNcmIntegralND;

void test_ncm_integral_nd_new_sinx (TestNcmIntegralND *test, gconstpointer pdata);
void test_ncm_integral_nd_free (TestNcmIntegralND *test, gconstpointer pdata);

void test_ncm_integral_nd_sinx_eval (TestNcmIntegralND *test, gconstpointer pdata);

void test_ncm_integral_nd_traps (TestNcmIntegralND *test, gconstpointer pdata);
void test_ncm_integral_nd_invalid_test (TestNcmIntegralND *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/ncm/integralnd/sinx/eval", TestNcmIntegralND, NULL,
              &test_ncm_integral_nd_new_sinx,
              &test_ncm_integral_nd_sinx_eval,
              &test_ncm_integral_nd_free);
  
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/integralnd/invalid/test/subprocess", TestNcmIntegralND, NULL,
              &test_ncm_integral_nd_new_sinx,
              &test_ncm_integral_nd_invalid_test,
              &test_ncm_integral_nd_free);
#endif
  g_test_run ();
}

static void
test_sin (NcmIntegralND * intnd, NcmVector * x, guint dim, guint npoints, guint fdim, NcmVector * fval)
{
   NcmVector * x_d = ncm_vector_new(dim);
   fval = ncm_vector_new(dim);
   ncm_vector_memcpy(x_d, x);
   ncm_vector_set_zero(fval);
}

static void
test_sin_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  guint dim_i, fdim_i;
  dim_i = 1;
  fdim_i = 1;
  fdim = &fdim_i;
  dim = &dim_i;
}

void
test_ncm_integral_nd_new_sinx (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmIntegralNDClass *self;
  self->integrand = &test_sin;
  self->get_dimensions 	  = &test_sin_dim;
  test->intnd  = NCM_INTEGRAL_ND (self);
  test->reltol = ncm_integral_nd_get_reltol (test->intnd);
 
  g_assert_cmpfloat (test->reltol, ==, NCM_INTEGRAL_ND_DEFAULT_RELTOL);
  
  g_assert_cmpfloat (ncm_integral_nd_get_abstol (test->intnd), ==, NCM_INTEGRAL_ND_DEFAULT_ABSTOL);
  
  g_assert_true (NCM_IS_INTEGRAL_ND (test->intnd));
}

void
_set_destroyed (gpointer b)
{
  gboolean *destroyed = b;
  
  *destroyed = TRUE;
}

void
test_ncm_integral_nd_free (TestNcmIntegralND *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_integral_nd_free, test->intnd);
}

#define NCM_INTEGRAL_ND_TESTCMP_ABS g_assert_cmpfloat (fabs (result), <=, prec)
#define NCM_INTEGRAL_ND_TESTCMP(d) g_assert_cmpfloat (fabs ((result - ((gdouble) d)) / result), <=, prec)

void
test_ncm_integral_nd_sinx_eval (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmVector * err = ncm_vector_new(1);
  NcmVector * xi = ncm_vector_new(1);
  NcmVector * xf = ncm_vector_new(1);
  NcmVector * res = ncm_vector_new(1);
  gdouble xi_i = 0.0;
  gdouble xf_i = 2.0 * ncm_c_pi ();
  ncm_vector_set (xf, 0, xf_i);
  ncm_vector_set (xi, 0, xi_i);
  const gdouble prec = 1.0e-11;
  gdouble result;

	
  ncm_integral_nd_set_reltol (test->intnd, prec);
  ncm_integral_nd_set_abstol (test->intnd, prec);
  ncm_integral_nd_eval (test->intnd, xi, xf, res, err);
  result = ncm_vector_get(res, 0);
  NCM_INTEGRAL_ND_TESTCMP_ABS;
  
  ncm_integral_nd_set_reltol (test->intnd, prec);
  ncm_integral_nd_set_abstol (test->intnd, 0.0);
  ncm_integral_nd_eval (test->intnd, xi, xf, res, err);
  result = ncm_vector_get(res, 0);
  NCM_INTEGRAL_ND_TESTCMP (2.0L);
}

void
test_ncm_integral_nd_traps (TestNcmIntegralND *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/integral_nd/invalid/test/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_integral_nd_invalid_test (TestNcmIntegralND *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}


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
void test_ncm_integral_nd_new_x_p_y (TestNcmIntegralND *test, gconstpointer pdata);
void test_ncm_integral_nd_new_acosx_by_p_cz (TestNcmIntegralND *test, gconstpointer pdata);

void test_ncm_integral_nd_free (TestNcmIntegralND *test, gconstpointer pdata);

void test_ncm_integral_nd_sinx_eval (TestNcmIntegralND *test, gconstpointer pdata);
void test_ncm_integral_nd_x_p_y_eval (TestNcmIntegralND *test, gconstpointer pdata);
void test_ncm_integral_nd_acosx_by_p_cz_eval (TestNcmIntegralND *test, gconstpointer pdata);

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

  g_test_add ("/ncm/integralnd/xpy/eval", TestNcmIntegralND, NULL,
              &test_ncm_integral_nd_new_x_p_y,
              &test_ncm_integral_nd_x_p_y_eval,
              &test_ncm_integral_nd_free);

  g_test_add ("/ncm/integralnd/acosx/eval", TestNcmIntegralND, NULL,
              &test_ncm_integral_nd_new_acosx_by_p_cz,
              &test_ncm_integral_nd_acosx_by_p_cz_eval,
              &test_ncm_integral_nd_free);
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/integralnd/invalid/test/subprocess", TestNcmIntegralND, NULL,
              &test_ncm_integral_nd_new_sinx,
              &test_ncm_integral_nd_invalid_test,
              &test_ncm_integral_nd_free);
#endif
  g_test_run ();
}

struct sin_data
{
  gdouble a;
};

struct x_p_y_data
{
  gdouble a;
};

struct acosx_by_p_cz_data
{
  gdouble a;
  gdouble b;
  gdouble c;
};



static void test_sin (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void test_sin_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

static void test_x_p_y (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void test_x_p_y_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

static void test_acosx_by_p_cz (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void test_acosx_by_p_cz_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

NCM_INTEGRAL_ND_DEFINE_TYPE (NCM, TEST_INT_SIN, NcmTestIntSin, ncm_test_int_sin, test_sin_dim, test_sin, struct sin_data)

NCM_INTEGRAL_ND_DEFINE_TYPE (NCM, TEST_INT_XPY, NcmTestIntXPY, ncm_test_int_x_p_y, test_x_p_y_dim, test_x_p_y, struct x_p_y_data)

NCM_INTEGRAL_ND_DEFINE_TYPE (NCM, TEST_INT_ACOSX_BY_P_CZ, NcmTestIntACosx, ncm_test_int_acosx_by_p_cz, test_acosx_by_p_cz_dim, test_acosx_by_p_cz, struct acosx_by_p_cz_data)


static void
test_sin (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcmTestIntSin *test_int_sin = NCM_TEST_INT_SIN (intnd);
  gint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble x_i = ncm_vector_get (x, i);

    ncm_vector_set (fval, i, sin (test_int_sin->data.a * x_i));
    printf ("%d % 22.15g % 22.15g\n", i, x_i, sin (test_int_sin->data.a * x_i));
  }
}

static void
test_sin_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 1;
  *fdim = 1;
}

void
test_ncm_integral_nd_new_sinx (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmTestIntSin *test_int_sin = g_object_new (ncm_test_int_sin_get_type (), NULL);

  test->intnd = NCM_INTEGRAL_ND (test_int_sin);

  test_int_sin->data.a = 3.0;
}


static void
test_x_p_y (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  gint i;
  for (i = 0; i < npoints; i++)
  {
    const gdouble x_i = ncm_vector_get (x, 2*i);
    const gdouble y_i = ncm_vector_get (x, 2*i + 1);
    ncm_vector_set (fval, i, x_i + y_i);
    printf ("%d % 22.15g % 22.15g % 22.15g\n", i, x_i, y_i, x_i + y_i);
  }
}

static void
test_x_p_y_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 2;
  *fdim = 1;
}

static void
test_acosx_by_p_cz (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  gint i;
  NcmTestIntACosx *test_int_acosx = NCM_TEST_INT_ACOSX_BY_P_CZ(intnd);
  for (i = 0; i < npoints; i++)
  {
    const gdouble x_i = ncm_vector_get (x, dim*i);
    const gdouble y_i = ncm_vector_get (x, 3*i + 1);
    const gdouble z_i = ncm_vector_get (x, 3*i + 2);
    ncm_vector_set (fval, i, test_int_acosx->data.a * sin(x_i));
    ncm_vector_set (fval, i+1, (test_int_acosx->data.b * y_i + test_int_acosx->data.c * z_i));
    printf ("%d % 22.15g % 22.15g % 22.15g \n", i, y_i, z_i, test_int_acosx->data.b * y_i + test_int_acosx->data.c * z_i);
  }
}

static void
test_acosx_by_p_cz_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 3;
  *fdim = 2;
}

void
test_ncm_integral_nd_new_x_p_y (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmTestIntXPY *test_int_x_p_y = g_object_new (ncm_test_int_x_p_y_get_type (), NULL);

  test->intnd = NCM_INTEGRAL_ND (test_int_x_p_y);

  test_int_x_p_y->data.a = 1.0;
}

void
test_ncm_integral_nd_new_acosx_by_p_cz (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmTestIntACosx *test_int_acosx = g_object_new (ncm_test_int_acosx_by_p_cz_get_type (), NULL);

  test->intnd = NCM_INTEGRAL_ND (test_int_acosx);

  test_int_acosx->data.a = 1.0;
  test_int_acosx->data.b = 1.0;
  test_int_acosx->data.c = 1.0;
}


void
test_ncm_integral_nd_sinx_eval (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmVector *err     = ncm_vector_new (1);
  NcmVector *xi      = ncm_vector_new (1);
  NcmVector *xf      = ncm_vector_new (1);
  NcmVector *res     = ncm_vector_new (1);
  gdouble xi_i       = 0.0;
  gdouble xf_i       = 2.0 * ncm_c_pi ();
  const gdouble prec = 1.0e-11;
  gdouble result;

  ncm_vector_set (xf, 0, xf_i);
  ncm_vector_set (xi, 0, xi_i);

  ncm_integral_nd_set_reltol (test->intnd, prec);
  ncm_integral_nd_set_abstol (test->intnd, prec);
  ncm_integral_nd_eval (test->intnd, xi, xf, res, err);

  result = ncm_vector_get (res, 0);
  g_assert_cmpfloat (fabs (result), <=, prec);
  ncm_vector_free (err);
  ncm_vector_free (xi);
  ncm_vector_free (xf);
  ncm_vector_free (res);
}
	
	

void
test_ncm_integral_nd_x_p_y_eval (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmVector *err     = ncm_vector_new (1);
  NcmVector *xi      = ncm_vector_new (2);
  NcmVector *xf      = ncm_vector_new (2);
  NcmVector *res     = ncm_vector_new (1);
  gdouble xi_0       = 0.0;
  gdouble xi_1 	     = 1.0;
  gdouble xf_0       = 2.0;
  gdouble xf_1	     = 5.0;
  const gdouble prec = 1.0e-11;
  gdouble result, error;

  ncm_vector_set (xi, 0, xi_0);
  ncm_vector_set (xi, 1, xi_1);
  ncm_vector_set (xf, 0, xf_0);
  ncm_vector_set (xf, 1, xf_1);
  ncm_integral_nd_set_reltol (test->intnd, prec);
  ncm_integral_nd_set_abstol (test->intnd, prec);
  ncm_integral_nd_eval (test->intnd, xi, xf, res, err);
  error = ncm_vector_get (err, 0);
  result = ncm_vector_get (res, 0);
  g_assert_cmpfloat (fabs (error), <=, prec);
  g_assert_cmpfloat (fabs(result - 32.0), <=, ncm_integral_nd_get_abstol (test->intnd));
  ncm_vector_free (err);
  ncm_vector_free (xi);
  ncm_vector_free (xf);
  ncm_vector_free (res);
}


void
test_ncm_integral_nd_acosx_by_p_cz_eval (TestNcmIntegralND *test, gconstpointer pdata)
{
  NcmVector *err     = ncm_vector_new (2);
  NcmVector *xi      = ncm_vector_new (3);
  NcmVector *xf      = ncm_vector_new (3);
  NcmVector *res     = ncm_vector_new (2);
  gdouble xi_0       = 0.0;
  gdouble xi_1       = 1.0;
  gdouble xi_2	     = 1.0;
  gdouble xf_0       = 2.0 * ncm_c_pi ();
  gdouble xf_1       = 2.0;
  gdouble xf_2       = 2.0;

  const gdouble prec = 1.0e-11;
  gdouble result_0, result_1, error_0, error_1;

  ncm_vector_set (xi, 0, xi_0);
  ncm_vector_set (xi, 1, xi_1);
  ncm_vector_set (xi, 2, xi_2);
  ncm_vector_set (xf, 0, xf_0);
  ncm_vector_set (xf, 1, xf_1);
  ncm_vector_set (xf, 2, xf_2);
  
  ncm_integral_nd_set_reltol (test->intnd, prec);
  ncm_integral_nd_set_abstol (test->intnd, prec);
  ncm_integral_nd_eval (test->intnd, xi, xf, res, err);
  
  error_0 = ncm_vector_get (err, 0);
  result_0 = ncm_vector_get (res, 0);
  
  error_1 = ncm_vector_get (err, 1);
  result_1 = ncm_vector_get (res, 1);
  
  g_assert_cmpfloat (fabs (error_0), <=, prec);
  g_assert_cmpfloat (fabs (error_1), <=, prec);

  g_assert_cmpfloat (fabs (result_0 - 0.0), <=, prec);
  g_assert_cmpfloat (fabs (result_1 / (2.0 * ncm_c_pi()) - 3.0), <=, prec);

  ncm_vector_free (err);
  ncm_vector_free (xi);
  ncm_vector_free (xf);
  ncm_vector_free (res);
}

void
test_ncm_integral_nd_free (TestNcmIntegralND *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_integral_nd_free, test->intnd);
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


/***************************************************************************
 *            test_ncm_ode_spline.c
 *
 *  Mon May 19 10:00:00 2026
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

#include <glib.h>

typedef struct _TestNcmOdeSpline
{
  NcmOdeSpline *os;
  NcmSpline *s;
  gdouble alpha;
} TestNcmOdeSpline;

static gdouble
_test_ode_dydx (gdouble y, gdouble x, gpointer userdata)
{
  TestNcmOdeSpline *test = (TestNcmOdeSpline *) userdata;

  return test->alpha * y;
}

void
test_ncm_ode_spline_new (TestNcmOdeSpline *test, gconstpointer pdata)
{
  test->alpha = 2.5;
  test->s     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  test->os    = ncm_ode_spline_new (test->s, _test_ode_dydx);

  ncm_ode_spline_set_reltol (test->os, 1.0e-10);

  g_assert_true (NCM_IS_ODE_SPLINE (test->os));
  g_assert_true (NCM_IS_SPLINE (test->s));
}

void
test_ncm_ode_spline_free (TestNcmOdeSpline *test, gconstpointer pdata)
{
  ncm_ode_spline_free (test->os);
  ncm_spline_free (test->s);
}

void
test_ncm_ode_spline_min_subdivisions_default (TestNcmOdeSpline *test, gconstpointer pdata)
{
  guint min_subdivisions = ncm_ode_spline_get_min_subdivisions (test->os);

  g_assert_cmpuint (min_subdivisions, ==, 0);
}

void
test_ncm_ode_spline_min_subdivisions_set_get (TestNcmOdeSpline *test, gconstpointer pdata)
{
  guint test_values[] = {0, 5, 10, 20, 100};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_values); i++)
  {
    ncm_ode_spline_set_min_subdivisions (test->os, test_values[i]);
    g_assert_cmpuint (ncm_ode_spline_get_min_subdivisions (test->os), ==, test_values[i]);
  }
}

void
test_ncm_ode_spline_min_subdivisions_effect (TestNcmOdeSpline *test, gconstpointer pdata)
{
  const gdouble x_0 = 0.0;
  const gdouble y_0 = 1.0;
  const gdouble x_1 = 5.0;
  guint min_subdivisions;
  NcmSpline *ss;
  NcmVector *xv;

  /* Test with min_subdivisions = 0 (default behavior) */
  ncm_ode_spline_set_min_subdivisions (test->os, 0);
  ncm_ode_spline_set_interval (test->os, y_0, x_0, x_1);
  ncm_ode_spline_prepare (test->os, test);

  ss = ncm_ode_spline_peek_spline (test->os);
  xv = ncm_spline_peek_xv (ss);

  /* With default behavior, we don't enforce a minimum */
  g_assert_true (ncm_vector_len (xv) >= 2);

  /* Test with min_subdivisions = 10 */
  min_subdivisions = 10;
  ncm_ode_spline_set_min_subdivisions (test->os, min_subdivisions);
  ncm_ode_spline_set_interval (test->os, y_0, x_0, x_1);
  ncm_ode_spline_prepare (test->os, test);

  ss = ncm_ode_spline_peek_spline (test->os);
  xv = ncm_spline_peek_xv (ss);

  /* With min_subdivisions, we should have at least that many points */
  g_assert_cmpuint (ncm_vector_len (xv), >=, min_subdivisions);

  /* Test with min_subdivisions = 50 */
  min_subdivisions = 50;
  ncm_ode_spline_set_min_subdivisions (test->os, min_subdivisions);
  ncm_ode_spline_set_interval (test->os, y_0, x_0, x_1);
  ncm_ode_spline_prepare (test->os, test);

  ss = ncm_ode_spline_peek_spline (test->os);
  xv = ncm_spline_peek_xv (ss);

  g_assert_cmpuint (ncm_vector_len (xv), >=, min_subdivisions);
}

void
test_ncm_ode_spline_min_subdivisions_prop (TestNcmOdeSpline *test, gconstpointer pdata)
{
  guint test_value = 42;
  guint retrieved_value;

  /* Test setting via GObject property */
  g_object_set (test->os, "min-subdivisions", test_value, NULL);
  g_object_get (test->os, "min-subdivisions", &retrieved_value, NULL);

  g_assert_cmpuint (retrieved_value, ==, test_value);
  g_assert_cmpuint (ncm_ode_spline_get_min_subdivisions (test->os), ==, test_value);
}

/* dy/dx = alpha y from y(0) = 1 is exp (alpha x); bound 10 times the measured 4.8e-8. */
void
test_ncm_ode_spline_values (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmSpline *ss;
  guint i;

  ncm_ode_spline_set_interval (test->os, 1.0, 0.0, 5.0);
  ncm_ode_spline_prepare (test->os, test);
  ss = ncm_ode_spline_peek_spline (test->os);

  for (i = 0; i <= 1000; i++)
  {
    const gdouble x = 5.0 * i / 1000.0;

    ncm_assert_cmpdouble_e (ncm_spline_eval (ss, x), ==, exp (test->alpha * x), 5.0e-7, 0.0);
  }
}

/* With yf set, the integration stops where y reaches it, x = ln (yf) / alpha; bound 10
 * times the measured 3.8e-9. */
void
test_ncm_ode_spline_yf (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmVector *xv;

  ncm_ode_spline_set_xi (test->os, 0.0);
  ncm_ode_spline_set_yi (test->os, 1.0);
  ncm_ode_spline_set_yf (test->os, exp (10.0));
  ncm_ode_spline_prepare (test->os, test);

  xv = ncm_spline_peek_xv (ncm_ode_spline_peek_spline (test->os));

  ncm_assert_cmpdouble_e (ncm_vector_get (xv, ncm_vector_len (xv) - 1), ==, 10.0 / test->alpha, 4.0e-8, 0.0);
  ncm_assert_cmpdouble_e (ncm_ode_spline_get_yf_attained (test->os), ==, exp (10.0), 1.0e-8, 0.0);
}

/* The yf stopping mode integrates toward increasing x also from a negative start. */
void
test_ncm_ode_spline_yf_negative (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmVector *xv;

  ncm_ode_spline_set_xi (test->os, -1.0);
  ncm_ode_spline_set_yi (test->os, 1.0);
  ncm_ode_spline_set_yf (test->os, exp (10.0));
  ncm_ode_spline_prepare (test->os, test);

  xv = ncm_spline_peek_xv (ncm_ode_spline_peek_spline (test->os));

  ncm_assert_cmpdouble_e (ncm_vector_get (xv, ncm_vector_len (xv) - 1), ==, -1.0 + 10.0 / test->alpha, 4.0e-8, 0.0);
}

/* The last knot is xf, also with a step limit that makes the final step short. */
void
test_ncm_ode_spline_ends_at_xf (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmVector *xv;

  ncm_ode_spline_set_min_subdivisions (test->os, 7);
  ncm_ode_spline_set_interval (test->os, 1.0, 0.0, 5.0);
  ncm_ode_spline_prepare (test->os, test);

  xv = ncm_spline_peek_xv (ncm_ode_spline_peek_spline (test->os));

  g_assert_cmpfloat (ncm_vector_get (xv, ncm_vector_len (xv) - 1), ==, 5.0);
}

static gdouble
_test_ode_d2ydx2 (gdouble y, gdouble x, gpointer userdata)
{
  TestNcmOdeSpline *test = (TestNcmOdeSpline *) userdata;

  return test->alpha * test->alpha * y;
}

/* With a NcmSplineHermite the knots carry y' and y'': from the callback, or from CVODE
 * with the difference quotient at the first knots. Measured at reltol 1e-10: 198 knots
 * against 2782 with the cubic, both sources of y'' give exp (alpha x) to 6.0e-10 (the
 * cubic to 4.8e-8), and CVODE's y'' is within 9.9e-5 of the exact one; the bounds are 10
 * times these. */
static void
_test_ode_spline_hermite (TestNcmOdeSpline *test, NcmOdeSplineD2ydx2 d2ydx2, const gdouble d2_bound)
{
  NcmOdeSpline *os = ncm_ode_spline_new_hermite (_test_ode_dydx, d2ydx2, 1.0, 0.0, 5.0);
  NcmSpline *ss    = ncm_ode_spline_peek_spline (os);
  NcmVector *xv, *d1v, *d2v;
  guint i;

  test->alpha = 2.5;
  ncm_ode_spline_set_reltol (os, 1.0e-10);
  ncm_ode_spline_prepare (os, test);

  g_assert_true (NCM_IS_SPLINE_HERMITE (ss));
  xv  = ncm_spline_peek_xv (ss);
  d1v = ncm_spline_hermite_peek_d1v (NCM_SPLINE_HERMITE (ss));
  d2v = ncm_spline_hermite_peek_d2v (NCM_SPLINE_HERMITE (ss));
  g_assert_cmpuint (ncm_vector_len (d1v), ==, ncm_vector_len (xv));
  g_assert_cmpuint (ncm_vector_len (d2v), ==, ncm_vector_len (xv));

  for (i = 0; i < ncm_vector_len (xv); i++)
  {
    const gdouble y = ncm_spline_eval (ss, ncm_vector_get (xv, i));

    g_assert_cmpfloat (ncm_vector_get (d1v, i), ==, test->alpha * y);
    g_assert_true (gsl_finite (ncm_vector_get (d2v, i)));
    ncm_assert_cmpdouble_e (ncm_vector_get (d2v, i), ==, test->alpha * test->alpha * y, d2_bound, 0.0);
  }

  for (i = 0; i <= 1000; i++)
  {
    const gdouble x = 5.0 * i / 1000.0;

    ncm_assert_cmpdouble_e (ncm_spline_eval (ss, x), ==, exp (test->alpha * x), 6.0e-9, 0.0);
  }

  g_assert_cmpuint (ncm_vector_len (xv), <, 400);

  ncm_ode_spline_free (os);
}

void
test_ncm_ode_spline_hermite_d2ydx2 (TestNcmOdeSpline *test, gconstpointer pdata)
{
  _test_ode_spline_hermite (test, _test_ode_d2ydx2, 0.0);
}

void
test_ncm_ode_spline_hermite_dky (TestNcmOdeSpline *test, gconstpointer pdata)
{
  _test_ode_spline_hermite (test, NULL, 1.0e-3);
}

/* The yf stop and the knot replacement at xf also carry the derivatives */
void
test_ncm_ode_spline_hermite_yf_and_xf (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmOdeSpline *os = ncm_ode_spline_new_hermite (_test_ode_dydx, _test_ode_d2ydx2, 1.0, 0.0, 5.0);
  NcmSpline *ss    = ncm_ode_spline_peek_spline (os);
  NcmVector *xv;

  test->alpha = 2.5;
  ncm_ode_spline_set_reltol (os, 1.0e-10);
  ncm_ode_spline_set_min_subdivisions (os, 7);
  ncm_ode_spline_prepare (os, test);
  xv = ncm_spline_peek_xv (ss);
  g_assert_cmpfloat (ncm_vector_get (xv, ncm_vector_len (xv) - 1), ==, 5.0);
  g_assert_cmpfloat (ncm_vector_get (ncm_spline_hermite_peek_d1v (NCM_SPLINE_HERMITE (ss)), ncm_vector_len (xv) - 1), ==,
                     test->alpha * ncm_vector_get (ncm_spline_peek_yv (ss), ncm_vector_len (xv) - 1));

  ncm_ode_spline_set_xf (os, GSL_NAN);
  ncm_ode_spline_set_yf (os, exp (10.0));
  ncm_ode_spline_prepare (os, test);
  xv = ncm_spline_peek_xv (ss);
  ncm_assert_cmpdouble_e (ncm_vector_get (xv, ncm_vector_len (xv) - 1), ==, 10.0 / test->alpha, 1.0e-7, 0.0);
  g_assert_cmpuint (ncm_vector_len (ncm_spline_hermite_peek_d2v (NCM_SPLINE_HERMITE (ss))), ==, ncm_vector_len (xv));

  ncm_ode_spline_free (os);
}

/* An xf within NCM_ODE_SPLINE_MIN_STEP of the knot before it: the last step replaces that
 * knot, with its derivatives, and the spline ends at xf with one knot fewer. A fixed initial
 * step keeps the knots independent of xf (CVODE's own initial step depends on it). */
void
test_ncm_ode_spline_hermite_replace_last (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmOdeSpline *os = ncm_ode_spline_new_hermite (_test_ode_dydx, _test_ode_d2ydx2, 1.0, 0.0, 5.0);
  NcmSpline *ss    = ncm_ode_spline_peek_spline (os);
  NcmVector *xv;
  gdouble x_prev, xf;
  guint n;

  test->alpha = 2.5;
  ncm_ode_spline_set_reltol (os, 1.0e-10);
  ncm_ode_spline_set_ini_step (os, 1.0e-3);
  ncm_ode_spline_prepare (os, test);
  xv     = ncm_spline_peek_xv (ss);
  n      = ncm_vector_len (xv);
  x_prev = ncm_vector_get (xv, n - 2);
  xf     = x_prev * (1.0 + 0.3 * NCM_ODE_SPLINE_MIN_STEP);

  ncm_ode_spline_set_xf (os, xf);
  ncm_ode_spline_prepare (os, test);
  xv = ncm_spline_peek_xv (ss);
  g_assert_cmpuint (ncm_vector_len (xv), ==, n - 1);
  g_assert_cmpfloat (ncm_vector_get (xv, n - 2), ==, xf);
  g_assert_cmpfloat (ncm_vector_get (ncm_spline_hermite_peek_d1v (NCM_SPLINE_HERMITE (ss)), n - 2), ==,
                     test->alpha * ncm_vector_get (ncm_spline_peek_yv (ss), n - 2));
  ncm_assert_cmpdouble_e (ncm_spline_eval (ss, xf), ==, exp (test->alpha * xf), 6.0e-9, 0.0);

  ncm_ode_spline_free (os);
}

/* max-order: zero selects the spline's default, 3 for a cubic and 7 for a Hermite spline;
 * a higher order takes fewer steps */
void
test_ncm_ode_spline_max_order (TestNcmOdeSpline *test, gconstpointer pdata)
{
  NcmOdeSpline *osh = ncm_ode_spline_new_hermite (_test_ode_dydx, _test_ode_d2ydx2, 1.0, 0.0, 5.0);
  guint n3, n7, retrieved;

  g_assert_cmpuint (ncm_ode_spline_get_max_order (test->os), ==, 0);
  g_assert_cmpuint (ncm_ode_spline_get_max_order (osh), ==, 0);

  ncm_ode_spline_set_reltol (osh, 1.0e-10);
  ncm_ode_spline_prepare (osh, test);
  n7 = ncm_spline_get_len (ncm_ode_spline_peek_spline (osh));

  ncm_ode_spline_set_max_order (osh, NCM_ODE_SPLINE_CUBIC_MAX_ORDER);
  ncm_ode_spline_prepare (osh, test);
  n3 = ncm_spline_get_len (ncm_ode_spline_peek_spline (osh));
  g_assert_cmpuint (n7, <, n3);

  ncm_ode_spline_set_interval (test->os, 1.0, 0.0, 5.0);
  ncm_ode_spline_prepare (test->os, test);
  g_assert_cmpuint (ncm_spline_get_len (ncm_ode_spline_peek_spline (test->os)), ==, n3);

  g_object_set (osh, "max-order", 5, NULL);
  g_object_get (osh, "max-order", &retrieved, NULL);
  g_assert_cmpuint (retrieved, ==, 5);

  {
    gpointer d2ydx2 = NULL;

    g_object_get (osh, "d2ydx2", &d2ydx2, NULL);
    g_assert_true (d2ydx2 == (gpointer) & _test_ode_d2ydx2);
  }
  ncm_ode_spline_prepare (osh, test);
  g_assert_cmpuint (ncm_spline_get_len (ncm_ode_spline_peek_spline (osh)), <, n3);

  ncm_ode_spline_free (osh);
}

static gdouble
_test_ode_dydx_blowup (gdouble y, gdouble x, gpointer userdata)
{
  return (x < 1.0) ? y : GSL_NAN;
}

/* A non-finite right-hand side stalls CVODE at x = 1, which aborts. */
void
test_ncm_ode_spline_failure (void)
{
  g_test_trap_subprocess ("/ncm/ode_spline/failure/subprocess", 60000000, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*step size underflowed*");
}

void
test_ncm_ode_spline_failure_subprocess (void)
{
  NcmSpline *s     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmOdeSpline *os = ncm_ode_spline_new (s, _test_ode_dydx_blowup);

  ncm_ode_spline_set_reltol (os, 1.0e-10);
  ncm_ode_spline_set_interval (os, 1.0, 0.0, 5.0);
  ncm_ode_spline_prepare (os, NULL);
}

/* Without stop-hnil the stalled integration stops, and the spline ends where it stalled,
 * within the knot spacing NCM_ODE_SPLINE_MIN_STEP. */
void
test_ncm_ode_spline_stall_stops (void)
{
  NcmSpline *s     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmOdeSpline *os = ncm_ode_spline_new (s, _test_ode_dydx_blowup);
  NcmVector *xv;

  g_object_set (os, "stop-hnil", FALSE, NULL);
  ncm_ode_spline_set_reltol (os, 1.0e-10);
  ncm_ode_spline_set_interval (os, 1.0, 0.0, 5.0);
  ncm_ode_spline_prepare (os, NULL);

  xv = ncm_spline_peek_xv (ncm_ode_spline_peek_spline (os));
  ncm_assert_cmpdouble_e (ncm_vector_get (xv, ncm_vector_len (xv) - 1), ==, 1.0, NCM_ODE_SPLINE_MIN_STEP, 0.0);

  ncm_ode_spline_free (os);
  ncm_spline_free (s);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/ode_spline/values", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_values,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/yf", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_yf,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/yf_negative", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_yf_negative,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/ends_at_xf", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_ends_at_xf,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/hermite/d2ydx2", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_hermite_d2ydx2,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/hermite/dky", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_hermite_dky,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/hermite/yf_and_xf", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_hermite_yf_and_xf,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/hermite/replace_last", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_hermite_replace_last,
              &test_ncm_ode_spline_free);
  g_test_add ("/ncm/ode_spline/max_order", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_max_order,
              &test_ncm_ode_spline_free);
  g_test_add_func ("/ncm/ode_spline/failure", &test_ncm_ode_spline_failure);
  g_test_add_func ("/ncm/ode_spline/failure/subprocess", &test_ncm_ode_spline_failure_subprocess);
  g_test_add_func ("/ncm/ode_spline/stall_stops", &test_ncm_ode_spline_stall_stops);
  g_test_add ("/ncm/ode_spline/min_subdivisions/default", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_min_subdivisions_default,
              &test_ncm_ode_spline_free);

  g_test_add ("/ncm/ode_spline/min_subdivisions/set_get", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_min_subdivisions_set_get,
              &test_ncm_ode_spline_free);

  g_test_add ("/ncm/ode_spline/min_subdivisions/effect", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_min_subdivisions_effect,
              &test_ncm_ode_spline_free);

  g_test_add ("/ncm/ode_spline/min_subdivisions/prop", TestNcmOdeSpline, NULL,
              &test_ncm_ode_spline_new,
              &test_ncm_ode_spline_min_subdivisions_prop,
              &test_ncm_ode_spline_free);

  g_test_run ();

  return 0;
}


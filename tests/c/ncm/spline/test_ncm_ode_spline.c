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


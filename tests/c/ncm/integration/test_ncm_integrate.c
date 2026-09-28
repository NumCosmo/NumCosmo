/***************************************************************************
 *            test_ncm_integrate.c
 *
 *  Sat September 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_integrate.c
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
#include <gsl/gsl_errno.h>

static gdouble
_f_osc (gdouble x, gpointer p)
{
  if (p != NULL)
    (*(guint *) p)++;

  return exp (sin (5.0 * x));
}

/* A GSL failure aborts instead of returning a result. */
static void
test_ncm_integrate_locked_failure (void)
{
  g_test_trap_subprocess ("/ncm/integrate/locked_a_b/failure/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_integral_locked_a_b*");
}

static void
test_ncm_integrate_locked_failure_subprocess (void)
{
  gsl_function F = {&_f_osc, NULL};
  gdouble res, err;

  /* A relative tolerance GSL rejects, with its own error handler off */
  gsl_set_error_handler_off ();
  ncm_integral_locked_a_b (&F, 0.0, 3.0, 0.0, 1.0e-20, &res, &err);
}

/* The cached integral from zero uses the tolerances of the cache: a loose cache takes
 * fewer evaluations and stays within its tolerance of a tight one. */
static void
test_ncm_integrate_cached_0_x (void)
{
  guint n_loose           = 0;
  guint n_tight           = 0;
  gsl_function F_loose    = {&_f_osc, &n_loose};
  gsl_function F_tight    = {&_f_osc, &n_tight};
  NcmFunctionCache *loose = ncm_function_cache_new (1, 0.0, 1.0e-3);
  NcmFunctionCache *tight = ncm_function_cache_new (1, 0.0, 1.0e-13);
  gdouble res_loose, res_tight, err;

  ncm_integral_cached_0_x (loose, &F_loose, 3.0, &res_loose, &err);
  ncm_integral_cached_0_x (tight, &F_tight, 3.0, &res_tight, &err);

  g_assert_cmpuint (n_loose, <, n_tight);
  ncm_assert_cmpdouble_e (res_loose, ==, res_tight, 1.0e-3, 0.0);

  ncm_function_cache_free (loose);
  ncm_function_cache_free (tight);
}

static gdouble
_f_x_p_y (gdouble x, gdouble y, gpointer p)
{
  return x + y;
}

/* Divonne integrates correctly and leaves the caller's peak list unchanged. */
static void
test_ncm_integrate_divonne_xgiven (void)
{
  NcmIntegrand2dim integ = {NULL, &_f_x_p_y};
  gdouble xgiven[2]      = {0.25, 1.5};
  gdouble res, err;

  ncm_integrate_2dim_divonne (&integ, 0.0, 0.0, 1.0, 2.0, 1.0e-8, 0.0, 1, 2, xgiven, &res, &err);

  g_assert_cmpfloat (xgiven[0], ==, 0.25);
  g_assert_cmpfloat (xgiven[1], ==, 1.5);
  ncm_assert_cmpdouble_e (res, ==, 3.0, 1.0e-8, 0.0);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/integrate/locked_a_b/failure", &test_ncm_integrate_locked_failure);
  g_test_add_func ("/ncm/integrate/locked_a_b/failure/subprocess", &test_ncm_integrate_locked_failure_subprocess);
  g_test_add_func ("/ncm/integrate/cached_0_x", &test_ncm_integrate_cached_0_x);
  g_test_add_func ("/ncm/integrate/divonne_xgiven", &test_ncm_integrate_divonne_xgiven);

  g_test_run ();
}


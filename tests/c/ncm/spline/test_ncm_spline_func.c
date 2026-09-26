/***************************************************************************
 *            test_ncm_spline_func.c
 *
 *  Sat September 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spline_func.c
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

static gdouble
_f (gdouble x, gpointer p)
{
  return exp (-x) * sin (3.0 * x) + 2.0;
}

static gdouble
_step (gdouble x, gpointer p)
{
  return (x < 1.0 / 3.0) ? 0.0 : 1.0;
}

static gdouble
_max_rel_err (NcmSpline *s, const gdouble xi, const gdouble xf)
{
  gdouble err = 0.0;
  guint i;

  for (i = 0; i <= 4000; i++)
  {
    const gdouble x = xi * pow (xf / xi, i / 4000.0);

    err = GSL_MAX (err, fabs (ncm_spline_eval (s, x) / _f (x, NULL) - 1.0));
  }

  return err;
}

/* Every adaptive type reaches the requested relative error away from the knots; the
 * largest ratio measured is 0.45 (sinhknot at 1e-8). */
static void
test_ncm_spline_func_adaptive (void)
{
  const NcmSplineFuncType types[] = {
    NCM_SPLINE_FUNCTION_SPLINE,
    NCM_SPLINE_FUNCTION_SPLINE_LNKNOT,
    NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT,
  };
  const gdouble rels[] = {1.0e-4, 1.0e-6, 1.0e-8};
  gsl_function F       = {&_f, NULL};
  guint t, r;

  for (t = 0; t < G_N_ELEMENTS (types); t++)
  {
    for (r = 0; r < G_N_ELEMENTS (rels); r++)
    {
      NcmSpline *s = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

      ncm_spline_set_func (s, types[t], &F, 0.1, 5.0, 0, rels[r]);

      g_assert_cmpfloat (_max_rel_err (s, 0.1, 5.0), <, rels[r]);

      ncm_spline_free (s);
    }
  }
}

/* The grids have the requested number of knots, uniform in x or in ln x. */
static void
test_ncm_spline_func_grid (void)
{
  gsl_function F = {&_f, NULL};
  NcmSpline *s   = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmVector *xv;
  guint i;

  ncm_spline_set_func_grid (s, NCM_SPLINE_FUNC_GRID_LINEAR, &F, 0.1, 5.0, 50);
  xv = ncm_spline_peek_xv (s);
  g_assert_cmpuint (ncm_vector_len (xv), ==, 50);

  for (i = 0; i < 50; i++)
  {
    ncm_assert_cmpdouble_e (ncm_vector_get (xv, i), ==, 0.1 + 4.9 * i / 49.0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_spline_eval (s, ncm_vector_get (xv, i)), ==, _f (ncm_vector_get (xv, i), NULL), 1.0e-15, 0.0);
  }

  ncm_spline_set_func_grid (s, NCM_SPLINE_FUNC_GRID_LOG, &F, 0.1, 5.0, 50);
  xv = ncm_spline_peek_xv (s);
  g_assert_cmpuint (ncm_vector_len (xv), ==, 50);

  for (i = 0; i < 50; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (xv, i), ==, 0.1 * pow (50.0, i / 49.0), 1.0e-14, 0.0);

  ncm_spline_free (s);
}

/* Exceeding max_nodes stops the refinement with a warning. */
static void
test_ncm_spline_func_max_nodes (void)
{
  g_test_trap_subprocess ("/ncm/spline_func/max_nodes/subprocess", 0, 0);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stderr ("*cannot achieve requested precision*");
}

static void
test_ncm_spline_func_max_nodes_subprocess (void)
{
  gsl_function F = {&_f, NULL};
  NcmSpline *s   = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  g_log_set_always_fatal (G_LOG_FATAL_MASK);
  ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.1, 5.0, 20, 1.0e-8);
  g_assert_cmpuint (ncm_spline_get_len (s), >, 20);

  ncm_spline_free (s);
}

/* A discontinuity drives the knot spacing to NCM_SPLINE_KNOT_DIFF_TOL, which aborts. */
static void
test_ncm_spline_func_discontinuous (void)
{
  g_test_trap_subprocess ("/ncm/spline_func/discontinuous/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*probably discontinuous*");
}

static void
test_ncm_spline_func_discontinuous_subprocess (void)
{
  gsl_function F = {&_step, NULL};
  NcmSpline *s   = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.1, 5.0, 0, 1.0e-6);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/spline_func/adaptive", &test_ncm_spline_func_adaptive);
  g_test_add_func ("/ncm/spline_func/grid", &test_ncm_spline_func_grid);
  g_test_add_func ("/ncm/spline_func/max_nodes", &test_ncm_spline_func_max_nodes);
  g_test_add_func ("/ncm/spline_func/max_nodes/subprocess", &test_ncm_spline_func_max_nodes_subprocess);
  g_test_add_func ("/ncm/spline_func/discontinuous", &test_ncm_spline_func_discontinuous);
  g_test_add_func ("/ncm/spline_func/discontinuous/subprocess", &test_ncm_spline_func_discontinuous_subprocess);

  g_test_run ();
}


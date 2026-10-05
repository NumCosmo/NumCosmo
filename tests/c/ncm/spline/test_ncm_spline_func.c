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

/* Counts the evaluations of the wrapped function. */
typedef struct _TestCountF
{
  gsl_function *F;
  guint n;
} TestCountF;

static gdouble
_count_f (gdouble x, gpointer p)
{
  TestCountF *c = p;

  c->n++;

  return GSL_FN_EVAL (c->F, x);
}

static gdouble
_linear (gdouble x, gpointer p)
{
  return 2.0 * x + 1.0;
}

/*
 * A line, a Gaussian at 0.55 and a pair of opposite Gaussians at 0.65 and 0.75. On [0, 1]
 * the rounds start from the knots 0, 0.2, ..., 1 and test the midpoints 0.1, 0.3, ..., 0.9:
 * the Gaussian is below 1e-24 at 0.4 and 0.6 and 1.9e-3 at 0.5, so [0.4, 0.6] alone fails
 * the first round; the pair is below 1e-17 at 0.6, 0.7 and 0.8, so the test at 0.7 passes
 * on the line and accepts [0.6, 0.7] and [0.7, 0.8] with the pair inside them. The
 * refinement around 0.55 reaches 0.6, the closure flags [0.6, 0.7] against that fine
 * neighbor, and its midpoint 0.65 fails by 1e8 times the tolerance.
 */
static gdouble
_hidden_pair (gdouble x, gpointer p)
{
  return 1.0 + x + exp (-gsl_pow_2 ((x - 0.55) / 0.02))
         + exp (-gsl_pow_2 ((x - 0.65) / 0.008)) - exp (-gsl_pow_2 ((x - 0.75) / 0.008));
}

static void
test_ncm_spline_func_closure_none (void)
{
  /* A cubic spline reproduces a line: every first-round midpoint passes, the mesh stays
   * uniform, and the closure round flags nothing. */
  gsl_function F  = {&_linear, NULL};
  TestCountF cf   = {&F, 0};
  gsl_function Fc = {&_count_f, &cf};
  NcmSpline *s    = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  guint rounds, flagged, failed;

  ncm_spline_set_func_full (s, NCM_SPLINE_FUNCTION_SPLINE, &Fc, 0.1, 5.0, 0, 1.0e-8, 0.0, &rounds, &flagged, &failed);

  g_assert_cmpuint (rounds, ==, 1);
  g_assert_cmpuint (flagged, ==, 0);
  g_assert_cmpuint (failed, ==, 0);
  g_assert_cmpuint (ncm_spline_get_len (s), ==, 11);
  g_assert_cmpuint (cf.n, ==, ncm_spline_get_len (s));

  ncm_spline_free (s);
}

static void
test_ncm_spline_func_closure_pass (void)
{
  /* A smooth function with a varying scale: the settled mesh has level boundaries, so
   * closure rounds flag intervals, and every flagged interval passes. */
  const NcmSplineFuncType types[] = {NCM_SPLINE_FUNCTION_SPLINE, NCM_SPLINE_FUNCTION_SPLINE_LNKNOT, NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT};
  guint t;

  for (t = 0; t < G_N_ELEMENTS (types); t++)
  {
    gsl_function F  = {&_f, NULL};
    TestCountF cf   = {&F, 0};
    gsl_function Fc = {&_count_f, &cf};
    NcmSpline *s    = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    guint rounds, flagged, failed;

    ncm_spline_set_func_full (s, types[t], &Fc, 0.1, 5.0, 0, 1.0e-8, 0.0, &rounds, &flagged, &failed);

    g_assert_cmpuint (rounds, >=, 1);
    g_assert_cmpuint (flagged, >=, 1);
    g_assert_cmpuint (failed, ==, 0);
    g_assert_cmpuint (cf.n, ==, ncm_spline_get_len (s));
    g_assert_cmpfloat (_max_rel_err (s, 0.1, 5.0), <, 1.0e-8);

    ncm_spline_free (s);
  }
}

static void
test_ncm_spline_func_closure_fail (void)
{
  gsl_function F    = {&_hidden_pair, NULL};
  TestCountF cf     = {&F, 0};
  gsl_function Fc   = {&_count_f, &cf};
  NcmSpline *s      = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  const gdouble rel = 1.0e-8;
  gdouble ratio     = 0.0;
  guint rounds, flagged, failed, i;

  ncm_spline_set_func_full (s, NCM_SPLINE_FUNCTION_SPLINE, &Fc, 0.0, 1.0, 0, rel, 1.0, &rounds, &flagged, &failed);

  /* Both halves fail, [0.6, 0.7] first and [0.7, 0.8] once its left neighbor is fine. */
  g_assert_cmpuint (failed, >=, 2);
  g_assert_cmpuint (flagged, >=, failed);
  g_assert_cmpuint (cf.n, ==, ncm_spline_get_len (s));

  for (i = 0; i <= 20000; i++)
  {
    const gdouble x = i / 20000.0;
    const gdouble y = _hidden_pair (x, NULL);

    ratio = GSL_MAX (ratio, fabs (ncm_spline_eval (s, x) - y) / (rel * (fabs (y) + 1.0)));
  }

  g_assert_cmpfloat (ratio, <, 1.0);

  ncm_spline_free (s);
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
  g_test_add_func ("/ncm/spline_func/closure/none", &test_ncm_spline_func_closure_none);
  g_test_add_func ("/ncm/spline_func/closure/pass", &test_ncm_spline_func_closure_pass);
  g_test_add_func ("/ncm/spline_func/closure/fail", &test_ncm_spline_func_closure_fail);

  g_test_run ();
}


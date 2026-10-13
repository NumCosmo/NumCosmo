/***************************************************************************
 *            test_ncm_function_sample_set.c
 *
 *  Wed September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_function_sample_set.c
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

/* F (x) = (x^2, exp (-x)) */
static void
_test_func (const gdouble x, NcmVector *y, gpointer user_data)
{
  ncm_vector_set (y, 0, x * x);
  ncm_vector_set (y, 1, exp (-x));
}

static NcmFunctionSampleSet *
_test_fss_new (const gdouble *x, const guint n)
{
  NcmFunctionSampleSet *fss = ncm_function_sample_set_new (2);
  guint i;

  for (i = 0; i < n; i++)
    ncm_function_sample_set_add_func (fss, x[i], &_test_func, NULL);

  return fss;
}

/* Knots strictly increasing, and the tracked range and maxima those of the samples. */
static void
_test_fss_check_order_and_tracking (NcmFunctionSampleSet *fss)
{
  NcmFunctionSampleSetIter iter_s;
  NcmFunctionSampleSetIter *iter = &iter_s;
  gdouble x_prev                 = GSL_NEGINF;
  gdouble max0                   = 0.0;
  gdouble max1                   = 0.0;
  gdouble x_first                = GSL_NAN;
  gdouble x_last                 = GSL_NAN;

  ncm_function_sample_set_iter_begin (fss, &iter);

  while (ncm_function_sample_set_iter_is_valid (iter))
  {
    const gdouble x = ncm_function_sample_set_iter_get_x (iter);
    NcmVector *y    = ncm_function_sample_set_iter_get_y (iter);

    g_assert_cmpfloat (x, >, x_prev);

    if (!gsl_finite (x_first))
      x_first = x;

    x_last = x;
    x_prev = x;
    max0   = GSL_MAX (max0, fabs (ncm_vector_get (y, 0)));
    max1   = GSL_MAX (max1, fabs (ncm_vector_get (y, 1)));

    ncm_function_sample_set_iter_next (iter);
  }

  g_assert_cmpfloat (ncm_function_sample_set_get_x_min (fss), ==, x_first);
  g_assert_cmpfloat (ncm_function_sample_set_get_x_max (fss), ==, x_last);
  g_assert_cmpfloat (ncm_function_sample_set_get_absmaxF (fss, 0, NULL), ==, max0);
  g_assert_cmpfloat (ncm_function_sample_set_get_absmaxF (fss, 1, NULL), ==, max1);
}

static void
test_ncm_function_sample_set_tracking (void)
{
  const gdouble x[]         = {2.0, 1.0, 3.0};
  NcmFunctionSampleSet *fss = _test_fss_new (x, G_N_ELEMENTS (x));
  NcmFunctionSampleSetIter iter_s, new_s;
  NcmFunctionSampleSetIter *iter = &iter_s;
  NcmFunctionSampleSetIter *new  = &new_s;
  gdouble x_at;

  /* add, insert_after and insert_before all keep the range and the maxima. */
  ncm_function_sample_set_iter_end (fss, &iter);
  ncm_function_sample_set_iter_insert_after_func (fss, iter, 4.0, &_test_func, NULL, &new);
  g_assert_true (ncm_function_sample_set_iter_get_new_point (new));

  ncm_function_sample_set_iter_begin (fss, &iter);
  ncm_function_sample_set_iter_insert_before_func (fss, iter, 0.5, &_test_func, NULL, &new);
  ncm_function_sample_set_add_old_func (fss, 2.5, &_test_func, NULL);

  g_assert_cmpuint (ncm_function_sample_set_get_nsamples (fss), ==, 6);
  _test_fss_check_order_and_tracking (fss);

  g_assert_cmpfloat (ncm_function_sample_set_get_absmaxF (fss, 0, &x_at), ==, 16.0);
  g_assert_cmpfloat (x_at, ==, 4.0);
  g_assert_cmpfloat (ncm_function_sample_set_get_absmaxF (fss, 1, &x_at), ==, exp (-0.5));
  g_assert_cmpfloat (x_at, ==, 0.5);

  /* prev from the first sample, and next from the last, leave the iterator invalid. */
  ncm_function_sample_set_iter_begin (fss, &iter);
  ncm_function_sample_set_iter_prev (iter);
  g_assert_false (ncm_function_sample_set_iter_is_valid (iter));
  ncm_function_sample_set_iter_end (fss, &iter);
  ncm_function_sample_set_iter_next (iter);
  g_assert_false (ncm_function_sample_set_iter_is_valid (iter));

  ncm_function_sample_set_free (fss);
}

static void
test_ncm_function_sample_set_expand_domain (void)
{
  /* exp (-x) falls below 1e-3 of the peak within [1, 1e3]: the right side converges. */
  const gdouble x[]         = {1.0, 1.5, 2.0};
  NcmFunctionSampleSet *fss = _test_fss_new (x, G_N_ELEMENTS (x));

  ncm_function_sample_set_expand_domain (fss, &_test_func, 0.5, 1.0e3, 0.2, 1.0e-3, 200, 3, NULL);

  g_assert_cmpfloat (ncm_function_sample_set_get_x_min (fss), ==, 0.5);
  g_assert_cmpfloat (ncm_function_sample_set_get_x_max (fss), >, 2.0);
  _test_fss_check_order_and_tracking (fss);

  ncm_function_sample_set_free (fss);
}

static void
test_ncm_function_sample_set_expand_domain_at_limits (void)
{
  /* Both sides already at their hard limits: nothing is added (a clamped proposal
   * used to duplicate the endpoint, which the spline then refuses). */
  const gdouble x[]         = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0};
  NcmFunctionSampleSet *fss = _test_fss_new (x, G_N_ELEMENTS (x));
  NcmSpline *s              = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmSplineVec *sv;

  ncm_function_sample_set_expand_domain (fss, &_test_func, 1.0, 10.0, 0.2, 1.0e-3, 10, 3, NULL);

  g_assert_cmpuint (ncm_function_sample_set_get_nsamples (fss), ==, 6);
  _test_fss_check_order_and_tracking (fss);

  sv = ncm_function_sample_set_to_spline_vec (fss, s);
  ncm_spline_vec_free (sv);

  ncm_spline_free (s);
  ncm_function_sample_set_free (fss);
}

static void
test_ncm_function_sample_set_expand_domain_nonpositive_subprocess (void)
{
  const gdouble x[]         = {-1.0, 2.0};
  NcmFunctionSampleSet *fss = _test_fss_new (x, G_N_ELEMENTS (x));

  ncm_function_sample_set_expand_domain (fss, &_test_func, -10.0, 10.0, 0.2, 1.0e-3, 10, 3, NULL);
}

static void
test_ncm_function_sample_set_expand_domain_empty_subprocess (void)
{
  NcmFunctionSampleSet *fss = ncm_function_sample_set_new (2);

  ncm_function_sample_set_expand_domain (fss, &_test_func, 0.1, 10.0, 0.2, 1.0e-3, 10, 3, NULL);
}

static void
test_ncm_function_sample_set_expand_domain_traps (void)
{
  g_test_trap_subprocess ("/ncm/function_sample_set/expand_domain/nonpositive/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the domain * must be positive*");

  g_test_trap_subprocess ("/ncm/function_sample_set/expand_domain/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the sample set is empty*");
}

typedef struct _TestMessageCount
{
  guint max_iter;
  guint at_precision;
} TestMessageCount;

static void
_test_count_messages (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  TestMessageCount *count = user_data;

  if (g_strrstr (message, "Max iterations") != NULL)
    count->max_iter++;

  if (g_strrstr (message, "reached machine precision") != NULL)
    count->at_precision++;
}

static void
test_ncm_function_sample_set_adaptive_midpoint_message (void)
{
  /*
   * The "Max iterations" message is for a refinement that did not converge. Run with
   * one more iteration allowed each time: at the first budget that converges, the loop
   * ends by exhausting it, and no message may be printed.
   */
  NcmSpline *s           = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  TestMessageCount count = {0, 0};
  guint handler_id;
  guint max_iter;

  handler_id = g_log_set_handler ("NUMCOSMO", G_LOG_LEVEL_MESSAGE, &_test_count_messages, &count);

  for (max_iter = 1; max_iter < 60; max_iter++)
  {
    const gdouble x[]         = {0.1, 1.0, 5.0, 10.0};
    NcmFunctionSampleSet *fss = _test_fss_new (x, G_N_ELEMENTS (x));
    gboolean converged;

    count.max_iter = 0;
    ncm_function_sample_set_mark_all_old (fss);
    ncm_function_sample_set_adaptive_midpoint (fss, &_test_func, 1.0e-8, 0.0, max_iter, 1, s, NULL);
    converged = ncm_function_sample_set_all_intervals_ok (fss, 1);
    ncm_function_sample_set_free (fss);

    if (converged)
    {
      g_assert_cmpuint (count.max_iter, ==, 0);
      break;
    }

    g_assert_cmpuint (count.max_iter, ==, 1);
  }

  g_assert_cmpuint (max_iter, <, 60);
  g_assert_cmpuint (count.at_precision, ==, 0);

  g_log_remove_handler ("NUMCOSMO", handler_id);
  ncm_spline_free (s);
}

/* A unit step at x = 0.5 in the first component. */
static void
_test_step (const gdouble x, NcmVector *y, gpointer user_data)
{
  ncm_vector_set (y, 0, (x < 0.5) ? 0.0 : 1.0);
  ncm_vector_set (y, 1, 0.0);
}

static void
test_ncm_function_sample_set_adaptive_midpoint_precision (void)
{
  /* No spline meets the tolerance across a jump: the intervals there shrink to machine
   * precision, are marked as passed, and the count is reported. */
  const gdouble x[]         = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  NcmFunctionSampleSet *fss = ncm_function_sample_set_new (2);
  NcmSpline *s              = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  TestMessageCount count    = {0, 0};
  guint handler_id, i;

  for (i = 0; i < G_N_ELEMENTS (x); i++)
    ncm_function_sample_set_add_func (fss, x[i], &_test_step, NULL);

  ncm_function_sample_set_mark_all_old (fss);

  handler_id = g_log_set_handler ("NUMCOSMO", G_LOG_LEVEL_MESSAGE, &_test_count_messages, &count);
  ncm_function_sample_set_adaptive_midpoint (fss, &_test_step, 1.0e-8, 1.0e-12, 500, 1, s, NULL);
  g_log_remove_handler ("NUMCOSMO", handler_id);

  g_assert_true (ncm_function_sample_set_all_intervals_ok (fss, 1));
  g_assert_cmpuint (count.at_precision, ==, 1);
  g_assert_cmpuint (count.max_iter, ==, 0);

  ncm_spline_free (s);
  ncm_function_sample_set_free (fss);
}

typedef gdouble (*TestEval) (const gdouble x);

static gdouble
_test_line_eval (const gdouble x)
{
  return 2.0 * x + 1.0;
}

static gdouble
_test_lorentzian_eval (const gdouble x)
{
  return 1.0 / (1.0 + x * x);
}

/*
 * A Gaussian at 3.5, seen by the first midpoint test of [3, 4], and a pair of opposite
 * Gaussians at 4.25 and 4.75, which vanish to rounding at 4, 4.5 and 5, the knots and
 * midpoint the first test of [4, 5] uses.
 */
#define TEST_HIDDEN_SIGMA 0.03

static gdouble
_test_hidden_pair_eval (const gdouble x)
{
  const gdouble g0 = exp (-0.5 * gsl_pow_2 ((x - 3.5) / TEST_HIDDEN_SIGMA));
  const gdouble g1 = exp (-0.5 * gsl_pow_2 ((x - 4.25) / TEST_HIDDEN_SIGMA));
  const gdouble g2 = exp (-0.5 * gsl_pow_2 ((x - 4.75) / TEST_HIDDEN_SIGMA));

  return g0 + g1 - g2;
}

static void
_test_closure_f (const gdouble x, NcmVector *y, gpointer user_data)
{
  TestEval eval = (TestEval) user_data;

  ncm_vector_set (y, 0, eval (x));
}

/*
 * Refines from seeds 0, 1, ..., 10, reports the closure, and returns the largest
 * |spline - F| / (reltol |F| + abstol) on a grid over the seed range.
 */
static gdouble
_test_closure (TestEval eval, const gdouble reltol, const gdouble abstol, guint *rounds, guint *flagged, guint *failed, guint *nsamples)
{
  const guint n_grid        = 20001;
  NcmFunctionSampleSet *fss = ncm_function_sample_set_new (1);
  NcmSpline *s              = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmVector *y              = ncm_vector_new (1);
  gdouble max_ratio         = 0.0;
  NcmSplineVec *sv;
  guint i;

  for (i = 0; i <= 10; i++)
    ncm_function_sample_set_add_func (fss, i, &_test_closure_f, eval);

  ncm_function_sample_set_mark_all_old (fss);
  ncm_function_sample_set_adaptive_midpoint_full (fss, &_test_closure_f, reltol, abstol, 1000, 1, s, eval,
                                                  rounds, flagged, failed);
  g_assert_true (ncm_function_sample_set_all_intervals_ok (fss, 1));
  *nsamples = ncm_function_sample_set_get_nsamples (fss);

  sv = ncm_function_sample_set_to_spline_vec (fss, s);

  for (i = 0; i < n_grid; i++)
  {
    const gdouble x  = 10.0 * i / (n_grid - 1.0);
    const gdouble fx = eval (x);

    ncm_spline_vec_eval (sv, x, y);
    max_ratio = GSL_MAX (max_ratio, fabs (ncm_vector_get (y, 0) - fx) / (reltol * fabs (fx) + abstol));
  }

  ncm_spline_vec_free (sv);
  ncm_vector_free (y);
  ncm_spline_free (s);
  ncm_function_sample_set_free (fss);

  return max_ratio;
}

static void
test_ncm_function_sample_set_adaptive_midpoint_closure_none (void)
{
  /* A cubic spline reproduces a line: every first-round midpoint passes, the set stays
   * uniform, and the closure round flags nothing. */
  guint rounds, flagged, failed, nsamples;

  g_assert_cmpfloat (_test_closure (&_test_line_eval, 1.0e-8, 0.0, &rounds, &flagged, &failed, &nsamples), <, 1.0);
  g_assert_cmpuint (rounds, ==, 1);
  g_assert_cmpuint (flagged, ==, 0);
  g_assert_cmpuint (failed, ==, 0);
  g_assert_cmpuint (nsamples, ==, 21);
}

static void
test_ncm_function_sample_set_adaptive_midpoint_closure_pass (void)
{
  /* A smooth function with a varying scale: the settled set has level boundaries, so
   * the closure round flags intervals, and every flagged interval passes. */
  guint rounds, flagged, failed, nsamples;

  g_assert_cmpfloat (_test_closure (&_test_lorentzian_eval, 1.0e-6, 0.0, &rounds, &flagged, &failed, &nsamples), <, 1.0);
  g_assert_cmpuint (rounds, ==, 1);
  g_assert_cmpuint (flagged, >=, 1);
  g_assert_cmpuint (failed, ==, 0);
}

static void
test_ncm_function_sample_set_adaptive_midpoint_closure_fail (void)
{
  /* The midpoint test at 4.5 accepts [4, 4.5] and [4.5, 5] with the pair inside them;
   * the closure flags [4, 4.5] against its refined left neighbor, its midpoint fails,
   * the failure flags [4.5, 5], and both are resolved to the requested tolerance. */
  guint rounds, flagged, failed, nsamples;

  g_assert_cmpfloat (_test_closure (&_test_hidden_pair_eval, 1.0e-6, 1.0e-8, &rounds, &flagged, &failed, &nsamples), <, 1.0);
  g_assert_cmpuint (rounds, >=, 2);
  g_assert_cmpuint (failed, >=, 1);
  g_assert_cmpuint (flagged, >=, failed);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/function_sample_set/tracking", &test_ncm_function_sample_set_tracking);
  g_test_add_func ("/ncm/function_sample_set/expand_domain", &test_ncm_function_sample_set_expand_domain);
  g_test_add_func ("/ncm/function_sample_set/expand_domain/at_limits", &test_ncm_function_sample_set_expand_domain_at_limits);
  g_test_add_func ("/ncm/function_sample_set/expand_domain/traps", &test_ncm_function_sample_set_expand_domain_traps);
  g_test_add_func ("/ncm/function_sample_set/expand_domain/nonpositive/subprocess", &test_ncm_function_sample_set_expand_domain_nonpositive_subprocess);
  g_test_add_func ("/ncm/function_sample_set/expand_domain/empty/subprocess", &test_ncm_function_sample_set_expand_domain_empty_subprocess);
  g_test_add_func ("/ncm/function_sample_set/adaptive_midpoint/message", &test_ncm_function_sample_set_adaptive_midpoint_message);
  g_test_add_func ("/ncm/function_sample_set/adaptive_midpoint/precision", &test_ncm_function_sample_set_adaptive_midpoint_precision);
  g_test_add_func ("/ncm/function_sample_set/adaptive_midpoint/closure/none", &test_ncm_function_sample_set_adaptive_midpoint_closure_none);
  g_test_add_func ("/ncm/function_sample_set/adaptive_midpoint/closure/pass", &test_ncm_function_sample_set_adaptive_midpoint_closure_pass);
  g_test_add_func ("/ncm/function_sample_set/adaptive_midpoint/closure/fail", &test_ncm_function_sample_set_adaptive_midpoint_closure_fail);

  g_test_run ();
}


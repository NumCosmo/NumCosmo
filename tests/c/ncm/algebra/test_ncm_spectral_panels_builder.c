/***************************************************************************
 *            test_ncm_spectral_panels_builder.c
 *
 *  Sat Oct 10 2026
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/* Every point F is called at, and whether it lay in the interval it was called for */
typedef struct _Calls
{
  GArray *x;
  NcmVector *breakpoints;
  guint outside;
} Calls;

static Calls *
_calls_new (NcmVector *breakpoints)
{
  Calls *calls = g_new0 (Calls, 1);

  calls->x           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  calls->breakpoints = breakpoints;

  return calls;
}

static void
_calls_free (Calls *calls)
{
  g_array_unref (calls->x);
  g_free (calls);
}

static void
_calls_record (Calls *calls, guint i, gdouble x)
{
  if (calls == NULL)
    return;

  g_array_append_val (calls->x, x);

  if (calls->breakpoints != NULL)
  {
    const guint n_bp = ncm_vector_len (calls->breakpoints);

    if (((i > 0) && (x < ncm_vector_get (calls->breakpoints, i - 1))) ||
        ((i < n_bp) && (x > ncm_vector_get (calls->breakpoints, i))))
      calls->outside++;
  }
}

/* exp(x) and sin(5x) */
static void
_smooth (gpointer user_data, guint i, gdouble x, NcmVector *fval)
{
  _calls_record (user_data, i, x);
  ncm_vector_set (fval, 0, exp (x));
  ncm_vector_set (fval, 1, sin (5.0 * x));
}

/* A Runge function and a steep tanh, which need bisection on [-1, 1] */
static void
_hard (gpointer user_data, guint i, gdouble x, NcmVector *fval)
{
  _calls_record (user_data, i, x);
  ncm_vector_set (fval, 0, 1.0 / (1.0 + 25.0 * x * x));
  ncm_vector_set (fval, 1, tanh (20.0 * x));
}

/* With breakpoints at 0.7 and 1: x^2 + i, which jumps by one at each, and |x - 1| */
static gdouble
_piecewise_exact (guint r, guint i, gdouble x)
{
  return (r == 0) ? x * x + i : ((i < 2) ? 1.0 - x : x - 1.0);
}

static void
_piecewise (gpointer user_data, guint i, gdouble x, NcmVector *fval)
{
  _calls_record (user_data, i, x);
  ncm_vector_set (fval, 0, _piecewise_exact (0, i, x));
  ncm_vector_set (fval, 1, _piecewise_exact (1, i, x));
}

static void
_linear (gpointer user_data, guint i, gdouble x, NcmVector *fval)
{
  _calls_record (user_data, i, x);
  ncm_vector_set (fval, 0, 2.0 * x - 1.0);
}

/* A step at 0.3 that is not given as a breakpoint */
static void
_undeclared_step (gpointer user_data, guint i, gdouble x, NcmVector *fval)
{
  ncm_vector_set (fval, 0, (x < 0.3) ? 0.0 : 1.0);
}

static gint
_cmp_gdouble (gconstpointer pa, gconstpointer pb)
{
  const gdouble a = *(const gdouble *) pa;
  const gdouble b = *(const gdouble *) pb;

  return (a > b) - (a < b);
}

static NcmVector *
_vector_new (guint n, const gdouble *v)
{
  NcmVector *vec = ncm_vector_new (n);
  guint i;

  for (i = 0; i < n; i++)
    ncm_vector_set (vec, i, v[i]);

  return vec;
}

/* The largest |p_r(x) - F_r(x)| over n + 1 points of [a, b], relative to the largest |F_r| */
static gdouble
_sup_rel_error (NcmSpectralPanels *sp, NcmSpectralPanelsF F, guint r, gdouble a, gdouble b, guint n)
{
  NcmVector *fval = ncm_vector_new (ncm_spectral_panels_get_fdim (sp));
  gdouble err = 0.0, peak = 0.0;
  guint j;

  for (j = 0; j <= n; j++)
  {
    const gdouble x = a + (b - a) * j / n;

    F (NULL, 0, x, fval);
    err  = MAX (err, fabs (ncm_spectral_panels_eval_row (sp, r, x) - ncm_vector_get (fval, r)));
    peak = MAX (peak, fabs (ncm_vector_get (fval, r)));
  }

  ncm_vector_free (fval);

  return err / peak;
}

static void
_assert_contiguous (NcmSpectralPanels *sp, gdouble a, gdouble b)
{
  const guint n = ncm_spectral_panels_get_n_panels (sp);
  gdouble last  = a;
  guint i;

  g_assert_cmpuint (n, >, 0);

  for (i = 0; i < n; i++)
  {
    gdouble pa, pb;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);
    g_assert_cmpfloat (pa, ==, last);
    g_assert_cmpfloat (pb, >, pa);
    last = pb;
  }

  g_assert_cmpfloat (last, ==, b);
}

static void
_assert_same_panels (NcmSpectralPanels *sp1, NcmSpectralPanels *sp2)
{
  const guint n = ncm_spectral_panels_get_n_panels (sp1);
  guint i, k;

  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp2), ==, n);

  for (i = 0; i < n; i++)
  {
    gdouble a1, b1, a2, b2;
    guint N1, N2;
    const gdouble *c1 = ncm_spectral_panels_peek_panel (sp1, i, &a1, &b1, &N1);
    const gdouble *c2 = ncm_spectral_panels_peek_panel (sp2, i, &a2, &b2, &N2);

    g_assert_cmpfloat (a1, ==, a2);
    g_assert_cmpfloat (b1, ==, b2);
    g_assert_cmpuint (N1, ==, N2);

    for (k = 0; k < ncm_spectral_panels_get_fdim (sp1) * N1; k++)
      g_assert_cmpfloat (c1[k], ==, c2[k]);
  }
}

static NcmSpectralPanelsBuilder *
_builder_new (gdouble min_panels_per_efold)
{
  NcmSpectralPanelsBuilder *builder = ncm_spectral_panels_builder_new ();

  ncm_spectral_panels_builder_set_min_panels_per_efold (builder, min_panels_per_efold);

  return builder;
}

static void
test_ncm_spectral_panels_builder_defaults (void)
{
  NcmSpectralPanelsBuilder *builder = ncm_spectral_panels_builder_new ();

  g_assert_cmpuint (ncm_spectral_panels_builder_get_level_min (builder), ==, 3);
  g_assert_cmpuint (ncm_spectral_panels_builder_get_level_cap (builder), ==, 5);
  g_assert_cmpfloat (ncm_spectral_panels_builder_get_min_panels_per_efold (builder), ==, 1.0);
  g_assert_cmpfloat (ncm_spectral_panels_builder_get_min_panel_fraction (builder), ==, 1.0e-6);

  ncm_spectral_panels_builder_set_level_min (builder, 2);
  ncm_spectral_panels_builder_set_level_cap (builder, 7);
  ncm_spectral_panels_builder_set_min_panel_fraction (builder, 1.0e-8);
  g_assert_cmpuint (ncm_spectral_panels_builder_get_level_min (builder), ==, 2);
  g_assert_cmpuint (ncm_spectral_panels_builder_get_level_cap (builder), ==, 7);
  g_assert_cmpfloat (ncm_spectral_panels_builder_get_min_panel_fraction (builder), ==, 1.0e-8);

  ncm_spectral_panels_builder_free (builder);
}

/* A smooth F is resolved to the relative tolerance of each row */
static void
test_ncm_spectral_panels_builder_smooth (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (2);
  const gdouble reltol              = 1.0e-10;

  ncm_spectral_panels_builder_build (builder, _smooth, NULL, 0.5, 4.0, NULL, NULL, NULL, reltol, 0.0, sp);

  _assert_contiguous (sp, 0.5, 4.0);
  g_assert_cmpfloat (_sup_rel_error (sp, _smooth, 0, 0.5, 4.0, 2000), <=, reltol);
  g_assert_cmpfloat (_sup_rel_error (sp, _smooth, 1, 0.5, 4.0, 2000), <=, reltol);

  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (builder);
}

/* A panel that does not converge up to level-cap is bisected, and the halves reuse the
 * values at their ends: no point is evaluated twice */
static void
test_ncm_spectral_panels_builder_bisect (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (2);
  Calls *calls                      = _calls_new (NULL);
  const gdouble reltol              = 1.0e-12;
  guint j;

  ncm_spectral_panels_builder_build (builder, _hard, calls, -1.0, 1.0, NULL, NULL, NULL, reltol, 0.0, sp);

  _assert_contiguous (sp, -1.0, 1.0);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), >, 1);
  g_assert_cmpfloat (_sup_rel_error (sp, _hard, 0, -1.0, 1.0, 4000), <=, reltol);
  g_assert_cmpfloat (_sup_rel_error (sp, _hard, 1, -1.0, 1.0, 4000), <=, reltol);

  g_array_sort (calls->x, _cmp_gdouble);

  for (j = 1; j < calls->x->len; j++)
    g_assert_cmpfloat (g_array_index (calls->x, gdouble, j), !=, g_array_index (calls->x, gdouble, j - 1));

  _calls_free (calls);
  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (builder);
}

/* No panel straddles a breakpoint, F is called inside the interval of each panel, and at a
 * breakpoint the panel set takes the value from the right */
static void
test_ncm_spectral_panels_builder_breakpoints (void)
{
  const gdouble bp_v[]              = {0.7, 1.0};
  NcmVector *breakpoints            = _vector_new (2, bp_v);
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (2);
  Calls *calls                      = _calls_new (breakpoints);
  const gdouble edges[]             = {0.0, 0.7, 1.0, 2.0};
  guint i, j, r;

  ncm_spectral_panels_builder_build (builder, _piecewise, calls, 0.0, 2.0, breakpoints, NULL, NULL, 1.0e-12, 0.0, sp);

  g_assert_cmpuint (calls->outside, ==, 0);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 3);

  for (i = 0; i < 3; i++)
  {
    gdouble pa, pb;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);
    g_assert_cmpfloat (pa, ==, edges[i]);
    g_assert_cmpfloat (pb, ==, edges[i + 1]);

    for (j = 0; j < 100; j++)
    {
      const gdouble x = pa + (pb - pa) * j / 100.0;

      for (r = 0; r < 2; r++)
        g_assert_cmpfloat (fabs (ncm_spectral_panels_eval_row (sp, r, x) - _piecewise_exact (r, i, x)), <=, 1.0e-14);
    }
  }

  g_assert_cmpfloat (fabs (ncm_spectral_panels_eval_row (sp, 0, 0.7) - (0.49 + 1.0)), <=, 1.0e-14);
  g_assert_cmpfloat (fabs (ncm_spectral_panels_eval_row (sp, 0, 1.0) - 3.0), <=, 1.0e-14);

  _calls_free (calls);
  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (builder);
  ncm_vector_free (breakpoints);
}

/* Known values on each side of the breakpoints replace two evaluations per breakpoint and
 * change nothing else */
static void
test_ncm_spectral_panels_builder_known_values (void)
{
  const gdouble bp_v[]              = {0.7, 1.0};
  NcmVector *breakpoints            = _vector_new (2, bp_v);
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (2);
  NcmSpectralPanels *sp_known       = ncm_spectral_panels_new (2);
  NcmMatrix *f_minus                = ncm_matrix_new (2, 2);
  NcmMatrix *f_plus                 = ncm_matrix_new (2, 2);
  Calls *calls                      = _calls_new (breakpoints);
  Calls *calls_known                = _calls_new (breakpoints);
  guint j, r;

  for (j = 0; j < 2; j++)
  {
    for (r = 0; r < 2; r++)
    {
      ncm_matrix_set (f_minus, j, r, _piecewise_exact (r, j, bp_v[j]));
      ncm_matrix_set (f_plus, j, r, _piecewise_exact (r, j + 1, bp_v[j]));
    }
  }

  ncm_spectral_panels_builder_build (builder, _piecewise, calls, 0.0, 2.0, breakpoints, NULL, NULL, 1.0e-12, 0.0, sp);
  ncm_spectral_panels_builder_build (builder, _piecewise, calls_known, 0.0, 2.0, breakpoints, f_minus, f_plus, 1.0e-12, 0.0, sp_known);

  g_assert_cmpuint (calls_known->x->len + 4, ==, calls->x->len);
  _assert_same_panels (sp, sp_known);

  _calls_free (calls_known);
  _calls_free (calls);
  ncm_matrix_free (f_plus);
  ncm_matrix_free (f_minus);
  ncm_spectral_panels_free (sp_known);
  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (builder);
  ncm_vector_free (breakpoints);
}

/* The initial grid is exp (j / p) inside the range; a grid point within a quarter of the
 * spacing of a breakpoint is dropped */
static void
test_ncm_spectral_panels_builder_grid (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (2.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);
  const gdouble bp_v[]              = {exp (1.05)};
  NcmVector *breakpoints            = _vector_new (1, bp_v);
  const gdouble grid[]              = {1.0, exp (0.5), exp (1.0), exp (1.5), exp (2.0), exp (2.5), exp (3.0)};
  const gdouble grid_bp[]           = {1.0, exp (0.5), exp (1.05), exp (1.5), exp (2.0), exp (2.5), exp (3.0)};
  guint i;

  ncm_spectral_panels_builder_build (builder, _linear, NULL, 1.0, exp (3.0), NULL, NULL, NULL, 1.0e-12, 0.0, sp);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 6);

  for (i = 0; i < 6; i++)
  {
    gdouble pa, pb;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);
    g_assert_cmpfloat (pa, ==, grid[i]);
    g_assert_cmpfloat (pb, ==, grid[i + 1]);
  }

  ncm_spectral_panels_builder_build (builder, _linear, NULL, 1.0, exp (3.0), breakpoints, NULL, NULL, 1.0e-12, 0.0, sp);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 6);

  for (i = 0; i < 6; i++)
  {
    gdouble pa, pb;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);
    g_assert_cmpfloat (pa, ==, grid_bp[i]);
    g_assert_cmpfloat (pb, ==, grid_bp[i + 1]);
  }

  ncm_vector_free (breakpoints);
  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (builder);
}

static gboolean
_accept_narrow (gpointer user_data, gdouble a, gdouble b, const gdouble *c, guint N, guint tda)
{
  return (b - a) <= 0.25;
}

/* A converged panel that accept refuses is bisected */
static void
test_ncm_spectral_panels_builder_accept (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);
  guint i;

  ncm_spectral_panels_builder_build_full (builder, _linear, _accept_narrow, NULL, 0.0, 1.0, NULL, NULL, NULL, 1.0e-12, 0.0, sp);

  _assert_contiguous (sp, 0.0, 1.0);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 4);

  for (i = 0; i < 4; i++)
  {
    gdouble pa, pb;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &N);
    g_assert_cmpfloat (pb - pa, ==, 0.25);
  }

  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (builder);
}

/* A builder and a panel set reused across builds, of another dimension too, give the
 * panels of a fresh pair */
static void
test_ncm_spectral_panels_builder_reuse (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanelsBuilder *fresh   = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (2);
  NcmSpectralPanels *sp_fresh       = ncm_spectral_panels_new (2);
  NcmSpectralPanels *sp1            = ncm_spectral_panels_new (1);
  NcmSpectralPanels *sp1_fresh      = ncm_spectral_panels_new (1);

  ncm_spectral_panels_builder_build (builder, _hard, NULL, -1.0, 1.0, NULL, NULL, NULL, 1.0e-12, 0.0, sp);
  ncm_spectral_panels_builder_build (builder, _linear, NULL, 0.0, 1.0, NULL, NULL, NULL, 1.0e-12, 0.0, sp1);
  ncm_spectral_panels_builder_build (builder, _smooth, NULL, 0.5, 4.0, NULL, NULL, NULL, 1.0e-10, 0.0, sp);

  ncm_spectral_panels_builder_build (fresh, _smooth, NULL, 0.5, 4.0, NULL, NULL, NULL, 1.0e-10, 0.0, sp_fresh);
  ncm_spectral_panels_builder_build (fresh, _linear, NULL, 0.0, 1.0, NULL, NULL, NULL, 1.0e-12, 0.0, sp1_fresh);

  _assert_same_panels (sp, sp_fresh);
  _assert_same_panels (sp1, sp1_fresh);

  ncm_spectral_panels_free (sp1_fresh);
  ncm_spectral_panels_free (sp1);
  ncm_spectral_panels_free (sp_fresh);
  ncm_spectral_panels_free (sp);
  ncm_spectral_panels_builder_free (fresh);
  ncm_spectral_panels_builder_free (builder);
}

typedef struct _ThreadBuild
{
  NcmSpectralPanelsBuilder *builder;
  NcmSpectralPanels *sp;
} ThreadBuild;

static gpointer
_thread_build (gpointer data)
{
  ThreadBuild *tb = data;
  guint k;

  for (k = 0; k < 20; k++)
    ncm_spectral_panels_builder_build (tb->builder, _hard, NULL, -1.0, 1.0, NULL, NULL, NULL, 1.0e-12, 0.0, tb->sp);

  return NULL;
}

/* Builds in two threads on one builder give the serial panels */
static void
test_ncm_spectral_panels_builder_threads (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *serial         = ncm_spectral_panels_new (2);
  ThreadBuild tb[2]                 = {
    {
      builder, ncm_spectral_panels_new (2)
    }, {
      builder, ncm_spectral_panels_new (2)
    }
  };
  GThread *threads[2];
  guint t;

  ncm_spectral_panels_builder_build (builder, _hard, NULL, -1.0, 1.0, NULL, NULL, NULL, 1.0e-12, 0.0, serial);

  for (t = 0; t < 2; t++)
    threads[t] = g_thread_new ("build", _thread_build, &tb[t]);

  for (t = 0; t < 2; t++)
  {
    g_thread_join (threads[t]);
    _assert_same_panels (serial, tb[t].sp);
    ncm_spectral_panels_free (tb[t].sp);
  }

  ncm_spectral_panels_free (serial);
  ncm_spectral_panels_builder_free (builder);
}

static void
test_ncm_spectral_panels_builder_undeclared_step (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_builder/undeclared_step/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*without converging*");
}

static void
test_ncm_spectral_panels_builder_undeclared_step_subprocess (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);

  ncm_spectral_panels_builder_build (builder, _undeclared_step, NULL, 0.1, 1.0, NULL, NULL, NULL, 1.0e-10, 0.0, sp);
}

static void
test_ncm_spectral_panels_builder_levels (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_builder/levels/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not below level-cap*");
}

static void
test_ncm_spectral_panels_builder_levels_subprocess (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);

  ncm_spectral_panels_builder_set_level_min (builder, 5);
  ncm_spectral_panels_builder_build (builder, _linear, NULL, 0.0, 1.0, NULL, NULL, NULL, 1.0e-10, 0.0, sp);
}

static void
test_ncm_spectral_panels_builder_grid_range (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_builder/grid_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*needs a > 0*");
}

static void
test_ncm_spectral_panels_builder_grid_range_subprocess (void)
{
  NcmSpectralPanelsBuilder *builder = _builder_new (1.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);

  ncm_spectral_panels_builder_build (builder, _linear, NULL, 0.0, 1.0, NULL, NULL, NULL, 1.0e-10, 0.0, sp);
}

static void
test_ncm_spectral_panels_builder_breakpoint_order (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_builder/breakpoint_order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not ascending inside*");
}

static void
test_ncm_spectral_panels_builder_breakpoint_order_subprocess (void)
{
  const gdouble bp_v[]              = {0.6, 0.4};
  NcmVector *breakpoints            = _vector_new (2, bp_v);
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);

  ncm_spectral_panels_builder_build (builder, _linear, NULL, 0.0, 1.0, breakpoints, NULL, NULL, 1.0e-10, 0.0, sp);
}

static void
test_ncm_spectral_panels_builder_known_values_pair (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels_builder/known_values_pair/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*given together*");
}

static void
test_ncm_spectral_panels_builder_known_values_pair_subprocess (void)
{
  const gdouble bp_v[]              = {0.5};
  NcmVector *breakpoints            = _vector_new (1, bp_v);
  NcmMatrix *f_minus                = ncm_matrix_new (1, 1);
  NcmSpectralPanelsBuilder *builder = _builder_new (0.0);
  NcmSpectralPanels *sp             = ncm_spectral_panels_new (1);

  ncm_spectral_panels_builder_build (builder, _linear, NULL, 0.0, 1.0, breakpoints, f_minus, NULL, 1.0e-10, 0.0, sp);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/spectral_panels_builder/defaults", &test_ncm_spectral_panels_builder_defaults);
  g_test_add_func ("/ncm/spectral_panels_builder/smooth", &test_ncm_spectral_panels_builder_smooth);
  g_test_add_func ("/ncm/spectral_panels_builder/bisect", &test_ncm_spectral_panels_builder_bisect);
  g_test_add_func ("/ncm/spectral_panels_builder/breakpoints", &test_ncm_spectral_panels_builder_breakpoints);
  g_test_add_func ("/ncm/spectral_panels_builder/known_values", &test_ncm_spectral_panels_builder_known_values);
  g_test_add_func ("/ncm/spectral_panels_builder/grid", &test_ncm_spectral_panels_builder_grid);
  g_test_add_func ("/ncm/spectral_panels_builder/accept", &test_ncm_spectral_panels_builder_accept);
  g_test_add_func ("/ncm/spectral_panels_builder/reuse", &test_ncm_spectral_panels_builder_reuse);
  g_test_add_func ("/ncm/spectral_panels_builder/threads", &test_ncm_spectral_panels_builder_threads);
  g_test_add_func ("/ncm/spectral_panels_builder/undeclared_step", &test_ncm_spectral_panels_builder_undeclared_step);
  g_test_add_func ("/ncm/spectral_panels_builder/undeclared_step/subprocess", &test_ncm_spectral_panels_builder_undeclared_step_subprocess);
  g_test_add_func ("/ncm/spectral_panels_builder/levels", &test_ncm_spectral_panels_builder_levels);
  g_test_add_func ("/ncm/spectral_panels_builder/levels/subprocess", &test_ncm_spectral_panels_builder_levels_subprocess);
  g_test_add_func ("/ncm/spectral_panels_builder/grid_range", &test_ncm_spectral_panels_builder_grid_range);
  g_test_add_func ("/ncm/spectral_panels_builder/grid_range/subprocess", &test_ncm_spectral_panels_builder_grid_range_subprocess);
  g_test_add_func ("/ncm/spectral_panels_builder/breakpoint_order", &test_ncm_spectral_panels_builder_breakpoint_order);
  g_test_add_func ("/ncm/spectral_panels_builder/breakpoint_order/subprocess", &test_ncm_spectral_panels_builder_breakpoint_order_subprocess);
  g_test_add_func ("/ncm/spectral_panels_builder/known_values_pair", &test_ncm_spectral_panels_builder_known_values_pair);
  g_test_add_func ("/ncm/spectral_panels_builder/known_values_pair/subprocess", &test_ncm_spectral_panels_builder_known_values_pair_subprocess);

  g_test_run ();

  return 0;
}


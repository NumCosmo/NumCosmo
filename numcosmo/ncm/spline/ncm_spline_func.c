/***************************************************************************
 *            ncm_spline_func.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

/**
 * NcmSplineFunc:
 *
 * Knot placement for a #NcmSpline interpolating a function.
 *
 * ncm_spline_set_func() and ncm_spline_set_func_scale() place the knots adaptively, and
 * ncm_spline_set_func_grid() on a fixed grid, according to #NcmSplineFuncType.
 *
 * The adaptive types, %NCM_SPLINE_FUNCTION_SPLINE, %NCM_SPLINE_FUNCTION_SPLINE_LNKNOT and
 * %NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, implement AutoKnots, [Vitenti et al.
 * (2025)](https://doi.org/10.1016/j.ascom.2025.100970), with the notation used there. They
 * start from $\max(m, 3)$ knots uniform in $x$, $\ln x$ or $\sinh^{-1} x$, where $m$ is
 * ncm_spline_min_size(). Each round visits every interval $[x_i, x_{i+1}]$ not yet
 * accepted, evaluates $f$ at its midpoint $\overline{x}_i$ in the same variable and
 * inserts $\overline{x}_i$ as a knot. Both new intervals are accepted when
 * $$|f(\overline{x}_i) - \hat{f}(\overline{x}_i)| \le \delta\,(|f(\overline{x}_i)| +
 * \varepsilon) \quad\text{and}\quad |\widetilde{\mathcal{I}}_i - \hat{\mathcal{I}}_i| \le
 * \delta\,(|\widetilde{\mathcal{I}}_i| + \varepsilon\,h_i),$$
 * where $\hat{f}$ is the spline on the knots of the previous round, $\hat{\mathcal{I}}_i$
 * its integral over the interval, $\widetilde{\mathcal{I}}_i$ the integral of the
 * quadratic through the three points (Simpson's rule for %NCM_SPLINE_FUNCTION_SPLINE),
 * $h_i = x_{i+1} - x_i$, $\delta$ the relative tolerance and $\varepsilon$ the scale (zero
 * for ncm_spline_set_func()). The rounds end when a round accepts every interval it
 * visits.
 *
 * A midpoint can pass by coincidence where the spline is still wrong, and the two intervals
 * it accepts are then left wider than the mesh around them once that mesh has settled. So
 * the settled mesh is closed: every interval wider in the knot variable than the mean of
 * its two neighbors (a border interval, than its one neighbor), by more than one part in
 * $10^8$ (bisection leaves exact factors of two), has its midpoint evaluated and tested as
 * above, the point kept as a knot; the
 * intervals that fail are split, the rounds resume, and the closure repeats until a closure
 * round flags nothing or every flagged interval passes. The closure removes errors from an
 * adaptive mesh and discovers no feature: a feature the starting knots and their midpoints
 * carry no evidence of is found by no test on those samples.
 *
 * The adaptive types stop with a warning when the number of knots exceeds `max_nodes`,
 * unlimited when zero, and abort when an interval becomes shorter than
 * %NCM_SPLINE_KNOT_DIFF_TOL relative to its midpoint, which indicates a discontinuity.
 *
 * <inlinegraphic fileref="spline_func_knots_evolution.png" format="PNG" scale="98" align="right"/>
 *
 * The figure shows three rounds of %NCM_SPLINE_FUNCTION_SPLINE, from $6$ knots (first
 * line) to $19$ (last line); the markers on each line are the midpoints tested in that
 * round, placed only in the intervals that failed the previous one.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/spline/ncm_spline_func.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"

typedef struct
{
  gdouble x;
  gdouble y;
  gint ok;
} _BIVec;

#define BIVEC_LIST_APPEND(dlist, xX, yY)              \
        do {                                          \
          _BIVec *bv = g_slice_new (_BIVec);          \
          bv->x = (xX); bv->y = (yY); bv->ok = FALSE; \
          dlist = g_list_append (dlist, bv);          \
        } while (FALSE)
#define BIVEC_LIST_INSERT_BEFORE(dlist, node, xX, yY)     \
        do {                                              \
          _BIVec *bv = g_slice_new (_BIVec);              \
          bv->x = (xX); bv->y = (yY); bv->ok = FALSE;     \
          dlist = g_list_insert_before (dlist, node, bv); \
        } while (FALSE)

#define BIVEC_LIST_X(dlist) (((_BIVec *) (dlist)->data)->x)
#define BIVEC_LIST_Y(dlist) (((_BIVec *) (dlist)->data)->y)
#define BIVEC_LIST_OK(dlist) (((_BIVec *) (dlist)->data)->ok)

/* Relative margin of the closure's width comparison: bisection leaves exact factors of two. */
#define NCM_SPLINE_FUNC_CLOSURE_MARGIN (1.0e-8)

static void
_BIVec_free (gpointer mem)
{
  g_slice_free (_BIVec, mem);
}

/* The knot variable of each adaptive type, u (x), and its inverse. */
static gdouble
_u_identity (gdouble x)
{
  return x;
}

static gdouble
_u_log (gdouble x)
{
  return log (x);
}

static gdouble
_u_exp (gdouble u)
{
  return exp (u);
}

static gdouble
_u_asinh (gdouble x)
{
  return asinh (x);
}

static gdouble
_u_sinh (gdouble u)
{
  return sinh (u);
}

/* The closure of the settled mesh: rounds run, intervals flagged, intervals that failed. */
typedef struct
{
  guint rounds;
  guint flagged;
  guint failed;
} _ClosureStats;

static guint _ncm_spline_func_flag_wide (GList *nodes, gdouble (*fwd) (gdouble));

/*
 * The adaptive placement on [xi, xf] with knots bisected in u = fwd (x), x = inv (u):
 * the rounds of midpoint tests, then the closure of the settled mesh, see NcmSplineFunc.
 */
static void
_ncm_spline_new_function_adaptive (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble f_scale, gdouble (*fwd) (gdouble), gdouble (*inv) (gdouble), _ClosureStats *stats)
{
  GArray *x_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *xt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *yt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GList *nodes = NULL, *wnodes = NULL;
  gsize n          = ncm_spline_min_size (s);
  const gdouble ui = fwd (xi);
  const gdouble uf = fwd (xf);
  gboolean closing = FALSE; /* the round under way tests the flagged intervals */
  guint n_flagged  = 0;
  guint i;

  n = (n < 3) ? 3 : n;

  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);
  g_assert (gsl_finite (ui) && gsl_finite (uf));
  g_assert_cmpfloat (f_scale, >=, 0.0);

  max_nodes = (max_nodes <= 0) ? G_MAXUINT64 : max_nodes;

  g_array_set_size (xt_array, n);
  g_array_set_size (yt_array, n);

  for (i = 0; i < n; i++)
  {
    const gdouble x = inv (ui + (uf - ui) / (n - 1.0) * i);
    const gdouble y = GSL_FN_EVAL (F, x);

    BIVEC_LIST_APPEND (nodes, x, y);
    BIVEC_LIST_OK (nodes) = 0;

    g_array_append_val (x_array, x);
    g_array_append_val (y_array, y);
    g_assert (gsl_finite (x));
    g_assert (gsl_finite (y));
  }

  ncm_spline_set_array (s, x_array, y_array, TRUE);

#define SWAP_PTR(a, b)                                    \
        do {                                              \
          const gpointer tmp = (b); (b) = (a); (a) = tmp; \
        } while (FALSE)

  while (TRUE)
  {
    gsize improves = 0;

    wnodes = nodes;
    g_array_set_size (xt_array, 0);
    g_array_set_size (yt_array, 0);

    do {
      g_array_append_val (xt_array, BIVEC_LIST_X (wnodes));
      g_array_append_val (yt_array, BIVEC_LIST_Y (wnodes));

      if (BIVEC_LIST_OK (wnodes) == 1)
      {
        continue;
      }
      else
      {
        const gdouble x0      = BIVEC_LIST_X (wnodes);
        const gdouble x1      = BIVEC_LIST_X (wnodes->next);
        const gdouble y0      = BIVEC_LIST_Y (wnodes);
        const gdouble y1      = BIVEC_LIST_Y (wnodes->next);
        const gdouble x       = inv (0.5 * (fwd (x0) + fwd (x1)));
        const gdouble y       = GSL_FN_EVAL (F, x);
        const gdouble ys      = ncm_spline_eval (s, x);
        const gdouble delta   = (x - 0.5 * (x1 + x0)) / (0.5 * (x1 - x0));
        const gdouble Iyc     = (x1 - x0) * (y1 * (3.0 - 2.0 / (1.0 - delta)) + y0 * (3.0 - 2.0 / (1.0 + delta)) + 4.0 * y / (1.0 - delta * delta)) / 6.0;
        const gdouble Iys     = ncm_spline_eval_integ (s, x0, x1);
        const gboolean test_p = fabs (y - ys)    <= rel_error * (fabs (y)   + f_scale);
        const gboolean test_I = fabs (Iyc - Iys) <= rel_error * (fabs (Iyc) + f_scale * (x1 - x0));

        if (fabs ((x - x0) / x) < NCM_SPLINE_KNOT_DIFF_TOL)
          g_error ("Tolerance of the difference between knots was reached. Interpolated function is probably discontinuous at x = (% 20.15g, % 20.15g, % 20.15g).\n"
                   "\tFunction value at f(x0) = % 22.15g, f(x) = % 22.15g and f(x1) = % 22.15g, cmp (%e, %e).",
                   x0, x, x1,
                   y0, y, y1,
                   fabs (y0 / y - 1.0),
                   fabs (y1 / y - 1.0));

        BIVEC_LIST_INSERT_BEFORE (nodes, wnodes->next, x, y);
        wnodes = g_list_next (wnodes);
        g_array_append_val (xt_array, BIVEC_LIST_X (wnodes));
        g_array_append_val (yt_array, BIVEC_LIST_Y (wnodes));
        BIVEC_LIST_OK (wnodes) = BIVEC_LIST_OK (wnodes->prev);

        if (test_p && test_I)
        {
          BIVEC_LIST_OK (wnodes->prev)++;
          BIVEC_LIST_OK (wnodes)++;
        }
        else
        {
          improves++;
        }
      }
    } while ((wnodes = g_list_next (wnodes)) && wnodes->next);

    if (wnodes != NULL)
    {
      g_array_append_val (xt_array, BIVEC_LIST_X (wnodes));
      g_array_append_val (yt_array, BIVEC_LIST_Y (wnodes));
    }

    SWAP_PTR (x_array, xt_array);
    SWAP_PTR (y_array, yt_array);

    ncm_spline_set_array (s, x_array, y_array, TRUE);

    if (x_array->len > max_nodes)
    {
      g_warning ("ncm_spline_set_func: cannot achieve requested precision with at most %zu nodes", max_nodes);
      break;
    }

    if (closing)
    {
      stats->rounds++;
      stats->flagged += n_flagged;
      stats->failed  += improves;

      if (improves == 0)
        break;

      closing = FALSE;
      continue;
    }

    if (improves == 0)
    {
      n_flagged = _ncm_spline_func_flag_wide (nodes, fwd);

      if (n_flagged == 0)
      {
        stats->rounds++;
        break;
      }

      closing = TRUE;
    }
  }

  g_list_free_full (nodes, _BIVec_free);

  g_array_unref (x_array);
  g_array_unref (xt_array);
  g_array_unref (y_array);
  g_array_unref (yt_array);
}

/*
 * Reopens, on a settled mesh, every interval wider in u = fwd (x) than the mean of its two
 * neighbors by more than NCM_SPLINE_FUNC_CLOSURE_MARGIN, a border interval than its one
 * neighbor, and returns how many it reopened.
 */
static guint
_ncm_spline_func_flag_wide (GList *nodes, gdouble (*fwd) (gdouble))
{
  guint n_flagged = 0;
  GList *w;

  for (w = nodes; (w != NULL) && (w->next != NULL); w = w->next)
  {
    const gdouble h = fwd (BIVEC_LIST_X (w->next)) - fwd (BIVEC_LIST_X (w));
    gdouble h_mean  = 0.0;

    if ((w->prev != NULL) && (w->next->next != NULL))
      h_mean = 0.5 * (fwd (BIVEC_LIST_X (w)) - fwd (BIVEC_LIST_X (w->prev)) + fwd (BIVEC_LIST_X (w->next->next)) - fwd (BIVEC_LIST_X (w->next)));
    else if (w->prev != NULL)
      h_mean = fwd (BIVEC_LIST_X (w)) - fwd (BIVEC_LIST_X (w->prev));
    else if (w->next->next != NULL)
      h_mean = fwd (BIVEC_LIST_X (w->next->next)) - fwd (BIVEC_LIST_X (w->next));
    else
      continue;

    if (h > h_mean * (1.0 + NCM_SPLINE_FUNC_CLOSURE_MARGIN))
    {
      BIVEC_LIST_OK (w) = 0;
      n_flagged++;
    }
  }

  return n_flagged;
}

/**
 * ncm_spline_set_func: (skip)
 * @s: a #NcmSpline
 * @ftype: a #NcmSplineFuncType, one of the adaptive types
 * @F: the function
 * @xi: the lower limit
 * @xf: the upper limit
 * @max_nodes: the maximum number of knots
 * @rel_error: the relative tolerance
 *
 * Places the knots of @s on [@xi, @xf] adaptively and prepares it, see #NcmSplineFunc,
 * with scale zero.
 */
void
ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error)
{
  ncm_spline_set_func_scale (s, ftype, F, xi, xf, max_nodes, rel_error, 0.0);
}

/**
 * ncm_spline_set_func_scale: (skip)
 * @s: a #NcmSpline
 * @ftype: a #NcmSplineFuncType, one of the adaptive types
 * @F: the function
 * @xi: the lower limit
 * @xf: the upper limit
 * @max_nodes: the maximum number of knots
 * @rel_error: the relative tolerance
 * @scale: the scale of the function values
 *
 * Places the knots of @s on [@xi, @xf] adaptively and prepares it, see #NcmSplineFunc;
 * the absolute tolerance is @rel_error times @scale.
 */
void
ncm_spline_set_func_scale (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble scale)
{
  ncm_spline_set_func_full (s, ftype, F, xi, xf, max_nodes, rel_error, scale, NULL, NULL, NULL);
}

/**
 * ncm_spline_set_func_full: (skip)
 * @s: a #NcmSpline
 * @ftype: a #NcmSplineFuncType, one of the adaptive types
 * @F: the function
 * @xi: the lower limit
 * @xf: the upper limit
 * @max_nodes: the maximum number of knots
 * @rel_error: the relative tolerance
 * @scale: the scale of the function values
 * @closure_rounds: (out) (allow-none): number of closure rounds
 * @closure_flagged: (out) (allow-none): number of intervals the closure rounds flagged
 * @closure_failed: (out) (allow-none): number of flagged intervals that failed
 *
 * Same as ncm_spline_set_func_scale(), reporting the closure of the settled mesh, see
 * #NcmSplineFunc. A closure that flagged nothing counts one round with nothing flagged;
 * one whose flagged intervals all passed counts their number with nothing failed.
 */
void
ncm_spline_set_func_full (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble scale, guint *closure_rounds, guint *closure_flagged, guint *closure_failed)
{
  _ClosureStats stats = {0, 0, 0};

  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);

  switch (ftype)
  {
    case NCM_SPLINE_FUNCTION_SPLINE:
      _ncm_spline_new_function_adaptive (s, F, xi, xf, max_nodes, rel_error, scale, &_u_identity, &_u_identity, &stats);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE_LNKNOT:
      g_assert_cmpfloat (xi, >, 0.0);
      _ncm_spline_new_function_adaptive (s, F, xi, xf, max_nodes, rel_error, scale, &_u_log, &_u_exp, &stats);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT:
      _ncm_spline_new_function_adaptive (s, F, xi, xf, max_nodes, rel_error, scale, &_u_asinh, &_u_sinh, &stats);
      break;
    default:
      g_assert_not_reached ();

      return;
  }

  if (closure_rounds != NULL)
    *closure_rounds = stats.rounds;

  if (closure_flagged != NULL)
    *closure_flagged = stats.flagged;

  if (closure_failed != NULL)
    *closure_failed = stats.failed;
}

/**
 * ncm_spline_set_func1:
 * @s: a #NcmSpline
 * @ftype: a #NcmSplineFuncType, one of the adaptive types
 * @F: (scope call): the function
 * @obj: (allow-none): the object passed to @F
 * @xi: the lower limit
 * @xf: the upper limit
 * @max_nodes: the maximum number of knots
 * @rel_error: the relative tolerance
 *
 * Same as ncm_spline_set_func(), for a #NcmSplineFuncF.
 */
void
ncm_spline_set_func1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
{
  gsl_function gslF = {(gdouble (*)(gdouble, gpointer)) F, obj};

  ncm_spline_set_func (s, ftype, &gslF, xi, xf, max_nodes, rel_error);
}

static void
_ncm_spline_new_function_grid_linear (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize nnodes)
{
  guint i;
  NcmVector *xv = ncm_vector_new (nnodes);
  NcmVector *yv = ncm_vector_new (nnodes);

  for (i = 0; i < nnodes; i++)
  {
    const gdouble x = xi + (xf - xi) / (nnodes - 1.0) * i;
    const gdouble y = GSL_FN_EVAL (F, x);

    g_assert (gsl_finite (x));
    g_assert (gsl_finite (y));

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, y);
  }

  ncm_spline_set (s, xv, yv, TRUE);

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return;
}

static void
_ncm_spline_new_function_grid_log (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize nnodes)
{
  guint i;
  const gdouble lnxi = log (xi);
  const gdouble lnxf = log (xf);
  NcmVector *xv      = ncm_vector_new (nnodes);
  NcmVector *yv      = ncm_vector_new (nnodes);

  g_assert (xi > 0.0);

  for (i = 0; i < nnodes; i++)
  {
    const gdouble x = exp (lnxi + (lnxf - lnxi) / (nnodes - 1.0) * i);
    const gdouble y = GSL_FN_EVAL (F, x);

    g_assert (gsl_finite (x));
    g_assert (gsl_finite (y));

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, y);
  }

  ncm_spline_set (s, xv, yv, TRUE);

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return;
}

/**
 * ncm_spline_set_func_grid: (skip)
 * @s: a #NcmSpline
 * @ftype: a #NcmSplineFuncType, %NCM_SPLINE_FUNC_GRID_LINEAR or %NCM_SPLINE_FUNC_GRID_LOG
 * @F: the function
 * @xi: the lower limit
 * @xf: the upper limit
 * @nnodes: the number of knots
 *
 * Sets @s to @nnodes knots uniform on [@xi, @xf] in $x$ or in $\ln x$, both limits
 * included, and prepares it. @nnodes must be larger than ncm_spline_min_size() and
 * smaller than %NCM_SPLINE_FUNC_DEFAULT_MAX_NODES; %NCM_SPLINE_FUNC_GRID_LOG requires
 * @xi > 0.
 */
void
ncm_spline_set_func_grid (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize nnodes)
{
  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);

  g_assert_cmpuint (nnodes, <, NCM_SPLINE_FUNC_DEFAULT_MAX_NODES);
  g_assert_cmpuint (nnodes, >, ncm_spline_min_size (s));

  switch (ftype)
  {
    case NCM_SPLINE_FUNC_GRID_LINEAR:
      _ncm_spline_new_function_grid_linear (s, F, xi, xf, nnodes);
      break;
    case NCM_SPLINE_FUNC_GRID_LOG:
      _ncm_spline_new_function_grid_log (s, F, xi, xf, nnodes);
      break;
    default:
      g_assert_not_reached ();

      return;
  }
}

/**
 * ncm_spline_set_func_grid1:
 * @s: a #NcmSpline
 * @ftype: a #NcmSplineFuncType, %NCM_SPLINE_FUNC_GRID_LINEAR or %NCM_SPLINE_FUNC_GRID_LOG
 * @F: (scope call): the function
 * @obj: (allow-none): the object passed to @F
 * @xi: the lower limit
 * @xf: the upper limit
 * @nnodes: the number of knots
 *
 * Same as ncm_spline_set_func_grid(), for a #NcmSplineFuncF.
 */
void
ncm_spline_set_func_grid1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize nnodes)
{
  gsl_function gslF = {(gdouble (*)(gdouble, gpointer)) F, obj};

  ncm_spline_set_func_grid (s, ftype, &gslF, xi, xf, nnodes);
}


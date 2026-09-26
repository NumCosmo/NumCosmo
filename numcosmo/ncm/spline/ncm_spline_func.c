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
 * visits. For %NCM_SPLINE_FUNCTION_SPLINE, the intervals with $h_i$ larger than the mean
 * plus `refine_ns` standard deviations are then reopened and the rounds resumed, `refine`
 * times, see ncm_spline_set_func_scale(); ncm_spline_set_func() uses one pass with one
 * standard deviation.
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
#include "ncm/stats/ncm_stats_vec.h"

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

static void
_BIVec_free (gpointer mem)
{
  g_slice_free (_BIVec, mem);
}

static void
ncm_spline_new_function_spline (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble f_scale, gint refine, gdouble refine_ns)
{
  GArray *x_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *xt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *yt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GList *nodes = NULL, *wnodes = NULL;
  NcmStatsVec *dx_stats = ncm_stats_vec_new (1, NCM_STATS_VEC_VAR, FALSE);
  gsize n               = ncm_spline_min_size (s);
  gdouble max_dx, min_dx;
  guint i;

  n = (n < 3) ? 3 : n;

  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);
  g_assert_cmpfloat (f_scale, >=, 0.0);

  max_nodes = (max_nodes <= 0) ? G_MAXUINT64 : max_nodes;

  g_array_set_size (xt_array, n);
  g_array_set_size (yt_array, n);

  for (i = 0; i < n; i++)
  {
    const gdouble x = xi + (xf - xi) / (n - 1.0) * i;
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

    ncm_stats_vec_reset (dx_stats, TRUE);

    max_dx = 0.0;
    min_dx = 1.0e300;

    do {
      const gdouble x0 = BIVEC_LIST_X (wnodes);
      const gdouble x1 = BIVEC_LIST_X (wnodes->next);
      const gdouble dx = x1 - x0;

      max_dx = MAX (max_dx, dx);
      min_dx = MIN (min_dx, dx);

      g_array_append_val (xt_array, BIVEC_LIST_X (wnodes));
      g_array_append_val (yt_array, BIVEC_LIST_Y (wnodes));

      ncm_stats_vec_set (dx_stats, 0, x1 - x0);
      ncm_stats_vec_update (dx_stats);

      if (BIVEC_LIST_OK (wnodes) == 1)
      {
        continue;
      }
      else
      {
        const gdouble y0      = BIVEC_LIST_Y (wnodes);
        const gdouble y1      = BIVEC_LIST_Y (wnodes->next);
        const gdouble x       = (x0 + x1) / 2.0;
        const gdouble y       = GSL_FN_EVAL (F, x);
        const gdouble ys      = ncm_spline_eval (s, x);
        const gdouble Iyc     = dx * (y1 + y0 + 4.0 * y) / 6.0;
        const gdouble Iys     = ncm_spline_eval_integ (s, x0, x1);
        const gboolean test_p = fabs (y - ys)    <= rel_error * (fabs (y)   + f_scale);
        const gboolean test_I = fabs (Iyc - Iys) <= rel_error * (fabs (Iyc) + f_scale * dx);

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
      g_warning ("ncm_spline_new_function_spline: cannot achieve requested precision with at most %zu nodes", max_nodes);
      break;
    }

    if (improves == 0)
    {
      if (refine < 1)
      {
        break;
      }
      else
      {
        const gdouble dx_mean = ncm_stats_vec_get_mean (dx_stats, 0);
        const gdouble dx_sd   = ncm_stats_vec_get_sd (dx_stats, 0);
        const gdouble dx_lim  = refine_ns * dx_sd + dx_mean;

        refine--;

        if (max_dx > dx_lim)
        {
          wnodes = nodes;

          do {
            const gdouble x0 = BIVEC_LIST_X (wnodes);
            const gdouble x1 = BIVEC_LIST_X (wnodes->next);
            const gdouble dx = x1 - x0;

            if (dx > dx_lim)
              BIVEC_LIST_OK (wnodes) = 0;
          } while ((wnodes = g_list_next (wnodes)) && wnodes->next);
        }
        else
        {
          break;
        }
      }
    }
  }

  g_list_free_full (nodes, _BIVec_free);

  g_array_unref (x_array);
  g_array_unref (xt_array);
  g_array_unref (y_array);
  g_array_unref (yt_array);

  ncm_stats_vec_clear (&dx_stats);

  return;
}

static void
ncm_spline_new_function_spline_lnknot (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, gdouble rel_error, const gdouble f_scale)
{
  GArray *x_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *xt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *yt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GList *nodes = NULL, *wnodes = NULL;
  gsize n = ncm_spline_min_size (s);
  guint i;
  const gdouble lnxi = log (xi);
  const gdouble lnxf = log (xf);

  n = (n < 3) ? 3 : n;

  max_nodes = (max_nodes <= 0) ? G_MAXUINT64 : max_nodes;

  g_assert (xi > 0.0 && xf > xi);
  g_assert_cmpfloat (f_scale, >=, 0.0);

  g_array_set_size (xt_array, n);
  g_array_set_size (yt_array, n);

  for (i = 0; i < n; i++)
  {
    gdouble x = exp (lnxi + (lnxf - lnxi) / (n - 1.0) * i);
    gdouble y = GSL_FN_EVAL (F, x);

    BIVEC_LIST_APPEND (nodes, x, y);
    BIVEC_LIST_OK (nodes) = 0;
    g_array_append_val (x_array, x);
    g_array_append_val (y_array, y);
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
        const gdouble lnx0    = log (x0);
        const gdouble lnx1    = log (x1);
        const gdouble y0      = BIVEC_LIST_Y (wnodes);
        const gdouble y1      = BIVEC_LIST_Y (wnodes->next);
        const gdouble lnx     = (lnx0 + lnx1) / 2.0;
        const gdouble x       = exp (lnx);
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
      g_warning ("ncm_spline_new_function_spline: cannot achieve requested precision with at most %zu nodes", max_nodes);
      break;
    }

    if (improves == 0)
      break;
  }

  g_list_free_full (nodes, _BIVec_free);

  g_array_unref (x_array);
  g_array_unref (xt_array);
  g_array_unref (y_array);
  g_array_unref (yt_array);

  return;
}

static void
ncm_spline_new_function_spline_sinhknot (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble f_scale)
{
  GArray *x_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *xt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *yt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GList *nodes = NULL, *wnodes = NULL;
  gsize n = ncm_spline_min_size (s);
  guint i;
  const gdouble axi = asinh (xi);
  const gdouble axf = asinh (xf);

  g_assert_cmpfloat (f_scale, >=, 0.0);

  n = (n < 3) ? 3 : n;

  max_nodes = (max_nodes <= 0) ? G_MAXUINT64 : max_nodes;

  g_array_set_size (xt_array, n);
  g_array_set_size (yt_array, n);

  for (i = 0; i < n; i++)
  {
    gdouble x = sinh (axi + (axf - axi) / (n - 1.0) * i);
    gdouble y = GSL_FN_EVAL (F, x);

    BIVEC_LIST_APPEND (nodes, x, y);
    BIVEC_LIST_OK (nodes) = 0;

    g_array_append_val (x_array, x);
    g_array_append_val (y_array, y);
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
        const gdouble ax0     = asinh (x0);
        const gdouble ax1     = asinh (x1);
        const gdouble y0      = BIVEC_LIST_Y (wnodes);
        const gdouble y1      = BIVEC_LIST_Y (wnodes->next);
        const gdouble ax      = (ax0 + ax1) / 2.0;
        const gdouble x       = sinh (ax);
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
      g_warning ("ncm_spline_new_function_spline: cannot achieve requested precision with at most %zu nodes", max_nodes);
      break;
    }

    if (improves == 0)
      break;
  }

  g_list_free_full (nodes, _BIVec_free);

  g_array_unref (x_array);
  g_array_unref (xt_array);
  g_array_unref (y_array);
  g_array_unref (yt_array);

  return;
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
 * with scale zero and one refinement pass.
 */
void
ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error)
{
  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);

  switch (ftype)
  {
    case NCM_SPLINE_FUNCTION_SPLINE:
      ncm_spline_new_function_spline (s, F, xi, xf, max_nodes, rel_error, 0.0, 1, 1.0);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE_LNKNOT:
      ncm_spline_new_function_spline_lnknot (s, F, xi, xf, max_nodes, rel_error, 0.0);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT:
      ncm_spline_new_function_spline_sinhknot (s, F, xi, xf, max_nodes, rel_error, 0.0);
      break;
    default:
      g_assert_not_reached ();

      return;
  }
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
 * @refine: the number of refinement passes
 * @refine_ns: the number of standard deviations above the mean spacing
 *
 * Places the knots of @s on [@xi, @xf] adaptively and prepares it, see #NcmSplineFunc;
 * the absolute tolerance is @rel_error times @scale. @refine and @refine_ns are used only
 * by %NCM_SPLINE_FUNCTION_SPLINE.
 */
void
ncm_spline_set_func_scale (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble scale, const gint refine, gdouble refine_ns)
{
  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);

  switch (ftype)
  {
    case NCM_SPLINE_FUNCTION_SPLINE:
      ncm_spline_new_function_spline (s, F, xi, xf, max_nodes, rel_error, scale, refine, refine_ns);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE_LNKNOT:
      ncm_spline_new_function_spline_lnknot (s, F, xi, xf, max_nodes, rel_error, scale);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT:
      ncm_spline_new_function_spline_sinhknot (s, F, xi, xf, max_nodes, rel_error, scale);
      break;
    default:
      g_assert_not_reached ();

      return;
  }
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


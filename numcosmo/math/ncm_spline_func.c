/***************************************************************************
 *            ncm_spline_func.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_spline_func
 * @title: NcmSplineFunc
 * @short_description: Automatic generation of the knots for a spline.
 * @stability: Stable
 * @include: numcosmo/math/ncm_spline_func.h
 *
 *
 * This function implements 4 different methods to automatically determine
 * the #NcmVector of knots, $\mathbf{x}$, of a #NcmSpline given a relative error between the function $f$
 * to be interpolated and the spline result $\hat{f}$.
 *
 * All available methods start with $n_0$ knots, $\mathbf{x}_0$, distributed  across [@xi, @xf], including both limiting points. 
 * The value of $n_0$ depends on the chosen interpolation method given by @s, e.g., #NcmSplineCubicNotaknot has $n_0 = 6$.
 *
 * The function $f$ is first interpolated at the $\mathbf{x}_0$ knots, producing the interpolated function $\hat{f}_0$. 
 * Next, the existing $n_0 - 1$ bins, $\Delta \mathbf{x}_0 = \mathbf{x}_0^{i+1} - \mathbf{x}_0^{i}$, are divided in half and 
 * test points are placed at those positions, $\overline{\mathbf{x}}_0 = \frac{\mathbf{x}_0^{i+1} + \mathbf{x}_0^{i}}{2}$.
 * The following tests are done for each one of the bins $\Delta \mathbf{x}_0$ separately:
 * \begin{equation*}
 *   \left| \frac{ \hat{f}_0(\overline{\mathbf{x}}_0) - f(\overline{\mathbf{x}}_0)}{f(\overline{\mathbf{x}}_0)} \right| < \mathrm{rel \\_ error}
 * \end{equation*}
 * and
 * \begin{equation*}
 *   \left| \frac{ \int_{\Delta \mathbf{x}_0} \hat{f}_0  - \int_{\Delta \mathbf{x}_0} f }{ \int_{\Delta \mathbf{x}_0} f } \right| < \mathrm{rel \\_ error}.
 * \end{equation*}
 * Where $\int_{\Delta \mathbf{x}_0} f$ is the integral of the input function $f$ evaluated using [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule)
 *\begin{equation*}
 * \int_{\Delta \mathbf{x}_0} f = \frac{\Delta \mathbf{x}_0}{6} \left[ f \left( \mathbf{x}_0^{i} \right) + 4 f \left( \overline{\mathbf{x}}_0 \right) + f \left( \mathbf{x}_0^{i+1} \right) \right]
 *\end{equation*}
 *
 * and the interpolated function $\hat{f}_0$ is integrated by applying ncm_spline_eval_integ(). 
 * If any bin passes those relations then, its associated test point $\overline{\mathbf{x}}_0$ together with both knots, 
 * are defined as a good representation of the function $f$ and this specific bin does not need to be refined anymore.
 * If not, then $\mathbf{x}_0 \cup \overline{\mathbf{x}}_0$ is splitted once again into two more symmetric test points around $\overline{\mathbf{x}}_0$.
 * A new set of knots is defined, $\mathbf{x}_1 = \mathbf{x}_0 \cup \overline{\mathbf{x}}_0$. 
 * The ones which did not pass the tests define another set of test points $\overline{\mathbf{x}}_1$ that lie between the $\mathbf{x}_0\cup \overline{\mathbf{x}}_0$. 
 * The new set of knots $\mathbf{x}_1$ are used to create a new interpolated function $\hat{f}_1$, always in the full range [@xi, @xf]. 
 * The same tests are performed as before, but now with $\hat{f}_1(\overline{\mathbf{x}}_1)$, $f(\overline{\mathbf{x}}_1)$ 
 * and the integral now has limits $\Delta \mathbf{x}_1 = \Delta \mathbf{x}_0/2$, but only for those bins that did not pass the previous test. 
 * Therefore, the tests are always applied on three knots at a time. 
 * This procedure is repeated until the desired accuracy is met across the whole range [@xi, @xf]. Note that it will most probably create a inhomogeneous set of knots. 
 *
 *
 *
 * <inlinegraphic fileref="spline_func_knots_evolution.png" format="PNG" scale="98" align="right"/>
 *
 *
 * The figure shows a schematically evolution of the methodology for choosing the knots.
 * It starts with 6 knots, $\mathbf{x}_0$ (black filled circles in the first line), used to create the interpolated function $\hat{f}_0$.
 * The $\overline{\mathbf{x}}_0$ test points are created (blue squares in the second line) and the first tests are performed in each one of the five bins. 
 * In this example, only the first and the fourth bins did not pass both tests. 
 * One new set of knots is created, $\mathbf{x}_1 = \mathbf{x}_0 \cup \overline{\mathbf{x}}_0$ (second line) and their new test points, 
 * $\overline{\mathbf{x}}_1$ (red diamonds in the third line). 
 * Note that $\overline{\mathbf{x}}_1$ are only placed in the middle of the previously bins that did not pass the tests.
 * The tests are done 4 times, one for each bin with a red diamond at it center, 
 * with width $\Delta \mathbf{x}_1 = \Delta \mathbf{x}_0/2$. Again, only two of them passed the tests. 
 * One new set of knots is created $\mathbf{x}_2 = \mathbf{x}_1 \cup \overline{\mathbf{x}}_1$ (third line) and their new test points, 
 * $\overline{\mathbf{x}}_2$ (black vertical ticks in the fourth line). 
 * The tests are done 4 times more, one for each bin with a black vertical tick at it center, with width $\Delta \mathbf{x}_2 = \Delta \mathbf{x}_1/2$. 
 * 
 * In this schematic example, the final set of knots is given by the last line 
 * $\mathbf{x} = \mathbf{x}_3 = \mathbf{x}_2 \cup \overline{\mathbf{x}}_2$, with 19 knots in total, 
 * also showing that the final distribution is not homogeneous. 
 * It is important to note that in all steps the interpolated function is created with all its knots: 
 * $\mathbf{x}_0 \rightarrow \hat{f}_0$, $\mathbf{x}_1 \rightarrow \hat{f}_1$, $\mathbf{x}_2 \rightarrow \hat{f}_2$, $\mathbf{x}_3 \rightarrow \hat{f}_3$. 
 * It is also worth noting that, in the description and example above, it was assumed a linear distribution of knots in each step, 
 * but there are other options listed at #NcmSplineFuncType. 
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_func.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_poly.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct { gdouble x;
                 gdouble y;
                 gint ok;
} _BIVec;

#define BIVEC_LIST_APPEND(dlist, xX, yY) \
  do { \
    _BIVec *bv = g_slice_new (_BIVec); \
    bv->x = (xX); bv->y = (yY); bv->ok = FALSE; \
    dlist = g_list_append (dlist, bv); \
  } while (FALSE)
#define BIVEC_LIST_INSERT_BEFORE(dlist, node, xX, yY) \
  do { \
    _BIVec *bv = g_slice_new (_BIVec); \
    bv->x = (xX); bv->y = (yY); bv->ok = FALSE; \
    dlist = g_list_insert_before (dlist, node, bv); \
  } while (FALSE)

#define BIVEC_LIST_X(dlist) (((_BIVec *) (dlist)->data)->x)
#define BIVEC_LIST_Y(dlist) (((_BIVec *) (dlist)->data)->y)
#define BIVEC_LIST_OK(dlist) (((_BIVec *) (dlist)->data)->ok)

static void
_test_and_eval_interior_4 (GList *nodes, gsl_function *F, const gdouble yinterp1, const gdouble yinterp2, const gdouble rel_error, const gdouble f_scale, gboolean ok)
{
  gint i;
  
  GList *wnodes[3];
  gdouble poly3_dd[4];
  gdouble poly3_x[4];
  gdouble poly3_y[4];
  gdouble step, sub_step;
  gboolean go_ok = FALSE;
  
  poly3_x[0] = BIVEC_LIST_X (nodes);
  poly3_x[3] = BIVEC_LIST_X (nodes->next);
  poly3_y[0] = BIVEC_LIST_Y (nodes);
  poly3_y[3] = BIVEC_LIST_Y (nodes->next);
  
  step       = (poly3_x[3] - poly3_x[0]) / 3.0;
  sub_step   = step / 3.0;
  poly3_x[1] = poly3_x[0] + step;
  poly3_x[2] = poly3_x[1] + step;
  poly3_y[1] = GSL_FN_EVAL (F, poly3_x[1]);
  poly3_y[2] = GSL_FN_EVAL (F, poly3_x[2]);
  
  BIVEC_LIST_INSERT_BEFORE (nodes, nodes->next, poly3_x[2], poly3_y[2]);
  BIVEC_LIST_INSERT_BEFORE (nodes, nodes->next, poly3_x[1], poly3_y[1]);
  
  gsl_poly_dd_init (poly3_dd, poly3_x, poly3_y, 4);
  
  if ((fabs (poly3_y[1] - yinterp1) < ((fabs (poly3_y[1]) + f_scale) * rel_error)) &&
      (fabs (poly3_y[2] - yinterp2) < ((fabs (poly3_y[2]) + f_scale) * rel_error)))
  {
    if (ok || ((sub_step / poly3_x[0]) < NCM_SPLINE_KNOT_DIFF_TOL))
      return;
    else
      go_ok = TRUE;
  }
  
  wnodes[0] = nodes;
  wnodes[1] = wnodes[0]->next;
  wnodes[2] = wnodes[1]->next;
  
  for (i = 0; i < 3; i++)
  {
    gdouble try_yinterp1, try_yinterp2;
    
    try_yinterp1 = gsl_poly_dd_eval (poly3_dd, poly3_x, 4, poly3_x[i] + sub_step);
    try_yinterp2 = gsl_poly_dd_eval (poly3_dd, poly3_x, 4, poly3_x[i] + 2.0 * sub_step);
    _test_and_eval_interior_4 (wnodes[i], F, try_yinterp1, try_yinterp2, rel_error, f_scale, go_ok);
  }
}

static void
_BIVec_free (gpointer mem)
{
  g_slice_free (_BIVec, mem);
}

static void
ncm_spline_new_function_4 (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble f_scale)
{
  gdouble poly3_dd[4];
  gdouble poly3_x[4];
  gdouble poly3_y[4];
  GList *nodes     = NULL;
  GList *wnodes[3] = {NULL, NULL, NULL};
  gint i;
  gdouble step     = (xf - xi) / 3.0;
  gdouble sub_step = step / 3.0;
  GArray *x_array;
  GArray *y_array;
  gsize n_elem;
  
  max_nodes = (max_nodes <= 0) ? G_MAXUINT64 : max_nodes;
  g_assert_cmpfloat (f_scale, >=, 0.0);
  
  for (i = 0; i < 4; i++)
  {
    poly3_x[i] = xi + i * step;
    poly3_y[i] = GSL_FN_EVAL (F, poly3_x[i]);
    BIVEC_LIST_APPEND (nodes, poly3_x[i], poly3_y[i]);
  }
  
  gsl_poly_dd_init (poly3_dd, poly3_x, poly3_y, 4);
  
  wnodes[0] = nodes;
  wnodes[1] = wnodes[0]->next;
  wnodes[2] = wnodes[1]->next;
  
  for (i = 0; i < 3; i++)
  {
    gdouble yinterp1, yinterp2;
    
    yinterp1 = gsl_poly_dd_eval (poly3_dd, poly3_x, 4, poly3_x[i] + sub_step);
    yinterp2 = gsl_poly_dd_eval (poly3_dd, poly3_x, 4, poly3_x[i] + 2.0 * sub_step);
    
    _test_and_eval_interior_4 (wnodes[i], F, yinterp1, yinterp2, rel_error, f_scale, FALSE);
  }
  
  n_elem  = g_list_length (nodes);
  x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_elem);
  y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_elem);
  
  wnodes[0] = nodes;
  
  do {
    g_array_append_val (x_array, BIVEC_LIST_X (nodes));
    g_array_append_val (y_array, BIVEC_LIST_Y (nodes));
    nodes = nodes->next;
  } while (nodes != NULL);
  
  g_list_free_full (wnodes[0], _BIVec_free);
  
  ncm_spline_set_array (s, x_array, y_array, TRUE);
  
  return;
}

static void
ncm_spline_new_function_spline (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble f_scale)
{
  GArray *x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *xt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *yt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GList *nodes = NULL, *wnodes = NULL;
  gsize n = ncm_spline_min_size (s);
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

#define SWAP_PTR(a, b) \
  do { \
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
        const gdouble x       = (x0 + x1) / 2.0;
        const gdouble y       = GSL_FN_EVAL (F, x);
        const gdouble ys      = ncm_spline_eval (s, x);
        const gdouble Iyc     = (x1 - x0) * (y1 + y0 + 4.0 * y) / 6.0;
        const gdouble Iys     = ncm_spline_eval_integ (s, x0, x1);
        const gboolean test_p = fabs (y - ys)    < rel_error * (fabs (y)   + f_scale);
        const gboolean test_I = fabs (Iyc - Iys) < rel_error * (fabs (Iyc) + f_scale * (x1 - x0));

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
      g_warning ("ncm_spline_new_function_spline: cannot archive requested precision with at most %zu nodes", max_nodes);
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
ncm_spline_new_function_spline_lnknot (NcmSpline *s, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, gdouble rel_error, const gdouble f_scale)
{
  GArray *x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
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
  
#define SWAP_PTR(a, b) \
  do { \
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
        const gboolean test_p = fabs (y - ys)    < rel_error * (fabs (y)   + f_scale);
        const gboolean test_I = fabs (Iyc - Iys) < rel_error * (fabs (Iyc) + f_scale * (x1 - x0));

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
      g_warning ("ncm_spline_new_function_spline: cannot archive requested precision with at most %zu nodes", max_nodes);
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
  GArray *x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  GArray *y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
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
  
#define SWAP_PTR(a, b) \
  do { \
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
        const gboolean test_p = fabs (y - ys)    < rel_error * (fabs (y)   + f_scale);
        const gboolean test_I = fabs (Iyc - Iys) < rel_error * (fabs (Iyc) + f_scale * (x1 - x0));

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
      g_warning ("ncm_spline_new_function_spline: cannot achive requested precision with at most %zu nodes", max_nodes);
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
 * @ftype: a #NcmSplineFuncType
 * @F: function to be approximated by spline functions
 * @xi: lower knot
 * @xf: upper knot
 * @max_nodes: maximum number of knots
 * @rel_error: relative error between the function to be interpolated and the spline result
 *
 * This function automatically determines the knots of @s in the interval [@xi, @xf] given a @ftype and @rel_error.
 */
void
ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error)
{
  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);
  
  switch (ftype)
  {
    case NCM_SPLINE_FUNCTION_4POINTS:
      ncm_spline_new_function_4 (s, F, xi, xf, max_nodes, rel_error, 0.0);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE:
      ncm_spline_new_function_spline (s, F, xi, xf, max_nodes, rel_error, 0.0);
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
 * @ftype: a #NcmSplineFuncType
 * @F: function to be approximated by spline functions
 * @xi: lower knot
 * @xf: upper knot
 * @max_nodes: maximum number of knots
 * @rel_error: relative error between the function to be interpolated and the spline result
 * @scale: scale of function, it is used to compute the absolute tolerance abstol = f_scale * rel_error.
 *
 * This function automatically determines the knots of @s in the interval [@xi, @xf] given a @ftype and @rel_error.
 */
void
ncm_spline_set_func_scale (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, const gdouble xi, const gdouble xf, gsize max_nodes, const gdouble rel_error, const gdouble scale)
{
  ncm_assert_cmpdouble_e (xf, >, xi, DBL_EPSILON, 0.0);

  switch (ftype)
  {
    case NCM_SPLINE_FUNCTION_4POINTS:
      ncm_spline_new_function_4 (s, F, xi, xf, max_nodes, rel_error, scale);
      break;
    case NCM_SPLINE_FUNCTION_SPLINE:
      ncm_spline_new_function_spline (s, F, xi, xf, max_nodes, rel_error, scale);
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
 * @s: a #NcmSpline.
 * @ftype: a #NcmSplineFuncType
 * @F: (scope call): function to be approximated by spline functions
 * @obj: (allow-none): #GObject used by the function @F  
 * @xi: lower knot
 * @xf: upper knot
 * @max_nodes: maximum number of knots
 * @rel_error: relative error between the function to be interpolated and the spline result
 *
 * This function automatically determines the knots of @s in the interval [@xi, @xf] given a @ftype and @rel_error.
 * 
 *   
 */
void 
ncm_spline_set_func1 (NcmSpline *s, NcmSplineFuncType ftype, NcmSplineFuncF F, GObject *obj, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
{
  gsl_function gslF = {(gdouble (*) (gdouble, gpointer))F, obj};
  ncm_spline_set_func (s, ftype, &gslF, xi, xf, max_nodes, rel_error);
}


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
 * @title: Spline Autoknots
 * @short_description: Automatic generation of the knots of a spline
 *
 * This set of functions implements 4 different methods to automatically determine
 * the #NcmVector of knots of a #NcmSpline given a relative error between the function
 * to be interpolated and the spline result.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_func.h"

#include <gsl/gsl_poly.h>

typedef struct { gdouble x; gdouble y; gint ok; } _BIVec;

#define BIVEC_LIST_APPEND(dlist,xX,yY) \
do { \
_BIVec *bv = g_slice_new (_BIVec); \
bv->x = (xX); bv->y = (yY); bv->ok = FALSE; \
dlist = g_list_append (dlist, bv); \
} while(FALSE)
#define BIVEC_LIST_INSERT_BEFORE(dlist,node,xX,yY) \
do { \
_BIVec *bv = g_slice_new (_BIVec); \
bv->x = (xX); bv->y = (yY); bv->ok = FALSE; \
dlist = g_list_insert_before (dlist, node, bv); \
} while(FALSE)

#define BIVEC_LIST_X(dlist) (((_BIVec *)(dlist)->data)->x)
#define BIVEC_LIST_Y(dlist) (((_BIVec *)(dlist)->data)->y)
#define BIVEC_LIST_OK(dlist) (((_BIVec *)(dlist)->data)->ok)
#define _NCM_SPLINE_MIN_DIST GSL_DBL_EPSILON

static void
_test_and_eval_interior_4 (GList *nodes, gsl_function *F, gdouble yinterp1, gdouble yinterp2, gdouble rel_error, gboolean ok)
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

  step = (poly3_x[3] - poly3_x[0]) / 3.0;
  sub_step = step / 3.0;
  poly3_x[1] = poly3_x[0] + step;
  poly3_x[2] = poly3_x[1] + step;
  poly3_y[1] = GSL_FN_EVAL (F, poly3_x[1]);
  poly3_y[2] = GSL_FN_EVAL (F, poly3_x[2]);

  //printf ("# A (% 15.10g % 15.10g)\n", BIVEC_LIST_X (nodes), BIVEC_LIST_X (nodes->next));

  BIVEC_LIST_INSERT_BEFORE (nodes, nodes->next, poly3_x[2], poly3_y[2]);
  BIVEC_LIST_INSERT_BEFORE (nodes, nodes->next, poly3_x[1], poly3_y[1]);

  //printf ("# B (% 15.10g % 15.10g % 15.10g % 15.10g)\n", BIVEC_LIST_X (nodes), BIVEC_LIST_X (nodes->next), BIVEC_LIST_X (nodes->next->next), BIVEC_LIST_X (nodes->next->next->next));

  gsl_poly_dd_init (poly3_dd, poly3_x, poly3_y, 4);

  //printf ("# CMP % 15.10g == % 15.10g | % 15.10g == % 15.10g\n", poly3_y[1], yinterp1, poly3_y[2], yinterp2);

  //if (gsl_fcmp (poly3_y[1], yinterp1, 1e-1) == 0 && gsl_fcmp (poly3_y[2], yinterp2, 1e-1) == 0)
  if (gsl_fcmp (poly3_y[1], yinterp1, rel_error) == 0 && gsl_fcmp (poly3_y[2], yinterp2, rel_error) == 0)
  {
	if (ok || ((sub_step / poly3_x[0]) < _NCM_SPLINE_MIN_DIST))
	{
	  //printf ("# Estou dentro     ok !!! [% 15.10g, % 15.10g, % 15.10g, % 15.10g] (%d, %d)\n", poly3_x[0], poly3_x[1], poly3_x[2], poly3_x[3], ok, go_ok);
	  return;
	}
	else
	  go_ok = TRUE;
  }
  //printf ("# Estou dentro not ok !!! [% 15.10g, % 15.10g, % 15.10g, % 15.10g] (%d, %d)\n", poly3_x[0], poly3_x[1], poly3_x[2], poly3_x[3], ok, go_ok);
  wnodes[0] = nodes;
  wnodes[1] = wnodes[0]->next;
  wnodes[2] = wnodes[1]->next;

  for (i = 0; i < 3; i++)
  {
	gdouble try_yinterp1, try_yinterp2;
	try_yinterp1 = gsl_poly_dd_eval (poly3_dd, poly3_x, 4, poly3_x[i] + sub_step);
	try_yinterp2 = gsl_poly_dd_eval (poly3_dd, poly3_x, 4, poly3_x[i] + 2.0 * sub_step);
	_test_and_eval_interior_4 (wnodes[i], F, try_yinterp1, try_yinterp2, rel_error, go_ok);
  }
}

static void
_BIVec_free (gpointer mem)
{
  g_slice_free (_BIVec, mem);
}

#if (GLIB_MAJOR_VERSION < 2) || (GLIB_MINOR_VERSION < 28)
static void g_list_free_full (GList *list, GDestroyNotify free_func)
{
  GList *first = g_list_first (list);
  GList *wl = first;
  while (wl != NULL)
  {
	free_func (wl->data);
	wl = g_list_next (wl);
  }
  g_list_free (first);
}
#endif /* (GLIB_MAJOR_VERSION < 2) || (GLIB_MINOR_VERSION < 28)  */

static void
ncm_spline_new_function_4 (NcmSpline *s, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
{
	gdouble poly3_dd[4];
	gdouble poly3_x[4];
	gdouble poly3_y[4];
	GList *nodes = NULL;
	GList *wnodes[3] = {NULL, NULL, NULL};
	gint i;
	gdouble step = (xf - xi) / 3.0;
	gdouble sub_step = step / 3.0;
	GArray *x_array;
	GArray *y_array;
	gsize n_elem;

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

		_test_and_eval_interior_4 (wnodes[i], F, yinterp1, yinterp2, rel_error, FALSE);
	}

	n_elem = g_list_length (nodes);
	x_array = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), n_elem);
	y_array = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), n_elem);

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
ncm_spline_new_function_2x2 (NcmSpline *s, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
{
	gdouble poly3_dd[4];
	gdouble poly3_x[4];
	gdouble poly3_y[4];
	GList *nodes = NULL;
	GList *wnodes[3] = {NULL, NULL, NULL};
	gint i;
	gdouble step = (xf - xi) / 3.0;
	gdouble sub_step = step / 3.0;
	GArray *x_array;
	GArray *y_array;
	gsize n_elem;

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

		_test_and_eval_interior_4 (wnodes[i], F, yinterp1, yinterp2, rel_error, FALSE);
	}

	n_elem = g_list_length (nodes);
	x_array = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), n_elem);
	y_array = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), n_elem);

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
ncm_spline_new_function_spline (NcmSpline *s, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
{
	GArray *x_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
	GArray *y_array  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
	GArray *xt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
	GArray *yt_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
	GList *nodes = NULL, *wnodes = NULL;
	gsize n = ncm_spline_min_size (s);
	guint i;

	g_assert (xf > xi);
	if (max_nodes == 0)
		max_nodes = NCM_SPLINE_FUNC_DEFAULT_MAX_NODES;

	g_array_set_size (xt_array, n);
	g_array_set_size (yt_array, n);

	for (i = 0; i < n; i++)
	{
		gdouble x = xi + (xf - xi) / (n - 1.0) * i;
		gdouble y = GSL_FN_EVAL (F, x);
		BIVEC_LIST_APPEND (nodes, x, y);
		BIVEC_LIST_OK (nodes) = 0;
		g_array_append_val (x_array, x);
		g_array_append_val (y_array, y);
	}
	ncm_spline_set_array (s, x_array, y_array, TRUE);

#define TEST_CMP(a,b) ((a) != 0.0 ? fabs(((b)-(a))/(a)) : fabs((b)-(a)))
#define SWAP_PTR(a,b) do { const gpointer tmp = (b); (b) = (a); (a) = tmp; } while (FALSE)

	while (TRUE)
	{
		gsize improves = 0;
		wnodes = nodes;
		g_array_set_size (xt_array, 0);
		g_array_set_size (yt_array, 0);
		do
		{
			g_array_append_val (xt_array, BIVEC_LIST_X(wnodes));
			g_array_append_val (yt_array, BIVEC_LIST_Y(wnodes));
			//printf ("[% 20.15g % 20.15g]\t", BIVEC_LIST_X(wnodes), BIVEC_LIST_X(wnodes->next));
			if (BIVEC_LIST_OK (wnodes) == 1)
			{
				//printf ("\n");
				continue;
			}
			else
			{
				const gdouble x0 = BIVEC_LIST_X(wnodes);
				const gdouble x1 = BIVEC_LIST_X(wnodes->next);
				const gdouble y0 = BIVEC_LIST_Y(wnodes);
				const gdouble y1 = BIVEC_LIST_Y(wnodes->next);
				const gdouble x = (x0 + x1) / 2.0;
				const gdouble y = GSL_FN_EVAL (F, x);
				const gdouble ys = ncm_spline_eval (s, x);
				const gdouble Iyc = (x1 - x0) * (y1 + y0 + 4.0 * y) / 6.0;
				const gdouble Iys = ncm_spline_eval_integ (s, x0, x1);
				const gboolean test_p = (TEST_CMP(y, ys) < rel_error);
				const gboolean test_I = (TEST_CMP(Iyc, Iys) < rel_error);
#ifdef _NCM_SPLINE_TEST_DIFF
				const gdouble dyc = (y1 - y0) / (x1 - x0);
				const gdouble dys = ncm_spline_eval_deriv (s, x);
				const gboolean test_d = (TEST_CMP(dyc, dys) < rel_error);
#endif /* _NCM_SPLINE_TEST_DIFF */

				if (fabs ((x-x0)/x) < NCM_SPLINE_KNOT_DIFF_TOL)
					g_error ("Tolerance of the difference between knots was reached. Interpolated function is probably discontinuous at % 20.15g.", x);

				BIVEC_LIST_INSERT_BEFORE (nodes, wnodes->next, x, y);
				wnodes = g_list_next (wnodes);
				g_array_append_val(xt_array, BIVEC_LIST_X(wnodes));
				g_array_append_val(yt_array, BIVEC_LIST_Y(wnodes));
				BIVEC_LIST_OK(wnodes) = BIVEC_LIST_OK(wnodes->prev);

				//printf ("improves %03zd (% 10.8g) [% 10.8g % 10.8g] [% 10.8g % 10.8g] ", improves, x, y, ys, Iyc, Iys);

#ifdef _NCM_SPLINE_TEST_DIFF
				if (test_p && test_I && test_d)
#else
					if (test_p && test_I)
#endif /* _NCM_SPLINE_TEST_DIFF */
				{
					BIVEC_LIST_OK(wnodes->prev)++;
					BIVEC_LIST_OK(wnodes)++;
				}
				else
					improves++;
				//printf ("OK %d %d\n", BIVEC_LIST_OK(wnodes->prev), BIVEC_LIST_OK(wnodes));
			}
		} while ((wnodes = g_list_next (wnodes)) && wnodes->next);
		g_array_append_val(xt_array, BIVEC_LIST_X(wnodes));
		g_array_append_val(yt_array, BIVEC_LIST_Y(wnodes));

		SWAP_PTR (x_array, xt_array);
		SWAP_PTR (y_array, yt_array);

		if (x_array->len > max_nodes)
			g_error ("ncm_spline_new_function_spline: cannot achive requested precision with at most %zu nodes", max_nodes);

		ncm_spline_set_array (s, x_array, y_array, TRUE);
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
ncm_spline_new_function_spline_lnknot (NcmSpline *s, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
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

  g_assert (xi > 0.0 && xf > xi);

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
	printf ("%d %g %g\n", i, x, y);
  }
  ncm_spline_set_array (s, x_array, y_array, TRUE);

#define TEST_CMP(a,b) ((a) != 0.0 ? fabs(((b)-(a))/(a)) : fabs((b)-(a)))
#define SWAP_PTR(a,b) do { const gpointer tmp = (b); (b) = (a); (a) = tmp; } while (FALSE)

  while (TRUE)
  {
	gsize improves = 0;
	wnodes = nodes;
	g_array_set_size (xt_array, 0);
	g_array_set_size (yt_array, 0);
	do
	{
	  g_array_append_val(xt_array, BIVEC_LIST_X(wnodes));
	  g_array_append_val(yt_array, BIVEC_LIST_Y(wnodes));
	  //printf ("[% 20.15g % 20.15g]\t", BIVEC_LIST_X(wnodes), BIVEC_LIST_X(wnodes->next));
	  if (BIVEC_LIST_OK (wnodes) == 1)
	  {
		//printf ("\n");
		continue;
	  }
	  else
	  {
		const gdouble x0 = BIVEC_LIST_X(wnodes);
		const gdouble x1 = BIVEC_LIST_X(wnodes->next);
		const gdouble lnx0 = log (x0);
		const gdouble lnx1 = log (x1);
		const gdouble y0 = BIVEC_LIST_Y(wnodes);
		const gdouble y1 = BIVEC_LIST_Y(wnodes->next);
		const gdouble lnx = (log(x0) + log(x1)) / 2.0;
		const gdouble x = exp (lnx);
		const gdouble y = GSL_FN_EVAL (F, x);
		const gdouble ys = ncm_spline_eval (s, x);
		const gdouble Iyc = (lnx1 - lnx0) * (x1 * y1 + x0 * y0 + 4.0 * x * y) / 6.0;
		const gdouble Iys = ncm_spline_eval_integ (s, x0, x1);
		const gboolean test_p = (TEST_CMP(y, ys) < rel_error);
		const gboolean test_I = (TEST_CMP(Iyc, Iys) < rel_error);
#ifdef _NCM_SPLINE_TEST_DIFF
		const gdouble dyc = (y1 - y0) / (x1 - x0);
		const gdouble dys = ncm_spline_eval_deriv (s, x);
		const gboolean test_d = (TEST_CMP(dyc, dys) < rel_error);
#endif /* _NCM_SPLINE_TEST_DIFF */

		BIVEC_LIST_INSERT_BEFORE (nodes, wnodes->next, x, y);
		wnodes = g_list_next (wnodes);
		g_array_append_val(xt_array, BIVEC_LIST_X(wnodes));
		g_array_append_val(yt_array, BIVEC_LIST_Y(wnodes));
		BIVEC_LIST_OK(wnodes) = BIVEC_LIST_OK(wnodes->prev);

		//printf ("improves %03zd (% 10.8g) [% 10.8g % 10.8g] [% 10.8g % 10.8g] ", improves, x, y, ys, Iyc, Iys);

#ifdef _NCM_SPLINE_TEST_DIFF
		if (test_p && test_I && test_d)
#else
		  if (test_p && test_I)
#endif /* _NCM_SPLINE_TEST_DIFF */
		{
		  BIVEC_LIST_OK(wnodes->prev)++;
		  BIVEC_LIST_OK(wnodes)++;
		}
		else
		  improves++;
		//printf ("OK %d %d\n", BIVEC_LIST_OK(wnodes->prev), BIVEC_LIST_OK(wnodes));
	  }
	} while ((wnodes = g_list_next (wnodes)) && wnodes->next);
	g_array_append_val(xt_array, BIVEC_LIST_X(wnodes));
	g_array_append_val(yt_array, BIVEC_LIST_Y(wnodes));

	SWAP_PTR (x_array, xt_array);
	SWAP_PTR (y_array, yt_array);

	ncm_spline_set_array (s, x_array, y_array, TRUE);
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
 * @s: a #NcmSpline.
 * @ftype: a #NcmSplineFuncType.
 * @F: function to be approximated by spline functions.
 * @xi: lower knot.
 * @xf: upper knot.
 * @max_nodes: maximum number of knots.
 * @rel_error: relative error between the function to be interpolated and the spline result.
 *
 * This function automatically determines the knots of @s in the interval [@xi, @xf] given a @ftype and @rel_error.
   */
void
ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error)
{
  switch (ftype)
  {
	case NCM_SPLINE_FUNCTION_4POINTS:
	  ncm_spline_new_function_4 (s, F, xi, xf, max_nodes, rel_error);
	  break;
	case NCM_SPLINE_FUNCTION_2x2POINTS:
	  g_assert_not_reached ();
	  ncm_spline_new_function_2x2 (s, F, xi, xf, max_nodes, rel_error);
	  break;
	case NCM_SPLINE_FUNCTION_SPLINE:
	  ncm_spline_new_function_spline (s, F, xi, xf, max_nodes, rel_error);
	  break;
	case NCM_SPLINE_FUNCTION_SPLINE_LNKNOT:
	  ncm_spline_new_function_spline_lnknot (s, F, xi, xf, max_nodes, rel_error);
	  break;
	default:
	  g_assert_not_reached ();
	  return;
  }
}

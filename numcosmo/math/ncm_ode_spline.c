/***************************************************************************
 *            ncm_ode_spline.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:ncm_ode_spline
 * @title: Spline
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>

typedef struct _NcmOdeSplineDydxData NcmOdeSplineDydxData;

/**
 * NcFunctionParams:
 *
 * FIXME
 */
struct _NcmOdeSplineDydxData
{
  NcmOdeSplineDydx dydx;
  gpointer userdata;
};

static gint
_ncm_ode_spline_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmOdeSplineDydxData *dydx_data = (NcmOdeSplineDydxData *) f_data;
  NV_Ith_S (ydot, 0) = dydx_data->dydx (NV_Ith_S (y, 0), x, dydx_data->userdata);
  return 0;
}

/**
 * ncm_ode_spline_new: (skip)
 * @s: a #NcmSpline
 * @dydx: a #NcmOdeSplineDydx
 * @userdata: FIXME
 * @yi: FIXME
 * @xi: FIXME
 * @xf: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmOdeSpline *
ncm_ode_spline_new (NcmSpline *s, NcmOdeSplineDydx dydx, gpointer userdata, gdouble yi, gdouble xi, gdouble xf)
{
  NcmOdeSpline *os = g_slice_new (NcmOdeSpline);
  N_Vector y = N_VNew_Serial(1);
  os->y = y;
  os->dydx = dydx;
  os->cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  os->y_array = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
  os->x_array = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
  os->yi = yi;
  os->xi = xi;
  os->xf = xf;

  NV_Ith_S(y, 0) = yi;
  CVodeInit (os->cvode, &_ncm_ode_spline_f, xi, y);
  os->s = ncm_spline_ref (s);

  os->ctrl = ncm_model_ctrl_new (NULL);

  return os;
}

/**
 * ncm_ode_spline_prepare:
 * @os: a #NcmOdeSpline
 * @userdata: FIXME
 *
 * FIXME
 */
void
ncm_ode_spline_prepare (NcmOdeSpline *os, gpointer userdata)
{
  NcmOdeSplineDydxData f_data = {os->dydx, userdata};
  gdouble x0;

  NV_Ith_S (os->y, 0) = os->yi;
  CVodeReInit (os->cvode, os->xi, os->y);

  g_array_set_size (os->x_array, 0);
  g_array_set_size (os->y_array, 0);

  g_array_append_val (os->x_array, os->xi);
  g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));

  CVodeSStolerances (os->cvode, NC_INT_ERROR, 1e-80); /* FIXME */
  CVodeSetMaxNumSteps (os->cvode, NC_INT_PARTITION);
  CVodeSetStopTime (os->cvode, os->xf);
  CVodeSetUserData (os->cvode, &f_data);

  while(1)
  {
	CVode (os->cvode, os->xf, os->y, &x0, CV_ONE_STEP);
	g_array_append_val (os->x_array, x0);
	g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));
	if (x0 == os->xf)
	  break;
  }

  ncm_spline_set_array (os->s, os->x_array, os->y_array, TRUE);
}

/**
 * ncm_ode_spline_free:
 * @os: a #NcmOdeSpline
 *
 * FIXME
 */
void
ncm_ode_spline_free (NcmOdeSpline *os)
{
  ncm_spline_free (os->s);
  g_array_unref (os->x_array);
  g_array_unref (os->y_array);

  CVodeFree (&os->cvode);
  N_VDestroy (os->y);

  ncm_model_ctrl_free (os->ctrl);

  g_slice_free (NcmOdeSpline, os);
}

/* Spline with a not a knot boundary conditions */
/* Got from gsl 1.15 and then adapted */

typedef struct __NcmSpline_state_t
{
  GArray *b;
  GArray *c;
  GArray *d;
  GArray *g;
  GArray *diag;
  GArray *offdiag;
} _NcmSpline_state_t;


static void *
_ncm_spline_gsl_notaknot_alloc (size_t size)
{
  _NcmSpline_state_t *state = g_slice_new (_NcmSpline_state_t);

  state->b = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), size);
  g_array_set_size (state->b, size);

  state->c = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), size);
  g_array_set_size (state->c, size);

  state->d = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), size);
  g_array_set_size (state->d, size);

  state->g = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), size);
  g_array_set_size (state->g, size);

  state->diag = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), size);
  g_array_set_size (state->diag, size);

  state->offdiag = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), size);
  g_array_set_size (state->offdiag, size);

  return state;
}

static gint
_ncm_spline_gsl_notaknot_init (void *vstate, const gdouble xa[], const gdouble ya[], size_t size)
{
  _NcmSpline_state_t *state = (_NcmSpline_state_t *) vstate;
  const size_t n = size - 1;
  const size_t sys_size = size - 2;
  const size_t nm1 = n - 1;
  const size_t nm2 = nm1 - 1;
  const size_t nm3 = nm2 - 1;
  size_t start_i = 0, pad_i = 0;
  size_t i;

  g_assert (sys_size > 1);

  for (i = 1; i <= nm3; i++)
  {
	const gdouble h_i   = xa[i + 1] - xa[i];
	const gdouble h_ip1 = xa[i + 2] - xa[i + 1];
	const gdouble ydiff_i   = ya[i + 1] - ya[i];
	const gdouble ydiff_ip1 = ya[i + 2] - ya[i + 1];
	const gdouble g_i = 1.0 / h_i;
	const gdouble g_ip1 = 1.0 / h_ip1;

	g_array_index(state->offdiag, gdouble, i) = h_ip1;
	g_array_index(state->diag, gdouble, i)    = 2.0 * (h_ip1 + h_i);
	g_array_index(state->g, gdouble, i)       = 3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i);
  }

  {
	const gdouble h_0 = xa[1] - xa[0];
	const gdouble h_1 = xa[2] - xa[1];
	const gdouble ydiff_0 = ya[1] - ya[0];
	const gdouble ydiff_1 = ya[2] - ya[1];
	if (fabs((h_0 - h_1) / h_0) < GSL_SQRT_DBL_EPSILON)
	{
	  start_i = 1;
	  g_array_index(state->c, gdouble, 1) = 3.0 * (ydiff_1 -  ydiff_0 * h_1 / h_0) / ((h_1 + h_0) * (2.0 * h_1 + h_0));
	  g_array_index(state->g, gdouble, 1) -= h_1 * g_array_index(state->c, gdouble, 1);
	}
	else
	{
	  g_array_index(state->offdiag, gdouble, 0) = h_1;
	  g_array_index(state->diag, gdouble, 0)    = h_1 * (2.0 * h_1 + h_0) / (h_1 - h_0);
	  g_array_index(state->g, gdouble, 0)       = 3.0 * h_1 * (ydiff_1 -  ydiff_0 * h_1 / h_0) / (h_1 * h_1 - h_0 * h_0);
	}
  }

  {
	const gdouble h_nm1   = xa[n] - xa[nm1];
	const gdouble h_nm2 = xa[nm1] - xa[nm2];
	const gdouble ydiff_nm1   = ya[n] - ya[nm1];
	const gdouble ydiff_nm2 = ya[nm1] - ya[nm2];
	if (fabs((h_nm2 - h_nm1) / h_nm2) < GSL_SQRT_DBL_EPSILON)
	{
	  pad_i = 1;
	  g_array_index(state->c, gdouble, nm1)      = 3.0 * (ydiff_nm1 * h_nm2 / h_nm1 - ydiff_nm2) / ((h_nm2 + h_nm1) * (2.0 * h_nm2 + h_nm1));
	  g_array_index(state->g, gdouble, nm1 - 2) -= h_nm2 * g_array_index (state->c, gdouble, nm1);
	}
	else
	{
	  g_array_index(state->offdiag, gdouble, nm2) = h_nm1;
	  g_array_index(state->diag, gdouble, nm2)    = h_nm2 * (2.0 * h_nm2 + h_nm1) / (h_nm2 - h_nm1);
	  g_array_index(state->g, gdouble, nm2)       = 3.0 * h_nm2 * (ydiff_nm1 * h_nm2 / h_nm1 - ydiff_nm2) / (h_nm2 * h_nm2 - h_nm1 * h_nm1);
	}
  }

  {
	gsize loc_sys_size = sys_size - start_i - pad_i;
	gsl_vector_view g_vec = gsl_vector_view_array(&g_array_index(state->g, gdouble, start_i), loc_sys_size);
	gsl_vector_view diag_vec = gsl_vector_view_array(&g_array_index(state->diag, gdouble, start_i), loc_sys_size);
	gsl_vector_view offdiag_vec = gsl_vector_view_array(&g_array_index(state->offdiag, gdouble, start_i), loc_sys_size - 1);
	gsl_vector_view solution_vec = gsl_vector_view_array (&g_array_index(state->c, gdouble, start_i + 1), loc_sys_size);

	gint status = gsl_linalg_solve_symm_tridiag (&diag_vec.vector,
	                                             &offdiag_vec.vector,
	                                             &g_vec.vector,
	                                             &solution_vec.vector);

	{
	  const gdouble h_0 = xa[1] - xa[0];
	  const gdouble h_1 = xa[2] - xa[1];
	  const gdouble h_nm1 = xa[n] - xa[nm1];
	  const gdouble h_nm2 = xa[nm1] - xa[nm2];
	  const gdouble c1   = g_array_index(state->c, gdouble, 1);
	  const gdouble c2   = g_array_index(state->c, gdouble, 2);
	  const gdouble cnm1 = g_array_index(state->c, gdouble, nm1);
	  const gdouble cnm2 = g_array_index(state->c, gdouble, nm2);

	  g_array_index (state->c, gdouble, 0) = c1 + h_0 * (c1 - c2) / h_1;
	  g_array_index (state->c, gdouble, n) = cnm1 + h_nm1 * (cnm1 - cnm2) / h_nm2;
	}

	for (i = 0; i < n; i++)
	{
	  const gdouble dx = xa[i + 1] - xa[i];
	  const gdouble dy = ya[i + 1] - ya[i];
	  const gdouble c_ip1 = g_array_index (state->c, gdouble, i + 1);
	  const gdouble c_i = g_array_index (state->c, gdouble, i);

	  g_array_index (state->b, gdouble, i) = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
	  g_array_index (state->d, gdouble, i) = (c_ip1 - c_i) / (3.0 * dx);
	}

	return status;
  }
}

static void
_ncm_spline_gsl_notaknot_free (void *vstate)
{
  _NcmSpline_state_t *state = (_NcmSpline_state_t *) vstate;
  g_array_free (state->b, TRUE);
  g_array_free (state->c, TRUE);
  g_array_free (state->d, TRUE);
  g_array_free (state->g, TRUE);
  g_array_free (state->diag, TRUE);
  g_array_free (state->offdiag, TRUE);
  g_slice_free (_NcmSpline_state_t, state);
}

/* function for common coefficient determination
 *
 static inline void
 _ncm_spline_notaknot_coeff_calc (const gdouble c_array[], gdouble dy, gdouble dx, size_t index, gdouble *b, gdouble *c, gdouble *d)
 {
   const gdouble c_i = c_array[index];
   const gdouble c_ip1 = c_array[index + 1];
   *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
   *c = c_i;
   *d = (c_ip1 - c_i) / (3.0 * dx);
   }*/

static gint
_ncm_spline_gsl_notaknot_eval (const void *vstate, const gdouble x_array[], const gdouble y_array[], size_t size, gdouble x, gsl_interp_accel *a, gdouble *y)
{
  const _NcmSpline_state_t *state = (const _NcmSpline_state_t *) vstate;

  size_t i;

  if (a != NULL)
	i = gsl_interp_accel_find (a, x_array, size, x);
  else
	i = gsl_interp_bsearch (x_array, x, 0, size - 1);

  {
	const gdouble delx = x - x_array[i];
	const gdouble b_i = g_array_index (state->b, gdouble, i);
	const gdouble c_i = g_array_index (state->c, gdouble, i);
	const gdouble d_i = g_array_index (state->d, gdouble, i);

	*y = y_array[i] + delx * (b_i + delx * (c_i + delx * d_i));
	return GSL_SUCCESS;
  }
}

static gint
_ncm_spline_gsl_notaknot_deriv (const void *vstate, const gdouble x_array[], const gdouble y_array[], size_t size, gdouble x, gsl_interp_accel * a, gdouble *dydx)
{
  const _NcmSpline_state_t *state = (const _NcmSpline_state_t *) vstate;
  size_t i;

  if (a != NULL)
	i = gsl_interp_accel_find (a, x_array, size, x);
  else
	i = gsl_interp_bsearch (x_array, x, 0, size - 1);

  {
	const gdouble delx = x - x_array[i];
	const gdouble b_i = g_array_index (state->b, gdouble, i);
	const gdouble c_i = g_array_index (state->c, gdouble, i);
	const gdouble d_i = g_array_index (state->d, gdouble, i);
	*dydx = b_i + delx * (2.0 * c_i + 3.0 * d_i * delx);
	return GSL_SUCCESS;
  }
}

static gint
_ncm_spline_gsl_notaknot_deriv2 (const void *vstate, const gdouble x_array[], const gdouble y_array[], size_t size, gdouble x, gsl_interp_accel *a, gdouble *y_pp)
{
  const _NcmSpline_state_t *state = (const _NcmSpline_state_t *) vstate;
  size_t i;

  if (a != 0)
	i = gsl_interp_accel_find (a, x_array, size, x);
  else
	i = gsl_interp_bsearch (x_array, x, 0, size - 1);

  {
	const gdouble delx = x - x_array[i];
	const gdouble c_i = g_array_index (state->c, gdouble, i);
	const gdouble d_i = g_array_index (state->d, gdouble, i);
	*y_pp = 2.0 * c_i + 6.0 * d_i * delx;
	return GSL_SUCCESS;
  }
}

static gint
_ncm_spline_gsl_notaknot_integ (const void *vstate, const gdouble x_array[], const gdouble y_array[], size_t size, gsl_interp_accel *acc, gdouble a, gdouble b, gdouble *result)
{
  const _NcmSpline_state_t *state = (const _NcmSpline_state_t *) vstate;

  size_t i, index_a, index_b;

  if (acc != NULL)
  {
	index_a = gsl_interp_accel_find (acc, x_array, size, a);
	index_b = gsl_interp_accel_find (acc, x_array, size, b);
  }
  else
  {
	index_a = gsl_interp_bsearch (x_array, a, 0, size - 1);
	index_b = gsl_interp_bsearch (x_array, b, 0, size - 1);
  }

  *result = 0.0;

  for(i = index_a; i <= index_b; i++)
  {
	const gdouble x_hi = x_array[i + 1];
	const gdouble x_lo = x_array[i];
	const gdouble y_lo = y_array[i];

	{
	  const gdouble b_i = g_array_index (state->b, gdouble, i);
	  const gdouble c_i = g_array_index (state->c, gdouble, i);
	  const gdouble d_i = g_array_index (state->d, gdouble, i);

	  if (i == index_a || i == index_b)
	  {
		gdouble x1 = (i == index_a) ? a : x_lo;
		gdouble x2 = (i == index_b) ? b : x_hi;
		*result += _ncm_spline_util_integ_eval (y_lo, b_i, c_i, d_i, x_lo, x1, x2);
	  }
	  else
	  {
		const gdouble dx = x_hi - x_lo;
		*result += dx * (y_lo + dx * (0.5 * b_i + dx * (c_i / 3.0 + 0.25 * d_i * dx)));
	  }
	}
  }

  return GSL_SUCCESS;
}

static const gsl_interp_type _ncm_spline_notaknot_type =
{
  "cspline_notaknot",
  6,
  &_ncm_spline_gsl_notaknot_alloc,
  &_ncm_spline_gsl_notaknot_init,
  &_ncm_spline_gsl_notaknot_eval,
  &_ncm_spline_gsl_notaknot_deriv,
  &_ncm_spline_gsl_notaknot_deriv2,
  &_ncm_spline_gsl_notaknot_integ,
  &_ncm_spline_gsl_notaknot_free
};

const gsl_interp_type * gsl_interp_nc_spline_notaknot = &_ncm_spline_notaknot_type;

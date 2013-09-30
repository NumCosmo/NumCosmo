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
 * @title: ODE Spline Interpolation
 * @short_description: Automatic generation of splines from ODE solutions
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_ode_spline.h"
#include "math/cvode_util.h"
#include "math/integral.h"

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <nvector/nvector_serial.h> 
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
  N_Vector y = N_VNew_Serial (1);

  NCM_UNUSED (userdata);
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
  gdouble x, x0;
  gint i = 0;

  NV_Ith_S (os->y, 0) = os->yi;
  CVodeReInit (os->cvode, os->xi, os->y);

  g_array_set_size (os->x_array, 0);
  g_array_set_size (os->y_array, 0);

  g_array_append_val (os->x_array, os->xi);
  g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));

  CVodeSStolerances (os->cvode, NCM_INTEGRAL_ERROR, 1e-80); /* FIXME */
  CVodeSetMaxNumSteps (os->cvode, NCM_INTEGRAL_PARTITION);
  CVodeSetStopTime (os->cvode, os->xf);
  CVodeSetUserData (os->cvode, &f_data);

  x0 = os->xi;
  while (TRUE)
  {
    gint flag = CVode (os->cvode, os->xf, os->y, &x, CV_ONE_STEP);
    if (!ncm_cvode_util_check_flag (&flag, "CVode", 1))
      g_error ("ncm_ode_spline_prepare: cannot integrate function.");
    if (x0 == x && i++ > 10)
      g_error ("ncm_ode_spline_prepare: cannot integrate function.");
    x0 = x;
    g_array_append_val (os->x_array, x);
    g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));
    if (x == os->xf)
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
  ncm_spline_clear (&os->s);
  if (os->x_array != NULL)
    g_array_unref (os->x_array);
  if (os->y_array != NULL)
    g_array_unref (os->y_array);

  if (os->cvode != NULL)
    CVodeFree (&os->cvode);
  if (os->y != NULL)
    N_VDestroy (os->y);

  ncm_model_ctrl_clear (&os->ctrl);

  g_slice_free (NcmOdeSpline, os);
}

/**
 * ncm_ode_spline_clear:
 * @os: a #NcmOdeSpline
 *
 * FIXME
 */
void
ncm_ode_spline_clear (NcmOdeSpline **os)
{
  g_clear_pointer (os, &ncm_ode_spline_free);
}

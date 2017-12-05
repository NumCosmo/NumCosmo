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
 * @title: NcmOdeSpline
 * @short_description: Automatic generation of splines from ODE solvers.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_ode_spline.h"
#include "math/integral.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_XI,
  PROP_XF,
  PROP_YI,
  PROP_YF,
  PROP_DYDX,
  PROP_SPLINE,
  PROP_STOP_HNIL,
};

G_DEFINE_TYPE (NcmOdeSpline, ncm_ode_spline, G_TYPE_OBJECT);

static void
ncm_ode_spline_init (NcmOdeSpline *os)
{
  os->cvode      = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
  os->cvode_init = FALSE;
  os->y          = N_VNew_Serial (1);
  os->y_array    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  os->x_array    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  os->xi         = GSL_NAN;
  os->xf         = GSL_NAN;
  os->yi         = GSL_NAN;
	os->yf         = GSL_NAN;
  os->reltol     = 0.0;
  os->abstol     = 0.0;
  os->dydx       = NULL;
  os->s          = NULL;
  os->s_init     = FALSE;
  os->hnil       = FALSE;
  os->stop_hnil  = FALSE;
  os->ctrl       = ncm_model_ctrl_new (NULL);  
}

static void
_ncm_ode_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmOdeSpline *os = NCM_ODE_SPLINE (object);
  g_return_if_fail (NCM_IS_ODE_SPLINE (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      ncm_ode_spline_set_reltol (os, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_ode_spline_set_abstol (os, g_value_get_double (value));
      break;
    case PROP_XI:
      ncm_ode_spline_set_xi (os, g_value_get_double (value));
      break;
    case PROP_XF:
      ncm_ode_spline_set_xf (os, g_value_get_double (value));
      break;
    case PROP_YI:
      ncm_ode_spline_set_yi (os, g_value_get_double (value));
      break;
    case PROP_YF:
      ncm_ode_spline_set_yf (os, g_value_get_double (value));
      break;
    case PROP_DYDX:
      os->dydx = (NcmOdeSplineDydx) g_value_get_pointer (value);
      break;
    case PROP_SPLINE:
      os->s = g_value_dup_object (value);
      break;
    case PROP_STOP_HNIL:
      os->stop_hnil = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_ode_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmOdeSpline *os = NCM_ODE_SPLINE (object);
  g_return_if_fail (NCM_IS_ODE_SPLINE (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, os->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, os->abstol);
      break;
    case PROP_XI:
      g_value_set_double (value, os->xi);
      break;
    case PROP_XF:
      g_value_set_double (value, os->xf);
      break;
    case PROP_YI:
      g_value_set_double (value, os->yi);
      break;
    case PROP_YF:
      g_value_set_double (value, os->yf);
      break;
    case PROP_DYDX:
      g_value_set_pointer (value, os->dydx);
      break;
    case PROP_SPLINE:
      g_value_set_object (value, os->s);
      break;
    case PROP_STOP_HNIL:
      g_value_set_boolean (value, os->stop_hnil);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_ode_spline_dispose (GObject *object)
{
  NcmOdeSpline *os = NCM_ODE_SPLINE (object);

  ncm_spline_clear (&os->s);
  g_clear_pointer (&os->x_array, g_array_unref);
  g_clear_pointer (&os->y_array, g_array_unref);  
  ncm_model_ctrl_clear (&os->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_ode_spline_parent_class)->dispose (object);
}

static void
_ncm_ode_spline_finalize (GObject *object)
{
  NcmOdeSpline *os = NCM_ODE_SPLINE (object);

  if (os->cvode != NULL)
  {
    CVodeFree (&os->cvode);
    os->cvode = NULL;
  }
  if (os->y != NULL)
  {
    N_VDestroy (os->y);
    os->y = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_ode_spline_parent_class)->finalize (object);
}

static void
ncm_ode_spline_class_init (NcmOdeSplineClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_ode_spline_set_property;
  object_class->get_property = &_ncm_ode_spline_get_property;
  object_class->dispose      = &_ncm_ode_spline_dispose;
  object_class->finalize     = &_ncm_ode_spline_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NCM_ODE_SPLINE_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abs",
                                                        NULL,
                                                        "Absolute tolerance",
                                                        0.0, 1.0, NCM_ODE_SPLINE_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_XI,
                                   g_param_spec_double ("xi",
                                                        NULL,
                                                        "Initial point",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_XF,
                                   g_param_spec_double ("xf",
                                                        NULL,
                                                        "Final point",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_YI,
                                   g_param_spec_double ("yi",
                                                        NULL,
                                                        "Initial Value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
                                   PROP_YF,
                                   g_param_spec_double ("yf",
                                                        NULL,
                                                        "Final Value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DYDX,
                                   g_param_spec_pointer ("dydx",
                                                         NULL,
                                                         "Pointer to the dydx function",
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SPLINE,
                                   g_param_spec_object ("spline",
                                                        NULL,
                                                        "Spline algorithm to be used",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_STOP_HNIL,
                                   g_param_spec_boolean ("stop-hnil",
                                                         NULL,
                                                         "Whether treat hnil as error",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

typedef struct _NcmOdeSplineDydxData NcmOdeSplineDydxData;

/**
 * NcFunctionParams:
 *
 * FIXME
 */
struct _NcmOdeSplineDydxData
{
	NcmOdeSpline *os;
  gpointer userdata;
};

static gint
_ncm_ode_spline_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmOdeSplineDydxData *dydx_data = (NcmOdeSplineDydxData *) f_data;
  NV_Ith_S (ydot, 0) = dydx_data->os->dydx (NV_Ith_S (y, 0), x, dydx_data->userdata);
  /*printf ("% 20.15g % 20.15g % 20.15g\n", x, NV_Ith_S (y, 0), NV_Ith_S (ydot, 0));*/
  return 0;
}

/**
 * ncm_ode_spline_new:
 * @s: a #NcmSpline
 * @dydx: (scope notified): a #NcmOdeSplineDydx
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmOdeSpline *
ncm_ode_spline_new (NcmSpline *s, NcmOdeSplineDydx dydx)
{
  NcmOdeSpline *os = g_object_new (NCM_TYPE_ODE_SPLINE,
                                   "dydx", dydx,
                                   "spline", s,
                                   NULL);
  return os;
}

/**
 * ncm_ode_spline_new_full:
 * @s: a #NcmSpline
 * @dydx: (scope notified): a #NcmOdeSplineDydx
 * @yi: FIXME
 * @xi: FIXME
 * @xf: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmOdeSpline *
ncm_ode_spline_new_full (NcmSpline *s, NcmOdeSplineDydx dydx, gdouble yi, gdouble xi, gdouble xf)
{
  NcmOdeSpline *os = g_object_new (NCM_TYPE_ODE_SPLINE,
                                   "xi", xi,
                                   "xf", xf,
                                   "yi", yi,
                                   "dydx", dydx,
                                   "spline", s,
                                   NULL);
  return os;
}

static gint 
_ncm_ode_spline_yf_root (realtype lambda, N_Vector y, realtype *gout, gpointer user_data)
{
  NcmOdeSplineDydxData *dydx_data = (NcmOdeSplineDydxData *) user_data;

	gout[0] = (NV_Ith_S (y, 0) - dydx_data->os->yf);

	return 0;
} 

/**
 * ncm_ode_spline_prepare:
 * @os: a #NcmOdeSpline
 * @userdata: (closure): FIXME
 *
 * FIXME
 */
void
ncm_ode_spline_prepare (NcmOdeSpline *os, gpointer userdata)
{
  NcmOdeSplineDydxData f_data = {os, userdata};
  gdouble x, x0;
  gint flag;

  NV_Ith_S (os->y, 0) = os->yi;
  
  if (!os->cvode_init)
  {
    flag = CVodeInit (os->cvode, &_ncm_ode_spline_f, os->xi, os->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    
    os->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (os->cvode, os->xi, os->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  g_array_set_size (os->x_array, 0);
  g_array_set_size (os->y_array, 0);

  g_array_append_val (os->x_array, os->xi);
  g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));

  flag = CVodeSStolerances (os->cvode, os->reltol, os->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );
  
  flag = CVodeSetMaxNumSteps (os->cvode, NCM_INTEGRAL_PARTITION);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
  
  flag = CVodeSetUserData (os->cvode, &f_data);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  x0 = os->xi;
  if (!gsl_finite (os->dydx (NV_Ith_S (os->y, 0), x0, f_data.userdata)))
    g_error ("ncm_ode_spline_prepare: not finite integrand.");

	if (!gsl_finite (os->yf))
	{
		g_assert (gsl_finite (os->xf));
		
		flag = CVodeSetStopTime (os->cvode, os->xf);
		NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

		os->hnil = FALSE;
		while (TRUE)
		{
			flag = CVode (os->cvode, os->xf, os->y, &x, CV_ONE_STEP);
			NCM_CVODE_CHECK (&flag, "ncm_ode_spline_prepare[CVode]", 1, );

			if (G_UNLIKELY (os->hnil))
			{
				if (os->stop_hnil)
				{
					g_error ("ncm_ode_spline_prepare: cannot integrate function %d.", flag);
				}
				else
					break;
			}

			/*if (((x0 != 0.0) ? fabs ((x - x0) / x0) : fabs (x - x0)) > NCM_ODE_SPLINE_MIN_STEP)*/
			if (x > x0 + fabs (x0) * NCM_ODE_SPLINE_MIN_STEP)
			{
				g_array_append_val (os->x_array, x);
				g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));
				/*printf ("% 20.15g % 20.15g\n", x, NV_Ith_S (os->y, 0));*/
				x0 = x;
			}

			if (x == os->xf)
				break;
		}
	}
	else
	{
		const gdouble xf = (os->xi != 0.0) ? os->xi * 2.0 : 1.0;
		gdouble last_y = GSL_NEGINF;
		
		flag = CVodeRootInit (os->cvode, 1, &_ncm_ode_spline_yf_root);
		NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );    
		
		os->hnil = FALSE;
		while (TRUE)
		{
			flag = CVode (os->cvode, xf, os->y, &x, CV_ONE_STEP);
			NCM_CVODE_CHECK (&flag, "ncm_ode_spline_prepare[CVode]", 1, );

			if (G_UNLIKELY (os->hnil))
			{
				if (os->stop_hnil)
				{
					g_error ("ncm_ode_spline_prepare: cannot integrate function %d.", flag);
				}
				else
					break;
			}

			if (x > x0 + fabs (x0) * NCM_ODE_SPLINE_MIN_STEP)
			{
				g_array_append_val (os->x_array, x);
				g_array_append_val (os->y_array, NV_Ith_S (os->y, 0));
				x0 = x;
			}

			/* Possible problem if the integrand stalls for a long time but resume growing afterwards! ATT */
			if (NV_Ith_S (os->y, 0) == last_y)
				break;
			last_y = NV_Ith_S (os->y, 0);
			
			if (flag == CV_ROOT_RETURN)
				break;
		}

		if (ncm_cmp (last_y, os->yf, 1.0e-2, 0.0) != 0)
			g_warning ("ncm_ode_spline_prepare: system has saturated at `% 22.15g' before attaining the required final value `% 22.15g'.",
			           last_y, os->yf);
	}
	
  ncm_spline_set_array (os->s, os->x_array, os->y_array, TRUE);
  os->s_init = TRUE;
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
  g_object_unref (os);
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
  g_clear_object (os);
}

/**
 * ncm_ode_spline_set_interval:
 * @os: a #NcmOdeSpline
 * @yi: FIXME
 * @xi: FIXME
 * @xf: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_interval (NcmOdeSpline *os, gdouble yi, gdouble xi, gdouble xf)
{
  g_assert_cmpfloat (xf, >, xi);
  
  ncm_ode_spline_set_xi (os, xi);
  ncm_ode_spline_set_xf (os, xf);
  ncm_ode_spline_set_yi (os, yi);
}

/**
 * ncm_ode_spline_set_reltol:
 * @os: a #NcmOdeSpline
 * @reltol: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_reltol (NcmOdeSpline *os, gdouble reltol)
{
  os->reltol = reltol;
}

/**
 * ncm_ode_spline_set_abstol:
 * @os: a #NcmOdeSpline
 * @abstol: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_abstol (NcmOdeSpline *os, gdouble abstol)
{
  os->abstol = abstol;
}

/**
 * ncm_ode_spline_set_xi:
 * @os: a #NcmOdeSpline
 * @xi: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_xi (NcmOdeSpline *os, gdouble xi)
{
  os->xi = xi;
}

/**
 * ncm_ode_spline_set_xf:
 * @os: a #NcmOdeSpline
 * @xf: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_xf (NcmOdeSpline *os, gdouble xf)
{
  os->xf = xf;
	if (gsl_finite (os->yf))
	{
		g_warning ("ncm_ode_spline_set_xf: setting xf when yf was also set, yf will take precedence.");
	}
}

/**
 * ncm_ode_spline_set_yi:
 * @os: a #NcmOdeSpline
 * @yi: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_yi (NcmOdeSpline *os, gdouble yi)
{
  os->yi = yi;
}

/**
 * ncm_ode_spline_set_yf:
 * @os: a #NcmOdeSpline
 * @yf: FIXME
 *
 * FIXME
 * 
 */
void 
ncm_ode_spline_set_yf (NcmOdeSpline *os, gdouble yf)
{
  os->yf = yf;
	if (gsl_finite (os->xf))
	{
		g_warning ("ncm_ode_spline_set_yf: setting yf when xf was also set, yf will take precedence.");
	}
}

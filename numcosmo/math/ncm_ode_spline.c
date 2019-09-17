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
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h> 
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h> 
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmOdeSplinePrivate
{
  gpointer cvode;
  SUNNonlinearSolver NLS;
  N_Vector y;
  GArray *y_array;
  GArray *x_array;
  gdouble xi;
  gdouble xf;
  gdouble yi;
	gdouble yf;
  gdouble reltol;
  gdouble abstol;
  NcmOdeSplineDydx dydx;
  gboolean s_init;
  gboolean cvode_init;
  gboolean hnil;
  gboolean stop_hnil;
  gboolean auto_abstol;
  gdouble ini_step;
  NcmModelCtrl *ctrl;
};

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
  PROP_AUTO_ABSTOL,
  PROP_INI_STEP,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmOdeSpline, ncm_ode_spline, G_TYPE_OBJECT);

typedef struct _NcmOdeSplineDydxData
{
	NcmOdeSpline *os;
  gpointer userdata;
} NcmOdeSplineDydxData;

static void
ncm_ode_spline_init (NcmOdeSpline *os)
{
  NcmOdeSplinePrivate * const self = os->priv = ncm_ode_spline_get_instance_private (os);

  os->spline        = NULL;
  self->cvode       = CVodeCreate (CV_ADAMS);
  self->cvode_init  = FALSE;
  self->y           = N_VNew_Serial (1);
  self->NLS         = NULL;
  self->y_array     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  self->x_array     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  self->xi          = GSL_NAN;
  self->xf          = GSL_NAN;
  self->yi          = GSL_NAN;
	self->yf          = GSL_NAN;
  self->reltol      = 0.0;
  self->abstol      = 0.0;
  self->dydx        = NULL;
  self->s_init      = FALSE;
  self->hnil        = FALSE;
  self->stop_hnil   = FALSE;
  self->auto_abstol = FALSE;
  self->ini_step    = 0.0;
  self->ctrl        = ncm_model_ctrl_new (NULL);

  self->NLS = SUNNonlinSol_FixedPoint (self->y, 0);
  NCM_CVODE_CHECK (self->NLS, "SUNNonlinSol_FixedPoint", 0, );
}

static void
_ncm_ode_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmOdeSpline *os = NCM_ODE_SPLINE (object);
  NcmOdeSplinePrivate * const self = os->priv;
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
      self->dydx = (NcmOdeSplineDydx) g_value_get_pointer (value);
      break;
    case PROP_SPLINE:
      ncm_spline_clear (&os->spline);
      os->spline = g_value_dup_object (value);
      break;
    case PROP_STOP_HNIL:
      self->stop_hnil = g_value_get_boolean (value);
      break;
    case PROP_AUTO_ABSTOL:
      ncm_ode_spline_auto_abstol (os, g_value_get_boolean (value));
      break;
    case PROP_INI_STEP:
      ncm_ode_spline_set_ini_step (os, g_value_get_double (value));
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
  NcmOdeSplinePrivate * const self = os->priv;
  g_return_if_fail (NCM_IS_ODE_SPLINE (object));

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, self->abstol);
      break;
    case PROP_XI:
      g_value_set_double (value, self->xi);
      break;
    case PROP_XF:
      g_value_set_double (value, self->xf);
      break;
    case PROP_YI:
      g_value_set_double (value, self->yi);
      break;
    case PROP_YF:
      g_value_set_double (value, self->yf);
      break;
    case PROP_DYDX:
      g_value_set_pointer (value, self->dydx);
      break;
    case PROP_SPLINE:
      g_value_set_object (value, os->spline);
      break;
    case PROP_STOP_HNIL:
      g_value_set_boolean (value, self->stop_hnil);
      break;
    case PROP_AUTO_ABSTOL:
      g_value_set_boolean (value, self->auto_abstol);
      break;
    case PROP_INI_STEP:
      g_value_set_boolean (value, ncm_ode_spline_get_ini_step (os));
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
  NcmOdeSplinePrivate * const self = os->priv;

  ncm_spline_clear (&os->spline);
  g_clear_pointer (&self->x_array, g_array_unref);
  g_clear_pointer (&self->y_array, g_array_unref);  
  ncm_model_ctrl_clear (&self->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_ode_spline_parent_class)->dispose (object);
}

static void
_ncm_ode_spline_finalize (GObject *object)
{
  NcmOdeSpline *os = NCM_ODE_SPLINE (object);
  NcmOdeSplinePrivate * const self = os->priv;

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode = NULL;
  }
  if (self->y != NULL)
  {
    N_VDestroy (self->y);
    self->y = NULL;
  }
  if (self->NLS != NULL)
  {
    SUNNonlinSolFree (self->NLS);
    self->NLS = NULL;
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
  g_object_class_install_property (object_class,
                                   PROP_AUTO_ABSTOL,
                                   g_param_spec_boolean ("auto-abstol",
                                                         NULL,
                                                         "Automatic abstol",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INI_STEP,
                                   g_param_spec_double ("ini-step",
                                                        NULL,
                                                        "Integration initial step size",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static gint
_ncm_ode_spline_f (realtype x, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmOdeSplineDydxData *dydx_data = (NcmOdeSplineDydxData *) f_data;
  NcmOdeSplinePrivate * const self = dydx_data->os->priv;
  NV_Ith_S (ydot, 0) = self->dydx (NV_Ith_S (y, 0), x, dydx_data->userdata);

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
  NcmOdeSplinePrivate * const self = dydx_data->os->priv;

	gout[0] = (NV_Ith_S (y, 0) - self->yf);

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
  NcmOdeSplinePrivate * const self = os->priv;
  NcmOdeSplineDydxData f_data = {os, userdata};
  gdouble x, x0;
  gint flag;

  NV_Ith_S (self->y, 0) = self->yi;

  if (self->auto_abstol)
  {
    self->abstol = fabs (self->dydx (NV_Ith_S (self->y, 0), self->xi, userdata) * self->reltol * NCM_ODE_SPLINE_MIN_STEP);
  }
  else if ((self->yi == 0.0) && (self->abstol == 0.0))
    g_error ("ncm_ode_spline_prepare: cannot integrate system where y_ini == 0.0 and abstol == 0.0.");
  
  if (!self->cvode_init)
  {
    flag = CVodeInit (self->cvode, &_ncm_ode_spline_f, self->xi, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetNonlinearSolver (self->cvode, self->NLS);
    NCM_CVODE_CHECK (&flag, "CVodeSetNonlinearSolver", 1, );
    
    self->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (self->cvode, self->xi, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
        
    flag = CVodeSetNonlinearSolver (self->cvode, self->NLS);
    NCM_CVODE_CHECK (&flag, "CVodeSetNonlinearSolver", 1, );
  }

  g_array_set_size (self->x_array, 0);
  g_array_set_size (self->y_array, 0);

  g_array_append_val (self->x_array, self->xi);
  g_array_append_val (self->y_array, NV_Ith_S (self->y, 0));

  flag = CVodeSStolerances (self->cvode, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );
  
  flag = CVodeSetMaxNumSteps (self->cvode, NCM_INTEGRAL_PARTITION);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
  
  flag = CVodeSetUserData (self->cvode, &f_data);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  if (self->ini_step > 0.0)
  {
    flag = CVodeSetInitStep (self->cvode, self->ini_step);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );  
  }

  x0 = self->xi;
  if (!gsl_finite (self->dydx (NV_Ith_S (self->y, 0), x0, f_data.userdata)))
    g_error ("ncm_ode_spline_prepare: not finite integrand at (% 22.15g, % 22.15g; % 22.15g).", 
             x0, 
             NV_Ith_S (self->y, 0), 
             self->dydx (NV_Ith_S (self->y, 0), x0, f_data.userdata));

	if (!gsl_finite (self->yf))
	{
		g_assert (gsl_finite (self->xf));
		
		flag = CVodeSetStopTime (self->cvode, self->xf);
		NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

		self->hnil = FALSE;
		while (TRUE)
		{
			flag = CVode (self->cvode, self->xf, self->y, &x, CV_ONE_STEP);
			NCM_CVODE_CHECK (&flag, "ncm_ode_spline_prepare[CVode]", 1, );

			if (G_UNLIKELY (self->hnil))
			{
				if (self->stop_hnil)
				{
					g_error ("ncm_ode_spline_prepare: cannot integrate function %d.", flag);
				}
				else
					break;
			}

			if (x > x0 + fabs (x0) * NCM_ODE_SPLINE_MIN_STEP)
			{
				g_array_append_val (self->x_array, x);
				g_array_append_val (self->y_array, NV_Ith_S (self->y, 0));
				x0 = x;
			}

			if (x == self->xf)
				break;
		}
	}
	else
	{
		const gdouble xf = (self->xi != 0.0) ? self->xi * 2.0 : 1.0;
		gdouble last_y = GSL_NEGINF;
		
		flag = CVodeRootInit (self->cvode, 1, &_ncm_ode_spline_yf_root);
		NCM_CVODE_CHECK (&flag, "CVodeRootInit", 1, );    
		
		self->hnil = FALSE;
		while (TRUE)
		{
			flag = CVode (self->cvode, xf, self->y, &x, CV_ONE_STEP);
			NCM_CVODE_CHECK (&flag, "ncm_ode_spline_prepare[CVode]", 1, );

			if (G_UNLIKELY (self->hnil))
			{
				if (self->stop_hnil)
				{
					g_error ("ncm_ode_spline_prepare: cannot integrate function %d.", flag);
				}
				else
					break;
			}

			if (x > x0 + fabs (x0) * NCM_ODE_SPLINE_MIN_STEP)
			{
				g_array_append_val (self->x_array, x);
				g_array_append_val (self->y_array, NV_Ith_S (self->y, 0));
				x0 = x;
			}

			/* Possible problem if the integrand stalls for a long time but resume growing afterwards! ATT */
			if (NV_Ith_S (self->y, 0) == last_y)
				break;
			last_y = NV_Ith_S (self->y, 0);
			
			if (flag == CV_ROOT_RETURN)
				break;
		}

		if (ncm_cmp (last_y, self->yf, 1.0e-2, 0.0) != 0)
			g_warning ("ncm_ode_spline_prepare: system has saturated at `% 22.15g' before attaining the required final value `% 22.15g'.",
			           last_y, self->yf);
	}
	
  ncm_spline_set_array (os->spline, self->x_array, self->y_array, TRUE);
  self->s_init = TRUE;
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
  NcmOdeSplinePrivate * const self = os->priv;
  self->reltol = reltol;
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
  NcmOdeSplinePrivate * const self = os->priv;
  self->abstol = abstol;
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
  NcmOdeSplinePrivate * const self = os->priv;
  self->xi = xi;
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
  NcmOdeSplinePrivate * const self = os->priv;
  self->xf = xf;
	if (gsl_finite (self->yf))
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
  NcmOdeSplinePrivate * const self = os->priv;
  self->yi = yi;
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
  NcmOdeSplinePrivate * const self = os->priv;
  self->yf = yf;
	if (gsl_finite (self->xf))
	{
		g_warning ("ncm_ode_spline_set_yf: setting yf when xf was also set, yf will take precedence.");
	}
}

/**
 * ncm_ode_spline_auto_abstol:
 * @os: a #NcmOdeSpline
 * @on: Whether to turn on the auto-abstol
 *
 * If @on is TRUE, the object uses the value of $\mathrm{d}y_i$ to estimate the 
 * abstol as $T_\mathrm{abs} = \dot{y}_i \mathrm{d}t_m T_\mathrm{rel}$,
 * where $T_\mathrm{rel}$ is the relative tolerance and $\mathrm{d}t_m$ is the 
 * minimum time step #NCM_ODE_SPLINE_MIN_STEP. Useful when computing integrals as ODEs.
 */
void 
ncm_ode_spline_auto_abstol (NcmOdeSpline *os, gboolean on)
{
  NcmOdeSplinePrivate * const self = os->priv;
  self->auto_abstol = on;
}

/**
 * ncm_ode_spline_set_ini_step:
 * @os: a #NcmOdeSpline
 * @ini_step: the initial step
 * 
 * Sets a guess for the initial step size. If @ini_step is
 * zero it uses the automatic determination based on the
 * tolerances.
 * 
 */
void 
ncm_ode_spline_set_ini_step (NcmOdeSpline *os, gdouble ini_step)
{
  NcmOdeSplinePrivate * const self = os->priv;
  self->ini_step = ini_step;
}

/**
 * ncm_ode_spline_get_ini_step:
 * @os: a #NcmOdeSpline
 * 
 * Gets the current guess for the initial step size. 
 * 
 * Returns: the current value of the initial guess (zero means disabled).
 */
gdouble
ncm_ode_spline_get_ini_step (NcmOdeSpline *os)
{
  NcmOdeSplinePrivate * const self = os->priv;
  return self->ini_step;
}

/**
 * ncm_ode_spline_peek_spline:
 * @os: a #NcmOdeSpline
 *
 * Peeks at the last prepared spline.
 * 
 * Returns: (transfer none): the last prepared spline.
 */

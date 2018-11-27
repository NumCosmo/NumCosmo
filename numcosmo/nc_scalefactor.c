/***************************************************************************
 *            nc_scalefactor.c
 *
 *  Wed Nov 12 14:46:09 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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
 * SECTION:nc_scalefactor
 * @title: NcScalefactor
 * @short_description: Scale factor as a function of the conformal time.
 *
 * Integrates the first order Friedmann equation, 
 * $$E^2 = \frac{\rho}{\rho_{\mathrm{crit}0}} + \Omega_{k0} x^2.$$ Where 
 * ${\mathrm{crit}0}$ is the critical density today [nc_hicosmo_crit_density()],  
 * $E = H / H_0$ is the dimensionless Hubble function [nc_hicosmo_E()] 
 * and $\Omega_{k0}$ is the curvature parameter today [nc_hicosmo_Omega_k0()].
 * 
 * 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_scalefactor.h"
#include "nc_enum_types.h"
#include "math/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <nvector/nvector_serial.h>
#if HAVE_SUNDIALS_MAJOR == 2
#define SUN_DENSE_ACCESS DENSE_ELEM
#elif HAVE_SUNDIALS_MAJOR == 3
#define SUN_DENSE_ACCESS SM_ELEMENT_D
#endif 

#endif /* NUMCOSMO_GIR_SCAN */


enum
{
  PROP_0,
  PROP_ZI,
  PROP_A0,
  PROP_A0_CONF_NORM,
  PROP_DIST,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_TTYPE,
};

G_DEFINE_TYPE (NcScalefactor, nc_scalefactor, G_TYPE_OBJECT);

static void
nc_scalefactor_init (NcScalefactor *a)
{
  a->ttype       = 0;
	a->a_eta       = ncm_spline_cubic_notaknot_new ();
	a->eta_a       = ncm_spline_cubic_notaknot_new ();
	a->t_eta       = ncm_spline_cubic_notaknot_new ();
	a->eta_t       = ncm_spline_cubic_notaknot_new ();
  a->spline_init = FALSE;
	a->dist        = NULL;
	a->ctrl        = ncm_model_ctrl_new (NULL);
	a->zf          = 0.0;
	a->a0          = 0.0;

  a->sets_conf_norm = FALSE;

  a->cvode       = CVodeCreate (CV_BDF, CV_NEWTON);
  NCM_CVODE_CHECK ((void *)a->cvode, "CVodeCreate", 0, );
  a->cvode_init  = FALSE;
  a->quad_init   = FALSE;

  a->reltol      = 0.0;
  a->abstol      = 0.0;

  a->y           = N_VNew_Serial (1);
  a->yQ          = N_VNew_Serial (1);
#if HAVE_SUNDIALS_MAJOR == 3
  a->A           = SUNDenseMatrix (1, 1);
  a->LS          = SUNDenseLinearSolver (a->y, a->A);
  
  NCM_CVODE_CHECK ((gpointer)a->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer)a->LS, "SUNDenseLinearSolver", 0, );
#endif

}

static void
nc_scalefactor_dispose (GObject *object)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  
  ncm_spline_clear (&a->a_eta);
  ncm_spline_clear (&a->eta_a);
  ncm_spline_clear (&a->t_eta);
  ncm_spline_clear (&a->eta_t);

	ncm_model_ctrl_clear (&a->ctrl);
	nc_distance_clear (&a->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_scalefactor_parent_class)->dispose (object);
}

static void
nc_scalefactor_finalize (GObject *object)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  
  CVodeFree (&a->cvode);
  N_VDestroy (a->y);
  N_VDestroy (a->yQ);

#if HAVE_SUNDIALS_MAJOR == 3
  if (a->A != NULL)
  {
    SUNMatDestroy (a->A);
    a->A = NULL;
  }
  if (a->LS != NULL)
  {
    SUNLinSolFree (a->LS);
    a->LS = NULL;
  }
#endif

  /* Chain up : end */
  G_OBJECT_CLASS (nc_scalefactor_parent_class)->finalize (object);
}

static void
nc_scalefactor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  g_return_if_fail (NC_IS_SCALEFACTOR (object));

  switch (prop_id)
  {
    case PROP_ZI:
      nc_scalefactor_set_zf (a, g_value_get_double (value));
      break;
    case PROP_A0:
      nc_scalefactor_set_a0 (a, g_value_get_double (value));
      break;
    case PROP_A0_CONF_NORM:
      nc_scalefactor_set_a0_conformal_normal (a, g_value_get_boolean (value));
      break;
    case PROP_DIST:
      a->dist = g_value_dup_object (value);
      break;
    case PROP_RELTOL:
      nc_scalefactor_set_reltol (a, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      nc_scalefactor_set_abstol (a, g_value_get_double (value));
      break;
    case PROP_TTYPE:
      nc_scalefactor_set_time_type (a, g_value_get_flags (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_scalefactor_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  g_return_if_fail (NC_IS_SCALEFACTOR (object));

  switch (prop_id)
  {
    case PROP_ZI:
      g_value_set_double (value, nc_scalefactor_get_zf (a));
      break;
    case PROP_A0:
      g_value_set_double (value, nc_scalefactor_get_a0 (a));
      break;
    case PROP_A0_CONF_NORM:
      g_value_set_boolean (value, a->sets_conf_norm);
      break;
    case PROP_DIST:
      g_value_set_object (value, a->dist);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_scalefactor_get_reltol (a));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, nc_scalefactor_get_abstol (a));
      break;
    case PROP_TTYPE:
      g_value_set_flags (value, nc_scalefactor_get_time_type (a));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_scalefactor_class_init (NcScalefactorClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = nc_scalefactor_set_property;
  object_class->get_property = nc_scalefactor_get_property;
  object_class->dispose      = nc_scalefactor_dispose;
  object_class->finalize     = nc_scalefactor_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Initial redshift",
                                                        0.0, G_MAXDOUBLE, NC_SCALEFACTOR_DEFAULT_ZF,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_A0,
                                   g_param_spec_double ("a0",
                                                        NULL,
                                                        "Scale factor today a_0",
                                                        0.0, G_MAXDOUBLE, NC_SCALEFACTOR_DEFAULT_A0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_A0_CONF_NORM,
                                   g_param_spec_boolean ("a0-conformal-normal",
                                                         NULL,
                                                         "Scale factor today a_0 from normalized curvature radius",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NC_SCALEFACTOR_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NC_SCALEFACTOR_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TTYPE,
                                   g_param_spec_flags ("time-type",
                                                        NULL,
                                                        "Different time type to be integrated",
                                                        NC_TYPE_SCALEFACTOR_TIME_TYPE, 0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static gint dz_deta_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);


#if HAVE_SUNDIALS_MAJOR == 2
static gint dz_deta_J (_NCM_SUNDIALS_INT_TYPE N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#elif HAVE_SUNDIALS_MAJOR == 3
static gint dz_deta_J (realtype lambda, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

static gint dt_deta (realtype eta, N_Vector y, N_Vector yQdot, gpointer user_data);

/**
 * nc_scalefactor_new:
 * @ttype: a #NcScalefactorTimeType
 * @zf: maximum redshift $z_f$
 * @dist: a #NcDistance
 *
 * Creates a new #NcScalefactor valid for calculations in
 * the $[0, z_f]$ interval.
 *
 * Returns: (transfer full): a new #NcScalefactor.
 */
NcScalefactor *
nc_scalefactor_new (NcScalefactorTimeType ttype, gdouble zf, NcDistance *dist)
{
  NcScalefactor *a;
  if (dist == NULL)
    dist = nc_distance_new (1.0);
  a = g_object_new (NC_TYPE_SCALEFACTOR,
                                     "time-type", ttype,
                                     "zf", zf,
                                     "dist", dist,
                                     NULL);
  nc_distance_free (dist);
  return a;
}

/**
 * nc_scalefactor_ref:
 * @a: a #NcScalefactor
 *
 * Increases the reference count of @a by one.
 *
 * Returns: (transfer full): @a.
 */
NcScalefactor *
nc_scalefactor_ref (NcScalefactor *a)
{
  return g_object_ref (a);
}

/**
 * nc_scalefactor_free:
 * @a: a #NcScalefactor
 *
 * Decreases the reference count of @a by one.
 * 
 */
void
nc_scalefactor_free (NcScalefactor *a)
{
  g_object_unref (a);
}

/**
 * nc_scalefactor_clear:
 * @a: a #NcScalefactor
 *
 * If *@a is different from NULL, decreases the reference 
 * count of *@a by one and sets *@a to NULL.
 * 
 */
void
nc_scalefactor_clear (NcScalefactor **a)
{
  g_clear_object (a);
}

static void _nc_scalefactor_init (NcScalefactor *a, NcHICosmo *cosmo);
static void _nc_scalefactor_calc_spline (NcScalefactor *a, NcHICosmo *cosmo);

/**
 * nc_scalefactor_prepare:
 * @a: a #NcScalefactor
 * @cosmo: FIXME
 *
 * FIXME
 *
 */
void
nc_scalefactor_prepare (NcScalefactor *a, NcHICosmo *cosmo)
{
  nc_distance_prepare_if_needed (a->dist, cosmo);
  
	_nc_scalefactor_init (a, cosmo);
	_nc_scalefactor_calc_spline (a, cosmo);
	ncm_model_ctrl_update (a->ctrl, NCM_MODEL (cosmo));
}

/**
 * nc_scalefactor_prepare_if_needed:
 * @a: a #NcScalefactor
 * @cosmo: FIXME
 *
 * FIXME
 *
 */
void
nc_scalefactor_prepare_if_needed (NcScalefactor *a, NcHICosmo *cosmo)
{
	if (ncm_model_ctrl_update (a->ctrl, NCM_MODEL (cosmo)))
		nc_scalefactor_prepare (a, cosmo);
}

/**
 * nc_scalefactor_set_zf:
 * @a: a #NcScalefactor
 * @zf: final redshift $z_f$
 * 
 * Sets the final redshift of the integration.
 * 
 */
void 
nc_scalefactor_set_zf (NcScalefactor *a, const gdouble zf)
{
  if (a->zf != zf)
  {
    a->zf = zf;
    ncm_model_ctrl_force_update (a->ctrl);
  }
}

/**
 * nc_scalefactor_set_a0:
 * @a: a #NcScalefactor
 * @a0: scale factor today $a_0$
 * 
 * Sets the value of the scale factor today.
 * 
 */
void 
nc_scalefactor_set_a0 (NcScalefactor *a, const gdouble a0)
{
  if (a->a0 != a0)
  {
    a->a0 = a0;
    ncm_model_ctrl_force_update (a->ctrl);
  }
}

/**
 * nc_scalefactor_set_reltol:
 * @a: a #NcScalefactor
 * @reltol: relative tolerance
 * 
 * Sets the relative tolerance of the integration.
 * 
 */
void 
nc_scalefactor_set_reltol (NcScalefactor *a, const gdouble reltol)
{
  if (a->reltol != reltol)
  {
    a->reltol = reltol;
    ncm_model_ctrl_force_update (a->ctrl);
  }
}

/**
 * nc_scalefactor_set_abstol:
 * @a: a #NcScalefactor
 * @abstol: absolute tolerance
 * 
 * Sets the absolute tolerance of the integration.
 * 
 */
void
nc_scalefactor_set_abstol (NcScalefactor *a, const gdouble abstol)
{
  if (a->abstol != abstol)
  {
    a->abstol = abstol;
    ncm_model_ctrl_force_update (a->ctrl);
  }
}

/**
 * nc_scalefactor_set_time_type:
 * @a: a #NcScalefactor
 * @ttype: a #NcScalefactorTimeType flag
 * 
 * Sets the which other time variables it should integrate.
 * 
 */
void 
nc_scalefactor_set_time_type (NcScalefactor *a, NcScalefactorTimeType ttype)
{
  if (a->ttype != ttype)
  {
    a->ttype = ttype;
    ncm_model_ctrl_force_update (a->ctrl);
  }
}

/**
 * nc_scalefactor_set_a0_conformal_normal:
 * @a: a #NcScalefactor
 * @enable: a boolean
 * 
 * When @enable is TRUE, it sets the value of the scale factor 
 * today $a_0$, assuming that the conformal hypersurface 
 *  the spatial hypersurface where ($a=1$) hascurvature
 * equals to 1Mpc, i.e., $1/\sqrt{K} = 1\,\mathrm{Mpc}$.
 * If @enable is FALSE it lets $a_0$ untouched.   * 
 * 
 */
void 
nc_scalefactor_set_a0_conformal_normal (NcScalefactor *a, gboolean enable) 
{
  if ((!a->sets_conf_norm && enable) || (a->sets_conf_norm && !enable))
  {
    a->sets_conf_norm = enable;
    ncm_model_ctrl_force_update (a->ctrl);
  }
}

/**
 * nc_scalefactor_get_zf:
 * @a: a #NcScalefactor
 * 
 * Gets the current final redshift $z_f$.
 * 
 * Returns: $z_f$.
 */
gdouble 
nc_scalefactor_get_zf (NcScalefactor *a)
{
  return a->zf;
}

/**
 * nc_scalefactor_get_a0:
 * @a: a #NcScalefactor
 * 
 * Gets the current value of the scale factor today $a_0$.
 * 
 * Returns: $a_0$.
 */
gdouble 
nc_scalefactor_get_a0 (NcScalefactor *a)
{
  return a->a0;
}

/**
 * nc_scalefactor_get_reltol:
 * @a: a #NcScalefactor
 * 
 * Gets the current relative tolerance.
 * 
 * Returns: reltol.
 */
gdouble 
nc_scalefactor_get_reltol (NcScalefactor *a)
{
  return a->reltol;
}

/**
 * nc_scalefactor_get_abstol:
 * @a: a #NcScalefactor
 * 
 * Gets the current absolute tolerance.
 * 
 * Returns: abstol.
 */
gdouble 
nc_scalefactor_get_abstol (NcScalefactor *a)
{
  return a->abstol;
}

/**
 * nc_scalefactor_get_time_type:
 * @a: a #NcScalefactor
 * 
 * Gets the current time type flag.
 * 
 * Returns: ttime.
 */
NcScalefactorTimeType 
nc_scalefactor_get_time_type (NcScalefactor *a)
{
  return a->ttype;
}

static void nc_scalefactor_init_cvode (NcScalefactor *a, NcHICosmo *cosmo);

static void
_nc_scalefactor_init (NcScalefactor *a, NcHICosmo *cosmo)
{
  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (cosmo);

  NV_Ith_S (a->y, 0) = a->zf;

  a->eta_i = nc_distance_conformal_time (a->dist, cosmo, a->zf);
  a->eta_f = nc_distance_conformal_time (a->dist, cosmo, 0.0);

  if (a->sets_conf_norm)
  {
    if (fabs (Omega_k0) < NC_SCALEFACTOR_OMEGA_K_ZERO)
    {
      a->a0 = 1.0;
    }
    else
    {
      const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
      a->a0 = RH_Mpc / sqrt (fabs (Omega_k0));
    }
  }
  nc_scalefactor_init_cvode (a, cosmo);

  return;
}

static void
nc_scalefactor_init_cvode (NcScalefactor *a, NcHICosmo *cosmo)
{
  gint flag;

  if (a->cvode_init)
  {
    flag = CVodeReInit (a->cvode, a->eta_i, a->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }
  else
  {
    flag = CVodeInit (a->cvode, &dz_deta_f, a->eta_i, a->y);
    NCM_CVODE_CHECK (&flag, "CVodeMalloc", 1, );
    a->cvode_init = TRUE;
  }

  flag = CVodeSStolerances (a->cvode, a->reltol, a->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

  flag = CVodeSetStopTime(a->cvode, a->eta_f);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetUserData (a->cvode, cosmo);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetMaxNumSteps (a->cvode, 50000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );


#if HAVE_SUNDIALS_MAJOR == 2
  flag = CVDense(a->cvode, 1);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (a->cvode, &dz_deta_J);
  NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );
#elif HAVE_SUNDIALS_MAJOR == 3
  flag = CVDlsSetLinearSolver (a->cvode, a->LS, a->A);
  NCM_CVODE_CHECK (&flag, "CVDlsSetLinearSolver", 1, );

  flag = CVDlsSetJacFn (a->cvode, &dz_deta_J);
  NCM_CVODE_CHECK (&flag, "CVDlsSetJacFn", 1, );
#endif
  
  return;
}

/**
 * nc_scalefactor_eval_z_eta:
 * @a: a #NcScalefactor
 * @eta: conformal time $\eta$
 *
 * Calculates the value of the redshift in $\eta$,
 * i.e., $z(\eta)$.
 *
 * Returns: $z(\eta)$.
 */
gdouble
nc_scalefactor_eval_z_eta (NcScalefactor *a, const gdouble eta)
{
  return -ncm_spline_eval (a->a_eta, eta);
}

/**
 * nc_scalefactor_eval_eta_z:
 * @a: a #NcScalefactor
 * @z: redshift $z$
 *
 * Calculates the value of the conformal time at $z$,
 * i.e., $\eta(z)$. 
 *
 * Returns: $\eta(z)$.
 */
gdouble
nc_scalefactor_eval_eta_z (NcScalefactor *a, const gdouble z)
{
  return ncm_spline_eval (a->eta_a, - z);
}

/**
 * nc_scalefactor_eval_eta_x:
 * @a: a #NcScalefactor
 * @x: redshift x variable $x = 1 + z$
 *
 * Calculates the value of the conformal time at $x$,
 * i.e., $\eta(z(x))$.
 *
 * Returns: $\eta(z(x))$.
 */
gdouble
nc_scalefactor_eval_eta_x (NcScalefactor *a, const gdouble x)
{
  return ncm_spline_eval (a->eta_a, -(x - 1.0));
}

/**
 * nc_scalefactor_eval_a_eta:
 * @a: a #NcScalefactor
 * @eta: conformal time $\eta$
 *
 * Calculates the value of the scale factor in $\eta$,
 * i.e., $a(\eta)$.
 *
 * Returns: $a(\eta)$.
 */
gdouble
nc_scalefactor_eval_a_eta (NcScalefactor *a, const gdouble eta)
{
  gdouble mz;
  mz = ncm_spline_eval (a->a_eta, eta);
  return a->a0 / (1.0 - mz);
}

/**
 * nc_scalefactor_eval_t_eta:
 * @a: a #NcScalefactor
 * @eta: conformal time $\eta$
 *
 * Calculates the value of the cosmic time at $\eta$,
 * i.e., $t(\eta)$.
 *
 * Returns: $t(\eta)$.
 */
gdouble 
nc_scalefactor_eval_t_eta (NcScalefactor *a, const gdouble eta)
{
  return ncm_spline_eval (a->t_eta, eta);
}

/**
 * nc_scalefactor_eval_eta_t:
 * @a: a #NcScalefactor
 * @t: cosmic time $t$
 *
 * Calculates the value of the conformal time at $t$,
 * i.e., $\eta(t)$.
 *
 * Returns: $\eta(t)$.
 */
gdouble 
nc_scalefactor_eval_eta_t (NcScalefactor *a, const gdouble t)
{
  return ncm_spline_eval (a->eta_t, t);
}

static void
_nc_scalefactor_calc_spline (NcScalefactor *a, NcHICosmo *cosmo)
{
  GArray *eta_a, *mz_a, *t_a = NULL;
  gdouble eta, mzi, last_eta;
  gint flag;
  gboolean int_t = (a->ttype & NC_SCALEFACTOR_TIME_TYPE_COSMIC);

	if (!ncm_spline_is_empty (a->a_eta))
	{
		eta_a = ncm_vector_get_array (a->a_eta->xv);
		mz_a  = ncm_vector_get_array (a->a_eta->yv);
		g_array_set_size (eta_a, 0);
		g_array_set_size (mz_a, 0);
	}
	else
	{
		eta_a = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
		mz_a  = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
	}

  last_eta = a->eta_i;
  mzi      = -a->zf;
  g_array_append_val (eta_a, a->eta_i);
  g_array_append_val (mz_a, mzi);

  if (int_t)
  {
    const gdouble ti = nc_distance_cosmic_time (a->dist, cosmo, a->zf);
    NV_Ith_S (a->yQ, 0) = ti;

    if (a->quad_init)
    {
      flag = CVodeQuadReInit (a->cvode, a->yQ);
      NCM_CVODE_CHECK (&flag, "_nc_scalefactor_calc_spline[CVodeQuadReInit]", 1, );
    }
    else
    {
      flag = CVodeQuadInit (a->cvode, &dt_deta, a->yQ);
      NCM_CVODE_CHECK (&flag, "_nc_scalefactor_calc_spline[CVodeQuadInit]", 1, );
      a->quad_init = TRUE;
    }

    if (!ncm_spline_is_empty (a->t_eta))
    {
      t_a = ncm_vector_get_array (a->t_eta->yv);
      g_array_set_size (t_a, 0);
    }
    else
      t_a = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);

    g_array_append_val (t_a, ti);
  }
  else if (a->quad_init)
  {
    CVodeQuadFree (a->cvode);
    a->quad_init = FALSE;
  }
  
  if (int_t)
  {
    while (TRUE)
    {
      gdouble eta_quad, t;
      
      flag = CVode (a->cvode, a->eta_f, a->y, &eta, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "_nc_scalefactor_calc_spline[CVode]", 1, );
      flag = CVodeGetQuad (a->cvode, &eta_quad, a->yQ);

      if ((eta - last_eta) / last_eta > NC_SCALEFACTOR_MIN_ETA_STEP)
      {
        mzi = - NV_Ith_S (a->y, 0);
        t   = NV_Ith_S (a->yQ, 0);

        g_array_append_val (eta_a, eta);
        g_array_append_val (mz_a, mzi);
        g_array_append_val (t_a, t);
        last_eta = eta;
      }
      
      if (a->eta_f == eta)
        break;
    }
  }
  else
  {
    while (TRUE)
    {
      flag = CVode (a->cvode, a->eta_f, a->y, &eta, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "_nc_scalefactor_calc_spline[CVode]", 1, );

      if ((eta - last_eta) / last_eta > NC_SCALEFACTOR_MIN_ETA_STEP)
      {
        mzi = - NV_Ith_S (a->y, 0);

        g_array_append_val (eta_a, eta);
        g_array_append_val (mz_a, mzi);

        last_eta = eta;
      }
      
      if (a->eta_f == eta)
        break;
    }
  }

	if (fabs (g_array_index (mz_a, gdouble, mz_a->len - 1)) < 1e-10)
	{
		g_array_index (mz_a, gdouble, mz_a->len - 1) = 0.0;
	}
	else
		g_error ("_nc_scalefactor_calc_spline today redshift must be zero not % 20.15g.",
		         fabs (g_array_index (mz_a, gdouble, mz_a->len - 1)));

	ncm_spline_set_array (a->a_eta, eta_a, mz_a, TRUE);
	ncm_spline_set_array (a->eta_a, mz_a, eta_a, TRUE);
  g_array_unref (eta_a);
  g_array_unref (mz_a);

  if (int_t)
  {
    ncm_spline_set_array (a->t_eta, eta_a, t_a, TRUE);
    ncm_spline_set_array (a->eta_t, t_a, eta_a, TRUE);
    g_array_unref (t_a);
  }
  
  return;
}

static gint
dz_deta_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *cosmo = NC_HICOSMO (f_data);
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble E = nc_hicosmo_E (cosmo, z);

  NCM_UNUSED (t);

  NV_Ith_S (ydot, 0) = - E;
  return 0;
}

static gint
#if HAVE_SUNDIALS_MAJOR == 2
dz_deta_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#elif HAVE_SUNDIALS_MAJOR == 3
dz_deta_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#endif
{
  NcHICosmo *cosmo = NC_HICOSMO (jac_data);
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble E = nc_hicosmo_E (cosmo, z);
  const gdouble dE2_dz = nc_hicosmo_dE2_dz (cosmo, z);

  NCM_UNUSED (t);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);
  
  SUN_DENSE_ACCESS (J, 0, 0) = - dE2_dz / (2.0 * E);

  return 0;
}

static gint 
dt_deta (realtype eta, N_Vector y, N_Vector yQdot, gpointer user_data)
{
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble x = 1.0 + z;

  NV_Ith_S (yQdot, 0) = 1.0 / x;

  return 0;
}


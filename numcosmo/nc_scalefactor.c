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
#include <cvode/cvode.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D
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
};

struct _NcScalefactorPrivate
{
  NcDistance *dist;
  NcmModelCtrl *ctrl;
  NcmSpline *a_eta;
  NcmSpline *eta_a;
  NcmSpline *t_eta;
  NcmSpline *eta_t;
  gdouble a0;
  gdouble zf;
  gdouble eta_i;
  gdouble eta_f;
  gdouble t_i;
  gboolean sets_conf_norm;
  gboolean spline_init;
  gboolean cvode_init;
  gboolean quad_init;
  gpointer cvode;
  gdouble reltol;
  gdouble abstol;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcScalefactor, nc_scalefactor, G_TYPE_OBJECT);

static void
nc_scalefactor_init (NcScalefactor *a)
{
  NcScalefactorPrivate * const self = a->priv = nc_scalefactor_get_instance_private (a);

	self->a_eta       = ncm_spline_cubic_notaknot_new ();
	self->eta_a       = ncm_spline_cubic_notaknot_new ();
	self->t_eta       = ncm_spline_cubic_notaknot_new ();
	self->eta_t       = ncm_spline_cubic_notaknot_new ();
  self->spline_init = FALSE;
	self->dist        = NULL;
	self->ctrl        = ncm_model_ctrl_new (NULL);

  self->a0          = 0.0;
	self->zf          = 0.0;
	self->eta_i       = 0.0;
	self->eta_f       = 0.0;
	self->t_i         = 0.0;

  self->sets_conf_norm = FALSE;

  self->cvode       = CVodeCreate (CV_BDF);
  NCM_CVODE_CHECK ((void *)self->cvode, "CVodeCreate", 0, );
  self->cvode_init  = FALSE;
  self->quad_init   = FALSE;

  self->reltol      = 0.0;
  self->abstol      = 0.0;

  self->y           = N_VNew_Serial (2);

  self->A           = SUNDenseMatrix (2, 2);
  self->LS          = SUNDenseLinearSolver (self->y, self->A);
  
  NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
}

static void
nc_scalefactor_dispose (GObject *object)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  NcScalefactorPrivate * const self = a->priv;
  
  ncm_spline_clear (&self->a_eta);
  ncm_spline_clear (&self->eta_a);
  ncm_spline_clear (&self->t_eta);
  ncm_spline_clear (&self->eta_t);

	ncm_model_ctrl_clear (&self->ctrl);
	nc_distance_clear (&self->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_scalefactor_parent_class)->dispose (object);
}

static void
nc_scalefactor_finalize (GObject *object)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  NcScalefactorPrivate * const self = a->priv;
  
  CVodeFree (&self->cvode);
  N_VDestroy (self->y);

  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }
  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (nc_scalefactor_parent_class)->finalize (object);
}

static void
nc_scalefactor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcScalefactor *a = NC_SCALEFACTOR (object);
  NcScalefactorPrivate * const self = a->priv;
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
      self->dist = g_value_dup_object (value);
      break;
    case PROP_RELTOL:
      nc_scalefactor_set_reltol (a, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      nc_scalefactor_set_abstol (a, g_value_get_double (value));
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
  NcScalefactorPrivate * const self = a->priv;
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
      g_value_set_boolean (value, self->sets_conf_norm);
      break;
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_scalefactor_get_reltol (a));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, nc_scalefactor_get_abstol (a));
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
}

static gint dz_deta_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint dz_deta_J (realtype lambda, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/**
 * nc_scalefactor_new:
 * @zf: maximum redshift $z_f$
 * @dist: a #NcDistance
 *
 * Creates a new #NcScalefactor valid for calculations in
 * the $[0, z_f]$ interval.
 *
 * Returns: (transfer full): a new #NcScalefactor.
 */
NcScalefactor *
nc_scalefactor_new (gdouble zf, NcDistance *dist)
{
  NcScalefactor *a;

  if (dist == NULL)
    dist = nc_distance_new (1.0);
  else
    dist = nc_distance_ref (dist);
  
  a = g_object_new (NC_TYPE_SCALEFACTOR,
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
  NcScalefactorPrivate * const self = a->priv;
  nc_distance_prepare_if_needed (self->dist, cosmo);
  
	_nc_scalefactor_init (a, cosmo);
	_nc_scalefactor_calc_spline (a, cosmo);
	ncm_model_ctrl_update (self->ctrl, NCM_MODEL (cosmo));
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
  NcScalefactorPrivate * const self = a->priv;
	if (ncm_model_ctrl_update (self->ctrl, NCM_MODEL (cosmo)))
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
  NcScalefactorPrivate * const self = a->priv;
  if (self->zf != zf)
  {
    self->zf = zf;
    ncm_model_ctrl_force_update (self->ctrl);
    if (self->dist != NULL)
      nc_distance_require_zf (self->dist, zf);
  }
}

/**
 * nc_scalefactor_require_zf:
 * @a: a #NcScalefactor
 * @zf: maximum redshift required
 *
 * Requires the final redshift of at least $z_f$ = @zf.
 *
 */
void 
nc_scalefactor_require_zf (NcScalefactor *a, const gdouble zf)
{
  NcScalefactorPrivate * const self = a->priv;
  if (zf > self->zf)
  {
    self->zf = zf;
    ncm_model_ctrl_force_update (self->ctrl);
    if (self->dist != NULL)
      nc_distance_require_zf (self->dist, zf);
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
  NcScalefactorPrivate * const self = a->priv;
  if (self->a0 != a0)
  {
    self->a0 = a0;
    ncm_model_ctrl_force_update (self->ctrl);
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
  NcScalefactorPrivate * const self = a->priv;
  if (self->reltol != reltol)
  {
    self->reltol = reltol;
    ncm_model_ctrl_force_update (self->ctrl);
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
  NcScalefactorPrivate * const self = a->priv;
  if (self->abstol != abstol)
  {
    self->abstol = abstol;
    ncm_model_ctrl_force_update (self->ctrl);
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
  NcScalefactorPrivate * const self = a->priv;
  if ((!self->sets_conf_norm && enable) || (self->sets_conf_norm && !enable))
  {
    self->sets_conf_norm = enable;
    ncm_model_ctrl_force_update (self->ctrl);
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
  NcScalefactorPrivate * const self = a->priv;
  return self->zf;
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
  NcScalefactorPrivate * const self = a->priv;
  return self->a0;
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
  NcScalefactorPrivate * const self = a->priv;
  return self->reltol;
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
  NcScalefactorPrivate * const self = a->priv;
  return self->abstol;
}

static void nc_scalefactor_init_cvode (NcScalefactor *a, NcHICosmo *cosmo);

static void
_nc_scalefactor_init (NcScalefactor *a, NcHICosmo *cosmo)
{
  NcScalefactorPrivate * const self = a->priv;
  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (cosmo);

  self->eta_i = nc_distance_conformal_time (self->dist, cosmo, self->zf);
  self->eta_f = nc_distance_conformal_time (self->dist, cosmo, 0.0);
  self->t_i   = nc_distance_cosmic_time (self->dist, cosmo, self->zf);

  NV_Ith_S (self->y, 0) = self->zf;
  NV_Ith_S (self->y, 1) = self->t_i;
  
  if (self->sets_conf_norm)
  {
    if (fabs (Omega_k0) < NC_SCALEFACTOR_OMEGA_K_ZERO)
    {
      self->a0 = 1.0;
    }
    else
    {
      const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo);
      self->a0 = RH_Mpc / sqrt (fabs (Omega_k0));
    }
  }
  nc_scalefactor_init_cvode (a, cosmo);

  return;
}

static void
nc_scalefactor_init_cvode (NcScalefactor *a, NcHICosmo *cosmo)
{
  NcScalefactorPrivate * const self = a->priv;
  gint flag;

  if (self->cvode_init)
  {
    flag = CVodeReInit (self->cvode, self->eta_i, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }
  else
  {
    flag = CVodeInit (self->cvode, &dz_deta_f, self->eta_i, self->y);
    NCM_CVODE_CHECK (&flag, "CVodeMalloc", 1, );
    self->cvode_init = TRUE;
  }

  flag = CVodeSStolerances (self->cvode, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

  flag = CVodeSetStopTime(self->cvode, self->eta_f);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetUserData (self->cvode, cosmo);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetMaxNumSteps (self->cvode, 50000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (self->cvode, &dz_deta_J);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
  
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
  NcScalefactorPrivate * const self = a->priv;
  return -ncm_spline_eval (self->a_eta, eta);
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
  NcScalefactorPrivate * const self = a->priv;
  return ncm_spline_eval (self->eta_a, - z);
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
  NcScalefactorPrivate * const self = a->priv;
  return ncm_spline_eval (self->eta_a, -(x - 1.0));
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
  NcScalefactorPrivate * const self = a->priv;
  gdouble mz;
  mz = ncm_spline_eval (self->a_eta, eta);
  return self->a0 / (1.0 - mz);
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
  NcScalefactorPrivate * const self = a->priv;
  return ncm_spline_eval (self->t_eta, eta);
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
  NcScalefactorPrivate * const self = a->priv;
  return ncm_spline_eval (self->eta_t, t);
}

static void
_nc_scalefactor_calc_spline (NcScalefactor *a, NcHICosmo *cosmo)
{
  NcScalefactorPrivate * const self = a->priv;
  GArray *eta_a, *mz_a, *t_a = NULL;
  gdouble mzi, last_eta;
  gint flag;

	if (!ncm_spline_is_empty (self->a_eta))
	{
		eta_a = ncm_vector_get_array (self->a_eta->xv);
		mz_a  = ncm_vector_get_array (self->a_eta->yv);
		g_array_set_size (eta_a, 0);
		g_array_set_size (mz_a, 0);
	}
	else
	{
		eta_a = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
		mz_a  = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
	}

  if (!ncm_spline_is_empty (self->t_eta))
  {
    t_a = ncm_vector_get_array (self->t_eta->yv);
    g_array_set_size (t_a, 0);
  }
  else
    t_a = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);

  
  last_eta = self->eta_i;
  mzi      = -self->zf;
  
  g_array_append_val (eta_a, self->eta_i);
  g_array_append_val (mz_a, mzi);
  g_array_append_val (t_a, self->t_i);

  while (TRUE)
  {
    gdouble eta;

    flag = CVode (self->cvode, self->eta_f, self->y, &eta, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "_nc_scalefactor_calc_spline[CVode]", 1, );

    if ((eta - last_eta) / last_eta > NC_SCALEFACTOR_MIN_ETA_STEP)
    {
      const gdouble mz = - NV_Ith_S (self->y, 0);
      const gdouble t  = NV_Ith_S (self->y, 1);

      g_array_append_val (eta_a, eta);
      g_array_append_val (mz_a, mz);
      g_array_append_val (t_a, t);
      last_eta = eta;
    }

    if (self->eta_f == eta)
      break;
  }
  
	if (fabs (g_array_index (mz_a, gdouble, mz_a->len - 1)) < 1e-10)
	{
		g_array_index (mz_a, gdouble, mz_a->len - 1) = 0.0;
	}
	else
		g_error ("_nc_scalefactor_calc_spline today redshift must be zero not % 20.15g.",
		         fabs (g_array_index (mz_a, gdouble, mz_a->len - 1)));

	ncm_spline_set_array (self->a_eta, eta_a, mz_a, TRUE);
	ncm_spline_set_array (self->eta_a, mz_a, eta_a, TRUE);
  g_array_unref (eta_a);
  g_array_unref (mz_a);

  ncm_spline_set_array (self->t_eta, eta_a, t_a, TRUE);
  ncm_spline_set_array (self->eta_t, t_a, eta_a, TRUE);
  g_array_unref (t_a);
  
  return;
}

static gint
dz_deta_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *cosmo = NC_HICOSMO (f_data);
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble x = 1.0 + z;
  const gdouble E = nc_hicosmo_E (cosmo, z);

  NV_Ith_S (ydot, 0) = - E;
  NV_Ith_S (ydot, 1) = 1.0 / x;
  return 0;
}

static gint
dz_deta_J (realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *cosmo     = NC_HICOSMO (jac_data);
  const gdouble z      = NV_Ith_S (y, 0);
  const gdouble x      = 1.0 + z;
  const gdouble E      = nc_hicosmo_E (cosmo, z);
  const gdouble dE2_dz = nc_hicosmo_dE2_dz (cosmo, z);
  
  SUN_DENSE_ACCESS (J, 0, 0) = - dE2_dz / (2.0 * E);
  SUN_DENSE_ACCESS (J, 0, 1) = 0.0;

  SUN_DENSE_ACCESS (J, 1, 0) = - 1.0 / (x * x);
  SUN_DENSE_ACCESS (J, 1, 1) = 0.0;
  
  return 0;
}

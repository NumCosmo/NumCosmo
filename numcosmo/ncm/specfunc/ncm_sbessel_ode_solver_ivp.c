/***************************************************************************
 *            ncm_sbessel_ode_solver_ivp.c
 *
 *  Mon Aug 24 2026
 *  Copyright 2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * NcmSBesselOdeSolverIVP:
 *
 * Explicit Runge--Kutta initial-value solver for the auxiliary spherical-Bessel
 * equation. The solver evolves the same polynomial form as
 * #NcmSBesselOdeSolver,
 * $$
 * x^2 u'' + [x^2 - \ell(\ell+1)]u = f(x)
 * $$
 * from a positive initial point $x_i$ with $u(x_i)=u'(x_i)=0$. This is the IVP counterpart of
 * #NcmSBesselOdeSolver and is intended as a proof-of-concept backend for
 * Levin-like integrations. The standard formulation evolves $u$ directly.
 * The Frobenius formulation evolves $u/(x/x_i)^{\ell+1}$, factoring out the
 * regular power-law solution in the non-oscillatory region.
 */


#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_sbessel_ode_solver_ivp.h"
#include "ncm/core/ncm_util.h"
#include "ncm_enum_types.h"

#include <float.h>
#include <math.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmSBesselOdeSolverIVPData
{
  NcmSBesselOdeSolverF f;
  gpointer user_data;
  NcmSBesselOdeSolverIVPMethod method;
  gdouble xi;
  guint ell;
} NcmSBesselOdeSolverIVPData;

typedef struct _NcmSBesselOdeSolverIVPPrivate
{
  gdouble reltol;
  gdouble abstol;
  NcmSBesselOdeSolverIVPMethod method;
  SUNContext sunctx;
  N_Vector y;
  gpointer arkode;
  NcmSBesselOdeSolverIVPData data;
} NcmSBesselOdeSolverIVPPrivate;

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_METHOD,
};

struct _NcmSBesselOdeSolverIVP
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSBesselOdeSolverIVP, ncm_sbessel_ode_solver_ivp, G_TYPE_OBJECT)

static int
_ncm_sbessel_ode_solver_ivp_rhs (sunrealtype x, N_Vector y, N_Vector ydot, gpointer user_data)
{
  NcmSBesselOdeSolverIVPData * const data = user_data;
  const gdouble ell                       = data->ell;
  const gdouble f                         = data->f != NULL ? data->f (data->user_data, x) : 0.0;
  const gdouble x2                        = x * x;

  switch (data->method)
  {
    case NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_STANDARD:
    {
      const gdouble u = NV_Ith_S (y, 0);

      NV_Ith_S (ydot, 0) = NV_Ith_S (y, 1);
      NV_Ith_S (ydot, 1) = ((ell * (ell + 1.0) - x2) * u + f) / x2;
      break;
    }
    case NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_FROBENIUS:
    {
      const gdouble p         = ell + 1.0;
      const gdouble scale_inv = pow (x / data->xi, -p);
      const gdouble v         = NV_Ith_S (y, 0);
      const gdouble dv        = NV_Ith_S (y, 1);

      NV_Ith_S (ydot, 0) = dv;
      NV_Ith_S (ydot, 1) = f * scale_inv / x2 - 2.0 * p * dv / x - v;
      break;
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }

  return isfinite (f) && isfinite (NV_Ith_S (ydot, 0)) && isfinite (NV_Ith_S (ydot, 1)) ? 0 : -1;
}

static void
ncm_sbessel_ode_solver_ivp_init (NcmSBesselOdeSolverIVP *solver)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  self->reltol = 1.0e-10;
  self->abstol = 1.0e-12;
  self->method = NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_STANDARD;
  self->arkode = NULL;

  if (SUNContext_Create (SUN_COMM_NULL, &self->sunctx))
    g_error ("ncm_sbessel_ode_solver_ivp_init: SUNContext_Create failed");

  self->y = N_VNew_Serial (2, self->sunctx);

  if (self->y == NULL)
    g_error ("ncm_sbessel_ode_solver_ivp_init: N_VNew_Serial failed");
}

static void
_ncm_sbessel_ode_solver_ivp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselOdeSolverIVP *solver = NCM_SBESSEL_ODE_SOLVER_IVP (object);

  switch (prop_id)
  {
    case PROP_RELTOL:
      ncm_sbessel_ode_solver_ivp_set_reltol (solver, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_sbessel_ode_solver_ivp_set_abstol (solver, g_value_get_double (value));
      break;
    case PROP_METHOD:
      ncm_sbessel_ode_solver_ivp_set_method (solver, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_ode_solver_ivp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselOdeSolverIVP *solver             = NCM_SBESSEL_ODE_SOLVER_IVP (object);
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  switch (prop_id)
  {
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, self->abstol);
      break;
    case PROP_METHOD:
      g_value_set_enum (value, self->method);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_ode_solver_ivp_finalize (GObject *object)
{
  NcmSBesselOdeSolverIVP *solver             = NCM_SBESSEL_ODE_SOLVER_IVP (object);
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  if (self->arkode != NULL)
    ERKStepFree (&self->arkode);

  if (self->y != NULL)
    N_VDestroy (self->y);

  SUNContext_Free (&self->sunctx);

  G_OBJECT_CLASS (ncm_sbessel_ode_solver_ivp_parent_class)->finalize (object);
}

static void
ncm_sbessel_ode_solver_ivp_class_init (NcmSBesselOdeSolverIVPClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_ode_solver_ivp_set_property;
  object_class->get_property = &_ncm_sbessel_ode_solver_ivp_get_property;
  object_class->finalize     = &_ncm_sbessel_ode_solver_ivp_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol", NULL, "Relative integration tolerance",
                                                        DBL_MIN, 1.0, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol", NULL, "Absolute integration tolerance",
                                                        DBL_MIN, 1.0, 1.0e-12,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method", NULL, "IVP numerical formulation",
                                                      NCM_TYPE_SBESSEL_ODE_SOLVER_IVP_METHOD,
                                                      NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_STANDARD,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sbessel_ode_solver_ivp_new:
 *
 * Creates a new #NcmSBesselOdeSolverIVP.
 *
 * Returns: (transfer full): a new #NcmSBesselOdeSolverIVP
 */
NcmSBesselOdeSolverIVP *
ncm_sbessel_ode_solver_ivp_new (void)
{
  return g_object_new (NCM_TYPE_SBESSEL_ODE_SOLVER_IVP, NULL);
}

/**
 * ncm_sbessel_ode_solver_ivp_ref:
 * @solver: a #NcmSBesselOdeSolverIVP
 *
 * Increases the reference count of @solver by one.
 *
 * Returns: (transfer full): the same @solver with increased reference count
 */
NcmSBesselOdeSolverIVP *
ncm_sbessel_ode_solver_ivp_ref (NcmSBesselOdeSolverIVP *solver)
{
  return g_object_ref (solver);
}

/**
 * ncm_sbessel_ode_solver_ivp_unref:
 * @solver: a #NcmSBesselOdeSolverIVP
 *
 * Decreases the reference count of @solver by one. If the reference count reaches
 * zero, the object is finalized and freed.
 */
void
ncm_sbessel_ode_solver_ivp_free (NcmSBesselOdeSolverIVP *solver)
{
  g_object_unref (solver);
}

/**
 * ncm_sbessel_ode_solver_ivp_clear:
 * @solver: a pointer to a #NcmSBesselOdeSolverIVP
 *
 * If @solver is not %NULL, decreases the reference count of the object by one and sets
 * @solver to %NULL. This is a convenience function for managing the lifetime of
 * #NcmSBesselOdeSolverIVP objects.
 *
 */
void
ncm_sbessel_ode_solver_ivp_clear (NcmSBesselOdeSolverIVP **solver)
{
  g_clear_object (solver);
}

/**
 * ncm_sbessel_ode_solver_ivp_set_reltol:
 * @solver: a #NcmSBesselOdeSolverIVP
 * @reltol: relative tolerance for the ODE solver
 *
 * Sets the relative tolerance for the ODE solver. The value must be in the range
 * (0, 1].
 *
 */
void
ncm_sbessel_ode_solver_ivp_set_reltol (NcmSBesselOdeSolverIVP *solver, gdouble reltol)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  g_return_if_fail (reltol > 0.0 && reltol <= 1.0);
  self->reltol = reltol;
}

/**
 * ncm_sbessel_ode_solver_ivp_get_reltol:
 * @solver: a #NcmSBesselOdeSolverIVP
 *
 * Gets the relative tolerance for the ODE solver.
 *
 * Returns: the relative tolerance
 */
gdouble
ncm_sbessel_ode_solver_ivp_get_reltol (NcmSBesselOdeSolverIVP *solver)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  return self->reltol;
}

/**
 * ncm_sbessel_ode_solver_ivp_set_abstol:
 * @solver: a #NcmSBesselOdeSolverIVP
 * @abstol: absolute tolerance for the ODE solver
 *
 * Sets the absolute tolerance for the ODE solver. The value must be in the range
 * (0, 1].
 *
 */
void
ncm_sbessel_ode_solver_ivp_set_abstol (NcmSBesselOdeSolverIVP *solver, gdouble abstol)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  g_return_if_fail (abstol > 0.0 && abstol <= 1.0);
  self->abstol = abstol;
}

/**
 * ncm_sbessel_ode_solver_ivp_get_abstol:
 * @solver: a #NcmSBesselOdeSolverIVP
 *
 * Gets the absolute tolerance for the ODE solver.
 *
 * Returns: the absolute tolerance
 */
gdouble
ncm_sbessel_ode_solver_ivp_get_abstol (NcmSBesselOdeSolverIVP *solver)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  return self->abstol;
}

/**
 * ncm_sbessel_ode_solver_ivp_set_method:
 * @solver: a #NcmSBesselOdeSolverIVP
 * @method: numerical formulation
 *
 * Sets the numerical formulation used to evolve the IVP.
 */
void
ncm_sbessel_ode_solver_ivp_set_method (NcmSBesselOdeSolverIVP *solver, NcmSBesselOdeSolverIVPMethod method)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  g_return_if_fail (method >= 0 && method < NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_LEN);
  self->method = method;
}

/**
 * ncm_sbessel_ode_solver_ivp_get_method:
 * @solver: a #NcmSBesselOdeSolverIVP
 *
 * Gets the numerical formulation used to evolve the IVP.
 *
 * Returns: the current #NcmSBesselOdeSolverIVPMethod
 */
NcmSBesselOdeSolverIVPMethod
ncm_sbessel_ode_solver_ivp_get_method (NcmSBesselOdeSolverIVP *solver)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);

  return self->method;
}

/**
 * ncm_sbessel_ode_solver_ivp_solve:
 * @solver: a #NcmSBesselOdeSolverIVP
 * @f: (nullable) (scope call): right-hand side $f(x)$ of the polynomial equation
 * @xi: positive initial point
 * @xf: final point, not smaller than @xi
 * @ell: angular momentum
 * @user_data: data passed to @f
 * @u: (out): value of $u(x_f)$
 * @du: (out): value of $u'(x_f)$
 *
 * Evolves the zero-initial-condition IVP from @xi to @xf. In the Frobenius
 * formulation, the internal state is $v = u/(x/x_i)^{\ell+1}$. Passing %NULL
 * for @f solves the homogeneous equation.
 */
void
ncm_sbessel_ode_solver_ivp_solve (NcmSBesselOdeSolverIVP *solver, NcmSBesselOdeSolverF f, gdouble xi, gdouble xf, guint ell, gpointer user_data, gdouble *u, gdouble *du)
{
  NcmSBesselOdeSolverIVPPrivate * const self = ncm_sbessel_ode_solver_ivp_get_instance_private (solver);
  sunrealtype x_attained                     = xi;
  gint flag;

  g_return_if_fail (NCM_IS_SBESSEL_ODE_SOLVER_IVP (solver));
  g_return_if_fail (xi > 0.0);
  g_return_if_fail (xf >= xi);
  g_return_if_fail (u != NULL);
  g_return_if_fail (du != NULL);

  self->data.f         = f;
  self->data.user_data = user_data;
  self->data.method    = self->method;
  self->data.xi        = xi;
  self->data.ell       = ell;

  NV_Ith_S (self->y, 0) = 0.0;
  NV_Ith_S (self->y, 1) = 0.0;

  if (xf == xi)
  {
    *u  = 0.0;
    *du = 0.0;

    return;
  }

  if (self->arkode == NULL)
  {
    self->arkode = ERKStepCreate (&_ncm_sbessel_ode_solver_ivp_rhs, xi, self->y, self->sunctx);

    if (self->arkode == NULL)
      g_error ("ncm_sbessel_ode_solver_ivp_solve: ERKStepCreate failed");
  }
  else
  {
    flag = ERKStepReInit (self->arkode, &_ncm_sbessel_ode_solver_ivp_rhs, xi, self->y);
    NCM_CVODE_CHECK (&flag, "ERKStepReInit", 1, );
  }

  flag = ERKStepSetUserData (self->arkode, &self->data);
  NCM_CVODE_CHECK (&flag, "ERKStepSetUserData", 1, );
  flag = ERKStepSStolerances (self->arkode, self->reltol, self->abstol);
  NCM_CVODE_CHECK (&flag, "ERKStepSStolerances", 1, );
  flag = ERKStepSetMaxNumSteps (self->arkode, 10000000);
  NCM_CVODE_CHECK (&flag, "ERKStepSetMaxNumSteps", 1, );
  flag = ERKStepEvolve (self->arkode, xf, self->y, &x_attained, ARK_NORMAL);
  NCM_CVODE_CHECK (&flag, "ERKStepEvolve", 1, );

  if (self->method == NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_FROBENIUS)
  {
    const gdouble p     = ell + 1.0;
    const gdouble scale = pow (xf / xi, p);
    const gdouble v     = NV_Ith_S (self->y, 0);
    const gdouble dv    = NV_Ith_S (self->y, 1);

    *u  = scale * v;
    *du = scale * (dv + p * v / xf);
  }
  else
  {
    *u  = NV_Ith_S (self->y, 0);
    *du = NV_Ith_S (self->y, 1);
  }
}

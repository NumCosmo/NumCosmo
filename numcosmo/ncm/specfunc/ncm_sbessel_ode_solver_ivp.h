/***************************************************************************
 *            ncm_sbessel_ode_solver_ivp.h
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
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NCM_SBESSEL_ODE_SOLVER_IVP_H_
#define _NCM_SBESSEL_ODE_SOLVER_IVP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.h>

G_BEGIN_DECLS

/**
 * NcmSBesselOdeSolverIVPMethod:
 * @NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_STANDARD: evolves $u$ and $u'$ directly.
 * @NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_FROBENIUS: evolves $u/(x/x_i)^{\ell+1}$
 *   and its derivative, removing the regular Frobenius power law.
 *
 * Numerical formulation used by #NcmSBesselOdeSolverIVP.
 */
typedef enum _NcmSBesselOdeSolverIVPMethod
{
  NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_STANDARD = 0,
  NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_FROBENIUS,
  /* < private > */
  NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_LEN, /*< skip >*/
} NcmSBesselOdeSolverIVPMethod;

#define NCM_TYPE_SBESSEL_ODE_SOLVER_IVP (ncm_sbessel_ode_solver_ivp_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselOdeSolverIVP, ncm_sbessel_ode_solver_ivp, NCM, SBESSEL_ODE_SOLVER_IVP, GObject)

NcmSBesselOdeSolverIVP *ncm_sbessel_ode_solver_ivp_new (void);
NcmSBesselOdeSolverIVP *ncm_sbessel_ode_solver_ivp_ref (NcmSBesselOdeSolverIVP *solver);
void ncm_sbessel_ode_solver_ivp_free (NcmSBesselOdeSolverIVP *solver);
void ncm_sbessel_ode_solver_ivp_clear (NcmSBesselOdeSolverIVP **solver);

void ncm_sbessel_ode_solver_ivp_set_reltol (NcmSBesselOdeSolverIVP *solver, gdouble reltol);
gdouble ncm_sbessel_ode_solver_ivp_get_reltol (NcmSBesselOdeSolverIVP *solver);
void ncm_sbessel_ode_solver_ivp_set_abstol (NcmSBesselOdeSolverIVP *solver, gdouble abstol);
gdouble ncm_sbessel_ode_solver_ivp_get_abstol (NcmSBesselOdeSolverIVP *solver);
void ncm_sbessel_ode_solver_ivp_set_method (NcmSBesselOdeSolverIVP *solver, NcmSBesselOdeSolverIVPMethod method);
NcmSBesselOdeSolverIVPMethod ncm_sbessel_ode_solver_ivp_get_method (NcmSBesselOdeSolverIVP *solver);

void ncm_sbessel_ode_solver_ivp_solve (NcmSBesselOdeSolverIVP *solver, NcmSBesselOdeSolverF f, gdouble xi, gdouble xf, guint ell, gpointer user_data, gdouble *u, gdouble *du);

G_END_DECLS

#endif /* _NCM_SBESSEL_ODE_SOLVER_IVP_H_ */


/***************************************************************************
 *            ncm_sbessel_ode_solver.h
 *
 *  Mon Jan 28 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
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
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NCM_SBESSEL_ODE_SOLVER_H_
#define _NCM_SBESSEL_ODE_SOLVER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_spectral.h>

G_BEGIN_DECLS

/**
 * NcmSBesselOdeSolverF:
 * @user_data: user data
 * @x: point to evaluate
 *
 * Function to be used as right-hand side of the ODE.
 *
 * Returns: the value of the function at @x
 */
typedef gdouble (*NcmSBesselOdeSolverF) (gpointer user_data, gdouble x);

#define NCM_TYPE_SBESSEL_ODE_SOLVER (ncm_sbessel_ode_solver_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselOdeSolver, ncm_sbessel_ode_solver, NCM, SBESSEL_ODE_SOLVER, GObject)

NcmSBesselOdeSolver *ncm_sbessel_ode_solver_new (gint l, gdouble a, gdouble b);
NcmSBesselOdeSolver *ncm_sbessel_ode_solver_ref (NcmSBesselOdeSolver *solver);

void ncm_sbessel_ode_solver_free (NcmSBesselOdeSolver *solver);
void ncm_sbessel_ode_solver_clear (NcmSBesselOdeSolver **solver);

void ncm_sbessel_ode_solver_set_l (NcmSBesselOdeSolver *solver, gint l);
gint ncm_sbessel_ode_solver_get_l (NcmSBesselOdeSolver *solver);

void ncm_sbessel_ode_solver_set_tolerance (NcmSBesselOdeSolver *solver, gdouble tol);
gdouble ncm_sbessel_ode_solver_get_tolerance (NcmSBesselOdeSolver *solver);

void ncm_sbessel_ode_solver_set_interval (NcmSBesselOdeSolver *solver, gdouble a, gdouble b);
void ncm_sbessel_ode_solver_get_interval (NcmSBesselOdeSolver *solver, gdouble *a, gdouble *b);

NcmVector *ncm_sbessel_ode_solver_solve (NcmSBesselOdeSolver *solver, NcmVector *rhs);
void ncm_sbessel_ode_solver_solve_endpoints (NcmSBesselOdeSolver *solver, NcmVector *rhs, gdouble *deriv_a, gdouble *deriv_b, gdouble *error);
NcmMatrix *ncm_sbessel_ode_solver_solve_batched (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin, guint n_l);
void ncm_sbessel_ode_solver_solve_endpoints_batched (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin, guint n_l, NcmMatrix *endpoints);

NcmMatrix *ncm_sbessel_ode_solver_get_operator_matrix (NcmSBesselOdeSolver *solver, gint nrows);
NcmMatrix *ncm_sbessel_ode_solver_get_operator_matrix_colmajor (NcmSBesselOdeSolver *solver, gint nrows);
NcmVector *ncm_sbessel_ode_solver_solve_dense (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint nrows);

NcmSpectral *ncm_sbessel_ode_solver_peek_spectral (NcmSBesselOdeSolver *solver);

G_END_DECLS

#endif /* _NCM_SBESSEL_ODE_SOLVER_H_ */


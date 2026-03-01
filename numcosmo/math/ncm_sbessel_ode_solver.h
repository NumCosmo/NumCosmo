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

/**
 * NcmSBesselOdeOperator:
 *
 * Opaque boxed type representing a spherical Bessel ODE operator.
 *
 * This structure encapsulates all problem-specific data for solving a spherical
 * Bessel ODE, including:
 * - Structural parameters: interval endpoints [a, b] and angular momentum range [ell_min, ell_max]
 * - Matrix storage: banded matrix rows and right-hand side vectors
 * - Diagonalization state: QR factorization and convergence information
 * - Tolerance: convergence criterion copied from the parent solver
 *
 * Operators are created from a #NcmSBesselOdeSolver via ncm_sbessel_ode_solver_create_operator()
 * and use reference counting for memory management. They can be reused for multiple solves
 * with the same parameters, or reset with different parameters using ncm_sbessel_ode_operator_reset().
 *
 * The operator supports both single and batched (multiple ell) solving modes, automatically
 * dispatching to optimized implementations based on the number of ell values.
 */
typedef struct _NcmSBesselOdeOperator NcmSBesselOdeOperator;

#define NCM_TYPE_SBESSEL_ODE_SOLVER (ncm_sbessel_ode_solver_get_type ())
#define NCM_TYPE_SBESSEL_OPERATOR (ncm_sbessel_ode_operator_get_type ())

G_DECLARE_FINAL_TYPE (NcmSBesselOdeSolver, ncm_sbessel_ode_solver, NCM, SBESSEL_ODE_SOLVER, GObject)
GType ncm_sbessel_ode_operator_get_type (void) G_GNUC_CONST;

NcmSBesselOdeSolver *ncm_sbessel_ode_solver_new (void);
NcmSBesselOdeSolver *ncm_sbessel_ode_solver_ref (NcmSBesselOdeSolver *solver);

void ncm_sbessel_ode_solver_free (NcmSBesselOdeSolver *solver);
void ncm_sbessel_ode_solver_clear (NcmSBesselOdeSolver **solver);

void ncm_sbessel_ode_solver_set_tolerance (NcmSBesselOdeSolver *solver, gdouble tol);
gdouble ncm_sbessel_ode_solver_get_tolerance (NcmSBesselOdeSolver *solver);

NcmMatrix *ncm_sbessel_ode_solver_get_operator_matrix (NcmSBesselOdeSolver *solver, const gdouble a, const gdouble b, guint ell, gint nrows);
NcmMatrix *ncm_sbessel_ode_solver_get_operator_matrix_colmajor (NcmSBesselOdeSolver *solver, const gdouble a, const gdouble b, guint ell, gint nrows);
NcmVector *ncm_sbessel_ode_solver_solve_dense (NcmSBesselOdeSolver *solver, const gdouble a, const gdouble b, guint ell, NcmVector *rhs, gint nrows);

NcmSpectral *ncm_sbessel_ode_solver_peek_spectral (NcmSBesselOdeSolver *solver);

NcmSBesselOdeOperator *ncm_sbessel_ode_solver_create_operator (NcmSBesselOdeSolver *solver, gdouble a, gdouble b, gint ell_min, gint ell_max);
NcmSBesselOdeOperator *ncm_sbessel_ode_operator_ref (NcmSBesselOdeOperator *op);
void ncm_sbessel_ode_operator_unref (NcmSBesselOdeOperator *op);
void ncm_sbessel_ode_operator_clear (NcmSBesselOdeOperator **op);

void ncm_sbessel_ode_operator_reset (NcmSBesselOdeOperator *op, gdouble a, gdouble b, gint ell_min, gint ell_max);

void ncm_sbessel_ode_operator_get_interval (NcmSBesselOdeOperator *op, gdouble *a, gdouble *b);
void ncm_sbessel_ode_operator_get_ell_range (NcmSBesselOdeOperator *op, gint *ell_min, gint *ell_max);
gdouble ncm_sbessel_ode_operator_get_tolerance (NcmSBesselOdeOperator *op);
glong ncm_sbessel_ode_operator_get_n_cols (NcmSBesselOdeOperator *op);

void ncm_sbessel_ode_operator_solve (NcmSBesselOdeOperator *op, GArray *rhs, GArray **solution, gsize *solution_len);
void ncm_sbessel_ode_operator_solve_endpoints (NcmSBesselOdeOperator *op, GArray *rhs, GArray **endpoints);


G_END_DECLS

#endif /* _NCM_SBESSEL_ODE_SOLVER_H_ */


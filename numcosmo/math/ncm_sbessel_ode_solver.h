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
 * @user_data: user-provided data pointer passed through from the caller
 * @x: evaluation point in the physical domain (not the mapped Chebyshev domain)
 *
 * Callback function type for the right-hand side forcing term $f(x)$ in the
 * modified spherical Bessel ODE. The function is evaluated at physical coordinates
 * $x \in [a,b]$ and should return $f(x)$ corresponding to the inhomogeneous term.
 *
 * Returns: the value $f(x)$ at the given point @x
 */
typedef gdouble (*NcmSBesselOdeSolverF) (gpointer user_data, gdouble x);

/**
 * NcmSBesselOdeOperator:
 *
 * Opaque boxed type for a configured spectral operator.
 *
 * Represents a fully configured two-point boundary value problem for the modified
 * spherical Bessel ODE over a specific interval $[a, b]$ and angular momentum range
 * $[\ell_{\min}, \ell_{\max}]$. Encapsulates:
 * - Problem parameters: interval endpoints, angular momentum range, and tolerance
 * - Discretized system: banded matrix representation and right-hand side storage
 * - Factorization state: adaptive QR decomposition and spectral truncation order
 *
 * Operators are created from a #NcmSBesselOdeSolver via ncm_sbessel_ode_solver_create_operator()
 * and managed via reference counting. Once created, an operator can be:
 * - Reused for multiple right-hand sides with the same $[a,b]$ and $[\ell_{\min}, \ell_{\max}]$
 * - Reconfigured for different parameters via ncm_sbessel_ode_operator_reset()
 * - Deallocated via ncm_sbessel_ode_operator_unref()
 *
 * Batched mode (multiple $\ell$ values) is automatically selected and optimized based on
 * $n_\ell = \ell_{\max} - \ell_{\min} + 1$.
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
gsize ncm_sbessel_ode_operator_get_operator_capacity (NcmSBesselOdeOperator *op);

void ncm_sbessel_ode_operator_solve (NcmSBesselOdeOperator *op, GArray *rhs, GArray **solution, gsize *solution_len);
void ncm_sbessel_ode_operator_solve_endpoints (NcmSBesselOdeOperator *op, GArray *rhs, GArray **endpoints);


G_END_DECLS

#endif /* _NCM_SBESSEL_ODE_SOLVER_H_ */


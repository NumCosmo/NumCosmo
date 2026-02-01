/***************************************************************************
 *            ncm_sbessel_ode_solver.c
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

/**
 * NcmSBesselOdeSolver:
 *
 * Spherical Bessel ODE solver using ultraspherical spectral methods.
 *
 * Solves the non-homogeneous spherical Bessel ODE mapped to the Chebyshev domain $[-1,1]$.
 * For an interval $[a,b]$, the mapping is $x = m + h\xi$ where $m = (a+b)/2$ and $h = (b-a)/2$.
 *
 * The ODE in the mapped coordinates is:
 * $$
 * \frac{(m+h\xi)^2}{h^2} \frac{d^2 y}{d\xi^2} + \frac{2(m+h\xi)}{h} \frac{dy}{d\xi} +
 * \left[(m+h\xi)^2 - l(l+1)\right] y = f(\xi)
 * $$
 * with Dirichlet boundary conditions $y(-1) = 0$ and $y(1) = 0$.
 *
 * The solver uses adaptive QR decomposition with the ultraspherical spectral method:
 * - The solution is expanded in Chebyshev polynomials $T_n(\xi)$
 * - Differential operators are represented in Gegenbauer $C^{(2)}_k$ basis ($\lambda=2$)
 * - Multiplication operators use exact polynomial multiplication formulas
 * - Adaptive QR determines the required series truncation automatically
 * - The matrix is built using specialized operator functions for maximum efficiency
 *
 * The method provides exponential convergence for smooth solutions and naturally
 * handles the coordinate transformation from the physical domain $[a,b]$ to $[-1,1]$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_ode_solver.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"
#include "ncm_cfg.h"

#include <math.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Operator bandwidth */
/* Number of rows to rotate lower bandwidth plus number of boundary conditions */
#define LOWER_BANDWIDTH 2
#define UPPER_BANDWIDTH 6
#define TOTAL_BANDWIDTH (LOWER_BANDWIDTH + UPPER_BANDWIDTH + 1)
#define ROWS_TO_ROTATE (LOWER_BANDWIDTH + 2)

/**
 * NcmSBesselOdeSolverRow:
 *
 * Internal structure for a matrix row.
 * Simple pointer to doubles with boundary condition coefficients.
 */
typedef struct _NcmSBesselOdeSolverRow
{
  gdouble *data;    /* Pointer to row data */
  gdouble bc_at_m1; /* Boundary condition row 1 coefficient */
  gdouble bc_at_p1; /* Boundary condition row 2 coefficient */
  glong col_index;  /* Column index of leftmost element in data */
} NcmSBesselOdeSolverRow;

typedef struct _NcmSBesselOdeSolverPrivate
{
  gint l;            /* Angular momentum parameter */
  gdouble tolerance; /* Convergence tolerance for adaptive QR */
  gint max_size;     /* Maximum matrix size */
  gdouble a;         /* Left endpoint of interval */
  gdouble b;         /* Right endpoint of interval */
  gdouble half_len;  /* (b-a)/2 - half length of interval */
  gdouble mid_point; /* (a+b)/2 - midpoint of interval */

  /* Matrix storage (adaptive) */
  GPtrArray *matrix_rows; /* Array of NcmSBesselOdeSolverRow* */

  /* Solution */
  NcmVector *solution; /* Result vector */
} NcmSBesselOdeSolverPrivate;

enum
{
  PROP_0,
  PROP_L,
  PROP_TOLERANCE,
  PROP_MAX_SIZE,
  PROP_INTERVAL,
};

struct _NcmSBesselOdeSolver
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSBesselOdeSolver, ncm_sbessel_ode_solver, G_TYPE_OBJECT)

/* Row operations */
static NcmSBesselOdeSolverRow *
_row_new (gdouble bc_at_m1, gdouble bc_at_p1, glong col_index)
{
  NcmSBesselOdeSolverRow *row = g_new (NcmSBesselOdeSolverRow, 1);

  row->data      = g_new0 (gdouble, TOTAL_BANDWIDTH);
  row->bc_at_m1  = bc_at_m1;
  row->bc_at_p1  = bc_at_p1;
  row->col_index = col_index;

  return row;
}

static void
_row_free (NcmSBesselOdeSolverRow *row)
{
  g_free (row->data);
  g_free (row);
}

/* Operator row computation */
static void _ncm_sbessel_compute_proj_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_x_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_x2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_d_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_x_d_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_d2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_x_d2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static void _ncm_sbessel_compute_x2_d2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff);
static gdouble _ncm_sbessel_bc_row (NcmSBesselOdeSolverRow *row, glong col_index);

static void
ncm_sbessel_ode_solver_init (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  self->l           = 0;
  self->tolerance   = NCM_SBESSEL_ODE_SOLVER_DEFAULT_TOLERANCE;
  self->max_size    = NCM_SBESSEL_ODE_SOLVER_DEFAULT_MAX_SIZE;
  self->a           = -1.0;
  self->b           = 1.0;
  self->half_len    = 1.0;
  self->mid_point   = 0.0;
  self->matrix_rows = g_ptr_array_new_with_free_func ((GDestroyNotify) _row_free);
  self->solution    = NULL;
}

static void
_ncm_sbessel_ode_solver_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselOdeSolver *solver = NCM_SBESSEL_ODE_SOLVER (object);

  g_return_if_fail (NCM_IS_SBESSEL_ODE_SOLVER (object));

  switch (prop_id)
  {
    case PROP_L:
      ncm_sbessel_ode_solver_set_l (solver, g_value_get_int (value));
      break;
    case PROP_TOLERANCE:
      ncm_sbessel_ode_solver_set_tolerance (solver, g_value_get_double (value));
      break;
    case PROP_MAX_SIZE:
      ncm_sbessel_ode_solver_set_max_size (solver, g_value_get_int (value));
      break;
    case PROP_INTERVAL:
    {
      NcmDTuple2 *interval = g_value_get_boxed (value);

      if (interval == NULL)
        g_error ("_ncm_sbessel_ode_solver_set_property: interval is NULL");

      ncm_sbessel_ode_solver_set_interval (solver, interval->elements[0], interval->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_ode_solver_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselOdeSolver *solver             = NCM_SBESSEL_ODE_SOLVER (object);
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  g_return_if_fail (NCM_IS_SBESSEL_ODE_SOLVER (object));

  switch (prop_id)
  {
    case PROP_L:
      g_value_set_int (value, self->l);
      break;
    case PROP_TOLERANCE:
      g_value_set_double (value, self->tolerance);
      break;
    case PROP_MAX_SIZE:
      g_value_set_int (value, self->max_size);
      break;
    case PROP_INTERVAL:
    {
      g_value_take_boxed (value, ncm_dtuple2_new (self->a, self->b));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_ode_solver_dispose (GObject *object)
{
  NcmSBesselOdeSolver *solver             = NCM_SBESSEL_ODE_SOLVER (object);
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  ncm_vector_clear (&self->solution);

  /* Chain up */
  G_OBJECT_CLASS (ncm_sbessel_ode_solver_parent_class)->dispose (object);
}

static void
_ncm_sbessel_ode_solver_finalize (GObject *object)
{
  /* Chain up */
  G_OBJECT_CLASS (ncm_sbessel_ode_solver_parent_class)->finalize (object);
}

static void
ncm_sbessel_ode_solver_class_init (NcmSBesselOdeSolverClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_ode_solver_set_property;
  object_class->get_property = &_ncm_sbessel_ode_solver_get_property;
  object_class->dispose      = &_ncm_sbessel_ode_solver_dispose;
  object_class->finalize     = &_ncm_sbessel_ode_solver_finalize;

  /**
   * NcmSBesselOdeSolver:l:
   *
   * Angular momentum quantum number l in the spherical Bessel ODE.
   */
  g_object_class_install_property (object_class,
                                   PROP_L,
                                   g_param_spec_int ("l",
                                                     NULL,
                                                     "Angular momentum parameter",
                                                     0, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselOdeSolver:tolerance:
   *
   * Convergence tolerance for adaptive QR decomposition.
   */
  g_object_class_install_property (object_class,
                                   PROP_TOLERANCE,
                                   g_param_spec_double ("tolerance",
                                                        NULL,
                                                        "Convergence tolerance",
                                                        0.0, 1.0, NCM_SBESSEL_ODE_SOLVER_DEFAULT_TOLERANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselOdeSolver:max-size:
   *
   * Maximum matrix size for adaptive algorithm.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_SIZE,
                                   g_param_spec_int ("max-size",
                                                     NULL,
                                                     "Maximum matrix size",
                                                     1, G_MAXINT, NCM_SBESSEL_ODE_SOLVER_DEFAULT_MAX_SIZE,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselOdeSolver:interval:
   *
   * Interval [a, b] where the ODE is solved.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTERVAL,
                                   g_param_spec_boxed ("interval",
                                                       NULL,
                                                       "Interval where the ODE is solved",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sbessel_ode_solver_new:
 * @l: angular momentum parameter
 * @a: left endpoint of interval
 * @b: right endpoint of interval
 *
 * Creates a new #NcmSBesselOdeSolver for the given angular momentum @l
 * and interval [a, b] where the ODE is solved.
 *
 * Returns: (transfer full): a new #NcmSBesselOdeSolver
 */
NcmSBesselOdeSolver *
ncm_sbessel_ode_solver_new (gint l, gdouble a, gdouble b)
{
  NcmDTuple2 *interval        = ncm_dtuple2_new (a, b);
  NcmSBesselOdeSolver *solver = g_object_new (NCM_TYPE_SBESSEL_ODE_SOLVER,
                                              "l", l,
                                              "interval", interval,
                                              NULL);

  ncm_dtuple2_free (interval);

  return solver;
}

/**
 * ncm_sbessel_ode_solver_ref:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Increases the reference count of @solver by one.
 *
 * Returns: (transfer full): @solver
 */
NcmSBesselOdeSolver *
ncm_sbessel_ode_solver_ref (NcmSBesselOdeSolver *solver)
{
  return g_object_ref (solver);
}

/**
 * ncm_sbessel_ode_solver_free:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Decreases the reference count of @solver by one.
 */
void
ncm_sbessel_ode_solver_free (NcmSBesselOdeSolver *solver)
{
  g_object_unref (solver);
}

/**
 * ncm_sbessel_ode_solver_clear:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Decreases the reference count of *@solver by one and sets *@solver to NULL.
 */
void
ncm_sbessel_ode_solver_clear (NcmSBesselOdeSolver **solver)
{
  g_clear_object (solver);
}

/**
 * ncm_sbessel_ode_solver_set_l:
 * @solver: a #NcmSBesselOdeSolver
 * @l: angular momentum parameter
 *
 * Sets the angular momentum parameter @l.
 */
void
ncm_sbessel_ode_solver_set_l (NcmSBesselOdeSolver *solver, gint l)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  if (self->l != l)
    self->l = l;
}

/**
 * ncm_sbessel_ode_solver_get_l:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Gets the angular momentum parameter l.
 *
 * Returns: the angular momentum parameter
 */
gint
ncm_sbessel_ode_solver_get_l (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return self->l;
}

/**
 * ncm_sbessel_ode_solver_set_tolerance:
 * @solver: a #NcmSBesselOdeSolver
 * @tol: convergence tolerance
 *
 * Sets the convergence tolerance for adaptive QR.
 */
void
ncm_sbessel_ode_solver_set_tolerance (NcmSBesselOdeSolver *solver, gdouble tol)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  self->tolerance = tol;
}

/**
 * ncm_sbessel_ode_solver_get_tolerance:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Gets the convergence tolerance.
 *
 * Returns: the convergence tolerance
 */
gdouble
ncm_sbessel_ode_solver_get_tolerance (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return self->tolerance;
}

/**
 * ncm_sbessel_ode_solver_set_max_size:
 * @solver: a #NcmSBesselOdeSolver
 * @max_size: maximum matrix size
 *
 * Sets the maximum matrix size for adaptive algorithm.
 */
void
ncm_sbessel_ode_solver_set_max_size (NcmSBesselOdeSolver *solver, gint max_size)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  self->max_size = max_size;
}

/**
 * ncm_sbessel_ode_solver_get_max_size:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Gets the maximum matrix size.
 *
 * Returns: the maximum matrix size
 */
gint
ncm_sbessel_ode_solver_get_max_size (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return self->max_size;
}

/**
 * ncm_sbessel_ode_solver_set_interval:
 * @solver: a #NcmSBesselOdeSolver
 * @a: left endpoint
 * @b: right endpoint
 *
 * Sets the interval [a, b] where the ODE is solved.
 * Computes and stores (b-a)/2 and (a+b)/2.
 */
void
ncm_sbessel_ode_solver_set_interval (NcmSBesselOdeSolver *solver, gdouble a, gdouble b)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  g_assert_cmpfloat (a, <, b);

  self->a         = a;
  self->b         = b;
  self->half_len  = (b - a) / 2.0;
  self->mid_point = (a + b) / 2.0;
}

/**
 * ncm_sbessel_ode_solver_get_interval:
 * @solver: a #NcmSBesselOdeSolver
 * @a: (out): left endpoint
 * @b: (out): right endpoint
 *
 * Gets the interval [a, b] where the ODE is solved.
 */
void
ncm_sbessel_ode_solver_get_interval (NcmSBesselOdeSolver *solver, gdouble *a, gdouble *b)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  *a = self->a;
  *b = self->b;
}

/**
 * ncm_sbessel_ode_solver_peek_solution:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Gets the solution vector (Chebyshev coefficients) from the last solve.
 * The returned vector is owned by @solver and should not be freed.
 *
 * Returns: (transfer none) (nullable): the solution vector or NULL if no solution computed
 */
NcmVector *
ncm_sbessel_ode_solver_peek_solution (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return self->solution;
}

/**
 * ncm_sbessel_ode_solver_get_solution_size:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Gets the size of the solution vector from the last solve.
 *
 * Returns: the solution size or 0 if no solution computed
 */
gint
ncm_sbessel_ode_solver_get_solution_size (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return (self->solution != NULL) ? (gint) ncm_vector_len (self->solution) : 0;
}

/**
 * _ncm_sbessel_ode_solver_row_add:
 * @row: row structure to update
 * @col: target column index in the global matrix
 * @value: value to add to the entry
 *
 * Helper function to add a value to a row entry at the specified column.
 * Handles offset calculation and bounds checking.
 *
 * The function computes the offset from row->col_index (leftmost column stored)
 * to the target column, validates it's within allocated space, and accumulates
 * the value. This allows building rows via linear combinations of operators.
 */
static void
_ncm_sbessel_ode_solver_row_add (NcmSBesselOdeSolverRow *row, glong col, gdouble value)
{
  const glong offset = col - row->col_index;

  /* Column must be at or to the right of col_index */
  g_assert_cmpint (offset, >=, 0);

  if (offset < TOTAL_BANDWIDTH)
    row->data[offset] += value;  /* Accumulate (supports += pattern) */
  else
    g_error ("_ncm_sbessel_ode_solver_row_add: offset %ld exceeds TOTAL_BANDWIDTH %d", offset, TOTAL_BANDWIDTH);
}

/**
 * _ncm_sbessel_compute_proj_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the projection operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients (ultraspherical basis with $\lambda=2$).
 * This is the identity operator expressed in different bases.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: Gegenbauer $C^{(2)}_k(x)$ basis coefficient (row k)
 *
 * Mathematical formula:
 * $$
 * C^{(2)}_k = \frac{1}{2} c_0 \delta_{k,0} +
 * \frac{c_k}{2(k+1)} - \frac{(k+2) c_{k+2}}{(k+1)(k+3)} + \frac{c_{k+4}}{2(k+3)}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k: coeff * 1/(2*(k+1))
 * - column k+2: coeff * -(k+2)/((k+1)*(k+3))
 * - column k+4: coeff * 1/(2*(k+3))
 * - For k=0 only: additional value coeff * 1/2 at column 0
 */
static void
_ncm_sbessel_compute_proj_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* General formula: three entries with rational function coefficients */
  const gdouble value_0 = coeff / (2.0 * (kd + 1.0));
  const gdouble value_1 = -coeff * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = coeff / (2.0 * (kd + 3.0));

  _ncm_sbessel_ode_solver_row_add (row, k, value_0);
  _ncm_sbessel_ode_solver_row_add (row, k + 2, value_1);
  _ncm_sbessel_ode_solver_row_add (row, k + 4, value_2);

  /* Special case: k=0 has an additional contribution */
  if (k == 0)
    _ncm_sbessel_ode_solver_row_add (row, 0, coeff * 0.5);
}

/**
 * _ncm_sbessel_compute_x_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the multiplication by x operator that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: Gegenbauer $C^{(2)}_k(x)$ basis coefficient for $x \cdot f(x)$ (row k)
 *
 * Mathematical formula:
 * $$
 * (x \cdot f)^{(2)}_k = \frac{\theta(k-1) c_{k-1}}{4(k+1)} - \frac{c_{k+1}}{4(k+3)} -
 * \frac{c_{k+3}}{4(k+1)} + \frac{c_{k+5}}{4(k+3)}
 * $$
 * plus special contributions: $\frac{c_1}{4}\delta_{k,0}$ and $\frac{c_0}{8}\delta_{k,1}$,
 * where $c_n$ are the input Chebyshev coefficients and $\theta$ is the Heaviside function.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - If k >= 1: column k-1: coeff * 1/(4*(k+1))
 * - column k+1: coeff * -1/(4*(k+3))
 * - column k+3: coeff * -1/(4*(k+1))
 * - column k+5: coeff * 1/(4*(k+3))
 * - For k=0 only: additional value coeff * 1/4 at column 1
 * - For k=1 only: additional value coeff * 1/8 at column 0
 */
static void
_ncm_sbessel_compute_x_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* General formula: up to four entries */
  if (k >= 1)
    _ncm_sbessel_ode_solver_row_add (row, k - 1, coeff / (4.0 * (kd + 1.0)));

  _ncm_sbessel_ode_solver_row_add (row, k + 1, -coeff / (4.0 * (kd + 3.0)));
  _ncm_sbessel_ode_solver_row_add (row, k + 3, -coeff / (4.0 * (kd + 1.0)));
  _ncm_sbessel_ode_solver_row_add (row, k + 5, coeff / (4.0 * (kd + 3.0)));

  /* Special cases */
  if (k == 0)
    _ncm_sbessel_ode_solver_row_add (row, 1, coeff * 0.25);  /* Additional 1/4 at column 1 */
  else if (k == 1)
    _ncm_sbessel_ode_solver_row_add (row, 0, coeff * 0.125);  /* Additional 1/8 at column 0 */
}

/**
 * _ncm_sbessel_compute_x2_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the multiplication by $x^2$ operator that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: Gegenbauer $C^{(2)}_k(x)$ basis coefficient for $x^2 \cdot f(x)$ (row k)
 *
 * Mathematical formula:
 * $$
 * (x^2 \cdot f)^{(2)}_k = \frac{\theta(k-2) c_{k-2}}{8(k+1)} + \frac{c_k}{4(k+1)(k+3)} -
 * \frac{(k+2) c_{k+2}}{4(k+1)(k+3)} - \frac{c_{k+4}}{4(k+1)(k+3)} + \frac{c_{k+6}}{8(k+3)}
 * $$
 * plus special contributions at k=0,1,2, where $c_n$ are the input Chebyshev
 * coefficients and $\theta$ is the Heaviside function.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - If k >= 2: column k-2: coeff * 1/(8*(k+1))
 * - column k: coeff * 1/(4*(k+1)*(k+3))
 * - column k+2: coeff * -(k+2)/(4*(k+1)*(k+3))
 * - column k+4: coeff * -1/(4*(k+1)*(k+3))
 * - column k+6: coeff * 1/(8*(k+3))
 * - For k=0 only: additional values coeff * 1/12 at column 0 and coeff * 1/8 at column 2
 * - For k=1 only: additional value coeff * 1/16 at column 1
 * - For k=2 only: additional value coeff * 1/24 at column 0
 */
static void
_ncm_sbessel_compute_x2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd  = (gdouble) k;
  const gdouble kp2 = kd + 2.0;

  /* General formula: up to five entries */
  if (k >= 2)
    _ncm_sbessel_ode_solver_row_add (row, k - 2, coeff / (8.0 * (kd + 1.0)));

  _ncm_sbessel_ode_solver_row_add (row, k, coeff / (4.0 * (kd + 1.0) * (kd + 3.0)));
  _ncm_sbessel_ode_solver_row_add (row, k + 2, -coeff * kp2 / (4.0 * (kd + 1.0) * (kd + 3.0)));
  _ncm_sbessel_ode_solver_row_add (row, k + 4, -coeff / (4.0 * (kd + 1.0) * (kd + 3.0)));
  _ncm_sbessel_ode_solver_row_add (row, k + 6, coeff / (8.0 * (kd + 3.0)));

  /* Special cases */
  if (k == 0)
  {
    _ncm_sbessel_ode_solver_row_add (row, 0, coeff / 12.0); /* Additional 1/12 at column 0 */
    _ncm_sbessel_ode_solver_row_add (row, 2, coeff / 8.0);  /* Additional 1/8 at column 2 */
  }
  else if (k == 1)
  {
    _ncm_sbessel_ode_solver_row_add (row, 1, coeff / 16.0); /* Additional 1/16 at column 1 */
  }
  else if (k == 2)
  {
    _ncm_sbessel_ode_solver_row_add (row, 0, coeff / 24.0); /* Additional 1/24 at column 0 */
  }
}

/**
 * _ncm_sbessel_compute_d_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the first derivative operator $\frac{d}{dx}$ that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, f' \rangle$ - projection of $f'$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, f' \rangle = c_{k+1} - c_{k+3}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k+1: coeff * 1.0
 * - column k+3: coeff * (-1.0)
 */
static void
_ncm_sbessel_compute_d_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  /* Two non-zero entries with opposite signs */
  _ncm_sbessel_ode_solver_row_add (row, k + 1, coeff);
  _ncm_sbessel_ode_solver_row_add (row, k + 3, -coeff);
}

/**
 * _ncm_sbessel_compute_x_d_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the $x \cdot \frac{d}{dx}$ operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, x \cdot f' \rangle$ - projection of $x \cdot f'$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, x \cdot f' \rangle = \frac{k c_k}{2(k+1)} + \frac{(k+2) c_{k+2}}{(k+1)(k+3)} -
 * \frac{(k+4) c_{k+4}}{2(k+3)}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k: coeff * k/(2*(k+1))
 * - column k+2: coeff * (k+2)/((k+1)*(k+3))
 * - column k+4: coeff * -(k+4)/(2*(k+3))
 */
static void
_ncm_sbessel_compute_x_d_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* Three non-zero entries with rational function coefficients */
  const gdouble value_0 = coeff * kd / (2.0 * (kd + 1.0));
  const gdouble value_1 = coeff * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = -coeff * (kd + 4.0) / (2.0 * (kd + 3.0));

  _ncm_sbessel_ode_solver_row_add (row, k, value_0);
  _ncm_sbessel_ode_solver_row_add (row, k + 2, value_1);
  _ncm_sbessel_ode_solver_row_add (row, k + 4, value_2);
}

/**
 * _ncm_sbessel_compute_d2_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the second derivative operator $\frac{d^2}{dx^2}$ that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, f'' \rangle$ - projection of $f''$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, f'' \rangle = 2(k+2) c_{k+2}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k+2: coeff * 2*(k+2) (single non-zero entry)
 */
static void
_ncm_sbessel_compute_d2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  _ncm_sbessel_ode_solver_row_add (row, k + 2, coeff * 2.0 * (kd + 2.0));
}

/**
 * _ncm_sbessel_compute_x_d2_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the $x \cdot \frac{d^2}{dx^2}$ operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, x \cdot f'' \rangle$ - projection of $x \cdot f''$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, x \cdot f'' \rangle = k c_{k+1} + (k+4) c_{k+3}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k+1: coeff * k
 * - column k+3: coeff * (k+4)
 */
static void
_ncm_sbessel_compute_x_d2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* Two non-zero entries in this row */
  _ncm_sbessel_ode_solver_row_add (row, k + 1, coeff * kd);
  _ncm_sbessel_ode_solver_row_add (row, k + 3, coeff * (kd + 4.0));
}

/**
 * _ncm_sbessel_compute_x2_d2_row:
 * @row: row structure to fill/update
 * @k: row index (0-based)
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the $x^2 \cdot \frac{d^2}{dx^2}$ operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, x^2 \cdot f'' \rangle$ - projection of $x^2 \cdot f''$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, x^2 \cdot f'' \rangle = \frac{k(k-1) c_k}{2(k+1)} +
 * \frac{(k+2)((k+2)^2-3) c_{k+2}}{(k+1)(k+3)} + \frac{(k+4)(k+5) c_{k+4}}{2(k+3)}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k: coeff * k*(k-1)/(2*(k+1))
 * - column k+2: coeff * (k+2)*((k+2)^2 - 3)/((k+1)*(k+3))
 * - column k+4: coeff * (k+4)*(k+5)/(2*(k+3))
 */
static void
_ncm_sbessel_compute_x2_d2_row (NcmSBesselOdeSolverRow *row, glong k, gdouble coeff)
{
  const gdouble kd  = (gdouble) k;
  const gdouble kp2 = kd + 2.0;

  /* Three non-zero entries with rational function coefficients */
  const gdouble value_0 = coeff * kd * (kd - 1.0) / (2.0 * (kd + 1.0));
  const gdouble value_1 = coeff * kp2 * (kp2 * kp2 - 3.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = coeff * (kd + 4.0) * (kd + 5.0) / (2.0 * (kd + 3.0));

  _ncm_sbessel_ode_solver_row_add (row, k, value_0);
  _ncm_sbessel_ode_solver_row_add (row, k + 2, value_1);
  _ncm_sbessel_ode_solver_row_add (row, k + 4, value_2);
}

static gdouble
_ncm_sbessel_bc_row (NcmSBesselOdeSolverRow *row, glong col_index)
{
  const gdouble bc_at_m1 = (col_index % 2) == 0 ? row->bc_at_m1 : -row->bc_at_m1;
  const gdouble bc_at_p1 = row->bc_at_p1;

  return (bc_at_m1 + bc_at_p1);
}

/**
 * _ncm_sbessel_create_row_bc_at_m1:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Creates boundary condition row for u(-1) = 0.
 * All data entries are zero, bc_at_m1 = 1.0, bc_at_p1 = 0.0.
 */
static NcmSBesselOdeSolverRow *
_ncm_sbessel_create_row_bc_at_m1 (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverRow *row = _row_new (1.0, 0.0, 0);

  /* All data entries remain zero (allocated with g_new0) */

  return row;
}

/**
 * _ncm_sbessel_create_row_bc_at_p1:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Creates boundary condition row for u(+1) = 0.
 * All data entries are zero, bc_at_m1 = 0.0, bc_at_p1 = 1.0.
 */
static NcmSBesselOdeSolverRow *
_ncm_sbessel_create_row_bc_at_p1 (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverRow *row = _row_new (0.0, 1.0, 0);

  /* All data entries remain zero (allocated with g_new0) */

  return row;
}

/**
 * _ncm_sbessel_create_row_operator:
 * @solver: a #NcmSBesselOdeSolver
 * @row_index: row index (0-based, corresponding to k in the Gegenbauer basis)
 *
 * Creates differential operator row for the given index k in the Gegenbauer $C^{(2)}_k$ basis.
 * Builds the full ODE operator using the formula:
 * $$\langle C^{(2)}_k, L[f] \rangle$$ where
 * $$L = \frac{(m+hx)^2}{h^2} \frac{d^2}{dx^2} + \frac{2(m+hx)}{h} \frac{d}{dx} +
 * (m+hx)^2 - l(l+1)$$
 *
 * The operator functions handle all special cases for low k values internally,
 * so this single function works for all k.
 */
static NcmSBesselOdeSolverRow *
_ncm_sbessel_create_row_operator (NcmSBesselOdeSolver *solver, glong row_index)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble m                         = self->mid_point;
  const gdouble h                         = self->half_len;
  const gdouble h2                        = h * h;
  const gdouble m2                        = m * m;
  const gint l                            = self->l;
  const glong k                           = row_index;
  const glong left_col                    = GSL_MAX (0, k - 2); /* Leftmost column index */
  NcmSBesselOdeSolverRow *row             = _row_new (0.0, 0.0, left_col);

  /* Second derivative term: (m^2/h^2) d^2 + (2m/h) x d^2 + x^2 d^2 */
  _ncm_sbessel_compute_d2_row (row, k, m2 / h2);
  _ncm_sbessel_compute_x_d2_row (row, k, 2.0 * m / h);
  _ncm_sbessel_compute_x2_d2_row (row, k, 1.0);

  /* First derivative term: (2m/h) d + 2 x d */
  _ncm_sbessel_compute_d_row (row, k, 2.0 * m / h);
  _ncm_sbessel_compute_x_d_row (row, k, 2.0);

  /* Identity term: (m^2 - l(l+1)) I + 2m h x + h^2 x^2 */
  _ncm_sbessel_compute_proj_row (row, k, m2 - (gdouble) (l * (l + 1)));
  _ncm_sbessel_compute_x_row (row, k, 2.0 * m * h);
  _ncm_sbessel_compute_x2_row (row, k, h2);

  return row;
}

/**
 * _ncm_sbessel_create_row:
 * @solver: a #NcmSBesselOdeSolver
 * @row_index: row index (first two rows are boundary conditions, rest are operators)
 *
 * Creates the appropriate row based on the row index.
 * Rows 0-1 are boundary conditions, rows >= 2 are differential operators.
 */
static NcmSBesselOdeSolverRow *
_ncm_sbessel_create_row (NcmSBesselOdeSolver *solver, glong row_index)
{
  if (row_index == 0)
    return _ncm_sbessel_create_row_bc_at_m1 (solver);
  else if (row_index == 1)
    return _ncm_sbessel_create_row_bc_at_p1 (solver);
  else
    return _ncm_sbessel_create_row_operator (solver, row_index - 2);
}

/**
 * _ncm_sbessel_apply_givens:
 * @solver: a #NcmSBesselOdeSolver
 * @row1: first row index (pivot row)
 * @row2: second row index (row to eliminate)
 * @c: right-hand side vector
 *
 * Applies Givens rotation to eliminate the entry at (row2, col=row1).
 * The rotation transforms the rows as:
 *   new_row1 = c*row1 + s*row2
 *   new_row2 = -s*row1 + c*row2
 * where c and s are chosen to zero out element (row2, row1).
 */
static void
_ncm_sbessel_apply_givens (NcmSBesselOdeSolver *solver, glong pivot_col, glong row1, glong row2, GArray *c)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  NcmSBesselOdeSolverRow *r1              = g_ptr_array_index (self->matrix_rows, row1);
  NcmSBesselOdeSolverRow *r2              = g_ptr_array_index (self->matrix_rows, row2);
  gdouble a_val                           = 0.0;
  gdouble b_val                           = 0.0;
  gdouble norm;

  g_assert_cmpuint (r1->col_index, ==, pivot_col); /* Pivot row aligned with its index */
  g_assert_cmpuint (r2->col_index, ==, pivot_col); /* Row to eliminate aligned with its index */

  a_val = r1->data[0] + _ncm_sbessel_bc_row (r1, pivot_col);
  b_val = r2->data[0] + _ncm_sbessel_bc_row (r2, pivot_col);

  norm = hypot (a_val, b_val);

  if (norm < 1.0e-100)
  {
    /* Entry already zero - just shift r2 without rotation */
    for (glong i = 1; i < TOTAL_BANDWIDTH; i++)
    {
      r2->data[i - 1] = r2->data[i];
    }

    r2->data[TOTAL_BANDWIDTH - 1] = 0.0; /* Last entry shifted out is now zero */
    r2->col_index++;

    return;
  }

  {
    /* Compute Givens rotation coefficients */
    const gdouble cos_theta = a_val / norm;
    const gdouble sin_theta = b_val / norm;
    const gdouble rhs1      = g_array_index (c, gdouble, row1);
    const gdouble rhs2      = g_array_index (c, gdouble, row2);
    const gdouble bc1_m1    = r1->bc_at_m1;
    const gdouble bc2_m1    = r2->bc_at_m1;
    const gdouble bc1_p1    = r1->bc_at_p1;
    const gdouble bc2_p1    = r2->bc_at_p1;
    glong i;

    g_array_index (c, gdouble, row1) = cos_theta * rhs1 + sin_theta * rhs2;
    g_array_index (c, gdouble, row2) = -sin_theta * rhs1 + cos_theta * rhs2;
    r1->bc_at_m1                     = cos_theta * bc1_m1 + sin_theta * bc2_m1;
    r2->bc_at_m1                     = -sin_theta * bc1_m1 + cos_theta * bc2_m1;
    r1->bc_at_p1                     = cos_theta * bc1_p1 + sin_theta * bc2_p1;
    r2->bc_at_p1                     = -sin_theta * bc1_p1 + cos_theta * bc2_p1;

    /* Updated leading entry in row1 */
    r1->data[0] = cos_theta * r1->data[0] + sin_theta * r2->data[0];

    /* Apply rotation to all matrix entries */
    for (i = 1; i < TOTAL_BANDWIDTH; i++)
    {
      const gdouble e1     = r1->data[i];
      const gdouble e2     = r2->data[i];
      const gdouble new_e1 = cos_theta * e1 + sin_theta * e2;
      const gdouble new_e2 = -sin_theta * e1 + cos_theta * e2;

      r1->data[i]     = new_e1;
      r2->data[i - 1] = new_e2;
    }

    r2->data[TOTAL_BANDWIDTH - 1] = 0.0; /* Last entry shifted out is now zero */
    r2->col_index++;                     /* Shift row2's column index right by 1 */
  }
}

/**
 * ncm_sbessel_ode_solver_solve:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 *
 * Solves the ODE using adaptive QR decomposition with ultraspherical spectral methods.
 * The algorithm grows the matrix size until convergence is achieved (error < tolerance).
 *
 * Returns: (transfer full): solution vector (Chebyshev coefficients)
 */
NcmVector *
ncm_sbessel_ode_solver_solve (NcmSBesselOdeSolver *solver, NcmVector *rhs)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const guint rhs_len                     = ncm_vector_len (rhs);
  GArray *c                               = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), rhs_len);
  glong col                               = 0;
  glong i, last_row_index;

  /* Clear previous solution */
  g_ptr_array_set_size (self->matrix_rows, 0);
  ncm_vector_clear (&self->solution);

  /* Add boundary condition rows */
  for (i = 0; i < ROWS_TO_ROTATE + 1; i++)
    g_ptr_array_add (self->matrix_rows, _ncm_sbessel_create_row (solver, i));

  last_row_index = ROWS_TO_ROTATE;

  /* Copy RHS to working array */
  for (i = 0; i < rhs_len; i++)
  {
    const gdouble rhs_i = ncm_vector_get (rhs, i);

    g_array_append_val (c, rhs_i);
  }

  for (i = 0; i < ROWS_TO_ROTATE; i++)
  {
    const gdouble zero = 0.0;

    g_array_append_val (c, zero);
  }

  for (col = 0; col < rhs_len; col++)
  {
    const glong last_row_for_col = GSL_MIN (ROWS_TO_ROTATE, rhs_len - col - 1);

    for (i = last_row_for_col; i > 0; i--)
    {
      const glong r1_index = col + i - 1;
      const glong r2_index = col + i;

      if (r2_index > last_row_index)
        g_ptr_array_add (self->matrix_rows, _ncm_sbessel_create_row (solver, ++last_row_index));

      _ncm_sbessel_apply_givens (solver, col, r1_index, r2_index, c);
    }
  }

  {
    /* Back substitution to solve R x = c */
    const glong n        = col; /* since col == rhs_len */
    gdouble acc_bc_at_m1 = 0.0;
    gdouble acc_bc_at_p1 = 0.0;
    gdouble *sol_ptr;
    glong row;

    self->solution = ncm_vector_new (n);
    sol_ptr        = ncm_vector_data (self->solution);

    for (row = n - 1; row >= 0; row--)
    {
      NcmSBesselOdeSolverRow *r = g_ptr_array_index (self->matrix_rows, row);
      gdouble sum               = g_array_index (c, gdouble, row);
      glong width               = GSL_MIN (TOTAL_BANDWIDTH, n - row);
      const gdouble diag        = r->data[0] + _ncm_sbessel_bc_row (r, row);
      gdouble sol;
      glong j;

      g_assert_cmpuint (r->col_index, ==, row); /* Banded matrix */

      for (j = 1; j < width; j++)
      {
        const glong col    = r->col_index + j;
        const gdouble a_ij = r->data[j];

        sum -= a_ij * sol_ptr[col];
      }

      sum -= acc_bc_at_m1 * r->bc_at_m1;
      sum -= acc_bc_at_p1 * r->bc_at_p1;

      sol          = sum / diag;
      sol_ptr[row] = sol;

      acc_bc_at_m1 += (row % 2 == 0 ? 1.0 : -1.0) * sol;
      acc_bc_at_p1 += sol;
    }
  }

  g_array_unref (c);

  return ncm_vector_ref (self->solution);
}

/**
 * _ncm_sbessel_ode_solver_fill_operator_matrix:
 * @solver: a #NcmSBesselOdeSolver
 * @nrows: number of rows
 * @ncols: number of columns
 * @data: pointer to matrix data
 * @colmajor: TRUE for column-major layout, FALSE for row-major
 *
 * Helper function to fill operator matrix data in either row-major or column-major format.
 * Extracts matrix entries from operator rows and fills them into the provided data array.
 */
static void
_ncm_sbessel_ode_solver_fill_operator_matrix (NcmSBesselOdeSolver *solver,
                                              gint nrows, gint ncols,
                                              gdouble *data, gboolean colmajor)
{
  /* Fill matrix rows */
  for (gint i = 0; i < nrows; i++)
  {
    NcmSBesselOdeSolverRow *row = _ncm_sbessel_create_row (solver, i);
    glong k;

    /* Handle boundary condition rows with infinite components */
    if (fabs (row->bc_at_m1) > 1.0e-100)
    {
      /* bc_at_m1 contributes (-1)^k at every column k */
      for (k = row->col_index; k < ncols; k++)
      {
        const gdouble sign  = (k % 2 == 0) ? 1.0 : -1.0;
        const gdouble value = row->bc_at_m1 * sign;
        const gint idx      = colmajor ? (k * nrows + i) : (i * ncols + k);

        data[idx] += value;
      }
    }

    if (fabs (row->bc_at_p1) > 1.0e-100)
    {
      /* bc_at_p1 contributes 1.0 at every column k */
      for (k = row->col_index; k < ncols; k++)
      {
        const gdouble value = row->bc_at_p1;
        const gint idx      = colmajor ? (k * nrows + i) : (i * ncols + k);

        data[idx] += value;
      }
    }

    /* Copy data entries from the operator part */
    for (glong j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col < ncols) && (fabs (row->data[j]) > 1.0e-100))
      {
        const gint idx = colmajor ? (col * nrows + i) : (i * ncols + col);

        data[idx] += row->data[j];
      }
    }

    _row_free (row);
  }
}

/**
 * ncm_sbessel_ode_solver_get_operator_matrix:
 * @solver: a #NcmSBesselOdeSolver
 * @nrows: number of rows to extract (including 2 boundary condition rows)
 *
 * Builds and returns a dense matrix representation of the differential operator
 * truncated to @nrows rows in row-major format. This is useful for validation,
 * testing against truth tables, and comparison with standard dense solvers.
 *
 * The matrix includes:
 * - Row 0: boundary condition u(-1) = 0
 * - Row 1: boundary condition u(+1) = 0
 * - Rows 2 to nrows-1: differential operator rows
 *
 * The matrix is square (nrows x nrows) for compatibility with standard solvers.
 *
 * Returns: (transfer full): dense matrix representation of the operator (row-major)
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_operator_matrix (NcmSBesselOdeSolver *solver, gint nrows)
{
  g_assert_cmpint (nrows, >, 2);

  /* Force square matrix */
  const gint ncols = nrows;

  /* Create matrix (row-major by default in NcmMatrix) */
  NcmMatrix *mat = ncm_matrix_new (nrows, ncols);

  ncm_matrix_set_zero (mat);

  /* Fill matrix data */
  _ncm_sbessel_ode_solver_fill_operator_matrix (solver, nrows, ncols,
                                                ncm_matrix_data (mat), FALSE);

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_operator_matrix_colmajor:
 * @solver: a #NcmSBesselOdeSolver
 * @nrows: number of rows to extract (including 2 boundary condition rows)
 *
 * Builds and returns a dense matrix representation of the differential operator
 * truncated to @nrows rows in column-major format (suitable for LAPACK/BLAS routines).
 *
 * The matrix includes:
 * - Row 0: boundary condition u(-1) = 0
 * - Row 1: boundary condition u(+1) = 0
 * - Rows 2 to nrows-1: differential operator rows
 *
 * The matrix is square (nrows x nrows) and stored in column-major order as expected
 * by LAPACK routines like dgesv.
 *
 * Returns: (transfer full): dense matrix representation of the operator (column-major)
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_operator_matrix_colmajor (NcmSBesselOdeSolver *solver, gint nrows)
{
  g_assert_cmpint (nrows, >, 2);

  /* Force square matrix */
  const gint ncols = nrows;

  /* Allocate column-major data */
  gdouble *data_colmajor = g_new0 (gdouble, nrows * ncols);

  /* Fill matrix data in column-major format */
  _ncm_sbessel_ode_solver_fill_operator_matrix (solver, nrows, ncols,
                                                data_colmajor, TRUE);

  /* Create matrix with column-major data (tda = nrows for column-major) */
  NcmMatrix *mat = ncm_matrix_new_data_malloc (data_colmajor, nrows, ncols);

  return mat;
}

/**
 * ncm_sbessel_ode_solver_solve_dense:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @nrows: size of the truncated system to solve
 *
 * Solves the ODE using a dense matrix representation with standard linear algebra.
 * This method is useful for validation and testing, providing a reference solution
 * to compare against the adaptive QR method.
 *
 * The system is truncated to @nrows equations (including 2 boundary conditions),
 * and solved using LAPACK's general linear solver (LU decomposition).
 *
 * Returns: (transfer full): solution vector (Chebyshev coefficients)
 */
NcmVector *
ncm_sbessel_ode_solver_solve_dense (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint nrows)
{
  /* Use column-major matrix for LAPACK */
  NcmMatrix *mat   = ncm_sbessel_ode_solver_get_operator_matrix_colmajor (solver, nrows);
  const gint ncols = ncm_matrix_ncols (mat);

  /* For LAPACK dgesv, we need a square matrix, so use minimum dimension */
  const gint n = GSL_MIN (nrows, ncols);

  /* Prepare RHS vector (truncate/pad if needed) */
  NcmVector *b = ncm_vector_new (n);

  for (gint i = 0; i < n; i++)
  {
    if (i < (gint) ncm_vector_len (rhs))
      ncm_vector_set (b, i, ncm_vector_get (rhs, i));
    else
      ncm_vector_set (b, i, 0.0);
  }

  /* Solve the system using LU decomposition */
  gint *ipiv = g_new (gint, n);
  gint ret   = ncm_lapack_dgesv (n, 1,
                                 ncm_matrix_data (mat), n,
                                 ipiv,
                                 ncm_vector_data (b), n);

  if (ret != 0)
    g_warning ("ncm_sbessel_ode_solver_solve_dense: LAPACK dgesv failed with code %d", ret);

  g_free (ipiv);
  ncm_matrix_free (mat);

  /* Solution is stored in b */
  return b;
}

/**
 * ncm_sbessel_ode_solver_compute_chebyshev_coeffs:
 * @F: (scope call): function to evaluate
 * @a: left endpoint
 * @b: right endpoint
 * @N: number of Chebyshev nodes
 * @user_data: user data for @F
 *
 * Computes Chebyshev coefficients of f(x) on [a,b] using FFTW DCT-I.
 * The function is sampled at Chebyshev nodes $x_k = (a+b)/2 - (b-a)/2\cos(k\pi/(N-1))$
 * and transformed using a Type-I discrete cosine transform.
 *
 * Returns: (transfer full): vector of Chebyshev coefficients
 */
NcmVector *
ncm_sbessel_ode_solver_compute_chebyshev_coeffs (NcmSBesselOdeSolverF F,
                                                 gdouble a, gdouble b,
                                                 guint N,
                                                 gpointer user_data)
{
  guint i;
  NcmVector *coeffs = ncm_vector_new (N);
  gdouble *f_vals   = fftw_malloc (sizeof (gdouble) * N);
  fftw_plan plan_r2r;
  const gdouble mid    = 0.5 * (a + b);
  const gdouble half_h = 0.5 * (b - a);

  /* Sample function at Chebyshev nodes */
  for (i = 0; i < N; i++)
  {
    const gdouble theta = M_PI * i / (N - 1);
    const gdouble x     = mid + half_h * cos (theta);

    f_vals[i] = F (user_data, x);
  }

  /* Create DCT-I plan and execute */
  ncm_cfg_load_fftw_wisdom ("ncm_sbessel_ode_solver");
  ncm_cfg_lock_plan_fftw ();
  plan_r2r = fftw_plan_r2r_1d (N, f_vals, ncm_vector_data (coeffs), FFTW_REDFT00, ncm_cfg_get_fftw_default_flag ());
  ncm_cfg_unlock_plan_fftw ();
  ncm_cfg_save_fftw_wisdom ("ncm_sbessel_ode_solver");

  fftw_execute (plan_r2r);

  /* Normalize coefficients */
  ncm_vector_set (coeffs, 0, ncm_vector_get (coeffs, 0) / ((N - 1.0) * 2.0));
  ncm_vector_set (coeffs, N - 1, ncm_vector_get (coeffs, N - 1) / ((N - 1.0) * 2.0));

  for (i = 1; i < N - 1; i++)
    ncm_vector_set (coeffs, i, ncm_vector_get (coeffs, i) / (N - 1.0));

  fftw_destroy_plan (plan_r2r);
  fftw_free (f_vals);

  return coeffs;
}

/**
 * ncm_sbessel_ode_solver_chebT_to_gegenbauer_lambda1:
 * @c: Chebyshev coefficients vector
 * @g: Gegenbauer $C^{(1)}_n$ coefficients vector (must have same length as @c, pre-allocated by caller)
 *
 * Converts Chebyshev $T_n$ coefficients to Gegenbauer $C^{(1)}_n$ coefficients ($\lambda=1$).
 * Uses the relationship: $T_n = \frac{1}{2}(C^{(1)}_n + C^{(1)}_{n-2})$ for $n \geq 2$.
 */
void
ncm_sbessel_ode_solver_chebT_to_gegenbauer_lambda1 (NcmVector *c, NcmVector *g)
{
  const guint N = ncm_vector_len (c);
  guint i;

  g_assert_cmpuint (ncm_vector_len (g), ==, N);

  if (N == 0)
    return;

  ncm_vector_set_zero (g);

  /* n = 0 case */
  ncm_vector_set (g, 0, ncm_vector_get (g, 0) + ncm_vector_get (c, 0));

  if (N == 1)
    return;

  /* n = 1 case */
  ncm_vector_set (g, 1, ncm_vector_get (g, 1) + ncm_vector_get (c, 1) * 0.5);

  /* n >= 2 */
  for (i = 2; i < N; i++)
  {
    const gdouble ci = ncm_vector_get (c, i);

    ncm_vector_set (g, i, ncm_vector_get (g, i) + 0.5 * ci);
    ncm_vector_set (g, i - 2, ncm_vector_get (g, i - 2) - 0.5 * ci);
  }
}

/**
 * ncm_sbessel_ode_solver_chebT_to_gegenbauer_lambda2:
 * @c: Chebyshev coefficients vector
 * @g: Gegenbauer $C^{(2)}_k$ coefficients vector (must have same length as @c, pre-allocated by caller)
 *
 * Converts Chebyshev $T_n$ coefficients to Gegenbauer $C^{(2)}_k$ coefficients ($\lambda=2$).
 *
 * Uses the projection formula:
 * $$g_k = \frac{1}{2} c_0 \delta_{k,0} + \frac{c_k}{2(k+1)} - \frac{(k+2) c_{k+2}}{(k+1)(k+3)} + \frac{c_{k+4}}{2(k+3)}$$
 * where $f(x) = \sum_n c_n T_n(x)$.
 */
void
ncm_sbessel_ode_solver_chebT_to_gegenbauer_lambda2 (NcmVector *c, NcmVector *g)
{
  const guint N = ncm_vector_len (c);
  guint k;

  g_assert_cmpuint (ncm_vector_len (g), ==, N);

  if (N == 0)
    return;

  /* Zero output vector */
  ncm_vector_set_zero (g);

  /* Apply projection formula for each k */
  for (k = 0; k < N; k++)
  {
    const gdouble kd = (gdouble) k;
    gdouble gk       = 0.0;

    /* Special case: k=0 has additional 1/2 * c[0] contribution */
    if (k == 0)
      gk += 0.5 * ncm_vector_get (c, 0);

    /* First term: c[k] / (2*(k+1)) */
    gk += ncm_vector_get (c, k) / (2.0 * (kd + 1.0));

    /* Second term: -(k+2) * c[k+2] / ((k+1)*(k+3)) */
    if (k + 2 < N)
      gk -= (kd + 2.0) * ncm_vector_get (c, k + 2) / ((kd + 1.0) * (kd + 3.0));

    /* Third term: c[k+4] / (2*(k+3)) */
    if (k + 4 < N)
      gk += ncm_vector_get (c, k + 4) / (2.0 * (kd + 3.0));

    ncm_vector_set (g, k, gk);
  }
}

/**
 * ncm_sbessel_ode_solver_gegenbauer_lambda1_eval:
 * @c: Gegenbauer $C^{(1)}_n$ coefficients vector
 * @x: point to evaluate
 *
 * Evaluates a Gegenbauer $C^{(1)}_n$ expansion at x using Clenshaw recurrence.
 * For $\lambda=1$, $C^{(1)}_n(x) = U_n(x)$ (Chebyshev polynomials of the second kind).
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(1)}_n(x)$
 */
gdouble
ncm_sbessel_ode_solver_gegenbauer_lambda1_eval (NcmVector *c, gdouble x)
{
  const guint N = ncm_vector_len (c);

  if (N == 0)
    return 0.0;

  /* Endpoint handling: C_n^{(1)}(+/-1) = (n+1)*(+/-1)^n */
  if (fabs (x - 1.0) < 1e-15)
  {
    gdouble sum = 0.0;
    guint n;

    for (n = 0; n < N; n++)
      sum += ncm_vector_get (c, n) * (gdouble) (n + 1);

    return sum;
  }

  if (fabs (x + 1.0) < 1e-15)
  {
    gdouble sum = 0.0;
    guint n;

    for (n = 0; n < N; n++)
      sum += ncm_vector_get (c, n) * ((n & 1) ? -(gdouble) (n + 1) : (gdouble) (n + 1));

    return sum;
  }

  {
    /* Stable recurrence for interior x */
    gdouble Cnm1 = 1.0; /* U_0 */
    gdouble sum  = ncm_vector_get (c, 0) * Cnm1;

    if (N == 1)
      return sum;

    gdouble Cn = 2.0 * x; /* U_1 */

    sum += ncm_vector_get (c, 1) * Cn;

    for (guint n = 1; n < N - 1; n++)
    {
      gdouble Cnp1 = 2.0 * x * Cn - Cnm1; /* U_{n+1} */

      sum += ncm_vector_get (c, n + 1) * Cnp1;
      Cnm1 = Cn;
      Cn   = Cnp1;
    }

    return sum;
  }
}

/**
 * ncm_sbessel_ode_solver_gegenbauer_lambda2_eval:
 * @c: Gegenbauer $C^{(2)}_n$ coefficients vector
 * @x: point to evaluate
 *
 * Evaluates a Gegenbauer $C^{(2)}_n$ expansion at x using Clenshaw recurrence.
 * For $\lambda=2$, the recurrence relation is:
 * $(n+1) C^{(2)}_{n+1}(x) = 2(n+2)x C^{(2)}_n(x) - (n+3) C^{(2)}_{n-1}(x)$
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(2)}_n(x)$
 */
gdouble
ncm_sbessel_ode_solver_gegenbauer_lambda2_eval (NcmVector *c, gdouble x)
{
  const guint N = ncm_vector_len (c);

  if (N == 0)
    return 0.0;

  /* Endpoint handling: C_n^{(2)}(+/-1) = binom(n+3,3)*(+/-1)^n = ((n+1)*(n+2)*(n+3)/6)*(+/-1)^n */
  if (fabs (x - 1.0) < 1e-15)
  {
    gdouble sum = 0.0;
    guint n;

    for (n = 0; n < N; n++)
      sum += ncm_vector_get (c, n) * (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

    return sum;
  }

  if (fabs (x + 1.0) < 1e-15)
  {
    gdouble sum = 0.0;
    guint n;

    for (n = 0; n < N; n++)
    {
      const gdouble val = (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

      sum += ncm_vector_get (c, n) * ((n & 1) ? -val : val);
    }

    return sum;
  }

  {
    /* Stable recurrence for interior x */
    gdouble Cnm1 = 1.0; /* C_0^{(2)} = 1 */
    gdouble sum  = ncm_vector_get (c, 0) * Cnm1;

    if (N == 1)
      return sum;

    gdouble Cn = 4.0 * x; /* C_1^{(2)} = 4x */

    sum += ncm_vector_get (c, 1) * Cn;

    for (guint n = 1; n < N - 1; n++)
    {
      /* (n+1) C_{n+1}^{(2)} = 2(n+2)x C_n^{(2)} - (n+3) C_{n-1}^{(2)} */
      gdouble Cnp1 = (2.0 * (gdouble) (n + 2) * x * Cn - (gdouble) (n + 3) * Cnm1) / (gdouble) (n + 1);

      sum += ncm_vector_get (c, n + 1) * Cnp1;
      Cnm1 = Cn;
      Cn   = Cnp1;
    }

    return sum;
  }
}

/**
 * ncm_sbessel_ode_solver_chebyshev_eval:
 * @a: Chebyshev coefficients vector
 * @t: point to evaluate in [-1,1]
 *
 * Evaluates a Chebyshev expansion $f(t) = \sum_{k=0}^{N-1} a_k T_k(t)$
 * using Clenshaw recurrence.
 *
 * Returns: the value of the Chebyshev expansion at t
 */
gdouble
ncm_sbessel_ode_solver_chebyshev_eval (NcmVector *a, gdouble t)
{
  const guint N = ncm_vector_len (a);

  if (N == 0)
    return 0.0;

  if (N == 1)
    return ncm_vector_get (a, 0);

  {
    gdouble b_kplus1 = 0.0;
    gdouble b_kplus2 = 0.0;
    gdouble two_t    = 2.0 * t;
    gint k;

    for (k = (gint) N - 1; k >= 1; k--)
    {
      gdouble b_k = two_t * b_kplus1 - b_kplus2 + ncm_vector_get (a, k);

      b_kplus2 = b_kplus1;
      b_kplus1 = b_k;
    }

    return t * b_kplus1 - b_kplus2 + ncm_vector_get (a, 0);
  }
}

/**
 * ncm_sbessel_ode_solver_chebyshev_deriv:
 * @a: Chebyshev coefficients vector (a_j multiplies T_j)
 * @t: point to evaluate in [-1,1]
 *
 * Evaluates the first derivative of a Chebyshev expansion at $t$.
 *
 * The Chebyshev series is
 * $$ f(t) = \sum_{j=0}^{N-1} a_j T_j(t) , $$
 * and its derivative can be written as
 * $$ f'(t) = \sum_{k=0}^{N-2} b_k T_k(t) . $$
 *
 * The derivative coefficients $b_k$ satisfy
 * $$ b_k = \sum_{j=k+1,k+3,\dots}^{N-1} 2 j a_j , \quad k \ge 1, $$
 * and
 * $$ b_0 = \sum_{j=1,3,5,\dots}^{N-1} j a_j . $$
 *
 * The derivative is evaluated using a fused backward recurrence and
 * the Clenshaw algorithm, without explicitly forming the coefficients $b_k$.
 *
 * Returns: the value of the derivative at $t$
 */
gdouble
ncm_sbessel_ode_solver_chebyshev_deriv (NcmVector *a, gdouble t)
{
  const gint N = ncm_vector_len (a);

  if (N <= 1)
    return 0.0;

  if (N == 2)
    return ncm_vector_get (a, 1);


  if (fabs (t - 1.0) < 1.0e-15)
  {
    /* ---- x = +1 ---- */
    gdouble d1  = 0.0;
    gdouble d2  = 0.0;
    gdouble sum = 0.0;

    for (gint k = N - 2; k >= 1; k--)
    {
      gdouble bk = d2 + 2.0 * (k + 1) * ncm_vector_get (a, k + 1);

      d2   = d1;
      d1   = bk;
      sum += bk;
    }

    /* b0 has the 1/2 factor */
    gdouble b0 = 0.5 * (d2 + 2.0 * ncm_vector_get (a, 1));

    sum += b0;

    return sum;
  }

  if (fabs (t + 1.0) < 1.0e-15)
  {
    /* ---- x = -1 ---- */
    gdouble d1  = 0.0;
    gdouble d2  = 0.0;
    gdouble sum = 0.0;

    /* start with (-1)^(N-2) */
    gdouble sign = ((N - 2) & 1) ? -1.0 : 1.0;

    for (gint k = N - 2; k >= 1; k--)
    {
      gdouble bk = d2 + 2.0 * (k + 1) * ncm_vector_get (a, k + 1);

      d2 = d1;
      d1 = bk;

      sum += sign * bk;
      sign = -sign;
    }

    /* b0 has sign +1 */
    gdouble b0 = 0.5 * (d2 + 2.0 * ncm_vector_get (a, 1));

    sum += b0;

    return sum;
  }

  {
    gdouble c1          = 0.0; /* Clenshaw state k+1 */
    gdouble c2          = 0.0; /* Clenshaw state k+2 */
    gdouble d1          = 0.0; /* recurrence helper */
    gdouble d2          = 0.0;
    const gdouble two_t = 2.0 * t;

    /* k = N-2 ... 1 */
    for (gint k = N - 2; k >= 1; k--)
    {
      /* build b[k] on the fly */
      gdouble bk = d2 + 2.0 * (k + 1) * ncm_vector_get (a, k + 1);

      /* update derivative recurrence */
      d2 = d1;
      d1 = bk;

      /* Clenshaw step */
      gdouble c0 = two_t * c1 - c2 + bk;

      c2 = c1;
      c1 = c0;
    }

    {
      /* k = 0 needs the 1/2 factor */
      gdouble b0 = 0.5 * (d2 + 2.0 * ncm_vector_get (a, 1));

      return t * c1 - c2 + b0;
    }
  }
}

/**
 * ncm_sbessel_ode_solver_get_proj_matrix:
 * @N: size of the matrix
 *
 * Returns the projection (identity) operator matrix that transforms Chebyshev $T_n$
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Returns: (transfer full): the projection operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_proj_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = k;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_proj_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_x_matrix:
 * @N: size of the matrix
 *
 * Returns the multiplication by $x$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x \cdot f(x)$.
 *
 * Returns: (transfer full): the $x$ operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_x_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = (k >= 1) ? (k - 1) : 0;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_x_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_x2_matrix:
 * @N: size of the matrix
 *
 * Returns the multiplication by $x^2$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x^2 \cdot f(x)$.
 *
 * Returns: (transfer full): the $x^2$ operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_x2_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = (k >= 2) ? (k - 2) : 0;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_x2_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_d_matrix:
 * @N: size of the matrix
 *
 * Returns the derivative operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $\frac{df}{dx}$.
 *
 * Returns: (transfer full): the derivative operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_d_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = k + 1;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_d_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_x_d_matrix:
 * @N: size of the matrix
 *
 * Returns the $x \cdot \frac{d}{dx}$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x \cdot \frac{df}{dx}$.
 *
 * Returns: (transfer full): the $x \cdot d$ operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_x_d_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = k;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_x_d_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_d2_matrix:
 * @N: size of the matrix
 *
 * Returns the second derivative operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $\frac{d^2f}{dx^2}$.
 *
 * Returns: (transfer full): the second derivative operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_d2_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = k + 2;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_d2_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_x_d2_matrix:
 * @N: size of the matrix
 *
 * Returns the $x \cdot \frac{d^2}{dx^2}$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x \cdot \frac{d^2f}{dx^2}$.
 *
 * Returns: (transfer full): the $x \cdot d^2$ operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_x_d2_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = k + 1;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_x_d2_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_x2_d2_matrix:
 * @N: size of the matrix
 *
 * Returns the $x^2 \cdot \frac{d^2}{dx^2}$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x^2 \cdot \frac{d^2f}{dx^2}$.
 *
 * Returns: (transfer full): the $x^2 \cdot d^2$ operator matrix
 */
NcmMatrix *
ncm_sbessel_ode_solver_get_x2_d2_matrix (guint N)
{
  NcmMatrix *mat = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col        = k;
    NcmSBesselOdeSolverRow *row = _row_new (0.0, 0.0, left_col);

    _ncm_sbessel_compute_x2_d2_row (row, k, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }

    _row_free (row);
  }

  return mat;
}


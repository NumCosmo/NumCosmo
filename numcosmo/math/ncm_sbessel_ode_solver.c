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
 * Modified spherical Bessel ODE solver using ultraspherical spectral methods.
 *
 * This object acts as a factory and configuration manager for solving modified spherical
 * Bessel ODEs. The actual solving is performed by #NcmSBesselOdeOperator instances created
 * from this solver, which encapsulate the problem-specific parameters (interval [a,b] and
 * angular momentum range [ell_min, ell_max]).
 *
 * # Problem Formulation
 *
 * Solves the modified spherical Bessel ODE for $u(x) = x \cdot j_\ell(x)$ (or $u(x) = x \cdot y_\ell(x)$)
 * mapped to the Chebyshev domain $[-1,1]$. For an interval $[a,b]$, the mapping is
 * $x = m + h\xi$ where $m = (a+b)/2$ and $h = (b-a)/2$.
 *
 * The ODE in the mapped coordinates is:
 * $$
 * \frac{1}{h^2} \frac{d^2 u}{d\xi^2} + \left[(m+h\xi)^2 - \ell(\ell+1)\right] u = f(\xi)
 * $$
 * with Dirichlet boundary conditions $u(-1) = 0$ and $u(1) = 0$.
 *
 * This formulation eliminates the first derivative term compared to the standard
 * spherical Bessel equation for $j_\ell(x)$. By solving for $u = x \cdot j_\ell(x)$ instead
 * of $j_\ell(x)$ directly, the equation becomes simpler and more numerically stable.
 *
 * # Numerical Method
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
 * # Usage Pattern
 *
 * 1. Create a solver with ncm_sbessel_ode_solver_new()
 * 2. Set desired tolerance with ncm_sbessel_ode_solver_set_tolerance()
 * 3. Create one or more operators with ncm_sbessel_ode_solver_create_operator()
 * 4. Use operators to solve problems via ncm_sbessel_ode_operator_solve() or ncm_sbessel_ode_operator_solve_endpoints()
 * 5. Operators can be reused with ncm_sbessel_ode_operator_reset() or destroyed with ncm_sbessel_ode_operator_unref()
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_ode_solver.h"
#include "math/ncm_sf_sbessel.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"
#include "ncm_cfg.h"

#include <math.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Operator bandwidth */
/* Number of rows to rotate lower bandwidth plus number of boundary conditions */
#define NUMBER_OF_BOUNDARY_CONDITIONS 2
#define LOWER_BANDWIDTH 2
#define UPPER_BANDWIDTH 6
#define TOTAL_BANDWIDTH (LOWER_BANDWIDTH + UPPER_BANDWIDTH + 1)
#define ROWS_TO_ROTATE (LOWER_BANDWIDTH + NUMBER_OF_BOUNDARY_CONDITIONS)

#define ALIGNMENT 64 /* 64-byte alignment for cache lines */

#define ROW_CORE_SIZE \
        (sizeof (gdouble) * (TOTAL_BANDWIDTH + 2) + sizeof (glong))

#define PADDED_BANDWIDTH \
        ((TOTAL_BANDWIDTH + 7) & ~7) /* round up to multiple of 8 */
#define OPERATOR_EXPONENT 1

/* Rotation parameter accessors for interleaved cos/sin storage */
#define ROTATION_COS(ptr) (ptr)[0]
#define ROTATION_SIN(ptr) (ptr)[1]

/**
 * NcmSBesselOdeSolverRow:
 *
 * Internal structure for a matrix row.
 * Simple pointer to doubles with boundary condition coefficients.
 */
typedef struct _NcmSBesselOdeSolverRow
{
  gdouble data[PADDED_BANDWIDTH] __attribute__ ((aligned (ALIGNMENT))); /* Pointer to row data */
  gdouble bc_at_m1;                                                     /* Boundary condition row 1 coefficient */
  gdouble bc_at_p1;                                                     /* Boundary condition row 2 coefficient */
  glong col_index;                                                      /* Column index of leftmost element in data */
} NcmSBesselOdeSolverRow __attribute__ ((aligned (ALIGNMENT)));

/**
 * _NcmSBesselOdeOperator:
 *
 * Opaque structure encapsulating operator-specific data for the spherical Bessel ODE.
 * Contains structural parameters, matrix storage, and diagonalization state.
 * Uses reference counting for memory management.
 */
struct _NcmSBesselOdeOperator
{
  /* Reference count */
  gint ref_count;

  /* Structural parameters */
  gdouble a;         /* Left endpoint of interval */
  gdouble b;         /* Right endpoint of interval */
  gdouble half_len;  /* (b-a)/2 - half length of interval */
  gdouble mid_point; /* (a+b)/2 - midpoint of interval */
  gint ell_min;      /* Minimum angular momentum */
  gint ell_max;      /* Maximum angular momentum */
  guint n_ell;       /* Number of ell values (ell_max - ell_min + 1) */
  gdouble tolerance; /* Convergence tolerance for adaptive QR */

  /* Matrix storage */
  NcmSBesselOdeSolverRow *matrix_rows; /* Aligned array of NcmSBesselOdeSolverRow */
  gsize matrix_rows_size;              /* Current size of matrix_rows */
  gsize matrix_rows_capacity;          /* Allocated capacity of matrix_rows */
  GArray *c;                           /* Array of gdouble for right-hand side */

  /* Factorization state */
  glong last_n_cols; /* Number of columns from last diagonalization (0 = no factorization) */

  /* Rotation storage - interleaved cos/sin pairs for cache efficiency
   *
   * For each column col, ROWS_TO_ROTATE Givens rotations are applied to eliminate
   * subdiagonal elements. The rotations are stored in application order:
   *   - rotation 0: eliminates row (col+ROWS_TO_ROTATE) using row (col+ROWS_TO_ROTATE-1)
   *   - rotation 1: eliminates row (col+ROWS_TO_ROTATE-1) using row (col+ROWS_TO_ROTATE-2)
   *   - ...
   *   - rotation (ROWS_TO_ROTATE-1): eliminates row (col+1) using row col
   *
   * Layout for single-ell (n_ell=1):
   *   rotation_params[2*i]   = cos for rotation i
   *   rotation_params[2*i+1] = sin for rotation i
   *   where i = col * ROWS_TO_ROTATE + rotation_within_col
   *
   * Layout for batched (n_ell>1):
   *   rotation_params[2*i]   = cos for rotation i
   *   rotation_params[2*i+1] = sin for rotation i
   *   where i = col * ROWS_TO_ROTATE * n_ell + rotation_within_col * n_ell + ell_idx
   */
  gdouble *rotation_params; /* Interleaved [cos0, sin0, cos1, sin1, ...] */
  gsize rotation_capacity;  /* Allocated capacity for rotation_params */

  /* Aligned arrays for temporary storage */
  gdouble *solution_batched;       /* Aligned array of gdouble for temporary solution storage in endpoint computation */
  gsize solution_batched_capacity; /* Allocated capacity of solution_batched */
  gdouble *acc_bc_at_m1;           /* Aligned array of gdouble for boundary condition accumulators at -1 */
  gdouble *acc_bc_at_p1;           /* Aligned array of gdouble for boundary condition accumulators at +1 */
  gsize acc_bc_capacity;           /* Allocated capacity of acc_bc_at_m1 and acc_bc_at_p1 */
};

typedef struct _NcmSBesselOdeSolverPrivate
{
  gdouble tolerance; /* Convergence tolerance for adaptive QR */

  /* Solution */
  NcmVector *solution; /* Result vector */

  /* Spectral methods context */
  NcmSpectral *spectral;

  /* Spherical Bessel array with cutoff */
  NcmSFSBesselArray *sba;
} NcmSBesselOdeSolverPrivate;

enum
{
  PROP_0,
  PROP_TOLERANCE,
};

struct _NcmSBesselOdeSolver
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSBesselOdeSolver, ncm_sbessel_ode_solver, G_TYPE_OBJECT)
G_DEFINE_BOXED_TYPE (NcmSBesselOdeOperator, ncm_sbessel_ode_operator, ncm_sbessel_ode_operator_ref, ncm_sbessel_ode_operator_unref)

/* Row operations */
static void
_row_reset (NcmSBesselOdeSolverRow *row, gdouble bc_at_m1, gdouble bc_at_p1, glong col_index)
{
  guint i;

  for (i = 0; i < TOTAL_BANDWIDTH; i++)
    row->data[i] = 0.0;

  row->bc_at_m1  = bc_at_m1;
  row->bc_at_p1  = bc_at_p1;
  row->col_index = col_index;
}

/* Operator row computation */
static gdouble _ncm_sbessel_bc_row (NcmSBesselOdeSolverRow *row, glong col_index);

static void
ncm_sbessel_ode_solver_init (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  self->tolerance = 0.0;
  self->solution  = NULL;
  self->sba       = ncm_sf_sbessel_array_new ();
  self->spectral  = ncm_spectral_new ();
}

static void
_ncm_sbessel_ode_solver_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselOdeSolver *solver = NCM_SBESSEL_ODE_SOLVER (object);

  g_return_if_fail (NCM_IS_SBESSEL_ODE_SOLVER (object));

  switch (prop_id)
  {
    case PROP_TOLERANCE:
      ncm_sbessel_ode_solver_set_tolerance (solver, g_value_get_double (value));
      break;
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
    case PROP_TOLERANCE:
      g_value_set_double (value, self->tolerance);
      break;
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
  ncm_spectral_clear (&self->spectral);
  ncm_sf_sbessel_array_clear (&self->sba);

  /* Chain up */
  G_OBJECT_CLASS (ncm_sbessel_ode_solver_parent_class)->dispose (object);
}

static void
_ncm_sbessel_ode_solver_finalize (GObject *object)
{
  /* Chain up */
  G_OBJECT_CLASS (ncm_sbessel_ode_solver_parent_class)->finalize (object);
}

/* Helper functions for aligned memory management in operators */
static void
_ensure_matrix_rows_capacity (NcmSBesselOdeOperator *op, gsize required_size)
{
  if (required_size > op->matrix_rows_capacity)
  {
    NcmSBesselOdeSolverRow *new_rows = NULL;
    gsize new_capacity               = required_size;

    if (new_capacity < 16)
      new_capacity = 16;


    if (posix_memalign ((void **) &new_rows, ALIGNMENT, new_capacity * sizeof (NcmSBesselOdeSolverRow)) != 0)
      g_error ("_ensure_matrix_rows_capacity: failed to allocate aligned memory");

    if (op->matrix_rows != NULL)
    {
      memcpy (new_rows, op->matrix_rows, op->matrix_rows_size * sizeof (NcmSBesselOdeSolverRow));
      free (op->matrix_rows);
    }

    op->matrix_rows          = new_rows;
    op->matrix_rows_capacity = new_capacity;
  }

  op->matrix_rows_size = required_size;
}

static void
_ensure_solution_batched_capacity (NcmSBesselOdeOperator *op, gsize required_size)
{
  if (required_size > op->solution_batched_capacity)
  {
    gdouble *new_array = NULL;
    gsize new_capacity = required_size;

    if (new_capacity < 16)
      new_capacity = 16;

    if (posix_memalign ((void **) &new_array, ALIGNMENT, new_capacity * sizeof (gdouble)) != 0)
      g_error ("_ensure_solution_batched_capacity: failed to allocate aligned memory");

    if (op->solution_batched != NULL)
      free (op->solution_batched);

    op->solution_batched          = new_array;
    op->solution_batched_capacity = new_capacity;
  }
}

static void
_ensure_acc_bc_capacity (NcmSBesselOdeOperator *op, gsize required_size)
{
  if (required_size > op->acc_bc_capacity)
  {
    gdouble *new_m1    = NULL;
    gdouble *new_p1    = NULL;
    gsize new_capacity = required_size;

    if (new_capacity < 16)
      new_capacity = 16;

    if (posix_memalign ((void **) &new_m1, ALIGNMENT, new_capacity * sizeof (gdouble)) != 0)
      g_error ("_ensure_acc_bc_capacity: failed to allocate aligned memory for acc_bc_at_m1");

    if (posix_memalign ((void **) &new_p1, ALIGNMENT, new_capacity * sizeof (gdouble)) != 0)
    {
      free (new_m1);
      g_error ("_ensure_acc_bc_capacity: failed to allocate aligned memory for acc_bc_at_p1");
    }

    if (op->acc_bc_at_m1 != NULL)
      free (op->acc_bc_at_m1);

    if (op->acc_bc_at_p1 != NULL)
      free (op->acc_bc_at_p1);

    op->acc_bc_at_m1    = new_m1;
    op->acc_bc_at_p1    = new_p1;
    op->acc_bc_capacity = new_capacity;
  }
}

static void
_ensure_rotation_capacity (NcmSBesselOdeOperator *op, gsize required_size)
{
  if (required_size > op->rotation_capacity)
  {
    gdouble *new_params = NULL;
    gsize new_capacity  = required_size;

    if (new_capacity < 16)
      new_capacity = 16;

    /* Allocate interleaved cos/sin pairs: need 2x space */
    if (posix_memalign ((void **) &new_params, ALIGNMENT, new_capacity * 2 * sizeof (gdouble)) != 0)
      g_error ("_ensure_rotation_capacity: failed to allocate aligned memory for rotation_params");

    if (op->rotation_params != NULL)
      free (op->rotation_params);

    op->rotation_params   = new_params;
    op->rotation_capacity = new_capacity;
  }
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
   * NcmSBesselOdeSolver:tolerance:
   *
   * Convergence tolerance for adaptive QR decomposition.
   */
  g_object_class_install_property (object_class,
                                   PROP_TOLERANCE,
                                   g_param_spec_double ("tolerance",
                                                        NULL,
                                                        "Convergence tolerance",
                                                        0.0, 1.0, DBL_EPSILON * 1.0e-1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sbessel_ode_solver_new:
 *
 * Creates a new #NcmSBesselOdeSolver.
 *
 * Returns: (transfer full): a new #NcmSBesselOdeSolver
 */
NcmSBesselOdeSolver *
ncm_sbessel_ode_solver_new (void)
{
  NcmSBesselOdeSolver *solver = g_object_new (NCM_TYPE_SBESSEL_ODE_SOLVER,
                                              NULL);

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
 * ncm_sbessel_ode_solver_create_operator:
 * @solver: a #NcmSBesselOdeSolver
 * @a: left endpoint of interval
 * @b: right endpoint of interval
 * @ell_min: minimum angular momentum
 * @ell_max: maximum angular momentum
 *
 * Creates a new #NcmSBesselOdeOperator with the specified structural parameters.
 * The operator encapsulates all problem-specific data including the interval,
 * angular momentum range, matrix storage, and diagonalization state.
 *
 * Returns: (transfer full): a new #NcmSBesselOdeOperator with refcount = 1
 */
NcmSBesselOdeOperator *
ncm_sbessel_ode_solver_create_operator (NcmSBesselOdeSolver *solver, gdouble a, gdouble b, gint ell_min, gint ell_max)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  NcmSBesselOdeOperator *op               = NULL;

  g_assert (solver != NULL);
  g_assert_cmpfloat (a, <, b);
  g_assert_cmpint (ell_min, <=, ell_max);
  g_assert_cmpint (ell_min, >=, 0);

  op = g_new0 (NcmSBesselOdeOperator, 1);

  /* Initialize reference count */
  op->ref_count = 1;

  /* Set structural parameters */
  op->a         = a;
  op->b         = b;
  op->half_len  = (b - a) / 2.0;
  op->mid_point = (a + b) / 2.0;
  op->ell_min   = ell_min;
  op->ell_max   = ell_max;
  op->n_ell     = (guint) (ell_max - ell_min + 1);
  op->tolerance = self->tolerance; /* Copy tolerance from solver */

  /* Initialize single-ell matrix storage */
  op->matrix_rows          = NULL;
  op->matrix_rows_size     = 0;
  op->matrix_rows_capacity = 0;
  op->c                    = g_array_new (FALSE, FALSE, sizeof (gdouble));

  /* Initialize factorization state */
  op->last_n_cols = 0;

  /* Initialize rotation storage */
  op->rotation_params   = NULL;
  op->rotation_capacity = 0;

  /* Initialize aligned arrays for temporary storage */
  op->solution_batched          = NULL;
  op->solution_batched_capacity = 0;
  op->acc_bc_at_m1              = NULL;
  op->acc_bc_at_p1              = NULL;
  op->acc_bc_capacity           = 0;

  return op;
}

/**
 * ncm_sbessel_ode_operator_ref:
 * @op: a #NcmSBesselOdeOperator
 *
 * Increases the reference count of @op by one.
 *
 * Returns: (transfer full): @op
 */
NcmSBesselOdeOperator *
ncm_sbessel_ode_operator_ref (NcmSBesselOdeOperator *op)
{
  g_assert (op != NULL);
  g_assert (op->ref_count > 0);

  g_atomic_int_inc (&op->ref_count);

  return op;
}

/**
 * ncm_sbessel_ode_operator_unref:
 * @op: a #NcmSBesselOdeOperator
 *
 * Decreases the reference count of @op by one.
 * When the reference count reaches zero, frees all allocated memory.
 */
void
ncm_sbessel_ode_operator_unref (NcmSBesselOdeOperator *op)
{
  g_assert (op != NULL);
  g_assert (op->ref_count > 0);

  if (g_atomic_int_dec_and_test (&op->ref_count))
  {
    /* Free matrix storage */
    if (op->matrix_rows != NULL)
      free (op->matrix_rows);

    g_clear_pointer (&op->c, g_array_unref);

    /* Free rotation storage */
    if (op->rotation_params != NULL)
      free (op->rotation_params);

    /* Free aligned arrays */
    if (op->solution_batched != NULL)
      free (op->solution_batched);

    if (op->acc_bc_at_m1 != NULL)
      free (op->acc_bc_at_m1);

    if (op->acc_bc_at_p1 != NULL)
      free (op->acc_bc_at_p1);

    /* Free the operator itself */
    g_free (op);
  }
}

/**
 * ncm_sbessel_ode_operator_clear:
 * @op: a #NcmSBesselOdeOperator
 *
 * Decreases the reference count of *@op by one and sets *@op to NULL.
 */
void
ncm_sbessel_ode_operator_clear (NcmSBesselOdeOperator **op)
{
  if (*op != NULL)
  {
    ncm_sbessel_ode_operator_unref (*op);
    *op = NULL;
  }
}

/**
 * ncm_sbessel_ode_operator_reset:
 * @op: a #NcmSBesselOdeOperator
 * @a: new left endpoint of interval
 * @b: new right endpoint of interval
 * @ell_min: new minimum angular momentum
 * @ell_max: new maximum angular momentum
 *
 * Resets the operator with new structural parameters without modifying the reference count.
 * Clears all factorization state, rotations, and convergence state, but preserves
 * allocated capacities for efficiency. After reset, the operator behaves exactly like
 * a freshly created one.
 */
void
ncm_sbessel_ode_operator_reset (NcmSBesselOdeOperator *op, gdouble a, gdouble b, gint ell_min, gint ell_max)
{
  g_assert (op != NULL);
  g_assert (op->ref_count > 0);
  g_assert_cmpfloat (a, <, b);
  g_assert_cmpint (ell_min, <=, ell_max);
  g_assert_cmpint (ell_min, >=, 0);

  /* Update structural parameters */
  op->a         = a;
  op->b         = b;
  op->half_len  = (b - a) / 2.0;
  op->mid_point = (a + b) / 2.0;
  op->ell_min   = ell_min;
  op->ell_max   = ell_max;
  op->n_ell     = (guint) (ell_max - ell_min + 1);

  /* Clear factorization state (preserve capacity) */
  op->matrix_rows_size = 0;
  op->last_n_cols      = 0;

  /* Clear RHS storage */
  if (op->c != NULL)
    g_array_set_size (op->c, 0);

  /* Note: We don't free allocated buffers to allow reuse */
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

static gdouble
_ncm_sbessel_bc_row (NcmSBesselOdeSolverRow *row, glong col_index)
{
  const gdouble bc_at_m1 = (col_index % 2) == 0 ? row->bc_at_m1 : -row->bc_at_m1;
  const gdouble bc_at_p1 = row->bc_at_p1;

  return (bc_at_m1 + bc_at_p1);
}

/**
 * _ncm_sbessel_create_row_bc_at_m1:
 * @op: a #NcmSBesselOdeOperator (unused in this function)
 * @row: output row
 *
 * Creates boundary condition row for u(-1) = 0.
 * All data entries are zero, bc_at_m1 = 1.0, bc_at_p1 = 0.0.
 */
static void
_ncm_sbessel_create_row_bc_at_m1 (NcmSBesselOdeOperator *op, NcmSBesselOdeSolverRow *row)
{
  G_GNUC_UNUSED NcmSBesselOdeOperator *unused_op = op;

  _row_reset (row, 1.0, 0.0, 0);
}

/**
 * _ncm_sbessel_create_row_bc_at_p1:
 * @op: a #NcmSBesselOdeOperator (unused in this function)
 * @row: output row
 *
 * Creates boundary condition row for u(+1) = 0.
 * All data entries are zero, bc_at_m1 = 0.0, bc_at_p1 = 1.0.
 */
static void
_ncm_sbessel_create_row_bc_at_p1 (NcmSBesselOdeOperator *op, NcmSBesselOdeSolverRow *row)
{
  G_GNUC_UNUSED NcmSBesselOdeOperator *unused_op = op;

  _row_reset (row, 0.0, 1.0, 0);
}

/**
 * _ncm_sbessel_create_row_operator:
 * @op: a #NcmSBesselOdeOperator (for interval parameters)
 * @row: output row
 * @row_index: row index (0-based, corresponding to k in the Gegenbauer basis)
 *
 * Creates differential operator row for the given index k in the Gegenbauer $C^{(2)}_k$ basis.
 * Builds the full ODE operator using the formula:
 * $$\langle C^{(2)}_k, L[f] \rangle$$ where
 * $$L = \frac{(m+hx)^2}{h^2} \frac{d^2}{dx^2} + \frac{2(m+hx)}{h} \frac{d}{dx} +
 * (m+hx)^2 - \ell(\ell+1)$$
 *
 * The operator functions handle all special cases for low k values internally,
 * so this single function works for all k.
 */
static void
_ncm_sbessel_create_row_operator (NcmSBesselOdeOperator *op, NcmSBesselOdeSolverRow *row, glong row_index)
{
  const gdouble m             = op->mid_point;
  const gdouble h             = op->half_len;
  const gdouble h2            = h * h;
  const gdouble m2            = m * m;
  const glong k               = row_index;
  const gdouble ell           = op->ell_min;
  const gdouble llp1          = ell * (ell + 1);
  const glong left_col        = GSL_MAX (0, k - 2); /* Leftmost column index */
  gdouble * restrict row_data = row->data;

  _row_reset (row, 0.0, 0.0, left_col);

  /* Compute offset once for all operators */
  const glong offset = k - row->col_index;

  /* Second derivative term: (m^2/h^2) d^2 + (2m/h) x d^2 + x^2 d^2 */
  ncm_spectral_compute_d2_row (row_data, k, offset, m2 / h2);
  ncm_spectral_compute_x_d2_row (row_data, k, offset, 2.0 * m / h);
  ncm_spectral_compute_x2_d2_row (row_data, k, offset, 1.0);

  /* First derivative term: (1.0 - o_e) * ((2m/h) d + 2 x d) */
#if OPERATOR_EXPONENT != 1
  ncm_spectral_compute_d_row (row_data, offset, 2.0 * (1.0 - OPERATOR_EXPONENT) * m / h);
  ncm_spectral_compute_x_d_row (row_data, k, offset, 2.0 * (1.0 - OPERATOR_EXPONENT));
#endif /* OPERATOR_EXPONENT != 1 */

  /* Identity term: (m^2 - ell(ell+1)) I + 2m h x + h^2 x^2 + o_e * (o_e - 1) */
  ncm_spectral_compute_proj_row (row_data, k, offset, m2 - llp1 + (OPERATOR_EXPONENT - 1.0) * OPERATOR_EXPONENT);
  ncm_spectral_compute_x_row (row_data, k, offset, 2.0 * m * h);
  ncm_spectral_compute_x2_row (row_data, k, offset, h2);
}

/**
 * _ncm_sbessel_create_row_operator_batched:
 * @op: a #NcmSBesselOdeOperator (for interval parameters)
 * @row: array of rows, one for each ell value
 * @row_index: row index (0-based, corresponding to k in the Gegenbauer basis)
 * @ell_min: minimum ell value
 * @n_ell: number of ell values to process
 *
 * Creates differential operator rows for multiple ell values at once using an optimized
 * incremental algorithm. The key insight is that operator rows for consecutive ell values
 * differ only in the ell(ell+1) term:
 *
 * - All derivative operators (d², x d², x² d², d, x d) are ell-independent
 * - Identity terms (x, x²) are ell-independent
 * - Only the projection term m² - ell(ell+1) depends on ell
 *
 * The optimization uses the recurrence: (ell+1)(ell+2) - ell(ell+1) = 2(ell+1), allowing
 * incremental updates by adding -2l for each successive ell value. This reduces
 * O(n_ell × operations) to O(operations + n_ell) complexity.
 *
 * Algorithm:
 * 1. Compute all ell-independent terms once into row[0]
 * 2. Add ell-dependent term (m² - ell_min(ell_min+1)) to row[0]
 * 3. For each subsequent row: memcpy row[0] and add incremental correction -2l
 */
static void
_ncm_sbessel_create_row_operator_batched (NcmSBesselOdeOperator *op, NcmSBesselOdeSolverRow *row, glong row_index, gint ell_min, guint n_ell)
{
  const gdouble m      = op->mid_point;
  const gdouble h      = op->half_len;
  const gdouble h2     = h * h;
  const gdouble m2     = m * m;
  const glong k        = row_index;
  const glong left_col = GSL_MAX (0, k - 2); /* Leftmost column index */
  gdouble ell          = (gdouble) ell_min;
  guint i;

  /* Initialize first row with all ell-independent terms */
  _row_reset (&row[0], 0.0, 0.0, left_col);

  /* Compute offset once */
  const glong offset = k - row[0].col_index;

  /* Second derivative term: (m^2/h^2) d^2 + (2m/h) x d^2 + x^2 d^2 - ell-independent */
  ncm_spectral_compute_d2_row (row[0].data, k, offset, m2 / h2);
  ncm_spectral_compute_x_d2_row (row[0].data, k, offset, 2.0 * m / h);
  ncm_spectral_compute_x2_d2_row (row[0].data, k, offset, 1.0);

  /* First derivative term: (2m/h) d + 2 x d - ell-independent */
#if OPERATOR_EXPONENT != 1
  ncm_spectral_compute_d_row (row[0].data, offset, 2.0 * (1.0 - OPERATOR_EXPONENT) * m / h);
  ncm_spectral_compute_x_d_row (row[0].data, k, offset, 2.0 * (1.0 - OPERATOR_EXPONENT));
#endif /* OPERATOR_EXPONENT != 1 */

  /* Identity term: 2m h x + h^2 x^2 - ell-independent part */
  ncm_spectral_compute_x_row (row[0].data, k, offset, 2.0 * m * h);
  ncm_spectral_compute_x2_row (row[0].data, k, offset, h2);

  /* Copy template from row[0] to all other rows, then add ell-dependent corrections in one pass */
  for (i = 1; i < n_ell; i++)
  {
    row[i] = row[0]; /* Struct copy is faster than memcpy for small structs and better for cache */
  }

  {
    const gdouble kd = (gdouble) k;

    /* General formula: three entries with rational function coefficients */
    const gdouble value_0 = 1.0 / (2.0 * (kd + 1.0));
    const gdouble value_1 = -1.0 * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
    const gdouble value_2 = 1.0 / (2.0 * (kd + 3.0));
    gdouble coeff         = m2 - ell * (ell + 1.0) + (OPERATOR_EXPONENT - 1.0) * OPERATOR_EXPONENT;

    for (i = 0; i < n_ell; i++)
    {
      gdouble * restrict row_data = row[i].data;

      row_data[offset]     += coeff * value_0;
      row_data[offset + 2] += coeff * value_1;
      row_data[offset + 4] += coeff * value_2;

      /* Special case: k=0 has an additional contribution */
      if (k == 0)
        row_data[offset] += coeff * 0.5;

      ell   += 1.0;
      coeff -= 2.0 * ell; /* Incremental update for next ell value */
    }
  }
}

/**
 * _ncm_sbessel_create_row:
 * @op: a #NcmSBesselOdeOperator
 * @row: output row
 * @row_index: row index (first two rows are boundary conditions, rest are operators)
 *
 * Creates the appropriate row based on the row index.
 * Rows 0-1 are boundary conditions, rows >= 2 are differential operators.
 */
static void
_ncm_sbessel_create_row (NcmSBesselOdeOperator *op, NcmSBesselOdeSolverRow *row, glong row_index)
{
  switch (row_index)
  {
    case 0:
      _ncm_sbessel_create_row_bc_at_m1 (op, row);
      break;

    case 1:
      _ncm_sbessel_create_row_bc_at_p1 (op, row);
      break;

    default:
      _ncm_sbessel_create_row_operator (op, row, row_index - 2);
      break;
  }
}

/**
 * _ncm_sbessel_create_row_batched:
 * @op: a #NcmSBesselOdeOperator
 * @row: array of rows, one for each ell value
 * @row_index: row index (first two rows are boundary conditions, rest are operators)
 * @ell_min: minimum ell value
 * @n_ell: number of ell values to process
 *
 * Creates the appropriate rows for multiple ell values at once.
 * The i-th element row[i] will contain the row for ell = ell_min + i.
 * Rows 0-1 are boundary conditions (same for all ell), rows >= 2 are differential operators.
 */
static void
_ncm_sbessel_create_row_batched (NcmSBesselOdeOperator *op, NcmSBesselOdeSolverRow *row, glong row_index, gint ell_min, guint n_ell)
{
  guint i;

  switch (row_index)
  {
    case 0:

      /* Boundary condition at -1: same for all ell values */
      for (i = 0; i < n_ell; i++)
        _ncm_sbessel_create_row_bc_at_m1 (op, &row[i]);

      break;

    case 1:

      /* Boundary condition at +1: same for all ell values */
      for (i = 0; i < n_ell; i++)
        _ncm_sbessel_create_row_bc_at_p1 (op, &row[i]);

      break;

    default:
      /* Differential operator: optimized batch computation */
      _ncm_sbessel_create_row_operator_batched (op, row, row_index - 2, ell_min, n_ell);
      break;
  }
}

__attribute__ ((optimize ("no-math-errno", "no-trapping-math")))

static inline double
_compute_inv_hypot (double a, double b)
{
  return 1.0 / sqrt (a * a + b * b);
}

/**
 * _ncm_sbessel_apply_givens:
 * @pivot_col: column index of the pivot
 * @r1: pointer to first row (pivot row)
 * @r2: pointer to second row (row to eliminate)
 * @c1: pointer to right-hand side element for row1
 * @c2: pointer to right-hand side element for row2
 * @rot_ptr: pointer to rotation parameter pair [cos, sin]
 *
 * Applies Givens rotation to eliminate the entry at (row2, col=pivot_col).
 * The rotation transforms the rows as:
 *   new_row1 = c*row1 + s*row2
 *   new_row2 = -s*row1 + c*row2
 * where c and s are chosen to zero out element (row2, pivot_col).
 * The rotation parameters are stored at rot_ptr[0] (cos) and rot_ptr[1] (sin).
 */
static inline __attribute__ ((hot)) void

_ncm_sbessel_apply_givens (glong pivot_col, NcmSBesselOdeSolverRow * restrict r1, NcmSBesselOdeSolverRow * restrict r2, gdouble * restrict c1, gdouble * restrict c2, gdouble * restrict rot_ptr)
{
  const gdouble a_val    = r1->data[0] + _ncm_sbessel_bc_row (r1, pivot_col);
  const gdouble b_val    = r2->data[0] + _ncm_sbessel_bc_row (r2, pivot_col);
  const gdouble inv_norm = _compute_inv_hypot (a_val, b_val);

  if (__builtin_expect ((inv_norm > 1.0e100), 0))
  {
    /* Entry already zero - just shift r2 without rotation */
    /* Store identity rotation */
    ROTATION_COS (rot_ptr) = 1.0;
    ROTATION_SIN (rot_ptr) = 0.0;

    /* Manually unrolled for TOTAL_BANDWIDTH = 9 */
    r2->data[0] = r2->data[1];
    r2->data[1] = r2->data[2];
    r2->data[2] = r2->data[3];
    r2->data[3] = r2->data[4];
    r2->data[4] = r2->data[5];
    r2->data[5] = r2->data[6];
    r2->data[6] = r2->data[7];
    r2->data[7] = r2->data[8];
    r2->data[8] = 0.0;

    r2->col_index++;

    return;
  }

  {
    /* Compute Givens rotation coefficients */
    #define FMA(a, b, c) ((a) * (b) + (c))
    /* #define FMA(a, b, c) fma (a, b, c) */
    const gdouble cos_theta    = a_val * inv_norm;
    const gdouble sin_theta    = b_val * inv_norm;
    const gdouble rhs1         = *c1;
    const gdouble rhs2         = *c2;
    const gdouble bc1_m1       = r1->bc_at_m1;
    const gdouble bc2_m1       = r2->bc_at_m1;
    const gdouble bc1_p1       = r1->bc_at_p1;
    const gdouble bc2_p1       = r2->bc_at_p1;
    gdouble * restrict r1_data = r1->data;
    gdouble * restrict r2_data = r2->data;
    const gdouble e1_0         = r1_data[0];
    const gdouble e1_1         = r1_data[1];
    const gdouble e1_2         = r1_data[2];
    const gdouble e1_3         = r1_data[3];
    const gdouble e1_4         = r1_data[4];
    const gdouble e1_5         = r1_data[5];
    const gdouble e1_6         = r1_data[6];
    const gdouble e1_7         = r1_data[7];
    const gdouble e1_8         = r1_data[8];
    const gdouble e2_0         = r2_data[0];
    const gdouble e2_1         = r2_data[1];
    const gdouble e2_2         = r2_data[2];
    const gdouble e2_3         = r2_data[3];
    const gdouble e2_4         = r2_data[4];
    const gdouble e2_5         = r2_data[5];
    const gdouble e2_6         = r2_data[6];
    const gdouble e2_7         = r2_data[7];
    const gdouble e2_8         = r2_data[8];

    /* Apply rotation to all matrix entries - manually unrolled for TOTAL_BANDWIDTH = 9 */

    *c1          = FMA (cos_theta, rhs1, sin_theta * rhs2);
    r1->bc_at_m1 = FMA (cos_theta, bc1_m1, sin_theta * bc2_m1);
    r1->bc_at_p1 = FMA (cos_theta, bc1_p1, sin_theta * bc2_p1);

    r1_data[0] = FMA (cos_theta, e1_0, sin_theta * e2_0);
    r1_data[1] = FMA (cos_theta, e1_1, sin_theta * e2_1);
    r1_data[2] = FMA (cos_theta, e1_2, sin_theta * e2_2);
    r1_data[3] = FMA (cos_theta, e1_3, sin_theta * e2_3);
    r1_data[4] = FMA (cos_theta, e1_4, sin_theta * e2_4);
    r1_data[5] = FMA (cos_theta, e1_5, sin_theta * e2_5);
    r1_data[6] = FMA (cos_theta, e1_6, sin_theta * e2_6);
    r1_data[7] = FMA (cos_theta, e1_7, sin_theta * e2_7);
    r1_data[8] = FMA (cos_theta, e1_8, sin_theta * e2_8);

    *c2          = FMA (-sin_theta, rhs1, cos_theta * rhs2);
    r2->bc_at_m1 = FMA (-sin_theta, bc1_m1, cos_theta * bc2_m1);
    r2->bc_at_p1 = FMA (-sin_theta, bc1_p1, cos_theta * bc2_p1);

    r2_data[0] = FMA (-sin_theta, e1_1, cos_theta * e2_1);
    r2_data[1] = FMA (-sin_theta, e1_2, cos_theta * e2_2);
    r2_data[2] = FMA (-sin_theta, e1_3, cos_theta * e2_3);
    r2_data[3] = FMA (-sin_theta, e1_4, cos_theta * e2_4);
    r2_data[4] = FMA (-sin_theta, e1_5, cos_theta * e2_5);
    r2_data[5] = FMA (-sin_theta, e1_6, cos_theta * e2_6);
    r2_data[6] = FMA (-sin_theta, e1_7, cos_theta * e2_7);
    r2_data[7] = FMA (-sin_theta, e1_8, cos_theta * e2_8);
    r2_data[8] = 0.0; /* Last entry shifted out is now zero */

    /* Store rotation parameters */
    ROTATION_COS (rot_ptr) = cos_theta;
    ROTATION_SIN (rot_ptr) = sin_theta;

    r2->col_index++; /* Shift row2's column index right by 1 */
  }
}

/**
 * _ncm_sbessel_check_convergence:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @col: current column being processed
 * @max_c_A: (inout): maximum coefficient-to-diagonal ratio
 * @quiet_cols: (inout): count of consecutive columns with small coefficients
 * @solution_order: (inout): current solution order
 * @max_solution_order: maximum allowed solution order
 *
 * Checks convergence and handles adaptive solution order increase.
 *
 * Returns: TRUE if converged (early stop condition met), FALSE otherwise
 */
static inline gboolean
_ncm_sbessel_check_convergence (NcmSBesselOdeOperator *op, glong col,
                                gdouble *max_c_A, guint *quiet_cols,
                                guint *solution_order, guint max_solution_order);

static inline gboolean
_ncm_sbessel_check_convergence_batched (NcmSBesselOdeOperator *op, glong col, guint n_ell,
                                        gdouble *max_c_A, guint *quiet_cols,
                                        guint *solution_order,
                                        guint max_solution_order);

static inline gboolean
_ncm_sbessel_check_convergence (NcmSBesselOdeOperator *op, glong col,
                                gdouble *max_c_A, guint *quiet_cols,
                                guint *solution_order, guint max_solution_order)
{
  NcmSBesselOdeSolverRow *row = &op->matrix_rows[col];
  const double diag           = row->data[0] + _ncm_sbessel_bc_row (row, col);
  const double c_col          = g_array_index (op->c, gdouble, col);
  const double Acol           = fabs (c_col / diag);

  if (Acol > *max_c_A)
  {
    *max_c_A    = Acol;
    *quiet_cols = 0;
  }
  else
  {
    /* Acol is small relative to established scale */
    if (Acol < *max_c_A * op->tolerance)
      (*quiet_cols)++;
    else
      *quiet_cols = 0;
  }

  if (*quiet_cols >= ROWS_TO_ROTATE + 1)
    return TRUE;  /* SAFE early stop */

  /* If we've reached the end without convergence, increase solution_order and continue */
  if ((col == *solution_order - ROWS_TO_ROTATE - 1) && (*solution_order < max_solution_order))
  {
    *solution_order *= 2;
    _ensure_matrix_rows_capacity (op, *solution_order);

    if (*solution_order + ROWS_TO_ROTATE > op->c->len)
      g_array_set_size (op->c, *solution_order + ROWS_TO_ROTATE);
  }

  return FALSE;
}

/**
 * _ncm_sbessel_apply_stored_rotation_to_rhs:
 * @rot_ptr: pointer to rotation parameter pair [cos, sin]
 * @c1: pointer to RHS element for row1
 * @c2: pointer to RHS element for row2
 *
 * Applies a stored Givens rotation to RHS elements only (no matrix operations).
 * This is much faster than the full Givens rotation as it only performs 4 FLOPs.
 * The rotation transforms:
 *   new_c1 = cos*c1 + sin*c2
 *   new_c2 = -sin*c1 + cos*c2
 */
static inline void
_ncm_sbessel_apply_stored_rotation_to_rhs (gdouble * restrict rot_ptr, gdouble * restrict c1, gdouble * restrict c2)
{
  const gdouble cos_theta = ROTATION_COS (rot_ptr);
  const gdouble sin_theta = ROTATION_SIN (rot_ptr);
  const gdouble rhs1      = *c1;
  const gdouble rhs2      = *c2;

  *c1 = cos_theta * rhs1 + sin_theta * rhs2;
  *c2 = -sin_theta * rhs1 + cos_theta * rhs2;
}

/**
 * _ncm_sbessel_apply_all_stored_rotations:
 * @op: a #NcmSBesselOdeOperator
 * @rhs: right-hand side vector
 * @max_col: maximum column to process (from last_n_cols)
 *
 * Applies all stored rotations to a new RHS, checking convergence after each column.
 * Returns the column where convergence occurred, or -1 if extension is needed.
 *
 * Returns: column where converged, or -1 if needs extension
 */
static glong
_ncm_sbessel_apply_all_stored_rotations (NcmSBesselOdeOperator *op, GArray *rhs, glong max_col)
{
  const gdouble *rhs_data = (gdouble *) rhs->data;
  const glong rhs_limit   = rhs->len - ROWS_TO_ROTATE;
  const glong first_limit = (max_col < rhs_limit) ? max_col : rhs_limit;
  gdouble max_c_A         = 0.0;
  guint quiet_cols        = 0;
  glong col;

  g_array_set_size (op->c, max_col + ROWS_TO_ROTATE);

  for (col = 0; col < ROWS_TO_ROTATE; col++)
    g_array_index (op->c, gdouble, col) = rhs_data[col];

  for (col = 0; col < first_limit; col++)
  {
    const glong new_row = col + ROWS_TO_ROTATE;
    glong i;

    g_array_index (op->c, gdouble, new_row) = rhs_data[new_row];

    for (i = ROWS_TO_ROTATE; i > 0; i--)
    {
      const glong rot_idx = col * ROWS_TO_ROTATE + (ROWS_TO_ROTATE - i);
      gdouble *rot_ptr    = &op->rotation_params[2 * rot_idx];
      gdouble *c1         = &g_array_index (op->c, gdouble, col + i - 1);
      gdouble *c2         = &g_array_index (op->c, gdouble, col + i);

      _ncm_sbessel_apply_stored_rotation_to_rhs (rot_ptr, c1, c2);
    }

    if (_ncm_sbessel_check_convergence (op, col, &max_c_A, &quiet_cols, NULL, 0))
      return col;
  }

  for (col = first_limit; col < max_col; col++)
  {
    const glong new_row = col + ROWS_TO_ROTATE;
    glong i;

    g_array_index (op->c, gdouble, new_row) = 0.0;

    for (i = ROWS_TO_ROTATE; i > 0; i--)
    {
      const glong rot_idx = col * ROWS_TO_ROTATE + (ROWS_TO_ROTATE - i);
      gdouble *rot_ptr    = &op->rotation_params[2 * rot_idx];
      gdouble *c1         = &g_array_index (op->c, gdouble, col + i - 1);
      gdouble *c2         = &g_array_index (op->c, gdouble, col + i);

      _ncm_sbessel_apply_stored_rotation_to_rhs (rot_ptr, c1, c2);
    }

    if (_ncm_sbessel_check_convergence (op, col, &max_c_A, &quiet_cols, NULL, 0))
      return col;
  }

  return -1;
}

/**
 * _ncm_sbessel_apply_all_stored_rotations_batched:
 * @op: a #NcmSBesselOdeOperator
 * @rhs: right-hand side vector
 * @max_col: maximum column to process (from last_n_cols)
 * @n_ell: number of ell values
 *
 * Batched version: applies all stored rotations to new RHS for all ell values.
 *
 * Returns: column where converged, or -1 if needs extension
 */
static glong
_ncm_sbessel_apply_all_stored_rotations_batched (NcmSBesselOdeOperator *op, GArray *rhs, glong max_col, guint n_ell)
{
  const gdouble *rhs_data = (gdouble *) rhs->data;
  const glong rhs_limit   = rhs->len - ROWS_TO_ROTATE;
  const glong first_limit = (max_col < rhs_limit) ? max_col : rhs_limit;
  gdouble max_c_A         = 0.0;
  guint quiet_cols        = 0;
  glong col;
  guint l_idx;

  g_array_set_size (op->c, (max_col + ROWS_TO_ROTATE) * n_ell);

  for (col = 0; col < ROWS_TO_ROTATE; col++)
  {
    #pragma omp simd
    for (l_idx = 0; l_idx < n_ell; l_idx++)
      g_array_index (op->c, gdouble, col * n_ell + l_idx) = rhs_data[col];
  }

  for (col = 0; col < first_limit; col++)
  {
    const glong new_row = col + ROWS_TO_ROTATE;
    glong i;

    #pragma omp simd
    for (l_idx = 0; l_idx < n_ell; l_idx++)
      g_array_index (op->c, gdouble, new_row * n_ell + l_idx) = rhs_data[new_row];

    for (i = ROWS_TO_ROTATE; i > 0; i--)
    {
      #pragma omp simd
      for (l_idx = 0; l_idx < n_ell; l_idx++)
      {
        const glong rot_idx = col * ROWS_TO_ROTATE * n_ell + (ROWS_TO_ROTATE - i) * n_ell + l_idx;
        gdouble *rot_ptr    = &op->rotation_params[2 * rot_idx];
        gdouble *c1         = &g_array_index (op->c, gdouble, (col + i - 1) * n_ell + l_idx);
        gdouble *c2         = &g_array_index (op->c, gdouble, (col + i) * n_ell + l_idx);

        _ncm_sbessel_apply_stored_rotation_to_rhs (rot_ptr, c1, c2);
      }
    }

    if (_ncm_sbessel_check_convergence_batched (op, col, n_ell, &max_c_A, &quiet_cols, NULL, 0))
      return col;
  }

  for (col = first_limit; col < max_col; col++)
  {
    const glong new_row = col + ROWS_TO_ROTATE;
    glong i;

    #pragma omp simd
    for (l_idx = 0; l_idx < n_ell; l_idx++)
      g_array_index (op->c, gdouble, new_row * n_ell + l_idx) = 0.0;

    for (i = ROWS_TO_ROTATE; i > 0; i--)
    {
      #pragma omp simd
      for (l_idx = 0; l_idx < n_ell; l_idx++)
      {
        const glong rot_idx = col * ROWS_TO_ROTATE * n_ell + (ROWS_TO_ROTATE - i) * n_ell + l_idx;
        gdouble *rot_ptr    = &op->rotation_params[2 * rot_idx];
        gdouble *c1         = &g_array_index (op->c, gdouble, (col + i - 1) * n_ell + l_idx);
        gdouble *c2         = &g_array_index (op->c, gdouble, (col + i) * n_ell + l_idx);

        _ncm_sbessel_apply_stored_rotation_to_rhs (rot_ptr, c1, c2);
      }
    }

    if (_ncm_sbessel_check_convergence_batched (op, col, n_ell, &max_c_A, &quiet_cols, NULL, 0))
      return col;
  }

  return -1;
}

/**
 * _ncm_sbessel_ode_solver_setup_initial_rows:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @solution_order: initial solution order
 *
 * Sets up the first ROWS_TO_ROTATE rows (boundary conditions) and initializes
 * the RHS vector. This function is called only when starting a fresh diagonalization.
 */
static void
_ncm_sbessel_ode_solver_setup_initial_rows (NcmSBesselOdeOperator *op, GArray *rhs, guint solution_order)
{
  const gdouble *rhs_data = (gdouble *) rhs->data;
  glong i;

  _ensure_matrix_rows_capacity (op, solution_order);
  g_array_set_size (op->c, solution_order + ROWS_TO_ROTATE);

  /* Add the first ROWS_TO_ROTATE rows (boundary conditions) */
  for (i = 0; i < ROWS_TO_ROTATE; i++)
  {
    NcmSBesselOdeSolverRow *row = &op->matrix_rows[i];

    g_array_index (op->c, gdouble, i) = rhs_data[i];
    _ncm_sbessel_create_row (op, row, i);
  }
}

/**
 * _ncm_sbessel_ode_solver_diagonalize:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 *
 * Diagonalizes the operator using adaptive QR decomposition.
 * This function applies Givens rotations to transform the system into upper triangular form
 * and applies the same rotations to the RHS vector. The transformed RHS is stored in
 * op->c and the upper triangular matrix is stored in op->matrix_rows.
 *
 * Returns: the effective number of columns used (may be less than rhs_len due to convergence)
 */
static glong
_ncm_sbessel_ode_solver_diagonalize (NcmSBesselOdeOperator *op, GArray *rhs)
{
  const guint rhs_len            = rhs->len;
  guint solution_order           = rhs_len * 2;
  const guint max_solution_order = 1 << 24;
  const gdouble *rhs_data        = (gdouble *) rhs->data;
  gboolean converged             = FALSE;
  const guint first_loop_len     = rhs_len - ROWS_TO_ROTATE;
  glong col                      = 0;
  gdouble max_c_A                = 0.0;
  guint quiet_cols               = 0;
  glong i;

  g_assert_cmpuint (rhs_len, >, ROWS_TO_ROTATE);

  /* Setup initial rows and RHS */
  _ncm_sbessel_ode_solver_setup_initial_rows (op, rhs, solution_order);

  for (col = 0; col < first_loop_len; col++)
  {
    glong new_row               = col + ROWS_TO_ROTATE;
    NcmSBesselOdeSolverRow *row = &op->matrix_rows[new_row];

    /* Ensure rotation storage capacity */
    _ensure_rotation_capacity (op, (col + 1) * ROWS_TO_ROTATE);

    _ncm_sbessel_create_row (op, row, new_row);
    g_array_index (op->c, gdouble, new_row) = rhs_data[new_row];

    for (i = ROWS_TO_ROTATE; i > 0; i--)
    {
      const glong r1_index       = col + i - 1;
      const glong r2_index       = col + i;
      const glong rot_idx        = col * ROWS_TO_ROTATE + (ROWS_TO_ROTATE - i);
      NcmSBesselOdeSolverRow *r1 = &op->matrix_rows[r1_index];
      NcmSBesselOdeSolverRow *r2 = &op->matrix_rows[r2_index];
      gdouble *c1                = &g_array_index (op->c, gdouble, r1_index);
      gdouble *c2                = &g_array_index (op->c, gdouble, r2_index);
      gdouble *rot_ptr           = &op->rotation_params[2 * rot_idx];

      _ncm_sbessel_apply_givens (col, r1, r2, c1, c2, rot_ptr);
    }

    /* Check convergence */
    if (_ncm_sbessel_check_convergence (op, col, &max_c_A, &quiet_cols,
                                        &solution_order, max_solution_order))
    {
      converged = TRUE;
      break; /* SAFE early stop */
    }
  }

  if (!converged)
  {
    for ( ; col < solution_order; col++)
    {
      glong new_row               = col + ROWS_TO_ROTATE;
      NcmSBesselOdeSolverRow *row = &op->matrix_rows[new_row];

      /* Ensure rotation storage capacity */
      _ensure_rotation_capacity (op, (col + 1) * ROWS_TO_ROTATE);

      _ncm_sbessel_create_row (op, row, new_row);
      g_array_index (op->c, gdouble, new_row) = 0.0;

      for (i = ROWS_TO_ROTATE; i > 0; i--)
      {
        const glong r1_index       = col + i - 1;
        const glong r2_index       = col + i;
        const glong rot_idx        = col * ROWS_TO_ROTATE + (ROWS_TO_ROTATE - i);
        NcmSBesselOdeSolverRow *r1 = &op->matrix_rows[r1_index];
        NcmSBesselOdeSolverRow *r2 = &op->matrix_rows[r2_index];
        gdouble *c1                = &g_array_index (op->c, gdouble, r1_index);
        gdouble *c2                = &g_array_index (op->c, gdouble, r2_index);
        gdouble *rot_ptr           = &op->rotation_params[2 * rot_idx];

        _ncm_sbessel_apply_givens (col, r1, r2, c1, c2, rot_ptr);
      }

      /* Check convergence */
      if (_ncm_sbessel_check_convergence (op, col, &max_c_A, &quiet_cols,
                                          &solution_order, max_solution_order))
        break;  /* SAFE early stop */
    }
  }

  /* Warn if we exhausted max_solution_order without convergence */
  if ((solution_order >= max_solution_order) && (quiet_cols < ROWS_TO_ROTATE + 1))
    g_warning ("_ncm_sbessel_ode_solver_diagonalize: "
               "reached max_solution_order=%u without convergence (quiet_cols=%u, needed %d). "
               "Results may be inaccurate.",
               max_solution_order, quiet_cols, ROWS_TO_ROTATE + 1);

  return col;
}

/**
 * _ncm_sbessel_ode_solver_build_solution:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @n_cols: number of columns in the solution (from diagonalization)
 * @solution: (out): output array for solution coefficients
 *
 * Builds the full Chebyshev coefficient solution by back-substitution on the
 * upper triangular system. Assumes _ncm_sbessel_ode_solver_diagonalize
 * has been called first.
 */
static void
_ncm_sbessel_ode_solver_build_solution (NcmSBesselOdeOperator *op, glong n_cols, GArray *solution)
{
  gdouble acc_bc_at_m1 = 0.0;
  gdouble acc_bc_at_p1 = 0.0;
  gdouble * restrict sol_ptr;
  glong row;

  g_array_set_size (solution, n_cols + TOTAL_BANDWIDTH);
  sol_ptr = (gdouble *) solution->data;
  memset (sol_ptr, 0, sizeof (gdouble) * (n_cols + TOTAL_BANDWIDTH));

  for (row = n_cols - 1; row >= 0; row--)
  {
    NcmSBesselOdeSolverRow *r = &op->matrix_rows[row];
    gdouble sum               = g_array_index (op->c, gdouble, row);
    const gdouble diag        = r->data[0] + _ncm_sbessel_bc_row (r, row);
    gdouble sol;
    glong j;

    g_assert_cmpuint (r->col_index, ==, row); /* Banded matrix */

    for (j = 1; j < TOTAL_BANDWIDTH; j++)
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

  g_array_set_size (solution, n_cols);
}

/**
 * _ncm_sbessel_ode_solver_compute_endpoints:
 * @solver: a #NcmSBesselOdeSolver
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @n_cols: number of columns in the solution (from diagonalization)
 * @endpoints: (out callee-allocates) (element-type gdouble): output array for endpoint derivatives and error estimate
 *
 * Computes endpoint derivatives u'(a) and u'(b) and error estimate directly from the
 * diagonalized system without building the full solution vector. This is much more
 * efficient when only endpoint information is needed, as it computes coefficients
 * on-the-fly during back-substitution and accumulates their contributions to the
 * derivatives without storing the full coefficient array.
 *
 * Assumes _ncm_sbessel_ode_solver_diagonalize has been called first.
 */
static void
_ncm_sbessel_ode_solver_compute_endpoints (NcmSBesselOdeOperator *op, glong n_cols, GArray *endpoints)
{
  const gdouble h        = op->half_len;
  gdouble acc_bc_at_m1   = 0.0;
  gdouble acc_bc_at_p1   = 0.0;
  gdouble deriv_at_m1    = 0.0; /* u'(-1) accumulator */
  gdouble deriv_at_p1    = 0.0; /* u'(+1) accumulator */
  gdouble error_estimate = 0.0; /* error accumulator */
  gdouble *sol_buf;
  glong row;

  /* Circular buffer for last TOTAL_BANDWIDTH coefficients - much smaller than full solution */
  sol_buf = g_new0 (gdouble, TOTAL_BANDWIDTH);

  /* Back-substitution: compute coefficients and accumulate derivative contributions */
  for (row = n_cols - 1; row >= 0; row--)
  {
    NcmSBesselOdeSolverRow *r = &op->matrix_rows[row];
    gdouble sum               = g_array_index (op->c, gdouble, row);
    const gdouble diag        = r->data[0] + _ncm_sbessel_bc_row (r, row);
    const gdouble row_sign    = (row % 2) == 0 ? 1.0 : -1.0;
    const gdouble kd          = (gdouble) row;
    const gdouble k_squared   = kd * kd;
    const glong buffer_pos    = row % TOTAL_BANDWIDTH; /* Circular buffer position */
    gdouble c_k;

    /* Precompute circular buffer positions for all coefficients we'll need */
    const glong buf_pos_1 = (row + 1) % TOTAL_BANDWIDTH;
    const glong buf_pos_2 = (row + 2) % TOTAL_BANDWIDTH;
    const glong buf_pos_3 = (row + 3) % TOTAL_BANDWIDTH;
    const glong buf_pos_4 = (row + 4) % TOTAL_BANDWIDTH;
    const glong buf_pos_5 = (row + 5) % TOTAL_BANDWIDTH;
    const glong buf_pos_6 = (row + 6) % TOTAL_BANDWIDTH;
    const glong buf_pos_7 = (row + 7) % TOTAL_BANDWIDTH;
    const glong buf_pos_8 = (row + 8) % TOTAL_BANDWIDTH;

    g_assert_cmpuint (r->col_index, ==, row); /* Banded matrix */

    /* Manually unrolled loop with precomputed buffer positions */
    sum -= r->data[1] * sol_buf[buf_pos_1];
    sum -= r->data[2] * sol_buf[buf_pos_2];
    sum -= r->data[3] * sol_buf[buf_pos_3];
    sum -= r->data[4] * sol_buf[buf_pos_4];
    sum -= r->data[5] * sol_buf[buf_pos_5];
    sum -= r->data[6] * sol_buf[buf_pos_6];
    sum -= r->data[7] * sol_buf[buf_pos_7];
    sum -= r->data[8] * sol_buf[buf_pos_8];
    sum -= acc_bc_at_m1 * r->bc_at_m1;
    sum -= acc_bc_at_p1 * r->bc_at_p1;

    c_k = sum / diag;

    /* Store coefficient in circular buffer for future back-substitution steps */
    sol_buf[buffer_pos] = c_k;

    /* Update boundary condition accumulators for next iteration */
    acc_bc_at_m1 += row_sign * c_k;
    acc_bc_at_p1 += c_k;

    /* Accumulate derivative contributions:
     * dy/dt|_{t=-1} = sum_k k^2 * (-1)^(k+1) * c_k
     * dy/dt|_{t=+1} = sum_k k^2 * c_k
     */
    deriv_at_m1    += k_squared * (-row_sign) * c_k; /* u'(-1) */
    deriv_at_p1    += k_squared * c_k;               /* u'(+1) */
    error_estimate += k_squared * fabs (c_k);        /* error estimate */
  }

  /* Convert from t-derivatives to x-derivatives and set output values */
  g_array_set_size (endpoints, 3);
  g_array_index (endpoints, gdouble, 0) = deriv_at_m1 / h;    /* u'(a) = u'(-1) / h */
  g_array_index (endpoints, gdouble, 1) = deriv_at_p1 / h;    /* u'(b) = u'(+1) / h */
  g_array_index (endpoints, gdouble, 2) = error_estimate / h; /* error estimate */

  g_free (sol_buf);
}

/**
 * _ncm_sbessel_apply_rotations_batched:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @col: current column being processed
 * @n_ell: number of ell values
 *
 * Applies Givens rotations to eliminate subdiagonal elements for all ell values.
 */
static inline void
_ncm_sbessel_apply_rotations_batched (NcmSBesselOdeOperator *op, glong col, guint n_ell)
{
  glong i;
  guint l_idx;

  #pragma omp simd

  for (i = ROWS_TO_ROTATE; i > 0; i--)
  {
    const glong r1_index = col + i - 1;
    const glong r2_index = col + i;

    /* Apply Givens rotations for all ell values */
    #pragma omp simd

    for (l_idx = 0; l_idx < n_ell; l_idx++)
    {
      const glong r1_idx_batch   = r1_index * n_ell + l_idx;
      const glong r2_idx_batch   = r2_index * n_ell + l_idx;
      const glong rot_idx        = col * ROWS_TO_ROTATE * n_ell + (ROWS_TO_ROTATE - i) * n_ell + l_idx;
      NcmSBesselOdeSolverRow *r1 = &op->matrix_rows[r1_idx_batch];
      NcmSBesselOdeSolverRow *r2 = &op->matrix_rows[r2_idx_batch];
      gdouble *c1                = &g_array_index (op->c, gdouble, r1_idx_batch);
      gdouble *c2                = &g_array_index (op->c, gdouble, r2_idx_batch);
      gdouble *rot_ptr           = &op->rotation_params[2 * rot_idx];

      _ncm_sbessel_apply_givens (col, r1, r2, c1, c2, rot_ptr);
    }
  }
}

/**
 * _ncm_sbessel_check_convergence_batched:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @col: current column being processed
 * @n_ell: number of ell values
 * @max_c_A: (inout): maximum coefficient-to-diagonal ratio
 * @quiet_cols: (inout): count of consecutive columns with small coefficients
 * @solution_order: (inout): current solution order
 * @total_rows: (inout): total number of rows
 * @max_solution_order: maximum allowed solution order
 *
 * Checks convergence across all ell values and handles adaptive solution order increase.
 *
 * Returns: TRUE if converged (early stop condition met), FALSE otherwise
 */
static inline gboolean
_ncm_sbessel_check_convergence_batched (NcmSBesselOdeOperator *op, glong col, guint n_ell,
                                        gdouble *max_c_A, guint *quiet_cols,
                                        guint *solution_order,
                                        guint max_solution_order)
{
  gdouble max_Acol = 0.0;
  guint l_idx;

  #pragma omp simd

  for (l_idx = 0; l_idx < n_ell; l_idx++)
  {
    const glong row_idx         = col * n_ell + l_idx;
    NcmSBesselOdeSolverRow *row = &op->matrix_rows[row_idx];
    const double diag           = row->data[0] + _ncm_sbessel_bc_row (row, col);
    const double c_col          = g_array_index (op->c, gdouble, row_idx);
    const double Acol           = fabs (c_col / diag);

    if (Acol > max_Acol)
      max_Acol = Acol;
  }

  if (max_Acol > *max_c_A)
  {
    *max_c_A    = max_Acol;
    *quiet_cols = 0;
  }
  else
  {
    /* max_Acol is small relative to established scale */
    if (max_Acol < *max_c_A * op->tolerance)
      (*quiet_cols)++;
    else
      *quiet_cols = 0;
  }

  if (*quiet_cols >= ROWS_TO_ROTATE + 1)
    return TRUE;  /* SAFE early stop */

  /* If we've reached the end without convergence, increase solution_order and continue */
  if ((col == *solution_order - ROWS_TO_ROTATE - 1) && (*solution_order < max_solution_order))
  {
    guint total_rows;

    *solution_order *= 2;
    total_rows       = *solution_order * n_ell;
    _ensure_matrix_rows_capacity (op, total_rows);

    if (total_rows + ROWS_TO_ROTATE * n_ell > op->c->len)
      g_array_set_size (op->c, total_rows + ROWS_TO_ROTATE * n_ell);
  }

  return FALSE;
}

/**
 * _ncm_sbessel_ode_operator_setup_initial_rows_batched:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @n_ell: number of ell values to process
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @solution_order: initial solution order
 *
 * Sets up the first ROWS_TO_ROTATE rows (boundary conditions) for all ell values
 * and initializes the RHS vectors. This function is called only when starting a
 * fresh diagonalization in batched mode.
 */
static void
_ncm_sbessel_ode_operator_setup_initial_rows_batched (NcmSBesselOdeOperator *op, const guint n_ell, GArray *rhs, guint solution_order)
{
  const gint ell_min      = op->ell_min;
  const guint total_rows  = solution_order * n_ell;
  const gdouble *rhs_data = (gdouble *) rhs->data;
  guint l_idx;
  glong i;

  _ensure_matrix_rows_capacity (op, total_rows);

  if (total_rows + ROWS_TO_ROTATE * n_ell > op->c->len)
    g_array_set_size (op->c, total_rows + ROWS_TO_ROTATE * n_ell);

  /* Add boundary condition rows for all ell values */
  for (i = 0; i < ROWS_TO_ROTATE; i++)
  {
    NcmSBesselOdeSolverRow *row = &op->matrix_rows[i * n_ell];

    _ncm_sbessel_create_row_batched (op, row, i, ell_min, n_ell);

    #pragma omp simd

    for (l_idx = 0; l_idx < n_ell; l_idx++)
    {
      const glong row_idx = i * n_ell + l_idx;

      g_array_index (op->c, gdouble, row_idx) = rhs_data[i];
    }
  }
}

/**
 * _ncm_sbessel_ode_operator_diagonalize_batched:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @n_ell: number of ell values to process
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 *
 * Diagonalizes the operator using adaptive QR decomposition for multiple ell values.
 * This function applies Givens rotations to transform the system into upper triangular form
 * and applies the same rotations to the RHS vectors. The transformed RHS is stored in
 * op->c and the upper triangular matrix is stored in op->matrix_rows.
 *
 * Returns: the effective number of columns used (may be less than rhs_len due to convergence)
 */
static glong
_ncm_sbessel_ode_operator_diagonalize_batched (NcmSBesselOdeOperator *op, const guint n_ell, GArray *rhs)
{
  const guint rhs_len            = rhs->len;
  guint solution_order           = rhs_len * 2;
  const guint max_solution_order = 1 << 24;
  const gint ell_min             = op->ell_min;
  const gdouble *rhs_data        = (gdouble *) rhs->data;
  glong col                      = 0;
  gdouble max_c_A                = 0.0;
  guint quiet_cols               = 0;
  const guint first_loop_len     = rhs_len - ROWS_TO_ROTATE;
  gboolean converged             = FALSE;
  guint l_idx;

  g_assert_cmpuint (n_ell, >, 0);
  g_assert_cmpuint (rhs_len, >, ROWS_TO_ROTATE);

  /* Setup initial rows and RHS */
  _ncm_sbessel_ode_operator_setup_initial_rows_batched (op, n_ell, rhs, solution_order);

  for (col = 0; col < first_loop_len; col++)
  {
    glong new_row               = col + ROWS_TO_ROTATE;
    NcmSBesselOdeSolverRow *row = &op->matrix_rows[new_row * n_ell];

    /* Ensure rotation storage capacity */
    _ensure_rotation_capacity (op, (col + 1) * ROWS_TO_ROTATE * n_ell);

    _ncm_sbessel_create_row_batched (op, row, new_row, ell_min, n_ell);

    #pragma omp simd

    for (l_idx = 0; l_idx < n_ell; l_idx++)
    {
      const glong row_idx = new_row * n_ell + l_idx;

      g_array_index (op->c, gdouble, row_idx) = rhs_data[new_row];
    }

    _ncm_sbessel_apply_rotations_batched (op, col, n_ell);

    /* Check convergence across all ell values */
    if (_ncm_sbessel_check_convergence_batched (op, col, n_ell, &max_c_A, &quiet_cols,
                                                &solution_order, max_solution_order))
    {
      converged = TRUE;
      break; /* SAFE early stop */
    }
  }

  if (!converged)
  {
    for ( ; col < solution_order; col++)
    {
      glong new_row               = col + ROWS_TO_ROTATE;
      NcmSBesselOdeSolverRow *row = &op->matrix_rows[new_row * n_ell];

      /* Ensure rotation storage capacity */
      _ensure_rotation_capacity (op, (col + 1) * ROWS_TO_ROTATE * n_ell);

      _ncm_sbessel_create_row_batched (op, row, new_row, ell_min, n_ell);

      #pragma omp simd

      for (l_idx = 0; l_idx < n_ell; l_idx++)
      {
        const glong row_idx = new_row * n_ell + l_idx;

        g_array_index (op->c, gdouble, row_idx) = 0.0;
      }

      _ncm_sbessel_apply_rotations_batched (op, col, n_ell);

      /* Check convergence across all ell values */
      if (_ncm_sbessel_check_convergence_batched (op, col, n_ell, &max_c_A, &quiet_cols,
                                                  &solution_order, max_solution_order))
        break;  /* SAFE early stop */
    }
  }

  /* Warn if we exhausted max_solution_order without convergence */
  if ((solution_order >= max_solution_order) && (quiet_cols < ROWS_TO_ROTATE + 1))
    g_warning ("_ncm_sbessel_ode_solver_diagonalize_batched: "
               "reached max_solution_order=%u without convergence (quiet_cols=%u, needed %d). "
               "Results may be inaccurate.",
               max_solution_order, quiet_cols, ROWS_TO_ROTATE + 1);


  return col;
}

/**
 * _ncm_sbessel_ode_operator_build_solution_batched:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @n_cols: number of columns in the solution (from diagonalization)
 * @n_ell: number of ell values
 * @solutions: output array of solution matrices
 *
 * Builds the full Chebyshev coefficient solution by back-substitution on the
 * upper triangular system. Assumes _ncm_sbessel_ode_solver_diagonalize_batched
 * has been called first.
 *
 * Returns: (transfer full): solution matrix where each row is the solution for one ell value
 */
static void
_ncm_sbessel_ode_operator_build_solution_batched (NcmSBesselOdeOperator *op, glong n_cols, guint n_ell, GArray *solutions)
{
  gdouble * restrict sol_data;
  glong row;
  guint l_idx;

  /* Ensure accumulator arrays have sufficient capacity */
  _ensure_acc_bc_capacity (op, n_ell);

  /* Zero out accumulators */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_ell; l_idx++)
  {
    op->acc_bc_at_m1[l_idx] = 0.0;
    op->acc_bc_at_p1[l_idx] = 0.0;
  }

  g_array_set_size (solutions, (n_cols + TOTAL_BANDWIDTH) * n_ell);
  sol_data = (gdouble *) solutions->data;
  memset (sol_data, 0, (n_cols + TOTAL_BANDWIDTH) * n_ell * sizeof (gdouble));

  /* Flipped loop order for better cache locality - process all ell values for each row */
  for (row = n_cols - 1; row >= 0; row--)
  {
    const glong row_base_idx = row * n_ell;
    const gdouble row_sign   = (row % 2) == 0 ? 1.0 : -1.0;

    #pragma omp simd

    for (l_idx = 0; l_idx < n_ell; l_idx++)
    {
      const glong row_idx                       = row_base_idx + l_idx;
      const NcmSBesselOdeSolverRow * restrict r = &op->matrix_rows[row_idx];
      const gdouble * restrict r_data           = r->data;
      const glong base_sol_idx                  = l_idx * n_cols;
      const gdouble * restrict sol_row          = &sol_data[base_sol_idx];
      gdouble sum                               = g_array_index (op->c, gdouble, row_idx);
      const gdouble diag                        = r_data[0] + _ncm_sbessel_bc_row ((NcmSBesselOdeSolverRow *) r, row);
      gdouble sol;

      g_assert_cmpuint (r->col_index, ==, row); /* Banded matrix */

      /* Manually unrolled loop for better performance */
      sum -= r_data[1] * sol_row[row + 1];
      sum -= r_data[2] * sol_row[row + 2];
      sum -= r_data[3] * sol_row[row + 3];
      sum -= r_data[4] * sol_row[row + 4];
      sum -= r_data[5] * sol_row[row + 5];
      sum -= r_data[6] * sol_row[row + 6];
      sum -= r_data[7] * sol_row[row + 7];
      sum -= r_data[8] * sol_row[row + 8];

      sum -= op->acc_bc_at_m1[l_idx] * r->bc_at_m1;
      sum -= op->acc_bc_at_p1[l_idx] * r->bc_at_p1;

      sol                          = sum / diag;
      sol_data[base_sol_idx + row] = sol;

      op->acc_bc_at_m1[l_idx] += row_sign * sol;
      op->acc_bc_at_p1[l_idx] += sol;
    }
  }

  g_array_set_size (solutions, (n_cols) * n_ell);
}

/**
 * _ncm_sbessel_ode_solver_compute_endpoints_batched:
 * @op: a #NcmSBesselOdeOperator (for matrix and RHS storage)
 * @n_cols: number of columns in the solution (from diagonalization)
 * @n_ell: number of ell values
 * @endpoints: (out): matrix with 3 columns per ell: [u'(a), u'(b), error]
 *
 * Computes endpoint derivatives u'(a) and u'(b) and error estimates directly from
 * the diagonalized system without building the full solution matrix. This is much more
 * efficient when only endpoint information is needed, as it computes coefficients
 * on-the-fly during back-substitution and accumulates their contributions to the
 * derivatives without storing the full coefficient array.
 *
 * This function writes the results into the provided @endpoints matrix, which should
 * have at least n_ell rows and exactly 3 columns. Each row corresponds to one ell value,
 * with columns for u'(a), u'(b), and error estimate.
 *
 * Assumes _ncm_sbessel_ode_solver_diagonalize_batched has been called first.
 *
 * Returns: (transfer full): matrix with 3 columns per ell: [u'(a), u'(b), error]
 */
static void
_ncm_sbessel_ode_operator_compute_endpoints_batched (NcmSBesselOdeOperator *op, glong n_cols, const guint n_ell, GArray *endpoints)
{
  const gdouble h = op->half_len;
  gdouble * restrict endp_data;
  glong row;
  guint l_idx;

  /* Circular buffer for last TOTAL_BANDWIDTH coefficients - much smaller than full solution */
  _ensure_solution_batched_capacity (op, TOTAL_BANDWIDTH * n_ell);
  memset (op->solution_batched, 0, sizeof (gdouble) * TOTAL_BANDWIDTH * n_ell);

  g_array_set_size (endpoints, n_ell * 3);

  /* Result matrix: n_ell rows x 3 columns [u'(a), u'(b), error] */
  endp_data = (gdouble *) endpoints->data;

  /* Initialize endpoint derivative accumulators and error */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_ell; l_idx++)
  {
    endp_data[l_idx * 3 + 0] = 0.0; /* u'(-1) accumulator */
    endp_data[l_idx * 3 + 1] = 0.0; /* u'(+1) accumulator */
    endp_data[l_idx * 3 + 2] = 0.0; /* error accumulator */
  }

  /* Ensure accumulator arrays have sufficient capacity */
  _ensure_acc_bc_capacity (op, n_ell);

  /* Zero out boundary condition accumulators */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_ell; l_idx++)
  {
    op->acc_bc_at_m1[l_idx] = 0.0;
    op->acc_bc_at_p1[l_idx] = 0.0;
  }

  /* Back-substitution: compute coefficients and accumulate derivative contributions */
  for (row = n_cols - 1; row >= 0; row--)
  {
    const glong row_base_idx = row * n_ell;
    const gdouble row_sign   = (row % 2) == 0 ? 1.0 : -1.0;
    const gdouble kd         = (gdouble) row;
    const gdouble k_squared  = kd * kd;
    const glong buffer_pos   = row % TOTAL_BANDWIDTH; /* Circular buffer position */
    /* Precompute circular buffer positions for all coefficients we'll need */
    const glong buf_pos_1 = (row + 1) % TOTAL_BANDWIDTH;
    const glong buf_pos_2 = (row + 2) % TOTAL_BANDWIDTH;
    const glong buf_pos_3 = (row + 3) % TOTAL_BANDWIDTH;
    const glong buf_pos_4 = (row + 4) % TOTAL_BANDWIDTH;
    const glong buf_pos_5 = (row + 5) % TOTAL_BANDWIDTH;
    const glong buf_pos_6 = (row + 6) % TOTAL_BANDWIDTH;
    const glong buf_pos_7 = (row + 7) % TOTAL_BANDWIDTH;
    const glong buf_pos_8 = (row + 8) % TOTAL_BANDWIDTH;

    #pragma omp simd

    for (l_idx = 0; l_idx < n_ell; l_idx++)
    {
      const glong row_idx                       = row_base_idx + l_idx;
      const NcmSBesselOdeSolverRow * restrict r = &op->matrix_rows[row_idx];
      const gdouble * restrict r_data           = r->data;
      const glong base_buffer_idx               = l_idx * TOTAL_BANDWIDTH;
      const gdouble * restrict sol_buf          = &op->solution_batched[base_buffer_idx];
      gdouble sum                               = g_array_index (op->c, gdouble, row_idx);
      const gdouble diag                        = r_data[0] + _ncm_sbessel_bc_row ((NcmSBesselOdeSolverRow *) r, row);
      gdouble c_k;

      g_assert_cmpuint (r->col_index, ==, row); /* Banded matrix */

      /* Manually unrolled loop with precomputed buffer positions */
      sum -= r_data[1] * sol_buf[buf_pos_1];
      sum -= r_data[2] * sol_buf[buf_pos_2];
      sum -= r_data[3] * sol_buf[buf_pos_3];
      sum -= r_data[4] * sol_buf[buf_pos_4];
      sum -= r_data[5] * sol_buf[buf_pos_5];
      sum -= r_data[6] * sol_buf[buf_pos_6];
      sum -= r_data[7] * sol_buf[buf_pos_7];
      sum -= r_data[8] * sol_buf[buf_pos_8];
      sum -= op->acc_bc_at_m1[l_idx] * r->bc_at_m1;
      sum -= op->acc_bc_at_p1[l_idx] * r->bc_at_p1;

      c_k = sum / diag;

      /* Store coefficient in circular buffer for future back-substitution steps */
      op->solution_batched[base_buffer_idx + buffer_pos] = c_k;

      /* Update boundary condition accumulators for next iteration */
      op->acc_bc_at_m1[l_idx] += row_sign * c_k;
      op->acc_bc_at_p1[l_idx] += c_k;

      /* Accumulate derivative contributions:
       * du/dt|_{t=-1} = sum_k k^2 * (-1)^(k+1) * c_k
       * du/dt|_{t=+1} = sum_k k^2 * c_k
       */
      endp_data[l_idx * 3 + 0] += k_squared * (-row_sign) * c_k; /* u'(-1) */
      endp_data[l_idx * 3 + 1] += k_squared * c_k;               /* u'(+1) */
      endp_data[l_idx * 3 + 2] += k_squared * fabs (c_k);        /* error estimate */
    }
  }

  /* Convert from t-derivatives to x-derivatives and finalize error */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_ell; l_idx++)
  {
    endp_data[l_idx * 3 + 0] /= h; /* u'(a) = u'(-1) / h */
    endp_data[l_idx * 3 + 1] /= h; /* u'(b) = u'(+1) / h */
    endp_data[l_idx * 3 + 2] /= h; /* error estimate */
  }
}

/**
 * _ncm_sbessel_ode_operator_solve_batched_internal:
 * @op: a #NcmSBesselOdeOperator
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @n_ell: number of ell values to solve for (ell = ell_min, ell_min+1, ..., ell_min+n_ell-1)
 * @solutions: array of solution matrices, one per ell value
 *
 * Internal batched solver implementation. Can be specialized at compile time
 * when n_ell is known at compile time for better optimization.
 *
 * This function uses the factored implementation: first diagonalizes the operator,
 * then builds the full solution.
 *
 * Returns: (transfer full): solution matrix where each row is the solution for one ell value
 */
void
_ncm_sbessel_ode_operator_solve_batched_internal (NcmSBesselOdeOperator *op, GArray *rhs, const guint n_ell, GArray *solutions)
{
  /* Step 1: Diagonalize the operator using QR decomposition */
  const glong n_cols = _ncm_sbessel_ode_operator_diagonalize_batched (op, n_ell, rhs);

  /* Step 2: Build the full solution by back-substitution */
  _ncm_sbessel_ode_operator_build_solution_batched (op, n_cols, n_ell, solutions);
}

/**
 * _ncm_sbessel_ode_operator_solve_endpoints_batched_internal:
 * @op: a #NcmSBesselOdeOperator
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @n_ell: number of ell values to solve for (ell = ell_min, ell_min+1, ..., ell_min+n_ell-1)
 * @solutions: array of solution matrices, one per ell value
 *
 * Internal batched solver for endpoint computations. Diagonalizes the operator and computes
 * endpoint derivatives and error estimates without building the full solution matrix.
 * Returns: (transfer full): matrix with 3 columns per ell: [u'(a), u'(b), error]
 */
void
_ncm_sbessel_ode_operator_solve_endpoints_batched_internal (NcmSBesselOdeOperator *op, GArray *rhs, const guint n_ell, GArray *solutions)
{
  /* Step 1: Diagonalize the operator using QR decomposition */
  const glong n_cols = _ncm_sbessel_ode_operator_diagonalize_batched (op, n_ell, rhs);

  /* Step 2: Compute endpoint derivatives and error estimates directly from diagonalized system */
  _ncm_sbessel_ode_operator_compute_endpoints_batched (op, n_cols, n_ell, solutions);
}

/* Specialized batched solvers for common sizes - enables better compiler optimizations */

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_operator_solve_batched_2 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_batched_internal (op, rhs, 2, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_operator_solve_batched_4 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_batched_internal (op, rhs, 4, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_operator_solve_batched_8 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_batched_internal (op, rhs, 8, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_operator_solve_batched_16 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_batched_internal (op, rhs, 16, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_operator_solve_batched_32 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_batched_internal (op, rhs, 32, solutions);
}

/* Endpoint computations */

static void
_ncm_sbessel_ode_operator_solve_endpoints_batched_2 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_endpoints_batched_internal (op, rhs, 2, solutions);
}

static void
_ncm_sbessel_ode_operator_solve_endpoints_batched_4 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_endpoints_batched_internal (op, rhs, 4, solutions);
}

static void
_ncm_sbessel_ode_operator_solve_endpoints_batched_8 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_endpoints_batched_internal (op, rhs, 8, solutions);
}

static void
_ncm_sbessel_ode_operator_solve_endpoints_batched_16 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_endpoints_batched_internal (op, rhs, 16, solutions);
}

static void
_ncm_sbessel_ode_operator_solve_endpoints_batched_32 (NcmSBesselOdeOperator *op, GArray *rhs, GArray *solutions)
{
  _ncm_sbessel_ode_operator_solve_endpoints_batched_internal (op, rhs, 32, solutions);
}

/**
 * ncm_sbessel_ode_operator_get_interval:
 * @op: a #NcmSBesselOdeOperator
 * @a: (out): left endpoint of the interval
 * @b: (out): right endpoint of the interval
 *
 * Gets the interval [a, b] on which the operator is defined.
 */
void
ncm_sbessel_ode_operator_get_interval (NcmSBesselOdeOperator *op, gdouble *a, gdouble *b)
{
  g_assert_nonnull (a);
  g_assert_nonnull (b);

  *a = op->mid_point - op->half_len;
  *b = op->mid_point + op->half_len;
}

/**
 * ncm_sbessel_ode_operator_get_ell_range:
 * @op: a #NcmSBesselOdeOperator
 * @ell_min: (out): minimum angular momentum value
 * @ell_max: (out): maximum angular momentum value
 *
 * Gets the range of angular momentum values [ell_min, ell_max] for which
 * the operator solves. For single-ell operators, ell_min == ell_max.
 */
void
ncm_sbessel_ode_operator_get_ell_range (NcmSBesselOdeOperator *op, gint *ell_min, gint *ell_max)
{
  g_assert_nonnull (ell_min);
  g_assert_nonnull (ell_max);

  *ell_min = op->ell_min;
  *ell_max = op->ell_max;
}

/**
 * ncm_sbessel_ode_operator_get_tolerance:
 * @op: a #NcmSBesselOdeOperator
 *
 * Gets the convergence tolerance used by the operator for adaptive QR decomposition.
 *
 * Returns: the convergence tolerance
 */
gdouble
ncm_sbessel_ode_operator_get_tolerance (NcmSBesselOdeOperator *op)
{
  return op->tolerance;
}

/**
 * ncm_sbessel_ode_operator_solve:
 * @op: a #NcmSBesselOdeOperator
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @solution: (out callee-allocates) (transfer full) (element-type gdouble): solution vector (Chebyshev coefficients)
 * @solution_len: (out): length of solution per ell value
 *
 * Solves the ODE using adaptive QR decomposition with ultraspherical spectral methods.
 * The algorithm grows the matrix size until convergence is achieved (error < tolerance).
 *
 * For operators with a single ell value, the @solution array has length @solution_len.
 * For operators with multiple ell values (n_ell), the @solution array has length @solution_len * n_ell,
 * with solutions stored contiguously: first solution at indices [0, solution_len), second at
 * [solution_len, 2*solution_len), etc.
 *
 * The @solution parameter supports two usage patterns:
 *
 * - **Python/convenience usage**: Pass `solution` pointing to NULL (`*solution == NULL`).
 *   A new #GArray will be allocated and returned. The `(out callee-allocates)` annotation
 *   ensures Python bindings automatically use this mode.
 *
 * - **C optimization**: Pass `solution` pointing to a pre-allocated #GArray (`*solution != NULL`).
 *   The existing array will be reused (cleared and refilled), avoiding repeated allocation/deallocation
 *   in performance-critical loops.
 *
 */
void
ncm_sbessel_ode_operator_solve (NcmSBesselOdeOperator *op, GArray *rhs, GArray **solution, gsize *solution_len)
{
  glong n_cols;

  if (*solution == NULL)
    *solution = g_array_new (FALSE, FALSE, sizeof (gdouble));

  if (op->n_ell == 1)
  {
    n_cols = _ncm_sbessel_ode_solver_diagonalize (op, rhs);

    _ncm_sbessel_ode_solver_build_solution (op, n_cols, *solution);
  }
  else
  {
    switch (op->n_ell)
    {
      case 2:
        _ncm_sbessel_ode_operator_solve_batched_2 (op, rhs, *solution);
        break;

      case 4:
        _ncm_sbessel_ode_operator_solve_batched_4 (op, rhs, *solution);
        break;

      case 8:
        _ncm_sbessel_ode_operator_solve_batched_8 (op, rhs, *solution);
        break;

      case 16:
        _ncm_sbessel_ode_operator_solve_batched_16 (op, rhs, *solution);
        break;

      case 32:
        _ncm_sbessel_ode_operator_solve_batched_32 (op, rhs, *solution);
        break;

      default:
        _ncm_sbessel_ode_operator_solve_batched_internal (op, rhs, op->n_ell, *solution);
        break;
    }

    n_cols = (*solution)->len / op->n_ell;
  }

  if (solution_len != NULL)
    *solution_len = n_cols;
}

/**
 * ncm_sbessel_ode_operator_solve_endpoints:
 * @op: a #NcmSBesselOdeOperator
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @endpoints: (out callee-allocates) (transfer full) (element-type gdouble): array with 3*n_ell elements storing [u'(a), u'(b), error] for each ell value
 *
 * Efficiently computes only the endpoint derivatives u'(a) and u'(b)
 * without building the full solution vector. This is much more efficient when only endpoint
 * information is needed (e.g., for integral computations via Green's identity).
 *
 * For an operator with n_ell angular momentum values (ell_min to ell_max), the output
 * array contains n_ell triplets. Each triplet at index i (where i = 0 to n_ell-1)
 * corresponds to ell = ell_min + i, with elements:
 * - endpoints[3*i + 0] = u'(a) for ell_min + i
 * - endpoints[3*i + 1] = u'(b) for ell_min + i
 * - endpoints[3*i + 2] = error estimate for ell_min + i
 *
 * The function first diagonalizes the operator using adaptive QR decomposition, then
 * performs back-substitution while accumulating the contributions to the endpoint
 * derivatives on-the-fly, avoiding the memory allocation and computation cost of the
 * full solution.
 *
 * The @endpoints parameter supports two usage patterns:
 *
 * - **Python/convenience usage**: Pass `endpoints` pointing to NULL (`*endpoints == NULL`).
 *   A new #GArray will be allocated and returned. The `(out callee-allocates)` annotation
 *   ensures Python bindings automatically use this mode.
 *
 * - **C optimization**: Pass `endpoints` pointing to a pre-allocated #GArray (`*endpoints != NULL`).
 *   The existing array will be reused (cleared and refilled), avoiding repeated allocation/deallocation
 *   in performance-critical loops.
 *
 */
void
ncm_sbessel_ode_operator_solve_endpoints (NcmSBesselOdeOperator *op, GArray *rhs, GArray **endpoints)
{
  if (*endpoints == NULL)
    *endpoints = g_array_new (FALSE, FALSE, sizeof (gdouble));

  if (op->n_ell == 1)
  {
    const glong n_cols = _ncm_sbessel_ode_solver_diagonalize (op, rhs);

    _ncm_sbessel_ode_solver_compute_endpoints (op, n_cols, *endpoints);

    return;
  }
  else
  {
    switch (op->n_ell)
    {
      case 2:
        _ncm_sbessel_ode_operator_solve_endpoints_batched_2 (op, rhs, *endpoints);
        break;

      case 4:
        _ncm_sbessel_ode_operator_solve_endpoints_batched_4 (op, rhs, *endpoints);
        break;

      case 8:
        _ncm_sbessel_ode_operator_solve_endpoints_batched_8 (op, rhs, *endpoints);
        break;

      case 16:
        _ncm_sbessel_ode_operator_solve_endpoints_batched_16 (op, rhs, *endpoints);
        break;

      case 32:
        _ncm_sbessel_ode_operator_solve_endpoints_batched_32 (op, rhs, *endpoints);
        break;

      default:
        _ncm_sbessel_ode_operator_solve_endpoints_batched_internal (op, rhs, op->n_ell, *endpoints);
        break;
    }
  }
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
                                              const gdouble a, const gdouble b, guint ell,
                                              guint nrows, guint ncols,
                                              gdouble *data, gboolean colmajor)
{
  NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, a, b, ell, ell);
  NcmSBesselOdeSolverRow row_mem;
  NcmSBesselOdeSolverRow *row = &row_mem;
  guint i;

  /* Fill matrix rows */
  for (i = 0; i < nrows; i++)
  {
    glong k;

    _ncm_sbessel_create_row (op, row, i);

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
  }

  ncm_sbessel_ode_operator_unref (op);
}

/**
 * ncm_sbessel_ode_solver_get_operator_matrix:
 * @solver: a #NcmSBesselOdeSolver
 * @a: left endpoint
 * @b: right endpoint
 * @ell: $\ell$ multipole order
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
ncm_sbessel_ode_solver_get_operator_matrix (NcmSBesselOdeSolver *solver, const gdouble a, const gdouble b, guint ell, gint nrows)
{
  g_assert_cmpint (nrows, >, 2);

  /* Force square matrix */
  const gint ncols = nrows;

  /* Create matrix (row-major by default in NcmMatrix) */
  NcmMatrix *mat = ncm_matrix_new (nrows, ncols);

  ncm_matrix_set_zero (mat);

  /* Fill matrix data */
  _ncm_sbessel_ode_solver_fill_operator_matrix (solver, a, b, ell, nrows, ncols,
                                                ncm_matrix_data (mat), FALSE);

  return mat;
}

/**
 * ncm_sbessel_ode_solver_get_operator_matrix_colmajor:
 * @solver: a #NcmSBesselOdeSolver
 * @a: left endpoint
 * @b: right endpoint
 * @ell: $\ell$ multipole order
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
ncm_sbessel_ode_solver_get_operator_matrix_colmajor (NcmSBesselOdeSolver *solver, const gdouble a, const gdouble b, guint ell, gint nrows)
{
  g_assert_cmpint (nrows, >, 2);
  {
    const gint ncols       = nrows;
    gdouble *data_colmajor = g_new0 (gdouble, nrows * ncols);
    NcmMatrix *mat         = ncm_matrix_new_data_malloc (data_colmajor, nrows, ncols);

    _ncm_sbessel_ode_solver_fill_operator_matrix (solver, a, b, ell, nrows, ncols,
                                                  data_colmajor, TRUE);

    return mat;
  }
}

/**
 * ncm_sbessel_ode_solver_solve_dense:
 * @solver: a #NcmSBesselOdeSolver
 * @a: left endpoint
 * @b: right endpoint
 * @ell: $\ell$ multipole order
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
ncm_sbessel_ode_solver_solve_dense (NcmSBesselOdeSolver *solver, const gdouble a, const gdouble b, guint ell, NcmVector *rhs, gint nrows)
{
  /* Use column-major matrix for LAPACK */
  NcmMatrix *mat   = ncm_sbessel_ode_solver_get_operator_matrix_colmajor (solver, a, b, ell, nrows);
  const gint ncols = ncm_matrix_ncols (mat);

  /* For LAPACK dgesv, we need a square matrix, so use minimum dimension */
  const gint n = GSL_MIN (nrows, ncols);

  /* Prepare RHS vector (truncate/pad if needed) */
  NcmVector *solution = ncm_vector_new (n);

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (solution), ==, 1);
  {
    const gdouble *rhs_data = ncm_vector_data (rhs);
    gdouble *b_data         = ncm_vector_data (solution);
    const gint rhs_len      = (gint) ncm_vector_len (rhs);
    const gint copy_len     = GSL_MIN (n, rhs_len);

    memcpy (b_data, rhs_data, copy_len * sizeof (gdouble));

    for (gint i = copy_len; i < n; i++)
    {
      b_data[i] = 0.0;
    }
  }

  /* Solve the system using LU decomposition */
  gint *ipiv = g_new (gint, n);
  gint ret   = ncm_lapack_dgesv (n, 1,
                                 ncm_matrix_data (mat), n,
                                 ipiv,
                                 ncm_vector_data (solution), n);

  if (ret != 0)
    g_warning ("ncm_sbessel_ode_solver_solve_dense: LAPACK dgesv failed with code %d", ret);

  g_free (ipiv);
  ncm_matrix_free (mat);

  /* Solution is stored in solution */
  return solution;
}

/**
 * ncm_sbessel_ode_solver_peek_spectral:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Returns a pointer to the internal spectral object used by the solver. This allows
 * the user to inspect the internal state of the solver. The spectral object contains
 * information about the Chebyshev nodes, weights, and Chebyshev polynomials.
 *
 * Returns: (transfer none): spectral object
 */
NcmSpectral *
ncm_sbessel_ode_solver_peek_spectral (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return self->spectral;
}


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
  gdouble data[TOTAL_BANDWIDTH]; /* Pointer to row data */
  gdouble bc_at_m1;              /* Boundary condition row 1 coefficient */
  gdouble bc_at_p1;              /* Boundary condition row 2 coefficient */
  glong col_index;               /* Column index of leftmost element in data */
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
  GArray *matrix_rows; /* Array of NcmSBesselOdeSolverRow* */
  GArray *c;           /* Array of gdouble for right-hand side */

  /* Batched matrix storage (adaptive) */
  GArray *matrix_rows_batched; /* Array of NcmSBesselOdeSolverRow for batched operations */
  GArray *c_batched;           /* Array of gdouble for batched right-hand side */
  GArray *acc_bc_at_m1;        /* Array of gdouble for batched boundary condition accumulators at -1 */
  GArray *acc_bc_at_p1;        /* Array of gdouble for batched boundary condition accumulators at +1 */

  /* Solution */
  NcmVector *solution; /* Result vector */

  /* Chebyshev coefficients computation cache */
  guint cheb_N_cached;     /* Cached N value */
  gdouble *cheb_f_vals;    /* Cached function values array */
  fftw_plan cheb_plan_r2r; /* Cached FFTW plan */

  /* Spherical Bessel array with cutoff */
  NcmSFSBesselArray *sba;
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
static void _ncm_sbessel_compute_proj_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_x_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_x2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_d_row (NcmSBesselOdeSolverRow *row, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_x_d_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_d2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_x_d2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static void _ncm_sbessel_compute_x2_d2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff);
static gdouble _ncm_sbessel_bc_row (NcmSBesselOdeSolverRow *row, glong col_index);

static void
ncm_sbessel_ode_solver_init (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  self->l                   = 0;
  self->tolerance           = NCM_SBESSEL_ODE_SOLVER_DEFAULT_TOLERANCE;
  self->max_size            = NCM_SBESSEL_ODE_SOLVER_DEFAULT_MAX_SIZE;
  self->a                   = -1.0;
  self->b                   = 1.0;
  self->half_len            = 1.0;
  self->mid_point           = 0.0;
  self->matrix_rows         = g_array_new (FALSE, FALSE, sizeof (NcmSBesselOdeSolverRow));
  self->c                   = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->matrix_rows_batched = g_array_new (FALSE, FALSE, sizeof (NcmSBesselOdeSolverRow));
  self->c_batched           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->acc_bc_at_m1        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->acc_bc_at_p1        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->solution            = NULL;

  /* Initialize Chebyshev cache */
  self->cheb_N_cached = 0;
  self->cheb_f_vals   = NULL;
  self->cheb_plan_r2r = NULL;
  self->sba           = ncm_sf_sbessel_array_new ();
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

  g_clear_pointer (&self->matrix_rows, g_array_unref);
  g_clear_pointer (&self->c, g_array_unref);
  g_clear_pointer (&self->matrix_rows_batched, g_array_unref);
  g_clear_pointer (&self->c_batched, g_array_unref);
  g_clear_pointer (&self->acc_bc_at_m1, g_array_unref);
  g_clear_pointer (&self->acc_bc_at_p1, g_array_unref);
  ncm_vector_clear (&self->solution);

  /* Clean up Chebyshev cache */
  if (self->cheb_plan_r2r != NULL)
  {
    fftw_destroy_plan (self->cheb_plan_r2r);
    self->cheb_plan_r2r = NULL;
  }

  if (self->cheb_f_vals != NULL)
  {
    fftw_free (self->cheb_f_vals);
    self->cheb_f_vals = NULL;
  }

  self->cheb_N_cached = 0;

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
_ncm_sbessel_compute_proj_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* General formula: three entries with rational function coefficients */
  const gdouble value_0 = coeff / (2.0 * (kd + 1.0));
  const gdouble value_1 = -coeff * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = coeff / (2.0 * (kd + 3.0));

  row->data[offset]     += value_0;
  row->data[offset + 2] += value_1;
  row->data[offset + 4] += value_2;

  /* Special case: k=0 has an additional contribution */
  if (k == 0)
    row->data[0] += coeff * 0.5;
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
_ncm_sbessel_compute_x_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* General formula: up to four entries */
  if (k >= 1)
    row->data[offset - 1] += coeff / (4.0 * (kd + 1.0));

  row->data[offset + 1] += -coeff / (4.0 * (kd + 3.0));
  row->data[offset + 3] += -coeff / (4.0 * (kd + 1.0));
  row->data[offset + 5] += coeff / (4.0 * (kd + 3.0));

  /* Special cases */
  if (k == 0)
    row->data[1 - row->col_index] += coeff * 0.25;  /* Additional 1/4 at column 1 */
  else if (k == 1)
    row->data[0 - row->col_index] += coeff * 0.125;  /* Additional 1/8 at column 0 */
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
_ncm_sbessel_compute_x2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd  = (gdouble) k;
  const gdouble kp2 = kd + 2.0;

  /* General formula: up to five entries */
  if (k >= 2)
    row->data[offset - 2] += coeff / (8.0 * (kd + 1.0));

  row->data[offset]     += coeff / (4.0 * (kd + 1.0) * (kd + 3.0));
  row->data[offset + 2] += -coeff * kp2 / (4.0 * (kd + 1.0) * (kd + 3.0));
  row->data[offset + 4] += -coeff / (4.0 * (kd + 1.0) * (kd + 3.0));
  row->data[offset + 6] += coeff / (8.0 * (kd + 3.0));

  /* Special cases */
  if (k == 0)
  {
    row->data[0 - row->col_index] += coeff / 12.0; /* Additional 1/12 at column 0 */
    row->data[2 - row->col_index] += coeff / 8.0;  /* Additional 1/8 at column 2 */
  }
  else if (k == 1)
  {
    row->data[1 - row->col_index] += coeff / 16.0; /* Additional 1/16 at column 1 */
  }
  else if (k == 2)
  {
    row->data[0 - row->col_index] += coeff / 24.0; /* Additional 1/24 at column 0 */
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
_ncm_sbessel_compute_d_row (NcmSBesselOdeSolverRow *row, glong offset, gdouble coeff)
{
  /* Two non-zero entries with opposite signs */
  row->data[offset + 1] += coeff;
  row->data[offset + 3] += -coeff;
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
_ncm_sbessel_compute_x_d_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* Three non-zero entries with rational function coefficients */
  const gdouble value_0 = coeff * kd / (2.0 * (kd + 1.0));
  const gdouble value_1 = coeff * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = -coeff * (kd + 4.0) / (2.0 * (kd + 3.0));

  row->data[offset]     += value_0;
  row->data[offset + 2] += value_1;
  row->data[offset + 4] += value_2;
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
_ncm_sbessel_compute_d2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  row->data[offset + 2] += coeff * 2.0 * (kd + 2.0);
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
_ncm_sbessel_compute_x_d2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* Two non-zero entries in this row */
  row->data[offset + 1] += coeff * kd;
  row->data[offset + 3] += coeff * (kd + 4.0);
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
_ncm_sbessel_compute_x2_d2_row (NcmSBesselOdeSolverRow *row, glong k, glong offset, gdouble coeff)
{
  const gdouble kd  = (gdouble) k;
  const gdouble kp2 = kd + 2.0;

  /* Three non-zero entries with rational function coefficients */
  const gdouble value_0 = coeff * kd * (kd - 1.0) / (2.0 * (kd + 1.0));
  const gdouble value_1 = coeff * kp2 * (kp2 * kp2 - 3.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = coeff * (kd + 4.0) * (kd + 5.0) / (2.0 * (kd + 3.0));

  row->data[offset]     += value_0;
  row->data[offset + 2] += value_1;
  row->data[offset + 4] += value_2;
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
static void
_ncm_sbessel_create_row_bc_at_m1 (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverRow *row)
{
  _row_reset (row, 1.0, 0.0, 0);
}

/**
 * _ncm_sbessel_create_row_bc_at_p1:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Creates boundary condition row for u(+1) = 0.
 * All data entries are zero, bc_at_m1 = 0.0, bc_at_p1 = 1.0.
 */
static void
_ncm_sbessel_create_row_bc_at_p1 (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverRow *row)
{
  _row_reset (row, 0.0, 1.0, 0);
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
static void
_ncm_sbessel_create_row_operator (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverRow *row, glong row_index)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble m                         = self->mid_point;
  const gdouble h                         = self->half_len;
  const gdouble h2                        = h * h;
  const gdouble m2                        = m * m;
  const gint l                            = self->l;
  const glong k                           = row_index;
  const gdouble llp1                      = (gdouble) l * (l + 1);
  const glong left_col                    = GSL_MAX (0, k - 2); /* Leftmost column index */

  _row_reset (row, 0.0, 0.0, left_col);

  /* Compute offset once for all operators */
  const glong offset = k - row->col_index;

  /* Second derivative term: (m^2/h^2) d^2 + (2m/h) x d^2 + x^2 d^2 */
  _ncm_sbessel_compute_d2_row (row, k, offset, m2 / h2);
  _ncm_sbessel_compute_x_d2_row (row, k, offset, 2.0 * m / h);
  _ncm_sbessel_compute_x2_d2_row (row, k, offset, 1.0);

  /* First derivative term: (2m/h) d + 2 x d */
  _ncm_sbessel_compute_d_row (row, offset, 2.0 * m / h);
  _ncm_sbessel_compute_x_d_row (row, k, offset, 2.0);

  /* Identity term: (m^2 - l(l+1)) I + 2m h x + h^2 x^2 */
  _ncm_sbessel_compute_proj_row (row, k, offset, m2 - llp1);
  _ncm_sbessel_compute_x_row (row, k, offset, 2.0 * m * h);
  _ncm_sbessel_compute_x2_row (row, k, offset, h2);
}

/**
 * _ncm_sbessel_create_row_operator_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @row: array of rows, one for each l value
 * @row_index: row index (0-based, corresponding to k in the Gegenbauer basis)
 * @lmin: minimum l value
 * @n_l: number of l values to process
 *
 * Creates differential operator rows for multiple l values at once using an optimized
 * incremental algorithm. The key insight is that operator rows for consecutive l values
 * differ only in the l(l+1) term:
 *
 * - All derivative operators (d², x d², x² d², d, x d) are l-independent
 * - Identity terms (x, x²) are l-independent
 * - Only the projection term m² - l(l+1) depends on l
 *
 * The optimization uses the recurrence: (l+1)(l+2) - l(l+1) = 2(l+1), allowing
 * incremental updates by adding -2l for each successive l value. This reduces
 * O(n_l × operations) to O(operations + n_l) complexity.
 *
 * Algorithm:
 * 1. Compute all l-independent terms once into row[0]
 * 2. Add l-dependent term (m² - lmin(lmin+1)) to row[0]
 * 3. For each subsequent row: memcpy row[0] and add incremental correction -2l
 */
static void
_ncm_sbessel_create_row_operator_batched (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverRow *row, glong row_index, gint lmin, guint n_l)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble m                         = self->mid_point;
  const gdouble h                         = self->half_len;
  const gdouble h2                        = h * h;
  const gdouble m2                        = m * m;
  const gdouble lmin_lmin_p_1             = (gdouble) lmin * (lmin + 1);
  const glong k                           = row_index;
  const glong left_col                    = GSL_MAX (0, k - 2); /* Leftmost column index */
  gdouble llp1                            = lmin_lmin_p_1;
  gdouble l                               = (gdouble) lmin;
  guint i;

  /* Initialize first row with all l-independent terms */
  _row_reset (&row[0], 0.0, 0.0, left_col);

  /* Compute offset once */
  const glong offset = k - row[0].col_index;

  /* Second derivative term: (m^2/h^2) d^2 + (2m/h) x d^2 + x^2 d^2 - l-independent */
  _ncm_sbessel_compute_d2_row (&row[0], k, offset, m2 / h2);
  _ncm_sbessel_compute_x_d2_row (&row[0], k, offset, 2.0 * m / h);
  _ncm_sbessel_compute_x2_d2_row (&row[0], k, offset, 1.0);

  /* First derivative term: (2m/h) d + 2 x d - l-independent */
  _ncm_sbessel_compute_d_row (&row[0], offset, 2.0 * m / h);
  _ncm_sbessel_compute_x_d_row (&row[0], k, offset, 2.0);

  /* Identity term: 2m h x + h^2 x^2 - l-independent part */
  _ncm_sbessel_compute_proj_row (&row[0], k, offset, m2 - llp1);
  _ncm_sbessel_compute_x_row (&row[0], k, offset, 2.0 * m * h);
  _ncm_sbessel_compute_x2_row (&row[0], k, offset, h2);

  /* Copy template to remaining rows and add incremental l-dependent correction */
  for (i = 1; i < n_l; i++)
  {
    l += 1.0;
    memcpy (&row[i], &row[i - 1], sizeof (NcmSBesselOdeSolverRow));
    _ncm_sbessel_compute_proj_row (&row[i], k, offset, -2.0 * l);
  }
}

/**
 * _ncm_sbessel_create_row:
 * @solver: a #NcmSBesselOdeSolver
 * @row_index: row index (first two rows are boundary conditions, rest are operators)
 *
 * Creates the appropriate row based on the row index.
 * Rows 0-1 are boundary conditions, rows >= 2 are differential operators.
 */
static void
_ncm_sbessel_create_row (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverRow *row, glong row_index)
{
  switch (row_index)
  {
    case 0:
      _ncm_sbessel_create_row_bc_at_m1 (solver, row);
      break;

    case 1:
      _ncm_sbessel_create_row_bc_at_p1 (solver, row);
      break;

    default:
      _ncm_sbessel_create_row_operator (solver, row, row_index - 2);
      break;
  }
}

/**
 * _ncm_sbessel_create_row_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @row: array of rows, one for each l value
 * @row_index: row index (first two rows are boundary conditions, rest are operators)
 * @lmin: minimum l value
 * @n_l: number of l values to process
 *
 * Creates the appropriate rows for multiple l values at once.
 * The i-th element row[i] will contain the row for l = lmin + i.
 * Rows 0-1 are boundary conditions (same for all l), rows >= 2 are differential operators.
 */
static void
_ncm_sbessel_create_row_batched (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverRow *row, glong row_index, gint lmin, guint n_l)
{
  guint i;

  switch (row_index)
  {
    case 0:

      /* Boundary condition at -1: same for all l values */
      for (i = 0; i < n_l; i++)
        _ncm_sbessel_create_row_bc_at_m1 (solver, &row[i]);

      break;

    case 1:

      /* Boundary condition at +1: same for all l values */
      for (i = 0; i < n_l; i++)
        _ncm_sbessel_create_row_bc_at_p1 (solver, &row[i]);

      break;

    default:
      /* Differential operator: optimized batch computation */
      _ncm_sbessel_create_row_operator_batched (solver, row, row_index - 2, lmin, n_l);
      break;
  }
}

__attribute__ ((optimize ("no-math-errno", "no-trapping-math")))

static inline gdouble
_compute_hypot (gdouble a, gdouble b)
{
  const gdouble abs_a = fabs (a);
  const gdouble abs_b = fabs (b);

  if (abs_a == 0.0)
    return abs_b;

  if (abs_b == 0.0)
    return abs_a;

  if (abs_a > abs_b)
  {
    const gdouble ratio = abs_b / abs_a;

    return abs_a * sqrt (1.0 + ratio * ratio);
  }
  else
  {
    const gdouble ratio = abs_a / abs_b;

    return abs_b * sqrt (1.0 + ratio * ratio);
  }
}

/**
 * _ncm_sbessel_apply_givens:
 * @solver: a #NcmSBesselOdeSolver
 * @pivot_col: column index of the pivot
 * @r1: pointer to first row (pivot row)
 * @r2: pointer to second row (row to eliminate)
 * @c1: pointer to right-hand side element for row1
 * @c2: pointer to right-hand side element for row2
 *
 * Applies Givens rotation to eliminate the entry at (row2, col=pivot_col).
 * The rotation transforms the rows as:
 *   new_row1 = c*row1 + s*row2
 *   new_row2 = -s*row1 + c*row2
 * where c and s are chosen to zero out element (row2, pivot_col).
 */
static inline __attribute__ ((hot)) void

_ncm_sbessel_apply_givens (NcmSBesselOdeSolver *solver, glong pivot_col, NcmSBesselOdeSolverRow * restrict r1, NcmSBesselOdeSolverRow * restrict r2, gdouble * restrict c1, gdouble * restrict c2)
{
  const gdouble a_val = r1->data[0] + _ncm_sbessel_bc_row (r1, pivot_col);
  const gdouble b_val = r2->data[0] + _ncm_sbessel_bc_row (r2, pivot_col);
  const gdouble norm  = _compute_hypot (a_val, b_val);

  if (__builtin_expect ((norm < 1.0e-100), 0))
  {
    /* Entry already zero - just shift r2 without rotation */
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
    const gdouble inv_norm     = 1.0 / norm;
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


    r2->col_index++; /* Shift row2's column index right by 1 */
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
  const gdouble *rhs_data                 = ncm_vector_data (rhs);
  glong col                               = 0;
  gdouble max_c_A                         = 0.0;
  guint quiet_cols                        = 0;
  glong i, last_row_index;

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);

  if (rhs_len > self->matrix_rows->len)
    g_array_set_size (self->matrix_rows, rhs_len);

  g_array_set_size (self->c, rhs_len + ROWS_TO_ROTATE);

  ncm_vector_clear (&self->solution);

  /* Add boundary condition rows */
  for (i = 0; i < ROWS_TO_ROTATE + 1; i++)
  {
    NcmSBesselOdeSolverRow *row = &g_array_index (self->matrix_rows, NcmSBesselOdeSolverRow, i);

    g_array_index (self->c, gdouble, i) = rhs_data[i];
    _ncm_sbessel_create_row (solver, row, i);
  }

  last_row_index = ROWS_TO_ROTATE;

  for (col = 0; col < rhs_len; col++)
  {
    const glong last_row_for_col = GSL_MIN (ROWS_TO_ROTATE, rhs_len - col - 1);

    for (i = last_row_for_col; i > 0; i--)
    {
      const glong r1_index = col + i - 1;
      const glong r2_index = col + i;

      if (r2_index > last_row_index)
      {
        NcmSBesselOdeSolverRow *row = &g_array_index (self->matrix_rows, NcmSBesselOdeSolverRow, ++last_row_index);

        _ncm_sbessel_create_row (solver, row, last_row_index);
        g_array_index (self->c, gdouble, last_row_index) = rhs_data[last_row_index];
      }

      {
        NcmSBesselOdeSolverRow *r1 = &g_array_index (self->matrix_rows, NcmSBesselOdeSolverRow, r1_index);
        NcmSBesselOdeSolverRow *r2 = &g_array_index (self->matrix_rows, NcmSBesselOdeSolverRow, r2_index);
        gdouble *c1                = &g_array_index (self->c, gdouble, r1_index);
        gdouble *c2                = &g_array_index (self->c, gdouble, r2_index);

        _ncm_sbessel_apply_givens (solver, col, r1, r2, c1, c2);
      }
    }

    {
      NcmSBesselOdeSolverRow *row =
        &g_array_index (self->matrix_rows, NcmSBesselOdeSolverRow, col);

      const double diag  = row->data[0] + _ncm_sbessel_bc_row (row, col);
      const double c_col = g_array_index (self->c, gdouble, col);
      const double Acol  = fabs (c_col / diag);

      if (Acol > max_c_A)
      {
        max_c_A    = Acol;
        quiet_cols = 0;
      }
      else
      {
        /* Acol is small relative to established scale */
        if (Acol < max_c_A * DBL_EPSILON * 1.0e-1)
          quiet_cols++;
        else
          quiet_cols = 0;
      }

      if (quiet_cols >= ROWS_TO_ROTATE + 1)
        break;  /* SAFE early stop */
    }
  }

  {
    /* Back substitution to solve R x = c */
    const glong n        = col; /* since col <= rhs_len */
    gdouble acc_bc_at_m1 = 0.0;
    gdouble acc_bc_at_p1 = 0.0;
    gdouble *sol_ptr;
    glong row;

    self->solution = ncm_vector_new (n);
    sol_ptr        = ncm_vector_data (self->solution);

    for (row = n - 1; row >= 0; row--)
    {
      NcmSBesselOdeSolverRow *r = &g_array_index (self->matrix_rows, NcmSBesselOdeSolverRow, row);
      gdouble sum               = g_array_index (self->c, gdouble, row);
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

  return ncm_vector_ref (self->solution);
}

/**
 * _ncm_sbessel_ode_solver_solve_batched_internal:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @lmin: minimum l value
 * @n_l: number of l values to solve for (l = lmin, lmin+1, ..., lmin+n_l-1)
 *
 * Internal batched solver implementation. Can be specialized at compile time
 * when n_l is known at compile time for better optimization.
 *
 * Returns: (transfer full): solution matrix where each row is the solution for one l value
 */
static NcmMatrix *
_ncm_sbessel_ode_solver_solve_batched_internal (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin, guint n_l)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const guint rhs_len                     = ncm_vector_len (rhs);
  const guint total_rows                  = rhs_len * n_l;
  const gdouble *rhs_data                 = ncm_vector_data (rhs);
  NcmMatrix *solution                     = NULL;
  glong col                               = 0;
  gdouble max_c_A                         = 0.0;
  guint quiet_cols                        = 0;
  glong i, last_row_index;
  guint l_idx;

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);
  g_assert_cmpuint (n_l, >, 0);

  /* Resize arrays only if necessary */
  if (total_rows > self->matrix_rows_batched->len)
    g_array_set_size (self->matrix_rows_batched, total_rows);

  if (total_rows + ROWS_TO_ROTATE * n_l > self->c_batched->len)
    g_array_set_size (self->c_batched, total_rows + ROWS_TO_ROTATE * n_l);

  /* Add boundary condition rows for all l values */
  for (i = 0; i < ROWS_TO_ROTATE + 1; i++)
  {
    NcmSBesselOdeSolverRow *row = &g_array_index (self->matrix_rows_batched, NcmSBesselOdeSolverRow, i * n_l);

    _ncm_sbessel_create_row_batched (solver, row, i, lmin, n_l);

    #pragma omp simd

    for (l_idx = 0; l_idx < n_l; l_idx++)
    {
      const glong row_idx = i * n_l + l_idx;

      g_array_index (self->c_batched, gdouble, row_idx) = rhs_data[i];
    }
  }

  last_row_index = ROWS_TO_ROTATE;

  for (col = 0; col < rhs_len; col++)
  {
    const glong last_row_for_col = GSL_MIN (ROWS_TO_ROTATE, rhs_len - col - 1);

    for (i = last_row_for_col; i > 0; i--)
    {
      const glong r1_index = col + i - 1;
      const glong r2_index = col + i;

      /* Create new rows if needed */
      if (r2_index > last_row_index)
      {
        NcmSBesselOdeSolverRow *row = &g_array_index (self->matrix_rows_batched, NcmSBesselOdeSolverRow, r2_index * n_l);

        _ncm_sbessel_create_row_batched (solver, row, r2_index, lmin, n_l);

        #pragma omp simd

        for (l_idx = 0; l_idx < n_l; l_idx++)
        {
          const glong row_idx = r2_index * n_l + l_idx;

          g_array_index (self->c_batched, gdouble, row_idx) = rhs_data[r2_index];
        }

        last_row_index = r2_index;
      }

      /* Apply Givens rotations for all l values */
      #pragma omp simd

      for (l_idx = 0; l_idx < n_l; l_idx++)
      {
        const glong r1_idx_batch   = r1_index * n_l + l_idx;
        const glong r2_idx_batch   = r2_index * n_l + l_idx;
        NcmSBesselOdeSolverRow *r1 = &g_array_index (self->matrix_rows_batched, NcmSBesselOdeSolverRow, r1_idx_batch);
        NcmSBesselOdeSolverRow *r2 = &g_array_index (self->matrix_rows_batched, NcmSBesselOdeSolverRow, r2_idx_batch);
        gdouble *c1                = &g_array_index (self->c_batched, gdouble, r1_idx_batch);
        gdouble *c2                = &g_array_index (self->c_batched, gdouble, r2_idx_batch);

        _ncm_sbessel_apply_givens (solver, col, r1, r2, c1, c2);
      }
    }

    /* Check convergence across all l values */
    {
      gdouble max_Acol = 0.0;

      #pragma omp simd

      for (l_idx = 0; l_idx < n_l; l_idx++)
      {
        const glong row_idx         = col * n_l + l_idx;
        NcmSBesselOdeSolverRow *row = &g_array_index (self->matrix_rows_batched, NcmSBesselOdeSolverRow, row_idx);
        const double diag           = row->data[0] + _ncm_sbessel_bc_row (row, col);
        const double c_col          = g_array_index (self->c_batched, gdouble, row_idx);
        const double Acol           = fabs (c_col / diag);

        if (Acol > max_Acol)
          max_Acol = Acol;
      }

      if (max_Acol > max_c_A)
      {
        max_c_A    = max_Acol;
        quiet_cols = 0;
      }
      else
      {
        /* max_Acol is small relative to established scale */
        if (max_Acol < max_c_A * DBL_EPSILON * 1.0e-1)
          quiet_cols++;
        else
          quiet_cols = 0;
      }

      if (quiet_cols >= ROWS_TO_ROTATE + 1)
        break;  /* SAFE early stop */
    }
  }

  /* Back substitution to solve R x = c for all l values */
  {
    const glong n = col; /* since col <= rhs_len */
    gdouble * restrict sol_data;
    glong row;

    /* Resize accumulator arrays if necessary */
    g_array_set_size (self->acc_bc_at_m1, n_l);
    g_array_set_size (self->acc_bc_at_p1, n_l);

    /* Zero out accumulators */
    #pragma omp simd

    for (l_idx = 0; l_idx < n_l; l_idx++)
    {
      g_array_index (self->acc_bc_at_m1, gdouble, l_idx) = 0.0;
      g_array_index (self->acc_bc_at_p1, gdouble, l_idx) = 0.0;
    }

    solution = ncm_matrix_new (n_l, n);
    sol_data = ncm_matrix_data (solution);

    /* Flipped loop order for better cache locality - process all l values for each row */
    for (row = n - 1; row >= 0; row--)
    {
      const glong row_base_idx = row * n_l;
      const gdouble row_sign   = (row % 2) == 0 ? 1.0 : -1.0;

      #pragma omp simd

      for (l_idx = 0; l_idx < n_l; l_idx++)
      {
        const glong row_idx       = row_base_idx + l_idx;
        NcmSBesselOdeSolverRow *r = &g_array_index (self->matrix_rows_batched, NcmSBesselOdeSolverRow, row_idx);
        gdouble sum               = g_array_index (self->c_batched, gdouble, row_idx);
        glong width               = GSL_MIN (TOTAL_BANDWIDTH, n - row);
        const gdouble diag        = r->data[0] + _ncm_sbessel_bc_row (r, row);
        gdouble sol;
        glong j;

        g_assert_cmpuint (r->col_index, ==, row); /* Banded matrix */

        for (j = 1; j < width; j++)
        {
          const glong col_idx = r->col_index + j;
          const gdouble a_ij  = r->data[j];

          sum -= a_ij * sol_data[l_idx * n + col_idx];
        }

        sum -= g_array_index (self->acc_bc_at_m1, gdouble, l_idx) * r->bc_at_m1;
        sum -= g_array_index (self->acc_bc_at_p1, gdouble, l_idx) * r->bc_at_p1;

        sol                       = sum / diag;
        sol_data[l_idx * n + row] = sol;

        g_array_index (self->acc_bc_at_m1, gdouble, l_idx) += row_sign * sol;
        g_array_index (self->acc_bc_at_p1, gdouble, l_idx) += sol;
      }
    }
  }

  return solution;
}

/* Specialized batched solvers for common sizes - enables better compiler optimizations */

static inline __attribute__ ((always_inline)) NcmMatrix *

_ncm_sbessel_ode_solver_solve_batched_2 (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin)
{
  return _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 2);
}

static inline __attribute__ ((always_inline)) NcmMatrix *

_ncm_sbessel_ode_solver_solve_batched_4 (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin)
{
  return _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 4);
}

static inline __attribute__ ((always_inline)) NcmMatrix *

_ncm_sbessel_ode_solver_solve_batched_8 (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin)
{
  return _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 8);
}

static inline __attribute__ ((always_inline)) NcmMatrix *

_ncm_sbessel_ode_solver_solve_batched_16 (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin)
{
  return _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 16);
}

static inline __attribute__ ((always_inline)) NcmMatrix *

_ncm_sbessel_ode_solver_solve_batched_32 (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin)
{
  return _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 32);
}

/**
 * ncm_sbessel_ode_solver_solve_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @lmin: minimum l value
 * @n_l: number of l values to solve for (l = lmin, lmin+1, ..., lmin+n_l-1)
 *
 * Solves the ODE for multiple l values simultaneously using batched operations.
 * Allocates rhs_len*n_l matrix rows and processes them in batches for efficiency.
 * Dispatches to specialized implementations for common sizes (8, 16, 32) for better
 * compiler optimizations including loop unrolling and vectorization.
 *
 * Returns: (transfer full): solution matrix where each row is the solution for one l value
 */
NcmMatrix *
ncm_sbessel_ode_solver_solve_batched (NcmSBesselOdeSolver *solver, NcmVector *rhs, gint lmin, guint n_l)
{
  /* Dispatch to specialized versions for common sizes */
  switch (n_l)
  {
    case 2:
      return _ncm_sbessel_ode_solver_solve_batched_2 (solver, rhs, lmin);

    case 4:
      return _ncm_sbessel_ode_solver_solve_batched_4 (solver, rhs, lmin);

    case 8:
      return _ncm_sbessel_ode_solver_solve_batched_8 (solver, rhs, lmin);

    case 16:
      return _ncm_sbessel_ode_solver_solve_batched_16 (solver, rhs, lmin);

    case 32:
      return _ncm_sbessel_ode_solver_solve_batched_32 (solver, rhs, lmin);

    default:
      return _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, n_l);
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
                                              guint nrows, guint ncols,
                                              gdouble *data, gboolean colmajor)
{
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  guint i;

  /* Fill matrix rows */
  for (i = 0; i < nrows; i++)
  {
    glong k;

    _ncm_sbessel_create_row (solver, row, i);

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

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (b), ==, 1);
  {
    const gdouble *rhs_data = ncm_vector_data (rhs);
    gdouble *b_data         = ncm_vector_data (b);
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
                                 ncm_vector_data (b), n);

  if (ret != 0)
    g_warning ("ncm_sbessel_ode_solver_solve_dense: LAPACK dgesv failed with code %d", ret);

  g_free (ipiv);
  ncm_matrix_free (mat);

  /* Solution is stored in b */
  return b;
}

/**
 * ncm_sbessel_ode_solver_integrate:
 * @solver: a #NcmSBesselOdeSolver
 * @F: (scope call): function to integrate against spherical Bessel function
 * @N: number of Chebyshev nodes for the approximation
 * @user_data: user data for @F
 *
 * Computes the integral $\int_a^b f(x) j_l(x) dx$ using Green's identity
 * with the spherical Bessel ODE. The method:
 *
 * 1. Computes Chebyshev coefficients of f(x) on [a,b]
 * 2. Converts to Gegenbauer $C^{(2)}$ basis
 * 3. Solves the ODE with homogeneous Dirichlet boundary conditions
 * 4. Evaluates derivatives at the endpoints
 * 5. Returns $b^2 j_l(b) y'(b) - a^2 j_l(a) y'(a)$
 *
 * By Green's identity, this equals $\int_a^b f(x) j_l(x) dx$ when y(x)
 * is the solution to the spherical Bessel ODE with right-hand side f(x).
 *
 * Returns: the value of the integral $\int_a^b f(x) j_l(x) dx$
 */
gdouble
ncm_sbessel_ode_solver_integrate (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverF F, guint N, gpointer user_data)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble a                         = self->a;
  const gdouble b                         = self->b;
  const gdouble h                         = self->half_len;
  const gint l                            = self->l;

  /* Step 1: Compute Chebyshev coefficients for f(x) */
  NcmVector *cheb_coeffs = ncm_sbessel_ode_solver_compute_chebyshev_coeffs (solver, F, a, b, N, user_data);

  /* Step 2: Convert to Gegenbauer C^(2) basis */
  NcmVector *gegen_coeffs = ncm_vector_new (N);

  ncm_sbessel_ode_solver_chebT_to_gegenbauer_lambda2 (cheb_coeffs, gegen_coeffs);

  /* Step 3: Set up RHS with homogeneous boundary conditions */
  NcmVector *rhs = ncm_vector_new (N + 2);

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (gegen_coeffs), ==, 1);
  {
    gdouble *rhs_data                = ncm_vector_data (rhs);
    const gdouble *gegen_coeffs_data = ncm_vector_data (gegen_coeffs);

    rhs_data[0] = 0.0; /* BC at x=a (t=-1) */
    rhs_data[1] = 0.0; /* BC at x=b (t=+1) */

    /* Copy Gegenbauer coefficients to RHS starting from index 2 */
    memcpy (&rhs_data[2], gegen_coeffs_data, N * sizeof (gdouble));
  }

  /* Step 4: Solve the ODE */
  NcmVector *solution = ncm_sbessel_ode_solver_solve (solver, rhs);

  /* Step 5: Compute derivatives at endpoints */
  /* dy/dx = (dy/dt) / h, where h = (b-a)/2 */
  const gdouble y_prime_at_minus1 = ncm_sbessel_ode_solver_chebyshev_deriv (solution, -1.0);
  const gdouble y_prime_at_plus1  = ncm_sbessel_ode_solver_chebyshev_deriv (solution, 1.0);
  const gdouble y_prime_a         = y_prime_at_minus1 / h;
  const gdouble y_prime_b         = y_prime_at_plus1 / h;

  /* Step 6: Evaluate j_l at endpoints */
  const gdouble j_l_a = gsl_sf_bessel_jl (l, a);
  const gdouble j_l_b = gsl_sf_bessel_jl (l, b);

  /* Step 7: Compute integral via Green's identity */
  const gdouble integral = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a;

  /* Clean up */
  ncm_vector_free (cheb_coeffs);
  ncm_vector_free (gegen_coeffs);
  ncm_vector_free (rhs);
  ncm_vector_free (solution);

  return integral;
}

/* Define Gaussian function */
typedef struct
{
  gdouble center;
  gdouble inv_std2;
  gdouble k;
} GaussianData;

gdouble
gaussian_func (gpointer user_data, gdouble x)
{
  GaussianData *gdata = (GaussianData *) user_data;
  const gdouble y     = x / gdata->k;
  const gdouble dx    = y - gdata->center;

  return exp (-dx * dx * gdata->inv_std2 * 0.5);
}

/**
 * ncm_sbessel_ode_solver_integrate_gaussian:
 * @solver: a #NcmSBesselOdeSolver
 * @center: center of the Gaussian
 * @std: standard deviation of the Gaussian
 * @k: scale factor
 * @N: number of Chebyshev nodes
 *
 * Computes the integral $\int_a^b \exp(-(x - \text{center})^2/(2\text{std}^2)) j_l(kx) dx$
 * using the spectral method.
 *
 * Returns: the value of the integral
 */
gdouble
ncm_sbessel_ode_solver_integrate_gaussian (NcmSBesselOdeSolver *solver, gdouble center, gdouble std, gdouble k, guint N)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble a_orig                    = self->a;
  const gdouble b_orig                    = self->b;
  const gdouble inv_std2                  = 1.0 / (std * std);
  GaussianData data                       = { center, inv_std2, k };

  /* Temporarily scale the interval by k for this integration */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig * k, b_orig * k);

  /* Integrate */
  const gdouble result = ncm_sbessel_ode_solver_integrate (solver, gaussian_func, N, &data);

  /* Restore original interval */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig, b_orig);

  return result;
}

/* Define rational function */
typedef struct
{
  gdouble center;
  gdouble inv_std;
  gdouble k;
} RationalData;

gdouble
rational_func (gpointer user_data, gdouble x)
{
  RationalData *rdata = (RationalData *) user_data;
  const gdouble y     = x / rdata->k;
  const gdouble dx    = (y - rdata->center) * rdata->inv_std;
  const gdouble denom = 1.0 + dx * dx;

  return (y * y) / (denom * denom * denom);
}

/**
 * ncm_sbessel_ode_solver_integrate_rational:
 * @solver: a #NcmSBesselOdeSolver
 * @center: center of the rational function
 * @std: width parameter
 * @k: scale factor
 * @N: number of Chebyshev nodes
 *
 * Computes the integral $\int_a^b \frac{x^2}{(1 + ((x-\text{center})/\text{std})^2)^4} j_l(kx) dx$
 * using the spectral method.
 *
 * Returns: the value of the integral
 */
gdouble
ncm_sbessel_ode_solver_integrate_rational (NcmSBesselOdeSolver *solver, gdouble center, gdouble std, gdouble k, guint N)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble a_orig                    = self->a;
  const gdouble b_orig                    = self->b;
  const gdouble inv_std                   = 1.0 / std;
  RationalData data                       = { center, inv_std, k };

  /* Temporarily scale the interval by k for this integration */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig * k, b_orig * k);

  /* Integrate */
  const gdouble result = ncm_sbessel_ode_solver_integrate (solver, rational_func, N, &data);

  /* Restore original interval */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig, b_orig);

  return result;
}

/**
 * ncm_sbessel_ode_solver_integrate_l_range:
 * @solver: a #NcmSBesselOdeSolver
 * @F: (scope call): function to integrate against spherical Bessel function
 * @N: number of Chebyshev nodes for the approximation
 * @lmin: minimum l value
 * @lmax: maximum l value
 * @user_data: user data for @F
 *
 * Computes the integrals $\int_a^b f(x) j_l(x) dx$ for each l from @lmin to @lmax.
 * This is more efficient than calling ncm_sbessel_ode_solver_integrate() multiple times
 * because the RHS (Chebyshev coefficients and Gegenbauer conversion) is computed only once.
 *
 * Returns: (transfer full): a #NcmVector with the integral values for each l from @lmin to @lmax
 */
NcmVector *
ncm_sbessel_ode_solver_integrate_l_range (NcmSBesselOdeSolver *solver, NcmSBesselOdeSolverF F, guint N, gint lmin, gint lmax, gpointer user_data)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble a                         = self->a;
  const gdouble b                         = self->b;
  const gdouble h                         = self->half_len;
  const gint l_orig                       = self->l;
  const gint my_lmax                      = GSL_MIN (lmax, ncm_sf_sbessel_array_eval_ell_cutoff (self->sba, b));
  const gint n_l                          = lmax - lmin + 1;
  gdouble *j_array_a, *j_array_b;

  j_array_a = g_new0 (gdouble, my_lmax + 1);
  j_array_b = g_new0 (gdouble, my_lmax + 1);

  ncm_sf_sbessel_array_eval (self->sba, my_lmax, a, j_array_a);
  ncm_sf_sbessel_array_eval (self->sba, my_lmax, b, j_array_b);

  /* Allocate result vector */
  NcmVector *result    = ncm_vector_new (n_l);
  gdouble *result_data = ncm_vector_data (result);

  ncm_vector_set_zero (result);

  /* Step 1: Compute Chebyshev coefficients for f(x) - done once */
  NcmVector *cheb_coeffs = ncm_sbessel_ode_solver_compute_chebyshev_coeffs (solver, F, a, b, N, user_data);

  /* Step 2: Convert to Gegenbauer C^(2) basis - done once */
  NcmVector *gegen_coeffs = ncm_vector_new (N);

  ncm_sbessel_ode_solver_chebT_to_gegenbauer_lambda2 (cheb_coeffs, gegen_coeffs);

  /* Step 3: Set up RHS with homogeneous boundary conditions - done once */
  NcmVector *rhs = ncm_vector_new (N);

  g_assert_cmpuint (ncm_vector_stride (rhs), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (gegen_coeffs), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (result), ==, 1);
  {
    gdouble *rhs_data                = ncm_vector_data (rhs);
    const gdouble *gegen_coeffs_data = ncm_vector_data (gegen_coeffs);

    rhs_data[0] = 0.0; /* BC at x=a (t=-1) */
    rhs_data[1] = 0.0; /* BC at x=b (t=+1) */

    /* Copy Gegenbauer coefficients to RHS starting from index 2 */
    memcpy (&rhs_data[2], gegen_coeffs_data, (N - 2) * sizeof (gdouble));
  }

  /* Step 4-7: Process l values in blocks for better cache locality */
  {
    const gint block_size = 8; /* Process this many l values at once - tune based on cache size */

    for (gint l_block_start = lmin; l_block_start <= my_lmax; l_block_start += block_size)
    {
      const gint l_block_end = GSL_MIN (l_block_start + block_size - 1, my_lmax);
      const guint n_l_block  = l_block_end - l_block_start + 1;

      /* Solve for all l values in this block simultaneously */
      NcmMatrix *solutions = ncm_sbessel_ode_solver_solve_batched (solver, rhs, l_block_start, n_l_block);

      /* Process each l in the block */
      for (guint i = 0; i < n_l_block; i++)
      {
        const gint l = l_block_start + i;

        /* Get solution row for this l (each row is one solution) */
        const guint sol_len = ncm_matrix_ncols (solutions);

        /* Compute derivatives at endpoints using solution row i */
        gdouble y_prime_at_minus1 = 0.0;
        gdouble y_prime_at_plus1  = 0.0;

        /* Compute Chebyshev derivative at t=-1 and t=+1 */
        for (guint k = 0; k < sol_len; k++)
        {
          const gdouble c_k = ncm_matrix_get (solutions, i, k);
          const gdouble kd  = (gdouble) k;

          /* Derivative at t=-1: sum_k k^2 * (-1)^(k+1) * c_k */
          y_prime_at_minus1 += kd * kd * ((k % 2 == 0) ? -1.0 : 1.0) * c_k;

          /* Derivative at t=+1: sum_k k^2 * c_k */
          y_prime_at_plus1 += kd * kd * c_k;
        }

        const gdouble y_prime_a = y_prime_at_minus1 / h;
        const gdouble y_prime_b = y_prime_at_plus1 / h;

        /* Evaluate j_l at endpoints */
        const gdouble j_l_a = j_array_a[l];
        const gdouble j_l_b = j_array_b[l];

        /* Compute integral via Green's identity */
        const gdouble integral = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a;

        /* Store result */
        result_data[l - lmin] = integral;
      }

      /* Clean up solutions for this block */
      ncm_matrix_free (solutions);
    }
  }

  /* Restore original l value */
  ncm_sbessel_ode_solver_set_l (solver, l_orig);

  /* Clean up shared resources */
  ncm_vector_free (cheb_coeffs);
  ncm_vector_free (gegen_coeffs);
  ncm_vector_free (rhs);

  return result;
}

/**
 * ncm_sbessel_ode_solver_integrate_gaussian_l_range:
 * @solver: a #NcmSBesselOdeSolver
 * @center: center of the Gaussian
 * @std: standard deviation of the Gaussian
 * @k: scale factor
 * @N: number of Chebyshev nodes
 * @lmin: minimum l value
 * @lmax: maximum l value
 *
 * Computes the integrals $\int_a^b \exp(-(x - \text{center})^2/(2\text{std}^2)) j_l(kx) dx$
 * for each l from @lmin to @lmax. This is more efficient than calling
 * ncm_sbessel_ode_solver_integrate_gaussian() multiple times because the RHS
 * is computed only once.
 *
 * Returns: (transfer full): a #NcmVector with the integral values for each l from @lmin to @lmax
 */
NcmVector *
ncm_sbessel_ode_solver_integrate_gaussian_l_range (NcmSBesselOdeSolver *solver, gdouble center, gdouble std, gdouble k, guint N, gint lmin, gint lmax)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble a_orig                    = self->a;
  const gdouble b_orig                    = self->b;
  const gdouble inv_std2                  = 1.0 / (std * std);
  GaussianData data                       = { center, inv_std2, k };

  /* Temporarily scale the interval by k for this integration */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig * k, b_orig * k);

  /* Integrate over l range */
  NcmVector *result = ncm_sbessel_ode_solver_integrate_l_range (solver, gaussian_func, N, lmin, lmax, &data);

  /* Restore original interval */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig, b_orig);

  return result;
}

/**
 * ncm_sbessel_ode_solver_integrate_rational_l_range:
 * @solver: a #NcmSBesselOdeSolver
 * @center: center of the rational function
 * @std: width parameter
 * @k: scale factor
 * @N: number of Chebyshev nodes
 * @lmin: minimum l value
 * @lmax: maximum l value
 *
 * Computes the integrals $\int_a^b \frac{x^2}{(1 + ((x-\text{center})/\text{std})^2)^4} j_l(kx) dx$
 * for each l from @lmin to @lmax. This is more efficient than calling
 * ncm_sbessel_ode_solver_integrate_rational() multiple times because the RHS
 * is computed only once.
 *
 * Returns: (transfer full): a #NcmVector with the integral values for each l from @lmin to @lmax
 */
NcmVector *
ncm_sbessel_ode_solver_integrate_rational_l_range (NcmSBesselOdeSolver *solver, gdouble center, gdouble std, gdouble k, guint N, gint lmin, gint lmax)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble a_orig                    = self->a;
  const gdouble b_orig                    = self->b;
  const gdouble inv_std                   = 1.0 / std;
  RationalData data                       = { center, inv_std, k };

  /* Temporarily scale the interval by k for this integration */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig * k, b_orig * k);

  /* Integrate over l range */
  NcmVector *result = ncm_sbessel_ode_solver_integrate_l_range (solver, rational_func, N, lmin, lmax, &data);

  /* Restore original interval */
  ncm_sbessel_ode_solver_set_interval (solver, a_orig, b_orig);

  return result;
}

/**
 * ncm_sbessel_ode_solver_compute_chebyshev_coeffs:
 * @solver: a #NcmSBesselOdeSolver
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
 * This method caches the FFTW plan and working arrays for efficiency. They are
 * only reallocated if N changes.
 *
 * Returns: (transfer full): vector of Chebyshev coefficients
 */
NcmVector *
ncm_sbessel_ode_solver_compute_chebyshev_coeffs (NcmSBesselOdeSolver *solver,
                                                 NcmSBesselOdeSolverF F,
                                                 gdouble a, gdouble b,
                                                 guint N,
                                                 gpointer user_data)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  guint i;
  NcmVector *coeffs    = ncm_vector_new (N);
  const gdouble mid    = 0.5 * (a + b);
  const gdouble half_h = 0.5 * (b - a);

  /* Reallocate and replan if N has changed */
  if (self->cheb_N_cached != N)
  {
    /* Clean up old resources */
    if (self->cheb_plan_r2r != NULL)
    {
      fftw_destroy_plan (self->cheb_plan_r2r);
      self->cheb_plan_r2r = NULL;
    }

    if (self->cheb_f_vals != NULL)
    {
      fftw_free (self->cheb_f_vals);
      self->cheb_f_vals = NULL;
    }

    /* Allocate new resources */
    self->cheb_f_vals = fftw_malloc (sizeof (gdouble) * N);

    /* Create new FFTW plan */
    ncm_cfg_load_fftw_wisdom ("ncm_sbessel_ode_solver");
    ncm_cfg_lock_plan_fftw ();
    self->cheb_plan_r2r = fftw_plan_r2r_1d (N, self->cheb_f_vals, ncm_vector_data (coeffs),
                                            FFTW_REDFT00, ncm_cfg_get_fftw_default_flag ());
    ncm_cfg_unlock_plan_fftw ();
    ncm_cfg_save_fftw_wisdom ("ncm_sbessel_ode_solver");

    self->cheb_N_cached = N;
  }

  /* Sample function at Chebyshev nodes */
  for (i = 0; i < N; i++)
  {
    const gdouble theta = M_PI * i / (N - 1);
    const gdouble x     = mid + half_h * cos (theta);

    self->cheb_f_vals[i] = F (user_data, x);
  }

  /* Execute FFTW plan (need to update output pointer for this execution) */
  fftw_execute_r2r (self->cheb_plan_r2r, self->cheb_f_vals, ncm_vector_data (coeffs));

  /* Normalize coefficients */
  g_assert_cmpuint (ncm_vector_stride (coeffs), ==, 1);
  {
    gdouble *coeffs_data = ncm_vector_data (coeffs);

    coeffs_data[0]     = coeffs_data[0] / ((N - 1.0) * 2.0);
    coeffs_data[N - 1] = coeffs_data[N - 1] / ((N - 1.0) * 2.0);

    for (i = 1; i < N - 1; i++)
      coeffs_data[i] = coeffs_data[i] / (N - 1.0);
  }

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
  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (g), ==, 1);

  if (N == 0)
    return;

  {
    const gdouble *c_data = ncm_vector_data (c);
    gdouble *g_data       = ncm_vector_data (g);

    memset (g_data, 0, N * sizeof (gdouble));

    /* n = 0 case */
    g_data[0] = c_data[0];

    if (N == 1)
      return;

    /* n = 1 case */
    g_data[1] = c_data[1] * 0.5;

    /* n >= 2 */
    for (i = 2; i < N; i++)
    {
      const gdouble ci = c_data[i];

      g_data[i]     += 0.5 * ci;
      g_data[i - 2] -= 0.5 * ci;
    }
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
  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (g), ==, 1);

  if (N == 0)
    return;

  {
    const gdouble *c_data = ncm_vector_data (c);
    gdouble *g_data       = ncm_vector_data (g);

    /* Zero output vector */
    memset (g_data, 0, N * sizeof (gdouble));

    /* Apply projection formula for each k */
    for (k = 0; k < N; k++)
    {
      const gdouble kd = (gdouble) k;
      gdouble gk       = 0.0;

      /* Special case: k=0 has additional 1/2 * c[0] contribution */
      if (k == 0)
        gk += 0.5 * c_data[0];

      /* First term: c[k] / (2*(k+1)) */
      gk += c_data[k] / (2.0 * (kd + 1.0));

      /* Second term: -(k+2) * c[k+2] / ((k+1)*(k+3)) */
      if (k + 2 < N)
        gk -= (kd + 2.0) * c_data[k + 2] / ((kd + 1.0) * (kd + 3.0));

      /* Third term: c[k+4] / (2*(k+3)) */
      if (k + 4 < N)
        gk += c_data[k + 4] / (2.0 * (kd + 3.0));

      g_data[k] = gk;
    }
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

  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);

  {
    const gdouble *c_data = ncm_vector_data (c);

    /* Endpoint handling: C_n^{(1)}(+/-1) = (n+1)*(+/-1)^n */
    if (fabs (x - 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * (gdouble) (n + 1);

      return sum;
    }

    if (fabs (x + 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * ((n & 1) ? -(gdouble) (n + 1) : (gdouble) (n + 1));

      return sum;
    }

    {
      /* Stable recurrence for interior x */
      gdouble Cnm1 = 1.0; /* U_0 */
      gdouble sum  = c_data[0] * Cnm1;

      if (N == 1)
        return sum;

      gdouble Cn = 2.0 * x; /* U_1 */

      sum += c_data[1] * Cn;

      for (guint n = 1; n < N - 1; n++)
      {
        gdouble Cnp1 = 2.0 * x * Cn - Cnm1; /* U_{n+1} */

        sum += c_data[n + 1] * Cnp1;
        Cnm1 = Cn;
        Cn   = Cnp1;
      }

      return sum;
    }
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

  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);

  {
    const gdouble *c_data = ncm_vector_data (c);

    /* Endpoint handling: C_n^{(2)}(+/-1) = binom(n+3,3)*(+/-1)^n = ((n+1)*(n+2)*(n+3)/6)*(+/-1)^n */
    if (fabs (x - 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

      return sum;
    }

    if (fabs (x + 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
      {
        const gdouble val = (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

        sum += c_data[n] * ((n & 1) ? -val : val);
      }

      return sum;
    }

    {
      /* Stable recurrence for interior x */
      gdouble Cnm1 = 1.0; /* C_0^{(2)} = 1 */
      gdouble sum  = c_data[0] * Cnm1;

      if (N == 1)
        return sum;

      gdouble Cn = 4.0 * x; /* C_1^{(2)} = 4x */

      sum += c_data[1] * Cn;

      for (guint n = 1; n < N - 1; n++)
      {
        /* (n+1) C_{n+1}^{(2)} = 2(n+2)x C_n^{(2)} - (n+3) C_{n-1}^{(2)} */
        gdouble Cnp1 = (2.0 * (gdouble) (n + 2) * x * Cn - (gdouble) (n + 3) * Cnm1) / (gdouble) (n + 1);

        sum += c_data[n + 1] * Cnp1;
        Cnm1 = Cn;
        Cn   = Cnp1;
      }

      return sum;
    }
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

  g_assert_cmpuint (ncm_vector_stride (a), ==, 1);

  {
    const gdouble *a_data = ncm_vector_data (a);

    if (N == 1)
      return a_data[0];

    {
      gdouble b_kplus1 = 0.0;
      gdouble b_kplus2 = 0.0;
      gdouble two_t    = 2.0 * t;
      gint k;

      for (k = (gint) N - 1; k >= 1; k--)
      {
        gdouble b_k = two_t * b_kplus1 - b_kplus2 + a_data[k];

        b_kplus2 = b_kplus1;
        b_kplus1 = b_k;
      }

      return t * b_kplus1 - b_kplus2 + a_data[0];
    }
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

  g_assert_cmpuint (ncm_vector_stride (a), ==, 1);

  {
    const gdouble *a_data = ncm_vector_data (a);

    if (N == 2)
      return a_data[1];


    if (fabs (t - 1.0) < 1.0e-15)
    {
      /* ---- x = +1 ---- */
      gdouble d1  = 0.0;
      gdouble d2  = 0.0;
      gdouble sum = 0.0;

      for (gint k = N - 2; k >= 1; k--)
      {
        gdouble bk = d2 + 2.0 * (k + 1) * a_data[k + 1];

        d2   = d1;
        d1   = bk;
        sum += bk;
      }

      /* b0 has the 1/2 factor */
      gdouble b0 = 0.5 * (d2 + 2.0 * a_data[1]);

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
        gdouble bk = d2 + 2.0 * (k + 1) * a_data[k + 1];

        d2 = d1;
        d1 = bk;

        sum += sign * bk;
        sign = -sign;
      }

      /* b0 has sign +1 */
      gdouble b0 = 0.5 * (d2 + 2.0 * a_data[1]);

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
        gdouble bk = d2 + 2.0 * (k + 1) * a_data[k + 1];

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
        gdouble b0 = 0.5 * (d2 + 2.0 * a_data[1]);

        return t * c1 - c2 + b0;
      }
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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = k;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_proj_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = (k >= 1) ? (k - 1) : 0;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_x_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = (k >= 2) ? (k - 2) : 0;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_x2_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = k + 1;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_d_row (row, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = k;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_x_d_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = k + 2;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_d2_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = k + 1;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_x_d2_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

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
  NcmSBesselOdeSolverRow *row = g_new0 (NcmSBesselOdeSolverRow, 1);
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  guint k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong left_col = k;

    _row_reset (row, 0.0, 0.0, left_col);

    _ncm_sbessel_compute_x2_d2_row (row, k, k - left_col, 1.0);

    for (guint j = 0; j < TOTAL_BANDWIDTH; j++)
    {
      const glong col = row->col_index + j;

      if ((col >= 0) && (col < (glong) N) && (fabs (row->data[j]) > 1.0e-100))
        ncm_matrix_set (mat, k, col, row->data[j]);
    }
  }

  g_free (row);

  return mat;
}


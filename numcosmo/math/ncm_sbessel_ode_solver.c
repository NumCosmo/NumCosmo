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

#define ALIGNMENT 64 /* 64-byte alignment for cache lines */

#define ROW_CORE_SIZE \
        (sizeof (gdouble) * (TOTAL_BANDWIDTH + 2) + sizeof (glong))

#define PADDED_BANDWIDTH \
        ((TOTAL_BANDWIDTH + 7) & ~7) /* round up to multiple of 8 */
#define OPERATOR_EXPONENT 1

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

typedef struct _NcmSBesselOdeSolverPrivate
{
  gint l;            /* Angular momentum parameter */
  gdouble tolerance; /* Convergence tolerance for adaptive QR */
  gdouble a;         /* Left endpoint of interval */
  gdouble b;         /* Right endpoint of interval */
  gdouble half_len;  /* (b-a)/2 - half length of interval */
  gdouble mid_point; /* (a+b)/2 - midpoint of interval */

  /* Matrix storage (adaptive) */
  NcmSBesselOdeSolverRow *matrix_rows; /* Aligned array of NcmSBesselOdeSolverRow */
  gsize matrix_rows_size;              /* Current size of matrix_rows */
  gsize matrix_rows_capacity;          /* Allocated capacity of matrix_rows */
  GArray *c;                           /* Array of gdouble for right-hand side */

  /* Batched matrix storage (adaptive) */
  NcmSBesselOdeSolverRow *matrix_rows_batched; /* Aligned array of NcmSBesselOdeSolverRow for batched operations */
  gsize matrix_rows_batched_size;              /* Current size of matrix_rows_batched */
  gsize matrix_rows_batched_capacity;          /* Allocated capacity of matrix_rows_batched */
  GArray *c_batched;                           /* Array of gdouble for batched right-hand side */

  /* Aligned arrays for batched operations */
  gdouble *solution_batched;       /* Aligned array of gdouble for temporary solution storage in endpoint computation */
  gsize solution_batched_capacity; /* Allocated capacity of solution_batched */
  gdouble *acc_bc_at_m1;           /* Aligned array of gdouble for batched boundary condition accumulators at -1 */
  gsize acc_bc_at_m1_capacity;     /* Allocated capacity of acc_bc_at_m1 */
  gdouble *acc_bc_at_p1;           /* Aligned array of gdouble for batched boundary condition accumulators at +1 */
  gsize acc_bc_at_p1_capacity;     /* Allocated capacity of acc_bc_at_p1 */

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
  PROP_L,
  PROP_TOLERANCE,
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
static gdouble _ncm_sbessel_bc_row (NcmSBesselOdeSolverRow *row, glong col_index);

static void
ncm_sbessel_ode_solver_init (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  self->l                            = 0;
  self->tolerance                    = 0.0;
  self->a                            = -1.0;
  self->b                            = 1.0;
  self->half_len                     = 1.0;
  self->mid_point                    = 0.0;
  self->matrix_rows                  = NULL;
  self->matrix_rows_size             = 0;
  self->matrix_rows_capacity         = 0;
  self->c                            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->matrix_rows_batched          = NULL;
  self->matrix_rows_batched_size     = 0;
  self->matrix_rows_batched_capacity = 0;
  self->c_batched                    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->solution_batched             = NULL;
  self->solution_batched_capacity    = 0;
  self->acc_bc_at_m1                 = NULL;
  self->acc_bc_at_m1_capacity        = 0;
  self->acc_bc_at_p1                 = NULL;
  self->acc_bc_at_p1_capacity        = 0;
  self->solution                     = NULL;
  self->sba                          = ncm_sf_sbessel_array_new ();
  self->spectral                     = ncm_spectral_new ();
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

  if (self->matrix_rows != NULL)
  {
    free (self->matrix_rows);
    self->matrix_rows          = NULL;
    self->matrix_rows_size     = 0;
    self->matrix_rows_capacity = 0;
  }

  g_clear_pointer (&self->c, g_array_unref);

  if (self->matrix_rows_batched != NULL)
  {
    free (self->matrix_rows_batched);
    self->matrix_rows_batched          = NULL;
    self->matrix_rows_batched_size     = 0;
    self->matrix_rows_batched_capacity = 0;
  }

  g_clear_pointer (&self->c_batched, g_array_unref);

  if (self->solution_batched != NULL)
  {
    free (self->solution_batched);
    self->solution_batched          = NULL;
    self->solution_batched_capacity = 0;
  }

  if (self->acc_bc_at_m1 != NULL)
  {
    free (self->acc_bc_at_m1);
    self->acc_bc_at_m1          = NULL;
    self->acc_bc_at_m1_capacity = 0;
  }

  if (self->acc_bc_at_p1 != NULL)
  {
    free (self->acc_bc_at_p1);
    self->acc_bc_at_p1          = NULL;
    self->acc_bc_at_p1_capacity = 0;
  }

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

/* Helper functions for aligned memory management */
static void
_ensure_matrix_rows_capacity (NcmSBesselOdeSolverPrivate *self, gsize required_size)
{
  if (required_size > self->matrix_rows_capacity)
  {
    NcmSBesselOdeSolverRow *new_rows = NULL;
    gsize new_capacity               = required_size;

    if (new_capacity < 16)
      new_capacity = 16;


    if (posix_memalign ((void **) &new_rows, ALIGNMENT, new_capacity * sizeof (NcmSBesselOdeSolverRow)) != 0)
      g_error ("_ensure_matrix_rows_capacity: failed to allocate aligned memory");

    if (self->matrix_rows != NULL)
    {
      memcpy (new_rows, self->matrix_rows, self->matrix_rows_size * sizeof (NcmSBesselOdeSolverRow));
      free (self->matrix_rows);
    }

    self->matrix_rows          = new_rows;
    self->matrix_rows_capacity = new_capacity;
  }

  self->matrix_rows_size = required_size;
}

static void
_ensure_matrix_rows_batched_capacity (NcmSBesselOdeSolverPrivate *self, gsize required_size)
{
  if (required_size > self->matrix_rows_batched_capacity)
  {
    NcmSBesselOdeSolverRow *new_rows = NULL;
    gsize new_capacity               = required_size;

    if (new_capacity < 16)
      new_capacity = 16;


    if (posix_memalign ((void **) &new_rows, ALIGNMENT, new_capacity * sizeof (NcmSBesselOdeSolverRow)) != 0)
      g_error ("_ensure_matrix_rows_batched_capacity: failed to allocate aligned memory");

    if (self->matrix_rows_batched != NULL)
    {
      memcpy (new_rows, self->matrix_rows_batched, self->matrix_rows_batched_size * sizeof (NcmSBesselOdeSolverRow));
      free (self->matrix_rows_batched);
    }

    self->matrix_rows_batched          = new_rows;
    self->matrix_rows_batched_capacity = new_capacity;
  }

  self->matrix_rows_batched_size = required_size;
}

static void
_ensure_solution_batched_capacity (NcmSBesselOdeSolverPrivate *self, gsize required_size)
{
  if (required_size > self->solution_batched_capacity)
  {
    gdouble *new_array = NULL;
    gsize new_capacity = required_size;

    if (new_capacity < 16)
      new_capacity = 16;

    if (posix_memalign ((void **) &new_array, ALIGNMENT, new_capacity * sizeof (gdouble)) != 0)
      g_error ("_ensure_solution_batched_capacity: failed to allocate aligned memory");

    if (self->solution_batched != NULL)
      free (self->solution_batched);

    self->solution_batched          = new_array;
    self->solution_batched_capacity = new_capacity;
  }
}

static void
_ensure_acc_bc_at_m1_capacity (NcmSBesselOdeSolverPrivate *self, gsize required_size)
{
  if (required_size > self->acc_bc_at_m1_capacity)
  {
    gdouble *new_array = NULL;
    gsize new_capacity = required_size;

    if (new_capacity < 16)
      new_capacity = 16;

    if (posix_memalign ((void **) &new_array, ALIGNMENT, new_capacity * sizeof (gdouble)) != 0)
      g_error ("_ensure_acc_bc_at_m1_capacity: failed to allocate aligned memory");

    if (self->acc_bc_at_m1 != NULL)
      free (self->acc_bc_at_m1);

    self->acc_bc_at_m1          = new_array;
    self->acc_bc_at_m1_capacity = new_capacity;
  }
}

static void
_ensure_acc_bc_at_p1_capacity (NcmSBesselOdeSolverPrivate *self, gsize required_size)
{
  if (required_size > self->acc_bc_at_p1_capacity)
  {
    gdouble *new_array = NULL;
    gsize new_capacity = required_size;

    if (new_capacity < 16)
      new_capacity = 16;

    if (posix_memalign ((void **) &new_array, ALIGNMENT, new_capacity * sizeof (gdouble)) != 0)
      g_error ("_ensure_acc_bc_at_p1_capacity: failed to allocate aligned memory");

    if (self->acc_bc_at_p1 != NULL)
      free (self->acc_bc_at_p1);

    self->acc_bc_at_p1          = new_array;
    self->acc_bc_at_p1_capacity = new_capacity;
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
                                                        0.0, 1.0, DBL_EPSILON * 1.0e-1,
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
  gdouble * restrict row_data             = row->data;

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

  /* Identity term: (m^2 - l(l+1)) I + 2m h x + h^2 x^2 + o_e * (o_e - 1) */
  ncm_spectral_compute_proj_row (row_data, k, offset, m2 - llp1 + (OPERATOR_EXPONENT - 1.0) * OPERATOR_EXPONENT);
  ncm_spectral_compute_x_row (row_data, k, offset, 2.0 * m * h);
  ncm_spectral_compute_x2_row (row_data, k, offset, h2);
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
  const glong k                           = row_index;
  const glong left_col                    = GSL_MAX (0, k - 2); /* Leftmost column index */
  gdouble l                               = (gdouble) lmin;
  guint i;

  /* Initialize first row with all l-independent terms */
  _row_reset (&row[0], 0.0, 0.0, left_col);

  /* Compute offset once */
  const glong offset = k - row[0].col_index;

  /* Second derivative term: (m^2/h^2) d^2 + (2m/h) x d^2 + x^2 d^2 - l-independent */
  ncm_spectral_compute_d2_row (row[0].data, k, offset, m2 / h2);
  ncm_spectral_compute_x_d2_row (row[0].data, k, offset, 2.0 * m / h);
  ncm_spectral_compute_x2_d2_row (row[0].data, k, offset, 1.0);

  /* First derivative term: (2m/h) d + 2 x d - l-independent */
#if OPERATOR_EXPONENT != 1
  ncm_spectral_compute_d_row (row[0].data, offset, 2.0 * (1.0 - OPERATOR_EXPONENT) * m / h);
  ncm_spectral_compute_x_d_row (row[0].data, k, offset, 2.0 * (1.0 - OPERATOR_EXPONENT));
#endif /* OPERATOR_EXPONENT != 1 */

  /* Identity term: 2m h x + h^2 x^2 - l-independent part */
  ncm_spectral_compute_x_row (row[0].data, k, offset, 2.0 * m * h);
  ncm_spectral_compute_x2_row (row[0].data, k, offset, h2);

  /* Copy template from row[0] to all other rows, then add l-dependent corrections in one pass */
  for (i = 1; i < n_l; i++)
  {
    row[i] = row[0]; /* Struct copy is faster than memcpy for small structs and better for cache */
  }

  {
    const gdouble kd = (gdouble) k;

    /* General formula: three entries with rational function coefficients */
    const gdouble value_0 = 1.0 / (2.0 * (kd + 1.0));
    const gdouble value_1 = -1.0 * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
    const gdouble value_2 = 1.0 / (2.0 * (kd + 3.0));
    gdouble coeff         = m2 - l * (l + 1.0) + (OPERATOR_EXPONENT - 1.0) * OPERATOR_EXPONENT;

    for (i = 0; i < n_l; i++)
    {
      gdouble * restrict row_data = row[i].data;

      row_data[offset]     += coeff * value_0;
      row_data[offset + 2] += coeff * value_1;
      row_data[offset + 4] += coeff * value_2;

      /* Special case: k=0 has an additional contribution */
      if (k == 0)
        row_data[offset] += coeff * 0.5;

      l     += 1.0;
      coeff -= 2.0 * l; /* Incremental update for next l value */
    }
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

static inline double
_compute_inv_hypot (double a, double b)
{
  return 1.0 / sqrt (a * a + b * b);
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
  const gdouble a_val    = r1->data[0] + _ncm_sbessel_bc_row (r1, pivot_col);
  const gdouble b_val    = r2->data[0] + _ncm_sbessel_bc_row (r2, pivot_col);
  const gdouble inv_norm = _compute_inv_hypot (a_val, b_val);

  if (__builtin_expect ((inv_norm > 1.0e100), 0))
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
 * _ncm_sbessel_ode_solver_diagonalize:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 *
 * Diagonalizes the operator using adaptive QR decomposition.
 * This function applies Givens rotations to transform the system into upper triangular form
 * and applies the same rotations to the RHS vector. The transformed RHS is stored in
 * self->c and the upper triangular matrix is stored in self->matrix_rows.
 *
 * Returns: the effective number of columns used (may be less than rhs_len due to convergence)
 */
static glong
_ncm_sbessel_ode_solver_diagonalize (NcmSBesselOdeSolver *solver, GArray *rhs)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const guint rhs_len                     = rhs->len;
  guint solution_order                    = rhs_len * 2;
  const guint max_solution_order          = 1 << 24;
  const gdouble *rhs_data                 = (gdouble *) rhs->data;
  glong col                               = 0;
  gdouble max_c_A                         = 0.0;
  guint quiet_cols                        = 0;
  glong i, last_row_index;

  _ensure_matrix_rows_capacity (self, solution_order);

  g_array_set_size (self->c, solution_order + ROWS_TO_ROTATE);

  /* Add boundary condition rows */
  for (i = 0; i < ROWS_TO_ROTATE + 1; i++)
  {
    NcmSBesselOdeSolverRow *row = &self->matrix_rows[i];

    g_array_index (self->c, gdouble, i) = rhs_data[i];
    _ncm_sbessel_create_row (solver, row, i);
  }

  last_row_index = ROWS_TO_ROTATE;

  for (col = 0; col < solution_order; col++)
  {
    const glong last_row_for_col = GSL_MIN (ROWS_TO_ROTATE, solution_order - col - 1);

    for (i = last_row_for_col; i > 0; i--)
    {
      const glong r1_index = col + i - 1;
      const glong r2_index = col + i;

      if (r2_index > last_row_index)
      {
        NcmSBesselOdeSolverRow *row = &self->matrix_rows[++last_row_index];

        _ncm_sbessel_create_row (solver, row, last_row_index);

        if (last_row_index < rhs_len)
          g_array_index (self->c, gdouble, last_row_index) = rhs_data[last_row_index];
        else
          g_array_index (self->c, gdouble, last_row_index) = 0.0;
      }

      {
        NcmSBesselOdeSolverRow *r1 = &self->matrix_rows[r1_index];
        NcmSBesselOdeSolverRow *r2 = &self->matrix_rows[r2_index];
        gdouble *c1                = &g_array_index (self->c, gdouble, r1_index);
        gdouble *c2                = &g_array_index (self->c, gdouble, r2_index);

        _ncm_sbessel_apply_givens (solver, col, r1, r2, c1, c2);
      }
    }

    {
      NcmSBesselOdeSolverRow *row = &self->matrix_rows[col];

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
        if (Acol < max_c_A * self->tolerance)
          quiet_cols++;
        else
          quiet_cols = 0;
      }

      if (quiet_cols >= ROWS_TO_ROTATE + 1)
        break;  /* SAFE early stop */

      /* If we've reached the end without convergence, increase solution_order and continue */
      if ((col == solution_order - ROWS_TO_ROTATE - 1) && (solution_order < max_solution_order))
      {
        solution_order *= 2;
        _ensure_matrix_rows_capacity (self, solution_order);

        if (solution_order + ROWS_TO_ROTATE > self->c->len)
          g_array_set_size (self->c, solution_order + ROWS_TO_ROTATE);
      }
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
 * @solver: a #NcmSBesselOdeSolver
 * @n_cols: number of columns in the solution (from diagonalization)
 *
 * Builds the full Chebyshev coefficient solution by back-substitution on the
 * upper triangular system. Assumes _ncm_sbessel_ode_solver_diagonalize
 * has been called first.
 *
 * Returns: (transfer full): solution vector (Chebyshev coefficients)
 */
static void
_ncm_sbessel_ode_solver_build_solution (NcmSBesselOdeSolver *solver, glong n_cols, GArray *solution)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  gdouble acc_bc_at_m1                    = 0.0;
  gdouble acc_bc_at_p1                    = 0.0;
  gdouble * restrict sol_ptr;
  glong row;

  g_array_set_size (solution, n_cols + TOTAL_BANDWIDTH);
  sol_ptr = (gdouble *) solution->data;
  memset (sol_ptr, 0, sizeof (gdouble) * (n_cols + TOTAL_BANDWIDTH));

  for (row = n_cols - 1; row >= 0; row--)
  {
    NcmSBesselOdeSolverRow *r = &self->matrix_rows[row];
    gdouble sum               = g_array_index (self->c, gdouble, row);
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
 * ncm_sbessel_ode_solver_solve:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @solution: (out callee-allocates) (transfer full) (element-type gdouble): solution vector (Chebyshev coefficients)
 *
 * Solves the ODE using adaptive QR decomposition with ultraspherical spectral methods.
 * The algorithm grows the matrix size until convergence is achieved (error < tolerance).
 *
 */
void
ncm_sbessel_ode_solver_solve (NcmSBesselOdeSolver *solver, GArray *rhs, GArray **solution)
{
  const glong n_cols = _ncm_sbessel_ode_solver_diagonalize (solver, rhs);

  if (solution != NULL)
    *solution = g_array_new (FALSE, FALSE, sizeof (gdouble));

  _ncm_sbessel_ode_solver_build_solution (solver, n_cols, *solution);
}

/**
 * _ncm_sbessel_ode_solver_compute_endpoints:
 * @solver: a #NcmSBesselOdeSolver
 * @n_cols: number of columns in the solution (from diagonalization)
 * @deriv_a: (out): derivative at point a, y'(a)
 * @deriv_b: (out): derivative at point b, y'(b)
 * @error: (out): error estimate
 *
 * Computes endpoint derivatives y'(a) and y'(b) and error estimate directly from the
 * diagonalized system without building the full solution vector. This is much more
 * efficient when only endpoint information is needed, as it computes coefficients
 * on-the-fly during back-substitution and accumulates their contributions to the
 * derivatives without storing the full coefficient array.
 *
 * Assumes _ncm_sbessel_ode_solver_diagonalize has been called first.
 */
static void
_ncm_sbessel_ode_solver_compute_endpoints (NcmSBesselOdeSolver *solver, glong n_cols, gdouble *deriv_a, gdouble *deriv_b, gdouble *error)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble h                         = self->half_len;
  gdouble acc_bc_at_m1                    = 0.0;
  gdouble acc_bc_at_p1                    = 0.0;
  gdouble deriv_at_m1                     = 0.0; /* y'(-1) accumulator */
  gdouble deriv_at_p1                     = 0.0; /* y'(+1) accumulator */
  gdouble error_estimate                  = 0.0; /* error accumulator */
  gdouble *sol_buf;
  glong row;

  /* Circular buffer for last TOTAL_BANDWIDTH coefficients - much smaller than full solution */
  sol_buf = g_new0 (gdouble, TOTAL_BANDWIDTH);

  /* Back-substitution: compute coefficients and accumulate derivative contributions */
  for (row = n_cols - 1; row >= 0; row--)
  {
    NcmSBesselOdeSolverRow *r = &self->matrix_rows[row];
    gdouble sum               = g_array_index (self->c, gdouble, row);
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
    deriv_at_m1    += k_squared * (-row_sign) * c_k; /* y'(-1) */
    deriv_at_p1    += k_squared * c_k;               /* y'(+1) */
    error_estimate += k_squared * fabs (c_k);        /* error estimate */
  }

  /* Convert from t-derivatives to x-derivatives and set output values */
  *deriv_a = deriv_at_m1 / h;    /* y'(a) = y'(-1) / h */
  *deriv_b = deriv_at_p1 / h;    /* y'(b) = y'(+1) / h */
  *error   = error_estimate / h; /* error estimate */

  g_free (sol_buf);
}

/**
 * ncm_sbessel_ode_solver_solve_endpoints:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @deriv_a: (out): derivative at point a, y'(a)
 * @deriv_b: (out): derivative at point b, y'(b)
 * @error: (out): error estimate
 *
 * Efficiently computes only the endpoint derivatives y'(a) and y'(b)
 * without building the full solution vector. This is much more efficient when only endpoint
 * information is needed (e.g., for integral computations via Green's identity).
 *
 * The function first diagonalizes the operator using adaptive QR decomposition, then
 * performs back-substitution while accumulating the contributions to the endpoint
 * derivatives on-the-fly, avoiding the memory allocation and computation cost of the
 * full solution.
 *
 * The computed values are returned via the output parameters.
 */
void
ncm_sbessel_ode_solver_solve_endpoints (NcmSBesselOdeSolver *solver, GArray *rhs, gdouble *deriv_a, gdouble *deriv_b, gdouble *error)
{
  const glong n_cols = _ncm_sbessel_ode_solver_diagonalize (solver, rhs);

  _ncm_sbessel_ode_solver_compute_endpoints (solver, n_cols, deriv_a, deriv_b, error);
}

/**
 * _ncm_sbessel_ode_solver_diagonalize_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @lmin: minimum l value
 * @n_l: number of l values to solve for (l = lmin, lmin+1, ..., lmin+n_l-1)
 *
 * Diagonalizes the operator using adaptive QR decomposition for multiple l values.
 * This function applies Givens rotations to transform the system into upper triangular form
 * and applies the same rotations to the RHS vectors. The transformed RHS is stored in
 * self->c_batched and the upper triangular matrix is stored in self->matrix_rows_batched.
 *
 * Returns: the effective number of columns used (may be less than rhs_len due to convergence)
 */
static glong
_ncm_sbessel_ode_solver_diagonalize_batched (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, guint n_l)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const guint rhs_len                     = rhs->len;
  guint solution_order                    = rhs_len * 2;
  const guint max_solution_order          = 1 << 24;
  guint total_rows                        = solution_order * n_l;
  const gdouble *rhs_data                 = (gdouble *) rhs->data;
  glong col                               = 0;
  gdouble max_c_A                         = 0.0;
  guint quiet_cols                        = 0;
  glong i, last_row_index;
  guint l_idx;

  g_assert_cmpuint (n_l, >, 0);

  /* Resize arrays only if necessary */
  _ensure_matrix_rows_batched_capacity (self, total_rows);

  if (total_rows + ROWS_TO_ROTATE * n_l > self->c_batched->len)
    g_array_set_size (self->c_batched, total_rows + ROWS_TO_ROTATE * n_l);

  /* Add boundary condition rows for all l values */
  for (i = 0; i < ROWS_TO_ROTATE + 1; i++)
  {
    NcmSBesselOdeSolverRow *row = &self->matrix_rows_batched[i * n_l];

    _ncm_sbessel_create_row_batched (solver, row, i, lmin, n_l);

    #pragma omp simd

    for (l_idx = 0; l_idx < n_l; l_idx++)
    {
      const glong row_idx = i * n_l + l_idx;

      g_array_index (self->c_batched, gdouble, row_idx) = rhs_data[i];
    }
  }

  last_row_index = ROWS_TO_ROTATE;

  for (col = 0; col < solution_order; col++)
  {
    const glong last_row_for_col = GSL_MIN (ROWS_TO_ROTATE, solution_order - col - 1);

    for (i = last_row_for_col; i > 0; i--)
    {
      const glong r1_index = col + i - 1;
      const glong r2_index = col + i;

      /* Create new rows if needed */
      if (r2_index > last_row_index)
      {
        NcmSBesselOdeSolverRow *row = &self->matrix_rows_batched[r2_index * n_l];

        _ncm_sbessel_create_row_batched (solver, row, r2_index, lmin, n_l);

        if (r2_index < rhs_len)
        {
          #pragma omp simd

          for (l_idx = 0; l_idx < n_l; l_idx++)
          {
            const glong row_idx = r2_index * n_l + l_idx;

            g_array_index (self->c_batched, gdouble, row_idx) = rhs_data[r2_index];
          }
        }
        else
        {
          #pragma omp simd

          for (l_idx = 0; l_idx < n_l; l_idx++)
          {
            const glong row_idx = r2_index * n_l + l_idx;

            g_array_index (self->c_batched, gdouble, row_idx) = 0.0;
          }
        }

        last_row_index = r2_index;
      }

      /* Apply Givens rotations for all l values */
      #pragma omp simd

      for (l_idx = 0; l_idx < n_l; l_idx++)
      {
        const glong r1_idx_batch   = r1_index * n_l + l_idx;
        const glong r2_idx_batch   = r2_index * n_l + l_idx;
        NcmSBesselOdeSolverRow *r1 = &self->matrix_rows_batched[r1_idx_batch];
        NcmSBesselOdeSolverRow *r2 = &self->matrix_rows_batched[r2_idx_batch];
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
        NcmSBesselOdeSolverRow *row = &self->matrix_rows_batched[row_idx];
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
        if (max_Acol < max_c_A * self->tolerance)
          quiet_cols++;
        else
          quiet_cols = 0;
      }

      if (quiet_cols >= ROWS_TO_ROTATE + 1)
        break;  /* SAFE early stop */

      /* If we've reached the end without convergence, increase solution_order and continue */
      if ((col == solution_order - ROWS_TO_ROTATE - 1) && (solution_order < max_solution_order))
      {
        solution_order *= 2;
        total_rows      = solution_order * n_l;
        _ensure_matrix_rows_batched_capacity (self, total_rows);

        if (total_rows + ROWS_TO_ROTATE * n_l > self->c_batched->len)
          g_array_set_size (self->c_batched, total_rows + ROWS_TO_ROTATE * n_l);
      }
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
 * _ncm_sbessel_ode_solver_build_solution_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @n_cols: number of columns in the solution (from diagonalization)
 * @n_l: number of l values
 * @solutions: output array of solution matrices
 *
 * Builds the full Chebyshev coefficient solution by back-substitution on the
 * upper triangular system. Assumes _ncm_sbessel_ode_solver_diagonalize_batched
 * has been called first.
 *
 * Returns: (transfer full): solution matrix where each row is the solution for one l value
 */
static void
_ncm_sbessel_ode_solver_build_solution_batched (NcmSBesselOdeSolver *solver, glong n_cols, guint n_l, GArray *solutions)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  gdouble * restrict sol_data;
  glong row;
  guint l_idx;

  /* Ensure accumulator arrays have sufficient capacity */
  _ensure_acc_bc_at_m1_capacity (self, n_l);
  _ensure_acc_bc_at_p1_capacity (self, n_l);

  /* Zero out accumulators */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_l; l_idx++)
  {
    self->acc_bc_at_m1[l_idx] = 0.0;
    self->acc_bc_at_p1[l_idx] = 0.0;
  }

  g_array_set_size (solutions, (n_cols + TOTAL_BANDWIDTH) * n_l);
  sol_data = (gdouble *) solutions->data;
  memset (sol_data, 0, (n_cols + TOTAL_BANDWIDTH) * n_l * sizeof (gdouble));

  /* Flipped loop order for better cache locality - process all l values for each row */
  for (row = n_cols - 1; row >= 0; row--)
  {
    const glong row_base_idx = row * n_l;
    const gdouble row_sign   = (row % 2) == 0 ? 1.0 : -1.0;

    #pragma omp simd

    for (l_idx = 0; l_idx < n_l; l_idx++)
    {
      const glong row_idx                       = row_base_idx + l_idx;
      const NcmSBesselOdeSolverRow * restrict r = &self->matrix_rows_batched[row_idx];
      const gdouble * restrict r_data           = r->data;
      const glong base_sol_idx                  = l_idx * n_cols;
      const gdouble * restrict sol_row          = &sol_data[base_sol_idx];
      gdouble sum                               = g_array_index (self->c_batched, gdouble, row_idx);
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

      sum -= self->acc_bc_at_m1[l_idx] * r->bc_at_m1;
      sum -= self->acc_bc_at_p1[l_idx] * r->bc_at_p1;

      sol                          = sum / diag;
      sol_data[base_sol_idx + row] = sol;

      self->acc_bc_at_m1[l_idx] += row_sign * sol;
      self->acc_bc_at_p1[l_idx] += sol;
    }
  }

  g_array_set_size (solutions, (n_cols) * n_l);
}

/**
 * _ncm_sbessel_ode_solver_compute_endpoints_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @n_cols: number of columns in the solution (from diagonalization)
 * @n_l: number of l values
 * @endpoints: (out): matrix with 3 columns per l: [y'(a), y'(b), error]
 *
 * Computes endpoint derivatives y'(a) and y'(b) and error estimates directly from
 * the diagonalized system without building the full solution matrix. This is much more
 * efficient when only endpoint information is needed, as it computes coefficients
 * on-the-fly during back-substitution and accumulates their contributions to the
 * derivatives without storing the full coefficient array.
 *
 * This function writes the results into the provided @endpoints matrix, which should
 * have at least n_l rows and exactly 3 columns. Each row corresponds to one l value,
 * with columns for y'(a), y'(b), and error estimate.
 *
 * Assumes _ncm_sbessel_ode_solver_diagonalize_batched has been called first.
 *
 * Returns: (transfer full): matrix with 3 columns per l: [y'(a), y'(b), error]
 */
static void
_ncm_sbessel_ode_solver_compute_endpoints_batched (NcmSBesselOdeSolver *solver, glong n_cols, guint n_l, GArray *endpoints)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);
  const gdouble h                         = self->half_len;
  gdouble * restrict endp_data;
  glong row;
  guint l_idx;

  /* Circular buffer for last TOTAL_BANDWIDTH coefficients - much smaller than full solution */
  _ensure_solution_batched_capacity (self, TOTAL_BANDWIDTH * n_l);
  memset (self->solution_batched, 0, sizeof (gdouble) * TOTAL_BANDWIDTH * n_l);

  g_array_set_size (endpoints, n_l * 3);

  /* Result matrix: n_l rows x 3 columns [y'(a), y'(b), error] */
  endp_data = (gdouble *) endpoints->data;

  /* Initialize endpoint derivative accumulators and error */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_l; l_idx++)
  {
    endp_data[l_idx * 3 + 0] = 0.0; /* y'(-1) accumulator */
    endp_data[l_idx * 3 + 1] = 0.0; /* y'(+1) accumulator */
    endp_data[l_idx * 3 + 2] = 0.0; /* error accumulator */
  }

  /* Ensure accumulator arrays have sufficient capacity */
  _ensure_acc_bc_at_m1_capacity (self, n_l);
  _ensure_acc_bc_at_p1_capacity (self, n_l);

  /* Zero out boundary condition accumulators */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_l; l_idx++)
  {
    self->acc_bc_at_m1[l_idx] = 0.0;
    self->acc_bc_at_p1[l_idx] = 0.0;
  }

  /* Back-substitution: compute coefficients and accumulate derivative contributions */
  for (row = n_cols - 1; row >= 0; row--)
  {
    const glong row_base_idx = row * n_l;
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

    for (l_idx = 0; l_idx < n_l; l_idx++)
    {
      const glong row_idx                       = row_base_idx + l_idx;
      const NcmSBesselOdeSolverRow * restrict r = &self->matrix_rows_batched[row_idx];
      const gdouble * restrict r_data           = r->data;
      const glong base_buffer_idx               = l_idx * TOTAL_BANDWIDTH;
      const gdouble * restrict sol_buf          = &self->solution_batched[base_buffer_idx];
      gdouble sum                               = g_array_index (self->c_batched, gdouble, row_idx);
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
      sum -= self->acc_bc_at_m1[l_idx] * r->bc_at_m1;
      sum -= self->acc_bc_at_p1[l_idx] * r->bc_at_p1;

      c_k = sum / diag;

      /* Store coefficient in circular buffer for future back-substitution steps */
      self->solution_batched[base_buffer_idx + buffer_pos] = c_k;

      /* Update boundary condition accumulators for next iteration */
      self->acc_bc_at_m1[l_idx] += row_sign * c_k;
      self->acc_bc_at_p1[l_idx] += c_k;

      /* Accumulate derivative contributions:
       * dy/dt|_{t=-1} = sum_k k^2 * (-1)^(k+1) * c_k
       * dy/dt|_{t=+1} = sum_k k^2 * c_k
       */
      endp_data[l_idx * 3 + 0] += k_squared * (-row_sign) * c_k; /* y'(-1) */
      endp_data[l_idx * 3 + 1] += k_squared * c_k;               /* y'(+1) */
      endp_data[l_idx * 3 + 2] += k_squared * fabs (c_k);        /* error estimate */
    }
  }

  /* Convert from t-derivatives to x-derivatives and finalize error */
  #pragma omp simd

  for (l_idx = 0; l_idx < n_l; l_idx++)
  {
    endp_data[l_idx * 3 + 0] /= h; /* y'(a) = y'(-1) / h */
    endp_data[l_idx * 3 + 1] /= h; /* y'(b) = y'(+1) / h */
    endp_data[l_idx * 3 + 2] /= h; /* error estimate */
  }
}

/**
 * _ncm_sbessel_ode_solver_solve_batched_internal:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: right-hand side vector (Chebyshev coefficients of f(x))
 * @lmin: minimum l value
 * @n_l: number of l values to solve for (l = lmin, lmin+1, ..., lmin+n_l-1)
 * @solutions: array of solution matrices, one per l value
 *
 * Internal batched solver implementation. Can be specialized at compile time
 * when n_l is known at compile time for better optimization.
 *
 * This function uses the factored implementation: first diagonalizes the operator,
 * then builds the full solution.
 *
 * Returns: (transfer full): solution matrix where each row is the solution for one l value
 */
void
_ncm_sbessel_ode_solver_solve_batched_internal (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, guint n_l, GArray *solutions)
{
  /* Step 1: Diagonalize the operator using QR decomposition */
  const glong n_cols = _ncm_sbessel_ode_solver_diagonalize_batched (solver, rhs, lmin, n_l);

  /* Step 2: Build the full solution by back-substitution */
  _ncm_sbessel_ode_solver_build_solution_batched (solver, n_cols, n_l, solutions);
}

/* Specialized batched solvers for common sizes - enables better compiler optimizations */

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_solver_solve_batched_2 (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, GArray *solutions)
{
  _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 2, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_solver_solve_batched_4 (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, GArray *solutions)
{
  _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 4, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_solver_solve_batched_8 (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, GArray *solutions)
{
  _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 8, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_solver_solve_batched_16 (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, GArray *solutions)
{
  _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 16, solutions);
}

static inline __attribute__ ((always_inline)) void

_ncm_sbessel_ode_solver_solve_batched_32 (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, GArray *solutions)
{
  _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, 32, solutions);
}

/**
 * ncm_sbessel_ode_solver_solve_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @lmin: minimum l value
 * @n_l: number of l values to solve for (l = lmin, lmin+1, ..., lmin+n_l-1)
 * @solutions: (out callee-allocates) (transfer full) (element-type gdouble): solution matrix where each row is the solution for one l value
 *
 * Solves the ODE for multiple l values simultaneously using batched operations.
 * Allocates rhs_len*n_l matrix rows and processes them in batches for efficiency.
 * Dispatches to specialized implementations for common sizes (8, 16, 32) for better
 * compiler optimizations including loop unrolling and vectorization.
 *
 */
void
ncm_sbessel_ode_solver_solve_batched (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, guint n_l, GArray **solutions)
{
  /* Dispatch to specialized versions for common sizes */
  if (*solutions == NULL)
    *solutions = g_array_new (FALSE, FALSE, sizeof (gdouble));

  switch (n_l)
  {
    case 2:
      _ncm_sbessel_ode_solver_solve_batched_2 (solver, rhs, lmin, *solutions);
      break;
    case 4:
      _ncm_sbessel_ode_solver_solve_batched_4 (solver, rhs, lmin, *solutions);
      break;
    case 8:
      _ncm_sbessel_ode_solver_solve_batched_8 (solver, rhs, lmin, *solutions);
      break;
    case 16:
      _ncm_sbessel_ode_solver_solve_batched_16 (solver, rhs, lmin, *solutions);
      break;
    case 32:
      _ncm_sbessel_ode_solver_solve_batched_32 (solver, rhs, lmin, *solutions);
      break;
    default:
      _ncm_sbessel_ode_solver_solve_batched_internal (solver, rhs, lmin, n_l, *solutions);
      break;
  }
}

static void
_ncm_sbessel_ode_solver_compute_endpoints_batched_4 (NcmSBesselOdeSolver *solver, glong n_cols, GArray *endpoints)
{
  _ncm_sbessel_ode_solver_compute_endpoints_batched (solver, n_cols, 4, endpoints);
}

static void
_ncm_sbessel_ode_solver_compute_endpoints_batched_8 (NcmSBesselOdeSolver *solver, glong n_cols, GArray *endpoints)
{
  _ncm_sbessel_ode_solver_compute_endpoints_batched (solver, n_cols, 8, endpoints);
}

static void
_ncm_sbessel_ode_solver_compute_endpoints_batched_16 (NcmSBesselOdeSolver *solver, glong n_cols, GArray *endpoints)
{
  _ncm_sbessel_ode_solver_compute_endpoints_batched (solver, n_cols, 16, endpoints);
}

static void
_ncm_sbessel_ode_solver_compute_endpoints_batched_32 (NcmSBesselOdeSolver *solver, glong n_cols, GArray *endpoints)
{
  _ncm_sbessel_ode_solver_compute_endpoints_batched (solver, n_cols, 32, endpoints);
}

/**
 * ncm_sbessel_ode_solver_solve_endpoints_batched:
 * @solver: a #NcmSBesselOdeSolver
 * @rhs: (element-type gdouble): right-hand side vector (Chebyshev coefficients of f(x))
 * @lmin: minimum l value
 * @n_l: number of l values to solve for (l = lmin, lmin+1, ..., lmin+n_l-1)
 * @endpoints: (out callee-allocates) (transfer full) (element-type gdouble): a matrix (n_l x 3) to store [y'(a), y'(b), error] for each l
 *
 * Efficiently computes only the endpoint derivatives y'(a) and y'(b) for multiple l values
 * without building the full solution matrix. This is much more efficient when only endpoint
 * information is needed (e.g., for integral computations via Green's identity).
 *
 * The result matrix must be pre-allocated with at least n_l rows and exactly 3
 * columns. Each row corresponds to one l value, with the first column storing y'(a),
 * the second column storing y'(b), and the third column storing the error estimate.
 *
 * The function first diagonalizes the operator using adaptive QR decomposition, then
 * performs back-substitution while accumulating the contributions to the endpoint
 * derivatives on-the-fly, avoiding the memory allocation and computation cost of the
 * full solution.
 */
void
ncm_sbessel_ode_solver_solve_endpoints_batched (NcmSBesselOdeSolver *solver, GArray *rhs, gint lmin, guint n_l, GArray **endpoints)
{
  const glong n_cols = _ncm_sbessel_ode_solver_diagonalize_batched (solver, rhs, lmin, n_l);

  if (*endpoints == NULL)
    *endpoints = g_array_new (FALSE, FALSE, sizeof (gdouble));

  switch (n_l)
  {
    case 4:
      _ncm_sbessel_ode_solver_compute_endpoints_batched_4 (solver, n_cols, *endpoints);
      break;
    case 8:
      _ncm_sbessel_ode_solver_compute_endpoints_batched_8 (solver, n_cols, *endpoints);
      break;
    case 16:
      _ncm_sbessel_ode_solver_compute_endpoints_batched_16 (solver, n_cols, *endpoints);
      break;
    case 32:
      _ncm_sbessel_ode_solver_compute_endpoints_batched_32 (solver, n_cols, *endpoints);
      break;
    default:
      _ncm_sbessel_ode_solver_compute_endpoints_batched (solver, n_cols, n_l, *endpoints);
      break;
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
  NcmSBesselOdeSolverRow row_mem;
  NcmSBesselOdeSolverRow *row = &row_mem;
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
  {
    const gint ncols       = nrows;
    gdouble *data_colmajor = g_new0 (gdouble, nrows * ncols);
    NcmMatrix *mat         = ncm_matrix_new_data_malloc (data_colmajor, nrows, ncols);

    _ncm_sbessel_ode_solver_fill_operator_matrix (solver, nrows, ncols,
                                                  data_colmajor, TRUE);

    return mat;
  }
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
 * ncm_sbessel_ode_solver_peek_spectral:
 * @solver: a #NcmSBesselOdeSolver
 *
 * Returns a pointer to the internal spectral object used by the solver. This allows
 * the user to inspect the internal state of the solver. The spectral object contains
 * information about the Chebyshev nodes, weights, and Chebyshev polynomials.
 *
 * Returns: (transfer full): spectral object
 */
NcmSpectral *
ncm_sbessel_ode_solver_peek_spectral (NcmSBesselOdeSolver *solver)
{
  NcmSBesselOdeSolverPrivate * const self = ncm_sbessel_ode_solver_get_instance_private (solver);

  return self->spectral;
}


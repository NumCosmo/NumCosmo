/***************************************************************************
 *            test_ncm_sbessel_ode_solver.c
 *
 *  Wed October 14 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sbessel_ode_solver.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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

/*
 * NcmSBesselOdeSolver and NcmSBesselOdeOperator: the entry points (solve, solve_values,
 * solve_endpoints, solve_dense) against each other and against an LU solve of the
 * operator matrix, zero right-hand sides, the accessors and the reuse of a
 * factorization. The right-hand side is the endpoint data u(a), u(b) followed by the
 * C^(2) coefficients of the forcing, x in [a, b] mapped to t in [-1, 1].
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <gsl/gsl_linalg.h>

void test_ncm_sbessel_ode_solver_solve_values (void);
void test_ncm_sbessel_ode_solver_batched_convergence (void);
void test_ncm_sbessel_ode_solver_qr_vs_lu (void);
void test_ncm_sbessel_ode_solver_dense_vs_lu (void);
void test_ncm_sbessel_ode_solver_zero_rhs (void);
void test_ncm_sbessel_ode_solver_zero_rhs_dense (void);
void test_ncm_sbessel_ode_solver_zero_rhs_endpoints (void);
void test_ncm_sbessel_ode_solver_zero_rhs_batched (void);
void test_ncm_sbessel_ode_solver_batched_vs_single (void);
void test_ncm_sbessel_ode_solver_batched_vs_single_exp (void);
void test_ncm_sbessel_ode_solver_endpoints_vs_full (void);
void test_ncm_sbessel_ode_solver_endpoints_finite (void);
void test_ncm_sbessel_ode_solver_endpoints_batched_vs_single (void);
void test_ncm_sbessel_ode_solver_endpoints_batched_vs_full (void);
void test_ncm_sbessel_ode_solver_accessors (void);
void test_ncm_sbessel_ode_solver_reconfigure_tolerance (void);
void test_ncm_sbessel_ode_solver_reuse (void);
void test_ncm_sbessel_ode_solver_reuse_correctness (void);
void test_ncm_sbessel_ode_solver_reuse_endpoints (void);
void test_ncm_sbessel_ode_solver_optimized_batched (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

  g_test_add_func ("/ncm/sbessel_ode_solver/solve_values", test_ncm_sbessel_ode_solver_solve_values);
  g_test_add_func ("/ncm/sbessel_ode_solver/batched_convergence", test_ncm_sbessel_ode_solver_batched_convergence);
  g_test_add_func ("/ncm/sbessel_ode_solver/qr_vs_lu", test_ncm_sbessel_ode_solver_qr_vs_lu);
  g_test_add_func ("/ncm/sbessel_ode_solver/dense_vs_lu", test_ncm_sbessel_ode_solver_dense_vs_lu);
  g_test_add_func ("/ncm/sbessel_ode_solver/zero_rhs/solve", test_ncm_sbessel_ode_solver_zero_rhs);
  g_test_add_func ("/ncm/sbessel_ode_solver/zero_rhs/dense", test_ncm_sbessel_ode_solver_zero_rhs_dense);
  g_test_add_func ("/ncm/sbessel_ode_solver/zero_rhs/endpoints", test_ncm_sbessel_ode_solver_zero_rhs_endpoints);
  g_test_add_func ("/ncm/sbessel_ode_solver/zero_rhs/batched", test_ncm_sbessel_ode_solver_zero_rhs_batched);
  g_test_add_func ("/ncm/sbessel_ode_solver/batched_vs_single", test_ncm_sbessel_ode_solver_batched_vs_single);
  g_test_add_func ("/ncm/sbessel_ode_solver/batched_vs_single_exp", test_ncm_sbessel_ode_solver_batched_vs_single_exp);
  g_test_add_func ("/ncm/sbessel_ode_solver/endpoints/vs_full", test_ncm_sbessel_ode_solver_endpoints_vs_full);
  g_test_add_func ("/ncm/sbessel_ode_solver/endpoints/finite", test_ncm_sbessel_ode_solver_endpoints_finite);
  g_test_add_func ("/ncm/sbessel_ode_solver/endpoints/batched_vs_single", test_ncm_sbessel_ode_solver_endpoints_batched_vs_single);
  g_test_add_func ("/ncm/sbessel_ode_solver/endpoints/batched_vs_full", test_ncm_sbessel_ode_solver_endpoints_batched_vs_full);
  g_test_add_func ("/ncm/sbessel_ode_solver/accessors", test_ncm_sbessel_ode_solver_accessors);
  g_test_add_func ("/ncm/sbessel_ode_solver/reconfigure_tolerance", test_ncm_sbessel_ode_solver_reconfigure_tolerance);
  g_test_add_func ("/ncm/sbessel_ode_solver/reuse/factorization", test_ncm_sbessel_ode_solver_reuse);
  g_test_add_func ("/ncm/sbessel_ode_solver/reuse/correctness", test_ncm_sbessel_ode_solver_reuse_correctness);
  g_test_add_func ("/ncm/sbessel_ode_solver/reuse/endpoints", test_ncm_sbessel_ode_solver_reuse_endpoints);
  g_test_add_func ("/ncm/sbessel_ode_solver/optimized_batched", test_ncm_sbessel_ode_solver_optimized_batched);

  g_test_run ();
}

#define TEST_A 1.0
#define TEST_B 20.0
#define TEST_H (0.5 * (TEST_B - TEST_A))
#define TEST_M (0.5 * (TEST_A + TEST_B))

/* |actual - desired| <= atol + rtol |desired|, the numpy assert_allclose criterion */
static void
_test_assert_close (gdouble actual, gdouble desired, gdouble rtol, gdouble atol)
{
  g_assert_true (gsl_finite (actual));
  g_assert_cmpfloat (fabs (actual - desired), <=, atol + rtol * fabs (desired));
}

/* A right-hand side of length n with zero endpoint data and the given leading coefficients */
static GArray *
_test_rhs_new (guint n, const gdouble *coeffs, guint n_coeffs)
{
  GArray *rhs = g_array_sized_new (FALSE, TRUE, sizeof (gdouble), n);
  guint i;

  g_array_set_size (rhs, n);

  for (i = 0; i < n_coeffs; i++)
    g_array_index (rhs, gdouble, 2 + i) = coeffs[i];

  return rhs;
}

/* Constant forcing: rhs[2] = c */
static GArray *
_test_rhs_const (guint n, gdouble c)
{
  return _test_rhs_new (n, &c, 1);
}

/* rhs[k] = 1 / (1 + k^2) for k >= 2 */
static GArray *
_test_rhs_decaying (guint n)
{
  GArray *rhs = _test_rhs_new (n, NULL, 0);
  guint k;

  for (k = 2; k < n; k++)
    g_array_index (rhs, gdouble, k) = 1.0 / (1.0 + k * k);

  return rhs;
}

typedef gdouble (*TestFx) (gdouble x);

static gdouble
_test_cheb_func (gpointer user_data, gdouble t)
{
  TestFx f = (TestFx) user_data;

  return f (TEST_M + TEST_H * t);
}

/* f (x) on [1, 20] as Chebyshev coefficients of order n - 2 converted to C^(2), after zero endpoint data */
static GArray *
_test_rhs_func (guint n, TestFx f)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *cheb          = NULL;
  GArray *gegen         = NULL;
  GArray *rhs;

  ncm_spectral_compute_chebyshev_coeffs (spectral, &_test_cheb_func, -1.0, 1.0, n - 2, &cheb, (gpointer) f);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (cheb, &gegen);
  rhs = _test_rhs_new (n, (const gdouble *) gegen->data, GSL_MIN (gegen->len, n - 2));

  g_array_unref (cheb);
  g_array_unref (gegen);
  ncm_spectral_free (spectral);

  return rhs;
}

/* Row @i of a batched solution of length @len per multipole, as a new array */
static GArray *
_test_row (GArray *sol, gsize len, guint i)
{
  GArray *row = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);

  g_array_append_vals (row, &g_array_index (sol, gdouble, i * len), len);

  return row;
}

static gdouble
_test_absmax (GArray *a)
{
  NcmVector *v = ncm_vector_new_array (a);
  gdouble absmin, absmax;

  ncm_vector_get_absminmax (v, &absmin, &absmax);
  ncm_vector_free (v);

  return absmax;
}

/* An operator on [a, b] for ell_min..ell_max, solved for @rhs */
static NcmSBesselOdeOperator *
_test_solve (NcmSBesselOdeSolver *solver, gdouble a, gdouble b, gint ell_min, gint ell_max, GArray *rhs, GArray **sol, gsize *len)
{
  NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, a, b, ell_min, ell_max);

  ncm_sbessel_ode_operator_solve (op, rhs, sol, len);

  return op;
}

/* (u(x0), u'(x0), u(x1), u'(x1)) from Chebyshev coefficients on [a, b] */
static void
_test_values (GArray *c, gdouble a, gdouble b, gdouble x0, gdouble x1, gdouble v[4])
{
  v[0] = ncm_spectral_chebyshev_eval (c, (2.0 * x0 - a - b) / (b - a));
  v[1] = ncm_spectral_chebyshev_deriv_x (c, a, b, x0);
  v[2] = ncm_spectral_chebyshev_eval (c, (2.0 * x1 - a - b) / (b - a));
  v[3] = ncm_spectral_chebyshev_deriv_x (c, a, b, x1);
}

/* solve_values against the Chebyshev evaluation of solve(), and against solve_endpoints at (a, b) */
void
test_ncm_sbessel_ode_solver_solve_values (void)
{
  const guint n_ell_a[] = {1, 8};
  const gdouble f[]     = {1.0, -0.2, 0.05, 0.03, -0.01};
  const gdouble xs[][2] = {
    {
      3.7, 14.2
    }, {
      TEST_A, TEST_B
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_new (65, f, G_N_ELEMENTS (f));
  guint c, i, j, l;

  for (c = 0; c < G_N_ELEMENTS (n_ell_a); c++)
  {
    const guint n_ell = n_ell_a[c];
    GArray *sol       = NULL;
    GArray *values    = NULL;
    GArray *endpoints = NULL;
    gsize len;
    NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, 3, 3 + n_ell - 1, rhs, &sol, &len);

    for (i = 0; i < G_N_ELEMENTS (xs); i++)
    {
      ncm_sbessel_ode_operator_solve_values (op, rhs, xs[i][0], xs[i][1], &values);
      g_assert_cmpuint (values->len, ==, 4 * n_ell);

      for (l = 0; l < n_ell; l++)
      {
        GArray *row = _test_row (sol, len, l);
        gdouble v[4];

        _test_values (row, TEST_A, TEST_B, xs[i][0], xs[i][1], v);

        for (j = 0; j < 4; j++)
          _test_assert_close (g_array_index (values, gdouble, 4 * l + j), v[j], 2.0e-11, 2.0e-12);

        g_array_unref (row);
      }
    }

    ncm_sbessel_ode_operator_solve_values (op, rhs, TEST_A, TEST_B, &values);
    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);

    for (l = 0; l < n_ell; l++)
    {
      _test_assert_close (g_array_index (values, gdouble, 4 * l + 1), g_array_index (endpoints, gdouble, 3 * l + 0), 2.0e-11, 0.0);
      _test_assert_close (g_array_index (values, gdouble, 4 * l + 3), g_array_index (endpoints, gdouble, 3 * l + 1), 2.0e-11, 0.0);
    }

    g_array_unref (sol);
    g_array_unref (values);
    g_array_unref (endpoints);
    ncm_sbessel_ode_operator_unref (op);
  }

  g_array_unref (rhs);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * A batch of ell 0..n_ell - 1 on [0.1, 10] at tolerance 1e-10 resolves every multipole:
 * each row agrees with its single-multipole solve, and solve_values and solve_endpoints
 * agree with the Chebyshev evaluation of the batched solution.
 */
void
test_ncm_sbessel_ode_solver_batched_convergence (void)
{
  const guint n_ell_a[]       = {1, 2, 3, 4, 7, 8, 16, 32, 64};
  const gdouble a             = 0.1;
  const gdouble b             = 10.0;
  const gdouble x0            = 1.3;
  const gdouble x1            = 8.7;
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_const (65, 1.0);
  guint c, l, j;

  ncm_sbessel_ode_solver_set_tolerance (solver, 1.0e-10);

  for (c = 0; c < G_N_ELEMENTS (n_ell_a); c++)
  {
    const guint n_ell = n_ell_a[c];
    GArray *sol       = NULL;
    GArray *values    = NULL;
    GArray *endpoints = NULL;
    gsize len;
    NcmSBesselOdeOperator *op = _test_solve (solver, a, b, 0, n_ell - 1, rhs, &sol, &len);

    for (l = 0; l < n_ell; l++)
    {
      GArray *single = NULL;
      gsize single_len;
      NcmSBesselOdeOperator *sop = _test_solve (solver, a, b, l, l, rhs, &single, &single_len);
      const gdouble scale        = _test_absmax (single);

      g_assert_cmpuint (len, >=, single_len);

      for (j = 0; j < len; j++)
      {
        const gdouble s = (j < single_len) ? g_array_index (single, gdouble, j) : 0.0;

        _test_assert_close (g_array_index (sol, gdouble, l * len + j), s, 5.0e-8, 5.0e-8 * scale);
      }

      g_array_unref (single);
      ncm_sbessel_ode_operator_unref (sop);
    }

    ncm_sbessel_ode_operator_solve_values (op, rhs, x0, x1, &values);
    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);

    for (l = 0; l < n_ell; l++)
    {
      GArray *row = _test_row (sol, len, l);
      gdouble v[4];

      _test_values (row, a, b, x0, x1, v);

      for (j = 0; j < 4; j++)
        _test_assert_close (g_array_index (values, gdouble, 4 * l + j), v[j], 2.0e-11, 2.0e-13);

      _test_assert_close (g_array_index (endpoints, gdouble, 3 * l + 0), ncm_spectral_chebyshev_deriv_x (row, a, b, a), 2.0e-11, 2.0e-13);
      _test_assert_close (g_array_index (endpoints, gdouble, 3 * l + 1), ncm_spectral_chebyshev_deriv_x (row, a, b, b), 2.0e-11, 2.0e-13);

      g_array_unref (row);
    }

    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, (glong) len);

    g_array_unref (sol);
    g_array_unref (values);
    g_array_unref (endpoints);
    ncm_sbessel_ode_operator_unref (op);
  }

  g_array_unref (rhs);
  ncm_sbessel_ode_solver_free (solver);
}

/* The solution of the n x n operator matrix for @rhs by LU (NcmMatrix has no general LU, so GSL's on its storage) */
static NcmVector *
_test_lu_solve (NcmSBesselOdeSolver *solver, guint ell, guint n, GArray *rhs)
{
  NcmMatrix *mat     = ncm_sbessel_ode_solver_get_operator_matrix (solver, TEST_A, TEST_B, ell, n);
  NcmVector *b       = ncm_vector_new_array (rhs);
  NcmVector *x       = ncm_vector_new (n);
  gsl_permutation *p = gsl_permutation_alloc (n);
  gint signum;

  g_assert_cmpuint (ncm_matrix_nrows (mat), ==, n);
  g_assert_cmpuint (ncm_matrix_ncols (mat), ==, n);
  g_assert_cmpuint (ncm_vector_len (b), ==, n);

  g_assert_cmpint (gsl_linalg_LU_decomp (ncm_matrix_gsl (mat), p, &signum), ==, GSL_SUCCESS);
  g_assert_cmpint (gsl_linalg_LU_solve (ncm_matrix_gsl (mat), p, ncm_vector_const_gsl (b), ncm_vector_gsl (x)), ==, GSL_SUCCESS);

  gsl_permutation_free (p);
  ncm_vector_free (b);
  ncm_matrix_free (mat);

  return x;
}

/*
 * The adaptive QR solve against an LU solve of the 2n x 2n operator matrix, ell 0..20 and
 * n = 32, 64, 128: sum (x_LU - x_QR)^2 / sum x_QR^2 < 1e-10 over the common length.
 */
void
test_ncm_sbessel_ode_solver_qr_vs_lu (void)
{
  const guint n_a[]           = {32, 64, 128};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint ell, c, j;

  for (ell = 0; ell <= 20; ell++)
  {
    for (c = 0; c < G_N_ELEMENTS (n_a); c++)
    {
      GArray *rhs  = _test_rhs_const (n_a[c], 1.0);
      GArray *rhs2 = _test_rhs_const (2 * n_a[c], 1.0);
      GArray *sol  = NULL;
      gsize len;
      NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, ell, ell, rhs, &sol, &len);
      NcmVector *lu             = _test_lu_solve (solver, ell, 2 * n_a[c], rhs2);
      const guint safe          = GSL_MIN (len, ncm_vector_len (lu));
      gdouble num               = 0.0;
      gdouble den               = 0.0;

      g_assert_cmpuint (len, >, 0);

      for (j = 0; j < safe; j++)
      {
        num += gsl_pow_2 (ncm_vector_get (lu, j) - g_array_index (sol, gdouble, j));
        den += gsl_pow_2 (g_array_index (sol, gdouble, j));
      }

      g_assert_cmpfloat (num / den, <, 1.0e-10);

      g_array_unref (rhs);
      g_array_unref (rhs2);
      g_array_unref (sol);
      ncm_vector_free (lu);
      ncm_sbessel_ode_operator_unref (op);
    }
  }

  ncm_sbessel_ode_solver_free (solver);
}

/* solve_dense (LAPACK dgesv on the column-major matrix) against GSL's LU, ell 0..20, n = 32, 64, 128 */
void
test_ncm_sbessel_ode_solver_dense_vs_lu (void)
{
  const guint n_a[]           = {32, 64, 128};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint ell, c, j;

  for (ell = 0; ell <= 20; ell++)
  {
    for (c = 0; c < G_N_ELEMENTS (n_a); c++)
    {
      GArray *rhs      = _test_rhs_const (n_a[c], 1.0);
      NcmVector *rhs_v = ncm_vector_new_array (rhs);
      NcmVector *lu    = _test_lu_solve (solver, ell, n_a[c], rhs);
      NcmVector *dense = ncm_sbessel_ode_solver_solve_dense (solver, TEST_A, TEST_B, ell, rhs_v, n_a[c]);

      g_assert_cmpuint (ncm_vector_len (dense), ==, n_a[c]);

      for (j = 0; j < n_a[c]; j++)
        _test_assert_close (ncm_vector_get (dense, j), ncm_vector_get (lu, j), 1.0e-10, 1.0e-12);

      g_array_unref (rhs);
      ncm_vector_free (lu);
      ncm_vector_free (rhs_v);
      ncm_vector_free (dense);
    }
  }

  ncm_sbessel_ode_solver_free (solver);
}

static const guint test_zero_ells[] = {0, 1, 5, 10, 20};
static const guint test_zero_ns[]   = {32, 64, 128};

/* A zero right-hand side, endpoint data included, gives a zero solution no longer than the input */
void
test_ncm_sbessel_ode_solver_zero_rhs (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint i, c, j;

  for (i = 0; i < G_N_ELEMENTS (test_zero_ells); i++)
  {
    for (c = 0; c < G_N_ELEMENTS (test_zero_ns); c++)
    {
      GArray *rhs = _test_rhs_new (test_zero_ns[c], NULL, 0);
      GArray *sol = NULL;
      gsize len;
      NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, test_zero_ells[i], test_zero_ells[i], rhs, &sol, &len);

      g_assert_cmpuint (len, <=, test_zero_ns[c]);

      for (j = 0; j < len; j++)
        g_assert_cmpfloat (fabs (g_array_index (sol, gdouble, j)), <=, 1.0e-15);

      g_array_unref (rhs);
      g_array_unref (sol);
      ncm_sbessel_ode_operator_unref (op);
    }
  }

  ncm_sbessel_ode_solver_free (solver);
}

void
test_ncm_sbessel_ode_solver_zero_rhs_dense (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint i, c, j;

  for (i = 0; i < G_N_ELEMENTS (test_zero_ells); i++)
  {
    for (c = 0; c < G_N_ELEMENTS (test_zero_ns); c++)
    {
      NcmVector *rhs = ncm_vector_new (test_zero_ns[c]);
      NcmVector *sol;

      ncm_vector_set_zero (rhs);
      sol = ncm_sbessel_ode_solver_solve_dense (solver, TEST_A, TEST_B, test_zero_ells[i], rhs, test_zero_ns[c]);

      g_assert_cmpuint (ncm_vector_len (sol), ==, test_zero_ns[c]);

      for (j = 0; j < test_zero_ns[c]; j++)
        g_assert_cmpfloat (fabs (ncm_vector_get (sol, j)), <=, 1.0e-15);

      ncm_vector_free (rhs);
      ncm_vector_free (sol);
    }
  }

  ncm_sbessel_ode_solver_free (solver);
}

void
test_ncm_sbessel_ode_solver_zero_rhs_endpoints (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint i, c;

  for (i = 0; i < G_N_ELEMENTS (test_zero_ells); i++)
  {
    for (c = 0; c < G_N_ELEMENTS (test_zero_ns); c++)
    {
      GArray *rhs               = _test_rhs_new (test_zero_ns[c], NULL, 0);
      GArray *endpoints         = NULL;
      NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, test_zero_ells[i], test_zero_ells[i]);

      ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);
      g_assert_cmpfloat (fabs (g_array_index (endpoints, gdouble, 0)), <=, 1.0e-15);
      g_assert_cmpfloat (fabs (g_array_index (endpoints, gdouble, 1)), <=, 1.0e-15);

      g_array_unref (rhs);
      g_array_unref (endpoints);
      ncm_sbessel_ode_operator_unref (op);
    }
  }

  ncm_sbessel_ode_solver_free (solver);
}

/* The same for a batch of ell 5..15, through solve and solve_endpoints */
void
test_ncm_sbessel_ode_solver_zero_rhs_batched (void)
{
  const guint n_ell           = 11;
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_new (128, NULL, 0);
  GArray *sol                 = NULL;
  GArray *endpoints           = NULL;
  gsize len;
  NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, 5, 15, rhs, &sol, &len);
  NcmSBesselOdeOperator *ep = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 15);
  guint j;

  g_assert_cmpuint (len, <=, 128);

  for (j = 0; j < n_ell * len; j++)
    g_assert_cmpfloat (fabs (g_array_index (sol, gdouble, j)), <=, 1.0e-15);

  ncm_sbessel_ode_operator_solve_endpoints (ep, rhs, &endpoints);

  for (j = 0; j < n_ell; j++)
  {
    g_assert_cmpfloat (fabs (g_array_index (endpoints, gdouble, 3 * j + 0)), <=, 1.0e-15);
    g_assert_cmpfloat (fabs (g_array_index (endpoints, gdouble, 3 * j + 1)), <=, 1.0e-15);
  }

  g_array_unref (rhs);
  g_array_unref (sol);
  g_array_unref (endpoints);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_operator_unref (ep);
  ncm_sbessel_ode_solver_free (solver);
}

/* A batch of ell 5..25 with rhs[k] = 1 / (1 + k^2) agrees with single-multipole solves on their common length */
void
test_ncm_sbessel_ode_solver_batched_vs_single (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_decaying (128);
  GArray *sol                 = NULL;
  gsize len;
  NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, 5, 25, rhs, &sol, &len);
  gint ell;
  guint j;

  for (ell = 5; ell <= 25; ell++)
  {
    GArray *single = NULL;
    gsize single_len;
    NcmSBesselOdeOperator *sop = _test_solve (solver, TEST_A, TEST_B, ell, ell, rhs, &single, &single_len);

    for (j = 0; j < GSL_MIN (len, single_len); j++)
      _test_assert_close (g_array_index (sol, gdouble, (ell - 5) * len + j), g_array_index (single, gdouble, j), 1.0e-12, 1.0e-14);

    g_array_unref (single);
    ncm_sbessel_ode_operator_unref (sop);
  }

  g_array_unref (rhs);
  g_array_unref (sol);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * A batched solve against single-multipole solves: the batch is at least as long as
 * each single solve, and each row agrees with the zero-padded single solution to
 * |d| <= rtol |s| + atol_factor max |s|.
 */
static void
_test_compare_batched_solve (NcmSBesselOdeSolver *solver, gint ell_min, gint ell_max, GArray *rhs, gdouble rtol, gdouble atol_factor)
{
  GArray *sol = NULL;
  gsize len;
  NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, ell_min, ell_max, rhs, &sol, &len);
  gint ell;
  guint j;

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    GArray *single = NULL;
    gsize single_len;
    NcmSBesselOdeOperator *sop = _test_solve (solver, TEST_A, TEST_B, ell, ell, rhs, &single, &single_len);
    const gdouble scale        = _test_absmax (single);

    g_assert_cmpuint (len, >=, single_len);

    for (j = 0; j < len; j++)
    {
      const gdouble s = (j < single_len) ? g_array_index (single, gdouble, j) : 0.0;

      _test_assert_close (g_array_index (sol, gdouble, (ell - ell_min) * len + j), s, rtol, atol_factor * scale);
    }

    g_array_unref (single);
    ncm_sbessel_ode_operator_unref (sop);
  }

  g_array_unref (sol);
  ncm_sbessel_ode_operator_unref (op);
}

/*
 * The same for solve_endpoints: (u'(a), u'(b)) and, when @with_error, the error estimate,
 * to |d| <= rtol |s| + atol_factor max (|u'(a)|, |u'(b)|).
 */
static void
_test_compare_batched_endpoints (NcmSBesselOdeSolver *solver, gint ell_min, gint ell_max, GArray *rhs, gdouble rtol, gdouble atol_factor, gboolean with_error)
{
  NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell_min, ell_max);
  GArray *batched           = NULL;
  gint ell;
  guint j;

  ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &batched);
  g_assert_cmpuint (batched->len, ==, 3 * (ell_max - ell_min + 1));

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    NcmSBesselOdeOperator *sop = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell, ell);
    GArray *single             = NULL;
    gdouble scale;

    ncm_sbessel_ode_operator_solve_endpoints (sop, rhs, &single);
    scale = GSL_MAX (fabs (g_array_index (single, gdouble, 0)), fabs (g_array_index (single, gdouble, 1)));

    for (j = 0; j < (with_error ? 3 : 2); j++)
      _test_assert_close (g_array_index (batched, gdouble, 3 * (ell - ell_min) + j), g_array_index (single, gdouble, j), rtol, atol_factor * scale);

    g_array_unref (single);
    ncm_sbessel_ode_operator_unref (sop);
  }

  g_array_unref (batched);
  ncm_sbessel_ode_operator_unref (op);
}

static gdouble
_test_exp (gdouble x)
{
  return exp (x);
}

static gdouble
_test_sin (gdouble x)
{
  return sin (x);
}

/* exp (x) forcing, ell 5..10, n = 128 */
void
test_ncm_sbessel_ode_solver_batched_vs_single_exp (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_func (128, &_test_exp);

  _test_compare_batched_solve (solver, 5, 10, rhs, 1.0e-12, 1.0e-12);

  g_array_unref (rhs);
  ncm_sbessel_ode_solver_free (solver);
}

/* solve_endpoints against the derivative of the full solve() at t = -1, +1, ell 0, 1, 5, 10, 20 */
void
test_ncm_sbessel_ode_solver_endpoints_vs_full (void)
{
  const guint ell_a[]         = {0, 1, 5, 10, 20};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_const (128, 1.0);
  guint i;

  for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
  {
    NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell_a[i], ell_a[i]);
    GArray *endpoints         = NULL;
    GArray *sol               = NULL;
    gsize len;

    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);
    ncm_sbessel_ode_operator_solve (op, rhs, &sol, &len);

    _test_assert_close (g_array_index (endpoints, gdouble, 0), ncm_spectral_chebyshev_deriv (sol, -1.0) / TEST_H, 1.0e-12, 1.0e-15);
    _test_assert_close (g_array_index (endpoints, gdouble, 1), ncm_spectral_chebyshev_deriv (sol, +1.0) / TEST_H, 1.0e-12, 1.0e-15);

    g_array_unref (endpoints);
    g_array_unref (sol);
    ncm_sbessel_ode_operator_unref (op);
  }

  g_array_unref (rhs);
  ncm_sbessel_ode_solver_free (solver);
}

static void
_test_assert_endpoints_finite (GArray *endpoints, guint n_ell)
{
  guint l;

  g_assert_cmpuint (endpoints->len, ==, 3 * n_ell);

  for (l = 0; l < n_ell; l++)
  {
    g_assert_true (gsl_finite (g_array_index (endpoints, gdouble, 3 * l + 0)));
    g_assert_true (gsl_finite (g_array_index (endpoints, gdouble, 3 * l + 1)));
    g_assert_true (gsl_finite (g_array_index (endpoints, gdouble, 3 * l + 2)));
    g_assert_cmpfloat (g_array_index (endpoints, gdouble, 3 * l + 2), >=, 0.0);
  }
}

/*
 * solve_endpoints returns finite derivatives and a non-negative error estimate: single
 * multipoles ell 0, 1, 5, 10 with rhs (m^2 + h^2 / 2, m h / 2), and a batch of ell 5..15
 * with rhs[k] = 1 / (1 + k^2).
 */
void
test_ncm_sbessel_ode_solver_endpoints_finite (void)
{
  const guint ell_a[]         = {0, 1, 5, 10};
  const gdouble f[]           = {TEST_M * TEST_M + TEST_H * TEST_H / 2.0, 0.25 * 2.0 * TEST_M * TEST_H};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_new (128, f, G_N_ELEMENTS (f));
  GArray *rhs_d               = _test_rhs_decaying (128);
  GArray *endpoints           = NULL;
  NcmSBesselOdeOperator *op;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
  {
    op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell_a[i], ell_a[i]);
    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);
    _test_assert_endpoints_finite (endpoints, 1);
    ncm_sbessel_ode_operator_unref (op);
  }

  op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 15);
  ncm_sbessel_ode_operator_solve_endpoints (op, rhs_d, &endpoints);
  _test_assert_endpoints_finite (endpoints, 11);
  ncm_sbessel_ode_operator_unref (op);

  g_array_unref (rhs);
  g_array_unref (rhs_d);
  g_array_unref (endpoints);
  ncm_sbessel_ode_solver_free (solver);
}

/* Batched solve_endpoints, ell 5..15, against single calls: all three entries */
void
test_ncm_sbessel_ode_solver_endpoints_batched_vs_single (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 15);
  GArray *rhs                 = _test_rhs_const (128, 1.0);
  GArray *batched             = NULL;
  gint ell;
  guint j;

  ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &batched);

  for (ell = 5; ell <= 15; ell++)
  {
    NcmSBesselOdeOperator *sop = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell, ell);
    GArray *single             = NULL;

    ncm_sbessel_ode_operator_solve_endpoints (sop, rhs, &single);

    for (j = 0; j < 3; j++)
      _test_assert_close (g_array_index (batched, gdouble, 3 * (ell - 5) + j), g_array_index (single, gdouble, j), 1.0e-12, 1.0e-15);

    g_array_unref (single);
    ncm_sbessel_ode_operator_unref (sop);
  }

  g_array_unref (rhs);
  g_array_unref (batched);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* Batched solve_endpoints, ell 5..15, against the derivative of the batched solve() */
void
test_ncm_sbessel_ode_solver_endpoints_batched_vs_full (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 15);
  GArray *rhs                 = _test_rhs_const (128, 1.0);
  GArray *endpoints           = NULL;
  GArray *sol                 = NULL;
  gsize len;
  guint l;

  ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);
  ncm_sbessel_ode_operator_solve (op, rhs, &sol, &len);

  for (l = 0; l < 11; l++)
  {
    GArray *row = _test_row (sol, len, l);

    _test_assert_close (g_array_index (endpoints, gdouble, 3 * l + 0), ncm_spectral_chebyshev_deriv (row, -1.0) / TEST_H, 1.0e-12, 1.0e-15);
    _test_assert_close (g_array_index (endpoints, gdouble, 3 * l + 1), ncm_spectral_chebyshev_deriv (row, +1.0) / TEST_H, 1.0e-12, 1.0e-15);

    g_array_unref (row);
  }

  g_array_unref (rhs);
  g_array_unref (endpoints);
  g_array_unref (sol);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

static void
_test_assert_accessors (NcmSBesselOdeOperator *op, gdouble a, gdouble b, gint ell_min, gint ell_max, gdouble tol)
{
  gdouble a_out, b_out;
  gint lo, hi;

  ncm_sbessel_ode_operator_get_interval (op, &a_out, &b_out);
  ncm_sbessel_ode_operator_get_ell_range (op, &lo, &hi);

  _test_assert_close (a_out, a, 1.0e-15, 1.0e-15);
  _test_assert_close (b_out, b, 1.0e-15, 1.0e-15);
  g_assert_cmpint (lo, ==, ell_min);
  g_assert_cmpint (hi, ==, ell_max);
  _test_assert_close (ncm_sbessel_ode_operator_get_tolerance (op), tol, 1.0e-15, 1.0e-15);
}

/*
 * The operator reports its interval, multipole range and tolerance (the solver's at
 * creation, default or set), and a solve leaves them unchanged.
 */
void
test_ncm_sbessel_ode_solver_accessors (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  const gdouble tol_default   = ncm_sbessel_ode_solver_get_tolerance (solver);
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 5);
  GArray *rhs                 = _test_rhs_const (64, 1.0);
  GArray *sol                 = NULL;
  gsize len;

  _test_assert_accessors (op, TEST_A, TEST_B, 5, 5, tol_default);
  ncm_sbessel_ode_operator_unref (op);

  op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 7, 7);
  _test_assert_accessors (op, TEST_A, TEST_B, 7, 7, tol_default);
  ncm_sbessel_ode_operator_unref (op);

  op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 15);
  _test_assert_accessors (op, TEST_A, TEST_B, 5, 15, tol_default);
  ncm_sbessel_ode_operator_unref (op);

  ncm_sbessel_ode_solver_set_tolerance (solver, 1.0e-10);
  op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 10, 10);
  _test_assert_accessors (op, TEST_A, TEST_B, 10, 10, 1.0e-10);
  ncm_sbessel_ode_operator_unref (op);

  ncm_sbessel_ode_solver_set_tolerance (solver, 1.0e-9);
  op = ncm_sbessel_ode_solver_create_operator (solver, 2.5, 18.3, 3, 12);
  _test_assert_accessors (op, 2.5, 18.3, 3, 12, 1.0e-9);
  ncm_sbessel_ode_operator_solve (op, rhs, &sol, &len);
  _test_assert_accessors (op, 2.5, 18.3, 3, 12, 1.0e-9);
  ncm_sbessel_ode_operator_unref (op);

  g_array_unref (rhs);
  g_array_unref (sol);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * reconfigure_operator takes the solver's current tolerance: after it, the reused
 * operator reports and uses the same tolerance as a new one (same column count).
 */
void
test_ncm_sbessel_ode_solver_reconfigure_tolerance (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op, *fresh;
  GArray *rhs = _test_rhs_const (64, 1.0);
  GArray *sol = NULL;
  gsize len;

  ncm_sbessel_ode_solver_set_tolerance (solver, 1.0e-6);
  op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 5);
  _test_assert_close (ncm_sbessel_ode_operator_get_tolerance (op), 1.0e-6, 1.0e-15, 0.0);

  ncm_sbessel_ode_solver_set_tolerance (solver, 1.0e-12);
  ncm_sbessel_ode_solver_reconfigure_operator (solver, op, TEST_A, TEST_B, 5, 5);
  _test_assert_close (ncm_sbessel_ode_operator_get_tolerance (op), 1.0e-12, 1.0e-15, 0.0);

  fresh = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 5);
  _test_assert_close (ncm_sbessel_ode_operator_get_tolerance (op), ncm_sbessel_ode_operator_get_tolerance (fresh), 1.0e-15, 0.0);

  ncm_sbessel_ode_operator_solve (op, rhs, &sol, &len);
  ncm_sbessel_ode_operator_solve (fresh, rhs, &sol, &len);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, ncm_sbessel_ode_operator_get_n_cols (fresh));

  g_array_unref (rhs);
  g_array_unref (sol);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_operator_unref (fresh);
  ncm_sbessel_ode_solver_free (solver);
}

/* At least one entry differs beyond numpy allclose (rtol 1e-10, atol 1e-8); equal lengths required */
static gboolean
_test_differ (GArray *x, GArray *y)
{
  gboolean differ = FALSE;
  guint j;

  g_assert_cmpuint (x->len, ==, y->len);

  for (j = 0; j < x->len; j++)
    differ = differ || (fabs (g_array_index (x, gdouble, j) - g_array_index (y, gdouble, j)) > 1.0e-8 + 1.0e-10 * fabs (g_array_index (y, gdouble, j)));

  return differ;
}

/*
 * The stored factorization: none on a new operator; a second right-hand side of the
 * same length keeps the column count and gives a different solution; a shorter one
 * keeps or lowers it; a longer one keeps or raises it; reconfigure_operator clears it.
 * For ell 5 and for the batch 2..8.
 */
void
test_ncm_sbessel_ode_solver_reuse (void)
{
  const gint ranges[][2] = {
    {
      5, 5
    }, {
      2, 8
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint r;

  for (r = 0; r < G_N_ELEMENTS (ranges); r++)
  {
    NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ranges[r][0], ranges[r][1]);
    GArray *rhs1              = _test_rhs_const (64, 1.0);
    GArray *rhs2              = _test_rhs_const (64, 2.0);
    GArray *rhs3              = _test_rhs_const (32, 3.0);
    GArray *rhs4              = _test_rhs_const (128, 1.0);
    GArray *sol1              = NULL;
    GArray *sol2              = NULL;
    GArray *sol               = NULL;
    gsize len;
    glong n_cols;

    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);

    ncm_sbessel_ode_operator_solve (op, rhs1, &sol1, &len);
    n_cols = ncm_sbessel_ode_operator_get_n_cols (op);
    g_assert_cmpint (n_cols, >, 0);

    ncm_sbessel_ode_operator_solve (op, rhs2, &sol2, &len);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols);
    g_assert_true (_test_differ (sol1, sol2));

    ncm_sbessel_ode_operator_solve (op, rhs3, &sol, &len);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), <=, n_cols);

    ncm_sbessel_ode_operator_solve (op, rhs4, &sol, &len);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), >=, n_cols);

    ncm_sbessel_ode_solver_reconfigure_operator (solver, op, TEST_A, TEST_B, ranges[r][0], ranges[r][1]);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);

    g_array_unref (rhs1);
    g_array_unref (rhs2);
    g_array_unref (rhs3);
    g_array_unref (rhs4);
    g_array_unref (sol1);
    g_array_unref (sol2);
    g_array_unref (sol);
    ncm_sbessel_ode_operator_unref (op);
  }

  ncm_sbessel_ode_solver_free (solver);
}

/* Two new operators and a reused one give the same solution, ell 7 */
void
test_ncm_sbessel_ode_solver_reuse_correctness (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op1  = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 7, 7);
  NcmSBesselOdeOperator *op2  = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 7, 7);
  GArray *rhs                 = _test_rhs_const (64, 1.0);
  GArray *sol1                = NULL;
  GArray *sol2                = NULL;
  GArray *sol3                = NULL;
  gsize len1, len2, len3;
  guint j;

  ncm_sbessel_ode_operator_solve (op1, rhs, &sol1, &len1);
  ncm_sbessel_ode_operator_solve (op2, rhs, &sol2, &len2);
  ncm_sbessel_ode_operator_solve (op1, rhs, &sol3, &len3);

  g_assert_cmpuint (len1, ==, len2);
  g_assert_cmpuint (len1, ==, len3);

  for (j = 0; j < len1; j++)
  {
    _test_assert_close (g_array_index (sol1, gdouble, j), g_array_index (sol2, gdouble, j), 1.0e-14, 1.0e-14);
    _test_assert_close (g_array_index (sol3, gdouble, j), g_array_index (sol1, gdouble, j), 1.0e-14, 1.0e-14);
  }

  g_array_unref (rhs);
  g_array_unref (sol1);
  g_array_unref (sol2);
  g_array_unref (sol3);
  ncm_sbessel_ode_operator_unref (op1);
  ncm_sbessel_ode_operator_unref (op2);
  ncm_sbessel_ode_solver_free (solver);
}

/* solve_endpoints stores the factorization and a second call reuses it with the same result, ell 0, 5, 10, 15 */
void
test_ncm_sbessel_ode_solver_reuse_endpoints (void)
{
  const guint ell_a[]         = {0, 5, 10, 15};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs                 = _test_rhs_const (64, 1.0);
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
  {
    NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell_a[i], ell_a[i]);
    GArray *e1                = NULL;
    GArray *e2                = NULL;
    glong n_cols;

    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);

    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &e1);
    n_cols = ncm_sbessel_ode_operator_get_n_cols (op);
    g_assert_cmpint (n_cols, >, 0);

    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &e2);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols);

    for (j = 0; j < 3; j++)
      _test_assert_close (g_array_index (e1, gdouble, j), g_array_index (e2, gdouble, j), 1.0e-14, 1.0e-14);

    g_array_unref (e1);
    g_array_unref (e2);
    ncm_sbessel_ode_operator_unref (op);
  }

  g_array_unref (rhs);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * The specialised batch sizes 2, 4, 8, 16, 32 and 64 (n = 64) against single-multipole
 * solves, for constant, linear (m, h / 4), quadratic (m^2 + h^2 / 2, m h / 2, h^2 / 6),
 * exp (x) and sin (x) forcings, through solve and solve_endpoints, and their reuse of
 * the factorization.
 */
void
test_ncm_sbessel_ode_solver_optimized_batched (void)
{
  const guint n_ell_a[]       = {2, 4, 8, 16, 32, 64};
  const gdouble lin[]         = {TEST_M, 0.25 * TEST_H};
  const gdouble quad[]        = {TEST_M * TEST_M + TEST_H * TEST_H / 2.0, 0.25 * 2.0 * TEST_M * TEST_H, TEST_H * TEST_H / 6.0};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  GArray *rhs_const           = _test_rhs_const (64, 1.0);
  GArray *rhs_const2          = _test_rhs_const (64, 2.0);
  GArray *rhs_lin             = _test_rhs_new (64, lin, G_N_ELEMENTS (lin));
  GArray *rhs_quad            = _test_rhs_new (64, quad, G_N_ELEMENTS (quad));
  GArray *rhs_exp             = _test_rhs_func (64, &_test_exp);
  GArray *rhs_sin             = _test_rhs_func (64, &_test_sin);
  guint c;

  for (c = 0; c < G_N_ELEMENTS (n_ell_a); c++)
  {
    const gint n_ell = n_ell_a[c];
    NcmSBesselOdeOperator *op;
    GArray *sol = NULL;
    gsize len;
    glong n_cols;

    _test_compare_batched_solve (solver, 5, 5 + n_ell - 1, rhs_const, 1.0e-12, 1.0e-12);
    _test_compare_batched_solve (solver, 3, 3 + n_ell - 1, rhs_lin, 1.0e-12, 1.0e-12);
    _test_compare_batched_solve (solver, 2, 2 + n_ell - 1, rhs_quad, 1.0e-12, 1.0e-12);
    _test_compare_batched_solve (solver, 1, 1 + n_ell - 1, rhs_exp, 1.0e-10, 1.0e-14);
    _test_compare_batched_solve (solver, 2, 2 + n_ell - 1, rhs_sin, 1.0e-10, 1.0e-13);

    _test_compare_batched_endpoints (solver, 3, 3 + n_ell - 1, rhs_const, 1.0e-11, 1.0e-11, TRUE);
    _test_compare_batched_endpoints (solver, 1, 1 + n_ell - 1, rhs_exp, 1.0e-11, 1.0e-11, FALSE);

    op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 2, 2 + n_ell - 1);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);
    ncm_sbessel_ode_operator_solve (op, rhs_const, &sol, &len);
    n_cols = ncm_sbessel_ode_operator_get_n_cols (op);
    g_assert_cmpint (n_cols, >, 0);
    ncm_sbessel_ode_operator_solve (op, rhs_const2, &sol, &len);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols);
    ncm_sbessel_ode_solver_reconfigure_operator (solver, op, TEST_A, TEST_B, 2, 2 + n_ell - 1);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);

    g_array_unref (sol);
    ncm_sbessel_ode_operator_unref (op);
  }

  g_array_unref (rhs_const);
  g_array_unref (rhs_const2);
  g_array_unref (rhs_lin);
  g_array_unref (rhs_quad);
  g_array_unref (rhs_exp);
  g_array_unref (rhs_sin);
  ncm_sbessel_ode_solver_free (solver);
}


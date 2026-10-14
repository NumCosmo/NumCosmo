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

#include "test_ncm_sbessel_arb_cases.h"

#define TEST_A 1.0

#define TEST_B 20.0

#define TEST_H (0.5 * (TEST_B - TEST_A))

#define TEST_M (0.5 * (TEST_A + TEST_B))

typedef gdouble (*TestFx) (gdouble x);

typedef enum _TestGreenForcing
{
  TEST_GREEN_ONE = 0,
  TEST_GREEN_X,
} TestGreenForcing;

/*
 * Stored rotations: a solve on an operator that already holds a factorization reuses
 * its rotations and extends them when the new right-hand side needs more columns. The
 * result must equal a solve on a new operator: same solution length and coefficients to
 * |reused - new| <= rtol |new| + rtol max |new| per multipole.
 */
typedef struct _TestCoeff
{
  guint k;
  gdouble value;
} TestCoeff;

/*
 * The tau constraint against Dirichlet data, ell 20, solver tolerance 1e-12. The forcing
 * is x F (x) with F a Gaussian centred on the panel [a, b] of width 0.35 h (h the
 * half-width), expanded to Chebyshev order 96 and converted to C^(2). The boundary
 * functional W(b) - W(a), W = x j_l u' - (x j_l)' u, is invariant under adding
 * homogeneous solutions, so both constraints give the same panel integral int F j_l dx.
 */
#define TEST_TAU_ELL 20

#define TEST_TAU_A 1.0e4

#define TEST_TAU_B 31622.776601683792 /* 10^4.5 */

typedef struct _TestBump
{
  gdouble a;
  gdouble b;
  gdouble center;
  gdouble width;
} TestBump;

/*
 * The tau working order and the roundoff bound it implies, on [1e4, 10^4.5] at ell 20:
 * the forcing is x F (x), F a Gaussian of width 0.25 h, expanded adaptively in x
 * (level 3, tolerance 1e-13) and converted to C^(2).
 */
typedef struct _TestXBump
{
  gdouble center;
  gdouble width;
} TestXBump;

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
void test_ncm_sbessel_ode_solver_homogeneous (void);
void test_ncm_sbessel_ode_solver_green_lu (void);
void test_ncm_sbessel_ode_solver_green_dense (void);
void test_ncm_sbessel_ode_solver_green_qr (void);
void test_ncm_sbessel_ode_solver_green_endpoints (void);
void test_ncm_sbessel_ode_solver_green_endpoints_batched (void);
void test_ncm_sbessel_ode_solver_forcing (void);
void test_ncm_sbessel_ode_solver_forcing_batched (void);
void test_ncm_sbessel_ode_solver_memory_accessors (void);
void test_ncm_sbessel_ode_solver_memory_grow (void);
void test_ncm_sbessel_ode_solver_memory_shrink (void);
void test_ncm_sbessel_ode_solver_memory_identical (void);
void test_ncm_sbessel_ode_solver_memory_same_length (void);
void test_ncm_sbessel_ode_solver_memory_increasing (void);
void test_ncm_sbessel_ode_solver_memory_varying (void);
void test_ncm_sbessel_ode_solver_memory_batched (void);
void test_ncm_sbessel_ode_solver_memory_endpoints (void);
void test_ncm_sbessel_ode_solver_memory_reset (void);
void test_ncm_sbessel_ode_solver_memory_stress (void);
void test_ncm_sbessel_ode_solver_rotations_identical (void);
void test_ncm_sbessel_ode_solver_rotations_extend (void);
void test_ncm_sbessel_ode_solver_rotations_sequence (void);
void test_ncm_sbessel_ode_solver_rotations_batched_extend (void);
void test_ncm_sbessel_ode_solver_rotations_leading_zeros (void);
void test_ncm_sbessel_ode_solver_rotations_convergence (void);
void test_ncm_sbessel_ode_solver_rotations_batched_sequence (void);
void test_ncm_sbessel_ode_solver_rotations_high_ell (void);
void test_ncm_sbessel_ode_solver_tau_deep_panel (void);
void test_ncm_sbessel_ode_solver_tau_narrow_panel (void);
void test_ncm_sbessel_ode_solver_tau_batched (void);
void test_ncm_sbessel_ode_solver_tau_shallow_panel (void);
void test_ncm_sbessel_ode_solver_tau_last_max_coeff (void);
void test_ncm_sbessel_ode_solver_tau_last_max_coeff_reconfigure (void);
void test_ncm_sbessel_ode_solver_tau_pinned (void);
void test_ncm_sbessel_ode_solver_tau_propagation (void);
void test_ncm_sbessel_ode_solver_floor_default (void);
void test_ncm_sbessel_ode_solver_floor_working_order (void);
void test_ncm_sbessel_ode_solver_floor_deriv_error (void);
void test_ncm_sbessel_ode_solver_floor_deriv_error_per_ell (void);
void test_ncm_sbessel_ode_solver_floor_matrix_constraint (void);
void test_ncm_sbessel_ode_solver_floor_matrix_band (void);
void test_ncm_sbessel_ode_solver_deriv_bound_never_below (void);
void test_ncm_sbessel_ode_solver_deriv_bound_separates (void);
void test_ncm_sbessel_ode_solver_conjugate_dirichlet_fails (void);
void test_ncm_sbessel_ode_solver_conjugate_pinned (void);
void test_ncm_sbessel_ode_solver_conjugate_away (void);

static const guint test_zero_ells[] = {0, 1, 5, 10, 20};

static const guint test_zero_ns[] = {32, 64, 128};

/*
 * Green's identity: with zero Dirichlet data and forcing g, the solution u satisfies
 * b j_l(b) u'(b) - a j_l(a) u'(a) = int_a^b g(x) j_l(x) / x dx. The forcing 1 (rhs[2] = 1)
 * gives int_1^20 j_l / x dx and the forcing x (rhs[2] = m, rhs[3] = h / 4) gives
 * int_1^20 j_l dx, certified by make_sbessel_arb_cases.py ode.
 */
static const TestArbCell test_arb_cells[] = {
  /* BEGIN make_sbessel_arb_cases.py ode */
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 0, 1.0, 0.5028396202159003},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 1, 1.0, 0.4607245247899321},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 2, 1.0, 0.3020747659379493},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 3, 1.0, 0.19574779143238527},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 4, 1.0, 0.1330169948034595},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 5, 1.0, 0.0956047082381122},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 6, 1.0, 0.07484428607483702},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 7, 1.0, 0.06312019765201392},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 8, 1.0, 0.05352439806530906},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 9, 1.0, 0.04344865430572876},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 10, 1.0, 0.03453229048880144},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 11, 1.0, 0.029113904722431434},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 12, 1.0, 0.027378313816413244},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 13, 1.0, 0.02736466382665843},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 14, 1.0, 0.026730174217619293},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 15, 1.0, 0.02423940889318163},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 16, 1.0, 0.02004830157463489},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 17, 1.0, 0.015109047820065719},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 18, 1.0, 0.01042813195282212},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 19, 1.0, 0.006635382719972067},
  {"green_inv_x", TEST_ARB_POWER, 1.0, 20.0, -1.0, 0.0, 0, FALSE, 20, 1.0, 0.003917880565087803},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 0, 1.0, 0.6021586306762569},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 1, 1.0, 0.7958237222715151},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 2, 1.0, 0.7800149436935394},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 3, 1.0, 0.7145501074182314},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 4, 1.0, 0.5902195963331576},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 5, 1.0, 0.48260284581290397},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 6, 1.0, 0.4614321942860765},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 7, 1.0, 0.49037287315997735},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 8, 1.0, 0.48537077049413235},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 9, 1.0, 0.41954189395027663},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 10, 1.0, 0.34015366131471403},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 11, 1.0, 0.3056362063145535},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 12, 1.0, 0.3294661473012089},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 13, 1.0, 0.3788216390957776},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 14, 1.0, 0.40937977601856873},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 15, 1.0, 0.39635341321518186},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 16, 1.0, 0.34204189967006177},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 17, 1.0, 0.2652405387477694},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 18, 1.0, 0.1867747740322384},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 19, 1.0, 0.12060034350664904},
  {"green_one", TEST_ARB_TOPHAT, 1.0, 20.0, 0.0, 0.0, 0, FALSE, 20, 1.0, 0.07200515204667218},
  {"rhs_quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, TRUE, 0, 1.0, 0.13222024405474772},
  {"rhs_quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, TRUE, 5, 1.0, 1.9425860132872468},
  {"rhs_quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, TRUE, 10, 1.0, 3.1777343624297174},
  {"rhs_quadratic", TEST_ARB_POWER, 1.0, 20.0, 2.0, 0.0, 0, TRUE, 15, 1.0, 6.583316123648764},
  {"rhs_cubic", TEST_ARB_POWER, 1.0, 20.0, 3.0, 0.0, 0, TRUE, 0, 1.0, -7.549864664479969},
  {"rhs_cubic", TEST_ARB_POWER, 1.0, 20.0, 3.0, 0.0, 0, TRUE, 5, 1.0, -3.7098496929237466},
  {"rhs_cubic", TEST_ARB_POWER, 1.0, 20.0, 3.0, 0.0, 0, TRUE, 10, 1.0, 25.794284875856388},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 0, 1.0, 274557.9565997802},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 2, 1.0, -405905.43815627106},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 3, 1.0, 747820.6721892703},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 4, 1.0, 679786.3761393836},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 5, 1.0, -443105.4051754701},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 6, 1.0, -943851.7547431279},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 7, 1.0, -182667.5296970764},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 8, 1.0, 826212.2069565651},
  {"rhs_exp", TEST_ARB_EXP, 1.0, 20.0, -1.0, 0.0, 0, TRUE, 10, 1.0, 57612.98468184247},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 0, 1.0, 0.6479721092843543},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 2, 1.0, 0.006110573663597096},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 3, 1.0, -0.08111944835354816},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 4, 1.0, -0.02449089352950558},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 5, 1.0, 0.02363007923070094},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 6, 1.0, 0.02066418140989152},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 7, 1.0, -0.0021544967305790973},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 8, 1.0, -0.012897440654211665},
  {"rhs_sin", TEST_ARB_SIN, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 10, 1.0, 0.0022078534135041633},
  {"rhs_cos", TEST_ARB_COS, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 0, 1.0, 0.0320600635222009},
  {"rhs_cos", TEST_ARB_COS, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 5, 1.0, 0.022739400291948583},
  {"rhs_cos", TEST_ARB_COS, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 10, 1.0, 0.008890728940125823},
  {"rhs_gauss", TEST_ARB_EXP2, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 0, 1.0, 0.08192634928161613},
  {"rhs_gauss", TEST_ARB_EXP2, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 5, 1.0, 4.873451437933564e-05},
  {"rhs_gauss", TEST_ARB_EXP2, 1.0, 20.0, 1.0, 0.0, 0, TRUE, 10, 1.0, 7.802205812987952e-10},
  {"rhs_lorentz", TEST_ARB_LORENTZ, 1.0, 20.0, 0.0, 1.0, 0, TRUE, 0, 1.0, 0.17909599346737495},
  {"rhs_lorentz", TEST_ARB_LORENTZ, 1.0, 20.0, 0.0, 1.0, 0, TRUE, 5, 1.0, 0.0037861935193706244},
  {"rhs_lorentz", TEST_ARB_LORENTZ, 1.0, 20.0, 0.0, 1.0, 0, TRUE, 10, 1.0, 0.0003449988561628076},
  {"tau_narrow", TEST_ARB_GAUSS, 100.0, 316.22776601683796, 208.11388300841898, 37.83985905294664, 0, FALSE, 20, 1.0, 3.470494682971913e-05},
  {"deriv_bound_30", TEST_ARB_GAUSS, 1000.0, 1030.0, 1015.0, 3.75, 0, FALSE, 20, 1.0, -3.8384799978730775e-06},
  {"deriv_bound_60", TEST_ARB_GAUSS, 1000.0, 1060.0, 1030.0, 7.5, 0, FALSE, 20, 1.0, 1.2305957763076221e-07},
  {"deriv_bound_100", TEST_ARB_GAUSS, 1000.0, 1100.0, 1050.0, 12.5, 0, FALSE, 20, 1.0, -2.530161116525624e-07},
  {"deriv_bound_150", TEST_ARB_GAUSS, 1000.0, 1150.0, 1075.0, 18.75, 0, FALSE, 20, 1.0, -2.2501000047625404e-07},
  {"deriv_bound_250", TEST_ARB_GAUSS, 1000.0, 1250.0, 1125.0, 31.25, 0, FALSE, 20, 1.0, -1.675615103467387e-07},
  {"deriv_bound_3000", TEST_ARB_GAUSS, 1000.0, 4000.0, 2500.0, 375.0, 0, FALSE, 20, 1.0, 1.824908439349761e-07},
  {"conjugate", TEST_ARB_GAUSS, 200.0, 225.13441585600583, 212.56720792800292, 4.39852277480102, 0, FALSE, 2, 1.0, -6.780495041869413e-05},
  {"conjugate_away", TEST_ARB_GAUSS, 200.0, 225.14441585600582, 212.5722079280029, 4.400272774801018, 0, FALSE, 2, 1.0, -6.719000261353115e-05},
  /* END make_sbessel_arb_cases.py ode */
};

static const guint test_green_ells[] = {0, 1, 5, 10, 20};

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
  g_test_add_func ("/ncm/sbessel_ode_solver/homogeneous", test_ncm_sbessel_ode_solver_homogeneous);
  g_test_add_func ("/ncm/sbessel_ode_solver/green/lu", test_ncm_sbessel_ode_solver_green_lu);
  g_test_add_func ("/ncm/sbessel_ode_solver/green/dense", test_ncm_sbessel_ode_solver_green_dense);
  g_test_add_func ("/ncm/sbessel_ode_solver/green/qr", test_ncm_sbessel_ode_solver_green_qr);
  g_test_add_func ("/ncm/sbessel_ode_solver/green/endpoints", test_ncm_sbessel_ode_solver_green_endpoints);
  g_test_add_func ("/ncm/sbessel_ode_solver/green/endpoints_batched", test_ncm_sbessel_ode_solver_green_endpoints_batched);
  g_test_add_func ("/ncm/sbessel_ode_solver/forcing", test_ncm_sbessel_ode_solver_forcing);
  g_test_add_func ("/ncm/sbessel_ode_solver/forcing/batched", test_ncm_sbessel_ode_solver_forcing_batched);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/accessors", test_ncm_sbessel_ode_solver_memory_accessors);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/grow", test_ncm_sbessel_ode_solver_memory_grow);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/shrink", test_ncm_sbessel_ode_solver_memory_shrink);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/identical", test_ncm_sbessel_ode_solver_memory_identical);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/same_length", test_ncm_sbessel_ode_solver_memory_same_length);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/increasing", test_ncm_sbessel_ode_solver_memory_increasing);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/varying", test_ncm_sbessel_ode_solver_memory_varying);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/batched", test_ncm_sbessel_ode_solver_memory_batched);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/endpoints", test_ncm_sbessel_ode_solver_memory_endpoints);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/reset", test_ncm_sbessel_ode_solver_memory_reset);
  g_test_add_func ("/ncm/sbessel_ode_solver/memory/stress", test_ncm_sbessel_ode_solver_memory_stress);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/identical", test_ncm_sbessel_ode_solver_rotations_identical);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/extend", test_ncm_sbessel_ode_solver_rotations_extend);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/sequence", test_ncm_sbessel_ode_solver_rotations_sequence);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/batched_extend", test_ncm_sbessel_ode_solver_rotations_batched_extend);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/leading_zeros", test_ncm_sbessel_ode_solver_rotations_leading_zeros);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/convergence", test_ncm_sbessel_ode_solver_rotations_convergence);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/batched_sequence", test_ncm_sbessel_ode_solver_rotations_batched_sequence);
  g_test_add_func ("/ncm/sbessel_ode_solver/rotations/high_ell", test_ncm_sbessel_ode_solver_rotations_high_ell);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/deep_panel", test_ncm_sbessel_ode_solver_tau_deep_panel);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/narrow_panel", test_ncm_sbessel_ode_solver_tau_narrow_panel);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/batched", test_ncm_sbessel_ode_solver_tau_batched);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/shallow_panel", test_ncm_sbessel_ode_solver_tau_shallow_panel);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/last_max_coeff", test_ncm_sbessel_ode_solver_tau_last_max_coeff);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/last_max_coeff_reconfigure", test_ncm_sbessel_ode_solver_tau_last_max_coeff_reconfigure);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/pinned", test_ncm_sbessel_ode_solver_tau_pinned);
  g_test_add_func ("/ncm/sbessel_ode_solver/tau/propagation", test_ncm_sbessel_ode_solver_tau_propagation);
  g_test_add_func ("/ncm/sbessel_ode_solver/floor/default", test_ncm_sbessel_ode_solver_floor_default);
  g_test_add_func ("/ncm/sbessel_ode_solver/floor/working_order", test_ncm_sbessel_ode_solver_floor_working_order);
  g_test_add_func ("/ncm/sbessel_ode_solver/floor/deriv_error", test_ncm_sbessel_ode_solver_floor_deriv_error);
  g_test_add_func ("/ncm/sbessel_ode_solver/floor/deriv_error_per_ell", test_ncm_sbessel_ode_solver_floor_deriv_error_per_ell);
  g_test_add_func ("/ncm/sbessel_ode_solver/floor/matrix_constraint", test_ncm_sbessel_ode_solver_floor_matrix_constraint);
  g_test_add_func ("/ncm/sbessel_ode_solver/floor/matrix_band", test_ncm_sbessel_ode_solver_floor_matrix_band);
  g_test_add_func ("/ncm/sbessel_ode_solver/deriv_bound/never_below", test_ncm_sbessel_ode_solver_deriv_bound_never_below);
  g_test_add_func ("/ncm/sbessel_ode_solver/deriv_bound/separates", test_ncm_sbessel_ode_solver_deriv_bound_separates);
  g_test_add_func ("/ncm/sbessel_ode_solver/conjugate/dirichlet_fails", test_ncm_sbessel_ode_solver_conjugate_dirichlet_fails);
  g_test_add_func ("/ncm/sbessel_ode_solver/conjugate/pinned", test_ncm_sbessel_ode_solver_conjugate_pinned);
  g_test_add_func ("/ncm/sbessel_ode_solver/conjugate/away", test_ncm_sbessel_ode_solver_conjugate_away);

  g_test_run ();
}

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

static gdouble
_test_cheb_func (gpointer user_data, gdouble t)
{
  TestFx f = (TestFx) user_data;

  return f (TEST_M + TEST_H * t);
}

/* F (t) on [-1, 1] as Chebyshev coefficients of order n - 2 converted to C^(2), after zero endpoint data */
static GArray *
_test_rhs_spectral (guint n, NcmSpectralF F, gpointer user_data)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *cheb          = NULL;
  GArray *gegen         = NULL;
  GArray *rhs;

  ncm_spectral_compute_chebyshev_coeffs (spectral, F, -1.0, 1.0, n - 2, &cheb, user_data);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (cheb, &gegen);
  rhs = _test_rhs_new (n, (const gdouble *) gegen->data, GSL_MIN (gegen->len, n - 2));

  g_array_unref (cheb);
  g_array_unref (gegen);
  ncm_spectral_free (spectral);

  return rhs;
}

/* f (x) on [1, 20] as the forcing */
static GArray *
_test_rhs_func (guint n, TestFx f)
{
  return _test_rhs_spectral (n, &_test_cheb_func, (gpointer) f);
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

static GArray *
_test_green_rhs (TestGreenForcing forcing)
{
  const gdouble f_x[] = {TEST_M, 0.25 * TEST_H};

  return (forcing == TEST_GREEN_ONE) ? _test_rhs_const (128, 1.0) : _test_rhs_new (128, f_x, G_N_ELEMENTS (f_x));
}

/* The certified int g j_l / x dx of the forcing at ell */
static gdouble
_test_green_value (TestGreenForcing forcing, guint ell)
{
  guint n;
  const TestArbCell *c = test_arb_cells_find (test_arb_cells, G_N_ELEMENTS (test_arb_cells),
                                              (forcing == TEST_GREEN_ONE) ? "green_inv_x" : "green_one", &n);

  g_assert_cmpuint (ell, <, n);
  g_assert_cmpuint (c[ell].ell, ==, ell);

  return c[ell].value;
}

/* b j_l(b) u'(b) - a j_l(a) u'(a) held to a certified integral, rtol 1e-8, atol 1e-20 */
static void
_test_green_assert (guint ell, gdouble du_a, gdouble du_b, gdouble value)
{
  const gdouble lhs = TEST_B * ncm_sf_sbessel (ell, TEST_B) * du_b - TEST_A * ncm_sf_sbessel (ell, TEST_A) * du_a;

  _test_assert_close (lhs, value, 1.0e-8, 1.0e-20);
}

static void
_test_green_check (TestGreenForcing forcing, guint ell, gdouble du_a, gdouble du_b)
{
  _test_green_assert (ell, du_a, du_b, _test_green_value (forcing, ell));
}

/* u'(a), u'(b) of Chebyshev coefficients on [1, 20] */
static void
_test_green_check_coeffs (TestGreenForcing forcing, guint ell, GArray *u)
{
  _test_green_check (forcing, ell, ncm_spectral_chebyshev_deriv (u, -1.0) / TEST_H, ncm_spectral_chebyshev_deriv (u, +1.0) / TEST_H);
}

/* The forcing g (x) of a certified row at x = m + h t: the row's shape without the 1 / x */
static gdouble
_test_cheb_cell (gpointer user_data, gdouble t)
{
  const TestArbCell *c = (const TestArbCell *) user_data;

  return test_arb_cell_shape (c, TEST_M + TEST_H * t);
}

/*
 * Checks one solution for the forcing of row @c: finite, not identically zero, u(a) =
 * u(b) = 0 to @bc_atol, and Green's identity against the row's certified
 * int_1^20 g j_l / x dx.
 */
static void
_test_forcing_check (const TestArbCell *c, GArray *u, gdouble bc_atol)
{
  guint j;
  gboolean nonzero = FALSE;

  for (j = 0; j < u->len; j++)
  {
    g_assert_true (gsl_finite (g_array_index (u, gdouble, j)));
    nonzero = nonzero || (fabs (g_array_index (u, gdouble, j)) > 1.0e-10);
  }

  g_assert_true (nonzero);
  g_assert_cmpfloat (fabs (ncm_spectral_chebyshev_eval (u, -1.0)), <=, bc_atol);
  g_assert_cmpfloat (fabs (ncm_spectral_chebyshev_eval (u, +1.0)), <=, bc_atol);

  _test_green_assert (c->ell, ncm_spectral_chebyshev_deriv (u, -1.0) / TEST_H, ncm_spectral_chebyshev_deriv (u, +1.0) / TEST_H, c->value);
}

/*
 * Memory and reuse of one operator across solves of different right-hand sides. Every
 * solve on the reused operator is also done on a new operator and the two agree to
 * |reused - new| <= 1e-12 |new| + 1e-12 max |new| (the shorter zero-padded); the stored
 * capacity (get_operator_size) never decreases and n_cols is positive after a solve.
 */
static GArray *
_test_rhs_unit (guint order, guint min_len, gdouble value)
{
  GArray *rhs = _test_rhs_new (GSL_MAX (min_len, order + 4), NULL, 0);

  g_array_index (rhs, gdouble, order) = value;

  return rhs;
}

static void
_test_assert_same_padded (GArray *x, gsize x_len, GArray *y, gsize y_len, guint n_ell)
{
  const gsize len = GSL_MAX (x_len, y_len);
  guint l, j;

  for (l = 0; l < n_ell; l++)
  {
    gdouble scale = 0.0;

    for (j = 0; j < y_len; j++)
      scale = GSL_MAX (scale, fabs (g_array_index (y, gdouble, l * y_len + j)));

    for (j = 0; j < len; j++)
    {
      const gdouble xv = (j < x_len) ? g_array_index (x, gdouble, l * x_len + j) : 0.0;
      const gdouble yv = (j < y_len) ? g_array_index (y, gdouble, l * y_len + j) : 0.0;

      g_assert_true (gsl_finite (xv));
      _test_assert_close (xv, yv, 1.0e-12, 1.0e-12 * scale);
    }
  }
}

/*
 * Solves @rhs on @op, checks the invariants (capacity >= *capacity, n_cols > 0, a
 * non-empty finite solution) and the agreement with a new operator; updates *capacity
 * and returns the solution (*len per multipole).
 */
static GArray *
_test_memory_solve (NcmSBesselOdeSolver *solver, NcmSBesselOdeOperator *op, GArray *rhs, gsize *capacity, gsize *len)
{
  GArray *sol   = NULL;
  GArray *fresh = NULL;
  gsize fresh_len;
  gint ell_min, ell_max;
  NcmSBesselOdeOperator *fop;

  ncm_sbessel_ode_operator_get_ell_range (op, &ell_min, &ell_max);
  ncm_sbessel_ode_operator_solve (op, rhs, &sol, len);

  g_assert_cmpuint (ncm_sbessel_ode_operator_get_operator_size (op), >=, *capacity);
  *capacity = ncm_sbessel_ode_operator_get_operator_size (op);
  g_assert_cmpuint (*capacity, >, 0);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), >, 0);
  g_assert_cmpuint (*len, >, 0);

  fop = _test_solve (solver, TEST_A, TEST_B, ell_min, ell_max, rhs, &fresh, &fresh_len);
  _test_assert_same_padded (sol, *len, fresh, fresh_len, ell_max - ell_min + 1);

  g_array_unref (fresh);
  ncm_sbessel_ode_operator_unref (fop);

  return sol;
}

/* Runs @n right-hand sides rhs[orders[i]] = values[i] (length at least @min_len) through one operator */
static void
_test_memory_sequence (gint ell_min, gint ell_max, const guint *orders, const gdouble *values, guint n, guint min_len)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, ell_min, ell_max);
  gsize capacity              = 0;
  guint i;

  for (i = 0; i < n; i++)
  {
    GArray *rhs = _test_rhs_unit (orders[i], min_len, (values != NULL) ? values[i] : 1.0);
    gsize len;
    GArray *sol = _test_memory_solve (solver, op, rhs, &capacity, &len);

    g_array_unref (sol);
    g_array_unref (rhs);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* A right-hand side of length n with zero endpoint data and the listed coefficients */
static GArray *
_test_rhs_sparse (guint n, const TestCoeff *c, guint n_c)
{
  GArray *rhs = _test_rhs_new (n, NULL, 0);
  guint i;

  for (i = 0; i < n_c; i++)
    g_array_index (rhs, gdouble, c[i].k) = c[i].value;

  return rhs;
}

static void
_test_assert_same (GArray *x, gsize x_len, GArray *y, gsize y_len, guint n_ell, gdouble rtol)
{
  guint l, j;

  g_assert_cmpuint (x_len, ==, y_len);

  for (l = 0; l < n_ell; l++)
  {
    gdouble scale = 0.0;

    for (j = 0; j < y_len; j++)
      scale = GSL_MAX (scale, fabs (g_array_index (y, gdouble, l * y_len + j)));

    for (j = 0; j < y_len; j++)
      _test_assert_close (g_array_index (x, gdouble, l * x_len + j), g_array_index (y, gdouble, l * y_len + j), rtol, rtol * scale);
  }
}

/* The largest |coefficient| of row l */
static gdouble
_test_row_absmax (GArray *sol, gsize len, guint l)
{
  GArray *row         = _test_row (sol, len, l);
  const gdouble scale = _test_absmax (row);

  g_array_unref (row);

  return scale;
}

/*
 * Solves @rhs on @op and on a new operator over the same interval and multipoles,
 * checks the agreement to @rtol and returns the reused solution; *n_cols_new receives
 * the new operator's n_cols.
 */
static GArray *
_test_reused_vs_new (NcmSBesselOdeSolver *solver, NcmSBesselOdeOperator *op, GArray *rhs, gdouble rtol, gsize *len, glong *n_cols_new)
{
  GArray *sol   = NULL;
  GArray *fresh = NULL;
  gsize fresh_len;
  gdouble a, b;
  gint ell_min, ell_max;
  NcmSBesselOdeOperator *fop;

  ncm_sbessel_ode_operator_get_interval (op, &a, &b);
  ncm_sbessel_ode_operator_get_ell_range (op, &ell_min, &ell_max);

  ncm_sbessel_ode_operator_solve (op, rhs, &sol, len);
  fop = _test_solve (solver, a, b, ell_min, ell_max, rhs, &fresh, &fresh_len);

  _test_assert_same (sol, *len, fresh, fresh_len, ell_max - ell_min + 1, rtol);

  if (n_cols_new != NULL)
    *n_cols_new = ncm_sbessel_ode_operator_get_n_cols (fop);

  g_array_unref (fresh);
  ncm_sbessel_ode_operator_unref (fop);

  return sol;
}

static gdouble
_test_bump_x_cheb (gpointer user_data, gdouble t)
{
  const TestBump *p = (const TestBump *) user_data;
  const gdouble x   = 0.5 * (p->a + p->b) + 0.5 * (p->b - p->a) * t;

  return x * exp (-0.5 * gsl_pow_2 ((x - p->center) / p->width));
}

/* x F (x) on [a, b], F the Gaussian (center, width), at Chebyshev order 96 */
static GArray *
_test_tau_rhs (const TestBump *p)
{
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *cheb          = NULL;
  GArray *gegen         = NULL;
  GArray *rhs;

  ncm_spectral_compute_chebyshev_coeffs (spectral, &_test_bump_x_cheb, -1.0, 1.0, 96, &cheb, (gpointer) p);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (cheb, &gegen);
  rhs = _test_rhs_new (gegen->len + 2, (const gdouble *) gegen->data, gegen->len);

  g_array_unref (cheb);
  g_array_unref (gegen);
  ncm_spectral_free (spectral);

  return rhs;
}

static TestBump
_test_tau_bump (gdouble a, gdouble b)
{
  TestBump p = {a, b, 0.5 * (a + b), 0.35 * 0.5 * (b - a)};

  return p;
}

/* j_l'(x), as scipy's spherical_jn(derivative=True) computes it */
static gdouble
_test_jl_deriv (guint ell, gdouble x)
{
  if (ell == 0)
    return -ncm_sf_sbessel (1, x);

  return ncm_sf_sbessel (ell - 1, x) - (ell + 1.0) / x * ncm_sf_sbessel (ell, x);
}

/* W(b) - W(a), W = x j_l u' - (j_l + x j_l') u, for the Chebyshev series @u on [a, b] */
static gdouble
_test_boundary_functional (GArray *u, gdouble a, gdouble b, guint ell)
{
  const gdouble h = 0.5 * (b - a);
  gdouble w[2];
  guint i;

  for (i = 0; i < 2; i++)
  {
    const gdouble t  = (i == 0) ? -1.0 : 1.0;
    const gdouble x  = (i == 0) ? a : b;
    const gdouble jl = ncm_sf_sbessel (ell, x);

    w[i] = x * jl * ncm_spectral_chebyshev_deriv (u, t) / h - (jl + x * _test_jl_deriv (ell, x)) * ncm_spectral_chebyshev_eval (u, t);
  }

  return w[1] - w[0];
}

/* A solver at tolerance 1e-12 with the given default constraint */
static NcmSBesselOdeSolver *
_test_tau_solver (NcmSBesselOdeConstraint constraint)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();

  ncm_sbessel_ode_solver_set_tolerance (solver, 1.0e-12);
  ncm_sbessel_ode_solver_set_default_constraint (solver, constraint);

  return solver;
}

/* Solves the bump on [a, b] for ell_min..ell_max; returns the operator */
static NcmSBesselOdeOperator *
_test_tau_solve (const TestBump *p, gint ell_min, gint ell_max, NcmSBesselOdeConstraint constraint, GArray **u, gsize *len)
{
  NcmSBesselOdeSolver *solver = _test_tau_solver (constraint);
  GArray *rhs                 = _test_tau_rhs (p);
  NcmSBesselOdeOperator *op   = _test_solve (solver, p->a, p->b, ell_min, ell_max, rhs, u, len);

  g_array_unref (rhs);
  ncm_sbessel_ode_solver_free (solver);

  return op;
}

/* max |a_j| of the tau solution over that of the Dirichlet one, at ell 2 on [a, b] */
static gdouble
_test_tau_contamination (gdouble a, gdouble b)
{
  const TestBump p = _test_tau_bump (a, b);
  GArray *u_d      = NULL;
  GArray *u_t      = NULL;
  gsize n;
  NcmSBesselOdeOperator *op_d = _test_tau_solve (&p, 2, 2, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET, &u_d, &n);
  NcmSBesselOdeOperator *op_t = _test_tau_solve (&p, 2, 2, NCM_SBESSEL_ODE_CONSTRAINT_TAU, &u_t, &n);
  const gdouble ratio         = _test_absmax (u_t) / _test_absmax (u_d);

  g_array_unref (u_d);
  g_array_unref (u_t);
  ncm_sbessel_ode_operator_unref (op_d);
  ncm_sbessel_ode_operator_unref (op_t);

  return ratio;
}

static gdouble
_test_xbump (gpointer user_data, gdouble x)
{
  const TestXBump *p = (const TestXBump *) user_data;

  return x * exp (-0.5 * gsl_pow_2 ((x - p->center) / p->width));
}

/* The right-hand side and, in *n_forcing, the number of Chebyshev coefficients of x F */
static GArray *
_test_floor_rhs (gdouble a, gdouble b, gdouble center, gdouble width, guint *n_forcing)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSpectral *spectral       = ncm_sbessel_ode_solver_peek_spectral (solver);
  TestXBump p                 = {center, width};
  GArray *cheb                = NULL;
  GArray *gegen               = NULL;
  GArray *rhs;

  ncm_spectral_compute_chebyshev_coeffs_adaptive (spectral, &_test_xbump, a, b, 3, 1.0e-13, &cheb, &p);
  ncm_spectral_chebT_to_gegenbauer_alpha2 (cheb, &gegen);
  g_assert_cmpuint (gegen->len, ==, cheb->len);
  rhs        = _test_rhs_new (gegen->len + 2, (const gdouble *) gegen->data, gegen->len);
  *n_forcing = cheb->len;

  g_array_unref (cheb);
  g_array_unref (gegen);
  ncm_sbessel_ode_solver_free (solver);

  return rhs;
}

static GArray *
_test_floor_rhs_default (guint *n_forcing)
{
  return _test_floor_rhs (TEST_TAU_A, TEST_TAU_B, 0.5 * (TEST_TAU_A + TEST_TAU_B), 0.25 * 0.5 * (TEST_TAU_B - TEST_TAU_A), n_forcing);
}

/* sum_j j^2 |a_j| / h of row l of a solution */
static gdouble
_test_weighted_sum (GArray *u, gsize n, guint l, gdouble h)
{
  gdouble s = 0.0;
  guint j;

  for (j = 0; j < n; j++)
    s += (gdouble) j * j * fabs (g_array_index (u, gdouble, l * n + j));

  return s / h;
}

/* The number of nonzero entries of row @i */
static guint
_test_row_nonzeros (NcmMatrix *m, guint i)
{
  guint count = 0;
  guint j;

  for (j = 0; j < ncm_matrix_ncols (m); j++)
    count += (ncm_matrix_get (m, i, j) != 0.0) ? 1 : 0;

  return count;
}

/*
 * The derivative-error bound: under the tau constraint the reported
 * get_last_deriv_error, scaled as eps deriv_error max (|x j_l|) / |I|, is never below
 * the realized relative error of the panel integral when that error exceeds 1e-9. The
 * panels [1e3, 1e3 + span] hold a Gaussian of width 0.25 h at ell 20; the integrals are
 * the certified deriv_bound_<span> cells.
 */
static void
_test_deriv_bound_solve (const TestArbCell *c, gdouble *value, gdouble *deriv_error)
{
  NcmSBesselOdeSolver *solver = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_TAU);
  guint n_forcing;
  GArray *rhs = _test_floor_rhs (c->a, c->b, c->p1, c->p2, &n_forcing);
  GArray *u   = NULL;
  gsize n;
  NcmSBesselOdeOperator *op = _test_solve (solver, c->a, c->b, TEST_TAU_ELL, TEST_TAU_ELL, rhs, &u, &n);

  *value       = _test_boundary_functional (u, c->a, c->b, TEST_TAU_ELL);
  *deriv_error = ncm_sbessel_ode_operator_get_last_deriv_error (op, 0);

  g_array_unref (u);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

static const TestArbCell *
_test_deriv_bound_cell (guint span)
{
  gchar *name          = g_strdup_printf ("deriv_bound_%u", span);
  const TestArbCell *c = test_arb_cells_find_ell (test_arb_cells, G_N_ELEMENTS (test_arb_cells), name, TEST_TAU_ELL);

  g_free (name);

  return c;
}

/*
 * Conjugate points: on [200, 200 + s], s the root of j_2(a) y_2(b) - j_2(b) y_2(a)
 * (certified row "conjugate"), the two-point Dirichlet problem at ell 2 has no unique
 * solution, and its error grows as the machine epsilon times the boundary condition
 * number. The pinned constraint imposes no endpoint condition, so it has no such
 * determinant. Solver tolerance 1e-12; the forcing is x F (x), F the row's Gaussian.
 */
static gdouble
_test_conjugate_integral (const TestArbCell *c, NcmSBesselOdeConstraint constraint, glong pin1, glong pin2)
{
  NcmSBesselOdeSolver *solver = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, c->a, c->b, c->ell, c->ell);
  guint n_forcing;
  GArray *rhs = _test_floor_rhs (c->a, c->b, c->p1, c->p2, &n_forcing);
  GArray *u   = NULL;
  gsize n;
  gdouble value;

  if (constraint == NCM_SBESSEL_ODE_CONSTRAINT_PINNED)
    ncm_sbessel_ode_operator_set_pinned_constraint (op, pin1, pin2);
  else
    ncm_sbessel_ode_operator_set_constraint (op, constraint);

  ncm_sbessel_ode_operator_solve (op, rhs, &u, &n);
  value = _test_boundary_functional (u, c->a, c->b, c->ell);

  g_array_unref (u);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);

  return value;
}

static const TestArbCell *
_test_conjugate_cell (const gchar *name)
{
  return test_arb_cells_find_ell (test_arb_cells, G_N_ELEMENTS (test_arb_cells), name, 2);
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

/*
 * With Dirichlet data a j_l(a), b j_l(b) and no forcing, the LU solution of the 128 x 128
 * operator matrix is u = x j_l(x): checked at seven interior points and at both ends,
 * ell 0..20, j_l from ncm_sf_sbessel().
 */
void
test_ncm_sbessel_ode_solver_homogeneous (void)
{
  const gdouble t_a[]         = {-0.8, -0.5, -0.2, 0.0, 0.2, 0.5, 0.8};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint ell, i;

  for (ell = 0; ell <= 20; ell++)
  {
    const gdouble u_a = TEST_A * ncm_sf_sbessel (ell, TEST_A);
    const gdouble u_b = TEST_B * ncm_sf_sbessel (ell, TEST_B);
    GArray *rhs       = _test_rhs_new (128, NULL, 0);
    NcmVector *x;
    GArray *u;

    g_array_index (rhs, gdouble, 0) = u_a;
    g_array_index (rhs, gdouble, 1) = u_b;
    x                               = _test_lu_solve (solver, ell, 128, rhs);
    u                               = ncm_vector_dup_array (x);

    for (i = 0; i < G_N_ELEMENTS (t_a); i++)
    {
      const gdouble xp = TEST_M + TEST_H * t_a[i];

      _test_assert_close (ncm_spectral_chebyshev_eval (u, t_a[i]), xp * ncm_sf_sbessel (ell, xp), 1.0e-10, 1.0e-15);
    }

    _test_assert_close (ncm_spectral_chebyshev_eval (u, -1.0), u_a, 1.0e-13, 1.0e-13);
    _test_assert_close (ncm_spectral_chebyshev_eval (u, +1.0), u_b, 1.0e-13, 1.0e-13);

    g_array_unref (rhs);
    g_array_unref (u);
    ncm_vector_free (x);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* LU solve of the 128 x 128 operator matrix, ell 0..20 */
void
test_ncm_sbessel_ode_solver_green_lu (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint f, ell;

  for (f = TEST_GREEN_ONE; f <= TEST_GREEN_X; f++)
  {
    GArray *rhs = _test_green_rhs (f);

    for (ell = 0; ell <= 20; ell++)
    {
      NcmVector *x = _test_lu_solve (solver, ell, 128, rhs);
      GArray *u    = ncm_vector_dup_array (x);

      _test_green_check_coeffs (f, ell, u);

      g_array_unref (u);
      ncm_vector_free (x);
    }

    g_array_unref (rhs);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* solve_dense, ell 0, 1, 5, 10, 20 */
void
test_ncm_sbessel_ode_solver_green_dense (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint f, i;

  for (f = TEST_GREEN_ONE; f <= TEST_GREEN_X; f++)
  {
    GArray *rhs      = _test_green_rhs (f);
    NcmVector *rhs_v = ncm_vector_new_array (rhs);

    for (i = 0; i < G_N_ELEMENTS (test_green_ells); i++)
    {
      NcmVector *x = ncm_sbessel_ode_solver_solve_dense (solver, TEST_A, TEST_B, test_green_ells[i], rhs_v, 128);
      GArray *u    = ncm_vector_dup_array (x);

      _test_green_check_coeffs (f, test_green_ells[i], u);

      g_array_unref (u);
      ncm_vector_free (x);
    }

    ncm_vector_free (rhs_v);
    g_array_unref (rhs);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* The adaptive QR solve, ell 0, 1, 5, 10, 20 */
void
test_ncm_sbessel_ode_solver_green_qr (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint f, i;

  for (f = TEST_GREEN_ONE; f <= TEST_GREEN_X; f++)
  {
    GArray *rhs = _test_green_rhs (f);

    for (i = 0; i < G_N_ELEMENTS (test_green_ells); i++)
    {
      GArray *u = NULL;
      gsize len;
      NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, test_green_ells[i], test_green_ells[i], rhs, &u, &len);

      _test_green_check_coeffs (f, test_green_ells[i], u);

      g_array_unref (u);
      ncm_sbessel_ode_operator_unref (op);
    }

    g_array_unref (rhs);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* solve_endpoints, ell 0, 1, 5, 10, 20 */
void
test_ncm_sbessel_ode_solver_green_endpoints (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint f, i;

  for (f = TEST_GREEN_ONE; f <= TEST_GREEN_X; f++)
  {
    GArray *rhs = _test_green_rhs (f);

    for (i = 0; i < G_N_ELEMENTS (test_green_ells); i++)
    {
      NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, test_green_ells[i], test_green_ells[i]);
      GArray *e                 = NULL;

      ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &e);
      _test_green_check (f, test_green_ells[i], g_array_index (e, gdouble, 0), g_array_index (e, gdouble, 1));

      g_array_unref (e);
      ncm_sbessel_ode_operator_unref (op);
    }

    g_array_unref (rhs);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* Batched solve_endpoints, ell 5..15 */
void
test_ncm_sbessel_ode_solver_green_endpoints_batched (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint f, ell;

  for (f = TEST_GREEN_ONE; f <= TEST_GREEN_X; f++)
  {
    GArray *rhs               = _test_green_rhs (f);
    NcmSBesselOdeOperator *op = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 15);
    GArray *e                 = NULL;

    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &e);

    for (ell = 5; ell <= 15; ell++)
      _test_green_check (f, ell, g_array_index (e, gdouble, 3 * (ell - 5) + 0), g_array_index (e, gdouble, 3 * (ell - 5) + 1));

    g_array_unref (e);
    g_array_unref (rhs);
    ncm_sbessel_ode_operator_unref (op);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * Each forcing of the certified table (x^2 and x^3 at ell 0, 5, 10 and 15 or 10; exp (x),
 * sin (x), cos (x), exp (-x^2) and 1 / (1 + x^2) at ell 0, 5, 10), expanded with
 * ncm_spectral_compute_chebyshev_coeffs() at order 126 and converted to C^(2), solved by
 * the adaptive QR for one multipole. The boundary tolerance is 1e-10 for the polynomials
 * and 1e-8 for the others.
 */
void
test_ncm_sbessel_ode_solver_forcing (void)
{
  const gchar *cases[]        = {"rhs_quadratic", "rhs_cubic", "rhs_exp", "rhs_sin", "rhs_cos", "rhs_gauss", "rhs_lorentz"};
  const guint ell_a[]         = {0, 5, 10, 15};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint i, l;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    const gboolean poly = (i < 2);
    guint n;
    const TestArbCell *c0 = test_arb_cells_find (test_arb_cells, G_N_ELEMENTS (test_arb_cells), cases[i], &n);
    GArray *rhs           = _test_rhs_spectral (128, &_test_cheb_cell, (gpointer) c0);

    for (l = 0; l < G_N_ELEMENTS (ell_a); l++)
    {
      const TestArbCell *c;
      NcmSBesselOdeOperator *op;
      GArray *u = NULL;
      gsize len;

      if ((ell_a[l] == 15) && (strcmp (cases[i], "rhs_quadratic") != 0))
        continue;

      c  = test_arb_cells_find_ell (test_arb_cells, G_N_ELEMENTS (test_arb_cells), cases[i], ell_a[l]);
      op = _test_solve (solver, TEST_A, TEST_B, c->ell, c->ell, rhs, &u, &len);

      _test_forcing_check (c, u, poly ? 1.0e-10 : 1.0e-8);

      g_array_unref (u);
      ncm_sbessel_ode_operator_unref (op);
    }

    g_array_unref (rhs);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* exp (x) and sin (x) through one batch of ell 2..8 */
void
test_ncm_sbessel_ode_solver_forcing_batched (void)
{
  const gchar *cases[]        = {"rhs_exp", "rhs_sin"};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  guint i, ell;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    guint n;
    const TestArbCell *c0 = test_arb_cells_find (test_arb_cells, G_N_ELEMENTS (test_arb_cells), cases[i], &n);
    GArray *rhs           = _test_rhs_spectral (128, &_test_cheb_cell, (gpointer) c0);
    GArray *sol           = NULL;
    gsize len;
    NcmSBesselOdeOperator *op = _test_solve (solver, TEST_A, TEST_B, 2, 8, rhs, &sol, &len);

    for (ell = 2; ell <= 8; ell++)
    {
      const TestArbCell *c = test_arb_cells_find_ell (test_arb_cells, G_N_ELEMENTS (test_arb_cells), cases[i], ell);
      GArray *u            = _test_row (sol, len, ell - 2);

      _test_forcing_check (c, u, 1.0e-8);

      g_array_unref (u);
    }

    g_array_unref (sol);
    g_array_unref (rhs);
    ncm_sbessel_ode_operator_unref (op);
  }

  ncm_mpsf_sbessel_free_cache ();
  ncm_sbessel_ode_solver_free (solver);
}

/* Nothing is stored before the first solve; ell 0..2 */
void
test_ncm_sbessel_ode_solver_memory_accessors (void)
{
  const guint orders[]        = {2};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 2);

  g_assert_cmpuint (ncm_sbessel_ode_operator_get_operator_size (op), ==, 0);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);

  _test_memory_sequence (0, 2, orders, NULL, 1, 64);
}

/* rhs[2] (length 32) then rhs[50] (length 128), ell 0 */
void
test_ncm_sbessel_ode_solver_memory_grow (void)
{
  const guint orders[]        = {2, 50};
  const guint lens[]          = {32, 128};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 0);
  gsize capacity              = 0;
  guint i;

  for (i = 0; i < 2; i++)
  {
    GArray *rhs = _test_rhs_unit (orders[i], lens[i], 1.0);
    gsize len;
    GArray *sol = _test_memory_solve (solver, op, rhs, &capacity, &len);

    g_array_unref (sol);
    g_array_unref (rhs);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* rhs[60] (length 128) then rhs[2] (length 32), ell 0 */
void
test_ncm_sbessel_ode_solver_memory_shrink (void)
{
  const guint orders[]        = {60, 2};
  const guint lens[]          = {128, 32};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 0);
  gsize capacity              = 0;
  guint i;

  for (i = 0; i < 2; i++)
  {
    GArray *rhs = _test_rhs_unit (orders[i], lens[i], 1.0);
    gsize len;
    GArray *sol = _test_memory_solve (solver, op, rhs, &capacity, &len);

    g_array_unref (sol);
    g_array_unref (rhs);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* The same right-hand side twice: capacity and n_cols unchanged, the solutions equal to 1e-14 */
void
test_ncm_sbessel_ode_solver_memory_identical (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 0);
  GArray *rhs                 = _test_rhs_unit (2, 64, 1.0);
  gsize capacity              = 0;
  gsize len1, len2, capacity1;
  glong n_cols1;
  GArray *sol1, *sol2;
  guint j;

  sol1      = _test_memory_solve (solver, op, rhs, &capacity, &len1);
  capacity1 = capacity;
  n_cols1   = ncm_sbessel_ode_operator_get_n_cols (op);
  sol2      = _test_memory_solve (solver, op, rhs, &capacity, &len2);

  g_assert_cmpuint (capacity, ==, capacity1);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols1);
  g_assert_cmpuint (len1, ==, len2);

  for (j = 0; j < len1; j++)
    _test_assert_close (g_array_index (sol2, gdouble, j), g_array_index (sol1, gdouble, j), 1.0e-14, 0.0);

  g_array_unref (sol1);
  g_array_unref (sol2);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* rhs[2] then rhs[30], both of length 64, ell 0 */
void
test_ncm_sbessel_ode_solver_memory_same_length (void)
{
  const guint orders[] = {2, 30};

  _test_memory_sequence (0, 0, orders, NULL, 2, 64);
}

/* Pairs (2, 10), (2, 30), (10, 30), (30, 50) at ell 5; the two solutions differ (numpy allclose defaults) */
void
test_ncm_sbessel_ode_solver_memory_increasing (void)
{
  const guint pairs[][2] = {
    {
      2, 10
    }, {
      2, 30
    }, {
      10, 30
    }, {
      30, 50
    }
  };
  guint p;

  for (p = 0; p < G_N_ELEMENTS (pairs); p++)
  {
    NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
    NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 5, 5);
    GArray *rhs1                = _test_rhs_unit (pairs[p][0], 64, 1.0);
    GArray *rhs2                = _test_rhs_unit (pairs[p][1], 64, 1.0);
    gsize capacity              = 0;
    gsize len1, len2;
    GArray *sol1    = _test_memory_solve (solver, op, rhs1, &capacity, &len1);
    GArray *sol2    = _test_memory_solve (solver, op, rhs2, &capacity, &len2);
    gboolean differ = FALSE;
    guint j;

    for (j = 0; j < GSL_MIN (len1, len2); j++)
      differ = differ || (fabs (g_array_index (sol1, gdouble, j) - g_array_index (sol2, gdouble, j)) > 1.0e-8 + 1.0e-5 * fabs (g_array_index (sol2, gdouble, j)));

    g_assert_true (differ);

    g_array_unref (sol1);
    g_array_unref (sol2);
    g_array_unref (rhs1);
    g_array_unref (rhs2);
    ncm_sbessel_ode_operator_unref (op);
    ncm_sbessel_ode_solver_free (solver);
  }
}

/* Order sequences (2, 10, 2), (2, 30, 10), (30, 10, 5), (5, 30, 5, 30), ell 0 */
void
test_ncm_sbessel_ode_solver_memory_varying (void)
{
  const guint s1[] = {2, 10, 2};
  const guint s2[] = {2, 30, 10};
  const guint s3[] = {30, 10, 5};
  const guint s4[] = {5, 30, 5, 30};

  _test_memory_sequence (0, 0, s1, NULL, G_N_ELEMENTS (s1), 64);
  _test_memory_sequence (0, 0, s2, NULL, G_N_ELEMENTS (s2), 64);
  _test_memory_sequence (0, 0, s3, NULL, G_N_ELEMENTS (s3), 64);
  _test_memory_sequence (0, 0, s4, NULL, G_N_ELEMENTS (s4), 64);
}

/* rhs[2] (length 64) then rhs[40] (length 128), batch ell 0..4 */
void
test_ncm_sbessel_ode_solver_memory_batched (void)
{
  const guint orders[]        = {2, 40};
  const guint lens[]          = {64, 128};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 4);
  gsize capacity              = 0;
  guint i;

  for (i = 0; i < 2; i++)
  {
    GArray *rhs = _test_rhs_unit (orders[i], lens[i], 1.0);
    gsize len;
    GArray *sol = _test_memory_solve (solver, op, rhs, &capacity, &len);

    g_array_unref (sol);
    g_array_unref (rhs);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* solve_endpoints with rhs[2] (length 32) then rhs[40] (length 128), ell 0: derivatives equal to a new operator's */
void
test_ncm_sbessel_ode_solver_memory_endpoints (void)
{
  const guint orders[]        = {2, 40};
  const guint lens[]          = {32, 128};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 0);
  gsize capacity              = 0;
  guint i, j;

  for (i = 0; i < 2; i++)
  {
    GArray *rhs                = _test_rhs_unit (orders[i], lens[i], 1.0);
    NcmSBesselOdeOperator *fop = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 0);
    GArray *e                  = NULL;
    GArray *fe                 = NULL;
    gdouble scale;

    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &e);
    g_assert_cmpuint (ncm_sbessel_ode_operator_get_operator_size (op), >=, capacity);
    capacity = ncm_sbessel_ode_operator_get_operator_size (op);
    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), >, 0);

    ncm_sbessel_ode_operator_solve_endpoints (fop, rhs, &fe);
    scale = GSL_MAX (fabs (g_array_index (fe, gdouble, 0)), fabs (g_array_index (fe, gdouble, 1)));

    for (j = 0; j < 2; j++)
      _test_assert_close (g_array_index (e, gdouble, j), g_array_index (fe, gdouble, j), 1.0e-12, 1.0e-12 * scale);

    g_array_unref (e);
    g_array_unref (fe);
    g_array_unref (rhs);
    ncm_sbessel_ode_operator_unref (fop);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* reconfigure_operator (to ell 1) drops the stored factorization */
void
test_ncm_sbessel_ode_solver_memory_reset (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_A, TEST_B, 0, 0);
  GArray *rhs                 = _test_rhs_unit (2, 64, 1.0);
  gsize capacity              = 0;
  gsize len;
  GArray *sol = _test_memory_solve (solver, op, rhs, &capacity, &len);

  ncm_sbessel_ode_solver_reconfigure_operator (solver, op, TEST_A, TEST_B, 1, 1);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, 0);

  g_array_unref (sol);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* Twenty solves, rhs[order_i] = 1 / (i + 1), the orders the Python test drew with numpy seed 21 */
void
test_ncm_sbessel_ode_solver_memory_stress (void)
{
  const guint orders[] = {5, 2, 20, 2, 2, 2, 15, 20, 10, 5, 10, 2, 400, 20, 300, 300, 5, 300, 400, 10};
  gdouble values[G_N_ELEMENTS (orders)];
  guint i;

  for (i = 0; i < G_N_ELEMENTS (orders); i++)
    values[i] = 1.0 / (i + 1.0);

  _test_memory_sequence (0, 0, orders, values, G_N_ELEMENTS (orders), 64);
}

/* The same right-hand side twice on [1, 10], ell 5: same length, n_cols and coefficients to 1e-15 */
void
test_ncm_sbessel_ode_solver_rotations_identical (void)
{
  const TestCoeff c[] = {
    {
      10, 1.0
    }, {
      15, 0.5
    }, {
      20, 0.25
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 10.0, 5, 5);
  GArray *rhs                 = _test_rhs_sparse (50, c, G_N_ELEMENTS (c));
  GArray *sol1                = NULL;
  GArray *sol2                = NULL;
  gsize len1, len2;
  glong n_cols1;
  guint j;

  ncm_sbessel_ode_operator_solve (op, rhs, &sol1, &len1);
  n_cols1 = ncm_sbessel_ode_operator_get_n_cols (op);
  ncm_sbessel_ode_operator_solve (op, rhs, &sol2, &len2);

  g_assert_cmpuint (len1, ==, len2);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols1);

  for (j = 0; j < len1; j++)
    _test_assert_close (g_array_index (sol2, gdouble, j), g_array_index (sol1, gdouble, j), 1.0e-15, 1.0e-15);

  g_array_unref (sol1);
  g_array_unref (sol2);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* A short then a longer right-hand side on [1, 10], ell 10: the extension converges at the new operator's column */
void
test_ncm_sbessel_ode_solver_rotations_extend (void)
{
  const TestCoeff c1[] = {
    {
      10, 1.0
    }
  };
  const TestCoeff c2[] = {
    {
      10, 1.0
    }, {
      15, 0.5
    }, {
      50, 0.1
    }, {
      60, 0.05
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 10.0, 10, 10);
  GArray *rhs1                = _test_rhs_sparse (40, c1, G_N_ELEMENTS (c1));
  GArray *rhs2                = _test_rhs_sparse (80, c2, G_N_ELEMENTS (c2));
  GArray *sol                 = NULL;
  gsize len;
  glong n_cols_new;

  ncm_sbessel_ode_operator_solve (op, rhs1, &sol, &len);
  g_array_unref (sol);

  sol = _test_reused_vs_new (solver, op, rhs2, 1.0e-12, &len, &n_cols_new);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols_new);
  g_assert_cmpfloat (_test_absmax (sol), >, 1.0e-10);

  g_array_unref (sol);
  g_array_unref (rhs1);
  g_array_unref (rhs2);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* Five single-coefficient right-hand sides in sequence on [1, 10], ell 8; consecutive solutions differ by more than 1e-8 */
void
test_ncm_sbessel_ode_solver_rotations_sequence (void)
{
  const guint lens[]  = {50, 50, 50, 60, 45};
  const TestCoeff c[] = {
    {
      10, 1.0
    }, {
      20, 0.8
    }, {
      15, 1.5
    }, {
      25, 1.2
    }, {
      12, 0.9
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 10.0, 8, 8);
  GArray *prev                = NULL;
  gsize prev_len              = 0;
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (c); i++)
  {
    GArray *rhs = _test_rhs_sparse (lens[i], &c[i], 1);
    gsize len;
    GArray *sol = _test_reused_vs_new (solver, op, rhs, 1.0e-12, &len, NULL);

    if (prev != NULL)
    {
      gdouble diff = 0.0;

      for (j = 0; j < GSL_MIN (len, prev_len); j++)
        diff = GSL_MAX (diff, fabs (g_array_index (sol, gdouble, j) - g_array_index (prev, gdouble, j)));

      g_assert_cmpfloat (diff, >, 1.0e-8);
      g_array_unref (prev);
    }

    prev     = sol;
    prev_len = len;
    g_array_unref (rhs);
  }

  g_array_unref (prev);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* The extension in a batch of ell 5..8 on [1, 10]; every multipole non-zero */
void
test_ncm_sbessel_ode_solver_rotations_batched_extend (void)
{
  const TestCoeff c1[] = {
    {
      12, 1.0
    }, {
      15, 0.3
    }
  };
  const TestCoeff c2[] = {
    {
      12, 1.0
    }, {
      15, 0.3
    }, {
      35, 0.15
    }, {
      50, 0.08
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 10.0, 5, 8);
  GArray *rhs1                = _test_rhs_sparse (40, c1, G_N_ELEMENTS (c1));
  GArray *rhs2                = _test_rhs_sparse (70, c2, G_N_ELEMENTS (c2));
  GArray *sol                 = NULL;
  gsize len;
  glong n_cols_new;
  guint l;

  ncm_sbessel_ode_operator_solve (op, rhs1, &sol, &len);
  g_array_unref (sol);

  sol = _test_reused_vs_new (solver, op, rhs2, 1.0e-12, &len, &n_cols_new);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols_new);

  for (l = 0; l < 4; l++)
    g_assert_cmpfloat (_test_row_absmax (sol, len, l), >, 1.0e-10);

  g_array_unref (sol);
  g_array_unref (rhs1);
  g_array_unref (rhs2);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* rhs[15] then only rhs[40] (leading zeros past the stored rotations) on [1, 10], ell 5: non-zero and equal to new */
void
test_ncm_sbessel_ode_solver_rotations_leading_zeros (void)
{
  const TestCoeff c1[] = {
    {
      15, 1.0
    }
  };
  const TestCoeff c2[] = {
    {
      40, 0.5
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 10.0, 5, 5);
  GArray *rhs1                = _test_rhs_sparse (50, c1, G_N_ELEMENTS (c1));
  GArray *rhs2                = _test_rhs_sparse (60, c2, G_N_ELEMENTS (c2));
  GArray *sol                 = NULL;
  gsize len;

  ncm_sbessel_ode_operator_solve (op, rhs1, &sol, &len);
  g_array_unref (sol);

  sol = _test_reused_vs_new (solver, op, rhs2, 1.0e-12, &len, NULL);
  g_assert_cmpfloat (_test_absmax (sol), >, 1.0e-10);

  g_array_unref (sol);
  g_array_unref (rhs1);
  g_array_unref (rhs2);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * A slightly longer second right-hand side on [1, 15], ell 12: the solution (length
 * included) equals the new operator's. get_n_cols counts the stored columns, which the
 * reused operator keeps from the first solve, so it is at least the new operator's.
 */
void
test_ncm_sbessel_ode_solver_rotations_convergence (void)
{
  const TestCoeff c1[] = {
    {
      18, 1.0
    }, {
      22, 0.1
    }
  };
  const TestCoeff c2[] = {
    {
      18, 1.0
    }, {
      22, 0.1
    }, {
      35, 0.02
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 15.0, 12, 12);
  GArray *rhs1                = _test_rhs_sparse (60, c1, G_N_ELEMENTS (c1));
  GArray *rhs2                = _test_rhs_sparse (70, c2, G_N_ELEMENTS (c2));
  GArray *sol                 = NULL;
  gsize len;
  glong n_cols_new;

  ncm_sbessel_ode_operator_solve (op, rhs1, &sol, &len);
  g_array_unref (sol);

  sol = _test_reused_vs_new (solver, op, rhs2, 1.0e-12, &len, &n_cols_new);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), >=, n_cols_new);

  g_array_unref (sol);
  g_array_unref (rhs1);
  g_array_unref (rhs2);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* Four right-hand sides of length 40, 50, 65 and 75 through one batch of ell 3..6 on [1, 12]; n_cols never decreases */
void
test_ncm_sbessel_ode_solver_rotations_batched_sequence (void)
{
  const guint lens[]          = {40, 50, 65, 75};
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 12.0, 3, 6);
  glong prev_n_cols           = 0;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (lens); i++)
  {
    const TestCoeff c[] = {
      {
        15, 1.0
      }, {
        20, 0.4
      }, {
        35, 0.15
      }, {
        55, 0.05
      }
    };
    const guint n_c = (lens[i] > 60) ? 4 : ((lens[i] > 40) ? 3 : 2);
    GArray *rhs     = _test_rhs_sparse (lens[i], c, n_c);
    gsize len;
    GArray *sol = _test_reused_vs_new (solver, op, rhs, 1.0e-12, &len, NULL);

    g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), >=, prev_n_cols);
    prev_n_cols = ncm_sbessel_ode_operator_get_n_cols (op);

    g_array_unref (sol);
    g_array_unref (rhs);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* ell 50 on [1, 20], 80 then 150 coefficients: equal to new to 1e-11, at least 20 columns added */
void
test_ncm_sbessel_ode_solver_rotations_high_ell (void)
{
  const TestCoeff c1[] = {
    {
      25, 1.0
    }, {
      30, 0.2
    }
  };
  const TestCoeff c2[] = {
    {
      25, 1.0
    }, {
      30, 0.2
    }, {
      60, 0.1
    }, {
      90, 0.05
    }
  };
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, 1.0, 20.0, 50, 50);
  GArray *rhs1                = _test_rhs_sparse (80, c1, G_N_ELEMENTS (c1));
  GArray *rhs2                = _test_rhs_sparse (150, c2, G_N_ELEMENTS (c2));
  GArray *sol                 = NULL;
  gsize len;
  glong n_cols_short, n_cols_new;

  ncm_sbessel_ode_operator_solve (op, rhs1, &sol, &len);
  n_cols_short = ncm_sbessel_ode_operator_get_n_cols (op);
  g_array_unref (sol);

  sol = _test_reused_vs_new (solver, op, rhs2, 1.0e-11, &len, &n_cols_new);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op), ==, n_cols_new);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_n_cols (op) - n_cols_short, >=, 20);

  g_array_unref (sol);
  g_array_unref (rhs1);
  g_array_unref (rhs2);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * Deep oscillatory panel [1e4, 10^4.5]: the same integral to 1e-9 with far fewer
 * columns, at most twice the forcing order and under 1/20 of Dirichlet's; the
 * operators report their constraint and column floor.
 */
void
test_ncm_sbessel_ode_solver_tau_deep_panel (void)
{
  const TestBump p = _test_tau_bump (TEST_TAU_A, TEST_TAU_B);
  GArray *rhs      = _test_tau_rhs (&p);
  GArray *u_d      = NULL;
  GArray *u_t      = NULL;
  gsize n_d, n_t;
  NcmSBesselOdeOperator *op_d = _test_tau_solve (&p, TEST_TAU_ELL, TEST_TAU_ELL, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET, &u_d, &n_d);
  NcmSBesselOdeOperator *op_t = _test_tau_solve (&p, TEST_TAU_ELL, TEST_TAU_ELL, NCM_SBESSEL_ODE_CONSTRAINT_TAU, &u_t, &n_t);
  const guint n_forcing       = rhs->len - 2;

  _test_assert_close (_test_boundary_functional (u_t, p.a, p.b, TEST_TAU_ELL), _test_boundary_functional (u_d, p.a, p.b, TEST_TAU_ELL), 1.0e-9, 0.0);
  g_assert_cmpuint (n_t, <=, 2 * n_forcing);
  g_assert_cmpuint (20 * n_t, <, n_d);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op_t), ==, 0);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op_d), >, 0);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_constraint (op_d), ==, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_constraint (op_t), ==, NCM_SBESSEL_ODE_CONSTRAINT_TAU);

  g_array_unref (rhs);
  g_array_unref (u_d);
  g_array_unref (u_t);
  ncm_sbessel_ode_operator_unref (op_d);
  ncm_sbessel_ode_operator_unref (op_t);
  ncm_mpsf_sbessel_free_cache ();
}

/* Narrow panel [100, 10^2.5]: the tau solve's integral against the certified tau_narrow cell, rtol 1e-9 */
void
test_ncm_sbessel_ode_solver_tau_narrow_panel (void)
{
  const TestArbCell *c = test_arb_cells_find_ell (test_arb_cells, G_N_ELEMENTS (test_arb_cells), "tau_narrow", TEST_TAU_ELL);
  const TestBump p     = {c->a, c->b, c->p1, c->p2};
  GArray *u            = NULL;
  gsize n;
  NcmSBesselOdeOperator *op = _test_tau_solve (&p, TEST_TAU_ELL, TEST_TAU_ELL, NCM_SBESSEL_ODE_CONSTRAINT_TAU, &u, &n);

  _test_assert_close (_test_boundary_functional (u, p.a, p.b, TEST_TAU_ELL), c->value, 1.0e-9, 0.0);

  g_array_unref (u);
  ncm_sbessel_ode_operator_unref (op);
  ncm_mpsf_sbessel_free_cache ();
}

/* A tau block of ell 20..27 against single solves: rtol 1e-9, atol 1e-13 max |single| */
void
test_ncm_sbessel_ode_solver_tau_batched (void)
{
  const TestBump p = _test_tau_bump (TEST_TAU_A, TEST_TAU_B);
  GArray *u_b      = NULL;
  gsize n_b;
  NcmSBesselOdeOperator *op = _test_tau_solve (&p, TEST_TAU_ELL, TEST_TAU_ELL + 7, NCM_SBESSEL_ODE_CONSTRAINT_TAU, &u_b, &n_b);
  guint l, j;

  for (l = 0; l < 8; l++)
  {
    GArray *u_s = NULL;
    gsize n_s;
    NcmSBesselOdeOperator *sop = _test_tau_solve (&p, TEST_TAU_ELL + l, TEST_TAU_ELL + l, NCM_SBESSEL_ODE_CONSTRAINT_TAU, &u_s, &n_s);
    const gdouble scale        = _test_absmax (u_s);

    for (j = 0; j < n_b; j++)
    {
      const gdouble s = (j < n_s) ? g_array_index (u_s, gdouble, j) : 0.0;

      _test_assert_close (g_array_index (u_b, gdouble, l * n_b + j), s, 1.0e-9, 1.0e-13 * scale);
    }

    g_array_unref (u_s);
    ncm_sbessel_ode_operator_unref (sop);
  }

  g_array_unref (u_b);
  ncm_sbessel_ode_operator_unref (op);
}

/*
 * On [3.162, 10] at ell 2 the panel holds about 4 oscillations against a working order
 * of about 32, so the tau solve admits homogeneous content (coefficient ratio above
 * 1e6); on [100, 10^2.5] it does not (below 10).
 */
void
test_ncm_sbessel_ode_solver_tau_shallow_panel (void)
{
  g_assert_cmpfloat (_test_tau_contamination (3.162, 10.0), >, 1.0e6);
  g_assert_cmpfloat (_test_tau_contamination (100.0, 316.22776601683796), <, 10.0);
}

/* get_last_max_coeff is the solution's own max |a_j|, for one multipole and a block of 8 */
void
test_ncm_sbessel_ode_solver_tau_last_max_coeff (void)
{
  const guint n_ell_a[] = {1, 8};
  const TestBump p      = _test_tau_bump (TEST_TAU_A, TEST_TAU_B);
  guint c, l;

  for (c = 0; c < G_N_ELEMENTS (n_ell_a); c++)
  {
    GArray *u = NULL;
    gsize n;
    NcmSBesselOdeOperator *op = _test_tau_solve (&p, TEST_TAU_ELL, TEST_TAU_ELL + n_ell_a[c] - 1, NCM_SBESSEL_ODE_CONSTRAINT_TAU, &u, &n);

    for (l = 0; l < n_ell_a[c]; l++)
      g_assert_cmpfloat (ncm_sbessel_ode_operator_get_last_max_coeff (op, l), ==, _test_row_absmax (u, n, l));

    g_array_unref (u);
    ncm_sbessel_ode_operator_unref (op);
  }
}

/* The same after a scalar solve and a reconfigure to a block of 8 (storage sized for one multipole first) */
void
test_ncm_sbessel_ode_solver_tau_last_max_coeff_reconfigure (void)
{
  const TestBump p            = _test_tau_bump (TEST_TAU_A, TEST_TAU_B);
  NcmSBesselOdeSolver *solver = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_TAU);
  GArray *rhs                 = _test_tau_rhs (&p);
  GArray *u                   = NULL;
  gsize n;
  NcmSBesselOdeOperator *op = _test_solve (solver, p.a, p.b, TEST_TAU_ELL, TEST_TAU_ELL, rhs, &u, &n);
  guint l;

  ncm_sbessel_ode_solver_reconfigure_operator (solver, op, p.a, p.b, TEST_TAU_ELL, TEST_TAU_ELL + 7);
  ncm_sbessel_ode_operator_solve (op, rhs, &u, &n);

  for (l = 0; l < 8; l++)
    g_assert_cmpfloat (ncm_sbessel_ode_operator_get_last_max_coeff (op, l), ==, _test_row_absmax (u, n, l));

  g_array_unref (u);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/*
 * Pinning coefficients 120 and 144 selects another member with the same integral (rtol
 * 1e-9) at exactly 145 columns; returning to Dirichlet restores its column floor.
 */
void
test_ncm_sbessel_ode_solver_tau_pinned (void)
{
  const TestBump p            = _test_tau_bump (TEST_TAU_A, TEST_TAU_B);
  NcmSBesselOdeSolver *solver = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  GArray *rhs                 = _test_tau_rhs (&p);
  GArray *u_d                 = NULL;
  GArray *u_p                 = NULL;
  gsize n_d, n_p;
  NcmSBesselOdeOperator *op = _test_solve (solver, p.a, p.b, TEST_TAU_ELL, TEST_TAU_ELL, rhs, &u_d, &n_d);
  const glong floor_cols    = ncm_sbessel_ode_operator_get_min_cols (op);
  glong pin1, pin2;

  ncm_sbessel_ode_operator_set_pinned_constraint (op, 120, 144);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_constraint (op), ==, NCM_SBESSEL_ODE_CONSTRAINT_PINNED);
  ncm_sbessel_ode_operator_get_pins (op, &pin1, &pin2);
  g_assert_cmpint (pin1, ==, 120);
  g_assert_cmpint (pin2, ==, 144);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op), ==, 145);

  ncm_sbessel_ode_operator_solve (op, rhs, &u_p, &n_p);
  g_assert_cmpuint (n_p, ==, 145);
  _test_assert_close (_test_boundary_functional (u_p, p.a, p.b, TEST_TAU_ELL), _test_boundary_functional (u_d, p.a, p.b, TEST_TAU_ELL), 1.0e-9, 0.0);

  ncm_sbessel_ode_operator_set_constraint (op, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_constraint (op), ==, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op), ==, floor_cols);

  g_array_unref (u_d);
  g_array_unref (u_p);
  g_array_unref (rhs);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
  ncm_mpsf_sbessel_free_cache ();
}

/* The solver default is Dirichlet; operators inherit the solver's constraint; per-operator setters override it */
void
test_ncm_sbessel_ode_solver_tau_propagation (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op;

  g_assert_cmpint (ncm_sbessel_ode_solver_get_default_constraint (solver), ==, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);

  ncm_sbessel_ode_solver_set_default_constraint (solver, NCM_SBESSEL_ODE_CONSTRAINT_TAU);
  op = ncm_sbessel_ode_solver_create_operator (solver, TEST_TAU_A, TEST_TAU_B, TEST_TAU_ELL, TEST_TAU_ELL);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_constraint (op), ==, NCM_SBESSEL_ODE_CONSTRAINT_TAU);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op), ==, 0);

  ncm_sbessel_ode_operator_set_constraint (op, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_constraint (op), ==, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op), >, 0);

  ncm_sbessel_ode_operator_set_min_cols (op, 7);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op), ==, 7);
  ncm_sbessel_ode_operator_set_min_cols (op, -1);
  g_assert_cmpint (ncm_sbessel_ode_operator_get_min_cols (op), >, 7);

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* The floor factor ships above one and reads back */
void
test_ncm_sbessel_ode_solver_floor_default (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();

  g_assert_cmpfloat (ncm_sbessel_ode_solver_get_tau_floor_factor (solver), >, 1.0);
  ncm_sbessel_ode_solver_set_tau_floor_factor (solver, 1.4);
  g_assert_cmpfloat (ncm_sbessel_ode_solver_get_tau_floor_factor (solver), ==, 1.4);

  ncm_sbessel_ode_solver_free (solver);
}

/* The tau solve keeps ceil (factor n_forcing) columns, factor 1.1, 1.5 and 2.0, as get_tau_constraint_order reports */
void
test_ncm_sbessel_ode_solver_floor_working_order (void)
{
  const gdouble factor_a[] = {1.1, 1.5, 2.0};
  guint n_forcing, i;
  GArray *rhs = _test_floor_rhs_default (&n_forcing);

  for (i = 0; i < G_N_ELEMENTS (factor_a); i++)
  {
    NcmSBesselOdeSolver *solver = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_TAU);
    NcmSBesselOdeOperator *op;
    GArray *u = NULL;
    gsize n;

    ncm_sbessel_ode_solver_set_tau_floor_factor (solver, factor_a[i]);
    op = _test_solve (solver, TEST_TAU_A, TEST_TAU_B, TEST_TAU_ELL, TEST_TAU_ELL, rhs, &u, &n);

    g_assert_cmpuint (n, ==, (gsize) ceil (factor_a[i] * n_forcing));
    g_assert_cmpint (ncm_sbessel_ode_operator_get_tau_constraint_order (op, rhs->len), ==, (glong) n);

    g_array_unref (u);
    ncm_sbessel_ode_operator_unref (op);
    ncm_sbessel_ode_solver_free (solver);
  }

  g_array_unref (rhs);
}

/* get_last_deriv_error is sum_j j^2 |a_j| / h after solve, solve_values and solve_endpoints, Dirichlet and tau */
void
test_ncm_sbessel_ode_solver_floor_deriv_error (void)
{
  const NcmSBesselOdeConstraint cons[] = {NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET, NCM_SBESSEL_ODE_CONSTRAINT_TAU};
  const gdouble h                      = 0.5 * (TEST_TAU_B - TEST_TAU_A);
  guint n_forcing, i;
  GArray *rhs = _test_floor_rhs_default (&n_forcing);

  for (i = 0; i < G_N_ELEMENTS (cons); i++)
  {
    NcmSBesselOdeSolver *solver = _test_tau_solver (cons[i]);
    GArray *u                   = NULL;
    GArray *values              = NULL;
    GArray *endpoints           = NULL;
    gsize n;
    NcmSBesselOdeOperator *op = _test_solve (solver, TEST_TAU_A, TEST_TAU_B, TEST_TAU_ELL, TEST_TAU_ELL, rhs, &u, &n);
    const gdouble expected    = _test_weighted_sum (u, n, 0, h);

    _test_assert_close (ncm_sbessel_ode_operator_get_last_deriv_error (op, 0), expected, 1.0e-13, 0.0);

    ncm_sbessel_ode_operator_solve_values (op, rhs, TEST_TAU_A, TEST_TAU_B, &values);
    _test_assert_close (ncm_sbessel_ode_operator_get_last_deriv_error (op, 0), expected, 1.0e-13, 0.0);

    ncm_sbessel_ode_operator_solve_endpoints (op, rhs, &endpoints);
    _test_assert_close (ncm_sbessel_ode_operator_get_last_deriv_error (op, 0), g_array_index (endpoints, gdouble, 2), 1.0e-13, 0.0);

    g_array_unref (u);
    g_array_unref (values);
    g_array_unref (endpoints);
    ncm_sbessel_ode_operator_unref (op);
    ncm_sbessel_ode_solver_free (solver);
  }

  g_array_unref (rhs);
}

/* Each member of a tau block of 1, 4 and 8 multipoles reports its own bound */
void
test_ncm_sbessel_ode_solver_floor_deriv_error_per_ell (void)
{
  const guint n_ell_a[] = {1, 4, 8};
  const gdouble h       = 0.5 * (TEST_TAU_B - TEST_TAU_A);
  guint n_forcing, c, l;
  GArray *rhs = _test_floor_rhs_default (&n_forcing);

  for (c = 0; c < G_N_ELEMENTS (n_ell_a); c++)
  {
    NcmSBesselOdeSolver *solver = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_TAU);
    GArray *u                   = NULL;
    gsize n;
    NcmSBesselOdeOperator *op = _test_solve (solver, TEST_TAU_A, TEST_TAU_B, TEST_TAU_ELL, TEST_TAU_ELL + n_ell_a[c] - 1, rhs, &u, &n);

    for (l = 0; l < n_ell_a[c]; l++)
      _test_assert_close (ncm_sbessel_ode_operator_get_last_deriv_error (op, l), _test_weighted_sum (u, n, l, h), 1.0e-13, 0.0);

    g_array_unref (u);
    ncm_sbessel_ode_operator_unref (op);
    ncm_sbessel_ode_solver_free (solver);
  }

  g_array_unref (rhs);
}

/* The 40-row dense matrix carries the operator's constraint rows: Dirichlet full, pinned one entry each, tau none */
void
test_ncm_sbessel_ode_solver_floor_matrix_constraint (void)
{
  const guint nrows                    = 40;
  NcmSBesselOdeSolver *solver          = _test_tau_solver (NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  NcmSBesselOdeOperator *op            = ncm_sbessel_ode_solver_create_operator (solver, TEST_TAU_A, TEST_TAU_B, TEST_TAU_ELL, TEST_TAU_ELL);
  const NcmSBesselOdeConstraint cons[] = {NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET, NCM_SBESSEL_ODE_CONSTRAINT_PINNED, NCM_SBESSEL_ODE_CONSTRAINT_TAU};
  const guint expected[]               = {nrows, 1, 0};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cons); i++)
  {
    NcmMatrix *m;

    if (cons[i] == NCM_SBESSEL_ODE_CONSTRAINT_PINNED)
      ncm_sbessel_ode_operator_set_pinned_constraint (op, 10, 20);
    else
      ncm_sbessel_ode_operator_set_constraint (op, cons[i]);

    m = ncm_sbessel_ode_operator_get_matrix (op, nrows);
    g_assert_cmpuint (_test_row_nonzeros (m, 0), ==, expected[i]);
    g_assert_cmpuint (_test_row_nonzeros (m, 1), ==, expected[i]);
    ncm_matrix_free (m);
  }

  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* Below the two constraint rows the Dirichlet and tau matrices are bitwise equal */
void
test_ncm_sbessel_ode_solver_floor_matrix_band (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeOperator *op   = ncm_sbessel_ode_solver_create_operator (solver, TEST_TAU_A, TEST_TAU_B, TEST_TAU_ELL, TEST_TAU_ELL);
  NcmMatrix *m_d, *m_t;
  guint i, j;

  ncm_sbessel_ode_operator_set_constraint (op, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET);
  m_d = ncm_sbessel_ode_operator_get_matrix (op, 40);
  ncm_sbessel_ode_operator_set_constraint (op, NCM_SBESSEL_ODE_CONSTRAINT_TAU);
  m_t = ncm_sbessel_ode_operator_get_matrix (op, 40);

  for (i = 2; i < 40; i++)
    for (j = 0; j < 40; j++)
      g_assert_cmpfloat (ncm_matrix_get (m_t, i, j), ==, ncm_matrix_get (m_d, i, j));

  ncm_matrix_free (m_d);
  ncm_matrix_free (m_t);
  ncm_sbessel_ode_operator_unref (op);
  ncm_sbessel_ode_solver_free (solver);
}

/* Spans 30, 60, 100, 150, 250 (the constraint fails) and 3000 (sound) */
void
test_ncm_sbessel_ode_solver_deriv_bound_never_below (void)
{
  const guint span_a[] = {30, 60, 100, 150, 250, 3000};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (span_a); i++)
  {
    const TestArbCell *c    = _test_deriv_bound_cell (span_a[i]);
    const gdouble amplitude = GSL_MAX (fabs (c->a * ncm_sf_sbessel (TEST_TAU_ELL, c->a)), fabs (c->b * ncm_sf_sbessel (TEST_TAU_ELL, c->b)));
    gdouble value, deriv_error, predicted, realized;

    _test_deriv_bound_solve (c, &value, &deriv_error);
    predicted = GSL_DBL_EPSILON * deriv_error * amplitude / fabs (c->value);
    realized  = fabs (value - c->value) / fabs (c->value);

    if (realized > 1.0e-9)
      g_assert_cmpfloat (predicted, >=, realized);
  }

  ncm_mpsf_sbessel_free_cache ();
}

/* A contaminated panel (span 30) reports a bound more than 1e6 times that of a sound one (span 3000) */
void
test_ncm_sbessel_ode_solver_deriv_bound_separates (void)
{
  gdouble value, deriv_sound, deriv_broken;

  _test_deriv_bound_solve (_test_deriv_bound_cell (3000), &value, &deriv_sound);
  _test_deriv_bound_solve (_test_deriv_bound_cell (30), &value, &deriv_broken);

  g_assert_cmpfloat (deriv_broken, >, 1.0e6 * deriv_sound);

  ncm_mpsf_sbessel_free_cache ();
}

/* Dirichlet on the conjugate span is wrong by more than 1e-6 */
void
test_ncm_sbessel_ode_solver_conjugate_dirichlet_fails (void)
{
  const TestArbCell *c = _test_conjugate_cell ("conjugate");
  const gdouble value  = _test_conjugate_integral (c, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET, 0, 0);

  g_assert_cmpfloat (fabs (value - c->value) / fabs (c->value), >, 1.0e-6);

  ncm_mpsf_sbessel_free_cache ();
}

/* Pinned at (0, 1), (2, 3) and (10, 11) on the same span: the certified value to 1e-9 */
void
test_ncm_sbessel_ode_solver_conjugate_pinned (void)
{
  const glong pins[][2] = {
    {
      0, 1
    }, {
      2, 3
    }, {
      10, 11
    }
  };
  const TestArbCell *c = _test_conjugate_cell ("conjugate");
  guint i;

  for (i = 0; i < G_N_ELEMENTS (pins); i++)
    _test_assert_close (_test_conjugate_integral (c, NCM_SBESSEL_ODE_CONSTRAINT_PINNED, pins[i][0], pins[i][1]), c->value, 1.0e-9, 0.0);

  ncm_mpsf_sbessel_free_cache ();
}

/* Dirichlet on [200, 200 + s + 0.01], off the conjugate span: the certified value to 1e-9 */
void
test_ncm_sbessel_ode_solver_conjugate_away (void)
{
  const TestArbCell *c = _test_conjugate_cell ("conjugate_away");

  _test_assert_close (_test_conjugate_integral (c, NCM_SBESSEL_ODE_CONSTRAINT_DIRICHLET, 0, 0), c->value, 1.0e-9, 0.0);

  ncm_mpsf_sbessel_free_cache ();
}


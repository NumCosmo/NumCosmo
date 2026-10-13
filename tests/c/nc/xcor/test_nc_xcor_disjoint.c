/***************************************************************************
 *            test_nc_xcor_disjoint.c
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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
 * Cross spectra of two kernels with disjoint radial support. Disjoint shells still
 * correlate: they are windows on the same 3D field, and the non-Limber kernel-space
 * method couples them through the outer k integral, which is the super-sample
 * covariance use case. Only the Limber tiers vanish, because there each multipole of a
 * kernel lives where (l + 1/2) / k is inside its own radial range, and disjoint shells
 * then have disjoint support in k. The tier, chosen per multipole, decides; the method
 * does not.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#include "test_nc_xcor_kquad_common.h"

/* Two shells in comoving distance [Mpc], far enough apart that nothing overlaps */
#define TEST_CHI1_LOWER 200.0
#define TEST_CHI1_UPPER 400.0
#define TEST_CHI2_LOWER 1500.0
#define TEST_CHI2_UPPER 1800.0

typedef struct _TestNcXcorDisjoint
{
  TestNcXcorKQuadEnv env;
  NcXcorKernel *k1;
  NcXcorKernel *k2;
} TestNcXcorDisjoint;

static void
test_nc_xcor_disjoint_new (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  test_nc_xcor_kquad_env_init (&test->env);
  test->k1 = NULL;
  test->k2 = NULL;
}

static void
test_nc_xcor_disjoint_free (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  nc_xcor_kernel_clear (&test->k1);
  nc_xcor_kernel_clear (&test->k2);
  test_nc_xcor_kquad_env_clear (&test->env);
}

/* The two shells, both in the given tier. */
static void
_test_nc_xcor_disjoint_kernels (TestNcXcorDisjoint *test, gint l_limber)
{
  test->k1 = test_nc_xcor_kquad_tophat (&test->env, TEST_CHI1_LOWER, TEST_CHI1_UPPER, l_limber);
  test->k2 = test_nc_xcor_kquad_tophat (&test->env, TEST_CHI2_LOWER, TEST_CHI2_UPPER, l_limber);
}

static NcXcor *
_test_nc_xcor_disjoint_xcor (TestNcXcorDisjoint *test, NcXcorMethod meth)
{
  NcXcor *xc = nc_xcor_new (test->env.dist, test->env.ps, meth);

  nc_xcor_prepare (xc, test->env.cosmo);

  return xc;
}

/*
 * Non-Limber, ell 0 to 2: finite and nonzero, bounded by the geometric mean of the
 * autos (Cauchy-Schwarz), and negative at ell = 0 for shells that do not touch. The
 * sign is set at the lowest k, so the kernels are converged: capped ones truncate that
 * end of the domain and get it wrong.
 */
static void
test_nc_xcor_disjoint_cross_nonzero (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  NcXcor *xc       = _test_nc_xcor_disjoint_xcor (test, NC_XCOR_METHOD_KERNEL_EXACT);
  NcmVector *cross = ncm_vector_new (3);
  NcmVector *auto1 = ncm_vector_new (3);
  NcmVector *auto2 = ncm_vector_new (3);
  guint i;

  test->k1 = test_nc_xcor_kquad_tophat_converged (&test->env, test->env.sbi, TEST_CHI1_LOWER, TEST_CHI1_UPPER, -1, 0.0, 0.0);
  test->k2 = test_nc_xcor_kquad_tophat_converged (&test->env, test->env.sbi, TEST_CHI2_LOWER, TEST_CHI2_UPPER, -1, 0.0, 0.0);

  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, 0, 2, cross);
  nc_xcor_compute (xc, test->k1, test->k1, test->env.cosmo, 0, 2, auto1);
  nc_xcor_compute (xc, test->k2, test->k2, test->env.cosmo, 0, 2, auto2);

  for (i = 0; i < 3; i++)
  {
    const gdouble c = ncm_vector_get (cross, i);

    g_assert_true (gsl_finite (c));
    g_assert_cmpfloat (c, !=, 0.0);
    g_assert_cmpfloat (fabs (c), <, sqrt (ncm_vector_get (auto1, i) * ncm_vector_get (auto2, i)));
  }

  g_assert_cmpfloat (ncm_vector_get (cross, 0), <, 0.0);

  ncm_vector_free (cross);
  ncm_vector_free (auto1);
  ncm_vector_free (auto2);
  nc_xcor_free (xc);
}

/* The solver reproduces nc_xcor_compute() for the non-Limber cross, ell 0 to 3. */
static void
test_nc_xcor_disjoint_solver (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  NcXcor *xc           = _test_nc_xcor_disjoint_xcor (test, NC_XCOR_METHOD_KERNEL_EXACT);
  NcXcorSolver *solver = nc_xcor_solver_new ();
  NcmVector *direct    = ncm_vector_new (4);
  NcmVector *solved;
  guint id1, id2, i;

  _test_nc_xcor_disjoint_kernels (test, -1);

  id1 = nc_xcor_solver_register_kernel (solver, test->k1);
  id2 = nc_xcor_solver_register_kernel (solver, test->k2);
  nc_xcor_solver_request_cl (solver, id1, id2, 0, 3);
  nc_xcor_solver_plan_blocks (solver, 8);
  nc_xcor_solver_solve (solver, xc, test->env.cosmo);
  solved = nc_xcor_solver_get_result (solver, 0);

  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, 0, 3, direct);

  for (i = 0; i < 4; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (solved, i), ==, ncm_vector_get (direct, i), 1.0e-6, 0.0);

  ncm_vector_free (direct);
  nc_xcor_solver_free (solver);
  nc_xcor_free (xc);
}

/*
 * Limber, ell 10 to 12: the cross is exactly zero, while the auto of the same kernel is
 * positive, so this is the tier's short-circuit and not an all-zero configuration.
 */
static void
test_nc_xcor_disjoint_limber_zero (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  NcXcor *xc       = _test_nc_xcor_disjoint_xcor (test, NC_XCOR_METHOD_KERNEL_EXACT);
  NcmVector *cross = ncm_vector_new (3);
  NcmVector *auto1 = ncm_vector_new (3);
  guint i;

  _test_nc_xcor_disjoint_kernels (test, 0);

  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, 10, 12, cross);
  nc_xcor_compute (xc, test->k1, test->k1, test->env.cosmo, 10, 12, auto1);

  for (i = 0; i < 3; i++)
  {
    g_assert_cmpfloat (ncm_vector_get (cross, i), ==, 0.0);
    g_assert_cmpfloat (ncm_vector_get (auto1, i), >, 0.0);
  }

  ncm_vector_free (cross);
  ncm_vector_free (auto1);
  nc_xcor_free (xc);
}

/*
 * l_limber = 6 and a request over ell 0 to 9: the Limber tail is exactly zero, and the
 * non-Limber head is nonzero and equal to a request that stops short of the threshold,
 * so the split does not perturb the multipoles below it. A request entirely at or
 * above the threshold is zero throughout.
 */
#define TEST_L_LIMBER 6
#define TEST_LMAX 9

static void
test_nc_xcor_disjoint_split_at_l_limber (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  NcXcor *xc            = _test_nc_xcor_disjoint_xcor (test, NC_XCOR_METHOD_KERNEL_EXACT);
  NcmVector *straddling = ncm_vector_new (TEST_LMAX + 1);
  NcmVector *head       = ncm_vector_new (TEST_L_LIMBER);
  NcmVector *tail       = ncm_vector_new (TEST_LMAX - TEST_L_LIMBER + 1);
  gboolean any_nonzero  = FALSE;
  guint i;

  _test_nc_xcor_disjoint_kernels (test, TEST_L_LIMBER);

  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, 0, TEST_LMAX, straddling);
  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, 0, TEST_L_LIMBER - 1, head);
  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, TEST_L_LIMBER, TEST_LMAX, tail);

  for (i = 0; i < TEST_L_LIMBER; i++)
  {
    const gdouble c = ncm_vector_get (straddling, i);

    g_assert_true (gsl_finite (c));
    any_nonzero = any_nonzero || (c != 0.0);
    ncm_assert_cmpdouble_e (c, ==, ncm_vector_get (head, i), 1.0e-9, 0.0);
  }

  g_assert_true (any_nonzero);

  for (i = TEST_L_LIMBER; i <= TEST_LMAX; i++)
  {
    g_assert_cmpfloat (ncm_vector_get (straddling, i), ==, 0.0);
    g_assert_cmpfloat (ncm_vector_get (tail, i - TEST_L_LIMBER), ==, 0.0);
  }

  ncm_vector_free (straddling);
  ncm_vector_free (head);
  ncm_vector_free (tail);
  nc_xcor_free (xc);
}

/* The redshift-space Limber method keeps its overlap short-circuit: zero at ell 10 to 12. */
static void
test_nc_xcor_disjoint_limber_z (TestNcXcorDisjoint *test, gconstpointer pdata)
{
  NcXcor *xc       = _test_nc_xcor_disjoint_xcor (test, NC_XCOR_METHOD_LIMBER_Z_GSL);
  NcmVector *cross = ncm_vector_new (3);
  guint i;

  _test_nc_xcor_disjoint_kernels (test, -1);

  nc_xcor_compute (xc, test->k1, test->k2, test->env.cosmo, 10, 12, cross);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (ncm_vector_get (cross, i), ==, 0.0);

  ncm_vector_free (cross);
  nc_xcor_free (xc);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorDisjoint *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"cross_nonzero",     &test_nc_xcor_disjoint_cross_nonzero    },
  {"solver",            &test_nc_xcor_disjoint_solver           },
  {"limber_zero",       &test_nc_xcor_disjoint_limber_zero      },
  {"split_at_l_limber", &test_nc_xcor_disjoint_split_at_l_limber},
  {"limber_z",          &test_nc_xcor_disjoint_limber_z         },
};

gint
main (gint argc, gchar *argv[])
{
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (test_checks); i++)
  {
    gchar *path = g_strdup_printf ("/nc/xcor/disjoint/%s", test_checks[i].name);

    g_test_add (path, TestNcXcorDisjoint, NULL,
                &test_nc_xcor_disjoint_new, test_checks[i].func, &test_nc_xcor_disjoint_free);

    g_free (path);
  }

  g_test_run ();
}


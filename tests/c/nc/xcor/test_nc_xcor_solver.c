/***************************************************************************
 *            test_nc_xcor_solver.c
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
 * #NcXcorSolver batches many spectra into one pass: kernels are registered, pairs are
 * requested, and the union of the requested multipole ranges is tiled into blocks that
 * are solved together. Registration and the block planner are pure bookkeeping and are
 * checked exactly here; solve() is run over a short band with Limber kernels, which is
 * cheap, and its agreement with nc_xcor_compute() on the same pair is what says the
 * batching does not change the answer.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_ZMAX 6.0

typedef struct _TestNcXcorSolver
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcXcorKernel *k1;
  NcXcorKernel *k2;
  NcXcorSolver *solver;
} TestNcXcorSolver;

static void
test_nc_xcor_solver_new (TestNcXcorSolver *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,      70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, 0.255);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B, 0.045);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W, -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", 0.0, NULL);

  test->cosmo = cosmo;
  test->dist  = nc_distance_new (TEST_ZMAX);
  test->ps    = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                       NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));

  nc_distance_prepare (test->dist, cosmo);
  ncm_powspec_prepare (test->ps, NCM_MODEL (cosmo));

  test->k1 = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (test->dist, test->ps, 1500.0, 300.0, 4.0));
  test->k2 = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (test->dist, test->ps, 1700.0, 350.0, 4.0));

  nc_xcor_kernel_prepare (test->k1, cosmo);
  nc_xcor_kernel_prepare (test->k2, cosmo);

  test->solver = nc_xcor_solver_new ();
  g_assert_true (NC_IS_XCOR_SOLVER (test->solver));
}

static void
test_nc_xcor_solver_free (TestNcXcorSolver *test, gconstpointer pdata)
{
  nc_xcor_kernel_free (test->k1);
  nc_xcor_kernel_free (test->k2);
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_solver_free, test->solver);
}

static void
test_nc_xcor_solver_register (TestNcXcorSolver *test, gconstpointer pdata)
{
  NcXcorSolver *solver2 = nc_xcor_solver_ref (test->solver);
  guint id1, id2;

  g_assert_cmpuint (nc_xcor_solver_get_n_kernels (test->solver), ==, 0);
  g_assert_cmpuint (nc_xcor_solver_get_n_requests (test->solver), ==, 0);

  id1 = nc_xcor_solver_register_kernel (test->solver, test->k1);
  id2 = nc_xcor_solver_register_kernel (test->solver, test->k2);

  /* The id is the handle every request is written in terms of, so it has to address the
   * kernel it was returned for. */
  g_assert_cmpuint (id1, !=, id2);
  g_assert_cmpuint (nc_xcor_solver_get_n_kernels (test->solver), ==, 2);
  g_assert_true (nc_xcor_solver_peek_kernel (test->solver, id1) == test->k1);
  g_assert_true (nc_xcor_solver_peek_kernel (test->solver, id2) == test->k2);

  nc_xcor_solver_clear (&solver2);
  g_assert_true (solver2 == NULL);
}

static void
test_nc_xcor_solver_requests (TestNcXcorSolver *test, gconstpointer pdata)
{
  const guint id1 = nc_xcor_solver_register_kernel (test->solver, test->k1);
  const guint id2 = nc_xcor_solver_register_kernel (test->solver, test->k2);
  guint a, b, lmin, lmax;

  nc_xcor_solver_request_cl (test->solver, id1, id2, 2, 12);
  nc_xcor_solver_request_cl (test->solver, id1, id1, 4, 20);

  g_assert_cmpuint (nc_xcor_solver_get_n_requests (test->solver), ==, 2);

  nc_xcor_solver_get_request (test->solver, 0, &a, &b, &lmin, &lmax);
  g_assert_cmpuint (a, ==, id1);
  g_assert_cmpuint (b, ==, id2);
  g_assert_cmpuint (lmin, ==, 2);
  g_assert_cmpuint (lmax, ==, 12);

  nc_xcor_solver_get_request (test->solver, 1, &a, &b, &lmin, &lmax);
  g_assert_cmpuint (a, ==, id1);
  g_assert_cmpuint (b, ==, id1);
  g_assert_cmpuint (lmin, ==, 4);
  g_assert_cmpuint (lmax, ==, 20);

  /* Clearing drops the requests and keeps the kernels: the registration is the expensive
   * half and a caller reusing a solver for a new set of spectra depends on that. */
  nc_xcor_solver_clear_requests (test->solver);
  g_assert_cmpuint (nc_xcor_solver_get_n_requests (test->solver), ==, 0);
  g_assert_cmpuint (nc_xcor_solver_get_n_kernels (test->solver), ==, 2);
}

static void
test_nc_xcor_solver_plan (TestNcXcorSolver *test, gconstpointer pdata)
{
  const guint id1        = nc_xcor_solver_register_kernel (test->solver, test->k1);
  const guint id2        = nc_xcor_solver_register_kernel (test->solver, test->k2);
  const guint block_size = 4;
  guint n_blocks, i, prev_lmax = 0;

  nc_xcor_solver_request_cl (test->solver, id1, id2, 2, 12);
  nc_xcor_solver_request_cl (test->solver, id1, id1, 8, 20);

  nc_xcor_solver_plan_blocks (test->solver, block_size);
  n_blocks = nc_xcor_solver_get_n_blocks (test->solver);

  g_assert_cmpuint (n_blocks, >, 0);

  for (i = 0; i < n_blocks; i++)
  {
    guint lmin, lmax;

    nc_xcor_solver_get_block (test->solver, i, &lmin, &lmax);

    g_assert_cmpuint (lmin, <=, lmax);
    g_assert_cmpuint (lmax - lmin + 1, <=, block_size);

    /* The blocks tile the union of the requested ranges contiguously. A gap would drop
     * multipoles silently; an overlap would compute them twice and disagree. */
    if (i == 0)
      g_assert_cmpuint (lmin, ==, 2);
    else
      g_assert_cmpuint (lmin, ==, prev_lmax + 1);

    prev_lmax = lmax;
  }

  g_assert_cmpuint (prev_lmax, ==, 20);

  /* Planning again replaces the previous blocks rather than appending to them. */
  nc_xcor_solver_plan_blocks (test->solver, block_size);
  g_assert_cmpuint (nc_xcor_solver_get_n_blocks (test->solver), ==, n_blocks);
}

/* A kernel evaluated in Limber mode below l_limber and non-Limber above it needs the
 * switch to fall on a block boundary, or one block would have to be both. */
static void
test_nc_xcor_solver_plan_l_limber (TestNcXcorSolver *test, gconstpointer pdata)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 32));
  NcXcorKernel *split       = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new_full (test->dist, test->ps,
                                                                                      1500.0, 300.0, 4.0, sbi));
  guint id, n_blocks, i;
  gboolean boundary_at_9 = FALSE;

  nc_xcor_kernel_set_l_limber (split, 9);
  nc_xcor_kernel_prepare (split, test->cosmo);

  id = nc_xcor_solver_register_kernel (test->solver, split);
  nc_xcor_solver_request_cl (test->solver, id, id, 2, 20);

  /* A block size that would happily span the threshold if nothing forced a break. */
  nc_xcor_solver_plan_blocks (test->solver, 16);
  n_blocks = nc_xcor_solver_get_n_blocks (test->solver);

  for (i = 0; i < n_blocks; i++)
  {
    guint lmin, lmax;

    nc_xcor_solver_get_block (test->solver, i, &lmin, &lmax);

    if (lmin == 9)
      boundary_at_9 = TRUE;

    /* No block may contain the threshold strictly inside it. */
    g_assert_false ((lmin < 9) && (lmax >= 9));
  }

  g_assert_true (boundary_at_9);

  nc_xcor_kernel_free (split);
  ncm_sbessel_integrator_free (sbi);
}

/* Batching is a scheduling change, not a numerical one: a spectrum solved in blocks must
 * be the spectrum computed directly. */
static void
test_nc_xcor_solver_solve (TestNcXcorSolver *test, gconstpointer pdata)
{
  const guint lmin  = 2;
  const guint lmax  = 11;
  const guint n_l   = lmax - lmin + 1;
  NcXcor *xc        = nc_xcor_new (test->dist, test->ps, NC_XCOR_METHOD_LIMBER_Z_GSL);
  NcmVector *direct = ncm_vector_new (n_l);
  const guint id1   = nc_xcor_solver_register_kernel (test->solver, test->k1);
  const guint id2   = nc_xcor_solver_register_kernel (test->solver, test->k2);
  NcmVector *batched;
  guint i;

  nc_xcor_prepare (xc, test->cosmo);

  nc_xcor_solver_request_cl (test->solver, id1, id2, lmin, lmax);
  nc_xcor_solver_plan_blocks (test->solver, 4);
  nc_xcor_solver_solve (test->solver, xc, test->cosmo);

  batched = nc_xcor_solver_get_result (test->solver, 0);

  g_assert_nonnull (batched);
  g_assert_cmpuint (ncm_vector_len (batched), ==, n_l);

  nc_xcor_compute (xc, test->k1, test->k2, test->cosmo, lmin, lmax, direct);

  for (i = 0; i < n_l; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (batched, i)));
    ncm_assert_cmpdouble_e (ncm_vector_get (batched, i), ==, ncm_vector_get (direct, i), 1.0e-9, 0.0);
  }

  ncm_vector_free (direct);
  nc_xcor_free (xc);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorSolver *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"register",       &test_nc_xcor_solver_register      },
  {"requests",       &test_nc_xcor_solver_requests      },
  {"plan",           &test_nc_xcor_solver_plan          },
  {"plan/l_limber",  &test_nc_xcor_solver_plan_l_limber },
  {"solve",          &test_nc_xcor_solver_solve         },
};

gint
main (gint argc, gchar *argv[])
{
  guint j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
  {
    gchar *path = g_strdup_printf ("/nc/xcor/solver/%s", test_checks[j].name);

    g_test_add (path, TestNcXcorSolver, NULL,
                &test_nc_xcor_solver_new, test_checks[j].func, &test_nc_xcor_solver_free);

    g_free (path);
  }

  g_test_run ();
}


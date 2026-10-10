/***************************************************************************
 *            test_nc_xcor_kquad_error.c
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
 * The error estimate of NC_XCOR_METHOD_KERNEL_EXACT on an auto spectrum, where the
 * positive integrand leaves nothing to amplify: it must bound the true error, taken
 * against the same kernel at reltol 1e-8, peak-epsilon at its floor of 1e-6 and
 * cheb-reltol 1e-10, and stay far below one, for both closures.
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

#define TEST_LMIN 2
#define TEST_LMAX 5
#define TEST_NELL (TEST_LMAX - TEST_LMIN + 1)

/**
 * TestCase:
 * @name: path component naming the case
 * @closure: how the kernel represents W_l(k)
 */
typedef struct _TestCase
{
  const gchar *name;
  NcXcorKernelClosure closure;
} TestCase;

static const TestCase test_cases[] = {
  {"spline",    NC_XCOR_KERNEL_CLOSURE_SPLINE   },
  {"chebyshev", NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV},
};

typedef struct _TestNcXcorKQuadError
{
  TestNcXcorKQuadEnv env;
  NcmSBesselIntegrator *sbi_tight;
  NcXcorKernel *kernel;
  NcXcorKernel *tight;
  NcXcor *xc;
} TestNcXcorKQuadError;

static void
test_nc_xcor_kquad_error_new (TestNcXcorKQuadError *test, gconstpointer pdata)
{
  const TestCase *tc = pdata;

  test_nc_xcor_kquad_env_init (&test->env);

  test->sbi_tight = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));
  ncm_sbessel_integrator_levin_set_cheb_reltol (NCM_SBESSEL_INTEGRATOR_LEVIN (test->sbi_tight), 1.0e-10);

  test->kernel = test_nc_xcor_kquad_tophat_converged (&test->env, test->env.sbi, 200.0, 400.0, -1, 0.0, 0.0);
  test->tight  = test_nc_xcor_kquad_tophat_converged (&test->env, test->sbi_tight, 200.0, 400.0, -1, 1.0e-8, 1.0e-6);

  test->xc = nc_xcor_new (test->env.dist, test->env.ps, NC_XCOR_METHOD_KERNEL_EXACT);
  nc_xcor_set_closure_type (test->xc, tc->closure);
  nc_xcor_prepare (test->xc, test->env.cosmo);
}

static void
test_nc_xcor_kquad_error_free (TestNcXcorKQuadError *test, gconstpointer pdata)
{
  nc_xcor_kernel_free (test->kernel);
  nc_xcor_kernel_free (test->tight);
  ncm_sbessel_integrator_free (test->sbi_tight);
  test_nc_xcor_kquad_env_clear (&test->env);

  NCM_TEST_FREE (nc_xcor_free, test->xc);
}

static void
test_nc_xcor_kquad_error_bounds_truth (TestNcXcorKQuadError *test, gconstpointer pdata)
{
  NcmVector *cl     = ncm_vector_new (TEST_NELL);
  NcmVector *cl_err = ncm_vector_new (TEST_NELL);
  NcmVector *truth  = ncm_vector_new (TEST_NELL);
  guint i;

  nc_xcor_compute_full (test->xc, test->kernel, NULL, test->env.cosmo, TEST_LMIN, TEST_LMAX, cl, cl_err);
  nc_xcor_compute (test->xc, test->tight, NULL, test->env.cosmo, TEST_LMIN, TEST_LMAX, truth);

  for (i = 0; i < TEST_NELL; i++)
  {
    const gdouble C        = ncm_vector_get (cl, i);
    const gdouble rel_est  = fabs (ncm_vector_get (cl_err, i) / C);
    const gdouble rel_true = fabs (C / ncm_vector_get (truth, i) - 1.0);

    g_test_message ("ell %u: rel_est %.3e, rel_true %.3e", TEST_LMIN + i, rel_est, rel_true);

    g_assert_cmpfloat (C, >, 0.0);
    g_assert_cmpfloat (rel_est, >=, rel_true);
    g_assert_cmpfloat (rel_est, <, 1.0e-2);
  }

  ncm_vector_free (cl);
  ncm_vector_free (cl_err);
  ncm_vector_free (truth);
}

gint
main (gint argc, gchar *argv[])
{
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (test_cases); i++)
  {
    gchar *path = g_strdup_printf ("/nc/xcor/kquad/error/%s/bounds_truth", test_cases[i].name);

    g_test_add (path, TestNcXcorKQuadError, &test_cases[i],
                &test_nc_xcor_kquad_error_new, &test_nc_xcor_kquad_error_bounds_truth, &test_nc_xcor_kquad_error_free);

    g_free (path);
  }

  g_test_run ();
}


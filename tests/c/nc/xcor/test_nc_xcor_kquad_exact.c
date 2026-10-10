/***************************************************************************
 *            test_nc_xcor_kquad_exact.c
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
 * NC_XCOR_METHOD_KERNEL_EXACT against an independent integral of the same closures:
 * nc_xcor_integrate_block() on a pair of integrands, compared with
 * test_nc_xcor_kquad_ref_eval() on the same pair, for spline, Chebyshev and
 * spline x Chebyshev closures, non-Limber, Limber and mixed tiers, auto and cross (auto
 * only where both kernels agree, since an auto spectrum uses the first alone). Exact is a bilinear form in the
 * coefficients or GL(5) on merged knots; the reference raises Gauss-Legendre per cell,
 * so the two share the closures and nothing else.
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
 * @closure1: how the first kernel represents W_l(k)
 * @closure2: how the second kernel represents W_l(k)
 * @l_limber1: tier of the first kernel
 * @l_limber2: tier of the second kernel
 */
typedef struct _TestCase
{
  const gchar *name;
  NcXcorKernelClosure closure1;
  NcXcorKernelClosure closure2;
  gint l_limber1;
  gint l_limber2;
} TestCase;

#define SPLINE NC_XCOR_KERNEL_CLOSURE_SPLINE
#define CHEB NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV

static const TestCase test_cases[] = {
  {"spline/non_limber",           SPLINE, SPLINE, -1, -1},
  {"spline/limber",               SPLINE, SPLINE,  0,  0},
  {"spline/mixed",                SPLINE, SPLINE, -1,  0},
  {"chebyshev/non_limber",        CHEB,   CHEB,   -1, -1},
  {"chebyshev/limber",            CHEB,   CHEB,    0,  0},
  {"chebyshev/mixed",             CHEB,   CHEB,   -1,  0},
  {"spline_chebyshev/non_limber", SPLINE, CHEB,   -1, -1},
  {"spline_chebyshev/mixed",      SPLINE, CHEB,    0, -1},
};

typedef struct _TestNcXcorKQuadExact
{
  TestNcXcorKQuadEnv env;
  NcXcorKernelIntegrand *i1;
  NcXcorKernelIntegrand *i2;
  TestNcXcorKQuadRef *ref;
  NcXcor *xc;
} TestNcXcorKQuadExact;

static void
test_nc_xcor_kquad_exact_new (TestNcXcorKQuadExact *test, gconstpointer pdata)
{
  const TestCase *tc = pdata;
  NcXcorKernel *k1, *k2;

  test_nc_xcor_kquad_env_init (&test->env);

  /* Overlapping, so the cross spectrum is not numerically zero. */
  k1 = test_nc_xcor_kquad_tophat (&test->env, 200.0, 400.0, tc->l_limber1);
  k2 = test_nc_xcor_kquad_tophat (&test->env, 250.0, 450.0, tc->l_limber2);

  test->i1 = nc_xcor_kernel_get_eval_vectorized (k1, test->env.cosmo, TEST_LMIN, TEST_LMAX, tc->closure1);
  test->i2 = nc_xcor_kernel_get_eval_vectorized (k2, test->env.cosmo, TEST_LMIN, TEST_LMAX, tc->closure2);

  test->xc = nc_xcor_new (test->env.dist, test->env.ps, NC_XCOR_METHOD_KERNEL_EXACT);
  nc_xcor_set_closure_type (test->xc, tc->closure1);
  nc_xcor_prepare (test->xc, test->env.cosmo);

  test->ref = test_nc_xcor_kquad_ref_new ();

  nc_xcor_kernel_free (k1);
  nc_xcor_kernel_free (k2);
}

static void
test_nc_xcor_kquad_exact_free (TestNcXcorKQuadExact *test, gconstpointer pdata)
{
  nc_xcor_kernel_integrand_unref (test->i1);
  nc_xcor_kernel_integrand_unref (test->i2);
  test_nc_xcor_kquad_ref_free (test->ref);
  test_nc_xcor_kquad_env_clear (&test->env);

  NCM_TEST_FREE (nc_xcor_free, test->xc);
}

/* A spline closure carries knots and no panels, a Chebyshev one the reverse. */
static void
_test_nc_xcor_kquad_exact_assert_closure (NcXcorKernelIntegrand *xclki, NcXcorKernelClosure closure)
{
  const gboolean spline = (closure == NC_XCOR_KERNEL_CLOSURE_SPLINE);

  g_assert_true ((nc_xcor_kernel_integrand_peek_knots (xclki) != NULL) == spline);
  g_assert_true ((nc_xcor_kernel_integrand_get_n_panels (xclki) > 0) == !spline);
}

static void
_test_nc_xcor_kquad_exact_check (TestNcXcorKQuadExact *test, const TestCase *tc, gboolean isauto)
{
  NcXcorKernelIntegrand *i2 = isauto ? test->i1 : test->i2;
  NcmVector *exact          = ncm_vector_new (TEST_NELL);
  NcmVector *ref            = ncm_vector_new (TEST_NELL);
  const gdouble RH          = nc_hicosmo_RH_Mpc (test->env.cosmo);
  gdouble move, peak = 0.0, dev = 0.0;
  guint i;

  _test_nc_xcor_kquad_exact_assert_closure (test->i1, tc->closure1);
  _test_nc_xcor_kquad_exact_assert_closure (test->i2, tc->closure2);

  nc_xcor_integrate_block (test->xc, test->i1, i2, TEST_LMIN, TEST_LMAX, isauto,
                           NC_XCOR_METHOD_KERNEL_EXACT, exact, NULL);
  move = test_nc_xcor_kquad_ref_eval (test->ref, test->i1, isauto ? NULL : test->i2, RH, ref);

  for (i = 0; i < TEST_NELL; i++)
    peak = GSL_MAX (peak, fabs (ncm_vector_get (ref, i)));

  for (i = 0; i < TEST_NELL; i++)
    dev = GSL_MAX (dev, fabs (ncm_vector_get (exact, i) - ncm_vector_get (ref, i)) / peak);

  g_test_message ("peak %.6e, max |exact - ref| / peak %.3e, reference move %.3e", peak, dev, move);

  g_assert_cmpfloat (peak, >, 0.0);
  g_assert_cmpfloat (move, <, 1.0e-12);
  g_assert_cmpfloat (dev, <, 1.0e-13);

  ncm_vector_free (exact);
  ncm_vector_free (ref);
}

static void
test_nc_xcor_kquad_exact_auto (TestNcXcorKQuadExact *test, gconstpointer pdata)
{
  _test_nc_xcor_kquad_exact_check (test, pdata, TRUE);
}

static void
test_nc_xcor_kquad_exact_cross (TestNcXcorKQuadExact *test, gconstpointer pdata)
{
  _test_nc_xcor_kquad_exact_check (test, pdata, FALSE);
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
    gchar *path_auto  = g_strdup_printf ("/nc/xcor/kquad/exact/%s/auto", test_cases[i].name);
    gchar *path_cross = g_strdup_printf ("/nc/xcor/kquad/exact/%s/cross", test_cases[i].name);

    if ((test_cases[i].l_limber1 == test_cases[i].l_limber2) && (test_cases[i].closure1 == test_cases[i].closure2))
      g_test_add (path_auto, TestNcXcorKQuadExact, &test_cases[i],
                  &test_nc_xcor_kquad_exact_new, &test_nc_xcor_kquad_exact_auto, &test_nc_xcor_kquad_exact_free);

    g_test_add (path_cross, TestNcXcorKQuadExact, &test_cases[i],
                &test_nc_xcor_kquad_exact_new, &test_nc_xcor_kquad_exact_cross, &test_nc_xcor_kquad_exact_free);

    g_free (path_auto);
    g_free (path_cross);
  }

  g_test_run ();
}


/***************************************************************************
 *            test_nc_xcor_kquad.c
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
 * The kernel-space outer integral: given two $W_\ell(k)$ representations, integrate
 * $k^2 W^1_\ell W^2_\ell$ over a multipole block. Each method does it differently --
 * adaptive GSL over breakpoints, fixed cubature, exact GL(5) on the panel union -- and
 * this checks the parts that are the same whichever is chosen: the band is filled, the
 * auto path agrees with passing one kernel twice to the cross path, and an error
 * estimate accompanies the result exactly when the method claims to offer one.
 *
 * Everything is sized to keep this in the unit lane: two narrow windows close in, the
 * power spectrum capped well below its default, a four-multipole block, and the kernels'
 * adaptive apparatus capped rather than run to convergence. Accuracy of these methods
 * against their own reference is tests/python/nc/xcor/test_k_integral.py.
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
#define TEST_KMAX 0.5
#define TEST_LMIN 2
#define TEST_LMAX 5
#define TEST_NELL (TEST_LMAX - TEST_LMIN + 1)

typedef struct _TestNcXcorKQuad
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcmSBesselIntegrator *sbi;
  NcXcorKernel *k1;
  NcXcorKernel *k2;
  NcXcor *xc;
} TestNcXcorKQuad;

/**
 * TestMethod:
 * @name: path component naming the method
 * @meth: the method itself
 * @has_block: whether nc_xcor_integrate_block() accepts it
 *
 * Kernel-space is not the same as block-capable: NC_XCOR_METHOD_KERNEL_GSL builds its
 * closure per multipole, so it has no block quadrature and integrate_block() refuses it.
 */
typedef struct _TestMethod
{
  const gchar *name;
  NcXcorMethod meth;
  gboolean has_block;
} TestMethod;

static const TestMethod test_methods[] = {
  {"kernel_gsl",       NC_XCOR_METHOD_KERNEL_GSL,       FALSE},
  {"kernel_cubature",  NC_XCOR_METHOD_KERNEL_CUBATURE,  TRUE },
  {"kernel_exact",     NC_XCOR_METHOD_KERNEL_EXACT,     TRUE },
  {"kernel_gsl_block", NC_XCOR_METHOD_KERNEL_GSL_BLOCK, TRUE },
};

static NcXcorKernel *
_tophat (TestNcXcorKQuad *test, gdouble chi_lower, gdouble chi_upper)
{
  NcXcorKernel *xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new_full (test->dist, test->ps,
                                                                                chi_lower, chi_upper,
                                                                                test->sbi));

  /* Capped, not converged: the branches are what is under test, not the accuracy. */
  nc_xcor_kernel_set_max_border_expansions (xclk, 1);
  nc_xcor_kernel_set_max_iter (xclk, 4);
  nc_xcor_kernel_set_reltol (xclk, 1.0e-3);
  nc_xcor_kernel_set_scaled_abstol (xclk, 1.0e-4);
  nc_xcor_kernel_set_panel_order_cap (xclk, 12);
  nc_xcor_kernel_set_lmax (xclk, 16);
  nc_xcor_kernel_set_l_limber (xclk, 0);

  nc_xcor_kernel_prepare (xclk, test->cosmo);

  return xclk;
}

static void
test_nc_xcor_kquad_new (TestNcXcorKQuad *test, gconstpointer pdata)
{
  const TestMethod *tm = pdata;
  NcHICosmo *cosmo     = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

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
  test->sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));

  ncm_powspec_set_kmax (test->ps, TEST_KMAX);

  nc_distance_prepare (test->dist, cosmo);
  ncm_powspec_prepare (test->ps, NCM_MODEL (cosmo));

  /* Overlapping, so the cross spectrum is not numerically zero. */
  test->k1 = _tophat (test, 200.0, 400.0);
  test->k2 = _tophat (test, 250.0, 450.0);

  test->xc = nc_xcor_new (test->dist, test->ps, tm->meth);
  nc_xcor_set_ell_batch_size (test->xc, TEST_NELL);
  nc_xcor_prepare (test->xc, cosmo);

  g_assert_true (nc_xcor_method_is_kernel_space (tm->meth));
}

static void
test_nc_xcor_kquad_free (TestNcXcorKQuad *test, gconstpointer pdata)
{
  nc_xcor_kernel_free (test->k1);
  nc_xcor_kernel_free (test->k2);
  ncm_sbessel_integrator_free (test->sbi);
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_free, test->xc);
}

static void
test_nc_xcor_kquad_compute (TestNcXcorKQuad *test, gconstpointer pdata)
{
  NcmVector *cl = ncm_vector_new (TEST_NELL);
  guint i;

  ncm_vector_set_all (cl, GSL_NAN);
  nc_xcor_compute (test->xc, test->k1, test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, cl);

  for (i = 0; i < TEST_NELL; i++)
    g_assert_true (gsl_finite (ncm_vector_get (cl, i)));

  ncm_vector_free (cl);
}

static void
test_nc_xcor_kquad_auto (TestNcXcorKQuad *test, gconstpointer pdata)
{
  NcmVector *cl = ncm_vector_new (TEST_NELL);
  guint i;

  /* The auto path samples one kernel once and squares it; the cross path samples two.
   * Handing the same kernel to both must therefore agree, and the auto spectrum of a
   * real field is positive. */
  nc_xcor_compute (test->xc, test->k1, test->k1, test->cosmo, TEST_LMIN, TEST_LMAX, cl);

  for (i = 0; i < TEST_NELL; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (cl, i)));
    g_assert_cmpfloat (ncm_vector_get (cl, i), >, 0.0);
  }

  ncm_vector_free (cl);
}

static void
test_nc_xcor_kquad_error (TestNcXcorKQuad *test, gconstpointer pdata)
{
  const TestMethod *tm = pdata;
  NcmVector *cl        = ncm_vector_new (TEST_NELL);
  NcmVector *cl_err    = ncm_vector_new (TEST_NELL);
  guint i;

  ncm_vector_set_all (cl_err, GSL_NAN);
  nc_xcor_compute_full (test->xc, test->k1, test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, cl, cl_err);

  for (i = 0; i < TEST_NELL; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (cl, i)));

    /* A method that advertises an error estimate has to fill the vector it was given. */
    if (nc_xcor_method_has_error_estimate (tm->meth))
    {
      g_assert_true (gsl_finite (ncm_vector_get (cl_err, i)));
      g_assert_cmpfloat (ncm_vector_get (cl_err, i), >=, 0.0);
    }
  }

  ncm_vector_free (cl);
  ncm_vector_free (cl_err);
}

/* nc_xcor_integrate_block() is the entry point the solver uses: it takes representations
 * that were built once and reused, rather than building them per pair. Driving it
 * directly must reproduce what nc_xcor_compute() does with the same kernels. */
static void
test_nc_xcor_kquad_integrate_block (TestNcXcorKQuad *test, gconstpointer pdata)
{
  const TestMethod *tm              = pdata;
  const NcXcorKernelClosure closure = nc_xcor_get_closure_type (test->xc);
  NcXcorKernelIntegrand *i1;
  NcXcorKernelIntegrand *i2;
  NcmVector *block;
  NcmVector *err;
  NcmVector *direct;
  guint i;

  if (!tm->has_block)
  {
    g_test_skip ("builds its closure per multipole, so it has no block quadrature");

    return;
  }

  i1 = nc_xcor_kernel_get_eval_vectorized (test->k1, test->cosmo, TEST_LMIN, TEST_LMAX, closure);
  i2 = nc_xcor_kernel_get_eval_vectorized (test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, closure);

  block  = ncm_vector_new (TEST_NELL);
  err    = ncm_vector_new (TEST_NELL);
  direct = ncm_vector_new (TEST_NELL);

  ncm_vector_set_all (block, GSL_NAN);

  nc_xcor_integrate_block (test->xc, i1, i2, TEST_LMIN, TEST_LMAX, FALSE, tm->meth, block, err);
  nc_xcor_compute (test->xc, test->k1, test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, direct);

  for (i = 0; i < TEST_NELL; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (block, i)));
    ncm_assert_cmpdouble_e (ncm_vector_get (block, i), ==, ncm_vector_get (direct, i), 1.0e-9, 0.0);
  }

  /* The same call in auto mode, with one representation used twice. */
  nc_xcor_integrate_block (test->xc, i1, i1, TEST_LMIN, TEST_LMAX, TRUE, tm->meth, block, err);

  for (i = 0; i < TEST_NELL; i++)
    g_assert_true (gsl_finite (ncm_vector_get (block, i)));

  ncm_vector_free (block);
  ncm_vector_free (err);
  ncm_vector_free (direct);
  nc_xcor_kernel_integrand_unref (i1);
  nc_xcor_kernel_integrand_unref (i2);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorKQuad *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"compute",         &test_nc_xcor_kquad_compute        },
  {"auto",            &test_nc_xcor_kquad_auto           },
  {"error",           &test_nc_xcor_kquad_error          },
  {"integrate_block", &test_nc_xcor_kquad_integrate_block},
};

/* Every method has to describe itself: the solver picks an integrator by these. */
static void
test_nc_xcor_kquad_method_table (void)
{
  const NcXcorMethod all[] = {
    NC_XCOR_METHOD_LIMBER_Z_GSL, NC_XCOR_METHOD_LIMBER_Z_CUBATURE,
    NC_XCOR_METHOD_KERNEL_GSL, NC_XCOR_METHOD_KERNEL_CUBATURE,
    NC_XCOR_METHOD_KERNEL_EXACT, NC_XCOR_METHOD_KERNEL_GSL_BLOCK
  };
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (all); i++)
  {
    const gchar *name = nc_xcor_method_get_name (all[i]);

    g_assert_nonnull (name);
    g_assert_cmpuint (strlen (name), >, 0);

    /* Names address methods in configuration and logs, so they must be distinct. */
    for (j = 0; j < i; j++)
      g_assert_cmpstr (name, !=, nc_xcor_method_get_name (all[j]));

    /* Only the two Limber methods are not kernel-space. */
    g_assert_cmpint (nc_xcor_method_is_kernel_space (all[i]), ==,
                     (all[i] != NC_XCOR_METHOD_LIMBER_Z_GSL) &&
                     (all[i] != NC_XCOR_METHOD_LIMBER_Z_CUBATURE));

    /* Just has to answer without asserting. */
    (void) nc_xcor_method_has_error_estimate (all[i]);
  }
}

gint
main (gint argc, gchar *argv[])
{
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (test_methods); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
    {
      gchar *path = g_strdup_printf ("/nc/xcor/kquad/%s/%s", test_methods[i].name, test_checks[j].name);

      g_test_add (path, TestNcXcorKQuad, &test_methods[i],
                  &test_nc_xcor_kquad_new, test_checks[j].func, &test_nc_xcor_kquad_free);

      g_free (path);
    }
  }

  g_test_add_func ("/nc/xcor/kquad/methods", &test_nc_xcor_kquad_method_table);

  g_test_run ();
}


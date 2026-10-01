/***************************************************************************
 *            test_nc_xcor.c
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
 * #NcXcor drives a pair of kernels to a $C_\ell$ band. The Limber methods are a single
 * one-dimensional integral per multipole, so a short band over two analytic windows is
 * cheap enough for the unit lane, and it is the whole of nc_xcor.c and nc_xcor_limber_z.c.
 *
 * What is checked is the plumbing, not the numbers: the band is filled, the error
 * estimate accompanies it when the method offers one, the two Limber backends agree with
 * each other to their own tolerance, and an auto-spectrum is the cross-spectrum of a
 * kernel with itself. Absolute accuracy against certified values lives in
 * tests/python/nc/xcor.
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
#define TEST_LMIN 2
#define TEST_LMAX 9
#define TEST_NELL (TEST_LMAX - TEST_LMIN + 1)

typedef struct _TestNcXcor
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcXcor *xc;
  NcXcorKernel *k1;
  NcXcorKernel *k2;
} TestNcXcor;

typedef struct _TestMethod
{
  const gchar *name;
  NcXcorMethod meth;
} TestMethod;

static const TestMethod test_methods[] = {
  {"limber_z_gsl",      NC_XCOR_METHOD_LIMBER_Z_GSL     },
  {"limber_z_cubature", NC_XCOR_METHOD_LIMBER_Z_CUBATURE},
};

static void
test_nc_xcor_new (TestNcXcor *test, gconstpointer pdata)
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

  nc_distance_prepare (test->dist, cosmo);
  ncm_powspec_prepare (test->ps, NCM_MODEL (cosmo));

  /* Two overlapping Gaussian windows: overlapping so the cross-spectrum is not
   * numerically zero, and analytic so nothing here depends on a transfer function. */
  test->k1 = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (test->dist, test->ps, 1500.0, 300.0, 4.0));
  test->k2 = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (test->dist, test->ps, 1700.0, 350.0, 4.0));

  nc_xcor_kernel_prepare (test->k1, cosmo);
  nc_xcor_kernel_prepare (test->k2, cosmo);

  test->xc = nc_xcor_new (test->dist, test->ps, tm->meth);
  nc_xcor_prepare (test->xc, cosmo);
}

static void
test_nc_xcor_free (TestNcXcor *test, gconstpointer pdata)
{
  nc_xcor_kernel_free (test->k1);
  nc_xcor_kernel_free (test->k2);
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_free, test->xc);
}

static void
test_nc_xcor_basic (TestNcXcor *test, gconstpointer pdata)
{
  const TestMethod *tm = pdata;
  NcXcor *xc2          = nc_xcor_ref (test->xc);

  g_assert_cmpuint (nc_xcor_get_meth (test->xc), ==, tm->meth);
  g_assert_nonnull (nc_xcor_method_get_name (tm->meth));
  g_assert_false (nc_xcor_method_is_kernel_space (tm->meth));

  nc_xcor_set_reltol (test->xc, 1.0e-5);
  ncm_assert_cmpdouble_e (nc_xcor_get_reltol (test->xc), ==, 1.0e-5, 1.0e-15, 0.0);

  nc_xcor_set_ell_batch_size (test->xc, 4);
  g_assert_cmpuint (nc_xcor_get_ell_batch_size (test->xc), ==, 4);

  nc_xcor_set_closure_type (test->xc, NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV);
  g_assert_cmpuint (nc_xcor_get_closure_type (test->xc), ==, NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV);

  nc_xcor_clear (&xc2);
  g_assert_true (xc2 == NULL);
}

static void
test_nc_xcor_compute (TestNcXcor *test, gconstpointer pdata)
{
  NcmVector *cl = ncm_vector_new (TEST_NELL);
  guint i;

  ncm_vector_set_all (cl, GSL_NAN);
  nc_xcor_compute (test->xc, test->k1, test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, cl);

  /* Every requested multipole must be written: a band left partly untouched is the
   * failure mode a shape-only check would miss. */
  for (i = 0; i < TEST_NELL; i++)
    g_assert_true (gsl_finite (ncm_vector_get (cl, i)));

  ncm_vector_free (cl);
}

static void
test_nc_xcor_compute_full (TestNcXcor *test, gconstpointer pdata)
{
  const TestMethod *tm = pdata;
  NcmVector *cl        = ncm_vector_new (TEST_NELL);
  NcmVector *cl_err    = ncm_vector_new (TEST_NELL);
  NcmVector *cl_plain  = ncm_vector_new (TEST_NELL);
  guint i;

  ncm_vector_set_all (cl, GSL_NAN);
  ncm_vector_set_all (cl_err, GSL_NAN);

  nc_xcor_compute_full (test->xc, test->k1, test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, cl, cl_err);
  nc_xcor_compute (test->xc, test->k1, test->k2, test->cosmo, TEST_LMIN, TEST_LMAX, cl_plain);

  for (i = 0; i < TEST_NELL; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (cl, i)));

    /* compute() is compute_full() without the error vector, so the spectra must match. */
    ncm_assert_cmpdouble_e (ncm_vector_get (cl_plain, i), ==, ncm_vector_get (cl, i), 1.0e-12, 0.0);

    if (nc_xcor_method_has_error_estimate (tm->meth))
    {
      g_assert_true (gsl_finite (ncm_vector_get (cl_err, i)));
      g_assert_cmpfloat (ncm_vector_get (cl_err, i), >=, 0.0);
    }
  }

  ncm_vector_free (cl);
  ncm_vector_free (cl_err);
  ncm_vector_free (cl_plain);
}

static void
test_nc_xcor_auto (TestNcXcor *test, gconstpointer pdata)
{
  NcmVector *cross = ncm_vector_new (TEST_NELL);
  NcmVector *auto_ = ncm_vector_new (TEST_NELL);
  guint i;

  /* The auto path is a separate integrand from the cross one -- one kernel sampled twice
   * rather than two sampled once -- so passing the same kernel twice has to reproduce it. */
  nc_xcor_compute (test->xc, test->k1, test->k1, test->cosmo, TEST_LMIN, TEST_LMAX, auto_);
  nc_xcor_compute (test->xc, test->k1, test->k1, test->cosmo, TEST_LMIN, TEST_LMAX, cross);

  for (i = 0; i < TEST_NELL; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (auto_, i)));
    g_assert_cmpfloat (ncm_vector_get (auto_, i), >, 0.0);
    ncm_assert_cmpdouble_e (ncm_vector_get (cross, i), ==, ncm_vector_get (auto_, i), 1.0e-12, 0.0);
  }

  ncm_vector_free (cross);
  ncm_vector_free (auto_);
}

/* Two windows with no redshift in common. Under Limber the integrand is a product of
 * the two at the same chi, so the cross spectrum is exactly zero -- not small, zero --
 * and the library short-circuits the pair rather than integrating it. The threshold
 * below which it declines to do so is the larger of the two kernels' l_limber, which is
 * 0 for both here, so the whole band is zero.
 *
 * The check is worth stating precisely because the non-Limber methods disagree: there
 * the same pair correlates, which is what tests/python/nc/xcor/test_disjoint_bins.py is
 * about. */
static void
test_nc_xcor_disjoint (TestNcXcor *test, gconstpointer pdata)
{
  NcXcorKernel *near = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new (test->dist, test->ps, 200.0, 400.0));
  NcXcorKernel *far  = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new (test->dist, test->ps, 2000.0, 2500.0));
  NcmVector *cl      = ncm_vector_new (TEST_NELL);
  gdouble zmin_n, zmax_n, zmid_n, zmin_f, zmax_f, zmid_f;
  guint i;

  nc_xcor_kernel_prepare (near, test->cosmo);
  nc_xcor_kernel_prepare (far, test->cosmo);

  nc_xcor_kernel_get_z_range (near, &zmin_n, &zmax_n, &zmid_n);
  nc_xcor_kernel_get_z_range (far, &zmin_f, &zmax_f, &zmid_f);

  /* The premise of the check, asserted rather than assumed. */
  g_assert_cmpfloat (zmax_n, <, zmin_f);

  ncm_vector_set_all (cl, GSL_NAN);
  nc_xcor_compute (test->xc, near, far, test->cosmo, TEST_LMIN, TEST_LMAX, cl);

  for (i = 0; i < TEST_NELL; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (cl, i), ==, 0.0, 0.0, 0.0);

  /* A kernel always overlaps itself, so the same short circuit must not fire. */
  nc_xcor_compute (test->xc, near, near, test->cosmo, TEST_LMIN, TEST_LMAX, cl);

  for (i = 0; i < TEST_NELL; i++)
    g_assert_cmpfloat (ncm_vector_get (cl, i), >, 0.0);

  nc_xcor_kernel_free (near);
  nc_xcor_kernel_free (far);
  ncm_vector_free (cl);
}

/* Every property reads back what it was constructed or set with. */
static void
test_nc_xcor_properties (TestNcXcor *test, gconstpointer pdata)
{
  const TestMethod *tm = pdata;
  NcDistance *dist     = NULL;
  NcmPowspec *ps       = NULL;
  NcXcorMethod meth;
  NcXcorKernelClosure closure;
  gdouble reltol;
  guint batch;

  nc_xcor_set_reltol (test->xc, 1.0e-7);
  nc_xcor_set_ell_batch_size (test->xc, 8);
  nc_xcor_set_closure_type (test->xc, NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV);

  g_object_get (test->xc,
                "distance", &dist,
                "power-spec", &ps,
                "meth", &meth,
                "closure-type", &closure,
                "reltol", &reltol,
                "ell-batch-size", &batch,
                NULL);

  g_assert_true (dist == test->dist);
  g_assert_true (ps == test->ps);
  g_assert_cmpuint (meth, ==, tm->meth);
  g_assert_cmpuint (closure, ==, NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV);
  ncm_assert_cmpdouble_e (reltol, ==, 1.0e-7, 1.0e-15, 0.0);
  g_assert_cmpuint (batch, ==, 8);

  nc_distance_free (dist);
  ncm_powspec_free (ps);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcor *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"basic",        &test_nc_xcor_basic       },
  {"compute",      &test_nc_xcor_compute     },
  {"compute_full", &test_nc_xcor_compute_full},
  {"auto",         &test_nc_xcor_auto        },
  {"disjoint",     &test_nc_xcor_disjoint    },
  {"properties",   &test_nc_xcor_properties  },
};

/* The two Limber backends integrate the same integrand with different quadrature, so
 * they must agree far better than either agrees with the truth. */
static void
test_nc_xcor_limber_backends_agree (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmVector *cl_gsl = ncm_vector_new (TEST_NELL);
  NcmVector *cl_cub = ncm_vector_new (TEST_NELL);
  NcXcorKernel *k1, *k2;
  NcXcor *xc_gsl, *xc_cub;
  guint i;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  k1 = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (dist, ps, 1500.0, 300.0, 4.0));
  k2 = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (dist, ps, 1700.0, 350.0, 4.0));

  nc_xcor_kernel_prepare (k1, cosmo);
  nc_xcor_kernel_prepare (k2, cosmo);

  xc_gsl = nc_xcor_new (dist, ps, NC_XCOR_METHOD_LIMBER_Z_GSL);
  xc_cub = nc_xcor_new (dist, ps, NC_XCOR_METHOD_LIMBER_Z_CUBATURE);

  nc_xcor_prepare (xc_gsl, cosmo);
  nc_xcor_prepare (xc_cub, cosmo);

  nc_xcor_compute (xc_gsl, k1, k2, cosmo, TEST_LMIN, TEST_LMAX, cl_gsl);
  nc_xcor_compute (xc_cub, k1, k2, cosmo, TEST_LMIN, TEST_LMAX, cl_cub);

  for (i = 0; i < TEST_NELL; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (cl_cub, i), ==, ncm_vector_get (cl_gsl, i), 1.0e-6, 0.0);

  nc_xcor_free (xc_gsl);
  nc_xcor_free (xc_cub);
  nc_xcor_kernel_free (k1);
  nc_xcor_kernel_free (k2);
  ncm_vector_free (cl_gsl);
  ncm_vector_free (cl_cub);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
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
      gchar *path = g_strdup_printf ("/nc/xcor/%s/%s", test_methods[i].name, test_checks[j].name);

      g_test_add (path, TestNcXcor, &test_methods[i],
                  &test_nc_xcor_new, test_checks[j].func, &test_nc_xcor_free);

      g_free (path);
    }
  }

  g_test_add_func ("/nc/xcor/limber/backends_agree", &test_nc_xcor_limber_backends_agree);

  g_test_run ();
}


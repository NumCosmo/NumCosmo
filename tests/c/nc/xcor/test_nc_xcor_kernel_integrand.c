/***************************************************************************
 *            test_nc_xcor_kernel_integrand.c
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
 * #NcXcorKernelIntegrand is what a kernel hands the outer k-integral: a callable
 * $W_\ell(k)$ plus the accessors the integrators use to place their panels. This checks
 * that surface -- range, panels, restriction, spectral form, component splitting -- for
 * both closures and for the Limber and non-Limber builders.
 *
 * The windows are deliberately small and the tolerances loose: what is under test is
 * that each accessor is wired to the representation actually built, not that the
 * representation has converged. Accuracy of $W_\ell(k)$ is checked against certified
 * values in tests/python/nc/xcor. Running the adaptive loop to convergence here would
 * cost far more and reach no line the capped path misses.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/* A narrow window close in, and a power spectrum cut well below its default kmax: the
 * Levin ODE's cost is linear in the top of y = k chi_max, and these bound it at a few
 * hundred rather than a few hundred thousand. */
#define TEST_CHI_LOWER 200.0
#define TEST_CHI_UPPER 400.0
#define TEST_KMAX 0.5
#define TEST_ZMAX 6.0
#define TEST_ELL 4

typedef struct _TestNcXcorKernelIntegrand
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcmSBesselIntegrator *sbi;
  NcXcorKernel *xclk;
} TestNcXcorKernelIntegrand;

typedef struct _TestCase
{
  const gchar *name;
  NcXcorKernelClosure closure;
  gboolean non_limber;
  gboolean multi_comp;
} TestCase;

static const TestCase test_cases[] = {
  {"limber/spline",         NC_XCOR_KERNEL_CLOSURE_SPLINE,    FALSE, FALSE},
  {"limber/chebyshev",      NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV, FALSE, FALSE},
  {"non_limber/spline",     NC_XCOR_KERNEL_CLOSURE_SPLINE,    TRUE,  FALSE},
  {"non_limber/chebyshev",  NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV, TRUE,  FALSE},

  /* Two disjoint bumps: the only case where the per-component accessors have more than
   * one component to split, and where the panel edges have to be merged across them.
   * This is also the case that lands the last point of the sampling grid on a k range
   * degenerate to within rounding -- see NC_XCOR_KERNEL_COMPONENT_K_RANGE_MIN_WIDTH. */
  {"multi/chebyshev",      NC_XCOR_KERNEL_CLOSURE_CHEBYSHEV, TRUE,  TRUE },
};

static NcXcorKernel *
_build_multi (TestNcXcorKernelIntegrand *test)
{
  NcmVector *chi_mean  = ncm_vector_new (2);
  NcmVector *chi_sigma = ncm_vector_new (2);
  NcmVector *weight    = ncm_vector_new (2);
  NcXcorKernel *xclk;

  ncm_vector_set (chi_mean, 0, 200.0);
  ncm_vector_set (chi_mean, 1, 500.0);
  ncm_vector_set (chi_sigma, 0, 30.0);
  ncm_vector_set (chi_sigma, 1, 30.0);
  ncm_vector_set (weight, 0, 1.0);
  ncm_vector_set (weight, 1, 0.6);

  xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_multi_new_full (test->dist, test->ps, chi_mean,
                                                                 chi_sigma, weight, 4.0, test->sbi));

  ncm_vector_free (chi_mean);
  ncm_vector_free (chi_sigma);
  ncm_vector_free (weight);

  return xclk;
}

static void
test_nc_xcor_kernel_integrand_new (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc = pdata;
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

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

  if (tc->multi_comp)
    test->xclk = _build_multi (test);
  else
    test->xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new_full (test->dist, test->ps,
                                                                          TEST_CHI_LOWER,
                                                                          TEST_CHI_UPPER,
                                                                          test->sbi));

  /* Cap the adaptive apparatus rather than let it converge: one border expansion, few
   * iterations, loose tolerances, a low panel order. Each branch is still entered. */
  nc_xcor_kernel_set_max_border_expansions (test->xclk, 1);
  nc_xcor_kernel_set_max_iter (test->xclk, 4);
  nc_xcor_kernel_set_reltol (test->xclk, 1.0e-3);
  nc_xcor_kernel_set_scaled_abstol (test->xclk, 1.0e-4);
  nc_xcor_kernel_set_panel_order_cap (test->xclk, 12);
  nc_xcor_kernel_set_lmax (test->xclk, 16);

  /* The threshold is "use Limber for l >= this": 0 means always, -1 means never. Both
   * settings need an integrator, which _new_full supplied. */
  nc_xcor_kernel_set_l_limber (test->xclk, tc->non_limber ? -1 : 0);

  nc_xcor_kernel_prepare (test->xclk, cosmo);
}

static void
test_nc_xcor_kernel_integrand_free (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  ncm_sbessel_integrator_free (test->sbi);
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_kernel_free, test->xclk);
}

static void
test_nc_xcor_kernel_integrand_eval (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc           = pdata;
  NcXcorKernelIntegrand *integ = nc_xcor_kernel_get_eval (test->xclk, test->cosmo, TEST_ELL, tc->closure);
  const guint len              = nc_xcor_kernel_integrand_get_len (integ);
  gdouble k_min, k_max;
  gdouble *W = g_new0 (gdouble, len);
  guint i, j;

  g_assert_cmpuint (len, >, 0);

  nc_xcor_kernel_integrand_get_range (integ, &k_min, &k_max);
  g_assert_true (gsl_finite (k_min));
  g_assert_true (gsl_finite (k_max));
  g_assert_cmpfloat (k_min, >, 0.0);
  g_assert_cmpfloat (k_min, <, k_max);

  for (j = 0; j <= 8; j++)
  {
    const gdouble k = k_min * pow (k_max / k_min, j / 8.0);
    GArray *arr;

    nc_xcor_kernel_integrand_eval (integ, k, W);

    for (i = 0; i < len; i++)
      g_assert_true (gsl_finite (W[i]));

    /* The array form allocates its own storage; it must agree entry for entry. */
    arr = nc_xcor_kernel_integrand_eval_array (integ, k);
    g_assert_cmpuint (arr->len, ==, len);

    for (i = 0; i < len; i++)
      ncm_assert_cmpdouble_e (g_array_index (arr, gdouble, i), ==, W[i], 1.0e-15, 0.0);

    g_array_unref (arr);
  }

  g_free (W);
  nc_xcor_kernel_integrand_unref (integ);
}

static void
test_nc_xcor_kernel_integrand_comps (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc           = pdata;
  NcXcorKernelIntegrand *integ = nc_xcor_kernel_get_eval (test->xclk, test->cosmo, TEST_ELL, tc->closure);
  const guint len              = nc_xcor_kernel_integrand_get_len (integ);
  gdouble *W                   = g_new0 (gdouble, len);
  gdouble *Wc                  = g_new0 (gdouble, len);
  gdouble k_min, k_max;
  guint i;

  nc_xcor_kernel_integrand_get_range (integ, &k_min, &k_max);

  for (i = 0; i < len; i++)
  {
    gdouble ck_min, ck_max;

    nc_xcor_kernel_integrand_get_range_comp (integ, i, &ck_min, &ck_max);

    g_assert_true (gsl_finite (ck_min));
    g_assert_true (gsl_finite (ck_max));
    g_assert_cmpfloat (ck_min, <, ck_max);

    /* Component ranges are what the block integrators union into the full range, so
     * none may stick out of it. */
    g_assert_cmpfloat (ck_min, >=, k_min);
    g_assert_cmpfloat (ck_max, <=, k_max);
  }

  /* Evaluating all components at once and one at a time must give the same thing: the
   * offset/len form is the hot path and the whole-array form is the reference. */
  nc_xcor_kernel_integrand_eval (integ, sqrt (k_min * k_max), W);

  for (i = 0; i < len; i++)
  {
    nc_xcor_kernel_integrand_eval_comps (integ, sqrt (k_min * k_max), i, 1, Wc);
    ncm_assert_cmpdouble_e (Wc[0], ==, W[i], 1.0e-15, 0.0);
  }

  g_free (W);
  g_free (Wc);
  nc_xcor_kernel_integrand_unref (integ);
}

/* Every accessor beyond eval/get_range is optional -- an integrand that does not offer
 * one answers 0 or FALSE rather than failing. What is not optional is that an integrator
 * can reach the representation by *some* route, and that whichever routes are offered
 * describe the same thing. */
static void
test_nc_xcor_kernel_integrand_panels (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc           = pdata;
  NcXcorKernelIntegrand *integ = nc_xcor_kernel_get_eval (test->xclk, test->cosmo, TEST_ELL, tc->closure);
  NcmMatrix *coeffs            = NULL;
  NcmVector *knots             = nc_xcor_kernel_integrand_peek_knots (integ);
  gdouble k_min, k_max, s_min, s_max;
  guint n_panels, i;
  gboolean has_spectral;

  nc_xcor_kernel_integrand_get_range (integ, &k_min, &k_max);
  n_panels     = nc_xcor_kernel_integrand_get_n_panels (integ);
  has_spectral = nc_xcor_kernel_integrand_peek_spectral (integ, &coeffs, &s_min, &s_max);

  g_assert_true (n_panels > 0 || has_spectral || knots != NULL);

  /* The panels tile the range in order and without gaps: the block integrators walk them
   * assuming exactly that, and a gap is silently lost integrand. */
  if (n_panels > 0)
  {
    gdouble prev_b = k_min;

    for (i = 0; i < n_panels; i++)
    {
      gdouble a, b;

      coeffs = NULL;
      nc_xcor_kernel_integrand_peek_panel (integ, i, &coeffs, &a, &b);

      g_assert_nonnull (coeffs);
      g_assert_cmpuint (ncm_matrix_nrows (coeffs), >, 0);
      g_assert_cmpfloat (a, <, b);
      ncm_assert_cmpdouble_e (a, ==, prev_b, 1.0e-12, 0.0);

      prev_b = b;
    }

    ncm_assert_cmpdouble_e (prev_b, ==, k_max, 1.0e-12, 0.0);
  }

  if (knots != NULL)
    g_assert_cmpuint (ncm_vector_len (knots), >, 1);

  /* Restricting to a sub-interval yields a series over it, when the representation
   * supports being restricted at all. */
  coeffs = NULL;

  if (nc_xcor_kernel_integrand_restrict (integ,
                                         k_min + 0.25 * (k_max - k_min),
                                         k_min + 0.75 * (k_max - k_min),
                                         &coeffs))
  {
    g_assert_nonnull (coeffs);
    g_assert_cmpuint (ncm_matrix_nrows (coeffs), >, 0);
  }

  nc_xcor_kernel_integrand_unref (integ);
}

/* The single-series spectral form exists exactly when the Chebyshev representation did
 * not have to split, so it is one panel spanning the whole range. */
static void
test_nc_xcor_kernel_integrand_spectral (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc           = pdata;
  NcXcorKernelIntegrand *integ = nc_xcor_kernel_get_eval (test->xclk, test->cosmo, TEST_ELL, tc->closure);
  NcmMatrix *coeffs            = NULL;
  gdouble k_min, k_max;

  if (nc_xcor_kernel_integrand_peek_spectral (integ, &coeffs, &k_min, &k_max))
  {
    g_assert_nonnull (coeffs);
    g_assert_cmpuint (ncm_matrix_nrows (coeffs), >, 0);
    g_assert_cmpuint (ncm_matrix_ncols (coeffs), >, 0);
    g_assert_cmpfloat (k_min, <, k_max);
    g_assert_cmpuint (nc_xcor_kernel_integrand_get_n_panels (integ), ==, 1);
  }
  else
  {
    g_assert_cmpuint (nc_xcor_kernel_integrand_get_n_panels (integ), !=, 1);
  }

  nc_xcor_kernel_integrand_unref (integ);
}

static void
test_nc_xcor_kernel_integrand_tolerances (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc           = pdata;
  NcXcorKernelIntegrand *integ = nc_xcor_kernel_get_eval (test->xclk, test->cosmo, TEST_ELL, tc->closure);
  NcXcorKernelIntegrand *ref   = nc_xcor_kernel_integrand_ref (integ);
  NcmMatrix *residuals         = ncm_matrix_new (2, 3);

  nc_xcor_kernel_integrand_set_tolerances (integ, 1.0e-7, 1.0e-5);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_integrand_get_reltol (integ), ==, 1.0e-7, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_integrand_get_scaled_abstol (integ), ==, 1.0e-5, 1.0e-15, 0.0);

  ncm_matrix_set_all (residuals, 1.0e-9);
  nc_xcor_kernel_integrand_set_residuals (integ, residuals);
  g_assert_true (nc_xcor_kernel_integrand_peek_residuals (integ) == residuals);

  nc_xcor_kernel_integrand_clear (&ref);
  g_assert_true (ref == NULL);

  ncm_matrix_free (residuals);
  nc_xcor_kernel_integrand_unref (integ);
}

/* The vectorized form builds one representation for a whole multipole block, with one
 * entry per multipole rather than per component.
 *
 * Only the shape is checked, not agreement with the single-l form. The fixture caps
 * max_iter at 4, so neither fit is converged -- the sampler says as much -- and two
 * unconverged fits of the same function are not approximations of it that can be
 * required to agree: on the two-bump kernel they differ by more than 100% at the same k.
 * Agreement between the block and single-l paths is a statement about accuracy and needs
 * converged fits, so it is checked in tests/python/nc/xcor at real tolerances. */
static void
test_nc_xcor_kernel_integrand_vectorized (TestNcXcorKernelIntegrand *test, gconstpointer pdata)
{
  const TestCase *tc           = pdata;
  NcXcorKernelIntegrand *block = nc_xcor_kernel_get_eval_vectorized (test->xclk, test->cosmo,
                                                                     TEST_ELL, TEST_ELL + 1, tc->closure);
  NcXcorKernelIntegrand *one = nc_xcor_kernel_get_eval (test->xclk, test->cosmo, TEST_ELL, tc->closure);
  const guint block_len      = nc_xcor_kernel_integrand_get_len (block);
  gdouble *Wb                = g_new0 (gdouble, block_len);
  gdouble k_min, k_max, b_min, b_max;
  guint i;

  g_assert_cmpuint (nc_xcor_kernel_integrand_get_len (one), ==, 1);
  g_assert_cmpuint (block_len, ==, 2);

  nc_xcor_kernel_integrand_get_range (one, &k_min, &k_max);
  nc_xcor_kernel_integrand_get_range (block, &b_min, &b_max);

  g_assert_cmpfloat (b_min, >, 0.0);
  g_assert_cmpfloat (b_min, <, b_max);

  nc_xcor_kernel_integrand_eval (block, sqrt (b_min * b_max), Wb);

  for (i = 0; i < block_len; i++)
    g_assert_true (gsl_finite (Wb[i]));

  g_free (Wb);
  nc_xcor_kernel_integrand_unref (block);
  nc_xcor_kernel_integrand_unref (one);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorKernelIntegrand *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"eval",       &test_nc_xcor_kernel_integrand_eval      },
  {"comps",      &test_nc_xcor_kernel_integrand_comps     },
  {"panels",     &test_nc_xcor_kernel_integrand_panels    },
  {"spectral",   &test_nc_xcor_kernel_integrand_spectral  },
  {"tolerances", &test_nc_xcor_kernel_integrand_tolerances},
  {"vectorized", &test_nc_xcor_kernel_integrand_vectorized},
};

gint
main (gint argc, gchar *argv[])
{
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (test_cases); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
    {
      gchar *path = g_strdup_printf ("/nc/xcor/kernel/integrand/%s/%s",
                                     test_cases[i].name, test_checks[j].name);

      g_test_add (path, TestNcXcorKernelIntegrand, &test_cases[i],
                  &test_nc_xcor_kernel_integrand_new,
                  test_checks[j].func,
                  &test_nc_xcor_kernel_integrand_free);

      g_free (path);
    }
  }

  g_test_run ();
}


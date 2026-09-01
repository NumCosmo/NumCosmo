/***************************************************************************
 *            test_nc_xcor_kernel_analytic.c
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
 * The seven analytic windows are library objects, not test scaffolding: they are
 * installed, registered and serializable, and every kernel-space check in the suite is
 * built on one. What is checked here is that contract -- construction, properties,
 * component supports, the window vanishing outside its own support, serialization -- not
 * whether the closed forms are right, which is what tests/python/nc/xcor covers against
 * certified values.
 *
 * The cosmology is NcHICosmoDEXcdm with an analytic power spectrum, so nothing here runs
 * CLASS or any transfer function.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/* Comoving distances (Mpc) the windows are placed at: well inside the distance object's
 * z <= 6 range for this cosmology, and narrow enough that a support grid is cheap. */
#define TEST_CHI_MEAN 1500.0
#define TEST_CHI_SIGMA 300.0
#define TEST_N_SIGMA 4.0
#define TEST_CHI_LOWER 500.0
#define TEST_CHI_UPPER 2500.0

#define TEST_GRID 64

typedef struct _TestNcXcorKernelAnalytic
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcXcorKernel *xclk;
} TestNcXcorKernelAnalytic;

typedef struct _TestProp
{
  const gchar *name;
  gdouble value;
} TestProp;

/**
 * TestKernel:
 * @name: path component naming the window
 * @build: builds the window on the fixture's distance and power spectrum
 * @n_comps: how many components the window is made of
 * @build_full: the _new_full form, taking an explicit sBessel integrator, or NULL
 * @check_getters: asserts this window's typed getters return the constructor arguments
 * @props: constructor arguments, by property name, terminated by a NULL name
 *
 * One analytic window's worth of what the shared checks need to know about it.
 */
typedef struct _TestKernel
{
  const gchar *name;

  NcXcorKernel *(*build) (NcDistance *dist, NcmPowspec *ps);
  NcXcorKernel *(*build_full) (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi);
  void (*check_getters) (NcXcorKernel *xclk);

  guint n_comps;
  const TestProp *props;
} TestKernel;

static NcXcorKernel *
_build_gauss (NcDistance *dist, NcmPowspec *ps)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (dist, ps, TEST_CHI_MEAN, TEST_CHI_SIGMA, TEST_N_SIGMA));
}

static NcXcorKernel *
_build_tophat (NcDistance *dist, NcmPowspec *ps)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new (dist, ps, TEST_CHI_LOWER, TEST_CHI_UPPER));
}

static NcXcorKernel *
_build_tophat_smooth (NcDistance *dist, NcmPowspec *ps)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_smooth_new (dist, ps, 1000.0, 2000.0, 150.0, 6.0));
}

static NcXcorKernel *
_build_student_t (NcDistance *dist, NcmPowspec *ps)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_student_t_new (dist, ps, TEST_CHI_MEAN, 200.0, 2.0, 6.0));
}

static NcXcorKernel *
_build_power_exp (NcDistance *dist, NcmPowspec *ps)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_power_exp_new (dist, ps, 1200.0, 2.0, 1.5, 50.0, 4000.0));
}

static NcXcorKernel *
_build_lensing (NcDistance *dist, NcmPowspec *ps)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_lensing_new (dist, ps, 50.0, 2000.0, 3000.0));
}

/* The multimodal window merges bumps whose n-sigma intervals overlap, so the same
 * constructor yields one component or two depending only on where the bumps sit. Both
 * are built here: that merge is the one piece of real logic in the analytic set. */
static NcXcorKernel *
_build_multi_full (NcDistance *dist, NcmPowspec *ps, gdouble mean0, gdouble mean1,
                   gdouble sigma0, gdouble sigma1)
{
  NcmVector *chi_mean  = ncm_vector_new (2);
  NcmVector *chi_sigma = ncm_vector_new (2);
  NcmVector *weight    = ncm_vector_new (2);
  NcXcorKernel *xclk;

  ncm_vector_set (chi_mean, 0, mean0);
  ncm_vector_set (chi_mean, 1, mean1);
  ncm_vector_set (chi_sigma, 0, sigma0);
  ncm_vector_set (chi_sigma, 1, sigma1);
  ncm_vector_set (weight, 0, 1.0);
  ncm_vector_set (weight, 1, 0.6);

  xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_multi_new (dist, ps, chi_mean, chi_sigma, weight, TEST_N_SIGMA));

  ncm_vector_free (chi_mean);
  ncm_vector_free (chi_sigma);
  ncm_vector_free (weight);

  return xclk;
}

static NcXcorKernel *
_build_multi (NcDistance *dist, NcmPowspec *ps)
{
  return _build_multi_full (dist, ps, 1000.0, 1600.0, 300.0, 300.0);
}

static NcXcorKernel *
_build_multi_disjoint (NcDistance *dist, NcmPowspec *ps)
{
  return _build_multi_full (dist, ps, 600.0, 2600.0, 100.0, 150.0);
}

static NcXcorKernel *
_build_gauss_full (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new_full (dist, ps, TEST_CHI_MEAN, TEST_CHI_SIGMA, TEST_N_SIGMA, sbi));
}

static void
_getters_gauss (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticGauss *k = NC_XCOR_KERNEL_ANALYTIC_GAUSS (xclk);

  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_gauss_get_chi_mean (k), ==, TEST_CHI_MEAN, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_gauss_get_chi_sigma (k), ==, TEST_CHI_SIGMA, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_gauss_get_n_sigma (k), ==, TEST_N_SIGMA, 1.0e-15, 0.0);
}

static NcXcorKernel *
_build_tophat_full (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_new_full (dist, ps, TEST_CHI_LOWER, TEST_CHI_UPPER, sbi));
}

static void
_getters_tophat (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticTophat *k = NC_XCOR_KERNEL_ANALYTIC_TOPHAT (xclk);

  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_tophat_get_chi_lower (k), ==, TEST_CHI_LOWER, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_tophat_get_chi_upper (k), ==, TEST_CHI_UPPER, 1.0e-15, 0.0);
}

static NcXcorKernel *
_build_tophat_smooth_full (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_tophat_smooth_new_full (dist, ps, 1000.0, 2000.0, 150.0, 6.0, sbi));
}

static void
_getters_tophat_smooth (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticTophatSmooth *k = NC_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH (xclk);

  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_tophat_smooth_get_chi_lower (k), ==, 1000.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_tophat_smooth_get_chi_upper (k), ==, 2000.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_tophat_smooth_get_chi_sigma (k), ==, 150.0, 1.0e-15, 0.0);
}

static NcXcorKernel *
_build_student_t_full (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_student_t_new_full (dist, ps, TEST_CHI_MEAN, 200.0, 2.0, 6.0, sbi));
}

static void
_getters_student_t (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticStudentT *k = NC_XCOR_KERNEL_ANALYTIC_STUDENT_T (xclk);

  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_student_t_get_chi_mean (k), ==, TEST_CHI_MEAN, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_student_t_get_chi_scale (k), ==, 200.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_student_t_get_nu (k), ==, 2.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_student_t_get_n_scale (k), ==, 6.0, 1.0e-15, 0.0);

  /* The truncated tail's missing mass: a diagnostic the window itself computes, so it
   * has to be a sane fraction rather than merely finite. */
  g_assert_cmpfloat (nc_xcor_kernel_analytic_student_t_get_tail_mass (k), >=, 0.0);
  g_assert_cmpfloat (nc_xcor_kernel_analytic_student_t_get_tail_mass (k), <, 1.0);
}

static NcXcorKernel *
_build_power_exp_full (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_power_exp_new_full (dist, ps, 1200.0, 2.0, 1.5, 50.0, 4000.0, sbi));
}

static void
_getters_power_exp (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticPowerExp *k = NC_XCOR_KERNEL_ANALYTIC_POWER_EXP (xclk);

  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_power_exp_get_chi_scale (k), ==, 1200.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_power_exp_get_alpha (k), ==, 2.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_power_exp_get_beta (k), ==, 1.5, 1.0e-15, 0.0);
}

static NcXcorKernel *
_build_lensing_full (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_analytic_lensing_new_full (dist, ps, 50.0, 2000.0, 3000.0, sbi));
}

static void
_getters_lensing (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticLensing *k = NC_XCOR_KERNEL_ANALYTIC_LENSING (xclk);

  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_lensing_get_chi_source_lower (k), ==, 2000.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_lensing_get_chi_source_upper (k), ==, 3000.0, 1.0e-15, 0.0);
}

static NcXcorKernel *
_build_multi_full_sbi (NcDistance *dist, NcmPowspec *ps, NcmSBesselIntegrator *sbi)
{
  NcmVector *chi_mean  = ncm_vector_new (2);
  NcmVector *chi_sigma = ncm_vector_new (2);
  NcmVector *weight    = ncm_vector_new (2);
  NcXcorKernel *xclk;

  ncm_vector_set (chi_mean, 0, 1000.0);
  ncm_vector_set (chi_mean, 1, 1600.0);
  ncm_vector_set (chi_sigma, 0, 300.0);
  ncm_vector_set (chi_sigma, 1, 300.0);
  ncm_vector_set (weight, 0, 1.0);
  ncm_vector_set (weight, 1, 0.6);

  xclk = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_multi_new_full (dist, ps, chi_mean, chi_sigma,
                                                                 weight, TEST_N_SIGMA, sbi));

  ncm_vector_free (chi_mean);
  ncm_vector_free (chi_sigma);
  ncm_vector_free (weight);

  return xclk;
}

static void
_getters_multi (NcXcorKernel *xclk)
{
  NcXcorKernelAnalyticMulti *k = NC_XCOR_KERNEL_ANALYTIC_MULTI (xclk);

  g_assert_cmpuint (nc_xcor_kernel_analytic_multi_get_n_bumps (k), ==, 2);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_analytic_multi_get_n_sigma (k), ==, TEST_N_SIGMA, 1.0e-15, 0.0);

  /* Bump count and vector length are set independently; a mismatch would index past the
   * end of these on the first evaluation. */
  g_assert_cmpuint (ncm_vector_len (nc_xcor_kernel_analytic_multi_peek_chi_mean (k)), ==, 2);
  g_assert_cmpuint (ncm_vector_len (nc_xcor_kernel_analytic_multi_peek_chi_sigma (k)), ==, 2);
  g_assert_cmpuint (ncm_vector_len (nc_xcor_kernel_analytic_multi_peek_weight (k)), ==, 2);
}

/* Uncrustify reflows nested array-in-struct initializers into something unreadable, so
 * each window's property list is its own named array and the descriptor points at it. */
static const TestProp _props_gauss[] = {
  {"chi-mean",  TEST_CHI_MEAN},
  {"chi-sigma", TEST_CHI_SIGMA},
  {"n-sigma",   TEST_N_SIGMA},
  {NULL,        0.0}
};

static const TestProp _props_tophat[] = {
  {"chi-lower", TEST_CHI_LOWER},
  {"chi-upper", TEST_CHI_UPPER},
  {NULL,        0.0}
};

static const TestProp _props_tophat_smooth[] = {
  {"chi-lower", 1000.0},
  {"chi-upper", 2000.0},
  {"chi-sigma", 150.0},
  {"n-sigma",   6.0},
  {NULL,        0.0}
};

static const TestProp _props_student_t[] = {
  {"chi-mean",  TEST_CHI_MEAN},
  {"chi-scale", 200.0},
  {"nu",        2.0},
  {"n-scale",   6.0},
  {NULL,        0.0}
};

static const TestProp _props_power_exp[] = {
  {"chi-scale", 1200.0},
  {"alpha",     2.0},
  {"beta",      1.5},
  {"chi-lower", 50.0},
  {"chi-upper", 4000.0},
  {NULL,        0.0}
};

static const TestProp _props_lensing[] = {
  {"chi-lower",        50.0},
  {"chi-source-lower", 2000.0},
  {"chi-source-upper", 3000.0},
  {NULL,               0.0}
};

static const TestProp _props_multi[] = {
  {"n-sigma", TEST_N_SIGMA},
  {NULL,      0.0}
};

static const TestKernel test_kernels[] = {
  {"gauss",          &_build_gauss,          &_build_gauss_full,         &_getters_gauss,         1, _props_gauss},
  {"tophat",         &_build_tophat,         &_build_tophat_full,        &_getters_tophat,        1, _props_tophat},
  {"tophat_smooth",  &_build_tophat_smooth,  &_build_tophat_smooth_full, &_getters_tophat_smooth, 1, _props_tophat_smooth},
  {"student_t",      &_build_student_t,      &_build_student_t_full,     &_getters_student_t,     1, _props_student_t},
  {"power_exp",      &_build_power_exp,      &_build_power_exp_full,     &_getters_power_exp,     1, _props_power_exp},

  /* Two components, one per branch of the lensing shape: they meet at the source-bin
   * edge with a kink no single Chebyshev panel fits. */
  {"lensing",        &_build_lensing,        &_build_lensing_full,       &_getters_lensing,       2, _props_lensing},
  {"multi",          &_build_multi,          &_build_multi_full_sbi,     &_getters_multi,         1, _props_multi},
  {"multi_disjoint", &_build_multi_disjoint, NULL,                       &_getters_multi,         2, _props_multi},
};

#define TEST_N_KERNELS (G_N_ELEMENTS (test_kernels))

static void
test_nc_xcor_kernel_analytic_new (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;
  NcHICosmo *cosmo     = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist     = nc_distance_new (6.0);
  NcmPowspec *ps       = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                                NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,      70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, 0.255);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B, 0.045);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W, -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", 0.0, NULL);

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  test->cosmo = cosmo;
  test->dist  = dist;
  test->ps    = ps;
  test->xclk  = tk->build (dist, ps);

  g_assert_true (NC_IS_XCOR_KERNEL (test->xclk));
  g_assert_true (NC_IS_XCOR_KERNEL_RADIAL (test->xclk));

  nc_xcor_kernel_prepare (test->xclk, cosmo);
}

static void
test_nc_xcor_kernel_analytic_free (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_kernel_free, test->xclk);
}

static void
test_nc_xcor_kernel_analytic_basic (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  NcXcorKernel *xclk2 = nc_xcor_kernel_ref (test->xclk);

  g_assert_true (nc_xcor_kernel_peek_dist (test->xclk) == test->dist);
  g_assert_true (nc_xcor_kernel_peek_powspec (test->xclk) == test->ps);

  g_assert_cmpuint (nc_xcor_kernel_obs_len (test->xclk), >, 0);

  nc_xcor_kernel_clear (&xclk2);
  g_assert_true (xclk2 == NULL);
}

static void
test_nc_xcor_kernel_analytic_props (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;
  guint i;

  for (i = 0; tk->props[i].name != NULL; i++)
  {
    gdouble value = 0.0;

    g_object_get (test->xclk, tk->props[i].name, &value, NULL);
    ncm_assert_cmpdouble_e (value, ==, tk->props[i].value, 1.0e-15, 0.0);
  }

  g_assert_cmpuint (i, >, 0);
}

static void
test_nc_xcor_kernel_analytic_support (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk       = pdata;
  NcXcorKernelRadial *radial = NC_XCOR_KERNEL_RADIAL (test->xclk);
  gdouble chi_min, chi_max;
  guint i;

  g_assert_cmpuint (nc_xcor_kernel_radial_get_n_comps (radial), ==, tk->n_comps);

  nc_xcor_kernel_radial_get_support (radial, &chi_min, &chi_max);

  g_assert_true (gsl_finite (chi_min));
  g_assert_true (gsl_finite (chi_max));
  g_assert_cmpfloat (chi_min, >=, 0.0);
  g_assert_cmpfloat (chi_min, <, chi_max);

  /* The total support is the union of the components', so every component must sit
   * inside it -- this is what the C_ell integrators use to place their panels. */
  for (i = 0; i < tk->n_comps; i++)
  {
    gdouble lo, hi;

    nc_xcor_kernel_radial_get_comp_support (radial, i, &lo, &hi);

    g_assert_true (gsl_finite (lo));
    g_assert_true (gsl_finite (hi));
    g_assert_cmpfloat (lo, <, hi);
    g_assert_cmpfloat (lo, >=, chi_min);
    g_assert_cmpfloat (hi, <=, chi_max);
  }
}

static void
test_nc_xcor_kernel_analytic_window (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk       = pdata;
  NcXcorKernelRadial *radial = NC_XCOR_KERNEL_RADIAL (test->xclk);
  gdouble chi_min, chi_max;
  gdouble mass = 0.0;
  guint i, j;

  nc_xcor_kernel_radial_get_support (radial, &chi_min, &chi_max);

  for (j = 0; j <= TEST_GRID; j++)
  {
    const gdouble chi = chi_min + (chi_max - chi_min) * j / (gdouble) TEST_GRID;
    const gdouble W   = nc_xcor_kernel_radial_eval_W (radial, chi);

    g_assert_true (gsl_finite (W));
    g_assert_cmpfloat (W, >=, 0.0);

    mass += W;
  }

  /* A window that is zero everywhere would pass every check above and integrate to
   * nothing downstream. */
  g_assert_cmpfloat (mass, >, 0.0);

  /* Outside its own support a component must be exactly zero: the integrators clamp to
   * these bounds and never sample beyond them, so a non-zero tail there is silently
   * dropped mass rather than an approximation error. */
  for (i = 0; i < tk->n_comps; i++)
  {
    gdouble lo, hi;

    nc_xcor_kernel_radial_get_comp_support (radial, i, &lo, &hi);

    g_assert_cmpfloat (nc_xcor_kernel_radial_eval_W_comp (radial, i, lo - 1.0e-2 * (hi - lo)), ==, 0.0);
    g_assert_cmpfloat (nc_xcor_kernel_radial_eval_W_comp (radial, i, hi + 1.0e-2 * (hi - lo)), ==, 0.0);
  }
}

static void
test_nc_xcor_kernel_analytic_z_range (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  gdouble zmin, zmax, zmid;

  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);

  g_assert_true (gsl_finite (zmin));
  g_assert_true (gsl_finite (zmax));
  g_assert_true (gsl_finite (zmid));
  g_assert_cmpfloat (zmin, >=, 0.0);
  g_assert_cmpfloat (zmin, <=, zmid);
  g_assert_cmpfloat (zmid, <=, zmax);
}

static void
test_nc_xcor_kernel_analytic_limber (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const gint l = 10;
  gdouble zmin, zmax, zmid;
  guint j;

  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);

  g_assert_true (gsl_finite (nc_xcor_kernel_eval_limber_z_prefactor (test->xclk, test->cosmo, l)));

  for (j = 1; j < 8; j++)
  {
    const gdouble z = zmin + (zmax - zmin) * j / 8.0;

    g_assert_true (gsl_finite (nc_xcor_kernel_eval_limber_z_full (test->xclk, test->cosmo, z, test->dist, l)));
  }
}

static void
test_nc_xcor_kernel_analytic_serialize (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk       = pdata;
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcXcorKernel *dup          = NC_XCOR_KERNEL (ncm_serialize_dup_obj (ser, G_OBJECT (test->xclk)));
  NcXcorKernelRadial *radial = NC_XCOR_KERNEL_RADIAL (test->xclk);
  NcXcorKernelRadial *rdup   = NC_XCOR_KERNEL_RADIAL (dup);
  gdouble chi_min, chi_max;
  guint i;

  g_assert_true (G_OBJECT_TYPE (dup) == G_OBJECT_TYPE (test->xclk));
  g_assert_cmpuint (nc_xcor_kernel_radial_get_n_comps (rdup), ==, tk->n_comps);

  for (i = 0; tk->props[i].name != NULL; i++)
  {
    gdouble value = 0.0;

    g_object_get (dup, tk->props[i].name, &value, NULL);
    ncm_assert_cmpdouble_e (value, ==, tk->props[i].value, 1.0e-15, 0.0);
  }

  /* The duplicate carries its own distance and power spectrum, so it must be prepared
   * before its window can be compared with the original's. */
  nc_xcor_kernel_prepare (dup, test->cosmo);

  nc_xcor_kernel_radial_get_support (radial, &chi_min, &chi_max);

  for (i = 0; i <= TEST_GRID; i++)
  {
    const gdouble chi = chi_min + (chi_max - chi_min) * i / (gdouble) TEST_GRID;

    ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_eval_W (rdup, chi), ==,
                            nc_xcor_kernel_radial_eval_W (radial, chi), 1.0e-14, 0.0);
  }

  nc_xcor_kernel_free (dup);
  ncm_serialize_free (ser);
}

static void
test_nc_xcor_kernel_analytic_getters (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;

  tk->check_getters (test->xclk);
}

/* _new_full() is _new() plus an explicit sBessel integrator, which is what moves the
 * kernel off the Limber tier. The window must not depend on that choice. */
static void
test_nc_xcor_kernel_analytic_new_full (TestNcXcorKernelAnalytic *test, gconstpointer pdata)
{
  const TestKernel *tk      = pdata;
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));
  NcXcorKernel *full;
  gdouble chi_min, chi_max;
  guint i;

  if (tk->build_full == NULL)
  {
    g_test_skip ("no _new_full form for this window");
    ncm_sbessel_integrator_free (sbi);

    return;
  }

  full = tk->build_full (test->dist, test->ps, sbi);

  g_assert_true (nc_xcor_kernel_peek_integrator (full) == sbi);
  nc_xcor_kernel_prepare (full, test->cosmo);

  nc_xcor_kernel_radial_get_support (NC_XCOR_KERNEL_RADIAL (test->xclk), &chi_min, &chi_max);

  for (i = 0; i <= TEST_GRID; i++)
  {
    const gdouble chi = chi_min + (chi_max - chi_min) * i / (gdouble) TEST_GRID;

    ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_eval_W (NC_XCOR_KERNEL_RADIAL (full), chi), ==,
                            nc_xcor_kernel_radial_eval_W (NC_XCOR_KERNEL_RADIAL (test->xclk), chi),
                            1.0e-14, 0.0);
  }

  nc_xcor_kernel_free (full);
  ncm_sbessel_integrator_free (sbi);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorKernelAnalytic *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"basic",     &test_nc_xcor_kernel_analytic_basic},
  {"props",     &test_nc_xcor_kernel_analytic_props},
  {"getters",   &test_nc_xcor_kernel_analytic_getters},
  {"new_full",  &test_nc_xcor_kernel_analytic_new_full},
  {"support",   &test_nc_xcor_kernel_analytic_support},
  {"window",    &test_nc_xcor_kernel_analytic_window},
  {"z_range",   &test_nc_xcor_kernel_analytic_z_range},
  {"limber",    &test_nc_xcor_kernel_analytic_limber},
  {"serialize", &test_nc_xcor_kernel_analytic_serialize},
};

gint
main (gint argc, gchar *argv[])
{
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < TEST_N_KERNELS; i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
    {
      gchar *path = g_strdup_printf ("/nc/xcor/kernel/analytic/%s/%s",
                                     test_kernels[i].name, test_checks[j].name);

      g_test_add (path, TestNcXcorKernelAnalytic, &test_kernels[i],
                  &test_nc_xcor_kernel_analytic_new,
                  test_checks[j].func,
                  &test_nc_xcor_kernel_analytic_free);

      g_free (path);
    }
  }

  g_test_run ();
}


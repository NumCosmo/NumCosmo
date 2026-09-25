/***************************************************************************
 *            test_ncm_util.c
 *
 *  Tue July 30 19:54:22 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_lambert.h>

void test_ncm_util_projected_radius (void);
void test_ncm_util_complex (void);
void test_ncm_util_gaussian_int (void);
void test_ncm_util_gaussian_int_rng_two_sides (void);
void test_ncm_util_gaussian_int_rng_one_side (void);
void test_ncm_util_gaussian_int_nonunit (void);
void test_ncm_util_lambert_W0_ln (void);
void test_ncm_util_elementary (void);
void test_ncm_util_exprel (void);
void test_ncm_util_sinh (void);
void test_ncm_util_mln_1mIexpzA_1pIexpmzA (void);
void test_ncm_util_cmp (void);
void test_ncm_util_rational_coarse_double (void);
void test_ncm_util_fact_size (void);
void test_ncm_util_basename_fits (void);
void test_ncm_util_function_params (void);
void test_ncm_util_smooth_trans (void);
void test_ncm_util_sky_geometry (void);
void test_ncm_util_error (void);
void test_ncm_util_error_traps (void);
void test_ncm_util_error_set_null_subprocess (void);
void test_ncm_util_error_pile_up_subprocess (void);
void test_ncm_util_error_forward_null_subprocess (void);
void test_ncm_util_function_params_invalid_subprocess (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/util/projected_radius", test_ncm_util_projected_radius);
  g_test_add_func ("/ncm/util/complex", test_ncm_util_complex);
  g_test_add_func ("/ncm/util/gaussian_integral/sigmas", test_ncm_util_gaussian_int);
  g_test_add_func ("/ncm/util/gaussian_integral/rng/two_sides", test_ncm_util_gaussian_int_rng_two_sides);
  g_test_add_func ("/ncm/util/gaussian_integral/rng/one_side", test_ncm_util_gaussian_int_rng_one_side);
  g_test_add_func ("/ncm/util/gaussian_integral/nonunit", test_ncm_util_gaussian_int_nonunit);
  g_test_add_func ("/ncm/util/lambert_W0_ln", test_ncm_util_lambert_W0_ln);
  g_test_add_func ("/ncm/util/elementary", test_ncm_util_elementary);
  g_test_add_func ("/ncm/util/exprel", test_ncm_util_exprel);
  g_test_add_func ("/ncm/util/sinh", test_ncm_util_sinh);
  g_test_add_func ("/ncm/util/mln_1mIexpzA_1pIexpmzA", test_ncm_util_mln_1mIexpzA_1pIexpmzA);
  g_test_add_func ("/ncm/util/cmp", test_ncm_util_cmp);
  g_test_add_func ("/ncm/util/rational_coarse_double", test_ncm_util_rational_coarse_double);
  g_test_add_func ("/ncm/util/fact_size", test_ncm_util_fact_size);
  g_test_add_func ("/ncm/util/basename_fits", test_ncm_util_basename_fits);
  g_test_add_func ("/ncm/util/function_params", test_ncm_util_function_params);
  g_test_add_func ("/ncm/util/smooth_trans", test_ncm_util_smooth_trans);
  g_test_add_func ("/ncm/util/sky_geometry", test_ncm_util_sky_geometry);
  g_test_add_func ("/ncm/util/error", test_ncm_util_error);
  g_test_add_func ("/ncm/util/error_traps", test_ncm_util_error_traps);
  g_test_add_func ("/ncm/util/error/set_null/subprocess", test_ncm_util_error_set_null_subprocess);
  g_test_add_func ("/ncm/util/error/pile_up/subprocess", test_ncm_util_error_pile_up_subprocess);
  g_test_add_func ("/ncm/util/error/forward_null/subprocess", test_ncm_util_error_forward_null_subprocess);
  g_test_add_func ("/ncm/util/function_params/invalid/subprocess", test_ncm_util_function_params_invalid_subprocess);

  g_test_run ();
}

void
test_ncm_util_projected_radius (void)
{
  NcmRNG *rng    = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  gdouble ntests = 1000;
  guint i;

  for (i = 0; i < ntests; i++)
  {
    gdouble d     = ncm_rng_uniform_gen (rng, 0.0, 100.0);
    gdouble theta = ncm_rng_uniform_gen (rng, 0.0, M_PI);

    g_assert_cmpfloat (ncm_util_projected_radius (theta, d), >=, 0);
  }

  ncm_rng_free (rng);
}

void
test_ncm_util_complex (void)
{
  {
    NcmComplex *c1 = ncm_complex_new ();

    ncm_complex_set (c1, 1.0, 2.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 2.0);

    ncm_complex_set_zero (c1);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 0.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 0.0);

    ncm_complex_free (c1);
  }

  {
    NcmComplex c1 = NCM_COMPLEX_INIT (1.0 + 2.0 * I);

    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);
  }

  {
    NcmComplex *c1 = ncm_complex_new ();
    NcmComplex *c2 = ncm_complex_new ();
    NcmComplex *c3 = ncm_complex_new ();

    ncm_complex_set (c1, 1.0, 2.0);
    ncm_complex_set (c2, 3.0, 4.0);
    ncm_complex_set (c3, 5.0, 6.0);

    ncm_complex_mul_real (c1, 3.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 3.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 6.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (3.0, 6.0));

    ncm_complex_res_add_mul (c1, c2, c3);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -6.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 44.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (-6.0, 44.0));

    ncm_complex_res_add_mul_real (c1, c2, 2.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 0.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 52.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (0.0, 52.0));

    ncm_complex_res_mul (c1, c2);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -208.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 156.0);
    g_assert_cmpfloat (ncm_complex_Abs (c1), ==, hypot (-208.0, 156.0));

    ncm_complex_free (c1);
    ncm_complex_free (c2);
    ncm_complex_free (c3);
  }

  {
    NcmComplex c1;
    complex double z;

    ncm_complex_set_c (&c1, 1.0 + 2.0 * I);
    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);
    g_assert_cmpfloat (ncm_complex_Abs (&c1), ==, hypot (1.0, 2.0));

    z = ncm_complex_c (&c1);
    g_assert_cmpfloat (creal (z), ==, 1.0);
    g_assert_cmpfloat (cimag (z), ==, 2.0);
  }

  {
    NcmComplex c1  = NCM_COMPLEX_INIT (1.0 + 2.0 * I);
    NcmComplex *c2 = ncm_complex_dup (&c1);

    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);
    g_assert_cmpfloat (ncm_complex_Abs (&c1), ==, hypot (1.0, 2.0));

    g_assert_cmpfloat (ncm_complex_Re (c2), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (c2), ==, 2.0);

    ncm_complex_free (c2);
  }

  {
    NcmComplex *c1 = ncm_complex_new ();

    ncm_complex_clear (&c1);
    g_assert_null (c1);
  }
}

void
test_ncm_util_gaussian_int (void)
{
  const gdouble tol = 1.0e-14;
  gdouble sigma_int[10];
  gdouble sign;
  gint i;

  for (i = 0; i < 10; i++)
  {
    sigma_int[i] = gsl_cdf_chisq_P ((i + 1.0) * (i + 1.0), 1);
  }

  ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (-100.0, +100.0), ==, +1.0, tol, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (+100.0, -100.0), ==, -1.0, tol, 0.0);

  ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (-100.0, +100.0, &sign), ==, 0.0, tol, 0.0);
  g_assert_cmpfloat (sign, ==, +1.0);
  ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (+100.0, -100.0, &sign), ==, 0.0, tol, 0.0);
  g_assert_cmpfloat (sign, ==, -1.0);

  for (i = 0; i < 10; i++)
  {
    const gdouble x = i + 1.0;
    gdouble logtol;

    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (-x, +x), ==, +sigma_int[i], tol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (+x, -x), ==, -sigma_int[i], tol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (0.0, +x), ==, +0.5 * sigma_int[i], tol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (0.0, -x), ==, -0.5 * sigma_int[i], tol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (-x, 0.0), ==, +0.5 * sigma_int[i], tol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (+x, 0.0), ==, -0.5 * sigma_int[i], tol, 0.0);

    logtol = tol / fabs (sigma_int[i] - 1.0);

    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (-x, +x, &sign), ==, log (sigma_int[i]), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, +1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (+x, -x, &sign), ==, log (sigma_int[i]), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, -1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (0.0, +x, &sign), ==, log (0.5 * sigma_int[i]), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, +1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (0.0, -x, &sign), ==, log (0.5 * sigma_int[i]), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, -1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (-x, 0.0, &sign), ==, log (0.5 * sigma_int[i]), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, +1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (+x, 0.0, &sign), ==, log (0.5 * sigma_int[i]), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, -1.0);
  }
}

void
test_ncm_util_gaussian_int_rng_two_sides (void)
{
  const gdouble tol = 1.0e-14;
  gdouble sign, logtol;
  gint i;

  for (i = 0; i < 10; i++)
  {
    const gdouble xl       = g_test_rand_double_range (-10.0, 0.0);
    const gdouble xu       = g_test_rand_double_range (0.0, 10.0);
    const gdouble symint_l = gsl_cdf_chisq_P (xl * xl, 1);
    const gdouble symint_u = gsl_cdf_chisq_P (xu * xu, 1);
    const gdouble symint_d = 0.5 * (symint_u - symint_l);
    const gdouble int_val  = symint_l + symint_d;

    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (xl, xu), ==, +int_val, tol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (xu, xl), ==, -int_val, tol, 0.0);

    logtol = tol / fabs (int_val - 1.0);

    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (xl, xu, &sign), ==, log (int_val), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, +1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (xu, xl, &sign), ==, log (int_val), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, -1.0);
  }
}

void
test_ncm_util_gaussian_int_rng_one_side (void)
{
  const gdouble tol = 1.0e-14;
  gdouble sign, logtol, ltol;
  gint i;

  for (i = 0; i < 10; i++)
  {
    const gdouble xl       = g_test_rand_double_range (3.0, 8.0);
    const gdouble xu       = g_test_rand_double_range (xl, 10.0);
    const gdouble symint_l = gsl_cdf_chisq_P (xl * xl, 1);
    const gdouble symint_u = gsl_cdf_chisq_P (xu * xu, 1);
    const gdouble int_val  = 0.5 * (symint_u - symint_l);

    ltol = tol / fabs (symint_u / symint_l - 1.0);

    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (xl, xu), ==, +int_val, ltol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (xu, xl), ==, -int_val, ltol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (-xl, -xu), ==, -int_val, ltol, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_normal_gaussian_integral (-xu, -xl), ==, +int_val, ltol, 0.0);

    logtol = ltol / fabs (int_val - 1.0);

    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (xl, xu, &sign), ==, log (int_val), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, +1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (xu, xl, &sign), ==, log (int_val), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, -1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (-xl, -xu, &sign), ==, log (int_val), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, -1.0);
    ncm_assert_cmpdouble_e (ncm_util_log_normal_gaussian_integral (-xu, -xl, &sign), ==, log (int_val), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, +1.0);
  }
}

void
test_ncm_util_gaussian_int_nonunit (void)
{
  const gdouble tol = 1.0e-14;
  gdouble sign, logtol;
  gint i;

  for (i = 0; i < 10; i++)
  {
    const gdouble mu      = g_test_rand_double_range (-1.0, 1.0);
    const gdouble sigma   = g_test_rand_double_range (0.5, 2.0);
    const gdouble xl      = g_test_rand_double_range (-10.0, 10.0);
    const gdouble xu      = g_test_rand_double_range (-10.0, 10.0);
    const gdouble nonunit = ncm_util_gaussian_integral (xl, xu, mu, sigma);
    const gdouble unit    = ncm_util_normal_gaussian_integral ((xl - mu) / sigma, (xu - mu) / sigma);

    ncm_assert_cmpdouble_e (nonunit, ==, unit, tol, 0.0);

    logtol = tol / fabs (fabs (unit) - 1.0);

    ncm_assert_cmpdouble_e (ncm_util_log_gaussian_integral (xl, xu, mu, sigma, &sign), ==, log (fabs (unit)), logtol, 0.0);
    g_assert_cmpfloat (sign, ==, GSL_SIGN (unit));
  }
}

void
test_ncm_util_lambert_W0_ln (void)
{
  const gdouble ln_y[] = {-30.0, -1.0, 0.0, 1.0, 10.0, 300.0, 700.0, 709.0, 709.78, 710.0, 1.0e3, 1.0e5, 1.0e10, 1.0e300};
  guint i;

  /* W e^W = y, written as W + ln(W) = ln(y) for y > 0 */
  for (i = 0; i < G_N_ELEMENTS (ln_y); i++)
  {
    const gdouble W = ncm_util_lambert_W0_ln (ln_y[i]);

    ncm_assert_cmpdouble_e (W + log (W), ==, ln_y[i], 1.0e-15, 1.0e-14);
  }

  /* Continuity at GSL_LOG_DBL_MAX - 1, where Newton's method takes over */
  ncm_assert_cmpdouble_e (ncm_util_lambert_W0_ln (nextafter (GSL_LOG_DBL_MAX - 1.0, 0.0)), ==,
                          ncm_util_lambert_W0_ln (GSL_LOG_DBL_MAX - 1.0), 1.0e-14, 0.0);

  ncm_assert_cmpdouble_e (ncm_util_lambert_W0_ln (0.0), ==, gsl_sf_lambert_W0 (1.0), 1.0e-15, 0.0);
}

void
test_ncm_util_elementary (void)
{
  const gdouble xs[] = {-0.7, -0.2, 0.3, 1.1, 2.5, 3.0};
  guint i;

  /* Well-conditioned points: the direct expressions */
  for (i = 0; i < G_N_ELEMENTS (xs); i++)
  {
    const gdouble x = xs[i];

    ncm_assert_cmpdouble_e (ncm_util_sqrt1px_m1 (x), ==, sqrt (1.0 + x) - 1.0, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_ln1pexpx (x), ==, log (1.0 + exp (x)), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1pcosx (sin (x), cos (x)), ==, 1.0 + cos (x), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1mcosx (sin (x), cos (x)), ==, 1.0 - cos (x), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1psinx (sin (x), cos (x)), ==, 1.0 + sin (x), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1msinx (sin (x), cos (x)), ==, 1.0 - sin (x), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_cos2x (sin (x), cos (x)), ==, cos (2.0 * x), 1.0e-14, 1.0e-15);
  }

  /* Points of cancellation: Taylor series of the exact function */
  {
    const gdouble h  = 1.0e-4;
    const gdouble h2 = h * h;

    ncm_assert_cmpdouble_e (ncm_util_sqrt1px_m1 (h), ==, h / 2.0 - h2 / 8.0 + h * h2 / 16.0 - 5.0 * h2 * h2 / 128.0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1mcosx (sin (h), cos (h)), ==, h2 / 2.0 - h2 * h2 / 24.0, 1.0e-15, 0.0);
  }

  /* References from mpmath at 40 digits for the double inputs M_PI - h and M_PI_2 - h */
  {
    const gdouble h = 1.0e-4;

    ncm_assert_cmpdouble_e (ncm_util_1pcosx (sin (M_PI - h), cos (M_PI - h)), ==, 4.9999999958666829219e-9, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1msinx (sin (M_PI_2 - h), cos (M_PI_2 - h)), ==, 4.9999999958383552275e-9, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_1psinx (sin (h - M_PI_2), cos (h - M_PI_2)), ==, 4.9999999958383552275e-9, 1.0e-15, 0.0);
  }

  /* ln(1 + e^x) = x + ln(1 + e^{-x}) for large x, no overflow */
  ncm_assert_cmpdouble_e (ncm_util_ln1pexpx (30.0), ==, 30.0 + log1p (exp (-30.0)), 1.0e-15, 0.0);
  g_assert_cmpfloat (ncm_util_ln1pexpx (1000.0), ==, 1000.0);
  ncm_assert_cmpdouble_e (ncm_util_ln1pexpx (-40.0), ==, exp (-40.0), 1.0e-15, 0.0);
}

void
test_ncm_util_exprel (void)
{
  /* {x, exprel, first, second and third derivatives}, from mpmath at 40 digits */
  const gdouble ref[][5] = {
    {-1.0, 0.6321205588285576784, 0.26424111765711535681, 0.16060279414278839202, 0.11392894125692285447},
    {1.0, 1.7182818284590452354, 1.0, 0.71828182845904523536, 0.56343634308190952928},
    {2.0, 3.1945280494653251136, 2.0972640247326625568, 1.5972640247326625568, 1.2986320123663312784},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (ref); i++)
  {
    const gdouble x = ref[i][0];

    ncm_assert_cmpdouble_e (ncm_exprel (x), ==, ref[i][1], 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_d1exprel (x), ==, ref[i][2], 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_d2exprel (x), ==, ref[i][3], 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_d3exprel (x), ==, ref[i][4], 1.0e-15, 0.0);
  }

  /* The n-th derivative at zero is 1 / (n + 1) */
  ncm_assert_cmpdouble_e (ncm_exprel (0.0), ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_d1exprel (0.0), ==, 1.0 / 2.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_d2exprel (0.0), ==, 1.0 / 3.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_d3exprel (0.0), ==, 1.0 / 4.0, 1.0e-15, 0.0);
}

void
test_ncm_util_sinh (void)
{
  const gdouble xs[] = {0.95, 1.5, 3.0, -2.0};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (xs); i++)
  {
    const gdouble x = xs[i];

    ncm_assert_cmpdouble_e (ncm_util_sinh1 (x), ==, sinh (x) / x, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_sinh3 (x), ==, 6.0 * (sinh (x) - x) / gsl_pow_3 (x), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_sinhx_m_xcoshx_x3 (x), ==, (sinh (x) - x * cosh (x)) / gsl_pow_3 (x), 1.0e-14, 0.0);
  }

  /* Taylor series at small x */
  {
    const gdouble x  = 1.0e-3;
    const gdouble x2 = x * x;

    ncm_assert_cmpdouble_e (ncm_util_sinh1 (x), ==, 1.0 + x2 / 6.0 + x2 * x2 / 120.0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_sinh3 (x), ==, 1.0 + x2 / 20.0 + x2 * x2 / 840.0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_sinhx_m_xcoshx_x3 (x), ==, -1.0 / 3.0 - x2 / 30.0 - x2 * x2 / 840.0, 1.0e-15, 0.0);
  }

  /* Continuity across |x| = 0.9, where the series is replaced */
  {
    const gdouble xm = nextafter (0.9, 0.0);

    ncm_assert_cmpdouble_e (ncm_util_sinh1 (xm), ==, ncm_util_sinh1 (0.9), 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_sinh3 (xm), ==, ncm_util_sinh3 (0.9), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_util_sinhx_m_xcoshx_x3 (xm), ==, ncm_util_sinhx_m_xcoshx_x3 (0.9), 1.0e-14, 0.0);
  }
}

void
test_ncm_util_mln_1mIexpzA_1pIexpmzA (void)
{
  /* {rho, theta, A}: e^{|rho|}|A| below and above 0.1, where the series is replaced */
  const gdouble points[][3] = {
    {0.3, 0.2, 0.05},
    {-1.2, 0.7, 0.02},
    {0.0, 1.0, 0.0999},
    {0.0, 1.0, 0.1001},
    {0.5, -0.4, 0.3},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (points); i++)
  {
    const gdouble rho       = points[i][0];
    const gdouble theta     = points[i][1];
    const gdouble A         = points[i][2];
    const complex double z  = rho + I * theta;
    const complex double z1 = z - clog ((1.0 - I * A * cexp (z)) / (1.0 + I * A * cexp (-z)));
    gdouble rho1, theta1;

    ncm_util_mln_1mIexpzA_1pIexpmzA (rho, theta, A, &rho1, &theta1);

    ncm_assert_cmpdouble_e (rho1, ==, creal (z1), 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (theta1, ==, cimag (z1), 1.0e-14, 1.0e-15);
  }
}

void
test_ncm_util_cmp (void)
{
  g_assert_cmpint (ncm_cmp (1.0, 1.0 + 1.0e-10, 1.0e-9, 0.0), ==, 0);
  g_assert_cmpint (ncm_cmp (1.0, 1.0 + 1.0e-8, 1.0e-9, 0.0), ==, -1);
  g_assert_cmpint (ncm_cmp (1.0 + 1.0e-8, 1.0, 1.0e-9, 0.0), ==, 1);
  g_assert_cmpint (ncm_cmp (1.0, 1.0 + 1.0e-8, 1.0e-9, 1.0e-7), ==, 0);
  g_assert_cmpint (ncm_cmp (0.0, 0.0, 0.0, 0.0), ==, 0);

  /* With one of them zero, reltol acts as an absolute tolerance */
  g_assert_cmpint (ncm_cmp (0.0, 1.0e-10, 1.0e-9, 0.0), ==, 0);
  g_assert_cmpint (ncm_cmp (0.0, 1.0e-8, 1.0e-9, 0.0), ==, -1);

  g_assert_cmpfloat (ncm_cmpdbl (2.0, 2.0), ==, 0.0);
  ncm_assert_cmpdouble_e (ncm_cmpdbl (1.0, 3.0), ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_cmpdbl (3.0, 1.0), ==, 1.0, 1.0e-15, 0.0);
}

void
test_ncm_util_rational_coarse_double (void)
{
  const gdouble xs[] = {0.0, 0.1, 1.0 / 3.0, M_PI, -2.5e10, 1.0e-20, 7.0};
  mpq_t q;
  guint i;

  mpq_init (q);

  for (i = 0; i < G_N_ELEMENTS (xs); i++)
  {
    ncm_rational_coarse_double (xs[i], q);
    ncm_assert_cmpdouble_e (mpq_get_d (q), ==, xs[i], 1.0e-15, 0.0);
  }

  /* Exact small rationals are recovered */
  ncm_rational_coarse_double (1.0 / 3.0, q);
  g_assert_cmpint (mpz_cmp_ui (mpq_numref (q), 1), ==, 0);
  g_assert_cmpint (mpz_cmp_ui (mpq_denref (q), 3), ==, 0);

  mpq_clear (q);

  {
    mpz_t a, b, c;

    ncm_mpz_inits (a, b, c, NULL);
    mpz_set_ui (c, 5);
    g_assert_cmpint (mpz_cmp_ui (a, 0), ==, 0);
    ncm_mpz_clears (a, b, c, NULL);
  }
}

void
test_ncm_util_fact_size (void)
{
  gulong n;

  /* n_f >= n and its only prime factors are 2, 3, 5 and 7 */
  for (n = 2; n <= 5000; n++)
  {
    gulong nf = ncm_util_fact_size (n);

    g_assert_cmpuint (nf, >=, n);

    while (nf % 2 == 0)
      nf /= 2;

    while (nf % 3 == 0)
      nf /= 3;

    while (nf % 5 == 0)
      nf /= 5;

    while (nf % 7 == 0)
      nf /= 7;

    g_assert_cmpuint (nf, ==, 1);
  }

  g_assert_cmpuint (ncm_util_fact_size (11), ==, 12);
  g_assert_cmpuint (ncm_util_fact_size (49), ==, 49);
}

void
test_ncm_util_basename_fits (void)
{
  const gchar *cases[][2] = {
    {"catalog.fits", "catalog"},
    {"dir/catalog.FIT", "dir/catalog"},
    {"catalog.fits.gz", "catalog.fits.gz"},
    {"catalog", "catalog"},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    gchar *base = ncm_util_basename_fits (cases[i][0]);

    g_assert_cmpstr (base, ==, cases[i][1]);
    g_free (base);
  }
}

void
test_ncm_util_function_params (void)
{
  gdouble *x = NULL;
  guint len  = 0;
  gchar *name;

  name = ncm_util_function_params ("gauss(1.5, -2e-3, 7)", &x, &len);
  g_assert_cmpstr (name, ==, "gauss");
  g_assert_cmpuint (len, ==, 3);
  g_assert_cmpfloat (x[0], ==, 1.5);
  g_assert_cmpfloat (x[1], ==, -2.0e-3);
  g_assert_cmpfloat (x[2], ==, 7.0);
  g_free (name);
  g_free (x);

  name = ncm_util_function_params ("NcHICosmoDEXcdm", &x, &len);
  g_assert_cmpstr (name, ==, "NcHICosmoDEXcdm");
  g_assert_cmpuint (len, ==, 0);
  g_assert_null (x);
  g_free (name);

  name = ncm_util_function_params ("1bad", &x, &len);
  g_assert_null (name);
}

void
test_ncm_util_function_params_invalid_subprocess (void)
{
  gdouble *x = NULL;
  guint len  = 0;

  ncm_util_function_params ("f(1.0, 2.0.0)", &x, &len);
}

void
test_ncm_util_smooth_trans (void)
{
  const gdouble f0 = 2.0;
  const gdouble f1 = -3.0;
  const gdouble z0 = 1.0;
  const gdouble dz = 0.5;
  gdouble theta0, theta1;

  ncm_assert_cmpdouble_e (ncm_util_smooth_trans (f0, f1, z0, dz, z0), ==, f0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_smooth_trans (f0, f1, z0, dz, z0 + dz), ==, f1, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_smooth_trans (f0, f1, z0, dz, z0 + 0.5 * dz), ==, 0.5 * (f0 + f1), 1.0e-15, 0.0);

  ncm_util_smooth_trans_get_theta (z0, dz, z0 + 0.3 * dz, &theta0, &theta1);
  ncm_assert_cmpdouble_e (theta0 + theta1, ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (f0 * theta0 + f1 * theta1, ==, ncm_util_smooth_trans (f0, f1, z0, dz, z0 + 0.3 * dz), 1.0e-15, 0.0);
}

void
test_ncm_util_sky_geometry (void)
{
  /* Separations along the equator and to the pole */
  ncm_assert_cmpdouble_e (ncm_util_great_circle_distance (0.0, 0.0, 90.0, 0.0), ==, 90.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_great_circle_distance (10.0, 20.0, 190.0, 20.0), ==, 140.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_great_circle_distance (35.0, -10.0, 200.0, 90.0), ==, 100.0, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_great_circle_distance (35.0, -10.0, 35.0, -10.0), ==, 0.0, 0.0, 1.0e-14);

  /* East of North */
  ncm_assert_cmpdouble_e (ncm_util_position_angle (0.0, 0.0, 0.0, 10.0), ==, 0.0, 0.0, 1.0e-15);
  ncm_assert_cmpdouble_e (ncm_util_position_angle (0.0, 0.0, 10.0, 0.0), ==, M_PI_2, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_util_position_angle (0.0, 0.0, -10.0, 0.0), ==, -M_PI_2, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (fabs (ncm_util_position_angle (0.0, 10.0, 0.0, 0.0)), ==, M_PI, 1.0e-15, 0.0);
}

void
test_ncm_util_error (void)
{
  const GQuark domain = g_quark_from_static_string ("test-ncm-util-error");
  GError *error       = NULL;
  GError *local_error = NULL;

  ncm_util_set_or_call_error (&error, domain, 3, "value %d", 7);
  g_assert_error (error, domain, 3);
  g_assert_cmpstr (error->message, ==, "value 7");
  g_clear_error (&error);

  /* Nothing to forward */
  ncm_util_forward_or_call_error (&error, NULL, "context");
  g_assert_null (error);

  ncm_util_set_or_call_error (&local_error, domain, 4, "inner");
  ncm_util_forward_or_call_error (&error, local_error, "outer %s", "call");
  g_assert_error (error, domain, 4);
  g_assert_cmpstr (error->message, ==, "outer call: inner");
  g_clear_error (&error);
}

void
test_ncm_util_error_set_null_subprocess (void)
{
  ncm_util_set_or_call_error (NULL, g_quark_from_static_string ("test"), 1, "fatal");
}

void
test_ncm_util_error_pile_up_subprocess (void)
{
  const GQuark domain = g_quark_from_static_string ("test");
  GError *error       = NULL;

  ncm_util_set_or_call_error (&error, domain, 1, "first");
  ncm_util_set_or_call_error (&error, domain, 2, "second");
}

void
test_ncm_util_error_forward_null_subprocess (void)
{
  GError *local_error = NULL;

  ncm_util_set_or_call_error (&local_error, g_quark_from_static_string ("test"), 1, "inner");
  ncm_util_forward_or_call_error (NULL, local_error, "outer");
}

void
test_ncm_util_error_traps (void)
{
  g_test_trap_subprocess ("/ncm/util/error/set_null/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/util/error/pile_up/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/util/error/forward_null/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/util/function_params/invalid/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}


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

void test_ncm_util_projected_radius (void);
void test_ncm_util_complex (void);
void test_ncm_util_gaussian_int (void);
void test_ncm_util_gaussian_int_rng_two_sides (void);
void test_ncm_util_gaussian_int_rng_one_side (void);
void test_ncm_util_gaussian_int_nonunit (void);

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

    ncm_complex_res_add_mul (c1, c2, c3);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -6.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 44.0);

    ncm_complex_res_add_mul_real (c1, c2, 2.0);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, 0.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 52.0);

    ncm_complex_res_mul (c1, c2);
    g_assert_cmpfloat (ncm_complex_Re (c1), ==, -208.0);
    g_assert_cmpfloat (ncm_complex_Im (c1), ==, 156.0);

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

    z = ncm_complex_c (&c1);
    g_assert_cmpfloat (creal (z), ==, 1.0);
    g_assert_cmpfloat (cimag (z), ==, 2.0);
  }

  {
    NcmComplex c1  = NCM_COMPLEX_INIT (1.0 + 2.0 * I);
    NcmComplex *c2 = ncm_complex_dup (&c1);

    g_assert_cmpfloat (ncm_complex_Re (&c1), ==, 1.0);
    g_assert_cmpfloat (ncm_complex_Im (&c1), ==, 2.0);

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
  gdouble sign, logtol;
  gint i;

  for (i = 0; i < 10; i++)
  {
    const gdouble xl       = g_test_rand_double_range (0.0, 0.0);
    const gdouble xu       = g_test_rand_double_range (xl, 10.0);
    const gdouble symint_l = gsl_cdf_chisq_P (xl * xl, 1);
    const gdouble symint_u = gsl_cdf_chisq_P (xu * xu, 1);
    const gdouble int_val  = 0.5 * (symint_u - symint_l);

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


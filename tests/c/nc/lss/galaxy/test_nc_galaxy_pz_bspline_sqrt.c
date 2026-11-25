/***************************************************************************
 *            test_nc_galaxy_pz_bspline_sqrt.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_pz_bspline_sqrt.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>

/*
 * The table of the tests: a Gaussian of mean TEST_PZ_MU and standard deviation
 * TEST_PZ_SIGMA at TEST_PZ_N redshifts from TEST_PZ_Z0 in steps of TEST_PZ_DZ, the step of
 * the HSC tables. The last redshift is TEST_PZ_ZF.
 */
#define TEST_PZ_MU 0.8
#define TEST_PZ_SIGMA 0.1
#define TEST_PZ_Z0 0.2
#define TEST_PZ_DZ 0.01
#define TEST_PZ_N 121
#define TEST_PZ_ZF (TEST_PZ_Z0 + TEST_PZ_DZ * (TEST_PZ_N - 1))

#define TEST_PZ_NPOINTS 1000
#define TEST_PZ_NDRAWS 20000

/*
 * The polynomial case of the moments: sqrt (p) is the degree-5 polynomial
 * _test_pz_poly_coeffs, tabulated at TEST_PZ_POLY_N non-uniform redshifts
 * z_j = TEST_PZ_POLY_Z0 + TEST_PZ_POLY_DZ (j / (TEST_PZ_POLY_N - 1))^1.5, and the moments are
 * compared with a TEST_PZ_POLY_GL-point Gauss-Legendre rule on the whole interval, exact
 * up to degree 2 TEST_PZ_POLY_GL - 1 >= 10 + TEST_PZ_POLY_KMAX.
 */
#define TEST_PZ_POLY_N 21
#define TEST_PZ_POLY_Z0 0.1
#define TEST_PZ_POLY_DZ 2.0
#define TEST_PZ_POLY_GL 40
#define TEST_PZ_POLY_KMAX 32

static const gdouble _test_pz_poly_coeffs[6] = {1.0, 0.5, 0.3, -0.1, 0.02, 0.003};

typedef struct _TestNcGalaxyPzBSplineSqrt
{
  NcmVector *zv;
  NcmVector *pv;
  NcGalaxyPzBSplineSqrt *gpzs;
} TestNcGalaxyPzBSplineSqrt;

/* A subclass that implements only the range, for the call to the moments it lacks */
#define TEST_TYPE_NC_GALAXY_PZ_NO_MOMENTS (test_nc_galaxy_pz_no_moments_get_type ())

G_DECLARE_FINAL_TYPE (TestNcGalaxyPzNoMoments, test_nc_galaxy_pz_no_moments, TEST, NC_GALAXY_PZ_NO_MOMENTS, NcGalaxyPz)

struct _TestNcGalaxyPzNoMoments
{
  NcGalaxyPz parent_instance;
};

G_DEFINE_TYPE (TestNcGalaxyPzNoMoments, test_nc_galaxy_pz_no_moments, NC_TYPE_GALAXY_PZ)

static void
_test_nc_galaxy_pz_no_moments_get_range (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max)
{
  *z_min = 0.0;
  *z_max = 1.0;
}

static void
test_nc_galaxy_pz_no_moments_init (TestNcGalaxyPzNoMoments *gpz)
{
}

static void
test_nc_galaxy_pz_no_moments_class_init (TestNcGalaxyPzNoMomentsClass *klass)
{
  NC_GALAXY_PZ_CLASS (klass)->get_range = &_test_nc_galaxy_pz_no_moments_get_range;
}

void test_nc_galaxy_pz_bspline_sqrt_new (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_free (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);

void test_nc_galaxy_pz_bspline_sqrt_range (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_eval (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_table (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_norm (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_gen (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_new_bspline (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_serialize (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);
void test_nc_galaxy_pz_bspline_sqrt_moments_norm (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata);

void test_nc_galaxy_pz_bspline_sqrt_negative (void);
void test_nc_galaxy_pz_bspline_sqrt_negative_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_few_points (void);
void test_nc_galaxy_pz_bspline_sqrt_few_points_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_zero (void);
void test_nc_galaxy_pz_bspline_sqrt_zero_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_order (void);
void test_nc_galaxy_pz_bspline_sqrt_order_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_missing (void);
void test_nc_galaxy_pz_bspline_sqrt_missing_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_poly (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_outside (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_outside_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_empty (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_empty_subprocess (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_missing (void);
void test_nc_galaxy_pz_bspline_sqrt_moments_missing_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/galaxy_pz_bspline_sqrt/range", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_range,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/eval", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_eval,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/table", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_table,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/norm", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_norm,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/gen", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_gen,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/new_bspline", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_new_bspline,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/serialize", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_serialize,
              &test_nc_galaxy_pz_bspline_sqrt_free);
  g_test_add ("/nc/galaxy_pz_bspline_sqrt/moments/norm", TestNcGalaxyPzBSplineSqrt, NULL,
              &test_nc_galaxy_pz_bspline_sqrt_new,
              &test_nc_galaxy_pz_bspline_sqrt_moments_norm,
              &test_nc_galaxy_pz_bspline_sqrt_free);

  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/negative", &test_nc_galaxy_pz_bspline_sqrt_negative);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/negative/subprocess", &test_nc_galaxy_pz_bspline_sqrt_negative_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/few_points", &test_nc_galaxy_pz_bspline_sqrt_few_points);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/few_points/subprocess", &test_nc_galaxy_pz_bspline_sqrt_few_points_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/zero", &test_nc_galaxy_pz_bspline_sqrt_zero);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/zero/subprocess", &test_nc_galaxy_pz_bspline_sqrt_zero_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/order", &test_nc_galaxy_pz_bspline_sqrt_order);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/order/subprocess", &test_nc_galaxy_pz_bspline_sqrt_order_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/missing", &test_nc_galaxy_pz_bspline_sqrt_missing);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/missing/subprocess", &test_nc_galaxy_pz_bspline_sqrt_missing_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/poly", &test_nc_galaxy_pz_bspline_sqrt_moments_poly);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/outside", &test_nc_galaxy_pz_bspline_sqrt_moments_outside);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/outside/subprocess", &test_nc_galaxy_pz_bspline_sqrt_moments_outside_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/empty", &test_nc_galaxy_pz_bspline_sqrt_moments_empty);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/empty/subprocess", &test_nc_galaxy_pz_bspline_sqrt_moments_empty_subprocess);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/missing", &test_nc_galaxy_pz_bspline_sqrt_moments_missing);
  g_test_add_func ("/nc/galaxy_pz_bspline_sqrt/moments/missing/subprocess", &test_nc_galaxy_pz_bspline_sqrt_moments_missing_subprocess);

  g_test_run ();
}

/* The Gaussian of the table at z */
static gdouble
_test_pz_gauss (const gdouble z)
{
  return exp (-0.5 * gsl_pow_2 ((z - TEST_PZ_MU) / TEST_PZ_SIGMA));
}

/* The polynomial sqrt (p) of the moments case at z, by Horner's rule */
static gdouble
_test_pz_poly (const gdouble z)
{
  gdouble q = 0.0;
  gint i;

  for (i = 5; i >= 0; i--)
    q = q * z + _test_pz_poly_coeffs[i];

  return q;
}

/*
 * The moments of the square of the polynomial on [a, b] for k <= n, by the
 * TEST_PZ_POLY_GL-point Gauss-Legendre rule on [a, b]: exact, and independent of the knots.
 */
static void
_test_pz_poly_moments (const gdouble a, const gdouble b, const guint n, gdouble *m)
{
  gsl_integration_glfixed_table *glt = gsl_integration_glfixed_table_alloc (TEST_PZ_POLY_GL);
  guint i, k;

  for (k = 0; k <= n; k++)
    m[k] = 0.0;

  for (i = 0; i < TEST_PZ_POLY_GL; i++)
  {
    gdouble zq, wq;

    gsl_integration_glfixed_point (a, b, i, &zq, &wq, glt);
    {
      const gdouble t  = (2.0 * zq - (a + b)) / (b - a);
      const gdouble pw = wq * gsl_pow_2 (_test_pz_poly (zq));

      for (k = 0; k <= n; k++)
        m[k] += pw * cos (k * acos (GSL_MAX (-1.0, GSL_MIN (1.0, t))));
    }
  }

  gsl_integration_glfixed_table_free (glt);
}

/* The n redshifts from TEST_PZ_Z0 in steps of TEST_PZ_DZ, and scale times the Gaussian at them */
static void
_test_pz_table (const guint n, const gdouble scale, NcmVector **zv, NcmVector **pv)
{
  guint j;

  *zv = ncm_vector_new (n);
  *pv = ncm_vector_new (n);

  for (j = 0; j < n; j++)
  {
    const gdouble z = TEST_PZ_Z0 + TEST_PZ_DZ * j;

    ncm_vector_set (*zv, j, z);
    ncm_vector_set (*pv, j, scale * _test_pz_gauss (z));
  }
}

void
test_nc_galaxy_pz_bspline_sqrt_new (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  _test_pz_table (TEST_PZ_N, 1.0, &test->zv, &test->pv);

  test->gpzs = nc_galaxy_pz_bspline_sqrt_new_from_table (test->zv, test->pv);

  g_assert_true (NC_IS_GALAXY_PZ_BSPLINE_SQRT (test->gpzs));
  g_assert_true (NC_IS_GALAXY_PZ (test->gpzs));
}

void
test_nc_galaxy_pz_bspline_sqrt_free (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  ncm_vector_free (test->zv);
  ncm_vector_free (test->pv);

  NCM_TEST_FREE (nc_galaxy_pz_free, NC_GALAXY_PZ (test->gpzs));
}

void
test_nc_galaxy_pz_bspline_sqrt_range (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  gdouble z_min, z_max;

  nc_galaxy_pz_get_range (NC_GALAXY_PZ (test->gpzs), &z_min, &z_max);

  g_assert_cmpfloat (z_min, ==, ncm_vector_get (test->zv, 0));
  g_assert_cmpfloat (z_max, ==, ncm_vector_get (test->zv, TEST_PZ_N - 1));
}

/*
 * The density is the square of the order-6 NcmSplineBSpline through (z_j, sqrt (p_j)),
 * and its logarithm twice the logarithm of the absolute value of that spline: equal
 * values, since both run the same evaluation on the same knots and control points.
 */
void
test_nc_galaxy_pz_bspline_sqrt_eval (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  NcGalaxyPz *gpz = NC_GALAXY_PZ (test->gpzs);
  NcmVector *sv   = ncm_vector_new (TEST_PZ_N);
  NcmRNG *rng     = ncm_rng_seeded_new (NULL, 20261117);
  NcmSpline *s;
  guint i;

  for (i = 0; i < TEST_PZ_N; i++)
    ncm_vector_set (sv, i, sqrt (ncm_vector_get (test->pv, i)));

  s = NCM_SPLINE (ncm_spline_bspline_new_full (NC_GALAXY_PZ_BSPLINE_SQRT_ORDER, test->zv, sv, TRUE));

  for (i = 0; i < TEST_PZ_NPOINTS; i++)
  {
    const gdouble z   = ncm_rng_uniform_gen (rng, TEST_PZ_Z0, TEST_PZ_ZF);
    const gdouble s_z = ncm_spline_eval (s, z);

    g_assert_cmpfloat (nc_galaxy_pz_eval_p (gpz, z), ==, s_z * s_z);
    g_assert_cmpfloat (nc_galaxy_pz_eval_ln_p (gpz, z), ==, 2.0 * log (fabs (s_z)));
  }

  ncm_rng_free (rng);
  ncm_spline_free (s);
  ncm_vector_free (sv);
}

/*
 * The spline interpolates sqrt (p_j), so the density at z_j is p_j up to the rounding of
 * the banded solve and of the square. The largest relative difference over the table
 * measured 8.9e-16; 1e-14 is asserted.
 */
void
test_nc_galaxy_pz_bspline_sqrt_table (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  NcGalaxyPz *gpz = NC_GALAXY_PZ (test->gpzs);
  guint j;

  for (j = 0; j < TEST_PZ_N; j++)
    ncm_assert_cmpdouble_e (nc_galaxy_pz_eval_p (gpz, ncm_vector_get (test->zv, j)), ==, ncm_vector_get (test->pv, j), 1.0e-14, 0.0);
}

/*
 * The normalization against the integral of the Gaussian over the range of the table.
 * The difference is the interpolation error of the spline of sqrt (p) with ten points per
 * standard deviation: it measured 1.6e-11 relative; 1e-10 is asserted.
 */
void
test_nc_galaxy_pz_bspline_sqrt_norm (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  const gdouble a     = (TEST_PZ_Z0 - TEST_PZ_MU) / (M_SQRT2 * TEST_PZ_SIGMA);
  const gdouble b     = (TEST_PZ_ZF - TEST_PZ_MU) / (M_SQRT2 * TEST_PZ_SIGMA);
  const gdouble exact = TEST_PZ_SIGMA * sqrt (0.5 * M_PI) * (gsl_sf_erf (b) - gsl_sf_erf (a));

  ncm_assert_cmpdouble_e (nc_galaxy_pz_get_norm (NC_GALAXY_PZ (test->gpzs)), ==, exact, 1.0e-10, 0.0);
}

/*
 * Draws from the density: every draw is in the range, the same seed gives the same draws,
 * and the mean and the variance of TEST_PZ_NDRAWS draws agree with those of the Gaussian
 * within five standard errors, sigma / sqrt (n) and sigma^2 sqrt (2 / n). The table cuts
 * the Gaussian at six standard deviations, which leaves the mean and lowers the variance
 * by 7.3e-8 of its value, against a tolerance of 5e-2 of it.
 */
void
test_nc_galaxy_pz_bspline_sqrt_gen (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  NcGalaxyPz *gpz        = NC_GALAXY_PZ (test->gpzs);
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, 1234);
  NcmRNG *rng2           = ncm_rng_seeded_new (NULL, 1234);
  const gdouble mean_tol = 5.0 * TEST_PZ_SIGMA / sqrt (TEST_PZ_NDRAWS);
  const gdouble var_tol  = 5.0 * TEST_PZ_SIGMA * TEST_PZ_SIGMA * sqrt (2.0 / TEST_PZ_NDRAWS);
  gdouble sum            = 0.0;
  gdouble sum2           = 0.0;
  guint i;

  for (i = 0; i < TEST_PZ_NDRAWS; i++)
  {
    const gdouble z = nc_galaxy_pz_gen (gpz, rng);

    g_assert_cmpfloat (z, >=, TEST_PZ_Z0);
    g_assert_cmpfloat (z, <=, TEST_PZ_ZF);
    g_assert_cmpfloat (z, ==, nc_galaxy_pz_gen (gpz, rng2));

    sum  += z - TEST_PZ_MU;
    sum2 += gsl_pow_2 (z - TEST_PZ_MU);
  }

  {
    const gdouble mean = sum / TEST_PZ_NDRAWS;
    const gdouble var  = sum2 / TEST_PZ_NDRAWS - mean * mean;

    ncm_assert_cmpdouble_e (mean, ==, 0.0, 0.0, mean_tol);
    ncm_assert_cmpdouble_e (var, ==, TEST_PZ_SIGMA * TEST_PZ_SIGMA, 0.0, var_tol);
  }

  ncm_rng_free (rng);
  ncm_rng_free (rng2);
}

/* An object built from the B-spline of another references it and gives the same density */
void
test_nc_galaxy_pz_bspline_sqrt_new_bspline (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  NcmBSpline *s               = nc_galaxy_pz_bspline_sqrt_peek_bspline (test->gpzs);
  NcGalaxyPzBSplineSqrt *gpzs = nc_galaxy_pz_bspline_sqrt_new (s);
  guint i;

  g_assert_true (NCM_IS_BSPLINE (s));
  g_assert_cmpuint (ncm_bspline_get_order (s), ==, NC_GALAXY_PZ_BSPLINE_SQRT_ORDER);
  g_assert_true (nc_galaxy_pz_bspline_sqrt_peek_bspline (gpzs) == s);

  g_assert_cmpfloat (nc_galaxy_pz_get_norm (NC_GALAXY_PZ (gpzs)), ==, nc_galaxy_pz_get_norm (NC_GALAXY_PZ (test->gpzs)));

  for (i = 0; i < 101; i++)
  {
    const gdouble z = TEST_PZ_Z0 + (TEST_PZ_ZF - TEST_PZ_Z0) * i / 100.0;

    g_assert_cmpfloat (nc_galaxy_pz_eval_p (NC_GALAXY_PZ (gpzs), z), ==, nc_galaxy_pz_eval_p (NC_GALAXY_PZ (test->gpzs), z));
  }

  nc_galaxy_pz_free (NC_GALAXY_PZ (gpzs));
}

void
test_nc_galaxy_pz_bspline_sqrt_serialize (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyPz *gpz   = NC_GALAXY_PZ (test->gpzs);
  NcGalaxyPz *dup   = NC_GALAXY_PZ (ncm_serialize_dup_obj (ser, G_OBJECT (gpz)));
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 4321);
  NcmRNG *rng2      = ncm_rng_seeded_new (NULL, 4321);
  gdouble z_min, z_max, dup_z_min, dup_z_max;
  guint i;

  g_assert_true (NC_IS_GALAXY_PZ_BSPLINE_SQRT (dup));
  g_assert_true (dup != gpz);

  nc_galaxy_pz_get_range (gpz, &z_min, &z_max);
  nc_galaxy_pz_get_range (dup, &dup_z_min, &dup_z_max);

  g_assert_cmpfloat (dup_z_min, ==, z_min);
  g_assert_cmpfloat (dup_z_max, ==, z_max);
  g_assert_cmpfloat (nc_galaxy_pz_get_norm (dup), ==, nc_galaxy_pz_get_norm (gpz));

  for (i = 0; i < 101; i++)
  {
    const gdouble z = z_min + (z_max - z_min) * i / 100.0;

    g_assert_cmpfloat (nc_galaxy_pz_eval_p (dup, z), ==, nc_galaxy_pz_eval_p (gpz, z));
    g_assert_cmpfloat (nc_galaxy_pz_eval_ln_p (dup, z), ==, nc_galaxy_pz_eval_ln_p (gpz, z));
    g_assert_cmpfloat (nc_galaxy_pz_gen (dup, rng), ==, nc_galaxy_pz_gen (gpz, rng2));
  }

  ncm_rng_free (rng);
  ncm_rng_free (rng2);
  nc_galaxy_pz_free (dup);
  ncm_serialize_free (ser);
}

void
test_nc_galaxy_pz_bspline_sqrt_negative (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/negative/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the p(z) table value 7 is -0.5 at z = 0.27; values must be non-negative*");
}

void
test_nc_galaxy_pz_bspline_sqrt_negative_subprocess (void)
{
  NcmVector *zv, *pv;

  _test_pz_table (TEST_PZ_N, 1.0, &zv, &pv);
  ncm_vector_set (pv, 7, -0.5);

  nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv);
}

void
test_nc_galaxy_pz_bspline_sqrt_few_points (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/few_points/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the p(z) table has 5 points, at least 6 are required*");
}

void
test_nc_galaxy_pz_bspline_sqrt_few_points_subprocess (void)
{
  NcmVector *zv, *pv;

  _test_pz_table (5, 1.0, &zv, &pv);

  nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv);
}

void
test_nc_galaxy_pz_bspline_sqrt_zero (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcGalaxyPzBSplineSqrt: the normalization is 0; it must be positive*");
}

void
test_nc_galaxy_pz_bspline_sqrt_zero_subprocess (void)
{
  NcmVector *zv, *pv;

  _test_pz_table (TEST_PZ_N, 0.0, &zv, &pv);

  nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv);
}

void
test_nc_galaxy_pz_bspline_sqrt_order (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcGalaxyPzBSplineSqrt: the B-spline must have order 6, got 4*");
}

void
test_nc_galaxy_pz_bspline_sqrt_order_subprocess (void)
{
  NcmVector *zv, *pv;
  NcmSplineBSpline *sbs;

  _test_pz_table (TEST_PZ_N, 1.0, &zv, &pv);
  sbs = ncm_spline_bspline_new_full (4, zv, pv, TRUE);

  nc_galaxy_pz_bspline_sqrt_new (ncm_spline_bspline_dup_bspline (sbs));
}

void
test_nc_galaxy_pz_bspline_sqrt_missing (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/missing/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcGalaxyPzBSplineSqrt: the B-spline is required*");
}

void
test_nc_galaxy_pz_bspline_sqrt_missing_subprocess (void)
{
  nc_galaxy_pz_bspline_sqrt_new (NULL);
}

/*
 * The mass m_0 of the whole range is the normalization, and masses add over two
 * intervals: both are exact integrals of the same s^2, equal to rounding. Measured 2.7e-15
 * and 2.4e-15 relative; 1e-14 is asserted.
 */
void
test_nc_galaxy_pz_bspline_sqrt_moments_norm (TestNcGalaxyPzBSplineSqrt *test, gconstpointer pdata)
{
  NcGalaxyPz *gpz = NC_GALAXY_PZ (test->gpzs);
  const gdouble c = 0.7543;
  GArray *m_all   = NULL;
  GArray *m_lo    = NULL;
  GArray *m_hi    = NULL;

  nc_galaxy_pz_cheb_moments (gpz, TEST_PZ_Z0, TEST_PZ_ZF, 0, &m_all);
  nc_galaxy_pz_cheb_moments (gpz, TEST_PZ_Z0, c, 0, &m_lo);
  nc_galaxy_pz_cheb_moments (gpz, c, TEST_PZ_ZF, 0, &m_hi);

  g_assert_cmpuint (m_all->len, ==, 1);
  ncm_assert_cmpdouble_e (g_array_index (m_all, gdouble, 0), ==, nc_galaxy_pz_get_norm (gpz), 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (g_array_index (m_lo, gdouble, 0) + g_array_index (m_hi, gdouble, 0), ==,
                          g_array_index (m_all, gdouble, 0), 1.0e-14, 0.0);

  g_array_unref (m_all);
  g_array_unref (m_lo);
  g_array_unref (m_hi);
}

/*
 * sqrt (p) a degree-5 polynomial on non-uniform knots: the B-spline reproduces it, and the
 * moments up to TEST_PZ_POLY_KMAX equal those of the Gauss-Legendre rule on the whole
 * interval, for the whole range, an interval inside one piece, and an interval across
 * pieces with ends between knots. The difference is rounding: measured 5.8e-15, 8.9e-14 and
 * 4.8e-15 of m_0 (the interval inside one piece has the least mass); 1e-12 of m_0 is
 * asserted. A reused array is resized to n + 1, also to a single value.
 */
void
test_nc_galaxy_pz_bspline_sqrt_moments_poly (void)
{
  const gdouble z_hi            = TEST_PZ_POLY_Z0 + TEST_PZ_POLY_DZ;
  const gdouble intervals[3][2] =
  {
    {TEST_PZ_POLY_Z0, z_hi},
    {TEST_PZ_POLY_Z0 + 0.4 * 0.01, TEST_PZ_POLY_Z0 + 0.6 * 0.01},
    {0.3172, 1.8341},
  };
  NcmVector *zv = ncm_vector_new (TEST_PZ_POLY_N);
  NcmVector *pv = ncm_vector_new (TEST_PZ_POLY_N);
  GArray *m     = NULL;
  gdouble m_ref[TEST_PZ_POLY_KMAX + 1];
  NcGalaxyPzBSplineSqrt *gpzs;
  guint j, r, k;

  for (j = 0; j < TEST_PZ_POLY_N; j++)
  {
    const gdouble z = TEST_PZ_POLY_Z0 + TEST_PZ_POLY_DZ * pow ((gdouble) j / (TEST_PZ_POLY_N - 1), 1.5);

    g_assert_cmpfloat (_test_pz_poly (z), >, 0.0);
    ncm_vector_set (zv, j, z);
    ncm_vector_set (pv, j, gsl_pow_2 (_test_pz_poly (z)));
  }

  gpzs = nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv);

  for (r = 0; r < 3; r++)
  {
    const gdouble a = intervals[r][0];
    const gdouble b = intervals[r][1];

    nc_galaxy_pz_cheb_moments (NC_GALAXY_PZ (gpzs), a, b, TEST_PZ_POLY_KMAX, &m);
    _test_pz_poly_moments (a, b, TEST_PZ_POLY_KMAX, m_ref);

    g_assert_cmpuint (m->len, ==, TEST_PZ_POLY_KMAX + 1);

    for (k = 0; k <= TEST_PZ_POLY_KMAX; k++)
      ncm_assert_cmpdouble_e (g_array_index (m, gdouble, k), ==, m_ref[k], 0.0, 1.0e-12 * m_ref[0]);
  }

  nc_galaxy_pz_cheb_moments (NC_GALAXY_PZ (gpzs), intervals[2][0], intervals[2][1], 0, &m);
  g_assert_cmpuint (m->len, ==, 1);
  ncm_assert_cmpdouble_e (g_array_index (m, gdouble, 0), ==, m_ref[0], 0.0, 1.0e-12 * m_ref[0]);

  g_array_unref (m);
  nc_galaxy_pz_free (NC_GALAXY_PZ (gpzs));
  ncm_vector_free (zv);
  ncm_vector_free (pv);
}

void
test_nc_galaxy_pz_bspline_sqrt_moments_outside (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/moments/outside/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the interval [0.1, 1] is empty or not inside the range [0.2, 1.4]*");
}

void
test_nc_galaxy_pz_bspline_sqrt_moments_outside_subprocess (void)
{
  NcmVector *zv, *pv;
  GArray *m = NULL;

  _test_pz_table (TEST_PZ_N, 1.0, &zv, &pv);
  nc_galaxy_pz_cheb_moments (NC_GALAXY_PZ (nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv)), 0.1, 1.0, 2, &m);
}

void
test_nc_galaxy_pz_bspline_sqrt_moments_empty (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/moments/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the interval [0.5, 0.5] is empty or not inside the range [0.2, 1.4]*");
}

void
test_nc_galaxy_pz_bspline_sqrt_moments_empty_subprocess (void)
{
  NcmVector *zv, *pv;
  GArray *m = NULL;

  _test_pz_table (TEST_PZ_N, 1.0, &zv, &pv);
  nc_galaxy_pz_cheb_moments (NC_GALAXY_PZ (nc_galaxy_pz_bspline_sqrt_new_from_table (zv, pv)), 0.5, 0.5, 2, &m);
}

void
test_nc_galaxy_pz_bspline_sqrt_moments_missing (void)
{
  g_test_trap_subprocess ("/nc/galaxy_pz_bspline_sqrt/moments/missing/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*TestNcGalaxyPzNoMoments does not implement the Chebyshev moments*");
}

void
test_nc_galaxy_pz_bspline_sqrt_moments_missing_subprocess (void)
{
  GArray *m = NULL;

  nc_galaxy_pz_cheb_moments (NC_GALAXY_PZ (g_object_new (TEST_TYPE_NC_GALAXY_PZ_NO_MOMENTS, NULL)), 0.2, 0.4, 2, &m);
}


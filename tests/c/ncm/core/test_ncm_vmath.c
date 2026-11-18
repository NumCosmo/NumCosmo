/***************************************************************************
 *            test_ncm_vmath.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_vmath.c
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
#include "ncm/core/ncm_vmath_private.h"

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/*
 * Tests of ncm_vmath_private.h: the six functions it declares, called in omp simd loops
 * of this file, against the same functions called one value at a time, which GCC does not
 * vectorize. Where the header declares nothing both are libm and the differences are zero.
 */

/* Values of each loop: not a multiple of a vector length, so the loops have a remainder */
#define TEST_NCM_VMATH_N 4099

/*
 * Largest difference accepted between a value of a loop and the value of libm, in
 * spacings of doubles at the libm value. It is the largest difference measured in
 * test_ncm_vmath_loop: 3 for sin and cos, 2 for exp and atan2, 1 for log and pow, the same
 * with the SSE2, AVX2 and AVX-512 variants of glibc 2.44.
 */
#define TEST_NCM_VMATH_TOL 3.0

typedef enum _TestNcmVMathFunc
{
  TEST_NCM_VMATH_EXP = 0,
  TEST_NCM_VMATH_LOG,
  TEST_NCM_VMATH_POW,
  TEST_NCM_VMATH_SIN,
  TEST_NCM_VMATH_COS,
  TEST_NCM_VMATH_ATAN2,
  TEST_NCM_VMATH_LEN,
} TestNcmVMathFunc;

typedef struct _TestNcmVMath
{
  gdouble *x;
  gdouble *y;
  gdouble *r;
} TestNcmVMath;

void test_ncm_vmath_new (TestNcmVMath *test, gconstpointer pdata);
void test_ncm_vmath_free (TestNcmVMath *test, gconstpointer pdata);
void test_ncm_vmath_loop (TestNcmVMath *test, gconstpointer pdata);
void test_ncm_vmath_special (TestNcmVMath *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/vmath/exp", TestNcmVMath, GINT_TO_POINTER (TEST_NCM_VMATH_EXP),
              &test_ncm_vmath_new, &test_ncm_vmath_loop, &test_ncm_vmath_free);
  g_test_add ("/ncm/vmath/log", TestNcmVMath, GINT_TO_POINTER (TEST_NCM_VMATH_LOG),
              &test_ncm_vmath_new, &test_ncm_vmath_loop, &test_ncm_vmath_free);
  g_test_add ("/ncm/vmath/pow", TestNcmVMath, GINT_TO_POINTER (TEST_NCM_VMATH_POW),
              &test_ncm_vmath_new, &test_ncm_vmath_loop, &test_ncm_vmath_free);
  g_test_add ("/ncm/vmath/sin", TestNcmVMath, GINT_TO_POINTER (TEST_NCM_VMATH_SIN),
              &test_ncm_vmath_new, &test_ncm_vmath_loop, &test_ncm_vmath_free);
  g_test_add ("/ncm/vmath/cos", TestNcmVMath, GINT_TO_POINTER (TEST_NCM_VMATH_COS),
              &test_ncm_vmath_new, &test_ncm_vmath_loop, &test_ncm_vmath_free);
  g_test_add ("/ncm/vmath/atan2", TestNcmVMath, GINT_TO_POINTER (TEST_NCM_VMATH_ATAN2),
              &test_ncm_vmath_new, &test_ncm_vmath_loop, &test_ncm_vmath_free);
  g_test_add ("/ncm/vmath/special", TestNcmVMath, NULL,
              &test_ncm_vmath_new, &test_ncm_vmath_special, &test_ncm_vmath_free);

  g_test_run ();
}

void
test_ncm_vmath_new (TestNcmVMath *test, gconstpointer pdata)
{
  test->x = g_new (gdouble, TEST_NCM_VMATH_N);
  test->y = g_new (gdouble, TEST_NCM_VMATH_N);
  test->r = g_new (gdouble, TEST_NCM_VMATH_N);
}

void
test_ncm_vmath_free (TestNcmVMath *test, gconstpointer pdata)
{
  g_free (test->x);
  g_free (test->y);
  g_free (test->r);
}

/* The six functions in loops that GCC vectorizes when the header declares their variants */
static void
_test_ncm_vmath_fill (TestNcmVMathFunc func, const gdouble *x, const gdouble *y, gdouble *r, const guint n)
{
  guint i;

  switch (func)
  {
    case TEST_NCM_VMATH_EXP:
      #pragma omp simd

      for (i = 0; i < n; i++)
        r[i] = exp (x[i]);

      break;
    case TEST_NCM_VMATH_LOG:
      #pragma omp simd

      for (i = 0; i < n; i++)
        r[i] = log (x[i]);

      break;
    case TEST_NCM_VMATH_POW:
      #pragma omp simd

      for (i = 0; i < n; i++)
        r[i] = pow (x[i], y[i]);

      break;
    case TEST_NCM_VMATH_SIN:
      #pragma omp simd

      for (i = 0; i < n; i++)
        r[i] = sin (x[i]);

      break;
    case TEST_NCM_VMATH_COS:
      #pragma omp simd

      for (i = 0; i < n; i++)
        r[i] = cos (x[i]);

      break;
    case TEST_NCM_VMATH_ATAN2:
      #pragma omp simd

      for (i = 0; i < n; i++)
        r[i] = atan2 (y[i], x[i]);

      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }
}

/* One value by libm: a call outside a loop is never replaced by a vector variant */
static G_GNUC_NO_INLINE gdouble
_test_ncm_vmath_one (TestNcmVMathFunc func, const gdouble x, const gdouble y)
{
  switch (func)
  {
    case TEST_NCM_VMATH_EXP:

      return exp (x);

    case TEST_NCM_VMATH_LOG:

      return log (x);

    case TEST_NCM_VMATH_POW:

      return pow (x, y);

    case TEST_NCM_VMATH_SIN:

      return sin (x);

    case TEST_NCM_VMATH_COS:

      return cos (x);

    case TEST_NCM_VMATH_ATAN2:

      return atan2 (y, x);

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/*
 * A function over the range where the library calls it in loops: exp on [-700, 700], log
 * and the base of pow on [1e-300, 1e300] (uniform in the logarithm), the exponent of pow
 * on [-3, 3] with the result kept finite, sin and cos on [-100, 100], atan2 on the square
 * [-1, 1]^2.
 */
void
test_ncm_vmath_loop (TestNcmVMath *test, gconstpointer pdata)
{
  const TestNcmVMathFunc func = GPOINTER_TO_INT (pdata);
  NcmRNG *rng                 = ncm_rng_seeded_new (NULL, 1234 + func);
  gdouble max_diff            = 0.0;
  guint i;

  for (i = 0; i < TEST_NCM_VMATH_N; i++)
  {
    switch (func)
    {
      case TEST_NCM_VMATH_EXP:
        test->x[i] = ncm_rng_uniform_gen (rng, -700.0, 700.0);
        test->y[i] = 0.0;
        break;
      case TEST_NCM_VMATH_LOG:
        test->x[i] = exp (ncm_rng_uniform_gen (rng, -690.0, 690.0));
        test->y[i] = 0.0;
        break;
      case TEST_NCM_VMATH_POW:
        test->x[i] = exp (ncm_rng_uniform_gen (rng, -230.0, 230.0));
        test->y[i] = ncm_rng_uniform_gen (rng, -3.0, 3.0);
        break;
      case TEST_NCM_VMATH_SIN:
      case TEST_NCM_VMATH_COS:
        test->x[i] = ncm_rng_uniform_gen (rng, -100.0, 100.0);
        test->y[i] = 0.0;
        break;
      case TEST_NCM_VMATH_ATAN2:
        test->x[i] = ncm_rng_uniform_gen (rng, -1.0, 1.0);
        test->y[i] = ncm_rng_uniform_gen (rng, -1.0, 1.0);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }
  }

  _test_ncm_vmath_fill (func, test->x, test->y, test->r, TEST_NCM_VMATH_N);

  for (i = 0; i < TEST_NCM_VMATH_N; i++)
  {
    const gdouble ref     = _test_ncm_vmath_one (func, test->x[i], test->y[i]);
    const gdouble spacing = nextafter (fabs (ref), GSL_POSINF) - fabs (ref);
    const gdouble diff    = fabs (test->r[i] - ref) / spacing;

    g_assert_true (gsl_finite (ref));
    max_diff = MAX (max_diff, diff);
  }

  g_test_message ("largest difference: %g spacings", max_diff);
  g_assert_cmpfloat (max_diff, <=, TEST_NCM_VMATH_TOL);

  ncm_rng_free (rng);
}

/*
 * The values at which libm returns zero, infinity or NaN, placed in a loop long enough
 * to be vectorized: the loop returns what libm returns for each of them.
 */
void
test_ncm_vmath_special (TestNcmVMath *test, gconstpointer pdata)
{
  const gdouble sx[] = {0.0, -0.0, 1.0, -1.0, GSL_POSINF, GSL_NEGINF, GSL_NAN, 710.0, -750.0, 1.0e-320, 1.0e308, -1.0e-320};
  const guint n_s    = G_N_ELEMENTS (sx);
  const guint n      = 64;
  TestNcmVMathFunc func;
  guint i, j;

  for (func = 0; func < TEST_NCM_VMATH_LEN; func++)
  {
    for (j = 0; j < n_s; j++)
    {
      for (i = 0; i < n; i++)
      {
        test->x[i] = sx[(i + j) % n_s];
        test->y[i] = sx[(i * 5 + j * 7 + 3) % n_s];
      }

      _test_ncm_vmath_fill (func, test->x, test->y, test->r, n);

      for (i = 0; i < n; i++)
      {
        const gdouble ref = _test_ncm_vmath_one (func, test->x[i], test->y[i]);

        if (gsl_isnan (ref))
        {
          g_assert_true (gsl_isnan (test->r[i]));
        }
        else if (gsl_isinf (ref) || (ref == 0.0))
        {
          g_assert_cmpfloat (test->r[i], ==, ref);
        }
        else
        {
          const gdouble spacing = nextafter (fabs (ref), GSL_POSINF) - fabs (ref);

          g_assert_cmpfloat (fabs (test->r[i] - ref) / spacing, <=, TEST_NCM_VMATH_TOL);
        }
      }
    }
  }
}


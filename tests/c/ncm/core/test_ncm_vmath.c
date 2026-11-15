/***************************************************************************
 *            test_ncm_vmath.c
 *
 *  Sat November 14 2026
 *  Copyright  2026  Caio Lima de Oliveira
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
#include "ncm/core/ncm_vmath_private.h"

#include <math.h>
#include <float.h>
#include <glib.h>

#define TEST_VMATH_N (200000)

/*
 * Built once more per x86-64 level (tests/c/meson.build), with NCM_VMATH_TEST_LEVEL 1, 3
 * or 4: the build must have generated the code of that level, and the tests are skipped
 * on a CPU without its instructions.
 */
#if defined (NCM_VMATH_TEST_LEVEL)
#  if (NCM_VMATH_TEST_LEVEL == 1) && defined (__AVX2__)
#    error "the x86-64 variant of test_ncm_vmath was built with AVX2"
#  elif (NCM_VMATH_TEST_LEVEL == 3) && (!defined (__AVX2__) || defined (__AVX512F__))
#    error "the x86-64-v3 variant of test_ncm_vmath was not built for AVX2 alone"
#  elif (NCM_VMATH_TEST_LEVEL == 4) && !defined (__AVX512F__)
#    error "the x86-64-v4 variant of test_ncm_vmath was not built with AVX-512"
#  endif
#  if (NCM_VMATH_TEST_LEVEL > 1) && !defined (NCM_VMATH_HAVE_FAST_EXP)
#    error "NCM_VMATH_HAVE_FAST_EXP is not defined with AVX2"
#  endif
#endif

/*
 * The CPU check and main() run before any test, so a variant built for an x86-64 level must
 * not execute an instruction of that level there: they are compiled for the baseline. The
 * tests themselves run only after the check passed.
 */
#if defined (NCM_VMATH_TEST_LEVEL) && (NCM_VMATH_TEST_LEVEL > 1)
#  define TEST_VMATH_BASELINE __attribute__ ((target ("arch=x86-64")))
#else
#  define TEST_VMATH_BASELINE
#endif

/*
 * A loop that vectorizes and a scalar call agree within a few ulp, not bit for bit. With
 * NCM_VMATH_HAVE_FAST_EXP the functions are the code of ncm_vmath_private.h, which the
 * compiler may contract (a * b + c into an fma) differently in vector and in scalar code
 * (GCC of Ubuntu 24.04 does): 2 ulp. Without it ncm_vmath_exp(), ncm_vmath_log() and
 * ncm_vmath_pow_pos() are those of libm, and the compiler may call the vector variants of
 * glibc (libmvec, within 4 ulp) in a loop it vectorizes: 4 ulp.
 */
#define TEST_VMATH_LOOP_ULP (2)
#define TEST_VMATH_LIBMVEC_ULP (4)

#ifdef NCM_VMATH_HAVE_FAST_EXP
#  define TEST_VMATH_ASSERT_LOOP(loop, scalar) g_assert_cmpuint (_ulp_diff ((loop), (scalar)), <=, TEST_VMATH_LOOP_ULP)
#else
#  define TEST_VMATH_ASSERT_LOOP(loop, scalar) g_assert_cmpuint (_ulp_diff ((loop), (scalar)), <=, TEST_VMATH_LIBMVEC_ULP)
#endif

static guint64 _ulp_diff (const gdouble a, const gdouble b);

TEST_VMATH_BASELINE static const gchar *
_test_ncm_vmath_missing_isa (void)
{
#if defined (NCM_VMATH_TEST_LEVEL) && (NCM_VMATH_TEST_LEVEL == 3)

  if (!__builtin_cpu_supports ("avx2") || !__builtin_cpu_supports ("fma"))
    return "the CPU does not have AVX2 and FMA";

#elif defined (NCM_VMATH_TEST_LEVEL) && (NCM_VMATH_TEST_LEVEL == 4)

  if (!__builtin_cpu_supports ("avx512f") || !__builtin_cpu_supports ("avx512vl") ||
      !__builtin_cpu_supports ("avx512bw") || !__builtin_cpu_supports ("avx512dq"))
    return "the CPU does not have AVX-512 (F, VL, BW and DQ)";

#endif

  return NULL;
}

TEST_VMATH_BASELINE static void
_test_ncm_vmath_skipped (gconstpointer missing)
{
  g_test_skip (missing);
}

static void test_ncm_vmath_select (void);
static void test_ncm_vmath_atan2_random (void);
static void test_ncm_vmath_atan2_special (void);
static void test_ncm_vmath_exp_random (void);
static void test_ncm_vmath_exp_special (void);
static void test_ncm_vmath_log_random (void);
static void test_ncm_vmath_log_special (void);
static void test_ncm_vmath_pow_pos_random (void);
static void test_ncm_vmath_sincos_random (void);
static void test_ncm_vmath_sincos_special (void);

TEST_VMATH_BASELINE gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  {
    const gchar *paths[]     = {"/ncm/vmath/select", "/ncm/vmath/atan2/random", "/ncm/vmath/atan2/special", "/ncm/vmath/exp/random", "/ncm/vmath/exp/special", "/ncm/vmath/log/random", "/ncm/vmath/log/special", "/ncm/vmath/pow_pos/random", "/ncm/vmath/sincos/random", "/ncm/vmath/sincos/special"};
    GTestFunc funcs[]        = {&test_ncm_vmath_select, &test_ncm_vmath_atan2_random, &test_ncm_vmath_atan2_special, &test_ncm_vmath_exp_random, &test_ncm_vmath_exp_special, &test_ncm_vmath_log_random, &test_ncm_vmath_log_special, &test_ncm_vmath_pow_pos_random, &test_ncm_vmath_sincos_random, &test_ncm_vmath_sincos_special};
    const gchar *missing_isa = _test_ncm_vmath_missing_isa ();
    guint i;

    for (i = 0; i < G_N_ELEMENTS (paths); i++)
    {
      if (missing_isa == NULL)
        g_test_add_func (paths[i], funcs[i]);
      else
        g_test_add_data_func (paths[i], missing_isa, &_test_ncm_vmath_skipped);
    }
  }

  g_test_run ();
}

/* Distance in ulp between two doubles of the same sign, 0 when both are NaN */
static guint64
_ulp_diff (const gdouble a, const gdouble b)
{
  gint64 ia, ib;

  if (isnan (a) && isnan (b))
    return 0;

  g_assert_true (!isnan (a) && !isnan (b));

  if (a == b)
    return 0;

  g_assert_cmpint (signbit (a) != 0, ==, signbit (b) != 0);

  memcpy (&ia, &a, sizeof (ia));
  memcpy (&ib, &b, sizeof (ib));

  return (ia > ib) ? (guint64) (ia - ib) : (guint64) (ib - ia);
}

/* Through loops that vectorize, as the functions are meant to be used */
static void
_atan2_array (const guint n, const gdouble *y, const gdouble *x, gdouble *out)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
    out[i] = ncm_vmath_atan2 (y[i], x[i]);
}

static void
_exp_array (const guint n, const gdouble *x, gdouble *out)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
    out[i] = ncm_vmath_exp (x[i]);
}

static void
_log_array (const guint n, const gdouble *x, gdouble *out)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
    out[i] = ncm_vmath_log (x[i]);
}

static void
_sincos_array (const guint n, const gdouble *x, gdouble *s, gdouble *c)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
    ncm_vmath_sincos (x[i], &s[i], &c[i]);
}

static void
_pow_pos_array (const guint n, const gdouble *x, const gdouble *y, gdouble *out)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
    out[i] = ncm_vmath_pow_pos (x[i], y[i]);
}

static void
test_ncm_vmath_select (void)
{
  const gdouble a = -0.0;
  const gdouble b = 3.5;
  gdouble r;

  g_assert_cmpfloat (ncm_vmath_select (TRUE, 1.0, 2.0), ==, 1.0);
  g_assert_cmpfloat (ncm_vmath_select (FALSE, 1.0, 2.0), ==, 2.0);
  g_assert_cmpfloat (ncm_vmath_select (7, 1.0, 2.0), ==, 1.0);

  /* The bits of the value chosen are kept */
  r = ncm_vmath_select (TRUE, a, b);
  g_assert_cmpfloat (r, ==, 0.0);
  g_assert_true (signbit (r));
  g_assert_true (isnan (ncm_vmath_select (TRUE, NAN, b)));
  g_assert_cmpfloat (ncm_vmath_select (FALSE, NAN, b), ==, b);
}

static void
test_ncm_vmath_atan2_random (void)
{
  NcmRNG *rng     = ncm_rng_seeded_new (NULL, 11);
  gdouble *y      = g_new (gdouble, TEST_VMATH_N);
  gdouble *x      = g_new (gdouble, TEST_VMATH_N);
  gdouble *out    = g_new (gdouble, TEST_VMATH_N);
  guint64 max_ulp = 0;
  guint i;

  for (i = 0; i < TEST_VMATH_N; i++)
  {
    if (i % 2 == 0)
    {
      /* Every quadrant, near the folds of the approximation */
      y[i] = ncm_rng_uniform_gen (rng, -2.0, 2.0);
      x[i] = ncm_rng_uniform_gen (rng, -2.0, 2.0);
    }
    else
    {
      /* Every scale and ratio */
      y[i] = ((i % 4 == 1) ? 1.0 : -1.0) * pow (10.0, ncm_rng_uniform_gen (rng, -300.0, 300.0));
      x[i] = ((i % 8 < 4) ? 1.0 : -1.0) * pow (10.0, ncm_rng_uniform_gen (rng, -300.0, 300.0));
    }
  }

  _atan2_array (TEST_VMATH_N, y, x, out);

  for (i = 0; i < TEST_VMATH_N; i++)
  {
    max_ulp = MAX (max_ulp, _ulp_diff (out[i], atan2 (y[i], x[i])));
    g_assert_cmpuint (_ulp_diff (out[i], ncm_vmath_atan2 (y[i], x[i])), <=, TEST_VMATH_LOOP_ULP);
  }

  g_assert_cmpuint (max_ulp, <=, 2);

  g_free (y);
  g_free (x);
  g_free (out);
  ncm_rng_free (rng);
}

static void
test_ncm_vmath_atan2_special (void)
{
  const gdouble v[] = {0.0, -0.0, 1.0, -1.0, 1.0e-310, -1.0e-310, DBL_MAX, -DBL_MAX, INFINITY, -INFINITY};
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (v); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (v); j++)
    {
      const gdouble y = v[i];
      const gdouble x = v[j];

      if (isinf (y) && isinf (x))
        g_assert_true (isnan (ncm_vmath_atan2 (y, x)));
      else
        g_assert_cmpuint (_ulp_diff (ncm_vmath_atan2 (y, x), atan2 (y, x)), <=, 2);
    }
  }

  /* Signed zeros, exactly */
  g_assert_cmpfloat (ncm_vmath_atan2 (0.0, -0.0), ==, M_PI);
  g_assert_cmpfloat (ncm_vmath_atan2 (-0.0, -0.0), ==, -M_PI);
  g_assert_true (signbit (ncm_vmath_atan2 (-0.0, 1.0)));
  g_assert_false (signbit (ncm_vmath_atan2 (0.0, 0.0)));

  g_assert_true (isnan (ncm_vmath_atan2 (NAN, 1.0)));
  g_assert_true (isnan (ncm_vmath_atan2 (1.0, NAN)));
}

static void
test_ncm_vmath_exp_random (void)
{
  const gdouble ranges[][2] = {
    {
      -746.0, 710.0
    }, {
      -30.0, 5.0
    }, {
      -1.0e-6, 1.0e-6
    }, {
      -745.2, -700.0
    }, {
      700.0, 709.8
    }
  };
  NcmRNG *rng  = ncm_rng_seeded_new (NULL, 13);
  gdouble *x   = g_new (gdouble, TEST_VMATH_N);
  gdouble *out = g_new (gdouble, TEST_VMATH_N);
  guint k, i;

  for (k = 0; k < G_N_ELEMENTS (ranges); k++)
  {
    guint64 max_ulp = 0;

    for (i = 0; i < TEST_VMATH_N; i++)
      x[i] = ncm_rng_uniform_gen (rng, ranges[k][0], ranges[k][1]);

    _exp_array (TEST_VMATH_N, x, out);

    for (i = 0; i < TEST_VMATH_N; i++)
    {
      max_ulp = MAX (max_ulp, _ulp_diff (out[i], exp (x[i])));
      TEST_VMATH_ASSERT_LOOP (out[i], ncm_vmath_exp (x[i]));
    }

#ifdef NCM_VMATH_HAVE_FAST_EXP
    g_assert_cmpuint (max_ulp, <=, 2);
#else
    g_assert_cmpuint (max_ulp, <=, TEST_VMATH_LIBMVEC_ULP);
#endif
  }

  g_free (x);
  g_free (out);
  ncm_rng_free (rng);
}

static void
test_ncm_vmath_exp_special (void)
{
  const gdouble v[] = {0.0, -0.0, 1.0, -1.0, 709.78, 709.79, 710.0, 800.0, -708.3, -708.4, -745.13, -745.14, -746.0, -800.0, 1.0e-300};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (v); i++)
    g_assert_cmpuint (_ulp_diff (ncm_vmath_exp (v[i]), exp (v[i])), <=, 2);

  g_assert_cmpfloat (ncm_vmath_exp (0.0), ==, 1.0);
  g_assert_cmpfloat (ncm_vmath_exp (-745.13), ==, exp (-745.13)); /* smallest subnormal */
  g_assert_cmpfloat (ncm_vmath_exp (-745.14), ==, 0.0);
  g_assert_cmpfloat (ncm_vmath_exp (-INFINITY), ==, 0.0);
  g_assert_true (isinf (ncm_vmath_exp (709.79)));
  g_assert_true (isinf (ncm_vmath_exp (INFINITY)));
  g_assert_true (isnan (ncm_vmath_exp (NAN)));
}

static void
test_ncm_vmath_log_random (void)
{
  /* log10 of x: every scale, subnormals, and close to 1 */
  const gdouble ranges[][2] = {
    {
      -300.0, 300.0
    }, {
      -1.0, 1.0
    }, {
      -1.0e-4, 1.0e-4
    }, {
      -323.0, -308.0
    }, {
      307.0, 308.2
    }
  };
  NcmRNG *rng  = ncm_rng_seeded_new (NULL, 17);
  gdouble *x   = g_new (gdouble, TEST_VMATH_N);
  gdouble *out = g_new (gdouble, TEST_VMATH_N);
  guint k, i;

  for (k = 0; k < G_N_ELEMENTS (ranges); k++)
  {
    guint64 max_ulp = 0;

    for (i = 0; i < TEST_VMATH_N; i++)
      x[i] = pow (10.0, ncm_rng_uniform_gen (rng, ranges[k][0], ranges[k][1]));

    _log_array (TEST_VMATH_N, x, out);

    for (i = 0; i < TEST_VMATH_N; i++)
    {
      max_ulp = MAX (max_ulp, _ulp_diff (out[i], log (x[i])));
      TEST_VMATH_ASSERT_LOOP (out[i], ncm_vmath_log (x[i]));
    }

#ifdef NCM_VMATH_HAVE_FAST_EXP
    g_assert_cmpuint (max_ulp, <=, 1);
#else
    g_assert_cmpuint (max_ulp, <=, TEST_VMATH_LIBMVEC_ULP);
#endif
  }

  g_free (x);
  g_free (out);
  ncm_rng_free (rng);
}

static void
test_ncm_vmath_log_special (void)
{
  const gdouble v[] = {1.0, 2.0, 0.5, M_SQRT1_2, M_SQRT2, DBL_MIN, 4.9e-324, 1.0e-310, DBL_MAX};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (v); i++)
    g_assert_cmpuint (_ulp_diff (ncm_vmath_log (v[i]), log (v[i])), <=, 1);

  g_assert_cmpfloat (ncm_vmath_log (1.0), ==, 0.0);
  g_assert_true (isinf (ncm_vmath_log (0.0)) && (ncm_vmath_log (0.0) < 0.0));
  g_assert_true (isinf (ncm_vmath_log (-0.0)) && (ncm_vmath_log (-0.0) < 0.0));
  g_assert_true (isinf (ncm_vmath_log (INFINITY)) && (ncm_vmath_log (INFINITY) > 0.0));
  g_assert_true (isnan (ncm_vmath_log (-1.0)));
  g_assert_true (isnan (ncm_vmath_log (-INFINITY)));
  g_assert_true (isnan (ncm_vmath_log (NAN)));
}

static void
test_ncm_vmath_pow_pos_random (void)
{
  NcmRNG *rng  = ncm_rng_seeded_new (NULL, 19);
  gdouble *x   = g_new (gdouble, TEST_VMATH_N);
  gdouble *y   = g_new (gdouble, TEST_VMATH_N);
  gdouble *out = g_new (gdouble, TEST_VMATH_N);
  guint i;

  /* The ranges of the integrals: bases in (0, 2], exponents of a few units */
  for (i = 0; i < TEST_VMATH_N; i++)
  {
    x[i] = pow (10.0, ncm_rng_uniform_gen (rng, -6.0, log10 (2.0)));
    y[i] = ncm_rng_uniform_gen (rng, -8.0, 8.0);
  }

  _pow_pos_array (TEST_VMATH_N, x, y, out);

  for (i = 0; i < TEST_VMATH_N; i++)
  {
    const gdouble ref = pow (x[i], y[i]);

    TEST_VMATH_ASSERT_LOOP (out[i], ncm_vmath_pow_pos (x[i], y[i]));
#ifdef NCM_VMATH_HAVE_FAST_EXP
    /* Within (2 + |y ln x|) ulp, with one more for the rounding of y ln x */
    g_assert_cmpfloat (fabs (out[i] / ref - 1.0), <=, (3.0 + fabs (y[i] * log (x[i]))) * DBL_EPSILON);
#else
    g_assert_cmpuint (_ulp_diff (out[i], ref), <=, TEST_VMATH_LIBMVEC_ULP);
#endif
  }

  g_free (x);
  g_free (y);
  g_free (out);
  ncm_rng_free (rng);
}

/*
 * Within 2 ulp of sin() and cos() where |value| > 1e-3 and within 2.3e-16 elsewhere, on
 * [-2 pi, 2 pi] and on ranges up to 1e5; the same values in and out of a vectorized loop.
 */
static void
test_ncm_vmath_sincos_random (void)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 17);
  gdouble *x        = g_new (gdouble, TEST_VMATH_N);
  gdouble *s        = g_new (gdouble, TEST_VMATH_N);
  gdouble *c        = g_new (gdouble, TEST_VMATH_N);
  const gdouble R[] = {2.0 * M_PI, 100.0, 1.0e5};
  guint k;

  for (k = 0; k < G_N_ELEMENTS (R); k++)
  {
    guint64 max_ulp = 0;
    gdouble max_abs = 0.0;
    guint i;

    for (i = 0; i < TEST_VMATH_N; i++)
      x[i] = ncm_rng_uniform_gen (rng, -R[k], R[k]);

    _sincos_array (TEST_VMATH_N, x, s, c);

    for (i = 0; i < TEST_VMATH_N; i++)
    {
      const gdouble sl = sin (x[i]);
      const gdouble cl = cos (x[i]);
      gdouble s1, c1;

      ncm_vmath_sincos (x[i], &s1, &c1);
      g_assert_cmpfloat (s[i], ==, s1);
      g_assert_cmpfloat (c[i], ==, c1);

      if (fabs (sl) > 1.0e-3)
        max_ulp = MAX (max_ulp, _ulp_diff (s[i], sl));

      if (fabs (cl) > 1.0e-3)
        max_ulp = MAX (max_ulp, _ulp_diff (c[i], cl));

      max_abs = MAX (max_abs, MAX (fabs (s[i] - sl), fabs (c[i] - cl)));
    }

    g_assert_cmpuint (max_ulp, <=, 2);
    g_assert_cmpfloat (max_abs, <=, 2.3e-16);
  }

  g_free (x);
  g_free (s);
  g_free (c);
  ncm_rng_free (rng);
}

/* Exact values at 0 and -0, near the quadrant boundaries, and NaN for NaN and infinities */
static void
test_ncm_vmath_sincos_special (void)
{
  const gdouble xs[] = {0.0, -0.0, M_PI_2, M_PI, 3.0 * M_PI_2, -M_PI_2, -M_PI, M_PI_4, -M_PI_4};
  gdouble s, c;
  guint i;

  ncm_vmath_sincos (0.0, &s, &c);
  g_assert_cmpfloat (s, ==, 0.0);
  g_assert_false (signbit (s));
  g_assert_cmpfloat (c, ==, 1.0);

  /* sin (-0) is +0, not -0 as in libm */
  ncm_vmath_sincos (-0.0, &s, &c);
  g_assert_cmpfloat (s, ==, 0.0);
  g_assert_cmpfloat (c, ==, 1.0);

  for (i = 0; i < G_N_ELEMENTS (xs); i++)
  {
    ncm_vmath_sincos (xs[i], &s, &c);
    g_assert_cmpfloat (fabs (s - sin (xs[i])), <=, 2.3e-16);
    g_assert_cmpfloat (fabs (c - cos (xs[i])), <=, 2.3e-16);
  }

  ncm_vmath_sincos (NAN, &s, &c);
  g_assert_true (isnan (s) && isnan (c));
  ncm_vmath_sincos (INFINITY, &s, &c);
  g_assert_true (isnan (s) && isnan (c));
  ncm_vmath_sincos (-INFINITY, &s, &c);
  g_assert_true (isnan (s) && isnan (c));
}


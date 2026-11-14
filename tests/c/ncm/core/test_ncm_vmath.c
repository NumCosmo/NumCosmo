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

static void test_ncm_vmath_select (void);
static void test_ncm_vmath_atan2_random (void);
static void test_ncm_vmath_atan2_special (void);
static void test_ncm_vmath_exp_random (void);
static void test_ncm_vmath_exp_special (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/vmath/select", &test_ncm_vmath_select);
  g_test_add_func ("/ncm/vmath/atan2/random", &test_ncm_vmath_atan2_random);
  g_test_add_func ("/ncm/vmath/atan2/special", &test_ncm_vmath_atan2_special);
  g_test_add_func ("/ncm/vmath/exp/random", &test_ncm_vmath_exp_random);
  g_test_add_func ("/ncm/vmath/exp/special", &test_ncm_vmath_exp_special);

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
    g_assert_cmpfloat (out[i], ==, ncm_vmath_atan2 (y[i], x[i]));
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
      g_assert_cmpfloat (out[i], ==, ncm_vmath_exp (x[i]));
    }

#ifdef NCM_VMATH_HAVE_FAST_EXP
    g_assert_cmpuint (max_ulp, <=, 2);
#else
    g_assert_cmpuint (max_ulp, ==, 0);
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


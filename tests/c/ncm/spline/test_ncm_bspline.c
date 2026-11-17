/***************************************************************************
 *            test_ncm_bspline.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_bspline.c
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

/* Points compared in each case of test_ncm_bspline_eval */
#define TEST_NCM_BSPLINE_NPOINTS 1000

void test_ncm_bspline_new (void);
void test_ncm_bspline_getters (void);
void test_ncm_bspline_eval (void);
void test_ncm_bspline_dup_is_a_copy (void);
void test_ncm_bspline_serialize (void);
void test_ncm_bspline_missing (void);
void test_ncm_bspline_missing_subprocess (void);
void test_ncm_bspline_stride (void);
void test_ncm_bspline_stride_subprocess (void);
void test_ncm_bspline_lengths (void);
void test_ncm_bspline_lengths_subprocess (void);
void test_ncm_bspline_decreasing (void);
void test_ncm_bspline_decreasing_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/bspline/new", &test_ncm_bspline_new);
  g_test_add_func ("/ncm/bspline/getters", &test_ncm_bspline_getters);
  g_test_add_func ("/ncm/bspline/eval", &test_ncm_bspline_eval);
  g_test_add_func ("/ncm/bspline/dup_is_a_copy", &test_ncm_bspline_dup_is_a_copy);
  g_test_add_func ("/ncm/bspline/serialize", &test_ncm_bspline_serialize);
  g_test_add_func ("/ncm/bspline/missing", &test_ncm_bspline_missing);
  g_test_add_func ("/ncm/bspline/missing/subprocess", &test_ncm_bspline_missing_subprocess);
  g_test_add_func ("/ncm/bspline/stride", &test_ncm_bspline_stride);
  g_test_add_func ("/ncm/bspline/stride/subprocess", &test_ncm_bspline_stride_subprocess);
  g_test_add_func ("/ncm/bspline/lengths", &test_ncm_bspline_lengths);
  g_test_add_func ("/ncm/bspline/lengths/subprocess", &test_ncm_bspline_lengths_subprocess);
  g_test_add_func ("/ncm/bspline/decreasing", &test_ncm_bspline_decreasing);
  g_test_add_func ("/ncm/bspline/decreasing/subprocess", &test_ncm_bspline_decreasing_subprocess);

  g_test_run ();
}

/*
 * A prepared NcmSplineBSpline of the given order through n samples of
 * exp (-x^2 / 2) cos (3 x) at points of [-2, 3] with unequal spacing, multiplied by scale.
 */
static NcmSplineBSpline *
_test_ncm_bspline_new_spline (const guint order, const guint n, const gdouble scale)
{
  NcmVector *xv = ncm_vector_new (n);
  NcmVector *yv = ncm_vector_new (n);
  NcmSplineBSpline *sbs;
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble u = i / (n - 1.0);
    const gdouble x = -2.0 + 5.0 * u * (0.5 + 0.5 * u);

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, scale * exp (-0.5 * x * x) * cos (3.0 * x));
  }

  sbs = ncm_spline_bspline_new_full (order, xv, yv, TRUE);

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return sbs;
}

/* A vector with the n values of d */
static NcmVector *
_test_ncm_bspline_vector (const gdouble *d, const guint n)
{
  NcmVector *v = ncm_vector_new (n);
  guint i;

  for (i = 0; i < n; i++)
    ncm_vector_set (v, i, d[i]);

  return v;
}

/*
 * The order-2 B-spline with knots (0, 0, 1, 2, 2) and control points (1, 3, 2) is the
 * polygon through (0, 1), (1, 3), (2, 2); outside [0, 2] it continues the edge segment.
 */
void
test_ncm_bspline_new (void)
{
  const gdouble t[5]        = {0.0, 0.0, 1.0, 2.0, 2.0};
  const gdouble c[3]        = {1.0, 3.0, 2.0};
  NcmVector *knots          = _test_ncm_bspline_vector (t, 5);
  NcmVector *control_points = _test_ncm_bspline_vector (c, 3);
  NcmBSpline *bs            = ncm_bspline_new (2, knots, control_points);

  g_assert_true (NCM_IS_BSPLINE (bs));

  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, 0.0), ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, 0.5), ==, 2.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, 1.0), ==, 3.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, 1.5), ==, 2.5, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, 2.0), ==, 2.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, -0.5), ==, 0.0, 0.0, 1.0e-15);
  ncm_assert_cmpdouble_e (ncm_bspline_eval (bs, 2.5), ==, 1.5, 1.0e-15, 0.0);

  ncm_vector_free (knots);
  ncm_vector_free (control_points);
  ncm_bspline_free (bs);
}

void
test_ncm_bspline_getters (void)
{
  const guint order         = 6;
  const guint n             = 40;
  NcmSplineBSpline *sbs     = _test_ncm_bspline_new_spline (order, n, 1.0);
  NcmBSpline *bs            = ncm_spline_bspline_dup_bspline (sbs);
  NcmVector *brk            = ncm_spline_bspline_get_breakpoints (sbs);
  NcmVector *knots          = ncm_bspline_peek_knots (bs);
  NcmVector *control_points = ncm_bspline_peek_control_points (bs);
  NcmBSpline *bs2           = ncm_bspline_new (order, knots, control_points);
  gdouble x_lo, x_hi;

  g_assert_cmpuint (ncm_bspline_get_order (bs), ==, order);
  g_assert_cmpuint (ncm_vector_len (knots), ==, n + order);
  g_assert_cmpuint (ncm_vector_len (control_points), ==, n);

  ncm_bspline_get_range (bs, &x_lo, &x_hi);
  g_assert_cmpfloat (x_lo, ==, ncm_vector_get (brk, 0));
  g_assert_cmpfloat (x_hi, ==, ncm_vector_get (brk, ncm_vector_len (brk) - 1));
  g_assert_cmpfloat (x_lo, ==, -2.0);
  g_assert_cmpfloat (x_hi, ==, 3.0);

  /* The vectors given at construction are referenced */
  g_assert_true (ncm_bspline_peek_knots (bs2) == knots);
  g_assert_true (ncm_bspline_peek_control_points (bs2) == control_points);

  ncm_bspline_free (bs2);
  ncm_vector_free (brk);
  ncm_bspline_free (bs);
  ncm_spline_free (NCM_SPLINE (sbs));
}

/*
 * NcmBSpline and the NcmSplineBSpline it was taken from run the same evaluation on the
 * same knots and control points, so their values are equal, not close: at random points
 * of the range and of a tenth of its length beyond each end, and at the breakpoints.
 */
void
test_ncm_bspline_eval (void)
{
  const guint n_cases[3] = {12, 57, 200};
  NcmRNG *rng            = ncm_rng_seeded_new (NULL, 20261117);
  guint order, j, i;

  for (order = 2; order <= NCM_SPLINE_BSPLINE_MAX_ORDER; order++)
  {
    for (j = 0; j < 3; j++)
    {
      NcmSplineBSpline *sbs = _test_ncm_bspline_new_spline (order, n_cases[j], 1.0);
      NcmBSpline *bs        = ncm_spline_bspline_dup_bspline (sbs);
      NcmVector *brk        = ncm_spline_bspline_get_breakpoints (sbs);
      gdouble x_lo, x_hi;

      ncm_bspline_get_range (bs, &x_lo, &x_hi);

      for (i = 0; i < TEST_NCM_BSPLINE_NPOINTS; i++)
      {
        const gdouble x = ncm_rng_uniform_gen (rng, x_lo - 0.1 * (x_hi - x_lo), x_hi + 0.1 * (x_hi - x_lo));

        g_assert_cmpfloat (ncm_bspline_eval (bs, x), ==, ncm_spline_eval (NCM_SPLINE (sbs), x));
      }

      for (i = 0; i < ncm_vector_len (brk); i++)
      {
        const gdouble x = ncm_vector_get (brk, i);

        g_assert_cmpfloat (ncm_bspline_eval (bs, x), ==, ncm_spline_eval (NCM_SPLINE (sbs), x));
      }

      ncm_vector_free (brk);
      ncm_bspline_free (bs);
      ncm_spline_free (NCM_SPLINE (sbs));
    }
  }

  ncm_rng_free (rng);
}

/* Setting other values in the NcmSplineBSpline does not change an NcmBSpline taken before */
void
test_ncm_bspline_dup_is_a_copy (void)
{
  const guint n         = 30;
  NcmSplineBSpline *sbs = _test_ncm_bspline_new_spline (4, n, 1.0);
  NcmSplineBSpline *ref = _test_ncm_bspline_new_spline (4, n, 1.0);
  NcmSplineBSpline *alt = _test_ncm_bspline_new_spline (4, n, 2.0);
  NcmBSpline *bs        = ncm_spline_bspline_dup_bspline (sbs);
  guint i;

  ncm_spline_set (NCM_SPLINE (sbs), ncm_spline_peek_xv (NCM_SPLINE (alt)), ncm_spline_peek_yv (NCM_SPLINE (alt)), TRUE);

  for (i = 0; i < 101; i++)
  {
    const gdouble x = -2.0 + 5.0 * i / 100.0;

    g_assert_cmpfloat (ncm_bspline_eval (bs, x), ==, ncm_spline_eval (NCM_SPLINE (ref), x));
    g_assert_cmpfloat (ncm_spline_eval (NCM_SPLINE (sbs), x), ==, ncm_spline_eval (NCM_SPLINE (alt), x));
  }

  g_assert_cmpfloat (ncm_bspline_eval (bs, 0.0), !=, ncm_spline_eval (NCM_SPLINE (sbs), 0.0));

  ncm_bspline_free (bs);
  ncm_spline_free (NCM_SPLINE (alt));
  ncm_spline_free (NCM_SPLINE (ref));
  ncm_spline_free (NCM_SPLINE (sbs));
}

void
test_ncm_bspline_serialize (void)
{
  NcmSplineBSpline *sbs = _test_ncm_bspline_new_spline (6, 57, 1.0);
  NcmBSpline *bs        = ncm_spline_bspline_dup_bspline (sbs);
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmBSpline *dup       = NCM_BSPLINE (ncm_serialize_dup_obj (ser, G_OBJECT (bs)));
  guint i;

  g_assert_true (NCM_IS_BSPLINE (dup));
  g_assert_true (dup != bs);
  g_assert_cmpuint (ncm_bspline_get_order (dup), ==, 6);
  g_assert_cmpuint (ncm_vector_len (ncm_bspline_peek_knots (dup)), ==, ncm_vector_len (ncm_bspline_peek_knots (bs)));
  g_assert_cmpuint (ncm_vector_len (ncm_bspline_peek_control_points (dup)), ==, ncm_vector_len (ncm_bspline_peek_control_points (bs)));

  for (i = 0; i < ncm_vector_len (ncm_bspline_peek_knots (bs)); i++)
    g_assert_cmpfloat (ncm_vector_get (ncm_bspline_peek_knots (dup), i), ==, ncm_vector_get (ncm_bspline_peek_knots (bs), i));

  for (i = 0; i < ncm_vector_len (ncm_bspline_peek_control_points (bs)); i++)
    g_assert_cmpfloat (ncm_vector_get (ncm_bspline_peek_control_points (dup), i), ==, ncm_vector_get (ncm_bspline_peek_control_points (bs), i));

  for (i = 0; i < 101; i++)
  {
    const gdouble x = -2.5 + 6.0 * i / 100.0;

    g_assert_cmpfloat (ncm_bspline_eval (dup, x), ==, ncm_bspline_eval (bs, x));
  }

  ncm_bspline_free (dup);
  ncm_serialize_free (ser);
  ncm_bspline_free (bs);
  ncm_spline_free (NCM_SPLINE (sbs));
}

void
test_ncm_bspline_missing (void)
{
  g_test_trap_subprocess ("/ncm/bspline/missing/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcmBSpline: both the knots and the control points are required*");
}

void
test_ncm_bspline_missing_subprocess (void)
{
  const gdouble t[5] = {0.0, 0.0, 1.0, 2.0, 2.0};
  NcmVector *knots   = _test_ncm_bspline_vector (t, 5);

  ncm_bspline_new (2, knots, NULL);
}

void
test_ncm_bspline_stride (void)
{
  g_test_trap_subprocess ("/ncm/bspline/stride/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcmBSpline: the knots and the control points must have stride one*");
}

void
test_ncm_bspline_stride_subprocess (void)
{
  const gdouble t[5]        = {0.0, 0.0, 1.0, 2.0, 2.0};
  gdouble c[6]              = {1.0, 0.0, 3.0, 0.0, 2.0, 0.0};
  NcmVector *knots          = _test_ncm_bspline_vector (t, 5);
  NcmVector *control_points = ncm_vector_new_data_static (c, 3, 2);

  ncm_bspline_new (2, knots, control_points);
}

void
test_ncm_bspline_lengths (void)
{
  g_test_trap_subprocess ("/ncm/bspline/lengths/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcmBSpline: order 2 needs at least 2 control points and as many knots as control points plus 2, got 3 control points and 4 knots*");
}

void
test_ncm_bspline_lengths_subprocess (void)
{
  const gdouble t[4]        = {0.0, 0.0, 1.0, 2.0};
  const gdouble c[3]        = {1.0, 3.0, 2.0};
  NcmVector *knots          = _test_ncm_bspline_vector (t, 4);
  NcmVector *control_points = _test_ncm_bspline_vector (c, 3);

  ncm_bspline_new (2, knots, control_points);
}

void
test_ncm_bspline_decreasing (void)
{
  g_test_trap_subprocess ("/ncm/bspline/decreasing/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcmBSpline: the knots must not decrease, knot 3 is*");
}

void
test_ncm_bspline_decreasing_subprocess (void)
{
  const gdouble t[5]        = {0.0, 0.0, 1.0, 0.5, 2.0};
  const gdouble c[3]        = {1.0, 3.0, 2.0};
  NcmVector *knots          = _test_ncm_bspline_vector (t, 5);
  NcmVector *control_points = _test_ncm_bspline_vector (c, 3);

  ncm_bspline_new (2, knots, control_points);
}


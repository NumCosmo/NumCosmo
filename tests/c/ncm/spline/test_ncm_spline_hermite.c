/***************************************************************************
 *            test_ncm_spline_hermite.c
 *
 *  Mon Oct 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

typedef struct
{
  gdouble (*f) (gdouble x);
  gdouble (*d1) (gdouble x);
  gdouble (*d2) (gdouble x);
} TestFunc;

static void test_ncm_spline_hermite_exact (void);
static void test_ncm_spline_hermite_convergence (void);
static void test_ncm_spline_hermite_knots (void);
static void test_ncm_spline_hermite_integ (void);
static void test_ncm_spline_hermite_serialize (void);
static void test_ncm_spline_hermite_copy_empty_reset (void);
static void test_ncm_spline_hermite_prepare_no_derivs_subprocess (void);
static void test_ncm_spline_hermite_prepare_no_derivs (void);
static void test_ncm_spline_hermite_derivs_no_knots_subprocess (void);
static void test_ncm_spline_hermite_derivs_no_knots (void);
static void test_ncm_spline_hermite_derivs_length_d1_subprocess (void);
static void test_ncm_spline_hermite_derivs_length_d2_subprocess (void);
static void test_ncm_spline_hermite_derivs_length_d1 (void);
static void test_ncm_spline_hermite_derivs_length_d2 (void);
static gdouble _poly (const gdouble x);
static gdouble _poly_d1 (const gdouble x);
static gdouble _poly_d2 (const gdouble x);
static gdouble _sin_d2 (const gdouble x);

/* A quintic, reproduced exactly on any knots */
static const gdouble _poly_c[6] = {1.0, 0.5, -2.0, 0.3, -0.7, 0.2};

static const TestFunc _poly_func = {_poly, _poly_d1, _poly_d2};

static const TestFunc _sin_func = {sin, cos, _sin_d2};

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/spline_hermite/exact", &test_ncm_spline_hermite_exact);
  g_test_add_func ("/ncm/spline_hermite/convergence", &test_ncm_spline_hermite_convergence);
  g_test_add_func ("/ncm/spline_hermite/knots", &test_ncm_spline_hermite_knots);
  g_test_add_func ("/ncm/spline_hermite/integ", &test_ncm_spline_hermite_integ);
  g_test_add_func ("/ncm/spline_hermite/serialize", &test_ncm_spline_hermite_serialize);
  g_test_add_func ("/ncm/spline_hermite/copy_empty_reset", &test_ncm_spline_hermite_copy_empty_reset);
  g_test_add_func ("/ncm/spline_hermite/prepare_no_derivs", &test_ncm_spline_hermite_prepare_no_derivs);
  g_test_add_func ("/ncm/spline_hermite/prepare_no_derivs/subprocess", &test_ncm_spline_hermite_prepare_no_derivs_subprocess);
  g_test_add_func ("/ncm/spline_hermite/derivs_no_knots", &test_ncm_spline_hermite_derivs_no_knots);
  g_test_add_func ("/ncm/spline_hermite/derivs_no_knots/subprocess", &test_ncm_spline_hermite_derivs_no_knots_subprocess);
  g_test_add_func ("/ncm/spline_hermite/derivs_length_d1", &test_ncm_spline_hermite_derivs_length_d1);
  g_test_add_func ("/ncm/spline_hermite/derivs_length_d1/subprocess", &test_ncm_spline_hermite_derivs_length_d1_subprocess);
  g_test_add_func ("/ncm/spline_hermite/derivs_length_d2", &test_ncm_spline_hermite_derivs_length_d2);
  g_test_add_func ("/ncm/spline_hermite/derivs_length_d2/subprocess", &test_ncm_spline_hermite_derivs_length_d2_subprocess);

  g_test_run ();

  return 0;
}

static gdouble
_poly (const gdouble x)
{
  return _poly_c[0] + x * (_poly_c[1] + x * (_poly_c[2] + x * (_poly_c[3] + x * (_poly_c[4] + x * _poly_c[5]))));
}

static gdouble
_poly_d1 (const gdouble x)
{
  return _poly_c[1] + x * (2.0 * _poly_c[2] + x * (3.0 * _poly_c[3] + x * (4.0 * _poly_c[4] + x * 5.0 * _poly_c[5])));
}

static gdouble
_poly_d2 (const gdouble x)
{
  return 2.0 * _poly_c[2] + x * (6.0 * _poly_c[3] + x * (12.0 * _poly_c[4] + x * 20.0 * _poly_c[5]));
}

static gdouble
_poly_integ (const gdouble x)
{
  return x * (_poly_c[0] + x * (_poly_c[1] / 2.0 + x * (_poly_c[2] / 3.0 + x * (_poly_c[3] / 4.0 + x * (_poly_c[4] / 5.0 + x * _poly_c[5] / 6.0)))));
}

static gdouble
_sin_d2 (const gdouble x)
{
  return -sin (x);
}

/* Knots x_i = xf (i / (n - 1))^p, uniform for p = 1 */
static NcmSplineHermite *
_new_on (const TestFunc *F, const guint n, const gdouble xf, const gdouble p, const gboolean init)
{
  NcmVector *xv  = ncm_vector_new (n);
  NcmVector *yv  = ncm_vector_new (n);
  NcmVector *d1v = ncm_vector_new (n);
  NcmVector *d2v = ncm_vector_new (n);
  NcmSplineHermite *sh;
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble x = xf * pow (i / (n - 1.0), p);

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, F->f (x));
    ncm_vector_set (d1v, i, F->d1 (x));
    ncm_vector_set (d2v, i, F->d2 (x));
  }

  sh = ncm_spline_hermite_new_full (xv, yv, d1v, d2v, init);

  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_vector_free (d1v);
  ncm_vector_free (d2v);

  return sh;
}

static gdouble
_max_midpoint_error (NcmSpline *s, const TestFunc *F)
{
  NcmVector *xv = ncm_spline_peek_xv (s);
  const guint n = ncm_vector_len (xv);
  gdouble err   = 0.0;
  guint i;

  for (i = 0; i + 1 < n; i++)
  {
    const gdouble xm = 0.5 * (ncm_vector_get (xv, i) + ncm_vector_get (xv, i + 1));

    err = GSL_MAX (err, fabs (ncm_spline_eval (s, xm) - F->f (xm)));
  }

  return err;
}

/* A quintic is reproduced with its derivatives and integral on non-uniform knots, and on two knots */
static void
test_ncm_spline_hermite_exact (void)
{
  const guint ns[] = {12, 2};
  guint k;

  for (k = 0; k < G_N_ELEMENTS (ns); k++)
  {
    NcmSplineHermite *sh = _new_on (&_poly_func, ns[k], 3.0, 1.5, TRUE);
    NcmSpline *s         = NCM_SPLINE (sh);
    guint i;

    g_assert_cmpstr (NCM_SPLINE_GET_CLASS (s)->name (s), ==, "NcmSplineHermite");

    for (i = 0; i <= 50; i++)
    {
      const gdouble x = 3.0 * i / 50.0;

      ncm_assert_cmpdouble_e (ncm_spline_eval (s, x), ==, _poly (x), 1.0e-12, 1.0e-12);
      ncm_assert_cmpdouble_e (ncm_spline_eval_deriv (s, x), ==, _poly_d1 (x), 1.0e-11, 1.0e-11);
      ncm_assert_cmpdouble_e (ncm_spline_eval_deriv2 (s, x), ==, _poly_d2 (x), 1.0e-10, 1.0e-10);
      /* The fifth derivative from data at spacing h carries roundoff of order eps / h^5 */
      ncm_assert_cmpdouble_e (ncm_spline_eval_deriv_nmax (s, x), ==, 120.0 * _poly_c[5], 1.0e-8, 1.0e-8);
      ncm_assert_cmpdouble_e (ncm_spline_eval_integ (s, 0.2, x), ==, _poly_integ (x) - _poly_integ (0.2), 1.0e-11, 1.0e-11);
    }

    ncm_spline_free (s);
  }
}

/* The error of a smooth function falls as h^6: halving h divides it by about 64 */
static void
test_ncm_spline_hermite_convergence (void)
{
  NcmSplineHermite *sh1 = _new_on (&_sin_func, 11, 2.0 * M_PI, 1.0, TRUE);
  NcmSplineHermite *sh2 = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, TRUE);
  const gdouble e1      = _max_midpoint_error (NCM_SPLINE (sh1), &_sin_func);
  const gdouble e2      = _max_midpoint_error (NCM_SPLINE (sh2), &_sin_func);

  g_assert_cmpfloat (e1, <, 1.0e-5);
  g_assert_cmpfloat (e2, <, 1.0e-7);
  g_assert_cmpfloat (e2, <, e1 / 32.0);

  ncm_spline_free (NCM_SPLINE (sh1));
  ncm_spline_free (NCM_SPLINE (sh2));
}

/* Values and derivatives at the knots are the given ones; the interval index is honored */
static void
test_ncm_spline_hermite_knots (void)
{
  NcmSplineHermite *sh = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, TRUE);
  NcmSpline *s         = NCM_SPLINE (sh);
  NcmVector *xv        = ncm_spline_peek_xv (s);
  const guint n        = ncm_vector_len (xv);
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble x = ncm_vector_get (xv, i);

    /* The last knot is the right end of its interval, matched to roundoff */
    ncm_assert_cmpdouble_e (ncm_spline_eval (s, x), ==, sin (x), 1.0e-14, 1.0e-14);
    ncm_assert_cmpdouble_e (ncm_spline_eval_deriv (s, x), ==, cos (x), 1.0e-13, 1.0e-13);
    ncm_assert_cmpdouble_e (ncm_spline_eval_deriv2 (s, x), ==, -sin (x), 1.0e-12, 1.0e-12);
  }

  for (i = 0; i + 1 < n; i++)
  {
    const gdouble xm = 0.5 * (ncm_vector_get (xv, i) + ncm_vector_get (xv, i + 1));

    g_assert_cmpfloat (ncm_spline_eval_idx (s, xm, i), ==, ncm_spline_eval (s, xm));
    g_assert_cmpfloat (ncm_spline_eval_deriv_idx (s, xm, i), ==, ncm_spline_eval_deriv (s, xm));
  }

  /* Outside the knots the boundary quintic continues */
  {
    const gdouble x = 2.0 * M_PI + 0.1;

    g_assert_cmpfloat (ncm_spline_eval (s, x), ==, ncm_spline_eval_idx (s, x, n - 2));
  }

  ncm_spline_free (s);
}

/* The integral over several intervals matches the function's, through both entry points */
static void
test_ncm_spline_hermite_integ (void)
{
  NcmSplineHermite *sh    = _new_on (&_sin_func, 41, 2.0 * M_PI, 1.0, TRUE);
  NcmSpline *s            = NCM_SPLINE (sh);
  const gdouble a         = 0.3;
  const gdouble b         = 5.1;
  const gdouble ref       = cos (a) - cos (b);
  const gdouble integ     = ncm_spline_eval_integ (s, a, b);
  const gdouble integ_idx = ncm_spline_eval_integ_idx (s, a, ncm_spline_get_index (s, a), b, ncm_spline_get_index (s, b));

  ncm_assert_cmpdouble_e (integ, ==, ref, 1.0e-9, 0.0);
  g_assert_cmpfloat (integ_idx, ==, integ);

  /* Within one interval */
  ncm_assert_cmpdouble_e (ncm_spline_eval_integ (s, 0.2, 0.3), ==, cos (0.2) - cos (0.3), 1.0e-9, 0.0);

  ncm_spline_free (s);
}

/* The derivative vectors are properties: a serialized copy evaluates identically */
static void
test_ncm_spline_hermite_serialize (void)
{
  NcmSplineHermite *sh = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, TRUE);
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpline *dup       = NCM_SPLINE (ncm_serialize_dup_obj (ser, G_OBJECT (sh)));
  NcmVector *d1v       = NULL;
  NcmVector *d2v       = NULL;
  guint i;

  g_assert_true (NCM_IS_SPLINE_HERMITE (dup));
  g_object_get (dup, "d1", &d1v, "d2", &d2v, NULL);
  g_assert_cmpuint (ncm_vector_len (d1v), ==, 21);
  g_assert_cmpuint (ncm_vector_len (d2v), ==, 21);
  g_assert_true (d1v == ncm_spline_hermite_peek_d1v (NCM_SPLINE_HERMITE (dup)));
  g_assert_true (d2v == ncm_spline_hermite_peek_d2v (NCM_SPLINE_HERMITE (dup)));

  ncm_spline_prepare (dup);

  for (i = 0; i <= 100; i++)
  {
    const gdouble x = 2.0 * M_PI * i / 100.0;

    g_assert_cmpfloat (ncm_spline_eval (dup, x), ==, ncm_spline_eval (NCM_SPLINE (sh), x));
  }

  ncm_vector_free (d1v);
  ncm_vector_free (d2v);
  ncm_spline_free (dup);
  ncm_serialize_free (ser);
  ncm_spline_free (NCM_SPLINE (sh));
}

/* copy_empty gives an empty spline of the type; knots of the same length keep the
 * derivatives, a different length drops them; NULL clears them */
static void
test_ncm_spline_hermite_copy_empty_reset (void)
{
  NcmSplineHermite *sh = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, TRUE);
  NcmSpline *s         = NCM_SPLINE (sh);
  NcmSpline *empty     = ncm_spline_copy_empty (s);
  NcmVector *yv        = ncm_vector_dup (ncm_spline_peek_yv (s));
  NcmVector *xv_short  = ncm_vector_new (5);
  NcmVector *yv_short  = ncm_vector_new (5);
  guint i;

  g_assert_true (NCM_IS_SPLINE_HERMITE (empty));
  g_assert_cmpuint (ncm_spline_get_len (empty), ==, 0);
  g_assert_null (ncm_spline_hermite_peek_d1v (NCM_SPLINE_HERMITE (empty)));
  g_assert_cmpuint (ncm_spline_min_size (empty), ==, 2);

  ncm_vector_scale (yv, 2.0);
  ncm_spline_set_yv (s, yv, TRUE);
  g_assert_nonnull (ncm_spline_hermite_peek_d1v (sh));
  g_assert_nonnull (ncm_spline_hermite_peek_d2v (sh));

  for (i = 0; i < 5; i++)
  {
    ncm_vector_set (xv_short, i, i);
    ncm_vector_set (yv_short, i, i * i);
  }

  ncm_spline_set (s, xv_short, yv_short, FALSE);
  g_assert_null (ncm_spline_hermite_peek_d1v (sh));
  g_assert_null (ncm_spline_hermite_peek_d2v (sh));

  ncm_spline_hermite_set_derivs (sh, xv_short, yv_short, TRUE);
  g_assert_nonnull (ncm_spline_hermite_peek_d1v (sh));
  ncm_spline_hermite_set_derivs (sh, NULL, NULL, FALSE);
  g_assert_null (ncm_spline_hermite_peek_d1v (sh));
  g_assert_null (ncm_spline_hermite_peek_d2v (sh));

  ncm_vector_free (yv);
  ncm_vector_free (xv_short);
  ncm_vector_free (yv_short);
  ncm_spline_free (empty);
  ncm_spline_free (s);
}

static void
test_ncm_spline_hermite_prepare_no_derivs_subprocess (void)
{
  NcmSplineHermite *sh = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, FALSE);

  ncm_spline_hermite_set_derivs (sh, NULL, NULL, FALSE);
  ncm_spline_prepare (NCM_SPLINE (sh));
}

static void
test_ncm_spline_hermite_prepare_no_derivs (void)
{
  g_test_trap_subprocess ("/ncm/spline_hermite/prepare_no_derivs/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*derivatives at the knots were not set*");
}

static void
test_ncm_spline_hermite_derivs_no_knots_subprocess (void)
{
  NcmSplineHermite *sh = ncm_spline_hermite_new ();
  NcmVector *v         = ncm_vector_new (3);

  ncm_spline_hermite_set_derivs (sh, v, v, FALSE);
}

static void
test_ncm_spline_hermite_derivs_no_knots (void)
{
  g_test_trap_subprocess ("/ncm/spline_hermite/derivs_no_knots/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*knots must be set first*");
}

static void
test_ncm_spline_hermite_derivs_length_d1_subprocess (void)
{
  NcmSplineHermite *sh = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, FALSE);
  NcmVector *v         = ncm_vector_new (3);

  ncm_spline_hermite_set_derivs (sh, v, ncm_spline_hermite_peek_d2v (sh), FALSE);
}

static void
test_ncm_spline_hermite_derivs_length_d2_subprocess (void)
{
  NcmSplineHermite *sh = _new_on (&_sin_func, 21, 2.0 * M_PI, 1.0, FALSE);
  NcmVector *v         = ncm_vector_new (3);

  ncm_spline_hermite_set_derivs (sh, ncm_spline_hermite_peek_d1v (sh), v, FALSE);
}

static void
test_ncm_spline_hermite_derivs_length_d1 (void)
{
  g_test_trap_subprocess ("/ncm/spline_hermite/derivs_length_d1/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*d1 has length 3, the spline has 21 knots*");
}

static void
test_ncm_spline_hermite_derivs_length_d2 (void)
{
  g_test_trap_subprocess ("/ncm/spline_hermite/derivs_length_d2/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*d2 has length 3, the spline has 21 knots*");
}


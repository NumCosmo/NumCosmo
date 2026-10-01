/***************************************************************************
 *            test_ncm_spline_bspline.c
 *
 *  Sat September 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spline_bspline.c
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
#include <gsl/gsl_sf_gamma.h>

/* The Gaussian of the accuracy table in the NcmSplineBSpline docs */
#define GAUSS_MU 1200.0
#define GAUSS_SIGMA 300.0
#define GAUSS_XF 2400.0

static gdouble
_gauss (const gdouble x)
{
  return exp (-0.5 * gsl_pow_2 ((x - GAUSS_MU) / GAUSS_SIGMA));
}

static void
_fill_gauss (NcmSpline *s, const guint n, const gboolean init)
{
  NcmVector *xv = ncm_vector_new (n);
  NcmVector *yv = ncm_vector_new (n);
  guint i;

  for (i = 0; i < n; i++)
  {
    const gdouble x = GAUSS_XF * i / (n - 1.0);

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, _gauss (x));
  }

  ncm_spline_set (s, xv, yv, init);

  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

static NcmSpline *
_new_gauss (const guint order, const guint n)
{
  NcmSpline *s = NCM_SPLINE (ncm_spline_bspline_new (order));

  _fill_gauss (s, n, TRUE);

  return s;
}

static gdouble
_max_midpoint_error (NcmSpline *s)
{
  NcmVector *xv = ncm_spline_peek_xv (s);
  const guint n = ncm_vector_len (xv);
  gdouble err   = 0.0;
  guint i;

  for (i = 0; i + 1 < n; i++)
  {
    const gdouble xm = 0.5 * (ncm_vector_get (xv, i) + ncm_vector_get (xv, i + 1));

    err = GSL_MAX (err, fabs (ncm_spline_eval (s, xm) - _gauss (xm)));
  }

  return err;
}

static void
test_ncm_spline_bspline_order (void)
{
  NcmSplineBSpline *sbs = ncm_spline_bspline_new (4);

  g_assert_cmpuint (ncm_spline_bspline_get_order (sbs), ==, 4);
  ncm_spline_bspline_set_order (sbs, 6);
  g_assert_cmpuint (ncm_spline_bspline_get_order (sbs), ==, 6);

  ncm_spline_free (NCM_SPLINE (sbs));
}

/* The spline reproduces the values at the knots, at every order. */
static void
test_ncm_spline_bspline_interpolates (void)
{
  const guint orders[] = {4, 6, 8};
  guint o;

  for (o = 0; o < G_N_ELEMENTS (orders); o++)
  {
    NcmSpline *s  = _new_gauss (orders[o], 200);
    NcmVector *xv = ncm_spline_peek_xv (s);
    guint i;

    for (i = 0; i < 200; i++)
    {
      const gdouble x = ncm_vector_get (xv, i);

      ncm_assert_cmpdouble_e (ncm_spline_eval (s, x), ==, _gauss (x), 0.0, 1.0e-14);
    }

    ncm_spline_free (s);
  }
}

/* Bounds just above the table in the NcmSplineBSpline docs. */
static void
test_ncm_spline_bspline_accuracy (void)
{
  const struct
  {
    guint order;
    guint n;
    gdouble err;
  } cases[] = {
    {4,  100, 1.0e-06},
    {6,  100, 1.0e-08},
    {8,  100, 1.0e-10},
    {10, 100, 1.0e-12},
    {4,  500, 1.0e-09},
    {6,  500, 1.0e-13},
    {4, 2000, 1.0e-11},
    {6, 2000, 1.0e-14},
    {8, 2000, 1.0e-14},
  };

  guint c;

  for (c = 0; c < G_N_ELEMENTS (cases); c++)
  {
    NcmSpline *s = _new_gauss (cases[c].order, cases[c].n);

    g_assert_cmpfloat (_max_midpoint_error (s), <, cases[c].err);

    ncm_spline_free (s);
  }

  /* The cubic spline stays above machine precision where order 8 reaches it */
  {
    NcmSpline *s4 = _new_gauss (4, 2000);

    g_assert_cmpfloat (_max_midpoint_error (s4), >, 1.0e-12);

    ncm_spline_free (s4);
  }
}

static void
test_ncm_spline_bspline_deriv_integ (void)
{
  NcmSpline *s     = _new_gauss (8, 1000);
  const gdouble x0 = 1350.0;
  const gdouble g  = _gauss (x0);
  const gdouble u  = (x0 - GAUSS_MU) / GAUSS_SIGMA;
  const gdouble d1 = -u / GAUSS_SIGMA * g;
  const gdouble d2 = (u * u - 1.0) / gsl_pow_2 (GAUSS_SIGMA) * g;
  const gdouble a  = GAUSS_SIGMA * M_SQRT2;
  const gdouble ex = GAUSS_SIGMA * sqrt (M_PI / 2.0) * (erf ((GAUSS_XF - GAUSS_MU) / a) - erf (-GAUSS_MU / a));

  ncm_assert_cmpdouble_e (ncm_spline_eval_deriv (s, x0), ==, d1, 1.0e-11, 0.0);
  ncm_assert_cmpdouble_e (ncm_spline_eval_deriv2 (s, x0), ==, d2, 1.0e-9, 0.0);
  ncm_assert_cmpdouble_e (ncm_spline_eval_integ (s, 0.0, GAUSS_XF), ==, ex, 1.0e-12, 0.0);

  ncm_spline_free (s);
}

/* Outside the knots the integral is that of the extrapolated edge polynomial, like the
 * value. Polynomials of degree order - 1 are reproduced, so the integrals are exact up to
 * rounding (measured 6.5e-13 at most). */
static void
test_ncm_spline_bspline_integ_extrap (void)
{
  const gdouble lims[][2] = {
    {
      -0.5, 1.5
    }, {
      1.2, 1.7
    }, {
      -1.0, -0.5
    }, {
      0.3, 0.7
    }, {
      1.5, -0.5
    }
  };
  guint order;

  for (order = 2; order <= 4; order++)
  {
    NcmVector *xv = ncm_vector_new (50);
    NcmVector *yv = ncm_vector_new (50);
    NcmSpline *s;
    guint i;

    for (i = 0; i < 50; i++)
    {
      const gdouble x = i / 49.0;

      ncm_vector_set (xv, i, x);
      ncm_vector_set (yv, i, gsl_pow_uint (x, order - 1) - 2.0 * x);
    }

    s = NCM_SPLINE (ncm_spline_bspline_new_full (order, xv, yv, TRUE));

    for (i = 0; i < G_N_ELEMENTS (lims); i++)
    {
      const gdouble a  = lims[i][0];
      const gdouble b  = lims[i][1];
      const gdouble ex = (gsl_pow_uint (b, order) - gsl_pow_uint (a, order)) / order - (b * b - a * a);

      ncm_assert_cmpdouble_e (ncm_spline_eval_integ (s, a, b), ==, ex, 0.0, 1.0e-11);
    }

    ncm_spline_free (s);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }
}

/* The top derivative is the degree's: (order - 1)! for x^(order - 1), and at orders 2
 * and 3 the first and second derivatives. The tolerance grows with the order: measured
 * 3e-10 at order 4 and 4e-6 at order 6 on 200 samples. */
static void
test_ncm_spline_bspline_deriv_nmax (void)
{
  const struct
  {
    guint order;
    gdouble reltol;
  } cases[] = {
    {4, 1.0e-8},
    {6, 1.0e-4},
  };

  const gdouble xs[] = {0.3, 0.9, 1.5};
  guint c, j;

  for (c = 0; c < G_N_ELEMENTS (cases); c++)
  {
    const guint deg = cases[c].order - 1;
    NcmVector *xv   = ncm_vector_new (200);
    NcmVector *yv   = ncm_vector_new (200);
    NcmSpline *s;
    guint i;

    for (i = 0; i < 200; i++)
    {
      const gdouble x = 2.0 * i / 199.0;

      ncm_vector_set (xv, i, x);
      ncm_vector_set (yv, i, gsl_pow_uint (x, deg));
    }

    s = NCM_SPLINE (ncm_spline_bspline_new_full (cases[c].order, xv, yv, TRUE));

    for (j = 0; j < G_N_ELEMENTS (xs); j++)
      ncm_assert_cmpdouble_e (ncm_spline_eval_deriv_nmax (s, xs[j]), ==, gsl_sf_fact (deg), cases[c].reltol, 0.0);

    ncm_spline_free (s);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }

  {
    NcmSpline *s2 = _new_gauss (2, 200);
    NcmSpline *s3 = _new_gauss (3, 200);

    for (j = 0; j < G_N_ELEMENTS (xs); j++)
    {
      const gdouble x = 1000.0 * xs[j];

      g_assert_cmpfloat (ncm_spline_eval_deriv_nmax (s2, x), ==, ncm_spline_eval_deriv (s2, x));
      g_assert_cmpfloat (ncm_spline_eval_deriv_nmax (s3, x), ==, ncm_spline_eval_deriv2 (s3, x));
    }

    ncm_spline_free (s2);
    ncm_spline_free (s3);
  }
}

/* A change of order rebuilds the workspace before the next preparation. */
static void
test_ncm_spline_bspline_reprepare (void)
{
  NcmSpline *s = _new_gauss (4, 500);
  gdouble err4, err8;

  err4 = _max_midpoint_error (s);
  ncm_spline_bspline_set_order (NCM_SPLINE_BSPLINE (s), 8);
  ncm_spline_prepare (s);
  err8 = _max_midpoint_error (s);

  g_assert_cmpfloat (err8, <, err4 / 100.0);

  ncm_spline_free (s);
}

static void
test_ncm_spline_bspline_serialize (void)
{
  NcmSpline *s      = _new_gauss (6, 200);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpline *dup    = NCM_SPLINE (ncm_serialize_dup_obj (ser, G_OBJECT (s)));

  ncm_spline_prepare (dup);

  g_assert_true (NCM_IS_SPLINE_BSPLINE (dup));
  g_assert_cmpuint (ncm_spline_bspline_get_order (NCM_SPLINE_BSPLINE (dup)), ==, 6);
  ncm_assert_cmpdouble_e (ncm_spline_eval (dup, 1234.5), ==, ncm_spline_eval (s, 1234.5), 1.0e-14, 0.0);

  ncm_spline_free (dup);
  ncm_serialize_free (ser);
  ncm_spline_free (s);
}

static void
test_ncm_spline_bspline_copy_empty (void)
{
  NcmSpline *s     = _new_gauss (6, 200);
  NcmSpline *empty = ncm_spline_copy_empty (s);

  g_assert_true (NCM_IS_SPLINE_BSPLINE (empty));
  g_assert_cmpuint (ncm_spline_bspline_get_order (NCM_SPLINE_BSPLINE (empty)), ==, 6);
  g_assert_cmpuint (ncm_spline_get_len (empty), ==, 0);

  /* Filling the copy leaves the original unchanged */
  _fill_gauss (empty, 50, TRUE);
  ncm_assert_cmpdouble_e (ncm_spline_eval (s, 700.0), ==, _gauss (700.0), 1.0e-9, 0.0);

  ncm_spline_free (empty);
  ncm_spline_free (s);

  /* The tolerances are copied, so the copy selects its own order */
  {
    NcmSpline *s_tol = NCM_SPLINE (ncm_spline_bspline_new_tol (1.0e-8, 1.0e-12));
    NcmSpline *e_tol = ncm_spline_copy_empty (s_tol);
    gdouble reltol, abstol;

    g_object_get (e_tol, "reltol", &reltol, "abstol", &abstol, NULL);
    g_assert_cmpfloat (reltol, ==, 1.0e-8);
    g_assert_cmpfloat (abstol, ==, 1.0e-12);

    _fill_gauss (e_tol, 100, TRUE);
    g_assert_cmpuint (ncm_spline_bspline_get_order (NCM_SPLINE_BSPLINE (e_tol)), ==, 6);

    ncm_spline_free (e_tol);
    ncm_spline_free (s_tol);
  }
}

/* The chosen order is the lowest even one meeting the tolerance, and the estimated error
 * agrees with the true one within a factor of two. */
static void
test_ncm_spline_bspline_tol (void)
{
  const struct
  {
    guint n;
    gdouble reltol;
    guint order;
  } cases[] = {
    { 100, 1.0e-06, 4},
    { 100, 1.0e-08, 6},
    { 100, 1.0e-10, 8},
    { 500, 1.0e-06, 4},
    { 500, 1.0e-10, 6},
    {2000, 1.0e-06, 4},
    {2000, 1.0e-14, 6},
  };

  guint c;

  for (c = 0; c < G_N_ELEMENTS (cases); c++)
  {
    NcmSplineBSpline *sbs = ncm_spline_bspline_new_tol (cases[c].reltol, 0.0);
    NcmSpline *s          = NCM_SPLINE (sbs);
    gdouble err;

    _fill_gauss (s, cases[c].n, TRUE);
    err = _max_midpoint_error (s);

    g_assert_cmpuint (ncm_spline_bspline_get_order (sbs), ==, cases[c].order);
    g_assert_cmpfloat (err, <=, cases[c].reltol);
    ncm_assert_cmpdouble_e (ncm_spline_bspline_get_achieved_error (sbs), ==, err, 0.5, 0.0);

    ncm_spline_free (s);
  }

  /* reltol = 0 keeps the order set and reports no estimate */
  {
    NcmSpline *s = _new_gauss (10, 2000);

    g_assert_cmpuint (ncm_spline_bspline_get_order (NCM_SPLINE_BSPLINE (s)), ==, 10);
    g_assert_cmpfloat (ncm_spline_bspline_get_achieved_error (NCM_SPLINE_BSPLINE (s)), ==, 0.0);

    ncm_spline_free (s);
  }
}

/* 100 samples reach about 3e-11 at the highest order, so 1e-13 aborts. */
static void
test_ncm_spline_bspline_tol_unreachable (void)
{
  g_test_trap_subprocess ("/ncm/spline_bspline/tol_unreachable/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*cannot support a requested interpolation error*");
}

static void
test_ncm_spline_bspline_tol_unreachable_subprocess (void)
{
  NcmSpline *s = NCM_SPLINE (ncm_spline_bspline_new_tol (1.0e-13, 0.0));

  _fill_gauss (s, 100, TRUE);
}

/* The min-size error names the spline by its order. */
static void
test_ncm_spline_bspline_min_size (void)
{
  g_test_trap_subprocess ("/ncm/spline_bspline/min_size/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*min size for [NcmSplineBSpline[order 8]] is 8*");
}

static void
test_ncm_spline_bspline_min_size_subprocess (void)
{
  NcmSpline *s = NCM_SPLINE (ncm_spline_bspline_new (8));

  _fill_gauss (s, 5, TRUE);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/spline_bspline/order", &test_ncm_spline_bspline_order);
  g_test_add_func ("/ncm/spline_bspline/interpolates", &test_ncm_spline_bspline_interpolates);
  g_test_add_func ("/ncm/spline_bspline/accuracy", &test_ncm_spline_bspline_accuracy);
  g_test_add_func ("/ncm/spline_bspline/deriv_integ", &test_ncm_spline_bspline_deriv_integ);
  g_test_add_func ("/ncm/spline_bspline/integ_extrap", &test_ncm_spline_bspline_integ_extrap);
  g_test_add_func ("/ncm/spline_bspline/deriv_nmax", &test_ncm_spline_bspline_deriv_nmax);
  g_test_add_func ("/ncm/spline_bspline/reprepare", &test_ncm_spline_bspline_reprepare);
  g_test_add_func ("/ncm/spline_bspline/serialize", &test_ncm_spline_bspline_serialize);
  g_test_add_func ("/ncm/spline_bspline/copy_empty", &test_ncm_spline_bspline_copy_empty);
  g_test_add_func ("/ncm/spline_bspline/tol", &test_ncm_spline_bspline_tol);
  g_test_add_func ("/ncm/spline_bspline/tol_unreachable", &test_ncm_spline_bspline_tol_unreachable);
  g_test_add_func ("/ncm/spline_bspline/tol_unreachable/subprocess", &test_ncm_spline_bspline_tol_unreachable_subprocess);
  g_test_add_func ("/ncm/spline_bspline/min_size", &test_ncm_spline_bspline_min_size);
  g_test_add_func ("/ncm/spline_bspline/min_size/subprocess", &test_ncm_spline_bspline_min_size_subprocess);

  g_test_run ();
}


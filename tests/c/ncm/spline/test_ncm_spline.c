/***************************************************************************
 *            test_ncm_spline.c
 *
 *  Mon April 23 23:15:46 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
#include <gsl/gsl_errno.h>

#include <glib.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>

typedef struct _TestNcmSpline
{
  NcmSpline *s_base;
  gsl_function F;
  guint nknots;
  gdouble xi;
  gdouble dx;
  gdouble prec;
  gdouble error;
  gdouble error_d1;
  gdouble error_d2;
  const gchar *name;
  gboolean deriv1;
  gboolean deriv2;
} TestNcmSpline;

gdouble
F_linear (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[0] + d[1] * x;
}

gdouble
F_linear_deriv (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[1];
}

gdouble
F_linear_deriv2 (gdouble x, gpointer p)
{
  return 0.0;
}

gdouble
F_linear_int_0x (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[0] * x + d[1] * x * x * 0.5;
}

gdouble
F_cubic (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[0] + (d[1] + (d[2] + d[3] * x) * x) * x;
}

gdouble
F_cubic_deriv (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[1] + (2.0 * d[2] + 3.0 * d[3] * x) * x;
}

gdouble
F_cubic_deriv2 (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return 2.0 * d[2] + 6.0 * d[3] * x;
}

gdouble
F_cubic_int_0x (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return (d[0] + (0.5 * d[1] + (d[2] / 3.0 + 0.25 * d[3] * x) * x) * x) * x;
}

gdouble
F_sin_poly (gdouble x, gpointer p)
{
  return sin (M_PI * x) * exp (x / 5.0) +
         (x * x / 2.0 + 5.0 * x + 0.1 * x * x * x + 1.0e-3 * gsl_pow_5 (x));
}

gdouble
F_sin_poly_deriv (gdouble x, gpointer p)
{
  return (M_PI * cos (M_PI * x) + sin (M_PI * x) / 5.0) * exp (x / 5.0) +
         (5.0 + x + 0.3 * x * x + 5.0e-3 * gsl_pow_4 (x));
}

gdouble
F_sin_poly_deriv2 (gdouble x, gpointer p)
{
  return (2.0 * M_PI * cos (M_PI * x) / 5.0 + (1.0 / 25.0 - M_PI * M_PI) * sin (M_PI * x)) * exp (x / 5.0) +
         (1.0 + 0.6 * x + 2.0e-2 * gsl_pow_3 (x));
}

void test_ncm_spline_cubic_notaknot_new_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_gsl_cspline_new_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_gsl_linear_new_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_gsl_steffen_new_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_new (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_new_array (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_new_data (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_copy_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_copy (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_set_type (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_set_type_serialize (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_serialize (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_idx (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_deriv (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_deriv_idx (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_deriv2 (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_int (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_integ_idx (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_get_index (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_get_index_stride2 (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_get_index_stride5 (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_get_index_acc (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_get_index_acc_stride2 (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_get_index_acc_stride5 (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_free_empty (TestNcmSpline *test, gconstpointer pdata);

void test_ncm_spline_invalid_vector_sizes (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_min_vector_sizes (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_x_vector (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_y_vector (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_xy_vector (TestNcmSpline *test, gconstpointer pdata);

void test_ncm_spline_invalid_array_sizes (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_min_array_sizes (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_x_array (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_y_array (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_invalid_xy_array (TestNcmSpline *test, gconstpointer pdata);

void test_ncm_spline_traps (TestNcmSpline *test, gconstpointer pdata);

typedef struct _TestNcmSplineFunc
{
  void (*func) (TestNcmSpline *test, gconstpointer pdata);

  const gchar *name;
} TestNcmSplineFunc;

TestNcmSplineFunc _test_ncm_spline_traps[] = {
  {&test_ncm_spline_invalid_vector_sizes,     "/vector/invalid/sizes/subprocess"},
  {&test_ncm_spline_invalid_min_vector_sizes, "/vector/invalid/min_sizes/subprocess"},
  {&test_ncm_spline_invalid_x_vector,         "/vector/invalid/x/subprocess"},
  {&test_ncm_spline_invalid_y_vector,         "/vector/invalid/y/subprocess"},
  {&test_ncm_spline_invalid_xy_vector,        "/vector/invalid/xy/subprocess"},
  {&test_ncm_spline_invalid_array_sizes,      "/array/invalid/sizes/subprocess"},
  {&test_ncm_spline_invalid_min_array_sizes,  "/array/invalid/min_sizes/subprocess"},
  {&test_ncm_spline_invalid_x_array,          "/array/invalid/x/subprocess"},
  {&test_ncm_spline_invalid_y_array,          "/array/invalid/y/subprocess"},
  {&test_ncm_spline_invalid_xy_array,         "/array/invalid/xy/subprocess"},
  {NULL, NULL}
};

TestNcmSplineFunc _test_ncm_spline_tests[] = {
  {&test_ncm_spline_new,                   "/new"},
  {&test_ncm_spline_new_array,             "/new/array"},
  {&test_ncm_spline_new_data,              "/new/data"},
  {&test_ncm_spline_copy_empty,            "/copy/empty"},
  {&test_ncm_spline_copy,                  "/copy"},
  {&test_ncm_spline_set_type,              "/set_type"},
  {&test_ncm_spline_set_type_serialize,    "/set_type/serialize"},
  {&test_ncm_spline_serialize,             "/serialize"},
  {&test_ncm_spline_eval,                  "/eval"},
  {&test_ncm_spline_eval_idx,              "/eval_idx"},
  {&test_ncm_spline_eval_deriv,            "/eval/deriv"},
  {&test_ncm_spline_eval_deriv_idx,        "/eval/deriv_idx"},
  {&test_ncm_spline_eval_deriv2,           "/eval/deriv2"},
  {&test_ncm_spline_eval_int,              "/int"},
  {&test_ncm_spline_eval_integ_idx,        "/int_idx"},
  {&test_ncm_spline_get_index,             "/get_index"},
  {&test_ncm_spline_get_index_stride2,     "/get_index/stride2"},
  {&test_ncm_spline_get_index_stride5,     "/get_index/stride5"},
  {&test_ncm_spline_get_index_acc,         "/get_index/acc"},
  {&test_ncm_spline_get_index_acc_stride2, "/get_index/acc/stride2"},
  {&test_ncm_spline_get_index_acc_stride5, "/get_index/acc/stride5"},
  {&test_ncm_spline_traps,                 "/traps"},
  {NULL, NULL}
};

void
_test_ncm_spline_add_tests (void ( *tnew ) (TestNcmSpline *test, gconstpointer pdata), void (*tfree) (TestNcmSpline *test, gconstpointer pdata), const gchar *name)
{
  guint i;

  for (i = 0; _test_ncm_spline_traps[i].func != NULL; i++)
  {
    gchar *tname = g_strdup_printf ("/ncm/%s%s", name, _test_ncm_spline_traps[i].name);

    g_test_add (tname, TestNcmSpline, NULL, tnew, _test_ncm_spline_traps[i].func, tfree);
    g_free (tname);
  }

  for (i = 0; _test_ncm_spline_tests[i].func != NULL; i++)
  {
    gchar *tname = g_strdup_printf ("/ncm/%s%s", name, _test_ncm_spline_tests[i].name);

    g_test_add (tname, TestNcmSpline, NULL, tnew, _test_ncm_spline_tests[i].func, tfree);
    g_free (tname);
  }
}

/* Not-a-knot reproduces a cubic exactly, also where it extrapolates the boundary cubics */
static void
test_ncm_spline_cubic_notaknot_extrapolates (void)
{
  NcmVector *xv = ncm_vector_new (11);
  NcmVector *yv = ncm_vector_new (11);
  NcmSpline *s;
  const gdouble xs[] = {-0.3, 0.37, 1.2};
  guint i;

  for (i = 0; i < 11; i++)
  {
    const gdouble x = 0.1 * i;

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, x * x * x - 2.0 * x);
  }

  s = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  for (i = 0; i < G_N_ELEMENTS (xs); i++)
  {
    const gdouble x = xs[i];

    g_assert_cmpfloat (fabs (ncm_spline_eval (s, x) - (x * x * x - 2.0 * x)), <, 2.0e-14);
    g_assert_cmpfloat (fabs (ncm_spline_eval_deriv (s, x) - (3.0 * x * x - 2.0)), <, 2.0e-13);
  }

  ncm_spline_free (s);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

/* With the GSL error handler off, NcmSplineGsl returns NaN outside the knots */
static void
test_ncm_spline_gsl_domain_nan (void)
{
  NcmVector *xv = ncm_vector_new (11);
  NcmVector *yv = ncm_vector_new (11);
  NcmSpline *s;
  guint i;

  for (i = 0; i < 11; i++)
  {
    ncm_vector_set (xv, i, 0.1 * i);
    ncm_vector_set (yv, i, 0.01 * i * i);
  }

  s = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_CSPLINE, xv, yv, TRUE));

  gsl_set_error_handler_off ();
  g_assert_true (isnan (ncm_spline_eval (s, 1.05)));
  g_assert_true (isnan (ncm_spline_eval (s, -0.05)));
  g_assert_true (isnan (ncm_spline_eval_integ (s, -0.1, 0.5)));
  ncm_cfg_enable_gsl_err_handler ();

  ncm_spline_free (s);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

/* Reversed limits give minus the integral, for every spline type */
static void
test_ncm_spline_integ_reversed (void)
{
  NcmVector *xv = ncm_vector_new (11);
  NcmVector *yv = ncm_vector_new (11);
  NcmSpline *splines[2];
  guint i;

  for (i = 0; i < 11; i++)
  {
    ncm_vector_set (xv, i, 0.1 * i);
    ncm_vector_set (yv, i, 0.01 * i * i);
  }

  splines[0] = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));
  splines[1] = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_CSPLINE, xv, yv, TRUE));

  for (i = 0; i < G_N_ELEMENTS (splines); i++)
  {
    const gdouble fwd = ncm_spline_eval_integ (splines[i], 0.15, 0.85);

    g_assert_cmpfloat (ncm_spline_eval_integ (splines[i], 0.85, 0.15), ==, -fwd);
    g_assert_cmpfloat (ncm_spline_eval_integ (splines[i], 0.18, 0.12), ==, -ncm_spline_eval_integ (splines[i], 0.12, 0.18));
  }

  {
    const gsize ia = ncm_spline_get_index (splines[0], 0.15);
    const gsize ib = ncm_spline_get_index (splines[0], 0.85);

    g_assert_cmpfloat (ncm_spline_eval_integ_idx (splines[0], 0.85, ib, 0.15, ia), ==,
                       -ncm_spline_eval_integ_idx (splines[0], 0.15, ia, 0.85, ib));
  }

  ncm_spline_free (splines[0]);
  ncm_spline_free (splines[1]);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

static NcmSpline *
_test_ncm_spline_cube (gdouble scale)
{
  NcmVector *xv = ncm_vector_new (11);
  NcmVector *yv = ncm_vector_new (11);
  NcmSpline *s;
  guint i;

  for (i = 0; i < 11; i++)
  {
    const gdouble x = 0.1 * i;

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, scale * x * x * x);
  }

  s = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));
  ncm_vector_free (xv);
  ncm_vector_free (yv);

  return s;
}

/* Not-a-knot reproduces x^3, so c(x) = 6x and N_p = (6^p / (p + 1))^(1/p) on [0, 1], at any scale */
static void
test_ncm_spline_curvature_norms (void)
{
  const gdouble ps[]    = {0.5, 1.0, 2.0};
  const gdouble scale[] = {1.0, 1.0e-9};
  guint i, j;

  for (j = 0; j < G_N_ELEMENTS (scale); j++)
  {
    NcmSpline *s = _test_ncm_spline_cube (scale[j]);
    NcmSpline *w = _test_ncm_spline_cube (0.0);

    ncm_vector_set_all (ncm_spline_peek_yv (w), 2.0);
    ncm_spline_prepare (w);

    for (i = 0; i < G_N_ELEMENTS (ps); i++)
    {
      const gdouble p     = ps[i];
      const gdouble exact = scale[j] * pow (pow (6.0, p) / (p + 1.0), 1.0 / p);
      const gdouble n_p   = ncm_spline_curvature_lp_norm (s, NCM_SPLINE_CURVATURE_D2, p, 0.0, 1.0);
      const gdouble n_w   = ncm_spline_curvature_weighted_lp_norm (s, NCM_SPLINE_CURVATURE_D2, p, w, 0.0, 1.0);

      g_assert_cmpfloat (fabs (n_p / exact - 1.0), <, 1.0e-13);
      g_assert_cmpfloat (fabs (n_w / exact - 1.0), <, 1.0e-13);
    }

    ncm_spline_free (s);
    ncm_spline_free (w);
  }
}

static void
test_ncm_spline_notaknot_repeated_knot (void)
{
  g_test_trap_subprocess ("/ncm/spline_cubic_notaknot/repeated_knot/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*strictly increasing*");
}

static void
test_ncm_spline_notaknot_repeated_knot_subprocess (void)
{
  NcmVector *xv = ncm_vector_new (8);
  NcmVector *yv = ncm_vector_new (8);
  guint i;

  for (i = 0; i < 8; i++)
  {
    ncm_vector_set (xv, i, (i == 0) ? 0.0 : 0.1 * (i - 1));
    ncm_vector_set (yv, i, 1.0 * i);
  }

  ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE);
}

/* After ncm_spline_cubic_d2_set_d2 the spline reproduces x^3 from its values and second derivatives */
static void
test_ncm_spline_cubic_d2_set_d2 (void)
{
  NcmVector *xv  = ncm_vector_new (11);
  NcmVector *yv  = ncm_vector_new (11);
  NcmVector *d2v = ncm_vector_new (11);
  NcmSplineCubicD2 *scd2;
  guint i;

  for (i = 0; i < 11; i++)
  {
    const gdouble x = 0.1 * i;

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, x * x);
    ncm_vector_set (d2v, i, 2.0);
  }

  scd2 = ncm_spline_cubic_d2_new (xv, yv, d2v, TRUE);

  for (i = 0; i < 11; i++)
  {
    const gdouble x = 0.1 * i;

    ncm_vector_set (yv, i, x * x * x);
    ncm_vector_set (d2v, i, 6.0 * x);
  }

  ncm_spline_cubic_d2_set_d2 (scd2, d2v, TRUE);

  for (i = 0; i < 10; i++)
  {
    const gdouble x = 0.1 * i + 0.05;

    g_assert_cmpfloat (fabs (ncm_spline_eval (NCM_SPLINE (scd2), x) - x * x * x), <, 1.0e-15);
  }

  ncm_spline_free (NCM_SPLINE (scd2));
  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_vector_free (d2v);
}

static void
test_ncm_spline_invalid_d2_len (void)
{
  g_test_trap_subprocess ("/ncm/spline_cubic_d2/invalid_d2_len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

static void
test_ncm_spline_invalid_d2_len_subprocess (void)
{
  NcmVector *xv  = ncm_vector_new (5);
  NcmVector *yv  = ncm_vector_new (5);
  NcmVector *d2v = ncm_vector_new (4);
  NcmSplineCubicD2 *scd2;

  ncm_vector_set_all (yv, 1.0);
  ncm_vector_set_all (d2v, 0.0);
  ncm_vector_set (xv, 0, 0.0);
  ncm_vector_set (xv, 1, 1.0);
  ncm_vector_set (xv, 2, 2.0);
  ncm_vector_set (xv, 3, 3.0);
  ncm_vector_set (xv, 4, 4.0);

  scd2 = ncm_spline_cubic_d2_new (xv, yv, yv, FALSE);
  ncm_spline_cubic_d2_set_d2 (scd2, d2v, FALSE);
}

static void
test_ncm_spline_gsl_unsupported_type (void)
{
  g_test_trap_subprocess ("/ncm/spline_gsl/unsupported_type/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not one of NcmSplineGslType*");
}

static void
test_ncm_spline_gsl_unsupported_type_subprocess (void)
{
  static const gsl_interp_type not_a_type = {"not-a-type", 2, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

  ncm_spline_gsl_new (&not_a_type);
}

/* Steffen stays between neighbouring knot values on monotone data; a natural cubic spline overshoots a step */
static void
test_ncm_spline_gsl_steffen_monotone (void)
{
  const guint nknots = 40;
  NcmVector *xv      = ncm_vector_new (nknots);
  NcmVector *yv      = ncm_vector_new (nknots);
  NcmSpline *steffen;
  NcmSpline *cspline;
  gdouble cspline_overshoot = 0.0;
  guint i, j;

  for (i = 0; i < nknots; i++)
  {
    ncm_vector_set (xv, i, i);
    ncm_vector_set (yv, i, (i < nknots / 2) ? 1.0e-6 * i : 1.0 + 1.0e-6 * i);
  }

  steffen = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_STEFFEN, xv, yv, TRUE));
  cspline = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_CSPLINE, xv, yv, TRUE));

  for (i = 0; i + 1 < nknots; i++)
  {
    const gdouble y_i   = ncm_vector_get (yv, i);
    const gdouble y_ip1 = ncm_vector_get (yv, i + 1);

    for (j = 0; j <= 10; j++)
    {
      const gdouble x = i + 0.1 * j;
      const gdouble s = ncm_spline_eval (steffen, x);
      const gdouble c = ncm_spline_eval (cspline, x);

      g_assert_cmpfloat (s, >=, y_i);
      g_assert_cmpfloat (s, <=, y_ip1);
      cspline_overshoot = GSL_MAX (cspline_overshoot, GSL_MAX (y_i - c, c - y_ip1));
    }
  }

  g_assert_cmpfloat (cspline_overshoot, >, 1.0e-2);
  g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (steffen)) == gsl_interp_steffen);

  {
    NcmSplineGsl *by_name = ncm_spline_gsl_new (gsl_interp_linear);

    ncm_spline_gsl_set_type_by_name (by_name, "steffen");
    g_assert_true (ncm_spline_gsl_get_type_id (by_name) == NCM_SPLINE_GSL_STEFFEN);
    ncm_spline_free (NCM_SPLINE (by_name));
  }

  ncm_spline_free (steffen);
  ncm_spline_free (cspline);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

/* The third derivative is constant on each interval, including for splines whose second derivative jumps at the knots */
static void
test_ncm_spline_gsl_deriv_nmax (void)
{
  const NcmSplineGslType types[] = {NCM_SPLINE_GSL_CSPLINE, NCM_SPLINE_GSL_AKIMA, NCM_SPLINE_GSL_STEFFEN};
  const guint nknots             = 21;
  NcmVector *xv                  = ncm_vector_new (nknots);
  NcmVector *yv                  = ncm_vector_new (nknots);
  guint t, i;

  for (i = 0; i < nknots; i++)
  {
    const gdouble x = 3.0 * i / (nknots - 1.0);

    ncm_vector_set (xv, i, x);
    ncm_vector_set (yv, i, sin (3.0 * x) + x * x * x);
  }

  for (t = 0; t < G_N_ELEMENTS (types); t++)
  {
    NcmSpline *s = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (types[t], xv, yv, TRUE));

    for (i = 0; i + 1 < nknots; i++)
    {
      const gdouble x_i = ncm_vector_get (xv, i);
      const gdouble dx  = ncm_vector_get (xv, i + 1) - x_i;
      const gdouble a   = x_i + 0.25 * dx;
      const gdouble b   = x_i + 0.75 * dx;
      const gdouble d3  = (ncm_spline_eval_deriv2 (s, b) - ncm_spline_eval_deriv2 (s, a)) / (b - a);

      ncm_assert_cmpdouble_e (ncm_spline_eval_deriv_nmax (s, x_i + 0.3 * dx), ==, d3, 1.0e-10, 1.0e-10);
    }

    ncm_spline_free (s);
  }

  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  _test_ncm_spline_add_tests (&test_ncm_spline_cubic_notaknot_new_empty,
                              &test_ncm_spline_free_empty,
                              "spline_cubic_notaknot");
  _test_ncm_spline_add_tests (&test_ncm_spline_gsl_cspline_new_empty,
                              &test_ncm_spline_free_empty,
                              "spline_gsl/cspline");
  _test_ncm_spline_add_tests (&test_ncm_spline_gsl_linear_new_empty,
                              &test_ncm_spline_free_empty,
                              "spline_gsl/linear");
  _test_ncm_spline_add_tests (&test_ncm_spline_gsl_steffen_new_empty,
                              &test_ncm_spline_free_empty,
                              "spline_gsl/steffen");
  g_test_add_func ("/ncm/spline_cubic_notaknot/extrapolates", &test_ncm_spline_cubic_notaknot_extrapolates);
  g_test_add_func ("/ncm/spline_gsl/domain_nan", &test_ncm_spline_gsl_domain_nan);
  g_test_add_func ("/ncm/spline/integ_reversed", &test_ncm_spline_integ_reversed);
  g_test_add_func ("/ncm/spline/curvature_norms", &test_ncm_spline_curvature_norms);
  g_test_add_func ("/ncm/spline_cubic_notaknot/repeated_knot", &test_ncm_spline_notaknot_repeated_knot);
  g_test_add_func ("/ncm/spline_cubic_notaknot/repeated_knot/subprocess", &test_ncm_spline_notaknot_repeated_knot_subprocess);
  g_test_add_func ("/ncm/spline_cubic_d2/set_d2", &test_ncm_spline_cubic_d2_set_d2);
  g_test_add_func ("/ncm/spline_cubic_d2/invalid_d2_len", &test_ncm_spline_invalid_d2_len);
  g_test_add_func ("/ncm/spline_cubic_d2/invalid_d2_len/subprocess", &test_ncm_spline_invalid_d2_len_subprocess);
  g_test_add_func ("/ncm/spline_gsl/unsupported_type", &test_ncm_spline_gsl_unsupported_type);
  g_test_add_func ("/ncm/spline_gsl/unsupported_type/subprocess", &test_ncm_spline_gsl_unsupported_type_subprocess);
  g_test_add_func ("/ncm/spline_gsl/steffen/monotone", &test_ncm_spline_gsl_steffen_monotone);
  g_test_add_func ("/ncm/spline_gsl/deriv_nmax", &test_ncm_spline_gsl_deriv_nmax);
  g_test_run ();
}

#define _NCM_SPLINE_TEST_NKNOTS 1000
#define _NCM_SPLINE_TEST_DX 0.003

void
test_ncm_spline_cubic_notaknot_new_empty (TestNcmSpline *test, gconstpointer pdata)
{
  test->name     = "spline_cubic_notaknot";
  test->deriv1   = TRUE;
  test->deriv2   = TRUE;
  test->nknots   = g_test_rand_int_range (_NCM_SPLINE_TEST_NKNOTS, 2 * _NCM_SPLINE_TEST_NKNOTS);
  test->dx       = _NCM_SPLINE_TEST_DX;
  test->xi       = 10.0 * GSL_SIGN (g_test_rand_double_range (-1.0, 1.0));
  test->prec     = 1.0e-5;
  test->error    = 1.0e-4;
  test->error_d1 = 5.0e-3;
  test->error_d2 = 1.0e-2;
  test->s_base   = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  g_assert_true (NCM_IS_SPLINE_CUBIC (test->s_base));
  g_assert_true (NCM_IS_SPLINE_CUBIC_NOTAKNOT (test->s_base));

  {
    NcmVector *xv = ncm_vector_new (test->nknots);
    NcmVector *yv = ncm_vector_new (test->nknots);
    guint i;

    for (i = 0; i < test->nknots; i++)
    {
      ncm_vector_set (xv, i, i * M_PI);
      ncm_vector_set (yv, i, i * M_PI_2);
    }

    ncm_spline_set (test->s_base, xv, yv, FALSE);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }
}

void
test_ncm_spline_gsl_cspline_new_empty (TestNcmSpline *test, gconstpointer pdata)
{
  test->name     = "spline_gsl/cspline";
  test->deriv1   = TRUE;
  test->deriv2   = FALSE;
  test->nknots   = g_test_rand_int_range (_NCM_SPLINE_TEST_NKNOTS, 2 * _NCM_SPLINE_TEST_NKNOTS);
  test->dx       = _NCM_SPLINE_TEST_DX;
  test->xi       = 10.0 * GSL_SIGN (g_test_rand_double_range (-1.0, 1.0));
  test->prec     = 1.0e-5;
  test->error    = 5.0e-3;
  test->error_d1 = 2.0e-1;
  test->error_d2 = 2.0e-1;
  test->s_base   = NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_cspline));
  g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
  {
    NcmVector *xv = ncm_vector_new (test->nknots);
    NcmVector *yv = ncm_vector_new (test->nknots);
    guint i;

    for (i = 0; i < test->nknots; i++)
    {
      ncm_vector_set (xv, i, i * M_PI);
      ncm_vector_set (yv, i, i * M_PI_2);
    }

    ncm_spline_set (test->s_base, xv, yv, FALSE);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }
}

void
test_ncm_spline_gsl_steffen_new_empty (TestNcmSpline *test, gconstpointer pdata)
{
  test->name     = "spline_gsl/steffen";
  test->deriv1   = TRUE;
  test->deriv2   = FALSE;
  test->nknots   = g_test_rand_int_range (_NCM_SPLINE_TEST_NKNOTS, 2 * _NCM_SPLINE_TEST_NKNOTS);
  test->dx       = _NCM_SPLINE_TEST_DX;
  test->xi       = 10.0 * GSL_SIGN (g_test_rand_double_range (-1.0, 1.0));
  test->prec     = 1.0e-5;
  test->error    = 5.0e-3;
  test->error_d1 = 2.0e-1;
  test->error_d2 = 2.0e-1;
  test->s_base   = NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_steffen));
  g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
  {
    NcmVector *xv = ncm_vector_new (test->nknots);
    NcmVector *yv = ncm_vector_new (test->nknots);
    guint i;

    for (i = 0; i < test->nknots; i++)
    {
      ncm_vector_set (xv, i, i * M_PI);
      ncm_vector_set (yv, i, i * M_PI_2);
    }

    ncm_spline_set (test->s_base, xv, yv, FALSE);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }
}

void
test_ncm_spline_gsl_linear_new_empty (TestNcmSpline *test, gconstpointer pdata)
{
  test->name     = "spline_gsl/linear";
  test->deriv1   = FALSE;
  test->deriv2   = FALSE;
  test->nknots   = g_test_rand_int_range (_NCM_SPLINE_TEST_NKNOTS, 2 * _NCM_SPLINE_TEST_NKNOTS);
  test->dx       = _NCM_SPLINE_TEST_DX;
  test->xi       = 10.0 * GSL_SIGN (g_test_rand_double_range (-1.0, 1.0));
  test->prec     = 1.0e-5;
  test->error    = 5.0e-3;
  test->error_d1 = 1.0;
  test->error_d2 = 1.0;
  test->s_base   = NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_linear));
  g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
  {
    NcmVector *xv = ncm_vector_new (test->nknots);
    NcmVector *yv = ncm_vector_new (test->nknots);
    guint i;

    for (i = 0; i < test->nknots; i++)
    {
      ncm_vector_set (xv, i, i * M_PI);
      ncm_vector_set (yv, i, i * M_PI_2);
    }

    ncm_spline_set (test->s_base, xv, yv, FALSE);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }
}

void
test_ncm_spline_free_empty (TestNcmSpline *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_spline_free, test->s_base);
}

void
test_ncm_spline_new_sanity (NcmSpline *s)
{
  g_assert_true (NCM_IS_SPLINE (s));
}

void
test_ncm_spline_new (TestNcmSpline *test, gconstpointer pdata)
{
  {
    NcmVector *x = ncm_vector_new (test->nknots);
    NcmVector *y = ncm_vector_new (test->nknots);
    NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert_true (!ncm_spline_is_init (s));

    ncm_spline_free (s);
    ncm_vector_free (x);
    ncm_vector_free (y);
  }

  {
    NcmVector *x = ncm_vector_new (test->nknots);
    NcmVector *y = ncm_vector_new (test->nknots);
    NcmSpline *s;
    gdouble d[2];
    guint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();

    for (i = 0; i < test->nknots; i++)
    {
      ncm_vector_set (x, i, test->xi + test->dx * i);
      ncm_vector_set (y, i, F_linear (ncm_vector_get (x, i), d));
    }

    s = ncm_spline_new (test->s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert_true (ncm_spline_is_init (s));

    ncm_spline_free (s);
    ncm_vector_free (x);
    ncm_vector_free (y);
  }
}

void
test_ncm_spline_new_array (TestNcmSpline *test, gconstpointer pdata)
{
  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);

    g_array_set_size (x, test->nknots);
    g_array_set_size (y, test->nknots);

    NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert_true (!ncm_spline_is_init (s));

    ncm_spline_free (s);
    g_array_unref (x);
    g_array_unref (y);
  }

  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
    NcmSpline *s;
    gdouble d[2];
    guint i;

    g_array_set_size (x, test->nknots);
    g_array_set_size (y, test->nknots);

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();

    for (i = 0; i < test->nknots; i++)
    {
      g_array_index (x, gdouble, i) = test->xi + test->dx * i;
      g_array_index (y, gdouble, i) = F_linear (g_array_index (x, gdouble, i), d);
    }

    s = ncm_spline_new_array (test->s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert_true (ncm_spline_is_init (s));

    ncm_spline_free (s);
    g_array_unref (x);
    g_array_unref (y);
  }
}

void
test_ncm_spline_new_data (TestNcmSpline *test, gconstpointer pdata)
{
  {
    gdouble x[test->nknots];
    gdouble y[test->nknots];
    NcmSpline *s = ncm_spline_new_data (test->s_base, x, y, test->nknots, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert_true (!ncm_spline_is_init (s));
    ncm_spline_free (s);
  }

  {
    gdouble x[test->nknots];
    gdouble y[test->nknots];
    NcmSpline *s;
    gdouble d[2];
    guint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();

    for (i = 0; i < test->nknots; i++)
    {
      x[i] = test->xi + test->dx * i;
      y[i] = F_linear (x[i], d);
    }

    s = ncm_spline_new_data (test->s_base, x, y, test->nknots, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert_true (ncm_spline_is_init (s));
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_copy_empty (TestNcmSpline *test, gconstpointer pdata)
{
  NcmSpline *s = ncm_spline_copy_empty (test->s_base);

  g_assert_true (G_TYPE_FROM_INSTANCE (s) == G_TYPE_FROM_INSTANCE (test->s_base));
  ncm_spline_free (s);
}

void
test_ncm_spline_copy (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *xv = ncm_vector_new (test->nknots);
  NcmVector *yv = ncm_vector_new (test->nknots);
  guint i;

  for (i = 0; i < test->nknots; i++)
  {
    ncm_vector_set (xv, i, i * M_PI);
    ncm_vector_set (yv, i, i * M_PI_2);
  }

  ncm_spline_set (test->s_base, xv, yv, FALSE);
  {
    NcmSpline *s         = ncm_spline_copy (test->s_base);
    NcmVector *s_xv      = ncm_spline_peek_xv (s);
    NcmVector *s_yv      = ncm_spline_peek_yv (s);
    NcmVector *s_base_xv = ncm_spline_peek_xv (test->s_base);
    NcmVector *s_base_yv = ncm_spline_peek_yv (test->s_base);

    g_assert_true (s_xv != s_base_xv && s_yv != s_base_yv);

    for (i = 0; i < test->nknots; i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (s_xv, i), ==, ncm_vector_get (s_base_xv, i));
      ncm_assert_cmpdouble (ncm_vector_get (s_yv, i), ==, ncm_vector_get (s_base_yv, i));
    }

    ncm_spline_free (s);
  }

  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

void
test_ncm_spline_serialize (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *xv = ncm_vector_new (test->nknots);
  NcmVector *yv = ncm_vector_new (test->nknots);
  guint i;

  for (i = 0; i < test->nknots; i++)
  {
    ncm_vector_set (xv, i, i * M_PI);
    ncm_vector_set (yv, i, i * M_PI_2);
  }

  ncm_spline_set (test->s_base, xv, yv, FALSE);
  {
    NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcmSpline *s         = NCM_SPLINE (ncm_serialize_dup_obj (ser, G_OBJECT (test->s_base)));
    NcmVector *s_xv      = ncm_spline_peek_xv (s);
    NcmVector *s_yv      = ncm_spline_peek_yv (s);
    NcmVector *s_base_xv = ncm_spline_peek_xv (test->s_base);
    NcmVector *s_base_yv = ncm_spline_peek_yv (test->s_base);

    ncm_spline_prepare (s);
    ncm_serialize_free (ser);

    g_assert_true (s_xv != s_base_xv && s_yv != s_base_yv);

    for (i = 0; i < test->nknots; i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (s_xv, i), ==, ncm_vector_get (s_base_xv, i));
      ncm_assert_cmpdouble (ncm_vector_get (s_yv, i), ==, ncm_vector_get (s_base_yv, i));
    }

    ncm_spline_free (s);
  }

  ncm_vector_free (xv);
  ncm_vector_free (yv);
}

void
test_ncm_spline_set_type (TestNcmSpline *test, gconstpointer pdata)
{
  if (!NCM_IS_SPLINE_GSL (test->s_base))
  {
    return;
  }
  else
  {
    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_linear);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_LINEAR);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_linear);

    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_polynomial);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_POLYNOMIAL);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_polynomial);

    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_cspline);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_CSPLINE);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_cspline);

    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_akima);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_AKIMA);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_akima);

    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_akima_periodic);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_AKIMA_PERIODIC);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_akima_periodic);

    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_cspline_periodic);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_CSPLINE_PERIODIC);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_cspline_periodic);

    ncm_spline_gsl_set_type (NCM_SPLINE_GSL (test->s_base), gsl_interp_steffen);
    g_assert_true (NCM_IS_SPLINE_GSL (test->s_base));
    g_assert_true (ncm_spline_gsl_get_type_id (NCM_SPLINE_GSL (test->s_base)) == NCM_SPLINE_GSL_STEFFEN);
    g_assert_true (ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base)) == gsl_interp_steffen);
  }
}

void
test_ncm_spline_set_type_serialize (TestNcmSpline *test, gconstpointer pdata)
{
  if (!NCM_IS_SPLINE_GSL (test->s_base))
  {
    return;
  }
  else
  {
    const gsl_interp_type *type = ncm_spline_gsl_get_gsl_type (NCM_SPLINE_GSL (test->s_base));
    NcmSerialize *ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcmSplineGsl *s             = NCM_SPLINE_GSL (ncm_serialize_dup_obj (ser, G_OBJECT (test->s_base)));

    g_assert_true (NCM_IS_SPLINE_GSL (s));
    g_assert_true (ncm_spline_gsl_get_gsl_type (s) == type);
    ncm_serialize_free (ser);
    NCM_TEST_FREE (ncm_spline_free, NCM_SPLINE (s));
  }
}

#define _TEST_EPSILON (1.00000001)

void
test_ncm_spline_eval (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  {
    NcmVector *x = ncm_vector_new (test->nknots);
    NcmVector *y = ncm_vector_new (test->nknots);
    NcmSpline *s;
    guint i;

    F.function = &F_linear;
    F.params   = d;

    ncm_vector_set (x, 0, test->xi);
    ncm_vector_set (y, 0, F_linear (ncm_vector_get (x, 0), d));

    ncm_vector_set (x, 1, test->xi + 2.0 * test->dx);
    ncm_vector_set (y, 1, F_linear (ncm_vector_get (x, 1), d));

    ncm_vector_set (x, 2, test->xi + 3.0 * test->dx);
    ncm_vector_set (y, 2, F_linear (ncm_vector_get (x, 2), d));

    for (i = 3; i < test->nknots - 1; i++)
    {
      ncm_vector_set (x, i, test->xi + test->dx * (i + 1.0));
      ncm_vector_set (y, i, F_linear (ncm_vector_get (x, i), d));
    }

    ncm_vector_set (x, test->nknots - 1, ncm_vector_get (x, test->nknots - 2) + 2.0 * test->dx * _TEST_EPSILON);
    ncm_vector_set (y, test->nknots - 1, F_linear (ncm_vector_get (x, test->nknots - 1), d));

    s = ncm_spline_new (test->s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert_true (ncm_spline_is_init (s) == TRUE);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      const gdouble xi = ncm_vector_get (x, 0);
      const gdouble xf = ncm_vector_get (x, test->nknots - 1);
      gdouble xval0    = xi + (xf - xi) / (2.0 * test->nknots - 1.0) * i;
      gdouble xval1    = GSL_MAX (xval0, xi);
      gdouble xval     = GSL_MIN (xval1, xf);
      gdouble f        = GSL_FN_EVAL (&F, xval);
      gdouble fs       = ncm_spline_eval (s, xval);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }

    ncm_spline_free (s);
    ncm_vector_free (x);
    ncm_vector_free (y);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function = &F_linear;
    F.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x  = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble f  = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }

    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function = &F_cubic;
    F.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x  = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble f  = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }

    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function = &F_sin_poly;
    F.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * ((test->nknots) / 100.0 - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x  = test->xi + test->dx * (test->nknots / 100.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble f  = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_idx (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  /* Only test for cubic splines which have eval_idx implemented */
  if (!NCM_IS_SPLINE_CUBIC (test->s_base))
  {
    g_test_skip ("eval_idx only implemented for cubic splines");

    return;
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function = &F_linear;
    F.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x    = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      guint idx    = ncm_spline_get_index (s, x);
      gdouble f    = GSL_FN_EVAL (&F, x);
      gdouble fs   = ncm_spline_eval (s, x);
      gdouble fsid = ncm_spline_eval_idx (s, x, idx);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
      ncm_assert_cmpdouble_e (fsid, ==, fs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function = &F_cubic;
    F.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x    = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      guint idx    = ncm_spline_get_index (s, x);
      gdouble f    = GSL_FN_EVAL (&F, x);
      gdouble fs   = ncm_spline_eval (s, x);
      gdouble fsid = ncm_spline_eval_idx (s, x, idx);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
      ncm_assert_cmpdouble_e (fsid, ==, fs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function = &F_sin_poly;
    F.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * ((test->nknots) / 100.0 - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x    = test->xi + test->dx * (test->nknots / 100.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      guint idx    = ncm_spline_get_index (s, x);
      gdouble f    = GSL_FN_EVAL (&F, x);
      gdouble fs   = ncm_spline_eval (s, x);
      gdouble fsid = ncm_spline_eval_idx (s, x, idx);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
      ncm_assert_cmpdouble_e (fsid, ==, fs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_deriv (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  gsl_function F_deriv;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function       = &F_linear;
    F.params         = d;
    F_deriv.function = &F_linear_deriv;
    F_deriv.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x   = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble df  = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error, 0.0);
    }

    ncm_spline_free (s);
  }

  if (test->deriv1)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function       = &F_cubic;
    F.params         = d;
    F_deriv.function = &F_cubic_deriv;
    F_deriv.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x   = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble df  = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error_d1, 0.0);
    }

    ncm_spline_free (s);
  }

  if (test->deriv1)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function       = &F_sin_poly;
    F.params         = d;
    F_deriv.function = &F_sin_poly_deriv;
    F_deriv.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots / 100.0 - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x   = test->xi + test->dx * (test->nknots / 100.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble df  = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error_d1, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_deriv_idx (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  gsl_function F_deriv;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  /* Only test for cubic splines which have deriv_idx implemented */
  if (!NCM_IS_SPLINE_CUBIC (test->s_base))
  {
    g_test_skip ("eval_deriv_idx only implemented for cubic splines");

    return;
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function       = &F_linear;
    F.params         = d;
    F_deriv.function = &F_linear_deriv;
    F_deriv.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x     = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      guint idx     = ncm_spline_get_index (s, x);
      gdouble df    = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs   = ncm_spline_eval_deriv (s, x);
      gdouble dfsid = ncm_spline_eval_deriv_idx (s, x, idx);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error, 0.0);
      ncm_assert_cmpdouble_e (dfsid, ==, dfs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }

  if (test->deriv1)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function       = &F_cubic;
    F.params         = d;
    F_deriv.function = &F_cubic_deriv;
    F_deriv.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x     = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      guint idx     = ncm_spline_get_index (s, x);
      gdouble df    = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs   = ncm_spline_eval_deriv (s, x);
      gdouble dfsid = ncm_spline_eval_deriv_idx (s, x, idx);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error_d1, 0.0);
      ncm_assert_cmpdouble_e (dfsid, ==, dfs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }

  if (test->deriv1)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function       = &F_sin_poly;
    F.params         = d;
    F_deriv.function = &F_sin_poly_deriv;
    F_deriv.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots / 100.0 - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x     = test->xi + test->dx * (test->nknots / 100.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      guint idx     = ncm_spline_get_index (s, x);
      gdouble df    = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs   = ncm_spline_eval_deriv (s, x);
      gdouble dfsid = ncm_spline_eval_deriv_idx (s, x, idx);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error_d1, 0.0);
      ncm_assert_cmpdouble_e (dfsid, ==, dfs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_deriv2 (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  gsl_function F_deriv2;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  if (test->deriv2)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function        = &F_linear;
    F.params          = d;
    F_deriv2.function = &F_linear_deriv2;
    F_deriv2.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x    = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble d2f  = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, test->error_d2, 0.0);
    }

    ncm_spline_free (s);
  }

  if (test->deriv2)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function        = &F_cubic;
    F.params          = d;
    F_deriv2.function = &F_cubic_deriv2;
    F_deriv2.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x    = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble d2f  = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, test->error_d2, 0.0);
    }

    ncm_spline_free (s);
  }

  if (FALSE)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function        = &F_sin_poly;
    F.params          = d;
    F_deriv2.function = &F_sin_poly_deriv2;
    F_deriv2.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots / 500.0 - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x    = test->xi + test->dx * (test->nknots / 500.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble d2f  = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, test->error_d2, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_int (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  gsl_function F_int;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function     = &F_linear;
    F.params       = d;
    F_int.function = &F_linear_int_0x;
    F_int.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x   = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble xi  = test->xi;
      gdouble If  = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs = ncm_spline_eval_integ (s, xi, x);

      ncm_assert_cmpdouble_e (Ifs, ==, If, test->error, 0.0);
    }

    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function     = &F_cubic;
    F.params       = d;
    F_int.function = &F_cubic_int_0x;
    F_int.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x   = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble xi  = test->xi;
      gdouble If  = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs = ncm_spline_eval_integ (s, xi, x);

      ncm_assert_cmpdouble_e (Ifs, ==, If, test->error, 0.0);
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_integ_idx (TestNcmSpline *test, gconstpointer pdata)
{
  gsl_function F;
  gsl_function F_int;
  guint i;
  gdouble d[4];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  /* Only test for cubic splines which have integ_idx implemented */
  if (!NCM_IS_SPLINE_CUBIC (test->s_base))
  {
    g_test_skip ("eval_integ_idx only implemented for cubic splines");

    return;
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function     = &F_linear;
    F.params       = d;
    F_int.function = &F_linear_int_0x;
    F_int.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x     = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble xi    = test->xi;
      guint idx_i   = ncm_spline_get_index (s, xi);
      guint idx_f   = ncm_spline_get_index (s, x);
      gdouble If    = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs   = ncm_spline_eval_integ (s, xi, x);
      gdouble Ifsid = ncm_spline_eval_integ_idx (s, xi, idx_i, x, idx_f);

      ncm_assert_cmpdouble_e (Ifs, ==, If, test->error, 0.0);
      ncm_assert_cmpdouble_e (Ifsid, ==, Ifs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);

    F.function     = &F_cubic;
    F.params       = d;
    F_int.function = &F_cubic_int_0x;
    F_int.params   = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, 0, test->prec);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x     = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble xi    = test->xi;
      guint idx_i   = ncm_spline_get_index (s, xi);
      guint idx_f   = ncm_spline_get_index (s, x);
      gdouble If    = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs   = ncm_spline_eval_integ (s, xi, x);
      gdouble Ifsid = ncm_spline_eval_integ_idx (s, xi, idx_i, x, idx_f);

      ncm_assert_cmpdouble_e (Ifs, ==, If, test->error, 0.0);
      ncm_assert_cmpdouble_e (Ifsid, ==, Ifs, 1.0e-15, 0.0);
    }

    ncm_spline_free (s);
  }
}

static void _test_ncm_spline_get_index_worker (NcmSpline *s, NcmVector *xv, NcmVector *yv, const guint size);

void
test_ncm_spline_get_index (TestNcmSpline *test, gconstpointer pdata)
{
  NcmSpline *s     = ncm_spline_copy (test->s_base);
  const guint size = 30000;
  NcmVector *xv    = ncm_vector_new (size);
  NcmVector *yv    = ncm_vector_new (size);

  _test_ncm_spline_get_index_worker (s, xv, yv, size);

  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_spline_free (s);
}

void
test_ncm_spline_get_index_stride2 (TestNcmSpline *test, gconstpointer pdata)
{
  if (NCM_IS_SPLINE_GSL (test->s_base))
  {
    g_test_skip ("GSL spline does not support stride");
  }
  else
  {
    NcmSpline *s       = ncm_spline_copy (test->s_base);
    const guint size   = 30000;
    NcmVector *xv_full = ncm_vector_new (2 * size);
    NcmVector *xv      = ncm_vector_get_subvector_stride (xv_full, 0, size, 2);
    NcmVector *yv      = ncm_vector_new (size);

    _test_ncm_spline_get_index_worker (s, xv, yv, size);

    ncm_vector_free (xv_full);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_get_index_stride5 (TestNcmSpline *test, gconstpointer pdata)
{
  if (NCM_IS_SPLINE_GSL (test->s_base))
  {
    g_test_skip ("GSL spline does not support stride");
  }
  else
  {
    NcmSpline *s       = ncm_spline_copy (test->s_base);
    const guint size   = 10000;
    NcmVector *xv_full = ncm_vector_new (5 * size);
    NcmVector *xv      = ncm_vector_get_subvector_stride (xv_full, 0, size, 5);
    NcmVector *yv      = ncm_vector_new (size);

    _test_ncm_spline_get_index_worker (s, xv, yv, size);

    ncm_vector_free (xv_full);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_get_index_acc (TestNcmSpline *test, gconstpointer pdata)
{
  NcmSpline *s     = ncm_spline_copy (test->s_base);
  const guint size = 30000;
  NcmVector *xv    = ncm_vector_new (size);
  NcmVector *yv    = ncm_vector_new (size);

  ncm_spline_acc (s, TRUE);
  _test_ncm_spline_get_index_worker (s, xv, yv, size);

  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_spline_free (s);
}

void
test_ncm_spline_get_index_acc_stride2 (TestNcmSpline *test, gconstpointer pdata)
{
  if (NCM_IS_SPLINE_GSL (test->s_base))
  {
    g_test_skip ("GSL spline does not support stride");
  }
  else
  {
    NcmSpline *s       = ncm_spline_copy (test->s_base);
    const guint size   = 30000;
    NcmVector *xv_full = ncm_vector_new (2 * size);
    NcmVector *xv      = ncm_vector_get_subvector_stride (xv_full, 0, size, 2);
    NcmVector *yv      = ncm_vector_new (size);

    ncm_spline_acc (s, TRUE);
    _test_ncm_spline_get_index_worker (s, xv, yv, size);

    ncm_vector_free (xv_full);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_get_index_acc_stride5 (TestNcmSpline *test, gconstpointer pdata)
{
  if (NCM_IS_SPLINE_GSL (test->s_base))
  {
    g_test_skip ("GSL spline does not support stride");
  }
  else
  {
    NcmSpline *s       = ncm_spline_copy (test->s_base);
    const guint size   = 10000;
    NcmVector *xv_full = ncm_vector_new (5 * size);
    NcmVector *xv      = ncm_vector_get_subvector_stride (xv_full, 0, size, 5);
    NcmVector *yv      = ncm_vector_new (size);

    ncm_spline_acc (s, TRUE);
    _test_ncm_spline_get_index_worker (s, xv, yv, size);

    ncm_vector_free (xv_full);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_spline_free (s);
  }
}

static void
_test_ncm_spline_get_index_worker (NcmSpline *s, NcmVector *xv, NcmVector *yv, const guint size)
{
  gdouble xi = -3.0123;
  gdouble d[4];
  guint i;

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  for (i = 0; i < size; i++)
  {
    ncm_vector_set (xv, i, xi);
    ncm_vector_set (yv, i, F_cubic (xi, d));
    xi += g_test_rand_double_range (1.0e-6, 0.1);
  }

  ncm_spline_set (s, xv, yv, TRUE);

  {
    const gdouble xmin = ncm_vector_get (xv, 0);
    const gdouble xmax = ncm_vector_get (xv, size - 1);

    g_assert_cmpuint (ncm_spline_get_index (s, xmin), ==, 0);
    g_assert_cmpuint (ncm_spline_get_index (s, xmax), ==, size - 2);

    g_assert_cmpuint (ncm_spline_get_index (s, xmin * (1.0 + GSL_DBL_EPSILON)), ==, 0);
    g_assert_cmpuint (ncm_spline_get_index (s, xmin * (1.0 - GSL_DBL_EPSILON)), ==, 0);

    g_assert_cmpuint (ncm_spline_get_index (s, xmax * (1.0 + GSL_DBL_EPSILON)), ==, size - 2);
    g_assert_cmpuint (ncm_spline_get_index (s, xmax * (1.0 - GSL_DBL_EPSILON)), ==, size - 2);

    for (i = 1; i < size - 1; i++)
    {
      gdouble x = ncm_vector_get (xv, i);
      guint idx;

      g_assert_cmpuint (ncm_spline_get_index (s, x), ==, i);

      idx = ncm_spline_get_index (s, x * (1.0 + GSL_DBL_EPSILON));
      g_assert_cmpfloat (ncm_vector_get (xv, idx), <=, x);
      g_assert_cmpfloat (ncm_vector_get (xv, idx + 1), >=, x);

      idx = ncm_spline_get_index (s, x * (1.0 - GSL_DBL_EPSILON));
      g_assert_cmpfloat (ncm_vector_get (xv, idx), <=, x);
      g_assert_cmpfloat (ncm_vector_get (xv, idx + 1), >=, x);
    }

    for (i = 0; i < size; i++)
    {
      const gdouble x = g_test_rand_double_range (xmin, xmax);
      const guint idx = ncm_spline_get_index (s, x);

      g_assert_cmpfloat (ncm_vector_get (xv, idx), <=, x);
      g_assert_cmpfloat (ncm_vector_get (xv, idx + 1), >=, x);
    }
  }
}

void
test_ncm_spline_invalid_vector_sizes (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *x = ncm_vector_new (test->nknots);
  NcmVector *y = ncm_vector_new (test->nknots + 1);
  NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_min_vector_sizes (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *x = ncm_vector_new (ncm_spline_min_size (test->s_base) - 1);
  NcmVector *y = ncm_vector_new (ncm_spline_min_size (test->s_base) - 1);
  NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_x_vector (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *x = NULL;
  NcmVector *y = ncm_vector_new (test->nknots);
  NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_y_vector (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *x = ncm_vector_new (test->nknots);
  NcmVector *y = NULL;
  NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_xy_vector (TestNcmSpline *test, gconstpointer pdata)
{
  NcmVector *x = NULL;
  NcmVector *y = NULL;
  NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_array_sizes (TestNcmSpline *test, gconstpointer pdata)
{
  GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
  GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots + 1);

  g_array_set_size (x, test->nknots);
  g_array_set_size (x, test->nknots + 1);

  NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_min_array_sizes (TestNcmSpline *test, gconstpointer pdata)
{
  GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
  GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);

  g_array_set_size (x, ncm_spline_min_size (test->s_base) - 1);
  g_array_set_size (y, ncm_spline_min_size (test->s_base) - 1);

  NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_x_array (TestNcmSpline *test, gconstpointer pdata)
{
  GArray *x = NULL;
  GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);

  g_array_set_size (y, test->nknots);

  NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_y_array (TestNcmSpline *test, gconstpointer pdata)
{
  GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
  GArray *y = NULL;

  g_array_set_size (x, test->nknots);

  NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_invalid_xy_array (TestNcmSpline *test, gconstpointer pdata)
{
  GArray *x    = NULL;
  GArray *y    = NULL;
  NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);

  ncm_spline_free (s);
}

void
test_ncm_spline_traps (TestNcmSpline *test, gconstpointer pdata)
{
  guint i;

  for (i = 0; _test_ncm_spline_traps[i].func != NULL; i++)
  {
    gchar *tname = g_strdup_printf ("/ncm/%s%s", test->name, _test_ncm_spline_traps[i].name);

    g_test_trap_subprocess (tname, 0, 0);
    g_test_trap_assert_failed ();
    g_free (tname);
  }
}


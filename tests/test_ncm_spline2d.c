/***************************************************************************
 *            test_ncm_spline2d.c
 *
 *  Thu April 27 00:40:15 2012
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

#include <glib.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>

typedef struct _TestNcmSpline2d
{
  NcmSpline2d *s2d_base;
  NcmSpline *s_base;
  gsl_function F;
  gdouble test_error;
} TestNcmSpline2d;

void test_ncm_spline2d_bicubic_notaknot_new_empty (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_gsl_cspline_new_empty (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_spline_new_empty (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_new (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_copy_empty (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_copy (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_serialize (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_serialize_init (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_no_stride (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_eval (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_eval_acc (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_deriv (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_eval_integ_dx (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_eval_integ_dy (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_eval_integ_dxdy (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_eval_integ_x_y_xy_spline (TestNcmSpline2d *test, gconstpointer pdata);
void test_ncm_spline2d_free_empty (TestNcmSpline2d *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  /* bi-cubic not-a-knot */

  g_test_add ("/ncm/spline2d_bicubic/notaknot/new", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_new,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/copy_empty", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_copy_empty,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/copy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_copy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/serialize", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_serialize,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/serialize/init", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_serialize_init,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/eval", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_eval,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/eval/acc", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_eval_acc,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/deriv", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_deriv,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/eval_integ_dx", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_eval_integ_dx,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/eval_integ_dy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_eval_integ_dy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/eval_integ_dxdy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_eval_integ_dxdy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_bicubic/notaknot/eval_integ_spline", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_bicubic_notaknot_new_empty,
              &test_ncm_spline2d_eval_integ_x_y_xy_spline,
              &test_ncm_spline2d_free_empty);

  /* gsl_cspline */

  g_test_add ("/ncm/spline2d_gsl/cspline/new", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_new,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/copy_empty", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_copy_empty,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/copy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_copy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/serialize", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_serialize,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/serialize/init", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_serialize_init,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/eval", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_eval,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/eval/acc", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_eval_acc,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/deriv", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_deriv,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/eval_integ_dx", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_eval_integ_dx,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/eval_integ_dy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_eval_integ_dy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_gsl/cspline/eval_integ_dxdy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_gsl_cspline_new_empty,
              &test_ncm_spline2d_eval_integ_dxdy,
              &test_ncm_spline2d_free_empty);

  g_test_add ("/ncm/spline2d_spline/notaknot/new", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_new,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/copy_empty", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_copy_empty,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/copy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_copy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/serialize", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_serialize,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/serialize/init", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_serialize_init,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/no_stride", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_no_stride,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/eval", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_eval,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/eval/acc", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_eval_acc,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/eval_integ_dx", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_eval_integ_dx,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/eval_integ_dy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_eval_integ_dy,
              &test_ncm_spline2d_free_empty);
  g_test_add ("/ncm/spline2d_spline/notaknot/eval_integ_dxdy", TestNcmSpline2d, NULL,
              &test_ncm_spline2d_spline_new_empty,
              &test_ncm_spline2d_eval_integ_dxdy,
              &test_ncm_spline2d_free_empty);

  g_test_run ();
}

static gdouble
F_linear (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return (d[0] + d[1] * x) * y;
}

static gdouble
F_linear_dx (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[1] * y;
}

static gdouble
F_linear_dy (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[0] + d[1] * x;
}

static gdouble
F_linear_dxdy (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return d[1];
}

static gdouble
F_linear_d2x (gdouble x, gdouble y, gpointer p)
{
  return 0.0;
}

static gdouble
F_linear_d2y (gdouble x, gdouble y, gpointer p)
{
  return 0.0;
}

static gdouble
F_linear_intx (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return (d[0] +  0.5 * d[1] * x) * x * y;
}

static gdouble
F_linear_inty (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return 0.5 * (d[0] + d[1] * x) * y * y;
}

static gdouble
F_linear_intxy (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;

  return 0.5 * (d[0] + 0.5 * d[1] * x) * x * y * y;
}

static gdouble
F_poly (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;

  return (d[0] + d[1] * x + d[2] * x2 + d[3] * x3) +
         (d[4] + d[0] * x + d[3] * x2 + d[1] * x3) * y +
         (d[2] + d[3] * x + d[4] * x2 + d[0] * x3) * y2 +
         (d[1] + d[2] * x + d[0] * x2 + d[4] * x3) * y3;
}

static gdouble
F_poly_intx (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble x4 = x3 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;

  return (d[0] * x + 0.5 * d[1] * x2 + d[2] * x3 / 3.0 + 0.25 * d[3] * x4) +
         (d[4] * x + 0.5 * d[0] * x2 + d[3] * x3 / 3.0 + 0.25 * d[1] * x4) * y +
         (d[2] * x + 0.5 * d[3] * x2 + d[4] * x3 / 3.0 + 0.25 * d[0] * x4) * y2 +
         (d[1] * x + 0.5 * d[2] * x2 + d[0] * x3 / 3.0 + 0.25 * d[4] * x4) * y3;
}

static gdouble
F_poly_inty (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  gdouble y4 = y3 * y;

  return (d[0] + d[1] * x + d[2] * x2 + d[3] * x3) * y +
         (d[4] + d[0] * x + d[3] * x2 + d[1] * x3) * 0.5 * y2 +
         (d[2] + d[3] * x + d[4] * x2 + d[0] * x3) * y3 / 3.0 +
         (d[1] + d[2] * x + d[0] * x2 + d[4] * x3) * y4 / 4.0;
}

static gdouble
F_poly_intxy (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *) p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble x4 = x3 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  gdouble y4 = y3 * y;

  return (d[0] * x + 0.5 * d[1] * x2 + d[2] * x3 / 3.0 + 0.25 * d[3] * x4) * y +
         (d[4] * x + 0.5 * d[0] * x2 + d[3] * x3 / 3.0 + 0.25 * d[1] * x4) * 0.5 * y2 +
         (d[2] * x + 0.5 * d[3] * x2 + d[4] * x3 / 3.0 + 0.25 * d[0] * x4) * y3 / 3.0 +
         (d[1] * x + 0.5 * d[2] * x2 + d[0] * x3 / 3.0 + 0.25 * d[4] * x4) * y4 / 4.0;
}

static gdouble
F_func (gdouble x, gdouble y, gpointer p)
{
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;

  return (16.0 +  15.0 * x +  14.0 * x2 +  13.0 * x3 +
          (12.0 +  11.0 * x +  10.0 * x2 +  9.0 * x3) * y +
          (8.0 + 7.0 * x + 6.0 * x2 + 5.0 * x3) * y2 +
          (4.0 + 3.0 * x + 2.0 * x2 + x3) * y3
         ) * cos (x * y * 0.01) * exp (0.001 * x * y);
}

static gdouble
F_func_x (gdouble x, gpointer p)
{
  return F_func (x, 15.0, p);
}

static gdouble
F_func_y (gdouble y, gpointer p)
{
  return F_func (55.0, y, p);
}

#define _NCM_SPLINE2D_TEST_NKNOTS_X 50
#define _NCM_SPLINE2D_TEST_NKNOTS_Y 25
#define _NCM_SPLINE2D_TEST_XI 10.0
#define _NCM_SPLINE2D_TEST_DX 0.3
#define _NCM_SPLINE2D_TEST_YI 7.0
#define _NCM_SPLINE2D_TEST_DY 0.15

void
test_ncm_spline2d_bicubic_notaknot_new_empty (TestNcmSpline2d *test, gconstpointer pdata)
{
  test->s_base     = NULL;
  test->s2d_base   = ncm_spline2d_bicubic_notaknot_new ();
  test->test_error = 1.0e-4;

  g_assert_true (NCM_IS_SPLINE2D_BICUBIC (test->s2d_base));
}

void
test_ncm_spline2d_gsl_cspline_new_empty (TestNcmSpline2d *test, gconstpointer pdata)
{
  test->s_base     = NULL;
  test->s2d_base   = ncm_spline2d_gsl_natural_new ();
  test->test_error = 1.0e-3;

  g_assert_true (NCM_IS_SPLINE2D_GSL (test->s2d_base));
}

void
test_ncm_spline2d_spline_new_empty (TestNcmSpline2d *test, gconstpointer pdata)
{
  test->s_base     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  test->s2d_base   = ncm_spline2d_spline_new (test->s_base);
  test->test_error = 1.0e-4;

  g_assert_true (NCM_IS_SPLINE2D_SPLINE (test->s2d_base));
}

void
test_ncm_spline2d_free_empty (TestNcmSpline2d *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_spline2d_free, test->s2d_base);

  if (test->s_base != NULL)
    NCM_TEST_FREE (ncm_spline_free, test->s_base);
}

void
test_ncm_spline2d_new_sanity (NcmSpline2d *s2d)
{
  NcmSpline2d *s2d_ref = ncm_spline2d_ref (s2d);

  ncm_spline2d_clear (&s2d_ref);

  g_assert_true (s2d_ref == NULL);
  g_assert_true (NCM_IS_SPLINE2D (s2d));
}

void
test_ncm_spline2d_new (TestNcmSpline2d *test, gconstpointer pdata)
{
  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    test_ncm_spline2d_new_sanity (s2d);
    g_assert_true (!ncm_spline2d_is_init (s2d));
    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d;
    gdouble d[2];
    gint i, j;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_vector_set (yv, j, _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_vector_set (xv, i, _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i);
        ncm_matrix_set (zm, j, i, F_linear (ncm_vector_get (xv, i), ncm_vector_get (yv, j), d));
      }
    }

    s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, TRUE);
    test_ncm_spline2d_new_sanity (s2d);
    g_assert_true (ncm_spline2d_is_init (s2d));
    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }
}

void
test_ncm_spline2d_copy_empty (TestNcmSpline2d *test, gconstpointer pdata)
{
  NcmSpline2d *s2d = ncm_spline2d_copy_empty (test->s2d_base);

  g_assert_true (G_TYPE_FROM_INSTANCE (s2d) == G_TYPE_FROM_INSTANCE (test->s2d_base));
  NCM_TEST_FREE (ncm_spline2d_free, s2d);
}

void
test_ncm_spline2d_copy (TestNcmSpline2d *test, gconstpointer pdata)
{
  NcmSpline2d *s2d = ncm_spline2d_copy_empty (test->s2d_base);
  NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
  NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);

  ncm_vector_set_all (xv, M_PI);
  ncm_vector_set_all (yv, M_PI_2);
  ncm_matrix_set_identity (zm);
  ncm_spline2d_set (s2d, xv, yv, zm, FALSE);

  {
    NcmSpline2d *s2d_cp = ncm_spline2d_copy (s2d);
    guint i, j;

    g_assert_true (ncm_spline2d_peek_xv (s2d_cp) != ncm_spline2d_peek_xv (s2d) && ncm_spline2d_peek_yv (s2d_cp) != ncm_spline2d_peek_yv (s2d) && ncm_spline2d_peek_zm (s2d_cp) != ncm_spline2d_peek_zm (s2d));

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_yv (s2d_cp), j), ==, ncm_vector_get (ncm_spline2d_peek_yv (s2d), j));

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_xv (s2d_cp), i), ==, ncm_vector_get (ncm_spline2d_peek_xv (s2d), i));
        ncm_assert_cmpdouble (ncm_matrix_get (ncm_spline2d_peek_zm (s2d_cp), j, i), ==, ncm_matrix_get (ncm_spline2d_peek_zm (s2d), j, i));
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d_cp);
  }
  NCM_TEST_FREE (ncm_spline2d_free, s2d);
  NCM_TEST_FREE (ncm_vector_free, xv);
  NCM_TEST_FREE (ncm_vector_free, yv);
  NCM_TEST_FREE (ncm_matrix_free, zm);
}

void
test_ncm_spline2d_serialize (TestNcmSpline2d *test, gconstpointer pdata)
{
  NcmSpline2d *s2d  = ncm_spline2d_copy_empty (test->s2d_base);
  NcmVector *xv     = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmVector *yv     = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
  NcmMatrix *zm     = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

  ncm_vector_set_all (xv, M_PI);
  ncm_vector_set_all (yv, M_PI_2);
  ncm_matrix_set_identity (zm);
  ncm_spline2d_set (s2d, xv, yv, zm, FALSE);

  {
    NcmSpline2d *s2d_cp = NCM_SPLINE2D (ncm_serialize_dup_obj (ser, G_OBJECT (s2d)));
    guint i, j;

    g_assert_true (ncm_spline2d_peek_xv (s2d_cp) != ncm_spline2d_peek_xv (s2d) && ncm_spline2d_peek_yv (s2d_cp) != ncm_spline2d_peek_yv (s2d) && ncm_spline2d_peek_zm (s2d_cp) != ncm_spline2d_peek_zm (s2d));

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_yv (s2d_cp), j), ==, ncm_vector_get (ncm_spline2d_peek_yv (s2d), j));

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_xv (s2d_cp), i), ==, ncm_vector_get (ncm_spline2d_peek_xv (s2d), i));
        ncm_assert_cmpdouble (ncm_matrix_get (ncm_spline2d_peek_zm (s2d_cp), j, i), ==, ncm_matrix_get (ncm_spline2d_peek_zm (s2d), j, i));
      }
    }

    NCM_TEST_FREE (ncm_serialize_free, ser);
    NCM_TEST_FREE (ncm_spline2d_free, s2d_cp);
  }

  NCM_TEST_FREE (ncm_spline2d_free, s2d);
  NCM_TEST_FREE (ncm_vector_free, xv);
  NCM_TEST_FREE (ncm_vector_free, yv);
  NCM_TEST_FREE (ncm_matrix_free, zm);
}

void
test_ncm_spline2d_serialize_init (TestNcmSpline2d *test, gconstpointer pdata)
{
  NcmSpline2d *s2d  = ncm_spline2d_copy_empty (test->s2d_base);
  NcmVector *xv     = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmVector *yv     = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
  NcmMatrix *zm     = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gint i;

  for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
    ncm_vector_set (xv, i, 1.0 * i);

  for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_Y; i++)
    ncm_vector_set (yv, i, 1.0e-2 * i);

  ncm_matrix_set_identity (zm);

  {
    NcmVector *xv_sub = ncm_vector_get_subvector (xv, 0, _NCM_SPLINE2D_TEST_NKNOTS_X - 1);
    NcmVector *yv_sub = ncm_vector_get_subvector (xv, 0, _NCM_SPLINE2D_TEST_NKNOTS_Y - 1);
    NcmMatrix *zm_sub = ncm_matrix_get_submatrix (zm, 0, 0, _NCM_SPLINE2D_TEST_NKNOTS_Y - 1, _NCM_SPLINE2D_TEST_NKNOTS_X - 1);

    ncm_spline2d_set (s2d, xv_sub, yv_sub, zm_sub, TRUE);
    ncm_vector_free (xv_sub);
    ncm_vector_free (yv_sub);
    ncm_matrix_free (zm_sub);
  }
  ncm_spline2d_set (s2d, xv, yv, zm, TRUE);

  {
    NcmSpline2d *s2d_cp = NCM_SPLINE2D (ncm_serialize_dup_obj (ser, G_OBJECT (s2d)));
    guint i, j;

    g_assert_true (ncm_spline2d_peek_xv (s2d_cp) != ncm_spline2d_peek_xv (s2d) && ncm_spline2d_peek_yv (s2d_cp) != ncm_spline2d_peek_yv (s2d) && ncm_spline2d_peek_zm (s2d_cp) != ncm_spline2d_peek_zm (s2d));

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_yv (s2d_cp), j), ==, ncm_vector_get (ncm_spline2d_peek_yv (s2d), j));

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_xv (s2d_cp), i), ==, ncm_vector_get (ncm_spline2d_peek_xv (s2d), i));
        ncm_assert_cmpdouble (ncm_matrix_get (ncm_spline2d_peek_zm (s2d_cp), j, i), ==, ncm_matrix_get (ncm_spline2d_peek_zm (s2d), j, i));
      }
    }

    NCM_TEST_FREE (ncm_serialize_free, ser);
    NCM_TEST_FREE (ncm_spline2d_free, s2d_cp);
  }

  NCM_TEST_FREE (ncm_spline2d_free, s2d);
  NCM_TEST_FREE (ncm_vector_free, xv);
  NCM_TEST_FREE (ncm_vector_free, yv);
  NCM_TEST_FREE (ncm_matrix_free, zm);
}

void
test_ncm_spline2d_no_stride (TestNcmSpline2d *test, gconstpointer pdata)
{
  NcmSpline2d *s2d = ncm_spline2d_copy_empty (test->s2d_base);
  gdouble xv_buf[2 * _NCM_SPLINE2D_TEST_NKNOTS_X];
  NcmVector *xv     = ncm_vector_new_full (xv_buf, _NCM_SPLINE2D_TEST_NKNOTS_X, 2, NULL, NULL);
  NcmVector *yv     = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
  NcmMatrix *zm     = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gint i;

  for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
    ncm_vector_set (xv, i, 1.0 * i);

  for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_Y; i++)
    ncm_vector_set (yv, i, 1.0e-2 * i);

  ncm_matrix_set_identity (zm);
  ncm_spline2d_set (s2d, xv, yv, zm, TRUE);

  {
    NcmSpline2d *s2d_cp = NCM_SPLINE2D (ncm_serialize_dup_obj (ser, G_OBJECT (s2d)));
    guint i, j;

    g_assert_true (ncm_spline2d_peek_xv (s2d_cp) != ncm_spline2d_peek_xv (s2d) && ncm_spline2d_peek_yv (s2d_cp) != ncm_spline2d_peek_yv (s2d) && ncm_spline2d_peek_zm (s2d_cp) != ncm_spline2d_peek_zm (s2d));

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_yv (s2d_cp), j), ==, ncm_vector_get (ncm_spline2d_peek_yv (s2d), j));

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (ncm_spline2d_peek_xv (s2d_cp), i), ==, ncm_vector_get (ncm_spline2d_peek_xv (s2d), i));
        ncm_assert_cmpdouble (ncm_matrix_get (ncm_spline2d_peek_zm (s2d_cp), j, i), ==, ncm_matrix_get (ncm_spline2d_peek_zm (s2d), j, i));
      }
    }

    NCM_TEST_FREE (ncm_serialize_free, ser);
    NCM_TEST_FREE (ncm_spline2d_free, s2d_cp);
  }

  NCM_TEST_FREE (ncm_spline2d_free, s2d);
  NCM_TEST_FREE (ncm_vector_free, xv);
  NCM_TEST_FREE (ncm_vector_free, yv);
  NCM_TEST_FREE (ncm_matrix_free, zm);
}

void
test_ncm_spline2d_eval (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_linear (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_poly (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmSpline2d *s2d = ncm_spline2d_copy_empty (test->s2d_base);
    gsl_function Fx, Fy;

    Fx.function = &F_func_x;
    Fx.params   = &d;
    Fy.function = &F_func_y;
    Fy.params   = &d;

    guint npx, npy;

    ncm_spline2d_set_function (s2d, NCM_SPLINE_FUNCTION_SPLINE, &Fx, &Fy,
                               _NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0),
                               _NCM_SPLINE2D_TEST_YI, _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0),
                               1.0e-7);

    npx = ncm_vector_len (ncm_spline2d_peek_xv (s2d));
    npy = ncm_vector_len (ncm_spline2d_peek_yv (s2d));

    for (j = 0; j < npy; j++)
    {
      gdouble y = ncm_vector_get (ncm_spline2d_peek_yv (s2d), j);

      for (i = 0; i < npx; i++)
      {
        gdouble x = ncm_vector_get (ncm_spline2d_peek_xv (s2d), i);
        gdouble f = F_func (x, y, d);

        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
        /*printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f); */
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * npy; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * npy - 1.0) * j;

      for (i = 0; i < 2 * npx; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * npx - 1.0) * i;
        gdouble f  = F_func (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
  }
}

void
test_ncm_spline2d_eval_acc (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    ncm_spline2d_use_acc (s2d, TRUE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_linear (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    ncm_spline2d_use_acc (s2d, TRUE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_poly (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmSpline2d *s2d = ncm_spline2d_copy_empty (test->s2d_base);
    gsl_function Fx, Fy;

    Fx.function = &F_func_x;
    Fx.params   = &d;
    Fy.function = &F_func_y;
    Fy.params   = &d;

    guint npx, npy;

    ncm_spline2d_set_function (s2d, NCM_SPLINE_FUNCTION_SPLINE, &Fx, &Fy,
                               _NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0),
                               _NCM_SPLINE2D_TEST_YI, _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0),
                               1.0e-7);

    ncm_spline2d_use_acc (s2d, TRUE);

    npx = ncm_vector_len (ncm_spline2d_peek_xv (s2d));
    npy = ncm_vector_len (ncm_spline2d_peek_yv (s2d));

    for (j = 0; j < npy; j++)
    {
      gdouble y = ncm_vector_get (ncm_spline2d_peek_yv (s2d), j);

      for (i = 0; i < npx; i++)
      {
        gdouble x = ncm_vector_get (ncm_spline2d_peek_xv (s2d), i);
        gdouble f = F_func (x, y, d);

        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
        /*printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f); */
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * npy; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * npy - 1.0) * j;

      for (i = 0; i < 2 * npx; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * npx - 1.0) * i;
        gdouble f  = F_func (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
  }
}

void
test_ncm_spline2d_deriv (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        {
          gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
          gdouble f  = F_linear_dx (x, y, d);
          gdouble fs = ncm_spline2d_deriv_dzdx (s2d, x, y);

          ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
        }
        {
          gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
          gdouble f  = F_linear_dy (x, y, d);
          gdouble fs = ncm_spline2d_deriv_dzdy (s2d, x, y);

          ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
        }
        {
          gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
          gdouble f  = F_linear_dxdy (x, y, d);
          gdouble fs = ncm_spline2d_deriv_d2zdxy (s2d, x, y);

          ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
        }
        {
          gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
          gdouble f  = F_linear_d2x (x, y, d);
          gdouble fs = ncm_spline2d_deriv_d2zdx2 (s2d, x, y);

          ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
        }
        {
          gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
          gdouble f  = F_linear_d2y (x, y, d);
          gdouble fs = ncm_spline2d_deriv_d2zdy2 (s2d, x, y);

          ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
        }
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }
}

void
test_ncm_spline2d_eval_integ_dx (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_linear_intx (x, y, d) - F_linear_intx (_NCM_SPLINE2D_TEST_XI, y, d);
        gdouble fs = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_poly_intx (x, y, d) - F_poly_intx (_NCM_SPLINE2D_TEST_XI, y, d);
        gdouble fs = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }
}

void
test_ncm_spline2d_eval_integ_dy (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_linear_inty (x, y, d) - F_linear_inty (x, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_poly_inty (x, y, d) - F_poly_inty (x, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }
}

void
test_ncm_spline2d_eval_integ_dxdy (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_linear_intxy (x, y, d) - F_linear_intxy (_NCM_SPLINE2D_TEST_XI, y, d) - F_linear_intxy (x, _NCM_SPLINE2D_TEST_YI, d) + F_linear_intxy (_NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    g_assert_true (gsl_finite (ncm_spline2dim_integ_total (s2d)));

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x  = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f  = F_poly_intxy (x, y, d) - F_poly_intxy (_NCM_SPLINE2D_TEST_XI, y, d) - F_poly_intxy (x, _NCM_SPLINE2D_TEST_YI, d) + F_poly_intxy (_NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test->test_error, 0.0);
      }
    }

    g_assert_true (gsl_finite (ncm_spline2dim_integ_total (s2d)));

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }
}

void
test_ncm_spline2d_eval_integ_x_y_xy_spline (TestNcmSpline2d *test, gconstpointer pdata)
{
  guint i, j;
  gdouble d[5];

  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x        = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f_int_x  = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);
        gdouble f_int_sx = ncm_spline2d_integ_dx_spline_val (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        gdouble f_int_y  = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble f_int_sy = ncm_spline2d_integ_dy_spline_val (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        gdouble f_int_xy = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble fs_x     = ncm_spline2d_integ_dxdy_spline_x (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble fs_y     = ncm_spline2d_integ_dxdy_spline_y (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (f_int_x, ==, f_int_sx, test->test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_y, ==, f_int_sy, test->test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_x, test->test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_y, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }

  {
    NcmVector *xv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv    = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm    = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (test->s2d_base, xv, yv, zm, FALSE);
    NcmSpline *sx;
    NcmSpline *sy;

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;

      ncm_vector_set (ncm_spline2d_peek_yv (s2d), j, y);

      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);

        ncm_vector_set (ncm_spline2d_peek_xv (s2d), i, x);
        ncm_matrix_set (ncm_spline2d_peek_zm (s2d), j, i, f);
      }
    }

    sx = ncm_spline2d_integ_dx_spline (s2d, _NCM_SPLINE2D_TEST_XI, ncm_vector_get (ncm_spline2d_peek_xv (s2d), _NCM_SPLINE2D_TEST_NKNOTS_X - 1));

    NCM_TEST_FREE (ncm_spline_free, sx);

    sy = ncm_spline2d_integ_dy_spline (s2d, _NCM_SPLINE2D_TEST_YI, ncm_vector_get (ncm_spline2d_peek_yv (s2d), _NCM_SPLINE2D_TEST_NKNOTS_Y - 1));

    NCM_TEST_FREE (ncm_spline_free, sy);

    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;

      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x        = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f_int_x  = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);
        gdouble f_int_sx = ncm_spline2d_integ_dx_spline_val (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        gdouble f_int_y  = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble f_int_sy = ncm_spline2d_integ_dy_spline_val (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        gdouble f_int_xy = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble fs_x     = ncm_spline2d_integ_dxdy_spline_x (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        gdouble fs_y = ncm_spline2d_integ_dxdy_spline_y (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (f_int_sx, ==, f_int_x, test->test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_sy, ==, f_int_y, test->test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_x, test->test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_y, test->test_error, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_spline2d_free, s2d);
    NCM_TEST_FREE (ncm_vector_free, xv);
    NCM_TEST_FREE (ncm_vector_free, yv);
    NCM_TEST_FREE (ncm_matrix_free, zm);
  }
}


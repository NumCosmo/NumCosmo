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

static gdouble
F_linear (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return  (d[0] +  d[1] * x) * y;
}

static gdouble
F_linear_intx (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return  (d[0] +  0.5 * d[1] * x) * x * y;
}

static gdouble
F_linear_inty (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return  0.5 * (d[0] + d[1] * x) * y * y;
}

static gdouble
F_linear_intxy (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return  0.5 * (d[0] + 0.5 * d[1] * x) * x * y * y;
}

static gdouble
F_poly (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  return  d[0] +  d[1] * x +  d[2] * x2 +  d[3] * x3 +
    (d[4] +  d[0] * x +  d[3] * x2 +  d[1] * x3 ) * y +
    (d[2] + d[3] * x + d[4] * x2 + d[0] * x3 ) * y2 +
    (d[1] + d[2] * x + d[0] * x2 + d[4] * x3 ) * y3;
}

static gdouble
F_poly_intx (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble x4 = x3 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  return  d[0] * x +  0.5 * d[1] * x2 +  d[2] * x3 / 3.0 + 0.25 * d[3] * x4 +
    (d[4] * x +  0.5 * d[0] * x2 +  d[3] * x3 / 3.0 + 0.25 * d[1] * x4 ) * y +
    (d[2] * x + 0.5 * d[3] * x2 + d[4] * x3 / 3.0 + 0.25 * d[0] * x4 ) * y2 +
    (d[1] * x + 0.5 * d[2] * x2 + d[0] * x3 / 3.0 + 0.25 * d[4] * x4 ) * y3;
}

static gdouble
F_poly_inty (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  gdouble y4 = y3 * y;
  return  (d[0] +  d[1] * x +  d[2] * x2 +  d[3] * x3) * y +
    (d[4] +  d[0] * x +  d[3] * x2 +  d[1] * x3 ) * 0.5 * y2 +
    (d[2] + d[3] * x + d[4] * x2 + d[0] * x3 ) * y3 / 3.0 +
    (d[1] + d[2] * x + d[0] * x2 + d[4] * x3 ) * y4 / 4.0;
}

static gdouble
F_poly_intxy (gdouble x, gdouble y, gpointer p)
{
  gdouble *d = (gdouble *)p;
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble x4 = x3 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  gdouble y4 = y3 * y;
  return  (d[0] * x + 0.5 * d[1] * x2 + d[2] * x3 /3.0 + 0.25 * d[3] * x4) * y +
    (d[4] * x + 0.5 * d[0] * x2 + d[3] * x3 /3.0 + 0.25 * d[1] * x4 ) * 0.5 * y2 +
    (d[2] * x + 0.5 * d[3] * x2 + d[4] * x3 /3.0 + 0.25 * d[0] * x4 ) * y3 / 3.0 +
    (d[1] * x + 0.5 * d[2] * x2 + d[0] * x3 /3.0 + 0.25 * d[4] * x4 ) * y4 / 4.0;
}

static gdouble
F_func (gdouble x, gdouble y, gpointer p)
{
  gdouble x2 = x * x;
  gdouble x3 = x2 * x;
  gdouble y2 = y * y;
  gdouble y3 = y2 * y;
  return (16.0 +  15.0 * x +  14.0 * x2 +  13.0 * x3 +
          (12.0 +  11.0 * x +  10.0 * x2 +  9.0 * x3 ) * y +
          (8.0 + 7.0 * x + 6.0 * x2 + 5.0 * x3 ) * y2 +
          (4.0 + 3.0 * x + 2.0 * x2 + x3 ) * y3
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

void test_ncm_spline2d_bicubic_notaknot_new_empty (void);
void test_ncm_spline2d_bicubic_spline_new_empty (void);
void test_ncm_spline2d_gsl_cspline_new_empty (void);
void test_ncm_spline2d_gsl_spline_new_empty (void);
void test_ncm_spline2d_spline_new_empty (void);
void test_ncm_spline2d_new (void);
void test_ncm_spline2d_copy_empty (void);
void test_ncm_spline2d_copy (void);
void test_ncm_spline2d_eval (void);
void test_ncm_spline2d_eval_integ_dx (void);
void test_ncm_spline2d_eval_integ_dy (void);
void test_ncm_spline2d_eval_integ_dxdy (void);
void test_ncm_spline2d_eval_integ_x_y_xy_spline (void);
void test_ncm_spline2d_free_empty (void);

NcmSpline2d *s2d_base = NULL;
NcmSpline *s_base = NULL;
gsl_function F;
gdouble test_error;

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/new_empty", &test_ncm_spline2d_bicubic_notaknot_new_empty);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/new", &test_ncm_spline2d_new);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/copy_empty", &test_ncm_spline2d_copy_empty);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/copy", &test_ncm_spline2d_copy);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/eval", &test_ncm_spline2d_eval);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/eval_integ_dx", &test_ncm_spline2d_eval_integ_dx);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/eval_integ_dy", &test_ncm_spline2d_eval_integ_dy);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/eval_integ_dxdy", &test_ncm_spline2d_eval_integ_dxdy);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/eval_integ_spline", &test_ncm_spline2d_eval_integ_x_y_xy_spline);
  g_test_add_func ("/ncm/spline2d_bicubic/notaknot/free_empty", &test_ncm_spline2d_free_empty);

  g_test_add_func ("/ncm/spline2d_gsl/cspline/new_empty", &test_ncm_spline2d_gsl_cspline_new_empty);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/new", &test_ncm_spline2d_new);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/copy_empty", &test_ncm_spline2d_copy_empty);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/copy", &test_ncm_spline2d_copy);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/eval", &test_ncm_spline2d_eval);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/eval_integ_dx", &test_ncm_spline2d_eval_integ_dx);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/eval_integ_dy", &test_ncm_spline2d_eval_integ_dy);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/eval_integ_dxdy", &test_ncm_spline2d_eval_integ_dxdy);
  g_test_add_func ("/ncm/spline2d_gsl/cspline/free_empty", &test_ncm_spline2d_free_empty);

  g_test_add_func ("/ncm/spline2d_spline/new_empty", &test_ncm_spline2d_spline_new_empty);
  g_test_add_func ("/ncm/spline2d_spline/new", &test_ncm_spline2d_new);
  g_test_add_func ("/ncm/spline2d_spline/copy_empty", &test_ncm_spline2d_copy_empty);
  g_test_add_func ("/ncm/spline2d_spline/copy", &test_ncm_spline2d_copy);
  g_test_add_func ("/ncm/spline2d_spline/eval", &test_ncm_spline2d_eval);
  g_test_add_func ("/ncm/spline2d_spline/eval_integ_dx", &test_ncm_spline2d_eval_integ_dx);
  g_test_add_func ("/ncm/spline2d_spline/eval_integ_dy", &test_ncm_spline2d_eval_integ_dy);
  g_test_add_func ("/ncm/spline2d_spline/eval_integ_dxdy", &test_ncm_spline2d_eval_integ_dxdy);
  g_test_add_func ("/ncm/spline2d_spline/free_empty", &test_ncm_spline2d_free_empty);

  g_test_run ();
}

#define _NCM_SPLINE2D_TEST_NKNOTS_X 50
#define _NCM_SPLINE2D_TEST_NKNOTS_Y 25
#define _NCM_SPLINE2D_TEST_XI 10.0
#define _NCM_SPLINE2D_TEST_DX 0.3
#define _NCM_SPLINE2D_TEST_YI 7.0
#define _NCM_SPLINE2D_TEST_DY 0.15

void
test_ncm_spline2d_bicubic_notaknot_new_empty (void)
{
  test_error = 1.0e-4;
  s2d_base = ncm_spline2d_bicubic_notaknot_new ();
  g_assert_true (NCM_IS_SPLINE2D_BICUBIC (s2d_base));
}

/*
 * void
 * test_ncm_spline2d_bicubic_spline_new_empty (void)
 * {
 *   s2d_base = ncm_spline2d_bicubic_new (s_base);
 *   g_assert_true (NCM_IS_SPLINE2D_BICUBIC (s2d_base));
 * }
 */
void
test_ncm_spline2d_gsl_cspline_new_empty (void)
{
  test_error = 1.0e-3;
  s2d_base = ncm_spline2d_gsl_natural_new ();
  g_assert_true (NCM_IS_SPLINE2D_GSL (s2d_base));
}
/*
 * void
 * test_ncm_spline2d_gsl_spline_new_empty (void)
 * {
 *   s2d_base = ncm_spline2d_gsl_new (s_base);
 *   g_assert_true (NCM_IS_SPLINE2D_GSL (s2d_base));
 * }
 */
void
test_ncm_spline2d_spline_new_empty (void)
{
  test_error = 1.0e-4;
  s_base = ncm_spline_cubic_notaknot_new ();
  s2d_base = ncm_spline2d_spline_new (s_base);
  g_assert_true (NCM_IS_SPLINE2D_SPLINE (s2d_base));
}

void
test_ncm_spline2d_free_empty (void)
{
  NCM_TEST_FREE (ncm_spline2d_free, s2d_base);
}

void
test_ncm_spline2d_new_sanity (NcmSpline2d *s2d)
{
  g_assert_true (NCM_IS_SPLINE2D (s2d));
}

void
test_ncm_spline2d_new (void)
{
  {
    NcmVector *x = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *y = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *z = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, x, y, z, FALSE);

    test_ncm_spline2d_new_sanity (s2d);
    g_assert_true (s2d->init == FALSE);
    ncm_spline2d_free (s2d);
  }

  {
    NcmVector *x = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *y = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *z = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d;
    gdouble d[2];
    gint i, j;
    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_vector_set (y, j, _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_vector_set (x, i, _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i);
        ncm_matrix_set (z, j, i, F_linear (ncm_vector_get (x, i), ncm_vector_get (y, j), d));
      }
    }

    s2d = ncm_spline2d_new (s2d_base, x, y, z, TRUE);
    test_ncm_spline2d_new_sanity (s2d);
    g_assert_true (s2d->init == TRUE);
    ncm_spline2d_free (s2d);
  }
}

void
test_ncm_spline2d_copy_empty (void)
{
  NcmSpline2d *s2d = ncm_spline2d_copy_empty (s2d_base);
  g_assert_true (G_TYPE_FROM_INSTANCE (s2d) == G_TYPE_FROM_INSTANCE (s2d_base));
  ncm_spline2d_free (s2d);
}

void
test_ncm_spline2d_copy (void)
{
  NcmSpline2d *s2d = ncm_spline2d_copy_empty (s2d_base);
  NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
  NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
  NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
  ncm_vector_set_all (xv, M_PI);
  ncm_vector_set_all (yv, M_PI_2);
  ncm_matrix_set_identity (zm);
  ncm_spline2d_set (s2d, xv, yv, zm, FALSE);

  {
    NcmSpline2d *s2d_cp = ncm_spline2d_copy (s2d);
    guint i, j;
    g_assert_true (s2d_cp->xv != s2d->xv && s2d_cp->yv != s2d->yv && s2d_cp->zm != s2d->zm);
    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (s2d_cp->yv, j), ==, ncm_vector_get (s2d->yv, j));
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (s2d_cp->xv, i), ==, ncm_vector_get (s2d->xv, i));
        ncm_assert_cmpdouble (ncm_matrix_get (s2d_cp->zm, j, i), ==, ncm_matrix_get (s2d->zm, j, i));
      }
    }

    ncm_spline2d_free (s2d_cp);
  }
  ncm_spline2d_free (s2d);
}


void
test_ncm_spline2d_eval (void)
{
  guint i, j;
  gdouble d[5];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);
    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_linear (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);
    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_poly (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

  {
    NcmSpline2d *s2d = ncm_spline2d_copy_empty (s2d_base);
    gsl_function Fx, Fy;
    Fx.function = &F_func_x;
    Fx.params = &d;
    Fy.function = &F_func_y;
    Fy.params = &d;
    guint npx, npy;

    ncm_spline2d_set_function (s2d, NCM_SPLINE_FUNCTION_SPLINE, &Fx, &Fy,
                               _NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0),
                               _NCM_SPLINE2D_TEST_YI, _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0),
                               1.0e-7);

    npx = ncm_vector_len (s2d->xv);
    npy = ncm_vector_len (s2d->yv);

    for (j = 0; j < npy; j++)
    {
      gdouble y = ncm_vector_get (s2d->yv, j);
      for (i = 0; i < npx; i++)
      {
        gdouble x = ncm_vector_get (s2d->xv, i);
        gdouble f = F_func (x, y, d);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * npy; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * npy - 1.0) * j;
      for (i = 0; i < 2 * npx; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * npx - 1.0) * i;
        gdouble f = F_func (x, y, d);
        gdouble fs = ncm_spline2d_eval (s2d, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

}

void
test_ncm_spline2d_eval_integ_dx (void)
{
  guint i, j;
  gdouble d[5];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }
    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_linear_intx (x, y, d) - F_linear_intx (_NCM_SPLINE2D_TEST_XI, y, d);
        gdouble fs = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);
    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_poly_intx (x, y, d) - F_poly_intx (_NCM_SPLINE2D_TEST_XI, y, d);
        gdouble fs = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }
}

void
test_ncm_spline2d_eval_integ_dy (void)
{
  guint i, j;
  gdouble d[5];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }
    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_linear_inty (x, y, d) - F_linear_inty (x, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);
    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 0; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 0; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_poly_inty (x, y, d) - F_poly_inty (x, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }
}

void
test_ncm_spline2d_eval_integ_dxdy (void)
{
  guint i, j;
  gdouble d[5];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }
    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_linear_intxy (x, y, d) - F_linear_intxy (_NCM_SPLINE2D_TEST_XI, y, d) - F_linear_intxy (x, _NCM_SPLINE2D_TEST_YI, d) + F_linear_intxy (_NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);
    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f = F_poly_intxy (x, y, d) - F_poly_intxy (_NCM_SPLINE2D_TEST_XI, y, d) - F_poly_intxy (x, _NCM_SPLINE2D_TEST_YI, d) + F_poly_intxy (_NCM_SPLINE2D_TEST_XI, _NCM_SPLINE2D_TEST_YI, d);
        gdouble fs = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (fs, ==, f, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }
}

void
test_ncm_spline2d_eval_integ_x_y_xy_spline (void)
{
  guint i, j;
  gdouble d[5];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();
  d[4] = g_test_rand_double ();

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_linear (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }
    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f_int_x = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);
        gdouble f_int_sx = ncm_spline2d_integ_dx_spline_val (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        gdouble f_int_y = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble f_int_sy = ncm_spline2d_integ_dy_spline_val (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        gdouble f_int_xy = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble fs_x = ncm_spline2d_integ_dxdy_spline_x (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble fs_y = ncm_spline2d_integ_dxdy_spline_y (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (f_int_x, ==, f_int_sx, test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_y, ==, f_int_sy, test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_x, test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_y, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }

  {
    NcmVector *xv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmVector *yv = ncm_vector_new (_NCM_SPLINE2D_TEST_NKNOTS_Y);
    NcmMatrix *zm = ncm_matrix_new (_NCM_SPLINE2D_TEST_NKNOTS_Y, _NCM_SPLINE2D_TEST_NKNOTS_X);
    NcmSpline2d *s2d = ncm_spline2d_new (s2d_base, xv, yv, zm, FALSE);
    NcmSpline *sx;
    NcmSpline *sy;

    for (j = 0; j < _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * j;
      ncm_vector_set (s2d->yv, j, y);
      for (i = 0; i < _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * i;
        gdouble f = F_poly (x, y, d);
        ncm_vector_set (s2d->xv, i, x);
        ncm_matrix_set (s2d->zm, j, i, f);
        //printf ("x = %.5g y = %.5g z = %.5g\n", x, y, f);
      }
    }

    sx = ncm_spline2d_integ_dx_spline (s2d, _NCM_SPLINE2D_TEST_XI, ncm_vector_get (s2d->xv, _NCM_SPLINE2D_TEST_NKNOTS_X - 1));

    NCM_TEST_FREE (ncm_spline_free, sx);

    sy = ncm_spline2d_integ_dy_spline (s2d, _NCM_SPLINE2D_TEST_YI, ncm_vector_get (s2d->yv, _NCM_SPLINE2D_TEST_NKNOTS_Y - 1));

    NCM_TEST_FREE (ncm_spline_free, sy);

    ncm_spline2d_prepare (s2d);

    for (j = 1; j < 2 * _NCM_SPLINE2D_TEST_NKNOTS_Y; j++)
    {
      gdouble y = _NCM_SPLINE2D_TEST_YI + _NCM_SPLINE2D_TEST_DY * (_NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_Y - 1.0) * j;
      for (i = 1; i < 2 * _NCM_SPLINE2D_TEST_NKNOTS_X; i++)
      {
        gdouble x = _NCM_SPLINE2D_TEST_XI + _NCM_SPLINE2D_TEST_DX * (_NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) / (2.0 * _NCM_SPLINE2D_TEST_NKNOTS_X - 1.0) * i;
        gdouble f_int_x = ncm_spline2d_integ_dx (s2d, _NCM_SPLINE2D_TEST_XI, x, y);
        gdouble f_int_sx = ncm_spline2d_integ_dx_spline_val (s2d, _NCM_SPLINE2D_TEST_XI, x, y);

        gdouble f_int_y = ncm_spline2d_integ_dy (s2d, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble f_int_sy = ncm_spline2d_integ_dy_spline_val (s2d, x, _NCM_SPLINE2D_TEST_YI, y);

        gdouble f_int_xy = ncm_spline2d_integ_dxdy (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);
        gdouble fs_x = ncm_spline2d_integ_dxdy_spline_x (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        gdouble fs_y = ncm_spline2d_integ_dxdy_spline_y (s2d, _NCM_SPLINE2D_TEST_XI, x, _NCM_SPLINE2D_TEST_YI, y);

        ncm_assert_cmpdouble_e (f_int_sx, ==, f_int_x, test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_sy, ==, f_int_y, test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_x, test_error, 0.0);
        ncm_assert_cmpdouble_e (f_int_xy, ==, fs_y, test_error, 0.0);
      }
    }
    ncm_spline2d_free (s2d);
  }
}

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

#include <glib.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>

gdouble
F_linear (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return d[0] + d[1] * x;
}

gdouble
F_linear_deriv (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *)p;
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
  gdouble *d = (gdouble *)p;
  return d[0] * x + d[1] * x * x * 0.5;
}

gdouble
F_cubic (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return d[0] + (d[1] + (d[2] + d[3] * x) * x) * x;
}

gdouble
F_cubic_deriv (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return d[1] + (2.0 * d[2] + 3.0 * d[3] * x) * x;
}

gdouble
F_cubic_deriv2 (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return 2.0 * d[2] + 6.0 * d[3] * x;
}

gdouble
F_cubic_int_0x (gdouble x, gpointer p)
{
  gdouble *d = (gdouble *)p;
  return (d[0] + (0.5 * d[1] + (d[2] / 3.0 + 0.25 * d[3] * x) * x) * x) * x;
}

gdouble
F_sin_poly (gdouble x, gpointer p)
{
  return sin (M_PI * x) * exp(x / 5.0) +
    (x * x / 2.0 + 5.0 * x + 0.1 * x * x * x + 1.0e-3 * gsl_pow_5 (x));
}

gdouble
F_sin_poly_deriv (gdouble x, gpointer p)
{
  return (M_PI * cos(M_PI * x) + sin (M_PI * x) / 5.0) * exp(x / 5.0) +
    (5.0 + x + 0.3 * x * x + 5.0e-3 * gsl_pow_4 (x));
}

gdouble
F_sin_poly_deriv2 (gdouble x, gpointer p)
{
  return (2.0 * M_PI * cos(M_PI * x) / 5.0 + (1.0/25.0 - M_PI * M_PI) * sin(M_PI * x)) * exp(x / 5.0) +
    (1.0 + 0.6 * x + 2.0e-2 * gsl_pow_3 (x));
}

void test_ncm_spline_cubic_notaknot_new_empty (void);
void test_ncm_spline_gsl_cspline_new_empty (void);
void test_ncm_spline_new (void);
void test_ncm_spline_new_array (void);
void test_ncm_spline_new_data (void);
void test_ncm_spline_copy_empty (void);
void test_ncm_spline_copy (void);
void test_ncm_spline_eval (void);
void test_ncm_spline_eval_deriv (void);
void test_ncm_spline_eval_deriv2 (void);
void test_ncm_spline_eval_int (void);
void test_ncm_spline_free_empty (void);

NcmSpline *s_base = NULL;
gsl_function F;

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/new_empty", &test_ncm_spline_cubic_notaknot_new_empty);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/new", &test_ncm_spline_new);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/new_array", &test_ncm_spline_new_array);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/new_data", &test_ncm_spline_new_data);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/copy_empty", &test_ncm_spline_copy_empty);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/copy", &test_ncm_spline_copy);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/eval", &test_ncm_spline_eval);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/eval_deriv", &test_ncm_spline_eval_deriv);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/eval_deriv2", &test_ncm_spline_eval_deriv2);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/eval_int", &test_ncm_spline_eval_int);
  g_test_add_func ("/numcosmo/ncm_spline_cubic_notaknot/free_empty", &test_ncm_spline_free_empty);

  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/new_empty", &test_ncm_spline_gsl_cspline_new_empty);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/new", &test_ncm_spline_new);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/new_array", &test_ncm_spline_new_array);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/new_data", &test_ncm_spline_new_data);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/copy_empty", &test_ncm_spline_copy_empty);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/copy", &test_ncm_spline_copy);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/eval", &test_ncm_spline_eval);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/eval_deriv", &test_ncm_spline_eval_deriv);
  //g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/eval_deriv2", &test_ncm_spline_eval_deriv2);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/eval_int", &test_ncm_spline_eval_int);
  g_test_add_func ("/numcosmo/ncm_spline_gsl/cspline/free_empty", &test_ncm_spline_free_empty);

  g_test_run ();
}

#define _NCM_SPLINE_TEST_NKNOTS 1000
#define _NCM_SPLINE_TEST_XI 10.0
#define _NCM_SPLINE_TEST_DX 0.3
#define _NCM_SPLINE_TEST_ERROR 1.0e-4

void
test_ncm_spline_cubic_notaknot_new_empty (void)
{
  s_base = ncm_spline_cubic_notaknot_new ();
  g_assert (NCM_IS_SPLINE_CUBIC (s_base));
  g_assert (NCM_IS_SPLINE_CUBIC_NOTAKNOT (s_base));
}

void
test_ncm_spline_gsl_cspline_new_empty (void)
{
  s_base = ncm_spline_gsl_new (gsl_interp_cspline);
  g_assert (NCM_IS_SPLINE_GSL (s_base));
}

void
test_ncm_spline_free_empty (void)
{
  ncm_spline_free (s_base);
  NCM_TEST_FAIL (ncm_spline_free (s_base));
  s_base = NULL;
}

void
test_ncm_spline_new_sanity (NcmSpline *s)
{
  g_assert (NCM_IS_SPLINE (s));
}

void
test_ncm_spline_new (void)
{
  NCM_TEST_FAIL (G_STMT_START
  {
    NcmVector *x = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmVector *y = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS + 1);
    NcmSpline *s = ncm_spline_new (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    NcmVector *x = ncm_vector_new (ncm_spline_min_size (s_base) - 1);
    NcmVector *y = ncm_vector_new (ncm_spline_min_size (s_base) - 1);
    NcmSpline *s = ncm_spline_new (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    NcmVector *x = NULL;
    NcmVector *y = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s = ncm_spline_new (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    NcmVector *x = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmVector *y = NULL;
    NcmSpline *s = ncm_spline_new (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    NcmVector *x = NULL;
    NcmVector *y = NULL;
    NcmSpline *s = ncm_spline_new (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  {
    NcmVector *x = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmVector *y = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s = ncm_spline_new (s_base, x, y, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert (s->init == FALSE);
    ncm_spline_free (s);
  }

  {
    NcmVector *x = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmVector *y = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s;
    gdouble d[2];
    gint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();
    for (i = 0; i < _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      ncm_vector_set (x, i, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * i);
      ncm_vector_set (y, i, F_linear (ncm_vector_get (x, i), d));
    }

    s = ncm_spline_new (s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert (s->init == TRUE);
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_new_array (void)
{
  NCM_TEST_FAIL (G_STMT_START
  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS + 1);
    g_array_set_size (x, _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (x, _NCM_SPLINE_TEST_NKNOTS + 1);
    NcmSpline *s = ncm_spline_new_array (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (x, ncm_spline_min_size (s_base) - 1);
    g_array_set_size (y, ncm_spline_min_size (s_base) - 1);
    NcmSpline *s = ncm_spline_new_array (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    GArray *x = NULL;
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (y, _NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s = ncm_spline_new_array (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    GArray *y = NULL;
    g_array_set_size (x, _NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s = ncm_spline_new_array (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  NCM_TEST_FAIL (G_STMT_START
  {
    GArray *x = NULL;
    GArray *y = NULL;
    NcmSpline *s = ncm_spline_new_array (s_base, x, y, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (x, _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (y, _NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s = ncm_spline_new_array (s_base, x, y, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert (s->init == FALSE);
    ncm_spline_free (s);
  }

  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (x, _NCM_SPLINE_TEST_NKNOTS);
    g_array_set_size (y, _NCM_SPLINE_TEST_NKNOTS);
    NcmSpline *s;
    gdouble d[2];
    gint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();
    for (i = 0; i < _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      g_array_index(x, gdouble, i) = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * i;
      g_array_index(y, gdouble, i) = F_linear (g_array_index(x, gdouble, i), d);
    }

    s = ncm_spline_new_array (s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert (s->init == TRUE);
    ncm_spline_free (s);
  }

}

void
test_ncm_spline_new_data (void)
{

  NCM_TEST_FAIL (G_STMT_START
  {
    gdouble x[_NCM_SPLINE_TEST_NKNOTS];
    gdouble y[_NCM_SPLINE_TEST_NKNOTS];
    NcmSpline *s = ncm_spline_new_data (s_base, x, y, ncm_spline_min_size (s_base) - 1, FALSE);
    ncm_spline_free (s);
  } G_STMT_END);

  {
    gdouble x[_NCM_SPLINE_TEST_NKNOTS];
    gdouble y[_NCM_SPLINE_TEST_NKNOTS];
    NcmSpline *s = ncm_spline_new_data (s_base, x, y, _NCM_SPLINE_TEST_NKNOTS, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert (s->init == FALSE);
    ncm_spline_free (s);
  }

  {
    gdouble x[_NCM_SPLINE_TEST_NKNOTS];
    gdouble y[_NCM_SPLINE_TEST_NKNOTS];
    NcmSpline *s;
    gdouble d[2];
    gint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();
    for (i = 0; i < _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      x[i] = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * i;
      y[i] = F_linear (x[i], d);
    }

    s = ncm_spline_new_data (s_base, x, y, _NCM_SPLINE_TEST_NKNOTS, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert (s->init == TRUE);
    ncm_spline_free (s);
  }

}

void
test_ncm_spline_copy_empty (void)
{
  NcmSpline *s = ncm_spline_copy_empty (s_base);
  g_assert (G_TYPE_FROM_INSTANCE (s) == G_TYPE_FROM_INSTANCE (s_base));
  ncm_spline_free (s);
}

void
test_ncm_spline_copy (void)
{
  NcmVector *xv = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
  NcmVector *yv = ncm_vector_new (_NCM_SPLINE_TEST_NKNOTS);
  guint i;

  for (i = 0; i < _NCM_SPLINE_TEST_NKNOTS; i++)
  {
    ncm_vector_set (xv, i, i * M_PI);
    ncm_vector_set (yv, i, i * M_PI_2);
  }

  ncm_spline_set (s_base, xv, yv, FALSE);
  {
    NcmSpline *s = ncm_spline_copy (s_base);
    g_assert (s->xv != s_base->xv && s->yv != s_base->yv);
    for (i = 0; i < _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (s->xv, i), ==, ncm_vector_get (s_base->xv, i));
      ncm_assert_cmpdouble (ncm_vector_get (s->yv, i), ==, ncm_vector_get (s_base->yv, i));
    }

    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval (void)
{
  gsl_function F;
  guint i;
  gdouble d[4];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_linear;
    F.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, _NCM_SPLINE_TEST_ERROR);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_cubic;
    F.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, _NCM_SPLINE_TEST_ERROR);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_sin_poly;
    F.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * ((_NCM_SPLINE_TEST_NKNOTS)/100.0 - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS / 100.0 - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, _NCM_SPLINE_TEST_ERROR * 50.0);
    }
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_deriv (void)
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
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_linear;
    F.params = d;
    F_deriv.function = &F_linear_deriv;
    F_deriv.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble df = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, _NCM_SPLINE_TEST_ERROR);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_cubic;
    F.params = d;
    F_deriv.function = &F_cubic_deriv;
    F_deriv.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble df = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, _NCM_SPLINE_TEST_ERROR * 50.0);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_sin_poly;
    F.params = d;
    F_deriv.function = &F_sin_poly_deriv;
    F_deriv.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS / 100.0 - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS / 100.0 - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble df = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, _NCM_SPLINE_TEST_ERROR * 50.0);
    }
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_deriv2 (void)
{
  gsl_function F;
  gsl_function F_deriv2;
  guint i;
  gdouble d[4];
  d[0] = g_test_rand_double ();
  d[1] = g_test_rand_double ();
  d[2] = g_test_rand_double ();
  d[3] = g_test_rand_double ();

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_linear;
    F.params = d;
    F_deriv2.function = &F_linear_deriv2;
    F_deriv2.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble d2f = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, _NCM_SPLINE_TEST_ERROR);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_cubic;
    F.params = d;
    F_deriv2.function = &F_cubic_deriv2;
    F_deriv2.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble d2f = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, _NCM_SPLINE_TEST_ERROR * 50.0);
    }
    ncm_spline_free (s);
  }

  if (FALSE)
  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_sin_poly;
    F.params = d;
    F_deriv2.function = &F_sin_poly_deriv2;
    F_deriv2.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS / 500.0 - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS / 500.0 - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble d2f = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, _NCM_SPLINE_TEST_ERROR * 100.0);
    }
    ncm_spline_free (s);
  }
}

void
test_ncm_spline_eval_int (void)
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
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_linear;
    F.params = d;
    F_int.function = &F_linear_int_0x;
    F_int.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble xi = _NCM_SPLINE_TEST_XI;
      gdouble If = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs = ncm_spline_eval_integ (s, xi, x);

      ncm_assert_cmpdouble_e (Ifs, ==, If, _NCM_SPLINE_TEST_ERROR);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (s_base);
    F.function = &F_cubic;
    F.params = d;
    F_int.function = &F_cubic_int_0x;
    F_int.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, _NCM_SPLINE_TEST_XI, _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1) , _NCM_SPLINE_TEST_NKNOTS, _NCM_SPLINE_TEST_ERROR);
    for (i = 0; i < 2 * _NCM_SPLINE_TEST_NKNOTS; i++)
    {
      gdouble x = _NCM_SPLINE_TEST_XI + _NCM_SPLINE_TEST_DX * (_NCM_SPLINE_TEST_NKNOTS - 1.0) / (2.0 * _NCM_SPLINE_TEST_NKNOTS - 1.0) * i;
      gdouble xi = _NCM_SPLINE_TEST_XI;
      gdouble If = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs = ncm_spline_eval_integ (s, xi, x);

      ncm_assert_cmpdouble_e (Ifs, ==, If, _NCM_SPLINE_TEST_ERROR);
    }
    ncm_spline_free (s);
  }
}

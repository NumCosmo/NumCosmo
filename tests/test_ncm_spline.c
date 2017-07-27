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
  gboolean deriv2;
} TestNcmSpline;

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

void test_ncm_spline_cubic_notaknot_new_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_gsl_cspline_new_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_new (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_new_array (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_new_data (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_copy_empty (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_copy (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_serialize (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_deriv (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_deriv2 (TestNcmSpline *test, gconstpointer pdata);
void test_ncm_spline_eval_int (TestNcmSpline *test, gconstpointer pdata);
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

typedef struct _TestNcmSplineFunc {
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
  {NULL}
};

TestNcmSplineFunc _test_ncm_spline_tests[] = {
  {&test_ncm_spline_new,         "/new"},
  {&test_ncm_spline_new_array,   "/new/array"},
  {&test_ncm_spline_new_data,    "/new/data"},
  {&test_ncm_spline_copy_empty,  "/copy/empty"},
  {&test_ncm_spline_copy,        "/copy"},
  {&test_ncm_spline_serialize,   "/serialize"},
  {&test_ncm_spline_eval,        "/eval"},
  {&test_ncm_spline_eval_deriv,  "/eval/deriv"},
  {&test_ncm_spline_eval_deriv2, "/eval/deriv2"},
  {&test_ncm_spline_eval_int,    "/int"},
  {&test_ncm_spline_traps,       "/traps"},
  {NULL}
};


void 
_test_ncm_spline_add_tests (void (*tnew) (TestNcmSpline *test, gconstpointer pdata), void (*tfree) (TestNcmSpline *test, gconstpointer pdata), const gchar *name)
{
  guint i;
  
#if GLIB_CHECK_VERSION(2,38,0)
  for (i = 0; _test_ncm_spline_traps[i].func != NULL; i++)
  {
    gchar *tname = g_strdup_printf ("/ncm/%s%s", name, _test_ncm_spline_traps[i].name);  
    g_test_add (tname, TestNcmSpline, NULL, tnew, _test_ncm_spline_traps[i].func, tfree);
    g_free (tname);
  }
#endif

  for (i = 0; _test_ncm_spline_tests[i].func != NULL; i++)
  {
    gchar *tname = g_strdup_printf ("/ncm/%s%s", name, _test_ncm_spline_tests[i].name);  
    g_test_add (tname, TestNcmSpline, NULL, tnew, _test_ncm_spline_tests[i].func, tfree);
    g_free (tname);
  }
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  _test_ncm_spline_add_tests (&test_ncm_spline_cubic_notaknot_new_empty,
                              &test_ncm_spline_free_empty,
                              "spline_cubic_notaknot");
  _test_ncm_spline_add_tests (&test_ncm_spline_gsl_cspline_new_empty,
                              &test_ncm_spline_free_empty,
                              "spline_gsl/cspline");  
  g_test_run ();
}

#define _NCM_SPLINE_TEST_NKNOTS 1000
#define _NCM_SPLINE_TEST_DX 0.003

void
test_ncm_spline_cubic_notaknot_new_empty (TestNcmSpline *test, gconstpointer pdata)
{
  test->name     = "spline_cubic_notaknot";
  test->deriv2   = TRUE;
  test->nknots   = g_test_rand_int_range (_NCM_SPLINE_TEST_NKNOTS, 2 * _NCM_SPLINE_TEST_NKNOTS);
  test->dx       = _NCM_SPLINE_TEST_DX;
  test->xi       = 10.0 * GSL_SIGN (g_test_rand_double_range (-1.0, 1.0));
  test->prec     = 1.0e-5;
  test->error    = 1.0e-4;
  test->error_d1 = 5.0e-3;
  test->error_d2 = 1.0e-2;
  test->s_base   = ncm_spline_cubic_notaknot_new ();
  g_assert (NCM_IS_SPLINE_CUBIC (test->s_base));
  g_assert (NCM_IS_SPLINE_CUBIC_NOTAKNOT (test->s_base));

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
  test->deriv2   = FALSE;
  test->nknots   = g_test_rand_int_range (_NCM_SPLINE_TEST_NKNOTS, 2 * _NCM_SPLINE_TEST_NKNOTS);
  test->dx       = _NCM_SPLINE_TEST_DX;
  test->xi       = 10.0 * GSL_SIGN (g_test_rand_double_range (-1.0, 1.0));
  test->prec     = 1.0e-5;
  test->error    = 5.0e-4;
  test->error_d1 = 1.0e-2;
  test->error_d2 = 1.0e-2;
  test->s_base   = ncm_spline_gsl_new (gsl_interp_cspline);
  g_assert (NCM_IS_SPLINE_GSL (test->s_base));
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
  g_assert (NCM_IS_SPLINE (s));
}

void
test_ncm_spline_new (TestNcmSpline *test, gconstpointer pdata)
{
  {
    NcmVector *x = ncm_vector_new (test->nknots);
    NcmVector *y = ncm_vector_new (test->nknots);
    NcmSpline *s = ncm_spline_new (test->s_base, x, y, FALSE);

    test_ncm_spline_new_sanity (s);
    g_assert (s->init == FALSE);
    ncm_spline_free (s);
  }

  {
    NcmVector *x = ncm_vector_new (test->nknots);
    NcmVector *y = ncm_vector_new (test->nknots);
    NcmSpline *s;
    gdouble d[2];
    gint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();
    for (i = 0; i < test->nknots; i++)
    {
      ncm_vector_set (x, i, test->xi + test->dx * i);
      ncm_vector_set (y, i, F_linear (ncm_vector_get (x, i), d));
    }

    s = ncm_spline_new (test->s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert (s->init == TRUE);
    ncm_spline_free (s);
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
    g_assert (s->init == FALSE);
    ncm_spline_free (s);
  }

  {
    GArray *x = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
    GArray *y = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nknots);
    g_array_set_size (x, test->nknots);
    g_array_set_size (y, test->nknots);
    NcmSpline *s;
    gdouble d[2];
    gint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();
    for (i = 0; i < test->nknots; i++)
    {
      g_array_index(x, gdouble, i) = test->xi + test->dx * i;
      g_array_index(y, gdouble, i) = F_linear (g_array_index(x, gdouble, i), d);
    }

    s = ncm_spline_new_array (test->s_base, x, y, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert (s->init == TRUE);
    ncm_spline_free (s);
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
    g_assert (s->init == FALSE);
    ncm_spline_free (s);
  }

  {
    gdouble x[test->nknots];
    gdouble y[test->nknots];
    NcmSpline *s;
    gdouble d[2];
    gint i;

    d[0] = g_test_rand_double ();
    d[1] = g_test_rand_double ();
    for (i = 0; i < test->nknots; i++)
    {
      x[i] = test->xi + test->dx * i;
      y[i] = F_linear (x[i], d);
    }

    s = ncm_spline_new_data (test->s_base, x, y, test->nknots, TRUE);
    test_ncm_spline_new_sanity (s);
    g_assert (s->init == TRUE);
    ncm_spline_free (s);
  }

}

void
test_ncm_spline_copy_empty (TestNcmSpline *test, gconstpointer pdata)
{
  NcmSpline *s = ncm_spline_copy_empty (test->s_base);
  g_assert (G_TYPE_FROM_INSTANCE (s) == G_TYPE_FROM_INSTANCE (test->s_base));
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
    NcmSpline *s = ncm_spline_copy (test->s_base);
    g_assert (s->xv != test->s_base->xv && s->yv != test->s_base->yv);
    for (i = 0; i < test->nknots; i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (s->xv, i), ==, ncm_vector_get (test->s_base->xv, i));
      ncm_assert_cmpdouble (ncm_vector_get (s->yv, i), ==, ncm_vector_get (test->s_base->yv, i));
    }

    ncm_spline_free (s);
  }
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
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcmSpline *s = NCM_SPLINE (ncm_serialize_dup_obj (ser, G_OBJECT (test->s_base)));

    ncm_spline_prepare (s);
    ncm_serialize_free (ser);
    
    g_assert (s->xv != test->s_base->xv && s->yv != test->s_base->yv);
    for (i = 0; i < test->nknots; i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (s->xv, i), ==, ncm_vector_get (test->s_base->xv, i));
      ncm_assert_cmpdouble (ncm_vector_get (s->yv, i), ==, ncm_vector_get (test->s_base->yv, i));
    }

    ncm_spline_free (s);
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
    F.function = &F_linear;
    F.params = d;
    NcmVector *x = ncm_vector_new (test->nknots);
    NcmVector *y = ncm_vector_new (test->nknots);
    NcmSpline *s;
    gint i;

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
    g_assert (s->init == TRUE);

    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble xval = ncm_vector_get (x, 0) + (ncm_vector_get (x, test->nknots - 1) - ncm_vector_get (x, 0)) / (2.0 * test->nknots - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, xval);
      gdouble fs = ncm_spline_eval (s, xval);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }
    
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_linear;
    F.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_cubic;
    F.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_sin_poly;
    F.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * ((test->nknots)/100.0 - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots / 100.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);

      ncm_assert_cmpdouble_e (fs, ==, f, test->error, 0.0);
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
    F.function = &F_linear;
    F.params = d;
    F_deriv.function = &F_linear_deriv;
    F_deriv.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble df = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error, 0.0);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_cubic;
    F.params = d;
    F_deriv.function = &F_cubic_deriv;
    F_deriv.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble df = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error_d1, 0.0);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_sin_poly;
    F.params = d;
    F_deriv.function = &F_sin_poly_deriv;
    F_deriv.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots / 100.0 - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots / 100.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble df = GSL_FN_EVAL (&F_deriv, x);
      gdouble dfs = ncm_spline_eval_deriv (s, x);

      ncm_assert_cmpdouble_e (dfs, ==, df, test->error_d1, 0.0);
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
    F.function = &F_linear;
    F.params = d;
    F_deriv2.function = &F_linear_deriv2;
    F_deriv2.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble d2f = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, test->error_d2, 0.0);
    }
    ncm_spline_free (s);
  }

  if (test->deriv2)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_cubic;
    F.params = d;
    F_deriv2.function = &F_cubic_deriv2;
    F_deriv2.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble d2f = GSL_FN_EVAL (&F_deriv2, x);
      gdouble d2fs = ncm_spline_eval_deriv2 (s, x);

      ncm_assert_cmpdouble_e (d2fs, ==, d2f, test->error_d2, 0.0);
    }
    ncm_spline_free (s);
  }

  if (FALSE)
  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_sin_poly;
    F.params = d;
    F_deriv2.function = &F_sin_poly_deriv2;
    F_deriv2.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots / 500.0 - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots / 500.0 - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble d2f = GSL_FN_EVAL (&F_deriv2, x);
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
    F.function = &F_linear;
    F.params = d;
    F_int.function = &F_linear_int_0x;
    F_int.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble xi = test->xi;
      gdouble If = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs = ncm_spline_eval_integ (s, xi, x);

      ncm_assert_cmpdouble_e (Ifs, ==, If, test->error, 0.0);
    }
    ncm_spline_free (s);
  }

  {
    NcmSpline *s = ncm_spline_copy (test->s_base);
    F.function = &F_cubic;
    F.params = d;
    F_int.function = &F_cubic_int_0x;
    F_int.params = d;
    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, test->xi, test->xi + (test->dx * (test->nknots - 1)) * _TEST_EPSILON, test->nknots, test->prec);
    for (i = 0; i < 2 * test->nknots; i++)
    {
      gdouble x = test->xi + test->dx * (test->nknots - 1.0) / (2.0 * test->nknots - 1.0) * i;
      gdouble xi = test->xi;
      gdouble If = GSL_FN_EVAL (&F_int, x) - GSL_FN_EVAL (&F_int, xi);
      gdouble Ifs = ncm_spline_eval_integ (s, xi, x);

      ncm_assert_cmpdouble_e (Ifs, ==, If, test->error, 0.0);
    }
    ncm_spline_free (s);
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
  GArray *x = NULL;
  GArray *y = NULL;
  NcmSpline *s = ncm_spline_new_array (test->s_base, x, y, FALSE);
  ncm_spline_free (s);
}

void 
test_ncm_spline_traps (TestNcmSpline *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  guint i;
  for (i = 0; _test_ncm_spline_traps[i].func != NULL; i++)
  {
    gchar *tname = g_strdup_printf ("/ncm/%s%s", test->name, _test_ncm_spline_traps[i].name);
    g_test_trap_subprocess (tname, 0, 0);
    g_test_trap_assert_failed ();
    g_free (tname);
  }
#endif
}

/***************************************************************************
 *            test_nc_window.c
 *
 *  Thu May 10 15:13:31 2012
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>


void test_nc_window_new_tophat (void);
void test_nc_window_new_gaussian (void);
void test_nc_window_eval_fourier_tophat (void);
void test_nc_window_eval_fourier_gaussian (void);
void test_nc_window_deriv_fourier_tophat (void);
void test_nc_window_deriv_fourier_gaussian (void);
void test_nc_window_eval_real_tophat (void); 
void test_nc_window_eval_real_gaussian (void);
void test_nc_window_free (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/window/tophat/new", &test_nc_window_new_tophat);
  g_test_add_func ("/nc/window/tophat/eval_fourier", &test_nc_window_eval_fourier_tophat);
  g_test_add_func ("/nc/window/tophat/deriv_fourier", &test_nc_window_deriv_fourier_tophat);
  g_test_add_func ("/nc/window/tophat/eval_real", &test_nc_window_eval_real_tophat);
  g_test_add_func ("/nc/window/tophat/free", &test_nc_window_free);

  g_test_add_func ("/nc/window/gaussian/new", &test_nc_window_new_gaussian);
  g_test_add_func ("/nc/window/gaussian/eval_fourier", &test_nc_window_eval_fourier_gaussian);
  g_test_add_func ("/nc/window/gaussian/deriv_fourier", &test_nc_window_deriv_fourier_gaussian);
  g_test_add_func ("/nc/window/gaussian/eval_real", &test_nc_window_eval_real_gaussian);
  g_test_add_func ("/nc/window/gaussian/free", &test_nc_window_free); 

  g_test_run ();
}

#define _NC_WINDOW_TEST_KMIN 0.01
#define _NC_WINDOW_TEST_KMAX 1.0e3
#define _NC_WINDOW_TEST_XMIN 0.01
#define _NC_WINDOW_TEST_XMAX 1.0e3

NcWindow *wf;

void
test_nc_window_new_tophat (void)
{
  wf = nc_window_tophat_new ();
  g_assert (NC_IS_WINDOW (wf));
  g_assert (NC_IS_WINDOW_TOPHAT (wf));

  test_nc_window_free ();

#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 30))
  wf = nc_window_new_from_name ("NcWindowTophat");
#else
  wf = nc_window_tophat_new ();
#endif
  g_assert (NC_IS_WINDOW (wf));
  g_assert (NC_IS_WINDOW_TOPHAT (wf));  
}

void
test_nc_window_new_gaussian (void)
{
  wf = nc_window_gaussian_new ();
  g_assert (NC_IS_WINDOW (wf));
  g_assert (NC_IS_WINDOW_GAUSSIAN (wf));

  test_nc_window_free();

#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 30))
  wf = nc_window_new_from_name ("NcWindowGaussian");
#else
  wf = nc_window_gaussian_new ();
#endif
  g_assert (NC_IS_WINDOW (wf));
  g_assert (NC_IS_WINDOW_GAUSSIAN (wf));  
}

void
test_nc_window_free (void)
{
  NCM_TEST_FREE (nc_window_free, wf);
}

void
test_nc_window_eval_fourier_tophat (void)
{
  guint i, j;

  {
    gdouble w1 = nc_window_eval_fourier (wf, 0.0, 5.0);
    ncm_assert_cmpdouble (w1, ==, 1.0);
  }

  for (j = 0; j < 30; j++)
  {
    gdouble R = 2.0 + 15.0 / 29.0 * j;
    for (i = 0; i < 50; i++)
    {
      gdouble k = _NC_WINDOW_TEST_KMIN + (_NC_WINDOW_TEST_KMAX - _NC_WINDOW_TEST_KMIN) / (49.0) * i;
      gdouble kR = k * R;
      gdouble w1 = nc_window_eval_fourier (wf, k, R);
      gdouble w2 = 3.0 * gsl_sf_bessel_j1 (kR) / kR;

      ncm_assert_cmpdouble_e (w1, ==, w2, 1e-8);
    }
  }
}

void
test_nc_window_deriv_fourier_tophat (void)
{
  guint i, j;

  for (j = 0; j < 30; j++)
  {
    gdouble R = 2.0 + 15.0 / 29.0 * j;
    for (i = 0; i < 50; i++)
    {
      gdouble k = _NC_WINDOW_TEST_KMIN + (_NC_WINDOW_TEST_KMAX - _NC_WINDOW_TEST_KMIN) / (49.0) * i;
      gdouble w1 = nc_window_deriv_fourier (wf, k, R);
      gdouble w2 = -3.0 * gsl_sf_bessel_j2 (k * R) / R;

      ncm_assert_cmpdouble_e (w1, ==, w2, 1e-11);
    }
  }
}

void
test_nc_window_eval_real_tophat (void)
{
  guint i, j;

  for (j = 0; j < 10; j++)
  {
    gdouble R = 2.0 + 15.0 / 9.0 * j;
    gdouble R3 = R * R * R;
    for (i = 0; i < 10; i++)
    {
      gdouble r = _NC_WINDOW_TEST_XMIN + (_NC_WINDOW_TEST_XMAX - _NC_WINDOW_TEST_XMIN) / (9.0) * i;
      gdouble w1 = nc_window_eval_realspace (wf, r, R);
      gdouble w2;
      if (r <= R)
        w2 = 3.0 / (4.0 * M_PI * R3);
      else
        w2 = 0.0;

      ncm_assert_cmpdouble_e (w1, ==, w2, 1e-11);
    }
  }
}

void
test_nc_window_eval_fourier_gaussian (void)
{
  guint i, j;

  for (j = 0; j < 30; j++)
  {
    gdouble R = 2.0 + 15.0 / 29.0 * j;
    for (i = 0; i < 50; i++)
    {
      gdouble k = _NC_WINDOW_TEST_KMIN + (_NC_WINDOW_TEST_KMAX - _NC_WINDOW_TEST_KMIN) / (49.0) * i;
      gdouble kR = k * R;
      gdouble kR2 = kR * kR;
      gdouble w1 = nc_window_eval_fourier (wf, k, R);
      gdouble w2 = exp (-kR2 / 2.0);

      ncm_assert_cmpdouble_e (w1, ==, w2, 1e-11);
    }
  }
}

void
test_nc_window_deriv_fourier_gaussian (void)
{
  guint i, j;

  for (j = 0; j < 20; j++)
  {
    gdouble R = 2.0 + 15.0 / 19.0 * j;
    for (i = 0; i < 20; i++)
    {
      gdouble k = _NC_WINDOW_TEST_KMIN + (_NC_WINDOW_TEST_KMAX - _NC_WINDOW_TEST_KMIN) / (19.0) * i;
      gdouble k2R = k * k * R; 
      gdouble kR2 = k2R * R;
      gdouble w1 = nc_window_deriv_fourier (wf, k, R);
      gdouble w2 = -k2R * exp (-kR2 / 2.0);

      ncm_assert_cmpdouble_e (w1, ==, w2, 1.0e-15);
    }
  }
}

void
test_nc_window_eval_real_gaussian (void)
{
  guint i, j;

  for (j = 0; j < 30; j++)
  {
    gdouble R = 2.0 + 15.0 / 29.0 * j;
    for (i = 0; i < 50; i++)
    {
      gdouble r = _NC_WINDOW_TEST_XMIN + (_NC_WINDOW_TEST_XMAX - _NC_WINDOW_TEST_XMIN) / (49.0) * i;
      gdouble r_R2 = r * r / (R * R);
      gdouble w1 = nc_window_eval_realspace (wf, r, R);
      gdouble w2 = 1.0 / gsl_pow_3 (sqrt (2 * M_PI * R * R)) * exp (-r_R2 / 2.0);

      ncm_assert_cmpdouble_e (w1, ==, w2, 1e-11);
    }
  }
}

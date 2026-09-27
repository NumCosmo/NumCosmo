/***************************************************************************
 *            test_ncm_sf_spherical_harmonics.c
 *
 *  Sun January 07 20:34:52 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <vitenti@uel.br>
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

#include <gsl/gsl_sf_legendre.h>

typedef struct _TestNcmSFSphericalHarmonics
{
  NcmSFSphericalHarmonics *spha;
  guint ntests;
} TestNcmSFSphericalHarmonics;

/*
 * Errors are scaled by the peak of |Y_l^m| over l at fixed m: relative errors are
 * meaningless at the zeros. The angles stay 0.01 away from the poles, where the step in m
 * after skipped orders loses accuracy (under review).
 */
#ifndef TEST_TOL
#define TEST_TOL (1.0e-11)
#endif
#define TEST_THETA_MIN (0.01)

static gdouble _test_max_err = 0.0;

static void test_ncm_sf_spherical_harmonics_new (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_free (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);

static void test_ncm_sf_spherical_harmonics_single_rec (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_rec2 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_rec4 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_recn (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);

static void test_ncm_sf_spherical_harmonics_array_single_rec (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_array_rec2 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_array_rec4 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_array_recn (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);

static void test_ncm_sf_spherical_harmonics_lmax (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);

static void test_ncm_sf_spherical_harmonics_traps (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);
static void test_ncm_sf_spherical_harmonics_invalid_test (TestNcmSFSphericalHarmonics *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/sf/spherical_harmonics/single_rec", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_single_rec,
              &test_ncm_sf_spherical_harmonics_free);
  g_test_add ("/ncm/sf/spherical_harmonics/rec2", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_rec2,
              &test_ncm_sf_spherical_harmonics_free);
  g_test_add ("/ncm/sf/spherical_harmonics/rec4", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_rec4,
              &test_ncm_sf_spherical_harmonics_free);
  g_test_add ("/ncm/sf/spherical_harmonics/recn", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_recn,
              &test_ncm_sf_spherical_harmonics_free);

  g_test_add ("/ncm/sf/spherical_harmonics/array/single_rec", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_array_single_rec,
              &test_ncm_sf_spherical_harmonics_free);
  g_test_add ("/ncm/sf/spherical_harmonics/array/rec2", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_array_rec2,
              &test_ncm_sf_spherical_harmonics_free);
  g_test_add ("/ncm/sf/spherical_harmonics/array/rec4", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_array_rec4,
              &test_ncm_sf_spherical_harmonics_free);
  g_test_add ("/ncm/sf/spherical_harmonics/array/recn", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_array_recn,
              &test_ncm_sf_spherical_harmonics_free);

  g_test_add ("/ncm/sf/spherical_harmonics/lmax", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_lmax,
              &test_ncm_sf_spherical_harmonics_free);

  g_test_add ("/ncm/sf/spherical_harmonics/traps", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_traps,
              &test_ncm_sf_spherical_harmonics_free);

  g_test_add ("/ncm/sf/spherical_harmonics/invalid/test/subprocess", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_invalid_test,
              &test_ncm_sf_spherical_harmonics_free);

  g_test_run ();

  g_test_message ("largest peak-scaled error: %.3e", _test_max_err);
}

static void
test_ncm_sf_spherical_harmonics_new (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint lmax = g_test_rand_int_range (512, 2048);

  test->spha = ncm_sf_spherical_harmonics_new (lmax);
}

static void
test_ncm_sf_spherical_harmonics_free (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha = test->spha;

  NCM_TEST_FREE (ncm_sf_spherical_harmonics_free, spha);
}

static gdouble *
_test_peaks (const gdouble *Yblm, const guint lmax)
{
  gdouble *peak = g_new0 (gdouble, lmax + 1);
  guint l, m;

  for (m = 0; m <= lmax; m++)
    for (l = m; l <= lmax; l++)
      peak[m] = GSL_MAX (peak[m], fabs (Yblm[gsl_sf_legendre_array_index (l, m)]));

  return peak;
}

static gboolean
_test_fails (const gdouble Ygsl, const gdouble Ync, const gdouble peak)
{
  const gdouble err = fabs (Ync - Ygsl) / peak;

  _test_max_err = GSL_MAX (_test_max_err, err);

  return err > TEST_TOL;
}

static void
test_ncm_sf_spherical_harmonics_single_rec (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha   = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta             = g_test_rand_double_range (TEST_THETA_MIN, M_PI - TEST_THETA_MIN);
  const gdouble x                 = cos (theta);
  const guint lmax                = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize               = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm                   = g_new (gdouble, asize);
  gdouble *peak;
  guint nerr = 0;
  gint l, m;

  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);

  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  peak = _test_peaks (Yblm, lmax);

  m = 0;
  l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

  while (TRUE)
  {
    while (TRUE)
    {
      gsize lm_index = gsl_sf_legendre_array_index (l, m);

      if (_test_fails (Yblm[lm_index], ncm_sf_spherical_harmonics_Y_get_lm (sphaY), peak[m]))
        nerr++;

      if (l < (gint) lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_l (sphaY);
        l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      }
      else
      {
        break;
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  g_free (Yblm);
  g_free (peak);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_rec2 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha   = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta             = g_test_rand_double_range (TEST_THETA_MIN, M_PI - TEST_THETA_MIN);
  const gdouble x                 = cos (theta);
  const guint lmax                = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize               = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm                   = g_new (gdouble, asize);
  gdouble *peak;
  gdouble Ylm[2];
  gint l, m, nerr = 0;

  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);

  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  peak = _test_peaks (Yblm, lmax);

  m = 0;

  while (TRUE)
  {
    while (TRUE)
    {
      gint j;

      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l + 2 > (gint) lmax)
        break;

      ncm_sf_spherical_harmonics_Y_next_l2 (sphaY, Ylm);

      for (j = 0; j < 2; j++)
      {
        if (_test_fails (Yblm[gsl_sf_legendre_array_index (l + j, m)], Ylm[j], peak[m]))
          nerr++;
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  g_free (Yblm);
  g_free (peak);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_rec4 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha   = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta             = g_test_rand_double_range (TEST_THETA_MIN, M_PI - TEST_THETA_MIN);
  const gdouble x                 = cos (theta);
  const guint lmax                = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize               = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm                   = g_new (gdouble, asize);
  gdouble *peak;
  gdouble Ylm[4];
  gint l, m, nerr = 0;

  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);

  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  peak = _test_peaks (Yblm, lmax);

  m = 0;

  while (TRUE)
  {
    while (TRUE)
    {
      gint j;

      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l + 4 > (gint) lmax)
        break;

      ncm_sf_spherical_harmonics_Y_next_l4 (sphaY, Ylm);

      for (j = 0; j < 4; j++)
      {
        if (_test_fails (Yblm[gsl_sf_legendre_array_index (l + j, m)], Ylm[j], peak[m]))
          nerr++;
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  g_free (Yblm);
  g_free (peak);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_recn (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint n                   = g_test_rand_int_range (3, 10);
  NcmSFSphericalHarmonics *spha   = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta             = g_test_rand_double_range (TEST_THETA_MIN, M_PI - TEST_THETA_MIN);
  const gdouble x                 = cos (theta);
  const guint lmax                = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize               = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm                   = g_new (gdouble, asize);
  gdouble *peak;
  gdouble *Ylm = g_new (gdouble, n + 2);
  guint nerr   = 0;
  gint l, m;

  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);

  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  peak = _test_peaks (Yblm, lmax);

  m = 0;

  while (TRUE)
  {
    while (TRUE)
    {
      guint j;

      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l + n + 2 > lmax)
        break;

      ncm_sf_spherical_harmonics_Y_next_l2pn (sphaY, Ylm, n);

      for (j = 0; j < n + 2; j++)
      {
        if (_test_fails (Yblm[gsl_sf_legendre_array_index (l + j, m)], Ylm[j], peak[m]))
          nerr++;
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  g_free (Yblm);
  g_free (Ylm);
  g_free (peak);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_array_single_rec (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint len                       = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha         = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta                        = g_new (gdouble, len);
  const guint lmax                      = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize                     = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm                        = g_new (gdouble *, len);
  gdouble **peak                        = g_new (gdouble *, len);
  const gdouble theta_b                 = g_test_rand_double_range (TEST_THETA_MIN / 0.9, (M_PI - TEST_THETA_MIN) / 1.1);
  guint nerr                            = 0;
  gint l, m;
  guint i;

  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);

    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
    peak[i] = _test_peaks (Yblm[i], lmax);
  }

  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);

  m = 0;
  l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

  while (TRUE)
  {
    while (TRUE)
    {
      gsize lm_index = gsl_sf_legendre_array_index (l, m);

      for (i = 0; i < len; i++)
      {
        const gdouble Ygsl = Yblm[i][lm_index];
        const gdouble Ync  = ncm_sf_spherical_harmonics_Y_array_get_lm (sphaYa, len, i);


        if (_test_fails (Ygsl, Ync, peak[i][m]))
          nerr++;
      }

      if (l < (gint) lmax)
      {
        ncm_sf_spherical_harmonics_Y_array_next_l (sphaYa, len);
        l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      }
      else
      {
        break;
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
    g_free (peak[i]);
  }

  g_free (Yblm);
  g_free (peak);
  g_free (theta);

  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_array_rec2 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint len                       = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha         = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta                        = g_new (gdouble, len);
  const guint lmax                      = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize                     = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm                        = g_new (gdouble *, len);
  gdouble **peak                        = g_new (gdouble *, len);
  gdouble *Ylm                          = g_new (gdouble, len * 2);
  const gdouble theta_b                 = g_test_rand_double_range (TEST_THETA_MIN / 0.9, (M_PI - TEST_THETA_MIN) / 1.1);
  guint nerr                            = 0;
  gint l, m;
  guint i;

  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);

    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
    peak[i] = _test_peaks (Yblm[i], lmax);
  }

  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);

  m = 0;

  while (TRUE)
  {
    while (TRUE)
    {
      guint j;

      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l + 2 > (gint) lmax)
        break;

      ncm_sf_spherical_harmonics_Y_array_next_l2 (sphaYa, len, Ylm);

      for (j = 0; j < 2; j++)
      {
        gsize lm_index = gsl_sf_legendre_array_index (l + j, m);

        for (i = 0; i < len; i++)
        {
          const gdouble Ygsl = Yblm[i][lm_index];
          const gdouble Ync  = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j, len)];


          if (_test_fails (Ygsl, Ync, peak[i][m]))
            nerr++;
        }
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
    g_free (peak[i]);
  }

  g_free (Yblm);
  g_free (peak);
  g_free (Ylm);
  g_free (theta);

  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_array_rec4 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint len                       = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha         = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta                        = g_new (gdouble, len);
  const guint lmax                      = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize                     = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm                        = g_new (gdouble *, len);
  gdouble **peak                        = g_new (gdouble *, len);
  gdouble *Ylm                          = g_new (gdouble, len * 4);
  const gdouble theta_b                 = g_test_rand_double_range (TEST_THETA_MIN / 0.9, (M_PI - TEST_THETA_MIN) / 1.1);
  guint nerr                            = 0;
  gint l, m;
  guint i;

  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);

    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
    peak[i] = _test_peaks (Yblm[i], lmax);
  }

  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);

  m = 0;

  while (TRUE)
  {
    while (TRUE)
    {
      gint j;

      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l + 4 > (gint) lmax)
        break;

      ncm_sf_spherical_harmonics_Y_array_next_l4 (sphaYa, len, Ylm);

      for (j = 0; j < 4; j++)
      {
        gsize lm_index = gsl_sf_legendre_array_index (l + j, m);

        for (i = 0; i < len; i++)
        {
          const gdouble Ygsl = Yblm[i][lm_index];
          const gdouble Ync  = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j, len)];


          if (_test_fails (Ygsl, Ync, peak[i][m]))
            nerr++;
        }
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
    g_free (peak[i]);
  }

  g_free (Yblm);
  g_free (peak);
  g_free (Ylm);
  g_free (theta);

  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_array_recn (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint n                         = g_test_rand_int_range (3, 10);
  const guint len                       = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha         = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta                        = g_new (gdouble, len);
  const guint lmax                      = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize                     = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm                        = g_new (gdouble *, len);
  gdouble **peak                        = g_new (gdouble *, len);
  gdouble *Ylm                          = g_new (gdouble, len * (n + 2));
  const gdouble theta_b                 = g_test_rand_double_range (TEST_THETA_MIN / 0.9, (M_PI - TEST_THETA_MIN) / 1.1);
  guint nerr                            = 0;
  gint l, m;
  guint i;

  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);

    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
    peak[i] = _test_peaks (Yblm[i], lmax);
  }

  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);

  m = 0;

  while (TRUE)
  {
    while (TRUE)
    {
      guint j;

      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l + n + 2 > lmax)
        break;

      ncm_sf_spherical_harmonics_Y_array_next_l2pn (sphaYa, len, Ylm, n);

      for (j = 0; j < n + 2; j++)
      {
        gsize lm_index = gsl_sf_legendre_array_index (l + j, m);

        for (i = 0; i < len; i++)
        {
          const gdouble Ygsl = Yblm[i][lm_index];
          const gdouble Ync  = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j, len)];


          if (_test_fails (Ygsl, Ync, peak[i][m]))
            nerr++;
        }
      }
    }

    if (m < (gint) lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

      if (l > (gint) lmax)
        break;
    }
    else
    {
      break;
    }
  }

  if (nerr > 0)
    g_error ("%u values off by more than %g of the peak, lmax %u.", nerr, TEST_TOL, lmax);

  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
    g_free (peak[i]);
  }

  g_free (Yblm);
  g_free (peak);
  g_free (Ylm);
  g_free (theta);

  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

/* Walks every (l, m) up to lmax at theta, in the order of the recursion */
static GArray *
_test_walk (NcmSFSphericalHarmonics *spha, const gdouble theta)
{
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gint lmax                 = ncm_sf_spherical_harmonics_get_lmax (spha);
  GArray *vals                    = g_array_new (FALSE, FALSE, sizeof (gdouble));

  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);

  while (TRUE)
  {
    while (TRUE)
    {
      const gdouble Ylm = ncm_sf_spherical_harmonics_Y_get_lm (sphaY);

      g_array_append_val (vals, Ylm);

      if (ncm_sf_spherical_harmonics_Y_get_l (sphaY) >= lmax)
        break;

      ncm_sf_spherical_harmonics_Y_next_l (sphaY);
    }

    if (ncm_sf_spherical_harmonics_Y_get_m (sphaY) >= lmax)
      break;

    ncm_sf_spherical_harmonics_Y_next_m (sphaY);

    if (ncm_sf_spherical_harmonics_Y_get_l (sphaY) > lmax)
      break;
  }

  ncm_sf_spherical_harmonics_Y_free (sphaY);

  return vals;
}

static void
test_ncm_sf_spherical_harmonics_lmax (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const gdouble theta             = 0.7;
  NcmSFSphericalHarmonics *spha   = ncm_sf_spherical_harmonics_new (0);
  NcmSFSphericalHarmonics *fresh  = ncm_sf_spherical_harmonics_new (60);
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  GArray *a, *b;
  guint i;

  /* lmax = 0 builds its tables: Y_0^0 and Y_1^0 */
  g_assert_cmpuint (ncm_sf_spherical_harmonics_get_lmax (spha), ==, 0);
  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);
  ncm_assert_cmpdouble_e (ncm_sf_spherical_harmonics_Y_get_lm (sphaY), ==, 0.5 / sqrt (M_PI), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_sf_spherical_harmonics_Y_get_lp1m (sphaY), ==, sqrt (3.0 / (4.0 * M_PI)) * cos (theta), 1.0e-15, 0.0);
  ncm_sf_spherical_harmonics_Y_free (sphaY);

  /* Growing and shrinking keeps the coefficients a fresh object computes */
  ncm_sf_spherical_harmonics_set_lmax (spha, 50);
  ncm_sf_spherical_harmonics_set_lmax (spha, 10);
  ncm_sf_spherical_harmonics_set_lmax (spha, 60);

  a = _test_walk (spha, theta);
  b = _test_walk (fresh, theta);

  g_assert_cmpuint (a->len, ==, b->len);
  g_assert_cmpuint (a->len, ==, 61 * 62 / 2);

  for (i = 0; i < a->len; i++)
    g_assert_cmpfloat (g_array_index (a, gdouble, i), ==, g_array_index (b, gdouble, i));

  g_array_unref (a);
  g_array_unref (b);
  ncm_sf_spherical_harmonics_free (spha);
  ncm_sf_spherical_harmonics_free (fresh);
}

static void
test_ncm_sf_spherical_harmonics_traps (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/sf/spherical_harmonics/invalid/test/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

static void
test_ncm_sf_spherical_harmonics_invalid_test (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}


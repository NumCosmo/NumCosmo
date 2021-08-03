/***************************************************************************
 *            test_ncm_sf_spherical_harmonics.c
 *
 *  Sun January 07 20:34:52 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <sandro@isoftware.com.br>
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

#define TEST_RELTOL (1.0e-7)

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
  
  g_test_add ("/ncm/sf/spherical_harmonics/traps", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_traps,
              &test_ncm_sf_spherical_harmonics_free);
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/sf/spherical_harmonics/invalid/test/subprocess", TestNcmSFSphericalHarmonics, NULL,
              &test_ncm_sf_spherical_harmonics_new,
              &test_ncm_sf_spherical_harmonics_invalid_test,
              &test_ncm_sf_spherical_harmonics_free);
#endif
  
  g_test_run ();
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

static void
test_ncm_sf_spherical_harmonics_single_rec (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta = g_test_rand_double_range (0.0, M_PI);
  const gdouble x = cos (theta);
  const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm = g_new (gdouble, asize);
  gint l, m, nerr = 0;
  
  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);
  
  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  
  m = 0;
  l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  
  while (TRUE)
  {
    while (TRUE)
    {
      gsize lm_index = gsl_sf_legendre_array_index (l, m);
      
      /*printf ("%6d %6d %6d % 22.15g % 22.15g % 22.15g %e\n", lmax, l, m, theta, Yblm[lm_index], ncm_sf_spherical_harmonics_Y_get_lm (sphaY), fabs (ncm_sf_spherical_harmonics_Y_get_lm (sphaY) / Yblm[lm_index] - 1.0));*/
      if (ncm_cmp (Yblm[lm_index], ncm_sf_spherical_harmonics_Y_get_lm (sphaY), TEST_RELTOL, 0.0) != 0)
        nerr++;
      
      if (l < lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_l (sphaY);
        l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      }
      else
      {
        break;
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  g_free (Yblm);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_rec2 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha   = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta             = g_test_rand_double_range (0.0, M_PI);
  const gdouble x                 = cos (theta);
  const guint lmax                = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize               = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm                   = g_new (gdouble, asize);
  gdouble Ylm[2];
  gint l, m, nerr = 0;
  
  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);
  
  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  
  m = 0;
  
  while (TRUE)
  {
    while (TRUE)
    {
      gint j;
      
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      
      if (l + 2 > lmax)
        break;
      
      ncm_sf_spherical_harmonics_Y_next_l2 (sphaY, Ylm);
      
      for (j = 0; j < 2; j++)
      {
        if (ncm_cmp (Yblm[gsl_sf_legendre_array_index (l + j, m)], Ylm[j], TEST_RELTOL, 0.0) != 0)
          nerr++;
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  g_free (Yblm);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_rec4 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  NcmSFSphericalHarmonics *spha   = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta             = g_test_rand_double_range (0.0, M_PI);
  const gdouble x                 = cos (theta);
  const guint lmax                = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize               = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm                   = g_new (gdouble, asize);
  gdouble Ylm[4];
  gint l, m, nerr = 0;
  
  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);
  
  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  
  m = 0;
  
  while (TRUE)
  {
    while (TRUE)
    {
      gint j;
      
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      
      if (l + 4 > lmax)
        break;
      
      ncm_sf_spherical_harmonics_Y_next_l4 (sphaY, Ylm);
      
      for (j = 0; j < 4; j++)
      {
        if (ncm_cmp (Yblm[gsl_sf_legendre_array_index (l + j, m)], Ylm[j], TEST_RELTOL, 0.0) != 0)
          nerr++;
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  g_free (Yblm);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_recn (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint n = g_test_rand_int_range (3, 10);
  NcmSFSphericalHarmonics *spha = test->spha;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gdouble theta = g_test_rand_double_range (0.0, M_PI);
  const gdouble x = cos (theta);
  const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize = gsl_sf_legendre_array_n (lmax);
  gdouble *Yblm = g_new (gdouble, asize);
  gdouble *Ylm = g_new (gdouble, n + 2);
  gint l, m, nerr = 0;
  
  ncm_sf_spherical_harmonics_start_rec (spha, sphaY, theta);
  
  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Yblm);
  
  m = 0;
  
  while (TRUE)
  {
    while (TRUE)
    {
      gint j;
      
      l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      
      if (l + n + 2 > lmax)
        break;
      
      ncm_sf_spherical_harmonics_Y_next_l2pn (sphaY, Ylm, n);
      
      for (j = 0; j < n + 2; j++)
      {
        if (ncm_cmp (Yblm[gsl_sf_legendre_array_index (l + j, m)], Ylm[j], TEST_RELTOL, 0.0) != 0)
          nerr++;
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  g_free (Yblm);
  g_free (Ylm);
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
test_ncm_sf_spherical_harmonics_array_single_rec (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint len = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta = g_new (gdouble, len);
  const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm = g_new (gdouble *, len);
  const gdouble theta_b = g_test_rand_double_range (0.0, M_PI);
  gint l, m, nerr = 0, i;
  
  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);
    
    if (theta[i] > M_PI)
      theta[i] -= M_PI;
    
    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
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
        
        /*printf ("[%d] %6d %6d %6d % 22.15g % 22.15g % 22.15g %e\n", i, lmax, l, m, theta[i] / M_PI, Ygsl, Ync, fabs (Ync / Ygsl - Ygsl / Ync));*/
        
        if (ncm_cmp (Ygsl, Ync, TEST_RELTOL, 0.0) != 0)
          nerr++;
      }
      
      if (l < lmax)
      {
        ncm_sf_spherical_harmonics_Y_array_next_l (sphaYa, len);
        l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      }
      else
      {
        break;
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
  }
  
  g_free (Yblm);
  
  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_array_rec2 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint len = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta = g_new (gdouble, len);
  const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm = g_new (gdouble *, len);
  gdouble *Ylm = g_new (gdouble, len * 2);
  const gdouble theta_b = g_test_rand_double_range (0.0, M_PI);
  gint l, m, nerr = 0, i;
  
  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);
    
    if (theta[i] > M_PI)
      theta[i] -= M_PI;
    
    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
  }
  
  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);
  
  m = 0;
  
  while (TRUE)
  {
    while (TRUE)
    {
      gint j;
      
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      
      if (l + 2 > lmax)
        break;
      
      ncm_sf_spherical_harmonics_Y_array_next_l2 (sphaYa, len, Ylm);
      
      for (j = 0; j < 2; j++)
      {
        gsize lm_index = gsl_sf_legendre_array_index (l + j, m);
        
        for (i = 0; i < len; i++)
        {
          const gdouble Ygsl = Yblm[i][lm_index];
          const gdouble Ync  = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j, len)];
          
          /*printf ("[%d] %6d %6d %6d % 22.15g % 22.15g % 22.15g %e\n", i, lmax, l, m, theta[i] / M_PI, Ygsl, Ync, fabs (Ync / Ygsl - Ygsl / Ync));*/
          
          if (ncm_cmp (Ygsl, Ync, TEST_RELTOL, 0.0) != 0)
            nerr++;
        }
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
  }
  
  g_free (Yblm);
  g_free (Ylm);
  
  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_array_rec4 (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint len = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta = g_new (gdouble, len);
  const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm = g_new (gdouble *, len);
  gdouble *Ylm = g_new (gdouble, len * 4);
  const gdouble theta_b = g_test_rand_double_range (0.0, M_PI);
  gint l, m, nerr = 0, i;
  
  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);
    
    if (theta[i] > M_PI)
      theta[i] -= M_PI;
    
    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
  }
  
  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);
  
  m = 0;
  
  while (TRUE)
  {
    while (TRUE)
    {
      gint j;
      
      l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      
      if (l + 4 > lmax)
        break;
      
      ncm_sf_spherical_harmonics_Y_array_next_l4 (sphaYa, len, Ylm);
      
      for (j = 0; j < 4; j++)
      {
        gsize lm_index = gsl_sf_legendre_array_index (l + j, m);
        
        for (i = 0; i < len; i++)
        {
          const gdouble Ygsl = Yblm[i][lm_index];
          const gdouble Ync  = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j, len)];
          
          /*printf ("[%d] %6d %6d %6d % 22.15g % 22.15g % 22.15g %e\n", i, lmax, l, m, theta[i] / M_PI, Ygsl, Ync, fabs (Ync / Ygsl - Ygsl / Ync));*/
          
          if (ncm_cmp (Ygsl, Ync, TEST_RELTOL, 0.0) != 0)
            nerr++;
        }
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
  }
  
  g_free (Yblm);
  g_free (Ylm);
  
  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_array_recn (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  const guint n = g_test_rand_int_range (3, 10);
  const guint len = g_test_rand_int_range (2, 6);
  NcmSFSphericalHarmonics *spha = test->spha;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (spha, len, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gdouble *theta = g_new (gdouble, len);
  const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
  const guint asize = gsl_sf_legendre_array_n (lmax);
  gdouble **Yblm = g_new (gdouble *, len);
  gdouble *Ylm = g_new (gdouble, len * (n + 2));
  const gdouble theta_b = g_test_rand_double_range (0.0, M_PI);
  gint l, m, nerr = 0, i;
  
  for (i = 0; i < len; i++)
  {
    Yblm[i]  = g_new (gdouble, asize);
    theta[i] = theta_b * g_test_rand_double_range (0.90, 1.1);
    
    if (theta[i] > M_PI)
      theta[i] -= M_PI;
    
    gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, cos (theta[i]), -1.0, Yblm[i]);
  }
  
  ncm_sf_spherical_harmonics_start_rec_array (spha, sphaYa, len, theta);
  
  m = 0;
  
  while (TRUE)
  {
    while (TRUE)
    {
      gint j;
      
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
          
          /*printf ("[%d] %6d %6d %6d % 22.15g % 22.15g % 22.15g %e\n", i, lmax, l, m, theta[i] / M_PI, Ygsl, Ync, fabs (Ync / Ygsl - Ygsl / Ync));*/
          
          if (ncm_cmp (Ygsl, Ync, TEST_RELTOL, 0.0) != 0)
            nerr++;
        }
      }
    }
    
    if (m < lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, len);
      m = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      
      if (l > lmax)
        break;
    }
    else
    {
      break;
    }
  }
  
  if (nerr > 5)
    g_error ("More than 5 failures `%d', lmax `%d'.", nerr, lmax);
  
  for (i = 0; i < len; i++)
  {
    g_free (Yblm[i]);
  }
  
  g_free (Yblm);
  g_free (Ylm);
  
  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
test_ncm_sf_spherical_harmonics_traps (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/sf/spherical_harmonics/invalid/test/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

static void
test_ncm_sf_spherical_harmonics_invalid_test (TestNcmSFSphericalHarmonics *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}


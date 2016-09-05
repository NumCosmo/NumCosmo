/***************************************************************************
 *            test_ncm_sphere_map_pix.c
 *
 *  Sun July 17 17:02:50 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <sandro@isoftware.com.br>
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

typedef struct _TestNcmSphereMapPix
{
  NcmSphereMapPix *pix;
  guint nside;
  gboolean try;
} TestNcmSphereMapPix;

void test_ncm_sphere_map_pix_new (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_free (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_sanity (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_angles (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_ring (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_pix2alm (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_pix2alm2pix (TestNcmSphereMapPix *test, gconstpointer pdata);

void test_ncm_sphere_map_pix_traps (TestNcmSphereMapPix *test, gconstpointer pdata);
void test_ncm_sphere_map_pix_invalid_nside (TestNcmSphereMapPix *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/sphere_map_pix/sanity", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_sanity,
              &test_ncm_sphere_map_pix_free);

  g_test_add ("/ncm/sphere_map_pix/angles", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_angles,
              &test_ncm_sphere_map_pix_free);

  g_test_add ("/ncm/sphere_map_pix/ring", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_ring,
              &test_ncm_sphere_map_pix_free);
  
  g_test_add ("/ncm/sphere_map_pix/pix2alm", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_pix2alm,
              &test_ncm_sphere_map_pix_free);

  g_test_add ("/ncm/sphere_map_pix/pix2alm2pix", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_pix2alm2pix,
              &test_ncm_sphere_map_pix_free);

  g_test_add ("/ncm/sphere_map_pix/traps", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_traps,
              &test_ncm_sphere_map_pix_free);

#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_add ("/ncm/sphere_map_pix/invalid/nside/subprocess", TestNcmSphereMapPix, NULL,
              &test_ncm_sphere_map_pix_new,
              &test_ncm_sphere_map_pix_invalid_nside,
              &test_ncm_sphere_map_pix_free);
#endif
  g_test_run ();
}

void
test_ncm_sphere_map_pix_new (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  test->nside = 64;//1 << g_test_rand_int_range (1, 8);
  
  test->pix   = ncm_sphere_map_pix_new (test->nside);
  test->try   = FALSE;

  g_assert (test->pix != NULL);
  g_assert (NCM_IS_SPHERE_MAP_PIX (test->pix));

  g_assert_cmpint (ncm_sphere_map_pix_get_nside (test->pix), ==, test->nside);
}

void
test_ncm_sphere_map_pix_free (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_sphere_map_pix_free, test->pix);
}

void
test_ncm_sphere_map_pix_sanity (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  g_assert (test->pix != NULL);

  g_assert_cmpint (ncm_sphere_map_pix_get_middle_size (test->pix) + 
                   2.0 * ncm_sphere_map_pix_get_cap_size (test->pix),
                   ==,
                   ncm_sphere_map_pix_get_npix (test->pix)
                   );  
  
  g_assert_cmpint (ncm_sphere_map_pix_get_nrings_middle (test->pix) + 
                   2.0 * ncm_sphere_map_pix_get_nrings_cap (test->pix),
                   ==,
                   ncm_sphere_map_pix_get_nrings (test->pix)
                   );

  {
    gint64 r_i;
    gint64 j = 0;
    
    for (r_i = 0; r_i < ncm_sphere_map_pix_get_nrings (test->pix); r_i++)
    {
      gint64 i = ncm_sphere_map_pix_get_ring_first_index (test->pix, r_i);

      g_assert_cmpint (i, ==, j);
      
      j += ncm_sphere_map_pix_get_ring_size (test->pix, r_i);
    }

    g_assert_cmpint (j, ==, ncm_sphere_map_pix_get_npix (test->pix));
  }
  
  ncm_sphere_map_pix_set_nside (test->pix, (1 << g_test_rand_int_range (1, 8)));

  if (!test->try)
  {
    test->try = TRUE;
    test_ncm_sphere_map_pix_sanity (test, pdata);
  }
}

void
test_ncm_sphere_map_pix_angles (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  const gint64 n_i = g_test_rand_int_range (200, 1000);
  const gint64 n_j = g_test_rand_int_range (200, 1000);
  gint64 i;

  g_assert (test->pix != NULL);

  for (i = 0; i < n_i; i++)
  {
    const gdouble theta_i = 2.0 * ncm_c_pi () / (n_i - 1.0) * i;
    gint64 j;
    for (j = 0; j < n_j; j++)
    {
      const gdouble phi_j = ncm_c_pi () / (n_j - 1.0) * j;

      gdouble theta_nest, phi_nest;
      gdouble theta_ring, phi_ring;
      
      gint64 nest_index, ring_index;
      
      ncm_sphere_map_pix_ang2pix_nest (test->pix, theta_i, phi_j, &nest_index);
      ncm_sphere_map_pix_ang2pix_ring (test->pix, theta_i, phi_j, &ring_index);

      ncm_sphere_map_pix_pix2ang_nest (test->pix, nest_index, &theta_nest, &phi_nest);
      ncm_sphere_map_pix_pix2ang_ring (test->pix, ring_index, &theta_ring, &phi_ring);

      g_assert_cmpfloat (theta_ring, ==, theta_nest);
      g_assert_cmpfloat (phi_ring, ==, phi_nest);

      g_assert_cmpint (ncm_sphere_map_pix_nest2ring (test->pix, nest_index), ==, ring_index);
      g_assert_cmpint (ncm_sphere_map_pix_ring2nest (test->pix, ring_index), ==, nest_index);
    }
  }
}

void
test_ncm_sphere_map_pix_ring (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  gint64 r_i;
  
  g_assert (test->pix != NULL);

  for (r_i = 0; r_i < ncm_sphere_map_pix_get_nrings (test->pix); r_i++)
  {
    const gint64 ring_fi  = ncm_sphere_map_pix_get_ring_first_index (test->pix, r_i);
    const gint64 r_i_size = ncm_sphere_map_pix_get_ring_size (test->pix, r_i);
    gint64 ir_i;
    gdouble last_phi = 0.0;

    for (ir_i = 0; ir_i < r_i_size; ir_i++)
    {
      const gint64 ring_index = ring_fi + ir_i;
      gdouble theta_i, phi_i;

      ncm_sphere_map_pix_pix2ang_ring (test->pix, ring_index, &theta_i, &phi_i);

      if (ir_i == 0)
        g_assert_cmpfloat (phi_i, >=, last_phi);
      else
      {
        g_assert_cmpfloat (phi_i, >, last_phi);
        ncm_assert_cmpdouble_e (phi_i - last_phi, ==, 2.0 * M_PI / r_i_size, 1.0e-10);
      }

      last_phi = phi_i;
      /*printf ("r_i %ld ir_i %ld ring_size %ld | theta % 20.15g phi % 20.15g\n", r_i, ir_i, r_i_size, theta_i, phi_i);*/
    }
  }
}

void
test_ncm_sphere_map_pix_pix2alm (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_new (NULL);
  const guint lmax = 1024;
  
  g_assert (test->pix != NULL);

  ncm_sphere_map_pix_add_noise (test->pix, 1.0, rng);

  ncm_sphere_map_pix_set_lmax (test->pix, lmax);

  ncm_sphere_map_pix_prepare_alm (test->pix);

  {
    guint l;
    for (l = 0; l <= lmax; l++)
    {
      const gdouble C_l = ncm_sphere_map_pix_get_Cl (test->pix, l);
      g_assert_cmpfloat (test->pix->npix * C_l / (2.0 * (l + 1.0) * (2.0 * M_PI)), >, 0.6);
      g_assert_cmpfloat (test->pix->npix * C_l / (2.0 * (l + 1.0) * (2.0 * M_PI)), <, 1.4);
    }
  }

  ncm_rng_free (rng);
}

void
test_ncm_sphere_map_pix_pix2alm2pix (TestNcmSphereMapPix *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_new (NULL);
  const guint lmax = 1024;
  
  g_assert (test->pix != NULL);

  ncm_sphere_map_pix_add_noise (test->pix, 1.0, rng);
 
  ncm_sphere_map_pix_set_lmax (test->pix, lmax);

  ncm_sphere_map_pix_prepare_alm (test->pix);

  ncm_sphere_map_pix_alm2map (test->pix);

  ncm_rng_free (rng);
}


void
test_ncm_sphere_map_pix_traps (TestNcmSphereMapPix *test, gconstpointer pdata)
{
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_trap_subprocess ("/ncm/sphere_map_pix/invalid/nside/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_sphere_map_pix_invalid_nside (TestNcmSphereMapPix *test, gconstpointer pdata)
{  
  ncm_sphere_map_pix_set_nside (test->pix, (1 << g_test_rand_int_range (1, 8)) + 1);
}

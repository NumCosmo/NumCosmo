/***************************************************************************
 *            test_ncm_sphere_map.c
 *
 *  Sun July 17 17:02:50 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <vitenti@uel.br>
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
#include <glib/gstdio.h>
#include <gsl/gsl_sf_legendre.h>

typedef struct _TestNcmSphereMap
{
  NcmSphereMap *pix;
  guint nside;
  gboolean try;
} TestNcmSphereMap;

void test_ncm_sphere_map_new (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_free (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_sanity (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_angles (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_ring (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_pix2alm (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_pix2alm2pix (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_properties (TestNcmSphereMap *test, gconstpointer pdata);

void test_ncm_sphere_map_traps (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_nside (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_order_roundtrip (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_cap_centres (void);
void test_ncm_sphere_map_fits_roundtrip (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_pixel (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_negative_pixel (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_ring (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_lmax_zero (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_alm_index (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_cross (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_Ctheta (void);
void test_ncm_sphere_map_healpy_pixels (void);
void test_ncm_sphere_map_healpy_ang2pix (void);
void test_ncm_sphere_map_healpy_transforms (void);
void test_ncm_sphere_map_healpy_fits (void);
void test_ncm_sphere_map_healpy_cross (void);
void test_ncm_sphere_map_vectors (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_update_Cl (void);
void test_ncm_sphere_map_fits_options (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_fits_headers (void);
void test_ncm_sphere_map_invalid_fits_explicit (void);
void test_ncm_sphere_map_fits_noorder_subprocess (void);
void test_ncm_sphere_map_invalid_fits_car (void);
void test_ncm_sphere_map_invalid_fits_short (void);
void test_ncm_sphere_map_invalid_fits_exists (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_noise (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_Cls_short (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_Cls_stale (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_get_pix (TestNcmSphereMap *test, gconstpointer pdata);
void test_ncm_sphere_map_invalid_set_map (TestNcmSphereMap *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/sphere_map/sanity", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_sanity,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/angles", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_angles,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/ring", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_ring,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/pix2alm", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_pix2alm,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/pix2alm2pix", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_pix2alm2pix,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/properties", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_properties,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/traps", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_traps,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/invalid/nside/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_nside,
              &test_ncm_sphere_map_free);

  g_test_add ("/ncm/sphere_map/order_roundtrip", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_order_roundtrip,
              &test_ncm_sphere_map_free);

  g_test_add_func ("/ncm/sphere_map/cap_centres", &test_ncm_sphere_map_cap_centres);

#ifdef HAVE_CFITSIO
  g_test_add ("/ncm/sphere_map/fits_roundtrip", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_fits_roundtrip,
              &test_ncm_sphere_map_free);
#endif /* HAVE_CFITSIO */

  g_test_add ("/ncm/sphere_map/invalid/pixel/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_pixel,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/negative_pixel/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_negative_pixel,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/ring/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_ring,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/lmax_zero/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_lmax_zero,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/alm_index/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_alm_index,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/cross/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_cross,
              &test_ncm_sphere_map_free);
  g_test_add_func ("/ncm/sphere_map/Ctheta", &test_ncm_sphere_map_Ctheta);
  g_test_add_func ("/ncm/sphere_map/healpy/pixels", &test_ncm_sphere_map_healpy_pixels);
  g_test_add_func ("/ncm/sphere_map/healpy/ang2pix", &test_ncm_sphere_map_healpy_ang2pix);
  g_test_add_func ("/ncm/sphere_map/healpy/transforms", &test_ncm_sphere_map_healpy_transforms);
  g_test_add_func ("/ncm/sphere_map/healpy/cross", &test_ncm_sphere_map_healpy_cross);
  g_test_add ("/ncm/sphere_map/vectors", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_vectors,
              &test_ncm_sphere_map_free);
  g_test_add_func ("/ncm/sphere_map/update_Cl", &test_ncm_sphere_map_update_Cl);
#ifdef HAVE_CFITSIO
  g_test_add_func ("/ncm/sphere_map/healpy/fits", &test_ncm_sphere_map_healpy_fits);
  g_test_add ("/ncm/sphere_map/fits_options", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_fits_options,
              &test_ncm_sphere_map_free);
  g_test_add_func ("/ncm/sphere_map/fits_headers", &test_ncm_sphere_map_fits_headers);
  g_test_add_func ("/ncm/sphere_map/fits_noorder/subprocess", &test_ncm_sphere_map_fits_noorder_subprocess);
  g_test_add_func ("/ncm/sphere_map/invalid/fits_explicit/subprocess", &test_ncm_sphere_map_invalid_fits_explicit);
  g_test_add_func ("/ncm/sphere_map/invalid/fits_car/subprocess", &test_ncm_sphere_map_invalid_fits_car);
  g_test_add_func ("/ncm/sphere_map/invalid/fits_short/subprocess", &test_ncm_sphere_map_invalid_fits_short);
  g_test_add ("/ncm/sphere_map/invalid/fits_exists/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_fits_exists,
              &test_ncm_sphere_map_free);
#endif /* HAVE_CFITSIO */
  g_test_add ("/ncm/sphere_map/noise", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_noise,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/Cls_short/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_Cls_short,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/Cls_stale/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_Cls_stale,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/get_pix/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_get_pix,
              &test_ncm_sphere_map_free);
  g_test_add ("/ncm/sphere_map/invalid/set_map/subprocess", TestNcmSphereMap, NULL,
              &test_ncm_sphere_map_new,
              &test_ncm_sphere_map_invalid_set_map,
              &test_ncm_sphere_map_free);

  g_test_run ();
}

void
test_ncm_sphere_map_new (TestNcmSphereMap *test, gconstpointer pdata)
{
  test->nside = 64; /*1 << g_test_rand_int_range (1, 8); */

  test->pix = ncm_sphere_map_new (test->nside);
  test->try = FALSE;

  g_assert_true (test->pix != NULL);
  g_assert_true (NCM_IS_SPHERE_MAP (test->pix));

  g_assert_cmpint (ncm_sphere_map_get_nside (test->pix), ==, test->nside);
}

void
test_ncm_sphere_map_free (TestNcmSphereMap *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_sphere_map_free, test->pix);
}

void
test_ncm_sphere_map_sanity (TestNcmSphereMap *test, gconstpointer pdata)
{
  g_assert_true (test->pix != NULL);

  g_assert_cmpint (ncm_sphere_map_get_middle_size (test->pix) +
                   2.0 * ncm_sphere_map_get_cap_size (test->pix),
                   ==,
                   ncm_sphere_map_get_npix (test->pix)
  );

  g_assert_cmpint (ncm_sphere_map_get_nrings_middle (test->pix) +
                   2.0 * ncm_sphere_map_get_nrings_cap (test->pix),
                   ==,
                   ncm_sphere_map_get_nrings (test->pix)
  );

  {
    gint64 r_i;
    gint64 j = 0;

    for (r_i = 0; r_i < ncm_sphere_map_get_nrings (test->pix); r_i++)
    {
      gint64 i = ncm_sphere_map_get_ring_first_index (test->pix, r_i);

      g_assert_cmpint (i, ==, j);

      j += ncm_sphere_map_get_ring_size (test->pix, r_i);
    }

    g_assert_cmpint (j, ==, ncm_sphere_map_get_npix (test->pix));
  }

  ncm_sphere_map_set_nside (test->pix, (1 << g_test_rand_int_range (1, 8)));

  if (!test->try)
  {
    test->try = TRUE;
    test_ncm_sphere_map_sanity (test, pdata);
  }
}

void
test_ncm_sphere_map_angles (TestNcmSphereMap *test, gconstpointer pdata)
{
  const gint64 n_i = g_test_rand_int_range (200, 1000);
  const gint64 n_j = g_test_rand_int_range (200, 1000);
  gint64 i;

  g_assert_true (test->pix != NULL);

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

      ncm_sphere_map_ang2pix_nest (test->pix, theta_i, phi_j, &nest_index);
      ncm_sphere_map_ang2pix_ring (test->pix, theta_i, phi_j, &ring_index);

      ncm_sphere_map_pix2ang_nest (test->pix, nest_index, &theta_nest, &phi_nest);
      ncm_sphere_map_pix2ang_ring (test->pix, ring_index, &theta_ring, &phi_ring);

      g_assert_cmpfloat (theta_ring, ==, theta_nest);
      g_assert_cmpfloat (phi_ring, ==, phi_nest);

      g_assert_cmpint (ncm_sphere_map_nest2ring (test->pix, nest_index), ==, ring_index);
      g_assert_cmpint (ncm_sphere_map_ring2nest (test->pix, ring_index), ==, nest_index);
    }
  }
}

void
test_ncm_sphere_map_ring (TestNcmSphereMap *test, gconstpointer pdata)
{
  gint64 r_i;

  g_assert_true (test->pix != NULL);

  for (r_i = 0; r_i < ncm_sphere_map_get_nrings (test->pix); r_i++)
  {
    const gint64 ring_fi  = ncm_sphere_map_get_ring_first_index (test->pix, r_i);
    const gint64 r_i_size = ncm_sphere_map_get_ring_size (test->pix, r_i);
    gint64 ir_i;
    gdouble last_phi = 0.0;

    for (ir_i = 0; ir_i < r_i_size; ir_i++)
    {
      const gint64 ring_index = ring_fi + ir_i;
      gdouble theta_i, phi_i;

      ncm_sphere_map_pix2ang_ring (test->pix, ring_index, &theta_i, &phi_i);

      if (ir_i == 0)
      {
        g_assert_cmpfloat (phi_i, >=, last_phi);
      }
      else
      {
        g_assert_cmpfloat (phi_i, >, last_phi);
        ncm_assert_cmpdouble_e (phi_i - last_phi, ==, 2.0 * M_PI / r_i_size, 1.0e-10, 0.0);
      }

      last_phi = phi_i;
      /*printf ("r_i %ld ir_i %ld ring_size %ld | theta % 20.15g phi % 20.15g\n", r_i, ir_i, r_i_size, theta_i, phi_i);*/
    }
  }
}

void
test_ncm_sphere_map_pix2alm (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint lmax = 1024;

  g_assert_true (test->pix != NULL);

  ncm_sphere_map_add_noise (test->pix, 1.0, rng);

  ncm_sphere_map_set_lmax (test->pix, lmax);

  ncm_sphere_map_prepare_alm (test->pix);

  {
    gdouble t = 0.0;
    guint l;

    for (l = 0; l <= lmax; l++)
    {
      const gdouble C_l  = ncm_sphere_map_get_Cl (test->pix, l);
      const gdouble NC_l = ncm_sphere_map_get_npix (test->pix) * C_l / (4.0 * M_PI);

      t += NC_l;

      if (l > 10)
      {
        g_assert_cmpfloat (NC_l, >, 1.0e-3);
        g_assert_cmpfloat (NC_l, <, 5.0);
      }
    }

    g_assert_cmpfloat (fabs (t / (lmax + 1.0) - 1.0), <, 1.0e-1);
    /*printf ("%u % 22.15e\n", lmax+1, fabs (t / (lmax + 1.0) - 1.0));*/
  }

  ncm_rng_free (rng);
}

void
test_ncm_sphere_map_pix2alm2pix (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  const guint lmax = 1024;

  g_assert_true (test->pix != NULL);

  ncm_sphere_map_add_noise (test->pix, 1.0, rng);

  ncm_sphere_map_set_lmax (test->pix, lmax);

  ncm_sphere_map_prepare_alm (test->pix);

  ncm_sphere_map_alm2map (test->pix);

  ncm_rng_free (rng);
}

void
test_ncm_sphere_map_properties (TestNcmSphereMap *test, gconstpointer pdata)
{
  const guint test_lmax = 512;
  const guint test_iter = 5;
  guint lmax_get, iter_get;
  NcmSphereMapOrder order_get;
  NcmSphereMapCoordSys coordsys_get;

  g_assert_true (test->pix != NULL);

  /* Test lmax property getter/setter and function */
  ncm_sphere_map_set_lmax (test->pix, test_lmax);
  lmax_get = ncm_sphere_map_get_lmax (test->pix);
  g_assert_cmpuint (lmax_get, ==, test_lmax);

  /* Test lmax property through GObject API */
  g_object_get (test->pix, "lmax", &lmax_get, NULL);
  g_assert_cmpuint (lmax_get, ==, test_lmax);

  /* No refinement unless asked */
  g_assert_cmpuint (ncm_sphere_map_get_iter (test->pix), ==, 0);

  /* Test iter property getter/setter */
  ncm_sphere_map_set_iter (test->pix, test_iter);
  iter_get = ncm_sphere_map_get_iter (test->pix);
  g_assert_cmpuint (iter_get, ==, test_iter);

  /* Test iter property through GObject API */
  g_object_get (test->pix, "iter", &iter_get, NULL);
  g_assert_cmpuint (iter_get, ==, test_iter);

  /* Test order property getter/setter */
  ncm_sphere_map_set_order (test->pix, NCM_SPHERE_MAP_ORDER_RING);
  order_get = ncm_sphere_map_get_order (test->pix);
  g_assert_cmpint (order_get, ==, NCM_SPHERE_MAP_ORDER_RING);

  /* Test order property through GObject API */
  g_object_get (test->pix, "order", &order_get, NULL);
  g_assert_cmpint (order_get, ==, NCM_SPHERE_MAP_ORDER_RING);

  ncm_sphere_map_set_order (test->pix, NCM_SPHERE_MAP_ORDER_NEST);
  order_get = ncm_sphere_map_get_order (test->pix);
  g_assert_cmpint (order_get, ==, NCM_SPHERE_MAP_ORDER_NEST);

  /* Test coordsys property getter/setter */
  ncm_sphere_map_set_coordsys (test->pix, NCM_SPHERE_MAP_COORD_SYS_GALACTIC);
  coordsys_get = ncm_sphere_map_get_coordsys (test->pix);
  g_assert_cmpint (coordsys_get, ==, NCM_SPHERE_MAP_COORD_SYS_GALACTIC);

  /* Test coordsys property through GObject API */
  g_object_get (test->pix, "coordsys", &coordsys_get, NULL);
  g_assert_cmpint (coordsys_get, ==, NCM_SPHERE_MAP_COORD_SYS_GALACTIC);

  ncm_sphere_map_set_coordsys (test->pix, NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC);
  coordsys_get = ncm_sphere_map_get_coordsys (test->pix);
  g_assert_cmpint (coordsys_get, ==, NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC);

  ncm_sphere_map_set_coordsys (test->pix, NCM_SPHERE_MAP_COORD_SYS_CELESTIAL);
  coordsys_get = ncm_sphere_map_get_coordsys (test->pix);
  g_assert_cmpint (coordsys_get, ==, NCM_SPHERE_MAP_COORD_SYS_CELESTIAL);
}

void
test_ncm_sphere_map_traps (TestNcmSphereMap *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/sphere_map/invalid/nside/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nside must be zero or a power of two*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/pixel/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_nest2ring: pixel index 49152 out of range [0, 49152)*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/negative_pixel/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_ring2nest: pixel index -1 out of range*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/ring/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_get_ring_size: ring index 255 out of range [0, 255)*");

  /* A zero lmax used to warn and leave the previous coefficients in place. */
  g_test_trap_subprocess ("/ncm/sphere_map/invalid/lmax_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_prepare_alm: lmax is zero*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/alm_index/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_get_alm: (l, m) = (11, 0) out of range, lmax = 10*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/cross/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the maps differ in lmax (10, 12)*");

  /* A short vector used to be read past its end. */
  g_test_trap_subprocess ("/ncm/sphere_map/invalid/Cls_short/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the vector has 5 values, lmax = 10 needs 11*");

  /* A new lmax zeroes the C_l; C(theta) used to come out zero without complaint. */
  g_test_trap_subprocess ("/ncm/sphere_map/invalid/Cls_stale/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_calc_Ctheta: no C_l*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/get_pix/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_map_get_pix: pixel index -1 out of range*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/set_map/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the array has 3 values, the map 49152 pixels*");
}

void
test_ncm_sphere_map_invalid_nside (TestNcmSphereMap *test, gconstpointer pdata)
{
  ncm_sphere_map_set_nside (test->pix, (1 << g_test_rand_int_range (1, 8)) + 1);
}

/* A change of ordering and back gives the same map bit by bit: the copy went through a
 * gfloat and rounded every pixel at 6e-8. */
void
test_ncm_sphere_map_order_roundtrip (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 11);
  const gint64 npix = ncm_sphere_map_get_npix (test->pix);
  GArray *map       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), npix);
  gint64 i;

  for (i = 0; i < npix; i++)
  {
    const gdouble v = 1.0 + 1.0e-3 * ncm_rng_uniform01_gen (rng);

    g_array_append_val (map, v);
  }

  ncm_sphere_map_set_map (test->pix, map);
  ncm_sphere_map_set_order (test->pix, NCM_SPHERE_MAP_ORDER_NEST);

  for (i = 0; i < npix; i++)
    g_assert_cmpfloat (ncm_sphere_map_get_pix (test->pix, ncm_sphere_map_ring2nest (test->pix, i)), ==, g_array_index (map, gdouble, i));

  ncm_sphere_map_set_order (test->pix, NCM_SPHERE_MAP_ORDER_RING);

  for (i = 0; i < npix; i++)
    g_assert_cmpfloat (ncm_sphere_map_get_pix (test->pix, i), ==, g_array_index (map, gdouble, i));

  g_array_unref (map);
  ncm_rng_free (rng);
}

void
test_ncm_sphere_map_invalid_pixel (TestNcmSphereMap *test, gconstpointer pdata)
{
  ncm_sphere_map_nest2ring (test->pix, ncm_sphere_map_get_npix (test->pix));
}

void
test_ncm_sphere_map_invalid_negative_pixel (TestNcmSphereMap *test, gconstpointer pdata)
{
  ncm_sphere_map_ring2nest (test->pix, -1);
}

void
test_ncm_sphere_map_invalid_ring (TestNcmSphereMap *test, gconstpointer pdata)
{
  ncm_sphere_map_get_ring_size (test->pix, ncm_sphere_map_get_nrings (test->pix));
}

/*
 * The first pixel of each polar-cap ring t (1 to nside - 1) sits at
 * theta = 2 asin (t / (sqrt(6) nside)) in the north and pi minus that in the south. The
 * centres came from acos (1 - t^2 / (3 nside^2)) and sqrt (1 - z^2), which lose precision
 * as nside^2 (3.6e-12 at nside 256, 5.8e-11 at 1024); both are now at rounding level.
 */
void
test_ncm_sphere_map_cap_centres (void)
{
  const gint64 nside = 256;
  NcmSphereMap *smap = ncm_sphere_map_new (nside);
  NcmTriVec *vec     = ncm_trivec_new ();
  gint64 t;

  for (t = 1; t < nside; t++)
  {
    const gdouble theta_n = 2.0 * asin (t / (sqrt (6.0) * nside));
    const gint64 north    = ncm_sphere_map_get_ring_first_index (smap, t - 1);
    const gint64 south    = ncm_sphere_map_get_ring_first_index (smap, ncm_sphere_map_get_nrings (smap) - t);
    gdouble theta, phi;

    ncm_sphere_map_pix2ang_ring (smap, north, &theta, &phi);
    ncm_assert_cmpdouble_e (theta, ==, theta_n, 1.0e-15, 0.0);

    ncm_sphere_map_pix2ang_ring (smap, south, &theta, &phi);
    ncm_assert_cmpdouble_e (theta, ==, M_PI - theta_n, 1.0e-15, 0.0);

    ncm_sphere_map_pix2vec_ring (smap, north, vec);
    ncm_assert_cmpdouble_e (hypot (vec->c[0], vec->c[1]), ==, sin (theta_n), 1.0e-15, 0.0);
  }

  ncm_trivec_free (vec);
  ncm_sphere_map_free (smap);
}

/* Save and load return the map bit by bit, with its ordering and coordinate system: the
 * column was single precision and changed every pixel by up to 6e-8. */
void
test_ncm_sphere_map_fits_roundtrip (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 5);
  const gint64 npix = ncm_sphere_map_get_npix (test->pix);
  GArray *map       = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), npix);
  gchar *dir        = g_dir_make_tmp ("ncm_sphere_map_XXXXXX", NULL);
  gchar *file       = g_build_filename (dir, "map.fits", NULL);
  NcmSphereMap *back;
  gint64 i;

  for (i = 0; i < npix; i++)
  {
    const gdouble v = ncm_rng_gaussian_gen (rng, 0.0, 1.0);

    g_array_append_val (map, v);
  }

  ncm_sphere_map_set_map (test->pix, map);
  ncm_sphere_map_set_order (test->pix, NCM_SPHERE_MAP_ORDER_NEST);
  ncm_sphere_map_set_coordsys (test->pix, NCM_SPHERE_MAP_COORD_SYS_GALACTIC);
  ncm_sphere_map_save_fits (test->pix, file, NULL, TRUE);

  back = ncm_sphere_map_new (1);
  ncm_sphere_map_load_fits (back, file, NULL);

  g_assert_cmpint (ncm_sphere_map_get_nside (back), ==, ncm_sphere_map_get_nside (test->pix));
  g_assert_cmpint (ncm_sphere_map_get_order (back), ==, NCM_SPHERE_MAP_ORDER_NEST);
  g_assert_cmpint (ncm_sphere_map_get_coordsys (back), ==, NCM_SPHERE_MAP_COORD_SYS_GALACTIC);

  for (i = 0; i < npix; i++)
    g_assert_cmpfloat (ncm_sphere_map_get_pix (back, i), ==, ncm_sphere_map_get_pix (test->pix, i));

  ncm_sphere_map_free (back);
  g_unlink (file);
  g_rmdir (dir);
  g_free (file);
  g_free (dir);
  g_array_unref (map);
  ncm_rng_free (rng);
}

void
test_ncm_sphere_map_invalid_lmax_zero (TestNcmSphereMap *test, gconstpointer pdata)
{
  ncm_sphere_map_prepare_alm (test->pix);
}

void
test_ncm_sphere_map_invalid_alm_index (TestNcmSphereMap *test, gconstpointer pdata)
{
  gdouble re, im;

  ncm_sphere_map_set_lmax (test->pix, 10);
  ncm_sphere_map_get_alm (test->pix, 11, 0, &re, &im);
}

void
test_ncm_sphere_map_invalid_cross (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmSphereMap *other = ncm_sphere_map_new (ncm_sphere_map_get_nside (test->pix));

  ncm_sphere_map_set_lmax (test->pix, 10);
  ncm_sphere_map_set_lmax (other, 12);
  ncm_sphere_map_compute_cross_Cl (test->pix, other);
}

/*
 * C(theta) = sum (2l + 1) / (4 pi) C_l P_l(cos theta), with P_l from GSL, for
 * C_l = 1 / (l + 1)^2 up to lmax 200; the spline at reltol 1e-8 was measured at 8.9e-10
 * of the peak.
 */
void
test_ncm_sphere_map_Ctheta (void)
{
  const guint lmax   = 200;
  NcmSphereMap *smap = ncm_sphere_map_new (16);
  NcmVector *Cls     = ncm_vector_new (lmax + 1);
  NcmSpline *Ctheta;
  gdouble peak = 0.0, maxdiff = 0.0;
  guint l, i;

  for (l = 0; l <= lmax; l++)
    ncm_vector_set (Cls, l, 1.0 / gsl_pow_2 (l + 1.0));

  ncm_sphere_map_set_lmax (smap, lmax);
  ncm_sphere_map_set_Cls (smap, Cls);
  Ctheta = ncm_sphere_map_calc_Ctheta (smap, 1.0e-8);

  for (i = 0; i <= 400; i++)
  {
    const gdouble theta = M_PI * i / 400.0;
    gdouble truth       = 0.0;

    for (l = 0; l <= lmax; l++)
      truth += (2.0 * l + 1.0) / (4.0 * M_PI) * ncm_vector_get (Cls, l) * gsl_sf_legendre_Pl (l, cos (theta));

    peak    = GSL_MAX (peak, fabs (truth));
    maxdiff = GSL_MAX (maxdiff, fabs (ncm_spline_eval (Ctheta, theta) - truth));
  }

  g_assert_cmpfloat (maxdiff, <, 1.0e-9 * peak);

  ncm_spline_free (Ctheta);
  ncm_vector_free (Cls);
  ncm_sphere_map_free (smap);
}

/* The added noise has zero mean and variance sd^2: over 49152 pixels, both within four
 * standard errors, for a fixed seed. */
void
test_ncm_sphere_map_noise (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 3);
  const gdouble sd  = 2.5;
  const gint64 npix = ncm_sphere_map_get_npix (test->pix);
  gdouble sum       = 0.0, sum2 = 0.0;
  gint64 i;

  ncm_sphere_map_clear_pixels (test->pix);
  ncm_sphere_map_add_noise (test->pix, sd, rng);

  for (i = 0; i < npix; i++)
  {
    const gdouble v = ncm_sphere_map_get_pix (test->pix, i);

    sum  += v;
    sum2 += v * v;
  }

  g_assert_cmpfloat (fabs (sum / npix), <, 4.0 * sd / sqrt (npix));
  g_assert_cmpfloat (fabs (sum2 / npix - sd * sd), <, 4.0 * sqrt (2.0) * sd * sd / sqrt (npix));

  ncm_rng_free (rng);
}

void
test_ncm_sphere_map_invalid_Cls_short (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmVector *Cls = ncm_vector_new (5);

  ncm_vector_set_zero (Cls);
  ncm_sphere_map_set_lmax (test->pix, 10);
  ncm_sphere_map_set_Cls (test->pix, Cls);
}

void
test_ncm_sphere_map_invalid_Cls_stale (TestNcmSphereMap *test, gconstpointer pdata)
{
  NcmVector *Cls = ncm_vector_new (11);

  ncm_vector_set_all (Cls, 1.0);
  ncm_sphere_map_set_lmax (test->pix, 10);
  ncm_sphere_map_set_Cls (test->pix, Cls);
  ncm_sphere_map_set_lmax (test->pix, 12);
  ncm_sphere_map_calc_Ctheta (test->pix, 1.0e-6);
}

void
test_ncm_sphere_map_invalid_get_pix (TestNcmSphereMap *test, gconstpointer pdata)
{
  ncm_sphere_map_get_pix (test->pix, -1);
}

void
test_ncm_sphere_map_invalid_set_map (TestNcmSphereMap *test, gconstpointer pdata)
{
  GArray *map = g_array_new (FALSE, TRUE, sizeof (gdouble));

  g_array_set_size (map, 3);
  ncm_sphere_map_set_map (test->pix, map);
}

/* The healpy truth tables, from tests/tools/make_sphere_healpy_truth_table.py. */
static NcmObjDictStr *
_test_healpy_table (const gchar *name)
{
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *file         = g_strdup_printf ("truth_tables/sphere/%s", name);
  gchar *path         = ncm_cfg_get_data_filename (file, TRUE);
  NcmObjDictStr *dict = ncm_serialize_dict_str_from_binfile (ser, path);

  g_free (path);
  g_free (file);
  ncm_serialize_free (ser);

  return dict;
}

static NcmMatrix *
_test_healpy_matrix (NcmObjDictStr *dict, const gchar *format, const gint64 n)
{
  gchar *key       = g_strdup_printf (format, n);
  NcmMatrix *table = NCM_MATRIX (ncm_obj_dict_str_peek (dict, key));

  g_assert_nonnull (table);
  g_free (key);

  return table;
}

/*
 * RING <-> NESTED and the pixel centres against healpy: every pixel at nside 1, 2, 4 and
 * 8, and 256 seeded pixels at nside 128 and 256. The indices are exact; the centres
 * agree to rounding.
 */
void
test_ncm_sphere_map_healpy_pixels (void)
{
  NcmObjDictStr *dict    = _test_healpy_table ("healpy_pixels.bin");
  const gint64 nside_a[] = {1, 2, 4, 8, 128, 256};
  guint j;

  for (j = 0; j < G_N_ELEMENTS (nside_a); j++)
  {
    const gint64 nside = nside_a[j];
    const gboolean all = (nside <= 8);
    NcmSphereMap *smap = ncm_sphere_map_new (nside);
    NcmMatrix *table   = _test_healpy_matrix (dict, "nside%" G_GINT64_FORMAT, nside);
    const guint offset = all ? 0 : 1;
    guint r;

    if (all)
      g_assert_cmpint (ncm_matrix_nrows (table), ==, ncm_sphere_map_get_npix (smap));

    for (r = 0; r < ncm_matrix_nrows (table); r++)
    {
      const gint64 i = all ? r : (gint64) ncm_matrix_get (table, r, 0);
      gdouble theta, phi;

      g_assert_cmpint (ncm_sphere_map_ring2nest (smap, i), ==, (gint64) ncm_matrix_get (table, r, offset + 0));
      g_assert_cmpint (ncm_sphere_map_nest2ring (smap, i), ==, (gint64) ncm_matrix_get (table, r, offset + 1));

      ncm_sphere_map_pix2ang_ring (smap, i, &theta, &phi);
      ncm_assert_cmpdouble_e (theta, ==, ncm_matrix_get (table, r, offset + 2), 1.0e-15, 1.0e-15);
      ncm_assert_cmpdouble_e (phi, ==, ncm_matrix_get (table, r, offset + 3), 1.0e-15, 1.0e-15);
    }

    ncm_sphere_map_free (smap);
  }

  ncm_obj_dict_str_unref (dict);
}

/*
 * The pixel containing a direction, RING and NESTED, against healpy at nside 1, 16 and
 * 256: 150 seeded directions and 50 at the poles, the edge of the polar caps, and phi at
 * 0, below 2 pi, negative and beyond 2 pi.
 */
void
test_ncm_sphere_map_healpy_ang2pix (void)
{
  NcmObjDictStr *dict    = _test_healpy_table ("healpy_ang2pix.bin");
  const gint64 nside_a[] = {1, 16, 256};
  guint j;

  for (j = 0; j < G_N_ELEMENTS (nside_a); j++)
  {
    NcmSphereMap *smap = ncm_sphere_map_new (nside_a[j]);
    NcmMatrix *table   = _test_healpy_matrix (dict, "nside%" G_GINT64_FORMAT, nside_a[j]);
    guint r;

    for (r = 0; r < ncm_matrix_nrows (table); r++)
    {
      const gdouble theta = ncm_matrix_get (table, r, 0);
      const gdouble phi   = ncm_matrix_get (table, r, 1);
      gint64 ring, nest;

      ncm_sphere_map_ang2pix_ring (smap, theta, phi, &ring);
      ncm_sphere_map_ang2pix_nest (smap, theta, phi, &nest);

      g_assert_cmpint (ring, ==, (gint64) ncm_matrix_get (table, r, 2));
      g_assert_cmpint (nest, ==, (gint64) ncm_matrix_get (table, r, 3));
    }

    ncm_sphere_map_free (smap);
  }

  ncm_obj_dict_str_unref (dict);
}

/*
 * map2alm (iter 0 and 3), C_l and alm2map against healpy on a seeded nside-8 map, at
 * lmax 23 (3 nside - 1) and 32 (4 nside, where m folds onto the ring Nyquist
 * frequency). healpy's a_lm run m-major in its own index order. Measured at the 1e-14
 * level of the largest coefficient or pixel.
 */
void
test_ncm_sphere_map_healpy_transforms (void)
{
  NcmObjDictStr *dict  = _test_healpy_table ("healpy_transforms.bin");
  NcmVector *map_v     = NCM_VECTOR (ncm_obj_dict_str_peek (dict, "map"));
  const guint lmax_a[] = {23, 32};
  const guint iter_a[] = {0, 3};
  GArray *map          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint j, k, i;

  for (i = 0; i < ncm_vector_len (map_v); i++)
  {
    const gdouble v = ncm_vector_get (map_v, i);

    g_array_append_val (map, v);
  }

  for (j = 0; j < G_N_ELEMENTS (lmax_a); j++)
  {
    const guint lmax = lmax_a[j];
    NcmMatrix *alm   = NULL;

    for (k = 0; k < G_N_ELEMENTS (iter_a); k++)
    {
      NcmSphereMap *smap = ncm_sphere_map_new (8);
      gchar *key         = g_strdup_printf ("alm_lmax%u_iter%u", lmax, iter_a[k]);
      gchar *cl_key      = g_strdup_printf ("cl_lmax%u_iter%u", lmax, iter_a[k]);
      NcmVector *cl      = NCM_VECTOR (ncm_obj_dict_str_peek (dict, cl_key));
      gdouble peak       = 0.0, maxdiff = 0.0;
      guint l, m, row = 0;

      alm = NCM_MATRIX (ncm_obj_dict_str_peek (dict, key));

      ncm_sphere_map_set_lmax (smap, lmax);
      ncm_sphere_map_set_iter (smap, iter_a[k]);
      ncm_sphere_map_set_map (smap, map);
      ncm_sphere_map_prepare_alm (smap);

      for (m = 0; m <= lmax; m++)
      {
        for (l = m; l <= lmax; l++, row++)
        {
          gdouble re, im;

          ncm_sphere_map_get_alm (smap, l, m, &re, &im);
          peak    = GSL_MAX (peak, hypot (ncm_matrix_get (alm, row, 0), ncm_matrix_get (alm, row, 1)));
          maxdiff = GSL_MAX (maxdiff, hypot (re - ncm_matrix_get (alm, row, 0), im - ncm_matrix_get (alm, row, 1)));
        }
      }

      g_assert_cmpuint (row, ==, ncm_matrix_nrows (alm));
      g_assert_cmpfloat (maxdiff, <, 1.0e-13 * peak);

      for (l = 0; l <= lmax; l++)
        ncm_assert_cmpdouble_e (ncm_sphere_map_get_Cl (smap, l), ==, ncm_vector_get (cl, l), 1.0e-12, 0.0);

      g_free (key);
      g_free (cl_key);
      ncm_sphere_map_free (smap);
    }

    /* alm2map of the iter-3 coefficients, the last ones read. */
    {
      NcmSphereMap *smap = ncm_sphere_map_new (8);
      gchar *key         = g_strdup_printf ("synth_lmax%u", lmax);
      NcmVector *synth   = NCM_VECTOR (ncm_obj_dict_str_peek (dict, key));
      gdouble peak       = 0.0, maxdiff = 0.0;
      guint l, m, row = 0;

      ncm_sphere_map_set_lmax (smap, lmax);

      for (m = 0; m <= lmax; m++)
        for (l = m; l <= lmax; l++, row++)
          ncm_sphere_map_set_alm (smap, l, m, ncm_matrix_get (alm, row, 0), ncm_matrix_get (alm, row, 1));

      ncm_sphere_map_alm2map (smap);

      for (i = 0; i < ncm_vector_len (synth); i++)
      {
        peak    = GSL_MAX (peak, fabs (ncm_vector_get (synth, i)));
        maxdiff = GSL_MAX (maxdiff, fabs (ncm_sphere_map_get_pix (smap, i) - ncm_vector_get (synth, i)));
      }

      g_assert_cmpfloat (maxdiff, <, 1.0e-13 * peak);

      g_free (key);
      ncm_sphere_map_free (smap);
    }
  }

  g_array_unref (map);
  ncm_obj_dict_str_unref (dict);
}

/* A NESTED celestial map written by healpy, 1024 pixels per row, loads exactly; its value
 * at NESTED pixel i is (i % 97) / 7 - 3. */
void
test_ncm_sphere_map_healpy_fits (void)
{
  gchar *path        = ncm_cfg_get_data_filename ("truth_tables/sphere/healpy_map_nest_nside16.fits", TRUE);
  NcmSphereMap *smap = ncm_sphere_map_new (1);
  gint64 i;

  ncm_sphere_map_load_fits (smap, path, NULL);

  g_assert_cmpint (ncm_sphere_map_get_nside (smap), ==, 16);
  g_assert_cmpint (ncm_sphere_map_get_order (smap), ==, NCM_SPHERE_MAP_ORDER_NEST);
  g_assert_cmpint (ncm_sphere_map_get_coordsys (smap), ==, NCM_SPHERE_MAP_COORD_SYS_CELESTIAL);

  for (i = 0; i < ncm_sphere_map_get_npix (smap); i++)
    g_assert_cmpfloat (ncm_sphere_map_get_pix (smap, i), ==, (i % 97) / 7.0 - 3.0);

  ncm_sphere_map_free (smap);
  g_free (path);
}

/*
 * Cross spectra against healpy's anafast (map, map2, iter = 3) at lmax 23 and 32; the
 * cross spectrum of a map with itself is its C_l, and swapping the maps changes nothing.
 */
void
test_ncm_sphere_map_healpy_cross (void)
{
  NcmObjDictStr *dict  = _test_healpy_table ("healpy_transforms.bin");
  NcmVector *map1_v    = NCM_VECTOR (ncm_obj_dict_str_peek (dict, "map"));
  NcmVector *map2_v    = NCM_VECTOR (ncm_obj_dict_str_peek (dict, "map2"));
  const guint lmax_a[] = {23, 32};
  GArray *map1         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *map2         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint j, i;

  for (i = 0; i < ncm_vector_len (map1_v); i++)
  {
    const gdouble v1 = ncm_vector_get (map1_v, i);
    const gdouble v2 = ncm_vector_get (map2_v, i);

    g_array_append_val (map1, v1);
    g_array_append_val (map2, v2);
  }

  for (j = 0; j < G_N_ELEMENTS (lmax_a); j++)
  {
    const guint lmax = lmax_a[j];
    NcmSphereMap *s1 = ncm_sphere_map_new (8);
    NcmSphereMap *s2 = ncm_sphere_map_new (8);
    gchar *key       = g_strdup_printf ("cross_lmax%u", lmax);
    NcmVector *truth = NCM_VECTOR (ncm_obj_dict_str_peek (dict, key));
    NcmVector *c12, *c21, *c11;
    gdouble peak = 0.0;
    guint l;

    ncm_sphere_map_set_lmax (s1, lmax);
    ncm_sphere_map_set_lmax (s2, lmax);
    ncm_sphere_map_set_iter (s1, 3);
    ncm_sphere_map_set_iter (s2, 3);
    ncm_sphere_map_set_map (s1, map1);
    ncm_sphere_map_set_map (s2, map2);
    ncm_sphere_map_prepare_alm (s1);
    ncm_sphere_map_prepare_alm (s2);

    c12 = ncm_sphere_map_compute_cross_Cl (s1, s2);
    c21 = ncm_sphere_map_compute_cross_Cl (s2, s1);
    c11 = ncm_sphere_map_compute_cross_Cl (s1, s1);

    for (l = 0; l <= lmax; l++)
      peak = GSL_MAX (peak, fabs (ncm_vector_get (truth, l)));

    for (l = 0; l <= lmax; l++)
    {
      /* Cross spectra change sign, so the scale is the largest one. */
      g_assert_cmpfloat (fabs (ncm_vector_get (c12, l) - ncm_vector_get (truth, l)), <, 1.0e-12 * peak);
      g_assert_cmpfloat (ncm_vector_get (c21, l), ==, ncm_vector_get (c12, l));
      ncm_assert_cmpdouble_e (ncm_vector_get (c11, l), ==, ncm_sphere_map_get_Cl (s1, l), 1.0e-14, 0.0);
    }

    ncm_vector_free (c12);
    ncm_vector_free (c21);
    ncm_vector_free (c11);
    g_free (key);
    ncm_sphere_map_free (s1);
    ncm_sphere_map_free (s2);
  }

  g_array_unref (map1);
  g_array_unref (map2);
  ncm_obj_dict_str_unref (dict);
}

/* pix2vec gives the unit vector of pix2ang's centre, in both orderings, and vec2pix
 * takes it back to the pixel. */
void
test_ncm_sphere_map_vectors (TestNcmSphereMap *test, gconstpointer pdata)
{
  const gint64 npix = ncm_sphere_map_get_npix (test->pix);
  NcmTriVec *vec    = ncm_trivec_new ();
  gint64 i;

  for (i = 0; i < npix; i += 7)
  {
    gdouble theta, phi;
    gint64 back;

    ncm_sphere_map_pix2vec_ring (test->pix, i, vec);
    ncm_sphere_map_pix2ang_ring (test->pix, i, &theta, &phi);
    ncm_assert_cmpdouble_e (ncm_trivec_norm (vec), ==, 1.0, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (vec->c[2], ==, cos (theta), 1.0e-15, 1.0e-15);
    ncm_assert_cmpdouble_e (atan2 (vec->c[1], vec->c[0]) + ((phi > M_PI) ? 2.0 * M_PI : 0.0), ==, phi, 1.0e-14, 1.0e-15);
    ncm_sphere_map_vec2pix_ring (test->pix, vec, &back);
    g_assert_cmpint (back, ==, i);

    ncm_sphere_map_pix2vec_nest (test->pix, i, vec);
    ncm_sphere_map_vec2pix_nest (test->pix, vec, &back);
    g_assert_cmpint (back, ==, i);
  }

  ncm_trivec_free (vec);
}

/* get_Cl returns the C_l of the last computation: set_alm leaves them until update_Cl. */
void
test_ncm_sphere_map_update_Cl (void)
{
  NcmSphereMap *smap = ncm_sphere_map_new (8);
  GArray *map        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble C1_before, re, im;
  gint64 i;

  for (i = 0; i < ncm_sphere_map_get_npix (smap); i++)
  {
    const gdouble v = sin (0.1 * i);

    g_array_append_val (map, v);
  }

  ncm_sphere_map_set_lmax (smap, 16);
  ncm_sphere_map_set_map (smap, map);
  ncm_sphere_map_prepare_alm (smap);
  C1_before = ncm_sphere_map_get_Cl (smap, 1);

  ncm_sphere_map_get_alm (smap, 1, 1, &re, &im);
  ncm_sphere_map_set_alm (smap, 1, 1, 2.0 * re, 2.0 * im);
  g_assert_cmpfloat (ncm_sphere_map_get_Cl (smap, 1), ==, C1_before);

  ncm_sphere_map_update_Cl (smap);
  {
    gdouble re0, im0;

    ncm_sphere_map_get_alm (smap, 1, 0, &re0, &im0);
    ncm_assert_cmpdouble_e (ncm_sphere_map_get_Cl (smap, 1), ==, (re0 * re0 + im0 * im0 + 2.0 * 4.0 * (re * re + im * im)) / 3.0, 1.0e-14, 0.0);
  }

  g_array_unref (map);
  ncm_sphere_map_free (smap);
}

/* A named column, and overwrite: a second save replaces the file only when asked. */
void
test_ncm_sphere_map_fits_options (TestNcmSphereMap *test, gconstpointer pdata)
{
  gchar *dir         = g_dir_make_tmp ("ncm_sphere_map_XXXXXX", NULL);
  gchar *file        = g_build_filename (dir, "map.fits", NULL);
  NcmSphereMap *back = ncm_sphere_map_new (1);
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, 9);
  gint64 i;

  ncm_sphere_map_add_noise (test->pix, 1.0, rng);
  ncm_sphere_map_save_fits (test->pix, file, "TEMPERATURE", TRUE);
  ncm_sphere_map_add_noise (test->pix, 1.0, rng);
  ncm_sphere_map_save_fits (test->pix, file, "TEMPERATURE", TRUE);

  ncm_sphere_map_load_fits (back, file, "TEMPERATURE");

  for (i = 0; i < ncm_sphere_map_get_npix (back); i++)
    g_assert_cmpfloat (ncm_sphere_map_get_pix (back, i), ==, ncm_sphere_map_get_pix (test->pix, i));

  g_unlink (file);
  g_rmdir (dir);
  g_free (file);
  g_free (dir);
  ncm_rng_free (rng);
  ncm_sphere_map_free (back);
}

/* The header cases: a file without ORDERING loads as RING, with a warning (it used to
 * abort on an unset buffer), and the rejections abort with their messages. */
void
test_ncm_sphere_map_fits_headers (void)
{
  g_test_trap_subprocess ("/ncm/sphere_map/fits_noorder/subprocess", 0, 0);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stderr ("*Could not find ORDERING in the fits file, assuming RING*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/fits_explicit/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*partial-sky map (INDXSCHM = EXPLICIT)*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/fits_car/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*has PIXTYPE `CAR', not HEALPIX*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/fits_short/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*holds 100 values (100 rows of 1)*");

  g_test_trap_subprocess ("/ncm/sphere_map/invalid/fits_exists/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_sphere_map_fits_noorder_subprocess (void)
{
  gchar *path        = ncm_cfg_get_data_filename ("truth_tables/sphere/fits_noorder.fits", TRUE);
  NcmSphereMap *smap = ncm_sphere_map_new (1);
  gint64 i;

  /* The missing ORDERING warns; only the load after it is under test. */
  g_log_set_always_fatal (G_LOG_FATAL_MASK);
  ncm_sphere_map_load_fits (smap, path, NULL);

  g_assert_cmpint (ncm_sphere_map_get_order (smap), ==, NCM_SPHERE_MAP_ORDER_RING);

  for (i = 0; i < ncm_sphere_map_get_npix (smap); i++)
    g_assert_cmpfloat (ncm_sphere_map_get_pix (smap, i), ==, (gdouble) i);

  ncm_sphere_map_free (smap);
  g_free (path);
}

static void
_test_load_fixture (const gchar *name)
{
  gchar *file        = g_strdup_printf ("truth_tables/sphere/%s", name);
  gchar *path        = ncm_cfg_get_data_filename (file, TRUE);
  NcmSphereMap *smap = ncm_sphere_map_new (1);

  ncm_sphere_map_load_fits (smap, path, NULL);
}

void
test_ncm_sphere_map_invalid_fits_explicit (void)
{
  _test_load_fixture ("fits_explicit.fits");
}

void
test_ncm_sphere_map_invalid_fits_car (void)
{
  _test_load_fixture ("fits_car.fits");
}

void
test_ncm_sphere_map_invalid_fits_short (void)
{
  _test_load_fixture ("fits_short.fits");
}

/* Without overwrite, saving onto an existing file is an error. */
void
test_ncm_sphere_map_invalid_fits_exists (TestNcmSphereMap *test, gconstpointer pdata)
{
  gchar *dir  = g_dir_make_tmp ("ncm_sphere_map_XXXXXX", NULL);
  gchar *file = g_build_filename (dir, "map.fits", NULL);

  ncm_sphere_map_save_fits (test->pix, file, NULL, FALSE);
  ncm_sphere_map_save_fits (test->pix, file, NULL, FALSE);
}


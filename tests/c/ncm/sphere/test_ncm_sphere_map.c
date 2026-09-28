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


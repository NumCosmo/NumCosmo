/***************************************************************************
 *            test_ncm_quaternion.c
 *
 *  Mon Sep 02 17:04:12 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

void test_ncm_trivec_new (void);
void test_ncm_trivec_new_full (void);
void test_ncm_trivec_new_full_c (void);
void test_ncm_trivec_new_sphere (void);
void test_ncm_trivec_new_astro_coord (void);
void test_ncm_trivec_new_astro_ra_dec (void);

void test_ncm_trivec_dup (void);
void test_ncm_trivec_memcpy (void);
void test_ncm_trivec_set_0 (void);
void test_ncm_trivec_scale (void);
void test_ncm_trivec_norm (void);
void test_ncm_trivec_dot (void);
void test_ncm_trivec_normalize (void);

void test_ncm_trivec_get_phi (void);
void test_ncm_trivec_set_spherical_coord (void);
void test_ncm_trivec_get_spherical_coord (void);

void test_ncm_trivec_set_astro_coord (void);
void test_ncm_trivec_get_astro_coord (void);

void test_ncm_trivec_set_astro_ra_dec (void);
void test_ncm_trivec_get_astro_ra_dec (void);

void test_ncm_quaternion_new (void);
void test_ncm_quaternion_new_from_vector (void);
void test_ncm_quaternion_new_from_data (void);

void test_ncm_quaternion_dup (void);
void test_ncm_quaternion_memcpy (void);
void test_ncm_quaternion_set_from_data (void);
void test_ncm_quaternion_set_random (void);
void test_ncm_quaternion_set_I (void);
void test_ncm_quaternion_set_0 (void);

void test_ncm_quaternion_norm (void);
void test_ncm_quaternion_normalize (void);
void test_ncm_quaternion_conjugate (void);

void test_ncm_quaternion_mul (void);
void test_ncm_quaternion_lmul (void);
void test_ncm_quaternion_rmul (void);

void test_ncm_quaternion_conjugate_u_mul (void);
void test_ncm_quaternion_conjugate_q_mul (void);

void test_ncm_quaternion_rotate_norm (void);
void test_ncm_quaternion_rotate_spherical (void);
void test_ncm_quaternion_rotate_spherical_onestep (void);
void test_ncm_quaternion_rotate_axis (void);
void test_ncm_quaternion_rotate_inverse (void);

void test_ncm_quaternion_set_to_rotate_to_x (void);
void test_ncm_quaternion_set_to_rotate_to_z (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/trivec/new", test_ncm_trivec_new);
  g_test_add_func ("/ncm/trivec/new_full", test_ncm_trivec_new_full);
  g_test_add_func ("/ncm/trivec/new_full_c", test_ncm_trivec_new_full_c);
  g_test_add_func ("/ncm/trivec/new_sphere", test_ncm_trivec_new_sphere);
  g_test_add_func ("/ncm/trivec/new_astro_coord", test_ncm_trivec_new_astro_coord);
  g_test_add_func ("/ncm/trivec/new_astro_ra_dec", test_ncm_trivec_new_astro_ra_dec);

  g_test_add_func ("/ncm/trivec/dup", test_ncm_trivec_dup);
  g_test_add_func ("/ncm/trivec/memcpy", test_ncm_trivec_memcpy);
  g_test_add_func ("/ncm/trivec/set_0", test_ncm_trivec_set_0);
  g_test_add_func ("/ncm/trivec/scale", test_ncm_trivec_scale);
  g_test_add_func ("/ncm/trivec/norm", test_ncm_trivec_norm);
  g_test_add_func ("/ncm/trivec/dot", test_ncm_trivec_dot);
  g_test_add_func ("/ncm/trivec/normalize", test_ncm_trivec_normalize);

  g_test_add_func ("/ncm/trivec/get_phi", test_ncm_trivec_get_phi);
  g_test_add_func ("/ncm/trivec/set_spherical_coord", test_ncm_trivec_set_spherical_coord);
  g_test_add_func ("/ncm/trivec/get_spherical_coord", test_ncm_trivec_get_spherical_coord);

  g_test_add_func ("/ncm/trivec/set_astro_coord", test_ncm_trivec_set_astro_coord);
  g_test_add_func ("/ncm/trivec/get_astro_coord", test_ncm_trivec_get_astro_coord);

  g_test_add_func ("/ncm/trivec/set_astro_ra_dec", test_ncm_trivec_set_astro_ra_dec);
  g_test_add_func ("/ncm/trivec/get_astro_ra_dec", test_ncm_trivec_get_astro_ra_dec);

  g_test_add_func ("/ncm/quaternion/new", test_ncm_quaternion_new);
  g_test_add_func ("/ncm/quaternion/new_from_vector", test_ncm_quaternion_new_from_vector);
  g_test_add_func ("/ncm/quaternion/new_from_data", test_ncm_quaternion_new_from_data);

  g_test_add_func ("/ncm/quaternion/dup", test_ncm_quaternion_dup);
  g_test_add_func ("/ncm/quaternion/memcpy", test_ncm_quaternion_memcpy);
  g_test_add_func ("/ncm/quaternion/set_from_data", test_ncm_quaternion_set_from_data);
  g_test_add_func ("/ncm/quaternion/set_random", test_ncm_quaternion_set_random);
  g_test_add_func ("/ncm/quaternion/set_I", test_ncm_quaternion_set_I);
  g_test_add_func ("/ncm/quaternion/set_0", test_ncm_quaternion_set_0);

  g_test_add_func ("/ncm/quaternion/norm", test_ncm_quaternion_norm);
  g_test_add_func ("/ncm/quaternion/normalize", test_ncm_quaternion_normalize);
  g_test_add_func ("/ncm/quaternion/conjugate", test_ncm_quaternion_conjugate);

  g_test_add_func ("/ncm/quaternion/mul", test_ncm_quaternion_mul);
  g_test_add_func ("/ncm/quaternion/lmul", test_ncm_quaternion_lmul);
  g_test_add_func ("/ncm/quaternion/rmul", test_ncm_quaternion_rmul);

  g_test_add_func ("/ncm/quaternion/conjugate_u_mul", test_ncm_quaternion_conjugate_u_mul);
  g_test_add_func ("/ncm/quaternion/conjugate_q_mul", test_ncm_quaternion_conjugate_q_mul);

  g_test_add_func ("/ncm/quaternion/rotate/norm", test_ncm_quaternion_rotate_norm);
  g_test_add_func ("/ncm/quaternion/rotate/spherical", test_ncm_quaternion_rotate_spherical);
  g_test_add_func ("/ncm/quaternion/rotate/spherical_onestep", test_ncm_quaternion_rotate_spherical_onestep);
  g_test_add_func ("/ncm/quaternion/rotate/axis", test_ncm_quaternion_rotate_axis);
  g_test_add_func ("/ncm/quaternion/rotate/inverse", test_ncm_quaternion_rotate_inverse);

  g_test_add_func ("/ncm/quaternion/set_to_rotate_to_x", test_ncm_quaternion_set_to_rotate_to_x);
  g_test_add_func ("/ncm/quaternion/set_to_rotate_to_z", test_ncm_quaternion_set_to_rotate_to_z);

  g_test_run ();
}

#define NTESTS 1000

void
test_ncm_trivec_new (void)
{
  NcmTriVec *v = ncm_trivec_new ();

  g_assert_cmpfloat (v->c[0], ==, 0.0);
  g_assert_cmpfloat (v->c[1], ==, 0.0);
  g_assert_cmpfloat (v->c[2], ==, 0.0);

  ncm_trivec_free (v);
}

void
test_ncm_trivec_new_full (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble c[3] = {
      g_test_rand_double_range (-100.0, 100.0),
      g_test_rand_double_range (-100.0, 100.0),
      g_test_rand_double_range (-100.0, 100.0)
    };
    NcmTriVec *v = ncm_trivec_new_full (c);

    ncm_assert_cmpdouble_e (v->c[0], ==, c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, c[2], reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_new_full_c (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v    = ncm_trivec_new_full_c (x, y, z);

    ncm_assert_cmpdouble_e (v->c[0], ==, x, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, y, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, z, reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_new_sphere (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble r     = g_test_rand_double_range (0.0, 100.0);
    const gdouble theta = g_test_rand_double_range (0.0, M_PI);
    const gdouble phi   = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmTriVec *v        = ncm_trivec_new_sphere (r, theta, phi);

    ncm_assert_cmpdouble_e (v->c[0], ==, r * sin (theta) * cos (phi), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, r * sin (theta) * sin (phi), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, r * cos (theta), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_new_astro_coord (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble r     = g_test_rand_double_range (0.0, 100.0);
    const gdouble delta = g_test_rand_double_range (-M_PI_2, M_PI_2);
    const gdouble alpha = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmTriVec *v        = ncm_trivec_new_astro_coord (r, delta, alpha);

    ncm_assert_cmpdouble_e (v->c[0], ==, r * cos (delta) * cos (alpha), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, r * cos (delta) * sin (alpha), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, r * sin (delta), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_new_astro_ra_dec (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble r   = g_test_rand_double_range (0.0, 100.0);
    const gdouble ra  = g_test_rand_double_range (0.0, 360.0);
    const gdouble dec = g_test_rand_double_range (-90.0, 90.0);
    NcmTriVec *v      = ncm_trivec_new_astro_ra_dec (r, ra, dec);

    ncm_assert_cmpdouble_e (v->c[0], ==, r * cos (M_PI * ra / 180.0) * cos (M_PI * dec / 180.0), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, r * sin (M_PI * ra / 180.0) * cos (M_PI * dec / 180.0), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, r * sin (M_PI * dec / 180.0), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_dup (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v    = ncm_trivec_new_full_c (x, y, z);
    NcmTriVec *u    = ncm_trivec_dup (v);

    ncm_assert_cmpdouble_e (u->c[0], ==, x, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[1], ==, y, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[2], ==, z, reltol, abstol);

    ncm_trivec_free (v);
    ncm_trivec_free (u);
  }
}

void
test_ncm_trivec_memcpy (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v    = ncm_trivec_new_full_c (x, y, z);
    NcmTriVec *u    = ncm_trivec_new ();

    ncm_trivec_memcpy (u, v);

    ncm_assert_cmpdouble_e (u->c[0], ==, x, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[1], ==, y, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[2], ==, z, reltol, abstol);

    ncm_trivec_free (v);
    ncm_trivec_free (u);
  }
}

void
test_ncm_trivec_set_0 (void)
{
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));

    ncm_trivec_set_0 (v);

    g_assert_cmpfloat (v->c[0], ==, 0.0);
    g_assert_cmpfloat (v->c[1], ==, 0.0);
    g_assert_cmpfloat (v->c[2], ==, 0.0);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_scale (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble scale = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v        = ncm_trivec_new_full_c (x, y, z);

    ncm_trivec_scale (v, scale);

    ncm_assert_cmpdouble_e (v->c[0], ==, x * scale, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, y * scale, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, z * scale, reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_norm (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z     = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v        = ncm_trivec_new_full_c (x, y, z);
    const gdouble norma = sqrt (x * x + y * y + z * z);

    ncm_assert_cmpdouble_e (ncm_trivec_norm (v), ==, norma, reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_dot (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble x1 = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y1 = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z1 = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v1    = ncm_trivec_new_full_c (x, y, z);
    NcmTriVec *v2    = ncm_trivec_new_full_c (x1, y1, z1);

    ncm_assert_cmpdouble_e (ncm_trivec_dot (v1, v2), ==, x * x1 + y * y1 + z * z1, reltol, abstol);

    ncm_trivec_free (v1);
    ncm_trivec_free (v2);
  }
}

void
test_ncm_trivec_normalize (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z     = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v        = ncm_trivec_new_full_c (x, y, z);
    const gdouble norma = sqrt (x * x + y * y + z * z);

    ncm_trivec_normalize (v);

    ncm_assert_cmpdouble_e (v->c[0], ==, x / norma, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, y / norma, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, z / norma, reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_get_phi (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    const gdouble phi = atan2 (v->c[1], v->c[0]);

    ncm_assert_cmpdouble_e (ncm_trivec_get_phi (v), ==, phi, reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_set_spherical_coord (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble r     = g_test_rand_double_range (0.0, 100.0);
    const gdouble theta = g_test_rand_double_range (0.0, M_PI);
    const gdouble phi   = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmTriVec *v        = ncm_trivec_new ();

    ncm_trivec_set_spherical_coord (v, r, theta, phi);

    ncm_assert_cmpdouble_e (v->c[0], ==, r * sin (theta) * cos (phi), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, r * sin (theta) * sin (phi), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, r * cos (theta), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_get_spherical_coord (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    gdouble r, theta, phi;

    ncm_trivec_get_spherical_coord (v, &r, &theta, &phi);

    ncm_assert_cmpdouble_e (r, ==, ncm_trivec_norm (v), reltol, abstol);
    ncm_assert_cmpdouble_e (theta, ==, acos (v->c[2] / r), reltol, abstol);
    ncm_assert_cmpdouble_e (phi, ==, ncm_trivec_get_phi (v), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_set_astro_coord (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble r     = g_test_rand_double_range (0.0, 100.0);
    const gdouble delta = g_test_rand_double_range (-M_PI_2, M_PI_2);
    const gdouble alpha = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmTriVec *v        = ncm_trivec_new ();

    ncm_trivec_set_astro_coord (v, r, delta, alpha);

    ncm_assert_cmpdouble_e (v->c[0], ==, r * cos (delta) * cos (alpha), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, r * cos (delta) * sin (alpha), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, r * sin (delta), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_get_astro_coord (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    gdouble r, delta, alpha;

    ncm_trivec_get_astro_coord (v, &r, &delta, &alpha);

    ncm_assert_cmpdouble_e (r, ==, ncm_trivec_norm (v), reltol, abstol);
    ncm_assert_cmpdouble_e (delta, ==, asin (v->c[2] / r), reltol, abstol);
    ncm_assert_cmpdouble_e (alpha, ==, ncm_trivec_get_phi (v), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_set_astro_ra_dec (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble r   = g_test_rand_double_range (0.0, 100.0);
    const gdouble ra  = g_test_rand_double_range (0.0, 360.0);
    const gdouble dec = g_test_rand_double_range (-90.0, 90.0);
    NcmTriVec *v      = ncm_trivec_new ();

    ncm_trivec_set_astro_ra_dec (v, r, ra, dec);

    ncm_assert_cmpdouble_e (v->c[0], ==, r * cos (M_PI * ra / 180.0) * cos (M_PI * dec / 180.0), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, r * sin (M_PI * ra / 180.0) * cos (M_PI * dec / 180.0), reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, r * sin (M_PI * dec / 180.0), reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_trivec_get_astro_ra_dec (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    gdouble r, ra, dec;

    ncm_trivec_get_astro_ra_dec (v, &r, &ra, &dec);

    ncm_assert_cmpdouble_e (r, ==, ncm_trivec_norm (v), reltol, abstol);
    ncm_assert_cmpdouble_e (dec, ==, asin (v->c[2] / r) * 180.0 / M_PI, reltol, abstol);
    ncm_assert_cmpdouble_e (ra, ==, ncm_trivec_get_phi (v) * 180.0 / M_PI, reltol, abstol);

    ncm_trivec_free (v);
  }
}

void
test_ncm_quaternion_new (void)
{
  NcmQuaternion *q = ncm_quaternion_new ();

  g_assert_cmpfloat (q->s, ==, 0.0);
  g_assert_cmpfloat (q->v.c[0], ==, 0.0);
  g_assert_cmpfloat (q->v.c[1], ==, 0.0);
  g_assert_cmpfloat (q->v.c[2], ==, 0.0);

  ncm_quaternion_free (q);
}

void
test_ncm_quaternion_new_from_vector (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z  = g_test_rand_double_range (-100.0, 100.0);
    NcmTriVec *v     = ncm_trivec_new_full_c (x, y, z);
    NcmQuaternion *q = ncm_quaternion_new_from_vector (v);

    ncm_assert_cmpdouble_e (q->s, ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[0], ==, x, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[1], ==, y, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[2], ==, z, reltol, abstol);

    ncm_quaternion_free (q);
    ncm_trivec_free (v);
  }
}

void
test_ncm_quaternion_new_from_data (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble theta = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmQuaternion *q    = ncm_quaternion_new_from_data (x, y, z, theta);
    const gdouble norma = sqrt (x * x + y * y + z * z);

    ncm_assert_cmpdouble_e (q->s, ==, cos (0.5 * theta), reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[0], ==, x * sin (0.5 * theta) / norma, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[1], ==, y * sin (0.5 * theta) / norma, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[2], ==, z * sin (0.5 * theta) / norma, reltol, abstol);

    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_dup (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble s  = g_test_rand_double_range (-100.0, 100.0);
    NcmQuaternion *q = ncm_quaternion_new_from_data (x, y, z, s);
    NcmQuaternion *r = ncm_quaternion_dup (q);

    ncm_assert_cmpdouble_e (r->s, ==, q->s, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, q->v.c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, q->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, q->v.c[2], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (r);
  }
}

void
test_ncm_quaternion_memcpy (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z  = g_test_rand_double_range (-100.0, 100.0);
    const gdouble s  = g_test_rand_double_range (-100.0, 100.0);
    NcmQuaternion *q = ncm_quaternion_new_from_data (x, y, z, s);
    NcmQuaternion *r = ncm_quaternion_new ();

    ncm_quaternion_memcpy (r, q);

    ncm_assert_cmpdouble_e (r->s, ==, q->s, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, q->v.c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, q->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, q->v.c[2], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (r);
  }
}

void
test_ncm_quaternion_set_from_data (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble x     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble y     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble z     = g_test_rand_double_range (-100.0, 100.0);
    const gdouble theta = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmQuaternion *q    = ncm_quaternion_new ();

    ncm_quaternion_set_from_data (q, x, y, z, theta);

    const gdouble norma = sqrt (x * x + y * y + z * z);

    ncm_assert_cmpdouble_e (q->s, ==, cos (0.5 * theta), reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[0], ==, x * sin (0.5 * theta) / norma, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[1], ==, y * sin (0.5 * theta) / norma, reltol, abstol);
    ncm_assert_cmpdouble_e (q->v.c[2], ==, z * sin (0.5 * theta) / norma, reltol, abstol);

    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_set_random (void)
{
  NcmRNG *rng          = ncm_rng_new (NULL);
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new ();

    ncm_quaternion_set_random (q, rng);

    ncm_assert_cmpdouble_e (q->s * q->s +
                            q->v.c[0] * q->v.c[0] +
                            q->v.c[1] * q->v.c[1] +
                            q->v.c[2] * q->v.c[2], ==, 1.0, reltol, abstol);

    ncm_quaternion_free (q);
  }

  ncm_rng_free (rng);
}

void
test_ncm_quaternion_set_I (void)
{
  NcmQuaternion *q     = ncm_quaternion_new ();
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;

  ncm_quaternion_set_I (q);

  ncm_assert_cmpdouble_e (q->s, ==, 1.0, reltol, abstol);
  ncm_assert_cmpdouble_e (q->v.c[0], ==, 0.0, reltol, abstol);
  ncm_assert_cmpdouble_e (q->v.c[1], ==, 0.0, reltol, abstol);
  ncm_assert_cmpdouble_e (q->v.c[2], ==, 0.0, reltol, abstol);

  ncm_quaternion_free (q);
}

void
test_ncm_quaternion_set_0 (void)
{
  NcmQuaternion *q     = ncm_quaternion_new ();
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;

  ncm_quaternion_set_0 (q);

  ncm_assert_cmpdouble_e (q->s, ==, 0.0, reltol, abstol);
  ncm_assert_cmpdouble_e (q->v.c[0], ==, 0.0, reltol, abstol);
  ncm_assert_cmpdouble_e (q->v.c[1], ==, 0.0, reltol, abstol);
  ncm_assert_cmpdouble_e (q->v.c[2], ==, 0.0, reltol, abstol);

  ncm_quaternion_free (q);
}

void
test_ncm_quaternion_norm (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));

    ncm_assert_cmpdouble_e (ncm_quaternion_norm (q), ==,
                            sqrt (q->s * q->s +
                                  q->v.c[0] * q->v.c[0] +
                                  q->v.c[1] * q->v.c[1] +
                                  q->v.c[2] * q->v.c[2]), reltol, abstol);

    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_normalize (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));

    ncm_quaternion_normalize (q);

    ncm_assert_cmpdouble_e (q->s * q->s +
                            q->v.c[0] * q->v.c[0] +
                            q->v.c[1] * q->v.c[1] +
                            q->v.c[2] * q->v.c[2], ==, 1.0, reltol, abstol);

    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_conjugate (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *r = ncm_quaternion_new ();

    ncm_quaternion_memcpy (r, q);
    ncm_quaternion_conjugate (r);

    ncm_assert_cmpdouble_e (r->s, ==, q->s, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, -q->v.c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, -q->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, -q->v.c[2], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (r);
  }
}

void
test_ncm_quaternion_mul (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *r = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *s = ncm_quaternion_new ();

    ncm_quaternion_mul (q, r, s);

    ncm_assert_cmpdouble_e (s->s, ==, q->s * r->s - ncm_trivec_dot (&q->v, &r->v), reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[0], ==, q->s * r->v.c[0] + r->s * q->v.c[0] + q->v.c[1] * r->v.c[2] - q->v.c[2] * r->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[1], ==, q->s * r->v.c[1] + r->s * q->v.c[1] + q->v.c[2] * r->v.c[0] - q->v.c[0] * r->v.c[2], reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[2], ==, q->s * r->v.c[2] + r->s * q->v.c[2] + q->v.c[0] * r->v.c[1] - q->v.c[1] * r->v.c[0], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (r);
    ncm_quaternion_free (s);
  }
}

void
test_ncm_quaternion_lmul (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *r = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *s = ncm_quaternion_new ();

    ncm_quaternion_memcpy (s, r);
    ncm_quaternion_lmul (s, q);

    ncm_assert_cmpdouble_e (s->s, ==, q->s * r->s - ncm_trivec_dot (&r->v, &q->v), reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[0], ==, q->s * r->v.c[0] + r->s * q->v.c[0] + q->v.c[1] * r->v.c[2] - q->v.c[2] * r->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[1], ==, q->s * r->v.c[1] + r->s * q->v.c[1] + q->v.c[2] * r->v.c[0] - q->v.c[0] * r->v.c[2], reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[2], ==, q->s * r->v.c[2] + r->s * q->v.c[2] + q->v.c[0] * r->v.c[1] - q->v.c[1] * r->v.c[0], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (r);
    ncm_quaternion_free (s);
  }
}

void
test_ncm_quaternion_rmul (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *r = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *s = ncm_quaternion_new ();

    ncm_quaternion_memcpy (s, q);
    ncm_quaternion_rmul (s, r);

    ncm_assert_cmpdouble_e (s->s, ==, q->s * r->s - ncm_trivec_dot (&q->v, &r->v), reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[0], ==, q->s * r->v.c[0] + r->s * q->v.c[0] + q->v.c[1] * r->v.c[2] - q->v.c[2] * r->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[1], ==, q->s * r->v.c[1] + r->s * q->v.c[1] + q->v.c[2] * r->v.c[0] - q->v.c[0] * r->v.c[2], reltol, abstol);
    ncm_assert_cmpdouble_e (s->v.c[2], ==, q->s * r->v.c[2] + r->s * q->v.c[2] + q->v.c[0] * r->v.c[1] - q->v.c[1] * r->v.c[0], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (r);
    ncm_quaternion_free (s);
  }
}

void
test_ncm_quaternion_conjugate_u_mul (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *u = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *r = ncm_quaternion_new ();
    NcmQuaternion *s = ncm_quaternion_new ();

    ncm_quaternion_conjugate_u_mul (q, q, r);

    ncm_assert_cmpdouble_e (r->s, ==, 1.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, 0.0, reltol, abstol);

    ncm_quaternion_conjugate_u_mul (u, u, r);

    ncm_assert_cmpdouble_e (r->s, ==, 1.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, 0.0, reltol, abstol);

    ncm_quaternion_conjugate_u_mul (q, u, r);
    ncm_quaternion_conjugate_u_mul (u, q, s);

    ncm_quaternion_conjugate (s);

    ncm_assert_cmpdouble_e (r->s, ==, s->s, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, s->v.c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, s->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, s->v.c[2], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (u);
    ncm_quaternion_free (r);
    ncm_quaternion_free (s);
  }
}

void
test_ncm_quaternion_conjugate_q_mul (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-15;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *u = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmQuaternion *r = ncm_quaternion_new ();
    NcmQuaternion *s = ncm_quaternion_new ();

    ncm_quaternion_conjugate_q_mul (q, q, r);

    ncm_assert_cmpdouble_e (r->s, ==, 1.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, 0.0, reltol, abstol);

    ncm_quaternion_conjugate_q_mul (u, u, r);

    ncm_assert_cmpdouble_e (r->s, ==, 1.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, 0.0, reltol, abstol);

    ncm_quaternion_conjugate_q_mul (q, u, r);
    ncm_quaternion_conjugate_q_mul (u, q, s);

    ncm_quaternion_conjugate (s);

    ncm_assert_cmpdouble_e (r->s, ==, s->s, reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[0], ==, s->v.c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[1], ==, s->v.c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (r->v.c[2], ==, s->v.c[2], reltol, abstol);

    ncm_quaternion_free (q);
    ncm_quaternion_free (u);
    ncm_quaternion_free (r);
    ncm_quaternion_free (s);
  }
}

void
test_ncm_quaternion_rotate_norm (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1e-12;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmQuaternion *q = ncm_quaternion_new_from_data (g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (-100.0, 100.0),
                                                     g_test_rand_double_range (0.0, 2.0 * M_PI));
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    NcmTriVec *w = ncm_trivec_new ();

    ncm_trivec_memcpy (w, v);
    ncm_quaternion_rotate (q, v);

    ncm_assert_cmpdouble_e (ncm_trivec_norm (v), ==, ncm_trivec_norm (w), reltol, abstol);

    ncm_quaternion_free (q);
    ncm_trivec_free (v);
    ncm_trivec_free (w);
  }
}

void
test_ncm_quaternion_rotate_spherical (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-10;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble psi   = g_test_rand_double_range (0.0, 2.0 * M_PI);
    const gdouble theta = g_test_rand_double_range (0.0, M_PI);
    const gdouble phi   = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmQuaternion *q0   = ncm_quaternion_new_from_data (0.0, 0.0, 1.0, psi);
    NcmQuaternion *q1   = ncm_quaternion_new_from_data (0.0, 1.0, 0.0, theta);
    NcmQuaternion *q2   = ncm_quaternion_new_from_data (0.0, 0.0, 1.0, phi);
    NcmTriVec *v        = ncm_trivec_new_full_c (0, 0, 1);
    gdouble r0, theta0, phi0;

    ncm_quaternion_rotate (q0, v);
    ncm_assert_cmpdouble_e (v->c[0], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, 1.0, reltol, abstol);

    ncm_quaternion_rotate (q1, v);

    ncm_trivec_get_spherical_coord (v, &r0, &theta0, &phi0);

    ncm_assert_cmpdouble_e (theta0, ==, theta, reltol, abstol);
    ncm_assert_cmpdouble_e (phi0, ==, 0.0, reltol, abstol);

    ncm_quaternion_rotate (q2, v);

    ncm_trivec_get_spherical_coord (v, &r0, &theta0, &phi0);

    ncm_assert_cmpdouble_e (theta0, ==, theta, reltol, abstol);
    ncm_assert_cmpdouble_e (ncm_c_radian_0_2pi (phi0), ==, phi, reltol, abstol);

    ncm_trivec_free (v);
    ncm_quaternion_free (q0);
    ncm_quaternion_free (q1);
    ncm_quaternion_free (q2);
  }
}

void
test_ncm_quaternion_rotate_spherical_onestep (void)
{
  const gdouble reltol = 1.0e-10;
  const gdouble abstol = 1.0e-10;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble psi   = g_test_rand_double_range (0.0, 2.0 * M_PI);
    const gdouble theta = g_test_rand_double_range (0.0, M_PI);
    const gdouble phi   = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmQuaternion *q0   = ncm_quaternion_new_from_data (0.0, 0.0, 1.0, psi);
    NcmQuaternion *q1   = ncm_quaternion_new_from_data (0.0, 1.0, 0.0, theta);
    NcmQuaternion *q2   = ncm_quaternion_new_from_data (0.0, 0.0, 1.0, phi);
    NcmTriVec *v        = ncm_trivec_new_full_c (0, 0, 1);
    gdouble r0, theta0, phi0;

    ncm_quaternion_rmul (q1, q0);
    ncm_quaternion_rmul (q2, q1);

    ncm_quaternion_rotate (q2, v);

    ncm_trivec_get_spherical_coord (v, &r0, &theta0, &phi0);

    ncm_assert_cmpdouble_e (theta0, ==, theta, reltol, abstol);
    ncm_assert_cmpdouble_e (ncm_c_radian_0_2pi (phi0), ==, phi, reltol, abstol);

    ncm_trivec_free (v);
    ncm_quaternion_free (q0);
    ncm_quaternion_free (q1);
    ncm_quaternion_free (q2);
  }
}

void
test_ncm_quaternion_rotate_axis (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-12;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble psi = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmTriVec *v      = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                               g_test_rand_double_range (-100.0, 100.0),
                                               g_test_rand_double_range (-100.0, 100.0));
    NcmQuaternion *q = ncm_quaternion_new_from_data (v->c[0], v->c[1], v->c[2], psi);
    NcmTriVec *w     = ncm_trivec_new ();

    ncm_trivec_memcpy (w, v);

    ncm_quaternion_rotate (q, w);

    ncm_assert_cmpdouble_e (v->c[0], ==, w->c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, w->c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, w->c[2], reltol, abstol);

    ncm_trivec_free (w);
    ncm_trivec_free (v);
    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_rotate_inverse (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-12;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    const gdouble psi = g_test_rand_double_range (0.0, 2.0 * M_PI);
    NcmTriVec *v      = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                               g_test_rand_double_range (-100.0, 100.0),
                                               g_test_rand_double_range (-100.0, 100.0));
    NcmQuaternion *q = ncm_quaternion_new_from_data (v->c[0], v->c[1], v->c[2], psi);
    NcmTriVec *w     = ncm_trivec_new ();

    ncm_trivec_memcpy (w, v);

    ncm_quaternion_rotate (q, w);
    ncm_quaternion_inv_rotate (q, w);

    ncm_assert_cmpdouble_e (v->c[0], ==, w->c[0], reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[1], ==, w->c[1], reltol, abstol);
    ncm_assert_cmpdouble_e (v->c[2], ==, w->c[2], reltol, abstol);

    ncm_trivec_free (v);
    ncm_trivec_free (w);
    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_set_to_rotate_to_x (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-10;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    NcmTriVec *w = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    NcmQuaternion *q = ncm_quaternion_new ();
    NcmTriVec *u     = ncm_trivec_new ();

    ncm_quaternion_set_to_rotate_to_x (q, v);
    ncm_trivec_memcpy (u, w);

    ncm_quaternion_rotate (q, u);
    /* Preserve the dot product */
    ncm_assert_cmpdouble_e (ncm_trivec_dot (v, w), ==, u->c[0] * ncm_trivec_norm (v), reltol, abstol);

    ncm_trivec_memcpy (u, v);
    ncm_quaternion_rotate (q, u);
    ncm_assert_cmpdouble_e (u->c[0], ==, ncm_trivec_norm (v), reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[2], ==, 0.0, reltol, abstol);

    ncm_trivec_free (v);
    ncm_trivec_free (w);
    ncm_trivec_free (u);
    ncm_quaternion_free (q);
  }
}

void
test_ncm_quaternion_set_to_rotate_to_z (void)
{
  const gdouble reltol = 1.0e-15;
  const gdouble abstol = 1.0e-10;
  gint i;

  for (i = 0; i < NTESTS; i++)
  {
    NcmTriVec *v = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    NcmTriVec *w = ncm_trivec_new_full_c (g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0),
                                          g_test_rand_double_range (-100.0, 100.0));
    NcmQuaternion *q = ncm_quaternion_new ();
    NcmTriVec *u     = ncm_trivec_new ();

    ncm_quaternion_set_to_rotate_to_z (q, v);
    ncm_trivec_memcpy (u, w);

    ncm_quaternion_rotate (q, u);
    /* Preserve the dot product */
    ncm_assert_cmpdouble_e (ncm_trivec_dot (v, w), ==, u->c[2] * ncm_trivec_norm (v), reltol, abstol);

    ncm_trivec_memcpy (u, v);
    ncm_quaternion_rotate (q, u);
    ncm_assert_cmpdouble_e (u->c[0], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[1], ==, 0.0, reltol, abstol);
    ncm_assert_cmpdouble_e (u->c[2], ==, ncm_trivec_norm (v), reltol, abstol);

    ncm_trivec_free (v);
    ncm_trivec_free (w);
    ncm_trivec_free (u);
    ncm_quaternion_free (q);
  }
}


/***************************************************************************
 *            test_nc_halo_position.c
 *
 *  Fri Aug 16 07:18:13 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <caiolimadeoliveira@pm.me>
 * Copyright (C) Sandro Dias Pinto Vitenti 2024 <vitenti@uel.br>
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

typedef struct _TestNcHaloPosition
{
  NcHaloPosition *hp;
} TestNcHaloPosition;

static void test_nc_halo_position_new (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_free (TestNcHaloPosition *test, gconstpointer pdata);

static void test_nc_halo_position_serialize (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_model_id (TestNcHaloPosition *test, gconstpointer pdata);

static void test_nc_halo_position_polar_angles (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_rotate (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_projected_radius (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_shared_distance (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_projected_radius_array (TestNcHaloPosition *test, gconstpointer pdata);
static void test_nc_halo_position_projected_radius_array_len (void);
static void test_nc_halo_position_projected_radius_array_len_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/halo_position/serialize", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_serialize,
              &test_nc_halo_position_free);

  g_test_add ("/nc/halo_position/model_id", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_model_id,
              &test_nc_halo_position_free);

  g_test_add ("/nc/halo_position/polar_angles", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_polar_angles,
              &test_nc_halo_position_free);

  g_test_add ("/nc/halo_position/rotate", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_rotate,
              &test_nc_halo_position_free);

  g_test_add ("/nc/halo_position/projected_radius", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_projected_radius,
              &test_nc_halo_position_free);

  g_test_add ("/nc/halo_position/shared_distance", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_shared_distance,
              &test_nc_halo_position_free);

  g_test_add ("/nc/halo_position/projected_radius_array", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_projected_radius_array,
              &test_nc_halo_position_free);

  g_test_add_func ("/nc/halo_position/projected_radius_array/len", &test_nc_halo_position_projected_radius_array_len);
  g_test_add_func ("/nc/halo_position/projected_radius_array/len/subprocess", &test_nc_halo_position_projected_radius_array_len_subprocess);

  g_test_run ();

  return 0;
}

static void
test_nc_halo_position_new (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcDistance *dist = nc_distance_new (1100.0);
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  test->hp = nc_halo_position_new (dist);

  g_assert_true (NC_IS_HALO_POSITION (test->hp));

  nc_halo_position_prepare (test->hp, cosmo);

  nc_distance_clear (&dist);
  nc_hicosmo_clear (&cosmo);
}

static void
test_nc_halo_position_free (TestNcHaloPosition *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_halo_position_free, test->hp);
}

static void
test_nc_halo_position_serialize (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcHaloPosition *hp     = test->hp;
  NcmSerialize *ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *hp_ser          = ncm_serialize_to_string (ser, G_OBJECT (hp), TRUE);
  NcHaloPosition *hp_dup = NC_HALO_POSITION (ncm_serialize_from_string (ser, hp_ser));
  NcHICosmo *cosmo       = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  gdouble ra, dec, theta, phi, r, theta_dup, phi_dup, r_dup;

  nc_halo_position_prepare_if_needed (hp_dup, cosmo);

  g_assert_true (NC_IS_HALO_POSITION (hp_dup));

  ra  = g_test_rand_double_range (0.0, 2.0 * M_PI);
  dec = g_test_rand_double_range (0.0, 2.0 * M_PI);

  nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);
  nc_halo_position_polar_angles (hp_dup, ra, dec, &theta_dup, &phi_dup);

  g_assert_cmpfloat (theta, ==, theta_dup);
  g_assert_cmpfloat (phi, ==, phi_dup);

  r     = nc_halo_position_projected_radius (hp, cosmo, theta);
  r_dup = nc_halo_position_projected_radius (hp_dup, cosmo, theta);

  g_assert_cmpfloat (r, ==, r_dup);

  ncm_serialize_clear (&ser);
  g_free (hp_ser);
  nc_halo_position_clear (&hp_dup);
  nc_hicosmo_clear (&cosmo);
}

static void
test_nc_halo_position_model_id (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcmMSet *model_set = ncm_mset_empty_new ();
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  ncm_mset_set (model_set, ncm_model_dup (NCM_MODEL (test->hp), ser), NULL);

  g_assert_true (NC_IS_HALO_POSITION (ncm_mset_peek (model_set, nc_halo_position_id ())));

  ncm_mset_clear (&model_set);
  ncm_serialize_clear (&ser);
}

static void
test_nc_halo_position_polar_angles (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcHaloPosition *hp = test->hp;
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  gdouble ra, dec, theta, phi;

  ra  = g_test_rand_double_range (0.0, 2.0 * M_PI);
  dec = g_test_rand_double_range (0.0, 2.0 * M_PI);

  nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);

  g_assert_cmpfloat (theta, >=, 0.0);
  g_assert_cmpfloat (theta, <=, M_PI);

  g_assert_cmpfloat (phi, >=, 0.0);
  g_assert_cmpfloat (phi, <=, 2.0 * M_PI);

  nc_hicosmo_clear (&cosmo);
}

/*
 * nc_halo_position_rotate() is a rotation that takes the halo to the positive z axis,
 * checked without nc_halo_position_polar_angles() (which is implemented with it): the
 * halo goes to z, the polar angle of a rotated position is the angular separation of
 * the spherical law of cosines, and lengths, the angle between two positions and their
 * handedness about the halo are kept. The angles of a rotated position are then those of
 * nc_halo_position_polar_angles().
 */
static void
test_nc_halo_position_rotate (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcHaloPosition *hp = test->hp;
  guint i;

  for (i = 0; i < 100; i++)
  {
    const gdouble ra_h  = g_test_rand_double_range (0.0, 360.0);
    const gdouble dec_h = g_test_rand_double_range (-80.0, 80.0);
    const gdouble ra    = ra_h + g_test_rand_double_range (-1.0, 1.0);
    const gdouble dec   = dec_h + g_test_rand_double_range (-1.0, 1.0);
    const gdouble ra2   = ra_h + g_test_rand_double_range (-1.0, 1.0);
    const gdouble dec2  = dec_h + g_test_rand_double_range (-1.0, 1.0);
    const gdouble d2r   = M_PI / 180.0;
    const gdouble cos_s = sin (dec_h * d2r) * sin (dec * d2r) + cos (dec_h * d2r) * cos (dec * d2r) * cos ((ra - ra_h) * d2r);
    gdouble theta, phi, r, theta_v, phi_v;
    NcmTriVec h, v, w, v0, w0;

    ncm_model_param_set_by_name (NCM_MODEL (hp), "ra", ra_h, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (hp), "dec", dec_h, NULL);

    /* The halo goes to the positive z axis */
    ncm_trivec_set_astro_ra_dec (&h, 1.0, ra_h, dec_h);
    nc_halo_position_rotate (hp, &h);
    g_assert_cmpfloat (fabs (h.c[0]), <, 1.0e-14);
    g_assert_cmpfloat (fabs (h.c[1]), <, 1.0e-14);
    g_assert_cmpfloat (fabs (h.c[2] - 1.0), <, 1.0e-14);

    /* The polar angle of a position is its angular separation from the halo */
    ncm_trivec_set_astro_ra_dec (&v0, 1.0, ra, dec);
    ncm_trivec_set_astro_ra_dec (&w0, 1.0, ra2, dec2);
    v = v0;
    w = w0;
    nc_halo_position_rotate (hp, &v);
    nc_halo_position_rotate (hp, &w);
    ncm_assert_cmpdouble_e (v.c[2], ==, cos_s, 1.0e-13, 0.0);

    /* A rotation: lengths, angles and handedness are kept */
    ncm_assert_cmpdouble_e (ncm_trivec_norm (&v), ==, 1.0, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_trivec_dot (&v, &w), ==, ncm_trivec_dot (&v0, &w0), 1.0e-13, 0.0);
    {
      const gdouble cross_z = v.c[0] * w.c[1] - v.c[1] * w.c[0];
      NcmTriVec hh;
      gdouble triple;

      ncm_trivec_set_astro_ra_dec (&hh, 1.0, ra_h, dec_h);
      triple = hh.c[0] * (v0.c[1] * w0.c[2] - v0.c[2] * w0.c[1]) -
               hh.c[1] * (v0.c[0] * w0.c[2] - v0.c[2] * w0.c[0]) +
               hh.c[2] * (v0.c[0] * w0.c[1] - v0.c[1] * w0.c[0]);
      ncm_assert_cmpdouble_e (cross_z, ==, triple, 1.0e-10, 1.0e-14);
    }

    /* The angles of a rotated position are those of nc_halo_position_polar_angles() */
    nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);
    ncm_trivec_get_spherical_coord (&v, &r, &theta_v, &phi_v);

    g_assert_cmpfloat (theta_v, ==, theta);
    g_assert_cmpfloat (phi_v, ==, phi);
  }
}

static void
test_nc_halo_position_projected_radius (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcHaloPosition *hp = test->hp;
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  gdouble ra, dec, theta, phi, r;

  ra  = g_test_rand_double_range (0.0, 2.0 * M_PI);
  dec = g_test_rand_double_range (0.0, 2.0 * M_PI);

  nc_halo_position_polar_angles (hp, ra, dec, &theta, &phi);

  r = nc_halo_position_projected_radius (hp, cosmo, theta);

  g_assert_cmpfloat (r, >=, 0.0);

  nc_hicosmo_clear (&cosmo);
}

static void
test_nc_halo_position_shared_distance (TestNcHaloPosition *test, gconstpointer pdata)
{
  /* nc_halo_position_prepare() has no local state: it only delegates to
   * nc_distance_prepare_if_needed(), which carries its own control. So it must
   * not mark itself current -- doing so lets prepare_if_needed() skip the
   * delegation, and a NcDistance shared with a second cosmology then stays
   * driven by that other cosmology. Regression for PR #318. */
  NcDistance *dist   = nc_distance_new (1100.0);
  NcHaloPosition *hp = nc_halo_position_new (dist);
  NcHICosmo *cosmo_a = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcHICosmo *cosmo_b = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  const gdouble z    = 1.0;
  gdouble d_a, d_after;

  ncm_model_param_set_by_name (NCM_MODEL (cosmo_a), "Omegac", 0.25, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo_b), "Omegac", 0.40, NULL);

  nc_halo_position_prepare (hp, cosmo_a);
  d_a = nc_distance_comoving (dist, cosmo_a, z);

  /* Another consumer drives the shared distance to the second cosmology. */
  nc_distance_prepare_if_needed (dist, cosmo_b);

  /* This call must restore it. */
  nc_halo_position_prepare_if_needed (hp, cosmo_a);
  d_after = nc_distance_comoving (dist, cosmo_a, z);

  ncm_assert_cmpdouble_e (d_after, ==, d_a, 1.0e-14, 0.0);

  nc_halo_position_free (hp);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo_a);
  nc_hicosmo_free (cosmo_b);
}

/* The array of projected radii equals the radius of each position, and a reused array is resized */
static void
test_nc_halo_position_projected_radius_array (TestNcHaloPosition *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  GArray *ra       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *dec      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *r        = NULL;
  const guint n    = 20;
  guint i;

  ncm_model_param_set (NCM_MODEL (test->hp), NC_HALO_POSITION_RA, 30.0);
  ncm_model_param_set (NCM_MODEL (test->hp), NC_HALO_POSITION_DEC, -5.0);
  ncm_model_param_set (NCM_MODEL (test->hp), NC_HALO_POSITION_Z, 0.3);
  nc_halo_position_prepare (test->hp, cosmo);

  for (i = 0; i < n; i++)
  {
    const gdouble ra_i  = 30.0 + g_test_rand_double_range (-0.5, 0.5);
    const gdouble dec_i = -5.0 + g_test_rand_double_range (-0.5, 0.5);

    g_array_append_val (ra, ra_i);
    g_array_append_val (dec, dec_i);
  }

  nc_halo_position_projected_radius_array (test->hp, cosmo, ra, dec, &r);
  g_assert_cmpuint (r->len, ==, n);

  for (i = 0; i < n; i++)
    g_assert_cmpfloat (g_array_index (r, gdouble, i), ==,
                       nc_halo_position_projected_radius_from_ra_dec (test->hp, cosmo, g_array_index (ra, gdouble, i), g_array_index (dec, gdouble, i)));

  g_array_set_size (ra, 3);
  g_array_set_size (dec, 3);
  nc_halo_position_projected_radius_array (test->hp, cosmo, ra, dec, &r);
  g_assert_cmpuint (r->len, ==, 3);

  g_array_unref (ra);
  g_array_unref (dec);
  g_array_unref (r);
  nc_hicosmo_clear (&cosmo);
}

/* Arrays of right ascensions and declinations of different lengths abort */
static void
test_nc_halo_position_projected_radius_array_len (void)
{
  g_test_trap_subprocess ("/nc/halo_position/projected_radius_array/len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*2 right ascensions and 1 declinations*");
}

static void
test_nc_halo_position_projected_radius_array_len_subprocess (void)
{
  NcDistance *dist   = nc_distance_new (1100.0);
  NcHaloPosition *hp = nc_halo_position_new (dist);
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  GArray *ra         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *dec        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *r          = NULL;
  const gdouble x    = 0.1;

  nc_halo_position_prepare (hp, cosmo);
  g_array_append_val (ra, x);
  g_array_append_val (ra, x);
  g_array_append_val (dec, x);

  nc_halo_position_projected_radius_array (hp, cosmo, ra, dec, &r);
}


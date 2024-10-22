/***************************************************************************
 *            test_nc_halo_position.c
 *
 *  Fri Aug 16 07:18:13 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <code.caio@limadeoliveira.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Caio Lima de Oliveira 2024 <code.caio@limadeoliveira.me>
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
static void test_nc_halo_position_projected_radius (TestNcHaloPosition *test, gconstpointer pdata);

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

  g_test_add ("/nc/halo_position/projected_radius", TestNcHaloPosition, NULL,
              &test_nc_halo_position_new,
              &test_nc_halo_position_projected_radius,
              &test_nc_halo_position_free);

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


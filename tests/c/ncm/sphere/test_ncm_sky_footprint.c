/***************************************************************************
 *            test_ncm_sky_footprint.c
 *
 *  Sun Jun 14 10:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sky_footprint.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include <gsl/gsl_math.h>

typedef struct _TestNcmSkyFootprint
{
  NcmSkyFootprintRectangular *rect;
} TestNcmSkyFootprint;

void test_ncm_sky_footprint_new (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_free (TestNcmSkyFootprint *test, gconstpointer pdata);

void test_ncm_sky_footprint_ref (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_limits (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_area (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_contains (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_density (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_gen (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_serialize (TestNcmSkyFootprint *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/sky_footprint/ref", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_ref, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/limits", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_limits, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/area", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_area, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/contains", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_contains, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/density", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_density, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/gen", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_gen, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/serialize", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_serialize, &test_ncm_sky_footprint_free);

  g_test_run ();

  return 0;
}

void
test_ncm_sky_footprint_new (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  test->rect = ncm_sky_footprint_rectangular_new (10.0, 40.0, -5.0, 25.0);

  g_assert_true (NCM_IS_SKY_FOOTPRINT_RECTANGULAR (test->rect));
  g_assert_true (NCM_IS_SKY_FOOTPRINT (test->rect));
}

void
test_ncm_sky_footprint_free (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_sky_footprint_rectangular_free, test->rect);
}

void
test_ncm_sky_footprint_ref (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSkyFootprintRectangular *rect_ref = ncm_sky_footprint_rectangular_ref (test->rect);

  g_assert_true (rect_ref == test->rect);

  ncm_sky_footprint_rectangular_clear (&rect_ref);
  g_assert_null (rect_ref);

  g_assert_true (NCM_IS_SKY_FOOTPRINT_RECTANGULAR (test->rect));
}

void
test_ncm_sky_footprint_limits (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  gdouble ra_min, ra_max, dec_min, dec_max;

  ncm_sky_footprint_rectangular_get_ra_lim (test->rect, &ra_min, &ra_max);
  ncm_sky_footprint_rectangular_get_dec_lim (test->rect, &dec_min, &dec_max);

  g_assert_cmpfloat (ra_min, ==, 10.0);
  g_assert_cmpfloat (ra_max, ==, 40.0);
  g_assert_cmpfloat (dec_min, ==, -5.0);
  g_assert_cmpfloat (dec_max, ==, 25.0);
}

void
test_ncm_sky_footprint_area (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  const gdouble expected = ncm_c_degree_to_radian (40.0 - 10.0) *
                           (sin (ncm_c_degree_to_radian (25.0)) - sin (ncm_c_degree_to_radian (-5.0)));

  g_assert_cmpfloat (ncm_sky_footprint_get_area (NCM_SKY_FOOTPRINT (test->rect)), ==, expected);
}

void
test_ncm_sky_footprint_contains (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSkyFootprint *fp = NCM_SKY_FOOTPRINT (test->rect);

  g_assert_true (ncm_sky_footprint_contains (fp, 25.0, 10.0));
  g_assert_false (ncm_sky_footprint_contains (fp, 100.0, 0.0));
  g_assert_false (ncm_sky_footprint_contains (fp, 25.0, 80.0));
}

void
test_ncm_sky_footprint_density (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSkyFootprint *fp     = NCM_SKY_FOOTPRINT (test->rect);
  const gdouble inside    = ncm_sky_footprint_density (fp, 25.0, 10.0);
  const gdouble ln_inside = ncm_sky_footprint_ln_density (fp, 25.0, 10.0);

  g_assert_cmpfloat (inside, >, 0.0);
  g_assert_cmpfloat (ln_inside, ==, log (inside));

  g_assert_cmpfloat (ncm_sky_footprint_density (fp, 100.0, 0.0), ==, 0.0);
  g_assert_true (gsl_isinf (ncm_sky_footprint_ln_density (fp, 100.0, 0.0)) == -1);
}

void
test_ncm_sky_footprint_gen (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSkyFootprint *fp = NCM_SKY_FOOTPRINT (test->rect);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 123);
  guint i;

  for (i = 0; i < 1000; i++)
  {
    gdouble ra, dec;

    ncm_sky_footprint_gen_ra_dec (fp, rng, &ra, &dec);

    g_assert_cmpfloat (ra, >=, 10.0);
    g_assert_cmpfloat (ra, <=, 40.0);
    g_assert_cmpfloat (dec, >=, -5.0);
    g_assert_cmpfloat (dec, <=, 25.0);
    g_assert_true (ncm_sky_footprint_contains (fp, ra, dec));
  }

  ncm_rng_free (rng);
}

void
test_ncm_sky_footprint_serialize (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSkyFootprintRectangular *dup;
  gdouble ra_min, ra_max, dec_min, dec_max;

  dup = NCM_SKY_FOOTPRINT_RECTANGULAR (ncm_serialize_dup_obj (ser, G_OBJECT (test->rect)));

  ncm_sky_footprint_rectangular_get_ra_lim (dup, &ra_min, &ra_max);
  ncm_sky_footprint_rectangular_get_dec_lim (dup, &dec_min, &dec_max);

  g_assert_cmpfloat (ra_min, ==, 10.0);
  g_assert_cmpfloat (ra_max, ==, 40.0);
  g_assert_cmpfloat (dec_min, ==, -5.0);
  g_assert_cmpfloat (dec_max, ==, 25.0);

  ncm_sky_footprint_rectangular_free (dup);
  ncm_serialize_free (ser);
}

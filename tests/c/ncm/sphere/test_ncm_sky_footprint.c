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
void test_ncm_sky_footprint_gen_moments (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_gen_truth_table (TestNcmSkyFootprint *test, gconstpointer pdata);
void test_ncm_sky_footprint_default (void);
void test_ncm_sky_footprint_traps (void);
void test_ncm_sky_footprint_invalid_ra_span (void);
void test_ncm_sky_footprint_invalid_ra_wide (void);
void test_ncm_sky_footprint_invalid_dec_range (void);
void test_ncm_sky_footprint_invalid_dec_order (void);

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
  g_test_add ("/ncm/sky_footprint/gen/moments", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_gen_moments, &test_ncm_sky_footprint_free);
  g_test_add ("/ncm/sky_footprint/gen/truth_table", TestNcmSkyFootprint, NULL,
              &test_ncm_sky_footprint_new, &test_ncm_sky_footprint_gen_truth_table, &test_ncm_sky_footprint_free);
  g_test_add_func ("/ncm/sky_footprint/default", &test_ncm_sky_footprint_default);
  g_test_add_func ("/ncm/sky_footprint/traps", &test_ncm_sky_footprint_traps);
  g_test_add_func ("/ncm/sky_footprint/invalid/ra_span/subprocess", &test_ncm_sky_footprint_invalid_ra_span);
  g_test_add_func ("/ncm/sky_footprint/invalid/ra_wide/subprocess", &test_ncm_sky_footprint_invalid_ra_wide);
  g_test_add_func ("/ncm/sky_footprint/invalid/dec_range/subprocess", &test_ncm_sky_footprint_invalid_dec_range);
  g_test_add_func ("/ncm/sky_footprint/invalid/dec_order/subprocess", &test_ncm_sky_footprint_invalid_dec_order);

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

  const gdouble dsin    = sin (ncm_c_degree_to_radian (25.0)) - sin (ncm_c_degree_to_radian (-5.0));
  const gdouble dec_a[] = {-5.0, 0.0, 10.0, 24.9};
  guint i;

  g_assert_cmpfloat (ln_inside, ==, log (inside));

  /* The density in dra ddec (degrees): cos(dec) / (dra (180 / pi) dsin), which integrates to
   * one over the rectangle. */
  for (i = 0; i < G_N_ELEMENTS (dec_a); i++)
  {
    const gdouble truth = cos (ncm_c_degree_to_radian (dec_a[i])) / ((40.0 - 10.0) * (180.0 / M_PI) * dsin);

    ncm_assert_cmpdouble_e (ncm_sky_footprint_density (fp, 25.0, dec_a[i]), ==, truth, 1.0e-15, 0.0);
  }

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

/*
 * Uniform on the sphere: ra uniform in [10, 40] and sin(dec) uniform in
 * [sin(-5), sin(25)]. The sample means and variances of both, from 200000 draws of a fixed
 * seed, must lie within 4 standard errors of the uniform values.
 */
void
test_ncm_sky_footprint_gen_moments (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSkyFootprint *fp = NCM_SKY_FOOTPRINT (test->rect);
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 20260927);
  const guint n       = 200000;
  const gdouble lo[2] = {10.0, sin (ncm_c_degree_to_radian (-5.0))};
  const gdouble hi[2] = {40.0, sin (ncm_c_degree_to_radian (25.0))};
  gdouble sum[2]      = {0.0, 0.0};
  gdouble sum2[2]     = {0.0, 0.0};
  guint i, c;

  for (i = 0; i < n; i++)
  {
    gdouble ra, dec, v[2];

    ncm_sky_footprint_gen_ra_dec (fp, rng, &ra, &dec);
    v[0] = ra;
    v[1] = sin (ncm_c_degree_to_radian (dec));

    for (c = 0; c < 2; c++)
    {
      const gdouble x = (v[c] - lo[c]) / (hi[c] - lo[c]);

      sum[c]  += x;
      sum2[c] += x * x;
    }
  }

  for (c = 0; c < 2; c++)
  {
    const gdouble mean = sum[c] / n;
    const gdouble var  = sum2[c] / n - mean * mean;

    /* Uniform on [0, 1]: mean 1/2, variance 1/12, standard errors sqrt(1/12/n) and
     * sqrt((1/80 - 1/144)/n). */
    g_assert_cmpfloat (fabs (mean - 0.5), <, 4.0 * sqrt (1.0 / 12.0 / n));
    g_assert_cmpfloat (fabs (var - 1.0 / 12.0), <, 4.0 * sqrt ((1.0 / 80.0 - 1.0 / 144.0) / n));
  }

  ncm_rng_free (rng);
}

/*
 * The draws of seed 123 against a stored reference (NcmMatrix N x 2 of ra, dec), to a
 * tolerance that survives the sub-ULP drift of asin across libm builds but catches any
 * change in the sampling.
 */
void
test_ncm_sky_footprint_gen_truth_table (TestNcmSkyFootprint *test, gconstpointer pdata)
{
  NcmSkyFootprint *fp = NCM_SKY_FOOTPRINT (test->rect);
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *path         = ncm_cfg_get_data_filename ("truth_tables/sphere/ncm_sky_footprint_rect_seed123.bin", TRUE);
  NcmMatrix *truth    = NCM_MATRIX (ncm_serialize_from_binfile (ser, path));
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 123);
  guint i;

  g_assert_cmpuint (ncm_matrix_ncols (truth), ==, 2);
  g_assert_cmpuint (ncm_matrix_nrows (truth), ==, 5000);

  for (i = 0; i < ncm_matrix_nrows (truth); i++)
  {
    gdouble ra, dec;

    ncm_sky_footprint_gen_ra_dec (fp, rng, &ra, &dec);
    g_assert_cmpfloat (fabs (ra - ncm_matrix_get (truth, i, 0)), <=, 1.0e-12 + 1.0e-9 * fabs (ncm_matrix_get (truth, i, 0)));
    g_assert_cmpfloat (fabs (dec - ncm_matrix_get (truth, i, 1)), <=, 1.0e-12 + 1.0e-9 * fabs (ncm_matrix_get (truth, i, 1)));
  }

  ncm_rng_free (rng);
  ncm_matrix_free (truth);
  g_free (path);
  ncm_serialize_free (ser);
}

/* A limit not given spans the whole sphere in that coordinate. */
void
test_ncm_sky_footprint_default (void)
{
  NcmSkyFootprintRectangular *full = g_object_new (NCM_TYPE_SKY_FOOTPRINT_RECTANGULAR, NULL);
  NcmDTuple2 ra_lim                = NCM_DTUPLE2_STATIC_INIT (10.0, 40.0);
  NcmSkyFootprintRectangular *band = g_object_new (NCM_TYPE_SKY_FOOTPRINT_RECTANGULAR, "ra-lim", &ra_lim, NULL);
  gdouble ra_min, ra_max, dec_min, dec_max;

  ncm_sky_footprint_rectangular_get_ra_lim (full, &ra_min, &ra_max);
  ncm_sky_footprint_rectangular_get_dec_lim (full, &dec_min, &dec_max);
  g_assert_cmpfloat (ra_min, ==, 0.0);
  g_assert_cmpfloat (ra_max, ==, 360.0);
  g_assert_cmpfloat (dec_min, ==, -90.0);
  g_assert_cmpfloat (dec_max, ==, 90.0);
  ncm_assert_cmpdouble_e (ncm_sky_footprint_get_area (NCM_SKY_FOOTPRINT (full)), ==, 4.0 * M_PI, 1.0e-15, 0.0);

  ncm_sky_footprint_rectangular_get_dec_lim (band, &dec_min, &dec_max);
  g_assert_cmpfloat (dec_min, ==, -90.0);
  g_assert_cmpfloat (dec_max, ==, 90.0);
  ncm_assert_cmpdouble_e (ncm_sky_footprint_get_area (NCM_SKY_FOOTPRINT (band)), ==, 2.0 * ncm_c_degree_to_radian (30.0), 1.0e-15, 0.0);

  ncm_sky_footprint_rectangular_free (full);
  ncm_sky_footprint_rectangular_free (band);
}

void
test_ncm_sky_footprint_traps (void)
{
  g_test_trap_subprocess ("/ncm/sky_footprint/invalid/ra_span/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*span must be in (0, 360]*");

  g_test_trap_subprocess ("/ncm/sky_footprint/invalid/ra_wide/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*span must be in (0, 360]*");

  g_test_trap_subprocess ("/ncm/sky_footprint/invalid/dec_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*-90 <= min < max <= 90*");

  g_test_trap_subprocess ("/ncm/sky_footprint/invalid/dec_order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*-90 <= min < max <= 90*");
}

void
test_ncm_sky_footprint_invalid_ra_span (void)
{
  ncm_sky_footprint_rectangular_new (40.0, 10.0, -5.0, 25.0);
}

void
test_ncm_sky_footprint_invalid_ra_wide (void)
{
  ncm_sky_footprint_rectangular_new (-10.0, 360.0, -5.0, 25.0);
}

/* dec_max = 100 folds sin(dec) back and gave a wrong area and density. */
void
test_ncm_sky_footprint_invalid_dec_range (void)
{
  ncm_sky_footprint_rectangular_new (10.0, 40.0, -5.0, 100.0);
}

void
test_ncm_sky_footprint_invalid_dec_order (void)
{
  ncm_sky_footprint_rectangular_new (10.0, 40.0, 25.0, 25.0);
}


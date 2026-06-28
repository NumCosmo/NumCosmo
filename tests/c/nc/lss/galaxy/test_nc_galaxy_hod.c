/***************************************************************************
 *            test_nc_galaxy_hod.c
 *
 *  Sun Jun 14 12:30 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_hod.c
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

typedef struct _TestNcGalaxyHOD
{
  NcGalaxyHODZheng07 *zheng07;
} TestNcGalaxyHOD;

void test_nc_galaxy_hod_new (TestNcGalaxyHOD *test, gconstpointer pdata);
void test_nc_galaxy_hod_free (TestNcGalaxyHOD *test, gconstpointer pdata);

void test_nc_galaxy_hod_ref (TestNcGalaxyHOD *test, gconstpointer pdata);
void test_nc_galaxy_hod_means (TestNcGalaxyHOD *test, gconstpointer pdata);
void test_nc_galaxy_hod_gen (TestNcGalaxyHOD *test, gconstpointer pdata);
void test_nc_galaxy_hod_deterministic (TestNcGalaxyHOD *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/galaxy_hod/ref", TestNcGalaxyHOD, NULL,
              &test_nc_galaxy_hod_new, &test_nc_galaxy_hod_ref, &test_nc_galaxy_hod_free);
  g_test_add ("/nc/galaxy_hod/means", TestNcGalaxyHOD, NULL,
              &test_nc_galaxy_hod_new, &test_nc_galaxy_hod_means, &test_nc_galaxy_hod_free);
  g_test_add ("/nc/galaxy_hod/gen", TestNcGalaxyHOD, NULL,
              &test_nc_galaxy_hod_new, &test_nc_galaxy_hod_gen, &test_nc_galaxy_hod_free);
  g_test_add ("/nc/galaxy_hod/deterministic", TestNcGalaxyHOD, NULL,
              &test_nc_galaxy_hod_new, &test_nc_galaxy_hod_deterministic, &test_nc_galaxy_hod_free);

  g_test_run ();

  return 0;
}

void
test_nc_galaxy_hod_new (TestNcGalaxyHOD *test, gconstpointer pdata)
{
  test->zheng07 = nc_galaxy_hod_zheng07_new ();

  g_assert_true (NC_IS_GALAXY_HOD_ZHENG07 (test->zheng07));
  g_assert_true (NC_IS_GALAXY_HOD (test->zheng07));
}

void
test_nc_galaxy_hod_free (TestNcGalaxyHOD *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_galaxy_hod_zheng07_free, test->zheng07);
}

void
test_nc_galaxy_hod_ref (TestNcGalaxyHOD *test, gconstpointer pdata)
{
  NcGalaxyHODZheng07 *ref = nc_galaxy_hod_zheng07_ref (test->zheng07);

  g_assert_true (ref == test->zheng07);

  nc_galaxy_hod_zheng07_clear (&ref);
  g_assert_null (ref);

  g_assert_true (NC_IS_GALAXY_HOD_ZHENG07 (test->zheng07));
}

void
test_nc_galaxy_hod_means (TestNcGalaxyHOD *test, gconstpointer pdata)
{
  NcGalaxyHOD *hod = NC_GALAXY_HOD (test->zheng07);
  /* Defaults: logMmin=12.72, sigma=0.26, logM0=12.7, logM1=13.93, alpha=1.15. */
  const gdouble log10_m = 13.5;
  const gdouble lnM     = log10_m * M_LN10;
  const gdouble arg     = (log10_m - 12.72) / 0.26;
  const gdouble exp_cen = 0.5 + ncm_util_normal_gaussian_integral (0.0, arg * M_SQRT2);
  const gdouble diff    = pow (10.0, log10_m) - pow (10.0, 12.7);
  const gdouble exp_sat = pow (diff / pow (10.0, 13.93), 1.15);

  g_assert_cmpfloat (fabs (nc_galaxy_hod_mean_n_central (hod, lnM) - exp_cen), <, 1.0e-12);
  g_assert_cmpfloat (fabs (nc_galaxy_hod_mean_n_satellite (hod, lnM) - exp_sat), <, 1.0e-12);

  /* Below the satellite cutoff there are no satellites. */
  g_assert_cmpfloat (nc_galaxy_hod_mean_n_satellite (hod, 11.0 * M_LN10), ==, 0.0);
}

void
test_nc_galaxy_hod_gen (TestNcGalaxyHOD *test, gconstpointer pdata)
{
  NcGalaxyHOD *hod = NC_GALAXY_HOD (test->zheng07);
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, 42);
  const gdouble lnM = 14.0 * M_LN10;
  guint i;

  for (i = 0; i < 2000; i++)
  {
    gint n_cen, n_sat;

    nc_galaxy_hod_gen (hod, lnM, rng, &n_cen, &n_sat);

    g_assert_true (n_cen == 0 || n_cen == 1);
    g_assert_cmpint (n_sat, >=, 0);

    if (n_cen == 0)
      g_assert_cmpint (n_sat, ==, 0);
  }

  ncm_rng_free (rng);
}

void
test_nc_galaxy_hod_deterministic (TestNcGalaxyHOD *test, gconstpointer pdata)
{
  NcGalaxyHOD *hod = NC_GALAXY_HOD (test->zheng07);
  NcmRNG *rng      = ncm_rng_seeded_new (NULL, 1);
  gint low_cen, high_cen, dummy;

  nc_galaxy_hod_set_stochastic_central (hod, FALSE);
  g_assert_false (nc_galaxy_hod_get_stochastic_central (hod));

  nc_galaxy_hod_gen (hod, 11.0 * M_LN10, rng, &low_cen, &dummy);
  nc_galaxy_hod_gen (hod, 14.0 * M_LN10, rng, &high_cen, &dummy);

  g_assert_cmpint (low_cen, ==, 0);
  g_assert_cmpint (high_cen, ==, 1);

  ncm_rng_free (rng);
}

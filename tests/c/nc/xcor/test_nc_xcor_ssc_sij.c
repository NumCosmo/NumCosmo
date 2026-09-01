/***************************************************************************
 *            test_nc_xcor_ssc_sij.c
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

/*
 * #NcXcorSSCSij's configuration surface. The two footprint descriptions -- a survey area
 * and a mask spectrum -- are mutually exclusive, and which one is set decides how f_sky
 * is derived; that exclusion is the only real logic here and is what these check. The
 * S_ij matrix itself is pinned against the Python implementation it replaced, in
 * tests/python/nc/xcor/test_ssc_sij.py.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_ZMAX 6.0
#define TEST_NBINS 3

typedef struct _TestNcXcorSSCSij
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcXcorSSCSij *ssc;
} TestNcXcorSSCSij;

static void
test_nc_xcor_ssc_sij_new (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcmVector *z_edges = ncm_vector_new (TEST_NBINS + 1);
  guint i;

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,      70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, 0.255);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B, 0.045);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W, -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", 0.0, NULL);

  for (i = 0; i <= TEST_NBINS; i++)
    ncm_vector_set (z_edges, i, 0.2 + 0.3 * i);

  test->cosmo = cosmo;
  test->dist  = nc_distance_new (TEST_ZMAX);
  test->ps    = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                       NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));

  nc_distance_prepare (test->dist, cosmo);
  ncm_powspec_prepare (test->ps, NCM_MODEL (cosmo));

  test->ssc = nc_xcor_ssc_sij_new (test->dist, test->ps, z_edges);

  g_assert_true (NC_IS_XCOR_SSC_SIJ (test->ssc));

  ncm_vector_free (z_edges);
}

static void
test_nc_xcor_ssc_sij_free (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_ssc_sij_free, test->ssc);
}

static void
test_nc_xcor_ssc_sij_basic (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  NcXcorSSCSij *ssc2 = nc_xcor_ssc_sij_ref (test->ssc);

  g_assert_cmpuint (nc_xcor_ssc_sij_get_nbins (test->ssc), ==, TEST_NBINS);

  /* Constructed without a footprint, the object still holds one: the full-sky mask,
   * which is the monopole alone -- hence lmax 0 and f_sky 1. */
  g_assert_nonnull (nc_xcor_ssc_sij_peek_mask_cl (test->ssc));
  g_assert_cmpuint (nc_xcor_ssc_sij_get_lmax (test->ssc), ==, 0);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_area (test->ssc), ==, 0.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_fsky (test->ssc), ==, 1.0, 1.0e-12, 0.0);

  nc_xcor_ssc_sij_clear (&ssc2);
  g_assert_true (ssc2 == NULL);
}

static void
test_nc_xcor_ssc_sij_knobs (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  nc_xcor_ssc_sij_set_block_size (test->ssc, 7);
  g_assert_cmpuint (nc_xcor_ssc_sij_get_block_size (test->ssc), ==, 7);

  nc_xcor_ssc_sij_set_reltol (test->ssc, 1.0e-5);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_reltol (test->ssc), ==, 1.0e-5, 1.0e-15, 0.0);

  nc_xcor_ssc_sij_set_scaled_abstol (test->ssc, 1.0e-4);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_scaled_abstol (test->ssc), ==, 1.0e-4, 1.0e-15, 0.0);

  nc_xcor_ssc_sij_set_method (test->ssc, NC_XCOR_METHOD_LIMBER_Z_CUBATURE);
  g_assert_cmpuint (nc_xcor_ssc_sij_get_method (test->ssc), ==, NC_XCOR_METHOD_LIMBER_Z_CUBATURE);

  nc_xcor_ssc_sij_set_method (test->ssc, NC_XCOR_METHOD_LIMBER_Z_GSL);
  g_assert_cmpuint (nc_xcor_ssc_sij_get_method (test->ssc), ==, NC_XCOR_METHOD_LIMBER_Z_GSL);
}

/* Area and mask are two descriptions of the same footprint, and the object keeps only
 * the one set last -- setting either resets the other, so f_sky never depends on which
 * of two stored values happens to be read. The mask is never dropped to NULL; it falls
 * back to the full-sky monopole. */
static void
test_nc_xcor_ssc_sij_footprint (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  const gdouble sqd_fullsky = 4.0 * ncm_c_pi () * gsl_pow_2 (180.0 / ncm_c_pi ());
  NcmVector *quarter_sky    = ncm_vector_new (4);
  guint i;

  /* A quarter-sky mask: f_sky is sqrt (C_0 / 4pi), so the monopole is 4pi f_sky^2. */
  ncm_vector_set (quarter_sky, 0, 4.0 * ncm_c_pi () * gsl_pow_2 (0.25));

  for (i = 1; i < 4; i++)
    ncm_vector_set (quarter_sky, i, 0.0);

  nc_xcor_ssc_sij_set_area (test->ssc, 0.25 * sqd_fullsky);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_area (test->ssc), ==, 0.25 * sqd_fullsky, 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_fsky (test->ssc), ==, 0.25, 1.0e-12, 0.0);
  g_assert_cmpuint (nc_xcor_ssc_sij_get_lmax (test->ssc), ==, 0);

  /* A mask takes over: the area is zeroed, and f_sky now comes from the monopole. Both
   * routes describing a quarter sky must land on the same f_sky. */
  nc_xcor_ssc_sij_set_mask_cl (test->ssc, quarter_sky);
  g_assert_true (nc_xcor_ssc_sij_peek_mask_cl (test->ssc) == quarter_sky);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_area (test->ssc), ==, 0.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_fsky (test->ssc), ==, 0.25, 1.0e-12, 0.0);
  g_assert_cmpuint (nc_xcor_ssc_sij_get_lmax (test->ssc), ==, 3);

  /* And an area set afterwards takes over again, putting the mask back to full sky. */
  nc_xcor_ssc_sij_set_area (test->ssc, 0.5 * sqd_fullsky);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_fsky (test->ssc), ==, 0.5, 1.0e-12, 0.0);
  g_assert_cmpuint (nc_xcor_ssc_sij_get_lmax (test->ssc), ==, 0);

  /* Clearing the mask restores the full-sky monopole rather than leaving none. */
  nc_xcor_ssc_sij_set_mask_cl (test->ssc, NULL);
  g_assert_nonnull (nc_xcor_ssc_sij_peek_mask_cl (test->ssc));
  g_assert_cmpuint (nc_xcor_ssc_sij_get_lmax (test->ssc), ==, 0);

  /* Zero area means the plain full-sky matrix, not an empty sky. */
  nc_xcor_ssc_sij_set_area (test->ssc, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_area (test->ssc), ==, 0.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_ssc_sij_get_fsky (test->ssc), ==, 1.0, 1.0e-12, 0.0);

  ncm_vector_free (quarter_sky);
}

static void
test_nc_xcor_ssc_sij_traps (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/nc/xcor/ssc_sij/invalid/area/negative/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/nc/xcor/ssc_sij/invalid/area/fullsky/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

static void
test_nc_xcor_ssc_sij_invalid_area_negative (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  nc_xcor_ssc_sij_set_area (test->ssc, -1.0);
}

static void
test_nc_xcor_ssc_sij_invalid_area_fullsky (TestNcXcorSSCSij *test, gconstpointer pdata)
{
  nc_xcor_ssc_sij_set_area (test->ssc, 4.0 * ncm_c_pi () * gsl_pow_2 (180.0 / ncm_c_pi ()));
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorSSCSij *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"basic",     &test_nc_xcor_ssc_sij_basic    },
  {"knobs",     &test_nc_xcor_ssc_sij_knobs    },
  {"footprint", &test_nc_xcor_ssc_sij_footprint},
  {"traps",     &test_nc_xcor_ssc_sij_traps    },
};

gint
main (gint argc, gchar *argv[])
{
  guint j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
  {
    gchar *path = g_strdup_printf ("/nc/xcor/ssc_sij/%s", test_checks[j].name);

    g_test_add (path, TestNcXcorSSCSij, NULL,
                &test_nc_xcor_ssc_sij_new, test_checks[j].func, &test_nc_xcor_ssc_sij_free);

    g_free (path);
  }

  g_test_add ("/nc/xcor/ssc_sij/invalid/area/negative/subprocess", TestNcXcorSSCSij, NULL,
              &test_nc_xcor_ssc_sij_new, &test_nc_xcor_ssc_sij_invalid_area_negative,
              &test_nc_xcor_ssc_sij_free);
  g_test_add ("/nc/xcor/ssc_sij/invalid/area/fullsky/subprocess", TestNcXcorSSCSij, NULL,
              &test_nc_xcor_ssc_sij_new, &test_nc_xcor_ssc_sij_invalid_area_fullsky,
              &test_nc_xcor_ssc_sij_free);

  g_test_run ();
}


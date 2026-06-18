/***************************************************************************
 *            test_nc_halo_catalog_member_generator.c
 *
 *  Sun Jun 14 14:30 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_halo_catalog_member_generator.c
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

typedef struct _TestNcHaloCatalogMemberGenerator
{
  NcGalaxyHOD *hod;
  NcHaloCatalogMemberGenerator *memgen;
  NcmMSet *mset;
} TestNcHaloCatalogMemberGenerator;

void test_nc_halo_catalog_member_generator_new (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_member_generator_free (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);

void test_nc_halo_catalog_member_generator_ref (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_member_generator_hod (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_member_generator_distance (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_member_generator_generate (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);
void test_nc_halo_catalog_member_generator_missing_column (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/halo_catalog_member_generator/ref", TestNcHaloCatalogMemberGenerator, NULL,
              &test_nc_halo_catalog_member_generator_new,
              &test_nc_halo_catalog_member_generator_ref,
              &test_nc_halo_catalog_member_generator_free);

  g_test_add ("/nc/halo_catalog_member_generator/hod", TestNcHaloCatalogMemberGenerator, NULL,
              &test_nc_halo_catalog_member_generator_new,
              &test_nc_halo_catalog_member_generator_hod,
              &test_nc_halo_catalog_member_generator_free);

  g_test_add ("/nc/halo_catalog_member_generator/distance", TestNcHaloCatalogMemberGenerator, NULL,
              &test_nc_halo_catalog_member_generator_new,
              &test_nc_halo_catalog_member_generator_distance,
              &test_nc_halo_catalog_member_generator_free);

  g_test_add ("/nc/halo_catalog_member_generator/generate", TestNcHaloCatalogMemberGenerator, NULL,
              &test_nc_halo_catalog_member_generator_new,
              &test_nc_halo_catalog_member_generator_generate,
              &test_nc_halo_catalog_member_generator_free);

  g_test_add ("/nc/halo_catalog_member_generator/missing_column", TestNcHaloCatalogMemberGenerator, NULL,
              &test_nc_halo_catalog_member_generator_new,
              &test_nc_halo_catalog_member_generator_missing_column,
              &test_nc_halo_catalog_member_generator_free);

  g_test_run ();

  return 0;
}

static NcHaloCatalog *
_host_catalog (void)
{
  const gchar *col_names[] = {"ra", "dec", "z_true", "lnM_true", "r_Delta", NULL};
  NcHaloCatalog *host      = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_HALO, NULL, NULL, 1, (GStrv) col_names, NULL, 0);
  NcmCatalog *hostc        = NCM_CATALOG (host);

  ncm_catalog_set (hostc, "ra", 0, 100.0, NULL);
  ncm_catalog_set (hostc, "dec", 0, 20.0, NULL);
  ncm_catalog_set (hostc, "z_true", 0, 0.3, NULL);
  ncm_catalog_set (hostc, "lnM_true", 0, 15.0 * M_LN10, NULL);
  ncm_catalog_set (hostc, "r_Delta", 0, 2.0, NULL);

  return host;
}

void
test_nc_halo_catalog_member_generator_new (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  test->hod    = NC_GALAXY_HOD (nc_galaxy_hod_zheng07_new ());
  test->memgen = nc_halo_catalog_member_generator_new (test->hod);
  test->mset   = ncm_mset_new (cosmo, NULL, NULL);

  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "H0", 70.0, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegac", 0.25, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegax", 0.70, NULL);

  g_assert_true (NC_IS_HALO_CATALOG_MEMBER_GENERATOR (test->memgen));

  nc_hicosmo_free (cosmo);
}

void
test_nc_halo_catalog_member_generator_free (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  ncm_mset_free (test->mset);
  nc_galaxy_hod_free (test->hod);
  NCM_TEST_FREE (nc_halo_catalog_member_generator_free, test->memgen);
}

void
test_nc_halo_catalog_member_generator_ref (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  NcHaloCatalogMemberGenerator *ref = nc_halo_catalog_member_generator_ref (test->memgen);

  g_assert_true (ref == test->memgen);

  nc_halo_catalog_member_generator_clear (&ref);
  g_assert_null (ref);

  g_assert_true (NC_IS_HALO_CATALOG_MEMBER_GENERATOR (test->memgen));
}

void
test_nc_halo_catalog_member_generator_hod (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  g_assert_true (nc_halo_catalog_member_generator_peek_hod (test->memgen) == test->hod);
}

void
test_nc_halo_catalog_member_generator_distance (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  NcDistance *dist = nc_distance_new (3.0);

  g_assert_null (nc_halo_catalog_member_generator_peek_distance (test->memgen));

  nc_halo_catalog_member_generator_set_distance (test->memgen, dist);
  g_assert_true (nc_halo_catalog_member_generator_peek_distance (test->memgen) == dist);

  nc_halo_catalog_member_generator_set_distance (test->memgen, NULL);
  g_assert_null (nc_halo_catalog_member_generator_peek_distance (test->memgen));

  nc_distance_free (dist);
}

void
test_nc_halo_catalog_member_generator_generate (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  NcHaloCatalog *host = _host_catalog ();
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 7);
  GError *error       = NULL;
  NcHaloCatalog *members;
  NcmCatalog *memc;
  guint n_central = 0;
  guint i;

  nc_galaxy_hod_set_stochastic_central (test->hod, FALSE);

  members = nc_halo_catalog_member_generator_generate (test->memgen, host, test->mset, rng, &error);
  g_assert_no_error (error);
  g_assert_nonnull (members);

  g_assert_cmpint (nc_halo_catalog_get_kind (members), ==, NC_HALO_CATALOG_KIND_MEMBER);

  memc = NCM_CATALOG (members);
  g_assert_cmpuint (ncm_catalog_len (memc), >, 0);

  for (i = 0; i < ncm_catalog_len (memc); i++)
  {
    const gboolean is_central = ncm_catalog_get_bool (memc, "is_central", i, NULL);

    /* Every member links back to the only host. */
    g_assert_cmpint (nc_halo_catalog_get_parent_id (members, i, NULL), ==, 0);

    if (is_central)
    {
      n_central++;
      /* The central sits exactly at the host position. */
      g_assert_cmpfloat (ncm_catalog_get (memc, "ra", i, NULL), ==, 100.0);
      g_assert_cmpfloat (ncm_catalog_get (memc, "dec", i, NULL), ==, 20.0);
      g_assert_cmpfloat (ncm_catalog_get (memc, "z", i, NULL), ==, 0.3);
    }
  }

  /* A massive deterministic-central host yields exactly one central. */
  g_assert_cmpuint (n_central, ==, 1);

  nc_halo_catalog_free (members);
  ncm_rng_free (rng);
  nc_halo_catalog_free (host);
}

void
test_nc_halo_catalog_member_generator_missing_column (TestNcHaloCatalogMemberGenerator *test, gconstpointer pdata)
{
  const gchar *col_names[] = {"ra", "dec", "z_true", NULL};
  NcHaloCatalog *bad       = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_HALO, NULL, NULL, 1, (GStrv) col_names, NULL, 0);
  NcmRNG *rng              = ncm_rng_seeded_new (NULL, 1);
  GError *error            = NULL;
  NcHaloCatalog *members;

  members = nc_halo_catalog_member_generator_generate (test->memgen, bad, test->mset, rng, &error);

  g_assert_null (members);
  g_assert_error (error, NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR, NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR_MISSING_COLUMN);
  g_clear_error (&error);

  ncm_rng_free (rng);
  nc_halo_catalog_free (bad);
}

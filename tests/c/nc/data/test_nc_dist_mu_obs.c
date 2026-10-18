/***************************************************************************
 *            test_nc_dist_mu_obs.c
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_dist_mu_obs.c
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
#include <string.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcDistMuObs
{
  NcDistMuObs *obs;
} TestNcDistMuObs;

void test_nc_dist_mu_obs_new (TestNcDistMuObs *test, gconstpointer pdata);
void test_nc_dist_mu_obs_free (TestNcDistMuObs *test, gconstpointer pdata);

void test_nc_dist_mu_obs_columns (TestNcDistMuObs *test, gconstpointer pdata);
void test_nc_dist_mu_obs_direction (TestNcDistMuObs *test, gconstpointer pdata);
void test_nc_dist_mu_obs_serialize (TestNcDistMuObs *test, gconstpointer pdata);
void test_nc_dist_mu_obs_cf4_snia (TestNcDistMuObs *test, gconstpointer pdata);
void test_nc_dist_mu_obs_missing_column (TestNcDistMuObs *test, gconstpointer pdata);
void test_nc_dist_mu_obs_pitfalls (TestNcDistMuObs *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/dist_mu_obs/columns", TestNcDistMuObs, NULL,
              &test_nc_dist_mu_obs_new,
              &test_nc_dist_mu_obs_columns,
              &test_nc_dist_mu_obs_free);

  g_test_add ("/nc/dist_mu_obs/direction", TestNcDistMuObs, NULL,
              &test_nc_dist_mu_obs_new,
              &test_nc_dist_mu_obs_direction,
              &test_nc_dist_mu_obs_free);

  g_test_add ("/nc/dist_mu_obs/serialize", TestNcDistMuObs, NULL,
              &test_nc_dist_mu_obs_new,
              &test_nc_dist_mu_obs_serialize,
              &test_nc_dist_mu_obs_free);

  g_test_add ("/nc/dist_mu_obs/cf4_snia", TestNcDistMuObs, NULL,
              &test_nc_dist_mu_obs_new,
              &test_nc_dist_mu_obs_cf4_snia,
              &test_nc_dist_mu_obs_free);

  g_test_add ("/nc/dist_mu_obs/pitfalls", TestNcDistMuObs, NULL,
              &test_nc_dist_mu_obs_new,
              &test_nc_dist_mu_obs_pitfalls,
              &test_nc_dist_mu_obs_free);

  g_test_add ("/nc/dist_mu_obs/missing_column/subprocess", TestNcDistMuObs, NULL,
              &test_nc_dist_mu_obs_new,
              &test_nc_dist_mu_obs_missing_column,
              &test_nc_dist_mu_obs_free);

  g_test_run ();
}

void
test_nc_dist_mu_obs_new (TestNcDistMuObs *test, gconstpointer pdata)
{
  const gchar *cols[] = {"pgc", NC_DIST_MU_OBS_COL_RA, NC_DIST_MU_OBS_COL_DEC, NC_DIST_MU_OBS_COL_Z, NC_DIST_MU_OBS_COL_MU, NC_DIST_MU_OBS_COL_SIGMA_MU, "DM_av", NULL};
  NcmCatalog *catalog;
  guint i;

  test->obs = nc_dist_mu_obs_new (4, (GStrv) cols);
  catalog   = NCM_CATALOG (test->obs);

  g_assert_true (NC_IS_DIST_MU_OBS (test->obs));

  for (i = 0; i < 4; i++)
  {
    ncm_catalog_set (catalog, "pgc", i, 100.0 + i, NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_RA, i, 10.0 + 100.0 * i, NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_DEC, i, -60.0 + 40.0 * i, NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_Z, i, 0.01 * (i + 1), NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_MU, i, 33.0 + 0.5 * i, NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_SIGMA_MU, i, 0.1 + 0.01 * i, NULL);
    ncm_catalog_set (catalog, "DM_av", i, 33.1 + 0.5 * i, NULL);
  }
}

void
test_nc_dist_mu_obs_free (TestNcDistMuObs *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_dist_mu_obs_free, test->obs);
}

void
test_nc_dist_mu_obs_columns (TestNcDistMuObs *test, gconstpointer pdata)
{
  NcmCatalog *catalog = NCM_CATALOG (test->obs);
  NcDistMuObs *basic  = nc_dist_mu_obs_new (2, NULL);

  g_assert_cmpuint (ncm_catalog_len (catalog), ==, 4);
  g_assert_cmpuint (ncm_catalog_ncols (catalog), ==, 7);
  g_assert_true (ncm_catalog_has_column (catalog, "DM_av"));
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, "DM_av", 2, NULL), ==, 34.1);
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, 3, NULL), ==, 0.04);

  g_assert_cmpuint (ncm_catalog_len (NCM_CATALOG (basic)), ==, 2);
  g_assert_cmpuint (ncm_catalog_ncols (NCM_CATALOG (basic)), ==, 5);
  g_assert_true (ncm_catalog_has_column (NCM_CATALOG (basic), NC_DIST_MU_OBS_COL_SIGMA_MU));

  nc_dist_mu_obs_free (basic);

  /* Typed columns */
  {
    const gchar *tcols[]              = {"id", NC_DIST_MU_OBS_COL_RA, NC_DIST_MU_OBS_COL_DEC, NC_DIST_MU_OBS_COL_Z, NC_DIST_MU_OBS_COL_MU, NC_DIST_MU_OBS_COL_SIGMA_MU, "flag", NULL};
    const NcmCatalogColType ttypes[7] = {NCM_CATALOG_COL_TYPE_INT, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_BOOL};
    NcDistMuObs *typed                = nc_dist_mu_obs_new_full (3, (GStrv) tcols, ttypes, 7);

    g_assert_cmpint (ncm_catalog_get_col_type (NCM_CATALOG (typed), "id", NULL), ==, NCM_CATALOG_COL_TYPE_INT);
    g_assert_cmpint (ncm_catalog_get_col_type (NCM_CATALOG (typed), "flag", NULL), ==, NCM_CATALOG_COL_TYPE_BOOL);
    g_assert_cmpint (ncm_catalog_get_col_type (NCM_CATALOG (typed), NC_DIST_MU_OBS_COL_MU, NULL), ==, NCM_CATALOG_COL_TYPE_DOUBLE);

    nc_dist_mu_obs_free (typed);
  }
}

void
test_nc_dist_mu_obs_direction (TestNcDistMuObs *test, gconstpointer pdata)
{
  NcmCatalog *catalog = NCM_CATALOG (test->obs);
  NcmTriVec n;
  guint i;

  for (i = 0; i < 4; i++)
  {
    const gdouble ra  = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_RA, i, NULL);
    const gdouble dec = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_DEC, i, NULL);
    const gdouble a   = ncm_c_degree_to_radian (ra);
    const gdouble d   = ncm_c_degree_to_radian (dec);

    nc_dist_mu_obs_get_direction (test->obs, i, &n);

    ncm_assert_cmpdouble_e (n.c[0], ==, cos (d) * cos (a), 1.0e-15, 1.0e-15);
    ncm_assert_cmpdouble_e (n.c[1], ==, cos (d) * sin (a), 1.0e-15, 1.0e-15);
    ncm_assert_cmpdouble_e (n.c[2], ==, sin (d), 1.0e-15, 1.0e-15);
    ncm_assert_cmpdouble_e (ncm_trivec_norm (&n), ==, 1.0, 1.0e-15, 0.0);
  }

  /* Writing a direction and reading it back */
  {
    NcmTriVec m = {
      {
        0.48, 0.6, 0.64
      }
    };

    nc_dist_mu_obs_set_direction (test->obs, 1, &m);
    nc_dist_mu_obs_get_direction (test->obs, 1, &n);

    ncm_assert_cmpdouble_e (n.c[0], ==, 0.48, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (n.c[1], ==, 0.6, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (n.c[2], ==, 0.64, 1.0e-14, 0.0);
  }
}

void
test_nc_dist_mu_obs_serialize (TestNcDistMuObs *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcDistMuObs *dup  = NC_DIST_MU_OBS (ncm_serialize_dup_obj (ser, G_OBJECT (test->obs)));
  NcmCatalog *a     = NCM_CATALOG (test->obs);
  NcmCatalog *b     = NCM_CATALOG (dup);
  guint i;

  g_assert_true (NC_IS_DIST_MU_OBS (dup));
  g_assert_cmpuint (ncm_catalog_len (b), ==, ncm_catalog_len (a));
  g_assert_cmpuint (ncm_catalog_ncols (b), ==, ncm_catalog_ncols (a));

  for (i = 0; i < ncm_catalog_len (a); i++)
  {
    ncm_assert_cmpdouble (ncm_catalog_get (b, NC_DIST_MU_OBS_COL_MU, i, NULL), ==, ncm_catalog_get (a, NC_DIST_MU_OBS_COL_MU, i, NULL));
    ncm_assert_cmpdouble (ncm_catalog_get (b, "DM_av", i, NULL), ==, ncm_catalog_get (a, "DM_av", i, NULL));
    ncm_assert_cmpdouble (ncm_catalog_get (b, NC_DIST_MU_OBS_COL_DEC, i, NULL), ==, ncm_catalog_get (a, NC_DIST_MU_OBS_COL_DEC, i, NULL));
  }

  nc_dist_mu_obs_free (dup);
  ncm_serialize_free (ser);
}

void
test_nc_dist_mu_obs_cf4_snia (TestNcDistMuObs *test, gconstpointer pdata)
{
  NcDistMuObs *cf4    = nc_dist_mu_obs_new_from_catalog_id (NC_DIST_MU_OBS_CATALOG_CF4_SNIA);
  NcmCatalog *catalog = NCM_CATALOG (cf4);
  gchar *source       = NULL;
  NcmTriVec n;
  guint i;

  g_assert_true (NC_IS_DIST_MU_OBS (cf4));
  g_assert_cmpuint (ncm_catalog_len (catalog), ==, 1008);
  g_assert_true (ncm_catalog_has_column (catalog, "DM_av"));
  g_assert_true (ncm_catalog_has_column (catalog, "Vcmb"));
  g_assert_cmpint (ncm_catalog_get_col_type (catalog, "pgc", NULL), ==, NCM_CATALOG_COL_TYPE_INT);

  g_assert_true (ncm_var_dict_get_string (ncm_catalog_peek_meta (catalog), "source", &source));
  g_assert_nonnull (strstr (source, "CF4 SNIa"));
  g_free (source);

  /* First row of the EDD export: pgc 250, RA 0.8958, DE 23.2008, Vcmb 6914, DM_zp 34.912, eDM 0.13, DM_av 35.05 */
  g_assert_cmpint (ncm_catalog_get_int (catalog, "pgc", 0, NULL), ==, 250);
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_RA, 0, NULL), ==, 0.8958);
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_DEC, 0, NULL), ==, 23.2008);
  ncm_assert_cmpdouble_e (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, 0, NULL), ==, 6914.0 / 299792.458, 1.0e-15, 0.0);
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_MU, 0, NULL), ==, 34.912);
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_SIGMA_MU, 0, NULL), ==, 0.13);
  ncm_assert_cmpdouble (ncm_catalog_get (catalog, "DM_av", 0, NULL), ==, 35.05);

  for (i = 0; i < ncm_catalog_len (catalog); i++)
  {
    nc_dist_mu_obs_get_direction (cf4, i, &n);
    ncm_assert_cmpdouble_e (ncm_trivec_norm (&n), ==, 1.0, 1.0e-15, 0.0);
    g_assert_cmpfloat (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_SIGMA_MU, i, NULL), >, 0.0);
    g_assert_cmpfloat (ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, i, NULL), >, 0.0);
  }

  nc_dist_mu_obs_free (cf4);
}

void
test_nc_dist_mu_obs_pitfalls (TestNcDistMuObs *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/nc/dist_mu_obs/missing_column/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*column `sigma_mu' is required*");
}

void
test_nc_dist_mu_obs_missing_column (TestNcDistMuObs *test, gconstpointer pdata)
{
  const gchar *cols[] = {NC_DIST_MU_OBS_COL_RA, NC_DIST_MU_OBS_COL_DEC, NC_DIST_MU_OBS_COL_Z, NC_DIST_MU_OBS_COL_MU, NULL};
  NcDistMuObs *obs    = nc_dist_mu_obs_new (2, (GStrv) cols);

  nc_dist_mu_obs_free (obs);
}


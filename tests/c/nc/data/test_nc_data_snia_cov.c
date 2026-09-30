/***************************************************************************
 *            test_nc_data_snia_cov.c
 *
 *  Tue Sep 29 10:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

/* Catalogs of TEST_NSNIA supernovae built in memory; no catalog file is read. */
#define TEST_NSNIA 3

void test_nc_data_snia_cov_lazy_cov_full (gconstpointer pdata);
void test_nc_data_snia_cov_v2_ignores_cov_full (void);
void test_nc_data_snia_cov_v2_no_lightcurve (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();


  g_test_add_data_func ("/nc/data_snia_cov/v0/lazy_cov_full", GUINT_TO_POINTER (0), &test_nc_data_snia_cov_lazy_cov_full);
  g_test_add_data_func ("/nc/data_snia_cov/v1/lazy_cov_full", GUINT_TO_POINTER (1), &test_nc_data_snia_cov_lazy_cov_full);
  g_test_add_func ("/nc/data_snia_cov/v2/ignores_cov_full", &test_nc_data_snia_cov_v2_ignores_cov_full);
  g_test_add_func ("/nc/data_snia_cov/v2/no_lightcurve", &test_nc_data_snia_cov_v2_no_lightcurve);

  g_test_run ();
}

/* A symmetric, diagonally dominant (so positive definite) 3n x 3n light-curve covariance. */
static NcmMatrix *
_test_nc_data_snia_cov_lightcurve_cov (void)
{
  const guint tn = 3 * TEST_NSNIA;
  NcmMatrix *cov = ncm_matrix_new (tn, tn);
  guint i, j;

  for (i = 0; i < tn; i++)
    for (j = 0; j < tn; j++)
      ncm_matrix_set (cov, i, j, (i == j) ? 2.0 : 0.1 / (1.0 + abs ((gint) i - (gint) j)));

  return cov;
}

static void
_test_nc_data_snia_cov_assert_equal (NcmMatrix *a, NcmMatrix *b)
{
  guint i, j;

  g_assert_cmpuint (ncm_matrix_nrows (a), ==, ncm_matrix_nrows (b));
  g_assert_cmpuint (ncm_matrix_ncols (a), ==, ncm_matrix_ncols (b));

  for (i = 0; i < ncm_matrix_nrows (a); i++)
    for (j = 0; j < ncm_matrix_ncols (a); j++)
      g_assert_cmpfloat (ncm_matrix_get (a, i, j), ==, ncm_matrix_get (b, i, j));
}

/* Versions 0 and 1 allocate the light-curve covariance on the first write and free it on resize */
void
test_nc_data_snia_cov_lazy_cov_full (gconstpointer pdata)
{
  const guint cat_version = GPOINTER_TO_UINT (pdata);
  const guint n           = TEST_NSNIA;
  NcDataSNIACov *snia_cov = nc_data_snia_cov_new (FALSE, cat_version);
  NcmMatrix *cov          = _test_nc_data_snia_cov_lightcurve_cov ();
  NcmVector *packed;
  guint i, j, ij;

  ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov), n);
  g_assert_null (nc_data_snia_cov_peek_cov_full (snia_cov));
  g_assert_null (nc_data_snia_cov_peek_cov_packed (snia_cov));

  nc_data_snia_cov_set_cov_full (snia_cov, cov);
  _test_nc_data_snia_cov_assert_equal (nc_data_snia_cov_peek_cov_full (snia_cov), cov);

  packed = nc_data_snia_cov_peek_cov_packed (snia_cov);
  g_assert_cmpuint (ncm_vector_len (packed), ==, NC_DATA_SNIA_COV_ORDER_LENGTH * n * (n + 1) / 2);

  /* Each packed entry averages (i, j) and (j, i) of its block. */
  ij = 0;

  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      const guint bi[NC_DATA_SNIA_COV_ORDER_LENGTH] = {0, 0, 0, 1, 1, 2};
      const guint bj[NC_DATA_SNIA_COV_ORDER_LENGTH] = {0, 1, 2, 1, 2, 2};
      guint k;

      for (k = 0; k < NC_DATA_SNIA_COV_ORDER_LENGTH; k++)
      {
        const gdouble c_ij = ncm_matrix_get (cov, bi[k] * n + i, bj[k] * n + j);
        const gdouble c_ji = ncm_matrix_get (cov, bi[k] * n + j, bj[k] * n + i);

        g_assert_cmpfloat (ncm_vector_get (packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + k), ==, 0.5 * (c_ij + c_ji));
      }

      ij++;
    }
  }

  /* With a complete covariance the positive-definiteness check leaves the matrix unchanged. */
  g_object_set (snia_cov, "has-complete-cov", TRUE, NULL);
  nc_data_snia_cov_set_cov_full (snia_cov, cov);
  _test_nc_data_snia_cov_assert_equal (nc_data_snia_cov_peek_cov_full (snia_cov), cov);

  ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov), n + 1);
  g_assert_null (nc_data_snia_cov_peek_cov_full (snia_cov));
  g_assert_null (nc_data_snia_cov_peek_cov_packed (snia_cov));

  ncm_matrix_free (cov);
  ncm_data_free (NCM_DATA (snia_cov));
}

/* Version 2 has no light-curve covariance: a cov-full set, as in older serializations, is dropped */
void
test_nc_data_snia_cov_v2_ignores_cov_full (void)
{
  NcDataSNIACov *snia_cov = nc_data_snia_cov_new (FALSE, 2);
  NcmMatrix *cov          = _test_nc_data_snia_cov_lightcurve_cov ();
  NcmSerialize *ser       = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GVariant *var, *params, *cov_full_var;

  ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov), TEST_NSNIA);
  g_object_set (snia_cov, "cov-full", cov, NULL);

  g_assert_null (nc_data_snia_cov_peek_cov_full (snia_cov));
  g_assert_null (nc_data_snia_cov_peek_cov_packed (snia_cov));

  var          = ncm_serialize_to_variant (ser, G_OBJECT (snia_cov));
  params       = g_variant_get_child_value (var, 1);
  cov_full_var = g_variant_lookup_value (params, "cov-full", NULL);
  g_assert_null (cov_full_var);

  g_variant_unref (params);
  g_variant_unref (var);
  ncm_serialize_free (ser);
  ncm_matrix_free (cov);
  ncm_data_free (NCM_DATA (snia_cov));
}

/* The light-curve paths of a version 2 catalog abort instead of using a matrix that does not exist */
void
test_nc_data_snia_cov_v2_no_lightcurve (void)
{
  if (g_test_subprocess ())
  {
    NcDataSNIACov *snia_cov = nc_data_snia_cov_new (FALSE, 2);
    NcDistance *dist        = nc_distance_new (2.0);
    NcHICosmo *cosmo        = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
    NcSNIADistCov *dcov     = nc_snia_dist_cov_new (dist, 1);
    NcmMSet *mset           = ncm_mset_new (cosmo, NULL, dcov, NULL);

    ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov), TEST_NSNIA);
    g_object_set (snia_cov, "has-complete-cov", TRUE, NULL);
    nc_data_snia_cov_estimate_width_colour (snia_cov, mset);

    return;
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*catalog version 2 has no light-curve covariance*");
}


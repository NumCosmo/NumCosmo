/***************************************************************************
 *            test_nc_data_hubble_cov.c
 *
 *  Wed Oct 21 13:26:04 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_data_hubble_cov.c
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

#define _TEST_LOUBSER2025_NP 3

typedef struct _TestNcDataHubbleCov
{
  NcDataHubbleCov *hubble_cov;
} TestNcDataHubbleCov;

/*
 * Loubser (2025), arXiv:2511.02730, table 1: z, H(z), statistical and systematic errors
 * [km/s/Mpc]; section 3.3: total-error correlations of the pairs (0, 1), (0, 2), (1, 2).
 */
static const gdouble _loubser2025_z[_TEST_LOUBSER2025_NP]          = {0.46, 0.67, 0.83};
static const gdouble _loubser2025_H[_TEST_LOUBSER2025_NP]          = {88.48, 119.45, 108.28};
static const gdouble _loubser2025_sigma_stat[_TEST_LOUBSER2025_NP] = {0.57, 6.39, 10.07};
static const gdouble _loubser2025_sigma_sys[_TEST_LOUBSER2025_NP]  = {12.32, 16.64, 15.08};
static const gdouble _loubser2025_rho[_TEST_LOUBSER2025_NP]        = {0.932, 0.830, 0.776};

void test_nc_data_hubble_cov_new_loubser2025 (TestNcDataHubbleCov *test, gconstpointer pdata);
void test_nc_data_hubble_cov_free (TestNcDataHubbleCov *test, gconstpointer pdata);

void test_nc_data_hubble_cov_loubser2025_values (TestNcDataHubbleCov *test, gconstpointer pdata);
void test_nc_data_hubble_cov_loubser2025_m2lnL (TestNcDataHubbleCov *test, gconstpointer pdata);
void test_nc_data_hubble_cov_serialization (TestNcDataHubbleCov *test, gconstpointer pdata);

void test_nc_data_hubble_cov_api_set_size (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_hubble_cov/loubser2025/values", TestNcDataHubbleCov, NULL,
              &test_nc_data_hubble_cov_new_loubser2025,
              &test_nc_data_hubble_cov_loubser2025_values,
              &test_nc_data_hubble_cov_free);
  g_test_add ("/nc/data_hubble_cov/loubser2025/m2lnL", TestNcDataHubbleCov, NULL,
              &test_nc_data_hubble_cov_new_loubser2025,
              &test_nc_data_hubble_cov_loubser2025_m2lnL,
              &test_nc_data_hubble_cov_free);
  g_test_add ("/nc/data_hubble_cov/loubser2025/serialization", TestNcDataHubbleCov, NULL,
              &test_nc_data_hubble_cov_new_loubser2025,
              &test_nc_data_hubble_cov_serialization,
              &test_nc_data_hubble_cov_free);

  g_test_add_func ("/nc/data_hubble_cov/api/set_size", &test_nc_data_hubble_cov_api_set_size);

  g_test_run ();
}

void
test_nc_data_hubble_cov_new_loubser2025 (TestNcDataHubbleCov *test, gconstpointer pdata)
{
  test->hubble_cov = nc_data_hubble_cov_new_from_id (NC_DATA_HUBBLE_COV_LOUBSER2025);
  g_assert_true (NC_IS_DATA_HUBBLE_COV (test->hubble_cov));
}

void
test_nc_data_hubble_cov_free (TestNcDataHubbleCov *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (test->hubble_cov));
}

void
test_nc_data_hubble_cov_loubser2025_values (TestNcDataHubbleCov *test, gconstpointer pdata)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->hubble_cov);
  NcmVector *mean        = ncm_data_gauss_cov_peek_mean (gauss);
  NcmMatrix *cov         = ncm_data_gauss_cov_peek_cov (gauss);
  NcmVector *z           = NULL;
  guint i, j, k = 0;

  g_object_get (test->hubble_cov, "z", &z, NULL);

  g_assert_cmpuint (ncm_data_gauss_cov_get_size (gauss), ==, _TEST_LOUBSER2025_NP);
  g_assert_cmpuint (ncm_vector_len (z), ==, _TEST_LOUBSER2025_NP);

  for (i = 0; i < _TEST_LOUBSER2025_NP; i++)
  {
    ncm_assert_cmpdouble (ncm_vector_get (z, i), ==, _loubser2025_z[i]);
    ncm_assert_cmpdouble (ncm_vector_get (mean, i), ==, _loubser2025_H[i]);

    for (j = 0; j < _TEST_LOUBSER2025_NP; j++)
    {
      /* Statistical errors independent, systematic error fully correlated. */
      const gdouble stat2  = (i == j) ? gsl_pow_2 (_loubser2025_sigma_stat[i]) : 0.0;
      const gdouble cov_ij = stat2 + _loubser2025_sigma_sys[i] * _loubser2025_sigma_sys[j];

      ncm_assert_cmpdouble_e (ncm_matrix_get (cov, i, j), ==, cov_ij, 1.0e-15, 0.0);
    }
  }

  /* The paper quotes the correlations to three digits. */
  for (i = 0; i < _TEST_LOUBSER2025_NP; i++)
  {
    for (j = i + 1; j < _TEST_LOUBSER2025_NP; j++)
    {
      const gdouble rho = ncm_matrix_get (cov, i, j) / sqrt (ncm_matrix_get (cov, i, i) * ncm_matrix_get (cov, j, j));

      g_assert_cmpfloat (fabs (rho - _loubser2025_rho[k]), <, 1.0e-3);
      k++;
    }
  }

  ncm_vector_free (z);
}

void
test_nc_data_hubble_cov_loubser2025_m2lnL (TestNcDataHubbleCov *test, gconstpointer pdata)
{
  NcmData *data          = NCM_DATA (test->hubble_cov);
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (test->hubble_cov);
  NcmMatrix *cov         = ncm_data_gauss_cov_peek_cov (gauss);
  NcHICosmo *cosmo       = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcmMSet *mset          = ncm_mset_new (NCM_MODEL (cosmo), NULL, NULL);
  gdouble r[_TEST_LOUBSER2025_NP], inv[_TEST_LOUBSER2025_NP][_TEST_LOUBSER2025_NP];
  gdouble chi2 = 0.0, m2lnL, det;
  guint i, j;

  /* Cofactor inverse of the 3x3 covariance, independent of the library's Cholesky path. */
  {
    const gdouble a = ncm_matrix_get (cov, 0, 0), b = ncm_matrix_get (cov, 0, 1), c = ncm_matrix_get (cov, 0, 2);
    const gdouble d = ncm_matrix_get (cov, 1, 1), e = ncm_matrix_get (cov, 1, 2), f = ncm_matrix_get (cov, 2, 2);

    inv[0][0] = d * f - e * e;
    inv[0][1] = c * e - b * f;
    inv[0][2] = b * e - c * d;
    inv[1][1] = a * f - c * c;
    inv[1][2] = b * c - a * e;
    inv[2][2] = a * d - b * b;
    inv[1][0] = inv[0][1];
    inv[2][0] = inv[0][2];
    inv[2][1] = inv[1][2];
    det       = a * inv[0][0] + b * inv[0][1] + c * inv[0][2];
  }

  for (i = 0; i < _TEST_LOUBSER2025_NP; i++)
    r[i] = _loubser2025_H[i] - nc_hicosmo_H (cosmo, _loubser2025_z[i]);

  for (i = 0; i < _TEST_LOUBSER2025_NP; i++)
    for (j = 0; j < _TEST_LOUBSER2025_NP; j++)
      chi2 += r[i] * inv[i][j] * r[j] / det;

  ncm_data_prepare (data, mset);
  ncm_data_m2lnL_val (data, mset, &m2lnL);

  ncm_assert_cmpdouble_e (m2lnL, ==, chi2, 1.0e-12, 0.0);

  ncm_mset_free (mset);
  nc_hicosmo_free (cosmo);
}

void
test_nc_data_hubble_cov_serialization (TestNcDataHubbleCov *test, gconstpointer pdata)
{
  NcDataHubbleCov *dup    = NC_DATA_HUBBLE_COV (ncm_serialize_global_dup_obj (G_OBJECT (test->hubble_cov)));
  NcmDataGaussCov *gauss  = NCM_DATA_GAUSS_COV (test->hubble_cov);
  NcmDataGaussCov *dgauss = NCM_DATA_GAUSS_COV (dup);
  NcmVector *z = NULL, *dz = NULL;
  guint i;

  g_assert_true (dup != test->hubble_cov);

  g_object_get (test->hubble_cov, "z", &z, NULL);
  g_object_get (dup, "z", &dz, NULL);

  g_assert_cmpuint (ncm_vector_len (dz), ==, ncm_vector_len (z));

  for (i = 0; i < ncm_vector_len (z); i++)
  {
    ncm_assert_cmpdouble (ncm_vector_get (dz, i), ==, ncm_vector_get (z, i));
    ncm_assert_cmpdouble (ncm_vector_get (ncm_data_gauss_cov_peek_mean (dgauss), i), ==,
                          ncm_vector_get (ncm_data_gauss_cov_peek_mean (gauss), i));
  }

  g_assert_cmpfloat (ncm_matrix_cmp (ncm_data_gauss_cov_peek_cov (gauss), ncm_data_gauss_cov_peek_cov (dgauss), 0.0), ==, 0.0);

  ncm_vector_free (z);
  ncm_vector_free (dz);
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (dup));
}

void
test_nc_data_hubble_cov_api_set_size (void)
{
  NcDataHubbleCov *hubble_cov = nc_data_hubble_cov_new_empty ();
  NcmDataGaussCov *gauss      = NCM_DATA_GAUSS_COV (hubble_cov);
  NcmVector *z                = NULL;

  g_assert_cmpuint (ncm_data_gauss_cov_get_size (gauss), ==, 0);

  ncm_data_gauss_cov_set_size (gauss, 5);
  g_object_get (hubble_cov, "z", &z, NULL);
  g_assert_cmpuint (ncm_data_gauss_cov_get_size (gauss), ==, 5);
  g_assert_cmpuint (ncm_vector_len (z), ==, 5);
  ncm_vector_free (z);

  ncm_data_gauss_cov_set_size (gauss, 0);
  g_object_get (hubble_cov, "z", &z, NULL);
  g_assert_null (z);

  NCM_TEST_FREE (ncm_data_free, NCM_DATA (hubble_cov));
}


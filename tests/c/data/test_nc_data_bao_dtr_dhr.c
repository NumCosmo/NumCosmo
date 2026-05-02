/***************************************************************************
 *            test_nc_data_bao_dtr_dhr.c
 *
 *  Thu Apr 17 19:50:16 2025
 *  Copyright  2025  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * test_nc_data_bao_dtr_dhr.c
 * Copyright (C) 2025 Mariana Penna-Lima <pennalima@unb.br>
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

typedef struct _TestNcDataBaoDtrDHr
{
  NcDataBaoDtrDHr *dtdh;
  NcDataBaoId id;
} TestNcDataBaoDtrDHr;

void test_nc_data_bao_dtr_dhr_free (TestNcDataBaoDtrDHr *test, gconstpointer pdata);

void test_nc_data_bao_dtr_dhr_new_dr12_2016_dr16compatible (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_set_sample_dr12_2016_dr16compatible (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_new_dr16_lrg_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_set_sample_dr16_lrg_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_new_dr16_qso_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_set_sample_dr16_qso_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_new_desi_dr1_lym_2025 (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_set_sample_desi_dr1_lym_2025 (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_get_distance (TestNcDataBaoDtrDHr *test, gconstpointer pdata);
void test_nc_data_bao_dtr_dhr_get_redshift (TestNcDataBaoDtrDHr *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_bao_dtr_dhr/set_sample/dr12dr16compatible", TestNcDataBaoDtrDHr, NULL,
              &test_nc_data_bao_dtr_dhr_new_dr12_2016_dr16compatible,
              &test_nc_data_bao_dtr_dhr_set_sample_dr12_2016_dr16compatible,
              &test_nc_data_bao_dtr_dhr_free);
  g_test_add ("/nc/data_bao_dtr_dhr/set_sample/dr16lrg2021", TestNcDataBaoDtrDHr, NULL,
              &test_nc_data_bao_dtr_dhr_new_dr16_lrg_2021,
              &test_nc_data_bao_dtr_dhr_set_sample_dr16_lrg_2021,
              &test_nc_data_bao_dtr_dhr_free);
  g_test_add ("/nc/data_bao_dtr_dhr/set_sample/dr16qso2021", TestNcDataBaoDtrDHr, NULL,
              &test_nc_data_bao_dtr_dhr_new_dr16_qso_2021,
              &test_nc_data_bao_dtr_dhr_set_sample_dr16_qso_2021,
              &test_nc_data_bao_dtr_dhr_free);
  g_test_add ("/nc/data_bao_dtr_dhr/set_sample/desi2025", TestNcDataBaoDtrDHr, NULL,
              &test_nc_data_bao_dtr_dhr_new_desi_dr1_lym_2025,
              &test_nc_data_bao_dtr_dhr_set_sample_desi_dr1_lym_2025,
              &test_nc_data_bao_dtr_dhr_free);
  g_test_add ("/nc/data_bao_dtr_dhr/get_distance", TestNcDataBaoDtrDHr, NULL,
              &test_nc_data_bao_dtr_dhr_new_desi_dr1_lym_2025,
              &test_nc_data_bao_dtr_dhr_get_distance,
              &test_nc_data_bao_dtr_dhr_free);
  g_test_add ("/nc/data_bao_dtr_dhr/get_redshift", TestNcDataBaoDtrDHr, NULL,
              &test_nc_data_bao_dtr_dhr_new_desi_dr1_lym_2025,
              &test_nc_data_bao_dtr_dhr_get_redshift,
              &test_nc_data_bao_dtr_dhr_free);
  g_test_run ();
}

void
test_nc_data_bao_dtr_dhr_free (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcmData *dtdh = NCM_DATA (test->dtdh);

  NCM_TEST_FREE (ncm_data_free, dtdh);
}

/* SDSS DR12 2016 - DR16 compatible */

void
test_nc_data_bao_dtr_dhr_new_dr12_2016_dr16compatible (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE;
  NcDistance *dist = nc_distance_new (3.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_dtr_dhr_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->dtdh = NC_DATA_BAO_DTR_DHR (data);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dtr_dhr_set_sample_dr12_2016_dr16compatible (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoDtrDHr *dtdh      = test->dtdh;
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (dtdh);
  NcDataBaoId id             = NC_DATA_BAO_DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE;
  NcmVector *y               = ncm_data_gauss_cov_peek_mean (gauss_cov);
  NcmMatrix *cov_m           = ncm_data_gauss_cov_peek_cov (gauss_cov);
  NcmVector *x               = NULL;

  const gdouble z0 = 0.38;
  const gdouble z1 = 0.38;
  const gdouble z2 = 0.51;
  const gdouble z3 = 0.51;

  const gdouble bf0 = 1.023406e+01;
  const gdouble bf1 = 2.498058e+01;
  const gdouble bf2 = 1.336595e+01;
  const gdouble bf3 = 2.231656e+01;

  gdouble covar[16] = {
    2.860520e-02, -4.939281e-02,  1.489688e-02, -1.387079e-02,
    -4.939281e-02,  5.307187e-01, -2.423513e-02,  1.767087e-01,
    1.489688e-02, -2.423513e-02,  4.147534e-02, -4.873962e-02,
    -1.387079e-02,  1.767087e-01, -4.873962e-02, 3.268589e-01
  };

  NcmMatrix *cov = ncm_matrix_new_data_static (covar, 4, 4);

  g_assert_true (dtdh != NULL);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (dtdh));

  g_assert_cmpuint (test->id, ==, id);

  g_object_get (dtdh, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (x, 3), ==, z3);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 3), ==, bf3);

  ncm_matrix_cmp (cov_m, cov, 1e-10);

  ncm_matrix_free (cov);
  ncm_vector_free (x);
}

/* SDSS DR16 LRG 2021*/

void
test_nc_data_bao_dtr_dhr_new_dr16_lrg_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021;
  NcDistance *dist = nc_distance_new (3.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_dtr_dhr_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->dtdh = NC_DATA_BAO_DTR_DHR (data);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dtr_dhr_set_sample_dr16_lrg_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoDtrDHr *dtdh      = test->dtdh;
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (dtdh);
  NcDataBaoId id             = NC_DATA_BAO_DTR_DHR_SDSS_DR16_LRG_2021;
  NcmVector *y               = ncm_data_gauss_cov_peek_mean (gauss_cov);
  NcmMatrix *cov_m           = ncm_data_gauss_cov_peek_cov (gauss_cov);
  NcmVector *x               = NULL;

  const gdouble z0 = 0.698;
  const gdouble z1 = 0.698;

  const gdouble bf0 = 17.85823691865007;
  const gdouble bf1 = 19.32575373059217;

  gdouble covar[4] = {
    0.1076634008565565, -0.05831820341302727,
    -0.0583182034130273, 0.2838176386340292
  };

  NcmMatrix *cov = ncm_matrix_new_data_static (covar, 2, 2);

  g_assert_true (dtdh != NULL);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (dtdh));

  g_assert_cmpuint (test->id, ==, id);

  g_object_get (dtdh, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);

  ncm_matrix_cmp (cov_m, cov, 1e-10);

  ncm_matrix_free (cov);
  ncm_vector_free (x);
}

/* SDSS DR16 QSO 2021*/

void
test_nc_data_bao_dtr_dhr_new_dr16_qso_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021;
  NcDistance *dist = nc_distance_new (2.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_dtr_dhr_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->dtdh = NC_DATA_BAO_DTR_DHR (data);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dtr_dhr_set_sample_dr16_qso_2021 (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoDtrDHr *dtdh      = test->dtdh;
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (dtdh);
  NcDataBaoId id             = NC_DATA_BAO_DTR_DHR_SDSS_DR16_QSO_2021;
  NcmVector *y               = ncm_data_gauss_cov_peek_mean (gauss_cov);
  NcmMatrix *cov_m           = ncm_data_gauss_cov_peek_cov (gauss_cov);
  NcmVector *x               = NULL;

  const gdouble z0 = 1.48;
  const gdouble z1 = 1.48;

  const gdouble bf0 = 30.6876;
  const gdouble bf1 = 13.2609;

  gdouble covar[4] = {
    0.63731604, 0.1706891,
    0.1706891,  0.30468415
  };

  NcmMatrix *cov = ncm_matrix_new_data_static (covar, 2, 2);

  g_assert_true (dtdh != NULL);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (dtdh));

  g_assert_cmpuint (test->id, ==, id);

  g_object_get (dtdh, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);

  ncm_matrix_cmp (cov_m, cov, 1e-10);

  ncm_matrix_free (cov);
  ncm_vector_free (x);
}

/* DESI DR1 - Lyman alpha 2025 */

void
test_nc_data_bao_dtr_dhr_new_desi_dr1_lym_2025 (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_DTR_DHR_DESI_DR1_LYA_2025;
  NcDistance *dist = nc_distance_new (3.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_dtr_dhr_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->dtdh = NC_DATA_BAO_DTR_DHR (data);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dtr_dhr_set_sample_desi_dr1_lym_2025 (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDataBaoDtrDHr *dtdh      = test->dtdh;
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (dtdh);
  NcDataBaoId id             = NC_DATA_BAO_DTR_DHR_DESI_DR1_LYA_2025;
  NcmVector *y               = ncm_data_gauss_cov_peek_mean (gauss_cov);
  NcmMatrix *cov_m           = ncm_data_gauss_cov_peek_cov (gauss_cov);
  NcmVector *x               = NULL;

  const gdouble z0 = 2.33;
  const gdouble z1 = 2.33;

  const gdouble bf0 = 39.71;
  const gdouble bf1 = 8.52;

  gdouble covar[4] = {
    0.9025, -0.07752,
    -0.07752, 0.0289
  };

  NcmMatrix *cov = ncm_matrix_new_data_static (covar, 2, 2);

  g_assert_true (dtdh != NULL);
  g_assert_true (NC_IS_DATA_BAO_DTR_DHR (dtdh));

  g_assert_cmpuint (test->id, ==, id);

  g_object_get (dtdh, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);

  ncm_matrix_cmp (cov_m, cov, 1e-10);

  ncm_matrix_free (cov);
  ncm_vector_free (x);
}

void
test_nc_data_bao_dtr_dhr_get_distance (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcDistance *dist      = NULL;
  NcDataBaoDtrDHr *dtdh = test->dtdh;

  g_object_get (dtdh, "dist", &dist, NULL);
  g_assert_true (dist != NULL);
  g_assert (NC_IS_DISTANCE (dist));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dtr_dhr_get_redshift (TestNcDataBaoDtrDHr *test, gconstpointer pdata)
{
  NcmVector *x          = NULL;
  NcDataBaoDtrDHr *dtdh = test->dtdh;

  g_object_get (dtdh, "z", &x, NULL);
  g_assert_true (x != NULL);
  g_assert (NCM_IS_VECTOR (x));

  ncm_vector_free (x);
}


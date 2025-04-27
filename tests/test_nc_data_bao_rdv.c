/***************************************************************************
 *            test_nc_data_bao_rdv.c
 *
 *  Wed February 11 12:11:42 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * test_nc_data_bao_rdv.c
 * Copyright (C) 2015 Mariana Penna-Lima <pennalima@gmail.com>
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

typedef struct _TestNcDataBaoRDV
{
  NcDataBaoRDV *rdv;
  NcDataBaoId id;
} TestNcDataBaoRDV;

void test_nc_data_bao_rdv_free (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_percival2007 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_percival2007 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_percival2010 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_percival2010 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_beutler2011 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_beutler2011 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_padmanabhan2012 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_padmanabhan2012 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_anderson2012 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_anderson2012 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_blake2012 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_blake2012 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_kazin2014 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_kazin2014 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_desi_dr1_bgs_qso_2024 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_desi_dr1_bgs_qso_2024 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_new_desi_dr2_bgs_2025 (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_set_sample_desi_dr2_bgs_2025 (TestNcDataBaoRDV *test, gconstpointer pdata);

void test_nc_data_bao_rdv_get_distance (TestNcDataBaoRDV *test, gconstpointer pdata);
void test_nc_data_bao_rdv_get_redshift (TestNcDataBaoRDV *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_bao_rdv/set_sample/percival2007", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_percival2007,
              &test_nc_data_bao_rdv_set_sample_percival2007,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/percival2010", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_percival2010,
              &test_nc_data_bao_rdv_set_sample_percival2010,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/beutler2011", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_beutler2011,
              &test_nc_data_bao_rdv_set_sample_beutler2011,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/padmanabhan2012", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_padmanabhan2012,
              &test_nc_data_bao_rdv_set_sample_padmanabhan2012,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/anderson2012", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_anderson2012,
              &test_nc_data_bao_rdv_set_sample_anderson2012,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/blake2012", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_blake2012,
              &test_nc_data_bao_rdv_set_sample_blake2012,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/kazin2014", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_kazin2014,
              &test_nc_data_bao_rdv_set_sample_kazin2014,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/desi_dr1_bgs_qso_2024", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_desi_dr1_bgs_qso_2024,
              &test_nc_data_bao_rdv_set_sample_desi_dr1_bgs_qso_2024,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/set_sample/desi_dr2_bgs_2025", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_desi_dr2_bgs_2025,
              &test_nc_data_bao_rdv_set_sample_desi_dr2_bgs_2025,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/get_distance", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_desi_dr1_bgs_qso_2024,
              &test_nc_data_bao_rdv_get_distance,
              &test_nc_data_bao_rdv_free);
  g_test_add ("/nc/data_bao_rdv/get_redshift", TestNcDataBaoRDV, NULL,
              &test_nc_data_bao_rdv_new_desi_dr1_bgs_qso_2024,
              &test_nc_data_bao_rdv_get_redshift,
              &test_nc_data_bao_rdv_free);


  g_test_run ();
}

void
test_nc_data_bao_rdv_free (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *rdv = NCM_DATA (test->rdv);

  NCM_TEST_FREE (ncm_data_free, rdv);
}

/* Percival2007 */

void
test_nc_data_bao_rdv_new_percival2007 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_PERCIVAL2007;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_percival2007 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_PERCIVAL2007;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = TRUE;
  NcmVector *x        = NULL;

  const gdouble z0 = 0.2;
  const gdouble z1 = 0.35;

  const gdouble bf0 = 0.1980;
  const gdouble bf1 = 0.1094;

  const gdouble icov00 = 35059.0;
  const gdouble icov01 = -24031.0;
  const gdouble icov10 = -24031.0;
  const gdouble icov11 = 108300.0;

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 1), ==, icov11);

  ncm_vector_free (x);
}

/* Percival 2010 */

void
test_nc_data_bao_rdv_new_percival2010 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_RDV_PERCIVAL2010;
  NcDistance *dist = nc_distance_new (2.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_percival2010 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_PERCIVAL2010;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = TRUE;
  NcmVector *x        = NULL;

  const gdouble z0 = 0.2;
  const gdouble z1 = 0.35;

  const gdouble bf0 = 0.1905;
  const gdouble bf1 = 0.1097;

  const gdouble icov00 = 30124.0;
  const gdouble icov01 = -17227.0;
  const gdouble icov10 = -17227.0;
  const gdouble icov11 = 86977.0;

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (R_DV);
  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 1), ==, icov11);

  ncm_vector_free (x);
}

/* Beutler 2011 */

void
test_nc_data_bao_rdv_new_beutler2011 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_BEUTLER2011;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_beutler2011 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_BEUTLER2011;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = TRUE;
  NcmVector *x        = NULL;

  const gdouble z0  = 0.106;
  const gdouble bf0 = 0.336 / 1.027;

  const gdouble icov00 = 1.0 * 1.027 * 1.027 / (0.015 * 0.015);

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);

  ncm_vector_free (x);
}

/* Padmanabhan 2012 */

void
test_nc_data_bao_rdv_new_padmanabhan2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_PADMANABHAN2012;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_padmanabhan2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_PADMANABHAN2012;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = FALSE;
  NcmVector *x        = NULL;

  const gdouble z0  = 0.35;
  const gdouble bf0 = 8.88 * 1.025;

  const gdouble icov00 = 1.0 / (0.17 * 0.17 * 1.025 * 1.025);

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (!R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);

  ncm_vector_free (x);
}

/* Anderson 2012 */

void
test_nc_data_bao_rdv_new_anderson2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_ANDERSON2012;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_anderson2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_ANDERSON2012;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = FALSE;
  NcmVector *x        = NULL;

  const gdouble z0  = 0.57;
  const gdouble bf0 = 13.67;

  const gdouble icov00 = 1.0 / (0.22 * 0.22);

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (!R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);

  ncm_vector_free (x);
}

/* Blake 2012 */

void
test_nc_data_bao_rdv_new_blake2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_BLAKE2012;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_blake2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_BLAKE2012;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = TRUE;
  NcmVector *x        = NULL;

  const gdouble z0 = 0.44;
  const gdouble z1 = 0.60;
  const gdouble z2 = 0.73;

  const gdouble bf0 = 0.0916;
  const gdouble bf1 = 0.0726;
  const gdouble bf2 = 0.0592;

  const gdouble icov00 = 24532.1;
  const gdouble icov01 = -25137.7;
  const gdouble icov02 = 12099.1;
  const gdouble icov10 = -25137.7;
  const gdouble icov11 = 134598.4;
  const gdouble icov12 = -64783.9;
  const gdouble icov20 = 12099.1;
  const gdouble icov21 = -64783.9;
  const gdouble icov22 = 128837.6;

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 2), ==, icov02);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 1), ==, icov11);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 2), ==, icov12);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 2, 0), ==, icov20);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 2, 1), ==, icov21);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 2, 2), ==, icov22);

  ncm_vector_free (x);
}

/* Kazin 2014 */

void
test_nc_data_bao_rdv_new_kazin2014 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_KAZIN2014;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_kazin2014 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_KAZIN2014;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = FALSE;
  NcmVector *x        = NULL;

  const gdouble z0 = 0.44;
  const gdouble z1 = 0.60;
  const gdouble z2 = 0.73;

  const gdouble bf0 = 11.550;
  const gdouble bf1 = 14.945;
  const gdouble bf2 = 16.932;

  const gdouble icov00 =  4.8116;
  const gdouble icov01 = -2.4651;
  const gdouble icov02 = 1.0375;
  const gdouble icov10 = -2.4651;
  const gdouble icov11 = 3.7697;
  const gdouble icov12 = -1.5865;
  const gdouble icov20 = 1.0375;
  const gdouble icov21 = -1.5865;
  const gdouble icov22 = 3.6498;

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (!R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 2), ==, icov02);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 1), ==, icov11);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 2), ==, icov12);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 2, 0), ==, icov20);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 2, 1), ==, icov21);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 2, 2), ==, icov22);

  ncm_vector_free (x);
}

/* DESI DR1 BGS and QSO 2024 */

void
test_nc_data_bao_rdv_new_desi_dr1_bgs_qso_2024 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_desi_dr1_bgs_qso_2024 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_DESI_DR1_BGS_QSO_2024;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = FALSE;
  NcmVector *x        = NULL;

  const gdouble z0 = 0.30;
  const gdouble z1 = 1.49;

  const gdouble bf0 = 12.56;
  const gdouble bf1 = 1.542;

  const gdouble icov00 = 44.444444444;
  const gdouble icov01 = 0.0;
  const gdouble icov10 = 0.0;
  const gdouble icov11 = 2.227667632;

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (!R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 1, 1), ==, icov11);

  ncm_vector_free (x);
}

/* DESI DR2 - BGS - 2025 */

void
test_nc_data_bao_rdv_new_desi_dr2_bgs_2025 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id   = NC_DATA_BAO_RDV_DESI_DR2_BGS_2025;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data     = NCM_DATA (nc_data_bao_rdv_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert_true (NC_IS_DATA_BAO_RDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_desi_dr2_bgs_2025 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv   = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id      = NC_DATA_BAO_RDV_DESI_DR2_BGS_2025;
  NcmVector *y        = ncm_data_gauss_peek_mean (gauss);
  NcmMatrix *inv_cov  = ncm_data_gauss_peek_inv_cov (gauss);
  gboolean R_DV       = FALSE;
  NcmVector *x        = NULL;

  const gdouble z0 = 0.295;

  const gdouble bf0 = 7.942;

  const gdouble icov00 = 1.0 / (0.075 * 0.075);

  g_assert_true (rdv != NULL);
  g_assert_true (NC_IS_DATA_BAO_RDV (rdv));

  g_assert_cmpuint (test->id, ==, id);
  g_assert_true (!R_DV);

  g_object_get (rdv, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (inv_cov, 0, 0), ==, icov00);

  ncm_vector_free (x);
}

void
test_nc_data_bao_rdv_get_distance (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDistance *dist  = NULL;
  NcDataBaoRDV *rdv = test->rdv;

  g_object_get (rdv, "dist", &dist, NULL);
  g_assert_true (dist != NULL);
  g_assert (NC_IS_DISTANCE (dist));

  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_get_redshift (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmVector *x      = NULL;
  NcDataBaoRDV *rdv = test->rdv;

  g_object_get (rdv, "z", &x, NULL);
  g_assert_true (x != NULL);
  g_assert (NCM_IS_VECTOR (x));

  ncm_vector_free (x);
}


/***************************************************************************
 *            test_nc_data_bao_dvr_dtdh.c
 *
 *  Wed Apr 16 14:23:38 2025
 *  Copyright  2025  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * test_nc_data_bao_dvr_dtdh.c
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

typedef struct _TestNcDataBaoDvrDtDh
{
  NcDataBaoDvrDtDh *dvdtdh;
  NcDataBaoId id;
} TestNcDataBaoDvrDtDh;

void test_nc_data_bao_dvr_dtdh_free (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);

void test_nc_data_bao_dvr_dtdh_new_desi_dr1_lrg_elg_2024 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);
void test_nc_data_bao_dvr_dtdh_set_sample_desi_dr1_lrg_elg_2024 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);

void test_nc_data_bao_dvr_dtdh_new_desi_dr2_lrg_elg_qso_lya_2025 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);
void test_nc_data_bao_dvr_dtdh_set_sample_desi_dr2_lrg_elg_qso_lya_2025 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);

void test_nc_data_bao_dvr_dtdh_get_distance (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);
void test_nc_data_bao_dvr_dtdh_get_redshift (TestNcDataBaoDvrDtDh *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_bao_dvr_dtdh/set_sample/desi2024", TestNcDataBaoDvrDtDh, NULL,
              &test_nc_data_bao_dvr_dtdh_new_desi_dr1_lrg_elg_2024,
              &test_nc_data_bao_dvr_dtdh_set_sample_desi_dr1_lrg_elg_2024,
              &test_nc_data_bao_dvr_dtdh_free);
  g_test_add ("/nc/data_bao_dvr_dtdh/set_sample/desi2025", TestNcDataBaoDvrDtDh, NULL,
              &test_nc_data_bao_dvr_dtdh_new_desi_dr2_lrg_elg_qso_lya_2025,
              &test_nc_data_bao_dvr_dtdh_set_sample_desi_dr2_lrg_elg_qso_lya_2025,
              &test_nc_data_bao_dvr_dtdh_free);
  g_test_add ("/nc/data_bao_dvr_dtdh/get_distance", TestNcDataBaoDvrDtDh, NULL,
              &test_nc_data_bao_dvr_dtdh_new_desi_dr1_lrg_elg_2024,
              &test_nc_data_bao_dvr_dtdh_get_distance,
              &test_nc_data_bao_dvr_dtdh_free);
  g_test_add ("/nc/data_bao_dvr_dtdh/get_redshift", TestNcDataBaoDvrDtDh, NULL,
              &test_nc_data_bao_dvr_dtdh_new_desi_dr1_lrg_elg_2024,
              &test_nc_data_bao_dvr_dtdh_get_redshift,
              &test_nc_data_bao_dvr_dtdh_free);

  g_test_run ();
}

void
test_nc_data_bao_dvr_dtdh_free (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcmData *dvdtdh = NCM_DATA (test->dvdtdh);

  NCM_TEST_FREE (ncm_data_free, dvdtdh);
}

/* DESI DR1 2024 */

void
test_nc_data_bao_dvr_dtdh_new_desi_dr1_lrg_elg_2024 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024;
  NcDistance *dist = nc_distance_new (3.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_dvr_dtdh_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->dvdtdh = NC_DATA_BAO_DVR_DTDH (data);
  g_assert_true (NC_IS_DATA_BAO_DVR_DTDH (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dvr_dtdh_set_sample_desi_dr1_lrg_elg_2024 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcDataBaoDvrDtDh *dvdtdh   = test->dvdtdh;
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (dvdtdh);
  NcDataBaoId id             = NC_DATA_BAO_DVR_DTDH_DESI_DR1_2024;
  NcmVector *y               = ncm_data_gauss_cov_peek_mean (gauss_cov);
  NcmMatrix *cov_m           = ncm_data_gauss_cov_peek_cov (gauss_cov);
  NcmVector *x               = NULL;

  const gdouble z0 = 0.51;
  const gdouble z1 = 0.51;
  const gdouble z2 = 0.71;
  const gdouble z3 = 0.71;
  const gdouble z4 = 0.93;
  const gdouble z5 = 0.93;
  const gdouble z6 = 1.32;
  const gdouble z7 = 1.32;

  const gdouble bf0 = 12.56;
  const gdouble bf1 = 1.542;
  const gdouble bf2 = 15.90;
  const gdouble bf3 = 1.193;
  const gdouble bf4 = 19.86;
  const gdouble bf5 = 0.824;
  const gdouble bf6 = 24.13;
  const gdouble bf7 = 0.498;


  gdouble covar[64] = {
    0.0225, -0.00420525, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.00420525, 0.003969, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.04, -0.004116, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.004116, 0.002401, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0289, -0.00145486, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.00145486, 0.000484, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1296, -0.00367632,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00367632, 0.000529
  };

  NcmMatrix *cov = ncm_matrix_new_data_static (covar, 8, 8);

  g_assert_true (dvdtdh != NULL);
  g_assert_true (NC_IS_DATA_BAO_DVR_DTDH (dvdtdh));

  g_assert_cmpuint (test->id, ==, id);

  g_object_get (dvdtdh, "z", &x, NULL);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (x, 3), ==, z3);
  ncm_assert_cmpdouble (ncm_vector_get (x, 4), ==, z4);
  ncm_assert_cmpdouble (ncm_vector_get (x, 5), ==, z5);
  ncm_assert_cmpdouble (ncm_vector_get (x, 6), ==, z6);
  ncm_assert_cmpdouble (ncm_vector_get (x, 7), ==, z7);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 3), ==, bf3);
  ncm_assert_cmpdouble (ncm_vector_get (y, 4), ==, bf4);
  ncm_assert_cmpdouble (ncm_vector_get (y, 5), ==, bf5);
  ncm_assert_cmpdouble (ncm_vector_get (y, 6), ==, bf6);
  ncm_assert_cmpdouble (ncm_vector_get (y, 7), ==, bf7);

  ncm_matrix_cmp (cov_m, cov, 1e-10);

  ncm_matrix_free (cov);
  ncm_vector_free (x);
}

/* DESI DR2 2025 */

void
test_nc_data_bao_dvr_dtdh_new_desi_dr2_lrg_elg_qso_lya_2025 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcDataBaoId id   = NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025;
  NcDistance *dist = nc_distance_new (3.0);
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_bao_dvr_dtdh_new_from_id (dist, id));
  g_assert_true (data != NULL);
  test->dvdtdh = NC_DATA_BAO_DVR_DTDH (data);
  g_assert_true (NC_IS_DATA_BAO_DVR_DTDH (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dvr_dtdh_set_sample_desi_dr2_lrg_elg_qso_lya_2025 (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcDataBaoDvrDtDh *dvdtdh   = test->dvdtdh;
  NcmDataGaussCov *gauss_cov = NCM_DATA_GAUSS_COV (dvdtdh);
  NcDataBaoId id             = NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025;
  NcmVector *y               = ncm_data_gauss_cov_peek_mean (gauss_cov);
  NcmMatrix *cov_m           = ncm_data_gauss_cov_peek_cov (gauss_cov);
  NcmVector *x               = NULL;


  gdouble zs[12]  = {0.51, 0.51, 0.706, 0.706, 0.934, 0.934, 1.321, 1.321, 1.484, 1.484, 2.33, 2.33};
  gdouble bfs[12] = {12.72, 0.622, 16.05, 0.892, 19.721, 1.223, 24.252, 1.948, 26.055, 2.386, 31.267, 4.518};

  const gdouble a11   =  0.099 * 0.099;
  const gdouble a12   =  0.099 * 0.017 * 0.05;
  const gdouble a22   =  0.017 * 0.017;
  const gdouble a33   =  0.11 * 0.11;
  const gdouble a34   =  0.11 * 0.021 * (-0.018);
  const gdouble a44   =  0.021 * 0.021;
  const gdouble a55   =  0.091 * 0.091;
  const gdouble a56   =  0.091 * 0.019 * 0.056;
  const gdouble a66   =  0.019 * 0.019;
  const gdouble a77   =  0.174 * 0.174;
  const gdouble a78   =  0.174 * 0.045 * 0.202;
  const gdouble a88   =  0.045 * 0.045;
  const gdouble a99   =  0.398 * 0.398;
  const gdouble a910  =  0.398 * 0.136 * 0.044;
  const gdouble a1010 =  0.136 * 0.136;
  const gdouble a1111 =  0.256 * 0.256;
  const gdouble a1112 =  0.256 * 0.097 * 0.574;
  const gdouble a1212 =  0.097 * 0.097;

  gdouble covar[144] = {
    a11, a12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    a12, a22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, a33, a34, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, a34, a44, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, a55, a56, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, a56, a66, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a77, a78, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a78, a88, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a99, a910, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a910, a1010, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a1111, a1112,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a1112, a1212
  };

  NcmVector *z   = ncm_vector_new_data_static (zs, 12, 1);
  NcmVector *bf  = ncm_vector_new_data_static (bfs, 12, 1);
  NcmMatrix *cov = ncm_matrix_new_data_static (covar, 12, 12);

  g_assert_true (dvdtdh != NULL);
  g_assert_true (NC_IS_DATA_BAO_DVR_DTDH (dvdtdh));

  g_assert_cmpuint (test->id, ==, id);

  g_object_get (dvdtdh, "z", &x, NULL);


  ncm_vector_cmp2 (z, x, 1e-10, 1e-10);
  ncm_vector_cmp2 (y, bf, 1e-10, 1e-10);
  ncm_matrix_cmp (cov_m, cov, 1e-10);

  ncm_matrix_free (cov);
  ncm_vector_free (z);
  ncm_vector_free (bf);
  ncm_vector_free (x);
}

void
test_nc_data_bao_dvr_dtdh_get_distance (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcDistance *dist         = NULL;
  NcDataBaoDvrDtDh *dvdtdh = test->dvdtdh;

  g_object_get (dvdtdh, "dist", &dist, NULL);
  g_assert_true (dist != NULL);
  g_assert (NC_IS_DISTANCE (dist));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dvr_dtdh_get_redshift (TestNcDataBaoDvrDtDh *test, gconstpointer pdata)
{
  NcmVector *x             = NULL;
  NcDataBaoDvrDtDh *dvdtdh = test->dvdtdh;

  g_object_get (dvdtdh, "z", &x, NULL);
  g_assert_true (x != NULL);
  g_assert (NCM_IS_VECTOR (x));

  ncm_vector_free (x);
}


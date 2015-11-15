/***************************************************************************
 *            test_nc_data_bao_dvdv.c
 *
 *  Thu February 12 21:09:27 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2015 <pennalima@gmail.com>
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

typedef struct _TestNcDataBaoDVDV
{
  NcDataBaoDVDV *dvdv;
  NcDataBaoId id;
} TestNcDataBaoDVDV;

void test_nc_data_bao_dvdv_free (TestNcDataBaoDVDV *test, gconstpointer pdata);

void test_nc_data_bao_dvdv_new_percival2007 (TestNcDataBaoDVDV *test, gconstpointer pdata);
void test_nc_data_bao_dvdv_set_sample_percival2007 (TestNcDataBaoDVDV *test, gconstpointer pdata);

void test_nc_data_bao_dvdv_new_percival2010 (TestNcDataBaoDVDV *test, gconstpointer pdata);
void test_nc_data_bao_dvdv_set_sample_percival2010 (TestNcDataBaoDVDV *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/numcosmo/nc_data_bao_dvdv/set_sample/percival2007", TestNcDataBaoDVDV, NULL,
              &test_nc_data_bao_dvdv_new_percival2007,
              &test_nc_data_bao_dvdv_set_sample_percival2007,
              &test_nc_data_bao_dvdv_free);
  g_test_add ("/numcosmo/nc_data_bao_dvdv/set_sample/percival2010", TestNcDataBaoDVDV, NULL,
              &test_nc_data_bao_dvdv_new_percival2010,
              &test_nc_data_bao_dvdv_set_sample_percival2010,
              &test_nc_data_bao_dvdv_free);

  g_test_run ();
}

void
test_nc_data_bao_dvdv_free (TestNcDataBaoDVDV *test, gconstpointer pdata)
{
  NcmData *dvdv = NCM_DATA(test->dvdv);
  NCM_TEST_FREE (ncm_data_free, dvdv);
}

/* Percival2007 */

void
test_nc_data_bao_dvdv_new_percival2007 (TestNcDataBaoDVDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_DVDV_PERCIVAL2007;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = NCM_DATA (nc_data_bao_dvdv_new_from_id (dist, id));
  g_assert (data != NULL);
  test->dvdv = NC_DATA_BAO_DVDV (data);
  g_assert (NC_IS_DATA_BAO_DVDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dvdv_set_sample_percival2007 (TestNcDataBaoDVDV *test, gconstpointer pdata)
{
  NcDataBaoDVDV *dvdv = test->dvdv;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (dvdv);
  NcDataBaoId id = NC_DATA_BAO_DVDV_PERCIVAL2007;

  const gdouble bf0 = 1.812;
  const gdouble sigma = 0.060;

  g_assert (dvdv != NULL);
  g_assert (NC_IS_DATA_BAO_DVDV (dvdv));

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (diag->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (diag->sigma, 0), ==, sigma);
}

/* Percival 2010 */

void
test_nc_data_bao_dvdv_new_percival2010 (TestNcDataBaoDVDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_DVDV_PERCIVAL2010;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = NCM_DATA (nc_data_bao_dvdv_new_from_id (dist, id));
  g_assert (data != NULL);
  test->dvdv = NC_DATA_BAO_DVDV (data);
  g_assert (NC_IS_DATA_BAO_DVDV (data));

  nc_distance_free (dist);
}

void
test_nc_data_bao_dvdv_set_sample_percival2010 (TestNcDataBaoDVDV *test, gconstpointer pdata)
{
  NcDataBaoDVDV *dvdv = test->dvdv;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (dvdv);
  NcDataBaoId id = NC_DATA_BAO_DVDV_PERCIVAL2010;

  const gdouble bf0 = 1.736;
  const gdouble sigma = 0.065;

  g_assert (dvdv != NULL);
  g_assert (NC_IS_DATA_BAO_DVDV (dvdv));

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (diag->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (diag->sigma, 0), ==, sigma);
}


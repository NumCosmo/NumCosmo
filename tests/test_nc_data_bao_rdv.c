/***************************************************************************
 *            test_nc_data_bao_rdv.c
 *
 *  Wed February 11 12:11:42 2015
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

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/percival2007", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_percival2007, 
              &test_nc_data_bao_rdv_set_sample_percival2007, 
              &test_nc_data_bao_rdv_free);
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/percival2010", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_percival2010, 
              &test_nc_data_bao_rdv_set_sample_percival2010, 
              &test_nc_data_bao_rdv_free);
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/beutler2011", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_beutler2011, 
              &test_nc_data_bao_rdv_set_sample_beutler2011, 
              &test_nc_data_bao_rdv_free);
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/padmanabhan2012", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_padmanabhan2012, 
              &test_nc_data_bao_rdv_set_sample_padmanabhan2012, 
              &test_nc_data_bao_rdv_free);
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/anderson2012", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_anderson2012, 
              &test_nc_data_bao_rdv_set_sample_anderson2012, 
              &test_nc_data_bao_rdv_free); 
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/blake2012", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_blake2012, 
              &test_nc_data_bao_rdv_set_sample_blake2012, 
              &test_nc_data_bao_rdv_free); 
  g_test_add ("/numcosmo/nc_data_bao_rdv/set_sample/kazin2014", TestNcDataBaoRDV, NULL, 
              &test_nc_data_bao_rdv_new_kazin2014, 
              &test_nc_data_bao_rdv_set_sample_kazin2014, 
              &test_nc_data_bao_rdv_free); 

  g_test_run ();
}

void _set_destroyed (gpointer b) { gboolean *destroyed = b; *destroyed = TRUE; }

void
test_nc_data_bao_rdv_free (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *rdv = NCM_DATA(test->rdv);
  gboolean destroyed = FALSE;
  g_object_set_data_full (G_OBJECT (rdv), "test-destroy", &destroyed, _set_destroyed);
  ncm_data_free (rdv);
  g_assert (destroyed);
}

/* Percival2007 */

void
test_nc_data_bao_rdv_new_percival2007 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_PERCIVAL2007;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_percival2007 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_PERCIVAL2007;
  gboolean R_DV = TRUE;
  
  const gdouble z0 = 0.2;
  const gdouble z1 = 0.35; 

  const gdouble bf0 = 0.1980;
  const gdouble bf1 = 0.1094;

  const gdouble icov00 = 35059.0;
  const gdouble icov01 = -24031.0;
  const gdouble icov10 = -24031.0;
  const gdouble icov11 = 108300.0;
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 1), ==, icov11);
}

/* Percival 2010 */

void
test_nc_data_bao_rdv_new_percival2010 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_PERCIVAL2010;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_percival2010 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_PERCIVAL2010;
  gboolean R_DV = TRUE;
  
  const gdouble z0 = 0.2;
  const gdouble z1 = 0.35; 

  const gdouble bf0 = 0.1905; 
  const gdouble bf1 = 0.1097;

  const gdouble icov00 = 30124.0; 
  const gdouble icov01 = -17227.0;
  const gdouble icov10 = -17227.0;
  const gdouble icov11 = 86977.0;
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 1), ==, icov11);
}

/* Beutler 2011 */

void
test_nc_data_bao_rdv_new_beutler2011 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_BEUTLER2011;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_beutler2011 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_BEUTLER2011;
  gboolean R_DV = TRUE;
  
  const gdouble z0 = 0.106;
  const gdouble bf0 = 0.336 / 1.027; 

  const gdouble icov00 = 1.0 * 1.027 * 1.027 / (0.015 * 0.015); 
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
}

/* Padmanabhan 2012 */

void
test_nc_data_bao_rdv_new_padmanabhan2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_PADMANABHAN2012;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_padmanabhan2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_PADMANABHAN2012;
  gboolean R_DV = FALSE;

  const gdouble z0 = 0.35;
  const gdouble bf0 = 8.88 * 1.025; 

  const gdouble icov00 = 1.0 / (0.17 * 0.17 * 1.025 * 1.025); 
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (!R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
}

/* Anderson 2012 */

void
test_nc_data_bao_rdv_new_anderson2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_ANDERSON2012;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_anderson2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_ANDERSON2012;
  gboolean R_DV = FALSE;

  const gdouble z0 = 0.57;
  const gdouble bf0 = 13.67; 

  const gdouble icov00 = 1.0 / (0.22 * 0.22); 
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (!R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
}

/* Blake 2012 */

void
test_nc_data_bao_rdv_new_blake2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_BLAKE2012;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_blake2012 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_BLAKE2012;
  gboolean R_DV = TRUE;
  
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
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 2), ==, icov02);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 1), ==, icov11);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 2), ==, icov12);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 2, 0), ==, icov20);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 2, 1), ==, icov21);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 2, 2), ==, icov22);
}

/* Kazin 2014 */

void
test_nc_data_bao_rdv_new_kazin2014 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcmData *data;
  NcDataBaoId id = NC_DATA_BAO_RDV_KAZIN2014;
  NcDistance *dist = nc_distance_new (2.0);

  test->id = id;
  data = nc_data_bao_rdv_new (dist, id);
  g_assert (data != NULL);
  test->rdv = NC_DATA_BAO_RDV (data);
  g_assert (NC_IS_DATA_BAO_RDV (data));
 
  nc_distance_free (dist);
}

void
test_nc_data_bao_rdv_set_sample_kazin2014 (TestNcDataBaoRDV *test, gconstpointer pdata)
{
  NcDataBaoRDV *rdv = test->rdv;
  NcmDataGauss *gauss = NCM_DATA_GAUSS (rdv);
  NcDataBaoId id = NC_DATA_BAO_RDV_KAZIN2014;
  gboolean R_DV = FALSE;
  
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
  
  g_assert (rdv != NULL);
  g_assert (NC_IS_DATA_BAO_RDV (rdv));
  
  g_assert_cmpuint (test->id, ==, id);
  g_assert (!R_DV);

  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (rdv->x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (gauss->y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 0), ==, icov00);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 1), ==, icov01);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 0, 2), ==, icov02);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 0), ==, icov10);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 1), ==, icov11);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 1, 2), ==, icov12);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 2, 0), ==, icov20);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 2, 1), ==, icov21);
  ncm_assert_cmpdouble (ncm_matrix_get (gauss->inv_cov, 2, 2), ==, icov22);
}


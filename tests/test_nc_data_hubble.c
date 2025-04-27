/***************************************************************************
 *            test_nc_data_hubble.c
 *
 *  Sat Apr 19 16:24:05 2025
 *  Copyright  2025  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * test_nc_data_hubble.c
 * Copyright (C) 2025 Mariana Penna-Lima <pennalima@unb.br>
 *
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

typedef struct _TestNcDataHubble
{
  NcDataHubble *hubble;
  NcDataHubbleId id;
} TestNcDataHubble;

void test_nc_data_hubble_free (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_simon2005 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_simon2005 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_cabre (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_cabre (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_stern2009 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_stern2009 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_moresco2012_bc03 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_moresco2012_bc03 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_moresco2012_mastro (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_moresco2012_mastro (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_moresco2015 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_moresco2015 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_moresco2016_bc03 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_moresco2016_bc03 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_moresco2016_mastro (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_moresco2016_mastro (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_busca2013 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_busca2013 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_riess2008 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_riess2008 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_zhang2012 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_zhang2012 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_riess2016 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_riess2016 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_ratsimbazafy2017 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_ratsimbazafy2017 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_gomez2018 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_gomez2018 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_riess2018 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_riess2018 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_borghi2022 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_borghi2022 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_jiao2023 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_jiao2023 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_jimenez2023 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_jimenez2023 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_new_tomasetti2023 (TestNcDataHubble *test, gconstpointer pdata);
void test_nc_data_hubble_set_sample_tomasetti2023 (TestNcDataHubble *test, gconstpointer pdata);

void test_nc_data_hubble_get_redshift (TestNcDataHubble *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/data_hubble/set_sample/simon2005", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_simon2005,
              &test_nc_data_hubble_set_sample_simon2005,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/cabre", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_cabre,
              &test_nc_data_hubble_set_sample_cabre,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/stern2009", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_stern2009,
              &test_nc_data_hubble_set_sample_stern2009,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/moresco2012bc03", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_moresco2012_bc03,
              &test_nc_data_hubble_set_sample_moresco2012_bc03,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/moresco2012mastro", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_moresco2012_mastro,
              &test_nc_data_hubble_set_sample_moresco2012_mastro,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/moresco2015", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_moresco2015,
              &test_nc_data_hubble_set_sample_moresco2015,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/moresco2016bc03", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_moresco2016_bc03,
              &test_nc_data_hubble_set_sample_moresco2016_bc03,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/moresco2016mastro", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_moresco2016_mastro,
              &test_nc_data_hubble_set_sample_moresco2016_mastro,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/busca2013", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_busca2013,
              &test_nc_data_hubble_set_sample_busca2013,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/riess2008", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_riess2008,
              &test_nc_data_hubble_set_sample_riess2008,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/zhang2012", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_zhang2012,
              &test_nc_data_hubble_set_sample_zhang2012,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/riess2016", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_riess2016,
              &test_nc_data_hubble_set_sample_riess2016,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/ratsimbazafy2017", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_ratsimbazafy2017,
              &test_nc_data_hubble_set_sample_ratsimbazafy2017,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/gomez12018", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_gomez2018,
              &test_nc_data_hubble_set_sample_gomez2018,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/riess2018", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_riess2018,
              &test_nc_data_hubble_set_sample_riess2018,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/borghi2022", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_borghi2022,
              &test_nc_data_hubble_set_sample_borghi2022,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/jiao2023", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_jiao2023,
              &test_nc_data_hubble_set_sample_jiao2023,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/jimenez2023", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_jimenez2023,
              &test_nc_data_hubble_set_sample_jimenez2023,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/set_sample/tomasetti2023", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_tomasetti2023,
              &test_nc_data_hubble_set_sample_tomasetti2023,
              &test_nc_data_hubble_free);
  g_test_add ("/nc/data_hubble/get_redshift", TestNcDataHubble, NULL,
              &test_nc_data_hubble_new_simon2005,
              &test_nc_data_hubble_get_redshift,
              &test_nc_data_hubble_free);

  g_test_run ();
}

void
test_nc_data_hubble_free (TestNcDataHubble *test, gconstpointer pdata)
{
  NcmData *hubble = NCM_DATA (test->hubble);

  NCM_TEST_FREE (ncm_data_free, hubble);
}

/* Simon 2005 */

void
test_nc_data_hubble_new_simon2005 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_SIMON2005;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_simon2005 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_SIMON2005;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  gdouble redshifts[9] = {
    0.09, 0.17, 0.27, 0.4, 0.88, 1.3, 1.43, 1.53, 1.75
  };

  gdouble bf[9] = {
    69.0, 83.0, 70.0, 87.0, 117.0, 168.0, 177.0, 140.0, 202.0
  };

  gdouble sigma[9] = {
    12.0, 8.3, 14.0, 17.4, 23.4, 13.4, 14.2, 14.0, 40.4
  };

  NcmVector *zs     = ncm_vector_new_data_static (redshifts, 9, 1);
  NcmVector *bfs    = ncm_vector_new_data_static (bf, 9, 1);
  NcmVector *sigmas = ncm_vector_new_data_static (sigma, 9, 1);

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_vector_cmp2 (zs, x, 1e-10, 0.0);
  ncm_vector_cmp2 (bfs, y, 1e-10, 0.0);
  ncm_vector_cmp2 (sigmas, sigma_v, 1e-10, 0.0);

  ncm_vector_free (x);
  ncm_vector_free (zs);
  ncm_vector_free (bfs);
  ncm_vector_free (sigmas);
}

/* Cabre */

void
test_nc_data_hubble_new_cabre (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_CABRE;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_cabre (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_CABRE;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 0.24;
  const gdouble z1 = 0.43;

  const gdouble bf0 = 83.2;
  const gdouble bf1 = 90.3;

  const gdouble sigma0 = 3.1;
  const gdouble sigma1 = 3.5;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);

  ncm_vector_free (x);
}

/* Stern 2009 */

void
test_nc_data_hubble_new_stern2009 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_STERN2009;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_stern2009 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_STERN2009;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  gdouble redshifts[11] = {
    0.1, 0.17, 0.27, 0.4, 0.48, 0.88, 0.9, 1.3, 1.43, 1.53, 1.75
  };

  gdouble bf[11] = {
    69.0, 83.0, 77.0, 95.0, 97.0, 90.0, 117.0, 168.0, 177.0, 140.0, 202.0
  };

  gdouble sigma[11] = {
    12.0, 8.0, 14.0, 17.0, 60.0, 40.0, 23.0, 17.0, 18.0, 14.0, 40.0
  };

  NcmVector *zs     = ncm_vector_new_data_static (redshifts, 11, 1);
  NcmVector *bfs    = ncm_vector_new_data_static (bf, 11, 1);
  NcmVector *sigmas = ncm_vector_new_data_static (sigma, 11, 1);

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_vector_cmp2 (zs, x, 1e-10, 0.0);
  ncm_vector_cmp2 (bfs, y, 1e-10, 0.0);
  ncm_vector_cmp2 (sigmas, sigma_v, 1e-10, 0.0);

  ncm_vector_free (x);
  ncm_vector_free (zs);
  ncm_vector_free (bfs);
  ncm_vector_free (sigmas);
}

/* Moresco 2012 BC03 */

void
test_nc_data_hubble_new_moresco2012_bc03 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_MORESCO2012_BC03;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_moresco2012_bc03 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_MORESCO2012_BC03;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 0.1791;
  const gdouble z1 = 0.1993;
  const gdouble z2 = 0.3519;
  const gdouble z3 = 0.5929;
  const gdouble z4 = 0.6797;
  const gdouble z5 = 0.7812;
  const gdouble z6 = 0.8754;
  const gdouble z7 = 1.037;

  const gdouble bf0 = 75.0;
  const gdouble bf1 = 75.0;
  const gdouble bf2 = 83.0;
  const gdouble bf3 = 104.0;
  const gdouble bf4 = 92.0;
  const gdouble bf5 = 105.0;
  const gdouble bf6 = 125.0;
  const gdouble bf7 = 154.0;

  const gdouble sigma0 = 4.0;
  const gdouble sigma1 = 5.0;
  const gdouble sigma2 = 14.0;
  const gdouble sigma3 = 13.0;
  const gdouble sigma4 = 8.0;
  const gdouble sigma5 = 12.0;
  const gdouble sigma6 = 17.0;
  const gdouble sigma7 = 20.0;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

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

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 2), ==, sigma2);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 3), ==, sigma3);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 4), ==, sigma4);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 5), ==, sigma5);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 6), ==, sigma6);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 7), ==, sigma7);

  ncm_vector_free (x);
}

/* Moresco 2012 MASTRO */

void
test_nc_data_hubble_new_moresco2012_mastro (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_MORESCO2012_MASTRO;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_moresco2012_mastro (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_MORESCO2012_MASTRO;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 0.1791;
  const gdouble z1 = 0.1993;
  const gdouble z2 = 0.3519;
  const gdouble z3 = 0.5929;
  const gdouble z4 = 0.6797;
  const gdouble z5 = 0.7812;
  const gdouble z6 = 0.8754;
  const gdouble z7 = 1.037;

  const gdouble bf0 = 81.0;
  const gdouble bf1 = 81.0;
  const gdouble bf2 = 88.0;
  const gdouble bf3 = 110.0;
  const gdouble bf4 = 98.0;
  const gdouble bf5 = 88.0;
  const gdouble bf6 = 124.0;
  const gdouble bf7 = 113.0;

  const gdouble sigma0 = 5.0;
  const gdouble sigma1 = 6.0;
  const gdouble sigma2 = 16.0;
  const gdouble sigma3 = 15.0;
  const gdouble sigma4 = 10.0;
  const gdouble sigma5 = 11.0;
  const gdouble sigma6 = 17.0;
  const gdouble sigma7 = 15.0;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

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

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 2), ==, sigma2);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 3), ==, sigma3);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 4), ==, sigma4);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 5), ==, sigma5);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 6), ==, sigma6);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 7), ==, sigma7);

  ncm_vector_free (x);
}

/* Moresco 2015 */

void
test_nc_data_hubble_new_moresco2015 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_MORESCO2015;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_moresco2015 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_MORESCO2015;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 1.363;
  const gdouble z1 = 1.965;

  const gdouble bf0 = 160.0;
  const gdouble bf1 = 186.5;

  const gdouble sigma0 = 33.6;
  const gdouble sigma1 = 50.4;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);

  ncm_vector_free (x);
}

/* Moresco 2016 BC03 */

void
test_nc_data_hubble_new_moresco2016_bc03 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_MORESCO2016_DR9_BC03;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_moresco2016_bc03 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_MORESCO2016_DR9_BC03;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 0.3802;
  const gdouble z1 = 0.4004;
  const gdouble z2 = 0.4247;
  const gdouble z3 = 0.4497;
  const gdouble z4 = 0.4783;

  const gdouble bf0 = 83.0;
  const gdouble bf1 = 77.0;
  const gdouble bf2 = 87.1;
  const gdouble bf3 = 92.8;
  const gdouble bf4 = 80.9;

  const gdouble sigma0 = 13.5;
  const gdouble sigma1 = 10.2;
  const gdouble sigma2 = 11.2;
  const gdouble sigma3 = 12.9;
  const gdouble sigma4 = 9.0;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (x, 3), ==, z3);
  ncm_assert_cmpdouble (ncm_vector_get (x, 4), ==, z4);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 3), ==, bf3);
  ncm_assert_cmpdouble (ncm_vector_get (y, 4), ==, bf4);

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 2), ==, sigma2);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 3), ==, sigma3);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 4), ==, sigma4);

  ncm_vector_free (x);
}

/* Moresco 2016 Mastro */

void
test_nc_data_hubble_new_moresco2016_mastro (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_MORESCO2016_DR9_MASTRO;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_moresco2016_mastro (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_MORESCO2016_DR9_MASTRO;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 0.3802;
  const gdouble z1 = 0.4004;
  const gdouble z2 = 0.4247;
  const gdouble z3 = 0.4497;
  const gdouble z4 = 0.4783;

  const gdouble bf0 = 89.3;
  const gdouble bf1 = 82.8;
  const gdouble bf2 = 93.7;
  const gdouble bf3 = 99.7;
  const gdouble bf4 = 86.6;

  const gdouble sigma0 = 14.1;
  const gdouble sigma1 = 10.6;
  const gdouble sigma2 = 11.7;
  const gdouble sigma3 = 13.4;
  const gdouble sigma4 = 8.7;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (x, 3), ==, z3);
  ncm_assert_cmpdouble (ncm_vector_get (x, 4), ==, z4);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 3), ==, bf3);
  ncm_assert_cmpdouble (ncm_vector_get (y, 4), ==, bf4);

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 2), ==, sigma2);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 3), ==, sigma3);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 4), ==, sigma4);

  ncm_vector_free (x);
}

/* BUSCA 2013 BAO WMAP */

void
test_nc_data_hubble_new_busca2013 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_BUSCA2013_BAO_WMAP;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_busca2013 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_BUSCA2013_BAO_WMAP;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 2.3;
  const gdouble bf0    = 223.74;
  const gdouble sigma0 = 7.92;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* RIESS 2008 */

void
test_nc_data_hubble_new_riess2008 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_RIESS2008_HST;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_riess2008 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_RIESS2008_HST;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.0;
  const gdouble bf0    = 73.8;
  const gdouble sigma0 = 2.4;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Zhang 2012 */

void
test_nc_data_hubble_new_zhang2012 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_ZHANG2012;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_zhang2012 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_ZHANG2012;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0 = 0.07;
  const gdouble z1 = 0.12;
  const gdouble z2 = 0.2;
  const gdouble z3 = 0.28;

  const gdouble bf0 = 69.0;
  const gdouble bf1 = 68.6;
  const gdouble bf2 = 72.9;
  const gdouble bf3 = 88.8;

  const gdouble sigma0 = 19.6;
  const gdouble sigma1 = 26.2;
  const gdouble sigma2 = 29.6;
  const gdouble sigma3 = 36.6;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (x, 1), ==, z1);
  ncm_assert_cmpdouble (ncm_vector_get (x, 2), ==, z2);
  ncm_assert_cmpdouble (ncm_vector_get (x, 3), ==, z3);

  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 1), ==, bf1);
  ncm_assert_cmpdouble (ncm_vector_get (y, 2), ==, bf2);
  ncm_assert_cmpdouble (ncm_vector_get (y, 3), ==, bf3);

  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 1), ==, sigma1);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 2), ==, sigma2);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 3), ==, sigma3);

  ncm_vector_free (x);
}

/* RIESS 2016 */

void
test_nc_data_hubble_new_riess2016 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_RIESS2016_HST_WFC3;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_riess2016 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_RIESS2016_HST_WFC3;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.0;
  const gdouble bf0    = 73.24;
  const gdouble sigma0 = 1.74;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Ratsimbazafy 2017 */

void
test_nc_data_hubble_new_ratsimbazafy2017 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_RATSIMBAZAFY2017;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_ratsimbazafy2017 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_RATSIMBAZAFY2017;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.47;
  const gdouble bf0    = 89.0;
  const gdouble sigma0 = 49.6;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Gomez Valent 2018 */

void
test_nc_data_hubble_new_gomez2018 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_GOMEZ_VALENT_COMP2018;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_gomez2018 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_GOMEZ_VALENT_COMP2018;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  gdouble redshifts[31] = {
    0.0700, 0.0900, 0.1200, 0.1700, 0.1791, 0.1993, 0.2000, 0.2700, 0.2800, 0.3519, 0.3802,
    0.4000, 0.4004, 0.4247, 0.4497, 0.4700, 0.4783, 0.4800, 0.5929, 0.6797, 0.7812, 0.8754,
    0.8800, 0.9000, 1.037, 1.3000, 1.3630, 1.4300, 1.5300, 1.7500, 1.965
  };

  gdouble bf[31] = {
    69.0,   69.0,   68.6,   83.0,   75.0,   75.0,   72.9,   77.0,   88.8,   83.0,   83.0,
    95.0,   77.0,   87.1,   92.8,   89.0,   80.9, 97.0,  104.0,   92.0,  105.0,  125.0,
    90.0,  117.0, 154.0,  168.0,  160.0,  177.0,  140.0,  202.0, 186.5
  };

  gdouble sigma[31] = {
    19.6,   12.0,   26.2,    8.0,    4.0,    5.0,   29.6,   14.0,   36.6,   14.0,   13.5,
    17.0,   10.2,   11.2,   12.9,   49.6,    9.0, 62.0,   13.0,    8.0,   12.0,   17.0,
    40.0,   23.0,  20.0,   17.0,   33.6,   18.0,   14.0,   40.0,  50.4
  };

  NcmVector *zs     = ncm_vector_new_data_static (redshifts, 31, 1);
  NcmVector *bfs    = ncm_vector_new_data_static (bf, 31, 1);
  NcmVector *sigmas = ncm_vector_new_data_static (sigma, 31, 1);

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_vector_cmp2 (zs, x, 1e-10, 0.0);
  ncm_vector_cmp2 (bfs, y, 1e-10, 0.0);
  ncm_vector_cmp2 (sigmas, sigma_v, 1e-10, 0.0);

  ncm_vector_free (x);
  ncm_vector_free (zs);
  ncm_vector_free (bfs);
  ncm_vector_free (sigmas);
}

/* RIESS 2018 */

void
test_nc_data_hubble_new_riess2018 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_RIESS2018;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_riess2018 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_RIESS2018;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.0;
  const gdouble bf0    = 73.45;
  const gdouble sigma0 = 1.66;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Borghi 2022 */

void
test_nc_data_hubble_new_borghi2022 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_BORGHI2022;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_borghi2022 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_BORGHI2022;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.75;
  const gdouble bf0    = 98.8;
  const gdouble sigma0 = 33.6;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Jiao 2023 */

void
test_nc_data_hubble_new_jiao2023 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_JIAO2023;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_jiao2023 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_JIAO2023;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.8;
  const gdouble bf0    = 113.1;
  const gdouble sigma0 = 25.22;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Jimenez 2023 */

void
test_nc_data_hubble_new_jimenez2023 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_JIMENEZ2023;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_jimenez2023 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_JIMENEZ2023;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 0.75;
  const gdouble bf0    = 105.0;
  const gdouble sigma0 = 10.76;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

/* Tomasetti 2023 */

void
test_nc_data_hubble_new_tomasetti2023 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubbleId id = NC_DATA_HUBBLE_TOMASETTI2023;
  NcmData *data;

  test->id = id;
  data     = NCM_DATA (nc_data_hubble_new_from_id (id));
  g_assert_true (data != NULL);
  test->hubble = NC_DATA_HUBBLE (data);
  g_assert_true (NC_IS_DATA_HUBBLE (data));
}

void
test_nc_data_hubble_set_sample_tomasetti2023 (TestNcDataHubble *test, gconstpointer pdata)
{
  NcDataHubble *hubble   = test->hubble;
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (hubble);
  NcDataHubbleId id      = NC_DATA_HUBBLE_TOMASETTI2023;
  NcmVector *y           = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *sigma_v     = ncm_data_gauss_diag_peek_std (diag);
  NcmVector *x           = NULL;

  const gdouble z0     = 1.26;
  const gdouble bf0    = 135.0;
  const gdouble sigma0 = 65.0;

  g_assert_true (hubble != NULL);
  g_assert_true (NC_IS_DATA_HUBBLE (hubble));
  g_object_get (hubble, "z", &x, NULL);

  g_assert_cmpuint (test->id, ==, id);

  ncm_assert_cmpdouble (ncm_vector_get (x, 0), ==, z0);
  ncm_assert_cmpdouble (ncm_vector_get (y, 0), ==, bf0);
  ncm_assert_cmpdouble (ncm_vector_get (sigma_v, 0), ==, sigma0);

  ncm_vector_free (x);
}

void
test_nc_data_hubble_get_redshift (TestNcDataHubble *test, gconstpointer pdata)
{
  NcmVector *x         = NULL;
  NcDataHubble *hubble = test->hubble;

  g_object_get (hubble, "z", &x, NULL);
  g_assert_true (x != NULL);
  g_assert (NCM_IS_VECTOR (x));

  ncm_vector_free (x);
}


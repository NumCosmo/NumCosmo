/***************************************************************************
 *            test_nc_hicosmo_qgw.c
 *
 *  Tue Feb 2 16:05:43 2024
 *  Copyright  2024  Eduardo Barroso
 *  <pennalima@gmail.com>
 ***************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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

typedef struct _TestNcHIcosmoQGW
{
  NcHICosmo *cosmo;
} TestNcHIcosmoQGW;

void test_nc_hicosmo_qgw_new (TestNcHIcosmoQGW *test, gconstpointer pdata);
void test_nc_hicosmo_qgw_units (TestNcHIcosmoQGW *test, gconstpointer pdata);
void test_nc_hicosmo_qgw_free (TestNcHIcosmoQGW *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/hicosmo_qgw/units", TestNcHIcosmoQGW, NULL,
              &test_nc_hicosmo_qgw_new,
              &test_nc_hicosmo_qgw_units,
              &test_nc_hicosmo_qgw_free);

  g_test_run ();
}

void
test_nc_hicosmo_qgw_free (TestNcHIcosmoQGW *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_hicosmo_free, test->cosmo);
}

void
test_nc_hicosmo_qgw_new (TestNcHIcosmoQGW *test, gconstpointer pdata)
{
  test->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  g_assert_true (test->cosmo != NULL);

  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_H0,       70.0);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_C,   0.255);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X,   0.7);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.7245);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_B,   0.045);
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_XCDM_W,   -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->cosmo), "Omegak", 0.0, NULL);
}

void
test_nc_hicosmo_qgw_units (TestNcHIcosmoQGW *test, gconstpointer pdata)
{
  /* NcHICosmo *cosmo = test->cosmo; */
  /* gdouble units    = 1.0; */
  /* g_assert_true(NC_IS_HICOSMO_QGW(cosmo)); */
  /* nc_hicosmo_qgw_set_units(cosmo, units); */
}


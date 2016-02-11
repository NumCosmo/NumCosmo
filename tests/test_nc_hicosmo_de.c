/***************************************************************************
 *            test_nc_recomb.c
 *
 *  Wed November 14 11:38:56 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

typedef struct _TestNcHICosmoDE
{
  NcHICosmo *cosmo;
} TestNcHICosmoDE;

void test_nc_hicosmo_de_xcdm_new (TestNcHICosmoDE *test, gconstpointer pdata);
void test_nc_hicosmo_de_free (TestNcHICosmoDE *test, gconstpointer pdata);

void test_nc_hicosmo_de_omega_x2omega_k (TestNcHICosmoDE *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/numcosmo/nc_hicosmo_de/omega_x2omega_k", TestNcHICosmoDE, NULL,
              &test_nc_hicosmo_de_xcdm_new,
              &test_nc_hicosmo_de_omega_x2omega_k,
              &test_nc_hicosmo_de_free);

  g_test_run ();
}

void
test_nc_hicosmo_de_xcdm_new (TestNcHICosmoDE *test, gconstpointer pdata)
{
  test->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  g_assert (test->cosmo != NULL);
  g_assert (NC_IS_HICOSMO_DE (test->cosmo));
  g_assert (NC_IS_HICOSMO_DE_XCDM (test->cosmo));
}

void
test_nc_hicosmo_de_free (TestNcHICosmoDE *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = test->cosmo;
  NCM_TEST_FREE (nc_hicosmo_free, cosmo);
}

void
test_nc_hicosmo_de_omega_x2omega_k (TestNcHICosmoDE *test, gconstpointer pdata)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (test->cosmo);
  nc_hicosmo_de_omega_x2omega_k (cosmo_de);

  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    const gdouble pOmega_k0 = ncm_model_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X);
    ncm_assert_cmpdouble_e (Omega_k0, ==, pOmega_k0, 1.0e-7);
  }
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    const gdouble pOmega_k0 = ncm_model_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X);
    ncm_assert_cmpdouble_e (Omega_k0, ==, pOmega_k0, 1.0e-7);
  }
  ncm_model_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    const gdouble pOmega_k0 = ncm_model_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X);
    ncm_assert_cmpdouble_e (Omega_k0, ==, pOmega_k0, 1.0e-7);
  }
  ncm_model_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0, 3.0);
  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    ncm_assert_cmpdouble_e (Omega_k0, ==, 0.0, 1.0e-7);
  }
}

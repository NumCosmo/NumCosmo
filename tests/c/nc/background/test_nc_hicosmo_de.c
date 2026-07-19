/***************************************************************************
 *            test_nc_recomb.c
 *
 *  Wed November 14 11:38:56 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

void test_nc_hicosmo_de_xcdm_new_full (void);
void test_nc_hicosmo_lcdm_new_full (void);
void test_nc_hicosmo_de_cpl_new_full (void);
void test_nc_hicosmo_de_jbp_new_full (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/hicosmo_de/omega_x2omega_k", TestNcHICosmoDE, NULL,
              &test_nc_hicosmo_de_xcdm_new,
              &test_nc_hicosmo_de_omega_x2omega_k,
              &test_nc_hicosmo_de_free);

  g_test_add_func ("/nc/hicosmo_de/xcdm/new_full", &test_nc_hicosmo_de_xcdm_new_full);
  g_test_add_func ("/nc/hicosmo_de/lcdm/new_full", &test_nc_hicosmo_lcdm_new_full);
  g_test_add_func ("/nc/hicosmo_de/cpl/new_full", &test_nc_hicosmo_de_cpl_new_full);
  g_test_add_func ("/nc/hicosmo_de/jbp/new_full", &test_nc_hicosmo_de_jbp_new_full);

  g_test_run ();
}

void
test_nc_hicosmo_de_xcdm_new (TestNcHICosmoDE *test, gconstpointer pdata)
{
  test->cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  g_assert_true (test->cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE (test->cosmo));
  g_assert_true (NC_IS_HICOSMO_DE_XCDM (test->cosmo));
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

  nc_hicosmo_de_omega_x2omega_k (cosmo_de, NULL);

  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    const gdouble pOmega_k0 = ncm_model_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X);

    ncm_assert_cmpdouble_e (Omega_k0, ==, pOmega_k0, 1.0e-7, 0.0);
  }
  ncm_model_orig_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    const gdouble pOmega_k0 = ncm_model_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X);

    ncm_assert_cmpdouble_e (Omega_k0, ==, pOmega_k0, 1.0e-7, 0.0);
  }
  ncm_model_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
  {
    const gdouble Omega_k0  = nc_hicosmo_Omega_k0 (test->cosmo);
    const gdouble pOmega_k0 = ncm_model_param_get (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_OMEGA_X);

    ncm_assert_cmpdouble_e (Omega_k0, ==, pOmega_k0, 1.0e-7, 0.0);
  }
  ncm_model_param_set (NCM_MODEL (test->cosmo), NC_HICOSMO_DE_T_GAMMA0, 3.0);
  {
    const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (test->cosmo);

    ncm_assert_cmpdouble_e (Omega_k0, ==, 0.0, 1.0e-7, 0.0);
  }
}

void
test_nc_hicosmo_de_xcdm_new_full (void)
{
  NcHIReion *reion       = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim         = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmoDEXcdm *cosmo = nc_hicosmo_de_xcdm_new_full (reion, prim);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE_XCDM (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  /* Both submodels are optional -- NULL/NULL must behave like the plain
   * zero-arg constructor. */
  cosmo = nc_hicosmo_de_xcdm_new_full (NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_lcdm_new_full (void)
{
  NcHIReion *reion     = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim       = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmoLCDM *cosmo = nc_hicosmo_lcdm_new_full (reion, prim);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_LCDM (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  cosmo = nc_hicosmo_lcdm_new_full (NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_de_cpl_new_full (void)
{
  NcHIReion *reion      = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim        = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmoDECpl *cosmo = nc_hicosmo_de_cpl_new_full (reion, prim);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE_CPL (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  cosmo = nc_hicosmo_de_cpl_new_full (NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_de_jbp_new_full (void)
{
  NcHIReion *reion      = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim        = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmoDEJbp *cosmo = nc_hicosmo_de_jbp_new_full (reion, prim);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE_JBP (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  cosmo = nc_hicosmo_de_jbp_new_full (NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}


/***************************************************************************
 *            test_nc_cbe.c
 *
 *  Thu January 05 19:23:54 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2017 <sandro@isoftware.com.br>
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

typedef struct _TestNcCBE
{
  NcCBE *cbe;
  NcHICosmo *cosmo;
  guint ntests;
} TestNcCBE;

static void test_nc_cbe_lcdm_new (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_xcdm_new (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_mnu_lcdm_new (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_mnu_xcdm_new (TestNcCBE *test, gconstpointer pdata);

static void test_nc_cbe_flat_lcdm_new (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_flat_xcdm_new (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_flat_mnu_lcdm_new (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_flat_mnu_xcdm_new (TestNcCBE *test, gconstpointer pdata);

static void test_nc_cbe_pad_new (TestNcCBE *test, gconstpointer pdata);

static void test_nc_cbe_free (TestNcCBE *test, gconstpointer pdata);

static void test_nc_cbe_sanity (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_compare_bg (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_serialize (TestNcCBE *test, gconstpointer pdata);
static void test_nc_cbe_calc_ps (TestNcCBE *test, gconstpointer pdata);

static void test_nc_cbe_traps (TestNcCBE *test, gconstpointer pdata);

/*static void test_nc_cbe_invalid_model (TestNcCBE *test, gconstpointer pdata);*/

typedef struct _TestNcCBEFunc
{
  void (*func) (TestNcCBE *, gconstpointer);
  const gchar *name;
  gconstpointer pdata;
} TestNcCBEFunc;

#define TEST_NC_CBE_MODEL_LEN 8
TestNcCBEFunc models[TEST_NC_CBE_MODEL_LEN] =
{
  {test_nc_cbe_lcdm_new,          "lcdm",          NULL},
  {test_nc_cbe_xcdm_new,          "xcdm",          NULL},
  {test_nc_cbe_mnu_lcdm_new,      "lcdm/mnu",      NULL},
  {test_nc_cbe_mnu_xcdm_new,      "xcdm/mnu",      NULL},
  {test_nc_cbe_flat_lcdm_new,     "lcdm/flat",     NULL},
  {test_nc_cbe_flat_xcdm_new,     "xcdm/flat",     NULL},
  {test_nc_cbe_flat_mnu_lcdm_new, "lcdm/flat/mnu", NULL},
  {test_nc_cbe_flat_mnu_xcdm_new, "xcdm/flat/mnu", NULL},
};

#define TEST_NC_CBE_TEST_LEN 4
TestNcCBEFunc tests[TEST_NC_CBE_TEST_LEN] =
{
  {test_nc_cbe_compare_bg, "compare_bg", NULL},
  {test_nc_cbe_sanity,     "sanity",     NULL},
  {test_nc_cbe_serialize,  "serialize",  NULL},
  {test_nc_cbe_calc_ps,    "calc_ps",    NULL},
};

gint
main (gint argc, gchar *argv[])
{
  gint i, j;
  
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  
  for (i = 0; i < TEST_NC_CBE_MODEL_LEN; i++)
  {
    for (j = 0; j < TEST_NC_CBE_TEST_LEN; j++)
    {
      gchar *path = g_strdup_printf ("/nc/cbe/%s/%s", models[i].name, tests[j].name);
      
      g_test_add (path, TestNcCBE, models[i].pdata, models[i].func, tests[j].func, &test_nc_cbe_free);
      
      g_free (path);
      
    }
  }
  
  g_test_add ("/nc/cbe/traps", TestNcCBE, NULL,
              &test_nc_cbe_lcdm_new,
              &test_nc_cbe_traps,
              &test_nc_cbe_free);
  
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_add ("/nc/cbe/invalid/model/subprocess", TestNcCBE, NULL,
              &test_nc_cbe_pad_new,
              &test_nc_cbe_compare_bg,
              &test_nc_cbe_free);
#endif
  g_test_run ();
}

void
test_nc_cbe_lcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.0>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_xcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.1>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_mnu_lcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.0>, 'massnu-length' : <1>, 'massnu' : <[0.6]>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_mnu_xcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.1>, 'massnu-length' : <1>, 'massnu' : <[0.6]>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_flat_lcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.0>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));
  ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_flat_xcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.1>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));
  ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_flat_mnu_lcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.0>, 'massnu-length' : <1>, 'massnu' : <[0.6]>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));
  ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_flat_mnu_xcdm_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEXcdm{'w' : <-1.1>, 'massnu-length' : <1>, 'massnu' : <[0.6]>}");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));
  ncm_model_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.0);
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_pad_new (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = nc_cbe_new ();
  NcHICosmo *cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO_DE, "NcHICosmoDEPad");
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim  *prim  = NC_HIPRIM  (nc_hiprim_power_law_new ());
  
  test->cbe    = cbe;
  test->cosmo  = cosmo;
  test->ntests = 1000;
  
  g_assert_true (cbe != NULL);
  g_assert_true (NC_IS_CBE (cbe));
  
  g_assert_true (NC_IS_HICOSMO (cosmo));
  g_assert_true (NC_IS_HIREION (reion));
  g_assert_true (NC_IS_HIPRIM  (prim));
  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
}

void
test_nc_cbe_free (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = test->cbe;
  NcHICosmo *cosmo = test->cosmo;
  
  NCM_TEST_FREE (nc_cbe_free, cbe);
  NCM_TEST_FREE (nc_hicosmo_free, cosmo);
}

static void
test_nc_cbe_sanity (TestNcCBE *test, gconstpointer pdata)
{
  nc_cbe_ref (test->cbe);
  nc_cbe_free (test->cbe);

  g_assert_true (NC_IS_CBE (test->cbe));
}

void
test_nc_cbe_compare_bg (TestNcCBE *test, gconstpointer pdata)
{
  NcCBE *cbe       = test->cbe;
  NcHICosmo *cosmo = test->cosmo;
  
  /*guint ntests     = test->ntests;*/
  
  nc_cbe_prepare (cbe, cosmo);
  {
    const gdouble err = nc_cbe_compare_bg (cbe, cosmo, FALSE);
    
    g_assert_cmpfloat (err, <, 1.0e-4);
  }
}

static void
test_nc_cbe_serialize (TestNcCBE *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  gchar *cbe_ser = ncm_serialize_to_string (ser, G_OBJECT (test->cbe), TRUE);

  nc_cbe_free (test->cbe);

  test->cbe = NC_CBE (ncm_serialize_from_string (ser, cbe_ser));

  g_free (cbe_ser);

  g_assert_true (NC_IS_CBE (test->cbe));

  test_nc_cbe_compare_bg (test, pdata);
}

static void
test_nc_cbe_calc_ps (TestNcCBE *test, gconstpointer pdata)
{
  nc_cbe_set_calc_transfer (test->cbe, TRUE);
  nc_cbe_set_max_matter_pk_z (test->cbe, 1.0);

  nc_cbe_prepare (test->cbe, test->cosmo);
  nc_cbe_prepare_if_needed (test->cbe, test->cosmo);

  ncm_model_state_mark_outdated (NCM_MODEL (test->cosmo));
  nc_cbe_prepare_if_needed (test->cbe, test->cosmo);

  {
    NcmSpline2d *ps = nc_cbe_get_matter_ps (test->cbe);
    ncm_spline2d_free (ps);
  }

}

void
test_nc_cbe_traps (TestNcCBE *test, gconstpointer pdata)
{
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_trap_subprocess ("/nc/cbe/invalid/model/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}


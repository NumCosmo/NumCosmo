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
void test_nc_hicosmo_bbn_construct_fixed (void);
void test_nc_hicosmo_bbn_construct_fixed_subprocess (void);
void test_nc_hicosmo_late_attach (void);
void test_nc_hicosmo_lcdm_new_full (void);
void test_nc_hicosmo_de_cpl_new_full (void);
void test_nc_hicosmo_de_jbp_new_full (void);

void test_nc_hicosmo_de_param_by_name_submodel (void);
void test_nc_hicosmo_de_reparam_old_layout (void);
void test_nc_hicosmo_de_reparam_old_layout_refused (void);

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
  g_test_add_func ("/nc/hicosmo_de/bbn/construct_fixed", &test_nc_hicosmo_bbn_construct_fixed);
  g_test_add_func ("/nc/hicosmo_de/bbn/construct_fixed/subprocess", &test_nc_hicosmo_bbn_construct_fixed_subprocess);
  g_test_add_func ("/nc/hicosmo_de/late_attach", &test_nc_hicosmo_late_attach);
  g_test_add_func ("/nc/hicosmo_de/lcdm/new_full", &test_nc_hicosmo_lcdm_new_full);
  g_test_add_func ("/nc/hicosmo_de/cpl/new_full", &test_nc_hicosmo_de_cpl_new_full);
  g_test_add_func ("/nc/hicosmo_de/jbp/new_full", &test_nc_hicosmo_de_jbp_new_full);

  g_test_add_func ("/nc/hicosmo_de/param_by_name/submodel", &test_nc_hicosmo_de_param_by_name_submodel);
  g_test_add_func ("/nc/hicosmo_de/reparam/old_layout", &test_nc_hicosmo_de_reparam_old_layout);
  g_test_add_func ("/nc/hicosmo_de/reparam/old_layout/refused", &test_nc_hicosmo_de_reparam_old_layout_refused);

  g_test_run ();
}

/*
 * An NcHICosmoDEReparamOk of length whose params-desc describes parameter key alone,
 * as a stored params-desc does: it replaces the descriptions the constructor set.
 */
static NcmReparam *
_test_nc_hicosmo_de_reparam_with_desc (guint length, gint key)
{
  NcmReparam *reparam = NCM_REPARAM (nc_hicosmo_de_reparam_ok_new (length));
  NcmObjDictInt *desc = ncm_obj_dict_int_new ();
  NcmSParam *sp       = ncm_sparam_new ("Bogus", "B", 0.0, 1.0, 0.1, 0.0, 0.5, NCM_PARAM_TYPE_FIXED);

  ncm_obj_dict_int_add (desc, key, G_OBJECT (sp));
  g_object_set (reparam, "params-desc", desc, NULL);

  ncm_sparam_free (sp);
  ncm_obj_dict_int_unref (desc);

  return reparam;
}

/* Setting the reparam of length describing key on an Xcdm fails with reason, leaving none. */
static void
_test_nc_hicosmo_de_reparam_refused (guint length, gint key, const gchar *reason)
{
  NcHICosmo *cosmo    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcmReparam *reparam = _test_nc_hicosmo_de_reparam_with_desc (length, key);
  GError *error       = NULL;
  gchar *message      = g_strdup_printf ("ncm_model_set_reparam: reparametrization `NcHICosmoDEReparamOk' "
                                         "is for %u parameters but model `NcHICosmoDEXcdm' has 7, and %s.",
                                         length, reason);

  ncm_model_set_reparam (NCM_MODEL (cosmo), reparam, &error);
  g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_REPARAM_INCOMPATIBLE);
  g_assert_cmpstr (error->message, ==, message);
  g_assert_null (ncm_model_peek_reparam (NCM_MODEL (cosmo)));

  g_free (message);
  g_clear_error (&error);
  ncm_reparam_free (reparam);
  nc_hicosmo_free (cosmo);
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
  NcBBN *bbn             = NC_BBN (nc_bbn_parametrized_new ());
  NcHICosmoDEXcdm *cosmo = nc_hicosmo_de_xcdm_new_full (reion, prim, bbn);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE_XCDM (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);
  g_assert_true (nc_hicosmo_peek_bbn (NC_HICOSMO (cosmo)) == bbn);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);
  nc_bbn_free (bbn);

  /* All submodels are optional -- NULL/NULL/NULL must behave like the plain
   * zero-arg constructor: reion and prim stay unset, bbn gets the default
   * #NcBBNParthenope. */
  cosmo = nc_hicosmo_de_xcdm_new_full (NULL, NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (NC_IS_BBN_PARTHENOPE (nc_hicosmo_peek_bbn (NC_HICOSMO (cosmo))));
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_bbn_construct_fixed_subprocess (void)
{
  NcHICosmoDEXcdm *cosmo = nc_hicosmo_de_xcdm_new_full (NULL, NULL, NULL);
  NcBBN *bbn             = NC_BBN (nc_bbn_parametrized_new ());

  /* bbn is a typed slot: replacing the default after construction must abort. */
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (bbn));

  nc_bbn_free (bbn);
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_bbn_construct_fixed (void)
{
  g_test_trap_subprocess ("/nc/hicosmo_de/bbn/construct_fixed/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*construction-fixed*");
}

void
test_nc_hicosmo_late_attach (void)
{
  NcHICosmoDEXcdm *cosmo = nc_hicosmo_de_xcdm_new_full (NULL, NULL, NULL);
  NcHIReion *reion       = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim         = NC_HIPRIM (nc_hiprim_power_law_new ());

  /* Deprecated 0.27 path: attaching after construction fills an empty slot. */
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));

  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (ncm_model_peek_host (NCM_MODEL (prim)) == NCM_MODEL (cosmo));
  g_assert_true (ncm_model_peek_host (NCM_MODEL (reion)) == NCM_MODEL (cosmo));

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);
}

void
test_nc_hicosmo_lcdm_new_full (void)
{
  NcHIReion *reion     = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim       = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmoLCDM *cosmo = nc_hicosmo_lcdm_new_full (reion, prim, NULL);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_LCDM (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  cosmo = nc_hicosmo_lcdm_new_full (NULL, NULL, NULL);
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
  NcHICosmoDECpl *cosmo = nc_hicosmo_de_cpl_new_full (reion, prim, NULL);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE_CPL (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  cosmo = nc_hicosmo_de_cpl_new_full (NULL, NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  nc_hicosmo_free (NC_HICOSMO (cosmo));

  /* The plain constructor delegates to _new_full (NULL, NULL, NULL). */
  cosmo = nc_hicosmo_de_cpl_new ();
  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_BBN_PARTHENOPE (nc_hicosmo_peek_bbn (NC_HICOSMO (cosmo))));
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_de_jbp_new_full (void)
{
  NcHIReion *reion      = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim        = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmoDEJbp *cosmo = nc_hicosmo_de_jbp_new_full (reion, prim, NULL);

  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_HICOSMO_DE_JBP (cosmo));
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == reion);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == prim);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
  nc_hireion_free (reion);
  nc_hiprim_free (prim);

  cosmo = nc_hicosmo_de_jbp_new_full (NULL, NULL, NULL);
  g_assert_true (cosmo != NULL);
  g_assert_true (nc_hicosmo_peek_reion (NC_HICOSMO (cosmo)) == NULL);
  g_assert_true (nc_hicosmo_peek_prim (NC_HICOSMO (cosmo)) == NULL);
  nc_hicosmo_free (NC_HICOSMO (cosmo));

  /* The plain constructor delegates to _new_full (NULL, NULL, NULL). */
  cosmo = nc_hicosmo_de_jbp_new ();
  g_assert_true (cosmo != NULL);
  g_assert_true (NC_IS_BBN_PARTHENOPE (nc_hicosmo_peek_bbn (NC_HICOSMO (cosmo))));
  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

void
test_nc_hicosmo_de_param_by_name_submodel (void)
{
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim   = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcmModel *model  = NCM_MODEL (cosmo);
  GError *error    = NULL;

  /* Host's own params are unaffected -- resolved on the host, no fallback
   * needed. */
  ncm_assert_cmpdouble_e (ncm_model___getitem__ (model, "H0", &error), ==,
                          ncm_model_orig_param_get (model, NC_HICOSMO_DE_H0), 1.0e-15, 0.0);
  g_assert_no_error (error);

  /* Bare name falls back into the (unique) submodel that has it. */
  g_assert_cmpfloat (ncm_model___getitem__ (model, "z_re", &error), ==, 13.0);
  g_assert_no_error (error);

  ncm_model___setitem__ (model, "z_re", 11.0, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (ncm_model_orig_param_get (NCM_MODEL (reion), NC_HIREION_CAMB_HII_HEII_Z), ==, 11.0);
  g_assert_cmpfloat (ncm_model___getitem__ (model, "z_re", &error), ==, 11.0);
  g_assert_no_error (error);

  g_assert_cmpfloat (ncm_model___getitem__ (model, "n_SA", &error), ==,
                     ncm_model_orig_param_get (NCM_MODEL (prim), NC_HIPRIM_POWER_LAW_N_SA));
  g_assert_no_error (error);

  /* Qualified "slot:param" form always works and matches the bare form. */
  g_assert_cmpfloat (ncm_model___getitem__ (model, "reion:z_re", &error), ==, 11.0);
  g_assert_no_error (error);

  ncm_model___setitem__ (model, "reion:z_re", 8.5, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (ncm_model___getitem__ (model, "z_re", &error), ==, 8.5);
  g_assert_no_error (error);

  g_assert_cmpfloat (ncm_model___getitem__ (model, "prim:n_SA", &error), ==,
                     ncm_model_orig_param_get (NCM_MODEL (prim), NC_HIPRIM_POWER_LAW_N_SA));
  g_assert_no_error (error);

  /* Unknown parameter name (host and submodels): clear NOT_FOUND error. */
  {
    gdouble val = ncm_model___getitem__ (model, "does_not_exist", &error);

    g_assert_true (isnan (val));
    g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND);
    g_clear_error (&error);
  }

  /* Qualified form with an unknown slot name: clear NOT_FOUND error. */
  {
    gdouble val = ncm_model___getitem__ (model, "nosuchslot:something", &error);

    g_assert_true (isnan (val));
    g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND);
    g_clear_error (&error);
  }

  /* Same behavior through the underlying by-name/desc C API, not just the
   * dict dunders. */
  g_assert_cmpfloat (ncm_model_param_get_by_name (model, "z_re", &error), ==, 8.5);
  g_assert_no_error (error);

  ncm_model_param_set_by_name (model, "z_re", 9.5, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (ncm_model_param_get_by_name (model, "z_re", &error), ==, 9.5);
  g_assert_no_error (error);

  {
    GHashTable *desc = ncm_model_param_get_desc (model, "z_re", &error);

    g_assert_no_error (error);
    g_assert_true (desc != NULL);
    g_hash_table_unref (desc);
  }

  /* Also settable via the qualified form. */
  {
    GHashTable *set_desc = g_hash_table_new (g_str_hash, g_str_equal);
    GValue *value        = g_new0 (GValue, 1);

    g_value_init (value, G_TYPE_BOOLEAN);
    g_value_set_boolean (value, TRUE);
    g_hash_table_insert (set_desc, "fit", value);

    ncm_model_param_set_desc (model, "reion:z_re", set_desc, &error);
    g_assert_no_error (error);

    g_hash_table_unref (set_desc);
    g_value_unset (value);
    g_free (value);
  }

  nc_hicosmo_free (cosmo);
  nc_hireion_free (reion);
  nc_hiprim_free (prim);
}

/*
 * NcHICosmoDE had Yp at index 4 (8 parameters); a reparam of length 8 was written before
 * the removal: its keys above 4 move down by one. The CMB case is the 2024 Planck one,
 * omegac at 1 and omegab at 6, which must land on 1 and 5 and leave w at 6.
 */
void
test_nc_hicosmo_de_reparam_old_layout (void)
{
  {
    NcHICosmo *cosmo    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
    NcmReparam *reparam = _test_nc_hicosmo_de_reparam_with_desc (8, 7);
    GError *error       = NULL;

    ncm_model_set_reparam (NCM_MODEL (cosmo), reparam, &error);
    g_assert_no_error (error);
    g_assert_cmpuint (ncm_model_len (NCM_MODEL (cosmo)), ==, 7);
    g_assert_cmpuint (ncm_reparam_get_length (ncm_model_peek_reparam (NCM_MODEL (cosmo))), ==, 7);
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 6), ==, "Bogus");
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 5), ==, "Omegab");
    /* The stored params-desc replaced the constructor's description of Omegak. */
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 2), ==, "Omegax");

    ncm_reparam_free (reparam);
    nc_hicosmo_free (cosmo);
  }

  {
    NcHICosmo *cosmo    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
    NcmReparam *reparam = NCM_REPARAM (nc_hicosmo_de_reparam_ok_new (8));
    GError *error       = NULL;

    g_object_set (reparam, "params-desc", NULL, NULL);
    ncm_model_set_reparam (NCM_MODEL (cosmo), reparam, &error);
    g_assert_no_error (error);
    g_assert_cmpuint (ncm_model_len (NCM_MODEL (cosmo)), ==, 7);
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 2), ==, "Omegax");

    ncm_reparam_free (reparam);
    nc_hicosmo_free (cosmo);
  }

  {
    NcHICosmo *cosmo    = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
    NcmReparam *reparam = NCM_REPARAM (nc_hicosmo_de_reparam_cmb_new (8));
    NcmObjDictInt *desc = ncm_obj_dict_int_new ();
    NcmSParam *sp_c     = ncm_sparam_new ("omegac", "\\omega_{c0}", 0.0075, 0.25, 0.01, 0.0, 0.12, NCM_PARAM_TYPE_FREE);
    NcmSParam *sp_b     = ncm_sparam_new ("omegab", "\\omega_{b0}", 0.005, 0.04, 0.01, 0.0, 0.022, NCM_PARAM_TYPE_FREE);
    GError *error       = NULL;

    ncm_obj_dict_int_add (desc, 1, G_OBJECT (sp_c));
    ncm_obj_dict_int_add (desc, 6, G_OBJECT (sp_b));
    g_object_set (reparam, "params-desc", desc, NULL);

    ncm_model_set_reparam (NCM_MODEL (cosmo), reparam, &error);
    g_assert_no_error (error);
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 1), ==, "omegac");
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 5), ==, "omegab");
    g_assert_cmpfloat (ncm_model_param_get_upper_bound (NCM_MODEL (cosmo), 5), ==, 0.04);
    g_assert_cmpstr (ncm_model_param_name (NCM_MODEL (cosmo), 6), ==, "w");

    ncm_sparam_free (sp_c);
    ncm_sparam_free (sp_b);
    ncm_obj_dict_int_unref (desc);
    ncm_reparam_free (reparam);
    nc_hicosmo_free (cosmo);
  }
}

/* The keys the resize refuses: Yp itself, a negative key, and both failures of a length
 * that matches neither layout. */
void
test_nc_hicosmo_de_reparam_old_layout_refused (void)
{
  _test_nc_hicosmo_de_reparam_refused (8, 4, "it reparametrizes parameter 4, which was removed from the model");
  _test_nc_hicosmo_de_reparam_refused (8, -1, "it reparametrizes parameter -1, which no longer exists");
  _test_nc_hicosmo_de_reparam_refused (9, 8, "it reparametrizes parameter 8, which no longer exists");
  _test_nc_hicosmo_de_reparam_refused (9, 5, "it reparametrizes parameter 5, at or after the removed parameter 4, "
                                             "and its length is not the 8 the model had before the removal");
}


/***************************************************************************
 *            test_ncm_object_serialization.c
 *
 *  Wed June 20 21:56:09 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <vitenti@uel.br> <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

typedef struct _TestNcmSerialize
{
  NcmSerialize *ser;
} TestNcmSerialize;

static void test_ncm_serialize_new (TestNcmSerialize *test, gconstpointer pdata);
void test_ncm_serialize_new_noclean_dup (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_free (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_global_from_string_plain (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_from_string_params (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_from_string_nest (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_from_string_plain_named (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_string_params_named (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_string_nest_named (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_from_string_nest_samename (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_to_file_from_file (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_to_binfile_from_binfile (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_reset_autosave_only (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_traps (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_invalid_from_string_syntax (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_global_invalid_from_string_nonexist (TestNcmSerialize *test, gconstpointer pdata);

static void test_ncm_serialize_invalid_from_string_mnames (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_invalid_from_string_samename (TestNcmSerialize *test, gconstpointer pdata);
static void test_ncm_serialize_invalid_from_string_wrongorder (TestNcmSerialize *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/ncm/serialize/global/from_string/plain", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_from_string_plain,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/global/from_string/params", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_from_string_params,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/global/from_string/nest", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_from_string_nest,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/plain/named", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_plain_named,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/params/named", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_params_named,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/nested/named", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_nest_named,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/nested/samename", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_from_string_nest_samename,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/from_string/plain/reset/autosave_only", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_reset_autosave_only,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/to_file/from_file", TestNcmSerialize, NULL,
              &test_ncm_serialize_new_noclean_dup,
              &test_ncm_serialize_to_file_from_file,
              &test_ncm_serialize_free);
  g_test_add ("/ncm/serialize/to_binfile/from_binfile", TestNcmSerialize, NULL,
              &test_ncm_serialize_new_noclean_dup,
              &test_ncm_serialize_to_binfile_from_binfile,
              &test_ncm_serialize_free);
  
  g_test_add ("/ncm/serialize/traps", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_traps,
              &test_ncm_serialize_free);
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_add ("/ncm/serialize/invalid/from_string/syntax/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_invalid_from_string_syntax,
              &test_ncm_serialize_free);
  
  g_test_add ("/ncm/serialize/invalid/from_string/nonexistent/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_global_invalid_from_string_nonexist,
              &test_ncm_serialize_free);
  
  g_test_add ("/ncm/serialize/invalid/from_string/multiplenames/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_invalid_from_string_mnames,
              &test_ncm_serialize_free);
  
  g_test_add ("/ncm/serialize/invalid/from_string/samename/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_invalid_from_string_samename,
              &test_ncm_serialize_free);
  
  g_test_add ("/ncm/serialize/invalid/from_string/wrongorder/subprocess", TestNcmSerialize, NULL,
              &test_ncm_serialize_new,
              &test_ncm_serialize_invalid_from_string_wrongorder,
              &test_ncm_serialize_free);
#endif
  
  g_test_run ();
}

void
test_ncm_serialize_new (TestNcmSerialize *test, gconstpointer pdata)
{
  test->ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
}

void
test_ncm_serialize_new_noclean_dup (TestNcmSerialize *test, gconstpointer pdata)
{
  test->ser = ncm_serialize_new (0);
}

void
test_ncm_serialize_free (TestNcmSerialize *test, gconstpointer pdata)
{
  NcmSerialize *ser = test->ser;
  
  NCM_TEST_FREE (ncm_serialize_free, ser);
}

void
test_ncm_serialize_global_from_string_plain (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj       = ncm_serialize_global_from_string ("NcHICosmoLCDM");
  NcHICosmo *hic     = NC_HICOSMO (obj);
  gchar *obj_ser     = ncm_serialize_global_to_string (obj, TRUE);
  GObject *obj_new   = ncm_serialize_global_from_string (obj_ser);
  gchar *obj_new_ser = ncm_serialize_global_to_string (obj_new, TRUE);
  
  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  
  obj_ser     = ncm_serialize_global_to_string (obj, FALSE);
  obj_new_ser = ncm_serialize_global_to_string (obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);
  
  NCM_TEST_FREE (nc_hicosmo_free, hic);
}

void
test_ncm_serialize_global_from_string_params (TestNcmSerialize *test, gconstpointer pdata)
{
  NcmModel *m = NCM_MODEL (ncm_serialize_global_from_string ("NcHICosmoLCDM{'H0':<12.3>,'Omegac':<0.2>}"));
  
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);
  
  NCM_TEST_FREE (ncm_model_free, m);
  
  m = NCM_MODEL (ncm_serialize_global_from_string ("{'NcHICosmoLCDM', {'H0':<12.3>,'Omegac':<0.2>}}"));
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);
  
  NCM_TEST_FREE (ncm_model_free, m);
}

void
test_ncm_serialize_global_from_string_nest (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_global_from_string ("NcPowspecMLTransfer{'transfer':<{'NcTransferFuncEH',@a{sv} {}}>}");
  
  NcPowspecML *ps_ml = NC_POWSPEC_ML (obj);
  
  gchar *obj_ser     = ncm_serialize_global_to_string (obj, TRUE);
  GObject *obj_new   = ncm_serialize_global_from_string (obj_ser);
  gchar *obj_new_ser = ncm_serialize_global_to_string (obj_new, TRUE);
  
  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  
  g_free (obj_ser);
  g_free (obj_new_ser);
  
  obj_ser     = ncm_serialize_global_to_string (obj, FALSE);
  obj_new_ser = ncm_serialize_global_to_string (obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);
  
  g_assert_true (NC_IS_TRANSFER_FUNC_EH (NC_POWSPEC_ML_TRANSFER (ps_ml)->tf));
  
  NCM_TEST_FREE (nc_powspec_ml_free, ps_ml);
}

static void
test_ncm_serialize_from_string_plain_named (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj       = ncm_serialize_from_string (test->ser, "NcHICosmoLCDM[myname:is]");
  NcHICosmo *hic     = NC_HICOSMO (obj);
  gchar *obj_ser     = ncm_serialize_to_string (test->ser, obj, TRUE);
  GObject *obj_new   = ncm_serialize_from_string (test->ser, obj_ser);
  gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);
  
  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  
  obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);
  
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);
  
  ncm_serialize_clear_instances (test->ser, FALSE);
  
  NCM_TEST_FREE (nc_hicosmo_free, hic);
}

static void
test_ncm_serialize_reset_autosave_only (TestNcmSerialize *test, gconstpointer pdata)
{
  /* Deserialized from string, since it is named it will be automatically added to the saved instances. */
  GObject *obj1 = ncm_serialize_from_string (test->ser, "NcHICosmoLCDM[myname:is]");
  /* Deserialized from string, since it is *not* named it will *not* be automatically added to the saved instances.*/
  GObject *obj2 = ncm_serialize_from_string (test->ser, "NcPowspecMLTransfer{'transfer':<{'NcTransferFuncEH',@a{sv} {}}>}");
  
  /* Serializing from object */
  gchar *obj_ser1 = ncm_serialize_to_string (test->ser, obj1, FALSE);
  /* Serializing from object, all objects will receive a name */
  gchar *obj_ser2 = ncm_serialize_to_string (test->ser, obj2, FALSE);
  
  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), ==, 1);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), >=, 2);
  
  {
    GObject *obj1_ref2 = ncm_serialize_from_string (test->ser, obj_ser1);
    GObject *obj2_ref2 = ncm_serialize_from_string (test->ser, obj_ser2);
    
    g_assert_true (obj1_ref2 == obj1);
    g_assert_true (obj2_ref2 != obj2);
    
    g_object_unref (obj1_ref2);
    g_object_unref (obj2_ref2);
    
    g_free (obj_ser1);
    g_free (obj_ser2);
  }
  
  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), >, 1);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), >=, 2);
  
  ncm_serialize_reset (test->ser, TRUE);
  
  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), ==, 1);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), ==, 0);
  
  ncm_serialize_reset (test->ser, FALSE);
  
  g_assert_cmpuint (ncm_serialize_count_instances (test->ser), ==, 0);
  g_assert_cmpuint (ncm_serialize_count_saved_serializations (test->ser), ==, 0);
  
  NCM_TEST_FREE (g_object_unref, obj1);
  NCM_TEST_FREE (g_object_unref, obj2);
}

void
test_ncm_serialize_from_string_params_named (TestNcmSerialize *test, gconstpointer pdata)
{
  NcmModel *m = NCM_MODEL (ncm_serialize_from_string (test->ser, "NcHICosmoLCDM[T0]{'H0':<12.3>,'Omegac':<0.2>}"));
  
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);
  
  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (ncm_model_free, m);
  
  m = NCM_MODEL (ncm_serialize_from_string (test->ser, "{'NcHICosmoLCDM[T1]', {'H0':<12.3>,'Omegac':<0.2>}}"));
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);
  
  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (ncm_model_free, m);
}

void
test_ncm_serialize_from_string_nest_named (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser, "NcPowspecMLTransfer[T0]{'transfer':<{'NcTransferFuncEH[T1]',@a{sv} {}}>}");
  
  NcPowspecML *ps_ml = NC_POWSPEC_ML (obj);
  
  gchar *obj_ser     = ncm_serialize_to_string (test->ser, obj, TRUE);
  GObject *obj_new   = ncm_serialize_from_string (test->ser, obj_ser);
  gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);
  
  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  
  g_free (obj_ser);
  g_free (obj_new_ser);
  
  obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);
  
  g_assert_true (NC_IS_TRANSFER_FUNC_EH (NC_POWSPEC_ML_TRANSFER (ps_ml)->tf));
  
  g_assert_cmphex (GPOINTER_TO_SIZE (ncm_serialize_peek_by_name (test->ser, "T0")), ==, GPOINTER_TO_SIZE (ps_ml));
  g_assert_cmphex (GPOINTER_TO_SIZE (ncm_serialize_peek_by_name (test->ser, "T1")), ==, GPOINTER_TO_SIZE (NC_POWSPEC_ML_TRANSFER (ps_ml)->tf));
  
  g_assert_cmpstr (ncm_serialize_peek_name (test->ser, NC_POWSPEC_ML_TRANSFER (ps_ml)->tf), ==, "T1");
  g_assert_cmpstr (ncm_serialize_peek_name (test->ser, ps_ml), ==, "T0");
  
  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (nc_powspec_ml_free, ps_ml);
}

void
test_ncm_serialize_from_string_nest_samename (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser,
                                            "NcmSplineCubicNotaknot{'length':<6>, "
                                            "'x' : <{'NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>}}>, "
                                            "'y' : <{'NcmVector[T0]', @a{sv} {}}>}");
  
  NcmSpline *s = NCM_SPLINE (obj);
  gchar *obj_ser, *obj_new_ser;
  GObject *obj_new;
  
  ncm_serialize_clear_instances (test->ser, FALSE);
  
  obj_ser = ncm_serialize_to_string (test->ser, obj, TRUE);
  obj_new = ncm_serialize_from_string (test->ser, obj_ser);
  
  ncm_serialize_reset (test->ser, FALSE);
  
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);
  
  g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  
  g_free (obj_ser);
  g_free (obj_new_ser);
  
  ncm_serialize_reset (test->ser, FALSE);
  obj_ser = ncm_serialize_to_string (test->ser, obj, FALSE);
  
  ncm_serialize_reset (test->ser, FALSE);
  obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);
  
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);
  
  ncm_serialize_clear_instances (test->ser, FALSE);
  NCM_TEST_FREE (ncm_spline_free, s);
}

static void
test_ncm_serialize_to_file_from_file (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj   = ncm_serialize_from_string (test->ser, "NcHICosmoDEXcdm{'w':<-2.0>}");
  NcHICosmo *hic = NC_HICOSMO (obj);
  gchar *obj_ser = ncm_serialize_to_string (test->ser, obj, TRUE);
  
  ncm_serialize_to_file (test->ser, obj, "test-serialize-file.obj");
  
  {
    GObject *obj_new   = ncm_serialize_from_file (test->ser, "test-serialize-file.obj");
    gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);
    
    g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
    g_assert_cmpstr (obj_ser, ==, obj_new_ser);
    g_free (obj_ser);
    g_free (obj_new_ser);
    
    obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
    obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);
    
    g_assert_cmpstr (obj_ser, ==, obj_new_ser);
    
    g_free (obj_ser);
    g_free (obj_new_ser);
    g_object_unref (obj_new);
    
    ncm_serialize_clear_instances (test->ser, FALSE);
    
    NCM_TEST_FREE (nc_hicosmo_free, hic);
  }
}

static void
test_ncm_serialize_to_binfile_from_binfile (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj   = ncm_serialize_from_string (test->ser, "NcHICosmoDEXcdm{'w':<-2.0>}");
  NcHICosmo *hic = NC_HICOSMO (obj);
  gchar *obj_ser = ncm_serialize_to_string (test->ser, obj, TRUE);
  
  ncm_serialize_to_binfile (test->ser, obj, "test-serialize-binfile.obj");
  
  {
    GObject *obj_new   = ncm_serialize_from_binfile (test->ser, "test-serialize-binfile.obj");
    gchar *obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, TRUE);
    
    g_assert_true (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
    g_assert_cmpstr (obj_ser, ==, obj_new_ser);
    g_free (obj_ser);
    g_free (obj_new_ser);
    
    obj_ser     = ncm_serialize_to_string (test->ser, obj, FALSE);
    obj_new_ser = ncm_serialize_to_string (test->ser, obj_new, FALSE);
    
    g_assert_cmpstr (obj_ser, ==, obj_new_ser);
    
    g_free (obj_ser);
    g_free (obj_new_ser);
    g_object_unref (obj_new);
    
    ncm_serialize_clear_instances (test->ser, FALSE);
    
    NCM_TEST_FREE (nc_hicosmo_free, hic);
  }
}

void
test_ncm_serialize_traps (TestNcmSerialize *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION (2, 38, 0)
  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/syntax/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/nonexistent/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/multiplenames/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/samename/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  
  g_test_trap_subprocess ("/ncm/serialize/invalid/from_string/wrongorder/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_serialize_global_invalid_from_string_syntax (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_global_from_string ("NcPowspecMLTransfer Error{'transfer':<{'NcTransferFuncEH',@a{sv} {}}>}");
  
  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_global_invalid_from_string_nonexist (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_global_from_string ("NcPowspecMLErrorName   {'transfer':<{'NcTransferFuncEH',@a{sv} {}}>}");
  
  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_invalid_from_string_mnames (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser, "NcPowspecMLTransfer[T0]{'transfer':<{'NcTransferFuncEH[T0]', @a{sv} {}}>}");
  
  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_invalid_from_string_samename (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser,
                                            "NcmSplineCubicNotaknot{'length':<6>, "
                                            "'x' : <{'NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>}}>, "
                                            "'y' : <{'NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>}}>}");
  
  NCM_TEST_FREE (g_object_unref, obj);
}

void
test_ncm_serialize_invalid_from_string_wrongorder (TestNcmSerialize *test, gconstpointer pdata)
{
  GObject *obj = ncm_serialize_from_string (test->ser,
                                            "NcmSplineCubicNotaknot{'length':<6>, "
                                            "'x' : <{'NcmVector[T0]', @a{sv} {}}>, "
                                            "'y' : <{'NcmVector[T0]', {'values':<[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]>}}>}");
  
  NCM_TEST_FREE (g_object_unref, obj);
}


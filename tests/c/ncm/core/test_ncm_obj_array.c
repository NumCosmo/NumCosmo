/***************************************************************************
 *            ncm_obj_array.c
 *
 *  Fri March 23 14:26:42 2012
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/*
 * NcmObjArray
 */

typedef struct _TestNcmObjArray
{
  NcmObjArray *oa;
  guint ntests;
} TestNcmObjArray;

void test_ncm_obj_array_new (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_free (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_add (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_to_from_keyfile (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_to_from_yaml (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_to_from_yaml_file (TestNcmObjArray *test, gconstpointer pdata);

void test_ncm_obj_array_traps (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_invalid_add (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_invalid_set (TestNcmObjArray *test, gconstpointer pdata);

/*
 * NcmObjDictStr
 */

typedef struct _TestNcmObjDictStr
{
  NcmObjDictStr *ods;
  guint ntests;
} TestNcmObjDictStr;

void test_ncm_obj_dict_str_new (TestNcmObjDictStr *test, gconstpointer pdata);
void test_ncm_obj_dict_str_free (TestNcmObjDictStr *test, gconstpointer pdata);
void test_ncm_obj_dict_str_add (TestNcmObjDictStr *test, gconstpointer pdata);
void test_ncm_obj_dict_str_to_from_yaml (TestNcmObjDictStr *test, gconstpointer pdata);
void test_ncm_obj_dict_str_to_from_yaml_file (TestNcmObjDictStr *test, gconstpointer pdata);

void test_ncm_obj_dict_str_traps (TestNcmObjDictStr *test, gconstpointer pdata);
void test_ncm_obj_dict_str_invalid_add (TestNcmObjDictStr *test, gconstpointer pdata);
void test_ncm_obj_dict_str_invalid_set (TestNcmObjDictStr *test, gconstpointer pdata);


/*
 * NcmObjDictInt
 */

typedef struct _TestNcmObjDictInt
{
  NcmObjDictInt *odi;
  guint ntests;
} TestNcmObjDictInt;

void test_ncm_obj_dict_int_new (TestNcmObjDictInt *test, gconstpointer pdata);
void test_ncm_obj_dict_int_free (TestNcmObjDictInt *test, gconstpointer pdata);
void test_ncm_obj_dict_int_add (TestNcmObjDictInt *test, gconstpointer pdata);
void test_ncm_obj_dict_int_to_from_yaml (TestNcmObjDictInt *test, gconstpointer pdata);
void test_ncm_obj_dict_int_to_from_yaml_file (TestNcmObjDictInt *test, gconstpointer pdata);

void test_ncm_obj_dict_int_traps (TestNcmObjDictInt *test, gconstpointer pdata);
void test_ncm_obj_dict_int_invalid_add (TestNcmObjDictInt *test, gconstpointer pdata);
void test_ncm_obj_dict_int_invalid_set (TestNcmObjDictInt *test, gconstpointer pdata);


/*
 * NcmVarDict:
 */

typedef struct _TestNcmVarDict
{
  NcmVarDict *vd;
  guint ntests;
} TestNcmVarDict;

void test_ncm_var_dict_new (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_free (TestNcmVarDict *test, gconstpointer pdata);

void test_ncm_var_dict_has_key (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_string (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_int (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_double (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_boolean (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_int_array (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_double_array (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_boolean_array (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_variant (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_to_from (TestNcmVarDict *test, gconstpointer pdata);

NcmVarDict *_test_ncm_var_dict_to_from_variant (NcmVarDict *vd);
NcmVarDict *_test_ncm_var_dict_to_from_yaml (NcmVarDict *vd);
NcmVarDict *_test_ncm_var_dict_to_from_variant_binfile (NcmVarDict *vd);
NcmVarDict *_test_ncm_var_dict_to_from_variant_file (NcmVarDict *vd);
NcmVarDict *_test_ncm_var_dict_to_from_yaml_file (NcmVarDict *vd);

void test_ncm_var_dict_traps (TestNcmVarDict *test, gconstpointer pdata);
void test_ncm_var_dict_invalid_set (TestNcmVarDict *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /*
   * NcmObjArray
   */

  g_test_add ("/ncm/obj_array/add", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_add,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/to_from/keyfile", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_to_from_keyfile,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/to_from/yaml", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_to_from_yaml,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/to_from/yaml_file", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_to_from_yaml_file,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/traps", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_traps,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/invalid/add/subprocess", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_invalid_add,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/invalid/set/subprocess", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_invalid_set,
              &test_ncm_obj_array_free);

  /*
   * NcmObjDictStr
   */

  g_test_add ("/ncm/obj_dict_str/add", TestNcmObjDictStr, NULL,
              &test_ncm_obj_dict_str_new,
              &test_ncm_obj_dict_str_add,
              &test_ncm_obj_dict_str_free);

  g_test_add ("/ncm/obj_dict_str/to_from/yaml", TestNcmObjDictStr, NULL,
              &test_ncm_obj_dict_str_new,
              &test_ncm_obj_dict_str_to_from_yaml,
              &test_ncm_obj_dict_str_free);

  g_test_add ("/ncm/obj_dict_str/to_from/yaml_file", TestNcmObjDictStr, NULL,
              &test_ncm_obj_dict_str_new,
              &test_ncm_obj_dict_str_to_from_yaml_file,
              &test_ncm_obj_dict_str_free);

  g_test_add ("/ncm/obj_dict_str/traps", TestNcmObjDictStr, NULL,
              &test_ncm_obj_dict_str_new,
              &test_ncm_obj_dict_str_traps,
              &test_ncm_obj_dict_str_free);

  g_test_add ("/ncm/obj_dict_str/invalid/add/subprocess", TestNcmObjDictStr, NULL,
              &test_ncm_obj_dict_str_new,
              &test_ncm_obj_dict_str_invalid_add,
              &test_ncm_obj_dict_str_free);

  g_test_add ("/ncm/obj_dict_str/invalid/set/subprocess", TestNcmObjDictStr, NULL,
              &test_ncm_obj_dict_str_new,
              &test_ncm_obj_dict_str_invalid_set,
              &test_ncm_obj_dict_str_free);

  /*
   * NcmObjDictInt
   */

  g_test_add ("/ncm/obj_dict_int/add", TestNcmObjDictInt, NULL,
              &test_ncm_obj_dict_int_new,
              &test_ncm_obj_dict_int_add,
              &test_ncm_obj_dict_int_free);

  g_test_add ("/ncm/obj_dict_int/to_from/yaml", TestNcmObjDictInt, NULL,
              &test_ncm_obj_dict_int_new,
              &test_ncm_obj_dict_int_to_from_yaml,
              &test_ncm_obj_dict_int_free);

  g_test_add ("/ncm/obj_dict_int/to_from/yaml_file", TestNcmObjDictInt, NULL,
              &test_ncm_obj_dict_int_new,
              &test_ncm_obj_dict_int_to_from_yaml_file,
              &test_ncm_obj_dict_int_free);

  g_test_add ("/ncm/obj_dict_int/traps", TestNcmObjDictInt, NULL,
              &test_ncm_obj_dict_int_new,
              &test_ncm_obj_dict_int_traps,
              &test_ncm_obj_dict_int_free);

  g_test_add ("/ncm/obj_dict_int/invalid/add/subprocess", TestNcmObjDictInt, NULL,
              &test_ncm_obj_dict_int_new,
              &test_ncm_obj_dict_int_invalid_add,
              &test_ncm_obj_dict_int_free);

  g_test_add ("/ncm/obj_dict_int/invalid/set/subprocess", TestNcmObjDictInt, NULL,
              &test_ncm_obj_dict_int_new,
              &test_ncm_obj_dict_int_invalid_set,
              &test_ncm_obj_dict_int_free);

  /*
   * NcmVarDict
   */

  g_test_add ("/ncm/var_dict/has_key", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_has_key,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/string", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_string,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/int", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_int,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/double", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_double,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/bool", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_boolean,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/int_array", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_int_array,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/double_array", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_double_array,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/boolean_array", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_boolean_array,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/variant", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_variant,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/to_from/variant", TestNcmVarDict, &_test_ncm_var_dict_to_from_variant,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_to_from,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/to_from/yaml", TestNcmVarDict, &_test_ncm_var_dict_to_from_yaml,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_to_from,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/to_from/variant_binfile", TestNcmVarDict, &_test_ncm_var_dict_to_from_variant_binfile,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_to_from,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/to_from/variant_file", TestNcmVarDict, &_test_ncm_var_dict_to_from_variant_file,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_to_from,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/to_from/yaml_file", TestNcmVarDict, &_test_ncm_var_dict_to_from_yaml_file,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_to_from,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/traps", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_traps,
              &test_ncm_var_dict_free);

  g_test_add ("/ncm/var_dict/invalid/set/subprocess", TestNcmVarDict, NULL,
              &test_ncm_var_dict_new,
              &test_ncm_var_dict_invalid_set,
              &test_ncm_var_dict_free);

  g_test_run ();
}

/*
 * NcmObjArray
 */

void
test_ncm_obj_array_new (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa;

  test->ntests = 1000;
  oa           = test->oa = ncm_obj_array_new ();

  g_assert_true (oa != NULL);
}

void
test_ncm_obj_array_free (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;

  ncm_obj_array_unref (oa);
}

void
test_ncm_obj_array_add (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;

  g_assert_true (oa != NULL);

  {
    NcmVector *v     = ncm_vector_new (10);
    NcmMatrix *m     = ncm_matrix_new (5, 8);
    NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_lcdm_new ());
    NcDistance *dist = nc_distance_new (1.0);
    NcDataBaoA *baoa = nc_data_bao_a_new_from_id (dist, 0);

    ncm_vector_set_all (v, 1.2);
    ncm_matrix_set_all (m, 2.0);

    ncm_obj_array_add (oa, G_OBJECT (v));
    ncm_obj_array_add (oa, G_OBJECT (m));
    ncm_obj_array_add (oa, G_OBJECT (cosmo));

    ncm_obj_array_set (oa, 0, G_OBJECT (dist));
    ncm_obj_array_add (oa, G_OBJECT (dist));

    ncm_obj_array_add (oa, G_OBJECT (baoa));

    ncm_vector_free (v);
    ncm_matrix_free (m);
    nc_hicosmo_free (cosmo);
    nc_distance_free (dist);
    g_object_unref (baoa);
  }
}

void
test_ncm_obj_array_to_from_keyfile (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa   = test->oa;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_obj_array_saved_XXXXXX", NULL);
  gchar *filename   = g_strdup_printf ("%s/test_ncm_obj_array_saved.oa", tmp_dir);

  test_ncm_obj_array_add (test, pdata);

  ncm_serialize_array_to_key_file (ser, oa, "test_ncm_obj_array_saved.oa", TRUE);

  {
    NcmObjArray *oa_load = ncm_serialize_array_from_key_file (ser, "test_ncm_obj_array_saved.oa");
    guint i;

    g_assert_cmpuint (ncm_obj_array_len (oa), ==, ncm_obj_array_len (oa_load));

    for (i = 0; i < ncm_obj_array_len (oa); i++)
    {
      g_assert_true (ncm_obj_array_peek (oa, i) != NULL);
      g_assert_true (ncm_obj_array_peek (oa_load, i) != NULL);

      g_assert_true (G_IS_OBJECT (ncm_obj_array_peek (oa, i)));
      g_assert_true (G_IS_OBJECT (ncm_obj_array_peek (oa_load, i)));

      g_assert_true (ncm_obj_array_peek (oa, i) != ncm_obj_array_peek (oa_load, i));
    }

    ncm_obj_array_unref (oa_load);
  }

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_array_to_from_yaml (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa   = test->oa;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *yaml_str;

  test_ncm_obj_array_add (test, pdata);

  yaml_str = ncm_serialize_array_to_yaml (ser, oa);

  {
    NcmObjArray *oa_load = ncm_serialize_array_from_yaml (ser, yaml_str);
    guint i;

    g_assert_cmpuint (ncm_obj_array_len (oa), ==, ncm_obj_array_len (oa_load));

    for (i = 0; i < ncm_obj_array_len (oa); i++)
    {
      g_assert_true (ncm_obj_array_peek (oa, i) != NULL);
      g_assert_true (ncm_obj_array_peek (oa_load, i) != NULL);

      g_assert_true (G_IS_OBJECT (ncm_obj_array_peek (oa, i)));
      g_assert_true (G_IS_OBJECT (ncm_obj_array_peek (oa_load, i)));

      g_assert_true (ncm_obj_array_peek (oa, i) != ncm_obj_array_peek (oa_load, i));
    }

    ncm_obj_array_unref (oa_load);
  }

  g_free (yaml_str);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_array_to_from_yaml_file (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa   = test->oa;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_obj_array_saved_XXXXXX", NULL);
  gchar *filename   = g_strdup_printf ("%s/test_ncm_obj_array_saved.yaml", tmp_dir);

  test_ncm_obj_array_add (test, pdata);

  ncm_serialize_array_to_yaml_file (ser, oa, filename);

  {
    NcmObjArray *oa_load = ncm_serialize_array_from_yaml_file (ser, filename);
    guint i;

    g_assert_cmpuint (ncm_obj_array_len (oa), ==, ncm_obj_array_len (oa_load));

    for (i = 0; i < ncm_obj_array_len (oa); i++)
    {
      g_assert_true (ncm_obj_array_peek (oa, i) != NULL);
      g_assert_true (ncm_obj_array_peek (oa_load, i) != NULL);

      g_assert_true (G_IS_OBJECT (ncm_obj_array_peek (oa, i)));
      g_assert_true (G_IS_OBJECT (ncm_obj_array_peek (oa_load, i)));

      g_assert_true (ncm_obj_array_peek (oa, i) != ncm_obj_array_peek (oa_load, i));
    }

    ncm_obj_array_unref (oa_load);
  }

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_array_traps (TestNcmObjArray *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/obj_array/invalid/add/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/obj_array/invalid/set/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_obj_array_invalid_add (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;

  ncm_obj_array_add (oa, NULL);
}

void
test_ncm_obj_array_invalid_set (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;

  ncm_obj_array_add (oa, NULL);

  test_ncm_obj_array_add (test, pdata);

  {
    NcmVector *v = ncm_vector_new (10);

    ncm_obj_array_set (oa, 10, G_OBJECT (v));

    ncm_vector_free (v);
  }
}

/*
 * NcmObjDictStr
 */

void
test_ncm_obj_dict_str_new (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods;

  test->ntests = 1000;
  ods          = test->ods = ncm_obj_dict_str_new ();

  g_assert_true (ods != NULL);
}

void
test_ncm_obj_dict_str_free (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods = test->ods;

  ncm_obj_dict_str_unref (ods);
}

void
test_ncm_obj_dict_str_add (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods = test->ods;

  g_assert_true (ods != NULL);

  {
    NcmVector *v     = ncm_vector_new (10);
    NcmMatrix *m     = ncm_matrix_new (5, 8);
    NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
    NcDistance *dist = nc_distance_new (1.0);
    NcDataBaoA *baoa = nc_data_bao_a_new_from_id (dist, 0);

    /* This sets an object containing a DictInt into cosmo */
    nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);

    ncm_vector_set_all (v, 1.2);
    ncm_matrix_set_all (m, 2.0);

    ncm_obj_dict_str_add (ods, "v", G_OBJECT (v));
    ncm_obj_dict_str_add (ods, "m", G_OBJECT (m));
    ncm_obj_dict_str_add (ods, "cosmo", G_OBJECT (cosmo));

    ncm_obj_dict_str_set (ods, "dist", G_OBJECT (dist));
    ncm_obj_dict_str_add (ods, "dist", G_OBJECT (dist));

    ncm_obj_dict_str_add (ods, "baoa", G_OBJECT (baoa));

    g_assert_null (ncm_obj_dict_str_get (ods, "IdontExist"));
    g_assert_null (ncm_obj_dict_str_peek (ods, "IdontExist"));

    {
      GObject *obj = ncm_obj_dict_str_get (ods, "v");

      g_assert_nonnull (obj);
      g_assert_true (G_IS_OBJECT (obj));
      g_assert_true (obj == G_OBJECT (v));

      g_assert_true (obj == ncm_obj_dict_str_peek (ods, "v"));

      g_object_unref (obj);
    }

    g_assert_nonnull (ncm_obj_dict_str_peek (ods, "v"));

    ncm_vector_free (v);
    ncm_matrix_free (m);
    nc_hicosmo_free (cosmo);
    nc_distance_free (dist);
    g_object_unref (baoa);
  }
}

void
test_ncm_obj_dict_str_to_from_yaml (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods = test->ods;
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *yaml_str;

  test_ncm_obj_dict_str_add (test, pdata);

  yaml_str = ncm_serialize_dict_str_to_yaml (ser, ods);

  {
    NcmObjDictStr *ods_load = ncm_serialize_dict_str_from_yaml (ser, yaml_str);
    GHashTableIter iter;
    gchar *key;
    GObject *value;

    g_assert_cmpuint (ncm_obj_dict_str_len (ods), ==, ncm_obj_dict_str_len (ods_load));

    g_hash_table_iter_init (&iter, ods);

    while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
    {
      GObject *obj = ncm_obj_dict_str_get (ods_load, key);

      g_assert_nonnull (obj);
      g_assert_nonnull (ncm_obj_dict_str_peek (ods_load, key));
      g_assert_true (obj == ncm_obj_dict_str_peek (ods_load, key));
      g_assert_true (G_IS_OBJECT (obj));

      g_assert_true (ncm_obj_dict_str_peek (ods, key) != ncm_obj_dict_str_peek (ods_load, key));

      g_object_unref (obj);
    }

    ncm_obj_dict_str_unref (ods_load);
  }

  g_free (yaml_str);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_dict_str_to_from_yaml_file (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods = test->ods;
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *tmp_dir     = g_dir_make_tmp ("test_ncm_obj_dict_str_saved_XXXXXX", NULL);
  gchar *filename    = g_strdup_printf ("%s/test_ncm_obj_dict_str_saved.yaml", tmp_dir);

  test_ncm_obj_dict_str_add (test, pdata);

  ncm_serialize_dict_str_to_yaml_file (ser, ods, filename);

  {
    NcmObjDictStr *ods_load = ncm_serialize_dict_str_from_yaml_file (ser, filename);
    GHashTableIter iter;
    gchar *key;
    GObject *value;

    g_assert_cmpuint (ncm_obj_dict_str_len (ods), ==, ncm_obj_dict_str_len (ods_load));

    g_hash_table_iter_init (&iter, ods);

    while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
    {
      GObject *obj = ncm_obj_dict_str_get (ods_load, key);

      g_assert_nonnull (obj);
      g_assert_nonnull (ncm_obj_dict_str_peek (ods_load, key));
      g_assert_true (obj == ncm_obj_dict_str_peek (ods_load, key));
      g_assert_true (G_IS_OBJECT (obj));

      g_assert_true (ncm_obj_dict_str_peek (ods, key) != ncm_obj_dict_str_peek (ods_load, key));

      g_object_unref (obj);
    }

    ncm_obj_dict_str_unref (ods_load);
  }

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_dict_str_traps (TestNcmObjDictStr *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/obj_dict_str/invalid/add/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/obj_dict_str/invalid/set/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_obj_dict_str_invalid_add (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods = test->ods;

  ncm_obj_dict_str_add (ods, NULL, NULL);
}

void
test_ncm_obj_dict_str_invalid_set (TestNcmObjDictStr *test, gconstpointer pdata)
{
  NcmObjDictStr *ods = test->ods;

  ncm_obj_dict_str_set (ods, NULL, NULL);
}

/*
 * NcmObjDictInt
 */

void
test_ncm_obj_dict_int_new (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi;

  test->ntests = 1000;
  odi          = test->odi = ncm_obj_dict_int_new ();

  g_assert_true (odi != NULL);
}

void
test_ncm_obj_dict_int_free (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi = test->odi;

  ncm_obj_dict_int_unref (odi);
}

void
test_ncm_obj_dict_int_add (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi = test->odi;

  g_assert_true (odi != NULL);

  {
    NcmVector *v     = ncm_vector_new (10);
    NcmMatrix *m     = ncm_matrix_new (5, 8);
    NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
    NcDistance *dist = nc_distance_new (1.0);
    NcDataBaoA *baoa = nc_data_bao_a_new_from_id (dist, 0);

    /* This sets an object containing a DictInt into cosmo */
    nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);

    ncm_vector_set_all (v, 1.2);
    ncm_matrix_set_all (m, 2.0);

    ncm_obj_dict_int_add (odi, 0, G_OBJECT (v));
    ncm_obj_dict_int_add (odi, 123, G_OBJECT (m));
    ncm_obj_dict_int_add (odi, 233, G_OBJECT (cosmo));

    ncm_obj_dict_int_set (odi, 3000, G_OBJECT (dist));
    ncm_obj_dict_int_add (odi, 3000, G_OBJECT (dist));

    ncm_obj_dict_int_add (odi, 66, G_OBJECT (baoa));

    g_assert_null (ncm_obj_dict_int_get (odi, 5));
    g_assert_null (ncm_obj_dict_int_peek (odi, 5));

    {
      GObject *obj = ncm_obj_dict_int_get (odi, 0);

      g_assert_nonnull (obj);
      g_assert_true (G_IS_OBJECT (obj));
      g_assert_true (obj == G_OBJECT (v));

      g_assert_true (obj == ncm_obj_dict_int_peek (odi, 0));

      g_object_unref (obj);
    }

    g_assert_nonnull (ncm_obj_dict_int_peek (odi, 0));

    ncm_vector_free (v);
    ncm_matrix_free (m);
    nc_hicosmo_free (cosmo);
    nc_distance_free (dist);
    g_object_unref (baoa);
  }
}

void
test_ncm_obj_dict_int_to_from_yaml (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi = test->odi;
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *yaml_str;

  test_ncm_obj_dict_int_add (test, pdata);

  yaml_str = ncm_serialize_dict_int_to_yaml (ser, odi);

  {
    NcmObjDictInt *odi_load = ncm_serialize_dict_int_from_yaml (ser, yaml_str);
    GHashTableIter iter;
    gint *key;
    GObject *value;

    g_assert_cmpuint (ncm_obj_dict_int_len (odi), ==, ncm_obj_dict_int_len (odi_load));

    g_hash_table_iter_init (&iter, odi);

    while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
    {
      GObject *obj = ncm_obj_dict_int_get (odi_load, *key);

      g_assert_nonnull (obj);
      g_assert_nonnull (ncm_obj_dict_int_peek (odi_load, *key));
      g_assert_true (obj == ncm_obj_dict_int_peek (odi_load, *key));
      g_assert_true (G_IS_OBJECT (obj));

      g_assert_true (ncm_obj_dict_int_peek (odi, *key) != ncm_obj_dict_int_peek (odi_load, *key));

      g_object_unref (obj);
    }

    ncm_obj_dict_int_unref (odi_load);
  }

  g_free (yaml_str);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_dict_int_to_from_yaml_file (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi = test->odi;
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *tmp_dir     = g_dir_make_tmp ("test_ncm_obj_dict_int_saved_XXXXXX", NULL);
  gchar *filename    = g_strdup_printf ("%s/test_ncm_obj_dict_int_saved.yaml", tmp_dir);

  test_ncm_obj_dict_int_add (test, pdata);

  ncm_serialize_dict_int_to_yaml_file (ser, odi, filename);

  {
    NcmObjDictInt *odi_load = ncm_serialize_dict_int_from_yaml_file (ser, filename);
    GHashTableIter iter;
    gint *key;
    GObject *value;

    g_assert_cmpuint (ncm_obj_dict_int_len (odi), ==, ncm_obj_dict_int_len (odi_load));

    g_hash_table_iter_init (&iter, odi);

    while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
    {
      GObject *obj = ncm_obj_dict_int_get (odi_load, *key);

      g_assert_nonnull (obj);
      g_assert_nonnull (ncm_obj_dict_int_peek (odi_load, *key));
      g_assert_true (obj == ncm_obj_dict_int_peek (odi_load, *key));
      g_assert_true (G_IS_OBJECT (obj));

      g_assert_true (ncm_obj_dict_int_peek (odi, *key) != ncm_obj_dict_int_peek (odi_load, *key));

      g_object_unref (obj);
    }

    ncm_obj_dict_int_unref (odi_load);
  }

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_dict_int_traps (TestNcmObjDictInt *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/obj_dict_int/invalid/add/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/obj_dict_int/invalid/set/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_obj_dict_int_invalid_add (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi = test->odi;

  ncm_obj_dict_int_add (odi, 0, NULL);
}

void
test_ncm_obj_dict_int_invalid_set (TestNcmObjDictInt *test, gconstpointer pdata)
{
  NcmObjDictInt *odi = test->odi;

  ncm_obj_dict_int_set (odi, 0, NULL);
}

/*
 * NcmVarDict
 */

void
test_ncm_var_dict_new (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd;

  test->ntests = 1000;
  vd           = test->vd = ncm_var_dict_new ();

  g_assert_true (vd != NULL);
}

void
test_ncm_var_dict_free (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  ncm_var_dict_unref (vd);
}

void
test_ncm_var_dict_has_key (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    gchar *str = g_strdup ("Hello World!");

    ncm_var_dict_set_string (vd, "str", str);
    ncm_var_dict_set_double (vd, "d", 123.456);
    ncm_var_dict_set_boolean (vd, "b", TRUE);
    ncm_var_dict_set_int (vd, "i", 123);

    g_assert_true (ncm_var_dict_has_key (vd, "str"));
    g_assert_true (ncm_var_dict_has_key (vd, "d"));
    g_assert_true (ncm_var_dict_has_key (vd, "b"));
    g_assert_true (ncm_var_dict_has_key (vd, "i"));

    g_assert_false (ncm_var_dict_has_key (vd, "not_exist"));
    g_assert_false (ncm_var_dict_has_key (vd, "not_exist_either"));

    g_free (str);
  }
}

void
test_ncm_var_dict_string (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    gchar *str = g_strdup ("Hello World!");

    ncm_var_dict_set_string (vd, "str", str);

    {
      gchar *str0;

      g_assert_false (ncm_var_dict_get_string (vd, "not_exist", &str0));
      g_assert_true (ncm_var_dict_get_string (vd, "str", &str0));

      g_assert_cmpstr (str0, ==, str);

      g_free (str0);
    }

    g_free (str);
  }
}

void
test_ncm_var_dict_int (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    gint i = 123;

    ncm_var_dict_set_int (vd, "i", i);

    {
      gint i0;

      g_assert_false (ncm_var_dict_get_int (vd, "not_exist", &i0));
      g_assert_true (ncm_var_dict_get_int (vd, "i", &i0));

      g_assert_cmpint (i0, ==, i);
    }
  }
}

void
test_ncm_var_dict_double (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    gdouble d = 123.456;

    ncm_var_dict_set_double (vd, "d", d);

    {
      gdouble d0;

      g_assert_false (ncm_var_dict_get_double (vd, "not_exist", &d0));
      g_assert_true (ncm_var_dict_get_double (vd, "d", &d0));

      g_assert_cmpfloat (d0, ==, d);
    }
  }
}

void
test_ncm_var_dict_boolean (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    gboolean b = TRUE;

    ncm_var_dict_set_boolean (vd, "b", b);

    {
      gboolean b0;

      g_assert_false (ncm_var_dict_get_boolean (vd, "not_exist", &b0));
      g_assert_true (ncm_var_dict_get_boolean (vd, "b", &b0));

      g_assert_cmpint (b0, ==, b);
    }
  }
}

void
test_ncm_var_dict_int_array (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    GArray *iarr = g_array_new (FALSE, FALSE, sizeof (gint));
    gint i;

    for (i = 0; i < 5; i++)
    {
      gint i = g_test_rand_int ();

      g_array_append_val (iarr, i);
    }

    ncm_var_dict_set_int_array (vd, "iarr", iarr);

    {
      GArray *iarr0;
      guint i;

      g_assert_false (ncm_var_dict_get_int_array (vd, "not_exist", &iarr0));
      g_assert_true (ncm_var_dict_get_int_array (vd, "iarr", &iarr0));

      g_assert_cmpuint (iarr->len, ==, iarr0->len);

      for (i = 0; i < iarr->len; i++)
      {
        g_assert_cmpint (g_array_index (iarr, gint, i), ==, g_array_index (iarr0, gint, i));
      }

      g_array_unref (iarr0);
    }

    g_array_unref (iarr);
  }
}

void
test_ncm_var_dict_double_array (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    GArray *darr = g_array_new (FALSE, FALSE, sizeof (gdouble));
    guint i;

    for (i = 0; i < 5; i++)
    {
      gdouble d = g_test_rand_double ();

      g_array_append_val (darr, d);
    }

    ncm_var_dict_set_double_array (vd, "darr", darr);

    {
      GArray *darr0;

      g_assert_false (ncm_var_dict_get_double_array (vd, "not_exist", &darr0));
      g_assert_true (ncm_var_dict_get_double_array (vd, "darr", &darr0));

      g_assert_cmpuint (darr->len, ==, darr0->len);

      for (i = 0; i < darr->len; i++)
      {
        g_assert_cmpfloat (g_array_index (darr, gdouble, i), ==, g_array_index (darr0, gdouble, i));
      }

      g_array_unref (darr0);
    }

    g_array_unref (darr);
  }
}

void
test_ncm_var_dict_boolean_array (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    GArray *barr = g_array_new (FALSE, FALSE, sizeof (gboolean));
    guint i;

    for (i = 0; i < 5; i++)
    {
      gboolean b = g_test_rand_int () % 2 == 0;

      g_array_append_val (barr, b);
    }

    ncm_var_dict_set_boolean_array (vd, "barr", barr);

    {
      GArray *barr0;

      g_assert_false (ncm_var_dict_get_boolean_array (vd, "not_exist", &barr0));
      g_assert_true (ncm_var_dict_get_boolean_array (vd, "barr", &barr0));

      g_assert_cmpuint (barr->len, ==, barr0->len);

      for (i = 0; i < barr->len; i++)
      {
        g_assert_cmpint (g_array_index (barr, gboolean, i), ==, g_array_index (barr0, gboolean, i));
      }

      g_array_unref (barr0);
    }

    g_array_unref (barr);
  }
}

void
test_ncm_var_dict_variant (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  g_assert_true (vd != NULL);

  {
    GVariant *v = g_variant_new_string ("Hello World!");

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }

  {
    GVariant *v = g_variant_new_int32 (123);

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }

  {
    GVariant *v = g_variant_new_double (123.456);

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }

  {
    GVariant *v = g_variant_new_boolean (TRUE);

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }

  {
    gint array[] = { 2341, 46822, -223 };
    GVariant *v  = g_variant_new_fixed_array (G_VARIANT_TYPE_INT32, array, G_N_ELEMENTS (array), sizeof (gint32));

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }

  {
    gdouble array[] = { 2341.123, 46822.456, -223.789 };
    GVariant *v     = g_variant_new_fixed_array (G_VARIANT_TYPE_DOUBLE, array, G_N_ELEMENTS (array), sizeof (gdouble));

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }

  {
    gchar array[] = { TRUE, FALSE, TRUE };
    GVariant *v   = g_variant_new_fixed_array (G_VARIANT_TYPE_BOOLEAN, array, G_N_ELEMENTS (array), sizeof (gchar));

    ncm_var_dict_set_variant (vd, "v", v);

    {
      GVariant *v0;

      g_assert_false (ncm_var_dict_get_variant (vd, "not_exist", &v0));
      g_assert_true (ncm_var_dict_get_variant (vd, "v", &v0));

      g_assert_true (g_variant_equal (v, v0));

      g_variant_unref (v0);
    }
  }
}

NcmVarDict *
_test_ncm_var_dict_to_from_variant (NcmVarDict *vd)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  GVariant *v       = ncm_serialize_var_dict_to_variant (ser, vd);
  NcmVarDict *vd0   = ncm_serialize_var_dict_from_variant (ser, v);

  g_assert_nonnull (vd0);
  g_assert_nonnull (v);

  g_assert_true (vd != vd0);

  g_variant_unref (v);
  ncm_serialize_clear (&ser);

  return vd0;
}

NcmVarDict *
_test_ncm_var_dict_to_from_yaml (NcmVarDict *vd)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  gchar *yaml_str   = ncm_serialize_var_dict_to_yaml (ser, vd);
  NcmVarDict *vd0   = ncm_serialize_var_dict_from_yaml (ser, yaml_str);

  g_assert_nonnull (vd0);
  g_assert_nonnull (yaml_str);

  g_assert_true (vd != vd0);

  g_free (yaml_str);
  ncm_serialize_clear (&ser);

  return vd0;
}

NcmVarDict *
_test_ncm_var_dict_to_from_variant_binfile (NcmVarDict *vd)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmVarDict *vd0   = NULL;
  gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_var_dict_saved_XXXXXX", NULL);
  gchar *filename   = g_strdup_printf ("%s/test_ncm_var_dict_saved.yaml", tmp_dir);

  ncm_serialize_var_dict_to_variant_file (ser, vd, filename, TRUE);
  vd0 = ncm_serialize_var_dict_from_variant_file (ser, filename, TRUE);

  g_assert_nonnull (vd0);
  g_assert_true (vd != vd0);

  ncm_serialize_clear (&ser);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  return vd0;
}

NcmVarDict *
_test_ncm_var_dict_to_from_variant_file (NcmVarDict *vd)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmVarDict *vd0   = NULL;
  gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_var_dict_saved_XXXXXX", NULL);
  gchar *filename   = g_strdup_printf ("%s/test_ncm_var_dict_saved.yaml", tmp_dir);

  ncm_serialize_var_dict_to_variant_file (ser, vd, filename, FALSE);
  vd0 = ncm_serialize_var_dict_from_variant_file (ser, filename, FALSE);

  g_assert_nonnull (vd0);
  g_assert_true (vd != vd0);

  ncm_serialize_clear (&ser);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  return vd0;
}

NcmVarDict *
_test_ncm_var_dict_to_from_yaml_file (NcmVarDict *vd)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmVarDict *vd0   = NULL;
  gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_var_dict_saved_XXXXXX", NULL);
  gchar *filename   = g_strdup_printf ("%s/test_ncm_var_dict_saved.yaml", tmp_dir);

  ncm_serialize_var_dict_to_yaml_file (ser, vd, filename);
  vd0 = ncm_serialize_var_dict_from_yaml_file (ser, filename);

  g_assert_nonnull (vd0);
  g_assert_true (vd != vd0);

  ncm_serialize_clear (&ser);

  g_unlink (filename);
  g_rmdir (tmp_dir);

  g_free (filename);
  g_free (tmp_dir);

  return vd0;
}

void
test_ncm_var_dict_to_from (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd   = test->vd;
  GArray *iarr     = g_array_new (FALSE, FALSE, sizeof (gint));
  GArray *darr     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *barr     = g_array_new (FALSE, FALSE, sizeof (gboolean));
  guint n_ints     = g_test_rand_int_range (5, 10);
  guint n_doubles  = g_test_rand_int_range (5, 10);
  guint n_booleans = g_test_rand_int_range (5, 10);

  NcmVarDict *(*to_from) (NcmVarDict *) = pdata;

  guint i;

  for (i = 0; i < n_ints; i++)
  {
    gint i = g_test_rand_int ();

    g_array_append_val (iarr, i);
  }

  for (i = 0; i < n_doubles; i++)
  {
    gdouble d = g_test_rand_double ();

    g_array_append_val (darr, d);
  }

  for (i = 0; i < n_booleans; i++)
  {
    gboolean b = g_test_rand_int () % 2 == 0;

    g_array_append_val (barr, b);
  }

  ncm_var_dict_set_string (vd, "str", "Hello World!");
  ncm_var_dict_set_int (vd, "i", 123);
  ncm_var_dict_set_double (vd, "d", 123.456);
  ncm_var_dict_set_boolean (vd, "b", TRUE);
  ncm_var_dict_set_int_array (vd, "iarr", iarr);
  ncm_var_dict_set_double_array (vd, "darr", darr);
  ncm_var_dict_set_boolean_array (vd, "barr", barr);

  {
    NcmVarDict *vd0 = to_from (vd);

    g_assert_true (vd0 != NULL);

    g_assert_cmpuint (ncm_var_dict_len (vd), ==, ncm_var_dict_len (vd0));

    g_assert_true (ncm_var_dict_has_key (vd0, "str"));
    g_assert_true (ncm_var_dict_has_key (vd0, "i"));
    g_assert_true (ncm_var_dict_has_key (vd0, "d"));
    g_assert_true (ncm_var_dict_has_key (vd0, "b"));
    g_assert_true (ncm_var_dict_has_key (vd0, "iarr"));
    g_assert_true (ncm_var_dict_has_key (vd0, "darr"));
    g_assert_true (ncm_var_dict_has_key (vd0, "barr"));

    g_assert_false (ncm_var_dict_has_key (vd0, "not_exist"));
    g_assert_false (ncm_var_dict_has_key (vd0, "not_exist_either"));

    {
      gchar *str;
      gchar *str0;

      g_assert_true (ncm_var_dict_get_string (vd, "str", &str));
      g_assert_true (ncm_var_dict_get_string (vd0, "str", &str0));

      g_assert_cmpstr (str, ==, str0);

      g_free (str);
      g_free (str0);
    }

    {
      gint i;
      gint i0;

      g_assert_true (ncm_var_dict_get_int (vd, "i", &i));
      g_assert_true (ncm_var_dict_get_int (vd0, "i", &i0));

      g_assert_cmpint (i, ==, i0);
    }

    {
      gdouble d;
      gdouble d0;

      g_assert_true (ncm_var_dict_get_double (vd, "d", &d));
      g_assert_true (ncm_var_dict_get_double (vd0, "d", &d0));

      g_assert_cmpfloat (d, ==, d0);
    }

    {
      gboolean b;
      gboolean b0;

      g_assert_true (ncm_var_dict_get_boolean (vd, "b", &b));
      g_assert_true (ncm_var_dict_get_boolean (vd0, "b", &b0));

      g_assert_cmpint (b, ==, b0);
    }

    {
      GArray *iarr0;
      guint i;

      g_assert_true (ncm_var_dict_get_int_array (vd, "iarr", &iarr0));

      g_assert_cmpuint (iarr->len, ==, iarr0->len);

      for (i = 0; i < iarr->len; i++)
      {
        g_assert_cmpint (g_array_index (iarr, gint, i), ==, g_array_index (iarr0, gint, i));
      }

      g_array_unref (iarr0);
    }

    {
      GArray *darr0;
      guint i;

      g_assert_true (ncm_var_dict_get_double_array (vd, "darr", &darr0));

      g_assert_cmpuint (darr->len, ==, darr0->len);

      for (i = 0; i < darr->len; i++)
      {
        g_assert_cmpfloat (g_array_index (darr, gdouble, i), ==, g_array_index (darr0, gdouble, i));
      }

      g_array_unref (darr0);
    }

    {
      GArray *barr0;
      guint i;

      g_assert_true (ncm_var_dict_get_boolean_array (vd, "barr", &barr0));

      g_assert_cmpuint (barr->len, ==, barr0->len);

      for (i = 0; i < barr->len; i++)
      {
        g_assert_cmpint (g_array_index (barr, gboolean, i), ==, g_array_index (barr0, gboolean, i));
      }

      g_array_unref (barr0);
    }

    ncm_var_dict_clear (&vd0);
  }

  g_array_free (iarr, TRUE);
  g_array_free (darr, TRUE);
  g_array_free (barr, TRUE);
}

void
test_ncm_var_dict_traps (TestNcmVarDict *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/var_dict/invalid/set/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_var_dict_invalid_set (TestNcmVarDict *test, gconstpointer pdata)
{
  NcmVarDict *vd = test->vd;

  ncm_var_dict_set_string (vd, NULL, NULL);
}


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

  test_ncm_obj_array_add (test, pdata);

  ncm_serialize_array_to_yaml_file (ser, oa, "test_ncm_obj_array_saved.yaml");

  {
    NcmObjArray *oa_load = ncm_serialize_array_from_yaml_file (ser, "test_ncm_obj_array_saved.yaml");
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
    nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));

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

  test_ncm_obj_dict_str_add (test, pdata);

  ncm_serialize_dict_str_to_yaml_file (ser, ods, "test_ncm_obj_dict_str_saved.yaml");

  {
    NcmObjDictStr *ods_load = ncm_serialize_dict_str_from_yaml_file (ser, "test_ncm_obj_dict_str_saved.yaml");
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
    }

    ncm_obj_dict_str_unref (ods_load);
  }

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
    nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo));

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

  test_ncm_obj_dict_int_add (test, pdata);

  ncm_serialize_dict_int_to_yaml_file (ser, odi, "test_ncm_obj_dict_int_saved.yaml");

  {
    NcmObjDictInt *odi_load = ncm_serialize_dict_int_from_yaml_file (ser, "test_ncm_obj_dict_int_saved.yaml");
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
    }

    ncm_obj_dict_int_unref (odi_load);
  }

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


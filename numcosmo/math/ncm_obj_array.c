/***************************************************************************
 *            ncm_obj_array.c
 *
 *  Wed October 16 11:04:01 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_obj_array.c
 *
 * Copyright (C) 2013 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:ncm_obj_array
 * @title: NcmObjArray
 * @short_description: GObjects array with serialization support.
 *
 * A #NcmObjArray is a #GPtrArray that holds #GObject's. It is used to
 * store #GObject's that can be serialized to a #GVariant.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_obj_array.h"
#include "math/ncm_cfg.h"

G_DEFINE_BOXED_TYPE (NcmObjArray, ncm_obj_array, ncm_obj_array_ref, ncm_obj_array_unref)
G_DEFINE_BOXED_TYPE (NcmObjDictStr, ncm_obj_dict_str, ncm_obj_dict_str_ref, ncm_obj_dict_str_unref)
G_DEFINE_BOXED_TYPE (NcmObjDictInt, ncm_obj_dict_int, ncm_obj_dict_int_ref, ncm_obj_dict_int_unref)
G_DEFINE_BOXED_TYPE (NcmVarDict, ncm_var_dict, ncm_var_dict_ref, ncm_var_dict_unref)

/*
 * NcmObjArray
 */

/**
 * ncm_obj_array_new:
 *
 * Creates a new #NcmObjArray.
 *
 * Returns: (transfer full): a new #NcmObjArray.
 */
NcmObjArray *
ncm_obj_array_new ()
{
  GPtrArray *oa = g_ptr_array_new ();

  g_ptr_array_set_free_func (oa, g_object_unref);

  return (NcmObjArray *) oa;
}

/**
 * ncm_obj_array_sized_new:
 * @n: initial allocation size.
 *
 * Creates a new #NcmObjArray with @n elements preallocated.
 *
 * Returns: (transfer full): a new #NcmObjArray.
 */
NcmObjArray *
ncm_obj_array_sized_new (guint n)
{
  GPtrArray *oa = g_ptr_array_sized_new (n);

  g_ptr_array_set_free_func (oa, g_object_unref);

  return (NcmObjArray *) oa;
}

/**
 * ncm_obj_array_ref:
 * @oa: a #NcmObjArray.
 *
 * Increases the reference count of @oa by one.
 *
 * Returns: (transfer full): @oa.
 */
NcmObjArray *
ncm_obj_array_ref (NcmObjArray *oa)
{
  return (NcmObjArray *) g_ptr_array_ref ((GPtrArray *) oa);
}

/**
 * ncm_obj_array_unref:
 * @oa: a #NcmObjArray.
 *
 * Decreases the reference count of @oa by one. If the reference count
 * reaches zero, all objects in the array are unreferenced.
 *
 */
void
ncm_obj_array_unref (NcmObjArray *oa)
{
  g_ptr_array_unref ((GPtrArray *) oa);
}

/**
 * ncm_obj_array_clear:
 * @oa: a pointer to a #NcmObjArray.
 *
 * If *@oa is not %NULL, unreferences it and sets *@oa to %NULL.
 *
 */
void
ncm_obj_array_clear (NcmObjArray **oa)
{
  g_clear_pointer (oa, ncm_obj_array_unref);
}

/**
 * ncm_obj_array_add:
 * @oa: a #NcmObjArray.
 * @obj: a #GObject.
 *
 * Adds a #GObject to a #NcmObjArray.
 *
 */
void
ncm_obj_array_add (NcmObjArray *oa, GObject *obj)
{
  g_assert (obj != NULL);
  g_ptr_array_add ((GPtrArray *) oa, g_object_ref (obj));
}

/**
 * ncm_obj_array_set:
 * @oa: a #NcmObjArray.
 * @i: object index.
 * @obj: a #GObject.
 *
 * Sets a #GObject to a #NcmObjArray. If there is already a #GObject
 * at position @i, it is unreferenced.
 *
 */
void
ncm_obj_array_set (NcmObjArray *oa, guint i, GObject *obj)
{
  g_assert_cmpuint (i, <, oa->len);

  g_assert (obj != NULL);

  if (obj != g_ptr_array_index ((GPtrArray *) oa, i))
  {
    g_object_unref (g_ptr_array_index ((GPtrArray *) oa, i));
    g_ptr_array_index ((GPtrArray *) oa, i) = g_object_ref (obj);
  }
}

/**
 * ncm_obj_array_get:
 * @oa: a #NcmObjArray
 * @i: object index
 *
 * Gets a #GObject from a #NcmObjArray at position @i.
 *
 * Returns: (transfer full): the #GObject at position @i.
 */
GObject *
ncm_obj_array_get (NcmObjArray *oa, guint i)
{
  return g_object_ref (ncm_obj_array_peek (oa, i));
}

/**
 * ncm_obj_array_peek:
 * @oa: a #NcmObjArray.
 * @i: object index.
 *
 * Peeks a #GObject from a #NcmObjArray at position @i without increasing its reference
 * count.
 *
 * Returns: (transfer none): the #GObject at position @i.
 */
GObject *
ncm_obj_array_peek (NcmObjArray *oa, guint i)
{
  g_assert_cmpuint (i, <, oa->len);

  return g_ptr_array_index ((GPtrArray *) oa, i);
}

/**
 * ncm_obj_array_len:
 * @oa: a #NcmObjArray
 *
 * Gets the length of a #NcmObjArray.
 *
 * Returns: array length
 */
guint
ncm_obj_array_len (NcmObjArray *oa)
{
  return oa->len;
}

/*
 * NcmObjDictStr
 */

/**
 * ncm_obj_dict_str_new:
 *
 * Creates a new #NcmObjDictStr.
 *
 * Returns: (transfer full): a new #NcmObjDictStr.
 */
NcmObjDictStr *
ncm_obj_dict_str_new ()
{
  GHashTable *ods = g_hash_table_new_full (g_str_hash, g_str_equal, g_free, g_object_unref);

  return (NcmObjDictStr *) ods;
}

/**
 * ncm_obj_dict_str_ref:
 * @ods: a #NcmObjDictStr.
 *
 * Increases the reference count of @ods by one.
 *
 * Returns: (transfer full): @ods.
 */
NcmObjDictStr *
ncm_obj_dict_str_ref (NcmObjDictStr *ods)
{
  return (NcmObjDictStr *) g_hash_table_ref ((GHashTable *) ods);
}

/**
 * ncm_obj_dict_str_unref:
 * @ods: a #NcmObjDictStr.
 *
 * Decreases the reference count of @ods by one. If the reference count
 * reaches zero, all objects in the dictionary are unreferenced.
 *
 */
void
ncm_obj_dict_str_unref (NcmObjDictStr *ods)
{
  g_hash_table_unref ((GHashTable *) ods);
}

/**
 * ncm_obj_dict_str_clear:
 * @ods: a pointer to a #NcmObjDictStr.
 *
 * If *@ods is not %NULL, unreferences it and sets *@ods to %NULL.
 *
 */
void
ncm_obj_dict_str_clear (NcmObjDictStr **ods)
{
  g_clear_pointer (ods, ncm_obj_dict_str_unref);
}

/**
 * ncm_obj_dict_str_add:
 * @ods: a #NcmObjDictStr.
 * @key: a string.
 * @obj: a #GObject.
 *
 * Adds a #GObject to a #NcmObjDictStr.
 *
 */
void
ncm_obj_dict_str_add (NcmObjDictStr *ods, const gchar *key, GObject *obj)
{
  g_assert (key != NULL);
  g_assert (obj != NULL);

  g_hash_table_insert ((GHashTable *) ods, g_strdup (key), g_object_ref (obj));
}

/**
 * ncm_obj_dict_str_set:
 * @ods: a #NcmObjDictStr.
 * @key: a string.
 * @obj: a #GObject.
 *
 * Sets a #GObject to a #NcmObjDictStr. If there is already a #GObject
 * with key @key, it is unreferenced.
 *
 */
void
ncm_obj_dict_str_set (NcmObjDictStr *ods, const gchar *key, GObject *obj)
{
  g_assert (key != NULL);
  g_assert (obj != NULL);

  if (obj != g_hash_table_lookup ((GHashTable *) ods, key))
    g_hash_table_replace ((GHashTable *) ods, g_strdup (key), g_object_ref (obj));
}

/**
 * ncm_obj_dict_str_get:
 * @ods: a #NcmObjDictStr
 * @key: a string.
 *
 * Gets a #GObject from a #NcmObjDictStr with key @key.
 *
 * Returns: (transfer full): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_str_get (NcmObjDictStr *ods, const gchar *key)
{
  GObject *obj = ncm_obj_dict_str_peek (ods, key);

  if (obj != NULL)
    return g_object_ref (obj);

  return NULL;
}

/**
 * ncm_obj_dict_str_peek:
 * @ods: a #NcmObjDictStr.
 * @key: a string.
 *
 * Peeks a #GObject from a #NcmObjDictStr with key @key without increasing its reference
 * count.
 *
 * Returns: (transfer none): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_str_peek (NcmObjDictStr *ods, const gchar *key)
{
  g_assert (key != NULL);

  return g_hash_table_lookup ((GHashTable *) ods, key);
}

/**
 * ncm_obj_dict_str_len:
 * @ods: a #NcmObjDictStr
 *
 * Gets the length of a #NcmObjDictStr.
 *
 * Returns: dictionary length
 */
guint
ncm_obj_dict_str_len (NcmObjDictStr *ods)
{
  return g_hash_table_size ((GHashTable *) ods);
}

/**
 * ncm_obj_dict_str_keys:
 * @ods: a #NcmObjDictStr
 *
 * Gets the keys of a #NcmObjDictStr.
 *
 * Returns: (transfer container): the keys of a #NcmObjDictStr.
 */
GStrv
ncm_obj_dict_str_keys (NcmObjDictStr *ods)
{
  return (GStrv) g_hash_table_get_keys_as_array ((GHashTable *) ods, NULL);
}

/*
 * NcmObjDictInt
 */

/**
 * ncm_obj_dict_int_new:
 *
 * Creates a new #NcmObjDictInt.
 *
 * Returns: (transfer full): a new #NcmObjDictInt.
 */
NcmObjDictInt *
ncm_obj_dict_int_new ()
{
  GHashTable *odi = g_hash_table_new_full (g_int_hash, g_int_equal, g_free, g_object_unref);

  return (NcmObjDictInt *) odi;
}

/**
 * ncm_obj_dict_int_ref:
 * @odi: a #NcmObjDictInt.
 *
 * Increases the reference count of @odi by one.
 *
 * Returns: (transfer full): @odi.
 */
NcmObjDictInt *
ncm_obj_dict_int_ref (NcmObjDictInt *odi)
{
  return (NcmObjDictInt *) g_hash_table_ref ((GHashTable *) odi);
}

/**
 * ncm_obj_dict_int_unref:
 * @odi: a #NcmObjDictInt.
 *
 * Decreases the reference count of @odi by one. If the reference count
 * reaches zero, all objects in the dictionary are unreferenced.
 *
 */
void
ncm_obj_dict_int_unref (NcmObjDictInt *odi)
{
  g_hash_table_unref ((GHashTable *) odi);
}

/**
 * ncm_obj_dict_int_clear:
 * @odi: a pointer to a #NcmObjDictInt.
 *
 * If *@odi is not %NULL, unreferences it and sets *@odi to %NULL.
 *
 */
void
ncm_obj_dict_int_clear (NcmObjDictInt **odi)
{
  g_clear_pointer (odi, ncm_obj_dict_int_unref);
}

/**
 * ncm_obj_dict_int_add:
 * @odi: a #NcmObjDictInt.
 * @key: an integer.
 * @obj: a #GObject.
 *
 * Adds a #GObject to a #NcmObjDictInt.
 *
 */
void
ncm_obj_dict_int_add (NcmObjDictInt *odi, gint key, GObject *obj)
{
  g_assert (obj != NULL);

  g_hash_table_insert ((GHashTable *) odi, g_memdup2 (&key, sizeof (gint)), g_object_ref (obj));
}

/**
 * ncm_obj_dict_int_set:
 * @odi: a #NcmObjDictInt.
 * @key: an integer.
 * @obj: a #GObject.
 *
 * Sets a #GObject to a #NcmObjDictInt. If there is already a #GObject
 * with key @key, it is unreferenced.
 *
 */
void
ncm_obj_dict_int_set (NcmObjDictInt *odi, gint key, GObject *obj)
{
  g_assert (obj != NULL);

  if (obj != g_hash_table_lookup ((GHashTable *) odi, &key))
    g_hash_table_replace ((GHashTable *) odi, g_memdup2 (&key, sizeof (gint)), g_object_ref (obj));
}

/**
 * ncm_obj_dict_int_get:
 * @odi: a #NcmObjDictInt
 * @key: an integer.
 *
 * Gets a #GObject from a #NcmObjDictInt with key @key.
 *
 * Returns: (transfer full): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_int_get (NcmObjDictInt *odi, gint key)
{
  GObject *obj = ncm_obj_dict_int_peek (odi, key);

  if (obj != NULL)
    return g_object_ref (obj);

  return NULL;
}

/**
 * ncm_obj_dict_int_peek:
 * @odi: a #NcmObjDictInt.
 * @key: an integer.
 *
 * Peeks a #GObject from a #NcmObjDictInt with key @key without increasing its reference
 * count.
 *
 * Returns: (transfer none): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_int_peek (NcmObjDictInt *odi, gint key)
{
  return g_hash_table_lookup ((GHashTable *) odi, &key);
}

/**
 * ncm_obj_dict_int_len:
 * @odi: a #NcmObjDictInt
 *
 * Gets the length of a #NcmObjDictInt.
 *
 * Returns: dictionary length
 */
guint
ncm_obj_dict_int_len (NcmObjDictInt *odi)
{
  return g_hash_table_size ((GHashTable *) odi);
}

/**
 * ncm_obj_dict_int_keys:
 * @odi: a #NcmObjDictInt
 *
 * Gets the keys of a #NcmObjDictInt.
 *
 * Returns: (transfer full) (array) (element-type int): the keys of a #NcmObjDictInt.
 */
GArray *
ncm_obj_dict_int_keys (NcmObjDictInt *odi)
{
  GArray *keys = g_array_new (FALSE, FALSE, sizeof (gint));
  GHashTableIter iter;
  gint *key;

  g_hash_table_iter_init (&iter, (GHashTable *) odi);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, NULL))
    g_array_append_val (keys, *key);

  return keys;
}

/*
 * NcmVarDict
 */

/**
 * ncm_var_dict_new:
 *
 * Creates a new #NcmVarDict.
 *
 * Returns: (transfer full): a new #NcmVarDict.
 */
NcmVarDict *
ncm_var_dict_new ()
{
  GHashTable *vd = g_hash_table_new_full (g_str_hash, g_str_equal, g_free, (GDestroyNotify) g_variant_unref);

  return (NcmVarDict *) vd;
}

/**
 * ncm_var_dict_ref:
 * @vd: a #NcmVarDict.
 *
 * Increases the reference count of @vd by one.
 *
 * Returns: (transfer full): @vd.
 */
NcmVarDict *
ncm_var_dict_ref (NcmVarDict *vd)
{
  return (NcmVarDict *) g_hash_table_ref ((GHashTable *) vd);
}

/**
 * ncm_var_dict_unref:
 * @vd: a #NcmVarDict.
 *
 * Decreases the reference count of @vd by one. If the reference count
 * reaches zero, all objects in the dictionary are unreferenced.
 *
 */
void
ncm_var_dict_unref (NcmVarDict *vd)
{
  g_hash_table_unref ((GHashTable *) vd);
}

/**
 * ncm_var_dict_clear:
 * @vd: a pointer to a #NcmVarDict.
 *
 * If *@vd is not %NULL, unreferences it and sets *@vd to %NULL.
 *
 */
void
ncm_var_dict_clear (NcmVarDict **vd)
{
  g_clear_pointer (vd, ncm_var_dict_unref);
}

/**
 * ncm_var_dict_peek:
 * @vd: a #NcmVarDict.
 * @key: a string.
 *
 * Peeks a #GVariant from a #NcmVarDict with key @key without increasing its reference
 * count.
 *
 * Returns: (transfer none): the #GVariant with key @key.
 */
GVariant *
ncm_var_dict_peek (NcmVarDict *vd, const gchar *key)
{
  g_assert (key != NULL);

  return g_hash_table_lookup ((GHashTable *) vd, key);
}

/**
 * ncm_var_dict_set_string:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: a string
 *
 * Sets a string to a #NcmVarDict. If there is already a string
 * with key @key, it is unreferenced.
 *
 */
void
ncm_var_dict_set_string (NcmVarDict *vd, const gchar *key, const gchar *value)
{
  g_assert (key != NULL);
  g_assert (value != NULL);

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                       g_variant_ref_sink (g_variant_new_string (value)));
}

/**
 * ncm_var_dict_set_int:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: an integer
 *
 * Sets an integer to a #NcmVarDict. If there is already an integer
 * with key @key, it is unreferenced.
 *
 */
void
ncm_var_dict_set_int (NcmVarDict *vd, const gchar *key, gint value)
{
  g_assert (key != NULL);

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                       g_variant_ref_sink (g_variant_new_int32 (value)));
}

/**
 * ncm_var_dict_set_double:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: a double
 *
 * Sets a double to a #NcmVarDict. If there is already a double
 * with key @key, it is unreferenced.
 *
 */
void
ncm_var_dict_set_double (NcmVarDict *vd, const gchar *key, gdouble value)
{
  g_assert (key != NULL);

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                       g_variant_ref_sink (g_variant_new_double (value)));
}

/**
 * ncm_var_dict_set_boolean:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: a boolean
 *
 * Sets a boolean to a #NcmVarDict. If there is already a boolean
 * with key @key, it is unreferenced.
 */
void
ncm_var_dict_set_boolean (NcmVarDict *vd, const gchar *key, gboolean value)
{
  g_assert (key != NULL);

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                       g_variant_ref_sink (g_variant_new_boolean (value)));
}

/**
 * ncm_var_dict_set_int_array:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (array) (element-type int): a #GArray of integers
 *
 * Sets an array of integers to a #NcmVarDict. If there is already an array
 * of integers with key @key, it is unreferenced.
 *
 */
void
ncm_var_dict_set_int_array (NcmVarDict *vd, const gchar *key, GArray *value)
{
  g_assert (key != NULL);
  g_assert (value != NULL);

  g_assert (g_array_get_element_size (value) == sizeof (gint));

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                       g_variant_ref_sink (g_variant_new_fixed_array (G_VARIANT_TYPE_INT32,
                                                                      value->data,
                                                                      value->len,
                                                                      sizeof (gint))));
}

/**
 * ncm_var_dict_set_double_array:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (array) (element-type double): a #GArray of doubles
 *
 * Sets an array of doubles to a #NcmVarDict. If there is already an array
 * of doubles with key @key, it is unreferenced.
 */
void
ncm_var_dict_set_double_array (NcmVarDict *vd, const gchar *key, GArray *value)
{
  g_assert (key != NULL);
  g_assert (value != NULL);

  g_assert (g_array_get_element_size (value) == sizeof (gdouble));

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                       g_variant_ref_sink (g_variant_new_fixed_array (G_VARIANT_TYPE_DOUBLE,
                                                                      value->data,
                                                                      value->len,
                                                                      sizeof (gdouble))));
}

/**
 * ncm_var_dict_set_boolean_array:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (array) (element-type boolean): a #GArray of booleans
 *
 * Sets an array of booleans to a #NcmVarDict. If there is already an array
 * of booleans with key @key, it is unreferenced.
 */
void
ncm_var_dict_set_boolean_array (NcmVarDict *vd, const gchar *key, GArray *value)
{
  g_assert (key != NULL);
  g_assert (value != NULL);

  g_assert (g_array_get_element_size (value) == sizeof (gboolean));

  {
    GVariantBuilder builder;
    guint i;

    g_variant_builder_init (&builder, G_VARIANT_TYPE ("ab"));

    for (i = 0; i < value->len; i++)
    {
      gchar b = g_array_index (value, gboolean, i);

      g_variant_builder_add (&builder, "b", b);
    }

    g_hash_table_insert ((GHashTable *) vd, g_strdup (key),
                         g_variant_ref_sink (g_variant_builder_end (&builder)));
  }
}

/**
 * ncm_var_dict_set_variant:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: a #GVariant
 *
 * Sets a #GVariant to a #NcmVarDict. If there is already a #GVariant
 * with key @key, it is unreferenced.
 *
 * Valid #GVariant types are:
 * - %G_VARIANT_TYPE_STRING
 * - %G_VARIANT_TYPE_INT32
 * - %G_VARIANT_TYPE_DOUBLE
 * - %G_VARIANT_TYPE_BOOLEAN
 * - %G_VARIANT_TYPE_ARRAY (element-type int)
 * - %G_VARIANT_TYPE_ARRAY (element-type double)
 * - %G_VARIANT_TYPE_ARRAY (element-type boolean)
 *
 */
void
ncm_var_dict_set_variant (NcmVarDict *vd, const gchar *key, GVariant *value)
{
  const GVariantType *allowed_types[] = {
    G_VARIANT_TYPE_STRING,
    G_VARIANT_TYPE_INT32,
    G_VARIANT_TYPE_DOUBLE,
    G_VARIANT_TYPE_BOOLEAN,
    G_VARIANT_TYPE ("ai"),
    G_VARIANT_TYPE ("ad"),
    G_VARIANT_TYPE ("ab"),
    NULL
  };
  gboolean is_allowed = FALSE;
  guint i;

  g_assert (key != NULL);
  g_assert (value != NULL);

  for (i = 0; allowed_types[i] != NULL; i++)
  {
    if (g_variant_is_of_type (value, allowed_types[i]))
    {
      is_allowed = TRUE;
      break;
    }
  }

  if (!is_allowed)
    g_error ("ncm_var_dict_set_variant: Invalid GVariant type");

  g_hash_table_insert ((GHashTable *) vd, g_strdup (key), g_variant_ref_sink (value));
}

/**
 * ncm_var_dict_has_key:
 * @vd: a #NcmVarDict
 * @key: a string
 *
 * Checks if a #NcmVarDict has a key @key.
 *
 * Returns: whether the key @key was found.
 */
gboolean
ncm_var_dict_has_key (NcmVarDict *vd, const gchar *key)
{
  return g_hash_table_contains ((GHashTable *) vd, key);
}

/**
 * ncm_var_dict_get_string:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out) (transfer full): a string
 *
 * Gets a string from a #NcmVarDict with key @key.
 *
 * Returns: whether the string with key @key was found.
 */
gboolean
ncm_var_dict_get_string (NcmVarDict *vd, const gchar *key, gchar **value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *value = g_variant_dup_string (v, NULL);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_int:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out): an integer
 *
 * Gets an integer from a #NcmVarDict with key @key.
 *
 * Returns: whether the integer with key @key was found.
 */
gboolean
ncm_var_dict_get_int (NcmVarDict *vd, const gchar *key, gint *value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *value = g_variant_get_int32 (v);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_double:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out): a double
 *
 * Gets a double from a #NcmVarDict with key @key.
 *
 * Returns: whether the double with key @key was found.
 */
gboolean
ncm_var_dict_get_double (NcmVarDict *vd, const gchar *key, gdouble *value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *value = g_variant_get_double (v);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_boolean:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out): a boolean
 *
 * Gets a boolean from a #NcmVarDict with key @key.
 *
 * Returns: whether the boolean with key @key was found.
 */
gboolean
ncm_var_dict_get_boolean (NcmVarDict *vd, const gchar *key, gboolean *value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *value = g_variant_get_boolean (v);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_int_array:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out) (transfer full) (element-type int): an array of integers
 *
 * Gets an array of integers from a #NcmVarDict with key @key.
 *
 * Returns: whether the array of integers with key @key was found.
 */
gboolean
ncm_var_dict_get_int_array (NcmVarDict *vd, const gchar *key, GArray **value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    gsize len;
    const gint *data = (const gint *) g_variant_get_fixed_array (v, &len, sizeof (gint));

    *value = g_array_sized_new (FALSE, FALSE, sizeof (gint), len);

    g_array_append_vals (*value, data, len);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_double_array:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out) (transfer full) (element-type double): an array of doubles
 *
 * Gets an array of doubles from a #NcmVarDict with key @key.
 *
 * Returns: whether the array of doubles with key @key was found.
 */
gboolean
ncm_var_dict_get_double_array (NcmVarDict *vd, const gchar *key, GArray **value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    gsize len;
    const gdouble *data = (const gdouble *) g_variant_get_fixed_array (v, &len, sizeof (gdouble));

    *value = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);

    g_array_append_vals (*value, data, len);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_boolean_array:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out) (transfer full) (element-type boolean): an array of booleans
 *
 * Gets an array of booleans from a #NcmVarDict with key @key.
 *
 * Returns: whether the array of booleans with key @key was found.
 */
gboolean
ncm_var_dict_get_boolean_array (NcmVarDict *vd, const gchar *key, GArray **value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    gsize len;
    const gchar *data = (const gchar *) g_variant_get_fixed_array (v, &len, sizeof (gchar));
    guint i;

    *value = g_array_sized_new (FALSE, FALSE, sizeof (gboolean), len);

    for (i = 0; i < len; i++)
    {
      gboolean b = data[i];

      g_array_append_val (*value, b);
    }

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_variant:
 * @vd: a #NcmVarDict
 * @key: a string
 * @value: (out) (transfer full): a #GVariant
 *
 * Gets a #GVariant from a #NcmVarDict with key @key.
 *
 * Returns: whether the #GVariant with key @key was found.
 */
gboolean
ncm_var_dict_get_variant (NcmVarDict *vd, const gchar *key, GVariant **value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *value = g_variant_ref_sink (v);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_len:
 * @vd: a #NcmVarDict
 *
 * Gets the length of a #NcmVarDict.
 *
 * Returns: dictionary length
 */
guint
ncm_var_dict_len (NcmVarDict *vd)
{
  return g_hash_table_size ((GHashTable *) vd);
}

/**
 * ncm_var_dict_keys:
 * @vd: a #NcmVarDict
 *
 * Gets the keys of a #NcmVarDict.
 *
 * Returns: (transfer container): the keys of a #NcmVarDict.
 */
GStrv
ncm_var_dict_keys (NcmVarDict *vd)
{
  return (GStrv) g_hash_table_get_keys_as_array ((GHashTable *) vd, NULL);
}


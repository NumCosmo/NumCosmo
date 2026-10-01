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
 * NcmObjArray:
 *
 * Reference-counted array of #GObject, serializable by #NcmSerialize.
 *
 * A #GPtrArray that holds a reference to each element.
 */
/**
 * NcmObjDictStr:
 *
 * Reference-counted dictionary from strings to #GObject, serializable by
 * #NcmSerialize.
 *
 * A #GHashTable that holds a copy of each key and a reference to each value.
 */
/**
 * NcmObjDictInt:
 *
 * Reference-counted dictionary from integers to #GObject, serializable by
 * #NcmSerialize.
 *
 * A #GHashTable that holds a reference to each value.
 */
/**
 * NcmVarDict:
 *
 * Reference-counted dictionary from strings to #GVariant, serializable by
 * #NcmSerialize.
 *
 * The values are strings, 32-bit integers, doubles, booleans, arrays of these three
 * numeric types, and serialized #GObject or #NcmObjArray; ncm_var_dict_set_variant()
 * lists the #GVariant types. Each getter requires the value under the key to have
 * its type.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_obj_array.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"

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
 * Creates a new empty #NcmObjArray.
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
 * @n: number of elements to preallocate
 *
 * Creates a new empty #NcmObjArray with room for @n elements.
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
 * @oa: a #NcmObjArray
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
 * @oa: a #NcmObjArray
 *
 * Decreases the reference count of @oa by one. At zero, releases the references to the elements.
 */
void
ncm_obj_array_unref (NcmObjArray *oa)
{
  g_ptr_array_unref ((GPtrArray *) oa);
}

/**
 * ncm_obj_array_clear:
 * @oa: a #NcmObjArray
 *
 * If *@oa is not %NULL, decreases its reference count by one and sets *@oa to %NULL.
 */
void
ncm_obj_array_clear (NcmObjArray **oa)
{
  g_clear_pointer (oa, ncm_obj_array_unref);
}

/**
 * ncm_obj_array_add:
 * @oa: a #NcmObjArray
 * @obj: a #GObject
 *
 * Appends @obj to @oa, holding a reference to it.
 */
void
ncm_obj_array_add (NcmObjArray *oa, GObject *obj)
{
  g_assert (obj != NULL);
  g_ptr_array_add ((GPtrArray *) oa, g_object_ref (obj));
}

/**
 * ncm_obj_array_set:
 * @oa: a #NcmObjArray
 * @i: index
 * @obj: a #GObject
 *
 * Replaces the element at @i, which must exist, by @obj, holding a reference to it.
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
 * @i: index
 *
 * Returns: (transfer full): the element at @i.
 */
GObject *
ncm_obj_array_get (NcmObjArray *oa, guint i)
{
  return g_object_ref (ncm_obj_array_peek (oa, i));
}

/**
 * ncm_obj_array_peek:
 * @oa: a #NcmObjArray
 * @i: index
 *
 * Returns: (transfer none): the element at @i.
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
 * Returns: the number of elements of @oa.
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
 * Creates a new empty #NcmObjDictStr.
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
 * @ods: a #NcmObjDictStr
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
 * @ods: a #NcmObjDictStr
 *
 * Decreases the reference count of @ods by one. At zero, releases the keys and the references to the values.
 */
void
ncm_obj_dict_str_unref (NcmObjDictStr *ods)
{
  g_hash_table_unref ((GHashTable *) ods);
}

/**
 * ncm_obj_dict_str_clear:
 * @ods: a #NcmObjDictStr
 *
 * If *@ods is not %NULL, decreases its reference count by one and sets *@ods to %NULL.
 */
void
ncm_obj_dict_str_clear (NcmObjDictStr **ods)
{
  g_clear_pointer (ods, ncm_obj_dict_str_unref);
}

/**
 * ncm_obj_dict_str_add:
 * @ods: a #NcmObjDictStr
 * @key: the key
 * @obj: a #GObject
 *
 * Stores @obj under @key, holding a reference to it. A value already under @key is
 * replaced.
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
 * @ods: a #NcmObjDictStr
 * @key: the key
 * @obj: a #GObject
 *
 * Same as ncm_obj_dict_str_add().
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
 * @key: the key
 *
 * Returns: (transfer full) (nullable): the value under @key, or %NULL.
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
 * @ods: a #NcmObjDictStr
 * @key: the key
 *
 * Returns: (transfer none) (nullable): the value under @key, or %NULL.
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
 * Returns: the number of keys of @ods.
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
 * Returns: (transfer container): the keys of @ods, in no particular order.
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
 * Creates a new empty #NcmObjDictInt.
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
 * @odi: a #NcmObjDictInt
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
 * @odi: a #NcmObjDictInt
 *
 * Decreases the reference count of @odi by one. At zero, releases the references to the values.
 */
void
ncm_obj_dict_int_unref (NcmObjDictInt *odi)
{
  g_hash_table_unref ((GHashTable *) odi);
}

/**
 * ncm_obj_dict_int_clear:
 * @odi: a #NcmObjDictInt
 *
 * If *@odi is not %NULL, decreases its reference count by one and sets *@odi to %NULL.
 */
void
ncm_obj_dict_int_clear (NcmObjDictInt **odi)
{
  g_clear_pointer (odi, ncm_obj_dict_int_unref);
}

/**
 * ncm_obj_dict_int_add:
 * @odi: a #NcmObjDictInt
 * @key: the key
 * @obj: a #GObject
 *
 * Stores @obj under @key, holding a reference to it. A value already under @key is
 * replaced.
 */
void
ncm_obj_dict_int_add (NcmObjDictInt *odi, gint key, GObject *obj)
{
  g_assert (obj != NULL);

  g_hash_table_insert ((GHashTable *) odi, g_memdup2 (&key, sizeof (gint)), g_object_ref (obj));
}

/**
 * ncm_obj_dict_int_set:
 * @odi: a #NcmObjDictInt
 * @key: the key
 * @obj: a #GObject
 *
 * Same as ncm_obj_dict_int_add().
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
 * @key: the key
 *
 * Returns: (transfer full) (nullable): the value under @key, or %NULL.
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
 * @odi: a #NcmObjDictInt
 * @key: the key
 *
 * Returns: (transfer none) (nullable): the value under @key, or %NULL.
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
 * Returns: the number of keys of @odi.
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
 * Returns: (transfer full) (array) (element-type int): the keys of @odi, in no particular order.
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
 * Creates a new empty #NcmVarDict.
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
 * @vd: a #NcmVarDict
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
 * @vd: a #NcmVarDict
 *
 * Decreases the reference count of @vd by one. At zero, releases the keys and values.
 */
void
ncm_var_dict_unref (NcmVarDict *vd)
{
  g_hash_table_unref ((GHashTable *) vd);
}

/**
 * ncm_var_dict_clear:
 * @vd: a #NcmVarDict
 *
 * If *@vd is not %NULL, decreases its reference count by one and sets *@vd to %NULL.
 */
void
ncm_var_dict_clear (NcmVarDict **vd)
{
  g_clear_pointer (vd, ncm_var_dict_unref);
}

/* The value under @key, or NULL */
static GVariant *
ncm_var_dict_peek (NcmVarDict *vd, const gchar *key)
{
  g_assert (key != NULL);

  return g_hash_table_lookup ((GHashTable *) vd, key);
}

/* The value under @key, or NULL; aborts if it is not of @type */
static GVariant *
_ncm_var_dict_peek_type (NcmVarDict *vd, const gchar *key, const GVariantType *type, const gchar *func)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if ((v != NULL) && !g_variant_is_of_type (v, type))
    g_error ("%s: the value under `%s' has type `%s', not `%.*s'.", func, key,
             g_variant_get_type_string (v), (gint) g_variant_type_get_string_length (type),
             g_variant_type_peek_string (type));

  return v;
}

/**
 * ncm_var_dict_set_string:
 * @vd: a #NcmVarDict
 * @key: the key
 * @value: a string
 *
 * Stores @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: an integer
 *
 * Stores @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: a double
 *
 * Stores @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: a boolean
 *
 * Stores @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: (array) (element-type int): an array of integers
 *
 * Stores a copy of @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: (array) (element-type double): an array of doubles
 *
 * Stores a copy of @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: (array) (element-type boolean): an array of booleans
 *
 * Stores a copy of @value under @key. A value already under @key is replaced.
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
 * @key: the key
 * @value: a #GVariant
 *
 * Stores @value under @key. A value already under @key is replaced. Aborts unless
 * @value has one of the types
 *
 * - G_VARIANT_TYPE_STRING
 * - G_VARIANT_TYPE_INT32
 * - G_VARIANT_TYPE_DOUBLE
 * - G_VARIANT_TYPE_BOOLEAN
 * - "ai", "ad" and "ab"
 * - #NCM_SERIALIZE_OBJECT_TYPE, as stored by ncm_var_dict_set_object()
 * - #NCM_SERIALIZE_OBJECT_ARRAY_TYPE, as stored by ncm_var_dict_set_object_array()
 *
 * The last two are accepted because ncm_serialize_var_dict_from_variant() restores
 * every entry through this function.
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
    G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE),
    G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_ARRAY_TYPE),
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
 * ncm_var_dict_set_object:
 * @vd: a #NcmVarDict
 * @key: the key
 * @ser: a #NcmSerialize
 * @obj: a #GObject
 *
 * Serializes @obj with @ser and stores it under @key. A value already under @key is
 * replaced. @ser carries the named instances and the autosave and autoname settings,
 * so objects shared across a larger serialization stay shared.
 */
void
ncm_var_dict_set_object (NcmVarDict *vd, const gchar *key, NcmSerialize *ser, GObject *obj)
{
  GVariant *var;

  g_assert (key != NULL);
  g_assert (obj != NULL);

  var = ncm_serialize_to_variant (ser, obj);
  g_hash_table_insert ((GHashTable *) vd, g_strdup (key), g_variant_ref_sink (var));
}

/**
 * ncm_var_dict_set_object_array:
 * @vd: a #NcmVarDict
 * @key: the key
 * @ser: a #NcmSerialize
 * @oa: a #NcmObjArray
 *
 * Serializes @oa with @ser and stores it under @key, as ncm_var_dict_set_object().
 */
void
ncm_var_dict_set_object_array (NcmVarDict *vd, const gchar *key, NcmSerialize *ser, NcmObjArray *oa)
{
  GVariant *var;

  g_assert (key != NULL);
  g_assert (oa != NULL);

  var = ncm_serialize_array_to_variant (ser, oa);
  g_hash_table_insert ((GHashTable *) vd, g_strdup (key), g_variant_ref_sink (var));
}

/**
 * ncm_var_dict_has_key:
 * @vd: a #NcmVarDict
 * @key: the key
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_has_key (NcmVarDict *vd, const gchar *key)
{
  return g_hash_table_contains ((GHashTable *) vd, key);
}

/**
 * ncm_var_dict_get_string:
 * @vd: a #NcmVarDict
 * @key: the key
 * @value: (out) (transfer full): the string
 *
 * Gets the string under @key. @value is set only when @key is present. Aborts if the value
 * under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_string (NcmVarDict *vd, const gchar *key, gchar **value)
{
  GVariant *v = _ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE_STRING, "ncm_var_dict_get_string");

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
 * @key: the key
 * @value: (out): the integer
 *
 * Gets the integer under @key. @value is set only when @key is present. Aborts if the value
 * under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_int (NcmVarDict *vd, const gchar *key, gint *value)
{
  GVariant *v = _ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE_INT32, "ncm_var_dict_get_int");

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
 * @key: the key
 * @value: (out): the double
 *
 * Gets the double under @key; an integer is converted. @value is set only when @key is
 * present. Aborts if the value under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_double (NcmVarDict *vd, const gchar *key, gdouble *value)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v == NULL)
    return FALSE;

  if (g_variant_is_of_type (v, G_VARIANT_TYPE_INT32))
    *value = g_variant_get_int32 (v);
  else
    *value = g_variant_get_double (_ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE_DOUBLE, "ncm_var_dict_get_double"));

  return TRUE;
}

/**
 * ncm_var_dict_get_boolean:
 * @vd: a #NcmVarDict
 * @key: the key
 * @value: (out): the boolean
 *
 * Gets the boolean under @key. @value is set only when @key is present. Aborts if the value
 * under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_boolean (NcmVarDict *vd, const gchar *key, gboolean *value)
{
  GVariant *v = _ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE_BOOLEAN, "ncm_var_dict_get_boolean");

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
 * @key: the key
 * @value: (out) (transfer full) (element-type int): the array of integers
 *
 * Copies the array of integers under @key to a new #GArray. @value is set only when @key
 * is present. Aborts if the value under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_int_array (NcmVarDict *vd, const gchar *key, GArray **value)
{
  GVariant *v = _ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE ("ai"), "ncm_var_dict_get_int_array");

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
 * @key: the key
 * @value: (out) (transfer full) (element-type double): the array of doubles
 *
 * Copies the array of doubles under @key to a new #GArray. @value is set only when @key
 * is present. Aborts if the value under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_double_array (NcmVarDict *vd, const gchar *key, GArray **value)
{
  GVariant *v = _ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE ("ad"), "ncm_var_dict_get_double_array");

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
 * @key: the key
 * @value: (out) (transfer full) (element-type boolean): the array of booleans
 *
 * Copies the array of booleans under @key to a new #GArray. @value is set only when @key
 * is present. Aborts if the value under @key has another type.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_boolean_array (NcmVarDict *vd, const gchar *key, GArray **value)
{
  GVariant *v = _ncm_var_dict_peek_type (vd, key, G_VARIANT_TYPE ("ab"), "ncm_var_dict_get_boolean_array");

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
 * @key: the key
 * @value: (out) (transfer full): the #GVariant
 *
 * Gets the #GVariant under @key. @value is set only when @key is present.
 *
 * Returns: whether @key is present.
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
 * ncm_var_dict_get_object:
 * @vd: a #NcmVarDict
 * @key: the key
 * @ser: a #NcmSerialize
 * @obj: (out) (transfer full): the object
 *
 * Deserializes the object under @key with @ser, see ncm_var_dict_set_object(). @obj
 * is set only when @key is present.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_object (NcmVarDict *vd, const gchar *key, NcmSerialize *ser, GObject **obj)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *obj = ncm_serialize_from_variant (ser, v);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_get_object_array:
 * @vd: a #NcmVarDict
 * @key: the key
 * @ser: a #NcmSerialize
 * @oa: (out) (transfer full): the #NcmObjArray
 *
 * Deserializes the #NcmObjArray under @key with @ser, see ncm_var_dict_set_object().
 * @oa is set only when @key is present.
 *
 * Returns: whether @key is present.
 */
gboolean
ncm_var_dict_get_object_array (NcmVarDict *vd, const gchar *key, NcmSerialize *ser, NcmObjArray **oa)
{
  GVariant *v = ncm_var_dict_peek (vd, key);

  if (v != NULL)
  {
    *oa = ncm_serialize_array_from_variant (ser, v);

    return TRUE;
  }

  return FALSE;
}

/**
 * ncm_var_dict_len:
 * @vd: a #NcmVarDict
 *
 * Returns: the number of keys of @vd.
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
 * Returns: (transfer container): the keys of @vd, in no particular order.
 */
GStrv
ncm_var_dict_keys (NcmVarDict *vd)
{
  return (GStrv) g_hash_table_get_keys_as_array ((GHashTable *) vd, NULL);
}


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
  GHashTable *od = g_hash_table_new_full (g_str_hash, g_str_equal, g_free, g_object_unref);

  return (NcmObjDictStr *) od;
}

/**
 * ncm_obj_dict_str_ref:
 * @od: a #NcmObjDictStr.
 *
 * Increases the reference count of @od by one.
 *
 * Returns: (transfer full): @od.
 */
NcmObjDictStr *
ncm_obj_dict_str_ref (NcmObjDictStr *od)
{
  return (NcmObjDictStr *) g_hash_table_ref ((GHashTable *) od);
}

/**
 * ncm_obj_dict_str_unref:
 * @od: a #NcmObjDictStr.
 *
 * Decreases the reference count of @od by one. If the reference count
 * reaches zero, all objects in the dictionary are unreferenced.
 *
 */
void
ncm_obj_dict_str_unref (NcmObjDictStr *od)
{
  g_hash_table_unref ((GHashTable *) od);
}

/**
 * ncm_obj_dict_str_clear:
 * @od: a pointer to a #NcmObjDictStr.
 *
 * If *@od is not %NULL, unreferences it and sets *@od to %NULL.
 *
 */
void
ncm_obj_dict_str_clear (NcmObjDictStr **od)
{
  g_clear_pointer (od, ncm_obj_dict_str_unref);
}

/**
 * ncm_obj_dict_str_add:
 * @od: a #NcmObjDictStr.
 * @key: a string.
 * @obj: a #GObject.
 *
 * Adds a #GObject to a #NcmObjDictStr.
 *
 */
void
ncm_obj_dict_str_add (NcmObjDictStr *od, const gchar *key, GObject *obj)
{
  g_assert (key != NULL);
  g_assert (obj != NULL);

  g_hash_table_insert ((GHashTable *) od, g_strdup (key), g_object_ref (obj));
}

/**
 * ncm_obj_dict_str_set:
 * @od: a #NcmObjDictStr.
 * @key: a string.
 * @obj: a #GObject.
 *
 * Sets a #GObject to a #NcmObjDictStr. If there is already a #GObject
 * with key @key, it is unreferenced.
 *
 */
void
ncm_obj_dict_str_set (NcmObjDictStr *od, const gchar *key, GObject *obj)
{
  g_assert (key != NULL);
  g_assert (obj != NULL);

  if (obj != g_hash_table_lookup ((GHashTable *) od, key))
    g_hash_table_replace ((GHashTable *) od, g_strdup (key), g_object_ref (obj));
}

/**
 * ncm_obj_dict_str_get:
 * @od: a #NcmObjDictStr
 * @key: a string.
 *
 * Gets a #GObject from a #NcmObjDictStr with key @key.
 *
 * Returns: (transfer full): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_str_get (NcmObjDictStr *od, const gchar *key)
{
  GObject *obj = ncm_obj_dict_str_peek (od, key);

  if (obj != NULL)
    return g_object_ref (obj);

  return NULL;
}

/**
 * ncm_obj_dict_str_peek:
 * @od: a #NcmObjDictStr.
 * @key: a string.
 *
 * Peeks a #GObject from a #NcmObjDictStr with key @key without increasing its reference
 * count.
 *
 * Returns: (transfer none): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_str_peek (NcmObjDictStr *od, const gchar *key)
{
  g_assert (key != NULL);

  return g_hash_table_lookup ((GHashTable *) od, key);
}

/**
 * ncm_obj_dict_str_len:
 * @od: a #NcmObjDictStr
 *
 * Gets the length of a #NcmObjDictStr.
 *
 * Returns: dictionary length
 */
guint
ncm_obj_dict_str_len (NcmObjDictStr *od)
{
  return g_hash_table_size ((GHashTable *) od);
}

/**
 * ncm_obj_dict_str_keys:
 * @od: a #NcmObjDictStr
 *
 * Gets the keys of a #NcmObjDictStr.
 *
 * Returns: (transfer container): the keys of a #NcmObjDictStr.
 */
GStrv
ncm_obj_dict_str_keys (NcmObjDictStr *od)
{
  return (GStrv) g_hash_table_get_keys_as_array ((GHashTable *) od, NULL);
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
  GHashTable *od = g_hash_table_new_full (g_int_hash, g_int_equal, g_free, g_object_unref);

  return (NcmObjDictInt *) od;
}

/**
 * ncm_obj_dict_int_ref:
 * @od: a #NcmObjDictInt.
 *
 * Increases the reference count of @od by one.
 *
 * Returns: (transfer full): @od.
 */
NcmObjDictInt *
ncm_obj_dict_int_ref (NcmObjDictInt *od)
{
  return (NcmObjDictInt *) g_hash_table_ref ((GHashTable *) od);
}

/**
 * ncm_obj_dict_int_unref:
 * @od: a #NcmObjDictInt.
 *
 * Decreases the reference count of @od by one. If the reference count
 * reaches zero, all objects in the dictionary are unreferenced.
 *
 */
void
ncm_obj_dict_int_unref (NcmObjDictInt *od)
{
  g_hash_table_unref ((GHashTable *) od);
}

/**
 * ncm_obj_dict_int_clear:
 * @od: a pointer to a #NcmObjDictInt.
 *
 * If *@od is not %NULL, unreferences it and sets *@od to %NULL.
 *
 */
void
ncm_obj_dict_int_clear (NcmObjDictInt **od)
{
  g_clear_pointer (od, ncm_obj_dict_int_unref);
}

/**
 * ncm_obj_dict_int_add:
 * @od: a #NcmObjDictInt.
 * @key: an integer.
 * @obj: a #GObject.
 *
 * Adds a #GObject to a #NcmObjDictInt.
 *
 */
void
ncm_obj_dict_int_add (NcmObjDictInt *od, gint key, GObject *obj)
{
  g_assert (obj != NULL);

  g_hash_table_insert ((GHashTable *) od, g_memdup2 (&key, sizeof (gint)), g_object_ref (obj));
}

/**
 * ncm_obj_dict_int_set:
 * @od: a #NcmObjDictInt.
 * @key: an integer.
 * @obj: a #GObject.
 *
 * Sets a #GObject to a #NcmObjDictInt. If there is already a #GObject
 * with key @key, it is unreferenced.
 *
 */
void
ncm_obj_dict_int_set (NcmObjDictInt *od, gint key, GObject *obj)
{
  g_assert (obj != NULL);

  if (obj != g_hash_table_lookup ((GHashTable *) od, &key))
    g_hash_table_replace ((GHashTable *) od, g_memdup2 (&key, sizeof (gint)), g_object_ref (obj));
}

/**
 * ncm_obj_dict_int_get:
 * @od: a #NcmObjDictInt
 * @key: an integer.
 *
 * Gets a #GObject from a #NcmObjDictInt with key @key.
 *
 * Returns: (transfer full): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_int_get (NcmObjDictInt *od, gint key)
{
  GObject *obj = ncm_obj_dict_int_peek (od, key);

  if (obj != NULL)
    return g_object_ref (obj);

  return NULL;
}

/**
 * ncm_obj_dict_int_peek:
 * @od: a #NcmObjDictInt.
 * @key: an integer.
 *
 * Peeks a #GObject from a #NcmObjDictInt with key @key without increasing its reference
 * count.
 *
 * Returns: (transfer none): the #GObject with key @key.
 */
GObject *
ncm_obj_dict_int_peek (NcmObjDictInt *od, gint key)
{
  return g_hash_table_lookup ((GHashTable *) od, &key);
}

/**
 * ncm_obj_dict_int_len:
 * @od: a #NcmObjDictInt
 *
 * Gets the length of a #NcmObjDictInt.
 *
 * Returns: dictionary length
 */
guint
ncm_obj_dict_int_len (NcmObjDictInt *od)
{
  return g_hash_table_size ((GHashTable *) od);
}

/**
 * ncm_obj_dict_int_keys:
 * @od: a #NcmObjDictInt
 *
 * Gets the keys of a #NcmObjDictInt.
 *
 * Returns: (transfer full) (array) (element-type int): the keys of a #NcmObjDictInt.
 */
GArray *
ncm_obj_dict_int_keys (NcmObjDictInt *od)
{
  GArray *keys = g_array_new (FALSE, FALSE, sizeof (gint));
  GHashTableIter iter;
  gint *key;

  g_hash_table_iter_init (&iter, (GHashTable *) od);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, NULL))
    g_array_append_val (keys, *key);

  return keys;
}


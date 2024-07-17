/***************************************************************************
 *            ncm_serialize.c
 *
 *  Mon August 26 13:38:17 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_serialize.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/**
 * SECTION:ncm_serialize
 * @title: NcmSerialize
 * @short_description: Serialization, deserialization and duplication object.
 *
 * This object provides serialization, deserialization and duplication of objects.
 * The serialization process is based on the #GObject object system.
 *
 * Serialization is the process of converting an object into a stream of bytes to store the object
 * or transmit it to memory, a database, or a file. Its main purpose is to save the state of an object
 * in order to be able to recreate it when needed. The reverse process is called deserialization.
 *
 * One support for serialized data #GVariant.
 * The #GVariant is a type-safe, reference counted, immutable, and memory-efficient container for arbitrary data.
 * It is a generic container that can hold any type of data, including basic types such as integers and floating
 * point numbers, strings, and byte arrays, as well as more complex types such as tuples, dictionaries, and variants.
 * The #GVariant type system is designed to be extensible, so that new types can be added in the future.
 * A serialized #GVariant object can be stored in binary or text format.
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_obj_array.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_vector.h"
#include "math/ncm_matrix.h"
#include "ncm_enum_types.h"

#include <gio/gio.h>
#include <float.h>
#include <errno.h>

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_LIBFYAML
#include <libfyaml.h>
#endif /* HAVE_LIBFYAML */
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_OPTS,
};

struct _NcmSerialize
{
  /*< private >*/
  GObject parent_instance;
  GHashTable *name_ptr;
  GHashTable *ptr_name;
  GHashTable *saved_ptr_name;
  GHashTable *saved_name_ser;
  GRegex *is_named_regex;
  GRegex *parse_obj_regex;
  NcmSerializeOpt opts;
  guint autosave_count;
};

G_DEFINE_TYPE (NcmSerialize, ncm_serialize, G_TYPE_OBJECT)

static void
ncm_serialize_init (NcmSerialize *ser)
{
  GError *error = NULL;

  ser->name_ptr = g_hash_table_new_full (&g_str_hash, &g_str_equal, &g_free,
                                         &g_object_unref);
  ser->ptr_name = g_hash_table_new_full (&g_direct_hash, &g_direct_equal, &g_object_unref,
                                         &g_free);

  ser->saved_ptr_name = g_hash_table_new_full (&g_direct_hash, &g_direct_equal, &g_object_unref,
                                               &g_free);
  ser->saved_name_ser = g_hash_table_new_full (&g_str_hash, &g_str_equal, &g_free,
                                               (GDestroyNotify) & g_variant_unref);

  ser->is_named_regex  = g_regex_new ("^\\s*([A-Za-z][A-Za-z0-9\\+\\-\\_]+)\\s*\\[([A-Za-z0-9\\:]+)\\]\\s*$", 0, 0, &error);
  ser->parse_obj_regex = g_regex_new ("^\\s*([A-Za-z][A-Za-z0-9\\+\\-\\_]+\\s*(?:\\[[A-Za-z0-9\\:]+\\])?)\\s*([\\{]?.*[\\}]?)\\s*$", 0, 0, &error);
  ser->autosave_count  = 0;
}

static void
_ncm_serialize_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSerialize *ser = NCM_SERIALIZE (object);

  g_return_if_fail (NCM_IS_SERIALIZE (object));

  switch (prop_id)
  {
    case PROP_OPTS:
      g_value_set_flags (value, ser->opts);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_serialize_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSerialize *ser = NCM_SERIALIZE (object);

  g_return_if_fail (NCM_IS_SERIALIZE (object));

  switch (prop_id)
  {
    case PROP_OPTS:
      ser->opts = g_value_get_flags (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_serialize_finalize (GObject *object)
{
  NcmSerialize *ser = NCM_SERIALIZE (object);

  g_hash_table_remove_all (ser->name_ptr);
  g_hash_table_remove_all (ser->ptr_name);

  g_hash_table_remove_all (ser->saved_ptr_name);
  g_hash_table_remove_all (ser->saved_name_ser);

  g_hash_table_unref (ser->name_ptr);
  g_hash_table_unref (ser->ptr_name);

  g_hash_table_unref (ser->saved_ptr_name);
  g_hash_table_unref (ser->saved_name_ser);

  g_regex_unref (ser->is_named_regex);
  g_regex_unref (ser->parse_obj_regex);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_serialize_parent_class)->finalize (object);
}

static void
ncm_serialize_class_init (NcmSerializeClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_serialize_set_property;
  object_class->get_property = &_ncm_serialize_get_property;
  object_class->finalize     = &_ncm_serialize_finalize;

  /**
   * NcmSerialize:options:
   *
   * Serialization options.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OPTS,
                                   g_param_spec_flags ("options",
                                                       NULL,
                                                       "Serialization options",
                                                       NCM_TYPE_SERIALIZE_OPT, NCM_SERIALIZE_OPT_NONE,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_serialize_new:
 * @sopt: a set of options from #NcmSerializeOpt.
 *
 * Creates a new #NcmSerialize object.
 *
 * Returns: a new #NcmSerialize.
 */
NcmSerialize *
ncm_serialize_new (NcmSerializeOpt sopt)
{
  NcmSerialize *ser = g_object_new (NCM_TYPE_SERIALIZE,
                                    "options", sopt,
                                    NULL);

  return ser;
}

/**
 * ncm_serialize_ref:
 * @ser: a #NcmSerialize.
 *
 * Increases the reference count of @ser by one.
 *
 * Returns: (transfer full): @ser.
 */
NcmSerialize *
ncm_serialize_ref (NcmSerialize *ser)
{
  return g_object_ref (ser);
}

/**
 * ncm_serialize_free:
 * @ser: a #NcmSerialize.
 *
 * Decreases the reference count of @ser by one.
 *
 */
void
ncm_serialize_free (NcmSerialize *ser)
{
  g_object_unref (ser);
}

/**
 * ncm_serialize_unref:
 * @ser: a #NcmSerialize.
 *
 * Same as ncm_serialize_free().
 *
 */
void
ncm_serialize_unref (NcmSerialize *ser)
{
  ncm_serialize_free (ser);
}

/**
 * ncm_serialize_clear:
 * @ser: a #NcmSerialize.
 *
 * Decreases the reference count of *@ser by one, and sets *@ser to NULL.
 *
 */
void
ncm_serialize_clear (NcmSerialize **ser)
{
  g_clear_object (ser);
}

static gboolean
_ncm_serialize_reset_remove_autosaved_name_key (gpointer key, gpointer value, gpointer user_data)
{
  gchar *name = (gchar *) key;

  return g_regex_match_simple ("^S[\\d]+$", name, 0, 0);
}

static gboolean
_ncm_serialize_reset_remove_autosaved_name_val (gpointer key, gpointer value, gpointer user_data)
{
  gchar *name = (gchar *) value;

  return g_regex_match_simple ("^S[\\d]+$", name, 0, 0);
}

/**
 * ncm_serialize_reset:
 * @ser: a #NcmSerialize
 * @autosave_only: a boolean
 *
 * Releases all objects in @ser and erase all serialized
 * objects. If @autosave_only is TRUE it will release only
 * autosaved objects.
 *
 */
void
ncm_serialize_reset (NcmSerialize *ser, gboolean autosave_only)
{
  ncm_serialize_clear_instances (ser, autosave_only);

  if (autosave_only)
  {
    g_hash_table_foreach_remove (ser->saved_ptr_name,
                                 &_ncm_serialize_reset_remove_autosaved_name_val,
                                 NULL);
    g_hash_table_foreach_remove (ser->saved_name_ser,
                                 &_ncm_serialize_reset_remove_autosaved_name_key,
                                 NULL);
  }
  else
  {
    g_hash_table_remove_all (ser->saved_ptr_name);
    g_hash_table_remove_all (ser->saved_name_ser);
  }

  ser->autosave_count = 0;
}

/**
 * ncm_serialize_clear_instances:
 * @ser: a #NcmSerialize
 * @autosave_only: a boolean
 *
 * Releases all objects in @ser. If @autosave_only is TRUE
 * it will release only autosaved objects.
 *
 */
void
ncm_serialize_clear_instances (NcmSerialize *ser, gboolean autosave_only)
{
  if (autosave_only)
  {
    g_hash_table_foreach_remove (ser->ptr_name,
                                 &_ncm_serialize_reset_remove_autosaved_name_val,
                                 NULL);
    g_hash_table_foreach_remove (ser->name_ptr,
                                 _ncm_serialize_reset_remove_autosaved_name_key,
                                 NULL);
  }
  else
  {
    g_hash_table_remove_all (ser->name_ptr);
    g_hash_table_remove_all (ser->ptr_name);
  }
}

/**
 * ncm_serialize_log_stats:
 * @ser: a #NcmSerialize.
 *
 * Releases all objects in @ser.
 *
 */
void
ncm_serialize_log_stats (NcmSerialize *ser)
{
  g_message ("# NcmSerialize: autosaved object %u\n", ser->autosave_count);
  g_message ("# NcmSerialize: name_ptr         %u\n", g_hash_table_size (ser->name_ptr));
  g_message ("# NcmSerialize: ptr_name         %u\n", g_hash_table_size (ser->ptr_name));
  g_message ("# NcmSerialize: saved_ptr_name   %u\n", g_hash_table_size (ser->saved_ptr_name));
  g_message ("# NcmSerialize: saved_name_ser   %u\n", g_hash_table_size (ser->saved_name_ser));
}

/**
 * ncm_serialize_contain_instance:
 * @ser: a #NcmSerialize.
 * @obj: (type GObject): a #GObject.
 *
 * Checks if the #GObject instance @obj is contained in @ser.
 *
 * Returns: if @obj is already in @ser.
 */
gboolean
ncm_serialize_contain_instance (NcmSerialize *ser, gpointer obj)
{
  g_assert (G_IS_OBJECT (obj));

  return g_hash_table_lookup_extended (ser->ptr_name, obj, NULL, NULL);
}

/**
 * ncm_serialize_contain_name:
 * @ser: a #NcmSerialize.
 * @name: an instance name.
 *
 * Checks if there is an instance named @name in @ser.
 *
 * Returns: if there is instance named @name in @ser.
 */
gboolean
ncm_serialize_contain_name (NcmSerialize *ser, const gchar *name)
{
  g_assert (name != NULL);

  return g_hash_table_lookup_extended (ser->name_ptr, name, NULL, NULL);
}

/**
 * ncm_serialize_count_instances:
 * @ser: a #NcmSerialize.
 *
 * Counts the number of instances registered in @ser.
 *
 * Returns: the number of instances in @ser.
 */
guint
ncm_serialize_count_instances (NcmSerialize *ser)
{
  return g_hash_table_size (ser->name_ptr);
}

/**
 * ncm_serialize_count_saved_serializations:
 * @ser: a #NcmSerialize.
 *
 * Counts the number of instances registered in @ser.
 *
 * Returns: the number of instances in @ser.
 */
guint
ncm_serialize_count_saved_serializations (NcmSerialize *ser)
{
  return g_hash_table_size (ser->saved_name_ser);
}

/**
 * ncm_serialize_peek_by_name:
 * @ser: a #NcmSerialize.
 * @name: an instance name.
 *
 * Peeks the instance @name or null if there isn't a instance named @name.
 *
 * Returns: (transfer none) (type GObject): Gets the instance named @name or NULL.
 */
gpointer
ncm_serialize_peek_by_name (NcmSerialize *ser, const gchar *name)
{
  g_assert (name != NULL);

  return g_hash_table_lookup (ser->name_ptr, name);
}

/**
 * ncm_serialize_get_by_name:
 * @ser: a #NcmSerialize.
 * @name: an instance name.
 *
 * Gets a new reference for the instance @name or null if there isn't a instance named @name.
 *
 * Returns: (transfer full) (type GObject): Gets the instance named @name or NULL.
 */
gpointer
ncm_serialize_get_by_name (NcmSerialize *ser, const gchar *name)
{
  g_assert (name != NULL);
  {
    gpointer obj = NULL;

    if (g_hash_table_lookup_extended (ser->name_ptr, name, NULL, &obj))
      return g_object_ref (obj);
    else
      return NULL;
  }
}

/**
 * ncm_serialize_peek_name:
 * @ser: a #NcmSerialize.
 * @obj: (type GObject): a #GObject.
 *
 * Gets the named associated to the instance @obj, it is an error to call this function
 * when the @obj is not contained in @ser.
 *
 * Returns: (transfer none): the name of @obj.
 */
gchar *
ncm_serialize_peek_name (NcmSerialize *ser, gpointer obj)
{
  g_assert (G_IS_OBJECT (obj));

  if (!ncm_serialize_contain_instance (ser, obj))
  {
    g_error ("ncm_serialize_peek_name: Cannot peek name of object %p, it is not a named instance.", obj);

    return NULL;
  }
  else
  {
    return g_hash_table_lookup (ser->ptr_name, obj);
  }
}

/**
 * ncm_serialize_set:
 * @ser: a #NcmSerialize.
 * @obj: (type GObject): a #GObject.
 * @name: the @obj name.
 * @overwrite: whether to overwrite if there is already an object named @name.
 *
 * Adds the object @obj to @ser using @name.
 *
 */
void
ncm_serialize_set (NcmSerialize *ser, gpointer obj, const gchar *name, gboolean overwrite)
{
  g_assert (G_IS_OBJECT (obj));
  g_assert (name != NULL);
  {
    gboolean ni_has_key = ncm_serialize_contain_name (ser, name);

    if (!overwrite && ni_has_key)
    {
      gpointer ni_obj = ncm_serialize_get_by_name (ser, name);

      if (ni_obj != obj)
      {
        g_error ("ncm_serialize_set: named instance already present, set overwrite to true if you want it replaced.");
      }
      else
      {
        g_object_unref (ni_obj);

        return;
      }
    }

    g_hash_table_insert (ser->name_ptr,
                         g_strdup (name), g_object_ref (obj));
    g_hash_table_insert (ser->ptr_name,
                         g_object_ref (obj), g_strdup (name));
  }
}

/**
 * ncm_serialize_unset:
 * @ser: a #NcmSerialize.
 * @obj: (type GObject): a #GObject.
 *
 * Removes the object @obj to @ser using @name, it does nothing
 * if the instance @obj is not present in @ser.
 *
 */
void
ncm_serialize_unset (NcmSerialize *ser, gpointer obj)
{
  if (ncm_serialize_contain_instance (ser, obj))
  {
    const gchar *name = ncm_serialize_peek_name (ser, obj);

    g_hash_table_remove (ser->name_ptr, name);
    g_hash_table_remove (ser->ptr_name, obj);
  }
}

/**
 * ncm_serialize_remove_ser:
 * @ser: a #NcmSerialize.
 * @obj: (type GObject): a #GObject.
 *
 * Removes the object @obj to @ser using @name, it does nothing
 * if the instance @obj is not present in @ser.
 *
 */
void
ncm_serialize_remove_ser (NcmSerialize *ser, gpointer obj)
{
  gchar *saved_name = NULL;

  if (g_hash_table_lookup_extended (ser->saved_ptr_name, obj, NULL, (gpointer *) &saved_name))
  {
    g_hash_table_remove (ser->saved_ptr_name, obj);
    g_hash_table_remove (ser->saved_name_ser, saved_name);
  }
}

static void
_ncm_serialize_save_ser (NcmSerialize *ser, gchar *name, gpointer obj, GVariant *ser_var)
{
  g_assert (G_IS_OBJECT (obj));
  g_assert (name != NULL);
  g_assert (ser != NULL);

  if (g_hash_table_lookup_extended (ser->saved_name_ser, name, NULL, NULL))
    g_error ("_ncm_serialize_save_ser: instance already saved.");

  if (g_hash_table_lookup_extended (ser->saved_ptr_name, obj, NULL, NULL))
    g_error ("_ncm_serialize_save_ser: instance already saved.");

  /*printf ("Saving: ``%s'' <=> %s\n", name, g_variant_print (ser_var, TRUE)); */

  g_hash_table_insert (ser->saved_ptr_name,
                       g_object_ref (obj), g_strdup (name));
  g_hash_table_insert (ser->saved_name_ser,
                       g_strdup (name), g_variant_ref_sink (ser_var));
}

/**
 * ncm_serialize_is_named:
 * @ser: a #NcmSerialize.
 * @serobj: serialized object.
 * @name: (allow-none) (out) (transfer full): object name.
 *
 * Checks if @serobj is a named serialized object, if so sets its name in @name
 * and returns TRUE.
 *
 * Returns: whether @serobj is a named serialized object.
 */
gboolean
ncm_serialize_is_named (NcmSerialize *ser, const gchar *serobj, gchar **name)
{
  GError *error          = NULL;
  GVariant *var_obj      = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), serobj, NULL, NULL, &error);
  GMatchInfo *match_info = NULL;
  gchar *obj_name_str    = NULL;
  gboolean is_named      = FALSE;

  g_clear_error (&error);

  if (var_obj != NULL)
  {
    GVariant *obj_name = g_variant_get_child_value (var_obj, 0);

    obj_name_str = g_variant_dup_string (obj_name, NULL);

    g_variant_unref (obj_name);
    g_variant_unref (var_obj);
  }
  else if (g_regex_match (ser->parse_obj_regex, serobj, 0, &match_info))
  {
    obj_name_str = g_match_info_fetch (match_info, 1);
    g_match_info_free (match_info);
  }
  else
  {
    g_error ("ncm_serialize_is_named: cannot identify object (%s) in string '%s'.",
             NCM_SERIALIZE_OBJECT_TYPE,
             serobj);
  }

  *name = NULL;

  if (g_regex_match (ser->is_named_regex, obj_name_str, 0, &match_info))
  {
    if (name != NULL)
      *name = g_match_info_fetch (match_info, 2);

    is_named = TRUE;
  }

  g_match_info_free (match_info);
  g_free (obj_name_str);

  return is_named;
}

/**
 * ncm_serialize_set_property:
 * @ser: a #NcmSerialize.
 * @obj: a #GObject.
 * @prop_str: a string containing the parameters to set.
 *
 * Deserialize the set of object properties in @params and sets the @obj.
 *
 */
void
ncm_serialize_set_property (NcmSerialize *ser, GObject *obj, const gchar *prop_str)
{
  GError *error    = NULL;
  GVariant *params = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE),
                                      prop_str, NULL, NULL, &error);

  if (params == NULL)
    g_error ("ncm_serialize_set_property: cannot parse prop_str %s: %s.",
             prop_str, error->message);

  g_assert (g_variant_is_of_type (params, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE)));

  {
    GVariantIter *p_iter = g_variant_iter_new (params);
    GVariant *var        = NULL;
    guint i              = 0;

    while ((var = g_variant_iter_next_value (p_iter)))
    {
      GVariant *var_key = g_variant_get_child_value (var, 0);
      GVariant *var_val = g_variant_get_child_value (var, 1);
      GVariant *val     = g_variant_get_variant (var_val);
      const gchar *name = g_variant_get_string (var_key, NULL);
      GValue value      = G_VALUE_INIT;

      if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)))
      {
        GVariant *nest_obj_key    = g_variant_get_child_value (val, 0);
        GVariant *nest_obj_params = g_variant_get_child_value (val, 1);
        GObject *nest_obj         =
          ncm_serialize_from_name_params (ser,
                                          g_variant_get_string (nest_obj_key, NULL),
                                          nest_obj_params);

        g_value_init (&value, G_TYPE_OBJECT);
        g_value_take_object (&value, nest_obj);
        g_variant_unref (nest_obj_key);
        g_variant_unref (nest_obj_params);
      }
      else
      {
        g_dbus_gvariant_to_gvalue (val, &value);
      }

      g_object_set_property (obj, name, &value);

      i++;
      g_variant_unref (var_key);
      g_variant_unref (val);
      g_variant_unref (var_val);
      g_variant_unref (var);
      g_value_unset (&value);
    }

    g_variant_iter_free (p_iter);
  }
}

/**
 * ncm_serialize_set_property_from_key_file:
 * @ser: a #NcmSerialize.
 * @obj: a #GObject.
 * @prop_file: a GKeyFile file containing the parameters to set.
 *
 * Deserializes the set of object properties in @prop_file and sets the @obj.
 *
 */
void
ncm_serialize_set_property_from_key_file (NcmSerialize *ser, GObject *obj, const gchar *prop_file)
{
  GKeyFile *key_file = g_key_file_new ();
  GError *error      = NULL;
  GString *prop_ser  = g_string_sized_new (200);

  if (!g_key_file_load_from_file (key_file, prop_file, G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS, &error))
  {
    gchar *contents;
    gsize length = 0;

    /* If it does not posses a group add one. */
    g_clear_error (&error);

    if (!g_file_get_contents (prop_file, &contents, &length, &error))
    {
      g_error ("ncm_serialize_set_property_from_key_file: Invalid mset configuration file: %s %s", prop_file, error->message);

      return;
    }
    else
    {
      gchar *gcontents = g_strconcat ("[Precision Parameters]\n\n\n", contents, NULL);

      g_free (contents);

      if (!g_key_file_load_from_data (key_file, gcontents, strlen (gcontents), G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS, &error))
      {
        g_error ("ncm_serialize_set_property_from_key_file: Invalid mset configuration file: %s %s", prop_file, error->message);

        return;
      }

      g_free (gcontents);
    }
  }

  g_string_append_printf (prop_ser, "@a{sv} {");
  {
    gsize nkeys  = 0;
    gchar **keys = g_key_file_get_keys (key_file, "Precision Parameters", &nkeys, &error);
    guint i;

    if (error != NULL)
      g_error ("ncm_serialize_set_property_from_key_file: %s", error->message);

    for (i = 0; i < nkeys; i++)
    {
      gchar *propval = g_key_file_get_value (key_file, "Precision Parameters", keys[i], &error);

      if (error != NULL)
        g_error ("ncm_serialize_set_property_from_key_file: %s", error->message);

      g_string_append_printf (prop_ser, "\'%s\':<%s>", keys[i], propval);
      g_free (propval);

      if (i + 1 != nkeys)
        g_string_append (prop_ser, ", ");
    }

    g_string_append (prop_ser, "}");
    g_strfreev (keys);
  }

  ncm_serialize_set_property (ser, obj, prop_ser->str);
  g_key_file_unref (key_file);
  g_string_free (prop_ser, TRUE);
}

/**
 * ncm_serialize_from_variant:
 * @ser: a #NcmSerialize.
 * @var_obj: A #GVariant containing the serialized version of the object.
 *
 * Deserialize and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_variant (NcmSerialize *ser, GVariant *var_obj)
{
  g_assert (var_obj != NULL);
  g_assert (g_variant_is_of_type (var_obj, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)));

  {
    GVariant *obj_name = g_variant_get_child_value (var_obj, 0);
    GVariant *params   = g_variant_get_child_value (var_obj, 1);
    GObject *obj       = ncm_serialize_from_name_params (ser, g_variant_get_string (obj_name, NULL), params);

    g_variant_unref (obj_name);
    g_variant_unref (params);

    return obj;
  }
}

/**
 * ncm_serialize_array_from_variant:
 * @ser: a #NcmSerialize.
 * @var: a #GVariant containing an array of objects.
 *
 * Creates a new #NcmObjArray from a #GVariant.
 *
 * Returns: (transfer full): a new #NcmObjArray.
 */
NcmObjArray *
ncm_serialize_array_from_variant (NcmSerialize *ser, GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_ARRAY_TYPE)));
  {
    guint i, n = g_variant_n_children (var);
    NcmObjArray *oa = ncm_obj_array_sized_new (n);

    for (i = 0; i < n; i++)
    {
      GVariant *cvar = g_variant_get_child_value (var, i);
      GObject *cobj  = ncm_serialize_from_variant (ser, cvar);

      ncm_obj_array_add (oa, cobj);
      g_object_unref (cobj);
      g_variant_unref (cvar);
    }

    return oa;
  }
}

/**
 * ncm_serialize_dict_str_from_variant:
 * @ser: a #NcmSerialize
 * @var: a #GVariant containing a dictionary of string keys
 *
 * Creates a new #NcmObjDictStr from a #GVariant.
 *
 * Returns: (transfer full): a new #NcmObjDictStr.
 */
NcmObjDictStr *
ncm_serialize_dict_str_from_variant (NcmSerialize *ser, GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_STR_TYPE)));
  {
    guint i, n = g_variant_n_children (var);
    NcmObjDictStr *ods = ncm_obj_dict_str_new ();

    for (i = 0; i < n; i++)
    {
      GVariant *cvar = g_variant_get_child_value (var, i);
      GVariant *key  = g_variant_get_child_value (cvar, 0);
      GVariant *val  = g_variant_get_child_value (cvar, 1);
      GObject *cobj  = ncm_serialize_from_variant (ser, val);

      ncm_obj_dict_str_set (ods, g_variant_get_string (key, NULL), cobj);

      g_object_unref (cobj);

      g_variant_unref (key);
      g_variant_unref (val);
      g_variant_unref (cvar);
    }

    return ods;
  }
}

/**
 * ncm_serialize_dict_int_from_variant:
 * @ser: a #NcmSerialize
 * @var: a #GVariant containing a dictionary of integers keys
 *
 * Creates a new #NcmObjDictInt from a #GVariant.
 *
 * Returns: (transfer full): a new #NcmObjDictInt.
 */
NcmObjDictInt *
ncm_serialize_dict_int_from_variant (NcmSerialize *ser, GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_INT_TYPE)));
  {
    guint i, n = g_variant_n_children (var);
    NcmObjDictInt *odi = ncm_obj_dict_int_new ();

    for (i = 0; i < n; i++)
    {
      GVariant *cvar = g_variant_get_child_value (var, i);
      GVariant *key  = g_variant_get_child_value (cvar, 0);
      GVariant *val  = g_variant_get_child_value (cvar, 1);
      GObject *cobj  = ncm_serialize_from_variant (ser, val);

      ncm_obj_dict_int_set (odi, g_variant_get_int32 (key), cobj);

      g_object_unref (cobj);

      g_variant_unref (key);
      g_variant_unref (val);
      g_variant_unref (cvar);
    }

    return odi;
  }
}

/**
 * ncm_serialize_var_dict_from_variant:
 * @ser: a #NcmSerialize
 * @var: a #GVariant containing a dictionary of string and variants
 *
 * Creates a new #NcmVarDict from a #GVariant.
 *
 * Returns: (transfer full): a new #NcmObjDictStr.
 */
NcmVarDict *
ncm_serialize_var_dict_from_variant (NcmSerialize *ser, GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE)));
  {
    guint i, n = g_variant_n_children (var);
    NcmVarDict *vd = ncm_var_dict_new ();

    for (i = 0; i < n; i++)
    {
      GVariant *cvar    = g_variant_get_child_value (var, i);
      GVariant *key     = g_variant_get_child_value (cvar, 0);
      GVariant *val_var = g_variant_get_child_value (cvar, 1);
      GVariant *val     = g_variant_get_variant (val_var);

      ncm_var_dict_set_variant (vd, g_variant_get_string (key, NULL), val);

      g_variant_unref (key);
      g_variant_unref (val_var);
      g_variant_unref (val);
      g_variant_unref (cvar);
    }

    return vd;
  }
}

/**
 * ncm_serialize_from_string:
 * @ser: a #NcmSerialize.
 * @obj_ser: String containing the serialized version of the object.
 *
 * Parses the serialized and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_string (NcmSerialize *ser, const gchar *obj_ser)
{
  GMatchInfo *match_info = NULL;
  GObject *obj           = NULL;
  GError *error          = NULL;
  gchar *error_msg       = NULL;
  GVariant *var_obj      = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), obj_ser, NULL, NULL, &error);

  if (error != NULL)
    error_msg = g_strdup (error->message);

  g_clear_error (&error);

  if (var_obj != NULL)
  {
    obj = ncm_serialize_from_variant (ser, var_obj);
    g_variant_unref (var_obj);
  }
  else if (g_regex_match (ser->parse_obj_regex, obj_ser, 0, &match_info))
  {
    GVariant *params = NULL;
    gchar *obj_name  = g_match_info_fetch (match_info, 1);
    gchar *obj_prop  = g_match_info_fetch (match_info, 2);

    g_free (error_msg);

    if ((obj_prop != NULL) && (strlen (obj_prop) > 0))
    {
      params = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE),
                                obj_prop, NULL, NULL, &error);

      if (error != NULL)
      {
        g_error ("ncm_serialize_from_string: cannot parse object `%s' parameters `%s', error: %s.", obj_name, obj_prop, error->message);
        g_error_free (error);
      }
    }

    obj = ncm_serialize_from_name_params (ser, obj_name, params);

    if (params != NULL)
      g_variant_unref (params);

    g_free (obj_name);
    g_free (obj_prop);
    g_match_info_free (match_info);
  }
  else
  {
    g_match_info_free (match_info);
    g_error ("ncm_serialize_from_string: cannot identify object (%s) in string '%s', variant error `%s'.",
             NCM_SERIALIZE_OBJECT_TYPE,
             obj_ser,
             error_msg);
    g_free (error_msg);
  }

  return obj;
}

static const GVariantType *_ncm_serialize_gtype_to_gvariant_type (GType t);

#ifdef HAVE_LIBFYAML

GObject *
_ncm_serialize_from_node (NcmSerialize *ser, struct fy_node *root)
{
  GObject *obj = NULL;

  if (fy_node_is_mapping (root))
  {
    struct fy_node_pair *obj_type_params = fy_node_mapping_get_by_index (root, 0);
    struct fy_node *obj_type_str         = fy_node_pair_key (obj_type_params);
    struct fy_node *obj_params           = fy_node_pair_value (obj_type_params);
    struct fy_anchor *obj_type_anchor    = fy_node_get_anchor (obj_type_str);
    gchar *anchor_label                  = NULL;

    g_assert (fy_node_is_scalar (obj_type_str));

    if (obj_type_anchor != NULL)
    {
      const gchar *anchor_text = NULL;
      size_t len;

      anchor_text  = fy_anchor_get_text (obj_type_anchor, &len);
      anchor_label = g_strdup_printf ("%.*s", (gint) len, anchor_text);
    }

    if (fy_node_is_alias (obj_params))
    {
      const gchar *alias_label = fy_node_get_scalar0 (obj_params);

      if (!ncm_serialize_contain_name (ser, alias_label))
        g_error ("_ncm_serialize_from_node: object alias `%s' is not saved.", alias_label);

      obj = ncm_serialize_get_by_name (ser, alias_label);
    }
    else if (fy_node_is_mapping (obj_params))
    {
      const gchar *obj_type_str0 = fy_node_get_scalar0 (obj_type_str);
      GType obj_type             = g_type_from_name (obj_type_str0);
      guint n_properties, i;

      if (obj_type == G_TYPE_INVALID)
      {
        g_error ("_ncm_serialize_from_node: object type `%s' is not registered.", obj_type_str0);
      }
      else
      {
        GObjectClass *klass     = g_type_class_ref (obj_type);
        GParamSpec **prop       = g_object_class_list_properties (klass, &n_properties);
        guint n_properties_yaml = fy_node_mapping_item_count (obj_params);

        if (n_properties_yaml > n_properties)
        {
          g_error ("_ncm_serialize_from_node: object yaml `%s' has more properties than expected.", obj_type_str0);
        }
        else if (n_properties_yaml == 0)
        {
          obj = g_object_new (obj_type, NULL);
        }
        else
        {
          GValue *values      = g_new (GValue, n_properties_yaml);
          const gchar **names = g_new (const gchar *, n_properties_yaml);
          guint n_found       = 0;

          for (i = 0; i < n_properties_yaml; i++)
          {
            struct fy_node_pair *prop_yaml = fy_node_mapping_get_by_index (obj_params, i);
            struct fy_node *prop_name      = fy_node_pair_key (prop_yaml);
            struct fy_node *prop_val       = fy_node_pair_value (prop_yaml);
            gboolean found                 = FALSE;
            guint j;

            if (!fy_node_is_scalar (prop_name))
              g_error ("_ncm_serialize_from_node: object property name must be a scalar.");

            names[i] = fy_node_get_scalar0 (prop_name);

            for (j = 0; j < n_properties; j++)
            {
              GParamSpec *pspec = prop[j];

              if ((pspec->flags & G_PARAM_READWRITE) != G_PARAM_READWRITE)
                continue;

              if (strcmp (pspec->name, names[i]) == 0)
              {
                GValue lval = G_VALUE_INIT;

                g_value_init (&lval, pspec->value_type);

                if (g_type_is_a (pspec->value_type, NCM_TYPE_OBJ_ARRAY))
                {
                  if (!fy_node_is_sequence (prop_val))
                  {
                    g_error ("_ncm_serialize_from_node: object property of type NcmObjArray `%s' must be a sequence.", names[i]);
                  }
                  else
                  {
                    NcmObjArray *obj_array = ncm_obj_array_new ();
                    guint n_items          = fy_node_sequence_item_count (prop_val);
                    guint k;

                    for (k = 0; k < n_items; k++)
                    {
                      struct fy_node *item = fy_node_sequence_get_by_index (prop_val, k);
                      GObject *item_obj    = _ncm_serialize_from_node (ser, item);

                      ncm_obj_array_add (obj_array, item_obj);

                      g_object_unref (item_obj);
                    }

                    g_value_take_boxed (&lval, obj_array);
                  }
                }
                else if (g_type_is_a (pspec->value_type, NCM_TYPE_OBJ_DICT_STR))
                {
                  if (!fy_node_is_mapping (prop_val))
                  {
                    g_error ("_ncm_serialize_from_node: object property of type NcmObjDictStr `%s' must be a mapping.", names[i]);
                  }
                  else
                  {
                    NcmObjDictStr *obj_dict_str = ncm_obj_dict_str_new ();
                    guint n_items               = fy_node_mapping_item_count (prop_val);
                    guint k;

                    for (k = 0; k < n_items; k++)
                    {
                      struct fy_node_pair *item = fy_node_mapping_get_by_index (prop_val, k);
                      struct fy_node *item_key  = fy_node_pair_key (item);
                      struct fy_node *item_val  = fy_node_pair_value (item);
                      const gchar *key          = fy_node_get_scalar0 (item_key);
                      GObject *val_obj          = _ncm_serialize_from_node (ser, item_val);

                      ncm_obj_dict_str_set (obj_dict_str, key, val_obj);

                      g_object_unref (val_obj);
                    }

                    g_value_take_boxed (&lval, obj_dict_str);
                  }
                }
                else if (g_type_is_a (pspec->value_type, NCM_TYPE_OBJ_DICT_INT))
                {
                  if (!fy_node_is_mapping (prop_val))
                  {
                    g_error ("_ncm_serialize_from_node: object property of type NcmObjDictInt `%s' must be a mapping.", names[i]);
                  }
                  else
                  {
                    NcmObjDictInt *obj_dict_int = ncm_obj_dict_int_new ();
                    guint n_items               = fy_node_mapping_item_count (prop_val);
                    guint k;

                    for (k = 0; k < n_items; k++)
                    {
                      struct fy_node_pair *item = fy_node_mapping_get_by_index (prop_val, k);
                      struct fy_node *item_key  = fy_node_pair_key (item);
                      struct fy_node *item_val  = fy_node_pair_value (item);
                      GObject *val_obj          = _ncm_serialize_from_node (ser, item_val);
                      glong lkey;

                      errno = 0;
                      lkey  = strtol (fy_node_get_scalar0 (item_key), NULL, 10);

                      if (errno != 0)
                        g_error ("_ncm_serialize_from_node: object dict int yaml key must be an integer, received '%s'.",
                                 fy_node_get_scalar0 (item_key));

                      ncm_obj_dict_int_set (obj_dict_int, lkey, val_obj);

                      g_object_unref (val_obj);
                    }

                    g_value_take_boxed (&lval, obj_dict_int);
                  }
                }
                else if (g_type_is_a (pspec->value_type, NCM_TYPE_VAR_DICT))
                {
                  if (!fy_node_is_mapping (prop_val))
                  {
                    g_error ("_ncm_serialize_from_node: object property of type NcmVarDict `%s' must be a mapping.", names[i]);
                  }
                  else
                  {
                    NcmVarDict *var_dict = ncm_var_dict_new ();
                    guint n_items        = fy_node_mapping_item_count (prop_val);
                    guint k;

                    for (k = 0; k < n_items; k++)
                    {
                      struct fy_node_pair *item = fy_node_mapping_get_by_index (prop_val, k);
                      struct fy_node *item_key  = fy_node_pair_key (item);
                      struct fy_node *item_val  = fy_node_pair_value (item);
                      const gchar *key          = fy_node_get_scalar0 (item_key);
                      gchar *value_string       = fy_emit_node_to_string (item_val, FYECF_WIDTH_INF | FYECF_MODE_FLOW_ONELINE);
                      GError *error             = NULL;
                      GVariant *item_val_var    = g_variant_parse (NULL, value_string, NULL, NULL, &error);

                      if (error != NULL)
                        g_error ("_ncm_serialize_var_dict_from_yaml: cannot parse dictionary element `%s' value `%s': %s.",
                                 key, value_string, error->message);

                      ncm_var_dict_set_variant (var_dict, key, item_val_var);

                      g_variant_unref (item_val_var);
                      g_free (value_string);
                    }

                    g_value_take_boxed (&lval, var_dict);
                  }
                }
                else if (g_type_is_a (pspec->value_type, G_TYPE_OBJECT))
                {
                  if (fy_node_is_mapping (prop_val))
                  {
                    GObject *item_obj = _ncm_serialize_from_node (ser, prop_val);

                    g_value_take_object (&lval, item_obj);
                  }
                  else if (g_type_is_a (pspec->value_type, NCM_TYPE_VECTOR))
                  {
                    if (!fy_node_is_sequence (prop_val))
                    {
                      g_error ("_ncm_serialize_from_node: object property of type NcmVector `%s' "
                               "must be either a sequence of a proper NcmVector object", names[i]);
                    }
                    else
                    {
                      GError *error       = NULL;
                      gchar *prop_val_str = fy_emit_node_to_string (prop_val, FYECF_WIDTH_INF | FYECF_MODE_FLOW_ONELINE);
                      GVariant *var       = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_VECTOR_TYPE), prop_val_str, NULL, NULL, &error);
                      NcmVector *vector   = NULL;

                      if (error != NULL)
                        g_error ("_ncm_serialize_from_node: cannot parse property `%s' value `%s': %s.", names[i], prop_val_str, error->message);

                      vector = ncm_vector_new_variant (var);

                      g_value_take_object (&lval, vector);

                      g_free (prop_val_str);
                    }
                  }
                  else if (g_type_is_a (pspec->value_type, NCM_TYPE_MATRIX))
                  {
                    if (!fy_node_is_sequence (prop_val))
                    {
                      g_error ("_ncm_serialize_from_node: object property of type NcmMatrix `%s' "
                               "must be either a sequence of a proper NcmMatrix object", names[i]);
                    }
                    else
                    {
                      GError *error       = NULL;
                      gchar *prop_val_str = fy_emit_node_to_string (prop_val, FYECF_WIDTH_INF | FYECF_MODE_FLOW_ONELINE);
                      GVariant *var       = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_MATRIX_TYPE), prop_val_str, NULL, NULL, &error);
                      NcmMatrix *matrix   = NULL;

                      if (error != NULL)
                        g_error ("_ncm_serialize_from_node: cannot parse property `%s' value `%s': %s.", names[i], prop_val_str, error->message);

                      matrix = ncm_matrix_new_variant (var);

                      g_value_take_object (&lval, matrix);

                      g_free (prop_val_str);
                    }
                  }
                  else
                  {
                    g_error ("_ncm_serialize_from_node: object property of type `%s' must be a mapping.", names[i]);
                  }
                }
                else if (g_type_is_a (pspec->value_type, G_TYPE_VARIANT))
                {
                  GParamSpecVariant *pspec_var = G_PARAM_SPEC_VARIANT (pspec);

                  g_assert (G_IS_PARAM_SPEC_VARIANT (pspec));
                  {
                    GError *error       = NULL;
                    gchar *prop_val_str = fy_emit_node_to_string (prop_val, FYECF_WIDTH_INF | FYECF_MODE_FLOW_ONELINE);
                    GVariant *var       = g_variant_parse (pspec_var->type, prop_val_str, NULL, NULL, &error);

                    if (error != NULL)
                      g_error ("_ncm_serialize_from_node: cannot parse property `%s' value `%s': %s.", names[i], prop_val_str, error->message);

                    g_value_take_variant (&lval, var);
                    g_free (prop_val_str);
                  }
                }
                else if (g_type_is_a (pspec->value_type, NCM_TYPE_DTUPLE2))
                {
                  GError *error             = NULL;
                  const gchar *prop_val_str = fy_node_get_scalar0 (prop_val);
                  GVariant *var             = g_variant_parse (G_VARIANT_TYPE (NCM_DTUPLE2_TYPE), prop_val_str, NULL, NULL, &error);
                  NcmDTuple2 *dtuple2;

                  if (error != NULL)
                    g_error ("_ncm_serialize_from_node: cannot parse property `%s' value `%s': %s.", names[i], prop_val_str, error->message);

                  dtuple2 = ncm_dtuple2_new_from_variant (var);
                  g_value_take_boxed (&lval, dtuple2);
                }
                else if (g_type_is_a (pspec->value_type, NCM_TYPE_DTUPLE3))
                {
                  GError *error             = NULL;
                  const gchar *prop_val_str = fy_node_get_scalar0 (prop_val);
                  GVariant *var             = g_variant_parse (G_VARIANT_TYPE (NCM_DTUPLE3_TYPE), prop_val_str, NULL, NULL, &error);
                  NcmDTuple3 *dtuple3;

                  if (error != NULL)
                    g_error ("_ncm_serialize_from_node: cannot parse property `%s' value `%s': %s.", names[i], prop_val_str, error->message);

                  dtuple3 = ncm_dtuple3_new_from_variant (var);
                  g_value_take_boxed (&lval, dtuple3);
                }
                else if (g_type_is_a (pspec->value_type, G_TYPE_STRV))
                {
                  if (!fy_node_is_sequence (prop_val))
                  {
                    g_error ("_ncm_serialize_from_node: object property of type G_TYPE_STRV `%s' must be a sequence.", names[i]);
                  }
                  else
                  {
                    guint n_items = fy_node_sequence_item_count (prop_val);
                    gchar **strv  = g_new (gchar *, n_items + 1);
                    guint k;

                    for (k = 0; k < n_items; k++)
                    {
                      struct fy_node *item = fy_node_sequence_get_by_index (prop_val, k);

                      strv[k] = g_strdup (fy_node_get_scalar0 (item));
                    }

                    strv[n_items] = NULL;

                    g_value_take_boxed (&lval, strv);
                  }
                }
                else
                {
                  GError *error       = NULL;
                  gchar *prop_val_str = fy_emit_node_to_string (prop_val, FYECF_WIDTH_INF | FYECF_MODE_FLOW_ONELINE);
                  GVariant *var       = g_variant_parse (_ncm_serialize_gtype_to_gvariant_type (pspec->value_type), prop_val_str, NULL, NULL, &error);

                  if (error != NULL)
                    g_error ("_ncm_serialize_from_node: cannot parse property `%s' value `%s': %s.", names[i], prop_val_str, error->message);

                  g_dbus_gvariant_to_gvalue (var, &lval);
                  g_variant_unref (var);

                  g_free (prop_val_str);
                }

                values[i] = lval;

                found = TRUE;
                n_found++;
                break;
              }
            }

            if (!found)
              g_error ("_ncm_serialize_from_node: object do not have property `%s'.", names[i]);
          }

          if (n_found != n_properties_yaml)
            g_error ("_ncm_serialize_from_node: object `%s' has less properties than expected.", obj_type_str0);

          obj = g_object_new_with_properties (obj_type, n_properties_yaml, names, values);

          for (i = 0; i < n_properties_yaml; i++)
            g_value_unset (&values[i]);

          g_free (names);
          g_free (values);
        }

        g_free (prop);
        g_type_class_unref (klass);
      }
    }
    else
    {
      g_error ("_ncm_serialize_from_node: object parameters must be a mapping or an alias.");
    }

    if ((anchor_label != NULL) && (ser->opts & NCM_SERIALIZE_OPT_AUTOSAVE_SER))
    {
      ncm_serialize_set (ser, obj, anchor_label, FALSE);
      g_free (anchor_label);
    }
  }
  else
  {
    g_error ("_ncm_serialize_from_node: object must be a mapping.");
  }

  return obj;
}

#endif /* HAVE_LIBFYAML */

/**
 * ncm_serialize_from_yaml:
 * @ser: a #NcmSerialize
 * @yaml_obj: string containing the serialized version of the object in YAML format
 *
 * Parses the serialized string in @yaml_obj and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_yaml (NcmSerialize *ser, const gchar *yaml_obj)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_build_from_string (NULL, yaml_obj, strlen (yaml_obj));

  if (doc == NULL)
  {
    g_error ("ncm_serialize_from_yaml: cannot parse YAML object ###\n%s\n###.", yaml_obj);

    return NULL;
  }
  else
  {
    struct fy_node *root = fy_document_root (doc);

    if (!fy_node_is_mapping (root))
    {
      fy_document_destroy (doc);

      g_error ("ncm_serialize_from_yaml: object yaml root must be a mapping.");

      return NULL;
    }

    if (fy_node_mapping_item_count (root) > 1)
    {
      fy_document_destroy (doc);

      g_error ("ncm_serialize_from_yaml: object yaml root has more than one item.");

      return NULL;
    }

    {
      GObject *obj = _ncm_serialize_from_node (ser, root);

      fy_document_destroy (doc);

      return obj;
    }
  }

#else

  g_error ("ncm_serialize_from_yaml: libfyaml not available.");

  return NULL;

#endif /* HAVE_LIBFYAML */
}

#ifdef HAVE_LIBFYAML

static NcmObjArray *
_ncm_serialize_array_from_yaml_node (NcmSerialize *ser, struct fy_node *root)
{
  NcmObjArray *array = ncm_obj_array_new ();

  if (fy_node_is_sequence (root))
  {
    guint n_items = fy_node_sequence_item_count (root);
    guint i;

    for (i = 0; i < n_items; i++)
    {
      struct fy_node *item = fy_node_sequence_get_by_index (root, i);
      GObject *item_obj    = _ncm_serialize_from_node (ser, item);

      ncm_obj_array_add (array, item_obj);

      g_object_unref (item_obj);
    }
  }
  else
  {
    g_error ("_ncm_serialize_array_from_yaml_node: object array yaml root must be a sequence.");
  }

  return array;
}

#endif /* HAVE_LIBFYAML */

/**
 * ncm_serialize_array_from_yaml:
 * @ser: a #NcmSerialize
 * @yaml_obj: string containing the serialized version of the object in YAML format
 *
 * Parses the serialized string in @yaml_obj and returns an array of newly created objects.
 *
 * Returns: (transfer full): A new #NcmObjArray.
 */
NcmObjArray *
ncm_serialize_array_from_yaml (NcmSerialize *ser, const gchar *yaml_obj)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_build_from_string (NULL, yaml_obj, strlen (yaml_obj));

  if (doc == NULL)
  {
    g_error ("ncm_serialize_from_yaml: cannot parse YAML object ###\n%s\n###.", yaml_obj);

    return NULL;
  }
  else
  {
    struct fy_node *root = fy_document_root (doc);
    NcmObjArray *array   = _ncm_serialize_array_from_yaml_node (ser, root);

    fy_document_destroy (doc);

    return array;
  }

#else

  g_error ("ncm_serialize_array_from_yaml: libfyaml not available.");

  return NULL;

#endif /* HAVE_LIBFYAML */
}

/**
 * ncm_serialize_dict_str_from_yaml:
 * @ser: a #NcmSerialize
 * @yaml_obj: string containing the serialized version of the object in YAML format
 *
 * Parses the serialized string in @yaml_obj and returns a #NcmObjDictStr containing
 * the object names as keys and the serialized objects as values.
 *
 * Returns: (transfer full): A new #NcmObjDictStr.
 */
NcmObjDictStr *
ncm_serialize_dict_str_from_yaml (NcmSerialize *ser, const gchar *yaml_obj)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_build_from_string (NULL, yaml_obj, strlen (yaml_obj));

  if (doc == NULL)
  {
    g_error ("ncm_serialize_from_yaml: cannot parse YAML object ###\n%s\n###.", yaml_obj);

    return NULL;
  }
  else
  {
    struct fy_node *root = fy_document_root (doc);
    NcmObjDictStr *dict  = ncm_obj_dict_str_new ();

    if (fy_node_is_mapping (root))
    {
      guint n_items = fy_node_mapping_item_count (root);
      guint i;

      for (i = 0; i < n_items; i++)
      {
        struct fy_node_pair *item = fy_node_mapping_get_by_index (root, i);
        struct fy_node *item_key  = fy_node_pair_key (item);
        struct fy_node *item_val  = fy_node_pair_value (item);
        const gchar *key          = fy_node_get_scalar0 (item_key);
        GObject *val              = _ncm_serialize_from_node (ser, item_val);

        ncm_obj_dict_str_set (dict, key, val);

        g_object_unref (val);
      }
    }
    else
    {
      g_error ("_ncm_serialize_dict_str_from_yaml: object dict str yaml root must be a mapping.");
    }

    fy_document_destroy (doc);

    return dict;
  }

#else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_dict_str_from_yaml: libfyaml not available.");

  return NULL;

#endif /* HAVE_LIBFYAML */
}

/**
 * ncm_serialize_dict_int_from_yaml:
 * @ser: a #NcmSerialize
 * @yaml_obj: string containing the serialized version of the object in YAML format
 *
 * Parses the serialized string in @yaml_obj and returns a #NcmObjDictInt containing
 * the object names as keys and the serialized objects as values.
 *
 * Returns: (transfer full): A new #NcmObjDictInt.
 */
NcmObjDictInt *
ncm_serialize_dict_int_from_yaml (NcmSerialize *ser, const gchar *yaml_obj)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_build_from_string (NULL, yaml_obj, strlen (yaml_obj));

  if (doc == NULL)
  {
    g_error ("ncm_serialize_from_yaml: cannot parse YAML object ###\n%s\n###.", yaml_obj);

    return NULL;
  }
  else
  {
    struct fy_node *root = fy_document_root (doc);
    NcmObjDictInt *dict  = ncm_obj_dict_int_new ();

    if (fy_node_is_mapping (root))
    {
      guint n_items = fy_node_mapping_item_count (root);
      guint i;

      for (i = 0; i < n_items; i++)
      {
        struct fy_node_pair *item = fy_node_mapping_get_by_index (root, i);
        struct fy_node *item_key  = fy_node_pair_key (item);
        struct fy_node *item_val  = fy_node_pair_value (item);
        GObject *val              = _ncm_serialize_from_node (ser, item_val);
        glong lkey;

        errno = 0;
        lkey  = strtol (fy_node_get_scalar0 (item_key), NULL, 10);

        if (errno != 0)
          g_error ("_ncm_serialize_dict_int_from_yaml: object dict int yaml key must be an integer, received '%s'.",
                   fy_node_get_scalar0 (item_key));

        ncm_obj_dict_int_set (dict, lkey, val);

        g_object_unref (val);
      }
    }
    else
    {
      g_error ("_ncm_serialize_dict_int_from_yaml: object dict int yaml root must be a mapping.");
    }

    fy_document_destroy (doc);

    return dict;
  }

#else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_dict_int_from_yaml: libfyaml not available.");

  return NULL;

#endif /* HAVE_LIBFYAML */
}

/**
 * ncm_serialize_var_dict_from_yaml:
 * @ser: a #NcmSerialize
 * @yaml_obj: string containing the serialized version of the #NcmVarDict in YAML format
 *
 * Parses the serialized string in @yaml_obj and returns a #NcmVarDict containing
 * the object names as keys and the serialized objects as values.
 *
 * Returns: (transfer full): A new #NcmVarDict.
 */
NcmVarDict *
ncm_serialize_var_dict_from_yaml (NcmSerialize *ser, const gchar *yaml_obj)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_build_from_string (NULL, yaml_obj, strlen (yaml_obj));

  if (doc == NULL)
  {
    g_error ("ncm_serialize_from_yaml: cannot parse YAML object ###\n%s\n###.", yaml_obj);

    return NULL;
  }
  else
  {
    struct fy_node *root = fy_document_root (doc);
    NcmVarDict *dict     = ncm_var_dict_new ();

    if (fy_node_is_mapping (root))
    {
      guint n_items = fy_node_mapping_item_count (root);
      guint i;

      for (i = 0; i < n_items; i++)
      {
        struct fy_node_pair *item = fy_node_mapping_get_by_index (root, i);
        struct fy_node *item_key  = fy_node_pair_key (item);
        struct fy_node *item_val  = fy_node_pair_value (item);
        const gchar *key          = fy_node_get_scalar0 (item_key);
        gchar *value_string       = fy_emit_node_to_string (item_val, FYECF_WIDTH_INF | FYECF_MODE_FLOW_ONELINE);
        GError *error             = NULL;
        GVariant *item_val_var    = g_variant_parse (NULL, value_string, NULL, NULL, &error);

        if (error != NULL)
          g_error ("ncm_serialize_var_dict_from_yaml: cannot parse dictionary element `%s' value `%s': %s.", key, value_string, error->message);

        ncm_var_dict_set_variant (dict, key, item_val_var);

        g_variant_unref (item_val_var);
        g_free (value_string);
      }
    }
    else
    {
      g_error ("_ncm_serialize_var_dict_from_yaml: object var dict yaml root must be a mapping.");
    }

    fy_document_destroy (doc);

    return dict;
  }

  #else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_var_dict_from_yaml: libfyaml not available.");

  return NULL;

#endif /* HAVE_LIBFYAML */
}

/**
 * ncm_serialize_from_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the object
 *
 * Parses the serialized string in @filename and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_file (NcmSerialize *ser, const gchar *filename)
{
  GError *error = NULL;
  gchar *file   = NULL;
  gsize length  = 0;
  GObject *obj  = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  obj = ncm_serialize_from_string (ser, file);

  g_free (file);

  return obj;
}

/**
 * ncm_serialize_from_binfile:
 * @ser: a #NcmSerialize.
 * @filename: File containing the binary serialized version of the object.
 *
 * Parses the serialized binary data in @filename and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_binfile (NcmSerialize *ser, const gchar *filename)
{
  GError *error = NULL;
  gchar *file   = NULL;
  gsize length  = 0;
  GObject *obj  = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("_nc_data_snia_cov_load_matrix: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  {
    GVariant *obj_ser = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE),
                                                 file,
                                                 length,
                                                 TRUE,
                                                 g_free,
                                                 file
                                                );

    obj = ncm_serialize_from_variant (ser, obj_ser);

    g_variant_unref (obj_ser);
  }

  return obj;
}

/**
 * ncm_serialize_var_dict_from_variant_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the #NcmVarDict
 * @binary: Whether the file contains binary data or not
 *
 * Parses the serialized string in @filename and returns a #NcmVarDict containing
 * the object names as keys and the serialized objects as values.
 *
 * Returns: (transfer full): A new #NcmVarDict.
 */
NcmVarDict *
ncm_serialize_var_dict_from_variant_file (NcmSerialize *ser, const gchar *filename, gboolean binary)
{
  GError *error  = NULL;
  gchar *file    = NULL;
  gsize length   = 0;
  NcmVarDict *vd = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  if (binary)
  {
    GVariant *vd_ser = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE),
                                                file,
                                                length,
                                                TRUE,
                                                g_free,
                                                file
                                               );

    vd = ncm_serialize_var_dict_from_variant (ser, vd_ser);

    g_variant_unref (vd_ser);
  }
  else
  {
    GVariant *vd_ser = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE),
                                        file,
                                        NULL,
                                        NULL,
                                        &error);

    if (error != NULL)
      g_error ("ncm_serialize_var_dict_from_variant_file: cannot parse file %s: %s",
               filename, error->message);

    vd = ncm_serialize_var_dict_from_variant (ser, vd_ser);

    g_variant_unref (vd_ser);
    g_free (file);
  }

  return vd;
}

/**
 * ncm_serialize_from_yaml_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the object in YAML format
 *
 * Parses the YAML in @filename and returns the newly created object.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_yaml_file (NcmSerialize *ser, const gchar *filename)
{
  GError *error = NULL;
  gchar *file   = NULL;
  gsize length  = 0;
  GObject *obj  = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  obj = ncm_serialize_from_yaml (ser, file);

  g_free (file);

  return obj;
}

/**
 * ncm_serialize_array_from_key_file:
 * @ser: a #NcmSerialize
 * @filename: oa filename
 *
 * Loads a #NcmObjArray from a file using a #NcmSerialize and a #GKeyFile.
 *
 * Returns: (transfer full): a new #NcmObjArray.
 */
NcmObjArray *
ncm_serialize_array_from_key_file (NcmSerialize *ser, const gchar *filename)
{
  NcmObjArray *oa  = ncm_obj_array_new ();
  GKeyFile *oafile = g_key_file_new ();
  GError *error    = NULL;
  gchar **groups   = NULL;
  gsize ngroups    = 0;
  guint i;

  if (!g_key_file_load_from_file (oafile, filename, G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS, &error))
  {
    g_error ("ncm_obj_array_load: Invalid GObject array file: %s %s", filename, error->message);

    return NULL;
  }

  if (g_key_file_has_group (oafile, "NcmObjArray"))
  {
    g_key_file_remove_group (oafile, "NcmObjArray", &error);

    if (error != NULL)
      g_error ("ncm_obj_array_load: %s", error->message);
  }

  groups = g_key_file_get_groups (oafile, &ngroups);

  for (i = 0; i < ngroups; i++)
  {
    GString *obj_ser = g_string_sized_new (200);
    gchar **a_pos    = g_strsplit (groups[i], ":", 2);
    gchar *obj_name  = NULL;

    g_assert_cmpuint (g_strv_length (a_pos), ==, 2);

    if (!g_key_file_has_key (oafile, groups[i], NCM_SERIALIZE_OBJECT_ARRAY_OBJ_NAME_STR, &error))
    {
      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);

      g_error ("ncm_obj_array_load: Every group must contain the key `%s' containing the object name.", NCM_SERIALIZE_OBJECT_ARRAY_OBJ_NAME_STR);
    }
    else
    {
      obj_name = g_key_file_get_value (oafile, groups[i], NCM_SERIALIZE_OBJECT_ARRAY_OBJ_NAME_STR, &error);

      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);

      g_key_file_remove_key (oafile, groups[i], NCM_SERIALIZE_OBJECT_ARRAY_OBJ_NAME_STR, &error);

      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);
    }

    g_string_append_printf (obj_ser, "(\'%s\', @a{sv} {", obj_name);
    g_clear_pointer (&obj_name, g_free);

    {
      gsize nkeys  = 0;
      gchar **keys = g_key_file_get_keys (oafile, groups[i], &nkeys, &error);
      guint j;

      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);

      for (j = 0; j < nkeys; j++)
      {
        gchar *propval = g_key_file_get_value (oafile, groups[i], keys[j], &error);

        if (error != NULL)
          g_error ("ncm_obj_array_load: %s", error->message);

        g_string_append_printf (obj_ser, "\'%s\':<%s>", keys[j], propval);
        g_free (propval);

        if (j + 1 != nkeys)
          g_string_append (obj_ser, ", ");
      }

      g_string_append (obj_ser, "})");
      g_strfreev (keys);
    }

    {
      GObject *obj = ncm_serialize_from_string (ser, obj_ser->str);

      g_assert (G_IS_OBJECT (obj));

      ncm_obj_array_add (oa, obj);
      g_object_unref (obj);
    }
    g_string_free (obj_ser, TRUE);
    g_strfreev (a_pos);
  }

  g_key_file_unref (oafile);
  g_strfreev (groups);

  return oa;
}

/**
 * ncm_serialize_array_from_yaml_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the object in YAML format
 *
 * Parses the YAML in @filename and returns an array of newly created objects.
 *
 * Returns: (transfer full): A new #NcmObjArray.
 */
NcmObjArray *
ncm_serialize_array_from_yaml_file (NcmSerialize *ser, const gchar *filename)
{
  GError *error      = NULL;
  gchar *file        = NULL;
  gsize length       = 0;
  NcmObjArray *array = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  array = ncm_serialize_array_from_yaml (ser, file);

  g_free (file);

  return array;
}

/**
 * ncm_serialize_dict_str_from_yaml_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the object in YAML format
 *
 * Parses the YAML in @filename and returns a #NcmObjDictStr containing
 *
 * Returns: (transfer full): A new #NcmObjDictStr.
 */
NcmObjDictStr *
ncm_serialize_dict_str_from_yaml_file (NcmSerialize *ser, const gchar *filename)
{
  GError *error       = NULL;
  gchar *file         = NULL;
  gsize length        = 0;
  NcmObjDictStr *dict = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  dict = ncm_serialize_dict_str_from_yaml (ser, file);

  g_free (file);

  return dict;
}

/**
 * ncm_serialize_dict_int_from_yaml_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the object in YAML format
 *
 * Parses the YAML in @filename and returns a #NcmObjDictInt containing
 *
 * Returns: (transfer full): A new #NcmObjDictInt.
 */
NcmObjDictInt *
ncm_serialize_dict_int_from_yaml_file (NcmSerialize *ser, const gchar *filename)
{
  GError *error       = NULL;
  gchar *file         = NULL;
  gsize length        = 0;
  NcmObjDictInt *dict = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  dict = ncm_serialize_dict_int_from_yaml (ser, file);

  g_free (file);

  return dict;
}

/**
 * ncm_serialize_var_dict_from_yaml_file:
 * @ser: a #NcmSerialize
 * @filename: File containing the serialized version of the #NcmVarDict in YAML format
 *
 * Parses the YAML in @filename and returns a #NcmVarDict containing
 * the element names as keys and their values.
 *
 * Returns: (transfer full): A new #NcmVarDict.
 */
NcmVarDict *
ncm_serialize_var_dict_from_yaml_file (NcmSerialize *ser, const gchar *filename)
{
  GError *error    = NULL;
  gchar *file      = NULL;
  gsize length     = 0;
  NcmVarDict *dict = NULL;

  g_assert (filename != NULL);

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("ncm_serialize_from_file: cannot open file %s: %s",
             filename, error->message);

  g_assert_cmpint (length, >, 0);

  dict = ncm_serialize_var_dict_from_yaml (ser, file);

  g_free (file);

  return dict;
}

/**
 * ncm_serialize_from_name_params:
 * @ser: a #NcmSerialize.
 * @obj_name: string containing the object name.
 * @params: a #GVariant containing the object parameters.
 *
 * Parses the serialized parameters and returns the newly created object using them.
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_from_name_params (NcmSerialize *ser, const gchar *obj_name, GVariant *params)
{
  GMatchInfo *match_info = NULL;
  gsize nprop            = (params != NULL) ? g_variant_n_children (params) : 0;
  GObject *obj           = NULL;
  gchar *name            = NULL;
  GType gtype;

  if (g_regex_match (ser->is_named_regex, obj_name, 0, &match_info))
  {
    gchar *cobj_name = g_match_info_fetch (match_info, 1);

    name = g_match_info_fetch (match_info, 2);

    gtype = g_type_from_name (cobj_name);
    g_free (cobj_name);
  }
  else
  {
    gtype = g_type_from_name (obj_name);
  }

  g_match_info_free (match_info);

  if (name != NULL)
/*
 *   if (g_hash_table_lookup (ser->saved_name_ser, name) != NULL)
 *     g_error ("ncm_serialize_from_name_params: deserializing object named `%s' but it is present in the list of saved serializations.",
 *              name);
 *   else
 */
    if (ncm_serialize_contain_name (ser, name))
      obj = ncm_serialize_get_by_name (ser, name);


  if (obj != NULL)
  {
    if (nprop > 0)
    {
      g_error ("ncm_serialize_from_name_params: cannot create object `%s' with properties, it was found on saved/named instances.",
               obj_name);
    }
    else
    {
      g_clear_pointer (&name, g_free);

      return obj;
    }
  }

  if (gtype == 0)
    g_error ("ncm_serialize_from_name_params: object `%s' is not registered.", obj_name);

  g_assert (params == NULL || g_variant_is_of_type (params, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE)));

  if (params != NULL)
  {
    GVariantIter *p_iter  = g_variant_iter_new (params);
    const gchar **names   = g_new (const gchar *, nprop);
    GValue *values        = g_new (GValue, nprop);
    GVariant *var         = NULL;
    gboolean is_NcmVector = g_type_is_a (gtype, NCM_TYPE_VECTOR);
    gboolean is_NcmMatrix = g_type_is_a (gtype, NCM_TYPE_MATRIX);
    guint i               = 0;

    while ((var = g_variant_iter_next_value (p_iter)))
    {
      GVariant *var_key = g_variant_get_child_value (var, 0);
      GVariant *var_val = g_variant_get_child_value (var, 1);
      GVariant *val     = g_variant_get_variant (var_val);

      names[i] = g_variant_get_string (var_key, NULL);

      if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_DTUPLE2_TYPE)))
      {
        NcmDTuple2 *dtuple2 = ncm_dtuple2_new_from_variant (val);
        GValue lval         = G_VALUE_INIT;

        g_value_init (&lval, NCM_TYPE_DTUPLE2);
        values[i] = lval;
        g_value_take_boxed (&values[i], dtuple2);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_DTUPLE3_TYPE)))
      {
        NcmDTuple3 *dtuple3 = ncm_dtuple3_new_from_variant (val);
        GValue lval         = G_VALUE_INIT;

        g_value_init (&lval, NCM_TYPE_DTUPLE3);
        values[i] = lval;
        g_value_take_boxed (&values[i], dtuple3);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_VECTOR_TYPE)) && !is_NcmVector)
      {
        NcmVector *vec = ncm_vector_new_variant (val);
        GValue lval    = G_VALUE_INIT;

        g_value_init (&lval, G_TYPE_OBJECT);
        values[i] = lval;
        g_value_take_object (&values[i], vec);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_MATRIX_TYPE)) && !is_NcmMatrix)
      {
        NcmMatrix *mat = ncm_matrix_new_variant (val);
        GValue lval    = G_VALUE_INIT;

        g_value_init (&lval, G_TYPE_OBJECT);
        values[i] = lval;
        g_value_take_object (&values[i], mat);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_STRV_TYPE)))
      {
        gchar **strv = g_variant_dup_strv (val, NULL);
        GValue lval  = G_VALUE_INIT;

        g_value_init (&lval, G_TYPE_STRV);
        values[i] = lval;
        g_value_take_boxed (&values[i], strv);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_ARRAY_TYPE)))
      {
        NcmObjArray *oa = ncm_serialize_array_from_variant (ser, val);
        GValue lval     = G_VALUE_INIT;

        g_value_init (&lval, NCM_TYPE_OBJ_ARRAY);
        values[i] = lval;
        g_value_take_boxed (&values[i], oa);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_STR_TYPE)))
      {
        NcmObjDictStr *dict = ncm_serialize_dict_str_from_variant (ser, val);
        GValue lval         = G_VALUE_INIT;

        g_value_init (&lval, NCM_TYPE_OBJ_DICT_STR);
        values[i] = lval;
        g_value_take_boxed (&values[i], dict);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_INT_TYPE)))
      {
        NcmObjDictInt *dict = ncm_serialize_dict_int_from_variant (ser, val);
        GValue lval         = G_VALUE_INIT;

        g_value_init (&lval, NCM_TYPE_OBJ_DICT_INT);
        values[i] = lval;
        g_value_take_boxed (&values[i], dict);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE)))
      {
        NcmVarDict *dict = ncm_serialize_var_dict_from_variant (ser, val);
        GValue lval      = G_VALUE_INIT;

        g_value_init (&lval, NCM_TYPE_VAR_DICT);
        values[i] = lval;
        g_value_take_boxed (&values[i], dict);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)))
      {
        GVariant *nest_obj_key    = g_variant_get_child_value (val, 0);
        GVariant *nest_obj_params = g_variant_get_child_value (val, 1);
        GValue lval               = G_VALUE_INIT;
        GObject *nest_obj         =
          ncm_serialize_from_name_params (ser,
                                          g_variant_get_string (nest_obj_key, NULL),
                                          nest_obj_params);

        g_value_init (&lval, G_TYPE_OBJECT);
        values[i] = lval;

        g_value_take_object (&values[i], nest_obj);
        g_variant_unref (nest_obj_key);
        g_variant_unref (nest_obj_params);
      }
      else
      {
        g_dbus_gvariant_to_gvalue (val, &values[i]);
      }

      i++;
      g_variant_unref (var_key);
      g_variant_unref (var_val);
      g_variant_unref (val);
      g_variant_unref (var);
    }

    obj = g_object_new_with_properties (gtype, nprop, names, values);

    for (i = 0; i < nprop; i++)
      g_value_unset (&values[i]);

    g_free (names);
    g_free (values);
    g_variant_iter_free (p_iter);
  }
  else
  {
    obj = g_object_new (gtype, NULL);
  }

  if ((name != NULL) && (ser->opts & NCM_SERIALIZE_OPT_AUTOSAVE_SER))
  {
    ncm_serialize_set (ser, obj, name, FALSE);
    g_free (name);
  }

  return obj;
}

static const GVariantType *
_ncm_serialize_gtype_to_gvariant_type (GType t)
{
  switch (t)
  {
    case G_TYPE_CHAR:
    case G_TYPE_UCHAR:
      return G_VARIANT_TYPE_BYTE;

      break;
    case G_TYPE_BOOLEAN:
      return G_VARIANT_TYPE_BOOLEAN;

      break;
    case G_TYPE_INT:
    {
      switch (sizeof (gint))
      {
        case 2:
          return G_VARIANT_TYPE_INT16;

          break;
        case 4:
          return G_VARIANT_TYPE_INT32;

          break;
        case 8:
          return G_VARIANT_TYPE_INT64;

          break;
        default:
          g_error ("Unknown gint size %"G_GSIZE_FORMAT ".", sizeof (gint));
          break;
      }

      break;
    }
    case G_TYPE_UINT:
    {
      switch (sizeof (guint))
      {
        case 2:
          return G_VARIANT_TYPE_UINT16;

          break;
        case 4:
          return G_VARIANT_TYPE_UINT32;

          break;
        case 8:
          return G_VARIANT_TYPE_UINT64;

          break;
        default:
          g_error ("Unknown gint size %"G_GSIZE_FORMAT ".", sizeof (guint));
          break;
      }

      break;
    }
    case G_TYPE_LONG:
    {
      switch (sizeof (glong))
      {
        case 2:
          return G_VARIANT_TYPE_INT16;

          break;
        case 4:
          return G_VARIANT_TYPE_INT32;

          break;
        case 8:
          return G_VARIANT_TYPE_INT64;

          break;
        default:
          g_error ("Unknown gint size %"G_GSIZE_FORMAT ".", sizeof (glong));
          break;
      }

      break;
    }
    case G_TYPE_ULONG:
    {
      switch (sizeof (gulong))
      {
        case 2:
          return G_VARIANT_TYPE_UINT16;

          break;
        case 4:
          return G_VARIANT_TYPE_UINT32;

          break;
        case 8:
          return G_VARIANT_TYPE_UINT64;

          break;
        default:
          g_error ("Unknown gint size %"G_GSIZE_FORMAT ".", sizeof (gulong));
          break;
      }

      break;
    }
    case G_TYPE_INT64:
      return G_VARIANT_TYPE_INT64;

      break;
    case G_TYPE_UINT64:
      return G_VARIANT_TYPE_UINT64;

      break;
    case G_TYPE_FLOAT:
    case G_TYPE_DOUBLE:
      return G_VARIANT_TYPE_DOUBLE;

      break;
    case G_TYPE_STRING:
      return G_VARIANT_TYPE_STRING;

      break;
    case G_TYPE_VARIANT:
      return G_VARIANT_TYPE_VARIANT;

      break;
    default:
      return NULL;

      break;
  }
}

/**
 * ncm_serialize_gvalue_to_gvariant:
 * @ser: a #NcmSerialize.
 * @val: a #GValue.
 *
 * Converts a #GValue to a #GVariant.
 *
 * Returns: (transfer full): A #GVariant convertion of @val.
 */
GVariant *
ncm_serialize_gvalue_to_gvariant (NcmSerialize *ser, GValue *val)
{
  GType t                      = G_VALUE_TYPE (val);
  GType fund_t                 = G_TYPE_FUNDAMENTAL (t);
  const GVariantType *var_type = _ncm_serialize_gtype_to_gvariant_type (fund_t);
  GVariant *var                = NULL;

  if (var_type == NULL)
  {
    switch (fund_t)
    {
      case G_TYPE_OBJECT:
      {
        GObject *nest_obj = g_value_get_object (val);

        if (nest_obj != NULL)
          var = ncm_serialize_to_variant (ser, nest_obj);

        break;
      }
      case G_TYPE_ENUM:
        var = g_variant_ref_sink (g_variant_new_int32 (g_value_get_enum (val)));
        break;
      case G_TYPE_FLAGS:
        var = g_variant_ref_sink (g_variant_new_uint32 (g_value_get_flags (val)));
        break;
      case G_TYPE_BOXED:
      {
        if (g_type_is_a (t, NCM_TYPE_DTUPLE2))
        {
          NcmDTuple2 *dtuple2 = g_value_get_boxed (val);

          if (dtuple2 != NULL)
            var = ncm_dtuple2_serialize (dtuple2);
        }
        else if (g_type_is_a (t, NCM_TYPE_DTUPLE3))
        {
          NcmDTuple3 *dtuple3 = g_value_get_boxed (val);

          if (dtuple3 != NULL)
            var = ncm_dtuple3_serialize (dtuple3);
        }
        else if (g_type_is_a (t, NCM_TYPE_OBJ_ARRAY))
        {
          NcmObjArray *oa = g_value_get_boxed (val);

          if (oa != NULL)
            var = ncm_serialize_array_to_variant (ser, oa);
        }
        else if (g_type_is_a (t, NCM_TYPE_OBJ_DICT_STR))
        {
          NcmObjDictStr *dict = g_value_get_boxed (val);

          if (dict != NULL)
            var = ncm_serialize_dict_str_to_variant (ser, dict);
        }
        else if (g_type_is_a (t, NCM_TYPE_OBJ_DICT_INT))
        {
          NcmObjDictStr *dict = g_value_get_boxed (val);

          if (dict != NULL)
            var = ncm_serialize_dict_int_to_variant (ser, dict);
        }
        else if (g_type_is_a (t, NCM_TYPE_VAR_DICT))
        {
          NcmVarDict *dict = g_value_get_boxed (val);

          if (dict != NULL)
            var = ncm_serialize_var_dict_to_variant (ser, dict);
        }
        else if (g_type_is_a (t, G_TYPE_STRV))
        {
          const gchar * const *strv = g_value_get_boxed (val);

          if (strv != NULL)
            var = g_variant_ref_sink (g_variant_new_strv (strv, -1));
        }
        else
        {
          g_error ("Cannot convert GValue '%s' to GVariant.", g_type_name (t));
        }

        break;
      }
      default:
        g_error ("Cannot convert GValue '%s' to GVariant.", g_type_name (t));
        break;
    }
  }
  else if (t == G_TYPE_VARIANT)
  {
    var = g_value_dup_variant (val);
  }
  else if (t == G_TYPE_STRING)
  {
    const gchar *str = g_value_get_string (val);

    if (str != NULL)
      var = g_variant_ref_sink (g_variant_new_string (str));
  }
  else
  {
    var = g_dbus_gvalue_to_gvariant (val, var_type);
  }

  return var;
}

/**
 * ncm_serialize_to_variant:
 * @ser: a #NcmSerialize.
 * @obj: a #GObject.
 *
 * Serialize the @obj to a @GVariant representation.
 *
 * Returns: (transfer full): A #GVariant dictionary describing the @obj.
 */
GVariant *
ncm_serialize_to_variant (NcmSerialize *ser, GObject *obj)
{
  GVariant *ser_var;
  gchar *obj_name   = g_strdup (G_OBJECT_TYPE_NAME (obj));
  gchar *saved_name = NULL;

  if (ncm_serialize_contain_instance (ser, obj))
  {
    gchar *ni_name = ncm_serialize_peek_name (ser, obj);
    gchar *fname   = g_strdup_printf ("%s[%s]", obj_name, ni_name);

    /*printf ("# Found instante %p at ptr_name %s.\n", obj, fname);*/
    ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_TYPE, fname, NULL));
    g_free (fname);
  }
  else if (g_hash_table_lookup_extended (ser->saved_ptr_name, obj, NULL, (gpointer *) &saved_name))
  {
    gchar *fname = g_strdup_printf ("%s[%s]", obj_name, saved_name);

    /*printf ("# Found instante %p at saved_ptr_name %s.\n", obj, fname);*/
    ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_TYPE, fname, NULL));
    g_free (fname);
  }
  else
  {
    GObjectClass *klass = G_OBJECT_GET_CLASS (obj);
    guint n_properties, i;
    GParamSpec **prop = g_object_class_list_properties (klass, &n_properties);
    gchar *name       = NULL;

    if (ser->opts & NCM_SERIALIZE_OPT_AUTONAME_SER)
    {
      gchar *tmp = obj_name;

      name = g_strdup_printf (NCM_SERIALIZE_AUTOSAVE_NAME NCM_SERIALIZE_AUTOSAVE_NFORMAT, ser->autosave_count);

      obj_name = g_strdup_printf ("%s[%s]", tmp, name);
      g_free (tmp);

      ser->autosave_count++;
    }

    if (n_properties == 0)
    {
      g_free (prop);
      ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_TYPE, obj_name, NULL));
    }
    else
    {
      GVariantBuilder b;
      GVariant *params;

      g_variant_builder_init (&b, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE));

      for (i = 0; i < n_properties; i++)
      {
        GVariant *var = NULL;
        GValue val    = G_VALUE_INIT;

        if ((prop[i]->flags & G_PARAM_READWRITE) != G_PARAM_READWRITE)
          continue;

        g_value_init (&val, prop[i]->value_type);
        g_object_get_property (obj, prop[i]->name, &val);

        var = ncm_serialize_gvalue_to_gvariant (ser, &val);

        if (var == NULL)
        {
          g_value_unset (&val);
          continue;
        }

        g_variant_builder_add (&b, NCM_SERIALIZE_PROPERTY_TYPE, prop[i]->name, var);

        g_variant_unref (var);
        g_value_unset (&val);
      }

      params = g_variant_builder_end (&b);

      g_free (prop);

      ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_FORMAT, obj_name, params));
    }

    if (name != NULL)
    {
      _ncm_serialize_save_ser (ser, name, obj, ser_var);
      g_free (name);
    }
  }

  g_free (obj_name);

  return ser_var;
}

/**
 * ncm_serialize_array_to_variant:
 * @ser: a #NcmSerialize
 * @oa: a #NcmObjArray
 *
 * Serializes a #NcmObjArray to a #GVariant.
 *
 * Returns: (transfer full): the serialized #GVariant.
 */
GVariant *
ncm_serialize_array_to_variant (NcmSerialize *ser, NcmObjArray *oa)
{
  GVariantBuilder *builder;
  GVariant *var;
  guint i;

  builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_ARRAY_TYPE));

  for (i = 0; i < oa->len; i++)
  {
    GVariant *cvar = ncm_serialize_to_variant (ser, ncm_obj_array_peek (oa, i));

    g_variant_builder_add_value (builder, cvar);
    g_variant_unref (cvar);
  }

  var = g_variant_ref_sink (g_variant_builder_end (builder));

  g_variant_builder_unref (builder);

  return var;
}

/**
 * ncm_serialize_dict_str_to_variant:
 * @ser: a #NcmSerialize
 * @ods: a #NcmObjDictStr
 *
 * Serializes a #NcmObjDictStr to a #GVariant.
 *
 * Returns: (transfer full): the serialized #GVariant.
 */
GVariant *
ncm_serialize_dict_str_to_variant (NcmSerialize *ser, NcmObjDictStr *ods)
{
  GVariantBuilder *builder;
  GVariant *var;
  GHashTableIter iter;
  gchar *key;
  GObject *val;

  builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_STR_TYPE));

  g_hash_table_iter_init (&iter, ods);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &val))
  {
    GVariant *key_var = g_variant_new_string (key);
    GVariant *cvar    = ncm_serialize_to_variant (ser, val);
    GVariant *entry   = g_variant_new_dict_entry (key_var, cvar);

    g_variant_builder_add_value (builder, entry);
    g_variant_unref (cvar);
  }

  var = g_variant_ref_sink (g_variant_builder_end (builder));

  g_variant_builder_unref (builder);

  return var;
}

/**
 * ncm_serialize_dict_int_to_variant:
 * @ser: a #NcmSerialize
 * @odi: a #NcmObjDictInt
 *
 * Serializes a #NcmObjDictInt to a #GVariant.
 *
 * Returns: (transfer full): the serialized #GVariant.
 */
GVariant *
ncm_serialize_dict_int_to_variant (NcmSerialize *ser, NcmObjDictInt *odi)
{
  GVariantBuilder *builder;
  GVariant *var;
  GHashTableIter iter;
  gint *key;
  GObject *val;

  builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_INT_TYPE));

  g_hash_table_iter_init (&iter, odi);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &val))
  {
    GVariant *key_var = g_variant_new_int32 (*key);
    GVariant *cvar    = ncm_serialize_to_variant (ser, val);
    GVariant *entry   = g_variant_new_dict_entry (key_var, cvar);

    g_variant_builder_add_value (builder, entry);
    g_variant_unref (cvar);
  }

  var = g_variant_ref_sink (g_variant_builder_end (builder));

  g_variant_builder_unref (builder);

  return var;
}

/**
 * ncm_serialize_var_dict_to_variant:
 * @ser: a #NcmSerialize
 * @vd: a #NcmVarDict
 *
 * Serializes a #NcmVarDict to a #GVariant.
 *
 * Returns: (transfer full): the serialized #GVariant.
 */
GVariant *
ncm_serialize_var_dict_to_variant (NcmSerialize *ser, NcmVarDict *vd)
{
  GVariantBuilder builder;
  GVariant *var;
  GHashTableIter iter;
  gchar *key;
  GVariant *val;

  g_variant_builder_init (&builder, G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE));

  g_hash_table_iter_init (&iter, vd);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &val))
  {
    GVariant *key_var = g_variant_new_string (key);
    GVariant *val_var = g_variant_new_variant (val);
    GVariant *entry   = g_variant_new_dict_entry (key_var, val_var);

    g_variant_builder_add_value (&builder, entry);
  }

  var = g_variant_ref_sink (g_variant_builder_end (&builder));

  return var;
}

#ifdef HAVE_LIBFYAML
static struct fy_node *_ncm_serialize_to_yaml_node (NcmSerialize *ser, struct fy_document *doc, GVariant *var_obj);

#endif /* HAVE_LIBFYAML */

/**
 * ncm_serialize_variant_to_yaml:
 * @ser: a #NcmSerialize
 * @var_obj: a #GObject serialized to a #GVariant
 *
 * Converts a #GObject serialized to a #GVariant to a YAML string.
 *
 * Returns: A pointer to the YAML string representation of the @var_obj.
 */
gchar *
ncm_serialize_variant_to_yaml (NcmSerialize *ser, GVariant *var_obj)
{
#ifdef HAVE_LIBFYAML
  g_assert (var_obj != NULL);
  g_assert (g_variant_is_of_type (var_obj, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)));

  {
    struct fy_document *doc = fy_document_create (NULL);
    struct fy_node *root    = _ncm_serialize_to_yaml_node (ser, doc, var_obj);
    gchar *yaml             = NULL;

    fy_document_set_root (doc, root);
    yaml = fy_emit_document_to_string (doc, FYECF_DEFAULT | FYECF_WIDTH_INF);
    fy_document_destroy (doc);

    return yaml;
  }
#else
  g_error ("ncm_serialize_variant_to_yaml: libfyaml not available.");

  return NULL;

#endif /* HAVE_LIBFYAML */
}

#ifdef HAVE_LIBFYAML

static struct fy_node *

_ncm_serialize_to_yaml_node (NcmSerialize *ser, struct fy_document *doc, GVariant *var_obj)
{
  GVariant *obj_name_var     = g_variant_get_child_value (var_obj, 0);
  GVariant *params_var       = g_variant_get_child_value (var_obj, 1);
  struct fy_node *root       = fy_node_create_mapping (doc);
  const gchar *obj_full_name = g_variant_get_string (obj_name_var, NULL);
  guint n_properties         = g_variant_n_children (params_var);
  GMatchInfo *match_info     = NULL;
  gchar *name;
  gchar *anchor;
  GType gtype;

  g_assert (doc != NULL);
  g_assert (root != NULL);

  if (g_regex_match (ser->is_named_regex, obj_full_name, 0, &match_info))
  {
    name   = g_match_info_fetch (match_info, 1);
    anchor = g_match_info_fetch (match_info, 2);
  }
  else
  {
    name   = g_strdup (obj_full_name);
    anchor = NULL;
  }

  gtype = g_type_from_name (name);

  g_match_info_free (match_info);

  if (gtype == 0)
    g_error ("_ncm_serialize_to_yaml_node: object `%s' is not registered.", name);

  g_assert (g_variant_is_of_type (params_var, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE)));

  if ((n_properties == 0) && (anchor != NULL))
  {
    fy_node_mapping_append (root,
                            fy_node_create_scalar_copy (doc, name, FY_NT),
                            fy_node_create_alias_copy (doc, anchor, FY_NT));
  }
  else
  {
    struct fy_node *properties = fy_node_create_mapping (doc);
    struct fy_node *root_key   = fy_node_create_scalar_copy (doc, name, FY_NT);

    if (anchor)
    {
      gint rc = fy_node_set_anchor (root_key, g_strdup (anchor), FY_NT);

      g_assert (rc == 0);
    }

    fy_node_mapping_append (root,
                            root_key,
                            properties);

    {
      GVariantIter *p_iter = g_variant_iter_new (params_var);
      GVariant *var        = NULL;
      guint i;

      while ((var = g_variant_iter_next_value (p_iter)))
      {
        GVariant *var_key      = g_variant_get_child_value (var, 0);
        GVariant *var_val      = g_variant_get_child_value (var, 1);
        GVariant *val          = g_variant_get_variant (var_val);
        const gchar *prop_name = g_variant_get_string (var_key, NULL);
        struct fy_node *value  = NULL;

        if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_ARRAY_TYPE)))
        {
          const guint n = g_variant_n_children (val);

          value = fy_node_create_sequence (doc);

          for (i = 0; i < n; i++)
          {
            GVariant *cvar = g_variant_get_child_value (val, i);

            fy_node_sequence_append (value, _ncm_serialize_to_yaml_node (ser, doc, cvar));
            g_variant_unref (cvar);
          }
        }
        else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_STR_TYPE)))
        {
          const guint n = g_variant_n_children (val);

          value = fy_node_create_mapping (doc);

          for (i = 0; i < n; i++)
          {
            GVariant *cvar     = g_variant_get_child_value (val, i);
            GVariant *cvar_key = g_variant_get_child_value (cvar, 0);
            GVariant *cvar_val = g_variant_get_child_value (cvar, 1);

            fy_node_mapping_append (value,
                                    fy_node_create_scalar_copy (doc, g_variant_get_string (cvar_key, NULL), FY_NT),
                                    _ncm_serialize_to_yaml_node (ser, doc, cvar_val));

            g_variant_unref (cvar_key);
            g_variant_unref (cvar_val);
            g_variant_unref (cvar);
          }
        }
        else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_DICT_INT_TYPE)))
        {
          const guint n = g_variant_n_children (val);

          value = fy_node_create_mapping (doc);

          for (i = 0; i < n; i++)
          {
            GVariant *cvar     = g_variant_get_child_value (val, i);
            GVariant *cvar_key = g_variant_get_child_value (cvar, 0);
            GVariant *cvar_val = g_variant_get_child_value (cvar, 1);

            fy_node_mapping_append (value,
                                    fy_node_create_scalarf (doc, "%d", g_variant_get_int32 (cvar_key)),
                                    _ncm_serialize_to_yaml_node (ser, doc, cvar_val));

            g_variant_unref (cvar_key);
            g_variant_unref (cvar_val);
            g_variant_unref (cvar);
          }
        }
        else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_VAR_DICT_TYPE)))
        {
          const guint n = g_variant_n_children (val);

          value = fy_node_create_mapping (doc);

          for (i = 0; i < n; i++)
          {
            GVariant *cvar     = g_variant_get_child_value (val, i);
            GVariant *cvar_key = g_variant_get_child_value (cvar, 0);
            GVariant *cvar_val = g_variant_get_child_value (cvar, 1);
            GVariant *val      = g_variant_get_variant (cvar_val);

            fy_node_mapping_append (value,
                                    fy_node_create_scalar_copy (doc, g_variant_get_string (cvar_key, NULL), FY_NT),
                                    fy_node_build_from_malloc_string (doc, g_variant_print (val, FALSE), FY_NT));

            g_variant_unref (cvar_key);
            g_variant_unref (cvar_val);
            g_variant_unref (cvar);
            g_variant_unref (val);
          }
        }
        else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)))
        {
          value = _ncm_serialize_to_yaml_node (ser, doc, val);
        }
        else
        {
          gchar *str = g_variant_print (val, FALSE);

          value = fy_node_build_from_malloc_string (doc, str, FY_NT);
        }

        fy_node_mapping_append (properties,
                                fy_node_create_scalar_copy (doc, prop_name, FY_NT),
                                value);

        g_variant_unref (var_key);
        g_variant_unref (var_val);
        g_variant_unref (val);
        g_variant_unref (var);
      }

      g_variant_iter_free (p_iter);
    }
  }

  g_clear_pointer (&name, g_free);
  g_clear_pointer (&anchor, g_free);
  g_variant_unref (obj_name_var);
  g_variant_unref (params_var);

  return root;
}

#endif /* HAVE_LIBFYAML */

/**
 * ncm_serialize_to_string:
 * @ser: a #NcmSerialize.
 * @obj: a #GObject.
 * @valid_variant: whether to use a valid #GVariant representation.
 *
 * Serialize the object @obj to a string.
 *
 * Returns: (transfer full): A string containing the serialized version of @obj.
 */
gchar *
ncm_serialize_to_string (NcmSerialize *ser, GObject *obj, gboolean valid_variant)
{
  GVariant *ser_var = ncm_serialize_to_variant (ser, obj);
  gchar *serstr     = NULL;

  if (valid_variant)
  {
    serstr = g_variant_print (ser_var, TRUE);
  }
  else
  {
    GVariant *params = NULL;
    gchar *obj_name  = NULL;
    gchar *params_str;

    g_variant_get (ser_var, NCM_SERIALIZE_OBJECT_FORMAT, &obj_name, &params);

    if (g_variant_n_children (params) == 0)
    {
      serstr = obj_name;
    }
    else
    {
      params_str = g_variant_print (params, TRUE);
      serstr     = g_strdup_printf ("%s%s", obj_name, params_str);
      g_free (params_str);
      g_free (obj_name);
    }

    g_variant_unref (params);
  }

  g_variant_unref (ser_var);

  return serstr;
}

/**
 * ncm_serialize_to_yaml:
 * @ser: a #NcmSerialize
 * @obj: a #GObject
 *
 * Serialize the object @obj to a YAML string.
 *
 * Returns: (transfer full): A YAML string containing the serialized version of @obj.
 */
gchar *
ncm_serialize_to_yaml (NcmSerialize *ser, GObject *obj)
{
  GVariant *ser_var = ncm_serialize_to_variant (ser, obj);
  gchar *yaml_str   = ncm_serialize_variant_to_yaml (ser, ser_var);

  g_variant_unref (ser_var);

  return yaml_str;
}

/**
 * ncm_serialize_array_to_yaml:
 * @ser: a #NcmSerialize
 * @oa: a #NcmObjArray
 *
 * Serialize the #NcmObjArray @oa to a YAML string.
 *
 * Returns: (transfer full): A YAML string containing the serialized version of @oa.
 */
gchar *
ncm_serialize_array_to_yaml (NcmSerialize *ser, NcmObjArray *oa)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_create (NULL);
  struct fy_node *root    = fy_node_create_sequence (doc);
  gchar *yaml_str         = NULL;
  guint i;

  for (i = 0; i < ncm_obj_array_len (oa); i++)
  {
    GObject *obj      = ncm_obj_array_peek (oa, i);
    GVariant *obj_var = ncm_serialize_to_variant (ser, obj);

    fy_node_sequence_append (root, _ncm_serialize_to_yaml_node (ser, doc, obj_var));
    g_variant_unref (obj_var);
  }

  fy_document_set_root (doc, root);
  yaml_str = fy_emit_document_to_string (doc, FYECF_DEFAULT | FYECF_WIDTH_INF);
  fy_document_destroy (doc);

  return yaml_str;

#else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_array_to_yaml: libfyaml not available.");

  return NULL;

#endif
}

/**
 * ncm_serialize_dict_str_to_yaml:
 * @ser: a #NcmSerialize
 * @ods: a #NcmObjDictStr
 *
 * Serialize the #NcmObjDictStr @ods to a YAML string.
 *
 * Returns: (transfer full): A YAML string containing the serialized version of @ods.
 */
gchar *
ncm_serialize_dict_str_to_yaml (NcmSerialize *ser, NcmObjDictStr *ods)
{
#ifdef HAVE_LIBFYAML
  struct fy_document *doc = fy_document_create (NULL);
  struct fy_node *root    = fy_node_create_mapping (doc);
  gchar *yaml_str         = NULL;
  GHashTableIter iter;
  gchar *key;
  GObject *value;

  g_hash_table_iter_init (&iter, ods);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
  {
    GVariant *obj_var = ncm_serialize_to_variant (ser, value);

    fy_node_mapping_append (root,
                            fy_node_create_scalar_copy (doc, key, FY_NT),
                            _ncm_serialize_to_yaml_node (ser, doc, obj_var));
    g_variant_unref (obj_var);
  }

  fy_document_set_root (doc, root);
  yaml_str = fy_emit_document_to_string (doc, FYECF_DEFAULT | FYECF_WIDTH_INF);
  fy_document_destroy (doc);

  return yaml_str;

#else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_dict_str_to_yaml: libfyaml not available.");

  return NULL;

#endif
}

/**
 * ncm_serialize_dict_int_to_yaml:
 * @ser: a #NcmSerialize
 * @odi: a #NcmObjDictInt
 *
 * Serialize the #NcmObjDictInt @odi to a YAML string.
 *
 * Returns: (transfer full): A YAML string containing the serialized version of @odi.
 */
gchar *
ncm_serialize_dict_int_to_yaml (NcmSerialize *ser, NcmObjDictInt *odi)
{
#ifdef HAVE_LIBFYAML

  struct fy_document *doc = fy_document_create (NULL);
  struct fy_node *root    = fy_node_create_mapping (doc);
  gchar *yaml_str         = NULL;
  GHashTableIter iter;
  gint *key;
  GObject *value;

  g_hash_table_iter_init (&iter, odi);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
  {
    GVariant *obj_var = ncm_serialize_to_variant (ser, value);

    fy_node_mapping_append (root,
                            fy_node_create_scalarf (doc, "%d", *key),
                            _ncm_serialize_to_yaml_node (ser, doc, obj_var));
    g_variant_unref (obj_var);
  }

  fy_document_set_root (doc, root);
  yaml_str = fy_emit_document_to_string (doc, FYECF_DEFAULT | FYECF_WIDTH_INF);
  fy_document_destroy (doc);

  return yaml_str;

#else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_dict_int_to_yaml: libfyaml not available.");

  return NULL;

#endif
}

/**
 * ncm_serialize_var_dict_to_yaml:
 * @ser: a #NcmSerialize
 * @dict: a #NcmVarDict
 *
 * Serialize the #NcmVarDict @dict to a YAML string.
 *
 * Returns: (transfer full): A YAML string containing the serialized version of @dict.
 */
gchar *
ncm_serialize_var_dict_to_yaml (NcmSerialize *ser, NcmVarDict *dict)
{
#ifdef HAVE_LIBFYAML
  struct fy_document *doc = fy_document_create (NULL);
  struct fy_node *root    = fy_node_create_mapping (doc);
  gchar *yaml_str         = NULL;
  GHashTableIter iter;
  gchar *key;
  GVariant *value;

  g_hash_table_iter_init (&iter, dict);

  while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
  {
    fy_node_mapping_append (root,
                            fy_node_create_scalar_copy (doc, key, FY_NT),
                            fy_node_build_from_malloc_string (doc, g_variant_print (value, FALSE), FY_NT));
  }

  fy_document_set_root (doc, root);
  yaml_str = fy_emit_document_to_string (doc, FYECF_DEFAULT | FYECF_WIDTH_INF);
  fy_document_destroy (doc);

  return yaml_str;

#else /* HAVE_LIBFYAML */

  g_error ("ncm_serialize_var_dict_to_yaml: libfyaml not available.");

  return NULL;

#endif
}

/**
 * ncm_serialize_to_file:
 * @ser: a #NcmSerialize
 * @obj: a #GObject
 * @filename: File where to save the serialized version of the object
 *
 * Serializes @obj and saves the string in @filename.
 *
 */
void
ncm_serialize_to_file (NcmSerialize *ser, GObject *obj, const gchar *filename)
{
  GError *error  = NULL;
  gchar *obj_ser = ncm_serialize_to_string (ser, obj, TRUE);
  gsize length   = strlen (obj_ser);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, obj_ser, length, &error))
    g_error ("ncm_serialize_to_file: cannot save to file %s: %s",
             filename, error->message);

  g_free (obj_ser);
}

/**
 * ncm_serialize_to_binfile:
 * @ser: a #NcmSerialize
 * @obj: a #GObject
 * @filename: File where to save the serialized version of the object
 *
 * Serializes @obj and saves the binary in @filename.
 *
 */
void
ncm_serialize_to_binfile (NcmSerialize *ser, GObject *obj, const gchar *filename)
{
  GError *error     = NULL;
  GVariant *obj_ser = ncm_serialize_to_variant (ser, obj);
  gsize length      = g_variant_get_size (obj_ser);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, g_variant_get_data (obj_ser), length, &error))
    g_error ("ncm_serialize_to_binfile: cannot save to file %s: %s",
             filename, error->message);

  g_variant_unref (obj_ser);
}

/**
 * ncm_serialize_var_dict_to_variant_file:
 * @ser: a #NcmSerialize
 * @vd: a #NcmVarDict
 * @filename: File where to save the serialized version of the object
 * @binary: whether to save the variant in binary format
 *
 * Serializes @vd and saves the variant string in @filename.
 *
 */
void
ncm_serialize_var_dict_to_variant_file (NcmSerialize *ser, NcmVarDict *vd, const gchar *filename, gboolean binary)
{
  GError *error = NULL;
  GVariant *var = ncm_serialize_var_dict_to_variant (ser, vd);
  gsize length  = g_variant_get_size (var);

  g_assert (filename != NULL);

  if (binary)
  {
    if (!g_file_set_contents (filename, g_variant_get_data (var), length, &error))
      g_error ("ncm_serialize_var_dict_to_variant_file: cannot save to file %s: %s",
               filename, error->message);
  }
  else
  {
    gchar *variant_string = g_variant_print (var, TRUE);

    if (!g_file_set_contents (filename, variant_string, strlen (variant_string), &error))
      g_error ("ncm_serialize_var_dict_to_variant_file: cannot save to file %s: %s",
               filename, error->message);

    g_free (variant_string);
  }

  g_variant_unref (var);
}

/**
 * ncm_serialize_to_yaml_file:
 * @ser: a #NcmSerialize
 * @obj: a #GObject
 * @filename: File where to save the serialized version of the object
 *
 * Serializes @obj and saves the YAML string in @filename.
 *
 */
void
ncm_serialize_to_yaml_file (NcmSerialize *ser, GObject *obj, const gchar *filename)
{
  GError *error = NULL;
  gchar *yaml   = ncm_serialize_to_yaml (ser, obj);
  gsize length  = strlen (yaml);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, yaml, length, &error))
    g_error ("ncm_serialize_to_yaml_file: cannot save to file %s: %s",
             filename, error->message);

  g_free (yaml);
}

/**
 * ncm_serialize_array_to_key_file:
 * @ser: a #NcmSerialize
 * @oa: a #NcmObjArray
 * @filename: oa filename
 * @save_comment: whether to save comments
 *
 * Saves a #NcmObjArray to a file using a #NcmSerialize and a #GKeyFile.
 *
 */
void
ncm_serialize_array_to_key_file (NcmSerialize *ser, NcmObjArray *oa, const gchar *filename, gboolean save_comment)
{
  GKeyFile *oafile = g_key_file_new ();
  guint i;

  {
    GError *error  = NULL;
    gchar *oa_desc = ncm_cfg_string_to_comment ("Whether NcmObjArray is empty");

    g_key_file_set_boolean (oafile, "NcmObjArray", "empty", oa->len == 0 ? TRUE : FALSE);

    if (save_comment)
      if (!g_key_file_set_comment (oafile, "NcmObjArray", NULL, oa_desc, &error))
        g_error ("ncm_obj_array_save: %s", error->message);


    g_free (oa_desc);
  }

  for (i = 0; i < oa->len; i++)
  {
    GObject *go          = ncm_obj_array_peek (oa, i);
    GObjectClass *oclass = G_OBJECT_GET_CLASS (go);
    GError *error        = NULL;
    gchar *group         = g_strdup_printf (NCM_SERIALIZE_OBJECT_ARRAY_POS_STR ":%d", i);
    GVariant *go_var     = ncm_serialize_to_variant (ser, go);
    GVariant *params     = NULL;
    gchar *obj_name      = NULL;
    guint nparams;

    g_variant_get (go_var, NCM_SERIALIZE_OBJECT_FORMAT, &obj_name, &params);
    nparams = g_variant_n_children (params);

    g_key_file_set_value (oafile, group, NCM_SERIALIZE_OBJECT_ARRAY_OBJ_NAME_STR, obj_name);

    if (nparams != 0)
    {
      GVariantIter iter;
      GVariant *value;
      gchar *key;

      g_variant_iter_init (&iter, params);

      while (g_variant_iter_next (&iter, "{sv}", &key, &value))
      {
        GParamSpec *param_spec = g_object_class_find_property (oclass, key);
        gchar *param_str       = g_variant_print (value, TRUE);

        if (param_spec == NULL)
          g_error ("ncm_obj_array_save: property `%s' not found in object `%s'.", key, obj_name);

        g_key_file_set_value (oafile, group, key, param_str);

        if (save_comment)
        {
          const gchar *blurb = g_param_spec_get_blurb (param_spec);

          if ((blurb != NULL) && (blurb[0] != 0))
          {
            gchar *desc = ncm_cfg_string_to_comment (blurb);

            if (!g_key_file_set_comment (oafile, group, key, desc, &error))
              g_error ("ncm_obj_array_save: %s", error->message);

            g_free (desc);
          }
        }

        g_variant_unref (value);
        g_free (key);
        g_free (param_str);
      }
    }

    g_free (obj_name);
    g_variant_unref (params);
    g_variant_unref (go_var);

    g_free (group);
  }

  {
    GError *error  = NULL;
    gsize len      = 0;
    gchar *oa_data = g_key_file_to_data (oafile, &len, &error);

    if (error != NULL)
      g_error ("Error converting NcmObjArray to configuration file: %s", error->message);

    if (!g_file_set_contents (filename, oa_data, len, &error))
      g_error ("Error saving configuration file to disk: %s", error->message);

    g_free (oa_data);
    g_key_file_free (oafile);
  }
}

/**
 * ncm_serialize_array_to_yaml_file:
 * @ser: a #NcmSerialize
 * @oa: a #NcmObjArray
 * @filename: oa filename
 *
 * Saves a #NcmObjArray to a file using a #NcmSerialize and a YAML string.
 *
 */
void
ncm_serialize_array_to_yaml_file (NcmSerialize *ser, NcmObjArray *oa, const gchar *filename)
{
  GError *error = NULL;
  gchar *yaml   = ncm_serialize_array_to_yaml (ser, oa);
  gsize length  = strlen (yaml);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, yaml, length, &error))
    g_error ("ncm_serialize_array_to_yaml_file: cannot save to file %s: %s",
             filename, error->message);

  g_free (yaml);
}

/**
 * ncm_serialize_dict_str_to_yaml_file:
 * @ser: a #NcmSerialize
 * @ods: a #NcmObjDictStr
 * @filename: ods filename
 *
 * Saves a #NcmObjDictStr to a file using a #NcmSerialize and a YAML string.
 *
 */
void
ncm_serialize_dict_str_to_yaml_file (NcmSerialize *ser, NcmObjDictStr *ods, const gchar *filename)
{
  GError *error = NULL;
  gchar *yaml   = ncm_serialize_dict_str_to_yaml (ser, ods);
  gsize length  = strlen (yaml);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, yaml, length, &error))
    g_error ("ncm_serialize_dict_str_to_yaml_file: cannot save to file %s: %s",
             filename, error->message);

  g_free (yaml);
}

/**
 * ncm_serialize_dict_int_to_yaml_file:
 * @ser: a #NcmSerialize
 * @odi: a #NcmObjDictInt
 * @filename: odi filename
 *
 * Saves a #NcmObjDictInt to a file using a #NcmSerialize and a YAML string.
 *
 */
void
ncm_serialize_dict_int_to_yaml_file (NcmSerialize *ser, NcmObjDictInt *odi, const gchar *filename)
{
  GError *error = NULL;
  gchar *yaml   = ncm_serialize_dict_int_to_yaml (ser, odi);
  gsize length  = strlen (yaml);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, yaml, length, &error))
    g_error ("ncm_serialize_dict_int_to_yaml_file: cannot save to file %s: %s",
             filename, error->message);

  g_free (yaml);
}

/**
 * ncm_serialize_var_dict_to_yaml_file:
 * @ser: a #NcmSerialize
 * @vd: a #NcmVarDict
 * @filename: vd filename
 *
 * Saves a #NcmVarDict to a file using a #NcmSerialize and a YAML string.
 *
 *
 */
void
ncm_serialize_var_dict_to_yaml_file (NcmSerialize *ser, NcmVarDict *vd, const gchar *filename)
{
  GError *error = NULL;
  gchar *yaml   = ncm_serialize_var_dict_to_yaml (ser, vd);
  gsize length  = strlen (yaml);

  g_assert (filename != NULL);

  if (!g_file_set_contents (filename, yaml, length, &error))
    g_error ("ncm_serialize_var_dict_to_yaml_file: cannot save to file %s: %s",
             filename, error->message);

  g_free (yaml);
}

/**
 * ncm_serialize_dup_obj:
 * @ser: a #NcmSerialize.
 * @obj: a #GObject.
 *
 * Duplicates @obj by serializing and deserializing a new object.
 *
 * Returns: (transfer full): A duplicate of @obj.
 */
GObject *
ncm_serialize_dup_obj (NcmSerialize *ser, GObject *obj)
{
  GVariant *var = ncm_serialize_to_variant (ser, obj);
  GObject *dup  = ncm_serialize_from_variant (ser, var);

  g_variant_unref (var);

  return dup;
}

/**
 * ncm_serialize_dup_array:
 * @ser: a #NcmSerialize
 * @oa: a #NcmObjArray
 *
 * Duplicates a #NcmObjArray, all objects are duplicated.
 *
 * Returns: (transfer full): a new #NcmObjArray.
 */
NcmObjArray *
ncm_serialize_dup_array (NcmSerialize *ser, NcmObjArray *oa)
{
  GVariant *var    = ncm_serialize_array_to_variant (ser, oa);
  NcmObjArray *dup = ncm_serialize_array_from_variant (ser, var);

  g_variant_unref (var);

  return dup;
}

static gsize _global_init        = FALSE;
static NcmSerialize *_global_ser = NULL;

/**
 * ncm_serialize_global:
 *
 * Gets the global serialization object, instanciates it if necessary.
 *
 * Returns: (transfer full): The global #NcmSerialize.
 */
NcmSerialize *
ncm_serialize_global (void)
{
  if (g_once_init_enter (&_global_init))
  {
    _global_ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    g_once_init_leave (&_global_init, TRUE);
  }

  return ncm_serialize_ref (_global_ser);
}

/**
 * ncm_serialize_global_reset:
 * @autosave_only: a boolean
 *
 * Releases all objects in global #NcmSerialize and erase
 * all serialized objects.
 *
 */
void
ncm_serialize_global_reset (gboolean autosave_only)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_reset (ser, autosave_only);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_clear_instances:
 * @autosave_only: a boolean
 *
 * Releases all objects in global #NcmSerialize.
 *
 */
void
ncm_serialize_global_clear_instances (gboolean autosave_only)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_clear_instances (ser, autosave_only);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_log_stats:
 *
 * Releases all objects in global #NcmSerialize.
 *
 */
void
ncm_serialize_global_log_stats (void)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_log_stats (ser);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_contain_instance:
 * @obj: (type GObject): a #GObject.
 *
 * Global version of ncm_serialize_contain_instance().
 *
 * Returns: if @obj is already in @ser.
 */
gboolean
ncm_serialize_global_contain_instance (gpointer obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gboolean ret      = ncm_serialize_contain_instance (ser, obj);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_contain_name:
 * @name: an instance name.
 *
 * Global version of ncm_serialize_contain_name().
 *
 * Returns: if there is instance named @name in @ser.
 */
gboolean
ncm_serialize_global_contain_name (const gchar *name)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gboolean ret      = ncm_serialize_contain_name (ser, name);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_count_instances:
 *
 * Global version of ncm_serialize_count_instances().
 *
 * Returns: the number of instances in @ser.
 */
guint
ncm_serialize_global_count_instances (void)
{
  NcmSerialize *ser = ncm_serialize_global ();
  guint ret         = ncm_serialize_count_instances (ser);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_count_saved_serializations:
 *
 * Global version of ncm_serialize_count_saved_serializations().
 *
 * Returns: the number of instances in @ser.
 */
guint
ncm_serialize_global_count_saved_serializations (void)
{
  NcmSerialize *ser = ncm_serialize_global ();
  guint ret         = ncm_serialize_count_saved_serializations (ser);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_get_by_name:
 * @name: an instance name.
 *
 * Global version of ncm_serialize_get_by_name().
 *
 * Returns: (transfer full) (type GObject): Gets the instance named @name or NULL.
 */
gpointer
ncm_serialize_global_get_by_name (const gchar *name)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gpointer ret      = ncm_serialize_get_by_name (ser, name);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_peek_name:
 * @obj: (type GObject): a #GObject.
 *
 * Global version of ncm_serialize_peek_name().
 *
 * Returns: (transfer none): the name of @obj.
 */
gchar *
ncm_serialize_global_peek_name (gpointer obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gchar *ret        = ncm_serialize_peek_name (ser, obj);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_set:
 * @obj: (type GObject): a #GObject.
 * @name: the @obj name.
 * @overwrite: whether to overwrite if there is already an object named @name.
 *
 * Global version of ncm_serialize_set().
 *
 */
void
ncm_serialize_global_set (gpointer obj, const gchar *name, gboolean overwrite)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_set (ser, obj, name, overwrite);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_unset:
 * @obj: (type GObject): a #GObject.
 *
 * Global version of ncm_serialize_unset().
 *
 */
void
ncm_serialize_global_unset (gpointer obj)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_unset (ser, obj);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_remove_ser:
 * @obj: (type GObject): a #GObject.
 *
 * Global version of ncm_serialize_remove_ser().
 *
 */
void
ncm_serialize_global_remove_ser (gpointer obj)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_remove_ser (ser, obj);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_is_named:
 * @serobj: serialized object.
 * @name: (out) (transfer full): object name.
 *
 * Global version of ncm_serialize_is_named().
 *
 * Returns: whether @serobj is a named serialized object.
 */
gboolean
ncm_serialize_global_is_named (const gchar *serobj, gchar **name)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gboolean ret      = ncm_serialize_is_named (ser, serobj, name);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_set_property:
 * @obj: a #GObject.
 * @prop_str: a string containing the parameters to set.
 *
 * Global version of ncm_serialize_set_property().
 *
 */
void
ncm_serialize_global_set_property (GObject *obj, const gchar *prop_str)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_set_property (ser, obj, prop_str);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_set_property_from_key_file:
 * @obj: a #GObject.
 * @prop_file: a #GKeyFile containing the parameters to set.
 *
 * Global version of ncm_serialize_set_property().
 *
 */
void
ncm_serialize_global_set_property_from_key_file (GObject *obj, const gchar *prop_file)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_set_property_from_key_file (ser, obj, prop_file);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_from_variant:
 * @var_obj: A #GVariant containing the serialized version of the object.
 *
 * Global version of ncm_serialize_from_variant().
 *
 * Returns: (transfer full): a new #GObject deserialized from @var_obj.
 */
GObject *
ncm_serialize_global_from_variant (GVariant *var_obj)
{
  GObject *obj;
  NcmSerialize *ser = ncm_serialize_global ();

  obj = ncm_serialize_from_variant (ser, var_obj);
  ncm_serialize_unref (ser);

  return obj;
}

/**
 * ncm_serialize_global_from_string:
 * @obj_ser: String containing the serialized version of the object.
 *
 * Global version of ncm_serialize_from_string().
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_global_from_string (const gchar *obj_ser)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GObject *ret      = ncm_serialize_from_string (ser, obj_ser);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_from_yaml:
 * @yaml_obj: a string containing the serialized version of the object in YAML format
 *
 * Global version of ncm_serialize_from_yaml().
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_global_from_yaml (const gchar *yaml_obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GObject *ret      = ncm_serialize_from_yaml (ser, yaml_obj);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_from_file:
 * @filename: File containing the serialized version of the object.
 *
 * Global version of ncm_serialize_from_file().
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_global_from_file (const gchar *filename)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GObject *ret      = ncm_serialize_from_file (ser, filename);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_from_binfile:
 * @filename: File containing the serialized version of the object.
 *
 * Global version of ncm_serialize_from_binfile().
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_global_from_binfile (const gchar *filename)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GObject *ret      = ncm_serialize_from_binfile (ser, filename);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_from_yaml_file:
 * @filename: File containing the serialized version of the object.
 *
 * Global version of ncm_serialize_from_yaml_file().
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_global_from_yaml_file (const gchar *filename)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GObject *ret      = ncm_serialize_from_yaml_file (ser, filename);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_from_name_params:
 * @obj_name: string containing the object name.
 * @params: a #GVariant containing the object parameters.
 *
 * Global version of ncm_serialize_from_name_params().
 *
 * Returns: (transfer full): A new #GObject.
 */
GObject *
ncm_serialize_global_from_name_params (const gchar *obj_name, GVariant *params)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GObject *ret      = ncm_serialize_from_name_params (ser, obj_name, params);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_gvalue_to_gvariant:
 * @val: a #GValue.
 *
 * Global version of ncm_serialize_gvalue_to_gvariant().
 *
 * Returns: (transfer full): A #GVariant convertion of @val.
 */
GVariant *
ncm_serialize_global_gvalue_to_gvariant (GValue *val)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GVariant *ret     = ncm_serialize_gvalue_to_gvariant (ser, val);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_to_variant:
 * @obj: a #GObject.
 *
 * Global version of ncm_serialize_to_variant().
 *
 * Returns: (transfer full): A #GVariant dictionary describing the @obj.
 */
GVariant *
ncm_serialize_global_to_variant (GObject *obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GVariant *ret     = ncm_serialize_to_variant (ser, obj);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_to_string:
 * @obj: a #GObject.
 * @valid_variant: whether to use a valid #GVariant representation.
 *
 * Global version of ncm_serialize_to_string().
 *
 * Returns: (transfer full): A string containing the serialized version of @obj.
 */
gchar *
ncm_serialize_global_to_string (GObject *obj, gboolean valid_variant)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gchar *ret        = ncm_serialize_to_string (ser, obj, valid_variant);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_to_yaml:
 * @obj: a #GObject
 *
 * Global version of ncm_serialize_to_yaml().
 *
 * Returns: (transfer full): A string containing the serialized version of @obj.
 */
gchar *
ncm_serialize_global_to_yaml (GObject *obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gchar *ret        = ncm_serialize_to_yaml (ser, obj);

  ncm_serialize_unref (ser);

  return ret;
}

/**
 * ncm_serialize_global_to_file:
 * @obj: a #GObject.
 * @filename: File where to save the serialized version of the object
 *
 * Global version of ncm_serialize_to_file().
 *
 */
void
ncm_serialize_global_to_file (GObject *obj, const gchar *filename)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_to_file (ser, obj, filename);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_to_binfile:
 * @obj: a #GObject.
 * @filename: File where to save the serialized version of the object
 *
 * Global version of ncm_serialize_to_binfile().
 *
 */
void
ncm_serialize_global_to_binfile (GObject *obj, const gchar *filename)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_to_binfile (ser, obj, filename);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_to_yaml_file:
 * @obj: a #GObject.
 * @filename: File where to save the serialized version of the object
 *
 * Global version of ncm_serialize_to_yaml_file().
 *
 */
void
ncm_serialize_global_to_yaml_file (GObject *obj, const gchar *filename)
{
  NcmSerialize *ser = ncm_serialize_global ();

  ncm_serialize_to_yaml_file (ser, obj, filename);
  ncm_serialize_unref (ser);
}

/**
 * ncm_serialize_global_dup_obj:
 * @obj: a #GObject.
 *
 * Global version of ncm_serialize_dup_obj().
 *
 * Returns: (transfer full): A duplicate of @obj.
 */
GObject *
ncm_serialize_global_dup_obj (GObject *obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  GVariant *var     = ncm_serialize_to_variant (ser, obj);
  GObject *dup      = ncm_serialize_from_variant (ser, var);

  g_variant_unref (var);
  ncm_serialize_unref (ser);

  return dup;
}

/**
 * ncm_serialize_global_variant_to_yaml:
 * @var_obj: a #GObject serialized to a #GVariant
 *
 * Global version of ncm_serialize_variant_to_yaml().
 * Converts a #GObject serialized to a #GVariant to a YAML string.
 *
 * Returns: A pointer to the YAML string representation of the @var_obj.
 */
gchar *
ncm_serialize_global_variant_to_yaml (GVariant *var_obj)
{
  NcmSerialize *ser = ncm_serialize_global ();
  gchar *ret        = ncm_serialize_variant_to_yaml (ser, var_obj);

  ncm_serialize_unref (ser);

  return ret;
}


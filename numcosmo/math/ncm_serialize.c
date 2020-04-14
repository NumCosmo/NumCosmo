/***************************************************************************
 *            ncm_serialize.c
 *
 *  Mon August 26 13:38:17 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_serialize.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_obj_array.h"
#include "math/ncm_vector.h"
#include "math/ncm_matrix.h"
#include "ncm_enum_types.h"

#include <gio/gio.h>

enum
{
  PROP_0,
  PROP_OPTS,
};

G_DEFINE_TYPE (NcmSerialize, ncm_serialize, G_TYPE_OBJECT);

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
                                               (GDestroyNotify)&g_variant_unref);
  
  ser->is_named_regex  = g_regex_new ("^\\s*([A-Za-z][A-Za-z0-9\\+\\_]+)\\s*\\[([A-Za-z0-9\\:]+)\\]\\s*$", 0, 0, &error);
  ser->parse_obj_regex = g_regex_new ("^\\s*([A-Za-z][A-Za-z0-9\\+\\_]+\\s*(?:\\[[A-Za-z0-9\\:]+\\])?)\\s*([\\{]?.*[\\}]?)\\s*$", 0, 0, &error);
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
    return g_hash_table_lookup (ser->ptr_name, obj);
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
        g_error ("ncm_serialize_set: named instance already present, set overwrite to true if you want it replaced.");
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
  if (g_hash_table_lookup_extended (ser->saved_ptr_name, obj, NULL, (gpointer *)&saved_name))
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

  //printf ("Saving: ``%s'' <=> %s\n", name, g_variant_print (ser_var, TRUE));

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
    obj_name_str       = g_variant_dup_string (obj_name, NULL);
    
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
  GError *error = NULL;
  GVariant *params = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE),
                                      prop_str, NULL, NULL, &error);
  if (params == NULL)
    g_error ("ncm_serialize_set_property: cannot parse prop_str %s: %s.",
             prop_str, error->message);

  g_assert (g_variant_is_of_type (params, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE)));

  {
    GVariantIter *p_iter = g_variant_iter_new (params);
    GVariant *var = NULL;
    guint i = 0;

    while ((var = g_variant_iter_next_value (p_iter)))
    {
      GVariant *var_key = g_variant_get_child_value (var, 0);
      GVariant *var_val = g_variant_get_child_value (var, 1);
      GVariant *val = g_variant_get_variant (var_val);
      const gchar *name = g_variant_get_string (var_key, NULL);
      GValue value = G_VALUE_INIT;

      if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)))
      {
        GVariant *nest_obj_key = g_variant_get_child_value (val, 0);
        GVariant *nest_obj_params = g_variant_get_child_value (val, 1);
        GObject *nest_obj =
          ncm_serialize_from_name_params (ser,
                                          g_variant_get_string (nest_obj_key, NULL),
                                          nest_obj_params);
        g_value_init (&value, G_TYPE_OBJECT);
        g_value_take_object (&value, nest_obj);
        g_variant_unref (nest_obj_key);
        g_variant_unref (nest_obj_params);
      }
      else
        g_dbus_gvariant_to_gvalue (val, &value);

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
  GError *error = NULL;
  GString *prop_ser = g_string_sized_new (200);

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
    gsize nkeys = 0;
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
    gchar *obj_name = g_match_info_fetch (match_info, 1);
    gchar *obj_prop = g_match_info_fetch (match_info, 2);

    g_free (error_msg);
    if (obj_prop != NULL && strlen (obj_prop) > 0)
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
  gchar *file = NULL;
  gsize length = 0;
  GObject *obj = NULL;

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
    gtype = g_type_from_name (obj_name);

  g_match_info_free (match_info);

  if (name != NULL)
  {
/*
    if (g_hash_table_lookup (ser->saved_name_ser, name) != NULL)
      g_error ("ncm_serialize_from_name_params: deserializing object named `%s' but it is present in the list of saved serializations.",
               name);
    else 
*/
    if (ncm_serialize_contain_name (ser, name))
      obj = ncm_serialize_get_by_name (ser, name);
	}

  if (obj != NULL)
  {
    if (nprop > 0)
      g_error ("ncm_serialize_from_name_params: cannot create object `%s' with properties, it was found on saved/named instances.",
               obj_name);
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
      names[i]          = g_variant_get_string (var_key, NULL);

      if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_VECTOR_TYPE)) && !is_NcmVector)
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
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE)))
      {
        NcmObjArray *oa = ncm_obj_array_new_from_variant (ser, val);
        GValue lval     = G_VALUE_INIT;
        
        g_value_init (&lval, NCM_TYPE_OBJ_ARRAY);
        values[i] = lval;
        g_value_take_boxed (&values[i], oa);
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
        g_dbus_gvariant_to_gvalue (val, &values[i]);

      i++;
      g_variant_unref (var_key);
      g_variant_unref (var_val);
      g_variant_unref (val);
      g_variant_unref (var);
    }

#if GLIB_CHECK_VERSION(2,54,0)
    obj = g_object_new_with_properties (gtype, nprop, names, values);
#else
    {
      GParameter *gprop     = g_new (GParameter, nprop);
      for (i = 0; i < nprop; i++)
      {
        gprop[i].name  = names[i];
        gprop[i].value = values[i];
      }
      obj = g_object_newv (gtype, nprop, gprop);
      g_free (gprop);
    }
#endif /* GLIB_CHECK_VERSION(2,54,0) */

    for (i = 0; i < nprop; i++)
      g_value_unset (&values[i]);

    g_free (names);
    g_free (values);
    g_variant_iter_free (p_iter);
  }
  else
    obj = g_object_new (gtype, NULL);

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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT".", sizeof(gint));
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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT".", sizeof(guint));
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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT".", sizeof(glong));
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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT".", sizeof(gulong));
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
  GType t = G_VALUE_TYPE (val);
  GType fund_t = G_TYPE_FUNDAMENTAL (t);
  const GVariantType *var_type = _ncm_serialize_gtype_to_gvariant_type (fund_t);
  GVariant *var = NULL;

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
        if (g_type_is_a (t, NCM_TYPE_OBJ_ARRAY))
        {
          NcmObjArray *oa = g_value_get_boxed (val);
          if (oa != NULL)
            var = ncm_obj_array_ser (oa, ser);
        }
        else if (g_type_is_a (t, G_TYPE_STRV))
        {
          const gchar * const * strv = g_value_get_boxed (val);
          if (strv != NULL)
            var = g_variant_ref_sink (g_variant_new_strv (strv, -1));
        }
        else
          g_error ("Cannot convert GValue '%s' to GVariant.", g_type_name (t));
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
  else
    var = g_dbus_gvalue_to_gvariant (val, var_type);

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
  gchar *obj_name = g_strdup (G_OBJECT_TYPE_NAME (obj));
  gchar *saved_name = NULL;

  if (ncm_serialize_contain_instance (ser, obj))
  {
    gchar *ni_name = ncm_serialize_peek_name (ser, obj);
    gchar *fname = g_strdup_printf ("%s[%s]", obj_name, ni_name);
    /*printf ("# Found instante %p at ptr_name %s.\n", obj, fname);*/
    ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_TYPE, fname, NULL));
    g_free (fname);
  }
  else if (g_hash_table_lookup_extended (ser->saved_ptr_name, obj, NULL, (gpointer *)&saved_name))
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
    GParamSpec **prop   = g_object_class_list_properties (klass, &n_properties);
    gchar *name = NULL;

    if (ser->opts & NCM_SERIALIZE_OPT_AUTONAME_SER)
    {
      gchar *tmp  = obj_name;
      
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
        GValue val = G_VALUE_INIT;

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
  gchar *serstr = NULL;

  if (valid_variant)
    serstr = g_variant_print (ser_var, TRUE);
  else
  {
    GVariant *params = NULL;
    gchar *obj_name = NULL;
    gchar *params_str;
    g_variant_get (ser_var, NCM_SERIALIZE_OBJECT_FORMAT, &obj_name, &params);

    if (g_variant_n_children (params) == 0)
    {
      serstr = obj_name;
    }
    else
    {
      params_str = g_variant_print (params, TRUE);
      serstr = g_strdup_printf ("%s%s", obj_name, params_str);
      g_free (params_str);
      g_free (obj_name);
    }
    g_variant_unref (params);
  }

  g_variant_unref (ser_var);
  return serstr;
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
    g_error ("ncm_serialize_to_file: cannot save to file %s: %s",
             filename, error->message);

  g_variant_unref (obj_ser);
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

static gsize _global_init = FALSE;
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
  gboolean ret = ncm_serialize_contain_instance (ser, obj);
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
  gboolean ret = ncm_serialize_contain_name (ser, name);
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
  guint ret = ncm_serialize_count_instances (ser);
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
  guint ret = ncm_serialize_count_saved_serializations (ser);
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
  gpointer ret = ncm_serialize_get_by_name (ser, name);
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
  gchar *ret = ncm_serialize_peek_name (ser, obj);
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
  gboolean ret = ncm_serialize_is_named (ser, serobj, name);
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
  GObject *ret = ncm_serialize_from_string (ser, obj_ser);
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
  GObject *ret = ncm_serialize_from_file (ser, filename);
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
  GObject *ret = ncm_serialize_from_binfile (ser, filename);
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
  GObject *ret = ncm_serialize_from_name_params (ser, obj_name, params);
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
  GVariant *ret = ncm_serialize_gvalue_to_gvariant (ser, val);
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
  GVariant *ret = ncm_serialize_to_variant (ser, obj);
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
  gchar *ret = ncm_serialize_to_string (ser, obj, valid_variant);
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
  GVariant *var = ncm_serialize_to_variant (ser, obj);
  GObject *dup = ncm_serialize_from_variant (ser, var);
  g_variant_unref (var);
  ncm_serialize_unref (ser);
  return dup;
}

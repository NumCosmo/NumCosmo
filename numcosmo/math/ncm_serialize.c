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
 * @title: Serialization object.
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
#include "ncm_enum_types.h"
#include <gio/gio.h>

enum
{
  PROP_0,
  PROP_OPTS,
};

G_DEFINE_TYPE (NcmSerialize, ncm_serialize, G_TYPE_OBJECT);

GVariant *_ncm_serialize_to_variant (NcmSerialize *ser, GObject *obj);

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
  
  ser->is_named_regex  = g_regex_new ("^\\s*[A-Z][A-Za-z0-9]+\\s*\\[([A-Za-z0-9\\:]+)\\]\\s*$", 0, 0, &error);
  ser->parse_obj_regex = g_regex_new ("^\\s*([A-Z][A-Za-z0-9]+)\\s*([\\{]?.*[\\}]?)\\s*$", 0, 0, &error);
  ser->autosave_count = 0;
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

/**
 * ncm_serialize_reset:
 * @ser: a #NcmSerialize.
 * 
 * Releases all objects in @ser.
 * 
 */
void 
ncm_serialize_reset (NcmSerialize *ser)
{
  g_hash_table_remove_all (ser->name_ptr);
  g_hash_table_remove_all (ser->ptr_name);
/*
  g_hash_table_remove_all (ser->saved_ptr_name);
  g_hash_table_remove_all (ser->saved_name_ser);
  ser->autosave_count = 0;
*/
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
  return (g_hash_table_lookup (ser->ptr_name, obj) != NULL);
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
ncm_serialize_contain_name (NcmSerialize *ser, gchar *name)
{
  g_assert (name != NULL);
  return (g_hash_table_lookup (ser->name_ptr, name) != NULL);
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
 * ncm_serialize_get_by_name:
 * @ser: a #NcmSerialize.
 * @name: an instance name.
 * 
 * Gets a new reference for the instance @name or null if there isn't a instance named @name.
 * 
 * Returns: (transfer full) (type GObject): Gets the instance named @name or NULL.
 */
gpointer 
ncm_serialize_get_by_name (NcmSerialize *ser, gchar *name)
{
  g_assert (name != NULL);
  {
    GObject *obj = g_hash_table_lookup (ser->name_ptr, name);
    if (obj != NULL)
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
    g_error ("ncm_serialize_peek_name: Cannot peek name of object %p, it is not a named instance.", obj);
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
ncm_serialize_set (NcmSerialize *ser, gpointer obj, gchar *name, gboolean overwrite)
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

static void 
_ncm_serialize_save_ser (NcmSerialize *ser, gchar *name, gpointer obj, GVariant *ser_var)
{
  g_assert (G_IS_OBJECT (obj));
  g_assert (name != NULL);
  g_assert (ser != NULL);

  if (g_hash_table_lookup (ser->saved_name_ser, name) != NULL)
    g_error ("_ncm_serialize_save_ser: instance already saved.");
  if (g_hash_table_lookup (ser->saved_ptr_name, obj) != NULL)
    g_error ("_ncm_serialize_save_ser: instance already saved.");

  /* printf ("\n``%s'' <=> %s\n", name, g_variant_print (ser_var, TRUE)); */
  
  g_hash_table_insert (ser->saved_ptr_name, 
                       g_object_ref (obj), g_strdup (name));
  g_hash_table_insert (ser->saved_name_ser, 
                       g_strdup (name), g_variant_ref_sink (ser_var));
  
  ser->autosave_count++;
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
  GError *error = NULL;
  GMatchInfo *match_info = NULL;
  GVariant *var_obj = NULL;
  GVariant *obj_name = NULL;
  GVariant *params = NULL;
  
  var_obj = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), serobj, NULL, NULL, &error);
  g_clear_error (&error);

  if (var_obj != NULL)
  {
    obj_name = g_variant_get_child_value (var_obj, 0);
    params = g_variant_get_child_value (var_obj, 1);
    serobj = g_variant_get_string (obj_name, NULL);
  }

  if (g_regex_match (ser->is_named_regex, serobj, 0, &match_info))
  {
    if (name != NULL)
      *name = g_match_info_fetch (match_info, 1);
    g_match_info_free (match_info);

    if (params != NULL && g_variant_n_children (params) > 0)
      g_error ("ncm_serialize_is_named: Invalid named object ``%s'', named objects cannot have properties ``%s''.", 
               serobj, g_variant_print (params, TRUE));

    if (obj_name != NULL)
      g_variant_unref (obj_name);
    if (var_obj != NULL)
      g_variant_unref (var_obj);
    if (params != NULL)
      g_variant_unref (params);
    return TRUE;
  }
  else
  {
    if (obj_name != NULL)
      g_variant_unref (obj_name);
    if (var_obj != NULL)
      g_variant_unref (var_obj);
    if (params != NULL)
      g_variant_unref (params);
    if (name != NULL)
      *name = NULL;
    return FALSE;
  }
}

/* Stole from glib 2.36.1 */
/* GDBus - GLib D-Bus Library
 *
 * Copyright (C) 2008-2010 Red Hat, Inc.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * Author: David Zeuthen <davidz@redhat.com>
 */
#if ((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 30))

static void
g_dbus_gvariant_to_gvalue (GVariant *value, GValue *out_gvalue)
{
  const GVariantType *type;
  gchar **array;

  g_return_if_fail (value != NULL);
  g_return_if_fail (out_gvalue != NULL);

  memset (out_gvalue, '\0', sizeof (GValue));

  switch (g_variant_classify (value))
  {
    case G_VARIANT_CLASS_BOOLEAN:
      g_value_init (out_gvalue, G_TYPE_BOOLEAN);
      g_value_set_boolean (out_gvalue, g_variant_get_boolean (value));
      break;

    case G_VARIANT_CLASS_BYTE:
      g_value_init (out_gvalue, G_TYPE_UCHAR);
      g_value_set_uchar (out_gvalue, g_variant_get_byte (value));
      break;

    case G_VARIANT_CLASS_INT16:
      g_value_init (out_gvalue, G_TYPE_INT);
      g_value_set_int (out_gvalue, g_variant_get_int16 (value));
      break;

    case G_VARIANT_CLASS_UINT16:
      g_value_init (out_gvalue, G_TYPE_UINT);
      g_value_set_uint (out_gvalue, g_variant_get_uint16 (value));
      break;

    case G_VARIANT_CLASS_INT32:
      g_value_init (out_gvalue, G_TYPE_INT);
      g_value_set_int (out_gvalue, g_variant_get_int32 (value));
      break;

    case G_VARIANT_CLASS_UINT32:
      g_value_init (out_gvalue, G_TYPE_UINT);
      g_value_set_uint (out_gvalue, g_variant_get_uint32 (value));
      break;

    case G_VARIANT_CLASS_INT64:
      g_value_init (out_gvalue, G_TYPE_INT64);
      g_value_set_int64 (out_gvalue, g_variant_get_int64 (value));
      break;

    case G_VARIANT_CLASS_UINT64:
      g_value_init (out_gvalue, G_TYPE_UINT64);
      g_value_set_uint64 (out_gvalue, g_variant_get_uint64 (value));
      break;

    case G_VARIANT_CLASS_DOUBLE:
      g_value_init (out_gvalue, G_TYPE_DOUBLE);
      g_value_set_double (out_gvalue, g_variant_get_double (value));
      break;

    case G_VARIANT_CLASS_STRING:
      g_value_init (out_gvalue, G_TYPE_STRING);
      g_value_set_string (out_gvalue, g_variant_get_string (value, NULL));
      break;

    case G_VARIANT_CLASS_OBJECT_PATH:
      g_value_init (out_gvalue, G_TYPE_STRING);
      g_value_set_string (out_gvalue, g_variant_get_string (value, NULL));
      break;

    case G_VARIANT_CLASS_SIGNATURE:
      g_value_init (out_gvalue, G_TYPE_STRING);
      g_value_set_string (out_gvalue, g_variant_get_string (value, NULL));
      break;
      
    case G_VARIANT_CLASS_ARRAY:
    {
      type = g_variant_get_type (value);
      switch (g_variant_type_peek_string (type)[1])
      {
        case G_VARIANT_CLASS_BYTE:
          g_value_init (out_gvalue, G_TYPE_STRING);
          g_value_set_string (out_gvalue, g_variant_get_bytestring (value));
          break;
          
        case G_VARIANT_CLASS_STRING:
          g_value_init (out_gvalue, G_TYPE_STRV);
          array = g_variant_dup_strv (value, NULL);
          g_value_take_boxed (out_gvalue, array);
          break;
          
        case G_VARIANT_CLASS_ARRAY:
        {
          switch (g_variant_type_peek_string (type)[2])
          {
            case G_VARIANT_CLASS_BYTE:
              g_value_init (out_gvalue, G_TYPE_STRV);
              array = g_variant_dup_bytestring_array (value, NULL);
              g_value_take_boxed (out_gvalue, array);
              break;
              
            default:
              g_value_init (out_gvalue, G_TYPE_VARIANT);
              g_value_set_variant (out_gvalue, value);
              break;
          }
          break;
        }
        default:
          g_value_init (out_gvalue, G_TYPE_VARIANT);
          g_value_set_variant (out_gvalue, value);
          break;
      }
      break;
    }
    case G_VARIANT_CLASS_HANDLE:
    case G_VARIANT_CLASS_VARIANT:
    case G_VARIANT_CLASS_MAYBE:
    case G_VARIANT_CLASS_TUPLE:
    case G_VARIANT_CLASS_DICT_ENTRY:
      g_value_init (out_gvalue, G_TYPE_VARIANT);
      g_value_set_variant (out_gvalue, value);
      break;
  }
}

static GVariant *
g_dbus_gvalue_to_gvariant (const GValue *gvalue, const GVariantType *type)
{
  GVariant *ret;
  const gchar *s;
  const gchar * const *as;
  const gchar *empty_strv[1] = {NULL};

  g_return_val_if_fail (gvalue != NULL, NULL);
  g_return_val_if_fail (type != NULL, NULL);

  ret = NULL;

  if (G_VALUE_TYPE (gvalue) == G_TYPE_VARIANT)
  {
    ret = g_value_dup_variant (gvalue);
  }
  else
  {
    switch (g_variant_type_peek_string (type)[0])
    {
      case G_VARIANT_CLASS_BOOLEAN:
        ret = g_variant_ref_sink (g_variant_new_boolean (g_value_get_boolean (gvalue)));
        break;

      case G_VARIANT_CLASS_BYTE:
        ret = g_variant_ref_sink (g_variant_new_byte (g_value_get_uchar (gvalue)));
        break;

      case G_VARIANT_CLASS_INT16:
        ret = g_variant_ref_sink (g_variant_new_int16 (g_value_get_int (gvalue)));
        break;

      case G_VARIANT_CLASS_UINT16:
        ret = g_variant_ref_sink (g_variant_new_uint16 (g_value_get_uint (gvalue)));
        break;

      case G_VARIANT_CLASS_INT32:
        ret = g_variant_ref_sink (g_variant_new_int32 (g_value_get_int (gvalue)));
        break;

      case G_VARIANT_CLASS_UINT32:
        ret = g_variant_ref_sink (g_variant_new_uint32 (g_value_get_uint (gvalue)));
        break;

      case G_VARIANT_CLASS_INT64:
        ret = g_variant_ref_sink (g_variant_new_int64 (g_value_get_int64 (gvalue)));
        break;

      case G_VARIANT_CLASS_UINT64:
        ret = g_variant_ref_sink (g_variant_new_uint64 (g_value_get_uint64 (gvalue)));
        break;

      case G_VARIANT_CLASS_DOUBLE:
        ret = g_variant_ref_sink (g_variant_new_double (g_value_get_double (gvalue)));
        break;

      case G_VARIANT_CLASS_STRING:
        s = g_value_get_string (gvalue);
        if (s == NULL)
          s = "";
        ret = g_variant_ref_sink (g_variant_new_string (s));
        break;

      case G_VARIANT_CLASS_OBJECT_PATH:
        s = g_value_get_string (gvalue);
        if (s == NULL)
          s = "/";
        ret = g_variant_ref_sink (g_variant_new_object_path (s));
        break;

      case G_VARIANT_CLASS_SIGNATURE:
        s = g_value_get_string (gvalue);
        if (s == NULL)
          s = "";
        ret = g_variant_ref_sink (g_variant_new_signature (s));
        break;

      case G_VARIANT_CLASS_ARRAY:
      {
        switch (g_variant_type_peek_string (type)[1])
        {
          case G_VARIANT_CLASS_BYTE:
          {
            s = g_value_get_string (gvalue);
            if (s == NULL)
              s = "";
            ret = g_variant_ref_sink (g_variant_new_bytestring (s));
            break;
          }
          case G_VARIANT_CLASS_STRING:
          {
            as = g_value_get_boxed (gvalue);
            if (as == NULL)
              as = empty_strv;
            ret = g_variant_ref_sink (g_variant_new_strv (as, -1));
            break;
          }
          case G_VARIANT_CLASS_ARRAY:
          {
            switch (g_variant_type_peek_string (type)[2])
            {
              case G_VARIANT_CLASS_BYTE:
                as = g_value_get_boxed (gvalue);
                if (as == NULL)
                  as = empty_strv;
                ret = g_variant_ref_sink (g_variant_new_bytestring_array (as, -1));
                break;
              default:
                ret = g_value_dup_variant (gvalue);
                break;
            }
            break;
          }
          default:
            ret = g_value_dup_variant (gvalue);
            break;
        }
        break;
      }
      case G_VARIANT_CLASS_HANDLE:
      case G_VARIANT_CLASS_VARIANT:
      case G_VARIANT_CLASS_MAYBE:
      case G_VARIANT_CLASS_TUPLE:
      case G_VARIANT_CLASS_DICT_ENTRY:
        ret = g_value_dup_variant (gvalue);
        break;
    }
  }

  if (ret == NULL)
  {
    GVariant *untrusted_empty;
    untrusted_empty = g_variant_new_from_data (type, NULL, 0, FALSE, NULL, NULL);
    ret = g_variant_ref_sink (g_variant_get_normal_form (untrusted_empty));
    g_variant_unref (untrusted_empty);
  }

  g_assert (!g_variant_is_floating (ret));

  return ret;
}

#endif
/* Stole from glib 2.36.1 */

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
    }

    g_variant_iter_free (p_iter);
  }
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
    GVariant *params = g_variant_get_child_value (var_obj, 1);
    GObject *obj = ncm_serialize_from_name_params (ser, g_variant_get_string (obj_name, NULL), params);

    g_variant_unref (obj_name);
    g_variant_unref (params);
    return obj;
  }
}

static GObject *
_ncm_serialize_from_saved_or_named (NcmSerialize *ser, const gchar *obj_ser, GVariant *params)
{
  GObject *obj = NULL;
  gchar *ni_name = NULL;

  if (ncm_serialize_is_named (ser, obj_ser, &ni_name))
  {
    GVariant *var_obj;
    if (params != NULL && g_variant_n_children (params) > 0)
      g_error ("_ncm_serialize_from_saved_or_named: Invalid named object parameters ``%s'', named objects cannot have properties.", 
               g_variant_print (params, TRUE));
    
    if (ncm_serialize_contain_name (ser, ni_name))
    {
      obj = ncm_serialize_get_by_name (ser, ni_name);
      /* printf ("found ni: ``%s'' <=> %s\n", ni_name, obj_ser); */
      g_free (ni_name);
      return obj;
    }
    else if ((var_obj = g_hash_table_lookup (ser->saved_name_ser, ni_name)) != NULL)
    {
      obj = ncm_serialize_from_variant (ser, var_obj);
      if (ser->opts & NCM_SERIALIZE_OPT_AUTONAME_SER)
        ncm_serialize_set (ser, obj, ni_name, FALSE);
      /* printf ("found saved: ``%s'' <=> %s\n", ni_name, obj_ser); */
      g_free (ni_name);
      return obj;
    }
    else
      g_error ("_ncm_serialize_from_saved_or_named: Creating object from ``%s'', depending object ``%s'' not found.", 
               obj_ser, ni_name);    
  }
  
  return obj;
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
  GError *error = NULL;
  GMatchInfo *match_info = NULL;
  GVariant *var_obj = NULL;
  GObject *obj = _ncm_serialize_from_saved_or_named (ser, obj_ser, NULL);

  if (obj != NULL)
    return obj;

  var_obj = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), obj_ser, NULL, NULL, &error);
  g_clear_error (&error);
  if (var_obj != NULL)
  {
    obj = ncm_serialize_from_variant (ser, var_obj);
    g_variant_unref (var_obj);
  }
  else if (g_regex_match (ser->parse_obj_regex, obj_ser, 0, &match_info))
  {
    gchar *obj_name = g_match_info_fetch (match_info, 1);
    gchar *obj_prop = g_match_info_fetch (match_info, 2);
    GVariant *params = NULL;

    if (obj_prop != NULL && strlen (obj_prop) > 0)
    {
      params = g_variant_parse (G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE), 
                                obj_prop, NULL, NULL, &error);
      if (error != NULL)
      {
        g_error ("ncm_serialize_from_string: cannot parse object '%s' parameters '%s', error %s\n", obj_name, obj_prop, error->message);
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
    g_error ("ncm_serialize_from_string: cannot indentify object (%s) in string '%s'\n", NCM_SERIALIZE_OBJECT_TYPE, obj_ser);
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
  GType gtype = g_type_from_name (obj_name);
  GObject *obj = _ncm_serialize_from_saved_or_named (ser, obj_name, params);

  if (obj != NULL)
    return obj;

  if (gtype == 0)
    g_error ("ncm_serialize_from_name_params: object '%s' is not registered\n", obj_name);

  g_assert (params == NULL || g_variant_is_of_type (params, G_VARIANT_TYPE (NCM_SERIALIZE_PROPERTIES_TYPE)));

  if (params != NULL)
  {
    GVariantIter *p_iter = g_variant_iter_new (params);
    gsize nprop = g_variant_iter_n_children (p_iter);
    GParameter *gprop = g_new (GParameter, nprop);
    GVariant *var = NULL;
    guint i = 0;

    while ((var = g_variant_iter_next_value (p_iter)))
    {
      GVariant *var_key = g_variant_get_child_value (var, 0);
      GVariant *var_val = g_variant_get_child_value (var, 1);
      GVariant *val = g_variant_get_variant (var_val);
      gprop[i].name = g_variant_get_string (var_key, NULL);
      /*
      if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_VECTOR_TYPE)))
      {
        NcmVector *vec = ncm_vector_new_variant (var);
        GValue lval = G_VALUE_INIT;
        g_value_init (&lval, G_TYPE_OBJECT);
        gprop[i].value = lval;
        g_value_take_object (&gprop[i].value, vec);
        ncm_vector_free (vec);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_MATRIX_TYPE)))
      {
        NcmMatrix *mat = ncm_matrix_new_variant (var);
        GValue lval = G_VALUE_INIT;
        g_value_init (&lval, G_TYPE_OBJECT);
        gprop[i].value = lval;
        g_value_take_object (&gprop[i].value, mat);
        ncm_matrix_free (mat);        
      }
      else */
      if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE)))
      {
        NcmObjArray *oa = ncm_obj_array_new_from_variant (ser, val);
        GValue lval = G_VALUE_INIT;
        g_value_init (&lval, NCM_TYPE_OBJ_ARRAY);
        gprop[i].value = lval;
        g_value_take_boxed (&gprop[i].value, oa);
      }
      else if (g_variant_is_of_type (val, G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE)))
      {
        GVariant *nest_obj_key = g_variant_get_child_value (val, 0);
        GVariant *nest_obj_params = g_variant_get_child_value (val, 1);
        GValue lval = G_VALUE_INIT;
        GObject *nest_obj =
          ncm_serialize_from_name_params (ser,
                                          g_variant_get_string (nest_obj_key, NULL),
                                          nest_obj_params);
        g_value_init (&lval, G_TYPE_OBJECT);
        gprop[i].value = lval;
        g_value_take_object (&gprop[i].value, nest_obj);
        g_variant_unref (nest_obj_key);
        g_variant_unref (nest_obj_params);
      }
      else
        g_dbus_gvariant_to_gvalue (val, &gprop[i].value);

      i++;
      g_variant_unref (var_key);
      g_variant_unref (var_val);
      g_variant_unref (val);
      g_variant_unref (var);
    }

    obj = g_object_newv (gtype, nprop, gprop);
    for (i = 0; i < nprop; i++)
      g_value_unset (&gprop[i].value);

    g_free (gprop);
    g_variant_iter_free (p_iter);
  }
  else
    obj = g_object_new (gtype, NULL);

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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(gint));
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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(guint));
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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(glong));
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
          g_error ("Unknown gint size %"G_GSIZE_FORMAT"\n", sizeof(gulong));
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

GVariant *_ncm_obj_array_ser (NcmObjArray *oa, NcmSerialize *ser);


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
          var = _ncm_serialize_to_variant (ser, nest_obj);
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
            var = _ncm_obj_array_ser (oa, ser);
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
  else
    var = g_dbus_gvalue_to_gvariant (val, var_type);

  return var;
}

GVariant *
_ncm_serialize_to_variant (NcmSerialize *ser, GObject *obj)
{
  GVariant *ser_var;
  const gchar *obj_name = G_OBJECT_TYPE_NAME (obj);
  gchar *saved_name = NULL;
  
  if (ncm_serialize_contain_instance (ser, obj))
  {
    gchar *ni_name = ncm_serialize_peek_name (ser, obj);
    gchar *fname = g_strdup_printf ("%s[%s]", obj_name, ni_name);

    ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_TYPE, fname, NULL));
    g_free (fname);
  }
  else if ((saved_name = g_hash_table_lookup (ser->saved_ptr_name, obj)) != NULL)
  {
    gchar *fname = g_strdup_printf ("%s[%s]", obj_name, saved_name);
    ser_var = g_variant_ref_sink (g_variant_new (NCM_SERIALIZE_OBJECT_TYPE, fname, NULL));
    g_free (fname);
  }
  else
  {
    GObjectClass *klass = G_OBJECT_GET_CLASS (obj);
    guint n_properties, i;
    GParamSpec **prop = g_object_class_list_properties (klass, &n_properties);

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

    if (ser->opts & NCM_SERIALIZE_OPT_AUTOSAVE_SER)
    {
      gchar *name = g_strdup_printf (NCM_SERIALIZE_AUTOSAVE_NAME NCM_SERIALIZE_AUTOSAVE_NFORMAT, ser->autosave_count);
      _ncm_serialize_save_ser (ser, name, obj, ser_var);      
      g_free (name);
      g_variant_unref (ser_var);
      ser_var = _ncm_serialize_to_variant (ser, obj);
    }
  }
  return ser_var;
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
  if (ser->opts & NCM_SERIALIZE_OPT_AUTOSAVE_SER)
  {
    GVariant *cser_var = _ncm_serialize_to_variant (ser, obj);
    gchar *name = g_hash_table_lookup (ser->saved_ptr_name, obj);
    GVariant *ser_var = g_hash_table_lookup (ser->saved_name_ser, name);
  
    g_assert (name != NULL);
    g_assert (ser_var != NULL);

    g_variant_unref (cser_var);
    
    return g_variant_ref_sink (ser_var);
  }
  else
    return _ncm_serialize_to_variant (ser, obj);
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
    }
    g_free (obj_name);
    g_variant_unref (params);
  }
  g_variant_unref (ser_var);
  return serstr;
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
  GObject *dup = ncm_serialize_from_variant (ser, var);
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
 * 
 * Releases all objects in global #NcmSerialize.
 * 
 */
void 
ncm_serialize_global_reset (void)
{
  NcmSerialize *ser = ncm_serialize_global ();
  ncm_serialize_reset (ser);
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
ncm_serialize_global_contain_name (gchar *name)
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
 * ncm_serialize_global_get_by_name:
 * @name: an instance name.
 * 
 * Global version of ncm_serialize_get_by_name().
 * 
 * Returns: (transfer full) (type GObject): Gets the instance named @name or NULL.
 */
gpointer 
ncm_serialize_global_get_by_name (gchar *name)
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
ncm_serialize_global_set (gpointer obj, gchar *name, gboolean overwrite)
{
  NcmSerialize *ser = ncm_serialize_global ();
  ncm_serialize_set (ser, obj, name, overwrite);
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

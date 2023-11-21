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
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_obj_array.h"
#include "math/ncm_cfg.h"

G_DEFINE_BOXED_TYPE (NcmObjArray, ncm_obj_array, ncm_obj_array_ref, ncm_obj_array_unref)

/**
 * ncm_obj_array_new:
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_new ()
{
  GPtrArray *oa = g_ptr_array_new ();
  g_ptr_array_set_free_func (oa, g_object_unref);
  return (NcmObjArray *) oa;
}

/**
 * ncm_obj_array_new_from_variant:
 * @ser: a #NcmSerialize.
 * @var: a #GVariant containing an array of objects.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_new_from_variant (NcmSerialize *ser, GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE)));
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
 * ncm_obj_array_sized_new:
 * @n: initial allocation size.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
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
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_ref (NcmObjArray *oa)
{
  return (NcmObjArray *)g_ptr_array_ref ((GPtrArray *)oa);
}

/**
 * ncm_obj_array_unref:
 * @oa: a #NcmObjArray.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_unref (NcmObjArray *oa)
{
  g_ptr_array_unref ((GPtrArray *)oa);
}

/**
 * ncm_obj_array_clear:
 * @oa: a pointer to a #NcmObjArray.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_clear (NcmObjArray **oa)
{
  g_clear_pointer (oa, ncm_obj_array_unref);
}

/**
 * ncm_obj_array_ser:
 * @oa: a #NcmObjArray.
 * @ser: a #NcmSerialize.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
GVariant *
ncm_obj_array_ser (NcmObjArray *oa, NcmSerialize *ser)
{
  GVariantBuilder *builder;
  GVariant *var;
  guint i;

  builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE));

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
 * ncm_obj_array_dup:
 * @oa: a #NcmObjArray.
 * @ser: a #NcmSerialize.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_dup (NcmObjArray *oa, NcmSerialize *ser)
{
  GVariant *var    = ncm_obj_array_ser (oa, ser);
  NcmObjArray *dup = ncm_obj_array_new_from_variant (ser, var);
  g_variant_unref (var);
  return dup;
}

/**
 * ncm_obj_array_add:
 * @oa: a #NcmObjArray.
 * @obj: a #GObject.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_add (NcmObjArray *oa, GObject *obj)
{
  g_assert (obj != NULL);
  g_ptr_array_add ((GPtrArray *)oa, g_object_ref (obj));
}

/**
 * ncm_obj_array_set:
 * @oa: a #NcmObjArray.
 * @i: object index.
 * @obj: a #GObject.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_set (NcmObjArray *oa, guint i, GObject *obj)
{
  g_assert_cmpuint (i, <, oa->len);

  g_assert (obj != NULL);
  if (obj != g_ptr_array_index ((GPtrArray *)oa, i))
	{
		g_object_unref (g_ptr_array_index ((GPtrArray *)oa, i));
		g_ptr_array_index ((GPtrArray *)oa, i) = g_object_ref (obj);
	}
}

/**
 * ncm_obj_array_get:
 * @oa: a #NcmObjArray.
 * @i: object index.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
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
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
GObject *
ncm_obj_array_peek (NcmObjArray *oa, guint i)
{
  g_assert_cmpuint (i, <, oa->len);
  return g_ptr_array_index ((GPtrArray *)oa, i);
}

/**
 * ncm_obj_array_len:
 * @oa: a #NcmObjArray
 *
 * FIXME
 *
 * Returns: array length
 */
gint
ncm_obj_array_len (NcmObjArray *oa)
{
  return oa->len;
}

/**
 * ncm_obj_array_save:
 * @oa: a #NcmObjArray
 * @ser: a #NcmSerialize
 * @filename: FIXME
 * @save_comment: FIXME
 *
 * FIXME
 *
 */
void
ncm_obj_array_save (NcmObjArray *oa, NcmSerialize *ser, const gchar *filename, gboolean save_comment)
{
  GKeyFile *oafile = g_key_file_new ();
  guint i;

  {
    GError *error = NULL;
    gchar *oa_desc = ncm_cfg_string_to_comment ("Whether NcmObjArray is empty");

    g_key_file_set_boolean (oafile, "NcmObjArray", "empty", oa->len == 0 ? TRUE : FALSE);
    if (save_comment)
    {
      if (!g_key_file_set_comment (oafile, "NcmObjArray", NULL, oa_desc, &error))
        g_error ("ncm_obj_array_save: %s", error->message);
    }

    g_free (oa_desc);
  }

  for (i = 0; i < oa->len; i++)
  {
    GObject *go          = ncm_obj_array_peek (oa, i);
    GObjectClass *oclass = G_OBJECT_GET_CLASS (go);
    GError *error        = NULL;
    gchar *group         = g_strdup_printf (NCM_OBJ_ARRAY_POS_STR":%d", i);
    GVariant *go_var     = ncm_serialize_to_variant (ser, go);
    GVariant *params     = NULL;
    gchar *obj_name      = NULL;
    guint nparams;

    g_variant_get (go_var, "{s@a{sv}}", &obj_name, &params);
    nparams = g_variant_n_children (params);
    
    g_key_file_set_value (oafile, group, NCM_OBJ_ARRAY_OBJ_NAME_STR, obj_name);
    
    if (nparams != 0)
    {
      GVariantIter iter;
      GVariant *value;
      gchar *key;
      g_variant_iter_init (&iter, params);
      while (g_variant_iter_next (&iter, "{sv}", &key, &value))
      {
        GParamSpec *param_spec = g_object_class_find_property (oclass, key);
        gchar *param_str = g_variant_print (value, TRUE);
        
        if (param_spec == NULL)
          g_error ("ncm_obj_array_save: property `%s' not found in object `%s'.", key, obj_name);

        g_key_file_set_value (oafile, group, key, param_str);
        
        if (save_comment)
        {
          const gchar *blurb = g_param_spec_get_blurb (param_spec);
          if (blurb != NULL && blurb[0] != 0)
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
    GError *error = NULL;
    gsize len = 0;
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
 * ncm_obj_array_load:
 * @filename: oa filename
 * @ser: a #NcmSerialize
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_load (const gchar *filename, NcmSerialize *ser)
{
  NcmObjArray *oa = ncm_obj_array_new ();
  GKeyFile *oafile = g_key_file_new ();
  GError *error = NULL;
  gchar **groups = NULL;
  gsize ngroups = 0;
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

    if (!g_key_file_has_key (oafile, groups[i], NCM_OBJ_ARRAY_OBJ_NAME_STR, &error))
    {
      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);
      g_error ("ncm_obj_array_load: Every group must contain the key `%s' containing the object name.", NCM_OBJ_ARRAY_OBJ_NAME_STR);
    }
    else
    {
      obj_name = g_key_file_get_value (oafile, groups[i], NCM_OBJ_ARRAY_OBJ_NAME_STR, &error);
      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);      
      g_key_file_remove_key (oafile, groups[i], NCM_OBJ_ARRAY_OBJ_NAME_STR, &error);
      if (error != NULL)
        g_error ("ncm_obj_array_load: %s", error->message);
    }

    g_string_append_printf (obj_ser, "{\'%s\', @a{sv} {", obj_name);    

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

      g_string_append (obj_ser, "}}");
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

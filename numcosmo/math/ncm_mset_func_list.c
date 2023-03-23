/***************************************************************************
 *            ncm_mset_func_list.c
 *
 *  Mon August 08 17:29:34 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_func_list.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_mset_func_list
 * @title: NcmMSetFuncList
 * @short_description: NcmMSet Functions list.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_func_list.h"

enum
{
  PROP_0,
  PROP_FULL_NAME,
  PROP_OBJECT,
};

G_DEFINE_TYPE (NcmMSetFuncList, ncm_mset_func_list, NCM_TYPE_MSET_FUNC);

static void
ncm_mset_func_list_init (NcmMSetFuncList *flist)
{
  flist->obj_type = G_TYPE_NONE;
  flist->func     = NULL;
  flist->obj      = NULL;
}

static void _ncm_mset_func_list_init_from_full_name (NcmMSetFuncList *flist, const gchar *full_name);

static void
_ncm_mset_func_list_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetFuncList *flist = NCM_MSET_FUNC_LIST (object);
  NcmMSetFunc *func = NCM_MSET_FUNC (object);
  g_return_if_fail (NCM_IS_MSET_FUNC_LIST (object));

  switch (prop_id)
  {
    case PROP_FULL_NAME:
      _ncm_mset_func_list_init_from_full_name (flist, g_value_get_string (value));
      break;
    case PROP_OBJECT:
      flist->obj = g_value_dup_object (value);
      if (flist->obj != NULL)
      {
        g_assert (g_type_is_a (G_OBJECT_TYPE (flist->obj), flist->obj_type));
      }
      else if (flist->obj_type != G_TYPE_NONE)
      {
        g_error ("_ncm_mset_func_list_set_property: object %s:%s requires an object `%s'.",
                 func->ns,
                 func->name,
                 g_type_name (flist->obj_type));
      }
      
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_func_list_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetFuncList *flist = NCM_MSET_FUNC_LIST (object);
  NcmMSetFunc *func = NCM_MSET_FUNC (object);
  g_return_if_fail (NCM_IS_MSET_FUNC_LIST (object));

  switch (prop_id)
  {
    case PROP_FULL_NAME:
    {
      gchar *ns_name = g_strdup_printf ("%s:%s", func->ns, func->name);
      g_value_take_string (value, ns_name);
      break;
    }
    case PROP_OBJECT:
      g_value_set_object (value, flist->obj);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_func_list_dispose (GObject *object)
{
  NcmMSetFuncList *flist = NCM_MSET_FUNC_LIST (object);

  g_clear_object (&flist->obj);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_list_parent_class)->dispose (object);
}

static void
_ncm_mset_func_list_finalize (GObject *object)
{
  /*NcmMSetFuncList *flist = NCM_MSET_FUNC_LIST (object);*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_list_parent_class)->finalize (object);
}

static void _ncm_mset_func_list_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_mset_func_list_class_init (NcmMSetFuncListClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_ncm_mset_func_list_set_property;
  object_class->get_property = &_ncm_mset_func_list_get_property;
  object_class->dispose      = &_ncm_mset_func_list_dispose;
  object_class->finalize     = &_ncm_mset_func_list_finalize;

  g_object_class_install_property (object_class,
                                   PROP_FULL_NAME,
                                   g_param_spec_string ("full-name",
                                                        NULL,
                                                        "Namespace and function name",
                                                        "NULL",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_OBJECT,
                                   g_param_spec_object ("object",
                                                        NULL,
                                                        "object",
                                                        G_TYPE_OBJECT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->func_array       = g_array_new (TRUE, TRUE, sizeof (NcmMSetFuncListStruct));
  klass->ns_hash          = g_hash_table_new (g_str_hash, g_str_equal);
  func_class->eval        = &_ncm_mset_func_list_eval;
}

G_LOCK_DEFINE_STATIC (insert_lock);

static void 
_ncm_mset_func_list_init_from_full_name (NcmMSetFuncList *flist, const gchar *full_name)
{
  G_LOCK (insert_lock);
  {
    NcmMSetFunc *func                 = NCM_MSET_FUNC (flist);
    NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);

    gchar **ns_name       = g_strsplit (full_name, ":", 2);
    GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, ns_name[0]);

    if (g_strv_length (ns_name) != 2)
      g_error ("_ncm_mset_func_list_init_from_full_name: invalid full_name `%s'.", full_name);
    
    if (func_hash == NULL)
    {
      g_error ("_ncm_mset_func_list_init_from_full_name: namespace `%s' not found.", ns_name[0]);
    }
    else
    {
      gpointer fdata_i;
      if (g_hash_table_lookup_extended (func_hash, ns_name[1], NULL, &fdata_i))
      {
        NcmMSetFuncListStruct *fdata = &g_array_index (flist_class->func_array, NcmMSetFuncListStruct, GPOINTER_TO_INT (fdata_i));

        g_clear_pointer (&func->name, g_free);
        g_clear_pointer (&func->symbol, g_free);
        g_clear_pointer (&func->ns, g_free);
        g_clear_pointer (&func->desc, g_free);

        func->name      = g_strdup (fdata->name);
        func->symbol    = g_strdup (fdata->symbol);
        func->ns        = g_strdup (fdata->ns);
        func->desc      = g_strdup (fdata->desc);
        flist->obj_type = fdata->obj_type;
        flist->func     = fdata->func;

        func->nvar     = fdata->nvar;
        func->dim      = fdata->dim;
      }
      else
        g_error ("_ncm_mset_func_list_init_from_full_name: name `%s' not found in namespace `%s'.", ns_name[1], ns_name[0]);
    }

    g_strfreev (ns_name);

    g_type_class_unref (flist_class);
  }
  G_UNLOCK (insert_lock);
}

static void 
_ncm_mset_func_list_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmMSetFuncList *flist = NCM_MSET_FUNC_LIST (func);

  if (flist->obj_type != G_TYPE_NONE && flist->obj == NULL)
    g_error ("_ncm_mset_func_list_eval: calling function without object `%s'.", g_type_name (flist->obj_type));
  
  flist->func (flist, mset, x, res);
}

/**
 * ncm_mset_func_list_register:
 * @name: function name
 * @symbol: function symbol
 * @ns: namespace
 * @desc: function description
 * @obj_type: object type
 * @func: (scope notified): function pointer
 * @nvar: number of variables
 * @dim: function dimension
 * 
 * Register a new function in the NcmMSetFuncList class.
 * 
 */
void
ncm_mset_func_list_register (const gchar *name, const gchar *symbol, const gchar *ns, const gchar *desc, GType obj_type, NcmMSetFuncListN func, guint nvar, guint dim)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  NcmMSetFuncListStruct flist_item = {
    g_strdup (name), 
    g_strdup (symbol), 
    g_strdup (ns),
    g_strdup (desc), 
    obj_type,
    func,
    nvar,
    dim,
    0
  };

  G_LOCK (insert_lock);

  flist_item.pos = flist_class->func_array->len;
  
  g_array_append_val (flist_class->func_array, flist_item);

  {
    GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, ns);
    if (func_hash == NULL)
    {
      func_hash = g_hash_table_new (g_str_hash, g_str_equal);
      g_hash_table_insert (flist_class->ns_hash, flist_item.ns, func_hash);
    }
    
    g_hash_table_insert (func_hash, flist_item.name, GINT_TO_POINTER (flist_item.pos));
  }

  G_UNLOCK (insert_lock);
  
  g_type_class_unref (flist_class);
}

/**
 * ncm_mset_func_list_select:
 * @ns: (allow-none): namespace
 * @nvar: number of variables
 * @dim: function dimension
 * 
 * Selects the NcmMSetFuncListStruct array containing the function
 * in the namespace @ns with @nvar and @dim. If @ns is NULL then 
 * gets from all namespaces, @nvar and/or @dim equals to -1 selects
 * any value. The contained strings must not be freed.
 * 
 * Returns: (transfer container) (element-type NcmMSetFuncListStruct): NcmMSetFuncListStruct array.
 */
GArray *
ncm_mset_func_list_select (const gchar *ns, gint nvar, gint dim)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  GArray *s = g_array_new (TRUE, TRUE, sizeof (NcmMSetFuncListStruct));
  G_LOCK (insert_lock);

  {
    gint i;
    for (i = 0; i < flist_class->func_array->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (flist_class->func_array, NcmMSetFuncListStruct, i);
      gboolean in_ns   = (ns == NULL) || g_str_has_prefix (fdata->ns, ns);
      gboolean in_nvar = (nvar == -1) || nvar == fdata->nvar;
      gboolean in_dim  = (dim  == -1) || dim  == fdata->dim;
      
      if (in_ns && in_nvar && in_dim)
      {
        g_array_append_val (s, fdata[0]);
      }
    }
  }

  G_UNLOCK (insert_lock);
  
  g_type_class_unref (flist_class);

  return s;
}

/**
 * ncm_mset_func_list_new:
 * @full_name: function full name
 * @obj: (allow-none): FIXME
 *
 * FIXME 
 * 
 * Returns: (transfer full): newly created #NcmMSetFuncList.
 */
NcmMSetFuncList *
ncm_mset_func_list_new (const gchar *full_name, GObject *obj)
{
  NcmMSetFuncList *flist = g_object_new (NCM_TYPE_MSET_FUNC_LIST,
                                         "full-name", full_name,
                                         "object", obj,
                                         NULL);
  return flist;
}

/**
 * ncm_mset_func_list_new_ns_name:
 * @ns: function namespace
 * @name: function name
 * @obj: (allow-none): FIXME
 *
 * FIXME 
 * 
 * Returns: (transfer full): newly created #NcmMSetFuncList.
 */
NcmMSetFuncList *
ncm_mset_func_list_new_ns_name (const gchar *ns, const gchar *name, GObject *obj)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  const gchar *full_ns = NULL;

  G_LOCK (insert_lock);
  {
    gint i;
    for (i = 0; i < flist_class->func_array->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (flist_class->func_array, NcmMSetFuncListStruct, i);
      if (g_str_has_prefix (fdata->ns, ns))
      {
        GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, fdata->ns);
        if ((func_hash != NULL) && g_hash_table_lookup_extended (func_hash, name, NULL, NULL))
        {
          full_ns = fdata->ns;
          break;
        }
      }
    }
  }
  G_UNLOCK (insert_lock);
  g_type_class_unref (flist_class);

  if (full_ns == NULL)
    g_error ("ncm_mset_func_list_new_ns_name: namespace `%s' not found.", ns);

  {
  gchar *full_name = g_strdup_printf ("%s:%s", full_ns, name);
  NcmMSetFuncList *flist = g_object_new (NCM_TYPE_MSET_FUNC_LIST,
                                         "full-name", full_name,
                                         "object", obj,
                                         NULL);
  g_free (full_name);
  return flist;
  }
}

/**
 * ncm_mset_func_list_has_ns_name:
 * @ns: function namespace
 * @name: function name
 *
 * Check if function @name exists in @ns.
 * 
 * Returns: whether the function @name exists in @ns.
 */
gboolean
ncm_mset_func_list_has_ns_name (const gchar *ns, const gchar *name)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  gboolean has_func = FALSE;
  G_LOCK (insert_lock);

  {
    gint i;
    for (i = 0; i < flist_class->func_array->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (flist_class->func_array, NcmMSetFuncListStruct, i);
      if (g_str_has_prefix (fdata->ns, ns))
      {
        GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, fdata->ns);
        if ((func_hash != NULL) && g_hash_table_lookup_extended (func_hash, name, NULL, NULL))
          has_func = TRUE;
      }
    }
  }

  G_UNLOCK (insert_lock);
  
  g_type_class_unref (flist_class);

  return has_func;
}

/**
 * ncm_mset_func_list_has_full_name:
 * @full_name: function full name
 *
 * Check if function @full_name exists.
 * 
 * Returns: whether the function @full_name.
 */
gboolean
ncm_mset_func_list_has_full_name (const gchar *full_name)
{
  gchar **ns_name   = g_strsplit (full_name, ":", 2);
  if (g_strv_length (ns_name) == 2)
  {
    gboolean has_func = ncm_mset_func_list_has_ns_name (ns_name[0], ns_name[1]);

    g_strfreev (ns_name);

    return has_func;
  }
  else
    g_error ("ncm_mset_func_list_has_full_name: invalid full_name `%s'.", full_name);
  
  return FALSE;
}

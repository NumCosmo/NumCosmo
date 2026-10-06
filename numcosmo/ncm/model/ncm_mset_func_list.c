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
 * NcmMSetFuncList:
 *
 * A #NcmMSetFunc chosen by name from a registry of C functions.
 *
 * Model classes register their functions with ncm_mset_func_list_register(), each
 * under a namespace and a name, for instance "NcHICosmo:H". A #NcmMSetFuncList is
 * built from that full name with ncm_mset_func_list_new() and takes the metadata
 * of the registered function. Some functions also need an object of a given type,
 * passed as the #NcmMSetFuncList:object property.
 *
 * ncm_mset_func_list_select() lists the registered functions.
 * ncm_mset_func_list_new_ns_name() and ncm_mset_func_list_has_ns_name() also accept
 * a namespace prefix, such as "NcHICosmo" for a function registered in
 * "NcHICosmoDE"; an exact namespace wins, and a name found under several
 * namespaces with that prefix is ambiguous.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/model/ncm_mset_func_list.h"

enum
{
  PROP_0,
  PROP_FULL_NAME,
  PROP_OBJECT,
};

typedef struct _NcmMSetFuncListPrivate
{
  GType obj_type;
  NcmMSetFuncListN func;
  GObject *obj;
} NcmMSetFuncListPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSetFuncList, ncm_mset_func_list, NCM_TYPE_MSET_FUNC)
G_DEFINE_BOXED_TYPE (NcmMSetFuncListStruct, ncm_mset_func_list_struct, ncm_mset_func_list_struct_copy, ncm_mset_func_list_struct_free)

static void
ncm_mset_func_list_init (NcmMSetFuncList *flist)
{
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);

  self->obj_type = G_TYPE_NONE;
  self->func     = NULL;
  self->obj      = NULL;
}

static void _ncm_mset_func_list_init_from_full_name (NcmMSetFuncList *flist, const gchar *full_name);

static void
_ncm_mset_func_list_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetFuncList *flist              = NCM_MSET_FUNC_LIST (object);
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);
  NcmMSetFunc *func                   = NCM_MSET_FUNC (object);

  g_return_if_fail (NCM_IS_MSET_FUNC_LIST (object));

  switch (prop_id)
  {
    case PROP_FULL_NAME:
    {
      const gchar *full_name = g_value_get_string (value);

      if (full_name == NULL)
        g_error ("_ncm_mset_func_list_set_property: a NcmMSetFuncList needs the full name "
                 "`namespace:name' of a registered function.");

      _ncm_mset_func_list_init_from_full_name (flist, full_name);
      break;
    }
    case PROP_OBJECT:
      g_clear_object (&self->obj);
      self->obj = g_value_dup_object (value);

      if (self->obj != NULL)
      {
        if (!g_type_is_a (G_OBJECT_TYPE (self->obj), self->obj_type))
          g_error ("_ncm_mset_func_list_set_property: function `%s:%s' requires an object of type `%s', but got a `%s'.",
                   ncm_mset_func_peek_ns (func),
                   ncm_mset_func_peek_name (func),
                   g_type_name (self->obj_type),
                   G_OBJECT_TYPE_NAME (self->obj));
      }
      else if (self->obj_type != G_TYPE_NONE)
      {
        g_error ("_ncm_mset_func_list_set_property: function `%s:%s' requires an object of type `%s'.",
                 ncm_mset_func_peek_ns (func),
                 ncm_mset_func_peek_name (func),
                 g_type_name (self->obj_type));
      }

      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_func_list_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetFuncList *flist              = NCM_MSET_FUNC_LIST (object);
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);
  NcmMSetFunc *func                   = NCM_MSET_FUNC (object);

  g_return_if_fail (NCM_IS_MSET_FUNC_LIST (object));

  switch (prop_id)
  {
    case PROP_FULL_NAME:
    {
      gchar *ns_name = g_strdup_printf ("%s:%s",
                                        ncm_mset_func_peek_ns (func),
                                        ncm_mset_func_peek_name (func));

      g_value_take_string (value, ns_name);
      break;
    }
    case PROP_OBJECT:
      g_value_set_object (value, self->obj);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_func_list_dispose (GObject *object)
{
  NcmMSetFuncList *flist              = NCM_MSET_FUNC_LIST (object);
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);

  g_clear_object (&self->obj);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_list_parent_class)->dispose (object);
}

static void _ncm_mset_func_list_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_mset_func_list_class_init (NcmMSetFuncListClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_ncm_mset_func_list_set_property;
  object_class->get_property = &_ncm_mset_func_list_get_property;
  object_class->dispose      = &_ncm_mset_func_list_dispose;

  /**
   * NcmMSetFuncList:full-name:
   *
   * The full name `namespace:name' of the registered function, required at
   * construction.
   */
  g_object_class_install_property (object_class,
                                   PROP_FULL_NAME,
                                   g_param_spec_string ("full-name",
                                                        NULL,
                                                        "Namespace and function name",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetFuncList:object:
   *
   * The object the registered function needs, of the type given at registration,
   * or %NULL for a function that needs none.
   */
  g_object_class_install_property (object_class,
                                   PROP_OBJECT,
                                   g_param_spec_object ("object",
                                                        NULL,
                                                        "object",
                                                        G_TYPE_OBJECT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->func_array = g_array_new (TRUE, TRUE, sizeof (NcmMSetFuncListStruct));
  klass->ns_hash    = g_hash_table_new (g_str_hash, g_str_equal);
  func_class->eval  = &_ncm_mset_func_list_eval;
}

G_LOCK_DEFINE_STATIC (insert_lock);

static void
_ncm_mset_func_list_init_from_full_name (NcmMSetFuncList *flist, const gchar *full_name)
{
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);
  NcmMSetFunc *func                   = NCM_MSET_FUNC (flist);
  NcmMSetFuncListClass *flist_class   = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  gchar **ns_name                     = g_strsplit (full_name, ":", 2);

  if (g_strv_length (ns_name) != 2)
    g_error ("_ncm_mset_func_list_init_from_full_name: invalid full name `%s', expected `namespace:name'.", full_name);

  G_LOCK (insert_lock);
  {
    GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, ns_name[0]);
    gpointer fdata_i;

    if (func_hash == NULL)
      g_error ("_ncm_mset_func_list_init_from_full_name: namespace `%s' not found.", ns_name[0]);

    if (!g_hash_table_lookup_extended (func_hash, ns_name[1], NULL, &fdata_i))
      g_error ("_ncm_mset_func_list_init_from_full_name: name `%s' not found in namespace `%s'.", ns_name[1], ns_name[0]);

    {
      NcmMSetFuncListStruct *fdata = &g_array_index (flist_class->func_array, NcmMSetFuncListStruct, GPOINTER_TO_INT (fdata_i));

      ncm_mset_func_set_meta (func, fdata->name, fdata->symbol, fdata->ns, fdata->desc, fdata->nvar, fdata->dim);

      self->obj_type = fdata->obj_type;
      self->func     = fdata->func;
    }
  }
  G_UNLOCK (insert_lock);

  g_strfreev (ns_name);
  g_type_class_unref (flist_class);
}

static void
_ncm_mset_func_list_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmMSetFuncList *flist              = NCM_MSET_FUNC_LIST (func);
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);

  if ((self->obj_type != G_TYPE_NONE) && (self->obj == NULL))
    g_error ("_ncm_mset_func_list_eval: calling function without object `%s'.", g_type_name (self->obj_type));

  self->func (flist, mset, x, res);
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
  NcmMSetFuncListStruct flist_item  = {
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
 * @ns: (allow-none): namespace prefix
 * @nvar: number of variables, or -1 for any
 * @dim: function dimension, or -1 for any
 *
 * Lists the registered functions whose namespace starts with @ns, or all of them if
 * @ns is %NULL, with @nvar variables and dimension @dim. The strings of the
 * elements belong to the registry and must not be freed.
 *
 * Returns: (transfer container) (element-type NcmMSetFuncListStruct): the matching functions.
 */
GArray *
ncm_mset_func_list_select (const gchar *ns, gint nvar, gint dim)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  GArray *s                         = g_array_new (TRUE, TRUE, sizeof (NcmMSetFuncListStruct));

  G_LOCK (insert_lock);

  g_assert_cmpint (nvar, >=, -1);
  g_assert_cmpint (dim, >=, -1);

  {
    guint i;

    for (i = 0; i < flist_class->func_array->len; i++)
    {
      NcmMSetFuncListStruct *fdata = &g_array_index (flist_class->func_array, NcmMSetFuncListStruct, i);
      gboolean in_ns               = (ns == NULL) || g_str_has_prefix (fdata->ns, ns);
      gboolean in_nvar             = (nvar == -1) || (guint) nvar == fdata->nvar;
      gboolean in_dim              = (dim  == -1) || (guint) dim  == fdata->dim;

      if (in_ns && in_nvar && in_dim)
        g_array_append_val (s, fdata[0]);
    }
  }

  G_UNLOCK (insert_lock);

  g_type_class_unref (flist_class);

  return s;
}

/**
 * ncm_mset_func_list_struct_copy:
 * @fdata: a #NcmMSetFuncListStruct
 *
 * Copies @fdata. The copy shares the strings of @fdata, which belong to the
 * registry and live for the whole program.
 *
 * Returns: (transfer full): a copy of @fdata.
 */
NcmMSetFuncListStruct *
ncm_mset_func_list_struct_copy (const NcmMSetFuncListStruct *fdata)
{
  return g_memdup2 (fdata, sizeof (NcmMSetFuncListStruct));
}

/**
 * ncm_mset_func_list_struct_free:
 * @fdata: a #NcmMSetFuncListStruct
 *
 * Frees a copy made by ncm_mset_func_list_struct_copy(); the shared strings are
 * not freed.
 *
 */
void
ncm_mset_func_list_struct_free (NcmMSetFuncListStruct *fdata)
{
  g_free (fdata);
}

/**
 * ncm_mset_func_list_new:
 * @full_name: function full name
 * @obj: (allow-none): associated object
 *
 * Generates a new instance of #NcmMSetFuncList based on the provided @full_name. The
 * @full_name should adhere to the "namespace:name" format, aligning with a registered
 * function. The associated @obj must match the type as the registered object.
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
 * ncm_mset_func_list_ref:
 * @flist: a #NcmMSetFuncList
 *
 * Increases the reference count of @flist by one.
 *
 * Returns: (transfer full): @flist.
 */
NcmMSetFuncList *
ncm_mset_func_list_ref (NcmMSetFuncList *flist)
{
  return g_object_ref (flist);
}

/**
 * ncm_mset_func_list_free:
 * @flist: a #NcmMSetFuncList
 *
 * Decreases the reference count of @flist by one. If the reference count reaches
 * zero, @flist is freed.
 *
 */
void
ncm_mset_func_list_free (NcmMSetFuncList *flist)
{
  g_object_unref (flist);
}

/**
 * ncm_mset_func_list_clear:
 * @flist: a #NcmMSetFuncList
 *
 * If *@flist is not %NULL, decreases the reference count of *@flist by one and
 * sets *@flist to %NULL.
 *
 */
void
ncm_mset_func_list_clear (NcmMSetFuncList **flist)
{
  g_clear_object (flist);
}

static const gchar *_ncm_mset_func_list_find_ns (NcmMSetFuncListClass *flist_class, const gchar *ns, const gchar *name);

/**
 * ncm_mset_func_list_new_ns_name:
 * @ns: function namespace or namespace prefix
 * @name: function name
 * @obj: (allow-none): associated object
 *
 * Creates a new #NcmMSetFuncList for the function @name. If @ns has no function
 * @name, the function is looked up in the namespaces starting with @ns; it aborts
 * if none or several of them have it. The @obj must have the type given at
 * registration.
 *
 * Returns: (transfer full): newly created #NcmMSetFuncList.
 */
NcmMSetFuncList *
ncm_mset_func_list_new_ns_name (const gchar *ns, const gchar *name, GObject *obj)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  gchar *full_name                  = NULL;

  G_LOCK (insert_lock);
  {
    const gchar *full_ns = _ncm_mset_func_list_find_ns (flist_class, ns, name);

    if (full_ns == NULL)
      g_error ("ncm_mset_func_list_new_ns_name: function `%s' not found in namespace `%s' "
               "nor in any namespace starting with it.",
               name, ns);

    full_name = g_strdup_printf ("%s:%s", full_ns, name);
  }
  G_UNLOCK (insert_lock);
  g_type_class_unref (flist_class);

  {
    NcmMSetFuncList *flist = ncm_mset_func_list_new (full_name, obj);

    g_free (full_name);

    return flist;
  }
}

/**
 * ncm_mset_func_list_has_ns_name:
 * @ns: function namespace or namespace prefix
 * @name: function name
 *
 * Checks if function @name exists in @ns or, as in
 * ncm_mset_func_list_new_ns_name(), in a namespace starting with @ns. Aborts if
 * several namespaces starting with @ns have @name.
 *
 * Returns: whether the function @name exists.
 */
gboolean
ncm_mset_func_list_has_ns_name (const gchar *ns, const gchar *name)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  gboolean has_func;

  G_LOCK (insert_lock);
  has_func = (_ncm_mset_func_list_find_ns (flist_class, ns, name) != NULL);
  G_UNLOCK (insert_lock);

  g_type_class_unref (flist_class);

  return has_func;
}

/*
 * _ncm_mset_func_list_find_ns:
 * @flist_class: the #NcmMSetFuncListClass
 * @ns: namespace or namespace prefix
 * @name: function name
 *
 * Finds the namespace of function @name: @ns itself if it has @name, otherwise the
 * namespace starting with @ns that has @name. Aborts if several such namespaces
 * have @name. Must be called with insert_lock held.
 *
 * Returns: (transfer none) (nullable): the namespace, or %NULL if none has @name.
 */
static const gchar *
_ncm_mset_func_list_find_ns (NcmMSetFuncListClass *flist_class, const gchar *ns, const gchar *name)
{
  GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, ns);

  if ((func_hash != NULL) && g_hash_table_contains (func_hash, name))
    return ns;

  {
    GPtrArray *found     = g_ptr_array_new ();
    const gchar *full_ns = NULL;
    GHashTableIter iter;
    gpointer key, value;

    g_hash_table_iter_init (&iter, flist_class->ns_hash);

    while (g_hash_table_iter_next (&iter, &key, &value))
    {
      if (g_str_has_prefix (key, ns) && g_hash_table_contains (value, name))
        g_ptr_array_add (found, key);
    }

    if (found->len > 1)
    {
      GString *list = g_string_new (NULL);
      guint i;

      g_ptr_array_sort_values (found, (GCompareFunc) g_strcmp0);

      for (i = 0; i < found->len; i++)
        g_string_append_printf (list, "%s`%s'", (i > 0) ? ", " : "", (gchar *) g_ptr_array_index (found, i));

      g_error ("_ncm_mset_func_list_find_ns: function `%s' is ambiguous under namespace `%s', found in %s; "
               "use the full namespace.",
               name, ns, list->str);
    }

    if (found->len == 1)
      full_ns = g_ptr_array_index (found, 0);

    g_ptr_array_unref (found);

    return full_ns;
  }
}

/**
 * ncm_mset_func_list_has_full_name:
 * @full_name: function full name
 *
 * Checks if the function @full_name, of the form `namespace:name', is registered.
 * The namespace must match exactly, as in ncm_mset_func_list_new().
 *
 * Returns: whether the function @full_name exists.
 */
gboolean
ncm_mset_func_list_has_full_name (const gchar *full_name)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  gchar **ns_name                   = g_strsplit (full_name, ":", 2);
  gboolean has_func;

  if (g_strv_length (ns_name) != 2)
    g_error ("ncm_mset_func_list_has_full_name: invalid full name `%s', expected `namespace:name'.", full_name);

  G_LOCK (insert_lock);
  {
    GHashTable *func_hash = g_hash_table_lookup (flist_class->ns_hash, ns_name[0]);

    has_func = (func_hash != NULL) && g_hash_table_contains (func_hash, ns_name[1]);
  }
  G_UNLOCK (insert_lock);

  g_strfreev (ns_name);
  g_type_class_unref (flist_class);

  return has_func;
}

/**
 * ncm_mset_func_list_peek_obj:
 * @flist: #NcmMSetFuncList
 *
 * Gets the object associated with the function.
 *
 * Returns: (transfer none): contained object.
 */
GObject *
ncm_mset_func_list_peek_obj (NcmMSetFuncList *flist)
{
  NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);

  return self->obj;
}


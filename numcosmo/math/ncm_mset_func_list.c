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
 * NcmMSet Functions list.
 *
 * This object is a subclass of #NcmMSetFunc, designed to manage a list of functions. To
 * register these functions, the #ncm_mset_func_list_register function is employed.
 * Selection of functions is accomplished through the use of the
 * ncm_mset_func_list_select() function. Additionally, external objects have the
 * capability to register functions by utilizing the #ncm_mset_func_list_register
 * function.
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


typedef struct _NcmMSetFuncListPrivate
{
  /*< private >*/
  NcmMSetFunc parent_instance;
  GType obj_type;
  NcmMSetFuncListN func;
  GObject *obj;
} NcmMSetFuncListPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSetFuncList, ncm_mset_func_list, NCM_TYPE_MSET_FUNC)

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
      _ncm_mset_func_list_init_from_full_name (flist, g_value_get_string (value));
      break;
    case PROP_OBJECT:
      self->obj = g_value_dup_object (value);

      if (self->obj != NULL)
        g_assert (g_type_is_a (G_OBJECT_TYPE (self->obj), self->obj_type));
      else if (self->obj_type != G_TYPE_NONE)
        g_error ("_ncm_mset_func_list_set_property: object %s:%s requires an object `%s'.",
                 ncm_mset_func_peek_ns (func),
                 ncm_mset_func_peek_name (func),
                 g_type_name (self->obj_type));

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
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
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

  klass->func_array = g_array_new (TRUE, TRUE, sizeof (NcmMSetFuncListStruct));
  klass->ns_hash    = g_hash_table_new (g_str_hash, g_str_equal);
  func_class->eval  = &_ncm_mset_func_list_eval;
}

G_LOCK_DEFINE_STATIC (insert_lock);

static void
_ncm_mset_func_list_init_from_full_name (NcmMSetFuncList *flist, const gchar *full_name)
{
  G_LOCK (insert_lock);
  {
    NcmMSetFuncListPrivate * const self = ncm_mset_func_list_get_instance_private (flist);
    NcmMSetFunc *func                   = NCM_MSET_FUNC (flist);
    NcmMSetFuncListClass *flist_class   = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);

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

        ncm_mset_func_set_meta (func, fdata->name, fdata->symbol, fdata->ns, fdata->desc, fdata->nvar, fdata->dim);

        self->obj_type = fdata->obj_type;
        self->func     = fdata->func;
      }
      else
      {
        g_error ("_ncm_mset_func_list_init_from_full_name: name `%s' not found in namespace `%s'.", ns_name[1], ns_name[0]);
      }
    }

    g_strfreev (ns_name);

    g_type_class_unref (flist_class);
  }
  G_UNLOCK (insert_lock);
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
 * ncm_mset_func_list_new_ns_name:
 * @ns: function namespace
 * @name: function name
 * @obj: (allow-none): associated object
 *
 * Creates a new instance of #NcmMSetFuncList based on the provided @ns and @name.
 * The @obj must match the type as the registered object.
 *
 * Returns: (transfer full): newly created #NcmMSetFuncList.
 */
NcmMSetFuncList *
ncm_mset_func_list_new_ns_name (const gchar *ns, const gchar *name, GObject *obj)
{
  NcmMSetFuncListClass *flist_class = g_type_class_ref (NCM_TYPE_MSET_FUNC_LIST);
  const gchar *full_ns              = NULL;

  G_LOCK (insert_lock);
  {
    guint i;

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
    gchar *full_name       = g_strdup_printf ("%s:%s", full_ns, name);
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
  gboolean has_func                 = FALSE;

  G_LOCK (insert_lock);

  {
    guint i;

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
  gchar **ns_name = g_strsplit (full_name, ":", 2);

  if (g_strv_length (ns_name) == 2)
  {
    gboolean has_func = ncm_mset_func_list_has_ns_name (ns_name[0], ns_name[1]);

    g_strfreev (ns_name);

    return has_func;
  }
  else
  {
    g_error ("ncm_mset_func_list_has_full_name: invalid full_name `%s'.", full_name);
  }

  return FALSE;
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


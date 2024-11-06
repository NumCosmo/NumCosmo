/***************************************************************************
 *            ncm_mset.c
 *
 *  Fri May 25 09:37:54 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * SECTION:ncm_mset
 * @title: NcmMSet
 * @short_description: A set of different NcmModel objects.
 *
 * A #NcmMSet is a set of different #NcmModel objects. It is used to
 * represent a set of models that can be used to fit a data set.
 *
 * When the model class is created the class method ncm_mset_model_register_id()
 * must be used to register the model class. This function must be used once
 * and only once in the model class definition. Any subclasse of the model
 * class will inherit the model id. The same compilation unit must call the macro
 * NCM_MSET_MODEL_REGISTER_ID() for each model class that will be used in the #NcmMSet.
 * It should also include NCM_MSET_MODEL_DECLARE_ID() in the header file.
 *
 * Models can be stackable or not. If a model is stackable, the #NcmMSet can contain
 * more than one instance of the same model. If a model is not stackable, the #NcmMSet
 * can contain only one instance of the model.
 *
 * The model can be a submodel of another model. In this case, the
 * model class must set the main_model_id field to the model id of the
 * parent model. If it is the main model, the main_model_id must be
 * set to NCM_MSET_MODEL_MAIN().
 *
 * The #NcmMSet can be created empty or with a list of models. The stackable models can
 * be added to the #NcmMSet using the function ncm_mset_push() to add the model to the
 * end of the list or using ncm_mset_set_pos() to add the model in a specific position.
 * The non-stackable models can be added using ncm_mset_set(). For both
 * ncm_mset_set_pos() and ncm_mset_set() if there is already a model in the position
 * it will be replaced.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_obj_array.h"

enum
{
  PROP_0,
  PROP_VALID_MAP,
  PROP_MARRAY,
  PROP_FMAP,
  PROP_SIZE,
};

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-mset-error, ncm_mset_error) 
/* *INDENT-ON* */
typedef struct _NcmMSetPrivate
{
  /*< private >*/
  GObject parent_instance;
  NcmObjArray *model_array;
  GHashTable *mid_item_hash;
  GHashTable *model_item_hash;
  GHashTable *fpi_hash;
  GPtrArray *fullname_parray;
  GArray *pi_array;
  GArray *mid_array;
  gboolean valid_map;
  guint total_len;
  guint fparam_len;
  NcmVector *temp_fparams;
} NcmMSetPrivate;

typedef struct _NcmMSetItem
{
  NcmModelID mid;
  NcmModel *model;
  gboolean dup;
  gint added_total_params;
} NcmMSetItem;

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSet, ncm_mset, G_TYPE_OBJECT)
G_DEFINE_BOXED_TYPE (NcmMSetPIndex, ncm_mset_pindex, ncm_mset_pindex_dup, ncm_mset_pindex_free)

static gint
_int_sort (gconstpointer a, gconstpointer b, gpointer user_data)
{
  NcmMSetItem **item_a = (NcmMSetItem **) a;
  NcmMSetItem **item_b = (NcmMSetItem **) b;

  return item_a[0]->mid - item_b[0]->mid;
}

NcmMSetItem *
_ncm_mset_item_new (NcmModel *model, NcmModelID mid)
{
  NcmMSetItem *item = g_slice_new0 (NcmMSetItem);

  item->model              = ncm_model_ref (model);
  item->mid                = mid;
  item->dup                = FALSE;
  item->added_total_params = ncm_model_len (model);

  return item;
}

static void
_ncm_mset_item_free (NcmMSetItem *item)
{
  ncm_model_clear (&item->model);
  g_slice_free (NcmMSetItem, item);
}

static void
ncm_mset_init (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  self->model_array     = g_ptr_array_sized_new (20);
  self->model_item_hash = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, NULL);
  self->mid_item_hash   = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, NULL);
  self->fpi_hash        = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, (GDestroyNotify) g_array_unref);
  self->fullname_parray = g_ptr_array_sized_new (NCM_MSET_INIT_MARRAY);
  self->pi_array        = g_array_sized_new (FALSE, TRUE, sizeof (NcmMSetPIndex), 20);
  self->mid_array       = g_array_sized_new (FALSE, TRUE, sizeof (NcmModelID), 20);

  g_ptr_array_set_free_func (self->model_array, (GDestroyNotify) _ncm_mset_item_free);
  g_ptr_array_set_free_func (self->fullname_parray, g_free);

/* self->fpi_array[i] = g_array_sized_new (FALSE, TRUE, sizeof (gint), 10); */

  self->valid_map = FALSE;
  self->total_len = 0;
}

static void
_ncm_mset_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSet *mset               = NCM_MSET (object);
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_return_if_fail (NCM_IS_MSET (object));

  switch (prop_id)
  {
    case PROP_VALID_MAP:
      g_value_set_boolean (value, self->valid_map);
      break;
    case PROP_MARRAY:
    {
      NcmObjArray *oa = ncm_obj_array_new ();
      guint i;

      for (i = 0; i < self->model_array->len; i++)
      {
        NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

        if (!NCM_MODEL_GET_CLASS (item->model)->is_submodel)
          ncm_obj_array_add (oa, G_OBJECT (item->model));
      }

      g_value_take_boxed (value, oa);
      break;
    }
    case PROP_FMAP:
      g_value_take_boxed (value, ncm_mset_get_fmap (mset));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSet *mset               = NCM_MSET (object);
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_return_if_fail (NCM_IS_MSET (object));

  switch (prop_id)
  {
    case PROP_VALID_MAP:
    {
      const gboolean valid_map = g_value_get_boolean (value);

      if (valid_map)
      {
        if (!self->valid_map)
          ncm_mset_prepare_fparam_map (mset);
      }
      else
      {
        self->valid_map = FALSE;
      }

      break;
    }
    case PROP_MARRAY:
    {
      NcmObjArray *oa = (NcmObjArray *) g_value_get_boxed (value);

      if (oa != NULL)
      {
        guint i;

        for (i = 0; i < oa->len; i++)
        {
          NcmModel *model = NCM_MODEL (ncm_obj_array_peek (oa, i));

          if (NCM_MODEL_GET_CLASS (model)->is_submodel)
            g_error ("_ncm_mset_set_property: NcmMSet model array cannot contain submodels `%s'.",
                     G_OBJECT_TYPE_NAME (model));

          ncm_mset_push (mset, model, NULL);
        }
      }

      break;
    }
    case PROP_FMAP:
    {
      const gchar * const *fmap = g_value_get_boxed (value);

      if (fmap != NULL)
        ncm_mset_set_fmap (mset, fmap, FALSE, NULL);

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_dispose (GObject *object)
{
  NcmMSet *mset               = NCM_MSET (object);
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_clear_pointer (&self->model_array, g_ptr_array_unref);

  g_clear_pointer (&self->mid_item_hash, g_hash_table_unref);
  g_clear_pointer (&self->model_item_hash, g_hash_table_unref);

  g_clear_pointer (&self->fpi_hash, g_hash_table_unref);
  g_clear_pointer (&self->fullname_parray, g_ptr_array_unref);
  g_clear_pointer (&self->pi_array, g_array_unref);
  g_clear_pointer (&self->mid_array, g_array_unref);

  ncm_vector_clear (&self->temp_fparams);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_parent_class)->dispose (object);
}

static void
_ncm_mset_finalize (GObject *object)
{
  /* NcmMSet *mset               = NCM_MSET (object); */
  /* NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset); */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_parent_class)->finalize (object);
}

static void
ncm_mset_class_init (NcmMSetClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_mset_set_property;
  object_class->get_property = &_ncm_mset_get_property;
  object_class->dispose      = &_ncm_mset_dispose;
  object_class->finalize     = &_ncm_mset_finalize;

  g_object_class_install_property (object_class,
                                   PROP_VALID_MAP,
                                   g_param_spec_boolean ("valid-map", NULL, "Valid properties map",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MARRAY,
                                   g_param_spec_boxed ("model-array",
                                                       NULL,
                                                       "NcmModel array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FMAP,
                                   g_param_spec_boxed ("fmap",
                                                       NULL,
                                                       "Free params map",
                                                       G_TYPE_STRV,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  {
    GError *error = NULL;

    klass->ns_table         = g_hash_table_new (g_str_hash, g_str_equal);
    klass->model_desc_array = g_array_new (FALSE, TRUE, sizeof (NcmMSetModelDesc));
    klass->fullname_regex   = g_regex_new ("^\\s*([A-Z][A-Za-z]+)\\:?([0-9]+)?\\:([0-9A-Z\\-a-z_]+)\\s*$", G_REGEX_OPTIMIZE, 0, &error);

    if (error != NULL)
    {
      g_error ("ncm_mset_class_init: %s", error->message);
      g_error_free (error);
    }
  }
}

/**
 * ncm_mset_pindex_new:
 * @mid: Model id
 * @pid: Parameter id
 *
 * Creates a new #NcmMSetPIndex.
 *
 * Returns: (transfer full): a new #NcmMSetPIndex
 */
NcmMSetPIndex *
ncm_mset_pindex_new (NcmModelID mid, guint pid)
{
  NcmMSetPIndex *pi = g_slice_new0 (NcmMSetPIndex);

  pi->mid = mid;
  pi->pid = pid;

  return pi;
}

/**
 * ncm_mset_pindex_dup:
 * @pi: a #NcmMSetPIndex
 *
 * Duplicate a #NcmMSetPIndex.
 *
 * Returns: (transfer full): a new #NcmMSetPIndex with the same values of @pi
 */
NcmMSetPIndex *
ncm_mset_pindex_dup (NcmMSetPIndex *pi)
{
  NcmMSetPIndex *new_pi = ncm_mset_pindex_new (pi->mid, pi->pid);

  return new_pi;
}

/**
 * ncm_mset_pindex_free:
 * @pi: a #NcmMSetPIndex
 *
 * Free a #NcmMSetPIndex.
 *
 */
void
ncm_mset_pindex_free (NcmMSetPIndex *pi)
{
  g_slice_free (NcmMSetPIndex, pi);
}

G_LOCK_DEFINE_STATIC (last_model_id);

/**
 * ncm_mset_model_register_id: (skip)
 * @model_class: a #NcmModelClass
 * @ns: model namespace
 * @desc: short description
 * @long_desc: long description
 * @can_stack: whether the models can stack in a #NcmMSet
 * @main_model_id: main model id, use -1 if this is a main model
 *
 * Register a model class in the #NcmMSet. This function must be used once and only
 * once in the model class definition. Any subclasse of the model class will inherit
 * the model id. The same compilation unit must call the macro
 * NCM_MSET_MODEL_REGISTER_ID() for each model class that will be used in the #NcmMSet.
 * It should also include NCM_MSET_MODEL_DECLARE_ID() in the header file.
 *
 * If @can_stack is TRUE, the models can stack in a #NcmMSet. If @can_stack is FALSE,
 * the #NcmMSet can contain only one instance of the model class or any of its
 * subclasses.
 *
 * If @main_model_id is NCM_MSET_MODEL_MAIN(), this is a main model. If @main_model_id
 * is not NCM_MSET_MODEL_MAIN(), this must be the id of the main model. This is used
 * to define an hierarchy of models. For example, the model class #NcHIPrim is a
 * submodel of #NcHICosmo. The main model of #NcHIPrim is #NcHICosmo. Thus, each
 * instance of #NcHICosmo can contain one instance of #NcHIPrim.
 *
 */
void
ncm_mset_model_register_id (NcmModelClass *model_class, const gchar *ns, const gchar *desc, const gchar *long_desc, gboolean can_stack, NcmModelID main_model_id)
{
  if (model_class->model_id < 0)
  {
    static NcmModelID last_model_id = 0;
    NcmMSetClass *mset_class        = g_type_class_ref (NCM_TYPE_MSET);
    NcmMSetModelDesc *model_desc    = NULL;
    guint id;

    G_LOCK (last_model_id);

    model_class->can_stack = can_stack;

    if (main_model_id > -1)
    {
      model_class->is_submodel   = TRUE;
      model_class->main_model_id = main_model_id;
    }

    model_class->model_id = last_model_id * NCM_MSET_MAX_STACKSIZE;
    id                    = last_model_id;

    last_model_id++;
    g_array_set_size (mset_class->model_desc_array, last_model_id);

    model_desc       = &g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, id);
    model_desc->init = TRUE;

    if (ns == NULL)
      g_error ("ncm_mset_model_register_id: Cannot register model without a namespace.");

    if (desc == NULL)
      g_error ("ncm_mset_model_register_id: Cannot register model without a description.");

    model_desc->ns   = g_strdup (ns);
    model_desc->desc = g_strdup (desc);

    if (long_desc != NULL)
      model_desc->long_desc = g_strdup (long_desc);
    else
      model_desc->long_desc = NULL;

    if (g_hash_table_lookup (mset_class->ns_table, ns) != NULL)
      g_error ("ncm_mset_model_register_id: Model namespace <%s> already registered.", ns);

    g_hash_table_insert (mset_class->ns_table, model_desc->ns, GINT_TO_POINTER (model_class->model_id));

    G_UNLOCK (last_model_id);
    g_type_class_unref (mset_class);
  }
  else
  {
    g_error ("ncm_mset_model_register_id: This model or its parent is already registered, "
             "id = %d. This function must be use once and only in the defining model.",
             model_class->model_id);
  }

  return;
}

/**
 * ncm_mset_split_full_name:
 * @fullname: full name of a parameter
 * @model_ns: (out) (transfer full): model namespace
 * @stackpos_id: (out): stack position id
 * @pname: (out) (transfer full): parameter name
 * @error: a #GError
 *
 * Splits the @fullname into @model_ns, @stackpos_id and @pname. The @fullname
 * should be specified with the parameter full name "model:parameter_name"
 * or "model:stackposition:parameter_name".
 *
 * Returns: %TRUE if the @fullname is valid, %FALSE otherwise.
 */
gboolean
ncm_mset_split_full_name (const gchar *fullname, gchar **model_ns, guint *stackpos_id, gchar **pname, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, FALSE);
  {
    NcmMSetClass * const klass = g_type_class_ref (NCM_TYPE_MSET);
    GMatchInfo *match_info     = NULL;
    gboolean ret               = FALSE;


    if (g_regex_match (klass->fullname_regex, fullname, 0, &match_info))
    {
      gint nm           = g_match_info_get_match_count (match_info);
      gchar *stackpos_s = NULL;

      g_assert_cmpint (nm, ==, 4);

      *model_ns  = g_match_info_fetch (match_info, 1);
      stackpos_s = g_match_info_fetch (match_info, 2);
      *pname     = g_match_info_fetch (match_info, 3);

      if (*stackpos_s != '\0')
      {
        gchar *endptr = NULL;

        *stackpos_id = g_ascii_strtoll (stackpos_s, &endptr, 10);

        if ((*endptr != '\0') || (*stackpos_id >= NCM_MSET_MAX_STACKSIZE))
        {
          ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_INVALID,
                                      "ncm_mset_param_split_full_name: invalid stackpos number (%s >= %d).",
                                      stackpos_s, NCM_MSET_MAX_STACKSIZE);
          g_free (*model_ns);
          g_free (*pname);
          g_free (stackpos_s);

          return FALSE;
        }
      }
      else
      {
        *stackpos_id = 0;
      }

      g_free (stackpos_s);

      ret = TRUE;
    }

    g_match_info_free (match_info);
    g_type_class_unref (klass);

    return ret;
  }
}

/**
 * ncm_mset_empty_new:
 *
 * Creates a new empty #NcmMSet.
 *
 * Returns: (transfer full): a new empty #NcmMSet
 */
NcmMSet *
ncm_mset_empty_new (void)
{
  return g_object_new (NCM_TYPE_MSET, NULL);
}

/**
 * ncm_mset_new:
 * @model0: a #NcmModel
 * @...: a null terminated list of #NcmModel
 *
 * Creates a new #NcmMSet with the models passed as arguments.
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_new (gpointer model0, GError **error, ...)
{
  NcmMSet *mset = NULL;
  va_list ap;

  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  va_start (ap, error);
  mset = ncm_mset_newv (model0, ap, error);
  va_end (ap);

  NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

  return mset;
}

/**
 * ncm_mset_newv:
 * @model0: a #NcmModel
 * @ap: a va_list
 * @error: a #GError
 *
 * Creates a new #NcmMSet with the models passed as arguments.
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_newv (gpointer model0, va_list ap, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  {
    NcmMSet *mset   = ncm_mset_empty_new ();
    NcmModel *model = NULL;


    g_assert (model0 != NULL);
    g_assert (NCM_IS_MODEL (model0));

    ncm_mset_set (mset, model0, error);
    NCM_UTIL_ON_ERROR_RETURN (error, g_clear_object (&mset), NULL);

    while ((model = va_arg (ap, NcmModel *)) != NULL)
    {
      g_assert (NCM_IS_MODEL (model));
      ncm_mset_push (mset, model, error);
      NCM_UTIL_ON_ERROR_RETURN (error, g_clear_object (&mset), NULL);
    }

    return mset;
  }
}

/**
 * ncm_mset_new_array:
 * @model_array: (array) (element-type NcmModel): a #GPtrArray of #NcmModel
 * @error: a #GError
 *
 * Creates a new #NcmMSet with the models passed as arguments.
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_new_array (GPtrArray *model_array, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  {
    NcmMSet *mset = ncm_mset_empty_new ();
    guint i;

    for (i = 0; i < model_array->len; i++)
    {
      NcmModel *model = NCM_MODEL (g_ptr_array_index (model_array, i));

      g_assert (NCM_IS_MODEL (model));
      ncm_mset_push (mset, model, error);
      NCM_UTIL_ON_ERROR_RETURN (error, g_clear_object (&mset), NULL);
    }

    return mset;
  }
}

/**
 * ncm_mset_ref:
 * @mset: a #NcmMSet
 *
 * Increases the reference count of @mset by one.
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_ref (NcmMSet *mset)
{
  return g_object_ref (mset);
}

/**
 * ncm_mset_dup:
 * @mset: a #NcmMSet
 * @ser: a #NcmSerialize
 *
 * Duplicate a #NcmMSet using a #NcmSerialize object.
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_dup (NcmMSet *mset, NcmSerialize *ser)
{
  return NCM_MSET (ncm_serialize_dup_obj (ser, G_OBJECT (mset)));
}

/**
 * ncm_mset_shallow_copy:
 * @mset: a #NcmMSet
 * @error: a #GError
 *
 * Creates a new #NcmMSet with the same models of @mset.
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_shallow_copy (NcmMSet *mset, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  {
    NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
    NcmMSet *mset_sc            = ncm_mset_empty_new ();
    const guint nmodels         = ncm_mset_nmodels (mset);
    guint i;

    for (i = 0; i < nmodels; i++)
    {
      NcmModel *model = ncm_mset_peek_array_pos (mset, i);

      if (!ncm_model_is_submodel (model))
      {
        ncm_mset_push (mset_sc, model, error);
        NCM_UTIL_ON_ERROR_RETURN (error, g_clear_object (&mset_sc), NULL);
      }
    }

    if (self->valid_map)
      ncm_mset_prepare_fparam_map (mset_sc);

    return mset_sc;
  }
}

/**
 * ncm_mset_free:
 * @mset: a #NcmMSet
 *
 * Atomically decreases the reference count of @mset by one. If the reference count drops to 0,
 * all memory allocated by @mset is released.
 *
 */
void
ncm_mset_free (NcmMSet *mset)
{
  g_object_unref (mset);
}

/**
 * ncm_mset_clear:
 * @mset: a #NcmMSet
 *
 * If *@mse is not NULL, decreases the reference count of @mset by one. If the
 * reference count drops to 0, all memory allocated by @mset is released and *@mset is
 * set to NULL.
 *
 */
void
ncm_mset_clear (NcmMSet **mset)
{
  g_clear_object (mset);
}

/**
 * ncm_mset_peek:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 *
 * Peeks a #NcmModel from the #NcmMSet using the model id @mid.
 *
 * Returns: (transfer none): a #NcmModel with the model id @mid
 */
NcmModel *
ncm_mset_peek (NcmMSet *mset, NcmModelID mid)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  NcmMSetItem *item           = g_hash_table_lookup (self->mid_item_hash, GINT_TO_POINTER (mid));

  if (item != NULL)
    return item->model;
  else
    return NULL;
}

/**
 * ncm_mset_fetch:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 *
 * Fetches a #NcmModel from the #NcmMSet using the model id @mid.
 * If the model is not found, an error is set.
 *
 * Returns: (transfer none): a #NcmModel with the model id @mid
 */
NcmModel *
ncm_mset_fetch (NcmMSet *mset, NcmModelID mid, GError **error)
{
  NcmModel *model = ncm_mset_peek (mset, mid);

  if (model == NULL)
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_NOT_SET,
                                "ncm_mset_fetch: model with id %d not found.", mid);


  return model;
}

/**
 * ncm_mset_peek_pos:
 * @mset: a #NcmMSet
 * @base_mid: a #NcmModelID
 * @stackpos_id: position in the stack
 *
 * Peeks a #NcmModel from the #NcmMSet using the model id @base_mid and stack position
 * @stackpos_id. This function is useful when the model is stackable.
 *
 * Returns: (transfer none): a #NcmModel with the model id @base_mid + @stackpos_id
 */
NcmModel *
ncm_mset_peek_pos (NcmMSet *mset, NcmModelID base_mid, guint stackpos_id)
{
  NcmModelID mid = base_mid + stackpos_id;

  return ncm_mset_peek (mset, mid);
}

/**
 * ncm_mset_get:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 *
 * Gets a #NcmModel from the #NcmMSet using the model id @mid.
 *
 * Returns: (transfer full): a #NcmModel with the model id @mid.
 */
NcmModel *
ncm_mset_get (NcmMSet *mset, NcmModelID mid)
{
  NcmModel *model = ncm_mset_peek (mset, mid);

  if (model != NULL)
    return ncm_model_ref (model);
  else
    return NULL;
}

/**
 * ncm_mset_peek_array_pos:
 * @mset: a #NcmMSet
 * @i: array position
 *
 * Peeks a #NcmModel from the #NcmMSet using the array position @i.
 * The array position is not guaranteed to be the same for different
 * #NcmMSet objects. This function is useful to iterate over the models
 * in the #NcmMSet.
 *
 * Returns: (transfer none): a #NcmModel with the array position @i.
 */
NcmModel *
ncm_mset_peek_array_pos (NcmMSet *mset, guint i)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert_cmpuint (i, <, self->model_array->len);

  return ((NcmMSetItem *) g_ptr_array_index (self->model_array, i))->model;
}

/**
 * ncm_mset_peek_by_name:
 * @mset: a #NcmMSet
 * @name: model namespace
 * @error: a #GError
 *
 * Peeks a #NcmModel from the #NcmMSet using the model namespace @name.
 * The name may be specified with the parameter full name "model:stackposition".
 * If the stack position is not specified, the first model with the model namespace
 * @name will be returned.
 *
 * Returns: (transfer none): a #NcmModel with the model namespace @name.
 */
NcmModel *
ncm_mset_peek_by_name (NcmMSet *mset, const gchar *name, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  {
    gchar **ns_stackpos = g_strsplit (name, ":", 2);
    NcmModel *model     = NULL;

    if (ns_stackpos[1] != NULL)
    {
      gchar *endptr     = NULL;
      guint stackpos_id = 0;

      stackpos_id = g_ascii_strtoll (ns_stackpos[1], &endptr, 10);

      if ((*endptr != '\0') || (stackpos_id >= NCM_MSET_MAX_STACKSIZE))
      {
        ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_INVALID,
                                    "ncm_mset_peek_by_name: invalid stackpos number (%s).",
                                    ns_stackpos[1]);

        g_strfreev (ns_stackpos);

        return NULL;
      }

      model = ncm_mset_peek_pos (mset, ncm_mset_get_id_by_ns (ns_stackpos[0]), stackpos_id);
    }
    else
    {
      model = ncm_mset_peek (mset, ncm_mset_get_id_by_ns (name));
    }

    g_strfreev (ns_stackpos);

    return model;
  }
}

/**
 * ncm_mset_fetch_by_name:
 * @mset: a #NcmMSet
 * @name: model namespace
 * @error: a #GError
 *
 * Fetches a #NcmModel from the #NcmMSet using the model namespace @name.
 * The name may be specified with the parameter full name "model:stackposition".
 * If the stack position is not specified, the first model with the model namespace
 * @name will be returned.
 * If the model is not found, an error is set.
 *
 * Returns: (transfer none): a #NcmModel with the model namespace @name.
 */
NcmModel *
ncm_mset_fetch_by_name (NcmMSet *mset, const gchar *name, GError **error)
{
  NcmModel *model = ncm_mset_peek_by_name (mset, name, error);

  NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_fetch_by_name: ");

  if (model == NULL)
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_NOT_SET,
                                "ncm_mset_fetch_by_name: model with name `%s' not found.", name);

  return model;
}

/**
 * ncm_mset_get_mid_array_pos:
 * @mset: a #NcmMSet
 * @i: array position
 *
 * Gets the model id of the model in the array position @i.
 * The array position is not guaranteed to be the same for different
 * #NcmMSet objects. This function is useful to iterate over the models
 * in the #NcmMSet.
 *
 * Returns: a #NcmModelID with the array position @i.
 */
NcmModelID
ncm_mset_get_mid_array_pos (NcmMSet *mset, guint i)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert_cmpuint (i, <, self->model_array->len);

  return ((NcmMSetItem *) g_ptr_array_index (self->model_array, i))->mid;
}

/**
 * ncm_mset_remove:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 *
 * Removes a #NcmModel from the #NcmMSet using the model id @mid.
 *
 */
void
ncm_mset_remove (NcmMSet *mset, NcmModelID mid)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  NcmMSetItem *item           = g_hash_table_lookup (self->mid_item_hash, GINT_TO_POINTER (mid));

  if (item != NULL)
  {
    self->total_len -= item->added_total_params;
    self->valid_map  = FALSE;

    g_hash_table_remove (self->mid_item_hash, GINT_TO_POINTER (mid));

    if (!item->dup)
      g_hash_table_remove (self->model_item_hash, item->model);

    g_assert (g_ptr_array_remove (self->model_array, item));
  }
}

/**
 * ncm_mset_set:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 * @error: a #GError
 *
 * Sets a #NcmModel in the #NcmMSet. If there is already a model with the same
 * model id, it will be replaced. If it is a stackable model, it will be added
 * to the first position.
 *
 */
void
ncm_mset_set (NcmMSet *mset, NcmModel *model, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  ncm_mset_set_pos (mset, model, 0, error);
  NCM_UTIL_ON_ERROR_RETURN (error, , );
}

/**
 * ncm_mset_push:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 * @error: a #GError
 *
 * Pushes a #NcmModel to the end of the #NcmMSet. If the model is not stackable,
 * it will be added to the first position, if there is already a model with the
 * same model id an error will be raised.
 *
 */
void
ncm_mset_push (NcmMSet *mset, NcmModel *model, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  g_assert (model != NULL);
  {
    NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
    NcmModelID base_mid         = ncm_model_id (model);
    guint stackpos_id           = 0;

    while (TRUE)
    {
      NcmModelID mid = base_mid + stackpos_id;

      if (g_hash_table_lookup (self->mid_item_hash, GINT_TO_POINTER (mid)) == NULL)
      {
        ncm_mset_set_pos (mset, model, stackpos_id, error);
        NCM_UTIL_ON_ERROR_RETURN (error, , );

        break;
      }

      stackpos_id++;
    }
  }
}

static void _ncm_mset_set_pos_intern (NcmMSet *mset, NcmModel *model, guint stackpos_id, GError **error);

/**
 * ncm_mset_set_pos:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 * @stackpos_id: stack position
 * @error: a #GError
 *
 * Sets a #NcmModel in the #NcmMSet in the stack position @stackpos_id.
 * If there is already a model with the same model id, it will be replaced.
 * If it is a stackable model, it will be added to the position @stackpos_id.
 * If @stackpos_id is 0, it will be added to the first position.
 *
 * If the @model is not stackable, it will raise an error if @stackpos_id is
 * different from 0.
 *
 */
void
ncm_mset_set_pos (NcmMSet *mset, NcmModel *model, guint stackpos_id, GError **error)
{
  NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);

  g_return_if_fail (error == NULL || *error == NULL);

  if (model_class->is_submodel)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_SUBMODEL,
                                "ncm_mset_set_pos: cannot add model `%s' directly to a NcmMSet.\n"
                                "This model is a submodel of `%s' and it should "
                                "be added to it instead of directly to a NcmMSet",
                                G_OBJECT_TYPE_NAME (model),
                                ncm_mset_get_ns_by_id (model_class->main_model_id));

    return;
  }
  else
  {
    _ncm_mset_set_pos_intern (mset, model, stackpos_id, error);
    NCM_UTIL_ON_ERROR_RETURN (error, , );
  }
}

static void
_ncm_mset_set_pos_intern (NcmMSet *mset, NcmModel *model, guint stackpos_id, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  g_assert (model != NULL);
  {
    NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
    NcmModelID mid              = ncm_model_id (model) + stackpos_id;
    guint model_len             = ncm_model_len (model);
    NcmMSetItem *item           = _ncm_mset_item_new (model, mid);
    GArray *fpi_array;
    NcmMSetItem *item0;
    guint i;

    if ((stackpos_id > 0) && !(NCM_MODEL_GET_CLASS (model)->can_stack))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_NOT_STACKABLE,
                                  "ncm_mset_set_pos: cannot stack object in position %u NcmMSet, type `%s' is not stackable.",
                                  stackpos_id, G_OBJECT_TYPE_NAME (model));
      g_slice_free (NcmMSetItem, item);

      return;
    }

    ncm_mset_remove (mset, mid);

    if ((item0 = g_hash_table_lookup (self->model_item_hash, model)) != NULL)
    {
      if (item->mid > item0->mid)
      {
        item->dup                = TRUE;
        item->added_total_params = 0;
      }
      else
      {
        item0->dup                = TRUE;
        item0->added_total_params = 0;
        g_hash_table_insert (self->model_item_hash, model, item);
      }

      fpi_array = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (item0->mid));
      g_assert (fpi_array != NULL);
      g_array_ref (fpi_array);
    }
    else
    {
      g_hash_table_insert (self->model_item_hash, model, item);
      self->total_len += item->added_total_params;
      fpi_array        = g_array_new (FALSE, TRUE, sizeof (gint));
      g_array_set_size (fpi_array, model_len);

      for (i = 0; i < model_len; i++)
        g_array_index (fpi_array, gint, i) = -1;
    }

    g_hash_table_insert (self->fpi_hash, GINT_TO_POINTER (mid), fpi_array);
    g_hash_table_insert (self->mid_item_hash, GINT_TO_POINTER (mid), item);

    g_ptr_array_add (self->model_array, item);
    g_ptr_array_sort_with_data (self->model_array, _int_sort, NULL);

    if (self->valid_map)
      ncm_mset_prepare_fparam_map (mset);

    for (i = 0; i < ncm_model_get_submodel_len (model); i++)
    {
      NcmModel *submodel = ncm_model_peek_submodel (model, i);

      _ncm_mset_set_pos_intern (mset, submodel, 0, error);
      NCM_UTIL_ON_ERROR_RETURN (error, , );
    }
  }
}

/**
 * ncm_mset_exists:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 *
 * Tests whether a #NcmModel exists in the #NcmMSet.
 *
 * Returns: TRUE if @model exists in @mset, FALSE otherwise.
 */
gboolean
ncm_mset_exists (NcmMSet *mset, NcmModel *model)
{
  NcmModelID mid = ncm_model_id (model);

  if (ncm_mset_peek (mset, mid) != NULL)
    return TRUE;
  else
    return FALSE;
}

/**
 * ncm_mset_exists_pos:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 * @stackpos_id: stack position
 *
 * Checks whether a #NcmModel with the same id as @model and stack position
 * @stackpos_id exists in the #NcmMSet.
 *
 * Returns: TRUE if @model exists in @mset, FALSE otherwise.
 */
gboolean
ncm_mset_exists_pos (NcmMSet *mset, NcmModel *model, guint stackpos_id)
{
  NcmModelID mid = NCM_MSET_MID (ncm_model_id (model), stackpos_id);

  if (ncm_mset_peek (mset, mid) != NULL)
    return TRUE;
  else
    return FALSE;
}

/**
 * ncm_mset_is_subset:
 * @mset: a #NcmMSet
 * @sub_mset: a #NcmMSet
 *
 * Checks whether @sub_mset is a subset of @mset, that is,
 * whether all models in @sub_mset are also in @mset.
 *
 * Returns: TRUE if @sub_mset is a subset of @mset, FALSE otherwise.
 */
gboolean
ncm_mset_is_subset (NcmMSet *mset, NcmMSet *sub_mset)
{
  guint sub_nmodels = ncm_mset_nmodels (sub_mset);
  guint i;

  for (i = 0; i < sub_nmodels; i++)
  {
    NcmModelID mid      = ncm_mset_get_mid_array_pos (sub_mset, i);
    NcmModel *sub_model = ncm_mset_peek_array_pos (sub_mset, i);
    NcmModel *model     = ncm_mset_peek (mset, mid);

    if (sub_model != model)
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_cmp_all:
 * @mset0: a #NcmMSet
 * @mset1: a #NcmMSet
 *
 * Compares two #NcmMSet objects. It returns 0 if they contain the same
 * models, with the same number of parameters, in the same order. It
 * returns -1 if @mset0 has fewer models than @mset1, or if they have
 * the same number of models but @mset0 has a model with fewer
 * parameters than the corresponding model in @mset1. It returns 1 if
 * @mset0 has more models than @mset1, or if they have the same number
 * of models but @mset0 has a model with more parameters than the
 * corresponding model in @mset1.
 *
 * Returns: the comparison result.
 */
gint
ncm_mset_cmp_all (NcmMSet *mset0, NcmMSet *mset1)
{
  guint sub_nmodels0 = ncm_mset_nmodels (mset0);
  guint sub_nmodels1 = ncm_mset_nmodels (mset1);
  guint i;

  if (sub_nmodels0 != sub_nmodels1)
  {
    return sub_nmodels0 < sub_nmodels1 ? -1 : 1;
  }
  else
  {
    for (i = 0; i < sub_nmodels0; i++)
    {
      NcmModelID mid0  = ncm_mset_get_mid_array_pos (mset0, i);
      NcmModelID mid1  = ncm_mset_get_mid_array_pos (mset1, i);
      NcmModel *model0 = ncm_mset_peek_array_pos (mset0, i);
      NcmModel *model1 = ncm_mset_peek_array_pos (mset1, i);

      if (mid0 != mid1)
        return mid0 < mid1 ? -1 : 1;
      else if (!ncm_model_is_equal (model0, model1))
        return ncm_model_len (model0) < ncm_model_len (model1) ? -1 : 1;
    }
  }

  return 0;
}

/**
 * ncm_mset_get_id_by_type:
 * @model_type: a #GType
 *
 * Gets the model id for a model type in the GObject type system.
 *
 * Returns: the model id for @model_type, or -1 if @model_type is not a #NcmModel.
 */
NcmModelID
ncm_mset_get_id_by_type (GType model_type)
{
  g_assert (g_type_is_a (model_type, NCM_TYPE_MODEL));
  {
    NcmMSetClass *mset_class = g_type_class_ref (NCM_TYPE_MSET);
    const gchar *ns          = g_type_name (model_type);
    gpointer model_id;
    gboolean has = g_hash_table_lookup_extended (mset_class->ns_table, ns, NULL, &model_id);

    g_type_class_unref (mset_class);

    if (has)
      return GPOINTER_TO_INT (model_id);
    else
      return -1;
  }
}

/**
 * ncm_mset_get_id_by_ns:
 * @ns: a string containing the model namespace
 *
 * Gets the model id for a model namespace.
 *
 * Returns: the model id for @ns, or -1 if @ns is not a registered model namespace.
 */
NcmModelID
ncm_mset_get_id_by_ns (const gchar *ns)
{
  gpointer model_id;
  NcmMSetClass *mset_class = g_type_class_ref (NCM_TYPE_MSET);
  gboolean has             = g_hash_table_lookup_extended (mset_class->ns_table, ns, NULL, &model_id);

  g_type_class_unref (mset_class);

  if (has)
    return GPOINTER_TO_INT (model_id);
  else
    return -1;
}

/**
 * ncm_mset_get_ns_by_id:
 * @id: namespace id
 *
 * Returns: (transfer none): namespace for @id
 */
const gchar *
ncm_mset_get_ns_by_id (NcmModelID id)
{
  NcmMSetClass *mset_class = g_type_class_ref (NCM_TYPE_MSET);
  guint base_mid           = id / NCM_MSET_MAX_STACKSIZE;

  g_assert_cmpint (base_mid, <, mset_class->model_desc_array->len);
  {
    const gchar *ns = g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, base_mid).ns;

    g_type_class_unref (mset_class);

    return ns;
  }
}

/**
 * ncm_mset_get_type_by_id:
 * @id: namespace id
 *
 * Returns: GType of model @id
 */
GType
ncm_mset_get_type_by_id (NcmModelID id)
{
  NcmMSetClass *mset_class = g_type_class_ref (NCM_TYPE_MSET);
  guint base_mid           = id / NCM_MSET_MAX_STACKSIZE;

  g_assert_cmpint (base_mid, <, mset_class->model_desc_array->len);
  {
    const gchar *ns = g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, base_mid).ns;
    GType t         = g_type_from_name (ns);

    g_type_class_unref (mset_class);

    g_assert_cmpint (t, !=, 0);

    return t;
  }
}

/**
 * ncm_mset_prepare_fparam_map:
 * @mset: a #NcmMSet
 *
 * Computes the free parameters map for @mset. This function must be
 * called before any other function that uses the free parameters map.
 *
 */
void
ncm_mset_prepare_fparam_map (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  self->fparam_len = 0;
  g_array_set_size (self->pi_array, 0);
  g_array_set_size (self->mid_array, 0);
  g_ptr_array_set_size (self->fullname_parray, 0);

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      GArray *fpi_array = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (item->mid));
      gboolean has_free = FALSE;
      gint pid;

      g_assert (fpi_array != NULL);

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        if (ncm_model_param_get_ftype (item->model, pid) == NCM_PARAM_TYPE_FREE)
        {
          NcmMSetPIndex pi = {item->mid, pid};

          g_array_append_val (self->pi_array, pi);
          g_array_index (fpi_array, gint, pid) = self->fparam_len;
          self->fparam_len++;
          has_free = TRUE;
        }
        else
        {
          g_array_index (fpi_array, gint, pid) = -1;
        }
      }

      if (has_free)
        g_array_append_val (self->mid_array, item->mid);
    }
  }

  g_ptr_array_set_size (self->fullname_parray, self->fparam_len);

  ncm_vector_clear (&self->temp_fparams);

  if (self->fparam_len > 0)
    self->temp_fparams = ncm_vector_new (self->fparam_len);

  self->valid_map = TRUE;
}

/**
 * ncm_mset_fparam_map_valid:
 * @mset: a #NcmMSet
 *
 * Checks whether the free parameters map for @mset is valid.
 *
 * Returns: TRUE if the free parameters map for @mset is valid, FALSE otherwise.
 */
gboolean
ncm_mset_fparam_map_valid (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  return self->valid_map;
}

/**
 * ncm_mset_set_fmap:
 * @mset: a #NcmMSet
 * @fmap: (in) (array zero-terminated=1) (element-type utf8): an array of strings
 * @update_models: a boolean
 * @error: a #GError
 *
 * Sets the free parameters map for @mset. This function must be called
 * before any other function that uses the free parameters map. The @fmap
 * array must be zero-terminated and contain the full names of the free
 * parameters in @mset.
 *
 */
void
ncm_mset_set_fmap (NcmMSet *mset, const gchar * const *fmap, gboolean update_models, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  g_assert (fmap != NULL);
  {
    NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
    guint len                   = g_strv_length ((gchar **) fmap);
    guint i;

    for (i = 0; i < self->model_array->len; i++)
    {
      NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

      if (item->dup)
      {
        continue;
      }
      else
      {
        GArray *fpi_array = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (item->mid));
        gint pid;

        g_assert (fpi_array != NULL);

        for (pid = 0; pid < item->added_total_params; pid++)
          g_array_index (fpi_array, gint, pid) = -1;
      }
    }

    self->fparam_len = len;
    g_array_set_size (self->pi_array, 0);
    g_ptr_array_set_size (self->fullname_parray, 0);

    for (i = 0; i < len; i++)
    {
      NcmMSetPIndex *pi = ncm_mset_param_get_by_full_name (mset, fmap[i], error);

      NCM_UTIL_ON_ERROR_RETURN (error, , );

      if (pi == NULL)
      {
        ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_NOT_FOUND,
                                    "ncm_mset_set_fmap: cannot set fmap, invalid param `%s'.", fmap[i]);

        return;
      }
      else
      {
        GArray *fpi_array = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (pi->mid));

        g_array_append_val (self->pi_array, *pi);
        g_array_index (fpi_array, gint, pi->pid) = i;
      }

      ncm_mset_pindex_free (pi);
    }

    g_ptr_array_set_size (self->fullname_parray, self->fparam_len);
    self->valid_map = TRUE;

    if (update_models)
      ncm_mset_param_set_ftype_from_fmap (mset);
  }
}

/**
 * ncm_mset_get_fmap:
 * @mset: a #NcmMSet
 *
 * Gets the free parameters map for @mset. The returned array must be
 * freed with g_strfreev. It contains the full names of the free
 * parameters in @mset.
 *
 * Returns: (transfer full) (array zero-terminated=1) (element-type utf8): an array of strings
 */
gchar **
ncm_mset_get_fmap (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  if (self->valid_map)
  {
    guint len    = ncm_mset_fparam_len (mset);
    gchar **fmap = g_new0 (gchar *, len + 1);
    guint i;

    for (i = 0; i < len; i++)
    {
      fmap[i] = g_strdup (ncm_mset_fparam_full_name (mset, i));
    }

    return fmap;
  }
  else
  {
    return NULL;
  }
}

/**
 * ncm_mset_total_len:
 * @mset: a #NcmMSet
 *
 * Gets the total number of parameters in @mset (including fixed and
 * free parameters).
 *
 * Returns: Total number of parameters in @mset.
 */
guint
ncm_mset_total_len (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  return self->total_len;
}

/**
 * ncm_mset_fparam_len:
 * @mset: a #NcmMSet
 *
 * Gets the number of free parameters in @mset.
 *
 * Returns: Number of free parameters in @mset.
 */
guint
ncm_mset_fparam_len (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  return self->fparam_len;
}

/**
 * ncm_mset_max_param_name:
 * @mset: a #NcmMSet
 *
 * Gets the maximum length of the parameter names in @mset.
 * This function is useful to print the parameters in a pretty way.
 *
 * Returns: Maximum length of the parameter names in @mset.
 */
guint
ncm_mset_max_param_name (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint name_size             = 0;
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        const gchar *pname = ncm_model_param_name (item->model, pid);

        name_size = GSL_MAX (name_size, strlen (pname));
      }
    }
  }

  return name_size;
}

/**
 * ncm_mset_max_fparam_name:
 * @mset: a #NcmMSet
 *
 * Gets the maximum length of the free parameter names in @mset.
 * This function is useful to print the parameters in a pretty way.
 *
 * Returns: Maximum length of the free parameter names in @mset.
 */
guint
ncm_mset_max_fparam_name (NcmMSet *mset)
{
  const guint free_params_len = ncm_mset_fparams_len (mset);
  guint name_size = 0, i;

  for (i = 0; i < free_params_len; i++)
  {
    const gchar *pname = ncm_mset_fparam_name (mset, i);

    name_size = GSL_MAX (name_size, strlen (pname));
  }

  return name_size;
}

/**
 * ncm_mset_max_model_nick:
 * @mset: a #NcmMSet
 *
 * Gets the maximum length of the model nick in @mset.
 * This function is useful to print the models in a pretty way.
 *
 * Returns: Maximum length of the model nick in @mset.
 */
guint
ncm_mset_max_model_nick (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint nick_size             = 0;
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        const gchar *nick = ncm_model_nick (item->model);

        nick_size = GSL_MAX (nick_size, strlen (nick));
      }
    }
  }

  return nick_size;
}

/**
 * ncm_mset_nmodels:
 * @mset: a #NcmMSet
 *
 * Gets the number of models in @mset.
 *
 * Returns: Number of models in @mset.
 */
guint
ncm_mset_nmodels (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  return self->model_array->len;
}

/**
 * ncm_mset_pretty_log:
 * @mset: a #NcmMSet
 *
 * This function prints the contents of @mset. It prints the model
 * nick and parameters' names and their values indicating if they are
 * fixed or free.
 *
 */
void
ncm_mset_pretty_log (NcmMSet *mset)
{
  NcmMSetClass *mset_class    = NCM_MSET_GET_CLASS (mset);
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint name_size             = ncm_mset_max_param_name (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      NcmModelID mid_base = NCM_MSET_GET_BASE_MID (item->mid);
      GArray *fpi_array   = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (item->mid));
      gint pid;

      ncm_cfg_msg_sepa ();
      g_message ("# Model[%05d]:\n#   - %s : %s\n", item->mid, g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).ns, ncm_model_name (item->model));
      ncm_cfg_msg_sepa ();
      g_message ("# Model parameters\n");

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        const gchar *pname = ncm_model_param_name (item->model, pid);

        g_message ("#   - %*s[%02d]: % -20.15g [%s]\n", name_size, pname, pid,
                   ncm_model_param_get (item->model, pid),
                   (g_array_index (fpi_array, gint, pid) >= 0) ? "FREE" : "FIXED");
      }
    }
  }
}

/**
 * ncm_mset_params_pretty_print:
 * @mset: a #NcmMSet
 * @out: name of the file
 * @header: pointer to the command line
 *
 * This function print the command line (first line, commented), the model nick and parameters' names (second line, commented)
 * and their values indicating if they are fixed or free.
 *
 */
void
ncm_mset_params_pretty_print (NcmMSet *mset, FILE *out, const gchar *header)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint name_size             = GSL_MAX (ncm_mset_max_param_name (mset), 5);
  guint nick_size             = ncm_mset_max_model_nick (mset);
  guint i;

  if (header != NULL)
    fprintf (out, "# %s\n ", header);
  else
    fprintf (out, "#\n");

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      const gchar *nick = ncm_model_nick (item->model);
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        const gchar *pname = ncm_model_param_name (item->model, pid);

        fprintf (out, "%*s:%*s %*s % 5.5g\n", nick_size, nick, name_size, pname, name_size,
                 (ncm_model_param_get_ftype (item->model, pid) == NCM_PARAM_TYPE_FREE) ? "FREE" : "FIXED",
                 ncm_model_param_get (item->model, pid));
      }
    }
  }

  return;
}

/**
 * ncm_mset_params_log_vals:
 * @mset: a #NcmMSet
 *
 * Logs the values of the parameters in @mset.
 *
 */
void
ncm_mset_params_log_vals (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  g_message ("#");

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        g_message (" % -20.15g", ncm_model_param_get (item->model, pid));
      }
    }
  }

  g_message ("\n");
}

/**
 * ncm_mset_params_print_vals:
 * @mset: a #NcmMSet
 * @out: a #FILE handler
 *
 * Prints the values of the parameters in @mset in
 * the file @out.
 *
 */
void
ncm_mset_params_print_vals (NcmMSet *mset, FILE *out)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        fprintf (out, " % -20.15g", ncm_model_param_get (item->model, pid));
      }
    }
  }

  fprintf (out, "\n");
}

/**
 * ncm_mset_fparams_log_covar:
 * @mset: a #NcmMSet
 * @covar: a #NcmMatrix
 *
 * Logs the covariance matrix of the free parameters in @mset.
 * The covariance matrix is assumed to be in the same order as
 * the free parameters in @mset and must be square and with
 * the same size as the number of free parameters ncm_mset_fparam_len().
 *
 */
void
ncm_mset_fparams_log_covar (NcmMSet *mset, NcmMatrix *covar)
{
  const guint name_size       = ncm_mset_max_fparam_name (mset);
  const guint free_params_len = ncm_mset_fparam_len (mset);
  const gchar *box            = "---------------";
  guint i, j;

  g_assert (covar != NULL);

  ncm_cfg_msg_sepa ();
  g_message ("# NcmMSet parameters covariance matrix\n");
  g_message ("#                                            ");

  for (i = 0; i < name_size; i++)
    g_message (" ");

  for (i = 0; i < free_params_len; i++)
    i ? g_message ("%s", box) : g_message ("-%s", box);

  if (i)
    g_message ("\n");

  for (i = 0; i < free_params_len; i++)
  {
    const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (mset, i);
    const gchar *pname      = ncm_mset_fparam_name (mset, i);

    g_message ("# %*s[%05d:%02d] = % -12.4g +/- % -12.4g |",
               name_size, pname, pi->mid, pi->pid,
               ncm_mset_fparam_get (mset, i),
               sqrt (ncm_matrix_get (covar, i, i))
              );

    for (j = 0; j < free_params_len; j++)
    {
      g_message (" % -12.4g |", ncm_matrix_get (covar, i, j) / (sqrt (ncm_matrix_get (covar, i, i) * ncm_matrix_get (covar, j, j))));
    }

    g_message ("\n");
  }

  g_message ("#                                            ");

  for (i = 0; i < name_size; i++)
    g_message (" ");

  for (i = 0; i < free_params_len; i++)
    i ? g_message ("%s", box) : g_message ("-%s", box);

  if (i)
    g_message ("\n");

  return;
}

/**
 * ncm_mset_params_valid:
 * @mset: a #NcmMSet
 *
 * Check whenever all models in @mset have valid parameters.
 *
 * Returns: If TRUE all models have valid parameters.
 */
gboolean
ncm_mset_params_valid (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
      continue;
    else if (!ncm_model_params_valid (item->model))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_params_valid_bounds:
 * @mset: a #NcmMSet
 *
 * Check whenever the parameters respect the bounds.
 *
 * Returns: If TRUE the parameter respect the bounds.
 */
gboolean
ncm_mset_params_valid_bounds (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
      continue;
    else if (!ncm_model_params_valid_bounds (item->model))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_cmp:
 * @mset0: a #NcmMSet
 * @mset1: a #NcmMSet
 * @cmp_model: whether to compare if the models correspond to the same objects
 *
 * Compares @mset0 and @mset1 and returns TRUE if both contains the same models types.
 * If @cmp_model is TRUE compare also if the models correspond to the same objects types.
 *
 * Returns: TRUE if @mset0 == @mset1.
 */
gboolean
ncm_mset_cmp (NcmMSet *mset0, NcmMSet *mset1, gboolean cmp_model)
{
  NcmMSetPrivate * const self0 = ncm_mset_get_instance_private (mset0);
  NcmMSetPrivate * const self1 = ncm_mset_get_instance_private (mset1);
  guint i;

  if (self0->model_array->len != self1->model_array->len)
    return FALSE;

  for (i = 0; i < self0->model_array->len; i++)
  {
    NcmMSetItem *item0 = g_ptr_array_index (self0->model_array, i);
    NcmMSetItem *item1 = g_ptr_array_index (self1->model_array, i);

    if (item0->mid != item1->mid)
      return FALSE;

    if (cmp_model && !ncm_model_is_equal (item0->model, item1->model))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_param_set0:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 * @x: the value to set
 *
 * Sets the value of the parameter @pid in the model @mid to @x.
 * This function does not update the model parameters. It is useful
 * when the parameters are being updated in a loop and the model
 * parameters are updated only once, after the loop.
 *
 */
void
ncm_mset_param_set0 (NcmMSet *mset, NcmModelID mid, guint pid, const gdouble x)
{
  ncm_model_param_set0 (ncm_mset_peek (mset, mid), pid, x);
}

/**
 * ncm_mset_param_set:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 * @x: the value to set
 *
 * Sets the value of the parameter @pid in the model @mid to @x.
 * This function updates the model parameters.
 *
 */
void
ncm_mset_param_set (NcmMSet *mset, NcmModelID mid, guint pid, const gdouble x)
{
  ncm_model_param_set (ncm_mset_peek (mset, mid), pid, x);
}

/**
 * ncm_mset_param_get:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the value of the parameter @pid in the model @mid.
 *
 * Returns: the value of the parameter @pid in the model @mid.
 */
gdouble
ncm_mset_param_get (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_orig_param_get:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the value of the original parameter @pid in the model @mid.
 * That is the value of the parameter before any reparametrization.
 *
 * Returns: the value of the original parameter @pid in the model @mid.
 */
gdouble
ncm_mset_orig_param_get (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_orig_param_get (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_name:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the name of the parameter @pid in the model @mid.
 *
 * Returns: the name of the parameter @pid in the model @mid.
 */
const gchar *
ncm_mset_param_name (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_name (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_symbol:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the symbol of the parameter @pid in the model @mid. The
 * parameter symbol is a string that represents the parameter
 * using LaTeX symbols.
 *
 * Returns: the symbol of the parameter @pid in the model @mid.
 */
const gchar *
ncm_mset_param_symbol (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_symbol (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_scale:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the scale of the parameter @pid in the model @mid.
 * This scale is a value that is used as a starting guess
 * for the variation of the parameter in a statistical
 * analysis.
 *
 * Returns: the scale of the parameter @pid in the model @mid.
 */
gdouble
ncm_mset_param_get_scale (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_scale (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_set_scale:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 * @scale: new scale
 *
 * Sets the scale of the parameter @pid in the model @mid to @scale.
 * This scale is a value that is used as a starting guess
 * for the variation of the parameter in a statistical
 * analysis.
 *
 */
void
ncm_mset_param_set_scale (NcmMSet *mset, NcmModelID mid, guint pid, gdouble scale)
{
  ncm_model_param_set_scale (ncm_mset_peek (mset, mid), pid, scale);
}

/**
 * ncm_mset_param_get_lower_bound:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the lower bound of the parameter @pid in the model @mid.
 *
 * Returns: the lower bound of the parameter @pid in the model @mid.
 */
gdouble
ncm_mset_param_get_lower_bound (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_lower_bound (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_upper_bound:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the upper bound of the parameter @pid in the model @mid.
 *
 * Returns: the upper bound of the parameter @pid in the model @mid.
 */
gdouble
ncm_mset_param_get_upper_bound (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_upper_bound (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_get_abstol:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 *
 * Gets the absolute tolerance of the parameter @pid in the model @mid.
 *
 * Returns: the absolute tolerance of the parameter @pid in the model @mid.
 */
gdouble
ncm_mset_param_get_abstol (NcmMSet *mset, NcmModelID mid, guint pid)
{
  return ncm_model_param_get_abstol (ncm_mset_peek (mset, mid), pid);
}

/**
 * ncm_mset_param_set_ftype:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 * @ftype: a #NcmParamType
 *
 * Sets the type of the parameter @pid in the model @mid to @ftype.
 * The parameter type can be fixed (#NCM_PARAM_TYPE_FIXED) or free (#NCM_PARAM_TYPE_FREE).
 *
 */
void
ncm_mset_param_set_ftype (NcmMSet *mset, NcmModelID mid, guint pid, NcmParamType ftype)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  NcmModel *model             = ncm_mset_peek (mset, mid);

  g_assert (model != NULL);
  ncm_model_param_set_ftype (model, pid, ftype);

  if (self->valid_map)
    ncm_mset_prepare_fparam_map (mset);
}

/**
 * ncm_mset_param_set_all_ftype:
 * @mset: a #NcmMSet
 * @ftype: a #NcmParamType
 *
 * Set all parameters of all models to @ftype.
 *
 */
void
ncm_mset_param_set_all_ftype (NcmMSet *mset, NcmParamType ftype)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_model_param_set_ftype (item->model, pid, ftype);
      }
    }
  }

  if (self->valid_map)
    ncm_mset_prepare_fparam_map (mset);
}

/**
 * ncm_mset_param_set_mid_ftype:
 * @mset: a #NcmMSet
 * @mid: model id
 * @ftype: a #NcmParamType
 *
 * Set all parameters of @mid model to @ftype.
 *
 */
void
ncm_mset_param_set_mid_ftype (NcmMSet *mset, NcmModelID mid, NcmParamType ftype)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup || (item->mid != mid))
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_model_param_set_ftype (item->model, pid, ftype);
      }
    }
  }

  if (self->valid_map)
    ncm_mset_prepare_fparam_map (mset);
}

/**
 * ncm_mset_param_set_all_but_mid_ftype:
 * @mset: a #NcmMSet
 * @mid: model id
 * @ftype: a #NcmParamType
 *
 * Set all parameters of all models but @mid to @ftype.
 *
 */
void
ncm_mset_param_set_all_but_mid_ftype (NcmMSet *mset, NcmModelID mid, NcmParamType ftype)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup || (item->mid == mid))
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_model_param_set_ftype (item->model, pid, ftype);
      }
    }
  }

  if (self->valid_map)
    ncm_mset_prepare_fparam_map (mset);
}

/**
 * ncm_mset_param_set_ftype_from_fmap:
 * @mset: a #NcmMSet
 *
 * Set all parameters of all models inside @mset in order
 * to reflect the current fmap.
 *
 */
void
ncm_mset_param_set_ftype_from_fmap (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint i;

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup)
    {
      continue;
    }
    else
    {
      GArray *fpi_array = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (item->mid));

      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        if (g_array_index (fpi_array, gint, pid) == -1)
          continue;

        ncm_model_param_set_ftype (item->model, pid, NCM_PARAM_TYPE_FREE);
      }
    }
  }

  if (self->valid_map)
    ncm_mset_prepare_fparam_map (mset);
}

/**
 * ncm_mset_param_set_vector:
 * @mset: a #NcmMSet
 * @params: a #NcmVector
 *
 * Sets the models parameters using values from the #NcmVector @params.
 *
 */
void
ncm_mset_param_set_vector (NcmMSet *mset, NcmVector *params)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint j                     = 0;
  guint i;

  g_assert_cmpuint (ncm_vector_len (params), ==, ncm_mset_total_len (mset));

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup || (item->added_total_params == 0))
    {
      continue;
    }
    else
    {
      NcmVector *mvec = ncm_vector_get_subvector (params, j, item->added_total_params);

      ncm_model_params_set_vector (item->model, mvec);
      ncm_vector_free (mvec);
      j += item->added_total_params;
    }
  }
}

/**
 * ncm_mset_param_get_vector:
 * @mset: a #NcmMSet
 * @params: a #NcmVector
 *
 * Sets the compontents of @params using the models parameters.
 *
 */
void
ncm_mset_param_get_vector (NcmMSet *mset, NcmVector *params)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint j                     = 0;
  guint i;

  g_assert (ncm_vector_len (params) == ncm_mset_total_len (mset));

  for (i = 0; i < self->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (self->model_array, i);

    if (item->dup || (item->added_total_params == 0))
    {
      continue;
    }
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_vector_set (params, j++, ncm_model_param_get (item->model, pid));
      }
    }
  }
}

/**
 * ncm_mset_param_set_mset:
 * @mset_dest: a #NcmMSet
 * @mset_src: a #NcmMSet
 *
 * Copy parameters from @mset_src to @mset_dest, both #NcmMSet must
 * be compatible.
 *
 */
void
ncm_mset_param_set_mset (NcmMSet *mset_dest, NcmMSet *mset_src)
{
  g_assert_cmpint (ncm_mset_cmp_all (mset_dest, mset_src), ==, 0);
  {
    const guint n     = ncm_mset_total_len (mset_dest);
    NcmVector *params = ncm_vector_new (n);

    ncm_mset_param_get_vector (mset_src, params);
    ncm_mset_param_set_vector (mset_dest, params);

    ncm_vector_free (params);
  }
}

/**
 * ncm_mset_param_get_ftype:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 * @pid: parameter id
 * @error: a #GError
 *
 * Gets the type #NcmParamType of the parameter @pid in the model @mid.
 *
 * Returns: the type #NcmParamType of the parameter @pid in the model @mid.
 */
NcmParamType
ncm_mset_param_get_ftype (NcmMSet *mset, NcmModelID mid, guint pid, GError **error)
{
  NcmModel *model = ncm_mset_peek (mset, mid);

  if (model == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_NOT_SET,
                                "ncm_mset_param_get_ftype: cannot get ftype of model-id %d, "
                                "model not set.", mid);

    return -1;
  }

  return ncm_model_param_get_ftype (model, pid);
}

/**
 * ncm_mset_param_set_pi:
 * @mset: a #NcmMSet
 * @pi: (array length=n) (element-type NcmMSetPIndex): a #NcmMSetPIndex array
 * @x: (array length=n) (element-type double): values to be set
 * @n: number of parameters to set
 *
 * Sets the values of the parameters in @mset using the values in @x.
 * The parameters are identified by the #NcmMSetPIndex @pi. This array
 * and @x must have the same size @n.
 *
 */
void
ncm_mset_param_set_pi (NcmMSet *mset, NcmMSetPIndex *pi, const gdouble *x, guint n)
{
  guint fpi;

  for (fpi = 0; fpi < n; fpi++)
    ncm_mset_param_set (mset, pi[fpi].mid, pi[fpi].pid, x[fpi]);
}

/**
 * ncm_mset_param_get_pi:
 * @mset: a #NcmMSet
 * @pi: (array length=n) (element-type NcmMSetPIndex): a #NcmMSetPIndex array
 * @x: (array length=n) (element-type double): array to store the values
 * @n: number of parameters to get
 *
 * Gets the values of the parameters in @mset and stores them in @x.
 * The parameters are identified by the #NcmMSetPIndex @pi. This array
 * and @x must have the same size @n.
 *
 */
void
ncm_mset_param_get_pi (NcmMSet *mset, NcmMSetPIndex *pi, gdouble *x, guint n)
{
  guint fpi;

  for (fpi = 0; fpi < n; fpi++)
    x[fpi] = ncm_mset_param_get (mset, pi[fpi].mid, pi[fpi].pid);
}

/**
 * ncm_mset_fparams_get_vector:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 *
 * Gets the free parameters of @mset and stores them in @x.
 * The size of @x must be equal to the number of free parameters
 * in @mset ncm_mset_fparams_len().
 *
 */
void
ncm_mset_fparams_get_vector (NcmMSet *mset, NcmVector *x)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint fpi;

  for (fpi = 0; fpi < self->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, fpi);

    ncm_vector_set (x, fpi, ncm_mset_param_get (mset, pi.mid, pi.pid));
  }
}

/**
 * ncm_mset_fparams_get_vector_offset:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 * @offset: starting index
 *
 * Gets the free parameters of @mset and stores them in @x starting
 * at @offset. The size of @x must be equal to the number of free
 * parameters in @mset ncm_mset_fparams_len() plus @offset.
 *
 */
void
ncm_mset_fparams_get_vector_offset (NcmMSet *mset, NcmVector *x, guint offset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint fpi;

  for (fpi = 0; fpi < self->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, fpi);

    ncm_vector_set (x, fpi + offset, ncm_mset_param_get (mset, pi.mid, pi.pid));
  }
}

/**
 * ncm_mset_fparams_set_vector:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 *
 * Sets the free parameters of @mset using the values of @x.
 * The size of @x must be equal to the number of free parameters
 * in @mset ncm_mset_fparams_len().
 *
 */
void
ncm_mset_fparams_set_vector (NcmMSet *mset, const NcmVector *x)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint fpi;

  for (fpi = 0; fpi < self->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, fpi);

    ncm_mset_param_set0 (mset, pi.mid, pi.pid, ncm_vector_get (x, fpi));
  }

  for (fpi = 0; fpi < self->mid_array->len; fpi++)
  {
    NcmModel *model = ncm_mset_peek (mset, g_array_index (self->mid_array, NcmModelID, fpi));

    ncm_model_params_update (model);
  }
}

/**
 * ncm_mset_fparams_set_vector_offset:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 * @offset: starting index
 *
 * Set the free parameters of @mset using the values of @x starting
 * at @offset.
 *
 */
void
ncm_mset_fparams_set_vector_offset (NcmMSet *mset, const NcmVector *x, guint offset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint fpi;

  for (fpi = 0; fpi < self->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, fpi);

    ncm_mset_param_set0 (mset, pi.mid, pi.pid, ncm_vector_get (x, fpi + offset));
  }

  for (fpi = 0; fpi < self->mid_array->len; fpi++)
  {
    NcmModel *model = ncm_mset_peek (mset, g_array_index (self->mid_array, NcmModelID, fpi));

    ncm_model_params_update (model);
  }
}

/**
 * ncm_mset_fparams_set_array:
 * @mset: a #NcmMSet
 * @x: (array) (element-type double): array with the values
 *
 * Sets the free parameters of @mset using the values of @x.
 * The size of @x must be equal to the number of free parameters
 * in @mset ncm_mset_fparams_len(). Otherwise the behaviour is
 * undefined.
 *
 */
void
ncm_mset_fparams_set_array (NcmMSet *mset, const gdouble *x)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint fpi;

  for (fpi = 0; fpi < self->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, fpi);

    ncm_mset_param_set0 (mset, pi.mid, pi.pid, x[fpi]);
  }

  for (fpi = 0; fpi < self->mid_array->len; fpi++)
  {
    NcmModel *model = ncm_mset_peek (mset, g_array_index (self->mid_array, NcmModelID, fpi));

    ncm_model_params_update (model);
  }
}

/**
 * ncm_mset_fparams_set_gsl_vector: (skip)
 * @mset: a #NcmMSet.
 * @x: a #gsl_vector.
 *
 * Sets the free parameters of @mset using the values of @x.
 * The size of @x must be equal to the number of free parameters
 * in @mset ncm_mset_fparams_len().
 *
 */
void
ncm_mset_fparams_set_gsl_vector (NcmMSet *mset, const gsl_vector *x)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint fpi;

  for (fpi = 0; fpi < self->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, fpi);

    ncm_mset_param_set0 (mset, pi.mid, pi.pid, gsl_vector_get (x, fpi));
  }

  for (fpi = 0; fpi < self->mid_array->len; fpi++)
  {
    NcmModel *model = ncm_mset_peek (mset, g_array_index (self->mid_array, NcmModelID, fpi));

    ncm_model_params_update (model);
  }
}

/**
 * ncm_mset_fparams_len:
 * @mset: a #NcmMSet
 *
 * Gets the number of free parameters in @mset.
 *
 * Returns: the number of free parameters in @mset.
 */
guint
ncm_mset_fparams_len (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map);

  return self->fparam_len;
}

/**
 * ncm_mset_fparam_name:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the name of the @n-th free parameter.
 *
 * Returns: (transfer none): the name of the @n-th free parameter.
 */
const gchar *
ncm_mset_fparam_name (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_name (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_symbol:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the symbol of the @n-th free parameter.
 *
 * Returns: (transfer none): the symbol of the @n-th free parameter.
 */
const gchar *
ncm_mset_fparam_symbol (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_symbol (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_full_name:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the full name of the @n-th free parameter. That is
 * the name of the model, the stack position and the parameter
 * name.
 *
 * Returns: (transfer none): the full name of the @n-th free parameter.
 */
const gchar *
ncm_mset_fparam_full_name (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  gchar *fullname;

  g_assert (self->valid_map && n < self->fparam_len);

  if ((fullname = g_ptr_array_index (self->fullname_parray, n)) != NULL)
  {
    return fullname;
  }
  else
  {
    NcmMSetPIndex pi        = g_array_index (self->pi_array, NcmMSetPIndex, n);
    const gchar *model_ns   = ncm_mset_get_ns_by_id (pi.mid); /* ncm_model_nick (ncm_mset_peek (mset, pi.mid));*/
    const gchar *pname      = ncm_mset_param_name (mset, pi.mid, pi.pid);
    const guint stackpos_id = pi.mid % NCM_MSET_MAX_STACKSIZE;

    if (stackpos_id > 0)
      fullname = g_strdup_printf ("%s:%02u:%s", model_ns, stackpos_id, pname);
    else
      fullname = g_strdup_printf ("%s:%s", model_ns, pname);

    g_ptr_array_index (self->fullname_parray, n) = fullname;
  }

  return fullname;
}

/**
 * ncm_mset_param_get_by_full_name:
 * @mset: a #NcmMSet
 * @fullname: param's full name
 * @error: a #GError
 *
 * Gets the #NcmMSetPIndex of the parameter identified by @fullname.
 * The @fullname must be in the form "model:stackpos:param_name" when
 * the model has a stack or "model:param_name" when the model has no
 * stack.
 *
 * Returns: (transfer full): the #NcmMSetPIndex of the parameter identified by @fullname.
 */
NcmMSetPIndex *
ncm_mset_param_get_by_full_name (NcmMSet *mset, const gchar *fullname, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  {
    NcmMSetPIndex *pi        = NULL;
    gchar *model_ns          = NULL;
    gchar *pname             = NULL;
    guint stackpos_id        = 0;
    gboolean full_name_found = FALSE;

    full_name_found = ncm_mset_split_full_name (fullname, &model_ns, &stackpos_id, &pname, error);
    NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

    if (full_name_found)
    {
      guint pid     = 0;
      gchar *endptr = NULL;
      NcmModelID mid;
      NcmModel *model;

      mid = ncm_mset_get_id_by_ns (model_ns);

      if (mid < 0)
      {
        ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_NOT_FOUND,
                                    "ncm_mset_param_get_by_full_name: namespace `%s' not found.",
                                    model_ns);
        g_free (pname);
        g_free (model_ns);

        return NULL;
      }

      mid  += stackpos_id;
      model = ncm_mset_peek (mset, mid);

      if (model != NULL)
      {
        endptr = NULL;
        pid    = g_ascii_strtoll (pname, &endptr, 10);

        if (*endptr != '\0')
        {
          gboolean found = ncm_model_param_index_from_name (model, pname, &pid, error);

          NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

          if (found)
            pi = ncm_mset_pindex_new (mid, pid);
        }
        else if (pid < ncm_model_len (model))
        {
          pi = ncm_mset_pindex_new (mid, pid);
        }
      }

      g_free (pname);
      g_free (model_ns);
    }
    else
    {
      ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_INVALID,
                                  "ncm_mset_param_get_by_full_name: invalid full name `%s'.",  fullname);

      return NULL;
    }

    return pi;
  }
}

/**
 * ncm_mset_fparam_get_scale:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the scale of the @n-th free parameter.
 *
 * Returns: the scale of the @n-th free parameter.
 */
gdouble
ncm_mset_fparam_get_scale (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_get_scale (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_lower_bound:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the lower bound of the @n-th free parameter.
 *
 * Returns: the lower bound of the @n-th free parameter.
 */
gdouble
ncm_mset_fparam_get_lower_bound (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_get_lower_bound (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_upper_bound:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the upper bound of the @n-th free parameter.
 *
 * Returns: the upper bound of the @n-th free parameter.
 */
gdouble
ncm_mset_fparam_get_upper_bound (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_get_upper_bound (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_bound_matrix:
 * @mset: a #NcmMSet
 *
 * Gets a matrix with the lower and upper bounds of all free parameters.
 * The returned matrix has two columns, the first column contains the
 * lower bounds and the second column contains the upper bounds.
 * The number of rows is equal to the number of free parameters in @mset.
 *
 * Returns: (transfer full): a matrix with the lower and upper bounds of all free parameters.
 */
NcmMatrix *
ncm_mset_fparam_get_bound_matrix (NcmMSet *mset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  NcmMatrix *bounds           = ncm_matrix_new (self->fparam_len, 2);
  guint i;

  g_assert (self->valid_map);

  for (i = 0; i < self->fparam_len; i++)
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, i);
    const gdouble ub_i     = ncm_mset_param_get_upper_bound (mset, pi.mid, pi.pid);
    const gdouble lb_i     = ncm_mset_param_get_lower_bound (mset, pi.mid, pi.pid);

    ncm_matrix_set (bounds, i, 0, lb_i);
    ncm_matrix_set (bounds, i, 1, ub_i);
  }

  return bounds;
}

/**
 * ncm_mset_fparam_get_abstol:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the absolute tolerance of the @n-th free parameter.
 *
 * Returns: the absolute tolerance of the @n-th free parameter.
 */
gdouble
ncm_mset_fparam_get_abstol (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_get_abstol (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_set_scale:
 * @mset: a #NcmMSet
 * @n: free parameter index
 * @scale: new scale
 *
 * Sets the scale of the @n-th free parameter to @scale.
 *
 */
void
ncm_mset_fparam_set_scale (NcmMSet *mset, guint n, gdouble scale)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    ncm_mset_param_set_scale (mset, pi.mid, pi.pid, scale);
  }
}

/**
 * ncm_mset_fparam_valid_bounds:
 * @mset: a #NcmMSet
 * @theta: free parameters vector
 *
 * Checks if the values of @theta respect the parameter bounds.
 *
 * Returns: whether @theta contain values respecting the parameter bounds.
 */
gboolean
ncm_mset_fparam_valid_bounds (NcmMSet *mset, NcmVector *theta)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map);
  g_assert_cmpuint (ncm_vector_len (theta), ==, self->fparam_len);
  {
    guint i;

    for (i = 0; i < self->fparam_len; i++)
    {
      const gdouble lb  = ncm_mset_fparam_get_lower_bound (mset, i);
      const gdouble ub  = ncm_mset_fparam_get_upper_bound (mset, i);
      const gdouble val = ncm_vector_get (theta, i);

      if ((val < lb) || (val > ub))
        return FALSE;
    }

    return TRUE;
  }
}

/**
 * ncm_mset_fparam_valid_bounds_offset:
 * @mset: a #NcmMSet
 * @theta: free parameters vector
 * @offset: starting index
 *
 * Checks if the values of @theta respect the parameter bounds.
 * The values are checked starting at @offset.
 *
 * Returns: whether @theta contain values respecting the parameter bounds.
 */
gboolean
ncm_mset_fparam_valid_bounds_offset (NcmMSet *mset, NcmVector *theta, guint offset)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map);
  g_assert_cmpuint (ncm_vector_len (theta), ==, self->fparam_len + offset);
  {
    guint i;

    for (i = 0; i < self->fparam_len; i++)
    {
      const gdouble lb  = ncm_mset_fparam_get_lower_bound (mset, i);
      const gdouble ub  = ncm_mset_fparam_get_upper_bound (mset, i);
      const gdouble val = ncm_vector_get (theta, i + offset);

      if ((val < lb) || (val > ub))
        return FALSE;
    }

    return TRUE;
  }
}

/**
 * ncm_mset_fparam_validate_all:
 * @mset: a #NcmMSet
 * @theta: free parameters vector
 *
 * Checks if the values of @theta respect all requirements.
 *
 * Returns: whether @theta contain values respecting all requirements.
 */
gboolean
ncm_mset_fparam_validate_all (NcmMSet *mset, NcmVector *theta)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  gboolean valid;

  g_assert (self->valid_map);
  g_assert_cmpuint (ncm_vector_len (theta), ==, self->fparam_len);

  ncm_mset_fparams_get_vector (mset, self->temp_fparams);
  ncm_mset_fparams_set_vector (mset, theta);

  valid = ncm_mset_params_valid (mset) && ncm_mset_params_valid_bounds (mset);
  ncm_mset_fparams_set_vector (mset, self->temp_fparams);

  return valid;
}

/**
 * ncm_mset_fparam_get:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the value of the @n-th free parameter.
 *
 * Returns: the value of the @n-th free parameter.
 */
gdouble
ncm_mset_fparam_get (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_get (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_set:
 * @mset: a #NcmMSet
 * @n: free parameter index
 * @x: new value
 *
 * Sets the value of the @n-th free parameter to @x.
 *
 */
void
ncm_mset_fparam_set (NcmMSet *mset, guint n, const gdouble x)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (self->pi_array, NcmMSetPIndex, n);

    return ncm_mset_param_set (mset, pi.mid, pi.pid, x);
  }
}

/**
 * ncm_mset_fparam_get_pi:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * Gets the #NcmMSetPIndex of the @n-th free parameter.
 *
 * Returns: (transfer none): the #NcmMSetPIndex of the @n-th free parameter.
 */
const NcmMSetPIndex *
ncm_mset_fparam_get_pi (NcmMSet *mset, guint n)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map && n < self->fparam_len);

  return &(g_array_index (self->pi_array, NcmMSetPIndex, n));
}

/**
 * ncm_mset_fparam_get_fpi:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 * @pid: parameter id
 *
 * Gets the free parameter index of the parameter @pid in the model @mid.
 *
 * Returns: the free parameter index of the parameter @pid in the model @mid.
 */
gint
ncm_mset_fparam_get_fpi (NcmMSet *mset, NcmModelID mid, guint pid)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);

  g_assert (self->valid_map);

  {
    GArray *fpi_array = g_hash_table_lookup (self->fpi_hash, GINT_TO_POINTER (mid));

    g_assert (fpi_array != NULL);

    return g_array_index (fpi_array, gint, pid);
  }
}

/**
 * ncm_mset_fparam_get_pi_by_name:
 * @mset: a #NcmMSet
 * @name: parameter name
 * @error: a #GError
 *
 * Gets the #NcmMSetPIndex of the parameter identified by @name.
 * The name can be the parameter name or the full name.
 *
 * Returns: (transfer none): the #NcmMSetPIndex of the parameter identified by @name.
 */
const NcmMSetPIndex *
ncm_mset_fparam_get_pi_by_name (NcmMSet *mset, const gchar *name, GError **error)
{
  NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
  guint match                 = 0;
  guint match_i               = 0;
  guint i;

  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  g_assert (self->valid_map);

  for (i = 0; i < self->fparam_len; i++)
  {
    const gchar *name_i = ncm_mset_fparam_name (mset, i);

    if (strcmp (name, name_i) == 0)
    {
      match++;
      match_i = i;
    }
  }

  if (match > 1)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_PARAM_NAME_AMBIGUOUS,
                                "ncm_mset_fparam_get_pi_by_name: more than one [%u] "
                                "parameters with the same name %s, "
                                "use the full name to avoid ambiguities.",
                                match, name);

    return NULL;
  }
  else if (match == 1)
  {
    return ncm_mset_fparam_get_pi (mset, match_i);
  }
  else
  {
    for (i = 0; i < self->fparam_len; i++)
    {
      const gchar *name_i = ncm_mset_fparam_full_name (mset, i);

      if (strcmp (name, name_i) == 0)
      {
        match++;
        match_i = i;
      }
    }

    g_assert_cmpuint (match, <=, 1);

    if (match == 0)
      return NULL;
    else
      return ncm_mset_fparam_get_pi (mset, match_i);
  }
}

/**
 * ncm_mset_save:
 * @mset: a #NcmMSet
 * @ser: a #NcmSerialize
 * @filename: a filename
 * @save_comment: whether to save comments
 * @error: a #GError
 *
 * Saves the #NcmMSet to a file using #GKeyFile.
 *
 */
void
ncm_mset_save (NcmMSet *mset, NcmSerialize *ser, const gchar *filename, gboolean save_comment, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  {
    NcmMSetPrivate * const self = ncm_mset_get_instance_private (mset);
    NcmMSetClass *mset_class    = NCM_MSET_GET_CLASS (mset);
    GKeyFile *msetfile          = g_key_file_new ();
    guint i;

    {
      GError *local_error        = NULL;
      gchar *mset_desc           = ncm_cfg_string_to_comment ("NcmMSet");
      gchar *mset_valid_map_desc = ncm_cfg_string_to_comment ("valid-map property");

      g_key_file_set_boolean (msetfile, "NcmMSet", "valid_map", self->valid_map);

      if (save_comment)
      {
        if (!g_key_file_set_comment (msetfile, "NcmMSet", NULL, mset_desc, &local_error))
        {
          ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: ");

          NCM_UTIL_ON_ERROR_RETURN (error, , );
        }

        if (!g_key_file_set_comment (msetfile, "NcmMSet", "valid_map", mset_valid_map_desc, &local_error))
        {
          ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: ");

          NCM_UTIL_ON_ERROR_RETURN (error, , );
        }
      }

      g_free (mset_desc);
      g_free (mset_valid_map_desc);
    }

    for (i = 0; i < self->model_array->len; i++)
    {
      NcmMSetItem *item       = g_ptr_array_index (self->model_array, i);
      NcmModelID mid_base     = item->mid / NCM_MSET_MAX_STACKSIZE;
      const guint stackpos_id = item->mid % NCM_MSET_MAX_STACKSIZE;
      const gchar *ns         = g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).ns;
      GError *local_error     = NULL;
      gchar *group;

      if (stackpos_id > 0)
        group = g_strdup_printf ("%s:%02u", ns, stackpos_id);
      else
        group = g_strdup_printf ("%s", ns);

      if (item->dup)
      {
        NcmMSetItem *item0 = g_hash_table_lookup (self->model_item_hash, item->model);

        g_assert (item0 != NULL);

        g_key_file_set_uint64 (msetfile, group, ns, item0->mid % NCM_MSET_MAX_STACKSIZE);

        if (save_comment)
        {
          gchar *model_desc = ncm_cfg_string_to_comment (g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).desc);

          if (!g_key_file_set_comment (msetfile, group, NULL, model_desc, &local_error))
          {
            ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: ");

            NCM_UTIL_ON_ERROR_RETURN (error, , );
          }

          g_free (model_desc);
        }
      }
      else
      {
        NcmModelClass *model_class = NCM_MODEL_GET_CLASS (item->model);
        GObjectClass *oclass       = G_OBJECT_CLASS (model_class);
        GVariant *model_var        = ncm_serialize_to_variant (ser, G_OBJECT (item->model));
        guint nsubmodels           = ncm_model_get_submodel_len (item->model);
        GVariant *params           = NULL;
        gchar *obj_name            = NULL;
        guint nparams, j;

        g_variant_get (model_var, NCM_SERIALIZE_OBJECT_FORMAT, &obj_name, &params);
        nparams = g_variant_n_children (params);
        g_key_file_set_value (msetfile, group, ns, obj_name);

        if (save_comment)
        {
          gchar *model_desc = ncm_cfg_string_to_comment (g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).desc);

          if (!g_key_file_set_comment (msetfile, group, NULL, model_desc, &local_error))
          {
            ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: ");

            NCM_UTIL_ON_ERROR_RETURN (error, , );
          }

          g_free (model_desc);
        }

        for (j = 0; j < nsubmodels; j++)
        {
          NcmModel *submodel = ncm_model_peek_submodel (item->model, j);

          ncm_serialize_unset (ser, submodel);
          ncm_serialize_remove_ser (ser, submodel);
        }

        if (nparams != 0)
        {
          GVariantIter iter;
          GVariant *value;
          gchar *key;

          g_variant_iter_init (&iter, params);

          while (g_variant_iter_next (&iter, "{sv}", &key, &value))
          {
            if (strcmp (key, "submodel-array") == 0)
            {
              g_variant_unref (value);
              g_free (key);
              continue;
            }
            else
            {
              GParamSpec *param_spec = g_object_class_find_property (oclass, key);
              gchar *param_str       = g_variant_print (value, TRUE);

              if (param_spec == NULL)
              {
                ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_PROPERTY_NOT_FOUND,
                                            "ncm_mset_save: property `%s' not found in object `%s'.", key, obj_name);

                g_free (obj_name);
                g_free (group);
                g_variant_unref (params);
                g_variant_unref (model_var);

                return;
              }

              g_key_file_set_value (msetfile, group, key, param_str);

              if (save_comment)
              {
                const gchar *blurb = g_param_spec_get_blurb (param_spec);

                if ((blurb != NULL) && (blurb[0] != 0))
                {
                  gchar *desc = ncm_cfg_string_to_comment (blurb);

                  if (!g_key_file_set_comment (msetfile, group, key, desc, &local_error))
                  {
                    ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: ");

                    NCM_UTIL_ON_ERROR_RETURN (error, , );
                  }

                  g_free (desc);
                }
              }

              g_variant_unref (value);
              g_free (key);
              g_free (param_str);
            }
          }
        }

        g_free (obj_name);
        g_variant_unref (params);
        g_variant_unref (model_var);
      }

      g_free (group);
    }

    {
      GError *local_error = NULL;
      gsize len           = 0;
      gchar *mset_data    = g_key_file_to_data (msetfile, &len, &local_error);

      ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: Error converting NcmMSet to configuration file: ");

      NCM_UTIL_ON_ERROR_RETURN (error,
                                g_free (mset_data);
                                g_key_file_free (msetfile), );

      if (!g_file_set_contents (filename, mset_data, len, &local_error))
      {
        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_save: Error saving configuration file to disk: ");

        NCM_UTIL_ON_ERROR_RETURN (error, g_free (mset_data);
                                  g_key_file_free (msetfile), );
      }

      g_free (mset_data);
      g_key_file_free (msetfile);
    }
  }
}

/**
 * ncm_mset_load: (constructor)
 * @filename: mset filename
 * @ser: a #NcmSerialize
 * @error: a #GError
 *
 * Loads a #NcmMSet from a configuration file using the #NcmSerialize
 * object @ser. The file must be in the same format as the one
 * generated by ncm_mset_save().
 *
 * Returns: (transfer full): the loaded #NcmMSet.
 */
NcmMSet *
ncm_mset_load (const gchar *filename, NcmSerialize *ser, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);
  {
    NcmMSet *mset             = ncm_mset_empty_new ();
    GKeyFile *msetfile        = g_key_file_new ();
    gchar **groups            = NULL;
    gsize ngroups             = 0;
    gboolean valid_map        = FALSE;
    GPtrArray *submodel_array = g_ptr_array_new ();
    GError *local_error       = NULL;
    guint i;

    g_ptr_array_set_free_func (submodel_array, (GDestroyNotify) ncm_model_free);

    if (!g_key_file_load_from_file (msetfile, filename, G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS, &local_error))
    {
      ncm_util_forward_or_call_error (error, local_error,
                                      "ncm_mset_load: Invalid mset configuration file: %s: ",
                                      filename);

      return NULL;
    }

    if (g_key_file_has_group (msetfile, "NcmMSet"))
    {
      if (g_key_file_has_key (msetfile, "NcmMSet", "valid_map", &local_error))
      {
        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
        NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

        valid_map = g_key_file_get_boolean (msetfile, "NcmMSet", "valid_map", &local_error);
        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
        NCM_UTIL_ON_ERROR_RETURN (error, , NULL);
      }

      g_key_file_remove_group (msetfile, "NcmMSet", &local_error);
      ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
      NCM_UTIL_ON_ERROR_RETURN (error, , NULL);
    }

    groups = g_key_file_get_groups (msetfile, &ngroups);

    for (i = 0; i < ngroups; i++)
    {
      GString *obj_ser = g_string_sized_new (200);
      gchar *ns        = g_strdup (groups[i]);
      gchar *twopoints = g_strrstr (ns, ":");
      guint stackpos   = 0;

      if (twopoints != NULL)
      {
        *twopoints = '\0';
        stackpos   = atoi (++twopoints);
      }

      if (!g_key_file_has_key (msetfile, groups[i], ns, &local_error))
      {
        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
        NCM_UTIL_ON_ERROR_RETURN (error, g_free (ns), NULL);

        ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_KEY_FILE_INVALID,
                                    "ncm_mset_load: Every group must contain a key with same name "
                                    "indicating the object type `%s' `%s'.", groups[i], ns);

        return NULL;
      }

      {
        gchar *obj_type = g_key_file_get_value (msetfile, groups[i], ns, &local_error);

        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
        NCM_UTIL_ON_ERROR_RETURN (error, g_free (obj_type), NULL);

        if (strlen (obj_type) < 5)
        {
          gchar *endptr       = NULL;
          guint stackpos_orig = g_ascii_strtoll (obj_type, &endptr, 10);

          g_assert (endptr != NULL);

          if (*endptr == '\0')
          {
            NcmModelID id = ncm_mset_get_id_by_ns (ns);

            g_assert_cmpint (id, >, -1);
            {
              NcmModel *model = ncm_mset_peek_pos (mset, id, stackpos_orig);

              g_assert (model != NULL);
              ncm_mset_set_pos (mset, model, stackpos, error);

              NCM_UTIL_ON_ERROR_FORWARD (error, g_free (obj_type), NULL, "ncm_mset_load: ");
            }

            g_free (ns);
            g_free (obj_type);
            g_string_free (obj_ser, TRUE);
            continue;
          }
        }

        g_string_append_printf (obj_ser, "(\'%s\', @a{sv} {", obj_type);
        g_free (obj_type);
        g_key_file_remove_key (msetfile, groups[i], ns, &local_error);
        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
        NCM_UTIL_ON_ERROR_RETURN (error, , NULL);
      }

      {
        gsize nkeys  = 0;
        gchar **keys = g_key_file_get_keys (msetfile, groups[i], &nkeys, &local_error);
        guint j;

        ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
        NCM_UTIL_ON_ERROR_RETURN (error, , NULL);


        for (j = 0; j < nkeys; j++)
        {
          if (strcmp (keys[j], "submodel-array") == 0)
          {
            ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_KEY_FILE_INVALID,
                                        "ncm_mset_load: serialized version of mset models cannot "
                                        "contain the submodel-array property.");

            return NULL;
          }
          else
          {
            gchar *propval = g_key_file_get_value (msetfile, groups[i], keys[j], &local_error);

            ncm_util_forward_or_call_error (error, local_error, "ncm_mset_load: ");
            NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

            g_string_append_printf (obj_ser, "\'%s\':<%s>", keys[j], propval);
            g_free (propval);

            if (j + 1 != nkeys)
              g_string_append (obj_ser, ", ");
          }
        }

        g_string_append (obj_ser, "})");
        g_strfreev (keys);
      }

      {
        GObject *obj = ncm_serialize_from_string (ser, obj_ser->str);

        g_assert (NCM_IS_MODEL (obj));

        if (ncm_mset_exists_pos (mset, NCM_MODEL (obj), stackpos))
        {
          ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_ALREADY_SET,
                                      "ncm_mset_load: a model ``%s'' already exists in NcmMSet.", obj_ser->str);
          g_object_unref (obj);
          g_string_free (obj_ser, TRUE);
          g_free (ns);

          return NULL;
        }

        if (ncm_model_is_submodel (NCM_MODEL (obj)))
        {
          NcmModel *submodel = ncm_model_ref (NCM_MODEL (obj));

          g_ptr_array_add (submodel_array, submodel);
        }
        else
        {
          ncm_mset_set_pos (mset, NCM_MODEL (obj), stackpos, error);
          NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_load: ");
        }

        ncm_model_free (NCM_MODEL (obj));
      }
      g_string_free (obj_ser, TRUE);
      g_free (ns);
    }

    for (i = 0; i < submodel_array->len; i++)
    {
      NcmModel *submodel  = g_ptr_array_index (submodel_array, i);
      NcmModelID mid      = ncm_model_main_model (submodel);
      NcmModel *mainmodel = ncm_mset_peek (mset, mid);

      g_assert_cmpint (mid, >=, 0);

      if (mainmodel == NULL)
      {
        ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MAIN_MODEL_NOT_FOUND,
                                    "ncm_mset_load: cannot add submodel `%s', main model `%s' not found.",
                                    G_OBJECT_TYPE_NAME (submodel),
                                    ncm_mset_get_ns_by_id (mid));

        g_ptr_array_unref (submodel_array);
        g_key_file_unref (msetfile);
        g_strfreev (groups);

        return NULL;
      }
      else
      {
        ncm_model_add_submodel (mainmodel, submodel);
        _ncm_mset_set_pos_intern (mset, submodel, 0, error);
        NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset_load: ");
      }
    }

    g_key_file_unref (msetfile);
    g_strfreev (groups);

    if (valid_map)
      ncm_mset_prepare_fparam_map (mset);

    g_ptr_array_unref (submodel_array);

    return mset;
  }
}

/**
 * ncm_mset___getitem__:
 * @mset: a #NcmMSet
 * @model_id: a #GValue
 * @error: a #GError
 *
 * Gets the model with the id @model_id.
 * This method is used to implement the Python __getitem__ method.
 * The parameter @model_id is a #GValue which can be an integer
 * containing the model id or a string containing the model namespace.
 *
 * Returns: (transfer none): the model with the id @model_id.
 */
NcmModel *
ncm_mset___getitem__ (NcmMSet *mset, GValue *model_id, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  if (G_VALUE_HOLDS_INT (model_id))
  {
    NcmModelID mid  = g_value_get_int (model_id);
    NcmModel *model = ncm_mset_fetch (mset, mid, error);

    NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset___getitem__: ");

    return model;
  }
  else if (G_VALUE_HOLDS_STRING (model_id))
  {
    const gchar *ns = g_value_get_string (model_id);
    NcmModel *model = ncm_mset_fetch_by_name (mset, ns, error);

    NCM_UTIL_ON_ERROR_FORWARD (error, , NULL, "ncm_mset___getitem__: ");

    return model;
  }
  else
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_INVALID_ID,
                                "ncm_mset___getitem__: invalid argument type.");

    return NULL;
  }
}

/**
 * ncm_mset___setitem__:
 * @mset: a #NcmMSet
 * @model_id: a #GValue
 * @model: a #NcmModel
 * @error: a #GError
 *
 * Sets the model with the id @model_id to @model.
 * This method is used to implement the Python __setitem__ method.
 * The parameter @model_id is a #GValue which can be an integer
 * containing the model id or a string containing the model namespace.
 *
 */
void
ncm_mset___setitem__ (NcmMSet *mset, GValue *model_id, NcmModel *model, GError **error)
{
  NcmModelID mid = 0;
  guint stackpos = 0;

  g_return_if_fail (error == NULL || *error == NULL);

  if (G_VALUE_HOLDS_INT (model_id))
  {
    mid = g_value_get_int (model_id);
  }
  else if (G_VALUE_HOLDS_STRING (model_id))
  {
    const gchar *ns     = g_value_get_string (model_id);
    gchar **ns_stackpos = g_strsplit (ns, ":", 2);
    guint nelem         = g_strv_length (ns_stackpos);

    if (nelem == 1)
    {
      stackpos = 0;
    }
    else
    {
      gchar *endptr = NULL;

      stackpos = g_ascii_strtoll (ns_stackpos[1], &endptr, 10);

      if (*endptr != '\0')
      {
        ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_INVALID,
                                    "ncm_mset___setitem__: invalid namespace `%s'.", ns);
        g_strfreev (ns_stackpos);

        return;
      }
    }

    mid = ncm_mset_get_id_by_ns (ns_stackpos[0]);

    g_strfreev (ns_stackpos);

    if (mid < 0)
    {
      ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_NOT_FOUND,
                                  "ncm_mset___setitem__: namespace `%s' not found.", ns);

      return;
    }
  }
  else
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_INVALID_ID,
                                "ncm_mset___setitem__: invalid argument type.");

    return;
  }

  if (ncm_model_id (model) != mid)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_ID_MISMATCH,
                                "ncm_mset___setitem__: model id mismatch, expected %d, got %d.",
                                mid, ncm_model_id (model));

    return;
  }

  if (stackpos > 0)
    ncm_mset_set_pos (mset, model, stackpos, error);
  else
    ncm_mset_set (mset, model, error);

  NCM_UTIL_ON_ERROR_RETURN (error, , );
}


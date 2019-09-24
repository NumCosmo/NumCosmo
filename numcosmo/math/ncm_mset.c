/***************************************************************************
 *            ncm_mset.c
 *
 *  Fri May 25 09:37:54 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_obj_array.h"

enum
{
  PROP_0,
  PROP_VALID_MAP,
  PROP_MARRAY,
  PROP_FMAP,
  PROP_SIZE,
};

typedef struct _NcmMSetItem
{
  NcmModelID mid;
  NcmModel *model;
  gboolean dup;
  gint added_total_params;
} NcmMSetItem;

G_DEFINE_TYPE (NcmMSet, ncm_mset, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcmMSetPIndex, ncm_mset_pindex, ncm_mset_pindex_dup, ncm_mset_pindex_free);

static gint
_int_sort (gconstpointer a, gconstpointer b, gpointer user_data)
{
  NcmMSetItem **item_a = (NcmMSetItem **)a;
  NcmMSetItem **item_b = (NcmMSetItem **)b;

  return item_a[0]->mid - item_b[0]->mid;
}

NcmMSetItem *
_ncm_mset_item_new (NcmModel *model, NcmModelID mid)
{
  NcmMSetItem *item = g_slice_new0 (NcmMSetItem);
  item->model = ncm_model_ref (model);
  item->mid   = mid;
  item->dup   = FALSE;
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
  GError *error = NULL;

  mset->model_array     = g_ptr_array_sized_new (20);
  mset->model_item_hash = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, NULL);
  mset->mid_item_hash   = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, NULL);
  mset->fpi_hash        = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, (GDestroyNotify)g_array_unref);
  mset->fullname_parray = g_ptr_array_sized_new (NCM_MSET_INIT_MARRAY);
  mset->pi_array        = g_array_sized_new (FALSE, TRUE, sizeof (NcmMSetPIndex), 20);

  g_ptr_array_set_free_func (mset->model_array, (GDestroyNotify)_ncm_mset_item_free);
  g_ptr_array_set_free_func (mset->fullname_parray, g_free);

/* mset->fpi_array[i] = g_array_sized_new (FALSE, TRUE, sizeof (gint), 10); */

  mset->fullname_regex = g_regex_new ("^\\s*([A-Z][A-Za-z]+)\\:?([0-9]+)?\\:([0-9A-Z\\-a-z_]+)\\s*$", G_REGEX_OPTIMIZE, 0, &error);

  mset->valid_map = FALSE;
  mset->total_len = 0;
}

static void
_ncm_mset_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSet *mset = NCM_MSET (object);
  g_return_if_fail (NCM_IS_MSET (object));

  switch (prop_id)
  {
    case PROP_VALID_MAP:
      g_value_set_boolean (value, mset->valid_map);
      break;
    case PROP_MARRAY:
    {
      NcmObjArray *oa = ncm_obj_array_new ();
      guint i;
      for (i = 0; i < mset->model_array->len; i++)
      {
        NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
        if (!NCM_MODEL_GET_CLASS (item->model)->is_submodel)
          ncm_obj_array_add (oa, G_OBJECT (item->model));
      }
      g_value_take_boxed (value, oa);
      break;
    }
    case PROP_FMAP:
      g_value_take_boxed (value, ncm_mset_get_fmap (mset));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSet *mset = NCM_MSET (object);
  g_return_if_fail (NCM_IS_MSET (object));

  switch (prop_id)
  {
    case PROP_VALID_MAP:
    {
      const gboolean valid_map = g_value_get_boolean (value);
      if (valid_map)
      {
        if (!mset->valid_map)
          ncm_mset_prepare_fparam_map (mset);
      }
      else
        mset->valid_map = FALSE;
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
          ncm_mset_push (mset, model);
        }
      }
      break;
    }
    case PROP_FMAP:
    {
      const gchar * const *fmap = g_value_get_boxed (value);
      if (fmap != NULL)
        ncm_mset_set_fmap (mset, fmap, FALSE);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_dispose (GObject *object)
{
  NcmMSet *mset = NCM_MSET (object);

  g_clear_pointer (&mset->model_array, g_ptr_array_unref);

  g_clear_pointer (&mset->mid_item_hash, g_hash_table_unref);
  g_clear_pointer (&mset->model_item_hash, g_hash_table_unref);

  g_clear_pointer (&mset->fpi_hash, g_hash_table_unref);
  g_clear_pointer (&mset->fullname_parray, g_ptr_array_unref);
  g_clear_pointer (&mset->pi_array, g_array_unref);

  ncm_vector_clear (&mset->temp_fparams);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_parent_class)->dispose (object);
}

static void
_ncm_mset_finalize (GObject *object)
{
  NcmMSet *mset = NCM_MSET (object);

  g_regex_unref (mset->fullname_regex);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_parent_class)->finalize (object);
}

static void
ncm_mset_class_init (NcmMSetClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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

  klass->ns_table         = g_hash_table_new (g_str_hash, g_str_equal);
  klass->model_desc_array = g_array_new (FALSE, TRUE, sizeof (NcmMSetModelDesc));
}


/**
 * ncm_mset_pindex_new:
 * @mid: Model id
 * @pid: Parameter id
 *
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: (transfer full): FIXME
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
 * FIXME
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
 * FIXME
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
      g_error ("Cannot register model without a namespace.");
    if (desc == NULL)
      g_error ("Cannot register model without a description.");

    model_desc->ns   = g_strdup (ns);
    model_desc->desc = g_strdup (desc);

    if (long_desc != NULL)
      model_desc->long_desc = g_strdup (long_desc);
    else
      model_desc->long_desc = NULL;

    if (g_hash_table_lookup (mset_class->ns_table, ns) != NULL)
      g_error ("Model namespace <%s> already registered.", ns);

    g_hash_table_insert (mset_class->ns_table, model_desc->ns, GINT_TO_POINTER (model_class->model_id));

    G_UNLOCK (last_model_id);
    g_type_class_unref (mset_class);
  }
  else
  {
    g_error ("This model or its parent is already registred, id = %d. This function must be use once and only in the defining model.", model_class->model_id);
  }

  return;
}

/**
 * ncm_mset_empty_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSet *
ncm_mset_empty_new (void)
{
  return g_object_new (NCM_TYPE_MSET, NULL);
}

/**
 * ncm_mset_new:
 * @model0: a #NcmModel
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSet *
ncm_mset_new (gpointer model0, ...)
{
  NcmMSet *mset;
  va_list ap;

  va_start (ap, model0);
  mset = ncm_mset_newv (model0, ap);
  va_end (ap);

  return mset;
}

/**
 * ncm_mset_newv:
 * @model0: a #NcmModel
 * @ap: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSet *
ncm_mset_newv (gpointer model0, va_list ap)
{
  NcmMSet *mset = ncm_mset_empty_new ();
  NcmModel *model = NULL;

  ncm_mset_set (mset, model0);

  while ((model = va_arg (ap, NcmModel *)) != NULL)
  {
    g_assert (NCM_IS_MODEL (model));
    ncm_mset_push (mset, model);
  }

  return mset;
}

/**
 * ncm_mset_new_array:
 * @model_array: (array) (element-type NcmModel): a #GPtrArray of #NcmModel.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSet *
ncm_mset_new_array (GPtrArray *model_array)
{
  NcmMSet *mset = ncm_mset_empty_new ();
  guint i;

  for (i = 0; i < model_array->len; i++)
  {
    NcmModel *model = NCM_MODEL (g_ptr_array_index (model_array, i));
    g_assert (NCM_IS_MODEL (model));
    ncm_mset_push (mset, model);
  }

  return mset;
}

/**
 * ncm_mset_ref:
 * @mset: a #NcmMSet
 *
 * FIXME
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
 * FIXME
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
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcmMSet
 */
NcmMSet *
ncm_mset_shallow_copy (NcmMSet *mset)
{
  NcmMSet *mset_sc = ncm_mset_empty_new ();
  const guint nmodels = ncm_mset_nmodels (mset);
  guint i;

  for (i = 0; i < nmodels; i++)
  {
    NcmModel *model = ncm_mset_peek_array_pos (mset, i);
    if (!ncm_model_is_submodel (model))
      ncm_mset_push (mset_sc, model);
  }
  if (mset->valid_map)
    ncm_mset_prepare_fparam_map (mset_sc);

  return mset_sc;
}

/**
 * ncm_mset_free:
 * @mset: a #NcmMSet
 *
 * FIXME
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
 * FIXME
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
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmModel *
ncm_mset_peek (NcmMSet *mset, NcmModelID mid)
{
  NcmMSetItem *item = g_hash_table_lookup (mset->mid_item_hash, GINT_TO_POINTER (mid));
  if (item != NULL)
    return item->model;
  else
    return NULL;
}

/**
 * ncm_mset_peek_pos:
 * @mset: a #NcmMSet
 * @base_mid: a #NcmModelID
 * @stackpos_id: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
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
 * FIXME
 *
 * Returns: (transfer full): FIXME
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
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmModel *
ncm_mset_peek_array_pos (NcmMSet *mset, guint i)
{
  g_assert_cmpuint (i, <, mset->model_array->len);
  return ((NcmMSetItem *)g_ptr_array_index (mset->model_array, i))->model;
}

/**
 * ncm_mset_get_mid_array_pos:
 * @mset: a #NcmMSet
 * @i: array position
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmModelID
ncm_mset_get_mid_array_pos (NcmMSet *mset, guint i)
{
  g_assert_cmpuint (i, <, mset->model_array->len);
  return ((NcmMSetItem *)g_ptr_array_index (mset->model_array, i))->mid;
}

/**
 * ncm_mset_remove:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 *
 * FIXME
 *
 */
void
ncm_mset_remove (NcmMSet *mset, NcmModelID mid)
{
  NcmMSetItem *item = g_hash_table_lookup (mset->mid_item_hash, GINT_TO_POINTER (mid));
  if (item != NULL)
  {
    mset->total_len -= item->added_total_params;
    mset->valid_map = FALSE;

    g_hash_table_remove (mset->mid_item_hash, GINT_TO_POINTER (mid));
    if (!item->dup)
      g_hash_table_remove (mset->model_item_hash, item->model);

    g_assert (g_ptr_array_remove (mset->model_array, item));
  }
}


/**
 * ncm_mset_set:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 *
 * FIXME
 *
 */
void
ncm_mset_set (NcmMSet *mset, NcmModel *model)
{
  ncm_mset_set_pos (mset, model, 0);
}

/**
 * ncm_mset_push:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 *
 * FIXME
 *
 */
void
ncm_mset_push (NcmMSet *mset, NcmModel *model)
{
  g_assert (model != NULL);
  {
    NcmModelID base_mid  = ncm_model_id (model);
    guint stackpos_id = 0;
    while (TRUE)
    {
      NcmModelID mid = base_mid + stackpos_id;
      if (g_hash_table_lookup (mset->mid_item_hash, GINT_TO_POINTER (mid)) == NULL)
      {
        ncm_mset_set_pos (mset, model, stackpos_id);
        break;
      }
      stackpos_id++;
    }
  }
}

static void _ncm_mset_set_pos_intern (NcmMSet *mset, NcmModel *model, guint stackpos_id);

/**
 * ncm_mset_set_pos:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 * @stackpos_id: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_set_pos (NcmMSet *mset, NcmModel *model, guint stackpos_id)
{
  NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);
  
  if (model_class->is_submodel)
    g_error ("ncm_mset_set_pos: cannot add model `%s' directly to a NcmMSet.\n"
             "                  This model is a submodel of `%s' and it should"
             " be added to it instead of directly to a NcmMSet", 
             G_OBJECT_TYPE_NAME (model), ncm_mset_get_ns_by_id (model_class->main_model_id));
  else
    _ncm_mset_set_pos_intern (mset, model, stackpos_id);
}
  
static void
_ncm_mset_set_pos_intern (NcmMSet *mset, NcmModel *model, guint stackpos_id)
{
  g_assert (model != NULL);
  {
    NcmModelID mid  = ncm_model_id (model) + stackpos_id;
    guint model_len = ncm_model_len (model);
    NcmMSetItem *item = _ncm_mset_item_new (model, mid);
    GArray *fpi_array;
    NcmMSetItem *item0;
    guint i;

    if (stackpos_id > 0 && !(NCM_MODEL_GET_CLASS (model)->can_stack))
      g_error ("ncm_mset_set_pos: cannot stack object in position %u NcmMSet, type `%s' not allowed.", 
               stackpos_id, G_OBJECT_TYPE_NAME (model));

    ncm_mset_remove (mset, mid);

    if ((item0 = g_hash_table_lookup (mset->model_item_hash, model)) != NULL)
    {
      if (item->mid > item0->mid)
      {
        item->dup = TRUE;
        item->added_total_params = 0;
      }
      else
      {
        item0->dup = TRUE;
        item0->added_total_params = 0;
        g_hash_table_insert (mset->model_item_hash, model, item);
      }
      fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (item0->mid));
      g_assert (fpi_array != NULL);
      g_array_ref (fpi_array);
    }
    else
    {
      g_hash_table_insert (mset->model_item_hash, model, item);
      mset->total_len += item->added_total_params;
      fpi_array = g_array_new (FALSE, TRUE, sizeof (gint));
      g_array_set_size (fpi_array, model_len);
      for (i = 0; i < model_len; i++)
        g_array_index (fpi_array, gint, i) = -1;
    }

    g_hash_table_insert (mset->fpi_hash, GINT_TO_POINTER (mid), fpi_array);
    g_hash_table_insert (mset->mid_item_hash, GINT_TO_POINTER (mid), item);

    g_ptr_array_add (mset->model_array, item);
    g_ptr_array_sort_with_data (mset->model_array, _int_sort, NULL);

    if (mset->valid_map)
      ncm_mset_prepare_fparam_map (mset);

    for (i = 0; i < ncm_model_get_submodel_len (model); i++)
    {
      NcmModel *submodel = ncm_model_peek_submodel (model, i);
      _ncm_mset_set_pos_intern (mset, submodel, 0);
    }
  }
}

/**
 * ncm_mset_exists:
 * @mset: a #NcmMSet
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
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
 * @stackpos_id: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
 */
gint
ncm_mset_cmp_all (NcmMSet *mset0, NcmMSet *mset1)
{
  guint sub_nmodels0 = ncm_mset_nmodels (mset0);
  guint sub_nmodels1 = ncm_mset_nmodels (mset1);
  guint i;

  if (sub_nmodels0 != sub_nmodels1)
    return sub_nmodels0 < sub_nmodels1 ? -1 : 1;
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
 * @model_type: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmModelID
ncm_mset_get_id_by_type (GType model_type)
{
  g_assert (g_type_is_a (model_type, NCM_TYPE_MODEL));
  {
    NcmMSetClass *mset_class = g_type_class_ref (NCM_TYPE_MSET);
    const gchar *ns = g_type_name (model_type);
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
 * @ns: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
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
    GType t = g_type_from_name (ns);
    g_type_class_unref (mset_class);

    g_assert_cmpint (t, !=, 0);
    
    return t;
  }
}

/**
 * ncm_mset_prepare_fparam_map:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 */
void
ncm_mset_prepare_fparam_map (NcmMSet *mset)
{
  guint i;

  mset->fparam_len = 0;
  g_array_set_size (mset->pi_array, 0);
  g_ptr_array_set_size (mset->fullname_parray, 0);

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    
    if (item->dup)
      continue;
    else
    {
      GArray *fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (item->mid));
      guint pid;

      g_assert (fpi_array != NULL);

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        if (ncm_model_param_get_ftype (item->model, pid) == NCM_PARAM_TYPE_FREE)
        {
          NcmMSetPIndex pi = {item->mid, pid};
          g_array_append_val (mset->pi_array, pi);
          g_array_index (fpi_array, gint, pid) = mset->fparam_len;
          mset->fparam_len++;
        }
        else
          g_array_index (fpi_array, gint, pid) = -1;
      }
    }
  }

  g_ptr_array_set_size (mset->fullname_parray, mset->fparam_len);

  ncm_vector_clear (&mset->temp_fparams);
  if (mset->fparam_len > 0)
  {
    mset->temp_fparams = ncm_vector_new (mset->fparam_len);
  }
    
  mset->valid_map = TRUE;
}

/**
 * ncm_mset_set_fmap:
 * @mset: a #NcmMSet
 * @fmap: (in) (array zero-terminated=1) (element-type utf8): an array of strings
 * @update_models: a boolean
 *
 * FIXME
 *
 */
void
ncm_mset_set_fmap (NcmMSet *mset, const gchar *const *fmap, gboolean update_models)
{
  g_assert (fmap != NULL);
  {
    guint len = g_strv_length ((gchar **)fmap);
    guint i;
    for (i = 0; i < mset->model_array->len; i++)
    {
      NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
      if (item->dup)
        continue;
      else
      {
        GArray *fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (item->mid));
        guint pid;
        g_assert (fpi_array != NULL);
        for (pid = 0; pid < item->added_total_params; pid++)
          g_array_index (fpi_array, gint, pid) = -1;
      }
    }

    mset->fparam_len = len;
    g_array_set_size (mset->pi_array, 0);
    g_ptr_array_set_size (mset->fullname_parray, 0);

    for (i = 0; i < len; i++)
    {
      NcmMSetPIndex *pi = ncm_mset_param_get_by_full_name (mset, fmap[i]);
      if (pi == NULL)
      {
        g_error ("ncm_mset_set_fmap: cannot set fmap, invalid param `%s'.", fmap[i]);
      }
      else
      {
        GArray *fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (pi->mid));
        g_array_append_val (mset->pi_array, *pi);
        g_array_index (fpi_array, gint, pi->pid) = i;
      }
      ncm_mset_pindex_free (pi);
    }

    g_ptr_array_set_size (mset->fullname_parray, mset->fparam_len);
    mset->valid_map = TRUE;
    if (update_models)
      ncm_mset_param_set_ftype_from_fmap (mset);
  }
}


/**
 * ncm_mset_get_fmap:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: (transfer full) (array zero-terminated=1) (element-type utf8): an array of strings
 */
gchar **
ncm_mset_get_fmap (NcmMSet *mset)
{
  if (mset->valid_map)
  {
    guint len = ncm_mset_fparam_len (mset);
    gchar **fmap = g_new0 (gchar *, len + 1);
    guint i;
    for (i = 0; i < len; i++)
    {
      fmap[i] = g_strdup (ncm_mset_fparam_full_name (mset, i));
    }
    return fmap;
  }
  else
    return NULL;
}

/**
 * ncm_mset_total_len:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_total_len (NcmMSet *mset)
{
  return mset->total_len;
}

/**
 * ncm_mset_fparam_len:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_fparam_len (NcmMSet *mset)
{
  return mset->fparam_len;
}

/**
 * ncm_mset_max_param_name:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_max_param_name (NcmMSet *mset)
{
  guint name_size = 0;
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      guint pid;
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_max_model_nick (NcmMSet *mset)
{
  guint nick_size = 0;
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      guint pid;

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
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_nmodels (NcmMSet *mset)
{
  return mset->model_array->len;
}

/**
 * ncm_mset_pretty_log:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 */
void
ncm_mset_pretty_log (NcmMSet *mset)
{
  NcmMSetClass *mset_class = NCM_MSET_GET_CLASS (mset);
  guint name_size = ncm_mset_max_param_name (mset);
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      NcmModelID mid_base = NCM_MSET_GET_BASE_MID (item->mid);
      GArray *fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (item->mid));
      guint pid;

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
  guint name_size = GSL_MAX (ncm_mset_max_param_name (mset), 5);
  guint nick_size = ncm_mset_max_model_nick (mset);
  guint i;

  if (header != NULL)
    fprintf (out, "# %s\n ", header);
  else
    fprintf (out, "#\n");

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      const gchar *nick = ncm_model_nick (item->model);
      guint pid;

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
 * FIXME
 *
 */
void
ncm_mset_params_log_vals (NcmMSet *mset)
{
  guint i;

  g_message ("#");
  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      guint pid;
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
 * @out: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_params_print_vals (NcmMSet *mset, FILE *out)
{
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      guint pid;
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
 * FIXME
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
  for (i = 0; i < name_size; i++) g_message (" ");

  for (i = 0; i < free_params_len; i++)
    i ? g_message ("%s", box) : g_message ("-%s",box);
  if (i)
    g_message ("\n");

  for (i = 0; i < free_params_len; i++)
  {
    const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (mset, i);
    const gchar *pname = ncm_mset_fparam_name (mset, i);
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
  for (i = 0; i < name_size; i++) g_message (" ");
  for (i = 0; i < free_params_len; i++)
    i ? g_message ("%s", box) : g_message ("-%s",box);
  if (i)
    g_message ("\n");

  return;
}

/**
 * ncm_mset_params_valid:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_mset_params_valid (NcmMSet *mset)
{
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      if (!ncm_model_params_valid (item->model))
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_mset_params_bounds:
 * @mset: a #NcmMSet
 *
 * Check whenever the parameters respect the bounds.
 *
 * Returns: If TRUE the parameter respect the bounds.
 */
gboolean
ncm_mset_params_valid_bounds (NcmMSet *mset)
{
  guint i;
  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      if (!ncm_model_params_valid_bounds (item->model))
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_mset_cmp:
 * @mset0: a #NcmMSet
 * @mset1: a #NcmMSet
 * @cmp_model: whether to compare if the models correspond to the same objects
 *
 * Compares @mset0 and @mset1 and returns TRUE if both coitains the same models types.
 * If @cmp_model is TRUE compare also if the models correspond to the same objects.
 *
 * Returns: TRUE if @mset0 == @mset1.
 */
gboolean
ncm_mset_cmp (NcmMSet *mset0, NcmMSet *mset1, gboolean cmp_model)
{
  guint i;

  if (ncm_mset_total_len (mset0) != ncm_mset_total_len (mset0))
    return FALSE;

  if (mset0->model_array->len != mset1->model_array->len)
    return FALSE;

  for (i = 0; i < mset0->model_array->len; i++)
  {
    NcmMSetItem *item0 = g_ptr_array_index (mset0->model_array, i);
    NcmMSetItem *item1 = g_ptr_array_index (mset1->model_array, i);

    if (item0->mid != item1->mid)
      return FALSE;

    if (cmp_model && !ncm_model_is_equal (item0->model, item1->model))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_param_set:
 * @mset: a #NcmMSet
 * @mid: model id
 * @pid: parameter id
 * @x: FIXME
 *
 * FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 * Returns: FIXME
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
 * FIXME
 *
 */
void
ncm_mset_param_set_ftype (NcmMSet *mset, NcmModelID mid, guint pid, NcmParamType ftype)
{
  NcmModel *model = ncm_mset_peek (mset, mid);
  g_assert (model != NULL);
  ncm_model_param_set_ftype (model, pid, ftype);
  if (mset->valid_map)
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
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_model_param_set_ftype (item->model, pid, ftype);
      }
    }
  }
  if (mset->valid_map)
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
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup || (item->mid != mid))
      continue;
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_model_param_set_ftype (item->model, pid, ftype);
      }
    }
  }
  if (mset->valid_map)
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
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup || (item->mid == mid))
      continue;
    else
    {
      gint pid;

      for (pid = 0; pid < item->added_total_params; pid++)
      {
        ncm_model_param_set_ftype (item->model, pid, ftype);
      }
    }
  }
  if (mset->valid_map)
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
  guint i;

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup)
      continue;
    else
    {
      GArray *fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (item->mid));

      gint pid;
      for (pid = 0; pid < item->added_total_params; pid++)
      {
        if (g_array_index (fpi_array, gint, pid) == - 1)
          continue;
        ncm_model_param_set_ftype (item->model, pid, NCM_PARAM_TYPE_FREE);
      }
    }
  }
  if (mset->valid_map)
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
  guint i;
  guint j = 0;

  g_assert_cmpuint (ncm_vector_len (params), ==, ncm_mset_total_len (mset));

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup || (item->added_total_params == 0))
      continue;
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
  guint j = 0;
  gint i;

  g_assert (ncm_vector_len (params) == ncm_mset_total_len (mset));

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item = g_ptr_array_index (mset->model_array, i);
    if (item->dup || (item->added_total_params == 0))
      continue;
    else
    {
      guint pid;
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
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmParamType
ncm_mset_param_get_ftype (NcmMSet *mset, NcmModelID mid, guint pid)
{
  NcmModel *model = ncm_mset_peek (mset, mid);
  if (model == NULL)
    g_error ("ncm_mset_param_get_ftype: cannot get ftype of mode %d, model not set.", mid);
  return ncm_model_param_get_ftype (model, pid);
}

/**
 * ncm_mset_param_set_pi:
 * @mset: a #NcmMSet
 * @pi: a #NcmMSetPIndex
 * @x: FIXME
 * @n: FIXME
 *
 * FIXME
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
 * @pi: a #NcmMSetPIndex
 * @x: FIXME
 * @n: FIXME
 *
 * FIXME
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
 * FIXME
 *
 */
void
ncm_mset_fparams_get_vector (NcmMSet *mset, NcmVector *x)
{
  guint fpi;

  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_vector_set (x, fpi, ncm_mset_param_get (mset, pi.mid, pi.pid));
  }
}

/**
 * ncm_mset_fparams_get_vector_offset:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 * @offset: starting index
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_get_vector_offset (NcmMSet *mset, NcmVector *x, guint offset)
{
  guint fpi;

  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_vector_set (x, fpi + offset, ncm_mset_param_get (mset, pi.mid, pi.pid));
  }
}

/**
 * ncm_mset_fparams_set_vector:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_vector (NcmMSet *mset, const NcmVector *x)
{
  guint fpi;

  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, ncm_vector_get (x, fpi));
  }
}

/**
 * ncm_mset_fparams_set_vector_offset:
 * @mset: a #NcmMSet
 * @x: a #NcmVector
 * @offset: starting index
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_vector_offset (NcmMSet *mset, const NcmVector *x, guint offset)
{
  guint fpi;

  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, ncm_vector_get (x, fpi + offset));
  }
}

/**
 * ncm_mset_fparams_set_array:
 * @mset: a #NcmMSet
 * @x: (array) (element-type double): FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_array (NcmMSet *mset, const gdouble *x)
{
  guint fpi;
  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, x[fpi]);
  }
}

/**
 * ncm_mset_fparams_set_gsl_vector: (skip)
 * @mset: a #NcmMSet.
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparams_set_gsl_vector (NcmMSet *mset, const gsl_vector *x)
{
  guint fpi;
  for (fpi = 0; fpi < mset->fparam_len; fpi++)
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, fpi);
    ncm_mset_param_set (mset, pi.mid, pi.pid, gsl_vector_get (x, fpi));
  }
}

/**
 * ncm_mset_fparams_len:
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_mset_fparams_len (NcmMSet *mset)
{
  g_assert (mset->valid_map);
  return mset->fparam_len;
}

/**
 * ncm_mset_fparam_name:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const gchar *
ncm_mset_fparam_name (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_name (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_symbol:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const gchar *
ncm_mset_fparam_symbol (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_symbol (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_full_name:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const gchar *
ncm_mset_fparam_full_name (NcmMSet *mset, guint n)
{
  gchar *fullname;
  g_assert (mset->valid_map && n < mset->fparam_len);
  if ((fullname = g_ptr_array_index (mset->fullname_parray, n)) != NULL)
    return fullname;
  else
  {
    NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    const gchar *model_ns = ncm_mset_get_ns_by_id (pi.mid); /* ncm_model_nick (ncm_mset_peek (mset, pi.mid));*/
    const gchar *pname = ncm_mset_param_name (mset, pi.mid, pi.pid);
    const guint stackpos_id = pi.mid % NCM_MSET_MAX_STACKSIZE;

    if (stackpos_id > 0)
      fullname = g_strdup_printf ("%s:%02u:%s", model_ns, stackpos_id, pname);
    else
      fullname = g_strdup_printf ("%s:%s", model_ns, pname);

    g_ptr_array_index (mset->fullname_parray, n) = fullname;
  }
  return fullname;
}

/**
 * ncm_mset_param_get_by_full_name:
 * @mset: a #NcmMSet
 * @fullname: param's full name
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetPIndex *
ncm_mset_param_get_by_full_name (NcmMSet *mset, const gchar *fullname)
{
  GMatchInfo *match_info = NULL;
  NcmMSetPIndex *pi = NULL;

  if (g_regex_match (mset->fullname_regex, fullname, 0, &match_info))
  {
    gint nm = g_match_info_get_match_count (match_info);
    gchar *ns         = NULL;
    gchar *stackpos_s = NULL;
    gchar *pid_s      = NULL;
    gchar *endptr     = NULL;
    guint pid      = 0;
    guint stackpos = 0;
    NcmModelID mid;
    NcmModel *model;

    g_assert_cmpint (nm, ==, 4);

    ns         = g_match_info_fetch (match_info, 1);
    stackpos_s = g_match_info_fetch (match_info, 2);
    pid_s      = g_match_info_fetch (match_info, 3);

    mid = ncm_mset_get_id_by_ns (ns);
    if (mid < 0)
      g_error ("ncm_mset_param_get_by_full_name: namespace `%s' not found.", ns);

    if (*stackpos_s != '\0')
    {
      stackpos = g_ascii_strtoll (stackpos_s, &endptr, 10);
      if (*endptr != '\0')
        g_error ("ncm_mset_param_get_by_full_name: invalid stackpos number `%s'.", stackpos_s);
    }
    g_free (stackpos_s);

    mid += stackpos;
    model = ncm_mset_peek (mset, mid);
    if (model != NULL)
    {
      endptr = NULL;
      pid = g_ascii_strtoll (pid_s, &endptr, 10);

      if (*endptr != '\0')
      {
        if (ncm_model_param_index_from_name (model, pid_s, &pid))
          pi = ncm_mset_pindex_new (mid, pid);
      }
      else if (pid < ncm_model_len (model))
      {
        pi = ncm_mset_pindex_new (mid, pid);
      }
    }
    g_free (pid_s);
    g_free (ns);
  }
  else
    g_error ("ncm_mset_param_get_by_full_name: invalid full name `%s'.",  fullname);

  g_match_info_free (match_info);

  return pi;
}

/**
 * ncm_mset_fparam_get_scale:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_scale (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_scale (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_lower_bound:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_lower_bound (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_lower_bound (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_upper_bound:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_upper_bound (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_upper_bound (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_get_abstol:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get_abstol (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get_abstol (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_set_scale:
 * @mset: a #NcmMSet
 * @n: free parameter index
 * @scale: new scale
 *
 * FIXME
 *
 */
void
ncm_mset_fparam_set_scale (NcmMSet *mset, guint n, gdouble scale)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    ncm_mset_param_set_scale (mset, pi.mid, pi.pid, scale);
  }
}

/**
 * ncm_mset_fparam_valid_bounds:
 * @mset: a #NcmMSet
 * @theta: free parameters vector
 *
 * FIXME
 *
 * Returns: whether @theta contain values respecting the parameter bounds.
 */
gboolean
ncm_mset_fparam_valid_bounds (NcmMSet *mset, NcmVector *theta)
{
  g_assert (mset->valid_map);
  g_assert_cmpuint (ncm_vector_len (theta), ==, mset->fparam_len);
  {
    guint i;
    for (i = 0; i < mset->fparam_len; i++)
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
 * FIXME
 *
 * Returns: whether @theta contain values respecting the parameter bounds.
 */
gboolean
ncm_mset_fparam_valid_bounds_offset (NcmMSet *mset, NcmVector *theta, guint offset)
{
  g_assert (mset->valid_map);
  g_assert_cmpuint (ncm_vector_len (theta), ==, mset->fparam_len + offset);
  {
    guint i;
    for (i = 0; i < mset->fparam_len; i++)
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
 * FIXME
 *
 * Returns: whether @theta contain values respecting all requirements.
 */
gboolean
ncm_mset_fparam_validate_all (NcmMSet *mset, NcmVector *theta)
{
  gboolean valid;
  g_assert (mset->valid_map);
  g_assert_cmpuint (ncm_vector_len (theta), ==, mset->fparam_len);

  ncm_mset_fparams_get_vector (mset, mset->temp_fparams);
  ncm_mset_fparams_set_vector (mset, theta);

  valid = ncm_mset_params_valid (mset) && ncm_mset_params_valid_bounds (mset);
  ncm_mset_fparams_set_vector (mset, mset->temp_fparams);
  
  return valid;
}

/**
 * ncm_mset_fparam_get:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_fparam_get (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_get (mset, pi.mid, pi.pid);
  }
}

/**
 * ncm_mset_fparam_set:
 * @mset: a #NcmMSet
 * @n: free parameter index
 * @x: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_fparam_set (NcmMSet *mset, guint n, const gdouble x)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  {
    const NcmMSetPIndex pi = g_array_index (mset->pi_array, NcmMSetPIndex, n);
    return ncm_mset_param_set (mset, pi.mid, pi.pid, x);
  }
}

/**
 * ncm_mset_fparam_get_pi:
 * @mset: a #NcmMSet
 * @n: free parameter index
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const NcmMSetPIndex *
ncm_mset_fparam_get_pi (NcmMSet *mset, guint n)
{
  g_assert (mset->valid_map && n < mset->fparam_len);
  return &(g_array_index (mset->pi_array, NcmMSetPIndex, n));
}

/**
 * ncm_mset_fparam_get_fpi:
 * @mset: a #NcmMSet
 * @mid: a #NcmModelID
 * @pid: parameter id
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint
ncm_mset_fparam_get_fpi (NcmMSet *mset, NcmModelID mid, guint pid)
{
  g_assert (mset->valid_map);

  {
    GArray *fpi_array = g_hash_table_lookup (mset->fpi_hash, GINT_TO_POINTER (mid));
    g_assert (fpi_array != NULL);
    return g_array_index (fpi_array, gint, pid);
  }
}

/**
 * ncm_mset_fparam_get_pi_by_name:
 * @mset: a #NcmMSet
 * @name: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const NcmMSetPIndex *
ncm_mset_fparam_get_pi_by_name (NcmMSet *mset, const gchar *name)
{
  guint match = 0;
  guint match_i = 0;
  guint i;
  g_assert (mset->valid_map);
  for (i = 0; i < mset->fparam_len; i++)
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
    g_warning ("ncm_mset_fparam_get_pi_by_name: more than one [%u] parameters with the same name %s, use the full name to avoid ambiguities, returning the last match.",
               match, name);
    return NULL;
  }
  else if (match == 1)
  {
    return ncm_mset_fparam_get_pi (mset, match_i);
  }
  else
  {
    for (i = 0; i < mset->fparam_len; i++)
    {
      const gchar *name_i = ncm_mset_fparam_full_name (mset, i);
      if (strcmp (name, name_i) == 0)
      {
        match++;
        match_i = i;
      }
    }
    if (match == 0)
      return NULL;
    else if (match > 1)
    {
      g_error ("ncm_mset_fparam_get_pi_by_name: more than one full names [%u] match %s.", match, name);
      return NULL;
    }
    else
      return ncm_mset_fparam_get_pi (mset, match_i);
  }
}

/**
 * ncm_mset_save:
 * @mset: a #NcmMSet
 * @ser: a #NcmSerialize
 * @filename: FIXME
 * @save_comment: FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_save (NcmMSet *mset, NcmSerialize *ser, const gchar *filename, gboolean save_comment)
{
  NcmMSetClass *mset_class = NCM_MSET_GET_CLASS (mset);
  GKeyFile *msetfile = g_key_file_new ();
  guint i;

  {
    GError *error = NULL;
    gchar *mset_desc = ncm_cfg_string_to_comment ("NcmMSet");
    gchar *mset_valid_map_desc = ncm_cfg_string_to_comment ("valid-map property");

    g_key_file_set_boolean (msetfile, "NcmMSet", "valid_map", mset->valid_map);

    if (save_comment)
    {
      if (!g_key_file_set_comment (msetfile, "NcmMSet", NULL, mset_desc, &error))
        g_error ("ncm_mset_save: %s", error->message);
      if (!g_key_file_set_comment (msetfile, "NcmMSet", "valid_map", mset_valid_map_desc, &error))
        g_error ("ncm_mset_save: %s", error->message);
    }

    g_free (mset_desc);
    g_free (mset_valid_map_desc);
  }

  for (i = 0; i < mset->model_array->len; i++)
  {
    NcmMSetItem *item       = g_ptr_array_index (mset->model_array, i);
    NcmModelID mid_base     = item->mid / NCM_MSET_MAX_STACKSIZE;
    const guint stackpos_id = item->mid % NCM_MSET_MAX_STACKSIZE;
    const gchar *ns = g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).ns;
    GError *error = NULL;
    gchar *group;

    if (stackpos_id > 0)
      group = g_strdup_printf ("%s:%02u", ns, stackpos_id);
    else
      group = g_strdup_printf ("%s", ns);

    if (item->dup)
    {
      NcmMSetItem *item0 = g_hash_table_lookup (mset->model_item_hash, item->model);
      g_assert (item0 != NULL);

      g_key_file_set_uint64 (msetfile, group, ns, item0->mid % NCM_MSET_MAX_STACKSIZE);
      if (save_comment)
      {
        gchar *model_desc = ncm_cfg_string_to_comment (g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).desc);
        if (!g_key_file_set_comment (msetfile, group, NULL, model_desc, &error))
          g_error ("ncm_mset_save: %s", error->message);
        g_free (model_desc);
      }
    }
    else
    {
      NcmModelClass *model_class = NCM_MODEL_GET_CLASS (item->model);
      GObjectClass *oclass       = G_OBJECT_CLASS (model_class);
      GVariant *model_var        = ncm_serialize_to_variant (ser, G_OBJECT (item->model));
      guint nsubmodels           = ncm_model_get_submodel_len (item->model);
      GVariant *params = NULL;
      gchar *obj_name = NULL;
      guint nparams, j;

      g_variant_get (model_var, "{s@a{sv}}", &obj_name, &params);
      nparams = g_variant_n_children (params);
      g_key_file_set_value (msetfile, group, ns, obj_name);

      if (save_comment)
      {
        gchar *model_desc = ncm_cfg_string_to_comment (g_array_index (mset_class->model_desc_array, NcmMSetModelDesc, mid_base).desc);
        if (!g_key_file_set_comment (msetfile, group, NULL, model_desc, &error))
          g_error ("ncm_mset_save: %s", error->message);
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
            gchar *param_str = g_variant_print (value, TRUE);
            if (param_spec == NULL)
              g_error ("ncm_mset_save: property `%s' not found in object `%s'.", key, obj_name);

            g_key_file_set_value (msetfile, group, key, param_str);
            if (save_comment)
            {
              const gchar *blurb = g_param_spec_get_blurb (param_spec);
              if (blurb != NULL && blurb[0] != 0)
              {
                gchar *desc = ncm_cfg_string_to_comment (blurb);
                if (!g_key_file_set_comment (msetfile, group, key, desc, &error))
                  g_error ("ncm_mset_save: %s", error->message);
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
    GError *error = NULL;
    gsize len = 0;
    gchar *mset_data = g_key_file_to_data (msetfile, &len, &error);
    if (error != NULL)
      g_error ("Error converting NcmMSet to configuration file: %s", error->message);
    if (!g_file_set_contents (filename, mset_data, len, &error))
      g_error ("Error saving configuration file to disk: %s", error->message);
    g_free (mset_data);
    g_key_file_free (msetfile);
  }
}


/**
 * ncm_mset_load:
 * @filename: mset filename
 * @ser: a #NcmSerialize
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSet *
ncm_mset_load (const gchar *filename, NcmSerialize *ser)
{
  NcmMSet *mset = ncm_mset_empty_new ();
  GKeyFile *msetfile = g_key_file_new ();
  GError *error = NULL;
  gchar **groups = NULL;
  gsize ngroups = 0;
  gboolean valid_map = FALSE;
  GPtrArray *submodel_array = g_ptr_array_new ();
  guint i;

  g_ptr_array_set_free_func (submodel_array, (GDestroyNotify) ncm_model_free);

  if (!g_key_file_load_from_file (msetfile, filename, G_KEY_FILE_KEEP_COMMENTS | G_KEY_FILE_KEEP_TRANSLATIONS, &error))
  {
    g_error ("ncm_mset_load: Invalid mset configuration file: %s %s", filename, error->message);
    return NULL;
  }

  if (g_key_file_has_group (msetfile, "NcmMSet"))
  {
    if (g_key_file_has_key (msetfile, "NcmMSet", "valid_map", &error))
    {
      if (error != NULL)
        g_error ("ncm_mset_load: %s", error->message);
      valid_map = g_key_file_get_boolean (msetfile, "NcmMSet", "valid_map", &error);
      if (error != NULL)
        g_error ("ncm_mset_load: %s", error->message);
    }
    g_key_file_remove_group (msetfile, "NcmMSet", &error);
    if (error != NULL)
      g_error ("ncm_mset_load: %s", error->message);
  }

  groups = g_key_file_get_groups (msetfile, &ngroups);
  for (i = 0; i < ngroups; i++)
  {
    GString *obj_ser = g_string_sized_new (200);
    gchar *ns = g_strdup (groups[i]);
    gchar *twopoints = g_strrstr (ns, ":");
    guint stackpos = 0;

    if (twopoints != NULL)
    {
      *twopoints = '\0';
      stackpos = atoi (++twopoints);
    }

    if (!g_key_file_has_key (msetfile, groups[i], ns, &error))
    {
      if (error != NULL)
        g_error ("ncm_mset_load: %s", error->message);
      g_error ("ncm_mset_load: Every group must contain a key with same name indicating the object type `%s' `%s'.", groups[i], ns);
    }

    {
      gchar *obj_type = g_key_file_get_value (msetfile, groups[i], ns, &error);

      if (strlen (obj_type) < 5)
      {
        gchar *endptr = NULL;
        guint stackpos_orig = g_ascii_strtoll (obj_type, &endptr, 10);

        g_assert (endptr != NULL);

        if (*endptr == '\0')
        {
          NcmModelID id = ncm_mset_get_id_by_ns (ns);
          g_assert_cmpint (id, >, -1);
          {
            NcmModel *model = ncm_mset_peek_pos (mset, id, stackpos_orig);
            g_assert (model != NULL);
            ncm_mset_set_pos (mset, model, stackpos);
          }
          continue;
        }
      }

      g_string_append_printf (obj_ser, "{\'%s\', @a{sv} {", obj_type);
      g_free (obj_type);
      g_key_file_remove_key (msetfile, groups[i], ns, &error);
      if (error != NULL)
        g_error ("ncm_mset_load: %s", error->message);
    }

    {
      gsize nkeys = 0;
      gchar **keys = g_key_file_get_keys (msetfile, groups[i], &nkeys, &error);
      guint j;
      if (error != NULL)
        g_error ("ncm_mset_load: %s", error->message);
      for (j = 0; j < nkeys; j++)
      {
        if (strcmp (keys[j], "submodel-array") == 0)
        {
          g_error ("ncm_mset_load: serialized version of mset models cannot contain the submodel-array property.");
        }
        else
        {
          gchar *propval = g_key_file_get_value (msetfile, groups[i], keys[j], &error);
          if (error != NULL)
            g_error ("ncm_mset_load: %s", error->message);
          g_string_append_printf (obj_ser, "\'%s\':<%s>", keys[j], propval);
          g_free (propval);
          if (j + 1 != nkeys)
            g_string_append (obj_ser, ", ");
        }
      }
      g_string_append (obj_ser, "}}");
      g_strfreev (keys);
    }

    {
      GObject *obj = ncm_serialize_from_string (ser, obj_ser->str);
      g_assert (NCM_IS_MODEL (obj));
      if (ncm_mset_exists_pos (mset, NCM_MODEL (obj), stackpos))
        g_error ("ncm_mset_load: a model ``%s'' already exists in NcmMSet.", obj_ser->str);
      if (ncm_model_is_submodel (NCM_MODEL (obj)))
      {
        NcmModel *submodel = ncm_model_ref (NCM_MODEL (obj));
        g_ptr_array_add (submodel_array, submodel);
      }
      else
      {
        ncm_mset_set_pos (mset, NCM_MODEL (obj), stackpos);
      }
      ncm_model_free (NCM_MODEL (obj));
    }
    g_string_free (obj_ser, TRUE);
    g_free (ns);
  }

  for (i = 0; i < submodel_array->len; i++)
  {
    NcmModel *submodel = g_ptr_array_index (submodel_array, i);
    NcmModelID mid = ncm_model_main_model (submodel);
    NcmModel *mainmodel = ncm_mset_peek (mset, mid);

    g_assert_cmpint (mid, >=, 0);
    if (mainmodel == NULL)
    {
      g_error ("ncm_mset_load: cannot add submodel `%s', main model `%s' not found.", 
               G_OBJECT_TYPE_NAME (submodel),
               ncm_mset_get_ns_by_id (mid));
    }
    else
    {
      ncm_model_add_submodel (mainmodel, submodel);
      _ncm_mset_set_pos_intern (mset, submodel, 0);
    }
  }
  
  g_key_file_unref (msetfile);
  g_strfreev (groups);
    
  if (valid_map)
    ncm_mset_prepare_fparam_map (mset);

  g_ptr_array_unref (submodel_array);
  return mset;
}

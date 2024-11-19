/***************************************************************************
 *            ncm_model.c
 *
 *  Fri February 24 21:18:21 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * SECTION:ncm_model
 * @title: NcmModel
 * @short_description: Abstract class for implementing models.
 *
 * The #NcmModel abstract class represents a general model. This object serves
 * for two general objectives. First, all the numerical properties (doubles), i.e.,
 * parameters, are implemented by the class functions described below, this
 * allows the implementation of a general statistical analyses based on these
 * models. Second, each child of NcmModel can register itself as a model type.
 * This allows multiples models types to be used simultaneously.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model.h"
#include "math/ncm_mset.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_obj_array.h"

typedef struct _NcmModelPrivate
{
  /*< private >*/
  GObject parent_instance;
  NcmReparam *reparam;
  NcmObjArray *sparams;
  NcmVector *params;
  NcmVector *p;
  GArray *sparam_modified;
  GArray *vparam_pos;
  GArray *vparam_len;
  GArray *ptypes;
  GHashTable *sparams_name_id;
  GHashTable *submodel_mid_pos;
  GPtrArray *submodel_array;
  guint total_len;
  guint64 pkey;
  guint64 skey;
  guint64 slkey[NCM_MODEL_MAX_STATES];
  gdouble *params_ptr;
} NcmModelPrivate;

enum
{
  PROP_0,
  PROP_NAME,
  PROP_NICK,
  PROP_SPARAMS_LEN,
  PROP_VPARAMS_LEN,
  PROP_IMPLEMENTATION,
  PROP_SPARAM_ARRAY,
  PROP_PTYPES,
  PROP_REPARAM,
  PROP_SUBMODEL_ARRAY,
};

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-model-error, ncm_model_error) 
/* *INDENT-ON* */

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmModel, ncm_model, G_TYPE_OBJECT)

static void
ncm_model_init (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gint i;

  self->sparams         = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_sparam_free);
  self->sparams_name_id = g_hash_table_new_full (&g_str_hash, &g_str_equal, &g_free, NULL);
  self->sparam_modified = g_array_new (FALSE, TRUE, sizeof (gboolean));
  self->params          = NULL;
  self->params_ptr      = NULL;
  self->p               = NULL;
  self->vparam_len      = g_array_sized_new (TRUE, TRUE, sizeof (guint), 0);
  self->vparam_pos      = g_array_sized_new (TRUE, TRUE, sizeof (guint), 0);
  self->pkey            = 1;
  self->skey            = 0;

  for (i = 0; i < NCM_MODEL_MAX_STATES; i++)
    self->slkey[i] = 0;

  self->reparam          = NULL;
  self->ptypes           = g_array_new (FALSE, TRUE, sizeof (NcmParamType));
  self->submodel_array   = g_ptr_array_new ();
  self->submodel_mid_pos = g_hash_table_new_full (g_direct_hash, g_direct_equal, NULL, NULL);

  g_ptr_array_set_free_func (self->submodel_array, (GDestroyNotify) ncm_model_free);
}

static void
_ncm_model_set_sparams (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  NcmModelClass *model_class   = NCM_MODEL_GET_CLASS (model);
  guint i;

  g_hash_table_remove_all (self->sparams_name_id);
  g_ptr_array_set_size (self->sparams, 0);
  g_ptr_array_set_size (self->sparams, self->total_len);
  g_array_set_size (self->sparam_modified, self->total_len);

  for (i = 0; i < model_class->sparam_len; i++)
  {
    NcmSParam *sp = g_ptr_array_index (model_class->sparam, i);

    g_array_index (self->ptypes, NcmParamType, i)      = NCM_PARAM_TYPE_FIXED;
    g_ptr_array_index (self->sparams, i)               = ncm_sparam_copy (sp);
    g_array_index (self->sparam_modified, gboolean, i) = FALSE;
    g_hash_table_insert (self->sparams_name_id, g_strdup (ncm_sparam_name (sp)), GUINT_TO_POINTER (i));
  }

  for (i = 0; i < model_class->vparam_len; i++)
  {
    const guint len = g_array_index (self->vparam_len, guint, i);
    const guint pos = g_array_index (self->vparam_pos, guint, i);
    NcmVParam *vp   = ncm_vparam_copy (g_ptr_array_index (model_class->vparam, i));
    guint j;

    ncm_vparam_set_len (vp, len);

    for (j = 0; j < len; j++)
    {
      const guint n = pos + j;
      NcmSParam *sp = ncm_vparam_peek_sparam (vp, j);

      g_array_index (self->ptypes, NcmParamType, n)      = NCM_PARAM_TYPE_FIXED;
      g_ptr_array_index (self->sparams, n)               = ncm_sparam_ref (sp);
      g_array_index (self->sparam_modified, gboolean, n) = FALSE;
      g_hash_table_insert (self->sparams_name_id, g_strdup (ncm_sparam_name (sp)), GUINT_TO_POINTER (n));
    }

    ncm_vparam_free (vp);
  }
}

static void
_ncm_model_set_sparams_from_dict (NcmModel *model, NcmObjDictInt *modified_sparams)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  if (modified_sparams != NULL)
  {
    GHashTableIter iter;
    gint *key;
    NcmSParam *value;

    g_hash_table_iter_init (&iter, modified_sparams);

    while (g_hash_table_iter_next (&iter, (gpointer *) &key, (gpointer *) &value))
    {
      const guint n         = *key;
      NcmSParam *current_sp = g_ptr_array_index (self->sparams, n);

      if (n >= self->total_len)
        g_error ("_ncm_model_set_sparams_from_dict: parameter %u is out of range (0-%u)", n, self->total_len - 1);

      g_assert_nonnull (current_sp);

      g_hash_table_remove (self->sparams_name_id, ncm_sparam_name (current_sp));
      ncm_sparam_clear ((NcmSParam **) &g_ptr_array_index (self->sparams, n));

      g_array_index (self->sparam_modified, gboolean, n) = TRUE;
      g_ptr_array_index (self->sparams, n)               = ncm_sparam_copy (value);
      g_hash_table_insert (self->sparams_name_id, g_strdup (ncm_sparam_name (value)), GUINT_TO_POINTER (n));
    }
  }
}

static void
_ncm_model_sparams_remove_reparam (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  if (self->reparam != NULL)
  {
    ncm_reparam_clear (&self->reparam);
    self->p = ncm_vector_ref (self->params);
  }
}

static void
_ncm_model_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_model_parent_class)->constructed (object);
  {
    NcmModel *model              = NCM_MODEL (object);
    NcmModelPrivate * const self = ncm_model_get_instance_private (model);
    NcmModelClass *model_class   = NCM_MODEL_GET_CLASS (model);
    guint i;

    self->total_len = model_class->sparam_len;

    for (i = 0; i < model_class->vparam_len; i++)
    {
      g_array_index (self->vparam_pos, guint, i) = self->total_len;
      self->total_len                           += g_array_index (self->vparam_len, guint, i);
    }

    self->params     = ncm_vector_new (self->total_len == 0 ? 1 : self->total_len);
    self->params_ptr = ncm_vector_data (self->params);
    self->p          = ncm_vector_ref (self->params);
    g_array_set_size (self->ptypes, self->total_len);
    _ncm_model_set_sparams (model);
    ncm_model_params_set_default (model);
  }
}

static void
_ncm_model_dispose (GObject *object)
{
  NcmModel *model              = NCM_MODEL (object);
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_clear (&self->params);
  ncm_vector_clear (&self->p);

  ncm_reparam_clear (&self->reparam);

  g_clear_pointer (&self->vparam_len,       g_array_unref);
  g_clear_pointer (&self->vparam_pos,       g_array_unref);
  g_clear_pointer (&self->ptypes,           g_array_unref);
  g_clear_pointer (&self->sparams,          g_ptr_array_unref);
  g_clear_pointer (&self->sparams_name_id,  g_hash_table_unref);
  g_clear_pointer (&self->sparam_modified,  g_array_unref);

  g_clear_pointer (&self->submodel_array,   g_ptr_array_unref);
  g_clear_pointer (&self->submodel_mid_pos, g_hash_table_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_parent_class)->dispose (object);
}

static void
_ncm_model_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_parent_class)->finalize (object);
}

static NcmSParam *_ncm_model_param_peek_desc (NcmModel *model, guint n, gboolean *is_original);
static NcmSParam *_ncm_model_orig_param_peek_desc (NcmModel *model, guint n);

static void
_ncm_model_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModel *model = NCM_MODEL (object);

  g_return_if_fail (NCM_IS_MODEL (object));

  switch (prop_id)
  {
    case PROP_SPARAM_ARRAY:
      _ncm_model_set_sparams_from_dict (model, g_value_get_boxed (value));
      break;
    case PROP_REPARAM:
      ncm_model_set_reparam (model, g_value_get_object (value), NULL);
      break;
    case PROP_SUBMODEL_ARRAY:
    {
      NcmObjArray *oa = (NcmObjArray *) g_value_get_boxed (value);

      if (oa != NULL)
      {
        guint i;

        for (i = 0; i < oa->len; i++)
        {
          NcmModel *submodel = NCM_MODEL (ncm_obj_array_peek (oa, i));

          if (!NCM_MODEL_GET_CLASS (submodel)->is_submodel)
            g_error ("_ncm_model_set_property: NcmModel submodel array can only contain submodels `%s'.",
                     G_OBJECT_TYPE_NAME (submodel));

          ncm_model_add_submodel (model, submodel);
        }
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_model_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModel *model              = NCM_MODEL (object);
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  NcmModelClass *model_class   = NCM_MODEL_GET_CLASS (object);

  g_return_if_fail (NCM_IS_MODEL (object));

  switch (prop_id)
  {
    case PROP_NAME:
      g_value_set_string (value, model_class->name);
      break;
    case PROP_NICK:
      g_value_set_string (value, model_class->nick);
      break;
    case PROP_SPARAMS_LEN:
      g_value_set_uint (value, model_class->sparam_len);
      break;
    case PROP_VPARAMS_LEN:
      g_value_set_uint (value, model_class->vparam_len);
      break;
    case PROP_IMPLEMENTATION:
      g_value_set_uint64 (value, model_class->impl_flag);
      break;
    case PROP_SPARAM_ARRAY:
    {
      NcmObjDictInt *modified_sparams = ncm_obj_dict_int_new ();
      guint i;

      for (i = 0; i < self->sparam_modified->len; i++)
      {
        if (g_array_index (self->sparam_modified, gboolean, i))
          ncm_obj_dict_int_add (modified_sparams, i,
                                G_OBJECT (_ncm_model_orig_param_peek_desc (model, i)));
      }

      g_value_take_boxed (value, modified_sparams);
      break;
    }
    case PROP_REPARAM:
      g_value_set_object (value, self->reparam);
      break;
    case PROP_PTYPES:
      g_value_set_boxed (value, self->ptypes);
      break;
    case PROP_SUBMODEL_ARRAY:
    {
      NcmObjArray *oa = ncm_obj_array_new ();
      guint i;

      for (i = 0; i < self->submodel_array->len; i++)
      {
        NcmModel *submodel = g_ptr_array_index (self->submodel_array, i);

        ncm_obj_array_add (oa, G_OBJECT (submodel));
      }

      g_value_take_boxed (value, oa);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gboolean
_ncm_model_valid (NcmModel *model)
{
  NCM_UNUSED (model);

  return TRUE;
}

static void _ncm_model_add_submodel (NcmModel *model, NcmModel *submodel);

static void
ncm_model_class_init (NcmModelClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_model_constructed;
  object_class->set_property = &_ncm_model_set_property;
  object_class->get_property = &_ncm_model_get_property;
  object_class->dispose      = &_ncm_model_dispose;
  object_class->finalize     = &_ncm_model_finalize;

  klass->valid        = &_ncm_model_valid;
  klass->set_property = NULL;
  klass->get_property = NULL;

  klass->model_id          = -1;
  klass->can_stack         = FALSE;
  klass->main_model_id     = -1;
  klass->is_submodel       = FALSE;
  klass->name              = NULL;
  klass->nick              = NULL;
  klass->nonparam_prop_len = 0;
  klass->sparam_len        = 0;
  klass->vparam_len        = 0;
  klass->sparam            = NULL;
  klass->vparam            = NULL;

  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Model's name",
                                                        NULL,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NICK,
                                   g_param_spec_string ("nick",
                                                        NULL,
                                                        "Model's nick",
                                                        NULL,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SPARAMS_LEN,
                                   g_param_spec_uint ("scalar-params-len",
                                                      NULL,
                                                      "Number of scalar parameters",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_VPARAMS_LEN,
                                   g_param_spec_uint ("vector-params-len",
                                                      NULL,
                                                      "Number of vector parameters",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_IMPLEMENTATION,
                                   g_param_spec_uint64  ("implementation",
                                                         NULL,
                                                         "Bitwise specification of functions implementation",
                                                         0, G_MAXUINT64, 0,
                                                         G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PTYPES,
                                   g_param_spec_boxed  ("params-types",
                                                        NULL,
                                                        "Parameters' types",
                                                        G_TYPE_ARRAY,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SPARAM_ARRAY,
                                   g_param_spec_boxed ("sparam-array",
                                                       NULL,
                                                       "NcmModel array of NcmSParam",
                                                       NCM_TYPE_OBJ_DICT_INT,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_REPARAM,
                                   g_param_spec_object  ("reparam",
                                                         NULL,
                                                         "Model reparametrization",
                                                         NCM_TYPE_REPARAM,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SUBMODEL_ARRAY,
                                   g_param_spec_boxed ("submodel-array",
                                                       NULL,
                                                       "NcmModel array of submodels",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  klass->add_submodel = &_ncm_model_add_submodel;
}

/*
 * ncm_model_class_get_property:
 * @object: a GObject descending from NcmModel
 * @prop_id: the gobject property id
 * @value: a GValue
 * @pspec: a GParamSpec
 *
 * get_property function, it should be used only when implementing NcmModels
 * in binded languages.
 *
 */
static void
ncm_model_class_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModel *model              = NCM_MODEL (object);
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (g_type_class_peek_static (pspec->owner_type));
  const guint sparam_id        = prop_id       - model_class->nonparam_prop_len + model_class->parent_sparam_len;
  const guint vparam_id        = sparam_id     - model_class->sparam_len        + model_class->parent_vparam_len;
  const guint vparam_len_id    = vparam_id     - model_class->vparam_len        + model_class->parent_vparam_len;
  const guint sparam_fit_id    = vparam_len_id - model_class->vparam_len        + model_class->parent_sparam_len;
  const guint vparam_fit_id    = sparam_fit_id - model_class->sparam_len        + model_class->parent_vparam_len;

  if ((prop_id < model_class->nonparam_prop_len) && model_class->get_property)
  {
    model_class->get_property (object, prop_id, value, pspec);
  }
  else if (sparam_id < model_class->sparam_len)
  {
    g_value_set_double (value, ncm_model_orig_param_get (model, sparam_id));
  }
  else if (vparam_id < model_class->vparam_len)
  {
    NcmVector *vp = ncm_model_orig_vparam_get_vector (model, vparam_id);

    g_value_take_object (value, vp);
  }
  else if (vparam_len_id < model_class->vparam_len)
  {
    g_value_set_uint (value, ncm_model_vparam_len (model, vparam_len_id));
  }
  else if (sparam_fit_id < model_class->sparam_len)
  {
    g_value_set_boolean (value, ncm_model_param_get_ftype (model, sparam_fit_id) == NCM_PARAM_TYPE_FREE ? TRUE : FALSE);
  }
  else if (vparam_fit_id < model_class->vparam_len)
  {
    gsize n = g_array_index (self->vparam_len, guint, vparam_fit_id);
    GVariantBuilder builder;
    GVariant *var;
    guint i;

    g_variant_builder_init (&builder, G_VARIANT_TYPE ("ab"));

    for (i = 0; i < n; i++)
    {
      guint pid      = ncm_model_vparam_index (model, vparam_fit_id, i);
      gboolean tofit = ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE ? TRUE : FALSE;

      g_variant_builder_add (&builder, "b", tofit);
    }

    var = g_variant_builder_end (&builder);
    g_variant_ref_sink (var);

    g_value_take_variant (value, var);
  }
  else
  {
    g_assert_not_reached ();
  }
}

/*
 * ncm_model_class_set_property:
 * @object: a GObject descending from NcmModel
 * @prop_id: the gobject property id
 * @value: a GValue
 * @pspec: a GParamSpec
 *
 * get_property function, it should be used only when implementing NcmModels
 * in binded languages.
 *
 */
static void
ncm_model_class_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModel *model              = NCM_MODEL (object);
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (g_type_class_peek_static (pspec->owner_type));
  const guint sparam_id        = prop_id       - model_class->nonparam_prop_len + model_class->parent_sparam_len;
  const guint vparam_id        = sparam_id     - model_class->sparam_len        + model_class->parent_vparam_len;
  const guint vparam_len_id    = vparam_id     - model_class->vparam_len        + model_class->parent_vparam_len;
  const guint sparam_fit_id    = vparam_len_id - model_class->vparam_len        + model_class->parent_sparam_len;
  const guint vparam_fit_id    = sparam_fit_id - model_class->sparam_len        + model_class->parent_vparam_len;

  /*printf ("[%u %u] [%u %u] [%u %u] [%u %u] [%u %u] [%u %u]\n", prop_id, model_class->nonparam_prop_len, sparam_id, model_class->sparam_len, vparam_id, model_class->vparam_len, vparam_len_id, model_class->vparam_len, sparam_fit_id, model_class->sparam_len, vparam_fit_id, model_class->vparam_len);*/

  if ((prop_id < model_class->nonparam_prop_len) && model_class->set_property)
  {
    model_class->set_property (object, prop_id, value, pspec);
  }
  else if (sparam_id < model_class->sparam_len)
  {
    gdouble val = g_value_get_double (value);

    ncm_model_orig_param_set (model, sparam_id, val);
  }
  else if (vparam_id < model_class->vparam_len)
  {
    NcmVector *vals = g_value_get_object (value);
    guint n         = ncm_vector_len (vals);

    if (n != g_array_index (self->vparam_len, guint, vparam_id))
      g_error ("set_property: cannot set value of vector parameter, vector contains %u elements but vparam dimension is %u",
               n, ncm_model_vparam_len (model, vparam_id));

    ncm_model_orig_vparam_set_vector (model, vparam_id, vals);
  }
  else if (vparam_len_id < model_class->vparam_len)
  {
    NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);
    guint psize                = g_value_get_uint (value);

    if (self->vparam_len->len == 0)
    {
      g_array_set_size (self->vparam_len, model_class->vparam_len);
      g_array_set_size (self->vparam_pos, model_class->vparam_len);
    }
    else
    {
      g_assert_cmpuint (self->vparam_len->len, ==, model_class->vparam_len);
      g_assert_cmpuint (self->vparam_pos->len, ==, model_class->vparam_len);
    }

    g_array_index (self->vparam_len, guint, vparam_len_id) = psize;
  }
  else if (sparam_fit_id < model_class->sparam_len)
  {
    gboolean tofit = g_value_get_boolean (value);

    ncm_model_param_set_ftype (model, sparam_fit_id, tofit ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);
  }
  else if (vparam_fit_id < model_class->vparam_len)
  {
    GVariant *var = g_value_get_variant (value);
    gsize n       = g_variant_n_children (var);
    gsize nv      = g_array_index (self->vparam_len, guint, vparam_fit_id);
    guint i;

    if (n == 1)
    {
      gboolean tofit;
      GVariant *varc = g_variant_get_child_value (var, 0);

      if (g_variant_is_of_type (varc, G_VARIANT_TYPE ("b")))
        tofit = g_variant_get_boolean (varc);
      else if (g_variant_is_of_type (varc, G_VARIANT_TYPE ("i")))
        tofit = g_variant_get_int32 (varc) != 0;
      else
        g_error ("set_property: Cannot convert `%s' variant to an array of booleans", g_variant_get_type_string (varc));

      g_variant_unref (varc);

      for (i = 0; i < nv; i++)
      {
        guint pid = ncm_model_vparam_index (model, vparam_fit_id, i);

        ncm_model_param_set_ftype (model, pid, tofit ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);
      }
    }
    else if (n != nv)
    {
      g_error ("set_property: cannot set fit type of vector parameter, variant contains %zu children but vector dimension is %u", n, g_array_index (self->vparam_len, guint, vparam_fit_id));
    }
    else
    {
      if (g_variant_is_of_type (var, G_VARIANT_TYPE ("ab")))
      {
        for (i = 0; i < n; i++)
        {
          guint pid = ncm_model_vparam_index (model, vparam_fit_id, i);
          gboolean tofit;

          g_variant_get_child (var, i, "b", &tofit);
          ncm_model_param_set_ftype (model, pid, tofit ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);
        }
      }
      else if (g_variant_is_of_type (var, G_VARIANT_TYPE ("ai")))
      {
        for (i = 0; i < n; i++)
        {
          guint pid = ncm_model_vparam_index (model, vparam_fit_id, i);
          gint tofit;

          g_variant_get_child (var, i, "i", &tofit);
          ncm_model_param_set_ftype (model, pid, tofit ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);
        }
      }
      else
      {
        g_error ("set_property: Cannot convert `%s' variant to an array of booleans", g_variant_get_type_string (var));
      }
    }
  }
  else
  {
    g_assert_not_reached ();
  }
}

/**
 * ncm_model_class_add_params:
 * @model_class: a #NcmModelClass
 * @sparam_len: number of scalar paramters
 * @vparam_len: number of vector parameters
 * @nonparam_prop_len: number of properties
 *
 * Class function to be used when implementing NcmModels, it defines the number
 * of scalar and vector parameters and the number of properties of the model.
 *
 */
void
ncm_model_class_add_params (NcmModelClass *model_class, guint sparam_len, guint vparam_len, guint nonparam_prop_len)
{
  GObjectClass *object_class = G_OBJECT_CLASS (model_class);

  object_class->set_property     = &ncm_model_class_set_property;
  object_class->get_property     = &ncm_model_class_get_property;
  model_class->parent_sparam_len = model_class->sparam_len;
  model_class->parent_vparam_len = model_class->vparam_len;
  model_class->sparam_len       += sparam_len;
  model_class->vparam_len       += vparam_len;
  model_class->nonparam_prop_len = nonparam_prop_len;

  if (model_class->sparam_len > 0)
  {
    if (model_class->sparam == NULL)
    {
      model_class->sparam = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_sparam_free);
      g_ptr_array_set_size (model_class->sparam, model_class->sparam_len);
    }
    else
    {
      GPtrArray *sparam = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_sparam_free);
      guint i;

      g_ptr_array_set_size (sparam, model_class->sparam_len);

      /* Copy all parent params info */
      for (i = 0; i < model_class->parent_sparam_len; i++)
        g_ptr_array_index (sparam, i) = ncm_sparam_copy (g_ptr_array_index (model_class->sparam, i));

      model_class->sparam = sparam;
    }
  }

  if (model_class->vparam_len > 0)
  {
    if (model_class->vparam == NULL)
    {
      model_class->vparam = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vparam_free);
      g_ptr_array_set_size (model_class->vparam, model_class->vparam_len);
    }
    else
    {
      GPtrArray *vparam = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_vparam_free);
      guint i;

      g_ptr_array_set_size (vparam, model_class->vparam_len);

      /* Copy all parent params info */
      for (i = 0; i < model_class->parent_vparam_len; i++)
        g_ptr_array_index (vparam, i) = ncm_vparam_copy (g_ptr_array_index (model_class->vparam, i));

      model_class->vparam = vparam;
    }
  }
}

/**
 * ncm_model_class_set_name_nick:
 * @model_class: a #NcmModelClass
 * @name: name and/or very short description of the model
 * @nick: model nickname
 *
 * Attributes @name and @nick, respectively, as the name and nickname of the model.
 *
 */
void
ncm_model_class_set_name_nick (NcmModelClass *model_class, const gchar *name, const gchar *nick)
{
  /*g_clear_pointer (&model_class->name, g_free);*/
  /*g_clear_pointer (&model_class->nick, g_free);*/
  model_class->name = g_strdup (name);
  model_class->nick = g_strdup (nick);
}

/**
 * ncm_model_class_set_sparam_obj:
 * @model_class: a #NcmModelClass
 * @sparam_id: id of the scalar parameter
 * @sparam: a #NcmSParam
 *
 * Sets the @sparam as the @sparam_id-th scalar parameter of the model.
 *
 */
void
ncm_model_class_set_sparam_obj (NcmModelClass *model_class, guint sparam_id, NcmSParam *sparam)
{
  GObjectClass *object_class = G_OBJECT_CLASS (model_class);
  const guint prop_id        = sparam_id - model_class->parent_sparam_len + model_class->nonparam_prop_len;
  const guint prop_fit_id    = prop_id + (model_class->sparam_len - model_class->parent_sparam_len) + 2 * (model_class->vparam_len - model_class->parent_vparam_len);

  if (sparam_id >= model_class->sparam_len)
    g_error ("ncm_model_class_set_sparam: cannot set parameter `%s` in model ``%s''. "
             "Parameter id %u is out of range (0-%u)",
             ncm_sparam_name (sparam), model_class->name, sparam_id + 1, model_class->sparam_len);

  g_assert_cmpint (prop_id, >, 0);

  if (g_ptr_array_index (model_class->sparam, sparam_id) != NULL)
    g_error ("Scalar Parameter: %u is already set.", sparam_id);

  g_ptr_array_index (model_class->sparam, sparam_id) = ncm_sparam_ref (sparam);

  g_object_class_install_property (object_class, prop_id,
                                   g_param_spec_double (ncm_sparam_name (sparam), NULL, ncm_sparam_symbol (sparam),
 /*ncm_sparam_get_lower_bound (sparam), ncm_sparam_get_upper_bound (sparam),*/ /*old behavior*/
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE,
                                                        ncm_sparam_get_default_value (sparam),
                                                        G_PARAM_READWRITE));

  {
    gchar *param_fit_name   = g_strdup_printf ("%s-fit", ncm_sparam_name (sparam));
    gchar *param_fit_symbol = g_strdup_printf ("%s:fit", ncm_sparam_symbol (sparam));

    g_object_class_install_property (object_class, prop_fit_id,
                                     g_param_spec_boolean (param_fit_name, NULL, param_fit_symbol,
                                                           ncm_sparam_get_fit_type (sparam) == NCM_PARAM_TYPE_FREE ? TRUE : FALSE,
                                                           G_PARAM_READWRITE));
    g_free (param_fit_name);
    g_free (param_fit_symbol);
  }
}

/**
 * ncm_model_class_set_vparam_obj:
 * @model_class: a #NcmModelClass
 * @vparam_id: id of the vector parameter
 * @vparam: a #NcmVParam
 *
 * Sets the @vparam as the @vparam_id-th vector parameter of the model.
 *
 */
void
ncm_model_class_set_vparam_obj (NcmModelClass *model_class, guint vparam_id, NcmVParam *vparam)
{
  GObjectClass *object_class = G_OBJECT_CLASS (model_class);
  const guint prop_id        = vparam_id + model_class->nonparam_prop_len - model_class->parent_vparam_len + (model_class->sparam_len - model_class->parent_sparam_len);
  const guint prop_len_id    = prop_id + (model_class->vparam_len - model_class->parent_vparam_len);
  const guint prop_fit_id    = prop_len_id + (model_class->vparam_len - model_class->parent_vparam_len) + (model_class->sparam_len - model_class->parent_sparam_len);

  if (vparam_id >= model_class->vparam_len)
    g_error ("ncm_model_class_set_vparam: setting parameter %u-th of %u (%s) parameters declared for model ``%s''.",
             vparam_id + 1, model_class->vparam_len, ncm_vparam_name (vparam), model_class->name);

  g_assert (prop_id > 0);
  g_assert (prop_len_id > 0);
  /*g_assert_cmpuint (default_length, >, 0);*/

  if (g_ptr_array_index (model_class->vparam, vparam_id) != NULL)
    g_error ("Vector Parameter: %u is already set.", vparam_id);

  g_ptr_array_index (model_class->vparam, vparam_id) = ncm_vparam_ref (vparam);
  g_object_class_install_property (object_class, prop_id,
                                   g_param_spec_object (ncm_vparam_name (vparam), NULL, ncm_vparam_symbol (vparam),
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE));
  {
    gchar *param_length_name   = g_strdup_printf ("%s-length", ncm_vparam_name (vparam));
    gchar *param_length_symbol = g_strdup_printf ("%s:length", ncm_vparam_symbol (vparam));
    gchar *param_fit_name      = g_strdup_printf ("%s-fit", ncm_vparam_name (vparam));
    gchar *param_fit_symbol    = g_strdup_printf ("%s:fit", ncm_vparam_symbol (vparam));

    g_object_class_install_property (object_class, prop_len_id,
                                     g_param_spec_uint (param_length_name,
                                                        NULL,
                                                        param_length_symbol,
                                                        0, G_MAXUINT, ncm_vparam_len (vparam),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));

    g_object_class_install_property (object_class, prop_fit_id,
                                     g_param_spec_variant (param_fit_name, NULL, param_fit_symbol,
                                                           G_VARIANT_TYPE ("ab"), NULL,
                                                           G_PARAM_READWRITE));
    g_free (param_length_name);
    g_free (param_length_symbol);
    g_free (param_fit_name);
    g_free (param_fit_symbol);
  }
}

/**
 * ncm_model_class_set_sparam:
 * @model_class: a #NcmModelClass
 * @sparam_id: id of the scalar parameter
 * @symbol: symbol of the scalar parameter
 * @name: name of the sacalar parameter
 * @lower_bound: lower-bound value
 * @upper_bound: upper-bound value
 * @scale: parameter scale
 * @abstol: parameter absolute tolerance
 * @default_value: default value
 * @ppt: a #NcmParamType
 *
 * Helper function to set a scalar parameter. It creates a #NcmSParam object
 * and calls ncm_model_class_set_sparam_obj().
 *
 */
void
ncm_model_class_set_sparam (NcmModelClass *model_class, guint sparam_id, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt)
{
  NcmSParam *sparam = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_value, ppt);

  ncm_model_class_set_sparam_obj (model_class, sparam_id, sparam);

  ncm_sparam_free (sparam);
}

/**
 * ncm_model_class_set_vparam:
 * @model_class: a #NcmModelClass
 * @vparam_id: id of the vector parameter
 * @default_length: default length of the vector parameter
 * @symbol: symbol of the vector parameter
 * @name: name of the vector parameter
 * @lower_bound: parameter lower bound
 * @upper_bound: parameter upper bound
 * @scale: parameter scale
 * @abstol: parameter absolute tolerance
 * @default_value: default value
 * @ppt: a #NcmParamType
 *
 * Helper function to set a vector parameter. It creates a #NcmVParam object
 * and calls ncm_model_class_set_vparam_obj().
 *
 */
void
ncm_model_class_set_vparam (NcmModelClass *model_class, guint vparam_id, guint default_length, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt)
{
  NcmVParam *vparam = ncm_vparam_full_new (default_length, name, symbol, lower_bound, upper_bound, scale, abstol, default_value, ppt);

  ncm_model_class_set_vparam_obj (model_class, vparam_id, vparam);

  ncm_vparam_free (vparam);
}

/**
 * ncm_model_class_check_params_info:
 * @model_class: a #NcmModelClass
 *
 * Class function to be used when implementing NcmModels, it checks if the
 * parameters information is correctly set. It must be called after all
 * parameters are set during the class initialization.
 *
 */
void
ncm_model_class_check_params_info (NcmModelClass *model_class)
{
  gulong i;
  guint total_params_len = model_class->sparam_len + model_class->vparam_len;

  if ((total_params_len == 0) && (model_class->nonparam_prop_len == 0))
    g_error ("Class size or params not initialized, call ncm_model_class_add_params.");

  for (i = 0; i < model_class->sparam_len; i++)
  {
    if (g_ptr_array_index (model_class->sparam, i) == NULL)
      g_error ("Class (%s) didn't initialized scalar parameter %lu/%u", model_class->name ? model_class->name : "no-name", i + 1, model_class->sparam_len);

    /* g_debug ("Model[%s][%s] id %lu\n", model_class->name, ((NcmSParam *)g_ptr_array_index (model_class->params_info, i))->name, i); */
  }

  for (i = 0; i < model_class->vparam_len; i++)
  {
    if (g_ptr_array_index (model_class->vparam, i) == NULL)
      g_error ("Class (%s) didn't initialized vector parameter %lu/%u", model_class->name ? model_class->name : "no-name", i + 1, model_class->vparam_len);

    /* g_debug ("Model[%s][%s] id %lu\n", model_class->name, ((NcmSParam *)g_ptr_array_index (model_class->params_info, i))->name, i); */
  }

  {
    GObjectClass *object_class = G_OBJECT_CLASS (model_class);

    if (object_class->set_property != &ncm_model_class_set_property)
      g_error ("Class (%s) is using object_class set_property, use model_class set_property instead.", model_class->name ? model_class->name : "no-name");

    if (object_class->get_property != &ncm_model_class_get_property)
      g_error ("Class (%s) is using object_class get_property, use model_class get_property instead.", model_class->name ? model_class->name : "no-name");
  }

  {
    if (model_class->nonparam_prop_len > 1)
    {
      if (model_class->set_property == NULL)
        g_error ("Class (%s) uses non parameter properties but does not set model_class->set_property.", model_class->name ? model_class->name : "no-name");

      if (model_class->get_property == NULL)
        g_error ("Class (%s) uses non parameter properties but does not set model_class->get_property.", model_class->name ? model_class->name : "no-name");
    }
  }
}

/**
 * ncm_model_class_add_impl_opts: (skip)
 * @model_class: a #NcmModelClass
 * @opt1: first option
 * @...: other options, must end with -1
 *
 * Class function to be used when implementing NcmModels, it defines the
 * implementation options of the model.
 *
 */
void
ncm_model_class_add_impl_opts (NcmModelClass *model_class, gint opt1, ...)
{
  gint opt_i;
  va_list ap;

  va_start (ap, opt1);

  model_class->impl_flag = model_class->impl_flag | NCM_MODEL_OPT2IMPL (opt1);

  while ((opt_i = va_arg (ap, gint)) != -1)
  {
    model_class->impl_flag = model_class->impl_flag | NCM_MODEL_OPT2IMPL (opt_i);
  }

  va_end (ap);
}

/**
 * ncm_model_class_add_impl_flag: (skip)
 * @model_class: a #NcmModelClass
 * @flag: implementation flag
 *
 * Class function to be used when implementing NcmModels, it defines the
 * implementation flags of the model.
 *
 */
void
ncm_model_class_add_impl_flag (NcmModelClass *model_class, guint64 flag)
{
  model_class->impl_flag = model_class->impl_flag | flag;
}

/**
 * ncm_model_dup:
 * @model: a #NcmModel
 * @ser: a #NcmSerialize
 *
 * Duplicates @model by serializing and deserializing it.
 *
 * Returns: (transfer full): a duplicate of @model.
 */
NcmModel *
ncm_model_dup (NcmModel *model, NcmSerialize *ser)
{
  return NCM_MODEL (ncm_serialize_dup_obj (ser, G_OBJECT (model)));
}

/**
 * ncm_model_ref:
 * @model: a #NcmModel
 *
 * Increments the reference count of @model by one.
 *
 * Returns: (transfer full): the same @model.
 */
NcmModel *
ncm_model_ref (NcmModel *model)
{
  return g_object_ref (model);
}

/**
 * ncm_model_free:
 * @model: a #NcmModel
 *
 * Atomically decrements the reference count of @model by one. If the reference count drops to 0,
 * all memory allocated by @model is released.
 *
 */
void
ncm_model_free (NcmModel *model)
{
  g_object_unref (model);
}

/**
 * ncm_model_clear:
 * @model: a #NcmModel
 *
 * Atomically decrements the reference count of @model by one. If the reference count drops to 0,
 * all memory allocated by @model is released. Set pointer to NULL.
 *
 */
void
ncm_model_clear (NcmModel **model)
{
  g_clear_object (model);
}

/**
 * ncm_model_set_reparam:
 * @model: a #NcmModel
 * @reparam: a #NcmReparam
 * @error: a #GError
 *
 * Sets the reparametrization of @model to @reparam.
 *
 */
void
ncm_model_set_reparam (NcmModel *model, NcmReparam *reparam, GError **error)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_return_if_fail (error == NULL || *error == NULL);

  if (reparam != NULL)
  {
    GType compat_type = ncm_reparam_get_compat_type (reparam);

    if (!g_type_is_a (G_OBJECT_TYPE (model), compat_type))
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_REPARAM_INCOMPATIBLE,
                                  "ncm_model_set_reparam: model `%s' is not compatible with the reparametrization `%s'",
                                  g_type_name (G_OBJECT_TYPE (model)), g_type_name (compat_type));

      return;
    }

    if (self->reparam != reparam)
    {
      ncm_reparam_clear (&self->reparam);
      self->reparam = ncm_reparam_ref (reparam);
    }

    if (self->p != ncm_reparam_peek_params (self->reparam))
    {
      ncm_vector_clear (&self->p);
      self->p = ncm_vector_ref (ncm_reparam_peek_params (self->reparam));
    }

    ncm_reparam_old2new (self->reparam, model);
  }
  else
  {
    _ncm_model_sparams_remove_reparam (model);
  }
}

/**
 * ncm_model_is_equal:
 * @model1: a #NcmModel
 * @model2: a #NcmModel
 *
 * Compares if model1 and model2 are the same, with same dimension and
 * reparametrization.
 *
 */
gboolean
ncm_model_is_equal (NcmModel *model1, NcmModel *model2)
{
  NcmModelPrivate * const self1 = ncm_model_get_instance_private (model1);
  NcmModelPrivate * const self2 = ncm_model_get_instance_private (model2);

  if (G_OBJECT_TYPE (model1) != G_OBJECT_TYPE (model2))
    return FALSE;

  if (self1->total_len != self2->total_len)
    return FALSE;

  if (self1->reparam)
  {
    if (self2->reparam == NULL)
      return FALSE;

    if (G_OBJECT_TYPE (self1->reparam) != G_OBJECT_TYPE (self2->reparam))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_model_get_reparam:
 * @model: a #NcmModel
 *
 * Gets the reparametrization of @model or NULL if it does not have one.
 *
 * Returns: (transfer full): the reparametrization of @model or NULL if it does not have one.
 */
NcmReparam *
ncm_model_get_reparam (NcmModel *model)
{
  NcmReparam *reparam;

  g_object_get (model, "reparam", &reparam, NULL);
  g_assert (NCM_IS_REPARAM (reparam));

  return reparam;
}

/**
 * ncm_model_params_set_default:
 * @model: a #NcmModel
 *
 * Sets the models parameters to their default values.
 *
 */
void
ncm_model_params_set_default (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
  {
    const NcmSParam *p = _ncm_model_param_peek_desc (model, i, NULL);

    ncm_vector_set (self->p, i, ncm_sparam_get_default_value (p));
  }

  ncm_model_params_update (model);
}

/**
 * ncm_model_params_save_as_default:
 * @model: a #NcmModel
 *
 * Saves the current parameters as the default values.
 *
 */
void
ncm_model_params_save_as_default (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
  {
    gboolean is_original;
    NcmSParam *p = _ncm_model_param_peek_desc (model, i, &is_original);

    ncm_sparam_set_default_value (p, ncm_vector_get (self->p, i));

    if (is_original)
      g_array_index (self->sparam_modified, gboolean, i) = TRUE;
  }
}

/**
 * ncm_model_params_copyto:
 * @model: a #NcmModel
 * @model_dest: a #NcmModel
 *
 * Copies the parameters of @model to @model_dest.
 *
 */
void
ncm_model_params_copyto (NcmModel *model, NcmModel *model_dest)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_assert (ncm_model_is_equal (model, model_dest));
  ncm_model_params_set_vector (model_dest, self->p);
}

/**
 * ncm_model_params_set_all:
 * @model: a #NcmModel
 * @...: a list of doubles
 *
 * Sets all parameters of @model to the values passed as arguments.
 * The number of arguments must be equal to the number of parameters.
 *
 */
void
ncm_model_params_set_all (NcmModel *model, ...)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  guint i;
  va_list ap;

  va_start (ap, model);

  for (i = 0; i < self->total_len; i++)
    ncm_vector_set (self->p, i, va_arg (ap, gdouble));

  va_end (ap);

  ncm_model_params_update (model);

  return;
}

/**
 * ncm_model_params_set_all_data:
 * @model: a #NcmModel
 * @data: an array of doubles
 *
 * Sets all parameters of @model to the values passed as arguments.
 * The size of the array must be equal to the number of parameters.
 *
 */
void
ncm_model_params_set_all_data (NcmModel *model, gdouble *data)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  guint i;

  for (i = 0; i < self->total_len; i++)
    ncm_vector_set (self->p, i, data[i]);

  ncm_model_params_update (model);

  return;
}

/**
 * ncm_model_params_set_vector:
 * @model: a #NcmModel
 * @v: a #NcmVector
 *
 * Sets all parameters of @model to the values of @v.
 * The size of @v must be equal to the number of parameters.
 *
 */
void
ncm_model_params_set_vector (NcmModel *model, NcmVector *v)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_memcpy (self->p, v);
  ncm_model_params_update (model);
}

/**
 * ncm_model_params_set_model:
 * @model: a #NcmModel
 * @model_src: a #NcmModel
 *
 * Sets all parameters of @model to the values of @model_src.
 *
 */
void
ncm_model_params_set_model (NcmModel *model, NcmModel *model_src)
{
  NcmModelPrivate * const self_src = ncm_model_get_instance_private (model_src);

  g_assert (ncm_model_is_equal (model, model_src));
  ncm_model_params_set_vector (model, self_src->p);
}

/**
 * ncm_model_params_print_all: (skip)
 * @model: a #NcmModel
 * @out: a file handle
 *
 * Prints all parameters of @model to @out.
 *
 */
void
ncm_model_params_print_all (NcmModel *model, FILE *out)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
    fprintf (out, "  % 22.16g", ncm_vector_get (self->p, i));

  fprintf (out, "\n");
  fflush (out);

  return;
}

/**
 * ncm_model_orig_params_log_all:
 * @model: a #NcmModel
 *
 * Logs all original parameters of @model. That is if there is a reparametrization
 * set, it return the values of the original parameters.
 *
 */
void
ncm_model_orig_params_log_all (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
    g_message ("  % 20.16g", ncm_vector_fast_get (self->params, i));

  g_message ("\n");

  return;
}

/**
 * ncm_model_params_log_all:
 * @model: a #NcmModel
 *
 * Logs all parameters of @model. It prints the values of the parameters
 * in the current reparametrization.
 *
 */
void
ncm_model_params_log_all (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
    g_message ("  % 22.16g", ncm_vector_get (self->p, i));

  g_message ("\n");

  return;
}

/**
 * ncm_model_params_get_all:
 * @model: a #NcmModel
 *
 * Creates a #NcmVector with all parameters of @model.
 *
 * Returns: (transfer full): a #NcmVector with all parameters of @model.
 */
NcmVector *
ncm_model_params_get_all (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return ncm_vector_dup (self->p);
}

/**
 * ncm_model_params_valid:
 * @model: a #NcmModel
 *
 * Check whenever the parameters are valid.
 *
 * Returns: TRUE if the parameter are valid.
 */
gboolean
ncm_model_params_valid (NcmModel *model)
{
  NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);

  return model_class->valid (model);
}

/**
 * ncm_model_params_valid_bounds:
 * @model: a #NcmModel
 *
 * Check whenever the parameters respect the bounds.
 *
 * Returns: if the parameter respect the bounds.
 */
gboolean
ncm_model_params_valid_bounds (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
  {
    const gdouble lb  = ncm_model_param_get_lower_bound (model, i);
    const gdouble ub  = ncm_model_param_get_upper_bound (model, i);
    const gdouble val = ncm_model_param_get (model, i);

    if ((val < lb) || (val > ub))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_model_id:
 * @model: a #NcmModel
 *
 * Gets the model id of @model.
 *
 * Returns: The model id of @model.
 */
NcmModelID
ncm_model_id (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->model_id;
}

/**
 * ncm_model_id_by_type:
 * @model_type: a GType
 * @error: a #GError
 *
 * Gets the model id of a model type. It is an error to call this function
 * with a type that is not a subclass of #NcmModel.
 *
 * Returns: The model id of @model_type.
 */
NcmModelID
ncm_model_id_by_type (GType model_type, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, 0);

  if (!g_type_is_a (model_type, NCM_TYPE_MODEL))
  {
    ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_INVALID_TYPE,
                                "ncm_model_id_by_type: type (%s) is not a %s",
                                g_type_name (model_type),
                                g_type_name (NCM_TYPE_MODEL));

    return 0;
  }
  else
  {
    NcmModelClass *model_class = NCM_MODEL_CLASS (g_type_class_ref (model_type));
    NcmModelID id              = model_class->model_id;

    g_type_class_unref (model_class);

    return id;
  }
}

/**
 * ncm_model_check_impl_flag:
 * @model: a #NcmModel
 * @impl: implementation flag
 *
 * Checks if the model implements the @impl flag.
 *
 * Returns: TRUE if the model implements the @impl flag.
 */
gboolean
ncm_model_check_impl_flag (NcmModel *model, guint64 impl)
{
  if (impl == 0)
    return TRUE;
  else
    return ((NCM_MODEL_GET_CLASS (model)->impl_flag & impl) == impl);
}

/**
 * ncm_model_check_impl_opt:
 * @model: a #NcmModel
 * @opt: implementation option
 *
 * Checks if the model implements the @opt option.
 *
 * Returns: TRUE if the model implements the @opt option.
 */
gboolean
ncm_model_check_impl_opt (NcmModel *model, gint opt)
{
  guint64 flag = NCM_MODEL_OPT2IMPL (opt);

  return ncm_model_check_impl_flag (model, flag);
}

/**
 * ncm_model_check_impl_opts:
 * @model: a #NcmModel
 * @opt1: first implementation option
 * @...: implementation options, must end with -1
 *
 * Checks if the model implements all the @opt1, @opt2, ... options.
 * The last argument must be -1.
 *
 * Returns: TRUE if the model implements all the @opt1, @opt2, ... options.
 */
gboolean
ncm_model_check_impl_opts (NcmModel *model, gint opt1, ...)
{
  guint64 flag = 0;
  gint opt_i;

  va_list ap;

  va_start (ap, opt1);

  flag = flag | NCM_MODEL_OPT2IMPL (opt1);

  while ((opt_i = va_arg (ap, gint)) != -1)
  {
    flag = flag | NCM_MODEL_OPT2IMPL (opt_i);
  }

  va_end (ap);

  return ncm_model_check_impl_flag (model, flag);
}

/**
 * ncm_model_len:
 * @model: a #NcmModel
 *
 * Count the total number of parameters of the model.
 *
 * Returns: the total number of parameters of the model.
 */
guint
ncm_model_len (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->total_len;
}

/**
 * ncm_model_state_is_update:
 * @model: a #NcmModel
 *
 * Check if the model is updated.
 *
 * Returns: TRUE if the model is updated.
 */
gboolean
ncm_model_state_is_update (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->pkey == self->skey;
}

/**
 * ncm_model_state_set_update:
 * @model: a #NcmModel
 *
 * Set the model as updated.
 *
 */
void
ncm_model_state_set_update (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  self->skey = self->pkey;
}

/**
 * ncm_model_lstate_is_update:
 * @model: a #NcmModel
 * @i: lstate index
 *
 * Check if the @i-th lstate is updated.
 * The parameter @i must be smaller than #NCM_MODEL_MAX_STATES.
 *
 * Returns: whether the @i-th lstate is updated.
 */
gboolean
ncm_model_lstate_is_update (NcmModel *model, guint i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->pkey == self->slkey[i];
}

/**
 * ncm_model_lstate_set_update:
 * @model: a #NcmModel
 * @i: lstate index
 *
 * Sets the @i-th lstate as updated. The parameter @i must be smaller than
 * #NCM_MODEL_MAX_STATES.
 *
 */
void
ncm_model_lstate_set_update (NcmModel *model, guint i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  self->slkey[i] = self->pkey;
}

/**
 * ncm_model_state_mark_outdated:
 * @model: a #NcmModel
 *
 * Set the model as outdated.
 *
 */
void
ncm_model_state_mark_outdated (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  self->pkey++;
}

/**
 * ncm_model_state_get_pkey:
 * @model: a #NcmModel
 *
 * Get the current pkey of the model.
 *
 * Returns: the current pkey of the model.
 */
guint64
ncm_model_state_get_pkey (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->pkey;
}

/**
 * ncm_model_sparam_len:
 * @model: a #NcmModel
 *
 * Count the number of scalar parameters of the model.
 *
 * Returns: the number of scalar parameters of the model.
 */
guint
ncm_model_sparam_len (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->sparam_len;
}

/**
 * ncm_model_vparam_array_len:
 * @model: a #NcmModel
 *
 * Count the number of vector parameters of the model.
 * Note that this function returns the number of vector parameters
 * of the model, not the length of the vector parameters.
 *
 * Returns: the number of vector parameters of the model.
 */
guint
ncm_model_vparam_array_len (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->vparam_len;
}

/**
 * ncm_model_name:
 * @model: a #NcmModel
 *
 * Get the name of the model.
 *
 * Returns: (transfer none): the name of the model.
 */
const gchar *
ncm_model_name (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->name;
}

/**
 * ncm_model_nick:
 * @model: a #NcmModel
 *
 * Get the nick of the model.
 *
 * Returns: (transfer none): the nick of the model.
 */
const gchar *
ncm_model_nick (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->nick;
}

/**
 * ncm_model_peek_reparam:
 * @model: a #NcmModel
 *
 * Peeks the current reparametrization of @model.
 *
 * Returns: (transfer none): the current reparametrization of @model or NULL if it does not have one.
 */
NcmReparam *
ncm_model_peek_reparam (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->reparam;
}

/**
 * ncm_model_param_finite:
 * @model: a #NcmModel
 * @i: parameter index
 *
 * Check if the @i-th parameter is finite.
 *
 * Returns: whether the @i-th parameter is finite.
 */
gboolean
ncm_model_param_finite (NcmModel *model, guint i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  NcmVector *params            = self->reparam ? ncm_reparam_peek_params (self->reparam) : self->params;

  return gsl_finite (ncm_vector_get (params, i));
}

/**
 * ncm_model_params_finite:
 * @model: a #NcmModel
 *
 * Check if all parameters are finite.
 *
 * Returns: whether all parameters are finite.
 */

gboolean
ncm_model_params_finite (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < ncm_model_len (model); i++)
  {
    if (!gsl_finite (ncm_vector_fast_get (self->params, i)))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_model_params_update:
 * @model: a #NcmModel
 *
 * Force the parameters to the update its internal flags and
 * update the original parameters if necessary.
 *
 */

void
ncm_model_params_update (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  self->pkey++;

  if (self->reparam)
    ncm_reparam_new2old (self->reparam, model);
}

/**
 * ncm_model_orig_params_update:
 * @model: a #NcmModel
 *
 * Update the new parameters. It causes an error to call this
 * function with a model without reparametrization.
 *
 */
void
ncm_model_orig_params_update (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  self->pkey++;

  if (self->reparam)
    ncm_reparam_old2new (self->reparam, model);
}

/**
 * ncm_model_orig_params_peek_vector:
 * @model: a #NcmModel
 *
 * Peeks the original parameters vector. This functions is provided for
 * reparametrization implementations and subclassing, do not use it in other contexts.
 *
 * The returned vector is the original parameters vector, that is, the parameters
 * before the reparametrization. It is guaranteed that the returned vector is
 * always the same and will stay valid until the model is destroyed.
 * The vector also always have stride 1, so it is safe to call ncm_vector_fast_get().
 *
 * Returns: (transfer none): the original parameters #NcmVector
 */
NcmVector *
ncm_model_orig_params_peek_vector (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->params;
}

/**
 * ncm_model_vparam_index:
 * @model: a #NcmModel
 * @n: vector parameter index
 * @i: vector component index
 *
 * Given a vector parameter index and a component index, returns the index of the
 * @i-th component of the @n-th vector in the full parameter vector.
 *
 * Returns: index of the i-th component of the n-th vector
 */
guint
ncm_model_vparam_index (NcmModel *model, guint n, guint i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return g_array_index (self->vparam_pos, guint, n) + i;
}

/**
 * ncm_model_vparam_len:
 * @model: a #NcmModel
 * @n: vector parameter index
 *
 * Given a vector parameter index, returns the length of the @n-th vector.
 *
 * Returns: length of the n-th vector
 */
guint
ncm_model_vparam_len (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return g_array_index (self->vparam_len, guint, n);
}

/**
 * ncm_model_set_vparam_len:
 * @model: a #NcmModel
 * @n: vector parameter index
 * @len: vector length
 *
 * Given a vector parameter index, sets the length of the @n-th vector to @len.
 * This function is provided for model implementations, do not use it in
 * other contexts. It will be removed in future versions.
 */
void
ncm_model_set_vparam_len (NcmModel *model, guint n, guint len)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_array_index (self->vparam_len, guint, n) = len;
}

/**
 * ncm_model_param_set0:
 * @model: a #NcmModel
 * @n: parameter index
 * @val: a double
 *
 * Sets the @n-th parameter of @model to @val. This function does not
 * update the model after setting the parameter. It is provided when
 * multiple parameters are set at once the model is updated only once.
 * ncm_model_params_update() must be called after setting all parameters.
 *
 */
void
ncm_model_param_set0 (NcmModel *model, guint n, gdouble val)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_set (self->p, n, val);
}

/**
 * ncm_model_param_set:
 * @model: a #NcmModel
 * @n: parameter index
 * @val: a double
 *
 * Sets the @n-th parameter of @model to @val.
 *
 */
void
ncm_model_param_set (NcmModel *model, guint n, gdouble val)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_set (self->p, n, val);
  ncm_model_params_update (model);

  return;
}

/**
 * ncm_model_param_set_default:
 * @model: a #NcmModel
 * @n: a parameter index
 *
 * Sets the @n-th parameter of @model to its default value.
 *
 */
void
ncm_model_param_set_default (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gboolean is_original;

  ncm_model_param_set (model, n, ncm_sparam_get_default_value (_ncm_model_param_peek_desc (model, n, &is_original)));

  if (is_original)
    g_array_index (self->sparam_modified, gboolean, n) = TRUE;
}

static NcmSParam *
_ncm_model_orig_param_peek_desc (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_assert_cmpuint (n, <, self->total_len);

  return g_ptr_array_index (self->sparams, n);
}

static NcmSParam *
_ncm_model_param_peek_desc (NcmModel *model, guint n, gboolean *is_original)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  NcmReparam *reparam          = ncm_model_peek_reparam (model);

  g_assert_cmpuint (n, <, self->total_len);

  if (is_original != NULL)
    *is_original = TRUE;

  if (reparam != NULL)
  {
    NcmSParam *sp = ncm_reparam_peek_param_desc (reparam, n);

    if (sp != NULL)
    {
      if (is_original != NULL)
        *is_original = FALSE;

      return sp;
    }
  }

  return _ncm_model_orig_param_peek_desc (model, n);
}

/**
 * ncm_model_param_get:
 * @model: a #NcmModel
 * @n: a parameter index
 *
 * Gets the @n-th parameter of @model.
 *
 * Returns: the @n-th parameter of @model.
 */
gdouble
ncm_model_param_get (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return ncm_vector_get (self->p, n);
}

/**
 * ncm_model_orig_param_set:
 * @model: a #NcmModel
 * @n: a parameter index
 * @val: a double
 *
 * Sets the @n-th original parameter of @model using the original
 * parametrization to @val.
 *
 */
void
ncm_model_orig_param_set (NcmModel *model, guint n, gdouble val)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_set (self->params, n, val);
  ncm_model_orig_params_update (model);

  return;
}

/**
 * ncm_model_orig_param_get:
 * @model: a #NcmModel
 * @n: a parameter index
 *
 * Gets the @n-th original parameter of @model using the original
 * parametrization.
 *
 * Returns: the @n-th original parameter of @model.
 */
gdouble
ncm_model_orig_param_get (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->params_ptr[n];
}

/**
 * ncm_model_vparam_set:
 * @model: a #NcmModel
 * @n: a vector parameter index
 * @i: a vector component index
 * @val: a double
 *
 * Sets the @i-th component of the @n-th vector parameter of @model to @val.
 *
 */
void
ncm_model_orig_vparam_set (NcmModel *model, guint n, guint i, gdouble val)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_set (self->params, ncm_model_vparam_index (model, n, i), val);
  ncm_model_orig_params_update (model);

  return;
}

/**
 * ncm_model_vparam_get:
 * @model: a #NcmModel
 * @n: a vector parameter index
 * @i: a vector component index
 *
 * Gets the @i-th component of the @n-th vector parameter of @model.
 *
 * Returns: the @i-th component of the @n-th vector parameter of @model.
 */
gdouble
ncm_model_orig_vparam_get (NcmModel *model, guint n, guint i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return ncm_vector_fast_get (self->params, ncm_model_vparam_index (model, n, i));
}

/**
 * ncm_model_orig_vparam_set_vector:
 * @model: a #NcmModel
 * @n: a vector parameter index
 * @val: a #NcmVector
 *
 * Sets the @n-th vector parameter of @model to @val.
 * The size of @val must be equal to the length of the @n-th vector parameter.
 *
 */
void
ncm_model_orig_vparam_set_vector (NcmModel *model, guint n, NcmVector *val)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  ncm_vector_memcpy2 (self->params, val,
                      ncm_model_vparam_index (model, n, 0), 0,
                      ncm_model_vparam_len (model, n));
  ncm_model_orig_params_update (model);
}

/**
 * ncm_model_orig_vparam_get_vector:
 * @model: a #NcmModel
 * @n: a vector parameter index
 *
 * Gets the @n-th vector parameter of @model.
 *
 * Returns: (transfer full): the @n-th vector parameter of @model.
 */
NcmVector *
ncm_model_orig_vparam_get_vector (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  const guint vparam_len       = ncm_model_vparam_len (model, n);

  if (vparam_len > 0)
  {
    NcmVector *val = ncm_vector_new (vparam_len);

    ncm_vector_memcpy2 (val, self->params,
                        0, ncm_model_vparam_index (model, n, 0),
                        ncm_model_vparam_len (model, n));

    return val;
  }
  else
  {
    return NULL;
  }
}

/**
 * ncm_model_orig_param_get_scale:
 * @model: a #NcmModel
 * @n: parameter index.
 *
 * Gets the scale of the original @n-th parameter.
 *
 * Returns: the scale of the original @n-th parameter.
 */
gdouble
ncm_model_orig_param_get_scale (NcmModel *model, guint n)
{
  const gdouble scale = ncm_sparam_get_scale (_ncm_model_orig_param_peek_desc (model, n));

  return scale;
}

/**
 * ncm_model_orig_param_get_lower_bound:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the lower bound of the original @n-th parameter.
 *
 * Returns: the lower bound of the original @n-th parameter.
 */
gdouble
ncm_model_orig_param_get_lower_bound (NcmModel *model, guint n)
{
  const gdouble lb = ncm_sparam_get_lower_bound (_ncm_model_orig_param_peek_desc (model, n));

  return lb;
}

/**
 * ncm_model_orig_param_get_upper_bound:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the upper bound of the original @n-th parameter.
 *
 * Returns: the upper bound of the original @n-th parameter.
 */
gdouble
ncm_model_orig_param_get_upper_bound (NcmModel *model, guint n)
{
  const gdouble ub = ncm_sparam_get_upper_bound (_ncm_model_orig_param_peek_desc (model, n));

  return ub;
}

/**
 * ncm_model_orig_param_get_abstol:
 * @model: a #NcmModel
 * @n: parameter index.
 *
 * Gets the absolute tolerance of the original @n-th parameter.
 *
 * Returns: the absolute tolerance of the original @n-th parameter.
 */
gdouble
ncm_model_orig_param_get_abstol (NcmModel *model, guint n)
{
  const gdouble abstol = ncm_sparam_get_absolute_tolerance (_ncm_model_orig_param_peek_desc (model, n));

  return abstol;
}

/**
 * ncm_model_param_get_scale:
 * @model: a #NcmModel
 * @n: parameter index.
 *
 * Gets the scale of the @n-th parameter.
 *
 * Returns: the scale of the @n-th parameter.
 */
gdouble
ncm_model_param_get_scale (NcmModel *model, guint n)
{
  const gdouble scale = ncm_sparam_get_scale (_ncm_model_param_peek_desc (model, n, NULL));

  return scale;
}

/**
 * ncm_model_param_get_lower_bound:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the lower bound of the @n-th parameter.
 *
 * Returns: the lower bound of the @n-th parameter.
 */
gdouble
ncm_model_param_get_lower_bound (NcmModel *model, guint n)
{
  const gdouble lb = ncm_sparam_get_lower_bound (_ncm_model_param_peek_desc (model, n, NULL));

  return lb;
}

/**
 * ncm_model_param_get_upper_bound:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the upper bound of the @n-th parameter.
 *
 * Returns: the upper bound of the @n-th parameter.
 */
gdouble
ncm_model_param_get_upper_bound (NcmModel *model, guint n)
{
  const gdouble ub = ncm_sparam_get_upper_bound (_ncm_model_param_peek_desc (model, n, NULL));

  return ub;
}

/**
 * ncm_model_param_get_abstol:
 * @model: a #NcmModel
 * @n: parameter index.
 *
 * Gets the absolute tolerance of the @n-th parameter.
 *
 * Returns: the absolute tolerance of the @n-th parameter.
 */
gdouble
ncm_model_param_get_abstol (NcmModel *model, guint n)
{
  const gdouble abstol = ncm_sparam_get_absolute_tolerance (_ncm_model_param_peek_desc (model, n, NULL));

  return abstol;
}

/**
 * ncm_model_param_get_ftype:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the fitting type of the @n-th parameter.
 *
 * Returns: the fitting type of the @n-th parameter.
 */
NcmParamType
ncm_model_param_get_ftype (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return g_array_index (self->ptypes, NcmParamType, n);
}

/**
 * ncm_model_param_set_scale:
 * @model: a #NcmModel
 * @n: parameter index
 * @scale: a double
 *
 * Sets @scale as the scale of the @n-th parameter.
 *
 */
void
ncm_model_param_set_scale (NcmModel *model, guint n, const gdouble scale)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gboolean is_original;

  ncm_sparam_set_scale (_ncm_model_param_peek_desc (model, n, &is_original), scale);

  if (is_original)
    g_array_index (self->sparam_modified, gboolean, n) = TRUE;

  ncm_model_state_mark_outdated (model);
}

/**
 * ncm_model_param_set_lower_bound:
 * @model: a #NcmModel
 * @n: parameter index
 * @lb: lower-bound value
 *
 * Sets @lb as the lower-bound value of the @n-th parameter.
 *
 */
void
ncm_model_param_set_lower_bound (NcmModel *model, guint n, const gdouble lb)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gboolean is_original;

  ncm_sparam_set_lower_bound (_ncm_model_param_peek_desc (model, n, &is_original), lb);

  if (is_original)
    g_array_index (self->sparam_modified, gboolean, n) = TRUE;

  ncm_model_state_mark_outdated (model);
}

/**
 * ncm_model_param_set_upper_bound:
 * @model: a #NcmModel
 * @n: parameter index
 * @ub: upper-bound value
 *
 * Sets @ub as the lower-bound value of the @n-th parameter.
 *
 */
void
ncm_model_param_set_upper_bound (NcmModel *model, guint n, const gdouble ub)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gboolean is_original;

  ncm_sparam_set_upper_bound (_ncm_model_param_peek_desc (model, n, &is_original), ub);

  if (is_original)
    g_array_index (self->sparam_modified, gboolean, n) = TRUE;

  ncm_model_state_mark_outdated (model);
}

/**
 * ncm_model_param_set_abstol:
 * @model: a #NcmModel
 * @n: parameter index
 * @abstol: the absolute tolerance
 *
 * Sets @abstol as the absolute tolerance of the @n-th parameter.
 *
 */
void
ncm_model_param_set_abstol (NcmModel *model, guint n, const gdouble abstol)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gboolean is_original;

  ncm_sparam_set_absolute_tolerance (_ncm_model_param_peek_desc (model, n, &is_original), abstol);

  if (is_original)
    g_array_index (self->sparam_modified, gboolean, n) = TRUE;

  ncm_model_state_mark_outdated (model);
}

/**
 * ncm_model_param_set_ftype:
 * @model: a #NcmModel
 * @n: parameter index
 * @ptype: a #NcmParamType
 *
 * Sets @ptype as #NcmParamType of the @n-th parameter.
 *
 */
void
ncm_model_param_set_ftype (NcmModel *model, guint n, const NcmParamType ptype)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_array_index (self->ptypes, NcmParamType, n) = ptype;
}

/**
 * ncm_model_params_set_default_ftype:
 * @model: a #NcmModel
 *
 * Sets all parameters #NcmParamType to their default values.
 *
 */
void
ncm_model_params_set_default_ftype (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  guint i;

  for (i = 0; i < self->total_len; i++)
  {
    const NcmSParam *p = _ncm_model_param_peek_desc (model, i, NULL);

    g_array_index (self->ptypes, NcmParamType, i) = ncm_sparam_get_fit_type (p);
  }
}

/**
 * ncm_model_param_get_desc:
 * @model: a #NcmModel
 * @param: parameter name
 * @error: a #GError
 *
 * Gets the description of the parameter @param. The output is a GHashTable which
 * contains the following keys:
 * - "name": the name of the parameter.
 * - "symbol": the symbol of the parameter.
 * - "scale": the scale of the parameter.
 * - "lower-bound": the lower bound of the parameter.
 * - "upper-bound": the upper bound of the parameter.
 * - "abstol": the absolute tolerance of the parameter.
 * - "fit": whether the parameter is a fitting parameter.
 * - "value": the current value of the parameter.
 *
 * Returns: (transfer full) (element-type utf8 GValue): the description of the parameter @param.
 */
GHashTable *
ncm_model_param_get_desc (NcmModel *model, gchar *param, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  {
    guint i;
    const gboolean has_param = ncm_model_param_index_from_name (model, param, &i, error);

    if (!has_param)
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
                                  "ncm_model_param_get_desc: model `%s' does not have a parameter called `%s'.",
                                  G_OBJECT_TYPE_NAME (model), param);

      return NULL;
    }
    else
    {
      GHashTable *desc = g_hash_table_new (g_str_hash, g_str_equal);

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_STRING);
        g_value_set_static_string (value, ncm_model_param_name (model, i));
        g_hash_table_insert (desc, g_strdup ("name"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_STRING);
        g_value_set_static_string (value, ncm_model_param_symbol (model, i));
        g_hash_table_insert (desc, g_strdup ("symbol"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_DOUBLE);
        g_value_set_double (value, ncm_model_param_get_scale (model, i));
        g_hash_table_insert (desc, g_strdup ("scale"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_DOUBLE);
        g_value_set_double (value, ncm_model_param_get_lower_bound (model, i));
        g_hash_table_insert (desc, g_strdup ("lower-bound"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_DOUBLE);
        g_value_set_double (value, ncm_model_param_get_upper_bound (model, i));
        g_hash_table_insert (desc, g_strdup ("upper-bound"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_DOUBLE);
        g_value_set_double (value, ncm_model_param_get_abstol (model, i));
        g_hash_table_insert (desc, g_strdup ("abstol"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_BOOLEAN);
        g_value_set_boolean (value, (ncm_model_param_get_ftype (model, i) == NCM_PARAM_TYPE_FREE) ? TRUE : FALSE);
        g_hash_table_insert (desc, g_strdup ("fit"), value);
      }

      {
        GValue *value = g_new0 (GValue, 1);

        g_value_init (value, G_TYPE_DOUBLE);
        g_value_set_double (value, ncm_model_param_get (model, i));
        g_hash_table_insert (desc, g_strdup ("value"), value);
      }

      return desc;
    }
  }
}

/**
 * ncm_model_param_set_desc:
 * @model: a #NcmModel
 * @param: parameter name
 * @desc: (element-type utf8 GValue): a GHashTable
 * @error: a #GError pointer
 *
 * Sets the description of the parameter @param. The input is a GHashTable which
 * may contain the following keys:
 *
 * - "name": the name of the parameter.
 * - "symbol": the symbol of the parameter.
 * - "scale": the scale of the parameter.
 * - "lower-bound": the lower bound of the parameter.
 * - "upper-bound": the upper bound of the parameter.
 * - "abstol": the absolute tolerance of the parameter.
 * - "fit": whether the parameter is a fitting parameter.
 * - "value": the current value of the parameter.
 *
 * Other keys are ignored.
 *
 */
void
ncm_model_param_set_desc (NcmModel *model, gchar *param, GHashTable *desc, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);

  {
    guint i;
    const gboolean has_param = ncm_model_param_index_from_name (model, param, &i, error);

    if (!has_param)
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
                                  "ncm_model_param_set_desc: model `%s' does not have a parameter called `%s'.",
                                  G_OBJECT_TYPE_NAME (model), param);

      return;
    }
    else
    {
      if (g_hash_table_lookup (desc, "scale"))
      {
        GValue *value = g_hash_table_lookup (desc, "scale");

        if (!G_VALUE_HOLDS_DOUBLE (value))
        {
          ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
                                      "ncm_model_param_set_desc: scale must be a double.");

          return;
        }

        ncm_model_param_set_scale (model, i, g_value_get_double (value));
      }

      if (g_hash_table_lookup (desc, "lower-bound"))
      {
        GValue *value = g_hash_table_lookup (desc, "lower-bound");

        if (!G_VALUE_HOLDS_DOUBLE (value))
        {
          ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
                                      "ncm_model_param_set_desc: lower-bound must be a double.");

          return;
        }

        ncm_model_param_set_lower_bound (model, i, g_value_get_double (value));
      }

      if (g_hash_table_lookup (desc, "upper-bound"))
      {
        GValue *value = g_hash_table_lookup (desc, "upper-bound");

        if (!G_VALUE_HOLDS_DOUBLE (value))
        {
          ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
                                      "ncm_model_param_set_desc: upper-bound must be a double.");

          return;
        }

        ncm_model_param_set_upper_bound (model, i, g_value_get_double (value));
      }

      if (g_hash_table_lookup (desc, "abstol"))
      {
        GValue *value = g_hash_table_lookup (desc, "abstol");

        if (!G_VALUE_HOLDS_DOUBLE (value))
        {
          ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
                                      "ncm_model_param_set_desc: abstol must be a double.");

          return;
        }

        ncm_model_param_set_abstol (model, i, g_value_get_double (value));
      }

      if (g_hash_table_lookup (desc, "fit"))
      {
        GValue *value = g_hash_table_lookup (desc, "fit");

        if (!G_VALUE_HOLDS_BOOLEAN (value))
        {
          ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
                                      "ncm_model_param_set_desc: fit must be a boolean.");

          return;
        }

        ncm_model_param_set_ftype (model, i, g_value_get_boolean (value) ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);
      }

      if (g_hash_table_lookup (desc, "value"))
      {
        GValue *value = g_hash_table_lookup (desc, "value");

        if (!G_VALUE_HOLDS_DOUBLE (value))
        {
          ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_TYPE,
                                      "ncm_model_param_set_desc: value must be a double.");

          return;
        }

        ncm_model_param_set (model, i, g_value_get_double (value));
      }

      return;
    }
  }
}

/**
 * ncm_model_orig_param_get_name:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the name of the @n-th original parameter of @model.
 *
 * Returns: (transfer none): the parameter name.
 */
const gchar *
ncm_model_orig_param_name (NcmModel *model, guint n)
{
  return ncm_sparam_name (_ncm_model_orig_param_peek_desc (model, n));
}

/**
 * ncm_model_param_get_name:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the name of the @n-th parameter of @model.
 *
 * Returns: (transfer none): the parameter name.
 */
const gchar *
ncm_model_param_name (NcmModel *model, guint n)
{
  return ncm_sparam_name (_ncm_model_param_peek_desc (model, n, NULL));
}

/**
 * ncm_model_orig_param_get_symbol:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the symbol of the @n-th original parameter of @model.
 *
 * Returns: (transfer none): the parameter symbol.
 */
const gchar *
ncm_model_orig_param_symbol (NcmModel *model, guint n)
{
  return ncm_sparam_symbol (_ncm_model_orig_param_peek_desc (model, n));
}

/**
 * ncm_model_param_names:
 * @model: a #NcmModel
 *
 * Gets an array containing the parameters names.
 *
 * Returns: (transfer container) (element-type utf8): an array containing the parameters names.
 */
GPtrArray *
ncm_model_param_names (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  GPtrArray *names             = g_ptr_array_new ();
  guint i;

  g_ptr_array_set_free_func (names, g_free);

  for (i = 0; i < self->sparams->len; i++)
  {
    gchar *name = g_strdup (ncm_sparam_name (_ncm_model_param_peek_desc (model, i, NULL)));

    g_ptr_array_add (names, name);
  }

  return names;
}

/**
 * ncm_model_param_get_symbol:
 * @model: a #NcmModel
 * @n: parameter index
 *
 * Gets the symbol of the @n-th parameter of @model.
 *
 * Returns: (transfer none): the parameter symbol.
 */
const gchar *
ncm_model_param_symbol (NcmModel *model, guint n)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_assert (n < self->total_len);

  return ncm_sparam_symbol (_ncm_model_param_peek_desc (model, n, NULL));
}

/**
 * ncm_model_orig_param_index_from_name:
 * @model: a #NcmModel
 * @param_name: parameter name
 * @i: (out): parameter index
 *
 * Looks for parameter named @param_name in the original parameters of @model
 * and puts its index in @i and returns TRUE if found.
 *
 * Returns: whether the parameter @param_name is found in the @model.
 */
gboolean
ncm_model_orig_param_index_from_name (NcmModel *model, const gchar *param_name, guint *i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gpointer param_id;
  gboolean found = g_hash_table_lookup_extended (self->sparams_name_id, param_name, NULL, &param_id);

  if (found)
    *i = GPOINTER_TO_UINT (param_id);
  else
    *i = -1;  /* Yup, I know. */

  return found;
}

/**
 * ncm_model_param_index_from_name:
 * @model: a #NcmModel
 * @param_name: parameter name
 * @i: (out): parameter index
 * @error: a #GError
 *
 * Looks for parameter named @param_name in @model and puts its index in @i
 * and returns TRUE if found.
 *
 * Returns: whether the parameter @param_name is found in the @model.
 */
gboolean
ncm_model_param_index_from_name (NcmModel *model, const gchar *param_name, guint *i, GError **error)
{
  NcmReparam *reparam = ncm_model_peek_reparam (model);

  g_return_val_if_fail (error == NULL || *error == NULL, FALSE);

  if (reparam != NULL)
  {
    if (ncm_reparam_index_from_name (reparam, param_name, i))
    {
      return TRUE;
    }
    else if (ncm_model_orig_param_index_from_name (model, param_name, i))
    {
      if (ncm_reparam_peek_param_desc (reparam, *i) != NULL)
      {
        ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_CHANGED,
                                    "ncm_model_param_index_from_name: parameter (%s) was changed by a NcmReparam, it is now named (%s).",
                                    param_name, ncm_sparam_name (ncm_reparam_peek_param_desc (reparam, *i)));

        return FALSE;
      }
      else
      {
        return TRUE;
      }
    }
    else
    {
      return FALSE;
    }
  }
  else
  {
    return ncm_model_orig_param_index_from_name (model, param_name, i);
  }
}

/**
 * ncm_model_param_set_by_name:
 * @model: a #NcmModel
 * @param_name: parameter name
 * @val: parameter value
 * @error: a #GError
 *
 * Sets the parameter value @val by @param_name.
 *
 */
void
ncm_model_param_set_by_name (NcmModel *model, const gchar *param_name, gdouble val, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  {
    guint i;
    const gboolean has_param = ncm_model_param_index_from_name (model, param_name, &i, error);

    if (!has_param)
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
                                  "ncm_model_param_set_by_name: model `%s' does not have a parameter called `%s'. "
                                  "Use the method ncm_model_param_index_from_name() to check if the parameter exists.",
                                  G_OBJECT_TYPE_NAME (model), param_name);

      return;
    }

    ncm_model_param_set (model, i, val);
  }
}

/**
 * ncm_model_orig_param_set_by_name:
 * @model: a #NcmModel
 * @param_name: parameter name
 * @val: parameter value
 * @error: a #GError
 *
 * Sets the parameter value @val by @param_name.
 *
 */
void
ncm_model_orig_param_set_by_name (NcmModel *model, const gchar *param_name, gdouble val, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  {
    guint i;
    const gboolean has_param = ncm_model_orig_param_index_from_name (model, param_name, &i);

    if (!has_param)
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_ORIG_PARAM_NAME_NOT_FOUND,
                                  "ncm_model_orig_param_set_by_name: model `%s' does not have a parameter called `%s'. "
                                  "Use the method ncm_model_orig_param_index_from_name() to check if the parameter exists.",
                                  G_OBJECT_TYPE_NAME (model), param_name);

      return;
    }

    ncm_model_orig_param_set (model, i, val);
  }
}

/**
 * ncm_model_param_get_by_name:
 * @model: a #NcmModel
 * @param_name: parameter name
 * @error: a #GError
 *
 * Gets the parameter value by @param_name.
 *
 * Returns: parameter value.
 */
gdouble
ncm_model_param_get_by_name (NcmModel *model, const gchar *param_name, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, GSL_NAN);
  {
    guint i;
    const gboolean has_param = ncm_model_param_index_from_name (model, param_name, &i, error);

    if (!has_param)
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
                                  "ncm_model_param_get_by_name: model `%s' does not have a parameter called `%s'. "
                                  "Use the method ncm_model_param_index_from_name() to check if the parameter exists.",
                                  G_OBJECT_TYPE_NAME (model), param_name);

      return GSL_NAN;
    }

    return ncm_model_param_get (model, i);
  }
}

/**
 * ncm_model_orig_param_get_by_name:
 * @model: a #NcmModel
 * @param_name: parameter name
 * @error: a #GError
 *
 * Gets the original parameter value by @param_name.
 *
 * Returns: parameter value.
 */
gdouble
ncm_model_orig_param_get_by_name (NcmModel *model, const gchar *param_name, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, GSL_NAN);
  {
    guint i;
    const gboolean has_param = ncm_model_orig_param_index_from_name (model, param_name, &i);

    if (!has_param)
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_ORIG_PARAM_NAME_NOT_FOUND,
                                  "ncm_model_orig_param_get_by_name: model `%s' does not have a parameter called `%s'. "
                                  "Use the method ncm_model_orig_param_index_from_name() to check if the parameter exists.",
                                  G_OBJECT_TYPE_NAME (model), param_name);

      return GSL_NAN;
    }

    return ncm_model_orig_param_get (model, i);
  }
}

/**
 * ncm_model_type_is_submodel:
 * @model_type: a GType
 *
 * Tests if @model_type is a submodel of other model class.
 *
 * Returns: TRUE if @model_type is a submodel of other model class.
 */
gboolean
ncm_model_type_is_submodel (GType model_type)
{
  g_assert (g_type_is_a (model_type, NCM_TYPE_MODEL));
  {
    NcmModelClass *model_class = g_type_class_ref (model_type);
    gboolean is_submodel       = model_class->is_submodel;

    g_type_class_unref (model_class);

    return is_submodel;
  }
}

/**
 * ncm_model_type_main_model:
 * @model_type: a GType
 *
 * If @model_type is a submodel returns the #NcmModelID of its
 * main model, otherwise returns -1.
 *
 * Returns: main model #NcmModelID or -1.
 */
NcmModelID
ncm_model_type_main_model (GType model_type)
{
  g_assert (g_type_is_a (model_type, NCM_TYPE_MODEL));
  {
    NcmModelClass *model_class = g_type_class_ref (model_type);
    NcmModelID main_model_id   = model_class->main_model_id;

    g_type_class_unref (model_class);

    return main_model_id;
  }
}

/**
 * ncm_model_is_submodel:
 * @model: a #NcmModel
 *
 * Tests if @model is a submodel of other model class.
 *
 * Returns: TRUE if @model is a submodel of other model class.
 */
gboolean
ncm_model_is_submodel (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->is_submodel;
}

/**
 * ncm_model_main_model:
 * @model: a #NcmModel
 *
 * If @model is a submodel returns the #NcmModelID of its
 * main model, otherwise returns -1.
 *
 * Returns: main model #NcmModelID or -1.
 */
NcmModelID
ncm_model_main_model (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->main_model_id;
}

/**
 * ncm_model_add_submodel: (virtual add_submodel)
 * @model: a #NcmModel
 * @submodel: a #NcmModel
 *
 * Adds the @submodel to the @model submodels.
 *
 */
void
ncm_model_add_submodel (NcmModel *model, NcmModel *submodel)
{
  NCM_MODEL_GET_CLASS (model)->add_submodel (model, submodel);
}

static void
_ncm_model_add_submodel (NcmModel *model, NcmModel *submodel)
{
  NcmModelPrivate * const self   = ncm_model_get_instance_private (model);
  NcmModelClass *submodel_class  = NCM_MODEL_GET_CLASS (submodel);
  const NcmModelID main_model_id = submodel_class->main_model_id;
  const gboolean is_submodel     = submodel_class->is_submodel;
  const NcmModelID submodel_mid  = ncm_model_id (submodel);
  gpointer pos_ptr, orig_key;

  g_assert (is_submodel);
  g_assert_cmpint (main_model_id, ==, ncm_model_id (model));

  if (g_hash_table_lookup_extended (self->submodel_mid_pos, GINT_TO_POINTER (submodel_mid), &orig_key, &pos_ptr))
  {
    const gint pos = GPOINTER_TO_INT (pos_ptr);

    g_assert_cmpint (pos, >, -1);
    g_assert_cmpint (pos, <, self->submodel_array->len);
    {
      NcmModel *old_submodel = g_ptr_array_index (self->submodel_array, pos);

      g_ptr_array_index (self->submodel_array, pos) = ncm_model_ref (submodel);
      ncm_model_free (old_submodel);
    }
  }
  else
  {
    gint pos = self->submodel_array->len;

    ncm_model_ref (submodel);
    g_ptr_array_add (self->submodel_array, submodel);
    g_hash_table_insert (self->submodel_mid_pos, GINT_TO_POINTER (submodel_mid), GINT_TO_POINTER (pos));
  }
}

/**
 * ncm_model_get_submodel_len:
 * @model: a #NcmModel
 *
 * Gets the number of submodels set in @model.
 *
 * Returns: the number of submodels set in @model.
 */
guint
ncm_model_get_submodel_len (NcmModel *model)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  return self->submodel_array->len;
}

/**
 * ncm_model_peek_submodel:
 * @model: a #NcmModel
 * @i: submodel position
 *
 * Gets the @i-th submodel.
 *
 * Returns: (transfer none): a #NcmModel.
 */
NcmModel *
ncm_model_peek_submodel (NcmModel *model, guint i)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);

  g_assert_cmpuint (i, <, ncm_model_get_submodel_len (model));

  return g_ptr_array_index (self->submodel_array, i);
}

/**
 * ncm_model_peek_submodel_by_mid:
 * @model: a #NcmModel
 * @mid: a #NcmModelID
 *
 * Gets the submodel if type #NcmModelID @mid.
 *
 * Returns: (transfer none): a #NcmModel.
 */
NcmModel *
ncm_model_peek_submodel_by_mid (NcmModel *model, NcmModelID mid)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gpointer pos_ptr, orig_key;

  if (g_hash_table_lookup_extended (self->submodel_mid_pos, GINT_TO_POINTER (mid), &orig_key, &pos_ptr))
  {
    gint pos = GPOINTER_TO_INT (pos_ptr);

    g_assert_cmpint (pos, >, -1);
    g_assert_cmpint (pos, <, self->submodel_array->len);

    return g_ptr_array_index (self->submodel_array, pos);
  }
  else
  {
    return NULL;
  }
}

/**
 * ncm_model_peek_submodel_pos_by_mid:
 * @model: a #NcmModel
 * @mid: a #NcmModelID
 *
 * Gets the submodel type #NcmModelID @mid position.
 *
 * Returns: the @mid position or -1 if not found.
 */
gint
ncm_model_peek_submodel_pos_by_mid (NcmModel *model, NcmModelID mid)
{
  NcmModelPrivate * const self = ncm_model_get_instance_private (model);
  gpointer pos_ptr, orig_key;

  if (g_hash_table_lookup_extended (self->submodel_mid_pos, GINT_TO_POINTER (mid), &orig_key, &pos_ptr))
  {
    gint pos = GPOINTER_TO_INT (pos_ptr);

    g_assert_cmpint (pos, >, -1);
    g_assert_cmpint (pos, <, self->submodel_array->len);

    return pos;
  }
  else
  {
    return -1;
  }
}

/**
 * ncm_model___getitem__:
 * @model: a #NcmModel
 * @param: parameter name
 * @error: a GError
 *
 * Gets the parameter by name.
 *
 * Returns: parameter value
 */
gdouble
ncm_model___getitem__ (NcmModel *model, gchar *param, GError **error)
{
  g_return_val_if_fail (error == NULL || *error == NULL, GSL_NAN);
  {
    guint i;
    gboolean exists = ncm_model_param_index_from_name (model, param, &i, error);

    if (error && *error)
      return GSL_NAN;

    if (exists)
    {
      return ncm_model_param_get (model, i);
    }
    else
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
                                  "Parameter named: %s does not exist in %s",
                                  param, G_OBJECT_TYPE_NAME (model));
      NCM_UTIL_ON_ERROR_RETURN (error, , GSL_NAN);
    }

    return GSL_NAN;
  }
}

/**
 * ncm_model___setitem__:
 * @model: a #NcmModel
 * @param: parameter name
 * @val: parameter value
 * @error: a pointer for GError
 *
 * Sets the parameter by name.
 *
 */
void
ncm_model___setitem__ (NcmModel *model, gchar *param, gdouble val, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  {
    guint i;
    gboolean exists = ncm_model_param_index_from_name (model, param, &i, error);

    if (error && *error)
      return;

    if (exists)
    {
      ncm_model_param_set (model, i, val);
    }
    else
    {
      ncm_util_set_or_call_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND,
                                  "Parameter named: %s does not exist in %s",
                                  param, G_OBJECT_TYPE_NAME (model));
      NCM_UTIL_ON_ERROR_RETURN (error, , );
    }
  }
}


/***************************************************************************
 *            ncm_model_builder.c
 *
 *  Fri November 06 12:18:27 2015
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_builder.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmModelBuilder:
 *
 * A #NcmModel builder.
 *
 * This model can be used to create runtime NcmModels. It is particularly useful to
 * create models in binded languages, e.g., python.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_builder.h"
#include "math/ncm_model.h"
#include "math/ncm_mset.h"
#include "math/ncm_obj_array.h"

struct _NcmModelBuilder
{
  /*< private >*/
  GObject parent_instance;
  gchar *name;
  gchar *desc;
  gchar *parent_type_string;
  GType ptype;
  GType type;
  NcmObjArray *sparams;
  NcmObjArray *vparams;
  gboolean stackable;
  gboolean created;
};

enum
{
  PROP_0,
  PROP_PARENT_TYPE_STRING,
  PROP_NAME,
  PROP_DESC,
  PROP_SPARAMS,
  PROP_VPARAMS,
  PROP_STACKABLE,
};

G_DEFINE_TYPE (NcmModelBuilder, ncm_model_builder, G_TYPE_OBJECT)

static void
ncm_model_builder_init (NcmModelBuilder *mb)
{
  mb->name               = NULL;
  mb->desc               = NULL;
  mb->parent_type_string = NULL;
  mb->ptype              = G_TYPE_INVALID;
  mb->type               = G_TYPE_INVALID;
  mb->sparams            = ncm_obj_array_new ();
  mb->vparams            = ncm_obj_array_new ();
  mb->created            = FALSE;
}

static void
ncm_model_builder_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModelBuilder *mb = NCM_MODEL_BUILDER (object);

  g_return_if_fail (NCM_IS_MODEL_BUILDER (object));

  switch (prop_id)
  {
    case PROP_PARENT_TYPE_STRING:
      g_assert_false (mb->created);
      mb->parent_type_string = g_value_dup_string (value);
      mb->ptype              = g_type_from_name (mb->parent_type_string);
      break;
    case PROP_NAME:
      g_assert_false (mb->created);
      g_clear_pointer (&mb->name, g_free);
      mb->name = g_value_dup_string (value);
      break;
    case PROP_DESC:
      g_assert_false (mb->created);
      g_clear_pointer (&mb->desc, g_free);
      mb->desc = g_value_dup_string (value);
      break;
    case PROP_SPARAMS:
    {
      NcmObjArray *sparams = g_value_get_boxed (value);

      if (sparams)
      {
        g_assert_false (mb->created);
        g_clear_pointer (&mb->sparams, ncm_obj_array_unref);
        mb->sparams = ncm_obj_array_ref (sparams);
      }

      break;
    }
    case PROP_VPARAMS:
    {
      NcmObjArray *vparams = g_value_get_boxed (value);

      if (vparams)
      {
        g_assert_false (mb->created);
        g_clear_pointer (&mb->vparams, ncm_obj_array_unref);
        mb->vparams = ncm_obj_array_ref (vparams);
      }

      break;
    }
    case PROP_STACKABLE:
      g_assert_false (mb->created);
      mb->stackable = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_model_builder_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModelBuilder *mb = NCM_MODEL_BUILDER (object);

  g_return_if_fail (NCM_IS_MODEL_BUILDER (object));

  switch (prop_id)
  {
    case PROP_PARENT_TYPE_STRING:
      g_value_set_string (value, mb->parent_type_string);
      break;
    case PROP_NAME:
      g_value_set_string (value, mb->name);
      break;
    case PROP_DESC:
      g_value_set_string (value, mb->desc);
      break;
    case PROP_SPARAMS:
      g_value_set_boxed (value, mb->sparams);
      break;
    case PROP_VPARAMS:
      g_value_set_boxed (value, mb->vparams);
      break;
    case PROP_STACKABLE:
      g_value_set_boolean (value, mb->stackable);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_model_builder_dispose (GObject *object)
{
  NcmModelBuilder *mb = NCM_MODEL_BUILDER (object);

  g_ptr_array_set_size (mb->sparams, 0);
  g_ptr_array_set_size (mb->vparams, 0);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_builder_parent_class)->dispose (object);
}

static void
ncm_model_builder_finalize (GObject *object)
{
  NcmModelBuilder *mb = NCM_MODEL_BUILDER (object);

  g_clear_pointer (&mb->sparams, ncm_obj_array_unref);
  g_clear_pointer (&mb->vparams, ncm_obj_array_unref);
  g_clear_pointer (&mb->name, g_free);
  g_clear_pointer (&mb->desc, g_free);
  g_clear_pointer (&mb->parent_type_string, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_builder_parent_class)->finalize (object);
}

static void
ncm_model_builder_class_init (NcmModelBuilderClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = ncm_model_builder_set_property;
  object_class->get_property = ncm_model_builder_get_property;
  object_class->dispose      = ncm_model_builder_dispose;
  object_class->finalize     = ncm_model_builder_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PARENT_TYPE_STRING,
                                   g_param_spec_string ("parent-type-string",
                                                        NULL,
                                                        "Parent type name",
                                                        "NcmModel",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Model's name",
                                                        "no-name",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DESC,
                                   g_param_spec_string ("description",
                                                        NULL,
                                                        "Model's description",
                                                        "no-description",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SPARAMS,
                                   g_param_spec_boxed ("sparams",
                                                       NULL,
                                                       "Scalar parameters",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_VPARAMS,
                                   g_param_spec_boxed ("vparams",
                                                       NULL,
                                                       "Vector parameters",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_STACKABLE,
                                   g_param_spec_boolean ("stackable",
                                                         NULL,
                                                         "Stackable",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_model_builder_new:
 * @ptype: Parent's type
 * @name: model's name
 * @desc: model's description
 *
 * Creates a new NcmModelBuilder object. This does not create the new class
 * after defining all parameters one should call ncm_model_builder_create()
 * to effectively define a new class.
 *
 * Returns: (transfer full): a new #NcmModelBuilder.
 */
NcmModelBuilder *
ncm_model_builder_new (GType ptype, const gchar *name, const gchar *desc)
{
  const gchar *parent_type_string = g_type_name (ptype);
  NcmModelBuilder *mb             = g_object_new (NCM_TYPE_MODEL_BUILDER,
                                                  "parent-type-string", parent_type_string,
                                                  "name",               name,
                                                  "description",        desc,
                                                  NULL);

  return mb;
}

/**
 * ncm_model_builder_ref:
 * @mb: a #NcmModelBuilder
 *
 * Increase reference count of @mb by one.
 *
 * Returns: (transfer full): @mb.
 */
NcmModelBuilder *
ncm_model_builder_ref (NcmModelBuilder *mb)
{
  return g_object_ref (mb);
}

/**
 * ncm_model_builder_add_sparam_obj:
 * @mb: a #NcmModelBuilder
 * @sparam: a #NcmSParam
 *
 * Adds the parameters described by @sparam to @mb.
 *
 */
void
ncm_model_builder_add_sparam_obj (NcmModelBuilder *mb, NcmSParam *sparam)
{
  g_assert (!mb->created);
  ncm_obj_array_add (mb->sparams, G_OBJECT (sparam));
}

/**
 * ncm_model_builder_add_vparam_obj:
 * @mb: a #NcmModelBuilder
 * @vparam: a #NcmVParam
 *
 * Adds the parameters described by @sparam to @mb.
 *
 */
void
ncm_model_builder_add_vparam_obj (NcmModelBuilder *mb, NcmVParam *vparam)
{
  g_assert (!mb->created);
  ncm_obj_array_add (mb->vparams, G_OBJECT (vparam));
}

/**
 * ncm_model_builder_add_sparam:
 * @mb: a #NcmModelBuilder
 * @symbol: symbol of the scalar parameter
 * @name: name of the sacalar parameter
 * @lower_bound: lower-bound value
 * @upper_bound: upper-bound value
 * @scale: parameter scale
 * @abstol: absolute tolerance
 * @default_value: default value
 * @ppt: a #NcmParamType
 *
 * Creates a new #NcmSParam from arguments and add it to @mb.
 *
 */
void
ncm_model_builder_add_sparam (NcmModelBuilder *mb, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt)
{
  NcmSParam *sparam = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_value, ppt);

  ncm_model_builder_add_sparam_obj (mb, sparam);

  ncm_sparam_free (sparam);
}

/**
 * ncm_model_builder_add_vparam:
 * @mb: a #NcmModelBuilder
 * @default_length: default length of the vector parameter
 * @symbol: symbol of the vector parameter
 * @name: name of the vector parameter
 * @lower_bound: parameter lower bound
 * @upper_bound: parameter upper bound
 * @scale: parameter scale
 * @abstol: absolute tolerance
 * @default_value: default value
 * @ppt: a #NcmParamType
 *
 * Creates a new #NcmVParam from arguments and add it to @mb.
 *
 */
void
ncm_model_builder_add_vparam (NcmModelBuilder *mb, guint default_length, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt)
{
  NcmVParam *vparam = ncm_vparam_full_new (default_length, name, symbol, lower_bound, upper_bound, scale, abstol, default_value, ppt);

  ncm_model_builder_add_vparam_obj (mb, vparam);

  ncm_vparam_free (vparam);
}

/**
 * ncm_model_builder_add_sparams:
 * @mb: a #NcmModelBuilder
 * @sparams: an array of #NcmSParam objects
 *
 * Adds all #NcmSParam objects in @sparams to @mb.
 *
 */
void
ncm_model_builder_add_sparams (NcmModelBuilder *mb, NcmObjArray *sparams)
{
  guint i;

  for (i = 0; i < sparams->len; i++)
  {
    NcmSParam *sparam = NCM_SPARAM (ncm_obj_array_get (sparams, i));

    g_assert (NCM_IS_SPARAM (sparam));
    ncm_model_builder_add_sparam_obj (mb, sparam);
  }
}

/**
 * ncm_model_builder_get_sparams:
 * @mb: a #NcmModelBuilder
 *
 * Returns: (transfer full): a #NcmObjArray containing all #NcmSParam objects in @mb.
 */
NcmObjArray *
ncm_model_builder_get_sparams (NcmModelBuilder *mb)
{
  NcmObjArray *oa = ncm_obj_array_new ();
  guint i;

  for (i = 0; i < mb->sparams->len; i++)
  {
    GObject *obj = ncm_obj_array_peek (mb->sparams, i);

    ncm_obj_array_add (oa, obj);
  }

  return oa;
}

void
_ncm_model_builder_class_init (gpointer g_class, gpointer class_data)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (g_class);
  NcmModelBuilder *mb        = NCM_MODEL_BUILDER (class_data);
  guint i;

  if (model_class->main_model_id == -1)
    ncm_mset_model_register_id (model_class, mb->name, mb->desc, NULL, mb->stackable, -1);

  ncm_model_class_set_name_nick (model_class, mb->name, mb->name);
  ncm_model_class_add_params (model_class, mb->sparams->len, mb->vparams->len, 1);

  for (i = 0; i < mb->sparams->len; i++)
  {
    NcmSParam *sparam = NCM_SPARAM (ncm_obj_array_peek (mb->sparams, i));

    ncm_model_class_set_sparam_obj (model_class, i, sparam);
  }

  for (i = 0; i < mb->vparams->len; i++)
  {
    NcmVParam *vparam = NCM_VPARAM (ncm_obj_array_peek (mb->vparams, i));

    ncm_model_class_set_vparam_obj (model_class, i, vparam);
  }

  ncm_model_class_check_params_info (model_class);
}

void
_ncm_model_builder_instance_init (GTypeInstance *instance, gpointer g_class)
{
}

/**
 * ncm_model_builder_create:
 * @mb: a #NcmModelBuilder
 *
 * Creates a new object type using the scalar and vector parameters defined
 * in @mb. If the object type was already created, this function just returns
 * the type.
 *
 * Returns: the new object type.
 */
GType
ncm_model_builder_create (NcmModelBuilder *mb)
{
  if (!mb->created)
  {
    GTypeQuery query = {0, };
    GTypeInfo info   = {0, };

    g_type_query (mb->ptype, &query);

    info.class_size     = query.class_size;
    info.base_init      = NULL;
    info.base_finalize  = NULL;
    info.class_init     = _ncm_model_builder_class_init;
    info.class_finalize = NULL;
    info.class_data     = ncm_model_builder_ref (mb);
    info.instance_size  = query.instance_size;
    info.n_preallocs    = 0;
    info.instance_init  = _ncm_model_builder_instance_init;
    info.value_table    = NULL;

    mb->type    = g_type_register_static (mb->ptype, mb->name, &info, 0);
    mb->created = TRUE;

    {
      GObjectClass *object_class = g_type_class_ref (mb->type);

      g_type_class_unref (object_class);
    }
  }

  return mb->type;
}


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
 * Builds a #NcmModel subclass at run time.
 *
 * A builder holds the parent type, the name and description of the new model and
 * its scalar and vector parameters; ncm_model_builder_create() registers the new
 * type. The builder is serializable, so a model defined in a language binding, for
 * instance Python, can be saved and rebuilt in another process.
 *
 * The parent must be #NcmModel or one of its subclasses. The new parameters come
 * after those of the parent. If the parent is not registered in #NcmMSet, the new
 * model is registered under its name; otherwise it shares the parent's model id,
 * as C subclasses do.
 *
 * Creating a type whose name is already registered returns the existing type if it
 * has the same parent and the same parameters, so the same builder may be created
 * again in one process; any other name clash aborts.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/model/ncm_model_builder.h"
#include "ncm/model/ncm_model.h"
#include "ncm/model/ncm_mset.h"
#include "ncm/core/ncm_obj_array.h"

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
  mb->stackable          = FALSE;
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
      mb->parent_type_string = g_value_dup_string (value);
      mb->ptype              = g_type_from_name (mb->parent_type_string);

      if (mb->ptype == 0)
        g_error ("ncm_model_builder_set_property: parent type `%s' is not registered.",
                 mb->parent_type_string);

      if (!g_type_is_a (mb->ptype, NCM_TYPE_MODEL))
        g_error ("ncm_model_builder_set_property: parent type `%s' is not a NcmModel.",
                 mb->parent_type_string);

      break;
    case PROP_NAME:
      mb->name = g_value_dup_string (value);
      break;
    case PROP_DESC:
      mb->desc = g_value_dup_string (value);
      break;
    case PROP_SPARAMS:
    {
      NcmObjArray *sparams = g_value_get_boxed (value);

      if (sparams != NULL)
      {
        g_clear_pointer (&mb->sparams, ncm_obj_array_unref);
        mb->sparams = ncm_obj_array_ref (sparams);
      }

      break;
    }
    case PROP_VPARAMS:
    {
      NcmObjArray *vparams = g_value_get_boxed (value);

      if (vparams != NULL)
      {
        g_clear_pointer (&mb->vparams, ncm_obj_array_unref);
        mb->vparams = ncm_obj_array_ref (vparams);
      }

      break;
    }
    case PROP_STACKABLE:
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

  /**
   * NcmModelBuilder:parent-type-string:
   *
   * The name of the parent type, #NcmModel or one of its subclasses.
   */
  g_object_class_install_property (object_class,
                                   PROP_PARENT_TYPE_STRING,
                                   g_param_spec_string ("parent-type-string",
                                                        NULL,
                                                        "Parent type name",
                                                        "NcmModel",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmModelBuilder:name:
   *
   * The name of the new type, also its #NcmMSet namespace when the parent is not
   * registered.
   */
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Model's name",
                                                        "no-name",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmModelBuilder:description:
   *
   * The description of the new model.
   */
  g_object_class_install_property (object_class,
                                   PROP_DESC,
                                   g_param_spec_string ("description",
                                                        NULL,
                                                        "Model's description",
                                                        "no-description",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmModelBuilder:sparams:
   *
   * The scalar parameters of the new model, a #NcmObjArray of #NcmSParam.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPARAMS,
                                   g_param_spec_boxed ("sparams",
                                                       NULL,
                                                       "Scalar parameters",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmModelBuilder:vparams:
   *
   * The vector parameters of the new model, a #NcmObjArray of #NcmVParam.
   */
  g_object_class_install_property (object_class,
                                   PROP_VPARAMS,
                                   g_param_spec_boxed ("vparams",
                                                       NULL,
                                                       "Vector parameters",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmModelBuilder:stackable:
   *
   * Whether several instances of the new model can be stacked in a #NcmMSet. Used
   * only when the new model is registered under its own name.
   */
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
 * @ptype: parent type, #NcmModel or one of its subclasses
 * @name: name of the new type
 * @desc: description of the new model
 *
 * Creates a new #NcmModelBuilder. Add the parameters and then call
 * ncm_model_builder_create() to register the new type.
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
 * Increases the reference count of @mb by one.
 *
 * Returns: (transfer full): @mb.
 */
NcmModelBuilder *
ncm_model_builder_ref (NcmModelBuilder *mb)
{
  return g_object_ref (mb);
}

/**
 * ncm_model_builder_free:
 * @mb: a #NcmModelBuilder
 *
 * Decreases the reference count of @mb by one. If the reference count reaches
 * zero, @mb is freed.
 *
 */
void
ncm_model_builder_free (NcmModelBuilder *mb)
{
  g_object_unref (mb);
}

/**
 * ncm_model_builder_clear:
 * @mb: a #NcmModelBuilder
 *
 * If *@mb is not %NULL, decreases the reference count of *@mb by one and sets
 * *@mb to %NULL.
 *
 */
void
ncm_model_builder_clear (NcmModelBuilder **mb)
{
  g_clear_object (mb);
}

/**
 * ncm_model_builder_add_sparam_obj:
 * @mb: a #NcmModelBuilder
 * @sparam: a #NcmSParam
 *
 * Adds the scalar parameter @sparam to @mb. Aborts if the type was already
 * created.
 *
 */
void
ncm_model_builder_add_sparam_obj (NcmModelBuilder *mb, NcmSParam *sparam)
{
  if (mb->created)
    g_error ("ncm_model_builder_add_sparam_obj: model `%s' was already created, "
             "cannot add the parameter `%s'.",
             mb->name, ncm_sparam_name (sparam));

  ncm_obj_array_add (mb->sparams, G_OBJECT (sparam));
}

/**
 * ncm_model_builder_add_vparam_obj:
 * @mb: a #NcmModelBuilder
 * @vparam: a #NcmVParam
 *
 * Adds the vector parameter @vparam to @mb. Aborts if the type was already
 * created.
 *
 */
void
ncm_model_builder_add_vparam_obj (NcmModelBuilder *mb, NcmVParam *vparam)
{
  if (mb->created)
    g_error ("ncm_model_builder_add_vparam_obj: model `%s' was already created, "
             "cannot add the parameter `%s'.",
             mb->name, ncm_vparam_name (vparam));

  ncm_obj_array_add (mb->vparams, G_OBJECT (vparam));
}

/**
 * ncm_model_builder_add_sparam:
 * @mb: a #NcmModelBuilder
 * @symbol: symbol of the scalar parameter
 * @name: name of the scalar parameter
 * @lower_bound: lower bound
 * @upper_bound: upper bound
 * @scale: parameter scale
 * @abstol: absolute tolerance
 * @default_value: default value
 * @ppt: a #NcmParamType
 *
 * Creates a new #NcmSParam from the arguments and adds it to @mb.
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
 * @lower_bound: lower bound
 * @upper_bound: upper bound
 * @scale: parameter scale
 * @abstol: absolute tolerance
 * @default_value: default value
 * @ppt: a #NcmParamType
 *
 * Creates a new #NcmVParam from the arguments and adds it to @mb.
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
 * @sparams: a #NcmObjArray of #NcmSParam
 *
 * Adds every #NcmSParam in @sparams to @mb.
 *
 */
void
ncm_model_builder_add_sparams (NcmModelBuilder *mb, NcmObjArray *sparams)
{
  guint i;

  for (i = 0; i < sparams->len; i++)
  {
    GObject *obj = ncm_obj_array_peek (sparams, i);

    if (!NCM_IS_SPARAM (obj))
      g_error ("ncm_model_builder_add_sparams: element %u is a `%s', not a NcmSParam.", i, G_OBJECT_TYPE_NAME (obj));

    ncm_model_builder_add_sparam_obj (mb, NCM_SPARAM (obj));
  }
}

/**
 * ncm_model_builder_add_vparams:
 * @mb: a #NcmModelBuilder
 * @vparams: a #NcmObjArray of #NcmVParam
 *
 * Adds every #NcmVParam in @vparams to @mb.
 *
 */
void
ncm_model_builder_add_vparams (NcmModelBuilder *mb, NcmObjArray *vparams)
{
  guint i;

  for (i = 0; i < vparams->len; i++)
  {
    GObject *obj = ncm_obj_array_peek (vparams, i);

    if (!NCM_IS_VPARAM (obj))
      g_error ("ncm_model_builder_add_vparams: element %u is a `%s', not a NcmVParam.", i, G_OBJECT_TYPE_NAME (obj));

    ncm_model_builder_add_vparam_obj (mb, NCM_VPARAM (obj));
  }
}

/**
 * ncm_model_builder_get_sparams:
 * @mb: a #NcmModelBuilder
 *
 * Gets the scalar parameters of @mb.
 *
 * Returns: (transfer full): a new #NcmObjArray with the #NcmSParam objects in @mb.
 */
NcmObjArray *
ncm_model_builder_get_sparams (NcmModelBuilder *mb)
{
  NcmObjArray *oa = ncm_obj_array_new ();
  guint i;

  for (i = 0; i < mb->sparams->len; i++)
    ncm_obj_array_add (oa, ncm_obj_array_peek (mb->sparams, i));

  return oa;
}

/**
 * ncm_model_builder_get_vparams:
 * @mb: a #NcmModelBuilder
 *
 * Gets the vector parameters of @mb.
 *
 * Returns: (transfer full): a new #NcmObjArray with the #NcmVParam objects in @mb.
 */
NcmObjArray *
ncm_model_builder_get_vparams (NcmModelBuilder *mb)
{
  NcmObjArray *oa = ncm_obj_array_new ();
  guint i;

  for (i = 0; i < mb->vparams->len; i++)
    ncm_obj_array_add (oa, ncm_obj_array_peek (mb->vparams, i));

  return oa;
}

static void _ncm_model_builder_class_init (gpointer g_class, gpointer class_data);
static void _ncm_model_builder_check_existing (NcmModelBuilder *mb, GType existing);

/**
 * ncm_model_builder_create:
 * @mb: a #NcmModelBuilder
 *
 * Registers the new type with the parameters of @mb and returns it. Later calls
 * return the same type. If a type with the name of @mb is already registered, it
 * is returned when it has the same parent and the same parameters; otherwise this
 * function aborts.
 *
 * Returns: the new type.
 */
GType
ncm_model_builder_create (NcmModelBuilder *mb)
{
  if (!mb->created)
  {
    const GType existing = g_type_from_name (mb->name);

    if (existing != 0)
    {
      _ncm_model_builder_check_existing (mb, existing);
      mb->type = existing;
    }
    else
    {
      GTypeQuery query = {0, };
      GTypeInfo info   = {0, };

      g_type_query (mb->ptype, &query);

      info.class_size    = query.class_size;
      info.class_init    = _ncm_model_builder_class_init;
      info.class_data    = ncm_model_builder_ref (mb);
      info.instance_size = query.instance_size;

      mb->type = g_type_register_static (mb->ptype, mb->name, &info, 0);

      g_type_class_unref (g_type_class_ref (mb->type));
    }

    mb->created = TRUE;
  }

  return mb->type;
}

/*
 * The class_init of the new type: the parameters of @mb come after those of the
 * parent, whose ids are absolute.
 */
static void
_ncm_model_builder_class_init (gpointer g_class, gpointer class_data)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (g_class);
  NcmModelBuilder *mb        = NCM_MODEL_BUILDER (class_data);
  guint i;

  if (model_class->model_id < 0)
    ncm_mset_model_register_id (model_class, mb->name, mb->desc, NULL, mb->stackable, -1);

  ncm_model_class_set_name_nick (model_class, mb->name, mb->name);
  ncm_model_class_add_params (model_class, mb->sparams->len, mb->vparams->len, 1);

  for (i = 0; i < mb->sparams->len; i++)
  {
    NcmSParam *sparam = NCM_SPARAM (ncm_obj_array_peek (mb->sparams, i));

    ncm_model_class_set_sparam_obj (model_class, model_class->parent_sparam_len + i, sparam);
  }

  for (i = 0; i < mb->vparams->len; i++)
  {
    NcmVParam *vparam = NCM_VPARAM (ncm_obj_array_peek (mb->vparams, i));

    ncm_model_class_set_vparam_obj (model_class, model_class->parent_vparam_len + i, vparam);
  }

  ncm_model_class_check_params_info (model_class);
}

static gboolean _ncm_model_builder_sparam_equal (const NcmSParam *a, const NcmSParam *b);

/*
 * Aborts unless @existing, a type already registered under the name of @mb, has
 * the parent and the parameters @mb would give it.
 */
static void
_ncm_model_builder_check_existing (NcmModelBuilder *mb, GType existing)
{
  NcmModelClass *model_class;
  guint i;

  if (g_type_parent (existing) != mb->ptype)
    g_error ("ncm_model_builder_create: a type named `%s' already exists with parent `%s', "
             "but this builder has parent `%s'.",
             mb->name, g_type_name (g_type_parent (existing)), g_type_name (mb->ptype));

  model_class = g_type_class_ref (existing);

  if ((model_class->sparam_len - model_class->parent_sparam_len != mb->sparams->len) ||
      (model_class->vparam_len - model_class->parent_vparam_len != mb->vparams->len))
    g_error ("ncm_model_builder_create: a type named `%s' already exists with %u scalar and %u vector "
             "parameter(s), but this builder has %u and %u.",
             mb->name,
             model_class->sparam_len - model_class->parent_sparam_len,
             model_class->vparam_len - model_class->parent_vparam_len,
             mb->sparams->len, mb->vparams->len);

  for (i = 0; i < mb->sparams->len; i++)
  {
    const NcmSParam *a = g_ptr_array_index (model_class->sparam, model_class->parent_sparam_len + i);
    const NcmSParam *b = NCM_SPARAM (ncm_obj_array_peek (mb->sparams, i));

    if (!_ncm_model_builder_sparam_equal (a, b))
      g_error ("ncm_model_builder_create: a type named `%s' already exists, but its scalar parameter %u "
               "(`%s' there, `%s' here) differs in name, symbol, bounds, scale, tolerance, default value "
               "or fit type.",
               mb->name, i, ncm_sparam_name (a), ncm_sparam_name (b));
  }

  for (i = 0; i < mb->vparams->len; i++)
  {
    NcmVParam *a   = g_ptr_array_index (model_class->vparam, model_class->parent_vparam_len + i);
    NcmVParam *b   = NCM_VPARAM (ncm_obj_array_peek (mb->vparams, i));
    gboolean equal = (g_strcmp0 (ncm_vparam_name (a), ncm_vparam_name (b)) == 0) &&
                     (g_strcmp0 (ncm_vparam_symbol (a), ncm_vparam_symbol (b)) == 0) &&
                     (ncm_vparam_len (a) == ncm_vparam_len (b));
    guint j;

    for (j = 0; equal && (j < ncm_vparam_len (a)); j++)
      equal = _ncm_model_builder_sparam_equal (ncm_vparam_peek_sparam (a, j), ncm_vparam_peek_sparam (b, j));

    if (!equal)
      g_error ("ncm_model_builder_create: a type named `%s' already exists, but its vector parameter %u "
               "(`%s' there, `%s' here) differs in name, symbol, length or components.",
               mb->name, i, ncm_vparam_name (a), ncm_vparam_name (b));
  }

  g_type_class_unref (model_class);
}

static gboolean
_ncm_model_builder_sparam_equal (const NcmSParam *a, const NcmSParam *b)
{
  return (g_strcmp0 (ncm_sparam_name (a), ncm_sparam_name (b)) == 0) &&
         (g_strcmp0 (ncm_sparam_symbol (a), ncm_sparam_symbol (b)) == 0) &&
         (ncm_sparam_get_lower_bound (a) == ncm_sparam_get_lower_bound (b)) &&
         (ncm_sparam_get_upper_bound (a) == ncm_sparam_get_upper_bound (b)) &&
         (ncm_sparam_get_scale (a) == ncm_sparam_get_scale (b)) &&
         (ncm_sparam_get_absolute_tolerance (a) == ncm_sparam_get_absolute_tolerance (b)) &&
         (ncm_sparam_get_default_value (a) == ncm_sparam_get_default_value (b)) &&
         (ncm_sparam_get_fit_type (a) == ncm_sparam_get_fit_type (b));
}


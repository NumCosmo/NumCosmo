/***************************************************************************
 *            ncm_reparam.c
 *
 *  Thu March 08 00:36:24 2012
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
 * SECTION:ncm_reparam
 * @title: NcmReparam
 * @short_description: Abstract class for model reparametrization.
 *
 * #NcmReparam is an abstract class for model reparametrization.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_reparam.h"
#include "math/ncm_model.h"
#include "math/ncm_cfg.h"
#include "math/ncm_sparam.h"
#include "math/ncm_vector.h"
#include "math/ncm_obj_array.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmReparamPrivate
{
  /*< private >*/
  GObject parent_instance;
  guint length;
  NcmVector *new_params;
  NcmObjDictInt *sparams;
  GHashTable *sparams_name_id;
  GType compat_type;
} NcmReparamPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmReparam, ncm_reparam, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_LEN,
  PROP_PARAMS_DESC,
  PROP_COMPAT_TYPE,
};

static void
ncm_reparam_init (NcmReparam *reparam)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  self->length          = 0;
  self->sparams         = ncm_obj_dict_int_new ();
  self->new_params      = NULL;
  self->sparams_name_id = g_hash_table_new_full (&g_str_hash, &g_str_equal, &g_free, NULL);
  self->compat_type     = G_TYPE_INVALID;
}

static void
_ncm_reparam_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmReparam *reparam            = NCM_REPARAM (object);
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  g_return_if_fail (NCM_IS_REPARAM (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, self->length);
      break;
    case PROP_PARAMS_DESC:
      g_value_set_boxed (value, self->sparams);
      break;
    case PROP_COMPAT_TYPE:
      g_value_set_string (value, g_type_name (self->compat_type));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_reparam_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmReparam *reparam            = NCM_REPARAM (object);
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  g_return_if_fail (NCM_IS_REPARAM (object));

  switch (prop_id)
  {
    case PROP_LEN:
      self->length = g_value_get_uint (value);
      break;
    case PROP_PARAMS_DESC:
      ncm_obj_dict_int_clear (&self->sparams);
      self->sparams = g_value_dup_boxed (value);
      break;
    case PROP_COMPAT_TYPE:
      self->compat_type = g_type_from_name (g_value_get_string (value));

      if (self->compat_type == G_TYPE_INVALID)
        g_error ("_ncm_reparam_set_property: GType `%s' unregistered or invalid.", g_value_get_string (value));

      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_reparam_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_reparam_parent_class)->constructed (object);
  {
    NcmReparam *reparam            = NCM_REPARAM (object);
    NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

    g_assert_cmpuint (self->length, >, 0);

    self->new_params = ncm_vector_new (self->length);
  }
}

static void
_ncm_reparam_finalize (GObject *object)
{
  NcmReparam *reparam            = NCM_REPARAM (object);
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  g_clear_pointer (&self->sparams_name_id, g_hash_table_unref);
  g_clear_pointer (&self->sparams, ncm_obj_dict_int_unref);
  ncm_vector_clear (&self->new_params);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_reparam_parent_class)->finalize (object);
}

static void
ncm_reparam_class_init (NcmReparamClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_reparam_set_property;
  object_class->get_property = &_ncm_reparam_get_property;
  object_class->constructed  = &_ncm_reparam_constructed;
  object_class->finalize     = &_ncm_reparam_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "System's length",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PARAMS_DESC,
                                   g_param_spec_boxed ("params-desc",
                                                       NULL,
                                                       "News parameter descriptions",
                                                       NCM_TYPE_OBJ_DICT_INT,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COMPAT_TYPE,
                                   g_param_spec_string ("compat-type",
                                                        NULL,
                                                        "Compatible type",
                                                        g_type_name (NCM_TYPE_MODEL),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->old2new = NULL;
  klass->new2old = NULL;
}

/**
 * ncm_reparam_ref:
 * @reparam: a #NcmReparam
 *
 * Increases the reference count of @reparam by one.
 *
 * Returns: (transfer full): the passed #NcmReparam.
 */
NcmReparam *
ncm_reparam_ref (NcmReparam *reparam)
{
  return NCM_REPARAM (g_object_ref (reparam));
}

/**
 * ncm_reparam_free:
 * @reparam: a #NcmReparam
 *
 * Decreases the reference count of @reparam by one. If the reference count
 * reaches zero, the #NcmReparam is freed.
 *
 */
void
ncm_reparam_free (NcmReparam *reparam)
{
  g_object_unref (reparam);
}

/**
 * ncm_reparam_clear:
 * @reparam: a #NcmReparam
 *
 * If *@reparam is not NULL, unrefs it and sets *@reparam to NULL.
 * If *@reparam is NULL, does nothing.
 *
 */
void
ncm_reparam_clear (NcmReparam **reparam)
{
  g_clear_object (reparam);
}

/**
 * ncm_reparam_set_compat_type:
 * @reparam: a #NcmReparam
 * @compat_type: a #GType
 *
 * Sets the compatible #GType for this reparametrization.
 *
 */
void
ncm_reparam_set_compat_type (NcmReparam *reparam, GType compat_type)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  self->compat_type = compat_type;
}

/**
 * ncm_reparam_get_compat_type:
 * @reparam: a #NcmReparam
 *
 * Gets the compatible #GType for this reparametrization.
 *
 * Returns: the compatible #GType for this reparametrization.
 */
GType
ncm_reparam_get_compat_type (NcmReparam *reparam)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  return self->compat_type;
}

/**
 * ncm_reparam_old2new: (virtual old2new)
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 *
 * Using the values set in the original parametrization update the
 * values of the new parametrization.
 *
 */
void
ncm_reparam_old2new (NcmReparam *reparam, NcmModel *model)
{
  NCM_REPARAM_GET_CLASS (reparam)->old2new (reparam, model);
}

/**
 * ncm_reparam_new2old: (virtual new2old)
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 *
 * Using the values set in the new parametrization update the
 * values of the original parametrization.
 *
 */
void
ncm_reparam_new2old (NcmReparam *reparam, NcmModel *model)
{
  NCM_REPARAM_GET_CLASS (reparam)->new2old (reparam, model);
}

/**
 * ncm_reparam_set_param_desc:
 * @reparam: a #NcmReparam
 * @i: index of the changed parameter.
 * @sp: #NcmSParam describing the new parameter.
 *
 * Change the @i-th parameter description using @sp.
 *
 */
void
ncm_reparam_set_param_desc (NcmReparam *reparam, guint i, NcmSParam *sp)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);
  NcmSParam *old_sp              = NCM_SPARAM (ncm_obj_dict_int_peek (self->sparams, i));

  g_assert (i < self->length);

  if (old_sp != NULL)
  {
    g_assert (g_hash_table_remove (self->sparams_name_id, ncm_sparam_name (old_sp)));
    ncm_sparam_clear (&old_sp);
  }

  g_hash_table_insert (self->sparams_name_id,
                       g_strdup (ncm_sparam_name (sp)),
                       GUINT_TO_POINTER (i));

  ncm_obj_dict_int_add (self->sparams, i, G_OBJECT (sp));
}

/**
 * ncm_reparam_peek_param_desc:
 * @reparam: a #NcmReparam
 * @i: index of the changed parameter.
 *
 * Peeks the @i-th parameter description.
 *
 * Returns: (transfer none): The @i-th parameter description.
 */
NcmSParam *
ncm_reparam_peek_param_desc (NcmReparam *reparam, guint i)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  g_assert_cmpuint (i, <, self->length);

  {
    GObject *sp = ncm_obj_dict_int_peek (self->sparams, i);

    if (sp != NULL)
      return NCM_SPARAM (sp);
    else
      return NULL;
  }
}

/**
 * ncm_reparam_get_param_desc:
 * @reparam: a #NcmReparam
 * @i: index of the changed parameter.
 *
 * Gets the @i-th parameter description.
 *
 * Returns: (transfer full): The @i-th parameter description.
 */
NcmSParam *
ncm_reparam_get_param_desc (NcmReparam *reparam, guint i)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  g_assert_cmpuint (i, <, self->length);

  {
    GObject *sp = ncm_obj_dict_int_peek (self->sparams, i);

    if (sp != NULL)
      return ncm_sparam_ref (NCM_SPARAM (sp));
    else
      return NULL;
  }
}

/**
 * ncm_reparam_set_param_desc_full:
 * @reparam: a #NcmReparam
 * @i: index of the changed parameter.
 * @name: #NcmSParam:name.
 * @symbol: #NcmSParam:symbol.
 * @lower_bound: value of #NcmSParam:lower-bound.
 * @upper_bound: value of #NcmSParam:upper-bound.
 * @scale: value of #NcmSParam:scale.
 * @abstol: value of #NcmSParam:absolute-tolerance.
 * @default_val: value of #NcmSParam:default-value.
 * @ftype: a #NcmParamType.
 *
 * Change the @i-th parameter description using the given values.
 *
 */
void
ncm_reparam_set_param_desc_full (NcmReparam *reparam, guint i, const gchar *name, const gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype)
{
  NcmSParam *sp = ncm_sparam_new (name, symbol, lower_bound, upper_bound,
                                  scale, abstol, default_val, ftype);

  ncm_reparam_set_param_desc (reparam, i, sp);


  ncm_sparam_free (sp);
}

/**
 * ncm_reparam_index_from_name:
 * @reparam: a #NcmReparam.
 * @param_name: parameter name.
 * @i: (out): parameter index.
 *
 * Looks for a parameter named @param_name and returns TRUE if found. If found
 * puts at @i its index.
 *
 * Returns: whenever the parameter is found.
 */
gboolean
ncm_reparam_index_from_name (NcmReparam *reparam, const gchar *param_name, guint *i)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);
  gpointer param_id;
  gboolean found = g_hash_table_lookup_extended (self->sparams_name_id, param_name, NULL, &param_id);

  if (found)
    *i = GPOINTER_TO_UINT (param_id);
  else
    *i = -1;  /* Yes, I known. */

  return found;
}

/**
 * ncm_reparam_get_length:
 * @reparam: a #NcmReparam
 *
 * Gets the number of parameters.
 *
 */
guint
ncm_reparam_get_length (NcmReparam *reparam)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  return self->length;
}

/**
 * ncm_reparam_peek_params:
 * @reparam: a #NcmReparam.
 *
 * Gets the #NcmVector containing the new parameters. This vector is owned by
 * the #NcmReparam and should not be freed. This method is used by #NcmModel
 * and subclasses to get the new parameters and should not be used by the
 * user. The pointer returned by this method is guaranteed to be valid until
 * the destruction of the #NcmReparam.
 *
 * Returns: (transfer none): a #NcmVector containing the new parameters.
 */
NcmVector *
ncm_reparam_peek_params (NcmReparam *reparam)
{
  NcmReparamPrivate * const self = ncm_reparam_get_instance_private (reparam);

  return self->new_params;
}


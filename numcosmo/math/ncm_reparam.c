/***************************************************************************
 *            ncm_reparam.c
 *
 *  Thu March 08 00:36:24 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * @title: Reparametrization Abstract Class
 * @short_description: Base class for model reparametrization
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_reparam.h"
#include "math/ncm_model.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_blas.h>

G_DEFINE_ABSTRACT_TYPE (NcmReparam, ncm_reparam, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_LEN,
  PROP_PARAMS_DESC,
};

static void
ncm_reparam_init (NcmReparam *reparam)
{
  reparam->length          = 0;
  reparam->sparams         = NULL;
  reparam->new_params      = NULL;
  reparam->sparams_name_id = g_hash_table_new_full (&g_str_hash, &g_str_equal, &g_free, NULL);
}

static void
_ncm_reparam_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_reparam_parent_class)->constructed (object);
  {
    NcmReparam *reparam = NCM_REPARAM (object);

    reparam->new_params = ncm_vector_new (reparam->length);
    reparam->sparams = g_ptr_array_sized_new (reparam->length);
    g_ptr_array_set_size (reparam->sparams, reparam->length);
  }
}

static void
_ncm_reparam_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmReparam *reparam = NCM_REPARAM (object);

  g_return_if_fail (NCM_IS_REPARAM (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, reparam->length);
      break;
    case PROP_PARAMS_DESC:
    {
      g_value_take_variant (value, ncm_reparam_get_params_desc_dict (reparam));
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_reparam_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmReparam *reparam = NCM_REPARAM (object);

  g_return_if_fail (NCM_IS_REPARAM (object));

  switch (prop_id)
  {
    case PROP_LEN:
      reparam->length = g_value_get_uint (value);
      break;
    case PROP_PARAMS_DESC:
      ncm_reparam_set_params_desc_dict (reparam, g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_reparam_finalize (GObject *object)
{
  NcmReparam *reparam = NCM_REPARAM (object);
  
  if (reparam->sparams_name_id != NULL)
  {
    g_hash_table_unref (reparam->sparams_name_id);
    reparam->sparams_name_id = NULL;
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_reparam_parent_class)->finalize (object);
}

static void
_ncm_reparam_copyto (NcmReparam *reparam, NcmReparam *reparam_dest)
{
  guint i;
  g_assert (G_OBJECT_TYPE (reparam) == G_OBJECT_TYPE (reparam_dest));
  if (reparam->length != reparam_dest->length)
  {
    ncm_vector_free (reparam_dest->new_params);
    reparam_dest->new_params = ncm_vector_new (reparam->length);
    g_ptr_array_set_size (reparam_dest->sparams, 0);
    g_ptr_array_set_size (reparam_dest->sparams, reparam->length);
  }

  ncm_vector_memcpy (reparam_dest->new_params, reparam->new_params);
  if (reparam->sparams)
  {
    for (i = 0; i < reparam->length; i++)
    {
      NcmSParam *pinfo = g_ptr_array_index (reparam->sparams, i);
      if (pinfo)
        g_ptr_array_index (reparam_dest->sparams, i) = ncm_sparam_copy (pinfo);
    }
  }
}

static void
ncm_reparam_class_init (NcmReparamClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  klass->copyto = &_ncm_reparam_copyto;
  klass->copy = NULL;
  klass->old2new = NULL;
  klass->new2old = NULL;
  klass->jac = NULL;

  object_class->set_property = &_ncm_reparam_set_property;
  object_class->get_property = &_ncm_reparam_get_property;
  object_class->constructed  = &_ncm_reparam_constructed;
  object_class->finalize     = &_ncm_reparam_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "System's length",
                                                      0.0, G_MAXUINT, 0.0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class, 
                                   PROP_PARAMS_DESC,
                                   g_param_spec_variant ("params-desc",
                                                         NULL, 
                                                         "News parameter descriptions", 
                                                         G_VARIANT_TYPE (NCM_REPARAM_PARAMS_DESC_DICT_TYPE), NULL, 
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  
}

/**
 * ncm_reparam_copyto:
 * @reparam: a #NcmReparam
 * @reparam_dest: a #NcmReparam
 *
 * FIXME
 */
void
ncm_reparam_copyto (NcmReparam *reparam, NcmReparam *reparam_dest)
{
  NCM_REPARAM_GET_CLASS (reparam)->copyto (reparam, reparam_dest);
}

/**
 * ncm_reparam_copy: (skip)
 * @reparam: a #NcmReparam
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmReparam *
ncm_reparam_copy (NcmReparam *reparam)
{
  if (!NCM_REPARAM_GET_CLASS (reparam)->copy)
  {
    g_error ("NcmReparam[%s] base class do not implement copy.", G_OBJECT_TYPE_NAME (reparam));
    return NULL;
  }
  else
    return NCM_REPARAM_GET_CLASS (reparam)->copy (reparam);
}

/**
 * ncm_reparam_ref:
 * @reparam: a #NcmReparam
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
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
 * FIXME
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
 * FIXME
 */
void
ncm_reparam_clear (NcmReparam **reparam)
{
  g_clear_object (reparam);
}

/**
 * ncm_reparam_old2new: (virtual old2new)
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 * @src: a #NcmVector
 * @dest: a #NcmVector
 *
 * FIXME
 */
void
ncm_reparam_old2new (NcmReparam *reparam, NcmModel *model, NcmVector *src, NcmVector *dest)
{
  NCM_REPARAM_GET_CLASS (reparam)->old2new (reparam, model, src, dest);
}

/**
 * ncm_reparam_new2old: (virtual new2old)
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 * @src: a #NcmVector
 * @dest: a #NcmVector
 *
 * FIXME
 */
void
ncm_reparam_new2old (NcmReparam *reparam, NcmModel *model, NcmVector *src, NcmVector *dest)
{
  NCM_REPARAM_GET_CLASS (reparam)->new2old (reparam, model, src, dest);
}

/**
 * ncm_reparam_jac: (virtual jac)
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 * @jac: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_reparam_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  NCM_REPARAM_GET_CLASS (reparam)->jac (reparam, model, jac);
}

/**
 * ncm_reparam_grad_old2new:
 * @reparam: a #NcmReparam
 * @model: FIXME
 * @jac: a #NcmMatrix
 * @old_grad: a #NcmVector
 * @new_grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_reparam_grad_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmVector *old_grad, NcmVector *new_grad)
{
  gint ret;
  NCM_REPARAM_GET_CLASS (reparam)->jac (reparam, model, jac);
  ret = gsl_blas_dgemv (CblasTrans, 1.0, ncm_matrix_gsl (jac), ncm_vector_gsl (old_grad), 0.0, ncm_vector_gsl (new_grad));
  NCM_TEST_GSL_RESULT("ncm_reparam_grad_old2new", ret);
}

/**
 * ncm_reparam_M_old2new:
 * @reparam: a #NcmReparam
 * @model: FIXME
 * @jac: a #NcmMatrix
 * @old_M: a #NcmMatrix
 * @new_M: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_reparam_M_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmMatrix *old_M, NcmMatrix *new_M)
{
  gint ret;
  NCM_REPARAM_GET_CLASS (reparam)->jac (reparam, model, jac);
  ret = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, ncm_matrix_gsl (old_M), ncm_matrix_gsl (jac), 0.0, ncm_matrix_gsl (new_M));
  NCM_TEST_GSL_RESULT("ncm_reparam_jac_old2new", ret);
}

/**
 * ncm_reparam_get_params_desc_dict: 
 * @reparam: a #NcmReparam.
 * 
 * Returns a #GVariant containing a dictionary describing the new parameters.
 * 
 * Returns: a #GVariant dictionary.
 */
GVariant *
ncm_reparam_get_params_desc_dict (NcmReparam *reparam)
{
  GVariantBuilder *builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_REPARAM_PARAMS_DESC_DICT_TYPE));
  
  guint i;
  for (i = 0; i < reparam->length; i++)
  {
    NcmSParam *sp_i = g_ptr_array_index (reparam->sparams, i); 
    if (sp_i != NULL)
    {
      GVariant *i_var = g_variant_new_uint32 (i);
      GVariant *sp_var_i = ncm_serialize_global_to_variant (G_OBJECT (sp_i));
      GVariant *entry_i = g_variant_new_dict_entry (i_var, sp_var_i);
      g_variant_builder_add_value (builder, entry_i);
    }
  }

  return g_variant_ref_sink (g_variant_builder_end (builder));
}

/**
 * ncm_reparam_set_params_desc_dict:
 * @reparam: a #NcmReparam.
 * @pdesc_dict: a #GVariant containing the new parameters descriptions.
 * 
 * Sets the new parameters descriptions using the information from @pdesc_dict.
 * 
 */
void 
ncm_reparam_set_params_desc_dict (NcmReparam *reparam, GVariant *pdesc_dict)
{
  if (pdesc_dict != NULL)
  {
    g_assert (g_variant_is_of_type (pdesc_dict, G_VARIANT_TYPE (NCM_REPARAM_PARAMS_DESC_DICT_TYPE)));
    {
      gsize n = g_variant_n_children (pdesc_dict);
      guint i;
      for (i = 0; i < n; i++)
      {
        GVariant *sp_var = NULL;
        guint j;
        NcmSParam *sp = NULL;
        g_variant_get_child (pdesc_dict, i, "{u@"NCM_SERIALIZE_OBJECT_TYPE"}", &j, &sp_var);
        sp = NCM_SPARAM (ncm_serialize_global_from_variant (sp_var));
        ncm_reparam_set_param_desc (reparam, j, sp);
        ncm_sparam_free (sp);
        g_variant_unref (sp_var);
      }
    }
  }
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
  NcmSParam **old_sp = (NcmSParam **) &g_ptr_array_index (reparam->sparams, i); 
  g_assert (i < reparam->sparams->len);
  
  if (*old_sp != NULL)
  {
    g_assert (g_hash_table_remove (reparam->sparams_name_id, ncm_sparam_name (*old_sp)));
    ncm_sparam_clear (old_sp);  
  }
  
  g_hash_table_insert (reparam->sparams_name_id, 
                       g_strdup (ncm_sparam_name (sp)), 
                       GUINT_TO_POINTER (i));

  g_ptr_array_index (reparam->sparams, i) = ncm_sparam_ref (sp);
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
  g_assert (i < reparam->sparams->len);
  return g_ptr_array_index (reparam->sparams, i);
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
  g_assert (i < reparam->sparams->len);
  if (g_ptr_array_index (reparam->sparams, i) != NULL)
    return ncm_sparam_ref (g_ptr_array_index (reparam->sparams, i));
  else
    return NULL;
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
 * FIXME
 * 
 */
void 
ncm_reparam_set_param_desc_full (NcmReparam *reparam, guint i, gchar *name, gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype)
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
ncm_reparam_index_from_name (NcmReparam *reparam, gchar *param_name, guint *i)
{
  gpointer param_id;
  gboolean found = g_hash_table_lookup_extended (reparam->sparams_name_id, param_name, NULL, &param_id);
  if (found)
    *i = GPOINTER_TO_UINT (param_id);
  else
    *i = -1; /* Yes, I known. */
  return found;
}

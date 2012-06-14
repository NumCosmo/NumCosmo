/***************************************************************************
 *            ncm_model.c
 *
 *  Fri February 24 21:18:21 2012
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
 * SECTION:ncm_model
 * @title: Generic Model Pure Virtual Class
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

enum
{
  PROP_0,
  PROP_NAME,
  PROP_NICK,
  PROP_SPARAMS_LEN,
  PROP_VPARAMS_LEN,
  PROP_PARAMS,
  PROP_IMPLEMENTATION,
  PROP_PTYPES,
  PROP_REPARAM,
};

G_DEFINE_ABSTRACT_TYPE (NcmModel, ncm_model, G_TYPE_OBJECT);

/**
 * ncm_model_copy: (skip)
 * @model: a #NcmModel.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmModel *
ncm_model_copy (NcmModel *model)
{
  return NCM_MODEL_GET_CLASS (model)->copy (model);
}

/**
 * ncm_model_copyto:
 * @model: a #NcmModel.
 * @model_dest: a #NcmModel.
 *
 * FIXME
 *
 */
void
ncm_model_copyto (NcmModel *model, NcmModel *model_dest)
{
  NCM_MODEL_GET_CLASS (model)->copyto (model, model_dest);
}

/**
 * ncm_model_free:
 * @model: a #NcmModel.
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

static void
ncm_model_init (NcmModel *model)
{
  NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);

  model->sparams = g_ptr_array_new_with_free_func ((GDestroyNotify)&ncm_sparam_free);
  model->sparams_name_id = g_hash_table_new_full (&g_str_hash, &g_str_equal, &g_free, NULL);
  model->params = NULL;
  model->p = NULL;

  model->vparam_len = g_array_sized_new (TRUE, TRUE, sizeof (guint), model_class->vparam_len);
  g_array_set_size (model->vparam_len, model_class->vparam_len);
  model->vparam_pos = g_array_sized_new (TRUE, TRUE, sizeof (guint), model_class->vparam_len);
  g_array_set_size (model->vparam_pos, model_class->vparam_len);

  model->pkey = 1;
  model->reparam = NULL;
  model->ptypes = g_array_new (FALSE, TRUE, sizeof (NcmParamType));
}

static void
_ncm_model_reset_sparams (NcmModel *model)
{
  NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);
  guint i;

  g_hash_table_remove_all (model->sparams_name_id);
  g_ptr_array_set_size (model->sparams, 0);
  g_ptr_array_set_size (model->sparams, model->total_len);

  for (i = 0; i < model_class->sparam_len; i++)
  {
	NcmSParam *sp = g_ptr_array_index (model_class->sparam, i);
	g_array_index (model->ptypes, NcmParamType, i) = ncm_sparam_get_fit_type (sp);
	g_ptr_array_index (model->sparams, i) = ncm_sparam_copy (sp);
	g_hash_table_insert (model->sparams_name_id, g_strdup (ncm_sparam_name (sp)), GUINT_TO_POINTER (i));
  }

  for (i = 0; i < model_class->vparam_len; i++)
  {
	const guint len = g_array_index (model->vparam_len, guint, i);
	const guint pos = g_array_index (model->vparam_pos, guint, i);
	NcmVParam *vp = ncm_vparam_copy (g_ptr_array_index (model_class->vparam, i));
	gint j;

	ncm_vparam_set_len (vp, len);

	for (j = 0; j < len; j++)
	{
	  const guint n = pos + j;
	  NcmSParam *sp = ncm_vparam_peek_sparam (vp, j);
	  g_array_index (model->ptypes, NcmParamType, n) = ncm_sparam_get_fit_type (sp);
	  g_ptr_array_index (model->sparams, n) = ncm_sparam_ref (sp);
	  g_hash_table_insert (model->sparams_name_id, g_strdup (ncm_sparam_name (sp)), GUINT_TO_POINTER (n));
	}

	ncm_vparam_free (vp);
  }
}

static void
_ncm_model_reset_sparams_from_reparam (NcmModel *model)
{
  guint i;

  g_hash_table_remove_all (model->sparams_name_id);
  g_ptr_array_set_size (model->sparams, 0);
  g_ptr_array_set_size (model->sparams, model->total_len);

  for (i = 0; i < model->total_len; i++)
  {
	NcmSParam *sp = g_ptr_array_index (model->reparam->sparams, i);
	if (sp == NULL)
	  sp = g_ptr_array_index (model->sparams, i);
	else
	{
	  ncm_sparam_free (g_ptr_array_index (model->sparams, i));
	  g_ptr_array_index (model->sparams, i) = ncm_sparam_copy (sp);
	}
	g_array_index (model->ptypes, NcmParamType, i) = ncm_sparam_get_fit_type (sp);
	g_hash_table_insert (model->sparams_name_id, g_strdup (ncm_sparam_name (sp)), GUINT_TO_POINTER (i));
  }
}

static void
_ncm_model_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_model_parent_class)->constructed (object);
  {
	NcmModel *model = NCM_MODEL (object);
	NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);
	gint i;

	model->total_len = model_class->sparam_len;
	for (i = 0; i < model_class->vparam_len; i++)
	{
	  g_array_index (model->vparam_pos, guint, i) = model->total_len;
	  model->total_len += g_array_index (model->vparam_len, guint, i);
	}
	model->params = ncm_vector_ref (ncm_vector_new (model->total_len));
	model->p = model->params;
	g_array_set_size (model->ptypes, model->total_len);
	_ncm_model_reset_sparams (model);
	ncm_model_params_set_default (model);
  }
}

static void
_ncm_model_dispose (GObject *object)
{
  NcmModel *model = NCM_MODEL (object);
  ncm_vector_free (model->params);
  model->params = NULL;
  model->p = NULL;

  if (model->reparam != NULL)
	ncm_reparam_free (model->reparam);

  g_array_unref (model->vparam_len);
  g_array_unref (model->vparam_pos);
  g_array_unref (model->ptypes);

  g_ptr_array_unref (model->sparams);

  g_hash_table_unref (model->sparams_name_id);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_parent_class)->dispose (object);
}

static void
_ncm_model_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_parent_class)->finalize (object);
}

static void
_ncm_model_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModel *model = NCM_MODEL (object);
  g_return_if_fail (NCM_IS_GMODEL (object));

  switch (prop_id)
  {
	case PROP_REPARAM:
	  ncm_model_set_reparam (model, g_value_get_object (value));
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
_ncm_model_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModel *model = NCM_MODEL (object);
  NcmModelClass *model_class = NCM_MODEL_GET_CLASS (object);

  g_return_if_fail (NCM_IS_GMODEL (object));

  switch (prop_id)
  {
	case PROP_NAME:
	  g_value_set_string (value, model_class->name);
	  break;
	case PROP_NICK:
	  g_value_set_string (value, model_class->nick);
	  break;
	case PROP_PARAMS:
	  g_value_set_object (value, ncm_vector_copy (model->params));
	  break;
	case PROP_SPARAMS_LEN:
	  g_value_set_uint (value, model_class->sparam_len);
	  break;
	case PROP_VPARAMS_LEN:
	  g_value_set_uint (value, model_class->vparam_len);
	  break;
	case PROP_IMPLEMENTATION:
	  g_value_set_ulong (value, model_class->impl);
	  break;
	case PROP_REPARAM:
	  g_value_set_object (value, model->reparam);
	  break;
	case PROP_PTYPES:
	  g_value_set_boxed (value, model->ptypes);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static NcmModel *
_ncm_model_copy (NcmModel *model)
{
  NcmModel *model_copy = g_object_new (G_OBJECT_TYPE (model), NULL);
  ncm_model_copyto (model, model_copy);
  return model_copy;
}

static void
_ncm_model_copyto (NcmModel *model, NcmModel *model_dest)
{
  guint i;
  g_assert (ncm_model_is_equal (model, model_dest));
  ncm_vector_memcpy (model_dest->p, model->p);
  for (i = 0; i < model->ptypes->len; i++)
	g_array_index (model_dest->ptypes, NcmParamType, i) = g_array_index (model->ptypes, NcmParamType, i);
  ncm_model_params_update (model_dest);
}

static void
ncm_model_class_init (NcmModelClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  object_class->constructed  = &_ncm_model_constructed;
  object_class->dispose      = &_ncm_model_dispose;
  object_class->finalize     = &_ncm_model_finalize;
  object_class->set_property = &_ncm_model_set_property;
  object_class->get_property = &_ncm_model_get_property;
  klass->copyto = &_ncm_model_copyto;
  klass->copy   = &_ncm_model_copy;
  klass->set_property = NULL;
  klass->get_property = NULL;

  klass->model_id          = -1;
  klass->name               = NULL;
  klass->nick               = NULL;
  klass->nonparam_prop_len  = 0;
  klass->sparam_len         = 0;
  klass->vparam_len         = 0;
  klass->sparam             = NULL;
  klass->vparam             = NULL;

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
                                   PROP_PARAMS,
                                   g_param_spec_object ("params",
                                                        NULL,
                                                        "Parameters vector",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_IMPLEMENTATION,
                                   g_param_spec_ulong  ("implementation",
                                                        NULL,
                                                        "Bitwise specification of functions implementation",
                                                        0, G_MAXULONG, 0,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PTYPES,
                                   g_param_spec_boxed  ("params-types",
                                                         NULL,
                                                         "Parameters' types",
                                                         G_TYPE_ARRAY,
                                                         G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_REPARAM,
                                   g_param_spec_object  ("reparam",
                                                         NULL,
                                                         "Model reparametrization",
                                                         NCM_TYPE_REPARAM,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

G_LOCK_DEFINE_STATIC (last_model_id);

void
ncm_model_register_id (NcmModelClass *model_class)
{
  if (model_class->model_id < 0)
  {
	static NcmModelID last_model_id = 0;
	G_LOCK (last_model_id);
	model_class->model_id = last_model_id++;
	G_UNLOCK (last_model_id);
  }
  else
  {
	g_error ("This model or its parent is already registred, id = %d. This function must be use once and only in the defining model.", model_class->model_id);
  }
  if (model_class->model_id > NCM_MODEL_MAX_ID)
	g_error ("Max model id was already attained. Increase by altering NCM_MODEL_MAX_ID.");

  return;
}

/**
 * ncm_model_class_get_property: (skip)
 * @object: a #NcmModelClass.
 * @prop_id: FIXME
 * @value: FIXME
 * @pspec: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_class_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModel *model = NCM_MODEL (object);
  NcmModelClass *model_class = NCM_MODEL_CLASS (g_type_class_peek_static (pspec->owner_type));
  const guint sparam_id = prop_id - model_class->nonparam_prop_len + model_class->parent_sparam_len;
  const guint vparam_id = sparam_id	- model_class->sparam_len + model_class->parent_vparam_len;
  const guint vparam_len_id = vparam_id	- model_class->vparam_len + model_class->parent_vparam_len;

  if (prop_id < model_class->nonparam_prop_len && model_class->get_property)
  {
	model_class->get_property (object, prop_id, value, pspec);
  }
  else if (sparam_id < model_class->sparam_len)
  {
	g_value_set_double (value, ncm_vector_get (model->params, sparam_id));
  }
  else if (vparam_id < model_class->vparam_len)
  {
	const guint len = g_array_index (model->vparam_len, guint, vparam_id);
	const guint pos = g_array_index (model->vparam_pos, guint, vparam_id);
	NcmVector *v = ncm_vector_get_subvector (model->params, pos, len);
	g_value_take_object (value, ncm_vector_ref (v));
  }
  else if (vparam_len_id < model_class->vparam_len)
  {
	g_value_set_uint (value, g_array_index (model->vparam_len, guint, vparam_len_id));
  }
  else
	g_assert_not_reached ();
}

/**
 * ncm_model_class_set_property: (skip)
 * @object: a #NcmModelClass.
 * @prop_id: FIXME
 * @value: FIXME
 * @pspec: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_class_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModel *model = NCM_MODEL (object);
  NcmModelClass *model_class = NCM_MODEL_CLASS (g_type_class_peek_static (pspec->owner_type));
  const guint sparam_id = prop_id - model_class->nonparam_prop_len + model_class->parent_sparam_len;
  const guint vparam_id = sparam_id	- model_class->sparam_len + model_class->parent_vparam_len;
  const guint vparam_len_id = vparam_id	- model_class->vparam_len + model_class->parent_vparam_len;

  if (prop_id < model_class->nonparam_prop_len && model_class->set_property)
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
	NcmVector *val = g_value_get_object (value);
	ncm_model_orig_vparam_set_vector (model, vparam_id, val);
  }
  else if (vparam_len_id < model_class->vparam_len)
  {
	guint psize = g_value_get_uint (value);
	g_array_index (model->vparam_len, guint, vparam_len_id) = psize;
  }
  else
	g_assert_not_reached ();
}

/**
 * ncm_model_class_add_params: (skip)
 * @model_class: a #NcmModelClass.
 * @sparam_len: FIXME
 * @vparam_len: FIXME
 * @nonparam_prop_len: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_class_add_params (NcmModelClass *model_class, guint sparam_len, guint vparam_len, guint nonparam_prop_len)
{
  model_class->parent_sparam_len = model_class->sparam_len;
  model_class->parent_vparam_len = model_class->vparam_len;
  model_class->sparam_len += sparam_len;
  model_class->vparam_len += vparam_len;
  model_class->nonparam_prop_len = nonparam_prop_len;

  if (model_class->sparam_len > 0)
  {
	if (model_class->sparam == NULL)
	{
	  model_class->sparam = g_ptr_array_new_with_free_func ((GDestroyNotify) &ncm_sparam_free);
	  g_ptr_array_set_size (model_class->sparam, model_class->sparam_len);
	}
	else
	{
	  GPtrArray *sparam = g_ptr_array_new_with_free_func ((GDestroyNotify) &ncm_sparam_free);
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
	  model_class->vparam = g_ptr_array_new_with_free_func ((GDestroyNotify) &ncm_vparam_free);
	  g_ptr_array_set_size (model_class->vparam, model_class->vparam_len);
	}
	else
	{
	  GPtrArray *vparam = g_ptr_array_new_with_free_func ((GDestroyNotify) &ncm_vparam_free);
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
 * @model_class: a #NcmModelClass.
 * @name: FIXME
 * @nick: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_class_set_name_nick (NcmModelClass *model_class, const gchar *name, const gchar *nick)
{
  g_free (model_class->name);
  g_free (model_class->nick);

  model_class->name = g_strdup (name);
  model_class->nick = g_strdup (nick);
}

/**
 * ncm_model_class_set_sparam:
 * @model_class: a #NcmModelClass.
 * @sparam_id: FIXME
 * @symbol: FIXME
 * @name: FIXME
 * @lower_bound: FIXME
 * @upper_bound: FIXME
 * @scale: FIXME
 * @abstol: FIXME
 * @default_value: FIXME
 * @ppt: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_class_set_sparam (NcmModelClass *model_class, guint sparam_id, gchar *symbol, gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt)
{
  GObjectClass* object_class = G_OBJECT_CLASS (model_class);
  const guint prop_id = sparam_id + model_class->nonparam_prop_len - model_class->parent_sparam_len;
  NcmSParam *sparam = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_value, ppt);

  g_assert (prop_id > 0);

  if (g_ptr_array_index (model_class->sparam, sparam_id) != NULL)
	g_error ("Scalar Parameter: %u is already set\n", sparam_id);

  g_ptr_array_index (model_class->sparam, sparam_id) = sparam;
  g_object_class_install_property (object_class, prop_id,
                                   g_param_spec_double (name, NULL, symbol,
                                                        lower_bound, upper_bound, default_value,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

/**
 * ncm_model_class_set_vparam:
 * @model_class: a #NcmModelClass
 * @vparam_id: FIXME
 * @default_length: FIXME
 * @symbol: FIXME
 * @name: FIXME
 * @lower_bound: FIXME
 * @upper_bound: FIXME
 * @scale: FIXME
 * @abstol: FIXME
 * @default_value: FIXME
 * @ppt: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_class_set_vparam (NcmModelClass *model_class, guint vparam_id, guint default_length, gchar *symbol, gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt)
{
  GObjectClass* object_class = G_OBJECT_CLASS (model_class);
  const guint prop_id = vparam_id + model_class->nonparam_prop_len - model_class->parent_vparam_len + (model_class->sparam_len - model_class->parent_sparam_len);
  const guint prop_len_id = prop_id + (model_class->vparam_len - model_class->parent_vparam_len);
  NcmVParam *vparam = ncm_vparam_full_new (default_length, name, symbol, lower_bound, upper_bound, scale, abstol, default_value, ppt);

  g_assert (prop_id > 0);
  g_assert (prop_len_id > 0);

  if (g_ptr_array_index (model_class->vparam, vparam_id) != NULL)
	g_error ("Vector Parameter: %u is already set\n", vparam_id);

  g_ptr_array_index (model_class->vparam, vparam_id) = vparam;
  g_object_class_install_property (object_class, prop_id,
                                   g_param_spec_object (name, NULL, symbol, NCM_TYPE_VECTOR, G_PARAM_READWRITE));

  g_object_class_install_property (object_class, prop_len_id,
                                   g_param_spec_uint (g_strdup_printf ("%s-length", name),
                                                      NULL,
                                                      g_strdup_printf ("Lenght of %s vector parameter", symbol),
                                                      0, G_MAXUINT, default_length,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY));
}

/**
 * ncm_model_class_check_params_info:
 * @model_class: a #NcmModelClass.
 *
 * FIXME
 *
 */
void
ncm_model_class_check_params_info (NcmModelClass *model_class)
{
  gulong i;
  guint total_params_len = model_class->sparam_len + model_class->vparam_len;
  if (!total_params_len)
  {
	g_error ("Class size or params not initialized, call ncm_model_class_add_params.");
  }

  for (i = 0; i < model_class->sparam_len; i++)
  {
	if (g_ptr_array_index (model_class->sparam, i) == NULL)
	{
	  g_error ("Class (%s) didn't initialized scalar parameter %lu/%u", model_class->name ? model_class->name : "no-name", i + 1, model_class->sparam_len);
	}
	/* g_debug ("Model[%s][%s] id %lu\n", model_class->name, ((NcmSParam *)g_ptr_array_index (model_class->params_info, i))->name, i); */
  }

  for (i = 0; i < model_class->vparam_len; i++)
  {
	if (g_ptr_array_index (model_class->vparam, i) == NULL)
	{
	  g_error ("Class (%s) didn't initialized vector parameter %lu/%u", model_class->name ? model_class->name : "no-name", i + 1, model_class->vparam_len);
	}
	/* g_debug ("Model[%s][%s] id %lu\n", model_class->name, ((NcmSParam *)g_ptr_array_index (model_class->params_info, i))->name, i); */
  }
}

static void _ncm_model_reset_sparams (NcmModel *model);
static void _ncm_model_reset_sparams_from_reparam (NcmModel *model);

/**
 * ncm_model_set_reparam:
 * @model: a #NcmModel.
 * @reparam: a #NcmReparam.
 *
 * FIXME
 *
 */
void
ncm_model_set_reparam (NcmModel *model, NcmReparam *reparam)
{
  if (reparam != NULL)
  {
	model->reparam = reparam;
	model->p = model->reparam->new_params;
	ncm_reparam_old2new (model->reparam, model, model->params, model->reparam->new_params);
	_ncm_model_reset_sparams_from_reparam (model);
  }
  else
  {
	model->reparam = NULL;
	model->p = model->params;
	_ncm_model_reset_sparams (model);
  }
}

/**
 * ncm_model_is_equal:
 * @model1: a #NcmModel.
 * @model2: a #NcmModel.
 *
 * Compares if model1 and model2 are the same,
 * with same dimension and reparametrization.
 *
 */
gboolean
ncm_model_is_equal (NcmModel *model1, NcmModel *model2)
{
  if (G_OBJECT_TYPE (model1) != G_OBJECT_TYPE (model2))
	return FALSE;
  if (ncm_vector_len (model1->params) != ncm_vector_len (model2->params))
	return FALSE;
  if (model1->reparam)
  {
	if (model2->reparam == NULL)
	  return FALSE;
	if (G_OBJECT_TYPE (model1->reparam) != G_OBJECT_TYPE (model2->reparam))
	  return FALSE;
  }
  return TRUE;
}

/**
 * ncm_model_get_reparam:
 * @model: a #NcmModel.
 *
 * FIXME
 *
 * Returns: FIXME
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
 * ncm_model_reparam_df:
 * @model: a #NcmModel.
 * @fv: a #NcmVector.
 * @v: a #NcmVector.
 *
 * FIXME
 *
 */
void
ncm_model_reparam_df (NcmModel *model, NcmVector *fv, NcmVector *v)
{
  g_assert (model->reparam);
  g_assert_not_reached ();
  ncm_reparam_grad_old2new (model->reparam, model, NULL, fv, v);
}

/**
 * ncm_model_reparam_J:
 * @model: a #NcmModel.
 * @fJ: a #NcmMatrix.
 * @J: a #NcmMatrix.
 *
 * FIXME
 *
 */
void
ncm_model_reparam_J (NcmModel *model, NcmMatrix *fJ, NcmMatrix *J)
{
  g_assert (model->reparam);
  g_assert_not_reached ();
  ncm_reparam_M_old2new (model->reparam, model, NULL, fJ, J);
}

/**
 * ncm_model_params_set_default:
 * @model: a #NcmModel.
 *
 * FIXME
 *
 */
void
ncm_model_params_set_default (NcmModel *model)
{
  gint i;

  for (i = 0; i < model->total_len; i++)
  {
	const NcmSParam *p = g_ptr_array_index (model->sparams, i);
	ncm_vector_set (model->params, i, ncm_sparam_get_default_value (p));
  }

  ncm_model_orig_params_update (model);
}

/**
 * ncm_model_params_save_as_default:
 * @model: a #NcmModel.
 *
 * FIXME
 *
 */
void
ncm_model_params_save_as_default (NcmModel *model)
{
  gint i;
  for (i = 0; i < model->total_len; i++)
  {
	NcmSParam *p = g_ptr_array_index (model->sparams, i);
	ncm_sparam_set_default_value (p, ncm_vector_get (model->params, i));
  }
}

/**
 * ncm_model_params_copyto:
 * @model: a #NcmModel.
 * @model_dest: a #NcmModel.
 *
 * FIXME
 *
 */
void
ncm_model_params_copyto (NcmModel *model, NcmModel *model_dest)
{
  g_assert (ncm_model_is_equal (model, model_dest));
  ncm_model_params_set_vector (model_dest, model->p);
}

/**
 * ncm_model_params_set_all:
 * @model: a #NcmModel.
 * @...: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_params_set_all (NcmModel *model, ...)
{
  guint i;
  va_list ap;
  va_start(ap, model);

  for (i = 0; i < ncm_vector_len (model->p); i++)
	ncm_vector_set (model->p, i, va_arg (ap, gdouble));

  va_end(ap);

  ncm_model_params_update (model);

  return;
}

/**
 * ncm_model_params_set_all_data:
 * @model: a #NcmModel.
 * @data: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_params_set_all_data (NcmModel *model, gdouble *data)
{
  guint i;

  for (i = 0; i < ncm_vector_len (model->p); i++)
	ncm_vector_set (model->p, i, data[i]);

  ncm_model_params_update (model);
  return;
}

/**
 * ncm_model_params_set_vector:
 * @model: a #NcmModel.
 * @v: a #NcmVector.
 *
 * FIXME
 *
 */
void
ncm_model_params_set_vector (NcmModel *model, NcmVector *v)
{
  ncm_vector_memcpy (model->p, v);
  ncm_model_params_update (model);
}

/**
 * ncm_model_params_set_model:
 * @model: a #NcmModel.
 * @model_src: a #NcmModel.
 *
 * FIXME
 *
 */
void
ncm_model_params_set_model (NcmModel *model, NcmModel *model_src)
{
  g_assert (ncm_model_is_equal (model, model_src));
  ncm_model_params_set_vector (model, model_src->p);
}

/**
 * ncm_model_params_print_all: (skip)
 * @model: a #NcmModel
 * @out: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_params_print_all (NcmModel *model, FILE *out)
{
  guint i;
  for (i = 0; i < ncm_vector_len (model->p); i++)
	fprintf (out, "  % 20.16g", ncm_vector_get (model->p, i));
  fprintf (out, "\n");
  fflush (out);
  return;
}

/**
 * ncm_model_params_log_all:
 * @model: a #NcmModel
 *
 * FIXME
 *
 */
void
ncm_model_params_log_all (NcmModel *model)
{
  guint i;
  for (i = 0; i < ncm_vector_len (model->p); i++)
	g_message ("  % 20.16g", ncm_vector_get (model->p, i));
  g_message ("\n");
  return;
}

/**
 * ncm_model_params_get_all:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
   */
NcmVector *
ncm_model_params_get_all (NcmModel *model)
{
  return ncm_vector_copy (model->p);
}

/**
 * ncm_model_ref:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
/**
 * ncm_model_id:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_impl:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_len:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_sparam_len:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_vparam_array_len:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_name:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
/**
 * ncm_model_nick:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
/**
 * ncm_model_peek_reparam:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
/**
 * ncm_model_param_finite:
 * @model: a #NcmModel
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_params_finite:
 * @model: a #NcmModel
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_params_update:
 * @model: a #NcmModel
 *
 * Force the parameters to the update its internal flags and
 * update the original parameters if necessary.
 *
 */
/**
 * ncm_model_orig_params_update:
 * @model: a #NcmModel
 *
 * Update the new parameters. It causes an error to call this
 * function with a model without reparametrization.
 *
 */
/**
 * ncm_model_vparam_index:
 * @model: a #NcmModel
 * @n: vector index
 * @i: vector component index
 *
 * FIXME
 *
 * Returns: index of the i-th component of the n-th vector
 */
/**
 * ncm_model_vparam_len:
 * @model: a #NcmModel
 * @n: vector index
 *
 * FIXME
 *
 * Returns: length of the n-th vector
 */
/**
 * ncm_model_param_set:
 * @model: a #NcmModel
 * @n: FIXME
 * @val: FIXME
 *
 * FIXME
 *
 */
/**
 * ncm_model_param_set_default:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 */
/**
 * ncm_model_param_get:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_orig_param_set:
 * @model: a #NcmModel
 * @n: FIXME
 * @val: FIXME
 *
 * FIXME
 *
 */
/**
 * ncm_model_orig_param_get:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_vparam_set:
 * @model: a #NcmModel
 * @n: FIXME
 * @i: FIXME
 * @val: FIXME
 *
 * FIXME
 *
 */
/**
 * ncm_model_vparam_get:
 * @model: a #NcmModel
 * @n: FIXME
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_orig_vparam_set_vector:
 * @model: a #NcmModel
 * @n: FIXME
 * @val: FIXME
 *
 * FIXME
 *
 */
/**
 * ncm_model_orig_vparam_get_vector:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
#ifndef NUMCOSMO_HAVE_INLINE
#define NUMCOSMO_HAVE_INLINE
#undef _NCM_MODEL_INLINE_H_
#define G_INLINE_FUNC
#include "ncm_model.h"
#undef _NCM_MODEL_INLINE_H_
#undef NUMCOSMO_HAVE_INLINE
#endif /* NUMCOSMO_HAVE_INLINE */

/**
 * ncm_model_param_get_scale:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_model_param_get_scale (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return ncm_sparam_get_scale (g_ptr_array_index (model->sparams, n));
}

/**
 * ncm_model_param_get_lower_bound:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_model_param_get_lower_bound (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return ncm_sparam_get_lower_bound (g_ptr_array_index (model->sparams, n));
}

/**
 * ncm_model_param_get_upper_bound:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_model_param_get_upper_bound (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return ncm_sparam_get_upper_bound (g_ptr_array_index (model->sparams, n));
}

/**
 * ncm_model_param_get_abstol:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_model_param_get_abstol (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return ncm_sparam_get_absolute_tolerance (g_ptr_array_index (model->sparams, n));
}

/**
 * ncm_model_param_get_ftype:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmParamType
ncm_model_param_get_ftype (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return g_array_index (model->ptypes, NcmParamType, n);
}

/**
 * ncm_model_param_set_scale:
 * @model: a #NcmModel
 * @n: FIXME
 * @scale: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_param_set_scale (NcmModel *model, guint n, const gdouble scale)
{
  g_assert (n < model->total_len);
  ncm_sparam_set_scale (g_ptr_array_index (model->sparams, n), scale);
}

/**
 * ncm_model_param_set_lower_bound:
 * @model: a #NcmModel
 * @n: FIXME
 * @lb: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_param_set_lower_bound (NcmModel *model, guint n, const gdouble lb)
{
  g_assert (n < model->total_len);
  ncm_sparam_set_lower_bound (g_ptr_array_index (model->sparams, n), lb);
}

/**
 * ncm_model_param_set_upper_bound:
 * @model: a #NcmModel
 * @n: FIXME
 * @ub: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_param_set_upper_bound (NcmModel *model, guint n, const gdouble ub)
{
  g_assert (n < model->total_len);
  ncm_sparam_set_upper_bound (g_ptr_array_index (model->sparams, n), ub);
}

/**
 * ncm_model_param_set_abstol:
 * @model: a #NcmModel
 * @n: FIXME
 * @abstol: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_param_set_abstol (NcmModel *model, guint n, const gdouble abstol)
{
  g_assert (n < model->total_len);
  ncm_sparam_set_absolute_tolerance (g_ptr_array_index (model->sparams, n), abstol);
}

/**
 * ncm_model_param_set_ftype:
 * @model: a #NcmModel
 * @n: FIXME
 * @ptype: FIXME
 *
 * FIXME
 *
 */
void
ncm_model_param_set_ftype (NcmModel *model, guint n, const NcmParamType ptype)
{
  g_assert (n < model->total_len);
  g_array_index (model->ptypes, NcmParamType, n) = ptype;
}

/**
 * ncm_model_param_get_name:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const gchar *
ncm_model_param_name (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return ncm_sparam_name (g_ptr_array_index (model->sparams, n));
}

/**
 * ncm_model_param_get_symbol:
 * @model: a #NcmModel
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
const gchar *
ncm_model_param_symbol (NcmModel *model, guint n)
{
  g_assert (n < model->total_len);
  return ncm_sparam_symbol (g_ptr_array_index (model->sparams, n));
}

/**
 * ncm_model_param_index_from_name:
 * @model: FIXME
 * @param_name: FIXME
 *
 * Returns: FIXME
 */
guint
ncm_model_param_index_from_name (NcmModel *model, gchar *param_name)
{
  gpointer param_id;
  gboolean found = g_hash_table_lookup_extended (model->sparams_name_id, param_name, NULL, &param_id);
  g_assert (found);
  return GPOINTER_TO_UINT (param_id);
}

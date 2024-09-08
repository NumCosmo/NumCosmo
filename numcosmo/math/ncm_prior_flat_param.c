/***************************************************************************
 *            ncm_prior_flat_param.c
 *
 *  Wed August 03 10:55:21 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_flat_param.c
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
 * SECTION:ncm_prior_flat_param
 * @title: NcmPriorFlatParam
 * @short_description: a flat prior on a sampling parameter.
 *
 * This object is a subclass of #NcmPriorFlat, tailored for defining a flat prior on a
 * sampling parameter. The prior is uniquely identified by the sampling parameter ID
 * and is characterized by user-specified lower and upper limits, along with the scale
 * of the prior.
 *
 * Users enjoy the flexibility to specify the parameter in various ways:
 * - Using the pair NcmModelID and the parameter pid.
 * - Providing a single NcmMSetPIndex.
 * - Supplying a string consisting of a parameter full name "model:parameter".
 *
 * These options offer diverse and convenient ways for parameter identification.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_flat_param.h"

enum
{
  PROP_0,
  PROP_MODEL_NS,
  PROP_STACK_POS,
  PROP_PARAM_NAME,
};

struct _NcmPriorFlatParam
{
  /*< private >*/
  NcmPriorFlat parent_instance;
  gchar *model_ns;
  gchar *param_name;
  guint stack_pos;
  NcmModelID mid;
};

G_DEFINE_TYPE (NcmPriorFlatParam, ncm_prior_flat_param, NCM_TYPE_PRIOR_FLAT)

static void
ncm_prior_flat_param_init (NcmPriorFlatParam *pfp)
{
}

static void
_ncm_prior_flat_param_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorFlatParam *pfp = NCM_PRIOR_FLAT_PARAM (object);

  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (object));

  switch (prop_id)
  {
    case PROP_MODEL_NS:
      ncm_prior_flat_param_set_model_ns (pfp, g_value_get_string (value), NULL);
      break;
    case PROP_STACK_POS:
      ncm_prior_flat_param_set_stack_pos (pfp, g_value_get_uint (value));
      break;
    case PROP_PARAM_NAME:
      ncm_prior_flat_param_set_param_name (pfp, g_value_get_string (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_flat_param_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorFlatParam *pfp = NCM_PRIOR_FLAT_PARAM (object);

  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (object));

  switch (prop_id)
  {
    case PROP_MODEL_NS:
      g_value_set_string (value, ncm_prior_flat_param_peek_model_ns (pfp));
      break;
    case PROP_STACK_POS:
      g_value_set_uint (value, ncm_prior_flat_param_get_stack_pos (pfp));
      break;
    case PROP_PARAM_NAME:
      g_value_set_string (value, ncm_prior_flat_param_peek_param_name (pfp));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_flat_param_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_flat_param_parent_class)->finalize (object);
}

static gdouble _ncm_prior_flat_param_mean (NcmPriorFlat *pf, NcmMSet *mset);

static void
ncm_prior_flat_param_class_init (NcmPriorFlatParamClass *klass)
{
  GObjectClass *object_class  = G_OBJECT_CLASS (klass);
  NcmPriorFlatClass *pf_class = NCM_PRIOR_FLAT_CLASS (klass);

  object_class->set_property = &_ncm_prior_flat_param_set_property;
  object_class->get_property = &_ncm_prior_flat_param_get_property;
  object_class->finalize     = &_ncm_prior_flat_param_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MODEL_NS,
                                   g_param_spec_string ("model-ns",
                                                        NULL,
                                                        "model namespace",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_STACK_POS,
                                   g_param_spec_uint ("stack-pos",
                                                      NULL,
                                                      "stack position",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PARAM_NAME,
                                   g_param_spec_string ("parameter-name",
                                                        NULL,
                                                        "parameter name",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  pf_class->mean = &_ncm_prior_flat_param_mean;
}

static gdouble
_ncm_prior_flat_param_mean (NcmPriorFlat *pf, NcmMSet *mset)
{
  NcmPriorFlatParam *pfp = NCM_PRIOR_FLAT_PARAM (pf);
  NcmModel *model        = ncm_mset_peek (mset, pfp->mid);

  return ncm_model_param_get_by_name (model, pfp->param_name, NULL);
}

/**
 * ncm_prior_flat_param_new:
 * @model: a #NcmModel
 * @pid: parameter id
 * @x_low: parameter lower limit
 * @x_upp: parameter upper limit
 * @scale: parameter scale
 *
 * Creates a new Flat prior for parameter @pid of @model.
 *
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_new (NcmModel *model, guint pid, gdouble x_low, gdouble x_upp, gdouble scale)
{
  const gchar *model_ns   = G_OBJECT_TYPE_NAME (model);
  const gchar *param_name = ncm_model_param_name (model, pid);

  NcmPriorFlatParam *pfp = g_object_new (NCM_TYPE_PRIOR_FLAT_PARAM,
                                         "x-low", x_low,
                                         "x-upp", x_upp,
                                         "scale", scale,
                                         "model-ns", model_ns,
                                         "parameter-name", param_name,
                                         NULL);

  return pfp;
}

/**
 * ncm_prior_flat_param_new_name:
 * @name: parameter name
 * @x_low: parameter lower limit
 * @x_upp: parameter upper limit
 * @scale: parameter scale
 * @error: a #GError
 *
 * Creates a new Flat prior for parameter named @name in @mset. See
 * ncm_mset_split_full_name() for details on the parameter name format.
 *
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_new_name (const gchar *name, gdouble x_low, gdouble x_upp, gdouble scale, GError **error)
{
  gchar *model_ns          = NULL;
  gchar *param_name        = NULL;
  guint stack_pos          = 0;
  gboolean full_name_found = FALSE;

  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  full_name_found = ncm_mset_split_full_name (name, &model_ns, &stack_pos, &param_name, error);
  NCM_UTIL_ON_ERROR_RETURN (error, , NULL);

  if (!full_name_found)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_INVALID,
                                "ncm_prior_flat_param_new_name: invalid parameter name `%s'.", name);
    NCM_UTIL_ON_ERROR_RETURN (error,
                              g_free (model_ns);
                              g_free (param_name), NULL);
  }

  {
    NcmPriorFlatParam *pfp = g_object_new (NCM_TYPE_PRIOR_FLAT_PARAM,
                                           "x-low", x_low,
                                           "x-upp", x_upp,
                                           "scale", scale,
                                           "model-ns", model_ns,
                                           "stack-pos", stack_pos,
                                           "parameter-name", param_name,
                                           NULL);

    g_free (model_ns);
    g_free (param_name);

    return pfp;
  }
}

/**
 * ncm_prior_flat_param_ref:
 * @pfp: a #NcmPriorFlatParam
 *
 * Increases the reference count of @pfp atomically.
 *
 * Returns: (transfer full): @pfp.
 */
NcmPriorFlatParam *
ncm_prior_flat_param_ref (NcmPriorFlatParam *pfp)
{
  return g_object_ref (pfp);
}

/**
 * ncm_prior_flat_param_free:
 * @pfp: a #NcmPriorFlatParam
 *
 * Decreases the reference count of @pfp atomically.
 *
 */
void
ncm_prior_flat_param_free (NcmPriorFlatParam *pfp)
{
  g_object_unref (pfp);
}

/**
 * ncm_prior_flat_param_clear:
 * @pfp: a #NcmPriorFlatParam
 *
 * Decreases the reference count of *@pfp and sets *@pfp to NULL.
 *
 */
void
ncm_prior_flat_param_clear (NcmPriorFlatParam **pfp)
{
  g_clear_object (pfp);
}

/**
 * ncm_prior_flat_param_set_model_ns:
 * @pfp: a #NcmPriorFlatParam
 * @model_ns: model namespace
 * @error: a #GError
 *
 * Sets the model namespace of @pfp to @model_ns.
 *
 */
void
ncm_prior_flat_param_set_model_ns (NcmPriorFlatParam *pfp, const gchar *model_ns, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (pfp));
  {
    GType model_type = g_type_from_name (model_ns); /* check if the model namespace is valid */

    if (model_type == G_TYPE_INVALID)
    {
      g_error ("ncm_prior_flat_param_set_model_ns: invalid model namespace `%s'.", model_ns);
    }
    else
    {
      NcmModelID mid = ncm_model_id_by_type (model_type, error);

      if (error && *error)
        return;

      g_clear_pointer (&pfp->model_ns, g_free);
      pfp->model_ns = g_strdup (model_ns);
      pfp->mid      = mid;
    }
  }
}

/**
 * ncm_prior_flat_param_set_stack_pos:
 * @pfp: a #NcmPriorFlatParam
 * @stack_pos: stack position
 *
 * Sets the stack position of @pfp to @stack_pos.
 *
 */
void
ncm_prior_flat_param_set_stack_pos (NcmPriorFlatParam *pfp, guint stack_pos)
{
  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (pfp));
  {
    pfp->stack_pos = stack_pos;
  }
}

/**
 * ncm_prior_flat_param_set_param_name:
 * @pfp: a #NcmPriorFlatParam
 * @param_name: parameter name
 *
 * Sets the parameter name of @pfp to @param_name.
 *
 */
void
ncm_prior_flat_param_set_param_name (NcmPriorFlatParam *pfp, const gchar *param_name)
{
  g_return_if_fail (NCM_IS_PRIOR_FLAT_PARAM (pfp));
  {
    g_clear_pointer (&pfp->param_name, g_free);

    pfp->param_name = g_strdup (param_name);
  }
}

/**
 * ncm_prior_flat_param_peek_model_ns:
 * @pfp: a #NcmPriorFlatParam
 *
 * Returns: (transfer none): the model namespace of @pfp.
 */
const gchar *
ncm_prior_flat_param_peek_model_ns (NcmPriorFlatParam *pfp)
{
  return pfp->model_ns;
}

/**
 * ncm_prior_flat_param_get_stack_pos:
 * @pfp: a #NcmPriorFlatParam
 *
 * Gets the stack position of the parameter.
 *
 * Returns: the stack position of @pfp.
 */
guint
ncm_prior_flat_param_get_stack_pos (NcmPriorFlatParam *pfp)
{
  return pfp->stack_pos;
}

/**
 * ncm_prior_flat_param_peek_param_name:
 * @pfp: a #NcmPriorFlatParam
 *
 * Returns: (transfer none): the parameter name of @pfp.
 */
const gchar *
ncm_prior_flat_param_peek_param_name (NcmPriorFlatParam *pfp)
{
  return pfp->param_name;
}


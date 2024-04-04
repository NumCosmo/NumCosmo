/***************************************************************************
 *            ncm_prior_gauss_param.c
 *
 *  Wed August 03 10:55:21 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_gauss_param.c
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
 * SECTION:ncm_prior_gauss_param
 * @title: NcmPriorGaussParam
 * @short_description: a gaussian prior on a parameter
 *
 * This object is a subclass of #NcmPriorGauss, precisely designed for a Gaussian
 * prior on a parameter. The prior relies on a parameter ID and
 * user-specified mean and standard deviation parameters.
 *
 * Users have flexibility in specifying the parameter in various ways:
 * - Using the pair NcmModelID and the parameter pid.
 * - Providing a single NcmMSetPIndex.
 * - Supplying a string consisting of a parameter full name "model:parameter".
 *
 * These options provide versatile and convenient methods for parameter identification.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_gauss_param.h"

enum
{
  PROP_0,
  PROP_MODEL_NS,
  PROP_STACK_POS,
  PROP_PARAM_NAME,
};

struct _NcmPriorGaussParam
{
  /*< private >*/
  NcmPriorGauss parent_instance;
  gchar *model_ns;
  gchar *param_name;
  guint stack_pos;
  NcmModelID mid;
};

G_DEFINE_TYPE (NcmPriorGaussParam, ncm_prior_gauss_param, NCM_TYPE_PRIOR_GAUSS)

static void
ncm_prior_gauss_param_init (NcmPriorGaussParam *pgp)
{
}

static void
_ncm_prior_gauss_param_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorGaussParam *pgp = NCM_PRIOR_GAUSS_PARAM (object);

  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (object));

  switch (prop_id)
  {
    case PROP_MODEL_NS:
      ncm_prior_gauss_param_set_model_ns (pgp, g_value_get_string (value));
      break;
    case PROP_STACK_POS:
      ncm_prior_gauss_param_set_stack_pos (pgp, g_value_get_uint (value));
      break;
    case PROP_PARAM_NAME:
      ncm_prior_gauss_param_set_param_name (pgp, g_value_get_string (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_gauss_param_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorGaussParam *pgp = NCM_PRIOR_GAUSS_PARAM (object);

  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (object));

  switch (prop_id)
  {
    case PROP_MODEL_NS:
      g_value_set_string (value, ncm_prior_gauss_param_peek_model_ns (pgp));
      break;
    case PROP_STACK_POS:
      g_value_set_uint (value, ncm_prior_gauss_param_get_stack_pos (pgp));
      break;
    case PROP_PARAM_NAME:
      g_value_set_string (value, ncm_prior_gauss_param_peek_param_name (pgp));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_gauss_param_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_gauss_param_parent_class)->finalize (object);
}

static gdouble _ncm_prior_gauss_param_mean (NcmPriorGauss *pg, NcmMSet *mset);

static void
ncm_prior_gauss_param_class_init (NcmPriorGaussParamClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmPriorGaussClass *pg_class = NCM_PRIOR_GAUSS_CLASS (klass);

  object_class->set_property = &_ncm_prior_gauss_param_set_property;
  object_class->get_property = &_ncm_prior_gauss_param_get_property;
  object_class->finalize     = &_ncm_prior_gauss_param_finalize;

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

  pg_class->mean = &_ncm_prior_gauss_param_mean;
}

static gdouble
_ncm_prior_gauss_param_mean (NcmPriorGauss *pg, NcmMSet *mset)
{
  NcmPriorGaussParam *pgp = NCM_PRIOR_GAUSS_PARAM (pg);
  NcmModel *model         = ncm_mset_peek (mset, pgp->mid);

  return ncm_model_param_get_by_name (model, pgp->param_name);
}

/**
 * ncm_prior_gauss_param_new:
 * @model: a #NcmModel
 * @pid: parameter id
 * @mu: mean
 * @sigma: standard deviation
 *
 * Creates a new Gaussian prior for parameter @pid of @model.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_new (NcmModel *model, guint pid, gdouble mu, gdouble sigma)
{
  const gchar *model_ns   = G_OBJECT_TYPE_NAME (model);
  const gchar *param_name = ncm_model_param_name (model, pid);
  NcmPriorGaussParam *pgp = g_object_new (NCM_TYPE_PRIOR_GAUSS_PARAM,
                                          "mu", mu,
                                          "sigma", sigma,
                                          "model-ns", model_ns,
                                          "parameter-name", param_name,
                                          NULL);

  return pgp;
}

/**
 * ncm_prior_gauss_param_new_name:
 * @name: parameter name
 * @mu: mean
 * @sigma: standard deviation
 *
 * Creates a new Gaussian prior for parameter named @name in @mset. See
 * ncm_mset_split_full_name() for details on the parameter name format.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_new_name (const gchar *name, gdouble mu, gdouble sigma)
{
  gchar *model_ns   = NULL;
  gchar *param_name = NULL;
  guint stack_pos   = 0;

  if (!ncm_mset_split_full_name (name, &model_ns, &stack_pos, &param_name))
  {
    g_error ("ncm_prior_gauss_param_new_name: invalid parameter name `%s'.", name);

    return NULL;
  }

  {
    NcmPriorGaussParam *pgp = g_object_new (NCM_TYPE_PRIOR_GAUSS_PARAM,
                                            "mu", mu,
                                            "sigma", sigma,
                                            "model-ns", model_ns,
                                            "stack-pos", stack_pos,
                                            "parameter-name", param_name,
                                            NULL);

    g_free (model_ns);
    g_free (param_name);

    return pgp;
  }
}

/**
 * ncm_prior_gauss_param_ref:
 * @pgp: a #NcmPriorGaussParam
 *
 * Increases the reference count of @pgp atomically.
 *
 * Returns: (transfer full): @pgp.
 */
NcmPriorGaussParam *
ncm_prior_gauss_param_ref (NcmPriorGaussParam *pgp)
{
  return g_object_ref (pgp);
}

/**
 * ncm_prior_gauss_param_free:
 * @pgp: a #NcmPriorGaussParam
 *
 * Decreases the reference count of @pgp atomically.
 *
 */
void
ncm_prior_gauss_param_free (NcmPriorGaussParam *pgp)
{
  g_object_unref (pgp);
}

/**
 * ncm_prior_gauss_param_clear:
 * @pgp: a #NcmPriorGaussParam
 *
 * Decreases the reference count of *@pgp and sets *@pgp to NULL.
 *
 */
void
ncm_prior_gauss_param_clear (NcmPriorGaussParam **pgp)
{
  g_clear_object (pgp);
}

/**
 * ncm_prior_gauss_param_set_model_ns:
 * @pgp: a #NcmPriorGaussParam
 * @model_ns: model namespace
 *
 * Sets the model namespace of @pgp to @model_ns.
 *
 */
void
ncm_prior_gauss_param_set_model_ns (NcmPriorGaussParam *pgp, const gchar *model_ns)
{
  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (pgp));
  {
    GType model_type = g_type_from_name (model_ns); /* check if the model namespace is valid */

    if (model_type == G_TYPE_INVALID)
    {
      g_error ("ncm_prior_gauss_param_set_model_ns: invalid model namespace `%s'.", model_ns);
    }
    else
    {
      NcmModelID mid = ncm_model_id_by_type (model_type);

      g_clear_pointer (&pgp->model_ns, g_free);
      pgp->model_ns = g_strdup (model_ns);
      pgp->mid      = mid;
    }
  }
}

/**
 * ncm_prior_gauss_param_set_stack_pos:
 * @pgp: a #NcmPriorGaussParam
 * @stack_pos: stack position
 *
 * Sets the stack position of @pgp to @stack_pos.
 *
 */
void
ncm_prior_gauss_param_set_stack_pos (NcmPriorGaussParam *pgp, guint stack_pos)
{
  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (pgp));
  {
    pgp->stack_pos = stack_pos;
  }
}

/**
 * ncm_prior_gauss_param_set_param_name:
 * @pgp: a #NcmPriorGaussParam
 * @param_name: parameter name
 *
 * Sets the parameter name of @pgp to @param_name.
 *
 */
void
ncm_prior_gauss_param_set_param_name (NcmPriorGaussParam *pgp, const gchar *param_name)
{
  g_return_if_fail (NCM_IS_PRIOR_GAUSS_PARAM (pgp));
  {
    g_clear_pointer (&pgp->param_name, g_free);

    pgp->param_name = g_strdup (param_name);
  }
}

/**
 * ncm_prior_gauss_param_peek_model_ns:
 * @pgp: a #NcmPriorGaussParam
 *
 * Returns: (transfer none): the model namespace of @pgp.
 */
const gchar *
ncm_prior_gauss_param_peek_model_ns (NcmPriorGaussParam *pgp)
{
  return pgp->model_ns;
}

/**
 * ncm_prior_gauss_param_get_stack_pos:
 * @pgp: a #NcmPriorGaussParam
 *
 * Gets the stack position of the parameter.
 *
 * Returns: the stack position of @pgp.
 */
guint
ncm_prior_gauss_param_get_stack_pos (NcmPriorGaussParam *pgp)
{
  return pgp->stack_pos;
}

/**
 * ncm_prior_gauss_param_peek_param_name:
 * @pgp: a #NcmPriorGaussParam
 *
 * Returns: (transfer none): the parameter name of @pgp.
 */
const gchar *
ncm_prior_gauss_param_peek_param_name (NcmPriorGaussParam *pgp)
{
  return pgp->param_name;
}


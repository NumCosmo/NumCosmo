/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_model_mvnd.c
 *
 *  Sun February 04 15:31:31 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_model_mvnd.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_model_mvnd
 * @title: NcmModelMVND
 * @short_description: Multivariate Normal Distribution mean model.
 *
 * Multivariate Normal distribution model of the mean.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_mvnd.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_DIM,
  PROP_SIZE,
};

struct _NcmModelMVNDPrivate
{
  gint dim;
  NcmVector *mu;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmModelMVND, ncm_model_mvnd, NCM_TYPE_MODEL);

static void
ncm_model_mvnd_init (NcmModelMVND *model_mvnd)
{
  model_mvnd->priv      = ncm_model_mvnd_get_instance_private (model_mvnd);
  model_mvnd->priv->dim = 0;
  model_mvnd->priv->mu  = NULL;
}

static void
_ncm_model_mvnd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModelMVND *model_mvnd = NCM_MODEL_MVND (object);
  g_return_if_fail (NCM_IS_MODEL_MVND (object));

  switch (prop_id)
  {
    case PROP_DIM:
      model_mvnd->priv->dim = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_mvnd_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModelMVND *model_mvnd = NCM_MODEL_MVND (object);
  g_return_if_fail (NCM_IS_MODEL_MVND (object));

  switch (prop_id)
  {
    case PROP_DIM:
      g_value_set_uint (value, model_mvnd->priv->dim);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_mvnd_dispose (GObject *object)
{
  NcmModelMVND *model_mvnd = NCM_MODEL_MVND (object);
  
  ncm_vector_clear (&model_mvnd->priv->mu);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_mvnd_parent_class)->dispose (object);
}

static void
_ncm_model_mvnd_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_mvnd_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (ncm_model_mvnd, NCM_TYPE_MODEL_MVND);

static void
ncm_model_mvnd_class_init (NcmModelMVNDClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_ncm_model_mvnd_set_property;
  model_class->get_property = &_ncm_model_mvnd_get_property;
  
  object_class->dispose      = &_ncm_model_mvnd_dispose;
  object_class->finalize     = &_ncm_model_mvnd_finalize;

  ncm_model_class_set_name_nick (model_class, "MVND", "NcmModelMVND");
  ncm_model_class_add_params (model_class, 0, NNCM_MODEL_MVND_VPARAM_LEN, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcmModelMVND",
                              "MVND",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_vparam (model_class, NCM_MODEL_MVND_MEAN, 1, "mu", "mu",
                              -50.0, 50.0, 1.0, 0.0, 0.0, NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);

  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dim",
                                                      NULL,
                                                      "Problem dimension",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_model_mvnd_new:
 * @dim: dimension of the MVND
 * 
 * Creates a new MVND mean model of @dim dimensions.
 * 
 * Returns: (transfer full): the newly created #NcmModelMVND
 */
NcmModelMVND *
ncm_model_mvnd_new (const guint dim)
{
  NcmModelMVND *gauss_mvnd = g_object_new (NCM_TYPE_MODEL_MVND,
                                           "dim", dim,
                                           "mu-length", dim,
                                           NULL);
  return gauss_mvnd;
}

/**
 * ncm_model_mvnd_ref:
 * @model_mvnd: a #NcmModelMVND
 * 
 * Increases the reference count of @model_mvnd by one.
 * 
 * Returns: (transfer full): @model_mvnd
 */
NcmModelMVND *
ncm_model_mvnd_ref (NcmModelMVND *model_mvnd)
{
  return g_object_ref (model_mvnd);
}

/**
 * ncm_model_mvnd_free:
 * @model_mvnd: a #NcmModelMVND
 * 
 * Decreases the reference count of @model_mvnd by one.
 * 
 */
void 
ncm_model_mvnd_free (NcmModelMVND *model_mvnd)
{
  g_object_unref (model_mvnd);
}

/**
 * ncm_model_mvnd_clear:
 * @model_mvnd: a #NcmModelMVND
 * 
 * If @model_mvnd is different from NULL, decreases the reference count of 
 * @model_mvnd by one and sets @model_mvnd to NULL.
 * 
 */
void 
ncm_model_mvnd_clear (NcmModelMVND **model_mvnd)
{
  g_clear_object (model_mvnd);
}

/**
 * ncm_model_mvnd_mean:
 * @model_mvnd: a #NcmModelMVND
 * @y: a #NcmVector
 * 
 * Copies into @y the mean vector.
 * 
 */
void 
ncm_model_mvnd_mean (NcmModelMVND *model_mvnd, NcmVector *y)
{
  if (model_mvnd->priv->mu == NULL)
  {
    NcmModel *model     = NCM_MODEL (model_mvnd);
    const guint mu_size = ncm_model_vparam_len (model, NCM_MODEL_MVND_MEAN);
    const guint mu_i    = ncm_model_vparam_index (model, NCM_MODEL_MVND_MEAN, 0);

    model_mvnd->priv->mu = ncm_vector_get_subvector (model->params, mu_i, mu_size);
  }

  ncm_vector_memcpy (y, model_mvnd->priv->mu);
}

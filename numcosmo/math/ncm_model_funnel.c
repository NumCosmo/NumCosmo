/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_model_funnel.c
 *
 *  Wed May 12 21:24:12 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_funnel.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_model_funnel
 * @title: NcmModelFunnel
 * @short_description: Multivariate Normal Distribution mean model.
 *
 * Multivariate Normal distribution model of the mean.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_funnel.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SIZE,
};

struct _NcmModelFunnelPrivate
{
  gint place_holder;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmModelFunnel, ncm_model_funnel, NCM_TYPE_MODEL)

static void
ncm_model_funnel_init (NcmModelFunnel *model_funnel)
{
  model_funnel->priv = ncm_model_funnel_get_instance_private (model_funnel);
}

static void
_ncm_model_funnel_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcmModelFunnel *model_funnel = NCM_MODEL_FUNNEL (object);*/
  g_return_if_fail (NCM_IS_MODEL_FUNNEL (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_funnel_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcmModelFunnel *model_funnel = NCM_MODEL_FUNNEL (object);*/
  g_return_if_fail (NCM_IS_MODEL_FUNNEL (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_funnel_dispose (GObject *object)
{
  /*NcmModelFunnel *model_funnel = NCM_MODEL_FUNNEL (object);*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_funnel_parent_class)->dispose (object);
}

static void
_ncm_model_funnel_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_funnel_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (ncm_model_funnel, NCM_TYPE_MODEL_FUNNEL);

static void
ncm_model_funnel_class_init (NcmModelFunnelClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_ncm_model_funnel_set_property;
  model_class->get_property = &_ncm_model_funnel_get_property;
  
  object_class->dispose      = &_ncm_model_funnel_dispose;
  object_class->finalize     = &_ncm_model_funnel_finalize;

  ncm_model_class_set_name_nick (model_class, "MFU", "NcmModelFunnel");
  ncm_model_class_add_params (model_class, NNCM_MODEL_FUNNEL_SPARAM_LEN, NNCM_MODEL_FUNNEL_VPARAM_LEN, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcmModelFunnel",
                              "MFU",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_sparam (model_class, NCM_MODEL_FUNNEL_NU, "\\nu", "nu",
                              -1000.0, 1000.0, 1.0e-1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_vparam (model_class, NCM_MODEL_FUNNEL_X, 9, "x", "x",
                              -1.0e6, 1.0e6, 1.0, 0.0, 0.0, NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);
}

/**
 * ncm_model_funnel_new:
 * @n: number of $x$ variables
 * 
 * Creates a new Funnel model.
 * 
 * Returns: (transfer full): the newly created #NcmModelFunnel
 */
NcmModelFunnel *
ncm_model_funnel_new (guint n)
{
  NcmModelFunnel *mfu = g_object_new (NCM_TYPE_MODEL_FUNNEL,
                                      "x-length", n,
                                      NULL);
  return mfu;
}

/**
 * ncm_model_funnel_ref:
 * @mfu: a #NcmModelFunnel
 * 
 * Increases the reference count of @mfu by one.
 * 
 * Returns: (transfer full): @mfu
 */
NcmModelFunnel *
ncm_model_funnel_ref (NcmModelFunnel *mfu)
{
  return g_object_ref (mfu);
}

/**
 * ncm_model_funnel_free:
 * @mfu: a #NcmModelFunnel
 * 
 * Decreases the reference count of @mfu by one.
 * 
 */
void 
ncm_model_funnel_free (NcmModelFunnel *mfu)
{
  g_object_unref (mfu);
}

/**
 * ncm_model_funnel_clear:
 * @mfu: a #NcmModelFunnel
 * 
 * If @mfu is different from NULL, decreases the reference count of
 * @mfu by one and sets @mfu to NULL.
 * 
 */
void 
ncm_model_funnel_clear (NcmModelFunnel **mfu)
{
  g_clear_object (mfu);
}

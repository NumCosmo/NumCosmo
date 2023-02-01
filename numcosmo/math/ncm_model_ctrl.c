/***************************************************************************
 *            ncm_model_ctrl.c
 *
 *  Mon February 27 12:10:09 2012
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
 * SECTION:ncm_model_ctrl
 * @title: NcmModelCtrl
 * @short_description: Control object for testing updates on model status.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_ctrl.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_MODEL,
};

G_DEFINE_TYPE (NcmModelCtrl, ncm_model_ctrl, G_TYPE_OBJECT);

static void
ncm_model_ctrl_init (NcmModelCtrl *ctrl)
{
  g_weak_ref_init (&ctrl->model_wr, NULL);
  ctrl->pkey = 0;

  ctrl->submodel_ctrl = g_ptr_array_new ();
  g_ptr_array_set_free_func (ctrl->submodel_ctrl, (GDestroyNotify) ncm_model_ctrl_free);

  ctrl->submodel_last_update = g_array_new (TRUE, TRUE, sizeof (gboolean));
}

static void
ncm_model_ctrl_dispose (GObject *object)
{
  NcmModelCtrl *ctrl = NCM_MODEL_CTRL (object);

  g_weak_ref_clear (&ctrl->model_wr);

  g_clear_pointer (&ctrl->submodel_ctrl, g_ptr_array_unref);
  g_clear_pointer (&ctrl->submodel_last_update, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_ctrl_parent_class)->dispose (object);
}

static void
ncm_model_ctrl_finalize (GObject *object)
{
  /*NcmModelCtrl *ctrl = NCM_MODEL_CTRL (object);*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_ctrl_parent_class)->finalize (object);
}

static void
ncm_model_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModelCtrl *ctrl = NCM_MODEL_CTRL (object);
  g_return_if_fail (NCM_IS_MODEL_CTRL (object));

  switch (prop_id)
  {
    case PROP_MODEL:
      ncm_model_ctrl_set_model (ctrl, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_model_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModelCtrl *ctrl = NCM_MODEL_CTRL (object);
  g_return_if_fail (NCM_IS_MODEL_CTRL (object));

  switch (prop_id)
  {
    case PROP_MODEL:
      g_value_take_object (value, g_weak_ref_get (&ctrl->model_wr));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_model_ctrl_class_init (NcmModelCtrlClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->set_property = ncm_model_set_property;
  object_class->get_property = ncm_model_get_property;
  object_class->dispose      = ncm_model_ctrl_dispose;
  object_class->finalize     = ncm_model_ctrl_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MODEL,
                                   g_param_spec_object ("model",
                                                        NULL,
                                                        "Last Model used",
                                                        NCM_TYPE_MODEL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_model_ctrl_new:
 * @model: (allow-none): FIXME
   *
 * FIXME
 *
 * Returns: FIXME
 */
NcmModelCtrl *
ncm_model_ctrl_new (NcmModel *model)
{
  if (model)
    return g_object_new (NCM_TYPE_MODEL_CTRL, "model", model, NULL);
  else
    return g_object_new (NCM_TYPE_MODEL_CTRL, NULL);
}

/**
 * ncm_model_ctrl_update:
 * @ctrl: a #NcmModelCtrl
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_ctrl_model_update:
 * @ctrl: a #NcmModelCtrl
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_model_ctrl_get_model:
 * @ctrl: a #NcmModelCtrl
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
/**
 * ncm_model_ctrl_model_last_update:
 * @ctrl: a #NcmModelCtrl
 * 
 * Checks if the main model was updated during the last call to
 * ncm_model_ctrl_update().
 * 
 * Returns: TRUE if the main model was updated.
 */
/**
 * ncm_model_ctrl_model_has_submodel:
 * @ctrl: a #NcmModelCtrl
 * @mid: a @NcmModelID
 * 
 * Checks if there is a submode inside ctrl model, it is an
 * error to call this function in an empty @ctrl.
 *
 * Returns: TRUE if there is a submodel with @mid inside the ctrl model.
 */
/**
 * ncm_model_ctrl_submodel_last_update:
 * @ctrl: a #NcmModelCtrl
 * @mid: a @NcmModelID
 *
 * Checks if the submodel @mid was updated during the last call to
 * ncm_model_ctrl_update().
 *
 * Returns: TRUE if the submodel with @mid inside the ctrl model was updated.
 */

/**
 * ncm_model_ctrl_set_model:
 * @ctrl: a #NcmModelCtrl
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_model_ctrl_set_model (NcmModelCtrl *ctrl, NcmModel *model)
{
  gboolean up = FALSE;
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);

  if (model != ctrl_model)
  {
    g_weak_ref_set (&ctrl->model_wr, model);
    ctrl->pkey  = model->pkey;
    up          = TRUE;
  }

  {
    const guint n = ncm_model_get_submodel_len (model);
    guint i;
    for (i = 0; i < n; i++)
    {
      NcmModel *submodel = ncm_model_peek_submodel (model, i);
      if (i >= ctrl->submodel_ctrl->len)
      {
        NcmModelCtrl *sub_ctrl = ncm_model_ctrl_new (submodel);
        g_ptr_array_add (ctrl->submodel_ctrl, sub_ctrl);
        up = TRUE;
      }
      else
      {
        NcmModelCtrl *sub_ctrl = g_ptr_array_index (ctrl->submodel_ctrl, i);
        if (ncm_model_ctrl_set_model (sub_ctrl, submodel))
        {
          up = TRUE;
        }
      }
    }
    g_ptr_array_set_size (ctrl->submodel_ctrl, n);
  }
  
  ncm_model_clear (&ctrl_model);
  return up;
}

/**
 * ncm_model_ctrl_has_model:
 * @ctrl: a #NcmModelCtrl
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_model_ctrl_has_model (NcmModelCtrl *ctrl, NcmModel *model)
{
  NcmModel *ctrl_model = ncm_model_ctrl_get_model (ctrl);
  gboolean is_model = (model == ctrl_model);

  ncm_model_clear (&ctrl_model);

  return is_model;
}

/**
 * ncm_model_ctrl_force_update:
 * @ctrl: a #NcmModelCtrl
 *
 * FIXME
 *
 */
void
ncm_model_ctrl_force_update (NcmModelCtrl *ctrl)
{
  g_weak_ref_set (&ctrl->model_wr, NULL);
  g_ptr_array_set_size (ctrl->submodel_ctrl, 0);
  g_array_set_size (ctrl->submodel_last_update, 0);
  return;
}

/**
 * ncm_model_ctrl_free:
 * @ctrl: a #NcmModelCtrl
 *
 * FIXME
 *
 */
void
ncm_model_ctrl_free (NcmModelCtrl *ctrl)
{
  g_object_unref (ctrl);
}

/**
 * ncm_model_ctrl_clear:
 * @ctrl: a #NcmModelCtrl
 *
 * FIXME
 *
 */
void
ncm_model_ctrl_clear (NcmModelCtrl **ctrl)
{
  g_clear_object (ctrl);
}
